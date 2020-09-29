!$Id$
!BOP
!
! !ROUTINE: prepoststep_ens --- Used-defined Pre/Poststep routine for PDAF
!
! !INTERFACE:
SUBROUTINE prepoststep_ens_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
! 
! The routine is called for global filters (e.g. SEIK)
! before the analysis and after the ensemble transformation.
! For local filters (e.g. LSEIK) the routine is called
! before and after the loop over all local analysis
! domains.
! The routine provides full access to the state 
! estimate and the state ensemble to the user.
! Thus, user-controlled pre- and poststep 
! operations can be performed here. For example 
! the forecast and the analysis states and ensemble
! covariance matrix can be analized, e.g. by 
! computing the estimated variances. 
! For the offline mode, this routine is the place
! in which the writing of the analysis ensemble
! can be performed.
!
! If a user considers to perform adjustments to the 
! estimates (e.g. for balances), this routine is 
! the right place for it.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code based on offline_1D
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: state_min_p, dim_state, dim_state_p, local_dims, epoch
  USE mod_parallel_pdaf, &
       ONLY: mype_filter, npes_filter, COMM_filter, MPI_DOUBLE_PRECISION, &
       MPIerr, MPIstatus, MPI_MODE_CREATE, MPI_OFFSET_KIND, &
       MPI_INFO_NULL, MPI_MODE_WRONLY, MPI_STATUS_IGNORE, MPI_SUM
  USE my_state_accessors, &
          ONLY: current_step

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step (not relevant for offline mode)
  INTEGER, INTENT(in) :: dim_p       ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens     ! Size of state ensemble
  INTEGER, INTENT(in) :: dim_ens_p   ! PE-local size of ensemble
  INTEGER, INTENT(in) :: dim_obs_p   ! PE-local dimension of observation vector
  REAL, INTENT(inout) :: state_p(dim_p) ! PE-local forecast/analysis state
  ! The array 'state_p' is not generally not initialized in the case of SEIK.
  ! It can be used freely here.
  REAL, INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) ! Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)      ! PE-local state ensemble
  INTEGER, INTENT(in) :: flag        ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_get_state      (as U_prepoststep)
! Called by: PDAF_X_update       (as U_prepoststep)
! Calls: MPI_send
! Calls: MPI_recv
!EOP

! *** local variables ***
  INTEGER :: i, j, member, domain      ! counters
  INTEGER, SAVE :: allocflag = 0       ! Flag for memory counting
  LOGICAL, SAVE :: firstio = .TRUE.    ! File output is peformed for first time?
  LOGICAL, SAVE :: firsttime = .TRUE.    ! Routine is called for first time?
  REAL :: invdim_ens                   ! Inverse ensemble size
  REAL :: invdim_ensm1                 ! Inverse of ensemble size minus 1
  REAL :: rmserror_est                 ! estimated RMS error
  REAL :: rmserror_est_p                 ! estimated RMS error
  REAL, ALLOCATABLE :: variance_p(:)     ! model state variances
  CHARACTER(len=5) :: ensstr          ! String for ensemble member
  CHARACTER(len=5) :: mpestr          ! String for ensemble member
  CHARACTER(len=5) :: epostr          ! String for ensemble member
  ! Variables for parallelization - global fields
  INTEGER :: offset   ! Row-offset according to domain decomposition
  REAL, ALLOCATABLE :: variance(:)     ! local variance
  REAL, ALLOCATABLE :: ens(:,:)       ! global ensemble
  REAL, ALLOCATABLE :: state(:)       ! global state vector
  REAL,ALLOCATABLE :: ens_p_tmp(:,:) ! Temporary ensemble for some PE-domain
  REAL,ALLOCATABLE :: state_p_tmp(:) ! Temporary state for some PE-domain
  integer(kind=MPI_OFFSET_KIND) :: disp
  INTEGER       :: file_id      ! MPI file handle
  INTEGER       :: ierr         ! MPI error handle

! **********************
! *** INITIALIZATION ***
! **********************

  IF (mype_filter == 0) THEN
     IF (firsttime) THEN
        WRITE (*, '(8x, a)') 'Analyze forecasted state ensemble'
     ELSE
        WRITE (*, '(8x, a)') 'Analyze and write assimilated state ensemble'
     END IF
  END IF
  ! Allocate fields
  ALLOCATE(variance_p(dim_p))

  ! Initialize numbers
  rmserror_est_p  = 0.0
  rmserror_est  = 0.0
  invdim_ens    = 1.0 / REAL(dim_ens)  
  invdim_ensm1  = 1.0 / REAL(dim_ens - 1)

! **************************************************************
! *** Perform prepoststep for SEIK with re-inititialization. ***
! *** The state and error information is completely in the   ***
! *** ensemble.                                              ***
! *** Also performed for SEIK without re-init at the initial ***
! *** time.                                                  ***
! **************************************************************

  ! *** Compute mean state
  IF (mype_filter == 0) WRITE (*, '(8x, a)') '--- compute ensemble mean'

  state_p = 0.0
  DO member = 1, dim_ens
     DO i = 1, dim_p
        state_p(i) = state_p(i) + ens_p(i, member)
     END DO
  END DO
  state_p(:) = invdim_ens * state_p(:)

  ! *** Compute sampled variances ***
  variance_p(:) = 0.0
  DO member = 1, dim_ens
     DO j = 1, dim_p
        variance_p(j) = variance_p(j) &
             + (ens_p(j, member) - state_p(j)) &
             * (ens_p(j, member) - state_p(j))
     END DO
  END DO
  variance_p(:) = invdim_ensm1 * variance_p(:)

! ************************************************************
! *** Compute RMS errors according to sampled covar matrix ***
! ************************************************************

  ! total estimated RMS error
  
  DO i = 1, dim_p
      rmserror_est_p = rmserror_est_p + variance_p(i)
  ENDDO

! collect values on rank 0
  
  call mpi_reduce( rmserror_est_p, rmserror_est, 1, &
    MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMM_filter, MPIerr )

! *****************
! *** Screen IO ***
! *****************

!  ! Output RMS errors given by sampled covar matrix
  IF (mype_filter == 0) THEN
    rmserror_est = SQRT(rmserror_est / dim_state)
     WRITE (*, '(12x, a, es12.4)') &
       'RMS error according to sampled variance: ', rmserror_est
  END IF
  
  DEALLOCATE(variance_p)
 
! *******************
! *** File output ***
! *******************

  notfirst: IF (.not. firsttime) THEN
    
    !disp = (state_min_p-1)*sizeof(ens_p(1,member))
    !
    !do member=1, dim_ens
		!	WRITE (ensstr, '(i5.5)') member
    !  call MPI_FILE_OPEN(COMM_filter, 'ens_'//TRIM(ensstr)//'_ana.txt', & 
    !                   MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
    !                   MPI_INFO_NULL, file_id, ierr) 
    !  call MPI_FILE_SET_VIEW(file_id, disp , MPI_DOUBLE_PRECISION, & 
    !                       MPI_DOUBLE_PRECISION, 'native', & 
    !                       MPI_INFO_NULL, ierr) 
    !  call MPI_FILE_WRITE(file_id, ens_p(1,member), dim_p, MPI_DOUBLE_PRECISION, & 
    !                    MPI_STATUS_IGNORE, ierr) 
    !  call MPI_FILE_CLOSE(file_id, ierr)   
    !end do

! write analysis state in ascii
	
    write(mpestr,'(i5.5)') mype_filter
    write(epostr,'(i5.5)') current_step
    if (step<0) then
      open(10, &
        file='for_state_rank'//TRIM(mpestr)//'_epoch'//TRIM(epostr)//'.txt', &
        form='formatted')
    else
      open(10, &
        file='ana_state_rank'//TRIM(mpestr)//'_epoch'//TRIM(epostr)//'.txt', &
        form='formatted')
    end if
    do i=1,dim_p
      write(10,"(I5, TR2, es12.4)") state_min_p -1 + i, state_p(i)
    end do
    close(10)

  END IF notfirst

  IF (firsttime) THEN
    
    write(mpestr,'(i5.5)') mype_filter
    write(epostr,'(i5.5)') current_step
    open(10, &
      file='initial_state_rank'//TRIM(mpestr)//'_epoch'//TRIM(epostr)//'.txt', &
      form='formatted')
    do i=1,dim_p
      write(10,"(I5, TR2, es12.4)") state_min_p -1 + i, state_p(i)
    end do
    close(10)
  
  END IF

! ********************
! *** finishing up ***
! ********************

  firsttime = .FALSE.

END SUBROUTINE prepoststep_ens_pdaf
