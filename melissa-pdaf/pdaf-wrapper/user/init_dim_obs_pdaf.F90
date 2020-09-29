!$Id$
!BOP
!
! !ROUTINE: init_dim_obs_pdaf --- Compute number of observations
!
! !INTERFACE:
SUBROUTINE init_dim_obs_pdaf(step, dim_obs_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/ETKF/ESTKF
!
! The routine is called at the beginning of each
! analysis step.  It has to initialize the size of 
! the observation vector according to the current 
! time step for the PE-local domain.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
      ONLY :  dim_state, &      ! dimension of full state
              obs_blk_size, &   ! block size for observations
              obs_prcnt, &      ! percentage of full state = number of observations
              state_min_p, &    ! minimum index state on PE
              state_max_p, &    ! maximum index state on PE
              obs_index_p, &       ! index array of observations
              obs_p, epoch
  USE mod_parallel_pdaf, &
        ONLY: mype_filter, COMM_filter, MPI_DOUBLE_PRECISION, &
        MPI_INFO_NULL, MPI_MODE_RDONLY, MPI_STATUS_IGNORE, &
        MPI_OFFSET_KIND, npes_filter, MPI_INTEGER
  USE mod_parallel_model, &
      ONLY: COMM_world, mpierr

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step       ! Current time step
  INTEGER, INTENT(out) :: dim_obs_p  ! Dimension of observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_seek_analysis    (as U_init_dim_obs)
! Called by: PDAF_seik_analysis, PDAF_seik_analysis_newT
! Called by: PDAF_enkf_analysis_rlm, PDAF_enkf_analysis_rsm
! Called by: PDAF_etkf_analysis, PDAF_etkf_analysis_T
! Called by: PDAF_estkf_analysis, PDAF_estkf_analysis_fixed
! Called by: PDAF_netf_analysis
!EOP

! *** Local variables
  INTEGER               :: i, j                     ! Counters
  INTEGER               :: num_reg    ! number of regions
  INTEGER               :: dim_obs    ! number of observations global
  INTEGER               :: stride     ! stride for the regions
  INTEGER               :: offset     ! offset for index in loop
  INTEGER               :: index_tmp  ! index dummy var
  INTEGER               :: cnt_obs_p  ! counter for observations on PE
  INTEGER               :: cnt_obs    ! counter for observations on PE
  INTEGER               :: file_id      ! MPI file handle
  INTEGER               :: ierr         ! MPI error handle
  character(len=5)      :: ensstr
  CHARACTER(len=512) :: dataset_path     ! pdaf path, load from environment variable
  CHARACTER(len=512) :: obs_file     ! pdaf path, load from environment variable
  INTEGER(KIND=MPI_OFFSET_KIND) :: disp
  integer, dimension(npes_filter) :: dim_obs_all
  logical :: obs_file_exists
! ****************************************
! *** Initialize observation dimension ***
! ****************************************

! compute total number of observations
  dim_obs = obs_prcnt * dim_state
  if (dim_obs .eq. 0) then
    dim_obs = 1
  end if
  
! compute number of regions
  num_reg = dim_obs / obs_blk_size
  if ( MODULO(dim_obs, obs_blk_size) .ne. 0 ) then
    num_reg = num_reg + 1
  end if

! compute stride for regions
  stride = dim_state / num_reg
  
! determine number of obs in pe
  dim_obs_p = 0
  cnt_obs = 0
  do i=1,num_reg
    offset = (i-1) * stride
    do j=1,obs_blk_size
      index_tmp = offset + j
      if ( (index_tmp .ge. state_min_p) .and. (index_tmp .le. state_max_p) ) then
        dim_obs_p = dim_obs_p + 1
      end if
      cnt_obs = cnt_obs + 1
      if ( cnt_obs .eq. dim_obs ) exit
      if ( index_tmp .eq. state_max_p ) exit
    end do
    if ( cnt_obs .eq. dim_obs ) exit
    if ( index_tmp .eq. state_max_p ) exit
  end do
  
  IF (ALLOCATED(obs_index_p)) DEALLOCATE(obs_index_p)
  ALLOCATE( obs_index_p(dim_obs_p) )
  
! assign indices to index array
  cnt_obs_p = 0
  do i=1,num_reg
    offset = (i-1) * stride
    do j=1,obs_blk_size
      index_tmp = offset + j
      if ( (index_tmp .ge. state_min_p) .and. (index_tmp .le. state_max_p) ) then
        cnt_obs_p = cnt_obs_p + 1
        obs_index_p(cnt_obs_p) = index_tmp - (state_min_p - 1)
      end if
      if ( cnt_obs_p .eq. dim_obs_p ) exit
    end do
    if ( cnt_obs_p .eq. dim_obs_p ) exit
  end do
    
  call mpi_allgather( dim_obs_p, 1, MPI_INTEGER, dim_obs_all, &
    1, MPI_INTEGER, COMM_filter, ierr )

  disp = 0

  do i = 1, mype_filter
    disp = disp + dim_obs_all(i) * sizeof( obs_p(1) )
  end do
  
  IF (ALLOCATED(obs_p)) DEALLOCATE(obs_p)
  ALLOCATE(obs_p(dim_obs_p))
 
  write(ensstr,"(i5.5)") step
  
  call get_environment_variable( 'DATASET_PATH', dataset_path )
  
  obs_file = TRIM(dataset_path)//'/obs-'//TRIM(ensstr)//'.dat'

  inquire(file=trim(obs_file), exist=obs_file_exists)

  if( .not. obs_file_exists ) then
     WRITE (*,"(4x,a,4x,a,a,a)") &
          'ERROR', &
          'observation file not found (', trim(obs_file), ')'
     CALL mpi_abort(COMM_world, MPIerr)
  end if

  call MPI_FILE_OPEN(COMM_filter, &
                   obs_file, &
                   MPI_MODE_RDONLY, & 
                   MPI_INFO_NULL, file_id, ierr) 
  call MPI_FILE_SET_VIEW(file_id, disp , MPI_DOUBLE_PRECISION, & 
                       MPI_DOUBLE_PRECISION, 'native', & 
                       MPI_INFO_NULL, ierr) 
  call MPI_FILE_READ(file_id, obs_p, dim_obs_p, MPI_DOUBLE_PRECISION, & 
                    MPI_STATUS_IGNORE, ierr) 
  call MPI_FILE_CLOSE(file_id, ierr)  
   
END SUBROUTINE init_dim_obs_pdaf

