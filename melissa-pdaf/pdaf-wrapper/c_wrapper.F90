module my_state_accessors
use iso_c_binding
implicit none
save
real(C_DOUBLE), POINTER :: distribute_state_to(:)
real(C_DOUBLE), POINTER :: collect_state_from(:)
integer(C_INT), POINTER :: index_map(:)
integer(C_INT), POINTER :: index_map_hidden(:)
integer :: current_step
end module


! TODO: take the dummy model or maybe even others to init parallel!
SUBROUTINE cwrapper_init_pdaf(param_dim_state, param_dim_state_p, param_ensemble_size, &
        param_comm_world, dim_index_map, param_index_map, &
dim_index_map_hidden, param_index_map_hidden) BIND(C,name='cwrapper_init_pdaf')
  USE iso_c_binding

  USE mod_assimilation, &
    ONLY: dim_state_p, dim_state, dim_ens, screen

  USE mod_parallel_pdaf, &
      ONLY: COMM_world, npes_filter

  !USE mod_model, &
  !    ONLY: init_nxy

  USE my_state_accessors, &
      ONLY: current_step, index_map, index_map_hidden
  IMPLICIT NONE

  INTEGER(kind=C_INT), intent(in) :: param_dim_state     ! Global state dimension
  INTEGER(kind=C_INT), intent(in) :: param_dim_state_p   ! Local state dimension
  INTEGER(kind=C_INT), intent(in) :: param_ensemble_size ! Ensemble size
  INTEGER(kind=C_INT), intent(in) :: param_comm_world    ! World communicator as given by the melissa_server
  INTEGER(kind=C_INT), INTENT(in) :: dim_index_map                   ! PE-local state dimension
  TYPE(C_PTR), INTENT(in), VALUE :: param_index_map
  INTEGER(kind=C_INT), INTENT(in) :: dim_index_map_hidden                   ! PE-local state dimension
  TYPE(C_PTR), INTENT(in), VALUE :: param_index_map_hidden

  print *, "Initing index_map s", dim_index_map, dim_index_map_hidden
  CALL C_F_POINTER( param_index_map, index_map,[dim_index_map])
  CALL C_F_POINTER( param_index_map_hidden, index_map_hidden,[dim_index_map_hidden])

  COMM_world = param_comm_world

  ! *** Define state dimension ***
  dim_state   = param_dim_state    ! Global state dimension
  dim_state_p = param_dim_state_p  ! Local state dimension

  dim_ens = param_ensemble_size

  current_step = 0


  ! Revise parallelization for ensemble assimilation
  screen = 3  ! lots of logs
  !                       dim_ens, 0 if initialized later.
  !                          log level (0 - no logs, 3 - lots of logs)
  CALL init_parallel_pdaf(0, screen)


  !call init_nxy(npes_filter)

  ! TODO: also init parallel
END SUBROUTINE

SUBROUTINE cwrapper_init_user(param_total_steps) BIND(C,name='cwrapper_init_user')
  USE iso_c_binding

  USE mod_model, &
    ONLY: total_steps


  IMPLICIT NONE

  INTEGER(kind=C_INT), intent(in) :: param_total_steps     ! total steps




  total_steps = param_total_steps


! *** Model specifications ***




! TODO: dim_state and the other things that require initialize must be parameters!,  see used variables in initialize.f90
  ! Initialize PDAF  ! TODO: dirty to call this here but init_pdaf depends on init_ens and init ens needs nx, ny and nx_p....
  CALL init_pdaf()

  ! deterministic repeatable experiments
  call srand(42)

END SUBROUTINE

SUBROUTINE cwrapper_init_ens_hidden(dim_p, dim_ens, member_id, hidden_state_p) &
        BIND(C,name='cwrapper_init_ens_hidden')
    USE iso_c_binding
    IMPLICIT NONE
    INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
    INTEGER, INTENT(in) :: dim_ens                 ! Size of ensemble
    INTEGER, INTENT(in) :: member_id
    REAL(C_DOUBLE), INTENT(out)   :: hidden_state_p(dim_p)            ! PE-local state ensemble of member_id

    EXTERNAL :: init_ens_hidden

    CALL init_ens_hidden(dim_p, dim_ens, member_id, hidden_state_p)
END SUBROUTINE


SUBROUTINE cwrapper_PDAF_deallocate() BIND(C,name='cwrapper_PDAF_deallocate')
  USE iso_c_binding

  IMPLICIT NONE

  CALL finalize_pdaf()

END SUBROUTINE




! called by get_state
SUBROUTINE my_distribute_state(dim_p, state_p)
  USE my_state_accessors, &
    only: distribute_state_to

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  ! local State vector

  !integer :: i


  !write(*,*) '--------- 5 cells distributing state: -----------'
  !do i=1, 5
        !write (*,*) state_p(i)
  !end do
  !write(*,*) '--------- end distributing state: -----------'

  distribute_state_to(:) = state_p(:)

END SUBROUTINE my_distribute_state

! called by put_state
SUBROUTINE my_collect_state(dim_p, state_p)

  USE my_state_accessors, &
       ONLY: collect_state_from

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  ! local State vector


  !integer :: i


  state_p(:) = collect_state_from(:)
  !write(*,*) '--------- 5 cells collect_from state: -----------'
  !do i=1, 5
        !write (*,*) state_p(i)
  !end do
  !write(*,*) '--------- end collect state: -----------'

END SUBROUTINE my_collect_state


FUNCTION cwrapper_PDAF_get_state(doexit, dim_state_analysis, state_analysis, status) &
    BIND(C, name='cwrapper_PDAF_get_state')

  use iso_c_binding
  USE my_state_accessors, &
  only: distribute_state_to


  IMPLICIT NONE

! Arguments:
  INTEGER(C_INT), intent(out) :: doexit
  INTEGER(C_INT), intent(in) :: dim_state_analysis
  TYPE(C_PTR), VALUE :: state_analysis
  INTEGER(C_INT), intent(out) :: status

  INTEGER(C_INT) :: cwrapper_PDAF_get_state

  ! External subroutines
  EXTERNAL :: my_collect_state, &   ! Routine to collect a state vector from model fields
       init_dim_obs_pdaf, &         ! Initialize Dimension Of Observation Vector
       obs_op_pdaf, &               ! Implementation of the Observation operator
       init_obs_pdaf, &             ! Routine to provide vector of measurements
       prepoststep_ens_pdaf, &      ! User supplied pre/poststep routine
       prodRinvA_pdaf, &            ! Provide product R^-1 A for some matrix A
       init_obsvar_pdaf, &          ! Initialize mean observation error variance
       next_observation_pdaf, &     ! Provide time step, model time, &
                                    ! and dimension of next observation
       my_distribute_state          ! Routine to distribute a state vector to model fields
  EXTERNAL :: init_n_domains_pdaf, &   ! Provide number of local analysis domains
       init_dim_l_pdaf, &              ! Initialize state dimension for local ana. domain
       init_dim_obs_l_pdaf,&           ! Initialize dim. of obs. vector for local ana. domain
       g2l_state_pdaf, &               ! Get state on local ana. domain from global state
       l2g_state_pdaf, &               ! Init global state from state on local analysis domain
       g2l_obs_pdaf, &                 ! Restrict a global obs. vector to local analysis domain
       init_obs_l_pdaf, &              ! Provide vector of measurements for local ana. domain
       prodRinvA_l_pdaf, &             ! Provide product R^-1 A for some local matrix A
       init_obsvar_l_pdaf, &           ! Initialize local mean observation error variance
       init_obs_f_pdaf, &              ! Provide full vector of measurements for PE-local domain
       obs_op_f_pdaf, &                ! Obs. operator for full obs. vector for PE-local domain
       init_dim_obs_f_pdaf             ! Get dimension of full obs. vector for PE-local domain


! local variables
  INTEGER :: nsteps    ! Number of time steps to be performed in current forecast
  REAL :: timenow      ! Current model time

  CALL C_F_POINTER( state_analysis, distribute_state_to,[dim_state_analysis])
  ! TODO: parameters depend on filtertype
  ! Ensemble-based filters (SEIK, EnKF, ETKF, LSEIK, LETKF, ESTKF, LESTKF)
  ! This call is for all ensemble-based filters, except for EKTF/LETKF
  CALL PDAF_get_state(nsteps, timenow, doexit, next_observation_pdaf, &
       my_distribute_state, prepoststep_ens_pdaf, status)

  cwrapper_PDAF_get_state = nsteps

END FUNCTION

SUBROUTINE cwrapper_PDAF_put_state(dim_state_background, state_background, status) &
  BIND(C, name='cwrapper_PDAF_put_state')

  USE iso_c_binding

  USE my_state_accessors, &
       ONLY: collect_state_from

  USE mod_assimilation, &      ! Variables for assimilation
       ONLY: filtertype

  IMPLICIT NONE

! Arguments
  INTEGER(C_INT) :: dim_state_background
  TYPE(C_PTR), VALUE :: state_background
  INTEGER(C_INT), intent(out) :: status

  ! External subroutines
  EXTERNAL :: my_collect_state, &   ! Routine to collect a state vector from model fields
       init_dim_obs_pdaf, &         ! Initialize Dimension Of Observation Vector
       obs_op_pdaf, &               ! Implementation of the Observation operator
       init_obs_pdaf, &             ! Routine to provide vector of measurements
       prepoststep_ens_pdaf, &      ! User supplied pre/poststep routine
       prodRinvA_pdaf, &            ! Provide product R^-1 A for some matrix A
       init_obsvar_pdaf, &          ! Initialize mean observation error variance
       next_observation_pdaf, &     ! Provide time step, model time, &
                                    ! and dimension of next observation
       my_distribute_state          ! Routine to distribute a state vector to model fields
  EXTERNAL :: init_n_domains_pdaf, &   ! Provide number of local analysis domains
       init_dim_l_pdaf, &              ! Initialize state dimension for local ana. domain
       init_dim_obs_l_pdaf,&           ! Initialize dim. of obs. vector for local ana. domain
       g2l_state_pdaf, &               ! Get state on local ana. domain from global state
       l2g_state_pdaf, &               ! Init global state from state on local analysis domain
       g2l_obs_pdaf, &                 ! Restrict a global obs. vector to local analysis domain
       init_obs_l_pdaf, &              ! Provide vector of measurements for local ana. domain
       prodRinvA_l_pdaf, &             ! Provide product R^-1 A for some local matrix A
       init_obsvar_l_pdaf, &           ! Initialize local mean observation error variance
       init_obs_f_pdaf, &              ! Provide full vector of measurements for PE-local domain
       obs_op_f_pdaf, &                ! Obs. operator for full obs. vector for PE-local domain
       init_dim_obs_f_pdaf             ! Get dimension of full obs. vector for PE-local domain

   ! EnKF:
   EXTERNAL :: add_obs_error_pdaf, init_obscovar_pdaf

  !collect_state_from => state_background
  CALL C_F_POINTER( state_background, collect_state_from,[dim_state_background])

  IF (filtertype == 6) THEN
      CALL PDAF_put_state_estkf(my_collect_state, init_dim_obs_pdaf, obs_op_pdaf, &
          init_obs_pdaf, prepoststep_ens_pdaf, prodRinvA_pdaf, init_obsvar_pdaf, status)
  ELSE IF (filtertype == 2) THEN
      CALL PDAF_put_state_enkf(my_collect_state, init_dim_obs_pdaf, obs_op_pdaf, &
          init_obs_pdaf, prepoststep_ens_pdaf, add_obs_error_pdaf, init_obscovar_pdaf, &
          status)
  END IF



  !IF (filtertype == 6) THEN
     !CALL PDAF_assimilate_estkf(collect_state_pdaf, distribute_state_pdaf, &
          !init_dim_obs_pdaf, obs_op_pdaf, init_obs_pdaf, prepoststep_ens_pdaf, &
          !prodRinvA_pdaf, init_obsvar_pdaf, next_observation_pdaf, status_pdaf)
  !ELSEIF (filtertype == 7) THEN
     !CALL PDAF_assimilate_lestkf(collect_state_pdaf, distribute_state_pdaf, &
          !init_dim_obs_f_pdaf, obs_op_f_pdaf, init_obs_f_pdaf, init_obs_l_pdaf, &
          !prepoststep_ens_pdaf, prodRinvA_l_pdaf, init_n_domains_pdaf, &
          !init_dim_l_pdaf, init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
          !g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, next_observation_pdaf, status_pdaf)
  !END IF

END SUBROUTINE

SUBROUTINE cwrapper_set_current_step(new_current_step) &
  BIND(C, name='cwrapper_set_current_step')

  USE iso_c_binding

  USE my_state_accessors, &
       ONLY: current_step

  IMPLICIT NONE

! Arguments
  INTEGER(C_INT) :: new_current_step

  write(*,"(4x,a,4x,a,i10)") "[DEBUG]",&
    "CURRENT_STEP IN C_WRAPPER: ", current_step

  ! PDAF starts counting at 0, the assimilator starts at one. TODO: do the same counting in melissa and pdaf?
  current_step = new_current_step

END SUBROUTINE

