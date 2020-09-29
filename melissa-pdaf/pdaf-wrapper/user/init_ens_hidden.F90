!BOP
!
! !ROUTINE: init_ens_hidden --- Initialize ensemble members hidden state, state that is not assimilated but needed to restart
! model...
!
! !INTERFACE:
SUBROUTINE init_ens_hidden(dim_p, dim_ens, member_id, hidden_state_p)
    ! !DESCRIPTION:
    ! User-supplied routine for Melissa-DA with PDAF.
    ! Used in the filters: SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
    !
    ! The routine is called when the filter is
    ! initialized in the constructor of the PDAF assimilator.  It has
    ! to initialize an ensemble of dim\_ens states.
    ! Typically, the ensemble will be directly read from files.
    !
    ! The routine is called by all filter processes and
    ! initializes the ensemble for the PE-local domain.
    !
    ! !REVISION HISTORY:
    ! 2013-02 - Lars Nerger - Initial code
    ! Later revisions - see svn log
    !
    ! !USES:

    IMPLICIT NONE

    ! !ARGUMENTS:
    INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
    INTEGER, INTENT(in) :: dim_ens                 ! Size of ensemble
    INTEGER, INTENT(in) :: member_id               ! id of the member that shall be initialized by this function. starts at 0 and
    ! goes up to dim_ens - 1
    REAL, INTENT(out)   :: hidden_state_p(dim_p)            ! PE-local state ensemble of member_id




! dummy, not needed for lorenz!

END SUBROUTINE init_ens_hidden
