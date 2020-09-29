!$Id: mod_parallel_model.F90 1411 2013-09-25 14:04:41Z lnerger $
!BOP
!
! !MODULE:
MODULE mod_parallel_model

! !DESCRIPTION:
! This modules provides variables for the MPI parallelization
! of the tutorial model to be shared between model-related routines.
!
! In addition, methods to initialize and finalize MPI are provided.
!
! !REVISION HISTORY:
! 2004-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
USE mod_parallel_pdaf, &
    ONLY: COMM_world
  IMPLICIT NONE
  SAVE

  INCLUDE 'mpif.h'

! !PUBLIC DATA MEMBERS:
  ! Basic variables for model state integrations
  INTEGER :: COMM_model  ! MPI communicator for model tasks
  INTEGER :: mype_model  ! Number of PEs in COMM_model
  INTEGER :: npes_model  ! PE rank in COMM_model
  INTEGER :: mype_world  ! Number of PEs in COMM_world
  INTEGER :: npes_world  ! PE rank in COMM_world
  INTEGER :: MPIerr      ! Error flag for MPI
!EOP

CONTAINS
!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_parallel - Initialize MPI
!
! !INTERFACE:
  SUBROUTINE init_parallel()

! !DESCRIPTION:
! Routine to initialize MPI, the number of PEs
! (npes\_world) and the rank of a PE (mype\_world).
! The model is executed within the scope of the
! communicator Comm_model. It is also initialized
! here together with its size (npes\_model) and
! the rank of a PE (mype\_model) within Comm_model.
!EOP

    IMPLICIT NONE

    INTEGER :: i  ! error flag

    ! Reuse mpi inited by server

    CALL MPI_Comm_Size(COMM_world,npes_world,i)
    CALL MPI_Comm_Rank(COMM_world,mype_world,i)

    ! Initialize model communicator, its size and the process rank
    ! Here the same as for COMM_world
    COMM_model = COMM_world   !pdaf will think there is one model. do w e really need this?....
    npes_model = npes_world
    mype_model = mype_world

  END SUBROUTINE init_parallel
!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: finalize_parallel - Finalize MPI
!
! !INTERFACE:
  SUBROUTINE finalize_parallel()

! !DESCRIPTION:
! Routine to finalize MPI
!EOP

    IMPLICIT NONE

    CALL  MPI_Barrier(COMM_world, MPIerr)
    ! Server will call MPI_Finalize

  END SUBROUTINE finalize_parallel
!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: abort_parallel - Abort MPI
!
! !INTERFACE:
  SUBROUTINE abort_parallel()

! !DESCRIPTION:
! Routine to abort MPI program
!EOP

    IMPLICIT NONE

    CALL  MPI_Abort(COMM_world, 1, MPIerr)

  END SUBROUTINE abort_parallel

END MODULE mod_parallel_model
