!$Id: mod_model.F90 1411 2013-09-25 14:04:41Z lnerger $
!BOP
!
! !MODULE:
MODULE mod_model

! !DESCRIPTION:
! This module provides variables needed for the
! 2-dimensional tutorial model without parallelization.
!
! !REVISION HISTORY:
! 2013-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
use iso_c_binding
  IMPLICIT NONE
  SAVE
!EOP


! *** Variables specific for 2D tutorial model ***

  INTEGER :: total_steps

END MODULE mod_model

