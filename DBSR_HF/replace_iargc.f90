!======================================================================
      Integer Function Iargc()
!======================================================================
!     replace iargc from FORTRAN-90
!     (this routine is absent in some last versions of gfortran)
!----------------------------------------------------------------------
      Implicit none
      iargc = command_argument_count()
      End  Function Iargc

