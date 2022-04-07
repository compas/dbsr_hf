!=======================================================================
      Module df_orbitals
!=======================================================================
!     contains description of atomic orbitals
!-----------------------------------------------------------------------
      Implicit none

      Integer :: nbf = 0                  ! number of orbitals
      Integer, allocatable :: nbs(:)      ! n-values
      Integer, allocatable :: kbs(:)      ! kappa-values
      Integer, allocatable :: lbs(:)      ! l-value
      Integer, allocatable :: jbs(:)      ! j-values
      Integer, allocatable :: ibs(:)      ! set numbers, not used here
      Integer, allocatable :: mbs(:)      ! number of splines
      Integer, allocatable :: iord(:)     ! order
      Real(8), allocatable :: qsum(:)     ! ocupation
      Real(8), allocatable :: dpm(:)      ! dicrepancy
      Logical, allocatable :: clsd(:)     ! closed shells
      Real(8), allocatable :: e(:,:)      ! orbital orthogonality
      Real(8), allocatable :: p(:,:,:)    ! radial functions
      Integer, allocatable :: iprm(:,:)   ! boundary conditions
      Character(5), allocatable :: ebs(:) ! spectroscopic notation

      End Module df_orbitals

!=======================================================================
      Subroutine alloc_df_orbitals(n)
!=======================================================================
!     allocates, deallocates or reallocate arrays in module DF_orbiyals
!-----------------------------------------------------------------------
      Use df_orbitals
      Implicit none
      Integer, intent(in) :: n
      if(allocated(nbs)) &
        Deallocate(nbs,kbs,lbs,jbs,ibs,mbs,ebs,e,iord,qsum,dpm,clsd)
        nbf=0
      if(n.le.0) return
      nbf=n
      Allocate(nbs(nbf),kbs(nbf),ibs(nbf),mbs(nbf),lbs(nbf),jbs(nbf), &
               ebs(nbf),e(nbf,nbf),iord(nbf),qsum(nbf),dpm(nbf),clsd(nbf))
      nbs=0; kbs=0; ibs=0; mbs=0; lbs=0; jbs=0; iord=0
      qsum=0.d0; dpm=0.d0; clsd = .false.

      End Subroutine alloc_DF_orbitals

!=======================================================================
      Subroutine alloc_DF_radial(ns)
!=======================================================================
!     allocates or deallocates radial one-electron functions
!-----------------------------------------------------------------------
      Use DBS_grid, only: ks
      Use df_orbitals
      Implicit none
      Integer :: ns
      if(allocated(p)) Deallocate(p,iprm)
      if(ns.le.0) Return
      Allocate(p(ns,2,nbf),iprm(ns+ns,nbf))

      End Subroutine alloc_DF_radial

!=======================================================================
      Integer Function Ifind_orb(n,k,iset)
!=======================================================================
!     find orbital (n,k,iset) in the list DF_orbitals
!----------------------------------------------------------------------
      Use df_orbitals
      Implicit none
      Integer, intent(in) :: n,k,iset
      Integer :: i
      Ifind_orb=0
      Do i=1,nbf
       if(n.ne.nbs(i)) Cycle
       if(k.ne.kbs(i)) Cycle
       if(iset.ne.ibs(i)) Cycle
       Ifind_orb = i
       Return
      End do
      End Function Ifind_orb

!=======================================================================
      Subroutine Get_pv(i,v)
!=======================================================================
!     get two-component function "i"  as one vector "v"
!-----------------------------------------------------------------------
      Use DBS_grid
      Use df_orbitals
      Implicit none
      Integer, intent(in) :: i
      Real(8), intent(out) :: v(ms)
      v(1:ns)=p(1:ns,1,i)
      v(ns+1:ms)=p(1:ns,2,i)
      End Subroutine Get_pv

!=======================================================================
      Subroutine Put_pv(i,v)
!=======================================================================
!     put two-component vector "v" in p-array, location "i"
!-----------------------------------------------------------------------
      Use DBS_grid
      Use DF_orbitals
      Implicit none
      Integer, intent(in) :: i
      Real(8), intent(in) :: v(ms)
      p(1:ns,1,i)=v(1:ns)
      p(1:ns,2,i)=v(ns+1:ms)
      End Subroutine Put_pv
