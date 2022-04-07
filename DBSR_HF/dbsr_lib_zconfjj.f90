!--------------------------------------------------------------
!     ZCONFJJ  LIBRARY
!--------------------------------------------------------------
!     contains  common routines dealing with configuration and
!     spectroscopic notation in jj-coupling
!--------------------------------------------------------------
!     MODULES:
!
!     mod_orb_jj.f90    -  list of one-electron orbitals
!     mod_symc.f90      -  description of configuration
!     mod_symt.f90      -  description of shell terms
!     mod_conf_jj.f90   -  description of atomic states
!     mod_boef_list.f90 -  matrix elements in uncouple nlmj-representation
!     mod_det_jq.f90    -  determinant shell's expentions
!--------------------------------------------------------------
!     ROUTINES:
!
!     cfp_jj.f90
!     coef_1conf.f90
!     EL_nljk.f90
!     Incode_cj.f90
!     Jdef_nc.f90
!     jterm.f90
!     mj_value.f90
!--------------------------------------------------------------



!=====================================================================
      Module orb_jj
!=====================================================================
!     defines the one-electron orbitals in jj-coupling
!---------------------------------------------------------------------
      Implicit none

! ... list of one-electron orbitals

      Integer :: mwf = 0        ! max. number of orbitals
      Integer :: nwf = 0        ! current number of orbitals
      Integer :: iwf = 2**10    ! initial prediction for mwf

      Integer, allocatable :: nef(:)   ! n-values
      Integer, allocatable :: kef(:)   ! kappa-value
      Integer, allocatable :: lef(:)   ! l-value
      Integer, allocatable :: jef(:)   ! j-value
      Integer, allocatable :: ief(:)   ! subset index
      Integer, allocatable :: ipef(:)  ! additional pointer

      Character(5), allocatable :: ELF(:) ! spectroscopic notation

! ... orbital orthogonality and AFTER conditions

      Integer :: JORT = 1
      Integer, allocatable :: IORT(:,:)

! ... possible set indexes for orthogonal subsets:

      Integer, parameter :: kset = 61
      Character(kset) :: ASET = &
       '123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
      Integer, parameter :: ksmax = kset*kset

! ... shift for the set indexes:

      Integer :: kshift = 0

      End Module orb_jj


!======================================================================
      Subroutine alloc_orb_jj(m)
!======================================================================
!     allocate, deallocate or reallocate orbital list in module orb_jj
!----------------------------------------------------------------------
      Use orb_jj

      Implicit none
      Integer, intent(in) :: m
      Integer :: i,j
      Integer, allocatable :: iarr(:)
      Integer, allocatable :: jarr(:,:)
      Character(5), allocatable :: carr(:)

      if(m.le.0) then
       if(allocated(nef)) Deallocate (nef,kef,lef,jef,ief,ipef,ELF,IORT)
       mwf = 0; nwf = 0
      elseif(.not.allocated(nef)) then
       mwf = m
       Allocate(nef(mwf),kef(mwf),lef(mwf),jef(mwf),ief(mwf),ipef(mwf), &
                ELF(mwf), IORT(mwf,mwf))
       nef=0; kef=0; lef=0; jef=0; ief=0; ipef=0; IORT=0
      elseif(m.le.mwf) then
       Return
      elseif(nwf.eq.0) then
       Deallocate (nef,kef,lef,jef,ief,ipef,ELF,IORT)
       mwf = m
       Allocate(nef(mwf),kef(mwf),lef(mwf),jef(mwf),ief(mwf),ipef(mwf), &
                ELF(mwf), IORT(mwf,mwf))
       nef=0; kef=0; lef=0; jef=0; ief=0; ipef=0; IORT=0
      else
       mwf=m
       Allocate(iarr(nwf))
       iarr(1:nwf)=nef(1:nwf); Deallocate(nef); Allocate(nef(m)); nef=0
       nef(1:nwf)=iarr(1:nwf)
       iarr(1:nwf)=kef(1:nwf); Deallocate(kef); Allocate(kef(m)); kef=0
       kef(1:nwf)=iarr(1:nwf)
       iarr(1:nwf)=lef(1:nwf); Deallocate(lef); Allocate(lef(m)); lef=0
       lef(1:nwf)=iarr(1:nwf)
       iarr(1:nwf)=jef(1:nwf); Deallocate(jef); Allocate(jef(m)); jef=0
       jef(1:nwf)=iarr(1:nwf)
       iarr(1:nwf)=ief(1:nwf); Deallocate(ief); Allocate(ief(m)); ief=0
       ief(1:nwf)=iarr(1:nwf)
       iarr(1:nwf)=ipef(1:nwf); Deallocate(ipef); Allocate(ipef(m)); ipef=0
       ipef(1:nwf)=iarr(1:nwf)
       Deallocate(iarr)

       Allocate(jarr(nwf,nwf))
       jarr(1:nwf,1:nwf)=IORT(1:nwf,1:nwf); Deallocate(IORT)
       Allocate(IORT(m,m)); IORT(1:nwf,1:nwf)=jarr(1:nwf,1:nwf)
       Deallocate(jarr)

       Allocate(carr(nwf))
       carr(1:nwf)=ELF(1:nwf); Deallocate(ELF); Allocate(ELF(m))
       ELF(1:nwf)=carr(1:nwf)
       Deallocate(carr)

      end if

      End Subroutine alloc_orb_jj

!=======================================================================
      Real(8) Function memory_orb_jj(m)
!=======================================================================
!     return requred memory for module orb_jj
!------------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: m
      memory_orb_jj = (6*4*m + 5*m + 4*m*m)/(1024d0*1024d0)
      End Function memory_orb_jj


!=======================================================================
      Integer Function Ifind_jjorb(n,k,ii,job)
!=======================================================================
!     find orbital (n,k,iset) in the list "orb_jj"
!     Options:
!     job = 0  -  no further actions
!     job = 1  -  stop if fail to find
!     job = 2  -  add new orbital
!------------------------------------------------------------------------
      USE orb_jj

      Implicit none
      Integer :: n,k,i,ii,j,job
      Character(5), external :: ELi

      Ifind_jjorb=0; i = ii + kshift

      Do j=1,nwf
       if(n.ne.nef(j)) Cycle
       if(k.ne.kef(j)) Cycle
       if(i.ge.0.and.i.ne.ief(j)) Cycle
       Ifind_jjorb=j
       Return
      End do
      if(job.eq.0) Return

      if(job.eq.1) then
       Write(*,'(a,a,5i5)') 'Ifind_jjorb can not find the orbital:',&
                            ' N,K,iset = ',n,k,i
       Stop
      end if

      if(nwf+1.gt.mwf) Call Alloc_orb_jj(mwf+iwf)
      nwf = nwf + 1
      nef(nwf)=n
      kef(nwf)=k
      lef(nwf) = k; if(k.lt.0) lef(nwf)=-k-1
      jef(nwf) = k+k-1; if(k.lt.0) jef(nwf)=-k-k-1
      ief(nwf)=i
      ELF(nwf)=ELi(n,k,i)
      Ifind_jjorb = nwf

      End Function Ifind_jjorb


!======================================================================
      Subroutine Read_bsw_orb_jj(nu)
!======================================================================
!     read only spectroscopic notation from bsw-file (unit nu)
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu
      Integer :: i,j,k,l,n
      Character(5) :: elw
      Integer, external :: Ifind_jjorb
      Real(8) :: p

      rewind(nu)
      read(nu) i ! ignore the rest of the record
    1 read(nu,end=2) elw
      Call EL_NLJK(elw,n,k,l,j,i)
      i = Ifind_jjorb(n,k,i,2)
      read(nu) p
      read(nu) p
      go to 1
    2 Close(nu)

      End Subroutine Read_bsw_orb_jj


!=======================================================================
      Subroutine Get_orb_jj(i,n,k,iset,ip)
!=======================================================================
!     get parameters for orbital "i" (for external calls)
!------------------------------------------------------------------------
      Use orb_jj

      Implicit none
      Integer, intent(in) :: i
      Integer, intent(out) :: n,k,iset,ip

      n = 0; k = 0;  iset=0; ip =0
      if(i.le.0.or.i.gt.nwf) Return
      n = nef(i)
      k = kef(i)
      iset = ief(i)
      ip = ipef(i)

      End Subroutine Get_orb_jj



!======================================================================
      MODULE symc_list
!======================================================================
!     Containes the "pure-configuration" list  (i.e., without terms)
!----------------------------------------------------------------------
      Implicit none

      Integer :: nsymc = 0       ! number of symmmetries
      Integer :: msymc = 0       ! current dimension of the list
      Integer :: isymc = 10000   ! initial dimension
      Integer :: jsymc = 10      ! average size of one symc.
      Integer :: ksymc = 0       ! max.dimension of all symc.s
      Integer :: lsymc = 0       ! last element

      Integer(1), allocatable :: JT_conf(:)
      Integer(1), allocatable :: no_conf(:)
      Integer,    allocatable :: ip_conf(:)
      Integer(1), allocatable :: iq_conf(:)
      Integer(1), allocatable :: kn_conf(:)

! ... IC_term(ic) - gives for given configuration 'ic' the range of corr.
!                   terms it1,it2 - first and final position in the
!                   ordered list of terms

      Integer, allocatable :: IC_term1(:)
      Integer, allocatable :: IC_term2(:)

! ... IC_need(:) - define the need of calc. between two config.s

      Integer, allocatable :: IC_need(:)

! ... JC_need(:) - define the need of calc. for the given config.

      Integer, allocatable :: JC_need(:)

      End MODULE symc_list


!======================================================================
      Subroutine alloc_symc(m)
!======================================================================
!     allocate, deallocate, or reallocate arrays in module symc_list
!----------------------------------------------------------------------
      Use symc_list

      Implicit none
      Integer, Intent(in) :: m
      Integer :: i
      Integer, allocatable :: iarr(:)

      if(m.le.0) then
       if(allocated(JT_conf)) Deallocate (JT_conf,no_conf,ip_conf, &
                                          iq_conf,kn_conf)
       msymc = 0; nsymc = 0; ksymc = 0; lsymc = 0
      elseif(.not.allocated(JT_conf)) then
       msymc = m; ksymc = msymc*jsymc; lsymc = 0
       Allocate(JT_conf(msymc),no_conf(msymc),ip_conf(msymc), &
                iq_conf(ksymc),kn_conf(ksymc))
      elseif(m.le.msymc) then
       Return
      elseif(nsymc.eq.0) then
       Deallocate (JT_conf,no_conf,ip_conf,iq_conf,kn_conf)
       msymc = m; ksymc = msymc*jsymc; lsymc = 0
       Allocate(JT_conf(msymc),no_conf(msymc),ip_conf(msymc), &
                iq_conf(ksymc),kn_conf(ksymc))
      else

       msymc=m; i=lsymc/nsymc+1; if(jsymc.lt.i) jsymc=i; ksymc=jsymc*msymc

       Allocate(iarr(lsymc))
       iarr(1:nsymc)=JT_conf(1:nsymc); Deallocate(JT_conf)
       Allocate(JT_conf(m)); JT_conf=0;  JT_conf(1:nsymc)=iarr(1:nsymc)
       iarr(1:nsymc)=no_conf(1:nsymc); Deallocate(no_conf)
       Allocate(no_conf(m)); no_conf=0;  no_conf(1:nsymc)=iarr(1:nsymc)
       iarr(1:nsymc)=ip_conf(1:nsymc); Deallocate(ip_conf)
       Allocate(ip_conf(m)); ip_conf=0;  ip_conf(1:nsymc)=iarr(1:nsymc)
       iarr(1:lsymc)=iq_conf(1:lsymc); Deallocate(iq_conf)
       Allocate(iq_conf(ksymc)); iq_conf=0;  iq_conf(1:lsymc)=iarr(1:lsymc)
       iarr(1:lsymc)=kn_conf(1:lsymc); Deallocate(kn_conf)
       Allocate(kn_conf(ksymc)); kn_conf=0;  kn_conf(1:lsymc)=iarr(1:lsymc)
       Deallocate(iarr)

      end if

      End Subroutine alloc_symc


!======================================================================
      Integer Function Iadd_symc (JT,no,iq,kn)
!======================================================================
!     add new conf.symmetry to the module symc_list
!----------------------------------------------------------------------
      Use symc_list

      Implicit none
      Integer :: no, JT, i,j,ip
      Integer :: iq(*),kn(*)

      Iadd_symc = 0
      if(no.le.0) Return
      if(msymc.le.0) Call Alloc_symc(isymc)

! ... check if the same symc. is already in the list:

      Do i=1,nsymc
       if(JT_conf(i).ne.JT) Cycle
       if(no_conf(i).ne.no) Cycle
       ip=ip_conf(i); Iadd_symc = i
       Do j=1,no; ip=ip+1
        if(iq(j).ne.iq_conf(ip)) then; Iadd_symc = 0; Exit; end if
        if(kn(j).ne.kn_conf(ip)) then; Iadd_symc = 0; Exit; end if
       End do
       if(Iadd_symc.ne.0) Return
      End do

! ... Add new symc.:

      if(nsymc.ge.msymc.or.lsymc+no.ge.ksymc) &
       Call Alloc_symc(msymc+isymc)

      nsymc=nsymc+1
      JT_conf(nsymc)=JT
      no_conf(nsymc)=no
      ip_conf(nsymc)=lsymc
      Do i=1,no; lsymc=lsymc+1
       iq_conf(lsymc)=iq(i)
       kn_conf(lsymc)=kn(i)
      End do
      Iadd_symc=nsymc

      End Function Iadd_symc


!======================================================================
      Subroutine Get_symc(ic,JT,no,nn,kn,ln,jn,iq,in)
!======================================================================
!     extract configuration 'ic' from symc_list
!----------------------------------------------------------------------
      Use symc_list

      Implicit none
      Integer :: ic,JT,no,i,ii,ip
      Integer, Dimension(*) :: nn,kn,ln,jn,iq,in
      Integer, External :: l_kappa, j_kappa

      if(ic.le.0.or.ic.gt.nsymc) Stop 'Get_symc: <ic> is out of range'

      JT = JT_conf(ic)
      no = no_conf(ic)
      ip = ip_conf(ic)
      Do i=1,no; ip = ip+1
       iq(i) = abs(iq_conf(ip))
       kn(i) = kn_conf(ip)
       ln(i) = l_kappa(kn(i))
       jn(i) = j_kappa(kn(i))
       nn(i) = i
       if(iq_conf(ip).lt.0) nn(i)=i-1
       in(i) = 0
      End do

      End Subroutine Get_symc


!======================================================================
      Integer Function Get_no(iconf)
!======================================================================
!     number of shells in configuration 'iconf'
!----------------------------------------------------------------------
      Use symc_list
      Implicit none
      Integer :: iconf
      Get_no = no_conf(iconf)
      End Function Get_no


!======================================================================
      Subroutine read_symc(nu)
!======================================================================
!     read symc_list from unit "nu"
!----------------------------------------------------------------------
      Use symc_list

      Implicit none
      Integer :: nu,i

      read(nu) nsymc,lsymc
      if(allocated(JT_conf)) Deallocate (JT_conf,no_conf,ip_conf, &
                                         iq_conf,kn_conf)
      Allocate(JT_conf(nsymc),no_conf(nsymc),ip_conf(nsymc), &
               iq_conf(lsymc),kn_conf(lsymc))
      read(nu) (JT_conf(i),i=1,nsymc)
      read(nu) (no_conf(i),i=1,nsymc)
      read(nu) (ip_conf(i),i=1,nsymc)
      read(nu) (iq_conf(i),i=1,lsymc)
      read(nu) (kn_conf(i),i=1,lsymc)
      msymc = nsymc; jsymc = lsymc/nsymc + 1; ksymc = lsymc

      End Subroutine read_symc


!======================================================================
      Subroutine dummy_symc(nu)
!======================================================================
!     dummy read symc_list from unit "nu"
!----------------------------------------------------------------------
      Use symc_list

      Implicit none
      Integer :: nu,i,j,n,l
      Integer(1) :: k

      read(nu) n,l
      read(nu) (k,i=1,n)
      read(nu) (k,i=1,n)
      read(nu) (j,i=1,n)
      read(nu) (k,i=1,l)
      read(nu) (k,i=1,l)

      End Subroutine dummy_symc


!======================================================================
      Subroutine write_symc(nu)
!======================================================================
!     write symc_list to unit "nu"
!----------------------------------------------------------------------
      Use symc_list

      Implicit none
      Integer :: nu,i

      write(nu) nsymc,lsymc
      write(nu) (JT_conf(i),i=1,nsymc)
      write(nu) (no_conf(i),i=1,nsymc)
      write(nu) (ip_conf(i),i=1,nsymc)
      write(nu) (iq_conf(i),i=1,lsymc)
      write(nu) (kn_conf(i),i=1,lsymc)
      Deallocate (JT_conf,no_conf,ip_conf,iq_conf,kn_conf)
      msymc = 0; jsymc = lsymc/nsymc + 1; ksymc = 0; lsymc=0

      End Subroutine write_symc


!======================================================================
      MODULE symt_list
!======================================================================
!     Containes the configuration angular symmetries (i.e. all terms)
!----------------------------------------------------------------------
      Implicit none

      Integer :: nsymt = 0       ! number of symmetries
      Integer :: msymt = 0       ! current dimension of the list
      Integer :: isymt = 10000   ! initial dimension
      Integer :: jsymt = 10      ! average size of one symmetry
      Integer :: ksymt = 0       ! dimension of all symt.s
      Integer :: lsymt = 0       ! last element

      Integer,    allocatable :: IT_conf(:)
      Integer,    allocatable :: ip_term(:)
      Integer(2), allocatable :: JS_term(:)
      Integer(1), allocatable :: VS_term(:)
      Integer(2), allocatable :: JI_term(:)

! ... IT_stat     - indicate existence the term in given case

      Integer(1), allocatable :: IT_stat(:)

! ... JP_term     - provides ordering of terms according to configurations

      Integer,    allocatable :: JP_term(:)

! ... IT_done(:) - pointer on the done calculation for specific
!                  operators and given terms

      Integer(1), allocatable :: IT_done(:)

! ... IT_need(:) - define the need of calc. for the given term.symmetry

      Integer(1), allocatable :: IT_need(:)

      Integer, parameter ::  mrecl = 100000 ! just to record IT_done

      End MODULE symt_list


!======================================================================
      Subroutine alloc_symt(m)
!======================================================================
!     allocate, deallocate or reallocate arrays in module symt_list
!----------------------------------------------------------------------
      Use symt_list

      Implicit none
      Integer, Intent(in) :: m
      Integer :: i
      Integer, Allocatable :: iarr(:)

      if(m.le.0) then
       if(allocated(IT_conf)) Deallocate (IT_conf,ip_term,&
                                          JS_term,VS_term,JI_term)
       msymt = 0; nsymt = 0; ksymt = 0; lsymt = 0
      elseif(.not.allocated(IT_conf)) then
       msymt = m; ksymt = msymt*jsymt; lsymt = 0
       Allocate(ip_term(msymt),IT_conf(msymt),&
                JS_term(ksymt),VS_term(ksymt),JI_term(ksymt))
      elseif(m.le.msymt) then
       Return
      elseif(nsymt.eq.0) then
       Deallocate (IT_conf,ip_term,JS_term,VS_term,JI_term)
       msymt = m; ksymt = msymt*jsymt; lsymt = 0
       Allocate(ip_term(msymt),IT_conf(msymt),&
                JS_term(ksymt),VS_term(ksymt),JI_term(ksymt))
      else
       msymt=m; i=lsymt/nsymt+1; if(jsymt.lt.i) jsymt=i; ksymt=jsymt*msymt
       Allocate(iarr(lsymt))
       iarr(1:nsymt)=IT_conf(1:nsymt); Deallocate(IT_conf)
       Allocate(IT_conf(m)); IT_conf=0;  IT_conf(1:nsymt)=iarr(1:nsymt)
       iarr(1:nsymt)=ip_term(1:nsymt); Deallocate(ip_term)
       Allocate(ip_term(m)); ip_term=0;  ip_term(1:nsymt)=iarr(1:nsymt)
       iarr(1:lsymt)=JS_term(1:lsymt); Deallocate(JS_term)
       Allocate(JS_term(ksymt)); JS_term=0;  JS_term(1:lsymt)=iarr(1:lsymt)
       iarr(1:lsymt)=VS_term(1:lsymt); Deallocate(VS_term)
       Allocate(VS_term(ksymt)); VS_term=0;  VS_term(1:lsymt)=iarr(1:lsymt)
       iarr(1:lsymt)=JI_term(1:lsymt); Deallocate(JI_term)
       Allocate(JI_term(ksymt)); JI_term=0;  JI_term(1:lsymt)=iarr(1:lsymt)
       Deallocate(iarr)
      end if

      End Subroutine alloc_symt


!======================================================================
      Integer Function Iadd_symt(iconf,no,Jshell,Vshell,Jintra)
!======================================================================
!     add new overlap conf.symmetry to symt_list
!----------------------------------------------------------------------
      Use symt_list

      Implicit none
      Integer :: no,iconf, i,j,ip
      Integer, Dimension(*) :: Jshell,Vshell,Jintra

      Iadd_symt = 0
      if(no.le.0) Return
      if(msymt.eq.0) Call Alloc_symt(isymt)

! ... check if the same symt. is already in the list:

      Do i=1,nsymt
       if(iabs(IT_conf(i)).ne.iconf) Cycle
       ip=ip_term(i); Iadd_symt = i
       Do j=1,no; ip = ip+1
        if(Jshell(j).ne.JS_term(ip)) then; Iadd_symt = 0; Exit; end if
        if(Vshell(j).ne.VS_term(ip)) then; Iadd_symt = 0; Exit; end if
        if(Jintra(j).ne.JI_term(ip)) then; Iadd_symt = 0; Exit; end if
       End do
       if(Iadd_symt.ne.0) Return
      End do

! ... Add new symt.:

      if(nsymt.ge.msymt.or.lsymt+no.ge.ksymt) &
         Call Alloc_symt(msymt+isymt)
      nsymt=nsymt+1
      ip_term(nsymt)=lsymt
      IT_conf(nsymt)=iconf
      Do i=1,no; lsymt=lsymt+1
       JS_term(lsymt)=Jshell(i)
       VS_term(lsymt)=Vshell(i)
       JI_term(lsymt)=Jintra(i)
      End do
      Iadd_symt=nsymt

      End Function Iadd_symt


!======================================================================
      Subroutine Get_symt(iterm,iconf,no,Jshell,Vshell,Jintra)
!======================================================================
!     extracts angular symmetry 'iterm'
!----------------------------------------------------------------------
      Use symt_list

      Implicit none
      Integer, Intent(in) :: iterm
      Integer, Intent(out) :: iconf,no
      Integer :: i,ip
      Integer, Dimension(*) :: Jshell,Vshell,Jintra
      Integer, External :: Get_no

      if(iterm.le.0.or.iterm.gt.nsymt)  then
        write(*,*) 'Get_symt: it =',iterm, &
                   ' is out of range, nsymt=',nsymt
        Stop
      end if

      iconf=it_conf(iterm)
      no = Get_no(iconf)
      ip = ip_term(iterm)
      Do i=1,no; ip=ip+1
       Jintra(i) = JI_term(ip)
       Jshell(i) = JS_term(ip)
       Vshell(i) = VS_term(ip)
      End do

      End Subroutine Get_symt


!======================================================================
      Subroutine Read_symt(nu)
!======================================================================
!     reads the list of angular symmetries from unit 'nu'
!----------------------------------------------------------------------
      Use symt_list

      Implicit none
      Integer :: nu,i

      read(nu) nsymt,lsymt
      if(Allocated(IT_conf)) Deallocate (IT_conf,ip_term, &
                                        JS_term,VS_term,JI_term)
      msymt = nsymt; jsymt=lsymt/nsymt+1; ksymt = lsymt
      Allocate(ip_term(msymt),IT_conf(msymt),&
               JS_term(ksymt),VS_term(ksymt),JI_term(ksymt))
      read(nu) (IT_conf(i),i=1,nsymt)
      read(nu) (ip_term(i),i=1,nsymt)
      read(nu) (JS_term(i),i=1,lsymt)
      read(nu) (VS_term(i),i=1,lsymt)
      read(nu) (JI_term(i),i=1,lsymt)

      End Subroutine Read_symt

!======================================================================
      Subroutine Dummy_symt(nu)
!======================================================================
!     skips the list of angular symmetries in unit 'nu'
!----------------------------------------------------------------------
      Implicit none
      Integer :: nu,i,j,n,l
      Integer(2) :: m
      Integer(1) :: k

      read(nu) n,l
      read(nu) (j,i=1,n)
      read(nu) (j,i=1,n)
      read(nu) (m,i=1,l)
      read(nu) (k,i=1,l)
      read(nu) (m,i=1,l)

      End Subroutine Dummy_symt

!======================================================================
      Subroutine Write_symt(nu)
!======================================================================
!     records the list of angular symmetries to unit 'nu'
!----------------------------------------------------------------------
      Use symt_list

      Implicit none
      Integer :: nu,i

      write(nu) nsymt,lsymt
      write(nu) (IT_conf(i),i=1,nsymt)
      write(nu) (ip_term(i),i=1,nsymt)
      write(nu) (JS_term(i),i=1,lsymt)
      write(nu) (VS_term(i),i=1,lsymt)
      write(nu) (JI_term(i),i=1,lsymt)

      Deallocate (IT_conf,ip_term,JS_term,VS_term,JI_term)
      msymt = 0; jsymt=lsymt/nsymt+1; ksymt = 0; lsymt=0

      End Subroutine Write_symt


!======================================================================
      Subroutine Write_done(nu)
!======================================================================
!     record the IT_done array to unit "nu"
!----------------------------------------------------------------------
      Use symt_list

      Implicit none
      Integer :: nu,i,i1,i2,n

      n = nsymt*(nsymt+1)/2
      Do i = 1,n
       if(IT_done(i).eq. 0) Stop 'Write_done: IT_done=0 - ?'
       if(IT_done(i).eq.-1) IT_done(i)=0
      End do

      write(nu) n
      i1=1; i2=mrecl; if(i2.gt.n) i2=n
      Do
       write(nu) (IT_done(i),i=i1,i2)
       i1=i1+mrecl; if(i1.gt.n) Exit
       i2=i2+mrecl; if(i2.gt.n) i2=n
      End do

      End Subroutine Write_done


!======================================================================
      Subroutine Read_done(nu)
!======================================================================
!     reads the IT_done array from unit "nu"
!----------------------------------------------------------------------
      Use symt_list

      Implicit none
      Integer :: nu,i,i1,i2,n

      read(nu) n
      i1=1; i2=mrecl; if(i2.gt.n) i2=n
      Do
       read(nu) (IT_done(i),i=i1,i2)
       i1=i1+mrecl; if(i1.gt.n) Exit
       i2=i2+mrecl; if(i2.gt.n) i2=n
      End do

      End Subroutine Read_done


!======================================================================
      Subroutine Dummy_done(nu)
!======================================================================
!     skips  the IT_done array in unit "nu"
!----------------------------------------------------------------------
      Use symt_list

      Implicit none
      Integer :: nu,i,i1,i2,n
      Integer(1) :: j

      read(nu) n
      i1 = 1; i2 = mrecl; if(i2.gt.n) i2=n
      Do
       read(nu) (j,i=i1,i2)
       i1=i1+mrecl; if(i1.gt.n) Exit
       i2=i2+mrecl; if(i2.gt.n) i2=n
      End do

      End Subroutine Dummy_done

!======================================================================
      Integer Function Get_iconf(iterm)
!======================================================================
!     provides "iconf" pointer for given "iterm" (for external calls)
!----------------------------------------------------------------------
      Use symt_list
      Implicit none
      Integer, intent(in) :: iterm
      Get_iconf=it_conf(iterm)
      End Function Get_iconf


!=====================================================================
      Module conf_jj
!=====================================================================
!     containes description of the configuration atomic states
!     and related information
!
!     The configuration atomic states (CAS) are defined by "iterm"
!     pointer on the "angular symmetry" and by list of specific
!     orbitals, which are placed in array ip_orb:
!     iterm     => IS_term(ic);
!     orbitals  => ip_orb(ip+1,...ip+no), where ip => IP_state(ic)
!     ip_orb(i) => pointer on orbital in common list of orbitals
!---------------------------------------------------------------------
      Implicit none

      Integer :: ne     = 0     !  number of electrons
      Integer :: parity = 0     !  parity of states (+1,-1)

! ... description of 1 conf.w.function:

      Integer, parameter :: msh = 31 ! max. number of shells behind core

      Integer :: no, Jtotal, iconf, iterm
      Integer, dimension(msh) :: nn,kn,ln,jn,iq,in,       &
                                 Jshell,Vshell,Jintra, &
                                 np_symc,np_symt, np_orb, np,mp

      Integer :: no1, Jtotal1, iconf1, iterm1
      Integer, dimension(msh) :: nn1,kn1,ln1,jn1,iq1,in1,       &
                                 Jshell1,Vshell1,Jintra1, &
                                 np_symc1,np_symt1, np_orb1, np1,mp1
      Integer :: no2, Jtotal2, iconf2, iterm2
      Integer, dimension(msh) :: nn2,kn2,ln2,jn2,iq2,in2,       &
                                 Jshell2,Vshell2,Jintra2, &
                                 np_symc2,np_symt2, np_orb2, np2,mp2

! ... storing configuration as a character string:

      Integer, parameter :: mas = 9*msh+3
      Character(mas) :: CONFIG, SHELLJ, INTRAJ, AS,BS,CS
      Integer :: ia

! ... core:

      Integer,parameter :: mcore = 50
      Character(250) :: core, closed
      Integer :: ncore = 0
      Integer :: nn_core(mcore),k_core(mcore),l_core(mcore),j_core(mcore)

! ... configuration list:

      Integer :: ncfg  = 0       !  current number of configurations
      Integer :: mcfg  = 0       !  max. number of configurations
      Integer :: icfg  = 50000   !  initial prediction of mcfg
      Integer :: jcfg  = 10      !  average number of shells
      Integer :: kcfg  = 0       !  max. dimension (mcfg*jcfg)
      Integer :: lcfg  = 0       !  last element

      Integer,    allocatable :: IP_state(:)
      Integer,    allocatable :: IS_term(:)
      Integer(2), allocatable :: IP_orb(:)

      Integer, allocatable :: IT_state1(:),IT_state2(:),IS_order(:)

      Integer :: J_min=-1, J_max=-1

! ... expansion coeficients:

      Real(8), allocatable :: WC(:)

! ... label assignments:

      Integer, parameter :: mlab = 64
      Character(mlab), allocatable :: LABEL(:)

! ... J-blocks:

      Integer :: njbl  = 0
      Integer, allocatable :: JJc (:),JTc1 (:),JTc2 (:),Jncfg (:), JTp(:)
      Integer :: njbl1 = 0
      Integer, allocatable :: JJ1c(:),JT1c1(:),JT1c2(:),Jncfg1(:), JTp1(:)
      Integer :: njbl2 = 0
      Integer, allocatable :: JJ2c(:),JT2c1(:),JT2c2(:),Jncfg2(:), JTp2(:)

      End Module conf_jj


!======================================================================
      Subroutine alloc_cfg(m)
!======================================================================
!     allocate, deallocate or reallocate the arrays in module conf_jj
!----------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Integer :: m,i
      Integer, allocatable :: iarr(:)
      Real(8), allocatable :: rarr(:)

      if(m.le.0) then
       if(allocated(ip_state)) Deallocate(ip_state,ip_orb,IS_term,WC)
       mcfg=0; ncfg = 0; lcfg = 0; ne = 0
      elseif(.not.allocated(ip_state)) then
       mcfg = m; kcfg = mcfg*jcfg
       Allocate(ip_state(mcfg),IS_term(mcfg),ip_orb(kcfg),WC(mcfg))
       ncfg = 0; lcfg = 0
      elseif(m.le.mcfg) then
       Return
      elseif(ncfg.le.0) then
       Deallocate (ip_state,ip_orb,IS_term,WC); ncfg = 0; lcfg = 0
       Allocate(ip_state(mcfg),IS_term(mcfg),ip_orb(kcfg),WC(mcfg))
      else
       mcfg=m; i=lcfg/ncfg+1; if(i.gt.jcfg) jcfg=i; kcfg=mcfg*jcfg
       Allocate(iarr(lcfg))
       iarr(1:ncfg)=IS_term(1:ncfg); Deallocate(IS_term)
       Allocate(IS_term(m)); IS_term=0; IS_term(1:ncfg)=iarr(1:ncfg)
       iarr(1:ncfg)=ip_state(1:ncfg); Deallocate(ip_state)
       Allocate(ip_state(m)); ip_state=0; ip_state(1:ncfg)=iarr(1:ncfg)
       iarr(1:lcfg)=ip_orb(1:lcfg); Deallocate(ip_orb)
       Allocate(ip_orb(kcfg)); ip_orb=0; ip_orb(1:lcfg)=iarr(1:lcfg)
       Deallocate(iarr)

       Allocate(rarr(ncfg))
       rarr(1:ncfg)=WC(1:ncfg); Deallocate(WC)
       Allocate(WC(m)); WC=0.d0; WC(1:ncfg)=rarr(1:ncfg)
       Deallocate(rarr)

      end if

      End Subroutine alloc_cfg


!======================================================================
      Integer Function Iadd_cfg_jj(job)
!======================================================================
!     add new or detect existing CAS in conf_jj list;
!     returns position of given CAS in the list
!     job  =  'add'     -  just add
!          =  'detect'  -  return -ic if exist
!          =   others   -  add if not exist
!----------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Character(*), intent(in) :: job
      Integer :: i,ic,ip
      Integer, external :: Ifind_jjorb, Iadd_symc, Iadd_symt

      Iadd_cfg_jj = 0
      if(no.le.0) Return

      if(mcfg.eq.0) Call Alloc_cfg(icfg)
      Jtotal = Jintra(no)
      iconf = Iadd_symc(Jtotal,no,iq,kn)
      iterm = Iadd_symt(iconf,no,Jshell,Vshell,Jintra)

      Do i = 1,no; np(i)=Ifind_jjorb(nn(i),kn(i),in(i),2); End do

! ... check if we already have such state:

      if(job.ne.'add') then
      Do ic = 1,ncfg
       if(IS_term(ic).ne.iterm) Cycle
       ip = ip_state(ic); Iadd_cfg_jj = ic; if(job.eq.'detect') Iadd_cfg_jj=-ic
       Do i = 1,no; ip=ip+1
        if(np(i).ne.IP_orb(ip)) then; Iadd_cfg_jj=0; Exit; end if
       End do
       if(Iadd_cfg_jj.ne.0) Return
      End do
      end if

      ncfg=ncfg+1
      if(ncfg.eq.mcfg.or.lcfg+no.gt.kcfg) Call Alloc_cfg(mcfg+icfg)
      IS_term(ncfg)=iterm
      ip_state(ncfg)=lcfg
      Do i=1,no; lcfg=lcfg+1; IP_orb(lcfg)=np(i); End do
      Iadd_cfg_jj = ncfg

      End Function Iadd_cfg_jj


!======================================================================
      Subroutine Get_cfg_jj(ic)
!======================================================================
!     extract the CAS "ic" from conf_jj list
!----------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Integer :: ic, i,j,ip, ii

      iterm=IS_term(ic)
      Call Get_symt(iterm,iconf,no,Jshell,Vshell,Jintra)
      Call Get_symc(iconf,JTOTAL,no,nn,kn,ln,jn,iq,in)
      ip = ip_state(ic)
      Do i=1,no; ip=ip+1
       j = IP_orb(ip)
       Call Get_orb_jj(j,nn(i),kn(i),in(i),ii)
      End do

      End Subroutine Get_cfg_jj


!=======================================================================
      Integer Function Iort_conf_jj(kk)
!=======================================================================
!     define "strong" orthogonality between config.1 and config.2
!     based on the difference in kappa-values
!----------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Integer :: i,j,k,kk

      Iort_conf_jj = 0
      np1=iq1; np2=iq2
      Do i=1,no1
       Do j=1,no2
        if(np2(j).eq.0) Cycle
        if(kn1(i).ne.kn2(j)) Cycle
        k=min(np1(i),np2(j)); np1(i)=np1(i)-k; np2(j)=np2(j)-k
        if(np1(i).eq.0) Exit
       End do
      End do
      k = SUM(np2(1:no2));  if(k.gt.kk) Iort_conf_jj = 1

      End Function Iort_conf_jj


!=======================================================================
      Subroutine Def_Jblocks
!=======================================================================
!     define the number of J-blocks in the conf_jj list, njbl,
!     and their attributes:
!     JJc   -  total angular momentum
!     Jncgf -  number of ACF in the block
!     JTc1  -  index of the first ACS in the block
!     JTc2  -  index of the lastt ACS in the block
!     JTp   -  parity of the block
!
!     remark: don't be confused - JJc and Jncfg first used as auxiliary
!             arrays to save the memory (re-write ?)
!----------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Integer :: i,j

      Allocate(JJc(ncfg),Jncfg(ncfg))
      Call Get_cfg_jj(1)
      njbl=1; JJc(1) = Jtotal; Jncfg(1)=1

      Do i=2,ncfg
       Call Get_cfg_jj(i)
       if(Jtotal.eq.JJc(njbl)) then
        Jncfg(njbl)=i
       else
        njbl=njbl+1
        JJc(njbl)=Jtotal
        Jncfg(njbl)=i
       end if
      End do

      Allocate(JTc1(njbl),JTc2(njbl))
      JTc1(1:njbl)=JJc(1:njbl)
      JTc2(1:njbl)=Jncfg(1:njbl)
      Deallocate(JJc,Jncfg)
      Allocate(JJc(njbl),Jncfg(njbl))
      JJc=JTc1; JTc1(1)=1; Jncfg(1)=JTc2(1)
      Do i = 2,njbl
       JTc1(i) = JTc2(i-1)+1
       Jncfg(i) = JTc2(i) - JTc2(i-1)
      End do

      Allocate(JTp(njbl))
      Do i = 1,njbl
       Call Get_cfg_jj(JTc1(i))
       j=SUM(ln(1:no)*iq(1:no)); JTp(i)=(-1)**j
      End do

      End Subroutine Def_Jblocks


!======================================================================
      Integer Function no_ic_jj (ic)
!======================================================================
!     number of shells in state "ic"
!----------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Integer, intent(in) :: ic
      Integer, external :: Get_iconf, Get_no

      iterm=IS_term(ic)
      iconf=Get_iconf(iterm)
      no_ic_jj = Get_no(iconf)

      End Function no_ic_jj


!======================================================================
      Subroutine Save_cfg_jj(i)
!======================================================================
!     save(restore) curent state in position i = (1|2)
!----------------------------------------------------------------------
      Use conf_jj

      Implicit none
      Integer, Intent(in) :: i

      Select Case(i)
      Case(1)
       no1=no; nn1=nn; kn1=kn; ln1=ln; jn1=jn; iq1=iq; in1=in
       Jshell1=Jshell; Vshell1=Vshell; Jintra1=Jintra
       Jtotal1=Jtotal; iconf1=iconf; iterm1=iterm
      Case(2)
       no2=no; nn2=nn; kn2=kn; ln2=ln; jn2=jn; iq2=iq; in2=in
       Jshell2=Jshell; Vshell2=Vshell; Jintra2=Jintra
       Jtotal2=Jtotal; iconf2=iconf; iterm2=iterm
      Case(-1)
       no=no1; nn=nn1; kn=kn1; ln=ln1; jn=jn1; iq=iq1; in=in1
       Jshell=Jshell1; Vshell=Vshell1; Jintra=Jintra1
       Jtotal=Jtotal1; iconf=iconf1; iterm=iterm1
      Case(-2)
       no=no2; nn=nn2; kn=kn2; ln=ln2; jn=jn2; iq=iq2; in=in2
       Jshell=Jshell2; Vshell=Vshell2; Jintra=Jintra2
       Jtotal=Jtotal2; iconf=iconf2; iterm=iterm2
      Case Default
       Stop 'Save_cfg_jj: ???'
      End Select

      End Subroutine Save_cfg_jj


!======================================================================
      MODULE boef_list
!======================================================================
!     Containes two-electron integrals for matrix elements
!     in uncouple nlmj-representation.
!     The coefficient in the list are recorded by blocks -
!     all integrals for all operators under concideration for
!     the given matrix element    < 1, 2 | O | 3, 4>
!     The operator O is defined by external procedure: me_jj
!     This list is introduced to decrease the number of calls for
!     nj-symbol subroutine.
!----------------------------------------------------------------------
      Implicit none

      Integer :: nboef = 0       ! number of integrals
      Integer :: mboef = 0       ! current dimension of list
      Integer :: iboef = 50000   ! initial dimension
      Integer :: jboef = 0       ! first new element

! ... IB_int(1:mboef) - integral indentifier
! ... boef  (1:mboef) - correspondent angular coefficient

      Integer, allocatable :: IB_int(:)
      Real(8), allocatable :: boef(:)

      Integer :: nblk  = 0       ! number of blocks
      Integer :: mblk  = 0       ! current dimension of list
      Integer :: iblk  = 5000    ! initial dimentsion
      Integer :: kblk  = 0       ! block under consideration

! ... identifiers of block:

      Integer, allocatable :: id1(:),id2(:),id3(:),id4(:)

! ... current identifiers:

      Integer :: jd1,jd2,jd3,jd4

! ... ncblk - pointer on the last element in the block
! ... ipblk - ordering pointer

      Integer, allocatable :: ipblk(:),ncblk(:)

! ... incoding basis parameters, see subroutine check_boef:

      Integer, parameter :: ib10  = 2**10
      Integer, parameter :: ib20  = 2**20

      End MODULE BOEF_list


!======================================================================
      Subroutine alloc_boef(mm)
!======================================================================
!     allocate, deallocate or reallocate coefficient arrays IB_int and
!     Boef in module boef_list
!----------------------------------------------------------------------
      Use boef_list

      Implicit none
      Integer, intent(in)  :: mm
      Integer              :: m
      Integer, allocatable :: ia(:)
      Real(8), allocatable :: ra(:)

      m = mm; if(mm.lt.0) m=iboef
      if(m.le.0) then
       if(allocated(IB_int)) Deallocate(IB_int,Boef)
       mboef=0; nboef=0
      elseif(.not.allocated(IB_int)) then
       Allocate(IB_int(m), Boef(m));  mboef=m; nboef=0
      elseif(m.le.mboef) then
       Return
      elseif(nboef.eq.0) then
       Deallocate(IB_int,Boef)
       Allocate(IB_int(m), Boef(m));  mboef=m
      else
       Allocate(ra(nboef))
       ra(1:nboef)=boef(1:nboef); Deallocate(boef)
       Allocate(boef(m)); boef(1:nboef)=ra(1:nboef)
       Deallocate(ra)
       Allocate(ia(nboef))
       ia(1:nboef)=IB_INT(1:nboef); Deallocate(IB_INT)
       Allocate(IB_INT(m)); IB_INT(1:nboef)=ia(1:nboef)
       Deallocate(ia)
       mboef=m
      end if

      End Subroutine alloc_boef


!======================================================================
      Subroutine alloc_blk(mm)
!======================================================================
!     allocate, deallocate or reallocate coefficient "block" arrays
!     in module boef_list
!----------------------------------------------------------------------
      Use boef_list

      Implicit none
      Integer, intent(in)  :: mm
      Integer              :: m
      Integer, allocatable :: ia(:)

      m = mm; if(mm.lt.0) m=iblk
      if(m.le.0) then
       if(allocated(ipblk)) Deallocate(ipblk,ncblk,id1,id2,id3,id4)
       mblk = 0; nblk = 0
      elseif(.not.allocated(ipblk)) then
       Allocate(ipblk(m),ncblk(0:m),id1(m),id2(m),id3(m),id4(m))
       mblk = m; nblk = 0; ncblk = 0
      elseif(m.le.mblk) then
       Return
      elseif(nblk.eq.0) then
       Deallocate(ipblk,ncblk,id1,id2,id3,id4)
       Allocate(ipblk(m),ncblk(0:m),id1(m),id2(m),id3(m),id4(m))
       mblk = m
      else
       Allocate(ia(0:nblk))
       ia(1:nblk)=ipblk(1:nblk); Deallocate(ipblk)
       Allocate(ipblk(m)); ipblk(1:nblk)=ia(1:nblk)
       ia(0:nblk)=ncblk(0:nblk); Deallocate(ncblk)
       Allocate(ncblk(0:m)); ncblk(0:nblk)=ia(0:nblk)
       ia(1:nblk)=id1(1:nblk); Deallocate(id1)
       Allocate(id1(m)); id1(1:nblk)=ia(1:nblk)
       ia(1:nblk)=id2(1:nblk); Deallocate(id2)
       Allocate(id2(m)); id2(1:nblk)=ia(1:nblk)
       ia(1:nblk)=id3(1:nblk); Deallocate(id3)
       Allocate(id3(m)); id3(1:nblk)=ia(1:nblk)
       ia(1:nblk)=id4(1:nblk); Deallocate(id4)
       Allocate(id4(m)); id4(1:nblk)=ia(1:nblk)
       Deallocate(ia)
       mblk = m
      end if

      End Subroutine alloc_blk


!=======================================================================
      Subroutine Iadd_boef(C,int)
!=======================================================================
!     add new integral to the coefficient list in module boef_list
!-----------------------------------------------------------------------
      Use boef_list

      Implicit none
      Integer, Intent(in) :: int
      Real(8), Intent(in) :: C
      Integer :: i

      if(mboef.eq.0.or.nboef.eq.mboef) Call Alloc_boef(mboef+iboef)

      nboef=nboef+1; Boef(nboef)=C; IB_int(nboef)=int

      End Subroutine Iadd_boef


!=======================================================================
      Subroutine Check_boef(l1,j1,m1,l2,j2,m2,l3,j3,m3,l4,j4,m4)
!=======================================================================
!     Check if we already have the m.e. for given orbitals,
!     otherwise - calculate them.
!     Procedure uses "packing" the orbitals parameters, and that
!     restricts the max. l to 2**10=1024.
!----------------------------------------------------------------------
      USE boef_list

      Implicit none
      Integer :: l1,j1,m1, l2,j2,m2, l3,j3,m3, l4,j4,m4
      Integer :: i1,i2,i3,i4, k,l,m,ipm

      if(mblk.eq.0) Call Alloc_blk(iblk)

! ... prepare indentifiers:

      jd1 = l1*ib20+j1*ib10+j1+m1
      jd2 = l2*ib20+j2*ib10+j2+m2
      jd3 = l3*ib20+j3*ib10+j3+m3
      jd4 = l4*ib20+j4*ib10+j4+m4

! ... look for the same case in the list:

      k=1; l = nblk
    1 if(k.gt.l) go to 2

      m=(k+l)/2; ipm=ipblk(m)

      if(jd1.lt.id1(ipm)) then;      l = m - 1
      elseif(jd1.gt.id1(ipm)) then;  k = m + 1
      else

      if(jd2.lt.id2(ipm)) then;      l = m - 1
      elseif(jd2.gt.id2(ipm)) then;  k = m + 1
      else

      if(jd3.lt.id3(ipm)) then;      l = m - 1
      elseif(jd3.gt.id3(ipm)) then;  k = m + 1
      else

      if(jd4.lt.id4(ipm)) then;      l = m - 1
      elseif(jd4.gt.id4(ipm)) then;  k = m + 1

      else
       kblk = ipm
       Return
      end if; end if; end if; end if

      go to 1
    2 Continue

! ... new block:

      jboef = nboef + 1
      Call me_jj(l1,j1,m1,l2,j2,m2,l3,j3,m3,l4,j4,m4,+1)
      Call me_jj(l1,j1,m1,l2,j2,m2,l4,j4,m4,l3,j3,m3,-1)

      nblk = nblk + 1;  ncblk(nblk) = nboef; kblk = nblk
      id1(nblk)=jd1; id2(nblk)=jd2; id3(nblk)=jd3; id4(nblk)=jd4

      if(k.eq.nblk) then
       ipblk(k)=nblk
      else
       Do m = nblk,k+1,-1; ipblk(m) = ipblk(m-1); End do
       ipblk(k)=nblk
      end if

! ... it is time for re-allocation:

      if(nblk.eq.mblk) Call Alloc_blk(mblk+iblk)

      End Subroutine Check_boef


!======================================================================
      Subroutine me_jj(l1,j1,m1,l2,j2,m2,l3,j3,m3,l4,j4,m4,ide)
!======================================================================
!     computes angular coefficients for the two-electron Coulomb
!     interaction in uncouple njm-representation;
!
!     should be called twice with interchange 3 <-> 4 and
!     ide = +1  ->  direct interaction
!     ide = -1  ->  exchange interaction
!
!     coefficients are stored in module "boef_list"
!     through call to Iadd_boef
!----------------------------------------------------------------------
      Implicit none
      Integer :: l1,j1,m1,l2,j2,m2,l3,j3,m3,l4,j4,m4,ide
      Integer :: KL,KM, KL1,KL2, KM1,KM2, k,m, kz,int, kk
      Real(8) :: C
      Real(8), external :: Z_3j2, Cjkj

      m = m1-m3; !if(m.ne.m4-m2) Return

! ... define the range of multipole indeces:

      KL1=IABS(j1-j3);  KL2=IABS(j2-j4);  KL=MAX0(KL1,KL2)/2
      KM1=     j1+j3 ;  KM2=     j2+j4 ;  KM=MIN0(KM1,KM2)/2
      if(KM.lt.KL) Return

      Do k = kl,km; kk=k+k

       C = Z_3j2(j1,-m1,kk,m1-m3,j3,m3) * Z_3j2(j2,-m2,kk,m2-m4,j4,m4)
       if(C.eq.0.d0) Cycle
       kz = (kk-m + j1-m1 + j2-m2)/2
       if(mod(kz,2).ne.0) C = -C
       C = C * ide

       C = C * Cjkj(j1,k,j3) * Cjkj(j2,k,j4)
       if(mod(k,2).eq.1) C=-C

       int = (k+1)*ide

       Call Iadd_boef(C,int)

      End do

      End Subroutine me_jj



!======================================================================
      Module det_jq
!======================================================================
!     contains data for expanding the j^q subshell wave function
!     over deterninant w.f.
!
!     Idet_jq  - identifier of subshell detetrminant w.f. as list of
!                the one-electron lsm-orbitals according to the common
!                list (-1/2, 1/2, -3/2, 3/2, ...)
!     MD_jq    - list of the correspondent total magnetic quantum numbers
!
!     Expansion coefficient for given determinant 'id' and term 'it'
!     is presented as  sqrt [ ID_jq(id,it) / JD_jq(id,it) ].
!
!     Routines:
!
!     DET_sh_jq (j,q,id,MJ,Idet) - gives the determinant 'id' and total MJ
!                                  for the subshell j^q
!
!     DETC_jq (j,q,it,id) - gives the expansion coefficient
!                           for the subshell 'jq', term 'it' and determinant 'id'
!----------------------------------------------------------------------
      Implicit none

! ... subshell [1/2]^2

      Integer, parameter :: kd_j1_q2 =  1
      Integer, parameter :: nt_j1_q2 =  1
      Integer :: Idet_j1_q2(2,kd_j1_q2)
      Integer :: MD_j1_q2(kd_j1_q2)
      Integer :: ID_j1_q2(nt_j1_q2,kd_j1_q2)
      Integer :: JD_j1_q2(nt_j1_q2,kd_j1_q2)

      Data Idet_j1_q2/ &
          1,    2  /
      Data MD_j1_q2/ &
          0  /
      Data ID_j1_q2/ &
         -1  /
      Data JD_j1_q2/ &
          1  /

! ... subshell [3/2]^2

      Integer, parameter :: kd_j3_q2 =  6
      Integer, parameter :: nt_j3_q2 =  2
      Integer :: Idet_j3_q2(2,kd_j3_q2)
      Integer :: MD_j3_q2(kd_j3_q2)
      Integer :: ID_j3_q2(nt_j3_q2,kd_j3_q2)
      Integer :: JD_j3_q2(nt_j3_q2,kd_j3_q2)

      Data Idet_j3_q2/ &
          1,    2, &
          1,    3, &
          1,    4, &
          2,    3, &
          2,    4, &
          3,    4  /
      Data MD_j3_q2/ &
          0,   -4,    2,   -2,    4,    0  /
      Data ID_j3_q2/ &
          1,   -1, &
          0,    1, &
          0,   -1, &
          0,    1, &
          0,   -1, &
         -1,   -1  /
      Data JD_j3_q2/ &
          2,    2, &
          1,    1, &
          1,    1, &
          1,    1, &
          1,    1, &
          2,    2  /

! ... subshell [3/2]^3

      Integer, parameter :: kd_j3_q3 =  4
      Integer, parameter :: nt_j3_q3 =  1
      Integer :: Idet_j3_q3(3,kd_j3_q3)
      Integer :: MD_j3_q3(kd_j3_q3)
      Integer :: ID_j3_q3(nt_j3_q3,kd_j3_q3)
      Integer :: JD_j3_q3(nt_j3_q3,kd_j3_q3)

      Data Idet_j3_q3/ &
          1,    2,    3, &
          1,    2,    4, &
          1,    3,    4, &
          2,    3,    4  /
      Data MD_j3_q3/ &
         -3,    3,   -1,    1  /
      Data ID_j3_q3/ &
          1, &
          1, &
         -1, &
         -1  /
      Data JD_j3_q3/ &
          1, &
          1, &
          1, &
          1  /
! ... subshell [3/2]^4

      Integer, parameter :: kd_j3_q4 =  1
      Integer, parameter :: nt_j3_q4 =  1
      Integer :: Idet_j3_q4(4,kd_j3_q4)
      Integer :: MD_j3_q4(kd_j3_q4)
      Integer :: ID_j3_q4(nt_j3_q4,kd_j3_q4)
      Integer :: JD_j3_q4(nt_j3_q4,kd_j3_q4)

      Data Idet_j3_q4/ &
          1,    2,    3,    4  /

      Data MD_j3_q4/ &
          0  /

      Data ID_j3_q4/ &
         -1  /

      Data JD_j3_q4/ &
          1  /

! ... subshell [5/2]^2

      Integer, parameter :: kd_j5_q2 = 15
      Integer, parameter :: nt_j5_q2 =  3

      Integer :: Idet_j5_q2(2,kd_j5_q2)
      Integer :: MD_j5_q2(kd_j5_q2)
      Integer :: ID_j5_q2(nt_j5_q2,kd_j5_q2)
      Integer :: JD_j5_q2(nt_j5_q2,kd_j5_q2)

      Data Idet_j5_q2/ &
          1,    2, &
          1,    3, &
          1,    4, &
          1,    5, &
          1,    6, &
          2,    3, &
          2,    4, &
          2,    5, &
          2,    6, &
          3,    4, &
          3,    5, &
          3,    6, &
          4,    5, &
          4,    6, &
          5,    6  /

      Data MD_j5_q2/ &
          0,   -4,    2,   -6,    4,   -2,    4,   -4,    6,    0, &
         -8,    2,   -2,    8,    0  /

      Data ID_j5_q2/ &
         -1,    8,   -2, &
          0,   -9,    5, &
          0,    2,   -5, &
          0,    0,    1, &
          0,   -5,   -9, &
          0,   -2,    5, &
          0,    9,   -5, &
          0,    5,    9, &
          0,    0,   -1, &
          1,   -1,   -9, &
          0,    0,    1, &
          0,   -5,   -2, &
          0,    5,    2, &
          0,    0,   -1, &
         -1,  -25,   -1  /

      Data JD_j5_q2/ &
          3,   21,    7, &
          1,   14,   14, &
          1,    7,    7, &
          1,    1,    1, &
          1,   14,   14, &
          1,    7,    7, &
          1,   14,   14, &
          1,   14,   14, &
          1,    1,    1, &
          3,   42,   14, &
          1,    1,    1, &
          1,    7,    7, &
          1,    7,    7, &
          1,    1,    1, &
          3,   42,   14  /

! ... subshell [5/2]^3

      Integer, parameter :: kd_j5_q3 = 20
      Integer, parameter :: nt_j5_q3 =  3

      Integer :: Idet_j5_q3(3,kd_j5_q3)
      Integer :: MD_j5_q3(kd_j5_q3)
      Integer :: ID_j5_q3(nt_j5_q3,kd_j5_q3)
      Integer :: JD_j5_q3(nt_j5_q3,kd_j5_q3)

      Data Idet_j5_q3/ &
          1,    2,    3, &
          1,    2,    4, &
          1,    2,    5, &
          1,    2,    6, &
          1,    3,    4, &
          1,    3,    5, &
          1,    3,    6, &
          1,    4,    5, &
          1,    4,    6, &
          1,    5,    6, &
          2,    3,    4, &
          2,    3,    5, &
          2,    3,    6, &
          2,    4,    5, &
          2,    4,    6, &
          2,    5,    6, &
          3,    4,    5, &
          3,    4,    6, &
          3,    5,    6, &
          4,    5,    6  /

      Data MD_j5_q3/ &
         -3,    3,   -5,    5,   -1,   -9,    1,   -3,    7,   -1, &
          1,   -7,    3,   -1,    9,    1,   -5,    5,   -3,    3  /

      Data ID_j5_q3/ &
          1,    8,    5, &
          1,   -8,    5, &
          1,    0,    1, &
          1,    0,    1, &
         -1,   -1,   -5, &
          0,    0,   -1, &
          0,   -5,   -2, &
          0,   -5,   16, &
          0,    0,    1, &
          1,   -1,   -5, &
         -1,    1,   -5, &
          0,    0,   -1, &
          0,   -5,  -16, &
          0,   -5,    2, &
          0,    0,    1, &
          1,    1,   -5, &
         -1,    0,    1, &
         -1,    0,    1, &
          1,   -8,   -5, &
          1,    8,   -5  /

      Data JD_j5_q3/ &
          2,   21,   42, &
          2,   21,   42, &
          2,    1,    2, &
          2,    1,    2, &
          2,    7,   14, &
          1,    1,    1, &
          1,    7,    7, &
          1,   21,   21, &
          1,    1,    1, &
          2,    7,   14, &
          2,    7,   14, &
          1,    1,    1, &
          1,   21,   21, &
          1,    7,    7, &
          1,    1,    1, &
          2,    7,   14, &
          2,    1,    2, &
          2,    1,    2, &
          2,   21,   42, &
          2,   21,   42  /

! ... subshell [5/2]^4

      Integer, parameter :: kd_j5_q4 = 15
      Integer, parameter :: nt_j5_q4 =  3

      Integer :: Idet_j5_q4(4,kd_j5_q4)
      Integer :: MD_j5_q4(kd_j5_q4)
      Integer :: ID_j5_q4(nt_j5_q4,kd_j5_q4)
      Integer :: JD_j5_q4(nt_j5_q4,kd_j5_q4)

      Data Idet_j5_q4/ &
          1,    2,    3,    4, &
          1,    2,    3,    5, &
          1,    2,    3,    6, &
          1,    2,    4,    5, &
          1,    2,    4,    6, &
          1,    2,    5,    6, &
          1,    3,    4,    5, &
          1,    3,    4,    6, &
          1,    3,    5,    6, &
          1,    4,    5,    6, &
          2,    3,    4,    5, &
          2,    3,    4,    6, &
          2,    3,    5,    6, &
          2,    4,    5,    6, &
          3,    4,    5,    6  /

      Data MD_j5_q4/ &
          0,   -8,    2,   -2,    8,    0,   -6,    4,   -4,    2, &
         -4,    6,   -2,    4,    0  /

      Data ID_j5_q4/ &
          1,  -25,   -1, &
          0,    0,    1, &
          0,   -5,   -2, &
          0,    5,    2, &
          0,    0,   -1, &
         -1,   -1,   -9, &
          0,    0,   -1, &
          0,    5,    9, &
          0,   -9,    5, &
          0,    2,   -5, &
          0,   -5,   -9, &
          0,    0,    1, &
          0,   -2,    5, &
          0,    9,   -5, &
          1,    8,   -2  /

      Data JD_j5_q4/ &
          3,   42,   14, &
          1,    1,    1, &
          1,    7,    7, &
          1,    7,    7, &
          1,    1,    1, &
          3,   42,   14, &
          1,    1,    1, &
          1,   14,   14, &
          1,   14,   14, &
          1,    7,    7, &
          1,   14,   14, &
          1,    1,    1, &
          1,    7,    7, &
          1,   14,   14, &
          3,   21,    7  /

! ... subshell [5/2]^5

      Integer, parameter :: kd_j5_q5 =  6
      Integer, parameter :: nt_j5_q5 =  1

      Integer :: Idet_j5_q5(5,kd_j5_q5)
      Integer :: MD_j5_q5(kd_j5_q5)
      Integer :: ID_j5_q5(nt_j5_q5,kd_j5_q5)
      Integer :: JD_j5_q5(nt_j5_q5,kd_j5_q5)

      Data Idet_j5_q5/ &
          1,    2,    3,    4,    5, &
          1,    2,    3,    4,    6, &
          1,    2,    3,    5,    6, &
          1,    2,    4,    5,    6, &
          1,    3,    4,    5,    6, &
          2,    3,    4,    5,    6  /

      Data MD_j5_q5/ &
         -5,    5,   -3,    3,   -1,    1  /

      Data ID_j5_q5/ &
         -1, &
         -1, &
          1, &
          1, &
         -1, &
         -1  /

      Data JD_j5_q5/ &
          1, &
          1, &
          1, &
          1, &
          1, &
          1  /

! ... subshell [5/2]^6

      Integer, parameter :: kd_j5_q6 =  1
      Integer, parameter :: nt_j5_q6 =  1

      Integer :: Idet_j5_q6(6,kd_j5_q6)
      Integer :: MD_j5_q6(kd_j5_q6)
      Integer :: ID_j5_q6(nt_j5_q6,kd_j5_q6)
      Integer :: JD_j5_q6(nt_j5_q6,kd_j5_q6)

      Data Idet_j5_q6/ &
          1,    2,    3,    4,    5,    6  /

      Data MD_j5_q6/ &
          0  /

      Data ID_j5_q6/ &
          1  /

      Data JD_j5_q6/ &
          1  /

! ... subshell [7/2]^2

      Integer, parameter :: kd_j7_q2 = 28
      Integer, parameter :: nt_j7_q2 =  4

      Integer :: Idet_j7_q2(2,kd_j7_q2)
      Integer :: MD_j7_q2(kd_j7_q2)
      Integer :: ID_j7_q2(nt_j7_q2,kd_j7_q2)
      Integer :: JD_j7_q2(nt_j7_q2,kd_j7_q2)

      Data Idet_j7_q2/ &
          1,    2, &
          1,    3, &
          1,    4, &
          1,    5, &
          1,    6, &
          1,    7, &
          1,    8, &
          2,    3, &
          2,    4, &
          2,    5, &
          2,    6, &
          2,    7, &
          2,    8, &
          3,    4, &
          3,    5, &
          3,    6, &
          3,    7, &
          3,    8, &
          4,    5, &
          4,    6, &
          4,    7, &
          4,    8, &
          5,    6, &
          5,    7, &
          5,    8, &
          6,    7, &
          6,    8, &
          7,    8  /

      Data MD_j7_q2/ &
          0,   -4,    2,   -6,    4,   -8,    6,   -2,    4,   -4, &
          6,   -6,    8,    0,   -8,    2,  -10,    4,   -2,    8, &
         -4,   10,    0,  -12,    2,   -2,   12,    0  /

      Data ID_j7_q2/ &
          1,  -25,   81,  -25, &
          0,   10,  -24,    7, &
          0,   -5,   27,  -35, &
          0,    0,   -4,    7, &
          0,    5,    1,   -7, &
          0,    0,    7,   15, &
          0,    0,   -7,   -4, &
          0,    5,  -27,   35, &
          0,  -10,   24,   -7, &
          0,   -5,   -1,    7, &
          0,    0,    4,   -7, &
          0,    0,    7,    4, &
          0,    0,   -7,  -15, &
         -1,    3,    9,  -27, &
          0,    0,  -15,    7, &
          0,    8,  -15,  -14, &
          0,    0,    0,    1, &
          0,   -1,  -15,   -5, &
          0,   -8,   15,   14, &
          0,    0,   15,   -7, &
          0,    1,   15,    5, &
          0,    0,    0,   -1, &
          1,    1, -169,  -25, &
          0,    0,    0,    1, &
          0,   -1,   -5,   -1, &
          0,    1,    5,    1, &
          0,    0,    0,   -1, &
         -1,   -7,   -7,   -1  /

      Data JD_j7_q2/ &
          4,   84,  308,  132, &
          1,   21,   77,   33, &
          1,   42,   77,   66, &
          1,    1,   11,   11, &
          1,   14,  154,   11, &
          1,    1,   22,   22, &
          1,    1,   11,   11, &
          1,   42,   77,   66, &
          1,   21,   77,   33, &
          1,   14,  154,   11, &
          1,    1,   11,   11, &
          1,    1,   11,   11, &
          1,    1,   22,   22, &
          4,   28,  308,   44, &
          1,    1,   22,   22, &
          1,   21,   77,   33, &
          1,    1,    1,    1, &
          1,    6,   22,   33, &
          1,   21,   77,   33, &
          1,    1,   22,   22, &
          1,    6,   22,   33, &
          1,    1,    1,    1, &
          4,   84,  308,  132, &
          1,    1,    1,    1, &
          1,    2,   11,   22, &
          1,    2,   11,   22, &
          1,    1,    1,    1, &
          4,   12,   44,  132  /

! ... subshell [7/2]^3

      Integer, parameter :: kd_j7_q3 = 56
      Integer, parameter :: nt_j7_q3 =  6

      Integer :: Idet_j7_q3(3,kd_j7_q3)
      Integer :: MD_j7_q3(kd_j7_q3)
      Integer :: ID_j7_q3(nt_j7_q3,kd_j7_q3)
      Integer :: JD_j7_q3(nt_j7_q3,kd_j7_q3)

      Data Idet_j7_q3/ &
          1,    2,    3, &
          1,    2,    4, &
          1,    2,    5, &
          1,    2,    6, &
          1,    2,    7, &
          1,    2,    8, &
          1,    3,    4, &
          1,    3,    5, &
          1,    3,    6, &
          1,    3,    7, &
          1,    3,    8, &
          1,    4,    5, &
          1,    4,    6, &
          1,    4,    7, &
          1,    4,    8, &
          1,    5,    6, &
          1,    5,    7, &
          1,    5,    8, &
          1,    6,    7, &
          1,    6,    8, &
          1,    7,    8, &
          2,    3,    4, &
          2,    3,    5, &
          2,    3,    6, &
          2,    3,    7, &
          2,    3,    8, &
          2,    4,    5, &
          2,    4,    6, &
          2,    4,    7, &
          2,    4,    8, &
          2,    5,    6, &
          2,    5,    7, &
          2,    5,    8, &
          2,    6,    7, &
          2,    6,    8, &
          2,    7,    8, &
          3,    4,    5, &
          3,    4,    6, &
          3,    4,    7, &
          3,    4,    8, &
          3,    5,    6, &
          3,    5,    7, &
          3,    5,    8, &
          3,    6,    7, &
          3,    6,    8, &
          3,    7,    8, &
          4,    5,    6, &
          4,    5,    7, &
          4,    5,    8, &
          4,    6,    7, &
          4,    6,    8, &
          4,    7,    8, &
          5,    6,    7, &
          5,    6,    8, &
          5,    7,    8, &
          6,    7,    8  /

      Data MD_j7_q3/ &
         -3,    3,   -5,    5,   -7,    7,   -1,   -9,    1,  -11, &
          3,   -3,    7,   -5,    9,   -1,  -13,    1,   -3,   11, &
         -1,    1,   -7,    3,   -9,    5,   -1,    9,   -3,   11, &
          1,  -11,    3,   -1,   13,    1,   -5,    5,   -7,    7, &
         -3,  -15,   -1,   -5,    9,   -3,    3,   -9,    5,    1, &
         15,    3,   -7,    7,   -5,    5  /

      Data ID_j7_q3/ &
          1,   -3,    3,   25,    3,    5, &
          1,   -3,   -3,  -25,    3,    5, &
          1,    0,    5,   25,    4,   25, &
          1,    0,   -5,  -25,    4,   25, &
          1,    0,    0, -175,    3,   20, &
          1,    0,    0,  175,    3,   20, &
         -1,    3,   -8,  -40,   -7,  -21, &
          0,    0,    0,  -80,   -4,   -1, &
          0,   -1,  -98,  -90,   -7,  -28, &
          0,    0,    0,    0,   -4,   -3, &
          0,    3,  -28,  -60,    5, -100, &
          0,    6,  -49,   15,   10,   64, &
          0,    0,    0,  -45,    5,    3, &
          0,    0,   -7, -315,   -5,  375, &
          0,    0,    0,   35,    4,   64, &
          1,   -8,   -9,  -55,    0, -175, &
          0,    0,    0,    0,    0,   -1, &
          0,    2,    7, -135,    2,  -64, &
          0,    1,   63,  -45,  -45,  243, &
          0,    0,    0,    0,   -3,    4, &
         -1,   -7,    1,  -35,    7,  -28, &
         -1,    3,    8,   40,   -7,  -21, &
          0,    0,    0,  -45,   -5,   -3, &
          0,   -6,  -49,   15,  -10,  -64, &
          0,    0,    0,   35,   -4,  -64, &
          0,    0,   -7, -315,    5, -375, &
          0,    1,  -98,  -90,    7,   28, &
          0,    0,    0,  -80,    4,    1, &
          0,   -3,  -28,  -60,   -5,  100, &
          0,    0,    0,    0,    4,    3, &
          1,   -8,    9,   55,    0, -175, &
          0,    0,    0,    0,    3,   -4, &
          0,   -1,   63,  -45,   45, -243, &
          0,   -2,    7, -135,   -2,   64, &
          0,    0,    0,    0,    0,    1, &
         -1,   -7,   -1,   35,    7,  -28, &
         -1,    0,   -5,  361,    1,   36, &
         -1,    0,    5, -361,    1,   36, &
         -1,    0,    0,   -7,  -12,   45, &
         -1,    0,    0,    7,  -12,   45, &
          1,    3,   -1, -477,    1,  -20, &
          0,    0,    0,    0,    0,   -1, &
          0,    1,    7,  -81,    5,   -5, &
          0,    0,   35,   63, -100,  192, &
          0,    0,    0,  -28,  -45,   20, &
         -1,    0,   25, -343,    8,   -5, &
          1,    3,    1,  477,    1,  -20, &
          0,    0,    0,  -28,   45,  -20, &
          0,    0,   35,   63,  100, -192, &
          0,   -1,    7,  -81,   -5,    5, &
          0,    0,    0,    0,    0,    1, &
         -1,    0,  -25,  343,    8,   -5, &
          1,    0,    0,  112,  -27,    5, &
          1,    0,    0, -112,  -27,    5, &
         -1,    0,    5,  -49,    1,   -1, &
         -1,    0,   -5,   49,    1,   -1  /

      Data JD_j7_q3/ &
          3,   14,   11, 3003,   22,  143, &
          3,   14,   11, 3003,   22,  143, &
          3,    1,   44, 1716,   11,  143, &
          3,    1,   44, 1716,   11,  143, &
          3,    1,    1,  429,   77,   91, &
          3,    1,    1,  429,   77,   91, &
          3,   70,   55, 3003,   22,  143, &
          1,    1,    1,  143,   11,   13, &
          1,   70,  165, 1001,   66,  143, &
          1,    1,    1,    1,    7,    7, &
          1,   10,  165,  143,  462, 1001, &
          1,   35,  660, 4004,   33,  143, &
          1,    1,    1,  143,   11,   13, &
          1,    1,  132,  572,  231, 1001, &
          1,    1,    1,  143,   77,   91, &
          3,  105,  110,  546,    1,  429, &
          1,    1,    1,    1,    1,    1, &
          1,   15,  110,  286,   11,  429, &
          1,   10,  220,  572,  154, 1001, &
          1,    1,    1,    1,    7,    7, &
          3,   30,  110,  858,   22,  429, &
          3,   70,   55, 3003,   22,  143, &
          1,    1,    1,  143,   11,   13, &
          1,   35,  660, 4004,   33,  143, &
          1,    1,    1,  143,   77,   91, &
          1,    1,  132,  572,  231, 1001, &
          1,   70,  165, 1001,   66,  143, &
          1,    1,    1,  143,   11,   13, &
          1,   10,  165,  143,  462, 1001, &
          1,    1,    1,    1,    7,    7, &
          3,  105,  110,  546,    1,  429, &
          1,    1,    1,    1,    7,    7, &
          1,   10,  220,  572,  154, 1001, &
          1,   15,  110,  286,   11,  429, &
          1,    1,    1,    1,    1,    1, &
          3,   30,  110,  858,   22,  429, &
          3,    1,   44, 1716,   11,  143, &
          3,    1,   44, 1716,   11,  143, &
          3,    1,    1,  429,   77,   91, &
          3,    1,    1,  429,   77,   91, &
          3,   14,  132, 1646,   66,  143, &
          1,    1,    1,    1,    1,    1, &
          1,    2,   66,  286,   66,  143, &
          1,    1,  132,  572,  231, 1001, &
          1,    1,    1,  143,   77,   91, &
          3,    1,  132, 1716,   33,  143, &
          3,   14,  132, 1646,   66,  143, &
          1,    1,    1,  143,   77,   91, &
          1,    1,  132,  572,  231, 1001, &
          1,    2,   66,  286,   66,  143, &
          1,    1,    1,    1,    1,    1, &
          3,    1,  132, 1716,   33,  143, &
          3,    1,    1,  429,   77,   91, &
          3,    1,    1,  429,   77,   91, &
          3,    1,   11,  429,   11,  143, &
          3,    1,   11,  429,   11,  143  /

! ... subshell [7/2]^4

      Integer, parameter :: kd_j7_q4 = 70
      Integer, parameter :: nt_j7_q4 =  8

      Integer :: Idet_j7_q4(4,kd_j7_q4)
      Integer :: MD_j7_q4(kd_j7_q4)
      Integer :: ID_j7_q4(nt_j7_q4,kd_j7_q4)
      Integer :: JD_j7_q4(nt_j7_q4,kd_j7_q4)

      Data Idet_j7_q4/ &
          1,    2,    3,    4, &
          1,    2,    3,    5, &
          1,    2,    3,    6, &
          1,    2,    3,    7, &
          1,    2,    3,    8, &
          1,    2,    4,    5, &
          1,    2,    4,    6, &
          1,    2,    4,    7, &
          1,    2,    4,    8, &
          1,    2,    5,    6, &
          1,    2,    5,    7, &
          1,    2,    5,    8, &
          1,    2,    6,    7, &
          1,    2,    6,    8, &
          1,    2,    7,    8, &
          1,    3,    4,    5, &
          1,    3,    4,    6, &
          1,    3,    4,    7, &
          1,    3,    4,    8, &
          1,    3,    5,    6, &
          1,    3,    5,    7, &
          1,    3,    5,    8, &
          1,    3,    6,    7, &
          1,    3,    6,    8, &
          1,    3,    7,    8, &
          1,    4,    5,    6, &
          1,    4,    5,    7, &
          1,    4,    5,    8, &
          1,    4,    6,    7, &
          1,    4,    6,    8, &
          1,    4,    7,    8, &
          1,    5,    6,    7, &
          1,    5,    6,    8, &
          1,    5,    7,    8, &
          1,    6,    7,    8, &
          2,    3,    4,    5, &
          2,    3,    4,    6, &
          2,    3,    4,    7, &
          2,    3,    4,    8, &
          2,    3,    5,    6, &
          2,    3,    5,    7, &
          2,    3,    5,    8, &
          2,    3,    6,    7, &
          2,    3,    6,    8, &
          2,    3,    7,    8, &
          2,    4,    5,    6, &
          2,    4,    5,    7, &
          2,    4,    5,    8, &
          2,    4,    6,    7, &
          2,    4,    6,    8, &
          2,    4,    7,    8, &
          2,    5,    6,    7, &
          2,    5,    6,    8, &
          2,    5,    7,    8, &
          2,    6,    7,    8, &
          3,    4,    5,    6, &
          3,    4,    5,    7, &
          3,    4,    5,    8, &
          3,    4,    6,    7, &
          3,    4,    6,    8, &
          3,    4,    7,    8, &
          3,    5,    6,    7, &
          3,    5,    6,    8, &
          3,    5,    7,    8, &
          3,    6,    7,    8, &
          4,    5,    6,    7, &
          4,    5,    6,    8, &
          4,    5,    7,    8, &
          4,    6,    7,    8, &
          5,    6,    7,    8  /

      Data MD_j7_q4/ &
          0,   -8,    2,  -10,    4,   -2,    8,   -4,   10,    0, &
        -12,    2,   -2,   12,    0,   -6,    4,   -8,    6,   -4, &
        -16,   -2,   -6,    8,   -4,    2,  -10,    4,    0,   14, &
          2,   -8,    6,   -6,    4,   -4,    6,   -6,    8,   -2, &
        -14,    0,   -4,   10,   -2,    4,   -8,    6,    2,   16, &
          4,   -6,    8,   -4,    6,    0,  -12,    2,   -2,   12, &
          0,  -10,    4,   -8,    2,   -4,   10,   -2,    8,    0  /

      Data ID_j7_q4/ &
         -1,    8,   -9,   -2,   -6,   45,    0,   -7, &
          0,    0,  -15,    7,    0,   -3,   -1,    1, &
          0,    4,  -15,   -7,   -3,   75,    0,   -7, &
          0,    0,    0,    1,    0,    0,   -5,    1, &
          0,   -1,  -15,   -5,   -3,    3,    3, -125, &
          0,   -4,   15,    7,    3,  -75,    0,    7, &
          0,    0,   15,   -7,    0,    3,   -1,   -1, &
          0,    1,   15,    5,    3,   -3,    3,  125, &
          0,    0,    0,   -1,    0,    0,   -5,   -1, &
          1,   -2,   -2,  -25,    3,   10,    0, -175, &
          0,    0,    0,    1,    0,    0,    0,    1, &
          0,   -1,   -5,   -1,   -1,   -1,    1,  -27, &
          0,    1,    5,    1,    1,    1,    1,   27, &
          0,    0,    0,   -1,    0,    0,    0,   -1, &
         -1,   -1,  -32,    2,   27,   -5,    0,  -56, &
          0,    0,    2,   -7,    0,    8,    3,   -1, &
          0,   -5,   -1,    7,   -5, -529,    1,   27, &
          0,    0,   -7,  -15,    0,   -7,    3, -135, &
          0,    0,    7,    2,    0,   -7,   -6,   25, &
          0,    5,  -12,    7,   20,   27,   -4,   25, &
          0,    0,    0,    0,    0,    0,    0,   -1, &
          0,    0,    0,    0,   20,    6,   -4,    5, &
          0,    0,    0,    0,    0,  -42,    2,  -27, &
          0,    0,    0,    0,    0,  -21,  -16,   20, &
          0,   -5,   12,   -7,   20,   27,   -4,   25, &
          0,   -5,   27,  -35,  125,   -3,   -7,  -35, &
          0,    0,    0,    0,    0,    0,    2,    5, &
          0,    0,    0,    0,    5,  -12,   16, -512, &
          0,    0,    0,    0,    5,    3,    1,   45, &
          0,    0,    0,    0,    0,    0,    0,   -1, &
          0,    5,  -27,   35,  125,   -3,   -7,  -35, &
          0,    0,    7,   15,    0,   -7,    3, -135, &
          0,    0,   -7,   -2,    0,   -7,   -6,   25, &
          0,    0,    2,   -7,    0,   -8,   -3,    1, &
          0,   -5,   -1,    7,    5,  529,   -1,  -27, &
          0,    5,    1,   -7,    5,  529,    1,  -27, &
          0,    0,   -2,    7,    0,   -8,    3,    1, &
          0,    0,   -7,   -2,    0,    7,   -6,  -25, &
          0,    0,    7,   15,    0,    7,    3,  135, &
          0,    5,  -27,   35, -125,    3,   -7,   35, &
          0,    0,    0,    0,    0,    0,    0,   -1, &
          0,    0,    0,    0,    5,    3,   -1,   45, &
          0,    0,    0,    0,    5,  -12,  -16, -512, &
          0,    0,    0,    0,    0,    0,   -2,    5, &
          0,   -5,   27,  -35, -125,    3,   -7,   35, &
          0,   -5,   12,   -7,  -20,  -27,   -4,  -25, &
          0,    0,    0,    0,    0,  -21,   16,   20, &
          0,    0,    0,    0,    0,  -42,   -2,  -27, &
          0,    0,    0,    0,   20,    6,    4,    5, &
          0,    0,    0,    0,    0,    0,    0,   -1, &
          0,    5,  -12,    7,  -20,  -27,   -4,  -25, &
          0,    0,    7,    2,    0,    7,   -6,  -25, &
          0,    0,   -7,  -15,    0,    7,    3,  135, &
          0,    5,    1,   -7,   -5, -529,   -1,   27, &
          0,    0,   -2,    7,    0,    8,   -3,   -1, &
         -1,    1,   32,   -2,   27,   -5,    0,  -56, &
          0,    0,    0,   -1,    0,    0,    0,    1, &
          0,    1,    5,    1,   -1,   -1,    1,  -27, &
          0,   -1,   -5,   -1,    1,    1,    1,   27, &
          0,    0,    0,    1,    0,    0,    0,   -1, &
          1,    2,    2,   25,    3,   10,    0, -175, &
          0,    0,    0,    1,    0,    0,    5,   -1, &
          0,   -1,  -15,   -5,    3,   -3,   -3,  125, &
          0,    0,   15,   -7,    0,   -3,   -1,    1, &
          0,   -4,   15,    7,   -3,   75,    0,   -7, &
          0,    1,   15,    5,   -3,    3,   -3, -125, &
          0,    0,    0,   -1,    0,    0,    5,    1, &
          0,    4,  -15,   -7,    3,  -75,    0,    7, &
          0,    0,  -15,    7,    0,    3,   -1,   -1, &
         -1,   -8,    9,    2,   -6,   45,    0,   -7  /

      Data JD_j7_q4/ &
          6,   21,  154,   33,   77,  182,    1,  858, &
          1,    1,   44,   44,    1,   13,    4,   52, &
          1,   21,  154,   33,   77,  182,    1,  143, &
          1,    1,    1,    2,    1,    1,   14,    7, &
          1,   12,   44,   66,   11,   52,   28, 2002, &
          1,   21,  154,   33,   77,  182,    1,  143, &
          1,    1,   44,   44,    1,   13,    4,   52, &
          1,   12,   44,   66,   11,   52,   28, 2002, &
          1,    1,    1,    2,    1,    1,   14,    7, &
          6,   21,   77,   66,  154,   91,    1,  858, &
          1,    1,    1,    2,    1,    1,    1,    2, &
          1,    4,   22,   44,   44,   26,    4,  143, &
          1,    4,   22,   44,   44,   26,    4,  143, &
          1,    1,    1,    2,    1,    1,    1,    2, &
          6,   42,   77,   33,  154,  182,    1,  429, &
          1,    1,   11,   22,    1,   65,   10,   13, &
          1,   28,  308,   22,   77, 1820,   20,  286, &
          1,    1,   44,   44,    1,   65,  140,  364, &
          1,    1,   22,   11,    1,  130,   35,   91, &
          1,   21,   77,   66,  231,  455,   15,  286, &
          1,    1,    1,    1,    1,    1,    1,    1, &
          1,    1,    1,    1,   33,   65,   15,  143, &
          1,    1,    1,    1,    1,   65,   35,   91, &
          1,    1,    1,    1,    1,   65,   35,   91, &
          1,   21,   77,   66,  231,  455,   15,  286, &
          1,   84,  154,  132,  924,  910,   60,  143, &
          1,    1,    1,    1,    1,    1,    7,    7, &
          1,    1,    1,    1,   33,   65,  105, 1001, &
          1,    1,    1,    1,   22,   26,    2,  286, &
          1,    1,    1,    1,    1,    1,    1,    1, &
          1,   84,  154,  132,  924,  910,   60,  143, &
          1,    1,   44,   44,    1,   65,  140,  364, &
          1,    1,   22,   11,    1,  130,   35,   91, &
          1,    1,   11,   22,    1,   65,   10,   13, &
          1,   28,  308,   22,   77, 1820,   20,  286, &
          1,   28,  308,   22,   77, 1820,   20,  286, &
          1,    1,   11,   22,    1,   65,   10,   13, &
          1,    1,   22,   11,    1,  130,   35,   91, &
          1,    1,   44,   44,    1,   65,  140,  364, &
          1,   84,  154,  132,  924,  910,   60,  143, &
          1,    1,    1,    1,    1,    1,    1,    1, &
          1,    1,    1,    1,   22,   26,    2,  286, &
          1,    1,    1,    1,   33,   65,  105, 1001, &
          1,    1,    1,    1,    1,    1,    7,    7, &
          1,   84,  154,  132,  924,  910,   60,  143, &
          1,   21,   77,   66,  231,  455,   15,  286, &
          1,    1,    1,    1,    1,   65,   35,   91, &
          1,    1,    1,    1,    1,   65,   35,   91, &
          1,    1,    1,    1,   33,   65,   15,  143, &
          1,    1,    1,    1,    1,    1,    1,    1, &
          1,   21,   77,   66,  231,  455,   15,  286, &
          1,    1,   22,   11,    1,  130,   35,   91, &
          1,    1,   44,   44,    1,   65,  140,  364, &
          1,   28,  308,   22,   77, 1820,   20,  286, &
          1,    1,   11,   22,    1,   65,   10,   13, &
          6,   42,   77,   33,  154,  182,    1,  429, &
          1,    1,    1,    2,    1,    1,    1,    2, &
          1,    4,   22,   44,   44,   26,    4,  143, &
          1,    4,   22,   44,   44,   26,    4,  143, &
          1,    1,    1,    2,    1,    1,    1,    2, &
          6,   21,   77,   66,  154,   91,    1,  858, &
          1,    1,    1,    2,    1,    1,   14,    7, &
          1,   12,   44,   66,   11,   52,   28, 2002, &
          1,    1,   44,   44,    1,   13,    4,   52, &
          1,   21,  154,   33,   77,  182,    1,  143, &
          1,   12,   44,   66,   11,   52,   28, 2002, &
          1,    1,    1,    2,    1,    1,   14,    7, &
          1,   21,  154,   33,   77,  182,    1,  143, &
          1,    1,   44,   44,    1,   13,    4,   52, &
          6,   21,  154,   33,   77,  182,    1,  858  /

! ... subshell [7/2]^5

      Integer, parameter :: kd_j7_q5 = 56
      Integer, parameter :: nt_j7_q5 =  6

      Integer :: Idet_j7_q5(5,kd_j7_q5)
      Integer :: MD_j7_q5(kd_j7_q5)
      Integer :: ID_j7_q5(nt_j7_q5,kd_j7_q5)
      Integer :: JD_j7_q5(nt_j7_q5,kd_j7_q5)

      Data Idet_j7_q5/ &
          1,    2,    3,    4,    5, &
          1,    2,    3,    4,    6, &
          1,    2,    3,    4,    7, &
          1,    2,    3,    4,    8, &
          1,    2,    3,    5,    6, &
          1,    2,    3,    5,    7, &
          1,    2,    3,    5,    8, &
          1,    2,    3,    6,    7, &
          1,    2,    3,    6,    8, &
          1,    2,    3,    7,    8, &
          1,    2,    4,    5,    6, &
          1,    2,    4,    5,    7, &
          1,    2,    4,    5,    8, &
          1,    2,    4,    6,    7, &
          1,    2,    4,    6,    8, &
          1,    2,    4,    7,    8, &
          1,    2,    5,    6,    7, &
          1,    2,    5,    6,    8, &
          1,    2,    5,    7,    8, &
          1,    2,    6,    7,    8, &
          1,    3,    4,    5,    6, &
          1,    3,    4,    5,    7, &
          1,    3,    4,    5,    8, &
          1,    3,    4,    6,    7, &
          1,    3,    4,    6,    8, &
          1,    3,    4,    7,    8, &
          1,    3,    5,    6,    7, &
          1,    3,    5,    6,    8, &
          1,    3,    5,    7,    8, &
          1,    3,    6,    7,    8, &
          1,    4,    5,    6,    7, &
          1,    4,    5,    6,    8, &
          1,    4,    5,    7,    8, &
          1,    4,    6,    7,    8, &
          1,    5,    6,    7,    8, &
          2,    3,    4,    5,    6, &
          2,    3,    4,    5,    7, &
          2,    3,    4,    5,    8, &
          2,    3,    4,    6,    7, &
          2,    3,    4,    6,    8, &
          2,    3,    4,    7,    8, &
          2,    3,    5,    6,    7, &
          2,    3,    5,    6,    8, &
          2,    3,    5,    7,    8, &
          2,    3,    6,    7,    8, &
          2,    4,    5,    6,    7, &
          2,    4,    5,    6,    8, &
          2,    4,    5,    7,    8, &
          2,    4,    6,    7,    8, &
          2,    5,    6,    7,    8, &
          3,    4,    5,    6,    7, &
          3,    4,    5,    6,    8, &
          3,    4,    5,    7,    8, &
          3,    4,    6,    7,    8, &
          3,    5,    6,    7,    8, &
          4,    5,    6,    7,    8  /

      Data MD_j7_q5/ &
         -5,    5,   -7,    7,   -3,  -15,   -1,   -5,    9,   -3, &
          3,   -9,    5,    1,   15,    3,   -7,    7,   -5,    5, &
         -1,  -13,    1,   -3,   11,   -1,  -11,    3,   -9,    1, &
         -5,    9,   -3,    7,   -1,    1,  -11,    3,   -1,   13, &
          1,   -9,    5,   -7,    3,   -3,   11,   -1,    9,    1, &
         -7,    7,   -5,    5,   -3,    3  /

      Data ID_j7_q5/ &
         -1,    0,   -5,   49,   -1,    1, &
         -1,    0,    5,  -49,   -1,    1, &
         -1,    0,    0,  112,  -27,    5, &
         -1,    0,    0, -112,  -27,    5, &
          1,    0,   25, -343,    8,   -5, &
          0,    0,    0,    0,    0,   -1, &
          0,    1,    7,  -81,    5,   -5, &
          0,    0,   35,   63, -100,  192, &
          0,    0,    0,  -28,  -45,   20, &
         -1,    3,   -1, -477,    1,  -20, &
          1,    0,  -25,  343,    8,   -5, &
          0,    0,    0,  -28,   45,  -20, &
          0,    0,   35,   63,  100, -192, &
          0,   -1,    7,  -81,   -5,    5, &
          0,    0,    0,    0,    0,    1, &
         -1,    3,    1,  477,    1,  -20, &
          1,    0,    0,   -7,  -12,   45, &
          1,    0,    0,    7,  -12,   45, &
         -1,    0,    5, -361,   -1,  -36, &
         -1,    0,   -5,  361,   -1,  -36, &
         -1,    7,   -1,   35,   -7,   28, &
          0,    0,    0,    0,    0,    1, &
          0,   -2,   -7,  135,   -2,   64, &
          0,   -1,  -63,   45,   45, -243, &
          0,    0,    0,    0,    3,   -4, &
          1,    8,    9,   55,    0,  175, &
          0,    0,    0,    0,   -4,   -3, &
          0,    3,  -28,  -60,    5, -100, &
          0,    0,    0,   80,    4,    1, &
          0,    1,   98,   90,    7,   28, &
          0,    0,   -7, -315,   -5,  375, &
          0,    0,    0,   35,    4,   64, &
          0,   -6,   49,  -15,  -10,  -64, &
          0,    0,    0,   45,   -5,   -3, &
         -1,   -3,    8,   40,    7,   21, &
         -1,    7,    1,  -35,   -7,   28, &
          0,    0,    0,    0,   -3,    4, &
          0,    1,  -63,   45,  -45,  243, &
          0,    2,   -7,  135,    2,  -64, &
          0,    0,    0,    0,    0,   -1, &
          1,    8,   -9,  -55,    0,  175, &
          0,    0,    0,   35,   -4,  -64, &
          0,    0,   -7, -315,    5, -375, &
          0,    0,    0,   45,    5,    3, &
          0,    6,   49,  -15,   10,   64, &
          0,   -3,  -28,  -60,   -5,  100, &
          0,    0,    0,    0,    4,    3, &
          0,   -1,   98,   90,   -7,  -28, &
          0,    0,    0,   80,   -4,   -1, &
         -1,   -3,   -8,  -40,    7,   21, &
         -1,    0,    0, -175,    3,   20, &
         -1,    0,    0,  175,    3,   20, &
          1,    0,   -5,  -25,   -4,  -25, &
          1,    0,    5,   25,   -4,  -25, &
         -1,   -3,    3,   25,    3,    5, &
         -1,   -3,   -3,  -25,    3,    5  /

      Data JD_j7_q5/ &
          3,    1,   11,  429,   11,  143, &
          3,    1,   11,  429,   11,  143, &
          3,    1,    1,  429,   77,   91, &
          3,    1,    1,  429,   77,   91, &
          3,    1,  132, 1716,   33,  143, &
          1,    1,    1,    1,    1,    1, &
          1,    2,   66,  286,   66,  143, &
          1,    1,  132,  572,  231, 1001, &
          1,    1,    1,  143,   77,   91, &
          3,   14,  132, 1646,   66,  143, &
          3,    1,  132, 1716,   33,  143, &
          1,    1,    1,  143,   77,   91, &
          1,    1,  132,  572,  231, 1001, &
          1,    2,   66,  286,   66,  143, &
          1,    1,    1,    1,    1,    1, &
          3,   14,  132, 1646,   66,  143, &
          3,    1,    1,  429,   77,   91, &
          3,    1,    1,  429,   77,   91, &
          3,    1,   44, 1716,   11,  143, &
          3,    1,   44, 1716,   11,  143, &
          3,   30,  110,  858,   22,  429, &
          1,    1,    1,    1,    1,    1, &
          1,   15,  110,  286,   11,  429, &
          1,   10,  220,  572,  154, 1001, &
          1,    1,    1,    1,    7,    7, &
          3,  105,  110,  546,    1,  429, &
          1,    1,    1,    1,    7,    7, &
          1,   10,  165,  143,  462, 1001, &
          1,    1,    1,  143,   11,   13, &
          1,   70,  165, 1001,   66,  143, &
          1,    1,  132,  572,  231, 1001, &
          1,    1,    1,  143,   77,   91, &
          1,   35,  660, 4004,   33,  143, &
          1,    1,    1,  143,   11,   13, &
          3,   70,   55, 3003,   22,  143, &
          3,   30,  110,  858,   22,  429, &
          1,    1,    1,    1,    7,    7, &
          1,   10,  220,  572,  154, 1001, &
          1,   15,  110,  286,   11,  429, &
          1,    1,    1,    1,    1,    1, &
          3,  105,  110,  546,    1,  429, &
          1,    1,    1,  143,   77,   91, &
          1,    1,  132,  572,  231, 1001, &
          1,    1,    1,  143,   11,   13, &
          1,   35,  660, 4004,   33,  143, &
          1,   10,  165,  143,  462, 1001, &
          1,    1,    1,    1,    7,    7, &
          1,   70,  165, 1001,   66,  143, &
          1,    1,    1,  143,   11,   13, &
          3,   70,   55, 3003,   22,  143, &
          3,    1,    1,  429,   77,   91, &
          3,    1,    1,  429,   77,   91, &
          3,    1,   44, 1716,   11,  143, &
          3,    1,   44, 1716,   11,  143, &
          3,   14,   11, 3003,   22,  143, &
          3,   14,   11, 3003,   22,  143  /

! ... subshell [7/2]^6

      Integer, parameter :: kd_j7_q6 = 28
      Integer, parameter :: nt_j7_q6 =  4

      Integer :: Idet_j7_q6(6,kd_j7_q6)
      Integer :: MD_j7_q6(kd_j7_q6)
      Integer :: ID_j7_q6(nt_j7_q6,kd_j7_q6)
      Integer :: JD_j7_q6(nt_j7_q6,kd_j7_q6)

      Data Idet_j7_q6/ &
          1,    2,    3,    4,    5,    6, &
          1,    2,    3,    4,    5,    7, &
          1,    2,    3,    4,    5,    8, &
          1,    2,    3,    4,    6,    7, &
          1,    2,    3,    4,    6,    8, &
          1,    2,    3,    4,    7,    8, &
          1,    2,    3,    5,    6,    7, &
          1,    2,    3,    5,    6,    8, &
          1,    2,    3,    5,    7,    8, &
          1,    2,    3,    6,    7,    8, &
          1,    2,    4,    5,    6,    7, &
          1,    2,    4,    5,    6,    8, &
          1,    2,    4,    5,    7,    8, &
          1,    2,    4,    6,    7,    8, &
          1,    2,    5,    6,    7,    8, &
          1,    3,    4,    5,    6,    7, &
          1,    3,    4,    5,    6,    8, &
          1,    3,    4,    5,    7,    8, &
          1,    3,    4,    6,    7,    8, &
          1,    3,    5,    6,    7,    8, &
          1,    4,    5,    6,    7,    8, &
          2,    3,    4,    5,    6,    7, &
          2,    3,    4,    5,    6,    8, &
          2,    3,    4,    5,    7,    8, &
          2,    3,    4,    6,    7,    8, &
          2,    3,    5,    6,    7,    8, &
          2,    4,    5,    6,    7,    8, &
          3,    4,    5,    6,    7,    8  /

      Data MD_j7_q6/ &
          0,  -12,    2,   -2,   12,    0,  -10,    4,   -8,    2, &
         -4,   10,   -2,    8,    0,   -8,    6,   -6,    4,   -4, &
          2,   -6,    8,   -4,    6,   -2,    4,    0  /

      Data ID_j7_q6/ &
         -1,    7,    7,    1, &
          0,    0,    0,   -1, &
          0,    1,    5,    1, &
          0,   -1,   -5,   -1, &
          0,    0,    0,    1, &
          1,   -1,  169,   25, &
          0,    0,    0,    1, &
          0,   -1,  -15,   -5, &
          0,    0,   15,   -7, &
          0,   -8,   15,   14, &
          0,    1,   15,    5, &
          0,    0,    0,   -1, &
          0,    8,  -15,  -14, &
          0,    0,  -15,    7, &
         -1,   -3,   -9,   27, &
          0,    0,   -7,  -15, &
          0,    0,    7,    4, &
          0,    0,   -4,    7, &
          0,    5,    1,   -7, &
          0,  -10,   24,   -7, &
          0,    5,  -27,   35, &
          0,    0,   -7,   -4, &
          0,    0,    7,   15, &
          0,   -5,   -1,    7, &
          0,    0,    4,   -7, &
          0,   -5,   27,  -35, &
          0,   10,  -24,    7, &
          1,   25,  -81,   25  /

      Data JD_j7_q6/ &
          4,   12,   44,  132, &
          1,    1,    1,    1, &
          1,    2,   11,   22, &
          1,    2,   11,   22, &
          1,    1,    1,    1, &
          4,   84,  308,  132, &
          1,    1,    1,    1, &
          1,    6,   22,   33, &
          1,    1,   22,   22, &
          1,   21,   77,   33, &
          1,    6,   22,   33, &
          1,    1,    1,    1, &
          1,   21,   77,   33, &
          1,    1,   22,   22, &
          4,   28,  308,   44, &
          1,    1,   22,   22, &
          1,    1,   11,   11, &
          1,    1,   11,   11, &
          1,   14,  154,   11, &
          1,   21,   77,   33, &
          1,   42,   77,   66, &
          1,    1,   11,   11, &
          1,    1,   22,   22, &
          1,   14,  154,   11, &
          1,    1,   11,   11, &
          1,   42,   77,   66, &
          1,   21,   77,   33, &
          4,   84,  308,  132  /

! ... subshell [7/2]^7

      Integer, parameter :: kd_j7_q7 =  8
      Integer, parameter :: nt_j7_q7 =  1

      Integer :: Idet_j7_q7(7,kd_j7_q7)
      Integer :: MD_j7_q7(kd_j7_q7)
      Integer :: ID_j7_q7(nt_j7_q7,kd_j7_q7)
      Integer :: JD_j7_q7(nt_j7_q7,kd_j7_q7)

      Data Idet_j7_q7/ &
          1,    2,    3,    4,    5,    6,    7, &
          1,    2,    3,    4,    5,    6,    8, &
          1,    2,    3,    4,    5,    7,    8, &
          1,    2,    3,    4,    6,    7,    8, &
          1,    2,    3,    5,    6,    7,    8, &
          1,    2,    4,    5,    6,    7,    8, &
          1,    3,    4,    5,    6,    7,    8, &
          2,    3,    4,    5,    6,    7,    8  /

      Data MD_j7_q7/ &
         -7,    7,   -5,    5,   -3,    3,   -1,    1  /

      Data ID_j7_q7/ &
         -1, &
         -1, &
          1, &
          1, &
         -1, &
         -1, &
          1, &
          1  /

      Data JD_j7_q7/ &
          1, &
          1, &
          1, &
          1, &
          1, &
          1, &
          1, &
          1  /

! ... subshell [7/2]^8

      Integer, parameter :: kd_j7_q8 =  1
      Integer, parameter :: nt_j7_q8 =  1

      Integer :: Idet_j7_q8(8,kd_j7_q8)
      Integer :: MD_j7_q8(kd_j7_q8)
      Integer :: ID_j7_q8(nt_j7_q8,kd_j7_q8)
      Integer :: JD_j7_q8(nt_j7_q8,kd_j7_q8)

      Data Idet_j7_q8/ &
          1,    2,    3,    4,    5,    6,    7,    8  /

      Data MD_j7_q8/ &
          0  /

      Data ID_j7_q8/ &
          1  /

      Data JD_j7_q8/ &
          1  /

! ... subshell [9/2]^2

      Integer, parameter :: kd_j9_q2 = 45
      Integer, parameter :: nt_j9_q2 =  5

      Integer :: Idet_j9_q2(2,kd_j9_q2)
      Integer :: MD_j9_q2(kd_j9_q2)
      Integer :: ID_j9_q2(nt_j9_q2,kd_j9_q2)
      Integer :: JD_j9_q2(nt_j9_q2,kd_j9_q2)

      Data Idet_j9_q2/ &
          1,    2, &
          1,    3, &
          1,    4, &
          1,    5, &
          1,    6, &
          1,    7, &
          1,    8, &
          1,    9, &
          1,   10, &
          2,    3, &
          2,    4, &
          2,    5, &
          2,    6, &
          2,    7, &
          2,    8, &
          2,    9, &
          2,   10, &
          3,    4, &
          3,    5, &
          3,    6, &
          3,    7, &
          3,    8, &
          3,    9, &
          3,   10, &
          4,    5, &
          4,    6, &
          4,    7, &
          4,    8, &
          4,    9, &
          4,   10, &
          5,    6, &
          5,    7, &
          5,    8, &
          5,    9, &
          5,   10, &
          6,    7, &
          6,    8, &
          6,    9, &
          6,   10, &
          7,    8, &
          7,    9, &
          7,   10, &
          8,    9, &
          8,   10, &
          9,   10  /

      Data MD_j9_q2/ &
          0,   -4,    2,   -6,    4,   -8,    6,  -10,    8,   -2, &
          4,   -4,    6,   -6,    8,   -8,   10,    0,   -8,    2, &
        -10,    4,  -12,    6,   -2,    8,   -4,   10,   -6,   12, &
          0,  -12,    2,  -14,    4,   -2,   12,   -4,   14,    0, &
        -16,    2,   -2,   16,    0  /

      Data ID_j9_q2/ &
         -1,    8, -162,   32,  -98, &
          0,  -25,   75,   -7,   21, &
          0,    2,  -27,   56, -294, &
          0,    0,   25,   -4,    6, &
          0,   -7,    7,    1,  -81, &
          0,    0,  -50,   -3,   81, &
          0,    0,   64,   -1,   -6, &
          0,    0,    0,    3,    2, &
          0,    0,  -18,  -15,   -5, &
          0,   -2,   27,  -56,  294, &
          0,   25,  -75,    7,  -21, &
          0,    7,   -7,   -1,   81, &
          0,    0,  -25,    4,   -6, &
          0,    0,  -64,    1,    6, &
          0,    0,   50,    3,  -81, &
          0,    0,   18,   15,    5, &
          0,    0,    0,   -3,   -2, &
          1,   -3,    9,    6, -392, &
          0,    0,   75,  -16,   12, &
          0,   -7,   42,   -4, -336, &
          0,    0,    0,   -2,    3, &
          0,    7,   21,  -25,  -75, &
          0,    0,    0,    3,    7, &
          0,    0,  -54,   -6,   -1, &
          0,    7,  -42,    4,  336, &
          0,    0,  -75,   16,  -12, &
          0,   -7,  -21,   25,   75, &
          0,    0,    0,    2,   -3, &
          0,    0,   54,    6,    1, &
          0,    0,    0,   -3,   -7, &
         -1,    1,  289,  -10,  -40, &
          0,    0,    0,   -7,    3, &
          0,    4,   -2,  -28,  -81, &
          0,    0,    0,    0,    1, &
          0,   -1,  -81,   -7,   -7, &
          0,   -4,    2,   28,   81, &
          0,    0,    0,    7,   -3, &
          0,    1,   81,    7,    7, &
          0,    0,    0,    0,   -1, &
          1,    2,  -22,  -11,  -49, &
          0,    0,    0,    0,    1, &
          0,   -4,  -72,   -7,   -4, &
          0,    4,   72,    7,    4, &
          0,    0,    0,    0,   -1, &
         -1,   -6, -162,   -3,    0  /

      Data JD_j9_q2/ &
          5,   33,  715,  165,  715, &
          1,   66,  286,   33,  143, &
          1,   33,  143,  165,  715, &
          1,    1,  143,   11,   13, &
          1,   22,  286,   11,  143, &
          1,    1,  143,  110,  130, &
          1,    1,  143,   11,   13, &
          1,    1,    1,    5,    5, &
          1,    1,  143,   22,   26, &
          1,   33,  143,  165,  715, &
          1,   66,  286,   33,  143, &
          1,   22,  286,   11,  143, &
          1,    1,  143,   11,   13, &
          1,    1,  143,   11,   13, &
          1,    1,  143,  110,  130, &
          1,    1,  143,   22,   26, &
          1,    1,    1,    5,    5, &
          5,   22, 1430,   55,  715, &
          1,    1,  143,   55,   65, &
          1,   33,  143,  165,  715, &
          1,    1,    1,    5,    5, &
          1,   33,  143,   66,  286, &
          1,    1,    1,   10,   10, &
          1,    1,  143,   11,   13, &
          1,   33,  143,  165,  715, &
          1,    1,  143,   55,   65, &
          1,   33,  143,   66,  286, &
          1,    1,    1,    5,    5, &
          1,    1,  143,   11,   13, &
          1,    1,    1,   10,   10, &
          5,   66, 1430,   33,  143, &
          1,    1,    1,   10,   10, &
          1,   11,  143,   55,  715, &
          1,    1,    1,    1,    1, &
          1,   11,  143,   22,  286, &
          1,   11,  143,   55,  715, &
          1,    1,    1,   10,   10, &
          1,   11,  143,   22,  286, &
          1,    1,    1,    1,    1, &
          5,   33,   65,   30, 1430, &
          1,    1,    1,    1,    1, &
          1,   11,  143,   55,  715, &
          1,   11,  143,   55,  715, &
          1,    1,    1,    1,    1, &
          5,   11,  715,  110,    1  /

! ... subshell [11/2]^2

      Integer, parameter :: kd_j11_q2 = 66
      Integer, parameter :: nt_j11_q2 =  6

      Integer :: Idet_j11_q2(2,kd_j11_q2)
      Integer :: MD_j11_q2(kd_j11_q2)
      Integer :: ID_j11_q2(nt_j11_q2,kd_j11_q2)
      Integer :: JD_j11_q2(nt_j11_q2,kd_j11_q2)

      Data Idet_j11_q2/ &
          1,    2, &
          1,    3, &
          1,    4, &
          1,    5, &
          1,    6, &
          1,    7, &
          1,    8, &
          1,    9, &
          1,   10, &
          1,   11, &
          1,   12, &
          2,    3, &
          2,    4, &
          2,    5, &
          2,    6, &
          2,    7, &
          2,    8, &
          2,    9, &
          2,   10, &
          2,   11, &
          2,   12, &
          3,    4, &
          3,    5, &
          3,    6, &
          3,    7, &
          3,    8, &
          3,    9, &
          3,   10, &
          3,   11, &
          3,   12, &
          4,    5, &
          4,    6, &
          4,    7, &
          4,    8, &
          4,    9, &
          4,   10, &
          4,   11, &
          4,   12, &
          5,    6, &
          5,    7, &
          5,    8, &
          5,    9, &
          5,   10, &
          5,   11, &
          5,   12, &
          6,    7, &
          6,    8, &
          6,    9, &
          6,   10, &
          6,   11, &
          6,   12, &
          7,    8, &
          7,    9, &
          7,   10, &
          7,   11, &
          7,   12, &
          8,    9, &
          8,   10, &
          8,   11, &
          8,   12, &
          9,   10, &
          9,   11, &
          9,   12, &
         10,   11, &
         10,   12, &
         11,   12  /

      Data MD_j11_q2/ &
          0,   -4,    2,   -6,    4,   -8,    6,  -10,    8,  -12, &
         10,   -2,    4,   -4,    6,   -6,    8,   -8,   10,  -10, &
         12,    0,   -8,    2,  -10,    4,  -12,    6,  -14,    8, &
         -2,    8,   -4,   10,   -6,   12,   -8,   14,    0,  -12, &
          2,  -14,    4,  -16,    6,   -2,   12,   -4,   14,   -6, &
         16,    0,  -16,    2,  -18,    4,   -2,   16,   -4,   18, &
          0,  -20,    2,   -2,   20,    0  /

      Data ID_j11_q2/ &
          1, -175,   28, -100,  991, -441, &
          0,   45,  -32,   36, -432,  462, &
          0,   -5,   16,  -40,  840, -879, &
          0,    0,  -14,   81,  -81,  231, &
          0,   40,   -9,   -1,  867, -862, &
          0,    0,  189,   -5, -125,  693, &
          0,    0,  -42,  128,   -1, -154, &
          0,    0,    0,   -8,    1,  154, &
          0,    0,  105,  169, -289, -385, &
          0,    0,    0,   11,   77,   70, &
          0,    0,    0,  -11,  -11,  -63, &
          0,    5,  -16,   40, -840,  879, &
          0,  -45,   32,  -36,  432, -462, &
          0,  -40,    9,    1, -867,  862, &
          0,    0,   14,  -81,   81, -231, &
          0,    0,   42, -128,    1,  154, &
          0,    0, -189,    5,  125, -693, &
          0,    0, -105, -169,  289,  385, &
          0,    0,    0,    8,   -1, -154, &
          0,    0,    0,   11,   11,   63, &
          0,    0,    0,  -11,  -77,  -70, &
         -1,  841,  -36,   -4,  804, -299, &
          0,    0, -245,  189, -189,  165, &
          0,  128, -250,   28,    3, -961, &
          0,    0,    0,    7,   -7,  275, &
          0, -216,  -15,  105, -875, -799, &
          0,    0,    0,  -35,   -5,  198, &
          0,    0,   54,    7, -343,  -99, &
          0,    0,    0,    0,   11,    8, &
          0,    0,   -3,  -35, -385,  -49, &
          0, -128,  250,  -28,   -3,  961, &
          0,    0,  245, -189,  189, -165, &
          0,  216,   15, -105,  875,  799, &
          0,    0,    0,   -7,    7, -275, &
          0,    0,  -54,   -7,  343,   99, &
          0,    0,    0,   35,    5, -198, &
          0,    0,    3,   35,  385,   49, &
          0,    0,    0,    0,  -11,   -8, &
          1, -289,  -13,  625, -325, -356, &
          0,    0,    0,   28,  -16,   55, &
          0, -243,  135,   21, -889, -757, &
          0,    0,    0,    0,   -8,   11, &
          0,  135,   39,  -21,  -91, -539, &
          0,    0,    0,    0,   11,   27, &
          0,    0,   -3,  -28, -154,   -4, &
          0,  243, -135,  -21,  889,  757, &
          0,    0,    0,  -28,   16,  -55, &
          0, -135,  -39,   21,   91,  539, &
          0,    0,    0,    0,    8,  -11, &
          0,    0,    3,   28,  154,    4, &
          0,    0,    0,    0,  -11,  -27, &
         -1,    0,   99,  -11, -420, -481, &
          0,    0,    0,    0,  -27,   11, &
          0,  320,    9,  -70, -750,  -88, &
          0,    0,    0,    0,    0,    1, &
          0,   -5,  -81,   -7,  -21,  -27, &
          0, -320,   -9,   70,  750,   88, &
          0,    0,    0,    0,   27,  -11, &
          0,    5,   81,    7,   21,   27, &
          0,    0,    0,    0,    0,   -1, &
          1,  625, -729, -961, -185,  -67, &
          0,    0,    0,    0,    0,    1, &
          0,  -25,  -45,   -7,   -6,    0, &
          0,   25,   45,    7,    6,    0, &
          0,    0,    0,    0,    0,   -1, &
         -1, -275,  -99,  -11,  -11,    0  /

      Data JD_j11_q2/ &
          6,  858,  143,  561, 6594, 4199, &
          1,  143,  143,  187, 2717, 4199, &
          1,  143,  143,  187, 2717, 2663, &
          1,    1,  143,  374,  247,  646, &
          1,  143,  143,  374, 5434, 1741, &
          1,    1,  572,  748,  988, 1292, &
          1,    1,  143,  561,  741,  323, &
          1,    1,    1,   17,   19,  323, &
          1,    1,  572,  748,  988, 1292, &
          1,    1,    1,  102,  114,  323, &
          1,    1,    1,   34,   19,  646, &
          1,  143,  143,  187, 2717, 2663, &
          1,  143,  143,  187, 2717, 4199, &
          1,  143,  143,  374, 5434, 1741, &
          1,    1,  143,  374,  247,  646, &
          1,    1,  143,  561,  741,  323, &
          1,    1,  572,  748,  988, 1292, &
          1,    1,  572,  748,  988, 1292, &
          1,    1,    1,   17,   19,  323, &
          1,    1,    1,   34,   19,  646, &
          1,    1,    1,  102,  114,  323, &
          6, 6006, 1001,  561, 4787,  620, &
          1,    1,  572,  748,  988, 1292, &
          1, 1001, 1001,  187, 2717, 2038, &
          1,    1,    1,   34,   19,  646, &
          1, 1001, 1001,  374, 5434, 2440, &
          1,    1,    1,  102,  114,  323, &
          1,    1,  143, 1122,  741,  646, &
          1,    1,    1,    1,   19,   19, &
          1,    1,   52,   68,  988, 1292, &
          1, 1001, 1001,  187, 2717, 2038, &
          1,    1,  572,  748,  988, 1292, &
          1, 1001, 1001,  374, 5434, 2440, &
          1,    1,    1,   34,   19,  646, &
          1,    1,  143, 1122,  741,  646, &
          1,    1,    1,  102,  114,  323, &
          1,    1,   52,   68,  988, 1292, &
          1,    1,    1,    1,   19,   19, &
          6, 6006,  308, 2244, 2508, 1063, &
          1,    1,    1,   51,   57,  323, &
          1, 1001, 1001,  374, 2283, 4281, &
          1,    1,    1,    1,   19,   19, &
          1, 1001,  154,  187,  209, 8398, &
          1,    1,    1,    1,   38,   38, &
          1,    1,   13,   51,  741,  323, &
          1, 1001, 1001,  374, 2283, 4281, &
          1,    1,    1,   51,   57,  323, &
          1, 1001,  154,  187,  209, 8398, &
          1,    1,    1,    1,   19,   19, &
          1,    1,   13,   51,  741,  323, &
          1,    1,    1,    1,   38,   38, &
          6,    1,  364,  204,  967, 6595, &
          1,    1,    1,    1,   38,   38, &
          1, 1001, 1001,  187, 2717, 4199, &
          1,    1,    1,    1,    1,    1, &
          1,   91,  182,   17,  247, 8398, &
          1, 1001, 1001,  187, 2717, 4199, &
          1,    1,    1,    1,   38,   38, &
          1,   91,  182,   17,  247, 8398, &
          1,    1,    1,    1,    1,    1, &
          6, 6006, 4004, 2244, 1621,13893, &
          1,    1,    1,    1,    1,    1, &
          1,   91,   91,   34,  247,    1, &
          1,   91,   91,   34,  247,    1, &
          1,    1,    1,    1,    1,    1, &
          6,  546,  364,  204, 2964,    1  /

      End Module det_jq


!====================================================================
      Subroutine  DET_sh_jq (j,q,id,MJ,Idet)
!====================================================================
!     gives the determinant id and MJ for the subshell j^q
!     Input:  jq - subshell
!             id - the determinant pointer
!     Output: Idet - description of the determinant as
!                    list of orbital positions in the common list
!             MJ - total azimutal quantum numbers (2*m+1)
!     Calls:  MJL_id, MS_id
!-----------------------------------------------------------------------
      Use det_jq

      Implicit none
      Integer, intent(in) :: j,q,id
      Integer, intent(out) :: MJ
      Integer, intent(out) :: Idet(q)
      Integer :: i,ii,k,nd,nq
      Integer, external :: MJ_value

      if(id.le.0.or.j.lt.0) Call DET_stop (j,q,id)
!-----------------------------------------------------------------------
! ... two-electron case for j>11, the only case we should compute:

      if(q.eq.2.and.j.gt.11) then
       nq=j+1; nd=nq*(nq-1)/2; if(id.gt.nd) Call DET_stop(j,q,id)
       k=0
       Do i=1,nq-1
        Do ii=i+1,nq
         k=k+1; if(k.eq.id) Exit
        End do
        if(k.eq.id) Exit
       End do
       Idet(1)=i; Idet(2)=j
       MJ = mj_value(id)
       Return
!----------------------------------------------------------------------

      elseif(q.eq.1) then                             ! j^1 - subshell:

       if(id.gt.j+1) Call DET_stop (j,q,id)
       MJ=MJ_value(id)
       Idet(1)=id
       Return

      else

      SELECT CASE(j*100+q)

      CASE(102)                                   !  [1/2]^2 - subshell
       if(id.gt.kd_j1_q2) Call DET_stop(j,q,id)
       MJ = MD_j1_q2(id)
       Idet(1:q)= Idet_j1_q2(1:q,id)
      CASE(302)                                   !  [3/2]^2 - subshell
       if(id.gt.kd_j3_q2) Call DET_stop(j,q,id)
       MJ = MD_j3_q2(id)
       Idet(1:q)= Idet_j3_q2(1:q,id)
      CASE(303)                                   !  [3/2]^3 - subshell
       if(id.gt.kd_j3_q3) Call DET_stop(j,q,id)
       MJ = MD_j3_q3(id)
       Idet(1:q)= Idet_j3_q3(1:q,id)
      CASE(304)                                   !  [3/2]^4 - subshell
       if(id.gt.kd_j3_q4) Call DET_stop(j,q,id)
       MJ = MD_j3_q4(id)
       Idet(1:q)= Idet_j3_q4(1:q,id)
      CASE(502)                                   !  [5/2]^2 - subshell
       if(id.gt.kd_j5_q2) Call DET_stop(j,q,id)
       MJ = MD_j5_q2(id)
       Idet(1:q)= Idet_j5_q2(1:q,id)
      CASE(503)                                   !  [5/2]^3 - subshell
       if(id.gt.kd_j5_q3) Call DET_stop(j,q,id)
       MJ = MD_j5_q3(id)
       Idet(1:q)= Idet_j5_q3(1:q,id)
      CASE(504)                                   !  [5/2]^4 - subshell
       if(id.gt.kd_j5_q4) Call DET_stop(j,q,id)
       MJ = MD_j5_q4(id)
       Idet(1:q)= Idet_j5_q4(1:q,id)
      CASE(505)                                   !  [5/2]^5 - subshell
       if(id.gt.kd_j5_q5) Call DET_stop(j,q,id)
       MJ = MD_j5_q5(id)
       Idet(1:q)= Idet_j5_q5(1:q,id)
      CASE(506)                                   !  [5/2]^6 - subshell
       if(id.gt.kd_j5_q6) Call DET_stop(j,q,id)
       MJ = MD_j5_q6(id)
       Idet(1:q)= Idet_j5_q6(1:q,id)
      CASE(702)                                   !  [7/2]^2 - subshell
       if(id.gt.kd_j7_q2) Call DET_stop(j,q,id)
       MJ = MD_j7_q2(id)
       Idet(1:q)= Idet_j7_q2(1:q,id)
      CASE(703)                                   !  [7/2]^3 - subshell
       if(id.gt.kd_j7_q3) Call DET_stop(j,q,id)
       MJ = MD_j7_q3(id)
       Idet(1:q)= Idet_j7_q3(1:q,id)
      CASE(704)                                   !  [7/2]^4 - subshell
       if(id.gt.kd_j7_q4) Call DET_stop(j,q,id)
       MJ = MD_j7_q4(id)
       Idet(1:q)= Idet_j7_q4(1:q,id)
      CASE(705)                                   !  [7/2]^5 - subshell
       if(id.gt.kd_j7_q5) Call DET_stop(j,q,id)
       MJ = MD_j7_q5(id)
       Idet(1:q)= Idet_j7_q5(1:q,id)
      CASE(706)                                   !  [7/2]^6 - subshell
       if(id.gt.kd_j7_q6) Call DET_stop(j,q,id)
       MJ = MD_j7_q6(id)
       Idet(1:q)= Idet_j7_q6(1:q,id)
      CASE(707)                                   !  [7/2]^7 - subshell
       if(id.gt.kd_j7_q7) Call DET_stop(j,q,id)
       MJ = MD_j7_q7(id)
       Idet(1:q)= Idet_j7_q7(1:q,id)
      CASE(708)                                   !  [7/2]^8 - subshell
       if(id.gt.kd_j7_q8) Call DET_stop(j,q,id)
       MJ = MD_j7_q8(id)
       Idet(1:q)= Idet_j7_q8(1:q,id)
      CASE(902)                                   !  [9/2]^2 - subshell
       if(id.gt.kd_j9_q2) Call DET_stop(j,q,id)
       MJ = MD_j9_q2(id)
       Idet(1:q)= Idet_j9_q2(1:q,id)
      CASE(1102)                                  ! [11/2]^2 - subshell
       if(id.gt.kd_j11_q2) Call DET_stop(j,q,id)
       MJ = MD_j11_q2(id)
       Idet(1:q)= Idet_j11_q2(1:q,id)
      Case default
        Call DET_stop(j,q,id)
      END SELECT

      end if

      End  Subroutine  DET_sh_jq


!======================================================================
      Subroutine DET_stop (j,q,id)
!======================================================================
!     stops execution of DET_sh_jq with debug information
!----------------------------------------------------------------------
      Integer, intent(in) :: j,q,id
      write(*,'(a,a,3i5)') 'DET_sh_jq: id is out the range: ',  &
                           'j,q,id = ',j,q,id
      Stop  'DET_sh_jq: id is out the range:'
      End Subroutine DET_stop


!======================================================================
      Real(8) Function  DETC_jq (j,q,it,id)
!======================================================================
!     gives the expansion coefficient for the given subshell 'jq',
!     term 'it' and determinant 'id'
!     Input:  jq - subshell
!             it - the term number
!             id - the determinant number
!     Calls:  CLEBSH2, MJ_value
!----------------------------------------------------------------------
      Use det_jq

      Implicit none
      Integer, Intent(in) :: j,q,it,id
      Integer :: i,i1,i2,k,mj,mj1,mj2,JT,JV,JW,JQ,nd,nq
      Integer, external :: MJ_value, Jterm
      Real(8) :: C
      Real(8), external :: Clebsh2

      if(it.le.0.or.id.le.0) Call DETC_stop(j,q,it,id)
!----------------------------------------------------------------------
! ... two-electron case for j > 11/2, the only case we should compute:

      if(q.eq.2.and.j.gt.11) then

       if(it.gt.(j+1)/2) Call DETC_stop(j,q,it,id)
       nq=j+1; nd=nq*(nq-1)/2; if(id.gt.nd) Call DET_stop(j,q,id)
       if(id.gt.nd) Call DETC_stop(j,q,it,id)
       k=0
       Do i1=1,nq-1;
        Do i2=i1+1,nq
         k=k+1
         if(k.eq.id) Exit
        End do
        if(k.eq.id) Exit
       End do
       mj1 = MJ_value(i1);  mj2 = MJ_value(i2); MJ = mj1+mj2

       i = Jterm(j,q,it,JT,JV,JW,JQ)                   !  ???
       DETC_jq = CLEBSH2(j,mj1,j,mj2,JT,MJ)
       Return

!----------------------------------------------------------------------

      elseif(q.eq.1) then                         !  one-electron shell

       if(it.ne.1.or.id.gt.j+1) Call DETC_stop(j,q,it,id)
       DETC_jq = 1.d0
       Return

      else

      SELECT CASE(j*100+q)

      CASE(102)                                    ! [1/2]^2 - subshell
       if(it.gt.nt_j1_q2) Call DETC_stop(j,q,it,id)
       if(id.gt.kd_j1_q2) Call DETC_stop(j,q,it,id)
       i1=ID_j1_q2(it,id)
       i2=JD_j1_q2(it,id)
      CASE(302)                                    ! [3/2]^2 - subshell
       if(it.gt.nt_j3_q2) Call DETC_stop(j,q,it,id)
       if(id.gt.kd_j3_q2) Call DETC_stop(j,q,it,id)
       i1=ID_j3_q2(it,id)
       i2=JD_j3_q2(it,id)
      CASE(303)                                    ! [3/2]^3 - subshell
       if(it.gt.nt_j3_q3) Call DETC_stop(j,q,it,id)
       if(id.gt.kd_j3_q3) Call DETC_stop(j,q,it,id)
       i1=ID_j3_q3(it,id)
       i2=JD_j3_q3(it,id)
      CASE(304)                                    ! [3/2]^4 - subshell
       if(it.gt.nt_j3_q4) Call DETC_stop(j,q,it,id)
       if(id.gt.kd_j3_q4) Call DETC_stop(j,q,it,id)
       i1=ID_j3_q4(it,id)
       i2=JD_j3_q4(it,id)
      CASE(502)                                    ! [5/2]^2 - subshell
       if(it.gt.nt_j5_q2) Call DETC_stop(j,q,it,id)
       if(id.gt.kd_j5_q2) Call DETC_stop(j,q,it,id)
       i1=ID_j5_q2(it,id)
       i2=JD_j5_q2(it,id)
      CASE(503)                                    ! [5/2]^3 - subshell
       if(it.gt.nt_j5_q3) Call DETC_stop(j,q,it,id)
       if(id.gt.kd_j5_q3) Call DETC_stop(j,q,it,id)
       i1=ID_j5_q3(it,id)
       i2=JD_j5_q3(it,id)
      CASE(504)                                    ! [5/2]^4 - subshell
       if(it.gt.nt_j5_q4) Call DETC_stop(j,q,it,id)
       if(id.gt.kd_j5_q4) Call DETC_stop(j,q,it,id)
       i1=ID_j5_q4(it,id)
       i2=JD_j5_q4(it,id)
      CASE(505)                                    ! [5/2]^5 - subshell
       if(it.gt.nt_j5_q5) Call DETC_stop(j,q,it,id)
       if(id.gt.kd_j5_q5) Call DETC_stop(j,q,it,id)
       i1=ID_j5_q5(it,id)
       i2=JD_j5_q5(it,id)
      CASE(506)                                    ! [5/2]^6 - subshell
       if(it.gt.nt_j5_q6) Call DETC_stop(j,q,it,id)
       if(id.gt.kd_j5_q6) Call DETC_stop(j,q,it,id)
       i1=ID_j5_q6(it,id)
       i2=JD_j5_q6(it,id)
      CASE(702)                                    ! [7/2]^2 - subshell
       if(it.gt.nt_j7_q2) Call DETC_stop(j,q,it,id)
       if(id.gt.kd_j7_q2) Call DETC_stop(j,q,it,id)
       i1=ID_j7_q2(it,id)
       i2=JD_j7_q2(it,id)
      CASE(703)                                    ! [7/2]^3 - subshell
       if(it.gt.nt_j7_q3) Call DETC_stop(j,q,it,id)
       if(id.gt.kd_j7_q3) Call DETC_stop(j,q,it,id)
       i1=ID_j7_q3(it,id)
       i2=JD_j7_q3(it,id)
      CASE(704)                                    ! [7/2]^4 - subshell
       if(it.gt.nt_j7_q4) Call DETC_stop(j,q,it,id)
       if(id.gt.kd_j7_q4) Call DETC_stop(j,q,it,id)
       i1=ID_j7_q4(it,id)
       i2=JD_j7_q4(it,id)
      CASE(705)                                    ! [7/2]^5 - subshell
       if(it.gt.nt_j7_q5) Call DETC_stop(j,q,it,id)
       if(id.gt.kd_j7_q5) Call DETC_stop(j,q,it,id)
       i1=ID_j7_q5(it,id)
       i2=JD_j7_q5(it,id)
      CASE(706)                                    ! [7/2]^6 - subshell
       if(it.gt.nt_j7_q6) Call DETC_stop(j,q,it,id)
       if(id.gt.kd_j7_q6) Call DETC_stop(j,q,it,id)
       i1=ID_j7_q6(it,id)
       i2=JD_j7_q6(it,id)
      CASE(707)                                    ! [7/2]^7 - subshell
       if(it.gt.nt_j7_q7) Call DETC_stop(j,q,it,id)
       if(id.gt.kd_j7_q7) Call DETC_stop(j,q,it,id)
       i1=ID_j7_q7(it,id)
       i2=JD_j7_q7(it,id)
      CASE(708)                                    ! [7/2]^8 - subshell
       if(it.gt.nt_j7_q8) Call DETC_stop(j,q,it,id)
       if(id.gt.kd_j7_q8) Call DETC_stop(j,q,it,id)
       i1=ID_j7_q8(it,id)
       i2=JD_j7_q8(it,id)
      CASE(902)                                    ! [9/2]^2 - subshell
       if(it.gt.nt_j9_q2) Call DETC_stop(j,q,it,id)
       if(id.gt.kd_j9_q2) Call DETC_stop(j,q,it,id)
       i1=ID_j9_q2(it,id)
       i2=JD_j9_q2(it,id)
      CASE(1102)                                    ![11/2]^2 - subshell
       if(it.gt.nt_j11_q2) Call DETC_stop(j,q,it,id)
       if(id.gt.kd_j11_q2) Call DETC_stop(j,q,it,id)
       i1=ID_j11_q2(it,id)
       i2=JD_j11_q2(it,id)

      CASE DEFAULT

       write(*,*) ' j,q = ',j,q
       Stop ' DETC_jq: jq - outside of possible values'

      END SELECT
      end if

      C=DBLE(i1)/DBLE(i2); C=sqrt(abs(C)); if(i1.lt.0) C=-C; DETC_jq=C

      End  Function  DETC_jq

!======================================================================
      Subroutine DETC_stop (j,q,it,id)
!======================================================================
!     stops execution of DETC_jq with debug information
!----------------------------------------------------------------------
      Integer, intent(in) :: j,q,it,id
      write(6,'(a,a,4i5)') 'DETC_jq: it(id) is out the range: ',  &
                           'j,q,it,id = ',j,q,it,id
      Stop  'DETC_jq: it(id) is out the range:'
      End Subroutine DETC_stop




!=====================================================================
      Real(8) Function cfp_jj(j,q,JP,VP,JD,VD)
!=====================================================================
!     fractional parentage coefficients for j^q subshells
!     j  - 2*j value
!     q  - the number of electrons in daughter subshell
!     JP - 2*J value of parent subshell
!     VP - seniority of parent term
!     JD - 2*J value of daughter subshell
!     VD - seniority of daughter term
!     coefficents are stored in form: +- sqrt(i1/i2)
!     subshell terms are defined in routine Jterm
!---------------------------------------------------------------------
      Implicit none

      Integer :: j,q,JP,VP,JD,VD, WP,QP, WD,QD, IP,ID
      Real(8) :: CN,CD,phase
      Real(8) :: zero = 0.d0, one = 1.d0
      Integer, external :: Jterm

      Integer, parameter :: n2_j3=2, n3_j3=1
      Integer :: num3_j3(n3_j3,n2_j3), norm3_j3(n3_j3)

      Integer, parameter :: n2_j5=3, n3_j5=3
      Integer :: num3_j5(n3_j5,n2_j5), norm3_j5(n3_j5)

      Integer, parameter :: n2_j7=4, n3_j7=6, n4_j7=8
      Integer :: num3_j7(n3_j7,n2_j7), norm3_j7(n3_j7)
      Integer :: num4_j7(n4_j7,n3_j7), norm4_j7(n4_j7)

! ... check of j and q values

      if(q.lt.1.or.q.gt.j+1) then
       write(*,*) 'cfp_jj: q is out of scope:',q
       Stop 'Stop in cfp_jj'
      end if

      if(j.lt.1.or.(j.ge.9.and.q.gt.2)) then
       write(*,*) 'cfp_jj: j is out of scope:',j
       Stop 'Stop in cfp_jj'
      end if

! ... term indexes:

      IP = Jterm (j,q-1,0,JP,VP,WP,QP)
      ID = Jterm (j,q,  0,JD,VD,WD,QD)

! .... trivial cases:

      if(q.le.2.or.q.eq.j+1) then
       cfp_jj = one
       Return
      end if

! ... select different subshells:

      phase = one
      Select case(j*100+q)
      Case(303)                                  !  3/2 ^ 3
       CN = num3_j3(ID,IP)
       CD = norm3_j3(ID)
       if(CN.lt.zero) phase=-phase
      Case(503)                                  !  5/2 ^ 3
       CN = num3_j5(ID,IP)
       CD = norm3_j5(ID)
       if(CN.lt.zero) phase=-phase
      Case(504)                                  !  5/2 ^ 4
       CN = num3_j5(IP,ID)
       CD = norm3_j5(IP)
       if(CN.lt.zero) phase=-phase
       CN = CN * ((7-q)*(1+JP))
       CD = CD * (q*(1+JD))
       phase = phase * (-1)**((JD-JP-VD+VP)/2-3)
      Case(505)                                  !  5/2 ^ 5
       CN = ((7-q)*(1+JP))
       CD = (q*(1+JD))
       phase =  (-1)**((JD-JP-VD+VP)/2-3)
      Case(703)                                  !  7/2 ^ 3
       CN = num3_j7(ID,IP)
       CD = norm3_j7(ID)
       if(CN.lt.zero) phase=-phase
      Case(704)                                  !  7/2 ^ 4
       CN = num4_j7(ID,IP)
       CD = norm4_j7(ID)
       if(CN.lt.zero) phase=-phase
      Case(705)                                  !  7/2 ^ 5
       CN = num4_j7(IP,ID)
       CD = norm4_j7(IP)
       if(CN.lt.zero) phase=-phase
       CN = CN * ((9-q)*(1+JP))
       CD = CD * (q*(1+JD))
       phase = phase * (-1)**((JD-JP-VD+VP)/2-3)
      Case(706)                                  !  7/2 ^ 6
       CN = num3_j7(IP,ID)
       CD = norm3_j7(IP)
       if(CN.lt.zero) phase=-phase
       CN = CN * ((9-q)*(1+JP))
       CD = CD * (q*(1+JD))
       phase = phase * (-1)**((JD-JP-VD+VP)/2-3)
      Case(707)                                  !  7/2 ^ 7
       CN = ((9-q)*(1+JP))
       CD = (q*(1+JD))
       phase =  (-1)**((JD-JP-VD+VP)/2-3)
      Case default
       write(*,*) 'cfp_jj: unknown subshell j^q: ',j,q
       Stop 'Stop in cfp_jj'
      End Select

      cfp_jj = sqrt(abs(CN)/CD) * phase
      Return

      DATA num3_j3 /  1, -5 /
      DATA norm3_j3/  6/

      DATA num3_j5 / -4, 0,  0, &
                      5,-5,  3, &
                      9, 2,-11 /
      DATA norm3_j5/ 18, 7, 14 /

      DATA num3_j7 /  9,   0,   0,   0,   0,   0, &
                     -5,   3, 121, 143, -55,   0, &
                     -9, -11,  12,-900,  39,   5, &
                    -13,   0, -65, 343, 104, -17  /
      DATA norm3_j7/ 36,  14, 198,1386, 198,  22  /

      DATA num4_j7 /  1,  280,  308, 1144,    0,    0,    0,   0, &
                      0,   54, -121,    0, -968,  169,  462,   0, &
                      0, -231,  -14,  195,  -77, 2366, -343,   0, &
                      0,  -65,  250, -245,-1755,   90, -945, 140, &
                      0, -210,   91,  624,  280, 2275,  650, 234, &
                      0,    0,  140,-1224,    0, -560,  680, 627  /
      DATA norm4_j7/  1,  840,  924, 3432, 3080, 5460, 3080,1001  /

      End Function cfp_jj



!======================================================================
      Subroutine coef_1conf(no,ln,jn,iq,Jshell,Vshell,Jintra,kmax,coefs)
!======================================================================
!     compute the angular coefficients for 1 atomic states
!----------------------------------------------------------------------
      Implicit none

! ... input-output:

      Integer, intent(in) :: no,ln(no),jn(no),iq(no), kmax, &
                             Jshell(no),Vshell(no),Jintra(no)
      Real(8), intent(out):: coefs(no,no,0:kmax)

! ... determinant expansion:

      Integer :: kdt
      Integer, allocatable :: IP_det(:,:)
      Real(8), allocatable :: C_det(:)
      Real(8) :: CC_det

! ... local variables:

      Integer             :: ne,Jtotal,i,j,k, jmax,kd1,kd2
      Integer             :: ip1(no),ip2(no)
      Integer, external   :: mj_value

! ... initialize arrays:

      coefs = 0.d0
      ne = SUM(iq(1:no))
      Jtotal = Jintra(no)
      ip1(1)=1; ip2(1)=iq(1)
      Do i=2,no
       ip1(i)=ip2(i-1)+1
       ip2(i)=ip2(i-1)+iq(i)
      End do

! ... determinant expansion:

      Call Det_expn_jj

! ... calculations:       total M - ???  Clebsh - ???

      Do kd1 = 1,kdt
      Do kd2 = kd1,kdt
        CC_det = C_det(kd1) * C_det(kd2)
        if(kd1.ne.kd2) CC_det = CC_det + CC_det
        Call me_det
      End do;  End do


CONTAINS

!======================================================================
      Subroutine Det_expn_jj
!======================================================================
!     procedure of exaustion of all possible determinants for given
!     configurations. The determinants and their coefficients
!     are recoded on unit 'nua'
!
!     Calls: Det_sh_jq, DETC_jq, Clebsh2, Ndets_jq, Jterm, mj_value
!----------------------------------------------------------------------
      Implicit none

      Integer :: i,j,k, mkdt, JW,JQ
      Integer, external :: Ndets_jq, Jterm, mj_value
      Real(8) :: C
      Real(8), external :: DETC_jq, Clebsh2

! ... shell values:
!     md(i)  - the max.number of det.'s for the i-th subshell
!     nd(i)  - determinant under consideration
!     ipn(i) - pointers on the orbitals of given shell
!     MJs    - shells MJ
!     MJi    - intermediate values MJ

      Integer :: md(no),nd(no),ipn(no), MJs(no),MJi(no), ips(no)
      Integer :: Idet(ne)

! ... prepare working arrays:

      k = 1; mkdt = 1
      Do i=1,no
       ipn(i)=k; k=k+iq(i); md(i)=Ndets_jq(jn(i),iq(i))
       mkdt = mkdt*md(i)
       ips(i) = Jterm(jn(i),iq(i),0,Jshell(i),Vshell(i),JW,JQ)
      End do

      if(allocated(C_det)) Deallocate(C_det); Allocate(C_det(mkdt))
      if(allocated(IP_det)) Deallocate(IP_det); Allocate(IP_det(ne,mkdt))

!--------------------------------------------------------------------
! ... exausting all possible determinants:

      kdt=0; i=1; nd(i)=1
    1 Call DET_sh_jq(jn(i),iq(i),nd(i),MJs(i),Idet(ipn(i)))

      if(i.eq.1) then
       MJi(1) = MJs(1)
      else
       MJi(i) = MJi(i-1)+MJs(i)
      end if
      if(i.lt.no) then;  i=i+1;  nd(i)=1;  go to 1; end if

       ! ... coefficient calculation:

       if(MJi(no).ne.Jtotal) go to 2
       C = DETC_jq(jn(1),iq(1),ips(1),nd(1))
       if(C.eq.0.d0) go to 2

       Do j=2,no
        C = C * DETC_jq(jn(j),iq(j),ips(j),nd(j))
        if(C.eq.0.d0) Exit
        C = C * Clebsh2(Jintra(j-1),MJi(j-1), &
                        Jshell(j  ),MJs(j  ), &
                        Jintra(j  ),MJi(j  ))
        if(C.eq.0.d0) Exit
       End do

       if(C.ne.0.d0) then
        kdt=kdt+1; IP_det(1:ne,kdt)=Idet(1:ne); C_det(kdt)=C
       end if

    2 nd(i)=nd(i)+1                ! selecting the next case

      if(nd(i).gt.md(i)) then
       if(i.eq.1) go to 3          ! to end
       i=i-1; go to 2
      end if
      go to 1

    3 Continue

      Do k=1,kdt; Do i=1,ne
       j=IP_det(i,k); IP_det(i,k)=mj_value(j)
      End do; End do

      End Subroutine Det_expn_jj


!======================================================================
      Subroutine me_det
!======================================================================
!     find the possible interaction orbitals in two determinants and
!     call the subroutines to calculate m.e. between possible
!     combinations of nj-orbitals (me_ee)
!----------------------------------------------------------------------
      Implicit none
      Integer :: i,i1,i2, j,j1,j2, k,kk,k1,k2, idif, jdif, io,jo
      Integer :: is,js, isym1,isym2, jsym1,jsym2
      Integer :: ii(ne,ne)

!----------------------------------------------------------------------
! ... the same determinants:

      if(kd1.eq.kd2) then
       Do i=1,no;  Do is = ip1(i),ip2(i)
       Do j=i,no;  Do js = ip1(j),ip2(j)
        if(js.gt.is) Call me_ee(i,j,is,js,is,js)
       End do; End do
       End do; End do
       Return
      end if

!----------------------------------------------------------------------
! ... check total orbitals differenc:

       ii = 0
       Do i=1,no; k=ip1(i); kk=ip2(i)
        Do i1=k,kk; Do i2=k,kk
         if(IP_det(i1,kd1).eq.IP_det(i2,kd2)) ii(i1,i2)=1
        End do; End do
       End do
       idif = ne - SUM(ii)
       if(idif.gt.2) Return
       if(idif.ne.2) Stop 'Det_me: idif <> 2 ???'

! ... find interaction orbitals:

       k = 1;  jdif=0
       Do i=1,no; k=ip1(i); kk=ip2(i); k1=k; k2=k
        idif = iq(i) - SUM(ii(k:kk,k:kk))
        if(idif.eq.0) Cycle

! ... first orbital:

        if(jdif.eq.0) then
         io = i
         Do i1=k,kk
          if(SUM(II(i1,k:kk)).eq.1) Cycle; isym1=i1; Exit
         End do
         Do i2=k,kk
          if(SUM(II(k:kk,i2)).eq.1) Cycle; isym2=i2; Exit
         End do
         jdif = 1
         idif = idif-1
         k1 = isym1+1
         k2 = isym2+1
        end if

	       if(idif.eq.0) Cycle

! ... second orbital:

        jo = i
        Do i1=k1,kk
         if(SUM(II(i1,k:kk)).eq.1) Cycle; jsym1=i1; Exit
        End do
        Do i2=k2,kk
         if(SUM(II(k:kk,i2)).eq.1) Cycle; jsym2=i2; Exit
        End do
        Exit
       End do

       Call me_ee(io,jo,isym1,jsym1,isym2,jsym2)

       End Subroutine me_det


!======================================================================
      SUBROUTINE me_ee (i,j,i1,j1,i2,j2)
!======================================================================
!     angular part of matrix elements between two det.w.f.
!     for two-electron operator
!     Calls: Check_boef
!----------------------------------------------------------------------
      USE boef_list

      Implicit none
      Integer, intent(in) :: i,j,i1,i2,j1,j2
      Integer :: k,kz,ib

      if(IP_det(i1,kd1)+IP_det(j1,kd1).ne. &
         IP_det(i2,kd2)+IP_det(j2,kd2)) Return

      Call Check_boef(ln(i),jn(i),IP_det(i1,kd1), &
                      ln(j),jn(j),IP_det(j1,kd1), &
                      ln(i),jn(i),IP_det(i2,kd2), &
                      ln(j),jn(j),IP_det(j2,kd2))

      kz = (-1)**(i1+i2+j1+j2)

      Do ib = ncblk(kblk-1)+1,ncblk(kblk)
       k = IB_int(ib)
       if(k.gt.0) then
         k = k - 1
         coefs(j,i,k) = coefs(j,i,k) + Boef(ib)*CC_det*kz
       elseif(i.eq.j) then
         k = -k - 1
         coefs(j,i,k) = coefs(j,i,k) + Boef(ib)*CC_det*kz
       else
         k = -k - 1
         coefs(i,j,k) = coefs(i,j,k) + Boef(ib)*CC_det*kz
       end if
      End do

      End Subroutine me_ee

      End Subroutine coef_1conf





!=============================================================================
      Subroutine coef_2conf_jj(no1,nn1,ln1,jn1,iq1,Jshell1,Vshell1,Jintra1,  &
                               no2,nn2,ln2,jn2,iq2,Jshell2,Vshell2,Jintra2,  &
                               mcoef,ncoef,icoefs,coefs)
!=============================================================================
!     compute the angular coefficients for 2 different atomic states
!-----------------------------------------------------------------------------
      Implicit none

! ... input-output:

      Integer, intent(in) :: no1,nn1(no1),ln1(no1),jn1(no1),iq1(no1), &
                             Jshell1(no1),Vshell1(no1),Jintra1(no1),  &
                             no2,nn2(no2),ln2(no2),jn2(no2),iq2(no2), &
                             Jshell2(no2),Vshell2(no2),Jintra2(no2),  &
                             mcoef
      Integer ::  ncoef,icoefs(5,mcoef)
      Real(8) ::  coefs(mcoef)
      Real(8), allocatable :: coef(:,:,:,:,:)

! ... determinant expansion:

      Integer :: kdt, kdt1, kdt2
      Integer, allocatable :: IP_det(:,:), IP_det1(:,:), IP_det2(:,:)
      Real(8), allocatable :: C_det(:), C_det1(:), C_det2(:)
      Real(8) :: CC_det, eps_C = 1.d-7

! ... local variables:

      Integer              :: ne, Jtotal, i,j,k, kmax,kd1,kd2, i1,i2,j1,j2
      Integer, allocatable :: ip1(:),ip2(:)
      Integer, external    :: mj_value

! ... initialize arrays:

      ncoef = 0
      ne = SUM(iq1(1:no1));  if(ne.ne.SUM(iq2(1:no2))) Return
      Jtotal = Jintra1(no1); if(Jtotal.ne.Jintra2(no2)) Return

! ... determinant expansions:

      Call Det_expn_jj (no1,ln1,jn1,iq1,Jshell1,Vshell1,Jintra1)
      kdt1=kdt

      Allocate(C_det1(kdt1), IP_det1(ne,kdt1))
      C_det1(1:kdt1) = C_det(1:kdt1)
      IP_det1(1:ne,1:kdt1) = IP_det(1:ne,1:kdt1)


      Call Det_expn_jj (no2,ln2,jn2,iq2,Jshell2,Vshell2,Jintra2)
      kdt2=kdt

      Allocate(C_det2(kdt2), IP_det2(ne,kdt2))
      C_det2(1:kdt2) = C_det(1:kdt2)
      IP_det2(1:ne,1:kdt2) = IP_det(1:ne,1:kdt2)

      Deallocate(C_det,IP_det)

      if(allocated(ip1)) Deallocate(ip1); Allocate(ip1(ne))
      if(allocated(ip2)) Deallocate(ip2); Allocate(ip2(ne))

      k=1; Do i=1,no1; ip1(k:k+iq1(i)-1)=i; k=k+iq1(i); End do
      k=1; Do i=1,no2; ip2(k:k+iq2(i)-1)=i; k=k+iq2(i); End do


! ... calculations:

      kmax = (maxval(jn1(1:no1)) + maxval(jn2(1:no2)))/2

      if(allocated(coef)) Deallocate(coef)
      Allocate(coef(no1,no1,no2,no2,-1:kmax))
      coef = 0.d0

      Do kd1 = 1,kdt1
      Do kd2 = 1,kdt2
        CC_det = C_det1(kd1) * C_det2(kd2);  Call me_det
      End do;  End do

! ... final

      ncoef = 0
      Do i1=1,no1; Do i2=1,no1
      Do j1=1,no2; Do j2=1,no2
      Do k=0,kmax

       if(abs(coef(i1,i2,j1,j2,k)).lt.eps_c) Cycle

       ncoef=ncoef+1
       if(ncoef.gt.mcoef) Stop 'Coef_ee_2conf: ncoef > mcoef'
       coefs(ncoef)=coef(i1,i2,j1,j2,k)
       icoefs(1,ncoef)=k
       icoefs(2,ncoef)=i1
       icoefs(3,ncoef)=i2
       icoefs(4,ncoef)=j1
       icoefs(5,ncoef)=j2

      End do
      End do; End do
      End do; End do

! ... one-electron integrals:

      Do i=1,no1; Do j=1,no2
       if(abs(coef(i,i,j,j,-1)).lt.eps_c) Cycle
       ncoef=ncoef+1
       if(ncoef.gt.mcoef) Stop 'Coef_ee_2conf: ncoef > mcoef'
       coefs(ncoef)=coef(i,i,j,j,-1)
       icoefs(1,ncoef)=-1
       icoefs(2,ncoef)=i
       icoefs(3,ncoef)=i
       icoefs(4,ncoef)=j
       icoefs(5,ncoef)=j
      End do; End do

CONTAINS

!======================================================================
      Subroutine Det_expn_jj(no, ln,jn,iq,Jshell,Vshell,Jintra)
!======================================================================
!     procedure of exaustion of all possible determinants for given
!     configurations. The determinants and their coefficients
!     are recoded on unit 'nua'
!
!     Calls: Det_sh_jq, DETC_jq, Clebsh2, Ndets_jq, Jterm, mj_value
!----------------------------------------------------------------------
      Implicit none

      Integer, intent(in) :: no,ln(no),jn(no),iq(no), &
                             Jshell(no),Vshell(no),Jintra(no)

      Integer :: i,j,k, mkdt, JW,JQ
      Integer, external :: Ndets_jq, Jterm, mj_value
      Real(8) :: C
      Real(8), external :: DETC_jq, Clebsh2

! ... shell values:
!     md(i)  - the max.number of det.'s for the i-th subshell
!     nd(i)  - determinant under consideration
!     ipn(i) - pointers on the orbitals of given shell
!     MJs    - shells MJ
!     MJi    - intermediate values MJ

      Integer :: md(no),nd(no),ipn(no), MJs(no),MJi(no), ips(no)
      Integer :: Idet(ne)

! ... prepare working arrays:

      k = 1; mkdt = 1
      Do i=1,no
       ipn(i)=k; k=k+iq(i); md(i)=Ndets_jq(jn(i),iq(i))
       mkdt = mkdt*md(i)
       ips(i) = Jterm(jn(i),iq(i),0,Jshell(i),Vshell(i),JW,JQ)
      End do

      if(allocated(C_det)) Deallocate(C_det); Allocate(C_det(mkdt))
      if(allocated(IP_det)) Deallocate(IP_det); Allocate(IP_det(ne,mkdt))

!--------------------------------------------------------------------
! ... exausting all possible determinants:

      kdt=0; i=1; nd(i)=1
    1 Call DET_sh_jq(jn(i),iq(i),nd(i),MJs(i),Idet(ipn(i)))

      if(i.eq.1) then
       MJi(1) = MJs(1)
      else
       MJi(i) = MJi(i-1)+MJs(i)
      end if
      if(i.lt.no) then;  i=i+1;  nd(i)=1;  go to 1; end if

       ! ... coefficient calculation:

       if(MJi(no).ne.Jtotal) go to 2
       C = DETC_jq(jn(1),iq(1),ips(1),nd(1))
       if(C.eq.0.d0) go to 2

       Do j=2,no
        C = C * DETC_jq(jn(j),iq(j),ips(j),nd(j))
        if(C.eq.0.d0) Exit
        C = C * Clebsh2(Jintra(j-1),MJi(j-1), &
                        Jshell(j  ),MJs(j  ), &
                        Jintra(j  ),MJi(j  ))
        if(C.eq.0.d0) Exit
       End do

       if(C.ne.0.d0) then
        kdt=kdt+1; IP_det(1:ne,kdt)=Idet(1:ne); C_det(kdt)=C
       end if

    2 nd(i)=nd(i)+1                ! selecting the next case

      if(nd(i).gt.md(i)) then
       if(i.eq.1) go to 3          ! to end
       i=i-1; go to 2
      end if
      go to 1

    3 Continue

      Do k=1,kdt; Do i=1,ne
       j=IP_det(i,k); IP_det(i,k)=mj_value(j)
      End do; End do

      End Subroutine Det_expn_jj


!======================================================================
      Subroutine me_det
!======================================================================
!     find the possible interaction orbitals in two determinants and
!     call the subroutines to calculate m.e. between possible
!     combinations of nj-orbitals (me_ee)
!----------------------------------------------------------------------
      Implicit none

      Integer :: i,i1,i2, j,j1,j2, idif, jdif
      Integer :: ii(ne),jj(ne)

!----------------------------------------------------------------------
! ... check total orbitals difference:

       ii = 0; jj=0
       Do i1=1,ne; j1=ip1(i1)
       Do i2=1,ne; j2=ip2(i2)
        if(nn1(j1).ne.nn2(j2)) Cycle
        if(ln1(j1).ne.ln2(j2)) Cycle
        if(jn1(j1).ne.jn2(j2)) Cycle
        if(IP_det1(i1,kd1).ne.IP_det2(i2,kd2)) Cycle
        ii(i1)=i2; jj(i2)=i1
       End do; End do

       idif = 0;  Do i=1,ne; if(ii(i).eq.0) idif=idif+1; End do
       jdif = 0;  Do i=1,ne; if(jj(i).eq.0) jdif=jdif+1; End do

       if(idif.ne.jdif) Stop 'me_det: idif <> jdif'
       if(idif.gt.2) Return

!----------------------------------------------------------------------
       Select case(idif)

       Case(0)

        Do i=1,ne; i1=ip1(i); i2=ip2(i)
         coef(i1,i1,i2,i2,-1) = coef(i1,i1,i2,i2,-1) + CC_det
        End do

        Do i1=1,ne-1;  j1=ii(i1)
        Do i2=i1+1,ne; j2=ii(i2)
         Call me_ee(i1,i2,min(j1,j2),max(j1,j2))
        End do; End do

       Case(1)

        Do i = 1,ne; if(ii(i).ne.0) Cycle; i1=i; Exit; End do
        Do j = 1,ne; if(jj(j).ne.0) Cycle; j1=j; Exit; End do

        i=ip1(i1); j=ip2(j1)
        if(IP_det1(i1,kd1).eq.IP_det2(j1,kd2).and. &
           jn1(i).eq.jn2(j).and.ln1(i).eq.ln2(j)) then
         coef(i,i,j,j,-1) = coef(i,i,j,j,-1) + CC_det * (-1)**(i1+j1)
        end if

        Do i2=1,ne; if(i2.eq.i1) Cycle; j2=ii(i2)
         Call me_ee(min(i1,i2),max(i1,i2),min(j1,j2),max(j1,j2))
        End do

       Case(2)

        Do i = 1,ne; if(ii(i).ne.0) Cycle; i1=i; Exit; End do
        Do i = ne,1,-1; if(ii(i).ne.0) Cycle; i2=i; Exit; End do

        Do j = 1,ne; if(jj(j).ne.0) Cycle; j1=j; Exit; End do
        Do j = ne,1,-1; if(jj(j).ne.0) Cycle; j2=j; Exit; End do

        Call me_ee(i1,i2,j1,j2)

       End Select

       End Subroutine me_det


!======================================================================
      SUBROUTINE me_ee (i1,j1,i2,j2)
!======================================================================
!     angular part of matrix elements in nljm-representation
!     for two-electron operator
!     Calls: Check_boef
!----------------------------------------------------------------------
      Use boef_list

      Implicit none
      Integer, intent(in) :: i1,i2,j1,j2
      Integer :: k,kz,ib, n1,n2,n3,n4

      if(IP_det1(i1,kd1)+IP_det1(j1,kd1).ne. &
         IP_det2(i2,kd2)+IP_det2(j2,kd2)) Return

      n1=ip1(i1); n2=ip1(j1); n3=ip2(i2); n4=ip2(j2)

      Call Check_boef(ln1(n1),jn1(n1),IP_det1(i1,kd1), &
                      ln1(n2),jn1(n2),IP_det1(j1,kd1), &
                      ln2(n3),jn2(n3),IP_det2(i2,kd2), &
                      ln2(n4),jn2(n4),IP_det2(j2,kd2))

      kz = (-1)**(i1+i2+j1+j2)

      Do ib = ncblk(kblk-1)+1,ncblk(kblk)
       k = IB_int(ib)
       if(k.gt.0) then
        k = k - 1;    if(k.gt.kmax) Cycle
        coef(n1,n2,n3,n4,k) = coef(n1,n2,n3,n4,k) + Boef(ib)*CC_det*kz
       else
        k = -k - 1;  if(k.gt.kmax) Cycle
        coef(n1,n2,n4,n3,k) = coef(n1,n2,n4,n3,k) + Boef(ib)*CC_det*kz
       end if
      End do

      End Subroutine me_ee

      End Subroutine coef_2conf_jj


!=======================================================================
      Subroutine EL_nljk(EL,n,kappa,l,j,k)
!=======================================================================
!     decodes the spectroscopic notation for electron orbital (n,l,j,k)
!
!     It is allowed following notations: 1s , 2s 3, 2p-30, 20p-3,
!     1s h, 20s h, kp , kp 1, kp-11, ns , ns 3, ns 33, ... .
!
!     Call:  LA, INDEX
!----------------------------------------------------------------------
      Use orb_jj, only: kset, ASET

      Implicit none
      Character(5), intent(in) :: EL
      Integer, intent(out) :: n,l,j,k,kappa
      Integer :: jj, k1,k2, n1,n2
      Integer, external :: LA, kappa_lj

      jj=0
      Do j=5,3,-1
       if(EL(j:j).eq.'-'.or.EL(j:j).eq.'+') then; jj=j; Exit; end if
      End do
      if(jj.eq.0) then
       Do j=5,3,-1
        if(EL(j:j).eq.' '.and.EL(j-1:j-1).ne.' ') then
         jj=j; Exit
        end if
       End do
      end if

      if(jj.eq.5) then

       k = 0
       l = LA(EL(4:4))
       n1 = INDEX(ASET,EL(3:3))
       n2 = INDEX(ASET,EL(2:2))
       n = n2*kset+n1

      elseif(jj.eq.4) then

       k = INDEX(ASET,EL(5:5))
       l = LA(EL(3:3))
       n1 = INDEX(ASET,EL(2:2))
       n2 = INDEX(ASET,EL(1:1))
       n = n2*kset+n1

      elseif(jj.eq.3) then

       k1 = INDEX(ASET,EL(5:5))
       k2 = INDEX(ASET,EL(4:4))
       k = k2*kset+k1
       l = LA(EL(2:2))
       n = INDEX(ASET,EL(1:1))

      else

       write(*,*) 'EL_NLJK: can not decode ',EL
       Stop ' '

      end if

      j = l+l+1; if(EL(jj:jj).eq.'-') j = l+l-1
      kappa = kappa_lj(l,j)

      End Subroutine EL_NLJK



!=======================================================================
      Character(5) Function ELi(n,kappa,k)
!=======================================================================
!     provides the spectroscopic notation for electron orbital (n,l,j,k)
!     set index k must be < 61*61 (see ASET for incoding k)
!-----------------------------------------------------------------------
      Use orb_jj, only: kset, ASET

      Implicit none
      Integer :: n,l,j,k, ll,jj, i,k1,k2,n1,n2, kappa
      Character(5) :: EL
      Character(1), external :: AL
      Integer, external :: l_kappa, j_kappa

      l = l_kappa(kappa)
      j = j_kappa(kappa)

      if(n.le.0.or.l.lt.0.or.j.le.0.or.k.lt.0) then
       write(*,'(a,4i5)') &
           ' ELj: parmeters are out of limits: n,l,j,k=',n,l,j,k
       Stop 'Stop in Elj'
      end if

      EL='    '; i=5

      if(k.lt.0) Stop 'ELi: set index < 0'
      if(k.gt.0) then
       if(k.le.kset) then
        EL(i:i)=ASET(k:k); i=i-1
       else
        k1=k/kset; k2=mod(k,kset);
        if(k2.eq.0) then; k1=k1-1; k2=kset; end if
        if(k1.gt.kset) Stop 'ELi: set index too big'
        EL(i:i)=ASET(k2:k2); i=i-1
        EL(i:i)=ASET(k1:k1); i=i-1
       end if
      end if

      if(j.eq.l+l-1) EL(i:i) = '-'; i=i-1

      EL(i:i)=AL(l,1);  i=i-1

      if(n.lt.0) Stop 'ELi: n < 0'
      if(n.gt.0) then
       if(n.le.kset) then
        EL(i:i)=ASET(n:n)
       else
        n1=n/kset; n2=mod(n,kset);
        if(n2.eq.0) then; n1=n1-1; n2=kset; end if
        if(n1.gt.kset) Stop 'ELi: n is too big'
        EL(i:i)=ASET(n2:n2); i=i-1
        EL(i:i)=ASET(n1:n1); i=i-1
       end if
      end if

      ELi = EL

      End Function ELi


!=======================================================================
      Character(5) Function ELj(n,l,j,k)
!=======================================================================
!     provides the spectroscopic notation for electron orbital (n,l,j,k)
!     set index k must be < 61*61 (see ASET for incoding k)
!-----------------------------------------------------------------------
      Use orb_jj, only: kset, ASET

      Implicit none
      Integer :: n,l,j,k, ll,jj, i,k1,k2,n1,n2, kappa
      Character(5) :: EL
      Character(1), external :: AL

      if(n.le.0.or.l.lt.0.or.j.le.0.or.k.lt.0) then
       write(*,'(a,4i5)') &
           ' ELj: parmeters are out of limits: n,l,j,k=',n,l,j,k
       Stop 'Stop in Elj'
      end if

      EL='    '; i=5

      if(k.lt.0) Stop 'ELi: set index < 0'
      if(k.gt.0) then
       if(k.le.kset) then
        EL(i:i)=ASET(k:k); i=i-1
       else
        k1=k/kset; k2=mod(k,kset);
        if(k2.eq.0) then; k1=k1-1; k2=kset; end if
        if(k1.gt.kset) Stop 'ELi: set index too big'
        EL(i:i)=ASET(k2:k2); i=i-1
        EL(i:i)=ASET(k1:k1); i=i-1
       end if
      end if

      if(j.eq.l+l-1) EL(i:i) = '-'; i=i-1

      EL(i:i)=AL(l,1);  i=i-1

      if(n.lt.0) Stop 'ELi: n < 0'
      if(n.gt.0) then
       if(n.le.kset) then
        EL(i:i)=ASET(n:n)
       else
        n1=n/kset; n2=mod(n,kset);
        if(n2.eq.0) then; n1=n1-1; n2=kset; end if
        if(n1.gt.kset) Stop 'ELi: n is too big'
        EL(i:i)=ASET(n2:n2); i=i-1
        EL(i:i)=ASET(n1:n1); i=i-1
       end if
      end if

      ELj = EL

      End Function ELj



!====================================================================
      Character FUNCTION AL(L,k)
!====================================================================
!     provides spectroscopic symbols for L values ( L <= 153 )
!     k=1 - small letters; k=2 - capital letters
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: L,K
      Integer :: I
      Character(21) :: AS, AB

      DATA AS/'spdfghiklmnoqrtuvwxyz'/
      DATA AB/'SPDFGHIKLMNOQRTUVWXYZ'/

      i = L+1; IF(k.eq.5.or.k.eq.6) i=(L-1)/2+1
      if(i.ge.1.and.i.le.21) then
       if(k.eq.1.or.k.eq.5) AL=AS(I:I)
       if(k.eq.2.or.k.eq.6) AL=AB(I:I)
      elseif(i.ge.22.and.i.le.153) then
       AL=CHAR(i+101)  ! from '{' and futher
      else
       write(*,*) 'L,k=',L,k
       Stop ' AL: L is out of range'
      end if

      End FUNCTION AL


!====================================================================
      Integer FUNCTION LA(a)
!====================================================================
!     gives the value of L from spetroscopic symbol "a"
!--------------------------------------------------------------------
      Implicit none
      Character, Intent(in) :: a
      Character(21) :: AS, AB
      Integer :: i

      DATA AS/'spdfghiklmnoqrtuvwxyz'/
      DATA AB/'SPDFGHIKLMNOQRTUVWXYZ'/

      i = INDEX (AS,a)
      if(i.eq.0) i = INDEX (AB,a)
      if(i.eq.0) i = ICHAR(a)-101
      LA = i-1

      End FUNCTION LA



!======================================================================
      Subroutine Incode_cj
!======================================================================
!     encodes the configuration from INTEGER format to the GRASP
!     c-file format (see module conf_jj for variables)
!----------------------------------------------------------------------
      Use conf_jj
      Implicit none
      Integer :: i,j,m,k,j1,j2,jj
      Character(5), external :: ELi

      CONFIG=' ';  SHELLJ=' ';   INTRAJ=' '

      m=-9; jj=0
      Do i=1,no
       m = m + 9
       CONFIG(m+1:m+5) = ELi(nn(i),kn(i),in(i))
       write(CONFIG(m+6:m+9),'(a1,i2,a1)') '(',iq(i),')'
       SHELLJ(m+1:m+9) = '        '
       INTRAJ(m+1:m+9) = '        '
       if(iq(i).eq.jn(i)+1) Cycle

       k = mod(Jshell(i),2)
       if(k.eq.0) then
        write(SHELLJ(m+1:m+9),'(i9)') Jshell(i)/2
       else
        write(SHELLJ(m+1:m+9),'(i7,a2)') Jshell(i),'/2'
       end if
       if(jn(i).eq.7.and.iq(i).eq.4.and.Jshell(i).gt.0)  &
         write(SHELLJ(m+1:m+5),'(i4,a1)') Vshell(i),';'

       if(i.eq.1.or.i.eq.no) Cycle

       if(SHELLJ(m+1:m+9).ne.'        ') then
        Do j=1,i-1; j1=(j-1)*9+1; j2=j*9
         if(SHELLJ(j1:j2).ne.'        ') jj=1
        End do
       end if

       !if(Jintra(i).eq.0) Cycle
       if(Jshell(i).eq.0) Cycle

       j1 = iabs(Jintra(i-1)-Jshell(i))
       j2 = iabs(Jintra(i-1)+Jshell(i))
       if(j1.eq.j2.and.jj.eq.0) Cycle
       k = mod(Jintra(i),2)
       if(k.eq.0) then
        write(INTRAJ(m+1:m+9),'(i9)') Jintra(i)/2
       else
        write(INTRAJ(m+1:m+9),'(i7,a2)') Jintra(i),'/2'
       end if
       jj = 1
       Call Clean_a(INTRAJ(m+1:m+9))
      End do

      k = mod(Jintra(no),2)
      if(k.eq.0) then
       write(INTRAJ(m+1:m+7),'(i7)') Jintra(no)/2
      else
       write(INTRAJ(m+1:m+7),'(i5,a2)') Jintra(no),'/2'
      end if

      k=0; Do i=1,no; k=k+ln(i)*iq(i); End do; parity=(-1)**k
      if(parity.eq.+1) INTRAJ(m+8:m+9)='+ '
      if(parity.eq.-1) INTRAJ(m+8:m+9)='- '
      Call Clean_a(INTRAJ(m+1:m+9))

      ia = no*9

      End Subroutine Incode_cj


!======================================================================
      Subroutine Decode_cj
!======================================================================
!     decode the configuration from c-file format to INTEGER format
!----------------------------------------------------------------------
      Use conf_jj
      Implicit none
      Integer :: i,j,k,m
      Character(9) :: bl = '        '

      ia=INDEX(CONFIG,')',BACK=.TRUE.); no=ia/9

      Vshell=0; Jshell=0; Jintra=0
      m=-9
      Do i=1,no
       m = m + 9

       Call EL_nljk(CONFIG(m+1:m+5),nn(i),kn(i),ln(i),jn(i),in(i))
       read(CONFIG(m+7:m+8),'(i2)') iq(i)

       if(SHELLJ(m+5:m+5).eq.';') then
        read(SHELLJ(m+4:m+4),'(i1)') Vshell(i)
        read(SHELLJ(m+9:m+9),'(i1)') J; Jshell(i) = 2*J
       elseif(SHELLJ(m+8:m+8).eq.'/') then
        read(SHELLJ(m+1:m+7),'(i7)') Jshell(i)
       elseif(SHELLJ(m+9:m+9).ne.' ') then
        read(SHELLJ(m+1:m+9),'(i9)') J; Jshell(i) = 2*J
       end if

       k = INDEX(INTRAJ(m+1:m+9),'/')
       if(i.eq.1) then
        Jintra(i) = Jshell(i)
       elseif(i.eq.no) then
        Cycle
       elseif(k.gt.0) then
        read(INTRAJ(m+1:m+k-1),*) Jintra(i)
       elseif(INTRAJ(m+1:m+9).ne.bl) then
        read(INTRAJ(m+1:m+9),*) J; Jintra(i) = 2*J
       else
       if(Jshell(i).eq.0) Jintra(i) = Jintra(i-1)
       if(Jintra(i-1).eq.0) Jintra(i) = Jshell(i)
       end if

      End do

      if(k.gt.0) then
       read(INTRAJ(m+1:m+k-1),*) Jintra(no)
      else
       k=INDEX(INTRAJ(m+1:m+9),'+')
       if(k.eq.0) k=INDEX(INTRAJ(m+1:m+9),'-')
       read(INTRAJ(m+1:m+k-1),*) J; Jintra(no) = 2*J
      end if

      Do i=1,no
       m = jn(i)+1
       if(Vshell(i).ne.0.or.iq(i).eq.m) Cycle
       if(iq(i).eq.1.or.iq(i).eq.m-1) then
         Vshell(i) = 1
       elseif(iq(i).eq.2.or.iq(i).eq.m-2) then
         if(Jshell(i).gt.0) Vshell(i) = 2
       elseif(iq(i).eq.3) then
         Vshell(i) = 3; if(Jshell(i).eq.jn(i)) Vshell(i)=1
       end if
      End do

      Jtotal = Jintra(no)

      End Subroutine Decode_cj


!======================================================================
      Subroutine Decode_confj
!======================================================================
!     decodes configuration from c-file format to INTEGER format
!----------------------------------------------------------------------
      Use conf_jj
      Implicit none
      Integer :: i,k,m

      no = INDEX(CONFIG,')',BACK=.TRUE.) / 9
      m=0
      Do i=1,no
       Call EL_nljk(CONFIG(m+1:m+5),nn(i),kn(i),ln(i),jn(i),in(i))
       read(CONFIG(m+7:m+8),'(i2)') iq(i)
       k = (2*ln(i)-jn(i))*(jn(i)+1)/2
       np_symc(i) = k*1000 + iq(i)
       m = m + 9
      End do

      End Subroutine Decode_confj

!======================================================================
      Subroutine Incode_confj
!======================================================================
!     incodes the configuration from INTEGER format to c-file format
!----------------------------------------------------------------------
      Use conf_jj
      Implicit none
      Integer :: i,m
      Character(5), External :: ELi

      m=0
      Do i=1,no
       if(iq(i).le.0) Cycle
       CONFIG(m+1:m+5) = ELi(nn(i),kn(i),in(i))
       write(CONFIG(m+6:m+9),'(a1,i2,a1)') '(',iq(i),')'
       m = m + 9
      End do

      ia = no*9

      End Subroutine Incode_confj


!======================================================================
      Subroutine Incode_confj1
!======================================================================
!     incodes the configuration (with index 1 in module conf_jj)
!     from INTEGER format to c-file format
!----------------------------------------------------------------------
      Use conf_jj
      Implicit none
      Integer :: i,m
      Character(5), External :: ELi

      m=0
      Do i=1,no1
       if(iq(i).le.0) Cycle
       CONFIG(m+1:m+5) = ELi(nn1(i),kn1(i),in1(i))
       write(CONFIG(m+6:m+9),'(a1,i2,a1)') '(',iq1(i),')'
       m = m + 9
      End do

      ia = no1*9

      End Subroutine Incode_confj1


!======================================================================
      Subroutine Decode_cjj(CONFIG,SHELLJ,INTRAJ,no,nn,kn,ln,jn,iq,in,&
                            Jshell,Vshell,Jintra)
!======================================================================
!     decode the configuration from c-file format to INTEGER format
!     Call: EL_nljk
!----------------------------------------------------------------------

      Implicit none
      Character(*), Intent(in) :: CONFIG,SHELLJ,INTRAJ
      Integer, Intent(out) :: no,nn(*),kn(*),ln(*),jn(*),iq(*),in(*),&
                              Jshell(*),Vshell(*),Jintra(*)
      Integer :: i,j,k,m

      m=INDEX(CONFIG,')',BACK=.TRUE.); no=m/9

      Vshell(1:no)=0; Jshell(1:no)=0; Jintra(1:no)=0
      m=-9
      Do i=1,no
       m = m + 9

       Call EL_nljk(CONFIG(m+1:m+5),nn(i),kn(i),ln(i),jn(i),in(i))
       read(CONFIG(m+7:m+8),'(i2)') iq(i)

       if(SHELLJ(m+5:m+5).eq.';') then
        read(SHELLJ(m+4:m+4),'(i1)') Vshell(i)
        read(SHELLJ(m+9:m+9),'(i1)') J; Jshell(i) = 2*J
       elseif(SHELLJ(m+8:m+8).eq.'/') then
        read(SHELLJ(m+1:m+7),'(i7)') Jshell(i)
       elseif(SHELLJ(m+9:m+9).ne.' ') then
        read(SHELLJ(m+1:m+9),'(i9)') J; Jshell(i) = 2*J
       end if

       k = INDEX(INTRAJ(m+1:m+9),'/')
       if(i.eq.1) then
        Jintra(i) = Jshell(i)
       elseif(i.eq.no) then
        Cycle
       elseif(k.gt.0) then
        read(INTRAJ(m+1:m+k-1),*) Jintra(i)
       elseif(LEN_TRIM(INTRAJ(m+1:m+9)).ne.0) then
        read(INTRAJ(m+1:m+9),*) J; Jintra(i) = 2*J
       else
       if(Jshell(i).eq.0) Jintra(i) = Jintra(i-1)
       if(Jintra(i-1).eq.0) Jintra(i) = Jshell(i)
       end if

      End do

      if(k.gt.0) then
       read(INTRAJ(m+1:m+k-1),*) Jintra(no)
      else
       k=INDEX(INTRAJ(m+1:m+9),'+')
       if(k.eq.0) k=INDEX(INTRAJ(m+1:m+9),'-')
       read(INTRAJ(m+1:m+k-1),*) J; Jintra(no) = 2*J
      end if

      Do i=1,no
       m = jn(i)+1
       if(Vshell(i).ne.0.or.iq(i).eq.m) Cycle
       if(iq(i).eq.1.or.iq(i).eq.m-1) then
         Vshell(i) = 1
       elseif(iq(i).eq.2.or.iq(i).eq.m-2) then
         if(Jshell(i).gt.0) Vshell(i) = 2
       elseif(iq(i).eq.3) then
         Vshell(i) = 3; if(Jshell(i).eq.jn(i)) Vshell(i)=1
       end if
      End do

      End Subroutine Decode_cjj


!======================================================================
      Integer Function Jdef_ncfg(nu)
!======================================================================
!     defines the number of configuration in c-file (unit nu)
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu
      Integer :: ncfg
      Character(6) :: AS

      rewind(nu)
      ncfg=0
    1 read(nu,'(a)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      if(AS(6:6).ne.'(') go to 1
      ncfg=ncfg+1
      go to 1
    2 rewind(nu)
      Jdef_ncfg=ncfg

      End Function Jdef_ncfg


!======================================================================
      Subroutine Jdef_kcfg(nu,ncfg,kcfg)
!======================================================================
!     defines the number of configuration in c-file (unit nu),
!     plus the total "length" of these configurations, kcfg
!----------------------------------------------------------------------
      Implicit none
      Integer :: i,nu,ncfg,kcfg
      Character(300) :: AS

      rewind(nu)
      ncfg=0; kcfg=0
    1 read(nu,'(a)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      if(AS(6:6).ne.'(') go to 1
      ncfg=ncfg+1
      i = INDEX(AS,')',BACK=.TRUE.)
      if((i/9)*9.ne.i) Stop 'Jdef_ncfg: problems with kcfg'
      kcfg = kcfg + i/9
      go to 1
    2 rewind(nu)

      End Subroutine Jdef_kcfg




!======================================================================
      Integer Function Jterm (j,q,k,JT,JV,JW,JQ)
!======================================================================
!     provides information about possible terms in the j^q subshell
!     for  j <= 9/2
!
!     Options:
!     k > 0  --> JT,JV,JQ,JM = k-th term of j^q-subshell
!     k = 0  --> Jterm = position of the (JT,JV) term in subshell list
!     k < 0  --> Jterm = number of terms in j^q-subshell
!
!     j, JT, JQ  -->  in 2*J representation
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: j,q,k
      Integer :: JT,JV,JQ,JW
      Integer :: qm,nterm,ip,i,ii,qq
      Integer :: term_list(3,65)

      Data term_list                                               &
       /1, 5, 2,  3, 3, 0,  3, 9, 0,                               & ! 5/2 ^ 3
        1, 7, 3,  3, 3, 1,  3, 5, 1,  3, 9, 1,  3,11, 1,  3,15, 1, & ! 7/2 ^ 3
        0, 0, 4,  2, 4, 2,  2, 8, 2,  2,12, 2,  4, 4, 0,  4, 8, 0, & ! 7/2 ^ 4
        4,10, 0,  4,16, 0,                                         &
        1, 9, 4,  3, 3, 2,  3, 5, 2,  3, 7, 2,  3, 9, 2,  3,11, 2, & ! 9/2 ^ 3
        3,13, 2,  3,15, 2,  3,17, 2,  3,21, 2,                     &
        0, 0, 5,  2, 4, 3,  2, 8, 3,  2,12, 3,  2,16, 3,  4, 0, 1, & ! 9/2 ^ 4
        4, 4, 1,  4, 6, 1,  4, 8,10,  4, 8,20,  2,10, 1,  4,12,10, &
        4,12,20,  4,14, 1,  4,16, 1,  4,18, 1,  4,20, 1,  4,24, 1, &
        1, 9, 4,  3, 3, 2,  3, 5, 2,  3, 7, 2,  3, 9, 2,  3,11, 2, & ! 9/2 ^ 5
        3,13, 2,  3,15, 2,  3,17, 2,  3,21, 2,  5, 1, 0,  5, 5, 0, &
        5, 7, 0,  5, 9, 0,  5,11, 0,  5,13, 0,  5,15, 0,  5,17, 0, &
        5,19, 0,  5,25, 0/

!       JW = JQ/10  -> additional seniority
!----------------------------------------------------------------------
! ... check input data:

      if(j.lt.0.or.mod(j,2).eq.0) then
       write(*,*) 'Jterm: unphysical j-value: j =',j
       Stop  'Stop in Jterm'
      end if

      qm = j + 1
      if(q.lt.0.or.q.gt.qm) then
       write(*,*) 'Jterm: number of electron is out of range: q =',q
       Stop  'Stop in Jterm'
      end if

      if(j.gt.9.and.q.gt.2) then
       write(*,*) 'Jterm: j^q out of scope: j,q =',j,q
       Stop  'Stop in Jterm'
      end if

      Jterm = 0

!----------------------------------------------------------------------
! ... the number of terms in j^q subshell:

      if(q.le.1.or.q.ge.qm-1) then
       nterm = 1
      elseif(q.eq.2.or.q.eq.qm-2) then
       nterm = qm/2
      else
       Select case(j*100 + q)
        Case(503);       nterm = 3; ip = 0    ! 5/2 ^ 3
        Case(703,705);   nterm = 6; ip = 3    ! 7/2 ^ 3,5
        Case(704);       nterm = 8; ip = 9    ! 7/2 ^ 4
        Case(903,907);   nterm =10; ip =17    ! 9/2 ^ 3
        Case(904,906);   nterm =18; ip =27    ! 9/2 ^ 4,6
        Case(905);       nterm =20; ip =45    ! 9/2 ^ 5
        Case default;
         write(*,*) 'Jterm: cannot find j^q subshell for j,q=',j,q
         Stop 'Stop in Jterm'
       End Select
      end if

      if(k.lt.0) then; Jterm = nterm; Return; end if
!----------------------------------------------------------------------
! ... position of the JT,JV term in the subshell list

      if(k.eq.0) then

      qq = min0(q,qm-q)
      Select case(qq)
       case(0);  JV=0
       case(1);  JV=1
       case(2);  JV=2
       case(3);  JV=3; if(j.eq.JT) JV=1
       case(4);  if(j.eq.7) then
                  if(JT.eq.12) JV=2
                  if(JT.eq.10) JV=4
                  if(JT.eq.16) JV=4
                 end if
                 if(j.eq.9.and.JT.ge.18) JV=4
       case(5);  if(j.eq.9) then
                  if(JT.eq. 1) JV=5
                  if(JT.eq. 3) JV=3
                  if(JT.eq.19) JV=5
                  if(JT.eq.21) JV=3
                  if(JT.eq.25) JV=5
                 end if
      End select
      if(JT.eq.0.and.j.lt.9) JV=0

       if(q.le.1.or.q.ge.j) then
        Jterm =  1
       elseif(q.eq.2.or.q.eq.qm-2) then
        Jterm =  JT/4+1
       else
        Do i=1,nterm; ii = ip + i
         if(JV.ne.term_list(1,ii).or.JT.ne.term_list(2,ii)) Cycle
         Jterm=i; Exit
        End do
      end if

      if(Jterm.eq.0) then
       write(*,*) 'Jterm:  incorect term j^q(JV,JT):',j,q,JV,JT
       write(*,*) '2j =',j
       write(*,*) 'q  =',q
       write(*,*) '2J =',JT
       write(*,*) 'v  =',JV
       Stop 'Stop in Jterm'
      end if

      Return
      end if

!----------------------------------------------------------------------
! ... k-th term of j^q subshell

      if(k.gt.nterm) then
       write(*,*) 'Jterm:  incorect input k =',k
       Stop 'Stop in Jterm'
      end if

      JW = 0
      if(q.eq.0.or.q.eq.qm) then                 !  j ^ 0
       JV = 0; JT = 0; JQ = qm/2
      elseif(q.eq.1.or.q.eq.qm-1) then           !  j ^ 1
       JV = 1; JT = j; JQ = qm/2-1
      elseif(q.eq.2.or.q.eq.qm-2) then           !  j ^ 2
       JV = 0; if(k.gt.1) JV = 2
       JT = (k-1) * 4
       JQ = qm/2; if(k.gt.1) JQ = JQ - 2
      else                                       !  j ^ q
       ii = ip + k
       JV = term_list(1,ii)
       JT = term_list(2,ii)
       JQ = term_list(3,ii)
       if(JQ.ge.10) then; JW = JQ/10; JQ = 1; end if
      end if

      End Function Jterm


!====================================================================
      Integer Function Index_mj (mj)
!====================================================================
!     index of jm-orbital in the common list: -1,+1,-3,+3,-5,+5, ...
!     mj -> in 2j-representation
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: mj
      Index_mj = iabs(mj)
      if(mj.gt.0) index_mj = index_mj + 1
      END Function Index_mj

!====================================================================
      Integer Function mj_value(i)
!====================================================================
!     mj value for orbital 'i' in the lit: -1,+1,-3,+3,-5,+5, ...
!     mj -> in 2j-representation
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: i
      mj_value = -i
      if(mod(i,2).eq.0) mj_value = i - 1
      END Function mj_value

!====================================================================
      Integer Function ndets_jq(j,q)
!====================================================================
!     number of det.s in subshells j^k (Newton's binom)
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: j,q
      Integer :: i
      Real(8) :: S
      if(q.gt.j+1) Stop 'ndets_jq:  q > q_max'
      S=1.d0
      Do i=q+1,j+1; S=S*i/(i-q); End do
      ndets_jq = S + 0.1d0
      End Function ndets_jq



!=====================================================================
      Subroutine Gen_jj_states (AF_inp,AF_out,jmin,jmax)
!======================================================================
!     preparation  the ASF list from a list of (nlj)^q configurations
!--------------------------------------------------------------------
      Use conf_jj; Use orb_jj

      Character(*), intent(in) :: AF_inp, AF_out
      Integer :: jmin, jmax, nu1,nu2, i,ii,ic
      Integer, external :: Icheck_file

! ... Check the files:

      if(Icheck_file(AF_inp).eq.0) Return
      Call Find_free_unit(nu1)
      open(nu1,file=AF_inp); rewind(nu1)
      Call Find_free_unit(nu2)
      open(nu2,file=AF_out); rewind(nu2)

! ... reload the header:

      Do
       read(nu1,'(a)',end=10) AS; write(nu2,'(a)') trim(AS)
       if(AS.eq.'CSF(s):') Exit
      End do

! ... check jmin and jmax:

      if(jmin.eq.-1.or.jmax.eq.-1) then
       i=0
       Do
        read(nu1,'(a)',end=10) AS
        if(AS(6:6).ne.'(') Cycle
        CONFIG = AS;  Call Decode_confj
        ii = SUM(iq(1:no)*jn(1:no)); if(ii.gt.i) i=ii
       End do
      end if
   10 Continue
      if(jmin.eq.-1) jmin=mod(i,2)
      if(jmax.eq.-1) jmax = i
      if(jmax.lt.jmin) jmax=jmin

! ... read configurations for each J-total:

      Do j = jmin,jmax,2;   J_min=j; J_max=j

      ncfg=0
      rewind(nu1)
    1 read(nu1,'(a)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      if(AS(6:6).ne.'(') go to 1

      CONFIG = AS;  Call Decode_confj

      Call Sum_Term

      go to 1
    2 Rewind(nu1)

      if(ncfg.eq.0) Cycle

      Do ic=1,ncfg
       Call Get_cfg_jj(ic)
       Call Incode_cj
       ii = no*9
       write(nu2,'(a)') CONFIG(1:ii)
       write(nu2,'(a)') SHELLJ(1:ii)
       write(nu2,'(9x,a)') INTRAJ(1:ii)
      End do
      if(j.lt.jmax) write(nu2,'(a)') ' *'
      if(j.eq.jmax) write(nu2,'(a)') '* '

      End do ! over j

      BACKSPACE(nu2);  write(nu2,'(a)') '* '


      End  Subroutine Gen_jj_states


!----------------------------------------------------------------------
      Subroutine Sum_Term
!----------------------------------------------------------------------
!     exhaustion of shell-terms
!----------------------------------------------------------------------
      Use conf_jj; Use orb_jj

      Implicit none
      Integer :: mt(msh),nt(msh)
      Integer :: i,ii,i1,i2, JT,JV,JW,JQ
      Integer, external :: Jterm

!     mt(i) - the number of term in shell i
!     nt(i) - the term inder consideration

      i1=1                     ! i1 - low  limit of shells
      i2=no                    ! i2 - high limit of shells
      Do i=i1,i2
       mt(i)=Jterm(jn(i),iq(i),-1,JT,JV,JW,JQ)
      End do

      i=i1                     ! first shell under consideration
      nt = 1

    1 Continue

      ii = Jterm(jn(i),iq(i),nt(i),Jshell(i),Vshell(i),JW,JQ)

      if(i.lt.i2) then
         i=i+1; nt(i)=1; go to 1
      else
         CALL Sum_Iterm
      end if

    2 nt(i)=nt(i)+1
      if(nt(i).gt.mt(i)) then
        if(i.eq.i1) go to 3
        i=i-1; go to 2
        end if
      go to 1

    3 Continue

      End Subroutine Sum_Term


!----------------------------------------------------------------------
      Subroutine Sum_Iterm
!----------------------------------------------------------------------
!     exhaustion of intermediate terms
!----------------------------------------------------------------------
      USE conf_jj; Use orb_jj

      Integer :: js_min(msh),js_max(msh)

      Jintra(1)=Jshell(1)
      if(no.eq.1) then
       if(Jshell(no).ge.J_min.and.Jshell(no).le.J_max)  ic=Iadd_cfg_jj('detect')
       Return
      end if

      i1=2                         ! i1 - low  limit
      i2=no                        ! i2 - high limit in array LS(...)

      i=i1
    1 j1=i-1; j2=i

      js_min(i)=IABS(Jintra(j1)-Jshell(j2))
      js_max(i)=     Jintra(j1)+Jshell(j2)
      Jintra(i)=js_min(i)

    2 if(i.lt.i2) then
       i=i+1; go to 1
      else

       if(Jintra(no).ge.J_min.and.Jintra(no).le.J_max)  ic=Iadd_cfg_jj('detect')

      end if

    3 if(Jintra(i).lt.js_max(i)) then
       Jintra(i)=Jintra(i)+2
       go to 2
      else
       if(i.eq.i1) go to 4
       i=i-1; go to 3
      end if

    4 Continue

      End  Subroutine Sum_Iterm






!====================================================================
      Integer Function Jparity(n,l,q)
!====================================================================
!     defines the parity for given configuration {l^q}
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: n
      Integer, intent(in) :: l(n),q(n)
      Integer :: m,i
      m=0;   Do i=1,n; m=m+l(i)*q(i);  End do
      Jparity=(-1)**m
      End Function Jparity




