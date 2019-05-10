!--------------------------------------------------------------
!     DBS LIBRARY
!--------------------------------------------------------------
!     contains routines for B-spline calculations in case
!     of different expansions for large and small components
!--------------------------------------------------------------
!     MODULES:
!
!     DBS_nuclear.f90
!     DBS_grid.f90
!     DBS_gauss.f90
!     DBS_debug.f90
!     DBS_dhl.f90
!     DBS_integrals.f90
!     DBS_moments.f90
!
!--------------------------------------------------------------
!     ROUTINES:
!
!     bdcwf.f90
!     bvalue2.f90
!     bvalue_bm.f90
!     bxv.f90
!     convert.f90
!     convol.f90
!     convol_sk.f90
!     def_Vnucl.f90
!     density.f90
!     mrk_pppp.f90
!     mrk_pqpq.f90
!     mrk_qpqp.f90
!     mrk_qqqq.f90
!     msk_ppqq.f90
!     msk_pqqp.f90
!     quadr.f90
!     rk.f90
!     sum_amb.f90
!     tvm.f90
!     update_hs.f90
!     zbsplvd.f90
!     zinty.f90
!--------------------------------------------------------------



!====================================================================
      Module DBS_nuclear
!====================================================================
!     nuclear parameters
!--------------------------------------------------------------------
      Use zconst

      Implicit none

      Real(8) :: atomic_number  = zero
      Real(8) :: atomic_weight  = zero

      Character(20) :: nuclear = 'Fermi'

!     nuclear = point   - point nuclear
!             = uniform - uniform destribution
!             = Fermi   - Fermi destribution
!
!     uniform distribution:
!     V(r) = -3/2 Z/R [1-r^2/(3R^2],      r <  R
!          = -Z/r                         r >  R

      Real(8) :: r_uniform = zero
      Real(8) :: ro_uniform = zero

! ... Fermi distribution:   rho(r) = rho'/[1+exp((r-c)/a)]

      Real(8) :: a_fermi = zero
      Real(8) :: c_fermi = zero
      Real(8) :: t_fermi = 2.30   ! in fermi
      Real(8) :: rrms    = zero
      Real(8) :: ro_fermi = zero

! ... nuclear potential in different basis:

      Real(8), allocatable :: VR_nucl(:,:)
      Real(8), allocatable :: VB_nucl(:,:)
      Real(8), allocatable :: VF_nucl(:,:)
      Real(8), allocatable :: ZR_nucl(:,:)

! ... Nuclear spin (I) (in units of h/2*pi):

      Real(8) :: I_nuclear = zero

! ... Nuclear dipole moment (in nuclear magnetons):

      Real(8) :: D_nuclear = zero

! ... Nuclear quadrupole moment (in barns):

      Real(8) :: Q_nuclear = zero

      End Module DBS_nuclear


!======================================================================
      Subroutine Read_nuclear(nu,z,atw)
!======================================================================
!     Set up the nuclear parameters in module DBS_nuclear
!     z,awt - data from the calling program
!     if they are zero - try to get them from file unit "nu"
!----------------------------------------------------------------------
      Use DBS_nuclear

      Implicit none
      Integer :: nu, i, an
      Real(8) :: apar,cpar, z,atw, rms, A
      Character(200) :: atom, core, conf


! ... read parameters from file:

      if(nu.gt.0) then

       Call Read_rpar(nu,'atomic_number',atomic_number)
       Call Read_rpar(nu,'atomic_weight',atomic_weight)
       Call Read_apar(nu,'nuclear',nuclear)
       Call Read_rpar(nu,'r_uniform',r_uniform)
       Call Read_rpar(nu,'a_fermi',a_fermi)
       Call Read_rpar(nu,'c_fermi',c_fermi)
       Call Read_rpar(nu,'t_fermi',t_fermi)
       Call Read_rpar(nu,'rrms',rrms)
       if(z.le.0.d0) z = atomic_number
       atw = atomic_weight

      end if

      if(z.le.0.d0) Stop 'Stop in Read_nuclear:  Z <= 0'
      an = NINT(z)
      if(an.lt.1.or.an.gt.104) Stop 'Stop in Read_nuclear:  Z out of range'
      Call Def_atom(an,atom,A,rms,core,conf)

      if(atw.le.0.d0) atw = A
      if(rrms.eq.0.d0) rrms = rms
      if(rrms.lt.0.d0) then
       if(an.le.90)  then
         rrms=(0.836d0*atw**(1.d0/3.d0)+0.570d0)
       else
         rrms=(0.77d0*atw**(1.d0/3.d0)+0.980d0)
       end if
      end if

      if(nuclear.eq.'Fermi'.and.rrms.lt.2.d0) then
       rrms = 2.d0
       write(*,'(70("_"))')
       write(*,*) 'rrms is change to 2.0 because GETCPR routine from GRASP'
       write(*,*) 'is nor working for small values of rrms'
       write(*,'(70("_")/)')
      end if  

      atomic_number = z
      atomic_weight = atw

! ... guess for nucler radius:

      if(nuclear.eq.'uniform'.and.r_uniform.eq.0.d0) then
       r_uniform = rrms * fermi_in_cm / bohr_radius_in_cm
      end if

! ... check the Fermi parameters:

      if(nuclear.eq.'Fermi') then

       apar = a_fermi / fermi_in_cm * bohr_radius_in_cm
       cpar = c_fermi / fermi_in_cm * bohr_radius_in_cm

       if(a_fermi.eq.0.d0) then
        apar = t_fermi/(four*LOG(three))
        a_fermi = apar * fermi_in_cm / bohr_radius_in_cm
       end if

       if(c_fermi.eq.0.d0) then
        CALL GETCPR(rrms,apar,cpar)
        c_fermi = cpar * fermi_in_cm / bohr_radius_in_cm
       end if
      end if

      End Subroutine Read_nuclear


!======================================================================
      Real(8) Function Z_nuclear(r)
!======================================================================
!     Nucleus charge distribution
!----------------------------------------------------------------------
      Use DBS_nuclear, b => r_uniform, a => a_fermi, c => c_fermi, &
                       Z => atomic_number
      Implicit none
      Real(8), intent(in) :: r
      Real(8) :: ro

      if(r.le.zero) then
       Z_nuclear = zero
      elseif(nuclear.eq.'point') then       ! point nucleus
       Z_nuclear = zero
       if(r.le.0.0219)  Z_nuclear = Z       ! after GRASP
      elseif(nuclear.eq.'uniform') then     ! uniform distribution
       ro_uniform = Z/(4*pi/3*b**3)
       Z_nuclear = zero
       if(r.le.b)  Z_nuclear = ro
      elseif(nuclear.eq.'Fermi') then       ! Fermi distribution
       if(ro_fermi.eq.zero) Call  DEF_ro_fermi
       Z_nuclear = ro_fermi/(one+exp((r-c)/a))
      else
        Stop 'Z_nuclear: unknown charge distribution '
      end if

      End Function Z_nuclear


!======================================================================
      Real(8) Function V_nuclear(r)
!======================================================================
!     Nucleus potential at point "r" for different charge distributions
!     F. A. Parpia and A. K. Mohanty
!     "Relativistic basis set calculations for atoms with Fermi nuclei"
!     Phys Rev A46,3735 (1992)
!----------------------------------------------------------------------
      Use DBS_nuclear, Z => atomic_number, b => r_uniform, &
                       a => a_fermi, c => c_fermi

      Implicit none
      Real(8), intent(in) :: r
      Real(8) :: V, N, ra,rc, ac,ac2,ac3, s,s2,s3, p
      Real(8), external :: SKFUN

      if(r.le.zero) then

       V = zero

      elseif(nuclear.eq.'point') then       ! point nucleus

       V = -Z / r

      elseif(nuclear.eq.'uniform') then     ! uniform distribution

       if(r.gt.b) then
        V = -Z / r
       else
        V = -(three*Z)/(two*b)*(one - r*r/(three*b*b))
       end if

      elseif(nuclear.eq.'Fermi') then       ! Fermi distribution

        ac = a/c; ac2=ac*ac; ac3=ac2*ac
        p = pi*pi*ac2; s = six*ac3*SKFUN(3,-c/a)
        N = one + p - s
        rc = r/c

       if(r.le.c) then
        V = -Z / r
        ra = (r-c)/a; s2 = ac2*SKFUN(2,ra); s3 = ac3*SKFUN(3,ra)
        V = (three - rc*rc + p)/two - three*s2 - s/rc + six/rc*s3
        V = -(Z*V)/(N*c)
       else
        ra = (c-r)/a; s2 = ac2*SKFUN(2,ra); s3 = ac3*SKFUN(3,ra)
        V = N + six*s3 + three*rc*s2
        V = -Z/r * (V/N)
       end if

      else

        Stop 'V_nuclear: unknown charge distribution '

      end if

      V_nuclear = V

      End Function V_nuclear


!=======================================================================
      Real(8) Function SKFUN (k,x)
!=======================================================================
!              infinity   (-1)^n e^nx
!      S (x) =   Sum      -----------
!       k        n=1          n^k
!======================================================================
      Implicit none
      Integer, intent(in) :: k
      Real(8), intent(in) :: x
      Integer :: n
      Real(8) :: eps, base, dnum, dk, delta

      SKFUN = 0.d0
      if(x.lt.-500d0) Return

      eps = 10.d0*epsilon(1.d0)
      BASE = -exp(x)
      DNUM = BASE
      SKFUN = BASE
      n = 1
      Do
       n = n + 1
       DNUM = DNUM*BASE
       dk = n; dk = dk ** k
       DELTA = DNUM/dk
       SKFUN = SKFUN+DELTA
       if(abs(DELTA/SKFUN).le.eps) Exit
      End do

      End Function SKFUN


!=====================================================================
      Real(8) Function ESTRMS (APARM,CPARM)
!=====================================================================
!     Determines the root mean square radius for a Fermi nucleus with
!     given  the parameters `c' (CPARM) and `a' (APARM).
!     Based on F. A. Parpia and A. K. Mohanty  Phys Rev A 46,3735(1992)
!     and W. R. Johnson, Lectures on Atomic Physics (2007)
!     Call(s) to: SKFUN.
!---------------------------------------------------------------------
      Use zconst

      Implicit none
      Real(8), intent(in) :: APARM,CPARM
      Real(8) :: A,P,X,DN,DD
      Real(8), external :: SKFUN

      A =  APARM/CPARM
      P =  PI*A
      X = -CPARM/APARM
      DN = 1.d0 + 10.d0/3.d0*P**2 + 7.d0/3.d0* P**4  &
               - 120.d0*A**5 * SKFUN(5,X)
      DD = 1.d0 + P**2 - 6.d0*A**3 * SKFUN(3,X)
      ESTRMS = CPARM * SQRT(3.d0*DN/5.d0/DD)

      End Function ESTRMS


!=====================================================================
      SUBROUTINE GETCPR (RRMS,APARM,CPARM)
!=====================================================================
!     Determines the parameter `c' (CPARM) for a Fermi nucleus,  given
!     the root mean square radius (RRMS) and the parameter `a' (APARM).
!     Call(s) to: ESTRMS.
!----------------------------------------------------------------------
      Use zconst

      Implicit none
      Real(8), intent(in) :: RRMS,APARM
      Real(8), intent(out) :: CPARM
      Real(8) :: ACCY,CPMIN,CPMAX,CPTRY,RMSTRY
      Real(8), external :: ESTRMS

! ... accuracy parameter

      ACCY = two*epsilon(one)

! ... bracket CPARM with a lower and upper limits:

      CPMIN = RRMS
      CPMAX = RRMS

    1 CPMIN = half*CPMIN

      if (ESTRMS(APARM,CPMIN) .GT. RRMS) GOTO 1

      CPMAX = RRMS
    2 CPMAX = two*CPMAX

      if (ESTRMS(APARM,CPMAX) .LT. RRMS) GOTO 2

! ... find CPARM by the method of bisection

    3 CPTRY = half*(CPMAX+CPMIN)

      RMSTRY = ESTRMS(APARM,CPTRY)

      if (RMSTRY .GT. RRMS) then
        CPMAX = CPTRY
      else
        CPMIN = CPTRY
      END if

      if ( ( (CPMAX-CPMIN)/(CPMAX+CPMIN) .GT. ACCY )   .AND. &
           ( ABS(RMSTRY-RRMS)/RRMS       .GT. ACCY ) )  GOTO 3

      CPARM = CPTRY

      End Subroutine GETCPR



!======================================================================
      Module DBS_grid
!======================================================================
!     the B-spline knot parameters
!----------------------------------------------------------------------
      Implicit none

      Integer :: grid_type = 1    ! type of grid

      Integer :: ks = 9           ! order of B-splines
      Integer :: ns = 0           ! number of splines
      Integer :: ms = 0           ! number of splines in (P,Q) basis
      Integer :: nv = 0           ! number of intervals ( = ns-ks+1 )

      Integer :: ml = 1           ! number of initial equal intervals
      Integer :: me = 0           ! number of intervals in the exponential region

      Real(8) :: hi   =  0.25d0   ! initial step parameter
      Real(8) :: he   =  0.25d0   ! exponential step factor
      Real(8) :: hmax =  1.d0     ! maximum step, t(ns+1) - t(ns)
      Real(8) :: rmax = 50.d0     ! input border radius
      Real(8) :: tmax =  0.d0     ! real border radius, t(ns+1)

! ... additional bases:

      Integer :: ksp=8, nsp=0     ! B-splines for p-functions
      Integer :: ksq=9, nsq=0     ! B-splines for q-functions

! ... knot sequence, t(1:ns+ks)

      Real(8), allocatable :: t(:)

! ... where to look for the grid:

      Character(80) :: AF_grid = 'knot.dat'

      End Module DBS_grid

!======================================================================
      SUBROUTINE def_grid (knot,name,z,atw)
!======================================================================
!     get input data for the grid and sets up the spline knots
!
!     knot -  if not empty, it defines file with grid parameters
!     name -  name of case; if is not empty, will be created name.knot,
!             file with grid parameters for the given case
!     z    -  nuclear charge
!     atw  -  atomic weight
!
!     if "knot" file is empty, program check file "name.knot";
!     if both files are absent, program creats grid file for given z;
!     if z=0, just example file created for z=1 and execution stops
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_nuclear

      Implicit none
      Character(*) :: knot, name
      Real(8) :: z,atw
      Integer :: nug=0
      Logical :: EX

! ... define the file name for input grid parameters:

      if(len_trim(knot).ne.0) then
       AF_grid = knot
       Inquire(file=AF_grid,exist=EX)
       if(EX) Call Find_free_unit(nug)
       if(INDEX(knot,'.bsw').eq.0) then
        if(nug.gt.0) Open(nug,file=AF_grid)
       else
        if(nug.gt.0) Open(nug,file=AF_grid,form='UNFORMATTED')
        grid_type = -2
       end if
      end if

      if(nug.eq.0.and.len_trim(name).ne.0) then
       AF_grid = trim(name)//'.knot'
       Call Find_free_unit(nug)
       Open(nug,file=AF_grid)
      end if

! ... create knot.dat as example:

      if(z.eq.0.and.nug.eq.0) then
       z = 1.d0; atw = 1.d0
       Call Read_nuclear(0,z,atw)
       Call mkgrid_01
       Call Write_knotdat
       Stop 'check knot.dat file first'
      end if

! ... check nuclear parameters first:

      if(grid_type.ne.-2) then
       Call Read_nuclear(nug,z,atw)
      else
       Call Read_nuclear(0,z,atw)
      end if

! ... read grid parameters from the file if any:

      if(nug.gt.0) then
       Call Read_ipar(nug,'grid_type',grid_type)
       Call Read_ipar(nug,'ns',ns)
       Call Read_ipar(nug,'ks',ks)
       Call Read_ipar(nug,'nv',nv)
       Call Read_ipar(nug,'ml',ml)
       Call Read_rpar(nug,'hi',hi)
       Call Read_rpar(nug,'he',he)
       Call Read_rpar(nug,'hmax',hmax)
       Call Read_rpar(nug,'rmax',rmax)
       Call Read_ipar(nug,'ksp',ksp)
       Call Read_ipar(nug,'ksq',ksq)
      end if

! ... read grid parameters from command line if any:

      Call Read_iarg('grid_type',grid_type)
      Call Read_iarg('ns',ns)
      Call Read_iarg('ks',ks)
      Call Read_iarg('nv',nv)
      Call Read_iarg('ml',ml)
      Call Read_rarg('hi',hi)
      Call Read_rarg('he',he)
      Call Read_rarg('hmax',hmax)
      Call Read_rarg('rmax',rmax)
      Call Read_iarg('ksp',ksp)
      Call Read_iarg('ksq',ksq)

      if(ksp.le.0.or.ksp.gt.ks) ksp=ks-1
      if(ksq.le.0.or.ksq.gt.ks) ksq=ks

! ... create the knot grid:

      Select case(grid_type)
       Case(1);      Call mkgrid_01
       Case(2);      Call mkgrid_02
       Case(-1);     Call read_grid(nug)
       Case(-2);     Call read_grid_bsw(nug)
       Case default; Stop 'Unknown grid_type'
      End Select

      Close(nug)

      ms = ns + ns
      nv = ns - ks + 1
      nsp=nv+ksp-1
      nsq=nv+ksq-1
      tmax = t(ns+1)

! ... record the grid parameters:

      if(grid_type.gt.0) Call Write_knotdat

      End Subroutine def_grid


!======================================================================
      Subroutine read_knot_dat
!======================================================================
      Character(40) :: knot = 'knot.dat', name = ' '
      Real(8) :: z = 0.d0, awt = 0.d0

      Call def_grid(knot,name,z,awt)

      End  Subroutine read_knot_dat


!======================================================================
      Subroutine read_grid (nu)
!======================================================================
!     read t-sequence from knot.dat;
!     it is option to work with old knot-sequences
!----------------------------------------------------------------------
      Use DBS_grid

      Implicit none
      Integer :: nu, i,j
      Integer, external :: Ifind_position

      if(ns.le.0)  Stop 'Stop in read_grid:  ns = 0'
      if(allocated(t)) Deallocate(t); Allocate(t(ns+ks))
      if(Ifind_position(nu,'grid points:').eq.0) &
       Stop 'Stop in read_grid:  can not find grid points'

      read(nu,*)
      read(nu,*)
      Do i=1,ns+ks; read(nu,*) j,t(j); End do

      End Subroutine read_grid


!======================================================================
      Subroutine read_grid_bsw (nu)
!======================================================================
!     read t-sequence from name.bsw unformated files;
!     it is option to work with old knot-sequences
!----------------------------------------------------------------------
      Use DBS_grid

      Implicit none
      Integer :: nu, i

      rewind(nu)
      read(nu) i,ns,ks
      if(allocated(t)) Deallocate(t); Allocate(t(ns+ks))
      rewind(nu)
      read(nu) i,ns,ks,t,ksp,ksq

      End Subroutine read_grid_bsw



!======================================================================
      Subroutine Write_knotdat
!======================================================================
!     create example of the knot.dat file
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_nuclear

      Implicit none
      Integer :: i, nug

      Call Find_free_unit(nug)
      Open(nug,file=AF_grid)
      rewind(nug)

      write(nug,'(a,i3,6x)') 'grid_type = ',grid_type
      write(nug,*)
      write(nug,'(a,i5,15x,a)') 'ks =',ks,'!  order  of splines'
      write(nug,'(a,i5,15x,a)') 'ns =',ns,'!  number of splines'
      write(nug,'(a,i5,15x,a)') 'nv =',nv,'!  number of intervals'
      write(nug,'(a,i5,15x,a)') 'ml =',ml,'!  number of initial equal-spaced intervals'
      write(nug,*)
      write(nug,'(a,i5,15x,a)') 'ksp=',ksp
      write(nug,'(a,i5,15x,a)') 'ksq=',ksq
      write(nug,*)
      write(nug,'(a,f16.8,2x,a)') 'hi =   ',hi,  '!  initial step for r*zg '
      write(nug,'(a,f16.8,2x,a)') 'he =   ',he,  '!  exponetial step size '
      write(nug,'(a,f16.8,2x,a)') 'hmax = ',hmax,'!  maximum step size for r, given'
      write(nug,'(a,f16.8,2x,a)') 'rmax = ',rmax,'!  maximum radius, given'
      write(nug,'(a,f16.8,2x,a)') 'tmax = ',tmax,'!  maximum radius, actual'
      write(nug,*)
      write(nug,'(a,f8.4)') 'atomic_number = ',atomic_number
      write(nug,'(a,f8.4)') 'atomic_weight = ',atomic_weight
      write(nug,'(a,f8.4)') 'rrms          = ',rrms
      write(nug,*)
      write(nug,'(a,a)') 'nuclear =  ',nuclear

      if(nuclear.eq.'Fermi') then
       write(nug,*)
       write(nug,'(a,d24.16,a)') 'c_fermi =   ',c_fermi,' !  in a.u.'
       write(nug,'(a,d24.16,a)') 'a_fermi =   ',a_fermi,' !  in a.u.'
      end if

      if(nuclear.eq.'uniform') then
       write(nug,*)
       write(nug,'(a,E16.8,a)') 'r_uniform = ',r_uniform
      end if

      write(nug,*)
      write(nug,'(a)') 'grid points:'
      write(nug,*)
      if(allocated(t)) then
       write(nug,'(i5,d24.16)') (i,t(i),i=1,ns+ks)
      End if
      write(nug,'(a)') '***'
      write(nug,*)
      write(nug,'(a)') 'possible grid-type:'
      write(nug,*)
      write(nug,'(a)') ' 1  -  semi-exponential, ns is derived from the parameters  '
      write(nug,'(a)') ' 2  -  exponential, for given ns and rmax '
      write(nug,'(a)') '-1  -  read from the given file '
      write(nug,'(a)') '-2  -  read from the bsw-file '

      Close(nug)

      End Subroutine Write_knotdat


!======================================================================
      Subroutine mkgrid_01
!======================================================================
!     sets up the knots for splines in semi-exponential grid:
!      1. first ml equal intervals "hi"  (this part depends on nuclear)
!      2. then exponantial grid with interval increasing as (1+he)
!      3. then again equal intervals "hmax" if any
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_nuclear, z => atomic_number

      Implicit none
      Integer :: i,ne
      Real(8) :: ti,tj,ts,ht,hs

! ... initial point:

      if(nuclear.eq.'point') then
       if(z.le.0.d0) Stop 'mkgrid_01: Z <= 0'
       hs=hi/z
      elseif(nuclear.eq.'Fermi') then
       ts = 4.d0 * log(3.d0) * a_fermi
       hs = ts*hi
      else
       hs = r_uniform*hi
      End if

      ti = hs*ml; ns = ks + ml

! ... "exponential" part of grid:

      ht = 1.d0 + he
      Do
       tj = ti*ht
       if(tj-ti.gt.hmax) Exit
       ns = ns + 1; ti = tj; ne = ns
       if(ti.gt.rmax) then; ns=ns-1; Exit; end if
      End do

! ... rest of interval with step = hmax

      if(ti.lt.rmax) then
       Do
        ti = ti + hmax
        if(ti.ge.rmax) Exit
        ns = ns + 1
       End do
      End if

      if(allocated(t)) Deallocate(t); Allocate(t(ns+ks))
      t = 0.d0

      Do i=ks+1,ks+ml;       t(i) = hs*(i-ks); End do
      Do i=ks+ml+1,ne;       t(i) = t(i-1) * ht; End do
      if(t(ne).lt.rmax) then
       Do i=ne+1,ns+1;       t(i) = t(i-1) + hmax; End do
      end if
      t(ns+2:ns+ks) = t(ns+1)
      tmax = t(ns+1)

      End Subroutine mkgrid_01


!======================================================================
      Subroutine mkgrid_02
!======================================================================
!     exponential grid for given ns and rmax;  he and hmax are derived.
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_nuclear, z => atomic_number

      Implicit none
      integer :: i
      Real(8) :: ts,ti

      ms=ns+ns; nv=ns-ks+1;  Allocate(t(ns+ks))

! ... initial point:

      if(nuclear.eq.'point') then
       ti=hi/z
      elseif(nuclear.eq.'Fermi') then
       ts = 4.d0 * log(3.d0) * a_fermi
       ti = ts*hi
      else
       ti = r_uniform*hi
      End if

! ... final point:

      tmax=rmax

      he = exp((log(tmax)-log(ti))/(ns-ks))

      t(1:ks)=0.d0; t(ks+1)=ti
      Do i=ks+2,ns+1; t(i)=t(i-1)*he; End do
      t(ns+1:ns+ks)=tmax

      hmax = t(ns+1)-t(ns)

      End Subroutine mkgrid_02





!======================================================================
     Module DBS_gauss
!======================================================================
!    defines the values of splines at the gaussian points for each
!    interval of a grid:  (1:nv+1,1:ks,1:ks)
!
!    bsp  (i,m,ith) - values of the i+ith-1 B-spline in interval i
!                     at gausian point m
!    bspd (i,m,ith) - corresponding values of first derivative
!
!    bspdd(i,m,ith) - corresponding values of second derivative
!
!    bsp(nv+1,1,.) and bspd(nv+1,1,.,.) - corresponding values at
!                                         last knot point (rmax)
!----------------------------------------------------------------------
     Implicit none

!... arrays for Gaussian-quadrature data (1:nv;1:ks)

     Real(8), allocatable :: gx(:)    ! gaussian points (1:ks)
     Real(8), allocatable :: gw(:)    ! gaussian points (1:ks)

     Real(8), allocatable :: gr (:,:) ! gr (i,m) - gaussian points m in the interval i
     Real(8), allocatable :: grm(:,:) ! grm(i,m) - reciprocal value of gr(i,m)
     Real(8), allocatable :: grw(:,:) ! grw(i,m) - corr. gaussian weights
     Real(8), allocatable :: ygw(:,:) ! working array

     Real(8), allocatable :: dbiatx(:,:,:) ! working array Used for call zbsplvd

!... standard B-splines: (keep for possible applications)

     Real(8), allocatable :: bsp(:,:,:)
     Real(8), allocatable :: bsq(:,:,:)
     Real(8), allocatable :: bspd(:,:,:)
     Real(8), allocatable :: bspdd(:,:,:)

!... basis for P-functions (ks=ksp):

     Real(8), allocatable, target :: pbsp(:,:,:)
     Real(8), allocatable :: pbsd(:,:,:)
     Real(8), allocatable :: dbip(:,:,:)
     Real(8), allocatable :: tp(:)
     Real(8), allocatable :: fpbs(:,:)
     Real(8), allocatable :: fpb_nucl(:,:)

!... basis for Q-functions (ks=ksq):

     Real(8), allocatable, target :: qbsp(:,:,:)
     Real(8), allocatable :: qbsd(:,:,:)
     Real(8), allocatable :: dbiq(:,:,:)
     Real(8), allocatable :: tq(:)
     Real(8), allocatable :: fqbs(:,:)
     Real(8), allocatable :: fqb_nucl(:,:)

!... mixed pq arrays:

     Real(8), allocatable :: fpqbs(:,:)
     Real(8), allocatable :: fqpbs(:,:)
     Real(8), allocatable :: fpqbsd(:,:)
     Real(8), allocatable :: fqpbsd(:,:)
     Real(8), allocatable :: fppqq(:,:)

     End Module DBS_gauss


!=======================================================================
      Real(8) Function memory_DBS_gauss(ns,ks)
!=======================================================================
!     return requred memory for module DBS_gauss (in Mb)
!------------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: ns,ks
      Integer :: nv,nw,ms
      nv = ns-ks+1; nw = nv + 1; ms = ns + ns
      memory_DBS_gauss = (nw*ks*ks*4*8 + nv*ks*4*8 + nv*ks*ks*8 +  &
                          nw*ks*ks*3*8 + ns*ns*2*8 +               &
                          nw*ks*ks*3*8 + ns*ns*2*8 +               &
                          ns*ns*4*8 + ms*ms*8 )/(1024d0*1024d0)
      End Function memory_DBS_gauss


!====================================================================
     Subroutine alloc_DBS_gauss
!====================================================================
! ... Allocate space and define the arrays in MODULE DBS_gauss
!--------------------------------------------------------------------
     Use DBS_grid
     Use DBS_gauss

     Implicit none
     Integer :: i,m,nw

     if(allocated(bsp)) &
      Deallocate(bsp,bsq,bspd,bspdd,gr,grm,grw,dbiatx)
     nw = nv+1
     Allocate(bsp (nw,ks,ks), bsq  (nw,ks,ks), &
              bspd(nw,ks,ks), bspdd(nw,ks,ks), &
              gr(nv,ks), grm(nv,ks), grw(nv,ks), ygw(nv,ks), &
              dbiatx(nv,ks,ks), gx(ks),gw(ks) )

     Call gauleg(0.d0,1.d0,gx,gw,ks)

     dbiatx = 0.d0
     Do m=1,ks
      gr (1:nv,m)= (t(1+ks:nv+ks)-t(ks:nv+ks-1))*gx(m)+t(ks:nv+ks-1)
      grm(1:nv,m)= 1.d0/gr(1:nv,m)
      grw(1:nv,m)= (t(1+ks:nv+ks)-t(ks:nv+ks-1))*gw(m)

      Call zbsplvd(ns,ks,nv,t,ks,nv,gr(1,m),3,dbiatx)

      bsp  (1:nv,m,1:ks) = dbiatx(1:nv,1:ks,1)
      bspd (1:nv,m,1:ks) = dbiatx(1:nv,1:ks,2)
      bspdd(1:nv,m,1:ks) = dbiatx(1:nv,1:ks,3)

      Do i = 1,ks
       bsq(1:nv,m,i) = bspd(1:nv,m,i) - grm(1:nv,m)*bsp(1:nv,m,i)
      End do

     End do

!... store values at the last knot

     Call zbsplvd(ns,ks,nv,t,ns,1,t(ns+1),3,dbiatx)

     bsp  (nw,:,:) = 0.d0
     bspd (nw,:,:) = 0.d0
     bspdd(nw,:,:) = 0.d0
     bsq  (nw,:,:) = 0.d0

     bsp  (nw,1,1:ks) = dbiatx(1,1:ks,1)
     bspd (nw,1,1:ks) = dbiatx(1,1:ks,2)
     bspdd(nw,1,1:ks) = dbiatx(1,1:ks,3)
     bsq  (nw,1,1:ks) = bspd(nw,1,1:ks) - bsp(nw,1,1:ks)/t(ns+1)

!... store also the values at the first knot

     Call zbsplvd(ns,ks,nv,t,ks,1,t(1),3,dbiatx)

     bsp  (nw,2,1:ks) = dbiatx(1,1:ks,1)
     bspd (nw,2,1:ks) = dbiatx(1,1:ks,2)
     bspdd(nw,2,1:ks) = dbiatx(1,1:ks,3)

     if(t(1).ne.0.d0) then
      bsq  (nw,2,1:ks) = bspd(nw,2,1:ks) - bsp(nw,2,1:ks)/t(1)
     else
      bsq  (nw,2,1:ks) = bspd(nw,2,1:ks) - bsp(nw,2,1:ks)/1.d-20
     end if

! ... case of different ks for p- and q-functions:

     Call Def_pbsp(ksp)
     Call Def_qbsp(ksq)
     Call Def_pq_arrays

     End Subroutine alloc_DBS_gauss


!====================================================================
     Subroutine Def_pbsp(k)
!====================================================================
! ... Allocate space and define the arrays specific for P-functions
!--------------------------------------------------------------------
     Use DBS_grid
     Use DBS_gauss

     Integer :: k, m, nw

     ksp=k; if(k.eq.0) ksp=ks; nsp=nv+ksp-1
     if(allocated(pbsp)) DeAllocate(pbsp); Allocate(pbsp(nv+1,ks,ksp))
     if(allocated(pbsd)) DeAllocate(pbsd); Allocate(pbsd(nv+1,ks,ksp))
     if(allocated(dbip)) DeAllocate(dbip); Allocate(dbip(nv,ksp,ksp))
     if(allocated(tp  )) DeAllocate(tp  ); Allocate(tp(nsp+ksp))
     if(allocated(fpbs)) DeAllocate(fpbs); Allocate(fpbs(ns,ns))
     if(allocated(fpb_nucl)) DeAllocate(fpb_nucl); Allocate(fpb_nucl(ns,ns))

     tp(1:ksp) = t(1)
     tp(ksp+1:nsp) = t(ks+1:ns)
     tp(nsp+1:nsp+ksp) = t(ns+1:ns+ksp)

     dbip = 0.d0
     Do m=1,ks
      Call zbsplvd(nsp,ksp,nv,tp,ksp,nv,gr(1,m),2,dbip)
      pbsp (1:nv,m,1:ksp) = dbip(1:nv,1:ksp,1)
      pbsd (1:nv,m,1:ksp) = dbip(1:nv,1:ksp,2)
     End do

     nw = nv + 1
     pbsp (nw,:,:) = 0.d0
     pbsd (nw,:,:) = 0.d0

!... store values at the last knot

     Call zbsplvd(nsp,ksp,nv,tp,nsp,1,tp(nsp+1),2,dbip)

     pbsp (nw,1,1:ksp) = dbip (1,1:ksp,1)
     pbsd (nw,1,1:ksp) = dbip (1,1:ksp,2)

     Call ZINTYm (nv,ks,ksp,ksp,pbsp,pbsp,grw,ns,fpbs)

     End Subroutine Def_pbsp


!====================================================================
     Subroutine Def_qbsp(k)
!====================================================================
! ... Allocate space and define the arrays specific for Q-functions
!--------------------------------------------------------------------
     Use DBS_grid
     Use DBS_gauss

     Integer :: k, m, nw

     ksq=k; if(k.eq.0) ksq=ks; nsq=nv+ksq-1

     if(allocated(qbsp)) DeAllocate(qbsp); Allocate(qbsp(nv+1,ks,ksq))
     if(allocated(qbsd)) DeAllocate(qbsd); Allocate(qbsd(nv+1,ks,ksq))
     if(allocated(dbiq)) DeAllocate(dbiq); Allocate(dbiq(nv,ksq,ksq))
     if(allocated(tq))   DeAllocate(tq);   Allocate(tq(nsq+ksq))
     if(allocated(fqbs)) DeAllocate(fqbs); Allocate(fqbs(ns,ns))
     if(allocated(fqb_nucl)) DeAllocate(fqb_nucl); Allocate(fqb_nucl(ns,ns))

     tq(1:ksq) = t(1)
     tq(ksq+1:nsq) = t(ks+1:ns)
     tq(nsq+1:nsq+ksq) = t(ns+1:ns+ksq)

     dbiq = 0.d0
     Do m=1,ks
      Call zbsplvd(nsq,ksq,nv,tq,ksq,nv,gr(1,m),2,dbiq)
      qbsp (1:nv,m,1:ksq) = dbiq(1:nv,1:ksq,1)
      qbsd (1:nv,m,1:ksq) = dbiq(1:nv,1:ksq,2)
     End do

     nw = nv + 1
     qbsp (nw,:,:) = 0.d0
     qbsd (nw,:,:) = 0.d0

!... store values at the last knot

     Call zbsplvd(nsq,ksq,nv,tq,nsq,1,tq(nsq+1),2,dbiq)

     qbsp (nw,1,1:ksq) = dbiq (1,1:ksq,1)
     qbsd (nw,1,1:ksq) = dbiq (1,1:ksq,2)

     Call ZINTYm (nv,ks,ksq,ksq,qbsp,qbsp,grw,ns,fqbs)

     End Subroutine Def_qbsp


!====================================================================
     Subroutine Def_pq_arrays
!====================================================================
! ... Allocate space and define the B-splines arrays
! ... connected  P- and Q-functions:
!--------------------------------------------------------------------
     Use DBS_grid
     Use DBS_gauss

     ygw = grm * grw
     if(Allocated(fpqbs)) DeAllocate(fpqbs); Allocate(fpqbs(ns,ns))
     Call ZINTYm (nv,ks,ksp,ksq,pbsp,qbsp,ygw,ns,fpqbs)
     if(Allocated(fqpbs)) DeAllocate(fqpbs); Allocate(fqpbs(ns,ns))
     Call ZINTYm (nv,ks,ksq,ksp,qbsp,pbsp,ygw,ns,fqpbs)

     if(Allocated(fpqbsd)) DeAllocate(fpqbsd); Allocate(fpqbsd(ns,ns))
     Call ZINTYm (nv,ks,ksp,ksq,pbsp,qbsd,grw,ns,fpqbsd)
     if(Allocated(fqpbsd)) DeAllocate(fqpbsd); Allocate(fqpbsd(ns,ns))
     Call ZINTYm (nv,ks,ksq,ksp,qbsp,pbsd,grw,ns,fqpbsd)

     if(Allocated(fppqq)) DeAllocate(fppqq); Allocate(fppqq(ms,ms))
     fppqq = 0.d0
     fppqq(1:ns,1:ns)=fpbs(1:ns,1:ns)
     fppqq(ns+1:ms,ns+1:ms)=fqbs(1:ns,1:ns)

     End Subroutine Def_pq_arrays


!======================================================================
      Module DBS_debug
!======================================================================
      Implicit none

      Real(8) :: time_convol = 0.d0
      Integer :: ic_convol =0

      Real(8) :: time_density = 0.d0
      Integer :: ic_density =0

      Real(8), external :: RRTC

      End Module DBS_debug





!====================================================================
      Module DBS_dhl_pq
!====================================================================
!     contains the spline representation of Dirac hamiltonian
!      <P|V_nucl|P>     c<P|-d/dr+k/r|Q>
!     c<Q|d/dr+k/r|P>   <Q|V_nucl-2c^2|Q>
!--------------------------------------------------------------------
      Implicit none

      Integer :: k_dhl = 0               ! k-value for HD-operator
      Real(8), allocatable :: dhl(:,:)

! ... Bloch operator:

      Integer :: mbloch  =   1    !  flag to include
      Real(8) :: RA = 0.d0        !  boder radius (=tmax)
      Real(8) :: RB = 0.d0        !  Bloch b-parameter
      Real(8) :: RN = 0.d0        !  Q(a)/P(a)
      Real(8) :: pnu = 1.d0       !  Q(a)/P(a)

      End Module DBS_dhl_pq

!======================================================================
      Subroutine alloc_dhl_pq(m)
!======================================================================
!     allocate or deallocate dhl array
!----------------------------------------------------------------------
      Use DBS_dhl_pq

      Integer :: m
      if(m.eq.0) then
       if(allocated(dhl)) Deallocate(dhl)
       k_dhl=0
      else
       if(allocated(dhl)) Deallocate(dhl)
       Allocate(dhl(m,m))
       k_dhl=0; dhl = 0.d0
      end if

      End Subroutine alloc_dhl_pq

!======================================================================
      Subroutine MAT_dhl_pq(k)
!======================================================================
!     Builds full matrix for the coulomb problem with kappa=k
!     from the symmetric forms for elementary operators
!----------------------------------------------------------------------
      USE zconst, only: c_au
      USE DBS_grid
      USE DBS_gauss
      USE DBS_dhl_pq

      Implicit none
      Real(8) :: sa,sb
      Integer :: k,i,j

      if(k_dhl.eq.k) Return
      if(.not.allocated(dhl)) Call Alloc_dhl_pq(ms)

! ... matrix of (1:1) = <B|V|B>
      dhl(1:ns,1:ns) = fpb_nucl
! ... matrix of (2:2) = <B|V|B> - 2*c^2*<B|B>
      sb = -2.d0 * c_au * c_au
      dhl(ns+1:ms,ns+1:ms) = fqb_nucl + sb*fqbs
! ... matrix of (2:1) = c <B | D+ | B>
      sa = c_au; sb = dble(k)*c_au
      dhl(ns+1:ms,1:ns) = sa*fqpbsd + sb*fqpbs
! ... matrix of (1:2) = c <B | D- | B>
      sa = -c_au; sb = dble(k)*c_au
      dhl(1:ns,ns+1:ms) = sa*fpqbsd + sb*fpqbs

      k_dhl = k

! ... symmetrize dhl by Bloch operator:

      if(mbloch.eq.1) then
       RN = (RB+k)/(2.d0*tmax*c_au)
       i = nsp; j = ns + nsq
       dhl(i,i) = dhl(i,i) - pnu*RN*c_au
       dhl(i,j) = dhl(i,j) + pnu*c_au
       dhl(j,i) = dhl(j,i) + (pnu-1.d0)*c_au
       dhl(j,j) = dhl(j,j) - (pnu-1.d0)/RN*c_au
      end if

      End Subroutine MAT_dhl_pq


!======================================================================
      Real(8) Function Vp_dhl(k1,v1,k2,v2)
!======================================================================
!     integral   <v1| dhl |v2>
!----------------------------------------------------------------------
      USE DBS_grid,   only: ms
      USE DBS_dhl_pq

      Implicit none
      Integer, Intent(in) :: k1,k2
      Real(8), Intent(in) :: v1(ms),v2(ms)
      Real(8) :: v(ms)
      Integer :: i,j

      Vp_dhl = 0.d0; if(k1.ne.k2) Return
      Call Mat_dhl_pq(k1)
      v = MATMUL(dhl,v2)
      Vp_dhl = SUM(v1*v)

      End Function Vp_dhl





!====================================================================
    Module DBS_moments
!====================================================================
!   contains moments defining as   <B_i|r^k|B_j>   over an interval.
!   These moments are used for calculation of two-electron integrals
!   according to the cell algorithm.
!--------------------------------------------------------------------
    Implicit none
    Integer :: kmk = -100, kk = 0
    Character(3) :: mtype = 'aaa'
    Real(8), allocatable :: rkd(:,:,:)
    Real(8), allocatable :: rkd1(:,:), rkd2(:,:), rkd3(:,:), rkd4(:,:)

    End Module DBS_moments

!====================================================================
    Subroutine alloc_DBS_moments
!====================================================================
!   allocate arrays in module DBS_moments and initialize the flags
!   mtype and kmk
!-------------------------------------------------------------------
    Use DBS_grid, only: nv,ks
    Use DBS_moments

    if(allocated(rkd)) Deallocate(rkd,rkd1,rkd2,rkd3,rkd4)
    kk = ks*ks
    Allocate(rkd(kk,kk,nv), rkd1(kk,nv),rkd2(kk,nv), &
                            rkd3(kk,nv),rkd4(kk,nv))
    rkd = 0.d0; rkd1 = 0.d0; rkd2 = 0.d0; rkd3 = 0.d0; rkd4 = 0.d0

    mtype='bbb';  kmk=-100

    End Subroutine alloc_DBS_moments


!====================================================================
    Subroutine dealloc_DBS_moments
!====================================================================
!   deallocate array in module DBS_moments and initialize the flags
!   mtype and kmk
!--------------------------------------------------------------------
    Use DBS_moments
    if(allocated(rkd)) Deallocate(rkd, rkd1,rkd2,rkd3,rkd4)
    mtype='aaa';  kmk=-100
    End Subroutine dealloc_DBS_moments


!=====================================================================
    Subroutine  moments (k,rkm,sym,dir)
!=====================================================================
!   Computes moments defining as <B_i|r^k|B_j> over an interval
!---------------------------------------------------------------------
!   on entry
!   --------
!      k      the power of moments
!      sym    's' or 'n' - symmetrical or non-symetrical case
!      div    'b' or 'd' - B_j or B'_j
!   on exit
!   -------
!      rkm    array of moments over all intervals.
!---------------------------------------------------------------------
      Use DBS_grid

      Implicit none
      Integer, intent(in) :: k
      Real(8), intent(out) :: rkm(ks*ks,nv)
      Character, intent(in) :: sym, dir
      ! local variables
      Integer :: iv, jk
      Real(8) :: hp
      Real(8) :: rv(ks*ks)

      jk = ks*(ks+1)/2; if(sym.eq.'n') jk = ks*ks

      if(me.eq.0) then     ! .. there is no exponential grid

       Do iv = 1,nv
        Call moment(iv,k,rv,sym,dir);  rkm(1:jk,iv) = rv(1:jk)
       End do

      else                 ! .. case of exponential grid

       ! .. the first non-exponential region
       DO iv = 1,ml+ks-1
        Call moment(iv,k,rv,sym,dir);  rkm(1:jk,iv) = rv(1:jk)
       End DO

       ! .. the exponential region -> using scaling law
       hp = (1.d0+he)**k; if(dir.eq.'b') hp = hp*(1.d0+he)
       DO iv=ml+ks,ml+me-ks+2
        rkm(1:jk,iv) = rkm(1:jk,iv-1)*hp
       End DO

       ! .. the last non-exponential region
       DO iv=ml+me-ks+3,nv
        Call moment(iv,k,rv,sym,dir);  rkm(1:jk,iv) = rv(1:jk)
       End DO

      End if

      End Subroutine moments


!=====================================================================
    Subroutine  moment (iv,k,rv,sym,dir)
!=====================================================================
!   Computes moment defining as <B_i|r^k|B_j> for given interval
!---------------------------------------------------------------------
!   on entry
!   --------
!      k      the power of moments
!      sym    's' or 'n' - symmetrical or non-symetrical case
!      div    'b' or 'd' - B_j or B'_j
!      iv     index of interval
!   on exit
!   -------
!      rv     array of moments for interval 'iv'
!---------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss

      Implicit none
      Integer, intent(in) :: k,iv
      Real(8), intent(out) :: rv(ks*ks)
      Character, intent(in) :: sym, dir
      ! local variables
      Integer :: i, j, ii, jj
      Real(8) :: bi(ks,ks),bj(ks,ks)

      gw = grw(iv,:)
      if( k > 0 ) gw = grw(iv,:) * gr (iv,:)**(+k)
      if( k < 0 ) gw = grw(iv,:) * grm(iv,:)**(-k)

      bi(:,:) = bsp(iv,:,:)
      if(dir.eq.'b') then
       Do j = 1,ks;  bj(:,j) = bsp(iv,:,j)*gw(:); End do
      else
       Do j = 1,ks;  bj(:,j) = bsq(iv,:,j)*gw(:); End do
      End if

      ii = 0
      Do i=1,ks; jj = 1; if(sym.eq.'s') jj = i
       Do j=jj,ks
        ii = ii + 1;  rv(ii) = SUM(bi(:,i)*bj(:,j))
       End do
      End do

      End Subroutine moment


!======================================================================
      Subroutine  moments_pp(k,kk,nv,rkm)
!======================================================================
!     Computes moments between P-functions
!     (cheking the scaling properties in future if lucky)
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: k,kk,nv
      Real(8), intent(out) :: rkm(kk,nv)
      Integer :: iv

      rkm = 0.d0
      Do iv = 1,nv; Call moment_pp(iv,k,rkm(1,iv)); End do

      End Subroutine moments_pp


!======================================================================
      Subroutine moment_pp(iv,k,rv)
!======================================================================
!     Computes moment between P-functions for interval "iv"
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss

      Implicit none
      Integer, intent(in) :: iv,k
      Real(8), intent(out) :: rv(*)
      Integer :: i, j, ii
      Real(8) :: b1(ks,ks),b2(ks,ks)

      gw = grw(iv,1:ks)
      if( k > 0 ) gw = grw(iv,1:ks) * gr (iv,1:ks)**(+k)
      if( k < 0 ) gw = grw(iv,1:ks) * grm(iv,1:ks)**(-k)

      b1(1:ks,1:ksp) = pbsp(iv,1:ks,1:ksp)
      Do j=1,ksp; b2(1:ks,j)=b1(1:ks,j)*gw(1:ks); End do

      ii = 0
      Do i=1,ksp; Do j=i,ksp
        ii=ii+1; rv(ii)=SUM(b1(1:ks,i)*b2(1:ks,j))
      End do; End do

      End Subroutine moment_pp


!======================================================================
      Subroutine  moments_qq(k,kk,nv,rkm)
!======================================================================
!     Computes moments between Q-functions
!     (cheking the scaling properties in future if lucky)
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: k,kk,nv
      Real(8), intent(out) :: rkm(kk,nv)
      Integer :: iv

      rkm = 0.d0
      Do iv = 1,nv; Call moment_qq(iv,k,rkm(1,iv)); End do

      End Subroutine moments_qq


!======================================================================
      Subroutine moment_qq(iv,k,rv)
!======================================================================
!     Computes moment between Q-functions for interval "iv"
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss

      Implicit none
      Integer, intent(in) :: iv,k
      Real(8), intent(out) :: rv(*)
      Integer :: i, j, ii
      Real(8) :: b1(ks,ks),b2(ks,ks)

      gw = grw(iv,1:ks)
      if( k > 0 ) gw = grw(iv,1:ks) * gr (iv,1:ks)**(+k)
      if( k < 0 ) gw = grw(iv,1:ks) * grm(iv,1:ks)**(-k)

      b1(1:ks,1:ksq) = qbsp(iv,1:ks,1:ksq)
      Do j=1,ksq; b2(1:ks,j)=b1(1:ks,j)*gw(1:ks); End do

      ii = 0
      Do i=1,ksq; Do j=i,ksq
        ii=ii+1; rv(ii)=SUM(b1(1:ks,i)*b2(1:ks,j))
      End do; End do

      End Subroutine moment_qq


!======================================================================
      Subroutine  moments_pq(k,kk,nv,rkm)
!======================================================================
!     Computes moments between P and Q-functions
!     (cheking the scaling properties in future if lucky)
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: k,kk,nv
      Real(8), intent(out) :: rkm(kk,nv)
      Integer :: iv

      rkm = 0.d0
      Do iv = 1,nv; Call moment_pq(iv,k,rkm(1,iv)); End do

      End Subroutine moments_pq


!======================================================================
      Subroutine moment_pq(iv,k,rv)
!======================================================================
!     Computes moment between P- and Q-functions for interval "iv"
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss

      Implicit none
      Integer, intent(in) :: iv,k
      Real(8), intent(out) :: rv(*)
      Integer :: i, j, ii
      Real(8) :: b1(ks,ks),b2(ks,ks)

      gw = grw(iv,1:ks)
      if( k > 0 ) gw = grw(iv,1:ks) * gr (iv,1:ks)**(+k)
      if( k < 0 ) gw = grw(iv,1:ks) * grm(iv,1:ks)**(-k)

      b1(1:ks,1:ksp) = pbsp(iv,1:ks,1:ksp)
      Do j=1,ksq; b2(1:ks,j)=qbsp(iv,1:ks,j)*gw(1:ks); End do

      ii = 0
      Do i=1,ksp; Do j=1,ksq
        ii=ii+1; rv(ii)=SUM(b1(1:ks,i)*b2(1:ks,j))
      End do; End do

      End Subroutine moment_pq


!======================================================================
      Subroutine  moments_qp(k,kk,nv,rkm)
!======================================================================
!     Computes moments between Q- and P-functions
!     (cheking the scaling properties in future if lucky)
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: k,kk,nv
      Real(8), intent(out) :: rkm(kk,nv)
      Integer :: iv

      rkm = 0.d0
      Do iv = 1,nv; Call moment_qp(iv,k,rkm(1,iv)); End do

      End Subroutine moments_qp


!======================================================================
      Subroutine moment_qp(iv,k,rv)
!======================================================================
!     Computes moment between Q- and P-functions for interval "iv"
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss

      Implicit none
      Integer, intent(in) :: iv,k
      Real(8), intent(out) :: rv(*)
      Integer :: i, j, ii
      Real(8) :: b1(ks,ks),b2(ks,ks)

      gw = grw(iv,1:ks)
      if( k > 0 ) gw = grw(iv,1:ks) * gr (iv,1:ks)**(+k)
      if( k < 0 ) gw = grw(iv,1:ks) * grm(iv,1:ks)**(-k)

      b1(1:ks,1:ksq) = qbsp(iv,1:ks,1:ksq)
      Do j=1,ksp; b2(1:ks,j)=pbsp(iv,1:ks,j)*gw(1:ks); End do

      ii = 0
      Do i=1,ksq; Do j=1,ksp
        ii=ii+1; rv(ii)=SUM(b1(1:ks,i)*b2(1:ks,j))
      End do; End do

      End Subroutine moment_qp



!====================================================================
      Module DBS_integrals
!====================================================================
!     contains the B-spline representation of two-electron integrals
!     rkb(i,j;i',j') in the symmetric column storage mode:
!                rkb(1:ns, 1:ns, 1:ks, 1:ks)
!     itype - character which indicates the type of integral
!     krk   - multipole index for the integral
!--------------------------------------------------------------------
      Implicit none

      Integer :: krk = -100
      Real(8), pointer :: rkb(:,:,:,:)
      Character(4) :: itype='aaaa'

      Integer :: kra_min = -100
      Integer :: kra_max = -100
      Integer :: ntype = 0
      Integer, allocatable :: irka(:,:)
      Real(8), allocatable, target :: rka(:,:,:,:,:,:)

      Real(8) :: memory_DBS_integrals = 0.d0

      End Module DBS_integrals


!====================================================================
      Subroutine alloc_DBS_integrals(ns,ks,kmin,kmax,ktype)
!====================================================================
! ... Allocates space for spline integrals:
!     1.  pppp
!     2.  qqqq
!     3.  pqpq
!     4.  qpqp
!--------------------------------------------------------------------
      Use DBS_integrals

      Implicit none
      Integer, intent(in) :: ns,ks,kmin,kmax,ktype

      if(allocated(irka)) Deallocate (irka,rka)
      Allocate(irka(kmin:kmax,ktype))
      Allocate(rka(ns,ns,ks,ks,kmin:kmax,ktype))
      kra_min = kmin
      kra_max = kmax
      irka = -100
      ntype = ktype

      if (associated(rkb)) Nullify(rkb)
      itype='bbbb';  krk=-100

      memory_DBS_integrals = ((kmax-kmin+1)*ntype*4 + &
       ns*ns*ks*ks*(kmax-kmin+1)*ktype*8)/(1024d0*1024d0)

      Call alloc_DBS_moments

      End Subroutine alloc_DBS_integrals


!====================================================================
      Subroutine dealloc_DBS_integrals
!====================================================================
! ... deAllocates arrays in module "spline_integrals"
!--------------------------------------------------------------------
      Use DBS_integrals

      if (associated(rkb)) nullify(rkb)
      itype='aaaa'; krk = -100
      if (allocated(irka)) Deallocate(irka,rka)
      kra_min = -100
      kra_max = -100
      ntype = 0
      memory_DBS_integrals = 0.d0

      End Subroutine dealloc_DBS_integrals





!====================================================================
      Module DBS_sk_integrals
!====================================================================
!     contains the B-spline representation of two-electron integrals
!     rkb(i,j;i',j') in the non-symmetric column storage mode:
!                rkb(1:ns, 1:ns, 1:2*ks-1, 1:2*ks-1)
!     itype - character which indicates the type of integral
!     krk   - multipole index for the integral
!--------------------------------------------------------------------
      Implicit none

      Integer :: krk = -100
      Real(8), pointer :: rkb(:,:,:,:)
      Character(4) :: itype='aaaa'

      Integer :: kra_min = -100
      Integer :: kra_max = -100
      Integer :: ntype = 0
      Integer, allocatable :: irka(:,:)
      Real(8), allocatable, target :: rka(:,:,:,:,:,:)

      Real(8) :: memory_DBS_sk_integrals = 0.d0

      End Module DBS_sk_integrals


!====================================================================
      Subroutine alloc_DBS_sk_integrals(ns,ks,kmin,kmax,ktype)
!====================================================================
! ... Allocates space for spline integrals:
!     1. ppqq
!     2. pqqp
!--------------------------------------------------------------------
      Use DBS_sk_integrals

      Implicit none
      Integer, intent(in) :: ns,ks,kmin,kmax,ktype

      if(allocated(irka)) Deallocate (irka,rka)
      Allocate(irka(kmin:kmax,ktype))
      Allocate(rka(ns,ns,2*ks-1,2*ks-1,kmin:kmax,ktype))
      kra_min = kmin
      kra_max = kmax
      irka = -100
      ntype = ktype

      if (associated(rkb)) Nullify(rkb)
      itype='bbbb';  krk=-100

      memory_DBS_sk_integrals = ((kmax-kmin+1)*ntype*4 + &
       ns*ns*ks*ks*(kmax-kmin+1)*ktype*8)/(1024d0*1024d0)

      Call alloc_DBS_moments

      End Subroutine alloc_DBS_sk_integrals


!====================================================================
      Subroutine dealloc_DBS_sk_integrals
!====================================================================
! ... deallocates arrays in module "DBS_sk_integrals"
!--------------------------------------------------------------------
      Use DBS_sk_integrals

      if (associated(rkb)) nullify(rkb)
      itype='aaaa'; krk = -100
      if (allocated(irka)) Deallocate(irka,rka)
      kra_min = -100
      kra_max = -100
      ntype = 0
      memory_DBS_sk_integrals = 0.d0

      End Subroutine dealloc_DBS_sk_integrals





!=======================================================================
    SUBROUTINE bdcwf_pq(n,k,z,Pcoef,Qcoef)
!=======================================================================
!   computes the spline expansion coefficients for hydrogenic radial
!   function for quantum numbers (n,k) and effective nuclear charge z.
!
!   Bases for P and Q may be different.
!
!   SUBROUTINE(S) called:  dcwf, gaussj
!-----------------------------------------------------------------------
!   on entry
!   --------
!       n,k     orbital quantum numbers
!       z       effective nuclear charge
!   on exit
!   -------
!       Pcoef   the spline coefficients for large component P(nl;r)
!       Qcoef   the spline coefficients for small component Q(nl;r)
!-----------------------------------------------------------------------
    Use DBS_grid,  only: ns,ks,nv, nsp,ksp, nsq,ksq
    Use DBS_gauss, only: gr,grw, pbsp,qbsp, fpbs,fqbs

    Implicit none
    Integer, intent(in)  :: n, k
    Real(8), intent(in)  :: z
    Real(8), intent(out) :: Pcoef(ns),Qcoef(ns)
    ! .. Local variables
    Integer :: i,iv,ip,m
    Real(8) :: E,pr(nv,ks),qr(nv,ks),a(ns,ns)

    ! .. obtain the values of the radial function
    ! .. at all the gaussian points
    ! .. and multiply by the gaussian weights

    Do m=1,ks
     Call dcwf(n,k,z,E,nv,gr(1,m),pr(1,m),qr(1,m))
     pr(:,m) = pr(:,m)*grw(:,m)
     qr(:,m) = qr(:,m)*grw(:,m)
    End do

    ! .. form the vector of inner products of the radial function and the
    ! .. spline basis functions

    Pcoef = 0.d0
    Do iv = 1,nv; Do ip = 1,ksp; i = iv+ip-1
     Pcoef(i) = Pcoef(i) + SUM(pr(iv,:)*pbsp(iv,:,ip))
    End do; End do

    Qcoef = 0.d0
    Do iv = 1,nv; Do ip = 1,ksq; i = iv+ip-1
     Qcoef(i) = Qcoef(i) + SUM(qr(iv,:)*qbsp(iv,:,ip))
    End do; End do

    ! .. solve the system of equations for coef

    a(1:nsp-1,1:nsp-1)=fpbs(2:nsp,2:nsp)
    Call gaussj (a,nsp-1,ns,Pcoef(2),1,ns)
    Pcoef(1)=0.d0

    a(1:nsq-1,1:nsq-1)=fqbs(2:nsq,2:nsq)
    Call gaussj (a,nsq-1,ns,Qcoef(2),1,ns)
    Qcoef(1)=0.d0

    End Subroutine bdcwf_pq



!=======================================================================
    Real(8) Function bvalu2 (t, bcoef, ns, ks, x, jderiv)
!=======================================================================
!   This routine is a modification of de Boor's bvalue routine, modified
!   to Return the value at right endpoint continuous from the left,
!   instead of 0. It assumes the usual knot multiplicity ks at the right
!   endpoint.
!   It calculates the value at  x of the jderiv-th derivative of spline
!   from b-representation. The spline is taken to be continuous from the
!   right.
!   SUBROUTINES contained: interv
!-----------------------------------------------------------------------
!   on entry
!   --------
!       t      knot sequence, of length  ns+ks, assumed nondecreasing.
!       bcoef  coefficient sequence, of length  ns .
!       ns     length of  bcoef, assumed positive.
!       ks     order of the spline .
!
!             . . . W A R N I N G . . .
!       The restriction  ks <= kmax (=15)  is imposed
!       arbitrarily by the parameter statement defining dimensions
!       for aj, dm, dp  below, but is  NEVER CHECKED.
!
!       x      the point at which to evaluate the spline.
!       jderiv integer giving the order of the derivative to be evaluated
!              ASSUMED to be zero or positive.
!
!   on exit
!   -------
!   bvalu2 - the value of the (jderiv)-th derivative of  f  at  x .
!
!   method
!   ------
!   the nontrivial knot interval  (t(i),t(i+1))  containing  x  is lo-
!   cated with the aid of  INTERV. the  ks  b-coeffs of  f  relevant for
!   this interval are then obtained from  bcoef (or taken to be zero if
!   not explicitly available) and are then differenced  jderiv  times to
!   obtain the b-coeffs of  (d**jderiv)f  relevant for that interval.
!   precisely, with  j = jderiv, we have from x.(12) of the text that
!
!           (d**j)f  =  sum ( bcoef(.,j)*b(.,ks-j,t) )
!
!   where
!                  / bcoef(.),                     ,  j = 0
!                  /
!   bcoef(.,j)  =  / bcoef(.,j-1) - bcoef(.-1,j-1)
!                  / ----------------------------- ,  j > 0
!                  /    (t(.+ks-j) - t(.))/(ks-j)
!
!   then, we use repeatedly the fact that
!
!      sum ( a(.)*b(.,m,t)(x) )  =  sum ( a(.,x)*b(.,m-1,t)(x) )
!   with
!                   (x - t(.))*a(.) + (t(.+m-1) - x)*a(.-1)
!      a(.,x)  =    ---------------------------------------
!                   (x - t(.))      + (t(.+m-1) - x)
!
!   to write  (d**j)f(x)  eventually as a linear combination of b-splines
!   of order  1 , and the coefficient for  b(i,1,t)(x)  must then
!   be the desired number  (d**j)f(x). (see x.(17)-(19) of text).
!-----------------------------------------------------------------------
    Implicit none
    Integer, Intent(in) :: ns, ks, jderiv
    Real(8), Intent(in) :: x, bcoef(ns), t(ns+ks)
    Real(8) :: aj(ks), dm(ks), dp(ks)
    Real(8) :: fkmj
    Integer :: i,mflag,km1,jcmin,imk,nmi,jcmax,jc,j,jj,kmj,ilo

    bvalu2 = 0.d0

    if (jderiv >= ks) Return

! ... find  i  such that 1 <= i < ns+ks  and  t(i) < t(i+1) and
! ... t(i) <= x < t(i+1) . if no such i can be found,  x  lies
! ... outside the support of  the spline  f  and  bvalu2 = 0.
! ... (the asymmetry in this choice of i makes f rightcontinuous)

    if( x /= t(ns+1) .OR. t(ns+1) /= t(ns+ks) ) then
      Call interv1(t,ns+ks,x, i, mflag)
      if (mflag /= 0) Return
    else
      i = ns
    end if

! ... if ks = 1 (and jderiv = 0), bvalu2 = bcoef(i).

    km1 = ks - 1
    if ( km1 <= 0 ) then
      bvalu2 = bcoef(i)
      Return
    end if

! ... store the ks b-spline coefficients relevant for the knot interval
! ... (t(i),t(i+1)) in aj(1),...,aj(ks); compute dm(j) = x - t(i+1-j),
! ... dp(j) = t(i+j) - x, j=1,...,ks-1. Set any of the aj not obtainable
! ... from input to zero. set any t.s not obtainable equal to t(1) or
! ... to t(ns+ks) appropriately.

    jcmin = 1
    imk = i - ks
    if (imk < 0) then
      jcmin = 1 - imk
      do j=1,i
        dm(j) = x - t(i+1-j)
      end do
      do j=i,km1
        aj(ks-j) = 0.
        dm(j) = dm(i)
      end do
    else
      do j=1,km1
        dm(j) = x - t(i+1-j)
      end do
    end if

    jcmax = ks
    nmi = ns - i
    if (nmi < 0) then
      jcmax = ks + nmi
      do j=1,jcmax
        dp(j) = t(i+j) - x
      end do
      do j=jcmax,km1
        aj(j+1) = 0.
        dp(j) = dp(jcmax)
      end do
    else
      do j=1,km1
        dp(j) = t(i+j) - x
      end do
    end if
      do jc=jcmin,jcmax
        aj(jc) = bcoef(imk + jc)
      end do

! ... difference the coefficients  jderiv  times.

    if (jderiv /= 0) then
      Do j=1,jderiv
        kmj = ks-j
        fkmj = kmj
        ilo = kmj
        Do jj=1,kmj
          aj(jj) = ((aj(jj+1) - aj(jj))/(dm(ilo) + dp(jj)))*fkmj
          ilo = ilo - 1
        End do
      End do
    end if

! ... compute value at  x  in (t(i),t(i+1)) of jderiv-th derivative,
! ... given its relevant b-spline coeffs in aj(1),...,aj(ks-jderiv).

    if (jderiv /= km1) then
      Do j=jderiv+1,km1
        kmj = ks-j
        ilo = kmj
        Do jj=1,kmj
          aj(jj) = (aj(jj+1)*dm(ilo) + aj(jj)*dp(jj))/(dm(ilo)+dp(jj))
          ilo = ilo - 1
        End do
      End do
    end if
    bvalu2 = aj(1)

    END FUNCTION bvalu2

!=======================================================================
      Subroutine INTERV1 ( xt, lxt, x, left, mflag )
!=======================================================================
!
!     Computes  left = max( i ; 1 <= i <= lxt  .and.  xt(i) <= x )
!     which is the interval containing x .
!
!     A reformatted version of the de Boor routine
!
!      on entry
!      --------
!      xt  a Real sequence, of length lxt, assumed to be nondecreasing
!      lxt number of terms in the sequence xt .
!      x   the point whose location with respect to the sequence xt is
!          to be determined.
!
!      on exit
!      -------
!      left, mflag   integers, whose values are
!         1     -1   if               x .lt.  xt(1)
!         i      0   if   xt(i)  .le. x .lt. xt(i+1)
!        lxt     1   if  xt(lxt) .le. x
!
!      In particular,  mflag = 0 is the 'usual' case.  mflag .ne. 0
!      indicates that  x  lies outside the halfopen interval
!      xt(1) .le. y .lt. xt(lxt) . the asymmetric treatment of the
!      interval is due to the decision to make all pp functions cont-
!      inuous from the right. (left - ?)
!
!      ...  m e t h o d  ...
!
!  The program is designed to be efficient in the common situation where
!  it is called repeatedly, with  x  taken from an increasing or decrea-
!  sing sequence. this will happen, e.g., when a pp function is to be
!  graphed. the first guess for  left  is therefore taken to be the val-
!  ue Returned at the previous Call and stored in the  l o c a l  varia-
!  ble  ilo . a first check ascertains that  ilo .lt. lxt (this is nec-
!  essary since the present call may have nothing to do with the previ-
!  ous call). then, if  xt(ilo) .le. x .lt. xt(ilo+1), we set  left =
!  ilo  and are done after just three comparisons.
!     otherwise, we repeatedly double the difference  istep = ihi - ilo
!  while also moving  ilo  and  ihi  in the direction of  x , until
!                      xt(ilo) .le. x .lt. xt(ihi) ,
!  after which we use bisection to get, in addition, ilo+1 = ihi .
!  left = ilo  is then Returned.
!-----------------------------------------------------------------------

      Integer :: left,lxt,mflag,ihi,ilo,istep,middle
      Real(8) :: x,xt(lxt)
      Data ilo /1/

      ihi = ilo + 1
      if (ihi .lt. lxt)                 go to 20
         if (x .ge. xt(lxt))            go to 110
         if (lxt .le. 1)                go to 90
         ilo = lxt - 1
         ihi = lxt

   20 if (x .ge. xt(ihi))               go to 40
      if (x .ge. xt(ilo))               go to 100

!     ... now x .lt. xt(ilo) . decrease  ilo  to capture  x

      istep = 1
   31    ihi = ilo
         ilo = ihi - istep
         if (ilo .le. 1)                go to 35
         if (x .ge. xt(ilo))            go to 50
         istep = istep*2
                                        go to 31
   35 ilo = 1
      if (x .lt. xt(1))                 go to 90
                                        go to 50

!     ... now x .ge. xt(ihi) . increase  ihi  to capture  x

   40 istep = 1
   41    ilo = ihi
         ihi = ilo + istep
         if (ihi .ge. lxt)              go to 45
         if (x .lt. xt(ihi))            go to 50
         istep = istep*2
                                        go to 41
   45 if (x .ge. xt(lxt))               go to 110
      ihi = lxt

!     ... now xt(ilo) .le. x .lt. xt(ihi) . narrow the interval

   50 middle = (ilo + ihi)/2
      if (middle .eq. ilo)              go to 100

!     note. it is assumed that middle = ilo in case ihi = ilo+1 .

      if (x .lt. xt(middle))            go to 53
         ilo = middle
                                        go to 50
   53    ihi = middle
                                        go to 50

!    ... set output and Return

   90 mflag = -1
      left = 1
                                        Return
  100 mflag = 0
      left = ilo
                                        Return
  110 mflag = 1
      left = lxt
                                        Return
      End Subroutine INTERV1



!======================================================================
      Subroutine BVALUE_bm(ksm,a,y,BSP)
!======================================================================
!     Computes the function y(r) = SUM(i)  a_i * B_i(r)
!     in all gausian points including border values
!----------------------------------------------------------------------
      Use DBS_grid

      Implicit none
      Integer :: ksm,i,j,iv,ith,m
      Real(8), intent(in)  :: a(*), bsp(nv+1,ks,ksm)
      Real(8), intent(out) :: y(*)

      Do iv = 1,nv                 ! over intervals
       Do m = 1,ks                 ! over gausian points
        Do ith = 1,ksm             ! over B-splines in given interval
         i = iv+ith-1              ! B-spline index
         j = 1 + (iv-1)*ks + m     ! radial point index
         y(j) = y(j) + a(i)*bsp(iv,m,ith)
        End do
       End do
      End do

      m = nv*ks+2                  ! last point
      Do i=1,ksm
       y(m) = y(m) + bsp(nv+1,1,i) * a(ns-ks+i)
      End do

      m = 1                        ! first point
      Do i=1,ksm
       y(m) = y(m) + bsp(nv+1,2,i) * a(i)
      End do

      End Subroutine BVALUE_bm



!=======================================================================
    SUBROUTINE bxv(k,n,b,v,y)
!=======================================================================
!   Computes   y = b * v    where b is a symmetric, banded matrix,
!   in lower-band storage mode,  and v, y are vectors.
!-----------------------------------------------------------------------
!
!   on entry
!   --------
!       k       the number of diagonals
!       n       the order of the matrix
!       b       the symmetric, banded matrix in column storage mode
!       v       vector
!
!   on exit
!   -------
!       y       y = b*v
!
!-----------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n, k
    REAL(KIND=8), DIMENSION(n), INTENT(IN) ::v
    REAL(KIND=8), DIMENSION(n,k), INTENT(IN) ::b
    REAL(KIND=8), DIMENSION(n), INTENT(out) :: y

    ! .. Local variables

    INTEGER :: i, j, jp

! ...   contribution from central diagonal (jp=k)

        Do i=1,n
          y(i) = b(i,k)*v(i)
        End do

! ...   off diagonal

        Do jp = 1,k-1
         Do i = k-jp+1,n
          j = i-k+jp
          y(i) = y(i) + b(i,jp)*v(j)             ! sub_diagonals
          y(j) = y(j) + b(i,jp)*v(i)             ! super-diagonals
         End do
        End do

  END SUBROUTINE bxv


!======================================================================
      Subroutine Convert_pq(nsw,ksw,tw,cw,nsv,ksv,cv,bsp,fbs,ll,mm)
!======================================================================
!     This programs converts B-spline orbital cw(1:nsw) given
!     on the grid tw to the B-spline orbital cv(1:nsv) with
!     grid defined in module DBS_grid
!     first ll and last mm splines are removed
!----------------------------------------------------------------------
      Use DBS_grid,   only: ns,ks,nv
      Use DBS_gauss,  only: gr, grw

      Implicit none
      Integer, intent(in)  :: nsw,ksw, nsv,ksv, ll,mm
      Real(8), intent(in)  :: tw(nsw+ksw), cw(nsw), &
                              bsp(nv+1,ks,ksv), fbs(ns,ns)
      Real(8), intent(out) :: cv(ns)
      ! local variables
      Integer :: i,j,ip,iv
      Real(8) :: a(ns,ns), ygw(nv,ks)
      Real(8), external :: bvalu2

! ... evaluate the function in gaussian points:

      Do i=1,nv; Do j=1,ks
       ygw(i,j) = bvalu2(tw,cw,nsw,ksw,gr(i,j),0) * grw(i,j)
      End do;  End do

! ... form the vector of inner products of the radial function
! ... and the spline basis functions:

      cv = 0.d0
      Do iv = 1,nv; Do ip = 1,ksv; i = iv+ip-1
       cv(i) = cv(i) + SUM(ygw(iv,:)*bsp(iv,:,ip))
      End do; End do

! ... B-spline overlap matrix:

      a(1:nsv-ll-mm,1:nsv-ll-mm) = fbs(ll+1:nsv-mm,ll+1:nsv-mm)

! ... solve the equation:  a cv = <B|cw>

      Call gaussj (a,nsv-ll-mm,ns,cv(ll+1),1,1)

      End Subroutine Convert_pq


!======================================================================
      Subroutine convol(ns,ks,a,d,icase,sym_i,sym_j)
!======================================================================
!     convolutes the rkb(i,j,i',j') array of spline integrals
!     with density matrix d(:,:)
!
!     results in array a(:,:)
!
!     icase =  1  - convolution other 2 and 4 variables, RK(.a;.b)
!              2  - convolution other 1 and 3 variables, RK(a.;b.)
!              3  - convolution other 2 and 3 variables, RK(.a;b.)
!              4  - convolution other 1 and 4 variables, RK(a.;.b)
!
!     sym_i  ->  symmetry in respect of i,i' variables ('s','l','n')
!     sym_j  ->  symmetry in respect of j,j' variables ('s','l','n')
!
!     combination of sym_i and sym_j leads to different represantation
!     for a and d:   a(ns,ks),  a(ns,2*ks+1),  a(ns,ns)
!----------------------------------------------------------------------
      Use DBS_integrals, only: rkb
      Use DBS_debug

      Implicit none
      Integer, intent(in) :: ns,ks,icase
      Character, intent(in) :: sym_i,sym_j
      Real(8), intent(in ) :: d(ns,*)
      Real(8), intent(out) :: a(ns,*)
      ! local variables
      Integer :: i,j, ip,jp, imin,imax, jmin,jmax, ii,jj
      Real(8) :: c,t1,t2

      t1 = RRTC()

      if(icase.le.2) a(1:ns,1:ks)=0.d0
      if(icase.gt.2) a(1:ns,1:ns)=0.d0

      Select case(icase)
!----------------------------------------------------------------------
      Case(1)                                          !  I( . a ; . b)

      if(sym_i.eq.'s'.and.sym_j.eq.'s') then

        do ip=1,ks
        do i=1,ns-ip+1;   c=0.d0
        do jp=1,ks
        do j=1,ns-jp+1;   c=c+d(j,jp)*rkb(i,j,ip,jp)
        end do; end do;   a(i,ip)=c
        end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'n') then

        do ip=1,ks+ks-1;  imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
        do i =imin,imax;  c=0.d0
        do jp=1,ks+ks-1;  jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
        do j =jmin,jmax;  c=c+d(j,jp)*rkb(i,j,ip,jp)
        end do; end do;   a(i,ip)=c
        end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'l') then

        do ip=1,ks
        do i=ks+1-ip,ns;  c=0.d0
        do jp=1,ks
        do j=ks+1-jp,ns;  c=c+d(j,jp)*rkb(i,j,ip,jp)
        end do; end do;   a(i,ip)=c
        end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'n') then

        do ip=1,ks
        do i=ks+1-ip,ns;  c=0.d0
        do jp=1,ks+ks-1;  jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
        do j =jmin,jmax;  c=c+d(j,jp)*rkb(i,j,ip,jp)
        end do; end do;   a(i,ip)=c
        end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'l') then

        do ip=1,ks+ks-1;  imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
        do i =imin,imax;  c=0.d0
        do jp=1,ks
        do j=ks+1-jp,ns;  c=c+d(j,jp)*rkb(i,j,ip,jp)
        end do; end do;   a(i,ip)=c
        end do; end do

      end if
!----------------------------------------------------------------------
      Case(2)                                         !  I( a . ; b . )

      if(sym_i.eq.'s'.and.sym_j.eq.'s') then

        do jp=1,ks
        do j=1,ns-jp+1;   c=0.d0
        do ip=1,ks
        do i=1,ns-ip+1;   c=c+d(i,ip)*rkb(i,j,ip,jp)
        end do; end do;   a(j,jp)=c
        end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'n') then

        do jp=1,ks+ks-1;  jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
        do j =jmin,jmax;  c=0.d0
        do ip=1,ks+ks-1;  imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
        do i =imin,imax;  c=c+d(i,ip)*rkb(i,j,ip,jp)
        end do; end do;   a(j,jp)=c
        end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'l') then

        do jp=1,ks
        do j=ks+1-jp,ns;  c=0.d0
        do ip=1,ks
        do i=ks+1-ip,ns;  c=c+d(i,ip)*rkb(i,j,ip,jp)
        end do; end do;   a(j,jp)=c
        end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'n') then

        do jp=1,ks+ks-1;  jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
        do j =jmin,jmax;  c=0.d0
        do ip=1,ks
        do i=ks+1-ip,ns;  c=c+d(i,ip)*rkb(i,j,ip,jp)
        end do; end do;   a(j,jp)=c
        end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'l') then

        do jp=1,ks
        do j=ks+1-jp,ns;  c=0.d0
        do ip=1,ks+ks-1;  imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
        do i =imin,imax;  c=c+d(i,ip)*rkb(i,j,ip,jp)
        end do; end do;   a(j,jp)=c
        end do; end do

      end if

!----------------------------------------------------------------------
      Case(3);  a(1:ns,1:ns) = 0.d0                    ! I( . a ; b . )

      if(sym_i.eq.'s'.and.sym_j.eq.'s') then

       do jp=1,ks
       do j =1,ns-jp+1; jj=j+jp-1
       do ip=1,ks
       do i =1,ns-ip+1; ii=i+ip-1
          c = rkb(i,j,ip,jp);      a( i,jj) = a( i,jj) + c*d(ii, j)
          if(ip.gt.1)              a(ii,jj) = a(ii,jj) + c*d( i, j)
          if(jp.gt.1)              a( i, j) = a( i, j) + c*d(ii,jj)
          if(ip.gt.1.and.jp.gt.1)  a(ii, j) = a(ii, j) + c*d( i,jj)
       end do;  end do
       end do;  end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'l') then

        do jp=1,ks
        do j=ks+1-jp,ns; jj=j+jp-ks
        do ip=1,ks
        do i=ks+1-ip,ns; ii=i+ip-ks
          c = rkb(i,j,ip,jp);       a( i,jj) = a( i,jj) + c*d(ii, j)
          if(ip.ne.ks)              a(ii,jj) = a(ii,jj) + c*d( i, j)
          if(jp.ne.ks)              a( i, j) = a( i, j) + c*d(ii,jj)
          if(ip.ne.ks.and.jp.ne.ks) a(ii, j) = a(ii, j) + c*d( i,jj)
        end do; end do
        end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'n') then

       do jp=1,ks+ks-1; jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
       do j =jmin,jmax; jj=j+jp-ks
       do ip=1,ks+ks-1; imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
       do i =imin,imax; ii=i+ip-ks
          c = rkb(i,j,ip,jp);      a( i,jj) = a( i,jj) + c*d(ii, j)
       end do; end do
       end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'l') then

        do jp=1,ks
        do j=ks+1-jp,ns; jj=j+jp-ks
        do ip=1,ks+ks-1; imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
        do i =imin,imax; ii=i+ip-ks
          c = rkb(i,j,ip,jp);      a( i,jj) = a( i,jj) + c*d(ii, j)
          if(jp.ne.ks)             a( i, j) = a( i, j) + c*d(ii,jj)
        end do; end do
        end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'n') then

        do jp=1,ks+ks-1; jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
        do j =jmin,jmax; jj=j+jp-ks
        do ip=1,ks
        do i=ks+1-ip,ns; ii=i+ip-ks
          c = rkb(i,j,ip,jp);       a( i,jj) = a( i,jj) + c*d(ii, j)
          if(ip.ne.ks)              a(ii,jj) = a(ii,jj) + c*d( i, j)
        end do; end do
        end do; end do

      end if
!----------------------------------------------------------------------
      Case(4);  a(1:ns,1:ns) = 0.d0                    ! I( a . ; . b )

      if(sym_i.eq.'s'.and.sym_j.eq.'s') then

       do jp=1,ks
       do j =1,ns-jp+1; jj=j+jp-1
       do ip=1,ks
       do i =1,ns-ip+1; ii=i+ip-1
          c = rkb(i,j,ip,jp);      a(ii, j) = a(ii, j) + c*d( i,jj)
          if(ip.gt.1)              a( i, j) = a( i, j) + c*d(ii,jj)
          if(jp.gt.1)              a(ii,jj) = a(ii,jj) + c*d( i, j)
          if(ip.gt.1.and.jp.gt.1)  a( i,jj) = a( i,jj) + c*d(ii, j)
       end do;  end do
       end do;  end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'n') then

       do jp=1,ks+ks-1; jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
       do j =jmin,jmax; jj=j+jp-ks
       do ip=1,ks+ks-1; imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
       do i =imin,imax; ii=i+ip-ks
          c = rkb(i,j,ip,jp);      a(ii, j) = a(ii, j) + c*d( i,jj)
       end do; end do
       end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'l') then

        do jp=1,ks
        do j=ks+1-jp,ns; jj=j+jp-ks
        do ip=1,ks
        do i=ks+1-ip,ns; ii=i+ip-ks
          c = rkb(i,j,ip,jp);       a(ii, j) = a(ii, j) + c*d( i,jj)
          if(ip.ne.ks)              a( i, j) = a( i, j) + c*d(ii,jj)
          if(jp.ne.ks)              a(ii,jj) = a(ii,jj) + c*d( i, j)
          if(ip.ne.ks.and.jp.ne.ks) a( i,jj) = a( i,jj) + c*d(ii, j)
        end do; end do
        end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'l') then

        do jp=1,ks
        do j=ks+1-jp,ns; jj=j+jp-ks
        do ip=1,ks+ks-1; imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
        do i =imin,imax; ii=i+ip-ks
          c = rkb(i,j,ip,jp);       a(ii, j) = a(ii, j) + c*d( i,jj)
          if(jp.ne.ks)              a(ii,jj) = a(ii,jj) + c*d( i, j)
        end do; end do
        end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'n') then

        do jp=1,ks+ks-1; jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
        do j =jmin,jmax; jj=j+jp-ks
        do ip=1,ks
        do i=ks+1-ip,ns; ii=i+ip-ks
          c = rkb(i,j,ip,jp);       a(ii, j) = a(ii, j) + c*d( i,jj)
          if(ip.ne.ks)              a( i, j) = a( i, j) + c*d(ii,jj)
        end do; end do
        end do; end do

      end if
!----------------------------------------------------------------------
      Case default

       Stop 'convol: unknown case'

      End Select

      t2 = RRTC()
      ic_convol = ic_convol + 1
      if(icase.le.2) time_convol=time_convol + (t2-t1)

      End Subroutine convol


!======================================================================
      Subroutine convol_sk(ns,ks,a,d,icase,sym_i,sym_j)
!======================================================================
!     convolutes the rkb(i,j,i',j') array of spline integrals
!     with density matrix d(:,:)
!
!     results in array a(:,:)
!
!     icase =  1  - convolution other 2 and 4 variables, RK(.a;.b)
!              2  - convolution other 1 and 3 variables, RK(a.;b.)
!              3  - convolution other 2 and 3 variables, RK(.a;b.)
!              4  - convolution other 1 and 4 variables, RK(a.;.b)
!
!     sym_i  ->  symmetry in respect of i,i' variables ('s','l','n')
!     sym_j  ->  symmetry in respect of j,j' variables ('s','l','n')
!
!     combination of sym_i and sym_j leads to different represantation
!     for a and d:   a(ns,ks),  a(ns,2*ks+1),  a(ns,ns)
!----------------------------------------------------------------------
      Use DBS_sk_integrals, only: rkb
      Use DBS_debug

      Implicit none
      Integer, intent(in) :: ns,ks,icase
      Character, intent(in) :: sym_i,sym_j
      Real(8), intent(in ) :: d(ns,*)
      Real(8), intent(out) :: a(ns,*)
      ! local variables
      Integer :: i,j, ip,jp, imin,imax, jmin,jmax, ii,jj
      Real(8) :: c,t1,t2

      t1 = RRTC()

      if(icase.le.2) a(1:ns,1:ks)=0.d0
      if(icase.gt.2) a(1:ns,1:ns)=0.d0

      Select case(icase)
!----------------------------------------------------------------------
      Case(1)                                          !  I( . a ; . b)

      if(sym_i.eq.'s'.and.sym_j.eq.'s') then

        do ip=1,ks
        do i=1,ns-ip+1;   c=0.d0
        do jp=1,ks
        do j=1,ns-jp+1;   c=c+d(j,jp)*rkb(i,j,ip,jp)
        end do; end do;   a(i,ip)=c
        end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'n') then

        do ip=1,ks+ks-1;  imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
        do i =imin,imax;  c=0.d0
        do jp=1,ks+ks-1;  jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
        do j =jmin,jmax;  c=c+d(j,jp)*rkb(i,j,ip,jp)
        end do; end do;   a(i,ip)=c
        end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'l') then

        do ip=1,ks
        do i=ks+1-ip,ns;  c=0.d0
        do jp=1,ks
        do j=ks+1-jp,ns;  c=c+d(j,jp)*rkb(i,j,ip,jp)
        end do; end do;   a(i,ip)=c
        end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'n') then

        do ip=1,ks
        do i=ks+1-ip,ns;  c=0.d0
        do jp=1,ks+ks-1;  jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
        do j =jmin,jmax;  c=c+d(j,jp)*rkb(i,j,ip,jp)
        end do; end do;   a(i,ip)=c
        end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'l') then

        do ip=1,ks+ks-1;  imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
        do i =imin,imax;  c=0.d0
        do jp=1,ks
        do j=ks+1-jp,ns;  c=c+d(j,jp)*rkb(i,j,ip,jp)
        end do; end do;   a(i,ip)=c
        end do; end do

      end if
!----------------------------------------------------------------------
      Case(2)                                         !  I( a . ; b . )

      if(sym_i.eq.'s'.and.sym_j.eq.'s') then

        do jp=1,ks
        do j=1,ns-jp+1;   c=0.d0
        do ip=1,ks
        do i=1,ns-ip+1;   c=c+d(i,ip)*rkb(i,j,ip,jp)
        end do; end do;   a(j,jp)=c
        end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'n') then

        do jp=1,ks+ks-1;  jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
        do j =jmin,jmax;  c=0.d0
        do ip=1,ks+ks-1;  imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
        do i =imin,imax;  c=c+d(i,ip)*rkb(i,j,ip,jp)
        end do; end do;   a(j,jp)=c
        end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'l') then

        do jp=1,ks
        do j=ks+1-jp,ns;  c=0.d0
        do ip=1,ks
        do i=ks+1-ip,ns;  c=c+d(i,ip)*rkb(i,j,ip,jp)
        end do; end do;   a(j,jp)=c
        end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'n') then

        do jp=1,ks+ks-1;  jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
        do j =jmin,jmax;  c=0.d0
        do ip=1,ks
        do i=ks+1-ip,ns;  c=c+d(i,ip)*rkb(i,j,ip,jp)
        end do; end do;   a(j,jp)=c
        end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'l') then

        do jp=1,ks
        do j=ks+1-jp,ns;  c=0.d0
        do ip=1,ks+ks-1;  imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
        do i =imin,imax;  c=c+d(i,ip)*rkb(i,j,ip,jp)
        end do; end do;   a(j,jp)=c
        end do; end do

      end if

!----------------------------------------------------------------------
      Case(3);  a(1:ns,1:ns) = 0.d0                    ! I( . a ; b . )

      if(sym_i.eq.'s'.and.sym_j.eq.'s') then

       do jp=1,ks
       do j =1,ns-jp+1; jj=j+jp-1
       do ip=1,ks
       do i =1,ns-ip+1; ii=i+ip-1
          c = rkb(i,j,ip,jp);      a( i,jj) = a( i,jj) + c*d(ii, j)
          if(ip.gt.1)              a(ii,jj) = a(ii,jj) + c*d( i, j)
          if(jp.gt.1)              a( i, j) = a( i, j) + c*d(ii,jj)
          if(ip.gt.1.and.jp.gt.1)  a(ii, j) = a(ii, j) + c*d( i,jj)
       end do;  end do
       end do;  end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'l') then

        do jp=1,ks
        do j=ks+1-jp,ns; jj=j+jp-ks
        do ip=1,ks
        do i=ks+1-ip,ns; ii=i+ip-ks
          c = rkb(i,j,ip,jp);       a( i,jj) = a( i,jj) + c*d(ii, j)
          if(ip.ne.ks)              a(ii,jj) = a(ii,jj) + c*d( i, j)
          if(jp.ne.ks)              a( i, j) = a( i, j) + c*d(ii,jj)
          if(ip.ne.ks.and.jp.ne.ks) a(ii, j) = a(ii, j) + c*d( i,jj)
        end do; end do
        end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'n') then

       do jp=1,ks+ks-1; jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
       do j =jmin,jmax; jj=j+jp-ks
       do ip=1,ks+ks-1; imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
       do i =imin,imax; ii=i+ip-ks
          c = rkb(i,j,ip,jp);      a( i,jj) = a( i,jj) + c*d(ii, j)
       end do; end do
       end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'l') then

        do jp=1,ks
        do j=ks+1-jp,ns; jj=j+jp-ks
        do ip=1,ks+ks-1; imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
        do i =imin,imax; ii=i+ip-ks
          c = rkb(i,j,ip,jp);      a( i,jj) = a( i,jj) + c*d(ii, j)
          if(jp.ne.ks)             a( i, j) = a( i, j) + c*d(ii,jj)
        end do; end do
        end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'n') then

        do jp=1,ks+ks-1; jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
        do j =jmin,jmax; jj=j+jp-ks
        do ip=1,ks
        do i=ks+1-ip,ns; ii=i+ip-ks
          c = rkb(i,j,ip,jp);       a( i,jj) = a( i,jj) + c*d(ii, j)
          if(ip.ne.ks)              a(ii,jj) = a(ii,jj) + c*d( i, j)
        end do; end do
        end do; end do

      end if
!----------------------------------------------------------------------
      Case(4);  a(1:ns,1:ns) = 0.d0                    ! I( a . ; . b )

      if(sym_i.eq.'s'.and.sym_j.eq.'s') then

       do jp=1,ks
       do j =1,ns-jp+1; jj=j+jp-1
       do ip=1,ks
       do i =1,ns-ip+1; ii=i+ip-1
          c = rkb(i,j,ip,jp);      a(ii, j) = a(ii, j) + c*d( i,jj)
          if(ip.gt.1)              a( i, j) = a( i, j) + c*d(ii,jj)
          if(jp.gt.1)              a(ii,jj) = a(ii,jj) + c*d( i, j)
          if(ip.gt.1.and.jp.gt.1)  a( i,jj) = a( i,jj) + c*d(ii, j)
       end do;  end do
       end do;  end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'n') then

       do jp=1,ks+ks-1; jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
       do j =jmin,jmax; jj=j+jp-ks
       do ip=1,ks+ks-1; imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
       do i =imin,imax; ii=i+ip-ks
          c = rkb(i,j,ip,jp);      a(ii, j) = a(ii, j) + c*d( i,jj)
       end do; end do
       end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'l') then

        do jp=1,ks
        do j=ks+1-jp,ns; jj=j+jp-ks
        do ip=1,ks
        do i=ks+1-ip,ns; ii=i+ip-ks
          c = rkb(i,j,ip,jp);       a(ii, j) = a(ii, j) + c*d( i,jj)
          if(ip.ne.ks)              a( i, j) = a( i, j) + c*d(ii,jj)
          if(jp.ne.ks)              a(ii,jj) = a(ii,jj) + c*d( i, j)
          if(ip.ne.ks.and.jp.ne.ks) a( i,jj) = a( i,jj) + c*d(ii, j)
        end do; end do
        end do; end do

      elseif(sym_i.eq.'n'.and.sym_j.eq.'l') then

        do jp=1,ks
        do j=ks+1-jp,ns; jj=j+jp-ks
        do ip=1,ks+ks-1; imin=max(1,1+ks-ip); imax=min(ns,ns+ks-ip)
        do i =imin,imax; ii=i+ip-ks
          c = rkb(i,j,ip,jp);       a(ii, j) = a(ii, j) + c*d( i,jj)
          if(jp.ne.ks)              a(ii,jj) = a(ii,jj) + c*d( i, j)
        end do; end do
        end do; end do

      elseif(sym_i.eq.'l'.and.sym_j.eq.'n') then

        do jp=1,ks+ks-1; jmin=max(1,1+ks-jp); jmax=min(ns,ns+ks-jp)
        do j =jmin,jmax; jj=j+jp-ks
        do ip=1,ks
        do i=ks+1-ip,ns; ii=i+ip-ks
          c = rkb(i,j,ip,jp);       a(ii, j) = a(ii, j) + c*d( i,jj)
          if(ip.ne.ks)              a( i, j) = a( i, j) + c*d(ii,jj)
        end do; end do
        end do; end do

      end if
!----------------------------------------------------------------------
      Case default

       Stop 'convol: unknown case'

      End Select

      t2 = RRTC()
      ic_convol = ic_convol + 1
      if(icase.le.2) time_convol=time_convol + (t2-t1)

      End Subroutine convol_sk


!======================================================================
      Subroutine Def_Vnucl
!======================================================================
!     define nuclear potential for p- and q-splines
!---------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DBS_nuclear

      Implicit none
      Integer :: i,j
      Real(8), external :: V_nuclear

      if(allocated(VR_nucl)) Deallocate(VR_nucl)
      Allocate(VR_nucl(nv,ks))

      Do i=1,nv
       Do j=1,ks
        VR_nucl(i,j) = V_nuclear(gr(i,j)) * grw(i,j)
       End do
      End do

      Call ZINTYm (nv,ks,ksp,ksp,pbsp,pbsp,VR_nucl,ns,fpb_nucl)
      Call ZINTYm (nv,ks,ksq,ksq,qbsp,qbsp,VR_nucl,ns,fqb_nucl)

      End Subroutine Def_Vnucl

!======================================================================
      Subroutine Def_ro_fermi
!======================================================================
!     define nuclear potential for p- and q-splines
!---------------------------------------------------------------------
      Use DBS_nuclear, Z => atomic_number, a => a_fermi, c => c_fermi
      Use DBS_grid
      Use DBS_gauss

      Implicit none
      Integer :: i,j
      Real(8) :: S

      Do i=1,nv
       Do j=1,ks
        ygw(i,j) = one/(one+exp((gr(i,j)-c)/a)) * grw(i,j) * gr(i,j) * gr(i,j)
       End do
      End do
      S =  SUM(ygw(:,:))

      ro_fermi = Z / (S * 4 * pi)

      End Subroutine Def_ro_fermi


!======================================================================
      Subroutine Def_Znucl
!======================================================================
!     define nuclear charge destribution in gausian points
!---------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DBS_nuclear

      Implicit none
      Integer :: i,j
      Real(8), external :: Z_nuclear

      if(allocated(ZR_nucl)) Deallocate(ZR_nucl)
      Allocate(ZR_nucl(nv,ks))

      Do i=1,nv
       Do j=1,ks
        ZR_nucl(i,j) = Z_nuclear(gr(i,j)) * grw(i,j)
       End do
      End do

      End Subroutine Def_Znucl



!======================================================================
      Subroutine density(ns,ks,d,p1,p2,type)
!======================================================================
!     d - density of two vectors p1,p2
!     type = 's' - symmetric upper-band column mode, or 'u'
!            'l' - symmetric lower-band column mode
!            'n' - non-symmetric band column mode
!            'x' - full matrix case, default
!----------------------------------------------------------------------
      Use DBS_debug

      Implicit none
      Character :: type
      Integer :: ns,ks, i,j,imin,imax
      Real(8) :: p1(ns),p2(ns)
      Real(8) :: d(ns,*)
      Real(8) :: t1,t2

      t1 = RRTC()

      if(type.eq.'s'.or.type.eq.'u') then

        d(1:ns,1:ks) = 0.d0
        do i=1,ns; d(i,1)=p1(i)*p2(i); end do
        do j=2,ks; do i=1,ns-j+1
          d(i,j)=p1(i)*p2(i+j-1) + p1(i+j-1)*p2(i)
        end do; end do

      elseif(type.eq.'l') then

        d(1:ns,1:ks) = 0.d0
        do i=1,ns; d(i,ks)= p1(i)*p2(i); end do
        do j=1,ks-1; do i=ks+1-j,ns
          d(i,j)=p1(i)*p2(i+j-ks) + p1(i+j-ks)*p2(i)
        end do; end do

      elseif(type.eq.'n') then

        d(1:ns,1:ks+ks-1) = 0.d0
        do j=1,ks+ks-1
          imin=max0( 1, 1+ks-j)
          imax=min0(ns,ns+ks-j)
          do i=imin,imax; d(i,j)=p1(i)*p2(i+j-ks); end do
        end do

      else

        d(1:ns,1:ns)=0.d0
        do i=1,ns; do j=1,ns
          d(i,j)=p1(i)*p2(j)
        end do; end do

      end if

      t2 = RRTC()
      ic_density=ic_density+1
      time_density = time_density + (t2-t1)

      End Subroutine density





!======================================================================
      Subroutine mrk_pppp(k)
!======================================================================
!     Defines matrix of Rk integrals in the B-spline basis
!     by cell algorithm between P-functions
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DBS_moments
      Use DBS_integrals

      Implicit none
      Integer, intent(in) :: k
      Integer :: i,j, ii,jj, iv,jv, ih,jh, ihp,jhp, ip,jp
      Integer, external :: Icheck_rka
      Real(8) :: c

! ... check the need of calculations

      if(itype == 'aaaa' .or. k < kra_min .or. k > kra_max) &
       Call alloc_DBS_integrals(ns,ks,0,k,4)

      if(itype == 'pppp' .and. krk == k) Return
      if (associated(rkb)) nullify(rkb)
      rkb => rka(:,:,:,:,k,1)
      if(irka(k,1) == 1) then
       krk=k; itype = 'pppp'
       Return
      end if

! ... compute the spline moments:

      Call moments_pp(  k   ,kk,nv,rkd1)
      Call moments_pp(-(k+1),kk,nv,rkd2)
      Call diag_pppp(k)

! ... generate the rkb array

      rkb=0.d0

      DO jv=1,nv;    jj = 0
      DO jh = 1,ksp; j  = jv  + jh - 1
      DO jhp=jh,ksp; jp = jhp - jh + 1
                     jj = jj  + 1

      DO iv=1,nv;    ii = 0
      DO ih=  1,ksp; i  = iv  + ih - 1
      DO ihp=ih,ksp; ip = ihp - ih + 1
                     ii = ii  + 1

          if     ( iv < jv ) then;   c = rkd1(ii,iv)*rkd2(jj,jv)
          else if( iv > jv ) then;   c = rkd1(jj,jv)*rkd2(ii,iv)
          else;                      c = rkd(ii,jj,iv)
          end if

          rkb(i,j,ip,jp) = rkb(i,j,ip,jp) +  c

      END DO;  END DO;  END DO
      END DO;  END DO;  END DO

      krk=k; itype = 'pppp'; irka(k,1)=1

      END Subroutine mrk_pppp


!======================================================================
      Subroutine diag_pppp(k)
!======================================================================
!     Controls the scaling propeties for diagonal B-spline Rk-interals
!     (not implimented yet)
!----------------------------------------------------------------------
      Use DBS_grid, only: nv

      Implicit none
      Integer :: k,iv

      Do iv=1,nv; Call triang_pppp(k,iv); End do

      END Subroutine diag_pppp


!======================================================================
    Subroutine triang_pppp (k,iv)
!======================================================================
!   Returns the two-dimensional array of B-spline integrals
!               <B_i B_j|r^k/r^(k+1)|B_i' B_j'>
!   over the given triangle diagonal cell
!
!   On entry   iv  -  index of the diagonal cell
!   --------
!
!   On exit    rkd(.,.,iv) - arrays of Rk B-spline integrals for given
!   --------                 interval iv in the reduced-dimension mode
!
!   Calls:   gauleg, zbsplvd
!----------------------------------------------------------------------
    Use DBS_grid
    Use DBS_gauss
    Use DBS_moments

    Implicit none
    Integer, intent(in) :: k,iv
    Integer :: i,j, ip,jp, ii,jj, m, left,  ik
    Real(8) :: xbase
    Real(8) :: x(ks),w(ks),bi(ks)
    Real(8) :: bspTmp(ks,ksp)
    Real(8) :: INTP(ksp,ksp,ks)
    Real(8) :: a(ksp*(ksp+1)/2,ksp*(ksp+1)/2)

    left=iv+ksp-1;  xbase=tp(left)

!.. setup the gaussian points

    Call gauleg(0.d0,1.d0,x,w,ks)

    DO m=1,ks

! .. the absolute coordinate at the new gaussian point

      gx(:) = (gr(iv,m)-xbase)*x(:) + xbase

! .. the bspline values at the new gaussian points

      DO i=1,ks
       Call zbsplvd(nsp,ksp,nv, tp,left, 1,gx(i),1,dbip)
       bspTmp(i,1:ksp)= dbip(1,1:ksp,1)
      END DO

! .. and the corresponding gaussian weights

      gw(:) = (gr(iv,m)-xbase)*w(:)
      if(k>1) then;            gx(:) = gw(:)*gx(:)**k
      else if(k==1) then;      gx(:) = gw(:)*gx(:)
      else if(k==0) then;      gx(:) = gw(:)
      end if

!            / r(iv,m)                             k
! .. INT =  |      bsp(iv,:,j)(r) bsp(iv,:,jp)(r) r  dr
!           / r_iv

      Do j=1,ksp
       gw(:) = gx(:)*bspTmp(:,j)
       Do jp=j,ksp; INTP(j,jp,m)= SUM(gw(:)*bspTmp(:,jp));  END DO
      End do

    END DO    !  over m

! .. second integration

    if(k/=0) then;   gx(:) = grw(iv,:)*grm(iv,:)**(k+1)
    else;            gx(:) = grw(iv,:)*grm(iv,:)
    end if

    ii = 0;  DO i=1,ksp;  DO ip=i,ksp;  ii = ii+1

             bi(:) = pbsp(iv,:,i)*pbsp(iv,:,ip)*gx(:)

    jj = 0;  DO j=1,ksp;  DO jp=j,ksp;  jj = jj+1

             a(ii,jj) =  SUM(bi(:)*INTP(j,jp,:))

             END DO; END DO
             END DO; END DO

    ik = ksp*(ksp+1)/2;  rkd(1:ik,1:ik,iv) = a + TRANSPOSE(a)

    END Subroutine triang_pppp


!======================================================================
      Subroutine mrk_pqpq(k)
!======================================================================
!     Defines matrix of Rk integrals in the B-spline basis
!     by cell algorithm for mixed PQ case
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DBS_moments
      Use DBS_integrals

      Implicit none
      Integer, intent(in) :: k
      Integer :: i,j, ii,jj, iv,jv, ih,jh, ihp,jhp, ip,jp
      Integer, external :: Icheck_rka
      Real(8) :: c

! ... check the need of calculations

      if(itype == 'aaaa' .or. k < kra_min .or. k > kra_max) &
       Call alloc_DBS_integrals(ns,ks,0,k,4)

      if(itype == 'pqpq' .and. krk == k) Return
      if (associated(rkb)) nullify(rkb)
      rkb => rka(:,:,:,:,k,3)
      if(irka(k,3) == 1) then
       krk=k; itype = 'pqpq'
       Return
      end if

! ... compute the spline moments:

      Call moments_pp(  k   ,kk,nv,rkd1)
      Call moments_qq(-(k+1),kk,nv,rkd2)
      Call moments_qq(  k   ,kk,nv,rkd3)
      Call moments_pp(-(k+1),kk,nv,rkd4)

      Call diag_pqpq(k)

! ... generate the rkb array

      rkb=0.d0

      DO jv = 1,nv;   jj = 0
      DO jh = 1,ksq;  j  = jv  + jh - 1
      DO jhp=jh,ksq;  jp = jhp - jh + 1
                      jj = jj  + 1

      DO iv = 1,nv;   ii = 0
      DO ih = 1,ksp;  i  = iv  + ih - 1
      DO ihp=ih,ksp;  ip = ihp - ih + 1
                      ii = ii  + 1

          if     ( iv < jv ) then;   c = rkd1(ii,iv)*rkd2(jj,jv)
          else if( iv > jv ) then;   c = rkd3(jj,jv)*rkd4(ii,iv)
          else;                      c = rkd(ii,jj,iv)
          end if

          rkb(i,j,ip,jp) = rkb(i,j,ip,jp) +  c

      END DO;  END DO;  END DO
      END DO;  END DO;  END DO

      krk=k; itype = 'pqpq'; irka(k,3) = 1

      END Subroutine mrk_pqpq


!======================================================================
      Subroutine diag_pqpq(k)
!======================================================================
!     Controls the scaling propeties for diagonal B-spline Rk-interals
!     (not implimented yet)
!----------------------------------------------------------------------
      Use DBS_grid, only: nv

      Implicit none
      Integer :: k,iv

      Do iv=1,nv; Call triang_pqpq(k,iv); End do

      END Subroutine diag_pqpq


!======================================================================
      Subroutine triang_pqpq(k,iv)
!======================================================================
!     Returns the two-dimensional array of B-spline integrals
!               <A_i B_j|r^k/r^(k+1)|C_i' D_j'>
!     over the given triangle diagonal cell
!
!     On entry   iv  -  index of the diagonal cell
!     --------
!
!     On exit    rkd(.,.,iv) - arrays of Rk B-spline integrals for given
!     --------                 interval iv in the reduced-dimension mode
!
!     Calls:   gauleg, zbsplvd
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DBS_moments

      Implicit none
      Integer :: iv,i,j, ip,jp, ii,jj, m, left, k
      Real(8) :: xbase
      Real(8) :: x(ks),w(ks), bp(ks),bq(ks)
      Real(8) :: bspTmp(ks,ksp)
      Real(8) :: bspTmq(ks,ksq)
      Real(8) :: INTP(ksp,ksp,ks)
      Real(8) :: INTQ(ksq,ksq,ks)

! ... setup the gaussian points

      Call gauleg(0.d0,1.d0,x,w,ks)

! ... first integration:

      left=iv+ks-1;  xbase=t(left)

      DO m=1,ks

! .. the absolute coordinate at the new gaussian point

      gx(:) = (gr(iv,m)-xbase)*x(:) + xbase

! .. the bspline values at the new gaussian points

      DO i=1,ks
       Call zbsplvd(nsp,ksp,nv, tp,iv+ksp-1, 1,gx(i),1,dbip)
       bspTmp(i,1:ksp)= dbip(1,1:ksp,1)
      END DO

      DO i=1,ks
       Call zbsplvd(nsq,ksq,nv, tq,iv+ksq-1, 1,gx(i),1,dbiq)
       bspTmq(i,1:ksq)= dbiq(1,1:ksq,1)
      END DO

! ... and the corresponding gaussian weights

      gw(:) = (gr(iv,m)-xbase)*w(:)
      gx(:) = gw(:)*gx(:)**k

!            / r(iv,m)                             k
! ... INT =  |      bsp(iv,:,j)(r) bsp(iv,:,jp)(r) r  dr
!           / r_iv

      Do j=1,ksp
       gw(:) = gx(:)*bspTmp(:,j)
       Do jp=j,ksp; INTP(j,jp,m)= SUM(gw(:)*bspTmp(:,jp));  END DO
      End do

      Do j=1,ksq
       gw(:) = gx(:)*bspTmq(:,j)
       Do jp=j,ksq; INTQ(j,jp,m)= SUM(gw(:)*bspTmq(:,jp));  END DO
      End do

      END DO !  over m

! ... second integration

      gx(:) = grw(iv,:)*grm(iv,:)**(k+1)

             rkd(:,:,iv) = 0.d0

      ii=0;  DO i=1,ksp;  DO ip=i,ksp; ii = ii+1
              bp(:) = pbsp(iv,:,i)*pbsp(iv,:,ip)*gx(:)
      jj=0;  DO j=1,ksq;  DO jp=j,ksq; jj = jj+1
              rkd(ii,jj,iv) = rkd(ii,jj,iv) + SUM(bp(:)*INTQ(j,jp,:))
             END DO; END DO
             END DO; END DO

      ii=0;  DO i=1,ksq;  DO ip=i,ksq; ii = ii+1
              bq(:) = qbsp(iv,:,i)*qbsp(iv,:,ip)*gx(:)
      jj=0;  DO j=1,ksp;  DO jp=j,ksp; jj = jj+1
              rkd(jj,ii,iv) = rkd(jj,ii,iv) + SUM(bq(:)*INTP(j,jp,:))
             END DO; END DO
             END DO; END DO

      End Subroutine triang_pqpq


!======================================================================
      Subroutine mrk_qpqp(k)
!======================================================================
!     Defines matrix of Rk integrals in the B-spline basis
!     by cell algorithm for mixed QP case.
!----------------------------------------------------------------------
      USE DBS_grid
      USE DBS_gauss
      USE DBS_moments
      USE DBS_integrals

      Implicit none
      Integer, intent(in) :: k
      Integer :: i,j, ii,jj, iv,jv, ih,jh, ihp,jhp, ip,jp
      Integer, external :: Icheck_rka
      Real(8) :: c

! ... check the need of calculations

      if(itype == 'aaaa' .or. k < kra_min .or. k > kra_max) &
       Call alloc_DBS_integrals(ns,ks,0,k,4)

      if(itype == 'qpqp' .and. krk == k) Return
      if (associated(rkb)) nullify(rkb)
      rkb => rka(:,:,:,:,k,4)
      if(irka(k,4) == 1) then
       krk=k; itype = 'qpqp'
       Return
      end if

! ... compute the spline moments:

      Call moments_pp(  k   ,kk,nv,rkd1)
      Call moments_qq(-(k+1),kk,nv,rkd2)
      Call moments_qq(  k   ,kk,nv,rkd3)
      Call moments_pp(-(k+1),kk,nv,rkd4)

      Call diag_qpqp(k)

! ... generate the rkb array

      rkb=0.d0

      DO jv = 1,nv;   jj = 0
      DO jh = 1,ksp;  j  = jv  + jh - 1
      DO jhp=jh,ksp;  jp = jhp - jh + 1
                      jj = jj  + 1

      DO iv=1,nv;     ii = 0
      DO ih=  1,ksq;  i  = iv  + ih - 1
      DO ihp=ih,ksq;  ip = ihp - ih + 1
                      ii = ii  + 1

          if     ( iv < jv ) then;   c = rkd3(ii,iv)*rkd4(jj,jv)
          else if( iv > jv ) then;   c = rkd1(jj,jv)*rkd2(ii,iv)
          else;                      c = rkd(ii,jj,iv)
          end if

          rkb(i,j,ip,jp) = rkb(i,j,ip,jp) +  c

      END DO;  END DO;  END DO
      END DO;  END DO;  END DO

      krk=k; itype = 'qpqp'; irka(k,4) = 1

      End Subroutine mrk_qpqp


!======================================================================
      Subroutine diag_qpqp(k)
!======================================================================
!     Controls the scaling propeties for diagonal B-spline Rk-interals
!     (not implimented yet)
!----------------------------------------------------------------------
      Use DBS_grid, only: nv

      Implicit none
      Integer :: k,iv

      Do iv=1,nv; Call triang_qpqp(k,iv); End do

      END Subroutine diag_qpqp


!======================================================================
      Subroutine triang_qpqp(k,iv)
!======================================================================
!     Returns the two-dimensional array of B-spline integrals
!               <Q_i P_j|r^k/r^(k+1)|Q_i' P_j'>
!     over the given triangle diagonal cell
!
!     On entry   iv  -  index of the diagonal cell
!     --------
!
!     On exit    rkd(.,.,iv) - arrays of Rk B-spline integrals for given
!     --------                 interval iv in the reduced-dimension mode
!
!     Calls:   gauleg, zbsplvd
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DBS_moments

      Implicit none
      Integer :: iv,i,j, ip,jp, ii,jj, m, left, k
      Real(8) :: xbase
      Real(8) :: x(ks),w(ks), bp(ks),bq(ks)
      Real(8) :: bspTmp(ks,ksp)
      Real(8) :: bspTmq(ks,ksq)
      Real(8) :: INTP(ksp,ksp,ks)
      Real(8) :: INTQ(ksq,ksq,ks)

! ... setup the gaussian points

      Call gauleg(0.d0,1.d0,x,w,ks)

! ... first integration:

      left=iv+ks-1;  xbase=t(left)

      DO m=1,ks

! .. the absolute coordinate at the new gaussian point

      gx(:) = (gr(iv,m)-xbase)*x(:) + xbase

! .. the bspline values at the new gaussian points

      DO i=1,ks
       Call zbsplvd(nsp,ksp,nv, tp,iv+ksp-1, 1,gx(i),1,dbip)
       bspTmp(i,1:ksp)= dbip(1,1:ksp,1)
      END DO

      DO i=1,ks
       Call zbsplvd(nsq,ksq,nv, tq,iv+ksq-1, 1,gx(i),1,dbiq)
       bspTmq(i,1:ksq)= dbiq(1,1:ksq,1)
      END DO

! ... and the corresponding gaussian weights

      gw(:) = (gr(iv,m)-xbase)*w(:)
      gx(:) = gw(:)*gx(:)**k

!            / r(iv,m)                             k
! ... INT =  |      bsp(iv,:,j)(r) bsp(iv,:,jp)(r) r  dr
!           / r_iv

      Do j=1,ksp
       gw(:) = gx(:)*bspTmp(:,j)
       Do jp=j,ksp; INTP(j,jp,m)= SUM(gw(:)*bspTmp(:,jp));  END DO
      End do

      Do j=1,ksq
       gw(:) = gx(:)*bspTmq(:,j)
       Do jp=j,ksq; INTQ(j,jp,m)= SUM(gw(:)*bspTmq(:,jp));  END DO
      End do

      END DO    !  over m

! ... second integration

      gx(:) = grw(iv,:)*grm(iv,:)**(k+1)

             rkd(:,:,iv) = 0.d0

      ii=0;  DO i=1,ksp;  DO ip=i,ksp; ii = ii+1
              bp(:) = pbsp(iv,:,i)*pbsp(iv,:,ip)*gx(:)
      jj=0;  DO j=1,ksq;  DO jp=j,ksq; jj = jj+1
              rkd(jj,ii,iv) = rkd(jj,ii,iv) + SUM(bp(:)*INTQ(j,jp,:))
             END DO; END DO
             END DO; END DO

      ii=0;  DO i=1,ksq;  DO ip=i,ksq; ii = ii+1
              bq(:) = qbsp(iv,:,i)*qbsp(iv,:,ip)*gx(:)
      jj=0;  DO j=1,ksp;  DO jp=j,ksp; jj = jj+1
              rkd(ii,jj,iv) = rkd(ii,jj,iv) + SUM(bq(:)*INTP(j,jp,:))
             END DO; END DO
             END DO; END DO

      End Subroutine triang_qpqp


!======================================================================
      Subroutine mrk_qqqq(k)
!======================================================================
!     Defines matrix of Rk integrals in the B-spline basis
!     by cell algorithm
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DBS_moments
      Use DBS_integrals

      Implicit none
      Integer, intent(in) :: k
      Integer :: i,j, ii,jj, iv,jv, ih,jh, ihp,jhp, ip,jp
      Integer, external :: Icheck_rka
      Real(8) :: c

! ... check the need of calculations

      if(itype == 'aaaa' .or. k < kra_min .or. k > kra_max) &
       Call alloc_DBS_integrals(ns,ks,0,k,4)

      if(itype == 'qqqq' .and. krk == k) Return
      if (associated(rkb)) nullify(rkb)
      rkb => rka(:,:,:,:,k,2)
      if(irka(k,2) == 1) then
       krk=k; itype = 'qqqq'
       Return
      end if

! ... compute the spline moments:

      Call moments_qq(  k   ,kk,nv,rkd1)
      Call moments_qq(-(k+1),kk,nv,rkd2)
      Call diag_qqqq(k)

! ... generate the rkb array

      rkb=0.d0

      DO jv=1,nv;    jj = 0
      DO jh = 1,ksq; j  = jv  + jh - 1
      DO jhp=jh,ksq; jp = jhp - jh + 1
                     jj = jj  + 1

      DO iv=1,nv;    ii = 0
      DO ih=  1,ksq; i  = iv  + ih - 1
      DO ihp=ih,ksq; ip = ihp - ih + 1
                     ii = ii  + 1

          if     ( iv < jv ) then;   c = rkd1(ii,iv)*rkd2(jj,jv)
          else if( iv > jv ) then;   c = rkd1(jj,jv)*rkd2(ii,iv)
          else;                      c = rkd(ii,jj,iv)
          end if

          rkb(i,j,ip,jp) = rkb(i,j,ip,jp) +  c

      END DO;  END DO;  END DO
      END DO;  END DO;  END DO

      krk=k; itype = 'qqqq'; irka(k,2) = 1

      END Subroutine mrk_qqqq


!======================================================================
      Subroutine diag_qqqq(k)
!======================================================================
!     Controls the scaling propeties for diagonal B-spline Rk-interals
!     (not implimented yet)
!----------------------------------------------------------------------
      Use DBS_grid, only: nv

      Implicit none
      Integer :: k,iv

      Do iv=1,nv; Call triang_qqqq(k,iv); End do

      END Subroutine diag_qqqq


!======================================================================
    Subroutine triang_qqqq (k,iv)
!======================================================================
!   Returns the two-dimensional array of B-spline integrals
!               <Q_i Q_j|r^k/r^(k+1)|Q_i' Q_j'>
!   over the given triangle diagonal cell
!
!   On entry   iv  -  index of the diagonal cell
!   --------
!
!   On exit    rkd(.,.,iv) - arrays of Rk B-spline integrals for given
!   --------                 interval iv in the reduced-dimension mode
!
!   Calls:   gauleg, zbsplvd
!----------------------------------------------------------------------
    Use DBS_grid
    Use DBS_gauss
    Use DBS_moments

    Implicit none
    Integer, intent(in) :: k,iv
    Integer :: i,j, ip,jp, ii,jj, m, left,  ik
    Real(8) :: xbase
    Real(8) :: x(ks),w(ks),bi(ks)
    Real(8) :: bspTmp(ks,ksq)
    Real(8) :: INTQ(ksq,ksq,ks)
    Real(8) :: a(ksq*(ksq+1)/2,ksq*(ksq+1)/2)

    left=iv+ksq-1;  xbase=tq(left)

!.. setup the gaussian points

    Call gauleg(0.d0,1.d0,x,w,ks)

    DO m=1,ks

! .. the absolute coordinate at the new gaussian point

      gx(:) = (gr(iv,m)-xbase)*x(:) + xbase

! .. the bspline values at the new gaussian points

      Do i=1,ks
       Call zbsplvd(nsq,ksq,nv, tq,left, 1,gx(i),1,dbiq)
       bspTmp(i,1:ksq)= dbiq(1,1:ksq,1)
      End do

! .. and the corresponding gaussian weights

      gw(:) = (gr(iv,m)-xbase)*w(:)

      IF(k>1) then;            gx(:) = gw(:)*gx(:)**k
      else IF(k==1) then;      gx(:) = gw(:)*gx(:)
      else IF(k==0) then;      gx(:) = gw(:)
      end if

!            / r(iv,m)                             k
! .. INT =  |      bsp(iv,:,j)(r) bsp(iv,:,jp)(r) r  dr
!           / r_iv

      Do j=1,ksq
       gw(:) = gx(:)*bspTmp(:,j)
       Do jp=j,ksq; INTQ(j,jp,m)= SUM(gw(:)*bspTmp(:,jp)); End do
      End do

    END DO   !  over m

! .. second integration

    IF(k/=0) then;   gx(:) = grw(iv,:)*grm(iv,:)**(k+1)
    else;            gx(:) = grw(iv,:)*grm(iv,:)
    end if

    ii = 0;  DO i=1,ksq;  DO ip=i,ksq;  ii = ii+1

             bi(:) = qbsp(iv,:,i)*qbsp(iv,:,ip)*gx(:)

    jj = 0;  DO j=1,ksq;  DO jp=j,ksq;  jj = jj+1

             a(ii,jj) =  SUM(bi(:)*INTQ(j,jp,:))

             END DO; END DO
             END DO; END DO

    ik = ksq*(ksq+1)/2;  rkd(1:ik,1:ik,iv) = a + TRANSPOSE(a)

    End Subroutine triang_qqqq


!======================================================================
      Subroutine msk_ppqq(k)
!======================================================================
!     Defines matrix of Sk integrals in the B-spline basis
!     by cell algorithm for PPQQ case.
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DBS_moments
      Use DBS_sk_integrals

      Implicit none
      Integer, intent(in) :: k
      Integer :: i,j, ii,jj, iv,jv, ih,jh, ihp,jhp, ip,jp
      Real(8) :: c

! ... check the need of calculations

      if(itype == 'aaaa' .or. k < kra_min .or. k > kra_max) &
       Call alloc_DBS_sk_integrals(ns,ks,0,k,2)

      if(itype == 'ppqq' .and. krk == k) Return
      if (associated(rkb)) nullify(rkb)
      rkb => rka(:,:,:,:,k,1)
      if(irka(k,1) == 1) then
       krk=k; itype = 'ppqq'
       Return
      end if

! ... compute the spline moments:

      Call moments_pq(  k   ,kk,nv,rkd1)
      Call moments_pq(-(k+1),kk,nv,rkd2)
      Call diag_ppqq(k)

! ... generate the rkb array

      rkb=0.d0

      DO jv = 1,nv;   jj = 0
      DO jh = 1,ksp;  j  = jv  + jh - 1
      DO jhp= 1,ksq;  jp = jhp - jh + ks
                      jj = jj + 1

      DO iv = 1,jv;   ii = 0
      DO ih = 1,ksp;  i  = iv  + ih - 1
      DO ihp= 1,ksq;  ip = ihp - ih + ks
                      ii = ii + 1

        IF( iv == jv ) THEN
          c =  rkd(ii,jj,iv)
        ELSE
          c =  rkd1(ii,iv)*rkd2(jj,jv)
        END IF

        rkb(i,j,ip,jp) = rkb(i,j,ip,jp) + c

      END DO;  END DO;  END DO
      END DO;  END DO;  END DO

      krk=k; itype = 'ppqq';  irka(k,1)=1

      End Subroutine msk_ppqq


!======================================================================
      Subroutine diag_ppqq(k)
!======================================================================
!     Controls the scaling propeties for diagonal B-spline Sk-interals
!     (not implimented)
!----------------------------------------------------------------------
      Use DBS_grid, only: nv

      Implicit none
      Integer :: k,iv

      Do iv=1,nv; Call triang_ppqq(k,iv); End do

      End Subroutine diag_ppqq


!======================================================================
    Subroutine triang_ppqq (k,iv)
!======================================================================
!   Returns the two-dimensional array of B-spline integrals
!               <B_i B_j|r^k/r^(k+1)|B_i' B_j'>
!   over the given triangle diagonal cell
!
!   On entry   iv  -  index of the diagonal cell
!   --------
!
!   On exit    rkd(.,.,iv) - arrays of Rk B-spline integrals for given
!   --------                 interval iv in the reduced-dimension mode
!
!   Calls:   gauleg, zbsplvd
!----------------------------------------------------------------------
    Use DBS_grid
    Use DBS_gauss
    Use DBS_moments

    Implicit none
    Integer, intent(in) :: k,iv
    Integer :: i,j, ip,jp, ii,jj, m, left
    REAL(8) :: xbase
    REAL(8) :: x(ks),w(ks),bi(ks)
    REAL(8) :: bspTmp(ks,ksp)
    REAL(8) :: bspTmq(ks,ksq)
    REAL(8) :: INTPQ(ksp,ksq,ks)

    left=iv+ks-1;  xbase=t(left)

!.. setup the gaussian points

    Call gauleg(0.d0,1.d0,x,w,ks)

    DO m=1,ks

! .. the absolute coordinate at the new gaussian point

      gx(:) = (gr(iv,m)-xbase)*x(:) + xbase

! .. the bspline values at the new gaussian points

      DO i=1,ks
       Call zbsplvd(nsp,ksp,nv, tp,iv+ksp-1, 1,gx(i),1,dbip)
       bspTmp(i,1:ksp)= dbip(1,1:ksp,1)
      END DO

      DO i=1,ks
       Call zbsplvd(nsq,ksq,nv, tq,iv+ksq-1, 1,gx(i),1,dbiq)
       bspTmq(i,1:ksq)= dbiq(1,1:ksq,1)
      END DO

! .. and the corresponding gaussian weights

      gw(:) = (gr(iv,m)-xbase)*w(:)

      gx(:) = gw(:)*gx(:)**k

!            / r(iv,m)                               k
! .. INT =  |      pbsp(iv,:,j)(r) qbsp(iv,:,jp)(r) r  dr
!           / r_iv

      Do j=1,ksp
       gw(:) = gx(:)*bspTmp(:,j)
       Do jp=1,ksq; INTPQ(j,jp,m)= SUM(gw(:)*bspTmq(:,jp));  END DO
      End do

    END DO     !  over m

! .. second integration

             gx(:) = grw(iv,:)*grm(iv,:)**(k+1)

             rkd(:,:,iv) = 0.d0

    ii = 0;  DO i=1,ksp;  DO ip=1,ksq;  ii = ii+1

             bi(:) = pbsp(iv,:,i)*qbsp(iv,:,ip)*gx(:)

    jj = 0;  DO j=1,ksp;  DO jp=1,ksq;  jj = jj+1

             rkd(jj,ii,iv) =  SUM(bi(:)*INTPQ(j,jp,:))

             END DO; END DO
             END DO; END DO

    End Subroutine triang_ppqq


!======================================================================
      Subroutine msk_pqqp(k)
!======================================================================
!     Defines matrix of Rk integrals in the B-spline basis
!     by cell algorithm
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DBS_moments
      Use DBS_sk_integrals

      Implicit none
      Integer, intent(in) :: k
      Integer :: i,j, ii,jj, iv,jv, ih,jh, ihp,jhp, ip,jp
      Real(8) :: c

! ... check the need of calculations

      if(itype == 'aaaa' .or. k < kra_min .or. k > kra_max) &
       Call alloc_DBS_sk_integrals(ns,ks,0,k,2)

      if(itype == 'pqqp' .and. krk == k) Return
      if (associated(rkb)) nullify(rkb)
      rkb => rka(:,:,:,:,k,2)
      if(irka(k,2) == 1) then
       krk=k; itype = 'pqqp'
       Return
      end if

! ... compute the spline moments:

      Call moments_pq(  k   ,kk,nv,rkd1)
      Call moments_qp(-(k+1),kk,nv,rkd2)
      Call diag_pqqp(k)

! ... generate the rkb array

      rkb=0.d0

      DO jv = 1,nv;   jj = 0
      DO jh = 1,ksq;  j  = jv  + jh - 1
      DO jhp= 1,ksp;  jp = jhp - jh + ks
                      jj = jj + 1

      DO iv = 1,jv;   ii = 0
      DO ih = 1,ksp;  i  = iv  + ih - 1
      DO ihp= 1,ksq;  ip = ihp - ih + ks
                      ii = ii + 1


        IF( iv == jv ) THEN
          c =  rkd(ii,jj,iv)
        ELSE
          c =  rkd1(ii,iv)*rkd2(jj,jv)
        END IF

        rkb(i,j,ip,jp) = rkb(i,j,ip,jp) + c

      END DO;  END DO;  END DO
      END DO;  END DO;  END DO

      krk=k; itype = 'pqqp';  irka(k,2)=1

      End Subroutine msk_pqqp


!======================================================================
      Subroutine diag_pqqp(k)
!======================================================================
!     Controls the scaling propeties for diagonal B-spline Rk-interals
!----------------------------------------------------------------------
      Use DBS_grid, only: nv

      Implicit none
      Integer :: k,iv

      Do iv=1,nv; Call triang_pqqp(k,iv); End do

      END Subroutine diag_pqqp


!======================================================================
    Subroutine triang_pqqp (k,iv)
!======================================================================
!   Returns the two-dimensional array of B-spline integrals
!               <P_i Q_j|r^k/r^(k+1)|P_i' Q_j'>
!   over the given triangle diagonal cell
!
!   On entry   iv  -  index of the diagonal cell
!   --------
!
!   On exit    rkd(.,.,iv) - arrays of Rk B-spline integrals for given
!   --------                 interval iv in the reduced-dimension mode
!
!   Calls:   gauleg, zbsplvd
!----------------------------------------------------------------------
    Use DBS_grid
    Use DBS_gauss
    Use DBS_moments

    Implicit none
    Integer, intent(in) :: k,iv
    Integer :: i,j, ip,jp, ii,jj, m, left
    Real(8) :: xbase
    Real(8) :: x(ks),w(ks),bi(ks)
    Real(8) :: bspTmp(ks,ksp)
    Real(8) :: bspTmq(ks,ksq)
    Real(8) :: INTpq(ksp,ksq,ks)

    left=iv+ks-1;  xbase=t(left)

!.. setup the gaussian points

    Call gauleg(0.d0,1.d0,x,w,ks)

    DO m=1,ks

! .. the absolute coordinate at the new gaussian point

      gx(:) = (gr(iv,m)-xbase)*x(:) + xbase

! .. the bspline values at the new gaussian points

      DO i=1,ks
       Call zbsplvd(nsp,ksp,nv, tp,iv+ksp-1, 1,gx(i),1,dbip)
       bspTmp(i,1:ksp)= dbip(1,1:ksp,1)
      END DO

      DO i=1,ks
       Call zbsplvd(nsq,ksq,nv, tq,iv+ksq-1, 1,gx(i),1,dbiq)
       bspTmq(i,1:ksq)= dbiq(1,1:ksq,1)
      END DO

! .. and the corresponding gaussian weights

      gw(:) = (gr(iv,m)-xbase)*w(:)
      gx(:) = gw(:)*gx(:)**k

!            / r(iv,m)                               k
! .. INT =  |      pbsp(iv,:,j)(r) qbsp(iv,:,jp)(r) r  dr
!           / r_iv

      Do j=1,ksp
       gw(:) = gx(:)*bspTmp(:,j)
       Do jp=1,ksq; INTpq(j,jp,m)= SUM(gw(:)*bspTmq(:,jp));  END DO
      End do

    END DO !  over m

! .. second integration

    gx(:) = grw(iv,:)*grm(iv,:)**(k+1)

             rkd(:,:,iv) = 0.d0

    ii = 0;  DO i=1,ksq;  DO ip=1,ksp;  ii = ii+1

             bi(:) = qbsp(iv,:,i)*pbsp(iv,:,ip)*gx(:)

    jj = 0;  DO j=1,ksp;  DO jp=1,ksq;  jj = jj+1

             rkd(jj,ii,iv) =  SUM(bi(:)*INTpq(j,jp,:))

             END DO; END DO
             END DO; END DO

    End Subroutine triang_pqqp


!======================================================================
      Real(8) Function quadr (fi,fj,m)
!======================================================================
!     Evaluates   <fi | r^m | fj>     with respect to r
!     where fi, fj  - two-component Dirac functions in B-spline basis
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss

      Implicit none
      Integer, intent(in) :: m
      Real(8), intent(in) :: fi(ns,2),fj(ns,2)
      Integer :: i,j, iv,ith,jth
      Real(8) :: pi,pj,qi,qj

      quadr = 0.d0

      Do iv=1,nv
       gx(:) = gr(iv,:)**m * grw(iv,:)
       Do ith=1,ksp;  i=iv+ith-1; pi=fi(i,1)
        gw(:) = gx(:)*pbsp(iv,:,ith)
        Do jth=1,ksp;  j=iv+jth-1; pj=fj(j,1)
         quadr = quadr + SUM(pbsp(iv,:,jth)*gw(:))*pi*pj
        End do
       End do
      End do

      Do iv=1,nv
       gx(:) = gr(iv,:)**m * grw(iv,:)
       Do ith=1,ksq;  i=iv+ith-1; qi=fi(i,2)
        gw(:) = gx(:)*qbsp(iv,:,ith)
        Do jth=1,ksq;  j=iv+jth-1; qj=fj(j,2)
         quadr = quadr + SUM(qbsp(iv,:,jth)*gw(:))*qi*qj
        End do
       End do
      End do

      End Function quadr



!======================================================================
      Real(8) Function rk (f1,f2,f3,f4,k)
!======================================================================
!     Returns  rk - integral, base on the assembling two-electron
!     B-spline integrals (see module DBS_integral)
!----------------------------------------------------------------------
      Use DBS_grid,         only: ns,ks

      Implicit none
      Integer, intent(in) :: k
      Real(8), intent(in) :: f1(ns,2),f2(ns,2),f3(ns,2),f4(ns,2)
      Real(8) :: dens1(ns,ks),dens2(ns,ks),dens3(ns,ks),dens4(ns,ks), &
                 conv(ns,ks), S
      Real(8), external :: SUM_AmB

      rk = 0.d0

      Call density (ns,ks,dens1,f1(1,1),f3(1,1),'s')
      Call density (ns,ks,dens2,f2(1,1),f4(1,1),'s')
      Call density (ns,ks,dens3,f1(1,2),f3(1,2),'s')
      Call density (ns,ks,dens4,f2(1,2),f4(1,2),'s')

      Call mrk_pppp(k)
      Call convol  (ns,ks,conv,dens1,2,'s','s')
      S = SUM_AmB(ns,ks,conv,dens2,'s')
      rk = rk + S

      Call mrk_qqqq(k)
      Call convol  (ns,ks,conv,dens3,2,'s','s')
      S = SUM_AmB(ns,ks,conv,dens4,'s')
      rk = rk + S

      Call mrk_qpqp(k)
      Call convol  (ns,ks,conv,dens3,2,'s','s')
      S = SUM_AmB(ns,ks,conv,dens2,'s')
      rk = rk + S

      Call mrk_pqpq(k)
      Call convol  (ns,ks,conv,dens1,2,'s','s')
      S = SUM_AmB(ns,ks,conv,dens4,'s')
      rk = rk + S

      End Function rk



!======================================================================
    Real(8) Function SUM_AmB(ns,ks,a,b,sym)
!======================================================================
!   Returns  SUM (a * b)  for banded and full matrixes
!   Possible array representation for matrixes a and b:
!   sym = 's' - symmetrical banded matrix, upper part
!   sym = 'l' - symmetrical banded matrix, lower part
!   sym = 'n' - non-symmetrical banded matrix
!   sym = 'x' - full matrix
!----------------------------------------------------------------------
    Implicit none
    Integer, intent(in) :: ns,ks
    Real(8), intent(in) :: a(ns,*),b(ns,*)
    Character(*), intent(in) :: sym
    Integer :: i,j, imin,imax
    Real(8) :: x

    x = 0.d0

    if(sym.eq.'x') then

     Do i = 1,ns;  Do j = 1,ns
      x = x + a(i,j)*b(i,j)
     End do; End do

    elseif(sym.eq.'s') then

     Do j = 1,ks; Do i = 1,ns-j+1
       x = x + a(i,j)*b(i,j)
     End do; End do

    elseif(sym.eq.'l') then

     Do j=1,ks;  Do i=ks+1-j,ns
      x = x + a(i,j)*b(i,j)
     End do; End do

    elseif(sym.eq.'n') then

     Do j = 1,ks+ks-1
      imin=max( 1, 1+ks-j)
      imax=min(ns,ns+ks-j)
     Do i = imin,imax
       x = x + a(i,j)*b(i,j)
     End do; End do

    else

     Stop ' SUM_AmB:  unknown symmetry '

    end if

    SUM_AmB = x

    End Function SUM_AmB


!======================================================================
      Subroutine TVM(k,p,q,TA,VA,MA)
!======================================================================
!     Calculate the expectation values for
!     TA - kinatic energy
!     VA - potential energy
!     MA - mass energy
!     for orbital (p,q) with kappa value "k"
!----------------------------------------------------------------------
      USE zconst
      USE DBS_nuclear, z => atomic_number
      USE DBS_grid
      USE DBS_gauss

      Implicit none
      Integer, intent(in) :: k
      Real(8), intent(in) :: p(ns),q(ns)
      Real(8), intent(out) :: TA,VA,MA
      Real(8) :: sa,sb,vv(ns),a(ns,ns)

! ... potential average:

      vv=MATMUL(fpb_nucl,p)
      VA = dot_product(p,vv)
      vv=MATMUL(fqb_nucl,q)
      VA = VA + dot_product(q,vv)

! ... mass average:

      sb = -2.d0 * c_au * c_au
      vv=MATMUL(fqbs,q)
      MA = sb*dot_product(q,vv)

! ... kinetic average:

! ... matrix of (2:1) = c <B | D+ | B>

      sa = c_au; sb = dble(k)*c_au
      a(1:ns,1:ns) = sa*fqpbsd + sb*fqpbs
      vv=MATMUL(a,p)
      TA = dot_product(q,vv)

! ... matrix of (1:2) = c <B | D- | B>

      sa = -c_au; sb = dble(k)*c_au
      a(1:ns,1:ns) = sa*fpqbsd + sb*fpqbs
      vv=MATMUL(a,q)
      TA = TA + dot_product(p,vv)

      End Subroutine TVM



!======================================================================
      Subroutine UPDATE_HS (ms,hm,ii,jj,ns,ks,d,sym)
!======================================================================
!     update 'small' channel block (ns,ns) in "big" channel block (ms,ms)
!     (ii,jj) = (1,1),(1,2),(2,1),(2,2) - define the "small" block
!     Possible array representation for input d-matrix:
!     sym = 's'  -->  symmetric banded upper-column storage mode
!     sym = 'n'  -->  non-symmetric band matrix
!     sym = 'x'  -->  non-symmetric full matrix
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: ii,jj,ms,ns,ks
      Real(8), intent(in) :: d(ns,*)
      Real(8), intent(inout) :: hm(ms,ms)
      Character(*), intent(in) :: sym
      Integer :: i,j, jp, imin,imax, i1,i2,j1,j2
      Real(8) :: y(ns,ns), yy(ns,ns)

      Select case(sym)

      Case('s')

       y = 0.d0
       do jp=1,ks;  do i=1,ns-jp+1;  j=i+jp-1
        y(i,j) = d(i,jp)
        y(j,i) = d(i,jp)
       end do; end do

      Case('n')

       y = 0.d0
       do jp = 1,ks+ks-1
         imin=max0( 1, 1 + ks-jp)
         imax=min0(ns,ns + ks-jp)
         do i = imin,imax;  j=i+jp-ks
           y(i,j) = d(i,jp)
         end do
       end do

      Case('x')

       y(1:ns,1:ns) = d(1:ns,1:ns)

      End Select

      Select Case (10*ii+jj)
       Case(11);  i1=1; i2=ns; j1=1; j2=ns
       Case(12);  i1=1; i2=ns; j1=ns+1;j2=ms
       Case(21);  i1=ns+1;i2=ms; j1=1; j2=ns
       Case(22);  i1=ns+1;i2=ms; j1=ns+1;j2=ms
       Case default; Stop 'unknown case in update_hs'
      End Select

      y = y/2.d0; yy = TRANSPOSE(y)
      hm(i1:i2,j1:j2) = hm(i1:i2,j1:j2) + y
      hm(j1:j2,i1:i2) = hm(j1:j2,i1:i2) + yy

      End Subroutine UPDATE_HS


!=======================================================================
   SUBROUTINE zbsplvd(ns,ks,nv, t, kg, ni, x, nderiv, dbiatx)
!=======================================================================
!
!  This routine calculates the values of the B-splines and their deriva-
!  tives, of order up to nderiv, that do not vanish at x(i), i=1,..ni
!  There are ks such B-splines at each point.
!
!  This routine is a vector version of bsplvd written by C. de Boor,
!  ``A Practical Guide to Splines".
!
!  subroutine contained: vbsplvb
!
!-----------------------------------------------------------------------
!  on entry
!  --------
!  t     the knot array, of length nt >=nv+2ks-1.  It is assumed
!        that t(i) < t(i+1) for each interval containing an x(i)
!        Division by zero will result otherwise (in vbsplvb).
!
!  kg    gives the beginning interval from which the B-splines are
!        evaluated at the Gaussian points.
!
!  ni    the number of intervals in which B-splines are to be evaluated
!        at all Gaussian points, its uplimit is nt.
!
!  x     the point array at which these values are sought,
!        one per interval, of length ni.
!
!  nderiv   an integer indicating that values of B-splines and their
!        derivatives up to but not including the  nderiv-th  are asked
!        for.
!
!  working area
!  ------------
!  w31   a three dimensional array, w31(i,j,m) (j=1,..,ks m=1,..,ks) con-
!        tains B-coeff.s of the derivatives of a certain order of the
!        ks B-splines of interest at point x(i)
!
!  w1,w2      one dimensional arrays
!
!  on return
!  ---------
!  dbiatx     a three dimensional array. its entry (i,j,m) contains
!        value of  (m-1)st  derivative of  (l-ks+j)-th B-spline of
!        order ks at point x(i) for knot sequence t, i=1..ni,
!        j=m..ks; m=1..nderiv;and l=kg..kg+ni-1
!
!  method
!  ------
!  values at x(i) of all the relevant B-splines of order ks,ks-1,...,
!  ks+1-nderiv  are generated via vbsplvb and stored temporarily
!  in dbiatx. then, the B-coeffs of the required derivatives of the
!  B-splines of interest are generated by differencing, each from the
!  preceding one of lower order, and combined with the values of B-
!  splines of corresponding order in dbiatx to produce the desired
!  values.
!----------------------------------------------------------------------
    Implicit none
    Integer, intent(in) :: ns,ks,nv, kg, ni, nderiv
    Real(8), intent(in) :: t(ns+ks)
    Real(8), intent(inout):: dbiatx(nv,ks,ks)
    Real(8), intent(in):: x(ni)
    ! local variables
    Real(8) :: fkpimm
    Real(8) :: w1(ni),w2(ni)
    Real(8) :: w31(ni,ks,ks)
    Real(8) :: deltar(ni,ks), deltal(ni,ks)
    Integer :: i, j, n, m, mhigh, kp1, jhigh, ideriv, &
               ldummy, kpimm, jlow, jpimid, il

    mhigh = MAX(MIN(nderiv,ks),1)   !mhigh is usually equal to nderiv.
    kp1 = ks+1
    jhigh = kp1-mhigh
    Call zbsplvb(kg,ni,x,jhigh,1,dbiatx)
    IF(mhigh == 1) Return

    ! ..the first row of dbiatx always contains the B-spline values
    ! ..for the current order. these are stored in row ks+1-current
    ! ..order before vbsplvb is called to put values for the next
    ! ..higher order on top of it. Vbsplvb only uses the first two dimensions

    ideriv = mhigh
    Do m = 2, mhigh
      jpimid = 1
      Do j = ideriv,ks
       dbiatx(1:ni,j,ideriv) = dbiatx(1:ni,jpimid,1)
       jpimid = jpimid+1
      End do
      ideriv = ideriv-1
      jhigh = kp1-ideriv
      Call zbsplvb(kg,ni,x,jhigh,2,dbiatx)
    End do

    ! at this point,  b(.,n-ks+i, ks+1-j)(x) is in dbiatx(.,i,j) for
    ! n=kg..kg+ni-1,i=j..ks,j=1..mhigh('='nderiv).in particular,the
    ! first row of  dbiatx  is already in final form. to obtain cor-
    ! responding derivatives of B-splines in subsequent rows, gene-
    ! rate their B-repr. by differencing, then evaluate at x(.).

    jlow = 1
    Do i = 1,ks
      w31(1:ni,jlow:ks,i) = 0.d0
      jlow = i
      w31(1:ni,i,i) = 1.d0
    End do

    ! at this point, w31(.,.,j) contains the B-coeffs for the j-th of the
    ! ks B-splines of interest here.

    Do m = 2,mhigh
      kpimm = kp1-m
      fkpimm = kpimm
      i = ks
      il = 0

      ! for j=1,...,ks, construct B-coeffs of  (m-1)st  derivative of
      ! B-splines from those for preceding derivative by differencing
      ! and store again in  w31(.,.,j). the fact that w31(i,j) = 0  for
      ! i < j is used.

      Do ldummy = 1, kpimm
       Do n = kg-il,ni+kg-il-1
        w1(n-kg+il+1) = fkpimm/(t(n+kpimm)-t(n))
       End do

        ! the assumption that t(n) < t(n+1) makes denominator
        ! in w1(1..ni) nonzero.

        Do j = 1,i
         w31(1:ni,i,j) = (w31(1:ni,i,j)-w31(1:ni,i-1,j))*w1(1:ni)
        End do
        il = il+1
        i = i-1
      End do

      ! for i=1,...,ks, combine B-coeffs a(.,.,i) with B-spline values
      ! stored in dbiatx(.,.,m) to get value of (m-1)st  derivative of
      ! i-th B-spline (of interest here) at x(.), and store in
      ! dbiatx(.,i,m). storage of this value over the value of a B-spline
      ! of order m there is safe since the remaining B-spline derivat-
      ! ive of the same order do not use this value due to the fact
      ! that  a(.,j,i) = 0  for j .lt. i .

      Do i = 1,ks
        w2(1:ni) = 0.d0
        jlow = MAX(i,m)
        Do j = jlow,ks
          w2(1:ni) = w2(1:ni) + w31(1:ni,j,i)*dbiatx(1:ni,j,m)
        End do
        dbiatx(1:ni,i,m) = w2(1:ni)
      End do

    End do

CONTAINS


!===================================================================
   Subroutine zbsplvb(kg, ni, x, jhigh, index, biatx)
!=====================================================================
!  This routine calculates the values of all possibly nonzero B-splines
!  at x(i) (i=1,..ni) of order
!               jout=max(jhigh,(j+1)*(index-1))
!  with knot sequence t.
!
!  This routine is a vector version of bsplvb written by C. de Boor,
!  "A Practical Guide to Splines", Chapter X, page 135
!
!  on entry
!  --------
!  t    -  knot sequence, of length nt=ns+ks, assumed to be nonde-
!          creasing, that is t(i) <= t(i+1)
!
!  jhigh-  choose jhigh=1 to get the B-spline values directly
!            by calling vbsplvb.
!
!  kg   -  gives the beginning interval from which the B-splines
!           are to be evaluated at Gaussin points.
!
!  ni   -  the number of intervals in which B-splines are
!            evaluated at all Gaussian points, its uplimit is nv.
!
!  x    -  the points at which the B-splines are to be evaluated,
!            its length is ni,
!
!  index-  integers which determine the order  jout = max(jhigh,
!            (j+1)*(index-1))  of the B-splines whose values at x(i)
!            are to be returned.  index is used to avoid recalcula-
!            tions when several columns of the triangular array of
!            B-spline values are needed (e.g., in vbsplvd ).
!            More precisely,
!                     if index = 1 ,
!            the calculation starts from scratch and the entire
!            triangular array of B-spline values of orders
!            1,2,...,jhigh is generated, order by order ,
!            i.e., column by column .
!                     if  index = 2 ,
!            only the B-spline values of order  j+1, j+2, ..., jout
!            are generated, the assumption being that  biatx,j,
!            deltal,deltar are, on  entry, as they were on exit at the
!            previous call. In particular, if  jhigh = 0, then
!            jout = j+1, i.e., just the next column of B-spline
!            values is generated.
!
!  working area
!  ------------
!  deltal, deltar: two dimensional arrays
!  term, saved:    one dimensional arrays.
!
!  on return
!  ---------
!  biatx.....two dimensional array, with biatx(j-k+1,i)(j=k..ni)
!        containing the value at x(j-k+1) of the polynomial of order
!        jout which agrees with the B-spline b(j-jout+1,jout,t) on
!        the interval (t(j),t(j+1)).
!
!  method
!  ------
!  The recurrence relation
!
!                       x - t(i)              t(i+j+1) - x
!     b(i,j+1)(x)  =  -----------b(i,j)(x) + ---------------b(i+1,j)(x)
!                     t(i+j)-t(i)            t(i+j+1)-t(i+1)
!
!  is used (repeatedly) to generate the (j+1)-vector  b(l-j,j+1)(x),
!  ...,b(l,j+1)(x)  from the j-vector  b(l-j+1,j)(x),...,
!  b(l,j)(x), storing the new values in  biatx  over the old. the
!  facts that
!            b(i,1) = 1  if  t(i) <= x < t(i+1)
!  and that
!            b(i,j)(x) = 0  unless  t(i) <= x < t(i+j)
!  are used. the particular organization of the calculations follows al-
!  gorithm (8) in chapter x of the text.
!-----------------------------------------------------------------------
    Implicit none
    Integer, intent(in):: kg, ni, jhigh, index
    Real(8), intent(in):: x(ni)
    Real(8), intent(inout):: biatx(nv,ks)
    ! .. Local variables
    Integer :: i, jp1, m
    Integer, save :: j=1
    Real(8) :: term(ni), saved(ni)

    if(index == 1) then
      j=1
      biatx(1:ni,1)=1.d0
      IF (j >= jhigh)  Return
    end if

    Do
      jp1=j+1
      saved(1:ni)=0.d0
      Do i=1,ni
        deltar(i,j) = t(i+kg-1+j) - x(i)
        deltal(i,j) = x(i) - t(i+kg-j)
      End do

      Do m=1,j
        Do i=1,ni
          term(i)=biatx(i,m)/(deltar(i,m)+deltal(i,jp1-m))
          biatx(i,m)=saved(i)+deltar(i,m)*term(i)
          saved(i)=deltal(i,jp1-m)*term(i)
        End do
      End do

      biatx(1:ni,jp1)=saved(1:ni)
      j=jp1
      if (j >= jhigh) Exit
    End do
    End Subroutine zbsplvb

    End Subroutine zbsplvd


!======================================================================
    Subroutine ZINTYM (nv,ks,ks1,ks2,B1,B2,ygr,mm,ym)
!======================================================================
!   Computes the array  ym(i,j) =  <B1_i | y(r) | B2 _j>
!
!   B1, B2 - two bases, spesified in gausian points ks
!            which may be different from ks1 or ks2
!
!   ygr - array of values of a specific function  y(r) at the gaussian
!         points of each interval, weighted by the gaussian weights
!----------------------------------------------------------------------
    Implicit none
    Integer, intent(in)  :: nv,ks,ks1,ks2,mm
    Real(8), intent(in)  :: B1(nv+1,ks,ks1),B2(nv+1,ks,ks2),ygr(nv,ks)
    Real(8), intent(out) :: ym(mm,mm)
    Integer :: iv,i,j,ith, jth

    ym = 0.d0
    Do iv = 1,nv                      ! over intervals
     Do ith = 1,ks1; i = iv+ith-1     ! over B1 splines
     Do jth = 1,ks2; j = iv+jth-1     ! over B2 splines
        ym(i,j) = ym(i,j) + SUM(ygr(iv,:)*b1(iv,:,ith)*b2(iv,:,jth))
     End do
     End do
    End do

    End Subroutine ZINTYM




