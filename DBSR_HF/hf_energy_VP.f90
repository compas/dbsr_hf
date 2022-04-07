!======================================================================
      Subroutine VP_energy(VPE)
!======================================================================
! ... vacuum-polarization corrections
!----------------------------------------------------------------------
      Use dbsr_hf
      Use DF_orbitals
      Use DBS_nuclear
      Use DBS_grid
      Use DBS_gauss

      Implicit none
      Integer :: i
      Real(8) :: VPE
      Real(8) :: vacpol(nv,ks),fp_nucl(ns,ns),fq_nucl(ns,ns),v(ns)

      VPE = 0.d0; if(mode_VP.eq.0) Return

      vacpol = 0.d0; fp_nucl=0.d0; fq_nucl=0.d0

      Call Def_Znucl

      Call set_vacpol_2 (vacpol)
      Call set_vacpol_4 (vacpol)

      vacpol = vacpol * grw

      Call ZINTYm (nv,ks,ksp,ksp,pbsp,pbsp,vacpol,ns,fp_nucl)
      Call ZINTYm (nv,ks,ksq,ksq,qbsp,qbsp,vacpol,ns,fq_nucl)

      VPE = 0.d0
      Do i=1,nbf
       v = matmul(fp_nucl,p(1:ns,1,i))
       VPE = VPE +  dot_product(v(1:ns),p(1:ns,1,i)) * qsum(i)
       v = matmul(fq_nucl,p(1:ns,2,i))
       VPE = VPE +  dot_product(v(1:ns),p(1:ns,2,i)) * qsum(i)
      End do

      End Subroutine VP_energy


!======================================================================
      Subroutine set_vacpol_2(vacpol)
!======================================================================
! Sets up the second-order vacuum polarization potential using
! equations (1) and (4) of L W Fullerton and  G A Rinker, Jr,
! Phys Rev A  13 (1976) 1283-1287.
!
! This routine is similar to vac2 from grasp92 [RCI92] which was
! written by F A Parpia (see also RATIP)
!
! Calls: vacpol_Kn
!--------------------------------------------------------------------
      Use zconst
      Use DBS_grid,    only: nv,ks
      Use DBS_gauss,   only: gr,ygw,grw
      Use DBS_nuclear, only: nuclear, atomic_number, ZR_nucl

      Implicit none
      Integer :: i,m, ii,mm
      Real(8) :: epsi, factor, twocv, x, xi, xm, xk, xp
      Real(8) :: vacpol(nv,ks)
      Real(kind=dp), external :: vacpol_Kn

      epsi  = 1.0e-20_dp
      twocv = c_au + c_au
      vacpol = zero

! ... Potential for a point nucleus, equation (1):
!                   2 Z
!     V_2(r) = -  -------- K_1(2cr)
!                 3 c pi r
! ... (this is also the asymptotic form for a finite nucleus)

      factor = -(two * atomic_number)/(three * pi * c_au)
      Do i=1,nv; Do m=1,ks
       x = twocv * gr(i,m)
       vacpol(i,m) = factor/gr(i,m) * vacpol_Kn(x,1)
       if (abs(vacpol(i,m)) < epsi) factor=0.d0
       if(factor.eq.0.d0) Exit
      End do; if(factor.eq.0.d0) Exit; End do

      if (nuclear.eq.'point') Return

! ... Potential for a finite nucleus: equation (4)
!                    2
!     V_2(r) = - -------- INT[r',0,inf; r' ro(r') {K_0(2c|r-r'|)-K_0(2c|r+r'|} ]
!                3 c^2 r

      factor = -two / (three * c_au**2)
      Do i=1,nv; Do m=1,ks
       xk = twocv * gr(i,m)
       ! set up integrand
       ygw = zero
       Do ii=1,nv; Do mm=1,ks
        xi = twocv * gr(ii,mm)
        xm = abs(xk - xi)
        xp = abs(xk + xi)
        ygw(ii,mm) = (vacpol_Kn(xm,0)-vacpol_Kn(xp,0))*ZR_nucl(ii,mm)
       End do; End do
       x = SUM(ygw*gr) * factor / gr(i,m)
       if(abs(x) < epsi) x=0.d0
       if(abs((x-vacpol(i,m))/x) < 1.d-5) x=0.d0
       vacpol(i,m) = x
       if(x.eq.0.d0) Exit
      End do; if(x.eq.0.d0) Exit; End do

      End Subroutine set_vacpol_2


!======================================================================
      Subroutine set_vacpol_4(vacpol)
!======================================================================
! Sets up the fourth-order vacuum polarization potential using
! equations (11) and (12) of L Wayne Fullerton and  G A Rinker, Jr,
! Phys  Rev  A 13 (1976) 1283-1287.
! This routine is similar to vac4 from GRASP92 [RCI92] which was
! written by F A Parpia.
!
! Calls: vacpol_Kn
!
! Remark: this potential is added to existing one.
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DBS_nuclear

      Implicit none
      Integer       :: i,m, ii,mm
      Real(kind=dp) :: epsi, factor, twocv, x, xi, xm, xk, xp
      Real(kind=dp) :: vacpol(nv,ks)
      Real(kind=dp), external :: vacpol_Lk

      epsi  = 1.0e-20_dp
      twocv = c_au + c_au

! ... Potential for a point nucleus: equation (12)
!                     Z
!     V_4(r) = -  ----------  L_1(2cr)
!                 c^2 pi^2 r
! ... (this is also the asymptotic form for a finite nucleus)

      factor = -atomic_number / (pi*c_au)**2

      Do i=1,nv; Do m=1,ks
       x = twocv * gr(i,m)
       vacpol(i,m) = vacpol(i,m) + factor/gr(i,m) * vacpol_Lk(x,1)
       if (abs(vacpol(i,m)) < epsi) factor=0.d0
       if(factor.eq.0.d0) Exit
      End do; if(factor.eq.0.d0) Exit; End do

      if (nuclear.eq.'point') Return

! ... Potential for finite nucleus: equation (11)
!                    1
!     V_4(r) = - -------- INT[r',0,inf; r' V(r') {L_0(2c|r-r'|)-L_0(2c|r+r'|} ]
!                c^3 pi r

      if (nuclear.eq.'point') Return

      factor = -one / (pi*c_au**3)
      Do i=1,nv; Do m=1,ks
       xk = twocv * gr(i,m)
       ! set up integrand
       ygw = zero
       Do ii=1,nv; Do mm=1,ks
        xi = twocv * gr(ii,mm)
        xm = abs(xk - xi)
        xp = abs(xk + xi)
        ygw(ii,mm) = (vacpol_Lk(xm,0)-vacpol_Lk(xp,0))*ZR_nucl(ii,mm)
       End do; End do
       x = SUM(ygw*gr) * factor / gr(i,m)
       if(abs(x) < epsi) x = 0.d0
       if(abs((x-vacpol(i,m))/x) < 1.0e-3_dp) x = 0.d0
       vacpol(i,m) = vacpol(i,m) + x
       if(x.eq.0.d0) Exit
      End do; if(x.eq.0.d0) Exit; End do

      End subroutine set_vacpol_4


!========================================================================
   function vacpol_Kn(x,n)                        result(Kn)
   !--------------------------------------------------------------------
   ! Evaluates the K_N(X) functions using the analytic functions defined
   ! in tables 1 and 3 of Fullerton and Rinker, Phys  Rev  A 13 (1976)
   ! 1283-1287.
   ! This routine is taken from RATIP and similar to ... from g
   ! [RCI92] which was written by F A Parpia;
   ! it has been adapted in RATIP to the Fortran90/95 standard.
   !--------------------------------------------------------------------
      Use zconst

      Implicit none

      integer, intent(in)       :: n
      real(kind=dp), intent(in) :: x
      real(kind=dp)             :: Kn
      !
      real(kind=dp), dimension(10,4), parameter :: p = reshape(source =   &
        (/  8.8357293375e-1_dp, -2.8259817381e-1_dp, -5.8904879578e-1_dp, &
            1.2500133434e-1_dp, -3.2729913852e-2_dp,  8.2888574511e-3_dp, &
           -1.0327765800e-5_dp,  6.3643668900e-5_dp,  0.0_dp,             &
            0.0_dp,                                                       &
           -7.1740181754e-1_dp,  1.1780972274_dp,    -3.7499963087e-1_dp, &
            1.3089675530e-1_dp, -3.8258286439e-2_dp, -2.4297287300e-5_dp, &
           -3.5920148670e-4_dp, -1.7170090700e-5_dp,  0.0_dp,             &
            0.0_dp,                                                       &
            9.9999999987e-1_dp,  1.9770200000e-8_dp, -7.5000050190e-1_dp, &
            7.8540306316e-1_dp, -3.4988601655e-1_dp,  6.4596333000e-5_dp, &
           -9.8189080747e-3_dp,  8.6513145800e-5_dp, -2.3969236620e-4_dp, &
            0.0_dp,                                                       &
            6.0000000002_dp,    -6.4305200000e-8_dp,  2.1049413000e-6_dp, &
           -2.6711271500e-5_dp, -1.3705236152e-1_dp, -6.3476104090e-4_dp, &
           -7.8739801501e-2_dp, -1.9641740173e-3_dp, -3.4752369349e-3_dp, &
           -7.3145316220e-4_dp  /), shape=(/10,4/))
      !
      real(kind=dp), dimension(2,4), parameter :: b = reshape(source =    &
        (/ -3.19999594323e+2_dp,    2.53900995981_dp,                     &
           -6.40514843293e+1_dp,    7.11722714285e-1_dp,                  &
            5.19010136460e+3_dp,    8.28495496200e+1_dp,                  &
            3.18150793824e+2_dp,    4.33898867347e+1_dp  /), shape=(/2,4/))
      !
      real(kind=dp), dimension(3,4), parameter :: c = reshape(source =    &
        (/ -3.19999594333e+2_dp,    2.53901020662_dp,                     &
            0.0_dp,                 6.40514843287e+1_dp,                  &
           -7.11722686403e-1_dp,    8.04220774800e-4_dp,                  &
            2.76805406060e+4_dp,   -3.27039477790e+2_dp,                  &
            0.0_dp,                 8.48402116837e+2_dp,                  &
           -2.56939867765e+1_dp,    3.20844906346e-1_dp  /), shape=(/3,4/))
      !
      real(kind=dp), dimension(5,4), parameter :: d = reshape(source =    &
        (/  5.018065179_dp,         7.1518912620e+1_dp,                   &
            2.116209929e+2_dp,      3.1403274780e+1_dp,                   &
           -1.0_dp,                 2.1723864090e+2_dp,                   &
            1.643364528e+3_dp,      2.1222445120e+3_dp,                   &
           -4.512004044e+1_dp,      1.0_dp,                               &
            8.540770444_dp,         6.0762427660e+1_dp,                   &
            9.714630584e+1_dp,      3.1549735930e+1_dp,                   &
            1.0_dp,                 5.9243015865e-1_dp,                   &
            2.0596312871_dp,        3.7785190424_dp,                      &
            3.5614853214_dp,        1.0_dp               /), shape=(/5,4/))
      !
      real(kind=dp), dimension(5,4), parameter :: e = reshape(source =    &
        (/  2.669207401_dp,         5.172549669e+1_dp,                    &
            2.969809720e+2_dp,      5.364324164e+2_dp,                    &
            1.535335924e+2_dp,      1.155589983e+2_dp,                    &
            1.292191441e+3_dp,      3.831198012e+3_dp,                    &
            2.904410075e+3_dp,      0.0_dp,                               &
            4.543392478_dp,         3.514920169e+1_dp,                    &
            6.019668656e+1_dp,      8.468839579_dp,                       &
            0.0_dp,                 3.1511867816e-1_dp,                   &
            3.473245222e-1_dp,      3.8791936870e-2_dp,                   &
           -1.3059741497e-3_dp,     0.0_dp               /), shape=(/5,4/))
      !
      integer, dimension(4), parameter :: np = (/ 8, 8, 9, 10 /)
      !
      integer       :: i, k, nn
      real(kind=dp) :: bsum, csum, dsum, esum, sum, x2, xm, xn
      logical       :: use_stop = .TRUE.
      !
      if (x == zero) goto 11
      if (use_stop   .and.                               &
         (n < 0   .or.   n == 2   .or.   n == 4   .or.  n > 5)) then
         print *, "Attempt to calculate FUNK (X,N) for N other than "// &
                  "0, 1, 3 and 5."
         stop     "vacpol_Kn(): program stop A."
      end if
      !
      select case(n-3)
      case(:-1)
         k  = n+1
         xn = one
      case(0)
         k  = n
         xn = one / (x**2)
      case(1:)
         k  = n-1
         xn = one / (x**4)
      case default
         stop "vacpol_Kn(): program stop B."
      end select
      if (x > one) goto 9
      !
      ! Calculate function for x < = 1
      nn  = np(k)
      sum = zero
      do  i = 1,nn
         sum = sum + p(i,k) * xn
         xn  = xn * x
      end do
      !
      x2   = x * x
      bsum = b(1,k) + x2*(b(2,k) + x2)
      csum = c(1,k) + x2*(c(2,k) + x2*c(3,k))
      !
      select case(k)
      case(1)
         bsum = bsum * x
      case(2,4)
      case(3)
         bsum = bsum * x2
      case default
         stop "vacpol_Kn(): program stop C."
      end select
      sum = sum+bsum*log (x)/csum
      !
      Kn = sum
      return
      !
      ! Calculate function for x > 1
    9 xn   = one
      dsum = zero
      esum = zero
      do  i = 1,5
         dsum = dsum + d(i,k) * xn
         esum = esum + e(i,k) * xn
         xn   = xn / x
      end do
      xm  = -x
      sum = dsum * exp(xm) / (esum*sqrt(x**3))
      !
      Kn  = sum
      return
      !
   11 if (use_stop   .and.   n /= 0) then
         print *, "Attempt to calculate Kn (0,N) for N > 0."
         stop     "vacpol_Kn(): program stop D."
      end if
      !
      Kn  = p(1,1)
      !
   end function vacpol_Kn



!========================================================================
   function vacpol_Lk(x,k)                        result(Lk)
   !--------------------------------------------------------------------
   ! Evaluates the LK(X) functions using the analytic functions defined
   ! in table 5  and equations (20) and  (21) of Fullerton and Rinker,
   ! Phys  Rev  A 13 (1976) 1283-1287.
   ! This routine is taken from RATIP and similar to ... from g
   ! [RCI92] which was written by F A Parpia;
   ! it has been adapted in RATIP to the Fortran90/95 standard.
   !--------------------------------------------------------------------
      Use zconst

      Implicit none

      integer, intent(in)       :: k
      real(kind=dp), intent(in) :: x
      real(kind=dp)             :: Lk
      !
      real(kind=dp), dimension(6,2), parameter :: f = reshape(source =      &
        (/ 2.008188_dp,    -2.397605_dp,    1.046471_dp,    -3.670660e-1_dp,&
           6.374000e-2_dp, -3.705800e-2_dp, 1.646407_dp,    -2.092942_dp,   &
           9.623100e-1_dp, -2.549600e-1_dp, 1.644040e-1_dp,  0.0_dp         &
                                                           /), shape=(/6,2/))
      !
      real(kind=dp), dimension(3,2), parameter :: g = reshape(source =      &
        (/   7.51198e-1_dp,  1.38889e-1_dp,  2.0886e-2_dp,                  &
             1.37691e-1_dp, -4.16667e-1_dp, -9.7486e-2_dp                   &
                                                           /), shape=(/3,2/))
      !
      real(kind=dp), dimension(2,2), parameter :: h = reshape(source =      &
        (/  -4.44444e-1_dp, -3.472e-3_dp,                                   &
             4.44444e-1_dp,  1.7361e-2_dp                  /), shape=(/2,2/))
      !
      real(kind=dp), parameter :: a = 2.2_dp,   b = -1.72_dp
      !
      integer       :: i, k1
      real(kind=dp) :: sum, sumg, sumh, x2, xm, xn
      logical       :: use_stop = .TRUE.
      !
      if (use_stop .and.  &
         (k < 0   .or.   k > 1)) then
         print *, "K must be either 0 or 1."
         stop     "vacpol_Lk(): program stop A."
      end if
      if (x  > two)  goto 3
      if (x == zero) goto 6
      !
      ! Use rational approximation for x < 2
      k1  = k + 1
      sum = zero
      xn  = one
      do  i = 1,6
         sum = sum + xn * f(i,k1)
         xn  = xn * x
      end do
      x2   = x * x
      sumg = g(1,k1) + x2 * (g(2,k1) + x2*g(3,k1))
      sumh = h(1,k1) + x2 * x2 * h(2,k1)
      xn   = log(x)
      sumg = xn*(sumg + xn*sumh)
      if (k == 0) goto 2
      sum  = sum + sumg
      goto 7
    2 sum  = sum + x * sumg
      goto 7
      !
    3 continue
      sum  = a + b / x
      if (k == 0) goto 4
      sum  = sum + (sum + b/x) / x
    4 sum  = sum / x
      xm   = -x
      sum  = sum * exp(xm)
      goto 7
    6 if (use_stop   .and.  k == 1) then
         print *, "Attempt to calculate function for zero argument "// &
                  "and K value of 1."
         stop     "vacpol_Lk(): program stop B."
      end if
      sum  = f(1,1)
    7 Lk   = sum
      !
   End function vacpol_Lk
