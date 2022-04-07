!======================================================================
      Subroutine Self_energy(SE)
!======================================================================
!     SE  =  Sum(orbitals) z^4/(pi*c^3*n^3) * F_QED(n,k,z)
!----------------------------------------------------------------------
      Use DBS_nuclear
      Use DBS_grid
      Use DBS_gauss
      Use dbsr_hf
      Use df_orbitals

      Implicit none
      Integer :: i,j
      Real(8) :: SE,S1,S2, ratio, r1,ar1,ar2,am1,am2,am3
      Real(8) :: Pcoef(ns),Qcoef(ns),fp_nucl(ns,ns),fq_nucl(ns,ns),v(ns)
      Real(8), external :: xAy, relci_qed_F

      SE = 0.d0; if(mode_SE.eq.0) Return

      Call Def_Znucl

      if(mode_SE.eq.1) then
       ZR_nucl = grw
      elseif(mode_SE.eq.3) then
       ZR_nucl = 0.d0
       Do i=1,nv; Do j=1,ks
        if(gr(i,j).le.0.0219) ZR_nucl(i,j) = grw(i,j)
       End do; End do
      end if

      Call ZINTYm (nv,ks,ksp,ksp,pbsp,pbsp,ZR_nucl,ns,fp_nucl)
      Call ZINTYm (nv,ks,ksq,ksq,qbsp,qbsp,ZR_nucl,ns,fq_nucl)

      SE = 0.d0
      Do i=1,nbf

       Call bdcwf_pq(nbs(i),kbs(i),z,Pcoef,Qcoef)

       S1 =  xAy(ns,ks,fp_nucl,'f',Pcoef,Pcoef) +  &
             xAy(ns,ks,fq_nucl,'f',Qcoef,Qcoef)
       S2 =  xAy(ns,ks,fp_nucl,'f',p(1,1,i),p(1,1,i)) +  &
             xAy(ns,ks,fq_nucl,'f',p(1,2,i),p(1,2,i))
       ratio = S2/S1
       if(debug.gt.0) &
       write(log,'(a,a,f15.5)') 'SE_ratio:  ',ebs(i),ratio
       Se = Se + ratio * relci_qed_F(nbs(i),kbs(i),z) * qsum(i) * &
                         z**4 / (pi*c_au**3*nbs(i)**3)
      End do

      End Subroutine Self_energy



!======================================================================
!  from RATIP
!======================================================================

   function relci_qed_F(n,kappa,Z)                             result(F)
   !--------------------------------------------------------------------
   ! Estimates the function  F (Z*\alpha) by using an interpolation of
   ! tabulated data from Mohr (1983) or Mohr and Kim (1992).
   !
   ! Calls: relci_qed_F_Mohr(), relci_qed_F_Mohr_Kim().
   !--------------------------------------------------------------------
      Use zconst
      Implicit none

      integer, intent(in)        :: n, kappa
      real(8), intent(in)  :: Z
      real(8)              :: F
      real(8), External    :: relci_qed_F_Mohr, relci_qed_F_Mohr_Kim

      if (n <= 2) then
         F = relci_qed_F_Mohr(n,kappa,z)
      else if (3 <= n   .and.   n <= 7) then
         F = relci_qed_F_Mohr_Kim(n,kappa,z)
      else
         F = zero
      end if
      !
   end function relci_qed_F
   !
   !
   function relci_qed_F_Klarsfeld(n,kappa,Z)                   result(F)
   !--------------------------------------------------------------------
   ! Estimates the function  F (Z*\alpha) by using a series expansion
   ! from S Klarsfeld and A Maquet, Physics Letters  43B (1973) 201,
   ! Eqs (1) and (2) and the table of Bethe logarithms. The
   ! vacuum-polarization contribution in Eq (2) is omitted.
   ! This procedure is adapted from RCI92 of GRASP92, written
   ! by Farid A Parpia, to the Fortran 95 standard.
   !
   ! This procedure is not used in the current version.
   !--------------------------------------------------------------------
      Use zconst
      Implicit none

      integer, intent(in)               :: n, kappa
      real(8), intent(in)         :: Z
      real(8)                     :: F

      Integer, External :: l_kappa

      real(kind=dp), dimension(36), parameter :: bethe = &
         (/ 2.9841285_dp,   2.8117699_dp,  -0.0300167_dp,   2.7676636_dp, &
           -0.0381902_dp,  -0.0052321_dp,   2.7498118_dp,  -0.0419549_dp, &
           -0.0067409_dp,  -0.0017337_dp,   2.7408237_dp,  -0.0440347_dp, &
           -0.0076008_dp,  -0.0022022_dp,  -0.0007721_dp,   2.7356642_dp, &
           -0.0453122_dp,  -0.0081472_dp,  -0.0025022_dp,  -0.0009628_dp, &
           -0.0004079_dp,   2.7324291_dp,  -0.0461552_dp,  -0.0085192_dp, &
           -0.0027091_dp,  -0.0010945_dp,  -0.0004997_dp,  -0.0002409_dp, &
            2.7302673_dp,  -0.0467413_dp,  -0.0087850_dp,  -0.0028591_dp, &
           -0.0011904_dp,  -0.0005665_dp,  -0.0002904_dp,  -0.0001539_dp /)
      !
      real(kind=dp), parameter :: C401 = 11.0_dp/24.0_dp,                 &
                               C402 = 3.0_dp/8.0_dp, ovlfac = 4.0_dp/3.0_dp
      !
      integer       :: l, loc
      real(kind=dp) :: bethel, factor, term
      !
      ! Ensure that the principal quantum number is in range
      if (n < 1   .or.   n > 8) then
         print *, "Principal quantum number,",n,", should be in the range 1-8."
         stop     "relci_qed_F_Klarsfeld(): program stop A."
      end if
      !
      l = l_kappa (kappa)
      if (l > n-1) then
         print *, "Kappa = ",kappa," is out of range for n = ",n,"."
         stop     "relci_qed_F_Klarsfeld(): program stop B."
      end if
      !
      ! Find the appropriate entry in the table
      loc    = (n*n-n)/2+l+1
      bethel = bethe(loc)
      !
      ! Determine the quantity in square brackets in eq.(1) of
      ! Klarsfeld and Maquet
      term = -bethel
      !
      if (kappa > 0) then
         term = term - c402 / (l*(l+l+one))
      else
         term = term + c402 / ((l+one)*(l+l+one))
         if (kappa == -1) then
            factor = log (Z/c_au)              !  ???
            factor = - (factor + factor)
            term   = term + factor + c401
         end if
      end if
      !
      F = ovlfac * term
      !
   end function relci_qed_F_Klarsfeld
   !
   !
   function relci_qed_F_Mohr(n,kappa,Z)                        result(F)
   !--------------------------------------------------------------------
   ! Computes the function  F (Z*\alpha) for the  1s  2s  2p-  2p
   ! symmetries by interpolating in, or extrapolating from, the table
   ! due to  P J Mohr. See  P J Mohr, At Data Nucl Data Tables 29
   ! (1983) 453.
   ! This procedure is adapted from RCI92 of GRASP92, written
   ! by Farid A Parpia, to the Fortran 95 standard.
   !
   ! Calls:
   !--------------------------------------------------------------------
      Use zconst
      Implicit none

      integer, intent(in) :: n, kappa
      real(8), intent(in) :: Z
      real(8)             :: F, value
      !
      ! Number of data points
      integer, parameter        :: numval = 12
      real(kind=dp), parameter  :: accy   = 1.0e-2_dp
      !
      ! 1s data
      real(kind=dp), dimension(numval), parameter :: val1s =     &
         (/   10.3168_dp,  4.6540_dp,   3.2460_dp,   2.5519_dp,  &
               2.1351_dp,  1.8644_dp,   1.6838_dp,   1.5675_dp,  &
               1.5032_dp,  1.4880_dp,   1.5317_dp,   1.6614_dp  /)
      ! 2s data
      real(kind=dp), dimension(numval), parameter :: val2s =     &
         (/   10.5468_dp,  4.8930_dp,   3.5063_dp,   2.8391_dp,  &
               2.4550_dp,  2.2244_dp,   2.0948_dp,   2.0435_dp,  &
               2.0650_dp,  2.1690_dp,   2.3870_dp,   2.7980_dp  /)
      ! 2p- data
      real(kind=dp), dimension(numval), parameter :: val2p1 =    &
         (/   -0.1264_dp, -0.1145_dp,  -0.0922_dp,  -0.0641_dp,  &
              -0.0308_dp,  0.0082_dp,   0.0549_dp,   0.1129_dp,  &
               0.1884_dp,  0.2934_dp,   0.4530_dp,   0.7250_dp  /)
      ! 2p data
      real(kind=dp), dimension(numval), parameter :: val2p3 =    &
         (/    0.1235_dp,  0.1303_dp,   0.1436_dp,   0.1604_dp,  &
               0.1794_dp,  0.1999_dp,   0.2215_dp,   0.2440_dp,  &
               0.2671_dp,  0.2906_dp,   0.3141_dp,   0.3367_dp  /)
      ! Z data
      real(kind=dp), dimension(numval), parameter :: arg =       &
         (/    1.0_dp,    10.0_dp,     20.0_dp,     30.0_dp,     &
              40.0_dp,    50.0_dp,     60.0_dp,     70.0_dp,     &
              80.0_dp,    90.0_dp,    100.0_dp,    110.0_dp     /)
      !
      ! Interpolate or issue error message as appropriate
      if (n == 1) then
         select case(kappa)
         case(-1)
            call interpolation_aitken(arg,val1s,numval,z,value,accy)
         case default
            print *, "Principal quantum number is ",n," and kappa",kappa
            stop     "relci_qed_F_Mohr(): program stop A."
         end select
      else if (n == 2) then
         select case(kappa)
         case(-1)
            call interpolation_aitken(arg,val2s, numval,z,value,accy)
         case(1)
            call interpolation_aitken(arg,val2p1,numval,z,value,accy)
         case(-2)
            call interpolation_aitken(arg,val2p3,numval,z,value,accy)
         case default
            print *, "Principal quantum number is ",n," and kappa",kappa
            stop     "relci_qed_F_Mohr(): program stop B."
         end select
      else
         print *, "Principal quantum number, ",n,"should be either 1 or 2."
         stop     "relci_qed_F_Mohr(): program stop B."
      end if
      !
      F = value
      !
   end function relci_qed_F_Mohr
   !
   !
   function relci_qed_F_Mohr_Kim(n,kappa,Z)                    result(F)
   !--------------------------------------------------------------------
   ! Computes the function  F (Z*\alpha) for the  1s  2s  2p-  2p
   ! symmetries by interpolating in, or extrapolating from, the table
   ! due to  P J Mohr and Y-K Kim. See P J Mohr and Y-K Kim,
   ! Phys. Rev. A45 (1992) 2723.
   !
   ! Since no values are given for Z = 1, these values were estimated
   ! below by (graphical) extrapolation.
   !
   ! Calls: interpolation_aitken().
   !--------------------------------------------------------------------
      Use zconst
      Implicit none

      integer, intent(in)         :: n, kappa
      real(8), intent(in)         :: Z
      real(8)                     :: F, value
      !
      ! Number of data points
      integer, parameter       :: numval = 12
      real(kind=dp), parameter :: accy   = 1.0e-2_dp
      !
      ! 3s data
      real(kind=dp), dimension(numval), parameter :: val3s =     &
         (/   10.5000_dp,  4.9524_dp,   3.5633_dp,   2.8940_dp,  &
               2.5083_dp,  2.2757_dp,   2.1431_dp,   2.0874_dp,  &
               2.1018_dp,  2.1935_dp,   2.3897_dp,   2.7609_dp  /)
      ! 3p- data
      real(kind=dp), dimension(numval), parameter :: val3p1 =    &
         (/   -0.1300_dp, -0.1021_dp,  -0.0760_dp,  -0.0430_dp,  &
               0.0041_dp,  0.0414_dp,   0.0956_dp,   0.1623_dp,  &
               0.2483_dp,  0.3660_dp,   0.5408_dp,   0.8322_dp  /)
      ! 3p data
      real(kind=dp), dimension(numval), parameter :: val3p3 =    &
         (/    0.1250_dp,  0.1421_dp,   0.1572_dp,   0.1761_dp,  &
               0.1977_dp,  0.2214_dp,   0.2470_dp,   0.2745_dp,  &
               0.3038_dp,  0.3350_dp,   0.3679_dp,   0.4020_dp  /)
      ! 3d- data
      real(kind=dp), dimension(numval), parameter :: val3d3 =    &
         (/   -0.0440_dp, -0.0428_dp,  -0.0420_dp,  -0.0410_dp,  &
              -0.0396_dp, -0.0378_dp,  -0.0353_dp,  -0.0321_dp,  &
              -0.0279_dp, -0.0225_dp,  -0.0154_dp,  -0.0062_dp  /)
      !
      ! 4s data
      real(kind=dp), dimension(numval), parameter :: val4s =     &
         (/   10.5000_dp,  4.9749_dp,   3.5834_dp,   2.9110_dp,  &
               2.5215_dp,  2.2842_dp,   2.1455_dp,   2.0814_dp,  &
               2.0840_dp,  2.1582_dp,   2.3262_dp,   2.6484_dp  /)
      ! 4p- data
      real(kind=dp), dimension(numval), parameter :: val4p1 =    &
         (/   -0.1200_dp, -0.0963_dp,  -0.0690_dp,  -0.0344_dp,  &
               0.0064_dp,  0.0538_dp,   0.1098_dp,   0.1780_dp,  &
               0.2649_dp,  0.3819_dp,   0.5525_dp,   0.8311_dp  /)
      ! 4p data
      real(kind=dp), dimension(numval), parameter :: val4p3 =    &
         (/    0.1250_dp,  0.1477_dp,   0.1630_dp,   0.1827_dp,  &
               0.2052_dp,  0.2299_dp,   0.2568_dp,   0.2858_dp,  &
               0.3170_dp,  0.3507_dp,   0.3868_dp,   0.4247_dp  /)
      ! 4d- data
      real(kind=dp), dimension(numval), parameter :: val4d3 =    &
         (/   -0.0410_dp, -0.0403_dp,  -0.0399_dp,  -0.0387_dp,  &
              -0.0371_dp, -0.0348_dp,  -0.0317_dp,  -0.0276_dp,  &
              -0.0222_dp, -0.0149_dp,  -0.0053_dp,   0.0074_dp  /)
      !
      ! 5s data
      real(kind=dp), dimension(numval), parameter :: val5s =     &
         (/   10.5000_dp,  4.9858_dp,   3.5923_dp,   2.9173_dp,  &
               2.5246_dp,  2.2833_dp,   2.1395_dp,   2.0686_dp,  &
               2.0619_dp,  2.1225_dp,   2.2696_dp,   2.5566_dp  /)
      ! 5p- data
      real(kind=dp), dimension(numval), parameter :: val5p1 =    &
         (/   -0.1200_dp, -0.0933_dp,  -0.0652_dp,  -0.0299_dp,  &
               0.0116_dp,  0.0597_dp,   0.1161_dp,   0.1843_dp,  &
               0.2703_dp,  0.3848_dp,   0.5497_dp,   0.8150_dp  /)
      ! 5p data
      real(kind=dp), dimension(numval), parameter :: val5p3 =    &
         (/    0.1300_dp,  0.1502_dp,   0.1662_dp,   0.1861_dp,  &
               0.2089_dp,  0.2341_dp,   0.2614_dp,   0.2910_dp,  &
               0.3229_dp,  0.3574_dp,   0.3946_dp,   0.4338_dp  /)
      ! 5d- data
      real(kind=dp), dimension(numval), parameter :: val5d3 =    &
         (/   -0.0405_dp, -0.0396_dp,  -0.0387_dp,  -0.0374_dp,  &
              -0.0356_dp, -0.0331_dp,  -0.0297_dp,  -0.0252_dp,  &
              -0.0190_dp, -0.0108_dp,   0.0001_dp,   0.0145_dp  /)
      ! Z data
      real(kind=dp), dimension(numval), parameter :: arg =       &
         (/    1.0_dp,    10.0_dp,     20.0_dp,     30.0_dp,     &
              40.0_dp,    50.0_dp,     60.0_dp,     70.0_dp,     &
              80.0_dp,    90.0_dp,    100.0_dp,    110.0_dp     /)
      !
      ! Interpolate or issue error message as appropriate
      if (n == 3) then
         select case(kappa)
         case(-1)
            call interpolation_aitken(arg,val3s,numval,z,value,accy)
         case(1)
            call interpolation_aitken(arg,val3p1,numval,z,value,accy)
         case(-2)
            call interpolation_aitken(arg,val3p3,numval,z,value,accy)
         case(2)
            call interpolation_aitken(arg,val3d3,numval,z,value,accy)
         case(-3)
            call interpolation_aitken(arg,val3d3,numval,z,value,accy)
            value = 0.9_dp * value
            value = zero
         case default
            print *, "Principal quantum number is ",n," and kappa",kappa
            stop     "relci_qed_F_Mohr_Kim(): program stop A."
         end select
      else if (n == 4) then
         select case(kappa)
         case(-1)
            call interpolation_aitken(arg,val4s,numval,z,value,accy)
         case(1)
            call interpolation_aitken(arg,val4p1,numval,z,value,accy)
         case(-2)
            call interpolation_aitken(arg,val4p3,numval,z,value,accy)
         case(2)
            call interpolation_aitken(arg,val4d3,numval,z,value,accy)
         case(-3)
            call interpolation_aitken(arg,val4d3,numval,z,value,accy)
            value = 0.9_dp * value
            value = zero
         case(3)
            value = zero
         case(-4)
            value = zero
         case default
            print *, "Principal quantum number is ",n," and kappa",kappa
            stop     "relci_qed_F_Mohr_Kim(): program stop B."
         end select
      else if (5 <= n  .and.   n <= 7) then
         select case(kappa)
         case(-1)
            call interpolation_aitken(arg,val5s,numval,z,value,accy)
         case(1)
            call interpolation_aitken(arg,val5p1,numval,z,value,accy)
         case(-2)
            call interpolation_aitken(arg,val5p3,numval,z,value,accy)
         case(2)
            call interpolation_aitken(arg,val5d3,numval,z,value,accy)
         case(-3)
            call interpolation_aitken(arg,val4d3,numval,z,value,accy)
            value = 0.9_dp * value
            value = zero
         case(3:6)
            value = zero
         case(-7:-4)
            value = zero
         case default
            print *, "Principal quantum number is ",n," and kappa",kappa
            stop     "relci_qed_F_Mohr_Kim(): program stop C."
         end select
      else
         print *, "Principal quantum number, ",n,"should be in the interval "//&
                  "3 <= n <= 7."
         stop     "relci_qed_F_Mohr(): program stop D."
      end if
      !
      F = value
      !
   end function relci_qed_F_Mohr_Kim


   subroutine interpolation_aitken(xarr,yarr,narr,xval,yval,accy)
   !--------------------------------------------------------------------
   ! This routine returns  yval  as the functions value of  xval
   ! by interpolating a pair of arrays xarr(1:narr), yarr(1:narr),
   ! that tabulate a  function.  Aitken's algorithm is used. See, for
   ! instance, F B Hildebrand, Introduction to Numerical Analysis,
   ! 2nd ed., McGraw-Hill, New York, NY, 1974. accy is the  desired
   ! accuracy of the estimate: a warning message is issued if this is
   ! not achieved.  A warning is also issued when the routine is
   ! extrapolating.
   ! This procedures is adapted to Fortran 90/95 from the routine
   ! interp() of GRASP92 which originally was written by F A Parpia.
   !--------------------------------------------------------------------

      Use zconst
      Implicit none

      integer,       intent(in)  :: narr
      real(kind=dp), intent(in)  :: xval
      real(kind=dp), intent(out) :: yval
      real(kind=dp), dimension(narr), intent(in) :: xarr, yarr
      real(kind=dp), intent(in)  :: accy
      !
      ! mxord is the maximum order of the interpolation
      integer, parameter                            :: mxord = 11
      logical, dimension(2*mxord+2)                 :: used
      real(kind=dp), dimension(mxord)               :: dx, x, est
      real(kind=dp), dimension((mxord*(mxord+1))/2) :: poly
      !
      logical       :: set
      integer       :: ibest, ilirok, ildiag, ilothr, irow, k, lhi, llo, llr, locnxt, nrsthi, nrstlo
      real(kind=dp) :: debe, debeb, diff, difft
      !
      ! Determine the nearest two XARR entries bounding XVAL
      if (xval < xarr(1)) then
         nrstlo = 1;   nrsthi = 1
         print *, "interpolation_aitken(): Extrapolating, not interpolating."
      elseif (xval > xarr(narr)) then
         nrstlo = narr;   nrsthi = narr
         print *, "interpolation_aitken(): Extrapolating, not interpolating."
      else
         k = 0
       1 k = k+1
         if (xarr(k) < xval) then
            nrstlo = k
            goto 1
         else
            nrsthi = k
         end if
      end if
      !
      ! Clear relevant piece of use-indicator array
      llo = max(nrstlo-mxord,   1)
      lhi = min(nrsthi+mxord,narr)
      llr = llo - 1
      do  k = llo,lhi
         used(k-llr) = .false.
      end do
      !
      ! Determine next-nearest XARR entry
      do  irow = 1,mxord
         llo = max(nrstlo-irow+1,   1)
         lhi = min(nrsthi+irow-1,narr)
         set = .false.
         do  k = llo,lhi
            if (.not.used(k-llr)) then
               if (.not.set) then
                  diff = xarr(k) - xval
                  locnxt = k
                  set = .true.
               else
                  difft = xarr(k) - xval
                  if (abs(difft) < abs(diff)) then
                     diff = difft
                     locnxt = k
                  end if
               end if
            end if
         end do
         used(locnxt-llr) = .true.
         x(irow)  = xarr(locnxt)
         dx(irow) = diff
         !
         ! Fill table for this row
         do  k = 1,irow
            ilirok = iloc(irow,k)
            if (k == 1) then
               poly(ilirok) = yarr(locnxt)
            else
               ildiag       = iloc(k-1,k-1)
               ilothr       = iloc(irow,k-1)
               !!x print *, "INTERPOLATION_AITKEN a"
               !!x print *, "k,irow,x(1:irow) = ",k,irow,x(1:irow)
               poly(ilirok) = (poly(ildiag)*dx(irow) - poly(ilothr)*dx(k-1)) &
                              / (x(irow)-x(k-1))
               !!x print *, "INTERPOLATION_AITKEN b"
            endif
         end do
         !
         ! Pick off the diagonal element
         ildiag    = iloc(irow,irow)
         est(irow) = poly(ildiag)
      end do
      !
      ! Now the estimate vector is filled in, so obtain the best estimate
      debeb = abs((est(2) - est(1)) / est(2))
      ibest = 2
      do  irow = 3,mxord
         debe = abs((est(irow) - est(irow-1)) / est(irow))
         if (debe < debeb) then
            debeb = debe
            ibest = irow
         end if
      end do
      yval = est(ibest)
      !
!zoi      if (present(accy)) then
!         if (debeb > accy) then
!            write(*,2) debeb, accy
!          2 format( "interpolation_aitken(): Accuracy of interpolation (", &
!                    e10.3,") is below input criterion (",e10.3,").")
!         end if
!zoi      end if
      !
      contains
         !
         function iloc (ind1,ind2)                           result(loc)
         !--------------------------------------------------------------
         ! This internal function dispenses with the need for a
         ! two-dimensional array for the interpolation. It replaces a
         ! statement function in the original code.
         !--------------------------------------------------------------
         !
         integer, intent(in) :: ind1, ind2
         integer             :: loc
         !
         loc = (ind1*(ind1-1)) / 2 + ind2
         !
         end function iloc
         !
   end subroutine interpolation_aitken
