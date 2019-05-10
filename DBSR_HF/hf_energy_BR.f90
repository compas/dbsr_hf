!======================================================================
      Subroutine Breit_energy (Ebreit)
!======================================================================
!     Breit corrections to Dirac-Coulomb energy
!     for the generalized configuration
!----------------------------------------------------------------------
      Use hf_energy
      Use dbsr_hf
      Use df_orbitals, only: nbf, kbs, ebs
      Use DBS_grid
      Use rk4_data      

      Implicit none
      Integer :: i,j,k, ki,kj, k1,k2,k3,k4, i1,i2,i3,i4, v
      Real(8) :: S(8), Ebreit
      Real(8), external :: SMU, sk_ppqq

      nrk = 0
      Do i = 1,nbf;  ki=kbs(i)
      Do j = 1,nbf;  kj=kbs(j) 
 
      if(i.ge.j) then
       k1=ki; k2=kj; k3=ki; k4=kj; i1=i; i2=j; i3=i; i4=j
      else
       k1=ki; k2=kj; k3=kj; k4=ki; i1=i; i2=j; i3=j; i4=i
      end if

      Do k = 0,kmax
       if(coef(i,j,k).eq.0.d0) Cycle
       Do v = k-1,k+1
        if(SMU(k1,k2,k3,k4,k,v,S).eq.0.d0) Cycle
        S = S * coef(i,j,k)
        Call Add_rk4_data(1,v,i1*ibi+i2,i3*ibi+i4,S(1)) 
        Call Add_rk4_data(1,v,i2*ibi+i1,i4*ibi+i3,S(2)) 
        Call Add_rk4_data(1,v,i3*ibi+i4,i1*ibi+i2,S(3)) 
        Call Add_rk4_data(1,v,i4*ibi+i3,i2*ibi+i1,S(4)) 
        Call Add_rk4_data(1,v,i1*ibi+i4,i3*ibi+i2,S(5)) 
        Call Add_rk4_data(1,v,i4*ibi+i1,i2*ibi+i3,S(6)) 
        Call Add_rk4_data(1,v,i3*ibi+i2,i1*ibi+i4,S(7)) 
        Call Add_rk4_data(1,v,i2*ibi+i3,i4*ibi+i1,S(8)) 
       End do
      End do

      End do; End do

      Ebreit = 0.d0
      Do v = 1,nrk; if(abs(crk(v)).lt.eps_c) Cycle
       k = kr2(v)
       i1 = kr3(v)/ibi; i2=mod(kr3(v),ibi);  i3 = kr4(v)/ibi; i4=mod(kr4(v),ibi) 
       Ebreit = Ebreit + crk(v) * sk_ppqq(i1,i2,i3,i4,k)
      End do

      End Subroutine Breit_energy



!=======================================================================
      Real(8) Function SMU (ka,kb,kc,kd,L,v,S)
!=======================================================================
! ... Computes the coefficient for Breit operator
! ... Grant and Pyper, J.Phys.B 9, 761 (1976)
! ... Grant book: p.345
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: ka,kb,kc,kd,L,v
      Real(8), intent(out) :: S(8)
      Integer :: i, la,lb,lc,ld, K,Kp, vv
      Integer, External :: ITRA, l_kappa
      Real(8) :: b,c,bp,cp

      SMU = 0.d0; S = 0.d0; if(L.lt.0) Return;  if(v.lt.0) Return 

      la=l_kappa(ka); lc=l_kappa(kc)
      if(mod(la+lc+v,2).eq.0) Return

      lb=l_kappa(kb); ld=l_kappa(kd)
      if(mod(lb+ld+v,2).eq.0) Return

      if(v.eq.L) then

       if(L.eq.0) Return
       SMU = -(ka+kc)*(kb+kd); SMU=SMU/L/(L+1); S=SMU

      else

       K=kc-ka; Kp=kd-kb
       
       if(v.eq.L-1) then

        bp=L+1; bp=bp/2/(L+L-1); cp=-(L-2); cp=cp/2/L/(L+L-1)

        S(1) = (L+K )*(bp+cp*KP)
        S(2) = (L+KP)*(bp+cp*K )
        S(3) = (L-K )*(bp-cp*KP)
        S(4) = (L-KP)*(bp-cp*K )
        S(5) =-(L+K )*(bp-cp*KP)
        S(6) =-(L-KP)*(bp+cp*K )
        S(7) =-(L-K )*(bp+cp*KP)
        S(8) =-(L+KP)*(bp-cp*K )

       elseif(v.eq.L+1) then

        i = -L-1
        b=L; b=b/2/(L+L+3); c=(L+3); c=c/(L+L+2)/(L+L+3)

        S(1) = (i+KP)*(b+c*K )
        S(2) = (i+K )*(b+c*KP)
        S(3) = (i-KP)*(b-c*K )
        S(4) = (i-K )*(b-c*KP)
        S(5) =-(i-KP)*(b+c*K )
        S(6) =-(i+K )*(b-c*KP)
        S(7) =-(i+KP)*(b-c*K )
        S(8) =-(i-K )*(b+c*KP)

       end if

        Do i=1,8; SMU = SMU + abs(S(i)); End do

      end if

      End Function SMU 


!======================================================================
      Subroutine Breit_coefs
!======================================================================
!     Breit-interaction angular coefficients
!----------------------------------------------------------------------
!     list of coefficients will contain:
!
! 1.  k + atype *1000;   atype=0 =>  sk_ppqq;  atype=1 => sk_pqqp
! 2   int = 1 -->  INT( i,j; i,j)
!     int = 2 -->  INT( j,i; j,j)
!     int = 3 -->  INT( i,j; j,i)
!     int = 4 -->  INT( j,i; i,j)
! 3.  i 
! 4.  j 
!----------------------------------------------------------------------
      Use dbsr_hf
      Use hf_energy
      Use df_orbitals, only: nbf,kbs

      Implicit none
      Integer :: i,j,k, ki,kj, v
      Real(8) :: S(8), C
      Real(8), external :: SMU

      Call alloc_rk4_data(0)

      Do i = 1,nbf;  ki=kbs(i)
      Do j = 1,nbf;  kj=kbs(j) 
 
      Do k = 0,kmax
       if(coef(i,j,k).eq.0.d0) Cycle
       Do v = k-1,k+1
         if(i.ge.j) then                             !  Rk(i,j,i,j)
          if(SMU(ki,kj,ki,kj,k,v,S).eq.0.d0) Cycle
          S = S * coef(i,j,k)
          C = Sum(S) 
          if(i.eq.j) then
           C = Sum(S) 
           Call Add_rk4_data(v,1,i,i,C) 
          else
           C = S(1)+S(3)+S(5)+S(7)
           Call Add_rk4_data(v,1,i,j,C) 
           C = S(2)+S(4)+S(6)+S(8)
           Call Add_rk4_data(v,2,i,j,C) 
          end if
         else                                       !  Rk(i,j,j,i)
          if(SMU(ki,kj,kj,ki,k,v,S).eq.0.d0) Cycle
          S = S * coef(i,j,k)
          C = S(1)+S(4)
          Call Add_rk4_data(v,3,i,j,C) 
          C = S(2)+S(3)
          Call Add_rk4_data(v,4,i,j,C) 
          C = S(5)+S(6)
          Call Add_rk4_data(v+1000,3,i,j,C) 
          C = S(7)+S(8)
          Call Add_rk4_data(v+1000,4,i,j,C) 
         end if  
       End do      ! v
      End do       ! k

      End do; End do  ! i,j

      End Subroutine Breit_coefs

!======================================================================
      Subroutine Breit_energy_coef (Ebreit)
!======================================================================
!     Breit corrections to Dirac-Coulomb energy
!     for the generalized configuration
!----------------------------------------------------------------------
      Use rk4_data      
      Use dbsr_hf,  only: kmax, eps_c
      Use DBS_grid, only: ns,ks

      Implicit none
      Integer :: i,j,k,m, int, atype, i1,i2,i3,i4
      Real(8) :: Ebreit
      Real(8), external :: sk_ppqq, sk_pqqp

      Call alloc_DBS_sk_integrals(ns,ks,0,kmax,2)

      Ebreit = 0.d0
      Do m = 1,nrk; if(abs(crk(m)).lt.eps_c) Cycle
       k = mod(kr1(m),1000); atype=kr1(m)/1000
       int = kr2(m); i = kr3(m); j = kr4(m)
       Select case(int)
        Case(1);   i1 = i; i2 = j;  i3 = i; i4 = j 
        Case(2);   i1 = j; i2 = i;  i3 = j; i4 = i 
        Case(3);   i1 = i; i2 = j;  i3 = j; i4 = i 
        Case(4);   i1 = j; i2 = i;  i3 = i; i4 = j 
       End Select
       if(atype.eq.0)  Ebreit = Ebreit + crk(m) * sk_ppqq(i1,i2,i3,i4,k)

       if(atype.eq.1)  Ebreit = Ebreit + crk(m) * sk_pqqp(i1,i2,i3,i4,k)
      End do

      End Subroutine Breit_energy_coef

