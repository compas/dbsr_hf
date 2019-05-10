!=======================================================================
      Module hf_energy
!=======================================================================
!     This module defines the energy expression as sum  of
!     one-electron       Sum_i  qsum(i) valueI (i)   
!     and two-electron   Sum_i,j,k  coef(i,j,k)  value(i,j,k)  
!     integrals.   value(i,j,k) ->  Fk(i,j) or Gk(i,j) integrals,
!     which is defined by subroutines a(i,j,k) and b(i,j,k). 
!     Note that coef(i,j,k) with i > = j refer to direct Fk-integrals
!     whereas coef(i,j,k) with i < j refer to exchage Gk-integrals. 
!----------------------------------------------------------------------
      Implicit none

      Real(8), allocatable :: coef (:,:,:), value(:,:,:), valueI(:)

      End Module hf_energy


!=======================================================================
      Real(8) Function a(i,j,k)
!=======================================================================
!     coefficient to the potential for electron i of y^k (i,j)
!-----------------------------------------------------------------------
      Use hf_energy
      Use df_orbitals, only: jbs

      Implicit none
      Integer, intent(in)  :: i,j,k
      a = 0.d0
      if(mod(k,2).ne.0) Return
      if(k > (jbs(i)+jbs(j))/2 ) Return
      if(k < 0) Return
      a = coef(max(i,j),min(i,j),k)

      End Function a


!=======================================================================
      Real(8) Function b(i,j,k) 
!=======================================================================
!     determine the coefficient of the y^k (i,j) p(j) term in the exchange
!     expression for electron i
!-----------------------------------------------------------------------
      Use hf_energy
      Use df_orbitals, only: lbs,jbs

      Implicit none
      Integer, Intent(IN) :: i,j,k

      b = 0.d0
      if(i.eq.j) Return
      if(mod(k+lbs(i)+lbs(j),2) /= 0) Return
      if( k < abs(jbs(i)-jbs(j))/2 ) Return
      if( k > (jbs(i)+jbs(j))/2 ) Return 
      b = coef(min(i,j),max(i,j),k)
    
      End Function b


!======================================================================
      Subroutine Def_energy_coef
!======================================================================
!     define angular coefficients for energy expression
!----------------------------------------------------------------------
      Use hf_energy
      Use dbsr_hf
      Use df_orbitals, only: nbf,jbs
      Integer :: i,j,k

      kmax = maxval(jbs(1:nbf))
      if(allocated(coef)) Deallocate(coef);  Allocate(coef(nbf,nbf,0:kmax))
      coef = 0.d0
      if(allocated(value)) Deallocate(value);  Allocate(value(nbf,nbf,0:kmax))
      value = 0.d0
      if(allocated(valueI)) Deallocate(valueI);  Allocate(valueI(nbf))
      valueI = 0.d0

      if(term.eq.'jj') then
       if(ncore.gt.0) Call av_energy_coef(ncore)
       Call jj_energy_coef(0)
      else
       Call av_energy_coef(nbf)
      end if

      if(mbreit.gt.0) Call Breit_coefs

      End Subroutine Def_energy_coef


!======================================================================
      Subroutine av_energy_coef(jj)
!======================================================================
!     Angular coefficients for energy in average-term approximation
!
!     Into shell, only direct interaction:  
!
!     q(q-1)/2 [ F0(a,a)  -  (2j+1)/2j  { j  k   j } ^2  Fk(a,a)  ]
!                                       {1/2 0 -1/2}
!
!     Between shells: direct only k=0
!
!     q_a*q_b/2 [ F0(a,b) -   {j_a  k j_b} ^2  Gk(a,b)   ]
!                             {1/2 0 -1/2}            
!
!     We devide on (2ja+1)(2jb+1) due to using Function Cjkj instead 3J-symbol 
!----------------------------------------------------------------------
      Use hf_energy
      Use df_orbitals

      Implicit none
      Integer, intent(in) :: jj
      Integer :: i,j, k
      Real(8) :: c
      Real(8), external :: Cjkj, WW 
	  
      if(jj.le.0) Return
      Do i = 1,nbf
      Do j = 1,i
       if(j.gt.jj) Exit
       if(i.eq.j) then
        c = WW(i,j)/2.d0       ! qsum(i)*(qsum(i)-1.d0)/2
        coef(i,i,0) = c
        Do k=1,(jbs(i)+jbs(i))/2
         coef(i,i,k)=-c*Cjkj(jbs(i),k,jbs(i))**2/(jbs(i)*(jbs(i)+1))
        End do
       else
        c = WW(i,j)            ! qsum(i)*qsum(j)
        coef(i,j,0)=c
        Do k=iabs(jbs(i)-jbs(j))/2,(jbs(i)+jbs(j))/2
         coef(j,i,k)=-c*Cjkj(jbs(i),k,jbs(j))**2/((jbs(i)+1)*(jbs(j)+1))
        End do
       end if
      End do; End do
	  
      End Subroutine av_energy_coef 


!======================================================================
      Real(8) Function WW(i,j)
!======================================================================
! ... generalized value of WC(ic)*WC(jc) - product of expansion
! ... coeficients for configurations ic and jc
!======================================================================
      Use dbsr_hf
      Use DF_orbitals
      
      Implicit none
      Integer, Intent(in) :: i,j 
      Integer :: ic

      if(nconf.eq.1) then
       if(i.eq.j) WW = qsum(i)*(qsum(i)-1.d0)       
       if(i.ne.j) WW = qsum(i)*qsum(j)
       Return
      end if

      WW = 0.d0
      if(i.eq.j) then
       Do ic=1,nconf
        WW = WW + iqconf(ic,i)*(iqconf(ic,i)-1)*weight(ic)
       End do
      else
       Do ic=1,nconf
        WW = WW + iqconf(ic,i)*iqconf(ic,j)*weight(ic)
       End do
      end if

      End Function WW


!======================================================================
      Subroutine update_int(ii)
!======================================================================
!     Update the integral list for orbital "ii"
!     if ii=0  -  all orbitals
!----------------------------------------------------------------------
      Use hf_energy
      Use df_orbitals, only: nbf, kbs, l => lbs, p

      Implicit none
      Integer, intent(in) :: ii
      Integer :: i,j,k
      Real(8), external :: Vp_dhl, a, b, hf_fk, hf_gk,  hf_fk1, hf_gk1

      Do i = 1,nbf

       if(ii.ne.0.and.i.ne.ii) Cycle

       valueI(i) = Vp_dhl (kbs(i),p(1,1,i),kbs(i),p(1,1,i))
 
       Do j = 1,i  ! nbf 
        Do  k = 0,2*min0(l(i),l(j)),2
         if(a(i,j,k).eq.0.d0) Cycle
         value(i,j,k) = hf_fk(i,j,k)
       End do
       End do 
 
       Do  j = 1,i-1 ! nbf 
        Do k = abs(l(i)-l(j)),l(i)+l(j),2
         if(b(i,j,k).eq.0.d0) Cycle
         value(j,i,k) = hf_gk(i,j,k)
        End do
       End do

      End do

      End Subroutine update_int


!======================================================================
      Subroutine energy
!======================================================================
!     Determine the total energy of the state, together with                                                                
!     one-electron and two-electron (direct and exchange) parts
!----------------------------------------------------------------------
      Use hf_energy
      Use dbsr_hf,     only: Etotal, E1body, E2body
      Use df_orbitals, only: nbf, kbs, l => lbs, qsum

      Implicit none
      Integer :: i,j,k
      Real(8) :: c, S
      Real(8), external :: a,b

! ... direct interaction and one-electron energies

      Etotal = 0.d0; E1body =0.d0; E2body =0.d0
      Do i = 1,nbf
       S = valueI(i)
       etotal = etotal + qsum(i)*S
       E1body = E1body + qsum(i)*S
       Do j = 1,i
        Do  k = 0,2*min0(l(i),l(j)),2
         c = a(i,j,k); if(c.eq.0.d0) Cycle
         S = value(i,j,k)                        
         etotal = etotal + c*S
         E2body = E2body + c*S
        End do
       End do 

! ... exchange interaction
 
       Do  j = 1,i-1
        Do k = abs(l(i)-l(j)),l(i)+l(j),2
         c = b(i,j,k); ; if(c.eq.0.d0) Cycle
         S = value(j,i,k)                         
         etotal = etotal + c*S
         E2body = E2body + c*S
        End do
       End do

      End do  !  over orbitals, i

      End Subroutine energy



!======================================================================
      Subroutine jj_energy_coef(jc)
!======================================================================
!     This routine gets agular coefficients for set of atomic states
!     in jj coupling (atomic states are taken from c-file, unit nuc)
!----------------------------------------------------------------------
      Use dbsr_hf
      Use hf_energy
      Use df_orbitals

      Implicit none
      Integer :: i,j,k, ic,jc, ii,jj, km
      Real(8) :: c
      Real(8), allocatable :: coefs(:,:,:)
      Integer, external :: Ifind_orb 

      rewind(nuc)
      ic = 0
      Do 
       read(nuc,'(a)') CONFIG
       if(CONFIG(6:6).ne.'(') Cycle
       read(nuc,'(a)') SHELLJ
       read(nuc,'(5x,a)') INTRAJ
       ic = ic + 1
       if(jc.ne.0.and.ic.ne.jc) Cycle

       Call Decode_cjj(CONFIG,SHELLJ,INTRAJ,no,nn,kn,ln,jn,iq,in,&
                                            Jshell,Vshell,Jintra)
       Allocate(coefs(no,no,0:kmax))

       Call coef_1conf(no,ln,jn,iq,Jshell,Vshell,Jintra,kmax,coefs)

       Do i=1,no;  ii = Ifind_orb(nn(i),kn(i),0)
       Do j=1,no;  jj = Ifind_orb(nn(j),kn(j),0)
       Do k=0,kmax
        C = coefs(i,j,k); if(abs(C).lt.eps_c) C=0.d0
        if(C.eq.0.d0) Cycle
        coef(ii,jj,k) = coef(ii,jj,k) + C * weight(ic)  
       End do; End do; End do
       Deallocate(coefs)

       if(jc.ne.0.and.ic.eq.jc) Exit
       if(ic.eq.nconf) Exit
      End do  ! over configurations
  
      End Subroutine jj_energy_coef 


