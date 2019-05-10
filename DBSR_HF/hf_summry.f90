!====================================================================== 
      Subroutine SUMMRY
!====================================================================== 
! ... the results of a calculation are summarized in log-file 
!----------------------------------------------------------------------
      Use DBS_grid
      Use DF_orbitals
      Use dbsr_hf

      Implicit none
      Integer :: i,j,m
      Real(8) :: en, r1,r2,rm1, EV,EK,EM, ratio, TA,VA,AM 
      Real(8), external :: quadr

      write(log,'(/84(''-''))')
      write(log,'(/a,a)') 'ATOM  ',ATOM 
      write(log,'(/a,a)') 'core: ',trim(adjustl(core)) 
      if(term.eq.'LS') then
       Call Clean_a(conf_LS)
       write(log,'(/a,a)') 'conf: ',trim(adjustl(conf_LS)) 
      elseif(term.eq.'jj') then
       Call Clean_a(configuration)
       write(log,'(/a,a)') 'conf: ',trim(adjustl(configuration)) 
      else
       write(log,'(/a,a)') 'conf: ',trim(adjustl(conf_AV)) 
      end if

! ... Compute Moments

      write(log,'(//2x,a,9x,a,10x,a,9x,a,6x,a,6x,a)') &
                'nl','e(nl)','dmp','ns','<r>','max_r'
      write(log,*)

      EK = 0.d0;  EM = 0.d0;  EV = 0.d0
      Do I = 1, NWF 
       Call TVM(kbs(i),p(1,1,i),p(1,2,i),TA,VA,AM)
       EK = EK + TA*qsum(i)
       EV = EV + VA*qsum(i)
       EM = EM + AM*qsum(i)
       r1  = QUADR(p(1,1,i),p(1,1,i),1) 
       r2  = QUADR(p(1,1,i),p(1,1,i),2)
       rm1 = QUADR(p(1,1,i),p(1,1,i),-1) 
       write(log,'(A5,f16.8,1PD13.2,I8,0P2f10.3)') &
        ebs(i), e(i,i), dpm(i), mbs(i), r1, t(mbs(i)+ks)  
      End do 
      
      write(log,'(/A,T25,1PD10.2)') 'Orbital convergence:', orb_diff
      write(log,'( A,T25,1PD10.2)') 'Energy  convergence:', scf_diff

      Call Conv_au (Z,ATW,au_cm,au_eV,0)
      if(debug.gt.0) then
       write(log,'(/a,f20.10)') 'a.u. --> cm-1 , au_cm = ', au_cm
       write(log,'(a,f20.10/)') 'a.u. --> eV   , au_eV = ', au_eV
      end if

      write(log,'(/a,T20,F24.15,a,T50,F25.15,a)')  &
        'Coulomb energy:', etotal,' au', etotal*au_eV, ' eV'

      r1 = EK + EV + EM
      EV = EV + E2body
      ratio = EV/EK

      if(debug.gt.0) then
       write(log,'( a,T20,F24.15)') 'Virial ratio:',ratio
       write(log,'(/a,T20,F24.15)') 'Kinetic:', EK
       write(log,'( a,T20,F24.15)') 'Potential:', EV
       write(log,'( a,T20,F24.15)') 'Mass:', EM
       write(log,'( a,T20,F24.15)') 'One-body:', E1body
       write(log,'( a,T20,F24.15)') 'Two-body:', E2body
      end if

      if(nconf.eq.1.and.mbreit.gt.0) then

       Call Breit_energy(E_breit)

       Call Self_energy(E_self)
       Call VP_energy(E_vacpol)                                  
       write(log,'(/a,T20,F24.15,T50,F24.15)') 'Breit:', E_breit,E_breit*au_eV
       write(log,'( a,T20,F24.15,T50,F24.15)') 'Self-energy:', E_self,E_self*au_eV
       write(log,'( a,T20,F24.15,T50,F24.15)') 'Vacuum polarization:', E_vacpol,E_vacpol*au_eV
       etotal = etotal+E_breit+E_self+E_vacpol
       write(log,'( a,T20,F24.15,T50,F24.15)') 'Sum:',etotal,etotal*au_eV

      else

       if(term.ne.'jj') Call Conf_energies
       if(term.eq.'jj') Call State_energies

      end if

      End Subroutine SUMMRY 




!====================================================================== 
      Subroutine Conf_energies
!====================================================================== 
! ... output of configuration energies 
! ... WARNING: weight and qsum are destroited here
!----------------------------------------------------------------------
      Use DBSR_hf
      Use DF_orbitals
      Use hf_energy, only: coef

      Implicit none
      Integer :: i,j,ic,jc,ii
      Real(8), allocatable :: EE_coul(:),EE_breit(:),EE_self(:),EE_vacpol(:), &
                              EE_tot(:)
      Integer, allocatable :: ipe(:)

      Allocate(EE_coul(nconf),EE_breit(nconf),EE_self(nconf), &
               EE_vacpol(nconf), EE_tot(nconf), ipe(nconf) )
      EE_coul=0.d0; EE_breit=0.d0; EE_self=0.d0; EE_vacpol=0.d0

! ... get all energies 

      Do ic=1,nconf
       qsum(:) = iqconf(ic,:)
       weight = 0.d0;  weight(ic) = 1.d0; coef = 0.d0
       Call av_energy_coef(nbf)

       Call ENERGY;    EE_coul(ic)=etotal

       if(mbreit.gt.0) then
        Call Breit_energy(EE_breit(ic))
        Call Self_energy(EE_self(ic))
        Call VP_energy(EE_vacpol(ic))
       end if

       EE_tot(ic) = EE_coul(ic) + EE_breit(ic) + EE_self(ic) + EE_vacpol(ic)

      End do ! over configurations

! ... sorting the energies:
 
      Call SORTR(nconf,EE_tot,ipe)

      write(log,'(/a/)') 'Configuration energies:'

! ... find biggest configuration:

      ii = 9
      Do ic=1,nconf
       i=0
       Do j=ncore+1,nwf
        if(iqconf(ic,j).eq.0) Cycle; i = i + 9
       End do
       if(i.gt.ii) ii=i
      End do

! ... output energies:

      Do jc=1,nconf;  ic = ipe(jc); etotal = EE_tot(ic)

       i=1
       Do j=ncore+1,nwf
        if(iqconf(ic,j).eq.0) Cycle
        write(configuration(i:),'(a5,a1,i2,a1)') ebs(j),'(',iqconf(ic,j),')'
        i = i + 9
       End do

       if(mbreit.eq.0) then 

        write(log,'(a,F24.15,a,3x,F24.15,a)') configuration(1:ii),etotal,' au', &
            etotal*au_eV,' eV'

       else  

        if(jc.eq.1) then
         SHELLJ=' '
         write(log,'(a,a16,3a13,a16)') &
         SHELLJ(1:ii),'Coulomb','Breit','Self','Vacpol','Sum'
        end if
 
        write(log,'(a,f16.4,3F13.4,f16.4)') configuration(1:ii),EE_coul(ic)*au_eV, &
         EE_breit(ic)*au_eV,EE_self(ic)*au_eV,EE_vacpol(ic)*au_eV,etotal*au_eV

       end if
      
      End do ! over configurations

      End Subroutine Conf_energies





!====================================================================
      Subroutine SORTR(n,S,IPT)
!--------------------------------------------------------------------
!     provide sorting pointer IPT for real array S
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in ) :: n
      Real(8), intent(in ) :: S(n) 
      Integer, intent(out) :: IPT(n) 
      Integer :: i,i1,j1,i2,j2
      Do i=1,n; IPT(i)=i;  End do
      Do i1=1,n-1;    j1=IPT(i1)
       Do i2=i1+1,n;  j2=IPT(i2)
        if(S(j1).gt.S(j2)) then; IPT(i2)=j1; j1=j2; end if
       End do
       IPT(i1)=j1
      End do
      End Subroutine SORTR

