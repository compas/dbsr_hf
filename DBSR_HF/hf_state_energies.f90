!====================================================================== 
      Subroutine State_energies
!====================================================================== 
! ... calculation and output of the state energies obtained in CI procedure
!----------------------------------------------------------------------
      Use DBSR_hf
      Use DF_orbitals
      Use hf_energy
      Use rk4_data
      Use DBS_grid, only: ns,ks                      

      Implicit none

      Integer, allocatable :: noc(:)
      Integer, allocatable :: nnc(:,:), knc(:,:),lnc(:,:),jnc(:,:),iqc(:,:),inc(:,:),ipc(:,:)
      Integer, allocatable :: Jshellc(:,:),Vshellc(:,:),Jintrac(:,:)
      Integer :: njbl
      Integer, allocatable :: JT(:), JP(:), JTc(:), ipe(:), index(:)
      Integer, parameter :: mlab = 200
      Character(mlab), allocatable :: Labels(:)
      Integer :: ilabel=0

      Integer, parameter :: mcoef = 1000
      Integer :: ncoef
      Real(8) :: coefs(mcoef)
      Integer :: icoefs(5,mcoef)

      Integer :: i,j, ic,jc, ii, jj, ij, ip, i1,i2,i3,i4, k,v,met, k1,k2,k3,k4, &
                 iset, info, nsol, ib, nc
      Integer, external :: Ifind_orb, Jdef_ncfg

      Real(8) :: S(8), C, E1
      Real(8), allocatable :: HC(:,:), HB(:,:), H(:,:), H_self(:), H_vacpol(:)
      Real(8), allocatable :: eval(:), EE_coul(:),EE_breit(:),EE_self(:),EE_vacpol(:), &
                              EE_tot(:)
      Real(8), external :: SMU, VP_dhl, hf_rk, sk_ppqq

      AF_cfg = trim(name)//BF_cfg
      Call Read_aarg('c',AF_cfg)
      inquire(file=AF_cfg,number=i)
      if(i.gt.0) Close(i)
      open(nuc,file=AF_cfg)
      rewind(nuc)      

      nconf = Jdef_ncfg(nuc)

      if(nconf.le.0) Return

! ... define all states from c-file:

      Allocate(noc(nconf), nnc(msh,nconf), knc(msh,nconf),lnc(msh,nconf), &
               jnc(msh,nconf), iqc(msh,nconf),inc(msh,nconf),ipc(msh,nconf), &
               Jshellc(msh,nconf),Vshellc(msh,nconf),Jintrac(msh,nconf), &
               JT(nconf), JP(nconf), JTc(0:nconf), Labels(nconf) )

      AF_cfg = trim(name)//BF_cfg
      Call Read_aarg('c',AF_cfg)
      inquire(file=AF_cfg,number=i)
      if(i.gt.0) Close(i)
      open(nuc,file=AF_cfg)
      rewind(nuc)      

      ic = 0;  njbl = 0; JTc = 0
      Do 
       read(nuc,'(a)') CONFIG
       if(CONFIG(6:6).ne.'(') Cycle
       read(nuc,'(a)') SHELLJ
       read(nuc,'(5x,a)') INTRAJ
       ic = ic + 1
       Call Decode_cjj(CONFIG,SHELLJ,INTRAJ, &
         noc(ic),nnc(1,ic),knc(1,ic),lnc(1,ic),jnc(1,ic),iqc(1,ic),inc(1,ic), &
         Jshellc(1,ic),Vshellc(1,ic),Jintrac(1,ic) )
       no = noc(ic); ip=0

       Do i=1,no
        ip = ip + lnc(i,ic)*iqc(i,ic)
       End do
       JT(ic) = Jintrac(noc(ic),ic)
       JP(ic) = (-1) ** ip

       if(ic.eq.1) then
        njbl = njbl + 1; JTc(njbl)=ic
       else
        if(JT(ic).eq.JT(ic-1).and.JP(ic).eq.JP(ic-1)) then
         JTc(njbl) = ic
        else
         njbl = njbl + 1; JTc(njbl)=ic
        end if
       end if 
       iset = 0
       Call labelc_jj (mlab,Labels(ic),iset, &
            no,nnc(1,ic),lnc(1,ic),knc(1,ic),jnc(1,ic),inc(1,ic),iqc(1,ic), &
            Jshellc(1,ic),Vshellc(1,ic),Jintrac(1,ic) )
       ii=LEN_TRIM(Labels(ic)); if(ii.gt.ilabel) ilabel=ii 

       ! ... add core:  
       if(ncore.gt.0) then
        noc(ic) = no + ncore
        nnc(ncore+1:ncore+no,ic) =  nnc(1:no,ic)
        lnc(ncore+1:ncore+no,ic) =  lnc(1:no,ic)
        knc(ncore+1:ncore+no,ic) =  knc(1:no,ic)
        jnc(ncore+1:ncore+no,ic) =  jnc(1:no,ic)
        inc(ncore+1:ncore+no,ic) =  inc(1:no,ic)
        iqc(ncore+1:ncore+no,ic) =  iqc(1:no,ic)
        nnc(1:ncore,ic) = n_core(1:ncore)
        lnc(1:ncore,ic) = l_core(1:ncore)
        knc(1:ncore,ic) = k_core(1:ncore)
        jnc(1:ncore,ic) = j_core(1:ncore)
        inc(1:ncore,ic) = 0
        iqc(1:ncore,ic) = j_core(1:ncore)+1
        Jshellc(ncore+1:ncore+no,ic) =  Jshellc(1:no,ic)
        Vshellc(ncore+1:ncore+no,ic) =  Vshellc(1:no,ic)
        Jintrac(ncore+1:ncore+no,ic) =  Jintrac(1:no,ic)
        Jshellc(1:ncore,ic)=0 
        Vshellc(1:ncore,ic)=0 
        Jintrac(1:ncore,ic)=0
       end if

       Do i=1,noc(ic)
        ipc(i,ic) = Ifind_orb(nnc(i,ic),knc(i,ic),inc(i,ic))
       End do


       if(ic.eq.nconf) Exit
      End do
      close(nuc)

! ... define all angular coefficients and store them in rk4_data module:

      Do ic=1,nconf
      Do jc=1,ic
       if(JT(ic).ne.JT(jc)) Cycle
       if(JP(ic).ne.JP(jc)) Cycle

       Call coef_2conf_jj(noc(ic),nnc(1,ic),lnc(1,ic),jnc(1,ic),iqc(1,ic), &
                          Jshellc(1,ic),Vshellc(1,ic),Jintrac(1,ic),  &
                          noc(jc),nnc(1,jc),lnc(1,jc),jnc(1,jc),iqc(1,jc), &
                          Jshellc(1,jc),Vshellc(1,jc),Jintrac(1,jc),  &
                          mcoef,ncoef,icoefs,coefs)
       
       if(ncoef.eq.0.d0) Cycle;  ij = ic*ibi+jc
       
       Do i = 1,ncoef; if(abs(coefs(i)).lt.eps_c) Cycle; met=2

        k = icoefs(1,i)
        if(k.lt.0) then; met=1; k=0; end if         
        i1=icoefs(2,i); i1=ipc(i1,ic); k1=kbs(i1)
        i2=icoefs(3,i); i2=ipc(i2,ic); k2=kbs(i2)
        i3=icoefs(4,i); i3=ipc(i3,jc); k3=kbs(i3)
        i4=icoefs(5,i); i4=ipc(i4,jc); k4=kbs(i4)

        Call Add_rk4_data(met*ibi+k,i1*ibi+i2,i3*ibi+i4,ij,coefs(i))

        if(met.ne.2.or.mbreit.eq.0) Cycle 

        Do v = k-1,k+1
         if(SMU(k1,k2,k3,k4,k,v,S).eq.0.d0) Cycle
         S = S * coefs(i);  met = 3*ibi+v
         Call Add_rk4_data(met,i1*ibi+i2,i3*ibi+i4,ij,S(1)) 
         Call Add_rk4_data(met,i2*ibi+i1,i4*ibi+i3,ij,S(2)) 
         Call Add_rk4_data(met,i3*ibi+i4,i1*ibi+i2,ij,S(3)) 
         Call Add_rk4_data(met,i4*ibi+i3,i2*ibi+i1,ij,S(4)) 
         Call Add_rk4_data(met,i1*ibi+i4,i3*ibi+i2,ij,S(5)) 
         Call Add_rk4_data(met,i4*ibi+i1,i2*ibi+i3,ij,S(6)) 
         Call Add_rk4_data(met,i3*ibi+i2,i1*ibi+i4,ij,S(7)) 
         Call Add_rk4_data(met,i2*ibi+i3,i4*ibi+i1,ij,S(8)) 
        End do                 
       End do         ! over coeffs, i
      End do; End do ! over states ic,jc

! ... define Hamiltonia matrix:

      Allocate(HC(nconf,nconf), HB(nconf,nconf), H_self(nconf), H_vacpol(nconf))
      HC = 0.d0; HB = 0.d0; H_self = 0.d0; H_vacpol = 0.d0

      ! one-electron integrals
 
      Do i=1,nrk;  met=kr1(i)/ibi; if(met.ne.1) Cycle
       i1 = kr2(i)/ibi;  i2 = kr3(i)/ibi  
       ic = kr4(i)/ibi;  jc = mod(kr4(i),ibi)  
       if(i1.eq.i2) then
        HC(ic,jc) = HC(ic,jc) + crk(i)*ValueI(i1)
       else
        HC(ic,jc) = HC(ic,jc) + crk(i)*Vp_dhl (kbs(i1),p(1,1,i1),kbs(i2),p(1,1,i2))
       end if
      End do

      ! two-electron Coulomb integrals: B-spline integrals should be allocated !!!

      Do i=1,nrk;  met=kr1(i)/ibi; if(met.ne.2) Cycle; k = mod(kr1(i),ibi)
       i1 = kr2(i)/ibi; i2 = mod(kr2(i),ibi) 
       i3 = kr3(i)/ibi; i4 = mod(kr3(i),ibi) 
       ic = kr4(i)/ibi; jc = mod(kr4(i),ibi)  

       if(mod(k+lbs(i1)+lbs(i3),2).ne.0) Cycle
       if(mod(k+lbs(i2)+lbs(i4),2).ne.0) Cycle

       if(ic.eq.jc.and.i1.eq.i3.and.i2.eq.i4) then            ! Fk integrals
         i3 = max(i1,i2); i4 = min(i1,i2)
         HC(ic,jc) = HC(ic,jc) + crk(i)*Value(i3,i4,k)
       elseif(ic.eq.jc.and.i1.eq.i4.and.i2.eq.i3) then        ! Gk integrals
         i3 = min(i1,i2); i4 = max(i1,i2)
         HC(ic,jc) = HC(ic,jc) + crk(i)*Value(i3,i4,k)
       else                                                   ! Rk integrals
         HC(ic,jc) = HC(ic,jc) + crk(i)*hf_rk(i1,i2,i3,i4,k)
       end if   
      End do

      ! two-electron Breit integrals: 

      if(mbreit.gt.0) then

       C = 0.d0
       Call Alloc_DBS_integrals(ns,2*ks-1,0,kmax+1,1)
       Do i=1,nrk;  met=kr1(i)/ibi; if(met.ne.3) Cycle; k = mod(kr1(i),ibi)
        ic = kr4(i)/ibi;  jc = kr4(i)/ibi  
        if(C.eq.0.d0) then
         i1 = kr2(i)/ibi; i2 = mod(kr2(i),ibi) 
         i3 = kr3(i)/ibi; i4 = mod(kr3(i),ibi) 
         C =  sk_ppqq(i1,i2,i3,i4,k)
        elseif(kr1(i).ne.kr1(i-1).or.kr2(i).ne.kr2(i-1).or. &
               kr3(i).ne.kr3(i-1))  then
         i1 = kr2(i)/ibi; i2 = mod(kr2(i),ibi) 
         i3 = kr3(i)/ibi; i4 = mod(kr3(i),ibi) 
         C =  sk_ppqq(i1,i2,i3,i4,k)
        end if
        HB(ic,jc) = HB(ic,jc) + crk(i)*C
       End do

       Do ic=1,nconf
        Call Self_energy(C);   H_self(ic) = C
        Call VP_energy(C);     H_vacpol(ic) = C
       End do

      end if

      Do i=1,nconf-1; Do j=i+1,nconf
        HC(i,j) = HC(j,i);  HB(i,j) = HB(i,j)   
      End do; End do                                                             

!----------------------------------------------------------------------
! ... diagonalizing the Hamiltonian matrix:

      Allocate(EE_coul(nconf),EE_breit(nconf),EE_self(nconf), &
               EE_vacpol(nconf), EE_tot(nconf) )
      EE_coul=0.d0; EE_breit=0.d0; EE_self=0.d0; EE_vacpol=0.d0

      Allocate(H(nconf,nconf),eval(nconf),ipe(nconf),index(nconf))

      AF = trim(name)//'.j'
      open(nuc,file=AF)
      write(nuc,'(a,i4)') 'ncfg = ',nconf
      write(nuc,'(a,i4)') 'nsol = ',nconf
      write(nuc,'(/a)') 'Solutions:'

      nsol = 0
      Do ib = 1,njbl
       i1 = JTc(ib-1)+1; i2=JTc(ib); nc = i2-i1+1

       H(1:nc,1:nc) = HC(i1:i2,i1:i2) + HB(i1:i2,i1:i2)
       Do i = 1,nc
        H(i,i) = H(i,i) + H_self(i+i1-1)
        H(i,i) = H(i,i) + H_vacpol(i+i1-1)
       End do

       Call LAP_DSYEV ('V','L',nc,nconf,H,eval,info)     

       ipe = 0
       Do k = 1,nc
        C=0.d0; ii=0                 ! find the leading component
        Do i=1,nc
         if(ipe(i).ne.0) Cycle
         if(abs(H(i,k)).le.C) Cycle
         C=abs(H(i,k)); ii=i
        End do
        ipe(ii)=k
        nsol = nsol + 1
        write(nuc,'(i8,2x,a)')  nsol, TRIM(Labels(ii+i1-1))
        write(nuc,'(f16.8,3i8,f16.4)')  eval(k), JT(ib), i1,i2, EVAL(k)*au_eV
        write(nuc,'(6f20.15)') H(1:nc,k)

        EE_tot(nsol) =  eval(k);  index(nsol)=ii+i1-1

        if(mbreit.eq.0) Cycle

         Do i = 1,nc; ii=i+i1-1
          EE_self(nsol) = EE_self(nsol) + H(i,k)*H(i,k)*H_self(ii)
          EE_vacpol(nsol) = EE_vacpol(nsol) + H(i,k)*H(i,k)*H_vacpol(ii)
         Do j = 1,nc; jj=j+i1-1
          EE_coul(nsol) = EE_coul(nsol) + H(i,k)*HC(ii,jj)*H(j,k)
          EE_breit(nsol) = EE_breit(nsol) +  H(i,k)*HB(ii,jj)*H(j,k)
         End do; End do

       End do

      End do  ! over J-blocks, ib

! ... output of state energies:


      Call SORTR(nconf,EE_tot,ipe)

      write(log,'(/a/)') 'Atomic state energies:'

      E1 = EE_tot(ipe(1))

      Do jc=1,nconf;  ic = ipe(jc); etotal = EE_tot(ic)
                      configuration = labels(index(ic)); ii = ilabel

       if(mbreit.eq.0) then 

        write(log,'(a,F25.15,a,3x,F14.5,a)') configuration(1:ii),etotal,' au', &
            (etotal-E1)*au_eV,' eV'

       else  
        SHELLJ = ' '
        if(jc.eq.1)  write(log,'(a,a16,3a13,a16)') &
                     SHELLJ(1:ii),'Coulomb','Breit','Self','Vacpol','Total'
 
        write(log,'(a,f16.4,3F13.4,f16.5)') configuration(1:ii),EE_coul(ic)*au_eV, &
         EE_breit(ic)*au_eV,EE_self(ic)*au_eV,EE_vacpol(ic)*au_eV,(etotal-E1)*au_eV
       end if
      
      End do ! over configurations, jc
 
      End Subroutine State_energies 



!======================================================================
      Subroutine labelc_jj (mlab,Lab,iset,no,nn,ln,kn,jn,in,iq, &
                                          Jshell,Vshell,Jintra)
!======================================================================
!     Lab  ->   packed configuration notation
!     iset = 0 -> all set numbers forced to 0
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: mlab,iset
      Integer, intent(in) :: no,nn(no),ln(no),kn(no),jn(no),in(no),iq(no), &
                             Jshell(no),Vshell(no),Jintra(no)
      Character(mlab), intent(out) ::  Lab
      Character(5) :: EL
      Character(5), external :: ELi
      Integer :: i,ii,j,jj,k,kk, JT,JV,JW,JQ
      Integer, external :: Jterm, Jparity
      Lab = ' '

      k=1
      Do i=1,no; if(iq(i).lt.1) Cycle

       kk = jn(i)+1
       if(k.eq.1.and.no-i.ge.2.and.iq(i).eq.kk) Cycle

       ! ... orbital

       ii=in(i); if(iset.eq.0) ii=0      
       EL=ELi(nn(i),kn(i),ii)
       Do j = 1,5
        if(EL(j:j).eq.' ') Cycle
        Lab(k:k)=EL(j:j); k=k+1
       End do

       ! ... number of electrons

       if(iq(i).gt.0.and.iq(i).le.9) then
        write(Lab(k:k),'(a)') '('; k=k+1 
        write(Lab(k:k),'(i1)') iq(i); k=k+1
        write(Lab(k:k),'(a)') ')'; k=k+1 
       elseif(iq(i).gt.9) then
        write(Lab(k:k),'(a)') '('; k=k+1 
        write(Lab(k:k+1),'(i2)') iq(i); k=k+2
        write(Lab(k:k),'(a)') ')'; k=k+1 
       end if
       Lab(k:k)='.'; k=k+1

       ! ... shell term
      
       if(Jterm (jn(i),iq(i),-1,JT,JV,JW,JQ).gt.1) then
        k=k-1;  Lab(k:k)='_'; k=k+1

        kk = mod(Jshell(i),2)
        if(kk.eq.0) then
         jj = Jshell(i)/2
         if(jj.le.9) then
          write(Lab(k:k),'(i1)') jj; k=k+1
         else
          write(Lab(k:k+1),'(i2)') jj; k=k+2
         end if
        else
         jj = Jshell(i)
         if(jj.le.9) then
          write(Lab(k:k),'(i1)') jj; k=k+1
         else
          write(Lab(k:k+1),'(i2)') jj; k=k+2
         end if
          write(Lab(k:k+1),'(a2)') '/2'; k=k+2
        end if

        Lab(k:k)='.'; k=k+1

       end if

       ! ... intermediate term

       if(no.gt.1.and.i.eq.1) Cycle
       kk = 0 
       if(i.gt.1) then
        if(iabs(Jshell(i)-Jintra(i-1)).lt.Jshell(i)+Jintra(i-1)) kk=1
       end if

       if(kk.eq.1.or.i.eq.no) then
        kk = mod(Jintra(i),2)
        if(kk.eq.0) then
         jj = Jintra(i)/2
         if(jj.le.9) then
          write(Lab(k:k),'(i1)') jj; k=k+1
         else
          write(Lab(k:k+1),'(i2)') jj; k=k+2
         end if
        else
         jj = Jintra(i)
         if(jj.le.9) then
          write(Lab(k:k),'(i1)') jj; k=k+1
         else
          write(Lab(k:k+1),'(i2)') jj; k=k+2
         end if
          write(Lab(k:k+1),'(a2)') '/2'; k=k+2
        end if
        if(i.ne.no) then; Lab(k:k)='_'; k=k+1; end if
       end if

      if(k-1.gt.mlab) write(*,*) 'Label_jj too long'

      End do

      i = Jparity(no,ln,iq); if(i.eq.-1) LAB(k:k)='*'

      End Subroutine labelc_jj

!======================================================================
      Subroutine LAP_DSYEV(job,UPLO,n,m,A,eval,info)
!======================================================================
!     Call LAPACK procedure DSYEV to obtain the eigenvalues (eval) and 
!     eigenvectors (A) for Real symmetric matrix A(n,n)
!     job = 'V' or 'N' - compute or not the eigenvectors
!     UPLO='U'  or 'L' -  used upper or lower triangle matrix
!---------------------------------------------------------------------
      Implicit none
      Character(1), intent(in) :: job,UPLO
      Integer, intent(in) :: n,m
      Integer, intent(out) :: info
      Real(8) :: A(m,m)
      Real(8) :: eval(m)
      Real(8) :: work(3*n-1)
      Integer :: lwork

      lwork = 3*n-1

      Call DSYEV(job,UPLO,n,A,m,eval,WORK,LWORK,INFO )

      if(INFO.ne.0) write(*,*) 'LAP_DSYEV: DSYEV(lapack) gives INFO = ',INFO

      End Subroutine LAP_DSYEV


