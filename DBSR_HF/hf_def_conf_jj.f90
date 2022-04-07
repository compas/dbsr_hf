!======================================================================
      Subroutine Def_conf_jj
!======================================================================
!     defines atomic states to be included into HF optimization
!     atomic states are geven in standard c-file
!----------------------------------------------------------------------
      Use dbsr_hf
      Use df_orbitals

      Implicit none
      Integer :: i,j,ic,io,is,k,iset
      Integer, external :: Jdef_ncfg, Ifind_position, Ifind_orb, &
                           Iadd_cfg_jj
      Character(5), external :: Eli

! ... check c-file:

      AF_cfg = trim(name)//BF_cfg
      Call Read_aarg('c',AF_cfg)
      if(Icheck_file(AF_cfg).eq.0) then
       write(log,'(/a/)') 'Stop in Def_conf_jj:'
       write(log,'( a/)') 'c-file is supposed to be provided for term=jj option'
       write(scr,'(/a/)') 'Stop in Def_conf_jj:'
       write(scr,'( a/)') 'c-file is supposed to be provided for term=jj option'
       Stop ' '
      end if

      open(nuc,file=AF_cfg)

! ... define number of atomic states:

      nconf = Jdef_ncfg(nuc)
      if(nconf.eq.0) then
       write(log,'(/a/)') 'ncfg=0 in c-file: nothing to do'
       write(scr,'(/a/)') 'ncfg=0 in c-file: nothing to do'
       Stop 'nconf=0 in c-file: nothing to do'
      end if

! ... define core:

      i=Ifind_position(nuc,'Core subshells:')
      read(nuc,'(/a)') core
      Call Def_core

! ... define one-electron orbitals:

      i=Ifind_position(nuc,'Peel subshells:')
      read(nuc,'(/a)') config
      i = LEN_TRIM(config)
      nwf = i/5; if(mod(i,5).ne.0) nwf=nwf+1; nwf = nwf + ncore

      Call alloc_DF_orbitals(nwf)

      Do i=1,ncore
       nbs(i)  = n_core(i)
       lbs(i)  = l_core(i)
       kbs(i)  = k_core(i)
       jbs(i)  = j_core(i)
       ebs(i)  = e_core(i)
       qsum(i) = jbs(i)+1
       ibs(i)  = 0
       clsd(i) = .TRUE.
      End do

      k = 0
      Do i=ncore+1,nwf
       read(config(k+1:k+5),'(a)') ebs(i)
       Call EL_NLJK(ebs(i),nbs(i),kbs(i),lbs(i),jbs(i),ibs(i))
       k = k + 5
      End do

      Do i=1,nwf; ebs(i)=ELi(nbs(i),kbs(i),0); End do

! ... atomic states paramters:

      Allocate(weight(nconf));  weight = 0.d0
      Call Read_ipar(inp,'eal',eal)
      Call Read_iarg('eal',eal)
      if(eal.eq.1) weight = 1.d0

      Allocate(iqconf(nconf,nwf)); iqconf = 0
      Do i=1,ncore;  iqconf(:,i) = qsum(i); End do

      i=Ifind_position(nuc,'CSF(s):')

      ic = 0
      Do
       read(nuc,'(a)') CONFIG

       if(CONFIG(6:6).ne.'(') Cycle
       read(nuc,'(a)') SHELLJ
       read(nuc,'(5x,a)') INTRAJ
       Call Decode_cjj(CONFIG,SHELLJ,INTRAJ,no,nn,kn,ln,jn,iq,in,&
                                            Jshell,Vshell,Jintra)
       ic = ic + 1

       Do is=1,no
        io=Ifind_orb(nn(is),kn(is),0)
        iqconf(ic,io) = iq(is)
       End do

       if(eal.eq.5) weight(ic) = JINTRA(no)+1
       if(eal.eq.9.and.LEN_TRIM(config(no*9+1:)).ne.0) &
        read(config(no*9+1:),*) weight(ic)

       if(ic.eq.nconf) Exit

      End do

      if(SUM(weight).eq.0.d0) weight = 1.d0
      weight = weight / SUM(weight)

! ... redefine qsum:

      qsum = 0.d0
      Do io=1,nwf
       Do ic=1,nconf
        qsum(io)=qsum(io) + iqconf(ic,io)*weight(ic)
       End do
      End do

      End Subroutine Def_conf_jj
