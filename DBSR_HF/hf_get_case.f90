!======================================================================
      Subroutine get_case
!======================================================================
!     get information about the problem to be solved
!     and how the spline methods to be applied
!----------------------------------------------------------------------
      Use dbsr_hf
      Use df_orbitals

      Implicit none
      Integer :: i,j

! ... name of case:

      Call Check_name

      AF_log = trim(name)//BF_log
      Call Read_aarg('log',AF_log)
      Open(log,file=AF_log)

! ... header:

      write(scr,'(/21x,a/21x,a/21x,a/)') &
      '=======================================',&
      ' B - S P L I N E  D I R A C - F O C K  ',&
      '======================================='
      write(scr,'(a,a)') 'name:  ',name

      write(log,'(/a/a/a)') &
      '========================================',&
      ' B - S P L I N E  D I R A C - F O C K   ',&
      '========================================'
      write(log,'(/a,a)') 'name:   ',name

! ... read input-data:

      Call Read_data

! ... initialize Lagrange multipliers:

      e = 0.d0
      Do i = 1,nwf
       Do j = 1, i-1
        if (kbs(i).ne.kbs(j)) Cycle
        e(i,j) = 1.d-5;  e(j,i) = 1.d-12
       End do
      End do

      End Subroutine get_case


!======================================================================
     Subroutine Check_name
!======================================================================
! .. read name of the case and check if name.inp exist;
! .. if exist - return, otherwise check if name is a true atom symbol
! .. and prepare file "name.inp";
! .. atomic number, an=.., can be supplied instead atomic symbol
!----------------------------------------------------------------------
     Use dbsr_hf
     Implicit none
     Real(8) :: a1,r1

     Call Read_name (name)
     if(len_trim(name).ne.0) Return

     Call Read_aarg('atom',atom)
     Call Read_iarg('an',an)
     Call Read_aarg('ion',ion)
     Call Read_iarg('ai',ai)

     if(len_trim(atom).gt.0) then
      Call Def_atom(an,atom,a1,r1,conf_AV,conf_LS)
     elseif(an.gt.0) then
      Call Def_atom(an,atom,a1,r1,conf_AV,conf_LS)
     end if

     if(len_trim(ion).gt.0) then
      Call Def_atom(ai,ion,a1,r1,conf_AV,conf_LS)
     elseif(ai.gt.0) then
      Call Def_atom(ai,ion,a1,r1,conf_AV,conf_LS)
     else
      ai = an
      Call Def_atom(ai,ion,a1,r1,conf_AV,conf_LS)
     end if

     name = atom
     if(len_trim(ion).gt.0.and.an.gt.0.and.atom.ne.ion) &
        write(name,'(a,a,i3.3)')  trim(ion),'_',an

     if(len_trim(name).eq.0) &
     Stop ' Please, indicate name of the case, e.g., atom symbol:  dbsr_hf Ne'

     End Subroutine Check_name


!======================================================================
      Subroutine Read_data
!======================================================================
!     This routine obtains information about the problem to be solved
!     and how the spline methods are to be applied
!----------------------------------------------------------------------
      Use dbsr_hf
      Use DBS_grid
      Use df_orbitals

      Implicit none
      Integer :: i,j, bn
      Real(8) :: a1,r1

      AF_dat = trim(name)//BF_dat
      Call Read_aarg('dat',AF_dat)
      Open(inp,file=AF_dat)

      Call Read_apar(inp,'atom',atom)
      Call Read_apar(inp,'ion',ion)
      Call Read_rpar(inp,'z',z)
      Call Read_rpar(inp,'atw',atw )
      Call Read_rpar(inp,'rms',rms )
      Call Read_string(inp,'core',core)
      Call Read_apar(inp,'conf',configuration)
      Call Read_apar(inp,'term',term)
      Call Read_string(inp,'varied',anit)

      Call Read_rpar(inp,'scf_tol',scf_tol)
      Call Read_rpar(inp,'orb_tol',orb_tol)
      Call Read_rpar(inp,'end_tol',end_tol)
      Call Read_ipar(inp,'max_it' ,max_it )
      Call Read_ipar(inp,'newton' ,newton )
      Call Read_ipar(inp,'rotate' ,rotate )
      Call Read_ipar(inp,'debug'  ,debug  )
      Call Read_ipar(inp,'ilzero' ,ilzero )
      Call Read_ipar(inp,'ibzero' ,ibzero )
      Call Read_ipar(inp,'mbreit' ,mbreit )
      Call Read_ipar(inp,'out_nl' ,out_nl )
      Call Read_ipar(inp,'out_w'  ,out_w  )
      Call Read_ipar(inp,'out_plot',out_plot)
      Call Read_ipar(inp,'mode_SE',mode_SE)
      Call Read_ipar(inp,'mode_VP',mode_VP)

      Call Read_apar(inp,'knot',knot)

!-----------------------------------------------------------------------------------

      Call Read_aarg('atom',atom)
      Call Read_iarg('an',an)
      Call Read_rarg('z',z)
      Call Read_aarg('ion',ion)
      Call Read_iarg('ai',ai)

      Call Read_rarg('z',z)
      Call Read_rarg('atw',atw)
      Call Read_rarg('rms',rms)
      Call Read_aarg('core',core)
      Call Read_aarg('conf',configuration)
      Call Read_aarg('term',term)
      Call Read_aarg('varied',anit)

      Call Read_rarg('scf_tol',scf_tol)
      Call Read_rarg('orb_tol',orb_tol)
      Call Read_rarg('end_tol',end_tol)
      Call Read_iarg('max_it' ,max_it )
      Call Read_iarg('newton' ,newton )
      Call Read_iarg('rotate' ,rotate )
      Call Read_iarg('debug'  ,debug  )
      Call Read_iarg('ilzero' ,ilzero )
      Call Read_iarg('ibzero' ,ibzero )
      Call Read_iarg('mbreit' ,mbreit )
      Call Read_iarg('out_nl' ,out_nl )
      Call Read_iarg('out_w'  ,out_w  )
      Call Read_iarg('out_plot',out_plot)
      Call Read_iarg('mode_SE',mode_SE)
      Call Read_iarg('mode_VP',mode_VP)

      Call Read_aarg('knot',knot)

!-----------------------------------------------------------------------------------

      if(len_trim(atom).eq.0) then
       if(an.gt.0) then
        Call Def_atom(an,atom,a1,r1,conf_AV,conf_LS)
       elseif(Z.gt.0.d0) then
        an = NINT(Z)
        Call Def_atom(an,atom,a1,r1,conf_AV,conf_LS)
       else
        atom = name
       end if
      end if

      Call Def_atom(an,atom,a1,r1,conf_AV,conf_LS)
      if(len_trim(atom).eq.0) &
        Stop 'Can not identify atom: indicate either atom=.. or an=..'

      if(Z.eq.0.d0) Z = an

      if(ai.eq.0.and.len_trim(ion).eq.0) then
       ion=atom; ai=an
      else
       Call Def_atom(ai,ion,a1,r1,conf_AV,conf_LS)
      end if

      if(atw.eq.0.d0)  Call Def_atom(an,atom,atw,r1,conf_AV,conf_LS)
      if(rms.eq.0.d0)  Call Def_atom(an,atom,a1,rms,conf_AV,conf_LS)

      if(len_trim(core).eq.0) &
         Call Def_atom(ai,ion,a1,r1,core,conf_LS)
      if(len_trim(configuration).eq.0) &
         Call Def_atom(ai,ion,a1,r1,conf_AV,configuration)

!-----------------------------------------------------------------------------------

      Call Def_core

      if(term.eq.'jj') then
       Call Def_conf_jj
      elseif(term.eq.'LS') then
       Call Def_conf_LS
       Call Def_conf
      else
       Call Def_conf
      end if

! ... define generalized configuration:

      i = 1
      Do j=ncore+1,nwf
       write(conf_AV(i:),'(a5,a1,f4.1,a1)') &
        ebs(j),'(',qsum(j),')'
       i = i + 11
      End do
      Call Clean_a(conf_AV)
      Call Clean_a(configuration)

      Call Def_nit(anit)

      write(log,'( /a,a)') 'ATOM    ',atom
      write(log,'(  a,a)') 'TERM    ',term
      write(log,'(a,f5.0)')'Z    ',z
      if(LEN_TRIM(core).ne.0) &
      write(log,'(/ a,a)') 'core:   ',trim(adjustl(core))
      write(log,'(  a,a)') 'conf:   ',trim(adjustl(conf_AV))
      write(log,'(/ a,a)') 'orbitals to optimize:  ',anit
      if(LEN_TRIM(core).ne.0) &
      write(scr,'(/a,a)') 'core:  ',trim(adjustl(core))
      if(len_trim(conf_AV).gt.0) &
      write(scr,'(a,a/)') 'conf:  ', trim(adjustl(conf_AV))

      if(nconf.gt.1) then
      if(term.eq.'LS')  write(log,'(/a,i3,a)') 'nconf =',nconf,&
	   '  - number of configurations to optimize, see conf-file'
      if(term.eq.'jj')  write(log,'(/a,i3,a)') 'nconf =',nconf,&
	   '  - atomic states to optimize, see c-file'

      if(eal.eq.1)  write(log,'(a,i3,a)') 'eal   =',eal,&
	   '  - weights of states are equal'
      if(eal.eq.5)  write(log,'(a,i3,a)') 'eal   =',eal,&
	   '  - weights of states are statistical'
      if(eal.eq.9)  write(log,'(a,i3,a)') 'eal   =',eal,&
	   '  - weights of states are chosen by user'
      end if

      write(log,'(//a/)') 'Running parameters:'
      write(log,'(a,1PE9.2,T25,a)') 'scf_tol = ',scf_tol, &
                '- energy convergence tolerance'
      write(log,'(a,1PE9.2,T25,a)') 'orb_tol = ',orb_tol, &
                '- orbital convergence tolerance'
      write(log,'(a,1PE9.2,T25,a)') 'end_tol = ',end_tol, &
                '- orbital tail cut-off'
      write(log,'(a,i3,T25,a)') 'max_it  = ',max_it, &
                '- max. number of iterations'
      write(log,'(a)')
      if(rotate.gt.0) &
      write(log,'(a,i2,T25,a)') 'rotate  = ',rotate, &
                '- use rotations of orbitals'

      Call Write_inp

      End Subroutine Read_data



!======================================================================
      Subroutine Def_core
!======================================================================
!     This routine defines the closed shells which will be considered as
!     a core (from the string "core").
!     DBSR_HF has no special treatment of core orbitals, however,
!     it is used for more compact notation of atomic configuration
!     and for consistency with other GRASP and DBSR programs.
!     INPUT:  string  "core"
!     OUTPUT: ncore, e_core,n_core,k_core,l_core(i),j_core
!                    all -> arrays (1:ncore)
!----------------------------------------------------------------------
      Use dbsr_hf
      Use atoms_par

      Implicit none
      Integer :: i,i1,i2,ii

      ncore = 0
      ii = LEN_TRIM (core)
      if(ii.eq.0) Return
      i1 = INDEX (core,'[')
      i2 = INDEX (core,']')
      if(i1*i2 /= 0) then    ! append default core (as in the periodic table)
       Select case(core(i1+1:i2-1))
        case ('He','1s'); core = trim(He)//' '//core(i2+1:ii)
        case ('Ne','2p'); core = trim(Ne)//' '//core(i2+1:ii)
        case ('Ar','3p'); core = trim(Ar)//' '//core(i2+1:ii)
        case ('Kr','4p'); core = trim(Kr)//' '//core(i2+1:ii)
        case ('Xe','5p'); core = trim(Xe)//' '//core(i2+1:ii)
        case ('Hg','6s'); core = trim(Hg)//' '//core(i2+1:ii)
        case ('Rn','6p'); core = trim(Rn)//' '//core(i2+1:ii)
        case default
         write(*,*) 'Default core must be one of ', &
                    'Ne, Be, Ne, Mg, Ar, Zn, Kr, Cd, Xe, Hg, Rn '
         write(*,*) 'Example:  core=[Ne]  or core=[2p]'
         Stop 'Stop in def_core routine: unable to define core orbitals'
       End Select
      end if

      Do i = 1,mcore
       read(core,*,IOSTAT=ii) e_core(1:i)
       if (ii /= 0) Exit
      End do
      ncore = i-1
      Do i = 1,ncore
       e_core(i) = adjustl(e_core(i))
       Call EL_NLJK(e_core(i),n_core(i),k_core(i),l_core(i),j_core(i),ii)
      End do
      core = ' '
      write(core,'(50a5)') (e_core(i),i=1,ncore)

      End Subroutine Def_core


!======================================================================
      Subroutine Def_conf
!======================================================================
!     Define configuration from the file  'name.conf' or from
!     the input:  conf = ...
!----------------------------------------------------------------------
      Use dbsr_hf
      Use df_orbitals

      Implicit none
      Integer :: i,j,start,i1,i2
      Character(5), external :: Eli

      Call Read_confs;  if(nconf.gt.0) Return

      Call Clean_a(configuration)
      nwf = ncore
      start = 1
      Do
       i = index(configuration(start:),')')
       if(i.eq.0) Exit
       nwf = nwf +1
       start = start+i
      End do

      if(nwf.eq.0) Stop 'Stop dbsr_hf: input number of orbitals, nwf=0'

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

      start = 1
      Do i=ncore+1,nwf
       i1 = index(configuration(start:),'(')+start-1
       i2 = index(configuration(start:),')')+start-1
       ebs(i) = configuration(start:i1-1)
       read(configuration(i1+1:i2-1),*) qsum(i)
       Call EL_NLJK(ebs(i),nbs(i),kbs(i),lbs(i),jbs(i),ibs(i))
       start = i2+1
       j = INT(qsum(i)); if(j.eq.jbs(i)+1) clsd(i)=.TRUE.
      End do

      Do i=1,nwf; ebs(i)=ELi(nbs(i),kbs(i),0); End do

      nconf=1
      Allocate(weight(nconf)); weight=1.d0
      Allocate(iqconf(nconf,nwf)); iqconf(nconf,1:nwf)=qsum(1:nwf)

      End Subroutine Def_conf


!======================================================================
      Subroutine Read_confs
!======================================================================
!     define configurations and their weights from conf-file;
!     if something wrong with conf-file - return to 1 configuration
!     from the input file
!----------------------------------------------------------------------
      Use dbsr_hf, iq_xx => iq
      Use df_orbitals

      Implicit none
      Character :: line*200, EL*5
      Integer :: i,j,m, n,l,k,iset,iq,io, start, i1,i2, ii,jj,iw
      Integer, external :: Ifind_orb, ndets_jq, Ifind_position
      Character(5), external :: Eli

! ... check conf-file:

      nconf = 0
      if(len_trim(name).ne.0) AF_conf = trim(name)//BF_conf
      Call Read_aarg('confs',AF_conf)


      if(Icheck_file(AF_conf).eq.0) Return

! ... define number of configurations:

      Open(nuc,file=AF_conf)
      rewind(nuc)
    1 read(nuc,'(a)',end=2) line
      if(INDEX(line,'(').eq.0) go to 1
      if(line(1:3).eq.'CSF') go to 1
      if(line(1:1).eq.'*') goto 2
      nconf = nconf + 1
      go to 1
    2 Continue
      if(nconf.eq.0) Return

! ... define core:

      i=Ifind_position(nuc,'Core subshells:')
      if(i.eq.0) Stop &
       'Stop in Reads_conf: no "Core subshells:" flag in the conf-file'
      read(nuc,'(/a)') core
      Call Def_core

! ... define one-electron orbitals:

      i=Ifind_position(nuc,'Peel subshells:')
      if(i.eq.0) Stop &
       'Stop in Reads_conf: no "Peel subshells:" flag in the conf-file'
      read(nuc,'(/a)') config

      m=1000; Allocate(ebs(m))
      Do i=1,m
       read(config,*,IOSTAT=ii) ebs(1:i)
       if (ii /= 0) Exit
      End do
      if(i.le.1) Stop 'Stop in Reads_conf: number of peel orbitals = 0'
      write(config,'(50a5)') (ebs(j),j=1,i-1)
      Deallocate(ebs)

      nwf = i -1 + ncore
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
       qsum(i) = 1
       k = k + 5
      End do

      Do i=1,nwf; ebs(i)=ELi(nbs(i),kbs(i),0); End do
      write(config,'(50a5)') (ebs(j),j=ncore+1,nwf)

! ... read configurations and define their weights:

      Call Read_ipar(inp,'eal',eal)
      Call Read_iarg('eal',eal)

      Allocate(weight(nconf));     weight = 1.d0
      Allocate(iqconf(nconf,nwf)); iqconf = 0
      Do i=1,ncore;  iqconf(:,i) = qsum(i); End do

      i = 0
      rewind(nuc)
   10 read(nuc,'(a)',end=20) line
      if(INDEX(line,'(').eq.0) go to 10
      if(line(1:3).eq.'CSF') go to 10
      if(line(1:1).eq.'*') goto 20

      i = i + 1;  iw = 1
      Call Clean_a(line)
      start = 1
      Do
       i1 = index(line(start:),'(')+start-1
       i2 = index(line(start:),')')+start-1
       if(i1.lt.start.or.i2.lt.start) Exit
       EL = line(start:i1-1)
       read(line(i1+1:i2-1),*) iq
       Call EL_NLJK(EL,n,k,l,jj,iset)
       io=Ifind_orb(n,k,iset)
       if(io.eq.0) Stop 'unknown orbital in conf-file'
       iqconf(i,io) = iq
       iw = iw * ndets_jq(jj,iq)
       start = i2+1
      End do
      if(eal.eq.5) weight(i)=iw
      if(eal.eq.9) read(line(start:),*) weight(i)

      go to 10
   20 Continue

      if(SUM(weight).eq.0.d0) weight = 1.d0
      weight = weight / SUM(weight)

! ... redefine qsum:

      qsum = 0.d0
      Do io=1,nwf
       Do i=1,nconf
        qsum(io)=qsum(io) + iqconf(i,io)*weight(i)
       End do
       i = INT(qsum(io)); if(i.eq.jbs(i)+1) clsd(io)=.TRUE.
      End do

      End Subroutine Read_confs


!======================================================================
      Subroutine Record_confs
!======================================================================
! ... record the existing configurations in conf-file
!----------------------------------------------------------------------
      Use dbsr_hf
      Use df_orbitals

      Implicit none
      Integer :: i,j,start

      AF_conf = trim(name)//trim(BF_conf)
      Call Read_aarg('confs',AF_conf)
      Open(nuc,file=AF_conf); rewind(nuc)
      write(nuc,'(a)') 'Core subshells:'
      write(nuc,'(20a5)') (ebs(i),i=1,ncore)
      write(nuc,'(a)') 'Peel subshells:'
      write(nuc,'(20a5)') (ebs(i),i=ncore+1,nwf)
      write(nuc,'(a)') 'CSF(s):'

      Do i = 1,nconf

      start = 1
      Do j=ncore+1,nwf;  if(iqconf(i,j).eq.0) Cycle
       write(configuration(start:),'(a5,a1,i2,a1)') &
        ebs(j),'(',iqconf(i,j),')'
       start = start + 9
      End do

      j = len_trim(configuration)
      if(j.le.72) then
       write(nuc,'(a,T73,F12.4)') configuration(1:j),weight(i)
      else
       write(nuc,'(a,F12.4)') configuration(1:j),weight(i)
      end if

      End do ! over configurations

      write(nuc,'(a)') '*'

      End Subroutine Record_confs



!======================================================================
      Subroutine Def_nit(string)
!======================================================================
!     This routine obtains information about the problem to be solved
!     and how the spline methods are to be applied
!----------------------------------------------------------------------
      Use dbsr_hf
      Use df_orbitals

      Implicit none
      Character(*) :: string
      Integer :: i,j,n,l,k,ip,start,iset
      Integer, external :: Ifind_orb
      Character(5) :: EL

      Call Clean_a(string)
      iord = 0

      if(string(1:3)=='ALL' .or. string(1:3)=='all' .or. &
         LEN_TRIM(string) == 0 ) then
       nit = nwf
       Do i=1,nwf; iord(i)=i; End do
      elseif(string(1:4)=='NONE' .or. string(1:4)=='none') then
       nit = 0
      elseif (INDEX(string,'n=') /= 0) then
       read(string(i+2:),*) n
       k=0; Do i=1,nwf; if(nbs(i).ne.n) Cycle; k=k+1; iord(k)=i; End do
       nit = k
      elseif (INDEX(string,'n>') /= 0) then
       read(string(i+2:),*) n
       k=0; Do i=1,nwf; if(nbs(i).le.n) Cycle; k=k+1; iord(k)=i; End do
       nit = k
      elseif (INDEX(string,'=') /= 0) then
       i = INDEX(string,'=')
       read(string(i+1:),*) nit
       k=0; Do i=nwf-nit+1,nwf; k=k+1; iord(k)=i; End do
      else
       start = 1; ip = 0
       Do
        i = index(string(start:),',')
        if (i /= 0 .or. LEN_TRIM(string(start:)) /= 0) then
         read(string(start:),*) EL
         Call EL_NLJK(EL,n,k,l,j,iset)
         j = Ifind_orb(n,k,iset)
         if(j.gt.0) then; ip=ip+1; iord(ip)=j; end if
        end if
        start = start + i
        if(i == 0 .or.  LEN_TRIM(string(start:)) == 0) Exit
       End do
       nit = ip
      end if

      Do j=1,nwf; i=iord(j)
       if(i.eq.0) Cycle
       if(qsum(i).eq.0.d0) iord(j)=0
      End do

      nit = 0;  Do j=1,nwf; if(iord(j).gt.0) nit=nit+1; End do

      End Subroutine Def_nit



!======================================================================
      Subroutine write_inp
!======================================================================
!     This routine prepare default input file
!----------------------------------------------------------------------
      Use dbsr_hf

      Implicit none

      AF_dat = trim(name)//BF_dat
      Call Read_aarg('dat',AF_dat)
      Open(inp,file=AF_dat)
      rewind(inp)

      write(inp,'( a/)')    'Main atomic parameters:'
      write(inp,'(a,a)')    'atom    =  ',trim(atom)
      if(z.eq.0.d0) z=an
      write(inp,'(a,f8.4)') 'z       = ' ,z
      write(inp,'(a,f8.4)') 'atw     = ' ,atw
      write(inp,'(a,f8.4)') 'rms     = ' ,rms
      if(ion.ne.atom) &
      write(inp,'(a,a)')    'ion     =  ',trim(ion)
      write(inp,'(a,a)')    'core    =  ',trim(adjustl(core))
      write(inp,'(a,a)')    'conf    =  ',trim(adjustl(configuration))
      write(inp,'(a,a)')    'varied  =  ',trim(anit)
      write(inp,'(a,a)')    'term    =  ',trim(term)
      write(inp,'(a,i1)')   'eal     =  ',eal

      write(inp,'(/a/)') 'Running parameters:'
      write(inp,'(a,i2,T40,a)') 'mbreit  = ',mbreit, &
                '- Breit interaction mode'
      write(inp,'(a,i2,T40,a)') 'mode_SE = ',mode_SE, &
                '- self-energy correction mode'
      write(inp,'(a,i2,T40,a)') 'mode_VP = ',mode_VP, &
                '- vacuum-polarization correction mode'
      write(inp,'(a,1PE9.2,T40,a)') 'scf_tol = ',scf_tol, &
                '- Energy convergence tolerance'
      write(inp,'(a,1PE9.2,T40,a)') 'orb_tol = ',orb_tol, &
                '- Orbital convergence tolerance'
      write(inp,'(a,1PE9.2,T40,a)') 'end_tol = ',end_tol, &
                '- Orbital tail cut-off'
      write(inp,'(a,i3,T40,a)') 'max_it  = ',max_it, &
                '- max. number of iterations'
      write(inp,'(a,i2,T40,a)') 'ilzero  = ',ilzero, &
                '- initial zero B-splines'
      write(inp,'(a,i2,T40,a)') 'ibzero  = ',ibzero, &
                '- minimum zero B-splines in the end'
      write(inp,'(a,i2,T40,a)') 'rotate  = ',rotate, &
                '- use orbital rotation'
      write(inp,'(a,i2,T40,a)') 'debug   = ',debug,  &
                '- additional debug output'

      write(inp,'(/a,a)')  'All parameters from input files ', &
                            'can be replaced from command line as:'
      write(inp,'(/a)') 'dbsr_hf name par1=value par2=value par3=value ... '

      write(inp,'(/a/)') 'Name-driven file-names and key-words for their re-definition:'

      write(inp,'(a,T15,a,T40,a)') 'name'//BF_dat, 'dat=...', '- input parameters'
      write(inp,'(a,T15,a,T40,a)') 'name'//BF_log, 'log=...', '- run output and summry'
      write(inp,'(a,T15,a,T40,a)') 'name'//BF_grid,'knot=...','- B-spline parameters'
      write(inp,'(a,T15,a,T40,a)') 'name'//BF_inp, 'inp=...', '- input wavefunctions if any'
      write(inp,'(a,T15,a,T40,a)') 'name'//BF_out, 'out=...', '- output wavefunctions'
      write(inp,'(a,T15,a,T40,a)') 'name'//BF_plt, 'plot=...','- plot wavefunctions'
      write(inp,'(a,T15,a,T40,a)') 'name'//BF_nl,  'nl=...',  '- output w.f. for whole spectrum of outer electron'
      write(inp,'(a,T15,a,T40,a)') 'name'//BF_rwf, 'w=...',   '- output w.f. in GRASP format'
      write(inp,'(a,T15,a,T40,a)') 'name'//BF_conf,'confs=..','- input relativistic configurations if any'
      write(inp,'(a,T15,a,T40,a)') 'name'//BF_cfg, 'c=...',   '- input jj atomic states if any'
      write(inp,'(a,T15,a,T40,a)') 'name'//BF_LS,  'LS=...',  '- input LS configurations if any'

      write(inp,'(/a/)') 'Call dbsr_hf as:'

      write(inp,'(a)') 'dbsr_hf name  -  name.inp is supposed to be created by user '

      write(inp,'(a)') 'dbsr_hf atom  -  choose the atom with given atomic symbol'

      write(inp,'(a)') 'dbsr_hf an=.. -  choose the atom with atomic number an'

      write(inp,'(a)') 'dbsr_hf name  an=...  ion=... -  choose needed ion for given an'


      write(inp,'(80("_"))')
      write(inp,'(/a/)') '! Additional information for input parameters:'

      write(inp,'(a)') '! term=AV  - optimize the configuration given as conf=...'
      write(inp,'(a)') '!            or set of configurations given in the name.conf file'
      write(inp,'(a)') '! term=LS  - optimize all jj-configurations consistent with LS '
      write(inp,'(a)') '!            configuration given as conf=... or set of configurations'
      write(inp,'(a)') '!            given in the name.LS file'
      write(inp,'(a)') '! term=jj  - optimize specific jj atomic states given in name.c file'
      write(inp,'(a)') '!            which is prepared by user or created from previous DBSR_HF run'
      write(inp,'(a)') '! eal      - indicates the mode for the state weights:'
      write(inp,'(a)') '!            =1 - equally weighted'
      write(inp,'(a)') '!            =5 - statistically weighed, default'
      write(inp,'(a)') '!            =9 - defined by user in .conf or .c files'
      write(inp,'(a)') '!            '
      write(inp,'(a)') '! varied   - possible options -> all, none, list of nl, =last, n=..., n>...'
      write(inp,'(a)') '!            '
      write(inp,'(a)') '! mbreit   - Breit interaction mode:'
      write(inp,'(a)') '!            =0 - no Breit interaction is included'
      write(inp,'(a)') '!            =1 - Breit interaction and QED correction is included '
      write(inp,'(a)') '!                 as first-order correction '
      write(inp,'(a)') '!            =2 - Breit interaction is included into orbital optimization'
      write(inp,'(a)') '!                 and QED as first order corrections'
      write(inp,'(a)') '!            '
      write(inp,'(a)') '! mode_SE  - self-energy correction mode:'
      write(inp,'(a)') '!            = 1 - hydrogen-based results'
      write(inp,'(a)') '!            = 2 - Welton concept'
      write(inp,'(a)') '!            = 3 - GRASP mode'
      write(inp,'(a)') '!            '
      write(inp,'(a)') '! out_w=1    - additional output of w.f. in GRASP format, name.w '
      write(inp,'(a)') '! out_nl=1   - additional output of name.nl with w.f. for outer electron'
      write(inp,'(a)') '! out_plot=1 - additional output in name.plot'
      write(inp,'(80("_"))')

      rewind(inp)

      End Subroutine write_inp

