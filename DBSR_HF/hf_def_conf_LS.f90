!======================================================================
      Subroutine Def_conf_LS
!======================================================================
!     Defines rel. configurations to be included into HF optimization
!     based on the LS input configurations
!----------------------------------------------------------------------
      Use dbsr_hf

      Implicit none
      Character(160) :: conf
      Integer :: iqsum(msh),iq_min(msh),iq_max(msh)
      Integer :: i,j,ii,m

      Integer, parameter :: morb=100
      Character(5) :: EL(morb), EL1, ELi
      Integer :: norb=0

! ... create file name.LS if absent:

      AF_LS = trim(name)//BF_LS
      Call Read_aarg('LS',AF_LS)
      if(Icheck_file(AF_LS).eq.0) then
       Call Decode_conf_jj(configuration,no,nn,kn,ln,jn,iq,in)
       if(no.eq.0) Stop 'Stop in Def_conf_LS: no = 0'
       Call Reduce_jj_LS(no,nn,kn,ln,jn,iq,in)
       Call Incode_conf_LS(configuration,no,nn,ln,iq,in)
       open(nus,file=AF_LS)
       write(nus,'(a)') atom
       write(nus,'(a)') trim(core)
       write(nus,'(a)') trim(configuration)
       write(nus,'(a)') '*'
       close(nus)
      end if

! ... generate all possible rel. configurations:

      open(nus,file=AF_LS)
      open(nua,file='LS')
      rewind(nus)
      read(nus,'(/a)') core
      nconf=0

      Do
       read(nus,'(a)',end=10,err=10) conf
       if(index(conf,'(').eq.0) Cycle
       Call Decode_conf_LS(conf,no,nn,ln,iq,in)
       nelc = sum(iq(1:no))

       Call Reduce_LS_jj (no,nn,kn,ln,jn,iq,in,iq_min,iq_max)
       iqsum(1:no) = iq(1:no)
       if(len_trim(conf_LS).eq.0) conf_LS=conf


! ... check new orbitals:

        Do i=1,no; if(iq_max(i).le.0) Cycle;  EL1 = ELi(nn(i),kn(i),in(i))
         m=0
         Do j=1,norb
          if(EL(j).ne.EL1) Cycle; m = 1; Exit
         End do
         if(m.eq.1) Cycle
         norb=norb+1; EL(norb)=EL1
        End do

        i=no; iq=iq_max
      1 ii=SUM(iq(1:no))

        if(ii.eq.nelc) then

         m = 0
         Do j=1,no; if(ln(j).eq.0) Cycle
          if(j.gt.1) then; if(ln(j).eq.ln(j-1)) Cycle; end if
          if(iq(j)+iq(j+1).eq.iqsum(j+1)) Cycle
          m = 1; Exit
         End do
         if(m.eq.0) Call Record_conf

        elseif(ii.lt.nelc.and.i.lt.no) then
         i=i+1;  iq(i)=iq_max(i); go to 1
        end if

      2 iq(i)=iq(i)-1

        if(iq(i).lt.iq_min(i)) then
         iq(i) = iq_min(i)
         if(i.eq.1) go to 3
         i=i-1; go to 2
        end if

        go to 1
      3 Continue

      End do
   10 if(nconf.eq.0) Stop 'Stop in Def_conf_LS: nconf=0'

      AF_conf = trim(name)//BF_conf
      Call Read_aarg('confs',AF_conf)
      open(nuc,file=AF_conf)
      write(nuc,'(a)') 'Core subshells:'
      write(nuc,'(a)') trim(core)
      write(nuc,'(a)') 'Peel subshells:'
      write(nuc,'(20a5)') (EL(i),i=1,norb)
      write(nuc,'(a)') 'CSF(s):'

      rewind(nua)
      Do i=1,nconf
       read(nua,'(a)') conf
       write(nuc,'(a)') trim(conf)
      End do
      write(nuc,'(a)') '*'
      Close(nuc)

      Close(nua,status='DELETE')
      Close(nus)

CONTAINS

!======================================================================
      Subroutine Record_conf
!======================================================================
! ... record the configuration in conf-file
!----------------------------------------------------------------------
      Implicit none
      Integer :: i
      Real(8) :: W
      Integer, external :: Ndets_jq

      Call Incode_conf_jj(conf,no,nn,kn,ln,jn,iq,in)

! ... statistical weight:

      W = 1.d0
      Do i=1,no; if(iq(i).le.0) Cycle
       W = W * Ndets_jq(jn(i),iq(i))
      End do

! ... record configuration

      i = len_trim(conf)
      if(i.le.72) then
       write(nua,'(a,T73,F12.4)') conf(1:i),W
      else
       write(nua,'(a,F12.4)') conf(1:i),W
      end if
      nconf = nconf + 1


      End Subroutine Record_conf

      End Subroutine Def_conf_LS


!=========================================================================
      Subroutine Decode_conf_jj(configuration,no,nn,kn,ln,jn,iq,in)
!=========================================================================
! ... decode the spectroscopic configuration into "integer" representation
!-------------------------------------------------------------------------
      Implicit none
      Character(*) :: configuration
      Integer :: nn(*),kn(*),ln(*),jn(*),iq(*),in(*)
      Integer :: no,start,i1,i2
      Character(5) :: EL

      Call Clean_a(configuration)
      no=0; start=1
      Do
       i1 = index(configuration(start:),'(')
       if(i1.eq.0) Exit
       i1=i1+start-1
       i2 = index(configuration(start:),')')+start-1
       no = no + 1
       read(configuration(i1+1:i2-1),*) iq(no)
       EL = configuration(start:i1-1)
       Call EL_NLJK(EL,nn(no),kn(no),ln(no),jn(no),in(no))
       start = i2+1
      End do

      End Subroutine Decode_conf_jj

!======================================================================
      Subroutine Incode_conf_jj(configuration,no,nn,kn,ln,jn,iq,in)
!======================================================================
!     incodes the configuration from INTEGER format to c-file format
!----------------------------------------------------------------------
      Implicit none
      Character(*) :: configuration
      Integer :: nn(*),kn(*),ln(*),jn(*),iq(*),in(*)
      Integer :: no,i,m
      Character(5), external :: ELi

      m=0; configuration=' '
      Do i=1,no
       if(iq(i).le.0) Cycle
       configuration(m+1:m+5) = ELi(nn(i),kn(i),in(i))
       write(configuration(m+6:m+9),'(a1,i2,a1)') '(',iq(i),')'
       m = m + 9
      End do

      End Subroutine Incode_conf_jj

!======================================================================
      Subroutine Reduce_jj_LS(no,nn,kn,ln,jn,iq,in)
!======================================================================
! ... convert jj- to LS-configuration
!----------------------------------------------------------------------
      Implicit none
      Integer :: no,i,n
      Integer :: nn(*),kn(*),ln(*),jn(*),iq(*),in(*)

      n=1
      Do i=2,no
       if(nn(n).eq.nn(i).and.ln(n).eq.ln(i)) then
        iq(n)=iq(n)+iq(i)
       else
        n=n+1
        nn(n)=nn(i); ln(n)=ln(i); iq(n)=iq(i); in(n)=in(i)
       end if
      End do
      no = n

      End Subroutine Reduce_jj_LS


!=========================================================================
      Subroutine Incode_conf_LS(configuration,no,nn,ln,iq,in)
!=========================================================================
! ... incode the spectroscopic configuration into "integer" representation
!-------------------------------------------------------------------------
      Implicit none
      Character(*) :: configuration
      Integer :: nn(*),ln(*),iq(*),in(*)
      Integer :: no,i,m
      Character(1), external :: AL

      m=1; configuration = ' '
      Do i=1,no
       if(in(i).eq.0) then
        write(configuration(m:),'(i3,a1)') nn(i),AL(ln(i),1)
       else
        write(configuration(m:),'(i2,a1,i1)') nn(i),AL(ln(i),1),in(i)
       end if
       write(configuration(m+4:),'(a1,i2,a1)') '(',iq(i),')'
       m=m+8
      End do

      End Subroutine Incode_conf_LS


!=========================================================================
      Subroutine Decode_conf_LS(configuration,no,nn,ln,iq,in)
!=========================================================================
! ... decode the spectroscopic configuration into "integer" representation
!-------------------------------------------------------------------------
      Implicit none
      Character(*) :: configuration
      Integer :: nn(*),ln(*),iq(*),in(*)
      Integer :: no,start,i1,i2,k,j
      Character(5) :: EL

      Call Clean_a(configuration)
      no=0; start=1
      Do
       i1 = index(configuration(start:),'(')
       if(i1.eq.0) Exit
       i1=i1+start-1
       i2 = index(configuration(start:),')')+start-1
       no = no + 1
       read(configuration(i1+1:i2-1),*) iq(no)
       EL = configuration(start:i1-1)
       Call EL_NLJK(EL,nn(no),k,ln(no),j,in(no))
       start = i2+1
      End do

      End Subroutine Decode_conf_LS


!======================================================================
      Subroutine Reduce_LS_jj(no,nn,kn,ln,jn,iq,in,iq_min,iq_max)
!======================================================================
! ... convert LS- to jj-configuration
!----------------------------------------------------------------------
      Implicit none
      Integer :: no,i,n,l,k,j,jj
      Integer :: nn(*),kn(*),ln(*),jn(*),iq(*),in(*),iq_min(*),iq_max(*)
      Integer :: n1(no),l1(no),q1(no),i1(no)

      n1(1:no) = nn(1:no)
      l1(1:no) = ln(1:no)
      q1(1:no) = iq(1:no)
      i1(1:no) = in(1:no)

      n=0
      Do i=1,no; l=l1(i)

       if(l.eq.0) then
        j=l+l+1; k=(l+l-j)*(j+1)/2
        n = n + 1
        nn(n) = n1(i); kn(n)=k; ln(n)=l; jn(n)=j; iq(n)=0; in(n)=i1(i)
        iq_min(n)=q1(i); iq_max(n)=q1(i)
        Cycle
       end if

       j=l+l-1; k=(l+l-j)*(j+1)/2; jj=l+l+1
       n = n + 1
       nn(n) = n1(i); kn(n)=k; ln(n)=l; jn(n)=j; iq(n)=0; in(n)=i1(i)
       iq_min(n)=max(0,q1(i)-jj-1); iq_max(n)=min(q1(i),j+1)

       j=l+l+1; k=(l+l-j)*(j+1)/2; jj=l+l-1
       n = n + 1
       nn(n) = n1(i); kn(n)=k; ln(n)=l; jn(n)=j; iq(n)=q1(i); in(n)=i1(i)
       iq_min(n)=max(0,q1(i)-jj-1); iq_max(n)=min(q1(i),j+1)

      End do
      no = n

      End Subroutine Reduce_LS_jj
