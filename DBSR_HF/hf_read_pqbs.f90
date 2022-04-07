!======================================================================
      Subroutine Read_pqbs(nu)
!======================================================================
!     read B-spline w.f. from bsw-file (unit nu)
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DF_orbitals
      Use dbsr_hf, iq_xx => iq

      Implicit none
      Integer, intent(in) :: nu
      Integer :: i,j,k,ip,iq,ka,l,n,m,itype,nsw,ksw,mw,kp,kq
      Character(5) :: elw
      Integer, external :: Ifind_orb
      Real(8), allocatable :: tw(:),pw(:),qw(:)
      Real(8) :: S

! ... read the written B-spline grid and check if it matches

      rewind(nu)
      read(nu) itype,nsw,ksw
      rewind(nu)
      allocate(tw(nsw+ksw),pw(nsw),qw(nsw))
      read(nu) itype,nsw,ksw,tw,kp,kq
      ip = ksw - kp
      iq = ksw - kq
      k=1
      if(ksw.ne.ks) k=0
      if(nsw.ne.ns) k=0
      if(ksp.ne.kp) k=0
      if(ksq.ne.kq) k=0
      Do i=1,ns+ks; if(k.eq.0) Exit
       if(abs(t(i)-tw(i)).lt.1.d-12) Cycle; k=0; Exit
      End do

! ... read radial functions and converte them if necessary

    1 read(nu,end=2) elw,mw, S
      pw=0.d0; read(nu) pw(1:mw)
      qw=0.d0; read(nu) qw(1:mw)
      Call EL_NLJK(elw,n,ka,l,j,i)
      m = Ifind_orb(n,ka,i)
      if(m.eq.0) go to 1    ! skip that orbital
      e(m,m) = S
      if(k.eq.1) then
       mbs(m)=mw
       p(:,1,m) = pw
       p(:,2,m) = qw
       write(log,'(a,a,a)') ebs(m),' - read from file ', trim(AF_inp)
      else
       Call Convert_pq(nsw-ip,kp,tw(1+ip),pw,nsp,ksp,p(1,1,m),pbsp,fpbs,1,0)
       Call Convert_pq(nsw-iq,kq,tw(1+iq),qw,nsq,ksq,p(1,2,m),qbsp,fqbs,1,0)
       mbs(m) = nsp-1
       write(log,'(a,a,a,a)') ebs(m),' - read from file ', trim(AF_inp),' and converted'
      end if

      go to 1
    2 Close(nu)

      if(allocated(tw)) deallocate(tw,pw,qw)

      End Subroutine Read_pqbs


!======================================================================
      Subroutine write_pqbs
!======================================================================
!     record B-spline w.f. to file (unit nuw)
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DF_orbitals
      Use dbsr_hf

      Implicit none
      Integer :: i

      AF_out = trim(name)//BF_out
      Call Read_ipar(inp,'out',AF_out)
      Call Read_iarg('out',AF_out)
      open(nuw,file=AF_out,form='UNFORMATTED')

      rewind(nuw)
      write(nuw) grid_type,ns,ks,t(1:ns+ks),ksp,ksq
      Do i=1,nbf
       write(nuw) ebs(i),mbs(i),e(i,i)
       write(nuw) p(1:mbs(i),1,i)
       write(nuw) p(1:mbs(i),2,i)
      End do
      Close(nuw)

      End Subroutine write_pqbs



!======================================================================
      Subroutine write_nl
!======================================================================
!     record B-spline w.f. to file (unit nuw)
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBSR_HF

      Implicit none
      Integer :: i

      if(out_nl.eq.0) Return
      AF_nl = trim(name)//BF_nl
      Call Read_ipar(inp,'nl',AF_nl)
      Call Read_iarg('nl',AF_nl)
      open(nuw,file=AF_nl,form='UNFORMATTED')

      rewind(nuw)
      write(nuw) grid_type,ns,ks,t(1:ns+ks),ksp,ksq
      Do i=1,nsol_nl
       write(nuw) e_nl(i)
       write(nuw) p_nl(1:ns,i)
       write(nuw) p_nl(ns+1:ns+ns,i)
      End do
      Close(nuw)

      End Subroutine write_nl

