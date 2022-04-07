!======================================================================
      Subroutine Read_Grasp(nu)
!======================================================================
!     Generation double (p,q) B-spline representation for orbitals
!     given in GRASP package format, w-files
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use DF_orbitals
      Use dbsr_hf

      Implicit none

      Real(8), allocatable :: R(:),PG(:),QG(:)
      Real(8), allocatable :: b3(:),c3(:),d3(:),b4(:),c4(:),d4(:)
      Real(8) :: pr(nv,ks),qr(nv,ks),a(ns,ns),Pcoef(ns),Qcoef(ns)
      Character(6) string
      Integer :: nu,n,kappa,npts,is,i,j,ii,np,nq,ip,iv
      Real(8) :: energy,p0,pn,PP,QQ,d,pm,qm
      Integer, external :: Ifind_orb
      Real(8), external :: SEVAL, bvalu2, QUADR

! ... check the file:

      read(nu) string
      if(string.ne.'G92RWF') Stop 'This is not a GRASP Radial File'

      write(log,'(/a/)') 'GRASP orbitals:'

! ... process the radial functions:

    1 read(nu,end=2)  n,kappa,energy,npts
      if(allocated(R)) Deallocate(R,PG,QG)
      Allocate(R(npts),PG(npts),QG(npts))
      READ(nu) p0,PG,QG
      READ(nu) R

! ... check if we need this orbital:

      is = n/1000; n = n - is*1000
      ii = Ifind_orb(n,kappa,is)
      if(ii.eq.0) go to 1
      e(ii,ii) = energy

! ... remove unphysical oscilations in the end:

      PM = maxval(abs(PG))
      Do i=npts,1,-1
       if(abs(PG(i))/PM.gt.end_tol) Exit; PG(i)=0.d0
      End do
      np = i

      QM = maxval(abs(QG))
      Do i=npts,1,-1
       if(abs(QG(i))/QM.gt.end_tol) Exit; QG(i)=0.d0
      End do
      nq = i

! ... interpolate function p(r) to cubic splines:

      if(allocated(b3)) Deallocate(b3,c3,d3,b4,c4,d4)
      Allocate(b3(np),c3(np),d3(np),b4(nq),c4(nq),d4(nq))
      Call Splin3(np,R,PG,b3,c3,d3)
      Call Splin3(nq,R,QG,b4,c4,d4)

! ... evaluate the function in gaussian points:

      pr = 0.d0; qr = 0.d0
      Do i=1,nv; Do j=1,ks
       if(gr(i,j).lt.R(np)) pr(i,j)=SEVAL(np,gr(i,j),R,PG,b3,c3,d3)*grw(i,j)
       if(gr(i,j).lt.R(nq)) qr(i,j)=SEVAL(nq,gr(i,j),R,QG,b4,c4,d4)*grw(i,j)
      End do; End do

! ... form the vector of inner products of the radial function and the
! ... spline basis functions

      Pcoef = 0.d0
      Do iv = 1,nv; Do ip = 1,ksp; i = iv+ip-1
       Pcoef(i) = Pcoef(i) + SUM(pr(iv,:)*pbsp(iv,:,ip))
      End do; End do

      Qcoef = 0.d0
      Do iv = 1,nv; Do ip = 1,ksq; i = iv+ip-1
       Qcoef(i) = Qcoef(i) + SUM(qr(iv,:)*qbsp(iv,:,ip))
      End do; End do

! ... solve the system of equations for coef

      a(1:nsp-1,1:nsp-1)=fpbs(2:nsp,2:nsp)
      Call gaussj (a,nsp-1,ns,Pcoef(2),1,ns)
      Pcoef(1)=0.d0

      a(1:nsq-1,1:nsq-1)=fqbs(2:nsq,2:nsq)
      Call gaussj (a,nsq-1,ns,Qcoef(2),1,ns)
      Qcoef(1)=0.d0

      p(:,1,ii) = Pcoef
      p(:,2,ii) = Qcoef
      Call Check_tails(ii)

! ... check the normalization

      PN = sqrt(QUADR(p(1,1,ii),p(1,1,ii),0))
      p(:,:,ii)=p(:,:,ii)/PN

! ... max.deviation from original orbital:

      pm = 0.d0
      Do j=2,np
       PP = bvalu2(tp,p(1,1,ii),nsp,ksp,R(j),0)
       d = abs(PP-PG(j));  if(d.gt.pm) pm = d
      End do

      qm = 0.d0
      Do j=2,nq
       QQ = bvalu2(tq,p(1,2,ii),nsq,ksq,R(j),0)
       d = abs(QQ-QG(j)); if(d.gt.qm) qm = d
      End do

      write(log,'(a5,4(a,E12.3))')  &
           ebs(ii), '  diff_p =',pm,'   diff_q =',qm,'   norm =',PN

      go to 1
    2 Close(nu)

      End Subroutine Read_GRASP


