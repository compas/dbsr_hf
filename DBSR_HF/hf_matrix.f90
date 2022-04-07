!=====================================================================
      Subroutine hf_matrix(i,hfm)
!=====================================================================
!     Set up the hf_matrix "hfm" for orbital i
!     Call(s):  a,b - angular coeficient routines
!-----------------------------------------------------------------------
      Use dbsr_hf                              !,     only: kmax,log
      Use DBS_grid,    only: ns,ks,ms
      Use DBS_dhl_pq,  only: dhl
      Use df_orbitals, only: kbs,nbf,qsum,p

      Implicit none
      Integer, intent(in) :: i
      Real(8), intent(out) :: hfm(ms,ms)
      Real(8) :: d(ns,ks),dd(ns,ks),x(ns,ns),xx(ns,ns)
      Real(8) :: c
      Real(8), external :: a,b
      Integer :: j,ii,jj,k,ip,met

      Real(8) :: t1,t2
      Real(8), external :: RRTC

      t1 = RRTC()

! ... one-electron integral contribution

      Call MAT_dhl_pq(kbs(i))
      hfm = qsum(i)*dhl

! ... two-electron contribution for direct and exchage integrals
! ... this contribution is devided on 4 parts: pppp, pqpq, qpqp, qqqq
! ... and added separately

      Do ii=1,2; Do jj=1,2   ! pppp, pqpq, qpqp, qqqq

      Do k=0,kmax

      Call mrk_aaaa(k,ii,jj)

! ... direct contribution

       dd = 0.d0; met=0
       Do j = 1,nbf
        c = a(i,j,k); if(c.eq.0.d0) cycle; if(i.eq.j) c=c+c
        Call density(ns,ks,d,p(1,jj,j),p(1,jj,j),'s')
        dd = dd + c*d; met=met+1
       End do
       if(met.gt.0) Call Convol(ns,ks,d,dd,1,'s','s')
       if(met.gt.0) Call Update_hs(ms,hfm,ii,ii,ns,ks,d,'s')

! ... exchange contribution

       xx = 0.d0; met=0
       Do j = 1,nbf
        c = b(i,j,k); if (c.eq.0.d0) cycle
        Call density(ns,ks,x,p(1,ii,j),p(1,jj,j),'x')
        xx = xx + x*c; met=met+1
       End do
       if(met.gt.0) Call Convol(ns,ks,x,xx,4,'s','s')
       if(met.gt.0) Call Update_hs(ms,hfm,ii,jj,ns,ks,x,'x')

      End do ! over k
      End do; End do ! over itype

      if(mbreit.eq.2) Call hf_matrix_breit(i,hfm)

      hfm = hfm / qsum(i)

      t2 = RRTC()
      time_hf_matrix = time_hf_matrix + t2-t1

      End Subroutine hf_matrix


!======================================================================
      Subroutine mrk_aaaa(k,ii,jj)
!======================================================================
      Integer, intent(in) :: k,ii,jj

      Select Case (10*ii+jj)
       Case(11);  Call mrk_pppp(k)   !  (. p p .)          1 1
       Case(12);  Call mrk_pqpq(k)   !  (. q p .)  jj ii   1 2
       Case(21);  Call mrk_qpqp(k)   !  (. p q .)  jj ii   2 1
       Case(22);  Call mrk_qqqq(k)   !  (. q q .)          2 2
       Case default; Stop 'unknown case in mrk_aaaa'
      End Select

      End Subroutine mrk_aaaa


!=====================================================================
      Subroutine hf_matrix_breit(i,hfm)
!=====================================================================
!     Set up the hf_matrix "hfm" for orbital i with Breit interaction
!-----------------------------------------------------------------------
      Use dbsr_hf
      Use rk4_data
      Use DBS_grid,    only: ns,ks,ms
      Use DBS_dhl_pq,  only: dhl
      Use df_orbitals, only: kbs,nbf,qsum,p

      Implicit none
      Integer, intent(in) :: i
      Real(8), intent(out) :: hfm(ms,ms)

      Character :: sym, sym_int(4)
      Integer :: itype,k,kk,met, int,jnt, ik, j, kint(4), ip1,ip2, jp1,jp2
      Real(8) :: x(ns,ns),xx(ns,ns)
      Real(8) :: c

      Real(8) :: t1,t2
      Real(8), external :: RRTC

      Data sym_int /'n','n','x','x'/
      Data kint /2,1,4,3/

      t1 = RRTC()

      Do itype=0,1
      Do k=0,kmax; kk = k + 1000*itype
       if(itype.eq.0)  Call msk_ppqq(k)
       if(itype.eq.1)  Call msk_pqqp(k)
      Do int = 1,4; sym = sym_int(int)

! ... this can be replaced by arrays

!---------------------------------------------------------------------------------------------
      Select case(itype)
       Case(0)
        Select case(int)
         Case(1);     ip1=1; ip2=2; jp1=1; jp2=2            !  I(ip1,jp1;ip2,jp1)   ( i . i .)
         Case(2);     ip1=1; ip2=2; jp1=1; jp2=2            !  I(jp1,ip1;jp2,ip2)   ( . i . i)
         Case(3);     ip1=1; ip2=2; jp1=2; jp2=1            !  I(ip1,jp2;jp1,ip2)   ( i . . i)
         Case(4);     ip1=2; ip2=1; jp1=1; jp2=2            !  I(jp1,ip2;ip1,jp2)   ( . i i .)
        End Select
       Case(1)
        Select case(int)
         Case(1);     ip1=1; ip2=2; jp1=2; jp2=1            !  I(ip1,jp1;ip2,jp2)   ( i . i .)
         Case(2);     ip1=2; ip2=1; jp1=1; jp2=2            !  I(jp1,ip1;jp2,ip2)   ( . i . i)
         Case(3);     ip1=1; ip2=1; jp1=2; jp2=2            !  I(ip1,jp2;jp1,ip2)   ( i . . i)
         Case(4);     ip1=2; ip2=2; jp1=1; jp2=1            !  I(jp1,ip2;ip1,jp2)   ( . i i .)
        End Select
      End Select
!----------------------------------------------------------------------------------------------

      xx = 0.d0; met=0
      Do ik = 1,nrk;  if(kk.ne.kr1(ik)) Cycle

       if(i.eq.kr3(ik)) then;          jnt = kr2(ik);       j = kr4(ik)
       elseif(i.eq.kr4(ik)) then;      jnt = kint(kr2(ik)); j = kr3(ik)
       else;                           Cycle
       end if
       if(jnt.ne.int) Cycle

       Call Density(ns,ks,x,p(1,jp1,j),p(1,jp2,j),sym)
       xx = xx + crk(ik)*x; met=met+1
      End do

      if(met.gt.0) Call Convol_sk(ns,ks,x,xx,int,sym,sym)
      if(met.gt.0) Call Update_hs(ms,hfm,ip1,ip2,ns,ks,x,sym)

      End do ! over int
      End do; End do ! over k and itype

      t2 = RRTC()
      time_hf_matrix_breit = time_hf_matrix_breit + t2-t1

      End Subroutine hf_matrix_breit


