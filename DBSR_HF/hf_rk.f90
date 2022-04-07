!======================================================================
      Real(8) Function hf_rk (i1,j1,i2,j2,k)
!======================================================================
!     Returns  hf_rk (i1, j1; i2, j2), base on the assembling two-electron
!     B-spline integrals (see module DBS_integral)
!----------------------------------------------------------------------
      Use DBS_grid,      only: ns,ks
      Use DF_orbitals,   only: p

      Implicit none
      Integer, intent(in) :: i1,j1,i2,j2,k
      Real(8) :: dens1(ns,ks),dens2(ns,ks),dens3(ns,ks),dens4(ns,ks), &
                 conv(ns,ks)
      Real(8), external :: SUM_AmB

      hf_rk = 0.d0

      Call density (ns,ks,dens1,p(1,1,i1),p(1,1,i2),'s')
      Call density (ns,ks,dens2,p(1,1,j1),p(1,1,j2),'s')
      Call density (ns,ks,dens3,p(1,2,i1),p(1,2,i2),'s')
      Call density (ns,ks,dens4,p(1,2,j1),p(1,2,j2),'s')

      Call mrk_pppp(k)
      Call convol  (ns,ks,conv,dens1,2,'s','s')
      hf_rk = hf_rk + SUM_AmB(ns,ks,conv,dens2,'s')

      Call mrk_qqqq(k)
      Call convol  (ns,ks,conv,dens3,2,'s','s')
      hf_rk = hf_rk + SUM_AmB(ns,ks,conv,dens4,'s')

      Call mrk_qpqp(k)
      Call convol  (ns,ks,conv,dens3,2,'s','s')
      hf_rk = hf_rk + SUM_AmB(ns,ks,conv,dens2,'s')

      Call mrk_pqpq(k)
      Call convol  (ns,ks,conv,dens1,2,'s','s')
      hf_rk = hf_rk + SUM_AmB(ns,ks,conv,dens4,'s')

      End Function hf_rk


!======================================================================
      Real(8) Function hf_fk (i,j,k)                  result(S)
!======================================================================
!     Returns  Fk (i, j) base on the assembling two-electron
!     B-spline integrals (see module DBS_integral)
!----------------------------------------------------------------------
      Use DBS_grid,      only: ns,ks
      Use df_orbitals,   only: p

      Implicit none
      Integer, intent(in) :: i,j,k
      Real(8) :: ppi(ns,ks),qqi(ns,ks),ppj(ns,ks),qqj(ns,ks), &
                 pppp(ns,ks),qqqq(ns,ks),qpqp(ns,ks),pqpq(ns,ks)
      Real(8), external :: SUM_AmB

      Call density (ns,ks,ppi,p(1,1,i),p(1,1,i),'s')
      Call density (ns,ks,ppj,p(1,1,j),p(1,1,j),'s')
      Call density (ns,ks,qqi,p(1,2,i),p(1,2,i),'s')
      Call density (ns,ks,qqj,p(1,2,j),p(1,2,j),'s')

      S = 0.d0

      Call mrk_pppp(k)
      Call convol (ns,ks,pppp,ppi,2,'s','s')
      S = S + SUM_AmB(ns,ks,pppp,ppj,'s')

      Call mrk_qqqq(k)
      Call convol (ns,ks,qqqq,qqi,2,'s','s')
      S = S + SUM_AmB(ns,ks,qqqq,qqj,'s')

      Call mrk_qpqp(k)
      Call convol (ns,ks,qpqp,qqi,2,'s','s')
      S = S + SUM_AmB(ns,ks,qpqp,ppj,'s')

      Call mrk_pqpq(k)
      Call convol (ns,ks,pqpq,ppi,2,'s','s')
      S = S + SUM_AmB(ns,ks,pqpq,qqj,'s')

      End Function hf_fk


!======================================================================
      Real(8) Function hf_gk (i,j,k)                  result(S)
!======================================================================
!     Exchange integral:  Rk(i,j;j,i)
!----------------------------------------------------------------------
      Use DBS_grid,      only: ns,ks
      Use df_orbitals,   only: p

      Implicit none
      Integer, intent(in) :: i,j,k
      Real(8) :: pp(ns,ks),qq(ns,ks),conv(ns,ks)
      Real(8), external :: SUM_AmB

      S = 0.d0

      Call density (ns,ks,pp,p(1,1,i),p(1,1,j),'s')
      Call density (ns,ks,qq,p(1,2,i),p(1,2,j),'s')

      Call mrk_pppp(k)
      Call convol (ns,ks,conv,pp,2,'s','s')
      S = S + SUM_AmB(ns,ks,conv,pp,'s')

      Call mrk_qqqq(k)
      Call convol (ns,ks,conv,qq,2,'s','s')
      S = S + SUM_AmB(ns,ks,conv,qq,'s')

      Call mrk_qpqp(k)
      Call convol (ns,ks,conv,qq,2,'s','s')
      S = S + 2.d0*SUM_AmB(ns,ks,conv,pp,'s')

      End Function hf_gk
