!======================================================================
      Subroutine rotate_ij(i,j)
!======================================================================
!     Determine the rotation parameter for a stationary solution for
!     the (i,j) orbitals, i < j
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBSR_hf
      Use df_orbitals, l=>lbs

      Implicit none
      Integer, intent(in) :: i,j
      Integer :: m,k
      Real(8) :: v(ns),w(ns),eps,dn,g,dg,  &
       Hii,Hij,Hjj, Fkii,Fkij,Fkjj,Gkij, RKij,RKji, FKim,Fkjm, GKim,Gkjm
      Real(8), external :: Vp_dhl, a, b, rk

      if(i.ge.j) Return
      if(kbs(i).ne.kbs(j)) Return
      if(iord(i).eq.0) Return
      if(iord(j).eq.0) Return
      if(clsd(i).and.clsd(j)) Return

      Hii = Vp_dhl(kbs(i),p(1,1,i),kbs(i),p(1,1,i))
      Hij = Vp_dhl(kbs(i),p(1,1,i),kbs(j),p(1,1,j))
      Hjj = Vp_dhl(kbs(j),p(1,1,j),kbs(j),p(1,1,j))
      G   = 2.d0*(qsum(j)-qsum(i))*Hij
      dG  = (qsum(i)-qsum(j))*(Hjj-Hii)

      Do k= 0,2*min0(l(i),l(j)),2
       RKij = rk(p(1,1,i),p(1,1,i),p(1,1,i),p(1,1,j),k)
       RKji = rk(p(1,1,i),p(1,1,j),p(1,1,j),p(1,1,j),k)
       G = G + 2.d0 * (a(i,j,k)+b(i,j,k)-2*a(i,i,k)) * RKij
       G = G - 2.d0 * (a(i,j,k)+b(i,j,k)-2*a(j,j,k)) * RKji
       FKii = rk(p(1,1,i),p(1,1,i),p(1,1,i),p(1,1,i),k)
       FKij = rk(p(1,1,i),p(1,1,j),p(1,1,i),p(1,1,j),k)
       FKjj = rk(p(1,1,j),p(1,1,j),p(1,1,j),p(1,1,j),k)
       GKij = rk(p(1,1,i),p(1,1,j),p(1,1,j),p(1,1,i),k)
       dG = dG + 2*a(i,i,k)*(FKij-FKii+2*GKij) + 2*a(j,j,k)*(FKij-FKjj+2*GKij)
       dG = dG + (a(i,j,k)+b(i,j,k)) * (FKii+Fkjj-2*FKij-4*GKij)
      End do

      Do m = 1,nbf
       if (m ==i .or. m == j) Cycle
       Do  k = 0,2*min0(l(i),l(j)),2
        if (a(i,m,k) == a(j,m,k)) Cycle
        g  = g + 2.d0*(a(j,m,k)-a(i,m,k))*rk(p(1,1,i),p(1,1,m),p(1,1,j),p(1,1,m),k)
        FKim = rk(p(1,1,i),p(1,1,m),p(1,1,i),p(1,1,m),k)
        FKjm = rk(p(1,1,j),p(1,1,m),p(1,1,j),p(1,1,m),k)
        dg =dg + (a(j,m,k)-a(i,m,k))*(FKim-FKjm)
       End do
       Do k = abs(l(i)-l(m)),l(i)+l(m),2
        if (b(i,m,k) == b(j,m,k)) Cycle
        g  = g + 2.d0*(b(j,m,k)-b(i,m,k))*rk(p(1,1,i),p(1,1,m),p(1,1,m),p(1,1,j),k)
        GKim = rk(p(1,1,i),p(1,1,m),p(1,1,m),p(1,1,i),k)
        GKjm = rk(p(1,1,j),p(1,1,m),p(1,1,m),p(1,1,j),k)
        dg =dg + (b(j,m,k)-b(i,m,k))*(GKim-GKjm)
       End do
      End do

      if ( abs(g) < 1.d-12 .and. abs(dg) < 1.d-12 ) then

       e(i,j) = 0.d0
       e(j,i) = 0.d0

      else

       eps = -g / (2*dg)
       dn = 1.d0/sqrt(1+eps*eps)
       v = (p(:,1,i)-eps*p(:,1,j))*dn
       w = (p(:,1,j)+eps*p(:,1,i))*dn
       p(:,1,i) = v
       p(:,1,j) = w
       v = (p(:,2,i)-eps*p(:,2,j))*dn
       w = (p(:,2,j)+eps*p(:,2,i))*dn
       p(:,2,i) = v
       p(:,2,j) = w

       Call Check_tails(i)
       Call Check_tails(j)

       if(debug.gt.0)  write(log,'(/a,a,a,a,3f16.8/)') &
        'Rotate ',ebs(i),ebs(j),' eps,g,dg = ', eps,g,dg

       e(i,j) = 1.d0
       e(j,i) = 1.d0

      end if

      End Subroutine rotate_ij




