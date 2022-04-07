!======================================================================
      Subroutine Check_tails(jo)
!======================================================================
!     nulify too small coefficients in the end for orbital jo
!     if jo=0 - check all orbitals
!     only large component is used to define small B-splines
!----------------------------------------------------------------------
      Use dbsr_hf
      Use DBS_grid
      Use df_orbitals

      Implicit none
      Integer :: io,jo,i
      Real(8) :: cm

      Do io=1,nbf; if(jo.ne.0.and.io.ne.jo) Cycle
       cm = maxval( abs( p(:,1,io) ) )
       Do i=nsp-1,1,-1
        mbs(io)=i
        if(abs(p(i,1,io))/cm.lt.end_tol) Cycle
        Exit
       End do
       p(i+1:ns,1,io)=0.d0
       p(i+1:ns,2,io)=0.d0
      End do
      Call Boundary_conditions

      End Subroutine Check_tails

