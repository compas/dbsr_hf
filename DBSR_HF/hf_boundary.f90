!====================================================================
      Subroutine Boundary_conditions
!====================================================================
! ... apply zero conditions at r=0 and on boundary r=a
! ... iprm(i) = 0 means deletion of B-spline "i"
!--------------------------------------------------------------------
      Use dbsr_hf
      Use DBS_nuclear
      Use DBS_grid
      Use df_orbitals

      Implicit none
      Integer :: i,j

      iprm=1
      if(nuclear.eq.'point'.or.ilzero.eq.0) then
       Do i=1,nbf
        iprm(1,i)=0;    j=nsp-ibzero+1; iprm(j:ns,i)=0
        iprm(ns+1,i)=0; j=nsq-ibzero+1; iprm(j+ns:ms,i)=0
       End do
      else
       Do i=1,nbf
        j=lbs(i)+1; if(j.gt.ksp-1) j=1; iprm(1:j,i)=0
        j=nsp-ibzero+1; iprm(j:ns,i)=0
        if(kbs(i).lt.0) j=lbs(i)+2
        if(kbs(i).ge.1) j=lbs(i)
        if(j.gt.ksq-1) j=1;  iprm(ns+1:ns+j,i)=0
        j=nsq-ibzero+1; iprm(j+ns:ms,i)=0
       End do
      end if

      End Subroutine Boundary_conditions


