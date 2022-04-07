!======================================================================
      Subroutine get_estimates
!======================================================================
!     Get initial estimates:
!     (1) read from bsw.inp
!     (2) screened hydrogenic
!----------------------------------------------------------------------
      Use DBS_grid
      Use DBSR_HF
      Use df_orbitals

      Implicit none
      Real(8) :: zz,ss,s
      Real(8), external :: QUADR
      Integer :: i, j,k, nnl
      Integer, allocatable :: nl(:)
      Real(8), allocatable :: qnl(:), snl(:)
      Integer, external :: ipointer

      write(log,'(//a/)') 'Initial estimations for orbitals:'

! ... Read AF_inp for the fixed orbitals or initial estimates

      AF_inp=trim(name)//BF_inp
      Call Read_apar(inp,'inp',AF_inp)
      Call Read_aarg('inp',AF_inp)

      if(Icheck_file(AF_inp).ne.0) then
       open(nuw,file=AF_inp,form='UNFORMATTED')
       i = INDEX(AF_inp,'.',BACK=.TRUE.)
       if(AF_inp(i:).eq.'.w')  Call Read_GRASP(nuw)
       if(AF_inp(i:).eq.'.bsw')  Call Read_pqbs(nuw)
      end if


      Allocate(nl(nwf),qnl(nwf),snl(nwf))

      nl=0; nnl = 0
      Do i = 1, nbf
       j = ipointer(nwf,nl,nbs(i)*1000+lbs(i))
       if(j.gt.0) then
        qnl(j) = qnl(j) + qsum(i)
       else
        nnl=nnl+1; nl(nnl) =  nbs(i)*1000+lbs(i)
                   qnl(nnl) = qsum(i)
       end if
      End do

      ss = 0.d0
      Do i = 1,nnl
       s =  max(0.d0,qnl(i)-1.d0)
       snl(i) = ss + s
       ss = ss + qnl(i)
      End do


      Do i = 1, nbf
       j = ipointer(nnl,nl,nbs(i)*1000+lbs(i))

       ss = snl(j); zz = z-ss

       if(mbs(i) /= 0) Cycle     ! We have an initial estimate

       Call bdcwf_pq(nbs(i),kbs(i),zz,p(1,1,i),p(1,2,i))
       mbs(i) = ns - 1
       Do j=nsp-1,nsp/3,-1
        if(abs(p(j,1,i)).gt.end_tol/1000) Exit   ! ???
        mbs(i) = j - 1
       End do

! ... set orthogonality constraints:

       Do j = 1,i-1
        if (kbs(i) /= kbs(j)) Cycle
        S=QUADR(p(1,1,i),p(1,1,j),0)
        if(abs(S).lt.Eps_ovl) Cycle
        p(:,:,i) = p(:,:,i) - s * p(:,:,j)
       End do
       S=QUADR(p(1,1,i),p(1,1,i),0)

       p(:,:,i) = p(:,:,i)/ sqrt(S)

       write(log,'(a,a,f5.2)') ebs(i),&
       ' - hydrogenic orbital with screening ',ss

      End do  ! over orbitals

      Call Check_tails(0)

      Call update_int(0)
      Call Energy

      End Subroutine get_estimates


