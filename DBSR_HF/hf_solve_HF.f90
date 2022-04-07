!=====================================================================
      Subroutine solve_HF
!=====================================================================
!     Solve the DF equations in turns
!---------------------------------------------------------------------
      Use DBS_grid
      Use dbsr_hf
      Use df_orbitals

      Implicit none
      Real(8) :: hfm(ms,ms), v(ms), et
      Real(8), external :: QUADR
      Integer :: ip, i,j, it

      Real(8) :: t1,t2,t3,t4 , S1,S2
      Real(8), external :: RRTC

      t1 = RRTC()

      et = etotal
      dpm = 0.d0
      Do it = 1,max_it
       write(log,'(//A,I6/A/)') 'Iteration ',it, &
                                '----------------'
       write(log,'(2x,a,9x,a,10x,a,9x,a/)')  'nl', 'e(nl)', 'dpm', 'ns'

! ... main iterations other orbitals:

       Do ip = 1,nwf; i=iord(ip); if(i.eq.0) Cycle

        Do j=i+1,nwf
         if(it.le.2) Cycle; if(rotate.eq.0) Cycle;  Call Rotate_ij(i,j)
        End do

        Call hf_matrix(i,hfm)

        ! .. diagonalize the hf matrix

        Call hf_eiv(i,hfm,v)


        dpm(i)=maxval(abs(p(1:ns,1,i)-v(1:ns)))/maxval(abs(p(:,1,i)))
        if(ip.eq.1) then
         if(dpm(i).lt.orb_tol) iord(i)=0
        else
         if(dpm(i).lt.orb_tol.and.iord(ip-1).eq.0) iord(i)=0
        end if

        write(log,'(A5,f16.8,1PD13.2,I8)')  ebs(i),e(i,i),dpm(i),mbs(i)


        Call Put_pv(i,v)

        Call Check_tails(i)

        t3=RRTC()
        Call Update_int(i)
        t4=RRTC()
        time_update_int=time_update_int+t4-t3

        ! .. remove tail zero

        if(debug.gt.0) then
         Do j = 1,nwf
          if(i.eq.j) Cycle
          if(e(i,j) < 1.d-10) Cycle
          if(kbs(i).ne.kbs(j)) Cycle
          write(log,'(a,a,a,f16.8)') &
           'Orthogonality ',ebs(i),ebs(j),QUADR(p(1,1,i),p(1,1,j),0)
         End do
        end if

       End do ! over functions

       Call Energy
       write(log,'(/a,F16.8)') 'Etotal = ',etotal

       orb_diff =  maxval(abs(dpm(1:nwf)))
       scf_diff =  abs(et-etotal)/abs(etotal)
       et = etotal
       write(log,*)
       write(log,'(A,T25,1PD10.2)') 'Maximum orbital diff. ', orb_diff
       write(log,'(A,T25,1PD10.2)') 'Energy difference     ', scf_diff
       write(scr,'(a,F20.12,i5,1P2E10.2)') &
        'etotal,it,E_diff,orb_diff = ',etotal,it,scf_diff,orb_diff

       Call Boundary_conditions

       if ( orb_diff < orb_tol .and. scf_diff  < scf_tol ) Exit

      End do ! over iterations

      if(debug.gt.0) &
      write(scr,'(a,T30,f10.2,a)') 'time_eiv:',time_hf_eiv,'  sec'
      if(debug.gt.0) &
      write(scr,'(a,T30,f10.2,a)') 'time_mat:',time_hf_matrix,'  sec'
      if(debug.gt.0) &
      write(scr,'(a,T30,f10.2,a)') 'time_int:',time_update_int,'  sec'

      End Subroutine solve_HF


!==================================================================
      Subroutine hf_eiv(i,hfm,v)
!==================================================================
!     Find the eigenvector v of hfm for the "positive-eneregy"
!     eigenvalue m=n-l. Each orthogonality condition to the lower
!     state reduces m - m - 1.  We supposed that nl orbitals are
!     ordered to their n-values.
!------------------------------------------------------------------
      Use DBS_grid
      Use DBS_gauss
      Use df_orbitals
      Use dbsr_hf
      Use zconst,  only: c_au

      Implicit none
      Integer, intent(in)    :: i
      Real(8), intent(inout) :: hfm(ms,ms)
      Real(8), intent(out)   :: v(ms)
      Real(8) :: aa(ms,ms),ss(ms,ms),w(3*ms),eval(ms),a(ms),s(ms),zz
      Integer :: j, jp, info, k,kk,m,mm, ipos(1)

      Real(8) :: t1,t2
      Real(8), external :: RRTC

      t1 = RRTC()

! ... apply orthogonality conditions for orbitals

      m = nbs(i)-lbs(i)
      Do j = 1,nwf
       if(i.eq.j) Cycle
       if(e(i,j) < 1.d-10) Cycle
       if(kbs(i).ne.kbs(j)) Cycle
       Call orthogonality(hfm,p(1,1,j))
       if(j.lt.i) m = m - 1
       if(debug.gt.0) write(log,'(a,a,a,a,i3)') &
        'Orthogonality ',ebs(i),ebs(j),' is applied, m =',m
      End do

! ... apply boundary conditions (delete extra B-splines)

      kk=0
      Do j=1,ms
       if(iprm(j,i).eq.0) Cycle; kk=kk+1
       k=0
       Do jp=1,ms
        if(iprm(jp,i).eq.0) Cycle
        k=k+1; a(k)=hfm(jp,j); s(k)=fppqq(jp,j)
       End do
       aa(1:k,kk)=a(1:k); ss(1:k,kk)=s(1:k)
      End do

! ... evaluates the eigenvalues and eigenvectors (LAPACK routine):

      Call dsygv(1,'V','L',kk,aa,ms,ss,ms,eval,w,3*ms,INFO)
      if (info /= 0) then
       write(scr,'(a,i6)') 'Error in Eigenvalue routine, dsygv', info
       Stop ' '
      end if

      mm=0; zz = -1.999*c_au**2
      Do j=1,kk; mm = j; if(eval(j).gt.zz) Exit; End do
      mm = m + mm - 1

! ... save all solutions if nl > 0:

      if(out_nl.gt.0.and.i.eq.nwf) then
       nsol_nl = kk - mm + 1
       if(.not.allocated(p_nl)) Allocate(p_nl(ms,nsol_nl),e_nl(nsol_nl))
       p_nl = 0.d0; e_nl=0.d0
       Do m=mm,kk
        a(1:ms) = aa(1:ms,m);  v=0.d0; k=0
        Do j=1,ms
         if(iprm(j,i).eq.0) Cycle; k=k+1; v(j)=a(k)
        End do

        if (v(ks) < 0.d0) v = -v

        p_nl(:,m-mm+1)=v(:)
        e_nl(m-mm+1) = eval(m)
       End do
      end if

! ... restore the solutions in original B-spline net:

      a(1:ms) = aa(1:ms,mm);  v=0.d0; k=0
      Do j=1,ms
       if(iprm(j,i).eq.0) Cycle; k=k+1; v(j)=a(k)
      End do

      ipos=maxloc(abs(v))
      if(v(ipos(1)).lt.0.d0) v=-v

!      if (v(ks) < 0.d0) v = -v


      e(i,i) = eval(mm)

      if(debug.gt.0) then
       write(log,'(a,5E15.5)') 'eval =',eval(mm-2:mm+2)
       write(log,'(a,2i5,E15.5)') 'we choose m,mm,e =',m,mm,eval(mm)
      end if

      t2 = RRTC()
      time_hf_eiv = time_hf_eiv + t2-t1

      End Subroutine hf_eiv

