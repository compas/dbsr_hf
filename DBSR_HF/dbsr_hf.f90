!======================================================================
!     Program DBSR_HF
!======================================================================
!                   C O P Y R I G H T -- 2015
!     Written by:   Oleg Zatsarinny and Charlotte Froese Fischer
!     email:        oleg_zoi@yahoo.com
!----------------------------------------------------------------------
!     This program computes the radial functions for Dirac-Fock
!     problem in B-spline basis. In addition to standard one-electron
!     case, the program allows also simultaineous calculation of
!     several configurations in different optimization schemes.
!     Another possibility - calculation whole Rydberg series in
!     frozen-core approximation.   
!
!     For short instructions, run dbsr_hf ? and look in name.inp.
!----------------------------------------------------------------------
      Use dbsr_hf
      Use DBS_debug

      Implicit none
      Real(8) :: t1,t2,t3,t4

      t1=RRTC()

! ... atomic paramters:

      Call Get_case

! ... set up the energy expression:

      Call Def_energy_coef

! ... B-spline parameters:

      Call Get_spline_param

! ... computing:

      t3=RRTC()

      Call Get_estimates

      t4=RRTC()
      if(debug.gt.0) &
      write(scr,'(a,T30,f10.2,a)') 'Get_estimates:',t4-t3,'  sec' 

      t3=RRTC()

      Call Solve_HF

      t4=RRTC()
      if(debug.gt.0) &
      write(scr,'(a,T30,f10.2,a)') 'Solve_HF:',t4-t3,'  sec' 

! ... output of results:

      if(term.ne.'jj') then
       Call record_confs
       AF_conf = trim(name)//trim(BF_conf)
       AF_cfg = trim(name)//trim(BF_cfg)
       close(nuc)
       Call gen_jj_states(AF_conf,AF_cfg,jmin,jmax)
      end if

      t3=RRTC()
      Call Write_pqbs
      Call Write_nl
      Call Plot_bsw
      Call Grasp_wfn
      Call Summry
      t4=RRTC()
      if(debug.gt.0) &
      write(scr,'(a,T30,f10.2,a)') 'Summry:',t4-t3,'  sec' 

      t2=RRTC()
      write(scr,'(/a,f10.2,a)') 'time:',t2-t1,'  sec' 
      write(log,'(/a,f10.2,a)') 'time:',t2-t1,'  sec' 

      if(debug.gt.0) &
      write(scr,'(a,T25,f10.2,a,i10)') 'Convol:',time_convol,'  sec',ic_convol 
      if(debug.gt.0) &
      write(scr,'(a,T25,f10.2,a,i10)') 'Density:',time_density,'  sec',ic_density 

      End ! Program DBSR_HF

