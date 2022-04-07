!======================================================================
      MODULE dbsr_hf
!======================================================================
!     main parameters for the DBSR_HF program
!----------------------------------------------------------------------
      Use zconst

      Implicit none

! ... input/output files and units:

      Integer, parameter :: ma = 80   ! limit for file-name length
      Character(ma) :: AF

      Integer :: inp = 5;   Character(ma) :: AF_dat  = 'name.inp'
                            Character(ma) :: BF_dat  = '.inp'
      Integer :: log = 3;   Character(ma) :: AF_log  = 'name.log'
                            Character(ma) :: BF_log  = '.log'
      Integer :: scr = 0
      Integer :: nuw = 11;  Character(ma) :: AF_inp  = 'name.bsw'
                            Character(ma) :: BF_inp  = '.bsw'
                            Character(ma) :: AF_out  = 'name.bsw'
                            Character(ma) :: BF_out  = '.bsw'
                            Character(ma) :: AF_rwf  = 'name.w'
                            Character(ma) :: BF_rwf  = '.w'
                            Character(ma) :: AF_nl   = 'name.nl'
                            Character(ma) :: BF_nl   = '.nl'
      Integer :: nup = 12;  Character(ma) :: AF_plt  = 'name.plot'
                            Character(ma) :: BF_plt  = '.plot'
      Integer :: nuc = 13;  Character(ma) :: AF_cfg  = 'name.c'
                            Character(ma) :: BF_cfg  = '.c'
                            Character(ma) :: AF_conf = 'name.conf'
                            Character(ma) :: BF_conf = '.conf'
      Integer :: nus = 14;  Character(ma) :: AF_LS   = 'name.LS'
                            Character(ma) :: BF_LS   = '.LS'
      Integer :: nug = 15;  Character(ma) :: AF_grid = 'name.knot'
                            Character(ma) :: BF_grid = '.knot'
      Integer :: nua = 21;  ! scratch file
! ... name of case:

      Character(ma) :: name = ' '
      Character(ma) :: knot = 'knot.dat'

! ... atomic parameters:

      Real(8) :: z = 0.d0
      Integer :: an = 0, ai = 0
      Real(8) :: atw = 0.d0
      Real(8) :: rms = 0.d0
      Real(8) :: Etotal, E1body, E2body, E_breit,E_self,E_vacpol

      Character(6)   :: atom = ' '
      Character(6)   :: ion  = ' '
      Character(160) :: configuration = ' ', conf_AV = ' ', conf_LS = ' '
      Character(6)   :: term = 'LS'
      Character(80)  :: anit = 'all'

      Integer :: nelc =  0
      Integer :: nconf = 0
      Integer :: eal = 5
      Real(8), allocatable :: weight(:)
      Integer, allocatable :: iqconf(:,:)

! ... convergence:

      Real(8) :: scf_tol=1.d-10,  scf_diff
      Real(8) :: orb_tol=1.d-07,  orb_diff
      Real(8) :: end_tol=1.d-07
      Real(8) :: eps_ovl=1.d-07
      Integer :: max_it = 75

! ... debug options:

      Integer :: debug  = 0
      Integer :: newton = 0
      Integer :: rotate = 0

! ... core

      Integer,parameter :: mcore = 50
      Character(250) :: core=' '
      Integer :: ncore = 0
      Integer :: n_core(mcore), k_core(mcore),l_core(mcore),j_core(mcore)
      Character(5) :: e_core(mcore)

! ... description of 1 conf.w.function:

      Integer, parameter :: msh = 31 ! max. number of shells behind core
      Integer :: no
      Integer, dimension(msh) :: nn,kn,ln,jn,iq,in,Jshell,Vshell,Jintra

! ... Storing configuration as character strings:

      Character(9*msh+9) :: CONFIG, SHELLJ, INTRAJ

! ... orbital variables that depend on nwf

      Integer :: nwf = 0
      Integer :: nit = 0
      Integer :: ilzero = 1
      Integer :: ibzero = 2
      Integer :: kmin = 0
      Integer :: kmax = 0
      Integer :: jmin = -1
      Integer :: jmax = -1

! ... output of Breit corrections:

      Integer :: mbreit =  0
      Integer :: mode_SE =  3
      Integer :: mode_VP =  1
      Real(8) :: eps_c = 1.d-6
      Integer :: ibi = 2**16

! ... solutions for Rydberg series:

      Integer :: out_nl = 0
      Integer :: nsol_nl = 0
      Real(8), allocatable :: p_nl(:,:), e_nl(:)

      Integer :: out_w = 0     ! output in the GRASP format
      Integer :: out_plot = 0  ! output in table form

! ... frequently called functions: (instead interface)

      Integer, external :: Icheck_file

! ... debuging time:

      Real(8) :: time_hf_eiv=0.d0, time_hf_matrix=0.d0, &
                 time_hf_matrix_breit=0.d0, time_update_int=0.d0
      Real(8) :: au_cm, au_eV

      End Module dbsr_hf

