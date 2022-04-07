!---------------------------------------------------------------
!     ZCOM LIBRARY
!---------------------------------------------------------------
!     contains common routines frequently used by other programs
!---------------------------------------------------------------
!     MODULES:
!
!     MOD_zconst.f90     -  physical constants
!     MOD_atoms_par.f90  -  atomic parameters and quantities
!---------------------------------------------------------------
!     ROUTINES:
!
!     ANG_dj_fact.f90
!     ANG_itra.f90
!     ANG_kappa.f90
!     ANG_t_ls_jk.f90
!     ANG_zcklm.f90
!     ANG_zclkl.f90
!     ANG_z_3j.f90
!     ANG_z_6j.f90
!     ANG_z_9j.f90
!     check_file.f90
!     clean_a.f90
!     full_mat.f90
!     FUN_alegf.f90
!     FUN_cgamma.f90
!     FUN_dirac.f90
!     gauss_LA.f90
!     gauss_points.f90
!     interv.f90
!     ipointer.f90
!     read_arg.f90
!     read_par.f90
!     rrtc.f90
!     splin3.f90
!---------------------------------------------------------------



!=======================================================================
      MODULE zconst
!=======================================================================
!     contains basic processor-dependent constants as well as
!     a larger variety of physical constants
!-----------------------------------------------------------------------

! ... kind-constants

      integer, parameter :: int4 = selected_int_kind(9)
      integer, parameter :: int2 = selected_int_kind(4)
      integer, parameter :: int1 = selected_int_kind(2)

      integer, parameter :: sp = kind(1.0)
      integer, parameter :: dp = selected_real_kind(2*precision(1.0_sp))
      integer, parameter :: qp = selected_real_kind(2*precision(1.0_dp))

! ... simple constants

      real(dp), parameter :: zero  = 0.0_dp,    &
                             one   = 1.0_dp,    &
                             two   = 2.0_dp,    &
                             three = 3.0_dp,    &
                             four  = 4.0_dp,    &
                             five  = 5.0_dp,    &
                             six   = 6.0_dp,    &
                             seven = 7.0_dp,    &
                             eight = 8.0_dp,    &
                             nine  = 9.0_dp,    &
                             ten   = 10.0_dp,   &
                             half  = one/two,   &
                             third = one/three

! ... some useful physical constants
! ... (based on NIST compilation 2001)

      real(dp), parameter ::                            &
       bohr_radius_in_cm      =    0.5291772083e-8_dp,  &
       ao_cm                  =    0.5291772083e-8_dp,  &
       one_over_alpha         =    137.03599976e+0_dp,  &
       alpha                  =    7.297352533e-3_dp,   &
       c_au                   =    137.03599976e+0_dp,  &
       c_vacuum_in_m_per_s    =    299792458e+0_dp,     &
       electron_charge_in_C   =    1.602176462e-19_dp,  &
       electron_mass_in_g     =    9.10938188e-28_dp,   &
       electron_mass_in_amu   =    5.485799110e-4_dp,   &
       proton_mass_in_amu     =    1.00727646688e+0_dp, &
       proton_electron        =    1836.1526675e+0_dp,  &
       electron_proton        =    5.446170232e-4_dp,   &
       hbar_in_J_s            =    1.054571596e-27_dp,  &
       hbar_in_eV_s           =    6.58211889e-16_dp,   &
       au_eV_inf              =    27.2113834e+0_dp,    &
       au_cm_inf              =    219471.62e+0_dp,     &
       Ry_eV_inf              =    13.60569175e+0_dp,   &
       Ry_cm_inf              =    109737.31568549_dp,  &
       time_au_in_sec         =    4.134138e+16_dp,     &
       time_au                =    2.4189e-17_dp,       &
       pi_a0_2                =    0.87973e+0_dp,       &
       fermi_in_cm            =    1.0e-13_dp,          &
       pi                     =    3.141592653589793238462643e+0_dp

! ... Standard input, output, error and scratch unit numbers

      integer(int4) :: istdi=5, istdo=6 ,istde=0

! ... maximum length for record:

      integer(int4), parameter ::  mrecl = 100000

! ... packing bases:

      Integer(int4), parameter :: ibf = 2**4  ! overlap factors
      Integer(int4), parameter :: ibd = 2**15 ! overlap determinants

      Integer(int4), parameter :: jb = 10     ! max. number of shells
      Integer(int4), parameter :: jb2=jb**2, jb4=jb**4, jb8=jb**8

      INTEGER(int4) :: ib1 = 2,     ib2 = 2**2,  ib3 = 2**3,  ib4 = 2**4, &
                    ib5 = 2**5,  ib6 = 2**6,  ib7 = 2**7,  ib8 = 2**8, &
                    ib9 = 2**9,  ib10= 2**10, ib15= 2**15, ib20= 2**20

      End module zconst


!======================================================================
      Module atoms_par
!======================================================================
!     atomic parameters according to periodic table
!----------------------------------------------------------------------
!     Nuclear radii are taken from
!     I. Angeli, K. P. Marinova, Atomic Data and Nuclear
!     Data Tables {\bf 99}, 69-95 (2013).
!     - - - - - - - - - - - - - - - - - - - - - - - - -
!     If no data for the element are available, we use
!     the empirical formula for Z < 91 from:
!     W.R. Johnson, G.Soff, Atomic Data and Nudear Data
!                   Tables {\bf 33}, p.405 (1985)
!     and for Z> 90 from:
!     I. Goidenko, Private communication.
!----------------------------------------------------------------------
      Implicit none
      Integer, parameter :: n_atoms = 104
      Type atomic
       INTEGER :: an
       CHARACTER(2) :: symbol
       CHARACTER(4) :: core
       CHARACTER(40) :: conf
       REAL(8) :: weight
       REAL(8) :: rrms
      End type atomic

      Type(atomic), dimension(n_atoms), parameter, public :: atoms = (/  &
       atomic(   1,  'H ',  '    ',  '1s(1)',                   1.d0, 0.8783d0), &
       atomic(   2,  'He',  '    ',  '1s(2)',                   4.d0, 1.6755d0), &
       atomic(   3,  'Li',  '[He]',  '2s(1)',                   6.d0, 2.5890d0), &
       atomic(   4,  'Be',  '[He]',  '2s(2)',                   9.d0, 2.5190d0), &
       atomic(   5,  'B ',  '[Be]',  '2p-(1)',                 10.d0, 2.4277d0), &
       atomic(   6,  'C ',  '[Be]',  '2p-(1)2p(1)',            12.d0, 2.4702d0), &
       atomic(   7,  'N ',  '[Be]',  '2p-(1)2p(2)',            14.d0, 2.5582d0), &
       atomic(   8,  'O ',  '[Be]',  '2p-(1)2p(3)',            16.d0, 2.6991d0), &
       atomic(   9,  'F ',  '[Be]',  '2p-(2)2p(3)',            19.d0, 2.8976d0), &
       atomic(  10,  'Ne',  '[Be]',  '2p-(2)2p(4)',            20.d0, 3.0055d0), &
       atomic(  11,  'Na',  '[Ne]',  '3s(1)',                  23.d0, 2.9936d0), &
       atomic(  12,  'Mg',  '[Ne]',  '3s(2)',                  24.d0, 3.0570d0), &
       atomic(  13,  'Al',  '[Mg]',  '3p-(1)',                 27.d0, 3.0610d0), &
       atomic(  14,  'Si',  '[Mg]',  '3p-(1)3p(1)',            28.d0, 3.1224d0), &
       atomic(  15,  'P ',  '[Mg]',  '3p-(1)3p(2)',            31.d0, 3.1889d0), &
       atomic(  16,  'S ',  '[Mg]',  '3p-(1)3p(3)',            32.d0, 3.2611d0), &
       atomic(  17,  'Cl',  '[Mg]',  '3p-(2)3p(3)',            35.d0, 3.3654d0), &
       atomic(  18,  'Ar',  '[Mg]',  '3p-(2)3p(4)',            38.d0, 3.4028d0), &
       atomic(  19,  'K ',  '[Ar]',  '4s(1)',                  39.d0, 3.4349d0), &
       atomic(  20,  'Ca',  '[Ar]',  '4s(2)',                  40.d0, 3.4776d0), &
       atomic(  21,  'Sc',  '[Ar]',  '3d-(1)4s(2)',            45.d0, 3.5459d0), &
       atomic(  22,  'Ti',  '[Ar]',  '3d-(2)4s(2)',            48.d0, 3.5921d0), &
       atomic(  23,  'V ',  '[Ar]',  '3d-(3)4s(2)',            51.d0, 3.6002d0), &
       atomic(  24,  'Cr',  '[Ar]',  '3d-(4)3d(1)4s(1)',       52.d0, 3.6452d0), &
       atomic(  25,  'Mn',  '[Ar]',  '3d-(4)3d(1)4s(2)',       55.d0, 3.7057d0), &
       atomic(  26,  'Fe',  '[Ar]',  '3d-(4)3d(2)4s(2)',       56.d0, 3.7377d0), &
       atomic(  27,  'Co',  '[Ar]',  '3d-(4)3d(3)4s(2)',       59.d0, 3.7875d0), &
       atomic(  28,  'Ni',  '[Ar]',  '3d-(4)3d(4)4s(2)',       60.d0, 3.8118d0), &
       atomic(  29,  'Cu',  '[Ar]',  '3d-(4)3d(6)4s(1)',       63.d0, 3.8823d0), &
       atomic(  30,  'Zn',  '[Ar]',  '3d-(4)3d(6)4s(2)',       66.d0, 3.9491d0), &
       atomic(  31,  'Ga',  '[Zn]',  '4p-(1)',                 69.d0, 3.9973d0), &
       atomic(  32,  'Ge',  '[Zn]',  '4p-(2)',                 74.d0, 4.0742d0), &
       atomic(  33,  'As',  '[Zn]',  '4p-(2)4p(1)',            75.d0, 4.0968d0), &
       atomic(  34,  'Se',  '[Zn]',  '4p-(2)4p(2)',            80.d0, 4.1400d0), &
       atomic(  35,  'Br',  '[Zn]',  '4p-(2)4p(3)',            79.d0, 4.1629d0), &
       atomic(  36,  'Kr',  '[Zn]',  '4p-(2)4p(4)',            86.d0, 4.1835d0), &
       atomic(  37,  'Rb',  '[Kr]',  '5s(1)',                  87.d0, 4.1989d0), &
       atomic(  38,  'Sr',  '[Kr]',  '5s(2)',                  88.d0, 4.2240d0), &
       atomic(  39,  'Y ',  '[Kr]',  '4d-(1)5s(2)',            89.d0, 4.2430d0), &
       atomic(  40,  'Zr',  '[Kr]',  '4d-(2)5s(2)',            90.d0, 4.2694d0), &
       atomic(  41,  'Nb',  '[Kr]',  '4d-(4)5s(1)',            93.d0, 4.3240d0), &
       atomic(  42,  'Mo',  '[Kr]',  '4d-(4)4d(1)5s(1)',       92.d0, 4.3151d0), &
       atomic(  43,  'Tc',  '[Kr]',  '4d-(4)4d(1)5s(2)',       97.d0, 0.0000d0), &  ! ???
       atomic(  44,  'Ru',  '[Kr]',  '4d-(4)4d(3)5s(1)',      104.d0, 4.5098d0), &
       atomic(  45,  'Rh',  '[Kr]',  '4d-(4)4d(4)5s(1)',      103.d0, 4.4945d0), &
       atomic(  46,  'Pd',  '[Kr]',  '4d-(4)4d(6)',           108.d0, 4.5563d0), &
       atomic(  47,  'Ag',  '[Kr]',  '4d-(4)4d(6)5s(1)',      109.d0, 4.5638d0), &
       atomic(  48,  'Cd',  '[Kr]',  '4d-(4)4d(6)5s(2)',      114.d0, 4.6087d0), &
       atomic(  49,  'In',  '[Cd]',  '5p-(1)',                115.d0, 4.6156d0), &
       atomic(  50,  'Sn',  '[Cd]',  '5p-(2)',                120.d0, 4.6519d0), &
       atomic(  51,  'Sb',  '[Cd]',  '5p-(2)5p(1)',           121.d0, 4.6802d0), &
       atomic(  52,  'Te',  '[Cd]',  '5p-(2)5p(2)',           130.d0, 4.7423d0), &
       atomic(  53,  'I ',  '[Cd]',  '5p-(2)5p(3)',           127.d0, 4.7500d0), &
       atomic(  54,  'Xe',  '[Cd]',  '5p-(2)5p(4)',           136.d0, 4.7964d0), &
       atomic(  55,  'Cs',  '[Xe]',  '6s(1)',                 133.d0, 4.8041d0), &
       atomic(  56,  'Ba',  '[Xe]',  '6s(2)',                 138.d0, 4.8378d0), &
       atomic(  57,  'La',  '[Xe]',  '5d-(1)6s(2)',           139.d0, 4.8550d0), &
       atomic(  58,  'Ce',  '[Xe]',  '4f-(1)5d-(1)6s(2)',     140.d0, 4.8771d0), &
       atomic(  59,  'Pr',  '[Xe]',  '4f-(3)6s(2)',           141.d0, 4.8919d0), &
       atomic(  60,  'Nd',  '[Xe]',  '4f-(4)6s(2)',           142.d0, 4.9123d0), &
       atomic(  61,  'Pm',  '[Xe]',  '4f-(5)6s(2)',           145.d0, 0.0000d0), &  !  ???
       atomic(  62,  'Sm',  '[Xe]',  '4f-(6)6s(2)',           144.d0, 4.9524d0), &
       atomic(  63,  'Eu',  '[Xe]',  '4f-(6)4f(1)6s(2)',      145.d0, 4.9663d0), &
       atomic(  64,  'Gd',  '[Xe]',  '4f-(6)4f(1)5d-(1)6s(2)',160.d0, 5.1734d0), &
       atomic(  65,  'Tb',  '[Xe]',  '4f-(6)4f(3)6s(2)',      159.d0, 5.0600d0), &
       atomic(  66,  'Dy',  '[Xe]',  '4f-(6)4f(4)6s(2)',      148.d0, 5.0455d0), &
       atomic(  67,  'Ho',  '[Xe]',  '4f-(6)4f(5)6s(2)',      165.d0, 5.2022d0), &
       atomic(  68,  'Er',  '[Xe]',  '4f-(6)4f(6)6s(2)',      170.d0, 5.2789d0), &
       atomic(  69,  'Tm',  '[Xe]',  '4f-(5)4f(8)6s(2)',      169.d0, 5.2256d0), &
       atomic(  70,  'Yb',  '[Xe]',  '4f-(6)4f(8)6s(2)',      176.d0, 5.3215d0), &
       atomic(  71,  'Lu',  '[4f]',  '5d-(1)6s(2)',           175.d0, 5.3700d0), &
       atomic(  72,  'Hf',  '[4f]',  '5d-(2)6s(2)',           178.d0, 5.3371d0), &
       atomic(  73,  'Ta',  '[4f]',  '5d-(3)6s(2)',           181.d0, 5.3507d0), &
       atomic(  74,  'W ',  '[4f]',  '5d-(4)6s(2)',           184.d0, 5.3658d0), &
       atomic(  75,  'Re',  '[4f]',  '5d-(4)5d(1)6s(2)',      185.d0, 5.3596d0), &
       atomic(  76,  'Os',  '[4f]',  '5d-(4)5d(2)6s(2)',      192.d0, 5.4126d0), &
       atomic(  77,  'Ir',  '[4f]',  '5d-(4)5d(3)6s(2)',      191.d0, 5.3968d0), &
       atomic(  78,  'Pt',  '[4f]',  '5d-(4)5d(5)6s(1)',      194.d0, 5.4236d0), &
       atomic(  79,  'Au',  '[4f]',  '5d-(4)5d(6)6s(1)',      197.d0, 5.4371d0), &
       atomic(  80,  'Hg',  '[4f]',  '5d-(4)5d(6)6s(2)',      198.d0, 5.4463d0), &
       atomic(  81,  'Tl',  '[Hg]',  '6p-(1)',                205.d0, 5.4759d0), &
       atomic(  82,  'Pb',  '[Hg]',  '6p-(2)',                208.d0, 5.5012d0), &
       atomic(  83,  'Bi',  '[Hg]',  '6p-(2)6p(1)',           209.d0, 5.5211d0), &
       atomic(  84,  'Po',  '[Hg]',  '6p-(2)6p(2)',           208.d0, 5.5584d0), &
       atomic(  85,  'At',  '[Hg]',  '6p-(2)6p(3)',           212.d0, 0.0000d0), &  ! ???
       atomic(  86,  'Rn',  '[Hg]',  '6p-(2)6p(4)',           212.d0, 5.5915d0), &
       atomic(  87,  'Fr',  '[Rn]',  '7s(1)',                 212.d0, 5.5915d0), &
       atomic(  88,  'Ra',  '[Rn]',  '7s(2)',                 214.d0, 5.6079d0), &
       atomic(  89,  'Ac',  '[Rn]',  '6d-(1)7s(2)',           227.d0, 0.0000d0), &  ! ???
       atomic(  90,  'Th',  '[Rn]',  '6d-(2)7s(2)',           232.d0, 5.7848d0), &
       atomic(  91,  'Pa',  '[Rn]',  '5f-(2)6d-(1)7s(2)',     231.d0, 0.0000d0), &  ! ???
       atomic(  92,  'U ',  '[Rn]',  '5f-(3)6d-(1)7s(2)',     238.d0, 5.8571d0), &
       atomic(  93,  'Np',  '[Rn]',  '5f-(4)6d-(1)7s(2)',     237.d0, 0.0000d0), &  ! ???
       atomic(  94,  'Pu',  '[Rn]',  '5f-(6)7s(2)',           239.d0, 5.8601d0), &
       atomic(  95,  'Am',  '[Rn]',  '5f-(6)5f(1)7s(2)',      243.d0, 5.9048d0), &
       atomic(  96,  'Cm',  '[Rn]',  '5f-(6)5f(1)6d-(1)7s(2)',244.d0, 5.8429d0), &
       atomic(  97,  'Bk',  '[Rn]',  '5f-(6)5f(3)7s(2)',      247.d0, 0.0000d0), &  ! ???  5.8160  where I took them?
       atomic(  98,  'Cf',  '[Rn]',  '5f-(6)5f(4)7s(2)',      251.d0, 0.0000d0), &  !      5.8440
       atomic(  99,  'Es',  '[Rn]',  '5f-(6)5f(5)7s(2)',      254.d0, 0.0000d0), &  !      5.8650
       atomic( 100,  'Fm',  '[Rn]',  '5f-(6)5f(6)7s(2)',      257.d0, 0.0000d0), &  !      5.8860
       atomic( 101,  'Md',  '[Rn]',  '5f-(6)5f(7)7s(2)',      258.d0, 0.0000d0), &  !      5.8930
       atomic( 102,  'No',  '[Rn]',  '5f-(6)5f(8)7s(2)',      259.d0, 0.0000d0), &  !      5.8860
       atomic( 103,  'Lr',  '[Rn]',  '5f-(6)5f(8)7s(2)7p-(1)',260.d0, 0.0000d0), &  !      5.9060
       atomic( 104,  'Rf',  '[Rn]',  '5f-(6)5f(8)6d(-2)7s(2)',265.d0, 0.0000d0) /)  !      5.8720

      Character(100) :: &
       He='1s', &
       Be='1s 2s', &
       Ne='1s 2s 2p- 2p', &
       Mg='1s 2s 2p- 2p 3s', &
       Ar='1s 2s 2p- 2p 3s 3p- 3p', &
       d3='1s 2s 2p- 2p 3s 3p- 3p 3d- 3d', &
       Zn='1s 2s 2p- 2p 3s 3p- 3p 3d- 3d 4s', &
       Kr='1s 2s 2p- 2p 3s 3p- 3p 3d- 3d 4s 4p- 4p', &
       d4='1s 2s 2p- 2p 3s 3p- 3p 3d- 3d 4s 4p- 4p 4d- 4d', &
       Cd='1s 2s 2p- 2p 3s 3p- 3p 3d- 3d 4s 4p- 4p 4d- 4d 5s', &
       Xe='1s 2s 2p- 2p 3s 3p- 3p 3d- 3d 4s 4p- 4p 4d- 4d 5s 5p- 5p', &
       f4='1s 2s 2p- 2p 3s 3p- 3p 3d- 3d 4s 4p- 4p 4d- 4d 4f- 4f 5s 5p- 5p', &
       d5='1s 2s 2p- 2p 3s 3p- 3p 3d- 3d 4s 4p- 4p 4d- 4d 4f- 4f 5s 5p- 5p 5d- 5d', &
       Hg='1s 2s 2p- 2p 3s 3p- 3p 3d- 3d 4s 4p- 4p 4d- 4d 4f- 4f 5s 5p- 5p 5d- 5d 6s', &
       Rn='1s 2s 2p- 2p 3s 3p- 3p 3d- 3d 4s 4p- 4p 4d- 4d 4f- 4f 5s 5p- 5p 5d- 5d 6s 6p- 6p', &
       f5='1s 2s 2p- 2p 3s 3p- 3p 3d- 3d 4s 4p- 4p 4d- 4d 4f- 4f 5s 5p- 5p 5d- 5d 5f- 5f 6s 6p- 6p'

      End Module atoms_par



!======================================================================
      Subroutine Def_atom(an,atom,atw,rms,core,conf)
!======================================================================
! ... define the ground state of atom
!----------------------------------------------------------------------
      Use atoms_par

      Integer :: an, i
      Real(8) :: atw, rms
      Character(*) :: atom,core,conf

! ... if atomic number = 0, first try to define atom by symbol:

      if(an.le.0.or.an.gt.n_atoms) then
       an = -1
       Do i = 1,n_atoms
        if(atom.ne.atoms(i)%symbol) Cycle
        an = i; Exit
       End do
       if(an.eq.-1) Return
      else
       atom = atoms(an)%symbol
      end if

! ... atomic weight:

      atw = atoms(an)%weight

! ... atomic radius:

      rms = atoms(an)%rrms
      if(rms.eq.0.d0) then
       if(an.le.90)  then
         rms=(0.836d0*atw**(1.d0/3.d0)+0.570d0)
       else
         rms=(0.77d0*atw**(1.d0/3.d0)+0.980d0)
       end if
      end if

! ... core label:

      core = atoms(an)%core

      if(core(1:1).eq.'[') then
       Select case(core(2:3))
         case ('He'); core = trim(He)
         case ('Be'); core = trim(Be)
         case ('Ne'); core = trim(Ne)
         case ('Mg'); core = trim(Mg)
         case ('Ar'); core = trim(Ar)
         case ('Zn'); core = trim(Zn)
         case ('Kr'); core = trim(Kr)
         case ('Cd'); core = trim(Cd)
         case ('Xe'); core = trim(Xe)
         case ('4f'); core = trim(f4)
         case ('Hg'); core = trim(Hg)
         case ('Rn'); core = trim(Rn)
         case ('3d'); core = trim(d3)
         case ('4d'); core = trim(d4)
         case ('5d'); core = trim(d5)
         CASE DEFAULT
         Stop 'Def_atom: problem with core'
       End Select
      end if

! ... configuration label:

      conf = atoms(an)%conf

      End Subroutine Def_atom


!======================================================================
      Real(8) Function DJ_fact(IL1,IL2,IS1,IS2,JOT1,JOT2,k)
!======================================================================
!     defines the J-dependence for reduced matrix elements of
!     k-pole electric and magnetic transition operator (type MA):
!
!     <LSJ||O[k]||L'S'J'> = DJ_fact *  <LS||O[k]||L'S'>
!
!     DJ_fact = (-1)^(L+S+J'+k) [J,J']^1/2  { L  S J },  S = S'
!                                           { J' k L'}
!
!     momenta are given in (2J+1)-representation
!======================================================================
      Implicit none
      Integer, intent(in) :: IL1,IL2,IS1,IS2,JOT1,JOT2,k
      Real(8), external :: Z_6j

      DJ_fact = 0.d0
      if(IS1.eq.IS2) DJ_fact = (-1)**((IS1+IL1+JOT2-3)/2+k)*  &
                               Dsqrt(1.d0*JOT1*JOT2)*         &
                               Z_6j(IL1,JOT1,IS1,JOT2,IL2,k+k+1)

      End Function DJ_fact


!====================================================================
      Real(8) Function DJM_fact(IL1,IL2,IS1,IS2,JOT1,JOT2,k)
!====================================================================
!     defines the J-dependence for reduced matrix elements of
!     k-pole magnetic operator MB:
!
!     <LSJ||MB[k]||L'S'J'> = DJM_fact *  <LS||MB[k]||L'S'>
!
!     DJM_fact =  [J,J',k]^1/2  { L  S  J }
!                               { L' S' J'}
!                               {k-1 1  k }
!
!     momenta are given in (2J+1)-representation
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: IL1,IL2,IS1,IS2,JOT1,JOT2,k
      Real(8), external :: Z_9j

      DJM_fact = Dsqrt(1.d0*JOT1*JOT2*(k+k+1))*         &
                 Z_9j(IL1,IS1,JOT1,IL2,IS2,JOT2,k+k-1,3,k+k+1)

      End Function DJM_fact


!====================================================================
      Integer Function ITRA (i1,i2,i3)
!====================================================================
!     check triangle relations for i1,i2,i3 momentums given
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: i1,i2,i3

      ITRA=1
      if(i1.gt.i2+i3.or.i2.gt.i1+i3.or.i3.gt.i1+i2) ITRA=0
      if(i1.lt.iabs(i2-i3).or.i2.lt.iabs(i1-i3).or. &
         i3.lt.iabs(i1-i2)) ITRA=0

      End function ITRA


!====================================================================
      Integer Function ITRI (i1,i2,i3)
!====================================================================
!     check triangle relations for i1,i2,i3 momentums given
!     in the (2J+1)-format
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: i1,i2,i3

      ITRI=1
      if(i1.gt.i2+i3-1.or.i2.gt.i1+i3-1.or.i3.gt.i1+i2-1) ITRI=0
      if(i1.lt.iabs(i2-i3)+1.or.i2.lt.iabs(i1-i3)+1.or. &
         i3.lt.iabs(i1-i2)+1) ITRI=0

      End Function ITRI




!=======================================================================
      Real(8) Function Cjkj (j1,k,j2)
!=======================================================================
!     Computes the relativistic reduced matrix element
!     (see Eq.15, Grant and Pyper, J.Phys.B9,761,1976)
!
!             k
!     (j1 || C  || j2)  =
!
!                                       (j1   k  j2 )
!         (-1)^(j1+1/2) sqrt([j1][j2])
!                                       (1/2  0 -1/2)
!
!     C(j,0,j) = sqrt([j])
!----------------------------------------------------------------------

      Implicit None

      Integer, Intent(in) :: j1,k,j2
      Real(8), External :: Z_3j2

      Cjkj = Z_3j2(j1,1,k+k,0,j2,-1) * &
            sqrt(real((j1+1)*(j2+1))) * (-1)**((j1+1)/2)

      END Function Cjkj


!=======================================================================
      Real(8) Function Ckap (kap1,k,kap2)
!=======================================================================
!     Computes the relativistic reduced matrix element
!                          k
!                (kap1 || C  || kap2)
!
!     see Cjkj as origin
!----------------------------------------------------------------------

      Implicit None

      Integer, Intent(in) :: kap1,k,kap2
      Integer, External :: j_kappa
      Real(8), External :: Cjkj
      Integer :: j1,j2

      j1 = j_kappa(kap1)
      j2 = j_kappa(kap1)
      Ckap = Cjkj(j1,k,j2)

      END Function Ckap



!=======================================================================
      Integer(4) Function kappa_lj(l,jj)
!=======================================================================
!     kappa = (l-j)*(2j+1);  jj = 2j
!-----------------------------------------------------------------------

      Integer(4) :: l,jj

      kappa_lj = (2*l-jj)*(jj+1)/2

      End Function kappa_lj


!=======================================================================
      Integer(4) Function l_kappa(kappa)
!=======================================================================

      Integer(4) :: kappa

      if(kappa.eq.0) Stop 'l_kappa: kappa=0'
      if(kappa.gt.0) then
       l_kappa =  kappa
      else
       l_kappa = -kappa-1
      end if

      End Function l_kappa


!=======================================================================
      Integer(4) Function j_kappa(kappa)
!=======================================================================
!     j_kappa = 2*j
!-----------------------------------------------------------------------

      Integer(4) :: kappa

      if(kappa.eq.0) Stop 'j_kappa: kappa=0'
      if(kappa.gt.0) then
       j_kappa =  kappa+kappa-1
      else
       j_kappa = -kappa-kappa-1
      end if

      End Function j_kappa




!====================================================================
      Real(8) Function T_LS_jk(l1,l2,s1,s2,L,S,js,jk,J)
!====================================================================
!     the recoupling coefficient from LS- to jK-coupling:
!
!     < (l1,l2)L,(s1,s2)S;J | (((l1,s1)js,l)K,s2)J >
!
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: l1,l2,s1,s2,L,S,js,jk,J
      Real(8), external :: Z_6j

      T_LS_jk = Z_6j(L,s1,jk,s2,J,S) * Z_6j(l2,l1,L,s1,jk,js) * &
                sqrt(1.d0*L*S*js*jk)* (-1)**((s2+J-l2-js)/2)

      End Function T_LS_jk

!====================================================================
      Real(8) Function T_LS_jj(l1,l2,s1,s2,L,S,j1,j2,J)
!====================================================================
!     recoupling coefficient from LS- to jj-coupling:
!
!     < (l1,l2)L,(s1,s2)S;J | ((l1,s1)j1,(l2,s2)j2;J >
!
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: l1,l2,s1,s2,L,S,j1,j2,J
      Integer :: k,k_min,k_max
      Real(8), external :: Z_6j

      T_LS_jj = 0.d0
      k_min = max(iabs(l2-j1),iabs(s2-J),iabs(L-s1))+1
      k_max = max(iabs(l2+j1),iabs(s2+J),iabs(L+s1))-1
      if(k_min.gt.k_max) Return
      Do k = k_min,k_max,2
       T_LS_jj = T_LS_jj + Z_6j(l2,s2,j2,J,j1,k) * &
                           Z_6j(L,S,J,s2,k,s1)   * &
                           Z_6j(l1,s1,j1,K,l2,L) * &
                           k * (-1)**(k-1)
      End do
      T_LS_jj = T_LS_jj * sqrt(1.d0*L*S*j1*j2)

      End Function T_LS_jj


!====================================================================
      Real(8) Function T_jj_jk(l2,s2,j1,j2,K,J)
!====================================================================
!     the recoupling coefficient from jK- to jj-coupling:
!
!     < (((l1,s1)j1,l2)K,s2)J | (((l1,s1)j1,(l2,s2)j2,J >
!
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: l2,s2,j1,j2,K,J
      Real(8), external ::  Z_6j2

      T_jj_jk = Z_6j2(j1,l2,K,s2,J,j2) * sqrt(1.d0*(K+1)*(j2+1)) &
                *(-1)**((j1+l2+s2+J)/2)

      End Function T_jj_jk


!====================================================================
      Real(8) Function ZCKLM (L1,M1,L2,M2,K)
!====================================================================
!     Condon-Shortly coefficients: ???
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: L1,M1,L2,M2,K
      Real(8), external :: CLEBSH
      Integer :: I,I1,I2, J,J1,J2, IM

      ZCKLM=0.0
      if(K.GT.L1+L2.OR.K.LT.IABS(L1-L2)) Return
      I=L1+L2+K
      if(I.NE.I/2*2) Return

      J = K+K+1
      J1= L1+L1+1
      J2= L2+L2+1
      I1=-M1-M1+1
      I2= M2+M2+1
      IM= I1+I2-1

      ZCKLM = SQRT(DBLE(J1*J2))/J*(-1)**M2  *  &
              CLEBSH(J1,1,J2,1,J,1)*CLEBSH(J1,I1,J2,I2,J,IM)

      End Function ZCKLM


!======================================================================
      Real(8) Function ZCLKL(K1,K2,K3)
!======================================================================
!     reduced matrix elements of spherical garmonics:
!
!     <k1||C[k2]||k3> = (-1)^k1 * sqtr[ (2*k1+1) * (2*k3+1) ] *
!                     3j(k1,0,k2,0,k3,0)
!
!---------------------------------------------------------------------
      Implicit none
      Integer, intent(In) :: k1,k2,k3
      Integer :: M, N
      Real(8), external :: ZCB

      M=K1+K2+K3
      N=M/2
      IF(N+N.NE.M) Return

      ZCLKL = ZCB(K1,K2,K3)
      ZCLKL = SQRT(ZCLKL*(K1+K1+1)*(K3+K3+1))
      IF(mod(N+K1,2).eq.1) ZCLKL=-ZCLKL

      End Function ZCLKL


!====================================================================
      Real(8) Function Z_3j (j1,m1,j2,m2,j3,m3)
!====================================================================
!     determines the value of the 3j-symbols without direct using of
!     factorials. The following expression for the 3j-symbols is used:
!         (A.P.JUCYS, A.A.BANDZAITIS, 1977)
!
!     3j{j1,m1,j2,m2,j3,m3} = delta(m1+m2,m3) * (2j3+1)^1/2 * {j1,j2,j3} *
!       sqrt[ (j1+m1)!*(j1-m1)!*(j2+m2)!*(j2-m2)!*(j3+m3)!*(j3-m3)! ]
!                         SUM(z) {   (-1)^z  /
!          [ z! *  (j1+j2-j3-z)! * (j1-m1-z)! * (j2-m2-z)! *
!                  (j3-j2+m1+z)! * (j3-j1-m2+z)! ] }
!
!     where {a,b,c}=sqrt[ (a+b-c)! * (a-b+c)! * (b-a+c)! / (a+b+c+1)! ]
!
!     If we introduce the auxiliary values a(i) and b(i)
!     (see below the text of program) then
!
!     3j =         (-1) ^ Sum[a(i)]
!          sqrt{ Pr[ (b(j)-a(i))! ] / [ Sum (b(j)-a(i))+1 ] }
!                  Sum(z) { (-1)^z  /
!          [  Pr[ (z-a(i))! ]  * Pr[ (b(j)-z)! ]   ] }
!
!     (below the moments are used in (2J+1)-representation)
!
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: j1,m1,j2,m2,j3,m3
      Integer :: i,i_max,k,kk,m,iz,iz_min,iz_max
      Real(8) :: x,y,z
      Integer :: a(3),b(3),J(16)

      Z_3j=0.0

      IF(M1+M2+M3-3.ne.0) Return    ! check of conservation rules
      J(1)= J1+J2-J3-1
      J(2)= J1-J2+J3-1
      J(3)= J2-J1+J3-1
      J(4)= J1+M1-2
      J(5)= J1-M1
      J(6)= J2-M2
      J(7)= J2+M2-2
      J(8)= J3+M3-2
      J(9)= J3-M3
      Do I=1,9
       IF(J(i).lt.0.or.mod(J(i),2).eq.1) Return
      End do

      a(1) = 0                         ! auxiliary values
      a(2) = (j2-j3-m1+1)/2
      a(3) = (j1-j3+m2-1)/2
      b(1) = (j1+j2-j3-1)/2
      b(2) = (j1-m1)/2
      b(3) = (j2+m2-2)/2

      IZ_min=MAX0(a(1),a(2),a(3))      ! limits of the sum
      IZ_max=MIN0(b(1),b(2),b(3))
      IF(IZ_max.LT.IZ_min) Return

      Do I=1,3                         ! constant factorial parameters
      Do K=1,3
       J(I+3*K-3)=b(i)-a(k)
      End do
      End do
      J(10)=(j1+j2+j3-3)/2+1

      Do I=1,3
       J(I+10)=IZ_min-a(i)               ! initial factorial parameters
       J(I+13)=b(i)-IZ_min               ! in the sum
      End do

      Z=0.0
      DO IZ=IZ_min,IZ_max                 ! summation

       I_max=0                            ! max. factorial
       Do I=1,16
        if(J(i).gt.I_max) I_max=J(i)
       End do

       Y=1.0
       DO I=2,I_max         ! estimation of one term in sum
        K=0                 ! K - the extent of the integer I in term
        DO M=1,9
         IF(J(M).GE.I) K=K+1
        End do
        IF(J(10).GE.I) K=K-1
        DO M=11,16
         IF(J(M).GE.I) K=K-2
        End do
        IF(K.EQ.0) Cycle

        X=DBLE(I)                   ! Y = Y * I ** K/2
        KK=IABS(K)/2
        IF(KK.GT.0) THEN
         DO M=1,KK
          IF(K.GT.0) Y=Y*X
          IF(K.LT.0) Y=Y/X
         END DO
        END IF
        IF(mod(K,2).EQ.+1) Y=Y*SQRT(X)
        IF(mod(K,2).EQ.-1) Y=Y/SQRT(X)
       End do

       IF(mod(IZ,2).eq.1) Y=-Y
       Z=Z+Y

       Do I=11,13                  ! new factorial parameters in sum
        J(I)=J(I)+1
       End do
       DO I=14,16
        J(I)=J(I)-1
       End do

      End do                       ! end of summation

      K=a(1)+a(2)+a(3)
      if(mod(k,2).ne.0) Z=-Z
      Z_3j=Z

      END Function Z_3j


!====================================================================
      Real(8) Function Z_3jj(j1,m1,j2,m2,j3,m3)
!====================================================================
!     3j-symbols for L-moments (without 2J+1 representation
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: j1,m1,j2,m2,j3,m3
      Real(8), external :: Z_3j

      Z_3jj=Z_3j(j1+j1+1,m1+m1+1,j2+j2+1,m2+m2+1,j3+j3+1,m3+m3+1)

      End Function Z_3jj


!====================================================================
      Real(8) Function Z_3j2(j1,m1,j2,m2,j3,m3)
!====================================================================
!     3j-symbols for moments in 2J representation
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: j1,m1,j2,m2,j3,m3
      Real(8), External :: Z_3j

      Z_3j2=Z_3j(j1+1,m1+1,j2+1,m2+1,j3+1,m3+1)

      End Function Z_3j2


!====================================================================
      Real(8) Function CLEBSH(J1,M1,J2,M2,J,M)
!====================================================================
!     determines the Clebsh-Gordon coefficients through the 3j-symbol
!     (the moments are used in (2J+1)-representation)
!
!     Call:  Z_3j
!--------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: J, M, J1, M1, J2, M2
      Real(8), external :: Z_3j

      Clebsh=(-1)**((j1-j2+m-1)/2)*sqrt(DBLE(J))*   &
             Z_3j(j1,m1,j2,m2,J,-m+2)

      END Function CLEBSH


!======================================================================
      Real(8) Function CLEBCH(L1,M1,L2,M2,L,M)
!======================================================================
!     determines the Clebsh-Gordon coefficients through the 3j-symbol
!     (the moments are in L-representation)
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: L, M, L1, M1, L2, M2
      Real(8), External :: Clebsh

      Clebch = Clebsh(l1+l1+1,m1+m1+1,l2+l2+1,m2+m2+1,l+l+1,m+m+1)

      End Function CLEBCH


!======================================================================
      Real(8) Function CLEBSH2(L1,M1,L2,M2,L,M)
!======================================================================
!     determines the Clebsh-Gordon coefficients through the 3j-symbol
!     (the moments are in 2J-representation)
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: L, M, L1, M1, L2, M2
      Real(8), External :: Clebsh

      Clebsh2 = Clebsh(l1+1,m1+1,l2+1,m2+1,l+1,m+1)

      End Function CLEBSH2


!======================================================================
      Real(8) Function Z_3j0 (K1,K2,K3)
!======================================================================
!     { k1 k2 k3}   -  3j-symbol with zero M-values
!     {  0  0  0}
!----------------------------------------------------------------------
      Implicit real(8) (A-H,O-Z)
      Integer, intent(in) :: k1,k2,k3

      Z_3j0=0.d0

      M=K1+K2+K3; N=M/2; IF(N+N.NE.M) RETURN; M=M+1

      N1=N-K1; N2=N-K2; N3=N-K3
      IF(N1.LT.0.OR.N2.LT.0.OR.N3.LT.0) RETURN
      M1=N1+N1; M2=N2+N2; M3=N3+N3

      A=1.d0
      DO I=2,M

      K=-1
      IF(M1.GE.I) K=K+1; IF(M2.GE.I) K=K+1; IF(M3.GE.I) K=K+1
      IF(N .GE.I) K=K+2
      IF(N1.GE.I) K=K-2; IF(N2.GE.I) K=K-2; IF(N3.GE.I) K=K-2

      IK=IABS(K); B=DBLE(I)
      IF(K.GT.0) THEN
       DO J=1,IK; A=A*B; END DO
      ELSEIF(K.LT.0) THEN
       DO J=1,IK; A=A/B; END DO
      END IF

      END DO

      Z_3j0=SQRT(A); IF(mod(N,2).eq.1) Z_3j0=-Z_3j0

      End Function Z_3j0

!======================================================================
      Real(8) Function ZCB(K1,K2,K3)
!======================================================================
!     CB =  3j(k1,0,k2,0,k3,0)**2
!
!     3j =  sqrt[ (2g-2k1)! (2g-2k2)! (2g-2k3)! / (2g+1)! ] *
!                   g! / [ (g-k1)! (g-k2)! (g-k3)! ],
!
!     where 2g = k1+k2+k3
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: k1,k2,k3
      Real(8) :: A,B,C
      Integer :: I,J, K,L, M,N, N1,N2,N3, M1,M2,M3

      ZCB=0.0;  M=K1+K2+K3;  N=M/2; if(N+N.NE.M) Return
      N1=N-K1; N2=N-K2; N3=N-K3
      IF(N1.LT.0.OR.N2.LT.0.OR.N3.LT.0) Return

      M1=N1+N1; M2=N2+N2; M3=N3+N3
      A = 1.d0/(M + 1)
      DO I = 1,N
       K = +1
       IF(M1.GE.I) K=K+1
       IF(M2.GE.I) K=K+1
       IF(M3.GE.I) K=K+1
       IF(N1.GE.I) K=K-2
       IF(N2.GE.I) K=K-2
       IF(N3.GE.I) K=K-2
       J = I + N
       L = -1
       IF(M1.GE.J) L=L+1
       IF(M2.GE.J) L=L+1
       IF(M3.GE.J) L=L+1
       B=I;  C=J;  A = A * B**K * C**L
      End do
      ZCB = A

      END FUNCTION ZCB


!======================================================================
      Real(8) Function CLB(K1,K2,K3)
!======================================================================
!     Clebsh(k1,0,k2,0;k3,0) = (-1)^(k1-k2) [k3]^1/2 3j(k1,0,k2,0;k3,0)
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: k1,k2,k3
      Real(8), external :: ZCB
      Real(8) :: A,B,C
      Integer :: I,J, K,L, M,N, N1,N2,N3, M1,M2,M3

      CLB = ZCB(k1,k2,k3)* (k3+k3+1)
      CLB = sqrt(CLB) * (-1)**(k1-k2)       ! sign?

      END FUNCTION CLB


!======================================================================
      Real(8) FUNCTION Z_6J (j1,j2,j3,j4,j5,j6)
!======================================================================
!     determination of 6j-symbols without direct using of factorials
!     accoding to formula:
!
!     6j{j1,j2,j3,j4,j5,j6) = {j1,j2,j3}*{j1,j5,j6}*{j4,j2,j3}*{j4,j5,j3}*
!                                SUM(z) {   (-1)^z * (z+1)!   /
!          [ (z-j1-j2-j3)! * (z-j1-j5-j6)! * (z-j4-j2-j3)! *(z-j4-j5-j3)! *
!              (j1+j2+j4+j5-z)! * (j1+j3+j4+j6-z)! * (j2+j3+j5+j6-z)! ]
!
!     where {a,b,c}=sqrt[ (a+b-c)! * (a-b+c)! * (b-a+c)! / (a+b+c+1)! ]
!
!     If we introduce the auxiliary values L(i)
!     (see below the text of program) then
!
!     6j = sqrt{ Pr(j=5,7,i=1,4) (L(j)-L(i))! / Pr(i=1,4) (L(i)+1)! }
!                Sum(z) { (-1)^z * (z+1)! /
!          [ Pr(i=1,4) (z-L(i))!  * Pr(j=5,7) (L(j)-z) ] }
!
!     (the momenta are used in (2J+1)-representation)
!
!----------------------------------------------------------------------

      Implicit none

      Integer(4), intent(in) :: j1,j2,j3,j4,j5,j6

      Integer(4) :: I, IZ_min, IZ_max, K, KK, M, IZ, I_max
      Integer(4) L(7),J(23)
      Real(8) :: X, R, C

      Z_6J = 0.0

      L(1)=j1+j2+j3-3                    ! auxiliary values
      L(2)=j1+j5+j6-3
      L(3)=j4+j2+j6-3
      L(4)=j4+j5+j3-3
      L(5)=j1+j2+j4+j5-4
      L(6)=j1+j3+j4+j6-4
      L(7)=j2+j3+j5+j6-4
      DO I=1,7
      IF(mod(L(I),2).eq.1) Return
       L(I)=L(I)/2
      END DO

      IZ_min=MAX0(L(1),L(2),L(3),L(4))   ! limits of the sum
      IZ_max=MIN0(L(5),L(6),L(7))
      IF(IZ_max.LT.IZ_min) Return

      Do I=1,4
       J(I)=L(I)+1
       Do K=5,7
        M=L(K)-L(I)
        IF(M.LT.0) Return                ! check of triangle rule
        J(4+3*I+K)=M
       End do
      End do
                                         ! initial factorial parameters
                                         ! in the sum
      Do I=5,8;  J(I)=IZ_min-L(I-4); End do
      Do I=9,11; J(I)=L(I-4)-IZ_min; End do

      C=0.0
      DO IZ=IZ_min,IZ_max           ! summation
       I_max=IZ+1

!      this limit for max. factorial follows from symmetry propeties:
!      let's a(i)=L(i) for i=1,4;  b(i)=L(j),j=5,7;
!      then  b(j)-a(i) <= a(k)  <= max[a(k)] = IZ_min;
!      also  a(j) <= b(j), then a(i)-a(j) <= b(i)-a(j) <= IZ_min;
!      and last (a(i)+1) < max(a(i))+1 < IZ_min+1;

       X=1.0
       DO I=2,I_max            ! estimation of one term in sum
        K=2                    ! K - the power of integer I in the term
        DO M=12,23; IF(J(M).GE.I) K=K+1;  End do
        DO M=1,4;   IF(J(M).GE.I) K=K-1;  End do
        DO M=5,11;  IF(J(M).GE.I) K=K-2;  End do
        IF(K.EQ.0) Cycle

        R=DBLE(I)               ! X = X * I ** K/2
        KK=IABS(K)/2
        IF(KK.GT.0) THEN
         DO M=1,KK
          IF(K.GT.0) X=X*R
          IF(K.LT.0) X=X/R
         END DO
        END IF
        IF(mod(K,2).EQ.+1) X=X*SQRT(R)
        IF(mod(K,2).EQ.-1) X=X/SQRT(R)
       End do

       IF(mod(IZ,2).eq.1) X=-X
       C=C+X
                                    ! new factorial parameters in sum
       Do I=5,8;  J(I)=J(I)+1; End do
       DO I=9,11; J(I)=J(I)-1; End do

      End do                        ! end of summation

      Z_6J=C

      END FUNCTION Z_6J


!====================================================================
      Real(8) FUNCTION Z_6jj(j1,j2,j3,j4,j5,j6)
!====================================================================

      IMPLICIT NONE

      Integer(4), intent(in) :: j1,j2,j3,j4,j5,j6
      Real(8), External :: Z_6j

      Z_6jj = Z_6j(j1+j1+1,j2+j2+1,j3+j3+1,j4+j4+1,j5+j5+1,j6+j6+1)

      End FUNCTION Z_6JJ

!====================================================================
      Real(8) FUNCTION Z_6j2(j1,j2,j3,j4,j5,j6)
!====================================================================

      IMPLICIT NONE

      Integer(4), intent(in) :: j1,j2,j3,j4,j5,j6
      Real(8), External :: Z_6j

      Z_6j2 = Z_6j(j1+1,j2+1,j3+1,j4+1,j5+1,j6+1)

      End FUNCTION Z_6j2


!======================================================================
      Real(8) FUNCTION Z_9j (j1,j2,j3,j4,j5,j6,j7,j8,j9)
!======================================================================
!
!     determination of 9j-symbols
!
! {j1 j2 j3}                           {j1 j2 j3} {j4 j5 j6} {j7 j8 j9}
! {j4 j5 j6} = SUM(j) (-1)^(2j) (2j+1)
! {j7 j8 j9}                           {j6 j9 j } {j2 j  j8} {j  j1 j4}
!
!
!     (the momenta are used in (2J+1)-representation)
!
!----------------------------------------------------------------------

      Implicit none

      Integer(4), intent(in) :: j1,j2,j3,j4,j5,j6,j7,j8,j9

      Integer(4) :: j,i1,i2
      Real(8), External :: Z_6j

      i1 = max(abs(j1-j9),abs(j2-j6),abs(j4-j8))+1
	  i2 = min(j1+j9,j2+j6,j4+j8)-1

	  Z_9j = 0.d0; if(i1.gt.i2) Return;  if(mod(i2-i1,2).ne.0) Return

      Do j = i1,i2,2
       Z_9j = Z_9j + j * (-1)**(j-1) * &
                     Z_6j(j1,j2,j3,j6,j9,j ) * &
                     Z_6j(j4,j5,j6,j2,j ,j8) * &
                     Z_6j(j7,j8,j9,j ,j1,j4)
      End do

      End FUNCTION Z_9j


!======================================================================
      Real(8) FUNCTION Z_9jj (j1,j2,j3,j4,j5,j6,j7,j8,j9)
!======================================================================

      Implicit none

      Integer(4), intent(in) :: j1,j2,j3,j4,j5,j6,j7,j8,j9
      Real(8), External :: Z_9j

      z_9jj = Z_9j (j1+j1+1,j2+j2+1,j3+j3+1,j4+j4+1,j5+j5+1,j6+j6+1, &
                    j7+j7+1,j8+j8+1,j9+j9+1)

      End FUNCTION Z_9jj


!======================================================================
      Subroutine Check_file(AF)
!======================================================================
!     stop with message if given file AF does not exist
!----------------------------------------------------------------------
      Character(*), Intent(in) :: AF
      Logical :: EX
      Inquire(file=AF,exist=EX)
      if(.not.EX) then
       write(*,*) ' can not find file  ',AF;  Stop
      end if
      End Subroutine Check_file


!======================================================================
      Integer Function Icheck_file(AF)
!======================================================================
!     check if the file AF exists
!----------------------------------------------------------------------
      Character(*), Intent(in) :: AF
      Logical :: EX
      Inquire(file=AF,exist=EX)
      Icheck_file = 1
      if(.not.EX) Icheck_file = 0
      End Function Icheck_file


!======================================================================
      Subroutine Find_free_unit(nu)
!======================================================================
!     provide free unit to open new file
!----------------------------------------------------------------------
      Implicit none
      Integer :: nu,i
      Logical :: connected

      nu = 0
      Do i=21,999
       INQUIRE(UNIT=i,OPENED=connected)
       if(connected) Cycle
       nu = i
       Exit
      End do
      if(nu.eq.0) Stop 'Find_free_unit: nu = 0'

      End Subroutine Find_free_unit


!======================================================================
      Subroutine Clean_a(a)
!======================================================================
!     remove blanks from the character string
!----------------------------------------------------------------------
      Character(*) :: a

      len = len_trim(a)
      m = 0
      Do i=1,len
       if(a(i:i).eq.' ') cycle
       m = m + 1
       a(m:m)=a(i:i)
      End do
      if(m.lt.len) a(m+1:len) = ' '

      End Subroutine Clean_a


!======================================================================
      Subroutine Conv_au (Z,AWT,au_cm,au_eV,ipri)
!======================================================================
!     "au_cm" is transformation from au to cm^-1.
!     It depends on the atomic weight AWT according to:
!        au_cm = 2*Ry_cm*AWT/(AWT+RATIO) cm-1/au,
!     where RATIO - the electron mass in au
!           Ry_cm - Rydberg constatnt in cm-1
!           Ry_eV - Rydberg constatnt in eV
!     These constants are taken from  http://physics.nist.gov (2001)
!----------------------------------------------------------------------
      Implicit none
      Real(8), intent(in) :: Z,AWT
      Real(8), intent(out) :: au_cm, au_ev
      Integer, intent(in) :: ipri
      Real(8) :: A
      Real(8), parameter :: RATIO = 5.485799110d-4
      Real(8), parameter :: Ry_cm = 109737.31568549d0
      Real(8), parameter :: Ry_ev = 13.6056917d0

! ... determine suitable atomic weight if any

      if(AWT.ne.0.d0) then
        A = AWT
      else
        IF ( Z .EQ. 1.) THEN
           A = 1.
         ELSE IF ( Z .GT. 10.) THEN
           A = 2*Z+1 + (Z-11)/2
         ELSE IF ( MOD(INT(Z),2) .EQ. 0 .OR. Z .EQ. 7. ) THEN
           A = 2*Z
         ELSE
           A = 2*Z+1
        END IF
      end if

      au_cm = 2*Ry_cm * A / (A + RATIO)
      au_eV = 2*Ry_ev * A / (A + RATIO)

      if(ipri.gt.0) then
       write(ipri,'(/a,f10.2)') 'Atomic weight = ',A
       write(ipri,'(a,f20.10)') &
            'a.u. --> cm-1 , au_cm = ', au_cm
       write(ipri,'(a,f20.10/)') &
            'a.u. --> eV   , au_eV = ', au_eV
      end if

      End Subroutine Conv_au


!======================================================================
      Subroutine Full_mat_sym(ns,ks,x,y,sym)
!======================================================================
!     Builds full matrix from different storage modes for banded matrix
!     (see p.328 in BSR description, CPC 174 (2006))
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: ns,ks
      Real(8), intent(in) :: x(ns,*)
      Real(8), intent(out) :: y(ns,ns)
      Character(*), intent(in) :: sym
      Integer ::  i,j, jp, imin,imax

      y = 0.d0

      Select case(sym)

      Case('s')   ! symmetric upper-column storage

       do jp=1,ks;  do i=1,ns-jp+1;  j=i+jp-1
        y(i,j) = x(i,jp)
        y(j,i) = x(i,jp)
       end do; end do

      Case('l')   ! symmetric lower-column storage

       do jp=1,ks;  do i=ks+1-jp,ns;  j=i+jp-ks
        y(i,j) = x(i,jp)
        y(j,i) = x(i,jp)
       end do; end do

      Case('n')   !  non-symmetric banded

       do jp = 1,ks+ks-1
         imin=max0( 1, 1 + ks-jp)
         imax=min0(ns,ns + ks-jp)
         do i = imin,imax;  j=i+jp-ks
           y(i,j) = x(i,jp)
         end do
       end do

      Case('x')   ! just full (for completeness)

        y(1:ns,1:ns) = x(1:ns,1:ns)

      End Select

      End Subroutine Full_mat_sym


!======================================================================
      Real(8) FUNCTION ALEGF (L,M,X,norm)
!======================================================================
! ... ASSOCIATED LEGENDRE POLYNOMIALS:
! ... NORM=0: UNNORMALISED LEGENDRE FUNCTIONS P_LM(x) (I.E. ORIGINAL)
! ....NORM=1: NORMALISED: * SQRT[ (2L+1)/2 * (L-M)!/(L+M)! ]
!----------------------------------------------------------------------
      Implicit real(8) (A-H,O-Z)
      Real(8) :: zero = 0.d0; one = 1.d0; two = 2.d0

      ALEGF=zero; MA=IABS(M); if(L.lt.MA) Return

      ALEGF=one
      FACT=one
      if(M.ne.0) then
       Do I=L-MA+1,L+MA; FACT=FACT*I; End do; FACT=one/FACT
      end if
      if(M.lt.0) ALEGF=FACT
      if(norm.eq.1.and.m.ge.0) ALEGF = ALEGF * sqrt(FACT*(L+L+1)/2)
      if(norm.eq.1.and.m.lt.0) ALEGF = ALEGF * sqrt((L+L+1)/(2*FACT))
      if(L.eq.0) Return

      if(abs(x).eq.one) then
       if(m.ne.0) ALEGF=zero; if(x.lt.zero) ALEGF=ALEGF*(-1)**L
       Return
      end if

      FACT=one; Do I=1,L+L-1,2; FACT=FACT*I; End do

      if(MA.eq.L) then
       ALEGF=ALEGF*FACT*(one-x*x)**(0.5d0*MA)
       Return
      end if
      if(MA.eq.L-1) then
       ALEGF=ALEGF*FACT*(one-x*x)**(0.5d0*MA)*x
       Return
      end if

! ... for m <= 1, use the recurrence relations in respect to l:
! ... P(l+1) = [(2l+1)P(l)-(l+m)p(l-1)] / (l-m+1)

      if(MA.le.1) then
       if(MA.eq.0) then
        P0=one;  P1=x
       else
        P0=zero; P1=DSQRT(one-x*x)
       end if
       DO N=1,L-1
        PN=((N+N+1)*x*P1-(N+MA)*P0)/(N-MA+1); P0=P1; P1=PN
       End do
        ALEGF=PN*ALEGF; Return
      end if

! ... use the recurrence relations in respect to m:
! ... P(m+2) = 2(m+1)x/sqrt(1-x^2)P(m+1)-(l-m)(l+m+1)P(m)

      z = x/sqrt(one-x*x)
      P2 = FACT*(one-x*x)**(0.5d0*L)
      P1 = P2 * z; z = z * two
      Do MM = L-2,MA,-1
       PM=((MM+1)*z*P1-P2)/((L-MM)*(L+MM+1)); P2=P1; P1=PM
      End do

      ALEGF=PM*ALEGF; Return

      END FUNCTION ALEGF

!=====================================================================
      Real(8) FUNCTION ALEGFM (L,M,X,norm)
!=====================================================================
!     sign correction for origional function ALEGF
!---------------------------------------------------------------------
      Implicit real(8) (A-H,O-Z)
      ALEGFM = ALEGF (L,M,X,norm)
      if(m.gt.0) ALEGFM = ALEGFM * (-1)**M
      End FUNCTION ALEGFM


!=======================================================================
      SUBROUTINE CGAMMA (ARGR,ARGI,RESR,RESI)
!=======================================================================
!   This subroutine returns in (RESr,RESi) the complex Gamma function  !
!   of the complex argument (ARGr,ARGi).                               !
!                                                                      !
!   Only RESR is nonzero if ARGI is zero.                              !
!                                                                      !
!   The  ARCTAN function required must return angles (in radians) in   !
!   the range  [0,2*pi).                                               !
!                                                                      !
!   Written by Farid A Parpia, at Oxford    Last update: 06 Oct 1992   !
!=======================================================================

      Implicit real(8) (A-H,O-Z)

!----------------------------------------------------------------------
!   These are the Bernoulli numbers B02, B04, ..., B14, expressed as
!   rational numbers. From Abramowitz and Stegun, p. 810.

      Real(8), SAVE :: PI, HLNTPI, eps, MAXEXP, MINEXP
      Real(8), DIMENSION(7), SAVE :: FN,FD
      DATA FN/ 1.0d0, -1.0d0,  1.0d0, -1.0d0,  5.0d0, -691.0d0, 7.0d0/
      DATA FD/ 6.0d0, 30.0d0, 42.0d0, 30.0d0, 66.0d0, 2730.0d0, 6.0d0/

!----------------------------------------------------------------------
      LOGICAL, SAVE :: FIRST, NEGARG
      DATA FIRST/.TRUE./
      Integer :: i

!   On the first entry to this routine, set up the constants required
!   for the reflection formula (cf. Abramowitz and Stegun 6.1.17) and
!   Stirling's approximation (cf. Abramowitz and Stegun 6.1.40).

      IF(FIRST) THEN

         PI = ACOS(-1.d0)
         HLNTPI = 0.5D0*LOG (PI+PI)
         eps = 2*EPSILON(1.d0)
         MAXEXP = MAXEXPONENT(1.d0) * LOG(2.d0)
         MINEXP = MINEXPONENT(1.d0) * LOG(2.d0)

         Do i = 1,7;  DI = i+i
           FN(I) = FN(I)/(FD(I)*DI*(DI-1.0d0))
         End do

         FIRST = .FALSE.

      END IF

!----------------------------------------------------------------------
! ... Cases where the argument is real

      IF (ARGI .EQ. 0.d0) THEN

! ... Cases where the argument is real and negative

         IF (ARGR .LE. 0.d0) THEN

! ... Stop with an error message if the argument is too near a pole

            IF (ABS (DBLE (NINT (ARGR))-ARGR) .LE. eps) &
               STOP 'CGAMMA: Argument too close to a pole.'

!   Otherwise use the reflection formula (Abramowitz and Stegun 6.1.17)
!   to ensure that the argument is suitable for Stirling's formula

            ARGUM = PI/(-ARGR*SIN(PI*ARGR))
            CLNGI = 0.d0
            IF (ARGUM .LT. 0.d0) THEN
              ARGUM = -ARGUM
              CLNGI = PI
            ENDIF
            FACNEG = LOG (ARGUM)
            ARGUR = -ARGR
            NEGARG = .TRUE.

! ...  Cases where the argument is real and positive

         ELSE

            CLNGI = 0.0d0
            ARGUR = ARGR
            NEGARG = .FALSE.

         END IF

!   Use Abramowitz and Stegun formula 6.1.15 to ensure that
!   the argument in Stirling's formula is greater than 10

         OVLFAC = 1.0d0
         Do
          IF (ARGUR .GE. 10.0d0) Exit
          OVLFAC = OVLFAC*ARGUR
          ARGUR = ARGUR+1.0d0
         End do

!   Now use Stirling's formula to compute Log (Gamma (ARGUM))

         CLNGR = (ARGUR-0.5d0)*LOG(ARGUR)-ARGUR+HLNTPI
         FAC = ARGUR
         OBASQ = 1.0d0/(ARGUR*ARGUR)
         Do I = 1,7
            FAC = FAC*OBASQ; CLNGR = CLNGR+FN(I)*FAC
         End do

!   Include the contributions from the recurrence and reflection
!   formulae

         CLNGR = CLNGR-LOG (OVLFAC)
         IF (NEGARG) CLNGR = FACNEG-CLNGR

      ELSE

!   Cases where the argument is complex

         ARGUR = ARGR
         ARGUI = ARGI
         ARGUI2 = ARGUI*ARGUI

!   Use the recurrence formula (Abramowitz and Stegun 6.1.15)
!   to ensure that the magnitude of the argument in Stirling's
!   formula is greater than 10

        OVLFR = 1.0d0
        OVLFI = 0.0d0
        Do
         ARGUM = SQRT (ARGUR*ARGUR+ARGUI2)
         IF (ARGUM .GE. 10.d0) Exit
         TERMR = OVLFR*ARGUR-OVLFI*ARGUI
         TERMI = OVLFR*ARGUI+OVLFI*ARGUR
         OVLFR = TERMR
         OVLFI = TERMI
         ARGUR = ARGUR+1.0d0
        End do

!   Now use Stirling's formula to compute Log (Gamma (ARGUM))

         ARGUR2 = ARGUR*ARGUR
         TERMR = 0.5d0*LOG (ARGUR2+ARGUI2)
         TERMI = ATAN2 (ARGUI,ARGUR)
         if(TERMI.lt.0.d0) TERMI = TERMI + PI + PI
         CLNGR = (ARGUR-0.5d0)*TERMR - ARGUI*TERMI-ARGUR+HLNTPI
         CLNGI = (ARGUR-0.5d0)*TERMI + ARGUI*TERMR-ARGUI
         FAC = (ARGUR2+ARGUI2)**(-2)
         OBASQR = (ARGUR2-ARGUI2)*FAC
         OBASQI = -2.0d0*ARGUR*ARGUI*FAC
         ZFACR = ARGUR
         ZFACI = ARGUI
         Do I = 1,7
            TERMR = ZFACR*OBASQR-ZFACI*OBASQI
            TERMI = ZFACR*OBASQI+ZFACI*OBASQR
            FAC = FN(I)
            CLNGR = CLNGR+TERMR*FAC
            CLNGI = CLNGI+TERMI*FAC
            ZFACR = TERMR
            ZFACI = TERMI
         End do

!   Add in the relevant pieces from the recurrence formula

         CLNGR = CLNGR - 0.5d0 * LOG (OVLFR*OVLFR+OVLFI*OVLFI)
         A = ATAN2(OVLFI,OVLFR); if(A.lt.0.d0) A = A + PI + PI
         CLNGI = CLNGI - ATAN2(OVLFI,OVLFR)

      END IF

!   Now exponentiate the complex Log Gamma function to get
!   the complex Gamma function

      IF (CLNGR.ge.MAXEXP.or.CLNGR.le.MINEXP) then
          write(*,*) CLNGR,MAXEXP,MINEXP
         STOP 'CGAMMA: Argument to exponential function out of range.'
      End if

      FAC = EXP(CLNGR)
      RESR = FAC*COS(CLNGI)
      RESI = FAC*SIN(CLNGI)

      END SUBROUTINE CGAMMA



!======================================================================
      FUNCTION RGAMMA(X)
!======================================================================
!     COMPUTES GAMMA(X) TO ABOUT 14 FIGURES
!     presented by Peter J. Mohr
!----------------------------------------------------------------------
      Use zconst

      Implicit none
      Real(dp), intent(in) :: X
      Real(dp) :: RGAMMA, Y,Z, A1,A2,A3,A4,A5,A6,A7, C,D,G

      Y=X
      G=one
    1 IF(Y.GT.seven)  GO TO 2
      G=G*Y
      Y=Y+one
      GO TO 1
    2 Z=-one/(Y*Y)

      C =.398942280401433e+0_dp
      D =.295506535947712e-1_dp

      A1=.282001658833287e+1_dp
      A2=.940005529444291e-1_dp
      A3=.268573008412655e-1_dp
      A4=.201429756309491e-1_dp
      A5=.284850160437664e-1_dp
      A6=.648894925920085e-1_dp
      A7=.216924352948683e+0_dp

      RGAMMA = (Y-half)*LOG(Y) - Y - LOG(C*G) + &
               (A1+Z*(A2+Z*(A3+Z*(A4+Z*(A5+Z*(A6+Z*(A7+Z)))))))*D/Y

      RGAMMA = EXP(RGAMMA)

      END FUNCTION RGAMMA


!=======================================================================
      SUBROUTINE DCWF (n,kappa,Z,E,NTP,R,P,Q)
!=======================================================================
!   This subroutine computes the  Dirac-Coulomb  bound-state orbital
!   radial wavefunction.
!
!   Input:
!
!      n          The (usual) principal quantum number
!      kappa      The relativistic angular quantum number
!      Z          The effective nuclear charge
!      R(1:NTP)   Radial grid
!
!   Output:
!
!      E          The Dirac-Coulomb Eigenenergy (E-mc^2)
!
!      P          r times the large component wavefunction of
!
!      Q          r times the small component wavefunction of
!
!   Call(s) to: CGAMMA
!
!   Written by Farid A Parpia, at Oxford    Last Update: 14 Oct 1992
!
!=======================================================================
!   For a given value of n > 0, there are 2*n-1 eigenfunctions:
!   n  with  kappa = -1,-2,...,-n
!   n-1 with kappa = 1,2,...,n-1
!
!   k = iabs(kappa);  gamma = sqrt[k^2-(alfa*Z)^2]
!
!   N = sqrt[n^2 - 2(n-k)(k-gamma)]
!
!   E(n,k) = c^2 / sqrt[1 + (alfa*Z)^2 / (gamma+n-k)^2 ]
!
!   x = 2*Z/N*r
!
!   P_nk(r) = sqrt[1+E(n,k)/c^2]  N(n,k)  e^-x/2  x^gamma *
!             [(N-kappa) F2  -  (n-k) F1]
!
!   Q_nk(r) = sqrt[1-E(n,k)/c^2]  N(n,k)  e^-x/2  x^gamma *
!             [(N-kappa) F2  +  (n-k) F1]
!
!   N(n,k) = 1/[N*G(2*gamma+1)]  *
!            sqrt [ [Z*G(2*gamma+1+n-k)] / [2(n-k)!(N-kappa)] ]
!
!   G(x) - GAMMA function
!
!   F2 = F(-n+k,2*gamma+1,x);  F1 = F(-n+k+1,2*gamma+1,x)
!
!   F(a,b,x) = 1 + a/b x^1/1!  + [a(a+1)]/[b(b+1)] x^2/2! + ...
!=======================================================================
      Use zconst, only: alpha, C => c_au

      Implicit real(8) (A-H,O-Z)
      Real(8) :: R(*), P(*), Q(*)
      Real(8) :: T1(n),T2(n)

! ... Check the input arguments:

      IF(n.le.0)     Stop 'DCWF: Principal quantum number < 0'
      IF(kappa.eq.0) Stop 'DCWF: Kappa quantum number = 0'
      IF(kappa.eq.n) Stop 'DCWF: Kappa quantum number = n'
      IF(kappa.gt.n) Stop 'DCWF: Kappa quantum number > n'
      IF(Z.le.0.d0)  Stop 'DCWF: Nuclear charge is too small, <= 0'
      IF(Z.gt.C)     Stop 'DCWF: Nuclear charge exceeds limit, c'

! ... Now determine all the parameters:

      fn = DBLE (n)
      fkappa = DBLE (kappa)
      k = IABS (kappa);  fk = DBLE (k)
      nr = n-k;  fnr = DBLE (nr)
      Za = Z*alpha
      gamma = SQRT (fk*fk-Za*Za)
      gg = gamma + gamma + 1.d0
      BIGN = SQRT (fn*fn-2.d0*fnr*(fk-gamma))
      EPS = 1.d0 /SQRT(1.d0+(Za/(gamma+fnr))**2)

! ... EPS is the total energy divided by C*C

      E = (1.d0-EPS)*C*C         !  E => E-mc^2

! ... normalization constant N(n,k):

      NRFAC=1; Do I=1,NR; NRFAC = NRFAC*I; End do

      ARGI = 0.d0
      ARGR = gg+FNR;    CALL CGAMMA (ARGR,ARGI,RGAMM1,DUMMY)
      ARGR = gg;        CALL CGAMMA (ARGR,ARGI,RGAMM2,DUMMY)

      FAC = - SQRT (RGAMM1)/(RGAMM2*SQRT (DBLE(NRFAC))) &
            * SQRT (Z/(2.d0*BIGN*BIGN*(BIGN-FKAPPA)))

! ... Ensure that the slope of the large-component function is
! ... positive at the origin:

      IF (KAPPA .GT. 0) FAC = -FAC

      FG = FAC * SQRT(1.d0+EPS)
      FF = FAC * SQRT(1.d0-EPS)

! ...  Now set up the coefficients of the confluent hypergeometric
! ...  functions  F(-NR+1,2*GAMMA+1;RHO)  and  F(-NR,2*GAMMA+1;RHO)
! ...  in the workspace arrays  TA  and  TB, respectively

      nt1=0; if(nr.gt.0) nt1=nr-1;  nt2=nr

      FAC = 1.d0;  FACN = FAC
      A   = -fnr;  AN1  = A + 1.d0;  AN2 = A
      B   =   gg;  BN   = B

      K = 0
    2 K = K+1
      FDEN = 1.d0/(FACN*BN)
      if (K .LE. nt1)      T1(K) = AN1*FDEN
      if (K .LE. nt2) then
                           T2(K) = AN2*FDEN
         A = A + 1.d0;     AN1  = AN1*(A+1.d0); AN2 = AN2*A
         B = B + 1.d0;     BN   = BN*B
         FAC = FAC + 1.d0; FACN = FACN*FAC
         go to 2
      end if

! ...  Now tabulate the function over the entire grid

      FAC = (Z+Z)/BIGN
      a1 = FNR; a2 = BIGN-FKAPPA
      Do i = 1,NTP
       x = FAC*R(i); y = x
       F1 = 1.d0; F2 = 1.d0
       K = 0
    3  K = K+1
       if (K .LE. nt1) F1 = F1+T1(K)*y
       if (K .LE. nt2) then;  F2 = F2+T2(K)*y; y=y*x; go to 3;  end if

       F1 = a1*F1;  F2 = a2*F2
       OVLFAC = EXP(-0.5d0*x) * x**gamma

       P(I) = FG*OVLFAC*(F1-F2)
       Q(I) = FF*OVLFAC*(F1+F2)

      End do

      END SUBROUTINE DCWF


!======================================================================
      Real(8) Function E_dcwf(n,k,Z)
!======================================================================
!     energy of the Dirac-Coulomb orbital nk:
!
!     E(n,k) = c^2 / sqrt[1 + (alfa*Z)^2 / (gamma+n-k)^2 ]  - c^2
!
!     gamma = sqrt[k^2-(alfa*Z)^2]
!----------------------------------------------------------------------
      Use zconst, only: alpha, C => c_au

      Implicit real(8) (A-H,O-Z)

      Za = Z*alpha

      gamma = SQRT (k*k-Za*Za)

      E = 1.d0 /SQRT(1.d0+(Za/(gamma+n-k))**2)

      E_dcwf = -(1.d0-E)*C*C

      End Function E_dcwf


!========================================================================
      SUBROUTINE DCME(n,kappa,z,ar1,ar2,am1,am2,am3)
!========================================================================
! ... provides radial moments <r^k> for Dirac-Coulomb wave functions
! ... expressions are taken from DRAKE HANDBOOK,2006
!------------------------------------------------------------------------
!     let  x = 2Z * r, then
!
!     <x^2> = 2*N^2*[(5*N^2-2*kappa^2)*R^2 + (1-gamma^2) - 3*kappa*R]
!
!     <x>   = -kappa + (3*N^2-kappa^2)*R
!
!     <x-1> = [n*gamma + (k-gamma)*k] / [2*gamma*N^3]
!
!     <x-2> = [kappa^2 * R] / [2*gamma^2*N^3*(2*gamma-sgn(kappa))]
!
!     <x-3> = [N^2 + 2*gamma^2*kappa^2 - 3* N^2*kappa*R] /
!             [4*N^5*gamma*(gamma^2-1)*(4*gamma^2-1)]
!
!     where  R = sqrt[1 - Z^2/(N^2*c^2)]
!            N = sqrt[n^2-2*(n-k)(k-gamma)]
!            gamma = sqrt[kappa^2-(Z/c)^2]
!            k = abs(kappa)
!------------------------------------------------------------------------
      USE zconst, ONLY: c_au

      Implicit none
      Integer, intent(in) :: n, kappa
      Real(8), intent(in) :: z
      Real(8), intent(out) :: ar1,ar2,am1,am2,am3
      Integer :: k, kk
      Real(8) :: Za, gamma, BigN, R, RR, bb, gg, x,xx,xxx

      k = iabs(kappa); kk = k*k
      Za = (z/c_au)**2
      gamma = sqrt(kk-Za); gg = gamma*gamma;
      BigN = sqrt(n*n-2*(n-k)*(k-gamma)); bb=BigN*BigN
      R = sqrt(1-Za/bb); RR = R*R

      ar2 = 4*bb*( (5*bb-2*kk)*RR + (1-gg) - 3*kappa*R )
      ar1 = -kappa + (3*bb-kk)*R
      am1 = (n*gamma + (k-gamma)*k) / (2*gamma*bb*BigN)
      am2 = kk*R / (2*gg*bb*BigN*(2*gamma-sign(1,kappa)))
      am3 = (bb + 2*gg*kk - 3*bb*kappa*R) / &
            (4*bb*bb*BigN*gamma*(gg-1)*(4*gg-1))

      x = 2*Z !/BigN;
      xx=x*x; xxx=xx*x

      ar2 = ar2 / xx
      ar1 = ar1 / x
      am1 = am1 * x
      am2 = am2 * xx
      am3 = am3 * xxx

      End SUBROUTINE DCME


!======================================================================
      Subroutine  Exp_dcwf(n,k,Z,ap,bp,cp)
!======================================================================
!     Expectation values for Dirac Coulomb w.f.
!     after notes of Mohr (???)
!----------------------------------------------------------------------
      Use zconst

      Implicit none
      Integer, intent(in) :: n,k
      Real(8), intent(in) :: Z
      Real(8) :: ap(-4:4),bp(-4:4),cp(-4:4)
      Real(8) :: la, nr, dn,dk, g, a, E, x,y, p
      Integer :: ip

      dn = DBLE(n); dk = DBLE(k)
      g = Z/c_au; la=sqrt(dk*dk-g*g); nr=dn-abs(dk)
      a = g/sqrt( (nr+la)**2+g*g )
      E = sqrt(one-a*a)

      ap(0)=one
      bp(0)=E
      cp(0)=dk * a**2 /g

      ap(-1) = a**3 / g**2 * (dk**2 /la + nr)
      bp(-1) = a*a / g
      cp(-1) = dk * a**3 / (g*la)

      ap(-2) = two*dk*a**3 / (g*la) * (two*dk*E-one) / (four*la*la-one)
      bp(-2) = two*a**3 / (g*la) * (two*la**2-dk*E) / (four*la*la-one)
      cp(-2) = two*a**3 / la * (two*dk*E-one) / (four*la*la-one)

      x = two*a**3 / la / (four*la*la-one)
      y = la*la-one

      ap(-3) = x * ( three*dk*E*(dk*E-one)/y - one)
      bp(-3) = x * E *( three*(la*la-dk*E)/y - one)
      cp(-3) = x /g * ( three*g*g*E*(dk*E-one)/y - a*a*dk)

      x = (four*la**2-three**2)
      y = E*bp(-3) - ap(-3)

      ap(-4) = two/three/x * &
               (four*g*dk*y + (-six*dk*E+four*g*g+nine)*cp(-3))
      bp(-4) = two/x * (two*g*y + (-three*E+two*dk)*cp(-3))
      ap(-4) = two/three/x * &
               ((four*dk**2-nine)*y + two*g*(two*dk-three*E)*cp(-3))

      ip=-3; p = DBLE(ip)
             x = (four*la*la-p*p)
             y = (E*bp(ip)-ap(ip))
      ap(ip-1) = two / p /x  * ( four*g*dk*y+ &
                (two*dk*p*E+four*g*g+p*p)*cp(ip) )
      bp(ip-1) = two / x * ( two*g*y + (two*dk+p*E)*cp(ip) )
      cp(ip-1) = two / p /x  * ( (four*dk*dk-p*p)*y + &
                 two*g*(two*dk+p*E)*cp(ip) )

      End Subroutine  Exp_dcwf




!----------------------------------------------------------------------
      SUBROUTINE gaussj(a,n,np,b,m,mp)
!----------------------------------------------------------------------
!     Linear equation solution by Gauss-Jordan elimination.
!     On input:
!     a(1:n,1:n) is an input matrix stored in an array of physical
!                dimensions np by np.
!     b(1:n,1:m) is an input matrix containing the m right-hand side
!                vectors, stored in an array np by mp.
!     On output:
!     a(1:n,1:n) is replaced by its matrix inverse
!     b(1:n,1:m) is replaced by the set of solution vectors.
!     (from numerical recipe,  www.nr.com)
!----------------------------------------------------------------------
      Integer :: m,mp,n,np
      Real(8) :: a(np,np),b(np,mp)
      Real(8) :: big,dum,pivinv, zero=0.d0, one=1.d0
      Integer :: i,icol,irow,j,k,l,ll
      Integer :: indxc(n),indxr(n),ipiv(n)

! ... The integer arrays ipiv, indxr, and indxc are used
! ... for bookkeeping on the pivoting.

      ipiv = 0

! ... This is the main loop over the columns to be reduced.

      Do i=1,n;   big=zero

! ... This is the outer loop of the search for a pivot element.

      do j=1,n;  if(ipiv(j).eq.1) Cycle
      do k=1,n;  if(ipiv(k).ne.0) Cycle
         if (abs(a(j,k)).ge.big) then
          big=abs(a(j,k)); irow=j; icol=k
         end if
      end do
      end do
      ipiv(icol)=ipiv(icol)+1

! ... We now have the pivot element, so we interchange rows, if needed,
! ... to put the pivot element on the diagonal. The columns are not
! ... physically interchanged, only relabeled:
! ... indxc(i), the column of the ith pivot element, is the ith column
! ...           that is reduced, while
! ... indxr(i)  is the row in which that pivot element was originally
! ...           located. If indxr(i) /= indxc(i) there is an implied
! ... column interchange. With this form of bookkeeping, the solution
! ... b’s will end up in the correct order, and the inverse matrix will
! ... be scrambled by columns.

      if (irow.ne.icol) then
       do l=1,n
        dum=a(irow,l); a(irow,l)=a(icol,l); a(icol,l)=dum
       end do
       do l=1,m
        dum=b(irow,l); b(irow,l)=b(icol,l); b(icol,l)=dum
       end do
      end if

! ... We are now ready to divide the pivot row by the pivot element,
! ... located at irow and icol.

      indxr(i)=irow
      indxc(i)=icol
      if (a(icol,icol).eq.zero) Stop 'singular matrix in gaussj'
      pivinv=one/a(icol,icol)
      a(icol,icol)=one
      do l=1,n; a(icol,l)=a(icol,l)*pivinv; end do
      do l=1,m; b(icol,l)=b(icol,l)*pivinv; end do

! ... Next, we reduce the rows, except for the pivot one, of course.

      do ll=1,n; if(ll.eq.icol) Cycle
       dum=a(ll,icol)
       a(ll,icol)=zero
       do l=1,n; a(ll,l)=a(ll,l)-a(icol,l)*dum; end do
       do l=1,m; b(ll,l)=b(ll,l)-b(icol,l)*dum; end do
      end do

     End do

! ... This is the end of the main loop over columns of the reduction.
! ... It only remains to unscramble the solution in view of the column
! ... interchanges. We do this by interchanging pairs of columns in
! ... the reverse order that the permutation was built up.

     do l=n,1,-1; if(indxr(l).eq.indxc(l)) Cycle
     do k=1,n
      dum=a(k,indxr(l)); a(k,indxr(l))=a(k,indxc(l)); a(k,indxc(l))=dum
     end do
     end do

     End Subroutine gaussj


!======================================================================
      SUBROUTINE gauleg(x1,x2,x,w,n)
!======================================================================
!     Given the lower and upper limits of integration x1 and x2,
!     and given n, this routine returns arrays x(1:n) and w(1:n)
!     of length n, containing the abscissas and weights of the
!     Gauss-Legendre n-point quadrature formula.
!     (from numerical recipe,  www.nr.com)
!-----------------------------------------------------------------------
      Integer n
      Double precision :: x1,x2,x(n),w(n)
      Double precision, parameter :: EPS=1.d-14 ! the relative precision.
      Integer :: i,j,m
      Double precision :: p1,p2,p3,pp,xl,xm,z,z1

      !  High precision is a good idea for this routine.

      m=(n+1)/2              !  The roots are symmetric in the interval,
      xm=0.5d0*(x2+x1)       !  so we only have to find half of them.
      xl=0.5d0*(x2-x1)

      Do i=1,m               !  Loop over the desired roots.

         z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))

      !  Starting with the above approximation to the ith root,
      !  we enter the main loop of refinement by Newton's method.

    1 Continue
        p1=1.d0
        p2=0.d0
        Do j=1,n             !  Loop up the recurrence relation to get the Leg
         p3=p2               !  endre polynomial evaluated at z.
         p2=p1
         p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
        End do
        !  p1 is now the desired Legendre polynomial.
        !  We next compute pp, its derivative, by a standard relation involving
        !  also p2, the polynomial of one lower order.
        pp=n*(z*p1-p2)/(z*z-1.d0)
        z1=z
        z=z1-p1/pp           !  Newton's method.
        if(abs(z-z1).gt.EPS) go to 1

        x(i)=xm-xl*z         !  Scale the root to the desired interval,
        x(n+1-i)=xm+xl*z     !  and put in its symmetric counterpart.
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)  !  Compute the weight
        w(n+1-i)=w(i)                    !  and its symmetric counterpart.
      End do

      End SUBROUTINE gauleg


!=======================================================================
    SUBROUTINE gauss(k,x,w)
!=======================================================================
!   Looks up the values of gaussian coordinates and guassian weights
!   for k-point gaussian quadrature over the interval [0, 1].
!
!   on entry:  k     the number of points in the quadrature
!   -------
!
!   on exit:  x(i)   Gaussian coordinates of the points
!   -------   w(i)   Gaussian weight to point x(i)
!
!   Restriction:     1<= k <= 15
!   -----------
!
!   This is an old version, see also gauleg.f90 for the general case.
!-----------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: k
    REAL(8), INTENT(OUT) :: x(k), w(k)

    IF (k < 1 .OR. k > 15)   &
       Stop 'Error in GAUSS: number of points is out of range '

    SELECT CASE ( k )
    CASE (1)
      x(1) = .5d0
      w(1) = 1.d0
    CASE (2)
      x(1) = .211324865405187d0
      x(2) = .788675134594813d0
      w(1) = .5d0
      w(2) = .5d0
    CASE (3)
      x(1) = .112701665379258d0
      x(2) = .5d0
      x(3) = .887298334620742d0
      w(1) = .277777777777778d0
      w(2) = .444444444444444d0
      w(3) = .277777777777778d0
    CASE (4)
      x(1) = .0694318442029737d0
      x(2) = .330009478207572d0
      x(3) = .669990521792428d0
      x(4) = .930568155797026d0
      w(1) = .173927422568727d0
      w(2) = .326072577431273d0
      w(3) = .326072577431273d0
      w(4) = .173927422568727d0
    CASE (5)
      x(1) = .046910077030668d0
      x(2) = .230765344947158d0
      x(3) = .5d0
      x(4) = .769234655052842d0
      x(5) = .953089922969332d0
      w(1) = .118463442528095d0
      w(2) = .239314335249683d0
      w(3) = .284444444444444d0
      w(4) = .239314335249683d0
      w(5) = .118463442528095d0
    CASE (6)
      x(1) = .033765242898424d0
      x(2) = .169395306766868d0
      x(3) = .380690406958402d0
      x(4) = .619309593041598d0
      x(5) = .830604693233132d0
      x(6) = .966234757101576d0
      w(1) = .0856622461895852d0
      w(2) = .180380786524069d0
      w(3) = .233956967286346d0
      w(4) = .233956967286346d0
      w(5) = .180380786524069d0
      w(6) = .0856622461895852d0
    CASE (7)
      x(1) = .0254460438286207d0
      x(2) = .129234407200303d0
      x(3) = .297077424311301d0
      x(4) = .5d0
      x(5) = .702922575688699d0
      x(6) = .870765592799697d0
      x(7) = .974553956171379d0
      w(1) = .0647424830844348d0
      w(2) = .139852695744638d0
      w(3) = .19091502525256d0
      w(4) = .208979591836735d0
      w(5) = .19091502525256d0
      w(6) = .139852695744638d0
      w(7) = .0647424830844348d0
    CASE (8)
      x(1) = .0198550717512319d0
      x(2) = .101666761293187d0
      x(3) = .237233795041835d0
      x(4) = .408282678752175d0
      x(5) = .591717321247825d0
      x(6) = .762766204958164d0
      x(7) = .898333238706813d0
      x(8) = .980144928248768d0
      w(1) = .0506142681451881d0
      w(2) = .111190517226687d0
      w(3) = .156853322938944d0
      w(4) = .181341891689181d0
      w(5) = .181341891689181d0
      w(6) = .156853322938944d0
      w(7) = .111190517226687d0
      w(8) = .0506142681451881d0
    CASE (9)
      x(1) = .015919880246187d0
      x(2) = .0819844463366821d0
      x(3) = .193314283649705d0
      x(4) = .337873288298095d0
      x(5) = .5d0
      x(6) = .662126711701904d0
      x(7) = .806685716350295d0
      x(8) = .918015553663318d0
      x(9) = .984080119753813d0
      w(1) = .0406371941807872d0
      w(2) = .0903240803474287d0
      w(3) = .130305348201468d0
      w(4) = .156173538520001d0
      w(5) = .16511967750063d0
      w(6) = .156173538520001d0
      w(7) = .130305348201468d0
      w(8) = .0903240803474287d0
      w(9) = .0406371941807872d0
    CASE (10)
      x(1) = .0130467357414141d0
      x(2) = .0674683166555077d0
      x(3) = .160295215850488d0
      x(4) = .283302302935376d0
      x(5) = .425562830509184d0
      x(6) = .574437169490816d0
      x(7) = .716697697064624d0
      x(8) = .839704784149512d0
      x(9) = .932531683344492d0
      x(10)= .986953264258586d0
      w(1) = .0333356721543441d0
      w(2) = .0747256745752903d0
      w(3) = .109543181257991d0
      w(4) = .134633359654998d0
      w(5) = .147762112357376d0
      w(6) = .147762112357376d0
      w(7) = .134633359654998d0
      w(8) = .109543181257991d0
      w(9) = .0747256745752903d0
      w(10)= .0333356721543441d0
    CASE (11)
      x(1) = .0108856709269715d0
      x(2) = .0564687001159523d0
      x(3) = .134923997212975d0
      x(4) = .240451935396594d0
      x(5) = .365228422023827d0
      x(6) = .5d0
      x(7) = .634771577976172d0
      x(8) = .759548064603406d0
      x(9) = .865076002787025d0
      x(10)= .943531299884048d0
      x(11)= .989114329073028d0
      w(1) = .0278342835580868d0
      w(2) = .0627901847324523d0
      w(3) = .0931451054638672d0
      w(4) = .116596882295995d0
      w(5) = .131402272255123d0
      w(6) = .13646254338895d0
      w(7) = .131402272255123d0
      w(8) = .116596882295995d0
      w(9) = .0931451054638672d0
      w(10)= .0627901847324523d0
      w(11)= .0278342835580868d0
    CASE (12)
      x(1) = .00921968287664038d0
      x(2) = .0479413718147626d0
      x(3) = .115048662902848d0
      x(4) = .206341022856691d0
      x(5) = .31608425050091d0
      x(6) = .437383295744266d0
      x(7) = .562616704255734d0
      x(8) = .68391574949909d0
      x(9) = .793658977143309d0
      x(10)= .884951337097152d0
      x(11)= .952058628185237d0
      x(12)= .99078031712336d0
      w(1) = .0235876681932559d0
      w(2) = .0534696629976592d0
      w(3) = .0800391642716731d0
      w(4) = .101583713361533d0
      w(5) = .116746268269177d0
      w(6) = .124573522906701d0
      w(7) = .124573522906701d0
      w(8) = .116746268269177d0
      w(9) = .101583713361533d0
      w(10)= .0800391642716731d0
      w(11)= .0534696629976592d0
      w(12)= .0235876681932559d0
    CASE (13)
      x(1) = .00790847264070593d0
      x(2) = .041200800388511d0
      x(3) = .099210954633345d0
      x(4) = .17882533027983d0
      x(5) = .275753624481777d0
      x(6) = .384770842022433d0
      x(7) = .5d0
      x(8) = .615229157977567d0
      x(9) = .724246375518223d0
      x(10)= .82117466972017d0
      x(11)= .900789045366655d0
      x(12)= .958799199611489d0
      x(13)= .992091527359294d0
      w(1) = .0202420023826579d0
      w(2) = .0460607499188642d0
      w(3) = .0694367551098937d0
      w(4) = .0890729903809729d0
      w(5) = .103908023768444d0
      w(6) = .113141590131449d0
      w(7) = .116275776615437d0
      w(8) = .113141590131449d0
      w(9) = .103908023768444d0
      w(10)= .0890729903809729d0
      w(11)= .0694367551098937d0
      w(12)= .0460607499188642d0
      w(13)= .0202420023826579d0
    CASE (14)
      x(1) = .00685809565159383d0
      x(2) = .0357825581682132d0
      x(3) = .0863993424651175d0
      x(4) = .156353547594157d0
      x(5) = .242375681820923d0
      x(6) = .340443815536055d0
      x(7) = .445972525646328d0
      x(8) = .554027474353672d0
      x(9) = .659556184463945d0
      x(10)= .757624318179077d0
      x(11)= .843646452405843d0
      x(12)= .913600657534882d0
      x(13)= .964217441831787d0
      x(14)= .993141904348406d0
      w(1) = .0175597301658759d0
      w(2) = .0400790435798801d0
      w(3) = .0607592853439516d0
      w(4) = .0786015835790968d0
      w(5) = .092769198738969d0
      w(6) = .102599231860648d0
      w(7) = .107631926731579d0
      w(8) = .107631926731579d0
      w(9) = .102599231860648d0
      w(10)= .092769198738969d0
      w(11)= .0786015835790968d0
      w(12)= .0607592853439516d0
      w(13)= .0400790435798801d0
      w(14)= .0175597301658759d0
    CASE (15)
      x(1) = .00600374098975728d0
      x(2) = .031363303799647d0
      x(3) = .0758967082947864d0
      x(4) = .137791134319915d0
      x(5) = .214513913695731d0
      x(6) = .302924326461218d0
      x(7) = .399402953001283d0
      x(8) = .5d0
      x(9) = .600597046998717d0
      x(10)= .697075673538782d0
      x(11)= .785486086304269d0
      x(12)= .862208865680085d0
      x(13)= .924103291705214d0
      x(14)= .968636696200353d0
      x(15)= .993996259010243d0
      w(1) = .0153766209980586d0
      w(2) = .0351830237440541d0
      w(3) = .053579610233586d0
      w(4) = .0697853389630772d0
      w(5) = .083134602908497d0
      w(6) = .0930805000077812d0
      w(7) = .0992157426635559d0
      w(8) = .101289120962781d0
      w(9) = .0992157426635559d0
      w(10)= .0930805000077812d0
      w(11)= .083134602908497d0
      w(12)= .0697853389630772d0
      w(13)= .053579610233586d0
      w(14)= .0351830237440541d0
      w(15)= .0153766209980586d0
   END SELECT
   END SUBROUTINE gauss


!=======================================================================
      Subroutine INTERV ( xt, lxt, x, left, mflag )
!=======================================================================
!
!     Computes  left = max( i ; 1 <= i <= lxt  .and.  xt(i) <= x )
!     which is the interval containing x .
!
!     A reformatted version of the de Boor routine
!
!      on entry
!      --------
!
!      xt  a real sequence, of length lxt, assumed to be nondecreasing
!      lxt number of terms in the sequence xt .
!      x   the point whose location with respect to the sequence xt is
!          to be determined.
!
!      on exit
!      -------
!      left, mflag   integers, whose values are
!         1     -1   if               x .lt.  xt(1)
!         i      0   if   xt(i)  .le. x .lt. xt(i+1)
!        lxt     1   if  xt(lxt) .le. x
!
!      In particular,  mflag = 0 is the 'usual' case.  mflag .ne. 0
!      indicates that  x  lies outside the halfopen interval
!      xt(1) .le. y .lt. xt(lxt) . the asymmetric treatment of the
!      interval is due to the decision to make all pp functions cont-
!      inuous from the right. (left - ?)
!
!      ...  m e t h o d  ...
!
!  The program is designed to be efficient in the common situation where
!  it is called repeatedly, with  x  taken from an increasing or decrea-
!  sing sequence. this will happen, e.g., when a pp function is to be
!  graphed. the first guess for  left  is therefore taken to be the val-
!  ue returned at the previous call and stored in the  l o c a l  varia-
!  ble  ilo . a first check ascertains that  ilo .lt. lxt (this is nec-
!  essary since the present call may have nothing to do with the previ-
!  ous call). then, if  xt(ilo) .le. x .lt. xt(ilo+1), we set  left =
!  ilo  and are done after just three comparisons.
!     otherwise, we repeatedly double the difference  istep = ihi - ilo
!  while also moving  ilo  and  ihi  in the direction of  x , until
!                      xt(ilo) .le. x .lt. xt(ihi) ,
!  after which we use bisection to get, in addition, ilo+1 = ihi .
!  left = ilo  is then returned.
!-----------------------------------------------------------------------

      integer left,lxt,mflag,   ihi,ilo,istep,middle
      double precision x,xt(lxt)
      data ilo /1/

      ihi = ilo + 1
      if (ihi .lt. lxt)                 go to 20
         if (x .ge. xt(lxt))            go to 110
         if (lxt .le. 1)                go to 90
         ilo = lxt - 1
         ihi = lxt

   20 if (x .ge. xt(ihi))               go to 40
      if (x .ge. xt(ilo))               go to 100

!     ... now x .lt. xt(ilo) . decrease  ilo  to capture  x

   30 istep = 1
   31    ihi = ilo
         ilo = ihi - istep
         if (ilo .le. 1)                go to 35
         if (x .ge. xt(ilo))            go to 50
         istep = istep*2
                                        go to 31
   35 ilo = 1
      if (x .lt. xt(1))                 go to 90
                                        go to 50

!     ... now x .ge. xt(ihi) . increase  ihi  to capture  x

   40 istep = 1
   41    ilo = ihi
         ihi = ilo + istep
         if (ihi .ge. lxt)              go to 45
         if (x .lt. xt(ihi))            go to 50
         istep = istep*2
                                        go to 41
   45 if (x .ge. xt(lxt))               go to 110
      ihi = lxt

!     ... now xt(ilo) .le. x .lt. xt(ihi) . narrow the interval

   50 middle = (ilo + ihi)/2
      if (middle .eq. ilo)              go to 100

!     note. it is assumed that middle = ilo in case ihi = ilo+1 .

      if (x .lt. xt(middle))            go to 53
         ilo = middle
                                        go to 50
   53    ihi = middle
                                        go to 50

!    ... set output and return

   90 mflag = -1
      left = 1
                                        return
  100 mflag = 0
      left = ilo
                                        return
  110 mflag = 1
      left = lxt
                                        return
      End Subroutine INTERV


!======================================================================
      Integer Function Ipointer (n,ii,i)
!======================================================================
!     Determines position of an element i in the integer array ii
!     A zero value indicates that the element is absent
!----------------------------------------------------------------------
      Integer, intent(in) :: n,ii(n),i

      ipointer = 0
      Do j = 1,n
        if (ii(j).eq.i) then; ipointer=j; Return; end if
      End do

      End Function Ipointer


!======================================================================
      Integer Function Rpointer (n,rr,r)
!======================================================================
!     Determines position of an element i in the set ii
!     A zero value indicates that element is not a member of the set
!----------------------------------------------------------------------
      Integer, intent(in) :: n
      Real(8), intent(in) :: r, rr(n)

      Rpointer = 0
      Do j = 1,n
        if (rr(j).eq.r) then; Rpointer=j; Return; end if
      End do

      End Function Rpointer


!======================================================================
      Integer Function Apointer (n,AA,a)
!======================================================================
!     Determines position of an element i in the set ii
!     A zero value indicates that element is not a member of the set
!----------------------------------------------------------------------
      Integer, intent(in)      :: n
      Character(*), intent(in) :: a,AA(n)

      Apointer = 0
      Do j = 1,n
        if (AA(j).eq.a) then; Apointer=j; Return; end if
      End do

      End Function Apointer


!======================================================================
      Subroutine Read_rarg(name,rvalue)
!======================================================================
!     read real argument as name=... from the command line
!----------------------------------------------------------------------
      Implicit none
      Character(*) :: name
      Real(8) :: rvalue
      Integer :: iarg,i,i1,i2,iname
      Character(80) :: AS
!      Integer, internal :: IARGC

      iarg = IARGC(); if(iarg.eq.0) Return
      iname=LEN_TRIM(name)
      Do i=1,iarg
       Call GETARG(i,AS)
       i1=INDEX(AS,'=')+1; i2=LEN_TRIM(AS)
       if(AS(1:i1-2).ne.name(1:iname)) Cycle
       read(AS(i1:i2),*) rvalue; Exit
      End do

      End Subroutine Read_rarg


!======================================================================
      Subroutine Read_iarg(name,ivalue)
!======================================================================
!     read integer argument as name=... from the command line
!----------------------------------------------------------------------
      Implicit none
      Character(*) :: name
      Integer :: ivalue
      Integer :: iarg,i,i1,i2,iname
      Character(80) :: AS
!      Integer, internal :: IARGC

      iarg = IARGC(); if(iarg.eq.0) Return
      iname=LEN_TRIM(name)
      Do i=1,iarg
       Call GETARG(i,AS)
       i1=INDEX(AS,'=')+1; i2=LEN_TRIM(AS)
       if(AS(1:i1-2).ne.name(1:iname)) Cycle
       read(AS(i1:i2),*) ivalue; Exit
      End do

      End Subroutine Read_iarg

!======================================================================
      Subroutine Read_aarg(name,avalue)
!======================================================================
!     read character argument as name=... from the command line
!----------------------------------------------------------------------
      Implicit none
      Character(*) :: name, avalue
      Integer :: iarg,i,i1,i2,iname
      Character(80) :: AS
!      Integer, internal :: IARGC

      iarg = IARGC(); if(iarg.eq.0) Return
      iname=LEN_TRIM(name)
      Do i=1,iarg
       Call GETARG(i,AS)
       i1=INDEX(AS,'=')+1; i2=LEN_TRIM(AS)
       if(AS(1:i1-2).ne.name(1:iname)) Cycle
       read(AS(i1:i2),'(a)') avalue; Exit
      End do

      End Subroutine Read_aarg


!======================================================================
      Subroutine Read_iarr(name,na,iarr)
!======================================================================
!     read integer arrray from string:  name=a,b,c-d,f...
!----------------------------------------------------------------------
      Implicit None
      Character(*) :: name
      Integer :: na,iarr(na)
      Integer :: iarg,ia,iname,i,i1,i2,j,j1,j2,k,k1,k2
      Character(180) :: AS
!      Integer, internal :: IARGC

      iarg = IARGC(); if(iarg.eq.0) Return
      iname=LEN_TRIM(name)
      k=0
      Do i=1,iarg
       Call GETARG(i,AS)
       i1=INDEX(AS,'=')+1; i2=LEN_TRIM(AS)
       if(AS(1:i1-2).ne.name(1:iname)) Cycle
       k=1; Exit
      End do
      if(k.eq.0) Return

      ia=0; j1=i1; ! iarr=0
      Do
       j2=INDEX(AS(j1:i2),',')
       if(j2.eq.0) then; j2=i2; else; j2=j2+j1-1; end if
       j=0 ! INDEX(AS(j1:j2),'-'); k = j+j1-1
       if(j.eq.0) then
        ia=ia+1; if(ia.gt.na) Stop 'Read_iarr: ia > na'
        read(AS(j1:j2),*) iarr(ia)
       else
        read(AS(j1:k-1),*) k1
        read(AS(k+1:j2-1),*) k2
        Do k=k1,k2
         ia=ia+1; if(ia.gt.na) Stop 'Read_iarr: ia > na'
         iarr(ia)=k
        End do
       end if
       j1=j2+1;
       if(j1.gt.i2) Exit
      End do

      End Subroutine Read_iarr


!======================================================================
      Subroutine Read_name(name)
!======================================================================
!     read "name" from the command line as first argument without "="
!----------------------------------------------------------------------
      Implicit none
      Character(*) :: name
      Integer :: iarg,i,i1,i2,iname
      Character(80) :: AS
!      Integer, internal :: IARGC

      iarg = IARGC(); if(iarg.eq.0) Return
      Do i=1,iarg
       Call GETARG(i,AS)
       if(INDEX(AS,'=').ne.0) Cycle
       name=AS
       Exit
      End do

      End Subroutine Read_name


!======================================================================
      Subroutine Read_ipar(nu,name,ivalue)
!======================================================================
!     read integer variable 'ivalue' with identifier 'name'
!     from unit 'nu', where the record like    name = ...
!     is supposed to exist; if absent - ivalue is not changed
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu
      Character(*), intent(in) :: name
      Integer :: ivalue
      Character(80) :: AS
      Integer :: i

      i=LEN_TRIM(name)
      rewind(nu)
    1 read(nu,'(a)',end=2) AS
      if(AS(1:i).ne.name) go to 1
      if(AS(i+1:i+1).ne.'='.and.AS(i+1:i+1).ne.' ') go to 1
      i=INDEX(AS,'=')+1
      if(i.eq.1) go to 1
      read(AS(i:),*) ivalue
    2 Continue

      End Subroutine Read_ipar

!======================================================================
      Subroutine Read_rpar(nu,name,rvalue)
!======================================================================
!     read real variable 'rvalue' with identifier 'name'
!     from unit 'nu', where the record begining with  name = ...
!     is supposed to exist; if absent - rvalue is not changed
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu
      Character(*), intent(in) :: name
      Real(8) :: rvalue
      Character(80) :: AS
      Integer :: i

      i=LEN_TRIM(name)
      rewind(nu)
    1 read(nu,'(a)',end=2) AS
      if(AS(1:i).ne.name) go to 1
      if(AS(i+1:i+1).ne.'='.and.AS(i+1:i+1).ne.' ') go to 1
      i=INDEX(AS,'=')+1
      if(i.eq.1) go to 1
      read(AS(i:),*,end=2,err=2) rvalue
    2 Continue

      End Subroutine Read_rpar


!======================================================================
      Subroutine Read_rval(nu,name,rvalue)
!======================================================================
!     read real variable 'rvalue' with identifier 'name'
!     from unit 'nu', where the record containing  name = ...
!     is supposed to exist; if absent - rvalue is not changed
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu
      Character(*), intent(in) :: name
      Real(8) :: rvalue
      Character(80) :: AS
      Integer :: i

      rewind(nu)
    1 read(nu,'(a)',end=2) AS
      i=INDEX(AS,name)
      if(i.eq.0) go to 1
      i=INDEX(AS,'=')+1
      if(i.eq.1) go to 1
      read(AS(i:),*) rvalue
    2 Continue

      End Subroutine Read_rval


!======================================================================
      Subroutine Read_apar(nu,name,avalue)
!======================================================================
!     read character variable 'avalue' with identifier 'name'
!     from unit 'nu', where the record begining with  name = ...
!     is supposed to exist; if absent - avalue is not changed
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu
      Character(*), intent(in) :: name
      Character(*) :: avalue
      Character(180) :: AS
      Integer :: i

      i=LEN_TRIM(name)
      rewind(nu)
    1 read(nu,'(a)',end=2) AS
      if(AS(1:i).ne.name(1:i)) go to 1
      if(AS(i+1:i+1).ne.'='.and.AS(i+1:i+1).ne.' ') go to 1
      i=INDEX(AS,'=')+1
      if(i.eq.1) go to 1
      read(AS(i:180),*) avalue
    2 Continue

      End Subroutine Read_apar


!======================================================================
      Subroutine Read_iarray(nu,name,nv,ivalue)
!======================================================================
!     read integer array 'ivalue(1:nv)' with identifier 'name'
!     from unit 'nu', where the record like    name = ...
!     is supposed to exist; if absent - ivalue is not changed
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu, nv
      Character(*), intent(in) :: name
      Integer, intent(out) :: ivalue(nv)
      Character(80) :: AS
      Integer :: i,j

      rewind(nu)
    1 read(nu,'(a)',end=2) AS
      i = INDEX(AS,name)
      if(i.eq.0) go to 1
      i=INDEX(AS,'=')+1
      if(i.eq.1) go to 1
      read(AS(i:),*) (ivalue(j),j=1,nv)
    2 Continue

      End Subroutine Read_iarray


!======================================================================
      Integer Function Ifind_position(nu,name)
!======================================================================
!     find position of line with "name" in the begining
!----------------------------------------------------------------------
      Implicit none
      Integer, Intent(in) :: nu
      Character(*), Intent(in) :: name
      Character(80) :: AS
      Integer :: i,j

      Ifind_position = 0
      i=LEN_TRIM(name); j=0
      rewind(nu)
    1 read(nu,'(a)',end=2) AS
      j=j+1
      if(AS(1:i).ne.name) go to 1
      Ifind_position = j
      Backspace(nu)
      Return
    2 rewind(nu)

      End Function Ifind_position


!======================================================================
      Integer Function Jfind_position(nu,name)
!======================================================================
!     find fisrt record in file containing "name"    not-finished yet?
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu
      Character(*), intent(in) :: name
      Character(80) :: AS
      Integer :: i,j

      Jfind_position = 0
      i=LEN_TRIM(name); j=0
      rewind(nu)
    1 read(nu,'(a)',end=2) AS
      j=j+1
      if(AS(1:i).ne.name) go to 1
    2 Jfind_position = j
      Backspace(nu)

      End Function Jfind_position


!======================================================================
      Subroutine Read_string(nu,name,avalue)
!======================================================================
!     read character variable 'avalue' with identifier 'name'
!     from unit 'nu', where the record begining with  name = ...
!     is supposed to exist; if absent - avalue is not changed.
!     Avalue may contains blanks.
!----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: nu
      Character(*), intent(in) :: name
      Character(*) :: avalue
      Character(180) :: AS
      Integer :: i,j

      i=LEN_TRIM(name)
      rewind(nu)
    1 read(nu,'(a)',end=2) AS
      if(AS(1:i).ne.name(1:i)) go to 1
      i=INDEX(AS,'=')+1
      j=LEN_TRIM(AS)
      avalue = AS(i:j)
    2 Continue

      End Subroutine Read_string



!======================================================================
      Real(8) Function RRTC()
!======================================================================
!     give the running time in seconds
!     (some old version to adjust for different compilers)
!----------------------------------------------------------------------
      Implicit none

!     CHARACTER(LEN=8) :: D
!     CHARACTER(LEN=10) :: T
!     INTEGER :: id,ih,im,is,ims

      Real :: TM(2), ETIME !, DTIME

! ... Power station Fortran 4.0:

!     Call DATE_AND_TIME(date=D,time=T)
!     read(D,'(6x,i2)') id
!     read(T,'(3i2,1x,i3)') ih,im,is,ims
!     RRTC = id*86400 + ih*3600 + im*60 + is
!     RRTC = RRTC + ims/1000.d0

! ... Digital Fortran 6.0:

      RRTC = ETIME(TM); ! RRTC = TM(1)

      End Function RRTC



!======================================================================
      SUBROUTINE SPLIN3 (N, X, Y, B, C, D)
!======================================================================
!     DEFINES COEFFICIENTS B(I), C(I), AND D(I), I=1,2,...,N
!     FOR A CUBIC INTERPOLATING SPLINE
!
!     S(X) = Y(I) + B(I)*(X-X(I)) + C(I)*(X-X(I))**2 + D(I)*(X-X(I))**3
!
!     FOR  X(I) .LE. X .LE. X(I+1)
!
! INPUT:
!
!     N = THE NUMBER OF DATA POINTS OR KNOTS (N=>2)
!     X = THE ABSCISSAS OF THE KNOTS IN STRICTLY INCREASING ORDER
!     Y = THE ORDINATES OF THE KNOTS
!
! OUTPUT:
!
!     B, C, D  -  ARRAYS OF SPLINE COEFFICIENTS AS DEFINED ABOVE.
!
!     Y(I) = S(X(I))
!     B(I) = S'(X(I))
!     C(I) = S''(X(I))/2
!     D(I) = S'''(X(I))/6  (DERIVATIVE FROM THE RIGHT)
!
!     ACCOMPANYING FUNCTION  'SEVAL'  CAN BE USED TO EVALUATE THE SPLINE.
!----------------------------------------------------------------------
      Implicit none
      Integer, Intent(in ) :: N
      Real(8), Intent(in ) :: X(n),Y(n)
      Real(8), Intent(out) :: B(n),C(n),D(n)
      Integer :: I,J,M
      Real(8) :: T

      IF ( N .LT. 2 ) RETURN

      IF ( N .LT. 3 ) THEN ! liniar interpolation for n = 2 :

       B(1) = (Y(2)-Y(1))/(X(2)-X(1)); C(1)=0.d0; D(1)=0.d0
       B(2) = B(1);  C(2)=0.d0; D(2)=0.d0
       RETURN

      END IF

      M = N-1

! ... SET UP TRIDIAGONAL SYSTEM

! ... B = DIAGONAL, D = OFFDIAGONAL, C = RIGHT HAND SIDE.

      D(1) = X(2) - X(1)
      C(2) = (Y(2) - Y(1))/D(1)
      Do I = 2, M
         D(I) = X(I+1) - X(I)
         B(I) = 2*(D(I-1) + D(I))
         C(I+1) = (Y(I+1) - Y(I))/D(I)
         C(I) = C(I+1) - C(I)
      End do

! ... END CONDITIONS. THIRD DERIVATIVES AT  X(1)  AND  X(N)
! ... OBTAINED FROM DIVIDED DIFFERENCES

      B(1) = -D(1)
      B(N) = -D(N-1)
      C(1) = 0.
      C(N) = 0.
      IF ( N .GT. 3 ) THEN
       C(1) = C(3)/(X(4)-X(2)) - C(2)/(X(3)-X(1))
       C(N) = C(N-1)/(X(N)-X(N-2)) - C(N-2)/(X(N-1)-X(N-3))
       C(1) = C(1)*D(1)**2/(X(4)-X(1))
       C(N) = -C(N)*D(N-1)**2/(X(N)-X(N-3))
      END IF

!     FORWARD ELIMINATION

      Do I = 2, N
         T = D(I-1)/B(I-1)
         B(I) = B(I) - T*D(I-1)
         C(I) = C(I) - T*C(I-1)
      End do

! ... BACK SUBSTITUTION

      C(N) = C(N)/B(N)
      Do J = 1, M
       I=N-J; C(I) = (C(I) - D(I)*C(I+1))/B(I)
      End do

! ... C(I) IS NOW THE SIGMA(I) OF THE TEXT

! ... COMPUTE POLYNOMIAL COEFFICIENTS

      B(N) = (Y(N) - Y(M))/D(M) + D(M)*(C(M) + 2.*C(N))
      Do I = 1, M
         B(I) = (Y(I+1) - Y(I))/D(I) - D(I)*(C(I+1) + 2.*C(I))
         D(I) = (C(I+1) - C(I))/D(I)
         C(I) = 3.*C(I)
      End do
      C(N) = 3.*C(N);  D(N) = D(N-1)

      End Subroutine SPLIN3


!=======================================================================
      Real(8) Function SEVAL(N,U,X,Y,B,C,D)
!=======================================================================
!     Evaluates the value of cubic spline for abscisa U
!-----------------------------------------------------------------------
      Implicit none
      Integer, intent(in) :: N
      Real(8), intent(in) :: U
      Real(8), intent(in) :: X(N),Y(N),B(N),C(N),D(N)
      Integer :: mflag, I
      Real(8) :: DX

      Call INTERV (X,N,U,I,mflag)

      DX=U-X(I); SEVAL=Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))

      End function SEVAL




!=====================================================================
      Real(8) Function xAy (n,k,array,type,x,y)
!=====================================================================
!
!     Returns   <x | array |y>  for different array types
!
!     The original size of array is n*n. The bandwidth of array is 2k-1
!     for the non-symmetric case, and is k for symmetric case.
!     The lower-band column storage mode is assumed.
!---------------------------------------------------------------------
      Integer, intent(in) :: n, k
      Character, intent(in) :: type
      Real(8), intent(in) :: array(n,*), x(n), y(n)
      Integer :: i

      xAy = 0.d0

      if ( type .eq. 'f') then
       Do i=1,n
        xAy = xAy + SUM(array(1:n,i)*x(1:n)) * y(i)
       End do
      else
       Stop  'xAy: not yet supported array type'
      end if

      End Function xAy



!====================================================================
      Module rk4_data
!====================================================================
!     contains a set of coefficients ordering accordint to
!     four indexes
!--------------------------------------------------------------------
      Implicit none

      Integer :: nrk = 0       ! current number of coefficients
      Integer :: mrk = 0       ! maximum dimension
      Integer :: irk = 2**16   ! initial dimension

! ... coefficients:

      Real(8), allocatable :: crk(:)

! ... their attributes:

      Integer, allocatable :: kr1(:),kr2(:),kr3(:),kr4(:)

      End Module rk4_data


!======================================================================
      Subroutine alloc_rk4_data(m)
!======================================================================
! ... allocate, deallocate, or reallocate the arrays in "rk4_data"
!----------------------------------------------------------------------
      Use rk4_data

      Implicit none
      Integer :: m
      Integer, allocatable :: iar(:)
      Real(8), allocatable :: rar(:)

      if(m.le.0) then
       if(allocated(crk)) Deallocate (crk,kr1,kr2,kr3,kr4)
       nrk = 0; mrk = 0
      elseif(.not.allocated(crk)) then
       mrk = m; nrk = 0
       Allocate(crk(mrk),kr1(mrk),kr2(mrk),kr3(mrk),kr4(mrk))
      elseif(m.le.mrk) then
       Return
      elseif(nrk.eq.0) then
       Deallocate (crk,kr1,kr2,kr3,kr4)
       mrk = m
       Allocate(crk(mrk),kr1(mrk),kr2(mrk),kr3(mrk),kr4(mrk))
      else
       Allocate(rar(nrk))
       rar=crk; Deallocate(crk); Allocate(crk(m)); crk(1:nrk)=rar
       Deallocate(rar)
       Allocate(iar(nrk))
       iar=kr1; Deallocate(kr1); Allocate(kr1(m)); kr1(1:nrk)=iar
       iar=kr2; Deallocate(kr2); Allocate(kr2(m)); kr2(1:nrk)=iar
       iar=kr3; Deallocate(kr3); Allocate(kr3(m)); kr3(1:nrk)=iar
       iar=kr4; Deallocate(kr4); Allocate(kr4(m)); kr4(1:nrk)=iar
       Deallocate(iar)
       mrk = m
       write(*,'(a,i10/a,f10.2,a)') ' Realloc_rk4_data: new dimension = ', m, &
        ' memory requred = ', 24.d0*m/(1024*1024),'  Mb'
      end if

      End Subroutine alloc_rk4_data


!======================================================================
      Subroutine Add_rk4_data(k1,k2,k3,k4,C)
!======================================================================
!     add new data to the list
!----------------------------------------------------------------------
      Use rk4_data

      Implicit none
      Integer, intent(in) ::  k1,k2,k3,k4
      Real(8), intent(in) :: C
      Integer :: i,k,l,m

      if(mrk.eq.0) Call alloc_rk4_data(irk)

! ... search position (k) for new integral

      k=1; l=nrk
    1 if(k.gt.l) go to 2
      m=(k+l)/2
      if    (k1.lt.kr1(m)) then;       l = m - 1
      elseif(k1.gt.kr1(m)) then;       k = m + 1
      else
       if    (k2.lt.kr2(m)) then;      l = m - 1
       elseif(k2.gt.kr2(m)) then;      k = m + 1
       else
        if    (k3.lt.kr3(m)) then;     l = m - 1
        elseif(k3.gt.kr3(m)) then;     k = m + 1
        else
         if    (k4.lt.kr4(m)) then;    l = m - 1
         elseif(k4.gt.kr4(m)) then;    k = m + 1
         else
          crk(m)=crk(m)+C
          Return ! the same integral
         end if
        end if
       end if
      end if
      go to 1
    2 Continue

! ... shift the rest data up:

      Do i=nrk,k,-1
       m = i + 1
       crk(m)=crk(i)
       kr1(m)=kr1(i); kr2(m)=kr2(i); kr3(m)=kr3(i); kr4(m)=kr4(i)
      End do

! ... add new integral:

      crk(k)=C; kr1(k)=k1; kr2(k)=k2; kr3(k)=k3; kr4(k)=k4; nrk=nrk+1
      if(nrk.eq.mrk) Call alloc_rk4_data(mrk+irk)

      End Subroutine Add_rk4_data


!======================================================================
      Subroutine Check_rk4_data
!======================================================================
!     check the same integrals,   ???
!----------------------------------------------------------------------
      Use rk4_data

      Implicit none
      Integer :: i,j

      Do i=2,nrk; j=i-1
       if(kr1(i).eq.kr1(j).and.kr2(i).eq.kr2(j).and.kr3(i).eq.kr3(j)) &
        kr4(i)=-kr4(i)
      End do

      End Subroutine Check_rk4_data

