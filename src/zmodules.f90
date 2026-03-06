
!                       MODULE definitions for Zelda

!                     Copyright (C) W.L. Morgan (1996)
!                        Kinema research & Software

!----------------------------------------------------------------------

      module PARAMETER_DEFS

!       NSP = number of different chemical species that the code can handle
!      NLEV = maximum number of states or levels allowed for each species
!     NMPTS = maximum number of points for momentum transfer cross section
!     NQPTS = maximum number of points for inelastic cross sections

      integer, parameter :: NSP   = 30
      integer, parameter :: NLEV  = 30
      integer, parameter :: NMPTS = 120
      integer, parameter :: NQPTS = 120

      end module PARAMETER_DEFS

!*--------------------------------------------------------------------72

      module EEDF_DIMENSIONS

!     NPTSF = number of points in finite-differenced electron energy
!             distribution function

!     make sure that NPTSF and nebinit are the same
      integer, parameter :: NPTSF=500, NPTSF1=NPTSF+1, NPTSF2=2*NPTSF1

      end module EEDF_DIMENSIONS

!*--------------------------------------------------------------------72

      module CONSTANTS

!     Some Constants:

      double precision, parameter :: Pi = 3.141592d0
      double precision, parameter :: Q0 = 1.0d-16         ! convert: 1 A^2->10^-16 cm^2
      double precision, parameter :: emass = 9.10956d-28  ! electron mass in gm
      double precision, parameter :: emassamu = 5.487d-4  !    "      "   "  amu
      double precision, parameter :: ergeV = 1.602192d-12 ! ergs/eV
      double precision, parameter :: EJeV = 1.602192d-19  ! Joules/eV
      double precision, parameter :: Bk1 = 8.6169d-5      ! Boltzmann's const in eV / K
      double precision, parameter :: a0 = 0.529d-8        ! Bohr radius in cm
      double precision, parameter :: Ryd = 13.6048        ! Rydberg constant in eV

      end module CONSTANTS

!*--------------------------------------------------------------------72

      module EFLAGS

      logical :: eeflag     = .FALSE.  ! no electron-electron collisions
      logical :: eiflag     = .FALSE.  ! no electron-ion      collisions
      logical :: recflag    = .FALSE.  ! no dissociative recombination
      logical :: fflag      = .FALSE.  ! use M-B dist. @ t=0 with Te = tel0
      logical :: esflag     = .FALSE.  ! no secondary e- from ionization
      logical :: ebflag     = .FALSE.  ! no external electron source
      logical :: ifpr       = .TRUE.   ! print f(E)
      logical :: ifpl       = .TRUE.   ! plot  f(E)
      logical :: rt_flag    = .FALSE.  ! do not write rate coeffs vs. t
      logical :: ft_flag    = .FALSE.  ! do not write f(E) vs. t
      logical :: enbl_flag  = .FALSE.  ! do not write energy balance fraction
      logical :: dv_flag    = .TRUE.   ! include all delta(v) vib transitions
      logical :: Cycle_flag = .FALSE.  ! do not calc for an integral number
                                       ! of AC cycles
      logical :: dt_flag    = .FALSE.  ! use default time step

      logical :: err_flag   = .FALSE.  ! error flag

      logical :: stand_alone = .TRUE.  ! stand-alone Boltzmann calculation

      logical :: ES_call    = .FALSE.
      logical :: FI_call    = .FALSE.
      logical :: REC_call   = .FALSE.
      logical :: EN_call    = .FALSE.

!       other flags:
      integer :: isec_flag  = 0        ! 1=> sec e- in bin 1; 2=> uniform

      character (len=2) :: enflag      ! 'DC', 'AC', or 'TB'

      end module EFLAGS

!*--------------------------------------------------------------------72

      module RUN_PARAMETERS

!       physical parameters:
      double precision    :: tgas     = 300.0      ! gas temperature in K
      double precision    :: press    =   1.0      ! gas pressure in Atm
      double precision    :: gdens    = 2.69d19   ! gas density (/cc)
      double precision    :: tel0     =   1.0      ! initial electron temperature in eV
      double precision    :: edens0   =   1.0      ! initial electron density (/cc)
      double precision    :: ion_density           ! ion density if different from Ne
      double precision    :: ave_Z    =   1.0      ! ave charge state of ions in plasma
      double precision    :: power    =   0.0      ! no photons
      double precision    :: photflux =   0.0      ! no photons
      double precision    :: ephoton  =   0.0      ! no photons

      double precision    :: ebyn0                 ! E/N or peak E/N for AC case
      double precision    :: freq                  ! AC frequency
      double precision    :: phase                 ! AC phase; 0 => E/N(t)=E/N(0)sin(ft)
      integer :: N_AC_Cycles           ! number of AC cycles

!       calculational parameters:
      integer :: ismax  = 100          ! integrate for 100 time steps
      integer :: iprint =   5          ! print every 5 time steps
      double precision    :: eps    = 1.d-6        ! convergence criterion
      double precision    :: dt     = 1.d-6        ! integration time step
      double precision    :: deltau
      double precision    :: umax

      end module RUN_PARAMETERS

!*--------------------------------------------------------------------72

      module FILENAMES_and_UNITS

!         file names:
        character (len=30)  :: run_file = ' '    ! run parameters
        character (len=30)  :: q_file   = ' '    ! cross sections
        character (len=30)  :: out_file = ' '    ! standard output file

        character (len=30)  :: ebyn_file    ! E/N(t)
        character (len=30)  :: fofe_file    ! f(E) @ t = 0
        character (len=30)  :: eb_file      ! external electron source
        character (len=30)  :: roft_file    ! rates vs. t
        character (len=30)  :: foft_file    ! f(E)  vs. t

!         unit numbers:
!           input:
        integer    :: aux_in   =   4   ! auxiliary input files
        integer    :: key_in   =   5   ! keyboard  input
        integer    :: std_in   =  15   ! primary   input file

!           output:
        integer    :: crt_out  =   6   ! screen  output
        integer    :: plt_out  =   7   ! plot    output file
        integer    :: t_out    =   8   ! time dependent output
        integer    :: f_out    =   9   ! f(E) vs. t output
        integer    :: mbr_plt  =  12   ! MB rates for plotting
        integer    :: mbr_out  =  13   ! MB rates output
        integer    :: err_out  =  14   ! error message file "errors.txt"
        integer    :: std_out  =  16   ! primary output file
        integer    :: acc_out  =  20   ! AC coefficients vs. t output
                                       ! "ac_coeff.out"

!         The following output is no longer used; the same information
!         appears in "roft_file" on unit t_out = 8:
        integer    :: acr_out  =  21   ! AC rates        vs. t output
                                       ! "ac_rates.out"


      end module FILENAMES_and_UNITS

!*--------------------------------------------------------------------72
!                The following are not used at present:
!*--------------------------------------------------------------------72

      module DISTRIBUTION_FUNCTION

      use EEDF_DIMENSIONS

      integer                                :: N_e_bins ! number of
                                                         ! energy bins

      type EEDF1                    ! EEDF class definition

        real              :: energy ! electron energy
        real              :: edf    ! f(E)

      end type EEDF1

      type (EEDF1) f0 (0:NPTSF1)    ! distribution function objects

      type EEDF2                    ! EEDF class definition

        double precision  :: energy ! electron energy
        double precision  :: edf    ! f(E)

      end type EEDF2

      type (EEDF2) f  (0:NPTSF1)

      end module DISTRIBUTION_FUNCTION

!*--------------------------------------------------------------------72

      module LARGE_ARRAYS

      use EEDF_DIMENSIONS    ! in the future use dynamic dimensioning

      double precision, dimension (NPTSF,NPTSF)   :: cf

      double precision, dimension (NPTSF,NPTSF)   :: aee

!      double precision, dimension (:,:), allocatable   :: cf

!      double precision, dimension (:,:), allocatable   :: aee

      end module LARGE_ARRAYS

!*--------------------------------------------------------------------72


