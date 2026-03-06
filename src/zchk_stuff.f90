      module FLY_PARAMETERS

!     make sure that NPTSF and nebinit are the same

      integer, parameter     :: numpp   = 4500 !3200 
      integer, parameter     :: nebinit = 500  
      integer, parameter     :: NMTout   = 500
      integer, parameter     :: mtran    = 100000 !50000
      integer, parameter     :: miso     = 49 !100
      integer, parameter     :: nlvi     = 170 !50
      end module FLY_PARAMETERS

!*--------------------------------------------------------------------72

      module POP_STUFF

!____POPSTUF

      use FLY_PARAMETERS

!  popR: pop() = THE ARRAY OF POPULATIONS THAT IS RETURNED FROM SOLVERS.  AT
!                ANY TIME THERE ARE NEQ POPULATIONS.
!                NOTE!  (neq .LE. numpop) DUE TO IP DEPRESSION.

!      double precision :: pop
!      common /popR/ pop(numpp)

!  popul8: popout() = THE ARRAY OF POPULATIONS INTO WHICH THE pop ARRAY IS MOVED
!                     AND THEN ADJUSTED: 1) PLACED IN CORRECT numpop LOCATION;
!                     2) SMALL NUMBERS SET TO 10^(-15) OF MAXIMUM; 3) NEGATIVE
!                     POPULATION ARE REPORTED

      double precision :: popout

      common /popul8/ popout(numpp)

!  compactI,R =  INDEXING INFORMATION FOR CALCULATED POPULATIONS ARRAY popout

      double precision :: efgrdp
      double precision :: evtoip
      double precision :: efromgr

      common /compactR/efromgr(numpp),evtoip(numpp),efgrdp(numpp)
      common /cindipd/indp2l(0:numpp),indl2p(0:numpp)
      common /cindipd2/indipd(numpp),indlev(100,0:100),indpop(0:numpp),indlev_z(100,0:100)

      integer:: nlev
      integer:: nlbn
      integer:: iptlev

      common/isocom/nlev(0:100),nlbn(0:100),iptlev(0:100)

      integer:: ntot
      integer:: isolev
      integer:: ionstg
      integer:: levl

      common/ilevls/ntot,isolev(numpp),ionstg(numpp),levl(numpp)

      double precision :: elev
      double precision :: glev
      double precision :: qlev
      double precision :: eltoip

      common/levels/elev(numpp),glev(numpp),qlev(numpp),eltoip(numpp)

      double precision, dimension(500) :: energyz
      
      common/edata/energyz
 
      integer:: ntran
      integer:: lupper
      integer:: llower 
      integer:: ntran_z ! mscho ddd 201027
      integer:: ltrans
      integer:: ltype
      common/rindx/ntran,lupper(mtran),llower(mtran),ltrans(mtran),ltype(mtran),ntran_z(mtran)
!     common/rindx/ntran,ltrans(mtran),ltype(mtran),llower_z(mtran),ntran_z(mtran)
      end module POP_STUFF

!*--------------------------------------------------------------------72

      module FofE_STUFF

      use FLY_PARAMETERS

!  velocity    = INITIAL ENERGY GRID, INTEGRATION WEIGHTS, AND CURRENT fe

!    energy()  = ENERGY GRID FOR ELECTRONS DISTRIBUTIONS
!    denergy() = WEIGHTS SUCH THAT SUM fe*sqrt(energy)*denergy OVER nebins = 1.
!    fe()      = ELECTRON ENERGY DISTRIBUTION DEFINED SO THAT Ne*fe*sqrt(energy)*de
!                IS THE NUMBER OF ELECTRONS IN ENERGY RANGE DE

      double precision  :: energy
      double precision  :: denergy
      double precision  :: fe
      

      common /velocity/ energy(0:nebinit),denergy(0:nebinit),fe(0:nebinit)

!  fetimesR,I = ELECTRON ENERGY DISTRIBUTION

!    fe0()    = ELECTRON DISTRIBUTION AT tlast, SEE COMMON/last/ IN timestuf
!    fe1()    = ELECTRON DISTRIBUTION AT tnext, SEE COMMON/next/ IN timestuf
!    feIN()   = ELECTRON DISTRIBUTION READ FROM fefile AT TIMES fetime()
!    feOut()  = CALCULATED ELECTRON DISTRIBUTION AT OUTPUT TIMES
!    fetime() = ARRAY OF TIMES AT WHICH feIN() DATA IS READ
!    ntimefe  = # of TIMES IN FILE fefile
!    nebins   = # OF ENERGIES IN ARRAY energy

      double precision  :: fe0
      double precision  :: fe1
      double precision  :: feIn
      double precision  :: feOut
      double precision  :: fetime

      common /fetimesR/ fe0(nebinit) ,fe1(nebinit), feIn(nebinit,NMTout),feOut(nebinit,NMTout),fetime(NMTout)

      common /fetimesI/ ntimefe,nebins

!  fefiler  = INFORMATION ON fe FILE OPTIONS

!    fefile = FILE WITH THE ELECTRON DISTRIBUTION AS A FUNCTION OF ENERGY & TIME
!    feflag = FLAG INDICATING METHOD OF GENERATING fe(): 'file', 'fefile'

      character (len = 20)   :: fefile
      character (len =  8)   :: feflag

      common /fefiler/ fefile,feflag


      double precision :: zqs, zqr

      common /xesource/zqs(nebinit),zqr(nebinit)

      double precision :: exthr

      common /xthr/ exthr

      end module FofE_STUFF

!*--------------------------------------------------------------------72

      module ATOMIC_and_PLASMA_STUFF

      use FLY_PARAMETERS

!   plasmar = INDEPENDENT VARIABLES OF THE CALCULATION (REALS)

!    t      = ELECTRON TEMPERATURE IN KELVIN
!    denow  = THE ELECTRON DENSITY (/CC)
!    dinow  = THE ION DENSITY (/CC) of popout variable: may not be the total density
!    atomz  = ATOMIC NUMBER OF THE SPECIES OF INTEREST
!    azp1   = ATOMZ+1
!    azm1   = ATOMZ-1
!    roott  = SQRT(tev)
!    boltzt = BOLTZMANN CONSTANT TIMES TEMPERATURE IN KELVIN (ERGS)
!    tev    = electron temperature in eV
!    debye  = THE IONIZATION POTENTIAL DEPRESSION OF ION WITH UNIT CHARGE (eV)
!    zave   = THE AVERAGE ION CHARGE, SUCH THAT zave*denow = ion density
!    c1     = 2.*(2.*pi*emass*boltzt/hp**2)**1.5 USED IN sAHA-bOLTZMANN RELATIONSHIP

      double precision  :: denow
      double precision  :: dinow
      double precision  :: atomz
      double precision  :: azp1
      double precision  :: azm1
      double precision  :: boltzt
      double precision  :: tev
      double precision  :: r_ion
      double precision  :: r_debye
      double precision  :: r_lattice
      double precision  :: zave
      double precision  :: c1

      common /plasmar/ denow,dinow,atomz,azp1,azm1,boltzt,tev,r_ion,r_debye,r_lattice,zave,c1

!  plasmai = INDEPENDENT VARIABLES OF THE CALCULATION (INTEGERS)

!    iatomz = atomz
!    iazp1  = azp1
!    iazm1  = azm1
!    ibias  -> ibias+1 IS THE INDEX OF THE FIRST IP OF THE ELEMENT atomz

      integer:: iatomz,iazp1,iazm1,ibias

      common /plasmai/ iatomz,iazp1,iazm1,ibias

      double precision :: eipz
      common /stage/ eipz(0:miso)

!    il        = NUMBER OF LOWEST ION STAGE TO BE EVER CONSIDERED (SET TO 1)
!    in        = NUMBER OF THE LAST ION STAGE TO BE TREATED AS A GROUND STATE ONLY
!    ilast     = NUMBER OF LOWEST ION STAGE AT CALCULATED TEMPERATURE & DENSITY.
!                ilast IS NOT ALWAYS EQUAL TO il DUE TO IP DEPRESSION!
!    ionow     = NUMBER OF ION STAGES AT CALCULATED TEMPERATURE AND DENSITY
!    neq       = TOTAL NUMBER OF STATES AT CALCULATED TEMPERATURE AND DENSITY. THIS IS
!                THE ORDER OF THE MATRIX SOLVE.
!    i1        = INDEX OF THE PART OF FILLED RATE MATRIX CONTAINING GROUND STATES ONLY
!    ire()     = NUMBER OF STATES REMOVED FROM EACH ION STAGE DUE TO IP DEPRESSION
!    ibo()     = NUMBER OF BOUND STATES REMAINING IN EACH ION STAGE AFTER THE IP
!                DEPRESSION HAS BEEN TAKEN INTO ACCOUNT: ibo() = nboun()-ire()
!    nboun()   = NUMBER OF STATES EXCLUDING AUTOIONIZING STATES:
!                nboun() = nlevel()-nauto()
!    irelast() = A STORAGE OF PREVIOUS ire(). uSED TO INDICATE THAT IF NUMBER OF
!                STATES CHANGES ONE MUST REDISTRIBUTE THE POPULATIONS AND START
!                THE CALCULATIONS WITH NEW neq.

      integer:: il,in,ilast,ionow,neq,i1
      integer:: ire,ibo,nboun,irelast,iraut,iaut

      common /inhab0/ il,in,ilast,ionow,neq,i1,ire(0:miso),ibo(0:miso),nboun(0:miso),irelast(0:miso),iraut(0:miso),iaut(0:miso)

!  tiinfo = THE ION TEMPERATURE INFORMATION
!
!   tiev = ION TEMPERATURE AT CURRENT TIME IN eV
!   tifixed =  FIXED ION TEMPERATURE SET ON INPUT (eV)
!   tibyte = ION TO ELECTRON TEMPERATURE RATIO TO DEFINE Ti FROM Te

      double precision  :: tiev,tifixed,tibyte

      common /tiinfo/ tiev,tifixed,tibyte

!  const = PHYSICAL CONSTANTS USED IN THE CALCULATIONS

      double precision  :: ao,boltz,speed,rh,pi,charge
      double precision  :: hbar,emass,amass,avogad,hp
      double precision  :: evtohz,rhtohz,rhtoev,evtoerg,hceverg

      common /const/ ao,boltz,speed,rh,pi,charge,hbar,emass,amass,avogad,hp,evtohz,rhtohz,rhtoev,evtoerg,hceverg


      character(len = 3):: flyflag
      common /ksh_flag/flyflag


      double precision :: Vcaa, Vchy, Vche, Vcli, Vchl
      common /vrates/Vcaa(nlvi,nlvi,0:miso),Vchy(25,25),Vche(31,31),Vcli(40,40),Vchl(70,70,3)	

      character (len=8) :: initflag,radflag,tiflag
      character (len=20):: initfile
      common /inits/ initflag, initfile, radflag, tiflag

      end module ATOMIC_and_PLASMA_STUFF


      module TIMING_STUFF

!     dtfly = time-increments dictated by FLYCHK
!     isfly = maximum number of time steps 

      double precision  :: tfirst
      double precision  :: tzelda
      double precision  :: dtfly
      integer :: isfly
      common /timefly/ tfirst,tzelda,dtfly,isfly

      double precision :: total_time
      common /timezeld/total_time

!     ee_scale: multiplicative factor for e-e collision rate.
!     Default = 1.0 (physical).
!     For Fig. 1(c) of Cho et al. PRE 109, 045207 (2024):
!       Set ee_scale ~ 100.0 and use the EEDF output of the <600fs run as
!       the initial fe input. The simulation time is then ~100x shorter
!       but linearly rescaled in the plot to show thermalization to ~20 ps.
!     This reproduces the non-equilibrium-to-Maxwellian relaxation figure.
      double precision :: ee_scale
      common /eescale/ ee_scale

      end module TIMING_STUFF

