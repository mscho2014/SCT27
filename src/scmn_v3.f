
      program scfly

c     This program is written by Hyun-Kyung Chung and Richard W. Lee.
c     It solves rate equations using a set ofsuper-configuration 
c     for a given plasma conditions in order to obtain 
c     charge state distributions and level population.
c     It requires one input file called "atomic.inp".
c
c     version 1.6-1.8 uses different photoionization c.x. from other versions.
c     temperature calculations for XFEL near threshold are very low. 
c
c*****************************
c*                           *
c*  d e c l a r a t i o n s  *
c*                           * 
c*****************************
c
      implicit real*8 (a-h,o-z)
      
      include 'mainstuf'
      include 'runstuf'
      include 'timestuf'
      include 'popstuf'
      include 'xc_stuf'

      character*1  iblank
      character*8  blank
      common /blkchr/ blank,iblank
      common /idatnum/ idat(2)
      common /eescale/ ee_scale
      common /en5sp/ en5sdata(84)
      common /exp1dat/ exp1data(336)
      common /exp1datI/ Iexp1dat(1)
      common /xexe1/ xex1data(336)
      common /xexe1I/ iex1data(1)
c 
      common /popwriteR/ bpop(numpp,NMTout),biso(0:100,NMTout),deninit
      common /paramnow/ trnow,tinow,zbarnow,totpnow,sizenow,rhonow
c
c     necessary for self-consistent Te calculations ! hydrodynamic expansion
      common /vhydro/dtev1,dtev2,zbar1,dedt1,dedt2,teint1,teint2,ehydro
      common /ihydro/ihydflag
      common /absorb/abb(100),abf(100),aff(100),att(100)

      dimension itestneg(numpp),barjtr(100)

      data ioo/6/, ioi/5/

      data jacobon/.false./

      ! initial setting, not a constant. It would be changed during calculation. 
      data ioutflag /2/
      data iatomz / 0/
      data nametime/'time0'/, nameinfo/ 'info0'/,trfile/'    '/
      data tr / 0.00/, trev / 0.00/, dilution/ 1.00/
      data tibyte/1.00/ , tiev/ 0.00/
      data tiflag/'off'/, initflag/'ss'/, radflag/'off'     /
      data fhflag/'off'/, feflag/'off'/, fefile/'                    '/
      data opacflag/'off'/
      data initfile/'    '/
      data r_lattice /0.d0/
      data ihydflag /0/

      data hp,evtohz,rhtohz,rhtoev,evtoerg,hceverg
     1    /6.626176d-27,2.417969d+14,3.289842d+15,1.3605804d+01,
     1     1.602d-12,7.7515d-20/
      data ao,boltz,speed,rh,pi,charge,hbar,emass,amass,avogad/
     1    5.2917706d-09,1.380662d-16,2.997924580d+10,1.097373d+05,
     2    3.141592653598,4.803242d-10,1.0545887d-27,9.109534d-28,
     3    1.660565d-24,6.022045d+23/
c
      character*25 filnam
      logical initdone, around

      data relflag/.true./

      common /codename/ version
      character*20 version 
      data version /'                    '/ 
      
c     For zelda run
      data zeldaflag/ .true./
      data isfly /20/ 

c     ee_scale: e-e collision rate multiplier for ZELDA.
c     Default 1.0 = physical. Set ~100.0 to reproduce PRE 2024 Fig.1(c)
c     (thermalization study with accelerated e-e collisions).
      ee_scale = 1.0d0

c*******************
c*  begin program  *
c*******************
c
c___print the version, and
c___get the input parameters from the user.
c
      write (ioo,5)
    5 format(' FLYCHK Super configuration with Te : Sep 2015')
      version = 'SCFLY-Te version 1.9 : SCT27 test version'
      write (ioo,10)
   10 format (1x)

      ! mscho 200511 : iteflag 1->0 (please confirm why : 201218 mscho)   
      iteflag=1    ! default =1 for self-consistent temperature calculations
                   ! add 0 at the line of 'z' input
      ibdr=1
      fbdr=1.d0
      isamp(1)=0  ! Gaunt Factor Fit Excitation
      isamp(2)=1  ! Brugess & Chidichimo Ionization
      nxprim=10   ! Up to n=10 states
      do i=1,12
         ratcon(i)=1.d0  ! rate control
      enddo
c      ratcon(10)=1.e-5  ! three-body recombination

   50 call getpar


c     call secnds to get the start time for this run

      SecZero = timing(0.0)

      write(filnam,801)nametime
 801  format('zb.',a20)
      open(23,file=filnam)

      write(filnam,805)nametime
 805  format('pp.',a20)
      open(22,file=filnam)

      write(filnam,803)nametime
 803  format('te.',a20)
      open(27,file=filnam)  ! write out tev as a function of t

      write(filnam,'("abs.",a20)')nametime
      open(28,file=filnam)
      write(28,'("iso iup ilo energy pop(lo) opl0  opl tot ")')

      write(filnam,'("tx.",a20)')nametime
      open(29,file=filnam)
      write(29,'("     time,       Te,      Ne,      Ni,    Size,",
     &     " hv(max),  J(max),",
     &     " Opacity[cm-1], BB,      BF,      FF")')

      write(filnam,'("ja.",a20)')nametime
      open(31,file=filnam)
      write(31,'("# of photon grid = ",i5,
     &     " at following energies")') ngrid
      write(31,'(10f11.3)')(xev(i),i=1,ngrid)
      write(31,'("time[s],    Te[eV],  Ne[cm-3],  Ni[cm-3],Size[cm]")')
      write(31,'("opacity[cm-1] at the photon energies")')


      write(filnam,'("jt.",a20)')nametime
      open(32,file=filnam)
      write(32,'(2i5)') ntimes,ngrid
      write(32,'(10f11.3)')(xev(i),i=1,ngrid)


      write(filnam,'("pm.",a20)')nametime
      open(98,file=filnam)

      write(filnam,'("fe.",a20)')nametime
      open(99,file=filnam)

      open(unit=77,file='rate_monitor',status='replace')
      close(77)
   


      open(unit=92,file='source_monitor',status='replace')
      write(92,*)'i source(i) sink(i)'
      close(92)
   
      open(unit=93,file='zelda.out',status='replace')
       write(93,*)'The code zelda output' 
      close(93)
     

      open(unit=79,file='ipd_monitor',status='replace')
      close(79)
      
      open(unit=78,file='pop_monitor',status='replace')
      close(78)

      open(unit=76,file='pe_monitor',status='replace')
       write(76,*)'eii,elev(j),eipz(iz),debye(iatomz-iz+1),elev(i)' 
       write(76,*) 'j, i, k, iz, field_e, xc_bcontxx'
      close(76)



      ! If you want to monitor some variables, then please input this one 
      ! open(unit=xx,file='xxx',status='old',position='append')
      !   write(xx,*)'xxxxxxx'
      ! close(xx) 


c***********************
c*  for a given ion:   *
c***********************
c
c___define element dependent variables.
c
      atomz  = iatomz

      iazp1  = iatomz+1
      azp1   = iazp1
      iazm1  = iatomz-1
      azm1   = iazm1
      ibias  = iatomz*iazm1/2


c ---------------------------------------------------------------
c     (1) determine the flags and indices for the calculations
c
c     (2) setup a level structure
c      
c     isos is the ***  input parameter ********
c     isomin = the lowest iso-sequence to have hydrgenic levels
c     isomax = the highest iso-sequence to have hydrgenic levels
c     All the other ions are ground-state only
c
c     (3) write out level energies and transition rates
c----------------------------------------------------------

      call atolev

c    define the minimum temperature allowable for this ion to stop overflows

      tmin  = eipz(iatomz)/300.
      tmin  = min(0.1d0,tmin)
      trmin = eipz(iatomz)/3.

c
c  Do the special case of the initial condition.  
c  There is a differnt method for each type of HISTORY option and 
c  for each type of INITIAL option.  
c  This leads to a large number of initialization paths
c
      tnext = timezero

c   define the limits of the write counter

      iolimit = 10
      if(numdens.ne.0)iolimit=min(numdens,10)

      if (initflag .eq. 'file') then
         initdone = .false.
         call loadinit(initdone)
         if (.not. initdone) then
            write(ioo,*)' $$$ Failure in INITLOAD try again'
            write(16,*)' Warning:Failure in INITLOAD try again'
            stop
         endif
      endif

      it=1 ! mscho explain: 1st time step 

      call tprep(it)

      write(16,'("")')
      write(ioo,'(" $$$ Working on # ",i3," of ",i3,": time ",1pe9.2,
     1     " with end at ",e9.2)') it,ntimes,tnext,timestop
      write(16,'(" Working on # ",i3," of ",i3,": time ",1pe9.2,
     1     " with end at ",e9.2)') it,ntimes,tnext,timestop

      ratioNZ = 1.00

          call worknlte_init(ioff,it,ratioNZ) !mscho 210218
          call depress(ioff) !mscho 
          call caleipd(S_E,P_F) !mscho

      if(ioff.eq.1) then
         write(ioo,'(" *** INITIAL TIME NOT working")')
         write(16,'(" Warning:INITIAL TIME NOT working")')
         stop
      endif
      
      deninit=totpnow

      do k=1,numpop
         itestneg(k) = 0.00
      enddo
      itestn = 0
      popmax = 1.00d-70
      do 4400 i=1,numpop
         popmax = max(popmax,popout(i))
         if (popout(i) .lt. 0.00) then
            itestn = itestn+1
            itestneg(itestn) = i
            go to 4400
         endif
 4400 continue

      popmin = 1.00d-70

      if (itestn .gt. 0) then
         write(16,4600) tnext,(itestneg(i),i=1,itestn)
 4600    format(/,' Warning:Time = ',1pe9.2,0p,
     1        ' Pop<0 for indices:',5i4,/,(5x,20i4))
      endif

      
      iolow = 1
      iocount = 1
      ratioNZ=1.00
      iohigh = iocount+iolow-1

      call getbpop(1,ratioNZ,popmin)


      if (iocount .eq. iolimit .or. ntimes .eq. 1)then
         ! 'this is writeout printed, but total tim 
         !           loop (ntimes) is 1. Very special case' 
         call writeout(iolow,iohigh)
         iolow = iohigh+1
         iocount = 0
      endif
     

c     write out atomic data

      tev=tevm(1)
      call wratmdat(sizenow,ratioNZ)

c ***************************************************
c   put calculated electron distribution fe -> feOut:
c
      if(feflag .ne.'off') then
         do i=1,nebins
            feOut(i,1) = fe(i)
         enddo
      endif

c ***************************************************
c     Initialize the next time-steps

      tlast = tnext
      t0 = tev
      ti0 = tiev
      rh0 = rhoAtT
      if(iruntyp1 .eq. 'nt'.and. percent.gt.0.d0) then
         totp0 =totpAtT
      else
         totp0 = deninit
      endif

      d0 = denow
      size0 = sizenow
      zb0 = zbAtT
      tr0 = trnow
      dlt0= dltAtT
      zmix0=zmixAtT
      perc0=percAtT

      if (radflag .eq. 'trfile') then
         do i=1,ngrid
            barj0(i) = barjAtT(i)  ! mscho : this is the variable for the radiation 
         enddo
      endif
      if(feflag.ne.'off') then
         fhot0=fhot
         do i=1,nebins
            fe0(i) = fe(i)
         enddo
      endif

c     
c************************************************
c*  loop over all times for the transient case  *
c************************************************
c     
c   find time for the intialization
      
      SecEnd(1) = timing(0.0)
      SecEndout = SecEnd(1)-SecZero

      write(16,'(" INITIALIZATION time",1pe10.2,/,1x,"Tev Ne ",
     1        2e10.2,3x,"delta t",0p,f9.3)')
     2        tnext,t0,d0,SecEndout
      if (iolevel .eq. 1 .and. initflag .eq. 'ss') then
         call subwraa(t0,d0,tnext)
      endif
 
      do i=1,iatomz
         irelast(i) = -1
      enddo

c -----------------------------------------------------
c      
c     time LOOP
c
c -----------------------------------------------------

      ratioNZ = 1.00
      istate = 2

      do ii=1, ntem
         fnem0(ii)=fnem(ii)
         tevm0(ii)=tevm(ii)
      enddo

      if(ngrid.ne.0)call abscrx(ratioNZ,0)

      write(32,'(1p,2e11.3)')times(1), 0.d0
      do i=1, ngrid
         barjtr(i)=barj1(i)*exp(-att(i)*size1)
      enddo
      write(32,'(1p,10e11.3)')(barjtr(i),i=1,ngrid)

c ... internal energy and total ion density
c     totin = internal energy [eV/cm3]
      totin=0.d0
      sumpp=0.d0
      do i = 1,numpop
         totin = totin + efgrdp(i)*popout(i)
         sumpp = sumpp + popout(i)
      enddo

c ... calculation of energies
c ... initial energy[eV/cm2]= (Te*(Ne+Ni)+Totin)*l[cm]

      tevbal=tevm0(1)
      if (tiflag .eq. 'off') tiev=tevbal

      ephabs=0.d0   ! dE/dt(t)
      dedt1=0.d0
      teint1=totin/sumpp
      tephab0=size0*(1.5*(denow*tevbal+sumpp*tiev)+totin) ! initial energy [eV/cm2]
      write(*,'("Initial Total Energy=(Te*Ne+Ti*Ni)*Size+Int.En:",
     &     1p,e10.2, " [eV/cm2]")')tephab0

      ehydro = 0.d0  ! hydrodynamic energy
      zbar1 = zbarnum(qqq)

      dtev1 = 0.d0
      dtev2 = 0.d0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do 200 it=2,ntimes ! Iteration loop starting point 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
         iocount = iocount+1

         call tprep(it)

         write(ioo,'(" $$$ MWorking on # ",i3," of ",i3,": time ",1pe9.2,
     1        " with end at ",e9.2)') it,ntimes,tnext,timestop
         write(16,'(" Working on # ",i3," of ",i3,": time ",1pe9.2,
     1        " with end at ",e9.2)') it,ntimes,tnext,timestop

c.... TIME INDEPENDENT SOLUTION ................................................         
         IF (ievolve .ne. 'td') THEN

            call worknlte(ioff,it,ratioNZ)

            ratioNZ=1.d0  ! popout is normalized properly 

            if(ioff.ne.0) then
               write(ioo,'(" *** STEADY-STATE NOT working at ",i5)'),it
               write(16,'(" Warning:STEADY-STATE NOT working at ",i5)')
     &              ,it
               stop
            endif

c.... TIME DEPENDENT SOLUTION ..........................................
         ELSE
           

            if(iteflag.eq.0)then
               write(*,*) '****** tdnlte start ********** '
               call tdnlte(ioff,it,ratioNZ)
              
            else
               write(*,*) ' ********** tdtcal start ********* '
               call tdtcal(ioff,it,ratioNZ) 

               if(ngrid.ne.0) call abscrx(ratioNZ,0)

               write(32,'(1p,2e11.3)')times(it), 0.d0
               do i=1, ngrid
                  barjtr(i)=barj1(i)*exp(-att(i)*size1)
               enddo

               write(32,'(1p,10e11.3)')(barjtr(i),i=1,ngrid)

            endif


c           popout is normalized to deninit at it=1
c           bpop will be normalized properly by ratioNZ  
            if(ioff.ne.0) then
               write(ioo,'(" *** TIME-DEPT NOT working at ",i5)'),it
               write(16,'(" Warning:TIME-DEPT NOT working at ",i5)'),it
               stop
            endif

         ENDIF
         
c     find computer time to advance this time step
         
         SecEnd(it) = timing(0.0)
         SecEndout = SecEnd(it)-SecEnd(it-1)


c     call the write out routine for diagnostic information
         
         write(16,'(" CALCULATION time ",1pe10.2,/,1x,"Tev Ne ",
     1        2e10.2,3x,"delta t",0p,f9.3)')
     2        tnext,tevm(1),denow,SecEndout
         if (iolevel .eq. 1) then
            if(iouttime(it) .gt. 0) then
               call subwraa(t1,denow,tnext)
            endif
         endif
         
c     test for negative populations and place in an array for printing
         
         do 420 k=1,100
            itestneg(k) = 0.00
 420     continue
         
         itestn = 0
         popmax = 1.00d-70
         do 440 i=1,numpop
            popmax = max(popmax,popout(i))
            if (popout(i) .lt. 0.00) then
               itestn = itestn+1
               itestneg(itestn) = i
               go to 440
            endif
 440     continue

         popmin = 1.00d-70
         
         if (itestn .gt. 0) then
            write(16,460) tnext,(itestneg(i),i=1,itestn)
 460        format(' Warning:Time = ',1pe9.2,0p,
     1           ' Pop<0 for indices:',5i4,/,(5x,20i4))
         endif
         
c     file array with values of variables used in calculations
!     This part generates "outfile" (cf - write(21,*))        
         iohigh = iocount+iolow-1
         call getbpop(it,ratioNZ,popmin)
         if (iocount .eq. iolimit .or. it .eq. ntimes)then
            call writeout(iolow,iohigh)
            iolow = iohigh+1
            iocount = 0
         endif

c ***************************************************
c   put calculated electron distribution fe -> feOut
c
         if(feflag.ne.'off') then
            do i=1,nebins
               feOut(i,iohigh) = fe1(i)
            enddo

c mscho 200506 : m-start : write out fe() 
            write(99,'("# ",1p,e12.4)')times(it)
            do i=1,nebins
               write(99,'(1p,10e12.4)')energy(i),fe1(i)
            enddo
            write(99,'(/)')
c mscho 200506 : m-finish 

            open(30,file='feout.d')
            nstep=3
            write(30,'(2i5)') nstep,nebins
            write(30,'(10e12.5)') (energy(i),i=1,nebins)
            write(30,'(2e12.5)') times(1)
            write(30,'(10e12.5)') (feOut(i,1),i=1,nebins)
            write(30,'(2e12.5)') times(it)
            write(30,'(10e12.5)') (fe1(i),i=1,nebins)
            write(30,'(2e12.5)') times(ntimes)
            write(30,'(10e12.5)') (fe1(i),i=1,nebins)
            close(30)
         endif
c     write out fe()
c     write out population density 
         open(unit=33,file='popout.d')
         do i=1,numpop
            write(33,8310) kname(indpop(i)),bpop(i,it)
 8310       format(1x,a8,1x,1pe12.4,9e12.4)
         enddo
         close(33)
c     write out rerun file 
c     rerun file should be as a history file 
         open(unit=33,file='rerun')
         write(33,910)times(1), times(it)
         write(33,911) iruntyp1, iruntyp2
         if(iruntyp1.eq.'ne')denrun1=enecalc(it)
         if(iruntyp1.eq.'nt')denrun1=totpcalc(it)
         if(iruntyp1.eq.'rho')denrun1=rhoout(it)
         if(iruntyp2.eq.' ')then
            write(33,912) times(1),teout(1), denrun1
            write(33,912) times(it),t1,denrun1
            write(33,912) times(it)*1.e3,t1,denrun1
         else
            if(iruntyp2.eq.'nt')denrun2=totpcalc(it)
            if(iruntyp2.eq.'rho')denrun2=rhoout(it)
            write(33,912) times(1),teout(1), denrun1,denrun2
            write(33,912) times(it),t1,denrun1,denrun2
            write(33,912) times(it)*1.e3,t1,denrun1,denrun2
         endif
         close(33)

 910     format('c  time=',1p,2e12.5)
 911     format('time      te',10x,a3,10x,a3)
 912     format(1p,e10.3,0p,2(1x,e12.4))
c ***************************************************
c....  Now place the 'next' variables in the 'last' variables

         tlast = tnext
         t0 = t1
         rh0 = rh1
         totp0 = totp1
         d0 = denow
         write(*,*)'size0, size1, sizenow', size0, size1, sizenow
         ! result is all zero (size0, size1, sizenow) - 210219 (mscho), why? 

         size0 = size1
         zb0 = zbarnum(qqq)
         ti0 = ti1
         tr0 = tr1
         dlt0= dlt1
         zmix0=zmix1
         perc0=perc1
         if (radflag .eq. 'trfile') then
            do i=1,ngrid
               barj0(i) = barj1(i)
            enddo
         endif
         do ii=1, ntem
            fnem0(ii)=fnem1(ii)
            tevm0(ii)=tevm1(ii)
         enddo
         if(feflag .ne. 'off') then
            fhot0=fhot1
            do i=1,nebins
               fe0(i)=fe1(i)
            enddo
         endif

      call wratmdat(sizenow,ratioNZ)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 200  continue ! time loop finished
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
 831  format(1x,a8,1x,1pe12.4,9e12.4)
c
      close(32)
      close(31)
      close(29)
      close(28)
      close(27)
      close(26)
      close(24)
      close(23)
      close(21)
      close(98)
      close(99)


      open(unit=94,file='fe5.outfile')
       
      !write(99,'("# ",1p,e12.4)')times(it)
      do i=1,nebins
         write(94,'(1p,10e12.4)')energy(i),fe1(i)
      enddo
      !write(99,'(/)')
      close(94)
c
c___print completion message with output file names.
c
      write (ioo,1000) nametime,nameinfo
      write (16,1000) nametime,nameinfo
 1000 format(/,' output to files: ',a20,'  and ',a20,//)

c     initialize all the parameters
      call initzero
      write(16,'("")')

      stop

      end

      subroutine initzero
      implicit real*8 (a-h,o-z)
      
      include 'mainstuf'
      include 'runstuf'
      include 'timestuf'
      include 'popstuf'
      include 'xc_stuf'

      iatomz=0
      numdens=0
      nametime='time0'
      nameinfo='info0'
      tr= 0.00
      trev=0.00
      dilution=1.00
      tibyte=1.00
      tiev=0.00
      tiflag='off'
      initflag='ss'
      ievolve='ss'
      radflag='off'
      fhflag='off'
      feflag='off'
      opacflag='off'
      itrout='off'
      fefile='                    '
      initfile='    '
      timename='    '
      trfile='    '
      fhotset = 0.00
      thotset = 0.00
      zbartemp=0.d0
      percent=0.d0
      zelse=0.d0
      elseznum=0.d0
      zmix=0.d0
      do k=1, NMTout
         eneout(k)=0.d0
         totpout(k)=0.d0
         zbout(k)=0.d0
         sizeout(k)=0.d0
         trout(k)=0.d0
         tiout(k)=0.d0
      enddo
      return
      end

c      real*8 function second(x)
c      real*4 secnds
c      second = secnds(x)
c      return
c      end

      real*8 function timing(x)
      integer secnds
      timing = secnds(x)
      return
      end

      subroutine getbpop(it,ratioNZ,popmin)
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'timestuf'
      include 'popstuf'
      include 'runstuf'

      common/NLTE_pw/pwrbb,pwrbf,pwrff,pwrtt

      common /popwriteR/ bpop(numpp,NMTout),biso(0:100,NMTout),deninit
      common /paramnow/ trnow,tinow,zbarnow,totpnow,sizenow,rhonow
      common /pwrrates/tpwrbb(NMTout),tpwrbf(NMTout),tpwrff(NMTout),
     &     tpwrtt(NMTout)
      dimension rlte(numpp)

      common /linmax/totmax
      character*30 filiso

      tev=tevm(1)

      trcalc(it) = trnow
      ticalc(it) = tinow
      zbarnow = zbarnum(xxxdum)
      zbcalc(it) = zbarnow
      enecalc(it) = denow
      totpcalc(it) = totpnow
      sizecalc(it) = sizenow
      
c.... Place the populations in the output population array bpop

c     bpop is in the energy-ascending order
c     ratioNZ is needed so that the sum of bpop(i)=totpnow

      sum=0.d0
      ptot=0.d0

      pmax=0.d0
      do ii=1, numpop
         ilv=indpop(ii)
         bpop(ii,it)=max(popout(ilv)*ratioNZ,popmin)
         sum=sum+bpop(ii,it)
         ptot=ptot+popout(ilv)
         pmax=max(pmax,bpop(ii,it))
      enddo
      dinow=ptot
      
      call caleipd(S_E,P_F)
      call pwrloss(ptot)

      const=ratioNZ/ptot ! will make pwrtot in units of [erg/s/atom]
      tpwrbb(it)=pwrbb*const
      tpwrbf(it)=pwrbf*const
      tpwrff(it)=pwrff*const
      tpwrtt(it)=pwrtt*const

      do iz=iatomz, 0, -1
         ilv0=iptlev(iz)-1
         sum=0.d0
         do i=1, nlev(iz)
            ilv=indlev(i,iz)
            ii=ilv0+i
            sum=sum+bpop(ii,it)
         enddo
         biso(iz,it)=sum/totpnow
      enddo


      if (ievolve .eq. 'td') then
         if(it.eq.1)then
            write(22,'("# case     time     Te        Ne",100i9)')
     &           (iz,iz=0,iatomz)
         endif
         write(22,819)it,tlast,teout(it),enecalc(it),
     &        (max(0.d0,biso(iatomz-iz,it)),iz=0,iatomz)

      else
         if(it.eq.1)then
            write(22,'("# case   steady     Te        Ne",100i9)')
     &           (iz,iz=0,iatomz)
         endif
         write(22,819)it,tlast,teout(it),enecalc(it),
     &        (max(0.d0,biso(iatomz-iz,it)),iz=0,iatomz)
      endif

      if(isomin.ne.1.and.biso(isomin,it).gt.1.e-4)then
         write(16,'(" Warning:Needs isomin <",i3," at ",i3)')
     &        isomin,it
      endif
      if(isomax.ne.iatomz.and.biso(isomax,it).gt.1.e-4)then
         write(16,'(" Warning:Needs isomax >",i3," at ",i3)')
     &        isomax,it
      endif

c     check for the double bumps.

      do iso=1, iatomz-1
         if((biso(iso,it).le.biso(iso-1,it)).and.
     &        (biso(iso,it).le.biso(iso+1,it)))then
            if(biso(iso,it).gt.1.e-6)then
               write(ioo,'("    Double bumps ",i5,10f10.6)')
     &              iso,biso(iso-1,it),biso(iso,it),biso(iso+1,it)
               write(16,'(" Warning:Double bumps ",i5,10f10.6)')
     &              iso,biso(iso-1,it),biso(iso,it),biso(iso+1,it)
            endif
         endif
      enddo

      write(16,'(" Info:ratioNZ,dilution ",10f10.6)') ratioNZ,dilution
      
 802  format(f10.1,10(1x,f10.4))
 803  format(1p,10(1x,e12.4))
 816  format(f12.4,1p,e13.4,2x,3i4,2x,2i2.2)
 817  format(i5,2x,1p,e13.4)
 818  format('#',i5,1p,3e10.2,1x,100e9.2)
 819  format(i5,1p,3e10.2,1x,100e9.2)
 820  format(i5,5x,a5,1p,2e10.2,1x,100e9.2)
 821  format('ionz=', i3,1p,10(1x,e10.3))
 822  format(1p,10(1x,e10.3))

      return
      end

      subroutine writeout(iolow,iohigh)
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'timestuf'
      include 'runstuf'
      include 'popstuf'

      common /pwrrates/tpwrbb(NMTout),tpwrbf(NMTout),tpwrff(NMTout),
     &     tpwrtt(NMTout)
      common /popwriteR/ bpop(numpp,NMTout),biso(0:100,NMTout),deninit
      character*8 timelab,telab
      character*8 nelab,nilab,zblab,sizlab,tilab,trlab
      character*8 neclab,niclab,zbclab,sizclab,ticlab,trclab
      data timelab/'Time    '/, telab /'Te      '/
      data nelab  /'Ne in   '/, nilab /'Ni in   '/,zblab /'Zbar in '/
      data sizlab /'Size in '/, tilab /'Ti in   '/,trlab /'Tr in   '/
      data neclab /'Ne Calc '/,niclab /'Ni Calc '/,zbclab /'ZbarCalc'/
      data sizclab/'SizeCalc'/,ticlab /'Ti calc '/,trclab /'Tr Calc '/

      character*8 ntelab,mtelab,mfelab
      data ntelab /'Num Te  '/,mtelab /'Te Hot  '/,mfelab/'Frac Hot'/


      write (21,8305) iolow,iohigh,numpop
      do i=1,numpop
         write(21,8310)
     1        kname(indpop(i)),(bpop(i,k),k=iolow,iohigh)
      enddo
      
      write(21,8310) nelab,(eneout(k),k=iolow,iohigh)
      write(21,8310) nilab,(totpout(k),k=iolow,iohigh)
      write(21,8310) zblab,(zbout(k),k=iolow,iohigh)
      write(21,8310) sizlab,(sizeout(k),k=iolow,iohigh)
      write(21,8310) trlab,(trout(k),k=iolow,iohigh)
      write(21,8310) tilab,(tiout(k),k=iolow,iohigh)
      
      write(21,8310) neclab,(enecalc(k),k=iolow,iohigh)
      write(21,8310) niclab,(totpcalc(k),k=iolow,iohigh)
      write(21,8310) zbclab,(zbcalc(k),k=iolow,iohigh)
      write(21,8310) sizclab,(sizecalc(k),k=iolow,iohigh)
      write(21,8310) trclab,(trcalc(k),k=iolow,iohigh)
      write(21,8310) ticlab,(ticalc(k),k=iolow,iohigh)
      
      write(21,8310) telab,(teout(k),k=iolow,iohigh)
      write(21,8310) timelab,(times(k),k=iolow,iohigh)

      write(21,8310) ntelab,(1.d0*ntemout(k),k=iolow,iohigh)

      nmx=1
      do k=iolow,iohigh
         nmx=max(ntemout(k),nmx)
      enddo
      do i=1, nmx
         write(21,8310) mtelab,(tevmout(i,k),k=iolow,iohigh)
      enddo
      do i=1, nmx
         write(21,8310) mfelab,(fnemout(i,k),k=iolow,iohigh)
      enddo

      if (timename.eq.'grid') then
         if(iolow.eq.1)then
            write(23,'("# case       steady",
     &      "       Te          Ne          zbar",
     &      "       pwrbb        pwrbf        pwrff        pwrtt")')
         endif
         do k=iolow, iohigh
            write(23,804)k,times(k),teout(k),enecalc(k),zbcalc(k),
     &           tpwrbb(k),tpwrbf(k),tpwrff(k),tpwrtt(k)
         enddo
      else
         if(iolow.eq.1)then
            write(23,'("# case       time  ",
     &      "       Te          Ne          zbar",
     &      "       pwrbb        pwrbf        pwrff        pwrtt")')
         endif
         do k=iolow, iohigh
            write(23,804)k,times(k),teout(k),enecalc(k),zbcalc(k),
     &           tpwrbb(k),tpwrbf(k),tpwrff(k),tpwrtt(k)
         enddo
      endif
 802  format(f10.1,10(1x,f10.4))
 803  format(i7,7x,a5,1p,10(1x,e12.4))
 804  format(i7,1p,10(1x,e12.4))
 8305 format (3i10)
 8310 format(1x,a8,1x,1pe12.4,9e12.4)

      return
      end


      subroutine tprep(it)
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'runstuf'
      include 'timestuf'
      include 'popstuf'
      include 'xc_stuf'

      common /paramnow/ trnow,tinow,zbarnow,totpnow,sizenow,rhonow


      tnext = times(it)
      call AtTime(tnext)
      teout(it) = max(teAtT,tmin)
      tiout(it) = tiAtT
      if (tiflag .eq. 'file') then
         if (tiAtT .le. 0.00) then
            tiev = teout(it)
         else
            tiev = tiout(it)
         endif
      else if (tiflag .eq. 'fixed') then
         tiev = tifixed
      else if (tiflag .eq. 'off') then
         tiev = teout(it)
      else if (tiflag .eq. 'ti/te') then
         tiev = teout(it)*tibyte
      endif

      eneout(it) = eneAtT
      sizeout(it) = sizeAtT
      totpout(it) = totpAtT
      rhoout(it) = rhoAtT
      trout(it) = trAtT


      t1 = max(teAtT,tmin)
      ti1 = tiev
      tinow = ti1
      rh1 = rhoAtT
      totp1 = totpAtT
      d1 = eneAtT
      zb1 = zbAtT
         
      ntemout(it)=ntem
      do ii=1, ntem
         fnem1(ii)=fnem(ii)
         tevm1(ii)=tevm(ii)
         fnemout(ii,it)=fnem(ii)
         tevmout(ii,it)=tevm(ii)
      enddo

      if (opacflag .eq. 'file') then
         sizenow = SizeAtT
      else if (opacflag .eq. 'size') then
         sizenow = sizefix
      else
         sizenow = 0.00
      endif
      size1 = sizenow

      if (radflag .eq. 'file' .or. radflag .eq. 'trfile') then
         if (trout(it) .le. trmin) then
            trnow = 0.00
         else
            trnow = trout(it)
         endif
         if (radflag .eq. 'trfile') then
            do i=1,ngrid
               barj1(i) = barjAtT(i)
            enddo
            barjend = 0.00
            if (xbarj .gt. 0.00.and.trnow.gt.0.d0) then
               xbytrev = xev(ngrid)/trnow
               if (xbytrev .lt. 100.) then
                  barjend = barj1(ngrid)/xev(ngrid)**3*exp(xbytrev)
               endif
            endif
         endif
      else if (radflag .eq . 'fixed') then
         trnow = trfixed
      else
         trnow = 0.00
      endif
      tr1 = trnow
      dlt1= dltAtT

      if(feflag.ne.'off') then
         fhot1=fhot
         sum=0.d0
         esum=0.d0
c--------------- in case of fe(i) from ZELDA
         if(it.ne.1.and.zeldaflag)then
            do i=1,nebins
               fe(i)=fe0(i)
            enddo
         endif
c-----------------------
         do i=1,nebins
            fe1(i) = fe(i)
            sum=sum+fe(i)*sqrt(energy(i))*denergy(i)
            esum=esum+fe(i)*sqrt(energy(i))*denergy(i)
     &           *energy(i)
         enddo
         tevfout(it)=2./3.*esum/sum
         fevfout(it)=fhot
      endif

c     new addition for mixture option
      zmix = zmixAtT
      zmix1= zmix
      perc1= percAtT

      percent=perc1

      if(percent.lt.1.d0)then
         zelse = zmix1*(percent/(1.-percent))*amratio
      else
         zelse = 0.d0
      endif

c     use zbout for the mixture option
      zbout(it) = zelse

      return
      end

      subroutine worknlte(ioff,it,ratioNZ)
      implicit real*8 (a-h,o-z)

      include 'mainstuf'
      include 'runstuf'
      include 'timestuf'
      include 'popstuf'
      common /paramnow/ trnow,tinow,zbarnow,totpnow,sizenow,rhonow
      common /outrt/filnn
      character*10 filnn

      external steady,depress,lte

      ioff=0
      write(filnn,'("rate.",i3.3)')it

      znow = zb1
      denow = d1
      totpnow = totpAtT/(1.+percent/(1.-percent)*amratio)
      rhonow = rhoAtT
      sizenow = size1

      tev = t1
      boltzt = boltz*tev*11604.5
      c1 = 2.*(2.*pi*emass*boltzt/hp**2)**1.5

      tiev=ti1
      trev=tr1
      tr = trev*11604.5
      dilution=dlt1
      zmix= zmix1

      if (radflag .eq. 'trfile') then
         do i=1,ngrid
            barjnow(i) = barj1(i)
         enddo
         barjend = 0.00
         if (xbarj .gt. 0.00.and.trev.gt.0.d0) then
            xbytrev = xev(ngrid)/trev
            if (xbytrev .lt. 100.) then
               barjend = barjnow(ngrid)/xev(ngrid)**3
     &              *exp(xbytrev)
            endif
         endif
      endif
      

c it's important to start with the right zave
      zave=0.1d0
      do iz=2, iatomz
         if(tev.ge.eipz(iz)/2..and.tev.lt.eipz(iz-1)/2.) then
            zave=iatomz-real(iz)+1
            goto 10
         endif
      enddo
 10   continue

c     correction on 03/11/2011
c     always initialize

      tev=tevm(1)
      boltzt = boltz*tev*11604.5
      c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
      
      call depress(ioff)

c  HISTORY option NE with either RHO or NT

      if (iruntyp1 .eq. 'ne' .and.
     1   (iruntyp2 .eq. 'rho' .or. iruntyp2 .eq. 'nt')) then

         if (iruntyp2 .eq. 'nt') then
            totpnow = totpAtT/(1.+percent/(1.-percent)*amratio)
         else if (iruntyp2 .eq. 'rho' .and. elseznum .eq. 0.00
     1           .and. percent .ne. 0.00) then
            write(ioo,*)' *** Must give Z of other species'
            write(16,*)' Warning:Must give Z of other species'
            ioff=1
            return
         else
            totpnow =rhoAtT/atmass(iatomz)*(1.-percent)
         endif

         if (initflag .eq. 'file'.and.it.eq.1) then
            dentotal = 0.00
            do i=1,numpop
               dentotal = popinit(i)+dentotal
            enddo
            
            do i=1,numpop
               popout(i) = popinit(i)*totpnow/dentotal
            enddo

         else
            
            tev=tevm(1)
            boltzt = boltz*tev*11604.5
            c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
            
            !write(ioo,*)'depress one time : ms' !ms 191229 - test version * this is not working in my setting (191229- opacity LTE/rho-T plot grid)  
            call depress(ioff)
            if(ioff.eq.1)return

c     if there are no levels only the bare nucleus put population in bare
            if (ionow .eq. 0) then
               popout(indlev(1,ionow)) = totpnow
            else
               if (initflag .eq. 'lte'.or.ievolve.eq.'lte') then
                  call lte(totpnow)
                  zave = zbarneq(xxxdum)
                  tev=tevm(1)
                  boltzt = boltz*tev*11604.5
                  c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
                  call depress(ioff)
                  if(ioff.eq.1)return
                  call lte(totpnow)
               else if (initflag .eq. 'ss'.or.ievolve.eq.'ss') then
                  call avnow(0.d0)
                  call steady(totpnow)
                  zave = zbarneq(xxxdum)

                  tev=tevm(1)
                  boltzt = boltz*tev*11604.5
                  c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
                  call depress(ioff)
                  if(ioff.eq.1)return
                  call steady(totpnow)
c iteration
                  nxx=0
 281              zave2 = zbarneq(xxxdum)
                  xxx=abs(zave2-zave)
                  yyy=min(zave, zave2)
                  if(xxx.gt.yyy.or.xxx.gt.0.2d0)then
                     nxx=nxx+1
                     if(nxx.gt.20)then
                        write(ioo,'(" *** Fail to converge")')
                        write(16,'(" Warning:Fail to converge")')
                        goto 282
                     else if(nxx.gt.10) then
                        if(zave2-zave.lt.0.d0) then
                           !write(ioo,*)'mstest'  !mscho 191229 
                           zave = 0.6*zave+0.4*zave2
                        else
                           !write(ioo,*)'mstest' !mscho 191229
                           zave = 0.4*zave+0.6*zave2
                        endif
                     else
                        !write(ioo,*)'mstest' !mscho 191229
                        zave = 0.5*(zave+zave2)
                     endif
                     tev=tevm(1)
                     boltzt = boltz*tev*11604.5
                     c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
                     call depress(ioff) !ms : 
                     if(ioff.eq.1)return
                     if(ionow.eq.0)then
                        pop(neq) = totpnow !ms : pop(neq) is used in zbarneq(xxxdum) function (191227) 
                        goto 282
                     endif
                     call avnow(0.d0)
                     call steady(totpnow)
                     goto 281
                  else
                     zave=zave2
                     tev=tevm(1)
                     boltzt = boltz*tev*11604.5
                     c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
                     call depress(ioff)
                     if(ioff.eq.1)return
                     if(ionow.eq.0)then
                        pop(neq) = totpnow
                        goto 282
                     endif
                     call avnow(0.d0)
                     call steady(totpnow)
                     goto 282
                  endif
 282              continue
                  if (sizenow .ne. 0.00) then
                     call fill
                     zxxx=zbarnum(qqq)
                     ratioNZ =totpnow/qqq
                     do i=1,numpop
                        popout(i) = popout(i)*ratioNZ
                     enddo
                     call avnow(sizenow)
                     call steady(totpnow)
                  endif
               endif
               call fill
            endif
         endif
        
c.... Test the agreement between the NE and the initial NT
        
         calcNE = totpnow*(zelse+zbarnum(qqq))
         ratioNE = calcNE/denow
         write(ioo,'(" $$$  Ne(calc)/Ne(input) = ",1pe11.2)') ratioNE
         write(16,'(" $$$ Ne(calc)/Ne(input) = ",1pe11.2)') ratioNE
      endif
      


c     HISTORY option NE with neither RHO or NT
      
      if ((iruntyp1 .eq. 'ne') .and. (iruntyp2 .eq. ' ')) then
         
         if (initflag .eq. 'file'.and.it.eq.1) then
            dentotal = 0.00
            do 1199 i=1,numpop
               popout(i) = popinit(i)
               dentotal =popinit(i)+dentotal
 1199       continue
            calcNZ = denow/(zelse+zbarnum(qqq))
            if(calcNZ.gt.1.d30) calcNZ=dentotal
            ratioNZ = calcNZ/dentotal
            totpnow = calcNZ

            do 1198 i=1,numpop
               popout(i) = popinit(i)*ratioNZ
 1198       continue

         else
            
c     determine which levels exist


            totpnow = denow/(zelse+zave)
            tev=tevm(1)
            boltzt = boltz*tev*11604.5
            c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5

            call depress(ioff)
            if(ioff.eq.1)return
            if (ionow .eq. 0) then
               popout(indlev(1,ionow)) = denow/(zelse+atomz)
            else
               if (initflag .eq. 'lte'.or.ievolve.eq.'lte') then
                  call lte(totpnow)
                  zave = zbarneq(xxxdum)
                  totpnow = denow/(zelse+zave)

                  tev=tevm(1)
                  boltzt = boltz*tev*11604.5
                  c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
                  call depress(ioff)
                  if(ioff.eq.1)return
                  if (ionow .eq. 0) then
                     popout(indlev(1,ionow)) = denow/(zelse+atomz)
                  else
                     call lte(totpnow)
                  endif
                  nxx=0
 290              zave2=zbarneq(qqq)
                  xxx=abs(zave2-zave)
                  yyy=min(zave, zave2)
                  if(xxx.gt.yyy.or.xxx.gt.0.2d0)then
                     nxx=nxx+1
                     if(nxx.gt.20)then
                        write(ioo,'(" $$$  Fail to converge")')
                        write(16,'(" Warning:Fail to converge")')
                        goto 300
                     else if(nxx.gt.10) then
                        if(zave2-zave.lt.0.d0) then
                           zave = 0.6*zave+0.4*zave2
                        else
                           zave = 0.4*zave+0.6*zave2
                        endif
                     else
                        zave = 0.5*(zave+zave2)
                     endif
                     totpnow = denow/(zelse+zave)
                     tev=tevm(1)
                     boltzt = boltz*tev*11604.5
                     c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
                     call depress(ioff)
                     if(ioff.eq.1)return
                     if (ionow .eq. 0) then
                        popout(indlev(1,ionow)) = denow/(zelse+atomz)
                     else
                        call lte(totpnow)
                     endif
                     goto 290
                  else
                     zave=zave2
                     totpnow = denow/(zelse+zave)
                     tev=tevm(1)
                     boltzt = boltz*tev*11604.5
                     c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
                     call depress(ioff)
                     if(ioff.eq.1)return
                     if (ionow .eq. 0) then
                        popout(indlev(1,ionow)) = denow/(zelse+atomz)
                     else
                        call lte(totpnow)
                     endif
                     goto 300
                  endif
               else if (initflag .eq. 'ss'.or.ievolve.eq.'ss') then
                  call avnow(0.d0)
                  call steady(totpnow)
                  zave = zbarneq(xxxdum)
                  totpnow = denow/(zelse+zave)
                  tev=tevm(1)
                  boltzt = boltz*tev*11604.5
                  c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
                  call depress(ioff)
                  if(ioff.eq.1)return
                  if (ionow .eq. 0) then
                     popout(indlev(1,ionow)) = denow/(zelse+atomz)
                  else
                     call avnow(0.d0)
                     call steady(totpnow)
                  endif
                  nxx=0
 291              zave2 = zbarneq(xxxdum)
                  xxx=abs(zave2-zave)
                  yyy=min(zave, zave2)
                  if(xxx.gt.yyy.or.xxx.gt.0.2d0)then
                     nxx=nxx+1
                     if(nxx.gt.20)then
                        write(ioo,'(" $$$ Fail to converge")')
                        write(16,'(" Warning:Fail to converge")')
                        goto 300
                     else if(nxx.gt.10) then
                        if(zave2-zave.lt.0.d0) then
                           zave = 0.6*zave+0.4*zave2
                        else
                           zave = 0.4*zave+0.6*zave2
                        endif
                     else
                        zave = 0.5*(zave+zave2)
                     endif
                     totpnow = denow/(zelse+zave)
                     tev=tevm(1)
                     boltzt = boltz*tev*11604.5
                     c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
                     call depress(ioff)
                     if(ioff.eq.1)return
                     if (ionow .eq. 0) then
                        popout(indlev(1,ionow)) = denow/(zelse+atomz)
                     else
                        call avnow(0.d0)
                        call steady(totpnow)
                     endif
                     goto 291
                  else
                     zave=zave2
                     totpnow = denow/(zelse+zave)
                     tev=tevm(1)
                     boltzt = boltz*tev*11604.5
                     c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
                     call depress(ioff)
                     if(ioff.eq.1)return
                     if (ionow .eq. 0) then
                        popout(indlev(1,ionow)) = denow/(zelse+atomz)
                     else
                        call avnow(0.d0)
                        call steady(totpnow)
                     endif
                     goto 300
                  endif
               endif
 300           continue
               call fill

               totpnow =denow/(zelse+zbarnum(qqq))
               ratioNZ =totpnow/qqq
               do 20 i=1,numpop
                  popout(i) = popout(i)*ratioNZ
 20            continue
               if (sizenow .ne. 0.00) then
                  if (initflag .eq. 'ss'.or.ievolve.eq.'ss') then
                     call avnow(sizenow) ! ratioNZ taken care of above
                     call steady(totpnow)
                     call fill
                  endif
               endif
            endif
         endif
      endif
      


!!!!!!!! Absolutely, this part is working in the opacity simulation (grid : T and Ni input) - minsang mscho 191229 !!!!!!!!!!!!!!!!!!!!

c     HISTORY option RHO or NT only
      
      if ((iruntyp1 .eq. 'rho') .or. (iruntyp1 .eq. 'nt')) then
         
         if (initflag .eq. 'file'.and.it.eq.1) then
            dentotal = 0.00
            do 1197 i=1,numpop
               dentotal =popinit(i)+dentotal
 1197       continue

            if (iruntyp1 .eq. 'rho') then
               calcNZ = (1.-percent)*rhoAtT/atmass(iatomz)
            else if (iruntyp1 .eq. 'nt') then
               calcNZ = totpAtT/(1.+percent/(1.-percent)*amratio)
            endif

            ratioNZ = calcNZ/dentotal
            do 1196 i=1,numpop
               popout(i) =popinit(i)*ratioNZ
 1196       continue
            znow = zbarnum(qqqq)
            totpnow = qqqq
            denow = totpnow*(zelse+znow)
            zave = znow


         else ! mscho 191229 : this part working!! 

            if (iruntyp1 .eq. 'rho') then
               totpnow =(1.-percent)*rhoAtT/atmass(iatomz)
            else if (iruntyp1 .eq. 'nt') then
               totpnow = totpAtT/(1.+percent/(1.-percent)*amratio)
            endif

c....     Find initial e density by getting estimate of zbar
c....       Thomas-Fermi for lte or coronal for steady state
            denow = totpnow*iatomz/2.
            zave = iatomz/2. 
            ! mscho 191229  origin : denow = totpnow*iatomz/2.  & denow = iatomz/2
            ! mscho 191229 : I tested initial zave = iatomz/6, but difference is not meaningful. 
            ! mscho 191229 : definitely, average ionization depends on the next part  
 
            tev=tevm(1)
            call depress(ioff) ! mscho 191229 : this is 2nd line of IPD in the file 'outfileinfo' 
            if (ionow .eq. 0) then
               pop(neq) = totpnow
               znow=atomz
            endif
            call avnow(0.d0)

            if(ioff.eq.1)return
            cplgam=4.*pi*denow*r_debye**3
            if (initflag .eq. 'lte'.or.ievolve.eq.'lte') then
               call zbarTF(atomz,totpnow,tev,znow)
            else if (initflag .eq. 'ss'.or.ievolve.eq.'ss') then
               if(cplgam.lt.0.5) then
                  call zbarTF(atomz,totpnow,tev,znow)
               else
                  call coronal(totpnow,znow)
                  call fill
               endif
            endif
            denow = totpnow*(zelse+znow)
            zave = znow 
            tev=tevm(1)
            boltzt = boltz*tev*11604.5
            c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
            call depress(ioff)
            if(ioff.eq.1)return

c  if there are no levels only the bare nucleus put population in bare

            if (ionow .eq. 0) then
               popout(indlev(1,ionow)) = totpnow
            else
               if (initflag .eq. 'lte'.or.ievolve.eq.'lte') then
                  call lte(totpnow)
                  call fill
                  zave2=zbarnum(qqq) !mscho 191229 : original zbar set. Initial ne/ni = zbar 
                  denow = totpnow*(zelse+zave2)
                  if(zave2.lt.1.d-7.and.zave2.lt.zave)then
                     denow = totpnow*(zelse+zave)
                     !write(ioo,*)'test' !mscho 191229 test - no showing  
                     call depress(ioff)
                     if(ioff.eq.1)return
                     call lte(totpnow)
                     call fill
                     denow=totpnow*(zelse+zave)
                     return
                  endif
                  nxx=0

           !!!!!!! this part gives the iteration for ztot - mscho 191229 !!!!!!!!!
 2002             zave2=zbarneq(qqq) !mscho 191229 : zbar re-calculation. LSODE solver produces 'zbar'  
                  xxx=abs(zave2-zave) !mscho 191229 : compare zbarneq(qqq) with zbarTF(atomz, totpnow, tev, znow)
                  yyy=min(zave, zave2)
                  zzz=xxx/zave2
                  if(xxx.gt.yyy.or.zzz.gt.0.01d0)then !mscho 191229 : this means "If Gap is still large"
                     nxx=nxx+1 !mscho 191229 : nxx means iteration number 
                     if(nxx.gt.30)then !mscho 191229 - If # of iteration is over than 30, it gives divergence value (keep in mind LTE/opacity) 
                        write(ioo,'("  $$ Fail to converge Z1 & Z2 ",
     &                       2f10.4)')zave,zave2
                        write(ieout,'(" Warning:Fail to converge")')
                        zave=(zave+zave2)/2.
                        denow=totpnow*(zelse+zave)
                        call depress(ioff)
                        if(ioff.eq.1)return
                        call lte(totpnow)
                        call fill
                        return
                     else if(nxx.gt.20)then
                        if(zave2-zave.lt.0.d0) then
                           zave = 0.6*zave+0.4*zave2
                        else
                           zave = 0.4*zave+0.6*zave2
                        endif
                     else if(nxx.gt.10) then
                        if(zave2-zave.lt.0.d0) then
                           zave = 0.6*zave+0.4*zave2
                        else
                           zave = 0.4*zave+0.6*zave2
                        endif
                     else
                        !write(ioo,*)'mscho' !mscho-191229 this part is working 5 times (or a little bit more) in each grid condition 
                        zave = 0.5*(zave+zave2)
                     endif
                     denow=totpnow*(zelse+zave)
                     !write(ioo,*)'mscho' !mscho-191229 this part is also working 5 times in each grid condition 
                     call depress(ioff)
                     if(ioff.eq.1)return
                     call lte(totpnow)
                     call fill
                     goto 2002
                  else
                     if(zave2.gt.1.d-10)then
                        zave=zave2
                     else
                        zave=znow
                     endif
                     denow=totpnow*(zelse+zave)
                    !write(ioo,*)'mscho' !mscho-191229 test: this part is also working just 1 time!!  
                     call depress(ioff)
                     if(ioff.eq.1)return
                     call lte(totpnow)
                     call fill
                     zave2=zbarneq(qqq)
                     if(zave2.gt.0.d0)then
                        zave=zave2
                        denow=totpnow*(zelse+zave)
                     endif
                     return
                  endif
               else if (initflag .eq. 'ss'.or.ievolve.eq.'ss') then
                  call avnow(0.d0)
                  call steady(totpnow)
                  call fill
                  iter=0
 2000             znew = zbarnum(totdum)
                  iter=iter+1
                  if (abs((znew-zave)/zave) .gt. 0.01) then
                     if(iter.gt.5) then
                        zave=(zave+znew)/2.0
                     else
                        zave = znew
                     endif
                     denow = totpnow*(zelse+zave)
                     tev=tevm(1)
                     boltzt = boltz*tev*11604.5
                     c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
                     call depress(ioff)
                     if(ioff.eq.1)return
                     if (ionow .eq. 0) then
                        pop(neq) = totpnow
                     else
                        call steady(totpnow)
                     endif
                     call fill
                     denow = totpnow*(zelse+zbarnum(qqq))
                     goto 2000
                  endif
                  if (sizenow .ne. 0.00) then
                     call avnow(sizenow*ratioNZ)
                     call steady(totpnow)
                     call fill
                     denow = totpnow*(zelse+zbarnum(qqq))
                  endif
               endif
            endif
         endif
      endif
      return
      end

      ! This subroutine is the same with worknlte except IPD. 
      ! For no ipd application, all depress functions are erased. 
      subroutine worknlte_init(ioff,it,ratioNZ)
      implicit real*8 (a-h,o-z)

      include 'mainstuf'
      include 'runstuf'
      include 'timestuf'
      include 'popstuf'
      common /paramnow/ trnow,tinow,zbarnow,totpnow,sizenow,rhonow
      common /outrt/filnn
      character*10 filnn

      external steady,depress,lte

      ioff=0
      write(filnn,'("rate.",i3.3)')it

      znow = zb1
      denow = d1
      totpnow = totpAtT/(1.+percent/(1.-percent)*amratio)
      rhonow = rhoAtT
      sizenow = size1

      tev = t1
      boltzt = boltz*tev*11604.5
      c1 = 2.*(2.*pi*emass*boltzt/hp**2)**1.5

      tiev=ti1
      trev=tr1
      tr = trev*11604.5
      dilution=dlt1
      zmix= zmix1

      if (radflag .eq. 'trfile') then
         do i=1,ngrid
            barjnow(i) = barj1(i)
         enddo
         barjend = 0.00
         if (xbarj .gt. 0.00.and.trev.gt.0.d0) then
            xbytrev = xev(ngrid)/trev
            if (xbytrev .lt. 100.) then
               barjend = barjnow(ngrid)/xev(ngrid)**3
     &              *exp(xbytrev)
            endif
         endif
      endif
      

c it's important to start with the right zave
      zave=0.1d0
      do iz=2, iatomz
         if(tev.ge.eipz(iz)/2..and.tev.lt.eipz(iz-1)/2.) then
            zave=iatomz-real(iz)+1
            goto 10
         endif
      enddo
 10   continue

c     correction on 03/11/2011
c     always initialize

      tev=tevm(1)
      boltzt = boltz*tev*11604.5
      c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
      
      call depress(ioff)

c  HISTORY option NE with either RHO or NT

      if (iruntyp1 .eq. 'ne' .and.
     1   (iruntyp2 .eq. 'rho' .or. iruntyp2 .eq. 'nt')) then

         if (iruntyp2 .eq. 'nt') then
            totpnow = totpAtT/(1.+percent/(1.-percent)*amratio)
         else if (iruntyp2 .eq. 'rho' .and. elseznum .eq. 0.00
     1           .and. percent .ne. 0.00) then
            write(ioo,*)' *** Must give Z of other species'
            write(16,*)' Warning:Must give Z of other species'
            ioff=1
            return
         else
            totpnow =rhoAtT/atmass(iatomz)*(1.-percent)
         endif

         if (initflag .eq. 'file'.and.it.eq.1) then
            dentotal = 0.00
            do i=1,numpop
               dentotal = popinit(i)+dentotal
            enddo
            
            do i=1,numpop
               popout(i) = popinit(i)*totpnow/dentotal
            enddo

         else
            
            tev=tevm(1)
            boltzt = boltz*tev*11604.5
            c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
            
            call depress(ioff)
            if(ioff.eq.1)return

c     if there are no levels only the bare nucleus put population in bare
            if (ionow .eq. 0) then
               popout(indlev(1,ionow)) = totpnow
            else
               if (initflag .eq. 'lte'.or.ievolve.eq.'lte') then
                  call lte(totpnow)
                  zave = zbarneq(xxxdum)
                  tev=tevm(1)
                  boltzt = boltz*tev*11604.5
                  c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
                  call depress(ioff)
                  if(ioff.eq.1)return
                  call lte(totpnow)
               else if (initflag .eq. 'ss'.or.ievolve.eq.'ss') then
                  call avnow(0.d0)
                  call steady(totpnow)
                  zave = zbarneq(xxxdum)

                  tev=tevm(1)
                  boltzt = boltz*tev*11604.5
                  c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
                  call depress(ioff)
                  if(ioff.eq.1)return
                  call steady(totpnow)
c iteration
                  nxx=0
 281              zave2 = zbarneq(xxxdum)
                  xxx=abs(zave2-zave)
                  yyy=min(zave, zave2)
                  if(xxx.gt.yyy.or.xxx.gt.0.2d0)then
                     nxx=nxx+1
                     if(nxx.gt.20)then
                        write(ioo,'(" *** Fail to converge")')
                        write(16,'(" Warning:Fail to converge")')
                        goto 282
                     else if(nxx.gt.10) then
                        if(zave2-zave.lt.0.d0) then
                           !write(ioo,*)'mstest'  !mscho 191229 
                           zave = 0.6*zave+0.4*zave2
                        else
                           !write(ioo,*)'mstest' !mscho 191229
                           zave = 0.4*zave+0.6*zave2
                        endif
                     else
                        !write(ioo,*)'mstest' !mscho 191229
                        zave = 0.5*(zave+zave2)
                     endif
                     tev=tevm(1)
                     boltzt = boltz*tev*11604.5
                     c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
                     call depress(ioff)  
                     if(ioff.eq.1)return
                     if(ionow.eq.0)then
                        pop(neq) = totpnow !ms : pop(neq) is used in zbarneq(xxxdum) function (191227) 
                        goto 282
                     endif
                     call avnow(0.d0)
                     call steady(totpnow)
                     goto 281
                  else
                     zave=zave2
                     tev=tevm(1)
                     boltzt = boltz*tev*11604.5
                     c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
                     call depress(ioff)
                     if(ioff.eq.1)return
                     if(ionow.eq.0)then
                        pop(neq) = totpnow
                        goto 282
                     endif
                     call avnow(0.d0)
                     call steady(totpnow)
                     goto 282
                  endif
 282              continue
                  if (sizenow .ne. 0.00) then
                     call fill
                     zxxx=zbarnum(qqq)
                     ratioNZ =totpnow/qqq
                     do i=1,numpop
                        popout(i) = popout(i)*ratioNZ
                     enddo
                     call avnow(sizenow)
                     call steady(totpnow)
                  endif
               endif
               call fill
            endif
         endif
        
c.... Test the agreement between the NE and the initial NT
        
         calcNE = totpnow*(zelse+zbarnum(qqq))
         ratioNE = calcNE/denow
         write(ioo,'(" $$$  Ne(calc)/Ne(input) = ",1pe11.2)') ratioNE
         write(16,'(" $$$ Ne(calc)/Ne(input) = ",1pe11.2)') ratioNE
      endif
      


c     HISTORY option NE with neither RHO or NT
      
      if ((iruntyp1 .eq. 'ne') .and. (iruntyp2 .eq. ' ')) then
         
         if (initflag .eq. 'file'.and.it.eq.1) then
            dentotal = 0.00
            do 1199 i=1,numpop
               popout(i) = popinit(i)
               dentotal =popinit(i)+dentotal
 1199       continue
            calcNZ = denow/(zelse+zbarnum(qqq))
            if(calcNZ.gt.1.d30) calcNZ=dentotal
            ratioNZ = calcNZ/dentotal
            totpnow = calcNZ

            do 1198 i=1,numpop
               popout(i) = popinit(i)*ratioNZ
 1198       continue

         else
            
c     determine which levels exist


            totpnow = denow/(zelse+zave)
            tev=tevm(1)
            boltzt = boltz*tev*11604.5
            c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5

            call depress(ioff)
            if(ioff.eq.1)return
            if (ionow .eq. 0) then
               popout(indlev(1,ionow)) = denow/(zelse+atomz)
            else
               if (initflag .eq. 'lte'.or.ievolve.eq.'lte') then
                  call lte(totpnow)
                  zave = zbarneq(xxxdum)
                  totpnow = denow/(zelse+zave)

                  tev=tevm(1)
                  boltzt = boltz*tev*11604.5
                  c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
                  call depress(ioff)
                  if(ioff.eq.1)return
                  if (ionow .eq. 0) then
                     popout(indlev(1,ionow)) = denow/(zelse+atomz)
                  else
                     call lte(totpnow)
                  endif
                  nxx=0
 290              zave2=zbarneq(qqq)
                  xxx=abs(zave2-zave)
                  yyy=min(zave, zave2)
                  if(xxx.gt.yyy.or.xxx.gt.0.2d0)then
                     nxx=nxx+1
                     if(nxx.gt.20)then
                        write(ioo,'(" $$$  Fail to converge")')
                        write(16,'(" Warning:Fail to converge")')
                        goto 300
                     else if(nxx.gt.10) then
                        if(zave2-zave.lt.0.d0) then
                           zave = 0.6*zave+0.4*zave2
                        else
                           zave = 0.4*zave+0.6*zave2
                        endif
                     else
                        zave = 0.5*(zave+zave2)
                     endif
                     totpnow = denow/(zelse+zave)
                     tev=tevm(1)
                     boltzt = boltz*tev*11604.5
                     c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
                     call depress(ioff)
                     if(ioff.eq.1)return
                     if (ionow .eq. 0) then
                        popout(indlev(1,ionow)) = denow/(zelse+atomz)
                     else
                        call lte(totpnow)
                     endif
                     goto 290
                  else
                     zave=zave2
                     totpnow = denow/(zelse+zave)
                     tev=tevm(1)
                     boltzt = boltz*tev*11604.5
                     c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
                     call depress(ioff)
                     if(ioff.eq.1)return
                     if (ionow .eq. 0) then
                        popout(indlev(1,ionow)) = denow/(zelse+atomz)
                     else
                        call lte(totpnow)
                     endif
                     goto 300
                  endif
               else if (initflag .eq. 'ss'.or.ievolve.eq.'ss') then
                  call avnow(0.d0)
                  call steady(totpnow)
                  zave = zbarneq(xxxdum)
                  totpnow = denow/(zelse+zave)
                  tev=tevm(1)
                  boltzt = boltz*tev*11604.5
                  c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
                  call depress(ioff)
                  if(ioff.eq.1)return
                  if (ionow .eq. 0) then
                     popout(indlev(1,ionow)) = denow/(zelse+atomz)
                  else
                     call avnow(0.d0)
                     call steady(totpnow)
                  endif
                  nxx=0
 291              zave2 = zbarneq(xxxdum)
                  xxx=abs(zave2-zave)
                  yyy=min(zave, zave2)
                  if(xxx.gt.yyy.or.xxx.gt.0.2d0)then
                     nxx=nxx+1
                     if(nxx.gt.20)then
                        write(ioo,'(" $$$ Fail to converge")')
                        write(16,'(" Warning:Fail to converge")')
                        goto 300
                     else if(nxx.gt.10) then
                        if(zave2-zave.lt.0.d0) then
                           zave = 0.6*zave+0.4*zave2
                        else
                           zave = 0.4*zave+0.6*zave2
                        endif
                     else
                        zave = 0.5*(zave+zave2)
                     endif
                     totpnow = denow/(zelse+zave)
                     tev=tevm(1)
                     boltzt = boltz*tev*11604.5
                     c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
                     call depress(ioff)
                     if(ioff.eq.1)return
                     if (ionow .eq. 0) then
                        popout(indlev(1,ionow)) = denow/(zelse+atomz)
                     else
                        call avnow(0.d0)
                        call steady(totpnow)
                     endif
                     goto 291
                  else
                     zave=zave2
                     totpnow = denow/(zelse+zave)
                     tev=tevm(1)
                     boltzt = boltz*tev*11604.5
                     c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
                     call depress(ioff)
                     if(ioff.eq.1)return
                     if (ionow .eq. 0) then
                        popout(indlev(1,ionow)) = denow/(zelse+atomz)
                     else
                        call avnow(0.d0)
                        call steady(totpnow)
                     endif
                     goto 300
                  endif
               endif
 300           continue
               call fill

               totpnow =denow/(zelse+zbarnum(qqq))
               ratioNZ =totpnow/qqq
               do 20 i=1,numpop
                  popout(i) = popout(i)*ratioNZ
 20            continue
               if (sizenow .ne. 0.00) then
                  if (initflag .eq. 'ss'.or.ievolve.eq.'ss') then
                     call avnow(sizenow) ! ratioNZ taken care of above
                     call steady(totpnow)
                     call fill
                  endif
               endif
            endif
         endif
      endif
      


!!!!!!!! Absolutely, this part is working in the opacity simulation (grid : T and Ni input) - minsang mscho 191229 !!!!!!!!!!!!!!!!!!!!

c     HISTORY option RHO or NT only
      
      if ((iruntyp1 .eq. 'rho') .or. (iruntyp1 .eq. 'nt')) then
         
         if (initflag .eq. 'file'.and.it.eq.1) then
            dentotal = 0.00
            do 1197 i=1,numpop
               dentotal =popinit(i)+dentotal
 1197       continue

            if (iruntyp1 .eq. 'rho') then
               calcNZ = (1.-percent)*rhoAtT/atmass(iatomz)
            else if (iruntyp1 .eq. 'nt') then
               calcNZ = totpAtT/(1.+percent/(1.-percent)*amratio)
            endif

            ratioNZ = calcNZ/dentotal
            do 1196 i=1,numpop
               popout(i) =popinit(i)*ratioNZ
 1196       continue
            znow = zbarnum(qqqq)
            totpnow = qqqq
            denow = totpnow*(zelse+znow)
            zave = znow


         else ! mscho 191229 : this part working!! 

            if (iruntyp1 .eq. 'rho') then
               totpnow =(1.-percent)*rhoAtT/atmass(iatomz)
            else if (iruntyp1 .eq. 'nt') then
               totpnow = totpAtT/(1.+percent/(1.-percent)*amratio)
            endif

c....     Find initial e density by getting estimate of zbar
c....       Thomas-Fermi for lte or coronal for steady state
            denow = totpnow*iatomz/2.
            zave = iatomz/2. 
            ! mscho 191229  origin : denow = totpnow*iatomz/2.  & denow = iatomz/2
            ! mscho 191229 : I tested initial zave = iatomz/6, but difference is not meaningful. 
            ! mscho 191229 : definitely, average ionization depends on the next part  
 
            tev=tevm(1)
            call depress(ioff) ! mscho 191229 : this is 2nd line of IPD in the file 'outfileinfo' 
            if (ionow .eq. 0) then
               pop(neq) = totpnow
               znow=atomz
            endif
            call avnow(0.d0)

            if(ioff.eq.1)return
            cplgam=4.*pi*denow*r_debye**3
            if (initflag .eq. 'lte'.or.ievolve.eq.'lte') then
               call zbarTF(atomz,totpnow,tev,znow)
            else if (initflag .eq. 'ss'.or.ievolve.eq.'ss') then
               if(cplgam.lt.0.5) then
                  call zbarTF(atomz,totpnow,tev,znow)
               else
                  call coronal(totpnow,znow)
                  call fill
               endif
            endif
            denow = totpnow*(zelse+znow)
            zave = znow 
            tev=tevm(1)
            boltzt = boltz*tev*11604.5
            c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
            call depress(ioff)
            if(ioff.eq.1)return

c  if there are no levels only the bare nucleus put population in bare

            if (ionow .eq. 0) then
               popout(indlev(1,ionow)) = totpnow
            else
               if (initflag .eq. 'lte'.or.ievolve.eq.'lte') then
                  call lte(totpnow)
                  call fill
                  zave2=zbarnum(qqq) !mscho 191229 : original zbar set. Initial ne/ni = zbar 
                  denow = totpnow*(zelse+zave2)
                  if(zave2.lt.1.d-7.and.zave2.lt.zave)then
                     denow = totpnow*(zelse+zave)
                     !write(ioo,*)'test' !mscho 191229 test - no showing  
                     call depress(ioff)
                     if(ioff.eq.1)return
                     call lte(totpnow)
                     call fill
                     denow=totpnow*(zelse+zave)
                     return
                  endif
                  nxx=0

           !!!!!!! this part gives the iteration for ztot - mscho 191229 !!!!!!!!!
 2002             zave2=zbarneq(qqq) !mscho 191229 : zbar re-calculation. LSODE solver produces 'zbar'  
                  xxx=abs(zave2-zave) !mscho 191229 : compare zbarneq(qqq) with zbarTF(atomz, totpnow, tev, znow)
                  yyy=min(zave, zave2)
                  zzz=xxx/zave2
                  if(xxx.gt.yyy.or.zzz.gt.0.01d0)then !mscho 191229 : this means "If Gap is still large"
                     nxx=nxx+1 !mscho 191229 : nxx means iteration number 
                     if(nxx.gt.30)then !mscho 191229 - If # of iteration is over than 30, it gives divergence value (keep in mind LTE/opacity) 
                        write(ioo,'("  $$ Fail to converge Z1 & Z2 ",
     &                       2f10.4)')zave,zave2
                        write(ieout,'(" Warning:Fail to converge")')
                        zave=(zave+zave2)/2.
                        denow=totpnow*(zelse+zave)
                        call depress(ioff)
                        if(ioff.eq.1)return
                        call lte(totpnow)
                        call fill
                        return
                     else if(nxx.gt.20)then
                        if(zave2-zave.lt.0.d0) then
                           zave = 0.6*zave+0.4*zave2
                        else
                           zave = 0.4*zave+0.6*zave2
                        endif
                     else if(nxx.gt.10) then
                        if(zave2-zave.lt.0.d0) then
                           zave = 0.6*zave+0.4*zave2
                        else
                           zave = 0.4*zave+0.6*zave2
                        endif
                     else
                        !write(ioo,*)'mscho' !mscho-191229 this part is working 5 times (or a little bit more) in each grid condition 
                        zave = 0.5*(zave+zave2)
                     endif
                     denow=totpnow*(zelse+zave)
                     !write(ioo,*)'mscho' !mscho-191229 this part is also working 5 times in each grid condition 
                     call depress(ioff)
                     if(ioff.eq.1)return
                     call lte(totpnow)
                     call fill
                     goto 2002
                  else
                     if(zave2.gt.1.d-10)then
                        zave=zave2
                     else
                        zave=znow
                     endif
                     denow=totpnow*(zelse+zave)
                    !write(ioo,*)'mscho' !mscho-191229 test: this part is also working just 1 time!!  
                     call depress(ioff)
                     if(ioff.eq.1)return
                     call lte(totpnow)
                     call fill
                     zave2=zbarneq(qqq)
                     if(zave2.gt.0.d0)then
                        zave=zave2
                        denow=totpnow*(zelse+zave)
                     endif
                     return
                  endif
               else if (initflag .eq. 'ss'.or.ievolve.eq.'ss') then
                  call avnow(0.d0)
                  call steady(totpnow)
                  call fill
                  iter=0
 2000             znew = zbarnum(totdum)
                  iter=iter+1
                  if (abs((znew-zave)/zave) .gt. 0.01) then
                     if(iter.gt.5) then
                        zave=(zave+znew)/2.0
                     else
                        zave = znew
                     endif
                     denow = totpnow*(zelse+zave)
                     tev=tevm(1)
                     boltzt = boltz*tev*11604.5
                     c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
                     call depress(ioff)
                     if(ioff.eq.1)return
                     if (ionow .eq. 0) then
                        pop(neq) = totpnow
                     else
                        call steady(totpnow)
                     endif
                     call fill
                     denow = totpnow*(zelse+zbarnum(qqq))
                     goto 2000
                  endif
                  if (sizenow .ne. 0.00) then
                     call avnow(sizenow*ratioNZ)
                     call steady(totpnow)
                     call fill
                     denow = totpnow*(zelse+zbarnum(qqq))
                  endif
               endif
            endif
         endif
      endif
      return
      end
c


c
c           ******  *******  *******    ***    ******   *     *
c          *           *     *         *   *    *    *   *   *
c          *           *     *        *     *   *    *    * *
c           *****      *     ****     *******   *    *     *
c                *     *     *        *     *   *    *     *
c                *     *     *        *     *   *    *     *
c          ******      *     *******  *     *  ******      *
c
      subroutine steady(deninit)
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'timestuf'
      include 'runstuf'
      include 'popstuf'
      include 'xc_stuf'
      common /outrt/filnn
      character*10 filnn
c
c    call steady with temperature and density at time(it+1)
c
      do i=1,neq
         pop(i) = 0.00
      enddo

c     
      if(cxflag) then

         nefmax=1
         do i=1, nebins
            if(fe(i).gt.1.e-25)then
               nefmax=i
            endif
         enddo

         if(fhot.le.0.d0) then
            write(ioo,'(" *** Check fhot in -steady-: exiting")')
            write(16,'(" Warning:Check fhot in -steady- exiting")')
            return
         else
            dhot=fhot
            if(fhot.le.1.1d0)dhot = fhot*denow
            write(16,'(" CX : dhot =",1p,2e12.4)')dhot,denow
         endif
      else
         dhot=0.d0
      endif


      ntemold=ntem
      if(rtflag) then
         if(ntem.eq.0) then
            write(ioo,'(" *** Check ntem in -steady-: exiting")')
            write(16,'(" Warning:Check ntem in -steady-: exiting")')
            return
         endif
      else
         ntem=0
      endif

c
c___zero array a(i,j)
c
      do j=1,numpop
        do i=1,numpop
           a(i,j) = 0.00
        enddo
      enddo

      if(iolevel .eq. 2) then
      open(89,file=filnn)
      write(89,'(1p,e15.5,0p,f15.4,1p,2e15.5)')timer,tevm(1),dhot,denow
      do iso=1,iatomz
         write(89,'("enot",i5,1x,2f16.5)')iso, eipz(iso)
         do i=1, nlev(iso)
            ilv=indlev(i,iso)
            if(indl2p(ilv).ne.0)
     &           write(89,101)iso,i,kname(ilv),elev(ilv),glev(ilv),
     &           (noctot(ilv,ip),ip=1,nprim),mobtot(ilv)
         enddo
      enddo
      endif

      do iso=isomin,isomax
         if(cxflag) then
            call Vcaafill(iso)
         endif
      enddo


      if(radflag.ne.'off') goto 100

      do 10 itran=1, ntran
         llo=llower(itran)
         lup=lupper(itran)
         ind1=indl2p(llo)
         ind2=indl2p(lup)
!        write(*,*) 'sub- steady 1781 ** llo lup', llo, lup
         if(ind1*ind2.eq.0) goto 10
         ilo=levl(llo)
         iup=levl(lup)
         iso=isolev(llo)
         rxdn=0.d0
         rxup=0.d0

         if(ltype(itran).eq.1)then
            rxdn=ratcon(1)*aratexx(iup,ilo,iso)+dhot*Vcxx(iup,ilo,iso)
            rxup=dhot*Vcxx(ilo,iup,iso)
            do ite=1,ntem
               tev=tevm(ite)
               frc=fnem(ite)
               rxdn=rxdn+denow*frc*ratcon(5)*cxx(iup,ilo,iso)
               rxup=rxup+denow*frc*ratcon(4)*cxx(ilo,iup,iso)
               if(iolevel.eq.2) then
               write(89,900)iso,ilo,iso,iup,tev,frc,
     &              dhot*Vcxx(ilo,iup,iso),
     &              dhot*Vcxx(iup,ilo,iso),
     &              denow*frc*cxx(ilo,iup,iso),
     &              denow*frc*cxx(iup,ilo,iso),
     &              aratexx(iup,ilo,iso)
               endif
            enddo
            if(iolevel.eq.2.and.ntem.eq.0)then
               xxx=Vcxx(iup,ilo,iso)
               yyy=cxx(iup,ilo,iso)
               if(xxx.gt.2*yyy.or.xxx.lt.0.5*yyy)then
               write(89,900)iso,ilo,iso,iup,tev,xxx/yyy,
     &              dhot*Vcxx(ilo,iup,iso),
     &              denow*cxx(ilo,iup,iso),
     &              dhot*xxx,
     &              denow*yyy,
     &              aratexx(iup,ilo,iso)
               endif
            endif
         else if(ltype(itran).eq.2)then
            rxdn=dhot*dhot*Vbetaxx(ilo,iup,iso)
     &           +dhot*Valfxx(ilo,iup,iso)
            rxup=dhot*Vsxx(ilo,iup,iso)
            do ite=1,ntem
               tev=tevm(ite)
               frc=fnem(ite)
               rxdn=rxdn
     &              +denow*frc*denow*frc*ratcon(10)*betaxx(ilo,iup,iso)
     &              +denow*frc*ratcon(7)*phrec(ilo,iup,iso)
               rxup=rxup
     &              +denow*frc*ratcon(9)*sxx(ilo,iup,iso)
               if(iolevel.eq.2) then
               write(89,901)iso,ilo,iso-1,iup,tev,frc,
     &              dhot*Vsxx(ilo,iup,iso),
     &              dhot*dhot*Vbetaxx(ilo,iup,iso),
     &              denow*frc*sxx(ilo,iup,iso),
     &              denow*frc*denow*frc*betaxx(ilo,iup,iso),
     &              dhot*Valfxx(ilo,iup,iso),
     &              denow*frc*phrec(ilo,iup,iso)
               endif
            enddo
            if(iolevel.eq.2.and.ntem.eq.0)then
               xxx=Vsxx(ilo,iup,iso)
               yyy=sxx(ilo,iup,iso)
               if(xxx.gt.2*yyy.or.xxx.lt.0.5*yyy)then
                  write(89,901)iso,ilo,iso-1,iup,tev,xxx/yyy,
     &                 dhot*Vsxx(ilo,iup,iso),
     &                 dhot*dhot*Vbetaxx(ilo,iup,iso),
     &                 dhot*Valfxx(ilo,iup,iso)
               endif
            endif
         else if(ltype(itran).eq.3)then
            rxdn=dhot*Vcapexx(ilo,iup,iso)
            rxup=ratcon(11)*autoxx(ilo,iup,iso)
            do ite=1,ntem
               tev=tevm(ite)
               frc=fnem(ite)
               rxdn=rxdn
     &              +denow*frc*ratcon(12)*capexx(ilo,iup,iso)  
               if(iolevel.eq.2) then
               write(89,902)iso,ilo,iso-1,iup,tev,frc,
     &              autoxx(ilo,iup,iso),
     &              dhot*Vcapexx(ilo,iup,iso),
     &              denow*frc*capexx(ilo,iup,iso)
               endif
            enddo
            if(iolevel.eq.2.and.ntem.eq.0)then
               write(89,902)iso,ilo,iso-1,iup,tev,
     &              autoxx(ilo,iup,iso),
     &              dhot*Vcapexx(ilo,iup,iso),
     &              denow*capexx(ilo,iup,iso)
            endif
         endif
         a(ind2,ind1)=rxup
         a(ind1,ind2)=rxdn
         a(ind1,ind1)=a(ind1,ind1)-a(ind2,ind1)
         a(ind2,ind2)=a(ind2,ind2)-a(ind1,ind2)
 10   continue

      goto 200

 100  continue


      do 20 itran=1, ntran
         llo=llower(itran)
         lup=lupper(itran)
         ind1=indl2p(llo)
         ind2=indl2p(lup)
         if(ind1*ind2.eq.0) goto 20
         ilo=levl(llo)
         iup=levl(lup)
         iso=isolev(llo)
         rxdn=0.d0
         rxup=0.d0
         if(ltype(itran).eq.1)then
            rxdn=ratcon(1)*aratexx(iup,ilo,iso)+dhot*Vcxx(iup,ilo,iso)
     &                 +ratcon(2)*bradxx(iup,ilo,iso)
            rxup=dhot*Vcxx(ilo,iup,iso)
     &                 +ratcon(3)**bradxx(ilo,iup,iso)
            do ite=1,ntem
               tev=tevm(ite)
               frc=fnem(ite)
               rxdn=rxdn+denow*frc*ratcon(5)*cxx(iup,ilo,iso)
               rxup=rxup+denow*frc*ratcon(4)*cxx(ilo,iup,iso)
               if(iolevel.eq.2) then
               write(89,900)iso,ilo,iso,iup,tev,frc,
     &              dhot*Vcxx(ilo,iup,iso),
     &              dhot*Vcxx(iup,ilo,iso),
     &              denow*frc*cxx(ilo,iup,iso),
     &              denow*frc*cxx(iup,ilo,iso),
     &              aratexx(iup,ilo,iso),
     &              bradxx(iup,ilo,iso),bradxx(ilo,iup,iso)
               endif
            enddo
         else if(ltype(itran).eq.2)then
            rxdn=dhot*dhot*Vbetaxx(ilo,iup,iso)
     &           +dhot*Valfxx(ilo,iup,iso)
     &           +dhot*Vbstmxx(ilo,iup,iso)
            rxup=dhot*Vsxx(ilo,iup,iso)
     &           +ratcon(6)*phion(ilo,iup,iso)
            do ite=1,ntem
               tev=tevm(ite)
               frc=fnem(ite)
               rxdn=rxdn
     &              +denow*frc*denow*frc*ratcon(10)*betaxx(ilo,iup,iso)
     &              +denow*frc*ratcon(7)*phrec(ilo,iup,iso)
               rxup=rxup+denow*frc*ratcon(9)*sxx(ilo,iup,iso)
               if(iolevel.eq.2) then
               write(89,901)iso,ilo,iso-1,iup,tev,frc,
     &              dhot*Vsxx(ilo,iup,iso),
     &              dhot*dhot*Vbetaxx(ilo,iup,iso),
     &              denow*frc*sxx(ilo,iup,iso),
     &              denow*frc*denow*frc*betaxx(ilo,iup,iso),
     &              dhot*Vbstmxx(ilo,iup,iso),
     &              phion(ilo,iup,iso),
     &              denow*frc*phrec(ilo,iup,iso)
               endif
            enddo
         else if(ltype(itran).eq.3)then
            rxdn=dhot*Vcapexx(ilo,iup,iso)
            rxup=ratcon(11)*autoxx(ilo,iup,iso)
            do ite=1,ntem
               tev=tevm(ite)
               frc=fnem(ite)
               rxdn=rxdn
     &              +denow*frc*ratcon(12)*capexx(ilo,iup,iso)  
               if(iolevel.eq.2) then
               write(89,902)iso,ilo,iso-1,iup,tev,frc,
     &              dhot*Vcapexx(ilo,iup,iso),
     &              ratcon(11)*autoxx(ilo,iup,iso),
     &              denow*frc*capexx(ilo,iup,iso)
               endif
            enddo
         endif
         a(ind2,ind1)=rxup
         a(ind1,ind2)=rxdn
         a(ind1,ind1)=a(ind1,ind1)-a(ind2,ind1)
         a(ind2,ind2)=a(ind2,ind2)-a(ind1,ind2)
 20   continue
c
 200  continue
      close(89)

c
c   invert the matrix set up
c
      do 400 i=1,neq
         b(i) = 0.00
 400  continue


c....  substitute the last row for the number conservation equation
c      changing from neq to 1 for low-T high-N cases to work. 


      do i=1, neq
         if(a(i,i).eq.0.d0) then
           write(ioo,'(" Warning:Zero diagonal element at ",2i5,1x,a8)')
     &           i, indp2l(i),kname(indp2l(i))
            write(16,'(" Warning:Zero diagonal element at ",2i5,1x,a8)')
     &           i, indp2l(i),kname(indp2l(i))
         endif
      enddo

c ... level index for the ground state of dominant ions
      izmpr=iatomz-int(zave)
      ilv=indlev(1,izmpr)
      ipv=indl2p(ilv)

      if(ipv.eq.0)then
         do i=izmpr-1,1,-1
            ipv=indl2p(indlev(1,i))
            if(ipv.ne.0)goto 402
         enddo
      endif

 402  azmax=0.d0
      azmin=1.d30
      do 401 i=1,neq
         azmin=min(azmin,abs(a(i,i)))
         azmax=max(azmax,abs(a(i,i)))
         a(ipv,i)=1.d0
 401  continue
      b(ipv)=deninit
c


c    solve the matrix by call into solve
c
      call solve(neq)
      popmax = 1.00d-70
      do i=1,neq
         pop(i) = x(i)
         popmax = max(popmax,pop(i))
      enddo
      do i=1,neq
         if (pop(i) .lt. popmax*1.00d-70) then
            pop(i) = popmax*1.00d-70
         endif
      enddo


      ntem=ntemold

 101  format('elev', 2i4,1x,a5,1x,f14.4,1x,1p,e16.8,0p,1x,10i3,1x,i4)
 900  format('e',2(i3,i5),1p,10e10.2)
 901  format('i',2(i3,i5),1p,12e10.2)
 902  format('a',2(i3,i5),1p,10e10.2)
 909  format('a(i,j) for i=',a8,'j=',2(2x,a8))
 910  format(1p,10e10.2)

 899  format(i3,1x,1p,20e9.2)

      return
      end


c mscho 200504: m-start: the subroutine workVnlte is added for fe.test4 (1st subroutine) 
c**********************************************************
! mscho : I copied only History option NE with either RHO or NT 
! mscho : There are two options more (History option NE with neither RHO nor NT
! mscho :                         and History option only RHO or NT) 
! mscho : After test, two options should be put! (200504) 

! mscho part1

! mscho part2

c it's important to start with the right zave

c determine the iso-sequence ranges
c      isomin=1
c      isomax=iatomz

c  HISTORY option NE with either RHO or NT


! mscho part3

c     if there are no levels only the bare nucleus put population in bare

! mscho part5 

c.... Test the agreement between the NE and the initial NT


c mscho 200504: m-finish for the subroutine workVnlte


c mscho 200504: m-start for the subroutine Vsteady (2nd subroutine) 
c This subroutine is inside the subroutine workVnlte  

! mscho part2-1
c

c mscho 200504: m-finish : The subroutine Vsteady (2nd) 






c**********************************************************
c
c          ******   *     *  ******   *******
c           *    *  **    *   *    *     *
c           *    *  * *   *   *    *     *
c           *    *  *  *  *   *    *     *
c           *    *  *   * *   *    *     *
c           *    *  *    **   *    *     *
c          ******   *     *  ******      *
c
      subroutine dndt(neqnow,timeODE,poptime,ff)
c
c*****************************
c*                           *
c*  d e c l a r a t i o n s  *
c*                           *
c*****************************
c

      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'timestuf'
      include 'popstuf'
      include 'runstuf'
      include 'xc_stuf'
c
      dimension ff(numpp),poptime(numpp)

      common /outrt/filnn
      character*10 filnn

      data timeLAST/-1.00/
c
c*******************
c*  begin routine  *
c*******************
c
c___find brackets for the temperature and the density

c      if (timeODE .eq. timeLAST) then
c         go to 171
c      endif

      timeLAST = timeODE
      timer = timeODE
c     

      tev = (t0*(tnext-timer)+t1*(timer-tlast))/(tnext-tlast)
      boltzt = boltz*tev*11604.5
      c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
      
      trev = (tr0*(tnext-timer)+tr1*(timer-tlast))/(tnext-tlast)
      if(trev.le.trmin) trev=0.d0
      tr = trev*11604.5

      zmix = (zmix0*(tnext-timer)+zmix1*(timer-tlast))/(tnext-tlast)
      zelse = zmix*(percent/(1.-percent))*amratio


c... determine the Ne from the input Ne, or pops if have only Nt or Rho

      if (iruntyp1 .eq. 'ne') then
         dinow = (totp0*(tnext-timer)+totp1*(timer-tlast))/
     &        (tnext-tlast)
         denow = exp((log(d0)*(tnext-timer)+log(d1)*(timer-tlast))/
     1        (tnext-tlast))
         
      else if (iruntyp1 .eq. 'nt') then
         dinow=exp((log(totp0)*(tnext-timer)+log(totp1)*
     1        (timer-tlast))/(tnext-tlast))
         denow = exp((log(totp0)*(tnext-timer)+log(totp1)*
     1        (timer-tlast))/(tnext-tlast))*(1.-percent)
     2        *(zelse+zbarneq(qqdum)) ! note zelse=zmix*percent/(1-percent)
      else if (iruntyp1 .eq. 'rho') then
         dinow = exp((log(rh0)*(tnext-timer)+log(rh1)*(timer-tlast))
     1        /(tnext-tlast))*(1.-percent)/atmass(iatomz)
c     already zelse has amratio multiplied
         denow =exp((log(rh0)*(tnext-timer)+log(rh1)*(timer-tlast))
     1        /(tnext-tlast))*(1.-percent)/atmass(iatomz)*
     2        (zelse+zbarneq(qqdum))
      endif

      do ii=1, ntem
         fnem(ii)= exp((log(fnem0(ii))*(tnext-timer)
     1        +log(fnem1(ii))*(timer-tlast))/(tnext-tlast))
         tevm(ii)=(tevm0(ii)*(tnext-timer)+tevm1(ii)*(timer-tlast))
     1        /(tnext-tlast)
      enddo

      if(cxflag) then
         sum=0.d0
         do i=1,nebins
            fe(i) = exp((log(fe0(i))*(tnext-timer)
     1           +log(fe1(i))*(timer-tlast))/(tnext-tlast))
            sum = sum + sqrt(energy(i))*fe(i)*denergy(i)
         enddo
         do i=1, nebins
           fe(i)=fe(i)/sum
         enddo
         nefmax=1
         do i=1, nebins
            if(fe(i).gt.1.e-25)then
               nefmax=i
            endif
         enddo
         if(fhot.gt.0.d0) then
            fhot= exp((log(fhot0)*(tnext-timer)
     1           +  log(fhot1)*(timer-tlast))/(tnext-tlast))
            dhot=fhot
            if(fhot.le.1.1d0)dhot = fhot*denow
            write(ioo,'(" CX : dhot =",1p,2e12.4)')dhot,denow
            write(16,'(" CX : dhot =",1p,2e12.4)')dhot,denow
         else
            write(ioo,'(" *** Check fhot in -dndt-: exiting")')
            write(16,'(" Warning:Check fhot in -dndt-: exiting")')
            return
         endif
      else
         dhot=0.d0
      endif


      ntemold=ntem
      if(rtflag) then
         if(ntem.eq.0) then
            write(ioo,'(" *** Check ntem in -dndt-: exiting")')
            write(16,'(" Warning:Check ntem in -dndt-: exiting")')
            return
         endif
      else
         ntem=0
      endif


c
c___zero array a(i,j)
c
      do j=1,numpop
         do i=1,numpop
            a(i,j) = 0.00
         enddo
      enddo

      if(iolevel .eq. 2) then
      open(89,file=filnn)
      write(89,'(1p,e15.5,0p,f15.4,1p,2e15.5)')timer,tevm(1),dhot,denow
      do iso=1,iatomz
         write(89,'("enot",i5,1x,2f16.5)')iso, eipz(iso)
         do i=1, nlev(iso)
            ilv=indlev(i,iso)
            if(indl2p(ilv).ne.0)
     &           write(89,101)iso,i,kname(ilv),elev(ilv),glev(ilv),
     &           (noctot(ilv,ip),ip=1,nprim),mobtot(ilv)
         enddo
      enddo
      endif

      do iso=isomin,isomax
         if(cxflag) then
            call Vcaafill(iso)
         endif
      enddo

      if(radflag.ne.'off') goto 100

      do 10 itran=1, ntran
         llo=llower(itran)
         lup=lupper(itran)
         ind1=indl2p(llo)
         ind2=indl2p(lup)
         if(ind1*ind2.eq.0) goto 10
         ilo=levl(llo)
         iup=levl(lup)
         iso=isolev(llo)
         rxdn=0.d0
         rxup=0.d0
         if(ltype(itran).eq.1)then
            rxdn=ratcon(1)*aratexx(iup,ilo,iso)
     &           +dhot*Vcxx(iup,ilo,iso)
            rxup=dhot*Vcxx(ilo,iup,iso)
            do ite=1,ntem
               tev=tevm(ite)
               frc=fnem(ite)
               rxdn=rxdn+denow*frc*ratcon(5)*cxx(iup,ilo,iso)
               rxup=rxup+denow*frc*ratcon(4)*cxx(ilo,iup,iso)
               if(iolevel.eq.2) then
               write(89,900)iso,ilo,iso,iup,tev,frc,
     &              dhot*Vcxx(ilo,iup,iso),
     &              dhot*Vcxx(iup,ilo,iso),
     &              denow*frc*cxx(ilo,iup,iso),
     &              denow*frc*cxx(iup,ilo,iso),
     &              aratexx(iup,ilo,iso)
               endif
            enddo
            if(iolevel.eq.2.and.ntem.eq.0)then
               write(89,900)iso,ilo,iso,iup,tev,frc,
     &              dhot*Vcxx(ilo,iup,iso),
     &              dhot*Vcxx(iup,ilo,iso),
     &              aratexx(iup,ilo,iso)
            endif
         else if(ltype(itran).eq.2)then
            iso=isolev(llo)
            rxdn=dhot*dhot*Vbetaxx(ilo,iup,iso)
     &           +dhot*Valfxx(ilo,iup,iso)
            rxup=dhot*Vsxx(ilo,iup,iso)
            do ite=1,ntem
               tev=tevm(ite)
               frc=fnem(ite)
               rxdn=rxdn
     &              +denow*frc*denow*frc*ratcon(10)*betaxx(ilo,iup,iso)
     &              +denow*frc*ratcon(7)*phrec(ilo,iup,iso)
               rxup=rxup
     &              +denow*frc*ratcon(9)*sxx(ilo,iup,iso)
               if(iolevel.eq.2) then
               write(89,901)iso,ilo,iso-1,iup,tev,frc,
     &              dhot*Vsxx(ilo,iup,iso),
     &              dhot*dhot*Vbetaxx(ilo,iup,iso),
     &              denow*frc*sxx(ilo,iup,iso),
     &              denow*frc*denow*frc*betaxx(ilo,iup,iso),
     &              dhot*Valfxx(ilo,iup,iso),
     &              +denow*frc*phrec(ilo,iup,iso)
               endif
            enddo
            if(iolevel.eq.2.and.ntem.eq.0)then
               write(89,901)iso,ilo,iso-1,iup,tev,frc,
     &              dhot*Vsxx(ilo,iup,iso),
     &              dhot*dhot*Vbetaxx(ilo,iup,iso),
     &              dhot*Valfxx(ilo,iup,iso)
            endif
         else if(ltype(itran).eq.3)then
            rxdn=dhot*Vcapexx(ilo,iup,iso)
            rxup=ratcon(11)*autoxx(ilo,iup,iso)
            do ite=1,ntem
               tev=tevm(ite)
               frc=fnem(ite)
               rxdn=rxdn
     &              +denow*frc*ratcon(12)*capexx(ilo,iup,iso)
               if(iolevel.eq.2) then
               write(89,902)iso,ilo,iso-1,iup,tev,frc,
     &              dhot*Vcapexx(ilo,iup,iso),
     &              autoxx(ilo,iup,iso),
     &              denow*frc*capexx(ilo,iup,iso)
               endif
            enddo
            if(iolevel.eq.2.and.ntem.eq.0)then
               write(89,902)iso,ilo,iso-1,iup,tev,frc,
     &              dhot*Vcapexx(ilo,iup,iso),
     &              autoxx(ilo,iup,iso)
            endif
         endif
         a(ind2,ind1)=rxup
         a(ind1,ind2)=rxdn
         a(ind1,ind1)=a(ind1,ind1)-a(ind2,ind1)
         a(ind2,ind2)=a(ind2,ind2)-a(ind1,ind2)
 10   continue

      goto 200

 100  continue

      if (radflag .eq. 'trfile') then
         do 1 i=1,ngrid
            barjnow(i)=(barj0(i)*(tnext-timer)+barj1(i)*(timer-tlast))
     &           /(tnext-tlast)
 1       continue
         barjend = 0.00
         if (xbarj .gt. 0.00.and.trev.gt.0.d0) then
            xbytrev = xev(ngrid)/trev
            if (xbytrev .lt. 100.) then
               barjend = barjnow(ngrid)/xev(ngrid)**3*exp(xbytrev)
            endif
         endif
      endif

      do 20 itran=1, ntran
         llo=llower(itran)
         lup=lupper(itran)
         ind1=indl2p(llo)
         ind2=indl2p(lup)
         ilo=levl(llo)
         iup=levl(lup)
         iso=isolev(llo)
         if(ind1*ind2.eq.0) goto 20
         rxdn=0.d0
         rxup=0.d0
         if(ltype(itran).eq.1)then
            atrad=bradxx(ilo,iup,iso)
            rxup=dhot*Vcxx(ilo,iup,iso)
     &                 +ratcon(3)*atrad
            rxdn=ratcon(1)*aratexx(iup,ilo,iso)+dhot*Vcxx(iup,ilo,iso)
     &                 +ratcon(2)*atrad*glev(llo)/glev(lup)
            do ite=1,ntem
               tev=tevm(ite)
               frc=fnem(ite)
               rxdn=rxdn+denow*frc*ratcon(5)*cxx(iup,ilo,iso)
               rxup=rxup+denow*frc*ratcon(4)*cxx(ilo,iup,iso)
               if(iolevel.eq.2) then
               write(89,900)iso,ilo,iso,iup,tev,frc,
     &              dhot*Vcxx(ilo,iup,iso),
     &              dhot*Vcxx(iup,ilo,iso),
     &              denow*frc*cxx(ilo,iup,iso),
     &              denow*frc*cxx(iup,ilo,iso),
     &              aratexx(iup,ilo,iso),
     &              bradxx(iup,ilo,iso),bradxx(ilo,iup,iso)
               endif
           enddo
            if(iolevel.eq.2.and.ntem.eq.0)then
              write(89,900)iso,ilo,iso,iup,tev,frc,
     &             dhot*Vcxx(ilo,iup,iso),
     &             dhot*Vcxx(iup,ilo,iso),
     &             aratexx(iup,ilo,iso),
     &             bradxx(iup,ilo,iso),bradxx(ilo,iup,iso)
           endif
         else if(ltype(itran).eq.2)then
            rxdn=dhot*dhot*Vbetaxx(ilo,iup,iso)
     &           +dhot*Valfxx(ilo,iup,iso)
     &           +dhot*Vbstmxx(ilo,iup,iso)
            rxup=dhot*Vsxx(ilo,iup,iso)
     &           +ratcon(6)*phion(ilo,iup,iso)
            do ite=1,ntem
               tev=tevm(ite)
               frc=fnem(ite)
               rxdn=rxdn
     &              +denow*frc*denow*frc*ratcon(10)*betaxx(ilo,iup,iso)
     &              +denow*frc*ratcon(7)*phrec(ilo,iup,iso)
               rxup=rxup+denow*frc*ratcon(9)*sxx(ilo,iup,iso)
               if(iolevel.eq.2) then
               write(89,901)iso,ilo,iso-1,iup,tev,frc,
     &              dhot*Vsxx(ilo,iup,iso),
     &              dhot*dhot*Vbetaxx(ilo,iup,iso),
     &              denow*frc*sxx(ilo,iup,iso),
     &              denow*frc*denow*frc*betaxx(ilo,iup,iso),
     &              dhot*Vbstmxx(ilo,iup,iso),
     &              phion(ilo,iup,iso),   
     &              denow*frc*phrec(ilo,iup,iso)
               endif
            enddo
            if(iolevel.eq.2.and.ntem.eq.0)then
               write(89,901)iso,ilo,iso-1,iup,tev,frc,
     &              dhot*Vsxx(ilo,iup,iso),
     &              dhot*dhot*Vbetaxx(ilo,iup,iso),
     &              dhot*Vbstmxx(ilo,iup,iso),
     &              phion(ilo,iup,iso)   
            endif
         else if(ltype(itran).eq.3)then
            rxdn=dhot*Vcapexx(ilo,iup,iso)
            rxup=ratcon(11)*autoxx(ilo,iup,iso)
            do ite=1,ntem
               tev=tevm(ite)
               frc=fnem(ite)
               rxdn=rxdn
     &              +denow*frc*ratcon(12)*capexx(ilo,iup,iso)  
               if(iolevel.eq.2) then
               write(89,902)iso,ilo,iso-1,iup,tev,frc,
     &              dhot*Vcapexx(ilo,iup,iso),
     &              autoxx(ilo,iup,iso),
     &              denow*frc*capexx(ilo,iup,iso)
               endif
            enddo
            if(iolevel.eq.2.and.ntem.eq.0)then
               write(89,902)iso,ilo,iso-1,iup,tev,frc,
     &              dhot*Vcapexx(ilo,iup,iso),
     &              autoxx(ilo,iup,iso)
            endif
         endif
         a(ind2,ind1)=rxup
         a(ind1,ind2)=rxdn
         a(ind1,ind1)=a(ind1,ind1)-a(ind2,ind1)
         a(ind2,ind2)=a(ind2,ind2)-a(ind1,ind2)
 20   continue
c
 200  continue
      close(89)

c___equation a*poptime = dndt returns the solution
c

  171 continue

      do i = 1,neqnow
         ff(i) = 0.00
         do j=1,neqnow
            ff(i) = ff(i)+a(i,j)*poptime(j)
         enddo
      enddo
      ntem=ntemold

 101  format('elev', 2i4,1x,a5,1x,f14.4,1x,1p,e16.8,0p,1x,10i3,1x,i4)
 900  format('e',2(i3,i5),1p,10e10.2)
 901  format('i',2(i3,i5),1p,10e10.2)
 902  format('a',2(i3,i5),1p,10e10.2)
 903  format(i3,i5,1p,10(1x,e8.1),/,10(3x,10(1x,e8.1),/))

      return 
      end

c**********************************************************
      subroutine derv(neq,timer,poptime,ml,mu,pw,nrowpw)

      implicit real*8 (a-h,o-z)
      include 'timestuf'
      include 'popstuf'

      dimension poptime(neq),pw(nrowpw,neq)

c.... calculate the partial derivatives of the rates

      do 1 i=1,neq
        jl = max(1,i-ml)
        ju = min(neq,i+mu)
        do 2 j=jl,ju
           pw(i-j+mu+1,j) = a(i,j)
    2   continue
    1 continue

      return
      end

      subroutine jac(neq,timer,poptime,ml,mu,pw,nrowpw)

      implicit real*8 (a-h,o-z)
      include 'timestuf'
      include 'popstuf'

      dimension poptime(neq),pw(nrowpw,neq)

c.... calculate the partial derivatives of the rates

      do 1 i=1,neq
         do 2 j=1, neq
            pw(i,j) = a(i,j)
 2       continue
 1    continue

      return
      end

      subroutine tdnlte(ioff,it,ratioNZ)
      implicit real*8 (a-h,o-z)

      include 'mainstuf'
      include 'runstuf'
      include 'timestuf'
      include 'popstuf'
      include 'xc_stuf'


      common /popwriteR/ bpop(numpp,NMTout),biso(0:100,NMTout),deninit
      common /paramnow/ trnow,tinow,zbarnow,totpnow,sizenow,rhonow

      parameter(lrw = 22+9*numpp+numpp**2)
      parameter(liw = 20 + numpp)         
      common /odework/ rwork(lrw),rtol(numpp),iwork(liw)

      dimension popold(numpp)
      external lsode, depress,lte,dndt,jac,initial,Vndt
      logical initdone, around, nogood

      common /outrt/filnn
      character*10 filnn

      logical rzlflag, ipdflag


      data error,atol,itol,mf,itask/0.001,0.0001,3,21,4/

c      error = .001
c      atol = 0.001
c      itol = 3
c      mf = 21                   ! for full matrix
c      itask = 4                 ! don't overshoot in time

      rwork(1) = tnext          ! NOT to overshoot


      ioff=0
      write(filnn,'("rate.",i3.3)')it

      sizenow = size0 ! up to April 2010 
      sizenow = size1 

      zave = zbarnum(dddd)
      denow = d0
      tev = t0
      boltzt = boltz*tev*11604.5
      if (iruntyp1 .eq. 'ne') then
         if (iruntyp2 .eq. ' ') then
            if(zave+zelse.le.1.d-30) then
               totpnow = 0.01*d0
            else
               totpnow = d0/(zave+zelse)
            endif
         else if (iruntyp2 .eq. 'nt') then
            totpnow =totp0/(1.+amratio*percent/(1.-percent))
         else if (iruntyp2 .eq. 'rho') then
            totpnow=rh0/atmass(iatomz)*(1.-percent)
         endif
      else if (iruntyp1 .eq. 'nt') then
         totpnow = totp0/(1.+amratio*percent/(1.-percent))
      else if (iruntyp1 .eq. 'rho') then
         totpnow = rh0/atmass(iatomz)*(1.-percent)
      endif
      ! mscho explain: totpnow is ion density 


c determine the iso-sequence ranges
c there is a danger where zave is very different 
c when hot electrons are involved.
c      isomin=max(iatomz-int(zave)-20,1)
c      isomax=min(iatomz-int(zave)+20,iatomz)

c mscho's test (200504)
c mscho 200506 : m-start : zeldaflag on 
      ipdflag=.false.
      rzlflag=.false. 
      if(zeldaflag)then
         if(initflag.eq.'ss'.or.it.gt.2)then
            ipdflag=.false.
         else
            rzlflag=.false.
         endif
      endif
      !write(*,*)'2571 line scmn.f:zeldaflag, ipdflag, rzlflag, initflag'
      !write(*,*) zeldaflag, ipdflag, rzlflag, initflag 
      !mscho result: zeldaflag=T, ipd=T(it<=2),F(it>3), rzl=F,T, init=file  

c   mscho 200506 : modified under 5 lines 
c   mscho - if you want to go back, then just erase 'if loop' 
      if(ipdflag)then
         tev=tevm0(1)
         boltzt = boltz*tev*11604.5
         c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
         call avnow(sizenow*ratioNZ)
         call depress(ioff)
         if(ioff.eq.1)return
      endif 
c mscho 200506 : m-finish 
      
      !write(*,*)'2586 line in scmn.f: ionow'
      !write(*,*) ionow 
      ! mscho ionow = 11 for this step 

      if (ionow .eq. 0) then
         do i=1, numpop
            popout(i)=0.d0
         enddo
         popout(indlev(1,ionow)) = deninit
         goto 202
      endif

      irestart = 0
      do 196  i=1,iatomz
         if (ire(i)+iraut(i) .ne. irelast(i)) then
            irestart = 1
            go to 195
         endif
 196  continue
 195  continue

      do 194 i=1,iatomz
         irelast(i) = ire(i)+iraut(i)
 194  continue

      if(irestart .eq.1) then
         do 201 i=1,neq
            rtol(i) = error
 201     continue
         call initial
      endif
      
      iopt = 0
      istate = 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c mscho 200506 : m-start : add zeldaflag again (This is key part!!)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      cxflag=.false. ! mscho 200511 add for test
      if(cxflag)then
         sum=0.d0
         do i=1,nebins
            if(fe0(i).gt.1.d-25)then
               sum=sum+((fe1(i)-fe0(i))/fe0(i))**2
            else if(fe1(i).gt.1.d-25)then
               sum=sum+((fe1(i)-fe0(i))/fe1(i))**2
            endif
         enddo
         sum=sqrt(sum)/nebins
         if(sum.lt.200.d0)then
            call fVndt
            call lsode(Vndt,neq,pop,tlast,tnext,itol,rtol,atol
     1           ,itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)
         else
            call lsode(dndt,neq,pop,tlast,tnext,itol,rtol,atol
     1           ,itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)
         endif
      else if(zeldaflag)then
         tfirst=times(1)
         tzelda=tnext-tlast
         dtfly=(tnext-tlast)/(isfly*1.d0)
         do i=1,neq
            popold(i)=pop(i)
         enddo
         ntimer=0
    
         do i=1, nebins !mscho 200506 initialize for the zelda 
            fe(i)=fe0(i)
         enddo 
         
         if(mod(it,2).eq.0)then
            izf=0
         else
            izf=1
         endif

         if(rzlflag) call runzelda(it,izf)

         do i=1,nebins !mscho 200506 finalize for the zelda 
            fe(i)=fe1(i)
         enddo
         tevfout(it)=tev
         fevfout(it)=denow
         tevm1(1)=tev
         fnem1(1)=1.d0 
         
         t1=tev
         call lsode(dndt,neq,pop,tlast,tnext,itol,rtol,atol
     1        ,itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)


      else 
         call lsode(dndt,neq,pop,tlast,tnext,itol,rtol,atol
     1        ,itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)
      endif
      cxflag=.true. ! mscho 200511 add (for test)

      if (iouttime(it) .gt. 0) then
         write(16,210) tnext,tevm(1),denow
 210     format(///,' LSODE Diagnostics for time',1pe10.2,
     1        /, 1x,'Tev Ne ',2e10.2)
         write(16,211) istate,it
         write(16,212) (iwork(k),k=11,13)
 211     format(5x,' istate & it',i3,1x,i3)
 212     format(5x,' nst,nfe,nje',5i6)
      endif
      
      if (istate .ne. 2) then
         write(ioo,'(/," *** Lsode error. See INFO File.",/)')
         write(16,'(/," Warning:Lsode error. See INFO File.",/)')
         ipos = ichar('a')
         nameinit = "init"//nametime(1:4)
 2130    inquire(file=nameinit,exist=around)
         if (around) then
            nameinit = "init"//nametime(1:4)//char(ipos)
            ipos = ipos+1
            if (ipos .eq. 27) then
               write(ioo,*) " *** Trouble writing NAMEINIT file"
               write(16,*) " Warning:Trouble writing NAMEINIT file"
               go to 2131
            endif
            go to 2130
         endif
         open(unit=33,file=nameinit,status="new")
         do 214 i=1,numpop
            write(33,831) kname(indpop(i)),bpop(i,it-1)
 214     continue
         close(33)
         write(16,213) nameinit,it-1,times(it-1),teout(it-1)
 213     format(1x," NameInit ",a8," for step ",i2,
     1        " at time ",1pe12.2," & Te ",e12.2,
     2        " written",/)
         
 2131    continue
         if (it .gt. iolow) then
            call writeout(iolow,it-1)
         endif
         close(21)
         close(16)
         ioff=1
         return
      endif
      
      call fill
      
 202  continue

      if (iruntyp1 .eq. 'ne') then
         if (iruntyp2 .eq. ' ') then
            totpnow = eneout(it)/(zelse+zbarnum(qqq))
            denow = eneout(it)
         else if (iruntyp2 .eq. 'nt') then
            totpnow =totpout(it)/(1.+amratio*percent/(1.-percent))
            denow = eneout(it)
         else if (iruntyp2 .eq. 'rho') then
            totpnow=rhoout(it)/atmass(iatomz)*(1.-percent)
            denow = eneout(it)
         endif
      else if (iruntyp1 .eq. 'nt') then
         totpnow = totpout(it)/(1.+amratio*percent/(1.-percent))
         denow = totpnow*(zelse+zbarnum(xxxx))
      else if (iruntyp1 .eq. 'rho') then
         totpnow = rhoout(it)/atmass(iatomz)*(1.-percent)
         denow = totpnow*(zelse+zbarnum(xxxx))
      endif
c     deninit doesn't changet after the first step.
c     popout normalized to deninit always. 
      ratioNZ = totpnow/deninit

 831  format(1x,a8,1x,1pe12.4,9e12.4)

      return
      end


      subroutine tdtcal(ioff,it,ratioNZ)
      implicit real*8 (a-h,o-z)

      include 'mainstuf'
      include 'runstuf'
      include 'timestuf'
      include 'popstuf'
      include 'xc_stuf'


      common /popwriteR/ bpop(numpp,NMTout),biso(0:100,NMTout),deninit
      common /paramnow/ trnow,tinow,zbarnow,totpnow,sizenow,rhonow

      parameter(lrw = 22+9*numpp+numpp**2)
      parameter(liw = 20 + numpp)         
      common /odework/ rwork(lrw),rtol(numpp),iwork(liw)

      dimension popold(numpp)
      dimension  ff(numpp)

      external lsode, depress,lte,tedt,jac,initial
      logical initdone, around, nogood, stuck

      common /outrt/filnn
      character*10 filnn

      data error,atol,itol,mf,itask/0.001,0.0001,3,21,4/

c     necessary for self-consistent Te calculations
      common /absorb/abb(100),abf(100),aff(100),att(100)
      common /vhydro/dtev1,dtev2,zbar1,dedt1,dedt2,teint1,teint2,ehydro
      common /ihydro/ihydflag

      logical rzlflag, ipdflag 
      logical ipdflag_in, ipdflag_out
      logical cxflag_in 

c      error = .001
c      atol = 0.001
c      itol = 3
c      mf = 21                   ! for full matrix
c      itask = 4                 ! don't overshoot in time


      open(unit=79,file='ipd_monitor',status='old',position='append')
        write(79,'("ver1. it =",I2)') it
        do i = 1, 600
          write(79,*) kname(i),popout(i)
        enddo
      close(79)



      rwork(1) = tnext          ! NOT to overshoot

      ioff=0
      write(filnn,'("rate.",i3.3)')it
      rzlflag =.true. 
      sizenow = size0 ! up to April 2010 
      sizenow = size1 

      zave = zbarnum(dddd)
      denow = d0
      tev = t0
      boltzt = boltz*tev*11604.5
      if (iruntyp1 .eq. 'ne') then
         if (iruntyp2 .eq. ' ') then
            if(zave+zelse.le.1.d-30) then
               totpnow = 0.01*d0
            else
               totpnow = d0/(zave+zelse)
            endif
         else if (iruntyp2 .eq. 'nt') then
            totpnow =totp0/(1.+amratio*percent/(1.-percent))
         else if (iruntyp2 .eq. 'rho') then
            totpnow=rh0/atmass(iatomz)*(1.-percent)
         endif
      else if (iruntyp1 .eq. 'nt') then
         totpnow = totp0/(1.+amratio*percent/(1.-percent))
      else if (iruntyp1 .eq. 'rho') then
         totpnow = rh0/atmass(iatomz)*(1.-percent)
      endif


c     check
      do ii=1, ntem
         tevm1(ii)=max(tevbal,0.001)
      enddo

      jbmax=1
      do i=1,ngrid
         if(barj1(i).gt.barj1(jbmax))then
            jbmax=i
         endif
      enddo

      ipdflag_out=.false.
      ipdflag_in=.false.
      rzlflag=.true. 
      if(zeldaflag)then
        if(initflag.eq.'ss'.or.it.gt.0)then
            ipdflag_in=.false.
        else
            rzlflag=.false. ! no run-zelda in the 1st time step 
        endif
      endif 

c ---  start the second time loop ---
      
      itzmax=2  ! itzmax original = 20 (mscho change 210603)
                ! itzmax = 2 for fast calculation (test): it should be back to 20
                ! for convergence 
      dtx0=(tnext-tlast)/itzmax
      tlast0=tlast
      tnext0=tlast
      atmaxx=atmass(iatomz)/amass
      zave=zbarnum(qqdum)

      dtx=dtx0
      itz=0
      stuck = .false.

      ! This is the starting point of time loop in "tdtcal"
      ! Te converges in this loop 

      if(ipdflag_out)then
        tev=tevm0(1)
        boltzt = boltz*tev*11604.5
        c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
        call depress(ioff)
        if(ioff.eq.1)return

        if (tiflag .eq. 'off')then
           tiev=tevm0(1)
        else
           call AtTime(tnext)
           if (tiflag .eq. 'file') then
              if (tiAtT .le. 0.00) then
                 tiev = teout(it)
              else
                 tiev = tiAtT
              endif
           else if (tiflag .eq. 'fixed') then
              tiev = tifixed
           else if (tiflag .eq. 'off') then
              tiev = teout(it)
           else if (tiflag .eq. 'ti/te') then
              tiev = teout(it)*tibyte
           endif
           tiev=ti0
           ti1=tiev
         endif
      endif 


 100  itz=itz+1 
      dt0=dtx
      write(*,*)'tnext0 > tnext?:escape',tnext0, tnext
      if(tnext0.gt.tnext) goto 200

c     correction on 03/11/2011
c     avnow should be called after depress
c     avnow : calculates the escape factor reductions of the a-values 

      if(ipdflag_in)then
        tev = tevm0(1)
        boltzt = boltz*tev*11604.5
        c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
        call depress(ioff)

        if(ioff.eq.1)return

        if (tiflag .eq. 'off')then
           tiev=tevm0(1)
        else
           call AtTime(tnext)
           if (tiflag .eq. 'file') then
              if (tiAtT .le. 0.00) then
                 tiev = teout(it)
              else
                 tiev = tiAtT
              endif
           else if (tiflag .eq. 'fixed') then
              tiev = tifixed
           else if (tiflag .eq. 'off') then
              tiev = teout(it)
           else if (tiflag .eq. 'ti/te') then
              tiev = teout(it)*tibyte
           endif
           tiev=ti0
           ti1=tiev
        endif
        call avnow(sizenow*ratioNZ)

      else 
        tev = tevm0(1)
        boltzt = boltz*tev*11604.5
        c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5

        if(ioff.eq.1)return

        if (tiflag .eq. 'off')then
           tiev=tevm0(1)
        else
           call AtTime(tnext)
           if (tiflag .eq. 'file') then
              if (tiAtT .le. 0.00) then
                 tiev = teout(it)
              else
                 tiev = tiAtT
              endif
           else if (tiflag .eq. 'fixed') then
              tiev = tifixed
           else if (tiflag .eq. 'off') then
              tiev = teout(it)
           else if (tiflag .eq. 'ti/te') then
              tiev = teout(it)*tibyte
           endif
           tiev=ti0
           ti1=tiev
        endif
        call avnow(sizenow*ratioNZ)
      endif



      if (ionow .eq. 0) then
         do i=1, numpop
            popout(i)=0.d0
         enddo
         popout(indlev(1,ionow)) = deninit
         goto 202
      endif

      irestart = 0
      do 196  i=1,iatomz
         if (ire(i)+iraut(i) .ne. irelast(i)) then
            irestart = 1
            go to 195
         endif
 196  continue
 195  continue

      do 194 i=1,iatomz
         irelast(i) = ire(i)+iraut(i)
 194  continue

      if(irestart .eq.1) then
         do 201 i=1,neq
            rtol(i) = error
 201     continue
         call initial
      endif


c     compute the level energy 
      call caleipd(S_E,P_F)


      if(irestart.eq.1) then    ! IPD imposed

c     compute temperature due to instantaneous redistribution by IPD

         totin=0.d0
         sumpp=0.d0
         do i = 1,neq
            ilev = indp2l(i)
            totin = totin + efgrdp(ilev)*pop(i) ! total internal energy
            sumpp = sumpp + pop(i)
         enddo
         teinti= totin/sumpp
         zbar = zbarneq(xxxdum) ! after IPD

         Etot1 = 1.500000d0*(tiev+zbar1*tevbal) +teint1 ! E(kinetic)+E(int)before IPD

c     change Apr 2014
c         tevbal = (Etot1-teinti-1.5*tiev)/1.500000d0/zbar ! new Te by IPD 

         if(tiflag.eq.'off')then
            tevbal=(Etot1-teinti)/1.500000d0/(zbar+1.d0) ! new Te by IPD 
         else
            tevbal = (Etot1-teinti-1.5*tiev)/1.500000d0/zbar ! new Te by IPD 
         endif

         if(tevbal.le.0.d0)then ! hyun's change  to make zero tevbal changed'
            write(*,'("Te adjusted by IPD", 1p,10e10.2)')
     &           zbar1,zbar,Etot1,teint1,totin/sumpp
            tevbal=0.001
            dt0=dtx*0.2
         else if(tevbal.le.0.01)then
            dt0=dtx*0.2
         endif

         teint1 = teinti   ! internal energy changed
         if(abs(zbar-zbar1).gt.0.1) dt0=dtx*0.2

      endif

c     run to tnext0


      tnext0=tlast0+dt0


      if(tnext0.gt.tnext)then
         dt0=tnext-tlast0
         tnext0=tnext
         if(dt0.lt.dtx/100.)then
            goto 200
         else if(dt0.gt.dtx)then
            dt0=dtx
         endif
         tnext0=tlast0+dt0
      endif

      icalev = 0
      tephabs = tephab0
      ttlast = tlast0
      
      iopt = 0
      istate = 1
      rwork(1) = tnext0

c    set up temperature from the last tevbal
c     check
      do ii=1, ntem
         tevm0(ii)=max(tevbal,0.001)
         tevm1(ii)=max(tevbal,0.001)
      enddo
      if(tiflag .eq. 'off') tiev=tevm0(1)


      zbar1 = zbarneq(xxxdum)
      tevb1= tevbal
      hsize1=sizenow 


      csound0=sqrt(zbar1*1.666667*tevb1/atmaxx)*9.82428d5/1e12 !   in cm/ps
      csound1=0.d0             ! no motion; zero sound speed  
      if(ihydflag.gt.0)then
         csound1=csound0
      endif

      if(zeldaflag)then
          cxflag_in = .false.
          cxflag = .true.
      endif 

      
      if(cxflag_in)then
         sum=0.d0
         do i=1,nebins
            if(fe0(i).gt.1.d-25)then
               sum=sum+((fe1(i)-fe0(i))/fe0(i))**2
            else if(fe1(i).gt.1.d-25)then
               sum=sum+((fe1(i)-fe0(i))/fe1(i))**2
            endif
         enddo
         sum=sqrt(sum)/nebins
         if(sum.lt.200.d0)then
            call fVndt
c            call lsode(Vtdt,neq,pop,tlast0,tnext0,itol,rtol,atol
c     1           ,itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)
         else
            call lsode(tedt,neq,pop,tlast0,tnext0,itol,rtol,atol
     1           ,itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)
         endif
      else if(zeldaflag)then 
         tfirst=times(1)
         tzelda=tnext-tlast
         dtfly=(tnext-tlast)/(isfly*1.d0)
         do i=1,neq
            popold(i)=pop(i)
         enddo
         ntimer=0
         
         do i=1, nebins
            fe(i)=fe0(i)
         enddo

         if(mod(it,2).eq.0)then
            izf=0
         else
            izf=1
         endif 
         
         if(rzlflag) call runzelda(it,izf)
          
         do i=1,nebins
            fe(i)=fe1(i)
         enddo
         !tevfout(it)=tev
         !fevfout(it)=denow
         !tevm1(1)=tev
         !fnem1(1)=1.d0

         t1=tev
         call lsode(tedt,neq,pop,tlast0,tnext0,itol,rtol,atol
     1        ,itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)
        
         write(*,*)'tedt finish in the zeldaflag loop' 
      else
         call lsode(tedt,neq,pop,tlast0,tnext0,itol,rtol,atol
     1        ,itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)
      endif

c --- error message ----------------------------------------------

      if (istate .ne. 2) then
         write(ioo,'(/," *** Lsode error. See INFO File.",/)')
         write(16,'(/," Warning:Lsode error. See INFO File.",/)')
         ipos = ichar('a')
         nameinit = "init"//nametime(1:4)
 2130    inquire(file=nameinit,exist=around)
         if (around) then
            nameinit = "init"//nametime(1:4)//char(ipos)
            ipos = ipos+1
            if (ipos .eq. 27) then
               write(ioo,*) " *** Trouble writing NAMEINIT file"
               write(16,*) " Warning:Trouble writing NAMEINIT file"
               go to 2131
            endif
            go to 2130
         endif
         open(unit=33,file=nameinit,status="new")
         do 214 i=1,numpop
            write(33,831) kname(indpop(i)),bpop(i,it-1)
 214     continue
         close(33)
         write(16,213) nameinit,it-1,times(it-1),teout(it-1)
 213     format(1x," NameInit ",a8," for step ",i2,
     1        " at time ",1pe12.2," & Te ",e12.2,
     2        " written",/)
         
 2131    continue
         if (it .gt. iolow) then
            call writeout(iolow,it-1)
         endif
         close(21)
         close(16)
         ioff=1
         return
      endif
c --------------------------------   error message  done -----------

c      call tedt(neq,tnext0,pop,ff)

      call fill


c     write absorption rates for the 1st time step
      if(itz.eq.1)then
         iwabs=1 
         write(28,'("time = ",1p,e10.2)')tnext0
      else
         iwabs=0
      endif

      if(ngrid.ne.0)call abscrx(ratioNZ,iwabs)

      write(29,'(1p,2e10.3,3e9.2,0p,f9.1,1p,10e9.2)')tnext0,tevbal,
c     &     denow,totpnow,sizenow,xev(jbmax),barj1(jbmax),
     &     denow,totpnow,sizenow,xev(jbmax),barjnow(jbmax),
     &     att(jbmax),abb(jbmax),abf(jbmax),aff(jbmax)

      write(31,'(1p,5e11.3)')tnext0,tevbal,denow,totpnow,sizenow
      write(31,'(1p,10e11.3)')(att(i),i=1,ngrid)
      !write(*,*)'***mscho tevbal test:',tevbal
      
c --------------- energy balance ----------------------

      totin=0.d0
      sumpp=0.d0
      do i = 1,neq
         ilev = indp2l(i)
         totin = totin + efgrdp(ilev)*pop(i) ! total internal energy
         sumpp = sumpp + pop(i)
      enddo
      
      teint2=totin/sumpp        ! Eint per ion
      dedt2=ephabs/sumpp        ! dE/dt per ion
      zave = zbarneq(xxxx)

      zbar2 = zave

c ... two temperatures using dE/dt(tlast0) and dE/dt(tnext0)   
c     use them for initial and next temperatures

      dtev1= (dedt1*dt0 ! increment by radiation by dedt1 at tlast0
     &     - (teint2-teint1)    ! increment by internal energy
     &     - csound1*dt0/hsize1*(tevb1*zbar1+tiev) ) ! by linear expansion

      dtev2= (dedt2*dt0 ! increment by radiation by dedt2 at tnext0
     &     - (teint2-teint1)    ! increment by internal energy
     &     - csound1*dt0/hsize1*(tevb1*zbar1+tiev) ) ! by linear expansion

      !write(*,*)'dedt1,dedt2,teint2,teint1,csound1,dt0,hsize1,tiev:'
      !write(*,*) dedt1,dedt2,teint2,teint1,csound1,dt0,hsize1,tiev
      !hsize1 = 1.d-0
      dtev=( (dedt1+dedt2)/2.*dt0 ! increment by radiation
     &     - (teint2-teint1)    ! increment by internal energy
     &     - csound1*dt0/hsize1*(tevb1*zbar1+tiev) ) ! by linear expansion

      !write(*,*)'tevb1,zbar1,dtev,zbar2:',tevb1,' ',zbar1,' ',dtev,' ',zbar2
      ! result: tevb1 = 2.58e-2, zbar1=1e-69, dtev = NAN, zbar2=2e-314
      if(tiflag.eq.'off')then
         tevbal=(tevb1*1.5*(1.d0+zbar1)+dtev)/(1.5*(1.d0+zbar2))
         tiev=tevbal
      else
         tevbal=(tevb1*1.5*zbar1+dtev)/(1.5*zbar2)          
      endif
      !write(*,*)'tiflag and tevbal:',tevbal,' ',tiflag
      ! result: tiflag = off / tevbal = NAN 

      do ii=1, ntem
         tevm0(ii)=max(tevb1,0.001)
         tevm1(ii)=max(tevbal,0.001)
      enddo
      !write(*,*)'simple test tevbal',tevbal

c     constant should be  3/2 * (r-1) = 3/2 * 2/3 = 1.
c ... determine if the time step is appropriate
c     proposed dt1 to give  min(0.1 * Te, 1 eV) increment

      dtx=dtx0
      dtevid = min(0.01*tevb1,1.d0)  ! less than 1 eV
      dt1= (dtevid+(teint2-teint1))/
     &     ((dedt1+dedt2)/2.
     &     -csound1/hsize1*tevb1*(1.+zbar1))
      if(dt1.gt.dtx.and.tevbal.gt.1.) then
         if(dt1.gt.2.*dtx)then
            dtx=dtx*2.
         else
            dtx=dt1  
         endif
      endif

c ... change in Rho, Ni due to hydro motion
      if(ihydflag.gt.0)then
         sizenow=hsize1 + csound1*dt0 ! change in size
         topold=totpnow
         totpnow=totpnow*hsize1/sizenow ! change in Te
         rhonow =rhonow*hsize1/sizenow ! change in rho
      endif

c ... new energy balance calculations due to hydro motion
      tephabs=tephab0+(dedt1+dedt2)/2.*dt0*sumpp*sizenow ! in eV/cm2
      ehydro=ehydro+csound1*dt0*tevb1*(1.+zbar1)*topold ! in eV/cm2
      ekin=sizenow*(tiev*totpnow+tevbal*denow)*1.5d0 ! in eV/cm2
      eint=totpnow*sizenow*teint2 ! in eV/cm2
      eresid = tephabs-(ekin+eint+ehydro) ! in eV/cm2

      tephab0=tephabs


      write(*,'("t=",1p,e11.4,": Cal_Te=",1p,e9.2,",Ne=",1p,e9.2,
     &     ",Tot_E=",1p,e9.2,",Int_E=",1p,e9.2,"[eV/cm2]",
     &     ",dE(rad)/dt=",1p,e9.2,"[eV/cm2/s]",
     &     ",C(sound)=",1p,e9.2,"[cm/ps]",1p,4e9.2)')
     &     tnext0,tevbal,denow,tephab0,eint, ! in eV/cm2
     &     (dedt1+dedt2)/2.*sumpp*sizenow, ! in eV/cm2/s
     &     csound0 !   in cm/s

      write(27,'("t=",1p,e11.4,": Cal_Te=",1p,e9.2,",Ne=",1p,e9.2,
     &     ",Tot_E=",1p,e9.2,",Int_E=",1p,e9.2,"[eV/cm2]",
     &     ",dE(rad)/dt=",1p,e9.2,"[eV/cm2/s]",
     &     ",C(sound)=",1p,e9.2,"[cm/ps]",1p,4e9.2)')
     &     tnext0,tevbal,denow,tephab0,eint, ! in eV/cm2
     &     (dedt1+dedt2)/2.*sumpp*sizenow, ! in eV/cm2/s
     &     csound0 !   in cm/s

      if(tevbal.le.0.d0)then
         print *, 'Warning: negative temperature ', tevbal
         tevbal=0.001
      endif

      if(abs(dtev2-dtev1).gt.0.2*abs(dtev))dtx=dt0/2.

c ... Updating the new temperature

      do ii=1, ntem
         tevm0(ii)=max(tevm1(ii),0.001)
      enddo


      dedt1 = dedt2
      teint1=teint2
      zbar1 = zave

      tlast0=tnext0
      if(tnext0.lt.tnext) goto 100

 200  continue

c check
      t1 = tevbal
      tev = tevbal

      do i=1, ntem
         tevm1(i)=max(tevbal,0.001)
         tevm(i)=max(tevbal,0.001)
      enddo


      if(tiflag.eq.'off')then
         ti1 = tevbal
      else
         ti1 = tiev
      endif

c----- end  of  time loop ----------


      teout(it)=tevbal
      tevmout(1,it)=tevbal
      fnemout(1,it)=1.

      if (iouttime(it) .gt. 0) then
         write(16,210) tnext,tevm(1),denow
 210     format(///,' LSODE Diagnostics for time',1pe10.2,
     1        /, 1x,'Tev Ne ',2e10.2)
         write(16,211) istate,it
         write(16,212) (iwork(k),k=11,13)
 211     format(5x,' istate & it',i3,1x,i3)
 212     format(5x,' nst,nfe,nje',5i6)
      endif
      
 202  continue



c     with hydro motion, totpnow, rhonow, sizenow changes
      if(ihydflag.gt.0)then
         totp1=totpnow
         rh1=rhonow
         size1=sizenow
      else
         if (iruntyp1 .eq. 'ne') then
            if (iruntyp2 .eq. ' ') then
               totpnow = eneout(it)/(zelse+zbarnum(qqq))
               denow = eneout(it)
            else if (iruntyp2 .eq. 'nt') then
               totpnow =totpout(it)/(1.+amratio*percent/(1.-percent))
               denow = eneout(it)
            else if (iruntyp2 .eq. 'rho') then
               totpnow=rhoout(it)/atmass(iatomz)*(1.-percent)
               denow = eneout(it)
            endif
         else if (iruntyp1 .eq. 'nt') then
            totpnow = totpout(it)/(1.+amratio*percent/(1.-percent))
            denow = totpnow*(zelse+zbarnum(xxxx))
         else if (iruntyp1 .eq. 'rho') then
            totpnow = rhoout(it)/atmass(iatomz)*(1.-percent)
            denow = totpnow*(zelse+zbarnum(xxxx))
         endif
      endif

c     deninit doesn't changet after the first step.
c     pop and popout arrays are normalized to deninit always. 
c     
      ratioNZ = totpnow/deninit

 831  format(1x,a8,1x,1pe12.4,9e12.4)

      return
      end

      function efevfe (x)
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'xc_stuf'
      sum=0.d0
      esum=0.d0
      do i=1,nebins
         sum=sum+fe(i)*sqrt(energy(i))*denergy(i)
         esum=esum+fe(i)*sqrt(energy(i))*denergy(i)
     &        *energy(i)
      enddo
      efevfe=2./3.*esum/sum  ! averaged hot electron temperature

      return
      end

      integer function indxc(eii,fxxx)
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'xc_stuf'

      indxc=0  ! it's the necessary condition
      !fxxx=0.
      !write(*,*) '3380 line in scmn.f',nebins
      if(eii.le.0.d0.or.eii.gt.energy(nebins))return
      if(inftyp.eq.0)then
         indxc = INT(eii/denergy(2)) + 1
         !fxxx=fe(indxc)         ! correction 9/24/2007
         !write(*,*)'3382scmn:inftyp,indxc,fe(i)',inftyp,indxc,fe(indxc)
         return
      else 
         if(inftyp.eq.1)then
            ifix=int(log10(eii))
         else if(inftyp.eq.2)then
            ifix=int(log(eii))
         endif
         i0=indfx(ifix-1)
         i1=indfx(ifix)
         i2=indfx(ifix+1)
         write(*,*)'line 3394 in scmn.f: i0,i1,i2)',i0,i1,i2
         if(eii.ge.energy(i1))then
            do i=i1,i2-1
               if(eii.ge.energy(i).and.eii.lt.energy(i+1))then
                  indxc=i
                  fxxx=(fe(i+1)*(eii-energy(i))+
     &                 fe(i)*(energy(i+1)-eii))/
     &                 (energy(i+1)-energy(i))
                  return
               endif
            enddo
         else
            do i=i0,i1-1
               if(eii.gt.energy(i).and.eii.le.energy(i+1))then
                  indxc=i
                  fxxx=(fe(i+1)*(eii-energy(i))+
     &                 fe(i)*(energy(i+1)-eii))/
     &                 (energy(i+1)-energy(i))
                  return
               endif
            enddo
         endif
      endif
      return
      end
c


      subroutine fVndt
c
c*****************************
c*                           *
c*  d e c l a r a t i o n s  *
c*                           *
c*****************************
c

      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'timestuf'
      include 'popstuf'
      include 'runstuf'
      include 'xc_stuf'
c
c
c*******************
c*  begin routine  *
c*******************
c
c ************ CALCULATIONS ** AT TLAST ***************
      nefmax=1
      do i=1, nebins
         fe(i)=fe0(i)
         if(fe(i).gt.1.e-25)then
            nefmax=i
         endif
      enddo

      do iso=isomin,isomax
         call Vcaafill(iso)
      enddo

      do 10 itran=1, ntran
         llo=llower(itran)
         lup=lupper(itran)
         ind1=indl2p(llo)
         ind2=indl2p(lup)
         if(ind1*ind2.eq.0) goto 10
         n1=n1+1
         ilo=levl(llo)
         iup=levl(lup)
         iso=isolev(llo)
         if(ltype(itran).eq.1)then
            Xcxx(iup,ilo,iso)=log(max(Vcxx(iup,ilo,iso),1.d-70))
            Xcxx(ilo,iup,iso)=log(max(Vcxx(ilo,iup,iso),1.d-70))
         else if(ltype(itran).eq.2)then
            Xbetaxx(ilo,iup,iso)=log(max(Vbetaxx(ilo,iup,iso),1.d-70))
            Xalfxx(ilo,iup,iso)=log(max(Valfxx(ilo,iup,iso),1.d-70))
            Xsxx(ilo,iup,iso)=log(max(Vsxx(ilo,iup,iso),1.d-70))
            if(radflag.ne.'off')
     &      Xbstmxx(ilo,iup,iso)=log(max(Vbstmxx(ilo,iup,iso),1.d-70))
            Xbcontxx(ilo,iup,iso)=Vbcontxx(ilo,iup,iso) !mscho 200504 add (because of compile error) 
            ! mscho 200504 Xbcontxx is in the zchkmn.f in CT27 (so I follows that one)  

         else if(ltype(itran).eq.3)then
            Xcapexx(ilo,iup,iso)=log(max(Vcapexx(ilo,iup,iso),1.d-70))
         endif
 10   continue

c ************ CALCULATIONS ** AT TNEXT ***************
      nefmax=1
      do i=1, nebins
         fe(i)=fe1(i)
         if(fe(i).gt.1.e-25)then
            nefmax=i
         endif
      enddo
      do iso=isomin,isomax
         call Vcaafill(iso)
      enddo

      do 20 itran=1, ntran
         llo=llower(itran)
         lup=lupper(itran)
         ind1=indl2p(llo)
         ind2=indl2p(lup)
         if(ind1*ind2.eq.0) goto 20
         n1=n1+1
         ilo=levl(llo)
         iup=levl(lup)
         iso=isolev(llo)
         if(ltype(itran).eq.1)then
            Xcyy(iup,ilo,iso)=log(max(Vcxx(iup,ilo,iso),1.d-70))
            Xcyy(ilo,iup,iso)=log(max(Vcxx(ilo,iup,iso),1.d-70))

         else if(ltype(itran).eq.2)then
            Xbetayy(ilo,iup,iso)=log(max(Vbetaxx(ilo,iup,iso),1.d-70))
            Xalfyy(ilo,iup,iso)=log(max(Valfxx(ilo,iup,iso),1.d-70))
            Xsyy(ilo,iup,iso)=log(max(Vsxx(ilo,iup,iso),1.d-70))
            if(radflag.ne.'off')
     &      Xbstmyy(ilo,iup,iso)=log(max(Vbstmxx(ilo,iup,iso),1.d-70))
         else if(ltype(itran).eq.3)then
            Xcapeyy(ilo,iup,iso)=log(max(Vcapexx(ilo,iup,iso),1.d-70))
         endif

 20   continue

c******************************************************

      return 

      end

      subroutine Vndt(neqnow,timeODE,poptime,ff)
c
c*****************************
c*                           *
c*  d e c l a r a t i o n s  *
c*                           *
c*****************************
c

      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'timestuf'
      include 'popstuf'
      include 'runstuf'
      include 'xc_stuf'
c
      dimension ff(numpp),poptime(numpp)

      common /outrt/filnn
      character*10 filnn

c
c*******************
c*  begin routine  *
c*******************
c
c___find brackets for the temperature and the density

c      if (timeODE .eq. timeLAST) then
c         go to 171
c      endif

      timeLAST = timeODE
      timer = timeODE
      dt0=(timer-tlast)/(tnext-tlast)
      dt1=(tnext-timer)/(tnext-tlast)
c     

      tev = (t0*dt1+t1*dt0)
      boltzt = boltz*tev*11604.5
      c1=2.*(2.*pi*emass*boltzt/hp**2)**1.5
      
      trev = (tr0*dt1+tr1*dt0)
      if(trev.le.trmin) trev=0.d0
      tr = trev*11604.5

      zmix = (zmix0*(tnext-timer)+zmix1*(timer-tlast))/(tnext-tlast)
      zelse = zmix*(percent/(1.-percent))*amratio      
      
c... determine the Ne from the input Ne, or pops if have only Nt or Rho
      
      if (iruntyp1 .eq. 'ne') then
         denow = exp(log(d0)*dt1+log(d1)*dt0)
      else if (iruntyp1 .eq. 'nt') then
         denow = exp(log(totp0)*dt1+log(totp1)*dt0)
     1        *(1.-percent)
     2        *(zelse+zbarneq(qqdum))
      else if (iruntyp1 .eq. 'rho') then
         denow =exp(log(rh0)*dt1+log(rh1)*dt0)
     1        *(1.-percent)/atmass(iatomz)*
     2        (zelse+zbarneq(qqdum))
      endif
      
      do ii=1, ntem
         fnem(ii)= exp(log(fnem0(ii))*dt1+log(fnem1(ii))*dt0)
         tevm(ii)=(tevm0(ii)*dt1+tevm1(ii)*dt0)
      enddo

      if(cxflag) then
         sum=0.d0
         do i=1,nebins
            fe(i) = exp( log(fe0(i))*dt1+log(fe1(i))*dt0 )
            sum = sum + sqrt(energy(i))*fe(i)*denergy(i)
         enddo
         do i=1, nebins
           fe(i)=fe(i)/sum
         enddo
         nefmax=1
         do i=1, nebins
            if(fe(i).gt.1.e-25)then
               nefmax=i
            endif
         enddo
         if(fhot.gt.0.d0) then
            fhot= exp(log(fhot0)*dt1+ log(fhot1)*dt0)
            dhot=fhot
            if(fhot.le.1.1d0)dhot = fhot*denow
            write(16,'(" CX : dhot =",1p,2e12.4)')dhot,denow
         else
            write(16,'(" Warning:Check fhot in -dndt-: exiting")')
            return
         endif
      else
         dhot=0.d0
      endif


      ntemold=ntem
      if(rtflag) then
         if(ntem.eq.0) then
            write(ioo,'(" *** Check ntem in -dndt-: exiting")')
            write(16,'(" Warning:Check ntem in -dndt-: exiting")')
            return
         endif
      else
         ntem=0
      endif


c
c___zero array a(i,j)
c
      do j=1,numpop
         do i=1,numpop
            a(i,j) = 0.00
         enddo
      enddo

      if(iolevel.eq.2) then
      open(89,file=filnn)
      write(89,'(1p,e15.5,0p,f15.4,1p,2e15.5)')timer,tevm(1),dhot,denow
      do iso=1,iatomz
         write(89,'("enot",i5,1x,2f16.5)')iso, eipz(iso)
         do i=1, nlev(iso)
            ilv=indlev(i,iso)
            if(indl2p(ilv).ne.0)
     &      write(89,101)iso,i,kname(ilv),elev(ilv),glev(ilv),
     &           (noctot(ilv,ip),ip=1,nprim),mobtot(ilv)
         enddo
      enddo
      endif

      if(radflag.ne.'off') goto 100

      do 10 itran=1, ntran
         llo=llower(itran)
         lup=lupper(itran)
         ind1=indl2p(llo)
         ind2=indl2p(lup)
         if(ind1*ind2.eq.0) goto 10
         ilo=levl(llo)
         iup=levl(lup)
         iso=isolev(llo)
         rxdn=0.d0
         rxup=0.d0
         if(ltype(itran).eq.1)then
            rxdn=ratcon(1)*aratexx(iup,ilo,iso)+
     &           dhot*exp(Xcxx(iup,ilo,iso)*dt1+Xcyy(iup,ilo,iso)*dt0)
            rxup=dhot*exp(Xcxx(ilo,iup,iso)*dt1+Xcyy(ilo,iup,iso)*dt0)
            do ite=1,ntem
               tev=tevm(ite)
               frc=fnem(ite)
               rxdn=rxdn+denow*frc*ratcon(5)*cxx(iup,ilo,iso)
               rxup=rxup+denow*frc*ratcon(4)*cxx(ilo,iup,iso)
               if(iolevel.eq.2) then
               write(89,900)iso,ilo,iso,iup,tev,frc,
     &              dhot*Xcxx(ilo,iup,iso),
     &              dhot*Xcxx(iup,ilo,iso),
     &              denow*frc*cxx(ilo,iup,iso),
     &              denow*frc*cxx(iup,ilo,iso),
     &              aratexx(iup,ilo,iso)
               endif
            enddo
            if(iolevel.eq.2.and.ntem.eq.0)then
               write(89,900)iso,ilo,iso,iup,tev,frc,
     &              dhot*Xcxx(ilo,iup,iso),
     &              dhot*Xcxx(iup,ilo,iso),
     &              aratexx(iup,ilo,iso)
            endif
         else if(ltype(itran).eq.2)then
            iso=isolev(llo)
            rxdn=dhot*dhot*
     &           exp(Xbetaxx(ilo,iup,iso)*dt1+Xbetayy(ilo,iup,iso)*dt0)
     &           +dhot*
     &           exp(Xalfxx(ilo,iup,iso)*dt1+Xalfyy(ilo,iup,iso)*dt0)
            rxup=dhot*exp(Xsxx(ilo,iup,iso)*dt1+Xsyy(ilo,iup,iso)*dt0)
            do ite=1,ntem
               tev=tevm(ite)
               frc=fnem(ite)
               rxdn=rxdn
     &              +denow*frc*denow*frc*ratcon(10)*betaxx(ilo,iup,iso)
     &              +denow*frc*alfxx(ilo,iup,iso)
               rxup=rxup
     &              +denow*frc*ratcon(9)*sxx(ilo,iup,iso)
               if(iolevel.eq.2) then
               write(89,901)iso,ilo,iso-1,iup,tev,frc,
     &              dhot*Xsxx(ilo,iup,iso),
     &              dhot*dhot*Xbetaxx(ilo,iup,iso),
     &              denow*frc*sxx(ilo,iup,iso),
     &              denow*frc*denow*frc*betaxx(ilo,iup,iso),
     &              dhot*Xalfxx(ilo,iup,iso),
     &              denow*frc*alfxx(ilo,iup,iso)
               endif
            enddo
            if(iolevel.eq.2.and.ntem.eq.0)then
               write(89,901)iso,ilo,iso-1,iup,tev,frc,
     &              dhot*Xsxx(ilo,iup,iso),
     &              dhot*dhot*Xbetaxx(ilo,iup,iso),
     &              dhot*Xalfxx(ilo,iup,iso)
            endif
         else if(ltype(itran).eq.3)then
            rxdn=dhot*
     &           exp(Xcapexx(ilo,iup,iso)*dt1+Xcapeyy(ilo,iup,iso)*dt0)
            rxup=ratcon(11)*autoxx(ilo,iup,iso)
            do ite=1,ntem
               tev=tevm(ite)
               frc=fnem(ite)
               rxdn=rxdn
     &              +denow*frc*ratcon(12)*capexx(ilo,iup,iso)  
               if(iolevel.eq.2) then
               write(89,902)iso,ilo,iso-1,iup,tev,frc,
     &              dhot*Xcapexx(ilo,iup,iso),
     &              autoxx(ilo,iup,iso),
     &              denow*frc*capexx(ilo,iup,iso)
               endif
            enddo
            if(iolevel.eq.2.and.ntem.eq.0)then
               write(89,902)iso,ilo,iso-1,iup,tev,frc,
     &              dhot*Xcapexx(ilo,iup,iso),
     &              autoxx(ilo,iup,iso)
            endif
         endif
         a(ind2,ind1)=rxup
         a(ind1,ind2)=rxdn
         a(ind1,ind1)=a(ind1,ind1)-a(ind2,ind1)
         a(ind2,ind2)=a(ind2,ind2)-a(ind1,ind2)
 10   continue

      goto 200

 100  continue

      if (radflag .eq. 'trfile') then
         do 1 i=1,ngrid
            barjnow(i)=(barj0(i)*dt1+barj1(i)*dt0)
 1       continue
         barjend = 0.00
         if (xbarj .gt. 0.00.and.trev.gt.0.d0) then
            xbytrev = xev(ngrid)/trev
            if (xbytrev .lt. 100.) then
               barjend = barjnow(ngrid)/xev(ngrid)**3*exp(xbytrev)
            endif
         endif
      endif

      do 20 itran=1, ntran
         llo=llower(itran)
         lup=lupper(itran)
         ind1=indl2p(llo)
         ind2=indl2p(lup)
         if(ind1*ind2.eq.0) goto 20
         ilo=levl(llo)
         iup=levl(lup)
         iso=isolev(llo)
         rxdn=0.d0
         rxup=0.d0
         if(ltype(itran).eq.1)then
            rxdn=ratcon(1)*aratexx(iup,ilo,iso)+
     &           dhot*exp(Xcxx(iup,ilo,iso)*dt1+Xcyy(iup,ilo,iso)*dt0)
     &                 +ratcon(2)*bradxx(iup,ilo,iso)
            rxup=dhot*exp(Xcxx(ilo,iup,iso)*dt1+Xcyy(ilo,iup,iso)*dt0)
     &                 +ratcon(3)*bradxx(ilo,iup,iso)
            do ite=1,ntem
               tev=tevm(ite)
               frc=fnem(ite)
               rxdn=rxdn+denow*frc*ratcon(5)*cxx(iup,ilo,iso)
               rxup=rxup+denow*frc*ratcon(4)*cxx(ilo,iup,iso)
               if(iolevel.eq.2) then
               write(89,900)iso,ilo,iso,iup,tev,frc,
     &              dhot*Xcxx(ilo,iup,iso),
     &              dhot*Xcxx(iup,ilo,iso),
     &              denow*frc*cxx(ilo,iup,iso),
     &              denow*frc*cxx(iup,ilo,iso),
     &              aratexx(iup,ilo,iso),
     &              bradxx(iup,ilo,iso),bradxx(ilo,iup,iso)
               endif
           enddo
           if(iolevel.eq.2.and.ntem.eq.0)then
              write(89,900)iso,ilo,iso,iup,tev,frc,
     &             dhot*Xcxx(ilo,iup,iso),
     &             dhot*Xcxx(iup,ilo,iso),
     &             aratexx(iup,ilo,iso),
     &             bradxx(iup,ilo,iso),bradxx(ilo,iup,iso)
           endif
         else if(ltype(itran).eq.2)then
            rxdn=dhot*dhot*
     &           exp(Xbetaxx(ilo,iup,iso)*dt1+Xbetayy(ilo,iup,iso)*dt0)
     &           +dhot*
     &           exp(Xalfxx(ilo,iup,iso)*dt1+Xalfyy(ilo,iup,iso)*dt0)
     &           +dhot*
     &           exp(Xbstmxx(ilo,iup,iso)*dt1+Xbstmyy(ilo,iup,iso)*dt0)
            rxup=dhot*exp(Xsxx(ilo,iup,iso)*dt1+Xsyy(ilo,iup,iso)*dt0)
     &           +ratcon(6)*phion(ilo,iup,iso)
            do ite=1,ntem
               tev=tevm(ite)
               frc=fnem(ite)
               rxdn=rxdn
     &              +denow*frc*denow*frc*ratcon(10)*betaxx(ilo,iup,iso)
     &              +denow*frc*ratcon(7)*phrec(ilo,iup,iso)
               rxup=rxup+denow*frc*ratcon(9)*sxx(ilo,iup,iso)
               if(iolevel.eq.2) then
               write(89,901)iso,ilo,iso-1,iup,tev,frc,
     &              dhot*Xsxx(ilo,iup,iso),
     &              dhot*dhot*Xbetaxx(ilo,iup,iso),
     &              denow*frc*sxx(ilo,iup,iso),
     &              denow*frc*denow*frc*betaxx(ilo,iup,iso),
     &              dhot*Xbstmxx(ilo,iup,iso),
     &              +denow*frc*phrec(ilo,iup,iso)
               endif
            enddo
            if(iolevel.eq.2.and.ntem.eq.0)then
               write(89,901)iso,ilo,iso-1,iup,tev,frc,
     &              dhot*Xsxx(ilo,iup,iso),
     &              dhot*dhot*Xbetaxx(ilo,iup,iso),
     &              dhot*Xbstmxx(ilo,iup,iso),
     &              phion(ilo,iup,iso)   
            endif
         else if(ltype(itran).eq.3)then
            rxdn=dhot*
     &           exp(Xcapexx(ilo,iup,iso)*dt1+Xcapeyy(ilo,iup,iso)*dt0)
            rxup=ratcon(11)*autoxx(ilo,iup,iso)
            do ite=1,ntem
               tev=tevm(ite)
               frc=fnem(ite)
               rxdn=rxdn
     &              +denow*frc*ratcon(12)*capexx(ilo,iup,iso)  
               if(iolevel.eq.2) then
               write(89,902)iso,ilo,iso-1,iup,tev,frc,
     &              dhot*Xcapexx(ilo,iup,iso),
     &              autoxx(ilo,iup,iso),
     &              denow*frc*capexx(ilo,iup,iso)
               endif
            enddo
            if(iolevel.eq.2.and.ntem.eq.0)then
               write(89,902)iso,ilo,iso-1,iup,tev,frc,
     &              dhot*Xcapexx(ilo,iup,iso),
     &              autoxx(ilo,iup,iso)
            endif
         endif
         a(ind2,ind1)=rxup
         a(ind1,ind2)=rxdn
         a(ind1,ind1)=a(ind1,ind1)-a(ind2,ind1)
         a(ind2,ind2)=a(ind2,ind2)-a(ind1,ind2)
 20   continue
c
 200  continue
      close(89)

c___equation a*poptime = dndt returns the solution
c

  171 continue

      do i = 1,neqnow
         ff(i) = 0.00
         do j=1,neqnow
            ff(i) = ff(i)+a(i,j)*poptime(j)
         enddo

      enddo
      ntem=ntemold

 101  format('elev', 2i4,1x,a5,1x,f14.4,1x,1p,e16.8,0p,1x,10i3,1x,i4)
 900  format('e',2(i3,i5),1p,10e10.2)
 901  format('i',2(i3,i5),1p,10e10.2)
 902  format('a',2(i3,i5),1p,10e10.2)
 903  format(i3,i5,1p,10(1x,e8.1),/,10(3x,10(1x,e8.1),/))

      return 
      end

c
c
      subroutine tedt(neqnow,timeODE,poptime,ff)
c     this subroutine computes matrix elements for a given temperature
c     all the time-relevant quantities are compared on original time grid
c
c*****************************
c*                           *
c*  d e c l a r a t i o n s  *
c*                           *
c*****************************
c

      implicit double precision (a-h,o-z)
      include 'mainstuf'
      include 'timestuf'
      include 'popstuf'
      include 'runstuf'
      include 'xc_stuf'
c
      common /paramnow/ trnow,tinow,zbarnow,totpnow,sizenow,rhonow

      dimension ff(numpp),poptime(numpp)

      common /outrt/filnn
      character*10 filnn
c
c*******************
c*  begin routine  *
c*******************
c
c___find brackets for the temperature and the density

c      if (timeODE .eq. timeLAST) then
c         go to 171
c      endif

      if(zbarneq(qqdum).lt.0.d0)then
         print *,'Z < 0', tev, denow, zbarneq(qqdum),qqdum
      endif

      timeLAST = timeODE
      timer = timeODE

      
      trev = (tr0*(tnext-timer)+tr1*(timer-tlast))/(tnext-tlast)
      if(trev.le.trmin) trev=0.d0
      tr = trev*11604.5

      zmix = (zmix0*(tnext-timer)+zmix1*(timer-tlast))/(tnext-tlast)
      zelse = zmix*(percent/(1.-percent))

      
c... determine the Ne from the input Ne, or pops if have only Nt or Rho
      
      if (iruntyp1 .eq. 'ne') then
         dinow = (totp0*(tnext-timer)+totp1*(timer-tlast))/
     &        (tnext-tlast)
         denow = exp((log(d0)*(tnext-timer)+log(d1)*(timer-tlast))/
     1        (tnext-tlast))
         
      else if (iruntyp1 .eq. 'nt') then
         dinow=exp((log(totp0)*(tnext-timer)+log(totp1)*
     1        (timer-tlast))/(tnext-tlast))
         denow = exp((log(totp0)*(tnext-timer)+log(totp1)*
     1        (timer-tlast))/(tnext-tlast))*(1.-percent)
     2        *(zelse+zbarneq(qqdum)) ! note zelse=zmix*percent/(1-percent)
      else if (iruntyp1 .eq. 'rho') then
         dinow = exp((log(rh0)*(tnext-timer)+log(rh1)*(timer-tlast))
     1        /(tnext-tlast))*(1.-percent)/atmass(iatomz)
         denow =exp((log(rh0)*(tnext-timer)+log(rh1)*(timer-tlast))
     1        /(tnext-tlast))*(1.-percent)/atmass(iatomz)*
     2        (zelse+zbarneq(qqdum))
      endif

      do ii=1, ntem
         fnem(ii)= exp((log(fnem0(ii))*(tnext-timer)
     1        +log(fnem1(ii))*(timer-tlast))/(tnext-tlast))
         tevm(ii)=(tevm0(ii)*(tnext0-timer)+tevm1(ii)*(timer-tlast0))
     1        /(tnext0-tlast0)
      enddo
      
      if(cxflag) then
         sum=0.d0
         do i=1,nebins
            fe(i) = exp((log(fe0(i))*(tnext-timer)
     1           +log(fe1(i))*(timer-tlast))/(tnext-tlast))
            !sum = sum + sqrt(energy(i))*fe(i)*denergy(i)
            sum = sum + fe(i)*denergy(i)
         enddo
         do i=1, nebins
           fe(i)=fe(i)/sum
         enddo
         nefmax=1
         do i=1, nebins
            if(fe(i).gt.1.e-25)then
               nefmax=i
            endif
         enddo
         if(fhot.gt.0.d0) then
            fhot= exp((log(fhot0)*(tnext-timer)
     1           +  log(fhot1)*(timer-tlast))/(tnext-tlast))
            dhot=fhot
            if(fhot.le.1.1d0)dhot = fhot*denow
            !write(ioo,*)'fhot:',fhot,'fhot0',fhot0,'fhot1',fhot1
            !write(ioo,*)'tnext, timer, tlast',tnext,timer,tlast
            write(ioo,'(" !!CX : dhot =",1p,2e12.4)')dhot,denow
            write(16,'(" CX : dhot =",1p,2e12.4)')dhot,denow
         else
            write(ioo,'(" *** Check fhot in -dndt-: exiting")')
            write(16,'(" Warning:Check fhot in -dndt-: exiting")')
            return
         endif
      else
         dhot=0.d0
      endif

      ntemold=ntem
      if(rtflag) then
         if(ntem.eq.0) then
            write(ioo,'(" *** Check ntem in -dndt-: exiting")')
            write(16,'(" Warning:Check ntem in -dndt-: exiting")')
            return
         endif
      else
         ntem=0
      endif


c
c___zero array a(i,j)
c
      do j=1,numpop
         do i=1,numpop
            a(i,j) = 0.00
         enddo
      enddo

      if(iolevel .eq. 2) then
      open(89,file=filnn)
      write(89,'(1p,e15.5,0p,f15.4,1p,2e15.5)')timer,tevm(1),dhot,denow
      do iso=1,iatomz
         write(89,'("enot",i5,1x,2f16.5)')iso, eipz(iso)
         do i=1, nlev(iso)
            ilv=indlev(i,iso)
            if(indl2p(ilv).ne.0)
     &           write(89,101)iso,i,kname(ilv),elev(ilv),glev(ilv),
     &           (noctot(ilv,ip),ip=1,nprim),mobtot(ilv)
         enddo
      enddo
      endif

      do iso=isomin,isomax
         if(cxflag) then
            call Vcaafill(iso)
         endif
      enddo


      eebabs=0.d0 ! calculate energy gain by electron beam 
      ephabs=0.d0 ! calculate energy gain by photon absorption  

      if(radflag.ne.'off') goto 100

      do 10 itran=1, ntran
         llo=llower(itran)
         lup=lupper(itran)
         ind1=indl2p(llo)
         ind2=indl2p(lup)
         if(ind1*ind2.eq.0) goto 10
         ilo=levl(llo)
         iup=levl(lup)
         iso=isolev(llo)
         rxdn=0.d0
         rxup=0.d0
         if(ltype(itran).eq.1)then
            rxdn=ratcon(1)*aratexx(iup,ilo,iso)
     &           +dhot*Vcxx(iup,ilo,iso)
            rxup=dhot*Vcxx(ilo,iup,iso)
            do ite=1,ntem
               tev=tevm(ite)
               frc=fnem(ite)
               !rxdn=rxdn+denow*frc*ratcon(5)*cxx(iup,ilo,iso)
               !rxup=rxup+denow*frc*ratcon(4)*cxx(ilo,iup,iso)
            enddo
            det=elev(lup)-elev(llo)
            eebabs=eebabs+det*dhot*(
     &           Vcxx(iup,ilo,iso)*poptime(ind1) ! gain by electron beam
     &           -Vcxx(ilo,iup,iso)*poptime(ind2)) ! loss by electron beam
            ephabs=ephabs
     &           -ratcon(1)*det*aratexx(iup,ilo,iso)*poptime(ind2) !radiative loss
         else if(ltype(itran).eq.2)then
            iso=isolev(llo)
            rxdn=dhot*dhot*Vbetaxx(ilo,iup,iso)
     &           +dhot*Valfxx(ilo,iup,iso)
            rxup=dhot*Vsxx(ilo,iup,iso)
            do ite=1,ntem
               tev=tevm(ite)
               frc=fnem(ite)
!              rxdn=rxdn
!    &              +denow*frc*denow*frc*ratcon(10)*betaxx(ilo,iup,iso)
!    &              +denow*frc*ratcon(7)*phrec(ilo,iup,iso)
               if(poptime(ind2).gt.1.e-30)then
                  ephabs=ephabs-denow*ratcon(7)*estm*poptime(ind2) ! radiative loss
               endif
!              rxup=rxup
!    &              +denow*frc*ratcon(9)*sxx(ilo,iup,iso)
            enddo
            det=elev(lup)-elev(llo)
            if(poptime(ind1).gt.1.e-30)then
               eebabs=eebabs+det*dhot*(
     &              Vsxx(ilo,iup,iso)*poptime(ind1) ! gain by electron beam
     &              -Vbetaxx(ilo,iup,iso)*dhot*poptime(ind1)) ! energy loss
            endif
            if(poptime(ind2).gt.1.e-30)then
               eebabs=eebabs+det*dhot*(
     &              -Valfxx(ilo,iup,iso)*poptime(ind2))
            endif

         else if(ltype(itran).eq.3)then
            rxdn=dhot*Vcapexx(ilo,iup,iso)
            rxup=ratcon(11)*autoxx(ilo,iup,iso)
            do ite=1,ntem
               tev=tevm(ite)
               frc=fnem(ite)
!              rxdn=rxdn
!    &              +denow*frc*ratcon(12)*capexx(ilo,iup,iso)
            enddo
            det=elev(lup)-elev(llo)
            eebabs=eebabs+det*dhot*(
     &           Vcapexx(ilo,iup,iso)*poptime(ind1))
         endif
         a(ind2,ind1)=rxup
         a(ind1,ind2)=rxdn
         a(ind1,ind1)=a(ind1,ind1)-a(ind2,ind1)
         a(ind2,ind2)=a(ind2,ind2)-a(ind1,ind2)
 10   continue

      goto 200

 100  continue

      if (radflag .eq. 'trfile') then
         do 1 i=1,ngrid
            barjnow(i)=(barj0(i)*(tnext-timer)+barj1(i)*(timer-tlast))
     &           /(tnext-tlast)
 1       continue
         barjend = 0.00
         if (xbarj .gt. 0.00.and.trev.gt.0.d0) then
            xbytrev = xev(ngrid)/trev
            if (xbytrev .lt. 100.) then
               barjend = barjnow(ngrid)/xev(ngrid)**3*exp(xbytrev)
            endif
         endif
      endif


      do 20 itran=1, ntran
         llo=llower(itran)
         lup=lupper(itran)
         ind1=indl2p(llo)
         ind2=indl2p(lup)
         ilo=levl(llo)
         iup=levl(lup)
         iso=isolev(llo)
         if(ind1*ind2.eq.0) goto 20
         rxdn=0.d0
         rxup=0.d0
         if(ltype(itran).eq.1)then
            atrad=bradxx(ilo,iup,iso)
            rxup=dhot*Vcxx(ilo,iup,iso)
     &                 +ratcon(3)*atrad
            rxdn=ratcon(1)*aratexx(iup,ilo,iso)+dhot*Vcxx(iup,ilo,iso)
     &                 +ratcon(2)*atrad*glev(llo)/glev(lup)
            det=elev(lup)-elev(llo)
            if(poptime(ind1).gt.1.e-30.and.poptime(ind2).gt.1.e-30)then
               eebabs=eebabs+det*dhot*(
     &              Vcxx(iup,ilo,iso)*poptime(ind1) ! gain by electron beam
     &              -Vcxx(ilo,iup,iso)*poptime(ind2)) ! loss by electron beam
               ephabs=ephabs+ratcon(3)*etrad*poptime(ind1)
     &              -ratcon(2)*etrad*glev(llo)/glev(lup)*poptime(ind2)
     &              -ratcon(1)*det*aratexx(iup,ilo,iso)*poptime(ind2) 

c               write(89,900)iso,ilo,iso,iup,det, etrad,
c     &              -etrad*glev(llo)/glev(lup),
c     &              -det*aratexx(iup,ilo,iso),
c     &              poptime(ind1),poptime(ind2),ephabs
            endif
            do ite=1,ntem
               tev=tevm(ite)
               frc=fnem(ite)
!              rxdn=rxdn+denow*frc*ratcon(5)*cxx(iup,ilo,iso)
!              rxup=rxup+denow*frc*ratcon(4)*cxx(ilo,iup,iso)
               if(iolevel.eq.2) then
               write(89,900)iso,ilo,iso,iup,tev,frc,
     &              dhot*Vcxx(ilo,iup,iso),
     &              dhot*Vcxx(iup,ilo,iso),
     &              denow*frc*cxx(ilo,iup,iso),
     &              denow*frc*cxx(iup,ilo,iso),
     &              aratexx(iup,ilo,iso),
     &              bradxx(iup,ilo,iso),bradxx(ilo,iup,iso)
               endif
           enddo
         else if(ltype(itran).eq.2)then
            rxdn=dhot*dhot*Vbetaxx(ilo,iup,iso)
     &           +dhot*Valfxx(ilo,iup,iso)
     &           +dhot*Vbstmxx(ilo,iup,iso)
            rxup=dhot*Vsxx(ilo,iup,iso)
     &           +ratcon(6)*phion(ilo,iup,iso)
            
            if(poptime(ind1).gt.1.e-30)then
               det=elev(lup)-elev(llo)
               eebabs=eebabs+det*dhot*(
     &              Vsxx(ilo,iup,iso)*poptime(ind1) ! gain by electron beam
     &              -Vbetaxx(ilo,iup,iso)*dhot*poptime(ind1)) ! energy loss
               ephabs=ephabs
     &              +ratcon(6)*econt*poptime(ind1) ! radiative gain
c               write(89,901)iso,ilo,iso-1,iup,det,econt,
c     &              poptime(ind1),ephabs
            endif
            if(poptime(ind2).gt.1.e-30)then
               eebabs=eebabs+det*dhot*(
     &              -Valfxx(ilo,iup,iso)*poptime(ind2)
     &           -Vbstmxx(ilo,iup,iso)*poptime(ind2)) ! energy loss
            endif

            do ite=1,ntem
               tev=tevm(ite)
               frc=fnem(ite)
!              rxdn=rxdn
!    &              +denow*frc*denow*frc*ratcon(10)*betaxx(ilo,iup,iso)
!    &              +denow*frc*ratcon(7)*phrec(ilo,iup,iso)
!              rxup=rxup+denow*frc*ratcon(9)*sxx(ilo,iup,iso)
               if(poptime(ind2).gt.1.e-30)then
                  ephabs=ephabs-denow*ratcon(7)*estm*poptime(ind2) ! radiative loss
c                  write(89,901)iso,ilo,iso-1,iup,det,-denow*estm,
c     &                 poptime(ind2),ephabs
               endif
               if(iolevel.eq.2) then
               write(89,901)iso,ilo,iso-1,iup,tev,frc,
     &              dhot*Vsxx(ilo,iup,iso),
     &              dhot*dhot*Vbetaxx(ilo,iup,iso),
     &              denow*frc*sxx(ilo,iup,iso),
     &              denow*frc*denow*frc*betaxx(ilo,iup,iso),
     &              dhot*Vbstmxx(ilo,iup,iso),
     &              phion(ilo,iup,iso),   
     &              denow*frc*phrec(ilo,iup,iso)
               endif
            enddo
            if(iolevel.eq.2.and.ntem.eq.0)then
               write(89,901)iso,ilo,iso-1,iup,tev,frc,
     &              dhot*Vsxx(ilo,iup,iso),
     &              dhot*dhot*Vbetaxx(ilo,iup,iso),
     &              dhot*Vbstmxx(ilo,iup,iso),
     &              phion(ilo,iup,iso)   
            endif
         else if(ltype(itran).eq.3)then
            rxdn=dhot*Vcapexx(ilo,iup,iso)
            rxup=ratcon(11)*autoxx(ilo,iup,iso)
            if(poptime(ind1).gt.1.e-30)then
               eebabs=eebabs+det*dhot*(
     &              Vcapexx(ilo,iup,iso)*poptime(ind1))
            endif
            do ite=1,ntem
               tev=tevm(ite)
               frc=fnem(ite)
!              rxdn=rxdn
!    &              +denow*frc*ratcon(12)*capexx(ilo,iup,iso)  
               if(iolevel.eq.2) then
               write(89,902)iso,ilo,iso-1,iup,tev,frc,
     &              dhot*Vcapexx(ilo,iup,iso),
     &              autoxx(ilo,iup,iso),
     &              denow*frc*capexx(ilo,iup,iso)
               endif
            enddo
            if(iolevel.eq.2.and.ntem.eq.0)then
               write(89,902)iso,ilo,iso-1,iup,tev,frc,
     &              dhot*Vcapexx(ilo,iup,iso),
     &              autoxx(ilo,iup,iso)
            endif
         endif
         a(ind2,ind1)=rxup
         a(ind1,ind2)=rxdn
         a(ind1,ind1)=a(ind1,ind1)-a(ind2,ind1)
         a(ind2,ind2)=a(ind2,ind2)-a(ind1,ind2)
 20   continue
c
 200  continue
      close(89)

c___equation a*poptime = dndt returns the solution
c

  171 continue

      close(89)


      do i = 1,neqnow
         ff(i) = 0.00
         do j=1,neqnow
            ff(i) = ff(i)+a(i,j)*poptime(j)
         enddo
      enddo
      ntem=ntemold

      sumin = 0.d0
      sumne = 0.d0
      totin =0.d0
      sumpop=0.d0
      pop2=0.d0

      tev=tevm(1)
      do i = 1,neqnow
         ilev = indp2l(i)
         iso=ionstg(ilev)
         pop2=pop2+poptime(i)*iso**2
         sumpop=sumpop + poptime(i)
         sumin = sumin + ff(i)*efgrdp(ilev)            ! Ei*dYi/dt
         sumne = sumne + ff(i)*iso*1.5*tev     ! dNe/dt*3/2kT  
         totin = totin + efgrdp(ilev)*poptime(i)       ! total internal energy
      enddo

c     choice of internal energy should be looked at
      teint = totin
      tekin = 1.5*tev*denow
      totin2 =sumin 

c     adding free-free components
      zave=zbarneq(qqdum)
      ephab0=ephabs

      ephabf1=effabs(sumpop,pop2)  ! Kramer with gaunt factor
      ephabf2=effcoul(sumpop,pop2) ! Johnson and Dawson with Coulomb Log
      ephabf=ephabf1

c     total absorption rate

      ephabs=ephab0+ephabf

c      write(26,903)timer,tevm0(1),tevm1(1),denow,ephabs,ephabf1,
c     &     tephabs,teint,tekin,ephabs*(timer-ttlast)

 101  format('elev', 2i4,1x,a5,1x,f14.4,1x,1p,e16.8,0p,1x,10i3,1x,i4)

 900  format('e',2(i3,i5),1p,10e10.2)
 901  format('i',2(i3,i5),1p,10e10.2)
 902  format('a',2(i3,i5),1p,10e10.2)
 903  format(1p,e12.4,1x,2e11.4,11(1x,e8.1))
 904  format(i3,i5,1p,10(1x,e8.1),/,10(3x,10(1x,e8.1),/))

      return 
      end

      double precision function effabs(sumpop,pop2)
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'runstuf'
      include 'timestuf'
c
c___ calculate photoionization from hydrogenic state i
c
      effabs = 0.d0
      if(tev.le.1.e-5) return

c     compute free-free Gaunt factor from a simple fit to the
c     results of Karzas and Latter (Ap. J. Suppl. 6, 167 (1961))

      gam2lg=log10(13.6058*pop2/sumpop/tev)
      gbrem=1.+0.44*exp(-0.25*(gam2lg+0.25)**2)

c      gbrem=0.06  ! to make Al case consistent with cold exp.

      if (radflag .ne. 'trfile') then

c         effabs=denow*pop2*1.53d-25*tev**0.5  ! in ergs/cm3/s
         effabs=gbrem*denow*pop2*9.55d-14*sqrt(tev)  ! in eV/cm3/s

      else
c         field should be in [ergs/cm2/s/Hz]
c         constant includes: 4Pi*2.42e-37 from ff opacity

c         bc=denow*pop2*7.37d-22/tev**0.5 ! in ergs/cm3/s
         bc=gbrem*denow*pop2*4.5896d-10/tev**0.5  ! in eV/cm3/s

         ny=51
         ymin=xev(1)
         ylast=xev(ngrid)
         dely=(ylast-ymin)/real(ny-1)

         sum = 0.5*field(ymin)/ymin**3*(1.-exp(-ymin/tev))
         do k=2,ny-1
            y2 = ymin+dely*real(k-1)
            sum = sum+field(y2)/y2**3*(1.-exp(-y2/tev))
         enddo
         sum = sum+0.5*field(ylast)/ylast**3*(1.-exp(-ylast/tev))

         effabs=bc*sum*dely

      endif

c
      return
      end
c
      double precision function effcoul(sumpop,pop2)
c     inverse bremsstrahlung absorption rate by Dawson and Oberman
c     using Lee and More Coulomb Logarithm.
c     we take Coulomb log from CRETIN subroutine
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'runstuf'
      include 'timestuf'
c
c     twopi=6.283185
c     1eV=2.42e14 Hz
c     
      effcoul = 0.d0
      if(tev.le.0.01) return
      wpfe=5.64e4 *denow**0.5  ! rad/sec
      ewpfe=wpfe/6.283185/2.42e14 ! plasma frequency in eV


c     Johnson and Dawson [Phys. Fluids, 16, 722, (1973)]

      if (radflag .ne. 'trfile') then

      else
         bc=denow*pop2*2.538d-10/tev**1.5*coulog(pop2)  ! in eV/cm3/s
         ny=11
         ymin=xev(1)
         ylast=xev(ngrid)
         dely=(ylast-ymin)/real(ny-1)

         if(ewpfe.lt.ymin) then
            sum = 0.5*field(ymin)/ymin**2
     &           /sqrt(1.-ewpfe**2/ymin**2)
         endif
         if(ewpfe.lt.ylast) then
            sum=sum+ 0.5*field(ylast)/ylast**2
     &           /sqrt(1.-ewpfe**2/ylast**2)
         endif
         do k=2,ny-1
            y2 = ymin+dely*real(k-1)
            if(ewpfe.lt.y2) then
               sum = sum+field(y2)/y2**2
     &              /sqrt(1.-ewpfe**2/y2**2)
            endif
         enddo
         effcoul=bc*sum*dely

      endif
c
      return
      end
c
      function coulog(pop2)
      implicit double precision(a-h,o-z)
      include 'mainstuf'
      include 'runstuf'
      include 'timestuf'


      tfermi=(hp**2/8./emass)*(3.*denow/pi)**(2./3.)/evtoerg ! in eV
      tveff=sqrt(tev**2+tfermi**2)
      rdeby2=5.5205d5/(denow/tveff+pop2/tiev)  !Debye-Huckel length^2 cm2
      rion2=(0.75*zave/(pi*denow))**(2./3.)
      bmax2=max(rdeby2,rion2)
      bmincl=(zave*4.8d-8/tev)**2
      bminqm=2.51d-15/tev
      bmin2=max(bmincl,bminqm)
      coulog= max(0.5*log(1.d0+bmax2/bmin2),2.d0)
c      print *, coulog,tev,rdeby2,rion2,bmincl,bminqm
      return
      end
c------------------------------------------------------
c------- absorption cross-sections---------------------
c------------------------------------------------------
      subroutine abscrx(ratioNZ,iwabs)
      implicit real*8(a-h,o-z)
c     this routine returns opacity of radiation field
c     in units of cm-1

      include 'mainstuf'
      include 'popstuf'
      include 'timestuf'
      data  const/2.33333E-15/

      common /absorb/abb(100),abf(100),aff(100),att(100)
      common /phistuff/ cp0,cp1,cp2,cp3,xn,cep,cem
      dimension abfk(100),abfs(100)
      do i=1,ngrid
         abb(i)=0.d0
         abf(i)=0.d0
         aff(i)=0.d0
         att(i)=0.d0
      enddo

c ---- ff component ----

      zbr2=0.d0
      tpop=0.d0
      do kk = 0, iatomz
         sum=0.d0
         do i=1, nlev(kk)
            ilv=indlev(i,kk)
            sum=sum+popout(ilv)
         enddo
         zbr2=zbr2+(iatomz-kk)**2 * sum
         tpop=tpop+sum
      enddo

      cono=denow*zbr2*2.42e-37

      zbr2=zbr2/tpop            ! fractional population 

      do ii=1,ntem
         frc=fnem(ii)
         tev=tevm(ii)
         gff=1.d0+0.44*exp(-0.25*(0.25+log10(13.6*zbr2/tev))**2)
         do k=1,ngrid
            o = xev(k)
            aff(k)=aff(k)+cono*gff*frc*(1-exp(-o/tev))/sqrt(tev)/o**3
         enddo
      enddo


c---- bound-bound component

      cwd=4.6337d-5*sqrt(tiev/(2.*atomz))
      etr0=max(0.,(xev(1)-1000*cwd*xev(1)))


      do 100 iso=isomin,isomax

         do 11 iup=2, nlev(iso)
            ilv=indlev(iup,iso)
            if(indl2p(ilv).eq.0) goto 11
            pu=popout(ilv)
            gu=glev(ilv)

            do 12 ilo=1, iup-1
               ilw=indlev(ilo,iso)
               if(indl2p(ilw).eq.0) goto 12
               pl=popout(ilw)
               gl=glev(ilw)
               flu=fxx(ilo,iup,iso)
               if(flu.lt.1.e-10) goto 12
               etr=elev(ilv)-elev(ilw) ! transition energy in eV
               opl=0.0265*pl*flu*(1.d0-pu/pl*gl/gu)
               if(pl/tpop.lt.1.e-20) goto 12

               wd=cwd*etr
               avt=(gamrat(ilo,iso)+gamrat(iup,iso))/(wd*12.57*evtohz)
               do ii=1,ngrid
                  xx=(xev(ii)-etr)/wd
                  phi0=const/wd*voigt(xx,avt)
                  abb(ii) = abb(ii) + opl*phi0
               enddo
               if(iwabs.eq.1)
     &          write(28,'("bb",3i3,f12.4,1p,10e10.2)')iso,iup,ilo,
     &              etr,pl,0.0265*pl*flu,opl,abb(1)
 12         continue
 11      continue


c -----------bound-free component

         do 21 i=1, nlev(iso)
            ilo=indlev(i,iso)
            if(indl2p(ilo).eq.0) goto 21

            pl=popout(ilo)
            gl=glev(ilo)
            zeff=qlev(ilo)
            if(pl/tpop.lt.1.e-20.or.zeff.lt.1.e-30) goto 21

            do 22 j=1, nlev(iso-1)

               if(istrn(i,j,iso).eq.0)goto 22

               iup=indlev(j,iso-1)
               if(indl2p(iup).eq.0) goto 22

               eip=elev(iup)-elev(ilo)
               if(eip.lt.1.e-10)goto 22

               pn=petrn(i,j,iso)
               pu=popout(iup)
               gu=glev(iup)

               kramer=0
               ish=nish(i,j,iso)
               xn =cphilsh(5,ish,iso)
               if(ish.eq.0)then
                  kramer=1
               else if(xn.le.1.d-30)then
                  kramer=1
               endif

               if(kramer.eq.1)then
c ------------ Kramer---------------------
                  oxx=2.9140d-17*eip**2.5/zeff*pl*pn
                  do k=1,ngrid
                     o = xev(k)
                     if(o.gt.eip) then
                        do ii=1,ntem
                           frc=fnem(ii)
                           tev=tevm(ii)
                           olt=(1.66d-22*denow)*(pu/pl)*(gl/gu)/tev**1.5
                           abf(k)=abf(k)+frc*oxx/o**3
     &                          *(1.-olt*(pu/pl)*exp((eip-o)/tev))

                        enddo
                        if(k.eq.1.and.iwabs.eq.1)
     &                       write(28,'("bk",3i3,f12.4,1p,10e10.2)')
     &                       iso,i,j,
     &                       eip,pl,frc*oxx/o**3,frc*oxx/o**3*
     &                       (1.-olt*(pu/pl)*exp((eip-o)/tev)),abf(k)
                     endif
                  enddo
               else
c ------------ Scofield---------------------

                  ish=nish(i,j,iso)
                  cp0=cphilsh(1,ish,iso)
                  cp1=cphilsh(2,ish,iso)
                  cp2=cphilsh(3,ish,iso)
                  cp3=cphilsh(4,ish,iso)
                  xn =cphilsh(5,ish,iso)
                  cep=cphilsh(6,ish,iso) 
                  cem=cphilsh(7,ish,iso) 
                  conphi = 1.d-18*(13.606/cep)*xn*pl
                  do k=1,ngrid
                     o = xev(k)
                     if(o.gt.eip) then                  
                        xx=log(o/cep)
                        do ii=1,ntem
                           frc=fnem(ii)
                           tev=tevm(ii)
                           olt=(1.66d-22*denow)*(pu/pl)*(gl/gu)/tev**1.5
                           abf(k)=abf(k)+frc*conphi*phicx(xx)
     &                          *(1.-olt*exp((eip-o)/tev))
                        enddo
                        if(k.eq.1.and.iwabs.eq.1)
     &                       write(28,'("bs",3i3,f12.4,1p,10e10.2)')
     &                       iso,i,j,eip,pl,frc*conphi*phicx(xx),
     &                       frc*conphi*phicx(xx)*
     &                       (1.-olt*(pu/pl)*exp((eip-o)/tev)),
     &                       abf(k)
                     endif
                  enddo
               endif
 22         continue
 21      continue
 100  continue

      do k=1, ngrid
         abb(k)=ratioNZ*abb(k)
         abf(k)=ratioNZ*abf(k)
         att(k)=abb(k)+abf(k)+aff(k)
      enddo

      return
      end


c**********************************************************

