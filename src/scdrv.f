c  Modified for EK model
c
c           ******  *******  *******  ******     ***    ******
c          *        *           *      *    *   *   *    *    *
c          *        *           *      *    *  *     *   *    *
c          *   ***  ****        *      *****   *******   *****
c          *     *  *           *      *       *     *   *  *
c          *     *  *           *      *       *     *   *   *
c           ******  *******     *      *       *     *   *    *
c
C
C
      subroutine getpar
      implicit real*8 (a-h,o-z)
      SAVE
c
c___     this routine allows the user to specify the operating parameters,
c___     and control the execution of the code via the command lines input.
c___     available commands are:
c___
C  Z        Atomic number of species of interest
C  Initial  SS or LTE or FILE filename
C  History  filename NE
C           filename NE with RHO or NT
C           filename RHO or NT
C           GRID and RHO,NT,NE gets grid-type file
C  Fe       FILE: create the f(e) vs eV from history file
C           FEFILE filename for f(e) vs eV file
C  Fhot     FILE use the THOT FHOT columns in history file
C           FIXED # # indicates the fixed FHOT fraction & THOT (ev)
C           OFF turns off the FHOT calculation
C  Tr (ev)  # indicates fixed tr
C           OFF turns off the radiation
C           FILE: read tr from history file
C           TRFILE filename for user J vs eV file
C           DILUTION followed by dilution factor
C  Ti (ev)  # indicates fixed ti 
C           OFF yields ti = te 
C           FILE: read ti from history file
C           Ti/Te followed by ratio of Ti/Te factor
C  Opacity  SIZE # indicates fixed size
C           OFF yields optical thin limit
C           FILE: read size from history file
C  Mixture  Zbar %  [Atomic #] of other species')
C  Outfile  name: if blank then use default name, ie, timeZ
C  Evolve   SS or LTE or TD
C  Runfile  filename. Instructions read from the file
C  Examine  look at the current history file
C  Time     TimeBegin TimeEnd DeltaTime
C  Time     log TimeBegin TimeEnd DeltaTime
C  IO       Level start# stop# delta#
C  run      start execution using information provided
C  help     produces list of current info
C           INFO will cause a command list
C  end
c___
c___     only one command per input line may be supplied.
c
c
      include 'mainstuf'
      include 'runstuf'
      include 'timestuf'
      include 'popstuf'
      include 'xc_stuf'
c
      dimension temp(5)
      logical around,readfile,open,histfile
      character filnam*10
      character ichar*2, namep*3, nameo*4
      character*120 msg

      parameter(numcom=20)
      character*8 icom(numcom),icol(10)
      character*20 isym(35),isizefix,runfile
      integer lsym(35)
      integer filnamln

      integer minargs(numcom),maxargs(numcom),itemp(5)
      equivalence (temp,itemp)
c
      data numnums/' ','0','1','2','3','4','5','6','7','8','9','X'/
      data column/'time','te  ','ne  ','rho ','nt  ',
     1     'tr  ','zbar','size','ti  ','thot',
     2     'fhot','mte ','mfe ','dilu','zmix','perc'/
c
      data filnamln /File_Name_Length/
      data ichar /'00'/
      data nameflag / 0/
      data nameo / 'info'/, namep/ 'tfe'/
      data itrout/'off'/, timezero/0.00/, timestop/0.00/, deltat/0.00/
      data logtime/0/
      data itiout/'off'/
      data timename/'   '/,ievolve/'ss'/
      data runfile/'    '/
      data isizefix /'                    '/
      data ifirst /1/, iolevel/0/
      data readfile/.false./


      data minargs / 2,2,2,2,2, 2,3,1,2,3, 3,1,1,1,2, 2,2,2,3,2/
      data maxargs / 3,3,4,5,3, 3,4,2,2,4, 7,1,2,1,2, 2,5,4,5,5/

      data icom /'z'       ,'initial' ,'fe      ','tr      ','ti      ',
     1           'opacity','mixture' ,'outfile' ,'evolve'  ,'history' ,
     2           'time'   ,'run'     ,'help'    , 'end'    ,'examine' ,
     3           'runfile','io'      ,'fhot'    , 'isos'   , 'lattice'/

      common /codename/ version
      character*20 version 
c
c
c___prompt for input & accept it.
c
      iolevel=0
      readfile=.false.

c reading runfile from the command line 

      call getarg(1,msg)
      call getsymb (msg,120,isym,nsym,lsym,35)
      if (nsym .eq. -1) write (ioo,*)' $$$ input symbol > 20 characters'
      if (nsym .eq. 0) then
         write (ioo,*)' $$$  No input file: Exiting'
         stop
      endif
      runfile = isym(1)
      inquire(file=runfile,exist=around,opened=open)
      if (open) then
         write(ioo,'(" $$$ ",a20," is opened.Something wrong")') runfile
         stop
      else if (around) then
         call initzero
         open(unit=50,file=runfile,status='old')
         histfile = .false.
         call runread(histfile,nameflag)
         if (histfile) then
           readfile = .true.
         else if (.not. readfile) then
           timename = '   '
         endif
         close(unit=50,status="keep")
      else
         write(ioo,*)' **** runfile not around'
         stop
      endif

c initializing the run
  750 continue
      if (.not. readfile) then
         write(ioo,*)' $$$ Must have history file before running'
         stop
      endif

c---- if the % was specified for the mass not the number density
c         modify the zelse to account for this

      amratio = 1.00
      if (elseznum .gt. 0.00) then
         iother = elseznum
         amratio = atmass(iatomz)/atmass(iother)
         zelse = zbartemp*(percent/(1.-percent))*amratio
      else if (elseznum .le. 0.00 .and. iruntyp1 .eq. 'rho' .and.
     1         zelse .gt. 0.00) then
        write(ioo,*)' $$$ Specify atomic # of other species for runtype'
        write(ioo,*)'     rho.  To modify re-enter MIXTURE command'
        stop
      endif

      if (timename .eq. 'grid') then
        if (ievolve .eq. 'td') then 
          write(ioo,*)' $$$ For grid option evolve must be LTE of SS'
          ievolve = initflag
        endif
        if (ievolve .ne. initflag) then
          write(ioo,*)' $$$ For grid option INITIAL must equal EVOLVE'          
          initflag = ievolve
        endif
      endif
       
c
c---- create the running times for the problem
c
 7511 if (timezero .eq. 0.00 .and. timestop .eq. 0.00) then
         timezero = timein(1)
         timestop = timein(ntimein)
         ntimes = ntimein
         if (ntimein .gt. NMTout) then
           write(ioo,*)' $$$ # of times too large. Will use evenspacing'
           ntimes = NMTout
           deltat = (timestop-timezero)/real(NMTout-1)
           do 7512 i=1,NMTout
             times(i) = timezero+deltat*real(i-1)
 7512      continue
         else
           do 755 i=1,ntimes
             times(i) = timein(i)
  755      continue
         endif

      else

         if (timezero.lt.timein(1).or.timestop.gt.timein(ntimein)) then
            if(timestop-timein(ntimein).gt.1.e-10*timestop)then
               write(ioo,*)' $$$ Start or stop time out of range'
               timezero = 0.00
               timestop = 0.00
               deltat = 0.00
               go to 7511
            endif
         endif

         if (deltat .eq. 0.00) then
            ncotime = 0
            do 756 i=1,ntimein
               if(timezero.ge.timein(i).and.timezero.lt.timein(i+1))then
                  ncotime = 1
                  times(1) = timezero
                  do 757 j=i+1,ntimein
                     if (timestop .ge. timein(j)) then
                        ncotime = ncotime+1
                        if (ncotime .gt. NMTout) then
                           write(ioo,*)' $$$ # of times too large #'
                           deltat = (timestop-timezero)/real(NMTout-1)
                           go to 7511
                        else
                           times(ncotime) = timein(j)
                        endif
                     endif
 757              continue
               endif
 756        continue
            if (times(ncotime) .lt. timestop) then
               ncotime = ncotime+1
               times(ncotime) = timestop
            endif
            ntimes = ncotime
            
         else if (deltat .ne. 0.00) then

            if (logtime .eq. 0) then
               ncotime = deltat
               deltat= (timestop-timezero)/(ncotime-1)
               ncotime = min(ncotime,NMTout)
               do 759 i=1,ncotime
                  times(i) = timezero+deltat*(i-1)
                  if (times(i) .ge. timestop) then
                     times(i) = timestop
                     ntimes = i
                     go to 758
                  endif
 759           continue
               ntimes = ncotime
 758           continue
            else if (logtime.eq.1) then
               ncotimes = int(deltat)
               if (ncotimes .gt. NMTout) then
                  ntimes = NMTout
                  deltat = real(ntimes)
               endif
               dellogt = (timestop/timezero)**(1./(deltat-1.))
               ntimes = deltat
               do 761 i=1,ntimes
                  times(i) = timezero*dellogt**(i-1)
                  if (times(i) .ge. timestop) then
                     times(i) = timestop
                     ntimes = i
                     go to 762
                  endif
 761           continue
 762           continue
            else if (logtime.eq.2) then
               ncotime=deltlg
               deltt=(timepause-timezero)/(ncotime-1)
               do i=1,ncotime
                  times(i) = timezero+deltt*(i-1)
               enddo
               ncotimes = int(deltat)
               if (ncotimes .gt. NMTout-ncotime) then
                  ntimes = NMTout
                  deltat = real(ntimes)
               endif
               dellogt=(timestop/timepause)**(1./ncotimes)
               ntimes=ncotimes+ncotime
               do 763 i=1,ncotimes
                  times(i+ncotime) = timepause*dellogt**(i)
                  if (times(i+ncotime) .ge. timestop) then
                     times(i+ncotime) = timestop
                     ntimes = i+ncotime
                     go to 764
                  endif
 763           continue
 764           continue
            else if (logtime.eq.3) then
               ncotime=int(deltlg)
               dellogt=(timepause/timezero)**(1./(ncotime-1.))
               do i=1,ncotime
                  times(i) = timezero*dellogt**(i-1)
               enddo
               ncotimes = deltat
               if (ncotimes .gt. NMTout-ncotime) then
                  ncotimes = NMTout-ncotime
               endif
               deltt=(timestop-timepause)/(ncotimes)
               ntimes=ncotimes+ncotime
               do 765 i=1,ncotimes
                  times(i+ncotime) = timepause+deltt*(i)
                  if (times(i+ncotime) .ge. timestop) then
                     times(i+ncotime) = timestop
                     ntimes = i+ncotime
                     go to 766
                  endif
 765           continue
 766           continue
            endif
         endif
      endif
c     
c___  size and IPD ok, flush and close previous files.
c
      if (ifirst .eq. 0) then
         close(unit=21,status='keep')
         close(unit=16,status='keep')
      else
         ifirst = 0
      endif
c
c___ If the FILE option, which names the output files, has not been used
c___ create the output file names.  The default is 'time' and 'info'
c___ with the ion number appended to the end to form the complete name.
c
  480 If (nameflag .eq. 0) then

         nametime='outfile.dat'
         nameinfo='info.dat'

      endif

c
c___destroy and create the files, as needed.
c

  495 open(unit=21,file=nametime)
      open(unit=16,file=nameinfo)


c
c___ test to see if trfile file must be read. If so go to readtr
c
      if (radflag .eq. 'trfile') then
         ireadit = 0
         call subrdtr(ireadit)
         if (ireadit .eq. 0) then
            write(ioo,*)' $$$ Difficulty reading TRFILE: radflg set off'
            write(16,*)' Exiting:Difficulty reading TRFILE'
            stop
            radflag = 'off'
            itrout = 'off'
         endif
      endif
c
c___ For Multi-temperature or Non-thermal components
c
      

      if(feflag .ne. 'off') then
         if(feflag(1:4).eq.'file') then
            readit = 0.0
            call subrdfe(readit)
            if (readit .eq. 0.00) then
               write(ioo,*)' $$$ Difficulty with feFILE: stop execution'
               write(16,*)' Exiting:Difficulty with feFILE: ',
     &              'stop execution'
               stop
            endif
         else if(feflag(1:4).eq.'func') then
c implement analytical form for fe(E) : not ready yet.
c           call subfunc(nfunc)
         endif
      endif

c
c___output the parameters for this calculation.
c
      if (sizefix .eq. 0.00 .or. opacflag .eq. 'off') then
         isizefix = '                    '
      else if (opacflag .eq. 'file') then
         isizefix = '                    '
      else
         write(isizefix,'(1pe8.2)') sizefix
      endif

      if(isomin.eq.0.and.isomax.eq.0) then
         isomin=1
         isomax=iatomz
      endif

      rtflag =.true.
      cxflag =.false. 
      if(feflag .ne. 'off') then
         cxflag =.true.
         if(feflag(5:8).eq.'only') rtflag=.false.
      endif


      write(21,802) iatomz,iteflag,initflag,initfile,itrout,dilution,
     1              itiout,tibyte,fhflag,fhotset,thotset,
     3              feflag,fefile,opacflag,isizefix,zbartemp,
     2              pctemp,elseznum,nametime,nameinfo,ievolve,
     2              timename,iruntyp1,iruntyp2,timezero,timestop,
     3              deltat,runfile,isomin,isomax,isamp(1),
     4              isamp(2),nxprim,r_lattice,version
  802 format('C',/,
     *'C VARIABLES  --      DESCRIPTION            - RUN VALUES',//,
     *'C  Z         -- atomic number: Te  option   - ',3x,i2,3x,i3/,
     *'C  Initial   -- choice of initial condition - ',3x,a4,2x,a20,/,
     *'C  Tr(ev) W  -- radiation temperature       - ',3x,a8,1pe9.2/,
     *'C  Ti(ev)    -- ion temperature  Ti/Te      - ',3x,a8,1pe9.2/,
     *'C  FHot      -- hot e-: Fhot & Thot(eV)     - ',3x,a8,1p,2e9.2/,
     *'C  Fe        -- e-beam: name of fe(E) file  - ',3x,a8,1x,a20/,
     *'C  Opacity   -- optical depth treatment     - ',3x,a4,1x,a10/,
     *'C  Mixture   -- Zbar  %  [Atomic #]         - ',3x,0p,3f8.2,/,
     *'C  Outfile   -- name of output files        - ',3x,a20,1x,a20,/,
     *'C  Evolve    -- evolve pop by SS, LTE or TD - ',3x,a3,/,
     *'C  History   -- definition of hydro data    - ',3x,a20,1x,2a3,/,
     *'C  Time      -- Time start stop [delta]     - ',3x,1p3e9.2,/,
     *'C  Runfile   -- name of file with input info- ',3x,a20,/,
     *'C  Isos      -- Min/Max iso for details     - ',3x,i10,3x,i10,/,
     *'C  Isamp     -- Collision data choice       - ',3x,i10,3x,i10,/,
     *'C  Nmax      -- Max. principal quantum no.  - ',3x,i10,/,
     *'C  LATTICE   -- LATTICE radius in cm        - ',3x,1p,e9.2,/,
     *'C  Version   -- ',a20,//)


c
c___go perform the requested calculation.
c
 100  return

c  error messages
 2250 write (ioo,2260)
 2260 format (' $$$ error creating filename. retry.')
      return
      
      end

      subroutine subrdfe(readit)
      implicit real*8 (a-h,o-z)

      include 'mainstuf'
      include 'timestuf'
      include 'runstuf'
      include 'xc_stuf'

      character*120 msg

      parameter(numcom=20)
      character*8 icom(numcom),icol(10)
      character*20 isym(35),isizefix,runfile
      integer lsym(35)

c
c   This routine reads the electron energy distribution in the units 
c
c                     n(v)dv = sqrt(e)*f(e)*de * fevfe(e)
c
c   Note that when n(v) is a maxwellian then 
c
c                f(e) = 2/sqrt(pi)*exp(-e/t)/t**1.5
c
c   the number of energy bins is nebins and all energys are in eV.
c


      open(unit=44,file=fefile,status='old')

      read(44,*) ntimefe,nebins
      igr=0
 20   read(44,'(a)')msg
      call getsymb (msg,120,isym,nsym,lsym,35)
      do i=1, nsym
         read(isym(i),'(bn,f20.0)',err=101) xxx
         igr=igr+1
         energy(igr)=xxx
      enddo
      if(igr.ne.nebins) goto 20

      do 30 j=1,ntimefe
         read(44,'(a)')msg
         call getsymb (msg,120,isym,nsym,lsym,35)
         read(isym(1),'(bn,f20.0)',err=101) xxx
         fetime(j)=xxx
         if(nsym.eq.1)then
            fevfe(j)=0.d0
         else
            read(isym(2),'(bn,f20.0)',err=101) xxx
            fevfe(j)=xxx
         endif
         igr=0
 40      read(44,'(a)')msg
         call getsymb (msg,120,isym,nsym,lsym,35)
         do i=1, nsym
            read(isym(i),'(bn,f20.0)',err=101) xxx
            igr=igr+1
            feIn(igr,j)=xxx
         enddo
         if(igr.ne.nebins) goto 40
   30 continue

      close(unit=44,status='keep')

c  Test if the times for the f(e) bracket the times for calculations
c  For grid and multe options, time is set up with an increment of 1e-12
c  Make fetime() include the timein() set up in -subgrid- and -submulte-

      if (fetime(ntimefe).lt.times(1) .or. 
     1     fetime(1).gt.times(ntimes)) then
         write(ioo,*)' $$$ fefile does not bracket calculation times'
         feflag = 'file'
         fefile = '    '
         stop
      endif
      
c  create and array of denergy weight such that the trapeziodal integral
c  is specified by SUM(i) f(i)*sqrt(E)*dE.  Also note that the first point
c  is uniques since we do not include the energy = 0, f(0) = 0.0 point but
c  use it to evaluate the weights.


c  energy(i) is defined at the center of the ith zone
c  denergy(i) is defined as the width of the ith zone

      icheck=0

      fe(0)=0.
      energy(0)=0.
      denergy(0)=0.
      denergy(1)= energy(1)+(energy(2)-energy(1))/2.d0
      sum = denergy(1)
      do i=2,nebins-1
         denergy(i) = (energy(i+1)-energy(i-1))/2.
         sum=sum + denergy(i)
	 if(abs(denergy(i)-denergy(i-1)).gt.1.e-5*denergy(i))
     &	 icheck=icheck+1
      enddo
      denergy(nebins)=2.*(energy(nebins)-energy(nebins-1))
     &               -denergy(nebins-1)

      inftyp=0
      if(icheck.gt.0)then
         write(16,'("CX: fe(i) is NOT evenly spaced")')
         infmin=10
         infmax=-10
         do i=1,nebins
            xfix=log10(energy(i))
            ifix=int(xfix)
            if(indfx(ifix).eq.0)indfx(ifix)=i
            infmin=min(infmin,ifix)
            infmax=max(infmax,ifix)
         enddo
         inftyp=1
         if(infmin.eq.infmax)then
            do i=1,nebins
               xfix=log(energy(i))
               ifix=int(xfix)
               if(indfx(ifix).eq.0)indfx(ifix)=i
               infmin=min(infmin,ifix)
               infmax=max(infmax,ifix)
            enddo
            inftyp=2
         endif
         do ifix=-10,10
            if(ifix.lt.infmin)then
               indfx(ifix)=1
            else if(ifix.gt.infmax)then
               indfx(ifix)=nebins
            endif
         enddo
         write(16,'("CX: inftype, infmin,infmax",3i6)')
     &        inftyp,infmin,infmax
      endif	


      if(feflag.eq.'filefe'.and. initflag.eq.'ss')then
         do j=1, ntimein
            tevnow=tein(j)
            t32=tevnow**1.5
            do i=1, nebins
               feIn(i,j)=1.12838*exp(-energy(i)/tevnow)/t32
            enddo
         enddo
      endif

c  create an array of temperatures from the integration of the f(e)

      do j=1,ntimefe
         dsum = 0.00
         eave = 0.00
         do k=1,nebins
            ck = feIn(k,j)*sqrt(energy(k))*denergy(k)
            dsum = dsum+ck
            eave = eave+ck*energy(k)
            feIn(k,j)=max(feIn(k,j),1.d-70)
         enddo
         tevfe(j) = 2./3.*eave/dsum
         write(16,'("CX : Thot & Nhot = ",1p,2e12.4)')tevfe(j),dsum


c  ensure normalization is maintained for intregration as it is performed
         if(abs(dsum-1.d0).gt.1.d-10) then
            write(16,'("CX: fe(e) normalized from ",i3,1p,e10.2)')j,dsum 
            if(fevfe(j).lt.1.d0)then
               fevfe(j)=1.d0
            endif
         endif
         do k=1,nebins
            feIn(k,j) = feIn(k,j)/dsum
         enddo
      enddo

      do j=1, ntimefe
         write(16,'("CX: info ",1p,3e12.4)')fetime(j),tevfe(j),fevfe(j)
      enddo

c     determine the minimum transition energy included in the model
c      exthr= (energy(2)+energy(3))/2.d0
      exthr= 0.d0

      readit = 1.0

 100  continue
      return

 101  write(ioo,'(" $$$ Error in reading Fe file")')
      write(16,'(" Warning:Error in reading Fe file")')
      ireadit = 0
      return
      end


      subroutine subrdtr(ireadit)
      implicit real*8 (a-h,o-z)


      include 'mainstuf'
      include 'timestuf'
      include 'runstuf'

      dimension trevin(NMTin)

      character*120 msg

      parameter(numcom=20)
      character*8 icom(numcom),icol(10)
      character*20 isym(35),isizefix,runfile
      integer lsym(35)

c   This routine reads the mean intensity to be used in the radiation
c   field calculations.  The number of energy point for the field is
c   ngrid and the number of times in the file is ntimtr.  The energy is
c   in eV and the mean intensity is in ergs/s/cm^2/Hz.  There must be
c

      open(unit=33,file=trfile,status='old')

c *** note that one should include 3rd number 
c     if the radiation field is extrapolated

      xbarj=0.d0
 10   read(33,'(a)')msg
      call getsymb (msg,120,isym,nsym,lsym,35)
      if(msg(1:1).eq.'c'.or.msg(1:1).eq.'#') goto 10
      read(isym(1),'(bn,f20.0)',err=101) xxx
      ntimtr=int(xxx)
      read(isym(2),'(bn,f20.0)',err=101) xxx
      ngrid=int(xxx)
      if(nsym.eq.3)then
         xbarj=1.d0
      endif

      igr=0
 20   read(33,'(a)')msg
      call getsymb (msg,120,isym,nsym,lsym,35)

      do i=1, nsym
         read(isym(i),'(bn,f20.0)',err=101) xxx
         igr=igr+1
         xev(igr)=xxx
      enddo
      if(igr.ne.ngrid) goto 20

      do 30 j=1,ntimtr
         read(33,*) rdtime(j),trevin(j)
         igr=0
 40      read(33,'(a)')msg
         call getsymb (msg,120,isym,nsym,lsym,35)
         do i=1, nsym
            read(isym(i),'(bn,f20.0)',err=101) xxx
            igr=igr+1
            barjin(igr,j)=xxx
         enddo
         if(igr.ne.ngrid) goto 40
   30 continue

      close(unit=33,status='keep')


c  Test if the times for the field bracket the ytimes for calculations
c  and place the values for the trevin in the array trin at the correct
c  times.  If there are no trevin calculate one.

      if (rdtime(ntimtr).lt.times(1).or.rdtime(1).gt.times(ntimes)) then
        write(ioo,*)' $$$ Radfile doesnt bracket calculation times'
        itrout = 'off'
        radflag = 'off'
        stop
      endif

c  now interpolate the trevin on to trin

      do 4 k=1,ntimein
         trin(k) = 0.00
         timenow = timein(k)
         do 5 j=2,ntimtr
            if (timenow .le. rdtime(j) .and.
     1          timenow .ge. rdtime(j-1)) then
               trin(k) = ((rdtime(j)-timenow)*trevin(j-1)+
     1                    (timenow-rdtime(j-1))*trevin(j))/
     2                   (rdtime(j)-rdtime(j-1))

               go to 4
            endif
    5    continue
    4 continue


c  ensure that there is a column of data dummied in for tr

      if (icolumn(6) .eq. 0) then
         column(ncolumn) = 'tr'
         icolumn(6) = 6
      endif

      delevtr = xev(2)-xev(1)
      ireadit = 1

 100  continue
      return

 101  write(ioo,'(" $$$ there is an error in reading radiation field")')
      write(16,'(" Warning:there is an error in reading ",
     &     "radiation field")')
      ireadit = 0
      return
      end

      subroutine JakeMake(tenowio,denowio,timeio)
      implicit real*8 (a-h,o-z)

      include 'mainstuf'
      include 'timestuf'
      include 'popstuf'
      include 'runstuf'

      data aijlow/1.00e-04/

c   if jacob is on then determine the maximum value of the a(i,j)
c   and the value of the elements that define the upper and lower bands
c   called aijmax, mu, ml respectively. The magnitude of the smallest
c   element is aijlow; thus no element below aijmax/aijlow is counted

      aijmax = 0.00
      do 900 i=1,neq
        do 900 j=1,neq
          aijmax = max(aijmax,abs(a(i,j)))
  900 continue
      smallaij = aijmax*aijlow
      munow = 0
      mlnow = 0
      do 905 i=1,neq
         do 906 j=1,neq

            absaij = abs(a(i,j))
            if (absaij .eq. 0.00) then
              jacob(i,j) = -9999
            else
              jacob(i,j) = int(log10(absaij/aijmax)+9.99999999)
            endif

           if (absaij .gt. smallaij) then
             if (j .gt. i .and. j-i .gt. munow) then
               munow = j-i
             endif
             if (i .gt. j .and. i-j .gt. mlnow) then
               mlnow = i-j
             endif
           endif

  906   continue
  905 continue

      if (iolevel .gt. 1) then

        write(16,'(//,"Jacobian at Te, Ne, time",1p,3e10.3,
     1                  " Aijmax ",e10.3," ml mu ",2i5,/)')
     2                  tenowio,denowio,timeio,aijmax,mlnow,munow

        do 902 i=1,neq
          do 901 j=1,neq
            if (jacob(i,j) .eq. -9999) then
              hjacob(j) = numnums(0)
            else if (jacob(i,j) .lt. 0) then
              hjacob(j) = "."
            else
              hjacob(j) = numnums(jacob(i,j)+1)
            endif
  901        continue
            write(16,'(1x,i4,3x,100a1)') i,(hjacob(j),j=1,neq)
  902   continue

      endif
      return
      end

      subroutine subfile(readfile)

      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'timestuf'
      include 'popstuf'
      include 'runstuf'

      logical readfile

      character*120 msg
      character*20 isym(35)
      dimension lsym(35)

      dimension hydrdat(NMTin,ncolumn)
      equivalence (hydrdat,timein(1))
c
c  read time dependent hydro data on unit 10:  first read header &
c       determine the variables in the file and where they are
c       the allowed variables are contained in the array 'column'.
c       In the file all columns not to be read in MUST be labeled
c       but not with one of the allowed variables
c       Blank lines and lines that start with 'c' or '*' are ignored


    1 read (10,'(a)',err=200) msg
      call getsymb(msg,120,isym,nsym,lsym,20)
      if (nsym .eq. 0) go to 1
      if (isym(1)(1:1) .eq. 'c' .or. isym(1)(1:1) .eq. '*') go to 1
      if (nsym .eq. -1)  then
         write (ioo,*)' $$$ Number of words > maximum'
         go to 1
      endif

      icount = 0
      do 2 i=1,ncolumn
         icolumn(i) = 0
         do j=1,NMTin
            hydrdat(j,i) = 0.00
         enddo
    2 continue

      nmte=0
      nmfe=0
      do 5 i=1,nsym
         do 4 j=1,ncolumn
            if (column(j) .eq. isym(i)) then
               if(j.eq.12)nmte=nmte+1
               if(j.eq.13)nmfe=nmfe+1
               if (icolumn(j) .eq. 0) then
                  icolumn(j) = i
                  icount=icount+1
               else if(j.ne.12.and.j.ne.13) then
                  write(ioo,*)' $$$ A Column appears more then once'
                  icount = 0
                  readfile = .false.
                  return
               endif
            else if(column(j) .eq. isym(i)(1:4))then
               if(j.eq.14)then
                  if (icolumn(j) .eq. 0) then
                     icolumn(j) = i
                     icount=icount+1
                  endif
               endif
            endif
    4    continue
    5 continue

c   Do tests for the consistency of the columns read in and
c      the calculational scheme chosen in HISTORY.

      if (icolumn(1) .eq. 0)  then
         write(ioo,*)' $$$ Time column must be given in HISTORY file'
         readfile = .false.
         return
      endif
      if (icolumn(2).eq.0.and.icolumn(12).eq.0)  then
         write(ioo,*)' $$$ Te column must be given in HISTORY file'
         readfile = .false.
         return
      endif
      if (icolumn(10)*icolumn(12) .ne. 0)  then
         write(ioo,*)' $$$ Use fhot for bi-maxwellian, mte for multe-Te'
         readfile = .false.
         return
      endif
      if (icolumn(10).ne.0.and.fhflag.eq.'off') then
         write(ioo,*)' $$$ fhot is not set: include "fhot file" '
         write(ioo,*)' $$$ or remove fhot data in file'
         readfile = .false.
         return
      endif
      if (icolumn(8).eq.0.and.opacflag.eq.'file') then
         write(ioo,*)' $$$ size is not set: include in history file'
         write(ioo,*)' $$$ or correct opacity option in runfile'
         readfile = .false.
         return
      endif
      if (icolumn(12).ne.0) then
         fhflag='multe'
      endif

      if (nmte.ne.nmfe) then
         write(ioo,*)
     &        ' $$$ Data points for Thot and Fhot are inconsistent'
         readfile = .false.
         return
      endif

c   now test if the history file is consistent with the iruntypes

      do 6 i=1,ncolumn
        if ((iruntyp1 .eq. column(i) .or. iruntyp2 .eq. column(i))
     1              .and. icolumn(i) .eq. 0) then
           write(ioo,*)' $$$ iruntype not found in history file'
           readfile = .false.
           return
        endif
    6 continue

c  Now read in all the information until EOF read

      ntimenow = 0
   10 continue
      read (10,'(a)',end=30,err=200) msg
      call getsymb(msg,120,isym,nsymow,lsym,20)
      if (nsymow .eq. 0) go to 10
      if (isym(1)(1:1) .eq. 'c' .or. isym(1)(1:1) .eq. '*') go to 10
      if (nsymow .eq. -1)  then
         write (ioo,*)' **** Number > 20 places.'
         go to 10
      endif
      
      ntimenow = ntimenow+1
      if (ntimenow .gt. NMTin) then 
         write(ioo,'(" $$$ Maximum # of times in history is ",I5)')
     1        ntimenow-1
         go to 200
      endif

      hydrdat(ntimenow,14)=dilution ! initial value
      hydrdat(ntimenow,15)=zmix ! initial value
      hydrdat(ntimenow,16)=percent ! initial value


      do 20 i=1,ncolumn
         if (icolumn(i) .ne. 0) then
            read(isym(icolumn(i)),'(bn,f20.0)') hydrdat(ntimenow,i)
            if(i.eq.16) hydrdat(ntimenow,i)=hydrdat(ntimenow,i)/100.
         endif
 20   continue

      ntemx=nmte
      nxx=icolumn(12)-1

      if (nsymow.ne.nsym) then
         if(icolumn(12).eq.0) then
            write(ioo,*)' $$$ # of columns .ne. # of column labels'
            readfile = .false.
            return
         else
            ntemx=(nsymow-nxx)/2
         endif
      endif

      if(fhflag.eq.'file') then
         ntemin(ntimenow)=2
         tevmin(1,ntimenow)=hydrdat(ntimenow,2)
         fnemin(1,ntimenow)=1.d0-hydrdat(ntimenow,11)
         tevmin(2,ntimenow)=hydrdat(ntimenow,10)
         fnemin(2,ntimenow)=hydrdat(ntimenow,11)
      else if(fhflag.eq.'fixed') then
         ntemin(ntimenow)=2
         tevmin(1,ntimenow)=hydrdat(ntimenow,2)
         fnemin(1,ntimenow)=1.d0-fhotset
         tevmin(2,ntimenow)=thotset
         fnemin(2,ntimenow)=fhotset
      else if(fhflag.eq.'multe') then
         ntemin(ntimenow)=ntemx
         do i=1, ntemx
            read(isym(nxx+i),'(bn,f20.0)') tevmin(i,ntimenow)
         enddo
         if(icolumn(2).eq.0)then
            icolumn(2)=icolumn(12)
            hydrdat(ntimenow,2)=tevmin(1,ntimenow)
         endif
         totf=0.d0
         do i=1, ntemx
            read(isym(nxx+ntemx+i),'(bn,f20.0)')fnemin(i,ntimenow)
            totf=totf+fnemin(i,ntimenow)
         enddo
         if(abs(totf-1.d0).gt.1.d-10)then
            write(ioo,'(1pe10.2," $$$ fractions renormalized ")')totf
            do j= 1, ntemx
               fnemin(i,ntimenow)=fnemin(i,ntimenow)/totf
            enddo
         endif
      else
         ntemin(ntimenow)=1
         tevmin(1,ntimenow)=hydrdat(ntimenow,2)
         fnemin(1,ntimenow)=1.d0
      endif

      go to 10

 30   ntimein = ntimenow
      close(unit=10, status='keep')
      readfile = .true.
      return


 200  write(ioo,*)' $$$ Read error. Is history file format correct?'
      close(unit=10, status='keep')
      readfile = .false.
      return

      end

      
      subroutine submulte(readfile,fromrun)

      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'timestuf'
      include 'popstuf'
      include 'runstuf'

      logical readfile,fromrun

      character*120 msg
      character*20 isym(35)
      dimension lsym(35)

      dimension hydrdat(NMTin,ncolumn)
      equivalence (hydrdat,timein(1))

      dimension ntemx(10),fnemx(10,10),tevmx(10,10)


      if (iatomz .eq. 0) then
        write(ioo,*)' $$$ For multe option must input Z before HISTORY'
        stop
      endif

      if (.not. fromrun) then
 1       write(ioo,2) iruntyp1
 2       format(1x,' MULTE option chosen must input:',/,
     1      10x,a3,' followed by min max delta',/,
     1      10x,'MTe  followed by  Te(1)  Te(2) ... Te(n) ',/,
     1      10x,'MFe  followed by  fNe(1) fNe(2)... fNe(n) ',/,
     1      10x,'in the order of fNe(1)>fNe(2)>... >fNe(n) ',/,
     1      10x,'END will stop the read in and return.'//)
      endif

      ithere = 0
      idhere = 0

  100 continue 

      if (fromrun) then
         read(50,'(a)',end=1000) msg
      else
         Write(ioo,'("  GRID OK: ",$)')
         if(batch) then
            read (18,'(a)') msg
         else
            read (ioi,'(a)') msg
         endif
      endif

      call getsymb (msg,120,isym,nsym,lsym,35)
      if (nsym .eq. -1) write (ioo,*)' $$$ input symbol > 20 characters'
      if (nsym .le. 0) then
         if (fromrun) then
            backspace(unit=50)
            go to 1000
         else          
            go to 100
         endif
      endif
      if (isym(1) .eq. 'end') then
         go to 200
      endif

      if (isym(1) .eq. iruntyp1) then
         if (nsym .lt. 4) then
            write(ioo,*)' $$$ Must have min max and delta on input line'
            if (fromrun) then
               backspace(unit=50)
               go to 1000
            else          
               go to 100
            endif
         endif
      else if (isym(1) .eq. 'mte' .or. isym(1) .eq. 'mfe') then
         if (nsym .le.2 ) then
            write(ioo,*)' $$$ Must have at least two temperatures'
            if (fromrun) then
               backspace(unit=50)
               go to 1000
            else          
               go to 100
            endif
         endif
      else
         write(ioo,'(" $$$ Must end first")')
         if (fromrun) then
            backspace(unit=50)
            go to 1000
         else          
            go to 100
         endif
      endif
c
c
c*******************************************
c*  number of multiple temperatures  N     *
c*******************************************
c

c*******************************************
c*  d e n s i t y   dmin dmax (10**) ddel  *
c*******************************************
c
      if (isym(1) .eq. iruntyp1) then
         read(isym(2),'(bn,f20.0)') denmin
         read(isym(3),'(bn,f20.0)') denmax
         if (isym(4) .eq. '10**' ) then
            read(isym(5),'(bn,f20.0)') dellog
            dendel=10.0**dellog
         else
            read(isym(4),'(bn,f20.0)') dendel
            if (dendel .le. 1.) then
               ndens = 1
               dellog = 1
               dendel = 10.
               idhere = 1
               go to 100
            else
               dellog=log10(dendel)
            endif
         endif

         ndens=max(1,int((log10(denmax)-log10(denmin))/dellog+1.4))
         if (ndens .gt. 10) then
            write(ioo,'(" $$$ The max # of ",a3,
     1           " points is 10. Try again")') iruntyp1
            if (fromrun) then
               backspace(unit=50)
               go to 1000
            else          
               go to 100
            endif
         endif
         idhere = 1
         go to 100
c
c*******************************************
c*  temperature     Te(1) Te(2) ... Te(N)   
c*******************************************
c
      else if (isym(1) .eq. 'mte') then

        ithere=ithere+1
        ntemx(ithere)=nsym-1
        do ii=2,nsym
           read(isym(ii),'(bn,f20.0)')tevmx(ii-1,ithere)
        enddo

        go to 100
c
c
c*******************************************
c*  d e n s i t y    fNe(1) fNe(2) ... fNe(N)   
c*******************************************
c
      else if (isym(1) .eq. 'mfe') then

         if(ithere.eq.0) then
            write(ioo,'("mte needs to be set up ")')
            goto 100
         else if(ntemx(ithere).ne.nsym-1) then
            write(ioo,'("mte and mfe points are different")')
            goto 1000
         endif

         totf=0.d0
         do ii=2,nsym
            read(isym(ii),'(bn,f20.0)')fnemx(ii-1,ithere)
            totf=totf+fnemx(ii-1,ithere)
         enddo
         if(abs(totf-1.d0).gt.1.d-10)then
            write(ioo, *) ' $$$ fractions are renormalized to unity'
            do j= 1, ntemx(ithere)
               fnemx(j,ithere)=fnemx(j,ithere)/totf
            enddo
         endif

         go to 100

      endif

  200 continue

      if (idhere .ne. 1 .or. ithere .eq. 0) then
        write(ioo,*)' $$$ Te and density data not determined. Returning'
        go to 1000
      endif

      do i=1,ncolumn
         icolumn(i) = 0
         do j=1,NMTin
            hydrdat(j,i) = 0.00
         enddo
      enddo

      icount = 3
      icolumn(1) = 1
      icolumn(2) = 2
      if (iruntyp1 .eq. 'ne') then
        icolumn(3) = 3
        ncolden = 3
      else if (iruntyp1 .eq. 'rho') then
        icolumn(4) = 3
        ncolden = 4
      else if (iruntyp1 .eq. 'nt') then
        icolumn(5) = 3
        ncolden = 5
      endif


c....  Create the time temperature and density arrays

      fhflag='multe'
      ntemps = ithere
      ntimein = ntemps*ndens
      tzero = 1.00e-12
      ncount=0
      do 300 i=1,ntemps
         tevnow = tevmin(1,i)
         ntem=ntemx(i)
         do 301 j=1,ndens
            ncount=ncount+1
            timein(ncount) = tzero*real(ncount-1)
            tein(ncount)=tevnow
            hydrdat(ncount,ncolden) = denmin*dendel**(j-1)
            ntemin(ncount)=ntem
            do k=1, ntem
               tevmin(k,ncount)=tevmx(k,i)
               fnemin(k,ncount)=fnemx(k,i)
            enddo
 301     continue
 300  continue
      numdens = ndens
      readfile = .true.

 1000 continue

      return
      end


      subroutine subgrid(readfile,fromrun)

      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'timestuf'
      include 'popstuf'
      include 'runstuf'
      include 'xc_stuf'


      logical readfile,fromrun

      character*120 msg
      character*20 isym(35)
      dimension lsym(35)

      dimension hydrdat(NMTin,ncolumn)
      equivalence (hydrdat,timein(1))

      if (iatomz .eq. 0) then
        write(ioo,*)' $$$ For grid option must input Z before HISTORY'
        go to 1000
      endif

      if (.not. fromrun) then
    1   write(ioo,2) iruntyp1,iruntyp1
    2   format(1x,' GRID option chosen must input:',/,
     1      10x,'Te  followed by min max delta',/,
     1      35x,'AND EITHER',/,
     1      10x,a3,' followed by min max delta',/,
     1      35x,'OR',/,
     1      10x,a3,' followed by min max 10** delta',/,
     1      10x,'END will stop the read in and return.'//)
      endif

      ithere = 0
      idhere = 0

  100 continue 

      if (fromrun) then
        if (idhere .eq. 1 .and. ithere .eq. 1) then
          go to 200
        else
          read(50,'(a)',end=1000) msg
        endif
      else
        Write(ioo,'("  GRID OK: ",$)')
        read (ioi,'(a)') msg
      endif

      call getsymb (msg,120,isym,nsym,lsym,35)
      if (nsym .eq. -1) write (ioo,*)' $$$ input symbol > 20 characters'
      if (nsym .le. 0) then
        if (fromrun) then
          backspace(unit=50)
          go to 1000
        else          
          go to 100
        endif
      endif

      if (isym(1) .eq. 'end') then
        go to 200
      endif

      if (isym(1) .eq. 'te' .or. isym(1) .eq. iruntyp1) then
        if (nsym .lt. 4) then
          write(ioo,*)' $$$ Must have min max and delta on input line'
          if (fromrun) then
            backspace(unit=50)
            go to 1000
          else          
            go to 100
          endif
        endif
      else
        write(ioo,'(" $$$ Must have: te ",a3," or end first")') iruntyp1
        if (fromrun) then
          backspace(unit=50)
          go to 1000
        else          
          go to 100
        endif
      endif
c
c
c*******************************************
c*  temperature     tmin tmax tdel         *
c*******************************************
c
      if (isym(1) .eq. 'te') then

        if (nsym .ne. 4) then
          write(ioo,*)' $$$ for te input only min max and delta'
          if (fromrun) then
            backspace(unit=50)
            go to 1000
          else          
            go to 100
          endif
        endif
        
        read(isym(2),'(bn,f20.0)') tempmin
        read(isym(3),'(bn,f20.0)') tempmax
        read(isym(4),'(bn,f20.0)') tempdel
        if (tempdel .le. 0.) then 
          tempdel = 0
          nttemp = 1
        else
          nttemp = max(1,int((tempmax-tempmin)/tempdel+1.4))
        endif

c
c____ check that the lower limit of the temperature is okay
c
        tevnow = tempmin-tempdel
        do 3000 it = 1,nttemp
          tevnow = tevnow+tempdel
          if ((iatomz .gt. 20) .and. 
     1         (tevnow .lt. tmin)) then
             write(ioo,*) 'Minimum temperature too low for Z'
             tempmin = tempmin+tempdel
             go to 3000
          else
            go to 3001
          endif
 3000   continue
          write(ioo,*)'The maximum Te is too low for this Z.'
          write(ioo,*)'Must be > ', int(tmin), 'eV, Try new range'
          go to 100
 3001   continue

        if (tempdel .le. 0.) then 
          ntemps = 1
        else
          ntemps = max(1,int((tempmax-tempmin)/tempdel+1.4))
        endif

        if (ntemps .gt. 10) then
          write(ioo,*)' $$$ The max # of te points is 10. Try again'
          if (fromrun) then
            backspace(unit=50)
            go to 1000
          else          
            go to 100
          endif
        endif
        ithere = 1
        go to 100
c
c
c*******************************************
c*  d e n s i t y   dmin dmax (10**) ddel  *
c*******************************************
c
      else if (isym(1) .eq. iruntyp1) then
        read(isym(2),'(bn,f20.0)') denmin
        read(isym(3),'(bn,f20.0)') denmax
        if (isym(4) .eq. '10**' ) then
          read(isym(5),'(bn,f20.0)') dellog
          dendel=10.0**dellog
        else
          read(isym(4),'(bn,f20.0)') dendel
          if (dendel .le. 1.) then
            ndens = 1
            dellog = 1
            dendel = 10.
            idhere = 1
            go to 100
          else
            dellog=log10(dendel)
          endif
        endif

        ndens=max(1,int((log10(denmax)-log10(denmin))/dellog+1.4))
        if (ndens .gt. 10) then
          write(ioo,'(" $$$ The max # of ",a3,
     1                " points is 10. Try again")') iruntyp1
          if (fromrun) then
            backspace(unit=50)
            go to 1000
          else          
            go to 100
          endif
        endif
        idhere = 1
        go to 100

      endif

  200 continue

      if (idhere .ne. 1 .or. ithere .ne. 1) then
        write(ioo,*)' $$$ Te and density data not determined. Returning'
        go to 1000
      endif

      do 201 i=1,ncolumn
         icolumn(i) = 0
         do 202 j=1,NMTin
            hydrdat(j,i) = 0.00
  202    continue
  201 continue

      icount = 3
      icolumn(1) = 1
      icolumn(2) = 2
      if (iruntyp1 .eq. 'ne') then
        icolumn(3) = 3
        ncolden = 3
      else if (iruntyp1 .eq. 'rho') then
        icolumn(4) = 3
        ncolden = 4
      else if (iruntyp1 .eq. 'nt') then
        icolumn(5) = 3
        ncolden = 5
      endif
      if (fhflag .eq. "fixed") then
        ncolthot = 10
        ncolfhot = 11
      endif

c....  Create the time temperature and density arrays

      ntimein = ntemps*ndens
      tzero = 1.00e-12
      ncount = 0
      do 300 i=1,ntemps
        tevnow = tempmin+tempdel*(i-1)
        do 301 j=1,ndens
          ncount = ncount+1
          timein(ncount) = tzero*real(ncount-1)
          tein(ncount) = tevnow
          hydrdat(ncount,ncolden) = denmin*dendel**(j-1)
          ntemin(ncount)=1
          tevmin(1,ncount)=tevnow
          fnemin(1,ncount)=1.d0
          if (fhflag .eq. 'fixed') then
             fnemin(1,ncount)=1.d0-fhotset
             ntemin(ncount)=2
             tevmin(2,ncount)=thotset
             fnemin(2,ncount)=fhotset
          endif
 301   continue
        hydrdat(ncount,14)=dilution
        hydrdat(ncount,15)=zmix
        hydrdat(ncount,16)=percent
 300  continue

      numdens = ndens
      readfile = .true.

 1000 continue

      return
      end



      subroutine runread(readfile,nameflag)
      implicit real*8 (a-h,o-z)
      SAVE
c
c___     this routine allows the user to specify the operating parameters,
c___     by using an input file specified with the runfile command.


      include 'mainstuf'
      include 'runstuf'
      include 'timestuf'
      include 'popstuf'
      include 'xc_stuf'
c
      dimension temp(5)

      logical around,readfile,open
      character*120 msg
      character filnam*10

      parameter(numcom=20)
      character*8 icom(numcom)
      character*20 isym(35)
      integer lsym(35)
      integer filnamln

      integer minargs(numcom),maxargs(numcom),itemp(5)

      equivalence (temp,itemp)
c
      data filnamln /File_Name_Length/
      data minargs / 2,2,2,2,2, 2,3,1,2,3, 3,1,1,1,1, 2,2,2,3,2/
      data maxargs / 3,3,4,5,3, 3,4,2,2,4, 7,1,2,1,1, 2,5,4,5,5/
      data icom /'z'       ,'initial' ,'fe      ','tr      ','ti      ',
     1 'opacity','mixture' ,'outfile' ,'evolve'  ,'history' ,'time'   ,
     2           'run'     ,'help'    , 'end'    ,'examine' ,'runfile',
     3           'io'      ,'fhot'    , 'isos'   , 'lattice'/
c
c
c___prompt for input & accept it.
c
  100 read (50,'(a)',end=1000) msg
      call getsymb (msg,120,isym,nsym,lsym,35)
      if(msg(1:1).eq.'c'.or.msg(1:1).eq.'*') go to 100
      if (nsym .eq. -1)  then
       write (ioo,'("  $$$ Input too long. Starts with ",a10)')msg(1:10)
       go to 1000
      endif
      if (nsym .eq. 0) go to 100
      if (isym(1)(1:1) .eq. 'c' .or. isym(1)(1:1) .eq. '*') then
        go to 100
      endif
c
c___now check for legality of command.
c
      do 130 ic=1,numcom
         if (isym(1) .eq. icom(ic)) go to 140
  130 continue

      write (ioo,'(" **** Unrecognizable command. Starts with",a10)')
     1      msg(1:10)
      go to 1000

c
c___are the minimum number of arguments for this command present?
c___if so, then go to the command handler...
c
  140 if (nsym .lt. minargs(ic)) then
         write (ioo,2060) ic, minargs(ic)
         write (16,2061) ic, minargs(ic)
 2060    format (' $$$ too few arguments for command #',i2,' (',i2,')')
 2061    format (' Warning:Too few arguments for command #',i2,' (',
     &        i2,')')
         go to 100
      endif

      if (nsym .gt. maxargs(ic)) then
         write (ioo,2062) ic,maxargs(ic)
         write (16,2063) ic,maxargs(ic)
 2062    format (' $$$ too many arguments for command #',i2,' (',i2,')')
 2063    format (' Warning:Too many arguments for command #',i2,' (',
     &        i2,')')
         go to 100
      endif


      go to (200,300,325,350,375,450,475,500,550,
     1     600,650,750,800,850,900,980,990,995,999,998),ic
c
c
c***************
c*  i o n   z  *
c***************
c
  200 read(isym(2),'(bn,f20.0)',err=2100) zion
      iatomz = zion
      isomin=1
      isomax=iatomz
      if(nsym.eq.3)iteflag=0
      go to 100
c
c
c****************************************
c*  initial condition calculation flag  *
c****************************************
c
  300 continue

      initfile = '        '
      if (isym(2) .eq. 'ss') then
         initflag = 'ss'
      else if (isym(2) .eq. 'lte') then
         initflag = 'lte'
      else if (isym(2) .eq. 'file') then
         initflag = 'file'
         if (nsym .ne. 3) then
            write(ioo,*)' $$$ Initial file name must be given.'
            write(16,*)' Warning:Initial file name must be given.'
            initfile = '   '
            stop
         else
            if (lsym(3) .gt. filnamln) then
               write(ioo,2315) lsym(3),filnamln
               write(16,2316) lsym(3),filnamln
               initflag = 'ss'
               go to 100
            endif

            initfile = isym(3)

c____ test if file = 'initfile' exists and if it is connected

            inquire(file=initfile,exist=around,opened=open)
            if (.not. around) then
               write(ioo,'(" $$$ ",a20, " not around")') initfile
               write(16,'(" Warning:",a20, " not around")') initfile
               initfile = '    '
               initflag = 'ss'
               stop
            else if (open) then
               write(ioo,'(" **** ",a20, " already opened")') initfile
               write(16,'(" Warning:",a20, " already opened")') initfile
               initfile = '    '
               initflag = 'ss'
               go to 100
            endif

c___ test to see if unit = 20 is connected. If it is disconnect
c         and save old file, then connect 'initfile' to it.

            inquire (unit=20,opened=open)
            if (open) then
               close(unit=20,status='keep')
            endif
         endif
      else
         write(ioo,*)' **** Initial options are: SS LTE FILE'
         write(16,*)' Warning:Initial options are: SS LTE FILE'
      endif

      go to 100
c
c
c*******************************************
c*  read in the electron distribution f(e) *
c*******************************************
c
  325 iplace = 2
      feflag ='off'

      if (isym(iplace) .eq. 'file'.or.isym(iplace).eq.'filefe') then
         if (nsym .eq. iplace) then
            write(ioo,*)' $$$ Name of file must come AFTER fefile'
            write(16,*)' Warning:Name of file must come AFTER fefile'
            stop
         endif
         if (lsym(iplace+1) .gt. filnamln) then
            write(ioo,2315) lsym(iplace+1),filnamln
            write(16,2316) lsym(iplace+1),filnamln
            go to 100
         endif
         iplace = iplace+1
         fefile = isym(iplace)
         

c____ test if file = fefile exists and if it is connected

         inquire(file=fefile,exist=around,opened=open)
         if (.not. around) then
            write(ioo,'(" $$$ ",a20, " not around")') fefile
            write(16,'(" Warning:",a20, " not around")') fefile
            fefile = '                   '  
            stop
         else if (open) then
            write(ioo,'(" *** ",a20, " already opened")') fefile
            write(16,'(" Warning:",a20, " already opened")') fefile
            go to 100
         endif

         if(isym(iplace-1).eq.'file')then
            feflag = 'file'
         else if(isym(iplace-1).eq.'filefe')then
            feflag='filefe'
         endif

         if(nsym.gt.iplace.and.isym(nsym).eq.'only') feflag = 'fileonly'

      else if(isym(iplace) .eq. 'func') then
         write(ioo,*)' **** fe is generated by analytical function'
         write(ioo,*)' **** not implemented yet'
      else
         write(ioo,*)' **** fe default is to create it from file Te'
      endif
      go to 100
c
c
c*******************************************
c*  read in the radiation temperature      *
c*******************************************
c
  350 itrout = 'off'
      radflag = 'off'
      fixhere = 0
      iplace = 1
  351 iplace = iplace+1
      If (isym(iplace) .eq. 'file') then
         if (radflag .ne. 'off') then
           write(ioo,*)' $$$ Incorrect tr input line'
           write(16,*)' Warning:Incorrect tr input line'
           radflag = 'off'
           itrout = 'off'
           stop
         endif
         radflag = 'file'
         itrout = 'file'
      else if (isym(iplace) .eq. 'trfile') then
         if (radflag .ne. 'off') then
           write(ioo,*)' $$$ Incorrect tr input line'
           write(16,*)' Warning:Incorrect tr input line'
           radflag = 'off'
           itrout = 'off'
           stop
         endif
         if (nsym .eq. iplace) then
            write(ioo,*)' $$$ Name of file must come AFTER trfile'
            write(16,*)' Warning:Name of file must come AFTER trfile'
            itrout = 'off'
            radflag = 'off'
            stop
         endif
         if (lsym(iplace+1) .gt. filnamln) then
            write(ioo,2315) lsym(3),filnamln
            write(16,2316) lsym(3),filnamln
            itrout = 'off'
            radflag = 'off'
            go to 100
         endif
         iplace = iplace+1
         trfile = isym(iplace)

c____ test if file = 'trfile' exists and if it is connected

         inquire(file=trfile,exist=around,opened=open)
         if (.not. around) then
            write(ioo,'(" $$$ ",a20, " not around")') trfile
            write(16,'(" Warning:",a20, " not around")') trfile
            trfile = '   '
            itrout = 'off'
            radflag = 'off'
            go to 100
         else if (open) then
            write(ioo,'(" $$$ ",a20, " already opened")') trfile
            write(16,'(" Warning",a20, " already opened")') trfile
            trfile = '    '
            itrout = 'off'
            radflag = 'off'
            go to 100
         endif
         radflag = 'trfile'
         itrout = trfile
      else if (isym(iplace) .eq. 'dilution') then
         if (nsym .eq. iplace) then
            write(ioo,*)' $$$ dilution factor must come AFTER dilution'
            write(16,*)' Warning:Dilution factor must come",
     &           " AFTER dilution'
            itrout = 'off'
            radflag = 'off'
            stop
         endif
         iplace = iplace+1
         read(isym(iplace),'(bn,f20.0)',err=2100) tempdil
         if (tempdil .le. 0. .or. tempdil .gt. 1.00) then
            write(ioo,*)' $$$ Dilution is lt 0 or gt 1.  Check.'
            write(16,*)' Warning:Dilution is lt 0 or gt 1.  Check.'
            itrout = 'off'
            radflag = 'off'
            stop
         endif
         dilution = tempdil
      else if (isym(iplace) .eq. 'off') then
         radflag = 'off'
         itrout = 'off'
      else
         if (fixhere .eq. 1 .or. radflag .ne. 'off') then
            write(ioo,*)' $$$ Fixed tr has 2 #s or line incorrect '
            write(16,*)' Warning:Fixed tr has 2 #s or line incorrect '
            radflag = 'off'
            itrout = 'off'
            stop
         endif
           
         read(isym(iplace),'(bn,f20.0)',err=2100) trdu
         if (trdu .gt. 1.e5 .or. trdu .lt. 0.00) then
            write(ioo,*)' $$$ Rad temperature can not be >100keV or <0.'
            write(16,*)' Warning:Rad temperature can not be >100keV', 
     &           ' or < 0.'
            radflag = 'off'
            itrout = 'off'
            stop
         endif
         fixhere = 1
         radflag = 'fixed'
         trfixed = trdu
         itrout = isym(iplace)
      endif
      if (iplace .lt. nsym) go to 351
      go to 100
c
c
c
c*******************************************
c*  read in the ion temperature            *
c*******************************************
c
  375 itiout = 'off'
      tiflag = 'off'
      tibyte = 1.00
      fixhere = 0
      iplace = 1
  376 iplace = iplace+1
      If (isym(iplace) .eq. 'file') then
         tiflag = 'file'
         itiout = 'file'
      else if (isym(iplace) .eq. 'ti/te') then
         if (nsym .eq. iplace) then
            write(ioo,*)' $$$ Ti/Te ratio must come AFTER ti/te'
            write(16,*)' Warning:Ti/Te ratio must come AFTER ti/te'
            itiout = 'off'
            tiflag = 'off'
            stop
         endif
         iplace = iplace+1
         read(isym(iplace),'(bn,f20.0)',err=2100) temprat
         if (temprat .le. 0.) then
            write(ioo,*)' $$$ Ti/Te is lt 0.  Check.'
            write(16,*)' Warning:Ti/Te is lt 0.  Check.'
            itiout = 'off'
            tiflag = 'off'
            stop
         endif
         tibyte = temprat
         itiout='ti/te'
         tiflag='ti/te'
      else if (isym(iplace) .eq. 'off') then
         tiflag = 'off'
         itiout = 'off'
      else
         if (fixhere .eq. 1) then
            write(ioo,*)' $$$ Fixed ti has # only and this is 2nd #'
            write(16,*)' Warning:Fixed ti has # only and this is 2nd #'
            tiflag = 'off'
            itiout = 'off'
            stop
         endif
           
         read(isym(iplace),'(bn,f20.0)',err=2100) tidu
         if (tidu .gt. 20000. .or. tidu .lt. 0.00) then
            write(ioo,*)' $$$ Ion temperature cant be >20keV or < 0.'
            write(16,*)' Warning:Ion temperature cant be >20keV', 
     &           ' or < 0.'
            tiflag = 'off'
            itiout = 'off'
            stop
         endif
         fixhere = 1
         tiflag = 'fixed'
         tifixed = tidu
         itiout = isym(iplace)
      endif
      if (iplace .lt. nsym .and. isym(iplace) .ne. 'off') go to 376
      go to 100

c
c
c***************
c   opacity    *
c***************
c

  450 continue

      if (isym(2) .eq. 'off') then
         opacflag = 'off'
      else if (isym(2) .eq. 'file') then
         opacflag = 'file'
      else if (isym(2) .eq. 'size') then
         if (nsym .lt. 3) then
            write(ioo,*)' $$$ Must give value of size when size used'
            write(16,*)' Warning:Must give value of size when size used'
            opacflag = 'off'
            stop
         endif

         read(isym(3),'(bn,f20.0)',err=2100) dummy
         if (dummy .lt.1.00e-06 .or. dummy .gt. 100.) then
            write(ioo,*)' $$$ Size value must be > 1e-6 and <100.'
            write(16,*)' Warning:Size value must be > 1e-6 and <100.'
            stop
         endif
         sizefix = dummy
         opacflag = 'size'
      else
         write(ioo,*)' $$$ Wrong parameters in opacity command'
         write(16,*)' Warning:Wrong parameters in opacity command'
         opacflag = 'off'
         stop
      endif

      go to 100
c
c**********************************************
c*  Mixture due to other species              *
c**********************************************
c
  475 continue
      read(isym(2),'(bn,f20.0)',err=2100) zbartemp
      read(isym(3),'(bn,f20.0)',err=2100) pctemp
      if (pctemp .ge. 100.) then
         write(ioo,*)' $$$ % of other species must be < 100'
         write(16,*)' Warning: % of other species must be < 100'
         stop
      endif
      if (nsym .eq. 4) then
         read(isym(4),'(bn,f20.0)',err=2100) elseznum
      endif

      percent =pctemp/100.

      zmix = zbartemp

      go to 100
c
c**********************************************
c*  File name of the output files             *
c**********************************************
c
  500 continue
      if (nsym .eq. 1) then
         nameflag = 0
      else
         if (lsym(2) .gt. filnamln) then
            write(ioo,2315) lsym(2),filnamln
            write(16,2316) lsym(2),filnamln
            go to 100
         endif
         nametime = isym(2)
         lchar = 16
         do 751 i=2,16
            if (isym(2)(i:i) .eq. ' ') then
               lchar = i-1
               go to 752
            endif
 751     continue
 752     continue
         nameinfo = isym(2)(1:lchar)//'info'
         nameflag = 1
      endif
      go to 100

c
c*************
c   evolve   *
c*************
c
  550 continue

      if (isym(2) .eq. 'ss') then
         ievolve = 'ss'
      else if (isym(2) .eq. 'lte') then
         ievolve = 'lte'
      else if (isym(2) .eq. 'td') then
         ievolve = 'td'
      else
         write(ioo,*)' $$$ Evolve options are: SS LTE TD. SS used'
         write(16,*)' Warning:Evolve options are: SS LTE TD. SS used'
         ievolve = 'ss'
      endif

      go to 100

c
c************
c  history  *
c************
c
  600 continue

      readfile = .false.

      if (isym(2) .eq. 'grid') then

        if (readfile) then
          write(ioo,*)' $$$ Can only have FILE or GRID  or MULTE'
          write(16,*)' Warning:Can only have FILE or GRID  or MULTE'
          timename = '   '
          stop
        endif
        timename = 'grid'
        readfile = .true.

      else if(isym(2) .eq. 'multe') then

        if (readfile) then
          write(ioo,*)' $$$ Can only have FILE or GRID or MULTE'
          write(16,*)' Warning:Can only have FILE or GRID or MULTE'
          timename = '   '
          stop
        endif
        timename = 'multe'
        readfile = .true.

      else 

        if (readfile) then
          write(ioo,*)' $$$ Can only have FILE or GRID  or MULTE'
          write(16,*)' Warning:Can only have FILE or GRID  or MULTE'
          timename = '   '
          stop
        endif
        if (lsym(2) .gt. filnamln) then
          write(ioo,2315) lsym(2),filnamln
          write(16,2316) lsym(2),filnamln
          go to 100
        endif
        timename = isym(2)

c____ test if file = 'timename' exists and if it is connected

        inquire(file=timename,exist=around,opened=open)
        if (.not. around) then
          write(ioo,'(" $$$ ",a20, " not around")') timename
          write(16,'(" Warning:",a20, " not around")') timename
          timename = '   '
          stop
        else if (open) then
          write(ioo,'(" **** ",a20, " already opened")') timename
          write(16,'(" Warning:",a20, " already opened")') timename
          timename = '    '
        endif
        readfile = .true.

      endif


      iplace = 2

      if (.not. readfile) then
        write(ioo,*)' $$$ FILENAME or GRID must be in HISTORY command!'
        write(16,*)' Warning:FILENAME or GRID must be in HISTORY ',
     &       'command!'
        stop
      endif
      iruntyp1 = ' '
      iruntyp2 = ' '

      do 602 i=3,nsym
        if (isym(i) .eq. 'ne') then
          if (iruntyp1 .eq. ' ') then
            iruntyp1 = 'ne'
          else
            iruntyp1 = ' '
            write(ioo,*)' $$$ Can not have NE as 2nd run type'
            write(16,*)' Warning:Can not have NE as 2nd run type'
            readfile = .false.
            timename = '   '
            stop
          endif
        else if (isym(i) .eq. 'rho'  .or. isym(i) .eq. 'nt') then
          if (iruntyp1 .eq. ' ') then
            iruntyp1 = isym(i)
          else if (iruntyp1 .eq. 'ne') then
            if (iruntyp2 .eq. ' ') then
              iruntyp2 = isym(i)
            else
              write(ioo,*)' $$$ Only 2 choices of runtype'
              write(16,*)' Warning:Only 2 choices of runtype'
              iruntyp2 = ' '
              iruntyp1 = ' '
              readfile = .false.
              timename = '    '
              stop
            endif
          endif
        else
          write(ioo,*)' $$$ Bad choice or order of runtypes'
          write(16,*)' Warning:Bad choice or order of runtypes'
          iruntyp1 = ' '
          iruntyp2 = ' '
          readfile = .false.
          timename = '    '
          stop
        endif
  602 continue

      if (iruntyp1 .eq. ' ') then
        write(ioo,*)' $$$ Must specify runtype on history input'
        write(16,*)' Warning:Must specify runtype on history input'
        iruntyp1 = ' '
        iruntyp2 = ' '
        timename = '    '
        stop
      endif

      if (timename .eq. 'grid') then
        readfile = .false.
        call subgrid(readfile,.true.)
        if (.not. readfile) then
          write(ioo,*)' $$$ Difficulty creating a grid type timefile'
          write(16,*)' Warning:Difficulty creating a grid type timefile'
          iruntyp1 = ' '
          iruntyp2 = ' '
          timename = '    '
          stop
        endif

      else if (timename .eq. 'multe') then

         readfile = .false.
         call submulte(readfile,.true.)
         if (.not. readfile) then
            write(ioo,*)' $$$ Difficulty creating a multe type timefile'
            write(16,*)' Warning:Difficulty creating a multe type ',
     &           'timefile'
            iruntyp1 = ' '
            iruntyp2 = ' '
            timename = '    '
            stop
         endif

      else

c___ test to see if unit = 10 is connected. If it is disconnect
c         and save old file, then connect 'timename' to it.

        open(unit=10,file=timename,status='old')

        readfile = .false.
        call subfile(readfile)


        if (.not. readfile) then
          write(ioo,'(" $$$ Difficulty reading file ",a20)') timename
          write(16,'(" Warning:Difficulty reading file ",a20)') timename
          iruntyp1 = ' '
          iruntyp2 = ' '
          timename = '    '
          close(unit=10,status='keep')
          stop
        endif

        close(unit=10,status='keep')
      endif

c
c____ Zero out any previous timestep info when the new file is read
c
      timezero = 0.00
      timestop = 0.00
      deltat = 0.00

      go to 100 
c
c***********
c*  Time   *
c***********
c
  650 continue
      deltat = 0.00
      logtime = 0
      if (isym(2) .eq. 'log') then
         ilow = 3
         if (nsym .ne. 5) then
           write(ioo,*)
     &           ' $$$ Must have t1, t2 & # per decade for log time'
           write(16,*)
     &           ' Warning:Must have t1, t2 & # per decade for log time'
            go to 100
         endif
         logtime = 1
      else if(isym(2) .eq. 'hybrid') then
         logtime = 2 !do linear first and log later
         ilow = 3
      else if(isym(2) .eq. 'pump') then
         logtime = 3 !do log first and linear later
         ilow = 3
      else 
         ilow = 2
      endif

      do 651 i=ilow,nsym
         read(isym(i),'(bn,f20.0)',err=2100) temp(i-ilow+1)
  651 continue

      if (temp(1) .ge. 0.00 .and. temp(1) .le. 1.00) then
         if (logtime .eq. 1 .and. temp(1) .eq. 0.00) then
           write(ioo,*)' $$$ Cannot have log time with initial time = 0'
           write(16,*)' Warning:Cannot have log time with initial ',
     &          'time = 0'
           logtime = 0
           stop
         endif
         timezero = temp(1)
      else
         write(ioo,*)' $$$ Time zero can not be > 1. or < 0.'
         write(16,*)' Warning:Time zero can not be > 1. or < 0.'
         stop
      endif

      if (temp(2) .gt. timezero .and. temp(2) .le. 1.00) then
         if(logtime .eq. 3) then
            timepause= temp(2)
            timestop = temp(3)
         else if(logtime .eq. 2) then
            timepause= temp(2)
            timestop = temp(3)
         else
            timestop = temp(2)
         endif
      else
         write(ioo,*)' $$$ Time stop can not be < timezero or > 1.0'
         write(16,*)' Warning:Time stop can not be < timezero or > 1.0'
         timezero = 0.00
         stop
      endif

      if (nsym .eq. 4 .and. logtime .eq. 0) then
         if (temp(3) .le. 0.00 .or. temp(3) .gt. NMTout) then
            write(ioo,'(" $$$ For time need ",i3,"> delta >0")')NMTout
            write(16,'(" Warning:For time need ",i3,"> delta >0")')
     &          NMTout
            deltat=NMTout
         else
            deltat = temp(3)
         endif
      else if (nsym .eq. 5) then
         if (temp(3) .le. 0.00 .or. temp(3) .gt. NMTout) then
           write(ioo,'(" $$$ For time need ",i3,"> delta >0")')NMTout
           write(16,'(" Warning:For time need ",i3,"> delta >0")')
     &          NMTout
            deltat=NMTout
         else
            deltat = temp(3)
         endif
      else if (nsym .eq. 6) then
         if (temp(4) .le. 0.00 .or. temp(4) .gt. NMTout) then
           write(ioo,'(" $$$ For time need ",i3,"> delta >0")')NMTout
           write(16,'(" Warning:For time need ",i3,"> delta >0")')
     &          NMTout
           deltat=NMTout
        else
           deltat=temp(4) 
        endif
      else if (nsym .eq. 7) then
         if (temp(4) .le. 0.00 .or. temp(4) .gt. NMTout) then
           temp(4) = NMTout/2
         else if (temp(5) .le. 0.00 .or. temp(5) .gt. NMTout) then
           temp(5) = NMTout/2
         endif
         deltlg = temp(4) 
         deltat = temp(5)
      endif

      go to 100

c
c***********
c*  r u n  *
c***********
c
c
  750 continue

      write(ioo,*)' $$$ RUN not valid inside RUNFILE. Ignored.'
      write(16,*)' Warning:RUN not valid inside RUNFILE. Ignored.'
      go to 100
c
c*************
c*  h e l p  *
c*************
c
  800 continue

      write(ioo,*)' $$$ HELP not valid inside RUNFILE. Ignored.'
      write(16,*)' Warning:HELP not valid inside RUNFILE. Ignored.'
      go to 100
c
c
c***********
c*  e n d  *
c***********
c
  850 continue
       return
c
c
c*******************
c*  e x a m i n e  *
c*******************
c

  900 continue

      write(ioo,*)' $$$ EXAMINE not valid inside RUNFILE. Ignored'
      write(16,*)' Warning:EXAMINE not valid inside RUNFILE. Ignored'
      go to 100

c*********************************
c                                *
c          RUNFILE               *
c                                *
c*********************************


  980 continue

      write(ioo,*)' $$$ RUNFILE not valid inside the RUNFILE. Ignored'
      write(16,*)' Warning:RUNFILE not valid inside the RUNFILE. ',
     &     'Ignored'
      go to 100

c*********************************
c                                *
c          IOlevel flags         *
c                                *
c*********************************


  990 continue
       read(isym(2),'(bn,f20.0)',err=2100) ziolevel
       if (ziolevel .ge. 0. .and. ziolevel .le. 6.) then
         iolevel = ziolevel
       else
         write(ioo,*)' $$$ IOlevel must be >= 0 and <= 6'
         write(16,*)' Warning:IOlevel must be >= 0 and <= 6'
       endif

       if (nsym .ge. 3) then
         read(isym(3),'(bn,f20.0)',err=2100) zio1st
         io1st = zio1st
         if (io1st .lt. 1. .or. io1st .gt. NMTout) then
           write(ioo,*)' $$$ io start # incorrect, set = to 1'
           write(16,*)' Warning:io start # incorrect, set = to 1'
           io1st = 1
         endif
         iostart = io1st
       endif

       if (nsym .ge. 4) then
         read(isym(4),'(bn,f20.0)',err=2100) zio2nd
         io2nd = zio2nd
         if (io2nd .le. io1st .or. io2nd .gt. NMTout) then
           write(ioo,*)' $$$ iostop incorrect,will be set to iostart+2'
           write(16,*)' Warning:iostop incorrect,will be set ',
     &          'to iostart+2'
           io2nd = iostart+2
         endif
         iostop = io2nd
       endif

       if (nsym .eq. 5) then
         read(isym(5),'(bn,f20.0)',err=2100) ziodel
         iodel = ziodel
         if (iodel .lt. 1 .or. iodel .gt. NMTout-1) then
           write(ioo,*)' $$$ iodelta incorrect, will be set to 1'
           write(16,*)' Warning:iodelta incorrect, will be set to 1'
           iodel = 1
         endif
         iodelta = iodel
       endif

       if (iostart .eq. 0) then
          iostart = 1
       endif
       if (iostop .eq. 0) then
          iostop = iostart+2
       endif
       if (iodelta .eq. 0) then
          iodelta = 1
       endif


       do 991 i=1,NMTout
          iouttime(i) = 0
  991  continue

       do 992 i=iostart,iostop,iodelta
          iouttime(i) = 1
 992   continue
       go to 100
c
c*******************************************
c*					   *
c*  read in the fhot fraction              *
c*                                         *
c*******************************************
c
 995  fhflag = 'off'
      fhotset = 0.00
      thotset = 0.00
      iplace = 1
      iplace = iplace+1
      If (isym(iplace) .eq. 'file') then
         if (iplace .ne. nsym) then
            write(ioo,*)' $$$ Fhot from FILE ,no other input allowed'
            write(16,*)' Warning:Fhot from FILE ,no other input allowed'
            go to 100
         endif           
         fhflag = 'file'
      else if (isym(iplace) .eq. 'fixed') then
         if (nsym .eq. iplace) then
            write(ioo,*)' $$$ FHOT & THOT ratio must come AFTER fixed'
            write(16,*)' Warning:FHOT & THOT ratio must come ',
     &           'AFTER fixed'
            go to 100
         endif
         iplace = iplace+1
         read(isym(iplace),'(bn,f20.0)',err=2100) tempfhot
         if (tempfhot .lt. 0.) then
            write(ioo,*)' $$$ fhot is lt 0.  Check.'
            write(16,*)' Warning:Fhot is lt 0.  Check.'
            go to 100
         endif
         fhotset = tempfhot
         iplace = iplace+1
         read(isym(iplace),'(bn,f20.0)',err=2100) tempthot
         if (tempthot .lt. 0.) then
            write(ioo,*)' $$$ Thot is lt 0.  Check.'
            write(16,*)' Warning:Thot is lt 0.  Check.'
            go to 100
         endif
         thotset = tempthot
         fhflag = 'fixed'
      endif
      go to 100
c
c***************
c*  i s o s    *
c***************
c
 999  read(isym(2),'(bn,f20.0)',err=2100) zsomin
      read(isym(3),'(bn,f20.0)',err=2100) zsomax
      isomin = int(zsomin)
      isomax = int(zsomax)
      if (isomax.gt.iatomz.or.isomin.gt.iatomz) then
         write(ioo,*)' **** check isos -> default: all detailed'
         write(16,*)' **** check isos -> default: all detailed'
         isomin=1
         isomax=iatomz
      endif
      if(nsym.gt.3)then
         read(isym(4),'(bn,f20.0)',err=2100) xprim!max.principal quantum number
         nxprim=int(xprim)
      endif

      go to 100
c
c
c***************
c*  Lattice radius    *
c***************
c
 998  read(isym(2),'(bn,f20.0)',err=2100) r_lattice
      print *, 'r_lattice', r_lattice
      go to 100
c
c
 1000 return
c
c*********************************
c*  e r r o r   m e s s a g e s  *
c*********************************
c
c

 2100 write (ioo,2110) i
      write (16,2111) i
 2110 format (' $$$ Argument # ',I2,' wrong. Input line ignored')
 2111 format (' Warning:Argument # ',I2,' wrong. Input line ignored')
      go to 100
c
 2315 format(' $$$ Input name is ',i2,' characters.',
     1       ' Limit is ',i2,'. input line ignored.')
 2316 format(' Warning:Input name is ',i2,' characters.',
     1       ' Limit is ',i2,'. input line ignored.')

      end

C
C
C************************************************************************
C*                                                                       
C*  ******  *******  *******   ******  *     *  *     *  ******          
C* *        *           *     *         *   *   **   **   *    *         
C* *        *           *     *          * *    * * * *   *    *         
C* *  ***  ****        *      *****      *     *  *  *   *****          
C* *     *  *           *           *     *     *     *   *    *         
C* *     *  *           *           *     *     *     *   *    *         
C*  ******  *******     *     ******      *     *     *  ******          
C*                                                                       
C************************************************************************
C
      subroutine getsymb (msg,msglen,isym,nsym,lsym,maxsym)
c
c*******************************************************************************
c
c__   this routine breaks a character string into space delimited
c__   symbols.  the symbols cannot be longer than 8 characters.
c__   the symbols are returned left justified.
c__
c__       msg       :  the character string to be broken.
c__
c__       msglen    :  the # of characters in the character string
c__
c__       isym     =:  the array to recieve the symbols.  (left justified)
c__
c__       nsym     =:  the # of symbols returned in 'isym'.  if = -1, then
c__                     a symbol was found that was > 20 characters.
c__
c__       lsym     =:  the array containing character length corresponding to
c__                    the isym array
c__
c__       maxsym    :  the maximum # of symbols to be returned in 'isym'.
c
c*******************************************************************************
c
c*****************************
c*                           *
c*  d e c l a r a t i o n s  *
c*                           *
c*****************************
c
      dimension lsym(35)
      character*1 itchar, iblank, comma, tab
      character*8  blank
      character*20 isym(35)
      character*120 msg
      data blank /'        '/
      data iblank / ' ' /
      data comma/','/
c
c*******************
c*  begin routine  *
c*******************
c
c   set initial conditions & clear symbol array.
c
      tab = char(09)
      iposcapa = ichar('A')
      iposcapz = ichar('Z')
      iposa = ichar('a')
      imove = iposcapa-iposa

      nsym = 0
      nclen = msglen
      do 10 i=1,maxsym
         isym(i) = blank
   10 continue
c
c   scan all characters in msg.
c
      i = 1
   20 continue
c
         itchar = msg(i:i)

         if (itchar .eq. iblank) go to 100
c
c   non-blank found, store as symbol until a blank.
c
         if (nsym .ge. maxsym) go to 200
         nsym = nsym+1
c
         k = 1
   50    if (k .gt. 20) go to 300
         ipos = ichar(itchar)
         if (ipos .ge. iposcapa .and. ipos .le. iposcapz) then
            ipos = ipos-imove
            itchar = char(ipos)
         endif

         isym(nsym)(k:k) = itchar
         k = k+1
         i = i+1
         if (i .gt. nclen) go to 100
         itchar = msg(i:i)
         if (itchar .ne. iblank) go to 50
         lsym(nsym) = k-1
c
  100 i = i+1
      if (i .le. nclen) go to 20
c
c   return to caller.
c
  200 return
c
c   too many words in message, set error.
c
  300 nsym = -1
      go to 200
c
c
      end

c
      subroutine AtTime(time)

      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'timestuf'
      include 'popstuf'
      include 'runstuf'
      include 'xc_stuf'
c
      dimension hydrdat(NMTin,ncolumn)
      equivalence (hydrdat,timein(1)) 

      dimension hydrAtT(ncolumn)
      equivalence (hydrAtT(1),TimeAtT)

      if (time .lt. timein(1) .or. time .gt. timein(ntimein)) then
         write(ioo,'(" **** time at ",e10.3," is not in range")') time
         write(16,'(" Exiting:Time at ",e10.3," is not in range")') time
         stop
      else
         do i=1,ncolumn
            hydrAtT(i) = 0.00
         enddo
      endif

      if(ntimein.eq.1) then
         hydrAtT(14)=dilution
         hydrAtT(15)=zmix
         hydrAtT(16)=percent
         do i=2, ncolumn
            if(icolumn(i).ne.0) then
               hydrAtT(i)=hydrdat(1,i)
            endif
         enddo
         do i=1,ngrid
            barjAtT(i)=barjin(i,1)
         enddo
         ntem=ntemin(ntimein)
         do ii=1, ntem
            fnem(ii)=fnemin(ii,1)
            tevm(ii)=tevmin(ii,1)
         enddo
         goto 100
      endif

      do 1 i=1,ntimein-1
        if (timein(i) .le. time .and. timein(i+1) .ge. time) then
          ip = i
          dt1 = time-timein(i)
          dt2 = timein(i+1)-time
          dt  = timein(i+1)-timein(i)
          go to 2
        endif
    1 continue
    2 continue

c do intepolations to get info at time. if the variable is zero do not interpolate

      hydrAtT(14)=dilution
      hydrAtT(15)=zmix
      hydrAtT(16)=percent
      do 3 i=2,ncolumn
         if (icolumn(i) .ne. 0) then
            xp1 = hydrdat(ip,i)
            xp2 = hydrdat(ip+1,i)
            inow= i
            if (xp1 .gt. 0.00 .and. xp2 .gt. 0.00) then
               If (column(i) .eq. 'ne ' .or.
     1              column(i) .eq. 'rho' .or.
     2             column(i) .eq. 'nt'  .or.
     3             column(i) .eq. 'dilu' ) then
                  hydrAtT(inow) = exp((log(xp1)*dt2+log(xp2)*dt1)/dt)
               else
                  hydrAtT(inow) = (xp1*dt2+xp2*dt1)/dt
               endif
            else
               hydrAtT(inow) = (xp1*dt2+xp2*dt1)/dt
            endif
         endif
 3    continue


c     temperatures are interpolated

      ntem1=ntemin(ip)
      ntem2=ntemin(ip+1)
      if(ntem1.ne.ntem2) then
         write(ioo,'(" **** Te interpolation error -AtTime-")')
         write(16,'(" Warning:Te interpolation error -AtTime-")')
      else if(ntem1.eq.1) then
         ntem=1
         fnem(1)=1.d0
         tevm(1)=teAtT
      else
         ntem=ntem1
         do ii=1, ntem
            xp1 = max(fnemin(ii,ip),1.d-30)
            xp2 = max(fnemin(ii,ip+1),1.d-30)
            fnem(ii)=exp((log(xp1)*dt2+log(xp2)*dt1)/dt)
            xp1 = max(tevmin(ii,ip),1.d-30)
            xp2 = max(tevmin(ii,ip+1),1.d-30)
            tevm(ii)= (xp1*dt2+xp2*dt1)/dt
         enddo
      endif

      if (radflag .eq. 'trfile') then
         ip = 0
         do 6 i=1,ntimtr-1
            if (rdtime(i) .le. time .and. rdtime(i+1) .ge. time) then
               ip = i
               dt1 = time-rdtime(i)
               dt2 = rdtime(i+1)-time
               dt  = rdtime(i+1)-rdtime(i)
               go to 7
            endif
 6       continue
 7       continue
         if (ip .eq. 0) then
            do 8 i=1,ngrid
               barjAtT(i) = 0.00
 8          continue
         else
            do 9 i=1,ngrid
               barjAtT(i) = (barjin(i,ip)*dt2+barjin(i,ip+1)*dt1)/dt
 9          continue
         endif
      endif

 100  fhot = 0.d0
      if(feflag .eq. 'off') return

      do i=1,ntimefe-1
         if (fetime(i) .le. time .and. fetime(i+1) .ge. time) then
           ip = i
           dt1 = time-fetime(i)
           dt2 = fetime(i+1)-time
           dt  = fetime(i+1)-fetime(i)
           go to 10
        endif
      enddo
 10   continue


      sum=0.d0
      do i=1,nebins
        fe(i) = exp(( log(feIn(i,ip))*dt2 + log(feIn(i,ip+1))*dt1 )/dt)
        sum = sum + sqrt(energy(i))*fe(i)*denergy(i)
      enddo

      if(sum.le.1.d-30)then
         print *, 'no fe()'
         fhot = 0.d0
      else
         do i=1, nebins
            fe(i)=fe(i)/sum
         enddo
         fhot1=fevfe(ip)
         fhot2=fevfe(ip+1)
         fhot = exp((log(fhot1)*dt2+log(fhot2)*dt1)/dt)
      endif


      return
      end



      subroutine fill

      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'timestuf'
      include 'popstuf'
c.... now pop array has all the populations that exist fill the
c      the total array popout

      do 1 i=1,numpop
         popout(i) = 0.d0
         ipop=indl2p(i)
         if(ipop.ne.0) then
            popout(i)=pop(ipop)
         endif
    1 continue

c.... put the existing populations in the correct population bin

      return
      end


      subroutine zbarTF(z,ddion,tev,zbar)

      implicit real*8 (a-h,o-z)

c
c
c
c    this function calculates an approximation to the thomas-fermi
c    LTE degree of ionization
c
c    z    - nuclear charge
c    tev  - temperature (ev)
c    rnen - electron number density (/cc)
c    zbar - # of free electrons per ion is returned in zbarTF
c    dion - ion number density (/cc)

c form tzero,r,tf
c
      dion = ddion
      fthrds = 4./3.
      tzero = tev/(z**fthrds)
      r = dion/z/6.022045e+23
      tf = tzero/(1.0+tzero)
c
c setup constants
c
      a1 = 0.003323467
      a2 = 0.97183224
      a3 = 9.26148e-05
      a4 = 3.1016524
      b0 = -1.762999
      b1 = 1.4317567
      b2 = 0.31546338
      c1 = -0.36666667
      c2 = 0.98333333
      alpha = 14.3139316
      beta = 0.66240046
c
c calculate a,b and c
c
      aa = (a1*(tzero**a2))+(a3*(tzero**a4))
      b = -exp(b0+(b1*tf)+(b2*(tf**7)))
      c = c2+(c1*tf)
c
c calculate q1 and thereby q
c
      q1 = aa*(r**b)
      q = (r**c)+(q1**c)
      cm = 1.0/c
      q = q**cm
c
c calculate x
c
      x = alpha*(q**beta)
c
c calculate zstar
c
      f = x/(1.0+x+(sqrt(1.0+(2.0*x))))
      rnen = f*z*dion

c... note that zbar = f*z

      zbar = f*z

      return
      end

      function zbarion(density,ismin,ismax)

      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'timestuf'
      include 'popstuf'
      dimension poptemp(numpp)

c     using population density of LSODE level index

      do 1 i=1,neq
        poptemp(i) =pop(i)
    1 continue
      do iz=0,iatomz
         popion(iz)=0.d0
      enddo
c.... now pop array has all the populations that exist find zbar

      partne=0.d0
      density=0.d0 
      do 2 i=1,neq
         ilev=indp2l(i)
         density = density+poptemp(i)
         partne = partne+poptemp(i)*ionstg(ilev)
         iso=isolev(ilev)
         popion(iso)=popion(iso)+poptemp(i)
 2    continue

      zbarion = partne/density

      pmax=0.d0
      imax=0
      do i=0, iatomz
         if(popion(i).gt.pmax)then
            pmax=popion(i)
            imax=i
         endif
      enddo
      do i=imax,0,-1
         if(popion(i).gt.1.e-10*pmax)then
            ismin=i
         endif
      enddo
      do i=imax,iatomz
         if(popion(i).gt.1.e-10*pmax)then
            ismax=i
         endif
      enddo
      return
      end

      function zbarneq(density)

      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'timestuf'
      include 'popstuf'
      dimension poptemp(numpp)

c     using population density of LSODE level index

      do 1 i=1,neq
        poptemp(i) =pop(i)
    1 continue
c.... now pop array has all the populations that exist find zbar

      partne=0.d0
      density=0.d0 
      do 2 i=1,neq
         ilev=indp2l(i)
         density = density+poptemp(i)
         partne = partne+poptemp(i)*ionstg(ilev)
 2    continue

      zbarneq = partne/density

      return
      end

      function zbarnum(density)

      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'timestuf'
      include 'popstuf'
      dimension poptemp(numpp)

c     using population density of original level index
      do 1 i=1,numpop
         poptemp(i) = popout(i)
    1 continue

      partne=0.d0
      density=0.d0 
      do 2 i=1,numpop
         density = density+poptemp(i)
         partne = partne+poptemp(i)*ionstg(i)
 2    continue

      zbarnum = partne/density

      return
      end
c
c
      subroutine subwraa(tenowio,denowio,tnextio)

      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'timestuf'
      include 'popstuf'
      include 'runstuf'
c
c___write diagnostic output
c
      poptot = 0.00

      if (ievolve .eq. 'td') then
         nmat = neq
      else
         nmat = neq-1
      endif

      do 11 i=1,neq
         b(i) = 0.00
         poptot = poptot+x(i)
         do 12 j=1,nmat
            b(i) = b(i)+a(j,i)
   12    continue
   11 continue
c
      do iz=iatomz,1,-1
         write(16,'("iz=",i3)')iz
         ibeg=indl2p(indlev(1,iz))
         izzz=max(iz-2,0)
         iend=indl2p(indlev(1,izzz))-1
         do i=ibeg,iend
            write(16,909)kname(indp2l(i)),
     &           kname(indp2l(ibeg)),kname(indp2l(iend))
            write(16,910)(a(i,j),j=ibeg,iend)
         enddo
      enddo
c     # depressed out

      write(16,901) (iz, ire(iz), iz=iatomz, 1, -1)
      write(16,904) (iz, iraut(iz), iz=iatomz, 1, -1)

c   sum of the verticals

      write(16,902) (b(i),i=1,neq)

c   total population and individual populations

      write(16,903) poptot,(pop(i),i=1,neq)


 900  format(a8,1p,10e10.2,10(/,8x,10e10.2))
 901  format(10(1x,'ire(',i2,')=',i2))
 904  format(10(1x,'ira(',i2,')=',i2))
 902  format(////,1x,"sumd",1p,10e10.2,/,(5x,10e10.2))
 903  format(/,1x,"pop ",1p,e10.2,/,(5x,10e10.2))
 909  format('a(i,j) for i=',a8,'j=',2(2x,a8))
 910  format(1p,10e10.2)

      return
      end

      subroutine loadinit(done)

      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'timestuf'
      include 'popstuf'
      include 'runstuf'

      character*120 msg
      character*20 isym(35)
      integer lsym(35)

      logical done,popgood

      done = .false.
      popgood = .false.

c... the information on the initial population file is opened on unit 20
c    the format is a level name and then a population number.
c    the level names must be in the existing level list

      do ilv=1,ntot
         popinit(ilv)=0.d0
         indipd(ilv)=1
      enddo

      open(unit=20,file=initfile,status='old')
   1  read(20,'(a)',end=100,err=200) msg
      call getsymb(msg,120,isym,nsym,lsym,35)
      if (nsym .eq. -1) then
         write(ioo,*)' **** Char string too long'
         write(16,*)' Warning:Char string too long'
         go to 300
      endif
      if (nsym .eq. 0) go to 1
      if (nsym .ne. 2) then
         write(ioo,*)' **** The init format is NAME # one per line'
         write(16,*)' Warning:The init format is NAME # one per line'
         go to 300
      endif

      read(isym(2),'(bn,f20.0)',err=201) temp
      if (temp .lt. 0.00) then
         write(ioo,*)' **** negative init pop read'
         write(16,*)' Warning:negative init pop read'
         go to 300
      else if (temp .gt. 0.00) then
         popgood = .true.
      else if (temp .eq. 0.00) then
         go to 1
      endif

      do 10 ilv=1,ntot
         if(kname(ilv).eq.isym(1))then
            popinit(ilv)=temp
            indipd(ilv)=1
            go to 1
         endif
 10   continue   
 2    continue
      write(ioo,'(" **** ",a6," is not in calculated array")') isym(1)
      write(16,'(" Warning:",a6," is not in calculated array")') isym(1)
      go to 300

  100 continue

      if (popgood) then
         done = .true.
      else
         write(ioo,*)' **** Initfile read successful, but total is 0'
         write(16,*)' Warning:Initfile read successful, but total is 0'
      endif

c   determine the smallest non-zero contribution

       if (popgood) then
         popmin = 1.00e+30
          do 110 j=1,numpop
             if (popinit(j) .gt. 0.00)  then
                popmin = min(popmin,popinit(j))
             endif
  110    continue
         popmin = popmin*1.00d-70
         do 120 j=1,numpop
            if (popinit(j) .eq. 0.00) then
               popinit(j) = popmin
            endif
  120    continue
      endif

      go to 300

  200 write(ioo,*)' **** Difficulty reading initfile'
      write(16,*)' Warning:Difficulty reading initfile'
      go to 300

  201 write(ioo,*)' **** Could not read an init pop into array'
      write(16,*)' Warning:Could not read an init pop into array'
      go to 300

  300 continue
      close(unit=20,status="keep")
      return
      end

c
c
c          *        *******  *******
c          *           *     *
c          *           *     *
c          *           *     ****
c          *           *     *
c          *           *     *
c          *******     *     *******
c
c
c
c
      subroutine lte(deninit)
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
c
      dimension dnew(0:miso)
c
c*******************
c*  begin routine  *
c*******************
      izlast=0
      do i=iatomz,1,-1
         ilv=indlev(1,i)
         if(indipd(ilv).ne.0) then
            izlast=i
            goto 10
         endif
      enddo
      do i=1, neq
         d(i)=0.d0
      enddo

 10   dnew(0)=1.d0
      ilev=indlev(1,0)
      ipop=indl2p(ilev)
      d(ipop)=1.d0

      do 1 i=1, izlast
         eionz = eipz(i)-debye(iatomz-i+1)
         expeion=exp(eionz/tev)
         ilv1=indlev(1,i)

         if(eionz.gt.1.d2*tev) then
            do k=0, i-1
               do j=1, nlev(k)
                  ilv=indlev(j,k)
                  if(indipd(ilv).ne.0) then
                     ipop=indl2p(ilv)
                     d(ipop)=0.d0
                  endif
               enddo
            enddo
            dnew(i)=1.d0
         else
            ilv1=indlev(1,i)
            ilv2=indlev(1,i-1)
            dnew(i)=denow/c1*glev(ilv1)/glev(ilv2)*dnew(i-1)*expeion
            if(dnew(i).gt.1.d30) then
               do k=0, i-1
                  do j=1, nlev(k)
                     ilv=indlev(j,k)
                     if(indipd(ilv).ne.0) then
                        ipop=indl2p(ilv)
                        d(ipop)=d(ipop)/dnew(i)
                     endif
                  enddo
               enddo
               dnew(i)=1.d0
            endif
         endif
         do 20 j=1, nlev(i)
            ilv=indlev(j,i)
            if(indipd(ilv).ne.0) then
               ipop=indl2p(ilv)
               eexct = (elev(ilv)-elev(ilv1))/tev
               d(ipop)=dnew(i)*exp(-eexct)*glev(ilv)/glev(ilv1)
            endif
 20      continue
 1    continue
 2    continue

c****************************************************
c*  fill level populations by total number density  *
c****************************************************
c
c
      sum = 0.0
      do 4 i=1,neq
         sum = sum+d(i)
    4 continue
c
      tot = deninit/sum
      do 5 i=1,neq
         pop(i) = tot*d(i)
    5 continue

      return
      end

      subroutine depress(iexit)

      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'timestuf'
      include 'popstuf'
      common /mixture/ zelse,percent,elseznum,amratio
      common /paramnow/ trnow,tinow,zbarnow,totpnow,sizenow,rhonow
      logical ipdexist
      dimension  indgrn(0:miso), ianorm(miso)

      ztot=zelse+zave

      r_debye = 743.3913*sqrt(tev/(denow*(1.+ztot)))

      if(r_lattice.gt.1.d-30)then

c  2012 July 5th ----------------- IPD model---------------
c  Contribution from Oxford University Justin Wark's group
c  Orlando Ciricosta, Sam Vinko
c  Ecker-Kroll model

         r_ion = (0.75*1.d0/(pi*(denow+totpnow)))**(1./3.)
         dencri=(0.75/pi)*(tev*evtoerg/(iatomz*charge)**2)**(3.)
         debye0=(charge**2)/(r_debye*evtoerg)
         debye1=(charge**2)/(r_ion*evtoerg)
         write(16,'("-IPD-"1p,10e10.2)')tev,denow,ztot,dencri,
     &        r_ion, r_debye, ccost, debye0,debye1,r_lattice
         
c         write(*,'("-IPD-"1p,10e10.2)')tev,denow,ztot,dencri,
c     &        r_ion, r_debye, ccost, debye0,debye1,r_lattice
         
         
      else


c  Stewart-Pyatt formalism on ion sphere radius
c      rion1 = (0.75*ztot/(pi*denow))**(1./3.)
         rion1 = (0.75*1.d0/(pi*denow))**(1./3.) ! for neutral
         r_ion = rion1
         if(ztot/denow.lt.1.d-30/totpnow) then
            rion2 = (0.75/(pi*totpnow))**(1./3.)
            r_ion=rion2
            write(16,'(" Warning: z/Ne is replaced by Ni",
     &      " in continnum lowering",1p,2e12.2)')totpnow,ztot/denow
         endif

         rfz = r_debye/r_ion
         debye0=2.16e-7/r_ion*((1.+rfz**3)**(2./3.)-rfz**2)
         debye1 = 1.45e-07/sqrt((r_ion/1.5)**2+r_debye**2)
         write(16,'("-IPD-"1p,10e10.2)')tev,denow,ztot,totpnow,
     &        r_ion, r_debye, debye0,debye1
c         write(*,'("-IPD-"1p,10e10.2)')tev,denow,ztot,totpnow,
c     &        r_ion, r_debye, debye0,debye1
         
      endif
      

c     Update energy levels
      do iso=1, iatomz
         do i=1,nlev(iso)
            ilv=indlev(i,iso)
            jlv=lvcont(ilv)
            evtoip(ilv)=ev2ip(ilv)-debye(iatomz-iso+1)
            elev(ilv)=elev(jlv)-evtoip(ilv)
         enddo
      enddo


c....   determine what states exist for this temperature and density

      ionow = 0
      indipd(indlev(1,0))=1 !for bare nuclei
      ibo(0)=1
      iaut(0)=0
      indgrn(0)=1

      do 1 i=1, iatomz 
         if(i.lt.isomin.or.i.gt.isomax) then
            ire(i)=nlbn(i)-1
            iraut(i)=nlev(i)-nlbn(i)
            do j=1,nlev(i)
               ilv=indlev(j,i)
               if (ipdexist(ilv,i).or.j.ne.1) then
                  indipd(ilv)=0
                  if(j.eq.1)ire(i)=nlbn(i)
               else
                  indipd(ilv)=1
               endif
            enddo
         else
            ire(i) = 0
            nboun(i)=nlbn(i)
            do j=1, nboun(i)
               ilv=indlev(j,i)
               indipd(ilv)=1
               if (ipdexist(ilv,i)) then
                  ire(i)=ire(i)+1
                  indipd(ilv)=0
               endif
            enddo
            ilv1=indlev(1,i)
            iraut(i)=0
            do j=nboun(i)+1, nlev(i)
               ilv=indlev(j,i)
               icont=lvcont(ilv)
               indipd(ilv)=1
               if(i.eq.isomin) then
                  if(indipd(ilv1).eq.0.or.ipdexist(ilv,i)) then 
                     iraut(i)=iraut(i)+1
                     indipd(ilv)=0
                  endif
               else
                  if(indipd(ilv1).eq.0) then ! if the ground state non-existent
                     iraut(i)=iraut(i)+1
                     indipd(ilv)=0
                  else if (ipdexist(ilv,i))then
                     iraut(i)=iraut(i)+1 
                     indipd(ilv)=0
                  else if (indipd(icont).eq.0) then 
                     iraut(i)=iraut(i)+1 ! if the continum state non-existent
                     indipd(ilv)=0
                  endif
               endif
            enddo
         endif
         ibo(i)=nboun(i)-ire(i)
         iaut(i)=nlev(i)-nboun(i)-iraut(i)
         if(ibo(i)+iaut(i).ne.0) then
            ionow = i
         endif
         indgrn(i)=indipd(indlev(1,i))
    1 continue

c check the anomalous cases 
c where the lower charge states surives while higher ones don't

 10   do i=1,isomax
         ianorm(i)=1
         if(indgrn(i).ne.0.and.indgrn(i-1).eq.0) then
            write(16,'("IPD at z ",i3,1p,4e13.2)')i,eipz(i-1),
     &           debye(iatomz-i+2),eipz(i),debye(iatomz-i+1)
            ianorm(i)=0
            write(16,'(4(10i2,1x,10i2,1x,10x))')(indgrn(j),j=1,isomax)
            write(16,'("IPD at z ",i3,1p,4e13.2)')i,eipz(i-1),
     &           debye(iatomz-i+2),eipz(i),debye(iatomz-i+1)
         endif
      enddo
      flag=0
      do i=1, isomax
         if(ianorm(i).eq.0) then
            flag=1
            ire(i)=0
            do j=1,nboun(i)
               ilv=indlev(j,i)
               ire(i)=ire(i)+1
               indipd(ilv)=0
            enddo
            iraut(i)=0
            do j=nboun(i)+1, nlev(i)
               ilv=indlev(j,i)
               iraut(i)=iraut(i)+1
               indipd(ilv)=0
            enddo
            ibo(i)=nboun(i)-ire(i)
            iaut(i)=nlev(i)-nboun(i)-iraut(i)
            indgrn(i)=indipd(indlev(1,i))
         endif
         if(ibo(i)+iaut(i).ne.0) then
            ionow = i
         endif
      enddo
      if(flag.ne.0) goto 10

      neq=0
      do i=0, iatomz
         neq=neq+ibo(i)+iaut(i)
      enddo


c setup indices
c Note that the indl2p is not in the monotonically ascending order.

      ipop=0
      do 11 i=iatomz,1,-1
         do 12 j=1, nlev(i)
            ilv=indlev(j,i)
            indl2p(ilv)=0
            if(indipd(ilv).ne.0) then
               ipop=ipop+1
               indp2l(ipop)=ilv
               indl2p(ilv)=ipop
            endif
 12      continue
 11   continue

      ipop=ipop+1
      ilv=indlev(1,0)
      indl2p(ilv)=ipop
      indp2l(ipop)=ilv


      if (ipop.ne.neq) then
         write(16,'(" Exiting:-depress- ",2i6)') ipop, neq
         iexit=1
      endif
      return
      end

      function debye(izspec)
      implicit real*8 (a-h,o-z)
      include 'mainstuf'


c     2012 July 5th ----------------- IPD model
      zz=real(izspec)
      if(r_ion.gt.0.d0) then
         if(r_lattice.gt.1.d-30) then
            debye0=(zz*charge**2)/(r_debye*evtoerg)
            debye1=(zz*charge**2)/(r_ion*evtoerg)
            if(denow.le.dencri) then
               debye = debye0
            else
               debye = r_lattice*debye1
            endif
         else
            rion=r_ion*zz**(1./3.)
            rfz = r_debye/rion
            debye0 = 2.16e-7/rion*((1.+rfz**3)**(2./3.)-rfz**2)
            debye1 = 1.45e-07/sqrt((rion/1.5)**2+r_debye**2)
            debye = zz*min(debye1,debye0)
         endif
      else
         debye=0.d0
      endif

      return
      end

      subroutine initial

c   determine the initial population for the next time step from the
c   last time step. This would be easy if the ipd did not change ; but
c   since the ipd changes then the # of states at the end of one time
c   step is not the same as the number or type of states at the beginning
c   of the next. This is due to the choice of picking the ipd from the
c   density and temperature initially.

      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'timestuf'
      include 'popstuf'

      dimension poptemp(numpp)

c   with the populations from the last time popout. first form an array,
c   pops, with 1's where states exist now and 0's where they do not exist
c.... redistribute the populations from popout into poptemp

      sum=0.d0
      do 1 i=1,numpop
         poptemp(i) = 0.00
         sum=sum+popout(i)
    1 continue

      poptemp(indlev(1,0)) =popout(indlev(1,0)) ! stripped

      ilv0=indlev(1,ionow)
      do i=1, iatomz
         ilv1=indlev(1,i)
         do j=1,nlev(i)
            ilev=indlev(j,i)
            ilon=lvcont(ilev) !continuum limit of ilev
            if(popout(ilev).gt.1.d-30) then
               if(indipd(ilev).ne.0) then
                  poptemp(ilev)=poptemp(ilev)+popout(ilev)  !existing pop
               else
                  if(indipd(ilon).ne.0) then
                     poptemp(ilon)=poptemp(ilon)+popout(ilev) !disappeared
                  else
                     !when the IPD suppresses a few ionization stages at once
                     if(indipd(ilv0).ne.0) then
                        poptemp(ilv0)=poptemp(ilv0)+popout(ilev) !disappeared
                     else
                        write(ioo,'(" **** Failure in -initial- ",
     &                       "due to IPD")')
                        write(16,'(" Warning:Failure in -initial- ",
     &                       "due to IPD")')
                     endif
                  endif
               endif
            else 
               if(indipd(ilev).ne.0.and.popout(ilon).gt.1.d-30) then
                  if(evtoip(ilev)/tev.le.10.)then
                     rlte=(denow*glev(ilev)/glev(ilon)/c1
     1                    *exp(evtoip(ilev)/tev))
                     ppp=popout(ilon)*rlte
                     if(ppp.lt.poptemp(ilon)) then
                        poptemp(ilon)=poptemp(ilon)-ppp
                        poptemp(ilev)=poptemp(ilev)+ppp
                     endif
                  endif
               endif
            endif
         enddo
      enddo

c    now place the array poptemp into the array pop for use in lsode
      psum=0.d0
      do 21 ipop=1,neq
         ilev=indp2l(ipop)
         pop(ipop)=poptemp(ilev)
         psum=psum+pop(ipop)
 21   continue

      if(abs(psum-sum).gt.max(1.d0,1.e-4*sum))then
         write(ioo,'(" Warning:Check pop sum in -initial- ",
     &        "due to IPD",1p,e12.2)')psum, sum
         write(16,'(" Warning:Check pop sum in -initial- ",
     &        "due to IPD",1p,e12.2)')psum, sum
         do 22 ipop=1,neq
            pop(ipop)=pop(ipop)*sum/psum
 22      continue
      endif

      dinow=psum

      return
      end

      subroutine coronal(deninit,zbar)

      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      external debye

      dimension d(100)

c....  iatomz+1 are the number of ion stages IP depression is ignored

      iazlow=0
      do 10 i=1, iatomz
         ebyt = (eipz(i)-debye(iatomz-i+1))/tev
         if (ebyt .lt.0.d0) then
            d(i) = 0.00
         else if (ebyt .lt. 20.00) then
            iazlow = i
            d(i) = 1.00
            go to 11
         else
            d(i) = 0.00
         endif
   10 continue
   11 continue


      mxatom=iatomz
      do i=1, iatomz
         ebyt = (eipz(i)-debye(iatomz-i+1))/tev
         if (ebyt .lt.0.d0) then
            mxatom=i-1
            goto 12
         endif
      enddo
 12   continue

c.... define the array d in which each element is the ratio of the population
c     of ion stage i to the fully stripped stage iatomz+1 (iazp1)

      if(iazlow.eq.0) then
         sum=0.d0
         goto 20
      endif
      sum = d(iazlow)
      do 1 i=iazlow+1,mxatom
         if(sxx(1,1,i).gt.0.d0) then
            d(i) = (alfxx(1,1,i)+dierec(i))/sxx(1,1,i)*d(i-1)
            sum = sum+d(i)
         else 
            d(i) = 0.d0
         endif
    1 continue

c... Sum is the ratio of the total density to pop(iazp1).
c     Find pop(iazp1)/deninit

      if(sum.gt.1.d-70)popiazp1 = 1./sum
         
c.... now sum up the electron density from each ion stage, ie, z*N

      sum = 0.d0
      do 2 i=iazlow, mxatom
         sum = sum+d(i)*popiazp1*real(iatomz-i+1)
    2 continue

 20   if(sum.eq.0.d0) then
         write(ioo,'(" $$$ check -coronal-; zbar=0.01")')
         write(16,'(" Warning:Check  -coronal- ; zbar=0.01")')
         sum=1.d-2
      endif
      zbar = sum

      return
      end
c
      logical function ipdexist (ind,nion)
c
c
c___  test for non-existant level in -nion- stage.
c
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'popstuf'
c
c___test for non-existance.
c
      ipdexist = evtoip(ind).le.0.d0

      return
      end

