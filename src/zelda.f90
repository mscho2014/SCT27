!       04/28/03 (in Berkeley):

!----------------------------------------------------------------------

!       11/08/02 (in Berkeley):

!         (1) Hyun's stuff included in DCOEFF & zhk_stuff.f90 replaces flystuff.f90
!         (2) problem with electron non-conservation fixed

!         (3) need to fix x2 jump in Ne after 1st time step
!         (4) need to make continuous approximation for small energy loss processes
!         (5) need to include ionization secondaries & 3-body recombination as sources/sinks
!         (6) use EnDegrad for fast electrons/beam component
!         (7) electron-electron collisions

!----------------------------------------------------------------------

!     module files:

!       'zmodules.f95'      ! include file containing modules
!       'zhk_stuff.f95'     ! include file containing commons from FLY - 11/08/02 from Hyun

!------------------------------------------------

      module MonteCarloSampling

        use EEDF_DIMENSIONS

        integer, parameter :: MaxNumberOfLevels = 3500

        type MCstuff

          logical :: MonteCarloFlag

          double precision, dimension (MaxNumberOfLevels,NPTSF) :: Rate

          double precision, dimension (MaxNumberOfLevels,NPTSF) :: DeltaE

          integer :: NEnergies

          double precision :: dE

          double precision, dimension (NPTSF) :: TotalRate, Energy

          integer, dimension (NPTSF) :: NumberOfLevels

        end type MCstuff

        type (MCstuff) :: MC

      end module MonteCarloSampling

!------------------------------------------------

      subroutine ZELDA(mcall)

      use FILENAMES_and_UNITS

      use EFLAGS


      implicit none
      
      integer ncall, mcall


!------------------------------------------------

      if (mcall == 0) open (unit=16, file='screen.out', status='unknown')
      ! zelda.out checks the results in zelda.f90 code  
!------------------------------------------------

      ncall = ncall + 1                 ! 1st pass

!      call STD_FILE_OPEN               ! open/create standard i/o files
!      call AUX_FILE_OPEN               ! create files for other output

!      write (crt_out,*) 'ZELDA call number ', ncall
!      write (crt_out,*) 'ZELDA call number ', mcall
      
      call SOLVE_BOLTZMANNS_EQN(mcall)

      return

      end subroutine ZELDA
      
!*--------------------------------------------------------------------72

      subroutine SOLVE_BOLTZMANNS_EQN(mcall)


      use EEDF_DIMENSIONS

      use FLY_PARAMETERS                                      ! FLY stuff

      use EFLAGS

      use FofE_STUFF                                          ! FLY stuff

      use POP_STUFF                                           ! FLY stuff

      use ATOMIC_and_PLASMA_STUFF                             ! FLY stuff

      use TIMING_STUFF                                        ! FLY stuff

      implicit none

      integer i, nmax, npts, ncall, mcall

      double precision:: fnorm, ave_engy

      double precision, dimension (NPTSF)    :: bfd

      double precision, dimension (0:NPTSF1) :: u
      double precision, dimension (0:NPTSF1) :: f

      double precision, dimension (0:NPTSF1) :: e0
      double precision, dimension (0:NPTSF1) :: f0

      integer NFdiffCall

!------------------------------------------------

!   Solution of linear equations by decomposition (DECOMP)
!   and backsolve (BKSOLVE)
!   Vector of solutions is returned in 'bfd'

      nmax = nebinit
      npts = nebins

      CALL FDIFF (bfd, npts, nmax, ncall)  ! key routine for solving Boltzmann's eqn

      DO i = 1, npts
         f(i) = bfd(i)
      ENDDO

!       value at zero energy:
      u(0) = 0.0
      f(0) = f(1) + 0.5 * ( f(1) - f(2) )       ! linear interpolation
!      f(0) = f(1) * dsqrt( f(1) / f(2) )        ! exponential  "

!------------------------------------------------

!   Normalize the distribution function and
!   calculate the mean energy

      CALL INTGRL (fnorm, ave_engy, f(0),npts, nmax)

      f(0) = fnorm * f(0)

!     The following sends stuff back to FlyChk:
      do i=1, npts
        fe1(i) = f(i)  ! fe1(i) is the final result of zelda.f90  
      enddo

      tev = dble(2. * ave_engy / 3.)

      return

      end subroutine SOLVE_BOLTZMANNS_EQN

!*--------------------------------------------------------------------72

      subroutine FDIFF (b, n, nmax, ncall)

!   Main routine for matrix calculation of distribution function

!   With capability of including electron-electron and free-free
!   processes

!   note that the quantity b(e) is sqrt(e)*f0(e)*Ne and has
!   units electrons/cc/eV

!------------------------------------------------

      use EEDF_DIMENSIONS

      use FILENAMES_and_UNITS

      use EFLAGS

      use ATOMIC_and_PLASMA_STUFF  ! FLY stuff

      use TIMING_STUFF             ! FLY stuff

      use POP_STUFF                ! FLY stuff

      use FofE_STUFF


!------------------------------------------------
      implicit none

      integer nmax
      integer mcall 

      double precision, dimension (0:NPTSF1) :: u
      double precision, dimension (0:NPTSF1) :: f

      double precision, dimension (0:NPTSF1) :: e0
      double precision, dimension (0:NPTSF1) :: f0

      double precision t1(NPTSF), t2(NPTSF),t3(NPTSF)
      double precision usqrt(NPTSF), Source(NPTSF), Sink(NPTSF)
      double precision S_mon(NPTSF) ! source monitor  

      double precision b(nmax)
      double precision bold(NPTSF), bolder(NPTSF)
      double precision :: du

      double precision :: cf0(NPTSF,NPTSF)

      integer ik
      DIMENSION ik(NPTSF)

!       set up for dynamic allocation of memory for large arrays:
!     dimension cf(:,:), aee(:,:)
!     allocatable :: cf, aee

      double precision cf(NPTSF,NPTSF), aee(NPTSF,NPTSF)

      integer ismax
      integer ierr_all
      integer NFdiffCall, ncall
      integer i, j, n, ifmax
      integer istep, ipf, iprint

      integer k ! mscho 200909

      double precision :: dt, eps
      double precision :: tlower, tupper,time
      double precision :: sum1, sum2, ebsum, bsum, bmax
      double precision :: edenst0, edenst, dnedt
      double precision :: bb, bfr, bfrmax, bfrave, ebfrmax
      double precision :: fnorm, ave_engy
      double precision :: cfmax, gdens
      double precision :: umax, edens, tel, telold
      double precision :: Src3b(NPTSF)
!------------------------------------------------

!      nmax = NPTSF

!------------------------------------------------

!       dynamic allocation of arrays:

!       dynamically allocate space for the coefficient array "cf":
!       (this now has to be passed as an argument)
      ierr_all = 0
!     allocate ( stat = ierr_all, cf(n,n) )
      if (ierr_all /= 0) then
        print *, 'ERROR in dynamic memory allocation: memory exhausted'
        write (std_out,*) 'ERROR in dynamic memory allocation for collision matrix: memory exhausted'
        err_flag = .TRUE.
        return
      endif

!       allocate space for e-e collision array if e-e collision flag is set:
      ierr_all = 0
      if (eeflag) then
!       allocate ( stat = ierr_all, aee(n,n) )
        if (ierr_all /= 0) then
          print *, 'ERROR in dynamic memory allocation: memory exhausted'
          write (std_out,*) 'ERROR in dynamic memory allocation for e-e collisions: memory exhausted'
          write (err_out,*) 'ERROR in dynamic memory allocation for e-e collisions: memory exhausted'
          err_flag = .TRUE.
          return
        endif
      endif

!------------------------------------------------

      du   = denergy(1)
      umax = energy(n)


! hyun's correction to reduce dt and get from FLY

      dt    = dtfly
      ismax = isfly  !isfly = 20 in zchkmn.f 
      iprint = 1
      eps    = 1.e-16

      write ( std_out,9001 ) n, umax, dt
      u = 0.d0

      DO i = 1 , n
          u(i)     = energy(i)
          usqrt(i) = sqrt(u(i))
      ENDDO

      edens = denow
      tel   = tev
      gdens = dinow  ! sum of popout(i)

      write (*,'(" Zelda input: Te,Ne,Ni,<Z> = ",f12.6,1p,3e12.4)') tev,denow,gdens,zave
      write (*,'(" Zelda input: du,umax,dt,t, = ",1p,4e12.4)')du,umax,real(dt),real(total_time+tfirst)

!hyun's commenting out FALSE and making it TRUE
!      eeflag = .FALSE.
       eeflag = .TRUE.
       eiflag = .FALSE.
       !cf = 0
!       compute matrix for electron-electron collisions:
      IF (eeflag) then

        CALL EECOLL1 (n, nmax, du, cf, aee)
        if (err_flag) return

      endif
   
      CALL DCOEFF (du, n, nmax, 0, ncall, cf)

      fflag = .TRUE.
      b = 0.d0

      IF (fflag) then
         do i=1, n
            b(i) = fe0(i) * sqrt(energy(i)) * edens
            !b(i) = fe0(i) * edens
         enddo
         if (err_flag) return
      ELSE
!        initial f(e) is Maxwell-Boltzmann:
         call MAXWELL_BOLTZMANN (b, n)
      ENDIF

      bsum = SUM (b)

      do i=1, n
         bold(i) = b(i) / sqrt(energy(i)) / bsum
         bolder(i) = bold(i)
      enddo
      bmax = maxval(bolder)

!      set up collisional matrix for decomposition:
      cfmax=0.d0
      do i=1, n
         if(dabs(cf(i,i)).ge.cfmax)then
            cfmax = dabs(cf(i,i))
            ifmax = i
         endif
      enddo
      write(*,'(" Zelda: 1 > dt * cfmax =  ",1p,2e12.2,i5)')dt*gdens*cfmax,cfmax*gdens,ifmax

! hyun's check--> OK
      cf = -(dt * gdens) * cf
      DO i = 1 , n
         cf(i,i) = 1.0d0 + cf(i,i)
      ENDDO


      CALL DECOMP (n, nmax, ik, cf)

      if (err_flag) return

!----------------------------------------------------------------------

       ebflag  = .TRUE.   ! source function comes from FlyChk
!      ebflag  = .FALSE.   ! source function comes from FlyChk

!----------------------------------------------------------------------

      NFdiffCall = NFdiffCall + 1


!------------------------------------------------
!---------- BEGINNING OF TIME LOOP --------------
!------------------------------------------------
      
!     external source of electrons included explicitly:

      IF (ebflag) then
         CALL ESOURCE (Source, b, ncall, n, nmax)
         if (err_flag) return

         bsum=0.d0
         do i=1, n
            bsum=bsum+zqs(i)-zqr(i)*fe(i)*denow
         enddo

         print *, 'extra source-ele', bsum

      ENDIF

      edenst = edens
      time   = 0.0d0
      tlower = total_time
!     BUG FIX: status='old' causes runtime error on first call (file doesn't exist yet).
      open(unit=92,file='source_monitor',status='unknown',position='append')
      write(92,*)'mcall=',mcall

      open(unit=93,file='zelda.out',status='old',position='append')


      Time_Loop: DO istep = 1 , ismax            ! begin time loop 
        !mscho 200326 : from now, ismax = 20 (=isfly in flychk code) 
         ipf    = mod(istep,iprint) !mscho 200326 : (1) mod(A,P) : computes the remainder of the division of A by P 
         time   = time + dble(dt) 
         total_time = total_time + dble(dt)
         telold = tel

! hyun's explicit take on Source and sinks

         if(ebflag)then

            b = b + dt * Source(1:n)
            do i=1, n
               Sink(i)=zqr(i)/usqrt(i)/denergy(i)*(dt*b(i))
               S_mon(i) = dt*Source(i)/usqrt(i)/denow !for source_monitor 
               if(istep==ismax)then
                write(92,'(I4," ",E12.4," ",E12.4)')i,S_mon(i),Sink(i)
               endif 
            enddo
            b = b - Sink(1:n)

         endif

!        backsolving matrix equation:
         CALL BKSOLVE (n, nmax, b, ik, cf)

         IF (eeflag) THEN
!           Quasi-implicit treatment of electron-electron collisions
           sum1 = DOT_PRODUCT (b, u(1:n))
           sum2 = SUM (b)
           tel = .66667 * sum1 / sum2
           tupper = tlower + dt
           call ee_Collisions (tlower, tupper, b) ! mscho 200326 - subroutine structure : ee_Collision > odeint > rkqs > rkck > derivs
           tlower = tupper
         ENDIF


         ebsum   = DOT_PRODUCT (b, u(1:n))
         bsum    = SUM (b)

!         new temperature
         tel     = .66667 * ebsum / bsum

!         new density
          edenst0 = edenst
          edenst  = du * bsum
          edens   = edenst
!         new density increase
          dnedt   = (edenst - edenst0) / dt

          write (93,'("Ncall,istep,t,Te,Ne,ave & max %chg & 
                  @ E = ",2i5,1p,e12.4,0p,f10.2,1p,3e12.4,0p,f10.2)') &
                  NFdiffCall, istep,total_time,tel,edenst, &
                  100.*real(bfrave),100.*real(bfrmax),ebfrmax  


      enddo Time_Loop                            ! end time loop

      close(92)
      close(93)

!-------------- END OF TIME LOOP ----------------
!------------------------------------------------

      ebsum   = DOT_PRODUCT (b, u(1:n))
      bsum    = SUM (b)
      
      tel     = ebsum / bsum * 2./ 3.


      write (*,'("*Ncall,t,Te,Ne,ave&max%chg@E",i5,1p,e12.4,0p,f14.6,1p,5e12.5)') &
            NFdiffCall,total_time,tel,edenst,100.*real(bfrave),100.*real(bfrmax),real(ebfrmax)

      write (std_out,'("Ncall,istep,t,Te,Ne = ",2i5,1p,e12.4,0p,f10.2,1p,e12.4)') &
            NFdiffCall, istep,total_time,tel,edenst

      bsum    = SUM (b)*denergy(1)

! making b-->f 

      b(1:n) = b(1:n) / usqrt(1:n) 
      ! mscho 200326 This means b(1:n) is the normalized distribution function (not DOS*distribution) 

      RETURN

 9001 FORMAT(/23x,'CALCULATIONAL PARAMETERS'/   &
               5x,'No. of Grid Points',         &
               5x,'Upper Integration Limit',    &
               5x,'Time Step'/                  &
              37x,'(eV)',17x,'(sec)'/           &
              11x,i3,20x,f8.3,13x,1p,e10.3)
 9002 FORMAT(/28x,'GAS PARAMETERS'/                 &
              17x,'Pressure',                       &
               5x,'Temperature',                    &
               5x,'Density'/                        &
              19x,'(atm)',10x,'(K)',10x,'(/cc)'/    &
              15x,1pe10.3,6x,0pf7.1,6x,1p,e10.3)
 9003 FORMAT(/25x,'ELECTRON PARAMETERS'/                &
              19x,'Density',5x,'Initial Temperature'/   &
              20x,'(/cc)',13x,'(eV)'/                   &
              18x,1p,e10.3,7x,e10.3)
 9004 FORMAT(/28x,'PHOTON PROCESS'/                           &
              11x,'Photon Energy',                            &
               5x,'Power Density',                            &
               5x,'Photon Flux'/                              &
              16x,'(eV)',11x,'(W/sq.cm)',7x,'(/sq.cm/sec)'/   &
              12x,1p,e10.3,8x,e10.3,7x,e10.3)
 9005 FORMAT(/16x,'ELECTRON-ELECTRON COLLISIONS INCLUDED')
 9006 FORMAT(/19x,'ELECTRON-ION COLLISIONS INCLUDED')
 9007 FORMAT(/24x,'RELAXATION CALCULATION'/                  &
               1x,'Cycle No.',5x,'Time',                     &
               5x,'Max. Fract. Change',15x,'Electron'/       &
              15x,'(sec)',30x,'Temperature',3x,'Density'/    &
              53x,'(eV)',8x,'(/cc)')
 9008 FORMAT(3x,i3,6x,1p,e10.3,6x,e10.3,12x,e10.3,3x,e10.3,3x,e10.3)
 9009 FORMAT(/21x,'DC ELECTRIC FIELD PARAMETERS'/             &
              22x,'E/N',19x,'E/P'/                            &
              19x,'(V-sq.cm)',12x,'(V/cm/Torr)')
 9010 FORMAT(18x,1p,e10.3,12x,e10.3)
 9011 format (3x,'E/N = ',1pe10.3/)
 9012 format (/1p,e12.4,5x,1p,e12.4/ (1x,1p,8e9.2))
 9013 FORMAT(/21x,'AC ELECTRIC FIELD PARAMETERS'/             &
              9x,'E/N',16x,'E/P',12x,'Freq',13x,'Phase'/      &
              6x,'(V-sq.cm)',9x,'(V/cm/Torr)',8x,'(Hz)',11x,'(Radians)')
 9014 FORMAT(5x,1p,e10.3,9x,1p,e10.3,6x,1p,e10.3,7x,1p,e10.3)
 9015 format (/27x,'TABULATED E/N(t)')

!*-----------------------------------------------
      CONTAINS
!*-----------------------------------------------

        subroutine MAXWELL_BOLTZMANN (b, n)

        implicit none

        integer n
        double precision :: b 
        dimension b(n)

        double precision :: b0, bexp

        b0   = 2. / (sqrt(Pi) * tel ** 1.5) * edens
        bexp = exp(-du / tel)
        b(1) = b0 * exp(-u(1) / tel)
        DO i = 2 , n
          b(i) = b0 * exp(-u(i) / tel)
        ENDDO
        DO i = 1 , n
          b(i) = usqrt(i) * b(i)
        ENDDO

        return
        end subroutine MAXWELL_BOLTZMANN

      END subroutine FDIFF


!*--------------------------------------------------------------------72

     subroutine ESOURCE (qs, b, ncall, n, ndim)

!    this routine can be used to provide an external source of electrons
!   (e.g., an electron beam); the units of q are /cc/eV/sec.
!    we read qs from FLYCHK variable zqs
!!!   this routine can be used to provide an external source of electrons
!!!   (e.g., an electron beam); the units of q are /cc/eV/sec.
!!

      use POP_STUFF

      use FILENAMES_and_UNITS

      use EFLAGS

      use FofE_STUFF

      use ATOMIC_and_PLASMA_STUFF

      implicit none

      integer ndim, ncall
      integer i, n
      integer itran, ietran, iso, llo, lup, ind1, ind2, ilo, iup
      double precision qs(ndim), b(ndim), ysum
      double precision Valfxx, Vbstmxx,Vbcontxx,autoxx,Vcapexx
      double precision N_upper, N_lower

!   qs(i) = no. of electrons/cm**3/sec/ev
!   n     = no. of grid points
!    the sink term should be in DCOEFF
!    qs(i) = zqs(i)/denergy(i) - zqr(i)/denergy(i)*fe(i)

      qs = 0.
      Energy_Loop: do i = 1, n
         qs(i) = zqs(i)/denergy(i)
      enddo Energy_Loop

      ysum=0.d0
      Transition: do itran=1, ntran

         if(ltype(itran).le.1)cycle Transition
         llo=llower(itran)
         lup=lupper(itran)
         ind1=indl2p(llo)
         ind2=indl2p(lup)
         if(ind1*ind2.eq.0) cycle Transition
         N_upper=popout(lup)
         if(N_upper.lt.1.e-30)cycle Transition
         N_lower=popout(llo)
         if(N_lower.lt.1.e-30)cycle Transition

         ilo=levl(llo)
         iup=levl(lup)
         iso=isolev(llo)

         if(ltype(itran).eq.2)then
            ysum=ysum-N_upper*Valfxx(ilo,iup,iso)*denow

            if(radflag.ne.'off') then
               ysum=ysum+N_lower*Vbcontxx(ilo,iup,iso)  &
                    -N_upper*Vbstmxx(ilo,iup,iso)*denow 
                        
            endif
         else if(ltype(itran).eq.3) then
            ysum=ysum+N_lower*autoxx(ilo,iup,iso) &
                     -N_upper*Vcapexx(ilo,iup,iso)*denow
         endif

      enddo Transition

      print *, 'extra source-ion', ysum

      return
      end subroutine ESOURCE

!*--------------------------------------------------------------------72
!               Monte Carlo Sampling of Energy Loss Time
!*--------------------------------------------------------------------72

      subroutine MonteCarlo

      use MonteCarloSampling

      use ATOMIC_and_PLASMA_STUFF

      implicit none

      integer iE, i
      integer nMC, NMCsamples, Ncoll

      double precision :: dt, t, eEnergy, CollProb, rn, cft
      double precision :: RANF

      MC%Rate = MC%Rate * 5.93094e-9 * dinow
      MC%TotalRate = MC%TotalRate * 5.93094e-9 * dinow

      iE = MC%NEnergies
      do i=1, MC%NumberOfLevels(iE)
       write (31,*) i,MC%Rate(i,iE), MC%DeltaE(i,iE)
      enddo

      write (*,*) 'dinow = ',dinow

      write (30,*) 'N Energies = ',MC%NEnergies

      do i=1, MC%NEnergies
        write (30,*) 'i,E,N Levels, Tot Rate = ',i,MC%Energy(i),MC%NumberOfLevels(i),MC%TotalRate(i)
      enddo

      NMCsamples = 100; dt = 2.e-18

      MC_loop: do nMC=1, NMCsamples

        t = 0.; eEnergy = MC%Energy(MC%NEnergies); iE = MC%NEnergies; NColl = 0

        TimeLoop: do while (.TRUE.)

          t = t + dt

          eEnergy = MC%Energy(iE)

          CollProb = 1. - exp(-MC%TotalRate(iE) * dt)

          if (CollProb == 0.) EXIT TimeLoop

          rn = RANF()

!          write (31,*) 'CollProb, rn = ',CollProb,rn

          if (rn > CollProb) then
            CYCLE TimeLoop
          endif

          cft = MC%TotalRate(iE)

          TransitionLoop: do i=1, MC%NumberOfLevels(iE)

            cft = cft - MC%Rate(i,iE)

            if (rn > cft * CollProb / MC%TotalRate(iE)) then

              write (31,*) 'i,iE,CollProb,cft,MC%Rate(i,iE),MC%DeltaE(i,iE) = ', &
                            i,iE,CollProb,cft,MC%Rate(i,iE),MC%DeltaE(i,iE)

              eEnergy = eEnergy - MC%DeltaE(i,iE)

              iE      = eEnergy / MC%dE

              NColl = NColl + 1

              write (33,*) t,eEnergy

              EXIT TransitionLoop

            endif

          enddo TransitionLoop

          if (iE < 1) EXIT TimeLoop

        enddo TimeLoop

        write (*,*)  't, NColl = ',t,NColl
        write (31,*) 't, NColl = ',t,NColl
        write (32,*)  t, NColl

      enddo MC_loop

      STOP 'MonteCarlo'

      return
      end subroutine MonteCarlo

!---------------------------------------------------------------------72

!      function RANF()
      double precision function RANF()

      implicit real*8 (A-H,O-Z)

!       Uniform random number generator from
!       Press, et al. "Numerical Recipes"

      DIMENSION R(97)
      PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=1.d0/m1)
      PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=1.d0/m2)
      PARAMETER (M3=243000,IA3=4561,IC3=51349)
      DATA IFF /0/
      DATA R/97*0./
      DATA IX1,IX2,IX3/0,0,0/
      DATA IDUM/-1/
!      save r,idum,iff,ix1,ix2,ix3,j

      IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
        IFF=1
        IX1=MOD(IC1-IDUM,M1)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX2=MOD(IX1,M2)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX3=MOD(IX1,M3)
        DO 11 J=1,97
          IX1=MOD(IA1*IX1+IC1,M1)
          IX2=MOD(IA2*IX2+IC2,M2)
          R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
11      CONTINUE
        IDUM=1
      ENDIF

      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      IX3=MOD(IA3*IX3+IC3,M3)
      J=1+(97*IX3)/M3
      IF(J.GT.97.OR.J.LT.1) STOP 'ranf error'
      ranf=R(J)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      RETURN
      END

!*--------------------------------------------------------------------72
!*--------------------------------------------------------------------72

      subroutine DECOMP (n, ndim, ip, cf)

!   matrix decomposition

!   input:
!         n    = order of matrix in a
!         ndim = dimension of array a
!         cf   = array containing matrix to be decomposed

!   output:
!         cf        = decomposed matrix
!         ip(k), k<n index of the k-th pivot row
!         ip(n)     = (-1)**(number of interchanges) or 0
!         ip(n)     = 0 if a is singular
!         determ(a) = ip(n) * a(1,1) * a(2,2) * ... * a(n,n)

      implicit none
      integer ndim
      integer n, i, j, k, l, kp1, ip
      dimension ip(ndim)

      double precision :: cf(ndim,ndim), t


      ip(n) = 1

      Main_Loop: do k = 1 , n

          if (k == n) call ERROR_CHECK (k, n, ndim, cf)
          kp1 = k + 1
          l   = k
          do i = kp1 , n
              if (dabs(cf(i,k)) > dabs(cf(l,k))) l = i
          enddo
          ip(k) = l
          if (l /= k) ip(n) = -ip(n)
          t       = cf(l,k)
          cf(l,k) = cf(k,k)
          cf(k,k) = t
          if (t == 0.) call ERROR_CHECK (k, n, ndim, cf)
          do i = kp1 , n
              cf(i,k) = -cf(i,k) / t
          enddo
          do j = kp1 , n
              t       = cf(l,j)
              cf(l,j) = cf(k,j)
              cf(k,j) = t
              if (t /= 0.) then
                do i = kp1 , n
                  cf(i,j) = cf(i,j) + cf(i,k) * t
                enddo
              endif
          enddo

      enddo Main_Loop

      return

      CONTAINS

        subroutine ERROR_CHECK (k, n, ndim, cf)

          use EFLAGS

          implicit none
          integer k, n, ndim, ip(ndim)
          double precision cf(ndim,ndim)

          if (cf(k,k) == 0.) then
            ip(n) = 0
!            call PRNT_ERR ('DECOMP  ', 'Matrix of coeffs is singular! ')
!            print *, 'DECOMP  ', 'Matrix of coeffs is singular! '
            err_flag = .TRUE.
            return
          endif

          return
        end subroutine ERROR_CHECK

      end subroutine DECOMP

!*--------------------------------------------------------------------72

      subroutine BKSOLVE (n, ndim, b, ip, cf)

!   backsolving of matrix decomposed by DECOMP

        implicit none

        integer n, ndim, ip
        integer i, k, l, kb, kp1, km1, nm1
        dimension ip(ndim)

        double precision :: b(ndim), cf(ndim,ndim), t

        if (n == 1) then
           b(1) = b(1) / cf(1,1)
           return
        endif

        nm1 = n - 1
        do k = 1 , nm1
           kp1  = k + 1
           l    = ip(k)
           t    = b(l)
           b(l) = b(k)
           b(k) = t
           do i = kp1 , n
              b(i) = b(i) + cf(i,k) * t
           enddo
        enddo

        do kb = 1 , nm1
           km1  = n - kb
           k    = km1 + 1
           b(k) = b(k) / cf(k,k)
           t    = -b(k)
           do i = 1 , km1
              b(i) = b(i) + cf(i,k) * t
           enddo
        enddo

        b(1) = b(1) / cf(1,1)

        return
      end subroutine BKSOLVE

!*--------------------------------------------------------------------72

      subroutine EECOLL1 (n, ndim, du, cf, aee)
      implicit none

!   computes the matrix used in electron-electron calculations

      integer n, ndim
      integer i, k, j, l 

      double precision cf(ndim,ndim), aee(ndim,ndim)
      double precision :: du
      double precision :: duby2, upk,  e1, e2, e3, e4, el

      double precision :: uplus, elk, aprime
      double precision :: z, dum1

!--------------------------------------------------------
!       original Rockwood finite difference scheme:
      uplus(z)    = dum1 + .25 / (z + duby2)
      elk(k,l)    = sqrt(float((l - 1) * (k + 1)) / float(l * k))
      aprime(k,l) = sqrt(cf(k,l) * cf(l - 1,k + 1) * elk(k,l))

!--------------------------------------------------------
!      if (n > NPTSF) then
!        call PRNT_ERR ('EECOLL1 ', 'n > NPTSF; increase dimension')
!        print *, 'EECOLL1 ', 'n > NPTSF; increase dimension'
!        err_flag = .TRUE.
!        return
!      endif

! BUG FIX: duby2 was set to 0 (commented out 0.5*du), which deviates from
!     the original Rockwood (1973) finite-difference scheme. Restored to 0.5*du.
!     With duby2=0, the uplus() function at low energies was slightly off.
      duby2 = 0.5d0 * du
      dum1  = 1.0 / du

      cf = 0.; aee = 0.

      e1 = 0.0
      e2 = du
      !write(*,*) 'mscho ****: du,e1,e2 plot : 1177 line in zelda.f90' 
      !write(*,*) du,e1,e2

      Outer_Loop: do k = 1 , n-1
          e1  = e1 + du
          e2  = e2 + du
          e3  = 1.0 / sqrt(e1)
          e4  = 1.0 / sqrt(e2)
          upk = uplus(e1)
          el  = du
          Inner_Loop: do l = 2 , n
              el = el + du

              if (k > l) then
                cf(k,l) = (e4 + e3) * (el * upk - 0.75)
                CYCLE
              endif

              if (k == l) then
                cf(k,l) = (e4 + e3) * (el * upk - 0.75) + e1 * upk / sqrt(el)
                CYCLE
              endif

              if (k < l) then
                if (k >= (l - 1)) then
                  cf(k,l) = e4 * (el * upk - 0.75) + (e2 + e1) * upk / sqrt(el)
                else
                  cf(k,l) = (e2 + e1) * upk / sqrt(el)
                endif
              endif

          enddo Inner_Loop
      enddo Outer_Loop

      do k = 1 , n-1
          do l = 2 , n
              aee(k,l) = aprime(k,l)
          enddo
      enddo

      return
      end subroutine EECOLL1

!*--------------------------------------------------------------------72

      subroutine INTERPOL (xi, yi, x, y, n, err_flag)

!       for general linear interpolation of y(x)
!       assume that x(i+1) > x(i) and that x(1) <= xi <= x(n)


      use FILENAMES_and_UNITS

      implicit none

      logical err_flag

      integer n, i, i1, i2
      double precision:: x(n), y(n)
      double precision :: xi, yi

!       check the limits:
      if (xi < x(1) .or. xi > x(n)) then
        yi = 0.
        print *, 'WARNING: subroutine INTERPOL'
        print *, ' argument outside of range of data'
        write (err_out,*) 'WARNING: subroutine INTERPOL argument'
        !outside of range of data'
        err_flag = .TRUE.
        return
      endif

!       perform the linear interpolation:
      do i=1, n
        if (x(i) > xi) then
          i2 = i
          i1 = i - 1
          EXIT
        endif
      enddo

      yi = y(i1) + (y(i2) - y(i1)) * (xi - x(i1)) / (x(i2) - x(i1))

      return
      end subroutine INTERPOL

!*--------------------------------------------------------------------72


      subroutine INTGRL (ynorm, ymean, y, n, ndim)

!   Calculates normalization constant for distribution function
!   using extended midpoint rule with first three terms modified
!   to start a staggered 3-point formula for the initial square
!   root behavior. Note that the y(i) are in
!   backwards order to sum over the small values first.
!   Also computes mean energy of distribution

!----------------------------------------------------------------------
      use FofE_STUFF, ONLY: energy, denergy

      implicit none

      integer:: n, ndim, k
      double precision:: ynorm, ymean,s1, s3
      double precision:: y(0:ndim)

!----------------------------------------------------------------------

      s1    = 0.0
      s3    = 0.0

      DO k = 1, n
          s1  = s1 + denergy(k)*y(k) !*sqrt(energy(k))
          s3  = s3 + denergy(k)*y(k)*(energy(k))**1.5
      ENDDO

      ynorm = 1.0 / s1
      ymean = s3 * ynorm

      DO k = 1 , n
          y(k) = ynorm * y(k)
      ENDDO

      RETURN
      END subroutine INTGRL

!*--------------------------------------------------------------------72

      subroutine PHOTON (n, ndim, cf)

!   computes the terms in coefficient matrix for free-free
!   processes

      use RUN_PARAMETERS

      use FofE_STUFF, ONLY: energy, denergy

      implicit none

      integer ndim, n, kep, i, k, imk, ipk
      double precision ::  cf(ndim,ndim)
      double precision :: Q0, const1, const2
      double precision :: akappa, ee, ep
      double precision :: du,  fn
      double precision :: z, qme, qm1, qmemk, qm2, ak

      data Q0 / 1.e-16 /
      data const1 /3.0664e-32/, const2 /5.930e+7/

!--------------------------------------------------------
      akappa(ee,ep,qm1) = (const1 / ep ** 2) * const2 * sqrt(ee + ep) * ( (ee + .5 * ep) / ep) * qm1

!--------------------------------------------------------

      du  = denergy(1)
      kep = NINT(ephoton / du)
      fn  = photflux * Q0

      do i = 1 , n
          z       = i * du - 0.5*du
          call QMOMNTM (z + .5 * ephoton,qme,qm2)
          ak      = akappa(z,ephoton,qme)
          cf(i,i) = cf(i,i) - fn * ak

          if (z > ephoton) then
            call QMOMNTM (z - .5 * ephoton,qmemk,qm2)
            cf(i,i) = cf(i,i) - fn * sqrt((z - ephoton) / z) * akappa(z - ephoton,ephoton,qmemk)
            imk = i - kep
            if (imk >= 1) then
              cf(i,imk) = cf(i,imk) + fn * akappa(imk*du,ephoton,qmemk)
            endif
          endif

          ipk = i + kep
          if (ipk <= n) then
            cf(i,ipk) = cf(i,ipk) + fn * sqrt(z / (z + ephoton)) * ak
          endif
      enddo

      return
      end subroutine PHOTON

!*--------------------------------------------------------------------72

      subroutine BPHOTON (arate, erate, f, u, n)
!   calculates absorption and emission rates for free-free processes

      use RUN_PARAMETERS

      implicit none
      integer n
      double precision u(n)
      double precision f(n)
      double precision :: arate, erate

      integer i, k
      double precision :: Q0, const1, const2
      double precision :: akappa, ee, ep
      double precision :: du,  fn
      double precision :: z, qme, qm1, qmemk, qm2, ak

      double precision :: ta, ta2, te, te2, suma, sume

      data Q0 / 1.e-16 /
      data const1 /3.0664e-32/, const2 /5.9308e+7/

!*--------------------------------------------------------------------72
      akappa(ee,ep,qm1) = (const1 / ep ** 2) * const2 * sqrt(ee + ep) * ( (ee + .5 * ep) / ep) * qm1
!*--------------------------------------------------------------------72

      if (photflux <= 0.0) then
        arate = 0.0
        erate = 0.0
        return
      endif

      du   = u(2)-u(1)
      ta   = 0.0
      suma = 0.0
      te   = 0.0
      sume = 0.0
      do i = 0 , n
          z    = u(i)
          call QMOMNTM (z + .5 * ephoton,qm1,qm2)
          ta2  = akappa(z,ephoton,qm1) * sqrt(z) * f(i)
          suma = suma + ta + ta2
          ta   = ta2

          if (z > ephoton) then
            call QMOMNTM (z - .5 * ephoton,qm1,qm2)
            te2  = akappa(z - ephoton,ephoton,qm1) * sqrt(z - ephoton) * f(i)
            sume = sume + te + te2
            te   = te2
          endif
      enddo

      fn    = ephoton * photflux * du * Q0
      arate = fn * suma
      erate = fn * sume

      return
      end subroutine BPHOTON

!*--------------------------------------------------------------------72

      subroutine QMOMNTM (x, theta, theta2)

!   interpolates the momentum transfer cross section

      use ATOMIC_and_PLASMA_STUFF, only: tev, denow, zave, pi, emass, amass, iatomz   ! FLY stuff

      implicit none
      double precision :: x, theta, theta2
      double precision :: tel, edens, ave_Z, atwt, QION
      double precision :: atmass ! function -zchkstr-

      tel   = tev
      edens = denow
      ave_Z = zave
      atwt  = atmass (iatomz)
     ! write(*,*) 'te, ne, zbar, atmass 1451 zelda.f90 *****'
     ! write(*,*) tel, edens, ave_Z, atwt !mscho 200910p 
      call Q_ELECTRON_ION (tel, edens, ave_Z, QION)

!      QION = 1.0

      theta  = QION                              ! Qm 10^-16 cm^2

      theta2 = (2 * emass / atwt) * QION         ! (2m/M)*Qm

      return

      CONTAINS

        subroutine Q_ELECTRON_ION (tel, edens, ave_Z, QION)

          implicit none
          double precision :: tel, edens, ave_Z, QION
          double precision ::  elambda, dlogelam

!             Coulomb log and collision frequency relations from
!             Rockwood (1973); originally from Spitzer:
!           gamma = 2.58e+9 * (3. * tel) * sqrt(tel / edens)

!             New Couloumb log relations from Zollweg & Liebermann (1987):
!             Note <Z> factors:
            elambda = sqrt( (743. * sqrt(tel / edens))**2 + 1. / (4.*3.141592 * edens/3.)**(0.66667) )  &
                      / (ave_Z * 4.8e-8 / tel)
            dlogelam = log( dsqrt( 1. + 1.4 * elambda**2) )

!             electron - ion momentum transfer cross section:
!           QION = 4.*651.4 * log(gamma) / (x + 1.e-6) ** 2
            QION = 4.*651.4 * ave_Z**2 * dlogelam / (x + 1.e-6) ** 2

          return

        end subroutine Q_ELECTRON_ION

      end subroutine QMOMNTM

!*--------------------------------------------------------------------72

      function GAMMLN (xx)

      implicit none

      integer j 

      double precision xx, GAMMLN

!     Gamma function

      double precision cof(6),stp,half,one,fpf,x,tmp,ser

!      dimension cof(6)

      data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,-1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/

      data half,one,fpf/0.5d0,1.0d0,5.5d0/

      x   = xx - one
      tmp = x + fpf
      tmp = (x + half) * log(tmp) - tmp
      ser = one
      do j=1, 6
        x   = x + one
        ser = ser + cof(j) / x
      enddo
      GAMMLN = tmp + log(stp * ser)

      return
      end function GAMMLN

! --------------------------------------------------------------------72

      subroutine PRINT_EEDF (f0, e0, n)

!   The following prints & plots the distribution function


      use EEDF_DIMENSIONS

      use FILENAMES_and_UNITS

      use EFLAGS

      implicit none

      integer n, i
      integer mcols, nrows, nrmc, nlines, nunit

      double precision  e0(0:n), f0(0:n)
      double precision  fmatrix(NPTSF2)

      ifpr = .FALSE.
!      ifpr = .TRUE.


      IF (ifpr) THEN
        write (std_out,'(18X,A)') ' D I S T R I B U T I O N     F U N C T I O N'
        write (std_out,'(4(1x,a6,a11))') 'eV','eV**3/2', 'eV','eV**3/2','eV','eV**3/2', 'eV','eV**3/2'
        write (std_out,'(1x,1p,8e9.2)') (e0(i),f0(i), i=0,n)
        write (std_out,'(4(1x,a6,a11))') 'eV','eV**3/2', 'eV','eV**3/2','eV','eV**3/2', 'eV','eV**3/2'
      ENDIF

      ifpl = .FALSE.
!      ifpl = .TRUE.

      IF (ifpl) then
        mcols  = 2
        nrows  = n + 1
        nrmc   = nrows * mcols
        nlines = 51
        nunit  = std_out
        CALL PRPLT (e0,f0,fmatrix,nrows,mcols,nrmc,nlines,nunit)
      endif

      call flush (nunit)

      return
      end subroutine PRINT_EEDF

!*--------------------------------------------------------------------72

      subroutine PRPLT (energy, fofe, fmatrix, nrows, mcols, nrmc, nlines, nunit)

      implicit none

!       printer plotting routine
!       original version obtained from Rubin Shuker in 1977


      integer nrows, nrmc, mcols, nlines, nunit
      integer ifmt, nll, kn, i, j, l, ll, jp, ix, my, iz

      double precision  energy(nrows), fofe(nrows)
      double precision ypr(11), fmatrix(nrmc)
      double precision xscal, ymax, ymin, yscal, xb, xpr
      double precision z

      character*1 blank, aster, out
      dimension out(101) 

!     dimension char(9)
!     data char/1h1,1h2,1h3,1h4,1h5,1h6,1h7,1h8,1h9/
      data blank/' '/
      data aster/'*'/

      write (nunit,10)
   10 format (a1)

      do i = 1,nrows
        fmatrix(i) = energy(i)
        if (fofe(i) < 0.0) fofe(i) = 1.e-9
        fmatrix(nrows+i) = log10(fofe(i))
      enddo
      ifmt = 0
      nll  = nlines

!       check plot length, if zero set to 51.
  100 if (nll /= 0) go to 120
  110 nll   = 51
  120 xscal = (fmatrix(nrows) - fmatrix(1)) / float(nll - 1)

      ymax  = 1.
      ymin  = -9.
      yscal = (ymax - ymin) / 58.

!       compute and print the top axis.
      write (nunit,280)
      ypr(1) = ymin
      do kn = 1,9
        ypr(kn+1) = ypr(kn) + yscal * 5.8
      enddo
      ypr(11) = ymax
      if (abs(ypr(1)-ypr(2))<.01 .or. abs(ymax)>9.9d6) go to 150
      ifmt = 1
      write (nunit,140) (nint(ypr(ix)),ix = 1,11)
  140 format (9x,11(4x,i2))
      go to 170
  150 write (nunit,160) (ypr(ix),ix = 1,11)
  160 format (15x,11e10.2)
  170 write (nunit,180)
  180 format (14x,'+     +     +     +     +     +     +     +     + &
              +     +')
      xb = fmatrix(1)
      l = 1
      my = mcols - 1

!       compute and print the plot, a line for each iteration.
      do 250 i = 1,nll

!     a fuzz factor of 1.0e-12 needs to be added to xpr.
        xpr = xb + (i-1) * xscal + 1.0d-12
        if (xpr < fmatrix(l)) go to 230

        do ix = 1,64
          out(ix) = blank
        enddo

        do 200 j = 1,my
          ll = l + j * nrows
          z = fmatrix(ll)
          if (z>1. .or. z<-9.) go to 200
          jp = int((z - ymin) / yscal) + 1
!          out(jp) = char(j)
          out(jp) = aster
  200   continue

        write (nunit,210) xpr,(out(iz),iz = 1,64)
  210   format (1x,f8.3,'>',6x,64a1)

  220   l = l + 1
        if (l > nrows) go to 260
        if (fmatrix(l) < xpr) go to 220
        if (fmatrix(l) >= xpr) go to 250
  230   write (nunit,240) xpr
  240   format (5x,f10.3,'>')

  250 continue

!     print the bottom axis.
  260 write (nunit,180)
      if (ifmt == 0) go to 290
      write (nunit,140) (nint(ypr(ix)),ix = 1,11)
      write (nunit,280)
  280 format (1x,'Energy (eV)',5x,'Log of Electron Energy Distribution &
              Fiunction (eV^-3/2)')
      return

  290 write (nunit,160) (ypr(ix),ix = 1,11)
      write (nunit,280)

      return
      end subroutine PRPLT

!*--------------------------------------------------------------------72

      subroutine BANNER (iunit)

        implicit none
        integer iunit

!       write ELENDIF banner to unit number iunit

      write (iunit,*)
      write (iunit,9000)

 9000 FORMAT(//31x,'*****************'/                                    &
               31x,'*               *'/                                    &
               31x,'* E L E N D I F *'/                                    &
               31x,"*       96      *"/                                    &
               31x,'*****************'//                                   &
               17x,"TIME-DEPENDENT SOLUTION OF BOLTZMANN'S EQUATION"/      &
               18x,'FOR THE ELECTRON ENERGY DISTRIBUTION FUNCTION'/        &
               26x,'IN A PARTIALLY IONIZED PLASMA'//                       &
               34x,'Version 2.0'//                                         &
         17x,'Copyright (C) W.L. Morgan, JILA Rpt. #19 (1979)'/            &
         17x,'Copyright (C) W.L. Morgan and B.M. Penetrante,'/             &
         23x,'Comp. Phys. Comm. 58:127-152 (1990)'/                        &
         5x,'Complete Revision: Copyright (C) W.L. Morgan,',               &
         ' Kinema Software (1993,1995)')

      return
      end subroutine BANNER

!*--------------------------------------------------------------------72



      module CollisionArrays

        use EEDF_DIMENSIONS

        double precision, dimension (NPTSF,NPTSF) :: cf, aee

        double precision, dimension (NPTSF) :: tim1, ti, tip1

        double precision :: h1

      end module CollisionArrays

!-------------------------------------------------------------

      subroutine ee_Collisions (tlower, tupper, ystart)

      use CollisionArrays

      implicit none

      INTEGER KMAXX,NVAR
      PARAMETER (KMAXX=200,NVAR=500)
      INTEGER i,kmax,kount,nbad,nok,nrhs, ncall

      double precision :: tlower, tupper
!      double precision :: dxsav,eps,hmin,h1,x1,x2,x,y,ystart(NVAR)
      double precision :: dxsav,eps,hmin,x1,x2,x,y,ystart(NVAR)

      common /count/ ncall

!      COMMON /path/ kmax,kount,dxsav,x(KMAXX),y(NMAX,KMAXX)
!      COMMON nrhs

      EXTERNAL derivs,rkqs

      data ncall/0/

      ncall = ncall + 1

      eps = 1.e-4; hmin = 1.e-22 ! mscho 200326 hmin is the time step controller 

      if (ncall == 1) h1 = 1.e-20

      write (16,*) 'odeint: tl,tu,dt = ',real(tlower),real(tupper),real(h1)

      call odeint(ystart,NVAR,tlower,tupper,eps,h1,hmin,nok,nbad,derivs,rkqs)

      write (16,*) 'odeint: nok,nbad = ',nok,nbad


      return
      end subroutine ee_Collisions

!*--------------------------------------------------------------------72

      subroutine derivs (t, rhs, qee)

      use ATOMIC_and_PLASMA_STUFF, only: tev, denow, pi         ! FLY stuff

      use CollisionArrays

      use TIMING_STUFF, only: ee_scale

      implicit none

      integer i, k, l, n

      double precision :: t
      double precision rhs(NPTSF),  qee(NPTSF)
      double precision :: tel, edens, elambda, dlogelam, alf, c, dt

      tel   = tev
      edens = denow

!-------------------------------------------------------------------------
!       Coulomb log and collision frequency relations from
!       Rockwood (1973); originally from Spitzer:
!     emv2 = 3. * tel

!     lambda = (ergev ** 1.5) * sqrt(tel / (4. * pi * edens * esu ** 2))
!              * mv2 / (2. * esu ** 2)
!     elambda = 2.58126e+9 * emv2 * sqrt(tel / edens)

!--------------------------------------------------------------------------
!       New Couloumb log relations from Zollweg & Liebermann (1987):
      elambda = sqrt( (743. * sqrt(tel / edens))**2 + 1. / (4.*Pi * edens/3.)**(0.66667) ) / (4.8e-8 / tel)
      dlogelam = log( sqrt( 1. + 1.4 * elambda**2) )
!---------------------------------------------------------------

!     alf = (2. * pi * esu ** 4/3.) * sqrt(2. / emass) * log(lambda)
!              / ergev ** 1.5
!     alf = 2.57569e-6 * log(elambda)
      alf = 2.57569e-6 * dlogelam

!     Apply ee_scale factor (from TIMING_STUFF module).
!     ee_scale=1.0 (default, physical). Set > 1 to artificially accelerate
!     e-e thermalization (e.g., for PRE 2024 Fig.1c reproduction).
      alf = alf * ee_scale

      tim1 = 0.; ti = 0.; tip1 = 0.

      n = NPTSF

      do l = 1 , n

          do k = 2 , n
              tim1(k) = tim1(k) + aee(k - 1,l) * rhs(l)
          enddo

          do k = 1 , n
              ti(k)   = ti(k)   - (aee(k,l) + aee(l,k)) * rhs(l)
          enddo

          do k = 1 , n-1
              tip1(k) = tip1(k) + aee(l,k + 1) * rhs(l)
          enddo

      enddo

!   explicit treatment of electron-electron collisions:

        qee = 0.

        qee(1) = qee(1) + alf * (ti(1) * rhs(1) + tip1(1) * rhs(2))
        qee(n) = qee(n) + alf * (tim1(n) * rhs(n - 1) + ti(n) * rhs(n))

        do i = 2 , n-1
            qee(i) = qee(i) + alf * (tim1(i) * rhs(i - 1) + ti(i) * rhs(i) +  tip1(i) * rhs(i + 1))
        enddo

!        rhs = rhs + dt * alf * qee

        return

      end subroutine derivs

!*--------------------------------------------------------------------72
!                          ODE Integrator
!*--------------------------------------------------------------------72

      SUBROUTINE odeint (ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,rkqs)
      implicit none

      INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
      double precision :: eps,h1,hmin,x1,x2,ystart(nvar),TINY
      EXTERNAL derivs,rkqs
      PARAMETER (MAXSTP=100000,NMAX=500,KMAXX=200,TINY=1.e-30)
      INTEGER i,kmax,kount,nstp

      double precision :: dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),yp(NMAX,KMAXX),yscal(NMAX)
      COMMON /path/ kmax,kount,dxsav,xp,yp

      x=x1
      h=sign(h1,x2-x1) !mscho 200326 : sign(x,y) = abs(x) if y>0 / sign(x,y) = -abs(x) if y<0  
      nok=0
      nbad=0
      kount=0
      do i=1,nvar
        y(i)=ystart(i)
      enddo

      if (kmax > 0) xsav=x-2.*dxsav

      OuterLoop: do nstp=1,MAXSTP

        call derivs(x,y,dydx)

        do i=1,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
        enddo

        if(kmax > 0)then
          if(abs(x-xsav) > abs(dxsav)) then
            if(kount < kmax-1)then
              kount=kount+1
              xp(kount)=x
              do i=1,nvar
                yp(i,kount)=y(i)
              enddo
              xsav=x
            endif
          endif
        endif

        if((x+h-x2)*(x+h-x1) > 0.) h=x2-x

        call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs) ! mscho 200326 the value "hnext" is the same as h (but sign can be different)  

        if(hdid == h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1) >= 0.)then
          do i=1,nvar
            ystart(i)=y(i)
          enddo

          h1 = h

          if(kmax /= 0)then
            kount=kount+1
            xp(kount)=x
            do i=1,nvar
              yp(i,kount)=y(i)
            enddo
          endif
          return
        endif

        if(abs(hnext) < hmin) then
        !   pause  'stepsize smaller than minimum in odeint' 
        !mscho 200325 - pause command itself gives "To resume execution, type : go" intrinsically 
        endif
        h=hnext

      enddo OuterLoop

      ! pause 'too many steps in odeint'!mscho 200325 

      return
      END SUBROUTINE odeint

      ! mscho 200326 subroutine structure : ee_Collision > odeint > rkqs > rkck > derivs 
      SUBROUTINE rkqs (y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
        implicit none

      INTEGER n,NMAX
      double precision :: eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      EXTERNAL derivs
      PARAMETER (NMAX=500)

!     USES derivs,rkck
      INTEGER i
      double precision :: errmax,h,htemp,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,PSHRNK,ERRCON
      PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4)

      h=htry

    1 call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)

      errmax=0.
      do i=1,n
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
      enddo
      errmax=errmax/eps
      if(errmax > 1.)then
        htemp=SAFETY*h*(errmax**PSHRNK)
        h=sign(max(abs(htemp),0.1*abs(h)),h)
        xnew=x+h
        if(xnew == x) pause 'stepsize underflow in rkqs' !mscho 200325
        goto 1
      else
        if(errmax > ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.*h
        endif
        hdid=h
        x=x+h
        do i=1,n
          y(i)=ytemp(i)
        enddo
        return
      endif

      END SUBROUTINE rkqs


      SUBROUTINE rkck (y,dydx,n,x,h,yout,yerr,derivs)
        implicit none

      INTEGER n,NMAX
      double precision :: h,x,dydx(n),y(n),yerr(n),yout(n)
      EXTERNAL derivs
      PARAMETER (NMAX=500)
!     USES derivs
      INTEGER i

      double precision :: ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX),  &
      ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,         &
      B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6

      PARAMETER (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,     &
      B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,           &
      B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.,          &
      B63=575./13824.,B64=44275./110592.,B65=253./4096.,C1=37./378.,    &
      C3=250./621.,C4=125./594.,C6=512./1771.,DC1=C1-2825./27648.,      &
      DC3=C3-18575./48384.,DC4=C4-13525./55296.,DC5=-277./14336.,       &
      DC6=C6-.25)

      do i=1,n
        ytemp(i)=y(i)+B21*h*dydx(i)
      enddo
      call derivs(x+A2*h,ytemp,ak2)
      do i=1,n
        ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
      enddo
      call derivs(x+A3*h,ytemp,ak3)
      do i=1,n
        ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
      enddo
      call derivs(x+A4*h,ytemp,ak4)
      do i=1,n
        ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
      enddo
      call derivs(x+A5*h,ytemp,ak5)
      do i=1,n
        ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+B65*ak5(i))
      enddo
      call derivs(x+A6*h,ytemp,ak6)
      do i=1,n
        yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
      enddo
      do i=1,n
        yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*ak6(i))
      enddo

      return
      END SUBROUTINE rkck

!*--------------------------------------------------------------------72
!*--------------------------------------------------------------------72

