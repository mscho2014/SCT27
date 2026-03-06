      subroutine DCOEFF (du, n, ndim, iflag, ncall, cf)


!   this subroutine computes the matrix of coefficients for
!   a solution by finite-differences

!     use MISC_PHYSICAL_QUANTITIES, only: alpha, tgev


!       11/08/02 (in Berkeley): Hyun's stuff included in DCOEFF & zhk_stuff.f90 replaces flystuff.f90


      use EEDF_DIMENSIONS

      use POP_STUFF

      use ATOMIC_and_PLASMA_STUFF

      use FofE_STUFF

!     use popstuf, ONLY:llower

      implicit none


      integer n, ndim, i, k, kmj, kpj, ncall, iflag, icheck
      integer itran, ietran, iso, llo, lup, ind1, ind2, ilo, iup
      !integer indxc
      integer j !mscho 200908 - for test
      
      double precision :: egain, eloss
      double precision :: dgain, dloss
      double precision :: cf(ndim,ndim)

      double precision  :: du, uplus, u
      double precision  :: N_lower
      double precision  :: N_upper
      double precision  :: etran, e

      double precision  :: q1 , q3
      double precision  :: qi1, qb1, qi2, qb2 

      double precision :: tgev, gdens


      double precision :: d1,d2,duby4,const1,const2
      double precision :: ai,aim1,bi,bip1,um1
      double precision :: qm1, qm2, qmum1, qm2um1
      double precision :: ak, bkp1, z, qel, qel2
      double precision :: alpha, qr1, qr2
      
      double precision :: xc_chy, xc_che, xc_cli, xc_chl, xc_caa
      double precision :: xc_sxx, xc_betaxx, xc_gsxx, xc_gbetaxx
      double precision :: debye
     
! check rates
      double precision :: esum, dsum
      double precision :: Vsxx, Vbetaxx
      integer   n1
      !double precision, dimension(500) :: energyz
      real*8 :: fxxa 
      real*8 :: fxa 
      integer, external :: indxc  
      !INTEGER FUNCTION indxc(fxx1,fxx2)
      !   real*8 :: fxx1, fxx2
      !END FUNCTION indxc 
      
!------------------------------------------------
!       original Rockwood finite difference scheme:
      ak(z,qel,qel2) = const1 * alpha * (d2 / (const2 * sqrt(z) * qel))  &
                       * (z + duby4) + .5 * d1 * const2 * sqrt(z) *      &
                       ((qel2 + qr1) * (-z) + (qel2 + qr2) * tgev *      &
                       (.5 + 2. * d1 * z))

      bkp1(z,qel,qel2) = const1 * alpha * (d2 / (const2 * sqrt(z) *      &
                         qel)) * (z - duby4) + .5 * d1 * const2 *        &
                         sqrt(z) * ((qel2 + qr1) * z + (qel2 + qr2) *    &
                         tgev * (-.5 + 2. * d1 * z))

!------------------------------------------------

      alpha = 0.   ! unless E/N > 0
      qr1 = 0.     ! rotational cross-sections
      qr2 = 0. 

      fxxa = 0. 

!     const1 = 2 e/m = 2.*(1.602192e-19*1.0e+7)**2/(1.602192e-12*9.10956e-28)
      const1 = 3.51761e+15
      duby4  = .25 * du

!     const2 = sqrt(2 e/m) = sqrt(2.*1.602192e-12/9.10956e-28)
      const2 = 5.93094e+7
!      write(*,*)'******* print iflag *****', iflag !mscho 200908
     ! write (*,*) '*******print tel ******', tel ! mscho 200512

      if (iflag == 0) then
         cf = 0.
         !write(*,*) 'iflag 000000' !mscho : iflag == 0  
      endif
!     write(*,*) cf ! mscho 200908
!      if (enflag /= 'DC') then
!        cfim1 = 0.; cfi = 0.; cfip1 = 0.
!      endif

!------------------------------------------------

!     BUG FIX (2024): tgev was hardcoded to 1.0 eV, ignoring actual ion temperature.
!     Now uses tiev (ion temperature in eV) passed from SCFLY.
!     For cold gas targets (e.g., Ne at 300K ~ 0.026 eV), this significantly
!     affects the elastic thermalization rate in the Rockwood FD scheme.
      tgev  = tiev

      gdens = dinow

      !write(*,*)'this is line 103 in zelda.f90: print out gdens(=dinow)'
      !write(*,*) gdens ! mscho 201102 OK (gdens = 6.00001e22) 
 
      d1   = 1.0 / du
      d2   = d1 ** 2

      !qm2  = 6.3240402667679558E-32
      qm2 = 0

      !mscho modification (201015) 
      !denergy(1) = 0.25
      !denergy(2) = 0.25

      !write(*,*) ' this is n, line 120 in zcoeff', n 
      !do i = 1, 3
      !  write(*,*) (cf(i,j),j=1,500)
        !write(*,*)'122 zcoeff:energy(i) denergy',i, energy(i),denergy(i)
        !energy(i)=energy(i+1) 
      !enddo 
      !write(*,*)'energy(i) in zcoeff 121 line',ncall, energyz(1), energyz(2)
      Energy_Loop_Elastic: do i = 1, n

          uplus = energy(i) + 0.5*denergy(i)
          !uplus = u(i) + 0.5*denergy(i)
          u     = energy(i)
          !u     = u(i)
          ! write(*,*) 'mscho uplus, u, qm1, qm2 in 108 line zco...f' 
          ! write(*,*) uplus, u, qm1, qm2 ! mscho 200910p - no problem (but different with ct27)  
          ! mscho  uplus:energy grid middle point / u: energy grid / qm1 & qm2 are distribution   
! ---------------------------------------------------------------------
!
!         elastic collision (tridiagonal) terms:
          !write(*,*) 'ok'
          call QMOMNTM (uplus, qm1, qm2)
          ai = dble(ak(uplus, qm1, qm2) / const2)
          if (i == n) ai = 0.0
          bip1 = dble(bkp1(uplus, qm1, qm2) / const2)
          if (i == 1) then
            bi   = 0.0d0
            aim1 = 0.0d0
           else
            um1  = uplus - du
            call QMOMNTM (um1, qmum1, qm2um1)
            aim1 = dble(ak(um1, qmum1, qm2um1) / const2)
            bi   = dble(bkp1(um1, qmum1, qm2um1) / const2)
          endif

!          ai = 0.; bi = 0.; aim1 = 0.; bip1 = 0.! ***** TESTING ENERGY BALANCE *****

! ---------------------------------------------------------------------
            cf(i,i)  =  - (ai + bi)
            if (i < n) cf(i,i+1) = bip1
            if (i > 1) cf(i,i-1) = aim1

! ---------------------------------------------------------------------
!         Other Energy Loss Processes

!         call dEdt_ee (u, dEdtee, PlasmonEnergy) ! energy loss from plasma electrons to conduction electrons (eV/sec)
!         iPlasmon = indxc(PlasmonEnergy)
!         if (i > 1) then
!           if (i > iPlasmon) then     ! discrete energy loss
!
!             cf(i,i)= cf(i,i)-(dEdtee/PlasmonEnergy)/(5.93094e-9 * gdens)  ! it gets multiplied again by 5.93094e-9 & gdens later
!             cf(i-iPlasmon,i) = cf(i-iPlasmon,i) + (dEdtee/PlasmonEnergy)/(5.93094e-9 * gdens)
!
!           else                       ! continuous energy loss
!
!             cf(i,i)   = cf(i,i)   - (dEdtee/du)/(5.93094e-9 * gdens)  ! it gets multiplied again by 5.93094e-9 & gdens later
!             cf(i-1,i) = cf(i-1,i) + (dEdtee/du)/(5.93094e-9 * gdens)
!
!           endif
!
!         endif

!         call dEdt_el (u, dEdtel)! energy loss from plasma electrons to lattice

!-----------------------------------------------------------------------------

      enddo Energy_Loop_Elastic   
      

      !do i=1,5
      !   write(*,*) (cf(i,k),k=1,500)
      !enddo  ! mscho 200910 - It is problem: it provides NAN value in the cf(2,k) (k=1,500) 
!-----------------------------------------------------------------------------
!                               Inelastic processes
!-----------------------------------------------------------------------------

      do iso=1,iatomz
         call Vcaafill(iso)
      enddo


      n1=0
      egain=0.
      eloss=0.
      dgain=0.
      dloss=0.

      do k=1,3
            !write(*,*)'fe(k) 206 zcoeff', fe(k)
            !fe(k)=fe(4)
      enddo

      Transition: do  itran=1, ntran

         llo=llower(itran)
         lup=lupper(itran)
         ind1=indl2p(llo)
         ind2=indl2p(lup)
         if(ind1*ind2.eq.0) cycle Transition

!        TESTING BB TRANSITIONS
         n1=n1+1
         N_upper=popout(lup)/gdens 
         if(N_upper.lt.1.e-30)cycle Transition
         N_lower=popout(llo)/gdens
         if(N_lower.lt.1.e-30)cycle Transition
         ilo=levl(llo)
         iup=levl(lup)
         iso=isolev(llo)


         if(ltype(itran).eq.1)then !For bound-bound(bb) transition 

            etran=elev(lup)-elev(llo) !energy level difference in specific b-b transition 
            ietran=indxc(etran)
! hyun's test
!            if(etran.le.denergy(1)) cycle Transition
            if(etran.le.exthr) cycle Transition
            if(ietran.le.0.or.ietran.gt.nebins)cycle Transition

            esum=0.d0
            dsum=0.d0
            do 100 k=ietran+1, n
                kmj=indxc(energy(k)-etran) !This is for collisional (de)excitation process 
                !kmj=indxc(20,0)
!               kmj=k-ietran
               if(kmj.le.0) goto 100
               if(iso.le.3.and.flyflag.ne.'off')then
                  if(flyflag.eq.'fly')then
                     if(iso.eq.1)then ! Collisional-CX define
                        q1 = xc_chy(k, ilo,iup) !Collisional CX for hydrogen-like atom 
                        q3 = xc_chy(k,iup,ilo)
                     else if(iso.eq.2)then
                        q1 = xc_che(k, ilo,iup) !Collisional CX for helium-like atom
                        q3 = xc_che(k,iup,ilo)
                     else 
                        q1 = xc_cli(k, ilo,iup) !Collisional CX for lithium-like atom 
                        q3 = xc_cli(k,iup,ilo)
                     endif
                  else if(flyflag.eq.'hul')then
                     q1 = xc_chl(k, ilo,iup,iso) ! Collisional CX for all atom in HULLAC atomic data
                     q3 = xc_chl(k, iup,ilo,iso)
                  endif
               else
                  q1 = xc_caa(k, ilo,iup,iso) ! Collisioonal CX for all atom (defined in advance) 
                  q3 = xc_caa(k, iup,ilo,iso)
               endif
               esum=esum+N_lower * q1 * fe(k)
               dsum=dsum+N_upper * q3 * fe(kmj)
               egain=egain-N_lower*q1*fe(k)*denow+N_upper*q3*fe(kmj)*denow

               !write(*,*)'test q1 mscho*** line 249 in zcoeff'! mscho: output is not shown in screen      
               q1 = q1 / (1.e-16*const2*energy(k)*denergy(k))
               !write(*,*) q1 !mscho: output is not shown in screen 
               cf(k,k) = cf(k,k) - N_lower * sqrt(energy(k)) * q1
               cf(kmj,k) = cf(kmj,k) + N_lower * sqrt(energy(k)) * q1
               
               q3 = q3 / (1.e-16*const2*energy(kmj)*denergy(kmj))
               cf(kmj,kmj) = cf(kmj,kmj) - N_upper * sqrt(energy(kmj)) * q3
               cf(k,kmj) = cf(k,kmj) + N_upper * sqrt(energy(kmj)) * q3
100         continue

            eloss=eloss-N_lower*Vcaa(ilo,iup,iso)*denow+N_upper*Vcaa(iup,ilo,iso)*denow

            ! check the balance : please confirm after 
            ! (1) Check the collisional excitation process
            if(abs(esum-N_lower*Vcaa(ilo,iup,iso)).gt.1.e-10*esum)then
              ! print *,'cx',iso,ilo,iup,esum,N_lower*Vcaa(ilo,iup,iso)
            endif
            ! (2) Check the collisional de-excitation process 
            if(abs(dsum-N_upper*Vcaa(iup,ilo,iso)).gt.1.e-10*dsum)then
              ! print *,'dx',iso,iup,ilo,dsum,N_upper*Vcaa(iup,ilo,iso)
            endif

         else if(ltype(itran).eq.2)then
            etran=elev(lup)-elev(llo)+eipz(iso)-debye(iatomz-iso+1)
            !write(*,*)'315 zcoeff:etran',etran
            if(etran.le.exthr) cycle Transition
            !ietran=indxc(180,0)
            fxa = REAL(etran)
            ietran=indxc(etran)

! Loss and gain by ionization
            esum=0.d0
            dsum=0.d0
            do 200 k=ietran, n
               e=energy(k)
               qi1=xc_sxx(k,ilo,iup,iso)/(1.e-16*const2*energy(k)*denergy(k))
               cf(k,k) = cf(k,k) - N_lower * sqrt(energy(k)) * qi1  

               dgain=dgain+N_lower*xc_sxx(k,ilo,iup,iso)*fe(k)*denow  
               esum=esum+N_lower*xc_sxx(k,ilo,iup,iso)*fe(k)  

               do 300 kmj=1, k-ietran+1
                  qi2=xc_gsxx(kmj,k,ilo,iup,iso)/(1.e-16*const2*energy(k)*denergy(k))
                  cf(kmj,k)= cf(kmj,k)+ N_lower*sqrt(energy(k)) * qi2
                  dsum=dsum+N_lower*xc_gsxx(kmj,k,ilo,iup,iso)*fe(k)  
 300           continue
200         continue

            if(abs(esum-dsum/2.).gt.1.e-10*dsum) then
               write(*,'("ion ",3i3,i4,1p,9e13.5)') &
                    iso,ilo,iup,ietran,esum,dsum/2.,Vsxx(ilo,iup,iso)*N_lower
            endif

            dloss = dloss+Vsxx(ilo,iup,iso)*N_lower*denow

! Loss and gain by 3-body recombination

            esum=0.d0
            dsum=0.d0
            do 400 k=1, n
               qb1 = xc_betaxx(k,ilo,iup,iso)/(1.e-16*const2*energy(k)*denergy(k))
               cf(k,k) = cf(k,k) - N_upper * sqrt(energy(k)) * qb1 * denow

               dgain=dgain-N_upper*xc_betaxx(k,ilo,iup,iso)*fe(k)*denow*denow/2.
               esum=esum+N_upper*xc_betaxx(k,ilo,iup,iso)*fe(k)


               do 500 kpj=k+ietran, nebins
                  qb2 = xc_gbetaxx(kpj,k,ilo,iup,iso)/(1.e-16*const2*energy(k)*denergy(k))
                  cf(kpj,k)= cf(kpj,k)+ N_upper*sqrt(energy(k)) * qb2 *denow
                  dsum=dsum+N_upper*xc_gbetaxx(kpj,k,ilo,iup,iso)*fe(k)
 500           continue
 400        continue

            if(abs(esum/2.-dsum).gt.1.e-10*dsum)then
               print *,'rec',iso,ilo,iup,esum/2,dsum,Vbetaxx(ilo,iup,iso)*N_upper
            endif

            dloss=dloss-N_upper*Vbetaxx(ilo,iup,iso)*denow*denow

         endif

      enddo Transition

!------------------------------------------------
!      DONE coefficients 
!------------------------------------------------

 999   cf = 5.93094e-9 * cf

      !do i=1,5
      !    write(*,*) (cf(i,j), j=1,500) !mscho 200908: Error: No inelastic collision
      !enddo 

      !do i=249,250
      !   write(*,*) (cf(i,j), j=1,10) !mscho 200908 
      !enddo 



!      print *, ' n1 = ', n1, popout(indp2l(9))
!      write(*, '(" egain-ele = ",1p,e15.8)') egain
!      write(*, '(" egain-ion = ",1p,e15.8)') eloss
!      write(*, '(" igain-ele = ",1p,e15.8)') dgain
!      write(*, '(" igain-ion = ",1p,e15.8)') dloss

      return

      end subroutine DCOEFF

!------------------------------------------------------------------------------
      subroutine dEdt_ee (u, dEdt, HbarOmegaPeV)     ! energy loss from plasma electrons to conduction electrons
        implicit none

      double precision :: u,dEdt, b0
      double precision :: Ne, Te, LambdaD, Ri, LambdaM, LogLambda, OmegaP, Ve, e, HbarOmegaP, HbarOmegaPeV

      Ne = 18.1d+22                         ! conduction electron density in Al (/cc)

      Te = 0.025                            ! I've taken this to be temperature of conduction electrons (eV)

      Ve     = 5.93e+7 * sqrt(u)            ! cm/sec

      OmegaP = 5.65e+4 * sqrt(Ne)           ! plasma frequency (/sec)

      HbarOmegaP = 1.055e-27 * OmegaP       ! erg

      HbarOmegaPeV = HbarOmegaP/1.6e-12     ! eV

      e      = 4.8e-10                      ! (erg-cm)^1/2

      if (u > HbarOmegaPeV) then

        LogLambda = log( 4.*u*1.6e-12 / HbarOmegaP)

        dEdt   = ( (Omegap * e)**2 / Ve ) * LogLambda / 1.6e-12        ! eV/sec; discrete

      else

        LambdaD = 742. * sqrt(Te/Ne)          ! Debye length (cm)

        Ri = (3 / (4*3.141592*Ne))**0.333333  ! mean ionic radius or ion sphere radius (cm)

        b0 = 7.2e-8 / u                       ! impact parameter

        LambdaM   = sqrt( LambdaD**2 + Ri**2 ) / b0

        LogLambda = log( 1 + 1.4 * LambdaM**2)     ! Coulomb logarithm based on Zollweg & Liebermann, JAP (1987)

        dEdt      = ( (Omegap * e)**2 / Ve ) * LogLambda / 1.6e-12        ! eV/sec; continuous

        dEdt      = dEdt * 0.63   ! use temporary scale factor to connect to discrete dEdt

      endif

!      write (20,*) u, dEdt

      return
      end  subroutine dEdt_ee

!------------------------------------------------------------------------------

