      subroutine gentrnx(iz)
      implicit real*8(a-h,o-z)
      include 'mainstuf'
      include 'flystuf'
      include 'popstuf'

      do i=2, nlev(iz)
         do j=1, i-1
            if(flyflag.eq.'fly') then
               if(iz.eq.1) then
                  oratexx(i,j,iz)=oratehy(i,j)
                  fxx(j,i,iz)=oratehy(i,j)*ghy(i)/ghy(j)
     &                 /(evhy(i)-evhy(j))**2/4.371d7
               else if(iz.eq.2)then
                  oratexx(i,j,iz)=oratehe(i,j)
                  fxx(j,i,iz)=oratehe(i,j)*ghe(i)/ghe(j)
     &                 /(evhe(i)-evhe(j))**2/4.371d7
               else if(iz.eq.3) then
                  oratexx(i,j,iz)=orateli(i,j)
                  fxx(j,i,iz)=orateli(i,j)*gli(i)/gli(j)
     &                 /(evli(i)-evli(j))**2/4.371d7
               else
                  oratexx(i,j,iz)=orateaa(i,j,iz)
                  fxx(j,i,iz)=faa(j,i,iz)
               endif
            else if(flyflag.eq.'hul')then
               if(iz.le.3)then
                  oratexx(i,j,iz)=oratehl(i,j,iz)                  
                  fxx(j,i,iz)=fhul(j,i,iz)
               else
                  oratexx(i,j,iz)=orateaa(i,j,iz)
                  fxx(j,i,iz)=faa(j,i,iz)
               endif
            else
               oratexx(i,j,iz)=orateaa(i,j,iz)
               fxx(j,i,iz)=faa(j,i,iz)
            endif
         enddo
      enddo

      return
      end
c---------------------------------------------------

c---------------------------------------------------
      function V2betaxx(i,k,iz)
      implicit real*8(a-h,o-z)
      include 'mainstuf'
      include 'popstuf'
      include 'xc_stuf'

      V2betaxx=0.d0
      if(.not.cxflag) return

      ntoti=indlev(i,iz)
      ntotj=indlev(k,iz-1)
      eii=eipz(iz)+elev(ntotj)-elev(ntoti)-debye(iatomz-iz+1)
      if(eii.le.exthr) return

      ieii = indxc(eii)
      pn=petrn(i,k,iz)
      qq=qlev(ntoti)-1
      bfact=0.25*((100.*qq+91.)/(4.*qq+3.))**0.5 - 1.25
      ci=1.4679d-22*glev(ntoti)/glev(ntotj)


      sum=0.d0
      do k2=1, nebins
         do j=k2+ieii,nebins
            e=energy(j)
            indx=indxc(e-energy(k2)-eii)
            if(e.gt.eii.and.indx.ge.1)then
               uu=e/eii
               etot=e-eii
               e2=energy(k2)/etot
               if(e2.le.1.d0)then
                  wfact=(log(e/eii))**(bfact*eii/e)
                  xctot=2.218d-6*pn/eii*log(e/eii)*wfact
                  qj=qdiff(etot,uu,e2)/2.d0
                  sum=sum+xctot*qj*
! both fe() are new
!     &                 fe(indx)*fe(k2)*denergy(k2)*denergy(j)! note denergy(j)
! use one old and one new fe()
     &                 fe0(indx)*fe(k2)*denergy(k2)*denergy(j)

               endif
            endif
         enddo
      enddo

      V2betaxx=ci*sum
      return
      end


c---------------------------------------------------
c
c---------------------------------------------------
c
c
      function capehe(i,k)
c
c___ rate of electron capture from hydrogenic ground state
c___     into autionizing state i
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'flystuf'
c
      ehe = evhe(i)-(evhy(k)+eipz(2)-debye(iatomz-1))
      if(ehe.le.0.) then
         capehe=0.d0
      else
         capehe=1.65615e-22/tev**1.5*ghe(i)/ghy(k)*exp(-ehe/tev)
     &        *autohe(i,k)
      endif
c
      return
      end
c
      function capeli(i,k)
c
c---- rate of electron capture from helium-like ground state
c----     into autionizing state i
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'flystuf'

      e= evli(i)-(evhe(k)+eipz(3)-debye(iatomz-2))
      if(e.le.0.) then
         capeli=0.d0
      else
         capeli=1.65615e-22/tev**1.5*gli(i)/ghe(k)*exp(-e/tev)
     &        *autoli(i,k)
      endif
      return
      end

      function capebe(i,j)
c
c---- rate of electron capture from the state (j, iz-1)
c----     into autionizing state (i,iz)
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'flystuf'
      include 'popstuf'

      iz=4
      jlv=indlev(j,iz-1)
      eii = evaa(i,iz)-elev(jlv)-eipz(iz)
     &     +debye(iatomz-iz+1)
      capebe = 1.65615d-22*exp(-eii/tev)/tev**1.5*
     &       gaa(i,iz)/glev(jlv)*autobe(i,j)
      return
      end
c
c

      function Vcapehe(i,j)
c
c___ rate of electron capture from hydrogenic ground state
c___     into autionizing state i
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'flystuf'
      include 'xc_stuf'
 
      Vcapehe=0.d0
      if(.not.cxflag) return
      ehe = evhe(i)-(evhy(j)+eipz(2)-debye(iatomz-1))
      glbygu = ghe(i)/ghy(j)
      k= indxc(ehe)
      if(k.gt.0.and.k.le.nebins) then
         Vcapehe=1.40e-22*autohe(i,j)*glbygu*energy(k)*fe(k)
      endif
      return
      end
c
      function Vcapeli(i,j)
c
c---- rate of electron capture from helium-like ground state
c----     into autionizing state i
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'flystuf'
      include 'xc_stuf'
 
      Vcapeli=0.d0
      if(.not.cxflag) return
      e= evli(i)-(evhe(j)+eipz(3)-debye(iatomz-2))
      k=indxc(e)
      if(k.gt.0.and.k.le.nebins) then
         Vcapeli=1.40e-22*autoli(i,j)*gli(i)/ghe(j)*energy(k)*fe(k)
      endif
      return
      end

      function Vcapebe(i,j)
c
c---- rate of electron capture from helium-like ground state
c----     into autionizing state i
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'flystuf'
      include 'xc_stuf'
      include 'popstuf'
c
      iz=4
      jlv=indlev(j,iz-1)
      Vcapebe=0.d0
      if(.not.cxflag) return

      eii = evaa(i,iz)-elev(jlv)-eipz(iz)
     &     +debye(iatomz-iz+1)

      glbygu = gaa(i,iz)/glev(jlv)
      k=indxc(eii)
      if(k.gt.0.and.k.le.nebins) then
         Vcapebe=1.40e-22*autobe(i,j)
     &        *glbygu *energy(k)*fe(k)
      endif

      return
      end

c
c
c **************************
      subroutine Vchyfill
c
c___ this calculates the collisional excitation or de-excitation
c___ from level i to level j.
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'flystuf'
      include 'xc_stuf'
c
c___ Fill the arrays of excitation and de-excitation cross-sections
c___ versus energy.  Noting that the 1/e dependence of the xc
c___ implied that, e.g., the collision strengths, e*xc(e) is energy
c___ independent and thus the excitation = de-excitation
c
      nhytot = 25
 
      do i=1,nhytot
        do j=1,nhytot
            Vchy(i,j) = 0.0
        enddo
      enddo
 
      do i=1,nhytot-1
        do 100 j=i+1,nhytot
          fij = fos(i,j)
          eji = evhy(j)-evhy(i)
          if(eji.le.exthr) goto 100
          ieji=indxc(eji)
          gij = ghy(i)/ghy(j)
          sumex=0.d0
          sumdx=0.d0
          do k=1,nebins
             e = energy(k)
             if (k-ieji.ge.1.and.e .ge. eji) then
                sumex=sumex+fe(k)*denergy(k)*2.823e-6*fij/eji
                sumdx=sumdx+fe(k-ieji)*denergy(k-ieji)*2.823e-6*fij/eji
             endif
          enddo 
          Vchy(i,j)=sumex
          Vchy(j,i)=sumdx
 100   continue
      enddo
c
      return
      end
c
c
      subroutine Vchefill
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'flystuf'
      include 'xc_stuf'

      common /dataa/ aaut(5,5),z2s2(5,5),coa(5,5),c1a(5,5),c2a(5,5)
      common /idatnum/ idatsam, idatboi
c
c___ calculate a constant.
c
      idsamp1 = idatsam+1
      iboiko = 1

      nbn0 = 25
      nbn = 13
      nau = 6
      nhetot = nau+nbn0
      zeff = atomz-1.
 
c   zero out the array for the helium-like collisions
 
      do i=1,nhetot
        do j=1,nhetot
            Vche(i,j) = 0.0
        enddo
      enddo
        
c
c___ calculate the collisional rates between helium like levels i to j
c
      do 2000 kku=2,nhetot
      do 2000 kkl=1,kku-1

        kki = kku
        if (kki .gt. nbn) kki = kki-nbn-nbn0
        kkj = kkl
        if (kkj .gt. nbn) kkj = kkj-nbn-nbn0
        ku=kki
        kl=kkj

        eu = evhe(kku)
        el = evhe(kkl)
        gu = ghe(kku)
        gl = ghe(kkl)
        delta = eu-el
        ieji=indxc(delta)
        if(delta.le.exthr)goto 2000

        do 100 k=1,nebins

           e = energy(k)
           ebyd = e/delta
           zbyd = ebyd+1.
           partex = 0.d0
           partde = 0.d0


           if(kkl.gt.nbn.and.kkl.le.nbn0)goto 5000
           if(kku.gt.nbn.and.kku.le.nbn0)goto 5000
c
c***********************************************************
c*                                                         *
c*  note:  for boiko data have all bound levels.           *
c*         for sampson data have detail for 5 levels.      *
c*                                                         *
c***********************************************************
c
c
c---- determine which set of rates to use
c
           if (ku.le.idatboi.and.iboiko.eq.1) go to 3000
c
c---- use sampsons rates - golden et.al. ap.j.supp.ser.,v45,p603,(1981)
c
           if (kl.gt.nbn) go to 109
           if (kl.gt.idatsam) go to 106
           go to (101,102,103,104,105) , kl
c
c___ treat autoionizing level collisions differently
c
 101       continue
c
c___ kl is the 1s state
c
           khl=1
           if (ku.gt.nbn) go to 5000
           if (ku.gt.idsamp1) go to 4000
           if (ku.eq.idsamp1) go to 111
           nd=ku/2
           partde = 1.5*xcsamex(zbyd,nd)
           if (ku-nd*2.gt.0) partde = -partde+2.*xcsamd(zbyd,nd)
           if (ebyd .ge. 1.) then
              partex = 1.5*xcsamex(ebyd,nd)
              if (ku-nd*2.gt.0) partex=-partex+2.*xcsamd(ebyd,nd)
           endif
           go to 1000
 111       continue
           partde = 2.*(xcsamd(zbyd,3)+xcsamd(zbyd,4)+xcsamd(zbyd,5))
           if (ebyd .ge. 1.) then      
              partex = 2.*(xcsamd(ebyd,3)+xcsamd(ebyd,4)+xcsamd(ebyd,5))
           endif
           go to 1000
c
c___ 2 triplet s state
c
 102       continue
           khl=2
           if (ku.gt.idsamp1) go to 4000
           if (ku.eq.idsamp1) go to 152
           if (ku-4) 122,132,142
 122       continue
           partde = .75*xcspnn1(5,zbyd)
           if (ebyd .ge. 1.) then
              partex = .75*xcspnn1(5,ebyd)
           endif
           go to 1000
 132       continue
           partde = xcspnno(1,5,zbyd)
           if (ebyd .ge. 1.) then
              partex = xcspnno(1,5,ebyd)
           endif
           go to 1000
 142       continue
           partde = .75*xcspnn1(1,zbyd)
           if (ebyd .ge. 1.) then
              partex = .75*xcspnn1(1,ebyd)
           endif
           go to 1000
 152       continue
           partde = 3.*(xcsamd(zbyd,6)+xcsamd(zbyd,7)+xcsamd(zbyd,8))
           if (ebyd .ge. 1.) then
              partex = 3.*(xcsamd(ebyd,6)+xcsamd(ebyd,7)+xcsamd(ebyd,8))
           endif
           go to 1000
c     
c___  2 singlet s state
c
 103       continue
           khl=2
           if (ku.gt.idsamp1) go to 4000
           if (ku.eq.idsamp1) go to 113
           if (ku-4) 123,123,133
 123       continue
           partde = .75*xcspnn1(1,zbyd)
           if (ebyd .ge. 1.) then
              partex = .75*xcspnn1(1,ebyd)
           endif
           go to 1000
 133       continue
           partde = xcspnno(1,4,zbyd)
           if (ebyd .ge. 1.) then
              partex = xcspnno(1,4,ebyd)
           endif
           go to 1000
 113       continue
           partde = (xcsamd(zbyd,6)+xcsamd(zbyd,7)+xcsamd(zbyd,8))
           if (ebyd .ge. 1.) then
              partex = (xcsamd(ebyd,6)+xcsamd(ebyd,7)+xcsamd(ebyd,8))
           endif
           go to 1000
c
c___ 2 triplet p state
c
 104       continue
           khl=2
           if (ku.gt.idsamp1) go to 4000
           if (ku.eq.idsamp1) go to 114
           partde = .75*xcspnn1(6,zbyd)
           if (ebyd .ge. 1.) then
              partex = .75*xcspnn1(6,ebyd)
           endif
           go to 1000
 114       continue
           partde = 3.*(xcsamd(zbyd,9)+xcsamd(zbyd,10)+xcsamd(zbyd,11))
           if (ebyd .ge. 1.) then
              partex = 3.*(xcsamd(ebyd,9)+xcsamd(ebyd,10)+
     &             xcsamd(ebyd,11))
           endif
           go to 1000
c     
c___ 2 singlet p state
c
 105       continue
           khl=2
           if (ku.gt.idsamp1) go to 4000
           if (ku.eq.idsamp1) go to 107
 107       continue
           partde = (xcsamd(zbyd,9)+xcsamd(zbyd,10)+xcsamd(zbyd,11))
           if (ebyd .ge. 1.) then
              partex = (xcsamd(ebyd,9)+xcsamd(ebyd,10)+xcsamd(ebyd,11))
           endif
           go to 1000
c
c---- use Vinogradovs rate for n = 1,2,3,4,5 - sov.j.quan.ele.,5,p630,(1978)
c
 3000      continue
           partde = xc_russian(kl,ku,zbyd)
           if (ebyd .ge. 1.) then
              partex = xc_russian(kl,ku,ebyd)
           endif
           go to 1000
c     
c___  levels treated as hydrogenic n .ge. 3
c
 106       continue
           khl=kl-3
           go to 4000
c
c___ use sampsons data for collisions between auto-ionizing levels
c___            goett etal. atomic and nuc. data tables v25, p185 (1980)
c     
 109       continue
           kal = i-nbn
           kau = j-i
           partde = coa(kau,kal)+1.66*z2s2(kau,kal)*log(zbyd)+
     1          c1a(kau,kal)/(aaut(kau,kal)+zbyd)+
     2          c2a(kau,kal)/(aaut(kau,kal)+zbyd)**2
           if (ebyd .ge. 1.) then
              partex = coa(kau,kal)+1.66*z2s2(kau,kal)*log(ebyd)+
     1             c1a(kau,kal)/(aaut(kau,kal)+ebyd)+
     2             c2a(kau,kal)/(aaut(kau,kal)+ebyd)**2
           endif
           go to 1000
           
 4000      continue
c
c___ use hydrogenic rates for bound levels above n = 3. Note that
c___ to keep the definition consistent with the collision strength
c___ used in this subroutine we multiply by fudge = 39.946*ghe(kkl)
c___ Note also that for the hydrogenic cross-section which depends
c___ on energy only as 1/e, the excitation and deexcitation strengths
c___ are the same
c
           if (ku.gt.nbn) go to 5000
           khu=ku-3
           fudge = 39.946*ghe(kkl)
           edel = 13.605*(1./real(khl)**2-1./real(khu)**2)
           partde = fos(khl,khu)/edel*fudge
           if (khl.eq.2) partde=partde*ghe(kku)/16.
           if (khl.eq.1) partde=partde*2.
           partex = partde
           go to 1000
c     
c___ for transitions from bound to auto use mewes formula
c___ with  a=.15, b=c=0.00, d=.28,
c___            mewe astron.astrophys., v20, p215 (1972)
c___ to ensure the definition of the collsion strength we must multiply
c___ by fudge = 198.10*ghe(kkl)
c
 5000      continue

           fudge = 198.10*ghe(kkl)
           osv = oratehe(kku,kkl)*1.3474e+21*ghe(kku)/(delta*evtohz)**2
           partde = (.15+.28*log(zbyd))/delta*osv*fudge
           if (ebyd .ge. 1.) then
              partex = (.15+.28*log(ebyd))/delta*osv*fudge
           endif
           go to 1000
           
c
c___ obtain final rates
c
 1000      continue
           
           if(k-ieji.ge.1.and.e.ge.delta)then
              Vche(kkl,kku)=Vche(kkl,kku)+fe(k)*denergy(k)
     &            *7.067e-8*partex/ghe(kkl)/zeff**2
           endif
           if(k+ieji.le.nebins)then
              Vche(kku,kkl)=Vche(kku,kkl)+fe(k)*denergy(k)
     &          *7.067e-8*partde/ghe(kku)/zeff**2
           endif
 100    continue
        
 2000 continue

      return
      end
c
c
      subroutine Vclifill
 
c
c  calculate the collision rates for li-like states.
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'flystuf'
      include 'xc_stuf'

      common /mewe/ f1(14,14)
      common /cli1comR/ sd1(6,2),sd2(6,2),se1(6,2),se2(6,2),
     1                  a(7),a1(7),b(7),c(7),d(7)
      common /cli1comI/ mm(19),ll(14)

      nlv=25
      nbn=19
      ndt=14
 
      do i=1,nlv
        do j=1,nlv
            Vcli(i,j)=0.0
        enddo
      enddo
c
c  bounded detailed for excitation from 2s and 2p up to 3d using the
c  data obtained in paper by Cochrane and McWhirter;physica scripta
c  vol.28,25-44,1983 refit by rwlee to obtain cross-section data in form
c  of Mewes parametrization
c
      xc=1.40e-5
      z=atomz
      n=0
 
c
c
      do 7 i=1,2
        do 6 j=1,5
          glbygu = gli(i)/gli(j)
          deli=abs(evli(j)-evli(i))
          ieji=indxc(deli)
          if(deli.le.exthr) goto 6
          if(ieji.le.0.or.ieji.gt.nebins) goto 6
          if (i.ge.j) go to 6
          if(i.eq.1)n=j-1
          if(i.eq.2)n=j+2
          if (deli .eq. 0.00) go to 6
          do k=1,nebins
            e = energy(k)
            de = denergy(k)
            ff = max(fli(i,j),1.d0)
            if (k-ieji.ge.1.and.y .ge. 1.) then
               y = e/deli
               Vcli(i,j) = Vcli(i,j)+fe(k)*xc/deli*ff*de*
     1               (a(n)+a1(n)/(z-2)+b(n)/y+c(n)/y**2+d(n)*log(y))
            endif
            if (k+ieji.le.nebins)then
               yp = energy(k+ieji)/deli
               Vcli(j,i) = Vcli(j,i)+fe(k)*xc/deli*ff*de*glbygu*
     1              (a(n)+a1(n)/(z-2)+b(n)/yp+c(n)/yp**2+d(n)*log(yp))
            endif
          enddo
6       continue
7     continue
c
c  These collisional excitation rates for 3s-3p,3s-3d,and 3p-3d
c  are calculated using the approximation in the paper
c  of Mewe Astron. & Astrophys.20,215-221(1972).
c
c  setup

      call mewe1

c  test if 3l will be used

      do 101 i=3,4
         do 102 j=4,5
            if (i .eq. j) go to 102
            deli = abs(evli(j)-evli(i))
            glbygu = gli(i)/gli(j)
            ieji=indxc(deli)
            if (deli .le. denergy(1)) go to 102
            if (ieji.le.0.or.ieji.gt.nebins) goto 102
c     
            do k=1,nebins
               e = energy(k)
               de = denergy(k)
               if (fli(i,j) .eq. 0.) then
                  if(k+ieji.le.nebins)then
                     Vcli(j,i) = Vcli(j,i)+
     1                 fe(k)*xc*f1(i,j)/deli*.15*de*glbygu
                  endif
                  if (k-ieji.ge.1.and.e.gt.deli) then
                     Vcli(i,j) = Vcli(j,i)+fe(k)*xc*f1(i,j)/deli*.15*de
                  endif
               else
                  if (k+ieji.le.nebins)then
                     yp = energy(k+ieji)/deli
                     Vcli(j,i) = Vcli(j,i)+fe(k)*xc*fli(i,j)/deli
     1                 *(.6+.28*log(yp))*de*glbygu
                  endif
                  if (k-ieji.ge.1.and.energy(k).gt.deli) then
                     y = energy(k)/deli
                     Vcli(i,j) = Vcli(i,j)+fe(k)*xc*fli(i,j)/deli
     1                    *(.6+.28*log(y))*de
                  endif
               endif
            enddo
 102     continue
 101  continue
c
c  using approximations from Mewe for the remaining transitions
c
c  collisional rates between detailed states
c
      do 250 i=1,ndt
        do 200 jj=1,ndt-5
          j=jj+5
          if (i .ge. j) go to 200
          deli = abs(evli(j)-evli(i))
          if (deli .le. denergy(1)) go to 200
          ieji=indxc(deli)
          if (ieji.le.0.or.ieji.gt.nebins) goto 200
          glbygu = gli(i)/gli(j)

c
c  test if transition is allowed
c
          if ( fli(i,j).ne.0.) then
            do k=1,nebins
              e = energy(k)
              de = denergy(k)
              y = e/deli
              if ( mm(i)-mm(j) .eq. 0.) then
                if (k-ieji.ge.1.and.y .gt. 1.) then
                   Vcli(i,j) = Vcli(i,j)+
     1                  fe(k)*xc*fli(i,j)/deli*(.6+.28*log(y))*de
                endif
                if (k+ieji.le.nebins)then
                   yp = energy(k+ieji)/deli
                   Vcli(j,i) = Vcli(j,i)+fe(k)*xc*fli(i,j)/deli
     1                  *(.6+.28*log(yp))*de*glbygu
                endif
              else
                if (k-ieji.ge.1.and.y .gt. 1.) then
                   Vcli(i,j) = Vcli(i,j)+fe(k)*xc*fli(i,j)/deli
     1                  *(.15+.28*log(y))*de
                endif
                if (k+ieji.le.nebins)then
                   yp = energy(k+ieji)/deli
                   Vcli(j,i) = Vcli(j,i)+fe(k)*xc*fli(i,j)/deli
     1                  *(.15+.28*log(yp))*de*glbygu
                endif
             endif
          enddo
c     
c  test for delta l > 2 or more
c
       elseif (abs(ll(i)-ll(j)) .le. 2) then
          do k=1,nebins
             e = energy(k)
             de = denergy(k)
             if (k-ieji.ge.1.and.e .gt. deli) then
                Vcli(i,j) = Vcli(i,j)+fe(k)*xc*f1(i,j)/deli*0.15*de
             endif
             if(k+ieji.le.nebins)then
                Vcli(j,i) = Vcli(j,i)+fe(k)*xc*f1(i,j)/deli*0.15*de
     &               *glbygu
             endif
          enddo
       endif
 200  continue
 250  continue
c
 103  continue
c     
c  for the case of bounded detailed to bounded non-detailed
c
      do 550 i=1,ndt
         do 501 j=ndt+1,nbn
            deli = abs(evli(j)-evli(i))
            if (deli .le. denergy(1)) go to 501
            glbygu = gli(i)/gli(j)
            ieji=indxc(deli)
            if (ieji.le.0.or.ieji.gt.nebins) goto 501

            do k=1,nebins
               e = energy(k)
               de = denergy(k)
               if (k-ieji.ge.1.and.e .gt. deli) then
                  Vcli(i,j)=Vcli(i,j)+fe(k)* xc*fli(i,j)/deli*0.38*de
               endif
               if(k+ieji.le.nebins)then
                  Vcli(j,i)=Vcli(j,i)+fe(k)*xc*fli(i,j)/deli*0.38*de
     &                 *glbygu
               endif
            enddo
 501     continue
 550  continue
c
c  for the case non-detailed to non-detailed
c
      do 650 i=ndt+1,nbn
        do 601 j=ndt+1,nbn
          if (i .ge. j) go to 601
          deli = abs(evli(j)-evli(i))
          if (deli .le. denergy(1)) go to 601
          ieji=indxc(deli)
          if (ieji.le.0.or.ieji.gt.nebins) goto 601
          glbygu = gli(i)/gli(j)

          do k=1,nebins
             e = energy(k)
             de = denergy(k)
             if (k-ieji.ge.1.and.e .gt. deli) then
                Vcli(i,j)=Vcli(i,j)+fe(k)* xc*fli(i,j)/deli*0.38*de
             endif
             if( k+ieji.le.nebins)then
                Vcli(j,i)=Vcli(j,i)+fe(k)*xc*fli(i,j)/deli*0.38*de
     &               *glbygu
             endif
          enddo
 601   continue
 650  continue
c
c  auto to atuo
c
c  test if auto states are used
c
900   if ((nlv-nbn) .ne. 6) go to 950
c
      ccc = 7.067e-8/(atomz-1.5)**2
      do 750 i=nbn+1,nlv
        do 701 j=nbn+1,nlv
          if (i .ge. j) go to 701
          deli = abs(evli(j)-evli(i))
          if (deli .le. denergy(1)) go to 701
          ieji=indxc(deli)
          if (ieji.le.0.or.ieji.gt.nebins) goto 701
          mi = i-nbn
          mj = j-nbn

          do 702 k=1,nebins
            e = energy(k)
            de = denergy(k)
            y = e/deli
            partex = 0.0
            partde=0.d0

            if(k+ieji.gt.nebins) goto 702
            yp = energy(k+ieji)/deli
            go to (300,400,500,600,700), mi
  300       continue
c
c___choose upper state.
c
            go to (304,305,306,307,308), (mj-1)
c
  304       partde = xcsliaut(z,1)
            if (y .ge. 1.) partex = xcsliaut(y,1)
            go to 1000
  305       partde = xcsliaut(yp,2)
            if (y .ge. 1.) partex = xcsliaut(y,2)
            go to 1000
  306       partde = xcsliaut(yp,3)
            if (y .ge. 1.) partex = xcsliaut(y,3)
            go to 1000
  307       partde = 0.0
            go to 1000
  308       partde = xcsliau1(yp,5)
            if (y .ge. 1.) partex = xcsliau1(y,5)
            go to 1000
c
c___1s2p (singlet p) 2s doublet p (mix)
c
  400       continue
c
c___choose upper state.
c
            go to (405,406,407,408), (mj-2)
c
  405       partde = 0.0
            go to 1000
  406       partde = xcsliaut(yp,7)
            if (y .ge. 1.) partex = xcsliaut(y,7)
            go to 1000
  407       partde = xcsliaut(yp,8)
            if (y .ge. 1.) partex = xcsliaut(y,8)
            go to 1000
  408       partde = xcsliaut(yp,9)
            if (y .ge. 1.) partex = xcsliaut(y,9)
            go to 1000
c
c___1s2p (triplet p) 2s doublet p (mix)
c
  500       continue
c
c___choose upper state.
c
            go to (506,507,508), (mj-3)
c
  506       partde = xcsliaut(yp,10)
            if (y .ge. 1.) partex = xcsliaut(y,10)
            go to 1000
  507       partde = xcsliaut(yp,11)
            if (y .ge. 1.) partex = xcsliaut(y,11)
            go to 1000
  508       partde = xcsliaut(yp,12)
            if (y .ge. 1.) partex = xcsliaut(y,12)
            go to 1000
c
c___1s2p2 doublet d
c
  600       continue
c
c___choose upper state.
c
            go to (607,608), (mj-4)
c
  607       partde = 0.0
            go to 1000
  608       partde = xcsliaut(yp,14)
            if (y .ge. 1.) partex = xcsliaut(y,14)
            go to 1000
c
c___1s2p2 doublet p
c
  700       continue
c
c___only one upper state.
c
            partde = 0.0
            go to 1000
c
c
 1000       continue
            if (k-ieji.ge.1.and.y .ge. 1.) then
               Vcli(i,j)=Vcli(i,j)+fe(k)* ccc*partex*de/gli(i)
            endif
            if(k+ieji.le.nebins)then
               Vcli(j,i)=Vcli(j,i)+fe(k)* ccc*partde*de/gli(j)
            endif

 702     continue

701     continue
750   continue
 
c
c   inner to auto
c
      ccc = 7.067e-8/(atomz-1.5)**2
      do 850 i=1,2
        do 800 j=nbn+1,nlv
          m = j-nbn
          deli = abs(evli(j)-evli(i))
          ieji=indxc(deli)
          if (ieji.le.0.or.ieji.gt.nebins) goto 800
          do k=1,nebins
            e = energy(k)
            de = denergy(k)
            if(k+ieji.le.nebins)then
               yp = energy(k+ieji)/deli
               part = sd1(m,i)*xcsamd(yp,1)+sd2(m,i)*xcsamd(yp,2)
     *              +se1(m,i)*xcsamex(yp,1)+se2(m,i)*xcsamex(yp,2)
               Vcli(j,i)=Vcli(j,i)+fe(k)* part*ccc*de*gli(j)
            endif
            if (k-ieji.ge.1.and.e .ge. deli) then
               y = e/deli
               part = sd1(m,i)*xcsamd(y,1)+sd2(m,i)*xcsamd(y,2)
     *              +se1(m,i)*xcsamex(y,1)+se2(m,i)*xcsamex(y,2)
               Vcli(i,j) = Vcli(i,j) + fe(k)* part*ccc*de*gli(i)
            endif
         enddo
 800  continue
 850  continue
 
 950  continue

      return
      end

c
c           ******    ***    *     *  ******
c          *         *   *   **   **   *    *
c          *        *     *  * * * *   *    *
c           *****   *******  *  *  *   *    *
c                *  *     *  *     *   *    *
c                *  *     *  *     *   *    *
c          ******   *     *  *     *  ******
c
      function xcsamd(y,i)
c
c---- direct hydrogenic cross section - Golden etal
c                       Ap.J.Supp.Ser.,V45,p603,(1981)
c
      implicit real*8 (a-h,o-z)
      common /gcgs/ ap(11),dp(11),c1(11),c2(11),a(11),ce1(11),
     1              ce2(11),ae(11)
      common /gcgsI/ mr(11)
c
      xcsamd=dp(i)+ap(i)*log(y)+c1(i)/(a(i)+y)+c2(i)/(a(i)+y)**2

      return
      end

c
c           ******    ***    *     *  *******  *     *
c          *         *   *   **   **  *         *   *
c          *        *     *  * * * *  *          * *
c           *****   *******  *  *  *  ****        *
c                *  *     *  *     *  *          * *
c                *  *     *  *     *  *         *   *
c          ******   *     *  *     *  *******  *     *
c
      function xcsamex(y,i)
c
c----exchange hydrogenic cross-section - see golden et al
c                          ap.j.supp.ser.,v45,p603,(1981)
c
      implicit real*8 (a-h,o-z)
      common /gcgs/ ap(11),dp(11),c1(11),c2(11),a(11),ce1(11),
     1              ce2(11),ae(11)
      common /gcgsI/ mr(11)
c
      m = mr(i)
      denom = ae(i)+y
 
      xcsamex = ce1(i)/denom**(m-1)+ce2(i)/denom**m
 
      return
      end
c
c           ******  ******   *     *  *     *   *
c          *         *    *  **    *  **    *  **
c          *         *    *  * *   *  * *   *   *
c           *****    *****   *  *  *  *  *  *   *
c                *   *       *   * *  *   * *   *
c                *   *       *    **  *    **   *
c          ******    *       *     *  *     *  ***
c
      function xcspnn1(m,y)
c
c___ sampson and parks ap.j.1974 delta spin change
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      common /sc2/ a2(9),b2(9),a3(3),b3(3),d3(3),c30(17),c31(17),
     1             c32(17),d2(9)
      common /sc2I/ mci(9)
c
      xcspnn1 = b2(m)/(a2(m)+y)**mci(m)*exp(-d2(m)*y)
c
      return
      end
c
c           ******  ******   *     *  *     *   *****
c          *         *    *  **    *  **    *  *     *
c          *         *    *  * *   *  * *   *  *     *
c           *****    *****   *  *  *  *  *  *  *     *
c                *   *       *   * *  *   * *  *     *
c                *   *       *    **  *    **  *     *
c          ******    *       *     *  *     *   *****
c
      function xcspnno(m,n,y)
c
c___sampson and parks ap.j. 1974 delta l change
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'

      common /sc2/ a2(9),b2(9),a3(3),b3(3),d3(3),c30(17),c31(17),
     1             c32(17),d2(9)
      common /sc2I/ mci(9)
c
      z=atomz
      t3=c30(n)+c31(n)/z+c32(n)/(z*z)
      xcspnno = b3(m)*y/(a3(m)+y)**2+d3(m)*log(a3(m)+y)+t3
c
      return
      end

c
c          ******   *     *   ******   ******  ***    ***    *     *
c           *    *  *     *  *        *         *    *   *   **    *
c           *    *  *     *  *        *         *   *     *  * *   *
c           *****   *     *   *****    *****    *   *******  *  *  *
c           *  *    *     *        *        *   *   *     *  *   * *
c           *   *   *     *        *        *   *   *     *  *    **
c           *    *   *****   ******   ******   ***  *     *  *     *
c
c
      function xc_russian(kl,ku,ebyd)
c
c---- calculate the collisional excitation or deexcitation from
c---- helium state kl to ku using the rates from Vinogradov
c---- refit to the forms of Sampson et al. so that we can obtain
c---- cross-sections.  Note that Russian as it is here returns the
c---- "z-scaled" collision strength as defined by Sampson.  Thus,
c---- the odd definition of the factors multiplying the bracket term
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'flystuf'
      dimension bkl(8,8),bfitkl(8,8),akl(8,8),ckl(8,8)
      data akl/
     1 1*0., 11.003, 1.9231, 3.9909, 33.333, 2.0000, 2.2727, 1.6129,
     1 2*0., 18982., 25.000, 1897.8, 2.2727, 1.8519, 1.0101,
     1 3*0.,1.6871e+6,25.00, 2.2727, 1.8519, 1.0101,
     1 4*0., 5522.4, 2.0000, 2.1739, 1.4706,
     1 5*0., 2.0000, 2.1739, 1.5873,
     1 6*0., 2.1277, 1.2195,
     1 7*0., 1.4286,
     1 8*0./
      data bkl/
     1 1*0., 3.4000, 6.0000, 20.880, 20.000, 73.400, 64.000, 16.000,
     1 2*0., 0.0470, 0.5400, 0.0640, 53.400, 40.000, 30.000,
     1 3*0., 8.54e-4,0.1700, 17.800, 13.400, 10.00,
     1 4*0., 0.2130, 201.60, 138.60, 77.800,
     1 5*0., 67.000, 46.000, 26.000,
     1 6*0., 700.00, 312.00,
     1 7*0., 420.00,
     1 8*0./
      data bfitkl/
     1 1*0., 0.0000, -1.200, 0.0000, -970.0, -1.300, -1.909, -.7355,
     1 2*0., 0.0000, -480.0, 0.0000, -1.909, -1.107, -.0071,
     1 3*0., 0.0000, -480.0, -1.909, -1.107, -.0071,
     1 4*0., 0.0000, -1.300, -1.761, -0.4706,
     1 5*0., -1.300, -1.761, -.7048,
     1 6*0., -1.691, -.1756,
     1 7*0., -.4286,
     1 8*0./
      data ckl/
     1 1*0., 2.1439, 0.3000, 0.8932, 29.000, 0.3000, 0.5000, 0.2000,
     1 2*0., 129.64, 19.000, 40.301, 0.5000, 0.3000, -.3000,
     1 3*0., 1231.0, 19.000, 0.5000, 0.3000, -.3000,
     1 4*0., 69.457, 0.3000, 0.5000, 0.0000,
     1 5*0., 0.3000, 0.5000, 0.2000,
     1 6*0., 0.5000, -.2000,
     1 7*0., 0.0000,
     1 8*0./

      xc_russian=0.d0
      eu = eipz(2)-debye(iatomz-1)-evhe(ku)
      el = eipz(2)-debye(iatomz-1)-evhe(kl)
      if(eu.le.0.d0.or.el.le.0.d0)return
      ehe = evhe(ku)-evhe(kl)
      if(ehe.le.0.d0) return
c
      if (bfitkl(ku,kl) .eq. 0.0) then
        bracket = bkl(ku,kl)*akl(ku,kl)/(ckl(ku,kl)+ebyd)**2
      else
        bracket = bkl(ku,kl)*(akl(ku,kl)+bfitkl(ku,kl)/
     1                      (ckl(ku,kl)+ebyd))
      endif

      xc_russian = (eu/el)**1.5*(13.6/ehe)/2.17*(atomz-1.)**2*bracket

      return
      end

c
c
      function xcsliau1(y,m)
c
      implicit real*8 (a-h,o-z)
      common /collali/ az2(16),b(16),c(16)
c
      xcsliau1 = az2(m)+c(m)*log(y)
c
      return
      end
c
c
c
c           ******  *        ***    ***    *     *  *******
c          *        *         *    *   *   *     *     *
c          *        *         *   *     *  *     *     *
c           *****   *         *   *******  *     *     *
c                *  *         *   *     *  *     *     *
c                *  *         *   *     *  *     *     *
c          ******   *******  ***  *     *   *****      *
c
c
c
c
      function xcsliaut(y,m)
c
      implicit real*8 (a-h,o-z)
      common /collali/ az2(16),b(16),c(16)
c
      xcsliaut = az2(m)/(b(m)+y)+c(m)*(log(y)+0.6931)
c
      return
      end
c
c
c
c           *****   *     *  *     *  *******  ***  *        *
c          *     *  *     *   *   *   *         *   *        *
c          *        *     *    * *    *         *   *        *
c          *        *******     *     ****      *   *        *
c          *        *     *     *     *         *   *        *
c          *     *  *     *     *     *         *   *        *
c           *****   *     *     *     *        ***  *******  *******
c
c
c
c
      function chy(iir,ijr)
c
c___ this calculates the collisional excitation or de-excitation
c___ from level i to level j.
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'flystuf'
c
      i=min(iir,ijr)
      j=max(ijr,iir)

      factor = 3.157e-06/(evhy(j)-evhy(i))*fos(i,j)/sqrt(tev)
c
c___  excitation rate
c
       if(iir.eq.i)then
          fexp=exp((evhy(i)-evhy(j))/tev)
          chy= factor*fexp
c
c___ de-excitation rate
c
       else if(iir.eq.j) then
          chy= factor*ghy(i)/ghy(j)
       endif

       return
       end
c
c
c
c           *****   *     *  *******  *******  ***  *        *
c          *     *  *     *  *        *         *   *        *
c          *        *     *  *        *         *   *        *
c          *        *******  ****     ****      *   *        *
c          *        *     *  *        *         *   *        *
c          *     *  *     *  *        *         *   *        *
c           *****   *     *  *******  *        ***  *******  *******
c
c
c
c
      function che(iir,ijr)
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'flystuf'

      common /dataa/ aaut(5,5),z2s2(5,5),coa(5,5),c1a(5,5),c2a(5,5)
      common /idatnum/ idatsam, idatboi
c
      data idatsam/5/, idatboi/8/
c
      data z2s2/.013,0.,0.,38.8,0.,161.7,.356,0.,.0015,0.,
     1          0.,.17,.0,0.,0.,53.8,4*0.,33.1,4*0./
      data aaut/33.7,8.5,48.,4.,4.2,5.2,8.9,30.5,11.5,0.,
     1          85.95,12.87,24.5,0.,0.,4.14,11.,3*0.,4.9,4*0./
      data coa/.076,.0082,.6,59.,.312,275.7,.4075,.0747,.0095,0.,
     1         .505,.31,.086,2*0.,114.8,14.75,3*0.,53.8,4*0./
      data c1a/-12.35,-.45,28.,-141.,.40,-110.6,-1.,-33.42,-1.35,0.,
     1     -227.7,3.68,-24.0,2*0.,726.6,-80.,3*0.,-261.9,4*0./
      data c2a/8351.,30.,-1359.,3200.,-.015
     1        ,22340.,1398.,18045.,209.,0.,343200.,3135.,4217.,2*0.,
     2         13332.,608.6,3*0.0,4415.,4*0.0/
c
c___ calculate a constant.
c

      idsamp1 = idatsam+1
      iboiko = 1

      nbn0= 25
      nbn = 13
      nau = 6

c
c___ calculate the collisional rates between helium like levels i to j
c
      kki = iir ! lev.index
      kkj = ijr ! lev. index

      if(iir.gt.nbn0)kki = kki+nbn-nbn0     !kin.index
      if(ijr.gt.nbn0)kkj = kkj+nbn-nbn0     !kin.index

      if (evhe(iir) .gt. evhe(ijr)) then
         kku = iir ! lev.index
         kkl = ijr
         ku = kki  ! kin.index
         kl = kkj
      else
         kku = ijr 
         kkl = iir
         ku = kkj
         kl = kki
      endif
      eu = evhe(kku)
      el = evhe(kkl)
      tz = (eu-el)/tev
      exptz = exp(-tz)

      if(iir.gt.nbn.and.iir.le.nbn0)goto 5000
      if(ijr.gt.nbn.and.ijr.le.nbn0)goto 5000
c
c***********************************************************
c*                                                         *
c*  note:  for boiko data have all bound levels.           *
c*         for sampson data have detail for 5 levels.      *
c*                                                         *
c***********************************************************
c
c
c---- determine which set of rates to use
c
      if (ku.le.idatboi.and.iboiko.eq.1) go to 3000
c
c---- use sampsons rates - golden et.al. ap.j.supp.ser.,v45,p603,(1981)
c
      ccc=2.17e-08/(atomz-1.)**2*sqrt(13.6058/tev)
      if (kl.gt.nbn) go to 109
      if (kl.gt.idatsam) go to 106
      go to (101,102,103,104,105) , kl
c
c___ treat autoionizing level collisions differently
c
 101  continue
c
c___ kl is the 1s state
c
      khl=1
      if (ku.gt.nbn) go to 5000
      if (ku.gt.idsamp1) go to 4000
      if (ku.eq.idsamp1) go to 111
      nd=ku/2
      part=1.5*samex(tz,nd)*ccc
      if (ku-nd*2.gt.0) part=-part+2.*samd(tz,nd)*ccc
      go to 1000
 111  continue
      part=2.*(samd(tz,3)+samd(tz,4)+samd(tz,5))*ccc
      go to 1000
c
c___ 2 triplet s state
c
 102  continue
      khl=2
      if (ku.gt.idsamp1) go to 4000
      if (ku.eq.idsamp1) go to 152
      if (ku-4) 122,132,142
 122  continue
      part=.75*spnn1(5)*ccc
      go to 1000
 132  continue
      part=spnno(1,5)*ccc
      go to 1000
 142  continue
      part=.75*spnn1(1)*ccc
      go to 1000
 152  continue
      part=3.*(samd(tz,6)+samd(tz,7)+samd(tz,8))*ccc
      go to 1000
c
c___ 2 singlet s state
c
 103  continue
      khl=2
      if (ku.gt.idsamp1) go to 4000
      if (ku.eq.idsamp1) go to 113
      if (ku-4) 123,123,133
 123  continue
      part=.75*spnn1(1)*ccc
      go to 1000
 133  continue
      part=spnno(1,4)*ccc
      go to 1000
 113  continue
      part=(samd(tz,6)+samd(tz,7)+samd(tz,8))*ccc
      go to 1000
c
c___ 2 triplet p state
c
 104  continue
      khl=2
      if (ku.gt.idsamp1) go to 4000
      if (ku.eq.idsamp1) go to 114
      part=.75*spnn1(6)*ccc
      go to 1000
 114  continue
      part=3.*(samd(tz,9)+samd(tz,10)+samd(tz,11))*ccc
      go to 1000
c
c___ 2 singlet p state
c
 105  continue
      khl=2
      if (ku.gt.idsamp1) go to 4000
      if (ku.eq.idsamp1) go to 107
 107  continue
      part=(samd(tz,9)+samd(tz,10)+samd(tz,11))*ccc
      go to 1000
c     
c---- use Vinogradovs rate for n = 1,2,3,4,5 - sov.j.quan.ele.,5,p630,(1978)
c
 3000 continue
      part = russian(kl,ku)
      go to 1000
c
c___ levels treated as hydrogenic n .ge. 3
c
 106  continue
      khl=kl-3
      go to 4000
c
c___ use sampsons data for collisions between auto-ionizing levels
c___            goett etal. atomic and nuc. data tables v25, p185 (1980)
c
 109  continue
      kal = kl-nbn
      kau = ku-kl
      expauttz = exp(aaut(kau,kal)*tz)
      expe = exptz/expauttz
      expa = aaut(kau,kal)*tz+tz
      part=ccc*(coa(kau,kal)*exptz+1.66*z2s2(kau,kal)*e1(tz,exptz)+
     1     tz*expauttz*(c1a(kau,kal)*e1(expa,expe)
     2     +c2a(kau,kal)/(aaut(kau,kal)+1.)*(expe-expa*e1(expa,expe))))
      go to 1000
      
 4000 continue
c
c___ use hydrogenic rates for bound levels above n = 3
c
      if (ku.gt.nbn) go to 5000
      khu=ku-3
      part=2.32e-07/sqrt(tev)*fos(khl,khu)/(atomz-1.)**2*ghe(kkl)
     1     *real(khu**2*khl**2)/(real(khu**2)-real(khl**2))
      if (khl.eq.2) part=part*ghe(kku)/16.
      if (khl.eq.1) part=part*2.
      go to 1000
c
c___ for transitions from bound to auto use mewes formula
c___ with  a=.15, b=c=0.00, d=.28,
c___            mewe astron.astrophys., v20, p215 (1972)
c
 5000 continue
      osv = oratehe(kku,kkl)*1.3474e+21*ghe(kku)/((eu-el)*evtohz)**2
      part = (.15+.28*xexpe1(tz)/tz)*1.578e-5/sqrt(tev)/(eu-el)*osv
      go to 1000
 
c
c___ obtain final rates
c     
 1000 continue

      if(iir.eq.kkl)then
c     
c___ excitation
c
         che = part*exptz/ghe(kkl)

      else if(iir.eq.kku) then
c
c___ de-excitation
c
         che = part/ghe(kku)
      endif
         
      return
      end
c
c
      function spnn1(m)
c
c___ sampson and parks ap.j.1974 delta spin change
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'flystuf'
      common /sc2/ a2(9),b2(9),a3(3),b3(3),d3(3),c30(17),c31(17),
     1             c32(17),d2(9)
      common /sc2I/ mci(9)
c
      x0=atomz**2*13.6058/tev
      z20=a2(m)*(x0+d2(m))
      xexpe1z2 = xexpe1(z20)
      mcinow = mci(m)
 
      if (mcinow .eq. 1) then
         top=b2(m)*x0/z20*xexpe1z2
         bot=a2(m)
      else
         top=b2(m)*x0*(1.-xexpe1z2)
         bot=a2(m)**2
      endif
      spnn1=top/bot
c
      return
      end
c
c           ******  ******   *     *  *     *   *****
c          *         *    *  **    *  **    *  *     *
c          *         *    *  * *   *  * *   *  *     *
c           *****    *****   *  *  *  *  *  *  *     *
c                *   *       *   * *  *   * *  *     *
c                *   *       *    **  *    **  *     *
c          ******    *       *     *  *     *   *****
c
      function spnno(m,n)
c
c___sampson and parks ap.j. 1974 delta l change
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'flystuf'
      common /sc2/ a2(9),b2(9),a3(3),b3(3),d3(3),c30(17),c31(17),
     1             c32(17),d2(9)
      common /sc2I/ mci(9)
c
      data a2/.335,.125,.12,.114,.27,.131,.24,.094,.05/
      data b2/.438,.07177,.1045,.3531,.6499,4.,.7239,2.377,5.92/
      data d2/0.,0.,0.,0.,0.,.17,0.,0.,.07/
      data mci/2,2,2,2,1,1,1,1,1/
      data a3/1.,1.,.75/
      data b3/480.,3000.,3000./
      data d3/72.,432.,540./
      data c30/214.3,1872.,2436.,165.7,193.4,1862.,2149.,3841.,2669.,
     1166.2,160.7,185.6,172.1,1641.,1239.,1612.,1919./
      data c31/-1744.,-12530.,-9132.,167.,-643.4,-5675.,-12330.,
     1-19540.,-9746.,-1796.,-1460.,-1231.,-1630.,-9040.,-6763.,-9795.,
     2-9290./
      data c32/3249.,27040.,15560.,-3169.,-409.4,4623.,29240.,57090.,
     116470.,2660.,1294.,574.2,1894.,10390.,1663.,10970.,7311./
c
      z=atomz
      z30=a3(m)*z*z*13.6058/tev
      b30x0=b3(m)*z*z*13.6058/tev
      xexpe1z3 = xexpe1(z30)
 
      t1=((b30x0+d3(m))/z30*xexpe1z3-b30x0*(1.-xexpe1z3))
      t2=d3(m)*log(a3(m))
      t3=c30(n)+c31(n)/z+c32(n)/(z*z)
      spnno=t1+t2+t3
c
      return
      end
c
      function russian(kl,ku)
c
c---- calculate the collisional excitation or deexcitation from
c---- helium state kl to ku using the rates from Vinogradov
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'flystuf'

      dimension bkl(8,8),xkl(8,8),phikl(8,8)
      data bkl/
     1 1*0.,3.4,6.0,20.88,20.,73.4,64.,16.,
     1 2*0.,.047,.54,.064,53.4,40.,30.,
     1 3*0.,8.54e-4,.17,17.8,13.4,10.,
     1 4*0.,.213,201.6,138.6,77.8,
     1 5*0.,67.,46.,26.,
     1 6*0.,700.,312.,
     1 7*0.,420.,
     1 8*0./
      data xkl/
     1 1*0.,.38,.52,.63,.03,.5,.44,.62,
     1 2*0.,.0008,.04,.029,.44,.54,.99,
     1 3*0.,.001,.04,.44,.54,.99,
     1 4*0.,.017,.5,.46,.68,
     1 5*0.,.5,.46,.63,
     1 6*0.,.47,.82,
     1 7*0.,.70,
     1 8*0./
      data phikl/
     1 1*0.,0.,1.,0.,1.,1.,1.,1.,
     1 2*0.,0.,1.,0.,1., 1.,1.,
     1 3*0.,0.,1.,1.,1.,1.,
     1 4*0.,0.,1.,1.,1.,
     1 5*0.,1.,1.,1.,
     1 6*0.,1.,1.,
     1 7*0.,1.,
     1 8*0./
c 
      eu = eipz(2)-debye(iatomz-1)-evhe(ku)
      el = eipz(2)-debye(iatomz-1)-evhe(kl)
      ehe = evhe(ku)-evhe(kl)
      ehet = ehe/tev
      pnow = sqrt(ehet)*(ehet+phikl(ku,kl))
c
      russian = 1.00e-8*(eu/el*13.6/ehe)**1.5*bkl(ku,kl)
     1        *pnow/(ehet+xkl(ku,kl))
      return
      end
c
c
      function samd(y,i)
c
c---- direct hydrogenic cross section - Golden etal
c                       Ap.J.Supp.Ser.,V45,p603,(1981)
c
      implicit real*8 (a-h,o-z)
      common /gcgs/ ap(11),dp(11),c1(11),c2(11),a(11),ce1(11),
     1              ce2(11),ae(11)
      common /gcgsI/ mr(11)
c
      data ap/0.000,2.220,0.000,.3560,0.000,0.000,
     1        12.52,0.000,1.174,0.000,60.12/
      data dp/.4439,.4474,.0881,.1986,.0678,5.316,
     1        -8.445,16.79,-.7230,17.46,5.517/
      data c1/-.2198,-1.256,-.04442,-.3481,-.09661,-.1898,
     1        14.76,5.845,.6093,-5.130,-10.99/
      data c2/4.5e-5,5.103,.04206,.736,.07271,-8.47,
     1       21.03,-469.1,5.758,7.413,311.1/
      data a/1.67,.70,.68,.43,-.14,2.10,
     1       .99,6.32,1.06,.54,1.42/
      data ce1/6.836e-3,1.056,1.414e-3,.2259,.01438,.05082,
     1         2.794,11.10,3.318,17.62,25.70/
      data ce2/.5483,3.212,.1090,.5057,.08622,1.21,
     1        13.26,11.66,10.68,22.55,85.96/
      data ae/.9,.8,.62,.54,.7,.75,1.68,1.2,2.,1.3,1.22/
      data mr/2,3,2,3,3,2,3,3,3,3,3/
c
      ay=a(i)*y
      ay1=ay+y
      xexpe1y = xexpe1(y)
      xexpe1ay = xexpe1(ay1)
      t1=dp(i)+ap(i)/y*xexpe1y
      t2=1./(1.+a(i))*(c1(i)*xexpe1ay+c2(i)*y*(1.-xexpe1ay))
c
      samd=t1+t2
      return
      end
c
c
c
      function samex(y,i)
c
c----exchange hydrogenic cross-section - see golden et al
c                          ap.j.supp.ser.,v45,p603,(1981)
c
      implicit real*8 (a-h,o-z)
      common /gcgs/ ap(11),dp(11),c1(11),c2(11),a(11),ce1(11),
     1              ce2(11),ae(11)
      common /gcgsI/ mr(11)
c
      ay1=(ae(i)+1.)*y
      mrnow = mr(i)
      xexpe1ay = xexpe1(ay1)
 
      if (mrnow .eq. 2) then
        t1=ce1(i)*xexpe1ay
        t2=ce2(i)*y*(1.-xexpe1ay)
      else
        t1=ce1(i)*y*(1.-xexpe1ay)
        t2=ce2(i)/(ae(i)+1.)*y*((1.-ay1)+ay1*xexpe1ay)*.5
      endif
 
c
      samex = (t1+t2)/(ae(i)+1.)
      return
      end
c
c
c
c           *****   *        ***  *******  ***  *        *
c          *     *  *         *   *         *   *        *
c          *        *         *   *         *   *        *
c          *        *         *   ****      *   *        *
c          *        *         *   *         *   *        *
c          *     *  *         *   *         *   *        *
c           *****   *******  ***  *        ***  *******  *******
c
c
c
c
c
      function cli(iir,ijr)
 
c
c  calculate the collision rates for li-like states.
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'flystuf'

      common /mewe/ f1(14,14)
      common /cli1comR/ sd1(6,2),sd2(6,2),se1(6,2),se2(6,2),
     1                  a(7),a1(7),b(7),c(7),d(7)
      common /cli1comI/ mm(19),ll(14)
c
      data a / 0.60,0.0450,-0.244, 0.151,-0.35,0.0485,-0.034/
      data a1/-1.50,0.0000, 0.000,-0.090, 0.20,0.0000, 0.000/
      data b / 0.27,0.0045, 0.250, 0.650, 0.25,0.0000, 0.238/
      data c / 5.00,0.6000, 4.000, 3.600, 7.50,0.0000, 3.000/
      data d / 0.00,3.5000, 0.000, 4.000, 0.00,0.0000, 0.000/
c
c  mm- principal quantum no. and ll- azimuthal quantum no.
c
      data sd1/1.80,4*0.0,.196,0.0,1.2,10.8,3*0.0/
      data sd2/0.0,0.5,3.6,3*0.0,.065,2*0.0,3.33,6.0,.601/
      data se1/6*0.0,0.0,1.8,-7.83,3*0.0/
      data se2/0.0,0.6,-2.61,3*0.0,4*0.0,-4.0,0.0/
 
 
      nlv=25
      nbn=19
      ndt=14
 
c
c  bounded detailed for excitation from 2s and 2p up to 3d using the
c  data obtained in paper by Cochrane and McWhirter;physica scripta
c  vol.28,25-44,1983
c  setup

      call mewe1

c
      xc=1.58e-5/sqrt(tev)
      z=atomz
c
      i=min(iir,ijr)
      j=max(iir,ijr)
      deli=abs(evli(j)-evli(i))
      y=deli/tev
      if (deli.eq.0.00.or.y.gt.450.d0)then
         cli=0.d0
         return
      endif

      cliup=0.d0
      if(j.le.5) then
         if(i.le.2) then
            if(i.eq.1)n=j-1
            if(i.eq.2)n=j+2
            if (fli(i,j).eq.0.) then
               ff=1.0
               cliup=abs(xc/deli*ff*exp(-y)
     1         *(a(n)+a1(n)/(z-2)+b(n)*log((1.+y*c(n))/(1.+y*d(n)))))
            else
               ff=fli(i,j)
               cliup=abs(xc/deli*ff*exp(-y)
     1         *(a(n)+a1(n)/(z-2)+b(n)*log((1.+y*c(n))/y )))
            endif
c
c  These collisional excitation rates for 3s-3p,3s-3d,and 3p-3d
c  are calculated using the approximation in the paper
c  of Mewe Astron. & Astrophys.20,215-221(1972).
c
         else if(i.le.4) then 
            if (fli(i,j) .eq. 0.) then
               cliup = xc*f1(i,j)/deli*exp(-y)*.15
            else
               cliup = xc*fli(i,j)/deli*exp(-y)
     1              *(.6+.28*(log((y+1.)/y)-0.4/(y+1.)**2))
            endif
         endif
c
c  using approximations from Mewe for the remaining transitions
c  collisional rates between detailed states
c
      else 

c     for the bounded detailed to bounded detailed
         if(i.le.ndt.and.j.le.ndt) then
            if ( fli(i,j).eq.0.) then
               if (abs(ll(i)-ll(j)).gt.2.) then
                  cliup = 0.0
               else
                  cliup = xc*f1(i,j)/deli*exp(-y)*0.15
               endif
            else
               if ( mm(i)-mm(j) .eq. 0.) then
                  cliup = xc*fli(i,j)/deli*exp(-y)
     1                 *(.6+.28*( log((y+1.)/y)-0.4/(y+1.)**2))
               else
                  cliup = xc*fli(i,j)/deli*exp(-y)
     1                 *(.15+.28*( log((y+1.)/y)-0.4/(y+1.)**2))
               endif
            endif
c
c  for the case of bounded detailed to bounded non-detailed
c
         else if(i.le.ndt.and.j.le.nbn) then
            cliup = xc*fli(i,j)/deli*exp(-y)*0.38
c
c  for the case non-detailed to non-detailed
c
         else if(j.le.nbn) then
            cliup = xc*fli(i,j)/deli*exp(-y)*0.38
c
c  auto to atuo
c
         else if(i.gt.nbn) then
            ccc = 2.17e-08*sqrt(13.605/tev)/(atomz-1.5)**2
            mi = i-nbn
            mj = j-nbn
            go to (300,400,500,600,700), mi
 300        continue
c
c___choose upper state.
c
            go to (304,305,306,307,308), (mj-1)
c
 304        part = sliaut(y,1)
            go to 1000
 305        part = sliaut(y,2)
            go to 1000
 306        part = sliaut(y,3)
            go to 1000
 307        part = 0.0
            go to 1000
 308        part = sliau1(y,5)
            go to 1000
c     
c___  1s2p (singlet p) 2s doublet p (mix)
c
 400        continue
c
c___choose upper state.
c
            go to (405,406,407,408), (mj-2)
c
 405        part = 0.0
            go to 1000
 406        part = sliaut(y,7)
            go to 1000
 407        part = sliaut(y,8)
            go to 1000
 408        part = sliaut(y,9)
            go to 1000
c
c___1s2p (triplet p) 2s doublet p (mix)
c
 500        continue
c
c___choose upper state.
c
            go to (506,507,508), (mj-3)
c
 506        part = sliaut(y,10)
            go to 1000
 507        part = sliaut(y,11)
            go to 1000
 508        part = sliaut(y,12)
            go to 1000
c
c___1s2p2 doublet d
c
 600        continue
c
c___choose upper state.
c
            go to (607,608), (mj-4)
c
 607        part = 0.0
            go to 1000
 608        part = sliaut(y,14)
            go to 1000
c
c___1s2p2 doublet p
c     
 700        continue
c
c___only one upper state.
c
            part = 0.0
            go to 1000
c
c
 1000       continue
            cliup=exp(-y)*part*ccc/gli(i)
c
c   inner to auto
c
         else if(i.le.2.and.j.gt.nbn) then
            m = j-nbn
            part= sd1(m,i)*samd(y,1)+sd2(m,i)*samd(y,2)
     *           +se1(m,i)*samex(y,1)+se2(m,i)*samex(y,2)
            ccc=2.17e-08*sqrt(13.605/tev)/(atomz-1.5)**2
            cliup=exp(-y)*part*ccc/gli(i)
         endif
      endif

c
c  now to defind deexcitational rate form detailed balance.
c
      if(iir.eq.i) then
         cli=cliup
      else if(iir.eq.j) then
         cli= gli(i)/gli(j)*exp(y)*cliup
      endif
c
      return
      end
 
c
c           ******  *        ***    ***    *     *   *
c          *        *         *    *   *   *     *  **
c          *        *         *   *     *  *     *   *
c           *****   *         *   *******  *     *   *
c                *  *         *   *     *  *     *   *
c                *  *         *   *     *  *     *   *
c          ******   *******  ***  *     *   *****   ***
c
c
c
c
      function sliau1 (y,m)
c
      implicit real*8 (a-h,o-z)
      common /collali/ az2(16),b(16),c(16)
c
      sliau1 = az2(m)+c(m)*xexpe1(y)/y
c
      return
      end
c
c
c
c           ******  *        ***    ***    *     *  *******
c          *        *         *    *   *   *     *     *
c          *        *         *   *     *  *     *     *
c           *****   *         *   *******  *     *     *
c                *  *         *   *     *  *     *     *
c                *  *         *   *     *  *     *     *
c          ******   *******  ***  *     *   *****      *
c
c
c
c
      function sliaut (y,m)
c
      implicit real*8 (a-h,o-z)
      common /collali/ az2(16),b(16),c(16)
c
      data az2/ 11.57, 235.4, 1935., 00.00, 0.865, 00.00, 50.26,744.6
     1         ,15.51, 2101.1,37.92, 1168., 0.00, 13280., 0.00, 0.0 /
      data b/ 0.81, .449, 340.825, 00.00, 00.00, 00.00, 1.29, 1.089,
     1        .718, 1.89, 1.64, 2.948, 00.00, 305.13, 00.00, 00.00 /
      data c/ 3.717, 105.7, 00.00, 00.00, -.050, 00.00, 11.58, 197.5,
     1        5.545, 329.1, 6.945, 141.9, 00.00, 00.00, 00.00, 00.00 /
c
c
      sliaut = az2(m)*xexpe1(y*(b(m)+1.))/(b(m)+1.)
     1         +c(m)*(xexpe1(y)/y+.6931)
c
      return
      end


c
c          *     *  *     *  *        *******  *     *  *******  *
c          *     *   *   *   *        *        *     *  *        *
c          *     *    * *    *        *         *   *   *        *
c          *******     *     *        ****      *   *   ****     *
c          *     *     *     *        *          * *    *        *
c          *     *     *     *        *          * *    *        *
c          *     *     *     *******  *******     *     *******  *******
c
c
c
c
      subroutine hylevel(nbn,nau,elv,glv,knam,nocion,morb0,gamrt)
c
c*****************************
c*                           *
c*  d e c l a r a t i o n s  *
c*                           *
c*****************************
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'flystuf'

      dimension elv(70),glv(70),nocion(70,nprim+1),morb0(70)
      dimension gamrt(70)
      character*8 knam(70)

      common /jhs/ e2p32(26),e2tp(26),e2sp(26),e2ts(26),e2ss(26)
      common /scofrh/ rhsc(26)
cc
c___ data from J.H. Scofield's tabulation LLNL Report
c             - UCID 16848  (July 1975)
c
      data rhsc/5*1.,3952056.8,5380081.2,7028385.,
     1               8897233. ,10986862.,13297660.,15829923.,18584114.,
     2               21560599.,24759891.,28182510.,31828989.,35699878.,
     3               39795768.,44117777.,48665353.,53440499.,58443646.,
     4               63675611.,69137168.,74829016./
      data  namhy /  'hy1    ',  'hy2    ',  'hy3    ',  'hy4    ',
     1               'hy5    ',  'hy6    ',  'hy7    ',  'hy8    ',
     2               'hy9    ',  'hy10   ',  'hy11   ',  'hy12   ',
     3               'hy13   ',  'hy14   ',  'hy15   ',  'hy16   ',
     4               'hy17   ',  'hy18   ',  'hy19   ',  'hy20   ',
     5               'hy21   ',  'hy22   ',  'hy23   ',  'hy24   ',
     6               'hy25   ' /
      data namhe /  'he1s   ',  'he2st  ',  'he2ss  ',  'he2pt  ',
     1               'he2ps  ',  'he3ps  ',  'he4ps  ',  'he5ps  ',
     2               'he6ps  ',  'he7ps  ',  'he8ps  ',  'he9ps  ',
     3               'he10ps ',  'he11ps ',  'he12ps ',  'he13ps ',
     4               'he14ps ',  'he15ps ',  'he16ps ',  'he17ps ',
     5               'he18ps ',  'he19ps ',  'he20ps ',  'he21ps ',
     6               'he22ps ',  '2s2s1s ',  '2s2p3p ',  '2p2p3p ',
     7               '2p2p1d ',  '2s2p1p ',  '2p2p1s ' /
      data namli/'li2s','li2p','li3s','li3p','li3d',
     1     'li4s','li4p','li4d','li4f',
     2     'li5s','li5p','li5d','li5f','li5g',
     3     'li6','li7' ,'li8' ,'li9' ,'li10',
     4     'op' ,'qr'  ,'st'  ,'jkl' ,'abcd','mn'/
cc
c*******************
c*  begin routine  *
c*******************
c
c
c___calculate energies and statistics.
c
      zt = max(2.0,atomz)
      rhz = rh/(1.0+1.0/(zt*atomz*1836.42))
      rhz = atomz*atomz*rhz
      if (iatomz .gt. 5) rhz=rhsc(iatomz)
      ccc = rhz*1.23981e-04
      eipz(1) = ccc
c
      nbn = 25
      nau= 0
      nlv = 10
      evhy(1) = 0.0
      ghy(1)  = 2.0
      evhy(2) = e2p32(iatomz)
      ghy(2)  = 8.0
      if (nbn .lt. 3) go to 200
      do 100 i=3,nbn
         nupper = i
         evhy(i) = ccc*(1.0-1.0/real(nupper)**2)
         ghy(i)  = 2.0*i*i
  100 continue
c
      do j=1,nbn
         elv(j)=evhy(j)
         glv(j)=ghy(j)
         knam(j)=namhy(j)
         do ii=1,nprim
            nocion(j,ii)=0
         enddo
         morb0(j)=j
         if(j.le.nprim)then
            nocion(j,j)=1
         else
            nocion(j,nprim+1)=j
         endif
      enddo
c
  200 continue

      gamrt(1)=0.d0
      do 300 j=2,nbn
         jm1 = j-1
         gamrt(j)=0.d0
         do 301 i=1,jm1
            oratehy(j,i) = 
     &           4.371e+07*(evhy(j)-evhy(i))**2*ghy(i)/ghy(j)*fos(i,j)
            gamrt(j)=gamrt(j)+oratehy(j,i)
  301    continue
  300 continue
c

      return
      end






c
c
c
c          *     *  *******  *        *******  *     *  *******  *
c          *     *  *        *        *        *     *  *        *
c          *     *  *        *        *         *   *   *        *
c          *******  ****     *        ****      *   *   ****     *
c          *     *  *        *        *          * *    *        *
c          *     *  *        *        *          * *    *        *
c          *     *  *******  *******  *******     *     *******  *******
c
c
c
c
      subroutine helevel(nbn,nau,elv,glv,knam,nocion,morb0,gamrt)
c
c*****************************
c*                           *
c*  d e c l a r a t i o n s  *
c*                           *
c*****************************
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'flystuf'

      dimension elv(70),glv(70),nocion(70,nprim+1),morb0(70)
      dimension gamrt(70)
      character*8 knam(70)

c
c____ Rates for the He-like n=2 states taken from GWF Drake Phys.Rev.
c             ,A3, p900,(1971) and Ap.J.,V163, p439 (1070)
c     Rates for the n>2 states taken from A. Bromage SRC Report
c             AL-R3 (1978)
c
      common /aautR/ etran(6,6),ztran(6)
      common /aautI/ nzlow(6)
      common /ggg/ gtemp(11)
      common /henergy/ wip(26)
      common /he34567/ w3p(26),w4p(26),w5p(26),w6p(26),w7p(26)

      common /jhs/ e2p32(26),e2tp(26),e2sp(26),e2ts(26),e2ss(26)

c
c___ data from J.H. Scofield's tabulation LLNL Report
c             - UCID 16848  (July 1975)
c
      data e2p32/ 0.0, 40.820, 91.8406, 163.07, 255.0,
     &          367.5, 500.3, 653.7,   827.5,  1021.9, 1236.9,
     &   1472.6,1728.9,2006.0,2303.9,2622.6, 2962.3,
     &   3322.9,3704.6,4107.4,4531.4,4976.7, 5443.4,
     &   5931.6,6441.4,6972.9/

      data e2ts/0.0,19.9819,59.0191,118.6,198.556,
     1          299.0,419.8,561.1,730.6,905.1,1107.8,1331.1,
     1          1575.0,1839.5,2124.6,2430.4,2756.9,3104.2,3472.3,
     2          3861.3,4271.2,4702.0,5154.0,5627.0,6121.2,6636.7/

      data e2ss/0.0,20.61523,60.9209,121.66,202.771,
     1          304.4,426.4,569.0,731.9,915.3,1119.3,1343.9,
     1          1589.0,1854.7,2141.1,2448.2,2776.0,3124.6,3494.0,
     2          3884.3,4295.5,4727.8,5181.1,5655.6,6151.2,6668.2/

      data e2tp/0.0,20.9636,61.2792,121.93,202.943,
     1          304.41,426.3,568.7,731.5,914.8,1118.7,1343.1,
     1          1588.1,1853.8,2140.1,2447.1,2775.0,3123.6,3493.0,
     2          3883.4,4315.4,4727.0,5180.6,5655.0,6150.8,6667.8/

      data e2sp/0.0,21.21748,62.2145,123.68,205.563,
     1          307.9,430.7,574.0,737.7,922.0,1126.8,1352.3,
     1          1598.3,1865.0,2152.5,2460.7,2789.7,3139.6,3510.5,
     2          3902.4,4294.7,4749.7,5205.2,5682.1,6180.5,6700.4/

      data w3p/ 0.0,0.186203,0.56175,1.13232,1.89818,
     1          2.85940, 4.01600, 5.36820, 6.91620,
     1          8.66030, 10.6007, 12.7377, 15.0718, 17.6032,
     2          20.3324, 23.2599, 26.3860, 29.7114, 33.2367,
     3          36.9623, 40.8890, 45.0175, 49.3485, 53.8827,
     4          58.6210, 63.5641/

      data w4p/ 0.0,0.191486,0.58283,1.179830,1.98275,
     1          2.99180, 4.20680, 5.62810, 7.25590,
     1          9.07050, 11.1321, 13.3811, 15.8378, 18.5027,
     2          21.3761, 24.4584, 27.7503, 31.2522, 34.9647,
     3          38.8884, 43.0241, 47.3722, 51.9337, 56.7093,
     4          61.6998, 66.9061/

      data w5p/0.0,0.193936,0.59263,1.201894,2.022,
     1         3.0529, 4.2950, 5.7485, 7.4132,
     1         9.2901,11.3784,13.6795,16.1928,18.9198,
     2        21.8598,25.0135,28.3826,31.9662,35.7652,
     3        39.7811,44.0132,48.4631,53.1307,58.0186,
     4        63.1258,68.4532/
      data w6p/0.0,0.195269,0.59788,1.213931,2.04336,
     1         3.0868, 4.3434, 5.8138, 7.4987,
     1         9.3982,11.5122,13.8416,16.3856,19.1457,
     2        22.1227,25.3152,28.7262,32.3541,36.1999,
     3         40.2650,44.5504,49.0551,53.7816,58.7292,
     4         63.9001,69.2937/
      data w7p/0.0,0.196073,0.60112,1.221135,2.0553,
     1         3.1061, 4.3724, 5.8533, 7.5504,
     1         9.4635,11.5929,13.9392,16.5017,19.2828,
     2         22.2808,25.4975,28.9327,32.5872,36.4620,
     3        40.5578,44.8738,49.4124,54.1736,59.1574,
     4        64.3671,69.8002/

c
      data wip/ 0.0,24.5876,75.638,153.893,259.368,
     1          392.08185,552.,739.34866,953.89,1195.7859,
     1          1465.09, 1761.8, 2085.7580, 2437.3323, 2816.4570,
     2          3223.1824, 3657.5695, 4119.6954,  4609.6, 5127.4,
     3          5673.2, 6246.9, 6843.9, 7479.0, 8137.5, 8824.4/
c
c___data for the auto-ionizing states from Vainshtein and Safronova
c       Atomic Data and Nuclear Data Tables, V21,49, (1978)
c
      data etran/9.878,9.934,9.977,10.013,10.042,10.068,
     1      10.074711,10.097775,10.127759,10.148617,10.168034,
     1      10.186273,
     1      10.06082,10.091267,10.116584,10.138633,10.158857,10.178247,
     1      10.036385,10.071438,10.100272,10.125705,10.149159,10.17190,
     1      10.114155,10.13736,10.1575,10.17591,10.19376,10.2114,
     1      10.187278,10.19751,10.20833,10.21971,10.23192,10.24519/
      data nzlow/5,2,4,5,3,5/
      data ztran/10.,12.,14.,16.,18.,20./
      data gtemp/1.,3.,1.,9.,3.,1.,9.,9.,5.,3.,1./

c*******************
c*  begin routine  *
c*******************
c
c___ load in the energy levels a values auto-ionization rates
c___ for helium like ions
c
c___ fill in the energy levels
c
      rzm1z4 = (azm1/atomz)**4
      evhe(1) = 0.0
c
c___ use Scofield data
c
      ncal = 5
      evhe(5) = e2sp(iatomz)
      evhe(4) = e2tp(iatomz)
      evhe(3) = e2ss(iatomz)
      evhe(2) = e2ts(iatomz)

      do 101 i=1,ncal
         ghe(i)=gtemp(i)
  101 continue

c---- place the n = 3 thru n =  7 energies thru weights

      evhe(6) = w3p(iatomz)*123.981
      ghe(6) = 36.
      evhe(7) = w4p(iatomz)*123.981
      ghe(7) = 64.
      evhe(8) = w5p(iatomz)*123.981
      ghe(8) = 100.
      evhe(9) = w6p(iatomz)*123.981
      ghe(9) = 144.
      evhe(10) = w7p(iatomz)*123.981
      ghe(10) = 196.

      nlv = 19
      nbn = 25
      nau = 6
      nlast = nbn+nau
      nb  = nbn-ncal
      if (nb .gt. 5) then
c
         do 12 i=6,nb
            if (iatomz .eq. 2) then
              nq = 4
            else if (iatomz .eq. 3) then
              nq = 1
            else
              nq = 0
            endif
            evhe(i+ncal)=wip(iatomz)*(1.-1./real(i+2+nq)**2)
            ghe(i+ncal)=4.*real(i+2)**2
   12    continue

         if (iatomz .eq. 2) then
           evhe(11) = 24.375
           evhe(12) = 24.419
           evhe(13) = 24.451
           evhe(14) = 24.475
           evhe(15) = 24.492
         endif

      endif

c___ put in auto-ionization level energies
c
      if (nau.ne.6) go to 100
      naum1 = nau-1
      do 13 i=1,naum1
         if (atomz.ge.ztran(i).and.atomz.le.ztran(i+1)) go to 14
   13 continue
      do 150 j=1,nau
         ghe(nbn+j)=gtemp(5+j)
         if (atomz .gt. ztran(nau)) then
            evhe(nbn+j)=atomz**2*etran(6,j)+evhe(nzlow(j))
         else
            if (iatomz .eq. 2) then
              evhe(nbn+1) = 57.87
              evhe(nbn+2) = 58.31
              evhe(nbn+3) = 59.67
              evhe(nbn+4) = 59.88
              evhe(nbn+5) = 60.13
              evhe(nbn+6) = 62.14
            else
              evhe(nbn+j)=atomz**2*etran(1,j)+evhe(nzlow(j))
            endif
         endif
  150 continue
      go to 100
   14 ii = i
c
      do 15 j=1,nau
         ghe(nbn+j)=gtemp(5+j)
         evhe(nbn+j)=atomz**2*((atomz-ztran(ii))*etran(ii+1,j)+
     1                         (ztran(ii+1)-atomz)*etran(ii,j))/
     2                         (ztran(ii+1)-ztran(ii))+evhe(nzlow(j))
   15 continue

c
c___ fill in the a rates and auto-ionization rates
c
  100 continue

      do j=1, nlast
         elv(j)=evhe(j)
         glv(j)=ghe(j)
         knam(j)=namhe(j)
         do ii=1,nprim
            nocion(j,ii)=0
         enddo
         nocion(j,1)=1
         morb0(j)=1
         if(j.eq.1)then
            nocion(j,1)=2
         else if(j.ge.2.and.j.le.5) then
            nocion(j,2)=1
            morb0(j)=2
         else if(j.ge.6.and.j.le.nbn)then
            morb0(j)=j-3
            if(morb0(j).le.nprim)then
               nocion(j,j-3)=1
            else
               nocion(j,nprim+1)=morb0(j)
            endif
         else if(j.gt.nbn)then
            nocion(j,1)=0
            nocion(j,2)=2
            morb0(j)=2
         endif
      enddo

      do 16 i=1,nlast
         do 16 j=1,nlast
            oratehe(i,j)=0.0
   16 continue
c
c____ The next set of avalues are taken from detailed calculations
c         of DWF Drake and are used in this code as averaged rates
c
      zn=atomz
      zn10=zn/10.
      z10sq=zn10**2

c____ 1s2s Triplet S - 1s1s Singlet S

      if (zn .eq. 2.) then
        oratehe(2,1) = 0.0
      else
        oratehe(2,1)=1.09e4*(zn10**10)*(1.524-.585/zn10+.061/z10sq)
       endif

c____ 1s2s Singlet S - 1s1s Singlet S

      if (zn .eq. 2.) then
        oratehe(3,1) = 0.0
      else
        oratehe(3,1)=1.0e7*(zn10**6)*(1.53-.59/zn10+.06/z10sq)
       endif

c_____ 1s2p Triplet P (j=1 & 2) - 1s1s singlet S

      if (zn .eq. 2.) then
        oratehe(4,1) = 0.0
      else
        a41j2=2.27e6*(zn10**8)*(1.637-.718/zn10+.081/z10sq)
        a41j1=5.43e9*(zn10**10)*(1.22-.23/zn10+.01/z10sq)
        oratehe(4,1)=(5.*a41j2+3.*a41j1)/9.
      endif

c_____ 1s2p Triplet P (j=1 & 2) - 1s2s Triplet S

      if (zn .eq. 2.) then
        oratehe(4,2) = 1.022e+07
      else
        fj2 = .685+.315/zn10
        fj1 = .91+.09/zn10
        if (zn10 .le. 1.) go to 31
           fj2 = 1.
           fj1 = 1.
        if (zn.le.15.) go to 31
           fj2 = .08+z10sq*.39
           fj1 = .6+zn10*.27
   31   continue
        a42j2=1.11e8*fj2*(zn10**1.666667)
        a42j1=1.06e8*fj1*(zn10**1.33333333)
        oratehe(4,2)=(5.*a42j2+3.*a42j1)/9.
      endif

c____ 1s2p Singlet P - 1s1s Singlet S

      if (zn .eq. 2.) then
        oratehe(5,1) = 1.799e+09
      else
        oratehe(5,1)=8.87e12*(zn10**4)*(1.435-.4785/zn10+.0435/z10sq)
      endif

c____ 1s2p Singlet P - 1s2s Singlet S
                          
      if (zn .eq. 2.) then
        oratehe(5,3) = 1.976e+06
      else
        oratehe(5,3)=3.20e5*(zn10**7)*(1.20-.205/zn10+.0083/z10sq)
      endif

c____ the n = 3 state is lumped .Thus detailed rates must be averaged

      if (zn .eq. 2) then
        oratehe(6,1) = 5.66e+08*3./36.
        oratehe(6,2) = 9.478e+06*9./36.
        oratehe(6,3) = 1.338e+07*3./36.

        oratehe(6,4) = 2.78e+07*3./36.
        oratehe(6,5) = 1.81e+07*1./36.

      else
        oratehe(6,1)=1.80e12*(zn10**4.47)*3./36.
        oratehe(6,2)=1.09e11*(zn10**4.60)*9./36.
        oratehe(6,3)=1.14e11*(zn10**4.51)*3./36.

        oratehe(6,4)=1.72e10*(zn10**4.23)*3./36.
        oratehe(6,4)=2.14e11*(zn10**4.41)*15./36.+oratehe(6,4)
        oratehe(6,4)=7.00e9*(zn10**6.24)*5./36.+oratehe(6,4)

        oratehe(6,5)=8.36e9*(zn10**6.27)*15./36.
        oratehe(6,5)=3.36e11*(zn10**4.35)*5./36.+oratehe(6,5)
      endif
c
c___ first the a values from singly to singly
c
      if (zn .eq. 2.) then
        oratehe(7,1) = 3./(4.*16.)*2.46e+8
        oratehe(8,1) = 3./(4.*25.)*1.28e+8
        oratehe(9,1) = 3./(4.*36.)*0.719e+8
        oratehe(10,1) = 3./(4.*49.)*0.507e+8
        oratehe(11,1) = 3./(4.*64.)*0.343e+8
        oratehe(12,1) = 3./(4.*81.)*0.237e+8
        oratehe(13,1) = 3./(4.*100.)*0.181e+8
        oratehe(14,1) = 3./(4.*121.)*0.130e+8
        oratehe(15,1) = 3./(4.*144.)*0.104e+8

c   use the fact that the A value scales as 1/(n**3 x delta energy)

        avalue = 0.104e+8*12.**3*evhe(15)
        do i=16,nbn-3
          oratehe(i,1) = 3./4.*avalue/((i-3.)**4*(evhe(i)-evhe(1)))
        enddo

        oratehe(7,2) = 9./(4.*16.)*0.0505e+08
        oratehe(8,2) = 9./(4.*25.)*0.0293e+08
        oratehe(9,2) = 9./(4.*36.)*0.0169e+08
        oratehe(10,2) = 9./(4.*49.)*0.0111e+08
        oratehe(11,2) = 9./(4.*64.)*0.0078e+08
        oratehe(12,2) = 9./(4.*81.)*0.0055e+08
        oratehe(13,2) = 9./(4.*100.)*0.00404e+08

c   use the fact that the A value scales as 1/(n**3 x delta energy)

        avalue = 0.00404e+8*10.**3*(evhe(13)-evhe(2))
        do i=14,nbn-3
          oratehe(i,2) = 9./4.*avalue/((i-3.)**4*(evhe(i)-evhe(2)))
        enddo

        oratehe(7,3) = 3./(4.*16.)*0.0717e+08
        oratehe(8,3) = 3./(4.*25.)*0.0376e+08
        oratehe(9,3) = 3./(4.*36.)*0.0239e+08
        oratehe(10,3) = 3./(4.*49.)*0.0130e+08
        oratehe(11,3) = 3./(4.*64.)*0.00901e+08
        oratehe(12,3) = 3./(4.*81.)*0.0065e+08
        oratehe(13,3) = 3./(4.*100.)*0.00490e+08

c   use the fact that the A value scales as 1/(n**3 x delta energy)

        avalue = 0.00490e+8*10.**3*(evhe(13)-evhe(3))
        do i=14,nbn-3
          oratehe(i,3) = 3./4.*avalue/((i-3.)**4*(evhe(i)-evhe(3)))
        enddo

        oratehe(7,4) = (3.*0.106e+08+15.*0.251e+08)/(4.*16.)
        oratehe(8,4) = (3.*0.0430e+08+15.*0.117e+08)/(4.*25.)
        oratehe(9,4) = (3.*0.0236e+08+15.*0.0589e+08)/(4.*36.)
        oratehe(10,4) = (3.*0.0150e+08+15.*0.0444e+08)/(4.*49.)
        oratehe(11,4) = (3.*0.0108e+08+15.*0.0261e+08)/(4.*64.)
        oratehe(12,4) = (3.*0.00800e+08+15.*0.0205e+08)/(4.*81.)
        oratehe(13,4) = (3.*0.00543e+08+15.*0.0131e+08)/(4.*100.)

c   use the fact that the A value scales as 1/(n**3 x delta energy)

        avalue =(3.*.00543e+8+15.*.0131e+8)/18.*1000.*(evhe(13)-evhe(4))
        do i=14,nbn-3
          oratehe(i,4) = 18./4.*avalue/((i-3.)**4*(evhe(i)-evhe(4)))
        enddo

        oratehe(7,5) = (0.0655e+08+5.*0.202e+08)/(4.*16.)
        oratehe(8,5) = (0.0313e+08+5.*0.0907e+08)/(4.*25.)
        oratehe(9,5) = (0.0176e+08+5.*0.0495e+08)/(4.*36.)
        oratehe(10,5) = (0.0109e+08+5.*0.0279e+08)/(4.*49.)
        oratehe(11,5) = (0.00718e+08+5.*0.0195e+08)/(4.*64.)
        oratehe(12,5) = (0.00400e+08+5.*0.0126e+08)/(4.*81.)
        oratehe(13,5) = (0.00303e+08+5.*0.00971e+08)/(4.*100.)

c   use the fact that the A value scales as 1/(n**3 x delta energy)

        avalue = (.00303e08+5.*.00971e+8)/6.*1000.*(evhe(13)-evhe(5))
        do i=14,nbn-3
          oratehe(i,5) = 6./4.*avalue/((i-3.)**4*(evhe(i)-evhe(5)))
        enddo
c
        do i=8,nbn
          iloup = i-1
          do j=7,iloup
            oratehe(i,j) = rzm1z4*oratehy(i-3,j-3)
          enddo
        enddo

      else
        nboo=nbn-3
        do 18 i=4,nboo
          oratehe(3+i,1)=rzm1z4*oratehy(i,1)*.5
          oratehe(3+i,2)=rzm1z4*oratehy(i,2)*3./16.
          oratehe(3+i,3)=rzm1z4*oratehy(i,2)/16.
          oratehe(3+i,4)=rzm1z4*oratehy(i,2)*9./16.
          oratehe(3+i,5)=rzm1z4*oratehy(i,2)*3./16.
          oratehe(3+i,6)=rzm1z4*oratehy(i,3)
   18   continue
c
        do 20 i=8,nbn
          iloup = i-1
          do 20 j=7,iloup
            oratehe(i,j) = rzm1z4*oratehy(i-3,j-3)
   20   continue
      endif
c
c___ fill in a rates for doubly excited levels
c
      oratehe(nbn+1,5)=.24*zn10**3.96*1.00e+13
      oratehe(nbn+4,5)=1.16*zn10**4.05*1.00e+13
      oratehe(nbn+6,5)=.900*zn10**4.11*1.00e+13
      oratehe(nbn+2,2)=.570*zn10**4.08*1.00e+13
      oratehe(nbn+5,3)=.590*zn10**4.06*1.00e+13
      oratehe(nbn+3,4)=1.15*zn10**4.07*1.00e+13
c
c___ auto-ionization rates
c
      
     
      do 19 i=1,nbn
         do j=1,10
            autohe(i,j)=0.0
         enddo
   19 continue
c
      autohe(nbn+1,1)=34.00e+13
      autohe(nbn+2,1)=1.36e+13
      autohe(nbn+3,1)=.0434*zn10**5.78*1.00e+13
      autohe(nbn+4,1)=37.40e+13
      autohe(nbn+5,1)=20.0e+13
      autohe(nbn+6,1)=1.65*zn10**.45*1.00e+13
c

      do 40 j=2,nbn+6
         jm1 = j-1
         gamrt(j)=0.d0
         do 41 i=1,jm1
            gamrt(j)=gamrt(j)+oratehe(j,i)
 41      continue
         if(j.gt.nbn)then
            gamrt(j)=gamrt(j)+autohe(j,1)
         endif
 40   continue

      return
      end
c
c
c
c          *        ***  *        *******  *     *  *******  *
c          *         *   *        *        *     *  *        *
c          *         *   *        *         *   *   *        *
c          *         *   *        ****      *   *   ****     *
c          *         *   *        *          * *    *        *
c          *         *   *        *          * *    *        *
c          *******  ***  *******  *******     *     *******  *******
c
c
c
c
c
      subroutine lilevel(nbn,nau,elv,glv,knam,nocion,morb0,gamrt)
c
c*****************************
c*                           *
c*  d e c l a r a t i o n s  *
c*                           *
c*****************************
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'flystuf'

      dimension elv(70),glv(70),nocion(70,nprim+1),morb0(70)
      dimension gamrt(70)
      character*8 knam(70)

c
c____ Lithium-like energy levels taken from Edlen, Physica Scripta
c              V19, 255 (1979); Vainshtein And Safronova Physica
c              Scripta, V31, 519 (1985)
c____ Lithium-like A values taken from Martin and Wiese,
c              Phys. Rev. V13,699 (1976) & J. Phys. Chem. Ref.
c              Data, V5, #3 (1976)
c
c____ Auto-ionizing energy levels from A.H. Gabriel MNRAS, V160, 99, (1972)
c     Autoionizing rates from Vainshtein and Safronova Atomic Data
c              and Nuclear data tables , V21, 49 (1978)
c
c____ All other rates calculated using the Klapisch Suite of Codes at LLNL
c              provided by B. Whitten and checked with scaling by L. Petway
c
      common /liauto/ etran(11,6),autor(7,6),arli(7,6)
      common/liauto2R/ ztran(6),ztran2(11)
      common/liauto2I/ nzlow(6)
      common /ggl/ gtemp(25)
      common /gbitli/ ggli(40)
      dimension ztem(11)
      data gtemp/2.,6.,2.,6.,10.,2.,6.,10.,14.,2.,6.,10.,14.,18.,
     1           72.,98.,128.,162.,200.,
     1           4.,6.,6.,14.,12.,4./
      data etran/ 55.42,283.84, 401.23, 539.05, 876.19, 1295.86,
     1                  1798.38,2382.88,3050.70,3803.10, 6563.29,
     2            55.42,299.98, 421.56, 563.04, 908.28, 1335.29,
     3                  1845.50,2437.70,3114.31,3873.20, 6664.19,
     4            55.42,303.95, 425.76, 567.94, 914.31, 1342.52,
     5                  1853.23,2447.32,3123.72,3885.34, 6676.39,
     6            52.45,298.39, 419.28, 560.74, 904.31, 1330.27,
     7                  1838.66,2430.05,3104.17,3863.54, 6650.25,
     8            52.08,299.61, 420.70, 562.27, 906.96, 1333.70,
     9                  1843.04,2435.78,3110.41,3872.00, 6657.77,
     +            52.08,305.37, 420.70, 569.76, 916.34, 1344.70,
     1                  1856.00,2450.22,3127.66,3890.22, 6683.59/
      data ztran/ 6., 8., 10., 14., 18., 26. /
      data ztran2/3.,6.,7.,8.,10.,12.,14.,16.,18.,20.,26./
      data nzlow/ 2, 1, 1, 2, 2, 2 /

c___ rates x e-13

      data autor/ 7.03, 8.97, 10.3, 11.7, 12.60, 13.0, 13.9,
     1            .478, .576, .585, .660, .7700, .828, 1.33,
     2            4.80, 6.21, 7.23, 8.32, 8.900, 9.26, 9.51,
     3            6.67, 9.46, 11.1, 13.0, 13.83, 14.2, 13.4,
     4           .0013,.0053, .017, .114, .4000, .685, 1.71,
     5            .495, .921, 1.19, 1.51, 1.700, 1.79, 1.98/
      data arli/ .0019, .0085, .0210, .1030, .3139, .4785, 1.570,
     1           .0930, .2950, .7200, 3.483, 10.23, 15.60, 44.60,
     2           .0100, .0320, .0780, .4240, 1.481, 2.257, 8.474,
     3           .0363, .1150, .2810, 1.388, 4.268, 6.505, 19.10,
     4           .0769, .2430, .5930, 2.900, 8.493, 12.95, 37.18,
     5           .0236, .0745, .1818, .8860, 2.759, 4.205, 13.82/
c

c*******************
c*  begin program  *
c*******************
c
c

      z = atomz
      nlv = 25
      nbn = 19
      nau = 6
      nlast = nbn+nau
      ndt=14


c__statistical weights
c
c  for bounded detailed state
c
      do 7 i=1,ndt
         gli(i) = gtemp(i)
    7 continue
c
c  for non detailed states
c
      if (ndt .eq. 0) then
         ggli(1)=1.0
         n=1
         goto 40
      else if (ndt .eq. 2) then
          ggli(1)=gli(1)/(gli(1)+gli(2))
          ggli(2)=gli(2)/(gli(1)+gli(2))
         n=3
         goto 40
      else if (ndt .eq. 5) then
          ggli(1)=gli(1)/(gli(1)+gli(2))
          ggli(2)=gli(2)/(gli(1)+gli(2))
          ggli(3)=gli(3)/(gli(3)+gli(4)+gli(5))
          ggli(4)=gli(4)/(gli(3)+gli(4)+gli(5))
          ggli(5)=gli(5)/(gli(3)+gli(4)+gli(5))
         n=4
         goto 40
      else if (ndt .eq. 9) then
          ggli(1)=gli(1)/(gli(1)+gli(2))
          ggli(2)=gli(2)/(gli(1)+gli(2))
          ggli(3)=gli(3)/(gli(3)+gli(4)+gli(5))
          ggli(4)=gli(4)/(gli(3)+gli(4)+gli(5))
          ggli(5)=gli(5)/(gli(3)+gli(4)+gli(5))
          ggli(6)=gli(6)/(gli(6)+gli(7)+gli(8)+gli(9))
          ggli(7)=gli(7)/(gli(6)+gli(7)+gli(8)+gli(9))
          ggli(8)=gli(8)/(gli(6)+gli(7)+gli(8)+gli(9))
          ggli(9)=gli(9)/(gli(6)+gli(7)+gli(8)+gli(9))
         n=5
         goto 40
      else if (ndt .eq. 14) then
          ggli(1)=gli(1)/(gli(1)+gli(2))
          ggli(2)=gli(2)/(gli(1)+gli(2))
          ggli(3)=gli(3)/(gli(3)+gli(4)+gli(5))
          ggli(4)=gli(4)/(gli(3)+gli(4)+gli(5))
          ggli(5)=gli(5)/(gli(3)+gli(4)+gli(5))
          ggli(6)=gli(6)/(gli(6)+gli(7)+gli(8)+gli(9))
          ggli(7)=gli(7)/(gli(6)+gli(7)+gli(8)+gli(9))
          ggli(8)=gli(8)/(gli(6)+gli(7)+gli(8)+gli(9))
          ggli(9)=gli(9)/(gli(6)+gli(7)+gli(8)+gli(9))
          ggli(10)=gli(10)/(gli(10)+gli(11)+gli(12)+gli(13)+gli(14))
          ggli(11)=gli(11)/(gli(10)+gli(11)+gli(12)+gli(13)+gli(14))
          ggli(12)=gli(12)/(gli(10)+gli(11)+gli(12)+gli(13)+gli(14))
          ggli(13)=gli(13)/(gli(10)+gli(11)+gli(12)+gli(13)+gli(14))
          ggli(14)=gli(14)/(gli(10)+gli(11)+gli(12)+gli(13)+gli(14))
         n=6
         goto 40
      endif
c
40      continue
c
      m=0
      k=n+nbn-ndt-1
      do 45 i=n,k
         m=m+1
         gli(ndt+m) = 2.*real(i)**2
         ggli(ndt+m)=1.0
45      continue
c
c  for auto states
c
      do 50 i=1,nau
          gli(nbn+i) = gtemp(19+i)
50      continue
c
c  calculate the enrgies
c
      call energyli(z)

c
c  setup the array for the oscillator strengths
c
      call fvalueli(z)
c
c---- calculate a-values between bounded states-------
c
      do 20 i=1,nlv
         do 22 j= 1,nlv
            orateli(i,j)=0.0
22         continue
20      continue
c
      confacto=1.0e-8/1.2399e-4
      const=6.67e+15
c
      do 26 i=1,nbn
         do 28 j=1,nbn
            if (i.eq.j) goto 28
            xx=(evli(j)-evli(i))*confacto
            xlambda=1.0/xx
            orateli(j,i)=const/xlambda**2*fli(i,j)*gli(i)/gli(j)
28         continue
26      continue
c
c___now fill autoionizing levels.
c
      if (nau .eq. 0) goto 111
      if (nau .ne. 6) call exit(1)
      if (z.lt.ztran(1))  then
         ii = 2
      else if (z .gt. ztran(6)) then
         ii = 6
      else
        do 1 i=2,6
           if (z .le. ztran(i)) go to 2
    1   continue
        call exit(1)
    2   ii = i
      endif
c
      do 9 i=2,11
        if(z .le. ztran2(i)) go to 10
   9  continue
      call exit(1)
   10 ii2=i
c
c___loop over auto levels.
c
      zlog=log(z)
      do 8 j=1,11
       ztem(j)=log(ztran2(j))
   8  continue
      do 3 i=1,nau
        eterp=(log(etran(ii2,i))*(zlog-ztem(ii2-1))+
     1   log(etran(ii2-1,i))*(ztem(ii2)-zlog))/
     1   (ztem(ii2)-ztem(ii2-1))

        evli(i+nbn) = evli(nzlow(i)) + exp(eterp)
    3 continue
c
c___fill autoionizing rates.
c
      do 4 i=1,nau
         do j=1, 30
            autoli(i+nbn,j)=0.d0
         enddo
         autoli(i+nbn,1) = autor(ii,i)*(z/ztran(ii))**
     1               (log10(autor(ii-1,i)/autor(ii,i))/
     2                log10(ztran(ii-1)/ztran(ii)))*1.00e+13
    4 continue
c
c___only have rates from auto-levels down.
c
      do 6 i=1,nau
         orateli(nbn+i,nzlow(i)) = arli(ii,i)*(z/ztran(ii))**
     1                           (log10(arli(ii-1,i)/arli(ii,i))/
     2                            log10(ztran(ii-1)/ztran(ii)))*1.e+13
    6 continue
c
 111  continue


      do j=1,nbn+nau
         elv(j)=evli(j)
         glv(j)=gli(j)
         knam(j)=namli(j)
         do ii=1,nprim
            nocion(j,ii)=0
         enddo
         nocion(j,1)=2
         if(j.le.2)then
            nocion(j,2)=1
            morb0(j)=2
         else if(j.le.5)then
            nocion(j,3)=1
            morb0(j)=3
         else if(j.le.9)then
            nocion(j,4)=1
            morb0(j)=4
         else if(j.le.14)then
            nocion(j,5)=1
            morb0(j)=5
         else if(j.le.nbn)then
            morb0(j)=j-9
            if(morb0(j).le.nprim)then
               nocion(j,j-9)=1
            else
               nocion(j,nprim+1)=morb0(j)
            endif
         else
            nocion(j,1)=1
            nocion(j,2)=2
            morb0(j)=2
         endif
      enddo

      do 30 j=2,nbn+nau
         jm1 = j-1
         gamrt(j)=0.d0
         do 31 i=1,jm1
            gamrt(j)=gamrt(j)+orateli(j,i)
 31      continue
         if(j.gt.nbn)then
            gamrt(j)=gamrt(j)+autoli(j,1)
         endif
 30   continue

      return
      end
 
      subroutine energyli(z)
 
c   This routine will calculate the energies for the Li-like ions
c   using the paper of Edlen  Physica Scripta 19, 255-266, 1979
c   the calculations are good for the levels 2 through 4 and for
c   elements up through z = 28. However the formulas for the ns
c   and np will work for all n. While the formulae should extra-
c   polate to higher z. Everthing is calculated in wavenumbers and
c   at the end converted to eV.
 
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'flystuf'

      common /ritz/ sa(26),pa(26)
      common /polar/ pnl(6),qnl(6)
      common /polarnl/ ppnl(12,12),qqnl(12,12)
      common /nstar/ t3s(28),t4s(28),t3p(28),t4p(28)
      common/en5sp/ d5s12(28),d5p12(28),d5p32(28)

      data pnl/ 60.2112, 28.5768, 4.0824, 15.3838, 2.3410, 0.55738/
      data qnl/ 0.22222, 0.27083, 0.02083, 0.29101, 0.02667, 0.00533/

      data t3s/ 0., 0., 0.14838, 0.13363, 0.12731, 0.12378, 0.12154,
     1         0.11998, 0.11884, 0.11798, 0.11730, 0.11676, 0.11632,
     2         0.11596, 0.11566, 0.11541, 0.11520, 0.11502, 0.11487,
     3         0.11475, 0.11464, 0.11455, 0.11448, 0.11443, 0.11438,
     4         0.11435, 0.11433, 0.11432/
 
      data t4s/ 0., 0., 0.07724, 0.07158, 0.06909, 0.06768, 0.06678,
     1         0.06615, 0.06568, 0.06533, 0.06506, 0.06484, 0.06466,
     2         0.06451, 0.06439, 0.06429, 0.06420, 0.06413, 0.06407,
     3         0.06402, 0.06397, 0.06394, 0.06391, 0.06389, 0.06387,
     4         0.06386, 0.06385, 0.06385/
 
      data t3p/ 0., 0., 0.11448, 0.11480, 0.11439, 0.11398, 0.11363,
     1         0.11335, 0.11313, 0.11295, 0.11281, 0.11270, 0.11260,
     2         0.11252, 0.11246, 0.11241, 0.11237, 0.11233, 0.11231,
     3         0.11229, 0.11228, 0.11227, 0.11227, 0.11227, 0.11228,
     4         0.11229, 0.11230, 0.11232/
 
      data t4p/ 0., 0., 0.06395, 0.06407, 0.06389, 0.06371, 0.06356,
     1         0.06345, 0.06335, 0.06328, 0.06322, 0.06317, 0.06313,
     2         0.06310, 0.06307, 0.06305, 0.06304, 0.06302, 0.06302,
     3         0.06302, 0.06301, 0.06300, 0.06300, 0.06301, 0.06301,
     4         0.06302, 0.06302, 0.06303/
 
      data d5s12/ 0.,0.,38300.,127335.,263156.,
     +            445221., 673721., 948521., 1269660., 1637210.,
     1           2051230., 2511820., 3019070., 3573100., 4174020.,
     2           4821980., 5517100., 6259560., 7049510., 7887130.,
     3           8772610., 9706140., 10688200., 11718300.,
     4           12797300., 13925300., 15102500., 16329300./

      data d5p12/ 0.,0.,39016.,128971.58,265719.7,
     +            448870., 678331., 954088., 1276190., 1644690.,
     1           2059670., 2521220., 3024440., 3584410., 4186330.,
     2           4835260., 5531370., 6274810., 7065750., 7904360.,
     3           8790840., 9725370., 10708200., 11739500.,
     4           12819600., 13948600., 15126900., 16354700./

      data d5p32/ 0.,0.,39016.1,128972.02,265720.2,
     +            448877., 678348., 954123., 1276257., 1644798.,
     1           2059840., 2521478., 3024814., 3584937., 4187052.,
     2           4836228., 5532641., 6276451., 7067836., 7906977.,
     3           8794084., 9729348., 10713032., 11745316.,
     4           12826546., 13956834., 15136595., 16366040./
 
      data sa/ 5*0.0,  .154325, .128559, .110367, .096886, .086542,
     1         .078397, .071856, .066523, .062127, .058473, .055419,
     2         .052857, .050708, .048908, .047407, .046166, .045153,
     3         .044340, .043708, .043237, .042913/
 
      data pa/ 5*0.0,  .038494, .033756, .030015, .027063, .024720,
     1         .022849, .021350, .020148, .019187, .018424, .017828,
     2         .017372, .017037, .016808, .016670, .016615, .016632,
     3         .016715, .016859, .017056, .017305/
c
      data ppnl/  14*0.0,60.2112,28.5768,15.3838,9.1392,5.8450,3.9548,
     1        2.7964, 2.0484, 1.5444 , 1.1928, 3*0.,4.0824, 2.3410,
     2        1.4336, .93273, 0.63788, .45427, .33443, 0.25308,
     3        .19600, 4*0.0 , 0.55738, .35840, .23931, 0.16623,
     4        .11960, .08867, 0.06745, .05244, 5*0.0 , 0.11404,
     5        .07914, .05620, 0.04101, .03070, .02351, 0.01837,
     6        6*0.0 , .03044, 0.02230, .01658, .01257, 0.00971,
     7        .00764, 7*0.0 , 0.00981, .00749, .00577, 0.00451,
     8        .00357, 8*0.0 , 0.00364, .00287, .00227, 0.00182,
     9        9*0.0 , .00151, 0.00122, .00099, 10*0.0, 0.00068,
     a        .00056, 11*0.0, 0.00056, 12*0.0/
 
      data qqnl/  14*0.0,.22222,.27083,.291010,.301470,.307620,.311550,
     1        .314220,.316110,.317510,.318560,3*0.,.020806,.026667,
     2        .029514, .031141, .032165, .032854, .033340, .033696,
     3        .033965, 4*0.0  , .005333, .006944, .007820, .008358,
     4        .008713, .008962, .009142, .009278, 5*0.0  , .001984,
     5        .002592, .002949, .003180, .003339, .003453, .003539,
     6        6*0.0    , .000907, .001181, .001351, .001466, .001548,
     7        .001609, 7*0.0  , .000473, .000612, .000703, .000766,
     8        .000813, 8*0.0  , .000271, .000348, .000400, .000434,
     9        9*0.0    , .000167, .000212, .000244, 10*0.0 , .000108,
     a        .000136, 11*0.0 , .000136, 12*0.0/ 
 
      dimension p5l1(4),p5l2(4),wnumber(19)
 
      nbn=19
      ndt=14
 
c   note that 's' is used as the screening factor, alf = fine structure
 
      alf = 1./(137.036)
 
c    calculate the Rydberg constant
 
      rz = 109737.318-60.200/(2.*z)
 
c   define the index for the atomic number iz
 
      if (z .le. 28) then
         iz = z
      else
         iz = 28
      endif
 
c   calculate the splitting of the 2p doublet P
 
      s = 1.7415+0.633/(z-0.80)
      zms = z-s
      split2p = 0.366076*zms**4+1.2183e-5*zms**6+4.3e-10*zms**8
 
c   calculate the relativistic energy contribution to 2s doublet S
 
      s = 1.2808
      zms = z-s
      rel2s = 0.456536*zms**4+1.2763e-5*zms**6+4.34e-10*zms**8
 
c   calculate the relativistic energy contribution to 2p doublet P
 
      s = 2.241
      zms = z-s
      rel2p = 0.21305*zms**4+0.466e-5*zms**6+1.48e-10*zms**8
 
c   calcluate the QED correction to 2s doublet S
 
      s = 1.60
      zms = z-s
      za = zms*alf
      el = log(za)
      qed2s = 0.0045246*zms**4*(-2.179-2.0*el+7.214*za-
     1        za**2*(3.*el**2+8.695*el+19.081))
 
c   calcluate the QED correction to 2p doublet P
 
      s = 2.
      zms = z-s
      za = zms*alf
      el = log(za)
      qed2p = .004525*zms**4*(0.0300-0.496*za**2*el)
 
c   calculate difference between 2S and center of gravity of 2P
 
      w2s2p = 0.141441*rz*(z-1.7025-0.768371/(z-0.975)*
     1        (1.-0.090333/(z-2.184))) - 0.21*(z-2.)**2+
     2        rel2s-rel2p-qed2s+qed2p
 
c   calculate the energies of the 2P 1/2 and 3/2
 
      w2p12 = w2s2p-2./3.*split2p
      w2p32 = w2s2p+1./3.*split2p
 
c   calculate the polarization energies to the 3d, 4d and 4f
 
      s = 0.3397+0.102/(z-0.4)
      a = 9.*(z-2.)**4/(z-s)**4
      ak = sqrt(a)*(0.2113*z+0.598-2.4/z)
 
      pol3d = a*pnl(1)*(1.+ak*qnl(1))
      pol4d = a*pnl(2)*(1.+ak*qnl(2))
      pol4f = a*pnl(3)*(1.+ak*qnl(3))
 
c  calculate the calcualted ionization potential
 
      wic = 0.25*rz*(z**2-3.18244*z+2.0038+0.20801/(z-1.3833))+
     1      rel2s-qed2s
 
 
c    calculate the TH formulae for the states
 
      s = 2.
      zms = z-2.
 
      th3d52 = rz*zms**2/9.+0.018036*zms**4+0.534e-7*zms**6
      th3d32 = rz*zms**2/9.+0.054108*zms**4+3.868e-7*zms**6
 
      th4d52 = rz*zms**2/16.+0.013315*zms**4
      th4d32 = rz*zms**2/16.+0.028533*zms**4
 
      th4f72 = rz*zms**2/16.+0.00571*zms**4
      th4f52 = rz*zms**2/16.+0.01331*zms**4
 
c   calculate the formulae for the nS and nP states using Ritz formula
 
      th3s12 = rz*zms**2*t3s(iz)
      th4s12 = rz*zms**2*t4s(iz)
 
      th3pcg = rz*zms**2*t3p(iz)
      th4pcg = rz*zms**2*t4p(iz)
 
c   calculate the splitting of the 3p and 4p levels
 
      th2pcg = wic-w2s2p
 
      split3p = 1.01*(th3pcg/th2pcg)**1.5*split2p
      split4p = 1.01*(th4pcg/th2pcg)**1.5*split2p
 
c   This is an add-on to include n=5 states; 5s12,5p12,
c   & 5p32 are taken directly from L.A. Vainshtein, et al,
c   Physica Scripta 31, 519-532,1985; the rest, 5nl are
c   calulated from formulas using the paper of Edlen
c   Physica Scripta 19, 255-266, 1979
 
 
c   CALCULATE NEEDED PARAMETERS
 
      xs=0.3397+0.102/(z-0.4)
      xaa=9.0*((z-2.)/(z-xs))**4
      xkk=(0.2113*z+0.598-2.4/z)*sqrt(xaa)
      xn=5.00
      xa=(1./137.036)**2
      y11=rz*((z-2.)/xn)**2
      y22=((z-2.)/xn)**4
 
      do 400 j=1,4
            p5l1(j)=0.
            p5l2(j)=0.
400      continue
 
 
 
      do 500 i=1,3
            xi=i
            xl=1.+xi
            y2=rz*xa*(xn/(xl+1.)-3./4.)
            y1=rz*xa*(xn/xl-3./4.)
            p1=xaa*pnl(i+3)*(1.+xkk*qnl(i+3))
 
            p5l1(i+1)=y11+y22*y1+p1
            p5l2(i+1)=y11+y22*y2+p1
 
500      continue
 
c   calculate the energies
 
      w3s12 = wic-th3s12
      w3p12 = wic-(th3pcg+2./3.*split3p)
      w3p32 = wic-(th3pcg-1./3.*split3p)
      w3d32 = wic-th3d32-pol3d
      w3d52 = wic-th3d52-pol3d
      w4s12 = wic-th4s12
      w4p12 = wic-(th4pcg+2./3.*split4p)
      w4p32 = wic-(th4pcg-1./3.*split4p)
      w4d32 = wic-th4d32-pol4d
      w4d52 = wic-th4d52-pol4d
      w4f52 = wic-th4f52-pol4f
      w4f72 = wic-th4f72-pol4f
      w5s12 = d5s12(iz)
      w5p12 = d5p12(iz)
      w5p32 = d5p32(iz)
      w5d32 = wic-p5l1(2)
      w5d52 = wic-p5l2(2)
      w5f52 = wic-p5l1(3)
      w5f72 = wic-p5l2(3)
      w5g72 = wic-p5l1(4)
      w5g92 = wic-p5l2(4)
c
c set up the statistical weights
c
      dg1=2. *  .5 + 1.
      dg3=2. * 1.5 + 1.
      dg5=2. * 2.5 + 1.
      dg7=2. * 3.5 + 1.
      dg9=2. * 4.5 + 1.
c
      ds=2.*(2. * 0.+1.)
      dp=2.*(2. * 1.+1.)
      dd=2.*(2. * 2.+1.)
      df=2.*(2. * 3.+1.)
      dg=2.*(2. * 4.+1.)
 
c calculate the energy for the multiplet
 
      w2s = 0.0
      w2p = (w2p12*dg1+w2p32*dg3)/(dg1+dg3)
      w3s = w3s12
      w3p = (w3p12*dg1+w3p32*dg3)/(dg1+dg3)
      w3d = (w3d32*dg3+w3d52*dg5)/(dg3+dg5)
      w4s = w4s12
      w4p = (w4p12*dg1+w4p32*dg3)/(dg1+dg3)
      w4d = (w4d32*dg3+w4d52*dg5)/(dg3+dg5)
      w4f = (w4f52*dg5+w4f72*dg7)/(dg5+dg7)
      w5s = w5s12
      w5p = (w5p12*dg1+w5p32*dg3)/(dg1+dg3)
      w5d = (w5d32*dg3+w5d52*dg5)/(dg3+dg5)
      w5f = (w5f52*dg5+w5f72*dg7)/(dg5+dg7)
      w5g = (w5g72*dg7+w5g92*dg9)/(dg7+dg9)
c
c____setup the array for the energies as if all possible detailed
c      levels will be used
c
      wnumber(1) = w2s
      wnumber(2) = w2p
      wnumber(3) = w3s
      wnumber(4) = w3p
      wnumber(5) = w3d
      wnumber(6) = w4s
      wnumber(7) = w4p
      wnumber(8) = w4d
      wnumber(9) = w4f
      wnumber(10) = w5s
      wnumber(11) = w5p
      wnumber(12) = w5d
      wnumber(13) = w5f
      wnumber(14) = w5g
 
c
c____ calculate the energy for the non-detailed levels
c
      if (ndt .eq. 0) then
 
         wnumber(1)=0.0
         wnumber(2)=(w3s*ds+w3p*dp+w3d*dd)/(ds+dp+dd)
         wnumber(3)=(w4s*ds+w4p*dp+w4d*dd+w4f*df)/(ds+dp+dd+df)
         wnumber(4)=(w5s*ds+w5p*dp+w5d*dd+w5f*df
     1             +w5g*dg)/(ds+dp+dd+df+dg)
         n=4
c
      else if (ndt .eq. 2) then
         wnumber(3)=(w3s*ds+w3p*dp+w3d*dd)/(ds+dp+dd)
         wnumber(4)=(w4s*ds+w4p*dp+w4d*dd+w4f*df)/(ds+dp+dd+df)
         wnumber(5)=(w5s*ds+w5p*dp+w5d*dd+w5f*df
     1             +w5g*dg)/(ds+dp+dd+df+dg)
         n=3
c
      else if (ndt .eq. 5) then
         wnumber(6)=(w4s*ds+w4p*dp+w4d*dd+w4f*df)/(ds+dp+dd+df)
         wnumber(7)=(w5s*ds+w5p*dp+w5d*dd+w5f*df
     1            +w5g*dg)/(ds+dp+dd+df+dg)
         n=2
 
      else if (ndt .eq. 9) then
         wnumber(10)=(w5s*ds+w5p*dp+w5d*dd+w5f*df
     1              +w5g*dg)/(ds+dp+dd+df+dg)
         n=1
c
      else if (ndt .eq. 14) then

         n=0
c
      endif
c
      mmm=nbn
      k = 6+nbn-ndt-1
      do 45 i=6,k
       n=n+1
       xn = i
       y11 = rz*((z-2.)/xn)**2
       y22 = ((z-2.)/xn)**4
       p1 = xaa*ppnl(i,2)*(1.0+xkk*qqnl(i,2))
       y2 = rz*xa*(xn/3.-.75)
       wnumber(ndt+n) = wic-y11-y22*y2-p1
45    continue
 
 
c   now convert the results to the units of eV
c
      do 1 i=1,mmm
         evli(i) = wnumber(i)*1.2399e-4
1      continue
 
      return
      end
 
      subroutine fvalueli(z)
 
c____ This subroutine calculate absorption oscilator
c     strength. The oscillator strength are obtain form a paper by
c     W.L.Wiese,J.Phys.Chem.Ref.Data,Vol.5,No.3,1976 and from the
c     Klapish codes that are avaliable at LLNL.
 
 
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'flystuf'

      common /cli1comR/ sd1(6,2),sd2(6,2),se1(6,2),se2(6,2),
     1                  a(7),a1(7),b(7),c(7),d(7)
      common /cli1comI/ mm(19),ll(14)
c
      data mm /2*2,3*3,4*4,5*5,6,7,8,9,10/
      data ll/0,1,0,1,2,0,1,2,3,0,1,2,3,4/
c
      nbn = 19
      ndt = 14
      nlv = 25
c
      do 8 i=1,25
        do 9 j=1,25
            fli(i,j)=0.0
9       continue
8     continue
c
      z = atomz
      z1 = 1./z

c
c  calculate the oscillator strengths between the detailed levels
c
c---2p 2s
      fli(1,2)=3.*(0.0014 + 0.3873*z1 +1.0766*z1**2)
c
c---3p 2s
      fli(1,4)=abs(3.*(0.1459-0.4421*z1-0.2737*z1**2))
c
c---4p 2s
      fli(1,7)=abs(3.*(0.0367-0.102*z1-0.0432*z1**2))
c
c---5p 2s
      fli(1,11)=abs(3.*(0.0135+0.0038*z1-0.1603*z1**2))
c
c---3s 2p
      fli(2,3)=0.0157+0.0282*z1+0.5625*z1**2
c
c---3d 2p
      fli(2,5)=abs(0.6581+0.635*z1-8.1674*z1**2+23.4758*z1**3)
c
c---4s 2p
      fli(2,6)=0.0028+0.0161*z1+0.0507*z1**2
c
c---4d 2p
      fli(2,8)=abs(0.1243-0.0158*z1)
c
c---5s 2p
      fli(2,10)=0.0012+0.0107*z1
c
c---5d 2p
      fli(2,12)=abs(0.0557-0.2829*z1+2.4187*z1**2-6.9541*z1**3)
c
c---3p 3s
      fli(3,4)=abs(3.*(-0.0099+0.9866*z1))
c
c---4p 3s
      fli(3,7)=abs(3.*(0.1643-0.5163*z1))
c
c---5p 3s
      fli(3,11)=abs(3.*(0.0392-0.08*z1))
c
c---3d 3p
      fli(4,5)=abs(10.*(-2.924e-4+0.0454*z1-0.0339*z1**2))
c
c---4s 3p
      fli(4,6)=0.0364+0.0712*z1+1.111*z1**2
c
c---4d 3p
      fli(4,8)=abs(0.6002-0.3755*z1)
c
c---5s 3p
      fli(4,10)=abs(0.0129-0.1063*z1+1.1571*z1**2-2.924*z1**3)
c
c---5d 3p
      fli(4,12)=abs(0.1431-0.0947*z1)
c
c---4p 3d
      fli(5,7)=6.*(0.0019+0.0047*z1+0.0046*z1**2)
c
c---4f 3d
      fli(5,9)=1.0114+0.0357*z1
c
c---5p 3d
      fli(5,11)=abs(6.*(5.492e-4-0.004*z1+0.0448*z1**2-0.1311*z1**3))
c
c---5f 3d
      fli(5,13)=abs(0.1574-0.0049*z1)
c
c---4p 4s
      fli(6,7)=3.*(0.0038+0.898*z1+2.4571*z1**2)
c
c---5p 4s
      fli(6,11)=abs(3.*(0.1789-0.5074*z1))
c
c---4d 4p
      fli(7,8)=abs(0.0703-1.1186*z1+15.8497*z1**2-45.7841*z1**3)
c
c---5s 4p
      fli(7,10)=0.0624+0.0599*z1+1.8323*z1**2
c
c---5d 4p
      fli(7,12)=abs(0.5931-0.5178*z1)
c
c---4f 4d
      fli(8,9)=abs(0.0026-1.474e-4*z+8.750e-6*z**2)
c
c---5p 4d
      fli(8,11)=6.*(0.0048+0.0114*z1)
c
c---5f 4d
      fli(8,13)=0.8828+0.0336*z1
c
c---5d 4f
      fli(9,12)=abs(15.*(5.957e-4+4.190e-5*z1))
c
c---5g 4f
      fli(9,14)=1.3421+0.0285*z1
c
c---5p 5s
      fli(10,11)=abs(3.*(-0.0176+1.732*z1))
c
c---5d 5p
      fli(11,12)=abs(0.0949-1.4671*z1+21.1873*z1**2-61.0494*z1**3)
c
c---5f 5d
      fli(12,13)=abs(0.0053-2.925e-4*z+1.625e-5*z**2)
c
c---5g 5f
      fli(13,14)=abs(4.889e-5-1.461e-5*z+3.432e-6*z**2)
c
c  calculate oscillator strengths between the detailed levels
c  to the non-detailed levels and between the non-detailed levels
c
      if (ndt .eq. 0) then
         n=1
      else if (ndt .eq. 2) then
         n=3
      else if (ndt .eq. 5) then
         n=4
      else if (ndt .eq. 9) then
         n=5
      else if (ndt .eq. 14) then
         n=6
      endif
 
c
c   detailed to non-detailed
c
      do 70 i=1,ndt
         do 68 j=1,nbn-ndt
           ii = mm(i)
           zii2 = ii*ii
           jj = n+j-1
c hyun's correction
c           fli(i,ndt+j) = fos(ii,jj)*gli(i)/2./zii2           
           fli(i,ndt+j) = fos(ii,jj)
68        continue
70    continue
 
c    non-detailed to non-detailed
 
      do 74 i=1,nbn-ndt
         do 72 j=1,nbn-ndt
            if (i .ge. j) go to 72
            ii = n+i-1
            jj = n+j-1
            fli(i+ndt,j+ndt) = fos(ii,jj)
72       continue
74    continue
c
c zero out the unused f-values
c
      do 11 i=1,nbn
        do 12 j=nbn+1,nlv
          fli(i,j) = 0.0
12      continue
11    continue
      return
      end
 
 
      subroutine mewe1
 
c This program creates arrays that are need to calculate collisional
c rates using the paper Mewe. What is needed is the principal quantum
c numbers for each detailed state and the multiplet oscilator strength
c for thoes detailed state transitions that are not allowed. For example
c the transition 2s-ns and 2s-nd uses the f-value for 2s-np.
 
      implicit real*8 (a-h,o-z)
      include 'flystuf'
      common /mewe/ f1(14,14)
 
c   the multiplet f-values
       do 20 i=1,14
          do 18 j=1,14
              f1(i,j)=0.0
18        continue
20     continue
 
       f1(1,6) = fli(1,7)
       f1(1,8) = fli(1,7)
       f1(1,10) = fli(1,11)
       f1(1,12) = fli(1,11)
       f1(2,7) = fli(2,8)
       f1(2,9) = fli(2,8)
       f1(2,11) = fli(2,12)
       f1(2,13) = fli(2,12)
       f1(3,5) = fli(3,4)
       f1(3,6) = fli(3,7)
       f1(3,8) = fli(3,7)
       f1(3,10) = fli(3,11)
       f1(3,12) = fli(3,11)
       f1(4,7) = fli(4,8)
       f1(4,9) = fli(4,8)
       f1(4,11) = fli(4,12)
       f1(4,13) = fli(4,12)
       f1(5,6) = fli(5,7)
       f1(5,8) = fli(5,7)
       f1(5,10) = fli(5,11)
       f1(5,12) = fli(5,11)
       f1(5,14) = fli(5,13)
       f1(6,8) = fli(6,7)
       f1(6,11) = fli(6,10)
       f1(6,12) = fli(6,10)
       f1(7,9) = fli(7,8)
       f1(7,11) = fli(7,10)
       f1(7,13) = fli(7,12)
       f1(8,10) = fli(8,11)
       f1(8,12) = fli(8,11)
       f1(8,14) = fli(8,13)
       f1(9,11) = fli(9,12)
       f1(9,13) = fli(9,12)
       f1(10,12) = fli(10,11)
       f1(11,13) = fli(11,12)
       f1(12,14) = fli(12,13)
       return
       end

      subroutine sthul(iso,epion)
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'popstuf'
      include 'flystuf'
!     common /gamrts/gamrat(2*nlvi,0:miso)

      dimension elv(70),glv(70),mcont(70),match(70),mnner(70)
      dimension tfit(5),ltran(70,70),lnner(70,2),gamrt(70)
      dimension epion(0:miso)
      character*8 knam(70)

      dimension nlvhul(3)
      data nlvhul/22,61,51/ 
c

      do ilo=1,70
         do iup=1,70
            ltran(ilo,iup)=0
            idcxhl(ilo,iup,iso)=0
         enddo
      enddo
      read(31,101)isoin,nlast,nbn,nau
      do ii=1, nlast
         read(31,101) isoin,i,mcont(i),match(i),mnner(i),
     &        knam(i),glv(i),elv(i)
      enddo
      do ii=2, nlast
         read(31,102)(fhul(i,ii,iso),i=1,ii-1)
         do i=1, ii-1
            if(fhul(i,ii,iso).gt.0.d0)then
               ntran=ntran+1
               ltran(i,ii)=ntran
               do ix=1,5
                  chul(ntran,ix)=0.d0
               enddo
            endif
         enddo
      enddo

      read(31,105)(ahul(i,iso),i=1,nlast)
 99   read(31,103,end=100)isoin,ilo,iup,(tfit(ix),ix=1,5)
      if(isoin.eq.0)goto 100
      if(ltran(ilo,iup).eq.0)then
         ntran=ntran+1
         ltran(ilo,iup)=ntran
      endif
      itran=ltran(ilo,iup)
      if(iso.ne.isoin) then
         write(ioo,'(" *** check hullac data files")')
         write(16,'(" Warning:Check hullac data files")')
      endif
      idcxhl(ilo,iup,isoin)=itran
      do ix=1,5
         chul(itran,ix)=tfit(ix)
      enddo
      goto 99
 100  continue

c done with reading
      nlast=nbn+nau
      nlev(iso)=nlast
      nlbn(iso)=nbn

c calculate a-rates and gamma rates
      do j=1,nlev(iso)
         gamrt(j)=ahul(j,iso)
         do i=1,j-1
            emn=elv(j)-elv(i)
            oratehl(j,i,iso)=4.34327d7*emn**2*fhul(i,j,iso)*
     &           glv(i)/glv(j)
            gamrt(j)=gamrt(j)+oratehl(j,i,iso)
         enddo
      enddo

      do i=1, nlast
         do k=1, 2
            lnner(i,k)=0
         enddo
         ntot=ntot+1
         efromgr(ntot)=epion(iso)+elv(i)
         elev(ntot)=elv(i)
         glev(ntot)=glv(i)
         qlev(ntot)=atomz-iso+1.
         isolev(ntot)=iso
         ionstg(ntot)=iatomz-iso
         levl(ntot)=i
         indlev(i,iso)=ntot
         gamrat(i,iso)=gamrt(i)
         kname(ntot)=knam(i)
         lmatch(ntot)=match(i)
         lvcont(ntot)=indlev(mcont(i),iso-1)
         if(mcont(i).ne.mnner(i))lnner(i,1)=mnner(i)
         evtoip(ntot)=efromgr(lvcont(ntot))-efromgr(ntot)
         icont=lvcont(ntot)-indlev(1,iso-1)+1
         mobtot(ntot)=morb(match(i),iso)
         do ip=1,nprim
            noctot(ntot,ip)=nocel(match(i),iso,ip)
         enddo
         noctot(ntot,nprim+1)=0
      enddo
      do il=1, nlast
         i=indlev(il,iso)
         if(iso.ge.1) then
            if(lnner(il,1).ne.0) then
               ii=indlev(lnner(il,1),iso-1)
               knam(1)(1:6)=kname(ii)(1:6)
               kk=1
               do j=1, ii-1
                  jj=indlev(j,iso-1)
                  knam(2)(1:6)=kname(jj)(1:6)
                  if(knam(1)(1:6).eq.knam(2)(1:6)) then
                     kk=kk+1
                     lnner(il,kk)=j
                  endif
               enddo
            endif
         endif
      enddo

c  set up excitation transitions
      do ii=2,nlev(iso)
         do i=1,ii-1
            if(ltran(i,ii).ne.0)then
               lupper(ltran(i,ii))=indlev(ii,iso)
               llower(ltran(i,ii))=indlev(i,iso)
               ltype(ltran(i,ii))=1
               write(*,*) 'line 3618 in zchkfly.f *********' ! mscho 201025 - no screen 
               write(*,*) lupper(ltran(i,ii)), indlev(ii,iso) ! mscho 201020 - no screen 
               ! mscho : it means that sthul is not loaded in the main code (201025) 
            endif
         enddo
      enddo

c  set up ionization transitions

      do i=1, nlev(iso)
         ii=indlev(i,iso)
         j=lvcont(ii)-indlev(1,iso-1)+1
         ntran=ntran+1
         lupper(ntran)=indlev(j,iso-1)
         llower(ntran)=ii
         write(*,*) '3633 line in zchkfly.f **** mscho' 
         write(*,*) llower(ntran) !mscho 
         ltype(ntran)=2
         petrn(i,j,iso)=1.
         if(i.eq.1.and.iso.eq.2)petrn(i,j,iso)=2.
         istrn(i,j,iso)=mobtot(ii)
         do jj=1, 2             !K-shell ionization 
            j=lnner(i,jj)     ! upper level index
            if(j.ne.0) then
               petrn(i,j,iso)=1.
               if(i.eq.1.and.iso.eq.2)petrn(i,j,iso)=2.
               istrn(i,j,iso)=1
               ntran=ntran+1
               lupper(ntran)=indlev(j,iso-1)
               llower(ntran)=ii
               write(*,*) '3646 line in zchkfly.f ******* mscho' 
               write(*,*) llower(ntran) ! mscho: It also is not working in the screen  
               ltype(ntran)=2
            endif
         enddo
      enddo
      do i=1, nlev(iso)
         if(ahul(i,iso).gt.0.d0)then
            ntran=ntran+1
            lupper(ntran)=indlev(1,iso-1)
            llower(ntran)=indlev(i,iso)
            write(*,*) ' wowowowowowowow mscho test 201025' ! mscho not working (no screen)  
            ltype(ntran)=3
         endif
      enddo

 101  format(5i5,1x,a8,1x,f10.2,2x,f15.6)
 102  format(1p,10e9.2)
 103  format(i3,2i6,1p,8e13.5)
 105  format(1p,10e9.2)
 1000 format (12i10)
 1010 format (1p,10e12.5)

      return
      end
c---------------------------------------------------
      function chl(iir,ijr,iz)
c
c___ this calculates the collisional excitation or de-excitation
c___ from level i to level j of an ion (with iz electrons)
c___ in the average-atom model.
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'popstuf'
      include 'flystuf'
      dimension tfit(5)

      i=min(iir,ijr)
      j=max(iir,ijr)
      ii=indlev(i,iz)
      jj=indlev(j,iz)
      gl=glev(ii)
      gu=glev(jj)

      chl=0.d0
      de=elev(jj)-elev(ii)
      if(de.le.0.d0) return
      ebt = de/tev
      fij = fhul(i,j,iz)
      if(fij.eq.0.d0) goto 100
      gaunt = 0.15 + 0.28*(log(1.+1./ebt)-0.4/(ebt+1.)**2)
      if(iir.eq.j) then
         chl=1.578d-5/de*fij*gaunt/sqrt(tev)*gl/gu  ! de-excitation
      else
         chl=1.578d-5/de*fij*gaunt/sqrt(tev)*exp(-ebt) ! excitation
      endif
 100  continue

      if(j.le.nlbn(iz))then
         ntrn=idcxhl(i,j,iz)
         if(ntrn.eq.0) goto 200
         do ix=1,5
            tfit(ix)=chul(ntrn,ix)
         enddo
         if(chl.gt.1.d-30.and.ebt.le.200.0)then
            chl=eval2(tev,tfit,gl,de)
            if(iir.eq.j) then
               chl=chl*exp(ebt)*gl/gu ! For de-excitation
            endif
         endif
      endif
 200  return
      end

c *********************************************************************
      function eval2(tev,a,gl,de)
      implicit real*8 (a-h,o-z)
      dimension a(5)
      data  tol/1.d-8/, mm/1/, kode/2/, nn/1/
      external exp1x
c
      c0=a(1)
      c1=a(2)
      c2=a(3)
      c3=a(4)
      a0=a(5)
      a1=a0+1.d0
      tx = de/tev

      eval2 = 8.01d-8/gl/tev**0.5*exp(-tx)
     &      *(c0*exp1x(tx) + (c1 + c3*tx/a1)+
     &     (c2-c3*tx)*tx*exp1x(tx*a1))
      if(eval2 .lt.0.d0)  eval2=0.
      return
      end
c
c *********************************************************************

c---------------------------------------------------
      subroutine Vchlfill(iz)
c
c___ this calculates the collisional excitation or de-excitation
c___ from level i to level j.
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'popstuf'
      include 'xc_stuf'
      include 'flystuf'
      dimension tfit(5)

      data conxx /7.099d-8/!concx*5.931d7
c      data concx /1.197d-15/
c___ Fill the arrays of excitation and de-excitation cross-sections
c___ versus energy.  Noting that the 1/e dependence of the xc
c___ implied that, e.g., the collision strengths, e*xc(e) is energy
c___ independent and thus the excitation = de-excitation
c
      
      do i=1,nlev(iz)
         do j=1,nlev(iz)
            Vchl(i,j,iz) = 0.0
         enddo
      enddo
      if(.not.cxflag) return

      naatotc=nlev(iz) ! mscho 200513 naatot -> naatotc (for type) 
      do 1000 i=1, naatotc-1
         ii=indlev(i,iz)
         gl=glev(ii)
         do 900 j=i+1,naatotc
            jj=indlev(j,iz)
            gu=glev(jj)
            eji = elev(jj)-elev(ii)
            if(eji.le.0.d0) goto 900
            fij=fhul(i,j,iz)
            if(fij.le.0.d0) goto 900
            ieji = indxc(eji)
            gij = gl/gu
            
            sumex=0.d0
            sumdx=0.d0
            if(j.le.nlbn(iz))then
               ntrn=idcxhl(i,j,iz)
               if(ntrn.ne.0) then
                  do ix=1,5
                     tfit(ix)=chul(ntrn,ix)
                  enddo
                  do k=1,nebins
                     e = energy(k)
                     de = denergy(k)
                     if (k-ieji .ge. 1 .and. e .gt. eji) then
                        xx=e/de
                        ocx=tfit(1)*log(xx)+tfit(2)+tfit(3)/(xx+tfit(5)) 
     &                       + tfit(4)/(xx+tfit(5))**2
c     crx=concx/gl/e*ocx*5.931E7*sqrt(e)*sqrt(e)*f(e)
                        crxx=conxx/gl*ocx
                        sumex = sumex + fe(k)*denergy(k)*crxx
                        sumdx = sumdx + fe(k-ieji)*denergy(k-ieji)
     &                       *gij*crxx
                     endif
                  enddo 
               else
                  do k=1,nebins
                     e = energy(k)
                     de = denergy(k)
                     if (k-ieji .ge. 1 .and. e .gt. eji) then
                        g=0.15+0.28*log(e/eji)
                        sumex=sumex+fe(k)*denergy(k)*1.40d-5*fij/eji*g
                        sumdx=sumdx+fe(k-ieji)*denergy(k-ieji)
     &                       *gij*1.40d-5*fij/eji*g
                     endif
                  enddo 
               endif
            else
               do k=1,nebins
                  e = energy(k)
                  if (k-ieji .ge. 1 .and. e .gt. eji) then
                     g=0.15+0.28*log(e/eji)
                     sumex=sumex+fe(k)*denergy(k)*1.40d-5*fij/eji*g
                     sumdx=sumdx+fe(k-ieji)*denergy(k-ieji)
     &                    *gij*1.40d-5*fij/eji*g
                  endif
               enddo 
            endif
 800        Vchl(i,j,iz)=sumex
            Vchl(j,i,iz)=sumdx
 900     continue
 1000 continue

      return
      end
c
c***********************************************************
      function capehl(i,iz)
c
c---- rate of electron capture from the state (j, iz-1)
c----     into autionizing state (i,iz)
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'popstuf'
      include 'flystuf'

      ntoti=indlev(i,iz)
      ntotj=indlev(1,iz)
      eii =elev(ntoti)-eipz(iz)+debye(iatomz-iz+1)
      capehl = 1.65615d-22*exp(-eii/tev)/tev**1.5*
     &     glev(ntoti)/glev(ntotj)*ahul(i,iz)
      return
      end

c
c***********************************************************
c
      function Vcapehl(i,iz)
c
c---- rate of electron capture from helium-like ground state
c----     into autionizing state i
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'xc_stuf'
      include 'popstuf'
      include 'flystuf'

c
      Vcapehl=0.d0
      if(.not.cxflag) return


      ntoti=indlev(i,iz)
      ntotj=indlev(1,iz)
      eii =elev(ntoti)-eipz(iz)+debye(iatomz-iz+1)
      capex = glev(ntoti)/glev(ntotj)*ahul(i,iz)
      k=indxc(eii)
      if(k.gt.0.and.k.le.nebins) then
         Vcapehl=1.40e-22*capex*energy(k)*fe(k)
      endif
      return
      end

