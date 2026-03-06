      subroutine runzelda(it,izf)
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'xc_stuf'
      include 'timestuf'
      include 'popstuf' 

      do kbin=1, nebins
         zqs(kbin)=0.d0
         zqr(kbin)=0.d0
      enddo

      call getzqs
      
      !write(*,*) 'mscho test 201021: line 14 in zchkxc.f'
      !write(*,*) llower(20),lupper(20)
      
      call Zelda(izf)

!      write(ioo,'(" *****: t, Te, Ne",1p,e12.5,0p,f13.5,1x,1p,e12.5)')
!    &     times(it),tevfout(it),fevfout(it)

      return
      end

      subroutine getzqs
      implicit double precision(a-h,o-z)

      include 'mainstuf'
      include 'xc_stuf'
      include 'popstuf'
      include 'runstuf'

      double precision N_upper, N_lower
      !write(*,*) 'line 34 in zchkxc.f' 
      !write(*,*) llower(20), lupper(20)
      xsum=0.d0
      ysum=0.d0
      do 100 itran=1, ntran
         !write(*,*) '40 in zchkxc.f:ltype(itran)',ltype(itran)
         if(ltype(itran).le.1)goto 100
         llo=llower(itran)
         lup=lupper(itran)
         ind1=indl2p(llo)
         ind2=indl2p(lup)         
         !write(*,*) llower(itran), lupper(itran) ! mscho 201022 - it is working well

         if(ind1*ind2.eq.0) goto 100
         N_upper=popout(lup)
         if(N_upper.lt.1.e-30)goto 100
         N_lower=popout(llo)
         if(N_lower.lt.1.e-30) goto 100

         ilo=levl(llo)
         iup=levl(lup)
         iso=isolev(llo)
      
         if(ltype(itran).eq.2)then
            do j=1, nebins
               zqr(j)=zqr(j)+N_upper*xc_alfxx(j,ilo,iup,iso)
               xsum=xsum-N_upper*xc_alfxx(j,ilo,iup,iso)*fe(j)
            enddo
            ysum=ysum+N_upper*Valfxx(ilo,iup,iso)*denow

            if(radflag.ne.'off') then
               do j=1, nebins
                  zqs(j)=zqs(j)+N_lower*xc_bcontxx(j,ilo,iup,iso)
                  xsum=xsum+N_lower*xc_bcontxx(j,ilo,iup,iso)
                  zqr(j)=zqr(j)+N_upper*xc_bstmxx(j,ilo,iup,iso)
                  xsum=xsum-N_upper*xc_bstmxx(j,ilo,iup,iso)*fe(j)

               enddo
            endif
            ysum=ysum-N_upper*Vbstmxx(ilo,iup,iso)*denow
            ysum=ysum+N_lower*Vbcontxx(ilo,iup,iso)

         else if(ltype(itran).eq.3)then
            eip=elev(llo)-elev(lup)-eipz(iso)+debye(iatomz-iso+1)
            j=indxc(eip)
            if(j.ge.1.and.j.le.nebins)then
               glbygu=glev(llo)/glev(lup)
               zqs(j)=zqs(j)+autoxx(ilo,iup,iso)*N_lower
               xsum=xsum+N_lower*autoxx(ilo,iup,iso)

               xxcape=1.4e-22*autoxx(ilo,iup,iso)*energy(j)*glbygu
               zqr(j)=zqr(j)+xxcape*N_upper
               xsum=xsum-N_upper*xxcape*fe(j)*denow
            endif
            ysum=ysum+N_lower*autoxx(ilo,iup,iso)
            ysum=ysum-N_upper*Vcapexx(ilo,iup,iso)
         endif
 100  continue

c      print *, 'extra source- ion', xsum,ysum

c      do ir=1, nlev(iz)
c         iir=indlev(ir,iz)
c         if(indl2p(iir).ne.0)then
c         do jr=1, nlev(iz-1)
c            ijr=indlev(jr,iz-1)
c            if(indl2p(ijr).ne.0.and.petrn(ir,jr,iz).gt.0.)then
c               do j=1, nebins
c                  zqr(j)=zqr(j)+popout(ijr)*xc_alfxx(j,ir,jr,iz)
c                  xsum=xsum+popout(ijr)*xc_alfxx(j,ir,jr,iz)*fe(j)
c               enddo
c               ysum=ysum+popout(ijr)*Valfxx(ir,jr,iz)
c               if(radflag.ne.'off') then
c                  do j=1, nebins
c                     zqs(j)=zqs(j)+popout(iir)*xc_bcontxx(j,ir,jr,iz)
c                     zqr(j)=zqr(j)+popout(ijr)*xc_bstmxx(j,ir,jr,iz)
c                  enddo
c               endif
c            endif
c         enddo
c         endif
c      enddo

c      print *, 'extra source- ion', iz,xsum*denow,ysum*denow

c autoionization and electron capture

c      do ir=1, nlev(iz)
c         iir=indlev(ir,iz)
c         if(indl2p(iir).ne.0)then
c         do jr=1, nlev(iz-1)
c            ijr=indlev(jr,iz-1)
c            if(indl2p(ijr).ne.0.and.autoxx(ir,jr,iz).gt.0.d0)then
c               eip=elev(iir)-elev(ijr)-eipz(iz)+debye(iatomz-iz+1)
c               j=indxc(eip)
c               if(j.ge.1)then
c               glbygu=glev(iir)/glev(ijr)
c               xxcape=1.4e-22*autoxx(ir,jr,iz)*energy(j)*glbygu
c               zqs(j)=zqs(j)+autoxx(ir,jr,iz)*popout(iir)
c               zqr(j)=zqr(j)+xxcape*popout(ijr)
c               endif
c            endif
c         enddo
c         endif
c      enddo

      return
      end
      

ccccccccccccccccccccccccccccc
c         
c        *******   ******   *        *****       ***   
c             *    *        *        *    *     *   *   
c            *     *        *        *     *   *     * 
c           *      ***      *        *     *   ******* 
c          *       *        *        *     *   *     * 
c         *        *        *        *    *    *     * 
c        *******   ******   * * **   *****     *     * 
c
      double precision function xc_alfxx(j,i,k,iz)
      implicit real*8 (a-h,o-z)
c
c___ this calculates the radiative recombination from ion stage
c___ i to i-1. this is for stages from hydrogen like to singly ionized
c___ Note that the charge on the ion i is i-1; i.e., i=1 is neutral
c
      include 'mainstuf'
      include 'popstuf'
      include 'xc_stuf'

      xc_alfxx=0.d0
      if(k.ne.1)return
      ntoti=indlev(i,iz)
      ntotj=indlev(k,iz-1)
      eii=eipz(iz)+elev(ntotj)-elev(ntoti)-debye(iatomz-iz+1)
      if(eii.le.exthr) return

      qq=qlev(ntoti)
      pn=petrn(i,k,iz)
      glbygu = glev(ntoti)/glev(ntotj)

      e = energy(j)
      de = denergy(j)
      xc_alfxx=1.692e-15*pn*glbygu*eii**2.5/qq/(eii+e)*de
      return
      end
c
c
       double precision function xc_gbetaxx(j,k2,i,k,iz)
c
c___ this calculates the three body recombination from ion stage
c___ i to i-1 for hydrogenic to singly ionized.
c___ Note that the charge on the ion i is i-1; i.e., i=1 is neutral
c___ the limits on the integral over xc_betali are from eip to inifinity
c

      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'popstuf'
      include 'xc_stuf'

c     gain by 3-body recombination at j from k2
c     energy(j)=energy(k2)+energy(new-e-)+eii

      xc_gbetaxx=0.d0
      ntoti=indlev(i,iz)
      ntotj=indlev(k,iz-1)
      eii=eipz(iz)+elev(ntotj)-elev(ntoti)-debye(iatomz-iz+1)
      if(eii.le.exthr) return

      ieii=indxc(eii)
      pn=petrn(i,k,iz)
      qq=qlev(ntoti)-1
      bfact=0.25*((100.*qq+91.)/(4.*qq+3.))**0.5 - 1.25
      ci=1.4679d-22*glev(ntoti)/glev(ntotj)

      e=energy(j)
      if(e.le.eii)return

      etot=e-eii
      uu=e/eii

      wfact=(log(e/eii))**(bfact*eii/e)
      xctot=2.218d-6*pn/eii*log(e/eii)*wfact      
c lotz     xctot=2.669e-6*pn/eii*log(e/eii)

      e2=energy(k2)/etot
      indx=indxc(e-energy(k2)-eii) ! new e-
      if(indx.ge.1.and.e2.le.1.d0)then
         qj=qdiff(etot,uu,e2)/2.d0
         xc_gbetaxx=ci*xctot*qj*fe(indx)*denergy(k2)*denergy(j)
      endif
      return
      end
c

c
       double precision function xc_betaxx(k2,i,k,iz)
c
c___ this calculates the three body recombination from ion stage
c___ i to i-1 for hydrogenic to singly ionized.
c___ Note that the charge on the ion i is i-1; i.e., i=1 is neutral
c___ the limits on the integral over xc_betali are from eip to inifinity
c

      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'popstuf'
      include 'xc_stuf'

c     loss by 3-body recombination at energy(j)
c     energy(j)=energy(k2)+energy(indx)+eii

      xc_betaxx=0.d0
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
      ie1=indxc(energy(1)+energy(k2)+eii)
      do j=1, nebins
         e=energy(j)
         indx=indxc(e-energy(k2)-eii)
         if(e.gt.energy(indx)+eii.and.indx.ge.1)then
            uu=e/eii
            etot=e-eii
            e2=energy(indx)/etot
c            if(e2.lt.1.d0)then
               wfact=(log(e/eii))**(bfact*eii/e)
               xctot=2.218d-6*pn/eii*log(e/eii)*wfact      
c lotz         xctot=2.669e-6*pn/eii*log(e/eii)
               qj=qdiff(etot,uu,e2)
               sum=sum+xctot*qj*fe(indx)*denergy(indx)*denergy(k2)
c            endif
         endif
      enddo
      xc_betaxx=ci*sum

      return
      end
c
c
c
c
      double precision function xc_bstmxx(j,i,k,iz)
c
c___ calculate the stimulated recombination from ion stage i to i-1
c___ Note that the charge on the ion i is i-1; i.e., i=1 is neutral
c___ A photon of energy e+IP takes an electron from bin e 
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'runstuf'
      include 'xc_stuf'
      include 'popstuf'

      xc_bstmxx =0.d0
      ntoti=indlev(i,iz)
      ntotj=indlev(k,iz-1)
      eii=eipz(iz)+elev(ntotj)-elev(ntoti)-debye(iatomz-iz+1)
      if(eii.le.exthr)return

      ebt = eii/tev
      glbygu = glev(ntoti)/glev(ntotj)
      qq= qlev(ntoti)
      pn = petrn(i,k,iz)

      e = energy(j)
      de = denergy(j)
      if (radflag .ne. 'trfile') then
         field_e = bbfield((e+eii)/trev)
      else
         field_e = field(e+eii)
      endif
      xc_bstmxx=8.167e-12*glbygu*eii**2.5/qq
     1        /(e+eii)**4*field_e*de*pn
      return
      end
c
c
c
c
      double precision function xc_caa(k,iir,ijr,iz)
c
c___ this calculates the collisional excitation or de-excitation
c___ from level i to level j.
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'xc_stuf'
c
c___ Fill the arrays of excitation and de-excitation cross-sections
c___ versus energy.  Noting that the 1/e dependence of the xc
c___ implied that, e.g., the collision strengths, e*xc(e) is energy
c___ independent and thus the excitation = de-excitation
c
      
      xc_caa= 0.d0
      i=min(iir,ijr)
      j=max(ijr,iir)

      eji = evaa(j,iz)-evaa(i,iz) ! energy difference
      if(eji.le.exthr)return
      ieji = indxc(eji)
      fij=faa(i,j,iz)
      if(fij.le.0.d0)return
      zeff = sqq(i,iz)
      gij = gaa(i,iz)/gaa(j,iz)

      de = denergy(k)
      e = energy(k)
      g = 0.15+0.28*log(e/eji)
      if(i.eq.iir)then
         xc_caa=1.40d-5*fij/eji*g*de
      else
         xc_caa=gij*1.40d-5*fij/eji*g*de
      endif
      
      if(npqe(i,j,iz).ne.0.and.isamp(1).eq.2) then
         n1=nloex(i,j,iz)       ! initial shell of lower level
         n2=nhiex(i,j,iz)       ! initial shell of lower level
         pm=npqe(i,j,iz)
         if(n1.le.3.and.(n2-n1).le.2) then
            ii=2*(n1-1)+(n2-n1)
            if(i.eq.iir)then
               xc_caa=pm*7.096e-8/(2*n1**2)*cstrsamp(e,ii)/zeff**2*de
            else
               xc_caa=gij*
     &              pm*7.096e-8/(2*n1**2)*cstrsamp(e,ii)/zeff**2*de
            endif
         endif
      endif

      return
      end
c
c
      double precision function xc_sxx(j,i,k,iz)
c
c___ this calculates the collisional ionization from state i
c___ to i+1 for all stages from neutral to helium like.
c___ Note that the charge on the ion i is i-1; i.e., i=1 is neutral
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'popstuf'
      include 'xc_stuf'

      data y1/3.5d5/,ems/5.11d5/,ems2/2.61121d5/

c
      xc_sxx=0.d0
      ntoti=indlev(i,iz)
      ntotj=indlev(k,iz-1)
      eii=eipz(iz)+elev(ntotj)-elev(ntoti)-debye(iatomz-iz+1)
      if(eii.le.exthr)return
      ieii = indxc(eii)

      pn=petrn(i,k,iz)
      qq=qlev(ntoti)-1
      bfact=0.25*((100.*qq+91.)/(4.*qq+3.))**0.5 - 1.25

   
      e = energy(j)
      if(e .gt. eii) then
         wfact=(log(e/eii))**(bfact*eii/e)
         uu=e/eii
         etot=e-eii
         qsum=0.d0
         do k2=1,j
            e2=energy(k2)/etot  ! secondary e- divided by total energy
            indx=indxc(e-energy(k2)-eii)
            if(e2.lt.1.and.indx.ge.1)then
               qj=qdiff(etot,uu,e2)/2.d0 
               qsum=qsum+qj*denergy(k2)
            endif
         enddo
         xc_sxx=qsum*2.218d-6*pn/eii*log(e/eii)*wfact*denergy(j)
      endif
      return
      end
c

      double precision function xc_gsxx(k2,j,i,k,iz)
c
c___ this calculates the collisional ionization from state i
c___ to i+1 for all stages from neutral to helium like.
c___ Note that the charge on the ion i is i-1; i.e., i=1 is neutral
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'popstuf'
      include 'xc_stuf'

      data y1/3.5d5/,ems/5.11d5/,ems2/2.61121d5/

c     e(j)=e(k2)+e(new e-)+eii
c     j=  initial bin of incoming e-
c     k2= final bin of outgoing e-

      xc_gsxx=0.d0
      ntoti=indlev(i,iz)
      ntotj=indlev(k,iz-1)
      eii=eipz(iz)+elev(ntotj)-elev(ntoti)-debye(iatomz-iz+1)
      if(eii.le.exthr)return
      if(energy(j).le.energy(k2)+eii)return

      ieii=indxc(eii)
      pn=petrn(i,k,iz)
      qq=qlev(ntoti)-1
      bfact=0.25*((100.*qq+91.)/(4.*qq+3.))**0.5 - 1.25

      e = energy(j)
      etot=e-eii
      e2=energy(k2)/etot
      indx=indxc(e-energy(k2)-eii)
      if(e2.lt.1.and.indx.ge.1)then
         uu=e/eii
         wfact=(log(e/eii))**(bfact*eii/e)
         xctot=2.218d-6*pn/eii*log(e/eii)*wfact
         qj=qdiff(etot,uu,e2)
         xc_gsxx=xctot*qj*denergy(k2)*denergy(j)
      endif
      return
      end
c

      
c
c     photoionization rates
c
      double precision function xc_scontaa(j,i,k,iz)
c
c     bbfiled is in the unit of erg/cm2/sec/Hz
c     use J. Scofield's photoionization data for inner-shell
c     K-shell photoionization cx are almost uniform for all ions
c     However continuum edge changes accordingly therefore
c     use the relevant continuum edge for each ion
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'runstuf'
      include 'timestuf'
!     include 'kinstuf'
      include 'xc_stuf'

      common /phistuff/ cp0,cp1,cp2,cp3,xn,cep,cem
      data fourpi /12.5664d0/

c
c___ calculate K-shell photoionization from hydrogenic state i
c       
c
      xc_scontaa=0.d0

c     use relevant continuum edge for integration
      eii = eipz(iz)+evaa(k,iz-1)-evaa(i,iz)
     &     -debye(iatomz-iz+1)
      if(eii.le.exthr)return


      ish=nish(i,k,iz)
      cp0=cphilsh(1,ish,iz)
      cp1=cphilsh(2,ish,iz)
      cp2=cphilsh(3,ish,iz)
      cp3=cphilsh(4,ish,iz)
      xn =cphilsh(5,ish,iz)
      pn =cphilsh(5,ish,iz)
      cep=cphilsh(6,ish,iz) 
      cem=cphilsh(7,ish,iz) 

      if(pn.eq.0.d0.or.eii.lt.0.9*cep)then
         xc_scontaa=xc_bcontxx(j,i,k,iz)
         return
      endif


      const = fourpi*1.d-18*(13.606/cep)*pn/hp

      ebt = cep  ! for a reproduction of Scofield's cx

      e=energy(j)
      de=denergy(j)*const
      xx=log(1.d0+(e/ebt))

      if (radflag .ne. 'trfile') then
         xc_scontaa=phicx(xx)/(e+eii)*bbfield((e+eii)/trev)*de
      else
         xc_scontaa=phicx(xx)/(e+eii)*field(e+eii)*de
      endif

      return
      end

      double precision function xc_bcontxx(j,i,k,iz)
c
c___ photoionization of ion species i using radiation field
c___ Note that the charge on the ion i is i-1; i.e., i=1 is neutral
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'runstuf'
      include 'popstuf'
      include 'xc_stuf'
c

      !open(unit=76,file='pe_monitor',status='old',position='append')
      xc_bcontxx=0.d0

      ntoti=indlev(i,iz)
      ntotj=indlev(k,iz-1)
      eii= elev(ntotj)+eipz(iz)-debye(iatomz-iz+1)-elev(ntoti)
      if(eii.le.exthr)return
      qq= qlev(ntoti)
      pn = petrn(i,k,iz)

      e = energy(j)
      de = denergy(j)

      if (radflag .ne. 'trfile') then
         field_e = bbfield((e+eii)/trev)
      else
         field_e = field(e+eii)
      endif
      xc_bcontxx= 5.5332e+10*eii**2.5/qq/(e+eii)**4
     1           *field_e*de*pn

      !if(iz == 10) then 
          !write(76,*)elev(ntotj),eipz(iz),debye(iatomz-iz+1),elev(ntoti)
          !write(76,*)eii, j, i, k, iz, field_e, xc_bcontxx

      !endif
      !close(76) 
      return
      end


c-----------------------------------------------------------
c
      function  xc_chy(k,ir,jr)
c
c___ this calculates the collisional excitation or de-excitation
c___ from level i to level j.
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'flystuf'
      include 'xc_stuf'
      
     
c___ Fill the arrays of excitation and de-excitation cross-sections
c___ versus energy.  Noting that the 1/e dependence of the xc
c___ implied that, e.g., the collision strengths, e*xc(e) is energy
c___ independent and thus the excitation = de-excitation
c
      xc_chy=0.d0

      i=min(ir,jr)
      j=max(ir,jr)
      
      

      fij = fos(i,j)
      eij = evhy(j)-evhy(i)
      if(eij.le.exthr)return

      ieji= indxc(eij)
      if(ieji.le.0.or.ieji.gt.nebins)return

      gij = ghy(i)/ghy(j)

      e = energy(k)
      de = denergy(k)
      if(i.eq.ir)then
         if (k-ieji.ge.1.and.e .ge. eij) then
            xc_chy=2.823e-6*fij/eij*de
         endif
      else
         if (k+ieji.le.nebins)then
            xc_chy=2.823e-6*fij/eij*de
         endif
      endif
c
      return
      end
c
c
      function xc_che(k,i,j)
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
 
      xc_che= 0.0
c
c___ calculate the collisional rates between helium like levels i to j
c
      kkl=min(i,j)
      kku=max(i,j)
      kkj = kkl
      kki = kku
      if (kki .gt. nbn) kki = kki-nbn-nbn0
      if (kkj .gt. nbn) kkj = kkj-nbn-nbn0
      ku=kki
      kl=kkj

      eu = evhe(kku)
      el = evhe(kkl)
      gu = ghe(kku)
      gl = ghe(kkl)
      delta = eu-el
      if(delta.le.exthr)return

      ieji=indxc(delta)
      if(ieji.le.0.or.ieji.gt.nebins)return

      e = energy(k)
      de = denergy(k)
      if(i.eq.kku.and.k+ieji.gt.nebins)return

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
 101  continue
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
 111  continue
      partde = 2.*(xcsamd(zbyd,3)+xcsamd(zbyd,4)+xcsamd(zbyd,5))
      if (ebyd .ge. 1.) then      
         partex = 2.*(xcsamd(ebyd,3)+xcsamd(ebyd,4)+xcsamd(ebyd,5))
      endif
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
      partde = .75*xcspnn1(5,zbyd)
      if (ebyd .ge. 1.) then
         partex = .75*xcspnn1(5,ebyd)
      endif
      go to 1000
 132  continue
      partde = xcspnno(1,5,zbyd)
      if (ebyd .ge. 1.) then
         partex = xcspnno(1,5,ebyd)
      endif
      go to 1000
 142  continue
      partde = .75*xcspnn1(1,zbyd)
      if (ebyd .ge. 1.) then
         partex = .75*xcspnn1(1,ebyd)
      endif
      go to 1000
 152  continue
      partde = 3.*(xcsamd(zbyd,6)+xcsamd(zbyd,7)+xcsamd(zbyd,8))
      if (ebyd .ge. 1.) then
         partex = 3.*(xcsamd(ebyd,6)+xcsamd(ebyd,7)+xcsamd(ebyd,8))
      endif
      go to 1000
c     
c___  2 singlet s state
c
 103  continue
      khl=2
      if (ku.gt.idsamp1) go to 4000
      if (ku.eq.idsamp1) go to 113
      if (ku-4) 123,123,133
 123  continue
      partde = .75*xcspnn1(1,zbyd)
      if (ebyd .ge. 1.) then
         partex = .75*xcspnn1(1,ebyd)
      endif
      go to 1000
 133  continue
      partde = xcspnno(1,4,zbyd)
      if (ebyd .ge. 1.) then
         partex = xcspnno(1,4,ebyd)
      endif
      go to 1000
 113  continue
      partde = (xcsamd(zbyd,6)+xcsamd(zbyd,7)+xcsamd(zbyd,8))
      if (ebyd .ge. 1.) then
         partex = (xcsamd(ebyd,6)+xcsamd(ebyd,7)+xcsamd(ebyd,8))
      endif
      go to 1000
c
c___ 2 triplet p state
c
 104  continue
      khl=2
      if (ku.gt.idsamp1) go to 4000
      if (ku.eq.idsamp1) go to 114
      partde = .75*xcspnn1(6,zbyd)
      if (ebyd .ge. 1.) then
         partex = .75*xcspnn1(6,ebyd)
      endif
      go to 1000
 114  continue
      partde = 3.*(xcsamd(zbyd,9)+xcsamd(zbyd,10)+xcsamd(zbyd,11))
      if (ebyd .ge. 1.) then
         partex = 3.*(xcsamd(ebyd,9)+xcsamd(ebyd,10)+
     &        xcsamd(ebyd,11))
      endif
      go to 1000
c     
c___ 2 singlet p state
c
 105  continue
      khl=2
      if (ku.gt.idsamp1) go to 4000
      if (ku.eq.idsamp1) go to 107
 107  continue
      partde = (xcsamd(zbyd,9)+xcsamd(zbyd,10)+xcsamd(zbyd,11))
      if (ebyd .ge. 1.) then
         partex = (xcsamd(ebyd,9)+xcsamd(ebyd,10)+xcsamd(ebyd,11))
      endif
      go to 1000
c
c---- use Vinogradovs rate for n = 1,2,3,4,5 - sov.j.quan.ele.,5,p630,(1978)
c
 3000 continue
      partde = xc_russian(kl,ku,zbyd)
      if (ebyd .ge. 1.) then
         partex = xc_russian(kl,ku,ebyd)
      endif
      go to 1000
c     
c___  levels treated as hydrogenic n .ge. 3
c     
 106  continue
      khl=kl-3
      go to 4000
c
c___ use sampsons data for collisions between auto-ionizing levels
c___            goett etal. atomic and nuc. data tables v25, p185 (1980)
c     
 109  continue
      kal = i-nbn
      kau = j-i
      partde = coa(kau,kal)+1.66*z2s2(kau,kal)*log(zbyd)+
     1     c1a(kau,kal)/(aaut(kau,kal)+zbyd)+
     2     c2a(kau,kal)/(aaut(kau,kal)+zbyd)**2
      if (ebyd .ge. 1.) then
         partex = coa(kau,kal)+1.66*z2s2(kau,kal)*log(ebyd)+
     1        c1a(kau,kal)/(aaut(kau,kal)+ebyd)+
     2        c2a(kau,kal)/(aaut(kau,kal)+ebyd)**2
      endif
      go to 1000
      
 4000 continue
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
 5000 continue
      
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
 1000 continue
           
      if (i.eq.kkl) then ! excitation 
         xc_che=7.067e-8*partex/ghe(kkl)/zeff**2*de
      else
         xc_che=7.067e-8*partde/ghe(kku)/zeff**2*de
      endif
      
 100  continue
      
 2000 continue

      return
      end
c
c
      function xc_cli(k,ir,jr)
 
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
      i=min(ir,jr)
      j=max(ir,jr)
      deli = abs(evli(j)-evli(i))
      if(deli.le.exthr)return

      ieji=indxc(deli)
      if(ieji.le.0.or.ieji.gt.nebins)return
      glbygu = gli(i)/gli(j)
      de = denergy(k)
      xc_cli=0.d0

      if(i.eq.j)return
c

      if((i.eq.1.or.i.eq.2).and.(j.ge.1.and.j.le.5))then
         if(i.eq.1)n=j-1
         if(i.eq.2)n=j+2
         if (deli .eq. 0.d0.or.e.le.deli) return
         ff = max(fli(i,j),1.d0)
         if(i.eq.ir)then        ! excitation
            e = energy(k)
            y = e/deli
            xc_cli=xc/deli*ff*de*
     1           (a(n)+a1(n)/(z-2)+b(n)/y+c(n)/y**2+d(n)*log(y))
            
         else
            e = energy(k+ieji)
            yp = e/deli
            xc_cli=xc/deli*ff*de*glbygu*
     1           (a(n)+a1(n)/(z-2)+b(n)/yp+c(n)/yp**2+d(n)*log(yp))
            
         endif
         return
      endif
c
c  These collisional excitation rates for 3s-3p,3s-3d,and 3p-3d
c  are calculated using the approximation in the paper
c  of Mewe Astron. & Astrophys.20,215-221(1972).
c
      

      if((i.eq.3.or.i.eq.4).and.(j.ge.4.and.j.le.5))then
         if (fli(i,j) .eq. 0.) then
            if(i.eq.ir)then     ! excitation
               xc_cli=xc*f1(i,j)/deli*.15*de
            else
               if(k+ieji.gt.nebins)return
               xc_cli=xc*f1(i,j)/deli*.15*de*glbygu
            endif
         else
            if(i.eq.ir)then
               e = energy(k)
               y = e/deli
               xc_cli=xc*fli(i,j)/deli
     1              *(.6+.28*log(y))*de
            else
               if(k+ieji.gt.nebins)return
               e = energy(k+ieji)
               yp = e/deli
               xc_cli=xc*fli(i,j)/deli
     1              *(.6+.28*log(yp))*de*glbygu
            endif
         endif
         return
      endif

c
c  using approximations from Mewe for the remaining transitions
c
c  collisional rates between detailed states
c
      if((i.ge.1.and.i.le.ndt).and.(j.ge.6.and.j.le.ndt))then
         if (i .ge. j)return
c
c  test if transition is allowed
c
         if ( fli(i,j).ne.0.) then
            if ( mm(i)-mm(j) .eq. 0.) then
               if(i.eq.ir)then
                  e = energy(k)
                  y = e/deli
                  xc_cli=xc*fli(i,j)/deli*(.6+.28*log(y))*de
               else
                  if(k+ieji.gt.nebins)return
                  e = energy(k+ieji)
                  yp= e/deli
                  xc_cli=xc*fli(i,j)/deli*(.6+.28*log(yp))*de*glbygu
               endif
            else
               if(i.eq.ir)then
                  e = energy(k)
                  y = e/deli
                  xc_cli=xc*fli(i,j)/deli*(.15+.28*log(y))*de
               else
                  if(k+ieji.gt.nebins)return
                  e = energy(k+ieji)
                  yp= e/deli
                  xc_cli=xc*fli(i,j)/deli*(.15+.28*log(yp))*de*glbygu
               endif
            endif
c     
c  test for delta l > 2 or more
c
         else if (abs(ll(i)-ll(j)) .le. 2) then
            if(i.eq.ir)then
               xc_cli=xc*f1(i,j)/deli*0.15*de
            else
               if(k+ieji.gt.nebins)return
               xc_cli=xc*f1(i,j)/deli*0.15*de*glbygu
            endif
         endif
         return
      endif
c
c  for the case of bounded detailed to bounded non-detailed
c

      if((i.ge.1.and.i.le.ndt).and.(j.ge.ndt+1.and.j.le.nbn))then
         if (i .ge. j)return
         if(i.eq.ir)then
            xc_cli= xc*fli(i,j)/deli*0.38*de
         else
            if(k+ieji.gt.nebins)return
            xc_cli=xc*fli(i,j)/deli*0.38*de*glbygu
         endif
         return
      endif
c
c  for the case non-detailed to non-detailed
c
      if((i.ge.ndt+1.and.i.le.nbn).and.(j.ge.ndt+1.and.j.le.nbn))then
         if (i .ge. j)return
          if(i.eq.ir)then
             xc_cli= xc*fli(i,j)/deli*0.38*de
          else
             if(k+ieji.gt.nebins)return
             xc_cli=xc*fli(i,j)/deli*0.38*de*glbygu
          endif
          return
       endif

c
c  auto to atuo
c
c  test if auto states are used
c
900   if ((nlv-nbn) .ne. 6) return
c
      ccc = 7.067e-8/(atomz-1.5)**2

      if((i.ge.nbn+1.and.i.le.nlv).and.(j.ge.nbn+1.and.j.le.nlv))then
         if (i .ge. j)return
         mi = i-nbn
         mj = j-nbn
         y=energy(k)/deli
         partex = 0.0
         if(i.eq.jr.and.k+ieji.gt.nebins)return
         yp = energy(k+ieji)/deli
         go to (300,400,500,600,700), mi
 300     continue
c
c___choose upper state.
c
         go to (304,305,306,307,308), (mj-1)
c
 304     partde = xcsliaut(z,1)
         if (y .ge. 1.) partex = xcsliaut(y,1)
         go to 1000
 305     partde = xcsliaut(yp,2)
         if (y .ge. 1.) partex = xcsliaut(y,2)
         go to 1000
 306     partde = xcsliaut(yp,3)
         if (y .ge. 1.) partex = xcsliaut(y,3)
         go to 1000
 307     partde = 0.0
         go to 1000
 308     partde = xcsliau1(yp,5)
         if (y .ge. 1.) partex = xcsliau1(y,5)
         go to 1000
c
c___1s2p (singlet p) 2s doublet p (mix)
c
 400     continue
c
c___choose upper state.
c
         go to (405,406,407,408), (mj-2)
c
 405     partde = 0.0
         go to 1000
 406     partde = xcsliaut(yp,7)
         if (y .ge. 1.) partex = xcsliaut(y,7)
         go to 1000
 407     partde = xcsliaut(yp,8)
         if (y .ge. 1.) partex = xcsliaut(y,8)
         go to 1000
 408     partde = xcsliaut(yp,9)
         if (y .ge. 1.) partex = xcsliaut(y,9)
         go to 1000
c
c___1s2p (triplet p) 2s doublet p (mix)
c
 500     continue
c
c___choose upper state.
c
         go to (506,507,508), (mj-3)
c
 506     partde = xcsliaut(yp,10)
         if (y .ge. 1.) partex = xcsliaut(y,10)
         go to 1000
 507     partde = xcsliaut(yp,11)
         if (y .ge. 1.) partex = xcsliaut(y,11)
         go to 1000
 508     partde = xcsliaut(yp,12)
         if (y .ge. 1.) partex = xcsliaut(y,12)
         go to 1000
c
c___1s2p2 doublet d
c
 600     continue
c
c___choose upper state.
c
         go to (607,608), (mj-4)
c
 607     partde = 0.0
         go to 1000
 608     partde = xcsliaut(yp,14)
         if (y .ge. 1.) partex = xcsliaut(y,14)
         go to 1000
c
c___1s2p2 doublet p
c
 700     continue
c
c___only one upper state.
c
         partde = 0.0
         go to 1000
c
c
 1000    continue
         if(i.eq.ir)then
            xc_cli=ccc*partex*de/gli(i)
         else
            if(k+ieji.gt.nebins)return
            xc_cli=ccc*partde*de/gli(j)
         endif
         return
      endif
 
c
c   inner to auto
c
      ccc = 7.067e-8/(atomz-1.5)**2
      if((i.eq.1.or.i.eq.2).and.(j.ge.nbn+1.and.j.le.nlv))then
         m = j-nbn
         if(i.eq.ir)then
            y=energy(k)/deli
            part = sd1(m,i)*xcsamd(y,1)+sd2(m,i)*xcsamd(y,2)
     *           +se1(m,i)*xcsamex(y,1)+se2(m,i)*xcsamex(y,2)
            xc_cli=part*ccc*de*gli(i)
         else
            if(k+ieji.gt.nebins)return         
            yp=energy(k+ieji)/deli
            part = sd1(m,i)*xcsamd(yp,1)+sd2(m,i)*xcsamd(yp,2)
     *           +se1(m,i)*xcsamex(yp,1)+se2(m,i)*xcsamex(yp,2)
            xc_cli=part*ccc*de*gli(j)
         endif
         return
      endif

      return
      end

c---------------------------------------------------
      function xc_chl(k,ir,jr,iz)
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

      naatot=nlev(iz)
      i=min(ir,jr)
      j=max(ir,jr)

      fij=fhul(i,j,iz)
      if(fij.le.0.d0) return

      ii=indlev(i,iz)
      jj=indlev(j,iz)
      gl=glev(ii)
      gu=glev(jj)
      eji = elev(jj)-elev(ii)
      if(eji.le.exthr)return


      ieji = indxc(eji)
      if(ieji.lt.1.or.ieji.gt.nebins)return

      gij = gl/gu
      if(i.eq.ir)then
         e = energy(k)
      else
         if(k+ieji.gt.nebins)return
         e = energy(k+ieji)
      endif
      de = denergy(k)

      ntrn=idcxhl(i,j,iz)
      if(j.le.nlbn(iz).and.ntrn.ne.0)then
         do ix=1,5
            tfit(ix)=chul(ntrn,ix)
         enddo
         xx=e/de
         ocx=tfit(1)*log(xx)+tfit(2)+tfit(3)/(xx+tfit(5)) 
     &        + tfit(4)/(xx+tfit(5))**2
         
         crxx=conxx/gl*ocx
         if(i.eq.ir)then
            xc_chl=crxx*de
         else
            xc_chl=gij*crxx*de
         endif
      else
         g=0.15+0.28*log(e/eji)
         if(i.eq.ir)then
            xc_chl=1.40d-5*fij/eji*g*de
         else
            xc_chl=gij*1.40d-5*fij/eji*g*de
         endif
      endif

      return
      end
c
c----------------------------------------------------------------
c

