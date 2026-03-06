      subroutine atolev

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

      parameter(mpq=10)
      common /sclvnn/nshell

      dimension epion(0:miso)
      dimension noc1(30),noc2(30)
      character*8 knam(40)

      character filnam*10
      logical around

      character*80  msg
       
c----------------------------------------------------------
c
c
      call readatom

      ngas=mgas
      call gauleg(-1.d0,1.d0,xgas,wgas,ngas)


      znuc=atomz
      call setfnn (znuc,fnn,ennn) ! needed for DR contributions

      
      epion(iatomz)=0.d0
      do iz=iatomz,1,-1
         epion(iz-1)=epion(iz)+eipz(iz)
      enddo


c     setup nlev, numpop, index

       nlev(0)=1
       nlbn(0)=1
       ntot=1
       efromgr(ntot)=epion(0)
       elev(ntot)=0.d0
       glev(ntot)=1.d0
       qlev(ntot)=atomz+1.
       isolev(ntot)=0
       ionstg(ntot)=iatomz
       levl(ntot)=1
       indlev(1,0)=ntot       
       kname(ntot)=namaa(1,0)
       mobtot(ntot)=0
       lmatch(ntot)=1
       do ip=1,nprim+1
          noctot(1,ip)=0
       enddo

       ntran=0 ! ntran is computed in sthul or at the end.

       imin=1
       do iso=imin,iatomz
          nlev(iso)=naabn(iso)+naaut(iso)
          nlbn(iso)=naabn(iso)
          do i=1, nlev(iso)
             ntot=ntot+1
             efromgr(ntot)=epion(iso)+evaa(i,iso)
             elev(ntot)=efromgr(ntot)-epion(0)
             glev(ntot)=gaa(i,iso)
             qlev(ntot)=sqq(i,iso)
             isolev(ntot)=iso
             ionstg(ntot)=iatomz-iso
             levl(ntot)=i
             indlev(i,iso)=ntot
             kname(ntot)=namaa(i,iso)
             lvcont(ntot)=indlev(1,iso-1)-1+ncont(i,iso)
             evtoip(ntot)=efromgr(lvcont(ntot))-efromgr(ntot)
             ev2ip(ntot)=evtoip(ntot)
             icont=lvcont(ntot)-indlev(1,iso-1)+1
             mobtot(ntot)=morb(i,iso)
             do ip=1,nprim
                noctot(ntot,ip)=nocel(i,iso,ip)
             enddo
             noctot(ntot,nprim+1)=0
             lmatch(ntot)=i
          enddo
       enddo
       numpop=ntot
       do i=0,iatomz
          nboun(i) = nlbn(i)
       enddo
       if(numpop.gt.numpp) then
          write(ioo,'(" Exiting:too many levels to solve ",i5,">",i5)')
     &         numpop,numpp
          write(16,'(" Exiting:too many levels to solve ",i5,">",i5)')
     &         numpop,numpp
          stop
       endif
      

c set up ionization trasition array  itrn(i,j,iz)
       imin=1
       do iz=imin, iatomz
          do i=1, nlev(iz)
             ilv1=indlev(i,iz)
             do ip=1,nprim
                noc1(ip)=noctot(ilv1,ip)
             enddo
             do j=1, nlev(iz-1)
                ilv2=indlev(j,iz-1)
                do ip=1,nprim
                   noc2(ip)=noctot(ilv2,ip)
                enddo
                petrn(i,j,iz)=0.
                istrn(i,j,iz)=0
                call getlnpqs(nprim,noc1,noc2,n1,n2)
                if(mobtot(ilv1).gt.nprim)then
                   if(mobtot(ilv2).eq.mobtot(ilv1))then
                      if(n1.gt.0.and.n2.eq.0) then
                         petrn(i,j,iz)=noc1(n1)*1.d0
                         istrn(i,j,iz)=n1
                      endif
                   else
                      if(n1.eq.0.and.n2.eq.0) then
                         petrn(i,j,iz)=1.
                         istrn(i,j,iz)=mobtot(ilv1)
                      endif
                   endif
                else
                   if(mobtot(ilv2).le.nprim)then
                      if(n1.gt.0.and.n2.eq.0) then
                         petrn(i,j,iz)=noc1(n1) *1.d0
                         istrn(i,j,iz)=n1 ! of the ionizing orbital
                      endif
                   endif
                endif
             enddo
          enddo
       enddo

c----------------------------------------------------------.
c     Output to files
c----------------------------------------------------------.
c     making a-rate arrays
      ipop=0
      do iz=iatomz,0,-1
         do il=1,nlev(iz)
            ipop=ipop+1
            indpop(ipop)=indlev(il,iz)
            if(il.eq.1)iptlev(iz)=ipop
         enddo
      enddo

      do iz=1, iatomz 
         nlv = nlev(iz)
         nbn = nlbn(iz)
         nau = nlev(iz)-nlbn(iz)
         write (21,1000) nlv,nbn,nau,iptlev(iz)-1,(iatomz-iz),1
         write(21,'("enot",2i5,1x,2f16.5)')iz,nprim,eipz(iz)
         do i=1, nlv
            ilv=indlev(i,iz)
            write(21,101)iz,i,kname(ilv),evaa(i,iz),glev(ilv),
     &           qlev(ilv),(noctot(ilv,ip),ip=1,nprim),mobtot(ilv)

            do j=i+1, nlv
               oratexx(j,i,iz)=orateaa(j,i,iz)
               fxx(i,j,iz)=faa(i,j,iz)
            enddo
         enddo
         do j=2,nlv
            write (21,1010) (oratexx(j,i,iz),i=1,j-1)
         enddo
         write(21,1010) (gamrat(i,iz),i=1,nlv)
      enddo
 101  format('elev', 2i4,1x,a8,1x,f14.4,1x,1p,e16.8,1x,e12.4,
     &     0p,1x,10i3,1x,i4)

c--------------------------------
c     set up transition lists
c--------------------------------

      imin=1
      n1=0
      n2=0
      n3=0
       do iso=imin, iatomz
          do 10 i=2,nlev(iso)
             lup=indlev(i,iso)
             do 20 j=1, i-1
                llo=indlev(j,iso)
                if(oratexx(i,j,iso).gt.0.d0)then
                   ntran=ntran+1
                   llower(ntran)=llo
                   lupper(ntran)=lup
                   ltype(ntran)=1
                   n1=n1+1
                endif
 20          continue
 10       continue

          do 50 j=1, nlev(iso-1)
             lup=indlev(j,iso-1)
             do 40 i=1, nlev(iso)
                llo=indlev(i,iso)
                if(petrn(i,j,iso).gt.0.)then
                   ntran=ntran+1
                   llower(ntran)=llo
                   lupper(ntran)=lup
                   ltype(ntran)=2
                   n2=n2+1
                endif
 40          continue
 50       continue

          do 80 j=1, nlev(iso-1)
             lup=indlev(j,iso-1)
             do 70 i=nlbn(iso)+1, nlev(iso)
                llo=indlev(i,iso)
                if(autoxx(i,j,iso).gt.0.d0)then
                   ntran=ntran+1
                   llower(ntran)=llo
                   lupper(ntran)=lup
                   ltype(ntran)=3
                   n3=n3+1
                endif
 70          continue
 80       continue
       enddo

       write(ioo,'(" # of levels and transitions used",i5,i8)')
     &      numpop, ntran
       write(16,'(" # of levels and transitions used",i5,i8)')
     &      numpop, ntran
       if(ntran.gt.mtran) then
          write(ioo,'(" *** ntran excceds mtran ",i8,i8)')
     &         ntran,mtran
          write(16,'(" Exiting:ntran excceds mtran ",i8,i8)')
     &         ntran,mtran
       endif

c     DO NOT PRINT OUT scfly data file

       return

c     write out energy levels and transitions
       write(filnam,'("scfly.",i2.2)')iatomz
       open(25, file=filnam)

       do i=1, iatomz
          write(25,'("# enot = ",i3,2i5,1x,f16.6)')
     &         i,nlev(i),nlbn(i),eipz(i)
          do j=1, nlev(i)
             ilev=indlev(j,i)
             ilv0=indlev(ncont(j,i),i-1)
             egr=elev(indlev(1,i))
             write(25,'(3i5,3x,a8,1x,f16.6,1p,4e12.5,1x,2i4,1x,30i2)')
     &            i,j,ilev,kname(ilev),
     &            elev(ilev)-egr,glev(ilev),qlev(ilev),
     &            (elev(ilv0)-elev(ilev)),gamrat(j,i),
     &            ncont(j,i),mobtot(ilev),
     &            (noctot(ilev,k),k=1,10)
          enddo
       enddo

       write(25,'("# ntran = ",i7)') ntran
       do i=1, ntran
          llo=llower(i)
          lup=lupper(i)
          iso=isolev(llo)
          ilo=levl(llo)
          iup=levl(lup)

          if(ltype(i).eq.1)then
             n1=nloex(ilo,iup,iso)
             n2=nhiex(ilo,iup,iso)
             itran=mpq*(n1-1)-((n1-1)*n1)/2+(n2-n1)
             if(itran.gt.0.and.itran.le.45)then
                write(25,'("b ",2(1x,i2,1x,i4),2i4,i5,1p,8e13.5)')
     &               iso,ilo,iso,iup,n1,n2,itran,
     &               fxx(ilo,iup,iso),oratexx(iup,ilo,iso),
     &               cfaa(ilo,iup,iso),(cxjj(itran,ix,iso),ix=1,5)
             else
                write(25,'("b ",2(1x,i2,1x,i4),2i4,i5,1p,8e13.5)')
     &               iso,ilo,iso,iup,n1,n2,itran,
     &               fxx(ilo,iup,iso),oratexx(iup,ilo,iso),
     &               cfaa(ilo,iup,iso),(0.d0,ix=1,5)
             endif
          else if(ltype(i).eq.2)then
             ish=istrn(ilo,iup,iso)
             pn =petrn(ilo,iup,iso)
             if(ish.gt.0.and.ish.lt.10) then
                write(25,'("f ",2(1x,i2,1x,i4),1x,i3,1x,1p,9e13.5)')
     &               iso,ilo,iso-1,iup,ish,pn,
     &               (cphilsh(ii,ish,iso),ii=1,8)
             else
                write(25,'("f ",2(1x,i2,1x,i4),1x,i3,1x,1p,9e13.5)')
     &               iso,ilo,iso-1,iup,ish,pn,(0.d0,ii=1,8)
             endif
          else if(ltype(i).eq.3)then
             write(25,'("a ",2(1x,i2,1x,i4),1x,1p,e13.5)')
     &            iso,ilo,iso-1,iup,autoxx(ilo,iup,iso)
          endif
       enddo

       close(25)


      return


 1000 format (12i10)
 1010 format (1p,10e12.5)
      
      end

      subroutine nmatch(noc,imatch,nlv,iz)
      implicit real*8(a-h,o-z)
      include 'mainstuf'
      dimension noc(70,nprim)
      dimension nocion(30),noclev(30),imatch(70)

      do 100 i=1,nlv
         do ip=1,nprim
            nocion(ip)=noc(i,ip)
         enddo
         imatch(i)=0
         do j=1,naabn(iz)+naaut(iz)
            do ip=1,nprim
               noclev(ip)=nocel(j,iz,ip)
            enddo
            call getlnpqs(nprim,noclev,nocion,n1,n2)
            if(n1.eq.0.and.n2.eq.0)then
               imatch(i)=j
            endif
         enddo
 100  continue

      return
      end

      subroutine wratmdat(size,rationNZ)
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'runstuf'
      include 'timestuf'
      include 'popstuf'
      include 'xc_stuf'

      open(24,file='atomic.data')
      write(24,908) tev, denow, trev

c energy levels
       do i=1, iatomz
          write(24,909) i, eipz(i)
          ilv0=indlev(1,i)
          egr=elev(ilv0)
          do j=1, nlev(i)
             ilev=indlev(j,i)
             write(24,910) i,j,kname(ilev),elev(ilev)-egr,glev(ilev),
     &            (noctot(ilev,k),k=1,10),mobtot(ilev)
          enddo
       enddo
       write(24,918)

c oscillator strength
       write(24,912)
       do i=1, iatomz
          do j=1, nlev(i)-1
             jlv=indlev(j,i)
             do k=j+1, nlev(i)
                klv=indlev(k,i)
                if(fxx(j,k,i).gt.0.d0.and. elev(klv).ne.elev(jlv))then
                    write(24,911)i,j,i,k,fxx(j,k,i),
     &               12398.d0/(elev(klv)-elev(jlv)),elev(klv)-elev(jlv)
                 endif
             enddo
          enddo
          write(24,*)
       enddo
       write(24,918)

c spontaneous emission

       write(24,913)
       if(radflag.ne.'off') call avnow(size*ratioNZ)

       do i=1, iatomz
          do j=2, nlev(i)
             do k=1, j
                if(oratexx(j,k,i).gt.1.d-30) then
                   if(radflag.ne.'off') then
                      write(24,911)i,j,i,k,
     &                  oratexx(j,k,i),bradxx(j,k,i),bradxx(k,j,i),
     &                     radint(j,k,i), radint(k,j,i)
                   else
                      write(24,911)i,j,i,k,
     &                  oratexx(j,k,i)
                   endif
                endif
             enddo
          enddo
       enddo

c photoionization
       write(24,914)
       do i=1, iatomz
          do k=1, nlev(i-1)
             klv=indlev(k,i-1)
             do j=1, nlev(i)
                jlv=indlev(j,i)
                eip=elev(klv)-elev(jlv)
                if(eip.gt.1.e-30.and.petrn(j,k,i).gt.0.) then
                   if(radflag.ne.'off') then
                      write(24,911)i,j,i-1,k,eip,phrec(j,k,i)*denow,
     &                     Valfxx(j,k,i)*denow,phion(j,k,i)
                   else
                      write(24,911)i,j,i-1,k,eip,phrec(j,k,i)*denow,
     &                     Valfxx(j,k,i)*denow
                   endif
                endif
             enddo
          enddo
       enddo

c collisional excitation
       write(24,915)
       do i=1, iatomz
          call Vcaafill(i)
          do k=1, nlev(i)-1
             do j=k+1, nlev(i)
                if(fxx(k,j,i).gt.0.d0) then
                   write(24,911)i,k,i,j,
     &                  cxx(k,j,i)*denow,cxx(j,k,i)*denow,
     &                  Vcxx(k,j,i)*denow,Vcxx(j,k,i)*denow
                endif
             enddo
          enddo
       enddo

c collisional ionization
       write(24,916)
       do i=1, iatomz
          k=1
          do j=1, nlev(i)
             jlv=indlev(j,i)
             eip=elev(klv)-elev(jlv)
             if(eip.gt.1.e-30.and.petrn(j,k,i).gt.0.) then
                write(24,911)i,j,i-1,k,
     &               sxx(j,k,i)*denow,betaxx(j,k,i)*denow**2,
     &               Vsxx(j,k,i)*denow,Vbetaxx(j,k,i)*denow**2
             endif
          enddo
          do j=1, nlev(i)
             jlv=indlev(j,i)
             do k=2, nlev(i-1)
                klv=indlev(k,i-1)
                eip=elev(klv)-elev(jlv)
                if(eip.gt.1.e-30.and.petrn(j,k,i).gt.0.) then
                   write(24,911)i,j,i-1,k,
     &               sxx(j,k,i)*denow,betaxx(j,k,i)*denow**2,
     &               Vsxx(j,k,i)*denow,Vbetaxx(j,k,i)*denow**2
                endif
             enddo
          enddo
       enddo

c autoionization
       write(24,917)
       do i=1, iatomz
          do j=1, nlev(i)
             jlv=indlev(j,i)
             do k=1, nlev(i-1)
                klv=indlev(k,i-1)
                eip=elev(jlv)-elev(klv)
                if(eip.gt.1.e-30.and.autoxx(j,k,i).gt.1.d-30) then
                   write(24,911)i,j,i-1,k,
     &                  autoxx(j,k,i),capexx(j,k,i)*denow,
     &                  Vcapexx(j,k,i)*denow
                endif
             enddo
          enddo
       enddo

       close(24)

      

c 905  format('C Atomic data')
 908  format('Te =',f10.2,' Ne= ',1p,e12.4, ' Tr =',f10.2)
 909  format('enot',i5,1x,2f16.5)
 910  format('elev', 2i4,1x,a8,1x,f12.4,1x,1p,e16.8,0p,1x,10i3,1x,i4)
 911  format('d  ', 4i5,1x,1p,10e12.3)
 912  format('data   phxs  ')
 913  format('  rate type:   photoexcitation')
 914  format('  rate type:   photoionization       ')
 915  format('  rate type:   collisional excitation')
 916  format('  rate type:   collisional ionization')
 917  format('  rate type:   augxs')
 918  format('end data')
 920  format(2i4,1p,10e10.2)
      return
      end
c
c
      subroutine avnow(size)
c
c____ this subroutine calculates the escape factor reductions of
c          the a-values.
c     Use voigt formalism and frequency-average and space-integrated
c     escape probability formalism 
c     correction on 2/18/2011

      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'timestuf'
      include 'popstuf'
      include 'runstuf'

      parameter(mphot=27,nphot=127)
      dimension evv(mphot),evf(nphot),phf(nphot),phv(nphot)
      dimension fvf(nphot),fev(nphot)
      dimension dhv(ngrid)

      data const/2.33333E-15/
      data evv /-70.0, -50., -30.,
     *      -15., -8.0, -4.0, -3.0, -2.0, 
     *       -1.5, -1., -0.5, -0.3, -0.15, 
     &     0., 0.15, 0.3, 0.50, 1., 1.500, 
     *     2.00, 3.00, 4.00, 8.00, 15.00, 
     &     30.0, 50.0, 70./ 

c    hceverg = (hc)**3/evtoerg**3.5*sqrt(amass)/(8*pi)
c    load working arrays with a-values

      external field

      iz1=isomin
      iz2=isomax

      do iz=iz1, iz2
         do i=1, nlev(iz)
            do j=1, i-1
               aratexx(i,j,iz)=oratexx(i,j,iz)
               escapexx(i,j,iz)=1.00
            enddo
         enddo
      enddo

c      goto 11  ! skipping the line inegration

c     if radiation field is not plankian, 
c     one needs to integrate over the line shape.


      if(radflag .eq. 'trfile'.and.xbarj.le.0.d0)then

         hvmin=xev(ngrid)
         hvmax=xev(1)
         do i=1, ngrid
            if(field(xev(i)).gt.0.d0)then
               hvmin=min(hvmin,xev(i))
               hvmax=max(hvmax,xev(i))
            endif
            if(i.eq.1)then
               dhv(i)=(xev(2)-xev(1))/2.*evtohz
            else if(i.eq.ngrid)then
               dhv(i)=(xev(ngrid)-xev(ngrid-1))/2.*evtohz
            else
               dhv(i)=(xev(i+1)-xev(i-1))/2.*evtohz
            endif
         enddo

         iwrad=1
         if(iwrad.eq.1)open(40,file='radint.d')

         cwd=4.6337d-5*sqrt(tiev/(2.*atomz))
         do iz=iz1, iz2
            do i=2, nlev(iz)
               ilv=indlev(i,iz)
               if(indl2p(ilv).ne.0)then
                  do j=1,i-1
                     jlv=indlev(j,iz)
                     radint(j,i,iz)=0.d0 ! integrated rates
                     eradint(j,i,iz)=0.d0 ! integrated energy rates
                     fij=fxx(j,i,iz)
                     if (fij .gt. 0.d0.and.indl2p(jlv).ne.0) then
                        eul = elev(ilv)-elev(jlv)
                        wd=cwd*eul
                        avt=(gamrat(i,iz)+gamrat(j,iz))/
     &                       (wd*12.57*evtohz)
                        f=max(avt, 1.d0)
                        evmn=eul+evv(1)*wd*f  ! minimum 
                        evmx=eul+evv(mphot)*wd*f  ! maximum 


                        if(evmn.gt.xev(ngrid).or.evmx.lt.xev(1))then
c             the line profile may not reach out to the radiation field
c             in that case, the rate should be included by picking up 
c             intensities at the radiation field grids:
c             energy grid is defined as the radiation field grid
                           esum=0.d0
                           fsum=0.d0
                           sum=0.d0
c    hyun test on cut-off 04/07/2011--- xcut*xcut =100.* avt*avt or xcut=100
c    the line profile behaves 1/x^2 to xcut and 1/x^4 above xcut. 
                           xcut=max(10.*avt, 100.)
                           cfar=const/wd*voigt(xcut,avt)*xcut**4

                           do jj=1, ngrid
                              xx=(xev(jj)-eul)/wd
                              if(abs(xx).le.xcut)then
                                 phix=const/wd*voigt(xx,avt)
                              else
                                 phix=cfar/xx**4
                              endif
                              esum=esum+phix*barjnow(jj)*dhv(jj)
     &                             *xev(jj)
                              fsum=fsum+phix*barjnow(jj)*dhv(jj)
                              sum=sum+phix*dhv(jj)
                           enddo
                           radint(j,i,iz)=fsum
                           eradint(j,i,iz)=esum
                           if(radint(j,i,iz).gt.0.d0.and.iwrad.eq.1)
     &                          write(40,'(3i3,1p,10e10.2)') iz,i,j,
     &                      eul,sum,radint(j,i,iz),eradint(j,i,iz),
     &                      xx,wd, avt*wd*12.57,avt,xcut    
                        else

                           npts=0
                           do jj=1,mphot
                              ex=eul+evv(jj)*wd*f ! normalized energy grid
                              if(ex.gt.1.d-7)then
                                 npts=npts+1
                                 fev(npts)=ex
                              endif
                           enddo
                           do jj=1,ngrid
                              npts=npts+1
                              fev(npts)=xev(jj) ! add radiation field grid
                           enddo
                           dx1=(xev(2)-xev(1))/100.
                           dx2=(xev(ngrid)-xev(ngrid-1))/100.
                           npts=npts+1
                           fev(npts)=xev(1)-dx1 ! add the left buffer
                           npts=npts+1
                           fev(npts)=xev(ngrid)+dx2 ! add the right buffer
                           call sort(npts,fev)
                           do jj=1,npts
                              ehv=fev(jj)-eul
                              phf(jj)=const/wd*voigt(ehv/wd,avt) !line profile at the energy
                              fvf(jj)=field(fev(jj)) ! radiation field at the given energy
                              evf(jj)=ehv*evtohz ! frequency grid from line-center
                           enddo
                           !calculating line-profile * energy spacing
                           phv(1)=phf(1)*(evf(2)-evf(1))/2.0
                           phv(npts)=phf(npts)*
     &                          (evf(npts)-evf(npts-1))/2.0
                           sum=phv(1)+phv(npts)
                           do jj=2,npts-1
                              phv(jj)=phf(jj)*(evf(jj+1)-evf(jj-1))/2.0
                              sum=sum+phv(jj)
                           enddo
                           !calculating the normalization constant and integration
                           fsum=phv(1)*fvf(1)+
     &                          phv(npts)*fvf(npts)
                           esum=phv(1)*fvf(1)*fev(1)+
     &                          phv(npts)*fvf(npts)*fev(npts)
                           do jj=2,npts-1
                              fsum=fsum+phv(jj)*fvf(jj)
                              esum=esum+phv(jj)*fvf(jj)*fev(jj)
                           enddo
                           radint(j,i,iz)=fsum/sum
                           eradint(j,i,iz)=esum/sum
                           if(radint(j,i,iz).gt.0.d0.and.iwrad.eq.1)
     &                          write(40,'(3i3,1p,10e10.2)') iz,i,j,eul,
     &                          sum,radint(j,i,iz),eradint(j,i,iz)

                        endif
                     endif
                  enddo
               endif
            enddo
         enddo

         if(iwrad.eq.1)close(40)

      endif


 11   if (size .eq. 0.00) return

c
      cwd=4.6337d-5*sqrt(tiev/(2.*atomz))
      do iz=iz1, iz2
         do i=2, nlev(iz)
            ilv=indlev(i,iz)
            pu=popout(ilv)
            if(indl2p(ilv).ne.0)then
               do 100 j=1,i-1
                  jlv=indlev(j,iz)
                  fij=fxx(j,i,iz)
                  pl=popout(jlv)
                  if (fij.gt.0.d0.and.indl2p(jlv).ne.0.and.
     &                 pl.gt.0.d0)then
                     eul = elev(ilv)-elev(jlv)
                     wd=cwd*eul
                     opl=0.0265*fij*pl*(1.0-pu/pl*glev(jlv)/glev(ilv))
                     if(opl.lt.1.d-70) goto 100
                     avt=(gamrat(i,iz)+gamrat(j,iz))/(wd*12.57*evtohz)
                     phi0=const/wd*voigt(0.d0,avt)
                     taun=opl*phi0*size
                     if(taun.gt.1.e-30) then
                        ef=espint(taun,avt)/taun
                     else
                        ef=1.d0
                     endif
                     aratexx(i,j,iz) = aratexx(i,j,iz)*ef
                     escapexx(i,j,iz) = ef
                  endif
 100           continue
            endif
         enddo
      enddo
c
c
      return
      end

      subroutine sort (npts,xin )
      implicit double precision (a-h,o-z)
c ... puts the elements of "xin" into ascending order
      dimension xin(npts)
      logical flag
      jt1=npts-1
      do 100 i=1, jt1
         etmp=xin(i)
         jmin=i
         ip1=i+1
         do 200 j=ip1, npts
            if (xin(j).gt.etmp) goto 200
            etmp=xin(j)
            jmin=j
 200     continue
         xin(jmin)=xin(i)
         xin(i)=etmp
 100  continue

      i=1
      do n=i+1,npts
         if(xin(n).gt.1.00001*xin(i))then
            i=i+1
            xin(i)=xin(n)
         endif
      enddo
      npts=i

      return
      end

      function espint(taun,avt)
      implicit real*8 (a-h,o-z)

      if(avt.eq.0.d0) then

c         if (taun .ge. 2.50) then
c            espint = 1.0/(taun*sqrt(pi*log(taun)))
c         else
c            espint = exp(-taun/1.73)
c         endif

         if ( taun .le. 0.01 ) then
            espint = taun * (1.-0.5*taun)
         else if ( taun .le. 5.18 ) then
            a = 0.675378 * taun + 0.756889
            espint = 2.32889 * ( atan(a) - 0.6478955 )
         else
            taulog = log( taun )
            espint = 0.209213 + 1.09369 * sqrt( taulog )
         endif

      else if ( avt .lt. 0.49 ) then

         tauc = 0.83 / ( avt*(1.+sqrt(avt)) )
         if ( taun .le. 0.01 ) then
            espint = taun * (1.-0.5*taun)
         else if ( taun.le.1. ) then
            espint = 0.666667 * log( 1.+1.5*taun )
         else if ( taun.le.tauc ) then
            espint = 0.61086 + 0.4*log( taun )
         else
            espint = 0.61086 + 0.4*log( tauc ) + 
     &               0.8 * ( sqrt( taun/tauc ) - 1. )
         endif

      else

         if ( taun .le. 0.01 ) then
            espint = taun * (1.-0.5*taun)
         else if ( taun.le.1. ) then
            espint = log( 1.+taun )
         else
            espint = sqrt( taun ) - 0.30685
         endif

      endif

      return
      end


c
c
c          ******   ******     ***    ******      * *      * *   
c           *    *   *    *   *   *    *    *   *     *  *     * 
c           *    *   *    *  *     *   *    *   *     *  *     * 
c           *****    *****   *******   *    *   *******  ******* 
c           *    *   *  *    *     *   *    *   *     *  *     * 
c           *    *   *   *   *     *   *    *   *     *  *     * 
c          ******    *    *  *     *  ******    *     *  *     * 
c
c
c
c
      function bradxx(iq,jq,iz)
c
      implicit real*8 (a-h,o-z)
c___ this calculates the photo-excitation 
c___ from hydrogenic state i to j.
c
      include 'mainstuf'
      include 'popstuf'
      include 'runstuf'

      common /paramnow/ trnow,tinow,zbarnow,totpnow,sizenow,rhonow
      data  const/2.33333E-15/

c
      bradxx=0.
      atrad=0.
      etrad=0.

      i=min(iq,jq)
      j=max(iq,jq)
      ilv1=indlev(i,iz)
      ilv2=indlev(j,iz)
      eaa=elev(ilv2)-elev(ilv1)
      if(abs(eaa).le.1.d-10) return
      fij=fxx(i,j,iz)
      if (radflag .ne. 'trfile') then
        if(trev.le.0.d0)return
        ebt=eaa/trev
        atrad = 2.08e+11/eaa*fij*bbfield(ebt)*escapexx(j,i,iz)
        etrad = atrad*eaa
      else 
c     radint is radiation field integrated over line profile
c     defined in avnow subroutine
         if(xbarj.gt.0.d0)then
            atrad = 2.08e+11/eaa*fij*field(eaa)*escapexx(j,i,iz)
            etrad = atrad*eaa
         else
            atrad = 2.08e+11/eaa*fij*radint(i,j,iz)*escapexx(j,i,iz)
            etrad = 2.08e+11/eaa*fij*eradint(i,j,iz)*escapexx(j,i,iz)
         endif
      endif

      if(i.eq.iq)then
         bradxx=atrad
      else
         bradxx=atrad*glev(ilv1)/glev(ilv2)
      endif

      return
      end

c
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                     *
*       EEEEEE  RRRRR   FFFFFF                                        *
*       EEEEEE  RR  RR  FFFFFF                                        *
*       EE      RR  RR  FF                                            *
*       EEEEEE  RRRRR   FFFFFF                                        *
*       EEEEEE  RRRR    FFFFFF                                        *
*       EE      RR RR   FF                                            *
*       EEEEEE  RR  RR  FF                                            *
*       EEEEEE  RR  RR  FF                                            *
*                                                                     *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      real*8 function erf(x)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  error function
c   erf(x)=2.0/sqrt(3.1415)*(integral of exp(-y**2) from y=0.0 to y=x)
c
c  reference: "handbook of mathematical functions" ,
c   abramowitz and stegun , national bureau of standards .
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit real*8 (a-h,o-z)

      dimension a(6)

      data(a(i),i=1,6) / 0.0705230784,0.0422820123,0.0092705272,
     10.0001520143,0.0002765672,0.0000430638 /

      f=1.0

      do 100 i=1,6
        f=f+a(i)*x**i
  100 continue

      erf=1.0-1.0/f**16

      return
      end
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                     *
*       EEEEEE  RRRRR   FFFFFF    CC                                  *
*       EEEEEE  RR  RR  FFFFFF   CCCC                                 *
*       EE      RR  RR  FF      CC  CC                                *
*       EEEEEE  RRRRR   FFFFFF  CC                                    *
*       EEEEEE  RRRR    FFFFFF  CC                                    *
*       EE      RR RR   FF      CC  CC                                *
*       EEEEEE  RR  RR  FF       CCCC                                 *
*       EEEEEE  RR  RR  FF        CC                                  *
*                                                                     *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      real*8 function erfc(y)

c     .. specifications for arguments

      real*8               y

c     .. specifications for local variables

      integer            isw,i

      dimension          p(5),q(3),p1(8),q1(7),p2(5),q2(4)

      real*8               p,q,p1,q1,p2,q2,xmin,xlarge,ssqpi,x,
     1                   res,xsq,xnum,xden,xi

c     .. coefficients for 0.0 .le. y .lt.
c     ..  477

      data               p(1)/-.44422647396874/,
     1                   p(2)/10.731707253648/,
     1                   p(3)/15.915606197771/,
     1                   p(4)/374.81624081284/,
     1                   p(5)/2.5612422994823d-02/
      data               q(1)/17.903143558843/,
     1                   q(2)/124.82892031581/,
     1                   q(3)/332.17224470532/

c     .. coefficients for .477 .le. y
c     ..  le. 4.0

      data               p1(1)/7.2117582508831/,
     1                   p1(2)/43.162227222057/,
     1                   p1(3)/152.98928504694/,
     1                   p1(4)/339.32081673434/,
     1                   p1(5)/451.91895371187/,
     1                   p1(6)/300.45926102016/,
     1                   p1(7)/-1.3686485738272d-07/,
     1                   p1(8)/.56419551747897/
      data               q1(1)/77.000152935229/,
     1                   q1(2)/277.58544474399/,
     1                   q1(3)/638.98026446563/,
     1                   q1(4)/931.35409485061/,
     1                   q1(5)/790.95092532790/,
     1                   q1(6)/300.45926095698/,
     1                   q1(7)/12.782727319629/

c     .. coefficients for 4.0 .lt. y

      data               p2(1)/-.22695659353969/,
     1                   p2(2)/-4.9473091062325d-02/,
     1                   p2(3)/-2.9961070770354d-03/,
     1                   p2(4)/-2.2319245973418d-02/,
     1                   p2(5)/-2.7866130860965d-01/
      data               q2(1)/1.0516751070679/,
     1                   q2(2)/.19130892610783/,
     1                   q2(3)/1.0620923052847d-02/,
     1                   q2(4)/1.9873320181714/

c     .. constants

      data               xmin/1.0e-8/,xlarge/5.6875/

c     .. erfc(xbig) .approx. setap

      data               xbig/25.90625/
      data               ssqpi/.56418958354776/

c     .. first executable statement

      x = y
      isw = 1

      if (x.ge.0.0e0) go to 5

      isw = -1
      x = -x

    5 if (x.lt..477e0) go to 10
      if (x.le.4.0e0) go to 30
      if (isw .gt. 0) go to 40
      if (x.lt.xlarge) go to 45

      res = 2.0e0
      go to 65

c     .. abs(y) .lt. .477, evaluate
c     .. approximation for erfc

   10 if (x.lt.xmin) go to 20

      xsq = x*x
      xnum = p(5)

      do 15 i = 1,4
        xnum = xnum*xsq+p(i)
   15 continue

      xden = ((q(1)+xsq)*xsq+q(2))*xsq+q(3)
      res = x*xnum/xden
      go to 25

   20 res = x*p(4)/q(3)

   25 if (isw.eq.-1) res = -res

      res = 1.0e0-res
      go to 65

c     ..  477 .le. abs(y) .le. 4.0
c     .. evaluate approximation for erfc

   30 xsq = x*x
      xnum = p1(7)*x+p1(8)
      xden = x+q1(7)

      do 35 i=1,6
        xnum = xnum*x+p1(i)
        xden = xden*x+q1(i)
   35 continue

      res = xnum/xden
      go to 55

c     .. 4.0 .lt. abs(y), evaluate
c     .. minimax approximation for erfc

   40 if (x.gt.xbig) go to 60

   45 xsq = x*x
      xi = 1.0e0/xsq
      xnum = p2(4)*xi+p2(5)
      xden = xi+q2(4)

      do 50 i = 1,3
        xnum = xnum*xi+p2(i)
        xden = xden*xi+q2(i)
   50 continue

      res = (ssqpi+xi*xnum/xden)/x

   55 res = res*exp(-xsq)

      if (isw.eq.-1) res = 2.0e0-res
      go to 65

   60 res = 0.0e0

   65 erfc = res

      return
      end
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                     *
*       EEEEEE  X    X  PPPPP   EEEEEE     1                          *
*       EEEEEE  X    X  PP  PP  EEEEEE    11                          *
*       EE       X  X   PP  PP  EE       111                          *
*       EEEEEE    XX    PPPPP   EEEEEE    11                          *
*       EEEEEE    XX    PP      EEEEEE    11                          *
*       EE       X  X   PP      EE        11                          *
*       EEEEEE  X    X  PP      EEEEEE    11                          *
*       EEEEEE  X    X  PP      EEEEEE   1111                         *
*                                                                     *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      real*8 function expe1 (px)

c     .. rev: 6/15/78

c     .. calculate  e1(px),  the exponential integral
c     .. e1(x) = integral from x to infinity of (exp(-t)/t)dt.

c     .. reference: "handbook of mathematical functions", abramowitz
c     .. and stegun, national bureau of standards

      implicit real*8 (a-h,o-z)

      if(px.ge.1.)go to 2
      if(px.gt.1.e-4)go to 1

      expe1 =   ( - dlog(px) - .57721566 + px )

      return

    1 expe1 = ( - dlog(px) - .57721566 + .99999193*px
     1         - .24991055*px**2 + .05519968*px**3 - .00976004*px**4
     1         + .00107857*px**5 )

      return

    2 expe1 = (px + 2.334733 + .250621/px) /
     1         (px**2 + 3.330657*px + 1.681534)
      expe1=expe1/exp(px)

      return
      end

c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c*                                                                     *
c*       EEEEEE  X    X  PPPPP   UU  UU  N    N  DDDD                  *
c*       EEEEEE  X    X  PP  PP  UU  UU  NN   N  DDDDD                 *
c*       EE       X  X   PP  PP  UU  UU  NNN  N  DD  DD                *
c*       EEEEEE    XX    PPPPP   UU  UU  N NN N  DD  DD                *
c*       EEEEEE    XX    PP      UU  UU  N NN N  DD  DD                *
c*       EE       X  X   PP      UU  UU  N  NNN  DD  DD                *
c*       EEEEEE  X    X  PP      UU  UU  N   NN  DDDDD                 *
c*       EEEEEE  X    X  PP       UUUU   N    N  DDDD                  *
c*                                                                     *
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      real*8 function expund (px,pchk)

      implicit real*8 (a-h,o-z)

c     .. prevent underflows in rate calculations by returning
c     .. 0. for exp(px) if px < pchk.

      if( px - pchk ) 1 , 1 , 2

    1 expund = exp(px)

      return

    2 expund = exp(px)

      return
      end



c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      function voigt (vv,a)
      implicit real*8 (a-h,o-z)

c ... this subroutine computes the Voigt function "h" as a function
c     of "a" and "v".  See Mihalas, "Stellar Atmospheres", 1978, for
c     definitions and details.
      dimension atab(11),vtab(11),htable(11,11)
      dimension table1(21),table2(11)

      data atab /  1.00E-01,1.58E-01,2.51E-01,3.98E-01,6.31E-01,
     2    1.00E+00,1.58E+00,2.51E+00,3.98E+00,6.31E+00,1.00E+01 /
      data vtab /  0.00E+00,5.00E-01,1.00E+00,1.50E+00,2.00E+00,
     2    2.50E+00,3.00E+00,3.50E+00,4.00E+00,4.50E+00,5.00E+00 /
      data htable /
     2  8.96E-01,8.44E-01,7.69E-01,6.72E-01,5.54E-01,4.28E-01,3.08E-01,
     3 2.10E-01,1.38E-01,8.83E-02,5.61E-02, 7.18E-01,6.85E-01,6.38E-01,
     4 5.72E-01,4.89E-01,3.91E-01,2.92E-01,2.04E-01,1.36E-01,8.78E-02,
     5 5.60E-02, 3.73E-01,3.74E-01,3.72E-01,3.63E-01,3.43E-01,3.05E-01,
     6 2.50E-01,1.88E-01,1.30E-01,8.63E-02,5.56E-02, 1.34E-01,1.48E-01,
     7 1.66E-01,1.87E-01,2.05E-01,2.12E-01,1.98E-01,1.65E-01,1.22E-01,
     8 8.39E-02,5.49E-02, 4.02E-02,5.18E-02,6.85E-02,9.07E-02,1.17E-01,
     9 1.40E-01,1.51E-01,1.40E-01,1.12E-01,8.07E-02,5.40E-02, 1.47E-02,
     1 2.19E-02,3.28E-02,4.86E-02,6.97E-02,9.38E-02,1.13E-01,1.17E-01,
     2 1.02E-01,7.69E-02,5.29E-02, 7.94E-03,1.25E-02,1.95E-02,3.01E-02,
     3 4.55E-02,6.53E-02,8.53E-02,9.64E-02,9.11E-02,7.27E-02,5.16E-02,
     4  5.34E-03,8.44E-03,1.33E-02,2.08E-02,3.21E-02,4.77E-02,6.58E-02,
     5 7.97E-02,8.09E-02,6.83E-02,5.01E-02, 3.92E-03,6.21E-03,9.81E-03,
     6 1.54E-02,2.40E-02,3.63E-02,5.18E-02,6.62E-02,7.16E-02,6.39E-02,
     7 4.85E-02, 3.02E-03,4.79E-03,7.57E-03,1.19E-02,1.86E-02,2.85E-02,
     8 4.17E-02,5.55E-02,6.33E-02,5.94E-02,4.69E-02, 2.41E-03,3.81E-03,
     9 6.03E-03,9.52E-03,1.49E-02,2.30E-02,3.42E-02,4.69E-02,5.59E-02,
     1 5.51E-02,4.51E-02 /

c ... "table1" and "table2" are used to compute Dawson's integral 
c     (which is used to compute the Voigt function for small "a"):
c                               _v
c       D.I.  =  exp(-v**2) * _/ dt exp(t**2)
c                            0
c
      data table1 / .00000, .09933, .19475, .28263, .35994,
     2      .42443, .47476, .51050, .53210, .54072, .53807,
     3      .52620, .50727, .48339, .45650, .42824, .39993,
     4      .37255, .34677, .32297, .30134 /
      data table2 / .50000, .50650, .51358, .52142, .53037,
     2      .54079, .55265, .56547, .57852, .59110, .60268 /

      data sqrtpi / 1.7725 /

      data xxxmax / 20. /


c ********************************************************************
c
c                           begin execution
c
c ********************************************************************

      v = abs( vv )

      if ( a .eq. 0. ) then

c ...    Doppler profile

         voigt = exp( -v*v )

      else if ( v.ge.5. .or. a.ge.10. ) then

c ...    Lorentzian profile in limit of large a,v

         voigt = a / sqrtpi / ( a*a + v*v )

      else if ( a.le.0.1 ) then

c ...    use Dawson's integral (see Mihalas) for small "a"

         if ( v .le. 2. ) then
            x = v
            ix = min ( 1. + x / 0.1 , xxxmax )
            f = x/0.1 - (ix-1)
            y = (1.-f) * table1(ix) + f * table1(ix+1)
            dawson = y
         else
            x = 1. / (v*v)
            ix = 1 + x / 0.025
            f = x/0.025 - (ix-1)
            y = (1.-f) * table2(ix) + f * table2(ix+1)
            dawson = y / v
         endif

         voigt = exp( -v*v ) + 2.*a/sqrtpi *
     &           ( 2.*v*dawson - 1. )

      else

c ...    use a table lookup

         aalog = log10( a )
         ia = aalog/0.2 + 6
         iv = v/0.5 + 1
         fa = aalog/0.2 - (ia-6)
         fv = v/0.5 - (iv-1)
         g1 = (1.-fa)*htable(ia,iv) + fa*htable(ia+1,iv)
         g2 = (1.-fa)*htable(ia,iv+1) + fa*htable(ia+1,iv+1)
         voigt = (1.-fv)*g1 + fv*g2

      endif
      return
      end

      subroutine gauleg(x1,x2,x,w,n)
c-----generate gauss-legendre points and weights
      implicit real*8 (a-h,o-z)
      parameter (eps=3.e-14)
      dimension x(n),w(n)
      data pi/3.141592654/
c
c     write(6,"('...in gauleg')")
      m=(n+1)/2
      xm=0.5*(x2+x1)
      xl=0.5*(x2-x1)
      do 12 i=1,m
         z=cos(pi*(i-0.25)/(n+0.5))
  1      continue
         p1=1.
         p2=0.
         do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.*j-1.)*z*p2-(j-1.)*p3)/j
 11      continue
         pp=n*(z*p1-p2)/(z*z-1.)
         z1=z
         z=z1-p1/pp
         if(abs(z-z1).gt.eps) goto 1
         x(i)=xm-xl*z
         x(n+1-i)=xm+xl*z
         w(i)=2.*xl/((1.-z*z)*pp*pp)
         w(n+1-i)=w(i)
 12   continue
      return
      end

      function atmass(nz)

      implicit real*8 (a-h,o-z)
      include 'mainstuf'
c
c   set atomic weight for atom with atomic number nz
c
      dimension atw(106)
      data atw/1.00,   4.00,   6.94,   9.01,  10.81,  12.01,  14.01,
     1        16.00,  19.00,  20.18,  22.99,  24.31,  26.98,  28.09,
     2        30.97,  32.06,  35.45,  39.95,  39.10,  40.08,  44.96,
     3        47.90,  50.94,  52.00,  54.94,  55.85,  58.93,  58.71,
     4        63.55,  65.37,  69.72,  72.59,  74.92,  78.96,  79.90,
     5        83.80,  85.47,  87.62,  88.91,  91.22,  92.91,  95.94,
     6        98.91, 101.07, 102.91, 106.40, 107.87, 112.40, 114.82,
     7       118.69, 121.75, 127.60, 126.90, 131.30, 132.90, 137.34,
     8       138.91, 140.12, 140.91, 144.24, 145.00, 150.35, 151.96,
     9       157.25, 158.92, 162.50, 164.93, 167.26, 168.93, 173.04,
     a       174.97, 178.49, 180.95, 183.85, 186.20, 190.20, 192.20,
     b       195.09, 196.97, 200.59, 204.37, 207.19, 208.98, 209.00,
     c       210.00, 222.00, 223.00, 226.00, 227.00, 232.04, 231.00,
     d       238.03, 237.00, 244.00, 243.00, 247.00, 247.00, 251.00,
     e       252.00, 257.00, 258.00, 259.00, 260.00, 261.00, 262.00,
     f       263.00/
c
      if (nz .le. 0) then
         atmass = 0.00
         write(ioo,*)' *** atomic mass sort for bad atomic number'
      else if (nz .lt. 106) then
         atmass = atw(nz)*amass
      else
         atmass = nz*2.50*amass
      endif

      return
      end


      subroutine setfnn (znuc,fnn,ennn)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   calculate the average oscillator strength and transition energy for
c
c   n-n transition . this is for ground state only .
c
c      zann(i) i=number of bounded electron
c 
c    See Post et al. Appendix of Atomic Data And Nuclear Data Tables 20, 397 (1977)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit real*8 (a-h,o-z)

      dimension zann(60),zbnn(60),zepsnn(60),zalpnn(60),fnn(100),
     1          ennn(100)

      data (zann(j),j=1,60) / 0., 0.,
     1 -.02, -.01, -.04, .04, -.01, -.01, .01, 0.,
     1  .01, .1, .25, .17, .08, -.02, -.14, -.27,
     1 -.29, -.3, -.3, -.29, -.27, -.24, -.2, -.14,
     1 -.08, 0., .97, 1.96, 1.92, 1.89, 1.86, 1.83,
     1 1.78, 1.73, 1.41, 1.05, .674, .261, -.167,
     1 -.637, -1.14, -1.67, -2.26, -2.88, -2.9,
     1 -2.83, -2.72, -2.61, -2.45, -2.27, -2.05, -1.81,
     1 -1.55, -1.25, -.94, -.63, -.31, 0. /
      data (zbnn(j),j=1,60) / 0., 0.,
     1 2., 4.33, 4.68, 3.6, 2.85, 2.4, 1.3, 0.,
     1 10.5, 20., 16.4, 24.7, 34.1, 44.8, 56.8, 70.3,
     1 68.6, 65.8, 61.9, 56.9, 50.7, 45.5, 34.4,
     1 24.2, 12.8, 0., -25.8, -53.1, -28.1, -2.97,
     1 23., 50.3, 79.1, 109.2, 141.7, 176., 213.9,
     1 256.6, 299.6, 348.2, 400.5, 456.6, 518.1,
     1 584.0, 571.9, 550.7, 524.9, 496.8, 463.4,
     1 427., 385.3, 339.8, 291.3, 237.4, 181.3, 122.9,
     1 62.2, 0. /
      data (zepsnn(j),j=1,28) / 0., 0.,
     1 2.04e-3, 4.49e-3, 6.8e-3, 8.16e-3, 6.80e-3,
     1 1.06e-2, 1.31e-2, 0., 2.3e-3, 3.88e-3,
     1 5.71e-3, 5.44e-3, 8.16e-3, 6.8e-3, 8.3e-3,
     1 4.32e-3, 5.11e-3, 6.04e-3, 7.08e-3,
     1 8.12e-3, 9.69e-3, 1.13e-2, 1.37e-2, 1.49e-2,
     1 1.71e-2, 0. /
      data (zepsnn(j),j=29,60) / 1.52e-5, 1.93e-5,
     1 5.62e-5, 1.20e-4, 2.13e-4, 3.34e-4, 4.66e-4,
     1 6.66e-4, 9.10e-4, 1.25e-3, 1.69e-3, 1.98e-3,
     1 3.03e-3, 3.92e-3, 5.19e-3, 6.40e-3, 8.05e-3,
     1 9.80e-3, 1.08e-2, 1.19e-2, 1.32e-2, 1.47e-2,
     1 1.59e-2, 1.81e-2, 1.92e-2, 2.10e-2, 2.33e-2,
     1 2.57e-2, 2.82e-2, 3.10e-2, 3.43e-3, 0. /
      data (zalpnn(j),j=1,60) / 1.0, 1.0,
     1 1.0, .93, .86, .80, .87, .77, .77, 1.0,
     1 .99, .90, .83, .87, .77, .82, .79, 1.06,
     1 1.03, 1., .97, .94, .91, .88, .85,
     1 .83, .80, 1.0, 2.46, 2.4, 2.14, 1.96,
     1 1.83, 1.72, 1.64, 1.57, 1.49, 1.41, 1.34,
     1 1.30, 1.19, 1.13, 1.06, 1.01, .95, .91,
     1 .89, .87, .85, .83, .82, .80, .78, .76,
     1 .74, .72, .70, .68, .66, 1.0 /

ccccccc  zero oscillator strength array and energy array

      do 5 jq=1,100
        fnn(jq)=0.0
        ennn(jq)=0.0
    5 continue

      nzn=znuc

      do 10 i=1,nzn

        if(i.gt.60) return

        ze=znuc-i
        fnn(i)=zann(i)+zbnn(i)/znuc
        alfa=zalpnn(i)
        ennn(i)=zepsnn(i)*ze**alfa*1.0e3  ! in eV

   10 continue

      return
      end



      subroutine getlnpqs(nprim,noclev,nocion,n1,n2)
      implicit real*8 (a-h,o-z)
      dimension noclev(30),nocion(30)
c ... compare shell occupations to identify principal quantum numbers
c ... for a potential transition

c ---------------------------------------------------------------------
      n1 = 0
      n2 = 0
      do 100 ip = 1,nprim
         ndp = noclev(ip)-nocion(ip)
         if (ndp.eq.1) then
            if (n1.eq.0) then
               n1 = ip
            else
               n1 = -1
            endif
         elseif (ndp.eq.-1) then
            if (n2.eq.0) then
               n2 = ip
            else
               n2 = -1
            endif
         elseif (ndp.ne.0) then
            n1 = -1
         endif
  100 continue
      return
      end

      subroutine readatom
      implicit double precision (a-h,o-z)
      include 'mainstuf'
      include 'popstuf'
      dimension cxin(10), indinp(nlvi,0:miso)
      dimension noc1(30), noc2(30)
      character*400 msg
      character*13 filnam
      logical around

      do i=1, iatomz
         do ish=1,10
            do ii=1,8
               cphilsh(ii,ish,i)=0.d0
            enddo
         enddo
      enddo

      if(iatomz.le.9) then
         write(namaa(1,0),'(i1.1,"+")')iatomz
      else if(iatomz.le.99) then
         write(namaa(1,0),'(i2.2,"+")')iatomz
      else
         write(namaa(1,0),'(i3.3,"+")')iatomz
      endif

c     bare nuclei
      ntot=1
      nlev(0)=1
      nlbn(0)=1
      elev(ntot)=0.d0
      glev(ntot)=1.d0
      qlev(ntot)=atomz+1.
      isolev(ntot)=0
      ionstg(ntot)=iatomz
      levl(ntot)=1
      indlev(1,0)=ntot       
      kname(ntot)=namaa(1,0)
      mobtot(ntot)=0
      lmatch(ntot)=1
      do ip=1,nprim+1
         noctot(1,ip)=0
      enddo


c     write out energy levels and transitions

      write(filnam,'("atomic.inp.",i2.2)')iatomz
      inquire(file=filnam,exist=around)
      if(.not.around)then
         write(*,'(" Please make ", a13," available")')filnam
         stop
      endif

      open(25, file=filnam)

      jlev=1
      indinp(1,0)=1
      do iso=1, iatomz
         read(25,'(9x,i3,2i5,1x,f16.6)')
     &        i,nlev(i),nlbn(i),eipz(i)
         jl=0
         jb=0
         do lv=1, nlev(i)
            read(25,'(3i5,3x,a8,1x,f16.6,1p,4e12.5,1x,2i4,1x,30i2)')
     &           i,j,ilev,kname(ilev),
     &           elev(ilev),glev(ilev),qlev(ilev),evtoip(ilev),gamx,
     &           ncont(j,i),mobtot(ilev),(noctot(ilev,k),k=1,10)
            indinp(j,i)=0
            if(mobtot(ilev).le.nxprim)then
               jl=jl+1
               jlev=jlev+1
               indinp(j,i)=jl
               indlev(jl,i)=jlev
               namaa(jl,i)=kname(ilev)
               evaa(jl,i)=elev(ilev)
               gaa(jl,i)=glev(ilev)
               sqq(jl,i)=qlev(ilev)
               gamrat(jl,i)=gamx
               morb(jl,i)=mobtot(ilev)
               lvcont(jlev)=indlev(1,i-1)+ncont(j,i)
               ionstg(jlev)=iatomz-i
               isolev(jlev)=i
               levl(jlev)=jl
               lmatch(jlev)=jl
               ncont(jl,i)=indinp(ncont(j,i),i-1)
               do ip=1, nprim
                  nocel(jl,i,ip)=noctot(ilev,ip)
               enddo
               if(ncont(jl,i).eq.1)jb=jb+1
c 2014-05-13
               if(lv.eq.1)then
                  ngr(i)=mobtot(ilev)
                  spnmax(i)=noctot(ilev,mobtot(ilev))
               endif

            endif
            naatot(i)=jl
            naabn(i)=jb
            naaut(i)=jl-jb
         enddo
       enddo


c     review of continuum state
       do iz=1, iatomz
          do i=1, nlev(iz)
             ncontlv=ncont(i,iz)
             do ip=1,nprim
                noc1(ip)=nocel(i,iz,ip)
             enddo
             do j=1, nlev(iz-1)
                do ip=1,nprim
                   noc2(ip)=nocel(j,iz-1,ip)
                enddo
                call getlnpqs(nprim,noc1,noc2,n1,n2)
                if(n1.gt.0.and.n2.eq.0) then
                   if(morb(i,iz).eq.n1.and.ncontlv.ne.j)then
                      ncont(i,iz)=j
                      print *,'change the contiuum state'
                  endif
               endif
            enddo
         enddo
      enddo
      

      ntran=0

      i=0
      m1=0
      m2=0
      m3=0
      read(25,'(10x,i7)') ntran
      
      do it=1, ntran
         read(25,'(a)',end=10)msg
         if(msg(1:1).eq.'b')then
            read(msg,'(2x,2(1x,i2,1x,i4),2i4,i8,1p,8e13.5)')
     &           iso,ilo,is2,iup,n1,n2,itran,
     &           faax, oratex, cfax, (cxin(ix),ix=1,5)
            if(indinp(ilo,iso)*indinp(iup,is2).gt.0)then
               i=i+1
               m1=m1+1
               ltype(i)=1
               faa(indinp(ilo,iso),indinp(iup,is2),iso)=faax
               orateaa(indinp(iup,is2),indinp(ilo,iso),iso)=oratex
               cfaa(indinp(ilo,iso),indinp(iup,is2),iso)=cfax
               if(itran.gt.0)then
                  do ix=1,5
                     cxjj(itran,ix,iso)=cxin(ix)
                  enddo
               endif
               nloex(indinp(ilo,iso),indinp(iup,is2),iso)=n1
               nhiex(indinp(ilo,iso),indinp(iup,is2),iso)=n2
            endif
         else if(msg(1:1).eq.'f')then
            read(msg,'(2x,2(1x,i2,1x,i4),1x,i3,1x,1p,9e13.5)')
     &           iso,ilo,is2,iup,ish,psh,(cxin(ii),ii=1,8)
            if(indinp(ilo,iso)*indinp(iup,is2).gt.0)then
               i=i+1
               m2=m2+1
               ltype(i)=2
               istrn(indinp(ilo,iso),indinp(iup,is2),iso)=ish
               petrn(indinp(ilo,iso),indinp(iup,is2),iso)=psh
               if(ish.gt.0.and.ish.le.10)then
                  do ii=1, 8
                     cphilsh(ii,ish,iso)=cxin(ii)
                  enddo
               endif
               npqi(indinp(ilo,iso),indinp(iup,is2),iso)=int(psh)
               nish(indinp(ilo,iso),indinp(iup,is2),iso)=ish
            endif
         else if(msg(1:1).eq.'a')then
            read(msg,'(2x,2(1x,i2,1x,i4),1x,1p,e13.5)')
     &           iso,ilo,is2,iup,autxx
            if(indinp(ilo,iso)*indinp(iup,is2).gt.0)then
               i=i+1
               m3=m3+1
               ltype(i)=3
               autoaa(indinp(ilo,iso),indinp(iup,is2),iso)=autxx
            endif
         else
            print *, 'error in readatom'
         endif
      enddo
      
 10   close(25)

      print *, '# of levels and transitions read from atomic.inp  ', 
     &     jlev, ntran

      ntran=i

      return
      end
