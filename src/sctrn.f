c          *******  ***  *******  *        ******
c          *         *   *        *         *    *
c          *         *   *        *         *    *
c          ****      *   ****     *         *    *
c          *         *   *        *         *    *
c          *         *   *        *         *    *
c          *        ***  *******  *******  ******
c
      function field(energy)
c
c___calculate the radiation field due to input field at energy
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'timestuf'
      include 'runstuf'
c
      field=0.d0
      if (energy .gt. xev(ngrid)) then
         if(xbarj.gt.0.d0) then
            y = energy/trev
            if (y .lt. 100.) then
               field = barjend*energy**3*exp(-y)*dilution
            endif
         endif
      else if (energy .lt. xev(1)) then
         if(xbarj.gt.0.d0) then
            field = barjnow(1)*energy/xev(1)*dilution
         endif
      else
         do i=1,ngrid-1
            if(energy.ge.xev(i).and.energy.lt.xev(i+1))then
               iplace=i
               bit = (energy-xev(iplace))/
     &              (xev(iplace+1)-xev(iplace))
               field = (barjnow(iplace)*(1.-bit)+
     1              barjnow(iplace+1)*bit)
     1              *dilution
               return
            endif
         enddo
      endif
c
      return
      end
c
c      *****    *******  ***  *******  *        *****
c      *    *   *         *   *        *        *    *
c      *    *   *         *   *        *        *    *
c      *****    ****      *   ****     *        *    *
c      *    *   *         *   *        *        *    *
c      *    *   *         *   *        *        *    *
c      *****    *        ***  *******  *******  *****
c
      function bbfield(y)
c
c___calculate the radiation field of a Plankian at energy = y*trev
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'timestuf'
      include 'runstuf'
c
      energy = trev*y
      if (y .lt. 0.001) then
         bbfield = 2.084e-04*trev**3*y**2*dilution
      else if (y .gt. 10.0) then
         bbfield = 2.084e-04*(energy)**3*exp(-y)*dilution
      else
         bbfield = 2.084e-04*(energy)**3/(exp(y)-1.0)*dilution
      endif

      return
      end
c
c
c
c          *******   *****    ******
c          *        *     *  *
c          *        *     *  *
c          ****     *     *   *****
c          *        *     *        *
c          *        *     *        *
c          *         *****   ******
c
c
c
c
      function fos(i,j)
c
c___ determine absorption oscillator strength
c___ if lower state is greater than 20 or upper
c___ state is greater than 25 use asymptotic form.
c
      implicit real*8 (a-h,o-z)
      common /fvalue/ fosc(290)
c
      data (fosc(ii),ii=1,60)/
     1.4162e+00,.7910e-01,.2899e-01,.1394e-01,.7799e-02,
     1.4814e-02,.3183e-02,.2216e-02,.1605e-02,.1201e-02,
     1.9214e-03,.7227e-03,.5774e-03,.4686e-03,.3856e-03,
     1.3211e-03,.2702e-03,.2296e-03,.1967e-03,.1690e-03,
     1.1470e-03,.1290e-03,.1140e-03,.1010e-03,.6407e+00,
     1.1193e+00,.4467e-01,.2209e-01,.1270e-01,.8036e-02,
     1.5429e-02,.3851e-02,.2835e-02,.2151e-02,.1672e-02,
     1.1326e-02,.1070e-02,.8764e-03,.7270e-03,.6099e-03,
     1.5167e-03,.4416e-03,.3810e-03,.3300e-03,.2880e-03,
     1.2534e-03,.2238e-03,.8421e+00,.1506e+00,.5584e-01,
     1.2768e-01,.1604e-01,.1023e-01,.6980e-02,.4996e-02,
     1.3711e-02,.2839e-02,.2224e-02,.1776e-02,.1443e-02/
      data (fosc(ii),ii= 61,120)/
     1.1188e-02,.9092e-03,.8361e-03,.7118e-03,.6020e-03,
     1.5100e-03,.4340e-03,.3710e-03,.3160e-03,.1038e+01,
     1.1793e+00,.6549e-01,.3230e-01,.1870e-01,.1196e-01,
     1.8187e-02,.5886e-02,.4393e-02,.3375e-02,.2656e-02,
     1.2131e-02,.1739e-02,.1439e-02,.1204e-02,.1019e-02,
     1.8410e-03,.7060e-03,.6000e-03,.5180e-03,.4410e-03,
     1.1231e+01,.2069e+00,.7448e-01,.3645e-01,.2104e-01,
     1.1344e-01,.9209e-02,.6631e-02,.4959e-02,.3821e-02,
     1.3014e-02,.2425e-02,.1984e-02,.1646e-02,.1382e-02,
     1.1170e-02,.1010e-02,.8570e-03,.8600e-03,.6480e-03,
     1.1424e+01,.2340e+00,.8315e-01,.4038e-01,.2320e-01,
     1.1479e-01,.1012e-01,.7289e-02,.5455e-02,.4207e-02/
      data (fosc(ii),ii=121,180)/
     1.3324e-02,.2679e-02,.2196e-02,.1825e-02,.1540e-02,
     1.1310e-02,.1120e-02,.9650e-03,.8400e-03,.1616e+01,
     1.2609e+00,.9163e-01,.4416e-01,.2525e-01,.1605e-01,
     1.1097e-01,.7891e-02,.5905e-02,.4556e-02,.3602e-02,
     1.2905e-02,.2383e-02,.1980e-02,.1660e-02,.1430e-02,
     1.1240e-02,.1080e-02,.1807e+01,.2876e+00,.1000e+00,
     1.4787e-01,.2724e-01,.1726e-01,.1177e-01,.8456e-02,
     1.6323e-02,.4877e-02,.3856e-02,.3122e-02,.2590e-02,
     1.2140e-02,.1790e-02,.1500e-02,.1280e-02,.1999e+01,
     1.3143e+00,.1083e+00,.5152e-01,.2918e-01,.1843e-01,
     1.1254e-01,.8995e-02,.6719e-02,.5180e-02,.4094e-02,
     1.3300e-02,.2700e-02,.2270e-02,.1920e-02,.1660e-02/
      data (fosc(ii),ii=181,240)/
     1.2190e+01,.3408e+00,.1166e+00,.5513e-01,.3109e-01,
     1.1958e-01,.1328e-01,.9515e-02,.7099e-02,.5468e-02,
     1.4270e-02,.3400e-02,.2810e-02,.2360e-02,.2010e-02,
     1.2381e+01,.3673e+00,.1248e+00,.5872e-01,.3298e-01,
     1.2070e-01,.1402e-01,.1002e-01,.7468e-02,.5520e-02,
     1.4400e-02,.3530e-02,.2900e-02,.2400e-02,.2572e+01,
     1.3938e+00,.1330e+00,.6228e-01,.3486e-01,.2182e-01,
     1.1474e-01,.1052e-01,.7390e-02,.5410e-02,.4020e-02,
     1.3040e-02,.2640e-02,.2763e+01,.4202e+00,.1412e+00,
     1.6584e-01,.3672e-01,.2292e-01,.1545e-01,.1060e-01,
     1.7510e-02,.5480e-02,.4060e-02,.3100e-02,.2954e+01,
     1.4467e+00,.1494e+00,.6938e-01,.3858e-01,.2402e-01/
      data (fosc(ii),ii=241,290)/
     1.1600e-01,.1160e-01,.8760e-02,.6700e-02,.5280e-02,
     1.3145e+01,.4731e+00,.1575e+00,.7292e-01,.4043e-01,
     1.2860e-01,.2130e-01,.1760e-01,.1230e-01,.1010e-01,
     1.3336e+01,.4995e+00,.1657e+00,.7644e-01,.4360e-01,
     1.3290e-01,.2610e-01,.2180e-01,.1960e-01,.3527e+01,
     1.5259e+00,.1738e+00,.8600e-01,.6650e-01,.5050e-01,
     1.3020e-01,.2990e-01,.3718e+01,.5523e+00,.1830e+00,
     1.9320e-01,.7930e-01,.3930e-01,.3190e-01,.3909e+01,
     1.5760e+00,.1910e+00,.8610e-01,.4320e-01,.3250e-01,
     1.4150e+01,.6000e+00,.1990e+00,.8900e-01,.4710e-01/
c
      if (i.gt.20) go to 20
      if (j.gt.25) go to 20
c
      iplace = 300-((25-i)*(26-i))/2+j-i
c
      fos=fosc(iplace)
      return
c
   20 ei = i
      ej = j
c
      fos=1.9603*(ej/(ej*ej-ei*ei))**3*ei
      return
      end
 
c
c           *****     * *      * *     *******  ***  *        *
c          *     *   *   *    *   *    *         *   *        *
c          *        *     *  *     *   *         *   *        *
c          *        *******  *******   ****      *   *        *
c          *        *     *  *     *   *         *   *        *
c          *     *  *     *  *     *   *         *   *        *
c           *****   *     *  *     *   *        ***  *******  *******
c
c
      function caa(iir,ijr,iz)
c
c___ this calculates the collisional excitation or de-excitation
c___ from level i to level j of an ion (with iz electrons)
c___ in the average-atom model.
c
c    collisional rate option
c    isamp(1)=1 Van Regemorter

      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'popstuf'
      parameter(mpq=10)

      caa = 0.d0

      i=min(iir,ijr)
      j=max(ijr,iir)
      ilv=indlev(i,iz)
      jlv=indlev(j,iz)

      ilv0=indlev(ncont(i,iz),iz-1)
c      eii = elev(ilv0)-elev(ilv)+debye(iatomz-iz+1)
c     correction on 2/28/2011
      eii = elev(ilv0)-elev(ilv)
      if(eii.le.0.d0) return

      ir=morb(i,iz)
      zeff = sqq(i,iz)
      fij = cfaa(i,j,iz)

      if(fij.eq.0.d0) return
      eji = elev(jlv)-elev(ilv)
      if(eji.le.0.d0) return

      ebt = eji/tev
      n1=nloex(i,j,iz)          ! initial shell of lower level
      n2=nhiex(i,j,iz)          ! final shell
      itr=mpq*(n1-1)-((n1-1)*n1)/2+n2-n1

      gaunt = 0.15 + 0.28*xexpe1(ebt)/ebt
      if(isamp(1).eq.0.and.n1*n2.ne.0)then
         if(n1.le.9.and.n2.le.10)gnew = gcexjj(ebt,itr,iz)
         if(gnew.gt.0.d0)gaunt=gnew
      endif

      if(iir.eq.i) then
         caa = 1.581d-5/eji*fij/sqrt(tev)*exp(-ebt)*gaunt
      else if(iir.eq.j) then
         caa = 1.581d-5/eji*fij/sqrt(tev)*glev(ilv)/glev(jlv)*gaunt
      endif


      return
      end

      function gcexjj(tx,i,iso)
      implicit double precision(a-h,o-z)
      common /jjcx/ cxjj(45,5,100) 

      c0=cxjj(i,1,iso)
      c1=cxjj(i,2,iso)
      c2=cxjj(i,3,iso)
      c3=cxjj(i,4,iso)
      a0=cxjj(i,5,iso)
      a1=a0+1.d0

      gcexjj=(c0*exp1x(tx) + (c1 + c3*tx/a1)+
     &     (c2-c3*tx)*tx*exp1x(tx*a1))
      return
      end
c
      function gcstjj(y,i,iso)
      implicit double precision(a-h,o-z)
      common /jjcx/ cxjj(45,5,100) 

      c0=cxjj(i,1,iso)
      c1=cxjj(i,2,iso)
      c2=cxjj(i,3,iso)
      c3=cxjj(i,4,iso)
      a0=cxjj(i,5,iso)

      gcstjj=c0*log(y)+c1+c2/(y+a0)+c3/(y+a0)**2
      return
      end


c           *****     ***    ******   *******    ***      ***   
c          *     *   *   *    *    *  *         *   *    *   *  
c          *        *     *   *    *  *        *     *  *     * 
c          *        *******   *****   ****     *******  ******* 
c          *        *     *   *       *        *     *  *     * 
c          *     *  *     *   *       *        *     *  *     * 
c           *****   *     *   *       *******  *     *  *     * 
c
c
c
c
      function capexx(i,j,iz)
c
c---- rate of electron capture from the state (j, iz-1)
c----     into autionizing state (i,iz)
c
      implicit real*8 (a-h,o-z)
      include 'mainstuf'
      include 'popstuf' 
      capexx=0.d0
      ntoti=indlev(i,iz)
      ntotj=indlev(j,iz-1)
      eii=elev(ntoti)-elev(ntotj)
      if(eii.ge.0.d0)then
         capexx = 1.656415d-22*exp(-eii/tev)/tev**1.5*
     &        glev(ntoti)/glev(ntotj)*autoaa(i,j,iz)
      endif

      return
      end

c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c*                                                                     *
c*       DDDD     IIII   EEEEEE  RRRRR   EEEEEE    CC                  *
c*       DDDDD     II    EEEEEE  RR  RR  EEEEEE   CCCC                 *
c*       DD  DD    II    EE      RR  RR  EE      CC  CC                *
c*       DD  DD    II    EEEEEE  RRRRR   EEEEEE  CC                    *
c*       DD  DD    II    EEEEEE  RRRR    EEEEEE  CC                    *
c*       DD  DD    II    EE      RR RR   EE      CC  CC                *
c*       DDDDD    IIII   EEEEEE  RR  RR  EEEEEE   CCCC                 *
c*       DDDD     IIII   EEEEEE  RR  RR  EEEEEE    CC                  *
c*                                                                     *
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      function dierec (iz)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      calculate dielectronic recombination rate to the ion iz from iz-1
c      in its ground state
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit real*8 (a-h,o-z)

      include 'mainstuf'

      external fos


      xne = denow
      nzn = atomz

      ze = real(iatomz-iz+1)    ! ionization state of the recombining ion iz-1
      ib = iz-1                 ! number of bound electrons of ion iz-1
      diet = 0.d0
      if(ib.gt.60) go to 1


      aa=1.0+0.015*ze**3/(ze+1.0)**2
      bb=sqrt(ze)*(ze+1.0)**2.5/sqrt(ze**2+13.4)

c original dca

      xntt=1.5084d+17*ze**6*sqrt(tev)/xne
      xntt=xntt**0.142857   
      ntt=xntt+1

c Griems'
c make executable with fly_n.out
c      xntt=1.9d+18*ze**6*sqrt(tev)/xne
c      xntt=xntt**(2./17.)
c      ntt=xntt+1

cccccccccccccccccccccccccc
c
c    n-m transition
c
cccccccccccccccccccccccccc

      nlte=amin0(ntt,nxprim)
      ddt=((ze+1.0)*xntt)**2*0.0015
      ddt=ddt/(1.0+ddt)
      nmax=ngr(ib)  ! caution for the index of ngr!!
      nmax1=nmax+1
      diemn=0.0

      if(nlte.lt.nmax1) go to 101


      do 100 j=nmax1,nlte

!         emn=evaa(j,ib)       
         emn=evaa(j-nmax+1,ib)    ! energy difference between j and g.s.
         emn=emn/(ze**2*13.605) ! 1/v_i^2 - 1/v_j^2
         y=(ze+1.0)*emn         ! (z+1)/(1/v_i^2 - 1/v_j^2)
         ebt=(ze+1.0)**2*emn/(7.35d-2*aa*tev)
         aay=0.5*sqrt(y)/(1.0+0.210*y+0.030*y**2)
         diemn=diemn+fos(nmax,j)*aay*exp(-ebt)
 100  continue

      diemn=diemn*ddt

 101  continue

cccccccccccccccccccccccc
c
c   n-n transition
c
cccccccccccccccccccccccc

      zeet=dmax1(ze,1.0d0)
      enn=ennn(ib)/(zeet**2*13.605)
      y=(ze+1.0)*enn
      ebt=(ze+1.0)**2*enn/(7.35d-2*aa*tev)
      aay=sqrt(y)/(1.0+0.105*y+0.015*y**2)
      ddt=xntt/200.0
      ddt=ddt/(1.0+ddt)
      dienn=fnn(ib)*aay*exp(-ebt)
c     
      dienn=dienn*ddt
c
c   sum of dn=0 and dn=!0 contributions
c
      if(ibdr.eq.1) then
         diet=2.4E-9/tev**1.5*bb*(diemn+dienn)*spnmax(ib)
      else if(ibdr.eq.2)then    
         diet=2.4E-9/tev**1.5*bb*(diemn)*spnmax(ib)
      else if(ibdr.eq.3) then
         diet=2.4E-9/tev**1.5*bb*(dienn)*spnmax(ib)
      else if(ibdr.eq.0) then
         diet=0.d0
      endif

    1 continue

c fbdr is a  DR multiplier
      dierec = diet*fbdr     
      return
      end

c
      function phicx(xx)
      implicit real*8 (a-h,o-z)      
      common /phistuff/ cp0,cp1,cp2,cp3,xn,cep,cem
      phicx=exp(cp0+xx*(cp1+xx*(cp2+xx*cp3)))
      return
      end

c----------------------------------------------
      real*8 function phrec(ir,jr,iz)
      implicit real*8 (a-h,o-z)

      include 'mainstuf'
      include 'runstuf'
      include 'timestuf'
      include 'popstuf'

      common /phistuff/ cp0,cp1,cp2,cp3,xn,cep,cem
      data fourpi /12.5664d0/
c
      dimension xt(mgas),ft(mgas),pt(mgas)
c
      estm=0.d0
      phrec=0.d0


      ntoti=indlev(ir,iz)
      ntotj=indlev(jr,iz-1)
      eip= elev(ntotj)-elev(ntoti)
      if(eip.le.0.d0) return

      debaln=1.656415d-22/tev**1.5*glev(ntoti)/glev(ntotj)

      ish=nish(ir,jr,iz)
      if(ish.ne.0)then
         cp0=cphilsh(1,ish,iz)
         cp1=cphilsh(2,ish,iz)
         cp2=cphilsh(3,ish,iz)
         cp3=cphilsh(4,ish,iz)
         xn =cphilsh(5,ish,iz)
         cep=cphilsh(6,ish,iz) 
         cem=cphilsh(7,ish,iz) 
         emax=cphilsh(8,ish,iz) 
      endif
      if(xn.le.1.d-30.or.ish.eq.0) then
         if(radflag.ne.'off')then
            phrec=alfxx(ir,jr,iz)  ! estm is computed within
            ealf=estm
            phrec=phrec+bstmxx(ir,jr,iz)  ! estm is computed within
            estm=ealf+estm
         else
            phrec=alfxx(ir,jr,iz)  ! estm is computed within
         endif
         return
      endif

c spontaneous recombination

      emin=eip !max(eip,cep)
      xm=(emax+emin)/2.
      xr=(emax-emin)/2.
      do i=1,ngas
         xt(i)=xr*xgas(i)+xm
      enddo
      call phifun(xt,ft,ngas,eip) ! cross-section in cm2

      valrec=0.
      do i=1,ngas
         if((xt(i)-emin)/tev.lt.200.)then
            prec=1.636e9*xt(i)**2*exp(-(xt(i)-emin)/tev) !(8Piv^2/c^2)=1.737e9
            valrec=valrec+wgas(i)*ft(i)*prec  
            estm=estm+wgas(i)*ft(i)*prec*xt(i)  
         endif
      enddo
      phrec=valrec*xr*2.4179e14*debaln
      estm=estm*xr*2.4179e14*debaln

      if(radflag.eq.'off') return

c stimulated recombination

      if(radflag.eq.'trfile') then
         emin=max(eip,xev(1)) !emin=max(cep,xev(1))
         emax=min(emax,xev(ngrid))
         xm=(emax+emin)/2.
         xr=(emax-emin)/2.
         if(xr.le.0.d0)return
      endif
      do i=1,ngas
         xt(i)=xr*xgas(i)+xm
         pt(i)=0.
      enddo
      call phifun(xt,ft,ngas,eip) ! cross-section in cm2
      call phtflux(xt,pt,ngas) ! 4Pi.Jv/hv in 1/cm2/s/Hz


c for stimulated recombination
c choose edge=min(cep,eip)    
      edge=min(cep,eip)
      valrec=0.d0
      xxx=0.d0
      do i=1,ngas
         if((xt(i)-edge)/tev.lt.200.)then
            prec=pt(i)*exp(-(xt(i)-edge)/tev)
            valrec=valrec+wgas(i)*ft(i)*prec
            xxx=xxx+xt(i)*wgas(i)*ft(i)*prec
         endif
      enddo

c total recombination

      phrec=phrec + valrec*xr*2.4179e14*debaln
      estm=estm + xxx*xr*2.4179e14*debaln

      return
      end
c----------------------------------------------
      real*8 function  phion(ir,jr,iz)
      implicit real*8 (a-h,o-z)

      include 'mainstuf'
      include 'runstuf'
      include 'timestuf'
      include 'popstuf'

      common /phistuff/ cp0,cp1,cp2,cp3,xn,cep,cem
      data fourpi /12.5664d0/
c
      dimension xt(mgas),ft(mgas),pt(mgas)
c
      econt=0.d0
      phion=0.d0

      ntoti=indlev(ir,iz)
      ntotj=indlev(jr,iz-1)
      eip= elev(ntotj)-elev(ntoti)
      if(eip.le.0.d0) return

      debaln=1.656415d-22/tev**1.5*exp(eip/tev)*glev(ntoti)/glev(ntotj)

      ish=nish(ir,jr,iz)
      if(ish.ne.0)then
         cp0=cphilsh(1,ish,iz)
         cp1=cphilsh(2,ish,iz)
         cp2=cphilsh(3,ish,iz)
         cp3=cphilsh(4,ish,iz)
         xn =cphilsh(5,ish,iz)
         cep=cphilsh(6,ish,iz) 
         cem=cphilsh(7,ish,iz) 
         emax=cphilsh(8,ish,iz) 
      endif
      if(xn.le.1.d-30.or.ish.eq.0) then
         phion=bcontxx(ir,jr,iz)
         return
      endif

      if(radflag.ne.'trfile') then
         if(trev.eq.0.d0) return
         emin=eip !max(eip,cep)
         xm=(emax+emin)/2.
         xr=(emax-emin)/2.
      else
         emin=max(eip,xev(1)) !emin=max(cep,xev(1))
         emax=min(emax,xev(ngrid))
         xm=(emax+emin)/2.
         xr=(emax-emin)/2.
c         if(xr.le.0.d0)phion=0.d0
         if(xr.le.0.d0) then
            phion=0.d0
            return
         endif
      endif
      do 10 i=1,ngas
         xt(i)=xr*xgas(i)+xm
         pt(i)=0.
 10   continue
c      call phifun(xt,ft,ngas)   ! cross-section in cm2
      call phifun(xt,ft,ngas,eip)   ! cross-section in cm2
      call phtflux(xt,pt,ngas) ! 4Pi.Jv/hv in 1/cm2/s/Hz

      valion=0.d0
      do 20 i=1,ngas
         valion=valion+wgas(i)*ft(i)*pt(i)
         econt=econt+xt(i)*wgas(i)*ft(i)*pt(i)
 20   continue

      phion=xr*valion*2.4179e14 ! in 1/s
      econt=xr*econt*2.4179e14

      return
      end


c--------------------------------
      subroutine phifun(xarr,farr,narr,eip)
      implicit real*8 (a-h,o-z)
      common /phistuff/ cp0,cp1,cp2,cp3,xn,cep,cem
      dimension xarr(narr),farr(narr)
      enemi=cep !min(eip,cep) 
      do 10 i=1,narr
         xx=max(log(xarr(i)/enemi),0.d0)
         yy=cp0+xx*(cp1+xx*(cp2+xx*cp3))
         farr(i)=1.0e-18*xn*exp(yy)*(13.606/enemi)
 10   continue
      return
      end
c------------------------------------------
      subroutine phtflux(xarr,farr,narr)
      implicit real*8 (a-h,o-z)
      include 'runstuf'
      common /tradi8/ tr,trev,trfixed,dilution,barjnow(100),
     1                barjend,trmin,xbarj,tmin
      dimension xarr(narr),farr(narr)
      if(radflag.ne.'trfile') then
         if(trev.gt.0.d0) then
            do i=1,narr
c            farr(i)=bbfield(xarr(i)/trev)*4*pi/xarr(i)/1.60219e-12
               farr(i)=bbfield(xarr(i)/trev)/xarr(i)*7.84325e12
            enddo
         else
            do i=1,narr
               farr(i)=0.d0
            enddo
         endif
      else 
         do i=1,narr
            farr(i)=field(xarr(i))/xarr(i)*7.84325e12
         enddo
      endif
      return
      end
c
c--------------------------------------------
c
c          *******   *
c          *        **
c          *         *
c          ****      *
c          *         *
c          *         *
c          *******  ***
c
c
c
c
      function e1(xt,expx)
 
      implicit real*8 (a-h,o-z)
      common /exp1dat/ e1data(334),dexp1,deoff
      common /exp1datI/ nexp1
 
      data nexp1/334/, dexp1/.075/, deoff/.025/
      data (e1data(i),i=1,334)/
     & 3.137E+00,1.823E+00, 1.333E+00, 1.044E+00, 8.471E-01, 7.024E-01,
     & 5.913E-01,5.034E-01, 4.323E-01, 3.738E-01, 3.250E-01, 2.840E-01,
     & 2.492E-01,2.194E-01, 1.937E-01, 1.716E-01, 1.523E-01, 1.355E-01,
     & 1.207E-01,1.078E-01, 9.638E-02, 8.631E-02, 7.740E-02, 6.949E-02,
     & 6.246E-02,5.620E-02, 5.062E-02, 4.564E-02, 4.118E-02, 3.719E-02,
     & 3.361E-02,3.040E-02, 2.751E-02, 2.491E-02, 2.258E-02, 2.047E-02,
     & 1.857E-02,1.686E-02, 1.531E-02, 1.391E-02, 1.264E-02, 1.149E-02,
     & 1.046E-02,9.517E-03, 8.664E-03, 7.891E-03, 7.189E-03, 6.552E-03,
     & 5.974E-03,5.448E-03, 4.970E-03, 4.535E-03, 4.139E-03, 3.779E-03,
     & 3.452E-03,3.153E-03, 2.881E-03, 2.633E-03, 2.407E-03, 2.201E-03,
     & 2.013E-03,1.841E-03, 1.684E-03, 1.541E-03, 1.411E-03, 1.291E-03,
     & 1.182E-03,1.083E-03, 9.919E-04, 9.086E-04, 8.325E-04, 7.629E-04,
     & 6.992E-04,6.409E-04, 5.876E-04, 5.388E-04, 4.941E-04, 4.532E-04,
     & 4.157E-04,3.814E-04, 3.499E-04, 3.211E-04, 2.947E-04, 2.705E-04,
     & 2.483E-04,2.279E-04, 2.093E-04, 1.922E-04, 1.765E-04, 1.621E-04,
     & 1.489E-04,1.368E-04, 1.257E-04, 1.155E-04, 1.061E-04, 9.752E-05,
     & 8.963E-05,8.239E-05, 7.573E-05, 6.962E-05, 6.401E-05, 5.886E-05,
     & 5.412E-05,4.977E-05, 4.578E-05, 4.210E-05, 3.873E-05, 3.563E-05,
     & 3.278E-05,3.015E-05, 2.775E-05, 2.553E-05, 2.349E-05, 2.162E-05,
     & 1.990E-05,1.832E-05, 1.686E-05, 1.552E-05, 1.429E-05, 1.315E-05,
     & 1.211E-05,1.115E-05, 1.027E-05, 9.454E-06, 8.706E-06, 8.018E-06,
     & 7.385E-06,6.802E-06, 6.265E-06, 5.771E-06, 5.316E-06, 4.897E-06,
     & 4.512E-06,4.157E-06, 3.830E-06, 3.529E-06, 3.252E-06, 2.997E-06,
     & 2.762E-06,2.545E-06, 2.346E-06, 2.162E-06, 1.993E-06, 1.837E-06,
     & 1.693E-06,1.561E-06, 1.439E-06, 1.326E-06, 1.223E-06, 1.127E-06,
     & 1.040E-06,9.585E-07, 8.838E-07, 8.150E-07, 7.515E-07, 6.931E-07,
     & 6.392E-07,5.895E-07, 5.436E-07, 5.014E-07, 4.625E-07, 4.266E-07,
     & 3.935E-07,3.630E-07, 3.348E-07, 3.089E-07, 2.849E-07, 2.629E-07,
     & 2.425E-07,2.238E-07, 2.065E-07, 1.905E-07, 1.758E-07, 1.622E-07,
     & 1.497E-07,1.381E-07, 1.274E-07, 1.176E-07, 1.085E-07, 1.002E-07,
     & 9.245E-08,8.532E-08, 7.875E-08, 7.268E-08, 6.709E-08, 6.193E-08,
     & 5.716E-08,5.276E-08, 4.871E-08, 4.496E-08, 4.151E-08, 3.832E-08,
     & 3.538E-08,3.266E-08, 3.015E-08, 2.784E-08, 2.570E-08, 2.373E-08,
     & 2.191E-08,2.023E-08, 1.868E-08, 1.725E-08, 1.593E-08, 1.471E-08,
     & 1.358E-08,1.255E-08, 1.159E-08, 1.070E-08, 9.881E-09, 9.126E-09,
     & 8.428E-09,7.784E-09, 7.190E-09, 6.640E-09, 6.133E-09, 5.665E-09,
     & 5.233E-09,4.834E-09, 4.465E-09, 4.124E-09, 3.810E-09, 3.519E-09,
     & 3.251E-09,3.003E-09, 2.775E-09, 2.563E-09, 2.368E-09, 2.188E-09,
     & 2.021E-09,1.867E-09, 1.725E-09, 1.594E-09, 1.473E-09, 1.361E-09,
     & 1.257E-09,1.162E-09, 1.074E-09, 9.920E-10, 9.167E-10, 8.471E-10,
     & 7.827E-10,7.233E-10, 6.684E-10, 6.177E-10, 5.708E-10, 5.275E-10,
     & 4.875E-10,4.505E-10, 4.164E-10, 3.848E-10, 3.556E-10, 3.287E-10,
     & 3.038E-10,2.808E-10, 2.595E-10, 2.399E-10, 2.217E-10, 2.049E-10,
     & 1.894E-10,1.751E-10, 1.618E-10, 1.496E-10, 1.383E-10, 1.278E-10,
     & 1.182E-10,1.092E-10, 1.010E-10, 9.334E-11, 8.628E-11, 7.976E-11,
     & 7.374E-11,6.817E-11, 6.302E-11, 5.826E-11, 5.386E-11, 4.980E-11,
     & 4.604E-11,4.257E-11, 3.935E-11, 3.639E-11, 3.364E-11, 3.110E-11,
     & 2.876E-11,2.659E-11, 2.459E-11, 2.273E-11, 2.102E-11, 1.944E-11,
     & 1.797E-11,1.662E-11, 1.537E-11, 1.421E-11, 1.314E-11, 1.215E-11,
     & 1.123E-11,1.039E-11, 9.607E-12, 8.884E-12, 8.216E-12, 7.598E-12,
     & 7.026E-12,6.498E-12, 6.009E-12, 5.557E-12, 5.139E-12, 4.753E-12,
     & 4.396E-12,4.065E-12, 3.760E-12, 3.477E-12, 3.216E-12, 2.974E-12,
     & 2.751E-12,2.544E-12, 2.353E-12, 2.177E-12, 2.013E-12, 1.862E-12,
     & 1.722E-12,1.593E-12, 1.473E-12, 1.363E-12, 1.261E-12, 1.166E-12,
     & 1.079E-12,9.977E-13, 9.229E-13, 8.537E-13, 7.897E-13, 7.305E-13,
     & 6.757E-13,6.251E-13, 5.782E-13, 5.349E-13/
 
      if (xt .lt. .1) then
        e1 = -log(xt)-0.57721+xt
      else if (xt .gt. 25.) then
         e1 = expx/xt*(1.-0.95/xt)
      else
         r = (xt-deoff)/dexp1+1.0
         intr = int(r)
         frac = r-real(intr)
         e1 = e1data(intr)*(1.-frac)+e1data(intr+1)*frac
      endif
 
      return
      end
 
      function xexpe1(x)
 
      implicit real*8 (a-h,o-z)
      common /xexe1/ exe1data(334),dexe1,dexe1off
      common /xexe1I/ nexe1
      data nexe1/334/, dexe1/.075/, dexe1off/.025/
      data (exe1data(i),i=1,334)/
     & 8.040E-02,2.015E-01, 2.780E-01, 3.352E-01, 3.810E-01, 4.191E-01,
     & 4.516E-01,4.799E-01, 5.047E-01, 5.269E-01, 5.468E-01, 5.648E-01,
     & 5.813E-01,5.963E-01, 6.102E-01, 6.231E-01, 6.350E-01, 6.461E-01,
     & 6.565E-01,6.662E-01, 6.754E-01, 6.840E-01, 6.921E-01, 6.998E-01,
     & 7.071E-01,7.140E-01, 7.205E-01, 7.268E-01, 7.327E-01, 7.384E-01,
     & 7.439E-01,7.491E-01, 7.540E-01, 7.588E-01, 7.634E-01, 7.678E-01,
     & 7.720E-01,7.761E-01, 7.800E-01, 7.838E-01, 7.875E-01, 7.910E-01,
     & 7.944E-01,7.977E-01, 8.008E-01, 8.039E-01, 8.069E-01, 8.098E-01,
     & 8.126E-01,8.153E-01, 8.179E-01, 8.205E-01, 8.230E-01, 8.254E-01,
     & 8.277E-01,8.300E-01, 8.322E-01, 8.344E-01, 8.365E-01, 8.385E-01,
     & 8.405E-01,8.425E-01, 8.444E-01, 8.462E-01, 8.481E-01, 8.498E-01,
     & 8.515E-01,8.532E-01, 8.549E-01, 8.565E-01, 8.581E-01, 8.596E-01,
     & 8.611E-01,8.626E-01, 8.640E-01, 8.654E-01, 8.668E-01, 8.681E-01,
     & 8.695E-01,8.708E-01, 8.720E-01, 8.733E-01, 8.745E-01, 8.757E-01,
     & 8.769E-01,8.780E-01, 8.791E-01, 8.802E-01, 8.813E-01, 8.824E-01,
     & 8.835E-01,8.845E-01, 8.855E-01, 8.865E-01, 8.875E-01, 8.884E-01,
     & 8.894E-01,8.903E-01, 8.912E-01, 8.921E-01, 8.930E-01, 8.938E-01,
     & 8.947E-01,8.955E-01, 8.964E-01, 8.972E-01, 8.980E-01, 8.988E-01,
     & 8.995E-01,9.003E-01, 9.010E-01, 9.018E-01, 9.025E-01, 9.032E-01,
     & 9.039E-01,9.046E-01, 9.053E-01, 9.060E-01, 9.067E-01, 9.073E-01,
     & 9.080E-01,9.086E-01, 9.092E-01, 9.099E-01, 9.105E-01, 9.111E-01,
     & 9.117E-01,9.123E-01, 9.128E-01, 9.134E-01, 9.140E-01, 9.145E-01,
     & 9.151E-01,9.156E-01, 9.162E-01, 9.167E-01, 9.172E-01, 9.177E-01,
     & 9.182E-01,9.188E-01, 9.192E-01, 9.197E-01, 9.202E-01, 9.207E-01,
     & 9.212E-01,9.216E-01, 9.221E-01, 9.226E-01, 9.230E-01, 9.235E-01,
     & 9.239E-01,9.243E-01, 9.248E-01, 9.252E-01, 9.256E-01, 9.260E-01,
     & 9.264E-01,9.269E-01, 9.273E-01, 9.277E-01, 9.280E-01, 9.284E-01,
     & 9.288E-01,9.292E-01, 9.296E-01, 9.299E-01, 9.303E-01, 9.307E-01,
     & 9.310E-01,9.314E-01, 9.318E-01, 9.321E-01, 9.325E-01, 9.328E-01,
     & 9.331E-01,9.335E-01, 9.338E-01, 9.341E-01, 9.345E-01, 9.348E-01,
     & 9.351E-01,9.354E-01, 9.357E-01, 9.360E-01, 9.363E-01, 9.367E-01,
     & 9.370E-01,9.373E-01, 9.375E-01, 9.378E-01, 9.381E-01, 9.384E-01,
     & 9.387E-01,9.390E-01, 9.393E-01, 9.395E-01, 9.398E-01, 9.401E-01,
     & 9.404E-01,9.406E-01, 9.409E-01, 9.412E-01, 9.414E-01, 9.417E-01,
     & 9.419E-01,9.422E-01, 9.424E-01, 9.427E-01, 9.429E-01, 9.432E-01,
     & 9.434E-01,9.437E-01, 9.439E-01, 9.441E-01, 9.444E-01, 9.446E-01,
     & 9.448E-01,9.451E-01, 9.453E-01, 9.455E-01, 9.457E-01, 9.460E-01,
     & 9.462E-01,9.464E-01, 9.466E-01, 9.468E-01, 9.470E-01, 9.472E-01,
     & 9.474E-01,9.477E-01, 9.479E-01, 9.481E-01, 9.483E-01, 9.485E-01,
     & 9.487E-01,9.489E-01, 9.491E-01, 9.493E-01, 9.495E-01, 9.496E-01,
     & 9.498E-01,9.500E-01, 9.502E-01, 9.504E-01, 9.506E-01, 9.508E-01,
     & 9.509E-01,9.511E-01, 9.513E-01, 9.515E-01, 9.517E-01, 9.518E-01,
     & 9.520E-01,9.522E-01, 9.524E-01, 9.525E-01, 9.527E-01, 9.529E-01,
     & 9.530E-01,9.532E-01, 9.534E-01, 9.535E-01, 9.537E-01, 9.538E-01,
     & 9.540E-01,9.542E-01, 9.543E-01, 9.545E-01, 9.546E-01, 9.548E-01,
     & 9.549E-01,9.551E-01, 9.552E-01, 9.554E-01, 9.555E-01, 9.557E-01,
     & 9.558E-01,9.560E-01, 9.561E-01, 9.563E-01, 9.564E-01, 9.566E-01,
     & 9.567E-01,9.568E-01, 9.570E-01, 9.571E-01, 9.573E-01, 9.574E-01,
     & 9.575E-01,9.577E-01, 9.578E-01, 9.579E-01, 9.581E-01, 9.582E-01,
     & 9.583E-01,9.585E-01, 9.586E-01, 9.587E-01, 9.588E-01, 9.590E-01,
     & 9.591E-01,9.592E-01, 9.593E-01, 9.595E-01, 9.596E-01, 9.597E-01,
     & 9.598E-01,9.600E-01, 9.601E-01, 9.602E-01, 9.603E-01, 9.604E-01,
     & 9.606E-01,9.607E-01, 9.608E-01, 9.609E-01, 9.610E-01, 9.611E-01,
     & 9.612E-01,9.614E-01, 9.615E-01, 9.616E-01, 9.617E-01, 9.618E-01,
     & 9.619E-01,9.620E-01, 9.621E-01, 9.622E-01, 9.623E-01, 9.624E-01,
     & 9.626E-01,9.627E-01, 9.628E-01, 9.629E-01/
 
      if (x .lt. .025) then
         xexpe1 = -x*(log(x)+.57721-x)
      else if (x .gt. 10.) then
         xexpe1 = 1.-.88/x
      else
        r = (x-dexe1off)/dexe1+1.0
        intr = int(r)
        frac = r-real(intr)
        xexpe1 = exe1data(intr)*(1.-frac)+exe1data(intr+1)*frac
      endif
 
      return
      end
 
c
      subroutine Vcaafill(iz)
c
c___ this calculates the collisional excitation or de-excitation
c___ from level i to level j.
c
      implicit real*8 (a-h,o-z)
      parameter(mpq=10)
      include 'mainstuf'
      include 'xc_stuf'
      include 'runstuf'
c
c___ Fill the arrays of excitation and de-excitation cross-sections
c___ versus energy.  Noting that the 1/e dependence of the xc
c___ implied that, e.g., the collision strengths, e*xc(e) is energy
c___ independent and thus the excitation = de-excitation
c

      do i=1,naatot(iz)
         do j=1,naatot(iz)
            Vcaa(i,j,iz) = 0.0
         enddo
      enddo
      if(.not.cxflag) return

      eip = eipz(iz)

      do 1000 i=1, naatot(iz)-1
         zeff = sqq(i,iz)
         do 900 j=i+1,naatot(iz)
            eji = evaa(j,iz)-evaa(i,iz) ! energy difference
            if(eji.le.exthr)goto 900
            fij=faa(i,j,iz)
            if(fij.le.0.d0) goto 900

            gij = gaa(i,iz)/gaa(j,iz)
            sumex=0.d0
            sumdx=0.d0

            do k=1,nefmax
               e = energy(k)
               if(e.gt.eji)then
                  g=0.15+0.28*log(e/eji)
                  sumex = sumex + fe(k)*denergy(k)
     &                 *1.40d-5*fij/eji*g
               endif
               g=0.15+0.28*log(e/eji+1.)
               sumdx = sumdx + fe(k)*denergy(k)
     &              *gij*1.40d-5*fij/eji*g
            enddo 

            Vcaa(i,j,iz)=sumex
            Vcaa(j,i,iz)=sumdx

            sum1=sumex

c ------   new gaunt factor ----------------
            n1=nloex(i,j,iz)    ! initial shell of lower level
            n2=nhiex(i,j,iz)    ! final shell
            if(isamp(1).eq.0.and.n1*n2.ne.0)then
               itr=mpq*(n1-1)-((n1-1)*n1)/2+n2-n1
               sumex=0.d0
               sumdx=0.d0
               do k=1,nefmax
                  e = energy(k)
                  y = e/eji
                  if(e.gt.eji)then
                     g=gcstjj(y,itr,iz)
                     sumex = sumex + fe(k)*denergy(k)
     &                    *1.40d-5*fij/eji*g
                  endif
                  g=gcstjj(y+1.,itr,iz)
                  sumdx = sumdx + fe(k)*denergy(k)
     &                 *gij*1.40d-5*fij/eji*g
               enddo 
               if(sumex.gt.0.d0)then
                  Vcaa(i,j,iz)=sumex
                  Vcaa(j,i,iz)=sumdx
               endif
            endif
c            if(sumex/sum1.gt.2.or.sumex/sum1.lt.0.5)
c     &           print *, iz, n1,n2, sumex/sum1
c-------------------------------------------------------
            if(npqe(i,j,iz).ne.0.and.isamp(1).eq.2) then
               pm=npqe(i,j,iz)
c               cxxx=7.096e-8/gaa(i,iz)/zeff**2
               cxxx=7.096e-8/zeff**2/(2.*n1**2)
               if(n1.le.3.and.(n2-n1).le.2) then
                  ii=2*(n1-1)+(n2-n1)
                  sumex=0.d0
                  sumdx=0.d0
                  do k=1,nefmax
                     e = energy(k)
                     y = e/eji
                     if(e.gt.eji)then
                        sumex = sumex + 
     &                       fe(k)*denergy(k)*cstrsamp(y,ii)*pm
                     endif
                     yp=y+1.
                     sumdx = sumdx + fe(k)*denergy(k)
     &                    *gij*cstrsamp(yp,ii)*pm
                  enddo
                  Vcaa(i,j,iz)=cxxx*sumex
                  Vcaa(j,i,iz)=cxxx*sumdx
               endif
            endif
c--------------------------------------------------------------
 900     continue
 1000 continue

      return
      end
c
cc* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      real*8 function cstrsamp(y,i)
c
c---- direct hydrogenic cross section - Golden etal
c                       Ap.J.Supp.Ser.,V45,p603,(1981)
c     cross-section(E)=8.7974e-17(pi*a0^2)/gi*(I_H/E)*cstrsamp/Zeff^2
c
      implicit real*8 (a-h,o-z)
      dimension ap(6),dp(6),c1(6),c2(6),a(6)
c

      data ap/ 2.220, 0.3560, 73.81, 10.18, 623.60,  76.23/
      data dp/ 0.8914,0.3546, 35.91, 18.25, 253.70, 168.90/
      data c1/-1.540,-0.4745,-14.27,-17.05, 351.00,-121.1/
      data c2/ 5.577, 0.7697,287.50, 34.29,2775.00, 304.4/
      data a/  0.75,  0.30,    1.18,  0.49,   1.42,   0.67/
c
      cstrsamp=ap(i)*log(y)+dp(i)+c1(i)/(y+a(i))+c2(i)/(y+a(i))**2
      return
      end

c---------------------------------------------------
      function cxx(i,j,iz)
      implicit real*8(a-h,o-z)
      include 'mainstuf'
      include 'popstuf'

      cxx=caa(i,j,iz)

      return
      end
c---------------------------------------------------
      function Vcxx(i,j,iz)
      implicit real*8(a-h,o-z)
      include 'mainstuf'
      include 'xc_stuf'

      Vcxx=Vcaa(i,j,iz)

      return
      end
c---------------------------------------------------
      function Vbetaxx(i,k,iz)
      implicit real*8(a-h,o-z)
      include 'mainstuf'
      include 'popstuf'
      include 'xc_stuf'

      Vbetaxx=0.d0
      if(.not.cxflag) return

      ntoti=indlev(i,iz)
      ntotj=indlev(k,iz-1)
      eii=elev(ntotj)-elev(ntoti)
      if(eii.le.exthr) return

      pn=petrn(i,k,iz)
      qq=qlev(ntoti)-1
      bfact=0.25*((100.*qq+91.)/(4.*qq+3.))**0.5 - 1.25
      ci=1.4679d-22*glev(ntoti)/glev(ntotj)

! Note that the 3-body recombination is done only with 
! the electrons on the grid. 
! No mixing with the thermal electrons and fe(E) electrons

      sum=0.d0
      do k2=1, nefmax
         do j=k2+1,nefmax
            e=energy(j)
            indx=indxc(e-energy(k2)-eii,fxxx)
            if(indx.ne.0)then
               uu=e/eii
               etot=e-eii
               e2=energy(k2)/etot
               if(e2.le.1.d0)then
                  wfact=(log(e/eii))**(bfact*eii/e)
                  xctot=2.218d-6*pn/eii*log(e/eii)*wfact
                  qj=qdiff(etot,uu,e2)/2.d0
                  sum=sum+xctot*qj*
     &                 fxxx*fe(k2)*denergy(k2)*denergy(j) ! note denergy(j)
               endif
            endif
         enddo
      enddo
      Vbetaxx=ci*sum
      return
      end
c---------------------------------------------------
      function Valfxx(i,k,iz)
      implicit real*8(a-h,o-z)
      include 'mainstuf'
      include 'popstuf'
      include 'xc_stuf'

      Valfxx=0.d0
      if(.not.cxflag) return
      if(k.ne.1)return

      ntoti=indlev(i,iz)
      ntotj=indlev(k,iz-1)
      eii=elev(ntotj)-elev(ntoti)
      if(eii.le.exthr) return

      qq=qlev(ntoti)
      pn=petrn(i,k,iz)
      glbygu = glev(ntoti)/glev(ntotj)

      sum = 0.d0
      do j=1,nefmax
         e = energy(j)
         de = denergy(j)
         sum = sum + fe(j)*
     &        1.692e-15*pn*glbygu*eii**2.5/qq/(eii+e)*de
      enddo

      Valfxx = sum

      return
      end

c---------------------------------------------------
      function Vsxx(i,k,iz)
      implicit real*8(a-h,o-z)
      include 'mainstuf'
      include 'popstuf'
      include 'runstuf'
      include 'xc_stuf'

      Vsxx=0.d0
      if(.not.cxflag) return

      ntoti=indlev(i,iz)
      ntotj=indlev(k,iz-1)
      eii=elev(ntotj)-elev(ntoti)
      if(eii.le.exthr)return

      pn=petrn(i,k,iz)
      qq=qlev(ntoti)-1
      bfact=0.25*((100.*qq+91.)/(4.*qq+3.))**0.5 - 1.25

      sum = 0.d0
      do j=1,nefmax
         e = energy(j)
         if(e .gt. eii) then
            wfact=(log(e/eii))**(bfact*eii/e)
            etot=e-eii
            uu=e/eii
            qsum=1.d0
            sum = sum+ qsum*2.218d-6*pn/eii*log(e/eii)*wfact
     &          *fe(j)*denergy(j)
         endif
      enddo

      Vsxx=sum
      return
      end
c
c---------------------------------------------------
      function betaxx(i,k,iz)
      implicit real*8(a-h,o-z)
      include 'mainstuf'
      include 'popstuf'

      betaxx = 0.d0

      ntoti=indlev(i,iz)
      ntotj=indlev(k,iz-1)
      eii=elev(ntotj)-elev(ntoti)
      ebt = eii/tev
      if(ebt.gt.500.d0)ebt=500.d0
      if(ebt.le.0.d0)return

      rc=sxx(i,k,iz)*exp(ebt)
      betaxx=1.656415d-22/tev**1.5*rc*glev(ntoti)/glev(ntotj)
      return
      end

c---------------------------------------------------
      function alfxx(i,k,iz)
      implicit real*8(a-h,o-z)
      include 'mainstuf'
      include 'popstuf'

c ... potential of the recombined level (i,iz+1)

      alfxx=0.
      estm=0.
c      if(k.ne.1)return

      ntoti=indlev(i,iz)
      ntotj=indlev(k,iz-1)
      eii=elev(ntotj)-elev(ntoti)
      if(eii.le.0.d0)return
      ebt = eii/tev
      if(ebt.gt.500.d0)ebt=500.d0
c
      qq=qlev(ntoti)
      pn=petrn(i,k,iz)
      glbygu = glev(ntoti)/glev(ntotj)
      alfxx=1.917d-15*pn*glbygu*sqrt(ebt)*xexpe1(ebt)*eii/qq
      estm=1.917d-15*pn*glbygu*(eii**2.5/tev**0.5)/qq

      return
      end

c---------------------------------------------------
c
      function sxx(i,k,iz)
      implicit real*8(a-h,o-z)
      parameter(cbar=2.3d0) !suggested value, cbar=2.77 for Lotz
      include 'mainstuf'
      include 'popstuf'
      data y1/3.5d5/,ems/5.11d5/,ems2/2.61121d5/

      sxx = 0.d0
      ntoti=indlev(i,iz)
      ntotj=indlev(k,iz-1)
      eii=elev(ntotj)-elev(ntoti)
      ebt = eii/tev
      if(ebt.gt.500.d0)ebt=500.d0
      if(ebt.le.0.d0) return

      pn=petrn(i,k,iz)
      qq=qlev(ntoti)-1.
      bfact=0.25*((100.*qq+91.)/(4.*qq+3.))**0.5 - 1.25
      wfact=log(1.+ 1./ebt)**(bfact/(1.+1./ebt))
      exebt=exp(-ebt)

      srel=0.d0
      if(relflag.and.(eii*100.+ems)/tev.le.100.) then
         sig1=(4.5d-14*log(y1/eii)/y1/eii)
         sig0=1.d-19/eii
         xxx=y1*(y1+2*ems)/eii/(eii+2.*ems)
         x0=(sig1-sig0)/log(xxx)
         x1=sig0-x0*log(eii*(eii+2.*ems))

         ex0=eii*1.d2
         exx=ex0+2*ems
         yy0=(ex0+tev)*exp(-ex0/tev)*log(ex0*exx)
     &        +(tev-2.*ems)*exp(2.*ems/tev)*expe1(exx/tev)
     &        +2.*tev*exp(-ex0/tev)
         yy1=(tev+ex0)*exp(-ex0/tev)
         srel=6.69d7/tev**0.5*(x0*yy0+x1*yy1)
      endif

c Burgess & Chidichimo
      if (isamp(2).eq.1) then
         sxx = 1.0891e-06*cbar/(sqrt(tev)*eii)*e1(ebt,exebt)*pn*wfact
      else if(isamp(2).eq.2) then
c  Lotz'
         sxx = 2.97e-6/(sqrt(tev)*eii)*e1(ebt,exebt)* pn
      else if(isamp(2).eq.3) then
c  Sampson
         ir=morb(i,iz)
         ze1=sqq(i,iz)
         rc=rcisamp(ir,eii,tev,ze1)
         sxx=rc*pn
      endif

      sxx=sxx+max(0.d0,srel*pn)

      return
      end
c
      real*8 function rcisamp (n,ei,tem,ze1)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    calculate the ionization rate using cross section fitted to the
c    results of sampson and golden(1972) (g. zimmerman)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit real*8 (a-h,o-z)

      dimension a(10),b(10),c(10),gbn(10)

      data(a(i),i=1,10)/3.1287,5.0048,7.0391,9.2107,11.4354,13.3434,
     115.2657,17.2635,19.3700,21.5106 /
      data(b(i),i=1,10)/-4.7764,-7.3133,-10.2440,-13.4870,-16.8678,
     1-19.7029,-22.5916,-25.6635,-28.9884,-32.4002 /
      data(c(i),i=1,10)/1.6476,2.3084,3.2048,4.2763,5.4323,6.3594,
     17.3259,8.4000,9.6183,10.8896 /
      data(gbn(i),i=1,10)/0.8675,0.932,0.952,0.960,0.965,0.969,0.972,
     10.975,0.978,0.981 /

      external erf, erfc, expe1, expund

      xn=n
      eit=ei/tem
      sqeit=sqrt(eit)
      f=sqrt(3.0)/(2.0*3.1415)*(128.0/9.0)*xn**3*gbn(n)
      an=a(n)/eit
      bn=b(n)/sqrt(eit)
      cn=c(n)
      r=5.88456d-9*sqrt(tem)*f/ze1**4

      if(sqeit.ge.9.0) go to 102

      eerf=1.0-erf(sqeit)

      if(eerf.le.1.0e-4)  eerf=erfc(sqeit)

      expe11=expe1(eit)
      go to 103

  102 continue

      eerf=0.0
      expe11=0.0

  103 continue

      eerf=0.8862269*eerf
      rcisamp=r*(an*expund(-eit,-80.d0)+2.0*bn*eerf+cn*expe11)
      rcisamp=rcisamp*eit*eit
      return
      end
c

      function autoxx(i,k,iz)
      implicit real*8(a-h,o-z)
      include 'mainstuf'
      include 'popstuf'

      autoxx=autoaa(i,k,iz)

      return
      end
c
c
      function Vcapexx(i,k,iz)
      implicit real*8(a-h,o-z)
      include 'mainstuf'
      include 'popstuf'
      include 'xc_stuf'

      Vcapexx=0.d0
      if(.not.cxflag) return

      ntoti=indlev(i,iz)
      ntotj=indlev(k,iz-1)
      eii=elev(ntoti)-elev(ntotj)

      j=indxc(eii,fxxx)
      glbygu=glev(ntoti)/glev(ntotj)
      Vcapexx=1.4679595906d-22*autoxx(i,k,iz)*glbygu*fxxx

      return
      end

c
c
      function Vbstmxx(i,k,iz)     
      implicit real*8(a-h,o-z)
      include 'mainstuf'
      include 'runstuf'
      include 'popstuf'
      include 'xc_stuf'

      Vbstmxx =0.d0
      if(.not.cxflag) return

      ntoti=indlev(i,iz)
      ntotj=indlev(k,iz-1)
      eii=elev(ntotj)-elev(ntoti)
      if(eii.le.exthr) return


      ebt = eii/tev
      glbygu = glev(ntoti)/glev(ntotj)
      qq= qlev(ntoti)
      pn = petrn(i,k,iz)

      sum=0.d0
      do j=1,nefmax
         e = energy(j)
         de = denergy(j)
         if (radflag .ne. 'trfile') then
            if(trev.eq.0.d0) return
            field_e = bbfield((e+eii)/trev)
         else
            field_e = field(e+eii)
         endif
         sum= sum+ fe(j)*
     1        8.167e-12*glbygu*eii**2.5/qq
     2        /(e+eii)**4*field_e*de*pn
      enddo
      Vbstmxx =sum

      return
      end
c
c
      function bstmxx(i,k,iz)     
      implicit real*8(a-h,o-z)
      include 'mainstuf'
      include 'runstuf'
      include 'timestuf'
      include 'popstuf'
c
      bstmxx=0.d0
      estm=0.d0

      ntoti=indlev(i,iz)
      ntotj=indlev(k,iz-1)
      eii= elev(ntotj)-elev(ntoti)
      if(eii.le.0.d0)return

      ebt = eii/tev
      glbygu = glev(ntoti)/glev(ntotj)
      qq=qlev(ntoti)
      pn=petrn(i,k,iz)
c
c     stimulated recombination 
c
      bc=9.162e-12/qq*glbygu*eii**2.5/tev**1.5*pn

      ymin=eii
      ylast=5.*eii
      ny=51
      dely=(ylast-ymin)/real(ny-1)

      if (radflag .ne. 'trfile') then

         dint = 0.5*exp(-(ymin-eii)/tev)*bbfield(ymin/trev)/ymin**4
         sum = dint
         esum= dint*ymin
         do 100 j=2,ny-1
            y2 = ymin+dely*real(j-1)
            dint = exp(-(y2-eii)/tev)*bbfield(y2/trev)/y2**4
            sum = sum + dint
            esum = esum + dint*y2
 100     continue
         dint=0.5*exp(-(ylast-eii)/tev)*bbfield(ylast/trev)/ylast**4
         sum = sum + dint
         esum= esum+ dint*ylast
      else

         if(xev(ngrid).lt.ymin) then
            return
         endif
         if(xev(1).gt.ymin) then
            ymin=xev(1)
            ylast=xev(ngrid)
         endif
         if(xev(ngrid).lt.ylast) then
            ylast=xev(ngrid)
         endif

         dely=(ylast-ymin)/real(ny-1)
         dint = 0.5*exp(-(ymin-eii)/tev)*field(ymin)/ymin**4
         sum = dint
         esum = dint*ymin
         do 200 j=2,ny-1
            y2 = ymin+dely*real(j-1)
            dint = exp(-(y2-eii)/tev)*field(y2)/y2**4
            sum = sum + dint
            esum = esum + dint * y2
 200     continue
         dint = 0.5*exp(-(ylast-eii)/tev)*field(ylast)/ylast**4
         sum = sum + dint
         esum= esum+ dint*ylast

      endif
c
      bstmxx = bc*sum*dely
      estm   = bc*esum*dely

      return
      end
c
      function bcontxx(i,k,iz)     
      implicit real*8(a-h,o-z)
      include 'mainstuf'
      include 'runstuf'
      include 'popstuf'
      include 'timestuf'
c
c     bbfiled is in the unit of erg/cm2/sec/Hz
c
      bcontxx=0.00
      econt = 0.d0

      ntoti=indlev(i,iz)
      ntotj=indlev(k,iz-1)
      eii= elev(ntotj)-elev(ntoti)
      if(eii.le.0.d0) return

      ebt = eii/tev
      qq= qlev(ntoti)
      pn = petrn(i,k,iz)

      if (radflag .ne. 'trfile') then

        tbyr = tev/trev
        bc=5.53e+10*ebt**2.5/qq/sqrt(tev)*pn
c
        ylast=max(10.*ebt,2.)
        ymin=ebt
        ny=51
        dely=(ylast-ymin)/real(ny-1)

        sum = 0.5*bbfield(ymin*tbyr)/ymin**4
        esum = 0.5*bbfield(ymin*tbyr)/ymin**3
        do 100 j=2,ny-1
          y2 = ymin+dely*real(j-1)
          sum = sum+bbfield(y2*tbyr)/y2**4
          esum= esum+bbfield(y2*tbyr)/y2**3
  100   continue
        sum = sum+bbfield(ylast*tbyr)/ylast**4*.5
        esum= esum+bbfield(ylast*tbyr)/ylast**3*.5

      else

         bc=5.53e+10*eii**2.5/qq*pn
         ymin=eii
         ylast=10.*eii
         
         if(xev(ngrid).lt.ymin) then
            return
         endif
         if(xev(1).gt.ymin) then
            ymin=xev(1)
            ylast=xev(ngrid)
         endif
         if(xev(ngrid).lt.ylast) then
            ylast=xev(ngrid)
         endif

         ny=51
         dely=(ylast-ymin)/real(ny-1)

         sum = 0.5*field(ymin)/ymin**4
         esum = 0.5*field(ymin)/ymin**3
         do 200 j=2,ny-1
            y2 = ymin+dely*real(j-1)
            sum = sum+field(y2)/y2**4
            esum = esum+field(y2)/y2**3
 200     continue
         sum = sum+field(ylast)/ylast**4*.5
         esum = esum+field(ylast)/ylast**3*.5                  
      endif

      bcontxx = bc*sum*dely
      econt   = bc*esum*dely
c
      return
      end
c
c***********************************************************
      function qdiff(etot,u,e2)
      implicit double precision(a-h,o-z)
! Secondary electron energy distribution  Clark et al. (1991)
      qdiff=(30.8+160.*(u**2-1.)*(e2-0.5)**4)/(14.4+u**2)/etot
! constant approximation
!      qdiff=2.d0/etot
      return
      end

c---------------------------------------------------

      function exp1x(x)

c ... exp(x) times
c ... Rational approximation to 1st exponential integral (within 5.e-5)
c ... (ref. Abramowitz and Stegun, p. 231

      real*8 x,x0,x1,a0,a1,a2,a3,a4,a5,b1,b2,b3,b4,exp1x

      data x0,x1 / 1.d-4, 1.d0 /

      data a0,a1,a2,a3,a4,a5
     &     / -.57721566d0, .99999193d0,-.24991055d0,
     &        .05519968d0,-.00976004d0, .00107857d0 /

      data b1,b2,b3,b4 / 2.334733d0,0.250621d0,3.330657d0,1.681534d0 /

      if (x.lt.x0) then
         exp1x = (1.d0 + x) * (-log(x) + a0 + x*(a1 + x*a2))
      elseif (x.lt.x1) then
         exp1x = exp(x) *
     &           (-log(x) + a0 + x*(a1 + x*(a2+ x*(a3+ x*(a4 + x*a5)))))
      else
         exp1x = (b2 + x*(b1 + x)) / (b4 + x*(b3 + x)) / x
      endif

      return
      end

c---------------------------------------------------
      subroutine caleipd(S_E,P_F)
      implicit real*8(a-h,o-z)
      include 'mainstuf'
      include 'popstuf'

c     setup internal energy consistent with IPD

      igr0=indlev(1,ionow)
      do 11 i=iatomz,1,-1 
         do 12 j=1, nlev(i)
            ilev=indlev(j,i)
            if(indipd(ilev).ne.0) then
               efgrdp(ilev)=elev(ilev)-elev(igr0)
            else
               efgrdp(ilev)=0.d0
            endif
 12      continue
 11   continue


      S_E=0.d0
      P_F=0.d0
      do i=1,numpop
         iso=ionstg(i)
         if(popout(i).gt.1.d-70)then
            S_E = S_E+popout(i)*efgrdp(i) ! specific internal energy
            if(indl2p(i).ne.0) then
               P_F = P_F+glev(i)*exp(-efgrdp(i)/tev) ! partition function
            endif
         endif
      enddo
      return
      end


      subroutine pwrloss(poptot)
      implicit real*8(a-h,o-z)
      include 'mainstuf'
      include 'popstuf'

      common/NLTE_pw/pwrbb,pwrbf,pwrff,pwrtt

      common /linmax/totmax

      iz=iatomz
      pden=denow
      sumbb=0.d0
      sumbf=0.d0
      totmax=0.d0

      tev=tevm(1)

      do kk = 1, iz
c   compute line emission 
         call embb(kk,poptot,totbb) ! totbb [eV/s/atom/Ne]
c   compute bound-free emission
         call embf(kk,poptot,totbf)
         sumbb=sumbb+totbb
         sumbf=sumbf+totbf
      enddo
c   compute free-free emission  ! [ergs/sec/cm-3]
      pwrbb=sumbb*denow*poptot*1.6022e-12
      pwrbf=sumbf*denow*poptot*1.6022e-12
      call emff(poptot,totff)
      pwrff=totff*denow*poptot*1.6022e-12
      pwrtt=pwrbb+pwrbf+pwrff
      return
      end

      subroutine embb(nion,popt,totrad)
      implicit real*8(a-h,o-z)

c     this routine returns the total bound-bound line emission rates
c     of the ionization stage "nion" at the case "inow"
c     in units of eV/s/atom/Ne

      include 'mainstuf'
      include 'popstuf'
      common /linmax/totmax
      common /emtrn/mtrn(50,100),etrn(13,50,100),ptrn(13,50,100)

      tot = 0.d0
      do iup=2, nlev(nion)
         ilv=indlev(iup,nion)
         pup=popout(ilv)
         do ilo=1, iup-1
            ilw=indlev(ilo,nion)
            aep=aratexx(iup,ilo,nion)! A-value weighted by escape probablity.
            etr=elev(ilv)-elev(ilw)! transition energy in eV
            tot=tot+etr*aep*pup          ! line radiation cooling rates eV/s
            totmax=max(totmax,tot)
         enddo
      enddo

      totrad=tot/popt/denow

      return
      end

      subroutine embf(nion,popt,totrad)
      implicit real*8(a-h,o-z)

c     this routine returns the total free-bound line emission rates
c     of the ionization stage "nion" at the case "inow"
c     we include the emission due to recombination from the next ground state.
c     in units of eV/s/atom/Ne

      include 'mainstuf'
      include 'popstuf'

      nlv = nlev(nion)
      ilv=indlev(1,nion-1)
      gup = glev(ilv)
      pup = popout(ilv)
      zion= (iatomz-nion+1)*1.d0
      const = 1.91e-15/tev**0.5*(pup/popt)/gup

      tot = 0.d0
      do ilo=1, nlbn(nion)
         ilw=indlev(ilo,nion)
         glo=glev(ilw)
         eip=elev(ilv)-elev(ilw)
         if(eip.gt.0.d0)then
            tot=tot+const*glo*eip**2.5/zion
         endif
      enddo
      totrad=tot
      return
      end


      subroutine emff(popt,totrad)
      implicit real*8(a-h,o-z)

c     this routine returns the total free-bound line emission rates
c     of the ionization stage "nion" at the case "inow"
c     we include the emission due to recombination from the next ground state.
c     in units of eV/s/atom/Ne

      include 'mainstuf'
      include 'popstuf'

      zbr2=0.d0
      do kk = 0, iatomz
         sum=0.d0
         do i=1, nlev(kk)
            ilv=indlev(i,kk)
            sum=sum+popout(ilv)
         enddo
         zbr2=zbr2+(iatomz-kk)**2 * sum
      enddo
      zbr2=zbr2/popt ! fractional population 

      totrad=0.d0
      do ii=1, ntem
         fevv=fnem(ii)
         tevv=tevm(ii)
c         gff=1.d0+0.44*exp(-0.25*(0.25+log10(13.6*zbr2/tevv)))
         gff=1.d0+0.44*exp(-0.25*(0.25+log10(13.6*zbr2/tevv))**2) ! correction 06Oct24
         totrad=totrad+fevv*tevv**0.5*gff
      enddo
      totrad=totrad*(9.55e-14*zbr2)

      return
      end

