

      SUBROUTINE LSODE (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     .                  ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
C
C  Declare arguments.
C
      implicit real*8 (a-h,o-z)
      EXTERNAL F, JAC
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
      REAL*8 Y, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW)
C
C  Declare externals (for linear system solver).
C
      external prepj, solsy
C
C  Declare local variables.
C
      integer illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
     1   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns
      integer icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     1   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, i1, i2, iflag, imxer, kgo, lf0,
     1   leniw, lenrw, lenwm, ml, mord, mu, mxhnl0, mxstp0
      real*8 rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real*8 atoli, ayi, big, ewti, h0, hmax, hmx, rh, rtoli,
     1   tcrit, tdist, tnext, tol, tolsf, tp, size, sum, w0,
     2   r1mach, vnorm
      dimension mord(2)
      logical ihit
c-----------------------------------------------------------------------
c the following internal common block contains
C (a) variables which are local to any subroutine but whose values must
c     be preserved between calls to the routine (own variables), and
c (b) variables which are communicated between subroutines.
c the structure of the block is as follows..  all real variables are
C listed first, followed by all integers.  within each type, the
C variables are grouped with those local to subroutine lsode first,
C then those local to subroutine stode, and finally those used
c for communication.  the block is declared in subroutines
c lsode, intdy, stode, prepj, and solsy.  groups of variables are
c replaced by dummy arrays in the common declarations in routines
c where those variables are not used.
c-----------------------------------------------------------------------
      common /ls0001R/ rowns(209),
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      common /ls0001I/
     2   illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
     3   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
c
      data  mord(1),mord(2)/12,5/, mxstp0/15000/, mxhnl0/10/
      data illin/0/, ntrep/0/
c-----------------------------------------------------------------------
c block a.
c this code block is executed on every call.
c it tests istate and itask for legality and branches appropriately.
c if istate .gt. 1 but the flag init shows that initialization has
c not yet been done, an error return occurs.
c if istate = 1 and tout = t, jump to block g and return immediately.
c-----------------------------------------------------------------------
C
C***FIRST EXECUTABLE STATEMENT  LSODE
      if (istate .lt. 1 .or. istate .gt. 3) go to 601
      if (itask .lt. 1 .or. itask .gt. 5) go to 602
      if (istate .eq. 1) go to 10
      if (init .eq. 0) go to 603
      if (istate .eq. 2) go to 200
      go to 20
 10   init = 0
      if (tout .eq. t) go to 430
 20   ntrep = 0
c-----------------------------------------------------------------------
c block b.
c the next code block is executed for the initial call (istate = 1),
c or for a continuation call with parameter changes (istate = 3).
c it contains checking of all inputs and various initializations.
c
c first check legality of the non-optional inputs neq, itol, iopt,
c mf, ml, and mu.
c-----------------------------------------------------------------------
      if (neq(1) .le. 0) go to 604
      if (istate .eq. 1) go to 25
      if (neq(1) .gt. n) go to 605
 25   n = neq(1)
      if (itol .lt. 1 .or. itol .gt. 4) go to 606
      if (iopt .lt. 0 .or. iopt .gt. 1) go to 607
      meth = mf/10
      miter = mf - 10*meth
      if (meth .lt. 1 .or. meth .gt. 2) go to 608
      if (miter .lt. 0 .or. miter .gt. 5) go to 608
      if (miter .le. 3) go to 30
      ml = iwork(1)
      mu = iwork(2)
      if (ml .lt. 0 .or. ml .ge. n) go to 609
      if (mu .lt. 0 .or. mu .ge. n) go to 610
 30   continue
c next process and check the optional inputs. --------------------------
      if (iopt .eq. 1) go to 40
      maxord = mord(meth)
      mxstep = mxstp0
      mxhnil = mxhnl0
      if (istate .eq. 1) h0 = 0.0d0
      hmxi = 0.0d0
      hmin = 0.0d0
      go to 60
 40   maxord = iwork(5)
      if (maxord .lt. 0) go to 611
      if (maxord .eq. 0) maxord = 100
      maxord = min0(maxord,mord(meth))
      mxstep = iwork(6)
      if (mxstep .lt. 0) go to 612
      if (mxstep .eq. 0) mxstep = mxstp0
      mxhnil = iwork(7)
      if (mxhnil .lt. 0) go to 613
      if (mxhnil .eq. 0) mxhnil = mxhnl0
      if (istate .ne. 1) go to 50
      h0 = rwork(5)
      if ((tout - t)*h0 .lt. 0.0d0) go to 614
 50   hmax = rwork(6)
      if (hmax .lt. 0.0d0) go to 615
      hmxi = 0.0d0
      if (hmax .gt. 0.0d0) hmxi = 1.0e0/hmax
      hmin = rwork(7)
      if (hmin .lt. 0.0d0) go to 616
c-----------------------------------------------------------------------
c set work array pointers and check lengths lrw and liw.
c pointers to segments of rwork and iwork are named by prefixing l to
c the name of the segment.  e.g., the segment yh starts at rwork(lyh).
c segments of rwork (in order) are denoted  yh, wm, ewt, savf, acor.
c-----------------------------------------------------------------------
 60   lyh = 21
      if (istate .eq. 1) nyh = n
      lwm = lyh + (maxord + 1)*nyh
      if (miter .eq. 0) lenwm = 0
      if (miter .eq. 1 .or. miter .eq. 2) lenwm = n*n + 2
      if (miter .eq. 3) lenwm = n + 2
      if (miter .ge. 4) lenwm = (2*ml + mu + 1)*n + 2
      lewt = lwm + lenwm
      lsavf = lewt + n
      lacor = lsavf + n
      lenrw = lacor + n - 1
      iwork(17) = lenrw
      liwm = 1
      leniw = 20 + n
      if (miter .eq. 0 .or. miter .eq. 3) leniw = 20
      iwork(18) = leniw
      if (lenrw .gt. lrw) go to 617
      if (leniw .gt. liw) go to 618
c check rtol and atol for legality. ------------------------------------
      rtoli = rtol(1)
      atoli = atol(1)
      do 70 i = 1,n
        if (itol .ge. 3) rtoli = rtol(i)
        if (itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
        if (rtoli .lt. 0.0d0) go to 619
        if (atoli .lt. 0.0d0) go to 620
 70     continue
      if (istate .eq. 1) go to 100
c if istate = 3, set flag to signal parameter changes to stode. --------
      jstart = -1
      if (nq .le. maxord) go to 90
c maxord was reduced below nq.  copy yh(*,maxord+2) into savf. ---------
      do 80 i = 1,n
 80     rwork(i+lsavf-1) = rwork(i+lwm-1)
c reload wm(1) = rwork(lwm), since lwm may have changed. ---------------
 90   if (miter .gt. 0) rwork(lwm) = sqrt(uround)
      if (n .eq. nyh) go to 200
c neq was reduced.  zero part of yh to avoid undefined references. -----
      i1 = lyh + l*nyh
      i2 = lyh + (maxord + 1)*nyh - 1
      if (i1 .gt. i2) go to 200
      do 95 i = i1,i2
 95     rwork(i) = 0.0d0
      go to 200
c-----------------------------------------------------------------------
c block c.
c the next block is for the initial call only (istate = 1).
c it contains all remaining initializations, the initial call to f,
c and the calculation of the initial step size.
c the error weights in ewt are inverted after being loaded.
c-----------------------------------------------------------------------
 100  uround = r1mach(4)
      tn = t
      if (itask .ne. 4 .and. itask .ne. 5) go to 110
      tcrit = rwork(1)
      if ((tcrit - tout)*(tout - t) .lt. 0.0d0) go to 625
      if (h0 .ne. 0.0d0 .and. (t + h0 - tcrit)*h0 .gt. 0.0d0)
     1   h0 = tcrit - t
 110  jstart = 0
      if (miter .gt. 0) rwork(lwm) = sqrt(uround)
      nhnil = 0
      nst = 0
      nje = 0
      nslast = 0
      hu = 0.0d0
      nqu = 0
      ccmax = 0.3e0
      maxcor = 3
      msbp = 20
      mxncf = 10
c initial call to f.  (lf0 points to yh(*,2).) -------------------------
      lf0 = lyh + nyh
      call f (neq, t, y, rwork(lf0))
      nfe = 1
c load the initial value vector in yh. ---------------------------------
      do 115 i = 1,n
 115    rwork(i+lyh-1) = y(i)
c load and invert the ewt array.  (h is temporarily set to 1.0.) -------
      nq = 1
      h = 1.0e0
      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      do 120 i = 1,n
        if (rwork(i+lewt-1) .le. 0.0d0) go to 621
 120    rwork(i+lewt-1) = 1.0e0/rwork(i+lewt-1)
c-----------------------------------------------------------------------
c the coding below computes the step size, h0, to be attempted on the
c first step, unless the user has supplied a value for this.
c first check that tout - t differs significantly from zero.
c a scalar tolerance quantity tol is computed, as max(rtol(i))
c if this is positive, or max(atol(i)/abs(y(i))) otherwise, adjusted
c so as to be between 100*uround and 1.0e-3.
c then the computed value h0 is given by..
c                                      neq
c   h0**2 = tol / ( w0**-2 + (1/neq) * sum ( f(i)/ywt(i) )**2  )
c                                       1
c where   w0     = max ( abs(t), abs(tout) ),
c         f(i)   = i-th component of initial value of f,
c         ywt(i) = ewt(i)/tol  (a weight for y(i)).
c the sign of h0 is inferred from the initial values of tout and t.
c-----------------------------------------------------------------------
      if (h0 .ne. 0.0d0) go to 180
      tdist = abs(tout - t)
      w0 = max(abs(t),abs(tout))
      if (tdist .lt. 2.0e0*uround*w0) go to 622
      tol = rtol(1)
      if (itol .le. 2) go to 140
      do 130 i = 1,n
 130    tol = max(tol,rtol(i))
 140  if (tol .gt. 0.0d0) go to 160
      atoli = atol(1)
      do 150 i = 1,n
        if (itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
        ayi = abs(y(i))
        if (ayi .ne. 0.0d0) tol = max(tol,atoli/ayi)
 150    continue
 160  tol = max(tol,100.0d0*uround)
      tol = min(tol,0.001e0)
      sum = vnorm (n, rwork(lf0), rwork(lewt))
      sum = 1.0e0/(tol*w0*w0) + tol*sum**2
      h0 = 1.0e0/sqrt(sum)
      h0 = min(h0,tdist)
      h0 = sign(h0,tout-t)
c adjust h0 if necessary to meet hmax bound. ---------------------------
 180  rh = abs(h0)*hmxi
      if (rh .gt. 1.0e0) h0 = h0/rh
c load h with h0 and scale yh(*,2) by h0. ------------------------------
      h = h0
      do 190 i = 1,n
 190    rwork(i+lf0-1) = h0*rwork(i+lf0-1)
      go to 270
c-----------------------------------------------------------------------
c block d.
c the next code block is for continuation calls only (istate = 2 or 3)
c and is to check stop conditions before taking a step.
c-----------------------------------------------------------------------
 200  nslast = nst
      go to (210, 250, 220, 230, 240), itask
 210  if ((tn - tout)*h .lt. 0.0d0) go to 250
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      if (iflag .ne. 0) go to 627
      t = tout
      go to 420
 220  tp = tn - hu*(1.0e0 + 100.0d0*uround)
      if ((tp - tout)*h .gt. 0.0d0) go to 623
      if ((tn - tout)*h .lt. 0.0d0) go to 250
      go to 400
 230  tcrit = rwork(1)
      if ((tn - tcrit)*h .gt. 0.0d0) go to 624
      if ((tcrit - tout)*h .lt. 0.0d0) go to 625
      if ((tn - tout)*h .lt. 0.0d0) go to 245
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      if (iflag .ne. 0) go to 627
      t = tout
      go to 420
 240  tcrit = rwork(1)
      if ((tn - tcrit)*h .gt. 0.0d0) go to 624
 245  hmx = abs(tn) + abs(h)
      ihit = abs(tn - tcrit) .le. 100.0d0*uround*hmx
      if (ihit) go to 400
      tnext = tn + h*(1.0e0 + 4.0e0*uround)
      if ((tnext - tcrit)*h .le. 0.0d0) go to 250
      h = (tcrit - tn)*(1.0e0 - 4.0e0*uround)
      if (istate .eq. 2) jstart = -2
c-----------------------------------------------------------------------
c block e.
c the next block is normally executed for all calls and contains
c the call to the one-step core integrator stode.
c
c this is a looping point for the integration steps.
c
c first check for too many steps being taken, update ewt (if not at
c start of problem), check for too much accuracy being requested, and
c check for h below the roundoff level in t.
c-----------------------------------------------------------------------
 250  continue
      if ((nst-nslast) .ge. mxstep) go to 500
      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      do 260 i = 1,n
        if (rwork(i+lewt-1) .le. 0.0d0) go to 510
 260    rwork(i+lewt-1) = 1.0e0/rwork(i+lewt-1)
 270  tolsf = uround*vnorm (n, rwork(lyh), rwork(lewt))
      if (tolsf .le. 1.0e0) go to 280
      tolsf = tolsf*2.0e0
      if (nst .eq. 0) go to 626
      go to 520
 280  if ((tn + h) .ne. tn) go to 290
      nhnil = nhnil + 1
      if (nhnil .gt. mxhnil) go to 290
      call xerrwv('lsode--  warning..internal t (=r1) and h (=r2) are',
     1   50, 101, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(
     1  '     such that in the machine, t + h = t on the next step  ',
     1   60, 101, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv('      (h = step size). solver will continue anyway',
     1   50, 101, 0, 0, 0, 0, 2, tn, h)
      if (nhnil .lt. mxhnil) go to 290
      call xerrwv('lsode--  above warning has been issued i1 times. ',
     1   50, 102, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv('      it will not be issued again for this problem',
     1   50, 102, 0, 1, mxhnil, 0, 0, 0.0d0, 0.0d0)
 290  continue
c-----------------------------------------------------------------------
c     call stode(neq,y,yh,nyh,yh,ewt,savf,acor,wm,iwm,f,jac,prepj,solsy)
c-----------------------------------------------------------------------
      call stode (neq, y, rwork(lyh), nyh, rwork(lyh), rwork(lewt),
     1   rwork(lsavf), rwork(lacor), rwork(lwm), iwork(liwm),
     2   f, jac, prepj, solsy)
      kgo = 1 - kflag
      go to (300, 530, 540), kgo
c-----------------------------------------------------------------------
c block f.
c the following block handles the case of a successful return from the
c core integrator (kflag = 0).  test for stop conditions.
c-----------------------------------------------------------------------
 300  init = 1
      go to (310, 400, 330, 340, 350), itask
c itask = 1.  if tout has been reached, interpolate. -------------------
 310  if ((tn - tout)*h .lt. 0.0d0) go to 250
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      t = tout
      go to 420
c itask = 3.  jump to exit if tout was reached. ------------------------
 330  if ((tn - tout)*h .ge. 0.0d0) go to 400
      go to 250
c itask = 4.  see if tout or tcrit was reached.  adjust h if necessary.
 340  if ((tn - tout)*h .lt. 0.0d0) go to 345
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      t = tout
      go to 420
 345  hmx = abs(tn) + abs(h)
      ihit = abs(tn - tcrit) .le. 100.0d0*uround*hmx
      if (ihit) go to 400
      tnext = tn + h*(1.0e0 + 4.0e0*uround)
      if ((tnext - tcrit)*h .le. 0.0d0) go to 250
      h = (tcrit - tn)*(1.0e0 - 4.0e0*uround)
      jstart = -2
      go to 250
c itask = 5.  see if tcrit was reached and jump to exit. ---------------
 350  hmx = abs(tn) + abs(h)
      ihit = abs(tn - tcrit) .le. 100.0d0*uround*hmx
c-----------------------------------------------------------------------
c block g.
c the following block handles all successful returns from lsode.
c if itask .ne. 1, y is loaded from yh and t is set accordingly.
c istate is set to 2, the illegal input counter is zeroed, and the
c optional outputs are loaded into the work arrays before returning.
c if istate = 1 and tout = t, there is a return with no action taken,
c except that if this has happened repeatedly, the run is terminated.
c-----------------------------------------------------------------------
 400  do 410 i = 1,n
 410    y(i) = rwork(i+lyh-1)
      t = tn
      if (itask .ne. 4 .and. itask .ne. 5) go to 420
      if (ihit) t = tcrit
 420  istate = 2
      illin = 0
      rwork(11) = hu
      rwork(12) = h
      rwork(13) = tn
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      return
c
 430  ntrep = ntrep + 1
      if (ntrep .lt. 5) return
      call xerrwv(
     1  'lsode--  repeated calls with istate = 1 and tout = t (=r1)  ',
     1   60, 301, 0, 0, 0, 0, 1, t, 0.0d0)
      go to 800
c-----------------------------------------------------------------------
c block h.
c the following block handles all unsuccessful returns other than
c those for illegal input.  first the error message routine is called.
c if there was an error test or convergence test failure, imxer is set.
c then y is loaded from yh, t is set to tn, and the illegal input
c counter illin is set to 0.  the optional outputs are loaded into
c the work arrays before returning.
c-----------------------------------------------------------------------
c the maximum number of steps was taken before reaching tout. ----------
 500  call xerrwv('lsode--  at current t (=r1), mxstep (=i1) steps   ',
     1   50, 201, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv('      taken on this call before reaching tout     ',
     1   50, 201, 0, 1, mxstep, 0, 1, tn, 0.0d0)
      istate = -1
      go to 580
c ewt(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  ewti = rwork(lewt+i-1)
      call xerrwv('lsode--  at t (=r1), ewt(i1) has become r2 .le. 0.',
     1   50, 202, 0, 1, i, 0, 2, tn, ewti)
      istate = -6
      go to 580
c too much accuracy requested for machine precision. -------------------
 520  call xerrwv('lsode--  at t (=r1), too much accuracy requested  ',
     1   50, 203, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv('      for precision of machine..  see tolsf (=r2) ',
     1   50, 203, 0, 0, 0, 0, 2, tn, tolsf)
      rwork(14) = tolsf
      istate = -2
      go to 580
c kflag = -1.  error test failed repeatedly or with abs(h) = hmin. -----
 530  call xerrwv('lsode--  at t(=r1) and step size h(=r2), the error',
     1   50, 204, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv('      test failed repeatedly or with abs(h) = hmin',
     1   50, 204, 0, 0, 0, 0, 2, tn, h)
      istate = -4
      go to 560
c kflag = -2.  convergence failed repeatedly or with abs(h) = hmin. ----
 540  call xerrwv('lsode--  at t (=r1) and step size h (=r2), the    ',
     1   50, 205, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv('      corrector convergence failed repeatedly     ',
     1   50, 205, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv('      or with abs(h) = hmin   ',
     1   30, 205, 0, 0, 0, 0, 2, tn, h)
      istate = -5
c compute imxer if relevant. -------------------------------------------
 560  big = 0.0d0
      imxer = 1
      do 570 i = 1,n
        size = abs(rwork(i+lacor-1)*rwork(i+lewt-1))
        if (big .ge. size) go to 570
        big = size
        imxer = i
 570    continue
      iwork(16) = imxer
c set y vector, t, illin, and optional outputs. ------------------------
 580  do 590 i = 1,n
 590    y(i) = rwork(i+lyh-1)
      t = tn
      illin = 0
      rwork(11) = hu
      rwork(12) = h
      rwork(13) = tn
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      return
c-----------------------------------------------------------------------
c block i.
c the following block handles all error returns due to illegal input
c (istate = -3), as detected before calling the core integrator.
c first the error message routine is called.  then if there have been
c 5 consecutive such returns just before this call to the solver,
c the run is halted.
c-----------------------------------------------------------------------
 601  call xerrwv('lsode--  istate (=i1) illegal ',
     1   30, 1, 0, 1, istate, 0, 0, 0.0d0, 0.0d0)
      go to 700
 602  call xerrwv('lsode--  itask (=i1) illegal  ',
     1   30, 2, 0, 1, itask, 0, 0, 0.0d0, 0.0d0)
      go to 700
 603  call xerrwv('lsode--  istate .gt. 1 but lsode not initialized  ',
     1   50, 3, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      go to 700
 604  call xerrwv('lsode--  neq (=i1) .lt. 1     ',
     1   30, 4, 0, 1, neq(1), 0, 0, 0.0d0, 0.0d0)
      go to 700
 605  call xerrwv('lsode--  istate = 3 and neq increased (i1 to i2)  ',
     1   50, 5, 0, 2, n, neq(1), 0, 0.0d0, 0.0d0)
      go to 700
 606  call xerrwv('lsode--  itol (=i1) illegal   ',
     1   30, 6, 0, 1, itol, 0, 0, 0.0d0, 0.0d0)
      go to 700
 607  call xerrwv('lsode--  iopt (=i1) illegal   ',
     1   30, 7, 0, 1, iopt, 0, 0, 0.0d0, 0.0d0)
      go to 700
 608  call xerrwv('lsode--  mf (=i1) illegal     ',
     1   30, 8, 0, 1, mf, 0, 0, 0.0d0, 0.0d0)
      go to 700
 609  call xerrwv('lsode--  ml (=i1) illegal.. .lt.0 or .ge.neq (=i2)',
     1   50, 9, 0, 2, ml, neq(1), 0, 0.0d0, 0.0d0)
      go to 700
 610  call xerrwv('lsode--  mu (=i1) illegal.. .lt.0 or .ge.neq (=i2)',
     1   50, 10, 0, 2, mu, neq(1), 0, 0.0d0, 0.0d0)
      go to 700
 611  call xerrwv('lsode--  maxord (=i1) .lt. 0  ',
     1   30, 11, 0, 1, maxord, 0, 0, 0.0d0, 0.0d0)
      go to 700
 612  call xerrwv('lsode--  mxstep (=i1) .lt. 0  ',
     1   30, 12, 0, 1, mxstep, 0, 0, 0.0d0, 0.0d0)
      go to 700
 613  call xerrwv('lsode--  mxhnil (=i1) .lt. 0  ',
     1   30, 13, 0, 1, mxhnil, 0, 0, 0.0d0, 0.0d0)
      go to 700
 614  call xerrwv('lsode--  tout (=r1) behind t (=r2)      ',
     1   40, 14, 0, 0, 0, 0, 2, tout, t)
      call xerrwv('      integration direction is given by h0 (=r1)  ',
     1   50, 14, 0, 0, 0, 0, 1, h0, 0.0d0)
      go to 700
 615  call xerrwv('lsode--  hmax (=r1) .lt. 0.0  ',
     1   30, 15, 0, 0, 0, 0, 1, hmax, 0.0d0)
      go to 700
 616  call xerrwv('lsode--  hmin (=r1) .lt. 0.0  ',
     1   30, 16, 0, 0, 0, 0, 1, hmin, 0.0d0)
      go to 700
 617  call xerrwv(
     1  'lsode--  rwork length needed, lenrw (=i1), exceeds lrw (=i2)',
     1   60, 17, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
      go to 700
 618  call xerrwv(
     1  'lsode--  iwork length needed, leniw (=i1), exceeds liw (=i2)',
     1   60, 18, 0, 2, leniw, liw, 0, 0.0d0, 0.0d0)
      go to 700
 619  call xerrwv('lsode--  rtol(i1) is r1 .lt. 0.0        ',
     1   40, 19, 0, 1, i, 0, 1, rtoli, 0.0d0)
      go to 700
 620  call xerrwv('lsode--  atol(i1) is r1 .lt. 0.0        ',
     1   40, 20, 0, 1, i, 0, 1, atoli, 0.0d0)
      go to 700
 621  ewti = rwork(lewt+i-1)
      call xerrwv('lsode--  ewt(i1) is r1 .le. 0.0         ',
     1   40, 21, 0, 1, i, 0, 1, ewti, 0.0d0)
      go to 700
 622  call xerrwv(
     1  'lsode--  tout (=r1) too close to t(=r2) to start integration',
     1   60, 22, 0, 0, 0, 0, 2, tout, t)
      go to 700
 623  call xerrwv(
     1  'lsode--  itask = i1 and tout (=r1) behind tcur - hu (= r2)  ',
     1   60, 23, 0, 1, itask, 0, 2, tout, tp)
      go to 700
 624  call xerrwv(
     1  'lsode--  itask = 4 or 5 and tcrit (=r1) behind tcur (=r2)   ',
     1   60, 24, 0, 0, 0, 0, 2, tcrit, tn)
      go to 700
 625  call xerrwv(
     1  'lsode--  itask = 4 or 5 and tcrit (=r1) behind tout (=r2)   ',
     1   60, 25, 0, 0, 0, 0, 2, tcrit, tout)
      go to 700
 626  call xerrwv('lsode--  at start of problem, too much accuracy   ',
     1   50, 26, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      call xerrwv(
     1  '      requested for precision of machine..  see tolsf (=r1) ',
     1   60, 26, 0, 0, 0, 0, 1, tolsf, 0.0d0)
      rwork(14) = tolsf
      go to 700
 627  call xerrwv('lsode--  trouble from intdy. itask = i1, tout = r1',
     1   50, 27, 0, 1, itask, 0, 1, tout, 0.0d0)
c
 700  if (illin .eq. 5) go to 710
      illin = illin + 1
      istate = -3
      return
 710  call xerrwv('lsode--  repeated occurrences of illegal input    ',
     1   50, 302, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
c
 800  call xerrwv('lsode--  run aborted.. apparent infinite loop     ',
     1   50, 303, 2, 0, 0, 0, 0, 0.0d0, 0.0d0)
      return
C----------------------- END OF SUBROUTINE LSODE -----------------------
      END


      SUBROUTINE CFODE (METH, ELCO, TESCO)
      implicit real*8 (a-h,o-z)
C
C  CFODE is called by the integrator routine to set coefficients
C  needed there.  The coefficients for the current method, as
C  given by the value of METH, are set for all orders and saved.
C  The maximum order assumed here is 12 if METH = 1 and 5 if METH = 2.
C  (A smaller value of the maximum order is also allowed.)
C  CFODE is called once at the beginning of the problem,
C  and is not called again unless and until METH is changed.
C
C  The ELCO array contains the basic method coefficients.
C  The coefficients el(i), 1 .le. i .le. nq+1, for the method of
C  order nq are stored in ELCO(i,nq).  They are given by a genetrating
C  polynomial, i.e.,
C      l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
C  For the implicit Adams methods, l(x) is given by
C      dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0.
C  For the BDF methods, l(x) is given by
C      l(x) = (x+1)*(x+2)* ... *(x+nq)/K,
C  where         K = factorial(nq)*(1 + 1/2 + ... + 1/nq).
C
C  The TESCO array contains test constants used for the
C  local error test and the selection of step size and/or order.
C  At order nq, TESCO(k,nq) is used for the selection of step
C  size at order nq - 1 if k = 1, at order nq if k = 2, and at order
C  nq + 1 if k = 3.
C

      integer meth
      integer i, ib, nq, nqm1, nqp1
      real*8 elco, tesco
      real*8 agamq, fnq, fnqm1, pc, pint, ragq,
     1   rqfac, rq1fac, tsign, xpin
      dimension elco(13,12), tesco(3,12)
      dimension pc(12)
c
C***FIRST EXECUTABLE STATEMENT  CFODE
      go to (100, 200), meth
c
 100  elco(1,1) = 1.0e0
      elco(2,1) = 1.0e0
      tesco(1,1) = 0.0d0
      tesco(2,1) = 2.0e0
      tesco(1,2) = 1.0e0
      tesco(3,12) = 0.0d0
      pc(1) = 1.0e0
      rqfac = 1.0e0
      do 140 nq = 2,12
c-----------------------------------------------------------------------
c the pc array will contain the coefficients of the polynomial
c     p(x) = (x+1)*(x+2)*...*(x+nq-1).
c initially, p(x) = 1.
c-----------------------------------------------------------------------
        rq1fac = rqfac
        rqfac = rqfac/float(nq)
        nqm1 = nq - 1
        fnqm1 = float(nqm1)
        nqp1 = nq + 1
c form coefficients of p(x)*(x+nq-1). ----------------------------------
        pc(nq) = 0.0d0
        do 110 ib = 1,nqm1
          i = nqp1 - ib
 110      pc(i) = pc(i-1) + fnqm1*pc(i)
        pc(1) = fnqm1*pc(1)
c compute integral, -1 to 0, of p(x) and x*p(x). -----------------------
        pint = pc(1)
        xpin = pc(1)/2.0e0
        tsign = 1.0e0
        do 120 i = 2,nq
          tsign = -tsign
          pint = pint + tsign*pc(i)/float(i)
 120      xpin = xpin + tsign*pc(i)/float(i+1)
c store coefficients in elco and tesco. --------------------------------
        elco(1,nq) = pint*rq1fac
        elco(2,nq) = 1.0e0
        do 130 i = 2,nq
 130      elco(i+1,nq) = rq1fac*pc(i)/float(i)
        agamq = rqfac*xpin
        ragq = 1.0e0/agamq
        tesco(2,nq) = ragq
        if (nq .lt. 12) tesco(1,nqp1) = ragq*rqfac/float(nqp1)
        tesco(3,nqm1) = ragq
 140    continue
      return
c
 200  pc(1) = 1.0e0
      rq1fac = 1.0e0
      do 230 nq = 1,5
c-----------------------------------------------------------------------
c the pc array will contain the coefficients of the polynomial
c     p(x) = (x+1)*(x+2)*...*(x+nq).
c initially, p(x) = 1.
c-----------------------------------------------------------------------
        fnq = float(nq)
        nqp1 = nq + 1
c form coefficients of p(x)*(x+nq). ------------------------------------
        pc(nqp1) = 0.0d0
        do 210 ib = 1,nq
          i = nq + 2 - ib
 210      pc(i) = pc(i-1) + fnq*pc(i)
        pc(1) = fnq*pc(1)
c store coefficients in elco and tesco. --------------------------------
        do 220 i = 1,nqp1
 220      elco(i,nq) = pc(i)/pc(2)
        elco(2,nq) = 1.0e0
        tesco(1,nq) = rq1fac
        tesco(2,nq) = float(nqp1)/elco(1,nq)
        tesco(3,nq) = float(nq+2)/elco(1,nq)
        rq1fac = rq1fac/fnq
 230    continue
      return
C----------------------- END OF SUBROUTINE CFODE -----------------------
      END
 
      SUBROUTINE EWSET (N, ITOL, RTOL, ATOL, YCUR, EWT)
      implicit real*8 (a-h,o-z)
C
C  This subroutine sets the error weight vector EWT according to
C      EWT(i) = RTOL(i)*ABS(YCUR(i)) + ATOL(i),  i = 1,...,N,
C  with the subscript on RTOL and/or ATOL possibly replaced by 1 above,
C  depending on the value of ITOL.
C
C
      integer n, itol
      integer i
      real*8 rtol, atol, ycur, ewt

      dimension rtol(1), atol(1), ycur(n), ewt(n)
c
C***FIRST EXECUTABLE STATEMENT  EWSET
      go to (10, 20, 30, 40), itol
 10   continue
      do 15 i = 1,n
 15     ewt(i) = rtol(1)*abs(ycur(i)) + atol(1)
      return
 20   continue
      do 25 i = 1,n
 25     ewt(i) = rtol(1)*abs(ycur(i)) + atol(i)
      return
 30   continue
      do 35 i = 1,n
 35     ewt(i) = rtol(i)*abs(ycur(i)) + atol(1)
      return
 40   continue
      do 45 i = 1,n
 45     ewt(i) = rtol(i)*abs(ycur(i)) + atol(i)
      return
C----------------------- END OF SUBROUTINE EWSET -----------------------
      END
 
      SUBROUTINE INTDY (T, K, YH, NYH, DKY, IFLAG)
      implicit real*8 (a-h,o-z)
C
C  INTDY computes interpolated values of the K-th derivative of the
C  dependent variable vector y, and stores it in DKY.  This routine
C  is called within the package with K = 0 and T = TOUT, but may
C  also be called by the user for any K up to the current order.
C  (See detailed instructions in the usage documentation.)
C
C  The computed values in DKY are gotten by interpolation using the
C  Nordsieck history array YH.  This array corresponds uniquely to a
C  vector-valued polynomial of degree NQCUR or less, and DKY is set
C  to the K-th derivative of this polynomial at T.
C  The formula for DKY is:
C               q
C   DKY(i)  =  sum  c(j,K) * (T - tn)**(j-K) * h**(-j) * YH(i,j+1)
C              j=K
C  where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, tn = TCUR, h = HCUR.
C  The quantities  nq = NQCUR, l = nq+1, N = NEQ, tn, and h are
C  communicated by COMMON.  The above sum is done in reverse order.
C  IFLAG is returned negative if either K or T is out of bounds.
C
      integer k, nyh, iflag
      integer iownd, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, ic, j, jb, jb2, jj, jj1, jp1
      real*8 t, yh, dky
      real*8 rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real*8 c, r, s, tp
      dimension yh(nyh,1), dky(1)
      common /ls0001R/ rowns(209),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      common /ls0001I/
     3   iownd(14), iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
C
C***FIRST EXECUTABLE STATEMENT  INTDY
      iflag = 0
      if (k .lt. 0 .or. k .gt. nq) go to 80
      tp = tn - hu -  100.0d0*uround*(tn + hu)
      if ((t-tp)*(t-tn) .gt. 0.0d0) go to 90
c
      s = (t - tn)/h
      ic = 1
      if (k .eq. 0) go to 15
      jj1 = l - k
      do 10 jj = jj1,nq
 10     ic = ic*jj
 15   c = float(ic)
      do 20 i = 1,n
 20     dky(i) = c*yh(i,l)
      if (k .eq. nq) go to 55
      jb2 = nq - k
      do 50 jb = 1,jb2
        j = nq - jb
        jp1 = j + 1
        ic = 1
        if (k .eq. 0) go to 35
        jj1 = jp1 - k
        do 30 jj = jj1,j
 30       ic = ic*jj
 35     c = float(ic)
        do 40 i = 1,n
 40       dky(i) = c*yh(i,jp1) + s*dky(i)
 50     continue
      if (k .eq. 0) return
 55   r = h**(-k)
      do 60 i = 1,n
 60     dky(i) = r*dky(i)
      return
c
 80   call xerrwv('ntdy--  k (=i1) illegal      ',
     1   30, 51, 0, 1, k, 0, 0, 0.0d0, 0.0d0)
      iflag = -1
      return
 90   call xerrwv('intdy--  t (=r1) illegal      ',
     1   30, 52, 0, 0, 0, 0, 1, t, 0.0d0)
      call xerrwv(
     1  '      t not in interval tcur - hu (= r1) to tcur (=r2)      ',
     1   60, 52, 0, 0, 0, 0, 2, tp, tn)
      iflag = -2
      return
C----------------------- END OF SUBROUTINE INTDY -----------------------
      END
 
      SUBROUTINE PREPJ(NEQ, Y, YH, NYH, EWT, FTEM, SAVF, WM, IWM, F,JAC)
      implicit real*8 (a-h,o-z)
C
C  PREPJ is called by STODE to compute and process the matrix
C  P = I - h*el(1)*J , where J is an approximation to the Jacobian.
C  Here J is computed by the user-supplied routine JAC if
C  MITER = 1 or 4, or by finite differencing if MITER = 2, 3, or 5.
C  If MITER = 3, a diagonal approximation to J is used.
C  J is stored in WM and replaced by P.  If MITER .ne. 3, P is then
C  subjected to LU decomposition in preparation for later solution
C  of linear systems with P as coefficient matrix.  This is done
C  by SGEFA if MITER = 1 or 2, and by SGBFA if MITER = 4 or 5.
C
C  In addition to variables described in STODE and LSODE prologues,
C  communication with PREPJ uses the following:
C  Y     = array containing predicted values on entry.
C  FTEM  = work array of length N (ACOR in STODE).
C  SAVF  = array containing f evaluated at predicted y.
C  WM    = real work space for matrices.  On output it contains the
C          inverse diagonal matrix if MITER = 3 and the LU decomposition
C          of P if MITER is 1, 2 , 4, or 5.
C          Storage of matrix elements starts at WM(3).
C          WM also contains the following matrix-related data:
C          WM(1) = SQRT(UROUND), used in numerical Jacobian increments.
C          WM(2) = h*EL0, saved for later use if MITER = 3.
C  IWM   = integer work space containing pivot information, starting at
C          IWM(21), if MITER is 1, 2, 4, or 5.  IWM also contains band
C          parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
C  EL0   = el(1) (input).
C  IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if
C          P matrix found to be singular.
C  JCUR  = output flag = 1 to indicate that the Jacobian matrix
C          (or approximation) is now current.
C  This routine also uses the COMMON variables EL0, H, TN, UROUND,
C  MITER, N, NFE, and NJE.
C
      external f, jac
      integer neq, nyh, iwm
      integer iownd, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, i1, i2, ier, ii, j, j1, jj, lenp,
     1   mba, mband, meb1, meband, ml, ml3, mu, np1
      real*8 y, yh, ewt, ftem, savf, wm
      real*8 rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real*8 con, di, fac, hl0, r, r0, srur, yi, yj, yjj,
     1   vnorm
      dimension neq(1), y(1), yh(nyh,1), ewt(1), ftem(1), savf(1),
     1   wm(1), iwm(1)
      common /ls0001R/ rowns(209),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      common /ls0001I/
     3   iownd(14), iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
C
C***FIRST EXECUTABLE STATEMENT  PREPJ
      ii2=2
      ii3=3  
      ii21=21
      nje = nje + 1
      ierpj = 0
      jcur = 1
      hl0 = h*el0
      go to (100, 200, 300, 400, 500), miter
c if miter = 1, call jac and multiply by scalar. -----------------------
 100  lenp = n*n
      do 110 i = 1,lenp
 110    wm(i+2) = 0.0d0
      call jac (neq, tn, y, 0, 0, wm(ii3), n)
      con = -hl0
      do 120 i = 1,lenp
 120    wm(i+2) = wm(i+2)*con
      go to 240
c if miter = 2, make n calls to f to approximate j. --------------------
 200  fac = vnorm (n, savf, ewt)
      r0 = 1000.0d0*abs(h)*uround*float(n)*fac
      if (r0 .eq. 0.0d0) r0 = 1.0e0
      srur = wm(1)
      j1 = 2
      do 230 j = 1,n
        yj = y(j)
        r = max(srur*abs(yj),r0/ewt(j))
        y(j) = y(j) + r
        fac = -hl0/r
        call f (neq, tn, y, ftem)
        do 220 i = 1,n
 220      wm(i+j1) = (ftem(i) - savf(i))*fac
        y(j) = yj
        j1 = j1 + n
 230    continue
      nfe = nfe + n
c add identity matrix. -------------------------------------------------
 240  j = 3
      np1 = n + 1
      do 250 i = 1,n
        wm(j) = wm(j) + 1.0e0
 250    j = j + np1
c do lu decomposition on p. --------------------------------------------
      call sgefa (wm(ii3), n, n, iwm(ii21), ier)
      if (ier .ne. 0) ierpj = 1
      return
c if miter = 3, construct a diagonal approximation to j and p. ---------
 300  wm(ii2) = hl0
      r = el0*0.1e0
      do 310 i = 1,n
 310    y(i) = y(i) + r*(h*savf(i) - yh(i,ii2))
      ii3=3        
      call f (neq, tn, y, wm(ii3))
      nfe = nfe + 1
      do 320 i = 1,n
        r0 = h*savf(i) - yh(i,ii2)
        di = 0.1e0*r0 - h*(wm(i+2) - savf(i))
        wm(i+2) = 1.0e0
        if (abs(r0) .lt. uround/ewt(i)) go to 320
        if (abs(di) .eq. 0.0d0) go to 330
        wm(i+2) = 0.1e0*r0/di
 320    continue
      return
 330  ierpj = 1
      return
c if miter = 4, call jac and multiply by scalar. -----------------------
 400  ml = iwm(1)
      mu = iwm(ii2)
      ml3 = ml + 3
      mband = ml + mu + 1
      meband = mband + ml
      lenp = meband*n
      do 410 i = 1,lenp
 410    wm(i+2) = 0.0d0
      call jac (neq, tn, y, ml, mu, wm(ml3), meband)
      con = -hl0
      do 420 i = 1,lenp
 420    wm(i+2) = wm(i+2)*con
      go to 570
c if miter = 5, make mband calls to f to approximate j. ----------------
 500  ml = iwm(1)
      mu = iwm(ii2)
      mband = ml + mu + 1
      mba = min0(mband,n)
      meband = mband + ml
      meb1 = meband - 1
      srur = wm(1)
      fac = vnorm (n, savf, ewt)
      r0 = 1000.0d0*abs(h)*uround*float(n)*fac
      if (r0 .eq. 0.0d0) r0 = 1.0e0
      do 560 j = 1,mba
        do 530 i = j,n,mband
          yi = y(i)
          r = max(srur*abs(yi),r0/ewt(i))
 530      y(i) = y(i) + r
        call f (neq, tn, y, ftem)
        do 550 jj = j,n,mband
          y(jj) = yh(jj,1)
          yjj = y(jj)
          r = max(srur*abs(yjj),r0/ewt(jj))
          fac = -hl0/r
          i1 = max0(jj-mu,1)
          i2 = min0(jj+ml,n)
          ii = jj*meb1 - ml + 2
          do 540 i = i1,i2
 540        wm(ii+i) = (ftem(i) - savf(i))*fac
 550      continue
 560    continue
      nfe = nfe + mba
c add identity matrix. -------------------------------------------------
 570  ii = mband + 2
      do 580 i = 1,n
        wm(ii) = wm(ii) + 1.0e0
 580    ii = ii + meband
c do lu decomposition of p. --------------------------------------------
      call sgbfa (wm(ii3), meband, n, ml, mu, iwm(ii21), ier)
      if (ier .ne. 0) ierpj = 1
      return
C----------------------- END OF SUBROUTINE PREPJ -----------------------
      END
 
      SUBROUTINE SOLSY (WM, IWM, X, TEM)
      implicit real*8 (a-h,o-z)
C
C  This routine manages the solution of the linear system arising from
C  a chord iteration.  It is called if MITER .ne. 0.
C  If MITER is 1 or 2, it calls SGESL to accomplish this.
C  If MITER = 3 it updates the coefficient h*EL0 in the diagonal
C  matrix, and then computes the solution.
C  If MITER is 4 or 5, it calls SGBSL.
C  Communication with SOLSY uses the following variables:
C  WM    = real work space containing the inverse diagonal matrix if
C          MITER = 3 and the LU decomposition of the matrix otherwise.
C          Storage of matrix elements starts at WM(3).
C          WM also contains the following matrix-related data:
C          WM(1) = SQRT(UROUND) (not used here),
C          WM(2) = HL0, the previous value of h*EL0, used if MITER = 3.
C  IWM   = integer work space containing pivot information, starting at
C          IWM(21), if MITER is 1, 2, 4, or 5.  IWM also contains band
C          parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
C  X     = the right-hand side vector on input, and the solution vector
C          on output, of length N.
C  TEM   = vector of work space of length N, not used in this version.
C  IERSL = output flag (in COMMON).  IERSL = 0 if no trouble occurred.
C          IERSL = 1 if a singular matrix arose with MITER = 3.
C  This routine also uses the COMMON variables EL0, H, MITER, and N.
C
      integer iwm
      integer iownd, iowns,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, meband, ml, mu
      real*8 wm, x, tem
      real*8 rowns,
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real*8 di, hl0, phl0, r
      dimension wm(1), iwm(1), x(1), tem(1)
      common /ls0001R/ rowns(209),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      common /ls0001I/
     3   iownd(14), iowns(6),
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
C
C***FIRST EXECUTABLE STATEMENT  SOLSY
      ii2=2
      ii3=3
      ii21=21
      iersl = 0
      go to (100, 100, 300, 400, 400), miter
 100  call sgesl (wm(ii3), n, n, iwm(ii21), x, 0)
      return
c
 300  phl0 = wm(ii2)
      hl0 = h*el0
      wm(ii2) = hl0
      if (hl0 .eq. phl0) go to 330
      r = hl0/phl0
      do 320 i = 1,n
        di = 1.0e0 - r*(1.0e0 - 1.0e0/wm(i+2))
        if (abs(di) .eq. 0.0d0) go to 390
 320    wm(i+2) = 1.0e0/di
 330  do 340 i = 1,n
 340    x(i) = wm(i+2)*x(i)
      return
 390  iersl = 1
      return
c
 400  ml = iwm(1)
      mu = iwm(ii2)
      meband = 2*ml + mu + 1
      call sgbsl (wm(ii3), meband, n, ml, mu, iwm(ii21), x, 0)
      return
C----------------------- END OF SUBROUTINE SOLSY -----------------------
      END
 
      SUBROUTINE SRCOM (RSAV, ISAV, JOB)
      implicit real*8 (a-h,o-z)
C
C  This routine saves or restores (depending on JOB) the contents of
C  the COMMON blocks LS0001 and EH0001, which are used internally
C  by one or more ODEPACK solvers.
C
C  RSAV = real array of length 218 or more.
C  ISAV = integer array of length 41 or more.
C  JOB  = flag indicating to save or restore the COMMON blocks:
C         JOB  = 1 if COMMON is to be saved (written to RSAV/ISAV)
C         JOB  = 2 if COMMON is to be restored (read from RSAV/ISAV)
C         A call with JOB = 2 presumes a prior call with JOB = 1.
C
      integer isav, job
      integer ieh, ils
      integer i, lenils, lenrls
      real*8 rsav,   rls
      dimension rsav(1), isav(1)
      common /ls0001R/ rls(218)
      common /ls0001I/ ils(39)
      common /eh0001/ ieh(3)
      data lenrls/218/, lenils/39/
c
C***FIRST EXECUTABLE STATEMENT  SRCOM
      if (job .eq. 2) go to 100
c
      do 10 i = 1,lenrls
 10     rsav(i) = rls(i)
      do 20 i = 1,lenils
 20     isav(i) = ils(i)
      isav(lenils+1) = ieh(1)
      isav(lenils+2) = ieh(2)
      return
c
 100  continue
      do 110 i = 1,lenrls
 110     rls(i) = rsav(i)
      do 120 i = 1,lenils
 120     ils(i) = isav(i)
      ieh(1) = isav(lenils+1)
      ieh(2) = isav(lenils+2)
      return
C----------------------- END OF SUBROUTINE SRCOM -----------------------
      END

      SUBROUTINE STODE (NEQ, Y, YH, NYH, YH1, EWT, SAVF, ACOR,
     1   WM, IWM, F, JAC, PJAC, SLVS)
      implicit real*8 (a-h,o-z)
C
C  STODE performs one step of the integration of an initial value
C  problem for a system of ordinary differential equations.
C  Note:  STODE is independent of the value of the iteration method
C  indicator MITER, when this is .ne. 0, and hence is independent
C  of the type of chord method used, or the Jacobian structure.
C  Communication with STODE is done with the following variables:
C
C  NEQ    = integer array containing problem size in NEQ(1), and
C           passed as the NEQ argument in all calls to F and JAC.
C  Y      = an array of length .ge. N used as the Y argument in
C           all calls to F and JAC.
C  YH     = an NYH by LMAX array containing the dependent variables
C           and their approximate scaled derivatives, where
C           LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
C           j-th derivative of y(i), scaled by h**j/factorial(j)
C           (j = 0,1,...,NQ).  on entry for the first step, the first
C           two columns of YH must be set from the initial values.
C  NYH    = a constant integer .ge. N, the first dimension of YH.
C  YH1    = a one-dimensional array occupying the same space as YH.
C  EWT    = an array of length N containing multiplicative weights
C           for local error measurements.  Local errors in Y(i) are
C           compared to 1.0/EWT(i) in various error tests.
C  SAVF   = an array of working storage, of length N.
C           Also used for input of YH(*,MAXORD+2) when JSTART = -1
C           and MAXORD .lt. the current order NQ.
C  ACOR   = a work array of length N, used for the accumulated
C           corrections.  On a successful return, ACOR(i) contains
C           the estimated one-step local error in Y(i).
C  WM,IWM = real and integer work arrays associated with matrix
C           operations in chord iteration (MITER .ne. 0).
C  PJAC   = name of routine to evaluate and preprocess Jacobian matrix
C           and P = I - h*el0*JAC, if a chord method is being used.
C  SLVS   = name of routine to solve linear system in chord iteration.
C  CCMAX  = maximum relative change in h*el0 before PJAC is called.
C  H      = the step size to be attempted on the next step.
C           H is altered by the error control algorithm during the
C           problem.  H can be either positive or negative, but its
C           sign must remain constant throughout the problem.
C  HMIN   = the minimum absolute value of the step size h to be used.
C  HMXI   = inverse of the maximum absolute value of h to be used.
C           HMXI = 0.0 is allowed and corresponds to an infinite hmax.
C           HMIN and HMXI may be changed at any time, but will not
C           take effect until the next change of h is considered.
C  TN     = the independent variable. TN is updated on each step taken.
C  JSTART = an integer used for input only, with the following
C           values and meanings:
C                0  perform the first step.
C            .gt.0  take a new step continuing from the last.
C               -1  take the next step with a new value of H, MAXORD,
C                     N, METH, MITER, and/or matrix parameters.
C               -2  take the next step with a new value of H,
C                     but with other inputs unchanged.
C           On return, JSTART is set to 1 to facilitate continuation.
C  KFLAG  = a completion code with the following meanings:
C                0  the step was succesful.
C               -1  the requested error could not be achieved.
C               -2  corrector convergence could not be achieved.
C               -3  fatal error in PJAC or SLVS.
C           A return with KFLAG = -1 or -2 means either
C           abs(H) = HMIN or 10 consecutive failures occurred.
C           On a return with KFLAG negative, the values of TN and
C           the YH array are as of the beginning of the last
C           step, and H is the last step size attempted.
C  MAXORD = the maximum order of integration method to be allowed.
C  MAXCOR = the maximum number of corrector iterations allowed.
C  MSBP   = maximum number of steps between PJAC calls (MITER .gt. 0).
C  MXNCF  = maximum number of convergence failures allowed.
C  METH/MITER = the method flags.  See description in driver.
C  N      = the number of first-order differential equations.
C  The values of CCMAX, H, HMIN, HMXI, TN, JSTART, KFLAG, MAXORD,
C  MAXCOR, MSBP, MXNCF, METH, MITER, and N are communicated via COMMON.
C
      external f, jac, pjac, slvs
      integer neq, nyh, iwm
      integer iownd, ialth, ipup, lmax, meo, nqnyh, nslp,
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, i1, iredo, iret, j, jb, m, ncf, newq
      real*8 y, yh, yh1, ewt, savf, acor, wm
      real*8 conit, crate, el, elco, hold, rmax, tesco,
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real*8 dcon, ddn, del, delp, dsm, dup, exdn, exsm, exup,
     1   r, rh, rhdn, rhsm, rhup, told, vnorm
      dimension neq(1), y(1), yh(nyh,1), yh1(1), ewt(1), savf(1),
     1   acor(1), wm(1), iwm(1)
      common /ls0001R/ conit, crate, el(13), elco(13,12),
     1   hold, rmax, tesco(3,12),
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      common /ls0001I/ iownd(14),
     3   ialth, ipup, lmax, meo, nqnyh, nslp,
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
C
C***FIRST EXECUTABLE STATEMENT  STODE
      kflag = 0
      told = tn
      ncf = 0
      ierpj = 0
      iersl = 0
      jcur = 0
      icf = 0
      delp = 0.0d0
      if (jstart .gt. 0) go to 200
      if (jstart .eq. -1) go to 100
      if (jstart .eq. -2) go to 160
c-----------------------------------------------------------------------
c on the first call, the order is set to 1, and other variables are
c initialized.  rmax is the maximum ratio by which h can be increased
c in a single step.  it is initially 1.e4 to compensate for the small
c initial h, but then is normally equal to 10.  if a failure
c occurs (in corrector convergence or error test), rmax is set at 2
c for the next increase.
c-----------------------------------------------------------------------
      lmax = maxord + 1
      nq = 1
      l = 2
      ialth = 2
      rmax = 10000.0d0
      rc = 0.0d0
      el0 = 1.0e0
      crate = 0.7e0
      hold = h
      meo = meth
      nslp = 0
      ipup = miter
      iret = 3
      go to 140
c-----------------------------------------------------------------------
c the following block handles preliminaries needed when jstart = -1.
c ipup is set to miter to force a matrix update.
c if an order increase is about to be considered (ialth = 1),
c ialth is reset to 2 to postpone consideration one more step.
c if the caller has changed meth, cfode is called to reset
c the coefficients of the method.
c if the caller has changed maxord to a value less than the current
c order nq, nq is reduced to maxord, and a new h chosen accordingly.
c if h is to be changed, yh must be rescaled.
c if h or meth is being changed, ialth is reset to l = nq + 1
c to prevent further changes in h for that many steps.
c-----------------------------------------------------------------------
 100  ipup = miter
      lmax = maxord + 1
      if (ialth .eq. 1) ialth = 2
      if (meth .eq. meo) go to 110
      call cfode (meth, elco, tesco)
      meo = meth
      if (nq .gt. maxord) go to 120
      ialth = l
      iret = 1
      go to 150
 110  if (nq .le. maxord) go to 160
 120  nq = maxord
      l = lmax
      do 125 i = 1,l
 125    el(i) = elco(i,nq)
      nqnyh = nq*nyh
      rc = rc*el(1)/el0
      el0 = el(1)
      conit = 0.5e0/float(nq+2)
      ddn = vnorm (n, savf, ewt)/tesco(1,l)
      exdn = 1.0e0/float(l)
      rhdn = 1.0e0/(1.3e0*ddn**exdn + 0.0000013e0)
      rh = min(rhdn,1.0e0)
      iredo = 3
      if (h .eq. hold) go to 170
      rh = min(rh,abs(h/hold))
      h = hold
      go to 175
c-----------------------------------------------------------------------
c cfode is called to get all the integration coefficients for the
c current meth.  then the el vector and related constants are reset
c whenever the order nq is changed, or at the start of the problem.
c-----------------------------------------------------------------------
 140  call cfode (meth, elco, tesco)
 150  do 155 i = 1,l
 155    el(i) = elco(i,nq)
      nqnyh = nq*nyh
      rc = rc*el(1)/el0
      el0 = el(1)
      conit = 0.5e0/float(nq+2)
      go to (160, 170, 200), iret
c-----------------------------------------------------------------------
c if h is being changed, the h ratio rh is checked against
c rmax, hmin, and hmxi, and the yh array rescaled.  ialth is set to
c l = nq + 1 to prevent a change of h for that many steps, unless
c forced by a convergence or error test failure.
c-----------------------------------------------------------------------
 160  if (h .eq. hold) go to 200
      rh = h/hold
      h = hold
      iredo = 3
      go to 175
 170  rh = max(rh,hmin/abs(h))
 175  rh = min(rh,rmax)
      rh = rh/max(1.0e0,abs(h)*hmxi*rh)
      r = 1.0e0
      do 180 j = 2,l
        r = r*rh
        do 180 i = 1,n
 180      yh(i,j) = yh(i,j)*r
      h = h*rh
      rc = rc*rh
      ialth = l
      if (iredo .eq. 0) go to 690
c-----------------------------------------------------------------------
c this section computes the predicted values by effectively
c multiplying the yh array by the pascal triangle matrix.
c rc is the ratio of new to old values of the coefficient  h*el(1).
c when rc differs from 1 by more than ccmax, ipup is set to miter
c to force pjac to be called, if a jacobian is involved.
c in any case, pjac is called at least every msbp steps.
c-----------------------------------------------------------------------
 200  if (abs(rc-1.0e0) .gt. ccmax) ipup = miter
      if (nst .ge. nslp+msbp) ipup = miter
      tn = tn + h
      i1 = nqnyh + 1
      do 215 jb = 1,nq
        i1 = i1 - nyh
cdir$ ivdep
        do 210 i = i1,nqnyh
 210      yh1(i) = yh1(i) + yh1(i+nyh)
 215    continue
c-----------------------------------------------------------------------
c up to maxcor corrector iterations are taken.  a convergence test is
c made on the r.m.s. norm of each correction, weighted by the error
c weight vector ewt.  the sum of the corrections is accumulated in the
c vector acor(i).  the yh array is not altered in the corrector loop.
c-----------------------------------------------------------------------
 220  m = 0
      do 230 i = 1,n
 230    y(i) = yh(i,1)
      call f (neq, tn, y, savf)
      nfe = nfe + 1
      if (ipup .le. 0) go to 250
c-----------------------------------------------------------------------
c if indicated, the matrix p = i - h*el(1)*j is reevaluated and
c preprocessed before starting the corrector iteration.  ipup is set
c to 0 as an indicator that this has been done.
c-----------------------------------------------------------------------
      call pjac (neq, y, yh, nyh, ewt, acor, savf, wm, iwm, f, jac)
      ipup = 0
      rc = 1.0e0
      nslp = nst
      crate = 0.7e0
      if (ierpj .ne. 0) go to 430
 250  do 260 i = 1,n
 260    acor(i) = 0.0d0
 270  if (miter .ne. 0) go to 350
c-----------------------------------------------------------------------
c in the case of functional iteration, update y directly from
c the result of the last function evaluation.
c-----------------------------------------------------------------------
      ii2=2
      do 290 i = 1,n
        savf(i) = h*savf(i) - yh(i,ii2)
 290    y(i) = savf(i) - acor(i)
      del = vnorm (n, y, ewt)
      do 300 i = 1,n
        y(i) = yh(i,1) + el(1)*savf(i)
 300    acor(i) = savf(i)
      go to 400
c-----------------------------------------------------------------------
c in the case of the chord method, compute the corrector error,
c and solve the linear system with that as right-hand side and
c p as coefficient matrix.
c-----------------------------------------------------------------------
 350  ii2=2
      do 360 i = 1,n
 360    y(i) = h*savf(i) - (yh(i,ii2) + acor(i))
      call slvs (wm, iwm, y, savf)
      if (iersl .lt. 0) go to 430
      if (iersl .gt. 0) go to 410
      del = vnorm (n, y, ewt)
      do 380 i = 1,n
        acor(i) = acor(i) + y(i)
 380    y(i) = yh(i,1) + el(1)*acor(i)
c-----------------------------------------------------------------------
c test for convergence.  if m.gt.0, an estimate of the convergence
c rate constant is stored in crate, and this is used in the test.
c-----------------------------------------------------------------------
 400  if (m .ne. 0) crate = max(0.2e0*crate,del/delp)
      dcon = del*min(1.0e0,1.5e0*crate)/(tesco(2,nq)*conit)
      if (dcon .le. 1.0e0) go to 450
      m = m + 1
      if (m .eq. maxcor) go to 410
      if (m .ge. 2 .and. del .gt. 2.0e0*delp) go to 410
      delp = del
      call f (neq, tn, y, savf)
      nfe = nfe + 1
      go to 270
c-----------------------------------------------------------------------
c the corrector iteration failed to converge.
c if miter .ne. 0 and the jacobian is out of date, pjac is called for
c the next try.  otherwise the yh array is retracted to its values
c before prediction, and h is reduced, if possible.  if h cannot be
c reduced or mxncf failures have occurred, exit with kflag = -2.
c-----------------------------------------------------------------------
 410  if (miter .eq. 0 .or. jcur .eq. 1) go to 430
      icf = 1
      ipup = miter
      go to 220
 430  icf = 2
      ncf = ncf + 1
      rmax = 2.0e0
      tn = told
      i1 = nqnyh + 1
      do 445 jb = 1,nq
        i1 = i1 - nyh
cdir$ ivdep
        do 440 i = i1,nqnyh
 440      yh1(i) = yh1(i) - yh1(i+nyh)
 445    continue
      if (ierpj .lt. 0 .or. iersl .lt. 0) go to 680
      if (abs(h) .le. hmin*1.00001e0) go to 670
      if (ncf .eq. mxncf) go to 670
      rh = 0.25e0
      ipup = miter
      iredo = 1
      go to 170
c-----------------------------------------------------------------------
c the corrector has converged.  jcur is set to 0
c to signal that the jacobian involved may need updating later.
c the local error test is made and control passes to statement 500
c if it fails.
c-----------------------------------------------------------------------
 450  jcur = 0
      if (m .eq. 0) dsm = del/tesco(2,nq)
      if (m .gt. 0) dsm = vnorm (n, acor, ewt)/tesco(2,nq)
      if (dsm .gt. 1.0e0) go to 500
c-----------------------------------------------------------------------
c after a successful step, update the yh array.
c consider changing h if ialth = 1.  otherwise decrease ialth by 1.
c if ialth is then 1 and nq .lt. maxord, then acor is saved for
c use in a possible order increase on the next step.
c if a change in h is considered, an increase or decrease in order
c by one is considered also.  a change in h is made only if it is by a
c factor of at least 1.1.  if not, ialth is set to 3 to prevent
c testing for that many steps.
c-----------------------------------------------------------------------
      kflag = 0
      iredo = 0
      nst = nst + 1
      hu = h
      nqu = nq
      do 470 j = 1,l
        do 470 i = 1,n
 470      yh(i,j) = yh(i,j) + el(j)*acor(i)
      ialth = ialth - 1
      if (ialth .eq. 0) go to 520
      if (ialth .gt. 1) go to 700
      if (l .eq. lmax) go to 700
      do 490 i = 1,n
 490    yh(i,lmax) = acor(i)
      go to 700
c-----------------------------------------------------------------------
c the error test failed.  kflag keeps track of multiple failures.
c restore tn and the yh array to their previous values, and prepare
c to try the step again.  compute the optimum step size for this or
c one lower order.  after 2 or more failures, h is forced to decrease
c by a factor of 0.2 or less.
c-----------------------------------------------------------------------
 500  kflag = kflag - 1
      tn = told
      i1 = nqnyh + 1
      do 515 jb = 1,nq
        i1 = i1 - nyh
cdir$ ivdep
        do 510 i = i1,nqnyh
 510      yh1(i) = yh1(i) - yh1(i+nyh)
 515    continue
      rmax = 2.0e0
      if (abs(h) .le. hmin*1.00001e0) go to 660
      if (kflag .le. -3) go to 640
      iredo = 2
      rhup = 0.0d0
      go to 540
c-----------------------------------------------------------------------
c regardless of the success or failure of the step, factors
c rhdn, rhsm, and rhup are computed, by which h could be multiplied
c at order nq - 1, order nq, or order nq + 1, respectively.
c in the case of failure, rhup = 0.0 to avoid an order increase.
c the largest of these is determined and the new order chosen
c accordingly.  if the order is to be increased, we compute one
c additional scaled derivative.
c-----------------------------------------------------------------------
 520  rhup = 0.0d0
      if (l .eq. lmax) go to 540
      do 530 i = 1,n
 530    savf(i) = acor(i) - yh(i,lmax)
      dup = vnorm (n, savf, ewt)/tesco(3,nq)
      exup = 1.0e0/float(l+1)
      rhup = 1.0e0/(1.4e0*dup**exup + 0.0000014e0)
 540  exsm = 1.0e0/float(l)
      rhsm = 1.0e0/(1.2e0*dsm**exsm + 0.0000012e0)
      rhdn = 0.0d0
      if (nq .eq. 1) go to 560
      ddn = vnorm (n, yh(1,l), ewt)/tesco(1,nq)
      exdn = 1.0e0/float(nq)
      rhdn = 1.0e0/(1.3e0*ddn**exdn + 0.0000013e0)
 560  if (rhsm .ge. rhup) go to 570
      if (rhup .gt. rhdn) go to 590
      go to 580
 570  if (rhsm .lt. rhdn) go to 580
      newq = nq
      rh = rhsm
      go to 620
 580  newq = nq - 1
      rh = rhdn
      if (kflag .lt. 0 .and. rh .gt. 1.0e0) rh = 1.0e0
      go to 620
 590  newq = l
      rh = rhup
      if (rh .lt. 1.1e0) go to 610
      r = el(l)/float(l)
      do 600 i = 1,n
 600    yh(i,newq+1) = acor(i)*r
      go to 630
 610  ialth = 3
      go to 700
 620  if ((kflag .eq. 0) .and. (rh .lt. 1.1e0)) go to 610
      if (kflag .le. -2) rh = min(rh,0.2e0)
c-----------------------------------------------------------------------
c if there is a change of order, reset nq, l, and the coefficients.
c in any case h is reset according to rh and the yh array is rescaled.
c then exit from 690 if the step was ok, or redo the step otherwise.
c-----------------------------------------------------------------------
      if (newq .eq. nq) go to 170
 630  nq = newq
      l = nq + 1
      iret = 2
      go to 150
c-----------------------------------------------------------------------
c control reaches this section if 3 or more failures have occured.
c if 10 failures have occurred, exit with kflag = -1.
c it is assumed that the derivatives that have accumulated in the
c yh array have errors of the wrong order.  hence the first
c derivative is recomputed, and the order is set to 1.  then
c h is reduced by a factor of 10, and the step is retried,
c until it succeeds or h reaches hmin.
c-----------------------------------------------------------------------
 640  if (kflag .eq. -10) go to 660
      rh = 0.1e0
      rh = max(hmin/abs(h),rh)
      h = h*rh
      do 645 i = 1,n
 645    y(i) = yh(i,1)
      call f (neq, tn, y, savf)
      nfe = nfe + 1
      ii2=2
      do 650 i = 1,n
 650    yh(i,ii2) = h*savf(i)
      ipup = miter
      ialth = 5
      if (nq .eq. 1) go to 200
      nq = 1
      l = 2
      iret = 3
      go to 150
c-----------------------------------------------------------------------
c all returns are made through this section.  h is saved in hold
c to allow the caller to change h on the next step.
c-----------------------------------------------------------------------
 660  kflag = -1
      go to 720
 670  kflag = -2
      go to 720
 680  kflag = -3
      go to 720
 690  rmax = 10.0d0
 700  r = 1.0e0/tesco(2,nqu)
      do 710 i = 1,n
 710    acor(i) = acor(i)*r
 720  hold = h
      jstart = 1
      return
C----------------------- END OF SUBROUTINE STODE -----------------------
      END
 
      REAL*8 FUNCTION VNORM (N, V, W)
C
C  This function routine computes the weighted root-mean-square norm
C  of the vector of length N contained in the array V, with weights
C  contained in the array W of length N:
C    VNORM = SQRT( (1/N) * SUM( V(i)*W(i) )**2 )
C
      integer n,   i
      real*8 v, w,sum
      dimension v(n), w(n)
C
C***FIRST EXECUTABLE STATEMENT  VNORM
      sum = 0.0d0
      do 10 i = 1,n
 10     sum = sum + (v(i)*w(i))**2
      vnorm = sqrt(sum/float(n))
      return
C----------------------- END OF FUNCTION VNORM -------------------------
      END

      real*8 function r1mach (idum)
      integer idum
c-----------------------------------------------------------------------
c this routine computes the unit roundoff of the machine in single
c precision.  this is defined as the smallest positive machine number
c u such that  1.0 + u .ne. 1.0.
c-----------------------------------------------------------------------
      real*8 u, comp
      u = 1.0
 10   u = u*0.5
      comp = 1.0 + u
      if (comp .ne. 1.0) go to 10
      r1mach = u*2.0
      return
c----------------------- end of function r1mach ------------------------
      end

      subroutine xerrwv (msg, nmes, nerr, level, ni, i1, i2, nr, r1, r2)
      implicit real*8 (a-h,o-z)
      integer nmes, nerr, level, ni, i1, i2, nr,
     1        i, lun, lunit0,lunit1, mesflg, ncpw, nch, nwds
      character msg(nmes)
c      dimension msg(nmes)
      include 'timestuf'
c-----------------------------------------------------------------------
c subroutines xerrwv, xsetf, and xsetun, as given here, constitute
c a simplified version of the slatec error handling package.
c written by a. c. hindmarsh at llnl.  version of march 30, 1987.
c this version is in real.
c
c all arguments are input arguments.
c
c msg    = the message (hollerith literal or integer array).
c nmes   = the length of msg (number of characters).
c nerr   = the error number (not used).
c level  = the error level..
c          0 or 1 means recoverable (control returns to caller).
c          2 means fatal (run is aborted--see note below).
c ni     = number of integers (0, 1, or 2) to be printed with message.
c i1,i2  = integers to be printed, depending on ni.
c nr     = number of reals (0, 1, or 2) to be printed with message.
c r1,r2  = reals to be printed, depending on nr.
c
c note..  this routine is machine-dependent and specialized for use
c in limited context, in the following ways..
c 1. the number of hollerith characters stored per word, denoted
c    by ncpw below, is a data-loaded constant.
c 2. the value of nmes is assumed to be at most 60.
c    (multi-line messages are generated by repeated calls.)
c 3. if level = 2, control passes to the statement   stop
c    to abort the run.  this statement may be machine-dependent.
c 4. r1 and r2 are assumed to be in real*8 and are printed
c    in d21.13 format.
c 5. the common block /eh0001/ below is data-loaded (a machine-
c    dependent feature) with default values.
c    this block is needed for proper retention of parameters used by
c    this routine which the user can reset by calling xsetf or xsetun.
c    the variables in this block are as follows..
c       mesflg = print control flag..
c                1 means print all messages (the default).
c                0 means no printing.
c       lunit  = logical unit number for messages.
c                the default is 6 (machine-dependent).
c-----------------------------------------------------------------------
c the following are instructions for installing this routine
c in different machine environments.
c
c to change the default output unit, change the data statement
c in the block data subprogram below.
c
c for a different number of characters per word, change the
c data statement setting ncpw below, and format 10.  alternatives for
c various computers are shown in comment cards.
c
c for a different run-abort command, change the statement following
c statement 100 at the end.
c-----------------------------------------------------------------------
      common /eh0001/ mesflg, lunit0, lunit1
c-----------------------------------------------------------------------
 
      data ncpw/4/
c-----------------------------------------------------------------------
      if (mesflg .eq. 0) go to 100
 
c get logical unit number. ---------------------------------------------

      lun = lunit1

c get number of words in message. --------------------------------------
c      nch = min0(nmes,60)
c      nwds = nch/ncpw
c      if (nch .ne. nwds*ncpw) nwds = nwds + 1
c write the message. ---------------------------------------------------
c      write (lun, 10) (msg(i),i=1,nwds)
      write(lun,*)'',msg
c-----------------------------------------------------------------------
c the following is valid for ncpw = 8.
  10  format(1x,15a4)
c-----------------------------------------------------------------------
      if (ni .eq. 1) write (lun, 20) i1
 20   format(6x,23hin above message,  i1 =,i10)
      if (ni .eq. 2) write (lun, 30) i1,i2
 30   format(6x,23hin above message,  i1 =,i10,3x,4hi2 =,i10)
      if (nr .eq. 1) write (lun, 40) r1
 40   format(6x,23hin above message,  r1 =,1p,e21.13)
      if (nr .eq. 2) write (lun, 50) r1,r2
 50   format(6x,15hin above,  r1 =,1p,e21.13,3x,4hr2 =,e21.13)
c abort the run if level = 2. ------------------------------------------
 100  if (level .ne. 2) return
      stop
c----------------------- end of subroutine xerrwv ----------------------
      end
      subroutine xsetf (mflag)
c
c this routine resets the print control flag mflag.
c
      integer mflag, mesflg, lunit0, lunit1
      common /eh0001/ mesflg, lunit0, lunit1
c
      if (mflag .eq. 0 .or. mflag .eq. 1) mesflg = mflag
      return
c----------------------- end of subroutine xsetf -----------------------
      end

      subroutine sgefa(a,lda,n,ipvt,info)
      implicit real*8 (a-h,o-z)
      integer lda,n,ipvt(1),info
      real*8 a(lda,1)
c
c     sgefa factors a real  matrix by gaussian elimination.
c
c     on entry
c
c        a       real (lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas raxpy,rscal,idamax
c
c     internal variables
c
      real*8  t
      integer idamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0/a(k,k)
            call rscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call raxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0) info = n
      return
      end

      subroutine sgesl(a,lda,n,ipvt,b,job)
      implicit real*8 (a-h,o-z)
      integer lda,n,ipvt(1),job
      real*8  a(lda,1),b(1)
c
c     dgesl solves the real  system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by sgefa.
c
c     on entry
c
c        a       real(lda, n)
c                the output from sgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       real(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if sgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call sgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call sgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and real*8 functions
c
c     blas raxpy,rdot
c
c     internal variables
c
      real*8  rdot,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call raxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call raxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = rdot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + rdot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end
      subroutine sgbfa(abd,lda,n,ml,mu,ipvt,info)
      implicit real*8 (a-h,o-z)
      integer lda,n,ml,mu,ipvt(1),info
      real*8 abd(lda,1)
c
c     rgbfa factors a real band matrix by elimination.
c
c     rgbfa is usually called by rgbco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c
c     on entry
c
c        abd     real(lda, n)
c                contains the matrix in band storage.  the columns
c                of the matrix are stored in the columns of  abd  and
c                the diagonals of the matrix are stored in rows
c                ml+1 through 2*ml+mu+1 of  abd .
c                see the comments below for details.
c
c        lda     integer
c                the leading dimension of the array  abd .
c                lda must be .ge. 2*ml + mu + 1 .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c                0 .le. ml .lt. n .
c
c        mu      integer
c                number of diagonals above the main diagonal.
c                0 .le. mu .lt. n .
c                more efficient if  ml .le. mu .
c     on return
c
c        abd     an upper triangular matrix in band storage and
c                the multipliers which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgbsl will divide by zero if
c                     called.  use  rcond  in dgbco for a reliable
c                     indication of singularity.
c
c     band storage
c
c           if  a  is a band matrix, the following program segment
c           will set up the input.
c
c                   ml = (band width below the diagonal)
c                   mu = (band width above the diagonal)
c                   m = ml + mu + 1
c                   do 20 j = 1, n
c                      i1 = max0(1, j-mu)
c                      i2 = min0(n, j+ml)
c                      do 10 i = i1, i2
c                         k = i - j + m
c                         abd(k,j) = a(i,j)
c                10    continue
c                20 continue
c
c           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
c           in addition, the first  ml  rows in  abd  are used for
c           elements generated during the triangularization.
c           the total number of rows needed in  abd  is  2*ml+mu+1 .
c           the  ml+mu by ml+mu  upper left triangle and the
c           ml by ml  lower right triangle are not referenced.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas raxpy,rscal,idamax
c     fortran max0,min0
c
c     internal variables
c
      real*8 t
      integer i,idamax,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1
c
c
      m = ml + mu + 1
      info = 0
c
c     zero initial fill-in columns
c
      j0 = mu + 2
      j1 = min0(n,m) - 1
      if (j1 .lt. j0) go to 30
      do 20 jz = j0, j1
         i0 = m + 1 - jz
         do 10 i = i0, ml
            abd(i,jz) = 0.0
   10    continue
   20 continue
   30 continue
      jz = j1
      ju = 0
c
c     gaussian elimination with partial pivoting
c
      nm1 = n - 1
      if (nm1 .lt. 1) go to 130
      do 120 k = 1, nm1
         kp1 = k + 1
c
c        zero next fill-in column
c
         jz = jz + 1
         if (jz .gt. n) go to 50
         if (ml .lt. 1) go to 50
            do 40 i = 1, ml
               abd(i,jz) = 0.0
   40       continue
   50    continue
c
c        find l = pivot index
c
         lm = min0(ml,n-k)
         l = idamax(lm+1,abd(m,k),1) + m - 1
         ipvt(k) = l + k - m
c
c        zero pivot implies this column already triangularized
c
         if (abd(l,k) .eq. 0.0) go to 100
c
c           interchange if necessary
c
            if (l .eq. m) go to 60
               t = abd(l,k)
               abd(l,k) = abd(m,k)
               abd(m,k) = t
   60       continue
c
c           compute multipliers
c
            t = -1.0/abd(m,k)
            call rscal(lm,t,abd(m+1,k),1)
c
c           row elimination with column indexing
c
            ju = min0(max0(ju,mu+ipvt(k)),n)
            mm = m
            if (ju .lt. kp1) go to 90
            do 80 j = kp1, ju
               l = l - 1
               mm = mm - 1
               t = abd(l,j)
               if (l .eq. mm) go to 70
                  abd(l,j) = abd(mm,j)
                  abd(mm,j) = t
   70          continue
               call raxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1)
   80       continue
   90       continue
         go to 110
  100    continue
            info = k
  110    continue
  120 continue
  130 continue
      ipvt(n) = n
      if (abd(m,n) .eq. 0.0d0) info = n
      return
      end

      subroutine sgbsl(abd,lda,n,ml,mu,ipvt,b,job)
      implicit real*8 (a-h,o-z)
      integer lda,n,ml,mu,ipvt(1),job
      real*8 abd(lda,1),b(1)
c
c     dgbsl solves the real band system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by rgbco or rgbfa.
c
c     on entry
c
c        abd     real(lda, n)
c                the output from dgbco or dgbfa.
c
c        lda     integer
c                the leading dimension of the array  abd .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c
c        mu      integer
c                number of diagonals above the main diagonal.
c
c        ipvt    integer(n)
c                the pivot vector from sgbco or sgbfa.
c
c        b       real*8(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b , where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if sgbco has set rcond .gt. 0.0
c        or dgbfa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call sgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call sgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas raxpy,rdot
c     fortran min0
c
c     internal variables
c
      real*8 rdot,t
      integer k,kb,l,la,lb,lm,m,nm1
c
      m = mu + ml + 1
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve l*y = b
c
         if (ml .eq. 0) go to 30
         if (nm1 .lt. 1) go to 30
            do 20 k = 1, nm1
               lm = min0(ml,n-k)
               l = ipvt(k)
               t = b(l)
               if (l .eq. k) go to 10
                  b(l) = b(k)
                  b(k) = t
   10          continue
               call raxpy(lm,t,abd(m+1,k),1,b(k+1),1)
   20       continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/abd(m,k)
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = -b(k)
            call raxpy(lm,t,abd(la,k),1,b(lb),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = rdot(lm,abd(la,k),1,b(lb),1)
            b(k) = (b(k) - t)/abd(m,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (ml .eq. 0) go to 90
         if (nm1 .lt. 1) go to 90
            do 80 kb = 1, nm1
               k = n - kb
               lm = min0(ml,n-k)
               b(k) = b(k) + rdot(lm,abd(m+1,k),1,b(k+1),1)
               l = ipvt(k)
               if (l .eq. k) go to 70
                  t = b(l)
                  b(l) = b(k)
                  b(k) = t
   70          continue
   80       continue
   90    continue
  100 continue
      return
      end
      subroutine raxpy(n,da,dx,incx,dy,incy)
      implicit real*8 (a-h,o-z)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      real*8 dx(1),dy(1),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end
      subroutine rscal(n,da,dx,incx)
      implicit real*8 (a-h,o-z)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      real*8 da,dx(1)
      integer i,incx,m,mp1,n,nincx
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
      real*8 function rdot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      real*8 dx(1),dy(1),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      rdot = 0.0
      dtemp = 0.0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      rdot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 rdot = dtemp
      return
      end

      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c
      real*8 dx(1),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n .lt. 1 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = abs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(abs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = abs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = abs(dx(1))
      do 30 i = 2,n
         if(abs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = abs(dx(i))
   30 continue
      return
      end
      block data lsd
c-----------------------------------------------------------------------
c this data subprogram loads variables into the internal common
c blocks used by the odepack solvers.  the variables are
c defined as follows..
c   illin  = counter for the number of consecutive times the package
c            was called with illegal input.  the run is stopped when
c            illin reaches 5.
c   ntrep  = counter for the number of consecutive times the package
c            was called with istate = 1 and tout = t.  the run is
c            stopped when ntrep reaches 5.
c   mesflg = flag to control printing of error messages.  1 means print,
c            0 means no printing.
c   lunit  = default value of logical unit number for printing of error
c            messages.
c-----------------------------------------------------------------------
      integer illin, iduma, ntrep, idumb, iowns, icomm, mesflg, lunit0,
     1        lunit1
      real*8 rowns, rcomm
      common /ls0001R/ rowns(209), rcomm(9)
      common /ls0001I/
     1   illin, iduma(10), ntrep, idumb(2), iowns(6), icomm(19)
      common /eh0001/ mesflg, lunit0, lunit1
c      data illin/0/, ntrep/0/
      data mesflg/1/, lunit0/16/, lunit1/16/
c
c----------------------- end of block data -----------------------------
      end
