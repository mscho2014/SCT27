c
c
c
c
c          ******   *******   *****
c           *    *  *        *     *
c           *    *  *        *
c           *    *  ****     *
c           *    *  *        *
c           *    *  *        *     *
c          ******   *******   *****
c
      subroutine dec (n, ndim, a, ip, ier)
c
      implicit real*8 (a-h,o-z)
      integer n, ndim, ip(n), ier
      dimension a(ndim,n)
 
c  matrix triangularization by gauss elimination with partial pivoting.
 
c  INPUT..
 
c     n = order of matrix.
c     ndim = declared first dimension of array  a.
c     a = matrix to be triangularized.
 
c  OUTPUT..
 
c     a(i,j), i.le.j = upper triangular factor, u .
c     a(i,j), i.gt.j = multipliers = lower triangular factor, i - l.
c     ip(k), k.lt.n = index of k-th pivot row.
c     ier = 0 if matrix a is nonsingular, or k if found to be
c           singular at stage k.
 
c      row interchanges are finished in u, only partly in l.
c      use  sol  to obtain solution of linear system.
c      if ier .ne. 0, a is singular, sol will divide by zero.
c
c  reference:  a. c. hindmarsh, l. j. sloan, and p. f. dubois,
c              dec/sol:  solution of dense systems of linear
c              algebraic equations,
c              lawrence livermore laboratory report ucid-30137, rev. 1,
c              december 1978.
c
      ier = 0
      if (n .eq. 1) go to 70
      nm1 = n - 1
      do 60 k = 1,nm1
        
        kp1 = k + 1
 
c  find the pivot in column k.  search rows k to n.
 
        m = k
        do 10 i = kp1,n
 10       if (abs(a(i,k)) .gt. abs(a(m,k))) m = i
        m = k
        ip(k) = m
 
c  interchange elements in rows k and m.
 
        t = a(m,k)
        if (m .eq. k) go to 20
        a(m,k) = a(k,k)
        a(k,k) = t
 20     if (t .eq. 0.00)  then
           write(16,'(" Warning: -dec- ",1p,10e11.3)')t, m, k
           go to 80
        endif
 
c  store multipliers in a(i,k), i = k+1,...,n.
 
        t = 1.d0/t

        do 30 i = kp1,n
 30       a(i,k) = -a(i,k)*t
 
c  apply multipliers to other columns of a.
 
        do 50 j = kp1,n
          t = a(m,j)
          a(m,j) = a(k,j)
          a(k,j) = t
          if (t .eq. 0.00) go to 50
          do 40 i = kp1,n
 40         a(i,j) = a(i,j) + a(i,k)*t
 50       continue
c
c
 60     continue
c
 70   k = n
      if (a(n,n) .eq. 0.00) go to 80
      return
 80   ier = k

      return
      end
c
c           ******   *****   *
c          *        *     *  *
c          *        *     *  *
c           *****   *     *  *
c                *  *     *  *
c                *  *     *  *
c          ******    *****   *******
c
      subroutine sol (n, ndim, a, b, ip)
 
c  solution of linear system a*x = b using output of dec.
c  input..
c     n = order of matrix.
c     ndim = declared first dimension of array  a.
c     a = triangularized matrix obtained from dec.
c     b = right hand side vector.
c     ip = pivot information vector obtained from dec.
c  do not use if dec has set ier .ne. 0.
c  output..
c     b = solution vector, x .
 
      implicit real*8 (a-h,o-z)
      dimension ip(n)
      dimension a(ndim,n), b(n)
 
      if (n .eq. 1) go to 50
      nm1 = n - 1
 
c  apply row permutations and multipliers to b. ------------------------
 
      do 20 k = 1,nm1
        kp1 = k + 1
        m = ip(k)
        t = b(m)
        b(m) = b(k)
        b(k) = t
        do 10 i = kp1,n
 10       b(i) = b(i) + a(i,k)*t
 20     continue
 
c  back solve. ---------------------------------------------------------
 
      do 40 kb = 1,nm1
        km1 = n - kb
        k = km1 + 1
        b(k) = b(k)/a(k,k)
        t = -b(k)
        do 30 i = 1,km1
 30       b(i) = b(i) + a(i,k)*t
 40     continue
 50   b(1) = b(1)/a(1,1)
      return
      end
c
c
c
c           ******   *****   *        *     *  *******
c          *        *     *  *        *     *  *
c          *        *     *  *         *   *   *
c           *****   *     *  *         *   *   ****
c                *  *     *  *          * *    *
c                *  *     *  *          * *    *
c          ******    *****   *******     *     *******
c
c
c
c
      subroutine solve(n)
c
c     dummy routine to drive dec/sol and not overwrite matrices
c
c____popstuf
 
      implicit real*8 (a-h,o-z)
      include 'timestuf'
      include 'popstuf'

      dimension ip(numpp)
      dimension ao(numpp,numpp),bo(numpp)
c

      do 1 i=1,n
        bo(i) = b(i)
        do 2 j=1,n
          ao(i,j) = a(i,j)
    2   continue
    1 continue


c      write(*,900)(ao(i,i),i=1,n)
c 900  format(1p,10e10.2)
      call dec (n,numpp,ao,ip,ier)

      if (ier .ne. 0)  then
         write(6,'(" Exiting:Trouble -solve-",i5)') ier
         write(16,'(" Exiting:Trouble -solve-",i5)') ier
         stop
      endif
      call sol (n,numpp,ao,bo,ip)
c
      do 4 i=1,n
        x(i) = bo(i)
   4  continue

c
      return
      end
 
