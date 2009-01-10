      integer function IGRAPHARPACKIDAMAX(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified to correct problem with negative increment, 8/21/90.
c
      double precision dx(1),dmax
      integer i,incx,ix,n
c
      IGRAPHARPACKIDAMAX = 0
      if( n .lt. 1 ) return
      IGRAPHARPACKIDAMAX = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      dmax = dabs(dx(ix))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         IGRAPHARPACKIDAMAX = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         IGRAPHARPACKIDAMAX = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
