c-------------------------------------------------------------------------------
      subroutine bfs(l,xsq,eps,fp)
c-------------------------------------------------------------------------------
c  this routine calculates spherical bessel function of the ist kind.
c  fp is equivalent to (r*dj/dr)/j
c  where r is radius and j is the sbf of order l and argument x=k*r
c  the technique employs the continued fraction approach
c  described in w. lentz's article in applied qptics, vol.15, #3, 1976
c-------------------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      save
      real*8 numer,nu
      if(xsq.le.0.d0) goto 10
      x=dsqrt(xsq)
      lp1=l+1
      rx=2.0d0/x
      nu=lp1-0.5d0
      rj=nu*rx
      rx=-rx
      denom=(nu+1.d0)*rx
      numer=denom+1.0d0/rj
      rj=rj*numer/denom
      nm1=1
    2 nm1=nm1+1
      rx=-rx
      a3=(nu+nm1)*rx
      denom=a3+1.d0/denom
      numer=a3+1.d0/numer
      ratio=numer/denom
      rj=rj*ratio
      if(dabs(dabs(ratio)-1.d0).gt.eps) goto 2
      fp=rj*x-lp1
      return
c  series solution
   10 f=1.d0
      fp=l
      a=1.d0
      b=l+l+1.d0
      c=2.d0
      d=l+2.d0
   15 a=-a*xsq/(c*(b+c))
      f=f+a
      fp=fp+a*d
      if(dabs(a*d).lt.eps) goto 20
      c=c+2.d0
      d=d+2.d0
      goto 15
   20 fp=fp/f
      return
      end
c-------------------------------------------------------------------------------
