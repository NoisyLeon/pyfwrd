      FUNCTION chebev(a,b,c,m,x)
      INTEGER m
      REAL*8 chebev,a,b,x,c(m)
      INTEGER j
      REAL*8 d,dd,sv,y,y2
      if ((x-a)*(x-b).gt.0.d0) pause 'x not in range in chebev'
      d=0.d0
      dd=0.d0
      y=(2.d0*x-a-b)/(b-a)
      y2=2.d0*y
      do 11 j=m,2,-1
        sv=d
        d=y2*d-dd+c(j)
        dd=sv
11    continue
      chebev=y*d-dd+0.5d0*c(1)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software W"..
