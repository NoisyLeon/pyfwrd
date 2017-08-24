      SUBROUTINE sphbes(n,x,sj,sy,sjp,syp)
      INTEGER n
      REAL*8 sj,sjp,sy,syp,x
CU    USES bessjy
      REAL*8 factor,order,rj,rjp,ry,ryp,RTPIO2
      PARAMETER (RTPIO2=1.253314137315500d0)
      if(n.lt.0.d0.or.x.le.0.d0)pause 'bad arguments in sphbes'
      order=n+0.5d0
      call bessjy(x,order,rj,ry,rjp,ryp)
      factor=RTPIO2/sqrt(x)
      sj=factor*rj
      sy=factor*ry
      sjp=factor*rjp-sj/(2.0d0*x)
      syp=factor*ryp-sy/(2.0d0*x)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software W"..
