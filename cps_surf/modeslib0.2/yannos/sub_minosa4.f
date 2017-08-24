      subroutine rprop(jf,jl,f)     
      use shanks
      use bits
      use eifx, a=>ar
      use param_modele
      use minos
      use yannos_flag
      implicit real*8(a-h,o-z)
c*** propagates soln ,f, for radial modes from jf to jl ***
      save
      real*8 nn
      dimension h(2,10),s(2),f(2)
      maxo1=maxo-1
      y=r(jf)
      vy=dsqrt((flam(jf)+2.d0*fmu(jf))/rho(jf))
      i=jf
      go to 50
   10 iq=i
      i=i+1
      x=y
      y=r(i)
      if(y.eq.x) goto 50
      qff=1.d0+xlam(iq)*fct
      qll=1.d0+qshear(iq)*fct
      qaa=1.d0+xa2(iq)*fct
      vx=vy
      vy=dsqrt((flam(i)+2.d0*fmu(i))/rho(i))
      q=dmax1(w/vx+1.d0/x,w/vy+1.d0/y)
      del=step(maxo)/q
      dxs=0.d0
   15 y=x+del
      if(y.gt.r(i)) y=r(i)
      dx=y-x
      if(dx.ne.dxs) call baylis(q,maxo1)
      dxs=dx
      s(1)=f(1)
      s(2)=f(2)
      do 40 ni=1,in
      z=x+c(ni)
      t=z-r(iq)
      ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
      gr=g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq)))
      ff=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq))))*qff
      if(ifanis.ne.0) goto 30
      nn=(lcon(iq)+t*(lspl(1,iq)+t*(lspl(2,iq)+t*lspl(3,iq))))*qll
      cc=ff+nn+nn
      aa=cc
      goto 35
   30 nn=(ncon(iq)+t*(nspl(1,iq)+t*(nspl(2,iq)+t*nspl(3,iq))))*qll
      cc=(ccon(iq)+t*(cspl(1,iq)+t*(cspl(2,iq)+t*cspl(3,iq))))*qaa
      aa=(acon(iq)+t*(aspl(1,iq)+t*(aspl(2,iq)+t*aspl(3,iq))))*qaa
   35 z=1.d0/z
      if ((.not.cancel_gravity).and.kg.eq.0) then
!approximation de cowling! 
         a21=-ro*wsq+4.d0*z*(z*(aa-nn-ff*ff/cc)-ro*(gr-ro/z))
      else
         a21=-ro*wsq+4.d0*z*(z*(aa-nn-ff*ff/cc)-ro*gr)
      endif
      h(1,ni)=(f(2)-2.d0*ff*z*f(1))/cc
      h(2,ni)=a21*f(1)+2.d0*z*f(2)*(ff/cc-1.d0)
   40 call rkdot(f,s,h,2,ni)
      if(knsw.ne.1) goto 45
!++++++++++++++++++++++++++++++++++yann
      select case (cond_limite)
      case(1)
!bord libre
         fp=a21*f(1)+2.d0*z*f(2)*(ff/cc-1.d0)
         call trknt(s(2),h(2,1),f(2),fp,x,y)
      case(2)
!bord rigide         
         fp=(f(2)-2.d0*ff*z*f(1))/cc
         call trknt(s(1),h(1,1),f(1),fp,x,y)
      end select
!+++++++++++++++++++++++++++++++++++fin yann
   45 x=y
      if(y.ne.r(i)) go to 15
   50 a(1,i)=f(1)
      a(2,i)=f(2)
      if(i.ne.jl) go to 10
      return
      end

      subroutine gauslv(r1,r2,iq,fint,nint)
c*** fifth order gauss-legendre integration ***
      implicit real*8(a-h,o-z)
      save
      dimension fint(*),vals(4),vals1(4),sum(4),w(2),x(2)
      data w,x/.478628670499366d0,.236926885056189d0,
     +         .538469310105683d0,.906179845938664d0/
      y1=.5d0*(r2+r1)
      y2=.5d0*(r2-r1)
      call intgds(y1,iq,vals)
      do 5 j=1,nint
    5 sum(j)=.568888888888889d0*vals(j)
      do 10 i=1,2
      t1=x(i)*y2
      call intgds(y1+t1,iq,vals)
      call intgds(y1-t1,iq,vals1)
      do 10 j=1,nint
   10 sum(j)=sum(j)+w(i)*(vals(j)+vals1(j))
      do 15 j=1,nint
   15 fint(j)=fint(j)+y2*sum(j)
      return
      end
!

      subroutine baylis(q,maxo1)
c    baylis returns the coefficients for rks integration.
c    see e. baylis shanks(1966 a. m. s.) and references therein for the
c    coefficients. the eight runge-kutta-shanks formulae are (1-1) (2-2)
c    (3-3) (4-4) (5-5) (6-6) (7-7) (8-10). for orders greater than 4 the
c    formulae are approximate rather than exact so incurring less roundoff.
!      use shanks, i=>in
      use shanks
      implicit real*8(a-h,o-z)
      save
      ds=q*dabs(dx)
      do 10 j=1,maxo1
      if(ds.gt.step(j)) go to 10
      in=j
      go to 15
   10 continue
      in=maxo
   15 c(1)=0.d0
      go to (1,2,3,4,5,6,7,8),in
    1 b(1)=dx
      return
    2 c(2)=dx
      b(1)=dx
      b(2)=.5d0*dx
      b(3)=1.d0
      return
    3 c(2)=.5d0*dx
      c(3)=dx
      b(1)=c(2)
      b(2)=-dx
      b(3)=-2.d0
      b(4)=.16666666666667d0*dx
      b(5)=4.d0
      b(6)=1.d0
      return
    4 c(2)=.01d0*dx
      c(3)=.6d0*dx
      c(4)=dx
      b(1)=c(2)
      b( 2)=-.17461224489790d+02*dx
      b( 3)=-.10343618513324d+01
      b( 4)= .59691275167780d+02*dx
      b( 5)=-.10140620414448d+01
      b( 6)= .30814908546230d-01
      b( 7)=-.25555555555556d+01*dx
      b( 8)=-.11165449632656d+01
      b( 9)=-.22568165070006d+00
      b(10)=-.49077733860351d-01
      return
    5 c( 2)= 1.1111111111111d-04*dx
      c( 3)= 3.0d-01*dx
      c( 4)= 7.5d-01*dx
      c( 5)= dx
      b( 1)=c(2)
      b( 2)=-.40470000000000d+03*dx
      b( 3)=-.10007412898443d+01
      b( 4)= .25301250000000d+04*dx
      b( 5)=-.10004446420631d+01
      b( 6)= .74107010523195d-03
      b( 7)=-.11494333333333d+05*dx
      b( 8)=-.10004929965491d+01
      b( 9)= .52629261224803d-03
      b(10)=-.12029545422812d-03
      b(11)= .92592592592593d-01*dx
      b(12)= .00000000000000d+00
      b(13)= .47619047619048d+01
      b(14)= .42666666666667d+01
      b(15)= .77142857142857d+00
      return
    6 c(2)=3.3333333333333d-03*dx
      c(3)=.2d0*dx
      c(4)=.6d0*dx
      c(5)=9.3333333333333d-01*dx
      c(6)=dx
      b( 1)=c(2)
      b( 2)=-.58000000000000d+01*dx
      b( 3)=-.10344827586207d+01
      b( 4)= .64600000000000d+02*dx
      b( 5)=-.10216718266254d+01
      b( 6)= .30959752321982d-01
      b( 7)=-.62975802469136d+03*dx
      b( 8)=-.10226149961576d+01
      b( 9)= .24906685695466d-01
      b(10)=-.37737402568887d-02
      b(11)=-.54275714285714d+04*dx
      b(12)=-.10225567867765d+01
      b(13)= .25375487829097d-01
      b(14)=-.31321559234596d-02
      b(15)= .12921040478749d-03
      b(16)= .53571428571429d-01*dx
      b(17)= .00000000000000d+00
      b(18)= .61868686868687d+01
      b(19)= .77777777777778d+01
      b(20)= .40909090909091d+01
      b(21)=-.38888888888889d+00
      return
    7 c(2)=5.2083333333333d-03*dx
      c(3)=1.6666666666667d-01*dx
      c(4)=.5d0*dx
      c(5)=dx
      c(6)=8.3333333333333d-01*dx
      c(7)=dx
      b( 1)=c(2)
      b( 2)=-.25000000000000d+01*dx
      b( 3)=-.10666666666667d+01
      b( 4)= .26166666666667d+02*dx
      b( 5)=-.10421204027121d+01
      b( 6)= .61228682966918d-01
      b( 7)=-.64500000000000d+03*dx
      b( 8)=-.10450612653163d+01
      b( 9)= .51262815703925d-01
      b(10)=-.77519379844961d-02
      b(11)=-.93549382716049d+02*dx
      b(12)=-.10450293206756d+01
      b(13)= .48394546673620d-01
      b(14)=-.11877268228307d-01
      b(15)=-.39590894094358d-03
      b(16)= .35111904761905d+03*dx
      b(17)=-.10446476812124d+01
      b(18)= .52479782656724d-01
      b(19)=-.71200922221468d-02
      b(20)=-.61029361904114d-03
      b(21)= .27463212856852d-02
      b(22)= .46666666666667d-01*dx
      b(23)= .57857142857143d+01
      b(24)= .78571428571429d+01
      b(25)= .00000000000000d+00
      b(26)= b(23)
      b(27)= .10000000000000d+01
      return
    8 c(2)=.14814814814815d0*dx
      c(3)=.22222222222222d0*dx
      c(4)=.33333333333333d0*dx
      c(5)= .5d0*dx
      c(6)=.66666666666667d0*dx
      c(7)=.16666666666667d0*dx
      c(8)=dx
      c(9)=.83333333333333d0*dx
      c(10)=dx
      b( 1)=c(2)
      b( 2)= .55555555555556d-01*dx
      b( 3)= .30000000000000d+01
      b( 4)= .83333333333333d-01*dx
      b( 5)= .00000000000000d+00
      b( 6)= .30000000000000d+01
      b( 7)= .12500000000000d+00*dx
      b( 8)= .00000000000000d+00
      b( 9)= .00000000000000d+00
      b(10)= .30000000000000d+01
      b(11)= .24074074074074d+00*dx
      b(12)= .00000000000000d+00
      b(13)=-.20769230769231d+01
      b(14)= .32307692307692d+01
      b(15)= .61538461538461d+00
      b(16)= .90046296296295d-01*dx
      b(17)= .00000000000000d+00
      b(18)=-.13881748071980d+00
      b(19)= .24832904884319d+01
      b(20)=-.21182519280206d+01
      b(21)= .62467866323908d+00
      b(22)=-.11550000000000d+02*dx
      b(23)=-.35064935064935d+00
      b(24)= .50389610389610d+01
      b(25)=-.28398268398268d+01
      b(26)= .52813852813853d+00
      b(27)=-.34632034632035d+01
      b(28)=-.44097222222222d+00*dx
      b(29)=-.14173228346457d+00
      b(30)= .53385826771654d+01
      b(31)=-.35905511811023d+01
      b(32)= .70866141732284d-01
      b(33)=-.45354330708661d+01
      b(34)=-.31496062992126d-01
      b(35)= .18060975609756d+01*dx
      b(36)=-.54692775151925d-01
      b(37)= .47967589466576d+01
      b(38)=-.22795408507765d+01
      b(39)= .48615800135044d-01
      b(40)=-.34031060094530d+01
      b(41)=-.40513166779204d-01
      b(42)= .48615800135044d+00
      b(43)= .48809523809524d-01*dx
      b(44)= .65853658536585d+00
      b(45)= .66341463414634d+01
      b(46)= .52682926829268d+01
      in=10
      return
      end
      subroutine steps(eps)
c*** computes 8 dimensionless step sizes for rks integration
      use shanks
      implicit real*8(a-h,o-z)
      save
      ps=dlog(eps)
      fac=1.d0
      do 2 n=1,8
      fn=n+1
      fac=fac*fn
      x=(dlog(fac)+ps)/fn
      x=dexp (x)
      s=x
      do 1 i=1,n
    1 s=x*dexp(-s/fn)
    2 step(n)=s
      return
      end
      subroutine drspln(i1,i2,x,y,q,f)
      implicit real*8(a-h,o-z)
      save
c   rspln computes cubic spline interpolation coefficients
c   for y(x) between grid points i1 and i2 saving them in q.  the
c   interpolation is continuous with continuous first and second
c   derivitives.  it agrees exactly with y at grid points and with the
c   three point first derivitives at both end points (i1 and i2).
c   x must be monotonic but if two successive values of x are equal
c   a discontinuity is assumed and seperate interpolation is done on
c   each strictly monotonic segment.  the arrays must be dimensioned at
c   least - x(i2), y(i2), q(3,i2), and f(3,i2).  f is working storage
c   for rspln.
c                                                     -rpb
      dimension x(*),y(*),q(3,1),f(3,1),yy(3)
      equivalence (yy(1),y0)
      data yy/3*0.d0/
      j1=i1+1
      y0=0.d0
c   bail out if there are less than two points total.
      if(i2-i1)13,17,8
 8    a0=x(j1-1)
c   search for discontinuities.
      do 3 i=j1,i2
      b0=a0
      a0=x(i)
      if(a0-b0)3,4,3
 3    continue
 17   j1=j1-1
      j2=i2-2
      go to 5
 4    j1=j1-1
      j2=i-3
c   see if there are enough points to interpolate (at least three).
 5    if(j2+1-j1)9,10,11
c   only two points.  use linear interpolation.
 10   j2=j2+2
      y0=(y(j2)-y(j1))/(x(j2)-x(j1))
      do 15 j=1,3
      q(j,j1)=yy(j)
 15   q(j,j2)=yy(j)
      go to 12
c   more than two points.  do spline interpolation.
 11   a0=0.d0
      h=x(j1+1)-x(j1)
      h2=x(j1+2)-x(j1)
      y0=h*h2*(h2-h)
      h=h*h
      h2=h2*h2
c   calculate derivitive at near end.
      b0=(y(j1)*(h-h2)+y(j1+1)*h2-y(j1+2)*h)/y0
      b1=b0
c   explicitly reduce banded matrix to an upper banded matrix.
      do 1 i=j1,j2
      h=x(i+1)-x(i)
      y0=y(i+1)-y(i)
      h2=h*h
      ha=h-a0
      h2a=h-2.d0*a0
      h3a=2.d0*h-3.*a0
      h2b=h2*b0
      q(1,i)=h2/ha
      q(2,i)=-ha/(h2a*h2)
      q(3,i)=-h*h2a/h3a
      f(1,i)=(y0-h*b0)/(h*ha)
      f(2,i)=(h2b-y0*(2.d0*h-a0))/(h*h2*h2a)
      f(3,i)=-(h2b-3.d0*y0*ha)/(h*h3a)
      a0=q(3,i)
 1    b0=f(3,i)
c   take care of last two rows.
      i=j2+1
      h=x(i+1)-x(i)
      y0=y(i+1)-y(i)
      h2=h*h
      ha=h-a0
      h2a=h*ha
      h2b=h2*b0-y0*(2.d0*h-a0)
      q(1,i)=h2/ha
      f(1,i)=(y0-h*b0)/h2a
      ha=x(j2)-x(i+1)
      y0=-h*ha*(ha+h)
      ha=ha*ha
c   calculate derivitive at far end.
      y0=(y(i+1)*(h2-ha)+y(i)*ha-y(j2)*h2)/y0
      q(3,i)=(y0*h2a+h2b)/(h*h2*(h-2.d0*a0))
      q(2,i)=f(1,i)-q(1,i)*q(3,i)
c   solve upper banded matrix by reverse iteration.
      do 2 j=j1,j2
      k=i-1
      q(1,i)=f(3,k)-q(3,k)*q(2,i)
      q(3,k)=f(2,k)-q(2,k)*q(1,i)
      q(2,k)=f(1,k)-q(1,k)*q(3,k)
 2    i=k
      q(1,i)=b1
c   fill in the last point with a linear extrapolation.
 9    j2=j2+2
      do 14 j=1,3
 14   q(j,j2)=yy(j)
c   see if this discontinuity is the last.
 12   if(j2-i2)6,13,13
c   no.  go back for more.
 6    j1=j2+2
      if(j1-i2)8,8,7
c   there is only one point left after the latest discontinuity.
 7    do 16 j=1,3
 16   q(j,i2)=yy(j)
c   fini.
 13   return
      end
      subroutine dsplin(n,x,y,q,f)
      implicit real*8(a-h,o-z)
      save
      dimension x(*),y(*),q(3,1),f(3,1),yy(3)
      equivalence (yy(1),y0)
      data yy/3*0.d0/
      a0=0.d0
      j2=n-2
      h=x(2)-x(1)
      h2=x(3)-x(1)
      y0=h*h2*(h2-h)
      h=h*h
      h2=h2*h2
      b0=(y(1)*(h-h2)+y(2)*h2-y(3)*h)/y0
      b1=b0
      do 5 i=1,j2
      h=x(i+1)-x(i)
      y0=y(i+1)-y(i)
      h2=h*h
      ha=h-a0
      h2a=h-2.d0*a0
      h3a=2.d0*h-3.*a0
      h2b=h2*b0
      q(1,i)=h2/ha
      q(2,i)=-ha/(h2a*h2)
      q(3,i)=-h*h2a/h3a
      f(1,i)=(y0-h*b0)/(h*ha)
      f(2,i)=(h2b-y0*(2.d0*h-a0))/(h*h2*h2a)
      f(3,i)=-(h2b-3.d0*y0*ha)/(h*h3a)
      a0=q(3,i)
    5 b0=f(3,i)
      i=j2+1
      h=x(i+1)-x(i)
      y0=y(i+1)-y(i)
      h2=h*h
      ha=h-a0
      h2a=h*ha
      h2b=h2*b0-y0*(2.d0*h-a0)
      q(1,i)=h2/ha
      f(1,i)=(y0-h*b0)/h2a
      ha=x(j2)-x(i+1)
      y0=-h*ha*(ha+h)
      ha=ha*ha
      y0=(y(i+1)*(h2-ha)+y(i)*ha-y(j2)*h2)/y0
      q(3,i)=(y0*h2a+h2b)/(h*h2*(h-2.d0*a0))
      q(2,i)=f(1,i)-q(1,i)*q(3,i)
      do 10 j=1,j2
      k=i-1
      q(1,i)=f(3,k)-q(3,k)*q(2,i)
      q(3,k)=f(2,k)-q(2,k)*q(1,i)
      q(2,k)=f(1,k)-q(1,k)*q(3,k)
   10 i=k
      q(1,i)=b1
      do 15 j=1,3
   15 q(j,n)=yy(j)
      return
      end
!

!
      subroutine rkdot(f,s,h,nvec,ni)
      use shanks
c*** performs dot product with rks coefficients ***
      implicit real*8(a-h,o-z)
      save
      dimension s(*),f(*),h(nvec,*)
      goto (1,2,3,4,5,6,7,8,9,10),ni
    1 do 21 j=1,nvec
   21 f(j)=s(j)+b(1)*h(j,1)
      return
    2 do 22 j=1,nvec
   22 f(j)=s(j)+b(2)*(h(j,1)+b(3)*h(j,2))
      return
    3 do 23 j=1,nvec
   23 f(j)=s(j)+b(4)*(h(j,1)+b(5)*h(j,2)+b(6)*h(j,3))
      return
    4 do 24 j=1,nvec
   24 f(j)=s(j)+b(7)*(h(j,1)+b(8)*h(j,2)+b(9)*h(j,3)+b(10)*h(j,4))
      return
    5 do 25 j=1,nvec
   25 f(j)=s(j)+b(11)*(h(j,1)+b(12)*h(j,2)+b(13)*h(j,3)+b(14)*h(j,4)+
     +b(15)*h(j,5))
      return
    6 do 26 j=1,nvec
   26 f(j)=s(j)+b(16)*(h(j,1)+b(17)*h(j,2)+b(18)*h(j,3)+b(19)*h(j,4)+
     +b(20)*h(j,5)+b(21)*h(j,6))
      return
    7 do 27 j=1,nvec
   27 f(j)=s(j)+b(22)*(h(j,1)+b(23)*h(j,3)+b(24)*h(j,4)+b(25)*h(j,5)+
     +b(26)*h(j,6)+b(27)*h(j,7))
      return
    8 do 28 j=1,nvec
   28 f(j)=s(j)+b(28)*(h(j,1)+b(29)*h(j,3)+b(30)*h(j,4)+b(31)*h(j,5)+
     +b(32)*h(j,6)+b(33)*h(j,7)+b(34)*h(j,8))
      return
    9 do 29 j=1,nvec
   29 f(j)=s(j)+b(35)*(h(j,1)+b(36)*h(j,3)+b(37)*h(j,4)+b(38)*h(j,5)+
     +b(39)*h(j,6)+b(40)*h(j,7)+b(41)*h(j,8)+b(42)*h(j,9))
      return
   10 do 30 j=1,nvec
   30 f(j)=s(j)+b(43)*(h(j,1)+h(j,10)+b(44)*(h(j,4)+h(j,6))+
     +b(45)*h(j,5)+b(46)*(h(j,7)+h(j,9)))
      return
      end
      subroutine intgds(rr,iq,vals)      
c*** interpolates integrands for normalisation,cg,q etc..for use with gauslv.
      use bits
      use eifx
      use param_modele
      use yannos_flag
      implicit real*8(a-h,o-z)
      save
      real*8  nn,ll
      dimension q(3),qp(3),vals(*)
      data d1,d2,d3,d4,d5,d6,d7/.111111111111111d0,
     + 0.066666666666667d0,0.666666666666667d0,1.333333333333333d0,
     + 2.666666666666667d0,3.333333333333333d0,5.333333333333333d0/
      t=rr-r(iq)
      hn=1.d0/(r(iq+1)-r(iq))
      hsq=hn*hn
      qff=1.d0+xlam(iq)*fct
      qll=1.d0+qshear(iq)*fct
      iq1=iq+1
      ifun=3
      if(jcom.ne.3) ifun=1
      do 10 i=1,ifun
      i2=2*i
      i1=i2-1
      a=((ar(i2,iq)+ar(i2,iq1))+2.d0*hn*(ar(i1,iq)-ar(i1,iq1)))*hsq
      b=-(2.d0*ar(i2,iq)+ar(i2,iq1))*hn-3.d0*(ar(i1,iq)-ar(i1,iq1))*hsq
      q(i)=(ar(i1,iq)+t*(ar(i2,iq)+t*(b+t*a)))/rr
   10 qp(i)=ar(i2,iq)+t*(2.d0*b+t*3.d0*a)
      rro=(rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq))))*rr
      gr=g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq)))
      ff=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq))))*qff
      ll=(lcon(iq)+t*(lspl(1,iq)+t*(lspl(2,iq)+t*lspl(3,iq))))*qll
      if(ifanis.ne.0) goto 15
      nn=ll
      cc=ff+ll+ll
      aa=cc
      goto 20
   15 qaa=1.d0+xa2(iq)*fct
      nn=(ncon(iq)+t*(nspl(1,iq)+t*(nspl(2,iq)+t*nspl(3,iq))))*qll
      cc=(ccon(iq)+t*(cspl(1,iq)+t*(cspl(2,iq)+t*cspl(3,iq))))*qaa
      aa=(acon(iq)+t*(aspl(1,iq)+t*(aspl(2,iq)+t*aspl(3,iq))))*qaa
   20 qrka=d1*(4.d0*(aa+ff-nn)+cc)
     1     *(qkappa(iq)+t*hn*(qkappa(iq1)-qkappa(iq)))
      qrmu=d2*(aa+cc-2.d0*ff+5.d0*nn+6.d0*ll)
     1     *(qshear(iq)+t*hn*(qshear(iq1)-qshear(iq)))
      if(jcom.ne.3) goto 25
      q1sq=q(1)*q(1)
      q2sq=q(2)*q(2)
      vals(1)=rr*rro*(q1sq+q2sq)
      fac=(fl+.5d0)/sfl3
      vals(2)=(sfl3*(ll*q1sq+aa*q2sq)+q(2)*((rro*gr+2.d0*(nn-aa-ll)+ff)
     +   *q(1)+rro*q(3)-ff*qp(1))+ll*qp(2)*q(1))*fac
     +   +.25d0*q(3)*(qp(3)+fl*q(3))
      t2=qrka+d7*qrmu
      t3=qrka+d4*qrmu
      t4=qrka+d6*qrmu
      t5=qrka-d5*qrmu
      t6=qrka-d3*qrmu
      vals(3)=.5d0*((fl3*qrmu+t2)*q1sq+(2.d0*qrmu+fl3*t3)*q2sq)
     1 -q(1)*sfl3*t4*q(2)+q(1)*(t5*qp(1)+sfl3*qrmu*qp(2))+q(2)*(-2.d0*
     2 qrmu*qp(2)-sfl3*t6*qp(1))+.5d0*(t3*qp(1)*qp(1)+qrmu*qp(2)*qp(2))
      if (cancel_gravity) then
      vals(4)=.5d0*((fl3*ll+4.d0*(aa-nn-ff)+cc)*q1sq+
     +(4.d0*ll-nn-nn+fl3*aa)*q2sq +fl*fl*.25d0*q(3)*q(3)+cc*qp(1)*qp(1)+
     +ll*qp(2)*qp(2)+.25d0*qp(3)*qp(3))+q(3)*(rro*sfl3*q(2)+fl*.25d0*qp
     +(3))+q(1)*(sfl3*(rro*gr+2.d0*(nn-aa-ll)+ff)*q(2)+rro*(qp(3)-q(3))+
     +(ff+ff-cc)*qp(1)+sfl3*ll*qp(2))-q(2)*(sfl3*ff*qp(1)+(ll+ll)*qp(2))
      else
      vals(4)=.5d0*((fl3*ll+4.d0*(rro*(rro-gr)+aa-nn-ff)+cc)*q1sq+
     +(4.d0*ll-nn-nn+fl3*aa)*q2sq +fl*fl*.25d0*q(3)*q(3)+cc*qp(1)*qp(1)+
     +ll*qp(2)*qp(2)+.25d0*qp(3)*qp(3))+q(3)*(rro*sfl3*q(2)+fl*.25d0*qp
     +(3))+q(1)*(sfl3*(rro*gr+2.d0*(nn-aa-ll)+ff)*q(2)+rro*(qp(3)-q(3))+
     +(ff+ff-cc)*qp(1)+sfl3*ll*qp(2))-q(2)*(sfl3*ff*qp(1)+(ll+ll)*qp(2))
      endif
      return
   25 q(1)=q(1)*rr
      vals(1)=rr*rro*q(1)*q(1)
      if(jcom.eq.1) goto 30
      vals(2)=nn*q(1)*q(1)
      t1=(rr*qp(1)-q(1))**2
      t2=(fl3-2.d0)*q(1)*q(1)
      vals(3)=(t1+t2)*qrmu
      vals(4)=t1*ll+t2*nn
      return
   30 t1=(rr*qp(1)+2.d0*q(1))**2
      t2=d4*(rr*qp(1)-q(1))**2
      vals(2)=t1*qrka+t2*qrmu
      if ((.not.cancel_gravity).and.kg.eq.0) then
!approximation de cowling! on rajoute le 4PIGrho^2U^2 manquant
         vals(3)=rr*qp(1)*(cc*rr*qp(1)+4.d0*ff*q(1))+4.d0*q(1)*q(1)
     +       *(aa-nn-rro*(gr-rro))
      else
         vals(3)=rr*qp(1)*(cc*rr*qp(1)+4.d0*ff*q(1))+4.d0*q(1)*q(1)
     +       *(aa-nn-rro*gr)
      endif
c$$$      vals(3)=rr*qp(1)*(cc*rr*qp(1)+4.d0*ff*q(1))+4.d0*q(1)*q(1)
c$$$     +    *(aa-nn-rro*gr)
      return
      end
      subroutine rps(i,a)
c*** radial mode start soln using sph bessel fns.
      use bits
      use param_modele
      implicit real*8(a-h,o-z)
      save
      dimension a(2)
      fla=flam(i)*(1.d0+xlam(i)*fct)
      sig=fla+2.d0*fmu(i)*(1.d0+qshear(i)*fct)
      zsq=r(i)*r(i)*rho(i)*(wsq+4.d0*g(i)/r(i))/sig
      call bfs(1,zsq,eps,fp)
      a(1)=r(i)
      a(2)=sig*fp+2.d0*fla
      return
      end
      subroutine tps(i,a)
c*** toroidal mode start soln using sph bessel fns.
      use bits
      use param_modele
      implicit real*8(a-h,o-z)
      save
      dimension a(2)
      fu=fmu(i)*(1.d0+qshear(i)*fct)
      zsq=r(i)*r(i)*wsq*rho(i)/fu
      call bfs(l,zsq,eps,fp)
      a(1)=r(i)
      a(2)=fu*(fp-1.d0)
      return
      end
      subroutine bfs(l,xsq,eps,fp)
c  this routine calculates spherical bessel function of the ist kind.
c  fp is equivalent to (r*dj/dr)/j
c  where r is radius and j is the sbf of order l and argument x=k*r
c  the technique employs the continued fraction approach
c  described in w. lentz's article in applied qptics, vol.15, #3, 1976
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
      subroutine remedy(ls)
c    obtains the eigenfunction of an awkward spheroidal mode by
c    integrating to the icb or the mcb.
      use bits
      use eifx
      use arem
      use rindex
      implicit real*8(a-h,o-z)      
      save
      dimension af(4,2),as(6,3),afr(4)
      print 900,ls
  900 format('in remedy with start level : ',i4)
      if(ls.gt.noc) return
      iexp=0
      do  k=1,2
         do  j=1,4
            af(j,k)=0.d0
         enddo
      enddo
      af(1,1)=1.d0
      if(kg.eq.1) af(2,2)=1.d0
      if(nsl.eq.n) goto 5
      do  i=nslp1,n
         do  k=1,3
            do  j=1,6
               a(j,k,i)=0.d0
            enddo
         enddo
      enddo
      call fprop(n,nslp1,af,iexp)
    5 call fsbdry(af,as,kg)
      do  k=1,3
         do  j=1,6
            a(j,k,nsl)=as(j,k)
         enddo
      enddo
      if(n.ne.nsl) call ortho(n,nsl,as,kg)
      call sprop(n,nsl,nocp1,as,iexp)
      call sfbdry(n,nocp1,as,af,kg)
      imtch=noc
      do i=1,4
         afr(i)=ar(i,noc)
      enddo
      if(ls.gt.nic) goto 15
      icomp=0
      call match(n,noc,kg,af,afr,icomp)
      if(icomp.eq.0) return
      call fprop(noc,nicp1,af,iexp)
      imtch=nic
      do  i=1,4
         afr(i)=ar(i,nicp1)
      enddo
   15 icomp=-1
      call match(n,imtch,kg,af,afr,icomp)
      return
      end
      subroutine fprop(jf,jl,f,iexp)
!      use shanks,sdum=>stepf,idum=>maxo 
      use shanks
      use bits
      use eifx
      use arem
      use rindex
      use param_modele
      use yannos_flag
c    fprop propagates the fundamental matrix f from jf to jl (a fluid region)
      implicit real*8(a-h,o-z)
      save
      dimension f(4,2),s(4,2),h(4,2,10)
      data econst/1048576.d0/
      kk=kg+1
      jj=2*kk
      jud=1
      if(jl.lt.jf) jud=-1
      y=r(jf)
      i=jf
      go to 80
   10 x=y
      y=r(i)
      if(y.eq.x) goto 80
      iq=min0(i,i-jud)
      qff=1.d0+xlam(iq)*fct
      zs=dmin1(x,y)
      xi=g(i)/y
      alfsq=(wsq+4.d0*rho(i)+xi-fl3*xi*xi/wsq)*rho(i)/flam(i)
      q=dmax1(sfl3/x,dsqrt(dabs(alfsq-fl3/(x*x)))+1.d0/zs)
      del=jud*step(8)/q
      dxs=0.d0
   15 y=x+del
      if(float(jud)*(y-r(i)).gt.0.d0) y=r(i)
      dx=y-x
      if(dx.ne.dxs) call baylis(q,7)
      dxs=dx
      do 20 k=1,kk
      do 20 j=1,jj
   20 s(j,k)=f(j,k)
      d=fl3/wsq
      do 40 ni=1,in
      z=x+c(ni)
      t=z-zs
      zr=1.d0/z
      ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
      flu=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq))))*qff
      gr=(g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq))))*zr
       if (cancel_gravity) then
          t21=0.0d0
       else
          t21=-4.d0*ro  !a26
       endif
      t12=d*zr*zr
      t11=(gr*d-1.d0)*zr
      s11=-ro*(wsq+4.d0*gr-gr*gr*d)
      c11=-t12/ro+1.d0/flu
      if(kg.eq.0) s11=s11-t21*ro
      if(kg.eq.0) goto 25
      t22=-fl*zr
      s22=ro*t12
      s12=ro*(t11+t22)
   25 do 70 k=1,kk
      if(kg.ne.0) goto 30
      h(1,k,ni)=t11*f(1,k)+c11*f(2,k)
      h(2,k,ni)=s11*f(1,k)-t11*f(2,k)
      goto 35
   30 h(1,k,ni)=t11*f(1,k)+t12*f(2,k)+c11*f(3,k)
      h(2,k,ni)=t21*f(1,k)+t22*f(2,k)+4.d0*f(4,k)
      h(3,k,ni)=s11*f(1,k)+s12*f(2,k)-t11*f(3,k)-t21*f(4,k)
      h(4,k,ni)=s12*f(1,k)+s22*f(2,k)-t12*f(3,k)-t22*f(4,k)
   35 do 70 j=1,jj
      go to (701,702,703,704,705,706,707,708,709,710),ni
  701 f(j,k)=s(j,k)+b(1)*h(j,k,1)
      go to 70
  702 f(j,k)=s(j,k)+b(2)*(h(j,k,1)+b(3)*h(j,k,2))
      go to 70
  703 f(j,k)=s(j,k)+b(4)*(h(j,k,1)+b(5)*h(j,k,2)+b(6)*h(j,k,3))
      go to 70
  704 f(j,k)=s(j,k)+b(7)*(h(j,k,1)+b(8)*h(j,k,2)+b(9)*h(j,k,3)+
     +b(10)*h(j,k,4))
      go to 70
  705 f(j,k)=s(j,k)+b(11)*(h(j,k,1)+b(12)*h(j,k,2)+b(13)*h(j,k,3)+
     +b(14)*h(j,k,4)+b(15)*h(j,k,5))
      go to 70
  706 f(j,k)=s(j,k)+b(16)*(h(j,k,1)+b(17)*h(j,k,2)+b(18)*h(j,k,3)+
     +b(19)*h(j,k,4)+b(20)*h(j,k,5)+b(21)*h(j,k,6))
      go to 70
  707 f(j,k)=s(j,k)+b(22)*(h(j,k,1)+b(23)*h(j,k,3)+b(24)*h(j,k,4)+
     +b(25)*h(j,k,5)+b(26)*h(j,k,6)+b(27)*h(j,k,7))
      go to 70
  708 f(j,k)=s(j,k)+b(28)*(h(j,k,1)+b(29)*h(j,k,3)+b(30)*h(j,k,4)+
     +b(31)*h(j,k,5)+b(32)*h(j,k,6)+b(33)*h(j,k,7)+b(34)*h(j,k,8))
      go to 70
  709 f(j,k)=s(j,k)+b(35)*(h(j,k,1)+b(36)*h(j,k,3)+b(37)*h(j,k,4)+
     +b(38)*h(j,k,5)+b(39)*h(j,k,6)+b(40)*h(j,k,7)+b(41)*h(j,k,8)+
     +b(42)*h(j,k,9))
      go to 70
  710 f(j,k)=s(j,k)+b(43)*(h(j,k,1)+h(j,k,10)+b(45)*h(j,k,5)+
     +b(44)*(h(j,k,4)+h(j,k,6))+b(46)*(h(j,k,7)+h(j,k,9)))
   70 continue
   40 continue
      x=y
      if(y.ne.r(i)) go to 15
   80 size=0.d0
      do 81 k=1,kk
      do 81 j=1,jj
   81 size=dmax1(size,dabs(f(j,k)))
   82 if(size.lt.1024.d0) goto 84
      do 83 k=1,kk
      do 83 j=1,jj
   83 f(j,k)=f(j,k)/econst
      size=size/econst
      iexp=iexp+20
      goto 82
   84 inorm(i)=iexp
      do 85 k=1,kk
      do 85 j=1,jj
   85 a(j,k,i)=f(j,k)
      if(i.eq.jl) return
      i=i+jud
      go to 10
      end
      subroutine sprop(li,jf,jl,f,iexp)
c    sprop propagates the fundamental matrix f from jf to jl (a solid region)
c    if iorth=1 the columns of f are orthogonalized at each level
c    except in regions of oscillatory p and s.
!      use shanks,sdum=>stepf,idum=>maxo 
      use shanks
      use bits
      use eifx
      use arem
      use rindex
      use param_modele
      use yannos_flag
      implicit real*8(a-h,o-z)
      save
      real*8 nn,ll
      dimension f(6,3),s(6,3),h(6,3,10)
      data econst/1048576.d0/
      kk=kg+2
      jj=2*kk
      jud=1
      if(jl.lt.jf) jud=-1
      y=r(jf)
      i=jf
      go to 80
   10 x=y
      y=r(i)
      if(x.eq.y) goto 80
      iq=min0(i,i-jud)
      qff=1.d0+xlam(iq)*fct
      qll=1.d0+qshear(iq)*fct
      qaa=1.d0+xa2(iq)*fct
      zs=dmin1(x,y)
      xi=g(i)/y
      vpsq=(flam(i)+2.d0*fmu(i))/rho(i)
      vssq=fmu(i)/rho(i)
      alfsq=(wsq+4.d0*rho(i)+xi)/vpsq
      betasq=wsq/vssq
      delsq=dsqrt((betasq-alfsq)**2+4.d0*fl3*xi*xi/(vssq*vpsq))
      fksq=.5d0*(alfsq+betasq+delsq)
      al=fl3/(x*x)
      jorth=1
      aq=fksq-delsq-al
      if(aq.gt.0.d0) jorth=0
      qs=dsqrt(dabs(fksq-al))+1.d0/zs
      qf=dsqrt(dabs(aq))+1.d0/zs
      q=dmax1(sfl3/x,qs,qf)
      del=jud*step(8)/q
      dxs=0.d0
   15 y=x+del
      if(float(jud)*(y-r(i)).gt.0.d0) y=r(i)
      dx=y-x
      if(dx.ne.dxs) call baylis(q,7)
      dxs=dx
      do 20 k=1,kk
      do 20 j=1,jj
   20 s(j,k)=f(j,k)
      do 50 ni=1,in
      z=x+c(ni)
      t=z-zs
      ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
      gr=g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq)))
      ff=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq))))*qff
      ll=(lcon(iq)+t*(lspl(1,iq)+t*(lspl(2,iq)+t*lspl(3,iq))))*qll
      if(ifanis.ne.0) goto 25
      nn=ll
      cc=ff+ll+ll
      aa=cc
      goto 30
   25 nn=(ncon(iq)+t*(nspl(1,iq)+t*(nspl(2,iq)+t*nspl(3,iq))))*qll
      cc=(ccon(iq)+t*(cspl(1,iq)+t*(cspl(2,iq)+t*cspl(3,iq))))*qaa
      aa=(acon(iq)+t*(aspl(1,iq)+t*(aspl(2,iq)+t*aspl(3,iq))))*qaa
   30 zr=1.d0/z
      sfl3z=sfl3*zr
      rogr=ro*gr
      c11=1.d0/cc
      c22=1.d0/ll
      dmg=aa-nn-ff*ff*c11
      zdmg=zr*dmg
      t11=-2.d0*ff*zr*c11+zr
      t12=sfl3z*ff*c11
      t21=-sfl3z
      t22=zr+zr
      s22=-ro*wsq
      s11=s22+4.d0*zr*(zdmg-rogr)
      s22=s22+zr*zr*(fl3*(dmg+nn)-nn-nn)
      s12=sfl3z*(rogr-zdmg-zdmg)
      if(kg.eq.0.and..not.cancel_gravity) s11=s11+4.d0*ro*ro
      if(kg.eq.0) goto 35
      t31=-4.d0*ro
      t33=-fl*zr
      s13=-fl1*zr*ro
      s23=ro*sfl3z
   35 do 70 k=1,kk
      if(kg.eq.1) goto 40
      h(1,k,ni)=t11*f(1,k)+t12*f(2,k)+c11*f(3,k)
      h(2,k,ni)=t21*f(1,k)+t22*f(2,k)+c22*f(4,k)
      h(3,k,ni)=s11*f(1,k)+s12*f(2,k)-t11*f(3,k)-t21*f(4,k)
      h(4,k,ni)=s12*f(1,k)+s22*f(2,k)-t12*f(3,k)-t22*f(4,k)
      goto 45
   40 h(1,k,ni)=t11*f(1,k)+t12*f(2,k)+c11*f(4,k)
      h(2,k,ni)=t21*f(1,k)+t22*f(2,k)+c22*f(5,k)
      h(3,k,ni)=t31*f(1,k)+t33*f(3,k)+4.d0*f(6,k)
      h(4,k,ni)=s11*f(1,k)+s12*f(2,k)+s13*f(3,k)-t11*f(4,k)-t21*f(5,k)
     +    -t31*f(6,k)
      h(5,k,ni)=s12*f(1,k)+s22*f(2,k)+s23*f(3,k)-t12*f(4,k)-t22*f(5,k)
      h(6,k,ni)=s13*f(1,k)+s23*f(2,k)-t33*f(6,k)
   45 do 70 j=1,jj
      go to (701,702,703,704,705,706,707,708,709,710),ni
  701 f(j,k)=s(j,k)+b(1)*h(j,k,1)
      go to 70
  702 f(j,k)=s(j,k)+b(2)*(h(j,k,1)+b(3)*h(j,k,2))
      go to 70
  703 f(j,k)=s(j,k)+b(4)*(h(j,k,1)+b(5)*h(j,k,2)+b(6)*h(j,k,3))
      go to 70
  704 f(j,k)=s(j,k)+b(7)*(h(j,k,1)+b(8)*h(j,k,2)+b(9)*h(j,k,3)+
     +b(10)*h(j,k,4))
      go to 70
  705 f(j,k)=s(j,k)+b(11)*(h(j,k,1)+b(12)*h(j,k,2)+b(13)*h(j,k,3)+
     +b(14)*h(j,k,4)+b(15)*h(j,k,5))
      go to 70
  706 f(j,k)=s(j,k)+b(16)*(h(j,k,1)+b(17)*h(j,k,2)+b(18)*h(j,k,3)+
     +b(19)*h(j,k,4)+b(20)*h(j,k,5)+b(21)*h(j,k,6))
      go to 70
  707 f(j,k)=s(j,k)+b(22)*(h(j,k,1)+b(23)*h(j,k,3)+b(24)*h(j,k,4)+
     +b(25)*h(j,k,5)+b(26)*h(j,k,6)+b(27)*h(j,k,7))
      go to 70
  708 f(j,k)=s(j,k)+b(28)*(h(j,k,1)+b(29)*h(j,k,3)+b(30)*h(j,k,4)+
     +b(31)*h(j,k,5)+b(32)*h(j,k,6)+b(33)*h(j,k,7)+b(34)*h(j,k,8))
      go to 70
  709 f(j,k)=s(j,k)+b(35)*(h(j,k,1)+b(36)*h(j,k,3)+b(37)*h(j,k,4)+
     +b(38)*h(j,k,5)+b(39)*h(j,k,6)+b(40)*h(j,k,7)+b(41)*h(j,k,8)+
     +b(42)*h(j,k,9))
      go to 70
  710 f(j,k)=s(j,k)+b(43)*(h(j,k,1)+h(j,k,10)+b(45)*h(j,k,5)+
     +b(44)*(h(j,k,4)+h(j,k,6))+b(46)*(h(j,k,7)+h(j,k,9)))
   70 continue
   50 continue
      x=y
      if(y.ne.r(i)) go to 15
   80 size=0.d0
      do 81 k=1,kk
      do 81 j=1,jj
   81 size=dmax1(size,dabs(f(j,k)))
   82 if(size.lt.1024.d0) goto 84
      do 83 k=1,kk
      do 83 j=1,jj
   83 f(j,k)=f(j,k)/econst
      size=size/econst
      iexp=iexp+20
      goto 82
   84 inorm(i)=iexp
      do 85 k=1,kk
      do 85 j=1,jj
   85 a(j,k,i)=f(j,k)
      if(jorth.eq.1) call ortho(li,i,f,kg)
      if(i.eq.jl) return
      i=i+jud
      go to 10
      end
      subroutine sfbdry(jf,jl,as,af,kg)
c*** the tangential traction scalar is forced to vanish at the solid
c*** side of a s/f boundary(level jl).a(j,3,i) is elliminated for
c*** i=jf...jl and af is loaded from a at level jl.
      use arem
      implicit real*8(a-h,o-z)
      save
      dimension as(6,*),af(4,*)
      n1=min0(jf,jl)
      n2=max0(jf,jl)
      if(kg.ne.0) goto 25
      i1=1
      i2=2
      if(dabs(as(4,2)).gt.dabs(as(4,1))) goto 10
      i1=2
      i2=1
   10 rat=-as(4,i1)/as(4,i2)
      do 15 i=n1,n2
      do 15 j=1,4
   15 a(j,1,i)=a(j,i1,i)+rat*a(j,i2,i)
      af(1,1)=a(1,1,jl)
      af(2,1)=a(3,1,jl)
      return
   25 ab53=dabs(as(5,3))
      do 30 k=1,2
      i1=k
      i2=3
      if(ab53.gt.dabs(as(5,k))) goto 35
      i1=3
      i2=k
   35 rat=-as(5,i1)/as(5,i2)
      do 40 i=n1,n2
      do 40 j=1,6
   40 a(j,k,i)=a(j,i1,i)+rat*a(j,i2,i)
      af(1,k)=a(1,k,jl)
      af(2,k)=a(3,k,jl)
      af(3,k)=a(4,k,jl)
   30 af(4,k)=a(6,k,jl)
      return
      end
      subroutine fsbdry(af,as,kg)
c    fsbdry creates solid fundamental matrix as from fluid fundamental matrix
c    af.it is presumed that fsbdry is used to cross a f/s boundary.
      implicit real*8(a-h,o-z)
      save
      dimension af(4,1),as(6,*)
      do 10 i=1,3
      do 10 j=1,6
   10 as(j,i)=0.d0
      if(kg.ne.0) goto 20
      as(1,1)=af(1,1)
      as(3,1)=af(2,1)
      as(2,2)=1.d0
      return
   20 do 25 k=1,2
      as(1,k)=af(1,k)
      as(3,k)=af(2,k)
      as(4,k)=af(3,k)
   25 as(6,k)=af(4,k)
      as(2,3)=1.d0
      return
      end
      subroutine match(n,j,kg,af,afr,icomp)
      use eifx
      use arem
      implicit real*8(a-h,o-z)
      save
      dimension af(4,*),afr(*),afi(4)
      k=j+2
      rms=0.d0
      fnor=0.d0
      if(kg.eq.1) go to 20
!modif yann?? sait pas trop pourquoi
      ctmp=af(1,1)**2+af(2,1)**2
      if (ctmp>0.0d0) then
         c=(af(1,1)*afr(1)+af(2,1)*afr(2))/ctmp
      else
         c=0.0d0
      endif
!fin modif yann
!      c=(af(1,1)*afr(1)+af(2,1)*afr(2))/(af(1,1)**2+af(2,1)**2)
      do 5 i=1,2
      afi(i)=af(i,1)*c
      rms=rms+(afi(i)-afr(i))**2
    5 fnor=fnor+afr(i)*afr(i)
      rms=dsqrt(rms/fnor)
      if(icomp.lt.0) goto 6
      if(rms.lt.1.d-3) goto 6
      icomp=1
      return
    6 idiff=inorm(j)-inorm(j+1)
      inorm(j+1)=inorm(j)
  999 format(4g20.10)
      do 10 i=k,n
      inorm(i)=inorm(i)+idiff
      do 10 jj=1,4
   10 ar(jj,i)=c*a(jj,1,i)
      return
   20 continue
      ctmp=(af(1,2)*af(3,1)-af(1,1)*af(3,2))
      if (ctmp==0.0d0) ctmp=1.d0
      a2=(af(3,1)*afr(1)-af(1,1)*afr(3))/ctmp
      ctmp=(af(1,1)*af(3,2)-af(3,1)*af(1,2))
      if (ctmp==0.0d0) ctmp=1.d0
      a1=(af(3,2)*afr(1)-af(1,2)*afr(3))/ctmp
      do 21 i=1,4
      afi(i)=a1*af(i,1)+a2*af(i,2)
      rms=rms+(afi(i)-afr(i))**2
   21 fnor=fnor+afr(i)*afr(i)
      rms=dsqrt(rms/fnor)
c      print 999,rms
      if(icomp.lt.0) goto 22
      if(rms.lt.1.d-3) goto 22
      icomp=1
      return
   22 idiff=inorm(j)-inorm(j+1)
      inorm(j+1)=inorm(j)
      do 25 i=k,n
      inorm(i)=inorm(i)+idiff
      do 25 jj=1,6
   25 ar(jj,i)=a1*a(jj,1,i)+a2*a(jj,2,i)
      return
      end
      subroutine ortho(li,lc,b,kg)
c    finds the orthogonal matrix v such that the columns of b*v are orthogonal
c   the array a is replaced by a*v for levels li - lc. array b is replaced
c   by b*v and is then ready fo entry to sprop at level lc.this is intended
c   to diminish the onset of degeneracy caused by rapid exponential growth
c    in the mantle for modes with deeply turning s and shallowly turning p.
      use arem
      implicit real*8(a-h,o-z)
      save
      dimension b(6,1),as(6,3)
      i1=min0(lc,li)
      i2=max0(lc,li)
      nc=kg+2
      nr=2*nc
      call svd(b,nr,nc)
      do 25 i=i1,i2
      do 20 j=1,nc
      do 20 k=1,nr
      as(k,j)=0.d0
      do 20 l=1,nc
   20 as(k,j)=as(k,j)+a(k,l,i)*b(l,j)
      do 25 j=1,nc
      do 25 k=1,nr
   25 a(k,j,i)=as(k,j)
      do 35 j=1,nc
      do 35 k=1,nr
   35 b(k,j)=a(k,j,lc)
      return
      end
      subroutine svd(a,mrow,ncol)
c    section i chapter 10 wilkenson and reinsch (1971 ,springer).
c    the matrix a is overwritten with v(ncol,ncol), the right side orthogonal
c    matrix in the svd decomposition. for use only in eos subs as ,to reduce
c    branching points, i have used the fact that ncol is lt mrow.
      implicit real*8(a-h,o-z)
      save
      dimension a(6,1),e(3),q(3)
      eps=1.5d-14
      tol=1.d-293
      g=0.d0
      x=0.d0
      do 60 i=1,ncol
      l=i+1
      e(i)=g
      s=0.d0
      do 10 j=i,mrow
   10 s=s+a(j,i)*a(j,i)
      if(s.gt.tol) go to 15
      q(i)=0.d0
      if(l.gt.ncol) goto 60
      go to 30
   15 q(i)=dsign(dsqrt(s),-a(i,i))
      h=a(i,i)*q(i)-s
      a(i,i)=a(i,i)-q(i)
      if(l.gt.ncol) go to 60
      do 25 j=l,ncol
      s=0.d0
      do 20 k=i,mrow
   20 s=s+a(k,i)*a(k,j)
      f=s/h
      do 25 k=i,mrow
   25 a(k,j)=a(k,j)+f*a(k,i)
   30 s=0.d0
      do 35 j=l,ncol
   35 s=s+a(i,j)*a(i,j)
      if(s.ge.tol)go to 40
      g=0.d0
      go to 60
   40 g=dsign(dsqrt(s),-a(i,l))
      h=a(i,l)*g-s
      a(i,l)=a(i,l)-g
      do 45 j=l,ncol
   45 e(j)=a(i,j)/h
      do 55 j=l,mrow
      s=0.d0
      do 50 k=l,ncol
   50 s=s+a(j,k)*a(i,k)
      do 55 k=l,ncol
   55 a(j,k)=a(j,k)+s*e(k)
   60 x=dmax1(dabs(q(i))+dabs(e(i)),x)
      goto 100
   75 if(g.eq.0.d0)go to 91
      h=a(i,l)*g
      do 80 j=l,ncol
   80 a(j,i)=a(i,j)/h
      do 90 j=l,ncol
      s=0.d0
      do 85 k=l,ncol
   85 s=s+a(i,k)*a(k,j)
      do 90 k=l,ncol
   90 a(k,j)=a(k,j)+s*a(k,i)
   91 do 95 j=l,ncol
      a(i,j)=0.d0
   95 a(j,i)=0.d0
  100 a(i,i)=1.d0
      g=e(i)
      l=i
      i=i-1
      if(i.ge.1)go to 75
      ep=eps*x
      k=ncol
  105 l=k
  110 if(dabs(e(l)).le.ep)go to 125
      if(dabs(q(l-1)).le.ep) go to 115
      l=l-1
      if(l.ge.1)go to 110
  115 c=0.d0
      s=1.d0
      do 120 i=l,k
      f=s*e(i)
      e(i)=c*e(i)
      if(dabs(f).le.ep)go to 125
      g=q(i)
      h=dsqrt(f*f+g*g)
      c=g/h
      s=-f/h
  120 q(i)=h
  125 z=q(k)
      if(l.eq.k)go to 145
      x=q(l)
      y=q(k-1)
      g=e(k-1)
      h=e(k)
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.d0*h*y)
      g=dsqrt(f*f+1.d0)
      f=((x-z)*(x+z)+h*(y/(f+dsign(g,f))-h))/x
      c=1.d0
      s=1.d0
      lp1=l+1
      do 140 i=lp1,k
      g=e(i)
      y=q(i)
      h=s*g
      g=c*g
      z=dsqrt(f*f+h*h)
      im1=i-1
      e(im1)=z
      c=f/z
      s=h/z
      f=s*g+c*x
      g=c*g-s*x
      h=s*y
      y=c*y
      do 130 j=1,ncol
      x=a(j,im1)
      z=a(j,i)
      a(j,im1)=c*x+s*z
  130 a(j,i)=c*z-s*x
      z=dsqrt(f*f+h*h)
      q(im1)=z
      c=f/z
      s=h/z
      f=s*y+c*g
  140 x=c*y-s*g
      e(l)=0.d0
      e(k)=f
      q(k)=x
      go to 105
  145 if(z.ge.0.d0)go to 155
      q(k)=-z
      do 150 j=1,ncol
  150 a(j,k)=-a(j,k)
  155 k=k-1
      if(k.ge.1)go to 105
      return
      end
c$$$c****************************************************
c$$$      subroutine zknt(s,sp,f,fp,x,y,ifsol,ibid)
c$$$      use bits
c$$$c*** given minor vector and derivs,constructs mode count ***
c$$$      implicit real*8(a-h,o-z)
c$$$      save
c$$$c      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
c$$$c     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
c$$$      dimension s(*),sp(*),f(*),fp(*),xs(4),val(4)
c$$$      if(ifsol.eq.0.and.kg.eq.0) goto 5
c$$$      y1=s(5)
c$$$      y2=f(5)
c$$$      y1p=sp(5)
c$$$      y2p=fp(5)
c$$$      t1=s(3)-s(4)
c$$$      t2=f(3)-f(4)
c$$$      t1p=sp(3)-sp(4)
c$$$      t2p=fp(3)-fp(4)
c$$$      goto 10
c$$$    5 y1=s(2)
c$$$      y2=f(2)
c$$$      y1p=sp(2)
c$$$      y2p=fp(2)
c$$$      t1=s(1)
c$$$      t2=f(1)
c$$$      t1p=sp(1)
c$$$      t2p=fp(1)
c$$$   10 h=y-x
c$$$      ns=0
c$$$      if(kount.ne.0) goto 15
c$$$      a1=y2-y1
c$$$      a2=0.d0
c$$$      a3=0.d0
c$$$      a22=0.d0
c$$$      a33=0.d0
c$$$      goto 50
c$$$   15 a1=h*y1p
c$$$      a2=-h*(2.d0*y1p+y2p)+3.d0*(y2-y1)
c$$$      a3=h*(y1p+y2p)-2.d0*(y2-y1)
c$$$      a33=3.d0*a3
c$$$      a22=2.d0*a2
c$$$      if(a3.ne.0.d0) goto 20
c$$$      if(a2.eq.0.d0) goto 50
c$$$      xs(2)=-a1/a22
c$$$      if(xs(2).ge.0.d0.and.xs(2).le.1.d0) ns=1
c$$$      goto 50
c$$$   20 disc=a2*a2-a1*a33
c$$$      if(disc) 50,25,30
c$$$   25 xs(2)=-a2/a33
c$$$      if(xs(2).ge.0.d0.and.xs(2).le.1.d0) ns=1
c$$$      goto 50
c$$$   30 disc=dsqrt(disc)
c$$$      tr1=(-a2+disc)/a33
c$$$      tr2=(-a2-disc)/a33
c$$$      if(dabs(a33).gt.dabs(a1)) goto 35
c$$$      fac=a1/a33
c$$$      tr1=fac/tr1
c$$$      tr2=fac/tr2
c$$$   35 if(tr1.lt.0.d0.or.tr1.gt.1.d0) goto 40
c$$$      xs(2)=tr1
c$$$      ns=1
c$$$   40 if(tr2.lt.0.d0.or.tr2.gt.1.d0) goto 50
c$$$      ns=ns+1
c$$$      xs(ns+1)=tr2
c$$$      if(ns.lt.2) goto 50
c$$$      if(tr2.ge.tr1) goto 50
c$$$      xs(2)=tr2
c$$$      xs(3)=tr1
c$$$   50 val(1)=y1
c$$$      xs(1)=0.d0
c$$$      ns2=ns+2
c$$$      val(ns2)=y2
c$$$      xs(ns2)=1.d0
c$$$      if(ns.eq.0) goto 60
c$$$      ns1=ns+1
c$$$      do 55 j=2,ns1
c$$$      t=xs(j)
c$$$   55 val(j)=y1+t*(a1+t*(a2+t*a3))
c$$$   60 ift=0
c$$$      do 100 j=2,ns2
c$$$      if(val(j-1)*val(j).gt.0.d0) goto 100
c$$$      if(val(j-1).ne.0.d0) goto 65
c$$$      tes=t1*a1
c$$$      goto 90
c$$$   65 rt1=0.5d0*(xs(j-1)+xs(j))
c$$$      rt=rt1
c$$$      do 70 i=1,5
c$$$      v=y1+rt*(a1+rt*(a2+rt*a3))
c$$$      vp=a1+rt*(a22+rt*a33)
c$$$      add=-v/vp
c$$$      rt=rt+add
c$$$      if(dabs(add).lt.1.d-5) goto 75
c$$$      if(dabs(rt-rt1).le..5d0) goto 70
c$$$      rt=rt1
c$$$      goto 75
c$$$   70 continue
c$$$   75 if(ift.ne.0) goto 85
c$$$      if(kount.ne.0) goto 80
c$$$      b1=t2-t1
c$$$      b2=0.d0
c$$$      b3=0.d0
c$$$      goto 85
c$$$   80 b1=h*t1p
c$$$      b2=-h*(2.d0*t1p+t2p)+3.d0*(t2-t1)
c$$$      b3=h*(t1p+t2p)-2.d0*(t2-t1)
c$$$      ift=1
c$$$   85 tes=t1+rt*(b1+rt*(b2+rt*b3))
c$$$      vp=a1+rt*(a22+rt*a33)
c$$$      tes=tes*vp
c$$$   90 if(tes.lt.0.d0) kount=1+kount
c$$$      if(tes.gt.0.d0) kount=kount-1
c$$$  100 continue
c$$$      return
c$$$      end
