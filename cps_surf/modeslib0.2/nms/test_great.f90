program test_great
  implicit none
  real*8 :: t1,t2,p1,p2,pi,conv,la1,lo1,la2,lo2,d,dd,g,a12,a21,a,b
  real*8, parameter ::  GEOC=0.993277d0
  pi = 3.141592653589793d0
  conv=pi/180.d0
!
  t1=90- 0003.30
  t2=90- 039.869
  p1= 0095.94
  p2=032.794
!
  la1=90.-t1
  lo1=p1
  la2=90.-t2
  lo2=p2
!
  call great(la1,lo1,la2,lo2,d,dd,g,a12,a21)
  print*,dd,a12,a21
!correction
  t1=pi/2-atan(GEOC*tan((pi/2-t1*conv)))
  t2=pi/2-atan(GEOC*tan((pi/2-t2*conv)))
  p1=p1*conv
  p2=p2*conv
  call euler(t1,p1,t2,p2,a,b,g)
  g=pi-g
  if (a<0) a=a+2*pi
  if (g<0) g=g+2*pi
  print*,b/conv,g/conv,a/conv

end program test_great

!----------------------------------------------------------------
      subroutine great(alat1,alon1,alat2,alon2,dist,disd,gc, &
       az12,az21)
! calculation of distance in km and in deg, length of great circle path,
! and azimuth
!  input (lat1,lon1) (lat2,lon2)
!  output	odist : distance in km
!		odisd : distance in degree
!		ogc   : length of great circle path
!		oz12  : azimuth of 2 at 1 measured clockwise from north
!		oz21  : azimuth of 1 at 2
!__________________________________________________________________
      implicit real*8 (a-h,o-z)
       ath=6378.140
      bth=6356.755
      pi = 3.141592653589793d0
      rad = pi/180.
      h = 1. - bth*bth/(ath*ath)
      p = h/(1. - h)
!print*,h
!p=1-0.993277
!p=0.d0
!h=p/(p+1)
!print*,h
      gr = alon1*rad
      tr = alat1*rad
      sintr =dsin(tr)
      costr =dcos(tr)
      if (sintr .eq. 0.) sintr = .00000100
      if (costr .eq. 0.) costr = .00000100
      r1 = ath/dsqrt(1. - h*sintr*sintr)
      z1 = r1*(1. - h)*sintr
      g = alon2*rad
      t = alat2*rad
      if (t .eq. 0.) t = .0000100
      sint =dsin(t)
      cost =dcos(t)
      r2 = ath/dsqrt(1. - h*sint*sint)
      dg = g - gr
      cosdg =dcos(dg)
      sindg =dsin(dg)
      dgr = gr - g
      dt = t - tr
      q = sint*costr/((1. + p)*cost*sintr) + h*r1*costr/(r2*cost)
      x = r2*cost*cosdg
      y = r2*cost*sindg
      z = r2*(1. - h)*sint
      az12 =datan2(sindg,(q - cosdg)*sintr)
      q = sintr*cost/(costr*sint*(1. + p)) + h*r2*cost/(r1*costr)
      az21 = datan2(dsin(dgr),sint*(q-dcos(dgr)))
      cos12 =dcos(az12)
      cta2 = costr*costr*cos12*cos12
      p0 = p*(cta2 + sintr*sintr)
      b0 = (r1/(1. + p0))*dsqrt(1. + p*cta2)
      e0 = p0/(1. + p0)
      gc = 2.*pi*b0*dsqrt(1. + p0)*(1. - e0*(.25 + e0*(3./64. &
                                               + 5.*e0/256.)))
      c0 = 1. + p0*(.25 - p0*(3./64. - 5.*p0/256.))
      c2 = p0*(-.125 + p0*(1./32. - 15.*p0/1024.))
      c4 = (-1./256. + 3.*p0/1024.)*p0*p0
      u0 =datan2(sintr,costr*cos12*dsqrt(1. + p0))
      u =datan2(r1*sintr + (1. + p0)*(z - z1),(x*cos12 - y*sintr* &
                                            dsin(az12))*dsqrt(1. + p0))
      disd = u - u0
      if (u .lt. u0) disd = pi + pi + disd
      dist = b0*(c0*( disd ) +c2*(dsin(u + u) -dsin(u0 + u0)) &
                            +c4*(dsin(4.*u) -dsin(4.*u0)))
      disd = disd/rad
      az12 = az12/rad
      az21 = az21/rad
      if (az12 .lt. 0.) az12 = 360. + az12
      if (az21 .lt. 0.) az21 = 360. + az21
    end subroutine great
!------------------------------------------------------------
      subroutine euler(t1,p1,t2,p2,a,b,g)
!_____________________________________________________________________
!   
!     calcule les 3 angles d'euler de la rotation definie par:
!     D(a b g)=D(-p1 -t1 0)oD(0 t2 p2)
!     (notation cf Edmonds)
!
!     Entrees:
!        angles p1 p2 dans [-pi,pi] (en fait peu importe mais en
!                                    radians)
!        angles t1 t2 DANS [ 0 ,pi]! (pour etre avec utilise avec 
!        la regle  de sommation d'harmoniques spheriques)
!     Sorties:
!        angles a et g dans [-pi,pi]
!        angle  b      dans [ 0 ,pi]
!_____________________________________________________________________
!
      IMPLICIT NONE
!
      real*8 t1,t2,p1,p2,a,b,g,ca,sa,cg,sg,sb,cb,pi
      real*8 rt1(3,3),rt2(3,3),rp1(3,3),rp2(3,3),prod1(3,3), &
           prod2(3,3),rot(3,3)
      integer i,j,k
!
      pi = 3.141592653589793d0
      do i=1,3
         do j=1,3
            rt1(i,j)  =0.0d0
            rt2(i,j)  =0.0d0
            rp1(i,j)  =0.0d0
            rp2(i,j)  =0.0d0
            prod1(i,j)=0.0d0
            prod2(i,j)=0.0d0
            rot(i,j)  =0.0d0
         enddo
      enddo
!      
      rt1(1,1)=dcos(t1)
      rt1(2,2)=1.0d0
      rt1(3,3)=rt1(1,1)
      rt1(1,3)=dsin(t1)
      rt1(3,1)=-rt1(1,3)
!
      rt2(1,1)=dcos(t2)
      rt2(2,2)=1.0d0
      rt2(3,3)=rt2(1,1)
      rt2(1,3)=-dsin(t2)
      rt2(3,1)=-rt2(1,3)
!
      rp1(1,1)=dcos(p1)
      rp1(2,2)=rp1(1,1)
      rp1(3,3)=1.0d0
      rp1(1,2)=-dsin(p1)
      rp1(2,1)=-rp1(1,2)
!
      rp2(1,1)=dcos(p2)
      rp2(2,2)=rp2(1,1)
      rp2(3,3)=1.0d0
      rp2(1,2)=dsin(p2)
      rp2(2,1)=-rp2(1,2)
!
      do i=1,3
         do j=1,3
            do k=1,3
               prod1(i,j)=prod1(i,j)+rp2(i,k)*rt2(k,j)
            enddo
         enddo
      enddo
      do i=1,3
         do j=1,3
            do k=1,3
               prod2(i,j)=prod2(i,j)+rp1(i,k)*prod1(k,j)
            enddo
         enddo
      enddo
      do i=1,3
         do j=1,3
            do k=1,3
               rot(i,j)=rot(i,j)+rt1(i,k)*prod2(k,j)
            enddo
         enddo
      enddo
!
      cb=rot(3,3)
      if (cb>=1.d0) then
         b=0.d0
      else if (cb<=-1.d0) then
         b=pi
      else
         b=dacos(cb)
      endif
      sb=dsin(b)
      if (abs(sb).le.1.0d-15) then
         a=p2-p1
         g=0.
      else 
         ca=rot(3,1)/sb
         sa=rot(3,2)/sb
         if (abs(ca-1.0d0).lt.1.0d-8) then
            a=0.0d0
         else 
            if (abs(ca+1.0d0).lt.1.0d-8) then
               a=pi
            else
               a=dacos(ca)
            endif
         endif
         if (sa.lt.(0.0d0)) a=-1.0d0*a
         cg=-rot(1,3)/sb
         sg=rot(2,3)/sb
         if (abs(cg-1.0d0).lt.1.0d-8) then
            g=0.0d0
         else
            if (abs(cg+1.0d0).lt.1.0d-8) then
               g=pi
            else
               g=dacos(cg)
            endif
         endif
         if (sg.lt.(0.0d0)) g=-1.0d0*g
      endif
!
      end subroutine euler
