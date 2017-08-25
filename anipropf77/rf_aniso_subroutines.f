c  rfgen - program to calculate transmission response of 
c  a stack of anisotropic layers to a plane wave incident from below,
c  and to compute P-S "receiver function".
c  built upon anirec.f 01/96, 04/96 VL
c  binary I/O taken out to facilitate use on Linux 09/99 VL
c  dynamic assignment of frequency sampling implemented
c  to make the code compatible with "recfunk" data analysis
c  software. 10/2001 VL
c
c  refft.o is an object file of the C-routine refft.c
c  compile this code with 
c      cc -c refft.c -o refft.o
c  f77 -o rfsyn rfsyn.f refft.o
c
c  for hexagonally symmetric media
c  reads fast axis orientation, constants A,B,C,D,E from file animodel
c  calculate quadratic eigenvalue problem based on the Christoffel matrix
c  see appendix of P. Shearer's thesis
c
c  read model, phase velocity of incident wave, P, SV, or SH
c
c  calc the eigenvector decomps for the layers
c  loop over frequency, calc reflection/transmission matrices
c  calc 3-comp transfer fct response at surface
c  find distortion of reference wavelet

        
      subroutine rf_aniso_subroutines(z_in,vp_in,vp2_in,vp4_in,
     x           vs_in,vs2_in,rho_in,theta_in,phig_in,nl,baz,
     x           Rrecfun,Trecfun,Timesave,Ntime) 

      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      
      real*8 baz
      integer nl,Ntime
      real*8 z_in(nl+1),vp_in(nl+1),vp2_in(nl+1),vp4_in(nl+1)
      real*8 vs_in(nl+1),vs2_in(nl+1),rho_in(nl+1)
      real*8 theta_in(nl+1),phig_in(nl+1)
      real*8 Rrecfun(Ntime),Trecfun(Ntime),Timesave(Ntime)
      
      complex*16 pp,u0,ee,z1,z0,xnu,e1,e2,zla,xl
      complex*16  rt,tt,rt0,trc,pfac,u,resp,zz
      real*4 cc4,cc5,zz4,frqq,amn,amx,amx1,amn1    
      common/stfff/w(3,101),t(3,3),ttl(3,3),s(3,3),stl(3,3),
     x             r(3,3),x(3),y(3)
      common/model/z(100),dz(100),rho(101),vp(101),vs(101),vp2(101),
     x             vp4(101),vs2(101),vss(101)
      common/model2/xmu(101),xla(101),xmu2(101),xla2(101),xla4(101)
      common/propag/xnu(6,101),xl(6,100),pfac(6,3),u(3,6)
      common/mstff/qq(6,6),wr(6),wi(6),zr(6,6),zi(6,6),iv(6),fv(6)
      common/pstff/pp(3),u0(3),ee(6,6,101),e1(6,6),e2(6,6),zla(6)
      common/rstff/rt(3,3,101),tt(3,3,101),rt0(3,3),trc(3,3)
      common/nstff/cc4(8200),zz4(8200),dat4(8200,3,3),ccc(8200),
     x		   cc5(8200)
      common/disper/resp(3,3,1024),frqq(1024)
      common/disper2/roota(101),rootb(101),jtrval(101),kroots(1024)
      common/evanes/ievan(10000)
      data pi/3.14159265358979d0/,eps/1.d-6/,tol/1.d-3/
      character*20 filein, outnm
c
c  we reduce the condition numbers of matrices by 
c  normalizing physical quantities to make them dimensionless
c  Ill use the normal-mode normalizations
c  which are a little peculiar for the crust, but what the hey!
      rbar=5.515d3
      ren=1.075190645d-3
      radi=6.371d6
      vbar=ren*radi
      con=rbar*vbar**2
      z1=dcmplx(1.d0,0.d0)
      z0=dcmplx(0.d0,0.d0)
c     read the frequency you want to work with,
c     and decide how many frequency points you need.
      nfrq=512
      gainfac = 16
c  gainfac above will compensate timeseries values for 
c the effects of spectral smoothing.
c  it includes a factor of 2 to remove the effect of the cos**2 taper
c and a halved ratio of freq. series length to the number of frequencies
c for nfrq=512 and npad = 8192 gainfac= 8

c      print *, 'enter maximum frequency'
c      read(*,*), frqmax
        frqmax=1.0
       if (frqmax.gt.1) then
	        nfrq=1024
	        gainfac=8
      endif

c      print *, 'frqmax', frqmax, 'nfrq', nfrq

c	read backazimuth (this only makes sence for anisotropic models)

 102  format(a)
      
      nlp=nl+1
      nlm=nl-1
      do i=1,nlp
c angle conventions for anisotropic tensor orientations
c numbers supplied in the input file are :
c phig - azimuth of symmetry axis, theta - tilt from vertical
c internal computations use phi measured
c counterclockwise from X (=radial=baz-180)

        theta=theta_in(i)
        phig=phig_in(i)

c-----correct for backazimuth
c        phi = baz-180-phig
c change to:		12/13/96 BUG DISCOVERED THE NIGHT BEFORE AGU...
        phi = phig-baz
c  end Change
        w(1,i)=dsind(theta)*dcosd(phi)
        w(2,i)=dsind(theta)*dsind(phi)
        w(3,i)=dcosd(theta)
c  read depth to ith interface, vp (m/sec), pk-to-pk cos(2th) relative P pert
c  pk-to-pk cos(4th) relative P pert, v_s,  pk-to-pk cos(2th) relative S pert
c  density (kg/m**3)
        z(i)=z_in(i)
        vp(i)=vp_in(i)
        vp2(i)=vp2_in(i)
        vp4(i)=vp4_in(i)
        vs(i)=vs_in(i)
        vs2(i)=vs2_in(i)
        rho(i)=rho_in(i)

c  recall that we interpret fractional values of b,c,e 
c  as peak-to-peak relative velocity perts.
c  therefore, e=0.02 is 2% pert to mu from slowest to fastest
        xmu(i)=rho(i)*vs(i)**2/con
        xmu2(i)=vs2(i)*xmu(i)
        xla(i)=rho(i)*vp(i)**2/con
        xla2(i)=vp2(i)*xla(i)
        xla4(i)=vp4(i)*xla(i)
        vs(i)=vs(i)/vbar
        vp(i)=vp(i)/vbar
        rho(i)=rho(i)/rbar
        z(i)=z(i)/radi
      end do
      close(7)
      do i=2,nl
        dz(i)=z(i)-z(i-1)
      end do
      dz(1)=z(1)
c  print the organ-pipe mode count for 1Hz
c  the lowest layer (nl+1) is taken as evanescent region.
      sn=0.
      pn=0.
      do i=1,nl
        sn=sn+(dz(i)/vs(i))/ren
        pn=pn+(dz(i)/vp(i))/ren
      end do
c      print *, 'organ-pipe mode count at 1 Hz in stack of layers: S & P'
c      print 104,sn,pn
c  104 format(2f10.1)
c  search for cmin, cmax is halfspace velocity
      cmin=vs(1)
      vss(1)=vs(1)
      do i=2,nlp
        if(cmin.gt.vs(i)) cmin=vs(i)
        vss(i)=vs(i)
      end do
  900 csmin=vs(nlp)*vbar/1000.
      cpmin=vp(nlp)*vbar/1000.
c      print *,'minimum phase velocity for S wave (km/sec)',csmin
c      print *,'minimum phase velocity for P wave (km/sec)',cpmin
c      print *,'enter phase velocity of incident wave (km/sec)'
c      read(*,*) cc
      cc=12.6  

      if(cc.le.0.d0) go to 950
      cc=cc*1000./vbar
c  calc the eigenvector decomps for the layers
c  need to identify upgoing P,SV,SH in the halfspace
      call matget(nl,cc)
c  loop over frequency, calc reflection/transmission matrices
c  calc 3-comp transfer fct response at surface
      df=frqmax/nfrq
      do jf=1,nfrq
        om=2.d0*pi*jf*df/ren
        frqq(jf)=jf*df
        call respget(nl,om,cc,resp(1,1,jf))
      end do

      npad=8192
      dt=1.d0/(npad*df)

      do i=1,npad
        zz4(i)=0.d0
      end do

c  zero the DC and Nyquist 
c  ick switches the sign of y and z components to transverse & vertical
      zz4(1)=0.
      zz4(2)=0.
c     do itype=1,3
c       iprint = 1
c -----here model only P-wave arrival    VL 01/96
c       ick=1
c       baseline=0.
c       do k=1,3
          cc4(1)=0.
          cc4(2)=0.
          cc5(1)=0.
          cc5(2)=0.
          cc4(1)=0.
          cc4(2)=0.

c       open(10,file = "RRF.aspec")
c        open(11,file = "fac")

          do jf=1,nfrq
c    figure out the value of the cosine-squared scaler
         	fac=cos(0.5d0*pi*(jf-1)/nfrq)**2
c
c           zz=ick*dcmplx(dble(zz4(2*jf+1)),dble(zz4(2*jf+2)))
c           zz=zz*resp(k,itype,jf)
	    zz= -resp(1,1,jf) / resp(3,1,jf)
            cc4(2*jf+1)=dreal(zz)*fac
            cc4(2*jf+2)=dimag(zz)*fac
c
c	ampsp =  cc4(2*jf+1)**2 + cc4(2*jf+2)**2
c	write (10,332) jf*df, ampsp
c	write (11,332) jf*df, fac
c	cc4 has spectra of P-SV reciever function
	    zz= resp(2,1,jf) / resp(3,1,jf)
            cc5(2*jf+1)=dreal(zz)*fac
            cc5(2*jf+2)=dimag(zz)*fac

c	ampsp = cc5(2*jf+1)**2 + cc5(2*jf+2)**2
c	write (11,332) jf*df, ampsp

c	cc5 has a spectra of P-SH reciever function
          end do


c	close(10)
c	close(11)
c  332   format(f10.4,1x,f11.8)

          do jf=2*nfrq+3,npad
            cc4(jf)=0.
	    cc5(jf)=0.
          end do


          call refft(cc4,npad,-1,-1)
          call refft(cc5,npad,-1,-1)
c	cc4 and cc5 have time domain reciever functions
c
c      now mult by 2 since cosine**2 taper is applied to spectral RF
c      also mult by gainfac to compensate for losing high freqs

	 do i=1,npad
	   cc4(i)=cc4(i)*gainfac
	   cc5(i)=cc5(i)*gainfac
	 end do
          amx=cc4(1)
          amn=cc4(1)
          do i=1,npad
            amx=amax1(amx,cc4(i))  
          end do
          amx1=cc5(1)
          amn1=cc5(1)
          do i=1,npad
            amx1=amax1(amx1,cc5(i))  
          end do
c	define new trace length - don't need a long tail of zeroes and junk
c	output only first 100 seconds
c	 npad = 1600 
c	output RFs in ascii form, appropriate for use, e.g., with 
c	GMT command pswiggle

c        open(10,file = "RRF.table")
c        open(11,file = "TRF.table")
c        do i=1,npad
c        timestep = (i-1)*dt
c        write (10,333) timestep, baz, cc4(i)
c	write (11,333) timestep, baz, cc5(i)
c        end do
c  333   format(f10.4,1x,f8.3,1x,f11.8)
	 

        ! output result 
        do i=1,Ntime
            Timesave(i)=(i-1)*dt 
            Rrecfun(i)=cc4(i)
            Trecfun(i)=cc5(i)
        end do 

  950 continue
c      stop
      end subroutine rf_aniso_subroutines
c
c -----------END OF MAIN CODE
c
c
      subroutine matget(nl,cc)
c  returns stress-displacement vectors for a stack of anisotropic layers
c  P waves may be evanescent, but S waves are oscillatory in the stack
c  the weirdness seen in the surface wave code should not appear in 
c  a receiver function code
c  however, the iev parameter is retained to avoid leaving timebombs 
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      complex*16 pp,u0,ee,pw,uw,pu,z1,z0,zz,xnu,eye,e1,e2,zla,rtm
      complex*16 pfac,u,xl
      common/stfff/w(3,101),t(3,3),ttl(3,3),s(3,3),stl(3,3),
     x                     			 r(3,3),x(3),y(3)
      common/model/z(100),dz(100),rho(101),vp(101),vs(101),vp2(101),
     x               vp4(101),vs2(101),vss(101)
      common/model2/xmu(101),xla(101),xmu2(101),xla2(101),xla4(101)
      common/propag/xnu(6,101),xl(6,100),pfac(6,3),u(3,6)
      common/defect/idfct(4,101),adf(2,101)
      common/rrt/rtm(6,6,100)
      common/mstff/qq(6,6),wr(6),wi(6),zr(6,6),zi(6,6),iv(6),fv(6)
      common/pstff/pp(3),u0(3),ee(6,6,101),e1(6,6),e2(6,6),zla(6)
      common/qstff/qi(6,6),xr(6),xi(6),yr(6),yi(6),ips(3)
      data pi/3.14159265358979d0/,eps/1.d-6/,tol/1.d-7/
c  set iev=1 ** should be superfluous, save for now
c  toggle to iev=0 if there is a purely propagating wave in the top layer n=1
      iev=1
      z1=dcmplx(1.d0,0.d0)
      z0=dcmplx(0.d0,0.d0)
      eye=dcmplx(0.d0,1.d0)
      rbar=5.515d3
      ren=1.075190645d-3
      radi=6.371d6
      vbar=radi*ren
      con=rbar*radi*radi*ren*ren
      nlp=nl+1
      nlm=nl-1
c  first calculate vertical wavenumbers and propagating waves for each layer
c   requires an eigenvector problem be solved
c  in general, the evanescent vertical wavenumbers have nonzero real parts 
c   complex exponential fct is used to avoid endless branching
c  horizontal slowness p_x
      px=1.d0/cc
      do n=1,nlp
        a=xla(n)
        b=xla2(n)
        c=xla4(n)
        d=xmu(n)
        e=xmu2(n)
c        print *,'a,b,c,d,e',a,b,c,d,e
        fact=8.d0*w(1,n)*w(1,n)*c+2.d0*e
        facs=16.d0*w(1,n)*w(3,n)*c
        facr=8.d0*w(3,n)*w(3,n)*c+2.d0*e
        do i=1,3
c  first the what-0-what tensor
          do j=1,3
            t(j,i)=fact*w(j,n)*w(i,n)
            s(j,i)=facs*w(j,n)*w(i,n)
            r(j,i)=facr*w(j,n)*w(i,n)
          end do
c  next the identity tensor - correct an error on 7/6/95
          t(i,i)=t(i,i)+d+e*(2.d0*w(1,n)*w(1,n)-1.d0)
          s(i,i)=s(i,i)+4.d0*e*w(1,n)*w(3,n)
          r(i,i)=r(i,i)+d+e*(2.d0*w(3,n)*w(3,n)-1.d0)
        end do 
c        print 101,(w(i,n),i=1,3)
c        print 101,fact,facs,facr
c        print *,'t,s,r'
c        print 101,((t(i,j),j=1,3),i=1,3)
c        print 101,((s(i,j),j=1,3),i=1,3)
c        print 101,((r(i,j),j=1,3),i=1,3)
        fac=b-4.d0*c-2.d0*e
c  next the what-0-xhat and what-0-zhat tensors
        do i=1,3
          t(1,i)=t(1,i)+fac*w(1,n)*w(i,n)
          t(i,1)=t(i,1)+fac*w(1,n)*w(i,n)
          s(1,i)=s(1,i)+fac*w(3,n)*w(i,n)
          s(i,1)=s(i,1)+fac*w(3,n)*w(i,n)
          s(3,i)=s(3,i)+fac*w(1,n)*w(i,n)
          s(i,3)=s(i,3)+fac*w(1,n)*w(i,n)
          r(3,i)=r(3,i)+fac*w(3,n)*w(i,n)
          r(i,3)=r(i,3)+fac*w(3,n)*w(i,n)
        end do
        fac=a-b+c-d+e
c  finally the xhat-0-xhat, zhat-0-zhat, xhat-0-zhat, zhat-0-xhat tensors
        t(1,1)=t(1,1)+fac
        s(3,1)=s(3,1)+fac
        s(1,3)=s(1,3)+fac
        r(3,3)=r(3,3)+fac
c  mult by horizontal slowness and calc the modified T-matrix 
        do i=1,3
          do j=1,3
            t(j,i)=t(j,i)*px*px
            s(j,i)=s(j,i)*px
          end do
          t(i,i)=t(i,i)-rho(n)
        end do
c  calculate R**(-1).S, R**(-1).T, using routine solve
        nn=3
        do i=1,3
          do j=1,3
            y(j)=s(j,i)
          end do
          call solve(nn,r,x,y)
          do j=1,3
            stl(j,i)=x(j)
          end do
          nn=-3
        end do
        do i=1,3
          do j=1,3
            y(j)=t(j,i)
          end do
          call solve(nn,r,x,y)
          do j=1,3
            ttl(j,i)=x(j)
          end do
        end do
c  fill the 6x6 Q-matrix
        do i=1,3
          do j=1,3
            qq(j,i)=-stl(j,i)
            qq(j,i+3)=-ttl(j,i)
            qq(j+3,i)=0.d0
            qq(j+3,i+3)=0.d0
          end do
          qq(i+3,i)=1.d0
        end do
c  solve eigenvalue problem for polarization vectors and vertical slownesses
c  matrix system is nonsymmetric real valued
c  solution from the eispack guide
        call balanc(6,6,qq,is1,is2,fv)
        call elmhes(6,6,is1,is2,qq,iv)
        call eltran(6,6,is1,is2,qq,iv,zr)
        call hqr2(6,6,is1,is2,qq,wr,wi,zr,ierr)
        if(ierr.ne.0) then
          print *, ierr,'   error!'
          stop
        endif      
        call balbak(6,6,is1,is2,fv,6,zr)
c        print *,'for layer',n
c        print *, 'for phase velocity',cc,'  the vertical slownesses are'
c        print 101,(wr(i),wi(i),i=1,6)
c        pause
  101 format(6g12.4)
c  eigenvector unpacking, see EISPACK guide, page 88
c  bad eigenvector order is flagged by wi(i)>0. for odd i
        iflag=0
        do i=1,6
          if(wi(i).eq.0.d0) then
            if(n.eq.1) iev=0
            do j=1,6
              zi(j,i)=0.d0
            end do
          elseif(wi(i).gt.0.d0) then
c  bad eigenvector order is flagged by wi(i)>0 for even i
            if((i/2)*2.eq.i) then
              iflag=iflag+1
              iv(iflag)=i
            endif
            do j=1,6
              zi(j,i)=zr(j,i+1)
            end do
          else
            do j=1,6
              zi(j,i)=-zi(j,i-1)
              zr(j,i)=zr(j,i-1)
            end do
          endif
c  normalize by the last three indices
          sum=0.d0
          do j=4,6
            sum=sum+zr(j,i)**2+zi(j,i)**2
          end do
          sum=dsqrt(sum)
          do j=1,6
            zr(j,i)=zr(j,i)/sum
            zi(j,i)=zi(j,i)/sum
          end do         
        end do
c  assemble the stress-displacement vectors
c  calculate the traction components, with i removed
        pp(1)=dcmplx(px,0.d0)
        pp(2)=z0
        do k=1,6
          pp(3)=dcmplx(wr(k),wi(k))
          do i=1,3
            u0(i)=dcmplx(zr(i+3,k),zi(i+3,k))
          end do
          pu=z0
          pw=z0
          uw=z0
          abcde=a-b+c-2.d0*d+2.d0*e
          bce=b-4.d0*c-4.d0*e
          de=d-e
          do i=1,3
            pu=pu+pp(i)*u0(i)
            pw=pw+pp(i)*w(i,n)
            uw=uw+u0(i)*w(i,n)
          end do
          do i=1,3
            e1(i,k)=u0(i)
            e1(i+3,k)=w(i,n)*(pu*w(3,n)*bce+8.d0*pw*uw*w(3,n)*c
     x                       +2.d0*(pw*u0(3)+uw*pp(3))*e)
            e1(i+3,k)=e1(i+3,k)+pp(i)*(u0(3)*de+2.d0*uw*w(3,n)*e)
            e1(i+3,k)=e1(i+3,k)+u0(i)*(pp(3)*de+2.d0*pw*w(3,n)*e)
          end do
          e1(6,k)=e1(6,k)+pu*abcde+pw*uw*bce
c  almost lastly, mult traction by i
          do i=1,3
            e1(i+3,k)=eye*e1(i+3,k)
          end do
        end do
c  reorder into upgoing and downgoing waves
c  we use the exp(-i*omega*t) convention with z increasing downward
c  so downgoing oscillatory waves have p_z>0, k_z real
c  downgoing evanescent waves have Im(p_z)>0
c  if the axis of symmetry is tilted, there are cases where a pair of 
c  near-horizontal plane waves will be both upgoing or both downgoing
c  since Chen's algorithm depends on a 3,3 split, we must adopt a kluge
c  similarly, there are cases where the EISPACK routines dont return
c  the vertical wavenumbers in ordered pairs, but mix them up a bit
c  this seems to cause problems, so a fix is necessary
c
c  first, test for bad eigenvector order, switch k-1->k+1, k->k-1, k+1->k
c   worst case is iflag=2, real,imag1+,imag1-,imag2+,imag2-,real
        if(iflag.gt.0) then
          do i=1,iflag
            k=iv(i)
            wrr=wr(k-1)
            wii=wi(k-1)
            wr(k-1)=wr(k)
            wi(k-1)=wi(k)
            wr(k)=wr(k+1)
            wi(k)=wi(k+1)
            wr(k+1)=wrr
            wi(k+1)=wii
            do j=1,6
              pu=e1(j,k-1)
              e1(j,k-1)=e1(j,k)
              e1(j,k)=e1(j,k+1)
              e1(j,k+1)=pu
            end do
          end do
        endif
c  second, divide into upgoing and downgoing waves
        isum=0
        do k=1,6
          iv(k)=0
          if(wi(k).eq.0.d0.and.wr(k).gt.0) iv(k)=1
          if(wi(k).gt.0.d0) iv(k)=1
          isum=isum+iv(k)
        end do
c  if up and downgoing cohorts are not equal, switch the sense of the
c  pure-oscillatory wave with smallest wavenumber
  140   continue
        if(isum.ne.3) then
          wr0=0.d0
          do k=1,6
            wr0=dmax1(wr0,dabs(wr(k)))
          end do
          do k=1,6
            if(wi(k).eq.0.d0) then
              if(dabs(wr(k)).lt.wr0) then
                wr0=dabs(wr(k))
                kk=k
              endif
            endif
          end do
          if(iv(kk).eq.0) then
            iv(kk)=1
          else
            iv(kk)=0
          endif
c  check that we have equal up/down cohorts
          isum=0
          do k=1,6
            isum=isum+iv(k)
          end do
          go to 140
        endif
        jdown=1
        jup=4
c        print *,'for layer',n,'  the vert wavenums are (0=up,1=dn)'
 1001 format(i2,2g15.6)
        do k=1,6
          if(iv(k).eq.1) then
            ki=jdown
            jdown=jdown+1
          else
            ki=jup
            jup=jup+1
          endif
          do i=1,6
            ee(i,ki,n)=e1(i,k)
          end do
c  incorporate the factor of i into the stored vertical slowness 
          xnu(ki,n)=dcmplx(-wi(k),wr(k))
        end do
c  OK, here's where we check whether two downgoing stress-disp vectors
c  are nearly parallel - we check the dotproducts of displacement components
        do i=1,4
          idfct(i,n)=0
          adf((i+1)/2,n)=0.d0
        end do
        do i=1,2
          do j=i+1,3
            r1=0.d0
            r2=0.d0
            zz=z0
            do k=1,3
              r1=r1+zabs(ee(k,i,n))**2
              r2=r2+zabs(ee(k,j,n))**2
              zz=zz+ee(k,j,n)*conjg(ee(k,i,n))
            end do
            qqq=1.d0-zabs(zz)/dsqrt(r1*r2)
            if(qqq.lt.tol) then
              ccc=cc*vbar
              idfct(1,n)=i
              idfct(2,n)=j
c              print 1008,'vert slownesses',xnu(i,n),' and',xnu(j,n)
c  we average eigenvalues (vert slownesses)
c  and solve for eigenvector in subroutine defective_rf
              xnu(i,n)=(xnu(i,n)+xnu(j,n))/2.d0
              xnu(j,n)=xnu(i,n)
c  calculate the extravector for defective_rf repeated eigenvalue
              call defective_rf(i,j,n,adf(1,n),a,b,c,d,e,px)
c              print *,i,j,n,ccc,qqq,adf(1,n)
            endif
          end do
        end do
 1008 format(a,2g15.6,a,2g15.6)
c  OK, here's where we check whether two upgoing stress-disp vectors
c  are nearly parallel - we check the dotproducts of displacement components
        do i=4,5
          do j=i+1,6
            r1=0.d0
            r2=0.d0
            zz=z0
            do k=1,3
              r1=r1+zabs(ee(k,i,n))**2
              r2=r2+zabs(ee(k,j,n))**2
              zz=zz+ee(k,j,n)*conjg(ee(k,i,n))
            end do
            qqq=1.d0-zabs(zz)/dsqrt(r1*r2)
            if(qqq.lt.tol) then
              ccc=cc*vbar
              idfct(3,n)=i
              idfct(4,n)=j
c              print 1008,'vert slownesses',xnu(i,n),' and',xnu(j,n)
c  we average the eigenvalues
              xnu(i,n)=(xnu(i,n)+xnu(j,n))/2.d0
              xnu(j,n)=xnu(i,n)
c  calculate the extravector for defective_rf repeated eigenvalue
c  as well as coefficient (adf) of linear (z-z0) term
              call defective_rf(i,j,n,adf(2,n),a,b,c,d,e,px)
c              print *,i,j,n,ccc,qqq,adf(2,n)
            endif
          end do
        end do
      end do
c  now, must identify which upgoing waves in the halfspace are P,SV,SH
c  crud, this goes back to array ee
c  3: SH is y-motion
c  2: SV is (-sqrt((1/vs)**2-p_x**2),0,-p_x)  ! recall that z points down
c  1: P is (p_x,0,-sqrt((1/vp)**2-p_x**2)
c  so we branch on size of u_y, and relative sign of u_x and u_z
c      print *,'in the halfspace:'
      do i=4,6
c        print *,'for i*k_z=',xnu(i,nlp),', the disp-stress vector is'
        do j=1,6
          xi(j)=dimag(ee(j,i,nlp))
          xr(j)=dreal(ee(j,i,nlp))
        end do
c        print 101,(xr(j),j=1,6),(xi(j),j=1,6)
      end do
      do i=4,6
        ips(i-3)=3
        if(zabs(ee(2,i,nlp)).lt.dsqrt(tol)) then  ! not SH
          test=dreal(ee(1,i,nlp))/dreal(ee(3,i,nlp))
          if(test.gt.0.d0) then
            ips(i-3)=2
          else
            ips(i-3)=1
          endif
        endif
      end do
c      print *,'wave types:',(ips(i),i=1,3)        
      return
      end
      subroutine respget(nl,om,cc,resp)
c  returns surface response for a stack of anisotropic layers
c  incident p,sv,sh waves with freq om and phase velocity cc
c  iev=1 if the waves are evanescent in the top layer, iev=0 otherwise
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      complex*16 pp,u0,ee,z1,z0,xnu,eye,e1,e2,zla,rtm
      complex*16 rt,tt,rt0,trc,xl,pfac,u,co,resp,ur
      common/stfff/w(3,101),t(3,3),ttl(3,3),s(3,3),stl(3,3),
     x                     			 r(3,3),x(3),y(3)
      common/model/z(100),dz(100),rho(101),vp(101),vs(101),vp2(101),
     x               vp4(101),vs2(101),vss(101)
      common/model2/xmu(101),xla(101),xmu2(101),xla2(101),xla4(101)
      common/propag/xnu(6,101),xl(6,100),pfac(6,3),u(3,6)
      common/defect/idfct(4,101),adf(2,101)
      common/rrt/rtm(6,6,100)
      common/mstff/qq(6,6),wr(6),wi(6),zr(6,6),zi(6,6),iv(6),fv(6)
      common/pstff/pp(3),u0(3),ee(6,6,101),e1(6,6),e2(6,6),zla(6)
      common/pstf2/co(6,101),ur(3)
      common/qstff/qi(6,6),xr(6),xi(6),yr(6),yi(6),ips(3)
      common/rstff/rt(3,3,101),tt(3,3,101),rt0(3,3),trc(3,3)
      dimension resp(3,3)
      data pi/3.14159265358979d0/,eps/1.d-6/,tol/1.d-7/
c  set iev=1
c  toggle to iev=0 if there is a purely propagating wave in the top layer n=1
      iev=1
      z1=dcmplx(1.d0,0.d0)
      z0=dcmplx(0.d0,0.d0)
      eye=dcmplx(0.d0,1.d0)
      rbar=5.515d3
      ren=1.075190645d-3
      radi=6.371d6
      vbar=radi*ren
      con=rbar*radi*radi*ren*ren
      nlp=nl+1
      nlm=nl-1
c  first calculate vertical wavenumbers and propagating waves for each layer
c    an eigenvector problem was solved in prior subroutine
c   results stored in array ee(6,6,101)
c  in general, the evanescent vertical wavenumbers have nonzero real parts 
c   complex exponential fct is used to avoid endless branching
c
c  calculate modified R/T coefficients
c  first calc the propagation factors
c  note that for dipping fast axes the upgoing and downgoing wavenumbers are
c  independent, so we must calc all to be safe
      do n=1,nl
        do k=1,3
          xl(k,n)=zexp(om*xnu(k,n)*dz(n))		! downgoing
          xl(k+3,n)=zexp(-om*xnu(k+3,n)*dz(n))	! upgoing
        end do
      end do
c      do i=1,6
c        print 1002,xnu(i,3),xl(i,3)
c      end do
 1002 format('i*k_z:',2g15.6,',  propfac is',2g15.6)
c  calculate modified R/T coefficients at each interface
      do n=1,nl
c  rearrange to e1: waves approaching and e2: waves leaving an interface
        do k=1,3
          do i=1,6
            e1(i,k)=ee(i,k,n+1)
            e2(i,k)=ee(i,k,n)
            e1(i,k+3)=-ee(i,k+3,n)
            e2(i,k+3)=-ee(i,k+3,n+1)
          end do
          zla(k)=xl(k,n)
          if(n.lt.nl) then
            zla(k+3)=xl(k+3,n+1)   
          else 
c  reference the upcoming wave amplitude to the top of halfspace
c  therefore no propagation factor, not valid for evanescent waves in halfspace
c  in surface wave code this is zero, so that upgoing evanescent waves vanish
            zla(k+3)=1.d0
          endif
        end do
c mult the columns of e2
        do k=1,6
          do i=1,6
            e2(i,k)=e2(i,k)*zla(k)
          end do
        end do
c  the possibility of defective_rf matrices must be contemplated here
c  k=1,2,3 columns are downgoing in nth layer
c  k=4,5,6 columns are upgoing in (n+1)th layer
c  the vector e2(.,k1) has already been multiplied by exponential factor zla
        if(idfct(1,n).ne.0) then
          k1=idfct(1,n)
          k2=idfct(2,n)
          do i=1,6
            e2(i,k2)=e2(i,k2)+adf(1,n)*dz(n)*e2(i,k1)
          end do
        endif
c  the sign change on dz is for upgoing waves
        if(idfct(3,n+1).ne.0) then
          k1=idfct(3,n+1)
          k2=idfct(4,n+1)
          do i=1,6
            e2(i,k2)=e2(i,k2)-adf(2,n+1)*dz(n+1)*e2(i,k1)
          end do
        endif
c  in order to use csolve to invert e1, must separate into real/imag parts
c  its clumsy, but im lazy
c  we calc e1^{-1}\cdot e2\cdot \Gamma one column at a time
        do k=1,6
          do i=1,6
            qq(i,k)=dreal(e1(i,k))
            qi(i,k)=dimag(e1(i,k))
          end do
        end do
        nn=6
        do k=1,6
          do i=1,6
            yr(i)=dreal(e2(i,k))
            yi(i)=dimag(e2(i,k))
          end do
          call csolve(nn,qq,qi,xr,xi,yr,yi)
          nn=-6
          do i=1,6
            rtm(i,k,n)=dcmplx(xr(i),xi(i))
          end do
        end do  
      end do
c  calc R_ud at the free surface
c  note that first two factors in Chen (20) dont collapse
c  mult by inv-matrix one column at a time
      do k=1,3
        do i=1,3
          rt0(i,k)=ee(i+3,k+3,1)*xl(k+3,1)
          s(i,k)=dreal(ee(i+3,k,1))
          t(i,k)=dimag(ee(i+3,k,1))
        end do
      end do
c  the possibility of defective_rf matrices must be contemplated here
c  these waves are upgoing in 1st layer
c  the sign change on dz is for upgoing waves, and xl(k1,1)=xl(k2,1)
      if(idfct(3,1).ne.0) then
        k1=idfct(3,1)
        k2=idfct(4,1)-3
        do i=1,3
          rt0(i,k2)=rt0(i,k2)-adf(2,1)*dz(1)*ee(i+3,k1,1)*xl(k1,1)
        end do
      endif
      nn=3
      do k=1,3
        do i=1,3
          yr(i)=dreal(rt0(i,k))
          yi(i)=dimag(rt0(i,k))
        end do
        call csolve(nn,s,t,xr,xi,yr,yi)
        nn=-3
        do i=1,3
          rt0(i,k)=-dcmplx(xr(i),xi(i))
        end do
      end do  
c  recursive calc of generalized R/T coefs:
c  in contrast to the surface-wave code, we start from the top layer and
c  iterate down to the halfspace
c  we also uses submatrices of generalized R/T matrix in different order
      do n=1,nl
c  first the generalized upward-transmission coef:
        do k=1,3
          do i=1,3
            trc(i,k)=z0
            if(n.gt.1) then
              do j=1,3
                trc(i,k)=trc(i,k)-rtm(i+3,j,n)*rt(j,k,n-1)
              end do
            else
c  use free-surface reflection matrix in top layer (interface "zero")
              do j=1,3
                trc(i,k)=trc(i,k)-rtm(i+3,j,n)*rt0(j,k)
              end do
            endif
          end do
          trc(k,k)=trc(k,k)+z1
        end do  
        do k=1,3
          do i=1,3
            s(i,k)=dreal(trc(i,k))
            t(i,k)=dimag(trc(i,k))
          end do
        end do
        nn=3
        do k=1,3
          do i=1,3
            yr(i)=dreal(rtm(i+3,k+3,n))
            yi(i)=dimag(rtm(i+3,k+3,n))
          end do
          call csolve(nn,s,t,xr,xi,yr,yi)
          nn=-3
          do i=1,3
            tt(i,k,n)=dcmplx(xr(i),xi(i))
          end do
        end do  
c  next the generalized reflection coef:
        do k=1,3
          do i=1,3
            trc(i,k)=z0
            if(n.gt.1) then
              do j=1,3
                trc(i,k)=trc(i,k)+rt(i,j,n-1)*tt(j,k,n)
              end do
            else
c  use free-surface reflection matrix in top layer (interface "zero")
              do j=1,3
                trc(i,k)=trc(i,k)+rt0(i,j)*tt(j,k,n)
              end do
            endif
          end do
        end do  
        do k=1,3
          do i=1,3
            rt(i,k,n)=rtm(i,k+3,n)
            do j=1,3
              rt(i,k,n)=rt(i,k,n)+rtm(i,j,n)*trc(j,k)
            end do
          end do
        end do  
      end do
 1001 format(6f14.6)
c      print *,'free-surface reflection'
c      print 1001,((rt0(i,j),j=1,3),i=1,3)
c      do n=1,nl
c        print *,'interface',n
c        print 1001,((rt(i,j,n),j=1,3),i=1,3)
c        print 1001,((tt(i,j,n),j=1,3),i=1,3)
c      end do
c  using the p,sv,sh identification, we propagate upward to the surface,
c  calculate the wave coefs in the top layer, then the particle displacement
      do iup=1,3
        do i=1,3
          co(i+3,nlp)=z0
        end do
        co(iup+3,nlp)=z1
c  from upgoing coefs in the n+1 layer, calculate
c  upgoing coefs in the nth layer, downgoing coefs in the n+1 layer
        do n=nl,1,-1
          do i=1,3
            co(i+3,n)=z0
            co(i,n+1)=z0
            do j=1,3
              co(i+3,n)=co(i+3,n)+tt(i,j,n)*co(j+3,n+1)
              co(i,n+1)=co(i,n+1)+rt(i,j,n)*co(j+3,n+1)
            end do
          end do
        end do
c  then downgoing coefs in the top layer:
        do i=1,3
          co(i,1)=z0
          do j=1,3
            co(i,1)=co(i,1)+rt0(i,j)*co(j+3,1)
          end do
        end do
c        print *,'upgoing coefs'
c        print 1001,((co(j+3,n),j=1,3),n=1,nlp)
c        print *,'downgoing coefs'
c        print 1001,((co(j,n),j=1,3),n=1,nlp)
c  calc the surface displacement
        h1=0.d0
        h2=z(1)
        do i=1,3
          ur(i)=z0
          do k=1,3
            ur(i)=ur(i)+co(k,1)*ee(i,k,1)
     x          +co(k+3,1)*ee(i,k+3,1)*(zexp(om*xnu(k+3,1)*(-h2)))
          end do
c  check for the xtra terms associated with defective_rf matrices
          if(idfct(1,1).ne.0) then
            ii=idfct(1,1)
            jj=idfct(2,1)
            ur(i)=ur(i)+co(jj,1)*adf(1,1)*ee(i,ii,1)*(-h1)
          endif
          if(idfct(3,1).ne.0) then
            ii=idfct(3,1)
            jj=idfct(4,1)
            ur(i)=ur(i)
     x +co(jj,1)*adf(2,1)*ee(i,ii,1)*(-h2)*(zexp(om*xnu(ii,1)*(-h2)))
          endif
        end do
c  copy the surface displacement into the response matrix
        do i=1,3
          resp(i,ips(iup))=ur(i)
        end do            
      end do     
      return
      end
      subroutine defective_rf(i,j,n,adf,a,b,c,d,e,px)
c  kluge for dealing with nearly defective_rf propagator matrices
c  in which the eigenvectors, 
c  which represent the particle motion of upgoing and downgoing waves
c  become nearly parallel.
c  in this case the solution for system of ODEs is
c  a_1 \bf_1 e^xnu*(z-z0) + a_2*(\bf_2 + adf*(z-z0)*\bf_1)e^xnu*(z-z0)
c
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      complex*16 pp,u0,ee,z1,z0,znu,xnu,e1,e2,zla,xl,u,pfac,eye
      complex*16 zq1,zq2,u1,u2,zq3,xee
      common/stfff/w(3,101),t(3,3),ttl(3,3),s(3,3),stl(3,3),
     x                     			 r(3,3),x(3),y(3)
      common/propag/xnu(6,101),xl(6,100),pfac(6,3),u(3,6)
      common/defect1/zq1(3,3),zq2(3,3),u1(3),u2(3),zq3(2,2),xee(3)
      common/defect2/edr(6),edi(6),qdr(5,5),qdi(5,5),ydr(6),ydi(6)
      common/defect3/q1r(3,3),q1i(3,3),q2r(3,3),q2i(3,3),fv2(3),fv3(3)
      common/mstff/qq(6,6),wr(6),wi(6),zr(6,6),zi(6,6),iv(6),fv(6)
      common/pstff/pp(3),u0(3),ee(6,6,101),e1(6,6),e2(6,6),zla(6)
      common/qstff/qi(6,6),xr(6),xi(6),yr(6),yi(6),ips(3)
      z1=dcmplx(1.d0,0.d0)
      z0=dcmplx(0.d0,0.d0)
      eye=dcmplx(0.d0,1.d0)
c  for the extravector, need to solve system of equations
c  based on original 6x6 Q matrix
c  the plane-wave solutions generalize to the form
c  u0*e^{i*nu*(z-z0)}  and  u1*e^{i*nu*(z-z0)} + adf* u0*(z-z0)*e^{i*nu*(z-z0)}
c  u1 is the solution to 
c  (\bTtil + nu*\bStil + nu^2*\bI).u1=i*adf*(\bStil + 2*nu*\bI).u0
c  in practice, we absorb the adf factor into u1, then normalize
c  (\bTtil + nu*\bStil + nu^2*\bI).(u1/adf)=i*(\bStil + 2*nu*\bI).u0
c  since nu is the known eigenvalue of u0, the solution is easier
c  form the matrices on either side
      znu=-eye*xnu(i,n)
      do ii=1,3
        do jj=1,3
          zq1(jj,ii)=dcmplx(ttl(jj,ii),0.d0)+znu*stl(jj,ii)
          zq2(jj,ii)=dcmplx(stl(jj,ii),0.d0)
        end do
        zq1(ii,ii)=zq1(ii,ii)+znu*znu
        zq2(ii,ii)=zq2(ii,ii)+2.d0*znu
      end do
c  we wish to find the eigenvector of the near-defective_rf matrix 
c  in the region where its eigenvectors are numerically unstable
c   we explicitly calculate the eigenvector with smallest right-eigenvalue of
c  (\bTtil + nu*\bStil + nu^2*\bI)=zq1
c  copy into real, imag matrices
      do ii=1,3
        do jj=1,3
          q1r(jj,ii)=dreal(zq1(jj,ii))
          q1i(jj,ii)=dimag(zq1(jj,ii))
        end do
      end do      
c  into eispack
      call cbal(3,3,q1r,q1i,low,igh,fv)
      call corth(3,3,low,igh,q1r,q1i,fv2,fv3)
      call comqr2(3,3,low,igh,fv2,fv3,q1r,q1i,wr,wi,q2r,q2i,ierr)
      if(ierr.ne.0) go to 400
      call cbabk2(3,3,low,igh,fv,3,q2r,q2i)
      amn=wr(1)**2+wi(1)**2
      ij=1
      do ii=2,3
        amm=wr(ii)**2+wi(ii)**2
        if(amm.lt.amn) then
          ij=ii
          amn=amm
        endif
      end do
      sum=0.d0
      do ii=1,3
        u0(ii)=dcmplx(q2r(ii,ij),q2i(ii,ij))
        sum=sum+zabs(u0(ii))**2
      end do
      sum=dsqrt(sum)
      do ii=1,3
        u0(ii)=u0(ii)/sum
      end do
c  assemble the ith stress-displacement vector
c  calculate the traction components, with i removed
      pp(1)=dcmplx(px,0.d0)
      pp(2)=z0
      pp(3)=znu
      pu=z0
      pw=z0
      uw=z0
      abcde=a-b+c-2.d0*d+2.d0*e
      bce=b-4.d0*c-4.d0*e
      de=d-e
      do ii=1,3
        pu=pu+pp(ii)*u0(ii)
        pw=pw+pp(ii)*w(ii,n)
        uw=uw+u0(ii)*w(ii,n)
      end do
      do ii=1,3
        ee(ii,i,n)= u0(ii)
        ee(ii+3,i,n)=w(ii,n)*(pu*w(3,n)*bce+8.d0*pw*uw*w(3,n)*c
     x                                  +2.d0*(pw*u0(3)+uw*pp(3))*e)
        ee(ii+3,i,n)=ee(ii+3,i,n)+pp(ii)*(u0(3)*de+2.d0*uw*w(3,n)*e)
        ee(ii+3,i,n)=ee(ii+3,i,n)+u0(ii)*(pp(3)*de+2.d0*pw*w(3,n)*e)
      end do
      ee(6,i,n)=ee(6,i,n)+pu*abcde+pw*uw*bce
c  almost lastly, mult traction by i
      do ii=1,3
        ee(ii+3,i,n)=eye*ee(ii+3,i,n)
      end do
c  extract u0 from ee(*,i,n) use it to calculate the additional traction terms
c  and store in ee(*,j,n)
c  additional traction terms involve gradient of (z-z0)
c  so can be calculated from standard formulas with \bk=zhat
c  we dont multiply by i
      pp(1)=z0
      pp(2)=z0
      pp(3)=z1
      pu=z0
      pw=z0
      uw=z0
      abcde=a-b+c-2.d0*d+2.d0*e
      bce=b-4.d0*c-4.d0*e
      de=d-e
      do ii=1,3
        u0(ii)=ee(ii,i,n)
        pu=pu+pp(ii)*u0(ii)
        pw=pw+pp(ii)*w(ii,n)
        uw=uw+u0(ii)*w(ii,n)
      end do
      do ii=1,3
        xee(ii)=w(ii,n)*(pu*w(3,n)*bce+8.d0*pw*uw*w(3,n)*c
     x            		        +2.d0*(pw*u0(3)+uw*pp(3))*e)
        xee(ii)=xee(ii)+pp(ii)*(u0(3)*de+2.d0*uw*w(3,n)*e)
        xee(ii)=xee(ii)+u0(ii)*(pp(3)*de+2.d0*pw*w(3,n)*e)
      end do
      xee(3)=xee(3)+pu*abcde+pw*uw*bce
c  extract u0 from ee(*,i,n), mult by i*(\bStil + 2*nu*\bI), replace in u0
      do ii=1,3
        u0(ii)=z0
        do jj=1,3
          u0(ii)=u0(ii)+zq2(ii,jj)*ee(jj,i,n)
        end do
        u0(ii)=eye*u0(ii)
      end do
 1002 format(3(2g14.6,3x))
c  for znu NOT an eigenvalue, 
c  but rather the average of closely-space eigenvalues
c  in this case, zq1 is nonsingular, and we just solve for u1 
      do ii=1,3
        yr(ii)=dreal(u0(ii))
        yi(ii)=dimag(u0(ii))
        do jj=1,3
          q1r(jj,ii)=dreal(zq1(jj,ii))
          q1i(jj,ii)=dimag(zq1(jj,ii))
        end do
      end do
      call csolve(3,q1r,q1i,xr,xi,yr,yi)
      do ii=1,3
        u1(ii)=dcmplx(xr(ii),xi(ii))
      end do
c   End, different tactic
c
c  normalize  
      sum=0.d0
      do ii=1,3
        sum=sum+zabs(u1(ii))**2
      end do
      sum=dsqrt(sum)
      do ii=1,3
        u1(ii)=u1(ii)/sum
      end do
c  adf is the normalization constant
      adf=1.d0/sum
c  calculate the traction
c  and place the new stress-displacement vector in column j
c  pp is the wavenumber vector, and first two components are already in place
      pp(1)=dcmplx(px,0.d0)
      pp(2)=z0
      pp(3)=znu
      pu=z0
      pw=z0
      uw=z0
      abcde=a-b+c-2.d0*d+2.d0*e
      bce=b-4.d0*c-4.d0*e
      de=d-e
      do ii=1,3
        pu=pu+pp(ii)*u1(ii)
        pw=pw+pp(ii)*w(ii,n)
        uw=uw+u1(ii)*w(ii,n)
      end do
      do ii=1,3
        ee(ii,j,n)=u1(ii)
        ee(ii+3,j,n)=w(ii,n)*(pu*w(3,n)*bce+8.d0*pw*uw*w(3,n)*c
     x                                  +2.d0*(pw*u1(3)+uw*pp(3))*e)
        ee(ii+3,j,n)=ee(ii+3,j,n)+pp(ii)*(u1(3)*de+2.d0*uw*w(3,n)*e)
        ee(ii+3,j,n)=ee(ii+3,j,n)+u1(ii)*(pp(3)*de+2.d0*pw*w(3,n)*e)
      end do
      ee(6,j,n)=ee(6,j,n)+pu*abcde+pw*uw*bce
c  almost lastly, mult traction by i 
c  and add extra traction from (z-z0) term (not mult by i)
c  TEST - mult xee by zero, see if it is important --- it IS important
      do ii=1,3
        ee(ii+3,j,n)=eye*ee(ii+3,j,n)+adf*xee(ii)
      end do
      return
 400  print *,'eispack error'
      stop
      end

                                                                 
c  program to test solve
c      implicit real*8 (a-h,o-z)
c      dimension a(100),x(10),y(10),aa(10,10)
c      call tnoua('n= ')
c      read(5,*) n
c      do 100 i=1,n
c      do 100 j=1,n
c      print 101, i,j
c  101 format('$enter a(',2i2,') ')
c      read(5,*) a(i+n*(j-1))
c      aa(i,j)=a(i+n*(j-1))
c  100 continue
c      do 200 i=1,n
c      print 102,i
c  102 format('$enter y(',i2,') ')
c      read(5,*) y(i)
c  200 continue
c      call solve(n,a,x,y)
c      print 103, (i,x(i),i=1,n)
c  103 format(' x(',i2,')=',e15.4)
cc  check by substitution
c      do 300 i=1,n
c      y(i)=0.d0
c      do 300 j=1,n
c  300 y(i)=y(i)+aa(i,j)*x(j)
c      call tnou('substitute')
c      print 103,(i,y(i),i=1,n)
c      stop
c      end

