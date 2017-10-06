c  aniprop - program to calculate propagating modes of anisotropic layer
c  writes a file of dispersion curves at evenly-space freq points
c
c reads a layered model like the following file (ignore the "c "s)
c
c  K&H SoCal model, deep crust horizontal anisotropy	TITLE
c  3							# OF LAYERS OVER HSPACE
c  45 45					THETA,PHI ORIENTATION ANGLES
c  4000 5500 0.06 0.00 3175 0.03 2600		FOR SYMMETRY AXIS
c  0 0					
c  27400 6300 0.00 0.00 3637 0.00 2800		DEPTH (M), VP (M/S), "B", "C",
c  90 45 					VS (M/S), "E" 
c  32400 6800 0.04 0.00 3925 0.02 2900		B,C,E ARE ANISOTROPIC PARAMETERS
c  0 0
c  60000 7800 0.00 0.00 4500 0.00 3200		NOTE: HSPACE MUST BE ISOTROPIC
c
c    revised to calculate dispersion curves: 7/7/95
c  revised to make rootfinder more focussed: 7/11/95
c  revised to make computations dimensionless: started 8/3/95, finished 8/14
c  bugfix for tilted axis of symmetry: 8/16/95
c  bugfix for group velocity: 8/22/95
c  bugfix for defective matrices: 8/22/95
c  bugfix for transition to evanescence in top layer: 9/10/95
c  multiple bugfixes in the rootfinder and in the test for 
c     horizontal waves in the surface layer  2/5/00
c  bugfix for underflow in routine grvel 3/8/00
c  
c
c  compile sequence 
c  eislib is the Eispack library
c  f77 -o aniprop aniprop.f /Users/jjpark/Ritz/eislib.a
c
c   compile sequence with Park-specific plotting programs for debugging
c  plotlib&jlib link with a complaint-that-I-ignore in Solaris 4.x
c  xf77 -o /Users/jjpark/bin/aniprop aniprop.f /Users/jjpark/Plotxy/plotlib.a /Users/jjpark/Ritz/eislib.a /Users/jjpark/Ritz/jlib.a
c
c  for hexagonally symmetric media
c  reads fast axis orientation, constants A,B,C,D,E from file animodel
c  calculate quadratic eigenvalue problem based on the Christoffel matrix
c  see appendix of P. Shearer's thesis
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      real*4 cc4,zz4,x4,frqq,uuu,x1,x2,dum4,u4,ddtt,time,ahead,evaleq
      character*80 name,title
      character*4 chead(158),ichan(3)
      complex*16 pp,u0,ee,z1,z0,zz,xnu,zzz,e1,e2,zla,xl,c,eigf
      complex*16  rt,tt,rt0,trc,pfac,u,eye,uu,ur,usz,usx,zz0
      character*5 outfile(3)
      common/stfff/w(3,101),t(3,3),ttl(3,3),s(3,3),stl(3,3),
     x                                           r(3,3),x(3),y(3)
c  max is 100 layers over hspace (would *you* parameterize 100 layers??)
      common/model/z(100),dz(100),rho(101),vp(101),vs(101),vp2(101),
     x               vp4(101),vs2(101),vss(101)
      common/model2/xmu(101),xla(101),xmu2(101),xla2(101),xla4(101)
      common/propag/xnu(6,101),xl(6,100),pfac(6,3),u(3,6)
      common/mstff/qq(6,6),wr(6),wi(6),zr(6,6),zi(6,6),iv(6),fv(6)
      common/pstff/pp(3),u0(3),ee(6,6,101),e1(6,6),e2(6,6),zla(6)
      common/rstff/rt(3,3,101),tt(3,3,101),rt0(3,3),trc(3,3)
      common/nstff/cc4(10000),zz4(10000),zzz(10000),ccc(10000),x4(10000)
      common/pstf2/c(6,101),eigf(1000,3)
      common/disper/cvel(4096,1616),gvel(4096,1616),nc(1616),frqq(8192)
      common/disper2/roota(1616),rootb(1616),jtrval(1616),kroots(4096)
      common/evanes/ievan(10000)
      common/vertwav/vertr(6),verti(6)
      common/source/ur(3),usx(3),usz(3),amom(3,3),uu(3)
      common/synthe/uuu(16384,3),dum4(16384,3),u4(16384)
      common/csynth/uur(16384),uui(16384)
      common/header/ahead(158)
      dimension iah(158)
      data pi/3.14159265358979d0/,eps/1.d-6/,tol/1.d-3/
      equivalence (iah,ahead),(chead,ahead)
      data outfile/'SYN.R','SYN.T','SYN.Z'/,ichan/'BHR ','BHT ','BHZ '/
c  we reduce the condition numbers of matrices by 
c  normalizing physical quantities to make them dimensionless
c  Ill use the normal-mode normalizations
c  which are a little peculiar for the crust, but what the hey!
      rbar=5.515d3
      radian=180.d0/pi
      ren=1.075190645d-3
      radi=6.371d6
      vbar=ren*radi
      con=rbar*vbar**2
c  ebar is energy normalization, for moment tensor
      ebar=con*radi**3
      z1=dcmplx(1.d0,0.d0)
      z0=dcmplx(0.d0,0.d0)
      eye=dcmplx(0.d0,1.d0)
c initialize cvel array to zero
      maxbr=1616
      do j=1,maxbr
        do i=1,4096
          cvel(i,j)=0.d0
          gvel(i,j)=0.d0
        end do
      end do  
      z1=dcmplx(1.d0,0.d0)
      z0=dcmplx(0.d0,0.d0)
      tolint=100.d0/vbar
      nfrq=500
      frqmax=0.5d0
c      nfrq=2000
c      frqmax=2.0d0
c      nfrq=200
c      frqmax=0.2d0
      ntry=1000
      mtry=5
      nfstart=4
 102  format(a)
      open(7,file='animodel',form='formatted')
      read(7,102) title
      print *,title
      read(7,*) nl
c  read in theta,phi in degrees - polar coords of fast axis
c  the coordinates are Z-down, X-propagation direction, Y toward the observer
c  seen from above, therefore, the azimuthal angle for W-hat is rotated CW from X toward Y
      nlp=nl+1
      nlm=nl-1
      do i=1,nlp
        read(7,*) theta,phi
        w(1,i)=dsin(theta/radian)*dcos(phi/radian)
        w(2,i)=dsin(theta/radian)*dsin(phi/radian)
        w(3,i)=dcos(theta/radian)
        print *,(w(j,i),j=1,3)
c  read depth to ith interface, vp (m/sec), pk-to-pk cos(2th) relative P pert
c  pk-to-pk cos(4th) relative P pert, v_s,  pk-to-pk cos(2th) relative S pert
c  density (kg/m**3)
        read(7,*) z(i),vp(i),vp2(i),vp4(i),vs(i),vs2(i),rho(i)
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
      print *, 'organ-pipe mode count at 1 Hz in stack of layers: S & P'
      print 104,sn,pn
  104 format(2f10.1)
c  search for cmin, cmax is halfspace velocity
      cmin=vs(1)
      vss(1)=vs(1)
      do i=2,nlp
        if(cmin.gt.vs(i)) cmin=vs(i)
        vss(i)=vs(i)
      end do
      cmax=vs(nlp)
      print *,'ntry',ntry
c  order interval velocities
      do i=1,nl
        do j=i+1,nlp
          if(vss(i).gt.vss(j)) then
            ssv=vss(i)
            vss(i)=vss(j)
            vss(j)=ssv
          endif
        end do
      end do
      cmin=vss(1)
      print *,cmin,cmax,' =cmin,cmax'
      print *,(vss(j),j=1,nlp+1)
c  keep a running table of the first 10 dispersion branvches
      open(7,file='out_cvel',form='formatted')
      open(8,file='out_gvel',form='formatted')
c  we search phase velocity via a focussed search.
c  first for smallest few frequencies, we search all [cmin,cmax]
c  and find the fundamental Rayleigh, maybe the fundamental Love wave
c  in succeeding frequencies, we use the solutions from the previous iteration 
c  to define intervals to search over.  The intervals are defined by 
c  group-velocity bounds, ranging from 0.9*cmin 
c  (for fundamental Rayleigh in the top layer) to cmax
c  some calculus shows that dc=(c/f)(1-c/U)df, with U=gvel<c
c  so dc+ = 0
c  and  dc- = (c/f)(c/cmin*(1./0.9)-1)
c  if the search intervals of two modes overlap, we coalesce the intervals 
c  which will densify the sampling
c  mtry grid points in each interval.  
c  Double-overlap has 2*mtry, triple has 3*mtry, und so weiter
      jj=0
      ddc=(cmax-cmin-2*eps)/(ntry-1)
      do sl=cmin+eps,cmax-eps,ddc
        jj=jj+1
        ccc(jj)=sl
        cc4(jj)=ccc(jj)*vbar
      end do
      ntry=jj
c  initialize overtone-branch counts
      do i=1,maxbr
        nc(i)=0
      end do
c  we step forward in frequency
      dfrq=frqmax/nfrq
c  we count backwards in frqq & cvel, but forwards in calculation
      do i=1,nfrq
        frqq(nfrq+1-i)=i*dfrq
      end do
      f=0.d0
      do ifrq=1,nfstart
        nm=0
        f=f+dfrq
ccccccccccccccccc test
c        f=0.25d0
ccccccccccccccccc
        om=2.d0*pi*f/ren
c        print *,ifrq,f,om 
c  first make an array of values to look for zeroes in
c  note that we arent trying to capture the fundamental Rayleigh wave
c   which usually has c<vsmin
        do jj=1,ntry
c          if(jj.gt.0.95*ntry) print *,jj
          cc=ccc(jj)
          call zzget(nl,cc,om,zz,iev)
          zzz(jj)=zz
c          if(jj.gt.0.95*ntry) print *,cmin,cc,cmax,zz
        end do
c        do i=1,ntry
c          zz4(i)=zabs(zzz(i))
c        end do
c        call plotit(cc4,zz4,dum,ntry,'Propagating Mode Condition',
c     x  'Phase velocity (km/sec)','Surface Traction',2,0,.05,0,21)
c        do i=1,ntry
c          x4(i)=i
c          zz4(i)=imag(zzz(i))
c        end do
c        call plotit(x4,zz4,dum,ntry,'rescaled',
c     x  'index count','Traction condition',2,0,0.05,0,0)
c        do i=1,ntry
c          zz4(i)=real(zzz(i))
c        end do
c        call plotit(x4,zz4,dum,ntry,'rescaled',
c     x  'index count','Traction condition',2,0,0,0,22)
c  zero in on normal modes
        do ij=2,ntry
c          print *,ij
          a1=dimag(zzz(ij-1))
          a2=dimag(zzz(ij))
          ar1=dreal(zzz(ij-1))
          ar2=dreal(zzz(ij))
c  we test on imag -- assumes that lowest-freq Rayleigh is not purely evanescent
          if(a1*a2.lt.0.d0.and.dabs(ar2).lt.1.d0) then
            ick=1
c  zoom in
            cc1=ccc(ij-1)
            cc2=ccc(ij)
            iter=0
            call zoom_in(nl,a1,a2,cc1,cc2,om,iter,cc,zz,ick)
            cccc=cc*vbar
c            print *,iter,cccc,zz
c  save if no flag - we count backwards in cvel
c  last check: horizontal waves in the top layer are spurious, so kill them
c  the test is based on the propagator in the top layer
c  a horizontal wave in the top layer will have complex propagtor = 1.d0
c  the arrays vertr and verti are filled in the subroutine zzget
            if(iter.ge.0) then
              iflag=0
              do kji=1,6
                test=dabs(vertr(kji)-1.d0)+dabs(verti(kji))
                if(test.lt.1.d-6) iflag=1
              end do
c              if(dabs(cc/cmin-1.d0).gt.1.d-10) then
              if(iflag.eq.0) then
c  get the group velocity via variational integrals
                call grvel(nl,cc,om,gv)
                nm=nm+1
c  these arrays are stored in units of m/sec
                cvel(nfrq+1-ifrq,nm)=cc*vbar
                gvel(nfrq+1-ifrq,nm)=gv*vbar
                nc(nm)=nc(nm)+1
              endif
            endif
          endif
        end do
        kroots(ifrq)=nm
      end do
  103 format(2g15.5)
      nroot=nm
      print *,nroot,' = nroot'
c  next we search for roots at higher frequencies
      do ifrq=nfstart+1,nfrq
c  set up the cvel search intervals -- uses f from previous iteration
c  we do not need to normalize cycle freq, since it appears as ratio dfrq/f
c
c  here we try to outsmart the algorithm, which may have lost roots
c  if nroot has decreased from a "long-term" average, we extrapolate from
c  the frequency before the "loss," looking to pick up the errant root
c  
        if((kroots(ifrq-3)+kroots(ifrq-2))/2.gt.kroots(ifrq-1)) then
          ddfrq=dfrq*2.d0/f
          iifrq=3
          nroot=kroots(ifrq-2)
      elseif((kroots(ifrq-4)+kroots(ifrq-3))/2.gt.kroots(ifrq-1)) then
          ddfrq=dfrq*3.d0/f
          iifrq=4
          nroot=kroots(ifrq-3)
        else
          ddfrq=dfrq/f
          iifrq=2
          nroot=kroots(ifrq-1)
        endif
        gfac=1.1d0+0.2d0*dexp(-(f/0.1)*dlog(2.d0))
        do kroot=1,nroot
          root=cvel(nfrq-ifrq+iifrq,kroot)/vbar
          gv=gvel(nfrq-ifrq+iifrq,kroot)/vbar  
          roota(kroot)=root*(1.d0+ddfrq*(1.d0-root/(0.8*gv)))
          rootb(kroot)=root*(1.d0+ddfrq*(1.d0-root/(1.2*gv)))
          jtrval(kroot)=mtry
c  at very low f, roota can undershoot cmin -- problem addressed by tolint
        end do
c  we put another interval near cmax to pick up a new overtone branch
        ntrval=nroot+1
        rootb(ntrval)=cmax-eps
c        roota(ntrval)=cmax-tolint
c  shrink the upper interval from 16% to 3% of cvel interval as f>>0.1 Hz
        propor=0.03+0.13*dexp(-(f/0.1)**2)
        roota(ntrval)=(1.d0-propor)*cmax+propor*cmin
        jtrval(ntrval)=mtry
c  coalesce search intervals
        itrval=1
c        print *,ntrval
          print *,cmin,roota(ntrval),rootb(ntrval),cmax
c        do i=1,ntrval
c          print *,roota(i),rootb(i)
c        end do
        do kroot=1,nroot
          if(rootb(itrval).gt.roota(kroot+1)) then
            jtrval(itrval)=jtrval(itrval)+mtry
c  need to take max and min because sometimes branches cross at finite angles
c  and the roota intervals pass through each other
            roota(itrval)=dmin1(roota(itrval),roota(kroot+1))
            rootb(itrval)=dmax1(rootb(itrval),rootb(kroot+1))
c  worst case: suppose the kroot+1 interval overlaps with kroot-1 interval also!
c  then we overlap farther back, delete an interval
c  and update the number of roots expected (in jtrval) 
  543       continue
            if(rootb(itrval-1).gt.roota(itrval)) then
              itrval=itrval-1
              jtrval(itrval)=jtrval(itrval)+jtrval(itrval+1)
              jtrval(itrval+1)=mtry
              roota(itrval)=dmin1(roota(itrval),roota(itrval+1))
              rootb(itrval)=dmax1(rootb(itrval),rootb(itrval+1))
              go to 543
            endif
          else           
            itrval=itrval+1
            roota(itrval)=roota(kroot+1)
            rootb(itrval)=rootb(kroot+1)
          endif
        end do
        ntrval=itrval
c  LAST CHECK -- sometimes the last few intervals overlap
c  because the last interval near cmax is meant to catch new overtone branches
c  this loop coalesces the overlapping intervals
  345   continue
        if(rootb(ntrval-1).gt.roota(ntrval)) then
          ntrval=ntrval-1
          jtrval(ntrval)=jtrval(ntrval)+jtrval(ntrval+1)
c  need to take max and min because sometimes branches cross at finite angles
c  and the roota intervals pass through each other
          roota(ntrval)=dmin1(roota(ntrval),roota(ntrval+1))
          rootb(ntrval)=dmax1(rootb(ntrval),rootb(ntrval+1))
          go to 345
        endif
c        print *,ntrval
c        do i=1,ntrval
c          print *,roota(i),rootb(i)
c        end do
        f=f+dfrq
        om=2.d0*pi*f/ren
        nm=0
        do nt=1,ntrval
c  in each interval, first make an array of values to look for zeroes in
          nr=jtrval(nt)
          nxroot=nr/mtry
          if(nt.eq.ntrval) nxroot=nxroot-1
  100     nmm=0
c  there is the special case of the purely evanescent Rayleigh wave
c  for which cvel < cmin and cvel=gvel
c  we test for this case, and short-circuit the rootfinder if found
          if(nt.eq.1) then
            cc=cvel(nfrq+2-ifrq,1)/vbar
            gv=gvel(nfrq+2-ifrq,1)/vbar
            if(cc.lt.cmin.and.(cc/gv-1.d0).lt.1.d-4) then
              nmm=nmm+1
              nm=nm+1
              cvel(nfrq+1-ifrq,nm)=cc*vbar
              gvel(nfrq+1-ifrq,nm)=gv*vbar
              nc(nm)=nc(nm)+1
              cc0=cc
              go to 747
            endif
          endif
          dr=(rootb(nt)-roota(nt))/(nr-1)
          do jj=1,nr
            cc=roota(nt)+(jj-1)*dr
            call zzget(nl,cc,om,zz,ievan(jj))
            ccc(jj)=cc
            zzz(jj)=zz
            cc4(jj)=cc*ren/1000.d0
          end do
c          do i=1,nr
c            zz4(i)=real(zzz(i))
c          end do
c          call plotit(cc4,zz4,dum,nr,'Propagating Mode Condition',
c     x       'Phase velocity (km/sec)','Surface Traction',2,0,0,0,0)
c          do i=1,nr
c            zz4(i)=imag(zzz(i))
c          end do
c          call plotit(cc4,zz4,dum,nr,'Propagating Mode Condition',
c     x     'Phase velocity (km/sec)','Surface Traction',2,0,.05,0,1)
c  zero in on normal modes
          cc0=0.d0
          do ij=2,nr
            if(ievan(ij).eq.1) then      ! P-wave evanescent in top layer
              ick=0  			 ! hunt for zero in real(zzz)
              a1=dreal(zzz(ij-1))
              a2=dreal(zzz(ij))
            else
              a1=dimag(zzz(ij-1))
              a2=dimag(zzz(ij))
              aa1=dreal(zzz(ij-1))
              aa2=dreal(zzz(ij))
c  in practice either imag or real part often approaches zero parabolically
c  therefore, it makes sense to choose the most auspicious pair
              if(a1*a2.lt.0.d0.and.dabs(aa2).le.1.d0) then
                ick=1   		! hunt for zero in imag(zzz)
              elseif(aa1*aa2.lt.0.d0.and.dabs(a2).le.1.d0) then
                ick=0			! hunt for zero in real(zzz)
                a1=aa1
                a2=aa2
              endif
            endif
            if(a1*a2.lt.0.d0.and.dabs(dreal(zzz(ij))).le.1.d0) then
c  zoom in
              cc1=ccc(ij-1)
              cc2=ccc(ij)
              iter=0
              call zoom_in(nl,a1,a2,cc1,cc2,om,iter,cc,zz,ick)
c  save if no flag, and if the root was not fortuitously at a grid point
              if(iter.ge.0.and.cc.ne.cc0) then
                call grvel(nl,cc,om,gv)
c                print *,(vertr(kji),kji=1,6)
c                print *,(verti(kji),kji=1,6)
c  last check: horizontal waves in the top layer will have cc \appox gv
c  these are spurious, so kill them
c  the test is based on the propagator in the top layer
c  a horizontal wave in the top layer will have complex propagtor = 1.d0
c  the arrays vertr and verti are filled in the subroutine zzget
                iflag=0
                do kji=1,6
                  test=dabs(vertr(kji)-1.d0)+dabs(verti(kji))
                  if(test.lt.1.d-6) iflag=1
                end do
c                if(dabs(cc/cmin-1.d0).gt.1.d-10) then
                if(iflag.eq.0) then
                  nmm=nmm+1
                  nm=nm+1
                  cvel(nfrq+1-ifrq,nm)=cc*vbar
                  gvel(nfrq+1-ifrq,nm)=gv*vbar
                  nc(nm)=nc(nm)+1
                  cc0=cc
                endif
              endif
            endif
          end do
          if(nmm.ne.nxroot) then
            print *,'at f=',f,'Hz'
            ra=roota(nt)*vbar
            rb=rootb(nt)*vbar
            print *,' in cvel interval',ra,rb
            print *,nxroot,' expected, but',nmm,' roots found'
            nr=2.0*nr
c  except, dont iterate forever on cmin
            if(cmin.gt.roota(nt).and.cmin.lt.rootb(nt)) then
              if(nr.gt.100) then  ! trigger the if-then loop to exit
                nr=9001
              endif
            endif
            if(nmm.lt.nxroot.and.nr.lt.9000) then ! run that play again
              nm=nm-nmm
              if(nmm.gt.0) then
                do i=1,nmm
                  nc(nm+i)=nc(nm+i)-1
                end do
              endif
c  expand the search area a bit
              droot=rootb(nt)-roota(nt)
              roota(nt)=roota(nt)-0.05d0*droot
              rootb(nt)=rootb(nt)+0.05d0*droot
              if(nt.gt.1) then
                roota(nt)=dmax1(roota(nt),rootb(nt-1))
              endif
              if(nt.lt.ntrval) then
                rootb(nt)=dmin1(rootb(nt),roota(nt+1))
              else
                rootb(nt)=dmin1(rootb(nt),cmax-eps)
              endif
              go to 100
            elseif(nmm.lt.nxroot) then
              print *,'failed!'
            endif
          endif
  747     continue   ! move to next interval if fundamental R is taken
        end do
        kroots(ifrq)=nm
        print *,kroots(ifrq),' = nroot,', ntrval,' = ntrval  at f=',f
        do j=1,nm
          zz4(j)=cvel(nfrq+1-ifrq,j)/1000.
        end do
        write(7,1010) (f,zz4(j),j=1,nm)
        do j=1,nm
          zz4(j)=gvel(nfrq+1-ifrq,j)/1000.
        end do
        write(8,1010) (f,zz4(j),j=1,nm)
c        if(ifrq.eq.(ifrq/50)*50) print *,nm,' modes at',f,' Hz'
      end do
      close(7)
      close(8)
  101 format(6g15.5)
  105 format('Modal Eigenfunction, f=',f3.1,' Hz,  c=',f7.4,' km/s')
c  after the cycle through frequencies, we plot the dispersion curves
c      x4(1)=0.
c      x4(2)=frqmax
c      x4(3)=vss(1)*0.89*vbar/1000.d0
c      x4(4)=vss(nlp)*1.01*vbar/1000.d0
      cmn=cmin*vbar
      cmx=cmax*vbar
      open(7,file='aniprop_out.dat',form='unformatted')
c  first write out the model
      write(7) title
      write(7) nl,maxbr
      do i=1,nlp
        write(7) (w(j,i),j=1,3),z(i),vp(i),vp2(i),
     x                               vp4(i),vs(i),vs2(i),rho(i)
      end do
c  then the phase and group velocities
      write(7) nfrq,frqmax,cmn,cmx,(nc(i),i=1,maxbr)
      print *,'output: nfrq,frqmax,cmn,cmx', nfrq,frqmax,cmn,cmx
      print *,(nc(i),i=1,10)
      write(6,1099)(nc(i),i=1,maxbr)
 1099 format(20i5)
      do imode=maxbr,1,-1
        if(nc(imode).gt.1) then
c          write(7) (cvel(i,imode),i=1,nc(imode))
          write(7) (cvel(i,imode),i=1,nfrq)
        endif
      end do
      do imode=maxbr,1,-1
        if(nc(imode).gt.1) then
c          write(7) (gvel(i,imode),i=1,nc(imode))
          write(7) (gvel(i,imode),i=1,nfrq)
        endif
      end do
      write(7) (kroots(i),i=1,nfrq)
      close(7)
c      do ii=1,nfrq
c        i=nfrq+1-ii
c        f=ii*dfrq
c        write(6,1010) f,(cvel(i,j),j=1,10)
c        write(6,1011) (gvel(i,j),j=1,10)
c      end do
 1010 format(f5.3,f10.6)
 1011 format(5x,10f10.3)
      call plotit_axes(x4(1),x4(2),x4(3),x4(4))
      do imode=maxbr,1,-1
        if(nc(imode).gt.1) then
          np=nc(imode)
          do i=1,np
            zz4(i)=gvel(i,imode)/1000.    ! copy into real*4 array for plotting
          end do
          call plotit(frqq,zz4,dum,np,title,'Frequency (Hz)',
     x        'Group Velocity (km/s)',2,0,0.05,2,1/imode)
        endif
      end do
      x4(3)=vss(1)*vbar*0.79/1000.d0
      call plotit_axes(x4(1),x4(2),x4(3),x4(4))
      do imode=maxbr,1,-1
        if(nc(imode).gt.1) then
          np=nc(imode)
          do i=1,np
            zz4(i)=cvel(i,imode)/1000.    ! copy into real*4 array for plotting
          end do
          call plotit(frqq,zz4,dum,np,title,'Frequency (Hz)',
     x        'Phase Velocity (km/s)',2,0,0.05,2,1/imode)
        endif
      end do
c  trapdoor for writing sac files to disk for later processing
c  read in a sac file for its header
c   99     print *,'enter SAC-format file for header info - bh? file'
c          read(5,101) name
c      name='split.0.0.bhz'
      name='/Users/jjpark/Courses/GEO557/T05C.BHZ'
      print *,'need a SAC-format data file '
      print *,'as template for SAC-format output'
      call sacread(name,u4,ierr)
c          if(ierr.ne.0) go to 99    ! in case of mistyped filename
      if(ierr.ne.0) then    ! in case of mistyped filename
        print *,'cant find the SAC-format file ',name
        stop
      endif
c     call plotit_axes(x4(1),x4(2),x4(3),x4(4))
c  after seeing the dispersion curves, lets sum up some synthetics
c  loop over frequency.
c  we assume that rhat=xhat, thetahat=yhat, and this simplifies the formulas
c  after summing the response in the freq domain, we inverse fft
      print *,'for EW strike-slip fault: 0,0,0,0,0,1'
      print *,'for NE-SW strike-slip fault: 1,-1,0,0,0,0'
      print *,'for NS striking 45-dip thrust/normal: 1,0,-1,0,0,0'
      print *,'for EW striking 45-dip thrust/normal: 0,1,-1,0,0,0'
      print *,'for EW striking 90-dip thrust/normal: 0,0,0,0,1,0'
      print *,'for NS striking 90-dip thrust/normal: 0,0,0,1,0,0'
      print *,'enter Mxx,Myy,Mzz,Mxy,Mxz,Myz, times 10**17 nt-m'
      read(5,*) amom(1,1),amom(2,2),amom(3,3),
     x                    amom(1,2),amom(1,3),amom(2,3)
c  for starts lets fix the source, phi=0
c      amom(1,1)=1.d17/ebar
c      amom(2,2)=-1.d17/ebar
c      amom(3,3)=0.d0
c      amom(1,2)=1.d17/ebar
c      amom(1,3)=0.d0
c      amom(2,3)=0.d0
      amom(2,1)=amom(1,2)
      amom(3,1)=amom(1,3)
      amom(3,2)=amom(2,3)
      do i=1,3
        do j=1,3
          amom(i,j)=amom(i,j)*1.d17/ebar
        end do
      end do
      xx=500000./radi   ! 500 km
      h=15000./radi     ! depth is 15 km to start
      print *,'enter downrange distance x and depth h in km'
      read(5,*) xx,h
      xx=xx*1000./radi
      h=h*1000./radi
  300 continue
      il=1
      do n=1,nlp
        if(h.gt.z(n)) il=il+1
      end do
      if(il.eq.1) then
        h1=0.d0
      else
        h1=z(il-1)
      endif
      h2=z(il)
      print *,il,h1,h2,' =layer #, bounds'
      df=frqmax/nfrq
      f=0.d0
      print *,'the maximum frequency of mode set is ',frqmax
      print *,'enter maximum frequency for synthetic '
      read(5,*) freqmx
      freqmx=dmin1(freqmx,frqmax)
      mfrq=freqmx/df      
 1009 format('frequency ',f7.4)
c      print *,'WE ARE ONLY INCLUDING FUNDAMENTAL MODES'
c  loop over frequencies
      do jf=1,mfrq
        do i=1,3
          uu(i)=z0
        end do
        f=f+df
        if(100*(jf/100).eq.jf) print 1009,f
        om=2.d0*pi*f/ren
c  loop over overtones at a certain frequency
        do iov=1,maxbr
c        do iov=1,2    
          cc=cvel(nfrq+1-jf,iov)/vbar
          gv=gvel(nfrq+1-jf,iov)/vbar
          if(cc.gt.tol) then  ! skip if no mode (cvel=0)
c          if(cc.gt.tol.and.gv.lt.gvtol) then  ! skip if no mode (cvel=0)
c          print *,cvel(jf,iov),gvel(jf,iov)
c  calc phase factor modulo 2pi
            akx=om/cc
            fac=akx*xx+pi/4.d0
            fac=fac-int(fac/(2.d0*pi))*2.d0*pi
            spread=dsqrt(2.d0/(pi*akx*xx))
            zz0=dcmplx(dcos(fac),dsin(fac))
c            print *,f,cc,om,akx,nl
            call zzget(nl,cc,om,zz,iev)
c  we use grvel subroutine to obtain the kinetic energy functional
c  which is the normalization factor
            call grvel1(nl,cc,om,ci1,gvx,gvy)
c            var2(jf,iov)=datan2d(gvy,gvx)
c            var4(jf,iov)=gvy/gvx
c            print *,f,cc,zz,gv,ci1
            cccc=cc*vbar
            if(zabs(zz).gt.0.1d0.and.cccc.gt.cmin) then
              print *,'is this root right?',f,cccc,zz
c             pause
            endif
c  calculate the eigenfunction at the surface and at source depth h
            call get_eigen(akx,f,nl,h,usx,usz,ur)
c  eigenmodes are calculated in a coordinate system with z increasing downward
c  and y increasing "out of" the board
c  to restore to x,y,z = radial,transverse,vertical, flip sign of y&z comps
            ur(2)=-ur(2)
            ur(3)=-ur(3)
            usx(2)=-usx(2)
            usx(3)=-usx(3)
            usz(2)=-usz(2)
            usz(3)=-usz(3)
c  calculate the SH amplitude at surface, relative to total amplitude
            sum=0.d0
            do i=1,3
              sum=sum+zabs(ur(i))**2
            end do
c            var1(jf,iov)=zabs(ur(2))/dsqrt(sum)
c  doubledot moment tensor with source strain (conjugated)
            zz=z0
            do i=1,3
              zz=zz+amom(1,i)*conjg(usx(i))
              zz=zz+amom(3,i)*conjg(usz(i))
            end do
            zz=spread*zz/(8.d0*cc*gv*ci1)
            zz=zz*zz0
            do i=1,3
              uu(i)=uu(i)+ur(i)*zz
            end do
          endif
        end do     
c  we pack the fourier coefficients into real array for inverse fft
        do i=1,3
          uuu(2*jf+1,i)=dreal(uu(i))    
          uuu(2*jf+2,i)=dimag(uu(i))    
          dum4(jf,i)=zabs(uu(i))
        end do
      end do
      print *,'done with frequency loop'
c  taper the spectrum near the freq cutoff
      nf1=3*mfrq/4
      do i=nf1,mfrq
        fac=dcos(0.5d0*pi*(i-nf1)/dfloat(mfrq-nf1))**2
        do j=1,3
          uuu(2*i+2,j)=uuu(2*i+2,j)*fac
          uuu(2*i+1,j)=uuu(2*i+1,j)*fac
          dum4(i,j)=dum4(i,j)*fac
        end do
      end do       
c  plot the spectrum
      npad=16384
      call plotit_axes(0.,0.,0.,0.)
      call manyplot(npad,1,npad,3,frqq,dum4,1.,
     x				'Spectrum','frq(Hz)',' ',)
      do i=1,3
        call plotit_axes(0.,0.,0.,0.)
        call plotit(frqq,dum4(1,i),dum,mfrq,'Spectrum','frq(Hz)',' ',
     x                                                 2,0,0.0,0,30+i)
      end do
c  clean up the details - zero at zero-freq, zero-pad, inverse fft
        do i=1,3
          uuu(1,i)=0.0
          uuu(2,i)=0.
          do j=mfrq*2+3,npad+2
            uuu(j,i)=0.
          end do
	  jj=0
          do j=1,npad,2
	    jj=jj+1
            uur(jj)=uuu(j,i)
	    uui(jj)=-uuu(j+1,i)
            uur(npad+2-jj)=uuu(j,i)
	    uui(npad+2-jj)=uuu(j+1,i)
          end do          
	  uur(npad/2+1)=0.
	  uui(npad/2+1)=0.
	  print *,'entering inverse fft',npad
          call fft2(uur,uui,npad)
	  print *,'leaving inverse fft',npad
c  return to displacement units (micro-meters)
          do j=1,npad
            uuu(j,i)=uur(j)*radi*1.d6/npad
          end do
	end do
      x4(1)=0.
      fny=(npad/2)*frqmax/dfloat(nfrq)
      dt=0.5/fny
      tmax=dt*npad
      print *,'a little farther'
c      call manyplot(npad,1,npad,3,x4,u,1.,title,'time(sec)','amp')
c      call plotit_axes(0.,0.,0.,0.)
c      do i=1,3
c        ii=(i-1)*1500
c        do j=1,1500
c          u4(ii+j)=uuu(j*3-2,i)
c        end do
c        x4(2)=3.*0.5/fny
c        call plotit(x4,uuu(1,i),dum,1500,title,'time(sec)',
c     x     'Displacement (microns)',1,0,0,0,30+i)
c      end do
c      call manyplot(1500,1,1500,3,x4,u4,1.,title,'time(sec)','amp')
c  for later analysis, write out the seismic records
c  but first splinefit it to 20 sps
      nscan=10000
      ahead(1)=0.05
      ahead(53)=270.  
c wave comes from the west from lat=0,lon=0 to lat=0, lon=lon1
c /bin/mv syn.z SYN.BHZ
c /bin/mv syn.r SYN.BHE
c /bin/mv syn.t SYN.BHN
      iah(80)=nscan
      ahead(36)=0.
      ahead(37)=0.
      ahead(39)=h*radi
      ahead(51)=xx*radi/1000.
      print *,ahead(39)
      ahead(32)=0.
c  lengths in code are normalized by earth radius 
c  mult by one deg/radian to get DELTA
      ahead(33)=xx*radian
      ahead(34)=0.
      print *,ahead(39)
      ahead(54)=ahead(33)
      ddtt=0.05/dt
      print *,'time interval in spline:',ddtt
      do j=1,3
        call splneq(npad,uuu(1,j),dum4(1,j))
        chead(151)=ichan(j)
        time=0.
        do i=1,nscan
          time=time+ddtt
          u4(i)=evaleq(time,npad,uuu(1,j),dum4(1,j),0,1.)
        end do
        name=outfile(j)
	print *,name
        call sacout(name,u4)
      end do
  220 print *,'max time = ',tmax
      print *,'enter the desired time interval '
      read(5,*) x1,x2
      i1=x1/dt+1
      i2=x2/dt+1
      np=i2-i1+1
      inc=np/1500+1
      x4(1)=(i1-1)*dt
      x4(2)=inc*dt
      ii=0
      do i=1,3
        do j=i1,i2,inc
          ii=ii+1
          u4(ii)=uuu(j,i)
        end do
      end do
      npp=ii/3
      call manyplot(npp,1,npp,3,x4,u4,1.,title,'time(sec)','amp')
c      do i=1,3
c        call plotit_axes(x1,x2,0.,0.)
c        call plotit(x4,uuu(1,i),dum,5000,title,'time(sec)',
c     x     'Displacement (microns)',1,0,0,0,30+i)
c      end do
      print *,'another interval (0=no)'
      read(5,*) ick
      if(ick.ne.0) go to 220
      print *,'enter new distance and depth in km (negative = no)'
      read(5,*) xx,h
      xx=xx*1000.d0/radi
      h=h*1000.d0/radi
      if(xx.gt.0.d0) go to 300


      stop
      end
      subroutine zoom_in(nl,a1,a2,cc1,cc2,om,iter,cc,zz,ick)
c  ick = 1 straddle a zero in imag(zz)
c  ick = 2 straddle a zero in real(zz)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      complex*16 zz
      data pi/3.14159265358979d0/,eps/1.d-8/,tol1/2.d-2/,tol2/1.d-5/
      do while(dmin1(dabs(a1),dabs(a2)).gt.eps.and.iter.lt.100)
        iter=iter+1
        if(dabs(a1)+dabs(a2).gt.tol1) then
c bisection - slow and steady
          frac=0.5d0
        else
c  secant, use when close to the root
          frac=dabs(a1)/(dabs(a1)+dabs(a2))
c  to avoid secant pivot for inauspicious curvature
          if(frac.lt.0.01d0) frac=dmax1(2.d0*frac,0.01d0)
          if(1.d0-frac.lt.0.01d0) frac=dmin1(2.d0*frac-1.d0,0.99d0)
        endif
        cc=cc1+frac*(cc2-cc1)
        call zzget(nl,cc,om,zz,iev)
        if(ick.eq.1) then
          a3=dimag(zz)
        else
          a3=dreal(zz)
        endif
        if(a3*a1.le.0.d0) then
          a2=a3
          cc2=cc
        else
          a1=a3
          cc1=cc
        endif
      end do
      if(dabs(a1).lt.dabs(a2)) then
        cc=cc1
      else
        cc=cc2 
      endif
c  check that the real part is zero
      if(zabs(zz).lt.tol2) then
        if(iter.ge.99) print *,cc,iter,zz
      else
c  flag to discard result, but test again to be certain 
c  if one a gridpoint was the root, there were no iterations, and zz is crap
c        print *,'failed?',cc,iter,zz
        call zzget(nl,cc,om,zz,iev)
c        print *,'recalculated',cc,zz
        if(dabs(dreal(zz))+dabs(dimag(zz)).gt.tol2) then
          iter=-1
c          print *,'failed again',cc,zz
        endif
      endif
      return
      end
      subroutine zzget(nl,cc,om,zz,iev)
c  returns determinant zz of a stack of anisotropic layers
c  iev=1 if the waves are evanescent in the top layer, iev=0 otherwise
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      complex*16 pp,u0,ee,pw,uw,pu,z1,z0,zz,xnu,eye,e1,e2,zla,rtm
      complex*16 rt,tt,rt0,trc,xl,pfac,u
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
      common/qstff/qi(6,6),xr(6),xi(6),yr(6),yi(6)
      common/rstff/rt(3,3,101),tt(3,3,101),rt0(3,3),trc(3,3)
      common/vertwav/vertr(6),verti(6)
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
c   requires an eigenvector problem be solved
c  in general, the evanescent vertical wavenumbers have nonzero real parts 
c   complex exponential fct is used to avoid endless branching
c  horizontal slowness p
      akx=om/cc
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
            t(j,i)=t(j,i)*akx*akx
            s(j,i)=s(j,i)*akx
          end do
          t(i,i)=t(i,i)-om*om*rho(n)
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
c  solve eigenvalue problem, nonsymmetric real valued
c  from the eispack guide
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
c        print *, 'for phase velocity',cc,'  the vertical wavenumbers are'
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
        pp(1)=dcmplx(akx,0.d0)
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
c  so downgoing oscillatory waves have k_z>0, k_z real
c  downgoing evanescent waves have Im(k_z)>0
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
c  incorporate the factor of i into the stored vertical wavenumber 
          xnu(ki,n)=dcmplx(-wi(k),wr(k))
        end do
c        do i=1,6
c          print *,'for i*k_z:',xnu(i,n),', the disp-stress vector is'
c          do j=1,6
c            xi(j)=dimag(ee(j,i,n))
c            xr(j)=dreal(ee(j,i,n))
c          end do
c          print 101,(xr(j),j=1,6),(xi(j),j=1,6)
c        end do
c        pause
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
c              print 1008,'vert wavenumbers',xnu(i,n),' and',xnu(j,n)
c  we average eigenvalues (vert wavenumbers)
c  and solve for eigenvector in subroutine defective
              xnu(i,n)=(xnu(i,n)+xnu(j,n))/2.d0
              xnu(j,n)=xnu(i,n)
c              zz=zz/(zabs(zz))
c              sq2=dsqrt(2.d0)
c              do jj=1,6
c                ee(jj,i,n)=(zz*ee(jj,i,n)+ee(jj,j,n))/sq2
c              end do
c  calculate the extravector for defective repeated eigenvalue
              call defective(i,j,n,adf(1,n),a,b,c,d,e,om,akx)
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
c              print 1008,'vert wavenumbers',xnu(i,n),' and',xnu(j,n)
c  we average the eigenvalues
              xnu(i,n)=(xnu(i,n)+xnu(j,n))/2.d0
              xnu(j,n)=xnu(i,n)
c             zz=zz/(zabs(zz))
c              sq2=dsqrt(2.d0)
c              do jj=1,6
c                ee(jj,i,n)=(zz*ee(jj,i,n)+ee(jj,j,n))/sq2
c              end do
c  calculate the extravector for defective repeated eigenvalue
c  as well as coefficient (adf) of linear (z-z0) term
              call defective(i,j,n,adf(2,n),a,b,c,d,e,om,akx)
c              print *,i,j,n,ccc,qqq,adf(2,n)
            endif
          end do
        end do
      end do
c  calculate modified R/T coefficients
c  first calc the propagation factors
c  note that for dipping fast axes the upgoing and downgoing wavenumbers are
c  independent, so we must calc all to be safe
      do n=1,nl
        do k=1,3
          xl(k,n)=zexp(xnu(k,n)*dz(n))		! downgoing
          xl(k+3,n)=zexp(-xnu(k+3,n)*dz(n))	! upgoing
        end do
      end do
c  ok, lets keep track of the propagation factors in the top layer
c  we wish to test whether we have a horizontally propagating S or P wave
c  this will have a phase velocity = S or P velocity in layer
c  and NOT be  a properly trapped resonance
      do k=1,6
        vertr(k)=dreal(xl(k,1))
        verti(k)=dimag(xl(k,1))
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
            zla(k+3)=0.d0
          endif
        end do
c mult the columns of e2
        do k=1,6
          do i=1,6
            e2(i,k)=e2(i,k)*zla(k)
          end do
        end do
c  the possibility of defective matrices must be contemplated here
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
c  recursive calc of generalized R/T coefs:
c  at the last interface:
      do k=1,3
        do i=1,3
          tt(i,k,nl)=rtm(i,k,nl)
          rt(i,k,nl)=rtm(i+3,k,nl)
        end do
      end do  
c  interate:
      do n=nlm,1,-1
c  first the generalized transmission coef:
        do k=1,3
          do i=1,3
            rt0(i,k)=z0
            do j=1,3
              rt0(i,k)=rt0(i,k)-rtm(i,j+3,n)*rt(j,k,n+1)
            end do
          end do
          rt0(k,k)=rt0(k,k)+z1
        end do  
        do k=1,3
          do i=1,3
            s(i,k)=dreal(rt0(i,k))
            t(i,k)=dimag(rt0(i,k))
          end do
        end do
        nn=3
        do k=1,3
          do i=1,3
            yr(i)=dreal(rtm(i,k,n))
            yi(i)=dimag(rtm(i,k,n))
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
            rt0(i,k)=z0
            do j=1,3
              rt0(i,k)=rt0(i,k)+rt(i,j,n+1)*tt(j,k,n)
            end do
          end do
        end do  
        do k=1,3
          do i=1,3
            rt(i,k,n)=rtm(i+3,k,n)
            do j=1,3
              rt(i,k,n)=rt(i,k,n)+rtm(i+3,j+3,n)*rt0(j,k)
            end do
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
c  the possibility of defective matrices must be contemplated here
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
c  calc the traction matrix
      do k=1,3
        do i=1,3
          trc(i,k)=z0
          do j=1,3
            trc(i,k)=trc(i,k)-rt0(i,j)*rt(j,k,1)
          end do
        end do
        trc(k,k)=trc(k,k)+z1
      end do  
c the determinant of the traction
      zz=trc(1,1)*(trc(2,2)*trc(3,3)-trc(3,2)*trc(2,3))
     x  +trc(1,2)*(trc(2,3)*trc(3,1)-trc(3,3)*trc(2,1))
     x  +trc(1,3)*(trc(2,1)*trc(3,2)-trc(3,1)*trc(2,2))
      return
      end
      subroutine defective(i,j,n,adf,a,b,c,d,e,om,akx)
c  patch for dealing with nearly defective propagator matrices
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
      common/qstff/qi(6,6),xr(6),xi(6),yr(6),yi(6)
      z1=dcmplx(1.d0,0.d0)
      z0=dcmplx(0.d0,0.d0)
      eye=dcmplx(0.d0,1.d0)
c  for the extravector, need to solve system of equations
c  based on original 6x6 Q matrix
c the plane-wave solutions generalize to the form
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
c  we wish to find the eigenvector of the near-defective matrix 
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
      pp(1)=dcmplx(akx,0.d0)
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
      pp(1)=dcmplx(akx,0.d0)
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
      subroutine grvel(nl,cc,om,gv)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      real*4 cc4,zz4,x4
      complex*16 pp,u0,ee,z1,z0,zz,xnu,zzz,e1,e2,zla,xl,eigf
      complex*16 rt,tt,rt0,trc,co,ci1,ci2,ci3,ci4,dci1,dci2,dci3,dci4
      complex*16 u,xxnu,eye,pfac,uu,wukk,wuk,dci31,dci32
      common/stfff/w(3,101),t(3,3),ttl(3,3),s(3,3),stl(3,3),
     x                                           r(3,3),x(3),y(3)
      common/model/z(100),dz(100),rho(101),vp(101),vs(101),vp2(101),
     x               vp4(101),vs2(101),vss(101)
      common/model2/xmu(101),xla(101),xmu2(101),xla2(101),xla4(101)
      common/propag/xnu(6,101),xl(6,100),pfac(6,3),u(3,6)
      common/mstff/qq(6,6),wr(6),wi(6),zr(6,6),zi(6,6),iv(6),fv(6)
      common/pstff/pp(3),u0(3),ee(6,6,101),e1(6,6),e2(6,6),zla(6)
      common/rstff/rt(3,3,101),tt(3,3,101),rt0(3,3),trc(3,3)
      common/nstff/cc4(10000),zz4(10000),zzz(10000),ccc(10000),x4(10000)
      common/pstf2/co(6,101),eigf(1000,3)
      common/cmstf/ci1,ci2,ci3,ci4,dci1,dci2,dci3,dci4,wukk,wuk,
     x             dci31,dci32,xxnu,eye,uu
      data pi/3.14159265358979d0/,eps/1.d-6/,tol/1.d-3/
      z1=dcmplx(1.d0,0.d0)
      z0=dcmplx(0.d0,0.d0)
      eye=dcmplx(0.d0,1.d0)
      ren=1.075190645d-3
      radi=6.371d6
      vbar=ren*radi
c      print *,nl,cc,om
      nlp=nl+1
c  update all the layer calculations
      call zzget(nl,cc,om,zz,iev)
c characteristics of the solution in the first layer
c            do k=1,6
c              print *,'for slowness:',xnu(k,1),', disp-stress vector'
c              print 103,(ee(j,k,1),j=1,6)
c            end do
  103 format(2g15.5)
c  calc the eigenfunction
c  play safe, since a 3x3 system of equations A.x=0 and detA=0
c  x is propto the normalized crossproduct of any two nonparallel rows of A
c  this condition gets dicey if we have degeneracy, because the 
c  determinant would have a double root, and all rows would be parallel
c   or when Rayleigh and Love decouple
c  so the safest route is to take the eigenvector of A with smallest eigenvalue
c  easiest variant: take Hermitian A^H.A, smallest (real) eigenvalue wins
c  in practice, gvel gets hybridized between dispersion branches 
c  when Love and Rayleigh couple
      do i=1,3
        do j=1,3
          rt0(j,i)=z0
          do k=1,3
            rt0(j,i)=rt0(j,i)+conjg(trc(k,j))*trc(k,i)
          end do
          s(j,i)=dreal(rt0(j,i))
          stl(j,i)=dimag(rt0(j,i))
        end do
      end do
      call htridi(3,3,s,stl,wi,fv,fv,wr)
      do i=1,3
        do j=1,3
          t(j,i)=0.d0
        end do
        t(i,i)=1.d0
      end do
      call tql2(3,3,wi,fv,t,ierr)
      if(ierr.ne.0) then
        print *,'tql2 error',ierr
        stop
      endif
      call htribk(3,3,s,stl,wr,3,t,ttl) 
      do i=1,3
        co(i,1)=dcmplx(t(i,1),ttl(i,1))
        co(i+3,nlp)=z0
      end do
      do n=1,nl
c  upgoing coefs in the nth layer, downgoing coefs in the n+1 layer
        do i=1,3
          co(i+3,n)=z0
          co(i,n+1)=z0
          do j=1,3
            co(i+3,n)=co(i+3,n)+rt(i,j,n)*co(j,n)
            co(i,n+1)=co(i,n+1)+tt(i,j,n)*co(j,n)
          end do
        end do
      end do       
c  calculate the group velocity
c  rather than integrate numerically, 
c  we use the Lagrangian interaction matrices 
c  for sums of plane waves in each layer
c  cross terms are complex valued
      ci1=z0
      ci2=z0
      ci3=z0
      ci4=z0
      do il=1,nlp
        if(il.eq.1) then
          h1=0.d0
        else
          h1=z(il-1)
        endif
        if(il.eq.nlp) then ! k=4,5,6 upgoing waves are zero
          nok=3
        else
          nok=6
        endif
        h2=z(il)
        a=xla(il)
        b=xla2(il)
        c=xla4(il)
        d=xmu(il)
        e=xmu2(il)
c  assemble the six polarization vectors and propagation factors
        do k=1,nok
c
c          pfac(k,1)=zexp(xnu(k,il)*h1)
c          if(il.lt.nlp) then
c            pfac(k,2)=zexp(xnu(k,il)*h2)
c          else
c            pfac(k,2)=z0
c          endif
c  normalize for the coefficient amplitudes in the Chen algorithm
c  upgoing coefs are referenced to layer bottom
c  downgoing coefs are referenced to layer top
c          pfac(k,3)=pfac(k,1+(k-1)/3)
c
C  START ALTERNATE COMPUTATION, HOPEFULLY MORE STABLE
c  normalize for the coefficient amplitudes in the Chen algorithm
c  during this computation, rather than dividing by extreme numbers later
c  upgoing coefs are referenced to layer bottom
c  downgoing coefs are referenced to layer top
c
c  we have a test to avoid underflows and overflows in the evenescent regime
c  if the decay factor in a layer is 10**(-50) or 10**(50) 
c  the eigenfunction in the layer cannot truly contribute to integrals
c  so we zero out the relevant pfac(.,1 or 2) pfac(.,3) propagation factors
c
c  pfac(k,3) has exponential negated in order to change  
c  division to multiplication in its use in the code below.  
c  This way we can set pfac(k,3)=0. when underflow is imminent
          test=dreal(xnu(k,il)*(h2-h1))
          iftest=1
          if(dabs(test).gt.100.d0) iftest=0
          if(k.le.3) then
            pfac(k,1)=z1
            if(il.lt.nlp.and.iftest.eq.1) then
              pfac(k,2)=zexp(xnu(k,il)*(h2-h1))
              pfac(k,3)=zexp(-xnu(k,il)*h1)
            else
              pfac(k,2)=z0
              pfac(k,3)=z0
            endif
          else
            if(iftest.eq.1) then
              pfac(k,1)=zexp(xnu(k,il)*(h1-h2))
              pfac(k,3)=zexp(-xnu(k,il)*h2)
            else
              pfac(k,1)=z0
              pfac(k,3)=z0
            endif
            if(il.lt.nlp) then
              pfac(k,2)=z1
            else
              pfac(k,2)=z0
            endif
          endif
C  END ALTERNATE COMPUTATION, HOPEFULLY MORE STABLE
          do i=1,3
            u(i,k)=co(k,il)*ee(i,k,il)
          end do
        end do
c  form matrix of vertical integrals -- check the coef convention!
c  kth wave is conjugated, kkth wave isnt
        do k=1,nok
c  note that we ADD conjg(xnu(k,il)) because xnu has a factor of i
          do kk=k,nok
            xxnu=xnu(kk,il)+conjg(xnu(k,il))
c            cvel=cc*vbar/1000.
c            if(cvel.lt.3.0d0.and.nlp-il.le.1) then
c              print *,k,kk
c              print *,'xxnu,xnu(kk),xnu(k)',xxnu,xnu(kk,il),xnu(k,il)
c              print *,'(pfac(k,i),i=1,3)',(pfac(k,i),i=1,3)
c              print *,'h1,h2',h1,h2
c              print *,'(pfac(kk,i),i=1,3)',(pfac(kk,i),i=1,3)
c            endif
            if(zabs(xxnu).lt.eps) then
c  true if propagating wave with zero vertical wavenumber (nu)
c  OR is integral of upgoing+downgoing evanescent waves with opposite nu
              xxnu=(h2-h1)/2.d0
c  normalize for the coefficient amplitudes in the Chen algorithm
c  upgoing coefs are referenced to layer bottom
c  downgoing coefs are referenced to layer top
c  taking a negative exponent in computation of pfac(k,3) allows us
c  to express this normalization as a multiplication
c  this way "zero=(0.,0.)" can be substituted in the far-evanescent region
              xxnu=xxnu*(pfac(kk,3)*conjg(pfac(k,3)))
            else
              xxnu=(pfac(kk,2)*conjg(pfac(k,2))
     x           -pfac(kk,1)*conjg(pfac(k,1)))/(2.d0*xxnu)
            endif
c  sum the plane wave contribution  
            uu=z0
            wukk=z0
            wuk=z0
            do i=1,3
              uu=uu+conjg(u(i,k))*u(i,kk)
              wuk=wuk+w(i,il)*conjg(u(i,k))
              wukk=wukk+w(i,il)*u(i,kk)
            end do
            wx=w(1,il)
            wz=w(3,il)
      dci1=rho(il)*uu
      dci2=(a-b+c-d+e)*conjg(u(1,k))*u(1,kk)
     x    +(d-e+2.d0*e*wx**2)*uu+(8.d0*c*wx**2+2.d0*e)*wuk*wukk
     x    +(b-4.d0*c-2.d0*e)*wx*(conjg(u(1,k))*wukk+wuk*u(1,kk))
      dci31=(a-b+c-2.d0*(d-e))*conjg(u(3,k))*u(1,kk)
     x    +2.d0*e*wx*wz*uu+8.d0*c*wx*wz*wuk*wukk
     x    +(b-4.d0*c-4.d0*e)*(wx*conjg(u(3,k))*wukk+wz*wuk*u(1,kk))
     x    +(d-e)*conjg(u(1,k))*u(3,kk)
     x    +2.d0*e*(wz*conjg(u(1,k))*wukk+wx*wuk*u(3,kk))
      dci32=(a-b+c-2.d0*(d-e))*conjg(u(1,k))*u(3,kk)
     x    +2.d0*e*wx*wz*uu+8.d0*c*wx*wz*wuk*wukk
     x    +(b-4.d0*c-4.d0*e)*(wz*conjg(u(1,k))*wukk+wx*wuk*u(3,kk))
     x    +(d-e)*conjg(u(3,k))*u(1,kk)
     x    +2.d0*e*(wx*conjg(u(3,k))*wukk+wz*wuk*u(1,kk))
c  note that we ADD conjg(xnu(k,il)) because xnu has a factor of i  (eye)
      dci3=eye*(conjg(xnu(k,il))*dci31-xnu(kk,il)*dci32)
      dci4=(a-b+c-d+e)*conjg(u(3,k))*u(3,kk)
     x    +(d-e+2.d0*e*wz**2)*uu+(8.d0*c*wz**2+2.d0*e)*wuk*wukk
     x    +(b-4.d0*c-2.d0*e)*wz*(conjg(u(3,k))*wukk+wuk*u(3,kk))
c  the extra factor of i in xnu cancels in this product
      dci4=dci4*xnu(kk,il)*conjg(xnu(k,il))
c      print *,k,kk,uu,xxnu,a,d
c      print *,'dci31',dci31
c      print *,'dci32',dci32
c      print *,'dci3',dci3
c      print *,'dci4',dci4
c     cvel=cc*vbar/1000.
c      if(cvel.lt.3.0d0.and.nlp-il.le.1) then
c        print *,k,kk
c        print *,ci1,ci2,ci3
c        print *,dci1,dci2,dci3
c        print *,xxnu,uu,wuk,wukk
c      endif
            if(k.eq.kk) then
              ci1=ci1+dci1*xxnu
              ci2=ci2+dci2*xxnu
              ci3=ci3+dci3*xxnu
              ci4=ci4+dci4*xxnu
            else
              ci1=ci1+2.d0*dreal(xxnu*dci1)
              ci2=ci2+2.d0*dreal(xxnu*dci2)
              ci3=ci3+2.d0*dreal(xxnu*dci3)
              ci4=ci4+2.d0*dreal(xxnu*dci4)
            endif
          end do
        end do
      end do
c      print *,'ci1',ci1
c      print *,'ci2',ci2
c      print *,'ci3',ci3
c      print *,'ci4',ci4
      gv=(ci2/cc+ci3/(2.d0*om))/ci1
      cvel=cc*vbar/1000.
      gvel=gv*vbar/1000.
c      if(cvel.lt.3.0d0) then
c        print *,'final',ci1,ci2,ci3
c        print *,cc,gv,om,vbar
c        print *,cvel,gvel
c      endif
c      print 106,cvel,gvel
  106 format('phase velocity',f9.5,'  group velocity',f9.5)
  101 format(6g15.5)
  105 format('Modal Eigenfunction, f=',f3.1,' Hz,  c=',f7.4,' km/s')
      return
      end
      subroutine grvel1(nl,cc,om,cci1,gvx,gvy)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      real*4 cc4,zz4,x4
      complex*16 pp,u0,ee,z1,z0,xnu,zzz,e1,e2,zla,xl,eigf
      complex*16 rt,tt,rt0,trc,co,cixz,dcixz,dc1,dc2
      complex*16 ci1,dci1,cixx,dcixx,cixy,dcixy,cizz,dcizz,ciyz,dciyz
      complex*16 u,xxnu,eye,pfac,uu,wukk,wuk
      common/stfff/w(3,101),t(3,3),ttl(3,3),s(3,3),stl(3,3),
     x                                           r(3,3),x(3),y(3)
      common/model/z(100),dz(100),rho(101),vp(101),vs(101),vp2(101),
     x               vp4(101),vs2(101),vss(101)
      common/model2/xmu(101),xla(101),xmu2(101),xla2(101),xla4(101)
      common/propag/xnu(6,101),xl(6,100),pfac(6,3),u(3,6)
      common/mstff/qq(6,6),wr(6),wi(6),zr(6,6),zi(6,6),iv(6),fv(6)
      common/pstff/pp(3),u0(3),ee(6,6,101),e1(6,6),e2(6,6),zla(6)
      common/rstff/rt(3,3,101),tt(3,3,101),rt0(3,3),trc(3,3)
      common/nstff/cc4(10000),zz4(10000),zzz(10000),ccc(10000),x4(10000)
      common/pstf2/co(6,101),eigf(1000,3)
      data pi/3.14159265358979d0/,eps/1.d-6/,tol/1.d-3/
      z1=dcmplx(1.d0,0.d0)
      z0=dcmplx(0.d0,0.d0)
      eye=dcmplx(0.d0,1.d0)
      ren=1.075190645d-3
      radi=6.371d6
      vbar=ren*radi
c      print *,nl,cc,om
      nlp=nl+1
c characteristics of the solution in the first layer
c            do k=1,6
c              print *,'for slowness:',xnu(k,1),', disp-stress vector'
c              print 103,(ee(j,k,1),j=1,6)
c            end do
  103 format(2g15.5)
c  calc the eigenfunction
c  play safe, since a 3x3 system of equations A.x=0 and detA=0
c  x is propto the normalized crossproduct of any two nonparallel rows of A
c  this condition gets dicey if we have degeneracy, because the 
c  determinant would have a double root, and all rows would be parallel
c   or when Rayleigh and Love decouple
c  so the safest route is to take the eigenvector of A with smallest eigenvalue
c  easiest variant: take Hermitian A^H.A, smallest (real) eigenvalue wins
      do i=1,3
        do j=1,3
          rt0(j,i)=z0
          do k=1,3
            rt0(j,i)=rt0(j,i)+conjg(trc(k,j))*trc(k,i)
          end do
          s(j,i)=dreal(rt0(j,i))
          stl(j,i)=dimag(rt0(j,i))
        end do
      end do
      call htridi(3,3,s,stl,wi,fv,fv,wr)
      do i=1,3
        do j=1,3
          t(j,i)=0.d0
        end do
        t(i,i)=1.d0
      end do
      call tql2(3,3,wi,fv,t,ierr)
      if(ierr.ne.0) then
        print *,'tql2 error',ierr
        stop
      endif
      call htribk(3,3,s,stl,wr,3,t,ttl) 
      do i=1,3
        co(i,1)=dcmplx(t(i,1),ttl(i,1))
        co(i+3,nlp)=z0
      end do
      do n=1,nl
c  upgoing coefs in the nth layer, downgoing coefs in the n+1 layer
        do i=1,3
          co(i+3,n)=z0
          co(i,n+1)=z0
          do j=1,3
            co(i+3,n)=co(i+3,n)+rt(i,j,n)*co(j,n)
            co(i,n+1)=co(i,n+1)+tt(i,j,n)*co(j,n)
          end do
        end do
      end do       
c  calculate the kinetic energy functional
c  rather than integrate numerically, 
c  we use the Lagrangian interaction matrices 
c  for sums of plane waves in each layer
c  cross terms are complex valued
      ci1=z0
      cixx=z0
      cixy=z0
      cixz=z0
      ciyz=z0
      cizz=z0
      do il=1,nlp
        if(il.eq.1) then
          h1=0.d0
        else
          h1=z(il-1)
        endif
        if(il.eq.nlp) then ! k=4,5,6 upgoing waves are zero
          nok=3
        else
          nok=6
        endif
        h2=z(il)
        a=xla(il)
        b=xla2(il)
        c=xla4(il)
        d=xmu(il)
        e=xmu2(il)
c  assemble the six polarization vectors and propagation factors
        do k=1,nok
C  START ALTERNATE COMPUTATION, HOPEFULLY MORE STABLE
c  normalize for the coefficient amplitudes in the Chen algorithm
c  during this computation, rather than dividing by extreme numbers later
c  upgoing coefs are referenced to layer bottom
c  downgoing coefs are referenced to layer top
c
c  we have a test to avoid underflows and overflows in the evenescent regime
c  if the decay factor in a layer is 10**(-50) or 10**(50) 
c  the eigenfunction in the layer cannot truly contribute to integrals
c  so we zero out the relevant pfac(.,1 or 2) pfac(.,3) propagation factors
c
c  pfac(k,3) has exponential negated in order to change  
c  division to multiplication in its use in the code below.  
c  This way we can set pfac(k,3)=0. when underflow is imminent
          test=dreal(xnu(k,il)*(h2-h1))
          iftest=1
          if(dabs(test).gt.100.d0) iftest=0
          if(k.le.3) then
            pfac(k,1)=z1
            if(il.lt.nlp.and.iftest.eq.1) then
              pfac(k,2)=zexp(xnu(k,il)*(h2-h1))
              pfac(k,3)=zexp(-xnu(k,il)*h1)
            else
              pfac(k,2)=z0
              pfac(k,3)=z0
            endif
          else
            if(iftest.eq.1) then
              pfac(k,1)=zexp(xnu(k,il)*(h1-h2))
              pfac(k,3)=zexp(-xnu(k,il)*h2)
            else
              pfac(k,1)=z0
              pfac(k,3)=z0
            endif
            if(il.lt.nlp) then
              pfac(k,2)=z1
            else
              pfac(k,2)=z0
            endif
          endif
C  END ALTERNATE COMPUTATION, HOPEFULLY MORE STABLE
          do i=1,3
            u(i,k)=co(k,il)*ee(i,k,il)
          end do
        end do
c  form matrix of vertical integrals -- check the coef convention!
c  kth wave is conjugated, kkth wave isnt
        do k=1,nok
c  note that we ADD conjg(xnu(k,il)) because xnu has a factor of i
          do kk=k,nok
            xxnu=xnu(kk,il)+conjg(xnu(k,il))
            if(zabs(xxnu).lt.eps) then
c  true if propagating wave with same vertical wavenumber (nu)
c  OR is integral of upgoing+downgoing evanescent waves with opposite nu
              xxnu=(h2-h1)/2.d0
c  normalize for the coefficient amplitudes in the Chen algorithm
c  upgoing coefs are referenced to layer bottom
c  downgoing coefs are referenced to layer top
c  taking a negative exponent in computation of pfac(k,3) allows us
c  to express this normalization as a multiplication
c  this way "zero=(0.,0.)" can be substituted in the far-evanescent region
              xxnu=xxnu*(pfac(kk,3)*conjg(pfac(k,3)))
            else
              xxnu=(pfac(kk,2)*conjg(pfac(k,2))
     x           -pfac(kk,1)*conjg(pfac(k,1)))/(2.d0*xxnu)
            endif
c  sum the plane wave contribution  
            uu=z0
            wukk=z0
            wuk=z0
            do i=1,3
              uu=uu+conjg(u(i,k))*u(i,kk)
              wuk=wuk+w(i,il)*conjg(u(i,k))
              wukk=wukk+w(i,il)*u(i,kk)
            end do
            wx=w(1,il)
            wy=w(2,il)
            wz=w(3,il)
      dci1=rho(il)*uu
      dcixx=(a-b+c-d+e)*conjg(u(1,k))*u(1,kk)
     x    +(d-e+2.d0*e*wx**2)*uu+(8.d0*c*wx**2+2.d0*e)*wuk*wukk
     x    +(b-4.d0*c-2.d0*e)*wx*(conjg(u(1,k))*wukk+wuk*u(1,kk))
      dcixy=(a-b+c-d+e)*(conjg(u(1,k))*u(2,kk)+conjg(u(2,k))*u(1,kk))
     x    +4.d0*e*wx*wy*uu+16.d0*c*wx*wy*wuk*wukk
     x    +(b-4.d0*c-2.d0*e)*wx*(conjg(u(2,k))*wukk+wuk*u(2,kk))
     x    +(b-4.d0*c-2.d0*e)*wy*(conjg(u(1,k))*wukk+wuk*u(1,kk))
      dc1=(a-b+c-2.d0*(d-e))*conjg(u(3,k))*u(1,kk)
     x    +2.d0*e*wx*wz*uu+8.d0*c*wx*wz*wuk*wukk
     x    +(b-4.d0*c-4.d0*e)*(wx*conjg(u(3,k))*wukk+wz*wuk*u(1,kk))
     x    +(d-e)*conjg(u(1,k))*u(3,kk)
     x    +2.d0*e*(wz*conjg(u(1,k))*wukk+wx*wuk*u(3,kk))
      dc2=(a-b+c-2.d0*(d-e))*conjg(u(1,k))*u(3,kk)
     x    +2.d0*e*wx*wz*uu+8.d0*c*wx*wz*wuk*wukk
     x    +(b-4.d0*c-4.d0*e)*(wz*conjg(u(1,k))*wukk+wx*wuk*u(3,kk))
     x    +(d-e)*conjg(u(3,k))*u(1,kk)
     x    +2.d0*e*(wx*conjg(u(3,k))*wukk+wz*wuk*u(1,kk))
c  note that we ADD conjg(xnu(k,il)) because xnu has a factor of i  (eye)
      dcixz=eye*(conjg(xnu(k,il))*dc1-xnu(kk,il)*dc2)
      dc1=(a-b+c-2.d0*(d-e))*conjg(u(3,k))*u(2,kk)
     x    +2.d0*e*wy*wz*uu+8.d0*c*wy*wz*wuk*wukk
     x    +(b-4.d0*c-4.d0*e)*(wy*conjg(u(3,k))*wukk+wz*wuk*u(2,kk))
     x    +(d-e)*conjg(u(2,k))*u(3,kk)
     x    +2.d0*e*(wz*conjg(u(2,k))*wukk+wy*wuk*u(3,kk))
      dc2=(a-b+c-2.d0*(d-e))*conjg(u(2,k))*u(3,kk)
     x    +2.d0*e*wy*wz*uu+8.d0*c*wy*wz*wuk*wukk
     x    +(b-4.d0*c-4.d0*e)*(wz*conjg(u(2,k))*wukk+wy*wuk*u(3,kk))
     x    +(d-e)*conjg(u(3,k))*u(2,kk)
     x    +2.d0*e*(wy*conjg(u(3,k))*wukk+wz*wuk*u(2,kk))
c  note that we ADD conjg(xnu(k,il)) because xnu has a factor of i  (eye)
      dciyz=eye*(conjg(xnu(k,il))*dc1-xnu(kk,il)*dc2)
c  dont really need cizz, just like we dont need ciyy
      dcizz=(a-b+c-d+e)*conjg(u(3,k))*u(3,kk)
     x    +(d-e+2.d0*e*wz**2)*uu+(8.d0*c*wz**2+2.d0*e)*wuk*wukk
     x    +(b-4.d0*c-2.d0*e)*wz*(conjg(u(3,k))*wukk+wuk*u(3,kk))
c  the extra factor of i in xnu cancels in this product
      dcizz=dci4*xnu(kk,il)*conjg(xnu(k,il))
            if(k.eq.kk) then
              ci1=ci1+dci1*xxnu
              cixx=cixx+dcixx*xxnu
              cixz=cixz+dcixz*xxnu
              cixy=cixy+dcixy*xxnu
              ciyz=ciyz+dciyz*xxnu
              cizz=cizz+dcizz*xxnu
            else
              ci1=ci1+2.d0*dreal(xxnu*dci1)
              cixx=cixx+2.d0*dreal(xxnu*dcixx)
              cixz=cixz+2.d0*dreal(xxnu*dcixz)
              cixy=cixy+2.d0*dreal(xxnu*dcixy)
              ciyz=ciyz+2.d0*dreal(xxnu*dciyz)
              cizz=cizz+2.d0*dreal(xxnu*dcizz)
            endif
          end do
        end do
      end do
      cci1=dreal(ci1)
      gvx=dreal((cixx/cc+cixz/(2.d0*om))/ci1)
      gvy=dreal((cixy/cc+ciyz/om)/(2.d0*ci1))
  101 format(6g15.5)
      return
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
      subroutine solve(nn,a,x,y)
c  solves the nxn system of equations a*x=y using gaussian elimination 
c  and partial pivoting
c  if n<0 the lu decomposition is already done
c  note that the matrix a is modified
      implicit real*8 (a-h,o-z)                                                 
      implicit integer*4 (i-n)                                                  
      common/solve_/ip(1000)
      dimension a(1),x(1),y(1)
      n=nn
      if(n.gt.0)then
        call lup(n,a,ip)
      else
        n=-n
      endif
c      print 101,((i,j,a(i,j),j=1,n),i=1,n)
c  101 format(' a(',2i2,')=',e15.5)
c      print 102,(ip(i),i=1,n)
c  102 format(10i5)
      call bcktr(n,a,x,y,ip)
      return
      end
      subroutine lup(n,a,ip)
c  finds lu decomp of a using partial pivoting         
c  output in c (upper triangle/unit lower triangle) and                      
c  pivoting sequence returned in ip(n)                                          
      implicit real*8 (a-h,o-z)                                                 
      implicit integer*4 (i-n)                                                  
      dimension a(n,1),ip(1)                             
      tol=1.d-17                                                                
      do 50 i=1,n
   50 ip(i)=i                                                                   
      nm1=n-1                                                                   
      if(n.eq.1) go to 700                                                      
      do 100 i=1,nm1                                                            
      aam=0.d0                                                                  
      do 200 j=i,n                                                              
      aai=a(ip(j),i)**2                                          
      if(aam.gt.aai) go to 200                                                  
      aam=aai                                                                    
      ipi=ip(j)                                                                 
      jm=j                                                                      
  200 continue                                                                  
      if(aam.lt.tol) go to 400                                                  
      ip(jm)=ip(i)                                                              
      ip(i)=ipi                                                                 
      i1=i+1                                                                    
      do 100 j=i1,n                                                             
      ipj=ip(j)                                                                 
c  if victim index is already zero, dont bother to rub it out                   
      tem=dabs(a(ipj,i))
      if(tem.lt.tol) go to 100                                                  
      b=(a(ipj,i)*a(ipi,i))/aam                             
      a(ipj,i)=b                                                                
      do 500 k=i1,n                                                             
      a(ipj,k)=a(ipj,k)-b*a(ipi,k)
  500 continue
  100 continue                                                                  
  700 continue              
      return                                                                    
  400 print 101,aam,i                                                           
  101 format('near-zero pivot ',e12.5,'  on column',i3)                         
      stop
      end                                                                       
      subroutine bcktr(n1,z,dr,er,ip)
c  performs backtransform on input vector er - 'y'
c  to find solution dr - 'x'
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension z(n1,1),ip(1),er(1),dr(1)
c  back transform with unit lower triangular matrix
      do 300 i=1,n1
  300 dr(i)=er(ip(i))
      if(n1.eq.1) go to 400
      do 310 i=2,n1
      i1=i-1
      do 310 j=1,i1
  310 dr(i)=dr(i)-z(ip(i),j)*dr(j)
  400 continue
c  back transform with upper triangular matrix
      do 320 ii=1,n1
      i=n1+1-ii
      ip1=i+1
      if(i.eq.n1) go to 320
      do 330 j=ip1,n1
  330 dr(i)=dr(i)-z(ip(i),j)*dr(j)
  320 dr(i)=dr(i)/z(ip(i),i)
      return
      end
      subroutine csolve(nn,a,ai,x,xi,y,yi)
c  solves the complex nxn system of equations a*x=y using gaussian elimination 
c  and partial pivoting
      implicit real*8 (a-h,o-z)                                                 
      implicit integer*4 (i-n)                                                  
      dimension ip(1000),a(1),x(1),y(1),ai(1),xi(1),yi(1)
      n=nn
      if(n.gt.0) then 
        call clup(n,a,ai,ip)
      else
        n=-n
      endif
      call cbcktr(n,a,ai,x,xi,y,yi,ip)
      return
      end
      subroutine clup(n,a,ai,ip)                                    
c  finds lu decomp of a+i*ai using partial pivoting 
c  pivoting sequence returned in ip(n)
      implicit real*8 (a-h,o-z)                                                 
      implicit integer*4 (i-n)                                                  
      dimension a(n,1),ai(n,1),ip(1)                             
      tol=1.d-17                                                                
c  initialize permutation vector                                                
      do 50 i=1,n
   50 ip(i)=i                                                                   
      nm1=n-1                                                                   
      if(n.eq.1) go to 700                                                      
      do 100 i=1,nm1                                                            
      aam=0.d0                                                                  
      do 200 j=i,n                                                              
      aai=a(ip(j),i)**2+ai(ip(j),i)**2                                          
      if(aam.gt.aai) go to 200                                                  
      aam=aai                                                                   
      ipi=ip(j)                                                                 
      jm=j                                                                      
  200 continue                                                                  
      if(aam.lt.tol) go to 400                                                  
      ip(jm)=ip(i)                                                              
      ip(i)=ipi                                                                 
      i1=i+1                                                                    
      do 100 j=i1,n                                                             
      ipj=ip(j)                                                                 
c  if victim index is already zero, dont bother to rub it out                   
      tem=dabs(a(ipj,i))+dabs(ai(ipj,i))                                        
      if(tem.lt.tol) go to 100                                                  
      b=(a(ipj,i)*a(ipi,i)+ai(ipj,i)*ai(ipi,i))/aam                             
      bi=(ai(ipj,i)*a(ipi,i)-a(ipj,i)*ai(ipi,i))/aam                            
      a(ipj,i)=b                                                                
      ai(ipj,i)=bi                                                              
      do 500 k=i1,n                                                             
      a(ipj,k)=a(ipj,k)-b*a(ipi,k)+bi*ai(ipi,k)                                 
  500 ai(ipj,k)=ai(ipj,k)-b*ai(ipi,k)-bi*a(ipi,k)                               
  100 continue                                                                  
  700 continue                                                                  
      return                                                                    
  400 print 101,aam,i                                                           
  101 format('near-zero pivot ',e12.5,'  on column',i3)                         
      stop                                                                 
      end                                                                       
      subroutine cbcktr(n1,z,zi,dr,di,er,ei,ip)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension z(n1,1),ip(1),er(1),ei(1),dr(1),di(1),zi(n1,1)
c  back transform with unit lower triangular matrix
      do 300 i=1,n1
      dr(i)=er(ip(i))
  300 di(i)=ei(ip(i))
      if(n1.eq.1) go to 400
      do 310 i=2,n1
      i1=i-1
      do 310 j=1,i1
      zkk1=z(ip(i),j)
      zkk2=zi(ip(i),j)
      dr(i)=dr(i)-zkk1*dr(j)+zkk2*di(j)
  310 di(i)=di(i)-zkk1*di(j)-zkk2*dr(j)
  400 continue
c  back transform with upper triangular matrix
      do 320 ii=1,n1
      i=n1+1-ii
      ip1=i+1
      uii=z(ip(i),i)**2+zi(ip(i),i)**2
      if(i.eq.n1) go to 340
      do 330 j=ip1,n1
      zkk1=z(ip(i),j)
      zkk2=zi(ip(i),j)
      dr(i)=dr(i)-zkk1*dr(j)+zkk2*di(j)
  330 di(i)=di(i)-zkk1*di(j)-zkk2*dr(j)
  340 dri=dr(i)
      dii=di(i)
      zkk1=z(ip(i),i)
      zkk2=zi(ip(i),i)
      dr(i)=(dri*zkk1+dii*zkk2)/uii
  320 di(i)=(dii*zkk1-dri*zkk2)/uii
      return
      end
      subroutine get_eigen(akx,f,nl,h,usx,usz,ur)
c  program will get the value of the eigenfct at the surface and source depths
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      complex*16 pp,u0,ee,z1,z0,xnu,zzz,e1,e2,zla,xl,co,eigf
      complex*16  rt,tt,rt0,trc,u,eye,pfac,usx,usz,ur
      common/stfff/w(3,101),t(3,3),ttl(3,3),s(3,3),stl(3,3),
     x                                           r(3,3),x(3),y(3)
      common/model/z(100),dz(100),rho(101),vp(101),vs(101),vp2(101),
     x               vp4(101),vs2(101),vss(101)
      common/model2/xmu(101),xla(101),xmu2(101),xla2(101),xla4(101)
      common/propag/xnu(6,101),xl(6,100),pfac(6,3),u(3,6)
      common/mstff/qq(6,6),wr(6),wi(6),zr(6,6),zi(6,6),iv(6),fv(6)
      common/mstff2/fv2(3),fv3(3)
      common/pstff/pp(3),u0(3),ee(6,6,101),e1(6,6),e2(6,6),zla(6)
      common/rstff/rt(3,3,101),tt(3,3,101),rt0(3,3),trc(3,3)
      common/pstf2/co(6,101),eigf(1000,3)
      common/defect/idfct(4,101),adf(2,101)
      dimension usx(1),usz(1),ur(1)
      data pi/3.14159265358979d0/,eps/1.d-8/,tol/5.d-4/
c  we reduce the condition numbers of matrices by 
c  normalizing physical quantities to make them dimensionless
c  Ill use the normal-mode normalizations
c  which are a little peculiar for the crust, but what the hey!
      rbar=5.515d3
      ren=1.075190645d-3
      radi=6.371d6
      vbar=radi*ren
      con=rbar*radi*radi*ren*ren
      z1=dcmplx(1.d0,0.d0)
      z0=dcmplx(0.d0,0.d0)
      eye=dcmplx(0.d0,1.d0)
      ntry=1000
      nlp=nl+1
c  calc the eigenfunction
c  play safe, since a 3x3 system of equations A.x=0 and detA=0
c  x is propto the normalized crossproduct of any two nonparallel rows of A
c  this condition gets dicey if we have degeneracy, because the 
c  determinant would have a double root, and all rows would be parallel
c   or when Rayleigh and Love decouple
c  so the safest route is to take the eigenvector of A with smallest eigenvalue
      do ii=1,3
        do jj=1,3
          s(jj,ii)=dreal(trc(jj,ii))
          stl(jj,ii)=dimag(trc(jj,ii))
        end do
      end do      
c  into eispack - fv,fv2,fv3 are dummy arrays
      call cbal(3,3,s,stl,low,igh,fv)
      call corth(3,3,low,igh,s,stl,fv2,fv3)
      call comqr2(3,3,low,igh,fv2,fv3,s,stl,wr,wi,t,ttl,ierr)
      if(ierr.ne.0) then
        print *,'comqr2 error',ierr
        stop
      endif
      call cbabk2(3,3,low,igh,fv,3,t,ttl)
c      print *,wr(1),wi(1)
      amn=wr(1)**2+wi(1)**2
      ij=1
      do ii=2,3
c        print *,wr(ii),wi(ii)
        amm=wr(ii)**2+wi(ii)**2
        if(amm.lt.amn) then
          ij=ii
          amn=amm
        endif
      end do
      sum=0.d0
      do ii=1,3
        u0(ii)=dcmplx(t(ii,ij),ttl(ii,ij))
        sum=sum+zabs(u0(ii))**2
      end do
      sum=dsqrt(sum)
      do i=1,3
        co(i,1)=u0(i)/sum
        co(i+3,nlp)=z0
      end do
c      print 101,(co(i,1),i=1,3)
      do n=1,nl
c  upgoing coefs in the nth layer, downgoing coefs in the n+1 layer
        do i=1,3
          co(i+3,n)=z0
          co(i,n+1)=z0
          do j=1,3
            co(i+3,n)=co(i+3,n)+rt(i,j,n)*co(j,n)
            co(i,n+1)=co(i,n+1)+tt(i,j,n)*co(j,n)
          end do
        end do
c        print 101,(co(i+3,n),i=1,3)
c        print 101,(co(i,n+1),i=1,3)
      end do       
c  first z=0, then z=h
      il=1
      h1=0.d0
      h2=z(1)
      do i=1,3
        ur(i)=z0
        do k=1,3
          ur(i)=ur(i)
     x        +co(k,il)*ee(i,k,il)*(zexp(xnu(k,il)*(-h1)))
     x      +co(k+3,il)*ee(i,k+3,il)*(zexp(xnu(k+3,il)*(-h2)))
        end do
c  check for the xtra terms associated with defective matrices
        if(idfct(1,il).ne.0) then
          ii=idfct(1,il)
          jj=idfct(2,il)
          ur(i)=ur(i)
     x +co(jj,il)*adf(1,il)*ee(i,ii,il)*(-h1)*(zexp(xnu(ii,il)*(-h1)))
        endif
        if(idfct(3,il).ne.0) then
          ii=idfct(3,il)
          jj=idfct(4,il)
          ur(i)=ur(i)
     x +co(jj,il)*adf(2,il)*ee(i,ii,il)*(-h2)*(zexp(xnu(ii,il)*(-h2)))
        endif
      end do
c  then z=h, ---> find the layer we are in
c  we calculate the vector gradient of eigenfct
      do n=1,nl
        if(h.gt.z(n)) il=il+1
      end do
      if(il.eq.1) then
        h1=0.d0
      else
        h1=z(il-1)
      endif
      h2=z(il)
c  zzz is horizontal wavenumber times i
c  the minus sign will arrive when we conjugate usz,usx in main program
      zzz=eye*akx
      do i=1,3
        usx(i)=z0
        usz(i)=z0
        do k=1,3
          usx(i)=usx(i)
     x      +zzz*co(k,il)*ee(i,k,il)*(zexp(xnu(k,il)*(h-h1)))
     x      +zzz*co(k+3,il)*ee(i,k+3,il)*(zexp(xnu(k+3,il)*(h-h2)))
          usz(i)=usz(i)
     x       +xnu(k,il)*co(k,il)*ee(i,k,il)*(zexp(xnu(k,il)*(h-h1)))
     x +xnu(k+3,il)*co(k+3,il)*ee(i,k+3,il)*(zexp(xnu(k+3,il)*(h-h2)))
        end do
c  check for the xtra terms associated with defective matrices
        if(idfct(1,il).ne.0) then
          ii=idfct(1,il)
          jj=idfct(2,il)
          usx(i)=usx(i)+zzz*co(jj,il)*adf(1,il)*ee(i,ii,il)
     x 				*(h-h1)*(zexp(xnu(ii,il)*(h-h1)))
          usz(i)=usz(i)+xnu(ii,il)*co(jj,il)*adf(1,il)*ee(i,ii,il)
     x 				*(h-h1)*(zexp(xnu(ii,il)*(h-h1)))
     x 	+co(jj,il)*adf(1,il)*ee(i,ii,il)*(zexp(xnu(ii,il)*(h-h1)))
        endif
        if(idfct(3,il).ne.0) then
          ii=idfct(3,il)
          jj=idfct(4,il)
          usx(i)=usx(i)+zzz*co(jj,il)*adf(2,il)*ee(i,ii,il)
     x 				*(h-h2)*(zexp(xnu(ii,il)*(h-h2)))
          usz(i)=usz(i)+xnu(ii,il)*co(jj,il)*adf(2,il)*ee(i,ii,il)
     x 				*(h-h2)*(zexp(xnu(ii,il)*(h-h2)))
     x 	+co(jj,il)*adf(2,il)*ee(i,ii,il)*(zexp(xnu(ii,il)*(h-h2)))
        endif
      end do
c      print *,zzz
c      print *,(ur(i),i=1,3)
c      print *,(usz(i),i=1,3)
c      print *,(usx(i),i=1,3)
c      pause
  106 format('phase velocity',f18.14,'  group velocity',f8.5)
  101 format(6g15.5)
  105 format('Modal Eigenfunct, f=',f5.3,' Hz,  c=',f7.4,' km/s')
      return
      end
      subroutine splneq(nn, u, s)
cz$$$$ calls no other routines
c  finds coeffs for a spline interpolation of equally spaced data
c  based on 'spline interpolation ona  digital computer' by r.f.thompson
c  nn  number of data points (may be negative - see d1,d2 below)
c  u  array of function values to be interpolated, assumed to samples at equal i
c     intervals of a twice differentiable function.
c  s  array to hold computed values of 2nd derivative of spline fit at sample
c     points.  these values are required by evaleq  to interpolate
c  if the user wishes to force specific values on the derivatives at the end
c  points, he should put h*du(1)/dx  nad  h*du(n)/dx  in  s(1),s(2), then call
c  splneq  with nn=-number of terms in series. h = sample spacing in x.
c  normally the derivatives are found by fitting a parabola through the
c  1st and last 3 points.
c  if the number of terms is between 1 and 3, straight-line interpolation is don
      implicit integer*4 (i-n)
      dimension u(3),s(3),a(13)
c
      n=iabs(nn)
      if (n.le.3) go to 5000
      d1=-0.5*u(3)  +2.0*u(2)  -1.5*u(1)
      dn= 0.5*u(n-2)-2.0*u(n-1)+1.5*u(n)
      if (nn.gt.0) go to 1000
      d1=s(1)
      dn=s(2) 
 1000 a(1)=2.0
      a(2)=3.5
      s(1)=u(2)-u(1)-d1
      s(2)=u(1)-2.0*u(2)+u(3)-0.5*s(1)
      n1=n-1
      do 3000 i=3,n1
      if (i.gt.13) go to 3000
      k=i
      a(k)=4.-1.0/a(k-1)
 3000 s(i)=u(i-1)-2.*u(i)+u(i+1)-s(i-1)/a(k-1)
      s(n)=u(n1)-u(n)+dn-s(n1)/a(k)
      s(n)=6.0*s(n)/(2.0-1.0/a(k))
      n2=n-2
c  compute 2nd derivatives by back-substitution
c  the array  a  tends to a constant (2+sqrt(3)) so only 13 elements are needed
      do 4000 j=1,n2
      i=n-j
      k=min0(i,k)
 4000 s(i)=(6.0*s(i)-s(i+1))/a(k)
      s(1)=3.0*s(1)-0.5*s(2)
      return
c  series too short for cubic spline.  fit straight lines.
 5000 do 5500 i=1,n
 5500 s(i)=0.0
      return
      end                                                               splneq
      function evaleq(y, nn, u, s,ick,h)
cz$$$$ calls no other routines
c  performs spline interpolation of equally spaced data.
c  based on 'spline interpolation on a digital computer'( by r.f.thompson.
c  evaluates a spline interpolate in a set of equally spaced samples.
c  the routine  splneq  should be called first, to establish the array  s .
c  y  the  coordinate at which interpolate is required, with y=1 for 1st
c     sample point, y=2 for the 2nd, etc.  if actual spacing is  h  and  x1 is
c     the 1st sample coordinate use  y = 1.0 + (x-x1)/h
c  nn  number of samples of function in original set.
c  u  array containing function samples.
c  s  array of normalized 2nd derivatives, computed by  splneq.  the derivatives
c     have been multiplied by h**2, where h is the sample spacing.
c  if  y  is out of the range (1,nn), the 1st or last sample value is used.
      dimension u(1),s(1)
      data z3/.333333333333/,z6/.16666666666666667/
c
      if (y.le.1.0) go to 1500
      if (y.ge.float(nn)) go to 2000
      k1=y
      k=k1+1
      dk=k-y
      dk1=y-k1
      if(ick.eq.1) go to 1000
      ff1=s(k1)*dk*dk*dk
      ff2=s(k)*dk1*dk1*dk1
      evaleq=(dk*(6.0*u(k1)-s(k1))+ dk1*(u(k)*6.0-s(k)) + ff1 +ff2)/6.0
      return
c  evaluate the first derivative
 1000 a1=dk1*(1.-dk1/2.)-z3
      a2=dk1*dk1/2.-z6
      a3=u(k)-u(k1)+a1*s(k1)+s(k)*a2
      evaleq=a3/h
      return
c  out of range.  supply constant values
 1500 evaleq=u(1)
      return
 2000 evaleq=u(nn)
      return
      end                                                               evaleq
c      subroutine sacread_noswap(name,a,ierr)
      subroutine sacread(name,a,ierr)
c  read a SAC-format file, dont swap the bytes
      real*4 ahead,ahd,ahhd
      character*80 name
      character*4 chead(158),chd,chhd
      character*1 chhead(4,158)
      common/header/ahead(158)
      dimension a(1)
      dimension iah(158)
      equivalence (ahead,iah),(ahead,chead),(ahead,chhead(1,1))
      equivalence (chd,ahd),(chhd,ahhd)
      ierr=0
      open(8,file=name,access='direct',recl=512,err=999)
c  read the 158-word header, can use its info 
      read(8,rec=1,err=998)(ahead(i),i=1,128)
      print *,(ahead(i),i=50,55)
      print *,(iah(70+i),i=1,10)
      nscan=iah(80)
      irec=2
      read(8,rec=2,err=998)(ahead(i+128),i=1,30),(a(j),j=1,98)
      if(chead(145).ne.' GRN'.and.chead(145).ne.' grn') then
        nword=nscan+158
        nrec=(nword-1)/128+1
        last=nword-128*(nrec-1)
      else       
        print *,'reading greens functions' 
        nword=6*nscan+158
        nrec=(nword-1)/128+1
        last=nword-128*(nrec-1)
      endif
      if(nrec.gt.3) then
        do irec=3,nrec-1
          read(8,rec=irec,err=998)(a((irec-3)*128+98+i),i=1,128)
        end do
      endif
c      print *,'rec=',nrec,'  Im trying to read the last partial record'
      read(8,rec=nrec,err=998)(a((nrec-3)*128+98+i),i=1,last)
      close(8)
      return
  999 print *,'open error'
      ierr=1
      return
  998 print *,'read error: reading',irec,' out of',nrec
      ierr=1
      print *,irec
      close(8)
      return
      end
      subroutine sacread_swap(name,a,ierr)
c      subroutine sacread(name,a,ierr)
c  read a SAC-format file -- swap bytes for Intel mac for all numbers, not characters
      real*4 ahead,ahd,ahhd
      character*80 name
      character*4 chead(158),chd,chhd
      character*1 chhead(4,158)
      common/header/ahead(158)
      dimension a(1)
      dimension iah(158)
      equivalence (ahead,iah),(ahead,chead),(ahead,chhead(1,1))
      equivalence (chd,ahd),(chhd,ahhd)
      ierr=0
      open(8,file=name,access='direct',recl=512,err=999)
c  read the 158-word header, can use its info 
      read(8,rec=1,err=998)(ahead(i),i=1,128)
c  swap bytes
      do i=1,120
        chd=chead(i)
	do j=1,4
	  chhead(j,i)=chd(5-j:5-j)
	end do
      end do
      print *,(ahead(i),i=50,55)
      print *,(iah(70+i),i=1,10)
      nscan=iah(80)
      irec=2
      read(8,rec=2,err=998)(ahead(i+128),i=1,30),(a(j),j=1,98)
c  swap bytes
      do i=1,98
        ahd=a(i)
	do j=1,4
	  chhd(j:j)=chd(5-j:5-j)
	end do
	a(i)=ahhd
      end do
c      print *,'rec=2'
      if(chead(145).ne.' GRN'.and.chead(145).ne.' grn') then
        nword=nscan+158
        nrec=(nword-1)/128+1
        last=nword-128*(nrec-1)
      else       
        print *,'reading greens functions' 
        nword=6*nscan+158
        nrec=(nword-1)/128+1
        last=nword-128*(nrec-1)
      endif
      if(nrec.gt.3) then
        do irec=3,nrec-1
c          print *,'rec=',irec
          read(8,rec=irec,err=998)(a((irec-3)*128+98+i),i=1,128)
c  swap bytes
          do i=1,128
            ahd=a((irec-3)*128+98+i)
       	    do j=1,4
	      chhd(j:j)=chd(5-j:5-j)
	    end do
	    a((irec-3)*128+98+i)=ahhd
          end do
        end do
      endif
      print *,'byteswapping read'
c  this here is a big fat kluge --- I cant read the last partial record
c  so I skip the last fragment and adjust the header to discard the points.
c      print *,'NOT READING ',last,' points in last record -- SACBUGFIX'
      read(8,rec=nrec,err=998)(a((nrec-3)*128+98+i),i=1,last)
c  swap bytes
      do i=1,last
        ahd=a((nrec-3)*128+98+i)
       	do j=1,4
	  chhd(j:j)=chd(5-j:5-j)
	end do
	a((nrec-3)*128+98+i)=ahhd
      end do
c      nscan=nscan-last
c      iah(80)=nscan
      close(8)
      return
  999 print *,'open error'
      ierr=1
      return
  998 print *,'read error: reading',irec,' out of',nrec
      ierr=1
      print *,irec
      close(8)
      return
      end
c      subroutine sacout_noswap(name,a)
      subroutine sacout(name,a)
c  write a SAC-format file
      real*4 ahead
      character*80 name
      character*4 chead(158)
      common/header/ahead(158)
      dimension a(1)
      dimension iah(158)
      equivalence (ahead,iah),(ahead,chead)
      open(8,file=name,access='direct',recl=512)
c  write the 158-word header
      write(8,rec=1)(ahead(i),i=1,128)
      print *,(ahead(i),i=50,55)
      print *,(iah(70+i),i=1,10)
      nscan=iah(80)
      if(chead(145).ne.' GRN'.and.chead(145).ne.' grn') then
        nword=nscan+158
        nrec=(nword-1)/128+1
      else       
        print *,'writing greens functions' 
        nword=6*nscan+158
        nrec=(nword-1)/128+1
      endif
      write(8,rec=2)(ahead(i+128),i=1,30),(a(j),j=1,98)
      do irec=3,nrec
        write(8,rec=irec)(a((irec-3)*128+98+i),i=1,128)
      end do
      print *,nrec,' records written'
      close(8)
      return
      end
c      subroutine sacout(name,a)
      subroutine sacout_swap(name,a)
c  write a SAC-format file
      real*4 ahead,ahd,ahhd
      character*80 name
      character*4 chead(158),chd,chhd
      character*1 chhead(4,158)
      common/header/ahead(158)
      dimension a(1)
      dimension iah(158),output(158)
      equivalence (ahead,iah),(ahead,chead),(ahead,chhead(1,1))
      equivalence (chd,ahd),(chhd,ahhd)
      open(8,file=name,access='direct',recl=512)
      print *,(ahead(i),i=50,55)
      print *,(iah(70+i),i=1,10)
c  write the 158-word header
c  swap bytes
      do i=1,120
        chd=chead(i)
	do j=1,4
	  chhd(j:j)=chd(5-j:5-j)
	end do
	output(i)=ahhd
      end do
      do i=121,128
        output(i)=ahead(i)
      end do
      write(8,rec=1)(output(i),i=1,128)
      nscan=iah(80)
      if(chead(145).ne.' GRN'.and.chead(145).ne.' grn') then
        nword=nscan+158
        nrec=(nword-1)/128+1
      else       
        print *,'writing greens functions' 
        nword=6*nscan+158
        nrec=(nword-1)/128+1
      endif
c  swap bytes
      do i=1,30
        output(i)=ahead(128+i)
      end do
      do i=1,98
        ahd=a(i)
	do j=1,4
	  chhd(j:j)=chd(5-j:5-j)
	end do
	output(30+i)=ahhd
      end do
      write(8,rec=2)(output(i),i=1,128)
      do irec=3,nrec
        iii=(irec-3)*128+98
        do i=1,128
          ahd=a(iii+i)
	  do j=1,4
	    chhd(j:j)=chd(5-j:5-j)
  	  end do
	  output(i)=ahhd
        end do
        write(8,rec=irec)(output(i),i=1,128)
      end do
      print *,nrec,'byteswapped records written'
      close(8)
      return
      end

      subroutine fft2(ar,ai,n)
c  fft routine with 2 real input arrays rather than complex
c  fft2 subroutine mults by exp(i\omega t)
c  OUTPUT: f=0 in ar(1), f=f_N in ar(n/2+1)
c
c  n is a power of two. 
c  For a real-valued time series, the FT of negative frequencies resides in the
c  output array beyond the Nyquist f_N with the usual conjugation:
c
c          ar(i)=ar(n+2-i), ai(i)=-ai(n+2-i)
c
c  fft2 is NOT a unitary transform, mults the series by sqrt(n)
c  the inverse FT can be effected by running fft2 on the conjugate of
c  the FFT-expansion, then taking the the conjugate of the output, 
c  and dividing thru by N. to wit:
c
c   assume Xr, Xi is in freq domain, xr, xi in the time domain
c
c   (Xr,Xi)=fft2(xr,xi,N)
c   (xr,-xi)=fft2(Xr,-Xi,N)/N
c
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension ar(1),ai(1)
      mex=dlog(dble(float(n)))/.693147d0
      nv2=n/2
      nm1=n-1
      j=1
      do 7 i=1,nm1
      if(i .ge. j) go to 5
      tr=ar(j)
      ti=ai(j)
      ar(j)=ar(i)
      ai(j)=ai(i)
      ar(i)=tr
      ai(i)=ti
   5  k=nv2
   6  if(k .ge. j) go to 7
      j=j-k
      k=k/2
      go to 6
   7  j=j+k
      pi=3.14159265358979d0
      do 20 l=1,mex
      le=2**l
      le1=le/2
      ur=1.0
      ui=0.
      wr=dcos(pi/le1 )
      wi=dsin (pi/le1)
      do 20 j=1,le1
      do 10 i=j,n,le
      ip=i+le1
      tr=ar(ip)*ur - ai(ip)*ui
      ti=ai(ip)*ur + ar(ip)*ui
      ar(ip)=ar(i)-tr
      ai(ip)=ai(i) - ti
      ar(i)=ar(i)+tr
      ai(i)=ai(i)+ti
  10  continue
      utemp=ur
      ur=ur*wr - ui*wi
      ui=ui*wr + utemp*wi
  20  continue
      return
      end
