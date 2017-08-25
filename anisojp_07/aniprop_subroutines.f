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
c  bugfix for defective_aniprop matrices: 8/22/95
c  bugfix for transition to evanescence in top layer: 9/10/95
c  
c
c  compile sequence 
c  eislib is the Eispack library
c  f77 -o aniprop -fast -native -O5 aniprop.f /data/d4/park/Ritz/eislib 
c
c   compile sequence with Park-specific plotting programs for debugging
c  plotlib&jlib link with a complaint-that-I-ignore in Solaris 4.x
c  xf77 -o aniprop -O aniprop.f /data/d4/park/Plotxy/plotlib /data/d4/park/Ritz/eislib /data/d4/park/Ritz/jlib
c
c  for hexagonally symmetric media
c  reads fast axis orientation, constants A,B,C,D,E from file animodel
c  calculate quadratic eigenvalue problem based on the Christoffel matrix
c  see appendix of P. Shearer's thesis
      
      subroutine aniprop_subroutines(z_in,vp_in,vp2_in,vp4_in,
     x           vs_in,vs2_in,rho_in,theta_in,phig_in,nl,baz,
     x           rphase,rgroup,lphase,lgroup,period,nfrq)

      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
     
      real*8 baz
      integer nl,nfrq 
      real*8 theta_in(nl+1),phig_in(nl+1),z_in(nl+1),vp_in(nl+1)
      real*8 vp2_in(nl+1),vp4_in(nl+1),vs_in(nl+1),vs2_in(nl+1)
      real*8 rho_in(nl+1)
      real*8 rphase(nfrq),rgroup(nfrq)
      real*8 lphase(nfrq),lgroup(nfrq),period(nfrq) 

c      real*4 cc4,zz4,x4
      character*80 name,title
      complex*16 pp,u0,ee,z1,z0,zz,xnu,zzz,e1,e2,zla,xl,c,eigf
      complex*16  rt,tt,rt0,trc,pfac,u
      real*4 cc4,zz4,x4,frqq
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
      common/disper/cvel(8192,202),gvel(8192,202),nc(202),frqq(8192)
      common/disper2/roota(202),rootb(202),jtrval(202),kroots(8192)
      common/evanes/ievan(10000)
      data pi/3.14159265358979d0/,eps/1.d-6/,tol/1.d-3/
c  we reduce the condition numbers of matrices by 
c  normalizing physical quantities to make them dimensionless
c  Ill use the normal-mode normalizations
c  which are a little peculiar for the crust, but what the hey!


      rbar=5.515d3
      ren=1.075190645d-3
      radi=6.371d6
      vbar=ren*radi
      con=rbar*vbar**2
c initialize cvel array to zero
      do j=1,maxbr
        do i=1,8192
          cvel(i,j)=0.d0
          gvel(i,j)=0.d0
        end do
      end do  
      z1=dcmplx(1.d0,0.d0)
      z0=dcmplx(0.d0,0.d0)
      tolint=100.d0/vbar
c      nfrq=1000
      frqmax=0.1d0
c      nfrq=2500
c      frqmax=2.5d0
      maxbr=202
      ntry=1000
      mtry=5
      nfstart=4
 102  format(a)

c  read in theta,phi in degrees - polar coords of fast axis
      nlp=nl+1
      nlm=nl-1
      do i=1,nlp

        theta=theta_in(i)
        phig=phig_in(i)

        phi=phig-baz

        w(1,i)=dsind(theta)*dcosd(phi)
        w(2,i)=dsind(theta)*dsind(phi)
        w(3,i)=dcosd(theta)

c        print *,(w(j,i),j=1,3)
c  read depth to ith interface, vp (m/sec), pk-to-pk cos(2th) relative P pert
c  pk-to-pk cos(4th) relative P pert, v_s,  pk-to-pk cos(2th) relative S pert
c  density (kg/m**3)
c        read(7,*) z(i),vp(i),vp2(i),vp4(i),vs(i),vs2(i),rho(i)
        z(i)  = z_in(i)
        vp(i) = vp_in(i)
        vp2(i)= vp2_in(i)
        vp4(i)= vp4_in(i)
        vs(i) = vs_in(i)
        vs2(i)= vs2_in(i)
        rho(i)= rho_in(i)

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
  104 format(2f10.1)
c  search for cmin, cmax is halfspace velocity
      cmin=vs(1)
      vss(1)=vs(1)
      do i=2,nlp
        if(cmin.gt.vs(i)) cmin=vs(i)
        vss(i)=vs(i)
      end do
      cmax=vs(nlp)
c      print *,'ntry',ntry
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
c      print *,cmin,cmax,' =cmin,cmax'
c      print *,(vss(j),j=1,nlp+1)
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
      ddc=(cmax-cmin-eps)/(ntry-1)
      do sl=cmin+eps,cmax,ddc
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
c  first make an array of values to look for zeroes in
c  note that we arent trying to capture the fundamental Rayleigh wave
c   which usually has c<vsmin
        do jj=1,ntry
          cc=ccc(jj)
          call zzget(nl,cc,om,zz,iev)
          zzz(jj)=zz
        end do
        do i=1,ntry
          zz4(i)=real(zzz(i))
        end do
c        call plotit(cc4,zz4,dum,ntry,'Propagating Mode Condition',
c     x  'Phase velocity (km/sec)','Surface Traction',2,0,0,0,0)
c        do i=1,ntry
c          zz4(i)=imag(zzz(i))
c        end do
c        call plotit(cc4,zz4,dum,ntry,'Propagating Mode Condition',
c     x  'Phase velocity (km/sec)','Surface Traction',2,0,.05,0,21)
c        do i=1,ntry
c          x4(i)=i
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
            if(iter.ge.0) then
c  get the group velocity via variational integrals
              call grvel(nl,cc,om,gv)
              nm=nm+1
c  these arrays are stored in units of m/sec
              cvel(nfrq+1-ifrq,nm)=cc*vbar
              gvel(nfrq+1-ifrq,nm)=gv*vbar
              nc(nm)=nc(nm)+1
  103 format(2g15.5)
            endif
          endif
        end do
        kroots(ifrq)=nm
      end do
      nroot=nm
c      print *,nroot,' = nroot'
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
        do kroot=1,nroot
          root=cvel(nfrq-ifrq+iifrq,kroot)/vbar
          gv=gvel(nfrq-ifrq+iifrq,kroot)/vbar
          roota(kroot)=root*(1.d0+ddfrq*(1.d0-root/(0.9*gv)))
          rootb(kroot)=root*(1.d0+ddfrq*(1.d0-root/(1.1*gv)))
c          roota(kroot)=root*(1.d0+ddfrq*(1.d0-root/(0.85*cmin)))
          rootb(kroot)=dmin1(root*1.00001d0,rootb(kroot))
          roota(kroot)=dmax1(roota(kroot),root-tolint)
          jtrval(kroot)=mtry
c  at very low f, roota can undershoot cmin -- problem addressed by tolint
        end do
c  we put another interval near cmax to pick up a new overtone branch
        ntrval=nroot+1
        rootb(ntrval)=cmax-eps
        roota(ntrval)=cmax-tolint
        jtrval(ntrval)=mtry
c  coalesce search intervals
        itrval=1
c        type *,ntrval
c        do i=1,ntrval
c          print *,roota(i),rootb(i)
c        end do
        do kroot=1,nroot
          if(rootb(itrval).gt.roota(itrval+1)) then
            ntrval=ntrval-1
            jtrval(itrval)=jtrval(itrval)+mtry
            rootb(itrval)=rootb(kroot+1)
            if(itrval.lt.nroot) then
              do i=itrval+1,nroot
                roota(i)=roota(i+1)
                rootb(i)=rootb(i+1)
              end do
            endif
          else           
            itrval=itrval+1
          endif
        end do
        ntrval=itrval
c        type *,ntrval
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
c  last check: horizontal waves in the top layer will have cc \appox gv
c  these are spurious, so kill them
                if(dabs(gv/cc-1.d0).gt.1.d-3) then
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
c            print *,'at f=',f,'Hz'
            ra=roota(nt)*vbar
            rb=rootb(nt)*vbar
c            print *,' in cvel interval',ra,rb
c            print *,nxroot,' expected, but',nmm,' roots found'
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
c              print *,'***!'
               nooutput=1
            endif
          endif
        end do
        kroots(ifrq)=nm
c        print *,kroots(ifrq),' = nroot,', ntrval,' = ntrval  at f=',f
c        if(ifrq.eq.(ifrq/50)*50) type *,nm,' modes at',f,' Hz'
      end do
  101 format(6g15.5)
  105 format('Modal Eigenfunction, f=',f3.1,' Hz,  c=',f7.4,' km/s')
c  after the cycle through frequencies, we plot the dispersion curves
      x4(1)=0.
      x4(2)=frqmax
      x4(3)=vss(1)*0.89*vbar/1000.d0
      x4(4)=vss(nlp)*1.01*vbar/1000.d0
      cmn=cmin*vbar
      cmx=cmax*vbar



c      open(7,file='aniprop_out.dat',form='formatted')
cc  first write out the model
cc      write(7) title
cc      write(7) nl
cc      do i=1,nlp
cc        write(7) (w(j,i),j=1,3),z(i),vp(i),vp2(i),
cc     x                               vp4(i),vs(i),vs2(i),rho(i)
cc      end do
cc  then the phase and group velocities
cc      write(7,*) nfrq,frqmax,cmn,cmx,(nc(i),i=1,maxbr)
c      print *,'output: nfrq,frqmax,cmn,cmx',nfrq,frqmax,cmn,cmx
c      print *,'maxbr:',maxbr
c      print *,(nc(i),i=1,10)
c      write(6,1099)(nc(i),i=1,maxbr)
c 1099 format(20i5)
cc      do imode=maxbr,1,-1
cc      do imode=2,1,-1
c       do imode=1,2
c        if(nc(imode).gt.1) then
cc          write(7) (cvel(i,imode),i=1,nc(imode))
c            do i=1,nfrq
c                frqq(i)=float(nfrq+1-i)*frqmax/float(nfrq)
c                write(7,*) 1.0/frqq(i),cvel(i,imode),gvel(i,imode)
c            end do 
c        endif
c      end do
cc      do imode=maxbr,1,-1
c      do imode=4,1,-1
c        if(nc(imode).gt.1) then
cc          write(7) (gvel(i,imode),i=1,nc(imode))
c          do i=1,nfrq
c             frqq(i)=float(nfrq+1-i)*frqmax/float(nfrq)
c             write(7,*) frqq(i),gvel(i,imode)
c          end do 
c        endif
c      end do
cc      write(7,*) (kroots(i),i=1,nfrq)
c      close(7)
      
      do i=1,nfrq
         frqq(i)=float(nfrq+1-i)*frqmax/float(nfrq)
         period(i)=1.0/frqq(i)
         
         rphase(i)=cvel(i,1)
         rgroup(i)=gvel(i,1)

         lphase(i)=cvel(i,2)
         lgroup(i)=gvel(i,2)
      end do 



c      call plotit_axes(x4(1),x4(2),x4(3),x4(4))
c      do imode=maxbr,1,-1
c        if(nc(imode).gt.1) then
c          np=nc(imode)
c          do i=1,np
c            zz4(i)=gvel(i,imode)/1000.    ! copy into real*4 array for plotting
c          end do
c          call plotit(frqq,zz4,dum,np,title,'Frequency (Hz)',
c     x        'Group Velocity (km/s)',2,0,0.05,2,1/imode)
c        endif
c      end do
c      x4(3)=vss(1)*vbar*0.79/1000.d0
c      call plotit_axes(x4(1),x4(2),x4(3),x4(4))
c      do imode=maxbr,1,-1
c        if(nc(imode).gt.1) then
c          np=nc(imode)
c          do i=1,np
c            zz4(i)=cvel(i,imode)/1000.    ! copy into real*4 array for plotting
c          end do
c          call plotit(frqq,gvel(1,imode),dum,np,title,'Frequency (Hz)',
c     x        'Phase Velocity (km/s)',2,0,0.05,2,1/imode)
c        endif
c      end do
c      stop
      end subroutine aniprop_subroutines

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
c        if(iter.ge.99) print *,cc,iter,zz
        if(iter.ge.99) nooutput=1
      else
c  flag to discard result, but test again to be certain 
c  if one a gridpoint was the root, there were no iterations, and zz is crap
c        print *,'failed?',cc,iter,zz
        call zzget(nl,cc,om,zz,iev)
c        type *,'recalculated',cc,zz
        if(dabs(dreal(zz))+dabs(dimag(zz)).gt.tol2) then
          iter=-1
c          type *,'failed again',cc,zz
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
c        type *,'a,b,c,d,e',a,b,c,d,e
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
c        type *,'t,s,r'
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
c        type *,'for layer',n
c        type *, 'for phase velocity',cc,'  the vertical wavenumbers are'
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
c        type *,'for layer',n,'  the vert wavenums are (0=up,1=dn)'
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
c          type *,'for i*k_z:',xnu(i,n),', the disp-stress vector is'
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
c  and solve for eigenvector in subroutine defective_aniprop
              xnu(i,n)=(xnu(i,n)+xnu(j,n))/2.d0
              xnu(j,n)=xnu(i,n)
c              zz=zz/(zabs(zz))
c              sq2=dsqrt(2.d0)
c              do jj=1,6
c                ee(jj,i,n)=(zz*ee(jj,i,n)+ee(jj,j,n))/sq2
c              end do
c  calculate the extravector for defective_aniprop repeated eigenvalue
              call defective_aniprop(i,j,n,adf(1,n),a,b,c,d,e,om,akx)
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
c  calculate the extravector for defective_aniprop repeated eigenvalue
c  as well as coefficient (adf) of linear (z-z0) term
              call defective_aniprop(i,j,n,adf(2,n),a,b,c,d,e,om,akx)
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
c  the possibility of defective_aniprop matrices must be contemplated here
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
c  the possibility of defective_aniprop matrices must be contemplated here
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
      subroutine defective_aniprop(i,j,n,adf,a,b,c,d,e,om,akx)
c  patch for dealing with nearly defective_aniprop propagator matrices
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
c  we wish to find the eigenvector of the near-defective_aniprop matrix 
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
c      real*4 cc4,zz4,x4
      complex*16 pp,u0,ee,z1,z0,zz,xnu,zzz,e1,e2,zla,xl,eigf
      complex*16 rt,tt,rt0,trc,co,ci1,ci2,ci3,ci4,dci1,dci2,dci3,dci4
      complex*16 u,xxnu,eye,pfac,uu,wukk,wuk,dci31,dci32
      real*4 cc4,zz4,x4
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
c      type *,nl,cc,om
      nlp=nl+1
c  update all the layer calculations
      call zzget(nl,cc,om,zz,iev)
c characteristics of the solution in the first layer
c            do k=1,6
c              type *,'for slowness:',xnu(k,1),', disp-stress vector'
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
          pfac(k,1)=zexp(xnu(k,il)*h1)
          if(il.lt.nlp) then
            pfac(k,2)=zexp(xnu(k,il)*h2)
          else
            pfac(k,2)=z0
          endif
c  normalize for the coefficient amplitudes in the Chen algorithm
c  upgoing coefs are referenced to layer bottom
c  downgoing coefs are referenced to layer top
          pfac(k,3)=pfac(k,1+(k-1)/3)
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
            else
              xxnu=(pfac(kk,2)*conjg(pfac(k,2))
     x           -pfac(kk,1)*conjg(pfac(k,1)))/(2.d0*xxnu)
            endif
c  normalize for the coefficient amplitudes in the Chen algorithm
c  upgoing coefs are referenced to layer bottom
c  downgoing coefs are referenced to layer top
            xxnu=xxnu/(pfac(kk,3)*conjg(pfac(k,3)))
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
c      type *,k,kk,uu,xxnu,a,d
c      type *,'dci31',dci31
c      type *,'dci32',dci32
c      type *,'dci3',dci3
c      type *,'dci4',dci4
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
c      type *,'ci1',ci1
c      type *,'ci2',ci2
c      type *,'ci3',ci3
c      type *,'ci4',ci4
      gv=(ci2/cc+ci3/(2.d0*om))/ci1
      cvel=cc*vbar/1000.
      gvel=gv*vbar/1000.
c      print 106,cvel,gvel
  106 format('phase velocity',f9.5,'  group velocity',f9.5)
  101 format(6g15.5)
  105 format('Modal Eigenfunction, f=',f3.1,' Hz,  c=',f7.4,' km/s')
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
                                                                
     

