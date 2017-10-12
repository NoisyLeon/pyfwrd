c Calculates a spread of seismograms for a given model.
c Requires readmodel.f,...
c Usage: seis-spread modelfile geometryfile phasefile arrivalfile

c####&

      subroutine raysum_interface(nlay, thick, rho, alpha, beta, 
     x      pct_a, pct_b, trend, plunge, strike, dip, isoflag,
     x      iphase_in, ntr, baz, slow, sta_dx, sta_dy,
     x      mults,nsamp,dt,width,align,shift,out_rot, phname_in,
     x      travel_time,amplitude,numph, Tr_cart, Tr_ph)

c Anyone who fails to start a Fortran program with this line
c should be severely beaten:
        implicit none
    
c Include constants:
        include 'params.h'

c Scratch variables:
        integer i,j,iargc
        character modname*(namelen),geomname*(namelen),phname*(namelen)
        character arrname*(namelen),tracename*(namelen)
        character iphname*(namelen)
        real amp_in
c Model parameters:
        integer nlay
        real thick(maxlay),rho(maxlay),alpha(maxlay),beta(maxlay)
        real pct_a(maxlay),pct_b(maxlay),trend(maxlay),plunge(maxlay)
        real strike(maxlay),dip(maxlay)
        logical isoflag(maxlay)
c Geometry parameters
        real baz(maxtr),slow(maxtr),sta_dx(maxtr),sta_dy(maxtr)
        integer ntr
c Phase parameters
        integer phaselist(maxseg,2,maxph),nseg(maxph),numph,iphase
c Arrivals
        real travel_time(maxph,maxtr),amplitude(3,maxph,maxtr)
c Traces
        real Tr_cart(3,maxsamp,maxtr),dt,width,shift
        real Tr_ph(3,maxsamp,maxtr)
        integer align,nsamp,mults,out_rot
c Input
        integer iphase_in
        character phname_in*(namelen)
        
c   aa is a list of rank-4 tensors (a_ijkl = c_ijkl/rho)
c   rot is a list of rotator matrices, used to rotate into the local
c   coordinate system of each interface.
        real aa(3,3,3,3,maxlay), rot(3,3,maxlay)
c   ar_list is aa, pre-rotated into reference frames for the interfaces
c   above (last index=2) and below (1) each respective layer.
        real ar_list(3,3,3,3,maxlay,2)
        
        
c Get filenames from command-line argument.

c        i = iargc()
c        if (i .lt. 5) then
c          write(*,*) 'Usage: seis-spread modelfile geometryfile ',
c     &               'phasefile arrivalfile tracefile [iphase]'
c          write(*,*) '-- Modelfile and geometryfile must already exist.'
c          write(*,*) '-- Other files will be overwritten.'
c          write(*,*) '-- iphase is P by default; may also be SV or SH.'
c          stop
c        end if
c        call getarg(1,modname)
c        write(*,*) 'Model is ',modname

c Determine initial phase
        iphase = iphase_in
        if (iphase .eq. 1) then
           iphname = 'P'
        else if (iphase .eq. 2) then
           iphname = 'SV'
        else if (iphase .eq. 3) then
           iphname = 'SH'
        else
            write (*,*) iphase,' is not a valid phase type.'
            write (*,*) 'Valid phases are 1-P, 2-SV and 3-SH.'
            stop
        end if
        write (*,*) 'Initial phase is ',iphname
        
c        if (i .eq. 5) then
c          iphase=1
c          iphname='P'
c        else
c          call getarg(6,iphname)
c          if (iphname .eq. 'P') then
c            iphase=1
c          else if (iphname .eq. 'SV') then
c            iphase=2
c          else if (iphname .eq. 'SH') then
c            iphase=3
c          else
c            write (*,*) iphname,' is not a valid phase type.'
c            write (*,*) 'Valid phases are P, SV and SH.'
c            stop
c          end if
c        end if
c       write (*,*) 'Initial phase is ',iphname
      
c Read in model      
c        call readmodel(modname,thick,rho,alpha,beta,isoflag,
c     &                 pct_a,pct_b,trend,plunge,strike,dip,nlay)
     
        do i=1,nlay
          strike(i)=strike(i)/180. * pi
          dip(i)=dip(i)/180. * pi
          trend(i)=trend(i)/180. * pi
          plunge(i)=plunge(i)/180. * pi
        end do
     
        call writemodel(6,thick,rho,alpha,beta,isoflag,
     &                  pct_a,pct_b,trend,plunge,strike,dip,nlay)
     
c Set up model for calculation, make rotators
        call buildmodel(aa,ar_list,rot,thick,rho,alpha,beta,isoflag,
     &                  pct_a,pct_b,trend,plunge,strike,dip,nlay)
c Read in geometry (desired traces)
c        call getarg(2,geomname)
c        write (*,*) 'Geometry is ',geomname
c        call readgeom(geomname,baz,slow,sta_dx,sta_dy,ntr)
        do i=1,ntr
          baz(i)=baz(i)/180. * pi
        end do
        call writegeom(6,baz,slow,sta_dx,sta_dy,ntr)

c Read in parameters from file 'raysum-params', if it exists
c        geomname='raysum-params'
c        call readparams(geomname,mults,nsamp,dt,width,align,
c     &                  shift,out_rot)

c Generate phase list
c        call getarg(3,phname)
        phname = phname_in
        numph=0
        if (mults .ne. 3) then
          call ph_direct(phaselist,nseg,numph,nlay,iphase)
        end if
        if (mults .eq. 1) then
          call ph_fsmults(phaselist,nseg,numph,nlay,1,iphase)
        end if
        if (mults .eq. 2) then
          do j=1,nlay-1
            call ph_fsmults(phaselist,nseg,numph,nlay,j,iphase)
          end do
        end if
c   Read phases from phase list file if mults == 3
        if (mults .eq. 3) then
          write(*,*) 'Reading phases from file ',phname
          call readphases(phname,phaselist,nseg,numph)
        end if
        
c        call printphases(phaselist,nseg,numph)
c   Write phases to phase list file if phname is not ''
        if (mults .ne. 3) then
            if (phname.ne.'') then
                open(unit=iounit1,file=phname,status='unknown')
                call writephases(iounit1,phaselist,nseg,numph)
                close(unit=iounit1)
                write(*,*) 'Phases written to ',phname
            end if
        end if

c        call getarg(4,arrname)
c        write(*,*) 'Arrivals will be written to ',arrname
c        open(unit=iounit1,file=arrname,status='unknown')
        
c        call getarg(5,tracename)
c        write(*,*) 'Traces will be written to ',tracename
c        open(unit=iounit2,file=tracename,status='unknown')
     
c        Perform calculation                   
        amp_in=1.
        call get_arrivals(travel_time,amplitude,thick,rho,isoflag,
     &       strike,dip,aa,ar_list,rot,baz,slow,sta_dx,sta_dy,
     &       phaselist,ntr,nseg,numph,nlay,amp_in)

c        Normalize arrivals
        if (iphase .eq. 1) then 
          call norm_arrivals(amplitude,baz,slow,alpha(1),beta(1),
     &                       rho(1),ntr,numph,1,1)
        end if

c        Write out arrivals
c        call writearrivals(iounit1,travel_time,amplitude,ntr,numph)
c        close(unit=iounit1)

c        Assemble traces
        call make_traces(travel_time,amplitude,ntr,numph,nsamp,
     &                   dt,width,align,shift,Tr_cart)
     
c        if (out_rot .eq. 0) then
c          call writetraces(iounit2,Tr_cart,ntr,nsamp,dt,align,shift)
c        else

        if (out_rot .eq. 1) then
c            Rotate to RTZ
            call rot_traces(Tr_cart,baz,ntr,nsamp,Tr_ph)
        else if (out_rot .eq. 2) then
c            Rotate to wavevector coordinates
            call fs_traces(Tr_cart,baz,slow,alpha(1),beta(1),
     &                     rho(1),ntr,nsamp,Tr_ph)
        end if
          
c          Write results
c          call writetraces(iounit2,Tr_ph,ntr,nsamp,dt,align,shift)
c        end if
c        close(unit=iounit2)
     
      end
      
      
      subroutine readparams(filename,mults,nsamp,dt,width,align,shift,
     &                      out_rot)
      
        implicit none
        include 'params.h'
        
        character filename*(namelen),buffer*(buffsize)
        integer mults,nsamp,align,ios,eof,out_rot
        real dt,width,shift
        
c          Default values
        mults=1
        nsamp=600
        dt=0.05
        width=1.
        align=1
        shift=5.
        out_rot=2
        
        open(unit=iounit1,status='old',file=filename,iostat=ios)
        
        if (ios .eq. 0) then
          write(*,*) 'Reading parameters from ',filename
          call getline(iounit1,buffer,eof)
          read (buffer,*) mults
          call getline(iounit1,buffer,eof)
          read (buffer,*) nsamp
          call getline(iounit1,buffer,eof)
          read (buffer,*) dt
          call getline(iounit1,buffer,eof)
          read (buffer,*) width
          call getline(iounit1,buffer,eof)
          read (buffer,*) align
          call getline(iounit1,buffer,eof)
          read (buffer,*) shift
          call getline(iounit1,buffer,eof)
          read (buffer,*) out_rot
        else
          write (*,*) 'Parameter file ',filename,' does not exist.'
          write (*,*) 'Writing default values.'
          close(unit=iounit1)
          open(unit=iounit1,status='unknown',file=filename)
          write(iounit1,*) '# Multiples: 0 for none, 1 for Moho, 2 ',
     &                      'for all first-order, 3 to read file'
          write(iounit1,*) mults
          write(iounit1,*) '# Number of samples per trace'
          write(iounit1,*) nsamp
          write(iounit1,*) '# Sample rate (seconds)'
          write(iounit1,*) dt
          write(iounit1,*) '# Gaussian pulse width (seconds)'
          write(iounit1,*) width
          write(iounit1,*) '# Alignment: 0 is none, 1 aligns on P'
          write(iounit1,*) align
          write(iounit1,*) '# Shift of traces -- t=0 at this time (sec)'
          write(iounit1,*) shift
          write(iounit1,*) '# Rotation to output: 0 is NS/EW/Z, ',
     &                     '1 is R/T/Z, 2 is P/SV/SH'
          write(iounit1,*) out_rot
        end if
        close(unit=iounit1)
                  
      end







