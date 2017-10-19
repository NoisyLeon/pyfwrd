# -*- coding: utf-8 -*-
"""
Module for synthetic receiver function computation,
    using anirec by Vadim Levin, theo by Takuo Shibutani, and raysum by Andrew Frederiksen

:Copyright:
    Author: Lili Feng
    Graduate Research Assistant
    CIEI, Department of Physics, University of Colorado Boulder
    email: lili.feng@colorado.edu
"""

import aniprop, theo, raysum
import numba
import numpy as np
import vmodel
import copy
import matplotlib.pyplot as plt

# define type of vmodel.model1d
model_type = numba.deferred_type()
model_type.define(vmodel.model1d.class_type.instance_type)

# @numba.jit(numba.float32[:](numba.float32[:], numba.int32, numba.float32, numba.float32))
def _phaseshift( x, nfft, DT, TSHIFT ):
    """Add a shift to the data into the freq domain, private function for _iter_deconvolve
    """
    Xf      = np.fft.fft(x)
    # phase shift in radians
    shift_i = round(TSHIFT/DT) # removed +1 from here.
    p       = np.arange(nfft)+1
    p       = 2*np.pi*shift_i/(nfft)*p
    # apply shift
    Xf      = Xf*(np.cos(p) - 1j*np.sin(p))
    # back into time
    x       = np.real( np.fft.ifft(Xf) )/np.cos(2*np.pi*shift_i/nfft)
    return x

# @numba.jit(numba.float32[:](numba.float32[:], numba.float32[:], numba.float32))
def _FreFilter(inW, FilterW, dt ):
    """Filter input array in frequency domain, private function for _iter_deconvolve
    """
    FinW    = np.fft.fft(inW)
    FinW    = FinW*FilterW*dt
    FilterdW= np.real(np.fft.ifft(FinW))
    return FilterdW

# @numba.jit(numba.float32[:](numba.float32, numba.int32, numba.float32))
def _gaussFilter( dt, nft, f0 ):
    """
    Compute a gaussian filter in the freq domain which is unit area in time domain
    private function for _iter_deconvolve
    ================================================================================
    Input:
    dt  - sampling time interval
    nft - number freq points
    f0  - width of filter
    
    Output:
    gauss  - Gaussian filter array (numpy)
    filter has the form: exp( - (0.5*w/f0)^2 ) the units of the filter are 1/s
    ================================================================================
    """
    df      = 1.0/(nft*dt)
    nft21   = int(0.5*nft + 1)
    # get frequencies
    f       = df*np.arange(nft21, dtype=np.float32)
    w       = 2*np.pi*f
    w       = w/f0
    kernel  = w**2
    # compute the gaussian filter
    gauss   = np.zeros(nft, dtype=np.float32)
    gauss[:nft21]   = np.exp( -0.25*kernel )/dt
    gauss[nft21:]   = np.flipud(gauss[1:nft21-1])
    return gauss

# @numba.jit(numba.float32[:](numba.float32[:], numba.float32[:], numba.float32, numba.int32, numba.int32, numba.float32, numba.float32, numba.float32))
def _iter_deconvolve(Ztr, RTtr, dt, npts, niter, tdel, f0, minderr):
    """
    Iterative deconvolution
    """
    RMS         = np.zeros(niter, dtype = np.float32)  # RMS errors
    nfft        = 2**(npts-1).bit_length() # number points in fourier transform
    P0          = np.zeros(nfft, dtype = np.float32)# predicted spikes
    # Resize and rename the numerator and denominator
    U0          = np.zeros(nfft, dtype = np.float32) #add zeros to the end
    W0          = np.zeros(nfft, dtype = np.float32)
    U0[:npts]   = RTtr  # clear UIN;
    W0[:npts]   = Ztr   # clear WIN;
    # get filter in Freq domain 
    gauss       = _gaussFilter( dt, nfft, f0 )
    # filter signals
    Wf0         = np.fft.fft(W0)
    FilteredU0  = _FreFilter(U0, gauss, dt )
    FilteredW0  = _FreFilter(W0, gauss, dt )
    R           = FilteredU0 #  residual numerator
    # Get power in numerator for error scaling
    powerU      = np.sum(FilteredU0**2)
    # Loop through iterations
    it          = 0
    sumsq_i     = 1
    d_error     = 100*powerU + minderr
    maxlag      = int(0.5*nfft)
    while( abs(d_error) > minderr  and  it < niter ):
        it          = it+1 # iteration advance
        #   Ligorria and Ammon method
        RW          = np.real(np.fft.ifft(np.fft.fft(R)*np.conj(np.fft.fft(FilteredW0))))
        sumW0       = np.sum(FilteredW0**2)
        RW          = RW/sumW0
        imax        = np.argmax(abs(RW[:maxlag]))
        amp         = RW[imax]/dt; # scale the max and get correct sign
        #   compute predicted deconvolution
        P0[imax]    = P0[imax] + amp  # get spike signal - predicted RF
        P           = _FreFilter(P0, gauss*Wf0, dt*dt ) # convolve with filter
        #   compute residual with filtered numerator
        R           = FilteredU0 - P
        sumsq       = np.sum(R**2)/powerU
        RMS[it-1]   = sumsq # scaled error
        d_error     = 100*(sumsq_i - sumsq)  # change in error 
        sumsq_i     = sumsq  # store rms for computing difference in next   
    # Compute final receiver function
    P   = _FreFilter(P0, gauss, dt )
    # Phase shift
    P   = _phaseshift(P, nfft, dt, tdel)
    # output first nt samples
    RFI = P[:npts]
    # output the rms values 
    RMS = RMS[:it]
    if it > 1:
        return RFI, (1.0-RMS[it-1])*100.0
    else:
        return RFI, (1.-d_error)*100.


class ref_solver(object):
    """
    An object solving for receiver function using anirec, theo, and raysum
    =====================================================================================================================
    ::: parameters :::
    model           - 1D Earth model object
    dt              - time interval
    dArr            - layer array (unit - km)
    =====================================================================================================================
    """
    def __init__(self, inmodel):
        if not isinstance(inmodel, vmodel.model1d):
            raise ValueError('Input model should be type of vmodel.model1d !')
        self.model  = inmodel
        self.dt     = 0.05
        self.rfrst  = []; self.rftst  = []
        self.bazArr  = np.array([])
        ##############################
        # input parameters for raysum
        ##############################
        self.mults  = 2 # Multiples: 0 for none, 1 for Moho, 2 for all first-order
        self.width  = 0.4 # Gaussian width
        self.align  = 1 # Alignment: 0 is none, 1 aligns on primary phase (P or S)
        self.shift  = 5. # Shift of traces -- t=0 at this time (sec)
        self.outrot = 1  # Rotation to output: 0 is NS/EW/Z, 1 is R/T/Z, 2 is P/SV/SH
        return
    
    def init_default(self, dh=1., nl=100):
        self.dArr   = np.ones(nl, dtype = np.float32)*np.float32(dh)
        return
    
    def init_default_2(self):
        self.dArr   = np.array([20.,  15.,  42.,  43.,  45.,  35.], dtype = np.float32)
        return
    
    def init_default_3(self):
        self.dArr   = np.array([30.,  170.], dtype = np.float32)
        return
    
    def solve_theo(self, t=25., slowness = 0.06, din = None):
        """
        Compute radial and transverse receiver function using theo
        ====================================================================================
        ::: input parameters :::
        t           - time length of output in sec
        slowness    - reference horizontal slowness (default - 0.06 s/km, 1./0.06=16.6667)
        din         - incident angle(default - None, unit - deg),
                        the value is determined from slowness if not assigned
        ::: output :::
        self.rf     - (radial) receiver function
        self.time   - time array
        ====================================================================================
        """
        if not self.model.is_iso():
            print 'WARNING: model is anisotropic!'
        dArr, rhoArr, AArr, CArr, FArr, LArr, NArr = self.model.get_layer_model(self.dArr, 100, 1.)
        hin         = np.zeros(100, dtype=np.float32)
        vsin        = np.zeros(100, dtype=np.float32)
        vpvs        = np.zeros(100, dtype=np.float32)
        qsin        = 1400.*np.ones(100, dtype=np.float32)
        qpin        = 600.*np.ones(100, dtype=np.float32)
        if dArr.size<100:
            nl      = dArr.size
        else:
            nl      = 100
        hin[:nl]    = dArr[:nl]
        vsin[:nl]   = (np.sqrt(LArr/rhoArr))[:nl]
        vpvs[:nl]   = (np.sqrt(CArr/LArr))[:nl]
        fs          = 1./self.dt
        ntimes      = int(t*fs)
        if din is None:
            din     = 180.*np.arcsin(vsin[nl-1]*vpvs[nl-1]*slowness)/np.pi
        # solve for receiver function using theo
        rx 	        = theo.theo(nl, vsin, hin, vpvs, qpin, qsin, fs, din, 2.5, 0.005, 0, ntimes)
        # receiver function (ONLY radial component)
        self.rf     = rx[:ntimes]
        self.time   = np.arange(ntimes, dtype=np.float32)*self.dt 
        return
    
    def solve_anirec(self, t=25., baz=0., savestream=True):
        """
        Compute radial and transverse receiver function using anirec
        Default maximum velocity is 16.6667, can be changed by modifying the source code
        c     phase velocity of incident wave, important! LF 
            cc=16.6667  (line 175)
        ===================================================================================================================================
        ::: input parameters :::
        t           - time length of output in sec
        baz         - back-azimuth of wave vector 
        ::: output :::
        self.rfr    - radial receiver function
        self.rft    - transverse receiver function
        self.time   - time array
        ===================================================================================================================================
        """
        self.model.aniprop_check_model()
        z, rho, vp0, vp2, vp4, vs0, vs2 = self.model.layer_aniprop_model(self.dArr, 200, 1.)
        if self.model.tilt:
            self.dip, self.strike = self.model.angles_aniprop_model(z)
        # frqmax  = 1.0;  nfrq = 512;          df = frqmax/nfrq
        # npad    = 8192; dt   = 1./(npad*df)
        ntimes  = int(t/self.dt)
        nl      = z.size - 1
        if self.model.tilt:
            theta               = self.dip
            phig                = np.zeros(nl+1, dtype=np.float32)
            # phig[self.dip>0.]   = self.strike[self.dip>0.] + 270.
            phig[self.dip>0.]   = self.strike[self.dip>0.] + 90.
            phig[phig>=360.]    = phig[phig>=360.] - 360.
        else:
            theta   = np.zeros(nl+1, dtype=np.float32)
            phig    = np.zeros(nl+1, dtype=np.float32)
        # # baz     = 180. + az
        # # if baz > 360.:
        # #     baz -= 360.
        z       *= 1000.
        rho     *= 1000.
        vp0     *= 1000.
        vs0     *= 1000.
        # solve for receiver function using anirec
        Rf,Tf,T     = aniprop.rf_aniso_interface(z,vp0,vp2,vp4,vs0,vs2,rho,theta,phig,nl,baz,ntimes)
        # radial component
        self.rfr    = Rf
        # transverse component
        self.rft    = Tf
        # time
        self.time   = T
        dt2         = T[1] - T[0]
        if self.dt  != dt2:
            raise ValueError('Incompatible time interval from anirec, change the source code!')
        if savestream:
            self.rfrst.append(Rf); self.rftst.append(Tf)
            self.bazArr  = np.append(self.bazArr, baz)
        return
    
    
    
    def solve_anirec_reproduce(self, baz=0., t=20.):
        """
        Compute radial and transverse receiver function using aniprop
        Default maximum velocity is 16.6667, can be changed by modifying the source code
        
        c     phase velocity of incident wave, important! LF 
            cc=16.6667  (line 175)
        ==============================================================================================
        ::: input parameters :::
        baz         - back-azimuth
        t           - time length of output in sec
        ::: output :::
        self.rfr    - radial receiver function
        self.rft    - transverse receiver function
        self.time   - time array
        ==============================================================================================
        """

        z       = np.array([30., 200.])
        rho     = np.array([2.7, 3.3])
        vp0     = np.array([6.5, 8.0])
        vs0     = np.array([3.6, 4.6])
        vs2     = np.array([0.05, 0])
        vp2     = np.array([0.05, 0])
        vp4     = np.array([0.0, 0])
        z       *= 1000.
        rho     *= 1000.
        vp0     *= 1000.
        vs0     *= 1000.
        
        # frqmax  = 1.0;  nfrq = 512;          df = frqmax/nfrq
        # npad    = 8192; dt   = 1./(npad*df)
        ntimes  = int(t/self.dt)
        nl      = z.size - 1
        theta   = np.zeros(nl+1, dtype=np.float32)
        phig    = np.zeros(nl+1, dtype=np.float32)
        theta[0]= 90.
        # solve for receiver function using anirec
        Rf,Tf,T     = aniprop.rf_aniso_interface(z,vp0,vp2,vp4,vs0,vs2,rho,theta,phig,nl,baz,ntimes)
        # radial component
        self.rfr    = Rf
        # transverse component
        self.rft    = Tf
        # time
        self.time   = T
        self.rfrst.append(Rf); self.rftst.append(Tf)
        self.bazArr  = np.append(self.bazArr, baz)
        return
    
    def solve_raysum(self, bazin=np.array([0.]), t=25., iphase=1, slowness=0.06, phfname='', ref=True):
        """
        Compute synthetics/receiver functions using raysum
        ===================================================================================================================================
        ::: input parameters :::
        bazin       - bazimuth array of wave vector
        t           - time length of output in sec
        iphase      - initial phase index (1 - P; 2 - SV; 3 - SH)
        slowness    - reference horizontal slowness (default - 0.06 s/km, 1./0.06=16.6667)
        phfname     - phase list file name for output(or input if self.mults == 3)
        ::: output :::
        self.tt     - travel time array (shape - (nphase, ntrace))
        self.amp    - amplitude array (shape - (3, nphase, ntrace))
        self.trENZ  - ENZ component traces (shape - (3, npts, ntrace))
        self.trROT  - rotated component traces (shape - (3, npts, ntrace)), will not be available if self.outrot = 0
        ===================================================================================================================================
        """
        dArr    = np.zeros(self.dArr.size+1, dtype=np.float32)
        dArr[1:]= self.dArr # first layer is zero
        din, rhoin, alphain, betain, dvpin, dvsin, isoin = self.model.layer_raysum_model(dArr, 15, 1.)
        nl      = din.size
        if nl > 14:
            raise ValueError('Maximum allowed number of layers is 15!')
        # initialize model arrays
        d       = np.zeros(15, dtype=np.float32)
        rho     = np.zeros(15, dtype=np.float32)
        alpha   = np.zeros(15, dtype=np.float32)
        beta    = np.zeros(15, dtype=np.float32)
        iso     = np.ones(15, dtype=np.int32)
        dvp     = np.zeros(15, dtype=np.float32)
        dvs     = np.zeros(15, dtype=np.float32)
        trend   = np.zeros(15, dtype=np.float32)
        plunge  = np.zeros(15, dtype=np.float32)
        strike  = np.zeros(15, dtype=np.float32)
        dip     = np.zeros(15, dtype=np.float32)
        #
        d[:nl]      = din[:]*1000.
        rho[:nl]    = rhoin[:]*1000.
        alpha[:nl]  = alphain[:]*1000.
        beta[:nl]   = betain[:]*1000.
        iso[:nl]    = isoin[:]
        dvp[:nl]    = dvpin[:]*100.
        dvs[:nl]    = dvsin[:]*100.
        # bottom half space
        nl          += 1
        d[nl-1]     = 0.
        rho[nl-1]   = rho[nl-2]
        alpha[nl-1] = alpha[nl-2]
        beta[nl-1]  = beta[nl-2]
        iso[nl-1]   = 1
        # topmost layer
        iso[0]      = 1
        dvp[0]      = 0.
        dvs[0]      = 0.
        
        if self.model.tilt:
            self.dip, self.strike = self.model.angles_raysum_model(din, 0)
            trend[:nl-1]   = self.strike[:]+270.; plunge[:nl-1] = 90. - self.dip[:] # double check
        if self.model.dipping:
            self.dipif, self.strikeif = self.model.angles_raysum_model(din, 1)
            dip[1:nl] = self.dipif[:]; strike[1:nl] = self.strikeif[:]
        # top most layer
        trend[0]    = 0.; plunge[0]=0.
        bazin       = np.asarray(bazin)
        ntr         = bazin.size
        baz         = np.zeros(200, dtype=np.float32);  baz[:ntr]   = bazin[:]
        slow        = np.zeros(200, dtype=np.float32);  slow[:ntr]  = slowness/1000. # s/km to s/m
        sta_dx      = np.zeros(200, dtype=np.float32)
        sta_dy      = np.zeros(200, dtype=np.float32)
        self.npts   = int(t/self.dt)
        t           += self.shift
        npts        = int(t/self.dt)
        self.nptsraysum = npts
        # Compute synthetics using raysum
        tt, amp, nphase, tr_cart, tr_ph = raysum.raysum_interface(nl, d, rho, alpha, beta, dvp, dvs, \
                    trend, plunge, strike, dip, iso, iphase,   ntr, baz, slow, sta_dx, sta_dy, \
                        self.mults, npts, self.dt, self.width, self.align, self.shift, self.outrot, phfname)
        self.tt     = tt[:nphase, :ntr]
        self.amp    = amp[:, :nphase, :ntr]
        self.trENZ  = tr_cart[:, :npts, :ntr]
        if self.outrot != 0:
            self.trROT  = tr_ph[:, :npts, :ntr]
        self.nphase = nphase; self.ntr = ntr
        self.time   = np.arange(self.npts, dtype=np.float32)*self.dt
        self.bazArr = bazin
        if ref: self.deconvolve_raysum()
        return
    
    def solve_raysum_normal(self, bazin=np.array([0.]), t=30., iphase=1, slowness=0.06, phfname='', ref=True):
        """
        Compute radial and transverse receiver function using raysum
        ===================================================================================================================================
        ::: input parameters :::
        bazin       - bazimuth array of wave vector
        t           - time length of output in sec
        iphase      - initial phase index (1 - P; 2 - SV; 3 - SH)
        slowness    - reference horizontal slowness (default - 0.06 s/km, 1./0.06=16.6667)
        phfname     - phase list file name for output(or input if self.mults == 3)
        ::: output :::
        self.tt     - travel time array (shape - (nphase, ntrace))
        self.amp    - amplitude array (shape - (3, nphase, ntrace))
        self.trENZ  - ENZ component traces (shape - (3, npts, ntrace))
        self.trROT  - rotated component traces (shape - (3, npts, ntrace)), will not be available if self.outrot = 0
        ===================================================================================================================================
        """
        din, rhoin, alphain, betain, dvpin, dvsin, isoin = self.model.layer_raysum_model(self.dArr, 15, 1.)
        nl      = din.size
        if nl > 14:
            raise ValueError('Maximum allowed number of layers is 15!')
        # initialize model arrays
        d       = np.zeros(15, dtype=np.float32)
        rho     = np.zeros(15, dtype=np.float32)
        alpha   = np.zeros(15, dtype=np.float32)
        beta    = np.zeros(15, dtype=np.float32)
        iso     = np.ones(15, dtype=np.int32)
        dvp     = np.zeros(15, dtype=np.float32)
        dvs     = np.zeros(15, dtype=np.float32)
        trend   = np.zeros(15, dtype=np.float32)
        plunge  = np.zeros(15, dtype=np.float32)
        strike  = np.zeros(15, dtype=np.float32)
        dip     = np.zeros(15, dtype=np.float32)
        #
        d[:nl]      = din[:]*1000.
        rho[:nl]    = rhoin[:]*1000.
        alpha[:nl]  = alphain[:]*1000.
        beta[:nl]   = betain[:]*1000.
        iso[:nl]   = isoin[:]
        dvp[:nl]   = dvpin[:]*100.
        dvs[:nl]   = dvsin[:]*100.
        
        nl          += 1
        d[nl-1]     = 0.
        rho[nl-1]   = rho[nl-2]
        alpha[nl-1] = alpha[nl-2]
        beta[nl-1]  = beta[nl-2]
        iso[nl-1]   = 1
        # 
        # iso[1:nl]   = isoin[:]
        # dvp[1:nl]   = dvpin[:]*100.
        # dvs[1:nl]   = dvsin[:]*100.
        
        if self.model.tilt:
            self.dip, self.strike = self.model.angles_raysum_model(din, 0)
            # trend[1:nl]   = self.strike[:]+270.; plunge[1:nl] = 90. - self.dip[:] # double check
            trend[:nl-1]   = self.strike[:]; plunge[:nl-1] = self.dip[:] # double check
        if self.model.dipping:
            self.dipif, self.strikeif = self.model.angles_raysum_model(din, 1)
            dip[1:nl] = self.dipif[:]; strike[1:nl] = self.strikeif[:]
        bazin   = np.asarray(bazin)
        ntr     = bazin.size
        baz     = np.zeros(200, dtype=np.float32); baz[:ntr] = bazin[:]
        slow    = np.zeros(200, dtype=np.float32); slow[:ntr]= slowness/1000. # s/km to s/m
        sta_dx  = np.zeros(200, dtype=np.float32)
        sta_dy  = np.zeros(200, dtype=np.float32)
        
        npts    = int(t/self.dt)
        
        tt, amp, nphase, tr_cart, tr_ph = raysum.raysum_interface(nl, d, rho, alpha, beta, dvp, dvs, \
                    trend, plunge, strike, dip, iso, iphase,   ntr, baz, slow, sta_dx, sta_dy, \
                        self.mults, npts, self.dt, self.width, self.align, self.shift, self.outrot, phfname)
        self.tt     = tt[:nphase, :ntr]
        self.amp    = amp[:, :nphase, :ntr]
        self.trENZ  = tr_cart[:, :npts, :ntr]
        if self.outrot != 0:
            self.trROT  = tr_ph[:, :npts, :ntr]
        self.nphase = nphase; self.ntr = ntr
        self.time   = np.arange(npts, dtype=np.float32)*self.dt
        self.bazArr = bazin
        if ref: self.deconvolve_raysum()
        return
    
    def deconvolve_raysum(self, tdel=0., f0 = 2.5, niter=200, minderr=0.0001):
        """
        Compute receiver function from raysum synthetics with iterative deconvolution algorithmn
        ========================================================================================================================
        ::: input parameters :::
        tdel        - phase delay
        f0          - Gaussian width factor
        niter       - number of maximum iteration
        minderr     - minimum misfit improvement, iteration will stop if improvement between two steps is smaller than minderr
        ::: output :::
        self.rfrst  - list for radial receiver functions
        self.rftst  - list for transverse receiver functions
        ========================================================================================================================
        """
        if self.outrot != 1:
            print 'No RTZ synthetics from raysum!'
            return
        nptsraysum      = self.nptsraysum
        for i in xrange(self.ntr):
            Ztr         = self.trROT[2,:,i]
            Rtr         = self.trROT[0,:,i]
            Ttr         = self.trROT[1,:,i]
            rfr, fitness= _iter_deconvolve(Ztr, Rtr, self.dt, nptsraysum, niter, tdel, f0, minderr)
            if fitness < 95.:
                print 'WARNING: fittness is',fitness,'for trace id:',i
            if np.all(Ttr == 0.):
                rft     = np.zeros(Ztr.size, np.float32)
            else:
                rft, fitness    = _iter_deconvolve(Ztr, Ttr, self.dt, nptsraysum, niter, tdel, f0, minderr)
                if fitness < 95.:
                    print 'WARNING: fittness is',fitness,'for trace id:',i
            self.rfrst.append(rfr[:self.npts]); self.rftst.append(rft[:self.npts])
        return

    
    def compute_diff_ps_time(self, slowness = 0.06, h=35., ptype='v', stype='v', vs0=None, vp0=None, vp1=None):
        """Compute the difference in arrival time between P and P-S converted waves.
        """
        dArr, rhoArr, AArr, CArr, FArr, LArr, NArr = self.model.get_layer_model(self.dArr, 100, 1.)
        if vs0 == None:
            if stype == 'v':
                vs0     = (np.sqrt(LArr/rhoArr))[0]
            elif stype == 'h':
                vs0     = (np.sqrt(NArr/rhoArr))[0]
            else:
                raise ValueError('Unexpected stype!')
        if vp0 == None or vp1 == None:
            if ptype == 'v':
                vp0     = (np.sqrt(CArr/rhoArr))[0]
                vp1     = (np.sqrt(CArr/rhoArr))[1]
            elif ptype == 'h':
                vp0     = (np.sqrt(AArr/rhoArr))[0]
                vp1     = (np.sqrt(AArr/rhoArr))[1]
            else:
                raise ValueError('Unexpected ptype!')
        phis0   = np.arcsin(vs0*slowness)
        phip0   = np.arcsin(vp0*slowness)
        phip1   = np.arcsin(vp1*slowness)
        dtps    = h* (np.tan(phip0) - np.tan(phis0))*np.sin(phip1)/vp1 + h/np.cos(phis0)/vs0 - h/np.cos(phip0)/vp0
        print dtps
        
    def plot_baz_rf(self, comp='T', showfig=True):
        ymax=361.
        ymin=-1.
        time    = self.time
        plt.figure()
        ax=plt.subplot()
        for i in xrange(self.bazArr.size):
            if comp=='R':
                yvalue  = self.rfrst[i]
            else:
                yvalue  = self.rftst[i]
            rfmax   = yvalue.max()
            yvalue  = -yvalue/rfmax*10.
            # yvalue[(time>1.8)*(time<3.2)] = 2.*yvalue[(time>1.8)*(time<3.2)]
            baz     = self.bazArr[i]
            ax.plot(time, yvalue+baz, '-k', lw=0.3)
            ax.fill_between(time, y2=baz, y1=yvalue+baz, where=yvalue>0, color='red', lw=0.01, interpolate=True)
            ax.fill_between(time, y2=baz, y1=yvalue+baz, where=yvalue<0, color='blue', lw=0.01, interpolate=True)
            plt.axis([0., 6, ymin, ymax])
            plt.xlabel('Time(sec)')
            plt.title(comp+' component')
            plt.gca().invert_yaxis()
        if showfig: plt.show()
    
    
    