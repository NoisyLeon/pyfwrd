# -*- coding: utf-8 -*-
"""
Module for synthetic shear wave splitting computation,
    using raysum by Andrew Frederiksen

The code is a python wrapper of the f77 code aniprop and theo

:Copyright:
    Author: Lili Feng
    Graduate Research Assistant
    CIEI, Department of Physics, University of Colorado Boulder
    email: lili.feng@colorado.edu
"""

import raysum
import numba
import numpy as np
import vmodel
import copy
import matplotlib.pyplot as plt
from obspy.signal.rotate import rotate_ne_rt 

# define type of vmodel.model1d
model_type = numba.deferred_type()
model_type.define(vmodel.model1d.class_type.instance_type)

class sws_solver(object):
    """
    An object solving for shear wave splitting using raysum
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
        self.dt     = 0.025
        self.delayt = np.array([])
        self.phi    = np.array([])
        self.bazArr = np.array([])
        ##############################
        # input parameters for raysum
        ##############################
        self.mults  = 0 # Multiples: 0 for none, 1 for Moho, 2 for all first-order
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
        # self.dArr   = np.array([20.,  15.,  42., 43.], dtype = np.float32)
        return
    
    def init_default_3(self):
        self.dArr   = np.array([30.,  170.], dtype = np.float32)
        return
    
    def solve_raysum(self, bazin=np.array([0.]), t=30., iphase=2, slowness=0.06, phfname='', sws=True):
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
        # nl          += 1
        # d[nl-1]     = 0.
        # rho[nl-1]   = rho[nl-2]
        # alpha[nl-1] = alpha[nl-2]
        # beta[nl-1]  = beta[nl-2]
        # iso[nl-1]   = 1
        # topmost layer
        iso[0]      = 1
        dvp[0]      = 0.
        dvs[0]      = 0.
        
        if self.model.tilt:
            self.dip, self.strike = self.model.angles_raysum_model(din, 0)
            # trend[:nl-1]   = self.strike[:]+270.; plunge[:nl-1] = 90. - self.dip[:] # double check
            trend[:nl]   = self.strike[:]+270.; plunge[:nl] = 90. - self.dip[:] # double check
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
        # Compute synthetics using raysum
        tt, amp, nphase, tr_cart, tr_ph = raysum.raysum_interface(nl, d, rho, alpha, beta, dvp, dvs, \
                    trend, plunge, strike, dip, iso, iphase,   ntr, baz, slow, sta_dx, sta_dy, \
                        self.mults, self.npts, self.dt, self.width, self.align, self.shift, self.outrot, phfname)
        self.tt     = tt[:nphase, :ntr]
        self.amp    = amp[:, :nphase, :ntr]
        self.trENZ  = tr_cart[:, :self.npts, :ntr]
        if self.outrot != 0:
            self.trROT  = tr_ph[:, :self.npts, :ntr]
        self.nphase = nphase; self.ntr = ntr
        self.time   = np.arange(self.npts, dtype=np.float32)*self.dt
        self.bazArr = bazin

        return
    
    def rotate(self, trid, angle):
        angle   = -angle/180.*np.pi
        x       = self.trROT[1, :, trid]
        y       = self.trROT[0, :, trid]
        
        compx   = x*np.cos(angle) - y* np.sin(angle)
        compy   = x*np.sin(angle) + y* np.cos(angle)
        return compx, compy

    
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
