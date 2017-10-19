# -*- coding: utf-8 -*-
"""
Module for synthetic shear wave splitting computation,
    using raysum by Andrew Frederiksen

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


def RickerIntSignal(dt, npts, fc, t0=None):
    """RickerInt signal defined in sw4 manual(p 18)
    """
    time    = np.arange(npts)*dt
    if t0 == None:
        Nshift  = np.int(1.35/fc/dt) - 1
        t0      = Nshift *dt
    else:
        Nshift  = np.int(t0/dt) - 1
        t0      = Nshift *dt
    return -(time - t0) *np.exp(- (np.pi*fc*(time-t0) )**2 )

# @numba.jit( numba.types.UniTuple(numba.float64[:,:], 3)(numba.int64[:], numba.float64[:], numba.float64[:], numba.float64[:], numba.float64[:] ))
# def _silverchan_mat(ind, SG, SH, phi_test, dt_test):
#     Ematrix = np.zeros([phi_test.size, dt_test.size], dtype=np.float64) # energy matrix
#     l1      = np.zeros([phi_test.size, dt_test.size], dtype=np.float64) # lambda 1
#     l2      = np.zeros([phi_test.size, dt_test.size], dtype=np.float64) # lambda 2
#     Nphi    = phi_test.size
#     Nt      = dt_test.size
#     for iphi in xrange(Nphi):
#         phi     = phi_test[iphi]
#         xphi    = SH*np.cos(phi) - SG* np.sin(phi)
#         yphi    = SH*np.sin(phi) + SG* np.cos(phi)
#         tmpfast = yphi[ind]
#         for idt  in xrange(Nt):
#             dt      = dt_test[idt]
#             tmpslow = np.take(xphi, ind+idt) # xphi[ind+idt]
#             SHc     = tmpslow*np.cos(-phi) - tmpfast* np.sin(-phi)
#             SGc     = tmpslow*np.sin(-phi) + tmpfast* np.cos(-phi)
#             Ematrix[iphi, idt] = np.linalg.norm(SHc) # Energy on transverse component
#             # Etr     = SHc**2
#             # ELst.append(np.linalg.norm(SHc))
#             # ELst.append(Etr.sum()/Etr.size)
#             SHe     = SHc - SHc.mean()
#             SGe     = SGc - SGc.mean()
#             # construct covariance matrix
#             covar   = np.cov(SHe, SGe)
#             # eigenvalue estimation
#             w, v    = np.linalg.eig(covar)
#             l1[iphi, idt]   = w[1]
#             l2[iphi, idt]   = w[0]
#     #         l1.append(w[1])
#     #         l2.append(w[0])
#     # Ematrix = np.array(ELst)
#     # l1      = np.array(l1)
#     # l2      = np.array(l2)
#     return Ematrix, l1, l2
# 
# fast_silverchan_mat= numba.jit( numba.types.UniTuple(numba.float32[:,:], 3)(numba.float32[:],\
#                     numba.float32[:], numba.float32[:], numba.float32[:], numba.float32[:] )) (_silverchan_mat)

class sws_solver(object):
    """
    An object solving for shear wave splitting using raysum
    =====================================================================================================================
    ::: parameters :::
    model           - 1D Earth model object
    dt              - time interval
    dArr            - layer array (unit - km)
    bazArr          - back-azimuth array (unit - deg)
    delayt          - delay time
    phi             - fast axis angle relative to RTZ system
    fa_phi          - fast axis angle relative to ENZ system
    Crc, Cev, Cme   - matrix storing correlation coefficients/eigenvalues/transverse energy
    =====================================================================================================================
    """
    def __init__(self, inmodel):
        if not isinstance(inmodel, vmodel.model1d):
            raise ValueError('Input model should be type of vmodel.model1d !')
        self.model  = inmodel
        self.dt     = 0.025
        self.bazArr = np.array([])
        ##############################
        # input parameters for raysum
        ##############################
        self.mults  = 0 # Multiples: 0 for none, 1 for Moho, 2 for all first-order
        self.width  = .4 # Gaussian width
        self.align  = 1 # Alignment: 0 is none, 1 aligns on primary phase (P or S)
        self.shift  = 10. # Shift of traces -- t=0 at this time (sec)
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
    
    def solve_raysum(self, bazin=np.array([0.]), t=50., iphase=2, slowness=0.06, phfname='', sws=True):
        """
        Compute synthetics using raysum
        ===================================================================================================================================
        ::: input parameters :::
        bazin       - back-azimuth array of wave vector
        t           - time length of output in sec
        iphase      - initial phase index (1 - P; 2 - SV; 3 - SH)
        slowness    - reference horizontal slowness (default - 0.06 s/km, 1./0.06=16.6667)
        phfname     - phase list file name for output(or input if self.mults == 3)
        ::: output :::
        self.tt     - travel time array (shape - (nphase, ntrace))
        self.amp    - amplitude array (shape - (3, nphase, ntrace))
        self.trNEZ  - ENZ component traces (shape - (3, npts, ntrace))
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
        d           = np.zeros(15, dtype=np.float32)
        rho         = np.zeros(15, dtype=np.float32)
        alpha       = np.zeros(15, dtype=np.float32)
        beta        = np.zeros(15, dtype=np.float32)
        iso         = np.ones(15, dtype=np.int32)
        dvp         = np.zeros(15, dtype=np.float32)
        dvs         = np.zeros(15, dtype=np.float32)
        trend       = np.zeros(15, dtype=np.float32)
        plunge      = np.zeros(15, dtype=np.float32)
        strike      = np.zeros(15, dtype=np.float32)
        dip         = np.zeros(15, dtype=np.float32)
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
            # trend[:nl]   = self.strike[:]+270.; plunge[:nl] = 90. - self.dip[:] # double check
        if self.model.dipping:
            self.dipif, self.strikeif = self.model.angles_raysum_model(din, 1)
            dip[1:nl] = self.dipif[:]; strike[1:nl] = self.strikeif[:]
        # top most layer
        trend[0]    = 0.; plunge[0]=0.
        trend[trend>360.] = trend[trend>360.] - 360.
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
        self.trNEZ  = tr_cart[:, :self.npts, :ntr]
        if self.outrot != 0:
            self.trROT  = tr_ph[:, :self.npts, :ntr]
        self.nphase = nphase; self.ntr = ntr
        self.time   = np.arange(self.npts, dtype=np.float32)*self.dt
        self.bazArr = bazin
        ###
        # store model parameters
        ###
        self.d          = d[:nl-1]
        self.beta       = beta[:nl-1]
        self.alpha      = alpha[:nl-1]
        self.slowness   = slowness
        return
    
    def rotate(self, trid, angle, dtype=1):
        """
        Rotate the synthetics clockwise by a given angle
        ======================================================================================================
        ::: input parameters :::
        trid        - trace index
        angle       - angle for rotation 
        dtype       - datatype for rotation (1 - RTZ before convolution; other - RTZ after convolution)
        ::: output :::
        compx, compy- rotated x/y components
        ======================================================================================================
        """
        angle   = angle/180.*np.pi
        if dtype == 1:
            x       = self.trSYNROT[1, :, trid]
            y       = self.trSYNROT[0, :, trid]
        else:
            x       = self.trROT[1, :, trid]
            y       = self.trROT[0, :, trid]
        # NOTE: Rotation matrix for a vector counter-clockwise with an angle
        # It is equivalent to rotate the coordinate system clockwise with the same angle
        compx   = x*np.cos(angle) - y* np.sin(angle)
        compy   = x*np.sin(angle) + y* np.cos(angle)
        return compx, compy

    def convolve(self, fc=2.5, rot=True):
        """
        Convolve synthetics from raysum with a Ricker integral signal
        ====================================================================================
        ::: input parameters :::
        fc              - central frequency of the Ricker integral signal
        rot             - convolve with RTZ (True) or NEZ (false) components
        ::: output :::
        self.trSYNROT   - convolved synthetic seismograms for RTZ system
        self.trSYN      - convolved synthetic seismograms for NEZ system
        self.stf        - source time function used for convolution
        ====================================================================================
        """
        fs          = 1./self.dt
        npts        = min(fs/fc*8, self.npts)
        Nshift      = np.int(1.35/fc/self.dt) - 1
        t0          = Nshift *self.dt
        stf         = RickerIntSignal(dt=self.dt, npts = npts, fc=fc, t0=t0)
        self.trSYN  = np.zeros(self.trROT.shape)
        if rot:     self.trSYNROT  = np.zeros(self.trROT.shape)
        for i in xrange(self.ntr):
            if rot:
                self.trSYNROT[0, :, i]  = np.convolve(stf, self.trROT[0, :, i])[Nshift:self.npts+Nshift]
                self.trSYNROT[1, :, i]  = np.convolve(stf, self.trROT[1, :, i])[Nshift:self.npts+Nshift]
                self.trSYNROT[2, :, i]  = np.convolve(stf, self.trROT[2, :, i])[Nshift:self.npts+Nshift]
            else:
                self.trSYN[0, :, i]     = np.convolve(stf, self.trNEZ[0, :, i])[Nshift:self.npts+Nshift]
                self.trSYN[1, :, i]     = np.convolve(stf, self.trNEZ[1, :, i])[Nshift:self.npts+Nshift]
                self.trSYN[2, :, i]     = np.convolve(stf, self.trNEZ[2, :, i])[Nshift:self.npts+Nshift]
        self.stf    = stf
        return
    
    def rotcorr(self, trid, maxtime = 1.0, twin=1.5, dphi=.1):
        """
            shear wave splitting using the Rotation-Correlation method
            (e.g. Bowman and Ando,1987)
        ====================================================================================
        ::: input parameters :::
        trid        - trace index
        maxtime     - maximum delay time (sec)
        twin        - time length for selected window +/- twin (sec)
        dphi        - fast angle interval (deg)
        ::: output :::
        delayt      - delay time (sec)
        fa_phi      - fast axis angle (deg)
        Cmatrix     - correlation coefficient matrix
        ====================================================================================
        """
        maxlags = int(np.ceil(maxtime/self.dt)) # only +-4 seconds relevant
        zerolag = int(np.ceil(twin/self.dt))
        phi_test= (np.mgrid[-90:90-dphi:dphi])/180*np.pi
        ind     = (self.time > (self.shift-twin))*(self.time < (self.shift+twin))
        Cmatrix = np.zeros([int(phi_test.size), int(2*maxlags+1)])
        SG      = self.trSYNROT[0, :, trid] # R component
        SH      = self.trSYNROT[1, :, trid] # T component
        for iphi in xrange(phi_test.size):
            phi     = phi_test[iphi]
            # test slow-fast seismograms
            xphi    = SH*np.cos(phi) - SG* np.sin(phi)
            yphi    = SH*np.sin(phi) + SG* np.cos(phi)
            corTr   = np.correlate(a=xphi[ind], v=yphi[ind], mode='full')/(np.linalg.norm(xphi[ind])*np.linalg.norm(yphi[ind]))
            L       = int((corTr.size - 1)/2)
            Cmatrix[iphi, :] = corTr[L-maxlags:L+maxlags+1]
        indmax  = Cmatrix.argmax()
        pmax    = int(np.floor(indmax/(2*maxlags+1)))
        tmax    = indmax - pmax*(2*maxlags+1)
        delayt  = (tmax-maxlags)*self.dt
        fa_phi  = phi_test[pmax]/np.pi*180.
        if delayt < 0.:
            print 'neg delay time'
            delayt  = -delayt
            fa_phi  += 90.
            if fa_phi > 90.:
                fa_phi -= 180.
        return delayt, fa_phi, Cmatrix
    
    def rotcorr_st(self, maxtime = 1.0, twin=1.5, dphi=.1):
        """
        shear wave splitting using the Rotation-Correlation method for the stream (all the traces)
        ============================================================================================
        ::: input parameters :::
        maxtime     - maximum delay time (sec)
        twin        - time length for selected window +/- twin (sec)
        dphi        - fast angle interval (deg)
        ============================================================================================
        """
        self.delayt = np.array([])
        self.phi    = np.array([])
        Crc         = []
        for trid in xrange(self.ntr):
            delayt, phi, Cmatrix = self.rotcorr(trid=trid, maxtime = maxtime, twin=twin, dphi=dphi)
            self.delayt = np.append(self.delayt, delayt)
            self.phi    = np.append(self.phi, phi)
            Crc.append(Cmatrix)
        # fast axis angle
        self.fa_phi                     = self.bazArr + self.phi
        self.fa_phi[self.fa_phi<0.]     = self.fa_phi[self.fa_phi<0.] + 180.
        self.fa_phi[self.fa_phi>180.]   = self.fa_phi[self.fa_phi>180.] - 180.
        self.Crc    = np.array(Crc)
        return
    
    def eigenvalue(self, trid, maxtime = 1.0, twin=1.5, dphi=.5, mtype=1):
        """
            shear wave splitting using the eigenvalue method Silver & Chan (1991)
        ====================================================================================
        ::: input parameters :::
        trid        - trace index
        maxtime     - maximum delay time (sec)
        twin        - time length for selected window +/- twin (sec)
        dphi        - fast angle interval (deg)
        mtype       - method type
                        1 - max(lambda1 / lambda2)
                        2 - min(lambda2)
                        3 - max(lambda1)
                        4 - min(lambda1 * lambda2)
                        5 - min(lambda2 / lambda1)
        ::: output :::
        delayt      - delay time (sec)
        phi         - fast axis angle (deg)
        C           - eigenvalue matrix
        ====================================================================================
        """
        maxlags = np.ceil(maxtime/self.dt) # only +-4 seconds relevant
        zerolag = np.ceil(twin/self.dt)
        phi_test= (np.mgrid[-90:90-dphi:dphi])/180*np.pi
        dt_test = np.mgrid[0.:maxtime:self.dt]
        ind     = np.where((self.time > (self.shift-twin))*(self.time < (self.shift+twin)))[0]
        SG      = self.trSYNROT[0, :, trid] # R component
        SH      = self.trSYNROT[1, :, trid] # T component
        l1      = np.zeros([phi_test.size, dt_test.size]) # lambda 1
        l2      = np.zeros([phi_test.size, dt_test.size]) # lambda 2
        for iphi in xrange(phi_test.size):
            phi     = phi_test[iphi]
            xphi    = SH*np.cos(phi) - SG* np.sin(phi)
            yphi    = SH*np.sin(phi) + SG* np.cos(phi)
            tmpfast = yphi[ind]
            for idt  in xrange(dt_test.size):
                dt      = dt_test[idt]
                tmpslow = xphi[ind+idt]
                SHc     = tmpslow*np.cos(-phi) - tmpfast* np.sin(-phi)
                SGc     = tmpslow*np.sin(-phi) + tmpfast* np.cos(-phi)
                # SHe     = SHc - SHc.mean()
                # SGe     = SGc - SGc.mean()
                # construct covariance matrix
                covar   = np.cov(SHc, SGc)
                # eigenvalue estimation
                w, v    = np.linalg.eig(covar)
                l1[iphi, idt]   = w[1]
                l2[iphi, idt]   = w[0]
        # if mtype == 0:
        #     print 'minimum energy'
        #     C               = Ematrix
        #     ind0            = C.argmin()
        #     self.maplabel   = 'transverse energy'
        #     print ind0
        if mtype == 1:
            print 'Eigenvalue: max(lambda1 / lambda2)'
            C               = l1/l2
            ind0            = C.argmax()
            self.maplabel   = r'$\lambda_1/\lambda_2$'
            print ind0
        elif mtype == 2:
            print 'Eigenvalue: min(lambda2)'
            C               = l2
            ind0            = C.argmin()
            self.maplabel   = r'$\lambda_2$'
            print ind0
        elif mtype == 3:
            print 'Eigenvalue: max(lambda1)'
            C               = l1
            ind0            = C.argmax()
            self.maplabel   = r'$\lambda_1$'
            print ind0
        elif mtype == 4:
            print 'Eigenvalue: min(lambda1 * lambda2)'
            C               = l1*l2
            ind0            = C.argmin()
            self.maplabel   = r'$\lambda_1*\lambda_2$'
            print ind0
        elif mtype == 5:
            print 'Eigenvalue: min(lambda2 / lambda1)'
            C               = l2/l1
            ind0            = C.argmin()
            self.maplabel   = r'$\lambda_2/\lambda_1$'
            print ind0
        else:
            raise ValueError('Unexpected value of mtype')
        pind    = int(np.floor(ind0/(dt_test.size)))
        tind    = int(ind0 - pind*(dt_test.size))
        delayt  = dt_test[tind]
        phi     = phi_test[pind]/np.pi*180.
        return delayt, phi, C
    
    def ev_st(self, maxtime = 1.0, twin=1.5, dphi=.1, mtype=1):
        """
        shear wave splitting using the eigenvalue method for the stream (all the traces)
        ====================================================================================
        ::: input parameters :::
        maxtime     - maximum delay time (sec)
        twin        - time length for selected window +/- twin (sec)
        dphi        - fast angle interval (deg)
        mtype       - method type
                        1 - max(lambda1 / lambda2)
                        2 - min(lambda2)
                        3 - max(lambda1)
                        4 - min(lambda1 * lambda2)
                        5 - min(lambda2 / lambda1)
        ====================================================================================
        """
        self.delayt = np.array([])
        self.phi    = np.array([])
        Cev         = []
        for trid in xrange(self.ntr):
            delayt, phi, Cmatrix = self.eigenvalue(trid=trid, maxtime = maxtime, twin=twin, dphi=dphi, mtype=mtype)
            self.delayt = np.append(self.delayt, delayt)
            self.phi    = np.append(self.phi, phi)
            Cev.append(Cmatrix)
        # fast axis angle
        self.fa_phi                     = self.bazArr + self.phi
        self.fa_phi[self.fa_phi<0.]     = self.fa_phi[self.fa_phi<0.] + 180.
        self.fa_phi[self.fa_phi>180.]   = self.fa_phi[self.fa_phi>180.] - 180.
        self.Cev    = np.array(Cev)
        return
    
    def minimum_energy(self, trid, maxtime = 1.0, twin=1.5, dphi=.1):
        """
        shear wave splitting using the minimum energy method by Silver & Chan (1991)
        ====================================================================================
        ::: input parameters :::
        trid        - trace index
        maxtime     - maximum delay time (sec)
        twin        - time length for selected window +/- twin (sec)
        dphi        - fast angle interval (deg)
        ::: output :::
        delayt      - delay time (sec)
        phi         - fast axis angle (deg)
        Ematrix     - energy matrix
        ====================================================================================
        """
        maxlags = np.ceil(maxtime/self.dt) # only +-4 seconds relevant
        zerolag = np.ceil(twin/self.dt)
        phi_test= (np.mgrid[-90:90-dphi:dphi])/180*np.pi
        dt_test = np.mgrid[0.:maxtime:self.dt]
        ind     = np.where((self.time > (self.shift-twin))*(self.time < (self.shift+twin)))[0]
        SG      = self.trSYNROT[0, :, trid] # R component
        SH      = self.trSYNROT[1, :, trid] # T component
        Ematrix = np.zeros([phi_test.size, dt_test.size]) # energy matrix
        xphi    = np.outer(SH, np.cos(phi_test)) - np.outer(SG, np.sin(phi_test) )
        yphi    = np.outer(SH, np.sin(phi_test)) + np.outer(SG, np.cos(phi_test) )
        tmpfast = yphi[ind, :]
        for idt  in xrange(dt_test.size):
            dt      = dt_test[idt]
            tmpslow = xphi[ind+idt, :]
            SHc     = tmpslow*np.cos(-phi_test) - tmpfast * np.sin(-phi_test)
            SGc     = tmpslow*np.sin(-phi_test) + tmpfast* np.cos(-phi_test)
            Ematrix[:, idt] = np.linalg.norm(SHc, axis=0) # Energy on transverse component
        ind0            = Ematrix.argmin()
        self.maplabel   = 'transverse energy'
        pind    = int(np.floor(ind0/(dt_test.size)))
        tind    = int(ind0 - pind*(dt_test.size))
        delayt  = dt_test[tind]
        phi     = phi_test[pind]/np.pi*180.
        return delayt, phi, Ematrix
        
    def me_st(self, maxtime = 1.0, twin=1.5, dphi=.1):
        """
        shear wave splitting using the minimum energy method for the stream (all the traces)
        ====================================================================================
        ::: input parameters :::
        maxtime     - maximum delay time (sec)
        twin        - time length for selected window +/- twin (sec)
        dphi        - fast angle interval (deg)
        ====================================================================================
        """
        self.delayt = np.array([])
        self.phi    = np.array([])
        Cme         = []
        for trid in xrange(self.ntr):
            delayt, fa_phi, Cmatrix = self.minimum_energy(trid=trid, maxtime = maxtime, twin=twin, dphi=dphi)
            self.delayt = np.append(self.delayt, delayt)
            self.phi    = np.append(self.phi, fa_phi)
            Cme.append(Cmatrix)
        # fast axis angle
        self.fa_phi                     = self.bazArr + self.phi
        self.fa_phi[self.fa_phi<0.]     = self.fa_phi[self.fa_phi<0.] + 180.
        self.fa_phi[self.fa_phi>180.]   = self.fa_phi[self.fa_phi>180.] - 180.
        self.Cme    = np.array(Cme)
        return
    
    
