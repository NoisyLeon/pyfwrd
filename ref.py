# -*- coding: utf-8 -*-
"""
Module for synthetic receiver function computation, using anirec by Vadim Levin, and theo by Takuo Shibutani

The code is a python wrapper of the f77 code aniprop and theo

:Copyright:
    Author: Lili Feng
    Graduate Research Assistant
    CIEI, Department of Physics, University of Colorado Boulder
    email: lili.feng@colorado.edu
"""

import aniprop, theo
import numba
import numpy as np
import vmodel
import copy
import matplotlib.pyplot as plt

# define type of vmodel.model1d
model_type = numba.deferred_type()
model_type.define(vmodel.model1d.class_type.instance_type)



class ref_solver(object):
    
    """
    An object solving for receiver function using anirec by Vadim Levin, and theo by T. Shibutani
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
        self.azArr  = np.array([])
        return
    
    def init_default(self, dh=1., nl=100):
        self.dArr   = np.ones(nl, dtype = np.float32)*np.float32(dh)
        return
    
    def init_default_2(self):
        self.dArr   = np.array([20.,  15.,  42.,  43.,  45.,  35.], dtype = np.float32)
        return
    
    def solve_anirec(self, az=0., t=20., savestream=True):
        """
        Compute radial and transverse receiver function using anirec
        Default maximum velocity is 16.6667, can be changed by modifying the source code
        
        c     phase velocity of incident wave, important! LF 
            cc=16.6667  (line 175)
        ===================================================================================================================================
        ::: input parameters :::
        az          - azimuth of wave vector (az = 180. - baz, baz is the azimuth look at event from station)
                    NOTE: the azimuth of wave vector at station is typically NOT the same as the azimuth looking at station from event
        t           - time length of output in sec
        ::: output :::
        self.rfr    - radial receiver function
        self.rft    - transverse receiver function
        self.time   - time array
        ===================================================================================================================================
        """
        self.model.aniprop_check_model()
        z, rho, vp0, vp2, vp4, vs0, vs2 = self.model.layer_aniprop_model(self.dArr, 200, 1.)
        z       *= 1000.
        rho     *= 1000.
        vp0     *= 1000.
        vs0     *= 1000.
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
        baz     = 180. + az
        if baz > 360.:
            baz -= 360.
        # solve for receiver function using aniprop
        Rf,Tf,T     = aniprop.rf_aniso_interface(z,vp0,vp2,vp4,vs0,vs2,rho,theta,phig,nl,baz,ntimes)
        # radial component
        self.rfr    = Rf
        # transverse component
        self.rft    = Tf
        # time
        self.time   = T
        if savestream:
            self.rfrst.append(Rf); self.rftst.append(Tf)
            self.azArr  = np.append(self.azArr, az)
        return
    
    
    
    def solve_aniprop_reproduce(self, az=0., t=20.):
        """
        Compute radial and transverse receiver function using aniprop
        Default maximum velocity is 16.6667, can be changed by modifying the source code
        
        c     phase velocity of incident wave, important! LF 
            cc=16.6667  (line 175)
        ==============================================================================================
        ::: input parameters :::
        az          - azimuth
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
        # solve for receiver function using aniprop
        Rf,Tf,T     = aniprop.rf_aniso_interface(z,vp0,vp2,vp4,vs0,vs2,rho,theta,phig,nl,az,ntimes)
        # radial component
        self.rfr    = Rf
        # transverse component
        self.rft    = Tf
        # time
        self.time   = T
        self.rfrst.append(Rf); self.rftst.append(Tf)
        self.azArr  = np.append(self.azArr, az)
        return
    
    def solve_theo(self, t=20., slowness = 0.06, din = None):
        """
        Compute radial and transverse receiver function using theo
        ====================================================================================
        ::: input parameters :::
        t           - time length of output in sec
        slowness    - reference horizontal slowness (default - 0.06, 1./0.06=16.6667)
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
        hin = np.zeros(100, dtype=np.float32)
        vsin= np.zeros(100, dtype=np.float32)
        vpvs= np.zeros(100, dtype=np.float32)
        qsin= 1400.*np.ones(100, dtype=np.float32)
        qpin= 600.*np.ones(100, dtype=np.float32)
        if dArr.size<100:
            nl          = dArr.size
        else:
            nl  = 100
        hin[:nl]  = dArr[:nl]
        vsin[:nl] = (np.sqrt(LArr/rhoArr))[:nl]
        vpvs[:nl] = (np.sqrt(CArr/LArr))[:nl]
        fs      = 1./self.dt
        ntimes  = int(t*fs)
        # # nd      = 1000
        if din is None:
            din = 180.*np.arcsin(vsin[nl-1]*vpvs[nl-1]*slowness)/np.pi
        # solve for receiver function using theo
        rx 	= theo.theo(nl, vsin, hin, vpvs, qpin, qsin, fs, din, 2.5, 0.005, 0, ntimes)
        # receiver function (ONLY radial component)
        self.rf     = rx[:ntimes]
        self.time   = np.arange(ntimes, dtype=np.float32)*self.dt 
        return
    
    def compute_diff_ps_time(self, slowness = 0.06, h=35.):
        """Compute the difference in arrival time between P and P-S converted waves.
        """
        dArr, rhoArr, AArr, CArr, FArr, LArr, NArr = self.model.get_layer_model(self.dArr, 100, 1.)
        vs0     = (np.sqrt(LArr/rhoArr))[0]
        vp0     = (np.sqrt(CArr/rhoArr))[0]
        vp1     = (np.sqrt(CArr/rhoArr))[1]
        phis0   = np.arcsin(vs0*slowness)
        phip0   = np.arcsin(vp0*slowness)
        phip1   = np.arcsin(vp1*slowness)
        dtps    = h* (np.tan(phip0) - np.tan(phis0))*np.sin(phip1)/vp1 + h/np.cos(phis0)/vs0 - h/np.cos(phip0)/vp0
        print dtps
        
    def plot_az_rf(self, comp='T', showfig=True):
        ymax=361.
        ymin=-1.
        time    = self.time
        plt.figure()
        ax=plt.subplot()
        for i in xrange(self.azArr.size):
            if comp=='R':
                yvalue  = self.rfrst[i]
            else:
                yvalue  = self.rftst[i]
            rfmax   = yvalue.max()
            yvalue  = -yvalue/rfmax*40.
            yvalue[(time>1.8)*(time<3.2)] = 2.*yvalue[(time>1.8)*(time<3.2)]
            azi     = self.azArr[i]
            ax.plot(time, yvalue+azi, '-k', lw=0.3)
            ax.fill_between(time, y2=azi, y1=yvalue+azi, where=yvalue>0, color='red', lw=0.01, interpolate=True)
            ax.fill_between(time, y2=azi, y1=yvalue+azi, where=yvalue<0, color='blue', lw=0.01, interpolate=True)
            plt.axis([0., 8, ymin, ymax])
            plt.xlabel('Time(sec)')
            plt.title(comp+' component')
            plt.gca().invert_yaxis()
        if showfig: plt.show()
    
    
    