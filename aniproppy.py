# -*- coding: utf-8 -*-
"""
Module for surface wave dispersion computation using reflectivity method by Jeff Park

The code is a python wrapper of the f77 code aniprop.f by Jeff Park
Numba is used for speeding up of the code.

:Copyright:
    Author: Lili Feng
    Graduate Research Assistant
    CIEI, Department of Physics, University of Colorado Boulder
    email: lili.feng@colorado.edu
"""

import aniprop
import numba
import numpy as np
import vmodel
import copy

# define type of vmodel.model1d
model_type = numba.deferred_type()
model_type.define(vmodel.model1d.class_type.instance_type)



class aniprop_solver(object):
    
    """
    An object solving for dispersion curves using aniprop by Jeffrey Park
    =====================================================================================================================
    ::: Parameters :::
    model           - 1D Earth model object
    nmodes          - number of modes  
    egn96           - solve for eigenfunctions/kernels or not (default - True)
    cmin/cmax       - minimum/maximum velocity for root searching (default - -1.)
    T/freq          - period/frequency array
    dArr            - layer array (unit - km)
    ilvry           - 1 - Love wave, 2 - Rayleigh wave
    =====================================================================================================================
    """
    def __init__(self, inmodel):
        if not isinstance(inmodel, vmodel.model1d):
            raise ValueError('Input model should be type of vmodel.model1d !')
        self.model  = inmodel
        self.nmodes = 1
        self.verbose= 0
        return
    
    def init_default(self, dh=1., nl=200):
        self.Tmin   = 5.
        self.dT     = 5.
        self.Nt     = 20
        self.Tmax   = self.Tmin + (self.Nt-1)*self.dT
        self.T      = np.arange(self.Nt)*self.dT + self.Tmin
        self.dArr   = np.ones(nl, dtype = np.float32)*np.float32(dh)
        return
    
    def init_default_2(self):
        self.Tmin   = 5.
        self.dT     = 5.
        self.Nt     = 20
        self.Tmax   = self.Tmin + (self.Nt-1)*self.dT
        self.T      = np.arange(self.Nt)*self.dT + self.Tmin
        self.dArr   = np.array([20.,  15.,  42.,  43.,  45.,  35.], dtype = np.float32)
        return
    
    def init_default_3(self):
        self.Tmin   = 5.
        self.dT     = 5.
        self.Nt     = 1
        self.Tmax   = self.Tmin + (self.Nt-1)*self.dT
        self.T      = np.arange(self.Nt)*self.dT + self.Tmin
        self.dArr   = np.array([20.,  15.,  42.,  43.,  45.,  35.], dtype = np.float32)
        return
    
    def solve(self):
        z, rho, vp0, vp2, vp4, vs0, vs2 = self.model.layer_aniprop_model(self.dArr, 200, 1.)
        z       *= 1000.
        rho     *= 1000.
        vp0     *= 1000.
        vs0     *= 1000.
        
        self.z  = z; self.rho = rho; self.vp0 = vp0
        nl      = z.size - 1
        theta   = np.zeros(nl+1, dtype=np.float32)
        phig    = np.zeros(nl+1, dtype=np.float32)
        baz     = 0.
        Rphase,Rgroup,Lphase,Lgroup,Period = aniprop.aniprop_interface(z,vp0,vp2,vp4,vs0,vs2,rho,theta,phig,nl,baz, self.Nt, self.Tmin, self.Tmax)
        
        
        self.CR = Rphase/1000.
        self.UR = Rgroup/1000.
        self.CL = Lphase/1000.
        self.UL = Lgroup/1000.
        self.T  = Period
        return
        
    