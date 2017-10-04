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
    T               - period array
    dArr            - layer array (unit - km)
    CR/UR           - Rayleigh wave phase/group velocity (unit - km/sec)
    CL/UL           - Love wave phase/group velocity (unit - km/sec)
    =====================================================================================================================
    """
    def __init__(self, inmodel):
        if not isinstance(inmodel, vmodel.model1d):
            raise ValueError('Input model should be type of vmodel.model1d !')
        self.model  = inmodel
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
        self.Tmin   = 10.
        self.dT     = 5.
        self.Nt     = 20
        self.Tmax   = self.Tmin + (self.Nt-1)*self.dT
        self.T      = np.arange(self.Nt)*self.dT + self.Tmin
        self.dArr   = np.array([20.,  15.,  42.,  43.,  45.,  35.], dtype = np.float32)
        return
    
    def solve_surf(self, baz):
        """
        Solve for Rayleigh and Love dispersion curves using reflectivity method
        """
        self.model.aniprop_check_model()
        if self.model.flat == 1:
            z, rho, vp0, vp2, vp4, vs0, vs2 = self.model.layer_aniprop_model(self.dArr, 200, 1.)
            self.dip, self.strike = self.model.angles_aniprop_model(z)
            self.z = z
            self.rho=rho
            self.vp0=vp0
            self.vp2=vp2
            self.vp4=vp4
            self.vs0=vs0
            
            # n = z.size
            # vs2= np.ones(n)*-0.2
            self.vs2=vs2
            
        else:
            zl, rhol, vp0l, vp2l, vp4l, vs0l, vs2l = vmodel.layer_aniprop_model_sph(inmodel = self.model, dArr = self.dArr, nl = 200, dh = 1., ilvry=1)
            zr, rhor, vp0r, vp2r, vp4r, vs0r, vs2r = vmodel.layer_aniprop_model_sph(inmodel = self.model, dArr = self.dArr, nl = 200, dh = 1., ilvry=2)
            self.zl = zl
            self.rhol=rhol
            self.vp0l=vp0l
            self.vp2l=vp2l
            self.vp4l=vp4l
            self.vs0l=vs0l
            self.vs2l=vs2l
            
            self.zr = zr
            self.rhor=rhor
            self.vp0r=vp0r
            self.vp2r=vp2r
            self.vp4r=vp4r
            self.vs0r=vs0r
            self.vs2r=vs2r
        if self.model.flat == 1:
            z       *= 1000.
            rho     *= 1000.
            vp0     *= 1000.
            vs0     *= 1000.
            ##########################################
            nl      = z.size - 1
            theta   = self.dip
            phig    = np.zeros(nl+1, dtype=np.float32)
            phig[self.dip>0.]  = self.strike[self.dip>0.] + 270.
            # phig[phig>=360.] = phig[phig>=360.] - 360.
            # theta   = np.zeros(nl+1, dtype=np.float32)
            # phig    = np.zeros(nl+1, dtype=np.float32)
            # baz     = 0.
            ###########################################
            print phig
            Rphase,Rgroup,Lphase,Lgroup,Period = aniprop.aniprop_interface(z,vp0,vp2,vp4,vs0,vs2,rho,theta,phig,nl,baz, self.Nt, self.Tmin, self.Tmax)
        else:
            zl      *= 1000.
            rhol    *= 1000.
            vp0l    *= 1000.
            vs0l    *= 1000.
            zr      *= 1000.
            rhor    *= 1000.
            vp0r    *= 1000.
            vs0r    *= 1000.
            ##########################################
            nl      = zl.size - 1
            theta   = np.zeros(nl+1, dtype=np.float32)
            phig    = np.zeros(nl+1, dtype=np.float32)
            baz     = 0.
            ###########################################
            Rphase0,Rgroup0,Lphase,Lgroup,Period = aniprop.aniprop_interface(zl,vp0l,vp2l,vp4l,vs0l,vs2l,rhol,\
                                                theta,phig,nl,baz, self.Nt, self.Tmin, self.Tmax)
            Rphase,Rgroup,Lphase0,Lgroup0,Period = aniprop.aniprop_interface(zr,vp0r,vp2r,vp4r,vs0r,vs2r,rhor,\
                                                theta,phig,nl,baz, self.Nt, self.Tmin, self.Tmax)
        self.CR = Rphase/1000.
        self.UR = Rgroup/1000.
        self.CL = Lphase/1000.
        self.UL = Lgroup/1000.
        self.T  = Period
        return
        
    