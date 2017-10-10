# -*- coding: utf-8 -*-
"""
Module for surface wave dispersion computation using fast_surf

The code is a python wrapper of the f77 code fast_surf

:Copyright:
    Author: Lili Feng
    Graduate Research Assistant
    CIEI, Department of Physics, University of Colorado Boulder
    email: lili.feng@colorado.edu
"""

import fast_surf
import numba
import numpy as np
import vmodel
import copy

# define type of vmodel.model1d
model_type = numba.deferred_type()
model_type.define(vmodel.model1d.class_type.instance_type)


#--------------------------------------------------------------------------------------------------
#- fundamental array manipulations
#--------------------------------------------------------------------------------------------------
@numba.jit(numba.float32[:](numba.float32, numba.float32, numba.float32))
def _get_array(xmin, xmax, dx):
    xlst= []
    Nx  = int((xmax - xmin)/dx + 1)
    for i in xrange(Nx): xlst.append(dx*i+xmin)
    return np.array(xlst, dtype=np.float32)

@numba.jit(numba.float32[:](numba.float32, numba.float32[:]))
def _value_divide_array(value, array):
    outArr  = np.zeros(array.size, dtype=np.float32)
    for i in xrange(array.size): outArr[i] = value/array[i]
    return outArr

@numba.jit(numba.float32[:](numba.float32, numba.float32[:]))
def _array_divide_value(value, array):
    outArr  = np.zeros(array.size, dtype=np.float32)
    for i in xrange(array.size): outArr[i] = array[i]/value
    return outArr

@numba.jit(numba.float32[:](numba.float32[:], numba.float32[:]))
def _merge_array(a1, a2):
    a3  = np.zeros(a1.size+a2.size, dtype=np.float32)
    Na1 = a1.size
    for i in xrange(a3.size):
        if i <= Na1-1:
            a3[i] = a1[i]
        else:
            a3[i] = a2[i-Na1]
    return a3

@numba.jit(numba.float32(numba.float32[:]))
def _abs_max_(array):
    mvalue=np.abs(array[0])
    for i in xrange(array.size):
        if np.abs(array[i]) > mvalue: mvalue = np.abs(array[i])
    return mvalue


class fsurf_solver(object):
    """
    An object solving for dispersion curves using fast_surf, faster than tcps_solver
    NOTE: Earth flattening transformation is included in the fast_surf
    =================================================================================
    ::: parameters :::
    model           - 1D Earth model object (should be isotropic)
    T/freq          - period/frequency array
    dArr            - layer array (unit - km)
    =================================================================================
    """
    def __init__(self, inmodel):
        if not isinstance(inmodel, vmodel.model1d):
            raise ValueError('Input model should be type of vmodel.model1d !')
        self.model  = inmodel
        if not self.model.is_iso():
            print 'WARNING: model is anisotropic!'
        return
    
    def init_default(self, dh=1., nl=200):
        Tmin        = 5.
        Tmax        = 100.
        dT          = 5.
        self.T      = _get_array(Tmin, Tmax, dT)
        self.dArr   = np.ones(nl, dtype = np.float32)*np.float32(dh)
        return
    
    def init_default_3(self, dh=1., nl=200):
        Tmin        = 5.
        Tmax        = 5.
        dT          = 5.
        self.T      = _get_array(Tmin, Tmax, dT)
        self.dArr   = np.ones(nl, dtype = np.float32)*np.float32(dh)
        return
    
    def init_default_2(self):
        Tmin        = 5.
        Tmax        = 100.
        dT          = 5.
        self.T      = _get_array(Tmin, Tmax, dT)
        self.dArr   = np.array([20.,  15.,  42.,  43.,  45.,  35.], dtype = np.float32)
        return
    
    # def init_dbase(self, T, c, rmin, dr, nmodes):
    #     
    #     self.T      = T
    #     self.c      = c
    #     self.omega  = _value_divide_array(np.float32(2.*np.pi), self.T)
    #     self.dr     = dr
    #     self.r      = _get_array(rmin, 6371000., self.dr)
    #     self.nmodes = nmodes
    #     return
    
    def solve_surf(self, ilvry=2, qs=600., checkiso=True):
        """
        Solve for P-SV motion dispersion curves
        ===================================================================================================
        ::: input parameters :::
        ilvry       - Love (1) or Rayleigh (2) wave computation
        qs          - S-wave quality factor (default - 600.)
        ::: output :::
        : model :
        dArr, vsArr, vpArr, rhoArr  - layerized model
        : dispersion :
        CR, UR      - Rayleigh wave phase/group velocity
        CL, UL      - Love wave phase/group velocity
        ===================================================================================================
        """
        try:
            vsin                = self.vsArr
            vpin                = self.vpArr
            rhoin               = self.rhoArr
            dArr                = self.dArr
        except AttributeError:
            dArr, rhoArr, AArr, CArr, FArr, LArr, NArr = self.model.get_layer_model(self.dArr, 200, 1.)
            vsin                = np.sqrt(LArr/rhoArr)
            vpin                = np.sqrt(CArr/rhoArr)
            rhoin               = rhoArr
            self.vsArr          = vsin
            self.vpArr          = vpin
            self.rhoArr         = rhoin
            self.dArr           = dArr
        qsinv               = 1./(np.ones(vsin.size, dtype=np.float32)*qs)
        nper                = self.T.size
        per                 = np.zeros(200, dtype=np.float32)
        per[:nper]          = self.T[:]
        if ilvry != 1 and ilvry !=2:
            raise ValueError('Unexpected ilvry value')
        #- root-finding algorithm using fast_surf, compute phase velocities
        (ur0,ul0,cr0,cl0)   = fast_surf.fast_surf(vsin.size, ilvry, vpin, vsin, rhoin, dArr, qsinv, per, nper)
        if ilvry == 2:
            self.CR = cr0[:nper]; self.UR  = ur0[:nper]
        else:
            self.CL = cl0[:nper]; self.UL  = ul0[:nper]
        return
    
    
    