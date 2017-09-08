# -*- coding: utf-8 -*-
"""
Module for surface wave dispersion computation using normal mode method by Robert Herrmann

The code is a python wrapper of the f77 code Computer Program in Seismology 
Numba is used for speeding up of the code.

:Copyright:
    Author: Lili Feng
    Graduate Research Assistant
    CIEI, Department of Physics, University of Colorado Boulder
    email: lili.feng@colorado.edu
"""

import tdisp96
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

spec = [('model',       model_type),
        ('dr',          numba.float32),
        ('cmin',          numba.float32),
        ('cmax',          numba.float32),
        ('dArr',        numba.float32[:]),
        ('T',           numba.float32[:]),
        ('c',           numba.float32[:]),
        ('Vph',         numba.float32[:, :]),
        ('Vgr',         numba.float32[:, :]),
        ('eArr',        numba.int32[:, :]),
        ('r1data',      numba.float32[:, :, :]),
        ('r2data',      numba.float32[:, :, :]),
        ('r3data',      numba.float32[:, :, :]),
        ('r4data',      numba.float32[:, :, :]),
        ('l1data',      numba.float32[:, :, :]),
        ('l2data',      numba.float32[:, :, :]),
        ('Kadata',      numba.float32[:, :, :]),
        ('Kcdata',      numba.float32[:, :, :]),
        ('Kfdata',      numba.float32[:, :, :]),
        ('Kldata',      numba.float32[:, :, :]),
        ('Kndata',      numba.float32[:, :, :]),
        ('Krho0data',   numba.float32[:, :, :]),
        ('Kvphdata',    numba.float32[:, :, :]),
        ('Kvpvdata',    numba.float32[:, :, :]),
        ('Kvshdata',    numba.float32[:, :, :]),
        ('Kvsvdata',    numba.float32[:, :, :]),
        ('Ketadata',    numba.float32[:, :, :]),
        ('Krhodata',    numba.float32[:, :, :]),
        ('freq',       numba.float32[:]),
        ('nmodes',      numba.int32)
        ]

# @numba.jitclass(spec)
class tcps_solver(object):
    
    def __init__(self, inmodel):
        self.model  = inmodel
        self.nmodes = 1
        return
    
    def init_default(self):
        Tmin   = 10.
        Tmax   = 50.
        dT     = 5.
        self.cmin   = -1.
        self.cmax   = -1.
        self.T      = _get_array(Tmin, Tmax, dT)
        self.freq   = _value_divide_array(1., self.T)
        dh          = 1.
        nl          = 200
        self.dArr   = np.ones(nl, dtype = np.float32)*np.float32(dh)
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
    
    def init_output(self, wavetype):
        Nt              = self.T.size
        Nz              = self.r.size
        Vph             = np.zeros(self.nmodes*Nt, np.float32)
        self.Vph        = Vph.reshape(self.nmodes, Nt)
        self.Vgr        = self.Vph.copy()
        eArr            = np.zeros(self.nmodes*Nt, np.int32)
        self.eArr       = eArr.reshape(self.nmodes, Nt)
        tempdata        = np.zeros(self.nmodes*Nt*Nz, np.float32)
        tempdata        = tempdata.reshape(self.nmodes, Nt, Nz)
        if wavetype==1:
            self.r1data     = tempdata.copy()
            self.r2data     = tempdata.copy()
            self.r3data     = tempdata.copy()
            self.r4data     = tempdata.copy()
        else:
            self.l1data     = tempdata.copy()
            self.l2data     = tempdata.copy()
        self.Kadata     = tempdata.copy()
        self.Kcdata     = tempdata.copy()
        self.Kfdata     = tempdata.copy()
        self.Kldata     = tempdata.copy()
        self.Kndata     = tempdata.copy()
        self.Krhodata   = tempdata.copy()
        self.Krho0data  = tempdata.copy()
        self.Kvphdata   = tempdata.copy()
        self.Kvpvdata   = tempdata.copy()
        self.Kvshdata   = tempdata.copy()
        self.Kvsvdata   = tempdata.copy()
        self.Ketadata   = tempdata.copy()
        return
    
    def solve_PSV(self):
        #- root-finding algorithm ---------------------------------------------------------------------
        # self.init_output(1)
        dArr, rhoArr, AArr, CArr, FArr, LArr, NArr = self.model.get_cps_model(self.dArr, 200, 1.)
        verby = 1
        nfval = self.freq.size
        if self.model.flat == 0:
            iflsph_in = 2
        else:
            iflsph_in = 0
        ilvry = 2
        c_out,d_out,TA_out,TC_out,TF_out,TL_out,TN_out,TRho_out = tdisp96.disprs(ilvry, 1., nfval, 1, verby, nfval, np.append(self.freq, np.zeros(2049-nfval)), self.cmin,self.cmax,\
               dArr, AArr,CArr,FArr,LArr,NArr,rhoArr,\
               dArr.size, iflsph_in, 0., self.nmodes, 0.5, 0.5)
        
    
    