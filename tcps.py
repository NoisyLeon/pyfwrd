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


spec = [('model',       model_type),
        ('dr',          numba.float32),
        ('r',           numba.float32[:]),
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
        ('omega',       numba.float32[:]),
        ('nmodes',      numba.int32)
        ]

@numba.jitclass(spec)
class tcps_solver(object):
    
    def __init__(self, inmodel):
        self.model  = inmodel
        
        self.nmodes = 1
        return
    
    