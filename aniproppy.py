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
    
    def __init__(self, inmodel):
        self.model  = inmodel
        self.dr     = 1000.
        self.nmodes = 1
        return