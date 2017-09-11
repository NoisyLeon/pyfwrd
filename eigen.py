# -*- coding: utf-8 -*-
"""
Module for surface wave eigenfunction and dispersion computation

The code is based on the surf package by Andreas Fichtner.
Numba is used for speeding up of the code.

:Copyright:
    Author: Lili Feng
    Graduate Research Assistant
    CIEI, Department of Physics, University of Colorado Boulder
    email: lili.feng@colorado.edu
"""

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

#--------------------------------------------------------------------------------------------------
#- fundamental functions for integrate_psv_alt
#--------------------------------------------------------------------------------------------------
@numba.jit(numba.float32(numba.float32, numba.float32, numba.float32, numba.float32))
def f1_psv_alt(C,L,r4,r5): return (r4 / L - r5 / C)

@numba.jit(numba.float32(numba.float32, numba.float32, numba.float32, numba.float32, numba.float32,
        numba.float32, numba.float32, numba.float32))
def f2_psv_alt(rho,A,C,F,omega,k,r4,r5): return (-omega**2 * rho * r4 + (omega**2 * rho - k**2 * (A - F**2 / C)) * r5)

@numba.jit(numba.float32(numba.float32, numba.float32, numba.float32, numba.float32, numba.float32))
def f3_psv_alt(C,F,k,r4,r5): return (k * r4 + k * F * r5 / C)

@numba.jit(numba.float32(numba.float32, numba.float32, numba.float32, numba.float32, numba.float32,
        numba.float32, numba.float32, numba.float32, numba.float32))
def f4_psv_alt(rho,A,C,F,omega,k,r1,r2,r3): return ((-omega**2 * rho + k**2 * (A - F**2 / C)) * r1 + r2 / C - 2 * k * F * r3 / C)

@numba.jit(numba.float32(numba.float32, numba.float32, numba.float32, numba.float32, numba.float32,
        numba.float32, numba.float32))
def f5_psv_alt(rho,L,omega,k,r1,r2,r3):
	return (omega**2 * rho * r1 - r2 / L - 2 * k * r3)

#--------------------------------------------------------------------------------------------------
#- numerical integration
#--------------------------------------------------------------------------------------------------

@numba.njit( numba.types.UniTuple(numba.float32[:], 5) (model_type, numba.float32[:], numba.float32, numba.float32, numba.float32) )
def integrate_psv_alt(model, r, dr, omega, k):
    """
    Integrate first-order P-SV system for a fixed angular frequency omega and a fixed wavenumber k.
    =================================================================================================
    Input Parameters:
    model   - input model1d object
    r       - radius array in m
    dr      - radius increment in m
    omege   - angular frequency in Hz
    k       - wave number in 1/m
    -------------------------------------------------------------------------------------------------
    Output:
    r1, ... -	variables of the alternative Rayleigh wave system
    =================================================================================================
    """
    #- initialization -----------------------------------------------------------------------------
    r1  = np.zeros(r.size, dtype=np.float32)
    r2  = np.zeros(r.size, dtype=np.float32)
    r3  = np.zeros(r.size, dtype=np.float32)
    r4  = np.zeros(r.size, dtype=np.float32)
    r5  = np.zeros(r.size, dtype=np.float32)
    rho, A, C, F, L, N = model.get_ind_Love_parameters_PSV(0)

    #- check if phase velocity is below S velocity ------------------------------------------------
    if (k**2 - (omega**2 * rho / L)) > 0.0:  
        #- set initial values
        r1[0] = 0.0	
        r2[0] = 1.0 
        r3[0] = 0.0
        r4[0] = 0.0
        r5[0] = 0.0
        #- integrate upwards with 4th-order Runge-Kutta--------------------------------------------
        for n in xrange(r.size-1):
            # # # dr = drArr[n]
            #- compute Runge-Kutta coeficients for l1 and l2
            rho, A, C, F, L, N = model.get_r_love_parameters_PSV(r[n])
            
            K1_1 = f1_psv_alt(C,L,r4[n],r5[n])
            K2_1 = f2_psv_alt(rho,A,C,F,omega,k,r4[n],r5[n])
            K3_1 = f3_psv_alt(C,F,k,r4[n],r5[n])
            K4_1 = f4_psv_alt(rho,A,C,F,omega,k,r1[n],r2[n],r3[n])
            K5_1 = f5_psv_alt(rho,L,omega,k,r1[n],r2[n],r3[n]) 
            
            rho, A, C, F, L, N = model.get_r_love_parameters_PSV(r[n]+dr/2.0)
            
            K1_2 = f1_psv_alt(C,L,r4[n]+0.5*K4_1*dr,r5[n]+0.5*K5_1*dr)
            K2_2 = f2_psv_alt(rho,A,C,F,omega,k,r4[n]+0.5*K4_1*dr,r5[n]+0.5*K5_1*dr)
            K3_2 = f3_psv_alt(C,F,k,r4[n]+0.5*K4_1*dr,r5[n]+0.5*K5_1*dr)
            K4_2 = f4_psv_alt(rho,A,C,F,omega,k,r1[n]+0.5*K1_1*dr,r2[n]+0.5*K2_1*dr,r3[n]+0.5*K3_1*dr)
            K5_2 = f5_psv_alt(rho,L,omega,k,r1[n]+0.5*K1_1*dr,r2[n]+0.5*K2_1*dr,r3[n]+0.5*K3_1*dr) 
            
            K1_3 = f1_psv_alt(C,L,r4[n]+0.5*K4_2*dr,r5[n]+0.5*K5_2*dr)
            K2_3 = f2_psv_alt(rho,A,C,F,omega,k,r4[n]+0.5*K4_2*dr,r5[n]+0.5*K5_2*dr)
            K3_3 = f3_psv_alt(C,F,k,r4[n]+0.5*K4_2*dr,r5[n]+0.5*K5_2*dr)
            K4_3 = f4_psv_alt(rho,A,C,F,omega,k,r1[n]+0.5*K1_2*dr,r2[n]+0.5*K2_2*dr,r3[n]+0.5*K3_2*dr)
            K5_3 = f5_psv_alt(rho,L,omega,k,r1[n]+0.5*K1_2*dr,r2[n]+0.5*K2_2*dr,r3[n]+0.5*K3_2*dr)
            
            rho, A, C, F, L, N = model.get_r_love_parameters_PSV(r[n]+dr)
            
            K1_4 = f1_psv_alt(C,L,r4[n]+K4_3*dr,r5[n]+K5_3*dr)
            K2_4 = f2_psv_alt(rho,A,C,F,omega,k,r4[n]+K4_3*dr,r5[n]+K5_3*dr)
            K3_4 = f3_psv_alt(C,F,k,r4[n]+K4_3*dr,r5[n]+K5_3*dr)
            K4_4 = f4_psv_alt(rho,A,C,F,omega,k,r1[n]+K1_3*dr,r2[n]+K2_3*dr,r3[n]+K3_3*dr)
            K5_4 = f5_psv_alt(rho,L,omega,k,r1[n]+K1_3*dr,r2[n]+K2_3*dr,r3[n]+K3_3*dr) 
            
            #- update
            r1[n + 1] = r1[n] + dr * (K1_1 + 2 * K1_2 + 2 * K1_3 + K1_4) / 6.0
            r2[n + 1] = r2[n] + dr * (K2_1 + 2 * K2_2 + 2 * K2_3 + K2_4) / 6.0
            r3[n + 1] = r3[n] + dr * (K3_1 + 2 * K3_2 + 2 * K3_3 + K3_4) / 6.0
            r4[n + 1] = r4[n] + dr * (K4_1 + 2 * K4_2 + 2 * K4_3 + K4_4) / 6.0
            r5[n + 1] = r5[n] + dr * (K5_1 + 2 * K5_2 + 2 * K5_3 + K5_4) / 6.0
            #- rescale to maximum to prevent overflow
            mm  = np.max(np.abs(r2))
            r1  = r1 / mm
            r2  = r2 / mm
            r3  = r3 / mm
            r4  = r4 / mm
            r5  = r5 / mm
    return r1, r2, r3, r4, r5

#--------------------------------------------------------------------------------------------------
#- fundamental functions for integrate_psv
#--------------------------------------------------------------------------------------------------
@numba.jit(numba.float32 (numba.float32, numba.float32, numba.float32, numba.float32, numba.float32) )
def f1_psv(C,F,k,r2,r3):
	return (r2 / C + k * F * r3 / C)

@numba.jit(numba.float32 (numba.float32, numba.float32, numba.float32, numba.float32, numba.float32) )
def f2_psv(rho,omega,k,r1,r4):
	return (-rho * omega**2 * r1 + k * r4)

@numba.jit(numba.float32 (numba.float32, numba.float32, numba.float32, numba.float32) )
def f3_psv(L,k,r1,r4):
	return (r4 / L - k * r1)

@numba.jit(numba.float32 (numba.float32, numba.float32, numba.float32, numba.float32, numba.float32,
         numba.float32, numba.float32, numba.float32) )
def f4_psv(rho,A,C,F,omega,k,r2,r3):
	return (-k * F * r2 / C + (k**2 * (A - F**2 / C) - rho * omega**2) * r3)

# 
@numba.jit( numba.types.UniTuple(numba.float32[:], 4) (model_type, numba.float32[:], numba.float32, numba.float32, numba.float32, numba.int32) )
def integrate_psv(model, r, dr, omega, k, initial_condition):
    """
    Integrate first-order P-SV system for a fixed angular frequency omega and a fixed wavenumber k.
    =================================================================================================
    Input Parameters:
    model               - input model1d object
    r                   - radius array in m
    dr                  - radius increment in m
    omege               - angular frequency in Hz
    k                   - wave number in 1/m
    initial_condition   - initial condition for eigenfunctions
    -------------------------------------------------------------------------------------------------
    Output:
    r1, ... -	variables of the alternative Rayleigh wave system
    =================================================================================================
    """
    #- initialization -----------------------------------------------------------------------------
    r1  = np.zeros(r.size, dtype=np.float32)
    r2  = np.zeros(r.size, dtype=np.float32)
    r3  = np.zeros(r.size, dtype=np.float32)
    r4  = np.zeros(r.size, dtype=np.float32)
    rho, A, C, F, L, N = model.get_ind_Love_parameters_PSV(0)
    
    #- check if phase velocity is below S velocity ------------------------------------------------
    if (1==1): #(k**2 - (omega**2 * rho / L)) > 0.0:
        #- set initial values
        if initial_condition == 1:
            r1[0] = 0.0	
            r2[0] = 1.0 
            r3[0] = 0.0
            r4[0] = 0.0
        else:
            r1[0] = 0.0	
            r2[0] = 0.0 
            r3[0] = 0.0
            r4[0] = 1.0
        #- integrate upwards with 4th-order Runge-Kutta--------------------------------------------
        for n in xrange(r.size-1):
            #- compute Runge-Kutta coeficients for l1 and l2
            rho, A, C, F, L, N = model.get_r_love_parameters_PSV(r[n])
            K1_1 = f1_psv(C,F,k,r2[n],r3[n])
            K2_1 = f2_psv(rho,omega,k,r1[n],r4[n])
            K3_1 = f3_psv(L,k,r1[n],r4[n])
            K4_1 = f4_psv(rho,A,C,F,omega,k,r2[n],r3[n])
    
            rho, A, C, F, L, N = model.get_r_love_parameters_PSV(r[n] + dr/2.0)
            K1_2 = f1_psv(C,F,k,r2[n]+0.5*K2_1*dr,r3[n]+0.5*K3_1*dr)
            K2_2 = f2_psv(rho,omega,k,r1[n]+0.5*K1_1*dr,r4[n]+0.5*K4_1*dr)
            K3_2 = f3_psv(L,k,r1[n]+0.5*K1_1*dr,r4[n]+0.5*K4_1*dr)
            K4_2 = f4_psv(rho,A,C,F,omega,k,r2[n]+0.5*K2_1*dr,r3[n]+0.5*K3_1*dr)
    
            K1_3 = f1_psv(C,F,k,r2[n]+0.5*K2_2*dr,r3[n]+0.5*K3_2*dr)
            K2_3 = f2_psv(rho,omega,k,r1[n]+0.5*K1_2*dr,r4[n]+0.5*K4_2*dr)
            K3_3 = f3_psv(L,k,r1[n]+0.5*K1_2*dr,r4[n]+0.5*K4_2*dr)
            K4_3 = f4_psv(rho,A,C,F,omega,k,r2[n]+0.5*K2_2*dr,r3[n]+0.5*K3_2*dr)
            
            rho, A, C, F, L, N = model.get_r_love_parameters_PSV(r[n] + dr)
            K1_4 = f1_psv(C,F,k,r2[n]+K2_3*dr,r3[n]+K3_3*dr)
            K2_4 = f2_psv(rho,omega,k,r1[n]+K1_3*dr,r4[n]+K4_3*dr)
            K3_4 = f3_psv(L,k,r1[n]+K1_3*dr,r4[n]+K4_3*dr)
            K4_4 = f4_psv(rho,A,C,F,omega,k,r2[n]+K2_3*dr,r3[n]+K3_3*dr)
    
            #- update
            r1[n + 1] = r1[n] + dr * (K1_1 + 2 * K1_2 + 2 * K1_3 + K1_4) / 6.0
            r2[n + 1] = r2[n] + dr * (K2_1 + 2 * K2_2 + 2 * K2_3 + K2_4) / 6.0
            r3[n + 1] = r3[n] + dr * (K3_1 + 2 * K3_2 + 2 * K3_3 + K3_4) / 6.0
            r4[n + 1] = r4[n] + dr * (K4_1 + 2 * K4_2 + 2 * K4_3 + K4_4) / 6.0
    
            #- rescale to maximum to prevent overflow
            if initial_condition == 1:
                mm = np.max(np.abs(r2))
            else:
                mm = np.max(np.abs(r4))
            r1 = r1 / mm
            r2 = r2 / mm
            r3 = r3 / mm
            r4 = r4 / mm
    return r1, r2, r3, r4

@numba.jit( numba.types.UniTuple(numba.float32, 3) (numba.float32[:], numba.float32[:],numba.float32[:], numba.float32[:],\
    numba.float32[:], numba.float32, numba.float32, numba.float32[:], numba.float32[:], numba.float32[:], numba.float32[:], \
    numba.float32[:], numba.float32[:]) )
def group_velocity_psv(r1, r2, r3, r4, r, k, phase_velocity, rho, A, C, F, L, N):
    """
    Compute Rayleigh wave group velocity for a given set of eigenfunctions (l1, l2), and
    phase velocity (phase_velocity). 
    
    U, I1, I3 = group_velocity_psv(r1, r2, r3, r4, r, phase_velocity, rho, A, C, F, L, N)
    """
    
    #- compute the integrals I1 and I3 ------------------------------------------------------------
    I1 = 0.0
    I3 = 0.0
    dr = r[1] - r[0]
    for n in xrange(r.size-1):
        I1 = I1 + rho[n] * (r1[n]**2 + r3[n]**2) 
        I3 = I3 + ((A[n]-F[n]**2/C[n])*r3[n]**2 + (r1[n]*r4[n]-F[n]*r2[n]*r3[n]/C[n])/k) 
    I1 = dr * I1
    I3 = dr * I3	
    #- return values ------------------------------------------------------------------------------
    U = I3 / (phase_velocity * I1)
    return U, I1, I3

@numba.jit( numba.types.UniTuple(numba.float32[:], 12) (numba.float32[:], numba.float32[:],numba.float32[:], numba.float32[:],\
    numba.float32[:], numba.float32, numba.float32, numba.float32, numba.float32[:], numba.float32[:], numba.float32[:], numba.float32[:], \
    numba.float32[:], numba.float32[:]) )
def kernels_psv(r, r1, r2, r3, r4, omega, k, I3, rho, A, C, F, L, N):
    """
    Compute and write sensitivity kernels for Love waves.
    
    kernels_psv(r, r1, r2, r3, r4, _omega, _k, I3, rho, A, C, F, L, N)
    """    
    #- initialisations ----------------------------------------------------------------------------
    
    K_rho_0 = np.zeros(r.size, dtype=np.float32)
    K_A_0   = np.zeros(r.size, dtype=np.float32)
    K_C_0   = np.zeros(r.size, dtype=np.float32)
    K_F_0   = np.zeros(r.size, dtype=np.float32)
    K_L_0   = np.zeros(r.size, dtype=np.float32)
    K_N_0   = np.zeros(r.size, dtype=np.float32)
    
    K_rho   = np.zeros(r.size, dtype=np.float32)
    K_vph   = np.zeros(r.size, dtype=np.float32)
    K_vpv   = np.zeros(r.size, dtype=np.float32)
    K_vsh   = np.zeros(r.size, dtype=np.float32)
    K_vsv   = np.zeros(r.size, dtype=np.float32)
    K_eta   = np.zeros(r.size, dtype=np.float32)
    #- phase velocity
    v = omega / k
    #- velocities and eta
    vpv = np.sqrt(A/rho)
    vph = np.sqrt(C/rho)
    vsh = np.sqrt(N/rho)
    vsv = np.sqrt(L/rho)
    eta = F / (A-np.float32(2.)*L)
    #- compute fundamental kernels ----------------------------------------------------------------
    if (I3!=0.0):
        K_rho_0 = -v**2 * omega * (r1**2 + r3**2)
        K_A_0 = omega * r3**2
        K_C_0 = v * (r2 + k*F*r3)**2 / (k*C**2)
        K_F_0 = -np.float32(2.)*v * (r2 + k*(F*r3)) * r3/C
        K_L_0 = (v/k) * (r4/L)**2
        K_N_0 = np.zeros(r.size, dtype=np.float32)
        
        K_rho_0 = K_rho_0 / (np.float32(2.)*k*I3)
        K_A_0 = K_A_0 / (np.float32(2.)*k*I3)
        K_C_0 = K_C_0 / (np.float32(2.)*k*I3)
        K_F_0 = K_F_0 / (np.float32(2.)*k*I3)
        K_L_0 = K_L_0 / (np.float32(2.)*k*I3)
        K_N_0 = K_N_0 / (np.float32(2.)*k*I3)
    #- compute relative kernels in velocity parametrisation ---------------------------------------
    K_vph   = np.float32(2.)*A*K_A_0 + np.float32(2.)*A*eta*K_F_0
    # return K_rho_0, K_A_0, K_C_0, K_F_0, K_L_0, K_N_0, K_rho, K_vph, K_vpv, K_vsh, K_vsv, K_eta
    K_vpv   = np.float32(2.)*C*K_C_0
    K_vsh   = np.float32(2.)*N*K_N_0
    K_vsv   = np.float32(2.)*L*K_L_0 - np.float32(4.)*L*eta*K_F_0
    K_rho   = rho*K_rho_0 + A*K_A_0 + C*K_C_0 + N*K_N_0 + L*K_L_0 + (A-np.float32(2.)*L)*eta*K_F_0 
    K_eta   = F*K_F_0
    #- convert to relative kernels ----------------------------------------------------------------
    K_rho_0 = rho*K_rho_0
    K_A_0   = A*K_A_0
    K_C_0   = C*K_C_0
    K_F_0   = F*K_F_0
    K_L_0   = L*K_L_0
    K_N_0   = N*K_N_0
    return K_rho_0, K_A_0, K_C_0, K_F_0, K_L_0, K_N_0, K_rho, K_vph, K_vpv, K_vsh, K_vsv, K_eta
    

###################################################################################
# functions for computation of SH wave
###################################################################################

#--------------------------------------------------------------------------------------------------
#- displacement function
#--------------------------------------------------------------------------------------------------
@numba.jit(numba.float32 (numba.float32, numba.float32) )
def f1_sh(L,l2):
    """
    Right-hand side of the displacement derivative for SH propagation.
    """
    return l2 / L
#--------------------------------------------------------------------------------------------------
#- stress function
#--------------------------------------------------------------------------------------------------
@numba.jit(numba.float32 (numba.float32, numba.float32, numba.float32, numba.float32, numba.float32) )
def f2_sh(N, rho, k, omega, l1):
    """
    Right-hand side of the stress derivative for SH propagation.
    """
    return (N * k**2 - rho * omega**2) * l1


#--------------------------------------------------------------------------------------------------
#- numerical integration
#--------------------------------------------------------------------------------------------------
@numba.jit( numba.types.UniTuple(numba.float32[:], 2) (model_type, numba.float32[:], numba.float32, numba.float32, numba.float32) )
def integrate_sh(model, r, dr, omega, k):
	"""
	Integrate first-order system for a fixed circular frequency omega and a fixed wavenumber k.
	l1, l2, r = integrate_sh(r_min, dr, omega, k, model)

	r_min:		minimum radius in m
	dr:			radius increment in m
	omega:		circular frequency in Hz
	k:			wave number in 1/m
	model:		Earth model, e.g. "PREM", "GUTENBERG", ... .

	l1, l2:		variables of the Love wave system
	r:			radius vector in m
	"""

	#- initialisation -----------------------------------------------------------------------------

	# r = np.arange(r_min, 6371000.0 + dr, dr, dtype=float)
	l1  = np.zeros(r.size, dtype=np.float32)
	l2  = np.zeros(r.size, dtype=np.float32)
	rho, A, C, F, L, N = model.get_ind_Love_parameters_SH(0)
	#- check if phase velocity is below S velocity ------------------------------------------------
	if (1==1): #(k**2 - (omega**2 * rho / L)) > 0.0:
		#- set initial values
		l2[0] = 1.0 
		l1[0] = 0.0

		#- integrate upwards with 4th-order Runge-Kutta--------------------------------------------

		for n in xrange(r.size-1):
			#- compute Runge-Kutta coeficients for l1 and l2
			rho, A, C, F, L, N = model.get_r_love_parameters_SH(r[n])
			K1_1 = f1_sh(L, l2[n])
			K2_1 = f2_sh(N, rho, k, omega, l1[n]) 

			rho, A, C, F, L, N = model.get_r_love_parameters_SH(r[n]+dr/2.0)
			K1_2 = f1_sh(L, l2[n] + K2_1 * dr / 2.0)
			K2_2 = f2_sh(N, rho, k, omega, l1[n] + K1_1 * dr / 2.0)

			K1_3 = f1_sh(L, l2[n] + K2_2 * dr / 2.0)
			K2_3 = f2_sh(N, rho, k, omega, l1[n] + K1_2 * dr / 2.0)

			rho, A, C, F, L, N = model.get_r_love_parameters_SH(r[n]+dr)
			K1_4 = f1_sh(L, l2[n] + K2_3 * dr)
			K2_4 = f2_sh(N, rho, k, omega, l1[n] + K1_3 * dr)

			#- update
			l1[n + 1] = l1[n] + dr * (K1_1 + 2 * K1_2 + 2 * K1_3 + K1_4) / 6.0
			l2[n + 1] = l2[n] + dr * (K2_1 + 2 * K2_2 + 2 * K2_3 + K2_4) / 6.0

			#- rescale to maximum to prevent overflow
			l1 = l1 / np.max(np.abs(l2))
			l2 = l2 / np.max(np.abs(l2))

	#- return -------------------------------------------------------------------------------------
	return l1, l2

@numba.jit( numba.types.UniTuple(numba.float32, 3) (numba.float32[:], numba.float32[:],numba.float32[:], \
        numba.float32, numba.float32[:], numba.float32[:]) )
def group_velocity_sh(l1, l2, r, phase_velocity, rho, N):
    """
    Compute Love wave group velocity for a given set of eigenfunctions (l1, l2), and
    phase velocity (phase_velocity). Further input: radius vector (r) and Earth model (rho, N).
    
    U, I1, I3 = group_velocity_sh(l1, l2, dr, _omega, phase_velocity, model)
    """
    #- compute the integrals I1 and I3 ------------------------------------------------------------
    I1 = 0.0
    I3 = 0.0
    dr = r[1] - r[0]
    for n in xrange(r.size-1):
        I1 = I1 + (rho[n] * l1[n]**2) 
        I3 = I3 + (N[n] * l1[n]**2) 
    I1 = dr * I1
    I3 = dr * I3	
    #- return values ------------------------------------------------------------------------------
    U = I3 / (phase_velocity * I1)
    
    return U, I1, I3

@numba.jit( numba.types.UniTuple(numba.float32[:], 12) (numba.float32[:], numba.float32[:],numba.float32[:],\
    numba.float32, numba.float32, numba.float32, numba.float32[:], numba.float32[:], numba.float32[:], numba.float32[:], \
    numba.float32[:], numba.float32[:]) )
def kernels_sh(r, l1, l2, omega, k, I3, rho, A, C, F, L, N):
    """
    Compute and write sensitivity kernels for Love waves.
    
    kernels_sh(r, l1, l2, _omega, _k, I3, L, write_output, output_directory, tag)
    """
    #- initialisations ----------------------------------------------------------------------------
    
    K_rho_0 = np.zeros(r.size, dtype=np.float32)
    K_A_0   = np.zeros(r.size, dtype=np.float32)
    K_C_0   = np.zeros(r.size, dtype=np.float32)
    K_F_0   = np.zeros(r.size, dtype=np.float32)
    K_L_0   = np.zeros(r.size, dtype=np.float32)
    K_N_0   = np.zeros(r.size, dtype=np.float32)
    
    K_rho   = np.zeros(r.size, dtype=np.float32)
    K_vph   = np.zeros(r.size, dtype=np.float32)
    K_vpv   = np.zeros(r.size, dtype=np.float32)
    K_vsh   = np.zeros(r.size, dtype=np.float32)
    K_vsv   = np.zeros(r.size, dtype=np.float32)
    K_eta   = np.zeros(r.size, dtype=np.float32)
    
    vpv = np.sqrt(A/rho)
    vph = np.sqrt(C/rho)
    vsh = np.sqrt(N/rho)
    vsv = np.sqrt(L/rho)
    eta = F / (A-np.float32(2.)*L)
    #- compute fundamental kernels ----------------------------------------------------------------
    
    K_rho_0 = -(omega**3 * l1**2) / (np.float32(2.) * k**3 * I3)
    K_A_0 = np.zeros(r.size, dtype=np.float32)
    K_C_0 = np.zeros(r.size, dtype=np.float32)
    K_F_0 = np.zeros(r.size, dtype=np.float32)
    K_L_0 = (omega * l1**2) / (np.float32(2.) * k**3 * I3 * L**2)
    K_N_0 = (omega * l1**2) / (np.float32(2.) * k * I3)
    
    #- compute relative kernels in velocity parametrisation ---------------------------------------
    
    K_vph = np.float32(2.)*A*K_A_0 + np.float32(2.)*A*eta*K_F_0
    K_vpv = np.float32(2.)*C*K_C_0
    K_vsh = np.float32(2.)*N*K_N_0
    K_vsv = np.float32(2.)*L*K_L_0 - np.float32(4.)*L*eta*K_F_0
    K_rho = rho*K_rho_0 + N*K_N_0 + L*K_L_0
    K_eta = F*K_F_0
    
    #- convert to relative kernels ----------------------------------------------------------------
    
    K_rho_0 = rho*K_rho_0
    K_L_0   = L*K_L_0
    K_N_0   = N*K_N_0
    return K_rho_0, K_A_0, K_C_0, K_F_0, K_L_0, K_N_0, K_rho, K_vph, K_vpv, K_vsh, K_vsv, K_eta

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
class eigen_solver(object):
    
    def __init__(self, inmodel):
        self.model  = inmodel
        self.dr     = 1000.
        self.nmodes = 1
        return
    
    def init_default(self):
        Tmin   = 5.
        Tmax   = 100.
        dT     = 2.
        cmin   = 3000.
        cmax   = 5000.
        dc     = 50.
        rmin   = 6171.*1000.
        self.T      = _get_array(Tmin, Tmax, dT)
        self.c      = _get_array(cmin, cmax, dc)
        self.omega  = _value_divide_array(2.*np.pi, self.T)
        self.r      = _get_array(rmin, 6371000., self.dr)
        return
        
    def init_dbase(self, T, c, rmin, dr, nmodes):
        
        self.T      = T
        self.c      = c
        self.omega  = _value_divide_array(np.float32(2.*np.pi), self.T)
        self.dr     = dr
        self.r      = _get_array(rmin, 6371000., self.dr)
        self.nmodes = nmodes
        return
    
    def init_output(self, wavetype):
        """
        Krho0data   - density kernel in terms of Love parameters
        Krhodata    - density kernel in terms of velocities
        """
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
        self.Krho0data  = tempdata.copy()
        self.Krhodata   = tempdata.copy()
        self.Kvphdata   = tempdata.copy()
        self.Kvpvdata   = tempdata.copy()
        self.Kvshdata   = tempdata.copy()
        self.Kvsvdata   = tempdata.copy()
        self.Ketadata   = tempdata.copy()
        return
    
    def solve_PSV(self):
        #- root-finding algorithm ---------------------------------------------------------------------
        self.init_output(1)

        rho = np.zeros(self.r.size, dtype=np.float32)
        A   = np.zeros(self.r.size, dtype=np.float32)
        C   = np.zeros(self.r.size, dtype=np.float32)
        F   = np.zeros(self.r.size, dtype=np.float32)
        L   = np.zeros(self.r.size, dtype=np.float32)
        N   = np.zeros(self.r.size, dtype=np.float32)
    
        for n in xrange(self.r.size): rho[n], A[n], C[n], F[n], L[n], N[n] = self.model.get_r_love_parameters_PSV(self.r[n])
        #- loop over angular frequencies
        for it in xrange(self.omega.size):
            omega       = self.omega[it]
            mode_count  = 0
            k           = _value_divide_array(omega, self.c)
            #- loop over trial wavenumbers
            r_left = 0.0
            for n in xrange(k.size):
                if self.nmodes == mode_count: break
                # - compute vertical wave functions using alternative system
                r1, r2, r3, r4, r5= integrate_psv_alt(self.model, self.r, self.dr, omega, k[n])
                r_right = r2[-1]
                #- check if there is a zero -----------------------------------------------------------
                if r_left * r_right < 0.0:
                    mode_count += 1
                    #- start bisection algorithm
                    rr_left     = r_left
                    rr_right    = r_right
                    k_left      = k[n-1]
                    k_right     = k[n]
                    for i in xrange(5):
                        k_new = (k_left * np.abs(rr_right) + k_right * np.abs(rr_left)) / (np.abs(rr_left) + np.abs(rr_right))
                        r1, r2, r3, r4, r5  = integrate_psv_alt(self.model, self.r, self.dr, omega, k_new)
                        rr = r2[len(r2)-1]
                        if rr * rr_left < 0.0:
                            k_right = k_new
                            rr_right = rr
                        elif rr * rr_right < 0.0:
                            k_left = k_new
                            rr_left = rr
                        else:
                            continue
                    #==================================================================================
                    #- compute final vertical wave functions and corresponding velocities and kernels -
                    #==================================================================================
                    
                    #- compute final vertical wave functions using the original first-order system
                    #- two independent solutions
                    r11, r21, r31, r41 = integrate_psv(self.model, self.r, self.dr, omega, k_new, 1)
                    r12, r22, r32, r42 = integrate_psv(self.model, self.r, self.dr, omega, k_new, 2)
                    #- determine their weights via boundary condition (weight q1 is set to 1)
                    q2 = -r21[-1] / r22[-1]
                    #- final solution with boundary condition
                    r1 = r11 + q2 * r12
                    r2 = r21 + q2 * r22
                    r3 = r31 + q2 * r32
                    r4 = r41 + q2 * r42
                    #- normalise to vertical displacement at the surface
                    mm = r1[-1]
                    r1 = r1 / mm
                    r2 = r2 / mm
                    r3 = r3 / mm
                    r4 = r4 / mm
                    # - store eigenfunction data
                    self.r1data[mode_count-1, it, :] = r1[:]
                    self.r2data[mode_count-1, it, :] = r2[:]
                    self.r3data[mode_count-1, it, :] = r3[:]
                    self.r4data[mode_count-1, it, :] = r4[:]
                    # -
                    self.eArr[mode_count-1, it] = 1
                    #- phase velocity
                    self.Vph[mode_count-1, it] = omega / k_new
                    #- group velocity
                    U, I1, I3 = group_velocity_psv(r1, r2, r3, r4, self.r, k_new, omega/k_new, rho, A, C, F, L, N)
                    self.Vgr[mode_count-1, it] = U
                    #- kernels
                    K_rho_0, K_A_0, K_C_0, K_F_0, K_L_0, K_N_0, K_rho, K_vph, K_vpv, K_vsh, K_vsv, K_eta = \
                            kernels_psv(self.r, r1, r2, r3, r4, omega, k_new, I3, rho, A, C, F, L, N)
                    self.Krho0data[mode_count-1, it, :] = K_rho_0[:]
                    self.Kadata[mode_count-1, it, :]    = K_A_0[:]
                    self.Kcdata[mode_count-1, it, :]    = K_C_0[:]
                    self.Kfdata[mode_count-1, it, :]    = K_F_0[:]
                    self.Kldata[mode_count-1, it, :]    = K_L_0[:]
                    self.Kndata[mode_count-1, it, :]    = K_N_0[:]
                    self.Krhodata[mode_count-1, it, :]  = K_rho[:]
                    self.Kvphdata[mode_count-1, it, :]  = K_vph[:]
                    self.Kvpvdata[mode_count-1, it, :]  = K_vpv[:]
                    self.Kvshdata[mode_count-1, it, :]  = K_vsh[:]
                    self.Kvsvdata[mode_count-1, it, :]  = K_vsv[:]
                    self.Ketadata[mode_count-1, it, :]  = K_eta[:]
                r_left = r_right # in the loop over k
        return
    
    def solve_SH(self):
        #- root-finding algorithm ---------------------------------------------------------------------
        self.init_output(2)
        rho = np.zeros(self.r.size, dtype=np.float32)
        A   = np.zeros(self.r.size, dtype=np.float32)
        C   = np.zeros(self.r.size, dtype=np.float32)
        F   = np.zeros(self.r.size, dtype=np.float32)
        L   = np.zeros(self.r.size, dtype=np.float32)
        N   = np.zeros(self.r.size, dtype=np.float32)
    
        for n in xrange(self.r.size): rho[n], A[n], C[n], F[n], L[n], N[n] = self.model.get_r_love_parameters_SH(self.r[n])
        #- loop over angular frequencies
        for it in xrange(self.omega.size):
            omega       = self.omega[it]
            mode_count  = 0
            k           = _value_divide_array(omega, self.c)
            #- loop over trial wavenumbers
            l_left = 0.0
            for n in xrange(k.size):
                if self.nmodes == mode_count: break
                #- compute vertical wave functions
                l1, l2 = integrate_sh(self.model, self.r, self.dr, omega, k[n])
                l_right = l2[-1]
                
                #- check if there is a zero -----------------------------------------------------------
                if l_left * l_right < 0.0:
                    mode_count += 1
                    #- start bisection algorithm
                    ll_left = l_left
                    ll_right= l_right
                    k_left  = k[n-1]
                    k_right = k[n]
                    for i in xrange(5):
                        k_new   = (k_left * np.abs(ll_right) + k_right * np.abs(ll_left)) / (np.abs(ll_left) + np.abs(ll_right))
                        l1, l2  = integrate_sh(self.model, self.r, self.dr, omega, k_new)
                        ll      = l2[-1]
                        if ll * ll_left < 0.0:
                            k_right = k_new
                            ll_right= ll
                        elif ll * ll_right < 0.0:
                            k_left  = k_new
                            ll_left = ll
                        else:
                            continue
                    #==================================================================================
                    #- compute final vertical wave functions and corresponding velocities and kernels -
                    #==================================================================================
                    #- stress and displacement functions 
                    l1, l2 = integrate_sh(self.model, self.r, self.dr, omega, k_new)
                    # -
                    self.eArr[mode_count-1, it] = 1
                    #- phase velocity
                    self.Vph[mode_count-1, it] = omega / k_new
                    #- group velocity
                    U, I1, I3 = group_velocity_sh(l1, l2, self.r, omega/k_new, rho, N)
                    self.Vgr[mode_count-1, it] = U
                    #- kernels
                    K_rho_0, K_A_0, K_C_0, K_F_0, K_L_0, K_N_0, K_rho, K_vph, K_vpv, K_vsh, K_vsv, K_eta = \
                        kernels_sh(self.r, l1, l2, omega, k_new, I3, rho, A, C, F, L, N)
                    self.Krho0data[mode_count-1, it, :] = K_rho_0[:]
                    self.Kadata[mode_count-1, it, :]    = K_A_0[:]
                    self.Kcdata[mode_count-1, it, :]    = K_C_0[:]
                    self.Kfdata[mode_count-1, it, :]    = K_F_0[:]
                    self.Kldata[mode_count-1, it, :]    = K_L_0[:]
                    self.Kndata[mode_count-1, it, :]    = K_N_0[:]
                    self.Krhodata[mode_count-1, it, :]  = K_rho[:]
                    self.Kvphdata[mode_count-1, it, :]  = K_vph[:]
                    self.Kvpvdata[mode_count-1, it, :]  = K_vpv[:]
                    self.Kvshdata[mode_count-1, it, :]  = K_vsh[:]
                    self.Kvsvdata[mode_count-1, it, :]  = K_vsv[:]
                    self.Ketadata[mode_count-1, it, :]  = K_eta[:]
                l_left = l_right # in the loop over k
        return
    