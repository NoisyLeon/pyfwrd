
import numba
import numpy as np
import vmodel
import asdf

# define type of vmodel.model1d
model_type = numba.deferred_type()
model_type.define(vmodel.model1d.class_type.instance_type)

#--------------------------------------------------------------------------------------------------
#- fundamental array manipulations
#--------------------------------------------------------------------------------------------------
@numba.jit(numba.float64[:](numba.float64, numba.float64, numba.float64))
def _get_array(xmin, xmax, dx):
    xlst= []
    Nx  = int((xmax - xmin)/dx + 1)
    for i in xrange(Nx): xlst.append(dx*i+xmin)
    return np.array(xlst, dtype=np.float64)

@numba.jit(numba.float64[:](numba.float64, numba.float64[:]))
def _value_divide_array(value, array):
    outArr  = np.zeros(array.size, dtype=np.float64)
    for i in xrange(array.size): outArr[i] = value/array[i]
    return outArr

@numba.jit(numba.float64[:](numba.float64, numba.float64[:]))
def _array_divide_value(value, array):
    outArr  = np.zeros(array.size, dtype=np.float64)
    for i in xrange(array.size): outArr[i] = array[i]/value
    return outArr

@numba.jit(numba.float64[:](numba.float64[:], numba.float64[:]))
def _merge_array(a1, a2):
    a3  = np.zeros(a1.size+a2.size, dtype=np.float64)
    Na1 = a1.size
    for i in xrange(a3.size):
        if i <= Na1-1:
            a3[i] = a1[i]
        else:
            a3[i] = a2[i-Na1]
    return a3

@numba.jit(numba.float64(numba.float64[:]))
def _abs_max_(array):
    mvalue=np.abs(array[0])
    for i in xrange(array.size):
        if np.abs(array[i]) > mvalue: mvalue = np.abs(array[i])
    return mvalue

#--------------------------------------------------------------------------------------------------
#- fundamental functions for integrate_psv_alt
#--------------------------------------------------------------------------------------------------
@numba.jit(numba.float64(numba.float64, numba.float64, numba.float64, numba.float64))
def f1_psv_alt(C,L,r4,r5): return (r4 / L - r5 / C)

@numba.jit(numba.float64(numba.float64, numba.float64, numba.float64, numba.float64, numba.float64,
        numba.float64, numba.float64, numba.float64))
def f2_psv_alt(rho,A,C,F,omega,k,r4,r5): return (-omega**2 * rho * r4 + (omega**2 * rho - k**2 * (A - F**2 / C)) * r5)

@numba.jit(numba.float64(numba.float64, numba.float64, numba.float64, numba.float64, numba.float64))
def f3_psv_alt(C,F,k,r4,r5): return (k * r4 + k * F * r5 / C)

@numba.jit(numba.float64(numba.float64, numba.float64, numba.float64, numba.float64, numba.float64,
        numba.float64, numba.float64, numba.float64, numba.float64))
def f4_psv_alt(rho,A,C,F,omega,k,r1,r2,r3): return ((-omega**2 * rho + k**2 * (A - F**2 / C)) * r1 + r2 / C - 2 * k * F * r3 / C)

@numba.jit(numba.float64(numba.float64, numba.float64, numba.float64, numba.float64, numba.float64,
        numba.float64, numba.float64))
def f5_psv_alt(rho,L,omega,k,r1,r2,r3):
	return (omega**2 * rho * r1 - r2 / L - 2 * k * r3)


#--------------------------------------------------------------------------------------------------
#- numerical integration
#--------------------------------------------------------------------------------------------------

    
@numba.njit( numba.types.UniTuple(numba.float64[:], 5) (model_type, numba.float64[:], numba.float64, numba.float64, numba.float64) )
def integrate_psv_alt(model, r, dr, omega, k):
# # # @numba.jit( numba.types.UniTuple(numba.float64[:], 5) (model_type, numba.float64[:], numba.float64[:], numba.float64, numba.float64) )
# # # def integrate_psv_alt(model, r, drArr, omega, k):
    """
    Integrate first-order system for a fixed circular frequency omega and a fixed wavenumber k.
    r1, r2, r3, r4, r5, r = integrate_psv_alt(r_min, dr, omega, k, model)
    
    r_min:		minimum radius in m
    dr:			radius increment in m
    omega:		circular frequency in Hz
    k:			wave number in 1/m
    model:		Earth model, e.g. "PREM", "GUTENBERG", ... .
    
    r1, ...:	variables of the alternative Rayleigh wave system
    r:			radius vector in m
    """
    
    #- initialisation -----------------------------------------------------------------------------
    
    # r   = _get_array(rmin, 6371000., dr)
    r1  = np.zeros(r.size, dtype=np.float64)
    r2  = np.zeros(r.size, dtype=np.float64)
    r3  = np.zeros(r.size, dtype=np.float64)
    r4  = np.zeros(r.size, dtype=np.float64)
    r5  = np.zeros(r.size, dtype=np.float64)
    
    rho, A, C, F, L, N = model.get_ind_Love_parameters(0)

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
            rho, A, C, F, L, N = model.get_r_love_parameters(r[n])
            
            K1_1 = f1_psv_alt(C,L,r4[n],r5[n])
            K2_1 = f2_psv_alt(rho,A,C,F,omega,k,r4[n],r5[n])
            K3_1 = f3_psv_alt(C,F,k,r4[n],r5[n])
            K4_1 = f4_psv_alt(rho,A,C,F,omega,k,r1[n],r2[n],r3[n])
            K5_1 = f5_psv_alt(rho,L,omega,k,r1[n],r2[n],r3[n]) 
        
            rho, A, C, F, L, N = model.get_r_love_parameters(r[n]+dr/2.0)
            
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
            
            rho, A, C, F, L, N = model.get_r_love_parameters(r[n]+dr)
            
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
            # # # # # mm  = _abs_max_(r2)
            mm  = np.max(np.abs(r2))
            # r1  = _array_divide_value(mm, r1)
            # r2  = _array_divide_value(mm, r2)
            # r3  = _array_divide_value(mm, r3)
            # r4  = _array_divide_value(mm, r4)
            # r5  = _array_divide_value(mm, r5)
            
            r1  = r1 / mm
            r2  = r2 / mm
            r3  = r3 / mm
            r4  = r4 / mm
            r5  = r5 / mm
    
    #- return -------------------------------------------------------------------------------------
    
    return r1, r2, r3, r4, r5

#--------------------------------------------------------------------------------------------------
#- fundamental functions for integrate_psv
#--------------------------------------------------------------------------------------------------
@numba.jit(numba.float64 (numba.float64, numba.float64, numba.float64, numba.float64, numba.float64) )
def f1_psv(C,F,k,r2,r3):
	return (r2 / C + k * F * r3 / C)

@numba.jit(numba.float64 (numba.float64, numba.float64, numba.float64, numba.float64, numba.float64) )
def f2_psv(rho,omega,k,r1,r4):
	return (-rho * omega**2 * r1 + k * r4)

@numba.jit(numba.float64 (numba.float64, numba.float64, numba.float64, numba.float64) )
def f3_psv(L,k,r1,r4):
	return (r4 / L - k * r1)

@numba.jit(numba.float64 (numba.float64, numba.float64, numba.float64, numba.float64, numba.float64,
         numba.float64, numba.float64, numba.float64) )
def f4_psv(rho,A,C,F,omega,k,r2,r3):
	return (-k * F * r2 / C + (k**2 * (A - F**2 / C) - rho * omega**2) * r3)

# 
@numba.jit( numba.types.UniTuple(numba.float64[:], 4) (model_type, numba.float64[:], numba.float64, numba.float64, numba.float64, numba.int32) )
def integrate_psv(model, r, dr, omega, k, initial_condition):

# @numba.jit( numba.types.UniTuple(numba.float64[:], 4) (model_type, numba.float64[:], numba.float64[:], numba.float64, numba.float64, numba.int32) )
# def integrate_psv(model, r, drArr, omega, k, initial_condition):
    """
    Integrate first-order system for a fixed circular frequency omega and a fixed wavenumber k.
    r1, r2, r3, r4, r = integrate_psv(r_min, dr, omega, k, model)
    
    r_min:		minimum radius in m
    dr:			radius increment in m
    omega:		circular frequency in Hz
    k:			wave number in 1/m
    model:		Earth model, e.g. "PREM", "GUTENBERG", ... .
    
    r1, ...:	variables of the Rayleigh wave system
    r:			radius vector in m
    """
    
    #- initialisation -----------------------------------------------------------------------------
    # initial_condition=1
    # r = np.arange(r_min, 6371000.0 + dr, dr, dtype=float)
    r1  = np.zeros(r.size, dtype=np.float64)
    r2  = np.zeros(r.size, dtype=np.float64)
    r3  = np.zeros(r.size, dtype=np.float64)
    r4  = np.zeros(r.size, dtype=np.float64)
    rho, A, C, F, L, N = model.get_ind_Love_parameters(0)
    
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
            # # # dr=drArr[n]
            #- compute Runge-Kutta coeficients for l1 and l2
            rho, A, C, F, L, N = model.get_r_love_parameters(r[n])
            K1_1 = f1_psv(C,F,k,r2[n],r3[n])
            K2_1 = f2_psv(rho,omega,k,r1[n],r4[n])
            K3_1 = f3_psv(L,k,r1[n],r4[n])
            K4_1 = f4_psv(rho,A,C,F,omega,k,r2[n],r3[n])
    
            rho, A, C, F, L, N = model.get_r_love_parameters(r[n] + dr/2.0)
            K1_2 = f1_psv(C,F,k,r2[n]+0.5*K2_1*dr,r3[n]+0.5*K3_1*dr)
            K2_2 = f2_psv(rho,omega,k,r1[n]+0.5*K1_1*dr,r4[n]+0.5*K4_1*dr)
            K3_2 = f3_psv(L,k,r1[n]+0.5*K1_1*dr,r4[n]+0.5*K4_1*dr)
            K4_2 = f4_psv(rho,A,C,F,omega,k,r2[n]+0.5*K2_1*dr,r3[n]+0.5*K3_1*dr)
    
            K1_3 = f1_psv(C,F,k,r2[n]+0.5*K2_2*dr,r3[n]+0.5*K3_2*dr)
            K2_3 = f2_psv(rho,omega,k,r1[n]+0.5*K1_2*dr,r4[n]+0.5*K4_2*dr)
            K3_3 = f3_psv(L,k,r1[n]+0.5*K1_2*dr,r4[n]+0.5*K4_2*dr)
            K4_3 = f4_psv(rho,A,C,F,omega,k,r2[n]+0.5*K2_2*dr,r3[n]+0.5*K3_2*dr)
            
            rho, A, C, F, L, N = model.get_r_love_parameters(r[n] + dr)
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
    
    #- return -------------------------------------------------------------------------------------
    
    return r1, r2, r3, r4

@numba.jit( numba.types.UniTuple(numba.float64, 3) (numba.float64[:], numba.float64[:],numba.float64[:], numba.float64[:],\
    numba.float64[:], numba.float64, numba.float64, numba.float64[:], numba.float64[:], numba.float64[:], numba.float64[:], \
    numba.float64[:], numba.float64[:]) )
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



###################################################################################
#
###################################################################################


spec = [('model', model_type),
        ('Tmin', numba.float64),
        ('Tmax', numba.float64),
        ('dT', numba.float64),
        ('cmin', numba.float64),
        ('cmax', numba.float64),
        ('dc', numba.float64),
        ('dr', numba.float64),
        ('rmin', numba.float64),
        ('r', numba.float64[:]),
        ('T', numba.float64[:]),
        ('c', numba.float64[:]),
        ('CR', numba.float64[:]),
        ('CL', numba.float64[:]),
        ('UR', numba.float64[:]),
        ('CL', numba.float64[:]),
        ('omega', numba.float64[:]),
        ('nmodes', numba.int32)
        ]

@numba.jitclass(spec)
class eigen_solver(object):
    
    def __init__(self, inmodel):
        self.model  = inmodel
        self.Tmin   = 10.
        self.Tmax   = 50.
        self.dT     = 5.
        self.cmin   = 3000.
        self.cmax   = 5000.
        self.dc     = 100.
        self.rmin   = 6171.*1000.
        self.dr     = 1000.
        self.nmodes = 1
        return
    
    # def initialization(self):
    #     self.T      = np.arange(self.Tmin,self.Tmax + self.dT, self.dT, dtype=np.float64)
    #     self.c      = np.arange(self.cmin,self.cmax + self.dc ,self.dc , dtype=np.float64)
    #     self.omega  = 2 * np.pi / T
    
    def initialization(self):
        
        self.T      = _get_array(self.Tmin, self.Tmax, self.dT)
        self.c      = _get_array(self.cmin, self.cmax, self.dc)
        self.omega  = _value_divide_array(2.*np.pi, self.T)
        self.r      = _get_array(self.rmin, 6371000., self.dr)
        # r1          = _get_array(6171000., 6371000., 1000.)
        # dr1         = np.ones(r1.size, dtype=np.float64)
        # r2          = _get_array(self.rmin, 6171000.-5000., 5000.)
        # dr2         = np.ones(r2.size, dtype=np.float64)*5.
        # self.r      = _merge_array(r2, r1)
        # ###
        # self.dr     = _merge_array(dr2, dr1)
        ###
    def assign(self, T, c, r, dr, nmodes):
        
        self.T      = T
        self.c      = c
        self.omega  = _value_divide_array(2.*np.pi, self.T)
        self.r      = r
        self.dr     = dr
        self.nmodes = nmodes
        
        return
        
    
    def solve_PSV(self):
        #- root-finding algorithm ---------------------------------------------------------------------
        mode = []
        phase_velocities = []
        group_velocities = []
        rho = np.zeros(self.r.size, dtype=np.float64)
        A   = np.zeros(self.r.size, dtype=np.float64)
        C   = np.zeros(self.r.size, dtype=np.float64)
        F   = np.zeros(self.r.size, dtype=np.float64)
        L   = np.zeros(self.r.size, dtype=np.float64)
        N   = np.zeros(self.r.size, dtype=np.float64)
    
        for n in xrange(self.r.size): rho[n], A[n], C[n], F[n], L[n], N[n] = self.model.get_r_love_parameters(self.r[n])
        #- loop over angular frequencies
        for omega in self.omega:
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
                    mode.append(mode_count)
                
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
                    
                    #- phase velocity
                    # periods.append(2*np.pi/_omega)
                    phase_velocities.append(omega / k_new)
                    #- group velocity
                    U, I1, I3 = group_velocity_psv(r1, r2, r3, r4, self.r, k_new, omega/k_new, rho, A, C, F, L, N)
                    group_velocities.append(U)
                    #- kernels
                    # kpsv.kernels_psv(r, r1, r2, r3, r4, _omega, k_new, I3, rho, A, C, F, L, N, write_output, output_directory, tag)
                r_left = r_right # in the loop over k
        self.CR = np.array(phase_velocities, dtype=np.float64)
        self.UR = np.array(group_velocities, dtype=np.float64)
        return
    
class eigenASDF(asdf.AsdfFile):
    """
    ASDF database for eigenfunction/dispersion computation and data storage
    Note that the ASDF here is Advanced Scientific Data Format, NOT Adaptable Seismic Data Format !
    """
    def init_dbase(self, inmodel=None, love=True, rayleigh=True, Tmin=10., Tmax=50., dT=5., cmin=2500., cmax=5000., dc=50.,
            nmodes=1, zmax=[500., 200.], zArr=None, dz=[5., 1.]):
        if isinstance(inmodel, vmodel.model1d):
            self.eigen  = eigen_solver(inmodel)
        else:
            m=vmodel.model1d()
            m.get_radius(400., 1.)
            m.model_prem()
            self.eigen  = eigen_solver(m)
        self.love       = love
        self.rayleigh   = rayleigh
        if zArr !=None:
            r = np.float64(6371000. - zArr*1000.)
        self.eigen.assign(T, c, r, dr, )
        
            
        
        
    
    