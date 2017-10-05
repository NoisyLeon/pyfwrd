# -*- coding: utf-8 -*-
"""
Module for surface wave dispersion computation using normal mode method by Robert Herrmann

The code is a python wrapper of the f77 code Computer Program in Seismology 
Numba is NOT used the object tcps_solver

:Copyright:
    Author: Lili Feng
    Graduate Research Assistant
    CIEI, Department of Physics, University of Colorado Boulder
    email: lili.feng@colorado.edu
"""

import tdisp96
import tregn96, tlegn96
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


class tcps_solver(object):
    """
    An object solving for dispersion curves, eigenfunctions and sensitivity kernels, use tdisp96 and tregn96/tlegn96
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
        self.egn96  = True
        self.tilt   = False
        return
    
    def init_default(self, dh=1., nl=200):
        Tmin        = 5.
        Tmax        = 100.
        dT          = 5.
        self.cmin   = -1.
        self.cmax   = -1.
        self.T      = _get_array(Tmin, Tmax, dT)
        self.freq   = _value_divide_array(1., self.T)
        self.dArr   = np.ones(nl, dtype = np.float32)*np.float32(dh)
        return
    
    def init_default_3(self, dh=1., nl=200):
        Tmin        = 5.
        Tmax        = 5.
        dT          = 5.
        self.cmin   = -1.
        self.cmax   = -1.
        self.T      = _get_array(Tmin, Tmax, dT)
        self.freq   = _value_divide_array(1., self.T)
        self.dArr   = np.ones(nl, dtype = np.float32)*np.float32(dh)
        return
    
    def init_default_2(self):
        Tmin        = 5.
        Tmax        = 100.
        dT          = 5.
        self.cmin   = -1.
        self.cmax   = -1.
        self.T      = _get_array(Tmin, Tmax, dT)
        self.freq   = _value_divide_array(1., self.T)
        self.dArr   = np.array([20.,  15.,  42.,  43.,  45.,  35.], dtype = np.float32)
        return
    
    def love2vel(self):
        """
        Love parameters to velocity parameters
        ::: Output :::
        ahArr   - vph (unit -km/s)
        avArr   - vpv (unit -km/s)
        bhArr   - vsh (unit -km/s)
        bvArr   - vsv (unit -km/s)
        nArr    - eta (dimensionless)
        """
        # # # if self.model.flat == 1:
        # # #     self.ahArr  = np.sqrt(self.AArr/self.rhoArr)
        # # #     self.avArr  = np.sqrt(self.CArr/self.rhoArr)
        # # #     self.bhArr  = np.sqrt(self.NArr/self.rhoArr)
        # # #     self.bvArr  = np.sqrt(self.LArr/self.rhoArr)
        # # #     self.nArr   = self.FArr/(self.AArr - np.float32(2.)* self.LArr)
        # # # else:
        # # #     self.ahArr  = np.sqrt(self.Asph/self.rhosph)
        # # #     self.avArr  = np.sqrt(self.Csph/self.rhosph)
        # # #     self.bhArr  = np.sqrt(self.Nsph/self.rhosph)
        # # #     self.bvArr  = np.sqrt(self.Lsph/self.rhosph)
        # # #     self.nArr   = self.Fsph/(self.Asph - np.float32(2.)* self.Lsph)
        self.ahArr  = np.sqrt(self.AArr/self.rhoArr)
        self.avArr  = np.sqrt(self.CArr/self.rhoArr)
        self.bhArr  = np.sqrt(self.NArr/self.rhoArr)
        self.bvArr  = np.sqrt(self.LArr/self.rhoArr)
        self.nArr   = self.FArr/(self.AArr - np.float32(2.)* self.LArr)
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
    
    def solve_PSV(self):
        """
        Solve P-SV motion dispersion curves, eigenfunctions and sensitivity kernels
        ===================================================================================================
        ::: Output :::
        : Model :
        dArr, AArr, CArr, FArr, LArr, NArr, rhoArr  - layerized model
        dsph, Asph, Csph, Fsph, Lsph, Nsph, rhosph  - flattening transformed Earth model
        : dispersion :
        Vph, Vgr    - phase/group velocity
        : eigenfunctions :
        uz, ur      - vertical/radial displacement functions
        tuz, tur    - vertical/radial stress functions
        duzdz, durdz- derivatives of vertical/radial displacement functions
        : velocity/density sensitivity kernels :
        dcdah/dcdav - vph/vpv kernel
        dcdbh/dcdbv - vsh/vsv kernel
        dcdn        - eta kernel
        dcdr        - density kernel
        : Love parameter/density sensitivity kernels, derived from the kernels above using chain rule :
        dcdA, dcdC, dcdF, dcdL, dcdN    - Love parameter kernels
        dcdrl                           - density kernel
        ===================================================================================================
        """
        #- root-finding algorithm using tdisp96, compute phase velocities ------------------------
        if self.model.tilt == 1:
            dArr, rhoArr, AArr, CArr, FArr, LArr, NArr, BcArr, BsArr, GcArr, GsArr, HcArr, HsArr, CcArr, CsArr =\
                        self.model.get_layer_tilt_model(self.dArr, 200, 1.)
            self.tilt   = True
            self.egn96  = True
            self.BcArr  = BcArr
            self.BsArr  = BsArr
            self.GcArr  = GcArr
            self.GsArr  = GsArr
            self.HcArr  = HcArr
            self.HsArr  = HsArr
            self.CcArr  = CcArr
            self.CsArr  = CsArr
        else:
            dArr, rhoArr, AArr, CArr, FArr, LArr, NArr = self.model.get_layer_model(self.dArr, 200, 1.)
        nfval   = self.freq.size
        if self.model.flat == 0:
            iflsph_in = 1
        else:
            iflsph_in = 0
        ilvry   = 2
        c_out,d_out,TA_out,TC_out,TF_out,TL_out,TN_out,TRho_out = tdisp96.disprs(ilvry, 1., nfval, 1, self.verbose, nfval, \
                np.append(self.freq, np.zeros(2049-nfval)), self.cmin,self.cmax, dArr, AArr,CArr,FArr,LArr,NArr,rhoArr, dArr.size,\
                iflsph_in, 0., self.nmodes, 0.5, 0.5)
        self.ilvry  = ilvry
        self.Vph    = c_out[:nfval]
        # Save flat transformed model is spherical flag is on
        if self.model.flat == 0:
            self.dsph   = d_out
            self.Asph   = TA_out
            self.Csph   = TC_out
            self.Lsph   = TL_out
            self.Fsph   = TF_out
            self.Nsph   = TN_out
            self.rhosph = TRho_out
        self.AArr   = AArr
        self.CArr   = CArr
        self.LArr   = LArr
        self.FArr   = FArr
        self.NArr   = NArr
        self.rhoArr = rhoArr
        self.love2vel()
        #- compute eigenfunction/kernels
        if self.egn96:
            hs_in       = 0.
            hr_in       = 0.
            ohr_in      = 0.
            ohs_in      = 0.
            refdep_in   = 0.
            dogam       = False # No attenuation
            nl_in       = dArr.size
            k           = 2.*np.pi/self.Vph/self.T
            k2d         = np.tile(k, (nl_in, 1))
            k2d         = k2d.T
            omega       = 2.*np.pi/self.T
            omega2d     = np.tile(omega, (nl_in, 1))
            omega2d     = omega2d.T
            if self.model.flat == 0:
                d_in    = d_out
                TA_in   = TA_out
                TC_in   = TC_out
                TF_in   = TF_out
                TL_in   = TL_out
                TN_in   = TN_out
                TRho_in = TRho_out
            else:
                d_in    = dArr
                TA_in   = AArr
                TC_in   = CArr
                TF_in   = FArr
                TL_in   = LArr
                TN_in   = NArr
                TRho_in = rhoArr
            qai_in      = np.ones(nl_in)*1.e6
            qbi_in      = np.ones(nl_in)*1.e6
            etapi_in    = np.zeros(nl_in)
            etasi_in    = np.zeros(nl_in)
            frefpi_in   = np.ones(nl_in)
            frefsi_in   = np.ones(nl_in)
            u_out, ur, tur, uz, tuz, dcdh,dcdav,dcdah,dcdbv,dcdbh,dcdn,dcdr = tregn96.tregn96(hs_in, hr_in, ohr_in, ohs_in,\
                refdep_in, dogam, nl_in, iflsph_in, d_in, TA_in, TC_in, TF_in, TL_in, TN_in, TRho_in, \
                qai_in,qbi_in,etapi_in,etasi_in, frefpi_in, frefsi_in, self.T.size, self.T, self.Vph)
            ######################################################
            # store output
            ######################################################
            self.Vgr    = u_out
            # eigenfunctions
            self.uz     = uz[:nfval,:nl_in]
            self.tuz    = tuz[:nfval,:nl_in]
            self.ur     = ur[:nfval,:nl_in]
            self.tur    = tur[:nfval,:nl_in]
            # sensitivity kernels for velocity parameters and density
            self.dcdah  = dcdah[:nfval,:nl_in]
            self.dcdav  = dcdav[:nfval,:nl_in]
            self.dcdbh  = dcdbh[:nfval,:nl_in]
            self.dcdbv  = dcdbv[:nfval,:nl_in]
            self.dcdr   = dcdr[:nfval,:nl_in]
            self.dcdn   = dcdn[:nfval,:nl_in]
            # Love parameters and density in the shape of nfval, nl_in
            # # A           = np.tile(TA_in, (nfval,1))
            # # C           = np.tile(TC_in, (nfval,1))
            # # F           = np.tile(TF_in, (nfval,1))
            # # L           = np.tile(TL_in, (nfval,1))
            # # N           = np.tile(TN_in, (nfval,1))
            # # rho         = np.tile(TRho_in, (nfval,1))
            # # d           = np.tile(d_in, (nfval,1))
            A           = np.tile(AArr, (nfval,1))
            C           = np.tile(CArr, (nfval,1))
            F           = np.tile(FArr, (nfval,1))
            L           = np.tile(LArr, (nfval,1))
            N           = np.tile(NArr, (nfval,1))
            rho         = np.tile(rhoArr, (nfval,1))
            d           = np.tile(dArr, (nfval,1))
            eta         = F/(A-2.*L)
            # derivative of eigenfunctions, $5.8 of R.Herrmann
            self.durdz  = 1./L*self.tur - k2d*self.uz
            self.duzdz  = k2d*F/C*self.ur + self.tuz/C
            # integrals for the Lagranian
            # # # I0          = np.sum( rho*(self.uz**2+self.ur**2)*d, axis=1)
            # # # I1          = np.sum( (L*(self.uz**2)+A*(self.ur**2))*d, axis=1)
            # # # I2          = np.sum( (L*self.uz*self.durdz - F*self.ur*self.duzdz)*d, axis=1)
            # # # I3          = np.sum( (L*(self.durdz**2)+C*(self.duzdz**2))*d, axis=1)
            # # U           = (k*I1+I2)/omega/I0
            # # # I02d        = np.tile(I0, (nl_in, 1))
            # # # I02d        = I02d.T
            # # # U2d         = np.tile(self.Vgr, (nl_in, 1))
            # # # U2d         = U2d.T
            # # # C2d         = np.tile(self.Vph, (nl_in, 1))
            # # # C2d         = C2d.T
            # sensitivity kernels for Love parameters and density, derived from velocity kernels using chain rule
            # NOTE: benchmark of love kernels predict somewhat different perturbed values compared with those using velocity kernel
            # Possible reason:
            #       By using the chain rule, infinitesimal perturbations for all the parameters is assumed.
            #       However, for the benckmark cases, the perturbation can be large.
            #       e.g. ah = sqrt(A/rho), dah = 0.5/sqrt(A*rho)*dA. But, if dA is large, dah = sqrt(A+dA/rho) - sqrt(A/rho) != 0.5/sqrt(A*rho)*dA
            #      
            self.dcdA   = 0.5/np.sqrt(A*rho)*self.dcdah - F/((A-2.*L)**2)*self.dcdn
            self.dcdC   = 0.5/np.sqrt(C*rho)*self.dcdav
            self.dcdF   = 1./(A-2.*L)*self.dcdn
            self.dcdL   = 0.5/np.sqrt(L*rho)*self.dcdbv + 2.*F/((A-2.*L)**2)*self.dcdn
            self.dcdN   = 0.5/np.sqrt(N*rho)*self.dcdbh
            # density kernel
            self.dcdrl  = -0.5*self.dcdah*np.sqrt(A/(rho**3)) - 0.5*self.dcdav*np.sqrt(C/(rho**3))\
                            -0.5*self.dcdbh*np.sqrt(N/(rho**3)) -0.5*self.dcdbv*np.sqrt(L/(rho**3))\
                                + self.dcdr
            # U2d3         = np.tile(U, (nl_in, 1))
            # U2d3         = U2d3.T
            # # # I0          = np.sum( rho*(self.uz**2+self.ur**2), axis=1)
            # # # I1          = np.sum( (L*(self.uz**2)+A*(self.ur**2)), axis=1)
            # # # I2          = np.sum( (L*self.uz*self.durdz - F*self.ur*self.duzdz), axis=1)
            # # # I3          = np.sum( (L*(self.durdz**2)+C*(self.duzdz**2)), axis=1)
            # self.I0 = I0
            # self.I1 = I1
            # self.I2 = I2
            # self.I3 = I3
            ##############################################
            # For benchmark purposes
            ##############################################
            # # derivative of eigenfunctions
            # self.durdz1 = -(self.ur[:,:-1] - self.ur[:,1:])/self.dArr[0]
            # self.duzdz1 = -(self.uz[:,:-1] - self.uz[:,1:])/self.dArr[0]
            # self.Vgr1   = (k*I1+I2)/omega/I0
            # 
            # # derived kernels from eigenfunctions, p 299 in Herrmann's notes
            # self.dcdah1 = eta*rho*np.sqrt(A/rho)*(self.ur**2 - 2.*eta/k2d*self.ur*self.duzdz)/U2d/I02d*d
            # # # self.dcdah2 = eta*rho*np.sqrt(A/rho)*(self.ur**2 - 2.*eta/k2d*self.ur*self.duzdz)/U2d3/I02d
            # self.dcdav1 = rho*np.sqrt(C/rho)*(1./k2d*self.duzdz)**2/U2d/I02d*d
            # self.dcdbv1 = rho*np.sqrt(L/rho)*( (self.uz)**2 + 2./k2d*self.uz*self.durdz + 4.*eta/k2d*self.ur*self.duzdz+\
            #             (1./k2d*self.durdz)**2 )/U2d/I02d*d
            # self.dcdr1  = -C2d**2/2./U2d/I02d*d*( (self.ur)**2 + (self.uz)**2)  \
            #                 + 1./2./U2d/I02d*d*A/rho*(self.ur)**2 + 1./2./U2d/I02d*d*C/rho*(1./k2d*self.duzdz)**2\
            #                 - 1./U2d/I02d*d/k2d*eta*(A/rho-2.*L/rho)*self.ur*self.duzdz \
            #                 + 1./2./U2d/I02d*d*L/rho*( (self.uz)**2 + 2./k2d*self.uz*self.durdz+ (1./k2d*self.durdz)**2)
            # self.dcdn1  = -1./U2d/I02d*d/k2d*F/eta*self.duzdz*self.ur
    
    def psv_azi_perturb(self, baz, theta4=False):
        """
        get azimuthal perturbation for tilted model
        """
        nfval   = self.freq.size
        Bc2d    = np.tile(self.BcArr, (nfval,1)); Bs2d    = np.tile(self.BsArr, (nfval,1))
        Gc2d    = np.tile(self.GcArr, (nfval,1)); Gs2d    = np.tile(self.GsArr, (nfval,1))
        Hc2d    = np.tile(self.HcArr, (nfval,1)); Hs2d    = np.tile(self.HsArr, (nfval,1))
        Cc2d    = np.tile(self.CcArr, (nfval,1)); Cs2d    = np.tile(self.CsArr, (nfval,1))
        dCrAA   = np.zeros(nfval, np.float32)
        dCrAA   += np.sum( (Bc2d * np.cos(2.*baz/180.*np.pi) + Bs2d * np.sin(2.*baz/180.*np.pi)) * self.dcdA, axis= 1)
        dCrAA   += np.sum( (Gc2d * np.cos(2.*baz/180.*np.pi) + Gs2d * np.sin(2.*baz/180.*np.pi)) * self.dcdL, axis= 1)
        dCrAA   += np.sum( (Hc2d * np.cos(2.*baz/180.*np.pi) + Hs2d * np.sin(2.*baz/180.*np.pi)) * self.dcdF, axis= 1)
        if theta4:
            # dCrAA   += np.sum( (Cc2d * np.cos(4.*baz/180.*np.pi) + Cs2d * np.sin(4.*baz/180.*np.pi)) * self.dcdA, axis= 1)
            dCrAA   += np.sum( (Cc2d * np.cos(4.*baz/180.*np.pi) + Cs2d * np.sin(4.*baz/180.*np.pi)) * self.dcdA, axis= 1)
            
        self.VphA=self.Vph + dCrAA
        return
    
    def psv_perturb_disp_vel(self, insolver):
        """
        Get dispersion curves for perturbed model, using velocity kernels
        """
        nfval   = self.freq.size
        dav     = np.tile( insolver.avArr- self.avArr, (nfval,1))
        dah     = np.tile( insolver.ahArr- self.ahArr, (nfval,1))
        dbv     = np.tile( insolver.bvArr- self.bvArr, (nfval,1))
        dn      = np.tile( insolver.nArr- self.nArr  , (nfval,1)) 
        dr      = np.tile( insolver.rhoArr- self.rhoArr, (nfval,1))
        dc      = np.zeros(nfval, np.float32)
        if (np.nonzero(dav)[0]).size!=0:
            print 'dav'
            dc  += np.sum( dav*self.dcdav, axis=1)
        if (np.nonzero(dah)[0]).size!=0:
            print 'dah'
            dc  += np.sum( dah*self.dcdah, axis=1)
        if (np.nonzero(dbv)[0]).size!=0:
            print 'dbv'
            dc  += np.sum( dbv*self.dcdbv, axis=1)
        if (np.nonzero(dr)[0]).size!=0:
            print 'dr'
            dc  += np.sum( dr*self.dcdr, axis=1)
        if (np.nonzero(dn)[0]).size!=0:
            print 'dn'
            dc  += np.sum( dn*self.dcdn, axis=1)
        self.Vph_pre    = self.Vph + dc
        
    def psv_perturb_disp_love(self, insolver):
        """
        Get dispersion curves for perturbed model, using Love kernels
        """
        nfval   = self.freq.size
        dA      = np.tile( insolver.AArr- self.AArr, (nfval,1))
        dC      = np.tile( insolver.CArr- self.CArr, (nfval,1))
        dF      = np.tile( insolver.FArr- self.FArr, (nfval,1))
        dL      = np.tile( insolver.LArr- self.LArr, (nfval,1))
        dN      = np.tile( insolver.NArr- self.NArr  , (nfval,1)) 
        dr      = np.tile( insolver.rhoArr- self.rhoArr, (nfval,1))
        dc      = np.zeros(nfval, np.float32)
        if (np.nonzero(dA)[0]).size!=0:
            print 'dA'
            dc  += np.sum( dA*self.dcdA, axis=1)
        if (np.nonzero(dC)[0]).size!=0:
            print 'dC'
            dc  += np.sum( dC*self.dcdC, axis=1)
        if (np.nonzero(dF)[0]).size!=0:
            print 'dF'
            dc  += np.sum( dF*self.dcdF, axis=1)
        if (np.nonzero(dL)[0]).size!=0:
            print 'dL'
            dc  += np.sum( dL*self.dcdL, axis=1)
        if (np.nonzero(dr)[0]).size!=0:
            print 'dr'
            dc  += np.sum( dr*self.dcdrl, axis=1)
        self.Vph_pre2    = self.Vph + dc
        
    def solve_SH(self):
        """
        Solve SH motion dispersion curves, eigenfunctions and sensitivity kernels
        ===================================================================================================
        ::: Output :::
        : Model :
        dArr, AArr, CArr, FArr, LArr, NArr, rhoArr  - layerized model
        dsph, Asph, Csph, Fsph, Lsph, Nsph, rhosph  - flattening transformed Earth model
        : dispersion :
        Vph, Vgr    - phase/group velocity
        : eigenfunctions :
        ut, tut     - transverse displacement/stress functions
        dutdz       - derivatives of transverse displacement functions
        : velocity/density sensitivity kernels :
        dcdah/dcdav - vph/vpv kernel
        dcdbh/dcdbv - vsh/vsv kernel
        dcdn        - eta kernel
        dcdr        - density kernel
        : Love parameter/density sensitivity kernels, derived from the kernels above using chain rule :
        dcdA, dcdC, dcdF, dcdL, dcdN    - Love parameter kernels
        dcdrl                           - density kernel
        ===================================================================================================
        """
        #- root-finding algorithm using tdisp96, compute phase velocities ------------------------
        if self.model.tilt == 1:
            dArr, rhoArr, AArr, CArr, FArr, LArr, NArr, BcArr, BsArr, GcArr, GsArr, HcArr, HsArr, CcArr, CsArr =\
                        self.model.get_layer_tilt_model(self.dArr, 200, 1.)
            self.tilt   = True
            self.egn96  = True
            self.BcArr  = BcArr
            self.BsArr  = BsArr
            self.GcArr  = GcArr
            self.GsArr  = GsArr
            self.HcArr  = HcArr
            self.HsArr  = HsArr
            self.CcArr  = CcArr
            self.CsArr  = CsArr
        else:
            dArr, rhoArr, AArr, CArr, FArr, LArr, NArr = self.model.get_layer_model(self.dArr, 200, 1.)
        nfval   = self.freq.size
        if self.model.flat == 0:
            iflsph_in = 1
        else:
            iflsph_in = 0
        ilvry   = 1
        c_out,d_out,TA_out,TC_out,TF_out,TL_out,TN_out,TRho_out = tdisp96.disprs(ilvry, 1., nfval, 1, self.verbose, nfval, \
                np.append(self.freq, np.zeros(2049-nfval)), self.cmin,self.cmax, dArr, AArr,CArr,FArr,LArr,NArr,rhoArr, dArr.size,\
                iflsph_in, 0., self.nmodes, 0.5, 0.5)
        self.ilvry  = ilvry
        self.Vph    = c_out[:nfval]
        # Save flat transformed model is spherical flag is on
        if self.model.flat == 0:
            self.dsph   = d_out
            self.Asph   = TA_out
            self.Csph   = TC_out
            self.Lsph   = TL_out
            self.Fsph   = TF_out
            self.Nsph   = TN_out
            self.rhosph = TRho_out
        self.AArr   = AArr
        self.CArr   = CArr
        self.LArr   = LArr
        self.FArr   = FArr
        self.NArr   = NArr
        self.rhoArr = rhoArr
        self.love2vel()
        #- compute eigenfunction/kernels
        if self.egn96:
            hs_in       = 0.
            hr_in       = 0.
            ohr_in      = 0.
            ohs_in      = 0.
            refdep_in   = 0.
            dogam       = False # No attenuation
            nl_in       = dArr.size
            k           = 2.*np.pi/self.Vph/self.T
            k2d         = np.tile(k, (nl_in, 1))
            k2d         = k2d.T
            omega       = 2.*np.pi/self.T
            omega2d     = np.tile(omega, (nl_in, 1))
            omega2d     = omega2d.T
            if self.model.flat == 0:
                d_in    = d_out
                TA_in   = TA_out
                TC_in   = TC_out
                TF_in   = TF_out
                TL_in   = TL_out
                TN_in   = TN_out
                TRho_in = TRho_out
            else:
                d_in    = dArr
                TA_in   = AArr
                TC_in   = CArr
                TF_in   = FArr
                TL_in   = LArr
                TN_in   = NArr
                TRho_in = rhoArr
            qai_in      = np.ones(nl_in)*1.e6
            qbi_in      = np.ones(nl_in)*1.e6
            etapi_in    = np.zeros(nl_in)
            etasi_in    = np.zeros(nl_in)
            frefpi_in   = np.ones(nl_in)
            frefsi_in   = np.ones(nl_in)
            u_out, ut, tut, dcdh,dcdav,dcdah,dcdbv,dcdbh,dcdn,dcdr = tlegn96.tlegn96(hs_in, hr_in, ohr_in, ohs_in,\
                refdep_in, dogam, nl_in, iflsph_in, d_in, TA_in, TC_in, TF_in, TL_in, TN_in, TRho_in, \
                qai_in,qbi_in,etapi_in,etasi_in, frefpi_in, frefsi_in, self.T.size, self.T, self.Vph)
            ######################################################
            # store output
            ######################################################
            self.Vgr    = u_out
            # eigenfunctions
            self.ut     = ut[:nfval,:nl_in]
            self.tut    = tut[:nfval,:nl_in]
            # sensitivity kernels for velocity parameters and density
            self.dcdah  = np.zeros((nfval, nl_in), dtype=np.float32)
            self.dcdav  = np.zeros((nfval, nl_in), dtype=np.float32)
            self.dcdbh  = dcdbh[:nfval,:nl_in]
            self.dcdbv  = dcdbv[:nfval,:nl_in]
            self.dcdr   = dcdr[:nfval,:nl_in]
            self.dcdn   = np.zeros((nfval, nl_in), dtype=np.float32)
            # Love parameters and density in the shape of nfval, nl_in
            # # A           = np.tile(TA_in, (nfval,1))
            # # C           = np.tile(TC_in, (nfval,1))
            # # F           = np.tile(TF_in, (nfval,1))
            # # L           = np.tile(TL_in, (nfval,1))
            # # N           = np.tile(TN_in, (nfval,1))
            # # rho         = np.tile(TRho_in, (nfval,1))
            # # d           = np.tile(d_in, (nfval,1))
            A           = np.tile(AArr, (nfval,1))
            C           = np.tile(CArr, (nfval,1))
            F           = np.tile(FArr, (nfval,1))
            L           = np.tile(LArr, (nfval,1))
            N           = np.tile(NArr, (nfval,1))
            rho         = np.tile(rhoArr, (nfval,1))
            d           = np.tile(dArr, (nfval,1))
            eta         = F/(A-2.*L)
            # derivative of eigenfunctions, $5.8 of R.Herrmann
            self.dutdz  = self.tut/L
            # integrals for the Lagranian
            I0          = np.sum( rho*(self.ut**2)*d, axis=1)
            I1          = np.sum( N*(self.ut**2)*d, axis=1)
            I2          = np.sum( L*(self.dutdz**2)*d, axis=1)
            # # U           = (k*I1)/omega/I0
            I02d        = np.tile(I0, (nl_in, 1))
            I02d        = I02d.T
            U2d         = np.tile(self.Vgr, (nl_in, 1))
            U2d         = U2d.T
            C2d         = np.tile(self.Vph, (nl_in, 1))
            C2d         = C2d.T
            # sensitivity kernels for Love parameters and density, derived from velocity kernels using chain rule
            # NOTE: benchmark of love kernels predict somewhat different perturbed values compared with those using velocity kernel
            # Possible reason:
            #       By using the chain rule, infinitesimal perturbations for all the parameters is assumed.
            #       However, for the benckmark cases, the perturbation can be large.
            #       e.g. ah = sqrt(A/rho), dah = 0.5/sqrt(A*rho)*dA. But, if dA is large, dah = sqrt(A+dA/rho) - sqrt(A/rho) != 0.5/sqrt(A*rho)*dA
            #      
            self.dcdA   = np.zeros((nfval, nl_in), dtype=np.float32)
            self.dcdC   = np.zeros((nfval, nl_in), dtype=np.float32)
            self.dcdF   = np.zeros((nfval, nl_in), dtype=np.float32)
            self.dcdL   = 0.5/np.sqrt(L*rho)*self.dcdbv + 2.*F/((A-2.*L)**2)*self.dcdn
            self.dcdN   = 0.5/np.sqrt(N*rho)*self.dcdbh
            # density kernel
            self.dcdrl  = -0.5*self.dcdah*np.sqrt(A/(rho**3)) - 0.5*self.dcdav*np.sqrt(C/(rho**3))\
                            -0.5*self.dcdbh*np.sqrt(N/(rho**3)) -0.5*self.dcdbv*np.sqrt(L/(rho**3))\
                                + self.dcdr
        #     # U2d3         = np.tile(U, (nl_in, 1))
        #     # U2d3         = U2d3.T
        #     # # # I0          = np.sum( rho*(self.uz**2+self.ur**2), axis=1)
        #     # # # I1          = np.sum( (L*(self.uz**2)+A*(self.ur**2)), axis=1)
        #     # # # I2          = np.sum( (L*self.uz*self.durdz - F*self.ur*self.duzdz), axis=1)
        #     # # # I3          = np.sum( (L*(self.durdz**2)+C*(self.duzdz**2)), axis=1)
        #     # self.I0 = I0
        #     # self.I1 = I1
        #     # self.I2 = I2
        #     # self.I3 = I3
            ##############################################
            # For benchmark purposes
            ##############################################
            # derivative of eigenfunctions
            self.dutdz1 = -(self.ut[:,:-1] - self.ut[:,1:])/self.dArr[0]
            self.Vgr1   = (k*I1)/omega/I0
            # derived kernels from eigenfunctions, p 299 in Herrmann's notes
            self.dcdbh1 = rho*np.sqrt(N/rho)*(self.ut**2)/U2d/I02d*d
            self.dcdbv1 = rho*np.sqrt(L/rho)*(1./k2d*self.dutdz)**2/U2d/I02d*d
            self.dcdr1  = -C2d**2/2./U2d/I02d*d*( (self.ut)**2 )  \
                            + 1./2./U2d/I02d*d*N/rho*(self.ut**2) \
                            + 1./2./U2d/I02d*d*L/rho*(1./k2d*self.dutdz)**2
            
    def sh_azi_perturb(self, baz, theta4=False):
        """
        get azimuthal perturbation for tilted model
        """
        nfval   = self.freq.size
        Gc2d    = np.tile(self.GcArr, (nfval,1)); Gs2d    = np.tile(self.GsArr, (nfval,1))
        Cc2d    = np.tile(self.CcArr, (nfval,1)); Cs2d    = np.tile(self.CsArr, (nfval,1))
        dClAA   = np.zeros(nfval, np.float32)
        dClAA   += np.sum( (- Gc2d * np.cos(2.*baz/180.*np.pi) - Gs2d * np.sin(2.*baz/180.*np.pi)) * self.dcdL, axis= 1)
        if theta4:
            dClAA   += np.sum( (- Cc2d * np.cos(4.*baz/180.*np.pi) - Cs2d * np.sin(4.*baz/180.*np.pi)) * self.dcdN, axis= 1)
        self.VphA=self.Vph + dClAA
        return
        
    def sh_perturb_disp_vel(self, insolver):
        """
        Get dispersion curves for perturbed model, using velocity kernels
        """
        nfval   = self.freq.size
        # dav     = np.tile( insolver.avArr- self.avArr, (nfval,1))
        dbh     = np.tile( insolver.bhArr- self.bhArr, (nfval,1))
        dbv     = np.tile( insolver.bvArr- self.bvArr, (nfval,1))
        # dn      = np.tile( insolver.nArr- self.nArr  , (nfval,1)) 
        dr      = np.tile( insolver.rhoArr- self.rhoArr, (nfval,1))
        dc      = np.zeros(nfval, np.float32)
        if (np.nonzero(dbv)[0]).size!=0:
            print 'dbv'
            dc  += np.sum( dbv*self.dcdbv, axis=1)
        if (np.nonzero(dbh)[0]).size!=0:
            print 'dbh'
            dc  += np.sum( dbh*self.dcdbh, axis=1)
        if (np.nonzero(dr)[0]).size!=0:
            print 'dr'
            dc  += np.sum( dr*self.dcdr, axis=1)
        self.Vph_pre    = self.Vph + dc
        
    def sh_perturb_disp_love(self, insolver):
        """
        Get dispersion curves for perturbed model, using Love kernels
        """
        nfval   = self.freq.size
        dA      = np.tile( insolver.AArr- self.AArr, (nfval,1))
        dC      = np.tile( insolver.CArr- self.CArr, (nfval,1))
        dF      = np.tile( insolver.FArr- self.FArr, (nfval,1))
        dL      = np.tile( insolver.LArr- self.LArr, (nfval,1))
        dN      = np.tile( insolver.NArr- self.NArr  , (nfval,1)) 
        dr      = np.tile( insolver.rhoArr- self.rhoArr, (nfval,1))
        dc      = np.zeros(nfval, np.float32)
        # if (np.nonzero(dA)[0]).size!=0:
        #     print 'dA'
        #     dc  += np.sum( dA*self.dcdA, axis=1)
        # if (np.nonzero(dC)[0]).size!=0:
        #     print 'dC'
        #     dc  += np.sum( dC*self.dcdC, axis=1)
        # if (np.nonzero(dF)[0]).size!=0:
        #     print 'dF'
        #     dc  += np.sum( dF*self.dcdF, axis=1)
        if (np.nonzero(dL)[0]).size!=0:
            print 'dL'
            dc  += np.sum( dL*self.dcdL, axis=1)
        if (np.nonzero(dN)[0]).size!=0:
            print 'dN'
            dc  += np.sum( dN*self.dcdN, axis=1)
        if (np.nonzero(dr)[0]).size!=0:
            print 'dr'
            dc  += np.sum( dr*self.dcdrl, axis=1)
        self.Vph_pre2    = self.Vph + dc
        
        
        
            
            
            
            
            
    
    
    
    
        
            
        
        
        
    
    