# -*- coding: utf-8 -*-
"""
Module for handling 1D velocity model objects.

Numba is used for speeding up of the code.

:Copyright:
    Author: Lili Feng
    Graduate Research Assistant
    CIEI, Department of Physics, University of Colorado Boulder
    email: lili.feng@colorado.edu
"""
import numpy as np
import numba

@numba.jit(numba.float32[:](numba.float32, numba.float32, numba.float32))
def _get_array(xmin, xmax, dx):
    xlst= []
    Nx  = int((xmax - xmin)/dx + 1)
    for i in xrange(Nx): xlst.append(dx*i+xmin)
    return np.array(xlst, dtype=np.float32)

def read_model(model, infname, unit=1000., isotropic=True,
        indz=0, indvpv=1, indvsv=2, indrho=3, indvph=4, indvsh=5, indeta=6, reverse=True):
    """
    Read model in txt format
    ===========================================================================================================
    Input Parameters:
    model                       - input model1d object
    infname                     - input txt file name
    unit                        - unit of input, default = 1000., means input has units of km
    isotropic                   - whether the input is isotrpic or not
    indz, indvpv, indvsv, indrho- column id(index) for depth, vpv, vsv, rho, vph, vsh, eta
    indvph, indvsh, indeta
    reverse                     - revert the arrays or not
    ===========================================================================================================
    """
    inArr   = np.loadtxt(infname, dtype=np.float32)
    z       = inArr[:, indz]
    radius  = (6371.-z)*unit
    rho     = inArr[:, indrho]*unit
    vpv     = inArr[:, indvpv]*unit
    vsv     = inArr[:, indvsv]*unit
    if isotropic:
        vph     = inArr[:, indvpv]*unit
        vsh     = inArr[:, indvsv]*unit
        eta     = np.ones(vph.size, dtype=np.float32)
    else:
        vph     = inArr[:, indvph]*unit
        vsh     = inArr[:, indvsh]*unit
    if reverse:
        vsv     = vsv[::-1]
        vsh     = vsh[::-1]
        vpv     = vpv[::-1]
        vph     = vph[::-1]
        eta     = eta[::-1]
        rho     = rho[::-1]
        radius  = radius[::-1]
    ind     = radius > 3700000.
    vsv     = vsv[ind]
    vsh     = vsh[ind]
    vpv     = vpv[ind]
    vph     = vph[ind]
    eta     = eta[ind]
    rho     = rho[ind]
    radius  = radius[ind]
    model.get_data_vel(vsv, vsh, vpv, vph, eta, rho, radius)
    return model

def read_axisem_bm(model, infname):
    """
    Read 1D block model from AxiSEM
    ===========================================================================================================
    Input Parameters:
    model       - input model1d object
    infname     - input txt file name
    ===========================================================================================================
    """
    with open(infname, 'rb') as f:
        f.readline()
        cline           = f.readline()
        cline           = cline.split()
        if cline[0] != 'NAME':
            raise ValueError('Unexpected header: '+cline[0])
        f.readline()
        cline           = f.readline()
        cline           = cline.split()
        if cline[0] != 'ANISOTROPIC':
            raise ValueError('Unexpected header: '+cline[0])
        anisotropic     = cline[1]
        if anisotropic == 'T': anisotropic = True
        elif anisotropic == 'F': anisotropic = False
        cline           = f.readline()
        cline           = cline.split()
        if cline[0] != 'UNITS':
            raise ValueError('Unexpected header: '+cline[0])
        if cline[1] == 'm': unit = 1.
        elif cline[1] == 'km': unit = 1000.
        cline           = f.readline()
        cline           = cline.split()
        if cline[0] != 'COLUMNS':
            raise ValueError('Unexpected header: '+cline[0])
        ind = {}
        i=0
        for hdrstr in cline[1:]:
            ind[hdrstr] = i
            i   += 1
        ###
        # Read model parameters
        ###
        z0 = 0.
        vsvArr  = np.array([], dtype=np.float32)
        vshArr  = np.array([], dtype=np.float32)
        vpvArr  = np.array([], dtype=np.float32)
        vphArr  = np.array([], dtype=np.float32)
        etaArr  = np.array([], dtype=np.float32)
        rhoArr  = np.array([], dtype=np.float32)
        rArr    = np.array([], dtype=np.float32)
        for line in f.readlines():
            cline   = line.split()
            if cline[0] == '#': continue
            r   = np.float32(cline[ ind['radius'] ])*unit
            vpv = np.float32(cline[ ind['vpv'] ])*unit
            vsv = np.float32(cline[ ind['vsv'] ])*unit
            rho = np.float32(cline[ ind['rho'] ])*unit
            if anisotropic:
                vph = np.float32(cline[ ind['vph'] ])*unit
                vsh = np.float32(cline[ ind['vsh'] ])*unit
                eta = np.float32(cline[ ind['eta'] ])*unit
            else:
                vph = vpv
                vsh = vsv
                eta = np.float32(1.)
            vsvArr  = np.append(vsvArr, vsv)
            vshArr  = np.append(vshArr, vsh)
            vpvArr  = np.append(vpvArr, vpv)
            vphArr  = np.append(vphArr, vph)
            etaArr  = np.append(etaArr, eta)
            rhoArr  = np.append(rhoArr, rho)
            rArr    = np.append(rArr, r)
    # revert array
    vsvArr  = vsvArr[::-1]
    vshArr  = vshArr[::-1]
    vpvArr  = vpvArr[::-1]
    vphArr  = vphArr[::-1]
    etaArr  = etaArr[::-1]
    rhoArr  = rhoArr[::-1]
    rArr    = rArr[::-1]
    # discard Earth core
    ind     = (rArr > 3700000.)
    vsvArr  = vsvArr[ind]
    vshArr  = vshArr[ind]
    vpvArr  = vpvArr[ind]
    vphArr  = vphArr[ind]
    etaArr  = etaArr[ind]
    rhoArr  = rhoArr[ind]
    rArr    = rArr[ind]
    # reassgin data type
    vsvArr  = vsvArr.astype(np.float32)
    vshArr  = vshArr.astype(np.float32)
    vpvArr  = vpvArr.astype(np.float32)
    vphArr  = vphArr.astype(np.float32)
    etaArr  = etaArr.astype(np.float32)
    rhoArr  = rhoArr.astype(np.float32)
    rArr    = rArr.astype(np.float32)
    model.get_data_vel(vsvArr, vshArr, vpvArr, vphArr, etaArr, rhoArr, rArr)        
    return model

spec = [('VsvArr', numba.float32[:]),
        ('VpvArr', numba.float32[:]),
        ('VshArr', numba.float32[:]),
        ('VphArr', numba.float32[:]),
        ('etaArr', numba.float32[:]),
        ('rhoArr', numba.float32[:]),
        ('rhoArrR', numba.float32[:]),
        ('rhoArrL', numba.float32[:]),
        ('rArr', numba.float32[:]),
        ('rArrS', numba.float32[:]),
        ('AArr', numba.float32[:]),
        ('CArr', numba.float32[:]),
        ('LArr', numba.float32[:]),
        ('FArr', numba.float32[:]),
        ('NArr', numba.float32[:]),
        ('AArrR', numba.float32[:]),
        ('CArrR', numba.float32[:]),
        ('LArrR', numba.float32[:]),
        ('FArrR', numba.float32[:]),
        ('NArrR', numba.float32[:]),
        ('AArrL', numba.float32[:]),
        ('CArrL', numba.float32[:]),
        ('LArrL', numba.float32[:]),
        ('FArrL', numba.float32[:]),
        ('NArrL', numba.float32[:]),
        ('dipArr', numba.float32[:]),
        ('strikeArr', numba.float32[:]),
        ('flat', numba.int32),
        ('rmin', numba.float32)
        ]

@numba.jitclass(spec)
class model1d(object):
    """
    Class defining a 1D Earth model
    ===============================================================================
    Parameters:
    VsvArr, VshArr, - Vsv, Vsh, Vpv, Vph velocity (unit - m/s)
    VpvArr, VphArr  
    rhoArr          - density (kg/m^3)
    etaArr          - eta(F/(A-2L)) dimentionless
    AArr, CArr, FArr- Love parameters (unit - Pa)
    LArr, NArr
    flat            - = 0 spherical Earth, = 1 flat Earth (default)
    ===============================================================================
    """
    def __init__(self):
        self.flat   = 1
        #
        return
    
    def get_data_vel(self, vsv, vsh, vpv, vph, eta, rho, radius):
        """
        Get model data given velocity arrays
        """
        self.rArr   = radius
        self.rhoArr = rho
        ###
        # Velocity
        ###
        self.VsvArr = vsv
        self.VshArr = vsh
        self.VpvArr = vpv
        self.VphArr = vph
        self.etaArr = eta
        ###
        # Love parameters
        ###
        self.AArr   = rho * vph**2
        self.CArr   = rho * vpv**2
        self.LArr   = rho * vsv**2
        self.FArr   = eta * (self.AArr - np.float32(2.)* self.LArr)
        self.NArr   = rho * vsh**2
        return
    
    def earth_flattening_kennett(self):
        """
        Kennett 2009, Seismic Wave Propagation in a Stratified Medium, p16
        """
        z           = np.float32(6371000.)*np.log(np.float32(6371000.)/self.rArr)
        ###
        # Velocity
        ###
        self.VsvArr = self.VsvArr*np.float32(6371000.)/self.rArr
        self.VshArr = self.VshArr*np.float32(6371000.)/self.rArr
        self.VpvArr = self.VpvArr*np.float32(6371000.)/self.rArr
        self.VphArr = self.VphArr*np.float32(6371000.)/self.rArr
        self.rhoArr = self.rhoArr*self.rArr/np.float32(6371000.)
        ###
        # Love parameters
        ###
        self.AArr   = self.rhoArr * self.VphArr**2
        self.CArr   = self.rhoArr * self.VpvArr**2
        self.LArr   = self.rhoArr * self.VsvArr**2
        self.FArr   = self.etaArr * (self.AArr - np.float32(2.)* self.LArr)
        self.NArr   = self.rhoArr * self.VshArr**2
        self.rArr   = np.float32(6371000.)-z
        
        return
    
    def earth_flattening(self):
        """
        Perform Earth flattening transformation for P-SV and SH waves
        Recoded from Robert Herrmann's Computer Program in Seismology,
        subrountine sphere in tdisp96.f
        """
        self.flat   = 0
        z           = np.float32(6371000.)*np.log(np.float32(6371000.)/self.rArr)
        tmp         = np.float32(6371000.)/self.rArr
        ###
        # P-SV, Love parameters
        ###
        self.AArrR  = self.AArr * tmp**(np.float32(-0.2750))
        self.CArrR  = self.CArr * tmp**(np.float32(-0.2750))
        self.LArrR  = self.LArr * tmp**(np.float32(-0.2750))
        self.FArrR  = self.FArr * tmp**(np.float32(-0.2750))
        self.NArrR  = self.NArr * tmp**(np.float32(-0.2750))
        self.rhoArrR= self.rhoArr * tmp**(np.float32(-2.2750))
        ###
        # SH, Love parameters
        ###
        self.AArrL  = self.AArr * tmp**(np.float32(-3.))
        self.CArrL  = self.CArr * tmp**(np.float32(-3.))
        self.LArrL  = self.LArr * tmp**(np.float32(-3.))
        self.FArrL  = self.FArr * tmp**(np.float32(-3.))
        self.NArrL  = self.NArr * tmp**(np.float32(-3.))
        self.rhoArrL= self.rhoArr * tmp**(np.float32(-5.))
        ###
        # radius
        ###
        self.rArrS  = np.float32(6371000.)-z
        self.rmin   = self.rArrS[0]
        return
        
    
    def get_radius(self, zmax, dz):
        """
        Get the radius array
        ===============================================================================
        Input Parameters:
        zmax        - maximum depth (unit - km)
        dz          - depth interval (unit - km)
        Output:
        self.rArr   - radius array (unit - m)
        ===============================================================================
        """
        rmin    = 6371.0 - zmax
        # # Nr      = int((6371. - rmin)/dz + 1)
        # # rlst    = []
        # # for i in xrange(Nr): rlst.append(1000.*(i+rmin)*dz)
        # # self.rArr   = np.array(rlst, dtype=np.float32)
        self.rArr   = _get_array(rmin*1000., 6371000., dz*1000.)
        return
    
    def model_prem(self):
        """
        Get 1D PREM model (Dziewonski & Anderson, PEPI 1981) for 
        a radius r in m. The reference frequency is 1 Hz. Crust continued into the ocean.
        """
        ALst=[]; CLst=[]; LLst=[]; FLst=[]; NLst=[]; rhoLst=[]
        vsvLst=[]; vpvLst=[]; vshLst=[]; vphLst=[]; etaLst=[]
        for r in self.rArr:
            #- normalised radius
            x = r / 6371000.0
            #- march through the various depth levels -----------------------------------------------------
            #- upper crust
            if (r >= 6356000.0):
                rho = 2.6
                vpv = 5.8
                vph = vpv
                vsv = 3.2
                vsh = vsv
                eta = 1.0
            #- lower crust
            elif (r >= 6346000.6) & (r < 6356000.0):
                rho = 2.9
                vpv = 6.8
                vph = vpv
                vsv = 3.9
                vsh = vsv
                eta = 1.0
            #- LID
            elif (r >= 6291000.0) & (r < 6346000.6):
                rho = 2.6910 + 0.6924 * x
                vpv = 0.8317 + 7.2180 * x
                vph = 3.5908 + 4.6172 * x
                vsv = 5.8582 - 1.4678 * x
                vsh = -1.0839 + 5.7176 * x
                eta = 3.3687 - 2.4778 * x
            #- LVZ
            elif (r >= 6151000.0) & (r < 6291000.0):
                rho = 2.6910 + 0.6924 * x
                vpv = 0.8317 + 7.2180 * x
                vph = 3.5908 + 4.6172 * x
                vsv = 5.8582 - 1.4678 * x
                vsh = -1.0839 + 5.7176 * x
                eta = 3.3687 - 2.4778 * x
            #- Transition zone 1
            elif (r >= 5971000.0) & (r < 6151000.0):
                rho = 7.1089 - 3.8045 * x
                vpv = 20.3926 - 12.2569 * x
                vph = vpv
                vsv = 8.9496 - 4.4597 * x
                vsh = vsv
                eta = 1.0
            #- Transition zone 2
            elif (r >= 5771000.0) & (r < 5971000.0):
                rho = 11.2494 - 8.0298 * x
                vpv = 39.7027 - 32.6166 * x
                vph = vpv
                vsv = 22.3512 - 18.5856 * x
                vsh = vsv
                eta = 1.0
            #- Transition zone 3
            elif (r >= 5701000.0) & (r < 5771000.0):
                rho = 5.3197 - 1.4836 * x
                vpv = 19.0957 - 9.8672 * x
                vph = vpv
                vsv = 9.9839 - 4.9324 * x
                vsh = vsv
                eta = 1.0
            #- Lower mantle 1
            elif (r >= 5600000.0) & (r < 5701000.0):
                rho = 7.9565 - 6.4761 * x + 5.5283 * x**2 - 3.0807 * x**3
                vpv = 29.2766 - 23.6027 * x + 5.5242 * x**2 - 2.5514 * x**3
                vph = vpv
                vsv = 22.3459 - 17.2473 * x - 2.0834 * x**2 + 0.9783 * x**3
                vsh = vsv
                eta = 1.0 
            #- Lower mantle 2
            elif (r >= 3630000.0) & (r < 5600000.0):
                rho = 7.9565 - 6.4761 * x + 5.5283 * x**2 - 3.0807 * x**3
                vpv = 24.9520 - 40.4673 * x + 51.4832 * x**2 - 26.6419 * x**3
                vph = vpv
                vsv = 11.1671 - 13.7818 * x + 17.4575 * x**2 - 9.2777 * x**3
                vsh = vsv
                eta = 1.0
            #- Lower mantle 3
            elif (r >= 3480000.0) & (r < 3630000.0):
                rho = 7.9565 - 6.4761 * x + 5.5283 * x**2 - 3.0807 * x**3
                vpv = 15.3891 - 5.3181 * x + 5.5242 * x**2 - 2.5514 * x**3
                vph = vpv
                vsv = 6.9254 + 1.4672 * x - 2.0834 * x**2 + 0.9783 * x**3
                vsh = vsv
                eta = 1.0
            #- Outer core
            elif (r >= 1221000.5) & (r < 3480000.0):
                rho = 12.5815 - 1.2638 * x - 3.6426 * x**2 - 5.5281 * x**3
                vpv = 11.0487 - 4.0362 * x + 4.8023 * x**2 - 13.5732 * x**3
                vph = vpv
                vsv = 0.0
                vsh = 0.0
                eta = 1.0
            #- Inner Core
            elif (r >= 0.0) & (r < 1221000.5):
                rho = 13.0885 - 8.8381 * x**2
                vpv = 11.2622 - 6.3640 * x**2
                vph = vpv
                vsv = 3.6678 - 4.4475 * x**2
                vsh = vsv
                eta = 1.0 
            #- convert to elastic parameters --------------------------------------------------------------
            rho = 1000.0 * rho
            vpv = 1000.0 * vpv
            vph = 1000.0 * vph
            vsv = 1000.0 * vsv
            vsh = 1000.0 * vsh
            
            vpvLst.append(vpv)
            vphLst.append(vph)
            vsvLst.append(vsv)
            vshLst.append(vsh)
            etaLst.append(eta)
            
            A = rho * vph**2
            C = rho * vpv**2
            N = rho * vsh**2
            L = rho * vsv**2
            F = eta * (A - 2 * L)
            
            ALst.append(A)
            CLst.append(C)
            LLst.append(L)
            FLst.append(F)
            NLst.append(N)
            rhoLst.append(rho)
        self.AArr   = np.array(ALst, dtype=np.float32)
        self.CArr   = np.array(CLst, dtype=np.float32)
        self.LArr   = np.array(LLst, dtype=np.float32)
        self.FArr   = np.array(FLst, dtype=np.float32)
        self.NArr   = np.array(NLst, dtype=np.float32)
        self.rhoArr = np.array(rhoLst, dtype=np.float32)
        self.VsvArr = np.array(vsvLst, dtype=np.float32)
        self.VshArr = np.array(vshLst, dtype=np.float32)
        self.VpvArr = np.array(vpvLst, dtype=np.float32)
        self.VphArr = np.array(vphLst, dtype=np.float32)
        self.etaArr = np.array(etaLst, dtype=np.float32)
        self.rmin   = self.rArr.min()
        return
    
    def model_ak135_cps(self):
        ak135_cps_arr = np.array([
                0.00000000e+00,   5.80000000e+00,   3.46000000e+00,   2.72000000e+00, 
                2.00000000e+01,   5.80000000e+00,   3.46000000e+00,   2.72000000e+00,
                2.00000000e+01,   6.50000000e+00,   3.85000000e+00,   2.92000000e+00,
                3.50000000e+01,   6.50000000e+00,   3.85000000e+00,   2.92000000e+00,
                3.50000000e+01,   8.04250000e+00,   4.48500000e+00,   3.33260000e+00,
                7.70000000e+01,   8.04250000e+00,   4.48500000e+00,   3.33260000e+00,
                7.70000000e+01,   8.04750000e+00,   4.49500000e+00,   3.35840000e+00,
                1.20000000e+02,   8.04750000e+00,   4.49500000e+00,   3.35840000e+00,
                1.20000000e+02,   8.11230000e+00,   4.50450000e+00,   3.39900000e+00,
                1.65000000e+02,   8.11230000e+00,   4.50450000e+00,   3.39900000e+00,
                1.65000000e+02,   8.23740000e+00,   4.51370000e+00,   3.34770000e+00,
                2.10000000e+02,   8.23740000e+00,   4.51370000e+00,   3.34770000e+00,
                2.10000000e+02,   8.39050000e+00,   4.56340000e+00,   3.34530000e+00,
                2.60000000e+02,   8.39050000e+00,   4.56340000e+00,   3.34530000e+00,
                2.60000000e+02,   8.57260000e+00,   4.65250000e+00,   3.38860000e+00,
                3.10000000e+02,   8.57260000e+00,   4.65250000e+00,   3.38860000e+00,
                3.10000000e+02,   8.75530000e+00,   4.73940000e+00,   3.43440000e+00,
                3.60000000e+02,   8.75530000e+00,   4.73940000e+00,   3.43440000e+00,
                3.60000000e+02,   8.93800000e+00,   4.82630000e+00,   3.48220000e+00,
                4.10000000e+02,   8.93800000e+00,   4.82630000e+00,   3.48220000e+00,
                4.10000000e+02,   9.44330000e+00,   5.13300000e+00,   3.92950000e+00,
                4.60000000e+02,   9.44330000e+00,   5.13300000e+00,   3.92950000e+00,
                4.60000000e+02,   9.61140000e+00,   5.23880000e+00,   3.92530000e+00,
                5.10000000e+02,   9.61140000e+00,   5.23880000e+00,   3.92530000e+00,
                5.10000000e+02,   9.77940000e+00,   5.34500000e+00,   3.92250000e+00,
                5.60000000e+02,   9.77940000e+00,   5.34500000e+00,   3.92250000e+00,
                5.60000000e+02,   9.94730000e+00,   5.45130000e+00,   3.92120000e+00,
                6.10000000e+02,   9.94730000e+00,   5.45130000e+00,   3.92120000e+00,
                6.10000000e+02,   1.01153000e+01,   5.55700000e+00,   3.92040000e+00,
                6.60000000e+02,   1.01153000e+01,   5.55700000e+00,   3.92040000e+00,
                6.60000000e+02,   1.08562000e+01,   6.02460000e+00,   4.26870000e+00,
                7.10000000e+02,   1.08562000e+01,   6.02460000e+00,   4.26870000e+00,
                7.10000000e+02,   1.09883000e+01,   6.14930000e+00,   4.32750000e+00,
                7.60000000e+02,   1.09883000e+01,   6.14930000e+00,   4.32750000e+00,
                7.60000000e+02,   1.10953000e+01,   6.22620000e+00,   4.38420000e+00,
                8.09500000e+02,   1.10953000e+01,   6.22620000e+00,   4.38420000e+00,
                8.09500000e+02,   1.11790000e+01,   6.26110000e+00,   4.43840000e+00,
                8.59000000e+02,   1.11790000e+01,   6.26110000e+00,   4.43840000e+00,
                8.59000000e+02,   1.12646000e+01,   6.29810000e+00,   4.49060000e+00,
                9.08500000e+02,   1.12646000e+01,   6.29810000e+00,   4.49060000e+00,
                9.08500000e+02,   1.13481000e+01,   6.33410000e+00,   4.54080000e+00,
                9.58000000e+02,   1.13481000e+01,   6.33410000e+00,   4.54080000e+00,
                9.58000000e+02,   1.14299000e+01,   6.36890000e+00,   4.57900000e+00,
                1.00750000e+03,   1.14299000e+01,   6.36890000e+00,   4.57900000e+00,
                1.00750000e+03,   1.15097000e+01,   6.40210000e+00,   4.60620000e+00,
                1.05700000e+03,   1.15097000e+01,   6.40210000e+00,   4.60620000e+00,
                1.05700000e+03,   1.15878000e+01,   6.43480000e+00,   4.63330000e+00,
                1.10650000e+03,   1.15878000e+01,   6.43480000e+00,   4.63330000e+00,
                1.10650000e+03,   1.16641000e+01,   6.46680000e+00,   4.66010000e+00,
                1.15600000e+03,   1.16641000e+01,   6.46680000e+00,   4.66010000e+00,
                1.15600000e+03,   1.17393000e+01,   6.49760000e+00,   4.68680000e+00,
                1.20550000e+03,   1.17393000e+01,   6.49760000e+00,   4.68680000e+00,
                1.20550000e+03,   1.18128000e+01,   6.52810000e+00,   4.71340000e+00,
                1.25500000e+03,   1.18128000e+01,   6.52810000e+00,   4.71340000e+00,
                1.25500000e+03,   1.18848000e+01,   6.55790000e+00,   4.73970000e+00,
                1.30450000e+03,   1.18848000e+01,   6.55790000e+00,   4.73970000e+00,
                1.30450000e+03,   1.19549000e+01,   6.58680000e+00,   4.76590000e+00,
                1.35400000e+03,   1.19549000e+01,   6.58680000e+00,   4.76590000e+00,
                1.35400000e+03,   1.20230000e+01,   6.61470000e+00,   4.79200000e+00,
                1.40350000e+03,   1.20230000e+01,   6.61470000e+00,   4.79200000e+00,
                1.40350000e+03,   1.20908000e+01,   6.64190000e+00,   4.81790000e+00,
                1.45300000e+03,   1.20908000e+01,   6.64190000e+00,   4.81790000e+00,
                1.45300000e+03,   1.21579000e+01,   6.66830000e+00,   4.84350000e+00,
                1.50250000e+03,   1.21579000e+01,   6.66830000e+00,   4.84350000e+00,
                1.50250000e+03,   1.22234000e+01,   6.69410000e+00,   4.86890000e+00,
                1.55200000e+03,   1.22234000e+01,   6.69410000e+00,   4.86890000e+00,
                1.55200000e+03,   1.22869000e+01,   6.71960000e+00,   4.89430000e+00,
                1.60150000e+03,   1.22869000e+01,   6.71960000e+00,   4.89430000e+00,
                1.60150000e+03,   1.23496000e+01,   6.74510000e+00,   4.91950000e+00,
                1.65100000e+03,   1.23496000e+01,   6.74510000e+00,   4.91950000e+00,
                1.65100000e+03,   1.24119000e+01,   6.76990000e+00,   4.94450000e+00,
                1.70050000e+03,   1.24119000e+01,   6.76990000e+00,   4.94450000e+00,
                1.70050000e+03,   1.24728000e+01,   6.79380000e+00,   4.96930000e+00,
                1.75000000e+03,   1.24728000e+01,   6.79380000e+00,   4.96930000e+00,
                1.75000000e+03,   1.25333000e+01,   6.81720000e+00,   4.99390000e+00,
                1.79950000e+03,   1.25333000e+01,   6.81720000e+00,   4.99390000e+00,
                1.79950000e+03,   1.25931000e+01,   6.84030000e+00,   5.01840000e+00,
                1.84900000e+03,   1.25931000e+01,   6.84030000e+00,   5.01840000e+00,
                1.84900000e+03,   1.26516000e+01,   6.86300000e+00,   5.04270000e+00,
                1.89850000e+03,   1.26516000e+01,   6.86300000e+00,   5.04270000e+00,
                1.89850000e+03,   1.27095000e+01,   6.88570000e+00,   5.06680000e+00,
                1.94800000e+03,   1.27095000e+01,   6.88570000e+00,   5.06680000e+00,
                1.94800000e+03,   1.27669000e+01,   6.90830000e+00,   5.09080000e+00,
                1.99750000e+03,   1.27669000e+01,   6.90830000e+00,   5.09080000e+00,
                1.99750000e+03,   1.28239000e+01,   6.93050000e+00,   5.11460000e+00,
                2.04700000e+03,   1.28239000e+01,   6.93050000e+00,   5.11460000e+00,
                2.04700000e+03,   1.28808000e+01,   6.95200000e+00,   5.13820000e+00,
                2.09650000e+03,   1.28808000e+01,   6.95200000e+00,   5.13820000e+00,
                2.09650000e+03,   1.29377000e+01,   6.97380000e+00,   5.16160000e+00,
                2.14600000e+03,   1.29377000e+01,   6.97380000e+00,   5.16160000e+00,
                2.14600000e+03,   1.29944000e+01,   6.99600000e+00,   5.18480000e+00,
                2.19550000e+03,   1.29944000e+01,   6.99600000e+00,   5.18480000e+00,
                2.19550000e+03,   1.30505000e+01,   7.01770000e+00,   5.20780000e+00,
                2.24500000e+03,   1.30505000e+01,   7.01770000e+00,   5.20780000e+00,
                2.24500000e+03,   1.31061000e+01,   7.03950000e+00,   5.23060000e+00,
                2.29450000e+03,   1.31061000e+01,   7.03950000e+00,   5.23060000e+00,
                2.29450000e+03,   1.31615000e+01,   7.06130000e+00,   5.25330000e+00,
                2.34400000e+03,   1.31615000e+01,   7.06130000e+00,   5.25330000e+00,
                2.34400000e+03,   1.32179000e+01,   7.08270000e+00,   5.27580000e+00,
                2.39350000e+03,   1.32179000e+01,   7.08270000e+00,   5.27580000e+00,
                2.39350000e+03,   1.32740000e+01,   7.10380000e+00,   5.29810000e+00,
                2.44300000e+03,   1.32740000e+01,   7.10380000e+00,   5.29810000e+00,
                2.44300000e+03,   1.33300000e+01,   7.12560000e+00,   5.32020000e+00,
                2.49250000e+03,   1.33300000e+01,   7.12560000e+00,   5.32020000e+00,
                2.49250000e+03,   1.33869000e+01,   7.14760000e+00,   5.34220000e+00,
                2.54200000e+03,   1.33869000e+01,   7.14760000e+00,   5.34220000e+00,
                2.54200000e+03,   1.34448000e+01,   7.16940000e+00,   5.36390000e+00,
                2.59150000e+03,   1.34448000e+01,   7.16940000e+00,   5.36390000e+00,
                2.59150000e+03,   1.35025000e+01,   7.19170000e+00,   5.38550000e+00,
                2.64000000e+03,   1.35025000e+01,   7.19170000e+00,   5.38550000e+00,
                2.64000000e+03,   1.35604000e+01,   7.21420000e+00,   5.40690000e+00,
                2.69000000e+03,   1.35604000e+01,   7.21420000e+00,   5.40690000e+00,
                2.69000000e+03,   1.36198000e+01,   7.23690000e+00,   5.42820000e+00,
                2.74000000e+03,   1.36198000e+01,   7.23690000e+00,   5.42820000e+00,
                2.74000000e+03,   1.36516000e+01,   7.25390000e+00,   5.70650000e+00,
                2.78966990e+03,   1.36516000e+01,   7.25390000e+00,   5.70650000e+00,
                2.78966990e+03,   1.36552000e+01,   7.26460000e+00,   5.73270000e+00,
                2.83933010e+03,   1.36552000e+01,   7.26460000e+00,   5.73270000e+00,
                2.83933010e+03,   1.36585000e+01,   7.27580000e+00,   5.75900000e+00,
                2.89150000e+03,   1.36585000e+01,   7.27580000e+00,   5.75900000e+00,
                2.89150000e+03,   8.01910000e+00,   0.00000000e+00,   9.95430000e+00,
                2.93933010e+03,   8.01910000e+00,   0.00000000e+00,   9.95430000e+00,
                2.93933010e+03,   8.08300000e+00,   0.00000000e+00,   1.00332000e+01,
                2.98965990e+03,   8.08300000e+00,   0.00000000e+00,   1.00332000e+01,
                2.98965990e+03,   8.17450000e+00,   0.00000000e+00,   1.01103000e+01,
                3.03999000e+03,   8.17450000e+00,   0.00000000e+00,   1.01103000e+01,
                3.03999000e+03,   8.26650000e+00,   0.00000000e+00,   1.01859000e+01,
                3.09032010e+03,   8.26650000e+00,   0.00000000e+00,   1.01859000e+01,
                3.09032010e+03,   8.35590000e+00,   0.00000000e+00,   1.02598000e+01,
                3.14065990e+03,   8.35590000e+00,   0.00000000e+00,   1.02598000e+01,
                3.14065990e+03,   8.44290000e+00,   0.00000000e+00,   1.03321000e+01,
                3.19099000e+03,   8.44290000e+00,   0.00000000e+00,   1.03321000e+01,
                3.19099000e+03,   8.52740000e+00,   0.00000000e+00,   1.04029000e+01,
                3.24132010e+03,   8.52740000e+00,   0.00000000e+00,   1.04029000e+01,
                3.24132010e+03,   8.60920000e+00,   0.00000000e+00,   1.04720000e+01,
                3.29164990e+03,   8.60920000e+00,   0.00000000e+00,   1.04720000e+01,
                3.29164990e+03,   8.68880000e+00,   0.00000000e+00,   1.05396000e+01,
                3.34198000e+03,   8.68880000e+00,   0.00000000e+00,   1.05396000e+01,
                3.34198000e+03,   8.76580000e+00,   0.00000000e+00,   1.06058000e+01,
                3.39231010e+03,   8.76580000e+00,   0.00000000e+00,   1.06058000e+01,
                3.39231010e+03,   8.83970000e+00,   0.00000000e+00,   1.06704000e+01,
                3.44263990e+03,   8.83970000e+00,   0.00000000e+00,   1.06704000e+01,
                3.44263990e+03,   8.91100000e+00,   0.00000000e+00,   1.07335000e+01,
                3.49297000e+03,   8.91100000e+00,   0.00000000e+00,   1.07335000e+01,
                3.49297000e+03,   8.97980000e+00,   0.00000000e+00,   1.07952000e+01,
                3.54330010e+03,   8.97980000e+00,   0.00000000e+00,   1.07952000e+01,
                3.54330010e+03,   9.04640000e+00,   0.00000000e+00,   1.08554000e+01,
                3.59363990e+03,   9.04640000e+00,   0.00000000e+00,   1.08554000e+01,
                3.59363990e+03,   9.11080000e+00,   0.00000000e+00,   1.09143000e+01,
                3.64397000e+03,   9.11080000e+00,   0.00000000e+00,   1.09143000e+01,
                3.64397000e+03,   9.17330000e+00,   0.00000000e+00,   1.09718000e+01,
                3.69430010e+03,   9.17330000e+00,   0.00000000e+00,   1.09718000e+01,
                3.69430010e+03,   9.23370000e+00,   0.00000000e+00,   1.10278000e+01,
                3.74462990e+03,   9.23370000e+00,   0.00000000e+00,   1.10278000e+01,
                3.74462990e+03,   9.29190000e+00,   0.00000000e+00,   1.10825000e+01,
                3.79496000e+03,   9.29190000e+00,   0.00000000e+00,   1.10825000e+01,
                3.79496000e+03,   9.34820000e+00,   0.00000000e+00,   1.11359000e+01,
                3.84529010e+03,   9.34820000e+00,   0.00000000e+00,   1.11359000e+01,
                3.84529010e+03,   9.40280000e+00,   0.00000000e+00,   1.11880000e+01,
                3.89562020e+03,   9.40280000e+00,   0.00000000e+00,   1.11880000e+01,
                3.89562020e+03,   9.45550000e+00,   0.00000000e+00,   1.12388000e+01,
                3.94595000e+03,   9.45550000e+00,   0.00000000e+00,   1.12388000e+01,
                3.94595000e+03,   9.50590000e+00,   0.00000000e+00,   1.12883000e+01,
                3.99628010e+03,   9.50590000e+00,   0.00000000e+00,   1.12883000e+01,
                3.99628010e+03,   9.55410000e+00,   0.00000000e+00,   1.13365000e+01,
                4.04662020e+03,   9.55410000e+00,   0.00000000e+00,   1.13365000e+01,
                4.04662020e+03,   9.60040000e+00,   0.00000000e+00,   1.13836000e+01,
                4.09695030e+03,   9.60040000e+00,   0.00000000e+00,   1.13836000e+01,
                4.09695030e+03,   9.64520000e+00,   0.00000000e+00,   1.14295000e+01,
                4.14727990e+03,   9.64520000e+00,   0.00000000e+00,   1.14295000e+01,
                4.14727990e+03,   9.68860000e+00,   0.00000000e+00,   1.14741000e+01,
                4.19761000e+03,   9.68860000e+00,   0.00000000e+00,   1.14741000e+01,
                4.19761000e+03,   9.73060000e+00,   0.00000000e+00,   1.15176000e+01,
                4.24794010e+03,   9.73060000e+00,   0.00000000e+00,   1.15176000e+01,
                4.24794010e+03,   9.77130000e+00,   0.00000000e+00,   1.15600000e+01,
                4.29827020e+03,   9.77130000e+00,   0.00000000e+00,   1.15600000e+01,
                4.29827020e+03,   9.81090000e+00,   0.00000000e+00,   1.16012000e+01,
                4.34860030e+03,   9.81090000e+00,   0.00000000e+00,   1.16012000e+01,
                4.34860030e+03,   9.84930000e+00,   0.00000000e+00,   1.16414000e+01,
                4.39893040e+03,   9.84930000e+00,   0.00000000e+00,   1.16414000e+01,
                4.39893040e+03,   9.88660000e+00,   0.00000000e+00,   1.16805000e+01,
                4.44926000e+03,   9.88660000e+00,   0.00000000e+00,   1.16805000e+01,
                4.44926000e+03,   9.92300000e+00,   0.00000000e+00,   1.17185000e+01,
                4.49960030e+03,   9.92300000e+00,   0.00000000e+00,   1.17185000e+01,
                4.49960030e+03,   9.95850000e+00,   0.00000000e+00,   1.17555000e+01,
                4.54993040e+03,   9.95850000e+00,   0.00000000e+00,   1.17555000e+01,
                4.54993040e+03,   9.99320000e+00,   0.00000000e+00,   1.17915000e+01,
                4.60026000e+03,   9.99320000e+00,   0.00000000e+00,   1.17915000e+01,
                4.60026000e+03,   1.00271000e+01,   0.00000000e+00,   1.18265000e+01,
                4.65059010e+03,   1.00271000e+01,   0.00000000e+00,   1.18265000e+01,
                4.65059010e+03,   1.00603000e+01,   0.00000000e+00,   1.18605000e+01,
                4.70092020e+03,   1.00603000e+01,   0.00000000e+00,   1.18605000e+01,
                4.70092020e+03,   1.00931000e+01,   0.00000000e+00,   1.18935000e+01,
                4.75125030e+03,   1.00931000e+01,   0.00000000e+00,   1.18935000e+01,
                4.75125030e+03,   1.01255000e+01,   0.00000000e+00,   1.19256000e+01,
                4.80158040e+03,   1.01255000e+01,   0.00000000e+00,   1.19256000e+01,
                4.80158040e+03,   1.01577000e+01,   0.00000000e+00,   1.19568000e+01,
                4.85191050e+03,   1.01577000e+01,   0.00000000e+00,   1.19568000e+01,
                4.85191050e+03,   1.01894000e+01,   0.00000000e+00,   1.19862000e+01,
                4.90224060e+03,   1.01894000e+01,   0.00000000e+00,   1.19862000e+01,
                4.90224060e+03,   1.02189000e+01,   0.00000000e+00,   1.20156000e+01,
                4.95258040e+03,   1.02189000e+01,   0.00000000e+00,   1.20156000e+01,
                4.95258040e+03,   1.02447000e+01,   0.00000000e+00,   1.20452000e+01,
                5.00291050e+03,   1.02447000e+01,   0.00000000e+00,   1.20452000e+01,
                5.00291050e+03,   1.02655000e+01,   0.00000000e+00,   1.20730000e+01,
                5.05324060e+03,   1.02655000e+01,   0.00000000e+00,   1.20730000e+01,
                5.05324060e+03,   1.02799000e+01,   0.00000000e+00,   1.21000000e+01,
                5.10357020e+03,   1.02872000e+01,   0.00000000e+00,   1.21262000e+01], dtype = np.float32)
        data    = ak135_cps_arr.reshape( ak135_cps_arr.size/4, 4)
        # assign model parameters
        z       = data[:, 0]
        radius  = (np.float32(6371.)-z)*np.float32(1000.)
        rho     = data[:, 3]*np.float32(1000.)
        vpv     = data[:, 1]*np.float32(1000.)
        vsv     = data[:, 2]*np.float32(1000.)
        vph     = vpv
        vsh     = vsv
        eta     = np.ones(vph.size, dtype=np.float32)
        # revert the array
        vsv     = vsv[::-1]
        vsh     = vsh[::-1]
        vpv     = vpv[::-1]
        vph     = vph[::-1]
        eta     = eta[::-1]
        rho     = rho[::-1]
        radius  = radius[::-1]
        # discard data with zero shear wave velocity
        ind     = (vsv!=0.)
        vsv     = vsv[ind]
        vsh     = vsh[ind]
        vpv     = vpv[ind]
        vph     = vph[ind]
        eta     = eta[ind]
        rho     = rho[ind]
        radius  = radius[ind]        
        # # # # assign zero shear velocity to 0.001 m/s
        # # # ind     = (vsv==0.)
        # # # vsv[ind]= np.float32(0.01)
        # # # vsh[ind]= np.float32(0.01)
        self.get_data_vel(vsv, vsh, vpv, vph, eta, rho, radius)
        self.rmin   = self.rArr.min()
        return
    
    def trim(self, zmax):
        """
        Trim the model given a maximum depth
        """
        rmin        = (np.float32(6371.)-zmax)*np.float32(1000.)
        ind         = self.rArr > rmin
        Nt          = (self.rArr[ind]).size
        N           = self.rArr.size
        ind[N-Nt-1] = True
        self.rArr   = self.rArr[ind]
        self.rhoArr = self.rhoArr[ind]
        # Love parameters
        self.AArr   = self.AArr[ind]
        self.CArr   = self.CArr[ind]
        self.FArr   = self.FArr[ind]
        self.LArr   = self.LArr[ind]
        self.NArr   = self.NArr[ind]
        # velocity
        self.VsvArr = self.VsvArr[ind]
        self.VshArr = self.VshArr[ind]
        self.VpvArr = self.VpvArr[ind]
        self.VphArr = self.VphArr[ind]
        self.etaArr = self.etaArr[ind]
        if self.flat == 0:
            self.earth_flattening()
        return
        
        
    def get_ind_Love_parameters(self, i):
        """
        Return Love paramaters and density given an index
        """
        return self.rhoArr[i], self.AArr[i], self.CArr[i], self.FArr[i], self.LArr[i], self.NArr[i]
    
    def get_ind_Love_parameters_PSV(self, i):
        """
        Return Love paramaters and density given an index, for P-SV waves
        """
        if self.flat == 0:
            return self.rhoArrR[i], self.AArrR[i], self.CArrR[i], self.FArrR[i], self.LArrR[i], self.NArrR[i]
        else:
            return self.get_ind_Love_parameters(i)
    
    def get_ind_Love_parameters_SH(self, i):
        """
        Return Love paramaters and density given an index, for SH waves
        """
        if self.flat == 0:
            return self.rhoArrL[i], self.AArrL[i], self.CArrL[i], self.FArrL[i], self.LArrL[i], self.NArrL[i]
        else:
            return self.get_ind_Love_parameters(i)
    
    def get_r_love_parameters(self, r):
        """
        Return Love paramaters and density given a radius
        """
        if r < self.rmin: raise ValueError('Required radius is out of the model range!')
        ind_l = -1; ind_r = -1
        for _r in self.rArr:
            if r < _r:
                ind_r += 1
                break
            if r == _r:
                ind_l += 1
                ind_r += 1
                break
            ind_l += 1
            ind_r += 1
        if ind_l == ind_r:
            if self.rArr[ind_l] == self.rArr[ind_l+1]: ind_l+=1
            return self.rhoArr[ind_l], self.AArr[ind_l], self.CArr[ind_l], \
                self.FArr[ind_l], self.LArr[ind_l], self.NArr[ind_l]
        r_left  = self.rArr[ind_l]
        r_right = self.rArr[ind_r]
        rhol    = self.rhoArr[ind_l]; rhor  =   self.rhoArr[ind_r]
        rho     = rhol + (r - r_left)*(rhor-rhol)/(r_right - r_left)
        Al      = self.AArr[ind_l]; Ar  =   self.AArr[ind_r]
        A       = Al + (r - r_left)*(Ar-Al)/(r_right - r_left)
        Cl      = self.CArr[ind_l]; Cr  =   self.CArr[ind_r]
        C       = Cl + (r - r_left)*(Cr-Cl)/(r_right - r_left)    
        Fl      = self.FArr[ind_l]; Fr  =   self.FArr[ind_r]
        F       = Fl + (r - r_left)*(Fr-Fl)/(r_right - r_left)
        Ll      = self.LArr[ind_l]; Lr  =   self.LArr[ind_r]
        L       = Ll + (r - r_left)*(Lr-Ll)/(r_right - r_left)
        Nl      = self.NArr[ind_l]; Nr  =   self.NArr[ind_r]
        N       = Nl + (r - r_left)*(Nr-Nl)/(r_right - r_left)
        return rho, A, C, F, L, N
    
    def get_r_love_parameters_PSV(self, r):
        """
        Return Love paramaters and density given a radius, for P-SV waves
        """
        if r < self.rmin: raise ValueError('Required radius is out of the model range!')
        if self.flat == 0:
            ind_l = -1; ind_r = -1
            for _r in self.rArrS:
                if r < _r:
                    ind_r += 1
                    break
                if r == _r:
                    ind_l += 1
                    ind_r += 1
                    break
                ind_l += 1
                ind_r += 1
            if ind_l == ind_r:
                if self.rArrS[ind_l] == self.rArrS[ind_l+1]: ind_l+=1
                return self.rhoArrR[ind_l], self.AArrR[ind_l], self.CArrR[ind_l], \
                    self.FArrR[ind_l], self.LArrR[ind_l], self.NArrR[ind_l]
            r_left  = self.rArrS[ind_l]
            r_right = self.rArrS[ind_r]
            rhol    = self.rhoArrR[ind_l]; rhor  =   self.rhoArrR[ind_r]
            rho     = rhol + (r - r_left)*(rhor-rhol)/(r_right - r_left)
            Al      = self.AArrR[ind_l]; Ar  =   self.AArrR[ind_r]
            A       = Al + (r - r_left)*(Ar-Al)/(r_right - r_left)
            Cl      = self.CArrR[ind_l]; Cr  =   self.CArrR[ind_r]
            C       = Cl + (r - r_left)*(Cr-Cl)/(r_right - r_left)    
            Fl      = self.FArrR[ind_l]; Fr  =   self.FArrR[ind_r]
            F       = Fl + (r - r_left)*(Fr-Fl)/(r_right - r_left)
            Ll      = self.LArrR[ind_l]; Lr  =   self.LArrR[ind_r]
            L       = Ll + (r - r_left)*(Lr-Ll)/(r_right - r_left)
            Nl      = self.NArrR[ind_l]; Nr  =   self.NArrR[ind_r]
            N       = Nl + (r - r_left)*(Nr-Nl)/(r_right - r_left)
            return rho, A, C, F, L, N
        else:
            return self.get_r_love_parameters(r)
    
    def get_r_love_parameters_SH(self, r):
        """
        Return Love paramaters and density given a radius, for SH waves
        """
        if r < self.rmin: raise ValueError('Required radius is out of the model range!')
        if self.flat == 0:
            ind_l = -1; ind_r = -1
            for _r in self.rArrS:
                if r < _r:
                    ind_r += 1
                    break
                if r == _r:
                    ind_l += 1
                    ind_r += 1
                    break
                ind_l += 1
                ind_r += 1
            if ind_l == ind_r:
                if self.rArrS[ind_l] == self.rArrS[ind_l+1]: ind_l+=1
                return self.rhoArrL[ind_l], self.AArrL[ind_l], self.CArrL[ind_l], \
                    self.FArrL[ind_l], self.LArrL[ind_l], self.NArrL[ind_l]
            r_left  = self.rArrS[ind_l]
            r_right = self.rArrS[ind_r]
            rhol    = self.rhoArrL[ind_l]; rhor  =   self.rhoArrL[ind_r]
            rho     = rhol + (r - r_left)*(rhor-rhol)/(r_right - r_left)
            Al      = self.AArrL[ind_l]; Ar  =   self.AArrL[ind_r]
            A       = Al + (r - r_left)*(Ar-Al)/(r_right - r_left)
            Cl      = self.CArrL[ind_l]; Cr  =   self.CArrL[ind_r]
            C       = Cl + (r - r_left)*(Cr-Cl)/(r_right - r_left)    
            Fl      = self.FArrL[ind_l]; Fr  =   self.FArrL[ind_r]
            F       = Fl + (r - r_left)*(Fr-Fl)/(r_right - r_left)
            Ll      = self.LArrL[ind_l]; Lr  =   self.LArrL[ind_r]
            L       = Ll + (r - r_left)*(Lr-Ll)/(r_right - r_left)
            Nl      = self.NArrL[ind_l]; Nr  =   self.NArrL[ind_r]
            N       = Nl + (r - r_left)*(Nr-Nl)/(r_right - r_left)
            return rho, A, C, F, L, N
        else:
            return self.get_r_love_parameters(r)
        
    def get_cps_model(self, dArr, nl, dh):
        """
        Get the layerized model for CPS
        Note that the unit is different from the default unit of the class
        ===================================================================
        ::: Input Parameters :::
        dArr            - numpy array of layer thickness (unit - km)
        nl              - number of layers
        dh              - thickness of each layer (unit - km)
        nl and dh will be used if and only if dArr.size = 0
        ::: Output :::
        dArr            - layer thickness array (unit - km)
        rhoArr          - density array (unit - g/cm^3)
        AArr, CArr, FArr- Love parameters (unit - GPa)
        LArr, NArr
        ===================================================================
        """
        if dArr.size==0:
            dh  *= 1000. 
            dArr= np.ones(nl, dtype = np.float32)*np.float32(dh)
        else:
            dArr    *= 1000.
            nl      = dArr.size
        ALst=[]; CLst=[]; LLst=[]; FLst=[]; NLst=[]; rhoLst=[]
        z0   = 0.
        z1   = dArr[0]
        for i in xrange(nl):
            r0  = 6371000.-z0
            r1  = 6371000.-z1
            rho0, A0, C0, F0, L0, N0  = self.get_r_love_parameters(r0)
            rho1, A1, C1, F1, L1, N1  = self.get_r_love_parameters(r1)
            # density is converted from kg/m^3 to g/cm^3
            rho = (rho0+rho1)/np.float32(1e3)/2.
            rhoLst.append(rho)
            # Love parameters are converted from Pa to GPa
            A   = (A0+A1)/1.e9/2.
            ALst.append(A)
            C   = (C0+C1)/1.e9/2.
            CLst.append(C)
            F   = (F0+F1)/1.e9/2.
            FLst.append(F)
            L   = (L0+L1)/1.e9/2.
            LLst.append(L)
            N   = (N0+N1)/1.e9/2.
            NLst.append(N)
            z0  += dArr[i]
            z1  += dArr[i]
        # layer thickness is converted from m to km
        dArr    /=1000.
        rhoArr  = np.array(rhoLst, dtype=np.float32)
        AArr    = np.array(ALst, dtype=np.float32)
        CArr    = np.array(CLst, dtype=np.float32)
        FArr    = np.array(FLst, dtype=np.float32)
        LArr    = np.array(LLst, dtype=np.float32)
        NArr    = np.array(NLst, dtype=np.float32)
        return dArr, rhoArr, AArr, CArr, FArr, LArr, NArr
    
        
    
    
        
