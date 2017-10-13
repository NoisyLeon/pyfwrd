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
import matplotlib.pyplot as plt
import tdisp96

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
    ::: input parameters :::
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
    ::: input parameters :::
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

def layer_aniprop_model_sph(inmodel, dArr, nl, dh, ilvry):
    """
    Get the flattening transformed layerized model
    ===================================================================
    ::: input parameters :::
    inmodel - input model1d object
    dArr    - numpy array of layer thickness (unit - km)
    nl      - number of layers 
    dh      - thickness of each layer (unit - km)
    ilvry   - 1 - Love wave; 2 - Rayleigh wave
    ===================================================================
    """
    # get the layerized model
    dArr, rhoArr, AArr, CArr, FArr, LArr, NArr = inmodel.get_layer_model(dArr, nl, dh)
    verby       = False
    nl_in       = rhoArr.size
    # Earth flattening transformation
    d_out,TA_out,TC_out,TF_out,TL_out,TN_out,TRho_out=tdisp96.flat2sphere(ilvry, 1., 1, 1, verby, \
                1, np.append(10., np.zeros(2048)), -1.,-1., dArr,AArr,CArr,FArr,LArr,NArr,rhoArr,nl_in, 0., 1, 0.5, 0.5)
    # convert layrized model to input model for aniprop
    zArr    = d_out.cumsum()   
    a       = (3.*TA_out + 3.*TC_out + 2.*TF_out + 4.*TL_out)/8.
    b       = 0.5*(TC_out - TA_out)
    c       = (TA_out + TC_out - 2.*TF_out - 4.*TL_out)/8.
    d       = 0.5*(TN_out + TL_out)
    e       = 0.5*(TL_out - TN_out)
    if np.any(e>0.000001):
        raise ValueError('aniprop does not accept L > N !')
    # total number of layers
    nl      = a.size
    # depth to each interface
    z       = zArr
    # model parameters
    rho     = TRho_out
    vp0     = np.sqrt(a/TRho_out)
    vp2     = b/a
    vp4     = c/a
    vs0     = np.sqrt(d/TRho_out)
    vs2     = e/d
    return z, rho, vp0, vp2, vp4, vs0, vs2


def layer_aniprop_model_sph_0(inmodel, dArr, nl, dh, ilvry):
    """
    Get the flattening transformed layerized model, add 0 to top, deprecated
    ===================================================================
    ::: input parameters ::::
    inmodel - input model1d object
    dArr    - numpy array of layer thickness (unit - km)
    nl      - number of layers 
    dh      - thickness of each layer (unit - km)
    ilvry   - 1 - Love wave; 2 - Rayleigh wave
    ===================================================================
    """
    # get the layerized model
    dArr, rhoArr, AArr, CArr, FArr, LArr, NArr = inmodel.get_layer_model(dArr, nl, dh)
    verby       = False
    nl_in       = rhoArr.size
    # Earth flattening transformation
    d_out,TA_out,TC_out,TF_out,TL_out,TN_out,TRho_out=tdisp96.flat2sphere(ilvry, 1., 1, 1, verby, \
                1, np.append(10., np.zeros(2048)), -1.,-1., dArr,AArr,CArr,FArr,LArr,NArr,rhoArr,nl_in, 0., 1, 0.5, 0.5)
    # convert layrized model to input model for aniprop
    zArr    = d_out.cumsum()   
    a       = (3.*TA_out + 3.*TC_out + 2.*TF_out + 4.*TL_out)/8.
    b       = 0.5*(TC_out - TA_out)
    c       = (TA_out + TC_out - 2.*TF_out - 4.*TL_out)/8.
    d       = 0.5*(TN_out + TL_out)
    e       = 0.5*(TL_out - TN_out)
    if np.any(e>0.0001):
        raise ValueError('aniprop does not accept L > N !')
    # total number of layers
    nl      = a.size
    # depth to each interface
    z       = np.zeros(nl+1, dtype=np.float32)
    z[1:]   = zArr
    z[0]    = 0.
    # model parameters
    rho     = np.zeros(nl+1, dtype=np.float32)
    rho[1:] = TRho_out
    rho[0]  = rho[1]
    
    vp0     = np.zeros(nl+1, dtype=np.float32)
    vp0[1:] = np.sqrt(a/TRho_out)
    vp0[0]  = vp0[1]
    
    vp2     = np.zeros(nl+1, dtype=np.float32)
    vp2[1:] = b/a
    vp2[0]  = vp2[1]
    
    vp4     = np.zeros(nl+1, dtype=np.float32)
    vp4[1:] = c/a
    vp4[0]  = vp4[1]
    
    vs0     = np.zeros(nl+1, dtype=np.float32)
    vs0[1:] = np.sqrt(d/TRho_out)
    vs0[0]  = vs0[1]
    
    vs2     = np.zeros(nl+1, dtype=np.float32)
    vs2[1:] = e/d
    vs2[0]  = vs2[1]
    return z, rho, vp0, vp2, vp4, vs0, vs2

def layer_aniprop_model_sph_2(inmodel, dArr, nl, dh, ilvry):
    """
    Get the flattening transformed layerized model, deprecated
    ===================================================================
    ::: input parameters :::
    inmodel - input model1d object
    dArr    - numpy array of layer thickness (unit - km)
    nl      - number of layers 
    dh      - thickness of each layer (unit - km)
    ilvry   - 1 - Love wave; 2 - Rayleigh wave
    ===================================================================
    """
    # get the layerized model
    dArr, rhoArr, AArr, CArr, FArr, LArr, NArr = inmodel.get_layer_model(dArr, nl, dh)
    verby       = False
    nl_in       = rhoArr.size
    # Earth flattening transformation
    d_out,TA_out,TC_out,TF_out,TL_out,TN_out,TRho_out=tdisp96.flat2sphere(ilvry, 1., 1, 1, verby, \
                1, np.append(10., np.zeros(2048)), -1.,-1., dArr,AArr,CArr,FArr,LArr,NArr,rhoArr,nl_in, 0., 1, 0.5, 0.5)
    # convert layrized model to input model for aniprop
    zArr    = d_out.cumsum()   
    a       = (3.*TA_out + 3.*TC_out + 2.*TF_out + 4.*TL_out)/8.
    b       = 0.5*(TC_out - TA_out)
    c       = (TA_out + TC_out - 2.*TF_out - 4.*TL_out)/8.
    d       = 0.5*(TN_out + TL_out)
    e       = 0.5*(TL_out - TN_out)
    if np.any(e>0.0001):
        raise ValueError('aniprop does not accept L > N !')
    # total number of layers
    nl      = a.size
    # depth to each interface
    z       = np.zeros(2*nl+1, dtype=np.float32)
    z[2:]   = (np.repeat(zArr, 2))[:-1]
    z[0]    = 0.; z[1] = 0.
    # model parameters
    rho     = np.zeros(2*nl+1, dtype=np.float32)
    rho[1:] = np.repeat(TRho_out, 2)
    rho[0]  = rho[1]
    
    vp0     = np.zeros(2*nl+1, dtype=np.float32)
    vp0[1:] = np.repeat( np.sqrt(a/TRho_out), 2)
    vp0[0]  = vp0[1]
    
    vp2     = np.zeros(2*nl+1, dtype=np.float32)
    vp2[1:] = np.repeat(b/a, 2)
    vp2[0]  = vp2[1]
    
    vp4     = np.zeros(2*nl+1, dtype=np.float32)
    vp4[1:] = np.repeat(c/a, 2)
    vp4[0]  = vp4[1]
    
    vs0     = np.zeros(2*nl+1, dtype=np.float32)
    vs0[1:] = np.repeat(np.sqrt(d/TRho_out), 2)
    vs0[0]  = vs0[1]
    
    vs2     = np.zeros(2*nl+1, dtype=np.float32)
    vs2[1:] = np.repeat(e/d, 2)
    vs2[0]  = vs2[1]
    return z, rho, vp0, vp2, vp4, vs0, vs2

def plot(model, dtype='vsv', unit='km', showfig=True):
    """
    Plot the model
    ===================================================================
    ::: input parameters :::
    inmodel - input model1d object
    dtype   - datatype for plotting
    unit    - km or m
    ===================================================================
    """
    dtype   = dtype.lower()
    if unit == 'km': factor= 1000.
    elif unit == 'm': factor = 1.
    else: raise ValueError('Unexpected unit: '+unit)
    depth   = (6371000. - model.rArr[::-1])/factor
    if dtype == 'rho':
        data    = model.rhoArr.copy()
    elif dtype == 'vsv':
        data    = model.VsvArr.copy()
    elif dtype == 'vsh':
        data    = model.VshArr.copy()
    elif dtype == 'vpv':
        data    = model.VpvArr.copy()
    elif dtype == 'vph':
        data    = model.VphArr.copy()
    elif dtype == 'eta':
        data    = model.etaArr.copy()
    data    = data[::-1]
    if unit == 'km':
        dunitdict = {'rho': r'$\rho (g/cm^3)$', 'vsv': 'vsv (km/s)', 'vsh': 'vsh (km/s)', 'vpv': 'vpv (km/s)',\
                    'vph': 'vph (km/s)', 'eta': r'$\eta$' }
    else:
        dunitdict = {'rho': r'$\rho (kg/m^3)$', 'vsv': 'vsv (m/s)', 'vsh': 'vsh (m/s)', 'vpv': 'vpv (m/s)',\
                    'vph': 'vph (m/s)', 'eta': r'$\eta$' }
    data = data/factor
    ax=plt.subplot()
    plt.plot(data, depth, 'o-', ms=10, lw=3)
    plt.xlabel(dunitdict[dtype], fontsize=30)
    plt.ylabel('Depth ('+unit+')', fontsize=30)
    plt.gca().invert_yaxis()
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    if showfig: plt.show()
    return
    
    
@numba.jit(numba.float32[:, :](numba.float32[:], numba.float32, numba.boolean))
def _rot2mat(axis, angle, is_normalized):
    ''' Rotation matrix for rotation angle `angle` around `axis`
    ===============================================================================
    :::Important Note:::
    The rotation matrix is defined for rotation of a tensor in a fixed coordinate
    The output rotation matrix generated by this function is the inverse of the
    rotation matrix in Bond's book(p12-13).
    ===============================================================================
    Input Parameters:
    axis            - 3 element sequence, vector specifying axis for rotation.
    angle           - scalar, angle of rotation in degree.
    is_normalized   - bool, optional
       True if `axis` is already normalized (has norm of 1).  Default False
    -----
    output  -   mat : array shape (3,3), rotation matrix for specified rotation
    ===============================================================================
    Notes
    -----
    From: http://en.wikipedia.org/wiki/Rotation_matrix#Axis_and_angle
    '''
    x, y, z = axis
    if not is_normalized:
        n = np.sqrt(x*x + y*y + z*z)
        x = x/n
        y = y/n
        z = z/n
    angle  = np.pi*angle/180.
    c = np.cos(angle); s = np.sin(angle); C = 1-c
    xs = x*s;   ys = y*s;   zs = z*s
    xC = x*C;   yC = y*C;   zC = z*C
    xyC = x*yC; yzC = y*zC; zxC = z*xC
    g       = np.zeros((3,3), np.float32)
    g[0,:]  = np.array([ x*xC+c,   xyC-zs,   zxC+ys ])
    g[1,:]  = np.array([ xyC+zs,   y*yC+c,   yzC-xs ])
    g[2,:]  = np.array([ zxC-ys,   yzC+xs,   z*zC+c ])
    # g[:,:]  = np.array([
    #         [ x*xC+c,   xyC-zs,   zxC+ys ],
    #         [ xyC+zs,   y*yC+c,   yzC-xs ],
    #         [ zxC-ys,   yzC+xs,   z*zC+c ]], dtype=np.float32)
    return g

@numba.jit(numba.float32[:,:](numba.float32[:], numba.float32))
def _bondmat(axis, angle):
    """
    Compute Bond Matrix for rotation of Voigt matrix (eq. 8.9 in Bond, 1943; eq. 1.54 in Carcione, 2014)
    :::Important Note:::
    The rotation matrix used for Bond matrix was originally defined for rotation of the coordinate system
    ,which is in the opposite direction for rotation of a tensor in a fixed coordinate.
    We define the rotation as the rotation for the tensor itself in a fixed coordinate,
    Therefore,
    M   = bondmat(axis, angle) should be equal to the Bond Matrix for an angle of opposite sign
    =========================================================================================================
    Input Parameters:
    axis            - 3 element sequence, vector specifying axis for rotation.
    angle           - scalar, angle of rotation in degree.
    -----
    output          - array shape (6,6), Bond matrix for rotation of Voigt matrix
    =========================================================================================================
    """
    g       = _rot2mat(axis, angle, False)
    M       = np.zeros((6,6), np.float32)
    M[0,:]  = np.array([g[0,0]**2, g[0,1]**2, g[0,2]**2, 2.*g[0,1]*g[0,2], 2.*g[0,2]*g[0,0], 2.*g[0,0]*g[0,1]])
    M[1,:]  = np.array([g[1,0]**2, g[1,1]**2, g[1,2]**2, 2.*g[1,1]*g[1,2], 2.*g[1,2]*g[1,0], 2.*g[1,0]*g[1,1]])
    M[2,:]  = np.array([g[2,0]**2, g[2,1]**2, g[2,2]**2, 2.*g[2,1]*g[2,2], 2.*g[2,2]*g[2,0], 2.*g[2,0]*g[2,1]])
    M[3,:]  = np.array(\
        [g[1,0]*g[2,0], g[1,1]*g[2,1], g[1,2]*g[2,2], g[1,1]*g[2,2]+g[1,2]*g[2,1], g[1,0]*g[2,2]+g[1,2]*g[2,0], g[1,1]*g[2,0]+g[1,0]*g[2,1]])
    M[4,:]  = np.array(\
        [g[2,0]*g[0,0], g[2,1]*g[0,1], g[2,2]*g[0,2], g[0,1]*g[2,2]+g[0,2]*g[2,1], g[0,2]*g[2,0]+g[0,0]*g[2,2], g[0,0]*g[2,1]+g[0,1]*g[2,0]])
    M[5,:]  = np.array(\
        [g[0,0]*g[1,0], g[0,1]*g[1,1], g[0,2]*g[1,2], g[0,1]*g[1,2]+g[0,2]*g[1,1], g[0,2]*g[1,0]+g[0,0]*g[1,2], g[0,0]*g[1,1]+g[0,1]*g[1,0]])
    # M[:,:]  = np.array([[g[0,0]**2, g[0,1]**2, g[0,2]**2, 2.*g[0,1]*g[0,2], 2.*g[0,2]*g[0,0], 2.*g[0,0]*g[0,1]],
    #                     [g[1,0]**2, g[1,1]**2, g[1,2]**2, 2.*g[1,1]*g[1,2], 2.*g[1,2]*g[1,0], 2.*g[1,0]*g[1,1]],
    #                     [g[2,0]**2, g[2,1]**2, g[2,2]**2, 2.*g[2,1]*g[2,2], 2.*g[2,2]*g[2,0], 2.*g[2,0]*g[2,1]],
    #     [g[1,0]*g[2,0], g[1,1]*g[2,1], g[1,2]*g[2,2], g[1,1]*g[2,2]+g[1,2]*g[2,1], g[1,0]*g[2,2]+g[1,2]*g[2,0], g[1,1]*g[2,0]+g[1,0]*g[2,1]],
    #     [g[2,0]*g[0,0], g[2,1]*g[0,1], g[2,2]*g[0,2], g[0,1]*g[2,2]+g[0,2]*g[2,1], g[0,2]*g[2,0]+g[0,0]*g[2,2], g[0,0]*g[2,1]+g[0,1]*g[2,0]],
    #     [g[0,0]*g[1,0], g[0,1]*g[1,1], g[0,2]*g[1,2], g[0,1]*g[1,2]+g[0,2]*g[1,1], g[0,2]*g[1,0]+g[0,0]*g[1,2], g[0,0]*g[1,1]+g[0,1]*g[1,0]]
    #     ], dtype=np.float32)
    return M

####################################################
# Predefine the parameters for the model2d object
####################################################
spec = [
        # velocity parameters
        ('VsvArr', numba.float32[:]),
        ('VpvArr', numba.float32[:]),
        ('VshArr', numba.float32[:]),
        ('VphArr', numba.float32[:]),
        ('etaArr', numba.float32[:]),
        ('rhoArr', numba.float32[:]),
        # Earth flattening transformed density for Rayleigh wave
        ('rhoArrR', numba.float32[:]),
        # Earth flattening transformed density for Love wave
        ('rhoArrL', numba.float32[:]),
        # radius array
        ('rArr', numba.float32[:]),
        # depth array
        ('zArr', numba.float32[:]),
        # Earth flattening transformed radius
        ('rArrS', numba.float32[:]),
        # Love parameters 
        ('AArr', numba.float32[:]),
        ('CArr', numba.float32[:]),
        ('LArr', numba.float32[:]),
        ('FArr', numba.float32[:]),
        ('NArr', numba.float32[:]),
        # Earth flattening transformed Love parameters for Rayleigh wave
        ('AArrR', numba.float32[:]),
        ('CArrR', numba.float32[:]),
        ('LArrR', numba.float32[:]),
        ('FArrR', numba.float32[:]),
        ('NArrR', numba.float32[:]),
        # Earth flattening transformed Love parameters for Love wave
        ('AArrL', numba.float32[:]),
        ('CArrL', numba.float32[:]),
        ('LArrL', numba.float32[:]),
        ('FArrL', numba.float32[:]),
        ('NArrL', numba.float32[:]),
        # effective Love parameters
        ('AArrE', numba.float32[:]),
        ('CArrE', numba.float32[:]),
        ('LArrE', numba.float32[:]),
        ('FArrE', numba.float32[:]),
        ('NArrE', numba.float32[:]),
        # 2 theta azimuthal term
        ('BcArr', numba.float32[:]),
        ('BsArr', numba.float32[:]),
        ('GcArr', numba.float32[:]),
        ('GsArr', numba.float32[:]),
        ('HcArr', numba.float32[:]),
        ('HsArr', numba.float32[:]),
        # 4-theta azimuthal terms
        ('CcArr', numba.float32[:]),
        ('CsArr', numba.float32[:]),
        # Dip/strike angles of anisotropy, see Xie et al.(2015) Fig. 1 for details.
        ('dipArr', numba.float32[:]),
        ('strikeArr', numba.float32[:]),
        # Dip/strike angles of interfaces, used for raysum
        ('dipifArr', numba.float32[:]),
        ('strikeifArr', numba.float32[:]),
        ('dipping', numba.boolean),
        # 4th order elastic tensor for tilted hexagonal symmetric media
        ('CijklArr', numba.float32[:, :, :, :, :]),
        # Voigt matrix
        ('CijArr', numba.float32[:, :, :]),
        # Voigt matrix for AA(azimuthally independent) part
        ('CijAA', numba.float32[:, :, :]),
        # Voigt matrix for ETI(effective transversely isotropic) part
        ('CijETI', numba.float32[:, :, :]),
        # flat Earth or not
        ('flat', numba.boolean),
        ('rmin', numba.float32),
        ('tilt', numba.boolean)
        ]

@numba.jitclass(spec)
class model1d(object):
    """
    An object for handling a 1D Earth model
    =====================================================================================================================
    ::: parameters :::
    VsvArr, VshArr, - Vsv, Vsh, Vpv, Vph velocity (unit - m/s)
    VpvArr, VphArr  
    rhoArr          - density (kg/m^3)
    etaArr          - eta(F/(A-2L)) dimensionless
    AArr, CArr, FArr- Love parameters (unit - Pa)
    LArr, NArr
    rArr            - radius array (unit - m), sorted from the rmin to rmax(6371000. m)
    zArr            - depth array (unit - km), sorted as rArr
    flat            - = 0 spherical Earth, = 1 flat Earth (default)
                        Note: different from CPS
    arrays with *E  - Love parameters for effective VTI tensor
    arrays with *R  - Love parameters and density arrays after Earth flattening transformation for PSV motion
    arrays with *L  - Love parameters and density arrays after Earth flattening transformation for SH motion
    rArrS           - radius array after Earth flattening transformation
    dipArr,strikeArr- dip/strike angles, used for tilted hexagonal symmetric media
    CijArr          - elastic tensor given rotational angles(dip, strike) (unit - Pa)
    CijAA           - azimuthally anisotropic elastic tensor (unit - Pa)
    =====================================================================================================================
    """
    def __init__(self):
        self.flat   = True
        return
    
    def get_depth(self):
        self.zArr   = (np.float32(6371000.) - self.rArr)/np.float32(1000.)
        return
    
    def get_data_vel(self, vsv, vsh, vpv, vph, eta, rho, radius):
        """
        Get model data given velocity/density/radius arrays
        """
        self.rArr   = radius
        self.rhoArr = rho
        if radius[-1] != 6371000.:
            raise ValueError('Last element of radius array should be 6371000. meter !')
        if np.any(vsv<500.) or np.any(vsh<500.) or np.any(vpv<500.) or np.any(vph<500.) or np.any(rho<500.):
            raise ValueError('Wrong unit for model parameters!')
        if np.any(radius< 10000.):
            raise ValueError('Wrong unit for radius!')
        ###
        # assign velocities
        ###
        self.VsvArr = vsv
        self.VshArr = vsh
        self.VpvArr = vpv
        self.VphArr = vph
        self.etaArr = eta
        ###
        # compute Love parameters
        ###
        self.AArr   = rho * vph**2
        self.CArr   = rho * vpv**2
        self.LArr   = rho * vsv**2
        self.FArr   = eta * (self.AArr - np.float32(2.)* self.LArr)
        self.NArr   = rho * vsh**2
        return
    
    def vel2love(self):
        """
        velocity parameters to Love parameters
        """
        self.AArr   = self.rhoArr * (self.VphArr)**2
        self.CArr   = self.rhoArr * (self.VpvArr)**2
        self.LArr   = self.rhoArr * (self.VsvArr)**2
        self.FArr   = self.etaArr * (self.AArr - np.float32(2.)* self.LArr)
        self.NArr   = self.rhoArr * (self.VshArr)**2
        return
        
    def love2vel(self):
        """
        Love parameters to velocity parameters
        """
        self.VphArr = np.sqrt(self.AArr/self.rhoArr)
        self.VpvArr = np.sqrt(self.CArr/self.rhoArr)
        self.VshArr = np.sqrt(self.NArr/self.rhoArr)
        self.VsvArr = np.sqrt(self.LArr/self.rhoArr)
        self.etaArr = self.FArr/(self.AArr - np.float32(2.)* self.LArr)
        return
        
    def earth_flattening_kennett(self):
        """
        Earth flattening transformation using Kennett's formulas.
        reference
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
        self.flat   = False
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
        ::: input parameters :::
        zmax        - maximum depth (unit - km)
        dz          - depth interval (unit - km)
        Output:
        self.rArr   - radius array (unit - m)
        ===============================================================================
        """
        rmin        = 6371.0 - zmax
        self.rArr   = _get_array(rmin*1000., 6371000., dz*1000.)
        return
    
    def model_prem(self):
        """
        PREM model (Dziewonski & Anderson, PEPI 1981)
        The reference frequency is 1 Hz. Crust continued into the ocean.
        """
        ALst=[]; CLst=[]; LLst=[]; FLst=[]; NLst=[]; rhoLst=[]
        vsvLst=[]; vpvLst=[]; vshLst=[]; vphLst=[]; etaLst=[]
        for r in self.rArr:
            #- normalized radius
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
        """
        ak135 model from CPS
        """
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
        data        = ak135_cps_arr.reshape( ak135_cps_arr.size/4, 4)
        # assign model parameters
        z           = data[:, 0]
        radius      = (np.float32(6371.)-z)*np.float32(1000.)
        rho         = data[:, 3]*np.float32(1000.)
        vpv         = data[:, 1]*np.float32(1000.)
        vsv         = data[:, 2]*np.float32(1000.)
        vph         = vpv
        vsh         = vsv
        eta         = np.ones(vph.size, dtype=np.float32)
        # revert the array
        vsv         = vsv[::-1]
        vsh         = vsh[::-1]
        vpv         = vpv[::-1]
        vph         = vph[::-1]
        eta         = eta[::-1]
        rho         = rho[::-1]
        radius      = radius[::-1]
        # discard data with zero shear wave velocity
        ind         = (vsv!=0.)
        vsv         = vsv[ind]
        vsh         = vsh[ind]
        vpv         = vpv[ind]
        vph         = vph[ind]
        eta         = eta[ind]
        rho         = rho[ind]
        radius      = radius[ind]        
        self.get_data_vel(vsv, vsh, vpv, vph, eta, rho, radius)
        self.rmin   = self.rArr.min()
        return
    
    def is_iso(self):
        """Check if the model is isotropic at each point.
        """
        tol = 1e-5
        for i in xrange(self.rArr.size):
            if abs(self.AArr[i] - self.CArr[i])> tol or abs(self.LArr[i] - self.NArr[i])> tol or abs(self.FArr[i] - (self.AArr[i]- 2.*self.LArr[i]) )> tol:
                return False
        return True
        # if np.allclose(self.AArr, self.CArr) and np.allclose(self.LArr, self.NArr) and np.allclose(self.FArr, (self.AArr- 2.*self.LArr)):
        #     return True
        # else:
        #     return False
    
    def trim_simple(self, zmax):
        """
        Trim the model given a maximum depth(unit - km), keep data points with depth smaller than zmax.
        """
        if zmax >= 6371.: raise ValueError('Input maximum depth should have a unit of km !')
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
        if not self.flat:
            self.earth_flattening()
        return
    
    def trim(self, zmax):
        """
        Trim the model given a maximum depth(unit - km)
        """
        if zmax >= 6371.:
            raise ValueError('Input maximum depth should have a unit of km !')
        rmin        = (np.float32(6371.)-zmax)*np.float32(1000.)
        ind         = np.where(self.rArr == rmin)[0]
        if ind.size == 0:
            ind2    = np.where(self.rArr > 6171000.0)[0]
            ind0    = ind2[0]-1
            rl      = self.rArr[ind0];      rr      = self.rArr[ind0+1]
            rhol    = self.rhoArr[ind0];    rhor    = self.rhoArr[ind0+1]
            Al      = self.AArr[ind0];      Ar      = self.AArr[ind0+1]
            Cl      = self.CArr[ind0];      Cr      = self.CArr[ind0+1]
            Fl      = self.FArr[ind0];      Fr      = self.FArr[ind0+1]
            Ll      = self.LArr[ind0];      Lr      = self.LArr[ind0+1]
            Nl      = self.NArr[ind0];      Nr      = self.NArr[ind0+1]
            vsvl    = self.VsvArr[ind0];    vsvr    = self.VsvArr[ind0+1]
            vshl    = self.VshArr[ind0];    vshr    = self.VshArr[ind0+1]
            vpvl    = self.VpvArr[ind0];    vpvr    = self.VpvArr[ind0+1]
            vphl    = self.VphArr[ind0];    vphr    = self.VphArr[ind0+1]
            etal    = self.etaArr[ind0];    etar    = self.etaArr[ind0+1]
            changel = True
        elif ind.size > 2:
            raise ValueError('More than TWO repeated radius value, check the model !')
        else:
            if ind.size == 2:
                ind0 = ind[1]
            else:
                ind0 = ind[0]
            changel = False
        self.rArr   = self.rArr[ind0:]
        self.rhoArr = self.rhoArr[ind0:]
        # Love parameters
        self.AArr   = self.AArr[ind0:]
        self.CArr   = self.CArr[ind0:]
        self.FArr   = self.FArr[ind0:]
        self.LArr   = self.LArr[ind0:]
        self.NArr   = self.NArr[ind0:]
        # velocity
        self.VsvArr = self.VsvArr[ind0:]
        self.VshArr = self.VshArr[ind0:]
        self.VpvArr = self.VpvArr[ind0:]
        self.VphArr = self.VphArr[ind0:]
        self.etaArr = self.etaArr[ind0:]
        # change the first element value, linear interpolation is applied
        if changel:
            self.rArr[0]    = rmin
            self.rhoArr[0]  = rhol + (rmin - rl)*(rhor-rhol)/(rr - rl)
            # Love parameters
            self.AArr[0]    = Al + (rmin - rl)*(Ar-Al)/(rr - rl)
            self.CArr[0]    = Cl + (rmin - rl)*(Cr-Cl)/(rr - rl)
            self.FArr[0]    = Fl + (rmin - rl)*(Fr-Fl)/(rr - rl)
            self.LArr[0]    = Ll + (rmin - rl)*(Lr-Ll)/(rr - rl)
            self.NArr[0]    = Nl + (rmin - rl)*(Nr-Nl)/(rr - rl)
            # velocity
            self.VsvArr[0]  = vsvl + (rmin - rl)*(vsvr-vsvl)/(rr - rl)
            self.VshArr[0]  = vshl + (rmin - rl)*(vshr-vshl)/(rr - rl)
            self.VpvArr[0]  = vpvl + (rmin - rl)*(vpvr-vpvl)/(rr - rl)
            self.VphArr[0]  = vphl + (rmin - rl)*(vphr-vphl)/(rr - rl)
            self.etaArr[0]  = etal + (rmin - rl)*(etar-etal)/(rr - rl)
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
        if not self.flat:
            return self.rhoArrR[i], self.AArrR[i], self.CArrR[i], self.FArrR[i], self.LArrR[i], self.NArrR[i]
        else:
            return self.get_ind_Love_parameters(i)
    
    def get_ind_Love_parameters_SH(self, i):
        """
        Return Love paramaters and density given an index, for SH waves
        """
        if not self.flat:
            return self.rhoArrL[i], self.AArrL[i], self.CArrL[i], self.FArrL[i], self.LArrL[i], self.NArrL[i]
        else:
            return self.get_ind_Love_parameters(i)
    
    def get_r_love_parameters(self, r):
        """
        Return Love paramaters and density given a radius
        NOTE : always yield the RIGHT value if repeated radius grid points appear
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
            if self.rArr[ind_l] == self.rArr[ind_l+1]: ind_l+=1 # the difference compared with get_r_love_parameters_left
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
    
    def get_r_love_parameters_left(self, r):
        """
        Return Love paramaters and density given a radius
        NOTE : always yield the LEFT value if repeated radius grid points appear
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
        NOTE : always yield the RIGHT value if repeated radius grid points appear
        """
        if r < self.rmin: raise ValueError('Required radius is out of the model range!')
        if not self.flat:
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
                if self.rArrS[ind_l] == self.rArrS[ind_l+1]: ind_l+=1 #  yield the RIGHT values
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
    
    def get_r_love_parameters_PSV_left(self, r):
        """
        Return Love paramaters and density given a radius, for P-SV waves
        NOTE : always yield the RIGHT value if repeated radius grid points appear
        """
        if r < self.rmin: raise ValueError('Required radius is out of the model range!')
        if not self.flat:
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
        NOTE : always yield the RIGHT value if repeated radius grid points appear
        """
        if r < self.rmin: raise ValueError('Required radius is out of the model range!')
        if not self.flat:
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

    #####################################################################
    # functions for layerized model as input for aniprop
    # The unit for the input/output uses is the same as CPS,
    # i.e. km/s, GPa, km etc.
    #####################################################################
    def is_layer_model(self):
        """
        Check if the model is a layerized one or not
        """
        r_inv   = self.rArr[::-1]
        if r_inv.size %2 !=0:
            return False
        self.vel2love()
        rho_inv = self.rhoArr[::-1]
        A_inv   = self.AArr[::-1]
        C_inv   = self.CArr[::-1]
        F_inv   = self.FArr[::-1]
        L_inv   = self.LArr[::-1]
        N_inv   = self.NArr[::-1]
        dip_inv = self.dipArr[::-1]
        s_inv   = self.strikeArr[::-1]
        dif_inv = self.dipifArr[::-1]
        sif_inv = self.strikeifArr[::-1]
        if s_inv.size == r_inv.size:
            tilt    = True
        else:
            tilt    = False
        if sif_inv.size == r_inv.size:
            dipping = True
        else:
            dipping = False
        
        for i in xrange(r_inv.size):
            if i == 0:
                if r_inv[i] != 6371000.:
                    return False
                else:
                    continue
            if i % 2 != 0: continue
            r0  = r_inv[i-1]; r1 = r_inv[i]
            if r0 != r1:
                return False
            A0  = A_inv[i-2]; A1 = A_inv[i-1]
            if A0 != A1:
                return False
            C0  = C_inv[i-2]; C1 = C_inv[i-1]
            if C0 != C1:
                return False
            F0  = F_inv[i-2]; F1 = F_inv[i-1]
            if F0 != F1:
                return False
            L0  = L_inv[i-2]; L1 = L_inv[i-1]
            if L0 != L1:
                return False
            N0  = N_inv[i-2]; N1 = N_inv[i-1]
            if N0 != N1:
                return False
            # check tilted angles of anisotropic axis
            if tilt: 
                d0  = dip_inv[i-2]; d1 = dip_inv[i-1]
                if d0 != d1:
                    return False
                s0  = s_inv[i-2]; s1 = s_inv[i-1]
                if s0 != s1:
                    return False
            # check dipping interface angles
            if dipping:
                dif0= dif_inv[i-2]; dif1 = dif_inv[i-1]
                if dif0 != dif1:
                    return False
                sif0  = sif_inv[i-2]; sif1 = sif_inv[i-1]
                if sif0 != sif1:
                    return False
        return True
    
    def add_perturb_layer(self, zmin, zmax, dtype, val, rel):
        """
        Add/perturb a layer given zmin/zmax
        ===============================================================================
        ::: input parameters ::::
        zmin/zmax   - min/max depth (unit -km)
        dtype       - datatype for the perturbation
                        0: vsv
                        1: vsh
                        2: vpv
                        3: vph
                        4: eta
                        5: rho
        val         - perturbed/absolute value (for absolute value, unit is km)
        rel         - relative or absolute perturbation (True : relative)
        ===============================================================================
        NOTE:
        Benchmarked cases:
        A. one layer
        1. one layer, both overlap
            m.add_perturb_layer(0., 20., 0, 1., False)
            m.add_perturb_layer(0., 20., 0, -.5, True)
        2. one layer, both overlap, rmax not at the surface
            m.add_perturb_layer(20., 35., 0, 1., False)
            m.add_perturb_layer(20., 35., 0, -.5, True)
        3. one layer, rmax overlap
            m.add_perturb_layer(0., 10., 0, 1., False)
            m.add_perturb_layer(0., 10., 0, -.5, True)
        4. one layer, rmax overlap, rmax not at the surface
            m.add_perturb_layer(20., 25., 0, 1., False)
            m.add_perturb_layer(20., 25., 0, -.5, True)
        5. one layer, rmin overlap
            m.add_perturb_layer(25., 35., 0, 1., False)
            m.add_perturb_layer(25., 35., 0, -.5, True)
        6. one layer, rmin overlap, rmin at the bottom
            m.add_perturb_layer(197, 200., 0, 1., False)
            m.add_perturb_layer(197, 200., 0, -.5, True)
        7. both overlap, , rmin at the bottom
            m.add_perturb_layer(165, 200., 0, 1., False)
            m.add_perturb_layer(165, 200., 0, -.5, True)
        B. multiple layer
        1. multiple layer, both overlap
            m.add_perturb_layer(0., 35, 0, 1., False)
            m.add_perturb_layer(0., 35, 0, -.5, True)
        2. multiple layer, both overlap, rmax not at the surface
            m.add_perturb_layer(20., 77, 0, 1., False)
            m.add_perturb_layer(20., 77, 0, -.5, True)
        3. multiple layer, rmax overlap
            m.add_perturb_layer(0., 75, 0, 1., False)
            m.add_perturb_layer(0., 75, 0, -.5, True)
        4. multiple layer, rmax overlap, rmax not at the surface
            m.add_perturb_layer(20., 75, 0, 1., False)
            m.add_perturb_layer(20., 75, 0, -.5, True)
        5. multiple layer, rmin overlap
            m.add_perturb_layer(13., 77, 0, 1., False)
            m.add_perturb_layer(13., 77, 0, -.5, True)
        6. multiple layer, rmin overlap, rmin at the bottom
            m.add_perturb_layer(13., 200., 0, 1., False)
            m.add_perturb_layer(13., 200., 0, -.5, True)
        7. multiple layer, both overlap, , rmin at the bottom
            m.add_perturb_layer(20., 200., 0, 1., False)
            m.add_perturb_layer(20., 200., 0, -.5, True)
        """
        if rel:
            if val <= -1.: raise ValueError('Unacceptable value for added layer!')
        else:
            if val < 0.: raise ValueError('Unacceptable value for added layer!')
            if dtype !=4:
                val = val *np.float32(1000.)
        if not self.is_layer_model():
            raise ValueError('The model is not a layerized one!')
        if zmin >= zmax:
            raise ValueError('Minimum depth should be smaller than maximum depth!')
        rmax    = (np.float32(6371.)-zmin)*np.float32(1000.)
        rmin    = (np.float32(6371.)-zmax)*np.float32(1000.)
        if self.rArr[-1] < rmin or self.rArr[0] > rmax or self.rArr[-1] < rmax:
            raise ValueError('Input depth range out of bound!')
        rLst    = []; vsvLst=[]; vshLst=[]; vpvLst=[]; vphLst=[]; etaLst=[]; rhoLst=[]
        add0    = False; add1   = False
        for i in xrange(self.rArr.size):
            r   = self.rArr[i]
            vsv = self.VsvArr[i]
            vsh = self.VshArr[i]
            vpv = self.VpvArr[i]
            vph = self.VphArr[i]
            eta = self.etaArr[i]
            rho = self.rhoArr[i]
            if r < rmin or ((r >= rmax) and add1):
                rLst.append(r)
                vsvLst.append(vsv)
                vshLst.append(vsh)
                vpvLst.append(vpv)
                vphLst.append(vph)
                etaLst.append(eta)
                rhoLst.append(rho)
                continue
            if r == rmin:
                if i == 0 or self.rArr[i+1] > r:
                    if dtype == 0:
                        if rel: vsv = (1.+val)*vsv
                        else:   vsv = val
                    if dtype == 1:
                        if rel: vsh = (1.+val)*vsh
                        else:   vsh = val
                    if dtype == 2:
                        if rel: vpv = (1.+val)*vpv
                        else:   vpv = val
                    if dtype == 3:
                        if rel: vph = (1.+val)*vph
                        else:   vph = val
                    if dtype == 4:
                        if rel: eta = (1.+val)*eta
                        else:   eta = val
                    if dtype == 5:
                        if rel: rho = (1.+val)*rho
                        else:   rho = val
                    add0 = True
                rLst.append(r)
                vsvLst.append(vsv)
                vshLst.append(vsh)
                vpvLst.append(vpv)
                vphLst.append(vph)
                etaLst.append(eta)
                rhoLst.append(rho)
                continue
            if r == rmax and add0:
                if add1: raise ValueError('DEBUG!')      
                if dtype == 0:
                    if rel: vsv = (1.+val)*vsv
                    else:   vsv = val
                if dtype == 1:
                    if rel: vsh = (1.+val)*vsh
                    else:   vsh = val
                if dtype == 2:
                    if rel: vpv = (1.+val)*vpv
                    else:   vpv = val
                if dtype == 3:
                    if rel: vph = (1.+val)*vph
                    else:   vph = val
                if dtype == 4:
                    if rel: eta = (1.+val)*eta
                    else:   eta = val
                if dtype == 5:
                    if rel: rho = (1.+val)*rho
                    else:   rho = val
                rLst.append(r)
                vsvLst.append(vsv)
                vshLst.append(vsh)
                vpvLst.append(vpv)
                vphLst.append(vph)
                etaLst.append(eta)
                rhoLst.append(rho)
                add1 = True
                continue
            # radius value is larger than rmin, but smaller/equal to rmax
            if r <= rmax: 
                if not add0:
                    # if rmin grid points has not been added
                    # left value for rmin
                    rLst.append(rmin)
                    vsvLst.append(vsv)
                    vshLst.append(vsh)
                    vpvLst.append(vpv)
                    vphLst.append(vph)
                    etaLst.append(eta)
                    rhoLst.append(rho)
                    if dtype == 0:
                        if rel: vsv = (1.+val)*vsv
                        else:   vsv = val
                    if dtype == 1:
                        if rel: vsh = (1.+val)*vsh
                        else:   vsh = val
                    if dtype == 2:
                        if rel: vpv = (1.+val)*vpv
                        else:   vpv = val
                    if dtype == 3:
                        if rel: vph = (1.+val)*vph
                        else:   vph = val
                    if dtype == 4:
                        if rel: eta = (1.+val)*eta
                        else:   eta = val
                    if dtype == 5:
                        if rel: rho = (1.+val)*rho
                        else:   rho = val
                    # right value for rmin
                    rLst.append(rmin)
                    vsvLst.append(vsv)
                    vshLst.append(vsh)
                    vpvLst.append(vpv)
                    vphLst.append(vph)
                    etaLst.append(eta)
                    rhoLst.append(rho)
                    add0    = True
                    if add1: raise ValueError('DEBUG!')
                    rLst.append(r)
                    vsvLst.append(vsv)
                    vshLst.append(vsh)
                    vpvLst.append(vpv)
                    vphLst.append(vph)
                    etaLst.append(eta)
                    rhoLst.append(rho)
                    if r == rmax: add1 = True
                    continue
                else:
                    if dtype == 0:
                        if rel: vsv = (1.+val)*vsv
                        else:   vsv = val
                    if dtype == 1:
                        if rel: vsh = (1.+val)*vsh
                        else:   vsh = val
                    if dtype == 2:
                        if rel: vpv = (1.+val)*vpv
                        else:   vpv = val
                    if dtype == 3:
                        if rel: vph = (1.+val)*vph
                        else:   vph = val
                    if dtype == 4:
                        if rel: eta = (1.+val)*eta
                        else:   eta = val
                    if dtype == 5:
                        if rel: rho = (1.+val)*rho
                        else:   rho = val
                    rLst.append(r)
                    vsvLst.append(vsv)
                    vshLst.append(vsh)
                    vpvLst.append(vpv)
                    vphLst.append(vph)
                    etaLst.append(eta)
                    rhoLst.append(rho)
                    if r == rmax: add1 = True
                    continue
            else: 
                if add1: raise ValueError('DEBUG!')                
                if not add0:
                    # left value for rmin
                    rLst.append(rmin)
                    vsvLst.append(vsv)
                    vshLst.append(vsh)
                    vpvLst.append(vpv)
                    vphLst.append(vph)
                    etaLst.append(eta)
                    rhoLst.append(rho)
                    if dtype == 0:
                        if rel: vsv = (1.+val)*vsv
                        else:   vsv = val
                    if dtype == 1:
                        if rel: vsh = (1.+val)*vsh
                        else:   vsh = val
                    if dtype == 2:
                        if rel: vpv = (1.+val)*vpv
                        else:   vpv = val
                    if dtype == 3:
                        if rel: vph = (1.+val)*vph
                        else:   vph = val
                    if dtype == 4:
                        if rel: eta = (1.+val)*eta
                        else:   eta = val
                    if dtype == 5:
                        if rel: rho = (1.+val)*rho
                        else:   rho = val
                    # right value for rmin
                    rLst.append(rmin)
                    vsvLst.append(vsv)
                    vshLst.append(vsh)
                    vpvLst.append(vpv)
                    vphLst.append(vph)
                    etaLst.append(eta)
                    rhoLst.append(rho)
                    add0    = True
                # add rmax values to the list
                add1 = True
                # left value for rmax
                rLst.append(rmax)
                vsvLst.append(vsvLst[-1])
                vshLst.append(vshLst[-1])
                vpvLst.append(vpvLst[-1])
                vphLst.append(vphLst[-1])
                etaLst.append(etaLst[-1])
                rhoLst.append(rhoLst[-1])
                # right value for rmax
                rLst.append(rmax)
                vsvLst.append(self.VsvArr[i])
                vshLst.append(self.VshArr[i])
                vpvLst.append(self.VpvArr[i])
                vphLst.append(self.VphArr[i])
                etaLst.append(self.etaArr[i])
                rhoLst.append(self.rhoArr[i])
                # left value for r
                rLst.append(r)
                vsvLst.append(self.VsvArr[i])
                vshLst.append(self.VshArr[i])
                vpvLst.append(self.VpvArr[i])
                vphLst.append(self.VphArr[i])
                etaLst.append(self.etaArr[i])
                rhoLst.append(self.rhoArr[i])
                continue
        self.rArr   = np.array(rLst, dtype=np.float32)
        self.zArr   = (np.float32(6371000.) - self.rArr)/np.float32(1000.)
        self.rhoArr = np.array(rhoLst, dtype=np.float32)
        self.VsvArr = np.array(vsvLst, dtype=np.float32)
        self.VshArr = np.array(vshLst, dtype=np.float32)
        self.VpvArr = np.array(vpvLst, dtype=np.float32)
        self.VphArr = np.array(vphLst, dtype=np.float32)
        self.etaArr = np.array(etaLst, dtype=np.float32)
        self.vel2love()
        self.get_depth()
        if not self.is_layer_model():
            raise ValueError('DEBUG: The model is no longer a layerized one!')
        return
    
    def add_perturb_layer_love(self, zmin, zmax, dtype, val, rel):
        """
        Add/perturb a layer given zmin/zmax
        ===============================================================================
        ::: input parameters ::::
        zmin/zmax   - min/max depth (unit -km)
        dtype       - datatype for the perturbation
                        0: A
                        1: C
                        2: F
                        3: L
                        4: N
                        5: rho
        val         - perturbed/absolute value (for absolute value, unit is GPa)
        rel         - relative or absolute perturbation (True : relative)
        ===============================================================================
        """
        if rel:
            if val <= -1.: raise ValueError('Unacceptable value for added layer!')
        else:
            if val < 0.: raise ValueError('Unacceptable value for added layer!')
            if dtype == 5:
                val = val *np.float32(1000.)
            else:
                val = val *np.float32(1e9)
        if not self.is_layer_model():
            raise ValueError('The model is not a layerized one!')
        if zmin >= zmax:
            raise ValueError('Minimum depth should be smaller than maximum depth!')
        rmax    = (np.float32(6371.)-zmin)*np.float32(1000.)
        rmin    = (np.float32(6371.)-zmax)*np.float32(1000.)
        if self.rArr[-1] < rmin or self.rArr[0] > rmax or self.rArr[-1] < rmax:
            raise ValueError('Input depth range out of bound!')
        rLst    = []; ALst=[]; CLst=[]; FLst=[]; LLst=[]; NLst=[]; rhoLst=[]
        add0    = False; add1   = False
        for i in xrange(self.rArr.size):
            r   = self.rArr[i]
            A   = self.AArr[i]
            C   = self.CArr[i]
            F   = self.FArr[i]
            L   = self.LArr[i]
            N   = self.NArr[i]
            rho = self.rhoArr[i]
            if r < rmin or ((r >= rmax) and add1):
                rLst.append(r)
                ALst.append(A)
                CLst.append(C)
                FLst.append(F)
                LLst.append(L)
                NLst.append(N)
                rhoLst.append(rho)
                continue
            if r == rmin:
                if i == 0 or self.rArr[i+1] > r:
                    if dtype == 0:
                        if rel: A = (1.+val)*A
                        else:   A = val
                    if dtype == 1:
                        if rel: C = (1.+val)*C
                        else:   C = val
                    if dtype == 2:
                        if rel: F = (1.+val)*F
                        else:   F = val
                    if dtype == 3:
                        if rel: L = (1.+val)*L
                        else:   L = val
                    if dtype == 4:
                        if rel: N = (1.+val)*N
                        else:   N = val
                    if dtype == 5:
                        if rel: rho = (1.+val)*rho
                        else:   rho = val
                    add0 = True
                rLst.append(r)
                ALst.append(A)
                CLst.append(C)
                FLst.append(F)
                LLst.append(L)
                NLst.append(N)
                rhoLst.append(rho)
                continue
            if r == rmax and add0:
                if add1: raise ValueError('DEBUG!')      
                if dtype == 0:
                    if rel: A = (1.+val)*A
                    else:   A = val
                if dtype == 1:
                    if rel: C = (1.+val)*C
                    else:   C = val
                if dtype == 2:
                    if rel: F = (1.+val)*F
                    else:   F = val
                if dtype == 3:
                    if rel: L = (1.+val)*L
                    else:   L = val
                if dtype == 4:
                    if rel: N = (1.+val)*N
                    else:   N = val
                if dtype == 5:
                    if rel: rho = (1.+val)*rho
                    else:   rho = val
                rLst.append(r)
                ALst.append(A)
                CLst.append(C)
                FLst.append(F)
                LLst.append(L)
                NLst.append(N)
                rhoLst.append(rho)
                add1 = True
                continue
            # radius value is larger than rmin, but smaller/equal to rmax
            if r <= rmax: 
                if not add0:
                    # if rmin grid points has not been added
                    # left value for rmin
                    rLst.append(rmin)
                    ALst.append(A)
                    CLst.append(C)
                    FLst.append(F)
                    LLst.append(L)
                    NLst.append(N)
                    rhoLst.append(rho)
                    if dtype == 0:
                        if rel: A = (1.+val)*A
                        else:   A = val
                    if dtype == 1:
                        if rel: C = (1.+val)*C
                        else:   C = val
                    if dtype == 2:
                        if rel: F = (1.+val)*F
                        else:   F = val
                    if dtype == 3:
                        if rel: L = (1.+val)*L
                        else:   L = val
                    if dtype == 4:
                        if rel: N = (1.+val)*N
                        else:   N = val
                    if dtype == 5:
                        if rel: rho = (1.+val)*rho
                        else:   rho = val
                    # right value for rmin
                    rLst.append(rmin)
                    ALst.append(A)
                    CLst.append(C)
                    FLst.append(F)
                    LLst.append(L)
                    NLst.append(N)
                    rhoLst.append(rho)
                    add0    = True
                    if add1: raise ValueError('DEBUG!')
                    rLst.append(r)
                    ALst.append(A)
                    CLst.append(C)
                    FLst.append(F)
                    LLst.append(L)
                    NLst.append(N)
                    rhoLst.append(rho)
                    if r == rmax: add1 = True
                    continue
                else:
                    if dtype == 0:
                        if rel: A = (1.+val)*A
                        else:   A = val
                    if dtype == 1:
                        if rel: C = (1.+val)*C
                        else:   C = val
                    if dtype == 2:
                        if rel: F = (1.+val)*F
                        else:   F = val
                    if dtype == 3:
                        if rel: L = (1.+val)*L
                        else:   L = val
                    if dtype == 4:
                        if rel: N = (1.+val)*N
                        else:   N = val
                    if dtype == 5:
                        if rel: rho = (1.+val)*rho
                        else:   rho = val
                    rLst.append(r)
                    ALst.append(A)
                    CLst.append(C)
                    FLst.append(F)
                    LLst.append(L)
                    NLst.append(N)
                    rhoLst.append(rho)
                    if r == rmax: add1 = True
                    continue
            else: 
                if add1: raise ValueError('DEBUG!')                
                if not add0:
                    # left value for rmin
                    rLst.append(rmin)
                    ALst.append(A)
                    CLst.append(C)
                    FLst.append(F)
                    LLst.append(L)
                    NLst.append(N)
                    rhoLst.append(rho)
                    if dtype == 0:
                        if rel: A = (1.+val)*A
                        else:   A = val
                    if dtype == 1:
                        if rel: C = (1.+val)*C
                        else:   C = val
                    if dtype == 2:
                        if rel: F = (1.+val)*F
                        else:   F = val
                    if dtype == 3:
                        if rel: L = (1.+val)*L
                        else:   L = val
                    if dtype == 4:
                        if rel: N = (1.+val)*N
                        else:   N = val
                    if dtype == 5:
                        if rel: rho = (1.+val)*rho
                        else:   rho = val
                    # right value for rmin
                    rLst.append(rmin)
                    ALst.append(A)
                    CLst.append(C)
                    FLst.append(F)
                    LLst.append(L)
                    NLst.append(N)
                    rhoLst.append(rho)
                    add0    = True
                # add rmax values to the list
                add1 = True
                # left value for rmax
                rLst.append(rmax)
                ALst.append(ALst[-1])
                CLst.append(CLst[-1])
                FLst.append(FLst[-1])
                LLst.append(LLst[-1])
                NLst.append(NLst[-1])
                rhoLst.append(rhoLst[-1])
                # right value for rmax
                rLst.append(rmax)
                ALst.append(self.AArr[i])
                CLst.append(self.CArr[i])
                FLst.append(self.FArr[i])
                LLst.append(self.LArr[i])
                NLst.append(self.NArr[i])
                rhoLst.append(self.rhoArr[i])
                # left value for r
                rLst.append(r)
                ALst.append(self.AArr[i])
                CLst.append(self.CArr[i])
                FLst.append(self.FArr[i])
                LLst.append(self.LArr[i])
                NLst.append(self.NArr[i])
                rhoLst.append(self.rhoArr[i])
                continue
        self.rArr   = np.array(rLst, dtype=np.float32)
        self.rhoArr = np.array(rhoLst, dtype=np.float32)
        self.AArr   = np.array(ALst, dtype=np.float32)
        self.CArr   = np.array(CLst, dtype=np.float32)
        self.FArr   = np.array(FLst, dtype=np.float32)
        self.LArr   = np.array(LLst, dtype=np.float32)
        self.NArr   = np.array(NLst, dtype=np.float32)
        self.love2vel()
        self.get_depth()
        if not self.is_layer_model():
            raise ValueError('DEBUG: The model is no longer a layerized one!')
        return
        
    def get_default_layer_model(self):
        """
        Get default layerized model, the model need to be a layerized one
        """
        if not self.is_layer_model():
            raise ValueError('The model is not a layerized one!')
        r_inv   = self.rArr[::-1]
        dArr    = np.zeros(r_inv.size/2, dtype=np.float32)
        ind_odd = np.arange(r_inv.size/2)*2 +1
        ind_even= np.arange(r_inv.size/2)*2 
        dArr    = (r_inv[ind_even] - r_inv[ind_odd])/np.float32(1000.)
        return self.get_layer_model(dArr, 1, 1.)
    
    def get_layer_model(self, dArr, nl, dh):
        """
        Get the layerized model for CPS
        Note: the unit is different from the default unit of the object
        ===================================================================
        ::: input parameters ::: 
        dArr            - numpy array of layer thickness (unit - km)
        nl              - number of layers 
        dh              - thickness of each layer (unit - km)
        nl and dh will be used if and only if dArr.size = 0
        ::: output :::
        dArr            - layer thickness array (unit - km)
        rhoArr          - density array (unit - g/cm^3)
        AArr, CArr, FArr- Love parameters (unit - GPa)
        LArr, NArr
        ===================================================================
        """
        if dArr.size==0:
            dh      *= 1000. 
            dArr    = np.ones(nl, dtype = np.float32)*np.float32(dh)
        else:
            dArr    *= 1000.
            nl      = dArr.size
        ALst    = []; CLst  = []; LLst  = []; FLst  = []; NLst  = []; rhoLst= []
        z0   = 0.
        z1   = dArr[0]
        for i in xrange(nl):
            r0  = 6371000.-z0
            r1  = 6371000.-z1
            rho0, A0, C0, F0, L0, N0  = self.get_r_love_parameters_left(r0) # bottom point value needs to use left value in radius array
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
            z1  += dArr[i+1]
        # layer thickness is converted from m to km
        dArr    /=1000.
        rhoArr  = np.array(rhoLst, dtype=np.float32)
        AArr    = np.array(ALst, dtype=np.float32)
        CArr    = np.array(CLst, dtype=np.float32)
        FArr    = np.array(FLst, dtype=np.float32)
        LArr    = np.array(LLst, dtype=np.float32)
        NArr    = np.array(NLst, dtype=np.float32)
        return dArr, rhoArr, AArr, CArr, FArr, LArr, NArr
    
    def get_layer_tilt_model(self, dArr, nl, dh):
        """
        Get the layerized tilted model for CPS
        Note: the unit is different from the default unit of the object
        ===================================================================
        ::: input parameters :::
        dArr                - numpy array of layer thickness (unit - km)
        nl                  - number of layers 
        dh                  - thickness of each layer (unit - km)
        nl and dh will be used if and only if dArr.size = 0
        ::: output :::
        dArr                - layer thickness array (unit - km)
        rhoArr              - density array (unit - g/cm^3)
        AArr, CArr, FArr    - effective Love parameters (unit - GPa)
        LArr, NArr
        BcArr, BsArr, GcArr - 2-theta azimuthal terms (unit - GPa)
        GsArr, HcArr, HsArr
        CcArr, CsArr        - 4-theta azimuthal terms (unit - GPa)
        ===================================================================
        """
        if not self.is_layer_model():
            print 'WARNING: Model is not layerized, unexpected error may occur!'
        if dArr.size==0:
            dh      *= 1000. 
            dArr    = np.ones(nl, dtype = np.float32)*np.float32(dh)
        else:
            dArr    *= 1000.
            nl      = dArr.size
        ALst    = []; CLst  = []; LLst  = []; FLst  = []; NLst  = []; rhoLst= []
        BcLst   = []; BsLst = []; GcLst = []; GsLst = []; HcLst = []; HsLst = []; CcLst = []; CsLst = []
        z0   = 0.
        z1   = dArr[0]
        for i in xrange(nl):
            r0  = 6371000.-z0
            r1  = 6371000.-z1
            rho0, A0, C0, F0, L0, N0, Bc0, Bs0, Gc0, Gs0, Hc0, Hs0, Cc0, Cs0\
                        =self.get_r_tilt_parameters(r0, True) # bottom point value needs to use left value in radius array
            rho1, A1, C1, F1, L1, N1, Bc1, Bs1, Gc1, Gs1, Hc1, Hs1, Cc1, Cs1\
                        =self.get_r_tilt_parameters(r1, False)
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
            # 2-theta azimuthal terms
            Bc  = (Bc0+Bc1)/1.e9/2.
            BcLst.append(Bc)
            Bs  = (Bs0+Bs1)/1.e9/2.
            BsLst.append(Bs)
            Gc  = (Gc0+Gc1)/1.e9/2.
            GcLst.append(Gc)
            Gs  = (Gs0+Gs1)/1.e9/2.
            GsLst.append(Gs)
            Hc  = (Hc0+Hc1)/1.e9/2.
            HcLst.append(Hc)
            Hs  = (Hs0+Hs1)/1.e9/2.
            HsLst.append(Hs)
            # 4-theta azimuthal terms
            Cc  = (Cc0+Cc1)/1.e9/2.
            CcLst.append(Cc)
            Cs  = (Cs0+Cs1)/1.e9/2.
            CsLst.append(Cs)
            z0  += dArr[i]
            z1  += dArr[i+1]
        # layer thickness is converted from m to km
        dArr    /=1000.
        rhoArr  = np.array(rhoLst, dtype=np.float32)
        AArr    = np.array(ALst, dtype=np.float32)
        CArr    = np.array(CLst, dtype=np.float32)
        FArr    = np.array(FLst, dtype=np.float32)
        LArr    = np.array(LLst, dtype=np.float32)
        NArr    = np.array(NLst, dtype=np.float32)
        BcArr   = np.array(BcLst, dtype=np.float32)
        BsArr   = np.array(BsLst, dtype=np.float32)
        GcArr   = np.array(GcLst, dtype=np.float32)
        GsArr   = np.array(GsLst, dtype=np.float32)
        HcArr   = np.array(HcLst, dtype=np.float32)
        HsArr   = np.array(HsLst, dtype=np.float32)
        CcArr   = np.array(CcLst, dtype=np.float32)
        CsArr   = np.array(CsLst, dtype=np.float32)
        return dArr, rhoArr, AArr, CArr, FArr, LArr, NArr, BcArr, BsArr, GcArr, GsArr, HcArr, HsArr, CcArr, CsArr
    
    #####################################################################
    # functions for layerized model as input for aniprop
    #####################################################################
    def layer_aniprop_model(self, dArr, nl, dh):
        """
        Get layrized model for aniprop
        ===================================================================================
        ::: input parameters :::
        dArr            - numpy array of layer thickness (unit - km)
        nl              - number of layers 
        dh              - thickness of each layer (unit - km)
        nl and dh will be used if and only if dArr.size = 0
        ::: output :::
        z               - depth array to interfaces (unit - km)
        rho             - density array (unit - g/cm^3)
        vp0, vp2, vp4   - P wave velocity and corresponding 2psi/4psi relative perturnation
        vs0, vs2        - S wave velocity and corresponding 2psi relative perturnation
        ===================================================================================
        """
        dArr, rhoArr, AArr, CArr, FArr, LArr, NArr = self.get_layer_model(dArr, nl, dh)
        zArr    = dArr.cumsum()
        # Get the a,b,c,d,e arrays
        a       = (3.*AArr + 3.*CArr + 2.*FArr + 4.*LArr)/8.
        b       = 0.5*(CArr - AArr)
        c       = (AArr + CArr - 2.*FArr - 4.*LArr)/8.
        d       = 0.5*(NArr + LArr)
        e       = 0.5*(LArr - NArr)
        # total number of laryers
        nl      = a.size
        z       = zArr
        rho     = rhoArr
        vp0     = np.sqrt(a/rhoArr)
        vp2     = b/a
        vp4     = c/a
        vs0     = np.sqrt(d/rhoArr)
        vs2     = e/d
        return z, rho, vp0, vp2, vp4, vs0, vs2
    
    def layer_aniprop_model_0(self, dArr, nl, dh):
        """
        Get layrized model for aniprop, add 0 to top, deprecated
        ===================================================================================
        ::: input parameters :::
        dArr            - numpy array of layer thickness (unit - km)
        nl              - number of layers 
        dh              - thickness of each layer (unit - km)
        nl and dh will be used if and only if dArr.size = 0
        ::: output :::
        z               - depth array to interfaces (unit - km)
        rho             - density array (unit - g/cm^3)
        vp0, vp2, vp4   - P wave velocity and corresponding 2psi/4psi relative perturnation
        vs0, vs2        - S wave velocity and corresponding 2psi relative perturnation
        ===================================================================================
        """
        dArr, rhoArr, AArr, CArr, FArr, LArr, NArr = self.get_layer_model(dArr, nl, dh)
        zArr    = dArr.cumsum()
        # Get the a,b,c,d,e arrays
        a       = (3.*AArr + 3.*CArr + 2.*FArr + 4.*LArr)/8.
        b       = 0.5*(CArr - AArr)
        c       = (AArr + CArr - 2.*FArr - 4.*LArr)/8.
        d       = 0.5*(NArr + LArr)
        e       = 0.5*(LArr - NArr)
        # total number of laryers
        nl      = a.size
        z       = np.zeros(nl+1, dtype=np.float32)
        z[1:]   = zArr
        z[0]    = 0.
        
        rho     = np.zeros(nl+1, dtype=np.float32)
        rho[1:] = rhoArr
        rho[0]  = rho[1]
        
        vp0     = np.zeros(nl+1, dtype=np.float32)
        vp0[1:] = np.sqrt(a/rhoArr)
        vp0[0]  = vp0[1]
        
        vp2     = np.zeros(nl+1, dtype=np.float32)
        vp2[1:] = b/a
        vp2[0]  = vp2[1]
        
        vp4     = np.zeros(nl+1, dtype=np.float32)
        vp4[1:] = c/a
        vp4[0]  = vp4[1]
        
        vs0     = np.zeros(nl+1, dtype=np.float32)
        vs0[1:] = np.sqrt(d/rhoArr)
        vs0[0]  = vs0[1]
        
        vs2     = np.zeros(nl+1, dtype=np.float32)
        vs2[1:] = e/d
        vs2[0]  = vs2[1]
        return z, rho, vp0, vp2, vp4, vs0, vs2
    
    def angles_aniprop_model(self, zArr):
        """
        Get dip/strike arrays given depth arrays
        ===================================================================================
        ::: input parameters :::
        zArr            - depth array (unit - km)
        ::: output :::
        dipArr/strikeArr- dip/strike arrays (unit - degree) 
        ===================================================================================
        """
        if not self.is_layer_model():
            print ('WARNING: Model is not layerized, unexpected error may occur!')
        dipLst  = []; strikeLst = []
        nl      = zArr.size
        for i in xrange(nl):
            r           = 6371000.-zArr[i]*1000.
            dip, strike = self.get_dip_strike(r, False)
            # Love parameters are converted from Pa to GPa
            dipLst.append(dip)
            strikeLst.append(strike)
        # layer thickness is converted from m to km
        dipArr      = np.array(dipLst, dtype=np.float32)
        strikeArr   = np.array(strikeLst, dtype=np.float32)
        return dipArr, strikeArr
    
    def aniprop_check_model(self):
        """check the model
        """
        if np.any(self.LArr > self.NArr):
            raise ValueError('aniprop does not accept L > N !')
        return
    
    ###############################################################################################
    # functions for raysum model
    ###############################################################################################
    def angles_raysum_model(self, dArr, atype):
        """
        Get dip/strike arrays given depth arrays
        =============================================================================================================
        ::: input parameters :::
        zArr            - depth array (unit - km)
        atype           - angle type (0 - anisotropic axis angles; others - dipping interface angles)
        ::: output :::
        dipArr/strikeArr- dip/strike arrays (unit - degree) 
        =============================================================================================================
        """
        zArr    = dArr.cumsum()
        if atype == 0:
            return self.angles_aniprop_model(zArr)
        if not self.is_layer_model():
            print ('WARNING: Model is not layerized, unexpected error may occur!')
        dipLst  = []; strikeLst = []
        nl      = zArr.size
        for i in xrange(nl):
            r           = 6371000.-zArr[i]*1000.
            dip, strike = self.get_dip_strike_interface(r, False)
            # Love parameters are converted from Pa to GPa
            dipLst.append(dip)
            strikeLst.append(strike)
        # layer thickness is converted from m to km
        dipArr      = np.array(dipLst, dtype=np.float32)
        strikeArr   = np.array(strikeLst, dtype=np.float32)
        return dipArr, strikeArr
    
    def layer_raysum_model(self, dArr, nl, dh):
        """
        Get layrized model for raysum
        ===================================================================================
        ::: input parameters :::
        dArr            - numpy array of layer thickness (unit - km)
        nl              - number of layers 
        dh              - thickness of each layer (unit - km)
        nl and dh will be used if and only if dArr.size = 0
        ::: output :::
        dArr            - thickness array (unit - km)
        rho             - density array (unit - g/cm^3)
        vp0, vp2, vp4   - P wave velocity and corresponding 2psi/4psi relative perturnation
        vs0, vs2        - S wave velocity and corresponding 2psi relative perturnation
        ===================================================================================
        """
        dArr, rhoArr, AArr, CArr, FArr, LArr, NArr = self.get_layer_model(dArr, nl, dh)
        vp      = (np.sqrt(CArr/rhoArr) + np.sqrt(AArr/rhoArr))/2.
        dvp     = (np.sqrt(CArr/rhoArr) - np.sqrt(AArr/rhoArr))/vp
        
        vs      = (np.sqrt(LArr/rhoArr) + np.sqrt(NArr/rhoArr))/2.
        dvs     = (np.sqrt(LArr/rhoArr) - np.sqrt(NArr/rhoArr))/vs
        iso     = np.zeros(dArr.size, dtype=np.int32)
        iso[(dvp==0.)*(dvs==0.)] = 1
        return dArr, rhoArr, vp, vs, dvp, dvs, iso
    
    ####################################################################################################
    # functions for tilted hexaganal symmetric model
    ####################################################################################################
    def init_dip_strike(self):
        """initialize dip/strike angle array
        """
        self.tilt       = True
        self.dipArr     = np.zeros(self.rArr.size, np.float32)
        self.strikeArr  = np.zeros(self.rArr.size, np.float32)
        return
    
    def is_stable(self):
        """
        Check that the elastic constants matrix is positive definite 
        That is, check that the structure is stable to small strains. This
        is done by finding the eigenvalues of the Voigt elastic stiffness matrix
        by diagonalization and checking that they are all positive.
        See Born & Huang, "Dynamical Theory of Crystal Lattices" (1954) page 141.
        """
        for i in xrange(self.rArr.size):
            A           = self.AArr[i]
            C           = self.CArr[i]
            F           = self.FArr[i]
            L           = self.LArr[i]
            N           = self.NArr[i]
            Cvoigt      = np.zeros((6,6), dtype=np.float32)
            Cvoigt[0,0] = A; Cvoigt[1,1] = A; Cvoigt[2,2]= C; Cvoigt[3,3]= L; Cvoigt[4,4] = L; Cvoigt[5,5] = N
            Cvoigt[0,1] = A-2.*N;   Cvoigt[1,0] = A-2.*N
            Cvoigt[0,2] = F;        Cvoigt[2,0] = F
            Cvoigt[1,2] = F;        Cvoigt[2,1] = F
            (eigenvalues, eigenvectors) = np.linalg.eig(Cvoigt)
            if not (eigenvalues.min() > 0.0):
                print self.rArr[i]/1000.,'km NOT stable!'
                return False
                # raise ValueError('Elastic tensor is not stable to small strains (Voigt matrix is not positive definite) at r='+rstr+' m')
        # if verbose: print 'Stability checked! Eigenvalues:', eigenvalues
        return True
    
    def init_etensor(self):
        """initialize elastic tensor arrays
        """
        # self.CijklArr   = np.zeros((3,3,3,3,self.rArr.size), dtype = np.float32)
        self.CijArr     = np.zeros((6,6,self.rArr.size), dtype = np.float32)
        self.CijAA      = np.zeros((6,6,self.rArr.size), dtype = np.float32)
        for i in xrange(self.rArr.size):
            A           = self.AArr[i]
            C           = self.CArr[i]
            F           = self.FArr[i]
            L           = self.LArr[i]
            N           = self.NArr[i]
            Cvoigt      = np.zeros((6,6), dtype=np.float32)
            Cvoigt[0,0] = A; Cvoigt[1,1] =A; Cvoigt[2,2]=C; Cvoigt[3,3]=L; Cvoigt[4,4]=L; Cvoigt[5,5]=N
            Cvoigt[0,1] = A-2.*N; Cvoigt[1,0] = A-2.*N
            Cvoigt[0,2] = F; Cvoigt[2,0] = F
            Cvoigt[1,2] = F; Cvoigt[2,1] = F
            self.CijArr[:, :, i] = Cvoigt
        # self.CijETI     = self.CijArr
        # self.voigt2Cijkl()
        return
    
    def init_tilt(self):
        """initialize tilted hexagonal symmetric model
        """
        self.init_dip_strike(); self.init_etensor()
        return
    
    # def voigt2Cijkl(self):
    #     """
    #     Convert Voigt matrix to 4th order elastic tensor 
    #     """
    #     m2t = np.zeros((3,3), dtype=np.int32)
    #     m2t[0,1] = 5; m2t[1,0] = 5
    #     m2t[0,2] = 4; m2t[2,0] = 4
    #     m2t[2,1] = 3; m2t[1,2] = 3
    #     m2t[1,1] = 1; m2t[2,2] = 2
    #     for ir in xrange(self.rArr.size):
    #         Cvoigt = self.CijArr[:,:,ir]
    #         for i in xrange(3):
    #             for j in xrange(3):
    #                 for k in xrange(3):
    #                     for l in xrange(3):
    #                         self.CijklArr[i, j, k, l, ir] = Cvoigt[m2t[i,j], m2t[k,l]]
    #     return
    # 
    # def Cijkl2voigt(self):
    #     """
    #     Convert 4th order elastic tensor to Voigt matrix
    #     """
    #     t2m = np.zeros((2,6), dtype=np.int32)
    #     t2m[0,1] = 1; t2m[0,2] = 2; t2m[0,3] = 1; t2m[0,4] = 2;
    #     t2m[1,1] = 1; t2m[1,2] = 2; t2m[1,3] = 2; t2m[1,5] = 1;
    #     for ir in xrange(self.rArr.size):
    #         Cijkl   = self.CijklArr[:,:,:,:, ir]
    #         for i in xrange(6):
    #             for j in xrange(6):
    #                 self.CijArr[i,j,ir] = Cijkl[t2m[0,i],t2m[1,i],t2m[0,j],t2m[1,j]]
    #     return
    
    def is_tilt_model(self, check):
        """check whether the model is a tilted hexaganal symmetric model or not
        """
        if (self.dipArr.size != self.rArr.size) or (self.strikeArr.size != self.rArr.size):
            if check:
                raise ValueError('Dip/Strike angles not specified!')
            else:
                self.init_dip_strike()
                return False
        return True
    
    def rot_dip_strike(self):
        """
        Rotate Voigt matrix using Bond matrix (eq. 1.58 in Carcione, 2014)
        Note that the rotation is the inverse of rotation of a coordinate system,
        thus the rotation matrix used to construct Bond matrix is the inverse of the
        rotation matrix in Bond's book (p12-13)
        """
        for i in xrange(self.rArr.size):
            if self.dipArr[i] == 0.: continue
            if self.dipArr[i] >90. or self.dipArr[i] < 0.: raise ValueError('Dip should be within [0., 90.]!')
            # rotation for dip, rotation axis is x (North)
            Mdip        = _bondmat(np.array([1.,0.,0.], dtype=np.float32), self.dipArr[i])
            Cvoigt      = np.dot(Mdip, self.CijArr[:,:, i])
            Cvoigt      = np.dot(Cvoigt, Mdip.T)
            # rotation for strike, rotation axis is z (downward)
            Mstrike     = _bondmat(np.array([0.,0.,1.], dtype=np.float32), self.strikeArr[i])
            Cvoigt      = np.dot(Mstrike, Cvoigt)
            Cvoigt      = np.dot(Cvoigt, Mstrike.T)
            self.CijArr[:,:, i] = Cvoigt
            # print(i)
        return
    
    def decompose(self):
        """
        Decompose the tilted elastic tensor into ETI and AA components
        """
        # initialize effective Love parameters
        self.AArrE  = self.AArr.copy()
        self.CArrE  = self.CArr.copy()
        self.FArrE  = self.FArr.copy()
        self.LArrE  = self.LArr.copy()
        self.NArrE  = self.NArr.copy()
        # initialize 2-theta azimuthal terms
        self.BcArr  = np.zeros(self.rArr.size, np.float32)
        self.BsArr  = np.zeros(self.rArr.size, np.float32)
        self.GcArr  = np.zeros(self.rArr.size, np.float32)
        self.GsArr  = np.zeros(self.rArr.size, np.float32)
        self.HcArr  = np.zeros(self.rArr.size, np.float32)
        self.HsArr  = np.zeros(self.rArr.size, np.float32)
        self.CcArr  = np.zeros(self.rArr.size, np.float32)
        self.CsArr  = np.zeros(self.rArr.size, np.float32)
        for i in xrange(self.rArr.size):
            if self.dipArr[i] == 0.: continue
            Cij             = self.CijArr[:,:, i]
            A               = 3./8.*(Cij[0,0] + Cij[1,1]) + Cij[0,1]/4. + Cij[5,5]/2.
            C               = Cij[2,2]
            F               = (Cij[0,2] + Cij[1,2])/2.
            L               = (Cij[3,3] + Cij[4,4])/2.
            N               = (Cij[0,0] + Cij[1,1])/8. - Cij[0,1]/4. + Cij[5,5]/2.
            CijETI          = np.zeros((6,6), np.float32)
            CijETI[0,:]     = np.array([A, A-2.*N, F, 0., 0., 0.])
            CijETI[1,:]     = np.array([A-2.*N, A, F, 0., 0., 0.])
            CijETI[2,:]     = np.array([F, F, C, 0., 0., 0.])
            CijETI[3,:]     = np.array([0., 0., 0., L, 0., 0.])
            CijETI[4,:]     = np.array([0., 0., 0., 0., L, 0.])
            CijETI[5,:]     = np.array([0., 0., 0., 0., 0., N])
            self.CijAA[:,:,i]   = self.CijArr[:,:, i] - CijETI
            CijAA               = self.CijAA[:,:,i]
            # # # CijAA               = self.CijArr[:,:, i]
            self.AArrE[i]       = A
            self.CArrE[i]       = C
            self.FArrE[i]       = F
            self.LArrE[i]       = L
            self.NArrE[i]       = N
            self.BcArr[i]       = (CijAA[0,0] - CijAA[1,1])/2.
            self.BsArr[i]       = CijAA[0,5] + CijAA[1,5]
            self.GcArr[i]       = (CijAA[4,4] - CijAA[3,3])/2.
            self.GsArr[i]       = CijAA[4,3]
            self.HcArr[i]       = (CijAA[0,2] - CijAA[1,2])/2.
            self.HsArr[i]       = CijAA[2,5]
            self.CcArr[i]       = (CijAA[0,0] + CijAA[1,1])/8. - CijAA[0,1]/4. - CijAA[5,5]/2.
            self.CsArr[i]       = (CijAA[0,5] - CijAA[1,5])/2.
        return

    
    # # # def get_hex_parameter(self, z, left):
    # # #     self.is_tilt_model(True)
    # # #     if not self.is_layer_model():
    # # #         print 'WARNING: the model is not layrized!'
    # # #     r   = np.float32( (6371.- z)*1000.)
    # # #     ind = np.where(self.rArr == r)[0]
    # # #     if ind.size == 0:
    # # #         if left:
    # # #             rho, A, C, F, L, N  = self.get_r_love_parameters_left(r)
    # # #             indds   = (np.where(self.rArr<r)[0])[-1]
    # # #             dip     = self.dipArr[indds]; strike     = self.strikeArr[indds]
    # # #         else:
    # # #             rho, A, C, F, L, N  = self.get_r_love_parameters(r)
    # # #             indds   = (np.where(self.rArr>r)[0])[0]
    # # #             dip     = self.dipArr[indds]; strike     = self.strikeArr[indds]
    # # #     else:
    # # #         if left:
    # # #             indds   = ind[0]
    # # #         else:
    # # #             indds   = ind[-1]
    # # #         rho     = self.rhoArr[indds]
    # # #         A       = self.AArr[indds]
    # # #         C       = self.CArr[indds]
    # # #         F       = self.FArr[indds]
    # # #         L       = self.LArr[indds]
    # # #         N       = self.NArr[indds]
    # # #         dip     = self.dipArr[indds]; strike     = self.strikeArr[indds]
    # # #     return rho, A, C, F, L, N, dip, strike
    
    def get_dip_strike(self, r, left):
        """
        Return dip/strike of anisotrpic axis given a radius
        ===================================================================================
        ::: input parameters :::
        r           - radius (unit - m)
        left        - yield the LEFT value if repeated radius grid points appear or not
        ::: output :::
        dip/strike  - dip/strike (unit - degree) 
        ===================================================================================
        """
        self.is_tilt_model(True)
        if not self.is_layer_model():
            print 'WARNING: the model is not layrized!'
        # r   = np.float32( (6371.- z)*1000.)
        ind = np.where(self.rArr == r)[0]
        if ind.size == 0:
            if left:
                indds   = (np.where(self.rArr<r)[0])[-1]
                dip     = self.dipArr[indds]; strike     = self.strikeArr[indds]
            else:
                indds   = (np.where(self.rArr>r)[0])[0]
                dip     = self.dipArr[indds]; strike     = self.strikeArr[indds]
        else:
            if left:
                indds   = ind[0]
            else:
                indds   = ind[-1]
            dip     = self.dipArr[indds]; strike     = self.strikeArr[indds]
        return dip, strike
    
    def get_r_tilt_parameters(self, r, left):
        """
        Return density, effective Love paramaters and 2-theta azimuthal terms given a radius
        left    - yield the LEFT value if repeated radius grid points appear or not
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
            if self.rArr[ind_l] == self.rArr[ind_l+1]:
                if not left:
                    ind_l+=1 # the difference compared with get_r_tilt_parameters_left
            return self.rhoArr[ind_l],\
                self.AArrE[ind_l], self.CArrE[ind_l], self.FArrE[ind_l], self.LArrE[ind_l], self.NArrE[ind_l],\
                self.BcArr[ind_l], self.BsArr[ind_l], self.GcArr[ind_l], self.GsArr[ind_l], self.HcArr[ind_l],\
                self.HsArr[ind_l], self.CcArr[ind_l], self.CsArr[ind_l]
        r_left  = self.rArr[ind_l]
        r_right = self.rArr[ind_r]
        
        rhol    = self.rhoArr[ind_l]; rhor  =   self.rhoArr[ind_r]
        rho     = rhol + (r - r_left)*(rhor-rhol)/(r_right - r_left)
        # effective Love parameters
        Al      = self.AArrE[ind_l]; Ar  =   self.AArrE[ind_r]
        A       = Al + (r - r_left)*(Ar-Al)/(r_right - r_left)
        
        Cl      = self.CArrE[ind_l]; Cr  =   self.CArrE[ind_r]
        C       = Cl + (r - r_left)*(Cr-Cl)/(r_right - r_left)    
        
        Fl      = self.FArrE[ind_l]; Fr  =   self.FArrE[ind_r]
        F       = Fl + (r - r_left)*(Fr-Fl)/(r_right - r_left)
        
        Ll      = self.LArrE[ind_l]; Lr  =   self.LArrE[ind_r]
        L       = Ll + (r - r_left)*(Lr-Ll)/(r_right - r_left)
        
        Nl      = self.NArrE[ind_l]; Nr  =   self.NArrE[ind_r]
        N       = Nl + (r - r_left)*(Nr-Nl)/(r_right - r_left)
        # 2-theta azimuthal terms
        Bcl     = self.BcArr[ind_l]; Bcr  =   self.BcArr[ind_r]
        Bc      = Bcl + (r - r_left)*(Bcr-Bcl)/(r_right - r_left)
        
        Bsl     = self.BsArr[ind_l]; Bsr  =   self.BsArr[ind_r]
        Bs      = Bsl + (r - r_left)*(Bsr-Bsl)/(r_right - r_left)
        
        Gcl     = self.GcArr[ind_l]; Gcr  =   self.GcArr[ind_r]
        Gc      = Gcl + (r - r_left)*(Gcr-Gcl)/(r_right - r_left)
        
        Gsl     = self.GsArr[ind_l]; Gsr  =   self.GsArr[ind_r]
        Gs      = Gsl + (r - r_left)*(Gsr-Gsl)/(r_right - r_left)
        
        Hcl     = self.HcArr[ind_l]; Hcr  =   self.HcArr[ind_r]
        Hc      = Hcl + (r - r_left)*(Hcr-Hcl)/(r_right - r_left)
        
        Hsl     = self.HsArr[ind_l]; Hsr  =   self.HsArr[ind_r]
        Hs      = Hsl + (r - r_left)*(Hsr-Hsl)/(r_right - r_left)
        ## 4-theta azimuthal terms
        Ccl     = self.CcArr[ind_l]; Ccr  =   self.CcArr[ind_r]
        Cc      = Ccl + (r - r_left)*(Ccr-Ccl)/(r_right - r_left)
        
        Csl     = self.CsArr[ind_l]; Csr  =   self.CsArr[ind_r]
        Cs      = Csl + (r - r_left)*(Csr-Csl)/(r_right - r_left)
        
        return rho, A, C, F, L, N, Bc, Bs, Gc, Gs, Hc, Hs, Cc, Cs
    
    #########################################################################
    # functions related to dipping interface
    #########################################################################
    def init_dip_strike_interface(self):
        """initialize dip/strike angle array
        """
        self.dipping    = True
        self.dipifArr   = np.zeros(self.rArr.size, np.float32)
        self.strikeifArr= np.zeros(self.rArr.size, np.float32)
        return
    
    def get_dip_strike_interface(self, r, left):
        """
        Return dip/strike of dipping interface given a radius
        ===================================================================================
        ::: input parameters :::
        r           - radius (unit - m)
        left        - yield the LEFT value if repeated radius grid points appear or not
        ::: output :::
        dip/strike  - dip/strike (unit - degree) 
        ===================================================================================
        """
        if not self.is_layer_model():
            print 'WARNING: the model is not layrized!'
        # r   = np.float32( (6371.- z)*1000.)
        ind = np.where(self.rArr == r)[0]
        if ind.size == 0:
            if left:
                indds   = (np.where(self.rArr<r)[0])[-1]
                dip     = self.dipifArr[indds]; strike     = self.strikeifArr[indds]
            else:
                indds   = (np.where(self.rArr>r)[0])[0]
                dip     = self.dipifArr[indds]; strike     = self.strikeifArr[indds]
        else:
            if left:
                indds   = ind[0]
            else:
                indds   = ind[-1]
            dip     = self.dipifArr[indds]; strike     = self.strikeifArr[indds]
        return dip, strike
    
    
    
        
    
    
    
    
    
