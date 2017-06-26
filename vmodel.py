
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
    # with open(infname, 'rb') as f:
    #     hdrstr  = f.readlines()
    #     lines   = hdrstr.split()
    #     if lines[0]=='#':
            
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
    ind = vsv!=0.
    vsv     = vsv[ind]
    vsh     = vsh[ind]
    vpv     = vpv[ind]
    vph     = vph[ind]
    eta     = eta[ind]
    rho     = rho[ind]
    radius  = radius[ind]
    
    model.get_data_vel(vsv, vsh, vpv, vph, eta, rho, radius)
    return model

spec = [('VsvArr', numba.float32[:]),
        ('VpvArr', numba.float32[:]),
        ('VshArr', numba.float32[:]),
        ('VphArr', numba.float32[:]),
        ('etaArr', numba.float32[:]),
        ('rhoArr', numba.float32[:]),
        ('rArr', numba.float32[:]),
        ('AArr', numba.float32[:]),
        ('CArr', numba.float32[:]),
        ('LArr', numba.float32[:]),
        ('FArr', numba.float32[:]),
        ('NArr', numba.float32[:])]

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
    ===============================================================================
    """
    def __init__(self):
        return
    
    def get_data_vel(self, vsv, vsh, vpv, vph, eta, rho, radius):
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
    
    def earth_flattening(self):
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
        # self.etaArr = eta
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
        Return rho, A, C, F, L, N for PREM (Dziewonski & Anderson, PEPI 1981) for 
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
        return
    
    # def model_ak135(self):
        
        
        
    def get_ind_Love_parameters(self, i):
        return self.rhoArr[i], self.AArr[i], self.CArr[i], self.FArr[i], self.LArr[i], self.NArr[i]
    
    def get_r_love_parameters(self, r):
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
        