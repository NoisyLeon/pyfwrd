
import numpy as np
import eigen
import vmodel
import matplotlib.pyplot as plt
import asdf
    
class eigenASDF(asdf.AsdfFile):
    """
    ASDF database for eigenfunction/dispersion computation and data storage
    Note that the ASDF here is Advanced Scientific Data Format, NOT Adaptable Seismic Data Format !
    """
    def init_eigensolver(self, inmodel=None, T=None, Tmin=10., Tmax=50., dT=5.,
            c=None, cmin=2500., cmax=5000., dc=50., nmodes=1, zmax=400., dz=1.):
        #####
        # initialize eigenfunction solver
        #####
        if isinstance(inmodel, vmodel.model1d):
            eigensolver  = eigen.eigen_solver(inmodel)
        else:
            m=vmodel.model1d()
            m.get_radius(zmax, dz)
            m.model_prem()
            eigensolver  = eigen.eigen_solver(m)
        if not isinstance(T, np.ndarray):
            T = np.arange(Tmin, Tmax+dT, dT, dtype=np.float32)
        if not isinstance(c, np.ndarray):
            c = np.arange(cmin, cmax+dc, dc, dtype=np.float32)
        rmin    = np.float32(6371000. - zmax*1000.)
        dr      = np.float32(dz*1000.)
        eigensolver.init_dbase(T, c, rmin, dr, np.int32(nmodes))
        # if newsolver: 
        #     self.eigensolver.init_dbase(T, c, rmin, dr, np.int32(nmodes))
        # else:
        #     r = np.arange(rmin, 6371000.+dr, dr, dtype=np.float32)
        #     if r.size == self.eigensolver.r.size and T.size == self.eigensolver.T.size and \
        #         c.size == self.eigensolver.c.size and nmodes == self.eigensolver.nmodes:
        #         if not (np.allclose(r, self.eigensolver.r) and np.allclose(T, self.eigensolver.T) \
        #                 and np.allclose(c, self.eigensolver.c) and nmodes == self.eigensolver.nmodes):
        #             self.eigensolver.init_dbase(T, c, rmin, dr, nmodes)
        #     else:
        #         self.eigensolver.init_dbase(T, c, rmin, dr, nmodes)
        return eigensolver
    
    
    def init_dbase(self, inmodel=None, love=True, rayleigh=True, T=None, Tmin=10., Tmax=50., dT=5.,
            c=None, cmin=2000., cmax=5500., dc=50., nmodes=1, zmax=200., dz=1.):
        self.love       = love
        self.rayleigh   = rayleigh
        if not isinstance(T, np.ndarray):
            self.T = np.arange(Tmin, Tmax+dT, dT, dtype=np.float32)
        if not isinstance(c, np.ndarray):
            self.c = np.arange(cmin, cmax+dc, dc, dtype=np.float32)
        if love:
            self.eigenL=self.init_eigensolver(inmodel=inmodel, T=self.T, c=self.c, nmodes=nmodes, zmax=zmax, dz=dz)
        if rayleigh:
            self.eigenR=self.init_eigensolver(inmodel=inmodel, T=self.T, c=self.c, nmodes=nmodes, zmax=zmax, dz=dz)
        ###
        # Initialize output trees
        ###
        self.tree  = {'ray':{}, 'love': {} }
        return
        
    def run(self, allroot=True, outfname=None):
        if self.love:
            self.eigenL.solve_SH()
            if allroot:
                if not np.all(self.eigenL.eArr):
                    raise ValueError('Some roots not found for Love wave!')
            lovetree = {'r': self.eigenL.r, 'nmodes':self.eigenL.nmodes, 'T': self.eigenL.T, 'root': self.eigenL.eArr,
                        'disp': {'C': self.eigenL.Vph, 'U': self.eigenL.Vgr},\
                        'egn': {'l1': self.eigenL.l1data, 'l2': self.eigenL.l2data},
                        'kernel':{'ka': self.eigenL.Kadata, 'kc': self.eigenL.Kcdata, 'kl': self.eigenL.Kldata,\
                        'kf': self.eigenL.Kfdata, 'kn': self.eigenL.Kndata, 'krho0': self.eigenL.Krho0data,\
                        'kvph': self.eigenL.Kvphdata, 'kvpv': self.eigenL.Kvpvdata, 'kvsh': self.eigenL.Kvshdata,\
                        'kvsv': self.eigenL.Kvsvdata,'keta': self.eigenL.Ketadata, 'krho': self.eigenL.Krhodata} }
            self.tree['love'].update(lovetree)
        if self.rayleigh:
            self.eigenR.solve_PSV()
            if allroot:
                if not np.all(self.eigenR.eArr):
                    raise ValueError('Some roots not found for Rayleigh wave!')
            raytree = {'r': self.eigenR.r, 'nmodes':self.eigenR.nmodes, 'T': self.eigenR.T, 'root': self.eigenR.eArr,
                        'disp': {'C': self.eigenR.Vph, 'U': self.eigenR.Vgr},\
                        'egn': {'r1': self.eigenR.r1data, 'r2': self.eigenR.r2data, 'r3': self.eigenR.r3data, 'r4': self.eigenR.r4data},
                        'kernel':{'ka': self.eigenR.Kadata, 'kc': self.eigenR.Kcdata, 'kl': self.eigenR.Kldata,\
                        'kf': self.eigenR.Kfdata, 'kn': self.eigenR.Kndata, 'krho0': self.eigenR.Krho0data,\
                        'kvph': self.eigenR.Kvphdata, 'kvpv': self.eigenR.Kvpvdata, 'kvsh': self.eigenR.Kvshdata,\
                        'kvsv': self.eigenR.Kvsvdata,'keta': self.eigenR.Ketadata, 'krho': self.eigenR.Krhodata} }
            self.tree['ray'].update(raytree)
        if outfname!=None: self.write_to(outfname)
        return
    
    def load(self, infname):
        """
        Load ASDF file
        """
        self.tree.update((asdf.AsdfFile.open(infname)).tree)
    
    def plot_disp(self, wavetype='ray', style='o', mode=0, showfig=True):
        wavetype=wavetype.lower()
        if wavetype == 'rayleigh': wavetype='ray'
        root    = self.tree[wavetype]['root'][mode, :]
        if not np.all(root):
            print 'WARNING: Some root NOT found!'
        ind     = root.astype(bool)
        T       = self.tree[wavetype]['T'][ind]
        if T.size == 0:
            print 'WARNING: No datapoint!'
            return
        C       = self.tree[wavetype]['disp']['C'][mode, ind]/1000.
        U       = self.tree[wavetype]['disp']['U'][mode, ind]/1000.
        fig, ax = plt.subplots()
        plt.plot(T, C, style, lw=3, ms=10, label='C')
        plt.plot(T, U, style, lw=3, ms=10, label='U')
        plt.xlabel( 'Period (sec)', fontsize=30)
        plt.ylabel('Velocity (km/s)', fontsize=30)
        plt.title('Dispersion Curve', fontsize=40)
        ax.tick_params(axis='x', labelsize=20)
        ax.tick_params(axis='y', labelsize=20)
        plt.legend(numpoints=1, fontsize=20, loc=0)
        if showfig: plt.show()
        return
        
        