import vmodel
import numpy as np

# m=vmodel.model1d(VsvArr=np.array([0], dtype=np.float32), VpvArr=np.array([0], dtype=np.float32),
            # VshArr=np.array([0],dtype=np.float32), VphArr=np.array([0],dtype=np.float32), etaArr=np.array([0],dtype=np.float32), rhoArr=np.array([0],dtype=np.float32), isotropic=True)
m=vmodel.model1d()
# m.get_radius(1000., 1.)
# m=vmodel.read_model(m, 'ak135.txt')
m.model_ak135_cps()
