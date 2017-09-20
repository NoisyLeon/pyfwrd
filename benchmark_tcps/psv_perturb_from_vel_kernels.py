import sys
sys.path.append('/home/leon/code/pysurf')

"""
Conclusion:

Flat Earth:
1. perturb top, 0-10 km
2. perturb mid, 30-40 km
3. perturb bottome 10 km
4. perturb all by 10%
    

Spherical Earth:
1. perturb top, 0-10 km.
    The perturbed dispersion curve using velocity difference BEFORE Earth flattening transformation
        is more close to directly computed curves from perturbed model.
    WHY???
    


"""
import eigen, tcps
import vmodel
import numpy as np
import matplotlib.pyplot as plt

m=vmodel.model1d()
m.model_ak135_cps()
m.flat=1

tcps1 = tcps.tcps_solver(m)
tcps1.init_default()
tcps1.solve_PSV()

# continue with 4 !!!
tcps2 = tcps.tcps_solver(m)
m.add_perturb_layer(30., 40., 4, 0.1, True)
tcps2.init_default()
# tcps2.verbose=1
tcps2.solve_PSV()

tcps1.psv_perturb_disp_vel(tcps2)


# v1 = tcps1.Vph_pre.copy()
# tcps1.model.flat=1
# tcps2.model.flat=1
# tcps1.love2vel()
# tcps2.love2vel()
# tcps1.psv_perturb_disp_vel(tcps2)
# v2 = tcps1.Vph_pre.copy()

plt.plot(tcps1.T, tcps1.Vph, 'o', ms=10)
plt.plot(tcps1.T, tcps1.Vph_pre, 'y^', ms=10)
plt.plot(tcps2.T, tcps2.Vph, 'kx', ms=15)

plt.show()

