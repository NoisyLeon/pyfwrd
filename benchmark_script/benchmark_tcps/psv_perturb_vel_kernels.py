import sys
sys.path.append('/home/leon/code/pysurf')

"""
Conclusion:

Flat Earth:
1. perturb top, 0-10 km, results are normal
2. perturb mid, 30-40 km,
 results are normal
3. perturb bottom, 180-200 km, results are normal
4. perturb all by 10%
    For m.add_perturb_layer(0, 200., 2, 0.02, True), m.add_perturb_layer(0, 200., 3, -0.02, True) and m.add_perturb_layer(0, 200., 4, -0.1, True) based on ak135 model (zmax=200 km),
        tdisp96 cannot find proper dispersion solution!
    But change zmax = 198. in the input for m.add_perturb_layer works!
    For perturbation in density, the dispersion curve has almost no changes. This can be explained by the density kernels.
    
Spherical Earth:

Important Note:
From Herrmann's email:
"These partials refer to the original spherical model, not the model listed in the SRDER.TXT or SLDER.TXT.
If you focus on the layer numbers everything will work."

1. perturb top, 0-10 km.
    For m.add_perturb_layer(0, 10., 0, 0.05, True), m.add_perturb_layer(0, 10., 2, 0.05, True),
        m.add_perturb_layer(0, 10., 5, +0.05/-0.05, True)
    The perturbed dispersion curve using velocity difference BEFORE Earth flattening transformation
        is more close to directly computed curves from perturbed model.
    m.add_perturb_layer(0, 10., 3, 0.05, True), only the two shorest periods that the spherical vel perturbation yield MORE accurate results
    m.add_perturb_layer(0, 10., 3, -0.05, True), only the two shorest periods  that the spherical vel perturbation yield LESS accurate results
    m.add_perturb_layer(0, 10., 4, -0.5, True), only the some intermediate periods that the spherical vel perturbation yield MORE accurate results
    m.add_perturb_layer(0, 10., 4, 0.5, True), only the two shorest periods that the spherical vel perturbation yield LESS accurate results
    WHY???
2. perturb mid, 30-40 km.
    similar to case 1, others are normal
3. perturb bottom, 180-200 km
    similar to case 1, others are normal
4. perturb all by 10%
    Similar to former case 1 for spherical and case 4 for flat Earth
    


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


tcps2 = tcps.tcps_solver(m)
m.add_perturb_layer_love(0, 10., 5, 0.1, True)
tcps2.init_default()
# tcps2.verbose=1
tcps2.solve_PSV()

tcps1.psv_perturb_disp_vel(tcps2)

# v1 = tcps1.C_pre.copy()
# tcps1.model.flat=1
# tcps2.model.flat=1
# tcps1.love2vel()
# tcps2.love2vel()
# tcps1.psv_perturb_disp_vel(tcps2)
# v2 = tcps1.C_pre.copy()

plt.plot(tcps1.T, tcps1.C, 'o', ms=10)
plt.plot(tcps1.T, tcps1.C_pre, 'y^', ms=10)
plt.plot(tcps2.T, tcps2.C, 'kx', ms=15)

# plt.plot(tcps1.T, tcps1.C, 'o', ms=10)
# plt.plot(tcps1.T, v1, 'y^', ms=10)
# plt.plot(tcps2.T, tcps2.C, 'kx', ms=15)

plt.show()

