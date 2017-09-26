import sys
sys.path.append('/home/leon/code/pysurf')

"""
Conclusion:

Flat Earth:
1. perturb top, 0-10 km
    vsv -   m.add_perturb_layer(0, 10., 0, 0.8, True), no changes
            m.add_perturb_layer(0, 10., 0, -0.8, True), no changes for peturbed disp, direct computed has large change in dispersion
    others parameters:   results are normal
2. perturb mid, 30-40 km, results are normal
    vsv -   m.add_perturb_layer(30, 40., 0, -0.8, True), no changes for peturbed disp, direct computed has large change in dispersion
    others parameters:   results are normal
3. perturb bottom, 180-200 km
    m.add_perturb_layer(180, 200., 1, 0.2, True) is problematic for direct computation, but m.add_perturb_layer(180, 198., 1, 0.2, True) works
4. perturb all
    problematic cases for direct computation:
        m.add_perturb_layer(0, 200., 0, -0.1, True), m.add_perturb_layer(0, 198., 0, 0.4, True), m.add_perturb_layer(0, 198., 1, 0.1, True)
        
    related normal cases for direct computation:
        m.add_perturb_layer(0, 198., 0, -0.1, True), m.add_perturb_layer(0, 200., 0, 0.1, True), m.add_perturb_layer(0, 198., 1, 0.02, True)
    
Spherical Earth:

Important Note:
From Herrmann's email:
"These partials refer to the original spherical model, not the model listed in the SRDER.TXT or SLDER.TXT.
If you focus on the layer numbers everything will work."

1. perturb top, 0-10 km.
    similar to case 1 in flat Earth
2. perturb mid, 30-40 km.
    similar to case 2 in flat Earth
3. perturb bottom, 180-200 km
    similar to case 3 in flat Earth
4. perturb all
    similar to case 4 in flat Earth
    


"""
import eigen, tcps
import vmodel
import numpy as np
import matplotlib.pyplot as plt

m=vmodel.model1d()
m.model_ak135_cps()
m.flat=0

tcps1 = tcps.tcps_solver(m)
tcps1.init_default()
tcps1.solve_SH()

m.add_perturb_layer(0, 10., 0, 0.1, True)
tcps2 = tcps.tcps_solver(m)
tcps2.init_default()
# tcps2.verbose=1
tcps2.solve_SH()

tcps1.sh_perturb_disp_vel(tcps2)

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

# plt.plot(tcps1.T, tcps1.Vph, 'o', ms=10)
# plt.plot(tcps1.T, v1, 'y^', ms=10)
# plt.plot(tcps2.T, tcps2.Vph, 'kx', ms=15)

plt.show()

