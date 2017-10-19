import sys
sys.path.append('/home/leon/code/pysurf')

"""
Conclusion:

Flat Earth:
1. perturb top, 0-10 km, results are normal
    
2. perturb mid, 30-40 km, results are normal
    
3. perturb bottom, 180-200 km, results are normal
    problematic cases for direct computation:
        m.add_perturb_layer_love(180, 200., 4, 0.1, True)
    related normal cases for direct computation:
        m.add_perturb_layer_love(180, 200., 4, 0.05, True)
4. perturb all
    problematic cases for direct computation:
       m.add_perturb_layer_love(0, 200., 3, 0.5, True), m.add_perturb_layer_love(0, 200., 3, -0.2, True), m.add_perturb_layer_love(0, 200., 4, 0.1, True)
        
    related normal cases for direct computation:
        m.add_perturb_layer_love(0, 200., 3, 0.2, True), m.add_perturb_layer_love(0, 198., 3, -0.2, True), m.add_perturb_layer_love(0, 200., 4, 0.03, True)
    
Spherical Earth:

Important Note:
From Herrmann's email:
"These partials refer to the original spherical model, not the model listed in the SRDER.TXT or SLDER.TXT.
If you focus on the layer numbers everything will work."

1. perturb top, 0-10 km, results are normal

2. perturb mid, 30-40 km, results are normal

3. perturb bottom, 180-200 km
    problematic cases for direct computation:
        m.add_perturb_layer_love(180, 200., 4, 0.3, True)
    related normal cases for direct computation:
        m.add_perturb_layer_love(180, 200., 4, 0.2, True)
4. perturb all
    similar to case 4 in flat Earth
    but, m.add_perturb_layer_love(0, 198, 4, 0.1, True) is normal

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
tcps1.solve_SH()
# 
# m.add_perturb_layer_love(0, 20., 0, -0.1, True)
# m.add_perturb_layer_love(0, 20., 1, +0.1, True)
# m.add_perturb_layer_love(0, 20., 3, -0.1, True)
# m.add_perturb_layer_love(0, 20., 4, +0.1, True)

# m.add_perturb_layer_love(0, 20., 0, -0.05, True)
# m.add_perturb_layer_love(0, 20., 1, +0.05, True)
# m.add_perturb_layer_love(0, 20., 3, -0.05, True)
# m.add_perturb_layer_love(0, 20., 4, +0.05, True)

m.add_perturb_layer_love(0, 20., 0, -0.02, True)
m.add_perturb_layer_love(0, 20., 1, +0.02, True)
m.add_perturb_layer_love(0, 20., 3, -0.02, True)
m.add_perturb_layer_love(0, 20., 4, +0.02, True)


tcps2 = tcps.tcps_solver(m)
tcps2.init_default()
# tcps2.verbose=1
tcps2.solve_SH()

tcps1.sh_perturb_disp_vel(tcps2)
tcps1.sh_perturb_disp_love(tcps2)


plt.figure()
ax=plt.subplot()
plt.plot(tcps2.T, tcps2.C, 'ro-', ms=10, label='CL, directly computed')
plt.plot(tcps2.T, tcps1.C_pre2, 'b^-', ms=8, label='CL, computed from kernels')
plt.plot(tcps2.T, tcps1.C, 'kx-', ms=10, label='CL, reference')
plt.legend(loc=0, fontsize=15)
plt.xlabel('Period (sec)', fontsize=30)
plt.ylabel('C (km/s)', fontsize=30)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
plt.title('Love wave', fontsize=35)
print np.abs(tcps2.C - tcps1.C_pre2)*1000.
plt.show()

