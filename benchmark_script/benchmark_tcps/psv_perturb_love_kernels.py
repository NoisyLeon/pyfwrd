import sys
sys.path.append('/home/leon/code/pysurf')

"""
Conclusion:

Flat Earth:
1. perturb top, 0-10 km, results are normal
2. perturb mid, 30-40 km, results are normal
    m.add_perturb_layer_love(30, 40., 0, -0.3, True), using velocity kernel is problematic. But Love kernel is correct.
3. perturb bottom, 180-200 km, results are normal
4. perturb all by 10%, results are normal
    
Spherical Earth:

Important Note:
From Herrmann's email:
"These partials refer to the original spherical model, not the model listed in the SRDER.TXT or SLDER.TXT.
If you focus on the layer numbers everything will work."

1. perturb top, 0-10 km, results are normal
2. perturb mid, 30-40 km, results are normal
    m.add_perturb_layer_love(30, 40., 0, -0.3, True), using velocity kernel is problematic. But Love kernel is correct.
    Possible reason:
    This case, A-2L = -0.20440674
3. perturb bottom, 180-200 km, results are normal
4. perturb all by 10%, results are normal

    


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
# m.add_perturb_layer_love(0, 10., 0, 0.1, True)
m.add_perturb_layer_love(30, 40., 0, -0.3, True)

m2 = vmodel.model1d()
m2.model_ak135_cps()
tcps2.init_default()
# tcps2.verbose=1
tcps2.solve_PSV()

tcps1.psv_perturb_disp_vel(tcps2)
tcps1.psv_perturb_disp_love(tcps2)


# v1 = tcps1.C_pre.copy()
# tcps1.model.flat=1
# tcps2.model.flat=1
# tcps1.love2vel()
# tcps2.love2vel()
# tcps1.psv_perturb_disp_vel(tcps2)
# v2 = tcps1.C_pre.copy()
plt.figure()
plt.plot(tcps1.T, tcps1.C, 'o', ms=10)
plt.plot(tcps1.T, tcps1.C_pre, 'y^', ms=10)
plt.plot(tcps2.T, tcps2.C, 'kx', ms=15)

plt.figure()
plt.plot(tcps1.T, tcps1.C, 'o', ms=10)
plt.plot(tcps1.T, tcps1.C_pre2, 'g^', ms=10)
plt.plot(tcps2.T, tcps2.C, 'kx', ms=15)

# plt.plot(tcps1.T, tcps1.C, 'o', ms=10)
# plt.plot(tcps1.T, v1, 'y^', ms=10)
# plt.plot(tcps2.T, tcps2.C, 'kx', ms=15)
print np.abs(tcps1.C_pre - tcps1.C_pre2)*1000.
plt.show()

# i=3
# plt.plot(tcps1.dArr.cumsum(), tcps1.dcdA[i, :], 'b-', lw=3, ms=10)
# plt.plot(tcps1.dArr.cumsum(), tcps1.dcdC[i, :], 'b--', lw=3,ms=10)
# plt.plot(tcps1.dArr.cumsum(), tcps1.dcdF[i, :], 'g-', lw=3,ms=10)
# plt.plot(tcps1.dArr.cumsum(), tcps1.dcdL[i, :], 'r--', lw=3,ms=10)
# plt.plot(tcps1.dArr.cumsum(), tcps1.dcdN[i, :], 'r-', lw=3,ms=10)
# plt.xlim(0, 100)
# plt.show()
