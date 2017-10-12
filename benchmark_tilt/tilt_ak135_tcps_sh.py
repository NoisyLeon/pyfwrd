import sys
sys.path.append('/home/leon/code/pysurf')

"""
Conclusion:

Flat Earth:
    T = 5 - 100. sec, dT = 5 sec
    Cases has been benchmarked:
    m.add_perturb_layer_love(0, 20., 0-2, +-0.1, True) (~ 0.05 %)
    m.add_perturb_layer_love(0, 20., 3, -0.3, True) (~ 0.1 %)
    m.add_perturb_layer_love(0, 20., 4, 0.1, True) (~ 0.05 %)
    m.add_perturb_layer_love(0, 20., 4, -0.1, True); m.add_perturb_layer_love(0, 20., 3, -0.3, True) (group ~ 0.1 %, phase ~0.03 %)
    m.add_perturb_layer_love(0, 20., 5, +-0.3, True) (~ 0.01 %, ~ 0.02 % for Love group)
    
    Cases failed !
    m.add_perturb_layer_love(0, 20., 3, 0.05, True), for val > 0.05, aniprop yields an error!
    m.add_perturb_layer_love(0, 20., 4, 0.5, True), T = 5 sec, aniprop yields seemly wrong result. OTHELS are consistent ( ~ 0.05 %)
    m.add_perturb_layer_love(0, 20., 4, 0.2, True), T = 55 sec, aniprop yields seemly wrong result. OTHELS are consistent ( ~ 0.05 %)
    
    Bug summary:
    When L > N, aniprop crashes ! However, both A > C and C < A works ! Added aniprop_check_model function to vmodel object.
    
Spherical Earth:
    T = 5 - 100. sec, dT = 5 sec
    Cases has been benchmarked:
    m.add_perturb_layer_love(0, 20., 0-2, +-0.2, True) (~ 0.15 %, ~ 0.25 % for Love group)
    m.add_perturb_layer_love(0, 20., 3, -0.2, True) (~ 0.15 %, ~ 0.25 % for Love group)
    m.add_perturb_layer_love(0, 20., 4, 0.2, True) (~ 0.15 %, ~ 0.25 % for Love group)
    m.add_perturb_layer_love(0, 20., 4, -0.1, True); m.add_perturb_layer_love(0, 20., 3, -0.2, True) (~ 0.15 %, ~ 0.25 % for Love group)
    m.add_perturb_layer_love(0, 20., 5, +-0.3, True) (~ 0.15 %, ~ 0.25 % for Love group)
    

    Cases failed !
    m.add_perturb_layer_love(0, 20., 3, -0.3, True) T = 5 sec, aniprop yields seemly wrong result. OTHELS are consistent ( ~ 0.05 %)
    m.add_perturb_layer_love(0, 20., 4, 0.3, True), T = 5 sec, aniprop yields seemly wrong result. OTHELS are consistent ( ~ 0.05 %)
    m.add_perturb_layer_love(0, 20., 4, -0.1, True); m.add_perturb_layer_love(0, 20., 3, -0.3, True) T = 5 sec, wrong results ( ~ 0.05 %)
"""
import eigen, tcps, aniproppy
import vmodel
import numpy as np
import matplotlib.pyplot as plt

m=vmodel.model1d()
m.model_ak135_cps()
m.flat=1
# 

# # model perturbation: isotropic -> VTI
# 
# # m.add_perturb_layer_love(0, 20., 4, -0.1, True)
# # m.add_perturb_layer_love(0, 20., 3, -0.3, True)
#
m.add_perturb_layer_love(0, 20., 0, -0.1, True)
m.add_perturb_layer_love(0, 20., 3, -0.1, True)
m.add_perturb_layer_love(0, 20., 4, 0.1, True)
tcpsL0 = tcps.tcps_solver(m)
tcpsL0.init_default()
tcpsL0.solve_SH()

m.init_tilt()

m.dipArr[-1] = 60; m.dipArr[-2] = 60
m.strikeArr[-1] = 0.; m.strikeArr[-2] = 0.

m.rot_dip_strike()
m.decompose()
# 
tcpsL1 = tcps.tcps_solver(m)
tcpsL1.init_default()
tcpsL1.solve_SH()
# 
CL1  = []
for baz in np.arange(360)*1.:
    tcpsL1.sh_azi_perturb(baz, True)
    CL1.append(tcpsL1.CA[1])

CL2  = []
for baz in np.arange(360)*1.:
    tcpsL1.sh_azi_perturb(baz)
    CL2.append(tcpsL1.CA[1])
    
m.init_tilt()
m.dipArr[-1] = 0; m.dipArr[-2] = 0
m.strikeArr[-1] = 0; m.strikeArr[-2] = 0

m.rot_dip_strike()
m.decompose()
# 
tcpsL2 = tcps.tcps_solver(m)
tcpsL2.init_default()
tcpsL2.solve_SH()
# 
CL3  = []
for baz in np.arange(360)*1.:
    tcpsL2.sh_azi_perturb(baz, True)
    CL3.append(tcpsL2.CA[1])

CL4  = []
for baz in np.arange(360)*1.:
    tcpsL2.sh_azi_perturb(baz)
    CL4.append(tcpsL2.CA[1])


CL1 = np.array(CL1)
CL2 = np.array(CL2)

CL3 = np.array(CL3)
CL4 = np.array(CL4)

plt.plot(np.arange(360)*1., CL1, 'ro', ms=5)
plt.plot(np.arange(360)*1., CL2, 'b^', ms=5)
# plt.plot(np.arange(360)*1., CL3, 'go', ms=5)
# plt.plot(np.arange(360)*1., CL4, 'k^', ms=5)
plt.plot(np.arange(360)*1., tcpsL1.C[1]*np.ones(360), 'x', ms=5)
plt.plot(np.arange(360)*1., tcpsL0.C[1]*np.ones(360), '-', ms=5)

plt.show()
