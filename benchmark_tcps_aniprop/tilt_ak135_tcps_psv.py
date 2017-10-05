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
    m.add_perturb_layer_love(0, 20., 4, 0.5, True), T = 5 sec, aniprop yields seemly wrong result. OTHERS are consistent ( ~ 0.05 %)
    m.add_perturb_layer_love(0, 20., 4, 0.2, True), T = 55 sec, aniprop yields seemly wrong result. OTHERS are consistent ( ~ 0.05 %)
    
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
    m.add_perturb_layer_love(0, 20., 3, -0.3, True) T = 5 sec, aniprop yields seemly wrong result. OTHERS are consistent ( ~ 0.05 %)
    m.add_perturb_layer_love(0, 20., 4, 0.3, True), T = 5 sec, aniprop yields seemly wrong result. OTHERS are consistent ( ~ 0.05 %)
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
tcpsR0 = tcps.tcps_solver(m)
tcpsR0.init_default()
tcpsR0.solve_PSV()

m.init_tilt()

m.dipArr[-1] = 45; m.dipArr[-2] = 45
m.strikeArr[-1] = 0.; m.strikeArr[-2] = 0.

m.rot_dip_strike()
m.decompose()
# 
tcpsR1 = tcps.tcps_solver(m)
tcpsR1.init_default()
tcpsR1.solve_PSV()
# 
CR1  = []
for baz in np.arange(360)*1.:
    tcpsR1.psv_azi_perturb(baz)
    CR1.append(tcpsR1.VphA[1])

CR2  = []
for baz in np.arange(360)*1.:
    tcpsR1.psv_azi_perturb_2theta(baz)
    CR2.append(tcpsR1.VphA[1])
    
m.init_tilt()
m.dipArr[-1] = 45.; m.dipArr[-2] = 45.
m.strikeArr[-1] = 90.; m.strikeArr[-2] = 90.

m.rot_dip_strike()
m.decompose()
# 
tcpsR2 = tcps.tcps_solver(m)
tcpsR2.init_default()
tcpsR2.solve_PSV()
# 
CR3  = []
for baz in np.arange(360)*1.:
    tcpsR2.psv_azi_perturb(baz)
    CR3.append(tcpsR2.VphA[1])

CR4  = []
for baz in np.arange(360)*1.:
    tcpsR2.psv_azi_perturb_2theta(baz)
    CR4.append(tcpsR2.VphA[1])


CR1 = np.array(CR1)
CR2 = np.array(CR2)

CR3 = np.array(CR3)
CR4 = np.array(CR4)

plt.plot(np.arange(360)*1., CR1, 'o', ms=5)
plt.plot(np.arange(360)*1., CR2, '^', ms=5)
plt.plot(np.arange(360)*1., CR3, 'o', ms=5)
plt.plot(np.arange(360)*1., CR4, '^', ms=5)
plt.plot(np.arange(360)*1., tcpsR2.Vph[1]*np.ones(360), 'x', ms=5)

plt.show()
