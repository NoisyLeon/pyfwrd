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

m1=vmodel.model1d()
m1.model_ak135_cps()
m1.flat=1
#
# m1.add_perturb_layer(0, 20., 0, 3.57, False)
# m1.add_perturb_layer(0, 20., 1, 3.74, False)
# m1.add_perturb_layer(0, 20., 2, 6.14, False)
# m1.add_perturb_layer(0, 20., 3, 6.52, False)
# m1.add_perturb_layer(0, 20., 4, 0.87, False)
# m1.add_perturb_layer(0, 20., 5, 2.79, False)

m1.add_perturb_layer(0, 20., 0, 3.494, False)
m1.add_perturb_layer(0, 20., 1, 3.702, False)
m1.add_perturb_layer(0, 20., 2, 5.94, False)
m1.add_perturb_layer(0, 20., 3, 6.28, False)
m1.add_perturb_layer(0, 20., 4, 0.82, False)
m1.add_perturb_layer(0, 20., 5, 2.73, False)

# m1.add_perturb_layer(0, 20., 0, 0.003, True)
# m1.add_perturb_layer(0, 20., 1, 0.02, True)



m1.init_tilt()

m1.dipArr[-1] = 34; m1.dipArr[-2] = 34
m1.strikeArr[-1] = 20; m1.strikeArr[-2] = 20

m1.rot_dip_strike()
m1.decompose()
###########################################
m2=vmodel.model1d()
m2.model_ak135_cps()
m2.flat=1
#
# m2.add_perturb_layer(0, 20., 0, 3.54, False)
# m2.add_perturb_layer(0, 20., 1, 3.72, False)
# m2.add_perturb_layer(0, 20., 2, 6.15, False)
# m2.add_perturb_layer(0, 20., 3, 6.47, False)
# m2.add_perturb_layer(0, 20., 4, 0.74, False)
# m2.add_perturb_layer(0, 20., 5, 2.79, False)

m2.add_perturb_layer(0, 20., 0, 3.45, False)
m2.add_perturb_layer(0, 20., 1, 3.61, False)
m2.add_perturb_layer(0, 20., 2, 6.06, False)
m2.add_perturb_layer(0, 20., 3, 6.24, False)
m2.add_perturb_layer(0, 20., 4, 0.72, False)
m2.add_perturb_layer(0, 20., 5, 2.73, False)


m2.init_tilt()
m2.dipArr[-1] = 27; m2.dipArr[-2] = 27
m2.strikeArr[-1] = 110; m2.strikeArr[-2] = 110

m2.rot_dip_strike()
m2.decompose()

# 
tcpsR1 = tcps.tcps_solver(m1)
tcpsR1.init_default()
tcpsR1.solve_PSV()

# 
tcpsR2 = tcps.tcps_solver(m2)
tcpsR2.init_default()
tcpsR2.solve_PSV()

# 
CR1  = []
for baz in np.arange(36)*10.:
    tcpsR1.psv_azi_perturb(baz, True)
    CR1.append(tcpsR1.CA[1])

CR2  = []
for baz in np.arange(36)*10.:
    tcpsR1.psv_azi_perturb(baz)
    CR2.append(tcpsR1.CA[1])
    


# 
CR3  = []
for baz in np.arange(36)*10.:
    tcpsR2.psv_azi_perturb(baz, True)
    CR3.append(tcpsR2.CA[1])

CR4  = []
for baz in np.arange(36)*10.:
    tcpsR2.psv_azi_perturb(baz)
    CR4.append(tcpsR2.CA[1])


CR1 = np.array(CR1)
CR2 = np.array(CR2)

CR3 = np.array(CR3)
CR4 = np.array(CR4)

plt.figure()
ax=plt.subplot()
plt.plot(np.arange(36)*10., CR2, 'ro', ms=10, label='group 1: dip=34, strike=20')
plt.plot(np.arange(36)*10., CR4, 'b^', ms=8, label='group 2: dip=27, strike=110')
# plt.plot(np.arange(36)*10., CR3, 'kx', ms=10, label='Montagner & Nataf, only theta-2 term')
plt.xlabel('Azimuth (deg)', fontsize=30)
plt.ylabel('C (km/s)', fontsize=30)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
plt.title('Rayleigh wave', fontsize=35)
plt.legend(loc=0, fontsize=15)
plt.show()
