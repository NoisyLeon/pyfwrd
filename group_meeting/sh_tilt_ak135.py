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

m.add_perturb_layer_love(0, 20., 0, -0.02, True)
m.add_perturb_layer_love(0, 20., 1, +0.02, True)
m.add_perturb_layer_love(0, 20., 3, -0.02, True)
m.add_perturb_layer_love(0, 20., 4, +0.02, True)


tcpsR0 = tcps.tcps_solver(m)
tcpsR0.init_default()
tcpsR0.solve_PSV()


m.init_tilt()
m.dipArr[-1] = 60; m.dipArr[-2] = 60
m.strikeArr[-1] = 0; m.strikeArr[-2] = 0

m.rot_dip_strike()
m.decompose()

tcpsL = tcps.tcps_solver(m)
tcpsL.init_default_2()
tcpsL.solve_SH()
# 
CL1  = []
for baz in np.arange(36)*10.:
    tcpsL.sh_azi_perturb(baz, True)
    CL1.append(tcpsL.CA[1])
CL1 = np.array(CL1)    

####################################
ani  = aniproppy.aniprop_solver(m)
ani.init_default_2()

try:
    CL2 = np.loadtxt('azi_data_0.02/CL_dip_'+str(m.dipArr[-1])+'_strike_'+str(m.strikeArr[-1])+'.txt')
except:
    CL2  = []
    for baz in np.arange(36)*10.:
        print baz
        ani.solve_surf(az=baz)
        CL2.append(ani.CL[1])
    CL2 = np.array(CL2)
    np.savetxt('azi_data_0.02/CL_dip_'+str(m.dipArr[-1])+'_strike_'+str(m.strikeArr[-1])+'.txt', CL2)
##################################################


CL3  = []
for baz in np.arange(36)*10.:
    tcpsL.sh_azi_perturb(baz)
    CL3.append(tcpsL.CA[1])
CL3 = np.array(CL3)   

plt.figure()
ax=plt.subplot()
plt.plot(np.arange(36)*10., CL1, 'ro', ms=10, label='Montagner & Nataf, both theta-2 and theta-4 term')
plt.plot(np.arange(36)*10., CL2, 'b^', ms=8, label='J.Park')
plt.plot(np.arange(36)*10., CL3, 'kx', ms=10, label='Montagner & Nataf, only theta-2 term')
# plt.plot(np.arange(36)*10., tcpsL.C[1]*np.ones(36, dtype=np.float32), 'g-', ms=8, label='Montagner & Nataf, average')
# plt.plot(np.arange(36)*10., CL2.mean()*np.ones(36, dtype=np.float32), 'y-', ms=8, label='J.Park, average')
# plt.plot(np.arange(36)*10., CL3, 'kx', ms=10, label='Montagner & Nataf, only theta-2 term')
plt.xlabel('Azimuth (deg)', fontsize=30)
plt.ylabel('C (km/s)', fontsize=30)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
plt.legend(loc=0, fontsize=20)
plt.title('Love wave', fontsize=35)


plt.show()
