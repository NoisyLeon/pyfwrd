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

# m.add_perturb_layer_love(0, 20, 0, -0.02, True)
# m.add_perturb_layer_love(0, 20, 1, 0.02, True)
# # m.add_perturb_layer_love(0, 20., 3, -0.05, True)
# m.add_perturb_layer_love(0, 20, 4, 0.02, True)


tcpsR0 = tcps.tcps_solver(m)
tcpsR0.init_default()
tcpsR0.solve_PSV()


m.init_tilt()
m.dipArr[-1] = 60; m.dipArr[-2] = 60
m.strikeArr[-1] = 0; m.strikeArr[-2] = 0

m.rot_dip_strike()
m.decompose()

tcpsR = tcps.tcps_solver(m)
tcpsR.init_default_2()
tcpsR.solve_PSV()
# 
CR1  = []
for baz in np.arange(36)*10.:
    tcpsR.psv_azi_perturb(baz, True)
    CR1.append(tcpsR.CA[1])
CR1 = np.array(CR1)    

####################################
ani  = aniproppy.aniprop_solver(m)
ani.init_default_2()

try:
    CR2 = np.loadtxt('azi_data_0.02/CR_dip_'+str(m.dipArr[-1])+'_strike_'+str(m.strikeArr[-1])+'.txt')
except:
    CR2  = []
    for baz in np.arange(36)*10.:
        print baz
        ani.solve_surf(az=baz)
        CR2.append(ani.CR[1])
    CR2 = np.array(CR2)
    np.savetxt('azi_data_0.02/CR_dip_'+str(m.dipArr[-1])+'_strike_'+str(m.strikeArr[-1])+'.txt', CR2)
##################################################


CR3  = []
for baz in np.arange(36)*10.:
    tcpsR.psv_azi_perturb(baz)
    CR3.append(tcpsR.CA[1])
CR3 = np.array(CR3)   

plt.figure()
ax=plt.subplot()
plt.plot(np.arange(36)*10., CR1, 'ro', ms=10, label='Montagner & Nataf, both theta-2 and theta-4 term')
plt.plot(np.arange(36)*10., CR2+0.0012, 'b^', ms=8, label='J.Park')
plt.plot(np.arange(36)*10., CR3, 'kx', ms=10, label='Montagner & Nataf, only theta-2 term')
plt.xlabel('Azimuth (deg)', fontsize=30)
plt.ylabel('C (km/s)', fontsize=30)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
plt.title('Rayleigh wave', fontsize=35)
# plt.show()


# plt.plot(np.arange(36)*10., tcpsR.C[1]*np.ones(36), 'gx', ms=5)

# cc = CR2[CR2<3.4]
# A = cc.max()-cc.min()
# A = CR2.max()-CR2.min()
# CR4 = np.cos(2.*np.arange(36)*10./180.*np.pi)*A/2 + (CR2.max()+CR2.min())/2
# # CR4 = np.cos(2.*np.arange(36)*10./180.*np.pi)*A/2 + (cc.max()+cc.min())/2
# plt.plot(np.arange(36)*10., CR4, 'g-', ms=5, label='J.Park, theta-2 fit')
# # 
# A = CR1.max()-CR1.min()
# CR5 = np.cos(2.*np.arange(36)*10./180.*np.pi)*A/2 + (CR1.max()+CR1.min())/2
# # # CR4 = np.cos(2.*np.arange(36)*10./180.*np.pi)*A/2 + (cc.max()+cc.min())/2
# plt.plot(np.arange(36)*10., CR5, 'y-', ms=5, label='Montagner & Nataf, theta-2 fit')
plt.legend(loc=0, fontsize=15)
plt.show()
