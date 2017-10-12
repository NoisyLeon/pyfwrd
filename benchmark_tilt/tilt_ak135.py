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
# m.add_perturb_layer_love(0, 20., 0, -0.05, True)
# m.add_perturb_layer_love(0, 20., 1, 0.05, True)
# m.add_perturb_layer_love(0, 20., 3, -0.05, True)
# m.add_perturb_layer_love(0, 20., 4, 0.05, True)

m.add_perturb_layer_love(0, 20., 0, -0.02, True)
m.add_perturb_layer_love(0, 20., 1, 0.02, True)
# m.add_perturb_layer_love(0, 20., 3, -0.05, True)
m.add_perturb_layer_love(0, 20., 4, 0.02, True)
tcpsR0 = tcps.tcps_solver(m)
tcpsR0.init_default()
tcpsR0.solve_PSV()


m.init_tilt()
m.dipArr[-1] = 90; m.dipArr[-2] = 90
m.strikeArr[-1] = 0; m.strikeArr[-2] = 0

m.rot_dip_strike()
# m.rot_dip_strike()
m.decompose()
#
# m.dipArr[-1] = 180; m.dipArr[-2] = 180
tcpsR = tcps.tcps_solver(m)
tcpsR.init_default_2()
tcpsR.solve_PSV()
# 
CR1  = []
for baz in np.arange(36)*10.:
    tcpsR.psv_azi_perturb(baz, True)
    CR1.append(tcpsR.CA[1])
CR1 = np.array(CR1)    

# # 
# plt.plot(np.arange(36)*10., CR1, '^', ms=5)
# plt.plot(np.arange(36)*10., tcpsR0.C[1]*np.ones(36), 'o', ms=5)
# plt.show()
# # 
# # tcpsL = tcps.tcps_solver(m)
# # tcpsL.init_default()
# # tcpsL.solve_SH()
# # #
# 
# # m.add_perturb_layer_love(0, 20., 0, -0.3, True)
ani  = aniproppy.aniprop_solver(m)
ani.init_default_2()
# # ani.init_default(nl=100, dh=2.)
# ani.solve_surf(baz=0.)
# print 'Start'
# ani.solve_surf(baz=80.)
# print 'End'
#
try:
    CR2 = np.loadtxt('azi_data/CR_dip_'+str(m.dipArr[-1])+'_strike_'+str(m.strikeArr[-1])+'.txt')
except:
    CR2  = []
    for baz in np.arange(36)*10.:
        print baz
        ani.solve_surf(az=baz)
        CR2.append(ani.CR[1])
    CR2 = np.array(CR2)
    np.savetxt('azi_data/CR_dip_'+str(m.dipArr[-1])+'_strike_'+str(m.strikeArr[-1])+'.txt', CR2)

CR3  = []
for baz in np.arange(36)*10.:
    tcpsR.psv_azi_perturb(baz)
    CR3.append(tcpsR.CA[1])
CR3 = np.array(CR3)   


plt.plot(np.arange(36)*10., CR1, 'ro', ms=5)
plt.plot(np.arange(36)*10., CR2, 'b^', ms=5)
plt.plot(np.arange(36)*10., CR3, 'kx', ms=5)


plt.plot(np.arange(36)*10., tcpsR.C[1]*np.ones(36), 'gx', ms=5)

cc = CR2[CR2<3.4]
# A = cc.max()-cc.min()
A = CR2.max()-CR2.min()
# CR4 = np.cos(2.*np.arange(36)*10./180.*np.pi)*A/2 + (CR2.max()+CR2.min())/2
# CR4 = np.cos(2.*np.arange(36)*10./180.*np.pi)*A/2 + (cc.max()+cc.min())/2
# plt.plot(np.arange(36)*10., CR4, '-', ms=5)
plt.show()
