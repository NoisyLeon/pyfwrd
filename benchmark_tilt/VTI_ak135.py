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

tcpsR0 = tcps.tcps_solver(m)
tcpsR0.init_default()
tcpsR0.solve_PSV()

tcpsL0 = tcps.tcps_solver(m)
tcpsL0.init_default()
tcpsL0.solve_SH()
# # model perturbation: isotropic -> VTI

# m.add_perturb_layer_love(0, 20., 4, -0.1, True)
# m.add_perturb_layer_love(0, 20., 3, -0.3, True)

m.add_perturb_layer_love(0, 20., 3, -0.3, True)
tcpsR = tcps.tcps_solver(m)
tcpsR.init_default()
tcpsR.solve_PSV()

tcpsL = tcps.tcps_solver(m)
tcpsL.init_default()
tcpsL.solve_SH()
# 
ani  = aniproppy.aniprop_solver(m)
ani.init_default_2()
# ani.init_default(nl=100, dh=2.)
ani.solve_surf()

dcr = (tcpsR.Vph - ani.CR)/tcpsR.Vph*100.
dur = (tcpsR.Vgr - ani.UR)/tcpsR.Vgr*100.
dcl = (tcpsL.Vph - ani.CL)/tcpsL.Vph*100.
dul = (tcpsL.Vgr - ani.UL)/tcpsL.Vgr*100.

plt.figure()
plt.plot(tcpsR.T, tcpsR.Vph, 'ro', ms=10, label='herrmann: CR VTI')
plt.plot(ani.T, ani.CR, 'b^', ms=10, label='Park: CR VTI')
plt.plot(tcpsR0.T, tcpsR0.Vph, 'kx', ms=10, label='herrmann: CR iso')
plt.legend(loc=0, fontsize=15)

plt.figure()
plt.plot(tcpsR.T, tcpsR.Vgr, 'ro', ms=10, label='herrmann: UR VTI')
plt.plot(ani.T, ani.UR, 'b^', ms=10, label='Park: UR VTI')
plt.plot(tcpsR0.T, tcpsR0.Vgr, 'kx', ms=10, label='herrmann: UR iso')
plt.legend(loc=0, fontsize=15)

plt.figure()
plt.plot(tcpsL.T, tcpsL.Vph, 'ro', ms=10, label='herrmann: CL VTI')
plt.plot(ani.T, ani.CL, 'b^', ms=10, label='Park: CL VTI')
plt.plot(tcpsL0.T, tcpsL0.Vph, 'kx', ms=10, label='herrmann: CL iso')
plt.legend(loc=0, fontsize=15)

plt.figure()
plt.plot(tcpsL.T, tcpsL.Vgr, 'ro', ms=10, label='herrmann: UL VTI')
plt.plot(ani.T, ani.UL, 'b^', ms=10, label='Park: UL VTI')
plt.plot(tcpsL0.T, tcpsL0.Vgr, 'kx', ms=10, label='herrmann: UL iso')
plt.legend(loc=0, fontsize=15)
plt.show()

print np.abs(dcr).max(), np.abs(dur).max(), np.abs(dcl).max(), np.abs(dul).max()