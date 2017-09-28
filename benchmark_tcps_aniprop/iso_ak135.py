import sys
sys.path.append('/home/leon/code/pysurf')

"""
Conclusion:

Flat Earth:
    T = 5 - 100. sec, dT = 5 sec
        Results are perfectly consistent (< ~ 0.02%) 

Spherical Earth:
    T = 5 - 100. sec, dT = 5 sec
        Results are consistent (< ~ 0.1%), result can be more consistent if nl is larger and dh is smaller. But this is computationally more expensive.

"""
import eigen, tcps, aniproppy
import vmodel
import numpy as np
import matplotlib.pyplot as plt

m=vmodel.model1d()
m.model_ak135_cps()
m.flat=0

tcpsR = tcps.tcps_solver(m)
tcpsR.init_default()
tcpsR.solve_PSV()

tcpsL = tcps.tcps_solver(m)
tcpsL.init_default()
tcpsL.solve_SH()
# 
ani  = aniproppy.aniprop_solver(m)
# ani.init_default(nl=50, dh=4.)
ani.init_default_2()
ani.solve_surf()


# 
# m.flat=1
# 
# tcpsRf = tcps.tcps_solver(m)
# tcpsRf.init_default(nl=50, dh=4)
# tcpsRf.solve_PSV()
# 
# 
# anif  = aniproppy.aniprop_solver(m)
# anif.init_default(nl=50, dh=4)
# anif.solve_surf()

dcr = (tcpsR.Vph - ani.CR)/tcpsR.Vph*100.
dur = (tcpsR.Vgr - ani.UR)/tcpsR.Vgr*100.
dcl = (tcpsL.Vph - ani.CL)/tcpsL.Vph*100.
dul = (tcpsL.Vgr - ani.UL)/tcpsL.Vgr*100.
print np.abs(dcr).max(), np.abs(dur).max(), np.abs(dcl).max(), np.abs(dul).max()