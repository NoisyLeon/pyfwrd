import sys
sys.path.append('/home/leon/code/pysurf')

"""
Conclusion:

Flat Earth:
1. T = 5 - 100. sec, dT = 5 sec
    Results are perfectly consistent. 

"""
import eigen, tcps, aniproppy
import vmodel
import numpy as np
import matplotlib.pyplot as plt

m=vmodel.model1d()
m.model_ak135_cps()
m.flat=0

tcpsR = tcps.tcps_solver(m)
tcpsR.init_default(nl=50, dh=4)
tcpsR.solve_PSV()


tcpsL = tcps.tcps_solver(m)
tcpsL.init_default(nl=50, dh=4)
tcpsL.solve_SH()
# 
ani  = aniproppy.aniprop_solver(m)
ani.init_default(nl=50, dh=4)
ani.solve_surf()

dcr = (tcpsR.Vph - ani.CR)/tcpsR.Vph*100.
dur = (tcpsR.Vgr - ani.UR)/tcpsR.Vgr*100.
dcl = (tcpsL.Vph - ani.CL)/tcpsL.Vph*100.
dul = (tcpsL.Vgr - ani.UL)/tcpsL.Vgr*100.