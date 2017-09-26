import sys
sys.path.append('/home/leon/code/pysurf')

"""
Conclusion:
A. Phase Velocity
1. Flat model
    eig : dr = 1000. m; rmin = 6071000 m
    tcps: h  = 1 km, zmax = 200 km
    Results : difference > 0.5 % when T > 55.
    
2. Spherical model
    eig : dr = 1000. m; rmin = 6071000 m
    tcps: h  = 1 km, zmax = 200 km
    Results : difference > 0.5 % when T > 55.

3. Spherical model
    coarse layer(tcps2.init_default_2()) of tcps is less accurate compared with eig
    finer layer of tcps is the most accurate
    
B. Group Velocity
    eigen is not accurate for short period, need correction !
"""
import eigen, tcps
import vmodel
import numpy as np
import matplotlib.pyplot as plt

m=vmodel.model1d()
m=vmodel.model1d()
m=vmodel.read_model(m, 'ak135.txt')
m.earth_flattening()
eig1 = eigen.eigen_solver(m)
eig1.init_default()
eig1.solve_PSV()


tcps1 = tcps.tcps_solver(m)
tcps1.init_default()
tcps1.solve_PSV()

tcps2 = tcps.tcps_solver(m)
tcps2.init_default_2()
tcps2.verbose=1
tcps2.solve_PSV()

plt.figure()
plt.plot((eig1.r), (eig1.r1data[0,0,:]), 'ro', ms=10)
plt.plot(tcps1.T, tcps1.Vph, 'b^', ms=10)

plt.figure()
plt.plot((eig1.T), (eig1.Vgr[0,:]/1000.), 'ro', ms=10)
plt.plot(tcps1.T, tcps1.Vgr, 'b^', ms=10)
plt.show()

