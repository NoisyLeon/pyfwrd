import anisojp
import numpy as np
z = np.arange(10., dtype=np.float64)*10*1000.
vp = np.array([5, 5, 5, 6, 6, 6, 7, 7, 7, 8])*1000.
vs = vp/1.8
vp2= 0.1 * np.ones(10.)
vs2=0.1 * np.ones(10.)

vp4 = 0.1 * np.ones(10.)
rho = vp/2.0
theta=np.zeros(10.)
phig=np.zeros(10.)
nl = 9
baz = 0
nperiod=10

# Rphase,Rgroup,Lphase,Lgroup,Period = anisojp.aniprop_interface(z,vp,vp2,vp4,vs,vs2,rho,theta,phig,nl,baz,nperiod)

# Rphase,Rgroup,Lphase,Lgroup,Period = anisojp.aniprop_subroutines(z,vp,vp2,vp4,vs,vs2,rho,theta,phig,nl,baz,nperiod)

# Rphase,Rgroup,Lphase,Lgroup,Period = anisojp.test(z,nl, 10)