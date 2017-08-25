import anisojp
import numpy as np
z = np.arange(10)*10
vp = np.array([5, 5, 5, 6, 6, 6, 7, 7, 7, 8])
vs = vp/1.8
vp2= 0.1 * np.ones(10)
vs2=0.1 * np.ones(10)

vp4 = 0.1 * np.ones(10)
rho = vp/2.0
theta=np.ones(10.)
phig=np.ones(10.)
nl = 10
baz = 0
nperiod=10

Rphase,Rgroup,Lphase,Lgroup,Period = anisojp.aniprop_interface(z,vp,vp2,vp4,vs,vs2,rho,theta,phig,nl,baz)