import sys
sys.path.append('/home/leon/code/pysurf')

"""

"""
import ref
import vmodel
import numpy as np
import matplotlib.pyplot as plt

m=vmodel.model1d()
m.model_ak135_cps()
m.flat=1
#

# m.add_perturb_layer(0, 30., 0, 3.494, False)
# m.add_perturb_layer(0, 30., 1, 3.702, False)
# m.add_perturb_layer(0, 30., 2, 5.94, False)
# m.add_perturb_layer(0, 30., 3, 6.28, False)
# m.add_perturb_layer(0, 30., 4, 0.82, False)
# m.add_perturb_layer(0, 30., 5, 2.7, False)
# 

m.add_perturb_layer(0, 200., 0, 4.6, False)
m.add_perturb_layer(0, 200., 1, 4.6, False)
m.add_perturb_layer(0, 200., 2, 8.0, False)
m.add_perturb_layer(0, 200., 3, 8.0, False)
m.add_perturb_layer(0, 200., 4, 1., False)
m.add_perturb_layer(0, 200., 5, 3.3, False)

m.add_perturb_layer(0, 30., 0, 3.2, False)
m.add_perturb_layer(0, 30., 1, 3.2, False)
m.add_perturb_layer(0, 30., 2, 5.8, False)
m.add_perturb_layer(0, 30., 3, 5.8, False)
m.add_perturb_layer(0, 30., 4, 1.0, False)
m.add_perturb_layer(0, 30., 5, 2.73, False)

# m.add_perturb_layer_love(0, 30., 0, 0.1, True)
# m.add_perturb_layer_love(0, 30., 1, -0.1, True)
m.add_perturb_layer_love(0, 30., 3, -0.05, True)
m.add_perturb_layer_love(0, 30., 4, 0.05, True)
# 
m.init_tilt()
m.dipArr[-1] = 90; m.dipArr[-2] = 90; m.dipArr[-3] = 90; m.dipArr[-4] = 90
# m.dipArr[-1] = 10; m.dipArr[-2] = 10; m.dipArr[-3] = 10; m.dipArr[-4] = 10
m.strikeArr[-1] = 0; m.strikeArr[-2] = 0

# 
# 
rsolver1  = ref.ref_solver(m)
rsolver1.init_default_3()
rsolver1.solve_anirec(baz=0.)
rsolver1.solve_anirec(baz=90.)


rsolver2  = ref.ref_solver(m)
rsolver2.init_default_3()
# rsolver2.dt = 0.025
rsolver2.solve_raysum(bazin=[0., 90.])
# 
#


m=vmodel.model1d()
m.model_ak135_cps()
m.flat=1
#

# m.add_perturb_layer(0, 30., 0, 3.494, False)
# m.add_perturb_layer(0, 30., 1, 3.702, False)
# m.add_perturb_layer(0, 30., 2, 5.94, False)
# m.add_perturb_layer(0, 30., 3, 6.28, False)
# m.add_perturb_layer(0, 30., 4, 0.82, False)
# m.add_perturb_layer(0, 30., 5, 2.7, False)
# 

m.add_perturb_layer(0, 200., 0, 4.6, False)
m.add_perturb_layer(0, 200., 1, 4.6, False)
m.add_perturb_layer(0, 200., 2, 8.0, False)
m.add_perturb_layer(0, 200., 3, 8.0, False)
m.add_perturb_layer(0, 200., 4, 1., False)
m.add_perturb_layer(0, 200., 5, 3.3, False)

m.add_perturb_layer(0, 30., 0, 3.2, False)
m.add_perturb_layer(0, 30., 1, 3.2, False)
m.add_perturb_layer(0, 30., 2, 5.8, False)
m.add_perturb_layer(0, 30., 3, 5.8, False)
m.add_perturb_layer(0, 30., 4, 1.0, False)
m.add_perturb_layer(0, 30., 5, 2.73, False)

rsolver3  = ref.ref_solver(m)
rsolver3.init_default_3()
rsolver3.solve_raysum()

ax=plt.subplot()
rf1 = rsolver1.rfrst[0]
ax.plot(rsolver1.time, rf1, 'r-', lw=3)


rf2  = rsolver2.rfrst[0]
rf2  = rf2/rf2.max()*rf1.max()
ax.plot(rsolver2.time, rf2, 'b--', lw=3)

rf4  = rsolver2.rfrst[1]
rf4  = rf4/rf4.max()*rf1.max()
ax.plot(rsolver2.time, rf4, 'go', lw=3)

rf3  = rsolver3.rfrst[0]
rf3  = rf3/rf3.max()*rf1.max()
ax.plot(rsolver1.time, rf3, 'ko', lw=3)

rf5  = rsolver1.rfrst[1]
rf5  = rf5/rf5.max()*rf1.max()
ax.plot(rsolver1.time, rf5, 'yo', lw=3)

plt.xlim([0., 20.])
plt.show()