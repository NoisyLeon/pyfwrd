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


# 

# m.add_perturb_layer(0, 200., 0, 4.6, False)
# m.add_perturb_layer(0, 200., 1, 4.6, False)
# m.add_perturb_layer(0, 200., 2, 8.0, False)
# m.add_perturb_layer(0, 200., 3, 8.0, False)
# m.add_perturb_layer(0, 200., 4, 1., False)
# m.add_perturb_layer(0, 200., 5, 3.3, False)
# 
# m.add_perturb_layer(0, 30., 0, 3.2, False)
# m.add_perturb_layer(0, 30., 1, 3.2, False)
# m.add_perturb_layer(0, 30., 2, 5.8, False)
# m.add_perturb_layer(0, 30., 3, 5.8, False)
# m.add_perturb_layer(0, 30., 4, 1.0, False)
# m.add_perturb_layer(0, 30., 5, 2.73, False)

# m.add_perturb_layer_love(0, 30., 0, 0.1, True)
# m.add_perturb_layer_love(0, 30., 1, -0.1, True)
m.add_perturb_layer_love(0, 20., 3, -0.05, True)
m.add_perturb_layer_love(0, 20., 4, 0.05, True)
# 
m.init_tilt()
m.dipArr[-1] = 90; m.dipArr[-2] = 90#; m.dipArr[-3] = 90; m.dipArr[-4] = 90
# m.dipArr[-1] = 10; m.dipArr[-2] = 10; m.dipArr[-3] = 10; m.dipArr[-4] = 10
m.strikeArr[-1] = 0; m.strikeArr[-2] = 0

# 
# 
rsolver1  = ref.ref_solver(m)
rsolver1.init_default_2()
rsolver1.solve_raysum(bazin=[0., 90])


ax=plt.subplot()

for i in xrange(rsolver1.bazArr.size):
    ax.plot(rsolver1.time, rsolver1.trROT[0,:,i], '--', lw=3)





plt.xlim([0., 20.])
plt.show()