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

# m.add_perturb_layer(0, 35., 0, 3.2, False)
# m.add_perturb_layer(0, 35., 1, 3.2, False)
# m.add_perturb_layer(0, 35., 2, 5.8, False)
# m.add_perturb_layer(0, 35., 3, 5.8, False)
# m.add_perturb_layer(0, 35., 4, 1.0, False)
# m.add_perturb_layer(0, 35., 5, 2.73, False)



rsolver1  = ref.ref_solver(m)
rsolver1.init_default_2()
rsolver1.solve_anirec()

rsolver2  = ref.ref_solver(m)
rsolver2.init_default_2()
rsolver2.solve_raysum()


ax=plt.subplot()
rf1 = rsolver1.rfrst[0]
ax.plot(rsolver1.time, rf1, 'r-', lw=3)
# tfill=rsolver.time[rsolver.rf>0]
# yfill=rsolver.rf[rsolver.rf>0]
# ax.fill_between(tfill, 0., yfill, color='red', linestyle='--', lw=0., alpha=0.3)

# ax.plot(rsolver.time, rsolver.rfr, 'b-')


rf2  = rsolver2.rfrst[0]
rf2  = rf2/rf2.max()*rf1.max()

ax.plot(rsolver1.time, rf2, 'b--', lw=3)

# tfill=rsolver.time[rsolver.rf<0]
# yfill=rsolver.rf[rsolver.rf<0]
# ax.fill_between(tfill, 0., yfill, color='blue', linestyle='--', lw=0.)

plt.show()