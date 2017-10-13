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



rsolver  = ref.ref_solver(m)
rsolver.init_default_2()
rsolver.solve_theo()
rsolver.solve_raysum()

ax=plt.subplot()
ax.plot(rsolver.time, rsolver.rf, 'k-')
tfill=rsolver.time[rsolver.rf>0]
yfill=rsolver.rf[rsolver.rf>0]
ax.fill_between(tfill, 0., yfill, color='red', linestyle='--', lw=0., alpha=0.3)

rfr  = rsolver.trROT[0,:,0]
ax.plot(rsolver.time, rfr, 'k--')
tfill=rsolver.time[rfr>0]
yfill=rfr[rfr>0]
ax.fill_between(tfill, 0., yfill, color='blue', linestyle='--', lw=0., alpha=0.3)

# tfill=rsolver.time[rsolver.rf<0]
# yfill=rsolver.rf[rsolver.rf<0]
# ax.fill_between(tfill, 0., yfill, color='blue', linestyle='--', lw=0.)

plt.show()