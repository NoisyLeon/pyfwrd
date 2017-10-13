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







rsolver  = ref.ref_solver(m)
for az in np.arange(12)*30.+15.:
    rsolver.solve_anirec_reproduce(az=az)
rsolver.plot_az_rf(comp='R')
# 
# ax=plt.subplot()
# ax.plot(rsolver.time, rsolver.rf, 'k-')
# tfill=rsolver.time[rsolver.rf>0]
# yfill=rsolver.rf[rsolver.rf>0]
# ax.fill_between(tfill, 0., yfill, color='red', linestyle='--', lw=0., alpha=0.3)
# 
# ax.plot(rsolver.time, rsolver.rfr, 'k--')
# tfill=rsolver.time[rsolver.rfr>0]
# yfill=rsolver.rfr[rsolver.rfr>0]
# ax.fill_between(tfill, 0., yfill, color='blue', linestyle='--', lw=0., alpha=0.3)
# 
# # tfill=rsolver.time[rsolver.rf<0]
# # yfill=rsolver.rf[rsolver.rf<0]
# # ax.fill_between(tfill, 0., yfill, color='blue', linestyle='--', lw=0.)
# 
# plt.show()