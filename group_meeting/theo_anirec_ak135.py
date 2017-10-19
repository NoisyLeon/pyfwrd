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
rsolver.solve_anirec()

ax=plt.subplot()
ax.plot(rsolver.time, rsolver.rf, 'g-', lw=4, label='theo, T. Shibutani')
time = rsolver.time
# ax.fill_between(time, y2=0., y1=rsolver.rf, where=rsolver.rf>0, color='green', lw=0.01, alpha=0.3, interpolate=True)
ax.plot(rsolver.time, rsolver.rfr, 'r--', lw=3, label='anirec, V.Levin & J. Park')
# ax.fill_between(time, y2=0., y1=rsolver.rfr, where=rsolver.rfr>0, color='red', lw=0.01, alpha=0.3, interpolate=True)

rsolver2  = ref.ref_solver(m)
rsolver2.init_default_2()
rsolver2.solve_raysum()

ax.plot(rsolver.time, rsolver2.rfrst[0]/rsolver2.rfrst[0].max()*rsolver.rf.max(), 'bo', lw=3, label='raysum, A. Fredriksen & M.Bostock')
# ax.fill_between(time, y2=0., y1=rsolver.rfr, where=rsolver.rfr>0, color='blue', lw=0.01, alpha=0.3, interpolate=True)

plt.xlabel('Time (sec)', fontsize=30)
# plt.ylabel('C (km/s)', fontsize=30)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
plt.title('Reference radial receiver function', fontsize=35)
plt.legend(loc=0, fontsize=20)
plt.xlim([0., 20.])
plt.show()