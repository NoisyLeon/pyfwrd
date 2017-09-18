import sys
sys.path.append('/home/leon/code/pysurf')

"""
Conclusion:





"""
import eigen, tcps
import vmodel
import numpy as np
import matplotlib.pyplot as plt

m=vmodel.model1d()
m.model_ak135_cps()


tcps1 = tcps.tcps_solver(m)
tcps1.init_default()
tcps1.solve_PSV()

# 
tcps2 = tcps.tcps_solver(m)
tcps2.init_default(nl=100., dh=2.)
tcps2.verbose=1
tcps2.solve_PSV()
# 
# 
i=5
# plt.plot(tcps1.dArr.cumsum(), np.abs(tcps1.dcdah[i, :]-tcps1.dcdah1[i, :])/tcps1.dcdah[i, :], 'ro-', ms=10)
# plt.plot(tcps1.dArr.cumsum(), np.abs(tcps1.dcdah[i, :]-tcps1.dcdah2[i, :])/tcps1.dcdah[i, :], 'bo-', ms=10)

plt.plot(tcps1.dArr.cumsum(), tcps1.dcdah[i, :], 'ro-', ms=10)
plt.plot(tcps1.dArr.cumsum(), tcps1.dcdah1[i, :], 'bo-', ms=10)
plt.plot(tcps2.dArr.cumsum(), tcps2.dcdah[i, :]/2, 'go-', ms=5)
plt.plot(tcps2.dArr.cumsum(), tcps2.dcdah1[i, :], 'y^-', ms=5)

# #     
# # for i in xrange(10):
# #     # plt.figure()
# #     plt.plot((6371000. - eig1.r[::-1])/1000., (eig1.Kvphdata[0,i,::-1]), 'ro-', ms=10)
# #     plt.plot(tcps1.dArr.cumsum(), -tcps1.Kvphdata[i, :]*8., 'b^-', ms=10)
# # # plt.plot(tcps2.dArr.cumsum(), tcps2.r1data[5, :], 'kx-', ms=10)
# # # 
# # # plt.figure()
# # # plt.plot((eig1.T), (eig1.Vgr[0,:]/1000.), 'ro', ms=10)
# # # plt.plot(tcps1.T, tcps1.Vgr, 'b^', ms=10)
plt.show()

