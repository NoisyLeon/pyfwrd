import sws
import vmodel
import numpy as np
import matplotlib.pyplot as plt


m=vmodel.model1d()
m.model_ak135_cps()
m.flat=1
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

m.init_tilt()
m.dipArr[-1] = 90; m.dipArr[-2] = 90; m.dipArr[-3] = 90; m.dipArr[-4] = 90
# # m.dipArr[-1] = 10; m.dipArr[-2] = 10; m.dipArr[-3] = 10; m.dipArr[-4] = 10
m.strikeArr[-1] = 0; m.strikeArr[-2] = 0

ssolver = sws.sws_solver(m)
ssolver.init_default_2()
ssolver.solve_raysum(bazin=[35., -66., 233.])

time = ssolver.time
# plt.plot(time, ssolver.trROT[0,:,0], 'b')
# plt.plot(time, ssolver.trROT[1,:,0], 'r')
# plt.plot(time, ssolver.trROT[2,:,0], 'k')
# c1, c2 = ssolver.rotate(trid=0, angle=45)
# plt.plot(time, c1, 'y')
# plt.plot(time, c2, 'g')
# plt.show()


ssolver.convolve()

# plt.plot(np.arange(ssolver.stf.size)*ssolver.dt, ssolver.stf, 'b')
# plt.plot(time, ssolver.trSYN[0,:,0], 'r')
# plt.plot(time, ssolver.trSYN[1,:,0], 'k')
# c1, c2 = ssolver.rotate(trid=0, angle=-34.8, dtype=2)
# c3, c4 = ssolver.rotate(trid=0, angle=45, dtype=2)
# plt.plot(time, c1, 'y')
# plt.plot(time, c2, 'g')
# 
# plt.plot(time, c3, 'y--')
# plt.plot(time, c4, 'g--')
# 
# plt.show()

# ssolver.rotcorr_st()
# ssolver.ev_st(mtype=3)
ssolver.me_st()

# determine Radial definition
# baz = (33.+180.)/180.*np.pi
# 
# Tp= np.cos(baz)*ssolver.trNEZ[1,:,0] - np.sin(baz)*ssolver.trNEZ[0,:,0]
# Rp= np.sin(baz)*ssolver.trNEZ[1,:,0] + np.cos(baz)*ssolver.trNEZ[0,:,0]
# 
# 
# plt.plot(time, ssolver.trROT[0,:,0], 'k')
# plt.plot(time, ssolver.trROT[1,:,0], 'r')
# 
# # plt.plot(time, np.sqrt(ssolver.trNEZ[1,:,0]**2 +ssolver.trNEZ[0,:,0]**2))
# 
# plt.plot(time, Rp, 'ko')
# plt.plot(time, Tp, 'ro')
# 
# plt.show()


# 
# delayt, phi, C = ssolver.eigenvalue(trid=0, mtype=1)
# delayt, phi, C = ssolver.eigenvalue(trid=0, mtype=2)
# delayt, phi, C = ssolver.eigenvalue(trid=0, mtype=3)
# delayt, phi, C = ssolver.eigenvalue(trid=0, mtype=4)
# delayt, phi, C = ssolver.eigenvalue(trid=0, mtype=5)