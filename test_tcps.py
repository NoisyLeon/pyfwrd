import tcps
import vmodel
import matplotlib.pyplot as plt

m=vmodel.model1d()
m.model_ak135_cps()
vs = 1.5

vp=0.9409+2.0947*vs-0.8206*vs**2+0.2683*vs**3-0.0251*vs**4 #Brocher Crust
rho=1.6612*vp-0.4721*vp**2+0.0671*vp**3-0.0043*vp**4+0.000106*vp**5

m.add_perturb_layer(0, 5., 0, vs, False)
m.add_perturb_layer(0, 5., 1, vs, False)
m.add_perturb_layer(0, 5., 2, vp, False)
m.add_perturb_layer(0, 5., 3, vp, False)
m.add_perturb_layer(0, 5., 5, rho, False)
m.add_perturb_layer(0, 5., 4, 1.0, False)
m.vel2love()

tcps1 = tcps.tcps_solver(m)
tcps1.init_default()
tcps1.solve_PSV()

m=vmodel.model1d()
m.model_ak135_cps()
vs = 1.5

vp=0.9409+2.0947*vs-0.8206*vs**2+0.2683*vs**3-0.0251*vs**4 #Brocher Crust
rho=1.6612*vp-0.4721*vp**2+0.0671*vp**3-0.0043*vp**4+0.000106*vp**5


m.add_perturb_layer(2., 5., 0, vs, False)
m.add_perturb_layer(2., 5., 1, vs, False)
m.add_perturb_layer(2., 5., 2, vp, False)
m.add_perturb_layer(2., 5., 3, vp, False)
m.add_perturb_layer(2., 5., 4, 1.0, False)
m.add_perturb_layer(2., 5., 5, rho, False)

m.add_perturb_layer(0., 2., 0, 0., False)
m.add_perturb_layer(0., 2., 1, 0., False)
m.add_perturb_layer(0., 2., 2, 1.5, False)
m.add_perturb_layer(0., 2., 3, 1.5, False)
m.add_perturb_layer(0., 2., 4, 1.0, False)
m.add_perturb_layer(0., 2., 5, 1.0, False)
m.vel2love()
tcps2 = tcps.tcps_solver(m)
tcps2.init_default()
tcps2.solve_PSV()

plt.plot(tcps1.T, tcps1.U, 'ro', ms=10)
plt.plot(tcps2.T, tcps2.U, 'bo', ms=10)
# plt.plot(tcps1.T, tcps1.C, 'rx', ms=10)
# plt.plot(tcps2.T, tcps2.C, 'bx', ms=10)
plt.show()
