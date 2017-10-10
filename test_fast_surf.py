import sys
sys.path.append('/home/leon/code/pysurf')

"""

"""
import fsurf, tcps, eigen
import vmodel
import numpy as np
import matplotlib.pyplot as plt

m=vmodel.model1d()
m.model_ak135_cps()
m.earth_flattening()
# 

# # model perturbation: isotropic -> VTI
# 
# # m.add_perturb_layer_love(0, 20., 4, -0.1, True)
# # m.add_perturb_layer_love(0, 20., 3, -0.3, True)
#
# m.add_perturb_layer_love(0, 20., 0, -0.05, True)
# m.add_perturb_layer_love(0, 20., 1, 0.05, True)
# m.add_perturb_layer_love(0, 20., 3, -0.05, True)
# m.add_perturb_layer_love(0, 20., 4, 0.05, True)

# m.add_perturb_layer_love(0, 20., 0, -0.02, True)
# m.add_perturb_layer_love(0, 20., 1, 0.02, True)
# # m.add_perturb_layer_love(0, 20., 3, -0.05, True)
# m.add_perturb_layer_love(0, 20., 4, 0.02, True)

ssolver  = fsurf.fsurf_solver(m)
ssolver.init_default_2()
ssolver.solve_surf()

tcpsR = tcps.tcps_solver(m)
tcpsR.init_default()
tcpsR.solve_PSV()

tcpsL = tcps.tcps_solver(m)
tcpsL.init_default()
tcpsL.solve_SH()


# eigR = eigen.eigen_solver(m)
# eigR.init_default()
# eigR.solve_PSV()
# 
# eigL = eigen.eigen_solver(m)
# eigL.init_default()
# eigL.solve_SH()


