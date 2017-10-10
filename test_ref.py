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


rsolver  = ref.ref_solver(m)
rsolver.init_default_2()


