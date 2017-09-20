import tcps
import vmodel
import numpy as np
import matplotlib.pyplot as plt

m=vmodel.model1d()
m.model_ak135_cps()
m.flat=1

tcps1 = tcps.tcps_solver(m)
tcps1.init_default()
tcps1.solve_PSV()

m=vmodel.model1d()
m.model_ak135_cps()
m.flat=0
tcps2 = tcps.tcps_solver(m)
tcps2.init_default()
tcps2.solve_PSV()

