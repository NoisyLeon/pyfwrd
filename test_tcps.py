import tcps
import vmodel

m=vmodel.model1d()
m=vmodel.model1d()
m.model_ak135_cps()
tcps1 = tcps.tcps_solver(m)
tcps1.init_default()
tcps1.solve_PSV()