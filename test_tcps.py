import tcps
import vmodel
m=vmodel.model1d()
# m.get_radius(4000., 1.)
# m.model_prem()
m=vmodel.model1d()
m=vmodel.read_model(m, 'ak135.txt')
tcps1 = tcps.tcps_solver(m)
tcps1.init_default()