import sys
sys.path.append('/home/leon/code/pysurf')

"""
Conclusion:

Flat Earth:
    T = 5 - 100. sec, dT = 5 sec
    Cases has been benchmarked:
    m.add_perturb_layer_love(0, 20., 0-2, +-0.1, True) (~ 0.05 %)
    m.add_perturb_layer_love(0, 20., 3, -0.3, True) (~ 0.1 %)
    m.add_perturb_layer_love(0, 20., 4, 0.1, True) (~ 0.05 %)
    m.add_perturb_layer_love(0, 20., 4, -0.1, True); m.add_perturb_layer_love(0, 20., 3, -0.3, True) (group ~ 0.1 %, phase ~0.03 %)
    m.add_perturb_layer_love(0, 20., 5, +-0.3, True) (~ 0.01 %, ~ 0.02 % for Love group)
    
    Cases failed !
    m.add_perturb_layer_love(0, 20., 3, 0.05, True), for val > 0.05, aniprop yields an error!
    m.add_perturb_layer_love(0, 20., 4, 0.5, True), T = 5 sec, aniprop yields seemly wrong result. OTHERS are consistent ( ~ 0.05 %)
    m.add_perturb_layer_love(0, 20., 4, 0.2, True), T = 55 sec, aniprop yields seemly wrong result. OTHERS are consistent ( ~ 0.05 %)
    
    Bug summary:
    When L > N, aniprop crashes ! However, both A > C and C < A works ! Added aniprop_check_model function to vmodel object.
    
Spherical Earth:
    T = 5 - 100. sec, dT = 5 sec
    Cases has been benchmarked:
    m.add_perturb_layer_love(0, 20., 0-2, +-0.2, True) (~ 0.15 %, ~ 0.25 % for Love group)
    m.add_perturb_layer_love(0, 20., 3, -0.2, True) (~ 0.15 %, ~ 0.25 % for Love group)
    m.add_perturb_layer_love(0, 20., 4, 0.2, True) (~ 0.15 %, ~ 0.25 % for Love group)
    m.add_perturb_layer_love(0, 20., 4, -0.1, True); m.add_perturb_layer_love(0, 20., 3, -0.2, True) (~ 0.15 %, ~ 0.25 % for Love group)
    m.add_perturb_layer_love(0, 20., 5, +-0.3, True) (~ 0.15 %, ~ 0.25 % for Love group)
    

    Cases failed !
    m.add_perturb_layer_love(0, 20., 3, -0.3, True) T = 5 sec, aniprop yields seemly wrong result. OTHERS are consistent ( ~ 0.05 %)
    m.add_perturb_layer_love(0, 20., 4, 0.3, True), T = 5 sec, aniprop yields seemly wrong result. OTHERS are consistent ( ~ 0.05 %)
    m.add_perturb_layer_love(0, 20., 4, -0.1, True); m.add_perturb_layer_love(0, 20., 3, -0.3, True) T = 5 sec, wrong results ( ~ 0.05 %)
"""
import eigen, tcps, aniproppy, ref
import vmodel
import numpy as np
import matplotlib.pyplot as plt

m1=vmodel.model1d()
m1.model_ak135_cps()
m1.flat=1


m1.add_perturb_layer(0, 20., 0, 3.494, False)
m1.add_perturb_layer(0, 20., 1, 3.702, False)
m1.add_perturb_layer(0, 20., 2, 5.94, False)
m1.add_perturb_layer(0, 20., 3, 6.28, False)
m1.add_perturb_layer(0, 20., 4, 0.82, False)
m1.add_perturb_layer(0, 20., 5, 2.73, False)

m1.init_tilt()

m1.dipArr[-1] = 34; m1.dipArr[-2] = 34
m1.strikeArr[-1] = 20; m1.strikeArr[-2] = 20

m1.rot_dip_strike()
m1.decompose()
###########################################
m2=vmodel.model1d()
m2.model_ak135_cps()
m2.flat=1

m2.add_perturb_layer(0, 20., 0, 3.45, False)
m2.add_perturb_layer(0, 20., 1, 3.61, False)
m2.add_perturb_layer(0, 20., 2, 6.06, False)
m2.add_perturb_layer(0, 20., 3, 6.24, False)
m2.add_perturb_layer(0, 20., 4, 0.72, False)
m2.add_perturb_layer(0, 20., 5, 2.73, False)


m2.init_tilt()
m2.dipArr[-1] = 27; m2.dipArr[-2] = 27
m2.strikeArr[-1] = 110; m2.strikeArr[-2] = 110

m2.rot_dip_strike()
m2.decompose()



rsolver1  = ref.ref_solver(m1)
rsolver1.init_default_2()
for baz in np.arange(12)*30.:
    rsolver1.solve_anirec(baz=baz)
rsolver1.plot_baz_rf(comp='T', showfig=False)


rsolver2  = ref.ref_solver(m2)
rsolver2.init_default_2()
for baz in np.arange(12)*30.:
    rsolver2.solve_anirec(baz=baz)
rsolver2.plot_baz_rf(comp='T')





