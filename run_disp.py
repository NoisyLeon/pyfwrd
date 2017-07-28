import database
import vmodel

m=vmodel.model1d()
m=vmodel.read_axisem_bm(m, '1dmodel_ak135f.txt')
m.earth_flattening()
# m.trim(500.)
dbase= database.eigenASDF()
dbase.init_dbase(inmodel=m, Tmax=40, dT=2., dc=10.,  nmodes=1, dz=1., zmax=1000.)
dbase.run(outfname='./ak135f_1s.asdf')
dbase.write_disp(outfname='benchmark_axisem/ak135f_ray_ph.txt')
dbase.write_disp(outfname='benchmark_axisem/ak135f_ray_gr.txt', dtype='U')
dbase.write_disp(outfname='benchmark_axisem/ak135f_love_ph.txt', wavetype='love')
dbase.write_disp(outfname='benchmark_axisem/ak135f_love_gr.txt', wavetype='love', dtype='U')
# dbase.init_dbase(nmodes=1)
# dbase.load('./ak135_cps.asdf')