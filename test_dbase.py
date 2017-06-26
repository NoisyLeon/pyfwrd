import database
import vmodel
m=vmodel.model1d()
m=vmodel.read_model(m, 'ak135.txt')
# m.earth_flattening()
dbase= database.eigenASDF()
dbase.init_dbase(inmodel=m, Tmax=100., dT=10., dc=50.,  nmodes=1, dz=1., zmax=1000.)
dbase.run(outfname='./ak135_cps.asdf')
dbase.write_disp(outfname='benchmark_cps/ak135_ray_ph.txt')
dbase.write_disp(outfname='benchmark_cps/ak135_ray_gr.txt', dtype='U')
dbase.write_disp(outfname='benchmark_cps/ak135_love_ph.txt', wavetype='love')
dbase.write_disp(outfname='benchmark_cps/ak135_love_gr.txt',wavetype='love', dtype='U')
# dbase.init_dbase(nmodes=1)
# dbase.load('./ak135_cps.asdf')