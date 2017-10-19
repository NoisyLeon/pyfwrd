import vmodel

m=vmodel.model1d()
m.model_ak135_cps()
m.flat=1

# vmodel.plot(m, showfig=False)
# vmodel.plot(m, dtype='vpv', showfig=True)


vmodel.plot(m, dtype='rho', showfig=True)