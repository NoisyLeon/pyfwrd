
# #   0.1000E+01                dt 
# #     7   -1   -2             npts    n1(user spcified or not, when < 0: user specified)  n2
# # out.mod                     mname
# #          T         T        dolove  dorayl
# #   0.0000E+00  0.0000E+00    hs      hr
# #     1                       nmodes
# #   0.5000E+01  0.5000E+01    faclov  facrayl
# #   0.10000000E+00
# #   0.66666670E-01
# #   0.50000001E-01
# #   0.39999999E-01
# #   0.33333335E-01
# #   0.28571429E-01
# #   0.25000000E-01

# input
# ilvry : love(1) or rayleigh(2)
# nn - npts
# mname - model name
# verby - screen output or not
# nfval - npts
# fval  - input frequencies (array)
# ccmin, ccmax - min/max phase vel

# output 
# iret


#    subroutine disprs(ilvry,dt,nn,iret,mname,
# 1          verby,nfval,fval,ccmin,ccmax)




import tegn96
import vmodel
import numpy as np
model=vmodel.Model1d(modelindex=2)
model.read('test_model.txt')
d_in    = model.HArr
TA_in   = model.rhoArr * model.VphArr**2
TC_in   = model.rhoArr * model.VpvArr**2
TL_in   = model.rhoArr * model.VsvArr**2
TF_in   = 1.0 * (TA_in - np.float32(2.)* TL_in)
TN_in   = model.rhoArr * model.VshArr**2
TRho_in = model.rhoArr
qai_in  = model.QpArr
qbi_in  = model.QsArr
etapi_in= model.etapArr
etasi_in= model.etasArr
frefpi_in= model.frefpArr
frefsi_in= model.frefsArr

hs_in=0.
hr_in=0.
ohr_in=0.
ohs_in=0.
refdep_in=0.
nl_in = model.HArr.size
iflsph_in = 0
mode_in=np.ones(7)
Nt_in=7
t_in = (np.arange(7)*5.+10.)
cp_in = np.array([ 3.50507307,  3.64215541,  3.78464317,  3.92029595,  4.03891897,
        4.13657808,  4.21457863])
tegn96.tregn96(hs_in, hr_in, ohr_in,ohs_in, refdep_in, nl_in, iflsph_in,\
               d_in,TA_in,TC_in,TF_in,TL_in,TN_in,TRho_in,\
               qai_in,qbi_in,etapi_in,etasi_in,frefpi_in,frefsi_in, mode_in, t_in, Nt_in, cp_in)