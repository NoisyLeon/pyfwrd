
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


import tdisp96
import vmodel_cps
import numpy as np
model=vmodel_cps.Model1d(modelindex=2)
model.read('tregn96_src/test_model.txt')
d_in    = model.HArr
TA_in   = model.rhoArr * model.VphArr**2
TC_in   = model.rhoArr * model.VpvArr**2
TL_in   = model.rhoArr * model.VsvArr**2
TF_in   = 1.0 * (TA_in - np.float32(2.)* TL_in)
TN_in   = model.rhoArr * model.VshArr**2
TRho_in = model.rhoArr


ilvry = 1
dt = 1.
npts = 7
iret =1
verby = True
nfval = 7
fval  = 1./(np.arange(7)*5.+10.)
fval  = np.append(fval, np.zeros(2049-nfval))
ccmin = -1.
ccmax = -1.
nl_in = TRho_in.size
iflsph_in = 1 # 0 flat, 1: love, 2: rayleigh
refdep_in = 0.
nmode = 1
# 
# 
c_out,d_out,TA_out,TC_out,TF_out,TL_out,TN_out,TRho_out=tdisp96.disprs(ilvry,dt,npts,iret,verby, nfval,fval,ccmin,ccmax,\
               d_in,TA_in,TC_in,TF_in,TL_in,TN_in,TRho_in,\
               nl_in, iflsph_in, refdep_in, nmode, 0.5, 0.5)


import tregn96
import tlegn96
import vmodel_cps
import numpy as np
model=vmodel_cps.Model1d(modelindex=2)
model.read('tregn96_src/test_model.txt')
d_in    = d_out
TA_in   = TA_out
TC_in   = TC_out
TL_in   = TL_out
TF_in   = TF_out
TN_in   = TN_out
TRho_in = TRho_out
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
iflsph_in = 1
mode_in=np.ones(7, dtype=np.int)
Nt_in=7
t_in = (np.arange(7)*5.+10.)
cp_in = c_out[:7]
dogam=False



# u_out,ur, tur, uz, tuz,dcdh,dcdav,dcdah,dcdbv,dcdbh,dcdn,dcdr=tregn96.tregn96(hs_in, hr_in, ohr_in,ohs_in, refdep_in,dogam, nl_in, iflsph_in,\
#                d_in,TA_in,TC_in,TF_in,TL_in,TN_in,TRho_in,\
#                qai_in,qbi_in,etapi_in,etasi_in,frefpi_in,frefsi_in, Nt_in, t_in, cp_in)
# 
u_out,ut, tut,dcdh,dcdav,dcdah,dcdbv,dcdbh,dcdn,dcdr=tlegn96.tlegn96(hs_in, hr_in, ohr_in,ohs_in, refdep_in,dogam, nl_in, iflsph_in,\
               d_in,TA_in,TC_in,TF_in,TL_in,TN_in,TRho_in,\
               qai_in,qbi_in,etapi_in,etasi_in,frefpi_in,frefsi_in, Nt_in, t_in, cp_in)


