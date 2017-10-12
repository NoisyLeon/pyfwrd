import raysum
import numpy as np
#thick rho  alph beta iso %P  %S  tr pl  st di
#  40000 2600 6540 3710  1  0   0    0  0   0  0
# 100000 3500 7700 4200  0  5.5 5.5  0 10   0  0
#      0 3500 8100 4500  1  0   0    0  0   0  10

thick   = np.zeros(15)
rho     = np.zeros(15)
alpha   = np.zeros(15)
beta    = np.zeros(15)

iso     = np.zeros(15)
dp      = np.zeros(15)
ds      = np.zeros(15)
trend   = np.zeros(15)
plunge  = np.zeros(15)
strike  = np.zeros(15)
dip     = np.zeros(15)

thick[:3]   = np.array([40, 100,0])*1000.
rho[:3]     = np.array([2.6, 3.5, 3.5]) * 1000.
alpha[:3]   = np.array([6.54, 7.7, 8.1]) * 1000.
beta[:3]   = np.array([3.71, 4.2, 4.5]) * 1000.

iso[:3]     = np.array([1,0,1])
dp[:3]      = np.array([0, 5.5, 0])
ds[:3]      = np.array([0, 5.5, 0])
trend[:3]   = np.zeros(3)
plunge[:3]  = np.array([0, 10, 0])
strike[:3]  = np.zeros(3)
dip[:3]     = np.array([0,0,10.])

iphase_in = 1
baz = np.zeros(200)
slow= np.zeros(200)
sta_dx=np.zeros(200)
sta_dy=np.zeros(200)
# baz[:12] = np.random.rand(12)*360.
baz[:12] = np.arange(12)*10.
slow[:12]=3e-5
ntr = 12


mults=1
nsamp=600
dt=0.05
width=1.
align=1
shift=5.
out_rot=0

phname_in=''

# tt[:nphase, :ntr]     
tt, amp, nphase, tr_cart, tr_ph = raysum.raysum_interface(3, thick, rho, alpha, beta, dp, ds, trend, plunge, strike, dip, iso, iphase_in,\
        ntr, baz, slow, sta_dx, sta_dy,\
        mults,nsamp,dt,width,align,shift,out_rot, phname_in)


