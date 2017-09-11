import numpy as np

ak135Arr = np.loadtxt('ak135_layer.txt')

hArr    = ak135Arr[:,0]
VsvArr     = ak135Arr[:,2]
VpvArr     = ak135Arr[:,1]
VshArr     = ak135Arr[:,2]
VphArr     = ak135Arr[:,1]
VpfArr     = np.sqrt ( (ak135Arr[:,1])**2 - 2.*((ak135Arr[:,2])**2) )
rhoArr     = ak135Arr[:,3]
zArr = hArr.cumsum()

zArr = np.repeat(zArr, 2)
VsArr = np.repeat(VsvArr, 2)
VpArr = np.repeat(VpvArr, 2)

zArr = np.append(0., zArr)
VsArr = np.append(VsvArr[0], VsArr)
VpArr = np.append(VpvArr[0], VpArr)
rhoArr= np.repeat(rhoArr, 2)
rhoArr= np.append(rhoArr[0], rhoArr)

outArr = np.append(zArr, VpArr)

outArr = np.append(outArr, VsArr)
outArr = np.append(outArr, rhoArr)

outArr=outArr.reshape( 4, VsArr.size,)
outArr=outArr.T
np.savetxt('ak135_cps.txt', outArr, fmt='%g')






