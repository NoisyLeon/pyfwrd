import numpy as np

inArr = np.loadtxt('ak135_dbase.txt')
H   = inArr[:, 0]
vpArr  = inArr[:, 1]
vsArr  = inArr[:, 2]
rhoArr = inArr[:, 3]

zbArr  = np.cumsum(H)
ztArr  = np.cumsum(H) - H
# z=np.array([])
outArr  = np.array([])
for i in xrange(H.size):
    # zb  = zbArr[i]
    zt  = ztArr[i]
    if i == 0:
        outArr = np.append(outArr, zt)
        outArr = np.append(outArr, vpArr[i])
        outArr = np.append(outArr, vsArr[i])
        outArr = np.append(outArr, rhoArr[i])
        continue
        # outArr = np.append(outArr, zt)
    if i==(H.size-1):
        outArr = np.append(outArr, zt)
        outArr = np.append(outArr, vpArr[i])
        outArr = np.append(outArr, vsArr[i])
        outArr = np.append(outArr, rhoArr[i])
        continue
    if vpArr[i-1] == vpArr[i] and rhoArr[i-1] == rhoArr[i] and vsArr[i-1] == vsArr[i]:
        outArr = np.append(outArr, zt)
        outArr = np.append(outArr, vpArr[i])
        outArr = np.append(outArr, vsArr[i])
        outArr = np.append(outArr, rhoArr[i])
        continue
    outArr = np.append(outArr, zt)
    outArr = np.append(outArr, vpArr[i-1])
    outArr = np.append(outArr, vsArr[i-1])
    outArr = np.append(outArr, rhoArr[i-1])
    outArr = np.append(outArr, zt)
    outArr = np.append(outArr, vpArr[i])
    outArr = np.append(outArr, vsArr[i])
    outArr = np.append(outArr, rhoArr[i])
# outArr = np.append()
outArr  = outArr.reshape(outArr.size/4, 4)
# outArr  = 
np.savetxt('ak135.txt',  outArr, fmt='%10.6f')