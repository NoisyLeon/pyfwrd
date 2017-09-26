import numpy as np

"""
Check list:
ak135 model, period: 10 - 100 sec, zmax = 1000 km
1. Rayleigh phase: checked
2. Love phase: checked
3. Rayleigh group: problematic
4. Love group: checked

ak135 model spherical model, period: 10 - 100 sec, zmax = 1000 km
1. Rayleigh phase: checked
2. Love phase: checked
3. Rayleigh group: 
4. Love group: checked

"""
# inArr1  = np.loadtxt('benchmark_cps/ak135_cps_ray_ph.txt')
# inArr2  = np.loadtxt('benchmark_cps/ak135_ray_ph.txt')
# 
# inArr1  = np.loadtxt('benchmark_cps/ak135_cps_ray_gr.txt')
# inArr2  = np.loadtxt('benchmark_cps/ak135_ray_gr.txt')
# # # 
inArr1  = np.loadtxt('benchmark_cps/ak135_cps_love_ph.txt')
inArr2  = np.loadtxt('benchmark_cps/ak135_love_ph.txt')
# # 
# # # 
# inArr1  = np.loadtxt('benchmark_cps/ak135_cps_love_gr.txt')
# inArr2  = np.loadtxt('benchmark_cps/ak135_love_gr.txt')
v1  = inArr1[:, 1]
v2  = inArr2[:, 1]

print (v1-v2)/v1*100.