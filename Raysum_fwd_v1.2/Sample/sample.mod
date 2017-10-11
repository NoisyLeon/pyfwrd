# Sample model file. Lines starting with '#' are ignored.
# Layers are listed from top to bottom. The bottom layer is
# assumed to be a half-space. Interface strike and dip apply
# to the upper interface of the layer.
#
# Format:
#     Column   Contents
#        1    Thickness (m)
#        2    Density (kg/m^3)
#        3    Average P-wave velocity (m/s)
#        4    Average S-wave velocity (m/s)
#        5    Isotropic-layer flag (1:isotropic, 0:anisotropic)
#        6    %P anisotropy
#        7    %S anisotropy (if 5 and 6 are zero, isotropic layer)
#        8    Trend of fast axis (degrees)
#        9    Plunge of fast axis (degrees)
#       10    Interface strike (degrees)
#       11    Interface dip (degrees)  
#     Note that the percentages of anisotropy are peak-to-peak
#     (the expressions used are from Farra et al. (1991))
#
# Layers: crust, anisotropic wedge, isotropic half-space.
#thick rho  alph beta iso %P  %S  tr pl  st di
 40000 2600 6540 3710  1  0   0    0  0   0  0
100000 3500 7700 4200  0  5.5 5.5  0 10   0  0
     0 3500 8100 4500  1  0   0    0  0   0  10
