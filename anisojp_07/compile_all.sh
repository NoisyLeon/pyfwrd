#f2py -h anisojp.pyf -m anisojp forward_solver.f90 refft.c aniprop_subroutines.f  eispack.f  rf_aniso_subroutines.f  util_subroutines.f
f2py --f77flags="-ffixed-line-length-none -O3" --f90flags="-O3" -c anisojp.pyf forward_solver.f90 refft.c aniprop_subroutines.f  eispack.f  rf_aniso_subroutines.f  util_subroutines.f --fcompiler=gfortran 
cp anisojp.so ..
#f2py -h anisojp.pyf -m anisojp refft.c aniprop_subroutines.f  eispack.f  rf_aniso_subroutines.f  util_subroutines.f
#f2py -c anisojp.pyf refft.c aniprop_subroutines.f  eispack.f  rf_aniso_subroutines.f  util_subroutines.f --f77flags=-ffixed-line-length-none --fcompiler=gfortran
