#f2py -h aniprop.pyf -m aniprop forward_solver.f90 refft.c aniprop_subroutines.f  eispack.f  rf_aniso_subroutines.f util_subroutines.f
f2py --f77flags="-ffixed-line-length-none -O3" --f90flags="-O3" -c aniprop.pyf forward_solver.f90 refft.c aniprop_subroutines.f  eispack.f  rf_aniso_subroutines.f  util_subroutines.f --fcompiler=gfortran 
cp aniprop.so ..
