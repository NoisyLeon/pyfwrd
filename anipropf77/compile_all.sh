#f2py -h anipropf77.pyf -m anipropf77 forward_solver.f90 refft.c aniprop_subroutines.f  eispack.f  rf_aniso_subroutines.f  util_subroutines.f
f2py --f77flags="-ffixed-line-length-none -O3" --f90flags="-O3" -c anipropf77.pyf forward_solver.f90 refft.c aniprop_subroutines.f  eispack.f  rf_aniso_subroutines.f  util_subroutines.f --fcompiler=gfortran 
cp anipropf77.so ..
