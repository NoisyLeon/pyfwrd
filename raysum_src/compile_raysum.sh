f2py -h raysum.pyf -m raysum buildmodel.f  eigenvec.f  eispack-cg.f  matrixops.f  misfit.f  phaselist.f  raysum.f  raysum_interface.f  readwrite.f trace.f params.h
f2py --f77flags="-ffixed-line-length-none -O3" --f90flags="-O3" -c raysum.pyf buildmodel.f  eigenvec.f  eispack-cg.f  matrixops.f  misfit.f  phaselist.f  raysum.f  raysum_interface.f  readwrite.f  trace.f params.h --fcompiler=gfortran 

