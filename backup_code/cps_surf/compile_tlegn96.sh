####################
## copy subs code
####################
cps_subs_dir='/home/leon/code/PROGRAMS.330/SUBS'

#cp $cps_subs_dir/f2csub.f .
#cp $cps_subs_dir/f96subf.f .
#cp $cps_subs_dir/igetmod.f .
cp $cps_subs_dir/tgetmod.f .
cp $cps_subs_dir/mnmarg.f .
cp $cps_subs_dir/mgtarg.f .
cp $cps_subs_dir/mchdep.f .
#cp $cps_subs_dir/solidf.f .
cp $cps_subs_dir/lgstr.f .
#cp $cps_subs_dir/grphsubf.f .


#f2py -h tegn96.pyf -m tegn96 tregn96_subroutine.f tlegn96_subroutine.f mnmarg.f mgtarg.f lgstr.f tio.f mchdep.f  tgetmod.f
#f2py --f77flags="-ffixed-line-length-none -O3" --f90flags="-O3" -c tegn96.pyf tregn96_subroutine.f tlegn96_subroutine.f mnmarg.f mgtarg.f lgstr.f tio.f mchdep.f tgetmod.f


f2py -h tlegn96.pyf -m tlegn96 tlegn96_subroutine.f mnmarg.f mgtarg.f lgstr.f tio.f mchdep.f  tgetmod.f
f2py --f77flags="-ffixed-line-length-none -O3" --f90flags="-O3" -c tlegn96.pyf tlegn96_subroutine.f mnmarg.f mgtarg.f lgstr.f tio.f mchdep.f tgetmod.f
