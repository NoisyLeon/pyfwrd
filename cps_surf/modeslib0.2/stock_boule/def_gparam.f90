module def_gparam
!module de parametres generaux. on peut ajouter, mais rien enlever
  integer, parameter :: P8 = selected_real_kind(8)
  integer, parameter :: DP = KIND(1.0D0)
  integer, parameter :: SP = KIND(1.0)
  integer, parameter :: I4 = selected_int_kind(9)
!
  real(P8),parameter :: pi=3.141592653589793_P8
  real(DP),parameter :: PI2=PI/2.0_DP,PI4=PI/4.0_DP
  real(DP),parameter :: deg2rad=PI/180._DP
  real(DP),parameter :: rad2deg=180._DP/PI
end module def_gparam
