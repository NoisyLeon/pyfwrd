!------------------------------------------------------------------------
module yannos_flag
!------------------------------------------------------------------------
  implicit none
!
  logical :: force_fmin           =.false.
  logical :: cancel_gravity       =.false.
  logical :: never_use_startlevel =.false.
  logical :: use_tref             =.true.
  logical :: check_modes          =.false.
  logical :: use_remedy           =.true.
  logical :: rescue               =.true.
  logical :: restart              =.false.
  logical :: force_systemic_search=.false.
  logical :: keep_bad_modes       =.false.
!
!  character(len=3) :: modout_format='ucb'
  character(len=3) :: modout_format='ipg'
!  character(len=3) :: modout_format='olm'
!
  doubleprecision:: seuil_ray=1.d-4
  integer :: l_startlevel=10
! 
  logical  :: use_startlevel
  contains

!-----------------------------------------------------------------------
  subroutine read_flags(unit)
!-----------------------------------------------------------------------
    implicit none
    integer, intent(in) :: unit
!
    read(unit,*)
    read(unit,*)  force_fmin 
    read(unit,*)
    read(unit,*)  cancel_gravity 
    read(unit,*)
    read(unit,*)  never_use_startlevel
    read(unit,*)
    read(unit,*)  use_tref            
    read(unit,*)
    read(unit,*)  check_modes         
    read(unit,*)
    read(unit,*)  use_remedy          
    read(unit,*)
    read(unit,*)  rescue              
    read(unit,*)
    read(unit,*)  restart             
    read(unit,*)
    read(unit,*)  force_systemic_search
    read(unit,*)
    read(unit,*)  keep_bad_modes 
    read(unit,*)
    read(unit,'(a3)')  modout_format
    read(unit,*)
    read(unit,*) seuil_ray
    read(unit,*)
    read(unit,*) l_startlevel
!-----------------------------------------------------------------------
  end subroutine read_flags
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
end module yannos_flag
!------------------------------------------------------------------------
