!---------------------------------------------------------------------
program nms
!there are 4 configuarations files:
!-source.dat         : contains sources info+frequency dans, time step
!                      etc
!-recepteur.dat      : contains receivers info
!
!-nms.dat           : contains the normal modes catalogue information 
!                     and some more info
!---------------------------------------------------------------------
  use global_main_param, only: init_global_main
  use modes
  implicit none
  integer, parameter :: logunit=11
!
  open(logunit,file='nms.log')
  call init_global_main(logunit)
  call init_modes(logunit)
  call get_and_sum()
  close(logunit)
!---------------------------------------------------------------------
end program nms
!---------------------------------------------------------------------
