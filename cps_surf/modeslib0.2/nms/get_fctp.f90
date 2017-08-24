!--------------------------------------------------------
program get_fctp
!--------------------------------------------------------
  use util_fctp
  implicit none
  character(len=100) :: file
  integer :: s,n,l,unit
  character(len=8) :: name

  unit=117
  print*,'Enter prefix of eigenfunction files:'
  read(*,'(a80)') file  
  print*,'Enter code (1=U,11,Up,2=V,22=Vp,3=W,33=Wp), n and l:'
  read*,s,n,l

  select case(s)
  case(1) 
     write(name,'(i3.3,"U",i3.3)') n,l
  case(11)
     write(name,'(i3.3,"Up",i3.3)') n,l
  case(2) 
     print*,'Warning, output will be multiplied by sqrt(l(l+1))'
     write(name,'(i3.3,"V",i3.3)') n,l
  case(22)
     print*,'Warning, output will be multiplied by sqrt(l(l+1))'
     write(name,'(i3.3,"Vp",i3.3)') n,l
  case(3) 
     print*,'Warning, output will be multiplied by sqrt(l(l+1))'
     write(name,'(i3.3,"W",i3.3)') n,l
  case(33)
     print*,'Warning, output will be multiplied by sqrt(l(l+1))'
     write(name,'(i3.3,"Wp",i3.3)') n,l
  end select
  open(unit,file=name,status='replace')
  if (s.eq.3.or.s.eq.33) then  
     call open_fctp1(file,'T')
     call write_modeT_ascii(unit,s,n,l)
  else
     call open_fctp1(file,'S') 
     call write_modeS_ascii(unit,s,n,l)
  endif
  close(unit)
  print*,'Ouput has been written in : ',name
!--------------------------------------------------------
end program get_fctp
!--------------------------------------------------------
