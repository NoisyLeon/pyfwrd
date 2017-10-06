!  Forward code to compute surface wave dispersion and receiver function 
!  in anisotropic models based on aniprop.f and rfgen.f from Yale
!
!  Copyright (C) 2012 University of Texas at Austin
!  
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!  
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

program main 
    implicit none 

    integer,parameter:: nperiod=400 
    integer,parameter:: ntime=300 
    
    double precision,dimension(:),allocatable::theta,phig,z,vp,vp2,vp4, & 
                                               vs,vs2,rho 
    
    double precision,dimension(nperiod)::Period,Rphase,Rgroup,Lphase,Lgroup

    double precision,dimension(ntime):: Rrecfun,Trecfun,Timesave

    integer:: i, ios, nlayer 
    double precision:: baz 
    character(len=256) :: arg 

    call getarg(1,arg); read(arg,*,iostat=ios) baz; 
    if (ios/=0) stop 'Error reading back azimuth'

    ! read in model 
    print*, '##### Read in model #####' 
    open(1001,file='ZMODEL.txt',status='unknown')
    read(1001,*) nlayer
    
    allocate(z(nlayer+1))
    allocate(vp(nlayer+1))
    allocate(vp2(nlayer+1))
    allocate(vp4(nlayer+1))
    allocate(vs(nlayer+1))
    allocate(vs2(nlayer+1))
    allocate(rho(nlayer+1))
    allocate(theta(nlayer+1))
    allocate(phig(nlayer+1))

    do i=1,nlayer+1 
        read(1001,*) z(i),vp(i),vp2(i),vp4(i),vs(i),vs2(i),rho(i),theta(i),phig(i)
    end do 
    close(1001)

    print*, '##### Start forward calculation #####' 
    call forward_solver(z,vp,vp2,vp4,vs,vs2,rho,theta,phig,nlayer,baz,&
                        Rphase,Rgroup,Lphase,Lgroup,Period,nperiod,&
                        Rrecfun,Trecfun,Timesave,ntime) 

    print*, '##### Write out results #####'
    open(1001,file='ZDISP_Rayleigh.txt',status='unknown')
    do i=1,nperiod
        write(1001,*) Period(i),Rphase(i),Rgroup(i)
    end do 
    close(1001)

    open(1001,file='ZDISP_Love.txt',status='unknown')
    do i=1,nperiod 
        write(1001,*) Period(i),Lphase(i),Lgroup(i)
    end do 
    close(1001)

    open(1001,file='ZRecFun.txt',status='unknown')
    do i=1,ntime
        write(1001,'(F10.3,2x,F12.8,2x,F12.8,2x,F8.3)') Timesave(i),Rrecfun(i),Trecfun(i),baz
    end do 
    close(1001)

    deallocate(z)
    deallocate(vp)
    deallocate(vp2)
    deallocate(vp4)
    deallocate(vs)
    deallocate(vs2)
    deallocate(rho)
    deallocate(theta)
    deallocate(phig)

end program main 

! 2014-11-16 Hejun Zhu !  
