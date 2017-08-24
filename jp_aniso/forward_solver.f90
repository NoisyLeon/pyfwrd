!  Interfaces for forward calculation subroutines 
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


!/*< Interface for forward calculation >*/
subroutine forward_solver(z,vp,vp2,vp4,vs,vs2,rho,theta,phig,nl,baz,&
                          Rphase,Rgroup,Lphase,Lgroup,Period,Nperiod,&
                          Rrecfun,Trecfun,Timesave,Ntime) 
    implicit none 

    ! Input 
    integer:: nl 
    double precision,dimension(nl+1):: theta,phig,z,vp,vp2,vp4,vs,vs2,rho 
    double precision:: baz 

    ! Output 
    integer:: Nperiod
    double precision,dimension(Nperiod):: Rphase,Rgroup,Lphase,Lgroup,Period 

    ! Output 
    integer:: Ntime 
    double precision,dimension(Ntime):: Rrecfun,Trecfun,Timesave
   
    ! calculate dispersion curve
    print*,'##### Calculate dispersion curve #####' 
    call aniprop_interface(z,vp,vp2,vp4,vs,vs2,rho,theta,phig,nl,baz,&
                           Rphase,Rgroup,Lphase,Lgroup,Period,Nperiod)


    ! calculate receiver function
    print*, '##### Calculate receiver function #####' 
    call rf_aniso_interface(z,vp,vp2,vp4,vs,vs2,rho,theta,phig,nl,baz,&
                            Rrecfun,Trecfun,Timesave,Ntime)

end subroutine forward_solver 


!/*< Interface for surface wave dispersion >*/
subroutine aniprop_interface(z,vp,vp2,vp4,vs,vs2,rho,theta,phig,nl,baz,&
                             Rphase,Rgroup,Lphase,Lgroup,Period,nperiod)
    implicit none 

    ! Input 
    integer:: nl 
    double precision,dimension(nl+1):: theta,phig,z,vp,vp2,vp4,vs,vs2,rho 
    double precision:: baz 

    ! Output 
    integer:: nperiod
    double precision,dimension(nperiod):: Rphase,Rgroup,Lphase,Lgroup,Period 

    ! Output from aniprop 
    integer,parameter:: nperiod_out=1000
    double precision,dimension(nperiod_out):: Rphase_out,Rgroup_out, & 
                                              Lphase_out,Lgroup_out,Period_out 

    integer:: i,j,id 
    double precision:: pmin,pmax,dp,dist,distmin 

    pmin=10.0
    pmax=200.0 
    dp=(pmax-pmin)/(nperiod-1)
    do i=1,nperiod
        Period(i)=pmin+(i-1)*dp
    end do 

    call aniprop_subroutines(z,vp,vp2,vp4,vs,vs2,rho,theta,phig,nl,baz,& 
                             Rphase_out,Rgroup_out,Lphase_out,Lgroup_out,&
                             Period_out,nperiod_out) 

    ! Interpolate result             
    do i=1,nperiod  
        distmin=10000000000000000.0 
        do j=1,nperiod_out 
            dist=(Period_out(j)-Period(i))*(Period_out(j)-Period(i))
            if (dist < distmin) then 
                distmin=dist 
                id=j 
            end if 
        end do 
        Rphase(i)=Rphase_out(id)
        Rgroup(i)=Rgroup_out(id)
        Lphase(i)=Lphase_out(id)
        Lgroup(i)=Lgroup_out(id)
    end do 

end subroutine aniprop_interface 


!/*< Inerface for receiver function >*/ 
subroutine rf_aniso_interface(z,vp,vp2,vp4,vs,vs2,rho,theta,phig,nl,baz,&
                              Rrecfun,Trecfun,Timesave,Ntime)
    implicit none 

    ! Input 
    integer:: nl
    double precision,dimension(nl+1):: theta,phig,z,vp,vp2,vp4,vs,vs2,rho
    double precision:: baz 

    ! Output 
    integer:: Ntime 
    double precision,dimension(Ntime):: Rrecfun,Trecfun,Timesave

    call rf_aniso_subroutines(z,vp,vp2,vp4,vs,vs2,rho,theta,phig,nl,baz,&
                              Rrecfun,Trecfun,Timesave,Ntime)

end subroutine rf_aniso_interface 
