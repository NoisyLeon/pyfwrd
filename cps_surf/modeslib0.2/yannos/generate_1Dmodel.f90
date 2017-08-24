!---------------------------------------------------------------------
program generate_1Dmodel
!---------------------------------------------------------------------
  implicit none
  character(len=20) :: filename
  integer :: nb_couches,nb_zones
  real, dimension(:), allocatable :: r,rho,vpv,vph,vsv,vsh,eta,qkappa,qshear
  integer :: i,i_zone,unit=134,ii,Di,dense,ocb,icb,oce,ic,nz,nbz,iz,nbp,n1,n2
  real :: Dr,c,rn,rm,dr1,dr2,a
  real, dimension(:), allocatable :: r_inter,rho_z,vp_z,vs_z,qmu_z
  integer, dimension(:), allocatable :: i_inter
  real :: r_min
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! input
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  print*,'Entrer le nom du model'
  read*,filename
  print*,'Entrer le nombre de couches du modele:'
  read*,nb_couches
  print*,'Entrer le nombre de parties:'
  read*,nb_zones
  allocate(rho_z(nb_zones),vp_z(nb_zones),vs_z(nb_zones) &
          ,r_inter(0:nb_zones),i_inter(0:nb_zones),qmu_z(nb_zones))
  print*,'Entrer le rayon de depart:'
  read*,r_min
  r_inter(0)=r_min  
  do nbz=1,nb_zones
     write(*,100) nbz 
     read*,r_inter(nbz)
     write(*,101) nbz
     read*,rho_z(nbz)
     write(*,102) nbz
     read*,vp_z(nbz)
     write(*,103) nbz
     read*,vs_z(nbz)
     write(*,104) nbz
     read*,qmu_z(nbz)
  enddo
  print*,'Entrer le type de densification au interfaces:'
  print*,' nulle: 1'
  print*,' soft : 2'
  print*,' Hard : 3'
  read*, dense
  print*,'Le rayon minimum a ete fixe a (en km):',r_min
!
100 format('enter le rayon (km) de la zone ',i2,' (en partant du centre)')
101 format('enter la masse volumique (kg/m3)  de la zone ',i2)
102 format('enter la vp (km/s) de la zone ',i2)
103 format('enter la vs (km/s) de la zone ',i2)
104 format('enter Qmu          de la zone ',i2)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! fin input
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  allocate(r(0:nb_couches),rho(nb_couches),vpv(nb_couches),vph(nb_couches) &
       ,vsv(nb_couches),vsh(nb_couches),eta(nb_couches),qkappa(nb_couches) &
       ,qshear(nb_couches))
!locate indice inter:
  rn=r_inter(nb_zones)
  r(0)=r_min
  i_inter(0)=0
  do iz=1,nb_zones
     i_inter(iz)=int(r_inter(iz)/rn*nb_couches)
  enddo
!rayon
  select case(dense)
!--------------------
  case(1)
!cas lineaire
!-------------------
     do nz=1,nb_zones
        nbp=i_inter(nz)-i_inter(nz-1)
        dr=(r_inter(nz)-r_inter(nz-1))/(nbp-1)
        print*,'dr=',dr
        do ic=1,nbp
           i=ic+i_inter(nz-1)
           r(i)=dr*(ic-1)+r(i_inter(nz-1))
        enddo
     enddo
  case (2)
!cas amplifier par 
     rm=0.05*rn
     C=100.
     do nz=1,nb_zones
        nbp=i_inter(nz)-i_inter(nz-1)
        a=r_inter(nz)-r_inter(nz-1)
        dr2=(2.*rm*C+a-2.*rm)/(nbp-1)
        dr1=dr2/C
        print*,'dr1=',dr1,'dr2=',dr2
        n1=int((rm+1.e-10)/dr1)+1
        n2=int((a-2*rm+1.e-10)/dr2)+1
        print*,'n1=',n1,'n2=',n2
        dr=dr1
        do ic=1,nbp
           i=ic+i_inter(nz-1)
           if (ic<n1.or.ic>n1+n2-1) then
              dr=dr1
           else
              dr=dr2
           endif
           if (ic==1) then
              r(i)=r(i-1)
           else
              r(i)=r(i-1)+dr
           endif
write(101,*) r(i),dr
        enddo
     enddo
  case default
     stop'cas de dense non prevu!'
  end select
  do iz=1,nb_zones
     do ii=i_inter(iz-1)+1,i_inter(iz)
        rho(ii)   =rho_z(iz)
        vpv(ii)   =vp_z (iz)*1000.
        vph(ii)   =vp_z (iz)*1000.
        vsv(ii)   =vs_z (iz)*1000.
        vsh(ii)   =vs_z (iz)*1000.
        eta(ii)   =1.0
        qkappa(ii)=1.E5
        qshear(ii)=Qmu_z(iz)
     enddo
  enddo
!
  nz=0
  icb=0
  ocb=0
  oce=0
  do i=1,nb_zones
     if (nz>=2) stop'+ de 2 couches d''eau : ca va pas aller'
     if (abs(vs_z(i))<1.E-10) then
        nz=nz+1
        if (nz/=1.and.i==nb_zones) then
           oce=i_inter(i)-i_inter(i-1)
        else
           icb=i_inter(i-1)
           ocb=i_inter(i)
        endif
     endif
  enddo
!
  r(:)=r(:)*1000.
  open(unit,file=filename)
  write(unit,*) filename
  write(unit,*) 1,1.,1
  write(unit,*) nb_couches,icb,ocb,oce
  do i=1,nb_couches
     write(unit,105)r(i),rho(i),vpv(i),vsv(i),qkappa(i),qshear(i),vph(i) &
                        ,vsh(i),eta(i)
  enddo
  close(unit)
 105  format(f8.0, 3f9.2, 2f9.1, 2f9.2, f9.5)
!---------------------------------------------------------------------
end program generate_1Dmodel
!---------------------------------------------------------------------
