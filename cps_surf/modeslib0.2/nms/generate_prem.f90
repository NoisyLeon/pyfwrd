!---------------------------------------------------------------------
program generate_prem
!---------------------------------------------------------------------
  implicit none
  character(len=20) :: filename
  integer :: nb_couches
  doubleprecision, dimension(:), allocatable :: r,rho,vpv,vph,vsv,vsh &
                                               ,eta,qkappa,qshear
  integer :: choix,nbinter,i,j,unit,icb,ocb,oce,nz,iz,ii,ideb,choix_derivees
  doubleprecision :: rmax,bid,ro,vp,vs,x0,eps,ra,vpht,vsht,Qmu,eta_a  &
                    ,rdeb,rop,vpp,vsp,ropp,vppp,vspp,Qkap,r1test
  doubleprecision, dimension(:), allocatable :: rinter,dr
  integer, dimension(:), allocatable :: iinter
  logical, parameter:: gravity=.false.

  logical  :: ocean=.false.,derivees,first_in_premF=.true.,cancel_att
!
  integer :: nbpi
  doubleprecision, dimension(:), allocatable :: xint,wint,a1p,a2p,a3p,a4p,a5p,qmup,rhop
  doubleprecision, dimension(:,:), allocatable :: bord
  unit=111
  eps=1.d-12
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! input
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  print*,'Entrer le nom du model'
  print*,'Enter the model output name '
  read*,filename
!  print*,'Voulez vous PREM      (1) (le vrai) '
  print*,'Do you want PREM      (1) (the real one) '
  print*,'PREM light            (2)'
  print*,'PREM very light       (3)'
  print*,'PREM super very light (4)'
  print*,'Modele a 4 layers     (5)'
  print*,'PREM avec aniso       (6)'
  print*,'PREM light 2          (7)'
  print*,'iapei                 (8)'
  print*,'PREM D" slow          (9)'
  print*,'PREM D" high density (10)'
  print*,'PREMc                (12)'
  print*,'PREM test            (13)'
  
  read*,choix
  if (choix==1) then
     print*,'Do you need ocean (1) or not (2):'
     read*,i
     if (i==1) then 
        ocean=.true.
     else
        ocean=.false.
     endif
  endif
  print*,'Enter the number of layers of the model '
  read*,nb_couches
  print*,'Do you want to cancel attenuation (yes=T,no=F)?'
  read*,cancel_att
!  print*,'Enter the starting radius'
!  read*,rdeb
!  print*,'Enter the maximum radius (0=max of prem==6371!)'
!  read*,rmax
rdeb=0.;rmax=0.
!  print*,'Do you need derivatives at the boundaries (for SEMM, yes)? (1) yes, (0) no.'
!  read*,choix_derivees
choix_derivees=0
  select case(choix_derivees)
  case(1)
     derivees=.true.
  case(0)
     derivees=.false.
  case default
     print*,'Choix pas prevu, je mets les derivees ... &
           &(a enlever a la main si necessaire)'
      derivees=.true.
  end select
!  if (choix==13) read*,r1test
     
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! fin input
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  select case(choix)
  case(1,11,112,113)
     call prem(0.5d0,bid,bid,bid,bid,gravity,nbi=nbinter)
  case(2)
     call prem_light(0.5d0,bid,bid,bid,bid,gravity,nbi=nbinter)
  case(7)
     call prem_light2(0.5d0,bid,bid,bid,gravity,nbi=nbinter)
  case(3)
     call prem_vlight(0.5d0,bid,bid,bid,gravity,nbi=nbinter)
  case(4)
     call prem_svlight(0.5d0,bid,bid,bid,gravity,nbi=nbinter)
  case(5)
     call modele_4couches(0.5d0,bid,bid,bid,gravity,nbi=nbinter)
  case(6)
     call prem_aniso(0.5d0,bid,bid,bid,bid,bid,bid,bid,nbi=nbinter)
  case(8)
     call emiask(0.5d0,bid,bid,bid,bid,bid,gravity,nbi=nbinter)
  case(9)
     call prem_d2low(0.5d0,bid,bid,bid,bid,gravity,nbi=nbinter)
  case(10)
     call prem_d2low2(0.5d0,bid,bid,bid,bid,gravity,nbi=nbinter)
  case(12)
     call premc(0.5d0,bid,bid,bid,bid,gravity,nbi=nbinter)
  case(13)
     call prem_test(0.5d0,bid,bid,bid,bid,gravity,nbi=nbinter)
  case default
!     stop 'Cas non prevu pour le type de modele!'
     stop 'Case not written yet'
  end select
  allocate(rinter(nbinter),iinter(nbinter),dr(nbinter))
  select case(choix)
  case(1,11,112,113)
     call prem(0.5d0,bid,bid,bid,bid,gravity,ri=rinter)
  case(2)
     call prem_light(0.5d0,bid,bid,bid,bid,gravity,ri=rinter)
  case(7)
     call prem_light2(0.5d0,bid,bid,bid,gravity,ri=rinter)
  case(3)
     call prem_vlight(0.5d0,bid,bid,bid,gravity,ri=rinter)
  case(4)
     call prem_svlight(0.5d0,bid,bid,bid,gravity,ri=rinter)
  case(5)
     call modele_4couches(0.5d0,bid,bid,bid,gravity,ri=rinter)
  case(6)
     call prem_aniso(0.5d0,bid,bid,bid,bid,bid,bid,bid,ri=rinter)
  case(8)
     call emiask(0.5d0,bid,bid,bid,bid,bid,gravity,ri=rinter)
  case(9)
     call prem_d2low(0.5d0,bid,bid,bid,bid,gravity,ri=rinter)
  case(10)
     call prem_d2low2(0.5d0,bid,bid,bid,bid,gravity,ri=rinter)
  case(12)
     call premc(0.5d0,bid,bid,bid,bid,gravity,ri=rinter)
  case(13)
     call prem_test(0.5d0,bid,bid,bid,bid,gravity,ri=rinter)
  end select

  ra=rinter(nbinter)
  if (rmax==0.d0) then
     rmax=rinter(nbinter)
  else
     j=0
     do i=1,nbinter
        if (rinter(i)<rmax) j=j+1
     enddo
     nbinter=j+1
     rinter(nbinter)=rmax
     do i=1,nbinter
        print*,i,rinter(i)
     enddo
  endif  
  if (rdeb<1.E-8) then
     ideb=1 
  else
     j=0
     do i=1,nbinter
        if (rinter(i)<=rdeb) j=j+1
     enddo
     ideb=j
  endif      
  print*,'ideb=',ideb
  if (ideb/=1) rinter(ideb)=rdeb
!
  do i=ideb,nbinter
     if (i==ideb) then
        iinter(ideb)=0
!test
!else if (i==4) then
!   iinter(i)=100+iinter(i-1)
     else

        iinter(i)=nint((rinter(i)-rinter(i-1))/(rmax-rdeb)*nb_couches)+iinter(i-1)
        if (iinter(i)<=iinter(i-1)+1)iinter(i)=iinter(i-1)+2
     endif
     if ( i > ideb ) &
     dr    (i)=(rinter(i)-rinter(i-1))/(iinter(i)-(iinter(i-1)+1))
  enddo
  if ( iinter(nbinter)/=nb_couches) then
     nb_couches=iinter(nbinter)
     print*,'The final number of layers is :',     nb_couches
  endif
  allocate(r(0:nb_couches),rho(nb_couches),vpv(nb_couches),vph(nb_couches) &
       ,vsv(nb_couches),vsh(nb_couches),eta(nb_couches),qkappa(nb_couches) &
       ,qshear(nb_couches))
  r=0.0d0;rho=0.0d0;vpv=0.0d0;vph=0.0d0;vsv=0.0d0;vsh=0.0d0
  eta=0.0d0;qkappa=0.0d0;qshear=0.0d0


  do iz=ideb,nbinter
     if (iz==ideb) then
        r(0)=rinter(iz)
     else
        do ii=iinter(iz-1)+1,iinter(iz)
           if (ii==iinter(iz-1)+1) then
               r(ii)=r(ii-1)
           else
              r(ii)=r(ii-1)+dr(iz)
           endif
           x0=r(ii)/ra
           if (ii==iinter(iz-1)+1) then
              x0=x0+eps
           else if (ii==iinter(iz)) then
              x0=x0-eps
           endif
           select case(choix)
           case(1)           
              call prem(x0,ro,vp,vs,Qmu,gravity)
           case(12)           
              call premc(x0,ro,vp,vs,Qmu,gravity)
           case(13)           
              call prem_test(x0,ro,vp,vs,Qmu,gravity)
           case(2)
              call prem_light(x0,ro,vp,vs,Qmu,gravity)
           case(7)
              call prem_light2(x0,ro,vp,vs,gravity)
           case(3)
              call prem_vlight(x0,ro,vp,vs,gravity)
           case(4)
              call prem_svlight(x0,ro,vp,vs,gravity)
           case(5)
              call modele_4couches(x0,ro,vp,vs,gravity)
           case(6)
              call prem_aniso(x0,ro,vp,vpht,vs,vsht,eta_a,Qmu)
              Qkap=1.e5
!Qmu=1.e5
!!$vpht=vp
!!$vsht=vs
!!$eta_a=1.d0
           case(8)
              call emiask(x0,ro,vp,vs,Qkap,Qmu,gravity)
              vpht=vp
              vsht=vs
              eta_a=1.D0
           case(9)           
              call prem_d2low(x0,ro,vp,vs,Qmu,gravity)
           case(10)           
              call prem_d2low2(x0,ro,vp,vs,Qmu,gravity)
           end select
           if (choix/=6.and.choix/=8.and.choix/=11.and.choix/=113) then
              vpht=vp
              vsht=vs
              eta_a=1.D0
              Qkap =1.E5
           endif
           if (choix==11.or.choix==112.or.choix==113) Qkap=1.e5
           if (choix/=6.and.choix/=8.and.choix/=1.and.choix/=9 &
               .and.choix/=10.and.choix/=2.and.choix/=11.and.  &
               choix/=12.and.choix/=112.and.choix/=113)   Qmu  =1.E5
           rho(ii)   =ro  *1000.d0
           vpv(ii)   =vp  *1000.d0
           vph(ii)   =vpht*1000.d0
           vsv(ii)   =vs  *1000.d0
           vsh(ii)   =vsht*1000.d0
           eta(ii)   =eta_a
           qkappa(ii)=Qkap
           qshear(ii)=Qmu
        enddo
     endif
     !           qshear(ii)=1.E5
  enddo
!
  nz=0
  icb=0
  ocb=0
  oce=0
  do i=nb_couches-1,1,-1
     if (abs(vsv(i+1))<1.d-10.and.abs(vsv(i))>1.d-10) then
        if (ocb/=0) then
           icb=i
        else
           oce=nb_couches-i
        endif
        nz=nz+1
     else if (abs(vsv(i+1))>1.d-10.and.abs(vsv(i))<1.d-10) then
        ocb=i
        nz=nz+1
     endif
  enddo
     if (nz>3) then
!        stop'Il y a trop de couches liquide'
     else if (nz==0) then
        if (abs(vsv(1))<1.d-10) then
!tout liquide
           ocb=nb_couches
           icb=0
           oce=0
        else if (abs(vsv(1))<1.d-10) then
!tout solide
           ocb=0
           icb=0
           oce=0
        endif
     else if (nz==1) then
        icb=nb_couches-oce
        ocb=nb_couches
        oce=0
     endif
!
     if (derivees) &
     call get_derive(choix,rdeb,rmax,ro,vp,vs,rop,vpp,vsp,ropp,vppp,vspp)
     if (cancel_att) then
        qkappa(:)=1e5
        qshear(:)=1e5
     endif
  r(:)=r(:)*1000.  
  open(unit,file=filename)
  write(unit,*) filename
  write(unit,*) 1,1.,1
  write(unit,*) nb_couches,icb,ocb,oce
  if (derivees) then 
     write(unit,*) ro  ,vp  ,vs
     write(unit,*) rop ,vpp ,vsp
     write(unit,*) ropp,vppp,vspp
  endif
  do i=1,nb_couches
     write(unit,105)r(i),rho(i),vpv(i),vsv(i),qkappa(i),qshear(i),vph(i) &
                        ,vsh(i),eta(i)
  enddo
  close(unit)
 105  format(f8.0, 3f9.2, 2f9.1, 2f9.2, f9.5)

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine get_derive(choix,rdeb,rmax,ro,vp,vs,rop,vpp,vsp,ropp,vppp,vspp)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  implicit none

  integer, intent(in) :: choix
  doubleprecision, intent(in) :: rdeb,rmax
  doubleprecision, intent(out):: ro,vp,vs,rop,vpp,vsp,ropp,vppp,vspp
!
  doubleprecision :: rdtn,eps=1d-5,x0,h,bid
  doubleprecision,dimension(3) :: ro_t,vp_t,vs_t,rd
  logical :: coq
  integer :: i
!
  h=eps*rmax
  if (rdeb > 1.d-8) then
     rdtn=rdeb
     coq=.true.
     h=-h
  else
     rdtn=rmax
     coq=.false.
  endif
  rd(1)=rdtn-eps*sign(h,1.d0)
  rd(2)=rdtn-h-eps*sign(h,1.d0)
  rd(3)=rdtn-2*h-eps*sign(h,1.d0)
  do i=1,3
     x0=rd(i)/ra
     select case(choix)        
     case(1)           
        call prem(x0,ro_t(i),vp_t(i),vs_t(i),bid,gravity)
     case(2)
        call prem_light(x0,ro_t(i),vp_t(i),vs_t(i),bid,gravity)
     case(7)
        call prem_light2(x0,ro_t(i),vp_t(i),vs_t(i),gravity)
     case(3)
        call prem_vlight(x0,ro_t(i),vp_t(i),vs_t(i),gravity)
     case(4)
        call prem_svlight(x0,ro_t(i),vp_t(i),vs_t(i),gravity)
     case(5)
        call modele_4couches(x0,ro_t(i),vp_t(i),vs_t(i),gravity)
     case(6)
        stop 'Pas prevu ....'
     end select
  enddo
  ro  = ro_t(1)                        *1000.d0
  rop =(ro_t(1)-  ro_t(2)        )/h
  ropp=(ro_t(1)-2*ro_t(2)+ro_t(3))/h**2/1000.d0
!
  vp  = vp_t(1)                        *1000.d0 
  vpp =(vp_t(1)-  vp_t(2)        )/h
  vppp=(vp_t(1)-2*vp_t(2)+vp_t(3))/h**2/1000.d0
!
  vs  = vs_t(1)                        *1000.d0 
  vsp =(vs_t(1)-  vs_t(2)        )/h
  vspp=(vs_t(1)-2*vs_t(2)+vs_t(3))/h**2/1000.d0

!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
end subroutine get_derive
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!!$  subroutine prem(x0,ro,vp,vs,Qmu,flag_g,GR,nbi,ri)
!!$!---------------------------------------
!!$!
!!$!   Given non-dimensionalized radius x0, prem returns
!!$!   density, ro, compressional velocity, vp, and shear 
!!$!   velocity, vs according to the PREM model (1s ref period)
!!$!     (Anderson and Dziewonski, 1981).  
!!$!
!!$!   Also returns gravity gr
!!$!
!!$!
!!$      integer 				:: i
!!$      doubleprecision    		:: x0,x,ro,vp,vs,Qmu
!!$      logical				:: flag_g
!!$      doubleprecision,optional   	:: GR
!!$      doubleprecision,parameter   	:: Rt = 6371000.d0,epsilon = 1.d-9
!!$      doubleprecision,parameter		:: bigg = 6.6723d-11
!!$      doubleprecision,dimension(14)	:: r,q
!!$      doubleprecision,dimension(13,4)	:: d,p,s
!!$      doubleprecision,dimension(13)	:: cumul
!!$      logical				:: pastrouve
!!$      integer, optional :: nbi
!!$      doubleprecision,optional,dimension(14) :: ri
!!$!
!!$      doubleprecision			:: CST,r1,r2,t1,t2,t3,t4
!!$      integer				:: j
!!$!
!!$   ro = 0.0d0; vp = 0.d0; vs = 0.d0;
!!$
!!$!
!!$   r(1) = 0.d0;     r(2) = 1221.5d0; r(3) = 3480.d0; r(4) = 3630.d0
!!$   r(5) = 5600.d0;  r(6) = 5701.d0;  r(7) = 5771.d0; 
!!$   r(8) = 5971.d0
!!$   r(9) = 6151.d0; r(10) = 6331.d0; r(11) = 6346.6d0
!!$  r(12) = 6363.d0; r(13) = 6368.d0; r(14) = 6371.d0
!!   r(8) = 5971.d0
!   r(9) = 6280.d0; r(10) = 6285.d0; r(11) = 6290d0
!   r(12) = 6295.d0; r(13) = 6301.d0; r(14) = 6371.d0
!!   r(9) = 6350.d0; r(10) = 6355.d0; r(11) = 6360d0
!!   r(12) = 6365.d0; r(13) = 6368.d0; r(14) = 6371.d0
!!$!
!!$  q(1)=84.6d0; q(2)=1.d5; q(3)=312.d0; q(4)=312.d0;q(5)=312.d0;
!!$  q(6)=312.d0; q(7)=143.d0; q(8)=143.d0;q(9)=80.d0;q(10)=80.d0;
!!$  q(11)=600.d0; q(12)=600.d0;q(13)=600.d0
!!$  if (present(nbi)) then
!!$     nbi=14
!!$  endif
!!$  if (present(ri)) then
!!$     ri(:)=r(:)
!!$  endif
!!$!
!!$  d(:,:) = 0.d0
!!$!
!!$  d(1,1) = 5.d0
!!$  d(2,1) =  5.d0
!!$  d(3,1) =  5.d0
!!$  d(4,1) =  5.d0
!!$  d(5,1) =  5.d0
!!$  d(6,1) =  5.d0
!!$  d(7,1) =  5.d0
!!$  d(8,1) =  5.d0
!!$  d(9,1) =  5.d0
!!$  d(10,1) = 5.d0
!!$  d(11,1) = 5.d0
!!$  d(12,1) = 5.d0
!!$  d(13,1) = 5.d0 
!!$!
!!$!----
!!$!
!!$  p(:,:) = 0.d0
!!$!
!!$  p(1,1) = 8.d0 
!!$  p(2,1) = 8.d0 
!!$  p(3,1) =8.d0 
!!$  p(4,1) =8.d0 
!!$  p(5,1) =8.d0 
!!$  p(6,1) =8.d0 
!!$  p(7,1) =8.d0 
!!$  p(8,1) =8.d0 
!!$  p(9,1) =8.d0 
!!$  p(10,1)= 8.d0 
!!$!  p(11,1)= 7.5d0 
!!$!  p(12,1)= 9.d0 
!!$!  p(13,1)=9.5d0 
!!$  p(11,1)= 8.d0
!!$  p(12,1)= 8.d0 
!!$  p(13,1)=8.d0 
!!$!
!!$!----
!!$!
!!$  s(:,:) = 0.d0
!!$!
!!$  s(1,1) = 6.d0
!!$  s(2,1) = 0.d0
!!$  s(3,1) = 6.d0
!!$  s(4,1) = 6.d0
!!$  s(5,1) = 6.d0
!!$  s(6,1) = 6.d0
!!$  s(7,1) = 6.d0
!!$!  s(8,1) = 7.d0
!!$  s(8,1) = 6.d0
!!$  s(9,1) = 6.d0
!!$  s(10,1) =6.d0
!!$!  s(11,1) =5.5d0
!!$!  s(12,1) =7.d0 
!!$!  s(13,1) =7.5d0
!!$  s(11,1) =6.d0
!!$  s(12,1) =6.d0 
!!$  s(13,1) =6.d0
!!$!
!!$!
!!$  r(:) = r(:) * 1000.d0
!!$  x    = x0   * Rt
!!$!
!!$  pastrouve = .true.
!!$!
!!$ do i = 1,13
!!$
!!$ if ( pastrouve.and.(x >= r(14)) ) then
!!$  x0 = 1.d0
!!$  ro = d(13,1) + x0*( d(13,2) + x0*( d(13,3) + x0*( d(13,4) )))
!!$  vp = p(13,1) + x0*( p(13,2) + x0*( p(13,3) + x0*( p(13,4) )))
!!$  vs = s(13,1) + x0*( s(13,2) + x0*( s(13,3) + x0*( s(13,4) )))
!!$  Qmu=q(13)
!!$  pastrouve = .false.
!!$ endif
!!$
!!$ if ( pastrouve.and.((x >= r(i)).and.(x < r(i+1))) ) then
!!$  ro = d(i,1) + x0*( d(i,2) + x0*( d(i,3) + x0*( d(i,4) )))
!!$  vp = p(i,1) + x0*( p(i,2) + x0*( p(i,3) + x0*( p(i,4) )))
!!$  vs = s(i,1) + x0*( s(i,2) + x0*( s(i,3) + x0*( s(i,4) )))
!!$  Qmu=q(i)
!!$  pastrouve = .false.
!!$ endif
!!$
!!$ enddo
!!$ vs=vs*(1.+0.1*cos(2.*3.14159*((x0-1.d0)*6371/30.+6371/15.)))
!!$
!!$  if (.not. flag_g) return
!!$
!!$!-- Gravity
!!$
!!$	CST = 16.d0*datan(1.d0)*bigg
!!$
!!$	do i = 1,13
!!$		t1 = d(i,1)/3.d0
!!$		t2 = d(i,2)/(Rt*4.d0)
!!$		t3 = d(i,3)/((Rt**2)*5.d0)		
!!$		t4 = d(i,4)/((Rt**3)*6.d0)
!!$		r2 = r(i+1)
!!$		r1 = r(i)
!!$
!!$		cumul(i) =  t1*(r2**3) + t2*(r2**4) + t3*(r2**5) + t4*(r2**6) - &
!!$			  ( t1*(r1**3) + t2*(r1**4) + t3*(r1**5) + t4*(r1**6) )
!!$	
!!$	enddo
!!$
!!$	if ( x > r(14) ) x = r(14)
!!$
!!$	do i =1,13
!!$ 	if ((x > r(i)).and.(x <= r(i+1))) then
!!$
!!$		GR = 0.d0
!!$		do j = 1,i-1
!!$			GR = GR + cumul(j)
!!$		enddo
!!$
!!$		t1 = d(i,1)/3.d0
!!$		t2 = d(i,2)/(Rt*4.d0)
!!$		t3 = d(i,3)/((Rt**2)*5.d0)		
!!$		t4 = d(i,4)/((Rt**3)*6.d0)
!!$		r2 = x
!!$		r1 = r(i)
!!$
!!$		GR = GR +   t1*(r2**3) + t2*(r2**4) + t3*(r2**5) + t4*(r2**6)	&
!!$			- ( t1*(r1**3) + t2*(r1**4) + t3*(r1**5) + t4*(r1**6) )
!!$		GR = GR * CST
!!$		if (x > epsilon) GR = GR/(x**2)
!!$		return
!!$	endif
!!$	enddo
!!$!
!!$  end subroutine prem
!!$
  subroutine prem(x0,ro,vp,vs,Qmu,flag_g,GR,nbi,ri)
!---------------------------------------
!
!   Given non-dimensionalized radius x0, prem returns
!   density, ro, compressional velocity, vp, and shear 
!   velocity, vs according to the PREM model (1s ref period)
!     (Anderson and Dziewonski, 1981).  
!
!   Also returns gravity gr
!
!
      integer 				:: i
      doubleprecision    		:: x0,x,ro,vp,vs,Qmu
      logical				:: flag_g
      doubleprecision,optional   	:: GR
      doubleprecision,parameter   	:: Rt = 6371000.d0,epsilon = 1.d-9
      doubleprecision,parameter		:: bigg = 6.6723d-11
      doubleprecision,dimension(14)	:: r,q
      doubleprecision,dimension(13,4)	:: d,p,s
      doubleprecision,dimension(13)	:: cumul
      logical				:: pastrouve
      integer, optional :: nbi
      doubleprecision,optional,dimension(14) :: ri
!
      doubleprecision			:: CST,r1,r2,t1,t2,t3,t4
      integer				:: j
!
   ro = 0.0d0; vp = 0.d0; vs = 0.d0;

!
   r(1) = 0.d0;     r(2) = 1221.5d0; r(3) = 3480.d0; r(4) = 3630.d0
   r(5) = 5600.d0;  r(6) = 5701.d0;  r(7) = 5771.d0; 
   r(8) = 5971.d0
   r(9) = 6151.d0; r(10) = 6291.d0; r(11) = 6346.6d0
  r(12) = 6356.d0; r(13) = 6368.d0; r(14) = 6371.d0
!   r(11) = 6360.d0
!  r(12) = 6367.d0; r(13) = 6369.d0; r(14) = 6371.d0
!prem
!  q(1)=84.6d0; q(2)=1.d5; q(3)=312.d0; q(4)=312.d0;q(5)=312.d0;
!  q(6)=312.d0; q(7)=143.d0; q(8)=143.d0;q(9)=80.d0;q(10)=80.d0;
!  q(11)=600.d0; q(12)=600.d0;q(13)=600.d0
!QL6:
  q(1)=84.6d0; q(2)=1.d5; q(3)=355.d0; q(4)=355.d0;q(5)=355.d0;
  q(6)=165.d0; q(7)=165.d0; q(8)=165.d0;q(9)=70.d0;q(10)=191.d0;
  q(11)=300.d0; q(12)=300.d0;q(13)=300.d0
  if (present(nbi)) then
     nbi=14
  endif
  if (present(ri)) then
     ri(:)=r(:)
  endif
!
  d(:,:) = 0.d0
!
  d(1,1) = 13.0885d0; 			  d(1,3) = -8.8381d0 
  d(2,1) = 12.5815d0; d(2,2) = -1.2638d0; d(2,3) = -3.6426d0; d(2,4) = -5.5281d0
  d(3,1) = 7.9565d0 ; d(3,2) = -6.4761;   d(3,3) = 5.5283d0;  d(3,4) = -3.0807d0
  d(4,1) = 7.9565d0 ; d(4,2) = -6.4761;   d(4,3) = 5.5283d0;  d(4,4) = -3.0807d0
  d(5,1) = 7.9565d0 ; d(5,2) = -6.4761;   d(5,3) = 5.5283d0;  d(5,4) = -3.0807d0
  d(6,1) = 5.3197d0 ; d(6,2) = -1.4836d0
  d(7,1) = 11.2494d0; d(7,2) = -8.0298d0
  d(8,1) = 7.1089d0 ; d(8,2) = -3.8045d0
  d(9,1) = 2.6910d0 ; d(9,2) = 0.6924d0
  d(10,1) = 2.6910d0; d(10,2) = 0.6924d0
  d(11,1) = 2.9d0  
  d(12,1) = 2.6d0  
  if (ocean) then
!
! oceanique
     d(13,1) = 1.02d0 
  else
!
! continental
     d(13,1) = d(12,1)
  endif
!!$!d(10,1) = 2.9d0
!!$!d(11,1) = 2.9d0
!!$!d(12,1) = 2.9d0
!!$!d(13,1) = 2.9d0
!
!----
!
  p(:,:) = 0.d0
!
  p(1,1) = 11.2622d0 ; p(1,3) = -6.3640d0
  p(2,1) = 11.0487d0 ; p(2,2) = -4.0362d0; p(2,3)  = 4.8023d0; p(2,4) = -13.5732d0
  p(3,1) = 15.3891d0 ; p(3,2) = -5.3181d0; p(3,3)  = 5.5242d0; p(3,4) = -2.5514d0
  p(4,1) = 24.952d0 ; p(4,2)  = -40.4673d0; p(4,3) = 51.4832d0; p(4,4) = -26.6419d0
  p(5,1) = 29.2766d0 ; p(5,2) = -23.6027d0; p(5,3) = 5.5242d0; p(5,4) = -2.5514d0
  p(6,1) = 19.0957d0 ; p(6,2)  = -9.8672d0
  p(7,1) = 39.7027d0 ; p(7,2)  = -32.6166d0
  p(8,1) = 20.3926d0 ; p(8,2)  = -12.2569d0
  p(9,1) = 4.1875d0 ; p(9,2)  = 3.9382d0
  p(10,1) = 4.1875d0 ; p(10,2) = 3.9382d0
  p(11,1) = 6.8d0 
  p(12,1) = 5.8d0
  if (ocean) then
!
! oceanique
     p(13,1) = 1.45d0 
  else
!
! continental
     p(13,1) = p(12,1)
  endif
!!$!  p(10,1) = 5.8d0
!!$!  p(11,1) = 5.8d0 
!!$!  p(12,1) = 5.8d0
!!$!  p(13,1) = 5.8d0 
!
!----
!
  s(:,:) = 0.d0
!
  s(1,1) = 3.6678d0; s(1,3) = -4.4475d0

  s(3,1) = 6.9254d0; s(3,2) = 1.4672d0; s(3,3) = -2.0834d0; s(3,4) = 0.9783d0
  s(4,1) = 11.1671d0; s(4,2) = -13.7818d0; s(4,3) = 17.4575d0; s(4,4) = -9.2777d0
  s(5,1) = 22.3459d0; s(5,2) = -17.2473d0; s(5,3) = -2.0834d0; s(5,4) = 0.9783d0
  s(6,1) = 9.9839; s(6,2) = -4.9324
  s(7,1) = 22.3512d0; s(7,2) = -18.5856d0 
  s(8,1) = 8.9496d0; s(8,2) = -4.4597
  s(9,1) = 2.1519d0; s(9,2) = 2.3481d0
  s(10,1) = 2.1519d0; s(10,2) = 2.3481d0
  s(11,1) = 3.9d0 
  s(12,1) = 3.2d0 
! oceanique (ne rien toucher)

!
! continental
  if (.not.ocean)  s(13,1) = s(12,1)
!!$
!!$!  s(11,1) = 3.0d0 
!!$!  s(12,1) = 4.8d0 
!!$!  s(13,1) = 3.9d0
!  s(10,1) = 2.1519d0; s(10,2) = 2.3481d0
!  s(11,1) = 2.1519d0; s(11,2) = 2.3481d0
!  s(12,1) = 2.1519d0; s(12,2) = 2.3481d0
!  s(13,1) = 2.1519d0; s(13,2) = 2.3481d0
!!$!  s(10,1) = 2.1519d0
!!$!  s(11,1) = 2.1519d0
!!$!  s(12,1) = 2.1519d0
!!$!  s(13,1) = 2.1519d0
!!$!
!
  r(:) = r(:) * 1000.d0
  x    = x0   * Rt
!
  pastrouve = .true.
!
 do i = 1,13

 if ( pastrouve.and.(x >= r(14)) ) then
  x0 = 1.d0
  ro = d(13,1) + x0*( d(13,2) + x0*( d(13,3) + x0*( d(13,4) )))
  vp = p(13,1) + x0*( p(13,2) + x0*( p(13,3) + x0*( p(13,4) )))
  vs = s(13,1) + x0*( s(13,2) + x0*( s(13,3) + x0*( s(13,4) )))
  Qmu=q(13)
  pastrouve = .false.
 endif

 if ( pastrouve.and.((x >= r(i)).and.(x < r(i+1))) ) then
  ro = d(i,1) + x0*( d(i,2) + x0*( d(i,3) + x0*( d(i,4) )))
  vp = p(i,1) + x0*( p(i,2) + x0*( p(i,3) + x0*( p(i,4) )))
  vs = s(i,1) + x0*( s(i,2) + x0*( s(i,3) + x0*( s(i,4) )))
  Qmu=q(i)
  pastrouve = .false.
 endif

 enddo

  if (.not. flag_g) return

!-- Gravity

	CST = 16.d0*datan(1.d0)*bigg

	do i = 1,13
		t1 = d(i,1)/3.d0
		t2 = d(i,2)/(Rt*4.d0)
		t3 = d(i,3)/((Rt**2)*5.d0)		
		t4 = d(i,4)/((Rt**3)*6.d0)
		r2 = r(i+1)
		r1 = r(i)

		cumul(i) =  t1*(r2**3) + t2*(r2**4) + t3*(r2**5) + t4*(r2**6) - &
			  ( t1*(r1**3) + t2*(r1**4) + t3*(r1**5) + t4*(r1**6) )
	
	enddo

	if ( x > r(14) ) x = r(14)

	do i =1,13
 	if ((x > r(i)).and.(x <= r(i+1))) then

		GR = 0.d0
		do j = 1,i-1
			GR = GR + cumul(j)
		enddo

		t1 = d(i,1)/3.d0
		t2 = d(i,2)/(Rt*4.d0)
		t3 = d(i,3)/((Rt**2)*5.d0)		
		t4 = d(i,4)/((Rt**3)*6.d0)
		r2 = x
		r1 = r(i)

		GR = GR +   t1*(r2**3) + t2*(r2**4) + t3*(r2**5) + t4*(r2**6)	&
			- ( t1*(r1**3) + t2*(r1**4) + t3*(r1**5) + t4*(r1**6) )
		GR = GR * CST
		if (x > epsilon) GR = GR/(x**2)
		return
	endif
	enddo
!
  end subroutine prem
!

  subroutine prem_d2low(x0,ro,vp,vs,Qmu,flag_g,GR,nbi,ri)
!---------------------------------------
! genere un modele prem avec la couche D" avec une ULVZ
!   Given non-dimensionalized radius x0, prem returns
!   density, ro, compressional velocity, vp, and shear 
!   velocity, vs according to the PREM model (1s ref period)
!     (Anderson and Dziewonski, 1981).  
!
!   Also returns gravity gr
!
!
      integer 				:: i
      doubleprecision    		:: x0,x,ro,vp,vs,Qmu
      logical				:: flag_g
      doubleprecision,optional   	:: GR
      doubleprecision,parameter   	:: Rt = 6371000.d0,epsilon = 1.d-9
      doubleprecision,parameter		:: bigg = 6.6723d-11
      doubleprecision,dimension(15)	:: r,q
      doubleprecision,dimension(14,4)	:: d,p,s
      doubleprecision,dimension(14)	:: cumul
      logical				:: pastrouve
      integer, optional :: nbi
      doubleprecision,optional,dimension(15) :: ri
!
      doubleprecision			:: CST,r1,r2,t1,t2,t3,t4
      integer				:: j
!
   ro = 0.0d0; vp = 0.d0; vs = 0.d0;
!
   r(1) = 0.d0;     r(2) = 1221.5d0; r(3) = 3480.d0; r(4)=3580.d0
   r(5) = 3630.d0
   r(6) = 5600.d0;  r(7) = 5701.d0;  r(8) = 5771.d0; r(9) = 5971.d0
   r(10) = 6151.d0;  r(11) = 6291.d0; r(12) = 6346.6d0
   r(13) = 6356.d0; r(14) = 6368.d0; r(15) = 6371.d0
!
  q(1)=84.6d0; q(2)=1.d5; q(3)=312.d0; q(4)=312.d0; q(5)=312.d0;q(6)=312.d0;
  q(7)=312.d0; q(8)=143.d0; q(9)=143.d0;q(10)=80.d0;q(11)=80.d0;
  q(12)=600.d0; q(13)=600.d0;q(14)=600.d0
!
  if (present(nbi)) then
     nbi=15
  endif
  if (present(ri)) then
     ri(:)=r(:)
  endif
!
  d(:,:) = 0.d0
!
  d(1,1) = 13.0885d0; 			  d(1,3) = -8.8381d0 
  d(2,1) = 12.5815d0; d(2,2) = -1.2638d0; d(2,3) = -3.6426d0; d(2,4) = -5.5281d0
!  d(3,1) = 7.9565d0 ; d(3,2) = -6.4761;   d(3,3) = 5.5283d0;  d(3,4) = -3.0807d0
!d high plateau
!  d(3,1) = 6.0d0
!d low palteau
  d(3,1) = 4.3d0
!d high gradiant
!  d(3,1) = 40.45d0  ; d(3,2) =-62.181d0
!d low gradiant
!  d(3,1) =-28.887d0  ; d(3,2) = 61.181d0
  d(4,1) = 7.9565d0 ; d(4,2) = -6.4761;   d(4,3) = 5.5283d0;  d(4,4) = -3.0807d0
  d(5,1) = 7.9565d0 ; d(5,2) = -6.4761;   d(5,3) = 5.5283d0;  d(5,4) = -3.0807d0
  d(6,1) = 7.9565d0 ; d(6,2) = -6.4761;   d(6,3) = 5.5283d0;  d(6,4) = -3.0807d0
  d(7,1) = 5.3197d0 ; d(7,2) = -1.4836d0
  d(8,1) = 11.2494d0; d(8,2) = -8.0298d0
  d(9,1) = 7.1089d0 ; d(9,2) = -3.8045d0
  d(10,1) = 2.6910d0 ; d(10,2) = 0.6924d0
  d(11,1) = 2.6910d0; d(11,2) = 0.6924d0
  d(12,1) = 2.9d0  
  d(13,1) = 2.6d0  
  if (ocean) then
!
! oceanique
     d(14,1) = 1.02d0 
  else
!
! continental
     d(14,1) = d(13,1)
  endif
!
!----
!
  p(:,:) = 0.d0
!
  p(1,1) = 11.2622d0 ; p(1,3) = -6.3640d0
  p(2,1) = 11.0487d0 ; p(2,2) = -4.0362d0;  p(2,3)  = 4.8023d0; p(2,4) = -13.5732d0
  p(3,1) = 15.3891d0 ; p(3,2) = -5.3181d0; p(3,3)  = 5.5242d0; p(3,4) = -2.5514d0
!  p(3,1) = -0.696d0  ; p(3,2) = 25.99d0 !vp faible pres de la cmb
!  p(3,1) = 20.69d0   ; p(3,2) = -12.61d0 !vp plus forte pres de la cmb
  p(4,1) = 15.3891d0 ; p(4,2) = -5.3181d0;  p(4,3)  = 5.5242d0; p(4,4) = -2.5514d0
  p(5,1) = 24.952d0  ; p(5,2)  = -40.4673d0;p(5,3) = 51.4832d0; p(5,4) = -26.6419d0
  p(6,1) = 29.2766d0 ; p(6,2) = -23.6027d0; p(6,3) = 5.5242d0;  p(6,4) = -2.5514d0
  p(7,1) = 19.0957d0 ; p(7,2)  = -9.8672d0
  p(8,1) = 39.7027d0 ; p(8,2)  = -32.6166d0
  p(9,1) = 20.3926d0 ; p(9,2)  = -12.2569d0
  p(10,1) = 4.1875d0 ; p(10,2)  = 3.9382d0
  p(11,1) = 4.1875d0 ; p(11,2) = 3.9382d0
  p(12,1) = 6.8d0 
  p(13,1) = 5.8d0
  if (ocean) then
!
! oceanique
     p(14,1) = 1.45d0 
  else
!
! continental
     p(14,1) = p(13,1)
  endif
!
!----
!
  s(:,:) = 0.d0
!
  s(1,1) = 3.6678d0; s(1,3) = -4.4475d0
! s(2.:) est nul et le reste
!prem
  s(3,1) = 6.9254d0; s(3,2) = 1.4672d0; s(3,3) = -2.0834d0; s(3,4) = 0.9783d0
!vs faible:
!  s(3,1) = -2.22d0; s(3,2) = 16.88d0
  s(4,1) = 6.9254d0; s(4,2) = 1.4672d0; s(4,3) = -2.0834d0; s(4,4) = 0.9783d0
  s(5,1) = 11.1671d0; s(5,2) = -13.7818d0; s(5,3) = 17.4575d0; s(5,4) = -9.2777d0
  s(6,1) = 22.3459d0; s(6,2) = -17.2473d0; s(6,3) = -2.0834d0; s(6,4) = 0.9783d0
  s(7,1) = 9.9839; s(7,2) = -4.9324
  s(8,1) = 22.3512d0; s(8,2) = -18.5856d0 
  s(9,1) = 8.9496d0; s(9,2) = -4.4597
  s(10,1) = 2.1519d0; s(10,2) = 2.3481d0
  s(11,1) = 2.1519d0; s(11,2) = 2.3481d0
  s(12,1) = 3.9d0 
  s(13,1) = 3.2d0 
!
! oceanique (ne rien toucher)

!
! continental
  if (.not.ocean)  s(14,1) = s(13,1)
!
!
  r(:) = r(:) * 1000.d0
  x    = x0   * Rt
!
  pastrouve = .true.
!
 do i = 1,14

 if ( pastrouve.and.(x >= r(15)) ) then
  x0 = 1.d0
  ro = d(13,1) + x0*( d(14,2) + x0*( d(14,3) + x0*( d(14,4) )))
  vp = p(13,1) + x0*( p(14,2) + x0*( p(14,3) + x0*( p(14,4) )))
  vs = s(13,1) + x0*( s(14,2) + x0*( s(14,3) + x0*( s(14,4) )))
  Qmu=q(13)
  pastrouve = .false.
 endif

 if ( pastrouve.and.((x >= r(i)).and.(x < r(i+1))) ) then
  ro = d(i,1) + x0*( d(i,2) + x0*( d(i,3) + x0*( d(i,4) )))
  vp = p(i,1) + x0*( p(i,2) + x0*( p(i,3) + x0*( p(i,4) )))
  vs = s(i,1) + x0*( s(i,2) + x0*( s(i,3) + x0*( s(i,4) )))
  Qmu=q(i)
  pastrouve = .false.
 endif

 enddo

  if (.not. flag_g) return

!-- Gravity

	CST = 16.d0*datan(1.d0)*bigg

	do i = 1,13
		t1 = d(i,1)/3.d0
		t2 = d(i,2)/(Rt*4.d0)
		t3 = d(i,3)/((Rt**2)*5.d0)		
		t4 = d(i,4)/((Rt**3)*6.d0)
		r2 = r(i+1)
		r1 = r(i)

		cumul(i) =  t1*(r2**3) + t2*(r2**4) + t3*(r2**5) + t4*(r2**6) - &
			  ( t1*(r1**3) + t2*(r1**4) + t3*(r1**5) + t4*(r1**6) )
	
	enddo

	if ( x > r(14) ) x = r(14)

	do i =1,13
 	if ((x > r(i)).and.(x <= r(i+1))) then

		GR = 0.d0
		do j = 1,i-1
			GR = GR + cumul(j)
		enddo

		t1 = d(i,1)/3.d0
		t2 = d(i,2)/(Rt*4.d0)
		t3 = d(i,3)/((Rt**2)*5.d0)		
		t4 = d(i,4)/((Rt**3)*6.d0)
		r2 = x
		r1 = r(i)

		GR = GR +   t1*(r2**3) + t2*(r2**4) + t3*(r2**5) + t4*(r2**6)	&
			- ( t1*(r1**3) + t2*(r1**4) + t3*(r1**5) + t4*(r1**6) )
		GR = GR * CST
		if (x > epsilon) GR = GR/(x**2)
		return
	endif
	enddo
!
  end subroutine prem_d2low
!


  subroutine prem_d2low2(x0,ro,vp,vs,Qmu,flag_g,GR,nbi,ri)
!---------------------------------------
! genere un modele prem avec la couche D" avec une ULVZ
!   Given non-dimensionalized radius x0, prem returns
!   density, ro, compressional velocity, vp, and shear 
!   velocity, vs according to the PREM model (1s ref period)
!     (Anderson and Dziewonski, 1981).  
!
!   Also returns gravity gr
!
!
      integer 				:: i
      doubleprecision    		:: x0,x,ro,vp,vs,Qmu
      logical				:: flag_g
      doubleprecision,optional   	:: GR
      doubleprecision,parameter   	:: Rt = 6371000.d0,epsilon = 1.d-9
      doubleprecision,parameter		:: bigg = 6.6723d-11
      doubleprecision,dimension(15)	:: r,q
      doubleprecision,dimension(14,4)	:: d,p,s
      doubleprecision,dimension(14)	:: cumul
      logical				:: pastrouve
      integer, optional :: nbi
      doubleprecision,optional,dimension(15) :: ri
!
      doubleprecision			:: CST,r1,r2,t1,t2,t3,t4
      integer				:: j
!
   ro = 0.0d0; vp = 0.d0; vs = 0.d0;
!
   r(1) = 0.d0;     r(2) = 1221.5d0; r(3) = 3480.d0; r(4)=3485.d0
   r(5) = 3630.d0
   r(6) = 5600.d0;  r(7) = 5701.d0;  r(8) = 5771.d0; r(9) = 5971.d0
   r(10) = 6151.d0;  r(11) = 6291.d0; r(12) = 6346.6d0
   r(13) = 6356.d0; r(14) = 6368.d0; r(15) = 6371.d0
!
  q(1)=84.6d0; q(2)=1.d5; q(3)=312.d0; q(4)=312.d0; q(5)=312.d0;q(6)=312.d0;
  q(7)=312.d0; q(8)=143.d0; q(9)=143.d0;q(10)=80.d0;q(11)=80.d0;
  q(12)=600.d0; q(13)=600.d0;q(14)=600.d0
!
  if (present(nbi)) then
     nbi=15
  endif
  if (present(ri)) then
     ri(:)=r(:)
  endif
!
  d(:,:) = 0.d0
!
  d(1,1) = 13.0885d0; 			  d(1,3) = -8.8381d0 
  d(2,1) = 12.5815d0; d(2,2) = -1.2638d0; d(2,3) = -3.6426d0; d(2,4) = -5.5281d0
!  d(3,1) = 7.9565d0 ; d(3,2) = -6.4761;   d(3,3) = 5.5283d0;  d(3,4) = -3.0807d0
!d high plateau
  d(3,1) = 7.6d0
!d low palteau
!  d(3,1) = 4.3d0
!d high gradiant
!  d(3,1) = 40.45d0  ; d(3,2) =-62.181d0
!d low gradiant
!  d(3,1) =-28.887d0  ; d(3,2) = 61.181d0
  d(4,1) = 7.9565d0 ; d(4,2) = -6.4761;   d(4,3) = 5.5283d0;  d(4,4) = -3.0807d0
  d(5,1) = 7.9565d0 ; d(5,2) = -6.4761;   d(5,3) = 5.5283d0;  d(5,4) = -3.0807d0
  d(6,1) = 7.9565d0 ; d(6,2) = -6.4761;   d(6,3) = 5.5283d0;  d(6,4) = -3.0807d0
  d(7,1) = 5.3197d0 ; d(7,2) = -1.4836d0
  d(8,1) = 11.2494d0; d(8,2) = -8.0298d0
  d(9,1) = 7.1089d0 ; d(9,2) = -3.8045d0
  d(10,1) = 2.6910d0 ; d(10,2) = 0.6924d0
  d(11,1) = 2.6910d0; d(11,2) = 0.6924d0
  d(12,1) = 2.9d0  
  d(13,1) = 2.6d0  
  if (ocean) then
!
! oceanique
     d(14,1) = 1.02d0 
  else
!
! continental
     d(14,1) = d(13,1)
  endif
!
!----
!
  p(:,:) = 0.d0
!
  p(1,1) = 11.2622d0 ; p(1,3) = -6.3640d0
  p(2,1) = 11.0487d0 ; p(2,2) = -4.0362d0;  p(2,3)  = 4.8023d0; p(2,4) = -13.5732d0
!  p(3,1) = 15.3891d0 ; p(3,2) = -5.3181d0; p(3,3)  = 5.5242d0; p(3,4) = -2.5514d0
!  p(3,1) = -0.696d0  ; p(3,2) = 25.99d0 !vp faible pres de la cmb
!  p(3,1) = 20.69d0   ; p(3,2) = -12.61d0 !vp plus forte pres de la cmb
  p(3,1) =9.75d0
  p(4,1) = 15.3891d0 ; p(4,2) = -5.3181d0;  p(4,3)  = 5.5242d0; p(4,4) = -2.5514d0
  p(5,1) = 24.952d0  ; p(5,2)  = -40.4673d0;p(5,3) = 51.4832d0; p(5,4) = -26.6419d0
  p(6,1) = 29.2766d0 ; p(6,2) = -23.6027d0; p(6,3) = 5.5242d0;  p(6,4) = -2.5514d0
  p(7,1) = 19.0957d0 ; p(7,2)  = -9.8672d0
  p(8,1) = 39.7027d0 ; p(8,2)  = -32.6166d0
  p(9,1) = 20.3926d0 ; p(9,2)  = -12.2569d0
  p(10,1) = 4.1875d0 ; p(10,2)  = 3.9382d0
  p(11,1) = 4.1875d0 ; p(11,2) = 3.9382d0
  p(12,1) = 6.8d0 
  p(13,1) = 5.8d0
  if (ocean) then
!
! oceanique
     p(14,1) = 1.45d0 
  else
!
! continental
     p(14,1) = p(13,1)
  endif
!
!----
!
  s(:,:) = 0.d0
!
  s(1,1) = 3.6678d0; s(1,3) = -4.4475d0
! s(2.:) est nul et le reste
!prem
!  s(3,1) = 6.9254d0; s(3,2) = 1.4672d0; s(3,3) = -2.0834d0; s(3,4) = 0.9783d0
!vs faible:
!  s(3,1) = -2.22d0; s(3,2) = 16.88d0
!  s(3,1) = 4.0d0
  s(3,1) = 10.0d0
  s(4,1) = 6.9254d0; s(4,2) = 1.4672d0; s(4,3) = -2.0834d0; s(4,4) = 0.9783d0
  s(5,1) = 11.1671d0; s(5,2) = -13.7818d0; s(5,3) = 17.4575d0; s(5,4) = -9.2777d0
  s(6,1) = 22.3459d0; s(6,2) = -17.2473d0; s(6,3) = -2.0834d0; s(6,4) = 0.9783d0
  s(7,1) = 9.9839; s(7,2) = -4.9324
  s(8,1) = 22.3512d0; s(8,2) = -18.5856d0 
  s(9,1) = 8.9496d0; s(9,2) = -4.4597
  s(10,1) = 2.1519d0; s(10,2) = 2.3481d0
  s(11,1) = 2.1519d0; s(11,2) = 2.3481d0
  s(12,1) = 3.9d0 
  s(13,1) = 3.2d0 
!
! oceanique (ne rien toucher)

!
! continental
  if (.not.ocean)  s(14,1) = s(13,1)
!
!
  r(:) = r(:) * 1000.d0
  x    = x0   * Rt
!
  pastrouve = .true.
!
 do i = 1,14

 if ( pastrouve.and.(x >= r(15)) ) then
  x0 = 1.d0
  ro = d(13,1) + x0*( d(14,2) + x0*( d(14,3) + x0*( d(14,4) )))
  vp = p(13,1) + x0*( p(14,2) + x0*( p(14,3) + x0*( p(14,4) )))
  vs = s(13,1) + x0*( s(14,2) + x0*( s(14,3) + x0*( s(14,4) )))
  Qmu=q(13)
  pastrouve = .false.
 endif

 if ( pastrouve.and.((x >= r(i)).and.(x < r(i+1))) ) then
  ro = d(i,1) + x0*( d(i,2) + x0*( d(i,3) + x0*( d(i,4) )))
  vp = p(i,1) + x0*( p(i,2) + x0*( p(i,3) + x0*( p(i,4) )))
  vs = s(i,1) + x0*( s(i,2) + x0*( s(i,3) + x0*( s(i,4) )))
  Qmu=q(i)
  pastrouve = .false.
 endif

 enddo

  if (.not. flag_g) return

!-- Gravity

	CST = 16.d0*datan(1.d0)*bigg

	do i = 1,13
		t1 = d(i,1)/3.d0
		t2 = d(i,2)/(Rt*4.d0)
		t3 = d(i,3)/((Rt**2)*5.d0)		
		t4 = d(i,4)/((Rt**3)*6.d0)
		r2 = r(i+1)
		r1 = r(i)

		cumul(i) =  t1*(r2**3) + t2*(r2**4) + t3*(r2**5) + t4*(r2**6) - &
			  ( t1*(r1**3) + t2*(r1**4) + t3*(r1**5) + t4*(r1**6) )
	
	enddo

	if ( x > r(14) ) x = r(14)

	do i =1,13
 	if ((x > r(i)).and.(x <= r(i+1))) then

		GR = 0.d0
		do j = 1,i-1
			GR = GR + cumul(j)
		enddo

		t1 = d(i,1)/3.d0
		t2 = d(i,2)/(Rt*4.d0)
		t3 = d(i,3)/((Rt**2)*5.d0)		
		t4 = d(i,4)/((Rt**3)*6.d0)
		r2 = x
		r1 = r(i)

		GR = GR +   t1*(r2**3) + t2*(r2**4) + t3*(r2**5) + t4*(r2**6)	&
			- ( t1*(r1**3) + t2*(r1**4) + t3*(r1**5) + t4*(r1**6) )
		GR = GR * CST
		if (x > epsilon) GR = GR/(x**2)
		return
	endif
	enddo
!
  end subroutine prem_d2low2
!



  subroutine prem_light2(x0,ro,vp,vs,flag_g,GR,nbi,ri)
!---------------------------------------
!
!   Given non-dimensionalized radius x0, prem returns
!   density, ro, compressional velocity, vp, and shear 
!   velocity, vs according to the PREM model (1s ref period)
!     (Anderson and Dziewonski, 1981).  
!
!   Also returns gravity gr
!
!
      integer 				:: i
      doubleprecision    		:: x0,x,ro,vp,vs
      logical				:: flag_g
      doubleprecision,optional   	:: GR
      doubleprecision,parameter   	:: Rt = 6371000.d0,epsilon = 1.d-9
      doubleprecision,parameter		:: bigg = 6.6723d-11
      doubleprecision,dimension(11)	:: r
      doubleprecision,dimension(10,4)	:: d,p,s
      doubleprecision,dimension(10)	:: cumul
      logical				:: pastrouve
      integer, optional :: nbi
      doubleprecision,optional,dimension(11) :: ri
!
      doubleprecision			:: CST,r1,r2,t1,t2,t3,t4
      integer				:: j
!
   ro = 0.0d0; vp = 0.d0; vs = 0.d0;

!
   r(1) = 0.d0;     r(2) = 1221.5d0; r(3) = 3480.d0; r(4) = 3630.d0
   r(5) = 5600.d0;  r(6) = 5701.d0;  r(7) = 5771.d0; r(8) = 5971.d0
   r(9) = 6151.d0; r(10) = 6291.d0; 
!r(11) = 6346.6d0
!  r(12) = 6356.d0; r(13) = 6368.d0; 
r(11) = 6371.d0
!
  if (present(nbi)) then
     nbi=11
  endif
  if (present(ri)) then
     ri(:)=r(:)
  endif
!
  d(:,:) = 0.d0
!
  d(1,1) = 13.0885d0; 			  d(1,3) = -8.8381d0 
  d(2,1) = 12.5815d0; d(2,2) = -1.2638d0; d(2,3) = -3.6426d0; d(2,4) = -5.5281d0
  d(3,1) = 7.9565d0 ; d(3,2) = -6.4761;   d(3,3) = 5.5283d0;  d(3,4) = -3.0807d0
  d(4,1) = 7.9565d0 ; d(4,2) = -6.4761;   d(4,3) = 5.5283d0;  d(4,4) = -3.0807d0
  d(5,1) = 7.9565d0 ; d(5,2) = -6.4761;   d(5,3) = 5.5283d0;  d(5,4) = -3.0807d0
  d(6,1) = 5.3197d0 ; d(6,2) = -1.4836d0
  d(7,1) = 11.2494d0; d(7,2) = -8.0298d0
  d(8,1) = 7.1089d0 ; d(8,2) = -3.8045d0
  d(9,1) = 2.6910d0 ; d(9,2) = 0.6924d0
  d(10,1) = 2.6910d0; d(10,2) = 0.6924d0
!  d(11,1) = 2.9d0  
!  d(12,1) = 2.6d0  
!  if (ocean) then
!
! oceanique
!     d(13,1) = 1.02d0 
!  else
!
! continental
!     d(13,1) = d(12,1)
!  endif
!
!----
!
  p(:,:) = 0.d0
!
  p(1,1) = 11.2622d0 ; p(1,3) = -6.3640d0
  p(2,1) = 11.0487d0 ; p(2,2) = -4.0362d0; p(2,3)  = 4.8023d0; p(2,4) = -13.5732d0
  p(3,1) = 15.3891d0 ; p(3,2) = -5.3181d0; p(3,3)  = 5.5242d0; p(3,4) = -2.5514d0
  p(4,1) = 24.952d0 ; p(4,2)  = -40.4673d0; p(4,3) = 51.4832d0; p(4,4) = -26.6419d0
  p(5,1) = 29.2766d0 ; p(5,2) = -23.6027d0; p(5,3) = 5.5242d0; p(5,4) = -2.5514d0
  p(6,1) = 19.0957d0 ; p(6,2)  = -9.8672d0
  p(7,1) = 39.7027d0 ; p(7,2)  = -32.6166d0
  p(8,1) = 20.3926d0 ; p(8,2)  = -12.2569d0
  p(9,1) = 4.1875d0 ; p(9,2)  = 3.9382d0
  p(10,1) = 4.1875d0 ; p(10,2) = 3.9382d0
!  p(11,1) = 6.8d0 
!  p(12,1) = 5.8d0
!  if (ocean) then
!
! oceanique
!     p(13,1) = 1.45d0 
!  else
!
! continental
!     p(13,1) = p(12,1)
!  endif
!
!----
!
  s(:,:) = 0.d0
!
  s(1,1) = 3.6678d0; s(1,3) = -4.4475d0

  s(3,1) = 6.9254d0; s(3,2) = 1.4672d0; s(3,3) = -2.0834d0; s(3,4) = 0.9783d0
  s(4,1) = 11.1671d0; s(4,2) = -13.7818d0; s(4,3) = 17.4575d0; s(4,4) = -9.2777d0
  s(5,1) = 22.3459d0; s(5,2) = -17.2473d0; s(5,3) = -2.0834d0; s(5,4) = 0.9783d0
  s(6,1) = 9.9839; s(6,2) = -4.9324
  s(7,1) = 22.3512d0; s(7,2) = -18.5856d0 
  s(8,1) = 8.9496d0; s(8,2) = -4.4597
  s(9,1) = 2.1519d0; s(9,2) = 2.3481d0
  s(10,1) = 2.1519d0; s(10,2) = 2.3481d0
!  s(11,1) = 3.9d0 
!  s(12,1) = 3.2d0 
!
! oceanique (ne rien toucher)

!
! continental
!  if (.not.ocean)  s(13,1) = s(12,1)
!
!
  r(:) = r(:) * 1000.d0
  x    = x0   * Rt
!
  pastrouve = .true.
!
 do i = 1,10

 if ( pastrouve.and.(x >= r(11)) ) then
  x0 = 1.d0
  ro = d(10,1) + x0*( d(10,2) + x0*( d(10,3) + x0*( d(10,4) )))
  vp = p(10,1) + x0*( p(10,2) + x0*( p(10,3) + x0*( p(10,4) )))
  vs = s(10,1) + x0*( s(10,2) + x0*( s(10,3) + x0*( s(10,4) )))
  pastrouve = .false.
 endif

 if ( pastrouve.and.((x >= r(i)).and.(x < r(i+1))) ) then
  ro = d(i,1) + x0*( d(i,2) + x0*( d(i,3) + x0*( d(i,4) )))
  vp = p(i,1) + x0*( p(i,2) + x0*( p(i,3) + x0*( p(i,4) )))
  vs = s(i,1) + x0*( s(i,2) + x0*( s(i,3) + x0*( s(i,4) )))
  pastrouve = .false.
 endif

 enddo

  if (.not. flag_g) return

!-- Gravity

	CST = 16.d0*datan(1.d0)*bigg

	do i = 1,10
		t1 = d(i,1)/3.d0
		t2 = d(i,2)/(Rt*4.d0)
		t3 = d(i,3)/((Rt**2)*5.d0)		
		t4 = d(i,4)/((Rt**3)*6.d0)
		r2 = r(i+1)
		r1 = r(i)

		cumul(i) =  t1*(r2**3) + t2*(r2**4) + t3*(r2**5) + t4*(r2**6) - &
			  ( t1*(r1**3) + t2*(r1**4) + t3*(r1**5) + t4*(r1**6) )
	
	enddo

	if ( x > r(11) ) x = r(11)

	do i =1,10
 	if ((x > r(i)).and.(x <= r(i+1))) then

		GR = 0.d0
		do j = 1,i-1
			GR = GR + cumul(j)
		enddo

		t1 = d(i,1)/3.d0
		t2 = d(i,2)/(Rt*4.d0)
		t3 = d(i,3)/((Rt**2)*5.d0)		
		t4 = d(i,4)/((Rt**3)*6.d0)
		r2 = x
		r1 = r(i)

		GR = GR +   t1*(r2**3) + t2*(r2**4) + t3*(r2**5) + t4*(r2**6)	&
			- ( t1*(r1**3) + t2*(r1**4) + t3*(r1**5) + t4*(r1**6) )
		GR = GR * CST
		if (x > epsilon) GR = GR/(x**2)
		return
	endif
	enddo
!
  end subroutine prem_light2
!


  subroutine prem_light(x0,ro,vp,vs,Qmu,flag_g,GR,nbi,ri)
!---------------------------------------
!
!   Given non-dimensionalized radius x0, prem returns
!   density, ro, compressional velocity, vp, and shear 
!   velocity, vs according to the PREM model (1s ref period)
!     (Anderson and Dziewonski, 1981).  
!
!   Also returns gravity gr
!
!
      integer 				:: i
      doubleprecision    		:: x0,x,ro,vp,vs,Qmu
      logical				:: flag_g
      doubleprecision,optional   	:: GR
      doubleprecision,parameter   	:: Rt = 6371000.d0,epsilon = 1.d-9
      doubleprecision,parameter		:: bigg = 6.6723d-11
      doubleprecision,dimension(7)	:: r,q
      doubleprecision,dimension(6,4)	:: d,p,s
      doubleprecision,dimension(6)	:: cumul
      logical				:: pastrouve
      integer, optional :: nbi
      doubleprecision,optional,dimension(7) :: ri
!
      doubleprecision			:: CST,r1,r2,t1,t2,t3,t4
      integer				:: j
!
   ro = 0.d0; vp = 0.d0; vs = 0.d0;

!
!!$   r(1 ) = 0.d0
!!$   r(2 ) = 1221.5d0
!!$   r(3 ) = 3480.d0
!!$   r(4 ) = 3630.d0
!!$   r(5 ) = 5600.d0
!!$   r(6 ) = 5701.d0
!!$   r(7 ) = 5771.d0
!!$   r(8 ) = 5971.d0
!!$   r(9 ) = 6151.d0
!!$   r(10) = 6291.d0
!!$   r(11) = 6346.6d0
!!$   r(12) = 6356.d0
!!$   r(13) = 6368.d0
!!$   r(14) = 6371.d0
   r(1 ) = 0.d0
   r(2 ) = 1221.5d0
   r(3 ) = 3480.d0
   r(4 ) = 5701.d0
   r(5 ) = 5971.d0
   r(6 ) = 6151.d0
   r(7 ) = 6371.d0
!attenuation
  q(1)=84.6d0; q(2)=1.d5; q(3)=312.d0; q(4)=312.d0;
  q(5)=143.d0; q(6)=600.d0
!
  if (present(nbi)) then
     nbi=7
  endif
  if (present(ri)) then
     ri(:)=r(:)
  endif
!
  d(:,:) = 0.d0
!
  d(1,1) = 13.0885d0; 			  d(1,3) = -8.8381d0 
  d(2,1) = 12.5815d0; d(2,2) = -1.2638d0; d(2,3) = -3.6426d0; d(2,4) = -5.5281d0
  d(3,1) = 7.9565d0 ; d(3,2) = -6.4761;   d(3,3) = 5.5283d0;  d(3,4) = -3.0807d0
  d(4,1) = 11.2494d0; d(4,2) = -8.0298d0
  d(5,1) = 7.1089d0 ; d(5,2) = -3.8045d0
  d(6,1) = 2.6910d0 ; d(6,2) = 0.6924d0
!
!----
!
  p(:,:) = 0.d0
!
  p(1,1) = 11.2622d0 ; p(1,3) = -6.3640d0
  p(2,1) = 11.0487d0 ; p(2,2) = -4.0362d0; p(2,3)  = 4.8023d0; p(2,4) = -13.5732d0
  p(3,1) = 24.952d0 ; p(3,2)  = -40.4673d0; p(3,3) = 51.4832d0; p(3,4) = -26.6419d0
  p(4,1) = 39.7027d0 ; p(4,2)  = -32.6166d0
  p(5,1) = 20.3926d0 ; p(5,2)  = -12.2569d0
  p(6,1) = 4.1875d0 ; p(6,2)  = 3.9382d0
!----
!
  s(:,:) = 0.d0
!
  s(1,1) = 3.6678d0; s(1,3) = -4.4475d0
  s(3,1) = 11.1671d0; s(3,2) = -13.7818d0; s(3,3) = 17.4575d0; s(3,4) = -9.2777d0
  s(4,1) = 22.3512d0; s(4,2) = -18.5856d0 
  s(5,1) = 8.9496d0; s(5,2) = -4.4597
  s(6,1) = 2.1519d0; s(6,2) = 2.3481d0
!
  r(:) = r(:) * 1000.d0
  x    = x0   * Rt
!
  pastrouve = .true.
!
 do i = 1,6

 if ( pastrouve.and.(x >= r(7)) ) then
  x0 = 1.d0
  ro = d(6,1) + x0*( d(6,2) + x0*( d(6,3) + x0*( d(6,4) )))
  vp = p(6,1) + x0*( p(6,2) + x0*( p(6,3) + x0*( p(6,4) )))
  vs = s(6,1) + x0*( s(6,2) + x0*( s(6,3) + x0*( s(6,4) )))
  Qmu=q(6)
  pastrouve = .false.  
 endif

 if ( pastrouve.and.((x >= r(i)).and.(x < r(i+1))) ) then
  ro = d(i,1) + x0*( d(i,2) + x0*( d(i,3) + x0*( d(i,4) )))
  vp = p(i,1) + x0*( p(i,2) + x0*( p(i,3) + x0*( p(i,4) )))
  vs = s(i,1) + x0*( s(i,2) + x0*( s(i,3) + x0*( s(i,4) )))
  Qmu=q(i)
  pastrouve = .false.
 endif

 enddo

  if (.not. flag_g) return

!-- Gravity

	CST = 16.d0*datan(1.d0)*bigg

	do i = 1,6
		t1 = d(i,1)/3.d0
		t2 = d(i,2)/(Rt*4.d0)
		t3 = d(i,3)/((Rt**2)*5.d0)		
		t4 = d(i,4)/((Rt**3)*6.d0)
		r2 = r(i+1)
		r1 = r(i)

		cumul(i) =  t1*(r2**3) + t2*(r2**4) + t3*(r2**5) + t4*(r2**6) - &
			  ( t1*(r1**3) + t2*(r1**4) + t3*(r1**5) + t4*(r1**6) )
	
	enddo

	if ( x > r(7) ) x = r(7)

	do i =1,6
 	if ((x > r(i)).and.(x <= r(i+1))) then

		GR = 0.d0
		do j = 1,i-1
			GR = GR + cumul(j)
		enddo

		t1 = d(i,1)/3.d0
		t2 = d(i,2)/(Rt*4.d0)
		t3 = d(i,3)/((Rt**2)*5.d0)		
		t4 = d(i,4)/((Rt**3)*6.d0)
		r2 = x
		r1 = r(i)

		GR = GR +   t1*(r2**3) + t2*(r2**4) + t3*(r2**5) + t4*(r2**6)	&
			- ( t1*(r1**3) + t2*(r1**4) + t3*(r1**5) + t4*(r1**6) )
		GR = GR * CST
		if (x > epsilon) GR = GR/(x**2)
		return
	endif
	enddo
!
  end subroutine prem_light

  subroutine prem_vlight(x0,ro,vp,vs,flag_g,GR,nbi,ri)
!---------------------------------------
!
!   Given non-dimensionalized radius x0, prem returns
!   density, ro, compressional velocity, vp, and shear 
!   velocity, vs according to the PREM model (1s ref period)
!     (Anderson and Dziewonski, 1981).  
!
!   Also returns gravity gr
!
!
      integer 				:: i
      doubleprecision    		:: x0,x,ro,vp,vs
      logical				:: flag_g
      doubleprecision,optional   	:: GR
      doubleprecision,parameter   	:: Rt = 6371000.d0,epsilon = 1.d-9
      doubleprecision,parameter		:: bigg = 6.6723d-11
      doubleprecision,dimension(7)	:: r
      doubleprecision,dimension(6,4)	:: d,p,s
      doubleprecision,dimension(6)	:: cumul
      logical				:: pastrouve
      integer, optional :: nbi
      doubleprecision,optional,dimension(7) :: ri
!
      doubleprecision			:: CST,r1,r2,t1,t2,t3,t4
      integer				:: j
!
   ro = 0.d0; vp = 0.d0; vs = 0.d0;

!
   r(1 ) = 0.d0
   r(2 ) = 1221.5d0
   r(3 ) = 3480.d0
   r(4 ) = 5701.d0
   r(5 ) = 5971.d0
   r(6 ) = 6151.d0
   r(7 ) = 6371.d0
!
  if (present(nbi)) then
     nbi=7
  endif
  if (present(ri)) then
     ri(:)=r(:)
  endif
!
  d(:,:) = 0.d0
!
  d(1,1) = 12.d0; 			  
  d(2,1) = 10.d0;
  d(3,1) = 7.9565d0 ; d(3,2) = -6.4761;   d(3,3) = 5.5283d0;  d(3,4) = -3.0807d0
  d(4,1) = 11.2494d0; d(4,2) = -8.0298d0
  d(5,1) = 7.1089d0 ; d(5,2) = -3.8045d0
  d(6,1) = 2.6910d0 ; d(6,2) = 0.6924d0
!
!----
!
  p(:,:) = 0.d0
!
  p(1,1) = 11.d0
  p(2,1) = 8.d0
  p(3,1) = 24.952d0 ; p(3,2)  = -40.4673d0; p(3,3) = 51.4832d0; p(3,4) = -26.6419d0
  p(4,1) = 39.7027d0 ; p(4,2)  = -32.6166d0
  p(5,1) = 20.3926d0 ; p(5,2)  = -12.2569d0
  p(6,1) = 4.1875d0 ; p(6,2)  = 3.9382d0
!----
!
  s(:,:) = 0.d0
!
  s(1,1) = 4.d0
  s(3,1) = 11.1671d0; s(3,2) = -13.7818d0; s(3,3) = 17.4575d0; s(3,4) = -9.2777d0
  s(4,1) = 22.3512d0; s(4,2) = -18.5856d0 
  s(5,1) = 8.9496d0; s(5,2) = -4.4597
  s(6,1) = 2.1519d0; s(6,2) = 2.3481d0
!
  r(:) = r(:) * 1000.d0
  x    = x0   * Rt
!
  pastrouve = .true.
!
 do i = 1,6

 if ( pastrouve.and.(x >= r(7)) ) then
  x0 = 1.d0
  ro = d(6,1) + x0*( d(6,2) + x0*( d(6,3) + x0*( d(6,4) )))
  vp = p(6,1) + x0*( p(6,2) + x0*( p(6,3) + x0*( p(6,4) )))
  vs = s(6,1) + x0*( s(6,2) + x0*( s(6,3) + x0*( s(6,4) )))
  pastrouve = .false.
 endif

 if ( pastrouve.and.((x >= r(i)).and.(x < r(i+1))) ) then
  ro = d(i,1) + x0*( d(i,2) + x0*( d(i,3) + x0*( d(i,4) )))
  vp = p(i,1) + x0*( p(i,2) + x0*( p(i,3) + x0*( p(i,4) )))
  vs = s(i,1) + x0*( s(i,2) + x0*( s(i,3) + x0*( s(i,4) )))
  pastrouve = .false.
 endif

 enddo

  if (.not. flag_g) return

!-- Gravity

	CST = 16.d0*datan(1.d0)*bigg

	do i = 1,6
		t1 = d(i,1)/3.d0
		t2 = d(i,2)/(Rt*4.d0)
		t3 = d(i,3)/((Rt**2)*5.d0)		
		t4 = d(i,4)/((Rt**3)*6.d0)
		r2 = r(i+1)
		r1 = r(i)

		cumul(i) =  t1*(r2**3) + t2*(r2**4) + t3*(r2**5) + t4*(r2**6) - &
			  ( t1*(r1**3) + t2*(r1**4) + t3*(r1**5) + t4*(r1**6) )
	
	enddo

	if ( x > r(7) ) x = r(7)

	do i =1,6
 	if ((x > r(i)).and.(x <= r(i+1))) then

		GR = 0.d0
		do j = 1,i-1
			GR = GR + cumul(j)
		enddo

		t1 = d(i,1)/3.d0
		t2 = d(i,2)/(Rt*4.d0)
		t3 = d(i,3)/((Rt**2)*5.d0)		
		t4 = d(i,4)/((Rt**3)*6.d0)
		r2 = x
		r1 = r(i)

		GR = GR +   t1*(r2**3) + t2*(r2**4) + t3*(r2**5) + t4*(r2**6)	&
			- ( t1*(r1**3) + t2*(r1**4) + t3*(r1**5) + t4*(r1**6) )
		GR = GR * CST
		if (x > epsilon) GR = GR/(x**2)
		return
	endif
	enddo
!
  end subroutine prem_vlight

  subroutine prem_svlight(x0,ro,vp,vs,flag_g,GR,nbi,ri)
!---------------------------------------
!
!   Given non-dimensionalized radius x0, prem returns
!   density, ro, compressional velocity, vp, and shear 
!   velocity, vs according to the PREM model (1s ref period)
!     (Anderson and Dziewonski, 1981).  
!
!   Also returns gravity gr
!
!
      integer 				:: i
      doubleprecision    		:: x0,x,ro,vp,vs
      logical				:: flag_g
      doubleprecision,optional   	:: GR
      doubleprecision,parameter   	:: Rt = 6371000.d0,epsilon = 1.d-9
      doubleprecision,parameter		:: bigg = 6.6723d-11
      doubleprecision,dimension(7)	:: r
      doubleprecision,dimension(6,4)	:: d,p,s
      doubleprecision,dimension(6)	:: cumul
      logical				:: pastrouve
      integer, optional :: nbi
      doubleprecision,optional,dimension(7) :: ri
!
      doubleprecision			:: CST,r1,r2,t1,t2,t3,t4
      integer				:: j
!
   ro = 0.d0; vp = 0.d0; vs = 0.d0;

!
   r(1 ) = 0.d0
   r(2 ) = 1221.5d0
   r(3 ) = 3480.d0
   r(4 ) = 5701.d0
   r(5 ) = 5971.d0
   r(6 ) = 6151.d0
   r(7 ) = 6371.d0
!
  if (present(nbi)) then
     nbi=7
  endif
  if (present(ri)) then
     ri(:)=r(:)
  endif
!
  d(:,:) = 0.d0
!
  d(1,1) = 12.d0; 			  
  d(2,1) = 10.d0;
  d(3,1) = 4.9d0  
  d(4,1) = 3.9d0 
  d(5,1) = 3.480d0  
  d(6,1) = 3.360d0  
!
!----
!
  p(:,:) = 0.d0
!
  p(1,1) = 11.d0
  p(2,1) = 8.d0
  p(3,1) = 12.4d0  
  p(4,1) = 9.90d0 
  p(5,1) = 8.7d0 
  p(6,1) = 8.00d0 
!----
!
  s(:,:) = 0.d0
!
  s(1,1) = 4.d0
  s(3,1) = 6.8d0
  s(4,1) = 5.35d0
  s(5,1) = 4.7d0
  s(6,1) = 4.45d0
!
  r(:) = r(:) * 1000.d0
  x    = x0   * Rt
!
  pastrouve = .true.
!
 do i = 1,6

 if ( pastrouve.and.(x >= r(7)) ) then
  x0 = 1.d0
  ro = d(6,1) + x0*( d(6,2) + x0*( d(6,3) + x0*( d(6,4) )))
  vp = p(6,1) + x0*( p(6,2) + x0*( p(6,3) + x0*( p(6,4) )))
  vs = s(6,1) + x0*( s(6,2) + x0*( s(6,3) + x0*( s(6,4) )))
  pastrouve = .false.
 endif

 if ( pastrouve.and.((x >= r(i)).and.(x < r(i+1))) ) then
  ro = d(i,1) + x0*( d(i,2) + x0*( d(i,3) + x0*( d(i,4) )))
  vp = p(i,1) + x0*( p(i,2) + x0*( p(i,3) + x0*( p(i,4) )))
  vs = s(i,1) + x0*( s(i,2) + x0*( s(i,3) + x0*( s(i,4) )))
  pastrouve = .false.
 endif

 enddo

  if (.not. flag_g) return

!-- Gravity

	CST = 16.d0*datan(1.d0)*bigg

	do i = 1,6
		t1 = d(i,1)/3.d0
		t2 = d(i,2)/(Rt*4.d0)
		t3 = d(i,3)/((Rt**2)*5.d0)		
		t4 = d(i,4)/((Rt**3)*6.d0)
		r2 = r(i+1)
		r1 = r(i)

		cumul(i) =  t1*(r2**3) + t2*(r2**4) + t3*(r2**5) + t4*(r2**6) - &
			  ( t1*(r1**3) + t2*(r1**4) + t3*(r1**5) + t4*(r1**6) )
	
	enddo

	if ( x > r(7) ) x = r(7)

	do i =1,6
 	if ((x > r(i)).and.(x <= r(i+1))) then

		GR = 0.d0
		do j = 1,i-1
			GR = GR + cumul(j)
		enddo

		t1 = d(i,1)/3.d0
		t2 = d(i,2)/(Rt*4.d0)
		t3 = d(i,3)/((Rt**2)*5.d0)		
		t4 = d(i,4)/((Rt**3)*6.d0)
		r2 = x
		r1 = r(i)

		GR = GR +   t1*(r2**3) + t2*(r2**4) + t3*(r2**5) + t4*(r2**6)	&
			- ( t1*(r1**3) + t2*(r1**4) + t3*(r1**5) + t4*(r1**6) )
		GR = GR * CST
		if (x > epsilon) GR = GR/(x**2)
		return
	endif
	enddo
!
  end subroutine prem_svlight
!
  subroutine modele_4couches(x0,ro,vp,vs,flag_g,GR,nbi,ri)
!---------------------------------------
!
      integer 				:: i
      doubleprecision    		:: x0,x,ro,vp,vs
      logical				:: flag_g
      doubleprecision,optional   	:: GR
      doubleprecision,parameter   	:: Rt = 6371000.d0,epsilon = 1.d-9
      doubleprecision,parameter		:: bigg = 6.6723d-11
      doubleprecision,dimension(5)	:: r
      doubleprecision,dimension(4,4)	:: d,p,s
      doubleprecision,dimension(4)	:: cumul
      logical				:: pastrouve
      integer, optional :: nbi
      doubleprecision,optional,dimension(5) :: ri
!
      doubleprecision			:: CST,r1,r2,t1,t2,t3,t4
      integer				:: j
!
   ro = 0.d0; vp = 0.d0; vs = 0.d0;

!
   r(1 ) = 0.d0
   r(2 ) = 1221.5d0
   r(3 ) = 3480.d0
   r(4 ) = 5406.d0
   r(5 ) = 6371.d0
!
  if (present(nbi)) then
     nbi=5
  endif
  if (present(ri)) then
     ri(:)=r(:)
  endif
!
  d(:,:) = 0.d0
!
  d(1,1) = 12.d0; 			  
  d(2,1) = 10.d0;
  d(3,1) = 4.9d0  
  d(4,1) = 3.9d0 
!
!----
!
  p(:,:) = 0.d0
!
  p(1,1) = 11.d0
  p(2,1) = 8.d0
  p(3,1) = 12.4d0  
  p(4,1) = 9.90d0 
!----
!
  s(:,:) = 0.d0
!
  s(1,1) = 4.d0
  s(3,1) = 6.8d0
  s(4,1) = 5.35d0
!
  r(:) = r(:) * 1000.d0
  x    = x0   * Rt
!
  pastrouve = .true.
!
 do i = 1,6

 if ( pastrouve.and.(x >= r(5)) ) then
  x0 = 1.d0
  ro = d(4,1) + x0*( d(4,2) + x0*( d(4,3) + x0*( d(4,4) )))
  vp = p(4,1) + x0*( p(4,2) + x0*( p(4,3) + x0*( p(4,4) )))
  vs = s(4,1) + x0*( s(4,2) + x0*( s(4,3) + x0*( s(4,4) )))
  pastrouve = .false.
 endif

 if ( pastrouve.and.((x >= r(i)).and.(x < r(i+1))) ) then
  ro = d(i,1) + x0*( d(i,2) + x0*( d(i,3) + x0*( d(i,4) )))
  vp = p(i,1) + x0*( p(i,2) + x0*( p(i,3) + x0*( p(i,4) )))
  vs = s(i,1) + x0*( s(i,2) + x0*( s(i,3) + x0*( s(i,4) )))
  pastrouve = .false.
 endif

 enddo

  if (.not. flag_g) return

!-- Gravity

	CST = 16.d0*datan(1.d0)*bigg

	do i = 1,4
		t1 = d(i,1)/3.d0
		t2 = d(i,2)/(Rt*4.d0)
		t3 = d(i,3)/((Rt**2)*5.d0)		
		t4 = d(i,4)/((Rt**3)*6.d0)
		r2 = r(i+1)
		r1 = r(i)

		cumul(i) =  t1*(r2**3) + t2*(r2**4) + t3*(r2**5) + t4*(r2**6) - &
			  ( t1*(r1**3) + t2*(r1**4) + t3*(r1**5) + t4*(r1**6) )
	
	enddo

	if ( x > r(5) ) x = r(5)

	do i =1,4
 	if ((x > r(i)).and.(x <= r(i+1))) then

		GR = 0.d0
		do j = 1,i-1
			GR = GR + cumul(j)
		enddo

		t1 = d(i,1)/3.d0
		t2 = d(i,2)/(Rt*4.d0)
		t3 = d(i,3)/((Rt**2)*5.d0)		
		t4 = d(i,4)/((Rt**3)*6.d0)
		r2 = x
		r1 = r(i)

		GR = GR +   t1*(r2**3) + t2*(r2**4) + t3*(r2**5) + t4*(r2**6)	&
			- ( t1*(r1**3) + t2*(r1**4) + t3*(r1**5) + t4*(r1**6) )
		GR = GR * CST
		if (x > epsilon) GR = GR/(x**2)
		return
	endif
	enddo
!
  end subroutine modele_4couches
!
!!$  subroutine prem_aniso(r_in,rho,vpv,vph,vsv,vsh,eta_aniso,Qmu,nbi,rio)
!!$
!!$  implicit none
!!$
!!$  integer, optional :: nbi
!!$  doubleprecision,optional,dimension(12) :: rio
!!$! given a radius r, gives the density rho,
!!$! speeds vp and vs, and the quality factor Qmu
!!$
!!$  double precision r_in,rho,Qmu,vpv,vph,vsv,vsh,eta_aniso,ri
!!$
!!$  double precision x
!!$
!!$! R_EARTH: radius of Earth (m)
!!$! RMOHO: radius of the Moho (m)
!!$! R220: radius of 220 km discontinuity (m)
!!$! R400: radius of 400 km discontinuity (m)
!!$! R600: radius of 600 km 2nd order discontinuity (m)
!!$! R670: radius of 670 km discontinuity (m)
!!$! R771: radius of 771 km 2nd order discontinuity (m)
!!$! RTOPDDOUBLEPRIME: radius of top of D" 2nd order discontinuity (m)
!!$! RCMB: radius of CMB (m)
!!$! RICB: radius of ICB (m)
!!$  double precision, parameter :: R_EARTH = 6371000.d0
!!$  double precision, parameter :: RMOHO = 6346600.d0
!!$  double precision, parameter :: R220 = 6151000.d0
!!$  double precision, parameter :: R400 = 5971000.d0
!!$  double precision, parameter :: R600 = 5771000.d0
!!$  double precision, parameter :: R670 = 5701000.d0
!!$  double precision, parameter :: R771 = 5600000.d0
!!$  double precision, parameter :: RTOPDDOUBLEPRIME = 3630000.d0
!!$  double precision, parameter :: RCMB = 3480000.d0
!!$  double precision, parameter :: RICB = 1221000.d0
!!$
!!$  if (present(nbi)) then
!!$     nbi=12
!!$  endif
!!$  if (present(rio)) then
!!$     rio(1)=0.0d0;rio(2)=RICB;rio(3)=RCMB;rio(4)=RTOPDDOUBLEPRIME;rio(5)=R771
!!$     rio(6)=R670;rio(7)=R600;rio(8)= R400;rio(9)=R220;rio(10)=6291000.d0;rio(11)=RMOHO;rio(12)=R_EARTH
!!$     rio(:)=rio(:)/1000.d0
!!$  endif
!!$! clear anisotropy for isotropic layers
!!$  eta_aniso = 1.d0
!!$
!!$! normalized radius
!!$  ri = r_in*R_EARTH
!!$  x = ri / R_EARTH
!!$  
!!$
!!$!
!!$!--- inner core
!!$!
!!$  if(ri >= 0.d0 .and. ri <= RICB) then
!!$    rho=13.0885d0-8.8381d0*x*x
!!$    vpv=11.2622d0-6.3640d0*x*x
!!$    vsv=3.6678d0-4.4475d0*x*x
!!$    vph=vpv
!!$    vsh=vsv
!!$    Qmu=84.6d0
!!$!
!!$!--- outer core
!!$!
!!$  else if(ri > RICB .and. ri <= RCMB) then
!!$    rho=12.5815d0-1.2638d0*x-3.6426d0*x*x-5.5281d0*x*x*x
!!$    vpv=11.0487d0-4.0362d0*x+4.8023d0*x*x-13.5732d0*x*x*x
!!$    vsv=0.0d0
!!$    vph=vpv
!!$    vsh=vsv
!!$    Qmu=1.d5
!!$!
!!$!--- D" at the base of the mantle
!!$!
!!$  else if(ri > RCMB .and. ri <= RTOPDDOUBLEPRIME) then
!!$    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
!!$    vpv=15.3891d0-5.3181d0*x+5.5242d0*x*x-2.5514d0*x*x*x
!!$    vsv=6.9254d0+1.4672d0*x-2.0834d0*x*x+0.9783d0*x*x*x
!!$    vph=vpv
!!$    vsh=vsv
!!$    Qmu=312.0d0
!!$!
!!$!--- mantle: from top of D" to d670
!!$!
!!$  else if(ri > RTOPDDOUBLEPRIME .and. ri <= R771) then
!!$    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
!!$    vpv=24.9520d0-40.4673d0*x+51.4832d0*x*x-26.6419d0*x*x*x
!!$    vsv=11.1671d0-13.7818d0*x+17.4575d0*x*x-9.2777d0*x*x*x
!!$    vph=vpv
!!$    vsh=vsv
!!$    Qmu=312.0d0
!!$  else if(ri > R771 .and. ri <= R670) then
!!$    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
!!$    vpv=29.2766d0-23.6027d0*x+5.5242d0*x*x-2.5514d0*x*x*x
!!$    vsv=22.3459d0-17.2473d0*x-2.0834d0*x*x+0.9783d0*x*x*x
!!$    vph=vpv
!!$    vsh=vsv
!!$    Qmu=312.0d0
!!$!
!!$!--- mantle: above d670
!!$!
!!$  else if(ri > R670 .and. ri <= R600) then
!!$    rho=5.3197d0-1.4836d0*x
!!$    vpv=19.0957d0-9.8672d0*x
!!$    vsv=9.9839d0-4.9324d0*x
!!$    vph=vpv
!!$    vsh=vsv
!!$!    Qmu=143.0d0
!!$    Qmu=312.0d0
!!$  else if(ri > R600 .and. ri <= R400) then
!!$    rho=11.2494d0-8.0298d0*x
!!$    vpv=39.7027d0-32.6166d0*x
!!$    vsv=22.3512d0-18.5856d0*x
!!$    vph=vpv
!!$    vsh=vsv
!!$!    Qmu=143.0d0
!!$    Qmu=156.0d0
!!$  else if(ri > R400 .and. ri <= R220) then
!!$    rho=7.1089d0-3.8045d0*x
!!$    vpv=20.3926d0-12.2569d0*x
!!$    vsv=8.9496d0-4.4597d0*x
!!$    vph=vpv
!!$    vsh=vsv
!!$!    Qmu=143.0d0
!!$    Qmu=156.0d0
!!$  else if(ri > R220 .and. ri <= 6291.0d3) then
!!$
!!$! anisotropy in PREM only above 220 km
!!$
!!$    rho=2.6910d0+0.6924d0*x
!!$    vpv=0.8317d0+7.2180d0*x
!!$    vph=3.5908d0+4.6172d0*x
!!$    vsv=5.8582d0-1.4678d0*x
!!$    vsh=-1.0839d0+5.7176d0*x
!!$    eta_aniso=3.3687d0-2.4778d0*x
!!$    Qmu=80.0d0
!!$!    Qmu=1.0d5
!!$
!!$  else
!!$! use PREM crust
!!$    if(ri > 6291.0d3 .and. ri <= RMOHO) then
!!$
!!$! anisotropy in PREM only above 220 km
!!$
!!$      rho=2.6910d0+0.6924d0*x
!!$      vpv=0.8317d0+7.2180d0*x
!!$      vph=3.5908d0+4.6172d0*x
!!$      vsv=5.8582d0-1.4678d0*x
!!$      vsh=-1.0839d0+5.7176d0*x
!!$      eta_aniso=3.3687d0-2.4778d0*x
!!$      Qmu=600.0d0
!!$!      Qmu=1.0d5
!!$
!!$! no anisotropy in the crust in PREM
!!$
!!$    else if(ri > RMOHO .and. ri <= 6356.0d3) then
!!$
!!$! same properties everywhere in PREM crust (only one layer in the crust)
!!$        rho=2.6d0
!!$        vpv=5.8d0
!!$        vsv=3.2d0
!!$        vph=vpv
!!$        vsh=vsv
!!$        Qmu=600.0d0
!!$!        Qmu=1.0d5
!!$
!!$    else if(ri > 6356.0d3 .and. ri <= 6368.0d3) then
!!$      rho=2.6d0
!!$      vpv=5.8d0
!!$      vsv=3.2d0
!!$      vph=vpv
!!$      vsh=vsv
!!$      Qmu=600.0d0
!!$!      Qmu=1.0d5
!!$    else if(ri > 6368.0d3 .and. ri <= R_EARTH) then
!!$      rho=2.6d0
!!$      vpv=5.8d0
!!$      vsv=3.2d0
!!$      vph=vpv
!!$      vsh=vsv
!!$      Qmu=600.0d0
!!$!      Qmu=1.0d5
!!$!
!!$!       replace water layer with crust
!!$!
!!$!        rho=1.02d0
!!$!        vpv=1.45d0
!!$!        vsv=0.0d0
!!$!        vph=vpv
!!$!        vsh=vsv
!!$!        Qmu=0.0d0
!!$    else
!!$      stop 'radius unreachable in subroutine prem'
!!$    endif
!!$  endif
!!$!test
!!$!           vph=vpv
!!$!           vsh=vsv
!!$!           eta_aniso=1.d0
!!$!fin test
!!$
!!$! return values in S.I. units
!!$
!!$!  rho=rho*1000.0d0
!!$!  vpv=vpv*1000.0d0
!!$!  vsv=vsv*1000.0d0
!!$!  vph=vph*1000.0d0
!!$!  vsh=vsh*1000.0d0
!!$
!!$  end subroutine prem_aniso
  subroutine prem_aniso(x0,ro,vpv,vph,vsv,vsh,eta,Qmu,nbi,ri)
!---------------------------------------
!
!   Given non-dimensionalized radius x0, prem returns
!   density, ro, compressional velocity, vp, and shear 
!   velocity, vs according to the PREM model (1s ref period)
!     (Anderson and Dziewonski, 1981).  
!
!   Also returns gravity gr
!
!
      integer 				:: i
      doubleprecision    		:: x0,x,ro,vpv,vph,vsv,vsh,eta
      doubleprecision   	:: Qmu
      doubleprecision,parameter   	:: Rt = 6371000.d0,epsilon = 1.d-9
      doubleprecision,parameter		:: bigg = 6.6723d-11
      doubleprecision,dimension(14)	:: r,q
      doubleprecision,dimension(13,4)	:: d,pv,ph,sv,sh,et
      doubleprecision,dimension(13)	:: cumul
      logical				:: pastrouve
!
      doubleprecision			:: CST,r1,r2,t1,t2,t3,t4
      integer				:: j
      integer, optional :: nbi
      doubleprecision,optional,dimension(14) :: ri
!

!
   r(1) = 0.d0;     r(2) = 1221.5d0; r(3) = 3480.d0; r(4) = 3630.d0
   r(5) = 5600.d0;  r(6) = 5701.d0;  r(7) = 5771.d0; r(8) = 5971.d0
   r(9) = 6151.d0; r(10) = 6291.d0; r(11) = 6346.6d0
  r(12) = 6356.d0; r(13) = 6368.d0; r(14) = 6371.d0
!
  if (present(nbi)) then
     nbi=14
  endif
  if (present(ri)) then
     ri(:)=r(:)
  endif
!
!  q(1)=84.6d0; q(2)=1.d5; q(3)=312.d0; q(4)=312.d0;q(5)=312.d0;
!  q(6)=312.d0; q(7)=156.d0; q(8)=156.d0;q(9)=70.d0;q(10)=191.d0;
!  q(11)=300.d0
!  q(6)=312.d0; q(7)=143.d0; q(8)=143.d0;q(9)=80.d0;q(10)=80.d0;
!  q(11)=600.d0; q(12)=600.d0;q(13)=600.d0
!QL6:
  q(1)=84.6d0; q(2)=1.d5; q(3)=355.d0; q(4)=355.d0;q(5)=355.d0;
  q(6)=165.d0; q(7)=165.d0; q(8)=165.d0;q(9)=70.d0;q(10)=191.d0;
  q(11)=300.d0; q(12)=300.d0;q(13)=300.d0
!
  d(:,:) = 0.d0
!
  d(1,1) = 13.0885d0; 			  d(1,3) = -8.8381d0 
  d(2,1) = 12.5815d0; d(2,2) = -1.2638d0; d(2,3) = -3.6426d0; d(2,4) = -5.5281d0
  d(3,1) = 7.9565d0 ; d(3,2) = -6.4761;   d(3,3) = 5.5283d0;  d(3,4) = -3.0807d0
  d(4,1) = 7.9565d0 ; d(4,2) = -6.4761;   d(4,3) = 5.5283d0;  d(4,4) = -3.0807d0
  d(5,1) = 7.9565d0 ; d(5,2) = -6.4761;   d(5,3) = 5.5283d0;  d(5,4) = -3.0807d0
  d(6,1) = 5.3197d0 ; d(6,2) = -1.4836d0
  d(7,1) = 11.2494d0; d(7,2) = -8.0298d0
  d(8,1) = 7.1089d0 ; d(8,2) = -3.8045d0
  d(9,1) = 2.6910d0 ; d(9,2) = 0.6924d0
  d(10,1) = 2.6910d0; d(10,2) = 0.6924d0
  d(11,1) = 2.9d0  
  d(12,1) = 2.6d0    
  d(13,1) = 2.6d0    
!  d(11,1) = 2.7d0  
!!  d(11,1) = 2.715d0  
!
!----
!
  pv(:,:) = 0.d0
!
  pv(1,1) = 11.2622d0 ; pv(1,3) = -6.3640d0
  pv(2,1) = 11.0487d0 ; pv(2,2) = -4.0362d0; pv(2,3)  = 4.8023d0; pv(2,4) = -13.5732d0
  pv(3,1) = 15.3891d0 ; pv(3,2) = -5.3181d0; pv(3,3)  = 5.5242d0; pv(3,4) = -2.5514d0
  pv(4,1) = 24.952d0 ; pv(4,2)  = -40.4673d0; pv(4,3) = 51.4832d0; pv(4,4) = -26.6419d0
  pv(5,1) = 29.2766d0 ; pv(5,2) = -23.6027d0; pv(5,3) = 5.5242d0; pv(5,4) = -2.5514d0
  pv(6,1) = 19.0957d0 ; pv(6,2)  = -9.8672d0
  pv(7,1) = 39.7027d0 ; pv(7,2)  = -32.6166d0
  pv(8,1) = 20.3926d0 ; pv(8,2)  = -12.2569d0
!  pv(9,1) = 4.1875d0 ; pv(9,2)  = 3.9382d0
!  pv(10,1) = 4.1875d0 ; pv(10,2) = 3.9382d0
  pv(9,1) =  0.8317d0; pv(9,2)  = 7.2180d0
  pv(10,1) =  0.8317d0 ; pv(10,2) =7.2180d0
  pv(11,1) = 6.8d0 
  pv(12,1) = 5.8d0
  pv(13,1) = pv(12,1)

!  pv(11,1) = 6.15d0 
!  pv(11,1) = 6.18d0 
!
!----
  ph(:,:)=pv(:,:)
  ph(9,1) = 3.5908d0 ; ph(9,2)  = 4.6172d0
  ph(10,1) = 3.5908d0 ; ph(10,2)  = 4.6172d0
!
  sv(:,:) = 0.d0
!
  sv(1,1) = 3.6678d0; sv(1,3) = -4.4475d0

  sv(3,1) = 6.9254d0; sv(3,2) = 1.4672d0; sv(3,3) = -2.0834d0; sv(3,4) = 0.9783d0
  sv(4,1) = 11.1671d0; sv(4,2) = -13.7818d0; sv(4,3) = 17.4575d0; sv(4,4) = -9.2777d0
  sv(5,1) = 22.3459d0; sv(5,2) = -17.2473d0; sv(5,3) = -2.0834d0; sv(5,4) = 0.9783d0
  sv(6,1) = 9.9839; sv(6,2) = -4.9324
  sv(7,1) = 22.3512d0; sv(7,2) = -18.5856d0 
  sv(8,1) = 8.9496d0; sv(8,2) = -4.4597
!  sv(9,1) = 2.1519d0; sv(9,2) = 2.3481d0
!  sv(10,1) = 2.1519d0; sv(10,2) = 2.3481d0
  sv(9,1) = 5.8582d0; sv(9,2) = -1.4678d0
  sv(10,1) = 5.8582d0; sv(10,2) = -1.4678d0
  sv(11,1) = 3.9d0 
  sv(12,1) = 3.2d0 
  sv(13,1) = sv(12,1)
!  sv(11,1) = 3.32d0 
!!  sv(11,1) = 3.47d0 
!
  sh(:,:)=sv(:,:)
!  sh(9,1) = -1.0839d0 ; sh(9,2)  = 5.7176d0
!  sh(10,1) = -1.0839d0 ; sh(10,2)  = 5.7176d0
!
  et(:,:)=0.0d0
  et(:,1)=1.d0
  et(9,1)=3.3687d0 ; et(9,2)=-2.4778d0
  et(10,1)=3.3687d0 ; et(10,2)=-2.4778d0
!
  r(:) = r(:) * 1000.d0
  x    = x0   * Rt
!
  pastrouve = .true.
!
 do i = 1,13
 if ( pastrouve.and.(x >= r(14)) ) then
  x0 = 1.d0
  ro = d(13,1) + x0*( d(13,2) + x0*( d(13,3) + x0*( d(13,4) )))
  vpv = pv(13,1) + x0*( pv(13,2) + x0*( pv(13,3) + x0*( pv(13,4) )))
  vsv = sv(13,1) + x0*( sv(13,2) + x0*( sv(13,3) + x0*( sv(13,4) )))
  vph = ph(13,1) + x0*( ph(13,2) + x0*( ph(13,3) + x0*( ph(13,4) )))
  vsh = sh(13,1) + x0*( sh(13,2) + x0*( sh(13,3) + x0*( sh(13,4) )))
  eta = et(13,1) + x0*( et(13,2) + x0*( et(13,3) + x0*( et(13,4) )))
  Qmu=q(13)
  pastrouve = .false.
 endif

 if ( pastrouve.and.((x >= r(i)).and.(x < r(i+1))) ) then
  ro  = d (i,1) + x0*( d (i,2) + x0*( d (i,3) + x0*( d (i,4) )))
  vpv = pv(i,1) + x0*( pv(i,2) + x0*( pv(i,3) + x0*( pv(i,4) )))
  vsv = sv(i,1) + x0*( sv(i,2) + x0*( sv(i,3) + x0*( sv(i,4) )))
  vph = ph(i,1) + x0*( ph(i,2) + x0*( ph(i,3) + x0*( ph(i,4) )))
  vsh = sh(i,1) + x0*( sh(i,2) + x0*( sh(i,3) + x0*( sh(i,4) )))
  eta = et(i,1) + x0*( et(i,2) + x0*( et(i,3) + x0*( et(i,4) )))
  Qmu=q(i)
  pastrouve = .false.
 endif

 enddo
!
  end subroutine prem_aniso
!

  subroutine emiask(x0,ro,vp,vs,Qalpha,Qbeta,GFLAG,GRAV,nbi,ri)
!------------------------------------------
!
!   Emiask returns model parameters for the IASPEI working model
!   (September 1990.1).
!   Given non-dimensionalized radius x0, emiasp returns
!   non-dimensionalized density, ro, compressional velocity, vp, and
!   shear velocity, vs.  Non-dimensionalization is according to the
!   scheme of Gilbert in program EOS:  x0 by a (the radius of the
!   Earth), ro by robar (the mean density of the Earth), and velocity
!   by a*sqrt(pi*G*robar) (where G is the universal gravitational
!   constant.
!
!---> Correction pour rho,vp,vs non normalises !
!
!
      doubleprecision    		:: x0,ro,vp,vs,Qbeta,Qalpha
      logical,optional			:: GFLAG
      doubleprecision,optional		:: GRAV
      integer, optional :: nbi
      doubleprecision,optional,dimension(12) :: ri
!
      doubleprecision    		:: xn,rn,vn,x,x1,Rt
      integer 				:: i
      save
      doubleprecision,dimension(12)	:: r
      doubleprecision,dimension(11,4)	:: d,p,s
      doubleprecision,dimension(11  )	:: Qm,Qe

      data r/0.      ,1217.1  ,3482.0  ,3631.  ,5611.   ,5711.   ,             &
      5961.   ,6161.   ,6251.   ,6336.   ,6351.    ,6371.    /

      data d/13.01219,12.58416, 6.8143 , 6.8143 , 6.8143 ,11.11978,            &
        7.15855, 7.15855, 7.15855,  2.92  , 2.72   ,                      &
              0.     ,-1.69929,-1.66273,-1.66273,-1.66273,-7.87054,            &
       -3.85999,-3.85999,-3.85999,2*0.,                                        &
             -8.45292,-1.94128,-1.18531,-1.18531,-1.18531,6*0.,                &
              0.     ,-7.11215,9*0./

      data p/11.12094,10.03904,14.49470,25.1486 ,25.969838,29.38896,           &
       30.78765,25.41389, 8.785412, 6.5   , 5.8   ,                      &
              0.   , 3.75665, -1.47089,-41.1538, -16.934118,-21.40656,         &
      -23.25415,-17.69722,-0.7495294, 2*0.,                                    &
             -4.09689,-13.67046, 0.0  ,51.9932,7*0.,                           &
              0.     , 0.      , 0.     ,-26.6083,7*0./  

      data s/ 3.56454, 0.      , 8.16616,12.9303 ,20.768902,17.70732,          &
       15.24213,5.750203, 6.706232, 3.75   , 3.36   ,                     &
              0.     , 0.      ,-1.58206,-21.2590,-16.531471,-13.50652,        &
      -11.08553,-1.274202,-2.248585, 2*0.,                                     &
             -3.45241, 0.      , 0.0    ,27.8988 ,7*0.,                        &
              0.     , 0.      , 0.     ,-14.1080,7*0./ 
      data Qm/11*300./

      data Qe/11*600./

      data xn,rn,vn/6371.,.18125793,.14598326/,i/1/
!
!
      if (present(nbi)) then
         nbi=12
      endif
      if (present(ri)) then
         ri(:)=r(:)
      endif
!
      Rt=r(12)
      x=dmax1(x0,0.0d0)
      x1=xn*x
 2    if(x1.ge.r(i)) go to 1
      i=i-1
      go to 2
 1    if(x1.le.r(i+1).or.i.ge.11) go to 3
      i=i+1
      if(i.lt.11) go to 1
 3    ro = d(i,1)+x*(d(i,2)+x*(d(i,3)+x*d(i,4)))
      vp = p(i,1)+x*(p(i,2)+x*(p(i,3)+x*p(i,4)))
      vs = s(i,1)+x*(s(i,2)+x*(s(i,3)+x*s(i,4)))
      Qbeta  = Qm(i)
      Qalpha = Qe(i)

      return
!
  end subroutine emiask  
  subroutine premc(x0,ro,vp,vs,Qmu,flag_g,GR,nbi,ri)
!---------------------------------------
!
!   Given non-dimensionalized radius x0, prem returns
!   density, ro, compressional velocity, vp, and shear 
!   velocity, vs according to the PREM model (1s ref period)
!     (Anderson and Dziewonski, 1981).  
!
!   Also returns gravity gr
!
!
      integer, optional :: nbi
      doubleprecision,optional,dimension(12) :: ri
      integer 				:: i
      doubleprecision    		:: x0,x,ro,vp,vs
      logical				:: flag_g
      doubleprecision,optional   	:: GR,Qmu
      doubleprecision,parameter   	:: Rt = 6371000.d0,epsilon = 1.d-9
      doubleprecision,parameter		:: bigg = 6.6723d-11
      doubleprecision,dimension(12)	:: r,q
      doubleprecision,dimension(11,4)	:: d,p,s
      doubleprecision,dimension(11)	:: cumul
      logical				:: pastrouve
!
      doubleprecision			:: CST,r1,r2,t1,t2,t3,t4
      integer				:: j
!
   ro = 0.0d0; vp = 0.d0; vs = 0.d0;

!
   r(1) = 0.d0;     r(2) = 1221.5d0; r(3) = 3480.d0; r(4) = 3630.d0
   r(5) = 5600.d0;  r(6) = 5701.d0;  r(7) = 5771.d0; r(8) = 5971.d0
   r(9) = 6151.d0; r(10) = 6291.d0; r(11) = 6346.6d0
   r(12) = 6371.d0
  if (present(nbi)) then
     nbi=12
  endif
  if (present(ri)) then
     ri(:)=r(:)
  endif
!
  q(1)=84.6d0; q(2)=1.d5; q(3)=312.d0; q(4)=312.d0;q(5)=312.d0;
  q(6)=312.d0; q(7)=156.d0; q(8)=156.d0;q(9)=70.d0;q(10)=191.d0;
  q(11)=300.d0
!
!  q(1)=84.6d0; q(2)=1.d5; q(3)=355.d0; q(4)=355.d0;q(5)=355.d0;
!  q(6)=165.d0; q(7)=165.d0; q(8)=165.d0;q(9)=70.d0;q(10)=191.d0;
!  q(11)=300.d0
!
  d(:,:) = 0.d0
!
  d(1,1) = 13.0885d0; 			  d(1,3) = -8.8381d0 
  d(2,1) = 12.5815d0; d(2,2) = -1.2638d0; d(2,3) = -3.6426d0; d(2,4) = -5.5281d0
  d(3,1) = 7.9565d0 ; d(3,2) = -6.4761;   d(3,3) = 5.5283d0;  d(3,4) = -3.0807d0
  d(4,1) = 7.9565d0 ; d(4,2) = -6.4761;   d(4,3) = 5.5283d0;  d(4,4) = -3.0807d0
  d(5,1) = 7.9565d0 ; d(5,2) = -6.4761;   d(5,3) = 5.5283d0;  d(5,4) = -3.0807d0
  d(6,1) = 5.3197d0 ; d(6,2) = -1.4836d0
  d(7,1) = 11.2494d0; d(7,2) = -8.0298d0
  d(8,1) = 7.1089d0 ; d(8,2) = -3.8045d0
  d(9,1) = 2.6910d0 ; d(9,2) = 0.6924d0
  d(10,1) = 2.6910d0; d(10,2) = 0.6924d0
  d(11,1) = 2.7d0  
!  d(11,1) = 2.715d0  
!
!----
!
  p(:,:) = 0.d0
!
  p(1,1) = 11.2622d0 ; p(1,3) = -6.3640d0
  p(2,1) = 11.0487d0 ; p(2,2) = -4.0362d0; p(2,3)  = 4.8023d0; p(2,4) = -13.5732d0
  p(3,1) = 15.3891d0 ; p(3,2) = -5.3181d0; p(3,3)  = 5.5242d0; p(3,4) = -2.5514d0
  p(4,1) = 24.952d0 ; p(4,2)  = -40.4673d0; p(4,3) = 51.4832d0; p(4,4) = -26.6419d0
  p(5,1) = 29.2766d0 ; p(5,2) = -23.6027d0; p(5,3) = 5.5242d0; p(5,4) = -2.5514d0
  p(6,1) = 19.0957d0 ; p(6,2)  = -9.8672d0
  p(7,1) = 39.7027d0 ; p(7,2)  = -32.6166d0
  p(8,1) = 20.3926d0 ; p(8,2)  = -12.2569d0
  p(9,1) = 4.1875d0 ; p(9,2)  = 3.9382d0
  p(10,1) = 4.1875d0 ; p(10,2) = 3.9382d0
  p(11,1) = 6.15d0 
!  p(11,1) = 6.18d0 
!
!----
!
  s(:,:) = 0.d0
!
  s(1,1) = 3.6678d0; s(1,3) = -4.4475d0

  s(3,1) = 6.9254d0; s(3,2) = 1.4672d0; s(3,3) = -2.0834d0; s(3,4) = 0.9783d0
  s(4,1) = 11.1671d0; s(4,2) = -13.7818d0; s(4,3) = 17.4575d0; s(4,4) = -9.2777d0
  s(5,1) = 22.3459d0; s(5,2) = -17.2473d0; s(5,3) = -2.0834d0; s(5,4) = 0.9783d0
  s(6,1) = 9.9839; s(6,2) = -4.9324
  s(7,1) = 22.3512d0; s(7,2) = -18.5856d0 
  s(8,1) = 8.9496d0; s(8,2) = -4.4597
  s(9,1) = 2.1519d0; s(9,2) = 2.3481d0
  s(10,1) = 2.1519d0; s(10,2) = 2.3481d0
  s(11,1) = 3.32d0 
!  s(11,1) = 3.47d0 
!
!
  r(:) = r(:) * 1000.d0
  x    = x0   * Rt
!
  pastrouve = .true.
!
 do i = 1,11
 if ( pastrouve.and.(x >= r(12)) ) then
  x0 = 1.d0
  ro = d(11,1) + x0*( d(11,2) + x0*( d(11,3) + x0*( d(11,4) )))
  vp = p(11,1) + x0*( p(11,2) + x0*( p(11,3) + x0*( p(11,4) )))
  vs = s(11,1) + x0*( s(11,2) + x0*( s(11,3) + x0*( s(11,4) )))
  Qmu=q(11)
  pastrouve = .false.
 endif

 if ( pastrouve.and.((x >= r(i)).and.(x < r(i+1))) ) then
  ro = d(i,1) + x0*( d(i,2) + x0*( d(i,3) + x0*( d(i,4) )))
  vp = p(i,1) + x0*( p(i,2) + x0*( p(i,3) + x0*( p(i,4) )))
  vs = s(i,1) + x0*( s(i,2) + x0*( s(i,3) + x0*( s(i,4) )))
  Qmu=q(i)
  pastrouve = .false.
 endif

 enddo
  if (.not. flag_g) return

!-- Gravity

	CST = 16.d0*datan(1.d0)*bigg

	do i = 1,11
		t1 = d(i,1)/3.d0
		t2 = d(i,2)/(Rt*4.d0)
		t3 = d(i,3)/((Rt**2)*5.d0)		
		t4 = d(i,4)/((Rt**3)*6.d0)
		r2 = r(i+1)
		r1 = r(i)

		cumul(i) =  t1*(r2**3) + t2*(r2**4) + t3*(r2**5) + t4*(r2**6) - &
			  ( t1*(r1**3) + t2*(r1**4) + t3*(r1**5) + t4*(r1**6) )
	
	enddo

	if ( x > r(12) ) x = r(12)

	do i =1,11
 	if ((x > r(i)).and.(x <= r(i+1))) then

		GR = 0.d0
		do j = 1,i-1
			GR = GR + cumul(j)
		enddo

		t1 = d(i,1)/3.d0
		t2 = d(i,2)/(Rt*4.d0)
		t3 = d(i,3)/((Rt**2)*5.d0)		
		t4 = d(i,4)/((Rt**3)*6.d0)
		r2 = x
		r1 = r(i)

		GR = GR +   t1*(r2**3) + t2*(r2**4) + t3*(r2**5) + t4*(r2**6)	&
			- ( t1*(r1**3) + t2*(r1**4) + t3*(r1**5) + t4*(r1**6) )
		GR = GR * CST
		if (x > epsilon) GR = GR/(x**2)
		return
	endif
	enddo
!
  end subroutine premc
  subroutine prem_test(x0,ro,vp,vs,Qmu,flag_g,GR,nbi,ri)
!---------------------------------------
!
!   Given non-dimensionalized radius x0, prem returns
!   density, ro, compressional velocity, vp, and shear 
!   velocity, vs according to the PREM model (1s ref period)
!     (Anderson and Dziewonski, 1981).  
!
!   Also returns gravity gr
!
!
      integer, optional :: nbi
      doubleprecision,optional,dimension(9) :: ri
      integer 				:: i
      doubleprecision    		:: x0,x,ro,vp,vs
      logical				:: flag_g
      doubleprecision,optional   	:: GR,Qmu
      doubleprecision,parameter   	:: Rt = 6371000.d0,epsilon = 1.d-9
      doubleprecision,parameter		:: bigg = 6.6723d-11
      doubleprecision,dimension(9)	:: r,q
      doubleprecision,dimension(8,4)	:: d,p,s
      doubleprecision,dimension(8)	:: cumul
      logical				:: pastrouve
!
      doubleprecision			:: CST,r1,r2,t1,t2,t3,t4
      integer				:: j
!
   ro = 0.0d0; vp = 0.d0; vs = 0.d0;

!
   r(1) = 0.d0;     r(2) = 1221.5d0; r(3) = 3480.d0; r(4) = 3630.d0
   r(5) = 5600.d0;  r(6) = 5701.d0;  r(7) = 5771.d0; 
   r(8) = 5971.d0
   r(9) = 6371.d0
  if (present(nbi)) then
     nbi=9
  endif
  if (present(ri)) then
     ri(:)=r(:)
  endif
!
!  q(1)=84.6d0; q(2)=1.d5; q(3)=312.d0; q(4)=312.d0;q(5)=312.d0;
!  q(6)=312.d0; q(7)=156.d0; q(8)=156.d0;q(9)=70.d0;q(10)=191.d0;
!  q(11)=300.d0
!
  q(1)=84.6d0; q(2)=1.d5; q(3)=355.d0; q(4)=355.d0;q(5)=355.d0;
  q(6)=165.d0; q(7)=165.d0; q(8)=300.d0; q(9)=300.d0
!
  d(:,:) = 0.d0
!
  d(1,1) = 13.0885d0; 			  d(1,3) = -8.8381d0 
  d(2,1) = 12.5815d0; d(2,2) = -1.2638d0; d(2,3) = -3.6426d0; d(2,4) = -5.5281d0
  d(3,1) = 7.9565d0 ; d(3,2) = -6.4761;   d(3,3) = 5.5283d0;  d(3,4) = -3.0807d0
  d(4,1) = 7.9565d0 ; d(4,2) = -6.4761;   d(4,3) = 5.5283d0;  d(4,4) = -3.0807d0
  d(5,1) = 7.9565d0 ; d(5,2) = -6.4761;   d(5,3) = 5.5283d0;  d(5,4) = -3.0807d0
  d(6,1) = 5.3197d0 ; d(6,2) = -1.4836d0
  d(7,1) = 11.2494d0; d(7,2) = -8.0298d0
  d(8,1) = 3.5d0
!
!----
!
  p(:,:) = 0.d0
!
  p(1,1) = 11.2622d0 ; p(1,3) = -6.3640d0
  p(2,1) = 11.0487d0 ; p(2,2) = -4.0362d0; p(2,3)  = 4.8023d0; p(2,4) = -13.5732d0
  p(3,1) = 15.3891d0 ; p(3,2) = -5.3181d0; p(3,3)  = 5.5242d0; p(3,4) = -2.5514d0
  p(4,1) = 24.952d0 ; p(4,2)  = -40.4673d0; p(4,3) = 51.4832d0; p(4,4) = -26.6419d0
  p(5,1) = 29.2766d0 ; p(5,2) = -23.6027d0; p(5,3) = 5.5242d0; p(5,4) = -2.5514d0
  p(6,1) = 19.0957d0 ; p(6,2)  = -9.8672d0
  p(7,1) = 39.7027d0 ; p(7,2)  = -32.6166d0
  p(8,1) = 8.55
!
!----
!
  s(:,:) = 0.d0
!
  s(1,1) = 3.6678d0; s(1,3) = -4.4475d0

  s(3,1) = 6.9254d0; s(3,2) = 1.4672d0; s(3,3) = -2.0834d0; s(3,4) = 0.9783d0
  s(4,1) = 11.1671d0; s(4,2) = -13.7818d0; s(4,3) = 17.4575d0; s(4,4) = -9.2777d0
  s(5,1) = 22.3459d0; s(5,2) = -17.2473d0; s(5,3) = -2.0834d0; s(5,4) = 0.9783d0
  s(6,1) = 9.9839; s(6,2) = -4.9324
  s(7,1) = 22.3512d0; s(7,2) = -18.5856d0 
  s(8,1) = 4.72
!
  r(:) = r(:) * 1000.d0
  x    = x0   * Rt
!
  pastrouve = .true.
!
 do i = 1,8
 if ( pastrouve.and.(x >= r(9)) ) then
  x0 = 1.d0
  ro = d(8,1) + x0*( d(8,2) + x0*( d(8,3) + x0*( d(8,4) )))
  vp = p(8,1) + x0*( p(8,2) + x0*( p(8,3) + x0*( p(8,4) )))
  vs = s(8,1) + x0*( s(8,2) + x0*( s(8,3) + x0*( s(8,4) )))
  Qmu=q(8)
  pastrouve = .false.
 endif

 if ( pastrouve.and.((x >= r(i)).and.(x < r(i+1))) ) then
  ro = d(i,1) + x0*( d(i,2) + x0*( d(i,3) + x0*( d(i,4) )))
  vp = p(i,1) + x0*( p(i,2) + x0*( p(i,3) + x0*( p(i,4) )))
  vs = s(i,1) + x0*( s(i,2) + x0*( s(i,3) + x0*( s(i,4) )))
  Qmu=q(i)
  pastrouve = .false.
 endif

 enddo
  if (.not. flag_g) return

!-- Gravity

	CST = 16.d0*datan(1.d0)*bigg

	do i = 1,8
		t1 = d(i,1)/3.d0
		t2 = d(i,2)/(Rt*4.d0)
		t3 = d(i,3)/((Rt**2)*5.d0)		
		t4 = d(i,4)/((Rt**3)*6.d0)
		r2 = r(i+1)
		r1 = r(i)

		cumul(i) =  t1*(r2**3) + t2*(r2**4) + t3*(r2**5) + t4*(r2**6) - &
			  ( t1*(r1**3) + t2*(r1**4) + t3*(r1**5) + t4*(r1**6) )
	
	enddo

	if ( x > r(9) ) x = r(9)

	do i =1,8
 	if ((x > r(i)).and.(x <= r(i+1))) then

		GR = 0.d0
		do j = 1,i-1
			GR = GR + cumul(j)
		enddo

		t1 = d(i,1)/3.d0
		t2 = d(i,2)/(Rt*4.d0)
		t3 = d(i,3)/((Rt**2)*5.d0)		
		t4 = d(i,4)/((Rt**3)*6.d0)
		r2 = x
		r1 = r(i)

		GR = GR +   t1*(r2**3) + t2*(r2**4) + t3*(r2**5) + t4*(r2**6)	&
			- ( t1*(r1**3) + t2*(r1**4) + t3*(r1**5) + t4*(r1**6) )
		GR = GR * CST
		if (x > epsilon) GR = GR/(x**2)
		return
	endif
	enddo
!
  end subroutine prem_test

!---------------------------------------------------------------------
end program generate_prem
!---------------------------------------------------------------------
