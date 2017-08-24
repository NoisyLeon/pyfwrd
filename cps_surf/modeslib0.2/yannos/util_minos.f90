!--------------------------------------------------------------------------
module minos
!--------------------------------------------------------------------------
private
public :: grav,startl,fsbm,sfbm,derms,dermf,trknt,tprop,spsm,fpsm  &
         ,sdepth,eifout,fprpmn,sprpmn,detqn,model,clean_minos_memory    &
         ,extension,write_fctpR,write_fctpL,check_eigenfreq2            &
         ,check_eigenfreq,zknt
logical :: first_in_startl=.true.,reallocation=.false.,allocation=.true.
integer :: nmx_backup=0,ncm_backup=0
!
contains
!###################################################################
!
! les routines suivantes sont celles de minos sans ces putains de goto
! (le prochain que je vois coder avec des goto, je le bute)
! yann Capdeville 28/10/98
!
!###################################################################
! 

!-------------------------------------------------------------------------
  subroutine model(mod)
!lecture du model d'entre, et calcul des constante, vitesses, des 
!numeros de couches etc..
!-------------------------------------------------------------------------
    use bits
    use rindex
    use param_modele
    use global_prameter
    use earth_modele
    use yannos_flag
!      
    implicit none
!
    type(modele), intent(in) :: mod
!
    character(len=80) :: title
    doubleprecision, dimension(:), allocatable :: vpv,vph,vsv,vsh,eta,wrk
!
    doubleprecision :: bigg=6.6723d-11,tau=1.d3,rhobar=5515.d0    
    doubleprecision :: rx,r1,r2,dr,rt,val,gn,vn2,rat
    integer :: ifdeck,ibidon,i,nreg,j,jj,knt,nn,ind,nlay,err
!
! model par couches:
!
       n        =mod%nbcou
       ncm=n
!       ncm=3000
       allocate(vsv(ncm),stat=err)
       allocate(eta(ncm),stat=err) 
       allocate(wrk(10*ncm),stat=err) 
       allocate(vsh(ncm),stat=err)
       allocate(vpv(ncm),stat=err)
       allocate(vph(ncm),stat=err)
       if (err/=0) stop ' model2 : pb a l''alloc'
!
!allocation des tableaux:
       if (allocation) call allocate_all()
!
       ifanis   =mod%ifanis
       tref     =mod%tref
       ifdeck   =mod%ifdeck
       nic      =mod%nic
       noc      =mod%noc
       r(1:n)     =mod%r(1:n)
       rho(1:n)   =mod%rho(1:n)
       vpv(1:n)   =mod%vpv(1:n)
       vsv(1:n)   =mod%vsv(1:n)
       qkappa(1:n)=mod%qkappa(1:n)
       qshear(1:n)=mod%qshear(1:n)
       vph(1:n)   =mod%vph(1:n)
       vsh(1:n)   =mod%vsh(1:n) 
       eta(1:n)   =mod%eta(1:n)
!
       ray(1:n)=r(1:n)
!
!------------------
    if(ifanis==0) then
!isotrope
       do  i=1,n
          vph(i)=vpv(i)
          vsh(i)=vsv(i)
          eta(i)=1.d0
       enddo
    endif
!
!***  normalise and spline ***
!
    rn=r(n)
    gn=pi*bigg*rhobar*rn
    vn2=gn*rn
    vn=dsqrt(vn2)
    wn=vn/rn
    do i=1,n
       r(i)=r(i)/rn
       if(i.gt.1) then
          if (dabs(r(i)-r(i-1)).lt.1.d-7) r(i)=r(i-1)
       endif
       if(qshear(i).gt.0.d0) qshear(i)=1.d0/qshear(i)
       if(qkappa(i).gt.0.d0) qkappa(i)=1.d0/qkappa(i)
       rho(i)=rho(i)/rhobar
       acon(i)=rho(i)*vph(i)*vph(i)/vn2
       ccon(i)=rho(i)*vpv(i)*vpv(i)/vn2
       lcon(i)=rho(i)*vsv(i)*vsv(i)/vn2
       ncon(i)=rho(i)*vsh(i)*vsh(i)/vn2
       fcon(i)=eta(i)*(acon(i)-2.d0*lcon(i))
       fmu(i)=(acon(i)+ccon(i)-2.d0*fcon(i)+5.d0*ncon(i)+6.d0*lcon(i))/15.d0
       flam(i)=(4.d0*(acon(i)+fcon(i)-ncon(i))+ccon(i))/9.d0-2.d0*fmu(i)/3.d0
       rat=4.d0*fmu(i)/(3.d0*(flam(i)+2.d0*fmu(i)))
       xlam(i)=((1.d0-rat)*qkappa(i)-.5d0*rat*qshear(i))/(1.d0-1.5d0*rat)
       xa2(i)=(1.d0-rat)*qkappa(i)+rat*qshear(i)
    enddo
    call drspln(1,n,r,rho,qro,wrk)
!     
!*** compute g *****
!    
    call grav(g,rho,qro,r,n)
!
    if (cancel_gravity) g(1:n)=0.0d0
!
    call drspln(1,n,r,g,qg,wrk)
    call drspln(1,n,r,fcon,fspl,wrk)
    call drspln(1,n,r,lcon,lspl,wrk)
    if(ifanis/=0) then
!anisotrope:
       call drspln(1,n,r,acon,aspl,wrk)
       call drspln(1,n,r,ccon,cspl,wrk)
       call drspln(1,n,r,ncon,nspl,wrk)
    endif
    nsl=n
    do while (vsv(nsl)<=0.d0) 
       nsl=nsl-1 
    enddo
    if (nsl==1.and.(vsv(nsl)<=0.d0))  nsl=0
    if (nsl==nic) nsl=n !ca veut dire qu'il n'y a qu'une couche d''eau
    nicp1=nic+1
    nocp1=noc+1
    nslp1=nsl+1
    tref=0.5d0*tref/pi
!
    deallocate(vpv,vph,vsv,vsh,eta,wrk)
!
!-------------------------------------------------------------------------
  end subroutine model
!-------------------------------------------------------------------------
!------------------------------------------------------------------
subroutine check_eigenfreq(l_in,wp,flag)
!------------------------------------------------------------------
  use bits, only: l,fl,fl1,fl2,fl3,sfl3,knsw,jcom,eps
  use shanks, only:maxo,stepf
  implicit none
!
  integer, intent(in) :: l_in
  doubleprecision, intent(in) :: wp
  logical, intent(out) :: flag
!
  integer, parameter :: sens=0,inss=5
  integer :: ibid
  doubleprecision :: w1,w2,w3,w4,det1,det2,det3,det4,rap
!
  l=l_in  
!
  fl=l
  fl1=fl+1.d0
  fl2=fl+fl1
  fl3=fl*fl1 
  sfl3=dsqrt(fl3) 
  knsw=0 !kount swich off
  if (l>4) then
     maxo=inss
     stepf=1.d0
  else
     maxo=inss
     stepf=0.25d0
  endif
  w1=wp*(1d0-1.d0   *dsqrt(eps*0.01))
  w2=wp*(1d0+1.d0   *dsqrt(eps*0.01))
  w3=wp*(1d0-2.d0 *dsqrt(eps*0.01))
  w4=wp*(1d0+2.d0 *dsqrt(eps*0.01))
  
  call detqn(w1,ibid,det1,sens)
  call detqn(w2,ibid,det2,sens)
  call detqn(w3,ibid,det3,sens)
  call detqn(w4,ibid,det4,sens)
  rap=abs((det2-det1)*(w4-w3)/(det4-det3)/(w2-w1))
!telle que sont choisies les wi le rapport precedant doit valoir 1
!si la frequence propre en est une et 2 si c''est seuelement
! un changement de signe foireux du determinant
  if ( (abs(rap-1.d0)>0.1d0)) then
     flag=.true.
     print*,'rap false=',rap,l
  else
     flag=.false.
     print*,'rap true=',rap,l
  endif
  
!------------------------------------------------------------------
end subroutine check_eigenfreq
!------------------------------------------------------------------

!------------------------------------------------------------------
subroutine check_eigenfreq2(l_in,wp,flag,wdiff)
!------------------------------------------------------------------
  use bits, only: l,fl,fl1,fl2,fl3,sfl3,knsw,jcom,eps
  use shanks, only:maxo,stepf
  implicit none
!
  integer, intent(in) :: l_in
  doubleprecision, intent(in) :: wp
  logical, intent(out) :: flag
  doubleprecision, intent(in), optional :: wdiff
!
  integer, parameter :: sens=0,inss=5
  integer :: ibid
  doubleprecision :: w1,w2,check1,check2,bid,epsloc
!
  l=l_in  
!
  fl=l
  fl1=fl+1.d0
  fl2=fl+fl1
  fl3=fl*fl1 
  sfl3=dsqrt(fl3) 
  knsw=0 !kount swich off
  if (l>4) then
     maxo=inss
     stepf=1.d0
  else
     maxo=inss
     stepf=0.25d0
  endif
  if (present(wdiff)) then
     epsloc=min(max(dsqrt(eps*0.01d0),abs(wdiff)*10.d0),1.d-4)
  else
     epsloc=dsqrt(eps*0.01d0)
  endif
  w1=wp*(1d0-epsloc)
  w2=wp*(1d0+epsloc)
  
  call detqn(w1,ibid,bid,sens,check=check1)
  call detqn(w2,ibid,bid,sens,check=check2)
  if (dsign(1.d0,check1)*dsign(1.d0,check2) > 0.0d0 ) then
!c'est un mode de la graine ou de Stonley:
     flag=.true.
!     print*,'check false=',dsign(1.d0,check1)*dsign(1.d0,check2),l
  else
!c'est un mode normal:
     flag=.false.
!     print*,'check true =',dsign(1.d0,check1)*dsign(1.d0,check2),l
  endif
  
!------------------------------------------------------------------
end subroutine check_eigenfreq2
!------------------------------------------------------------------


!----------------------------------------------------------------------
subroutine write_fctpR(l_in,f_in,iper,iinfo,ieigen,ilost)
!----------------------------------------------------------------------
  use bits, only: l,jcom,qinv,cg,vn,fl,fl1,fl2,fl3,sfl3,knsw,nord,wray,wn,wgrav,cond_limite
  use shanks, only:maxo,stepf
  use def_gparam
  use eifx,only: ar
  use yannos_flag
  implicit none
  integer, intent(in) :: l_in,iper,iinfo,ieigen,ilost
  real(DP), intent(in) :: f_in
!
  real(DP) :: wp,gcom,qmod,bid,tcom,WMHZ,wdiff
  integer :: sens,ibid
  logical :: flag
!initialisation de ar du module eifx:
  ar(:,:)=0.0d0
!
  wp=f_in*2._DP*PI
  l=l_in
  if (l==0) then
     jcom=1
  else
     jcom=3
  endif
!
  fl=l
  fl1=fl+1.d0
  fl2=fl+fl1
  fl3=fl*fl1 
  sfl3=dsqrt(fl3)   
  knsw=0 !kount swich off
  maxo=8 ! precision rks ?
  sens=1  
  stepf=1.d0
!
  call detqn(wp,ibid,bid,sens)
  gcom=vn*cg/1000.d0
  qmod=0.d0
  if(qinv.gt.0.d0) qmod=1.d0/qinv
  TCOM=2.D0*PI/wp
  WMHZ=1000.D0/TCOM
  wdiff=(wp-wray*wn)/wp
  call check_eigenfreq2(l,wp ,flag,wdiff)       
!  WRITE(iper,200) NORD,'S',L,wp,WMHZ,TCOM,GCOM,QMOD,WDIFF,flag
!on est obliger de rappeler detqn (pour faire autrement, il faudrait
!backuper ar et inorm
  if (abs(wdiff)<seuil_ray.or..not.flag) then
     call detqn(wp,ibid,bid,sens)
     WRITE(iper,200) NORD,'S',L,wp,WMHZ,TCOM,GCOM,QMOD,WDIFF,flag
     call modout(wp,qmod,gcom,wdiff,flag,ieigen,iinfo)
  else
     if (cond_limite==1) then
        print*,'force_remedy',wdiff,flag
        call detqn(wp,ibid,bid,sens,force_remedy=.true.)
     endif
     wdiff=(wp-wray*wn)/wp
     if (abs(wdiff)<seuil_ray.or..not.flag) then
        WRITE(iper,200) NORD,'S',L,wp,WMHZ,TCOM,GCOM,QMOD,WDIFF,flag
        call modout(wp,qmod,gcom,wdiff,flag,ieigen,iinfo)        
     else
        if (keep_bad_modes) then
!on met des zero
!           ar(:,:)=0.0d0
           call modout(wp,qmod,gcom,wdiff,flag,ieigen,iinfo)
           WRITE(iper,200) NORD,'S',L,wp,WMHZ,TCOM,GCOM,QMOD,WDIFF,flag
        endif
        WRITE(ilost,200) NORD,'S',L,wp,WMHZ,TCOM,GCOM,QMOD,WDIFF,flag
        if (.not.keep_bad_modes) print*,'Mode de la graine ou de Stoneley ellimine: ',NORD,'S',l,wdiff
     endif
  endif
200 FORMAT(I5,A2,I5,6G16.7,l7)
!----------------------------------------------------------------------
end subroutine write_fctpR
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine write_fctpL(l_in,f_in,iper,iinfo,ieigen)
!----------------------------------------------------------------------
  use bits, only: l,jcom,qinv,cg,vn,fl,fl1,fl2,fl3,sfl3,knsw,nord,wray,wn
  use shanks, only:maxo,stepf
  use def_gparam
  use yannos_flag
  implicit none
  integer, intent(in) :: l_in,iper,iinfo,ieigen
  real(DP), intent(in) :: f_in
!
  real(DP) :: wp,gcom,qmod,bid,tcom,WMHZ,wdiff
  integer :: sens,ibid
!
  wp=f_in*2._DP*PI
  l=l_in
  if (l/=0) then
     jcom=2
!
     fl=l
     fl1=fl+1.d0
     fl2=fl+fl1
     fl3=fl*fl1 
     sfl3=dsqrt(fl3)   
     knsw=0 !kount swich off
     maxo=8 ! precision rks ?
     sens=1  
     stepf=1.0d0
     call detqn(wp,ibid,bid,sens)
     gcom=vn*cg/1000.d0
     qmod=0.d0
     if(qinv.gt.0.d0) qmod=1.d0/qinv
     TCOM=2.D0*PI/wp
     WMHZ=1000.D0/TCOM
     wdiff=(wp-wray*wn)/wp
     if (abs(wdiff)<seuil_ray) then
        WRITE(iper,200) NORD,'T',L,wp,WMHZ,TCOM,GCOM,QMOD,WDIFF
        call modout(wp,qmod,gcom,wdiff,.false.,ieigen,iinfo)
     else
        print*,'La precision  mauvaise, je le vire',wdiff
     endif
  endif
200 FORMAT(I5,A2,I5,6G16.7)
!----------------------------------------------------------------------
end subroutine write_fctpL
!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine detqn(wdim,knt,det,ifeif,Ann,continu,check,force_remedy)
!input:
! wdim: frequence en Hz a laquelle on integer les mineurs
! ifeif: =0 on integer de bas en haut seulement pour recuperer
!           le determinant (et le mode kount)
!        /=0 on integer aussi de haut en bas pour sortir la solution
! output:
! knt= nombre de branches entre 0 et wdim (clculer si demande, cf common)
! det= determinant de la traction la surface normalise
!
! input par common:
! knsw =1 On compte le nombre de branches
!      /=0 on ne compte pas
!
!**** supevises the integration of the equations,it returns the value
!**** of the secular determinant as det and the count of zero crossings.
!----------------------------------------------------------------------
    use bits
    use eifx, a=>ar
    use rindex
    use param_modele
    use detqn_buffer
    use yannos_flag
    implicit none
!
    doubleprecision, intent(in) :: wdim
    integer, intent(in)  :: ifeif
    integer, intent(out) :: knt
    doubleprecision, intent(out) :: det
    doubleprecision, dimension(5),optional, intent(out) :: Ann
    integer, optional :: continu
    logical, optional ::    force_remedy
    doubleprecision, optional, intent(out) :: check
!
    doubleprecision, dimension(14) :: ass
    doubleprecision, dimension( 4) :: zi
    integer :: nvefm,nvesm,is,irem,nbakf,nbaks,nto,n2,nb    &
              ,ls1,i,ls,iexp,jexp,j,im
    doubleprecision :: r10,dnorm,asi1,asi2,ff,cc,rnrm,q,m3  &
                      ,detp,y4a2,ll,scale,rhobar=5515.d0,y2a2    
!
! on calcul les pas d'integration (plus la freq est elevee, plus le precision l''est)
!    call steps(eps/(1.d4*(wdim/2.d0/pi)**2))
    call steps(epsint/(1.d4*(wdim/2.d0/pi)))
!si on passe la premiere fois dans detqn, on alloue vf:
    call allocate_vf()
    vf(:)=0.0d0
    a(:,:)=0.0d0
    ass(:)=0.0d0
!
    iback=0 ! on va de bas en haut pour l'integration
    w=wdim/wn !on normalise la frequence d'entree
    if (abs(wdim)<1.d-30) then
       print*,'wdim=',wdim,w,wn
       print*,'l=',l
       stop 'detqn: la frequence d''entree &
                                  &ne doit pas etre nulle !!!'
    endif
    wsq=w*w
    iexp=0
    kount=0
    kg=0 !kg=0 on ne tient pas compte de la perturbation gravite    
    fct=0.d0
    if(tref.gt.0.d0.and.use_tref) &
             fct=2.d0*dlog(tref*wdim)/pi !utilisation de periode
                                                   !de ref
!    print*,'dans detqn: l=',l,' wdim=',w,kg
!
!
!--------------------------------------------------------------
    select case(jcom)
!--------------------------------------------------------------
    case(3)
! modes spheroiaux
!--------------------------------------------------------------
       if(wdim.le.wgrav) kg=1 !on tient compte de la gravite
       nvefm=2+kg*3 !nombre de mineurs dans le fluide (2, sans gravite
                    ! ou 5 avec gravite)
       nvesm=5+kg*9 !nombre de mineurs dans le solide (5, sans gravite
                    ! ou 14 avec gravite)
!on determine a quelle profondeur commencer l'integration
! sdepth = start depth
       if (use_startlevel) then
          call sdepth(wdim,ls)
          if (ls>nocp1) then
!             print*,'ls >nocp1!'
             ls=(nocp1+nicp1)/2
          endif
       else
          ls=2
       endif
!
       if (ls <= 2) then
          r10=4.5d-4*(fl+.5d0)/wdim
          if(r10 < r(2)) then
             r(1)=r10
             if (.not.cancel_gravity) then
                g(1)=rho(1)*r(1)*1.333333333333333d0
             else
                g(1)=0.
             endif
             ls=1
          endif
       endif


! nocp1 : numero de la couche outer core +1
! nicp1 : numero de la couche inner core +1
!++++++++++++++++++++++++++++++++++++++++++
       if(ls <= nocp1 .and.noc/=0) then
!++++++++++++++++++++++++++++++++++++++++++
!------------------------------------------
          if(ls <= nicp1 .and.nic/=0) then
!------------------------------------------
!on est en dessous la graine

!calcule des solution de depart avec des fonctions de bessel spheriques
! dans une region solide:
             call spsm(ls,nvesm,ass)
!propage les nineurs jusqu'au le noyau liquide (nic)
             call sprpmn(ls,nic,ass,vf,nvesm,iexp)
             r(1)=0.d0
             g(1)=0.d0
!convertie les mineurs d'une region solide vers une region liquides
             call sfbm(ass,kg,iback)
!------------------------------------------
          endif
!------------------------------------------
          is=max0(ls,nicp1)
!si le point de depart (sortie de sdepth)  est dans le liquide
!on calcule  les determinant de depart avec des fct de bessel
!dans une region liquide
          if(is.eq.ls) call fpsm(ls,nvefm,ass)
!propage les nineurs jusqu'au manteau (noc)
          call fprpmn(is,noc,ass,vf,nvefm,iexp)
          if (noc==n) then
!on est deja a la surface
             if (present(continu).and.present(Ann)) then
                Ann(:)=0.0d0
                Ann(1)=a(2,n)/a(1,n)*rhobar*vn**2/rn
             endif
             select case (cond_limite)
             case(1)
!surface libre
                det=a(2,nsl)/dsqrt(a(1,nsl)*a(1,nsl)+a(2,nsl)*a(2,nsl))
!bord rigide
             case(2)
                det=a(1,nsl)/dsqrt(a(1,nsl)*a(1,nsl)+a(2,nsl)*a(2,nsl))
             end select
             if (present(check)) then
                check=a(1,nsl)*a(2,nsl)/dsqrt(a(1,n)*a(1,n)+a(2,n)*a(2,n))
             endif
          else
!convertie les mineurs d'une region fluide vers une region solide
             call fsbm(ass,kg,iback)
          endif



!++++++++++++++++++++++++++++++++++++++++++
       endif
!++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++
       if (noc/=n) then
!++++++++++++++++++++++++++++++++++++++++++
          is=max0(ls,nocp1)
!si le point de depart (sortie de sdepth)  est dans manteau
!on calcule  les determinants de depart avec des fct de bessel
!dans une region solide:
          if(is.eq.ls) call spsm(ls,nvesm,ass)
!propage les nineurs jusqu'a  la surface solide (base de l'ocean
! si il existe) :
          call sprpmn(is,nsl,ass,vf,nvesm,iexp)
! si nsl=n ca veut dire qu'il n'y a pas d'ocean:
          if(nsl==n) then
             dnorm=a(1,nsl)*a(1,nsl)
             do  i=2,nvesm
                dnorm=dnorm+a(i,nsl)*a(i,nsl)
             enddo
!
             select case (cond_limite)
             case(1)
!surface libre
                det=a(5,nsl)/dsqrt(dnorm)
!bord rigide
             case(2)
                det=a(1,nsl)/dsqrt(dnorm) 
             end select
             if (present(check)) then
                check=a(1,nsl)*a(5,nsl)/dsqrt(dnorm) 
             endif
!
!sortie pour l'operateur de couplage
             if (present(Ann)) then             
                if (present(continu)) then
!A00:
                   Ann(1)= a(4,nsl)
!A10:
                   Ann(2)=-a(2,nsl)/sqrt(2.d0)
!A01:
                   Ann(3)=-a(2,nsl)/sqrt(2.d0)
!A11:
                   Ann(4)=-a(3,nsl)/2.d0
!
                   Ann(:)=-Ann/a(1,nsl)*rhobar*vn**2/rn
                else
!A00:
                   Ann(1)= a(4,nsl)
!A10:
                   Ann(2)=-a(2,nsl)/sqrt(2.d0)
!A01:
                   Ann(3)=-a(2,nsl)/sqrt(2.d0)
!A11:
                   Ann(4)=-a(3,nsl)/2.d0
!
                   m3=a(3,nsl)     
                endif
             endif
          else
! si nsl/=n ca veut dire qu'il y a un ocean:
! dans se cas on convertie les mineurs de solide vers liquide
             call sfbm(ass,kg,iback)
! on propage dans l'ocean:
             call fprpmn(nslp1,n,ass,vf,nvefm,iexp)
!
             select case (cond_limite)
             case(1)
!surface libre
                if(kg.eq.0) det=a(2,n)/dsqrt(a(1,n)*a(1,n)+a(2,n)*a(2,n))
                if(kg.ne.0) det=a(5,n)/dsqrt(a(1,n)**2+a(2,n)**2+a(3,n)**2+ &
                     a(4,n)**2+a(5,n)**2)
!bord rigide
             case(2)
                det=a(1,n)/dsqrt(a(1,n)*a(1,n)+a(2,n)*a(2,n))
             end select
             if (present(check)) then
                if (kg.eq.0) then
                   check=a(1,n)*a(2,n)/dsqrt(a(1,n)*a(1,n)+a(2,n)*a(2,n))
                else
                   check=a(1,n)*a(5,n)/dsqrt(a(1,n)*a(1,n)+a(5,n)*a(5,n))
                endif
             endif
          endif
!+++++++++++++++++++++++++++++++++++++++++++++
       endif
!+++++++++++++++++++++++++++++++++++++++++++++
! si la couche de depart est au dessus de la CMB:
       if(ls.gt.noc.and.noc/=0) det=-det
! si mode count:
       if(knsw==1) then
             if(ls.gt.noc.and.noc/=0) kount=kount-2
             if (nic==0.and.noc==0) kount=kount-1
             irem=mod(kount,2)
             if (noc/=0) then
                if(irem.eq.0.and.det.lt.0.d0) kount=kount+1
                if(irem.ne.0.and.det.gt.0.d0) kount=kount+1
             endif
             knt=kount
       endif
!-----------------------------------
       if(ifeif/=0) then
! on a demander le calule des fonctions propres (et non uniquement
! du determinant)
!-----------------------------------
          iback=1 ! on va de haut en bas pour l'integration
          jexp=0
          nbakf=1+kg*3
          nbaks=4+kg*10
          do  i=1,nbaks
             ass(i)=0.d0
          enddo
!+++++++++++++++++++++++++++++++++++++++++++++++++
          if(n==nsl.and.noc/=n) then 
! pas d'ocean, on fixe les determinant a la surface
!+++++++++++++++++++++++++++++++++++++++++++++++++
             if(kg/=0)  then
!gravite
                asi1=a(3,n)*a(3,n)+a(12,n)*a(12,n)
                asi2=a(4,n)*a(4,n)+a(11,n)*a(11,n)
                if(asi2.le.asi1) ass(1)=dsign(1.d0,a(3,n))
                if(asi2.gt.asi1) ass(2)=dsign(1.d0,a(2,n))
             else
!pas de gravite
                asi1=a(3,n)*a(3,n)
                asi2=a(4,n)*a(4,n)
                select case (cond_limite)
                case(1)
!pour surface libre
                   if(asi2.le.asi1) ass(1)=dsign(1.d0,a(3,n))
                   if(asi2.gt.asi1) ass(2)=dsign(1.d0,a(2,n))
                case(2)
!pour bord rigide:
                   if(asi2.le.asi1) ass(3)=dsign(1.d0,a(3,n))
                   if(asi2.gt.asi1) ass(4)=dsign(1.d0,a(2,n))
                end select
             endif
!+++++++++++++++++++++++++++++++++++++++++++++++++
          else
!depart dans un liquide
!+++++++++++++++++++++++++++++++++++++++++++++++++
! ocean:
             if(kg/=0)  then
!gravite:
                asi1=a(3,n)*a(3,n)
                asi2=a(4,n)*a(4,n)
                if(asi2.le.asi1) ass(1)=dsign(1.d0,a(3,n))
                if(asi2.gt.asi1) ass(2)=dsign(1.d0,a(2,n))
                if (n/=nsl) then
!propage dans l'ocean
                   call fprpmn(n,nslp1,ass,vf,nbakf,jexp)
!conversion liquide solide:
                   call fsbm(ass,kg,iback)
                endif
             else
! pas de gravite
                select case(cond_limite)
                !pour bord libre:
                case(1)
                   ass(1)=dsign(1.d0,a(1,n))
                case(2)
                   !pour bord rigide:
!                   ass(2)=dsign(1.d0,a(2,n))
                   ass(1)=dsign(1.d0,a(2,n))
                end select
                if (n/=nsl) then
!propage dans l'ocean
                   call fprpmn(n,nslp1,ass,vf,nbakf,jexp)
!conversion liquide solide:
                   call fsbm(ass,kg,iback)
                endif
             endif
!+++++++++++++++++++++++++++++++++++++++++++++++++
          endif
!+++++++++++++++++++++++++++++++++++++++++++++++++
! nto est le numero de couche du manteau ou un numero 
! plus haut si ca ne vaut pas la peine d'aller plus bas (sortie de sdepth)
!+++++++++++++++++++++++++++++++++++++++++++++++++
          nto=max0(ls,nocp1)
          if (noc/=n) then
!on part dans le solide
!+++++++++++++++++++++++++++++++++++++++++++++++++
! on porpage dans le solide
             call sprpmn(nsl,nto,ass,vf,nbaks,jexp)
             if (noc/=0.and.nto/=ls) then
! si on doit aller dans le liquide: conversion:
                call sfbm(ass,kg,iback)
             endif
!+++++++++++++++++++++++++++++++++++++++++++++++++
          endif
          if (noc/=0.and.nto/=ls) then
             nto=max0(ls,nicp1)
!on est dans le liquide
!+++++++++++++++++++++++++++++++++++++++++++++++++             
! on propage de liquide
             call fprpmn(noc,nto,ass,vf,nbakf,jexp)
             if(nto/=ls.and.nic/=0) then
! si on doit aller dans la graine on convertie liquide solide
                call fsbm(ass,kg,iback)
                nto=max0(ls,2)
! on propage dans le solide jusqu'au centre ou ls:
                call sprpmn(nic,nto,ass,vf,nbaks,jexp)
             endif
          endif
          if (use_remedy.and.(abs(det) > 5.d-4 .or. present(force_remedy) )) &
               call remedy(ls)
! normalisation etc .... du mode :
          call eifout(ls)
!sortie pour l'operateur de couplage: divistion par la
! derivee du determinant
          if (present(Ann)) then
             if (noc/=n) then
                scale=1.0d0/dsqrt(rhobar*rn**3)
                ll=lcon(nsl)*(1.d0+qshear(nsl)*fct)
                y4a2=(w*ll*(a(4,nsl)-a(3,nsl)+a(1,nsl))    &
                              *(scale*rhobar*vn**2/rn))**2
                detp=-2.d0*w*wn*m3/y4a2/rn**2
                Ann(:)=Ann(:)/detp
             else
!interface liquide:
                scale=1.0d0/sqrt(rhobar*rn**3)
                cc=ccon(nsl)*(1.d0+xlam(nsl)*fct)
                ff=fcon(nsl)*(1.d0+xa2(nsl)*fct)
                y2a2=(w*(cc*a(2,nsl)+ff*(2.d0*a(1,nsl)-sfl3*a(3,nsl))) &
                                   *(scale*rhobar*vn**2/rn))**2 
!A00:
                Ann(1)= y2a2/2.d0/w/wn*rn**2
                Ann(2:4)=0.d0
             endif
          endif
       endif
!---------------------------
    case (1)
!*** radial modes ***
!---------------------------
       if(wdim.le.wgrav) kg=1 !on tient compte de la gravite
       if (l/=0) stop 'Les modes radiaux n''existe que pour l=0!'
       ls=2
       call rps(ls,ass)
       call rprop(ls,n,ass)
       select case(cond_limite)
       case(1)
!bord libre:
          det=a(2,n)/dsqrt(a(1,n)*a(1,n)+a(2,n)*a(2,n))
       case(2)
!bord rigide:
          det=a(1,n)/dsqrt(a(1,n)*a(1,n)+a(2,n)*a(2,n))
       end select
       if (present(check)) then
          check=a(1,n)*a(2,n)/dsqrt(a(1,n)*a(1,n)+a(2,n)*a(2,n))
       endif
       knt=kount-1

!
       if (present(continu).and.present(Ann)) then
          Ann(:)=0.0d0
          Ann(1)=a(2,n)/a(1,n)*rhobar*vn**2/rn
       endif
!------------
       if(ifeif.eq.0) return
!------------
       a(1,1)=0.d0
       a(2,1)=0.d0
       do  i=ls,n
          ff=fcon(i)*(1.d0+xlam(i)*fct)
          cc=ccon(i)*(1.d0+xa2(i)*fct)
          a(2,i)=(a(2,i)-2.d0*ff*a(1,i)/r(i))/cc
       enddo
       zi(1)=0.d0
       zi(2)=0.d0
       zi(3)=0.d0
       do i=ls,n
          im=i-1
          if (im/=0) then
             if(r(i).ne.r(im)) call gauslv(r(im),r(i),im,zi,3)
          endif
       enddo
       rnrm=1.d0/(w*dsqrt(zi(1)))
       cg=0.d0
       qinv=zi(2)/(wsq*zi(1))
       wray=dsqrt(zi(3)/zi(1))
       do  i=2,n
          do  j=1,2
             a(j,i)=a(j,i)*rnrm
          enddo
       enddo
!sortie pour l'operateur de couplage
       if (present(Ann)) then
          scale=1.0d0/sqrt(rhobar*rn**3)
          cc=ccon(n)*(1.d0+xlam(n)*fct)
          ff=fcon(n)*(1.d0+xa2(n)*fct)
          y2a2=(w*(cc*a(2,n)+2.d0*ff*a(1,n))*(scale*rhobar*vn**2/rn))**2
!A00:
          Ann(1)= y2a2/2.d0/w/wn*rn**2
          Ann(2:4)=0.d0
       endif
!---------------------------
    case (2,4)
!*** toroidal modes *** pour l/=0 et quand il n'y a pas que de l'eau
!---------------------------
       if (l/=0.and.noc/=n) then
          nb=nocp1
          n2=nsl
          ass(1)=1.d0
          ass(2)=0.d0
          if(jcom==4) then
             nb=2
             a(1,1)=0.d0
             a(2,1)=0.d0
             n2=nic
          endif
!
          q=0.d0
          ls=nb
          call startl(ls,n2,fmu,ls,q)
          ls=max(ls,nocp1)
          ls=max(ls,2)
!
          if(ls.ne.nocp1) call tps(ls,ass)
          call tprop(ls,n2,ass)
!
          select case(cond_limite)
          case(1)
!bord libre:
             det=a(2,n2)/dsqrt(a(1,n2)*a(1,n2)+a(2,n2)*a(2,n2))
!test
!             det=(a(2,n2)+1.771918364456304d-002*(fl3-2.)*a(1,n2)/200.)/dsqrt(a(1,n2)*a(1,n2)+a(2,n2)*a(2,n2))
          case(2)
!bord rigide:
             det=a(1,n2)/dsqrt(a(1,n2)*a(1,n2)+a(2,n2)*a(2,n2))
          end select
          if (present(check)) then
             check=a(1,n2)*a(2,n2)/dsqrt(a(1,n2)*a(1,n2)+a(2,n2)*a(2,n2))
          endif
!
          if (present(continu).and.present(Ann)) then
             Ann(:)=0.0d0
             Ann(1)=a(2,n)/a(1,n)/2.0d0*rhobar*vn**2/rn
          endif
!
          if(ifeif/=0) then
             do  i=ls,n2
                a(2,i)=a(1,i)/r(i)+a(2,i)/(lcon(i)*(1.d0+qshear(i)*fct))
             enddo
             if(ls/=nb) then
                ls1=ls-1
                do  i=nb,ls1
                   a(1,i)=0.d0
                   a(2,i)=0.d0
                enddo
             endif
             do  i=1,4
                zi(i)=0.d0
             enddo
             do i=ls,n2
                im=i-1
                if (im/=0) then
                   if(r(i).ne.r(im)) call gauslv(r(im),r(i),im,zi,4)
                endif
             enddo
             rnrm=1.d0/(w*dsqrt(zi(1)))
!test
!if (l==30) then
!print*,'rnrm=',rnrm
!rnrm=9.512561442056486d-006
!endif
!fin test
             cg=(fl+0.5d0)*zi(2)/(w*zi(1))
             qinv=zi(3)/(wsq*zi(1))
             wray=dsqrt(zi(4)/zi(1))
             do  i=ls,n2
                do  j=1,2
                   a(j,i)=a(j,i)*rnrm
                enddo
             enddo
!sortie pour l'operateur de couplage
             if (present(Ann)) then
                scale=1.0d0/sqrt(rhobar*rn**3)
                ll=lcon(n2)*(1.d0+qshear(n2)*fct)
                y2a2=(w*ll*(a(2,n2)-a(1,n2))*(scale*rhobar*vn**2/rn))**2
!A5:
                Ann(1)= y2a2/4.d0/w/wn*rn**2
!                print*,'Ann(5)= ', Ann(1)
             endif
          else
             if(knsw==1) then
!test
!                knt=kount-1               

                knt=kount-1
                if(jcom/=4.and.l/=1)then
                   irem=mod(knt,2)
                   if(.not.((irem.eq.0.and.det.lt.0.d0).or. &
                        (irem.ne.0.and.det.gt.0.d0))) knt=knt+1
                endif
             endif
          endif
       else
!l=0
          knt=-2
          if (present(Ann))ann(1)=0.d0
       endif
!----------------------------------------------------------------------
    end select
!----------------------------------------------------------------------
!----------------------------------------------------------------------
  end subroutine detqn
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  subroutine sprpmn(jf,jl,f,h,nvesm,iexp)
! propage les mineurs dans une region solide de jf a jl
!----------------------------------------------------------------------
    use shanks
    use bits
    use eifx
    use param_modele
!
    implicit none
    integer, intent(in) :: jf,jl,nvesm
    doubleprecision, dimension(nvesm), intent(inout) :: f
    integer, intent(inout) :: iexp
    doubleprecision :: h(nvesm,1) ! a revoir
    
    doubleprecision ::  econst=1048576.d0
!
    doubleprecision, dimension(14) :: s,fp
    doubleprecision, dimension(6) :: rne
    doubleprecision:: x,y,qff,qll,qaa,zs,xi,vpsq,vssq,alfsq,betasq,delsq &
                     ,fksq,qt,q,del,dxs,size,t1,t2,t3,t4,z
    integer        :: maxo1,jud,i,iq,ni,ierr_count_yann,jj,j
    logical :: flag_first2
!
    maxo1=maxo-1 ! rapport a la precision d'integration (??)
    if(jl.lt.jf) then
       jud=-1 !jud : sens dans le quel on parcour le couche (increment quoi)
    else
       jud=1
    endif
!
    y=r(jf)
    x=y
    i=jf
!------------------------------------- 
    do while ( i /= jl + jud )
!on integre de jf a jl
!------------------------------------- 
!
       if(abs((y-x)/x)>1.0d-15) then
          iq=min0(i,i-jud)
          qff=1.d0+xlam(iq)*fct
          qll=1.d0+qshear(iq)*fct
          qaa=1.d0+xa2(iq)*fct
          zs=dmin1(x,y)
          xi=g(i)/y
          vpsq=(flam(i)+2.d0*fmu(i))/rho(i)
          vssq=fmu(i)/rho(i)
          alfsq=(wsq+4.d0*rho(i)+xi)/vpsq
          betasq=wsq/vssq
          delsq=dsqrt((betasq-alfsq)**2+4.d0*fl3*xi*xi/(vssq*vpsq))
          fksq=.5d0*(alfsq+betasq+delsq)-fl3/(x*x)
          qt=dsqrt(dabs(fksq))+dsqrt(dabs(fksq-delsq))+2.d0/zs
          q=(qt+float(kg)*sfl3/x)/stepf
          del=float(jud)*step(maxo)/q
          dxs=0.d0
          !
          flag_first2=.true.
!
          do while ( abs((y-r(i))/y) > 1.0d-15 .or. flag_first2 )
             flag_first2=.false.
             y=x+del
             if(float(jud)*(y-r(i)).gt.0.d0) y=r(i)
             dx=y-x
!si le pas a change, on calcule les coefs pour l'integration rks
!  (pour runge-kutta-shanks) (common/shanks/)
             if(dx.ne.dxs) call baylis(q,maxo1) !calcul aussi in
             dxs=dx
             do j=1,nvesm
                s(j)=f(j)
             enddo
             do  ni=1,in !in est dans le module shanks
                z=x+c(ni)
!calcule des drive des mineurs dans un solide
                call derms(iq,z,f,h(1,ni),0,qff,qll,qaa,nvesm)
!                  call derms(iq,z,f,h((ni-1)*nvesm+1),0,qff,qll,qaa
!     >                        ,nevsm)
!produit scalair avec les coefs rks
                call rkdot(f,s,h,nvesm,ni)
             enddo
             if (knsw == 1) then
! si on doit compter les modes
!calcule des drive des mineurs dans un solide
                call derms(iq,y,f,fp,1,qff,qll,qaa,nvesm)
!compte les modes:
                call zknt(s,h,f,fp,x,y,1,nvesm)
             endif
             x=y
          enddo
!fait une de manipe (de normalisation (???))
       endif
!
       size=dabs(f(1))
       do  j=2,nvesm
          size=dmax1(size,dabs(f(j)))
       enddo
       ierr_count_yann=0
       do while(size >= 1024.d0)
          do  j=1,nvesm
             f(j)=f(j)/econst 
          enddo
          size=size/econst
!econst vaut 2^20
! a chaque fois qu'on divise pas econst on incremente iexp de 20
          iexp=iexp+20
          ierr_count_yann=ierr_count_yann+1
          if (ierr_count_yann.gt.10000) stop 'Ooopss dans sprpmn'
       enddo
!--------------------------------------------
       if ( iback /=0 ) then
!si on doit sortir les modes on integre de haut en bas
!--------------------------------------------
          inorm(i)=inorm(i)+iexp
          if (kg/=0) then
!gravite, exctraction des fctps a prtir des mineurs
             t1=f(4)+f(8)
             t2=t1+f(4)
             t1=t1+f(8)
             t3=f(8)-f(4)
             rne(1)=ar(6,i)*f(10)-ar(14,i)*f(9)+ar(13,i)*t3     &
                  -ar(1,i)*f(7)-ar(7,i)*f(6)+ar(8,i)*f(5)       &
                  +ar(12,i)*f(3)-ar(2,i)*f(2)+ar(3,i)*f(1)
!
             rne(2)=ar(6,i)*f(13)+ar(14,i)*t2+ar(13,i)*f(12)    &
                  -ar(1,i)*f(11)-ar(9,i)*f(6)-ar(7,i)*f(5)      &
                  +ar(11,i)*f(3)-ar(4,i)*f(2)-ar(2,i)*f(1)
!
             rne(3)=ar(6,i)*f(14)-ar(7,i)*t1-ar(8,i)*f(12)      &
                  +ar(13,i)*f(11)-ar(9,i)*f(9)+ar(14,i)*f(7)    &
                  +ar(10,i)*f(3)+ar(11,i)*f(2)+ar(12,i)*f(1)
!
             rne(4)=ar(14,i)*f(14)+ar(7,i)*f(13)+ar(12,i)*f(12) &
                  -ar(2,i)*f(11)-ar(9,i)*f(10)-ar(11,i)*t3      &
                  +ar(4,i)*f(7)+ar(10,i)*f(5)+ar(5,i)*f(1)
!
             rne(5)=ar(13,i)*f(14)+ar(8,i)*f(13)-ar(12,i)*t2    &
                  -ar(3,i)*f(11)+ar(7,i)*f(10)-ar(11,i)*f(9)    &
                  -ar(2,i)*f(7)+ar(10,i)*f(6)+ar(5,i)*f(2)
!
             rne(6)=ar(1,i)*f(14)+ar(13,i)*f(13)-ar(2,i)*t1     &
                  -ar(3,i)*f(12)+ar(14,i)*f(10)-ar(4,i)*f(9)    &
                  -ar(11,i)*f(6)-ar(12,i)*f(5)+ar(5,i)*f(3)
          else
!pas de gravite, exctraction des fctps a prtir des mineurs
             rne(1)=-ar(1,i)*f(3)+ar(2,i)*f(2)-ar(3,i)*f(1)
             rne(2)=-ar(1,i)*f(4)+ar(4,i)*f(2)+ar(2,i)*f(1)
             rne(3)=-ar(2,i)*f(4)+ar(4,i)*f(3)-ar(5,i)*f(1)
             rne(4)=-ar(3,i)*f(4)-ar(2,i)*f(3)-ar(5,i)*f(2)
          endif
! on replace les fctp dans le buffer
          do  jj=1,nvesm
             ar(jj,i)=rne(jj)
          enddo
       else
! pas de sorties de modes,  on place le mineurs dans le buffer
          inorm(i)=iexp
          do  j=1,nvesm
             ar(j,i)=f(j)
          enddo
       endif
       i=i+jud
       x=y
       if (i/= jl+jud) y=r(i)
!------------------------------------- 
! fin do while qui dirige l'integration
    enddo
!------------------------------------- 
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
  end subroutine sprpmn
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
  subroutine fprpmn(jf,jl,f,h,nvefm,iexp)
! propage les mineurs dans une region liquide de jf a jl
!----------------------------------------------------------------------
    use shanks
    use bits
    use eifx 
    use param_modele
!
    implicit none
!*** propagate the minor vector in a fluid region from level jf to jl ***
    integer, intent(in) :: jf,jl,nvefm
    doubleprecision, dimension(nvefm), intent(inout) :: f
    integer, intent(inout) :: iexp
    doubleprecision :: h(nvefm,1) ! a revoir
!
    doubleprecision :: econst=1048576.d0
!
    doubleprecision, dimension(5) :: s,fp
    doubleprecision :: x,y,qff,zs,xi,alfsq,q,del,dxs,z,size,rne2,rne3
    integer   :: maxo1,jud,i,j,iq,ni,ierr_count_yann
    logical :: flag_first

!

    if(nvefm/=1) then
!
       maxo1=maxo-1           ! rapport a la precision d'integration (??)
!
       if(jl.lt.jf) then
          jud=-1              !jud : sens dans le quel on parcour le couche (increment quoi)
       else
          jud=1
       endif
!
       y=r(jf)
       x=y
       i=jf
!
!------------------------------------- 
       do while ( i /= jl + jud )
!on integre de jf a jl
!------------------------------------- 
          if(abs((y-x)/x)>1.0d-15) then
             iq=min0(i,i-jud)
             qff=1.d0+xlam(iq)*fct
             zs=dmin1(x,y)
             xi=g(i)/y
             alfsq=(wsq+4.d0*rho(i)+xi-fl3*xi*xi/wsq)*rho(i)/flam(i)
             q=(dsqrt(dabs(alfsq-fl3/(x*x)))+1.d0/zs+float(kg)*sfl3/x)/stepf
             del=float(jud)*step(maxo)/q
             dxs=0.d0
!
             flag_first=.true.
             do while ( abs((y-r(i))/y) > 1.0d-15 .or. flag_first )
                flag_first=.false.
                y=x+del
                if(float(jud)*(y-r(i)).gt.0.d0) y=r(i)
                dx=y-x
!si le pas a change, on calcule les coefs pour l'integration rks
!  (pour runge-kutta-shanks) (common/shanks/)
                if(dx.ne.dxs) call baylis(q,maxo1)
                dxs=dx
                do  j=1,nvefm
                   s(j)=f(j)
                enddo
                do  ni=1,in
                   z=x+c(ni)
!calcule des drive des mineurs dans un fluide
                   call dermf(iq,z,f,h(1,ni),0,qff,nvefm)
                   call rkdot(f,s,h,nvefm,ni)
                enddo
                if(knsw == 1)then
! si on doit compter les modes
!calcule des drive des mineurs dans un fluide
                   call dermf(iq,y,f,fp,1,qff,nvefm)
!compte les modes:
                   call zknt(s,h,f,fp,x,y,0,nvefm)
                endif
                x=y
             enddo
!
          endif
          size=dabs(f(1))
          do  j=2,nvefm
             size=dmax1(size,dabs(f(j)))
          enddo
          ierr_count_yann=0
          do while(size >= 1024.d0)
             do  j=1,nvefm
                f(j)=f(j)/econst 
             enddo
             size=size/econst
!econst vaut 2^20
! a chaque fois qu'on divise pas econst on incremente iexp de 20
             iexp=iexp+20
             ierr_count_yann=ierr_count_yann+1
             if (ierr_count_yann.gt.10000) stop 'Ooopss dans fprpmn'
          enddo
!
          if(iback/=0) then
             inorm(i)=inorm(i)+iexp
             rne2   =-ar(1,i)*f(4)+ar(4,i)*f(2)+ar(2,i)*f(1)
             ar(1,i)=-ar(1,i)*f(3)+ar(2,i)*f(2)-ar(3,i)*f(1)
             rne3   =-ar(2,i)*f(4)+ar(4,i)*f(3)-ar(5,i)*f(1)
             ar(4,i)=-ar(3,i)*f(4)-ar(2,i)*f(3)-ar(5,i)*f(2)
             ar(2,i)=rne2
             ar(3,i)=rne3
          else
             inorm(i)=iexp
             do  j=1,nvefm
                ar(j,i)=f(j)
             enddo
          endif
          i=i+jud
          x=y
          if (i/= jl+jud)  y=r(i)
!------------------------------------- 
! fin do while qui dirige l'integration
       enddo
!-------------------------------------        
    else
       do  i=jl,jf
          inorm(i)=inorm(i)+iexp
          do  j=1,2
             ar(j,i)=ar(j,i)*f(1)
          enddo
       enddo
    endif
!----------------------------------------------------------------------
  end subroutine fprpmn
!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine eifout(lsmin)
!massages spheroidal mode eigenfunctions before output ***
!lcmin : couche minimum (apres ls, les  fctps ne sont plus nulles)
! avant :
!    y1(:)=a(1,:)
!    y3(:)=a(2,:) 
!    y2(:)=a(3,:) 
!    y4(:)=a(4,:)  
!apres:
!    u (:)=a(1,:)
!    du(:)=a(2,:) 
!    v (:)=a(3,:) 
!    dv(:)=a(4,:)  
!+ normalisation
!
!----------------------------------------------------------------------
!++++++++++++++++++++++++++
    use bits
    use eifx, a=>ar
    use rindex
    use param_modele
!++++++++++++++++++++++++++
    implicit none
    integer, intent(in) :: lsmin
!
    doubleprecision ::  ll,ff,zr,sfl3z,d,v,p,iexp,ffi,rnorm,al
    doubleprecision, dimension(4) ::  zi
    integer :: i,j,i1,i2,iq,lsm1,ip,imax
!***************************
!----------------dans la graine  (si necessaire en fct de lsmin)
    i1=min0(nic,max0(2,lsmin))
    i2=nic
!
    if(i1/=i2) then
       do iq=i1,i2
          ff=fcon(iq)*(1.d0+xlam(iq)*fct)
          ll=lcon(iq)*(1.d0+qshear(iq)*fct)
          zr=1.d0/r(iq)
          sfl3z=sfl3*zr
          d=1.d0/(ccon(iq)*(1.d0+xa2(iq)*fct))
          v=a(2,iq)
          if(kg==0) then
!pas de gravite:
             a(2,iq)=(zr-2.d0*ff*d*zr)*a(1,iq)+sfl3z*ff*d*v+d*a(3,iq)
             a(4,iq)=-sfl3z*a(1,iq)+(zr+zr)*v+a(4,iq)/ll
             a(5,iq)=0.d0
             a(6,iq)=0.d0
          else
!gravite :
             a(2,iq)=(zr-2.d0*ff*d*zr)*a(1,iq)+sfl3z*ff*d*v+d*a(4,iq)
             a(4,iq)=-sfl3z*a(1,iq)+(zr+zr)*v+a(5,iq)/ll
             a(5,iq)=a(3,iq)
             a(6,iq)=4.d0*(a(6,iq)-rho(iq)*a(1,iq))-fl*zr*a(5,iq)
          endif
          a(3,iq)=v
       enddo
    endif
!
!----------------dans le noyau liquide (si necessaire en fct de lsmin)
    i1=min0(nsl,max0(lsmin,nocp1))
    i2=nsl
!
    if(i1/=i2) then
       do iq=i1,i2
          ff=fcon(iq)*(1.d0+xlam(iq)*fct)
          ll=lcon(iq)*(1.d0+qshear(iq)*fct)
          zr=1.d0/r(iq)
          sfl3z=sfl3*zr
          d=1.d0/(ccon(iq)*(1.d0+xa2(iq)*fct))
          v=a(2,iq)
          if(kg==0) then
!pas de gravite:
             a(2,iq)=(zr-2.d0*ff*d*zr)*a(1,iq)+sfl3z*ff*d*v+d*a(3,iq)
             a(4,iq)=-sfl3z*a(1,iq)+(zr+zr)*v+a(4,iq)/ll
             a(5,iq)=0.d0
             a(6,iq)=0.d0
          else
!gravite :
             a(2,iq)=(zr-2.d0*ff*d*zr)*a(1,iq)+sfl3z*ff*d*v+d*a(4,iq)
             a(4,iq)=-sfl3z*a(1,iq)+(zr+zr)*v+a(5,iq)/ll
             a(5,iq)=a(3,iq)
             a(6,iq)=4.d0*(a(6,iq)-rho(iq)*a(1,iq))-fl*zr*a(5,iq)
          endif
          a(3,iq)=v
       enddo
    endif
!
!!----------------dans le manteau (si necessaire )
    i1=min0(noc,max0(lsmin,nicp1))
    i2=noc
    if (i1 /= i2) then
       do  iq=i1,i2
          zr=1.d0/r(iq)
          sfl3z=sfl3*zr
          ffi=1.d0/(flam(iq)*(1.d0+xlam(iq)*fct))
          if(kg==0) then
!pas de gravite:
             p=a(2,iq)
             a(5,iq)=0.d0
             a(6,iq)=0.d0
          else
!gravite :
             p=a(3,iq)
             a(5,iq)=a(2,iq)
             a(6,iq)=4.d0*(a(4,iq)-rho(iq)*a(1,iq))-fl*zr*a(5,iq)
          endif
          a(3,iq)=sfl3z*(g(iq)*a(1,iq)-p/rho(iq)+a(5,iq))/wsq
          a(2,iq)=sfl3z*a(3,iq)-a(1,iq)*zr+p*ffi
          a(4,iq)=sfl3z*(a(1,iq)+p*(qro(1,iq)/(rho(iq)**2)+g(iq)*ffi)/wsq)
       enddo
    endif
    if(n/=nsl.and.i2/=n) then
!!----------------dans l'ocean (si necessaire )
       i1=nslp1
       i2=n
       do  iq=i1,i2
          zr=1.d0/r(iq)
          sfl3z=sfl3*zr
          ffi=1.d0/(flam(iq)*(1.d0+xlam(iq)*fct))
          if(kg==0) then
!pas de gravite:
             p=a(2,iq)
             a(5,iq)=0.d0
             a(6,iq)=0.d0
          else
!gravite :
             p=a(3,iq)
             a(5,iq)=a(2,iq)
             a(6,iq)=4.d0*(a(4,iq)-rho(iq)*a(1,iq))-fl*zr*a(5,iq)
          endif
          a(3,iq)=sfl3z*(g(iq)*a(1,iq)-p/rho(iq)+a(5,iq))/wsq
          a(2,iq)=sfl3z*a(3,iq)-a(1,iq)*zr+p*ffi
          a(4,iq)=sfl3z*(a(1,iq)+p*(qro(1,iq)/(rho(iq)**2)+g(iq)*ffi)/wsq)
       enddo
    endif
!
! renormalisation (avec les puissance de 2 contenues dans inorm)
    imax=0
    do  iq=lsmin,n
       imax=max0(inorm(iq),imax)
    enddo
    do  iq=lsmin,n
       iexp=inorm(iq)-imax
       al=0.d0
       if(iexp.ge.-80) al=2.d0**iexp
       do  j=1,6
          a(j,iq)=a(j,iq)*al
       enddo
    enddo
    lsm1=max0(1,lsmin-1)
    do  i=1,lsm1
       do  j=1,6
          a(j,i)=0.d0
       enddo
    enddo
!
    if(l<=1.and.lsmin<=2) then
       a(2,1)=1.5d0*a(1,2)/r(2)-.5d0*a(2,2)
       a(4,1)=1.5d0*a(3,2)/r(2)-.5d0*a(4,2)
    endif
    do j=1,4
       zi(j)=0.d0
    enddo
    i1=max0(lsmin,2)
    do  iq=i1,n
       ip=iq-1
!integration de gauss-legendre qu 5eme ordre de I1,I2,I3,I4 (voir intgds (c''est vraiment
! mal foutu):
       if(r(iq).ne.r(ip)) call gauslv(r(ip),r(iq),ip,zi,4)
    enddo
!I1=zi(1) ?
    cg=zi(2)/(w*zi(1)) 
    wray=dsqrt(2.d0*zi(4)/zi(1))
    qinv=2.d0*zi(3)/(wsq*zi(1))
    rnorm=1.d0/(w*dsqrt(zi(1)))
    do  iq=i1,n
       zr=1.d0/r(iq)
       a(1,iq)=a(1,iq)*zr
       a(2,iq)=(a(2,iq)-a(1,iq))*zr
       a(3,iq)=a(3,iq)*zr
       a(4,iq)=(a(4,iq)-a(3,iq))*zr
       a(5,iq)=a(5,iq)*zr
       a(6,iq)=(a(6,iq)-a(5,iq))*zr
       a(1,iq)=a(1,iq)*rnorm
       a(2,iq)=a(2,iq)*rnorm
       a(3,iq)=a(3,iq)*rnorm
       a(4,iq)=a(4,iq)*rnorm
       a(5,iq)=a(5,iq)*rnorm
       a(6,iq)=a(6,iq)*rnorm
    enddo
    if(lsmin<=2.and.l<=2) then 
       if(l/=2) then
          a(1,1)=a(1,2)-.5d0*a(2,2)*r(2)
          a(2,1)=0.d0
          a(3,1)=a(3,2)-.5d0*a(4,2)*r(2)
          a(4,1)=0.d0
          a(6,1)=1.5d0*a(5,2)/r(2)-.5d0*a(6,2)
       else
          a(2,1)=1.5d0*a(1,2)/r(2)-.5d0*a(2,2)
          a(4,1)=1.5d0*a(3,2)/r(2)-.5d0*a(4,2)
       endif
    endif
      
!----------------------------------------------------------------------
  end subroutine eifout
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine sdepth(wdim,ls)
! pour une frequnce de depart donnee, evalue le numero de couche a partir
! duquel il faut integer
!----------------------------------------------------------------------
    use bits
    use rindex
    use param_modele
    implicit none
!*** finds starting level,ls, for a given l and w ***
    doubleprecision, intent(in) :: wdim
    integer, intent(out) :: ls
!
    doubleprecision :: aw=-2.d-3,bw=2.25d-3,dw=1.28d-3
!
    doubleprecision :: q,wsoc,wsic
!
    q=0.d0
    w=wdim/wn !est dans le module bits
    wsoc=aw+dw*fl
    if(wdim <= wsoc) then
       call startl(nocp1,nsl,fmu,ls,q)
       if(ls.eq.nsl) ls=ls-1
       if(ls.gt.nocp1) return !hummm,, faudrait revoir cette facon de coder
    endif
    wsic=aw+bw*fl
    if(wdim <= wsic)  then
       call startl(nicp1,noc,flam,ls,q)
       if(ls.eq.noc) ls=ls-1
       if(ls.gt.nicp1) return
    endif
    call startl(2,nic,fmu,ls,q)
    if(ls.eq.nic) ls=ls-1
    return
!----------------------------------------------------------------------
  end subroutine sdepth
!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine fpsm(ls,nvefm,ass)
! calcule les solution des depart dans une region fluide (avec des fcts
! de Bessel spheriques
!----------------------------------------------------------------------
    use bits
    use param_modele
    use yannos_flag
    implicit none
!
    integer, intent(in) :: ls,nvefm
    doubleprecision, dimension(nvefm), intent(out) :: ass    
!
    doubleprecision :: x,fla,vpsq,xi,qsq,zsq,u,c1,c2,sum,fp
    integer :: i
!
    x=r(ls)
    fla=flam(ls)*(1.d0+xlam(ls)*fct)
    vpsq=fla/rho(ls)
    xi=g(ls)/x    
    qsq=(wsq+float(kg)*4.d0*rho(ls)+xi-fl3*xi*xi/wsq)/vpsq
    zsq=qsq*x*x
    call bfs(l,zsq,eps,fp)
    if(kg/=0) then
!gravite
       u=(fl-fp)/qsq
       c1=fl*g(ls)-wsq*x
       c2=fl2*c1*0.25d0/x-rho(ls)*fl
       ass(1)=-x*fl*vpsq-c1*u
       ass(2)=-x*fl*fla
       ass(3)=-fl*fl2*vpsq*0.25d0-u*c2
       ass(4)=x*fla*c1
       ass(5)=-x*fla*c2
    else
!pas gravite
       ass(1)=-(fl3*xi/wsq+fp)/qsq
       ass(2)=x*fla
    endif
!
    sum=ass(1)*ass(1)
    do  i=2,nvefm
       sum=sum+ass(i)*ass(i)
    enddo
    sum=1.d0/dsqrt(sum)
    if(ass(nvefm).lt.0.d0) sum=-sum
!
    do  i=1,nvefm
       ass(i)=ass(i)*sum
    enddo
!----------------------------------------------------------------------
  end subroutine fpsm
!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine spsm(ls,nvesm,ass)
! calcule les solution des depart dans une region solide (avec des fcts
! de Bessel spheriques
!
! nvesm: nombre de mineurs dans le solide (5, sans gravite
                    ! ou 14 avec gravite)
!----------------------------------------------------------------------
    use bits
    use param_modele
    use rindex
    use yannos_flag
    implicit none
!
    integer, intent(in) :: ls,nvesm
    doubleprecision, dimension(nvesm), intent(out) :: ass
!   
    doubleprecision, dimension(6,2) :: a
    doubleprecision, dimension(15)  :: e
    doubleprecision :: x,ro,fu,flu,vssq,vpsq,zeta,zi,alfsq,betasq &
                      ,delsq,fksq,qsq,zsq,b,c0,c1,c2,c3,sum,xi,fp,vs,vp,om
    integer :: i,j,k,jj,kk,ll,i1
!    
    x=r(ls)
    ro=rho(ls)
    fu=fmu(ls)*(1.d0+qshear(ls)*fct)
    flu=flam(ls)*(1.d0+xlam(ls)*fct)+2.d0*fu
    vssq=fu/ro  !vs^2
    vpsq=flu/ro !vp^2
    zeta=4.d0*ro
    xi=g(ls)/x
    alfsq=(wsq+float(kg)*zeta+xi)/vpsq
    betasq=wsq/vssq
    delsq=dsqrt((betasq-alfsq)**2+4.d0*fl3*xi*xi/(vpsq*vssq))
    fksq=.5d0*(alfsq+betasq+delsq)
!qsq est k^2 de la p244 Takeuchi & Saito
    qsq=fksq-delsq
    zsq=qsq*x*x
    if (cancel_gravity) then
       vp=sqrt(vpsq)
       vs=sqrt(vssq)
       om =sqrt(wsq )
       ass(:)=0.d0
       call spsm_nograve(om,l,r(ls),vs,vp,ro,ass)
    else
!b=-1/f de la p244 Takeuchi & Saito
       b=xi/(vssq*(betasq-qsq))
       do k=1,2
!solution pas 245 Takeuchi & Saito
          call bfs(l,zsq,eps,fp)
          a(1,k)=fl3*b+fp
          a(2,k)=1.d0+b+b*fp
          a(3,k)=-zsq
          a(4,k)=b*a(3,k)
          a(5,k)=1.d0
          a(6,k)=fl1-fl3*b
          if (k==1) then
             zsq=fksq*x*x
             b=-flu/(fu*fl3*b)
          endif
       enddo
!
       jj=3+2*kg
       kk=jj+1
       ll=0
!
       do  i=1,jj
          i1=i+1
          do  j=i1,kk
             ll=ll+1
             e(ll)=a(i,1)*a(j,2)-a(j,1)*a(i,2)
          enddo
       enddo
       if(kg==0) then
!pas de gravite
          ass(1)=x*x*e(1)
          ass(2)=fu*x*sfl3*(2.d0*e(1)-e(5))
          ass(3)=fu*x*(e(3)-2.d0*e(1))
          ass(4)=x*(flu*e(4)+4.d0*fu*e(1))
          ass(5)=fu*(flu*(e(6)+2.d0*e(4))+4.d0*fu*(fl3*(e(5)-e(1)) &
               -e(3)+2.d0*e(1)))
       else
!gravite
          c0=wsq-xi*fl
          c1=ro*fl+0.25d0*fl2*c0
          c2=2.d0*fu/x
          c3=c2*(fl-1.d0)
          ass(6)=x*x*(c0*e(1)-zeta*(fl*e(8)-e(4)))
          ass(14)=flu*(fl*e(6)-e(2))
          ass(13)=fu*sfl3*(fl*e(7)-e(3))
          ass(1)=x*(c1*e(1)-ro*(fl*e(9)-e(5)))
          ass(7)=x*flu*(c0*e(2)-zeta*fl*e(11))/sfl3+c2*sfl3*ass(6)
          ass(8)=x*fu*(c0*e(3)-zeta*fl*e(13))-c2*ass(6)
          ass(12)=(flu*fl*e(10)+2.d0*(ass(14)+sfl3*ass(13)))*fu/x
          ass(2)=flu*(c1*e(2)-ro*fl*e(12))/sfl3+c2*sfl3*ass(1)
          ass(3)=fu*(c1*e(3)-ro*fl*e(14))-c2*ass(1)
          ass(9)=(x*c0*ass(14)+sfl3*ass(7)-c3*fl*ass(6))/fl
          ass(11)=(sfl3*ass(12)+c3*(sfl3*ass(14)-fl*ass(13)))/fl
          ass(4)=(c1*ass(14)+sfl3*ass(2)-c3*fl*ass(1))/fl
          ass(10)=(x*c0*ass(11)-c3*(sfl3*ass(9)+fl*ass(7)))/sfl3
          ass(5)=(c1*ass(11)-c3*(sfl3*ass(4)+fl*ass(2)))/sfl3
       endif
    endif
    sum=ass(1)*ass(1)
    do  i=2,nvesm
       sum=sum+ass(i)*ass(i)
    enddo
    sum=1.d0/dsqrt(sum)
    if (cond_limite==1) then
       if(ass(5).lt.0.d0) sum=-sum
    else
       if(ass(1).lt.0.d0) sum=-sum
    endif
    do  i=1,nvesm
       ass(i)=ass(i)*sum
    enddo
!
!----------------------------------------------------------------------
  end subroutine spsm
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine tprop(jf,jl,f)
! propage les mineurs dans une region solide de jf a jl pour les modes
! toroidaux
!----------------------------------------------------------------------
    use shanks
    use bits
    use eifx,a=>ar
    use param_modele
    implicit none
!
    integer, intent(in) :: jf,jl
    doubleprecision, dimension(2), intent(inout) ::  f
    
    doubleprecision, dimension(2,10) :: h
    doubleprecision, dimension(2)    :: s
    doubleprecision :: nn,ll,fl3m2,y,vx,vy,x,qll,qx,qy,q,del,dxs,z,ro,fp,t
    integer :: maxo1,i,iq,ni
!
    fl3m2=fl3-2.d0
    maxo1=maxo-1
    y=r(jf)
    vy=fmu(jf)/rho(jf)
    i=jf
    a(1,i)=f(1)
    a(2,i)=f(2)
!--------------------
    do while (i/=jl)
!on integre de jf a jl
!--------------------
       iq=i
       i=i+1
       x=y
       y=r(i)
       if(dabs(y-x) > 1.d-18 )  then
          qll=1.d0+qshear(iq)*fct
          vx=vy
          vy=fmu(i)/rho(i)
          qx=1.d0/x+dsqrt(dabs(wsq/(vx)-fl3/(x*x)))
          qy=1.d0/y+dsqrt(dabs(wsq/(vy)-fl3/(y*y)))
          q=dmax1(qx,qy)
          del=step(maxo)/q
          dxs=0.d0
!--------------------
          do while (dabs(x-r(i))>1.d-18)
!--------------------
             y=x+del
             if(y.gt.r(i)) y=r(i)
             dx=y-x
!si le pas a change, on calcule les coefs pour l'integration rks
!  (pour runge-kutta-shanks) (common/shanks/)
             if(dx.ne.dxs) call baylis(q,maxo1)
             dxs=dx
             s(1)=f(1)
             s(2)=f(2)
             do ni=1,in
                z=x+c(ni)
                t=z-r(iq)
                ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
                ll=(lcon(iq)+t*(lspl(1,iq)+t*(lspl(2,iq)           &
                     +t*lspl(3,iq))))*qll
                nn=ll
                if(ifanis.ne.0) nn=(ncon(iq)+                      &
                     t*(nspl(1,iq)+t*(nspl(2,iq)+t*nspl(3,iq))))*qll
                z=1.d0/z
                h(1,ni)=z*f(1)+f(2)/ll
                h(2,ni)=(nn*fl3m2*z*z-ro*wsq)*f(1)-3.d0*z*f(2)
                call rkdot(f,s,h,2,ni)
             enddo
             if(knsw==1) then
!  on doit compter les modes:
                select case (cond_limite)
                case(1)
!             bord libre:
                   fp=(nn*fl3m2*z*z-ro*wsq)*f(1)-3.d0*z*f(2)
                   call trknt(s(2),h(2,1),f(2),fp,x,y)
                case(2)
!             bord rigide:
                   fp=z*f(1)+f(2)/ll
                   call trknt(s(1),h(1,1),f(1),fp,x,y)
                end select
             endif
             x=y
!--------------------
          enddo
!--------------------
       endif
!
       a(1,i)=f(1)
       a(2,i)=f(2)
!--------------------
    enddo
!--------------------
!----------------------------------------------------------------------
  end subroutine tprop
!----------------------------------------------------------------------
!--------------------------------------------------------------------------
  subroutine zknt(s,sp,f,fp,x,y,ifsol,nvec)
!construit le mode count a partir d'un vecteur de mineur et de sa derivee
! ifsol est un flag:
!      0:liquide
!      1:solide
!--------------------------------------------------------------------------
    use bits
    implicit none
    integer, intent(in) :: nvec,ifsol
    doubleprecision, dimension(nvec), intent(in) :: s,sp,f,fp
    doubleprecision, intent(in) :: x,y
!
    integer :: ns,ns1,ns2,ift,i,j
    doubleprecision, dimension(4) ::xs,val
    doubleprecision :: y1,y2,y1p,y2p,t1,t2,t1p,t2p,h,a1,a2,a3,a22,a33 &
                      ,disc,tr1,tr2,fac,tes,rt1,rt,v,vp,add,b1,b2,b3  &
                      ,t
    logical :: go_on
!
    if(ifsol/=0.or.kg/=0) then
!gravite ou partie solide
!
       select case (cond_limite)
       case(1)
! pour bord libre
          y1=s(5)
          y2=f(5)
          y1p=sp(5)
          y2p=fp(5)
          t1=s(3)-s(4)
          t2=f(3)-f(4)
          t1p=sp(3)-sp(4)
          t2p=fp(3)-fp(4)
!pour bord rigide
       case(2)
          y1=s(1)
          y2=f(1)
          y1p=sp(1)
          y2p=fp(1)
!
          t1=-(s(3)-s(4))
          t2=-(f(3)-f(4))
          t1p=-(sp(3)-sp(4))
          t2p=-(fp(3)-fp(4))
       end select
    else
!pas de gravite et liquite
       select case (cond_limite)
       case(1)
! pour bord libre
          y1=s(2)
          y2=f(2)
          y1p=sp(2)
          y2p=fp(2)
          t1=s(1)
          t2=f(1)
          t1p=sp(1)
          t2p=fp(1)
       case(2)
!pour bord rigide
          y1=s(1)
          y2=f(1)
          y1p=sp(1)
          y2p=fp(1)
          t1=-s(2)
          t2=-f(2)
          t1p=-sp(2)
          t2p=-fp(2)
       end select
    endif
    h=y-x
    ns=0
    if(kount==0) then
      a1=y2-y1
      a2=0.d0
      a3=0.d0
      a22=0.d0
      a33=0.d0
   else
      a1=h*y1p
      a2=-h*(2.d0*y1p+y2p)+3.d0*(y2-y1)
      a3=h*(y1p+y2p)-2.d0*(y2-y1)
      a33=3.d0*a3
      a22=2.d0*a2
      if(a3==0.d0) then
         if(a2/=0.d0) then
            xs(2)=-a1/a22
            if(xs(2).ge.0.d0.and.xs(2).le.1.d0) ns=1
         endif
     else
         disc=a2*a2-a1*a33
          if(disc==0.0d0) then
             xs(2)=-a2/a33
             if(xs(2).ge.0.d0.and.xs(2).le.1.d0) ns=1
          else if( disc>0.d0) then
             disc=dsqrt(disc)
             tr1=(-a2+disc)/a33
             tr2=(-a2-disc)/a33
             if(dabs(a33)<=dabs(a1)) then
                fac=a1/a33
                tr1=fac/tr1
                tr2=fac/tr2
             endif
             if(tr1>=0.d0 .and. tr1<=1.d0) then
                xs(2)=tr1
                ns=1
             endif
             if(tr2>=0.d0 .and. tr2<=1.d0) then
                ns=ns+1
                xs(ns+1)=tr2
                if(ns >= 2 .and. tr2 < tr1) then
                   xs(2)=tr2
                   xs(3)=tr1
                endif
             endif
          else
             xs(2)=-a2/a33
             if(xs(2).ge.0.d0.and.xs(2).le.1.d0) ns=1
          endif
       endif
    endif
    val(1)=y1
    xs(1)=0.d0
    ns2=ns+2
    val(ns2)=y2
    xs(ns2)=1.d0
    if(ns/=0) then
       ns1=ns+1
       do j=2,ns1
          t=xs(j)
          val(j)=y1+t*(a1+t*(a2+t*a3))
       enddo
    endif
    ift=0!
    do  j=2,ns2
       if(val(j-1)*val(j)<=0.d0) then
          if(val(j-1)==0.d0) then
             tes=t1*a1
          else
             rt1=0.5d0*(xs(j-1)+xs(j))
             rt=rt1
             i=1
             go_on=.true.
             do while(i<=5.and.go_on)
                v=y1+rt*(a1+rt*(a2+rt*a3))
                vp=a1+rt*(a22+rt*a33)
                add=-v/vp
                rt=rt+add
                i=i+1
                if (dabs(add) < 1.d-5 ) then
                   go_on=.false.
                else if (dabs(rt-rt1) <= 0.5d0 ) then
                   go_on=.true.
                else
                   go_on=.false.
                   rt=rt1
                endif
             enddo
             if(ift==0) then
                if(kount==0) then
                   b1=t2-t1
                   b2=0.d0
                   b3=0.d0
                else
                   b1=h*t1p
                   b2=-h*(2.d0*t1p+t2p)+3.d0*(t2-t1)
                   b3=h*(t1p+t2p)-2.d0*(t2-t1)
                   ift=1
                endif
             endif
             tes=t1+rt*(b1+rt*(b2+rt*b3))
             vp=a1+rt*(a22+rt*a33)
             tes=tes*vp
          endif
          if(tes.lt.0.d0) kount=kount+1
          if(tes.gt.0.d0) kount=kount-1
       endif
    enddo
!--------------------------------------------------------------------------
  end subroutine zknt
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
  subroutine trknt(y1,y1p,y2,y2p,x,y)
!construit le mode count a partir d'un vecteur de mineur et de sa derivee
! pour des modes radiaux et toroidaux
!--------------------------------------------------------------------------
    use bits
    implicit none
!
    doubleprecision, intent(in)   :: y1,y1p,y2,y2p,x,y
!
    integer :: ns,ns1,ns2,j
    doubleprecision :: h,a1,a2,a3,a33,a22,disc,tr1,tr2,t,fac
    doubleprecision, dimension(2) :: xs
    doubleprecision, dimension(4) :: val
!
    ns=0
    if(kount/=0) then
       h=y-x
       a1=h*y1p
       a2=-h*(2.d0*y1p+y2p)+3.d0*(y2-y1)
       a3=h*(y1p+y2p)-2.d0*(y2-y1)
       a33=3.d0*a3
       a22=2.d0*a2
!       if(abs(a3) <= 1.d-15) then
       if(a3 == 0.d0) then
          disc=a2*a2-a1*a33
!          if (dabs(disc)<1.d-15) then
          if (disc<0.d0) then
             xs(1)=-a2/a33
             if(xs(1).ge.0.d0.and.xs(1).le.1.d0) ns=1
          else if (disc>0.d0) then
             disc=dsqrt(disc)
             tr1=(-a2+disc)/a33
             tr2=(-a2-disc)/a33
             if(dabs(a33)<=dabs(a1)) then
                fac=a1/a33
                tr1=fac/tr1
                tr2=fac/tr2
             endif
             if(tr1>=0.d0.and.tr1<=1.d0) then
                xs(1)=tr1
                ns=1
             endif
             if(tr2>=0.d0.and.tr2<=1.d0) then
                ns=ns+1
                xs(ns)=tr2
                if(ns>=2.and.tr2<tr1) then
                   xs(1)=tr2
                   xs(2)=tr1
                endif
             endif
          endif
!       else if(abs(a2)>1.d-15) then
       else if(a2/=0.d0) then
          xs(1)=-a1/a22
          if(xs(1).ge.0.d0.and.xs(1).le.1.d0) ns=1
       endif
       if(ns/=0) then
          ns1=ns+1
          do  j=2,ns1
             t=xs(j-1)
             val(j)=y1+t*(a1+t*(a2+t*a3))
          enddo
       endif
    endif
    val(1)=y1
    ns2=ns+2
    val(ns2)=y2
!
    do  j=2,ns2
       if(val(j-1)*val(j).le.0.d0) then
          kount=kount+1
       endif
    enddo
    if(val(1).eq.0.d0) kount=kount-1
!
!--------------------------------------------------------------------------
    end subroutine trknt
!--------------------------------------------------------------------------
!---------------------------------------------------------------------------
  subroutine derms(iq,z,f,fp,iknt,qff,qll,qaa,nevsm)
! Calcules les derivees des mineurs (fp) dans un solid
!---------------------------------------------------------------------------
    use bits
    use param_modele
    use yannos_flag
    implicit none
!
    integer, intent(in) :: iq, iknt,nevsm
    doubleprecision, intent(in) :: z,qff,qll,qaa
    doubleprecision, dimension(nevsm), intent(in ) :: f
    doubleprecision, dimension(nevsm), intent(out) :: fp
    
    doubleprecision ::  nn,ll,ro,gr,ff,cc,aa,zr,sfl3z,dmg,zdmg,t &
                        ,rogr 
    doubleprecision, save :: c11,c22,t11,t12,t21,t22,s22,s11,s12 &
                            ,b11,b33,t31,t33,s23,b44,b55,b32 &
                            ,b42,b52,b313,b414,b919,b66,b99,t4   &
                            ,t5,s13,b914,b22
!-----------------------------
    if (iknt/=0) then
! mode count
!-----------------------------
       if (kg==0) then
! pas de gravite
          fp(3)=s22*f(1)-2.d0*t12*f(2)+b33*f(3)+c11*f(5)
          fp(4)=-s11*f(1)+2.d0*t21*f(2)-b33*f(4)-c22*f(5)
          fp(5)=-2.d0*s12*f(2)+s11*f(3)-s22*f(4)-b11*f(5)            
       else
! gravite
          fp(3)=s22*f(1)+b32*f(2)+b33*f(3)+c11*f(5)+b313*f(13)
          fp(4)=-s11*f(1)+b42*f(2)+b44*f(4)-c22*f(5)+b414*f(14)
          fp(5)=b52*f(2)+s11*f(3)-s22*f(4)+b55*f(5)-b313*f(11)  &
               +b414*f(12)
       endif
!-----------------------------
    else
! pas de mode count
!-----------------------------
       t=z-r(iq)
!
       if(abs(t)<=1.0d-15) then
!
!pas besoin d'interpolation
!
          ro=rho(iq)
          gr=g(iq)
          ff=fcon(iq)*qff
          ll=lcon(iq)*qll
          nn=ncon(iq)*qll
          cc=ccon(iq)*qaa
          aa=acon(iq)*qaa
       else
!
!interpolation
!
          ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
          gr=g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq)))
          ff=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq))))*qff
          ll=(lcon(iq)+t*(lspl(1,iq)+t*(lspl(2,iq)+t*lspl(3,iq))))*qll
!
          if(ifanis==0) then
!isotrope
             nn=ll
             cc=ff+ll+ll
             aa=cc
          else
!anisotrope
             nn=(ncon(iq)+t*(nspl(1,iq)+t*(nspl(2,iq)+t*nspl(3,iq)))) &
                  *qll
             cc=(ccon(iq)+t*(cspl(1,iq)+t*(cspl(2,iq)+t*cspl(3,iq)))) &
                  *qaa
             aa=(acon(iq)+t*(aspl(1,iq)+t*(aspl(2,iq)+t*aspl(3,iq)))) &
                  *qaa
          endif
       endif
       zr=1.d0/z
       sfl3z=sfl3*zr
       rogr=ro*gr
       c11=1.d0/cc            !1/C                   =a12
       c22=1.d0/ll            !1/L                   =a34
       dmg=aa-nn-ff*ff*c11    !A-F^2/C-N
       zdmg=zr*dmg            !(A-F^2/C-N)/r
       t11=-2.d0*ff*zr*c11+zr !(1-2F/C)/r             =a11
       t12=sfl3z*ff*c11       ![l(l+1/2)]^(1/2)/r*F/C =a13
       t21=-sfl3z             !-[l(l+1/2)]^(1/2)/r    =a24
       t22=zr+zr              !2/r                    =a33
       s22=-ro*wsq            !-p*w^2
       s11=s22+4.d0*zr*(zdmg-rogr)        !-p*w^2+4/r^2((A-F^2/C-N)-p*gr*r)=a23+...
       s22=s22+zr*zr*(fl3*(dmg+nn)-nn-nn) !-p*w^2+ ....   =a43
       s12=sfl3z*(rogr-zdmg-zdmg)         !
!-------------------------------------------
       if(kg==0) then      
!pas de gravite
!
          if (.not.cancel_gravity) s11=s11+4.d0*ro*ro
!-------------------------------------------
          if(iback.eq.1) then
!intergration de bas en haut
             fp(1)=t22*f(1)-t21*f(2)-c22*f(3)
             fp(2)=-t12*f(1)+t11*f(2)-c11*f(4)
             fp(3)=-s22*f(1)+s12*f(2)-t22*f(3)+t12*f(4)
             fp(4)=s12*f(1)-s11*f(2)+t21*f(3)-t11*f(4)
!-------------------------------------------
          else
!intergration de haut en bas
             b11=t11+t22
             b33=t11-t22
             fp(1)=b11*f(1)+c22*f(3)-c11*f(4)
             fp(2)=s12*f(1)-t21*f(3)+t12*f(4)
             fp(3)=s22*f(1)-2.d0*t12*f(2)+b33*f(3)+c11*f(5)
             fp(4)=-s11*f(1)+2.d0*t21*f(2)-b33*f(4)-c22*f(5)
             fp(5)=-2.d0*s12*f(2)+s11*f(3)-s22*f(4)-b11*f(5)
          endif
!-------------------------------------------
       else
!gravite
          t31=-4.d0*ro
          t33=-fl*zr
          s13=-fl1*zr*ro
          s23=ro*sfl3z
!-------------------------------------------
          if(iback/=1) then
!intergration de haut en bas
             b11=t11+t22-t33
             b33=t11-t22-t33
             b44=t22-t11-t33
             b55=-t11-t22-t33
             b32=-t12-t12
             b42=t21+t21
             b52=-s12-s12
             b313=-s23-s23
             b414=s13+s13
             b914=t31+t31
             fp(1)=b11*f(1)+c22*f(3)-c11*f(4)
             fp(2)=s12*f(1)-t33*f(2)-t21*f(3)+t12*f(4)-s13*f(13)    &
                  -s23*f(14)
             fp(6)=4.d0*f(1)-b55*f(6)+c22*f(8)-c11*f(9)
             fp(7)=4.d0*f(2)+s12*f(6)+t33*f(7)-t21*f(8)+t12*f(9)    &
                  -t31*f(13)
             fp(8)=4.d0*f(3)+s22*f(6)+b32*f(7)-b44*f(8)+c11*f(10)
             fp(9)=4.d0*f(4)-s11*f(6)+b42*f(7)-b33*f(9)-c22*f(10)   &
                  +b914*f(14)
             fp(10)=4.d0*f(5)+b52*f(7)+s11*f(8)-s22*f(9)-b11*f(10)  &
                  +b914*f(12)
             fp(11)=-t31*f(2)+s13*f(7)+s23*f(9)-t11*f(11)+t21*f(12) &
                  -s11*f(13)+s12*f(14)
             fp(12)=t31*f(3)+s23*f(7)-s13*f(8)+t12*f(11)-t22*f(12)  &
                  +s12*f(13)-s22*f(14)
             fp(13)=s23*f(6)-c11*f(11)+t11*f(13)-t12*f(14)
             fp(14)=-t31*f(1)+s13*f(6)-c22*f(12)-t21*f(13)+t22*f(14)
             fp(3)=s22*f(1)+b32*f(2)+b33*f(3)+c11*f(5)+b313*f(13)
             fp(4)=-s11*f(1)+b42*f(2)+b44*f(4)-c22*f(5)+b414*f(14)
             fp(5)=b52*f(2)+s11*f(3)-s22*f(4)+b55*f(5)-b313*f(11)   &
                    +b414*f(12)
!-------------------------------------------
          else
!intergration de bas en haut
             b11=t22+t33
             b22=t11+t33
             b33=t11+t22
             b55=t22-t33
             b66=t11-t33
             b99=t11-t22
             t4=f(4)+f(8)
             t5=t4+f(8)
             t4=t4+f(4)
             fp(1)=b11*f(1)-t21*f(2)-t31*f(3)-4.d0*f(5)+c22*f(7)
             fp(2)=-t12*f(1)+b22*f(2)-4.d0*f(6)+c11*f(11)
             fp(3)=b33*f(3)-c22*f(9)+c11*f(12)
             fp(4)=-s23*f(1)+s13*f(2)+t31*f(6)
             fp(5)=s13*f(3)+b55*f(5)-t21*f(6)-c22*f(10)
             fp(6)=s23*f(3)-t12*f(5)+b66*f(6)-c11*f(13)
             fp(7)=s22*f(1)-s12*f(2)-b55*f(7)+t31*f(9)+4.d0*f(10)  &
                  +t12*f(11)
             fp(8)=s23*f(1)-s12*f(3)-t21*f(9)+t12*f(12)
             fp(9)=s23*f(2)-s22*f(3)-t12*t5+b99*f(9)-c11*f(14)
             fp(10)=s23*(f(4)-f(8))-s22*f(5)+s12*f(6)+s13*f(9)     &
                  -b11*f(10)+t12*f(13)
             fp(11)=-s12*f(1)+s11*f(2)-t4*t31+t21*f(7)-b66*f(11)   &
                  +4.d0*f(13)
             fp(12)=-s13*f(1)+s11*f(3)+t21*t5-t31*f(5)-b99*f(12)   &
                  +c22*f(14)
             fp(13)=-t4*s13+s12*f(5)-s11*f(6)+t21*f(10)-s23*f(12)  &
                  -b22*f(13)
             fp(14)=s12*t5-s13*f(7)-s11*f(9)+t31*f(10)-s23*f(11)   &
                  +s22*f(12)-b33*f(14)
          endif
       endif
    endif
!--------------------------------------------------------------------------
  end subroutine derms
!--------------------------------------------------------------------------
      
!--------------------------------------------------------------------------
  subroutine dermf(iq,z,f,fp,iknt,qff,nvefm)
! Calcules les derivees des mineurs (fp) dans un liquide
!--------------------------------------------------------------------------
    use bits
    use param_modele
    use yannos_flag
    implicit none
!
    integer, intent(in) :: iq, iknt,nvefm
    doubleprecision, intent(in) :: z,qff
    doubleprecision, dimension(nvefm), intent(in)  :: f
    doubleprecision, dimension(nvefm), intent(out) :: fp
!
    doubleprecision ::  ro,gr,flu,t,zr
    doubleprecision, save :: t21,t11,s11,c11,t12,t22,s22,b11,s12,b33


!
!-----------------------------
    if(iknt/=0) then
! mode count
!-----------------------------
       if(kg==0) then
! pas de gravite
          fp(1)=t11*f(1)+c11*f(2)
          fp(2)=(s11-t21*ro)*f(1)-t11*f(2)
       else
! gravite
          fp(3)=s22*f(1)-(t12+t12)*f(2)+b33*f(3)+c11*f(5)
          fp(4)=-s11*f(1)+(t21+t21)*f(2)-b33*f(4)-4.d0*f(5)
          fp(5)=-(s12+s12)*f(2)+s11*f(3)-s22*f(4)-b11*f(5)
       endif
    else
       t=z-r(iq)
       if(abs(t) <=1.d-15) then
!pas besoin d'interpolation
          ro=rho(iq)
          flu=fcon(iq)*qff
          gr=g(iq)
       else
!interpolation
          ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
          flu=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq)))) &
               *qff
          gr=g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq)))
       endif
       if (cancel_gravity) then
          t21=0.0d0
       else
          t21=-4.d0*ro  !a26
       endif
       zr=1.d0/z
       t12=fl3*zr*zr/wsq  !l(l+1)/(r^2w^2)
       t11=gr*t12-zr   !a25 ?
       s11=ro*(gr*gr*t12-wsq)+t21*gr*zr !a21
       c11=-t12/ro+1.d0/flu  !a12
!-------------------------------------------
       if(kg==0)  then 
!pas de gravite
          fp(1)=t11*f(1)+c11*f(2)
          fp(2)=(s11-t21*ro)*f(1)-t11*f(2)
!-------------------------------------------
       else
!gravite
          t22=-fl*zr
          s22=ro*t12
          b11=t11+t22
          s12=ro*b11
          if(iback/=1) then
!intergration de haut en bas
             b33=t11-t22
             fp(1)=b11*f(1)+4.d0*f(3)-c11*f(4)
             fp(2)=s12*f(1)-t21*f(3)+t12*f(4)
             fp(3)=s22*f(1)-(t12+t12)*f(2)+b33*f(3)+c11*f(5)
             fp(4)=-s11*f(1)+(t21+t21)*f(2)-b33*f(4)-4.d0*f(5)
             fp(5)=-(s12+s12)*f(2)+s11*f(3)-s22*f(4)-b11*f(5)
          else
!intergration de bas en haut
             fp(1)=t22*f(1)-t21*f(2)-4.d0*f(3)
             fp(2)=-t12*f(1)+t11*f(2)-c11*f(4)
             fp(3)=-s22*f(1)+s12*f(2)-t22*f(3)+t12*f(4)
             fp(4)=s12*f(1)-s11*f(2)+t21*f(3)-t11*f(4)
          endif
       endif
    endif
!--------------------------------------------------------------------------
  end subroutine dermf
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
  subroutine sfbm(ass,kg,iback)
!convertit le vecteur de mineur a une discontinuite solide/fluide
!--------------------------------------------------------------------------
    implicit none
!    
    integer, intent(in) :: kg,iback
    doubleprecision, dimension(14), intent(inout) ::  ass
!
    integer :: j
    doubleprecision, dimension(14) :: as
! 
    do  j=1,14
       as(j)=ass(j)
       ass(j)=0.d0
    enddo
!
    if(iback/=1) then
!intergration de haut en bas
       if(kg==0) then
! pas de gravite
          ass(1)=as(3)
          ass(2)=as(5)
       else
! gravite
          ass(1)=as(8)
          ass(2)=-as(12)
          ass(3)=as(3)
          ass(4)=-as(10)
          ass(5)=as(5)
       endif
    else
!intergration de bas en haut
       if(kg==0) then
! pas de gravite
          ass(1)=-as(3)
       else
! gravite
          ass(1)=as(7)
          ass(2)=-as(9)
          ass(3)=-as(10)
          ass(4)=-as(14)
       endif
    endif
!--------------------------------------------------------------------------
  end subroutine sfbm
!--------------------------------------------------------------------------


!--------------------------------------------------------------------------
  subroutine fsbm(ass,kg,iback)
!convertit le vecteur de mineur a une discontinuite fluide/solide
!--------------------------------------------------------------------------
    implicit none
!    
    integer, intent(in) :: kg,iback
    doubleprecision, dimension(14), intent(inout) ::  ass
!
    integer :: j
    doubleprecision, dimension(14) :: as
!
    do  j=1,14
       as(j)=ass(j)
       ass(j)=0.d0
    enddo
!
    if(iback/=1) then
!intergration de haut en bas
       if(kg==0) then
! pas de gravite
          ass(1)=as(1)
          ass(4)=-as(2)
       else
! gravite
          ass(6)=as(1)
          ass(14)=as(2)
          ass(1)=as(3)
          ass(9)=as(4)
          ass(4)=-as(5)
       endif
    else
!intergration de bas en haut
       if(kg==0) then 
! pas de gravite
          ass(1)=-as(1)
       else
! gravite
          ass(1)=-as(1)
          ass(3)=-as(2)
          ass(5)=-as(3)
          ass(12)=as(4)
       endif
    endif
!--------------------------------------------------------------------------
    end subroutine fsbm
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine grav(g,rho,qro,r,n)
!*** given rho and spline coeffs,computes gravity ***
!--------------------------------------------------------------------------
  use global_prameter
  implicit none  
  save
  doubleprecision, dimension(3,ncm), intent(in ) :: qro
  doubleprecision, dimension(:)  , intent(in ) :: rho,r
  integer ,  intent(in ) :: n
  doubleprecision, dimension(:)  , intent(out) :: g
!
  doubleprecision :: del,rn2,trn,c1,c2,c3,c4,c5
  integer :: i,im1
!
  g(1)=0.d0
  do  i=2,n
     im1=i-1
     del=r(i)-r(im1)
     rn2=r(im1)*r(im1)
     trn=2.d0*r(im1)
     c1=rho(im1)*rn2
     c2=(qro(1,im1)*rn2+trn*rho(im1))*0.5d0
     c3=(qro(2,im1)*rn2+trn*qro(1,im1)+rho(im1))/3.d0
     c4=(qro(3,im1)*rn2+trn*qro(2,im1)+qro(1,im1))*.25d0
     c5=(trn*qro(3,im1)+qro(2,im1))*0.2d0
     g(i)=(g(im1)*rn2+4.d0*del*(c1+del*(c2+del*(c3+del*(c4+del* &
          (c5+del*qro(3,im1)/6.d0))))))/(r(i)*r(i))
  enddo
!--------------------------------------------------------------------------
end subroutine grav
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
  subroutine startl(jf,jl,v,ls,q)
!--------------------------------------------------------------------------
!*** finds start level between jf and jl using velocityv and ang. ord. l.
!*** upon entry q is the value of the exponent at r(jf) or at the turning
!*** point(q=0) depending on previous calls to startl. upon exit q is the
!*** value of the exponent at the starting level ls.
    use bits
    use rindex
    use param_modele
    use startl_buffer
!
    implicit none
!
    integer, intent(in) :: jf,jl
    doubleprecision, dimension(:), intent(in) :: v
    doubleprecision, intent(inout) :: q
    integer, intent(out) :: ls
!
    integer :: j,k,i
    doubleprecision :: pp
    logical :: go_on
!
!alloue rrlog et q si on est dans startl pour la premiere fois
    if(first_in_startl) then
       call allocate_startl_buffer()
       first_in_startl=.false.
       vertno=-dlog(eps)
       do  i=3,n
          rrlog(i)=.5d0*dlog(r(i)/r(i-1))
       enddo
    endif
!
    p(:)=0.0d0
    j=jf
    pp=1.0d0
    do  while(j<=jl .and. pp > 0.0d0)
       pp=fl3-wsq*r(j)*r(j)*rho(j)/v(j)
       if(pp > 0.d0) then
          p(j)=dsqrt(pp)
          j=j+1
       endif
    enddo
!
    p(j)=0.d0
!
    k=j
    j=j-1
!
    go_on=.true. 
!
    do while (j>jf .and. go_on)
       q=q+rrlog(k)*(p(j)+p(k))
       if(q < vertno) then
          k=j
          j=j-1
          go_on=.true.
       else
          go_on=.false.            
       endif
    enddo
!
    if(j.le.jf) then
       ls=jf
    else
       ls=j
    endif
!
!--------------------------------------------------------------------------
  end subroutine startl
!--------------------------------------------------------------------------



!----------------------------------------------------------------------------
    subroutine renum(mod)
! procedure de renumerotation des harmoniques pour tenir
! compte des modes de Stonley. Sur n_modes, si il y a
! n_sto modes de stonley, le premier trouve est place
! en numero n_modes et le dernier en n_modes-n_sto+1
!----------------------------------------------------------------------------
      use modes
      implicit none
!
      type(grp_modes), intent(inout) :: mod
!
      integer :: nb_modes,ncm
      type(harmonique) :: mod_tmp
      doubleprecision, dimension(:), pointer ::  r
      doubleprecision :: di,int,int2000,seuil
      integer ::  i,n,j,type,np1,nbstone,k
!
      nb_modes=mod%nbh
      r=>mod%r
      NCM=mod%ncm
      seuil=0.7d0
      if (mod%l>=80) seuil=0.8d0
!
      n=1
      nbstone=0
      do while(n <= nb_modes-nbstone)
         int=0.d0
!int2000 contient l'energie du mode en dessous de 2000km
         int2000=0.d0
         do i=2,NCM
            di=abs((r(i)-r(i-1))*mod%h(n)%n%u(i,1))
            int=int+di
write(200+n,*)r(i),mod%h(n)%n%u(i,1)          
!si la profondeur ets + grande que 2000km on calcul
! l'energie dans int2000
            if ((r(ncm)-r(i)).gt.2.d06) int2000=int2000+di
         enddo
         if (int.ge.1.e-3) then
            if ((int2000/int).gt.seuil) then 
               print*,'Mode de stonley pour n=',n, ' (',n+nbstone &
                     ,'), renumerotation ...'
               nbstone=nbstone+1
               mod_tmp%n=>mod%h(n)%n
               do k=n,nb_modes-nbstone
                  mod%h(k )%n=>mod%h(k+1)%n
               enddo
                mod%h(nb_modes-nbstone+1)%n=>mod_tmp%n
            else
               n=n+1
            endif
         endif
      enddo
      print*,nbstone,' modes de stonley trouves'
      mod%nb_stonley=nbstone
!------------------------------------------------------------------
end subroutine renum
!------------------------------------------------------------------

!------------------------------------------------------------------
subroutine clean_minos_memory
!------------------------------------------------------------------
  use mtab
  use global_prameter
  use detqn_buffer
  use startl_buffer
  use modout_buffer
  implicit none
  
  call deallocate_all
  call deallocate_vf
!  call deallocate_startl_buffer
  call deallocate_modout_buffer  
  nmx_backup=0
  ncm_backup=0
!------------------------------------------------------------------
end subroutine clean_minos_memory
!------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine extension(fichier1,fichier2,ext) 
!---------------------------------------------------------------------
!
  implicit none
!
  character(len=*) ::  fichier1,fichier2,ext
  integer  :: lenfichier1,lenext
!
  lenfichier1=index (fichier1,' ') -1  
  lenext=index (ext,' ') -1  
  if ((lenfichier1+lenext).gt.80) then
     print*,'il y a trop de caracter a ',fichier1
  endif
!
  fichier2(1:lenfichier1)=fichier1(1:lenfichier1)
  fichier2(lenfichier1+1:(lenfichier1+lenext))=ext(1:lenext)
!---------------------------------------------------------------------
end subroutine extension
!---------------------------------------------------------------------
!----------------------------------------------------------------------
subroutine modout(wcom,qmod,gcom,wdiff,flag,ioeig,iinfo)
!----------------------------------------------------------------------
  use bits
  use eifx, a=>ar
  use rindex
  use modout_buffer
  use layer
  use yannos_flag
!
  implicit none
!
  doubleprecision, intent(in) :: wcom,qmod,gcom,wdiff
  integer, intent(in) :: ioeig,iinfo
  logical, intent(in) :: flag
!
  doubleprecision :: ww,qq,gc,fll,fnn,bid=0.0d0
  integer :: nn,ll,i,j,nvec,indicerec,k,ii,iii
!
  call allocate_modout_buffer()
!
  if (modout_format == 'ucb') then
!sortie au format "berkeley"
!     fnn=float(nord)
     nn=nord 
     ll=l 
     ww=wcom 
     qq=qmod 
     gc=gcom 
     nvec=2*n
     if(jcom==3) nvec=6*n
! order: u du v dv phi dphi/dr (watch normalizations)
     do i=1,n
        buf(i)=a(1,i)
        j=i+n
        buf(j)=a(2,i)
        if(jcom==3) then
           j=j+n
           buf(j)=a(3,i)
           j=j+n
           buf(j)=a(4,i)
           j=j+n
           buf(j)=a(5,i)
           j=j+n
           buf(j)=a(6,i)
        endif
     enddo
     write(ioeig)nn,ll,ww,qq,gc,(buf(i),i=1,nvec)
  else if (modout_format == 'olm') then
!sortie au format "old minos"
     if (nb_lay/=1) stop 'Dans le module layer, il ne doit y avoir qu''une seule&
                         & couche pour etre compatible avec le format ''old minos'''
     nn=nord 
     ll=l 
     ww=wcom 
     qq=qmod 
     gc=gcom 
     if(jcom==3) then
! order: u du v dv phi dphi/dr (watch normalizations)
        write(ioeig)nn,ll,ww,qq,gc,wdiff,(real(a(1,i)),i=i_la1(1),n) &
                                        ,(real(a(2,i)),i=i_la1(1),n) &
                                        ,(real(a(3,i)),i=i_la1(1),n) &
                                        ,(real(a(4,i)),i=i_la1(1),n) 
     else if(jcom==1) then
        a(3,:)=0.0d0
        a(4,:)=0.0d0
        write(ioeig)nn,ll,ww,qq,gc,wdiff,(real(a(1,i)),i=i_la1(1),n) &
                                        ,(real(a(2,i)),i=i_la1(1),n) &
                                        ,(real(a(3,i)),i=i_la1(1),n) &
                                        ,(real(a(4,i)),i=i_la1(1),n) 
     else 
        write(ioeig)nn,ll,ww,qq,gc,wdiff,(real(a(1,i)),i=i_la1(1),n) &
                                        ,(real(a(2,i)),i=i_la1(1),n) 
     endif
  else
   !sortie au format yann
     nn=nord
     ll=l
     ww=wcom
     qq=qmod
     gc=gcom
     if(jcom.eq.3) then
        nvec=6*nbcou_lay
        ii=0
        do k=1,nb_lay
           do  i=i_la1(k),i_la2(k)
              ii=ii+1
              j=ii
              buf(j)=a(1,i)
              j=ii+nbcou_lay
              buf(j)=a(2,i)
              j=j+nbcou_lay
              buf(j)=a(3,i)
              j=j+nbcou_lay
              buf(j)=a(4,i)
              j=j+nbcou_lay
              buf(j)=a(5,i)
              j=j+nbcou_lay
              buf(j)=a(6,i)
           enddo
        enddo
     else if(jcom==2) then 
        nvec=2*nbcou_lay
        ii=0
        do  i=1,noc
           a(1,i)=0.d0
           a(2,i)=0.d0
        enddo
        do k=1,nb_lay
           do  i=i_la1(k),i_la2(k)
              ii=ii+1
              j=ii
              buf(j)=a(1,i)
              j=j+nbcou_lay
              buf(j)=a(2,i)
           enddo
        enddo
     else if (jcom==1) then
        nvec=6*nbcou_lay
        ii=0
        buf(:)=0.d0
        do k=1,nb_lay
           do  i=i_la1(k),i_la2(k)
              ii=ii+1
              j=ii
              buf(j)=a(1,i)
              j=j+nbcou_lay
              buf(j)=a(2,i)
           enddo
        enddo
     else
        stop 'modout : jcom /= de 1, 2 ou 3 n''est pas code!'
     endif
! normalisation des modes spheroidaux(yann 13/03/97)
     fll=dble(ll)
     if (jcom.eq.3.or.jcom.eq.1) then
        do i=1,nbcou_lay
           buf(i)=buf(i)*ww/wn
           buf(i+nbcou_lay)=buf(i+nbcou_lay)*ww/wn
           if (jcom.ne.1) then
              buf(i+2*nbcou_lay)=buf(i+2*nbcou_lay)*ww/wn/dsqrt(fll*(fll+1))
              buf(i+3*nbcou_lay)=buf(i+3*nbcou_lay)*ww/wn/dsqrt(fll*(fll+1))
           endif
           buf(i+4*nbcou_lay)=buf(i+4*nbcou_lay)*ww/wn
           buf(i+5*nbcou_lay)=buf(i+5*nbcou_lay)*ww/wn
        enddo
     endif
! normalisation des modes toroidaux 
     if (jcom.eq.2) then
        do i=1,nbcou_lay
           buf(i)=buf(i)*ww/wn/dsqrt(fll*(fll+1))
           buf(i+nbcou_lay)=buf(i+nbcou_lay)*ww/wn/sqrt(fll*(fll+1))
        enddo
     endif
! attention l'ordre de sortie va etre un peu tordu pour
! gagner du temps a la relecture.:
     if (jcom.eq.3.or.jcom.eq.1) then
        do i=1,nbcou_lay
           bufout(4*(i-1)+1)=sngl(buf(  nbcou_lay-i+1))
           bufout(4*(i-1)+2)=sngl(buf(2*nbcou_lay-i+1))
           bufout(4*(i-1)+3)=sngl(buf(3*nbcou_lay-i+1))
           bufout(4*(i-1)+4)=sngl(buf(4*nbcou_lay-i+1))
        enddo
        do i=1,nbcou_lay
           bufout((i-1)+1+4*nbcou_lay)=sngl(buf(5*nbcou_lay-i+1))
           bufout((i-1)+1+5*nbcou_lay)=sngl(buf(6*nbcou_lay-i+1))
        enddo
     else if (jcom.eq.2) then
        do i=1,nbcou_lay
           bufout(2*i-1)=sngl(buf(  nbcou_lay-i+1))
           bufout(2*i  )=sngl(buf(2*nbcou_lay-i+1))
        enddo
     endif
!
     inquire(ioeig,nextrec=indicerec)
     write(ioeig,rec=indicerec) nn,ll,ww,qq,gc,(bufout(i),i=1,nvec)
     write(iinfo,*) nn,ll,ww,qq,gc,indicerec,wdiff,flag
  endif
!----------------------------------------------------------------------
end subroutine modout
!----------------------------------------------------------------------
!--------------------------------------------------------------
subroutine spsm_nograve(w,l,r,vs,vp,rho,ass)
!--------------------------------------------------------------
  integer, intent(in) :: l
  doubleprecision, intent(in)  :: w,r,vs,vp,rho
  doubleprecision, dimension(6), intent(out) :: ass
!
  doubleprecision:: y11,y12,y13,y14,y21,y22,y23,y24
!
  call rayleigh(l,w,r,rho,vp,vs,y11,y12,y13,y14,y21,y22,y23,y24)  
!
! le difference de conention explique ces changements:
!
  ass(2)= (y11*y22-y12*y21)/sqrt(dble(l*(l+1)))
  ass(1)=  y11*y23-y13*y21
  ass(3)=  y11*y24-y14*y21
!
  ass(4)=-(y12*y23-y13*y22)
  ass(5)=  y12*y24-y14*y22
!
!  ass(6)= (y13*y24-y14*y23)*sqrt(dble(l*(l+1)))
  ass(6)=0.d0
!----------------------------------------------------------------------
end subroutine spsm_nograve
!--------------------------------------------------------------

!----------------------------------------------------------------------
  subroutine rayleigh(l,om,r,rho,vp,vs,y11,y12,y13,y14,y21,y22,y23,y24)
!----------------------------------------------------------------------
! vertion de la sub rayleigh pour x=om/v*r proche de zero
!
! ATTENTION, cette routine a deux solutions
! soit utilisation de bsf, soit dev limite
! si elle utilise le dev limite, il faut multiplier le sol par
!  x^n
!
! y11,y12,y13,y14 sont a multiplier par jl(om*r/vp)
! y21,y22,y23,y24 sont a multiplier par jl(om*r/vp)
!
! Calcul les solutions spheroidales (1 et 2 des pages 245 et 246 de TAKEUCHI 
! et SAITO 1972) pour une boule homogene 
! INPUT: 
! l entier
! r rayon en m, real(DP)
! om : pulsation
! vp vitesse des ondes p en m/s, real(DP)
! vs vitesse des ondes s en m/s, real(DP)
! lamb: lambda en, kg m^5s^-4 real(DP)
! mu : mu en  kg m^5s^-4  , real(DP)
! OUTPUT:
! y: real(DP) , les solutions 1 et 2 des pages 245 et 246 de TAKEUCHI et SAITO 1972
!----------------------------------------------------------------------
    use def_gparam
    implicit none
!
    integer , intent(in) ::  l
    real(DP), intent(in) ::  om,r,rho,vp,vs
    real(DP), intent(out)::  y11,y12,y13,y14,y21,y22,y23,y24
!
    integer l2
!
    real(DP) x,jl,jlp1,lp2mu,x2,f2mu,bidon,f2x,r2,lamb,mu
    logical, parameter :: liquide=.false.
!
    mu   =rho*vs**2
    lamb  =rho*vp**2-2.00_DP*mu
!
    x=om*r/vp 
    x2=x*x 
    f2x=2.0_DP*x
    f2mu=2.0_DP*mu
    lp2mu=lamb+f2mu 
    r2=r*r 
  !
!
! attention les sol sont multipliees par r par rapport a takeuchi & saito!
!
    jl=1.0_DP
    call bfs(l,x2,1.E-10_DP,jlp1)
    jlp1=(real(l,DP)-jlp1)/x
    if (liquide) then
       y11=(dble(l)*jl-x*jlp1)/r 
       y12=-lamb*x2*jl/r2
       y13=0.0_DP
       y14=0.0_DP
    else
       if (l.ne.0) then
          y11=(dble(l)*jl-x*jlp1)/r 
          y12=(-lp2mu*x2*jl+f2mu*(dble(l*(l-1))*jl+f2x*jlp1))/r2 
          y13=jl/r 
          y14=f2mu*(dble(l-1)*jl-x*jlp1)/r2 
       else
          y11=-x*jlp1/r
          y12=-lp2mu*(x/r)**2*jl+4*mu/r2*x*jlp1
       endif
    endif
!   
    x=om*r/vs 
    x2=x*x 
    f2x=2*x 
    l2=l*l 
!
    if (l.ne.0.and..not.liquide) then
!     
       jl=1.0_DP
       call bfs(l,x2,1.E-10_DP,jlp1)
       jlp1=(real(l,DP)-jlp1)/x
!
       y21=-dble(l*(l+1))*jl/r
       y22=f2mu*(dble(-l*(l2-1))*jl+dble(l*(l+1))*x*jlp1)/r2
       y23=(-dble((l+1))*jl+x*jlp1)/r
       y24=mu*(x2*jl-dble(2*(l2-1))*jl-f2x*jlp1)/r2

    else
       y21=0.0_DP
       y22=0.0_DP
       y23=0.0_DP
       y24=0.0_DP
    endif
!----------------------------------------------------------------------
  end subroutine rayleigh
!----------------------------------------------------------------------
!--------------------------------------------------------------------------
end module minos
!--------------------------------------------------------------------------
