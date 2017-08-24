!----------------------------------------------------------------------
module mtab
! tableaux de la table de frequences
!----------------------------------------------------------------------
  implicit none

  logical :: mtab_not_allocated=.true.
  doubleprecision, allocatable, dimension(:) :: we,de,wt
  integer        , allocatable, dimension(:) ::  ke
  doubleprecision, allocatable, dimension(:) :: wtry,bm,um
contains
!----------------------------------
  subroutine allocate_mtab(nmx)
!----------------------------------
    implicit none
    integer, intent(in) :: nmx
!
    allocate(wt(nmx))
    allocate(we(2*nmx),de(2*nmx),ke(2*nmx),wtry(nmx),bm(nmx),um(nmx))
!    allocate(wt(215))
!    allocate(we(430),de(430),ke(430),wtry(250),bm(250),um(250))
    wt(:)=0.0d0
    we(:)=0.0d0
    de(:)=0.0d0
    ke(:)=0.0d0
    wtry(:)=0.0d0
    bm(:)=0.0d0
    um(:)=0.0d0
    mtab_not_allocated=.false.
!----------------------------------
  end subroutine allocate_mtab
!----------------------------------
!
!----------------------------------
  subroutine deallocate_mtab
!----------------------------------
    implicit none
!
    deallocate(wt)
    deallocate(we,de,ke,wtry,bm,um)
!----------------------------------
  end subroutine deallocate_mtab
!----------------------------------   
!
!----------------------------------------------------------------------
end module mtab
!----------------------------------------------------------------------
!!$!----------------------------------------------------------------------
!!$module mtab
!!$! tableaux de la table de frequences
!!$!----------------------------------------------------------------------
!!$  implicit none 
!!$  
!!$  doubleprecision, dimension(430) :: we,de 
!!$  integer        , dimension(430) ::  ke 
!!$  doubleprecision, dimension(250) :: wtry,bm,um 
!!$!----------------------------------------------------------------------
!!$end module mtab 
!!$!----------------------------------------------------------------------
!---------------------------------------------------------------------- 
module shanks
! 
! param pour runge-futta-shanks
!----------------------------------------------------------------------
  implicit none
  doubleprecision, dimension(46) :: b
  doubleprecision, dimension(10) :: c
  doubleprecision, dimension( 8) :: step
  doubleprecision :: dx,stepf
  integer         :: in,maxo
!----------------------------------------------------------------------
end module shanks
!----------------------------------------------------------------------
!----------------------------------------------------------------------
module bits
!
! constantes etc...:
! pi    : pi
! rn    : rayon max  du modele de terre
! vn    : cste de normalisation des vitesses (sqrt(pi*G*rhobar*rn**2))
! wn    : cste de normalisation des frequences
! w     : frequence (pas la pulsation) actuelle
! wsq   : w^2
! wray  : (? voir detqn)
! qinv  : (? voir detqn)
! cg    : ? vitesse de groupe
! wgrav : frequqnce jusqu'a la quelle on titent compte de la gravite
! tref  : periode de ref du modele
! fct   : 2.d0*dlog(tref*wdim)/pi, utilisation de la periode de ref
! eps   : precision de la dicho (param d'entree) 
! epsint: precision de l'integration (param d'entree)
! fl    : l (ordre angulaire)
! fl1   : l+1
! fl2   : l(l+1)
! fl3   : l(l+1)
! sfl3  : (l(l+1))^(1/2)
!
! jcom  : 1 mode radial, 2 toroidal, 3spheroidale
! l     : ordre angulaire
! kg    : 0: on ne tient pas compte de la gravite, 1: on en tient compte
! kount : nombre d'harmoniques present avant la frequence courante (cf detqn)
! knsw  : 1 On compte le nombre de branches, 0 on ne compte pas (dans detqn)
! ifanis: param d'entre pour l'anisotropie (pour les modeles polynomiaux 
!                                           je crois)
!        =0 isotrope, =1,=2 anisotrope?
! iback : sens de l'integration (dans detqn) 1:on va de haut en bas
!                                            0:on va de bas en haut
! nord  : n order, numero de la solution (~branche avec stonley)
!cond_limite: 1 bord libre a la surface
!             2 bord rigide " "   " 
!             3 bord libre a la surface et bord rigide  a la base
!----------------------------------------------------------------------
  implicit none
  doubleprecision, parameter :: pi=3.14159265358979d0
  doubleprecision :: rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl, &
          fl1,fl2,fl3,sfl3,epsint,epsint_base
  integer         :: jcom,l,kg,kount,knsw,ifanis,iback,nord,cond_limite

  contains
!-------------------------------------------------------------------
  subroutine set_nogravity()
!-------------------------------------------------------------------
    
    implicit none
!
    kg=0
!
!-------------------------------------------------------------------
  end subroutine set_nogravity
!-------------------------------------------------------------------

!-------------------------------------------------------------------
  subroutine set_pres(eps_in)
!-------------------------------------------------------------------
    implicit none
    doubleprecision, intent(in) :: eps_in
!
    eps=eps_in
!
!-------------------------------------------------------------------
  end subroutine set_pres
!-------------------------------------------------------------------

!-------------------------------------------------------------------
  subroutine set_l(l_in)
!-------------------------------------------------------------------
    implicit none
    integer, intent(in) :: l_in
!
    l=l_in
!
!-------------------------------------------------------------------
  end subroutine set_l
!-------------------------------------------------------------------

!-------------------------------------------------------------------
  subroutine set_comp(j)
!-------------------------------------------------------------------
    implicit none 
    integer, intent(in) :: j
!
    jcom=j
!
!------------------------------------------------------------------- 
  end subroutine set_comp
!-------------------------------------------------------------------

!-------------------------------------------------------------------
  subroutine set_gravity()
!-------------------------------------------------------------------
    
    implicit none
!
    kg=1
!
!-------------------------------------------------------------------
  end subroutine set_gravity
!-------------------------------------------------------------------
!----------------------------------------------------------------------
end module bits
!----------------------------------------------------------------------
!----------------------------------------------------------------------
module eifx
! ar buffer ou sont stocker les mineurs puis les solutions
!----------------------------------------------------------------------
  implicit none
!
  integer, parameter :: nb_max_mineurs=14
! dimension(14,ncm):
  doubleprecision, dimension(:,:), allocatable :: ar,ar_backup
!dimension(ncm):
  integer,         dimension(:),   allocatable :: inorm,jjj
!----------------------------------------------------------------------
end module eifx
!----------------------------------------------------------------------

!----------------------------------------------------------------------
module arem
!----------------------------------------------------------------------
  implicit none
! dimension(6,3,ncm) :
  doubleprecision, dimension(:,:,:), allocatable :: a
!----------------------------------------------------------------------
end module arem
!----------------------------------------------------------------------

!----------------------------------------------------------------------
module rindex
!----------------------------------------------------------------------
! indexs de rayons:
!
! noc   : numero de la couche outer core
! nic   : numero de la couche inner core
! nocp1 : numero de la couche outer core +1
! nicp1 : numero de la couche inner core +1
! nsl   : numero de la  dernier couche avant l'ocean
! nslp1 : numero de la  dernier couche avant l'ocean+1
! n     : nombre de couches du model (pas d'ocean <=> n=nsl)
!----------------------------------------------------------------------
  implicit none
  integer :: nic,noc,nsl,nicp1,nocp1,nslp1,n
!----------------------------------------------------------------------
end module rindex
!----------------------------------------------------------------------

!----------------------------------------------------------------------
module param_modele
!----------------------------------------------------------------------
  implicit none
! dimension(ncm) :
  doubleprecision, dimension(:), allocatable   :: r,fmu,flam,qshear &
                     ,qkappa,xa2,xlam,rho,fcon,lcon,ccon,acon,ncon,g
  real, dimension(:), allocatable   :: ray
! dimension(3,ncm):
 doubleprecision, dimension(:,:), allocatable :: qro,qg,fspl,lspl  &
                    ,nspl,cspl,aspl
!  
!----------------------------------------------------------------------
end module param_modele
!----------------------------------------------------------------------

!----------------------------------------------------------------------
module global_prameter
!ncm : nombre de couches du modelw
!----------------------------------------------------------------------
  implicit none
  integer :: NCM
  logical, private :: allocated=.false.
contains
!----------------------------------------------------------------------
  subroutine allocate_all
!----------------------------------------------------------------------
    use  param_modele
    use eifx
    use arem
    implicit none
!
    if (.not.allocated) then
       allocate(r(ncm),fmu(ncm),flam(ncm),qshear(ncm),qkappa(ncm),xa2(ncm) &
               ,xlam(ncm),rho(ncm),fcon(ncm),lcon(ncm),ccon(ncm),acon(ncm) &
               ,ncon(ncm),g(ncm),ray(ncm))
       allocate(qro(3,ncm),qg(3,ncm),fspl(3,ncm),lspl(3,ncm),nspl(3,ncm)   &
               ,cspl(3,ncm),aspl(3,ncm))
       allocate(ar(nb_max_mineurs,ncm),ar_backup(nb_max_mineurs,ncm))
       allocate(inorm(ncm),jjj(ncm))       
       allocate(a(6,3,ncm))
       ar(:,:) =0.0d0       
       inorm(:)=0
       jjj  (:)=0
       a(:,:,:)=0.0d0
       allocated=.true.
    else
       stop 'allocate_all: deja fait'
    endif
!----------------------------------------------------------------------
  end subroutine allocate_all
!----------------------------------------------------------------------  
!----------------------------------------------------------------------
  subroutine deallocate_all
!----------------------------------------------------------------------
    use  param_modele
    use eifx
    use arem
    implicit none
!
    if (allocated) then
       deallocate(r,fmu,flam,qshear,qkappa,xa2,xlam,rho,fcon,lcon,ccon,acon &
               ,ncon,g,ray)
       deallocate(qro,qg,fspl,lspl,nspl,cspl,aspl)
       deallocate(ar,ar_backup)
       deallocate(inorm,jjj)       
       deallocate(a)
       allocated=.false.
    else
       stop 'deallocate_all: pas alloue?'
    endif
!----------------------------------------------------------------------
  end subroutine deallocate_all
!----------------------------------------------------------------------  
!----------------------------------------------------------------------
end module global_prameter
!----------------------------------------------------------------------

!----------------------------------------------------------------------
module detqn_buffer
!----------------------------------------------------------------------
  implicit none
  public
  doubleprecision, dimension(:), allocatable :: vf
  logical, private :: first_in_detqn=.true.
  
  contains
!----------------------------------------------------------------------
    subroutine allocate_vf
!----------------------------------------------------------------------
      use global_prameter
      implicit none
!
      if (first_in_detqn) then
         allocate(vf(ncm))
         first_in_detqn=.false.
      endif
!----------------------------------------------------------------------
    end subroutine allocate_vf
!----------------------------------------------------------------------

!----------------------------------------------------------------------
    subroutine deallocate_vf
!----------------------------------------------------------------------
      use global_prameter
      implicit none
!
      deallocate(vf)
      first_in_detqn=.true.
!----------------------------------------------------------------------
    end subroutine deallocate_vf
!----------------------------------------------------------------------
!----------------------------------------------------------------------
end module detqn_buffer
!----------------------------------------------------------------------

!----------------------------------------------------------------------
module startl_buffer
!----------------------------------------------------------------------
  implicit none
  public
  doubleprecision, dimension(:), allocatable :: rrlog,p
  doubleprecision :: vertno
  logical, private :: first_in_startl=.true.
  
  contains
!----------------------------------------------------------------------
    subroutine allocate_startl_buffer
!----------------------------------------------------------------------
      use global_prameter
      implicit none
!
      if (first_in_startl) then
         allocate(rrlog(0:ncm+1),p(0:ncm+1))
         rrlog(:)=0.0d0
         p(:)=0.0d0
!         first_in_startl=.false.
      endif
!----------------------------------------------------------------------
    end subroutine allocate_startl_buffer
!----------------------------------------------------------------------

!----------------------------------------------------------------------
    subroutine deallocate_startl_buffer
!----------------------------------------------------------------------
      use global_prameter
      implicit none
!
      deallocate(rrlog,p)
      first_in_startl=.true.
!----------------------------------------------------------------------
    end subroutine deallocate_startl_buffer
!----------------------------------------------------------------------

!----------------------------------------------------------------------
end module startl_buffer
!----------------------------------------------------------------------

!----------------------------------------------------------------------
module modout_buffer
!----------------------------------------------------------------------
  implicit none
  public
!
  logical, private :: first_in_modout=.true.
!
  doubleprecision, dimension(:), allocatable ::  buf
  real           , dimension(:), allocatable ::  bufout
!  
contains
!----------------------------------------------------------------------
  subroutine allocate_modout_buffer
!----------------------------------------------------------------------
    use global_prameter
    implicit none
!
    if (first_in_modout) then
       allocate(buf(6*ncm),bufout(6*ncm))
       first_in_modout=.false.
    endif
!----------------------------------------------------------------------
  end subroutine allocate_modout_buffer
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  subroutine deallocate_modout_buffer
!----------------------------------------------------------------------
    use global_prameter
    implicit none
!
    if (.not.first_in_modout) then
       deallocate(buf,bufout)
       first_in_modout=.true.
    endif
!----------------------------------------------------------------------
  end subroutine deallocate_modout_buffer
!----------------------------------------------------------------------
  
!----------------------------------------------------------------------
end module modout_buffer
!----------------------------------------------------------------------

!----------------------------------------------------------------------
module modes
!----------------------------------------------------------------------
  implicit none
!
!----------------------
  type mode
!----------------------
     integer :: n   !numero de ce mode
     doubleprecision :: w !frequence
     doubleprecision :: q !attenuation
     doubleprecision :: g !vitesse de groupe
!
!     integer :: NCM ! nombre de couches d'echant verticale des fctp
! deplacement (NCM,2) pour rayleigh,(NCM,1) pour love:
     doubleprecision, dimension(:,:), pointer :: u 
! derivee du deplacement  (NCM,2) pour rayleigh,(NCM,1) pour love          
     doubleprecision, dimension(:,:), pointer :: up 
!----------------------
  end type mode
!----------------------

!----------------------
  type harmonique
!----------------------
     type(mode), pointer :: n
!----------------------
  end type harmonique
!----------------------

!----------------------
  type grp_modes
!----------------------
     integer :: type ! 2 :love, 3 Rayleigh
     integer :: l   !ordre angulaire
!
     integer :: NCM ! nombre de couches d'echant verticale des fctp
! rayon non normalise (r(ncm)=6371000.,  si PREM, r(1)=0.)
! c a d que pour les tableaux contenant les modes u(NCM)
!est la surfvace et u(1) le centre de la terre
     doubleprecision, dimension(:), pointer :: r
!
     integer :: nbh ! nombre d'harmoniques presents dans le groupe
     integer :: nb_stonley ! nombre de modes de stonley presents dans le groupe
!tableau de modes:     
     type(mode), dimension(:), pointer :: tab_modes
!tableau de pointers sur les modes, pour manipe
     type(harmonique), dimension(:),pointer :: h
!----------------------
  end type grp_modes
!----------------------

!
contains
!----------------------------------------------------------------------
  subroutine allocate_mode(mod,nb_couches,type)
!----------------------------------------------------------------------
    implicit none
    type(mode), intent(inout) :: mod
    integer,  intent(in) :: nb_couches,type
    integer :: ier
!
    select case(type)
    case(1,2)
!love:
       allocate(mod%u(nb_couches,1),mod%up(nb_couches,1))
       mod%u (:,:)=0.0d0
       mod%up(:,:)=0.0d0       
    case(3)
!rayleigh
       allocate(mod%u (nb_couches,2),stat=ier)
       allocate(mod%up(nb_couches,2),stat=ier)
       mod%u (:,:)=0.0d0
       mod%up(:,:)=0.0d0
    case default
       stop 'allocate_mode: type de mode inconnu'
    end select
!----------------------------------------------------------------------
  end subroutine allocate_mode
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine deallocate_mode(mod)
!----------------------------------------------------------------------
  implicit none 
  type(mode), intent(inout) :: mod 
!
  deallocate(mod%u,mod%up) 
!---------------------------------------------------------------------- 
end subroutine deallocate_mode 
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  subroutine allocate_grp_modes(grp,nb_modes,nb_couches,type)
!----------------------------------------------------------------------
    implicit none
    type(grp_modes) , intent(inout) :: grp
    integer,  intent(in) :: nb_couches,type,nb_modes    
!tableau automatique
    integer :: i,err
!
    grp%NCM =nb_couches
    grp%type=type
    grp%nbh =nb_modes
    allocate(grp%r(nb_couches),stat=err)
    allocate(grp%tab_modes(nb_modes),stat=err)
    allocate(grp%h(nb_modes),stat=err)
    if (err/=0) stop 'allocate_grp_modes: pb d''allocation dynamique'
    do i=1,nb_modes
       call allocate_mode(grp%tab_modes(i),nb_couches,type)
       grp%h(i)%n=>grp%tab_modes(i)
    enddo
!----------------------------------------------------------------------
  end subroutine allocate_grp_modes
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  subroutine deallocate_grp_modes(grp)
!----------------------------------------------------------------------
    implicit none
    type(grp_modes) , intent(inout) :: grp
!
    integer :: i
!
    deallocate(grp%r)
    do i=1,grp%nbh
       call deallocate_mode(grp%tab_modes(i))
    enddo
    deallocate(grp%tab_modes)
    deallocate(grp%h)
!----------------------------------------------------------------------
  end subroutine deallocate_grp_modes
!----------------------------------------------------------------------

!----------------------------------------------------------------------
end module modes
!----------------------------------------------------------------------
