!--------------------------------------------------------------------------
module earth_modele
!--------------------------------------------------------------------------
  implicit none
!
!----------------------------------------------------------
  type mod_der
!----------------------------------------------------------
     doubleprecision :: vs,vsp,vspp,vp,vpp,vppp,ro,rop,ropp
!----------------------------------------------------------
  end type mod_der
!----------------------------------------------------------

!----------------------------
  type modele
!
! titre
! ifanis,tref,ifdeck
! nbcou,nic,noc,nbceau
! r(i),rho(i),vpv(i),vsv(i),qkappa(i),qshear(i),vph(i),vsh(i),eta(i)
!
! ifanis/=0 : anisotrope
! ifdeck : je sais pas c'etait la
! nbcou : nmbre de couches du modele
! nic :    numero de la couche inner core
! noc   : numero de la couche outer core
! nbceau : nombre de couches de l'ocean (0 = pas d'ocean)
! champs de dimension (nbcou)
! r : rayon en km
! rho : kg/m3
! vpv,vsv,vph,vsh: km/s
! qkappa,qshear : attenuation
! eta : F/(A-2L)
!----------------------------
     character(len=80) :: titre
     integer :: nbcou     
     integer :: ifanis,ifdeck,nic,noc,nbceau
     real    :: tref
     type(mod_der) :: der
     doubleprecision, dimension(:), pointer :: r,rho,vpv,vsv,qkappa,qshear,vph,vsh,eta
     logical :: allocated
!----------------------------
  end type modele
!----------------------------
  
!--------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------
    
!--------------------------------------------------------------------------
  subroutine  init_modele(mod)
!--------------------------------------------------------------------------
    implicit none
    type(modele), intent(out) :: mod
!
    mod%allocated=.false.
    mod%nbceau=-1
!--------------------------------------------------------------------------
  end subroutine init_modele
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
  subroutine allocate_modele(mod,nbcou)
!--------------------------------------------------------------------------
    implicit none
!
    type(modele), intent(inout) :: mod
    integer, optional,intent(in) :: nbcou
!
    if (present(nbcou)) mod%nbcou=nbcou
    if (mod%nbcou==-1) stop 'allocate_modele: le nombre de couches&
         & n''est pas connu'
    if (mod%allocated) then
       print*,'allocate_modele: le modele',mod%titre,' est deja alloue'
       stop
    endif
    allocate(mod%r(mod%nbcou),mod%rho(mod%nbcou),mod%vpv(mod%nbcou)         &
         ,mod%vsv(mod%nbcou),mod%qkappa(mod%nbcou),mod%qshear(mod%nbcou) &
         ,mod%vph(mod%nbcou),mod%vsh(mod%nbcou),mod%eta(mod%nbcou))
!
    mod%r     (:)=0.d0
    mod%rho   (:)=0.d0
    mod%vpv   (:)=0.d0
    mod%vsv   (:)=0.d0
    mod%qkappa(:)=0.d0
    mod%qshear(:)=0.d0
    mod%vph   (:)=0.d0
    mod%vsh   (:)=0.d0
    mod%eta   (:)=0.d0
!
    mod%allocated=.true.
!
!--------------------------------------------------------------------------
  end subroutine allocate_modele
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
  subroutine deallocate_modele(mod)
!--------------------------------------------------------------------------
    implicit none 
!
    type(modele), intent(inout) :: mod 
!
    if (.not.mod%allocated) then
       print*,'deallocate_modele: le modele',mod%titre,' n''est pas alloue??' 
       stop 
    endif 
    deallocate(mod%r,mod%rho,mod%vpv,mod%vsv,mod%qkappa,mod%qshear  &
              ,mod%vph,mod%vsh,mod%eta)
!
    mod%allocated=.false.
!
!--------------------------------------------------------------------------
  end subroutine deallocate_modele
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
  subroutine read_modele(name,mod)
!--------------------------------------------------------------------------
    implicit none
!
    character(*), intent(in) :: name
    type(modele), intent(out) :: mod
!
    integer :: i,unit=127,icode
!
    call init_modele(mod)
    open(unit,file=name,iostat=icode,status='old',action='read')
    if (icode>0) then
       print*,'probleme d''ouverture de',name
       stop 'read_modele: prob de d''ouverture'
    endif
!
    read(unit,100,iostat=icode) mod%titre
    read(unit,*  ,iostat=icode) mod%ifanis,mod%tref,mod%ifdeck
    read(unit,*  ,iostat=icode) mod%nbcou,mod%nic,mod%noc  &
                             ,mod%nbceau
!    read(unit,*  ,iostat=icode) mod%der%ro  ,mod%der%vp  ,mod%der%vs
!    read(unit,*  ,iostat=icode) mod%der%rop ,mod%der%vpp ,mod%der%vsp
!    read(unit,*  ,iostat=icode) mod%der%ropp,mod%der%vppp,mod%der%vspp
!
    call allocate_modele(mod)
!
    read(unit,105,iostat=icode) (mod%r(i),mod%rho(i),mod%vpv(i),mod%vsv(i) &
                   ,mod%qkappa(i),mod%qshear(i),mod%vph(i)       &
                   ,mod%vsh(i),mod%eta(i),i=1,mod%nbcou)
    if (icode>0) then
       print*,'probleme de lecture dans',name
       stop 'read_modele: prob de de lecture'
    endif
    close(unit)
!      
100 format(a80)
105 format(f8.0, 3f9.2, 2f9.1, 2f9.2, f9.5)     
!--------------------------------------------------------------------------
  end subroutine read_modele
!--------------------------------------------------------------------------


!--------------------------------------------------------------------------
  subroutine write_modele(name,mod)
!--------------------------------------------------------------------------
    character(*), intent(in) :: name
    type(modele), intent(in) :: mod
!
    integer :: i,unit=127,icode
!
    open(unit,file=name,iostat=icode,action='write')
    if (icode>0) then
       print*,'probleme d''ouverture de',name
       stop 'write_modele: prob de d''ouverture'
    endif
!
    write(unit,100,iostat=icode) mod%titre
    write(unit,*  ,iostat=icode) mod%ifanis,mod%tref,mod%ifdeck
    write(unit,*  ,iostat=icode) mod%nbcou,mod%nic,mod%noc  &
                             ,mod%nbceau
    write(unit,*  ,iostat=icode) mod%der%ro  ,mod%der%vp  ,mod%der%vs
    write(unit,*  ,iostat=icode) mod%der%rop ,mod%der%vpp ,mod%der%vsp
    write(unit,*  ,iostat=icode) mod%der%ropp,mod%der%vppp,mod%der%vspp
!
    write(unit,105,iostat=icode) (mod%r(i),mod%rho(i),mod%vpv(i),mod%vsv(i) &
                   ,mod%qkappa(i),mod%qshear(i),mod%vph(i)       &
                   ,mod%vsh(i),mod%eta(i),i=1,mod%nbcou)
    if (icode>0) then
       print*,'probleme de lecture dans',name
       stop 'write_modele: prob d''ecriture'
    endif
    close(unit)
!      
100 format(a80)
105 format(f8.0, 3f9.2, 2f9.1, 2f9.2, f9.5)     
!--------------------------------------------------------------------------
  end subroutine write_modele
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
  subroutine get_beta_value(mod,coq,vp,vs,r,rho,lamb,mu)
!--------------------------------------------------------------------------
    use def_gparam
    implicit none
    type(modele), intent(in) :: mod
    logical, intent(in) :: coq
    real(DP), intent(out) :: vp,vs,r,rho,lamb,mu
    integer :: n
!
    if (coq) then
       n=1
    else
       n  =mod%nbcou
    endif
    vp =mod%vpv(n)
    vs =mod%vsv(n)
    r  =mod%r  (n)
    rho=mod%rho(n)
!
    mu    =rho*vs**2
    lamb  =rho*vp**2-2.00_DP*mu       
!
!--------------------------------------------------------------------------
  end subroutine get_beta_value
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
  logical function modele_different(mod,mod1)
!--------------------------------------------------------------------------
    use def_gparam
    implicit none
    type(modele), intent(in) :: mod,mod1
!
    real :: epsi=1.E-9
    real:: mrho,mr,mvpv,mvph,mvsv,mvsh,mqkappa,mqshear,meta
!
    mrho   =maxval(abs(mod1%rho   (:)))
    mr     =maxval(abs(mod1%r     (:)))
    mvpv   =maxval(abs(mod1%vpv   (:)))
    mvph   =maxval(abs(mod1%vph   (:)))
    mvsv   =maxval(abs(mod1%vsv   (:)))
    mvsv   =max(mvsv,1.) !pour eviter les div par 0
    mvsh   =maxval(abs(mod1%vsh   (:)))
    mvsh   =max(mvsh,1.)
    mqshear=maxval(abs(mod1%qshear(:)))
    mqkappa=maxval(abs(mod1%qkappa(:)))
    meta   =maxval(abs(mod1%eta   (:)))
!
    if (mod1%nbcou /=mod%nbcou .or.mod1%ifanis/=mod%ifanis  .or.   &
        mod1%nic   /=mod%nic   .or.mod1%noc   /=mod%noc     .or.   &
        mod1%nbceau/=mod%nbceau.or.mod1%ifdeck/=mod%ifdeck  .or.   &
        abs(mod1%tref-mod%tref)>epsi) then
!
       modele_different=.true.
!
    else if (maxval(abs(mod1%rho   (:)-mod%rho   (:)))/mrho   >epsi .or.   &
             maxval(abs(mod1%r     (:)-mod%r     (:)))/mr     >epsi .or.   &
             maxval(abs(mod1%vpv   (:)-mod%vpv   (:)))/mvpv   >epsi .or.   &
             maxval(abs(mod1%vph   (:)-mod%vph   (:)))/mvph   >epsi .or.   &
             maxval(abs(mod1%vsv   (:)-mod%vsv   (:)))/mvsv   >epsi .or.   &
             maxval(abs(mod1%vsh   (:)-mod%vsh   (:)))/mvsh   >epsi .or.   &
             maxval(abs(mod1%qkappa(:)-mod%qkappa(:)))/mqkappa>epsi .or.   &
             maxval(abs(mod1%qshear(:)-mod%qshear(:)))/mqshear>epsi .or.   &
             maxval(abs(mod1%eta   (:)-mod%eta   (:)))/meta   >epsi ) then
!
       modele_different=.true.
!
    else
!
       modele_different=.false.
!
    endif     
!--------------------------------------------------------------------------
  end function modele_different
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
end module earth_modele
!--------------------------------------------------------------------------
