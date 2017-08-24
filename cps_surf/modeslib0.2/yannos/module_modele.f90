!--------------------------------------------------------------------------
module earth_modele
!--------------------------------------------------------------------------
  implicit none
!
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
     real, dimension(:), pointer :: r,rho,vpv,vsv,qkappa,qshear,vph,vsh,eta
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
! pour la version polynomial
    doubleprecision :: rx,dr,r1,r2,rt,val
    doubleprecision, parameter :: tau=1.0d3
    doubleprecision, dimension(5) :: wrk
    integer :: nreg,n,nn,knt,nlay,ind,j,jj
!
    integer :: i,unit=127,icode
!
!test
    call init_modele(mod)
    open(unit,file=name,iostat=icode,status='old',action='read')
    if (icode>0) then
       print*,'probleme d''ouverture de',name
       stop 'read_modele: prob de d''ouverture'
    endif
!
    read(unit,100,iostat=icode) mod%titre
    read(unit,*  ,iostat=icode) mod%ifanis,mod%tref,mod%ifdeck
!---------------------------------
    if (mod%ifdeck==0) then
       print*,'Modele polynomial'
!modele polynomial
!---------------------------------
!lit d'abord le fichier pour compter le mombre de couche!
       read(unit,*) nreg,mod%nic,mod%noc,rx
       n=0
       jj=5
       if(mod%ifanis.ne.0) jj=8
       do  nn=1,nreg
          read(unit,*) nlay,r1,r2
          do  i=1,nlay
             n=n+1
          enddo
          do  j=1,jj
             read(unit,110) (wrk(i),i=1,5)
          enddo
       enddo
       mod%nbcou=n
!on reroule le tout
       rewind unit
       
!cette fois, lecture pour de bon!
       
       call allocate_modele(mod)
       read(unit,100,iostat=icode) mod%titre
       read(unit,*  ,iostat=icode) mod%ifanis,mod%tref,mod%ifdeck
       read(unit,*) nreg,mod%nic,mod%noc,rx
       rx=rx*tau
       n=0
       knt=0
       jj=5
       if(mod%ifanis.ne.0) jj=8
       do  nn=1,nreg
          read(unit,*) nlay,r1,r2
          r1=r1*tau
          r2=r2*tau
          dr=(r2-r1)/float(nlay-1)
          do  i=1,nlay
             n=n+1
             mod%r(n)=r1+dr*float(i-1)
          enddo
          do  j=1,jj
             read(unit,110) (wrk(i),i=1,5)
             do  i=1,nlay
                ind=knt+i
                rt=mod%r(ind)/rx
                val=wrk(1)+rt*(wrk(2)+rt*(wrk(3)+rt*(wrk(4)+rt*wrk(5))))
                if(j.eq.1) mod%rho(ind)=val*tau
                if(j.eq.2) mod%vpv(ind)=val*tau
                if(j.eq.3) mod%vsv(ind)=val*tau
                if(j.eq.4) mod%qkappa(ind)=val
                if(j.eq.5) mod%qshear(ind)=val
                if(mod%ifanis/=0) then
                   if(j.eq.6) mod%vph(ind)=val*tau
                   if(j.eq.7) mod%vsh(ind)=val*tau
                   if(j.eq.8) mod%eta(ind)=val
                endif
             enddo
          enddo
          do i=1,nlay
             ind=knt+i
             print*,ind,mod%r(ind),mod%rho(ind),mod%vpv(ind),mod%vsv(ind)
          enddo
          knt=knt+nlay
          if(mod%ifanis==0) then
             do  i=1,n
                mod%vph(i)=mod%vpv(i)
                mod%vsh(i)=mod%vsv(i)
                mod%eta(i)=1.d0       
             enddo
          endif
       enddo
!---------------------------------
    else
!modele par couche
!---------------------------------
       read(unit,*  ,iostat=icode) mod%nbcou,mod%nic,mod%noc  &
            ,mod%nbceau
!
       call allocate_modele(mod)
!
       read(unit,105,iostat=icode) (mod%r(i),mod%rho(i),mod%vpv(i),mod%vsv(i) &
                      ,mod%qkappa(i),mod%qshear(i),mod%vph(i)       &
                      ,mod%vsh(i),mod%eta(i),i=1,mod%nbcou)
       if (mod%r(mod%nbcou)/=6371000.) then
          print*,'Warning, the earth radius is not 6371km'
          print*,'This program can run anyway, but this will cause problems'
          print*,'when using normal mode summation codes (I know, this is bad :-(  )' 
          print*,'If you went ti run anyway, edit module_modele.f90, remove the stop'
          print*,'and recompile'
          stop
       endif
       if (icode>0) then
          print*,'probleme de lecture dans',name
          print*,i,mod%nbcou
          stop 'read_modele: prob de de lecture'
       endif
!---------------------------------
    endif
!---------------------------------
    close(unit)
!      
100 format(a80)
105 format(f8.0, 3f9.2, 2f9.1, 2f9.2, f9.5)     
110 format(5f9.5)
!--------------------------------------------------------------------------
  end subroutine read_modele
!--------------------------------------------------------------------------

!!$!--------------------------------------------------------------------------
!!$  subroutine read_modele_MPI(name,mod,rang)
!!$!--------------------------------------------------------------------------
!!$    implicit none
!!$    include 'mpif.h'
!!$!
!!$    integer, intent(in) :: rang
!!$    character(len=100), intent(in) :: name
!!$    type(modele), intent(out) :: mod
!!$!
!!$    integer :: i,unit=127,icode,ier,nbv
!!$    logical :: flag
!!$!
!!$    call MPI_BARRIER(MPI_COMM_WORLD,ier)
!!$    call init_modele(mod)
!!$    if (rang==0) then
!!$       inquire(file=name,exist=flag)
!!$       if (.not.flag) then
!!$          print*,'read_modele_MPI: le modele n''existe pas :',name
!!$          stop
!!$       endif     
!!$       open(unit,file=name,status='old',action='read')
!!$!
!!$       read(unit,100,iostat=icode) mod%titre
!!$       read(unit,*  ,iostat=icode) mod%ifanis,mod%tref,mod%ifdeck
!!$       read(unit,*  ,iostat=icode) mod%nbcou,mod%nic,mod%noc  &
!!$            ,mod%nbceau
!!$!
!!$       call allocate_modele(mod)
!!$!
!!$       read(unit,105,iostat=icode) (mod%r(i),mod%rho(i),mod%vpv(i),mod%vsv(i) &
!!$            ,mod%qkappa(i),mod%qshear(i),mod%vph(i)       &
!!$            ,mod%vsh(i),mod%eta(i),i=1,mod%nbcou)
!!$    endif
!!$    call MPI_BCAST(icode,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
!!$    if (icode>0) then
!!$       if (rang==0)  print*,'probleme de lecture (ou d''ouverture) dans',name
!!$       stop 'read_modele_MPI: prob de de lecture'
!!$    endif
!!$!
!!$!broad cast du modele d0 vers tous les procs:
!!$!
!!$    call MPI_BCAST(mod%ifanis,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
!!$    call MPI_BCAST(mod%ifdeck,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
!!$    call MPI_BCAST(mod%nbcou ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
!!$    call MPI_BCAST(mod%nic   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
!!$    call MPI_BCAST(mod%noc   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
!!$    call MPI_BCAST(mod%nbceau,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
!!$    call MPI_BCAST(mod%tref  ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ier)
!!$!
!!$    if (rang/=0) call allocate_modele(mod)
!!$!
!!$    nbv=mod%nbcou
!!$    call MPI_BCAST(mod%r     ,nbv,MPI_REAL ,0,MPI_COMM_WORLD,ier)
!!$    call MPI_BCAST(mod%rho   ,nbv,MPI_REAL ,0,MPI_COMM_WORLD,ier) 
!!$    call MPI_BCAST(mod%vpv   ,nbv,MPI_REAL ,0,MPI_COMM_WORLD,ier)
!!$    call MPI_BCAST(mod%vsv   ,nbv,MPI_REAL ,0,MPI_COMM_WORLD,ier)
!!$    call MPI_BCAST(mod%qkappa,nbv,MPI_REAL ,0,MPI_COMM_WORLD,ier)
!!$    call MPI_BCAST(mod%qshear,nbv,MPI_REAL ,0,MPI_COMM_WORLD,ier)
!!$    call MPI_BCAST(mod%vph   ,nbv,MPI_REAL ,0,MPI_COMM_WORLD,ier)
!!$    call MPI_BCAST(mod%vsh   ,nbv,MPI_REAL ,0,MPI_COMM_WORLD,ier)
!!$    call MPI_BCAST(mod%eta   ,nbv,MPI_REAL ,0,MPI_COMM_WORLD,ier)
!!$    
!!$!      
!!$100 format(a80)
!!$105 format(f8.0, 3f9.2, 2f9.1, 2f9.2, f9.5)     
!!$!--------------------------------------------------------------------------
!!$  end subroutine read_modele_MPI
!!$!--------------------------------------------------------------------------

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
  subroutine get_beta_value(mod,vp,vs,r,rho,lamb,mu)
!--------------------------------------------------------------------------
    use def_gparam
    implicit none
    type(modele), intent(in) :: mod
    real(DP), intent(out) :: vp,vs,r,rho,lamb,mu
    integer :: n
!
    n  =mod%nbcou
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
    mvsh   =maxval(abs(mod1%vsh   (:)))
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
