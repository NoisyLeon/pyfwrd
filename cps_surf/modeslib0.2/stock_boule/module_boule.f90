!-------------------------------------------------------
module sol_ana
!--------------------------------------------------------
!
  use def_gparam
  implicit none
!
  public :: get_0, AA1lim, AA2lim, set_param, set_pres, init_tree_freq &
           ,read_freq ,allocate_freq,  deallocate_freq, write_freq     &
           ,open_freq, close_freq, det1, det2, love, rayleigh     &
           ,get_freq_param, filtre_coef, set_filtre,close_boule,get_0_sl &
           ,write_fctp,set_boule_param
!
  integer , dimension(:,:), allocatable, public:: n_zero
  real(DP), dimension(:,:), allocatable, public:: Al00f,Al11f,Al01f &
          ,Al10f,Al5f
  real(DP), dimension(:,:,:), allocatable, public :: f
!
  private
     integer, parameter  :: uARf =1009      &
                           ,uALf =5559      &
                           ,uinfo=5569

      real(DP) ::  r,vp,vs,rho,lamb,mu
      real(DP) :: df=0.1E-5_DP,fmax,fmin
      integer  :: ipf
      real(DP) :: eps=1.0E-10_DP,nbdiv=100.0_DP
      integer  :: Npre=100,Ipmax,Lmax,lmax_pres_l=0,lmax_pres_r=0
      real(DP) :: f_filtre_1,f_filtre_2,f_filtre_3,f_filtre_4

     character(len=100) ::  tree_path= &
          '/donnee/user3/ipg/capdevil/croute/flt/int_freq_ana ' &
          ,tree
     character(len=100) ::  finfo,fARf,fALf

      logical  :: param_set=.false., init=.false., ouvert=.false.       &
                 ,allocated_freq=.false., open_mode, taille_tab=.false. &
                 ,filtre_set=.false.,freq_read=.false.,oldminos=.false.
! info pour faire des modes propres d'une boule homogene:
!
      real(DP) :: dr
      real(DP), dimension(:), allocatable :: ray2
      integer :: NCM ! nombre de points pour une fctp
      integer :: NMAX ! nombre d'harmonique calcule
      logical :: surface_libre=.false. !si vrai on cherche les solutions
                                       !d'une boule libre
      integer :: unit_rd, unit_ld,unit_ri, unit_li! unite loghique des fichier rayleigh et love
      integer :: love_rec,rayleigh_rec,len_l,len_r
      character(len=105) :: file_directr,file_infor,file_directl,file_infol
      real(DP) :: scale ! parametres d'echelles pour coller a minos
      real(DP):: rhobar      
!
! pour l'integration de gauss
      real(DP), dimension(:), allocatable :: r_gauss,r_gauss2,w_gauss
      logical ::  gauss_param=.false.
!
!--------------------------------------------------------
contains
!--------------------------------------------------------

!--------------------------------------------------------
  subroutine allocate_freq(Lmax_in,Ipmax_in)
!--------------------------------------------------------
    implicit none
!
    integer, optional, intent(in):: Lmax_in,Ipmax_in
!
    integer :: ier
!
    if (present(Lmax_in).neqv.present(Ipmax_in)) &
       stop 'allocate_freq: si un des 2 parametres est present,&
           & les 2 doivent l''etre!'
    if (present(Lmax_in)) then
       Lmax=Lmax_in
       Ipmax=Ipmax_in
       taille_tab=.true.
       if (mod(lmax,2)==0) then
          print*,'allocate_freq: attention, lmax doit etre impair'
          lmax=lmax+1 !lmax doit etre impair
       endif
    endif
    if (.not.taille_tab) stop 'allocate_freq: Lmax et Ipmax ne sont pas&
                             & present, et ne sont pas initialises!'
    if (allocated_freq) then
       print*,'allocated_freq: attention les tableaux deja&
              & alloues, je ne fais rien!!'
    else
       allocate(n_zero(0:LMAX,2)   ,stat=ier)
       allocate(Al00f(IPMAX,0:Lmax),stat=ier)
       allocate(Al11f(IPMAX,0:Lmax),stat=ier)
       allocate(Al01f(IPMAX,0:Lmax),stat=ier)
       allocate(Al10f(IPMAX,0:Lmax),stat=ier)
       allocate(Al5f (IPMAX,0:Lmax),stat=ier)
       allocate(f(IPMAX,0:Lmax,2)  ,stat=ier)
       if (ier/=0) stop 'allocate_freq: pb d''alloc memoire'
       allocated_freq=.true.
       n_zero(:,:  )=0
       Al00f (:,:  )=0._DP
       Al11f (:,:  )=0._DP
       Al01f (:,:  )=0._DP
       Al10f (:,:  )=0._DP
       f     (:,:,:)=0._DP       
    endif
!
!--------------------------------------------------------
  end subroutine allocate_freq
!--------------------------------------------------------


!--------------------------------------------------------
  subroutine deallocate_freq
!--------------------------------------------------------
    implicit none
!
    if (.not.allocated_freq) then
       print*,'deallocate_freq: attention, les tableaux ne sont&
              & pas alloue'
    else
       deallocate(n_zero,Al00f,Al11f,Al01f,Al10f,Al5f,f)
       allocated_freq=.false.
       freq_read     =.false.
    endif
!
!--------------------------------------------------------
  end subroutine deallocate_freq
!--------------------------------------------------------
  
!--------------------------------------------------------
  subroutine write_freq(lmax_present,Ipmax_present)
!--------------------------------------------------------
  implicit none
  integer, intent(in) :: lmax_present,Ipmax_present

  integer :: l,ip
!
  if (.not.allocated_freq) stop 'write_freq: vous devez avoir alloue les&
                           & tableaux avant!'
  
  if (.not.ouvert) stop 'write_freq: vous devez avoir ouvert les&
                           & fichiers avant!'
  if (.not.open_mode) stop 'write_freq: vous devez avoir ouvert les&     
                           & fichiers en ecriture!'
  if (.not.param_set) print*,'write_freq: attention les parametres n''ont&
                             & pas etes fixes par set_param'
  write(uinfo,*)lmax_present,Ipmax_present
  write(uinfo,*) vp,vs,r,rho,lamb,mu,fmin,fmax                       &
                ,(n_zero(l,1),l=0,lmax_present)                      &
                ,(n_zero(l,2),l=0,lmax_present)
  do l=0,lmax_present
     write(uARf,*) l,n_zero(l,1)
     do ip=1,n_zero(l,1)
        write(uARf,*) f(ip,l,1),Al00f(ip,l),Al01f(ip,l),Al10f(ip,l)  &
                     ,Al11f(ip,l)
     enddo
     write(uALf,*) l,n_zero(l,2)
     do ip=1,n_zero(l,2)
        write(uALf,*) f(ip,l,2),Al5f(ip,l)
     enddo
  enddo
!--------------------------------------------------------
  end subroutine write_freq
!--------------------------------------------------------

!--------------------------------------------------------
  subroutine read_freq(fermeture_in)
!--------------------------------------------------------
    implicit none
!
    logical, optional,intent(in) :: fermeture_in
    logical  :: fermeture
    integer :: l,ip,l_verif,n_verif
!
    if (present(fermeture_in)) then
       fermeture=fermeture_in
    else
       fermeture=.true.  !valeur pas defaut
    endif
    
    if (.not.init)   call init_tree_freq()
    if (.not.ouvert) call open_freq()
    
    if (open_mode)  stop 'read_freq: vous devez avoir ouvert les&     
                           & fichiers en lecture!'    
    read(uinfo,*)lmax,Ipmax
    taille_tab=.true.
    if (.not.allocated_freq) then
       call allocate_freq
    else
       stop 'read_freq: les tableaux sont deja alloues avant de&
           & connaitre leurs tails'
    endif
    read(uinfo,*) vp,vs,r,rho,lamb,mu,fmin,fmax       &
         ,(n_zero(l,1),l=0,lmax)                      &
         ,(n_zero(l,2),l=0,lmax)
    param_set=.true.
    do l=0,lmax
       read(uARf,*) l_verif,n_verif
       if (n_verif /= n_zero(l,1)) stop 'read_freq: n_verif /= n_zero(l,1)??'
       do ip=1,n_zero(l,1)
          read(uARf,*) f(ip,l,1),Al00f(ip,l),Al01f(ip,l),Al10f(ip,l)  &
               ,Al11f(ip,l)
       enddo
       read(uALf,*) l_verif,n_verif
       if (n_verif /= n_zero(l,2)) stop 'read_freq: n_verif /= n_zero(l,2)??'
       do ip=1,n_zero(l,2)
          read(uALf,*) f(ip,l,2),Al5f(ip,l)
       enddo
    enddo
    freq_read=.true.
!    if (mod(lmax,2)==0) then
!       lmax=lmax+1 !lmax doit etre impair
!       print*,'read_freq: attention, lmax doit etre impair'
!    endif
    if (fermeture) call close_freq
!--------------------------------------------------------
  end subroutine read_freq
!--------------------------------------------------------

!--------------------------------------------------------
  subroutine close_freq
!--------------------------------------------------------
    implicit none
    close(uinfo)
    close(uARf)
    close(uALf)
    ouvert=.false.    
!--------------------------------------------------------
  end subroutine close_freq
!--------------------------------------------------------

!------------------------------------------------------------------------
    subroutine init_tree_freq(creat)
!------------------------------------------------------------------------
! initialise l'arborescence pour le stockage (creat=.true.)
! ou la lecture (creat=.false.) des fonctions de Legendre
!------------------------------------------------------------------------
!
      logical, optional, intent(in) :: creat
!------------------------------------------------------------------------
      logical           :: deja,creat_in
      integer           :: len
!
      if (.not.present(creat)) then
         creat_in=.false.   !valeur par defaut 
      else
         creat_in=creat
      endif
      open_mode=creat_in
      if (.not.init) then
         tree_path=adjustl(tree_path)
         len=LEN_TRIM(tree_path)+1
         tree(1:len)=tree_path(1:len)
         tree(len:)=' '
!
         inquire (file=tree,exist=deja)
         if (.not.deja) then
            if (creat_in) then
               call system('mkdir '//tree)
               write(*,*) 'creation de: ',tree
            else
               write(*,*) tree,' n''existe pas!'
               write(*,*) 'il n''a donc pas ete remplis par les fcts de legendre'
               write(*,*) 'il inutile de continuer: je stope.'
               stop 'init_tree_freq: tree n''existe pas'
            endif
         endif
         tree(len:len)='/'
!
         finfo=tree(1:len)//'info'
         fARf=tree(1:len)//'aRf'
         fALf=tree(1:len)//'aLf'
         print*,'finfo:',finfo
         print*,'fARf :',fARf
         print*,'fALf :',fALf
!
         init=.true.
      else
         print*,'init_tree_freq: l''arborescence a deja ete initialisee!'
      endif
!---------------------------------------------------------------------      
    end subroutine init_tree_freq
!---------------------------------------------------------------------
!----------------------------------------------------------------------------
    subroutine open_freq(creat)
!----------------------------------------------------------------------------
! ouvre les fichiers pour stocjage ou lecture de fqp + fctp.
! uA... sont les numeros d'unite logique des fichier
! creat est un logical comme dans init_tree
!----------------------------------------------------------------------------
      implicit none
!
      logical, optional, intent(in) :: creat
      logical                       :: creat_in
!----------------------------------------------------------------------------
      if (.not.present(creat)) then
         creat_in=.false.  !valeur par defaut
      else
         creat_in=creat
      endif
      if (.not.(creat_in.and.open_mode).and.(creat_in.or.open_mode)) then
         print*,'open_freq: ouverture dans un mode different que celui&
                 & de init_flt!'
         stop 'open_freq: conflit de mode de lecture/ecriture'
      endif
      if (ouvert)   stop  'open_freq: les fichiers sont DEJA ouverts'
      if (.not.init) stop 'open_freq: l''arborescence n''a pas ete initialise&
              & avec ini_tree_freq, impossible d''ouvrir les fichiers'
      if (creat_in) then
         open(uinfo,file=finfo,action='write',status='replace' &
              ,position='rewind',form='formatted')
         open(uARf,file=fARf,action='write',status='replace' &
              ,position='rewind',form='formatted')
         open(uALf,file=fALf,action='write',status='replace' &
              ,position='rewind',form='formatted')
      else
         open(uinfo,file=finfo,action='read',status='old' &
              ,position='rewind',form='formatted')
         open(uARf,file=fARf,action='read',status='old' &
              ,position='rewind',form='formatted')
         open(uALf,file=fALf,action='read',status='old' &
              ,position='rewind',form='formatted')
      endif
      ouvert=.true.
!---------------------------------------------------------------------
    end subroutine open_freq
!---------------------------------------------------------------------

!---------------------------------------------------------------------
  subroutine set_boule_param(r_in,vp_in,vs_in,rho_in,lamb_in,mu_in,fmin_in &
                      ,fmax_in,NCM_in,NMAX_in,LMAX_in,fileR,fileL,df_in)
!---------------------------------------------------------------------
    use def_gparam
    implicit none
!
    real(DP), intent(in) :: r_in,vp_in,vs_in,lamb_in,mu_in,fmin_in    & 
                           ,fmax_in,rho_in
    integer,  intent(in):: NCM_in,NMAX_in,LMAX_in
    real(DP), optional, intent(in) :: df_in   
    real(SP) :: test_sp
    character(len=100), intent(in) :: fileR,fileL
    integer :: len,len_dp,len_int,len_sp,i
!
    r=r_in
    if (r_in <= 0.0_DP) stop 'set_param: r doit etre >= 0!'
    vp=vp_in
    if (vp_in <= 0.0_DP) stop 'set_param: vp doit etre >= 0!'
    vs=vs_in
    if (vs_in <= 0.0_DP) stop 'set_param: vs doit etre >= 0!'
    rho=rho_in
    if (rho_in <= 0.0_DP) stop 'set_param: rho doit etre >= 0!'
    lamb=lamb_in
!    if (lamb_in <= 0.0_DP) stop 'set_param: lamb doit etre >= 0!'
    mu=mu_in
!    if (mu_in <= 0.0_DP) stop 'set_param: mu doit etre >= 0!'
    fmin=fmin_in     
    if (fmin_in <= 0.0_DP) stop 'set_param: fmin doit etre >= 0!'
    fmax=fmax_in
    if (fmax_in <= 0.0_DP) stop 'set_param: fmax doit etre >= 0!'
    if (fmax <= fmin)  stop 'set_param: fmin doit etre > fmax!'
    if (present(df_in)) df=df_in
    if (df <= 0.0_DP) stop 'set_param: df doit etre > 0!'
    NCM=NCM_in
    NMAX=NMAX_in
    Lmax=Lmax_in
    ipf=fmax/df+1
    dr=r/NCM 
    param_set=.true.
    surface_libre=.true.
    unit_rd=133
    unit_ri=134
    unit_ld=135
    unit_li=136
    call extension(fileR,file_directr,'.direct ')
    call extension(fileR,file_infor,'.info ')
    call extension(fileL,file_directl,'.direct ')
    call extension(fileL,file_infol,'.info ')
    if (oldminos) then
       print*,'**************************************************'
       print*,'********* Format de sorties oldminos **************'
       print*,'**************************************************'
       open(unit_rd,file=fileR,form='unformatted')
       open(unit_ld,file=fileL,form='unformatted')
    else
       inquire(iolength=len_dp)r
       inquire(iolength=len_sp)test_sp
       inquire(iolength=len_int)ipf
       len_r=len_int*2+3*len_dp+4*NCM*len_sp
       open(unit_rd,file=file_directr,access='DIRECT',recl=len_r)    
       len_l=len_int*2+3*len_dp+2*NCM*len_sp
       open(unit_ld,file=file_directl,access='DIRECT',recl=len_l)    
       open(unit_ri,file=file_infor)
       open(unit_li,file=file_infol)
    endif
    love_rec=1
    rayleigh_rec=1
    allocate(n_zero(0:LMAX,2))
    allocate(f(NMAX,0:Lmax,2))
    allocate(ray2(NCM))
    do i=1,NCM
!       ray2(i)=(r-(i-1)*dr)**2
       ray2(i)=(i*dr)**2
    enddo
    rhobar=5515.0_DP
    scale=1.0_DP/sqrt(rhobar*r**3)
!---------------------------------------------------------------------
  end subroutine set_boule_param
!---------------------------------------------------------------------

!---------------------------------------------------------------------
  subroutine set_param(r_in,vp_in,vs_in,rho_in,lamb_in,mu_in,fmin_in &
                      ,fmax_in,df_in)
!---------------------------------------------------------------------
    use def_gparam
    implicit none
!
    real(DP), intent(in) :: r_in,vp_in,vs_in,lamb_in,mu_in,fmin_in    & 
                           ,fmax_in,rho_in
    real(DP), optional, intent(in) :: df_in
!
    r=r_in
    if (r_in <= 0.0) stop 'set_param: r doit etre >= 0!'
    vp=vp_in
    if (vp_in <= 0.0) stop 'set_param: vp doit etre >= 0!'
    vs=vs_in
    if (vs_in <= 0.0) stop 'set_param: vs doit etre >= 0!'
    rho=rho_in
    if (rho_in <= 0.0) stop 'set_param: rho doit etre >= 0!'
    lamb=lamb_in
    if (lamb_in <= 0.0) stop 'set_param: lamb doit etre >= 0!'
    mu=mu_in
    if (mu_in <= 0.0) stop 'set_param: mu doit etre >= 0!'
    fmin=fmin_in     
    if (fmin_in <= 0.0) stop 'set_param: fmin doit etre >= 0!'
    fmax=fmax_in
    if (fmax_in <= 0.0) stop 'set_param: fmax doit etre >= 0!'
    if (fmax <= fmin)  stop 'set_param: fmin doit etre > fmax!'
    if (present(df_in)) df=df_in
    if (df <= 0.0) stop 'set_param: df doit etre > 0!'
    ipf=fmax/df+1
    param_set=.true.
!---------------------------------------------------------------------
  end subroutine set_param
!---------------------------------------------------------------------

!---------------------------------------------------------------------
  subroutine get_freq_param(lmax_out,Ipmax_out,r_out,vp_out, &
         vs_out,rho_out,lamb_out,mu_out,fmin_out,fmax_out,df_out)
!---------------------------------------------------------------------
    use def_gparam
    implicit none
!
    real(DP), optional,intent(out) :: r_out,vp_out,vs_out   &
          ,lamb_out,mu_out,fmin_out,fmax_out,rho_out,df_out
    integer, optional, intent(out) :: lmax_out,Ipmax_out
!
    if (present(lmax_out )) lmax_out =lmax
    if (present(Ipmax_out)) Ipmax_out=Ipmax
    if (present(r_out    )) r_out    =r
    if (present(vp_out   )) vp_out   =vp
    if (present(vs_out   )) vs_out   =vs
    if (present(rho_out  )) rho_out  =rho
    if (present(lamb_out )) lamb_out =lamb
    if (present(mu_out   )) mu_out   =mu
    if (present(fmin_out )) fmin_out =fmin
    if (present(fmax_out )) fmax_out =fmax
    if (present(df_out   )) df_out   =df
!---------------------------------------------------------------------
  end subroutine get_freq_param
!---------------------------------------------------------------------

!---------------------------------------------------------------------
  subroutine set_pres(eps_in,nbdiv_in,Npre_in)
!---------------------------------------------------------------------
    use def_gparam
    implicit none
!
    real(DP), optional, intent(in) :: eps_in, nbdiv_in
    integer , optional, intent(in) :: Npre_in
!
    if (present(eps_in))   eps  =eps_in
    if (present(nbdiv_in)) nbdiv=nbdiv_in
    if (present(eps_in))   Npre =Npre_in
!
!---------------------------------------------------------------------
  end subroutine set_pres
!---------------------------------------------------------------------

!---------------------------------------------------------------------
  subroutine get_0(l,det,fdeb,f_zero,n_zero)
!---------------------------------------------------------------------
    use def_gparam
    implicit none
!
    real(DP), external   ::  det
    real(DP), intent(in) :: fdeb
    integer , intent(in) :: l
    integer , intent(out):: n_zero
    real(DP), dimension(:), intent(out):: f_zero
!
    real(DP)  :: prod,val2,val1,f1,f2,ddf,f,detf,ap,am,detfm1
    integer   :: k,ipdeb,ip
!
    if (.not.param_set) stop 'get_0: il faut passer par set_param avant!'
    n_zero=0
    f_zero(:)=0.0_DP
!
! 1 on cherche les singularites:, c.a.d les zeros de det!
!
    detfm1=det(l,df)
    ipdeb=max(int(fdeb/df)-1,1)
    do ip=ipdeb,ipf
       f=(ip-1)*df
       detf=det(l,f)
!       print*,'detf=',detf
       prod=detfm1*detf
       if (prod.lt.0.0d0.and.(ip*df).lt.fmax.and.f.ge.fmin) &
            then
!     
!     il y a un zero entre ip-1 et ip, on cherche a 2^Npre pres
!     
          f1=df*(ip-2)        !f(ip)=(ip-1)*df
          f2=df*(ip-1)
          ddf=df/2.0d0
          f =f1+ddf
          val1=detfm1
          val2=detf
          detfm1=detf
          do k=1,Npre
             detf=det(l,f)
             if (val1*detf.le.0.0d0) then
                val2=detf
                f2=f
             else
                f1=f
                val1=detf
             endif
             ddf=ddf/2.0d0
             f=f1+ddf
          enddo
          if (f.lt.fmax.and.f.ge.fmin) then
             n_zero=n_zero+1
             f_zero(n_zero)=f
          endif
       else
          detfm1=detf
       endif
    enddo
    print*,'nombre de singularite(s):',n_zero

!--------------------------------------------------------
  end subroutine get_0
!--------------------------------------------------------


!---------------------------------------------------------------------
  subroutine get_0_sl(l,det,fdeb,f_zero,n_z,type)
! cherches le fonctions propres d'une boule homogene a surface libre
! type=0: rayleight
! type=1: love
!---------------------------------------------------------------------
    use def_gparam
    implicit none
!
    real(DP), external   :: det
    real(DP), intent(in) :: fdeb
    integer , intent(in) :: l,type
    integer , intent(out):: n_z
    real(DP), dimension(:), intent(out):: f_zero
!
    real(DP)  :: prod,val2,val1,f1,f2,ddf,f,detf,ap,am,detfm1
    integer   :: k,ipdeb,ip
!
    if (.not.param_set) stop 'get_0: il faut passer par set_param avant!'
    surface_libre=.true.
    n_z=0
    f_zero(:)=0.0_DP
!
! 1 on cherche les singularites:, c.a.d les zeros de det!
!
    ipdeb=max(int(fdeb/df)-1,1)
    detfm1=det(l,max((ipdeb-1)*df,df))
    do ip=ipdeb,ipf
       f=(ip-1)*df
       detf=det(l,f)
!       write(107,*) f,detf
       prod=detfm1*detf
       if (prod.lt.0.0d0.and.(ip*df).lt.fmax.and.f.ge.fmin.and.n_z<NMAX) &
            then
!     
!     il y a un zero entre ip-1 et ip, on cherche a 2^Npre pres
!     
          f1=df*(ip-2)        !f(ip)=(ip-1)*df
          f2=df*(ip-1)
          ddf=df/2.0d0
          f =f1+ddf
          val1=detfm1
          val2=detf
          detfm1=detf
          do k=1,Npre
             detf=det(l,f)
             if (val1*detf.le.0.0d0) then
                val2=detf
                f2=f
             else
                f1=f
                val1=detf
             endif
             ddf=ddf/2.0d0
             f=f1+ddf
          enddo
          if (f.lt.fmax.and.f.ge.fmin) then
             n_z=n_z+1
             f_zero(n_z)=f
          endif
       else
          detfm1=detf
       endif
    enddo
    print*,'nombre de singularite(s):',n_z
    if (n_z/=0) then
       if (type==1) lmax_pres_r=max(lmax_pres_r,l)
       if (type==2) lmax_pres_l=max(lmax_pres_l,l)
    endif

!--------------------------------------------------------
  end subroutine get_0_sl
!--------------------------------------------------------

!--------------------------------------------------------
  function det1(l,f) result(det1_out)
!--------------------------------------------------------
!     
    use def_gparam
    implicit none
!  
    integer   :: l
    real(DP)  :: f
!
    real(DP)  ::det1_out
!
    real(DP)  :: om,y11,y12,y13,y14,y21,y22,y23,y24,norme
!
    if (.not.param_set) stop 'det1: il faut passer par set_param avant!'
    om=2.d0*pi*f
    if (f.ne.0.0d0) then
       call rayleigh(l,om,y11,y12,y13,y14,y21,y22,y23,y24)  
       if (l.eq.0) then
          if (surface_libre) then
             norme=abs(y12+y11)        
!             norme=1.0_DP
             if (norme>1.0E-80_DP) then    
                det1_out=-y12/norme
             else
                det1_out=0.0_DP
             endif
          else
             det1_out=-y11
          endif
       else
          if (surface_libre) then
             norme=abs(y14*y22)+abs(y12*y24)
!             norme=1.0_DP
             if (norme>1.0E-80_DP) then
                det1_out=(y14*y22-y12*y24)/norme
             else
                det1_out=0.0_DP
             endif
          else
             det1_out=y13*y21-y11*y23
          endif
       endif
    else
       det1_out=0.0d0
    endif
!    
!-------------------------------------------------------- 
  end function det1
!--------------------------------------------------------
!--------------------------------------------------------
  function det2(l,f) result(det2_out)
!--------------------------------------------------------
!     
    use def_gparam
    implicit none
!
    integer  :: l
    real(DP) :: f
    real(DP) :: det2_out
!
    real(DP) :: om,y5,y6
!
    if (.not.param_set) stop 'det2: il faut passer par set_param avant!'
    om=2.d0*pi*f
    if (f.ne.0.0_DP.and.l.ne.0) then
       call love(l,om,y5,y6)
       if (surface_libre) then
          det2_out=y6
       else
          det2_out=y5
       endif
    else
       det2_out=0.0_DP
    endif
!
!--------------------------------------------------------     
  end function det2
!--------------------------------------------------------      

!--------------------------------------------------------      
  subroutine  AA1lim(l,f_in,A00,A11,A01,A10)
!--------------------------------------------------------      
!
    use def_gparam
    implicit none
!
    integer, intent(in)  :: l
    real(DP),intent(in)  :: f_in
    real(DP),intent(out) ::  A00,A11,A01,A10
!
    integer :: i,j
    real(DP):: y11,y12,y13,y14,y21,y22,y23,y24,y5,y6,om,Ol,wt
    real(DP):: det(4),detp,fi,ddf
!
    if (.not.param_set) stop 'AA1lim: il faut passer par set_param avant!'
!
    ddf=df/nbdiv
!
    call wtcoef(f_in,fmin,fmin,fmax,fmax,wt) !filtre passe bas en cos tapper
    if (wt.ne.0.0_DP) then
       if (f_in.ne.0.0_DP) then
          Ol=sqrt(real(l*(l+1),DP)/2.0_DP)
!
! on calcule la derivee du determinant en f:
!
          do i=1,4
             if (i.le.2) then
                fi=f_in+(i-3)*ddf
             else
                fi=f_in+(i-2)*ddf
             endif
             det(i)=det1(l,fi)
          enddo
          detp=(-det(4)+8.0_DP*(det(3)-det(2))+det(1)) &
               /12.0_DP/ddf
!
          om=2.*pi*f_in
          call rayleigh(l,om,y11,y12,y13,y14,y21,y22 &
             ,y23,y24)    
          if ( l /= 0) then
             A00=   (y13*y22-y23*y12)       /detp
             A11=   (y14*y21-y24*y11)/2.0_DP/detp
             A10=   (y21*y12-y11*y22)/2.0_DP/detp/Ol
             A01=Ol*(y13*y24-y23*y14)       /detp   
          else
             A00=-y12/detp 
             A11=0.
             A10=0.
             A01=0.
          endif
       endif
    endif
!--------------------------------------------------------            
  end subroutine AA1lim
!--------------------------------------------------------      

!--------------------------------------------------------      
  subroutine  AA2lim(l,f_in,A)
!--------------------------------------------------------      
!
    use def_gparam  
    implicit none
!
    real(DP), intent(in)  :: f_in
    integer , intent(in)  :: l
    real(DP), intent(out) :: A
!
    integer :: i,j
    real(DP) :: y5,y6,om,Ol,det(4),detp,fi,ddf,wt
!
    if (.not.param_set) stop 'AA2lim: il faut passer par set_param avant!'
!
    ddf=df/nbdiv
!
    call wtcoef(f_in,fmin,fmin,fmax,fmax,wt) !filtre passe bas en cos tapper
    if (wt.ne.0.0_DP) then
       if (f_in.ne.0.0_DP.and.l.ne.0) then
          Ol=sqrt(real((l*(l+1))/2.,DP))
!
! on calcule la derivee du determinant en f:
!
          do i=1,4
             if (i.le.2) then
                fi=f_in+(i-3)*ddf
             else
                fi=f_in+(i-2)*ddf
             endif
             det(i)=det2(l,fi)
          enddo
          detp=(-det(4)+8.0_DP*(det(3)-det(2))+det(1)) &
                /12.0_DP/ddf
!
          om=2.*pi*f_in
          call love(l,om,y5,y6)
          if (y5.gt.1.0D-2) then
             print*,'warning dans AA2lim, det=',y5
          endif
          A=y6/2.0_DP/detp
       else
          A=0.
       endif
    endif
!--------------------------------------------------------            
  end subroutine AA2lim
!--------------------------------------------------------      

!--------------------------------------------------------      
  subroutine mat_construct(l,f_in,D,Dm1,T,det)
!--------------------------------------------------------      
    use def_gparam
    implicit none
!     

    integer, intent(in)  :: l
    real(DP),intent(in)  :: f_in
    real(DP), dimension(2), intent(out) :: det
    complex(DP), dimension(3,-1:1), intent(out) ::  D,T
    complex(DP), dimension(-1:1,3), intent(out) ::  Dm1
    integer :: i,j
    real(DP) y11,y12,y13,y14,y21,y22,y23,y24,y5,y6,om,Ol
!
    if (.not.param_set) stop 'mat_construct: il faut passer par &
                            &set_param avant!'
!
    if (f_in.ne.0.0_DP) then
       om=2.*pi*f_in
       call rayleigh(l,om,y11,y12,y13,y14,y21,y22 &
            ,y23,y24)     
       call love(l,om,y5,y6)       
!
       Ol=sqrt(real((l*(l+1))/2.,DP))
!
       det(1)=-(y11*y23-y13*y21)
       det(2)= y5
!     
       D(1,-1)=cmplx(Ol*y13,0.0   ,DP )
       D(2,-1)=cmplx(Ol*y23,0.0   ,DP )
       D(3,-1)=cmplx(0.0   ,-Ol*y5,DP )
       D(1, 0)=cmplx(y11   ,0.0   ,DP )
       D(2, 0)=cmplx(y21   ,0.0   ,DP )
       D(3, 0)=cmplx(0.0   ,0.0   ,DP )
       D(1, 1)= D(1,-1)
       D(2, 1)= D(2,-1)
       D(3, 1)=-D(3,-1)
!     
       Dm1(-1,1)=cmplx( y21/2.0/det(1)/Ol,0.0  , DP)
       Dm1( 0,1)=cmplx(-y23/det(1),0.0         , DP)
       Dm1( 1,1)= Dm1(-1,1)
       Dm1(-1,2)=cmplx(-y11/2.0/det(1)/Ol,0.0  , DP)
       Dm1( 0,2)=cmplx(y13/det(1),0.0          , DP)
       Dm1( 1,2)= Dm1(-1,2)
       Dm1(-1,3)=cmplx(0.0,.5/Ol/y5            , DP)
       Dm1( 0,3)=cmplx(0.0,0.0                 , DP)
       Dm1( 1,3)=-Dm1(-1,3)
!
       T(1,-1)=cmplx(Ol*y14, 0.0               , DP)
       T(2,-1)=cmplx(Ol*y24, 0.0               , DP)
       T(3,-1)=cmplx(0.0   ,-y6                , DP)
       T(1, 0)=cmplx(y12   , 0.0               , DP)
       T(2, 0)=cmplx(y22   , 0.0               , DP)
       T(3 ,0)=cmplx(0.0   , 0.0               , DP )
       T(1, 1)=cmplx(Ol*y14, 0.0               , DP)
       T(2, 1)=cmplx(Ol*y24, 0.0               , DP)
       T(3, 1)=cmplx(0.0   , y6                , DP)
    else
       D(:,:)=0.0_DP
       T(:,:)=0.0_DP
    endif
!
!----------------------------------------------------------------------
  end subroutine mat_construct
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  subroutine rayleigh(l,om,y11,y12,y13,y14,y21,y22,y23,y24)
!----------------------------------------------------------------------
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
    real(DP), intent(in) ::  om
    real(DP), intent(out)::  y11,y12,y13,y14,y21,y22,y23,y24
!
    integer l2
!
    real(DP) x,jl,jlp1,lp2mu,x2,f2mu,bidon,f2x,r2 
!
    if (.not.param_set) stop 'rayleigh: il faut passer par set_param avant!'
!
    x=om*r/vp 
    x2=x*x 
    f2x=2.0_DP*x
    f2mu=2.0_DP*mu
    lp2mu=lamb+f2mu 
    r2=r*r 
    !
    if (x.ne.0.0_DP) then
       call sphbes(l  ,x,jl  ,bidon,bidon,bidon)       
       call sphbes(l+1,x,jlp1,bidon,bidon,bidon) 
!     detection des NaN: 
       if (jl  *0.0_DP.ne.0.0_DP) then
          jl  =0.0_DP
!            print*,'NaN dans sphbes' 
       endif
       if (jlp1*0.0_DP.ne.0.0_DP) then
          jlp1=0.0_DP
!            print*,'NaN dans sphbes' 
       endif
!
! attention les sol sont multipliees par r par rapport a takeuchi & saito!
!
       if (l.ne.0) then
          y11=(dble(l)*jl-x*jlp1)/r 
          y12=(-lp2mu*x2*jl+f2mu*(dble(l*(l-1))*jl+f2x*jlp1))/r2 
          y13=jl/r 
          y14=f2mu*(dble(l-1)*jl-x*jlp1)/r2 
       else
          y11=-x*jlp1/r
          y12=-lp2mu*(x/r)**2*jl+4*mu/r2*x*jlp1
       endif
    else   
       y11=0.0_DP 
       y12=0.0_DP
       y13=0.0_DP  
       y14=0.0_DP 
    endif
!   
    x=om*r/vs 
    x2=x*x 
    f2x=2*x 
    l2=l*l 
!
    if (x.ne.0.0_DP.and.l.ne.0) then
       call sphbes(l  ,x,jl  ,bidon,bidon,bidon)      
       call sphbes(l+1,x,jlp1,bidon,bidon,bidon)      
!     detection des NaN:
       if (jl  *0.0_DP.ne.0.0_DP) then
          jl  =0.0_DP
!            print*,'NaN dans sphbes'
       endif
       if (jlp1*0.0_DP.ne.0.0_DP) then
          jlp1=0.0_DP
!            print*,'NaN dans sphbes'
       endif
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

!----------------------------------------------------------------------
  subroutine love(l,om,y1,y2)
!----------------------------------------------------------------------
! Calcul le solution toroidale (page 243 de TAKEUCHI et SAITO 1972)
! pour une boule homogenes.
! INPUT: 
! l entier
! r rayon en m, real(DP)
! om : pulsation
! vs vitesse des ondes s en m/s, real(DP)
! mu : mu en  kg m^5s^-4  , real(DP)
! OUTPUT:
! y: real(DP) , la solution de la page 243 de TAKEUCHI et SAITO 1972
!----------------------------------------------------------------------
    use def_gparam
    implicit none
!
    real(DP), intent(in) :: om
    integer , intent(in) :: l
    real(DP), intent(out):: y1,y2
!
    real(DP)  :: x,jl,jlp1,bidon
!
    if (.not.param_set) stop 'love: il faut passer par set_param avant!'
    x=om*r/vs
!
    if (x.ne.0.0_DP.and.l.ne.0) then
       call sphbes(l  ,x,jl  ,bidon,bidon,bidon)      
       call sphbes(l+1,x,jlp1,bidon,bidon,bidon)      
!    detection des NaN:
       if (jl  *0.0_DP.ne.0.0_DP) jl  =0.0_DP
       if (jlp1*0.0_DP.ne.0.0_DP) jlp1=0.0_DP
!    
       y1=jl
       y2=mu/r*(dble(l-1)*jl-x*jlp1)
    else
       y1=0.0_DP
       y2=0.0_DP
    endif
!----------------------------------------------------------------------
  end subroutine love
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  subroutine R_mode(l,f,u,du,v,dv,flag)
!----------------------------------------------------------------------
    use def_gparam
    implicit none
!
    integer, intent(in) :: l
    real(DP), intent(in) :: f
    real(DP), dimension(:), intent(out) :: u,du,v,dv
    logical :: flag
!
    integer :: i
    real(DP) :: om,y11,y12,y13,y14,y21,y22,y23,y24,a,u1,u2,v1,v2,du1,dv1 &
               ,du2,dv2,ray,norme
!
    om=2.0_DP*pi*f
    call rayleigh(l,om,y11,y12,y13,y14,y21,y22 &
         ,y23,y24) 
    if (abs(y22)>1.E-17_DP) then
       a=-y12/y22
    else
       a=0.0_DP
    endif 
!
    do i=1,NCM
       ray=(i*dr)
       call rayleigh2(l,om,ray,u1,du1,v1,dv1,u2,du2,v2,dv2) 
       u(i) = u1+a* u2
       v(i) = v1+a* v2
       du(i)=du1+a*du2
       dv(i)=dv1+a*dv2
    enddo
!
    if (.not.gauss_param) then
       allocate(r_gauss(NCM),r_gauss2(NCM),w_gauss(NCM))
       call gauleg(0._DP,r,r_gauss,w_gauss,NCM)
       r_gauss2(:)=r_gauss(:)**2
       gauss_param=.true.
    endif
    call R_norme(l,om,NCM,norme)
    if (abs(norme) > 1.E-50_DP) then
       u(:) =u(:) /norme
       du(:)=du(:)/norme*r
       v(:) =v(:) /norme
       dv(:)=dv(:)/norme*r
    endif
    flag=.false.
    
!----------------------------------------------------------------------
  end subroutine R_mode
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  subroutine R_norme(l,om,nbp,norme)
! nbp : nombre de point sur lesquels on fait l'integration
!----------------------------------------------------------------------
    implicit none
    integer, intent(in) :: l,nbp
    real(DP), intent(in) :: om
    real(DP), intent(out) :: norme
!tableaux automatiques:
    real(DP), dimension(nbp) :: u,du,v,dv    
!
    real(DP) :: y,y12,y22,y21,y11,a,u1,du1,v1,dv1,u2,du2,v2,dv2
    integer :: i,ll
!
    if (.not.gauss_param) stop 'R_norme: integration de gauss non initialisee'
!
!on cherche a annuler la traction  a la surface:
    call rayleigh(l,om,y11,y12,y,y,y21,y22,y,y) 
    if (abs(y22)>1.E-17_DP) then
       a=-y12/y22
    else
       a=0.0_DP
    endif
!    if (abs(y21)>1.E-17_DP) then
!       a=-y11/y21
!    else
!       a=0.0_DP
!    endif  
    do i=1,nbp
       call rayleigh2(l,om,r_gauss(i),u1,du1,v1,dv1,u2,du2,v2,dv2) 
       u(i) = u1+a* u2
       v(i) = v1+a* v2
       du(i)=du1+a*du2
       dv(i)=dv1+a*dv2
    enddo  
    if (l==0) then
       norme=SUM(u(:)**2*r_gauss2(:)*w_gauss(:))
       norme=sqrt(rho*norme)*scale
    else
       ll=real(l*(l+1),DP)
       norme=SUM((u(:)**2+real(l*(l+1),DP)*v(:)**2)*r_gauss2(:)*w_gauss(:))
       norme=sqrt(rho*norme)*scale
    endif  
    
!----------------------------------------------------------------------
  end subroutine R_norme
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  subroutine L_mode(l,f,w,dw,flag)
!----------------------------------------------------------------------
    use def_gparam
    implicit none
!
    integer, intent(in) :: l
    real(DP), intent(in) :: f
    real(DP), dimension(:), intent(out) :: w,dw
    real(DP) :: ray
    logical :: flag
!
    integer  :: i
    real(DP) :: om,norme
!
    om=2.d0*pi*f
    do i=1,NCM
       ray=(i*dr)
       call love2(l,om,ray,w(i),dw(i))
    enddo
    if (.not.gauss_param) then
       allocate(r_gauss(NCM),r_gauss2(NCM),w_gauss(NCM))
       call gauleg(0._DP,r,r_gauss,w_gauss,NCM)
       r_gauss2(:)=r_gauss(:)**2
       gauss_param=.true.
    endif
    call L_norme(l,om,NCM,norme)
    if (abs(norme) > 1.E-50_DP) then
       w(:) =w(:) /norme
       dw(:)=dw(:)/norme*r
    endif
    flag=.false.
!----------------------------------------------------------------------
  end subroutine L_mode
!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine L_norme(l,om,nbp,norme)
! nbp : nombre de point sur lesquels on fait l'integration
!----------------------------------------------------------------------
    implicit none
    integer, intent(in) :: l,nbp
    real(DP), intent(in) :: om
    real(DP), intent(out) :: norme
!tableaux automatiques:
    real(DP), dimension(nbp) :: w,dw
!
    integer :: i,ll
!
    if (.not.gauss_param) stop 'L_norme: integration de gauss non initialisee'
!
    if (l/=0) then
       do i=1,nbp
          call love2(l,om,r_gauss(i),w(i),dw(i))
       enddo
       norme=SUM(w(:)**2*r_gauss2(:)*w_gauss(:))
       ll=real(l*(l+1),DP)
       norme=sqrt(real(l*(l+1),DP)*rho*norme)*scale
    else
       norme=0.0_DP
    endif
!----------------------------------------------------------------------
  end subroutine L_norme
!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine norme_Lmode(l,w,dw,flag)
!----------------------------------------------------------------------
    use def_gparam
    implicit none
    integer, intent(in)  :: l
    real(DP), dimension(:), intent(inout) :: w,dw
    logical, intent(out) :: flag
!
    real(DP) :: int
    integer  :: i
!
    if (l/=0) then
       int=(w(1)**2*ray2(1)+w(NCM)**2*ray2(NCM))/2.0_DP
       do i=2,NCM-1
          int=int+w(i)**2*ray2(i)
       enddo
       int=sqrt(real(l*(l+1),DP)*rho*int*dr)*scale
    else
       int=0.0_DP
    endif
    print*,'int toute conne pour L:',int
    if (abs(int) > 1.E-50_DP) then
       w(:) =w(:) /int
       dw(:)=dw(:)/int*r
       flag=.false.
    else
       flag=.true.
    endif

!----------------------------------------------------------------------
  end subroutine norme_Lmode
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  subroutine norme_Rmode(l,u,du,v,dv,flag)
!----------------------------------------------------------------------
    use def_gparam
    implicit none
    integer, intent(in)  :: l
    real(DP), dimension(:), intent(inout) :: u,du,v,dv
    logical, intent(out) :: flag
!
    real(DP) :: int,ll
    integer :: i
!
    int=0
    if (l==0) then
       int=(u(1)**2*ray2(1)+u(NCM)**2*ray2(NCM))/2.0_DP
       do i=2,NCM-1
          int=int+u(i)**2*ray2(i)
       enddo
       int=sqrt(rho*int*dr)*scale
    else
       ll=real(l*(l+1),DP)
       int=((u(1  )**2+ll*v(1  )**2)*ray2(1  ) &
            +(u(NCM)**2+ll*v(NCM)**2)*ray2(NCM))
       do i=2,NCM-1
          int=int+(u(i)**2+real(l*(l+1),DP)*v(i)**2)*ray2(i)
       enddo
       int=sqrt(rho*int*dr)*scale
    endif
    print*,'Pour int toute conne int=',int
   if (abs(int) > 1.E-50_DP) then
       u(:) =u(:) /int
       du(:)=du(:)/int*r
       v(:) =v(:) /int
       dv(:)=dv(:)/int*r
       flag=.false.
    else
       flag=.true.
    endif
!----------------------------------------------------------------------
  end subroutine norme_Rmode
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  subroutine norme_Rmode_gauss(l,u,du,v,dv,flag)
!----------------------------------------------------------------------
    use def_gparam
    implicit none
    integer, intent(in)  :: l
    real(DP), dimension(:), intent(inout) :: u,du,v,dv
    logical, intent(out) :: flag
!
    real(DP) :: int,ll
    integer :: i
!
    if (l==0) then
       int=SUM(u(:)**2*r_gauss2(:)*w_gauss(:))
       int=sqrt(rho*int)*scale
    else
       ll=real(l*(l+1),DP)
       int=SUM((u(:)**2+real(l*(l+1),DP)*v(:)**2)*r_gauss2(:)*w_gauss(:))
       int=sqrt(rho*int)*scale
    endif
    print*,'Pour int de gauss int=',int
   if (abs(int) > 1.E-50_DP) then
       u(:) =u(:) /int
       du(:)=du(:)/int*r
       v(:) =v(:) /int
       dv(:)=dv(:)/int*r
       flag=.false.
    else
       flag=.true.
    endif

!----------------------------------------------------------------------
  end subroutine norme_Rmode_gauss
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  subroutine rayleigh2(l,om,r_in,y11,y12,y13,y14,y21,y22,y23,y24)
!----------------------------------------------------------------------
! meme chose que rayleigh, mais y12, y14, y22 et y24 ne sont plus les 
! tractions mais la derivees par rapport a r
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
    real(DP), intent(in) ::  om,r_in
    real(DP), intent(out)::  y11,y12,y13,y14,y21,y22,y23,y24
!
    integer l2
!
    real(DP) x,jl,jlp1,bidon,dx,djl,djlp1
!
    if (.not.param_set) stop 'rayleigh: il faut passer par set_param avant!'
!
    x=om*r_in/vp
    dx=om/vp
    !
    if (x.ne.0.0_DP) then
       call sphbes(l  ,x,jl  ,bidon,djl  ,bidon)       
       call sphbes(l+1,x,jlp1,bidon,djlp1,bidon) 
!     detection des NaN: 
       if (jl  *0.0_DP.ne.0.0_DP) then
          jl  =0.0_DP
!            print*,'NaN dans sphbes' 
       endif
       if (jlp1*0.0_DP.ne.0.0_DP) then
          jlp1=0.0_DP
!            print*,'NaN dans sphbes' 
       endif
!
! attention les sol sont multipliees par r par rapport a takeuchi & saito!
!
       if (l.ne.0) then
          y11=(dble(l)*jl-x*jlp1)/r_in 
          y12=(dble(l)*djl*dx-x*djlp1*dx-dx*jlp1)/r_in-y11/r_in
          y13=jl/r_in 
          y14=djl*dx/r_in -y13/r_in
       else
          y11=-x*jlp1/r_in
          y12=(-dx*jlp1-x*dx*djlp1)/r_in-y11/r_in
       endif
    else   
       y11=0.0_DP 
       y12=0.0_DP
       y13=0.0_DP  
       y14=0.0_DP 
    endif
!   
    x=om*r_in/vs 
    dx=om/vs
!
    if (x.ne.0.0_DP.and.l.ne.0) then
       call sphbes(l  ,x,jl  ,bidon,djl  ,bidon)      
       call sphbes(l+1,x,jlp1,bidon,djlp1,bidon)      
!     detection des NaN:
       if (jl  *0.0_DP.ne.0.0_DP) then
          jl  =0.0_DP
!            print*,'NaN dans sphbes'
       endif
       if (jlp1*0.0_DP.ne.0.0_DP) then
          jlp1=0.0_DP
!            print*,'NaN dans sphbes'
       endif
!     
       y21=-dble(l*(l+1))*jl/r_in
       y22=-dble(l*(l+1))*djl*dx/r_in-y21/r_in
       y23=(-dble((l+1))*jl+x*jlp1)/r_in
       y24=(-dble((l+1))*djl*dx+dx*jlp1+x*dx*djlp1)/r_in-y23/r_in
    else
       y21=0.0_DP
       y22=0.0_DP
       y23=0.0_DP
       y24=0.0_DP
    endif
!----------------------------------------------------------------------
  end subroutine rayleigh2
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  subroutine love2(l,om,r_in,y1,y2)
!----------------------------------------------------------------------
! meme chose que love, mais y12,n'est  plus la 
! tractions mains le derivees pa rapport a r
! Calcul le solution toroidale (page 243 de TAKEUCHI et SAITO 1972)
! pour une boule homogenes.
! INPUT: 
! l entier
! r rayon en m, real(DP)
! om : pulsation
! vs vitesse des ondes s en m/s, real(DP)
! mu : mu en  kg m^5s^-4  , real(DP)
! OUTPUT:
! y: real(DP) , la solution de la page 243 de TAKEUCHI et SAITO 1972
!----------------------------------------------------------------------
    use def_gparam
    implicit none
!
    real(DP), intent(in) :: om
    integer , intent(in) :: l
    real(DP), intent(in)::    r_in
    real(DP), intent(out):: y1,y2
!
    real(DP)  :: x,jl,bidon,djl,dx
!
    if (.not.param_set) stop 'love: il faut passer par set_param avant!'
    x=om*r_in/vs
    dx=om/vs
!
    if (x.ne.0.0_DP.and.l.ne.0) then
       call sphbes(l  ,x,jl  ,bidon,djl,bidon)      
!    detection des NaN:
       if (jl  *0.0_DP.ne.0.0_DP) jl  =0.0_DP
!    
       y1=jl
       y2=dx*djl
    else
       y1=0.0_DP
       y2=0.0_DP
    endif
!----------------------------------------------------------------------
  end subroutine love2
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  function filtre_coef(f) result(wt)
!----------------------------------------------------------------------
    use def_gparam
    implicit none
!    
    real(DP) :: f,wt
!
    if (.not.filtre_set) stop 'filtre_coef: les frequences du filtre &
                             &n''ont pas ete fixees!'
    call wtcoef(f,f_filtre_1,f_filtre_2,f_filtre_3,f_filtre_4,wt)
!   
!----------------------------------------------------------------------
  end function filtre_coef
!----------------------------------------------------------------------'

!----------------------------------------------------------------------
  subroutine wtcoef(f,f1,f2,f3,f4,wt)
!----------------------------------------------------------------------
    use def_gparam
    implicit none
!
    real(DP), intent(in) ::  f,f1,f2,f3,f4
    real(DP), intent(out)::  wt
!
    if (f3.gt.f4) stop 'wtcoef: f3>f4 '
    if (f1.gt.f2) stop 'wtcoef: f1>f2 '
    if (f.ge.f4.or.f.lt.f1 ) wt=0.0_DP
    if (f.le.f3.and.f.ge.f2) wt=1.0_DP
    if (f.gt.f3.and.f.lt.f4) wt=0.5_DP*(1.0+cos(pi*(f-f3)/(f4-f3)))
    if (f.gt.f1.and.f.lt.f2) wt=0.5_DP*(1.0+cos(pi*(f-f2)/(f2-f1)))
!----------------------------------------------------------------------
  end subroutine wtcoef
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  subroutine set_filtre(f1_in,f2_in,f3_in,f4_in,lmax_freq)
!----------------------------------------------------------------------
    use def_gparam
    implicit none
!
    real(DP),intent(in) :: f1_in,f2_in,f3_in,f4_in
    integer ,optional, intent(out):: lmax_freq

    integer :: i,l,n
!
    f_filtre_1=f1_in
    f_filtre_2=f2_in
    f_filtre_3=f3_in
    f_filtre_4=f4_in

    lmax_freq=-1
    if (present(lmax_freq)) then
       if (.not.freq_read) then
          stop 'set_filtre: je ne peux pas &
               &determiner le lmax car les frequences n''ont ete lues'
       else
          do i=1,2
             do l=0,lmax
                do n=1,n_zero(l,i)
                   if (f(n,l,i) < f_filtre_4 .and. lmax_freq < l ) lmax_freq=l
                enddo
             enddo
          enddo
       endif
       if (mod(lmax_freq,2) == 0 )  lmax_freq=lmax_freq+1
       lmax=lmax_freq
       filtre_set=.true.
    else
       print*,'set_filtre: Attention, vous n''avez pas recalcule lmax,&
              & celui utilise risque d''etre trop grand'
    endif  
    if (freq_read) then
       if (f_filtre_1 < fmin) print*,'set_filtre: Attention, la &
                             &frequence min du filtre est < fmin presente'
       if (f_filtre_4 > fmax) print*,'set_filtre: Attention, la &
                             &frequence max du filtre est > fmax presente'
    endif
!    
!----------------------------------------------------------------------
  end subroutine set_filtre
!----------------------------------------------------------------------


!----------------------------------------------------------------------
subroutine write_fctp(f,l,n,type)
!----------------------------------------------------------------------
  use def_gparam
  implicit none
!
  real(DP), intent(in) :: f
  integer, intent(in)  :: l,n,type
  real(DP), dimension(NCM) :: u,v,du,dv,w,dw
  integer :: i
  real(DP) :: bid=0.0_DP
  logical :: flag
!
  select case(type)
  case(1)
!rayleigh:
     call R_mode(l,f,u,du,v,dv,flag)
     if (.not.flag) then
        if (oldminos) then
           write(unit_rd) n-1,l,2.0_DP*pi*f,bid,bid,bid,&
                (real( u(i),SP),i=1,NCM)  &
                ,(real(du(i),SP),i=1,NCM)  &
                ,(real( v(i),SP),i=1,NCM)  &
                ,(real(dv(i),SP),i=1,NCM)  
        else
           rayleigh_rec=rayleigh_rec+1
           write(unit_rd,rec=rayleigh_rec) n-1,l,2.0_DP*pi*f,0.0_DP,0.0_DP,&
                (real(u(i),SP),real(du(i),SP),real(v(i),SP),real(dv(i),SP),i=NCM,1,-1)  
           write(unit_ri,*) n-1,l,2.0_DP*pi*f,0.0_DP,0.0_DP,rayleigh_rec,0.d0,.false.
        endif
     endif
!test
!     if (n==1) then
!        print*,'ecriture de ',l
!        do i=1,NCM
!           write(l,*)i,u(i)
!        enddo
!     endif
!fin test
  case(2)
!love:
     call L_mode(l,f,w,dw,flag)
     if (.not.flag) then
        if (oldminos) then
           write(unit_ld) n-1,l,2.0_DP*pi*f,bid,bid,bid, &
                (real( w(i),SP),i=1,NCM)  &
                ,(real(dw(i),SP),i=1,NCM)  
        else
           love_rec=love_rec+1
           write(unit_ld,rec=love_rec) n-1,l,2.0_DP*pi*f,0.0_DP,0.0_DP,&
                (real(w(i),SP),real(dw(i),SP),i=NCM,1,-1)
           write(unit_li,*) n-1,l,2.0_DP*pi*f,0.0_DP,0.0_DP,love_rec,0.d0,.false.
        endif
     endif
  case default
     stop 'write_fctp: type doit valoir 1 (rayleigh) ou 2(love)'
  end select

!----------------------------------------------------------------------
end subroutine write_fctp
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine close_boule()
!----------------------------------------------------------------------
  use def_gparam
  implicit none
  integer :: i,l,n
  if (.not.oldminos) then 
  write(unit_rd,REC=1)len_r,lmax_pres_r,rayleigh_rec,0,NCM &
                    ,(real(1.0_DP-(i-1)*dr/r,SP),i=NCM,1,-1)
  write(unit_ld,REC=1)len_l,lmax_pres_l,love_rec    ,0,NCM &
                    ,(real(1.0_DP-(i-1)*dr/r,SP),i=NCM,1,-1)
  endif
  open(23,file='freqR')
  do l=0,Lmax_pres_r
     do n=1,n_zero(l,1)
        write(23,*)l,n-1,f(n,l,1)
     enddo
  enddo
  close(23)
  open(23,file='freqL')
  do l=0,Lmax_pres_l
     do n=1,n_zero(l,2)
        write(23,*)l,n-1,f(n,l,2)
     enddo
  enddo
  close(23)
  close(unit_rd)
  close(unit_ri)
  close(unit_ld)
  close(unit_li)
!----------------------------------------------------------------------
end subroutine close_boule
!----------------------------------------------------------------------
!---------------------------------------------------------------------
      SUBROUTINE extension(fichier1,fichier2,ext) 
!---------------------------------------------------------------------
!
      IMPLICIT NONE
!
      character(len=*), intent(in) ::  fichier1,ext
      character(len=*), intent(out):: fichier2
      integer lenfichier1,lenext
!
      lenfichier1=INDEX (fichier1,' ') -1  
      lenext=INDEX (ext,' ') -1  
      if ((lenfichier1+lenext).gt.80) then
         print*,'il y a trop de caracter a ',fichier1
      endif
!
      fichier2(1:lenfichier1)=fichier1(1:lenfichier1)
      fichier2(lenfichier1+1:(lenfichier1+lenext))=ext(1:lenext)
!---------------------------------------------------------------------
end SUBROUTINE extension
!---------------------------------------------------------------------
!----------------------------------------------------------------------
end module sol_ana
!----------------------------------------------------------------------






