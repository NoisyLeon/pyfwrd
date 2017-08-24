!------------------------------------------------------
module global_main_param
!------------------------------------------------------
  use def_gparam
  implicit none
  public:: recepteur_dat  &
         ,nms_dat,NBT,dt,t0,f1,f2,f3,f4,RA,f0    &
         ,eigenfileT,eigenfileS,amp_source   &
         ,get_NBE,init_global_main,GEOC,source_type &
         ,geocentric_correction,nbcomp,NBR,NBE,stations_name      &
         ,get_source_dat,sources_name,logunit,euler&
         ,coord_stations,get_receivers,load_source     &
         ,set_logunit_in_global_main,t_delay_sources &
         ,SEXP,M0,source_delay,comp,compr    &
         ,local_backup_dir,reset_sources_delay,get_Mnt &
         ,load_receivers,get_NBR,rotation ,stacking_sources,nmax_in,nmin
  private
!configuration files:
  character(len=30)            :: source_dat         ='sources.dat'
  character(len=30)            :: recepteur_dat      ='receivers.dat'
  character(len=30), parameter :: nms_dat            ='nms.dat'
  character(len=1), dimension(3), parameter  :: compr=(/'Z','R','T'/)
  character(len=1), dimension(3), parameter  :: comp=(/'Z','N','E'/)
!sources, stations, 
  integer :: NBE
  integer :: NBR
  integer, parameter :: nbcomp=3!nombre de composante (1== veritcal, 3 == all)
  integer :: iNNmax
  real(DP), dimension(:,:), allocatable :: coord_stations,Mtmp
  real(DP), dimension(:,:), allocatable :: coord_sources
  real(DP), dimension(:  ), allocatable :: t_delay_sources
  logical :: geoc_corr,rotation
  character(len=8), dimension(: ), allocatable :: sources_name
  logical :: sources_available=.false.,stacking_sources,source_delay
  doubleprecision, parameter :: SEXP=1.d18
  character(len=4), dimension(: ), allocatable :: stations_name
!temps:
  integer :: NBT
  real(DP):: dt,t0,amp_source
  real(DP):: f1,f2,f3,f4,f0 
  character(len=6), parameter :: source_type='heavis'
  logical :: local_backup
  character(len=100)  :: local_backup_dir
!frequence
  integer  :: NBF,NBF2 
   real(DP):: dtf,duree,deltaf,dtf2
  character(len=100) :: eigenfileT,eigenfileS
!
  real(DP), parameter :: RA=6371000.0, GEOC=0.993277d0
  integer :: nmax_in,nmin
!
!logical unit for log file:
  integer :: logunit
!pour les FFT
  integer :: coef_nbt
!
contains 
!---------------------------------------------------
  subroutine init_global_main(logunit_)
!---------------------------------------------------
    implicit none
    integer, intent(in) :: logunit_
!
    logunit=logunit_
!
    call load_nms()
    call init_rcv_src()
!---------------------------------------------------
  end subroutine init_global_main
!---------------------------------------------------
!---------------------------------------------------
  subroutine set_logunit_in_global_main(log)
!---------------------------------------------------
    implicit none
    integer, intent(in) :: log
    logunit=log
!---------------------------------------------------
  end subroutine set_logunit_in_global_main
!---------------------------------------------------
!-------------------------------------------------------------------------
    logical function geocentric_correction()
!-------------------------------------------------------------------------
      implicit none
      geocentric_correction=geoc_corr
!-------------------------------------------------------------------------
    end function geocentric_correction
!-------------------------------------------------------------------------
!---------------------------------------------------
  subroutine load_nms()
!---------------------------------------------------
    implicit none
    integer, parameter:: unit=13
    integer :: i
    
!
!    nbcomp=1
    open(unit,file=nms_dat,status='old')
    read(unit,*)
    read(unit,*) dt
    read(unit,*)
    read(unit,*) NBT
    read(unit,*)
    read(unit,*) t0
    read(unit,*)
    read(unit,*) amp_source
    read(unit,*) 
    read(unit,*) f1,f2,f3,f4
    read(unit,*)
    read(unit,*) nmin,nmax_in
    if (nmin<0) nmin=0
    if (nmax_in<0) nmax_in=-1
    read(unit,*)
    read(unit,*) geoc_corr
    read(unit,*)
    read(unit,*) rotation
    read(unit,*)
    read(unit,'(a)') eigenfileS
    read(unit,*)
    read(unit,'(a)') eigenfileT
    write(logunit,*) "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    write(logunit,*) "time step:",dt
    write(logunit,*) "Number of time step:",NBT
    write(logunit,*) "Source t0 time:",t0
    write(logunit,*) "source frequency band:",sngl(f1),sngl(f2),sngl(f3),sngl(f4)
    write(logunit,*) "Geocentric correction?:",geoc_corr 
    write(logunit,*) "Egeinfunctions files:"
    write(logunit,'(a)') eigenfileS
    write(logunit,'(a)') eigenfileT
!---------------------------------------------------
  end subroutine load_nms
!---------------------------------------------------

!---------------------------------------------------
  subroutine get_NBE(NBE_)
!---------------------------------------------------
    implicit none
    integer, intent(out) :: NBE_
    NBE_=NBE
!---------------------------------------------------
  end subroutine get_NBE
!---------------------------------------------------
!---------------------------------------------------
  subroutine get_NBR(NBR_)
!---------------------------------------------------
    implicit none
    integer, intent(out) :: NBR_
    NBR_=NBR
!---------------------------------------------------
  end subroutine get_NBR
!---------------------------------------------------
!============SOURCES:
!-------------------------------------------------------------------------
  subroutine init_rcv_src()
!-------------------------------------------------------------------------
    implicit none
    logical :: err
    call load_source()
    call load_receivers(err)
    if (err) STOP 'Error when loading the receivers file'
!-------------------------------------------------------------------------
  end subroutine init_rcv_src
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
  subroutine load_source()
!-------------------------------------------------------------------------
    implicit none
    character(len=8) :: name
    integer     :: i
!pour les corrections geocentrique 
    character(len=100) :: source_file
    integer :: units=527,j,ibid
!
!lecture du fichier de source:
    open(units,file=source_dat,status='old')
    read(units,*) 
    read(units,*) NBE
    read(units,*) 
    read(units,*) stacking_sources
    read(units,*) 
    allocate(coord_sources(3,NBE),Mtmp(6,NBE),sources_name(NBE),t_delay_sources(NBE)) 
    source_delay=.false.
    do i=1,NBE
        read(units,'(i2.2,1x,a8,1x,6(e9.2,1x),3(f6.2,1x),f7.2)') ibid,name, &
           (Mtmp (j,i),j=1,6),(coord_sources(j,i),j=1,3),t_delay_sources(i)
        sources_name(i)=name
        if (abs(t_delay_sources(i)) > 1.E-4_DP)  source_delay=.true.
!conversion latitude->colatitude      
       coord_sources(2,i)=90.d0-coord_sources(2,i)
!conversion profondeur (en km)->rayon (en m)
       coord_sources(1,i)=RA-coord_sources(1,i)*1000.d0
    enddo
    close(units)
    write(logunit,*) "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    write(logunit,*) "Number of sources to used:",NBE
    if (source_delay)  then
       write(logunit,*) " Sources delay mode (not all sources will triggered at t=0)"
    else
       write(logunit,*) " All triggered at t=0"
    endif
    if (rotation.and.stacking_sources) then
       print*,'Warning, can''t rotate component when stacking sources'
       rotation=.false.
    endif
!
!-------------------------------------------------------------------------
  end subroutine load_source
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
  subroutine get_Mnt(Mout,is)
!-------------------------------------------------------------------------
    implicit none
    real, dimension(:), intent(out):: Mout
    integer, intent(in) :: is
    Mout(:)=real(Mtmp(:,is))
!-------------------------------------------------------------------------
  end subroutine get_Mnt
!-------------------------------------------------------------------------


!-------------------------------------------------------------------------
  subroutine load_receivers(error)
!-------------------------------------------------------------------------
    implicit none
    logical, intent(out) :: error
    integer :: irec,unitrec
!    
!=========================================
!Lecture des donnees dans recepteur.dat
!=========================================
    unitrec=135
    open(unitrec, file=recepteur_dat,status='old',err=100)
    read(unitrec,*)      
    read(unitrec,*) NBR
    read(unitrec,*)      
    allocate(stations_name(NBR),coord_stations(2,NBR))
    do irec= 1,NBR
       read(unitrec,'(a4,1x,f8.3,1x,f8.3)') stations_name(irec) &
                ,coord_stations(1,irec),coord_stations(2,irec)
!latitude-> colatitude
       coord_stations(1,irec)=90.d0-coord_stations(1,irec)
    enddo

!Fermeture du fichier recepteur.dat
    close(unitrec)
    error=.false.
    write(logunit,*) "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    write(logunit,*) "Number of receivers:",NBR
    return
100 print*,'No recepteur.dat file found! ... skeeping receivers management..'     
    error=.true.
!-------------------------------------------------------------------------
  end subroutine load_receivers
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
  subroutine write_receivers(weight)
!-------------------------------------------------------------------------
    implicit none
    integer :: irec,unitrec
    real(DP), dimension(:), optional :: weight
!    
!=========================================
!Ecriture des donnees dans recepteur_out.dat
!=========================================
    unitrec=135
    open(unitrec, file='recepteur_out.dat')
    write(unitrec,'(a)')'#nombre de recepteurs:'
    write(unitrec,*) NBR
    write(unitrec,'(a)') '#nom (4 characters), latitude, longitude en degres &
                     &, 1/poids de la station'
    if (.not.present(weight) ) then
       do irec= 1,NBR
          write(unitrec,'(a4,1x,f8.3,1x,f8.3)') stations_name(irec) &
               ,90.d0-coord_stations(1,irec),coord_stations(2,irec)
       enddo
    else
       do irec= 1,NBR
          write(unitrec,'(a4,1x,f8.3,1x,f8.3,1x,f5.2)') stations_name(irec) &
               ,90.d0-coord_stations(1,irec),coord_stations(2,irec),weight(irec)
       enddo
    endif

!Fermeture du fichier recepteur_out.dat
    close(unitrec)
!cleaning
!-------------------------------------------------------------------------
  end subroutine write_receivers
!-------------------------------------------------------------------------
!-------------------------------------------------
  subroutine get_receivers(rec_coord)
!-------------------------------------------------
    use def_gparam
    implicit none
    real(DP), dimension(:,:), intent(out)  :: rec_coord
!
    rec_coord(:,:)=coord_stations(:,:)
    if (geoc_corr)                             &
         rec_coord(1,:)=rad2deg*(PI/2.d0-atan(GEOC*tan((PI/2.d0-  &
         rec_coord(1,:)/rad2deg))))
!on transforme les degres en radian
    rec_coord(:,:)=rec_coord(:,:)*deg2rad
!------------------------------------------------------------------------
  end subroutine get_receivers
!------------------------------------------------------------------------


!-------------------------------------------------------------------------
  doubleprecision function M0(M)
!calcul la norme de la matrice 3x3 symetrique
!| M(1) M(4) M(5) |
!|      M(2) M(6) |
!|           M(3) |
!-------------------------------------------------------------------------
    implicit none
    doubleprecision, dimension(6) ::M
    M0=sqrt(M(1)**2+M(2)**2+M(3)**2+2.d0*(M(4)**2+M(5)**2+M(6)**2))
!-------------------------------------------------------------------------
  end function M0
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    subroutine get_source_dat(coord_out,Mtmp_out,name_out)
!-------------------------------------------------------------------------
      implicit none
      doubleprecision, dimension(3,NBE), intent(out)  :: coord_out
      doubleprecision, dimension(6,NBE), intent(out)  :: Mtmp_out
      character(len=8),dimension(NBE), optional, intent(out) :: name_out
      integer :: i,j
!
      coord_out=coord_sources
      Mtmp_out =Mtmp
      if (geoc_corr)                             &
           coord_out(2,:)=rad2deg*(PI/2.d0-atan(GEOC*tan((PI/2.d0-  &
           coord_out(2,:)/rad2deg))))      
!passage en radian
      coord_out(2:3,:)=coord_out(2:3,:)*deg2rad
      if (present(name_out)) name_out(:)=sources_name(:)
!
!-------------------------------------------------------------------------
    end subroutine get_source_dat
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    subroutine reset_sources_delay
!-------------------------------------------------------------------------
      use def_gparam
      implicit none
      write(logunit,*)'Warning, forcing source delays to zero !!!!'
      t_delay_sources(:)=0.0_DP
!-------------------------------------------------------------------------
    end subroutine reset_sources_delay
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
  real(DP) function dist3D(pt1,pt2)
!pseudo distance dans le volume en radian entre pt1 et pt2
!-------------------------------------------------------------------------
    use def_gparam
    implicit none
    real(DP), dimension(3), intent(in) :: pt1,pt2
!
    real(DP) :: a,b,g
!
    call euler(pt1(2),pt1(3),pt2(2),pt2(3),a,b,g)
    dist3D=sqrt(b**2+((pt1(1)-pt2(1))/RA)**2)
!-------------------------------------------------------------------------
  end function dist3D
!-------------------------------------------------------------------------

      subroutine euler(t1,p1,t2,p2,a,b,g)
!_____________________________________________________________________
!   
!     calcule les 3 angles d'euler de la rotation definie par:
!     D(a b g)=D(-p1 -t1 0)oD(0 t2 p2)
!     (notation cf Edmonds)
!
!     Entrees:
!        angles p1 p2 dans [-pi,pi] (en fait peu importe mais en
!                                    radians)
!        angles t1 t2 DANS [ 0 ,pi]! (pour etre avec utilise avec 
!        la regle  de sommation d'harmoniques spheriques)
!     Sorties:
!        angles a et g dans [-pi,pi]
!        angle  b      dans [ 0 ,pi]
!_____________________________________________________________________
!
      IMPLICIT NONE
!
      real*8 t1,t2,p1,p2,a,b,g,ca,sa,cg,sg,sb,cb,pi
      real*8 rt1(3,3),rt2(3,3),rp1(3,3),rp2(3,3),prod1(3,3), &
           prod2(3,3),rot(3,3)
      integer i,j,k
!
      pi = 3.141592653589793d0
      do i=1,3
         do j=1,3
            rt1(i,j)  =0.0d0
            rt2(i,j)  =0.0d0
            rp1(i,j)  =0.0d0
            rp2(i,j)  =0.0d0
            prod1(i,j)=0.0d0
            prod2(i,j)=0.0d0
            rot(i,j)  =0.0d0
         enddo
      enddo
!      
      rt1(1,1)=dcos(t1)
      rt1(2,2)=1.0d0
      rt1(3,3)=rt1(1,1)
      rt1(1,3)=dsin(t1)
      rt1(3,1)=-rt1(1,3)
!
      rt2(1,1)=dcos(t2)
      rt2(2,2)=1.0d0
      rt2(3,3)=rt2(1,1)
      rt2(1,3)=-dsin(t2)
      rt2(3,1)=-rt2(1,3)
!
      rp1(1,1)=dcos(p1)
      rp1(2,2)=rp1(1,1)
      rp1(3,3)=1.0d0
      rp1(1,2)=-dsin(p1)
      rp1(2,1)=-rp1(1,2)
!
      rp2(1,1)=dcos(p2)
      rp2(2,2)=rp2(1,1)
      rp2(3,3)=1.0d0
      rp2(1,2)=dsin(p2)
      rp2(2,1)=-rp2(1,2)
!
      do i=1,3
         do j=1,3
            do k=1,3
               prod1(i,j)=prod1(i,j)+rp2(i,k)*rt2(k,j)
            enddo
         enddo
      enddo
      do i=1,3
         do j=1,3
            do k=1,3
               prod2(i,j)=prod2(i,j)+rp1(i,k)*prod1(k,j)
            enddo
         enddo
      enddo
      do i=1,3
         do j=1,3
            do k=1,3
               rot(i,j)=rot(i,j)+rt1(i,k)*prod2(k,j)
            enddo
         enddo
      enddo
!
      cb=rot(3,3)
      if (cb>=1.d0) then
         b=0.d0
      else if (cb<=-1.d0) then
         b=pi
      else
         b=dacos(cb)
      endif
      sb=dsin(b)
      if (abs(sb).le.1.0d-15) then
         a=p2-p1
         g=0.
      else 
         ca=rot(3,1)/sb
         sa=rot(3,2)/sb
         if (abs(ca-1.0d0).lt.1.0d-8) then
            a=0.0d0
         else 
            if (abs(ca+1.0d0).lt.1.0d-8) then
               a=pi
            else
               a=dacos(ca)
            endif
         endif
         if (sa.lt.(0.0d0)) a=-1.0d0*a
         cg=-rot(1,3)/sb
         sg=rot(2,3)/sb
         if (abs(cg-1.0d0).lt.1.0d-8) then
            g=0.0d0
         else
            if (abs(cg+1.0d0).lt.1.0d-8) then
               g=pi
            else
               g=dacos(cg)
            endif
         endif
         if (sg.lt.(0.0d0)) g=-1.0d0*g
      endif
!
      end subroutine euler
!-------------------------------------------------------------------------------
!------------------------------------------------------
end module global_main_param
!------------------------------------------------------
