!25/11/02: passage en F90 pour scat_sem .... Yahooo
!08/04/97:
!nouvelle version real*4 pour ww, q et g!mais real*4 pour les fctps
!    .(utilisable ave!les sorties de minosa3!)
!21/03/97 introduction des parametres NCM et NCM6
!la lecture des fctp est un peu bizard pour gagner du temps
!
!----------------------------------------------------------------------
module util_fctp
!----------------------------------------------------------------------
  use def_gparam
  use earth_modele
  implicit none
  public :: open_fctp,readS,readT,read_interpS,read_interpT,get_Lmax,get_Nmax &
           ,get_eigenw,get_1Dmodel,mode_existe,get_fmaxc,get_tref,euler,close_fctp &
           ,get_1Dmodelb,write_modeT_ascii, write_modeS_ascii,open_fctp1
  private
  integer   :: LLmaxS,NNmaxS,LLmaxT,NNmaxT,nkS,nkT,nb_discT,nb_discS
  integer,parameter   :: unitfctpS=155,unitfctpT=154
  real(DP), dimension(:,:), allocatable :: omtabS,qtabS,gctabS
  integer , dimension(:,:), allocatable :: irecS
  real(DP), dimension(:)  , allocatable :: rayS,discS_rad,discT_rad
  integer , dimension(:)  , allocatable :: discS_index,discT_index
  real(DP), dimension(:,:), allocatable :: omtabT,qtabT,gctabT
  integer , dimension(:,:), allocatable :: irecT
  real(DP), dimension(:)  , allocatable :: rayT
  type(modele) :: modele1D
!
  logical :: files_open=.false.,modele1D_read=.false. &
            ,files_openS=.false.,files_openT=.false.
!====================================================================
contains
!====================================================================
!-------------------------------------------------------------------- 
      SUBROUTINE open_fctp(fichierfctpS,fichierfctpT,model1Dfile,rad_recepteur)
!
!
!     fichierfctp: character*90, nom du fichier de fonction propre
!                  sans extension
!     unitfctp   : entier, unite logique du fichier de fctp desiree
!
!    irec(0:LLmax,0:NNmax) :tableaux d'entier des numeros de record
!                           des foctions propres
!    omtab(0:LLmax,0:NNmax):tableau de real*8,  w en rad/s
!    qtab(0:LLmax,0:NNmax) :tableau de real*8, attenuations
!    gctab(0:LLmax,0:NNmax):tableau de real*8, vitesses de groupe
!    nocean                :entier, nb de couches de l'ocean
!    nk                    :entier , nombre de couches du model
!    ray(NCM)              :tableau de real*4, couches du model
!--------------------------------------------------------------------
        use global_main_param, only: RA
        use earth_modele
        IMPLICIT NONE
!
      character(len=90), intent(in) :: fichierfctpS,fichierfctpT
      character(len=90), optional, intent(in) :: model1Dfile
      real(DP), optional, intent(out) :: rad_recepteur
!
      character(len=90) :: fichierdirect,fichierinfo,ext
      real(DP) :: ww,qq,gc,wdiff      
      integer  ::  LEN,unitinfo,indicepointeur,n,l,ibid,i,nocean,ier,k
      logical ::  flag
      real(SP), dimension(:)  , allocatable :: ray_tmp,ray_tmp2
!
!on Lit d'abord le fichier Spheroidaux
!
	unitinfo=73
      ext='.direct '
      call extension(fichierfctpS,fichierdirect,ext)
      ext='.info'
      call extension(fichierfctpS,fichierinfo,'.info ')
!
      inquire(file=fichierdirect,EXIST=flag)
      if (.not.flag) then
         print*,'le fichier',fichierdirect,'n''existe pas???'
         stop
      endif
      inquire(file=fichierinfo,EXIST=flag)
      if (.not.flag) then
         print*,'le fichier',fichierinfo,'n''existe pas???'
         stop
      endif
      call get_len(fichierdirect,LEN)
      open(unitfctpS,file=fichierdirect,status='old',access='direct' &
               ,recl=LEN)
      open(unitinfo,file=fichierinfo)
      read(unitfctpS,rec=1)ibid,ibid,ibid,nocean,nkS
      allocate(rayS(nkS),ray_tmp(nkS))
!reading raius and track for discontinuities
      read(unitfctpS,rec=1)ibid,ibid,ibid,nocean,nkS,(ray_tmp(i),i=1,nkS)
!km SP ->m DP
      rayS(:)=dble(ray_tmp(1:nkS))*RA 
      if (present(rad_recepteur)) rad_recepteur=rayS(nkS-nocean)-0.1
      nb_discS=0
      do i=1,nkS-1
         if (abs(rayS(i)-rayS(i+1))<1e-8)  nb_discS=nb_discS+1
      enddo  
      allocate(discS_index(0:nb_discS+1),discS_rad(0:nb_discS+1))
      discS_index(0)=0
      discS_index(nb_discS+1)=nkS
      discS_rad(0)  =rayS(1)
      discS_rad(nb_discS+1) =rayS(nkS)
      k=0
      do i=1,nkS-1
         if (abs(rayS(i)-rayS(i+1))<1e-8) then
            k=k+1
            discS_index(k)=i
            discS_rad  (k)=rayS(i)
         endif
      enddo  
      deallocate(ray_tmp)
!scaning for LLmax and NNmax:
      LLmaxS=0;NNmaxS=0
 10   read(unitinfo,*,END=99) n,l,ww,qq,gc,indicepointeur,wdiff,flag
         if (l>LLmaxS) LLmaxS=l
         if (n>NNmaxS) NNmaxS=n
         goto 10 
 99   continue
      allocate(irecS(0:LLmaxS,0:NNmaxS),omtabS(0:LLmaxS,0:NNmaxS)  &
              ,qtabS(0:LLmaxS,0:NNmaxS),gctabS(0:LLmaxS,0:NNmaxS), stat=ier  )
      irecS(:,:)=0
      if (ier/=0) then
         print*,'LLmaxS=',LLmaxS
         print*,'NNmaxS=',NNmaxS
         print*,(4.+8.*3.)*dble(LLmaxS+1)*dble(NNmaxS+1)/1024./1024.,' Mb'
         stop 'open_fctp: probleme d''allocation memoire'
      endif
      rewind(unitinfo)
!relecture
 11   read(unitinfo,*,END=98) n,l,ww,qq,gc,indicepointeur,wdiff,flag
         irecS (l,n)=indicepointeur
         omtabS(l,n)=ww
         qtabS (l,n)=qq
         gctabS(l,n)=gc
         if ((abs(wdiff).gt.1.d-9.and.flag).or.abs(wdiff).gt.1d-5) irecS(l,n)=0
      goto 11 
 98   continue     
      close(unitinfo)
      fichierdirect(:)=' '
      fichierinfo(:)=' '

!
!on Lit ensuite le fichier Toroidaux
!
      call extension(fichierfctpT,fichierdirect,'.direct ')
      call extension(fichierfctpT,fichierinfo,'.info ')
!
      inquire(file=fichierdirect,EXIST=flag)
      if (.not.flag) then
         print*,'le fichier',fichierdirect,'n''existe pas???'
         stop
      endif
      inquire(file=fichierinfo,EXIST=flag)
      if (.not.flag) then
         print*,'le fichier',fichierinfo,'n''existe pas???'
         stop
      endif
      call get_len(fichierdirect,LEN)
      open(unitfctpT,file=fichierdirect,status='old',access='direct' &
                ,recl=LEN)
      open(unitinfo,file=fichierinfo)
      read(unitfctpT,rec=1)ibid,ibid,ibid,nocean,nkT
      allocate(rayT(nkT),ray_tmp(nkT))
!reading raius and track for discontinuities
      read(unitfctpT,rec=1)ibid,ibid,ibid,nocean,nkT,(ray_tmp(i),i=1,nkT)
!km SP ->m DP
      rayT(:)=dble(ray_tmp(1:nkT))*RA
      nb_discT=0
      do i=1,nkT-1
         if (abs(rayT(i)-rayT(i+1))<1e-8)  nb_discT=nb_discT+1
      enddo  
      allocate(discT_index(0:nb_discT+1),discT_rad(0:nb_discT+1))
      discT_index(0)=0
      discT_index(nb_discT+1)=nkT
      discT_rad(0)  =rayT(1)
      discT_rad(nb_discT+1) =rayT(nkT)
      k=0
      do i=1,nkT-1
         if (abs(rayT(i)-rayT(i+1))<1e-8) then
            k=k+1
            discT_index(k)=i
            discT_rad  (k)=rayT(i)
         endif
      enddo  
      deallocate(ray_tmp)
!scaning for LLmax and NNmax:
      LLmaxT=0;NNmaxT=0
 110   read(unitinfo,*,END=199) n,l,ww,qq,gc,indicepointeur,wdiff,flag
         if (l>LLmaxT) LLmaxT=l
         if (n>NNmaxT) NNmaxT=n
         goto 110 
 199   continue
      allocate(irecT(0:LLmaxT,0:NNmaxT),omtabT(0:LLmaxT,0:NNmaxT))
      allocate(qtabT(0:LLmaxT,0:NNmaxT),gctabT(0:LLmaxT,0:NNmaxT))
      irecT(:,:)=0 
      rewind(unitinfo)
!relecture
 111   read(unitinfo,*,END=198) n,l,ww,qq,gc,indicepointeur,wdiff,flag
         irecT (l,n)=indicepointeur
         omtabT(l,n)=ww
         qtabT (l,n)=qq
         gctabT(l,n)=gc
!         if ((abs(wdiff).gt.1.d-9.and.flag).or.abs(wdiff).gt.1d-5) irecT(l,n)=0
      goto 111
 198   continue     
      close(unitinfo)
      fichierdirect(:)=' '
      fichierinfo(:)=' '
!
      if (present(model1Dfile)) then ! on lit aussi le model
         call read_modele(model1Dfile,modele1D)
         modele1D_read=.true.
      endif
!
      files_open=.true.
!----------------------------------------------------------------------------
      end subroutine open_fctp
!----------------------------------------------------------------------------

!-------------------------------------------------------------------- 
      SUBROUTINE open_fctp1(fichierfctp,type)
!
!
!     fichierfctp: character*90, nom du fichier de fonction propre
!                  sans extension
!     unitfctp   : entier, unite logique du fichier de fctp desiree
!
!    irec(0:LLmax,0:NNmax) :tableaux d'entier des numeros de record
!                           des foctions propres
!    omtab(0:LLmax,0:NNmax):tableau de real*8,  w en rad/s
!    qtab(0:LLmax,0:NNmax) :tableau de real*8, attenuations
!    gctab(0:LLmax,0:NNmax):tableau de real*8, vitesses de groupe
!    nocean                :entier, nb de couches de l'ocean
!    nk                    :entier , nombre de couches du model
!    ray(NCM)              :tableau de real*4, couches du model
!--------------------------------------------------------------------
        use global_main_param, only: RA
        use earth_modele
        IMPLICIT NONE
!
      character(len=90), intent(in) :: fichierfctp
      character(len=1), intent(in) :: type
!
      character(len=90) :: fichierdirect,fichierinfo,ext
      real(DP) :: ww,qq,gc,wdiff      
      integer  ::  LEN,unitinfo,indicepointeur,n,l,ibid,i,nocean,ier,k,unitfctp
      logical ::  flag
      real(SP), dimension(:)  , allocatable :: ray_tmp,ray_tmp2
!
!on Lit d'abord le fichier Spheroidaux
!
      unitfctp=118
      unitinfo=73
      if (type=='S') then
         ext='.direct '
         call extension(fichierfctp,fichierdirect,ext)
         ext='.info'
         call extension(fichierfctp,fichierinfo,'.info ')
!
         inquire(file=fichierdirect,EXIST=flag)
         if (.not.flag) then
            print*,'le fichier',fichierdirect,'n''existe pas???'
            stop
         endif
         inquire(file=fichierinfo,EXIST=flag)
         if (.not.flag) then
            print*,'le fichier',fichierinfo,'n''existe pas???'
            stop
         endif
         call get_len(fichierdirect,LEN)
         open(unitfctpS,file=fichierdirect,status='old',access='direct' &
           ,recl=LEN)
         open(unitinfo,file=fichierinfo)
         read(unitfctpS,rec=1)ibid,ibid,ibid,nocean,nkS
         allocate(rayS(nkS),ray_tmp(nkS))
!reading raius and track for discontinuities
         read(unitfctpS,rec=1)ibid,ibid,ibid,nocean,nkS,(ray_tmp(i),i=1,nkS)
!km SP ->m DP
         rayS(:)=dble(ray_tmp(1:nkS))*RA 
         nb_discS=0
         do i=1,nkS-1
            if (abs(rayS(i)-rayS(i+1))<1e-8)  nb_discS=nb_discS+1
         enddo
         allocate(discS_index(0:nb_discS+1),discS_rad(0:nb_discS+1))
         discS_index(0)=0
         discS_index(nb_discS+1)=nkS
         discS_rad(0)  =rayS(1)
         discS_rad(nb_discS+1) =rayS(nkS)
         k=0
         do i=1,nkS-1
            if (abs(rayS(i)-rayS(i+1))<1e-8) then
               k=k+1
               discS_index(k)=i
               discS_rad  (k)=rayS(i)
            endif
         enddo
         deallocate(ray_tmp)
!scaning for LLmax and NNmax:
         LLmaxS=0;NNmaxS=0
10       read(unitinfo,*,END=99) n,l,ww,qq,gc,indicepointeur,wdiff,flag
         if (l>LLmaxS) LLmaxS=l
         if (n>NNmaxS) NNmaxS=n
         goto 10 
 99   continue
      allocate(irecS(0:LLmaxS,0:NNmaxS),omtabS(0:LLmaxS,0:NNmaxS)  &
              ,qtabS(0:LLmaxS,0:NNmaxS),gctabS(0:LLmaxS,0:NNmaxS), stat=ier  )
      irecS(:,:)=0
      if (ier/=0) then
         print*,'LLmaxS=',LLmaxS
         print*,'NNmaxS=',NNmaxS
         print*,(4.+8.*3.)*dble(LLmaxS+1)*dble(NNmaxS+1)/1024./1024.,' Mb'
         stop 'open_fctp: probleme d''allocation memoire'
      endif
      rewind(unitinfo)
!relecture
 11   read(unitinfo,*,END=98) n,l,ww,qq,gc,indicepointeur,wdiff,flag
         irecS (l,n)=indicepointeur
         omtabS(l,n)=ww
         qtabS (l,n)=qq
         gctabS(l,n)=gc
         if ((abs(wdiff).gt.1.d-9.and.flag).or.abs(wdiff).gt.1d-5) irecS(l,n)=0
      goto 11 
 98   continue     
      close(unitinfo)
      fichierdirect(:)=' '
      fichierinfo(:)=' '
      files_openS=.true.
!-----------------------------
      else if (type=='T') then
!-----------------------------

!
!on Lit ensuite le fichier Toroidaux
!
      call extension(fichierfctp,fichierdirect,'.direct ')
      call extension(fichierfctp,fichierinfo,'.info ')
!
      inquire(file=fichierdirect,EXIST=flag)
      if (.not.flag) then
         print*,'le fichier',fichierdirect,'n''existe pas???'
         stop
      endif
      inquire(file=fichierinfo,EXIST=flag)
      if (.not.flag) then
         print*,'le fichier',fichierinfo,'n''existe pas???'
         stop
      endif
      call get_len(fichierdirect,LEN)
      open(unitfctpT,file=fichierdirect,status='old',access='direct' &
                ,recl=LEN)
      open(unitinfo,file=fichierinfo)
      read(unitfctpT,rec=1)ibid,ibid,ibid,nocean,nkT
      allocate(rayT(nkT),ray_tmp(nkT))
!reading raius and track for discontinuities
      read(unitfctpT,rec=1)ibid,ibid,ibid,nocean,nkT,(ray_tmp(i),i=1,nkT)
!km SP ->m DP
      rayT(:)=dble(ray_tmp(1:nkT))*RA
      nb_discT=0
      do i=1,nkT-1
         if (abs(rayT(i)-rayT(i+1))<1e-8)  nb_discT=nb_discT+1
      enddo  
      allocate(discT_index(0:nb_discT+1),discT_rad(0:nb_discT+1))
      discT_index(0)=0
      discT_index(nb_discT+1)=nkT
      discT_rad(0)  =rayT(1)
      discT_rad(nb_discT+1) =rayT(nkT)
      k=0
      do i=1,nkT-1
         if (abs(rayT(i)-rayT(i+1))<1e-8) then
            k=k+1
            discT_index(k)=i
            discT_rad  (k)=rayT(i)
         endif
      enddo  
      deallocate(ray_tmp)
!scaning for LLmax and NNmax:
      LLmaxT=0;NNmaxT=0
 110   read(unitinfo,*,END=199) n,l,ww,qq,gc,indicepointeur,wdiff,flag
         if (l>LLmaxT) LLmaxT=l
         if (n>NNmaxT) NNmaxT=n
         goto 110 
 199   continue
      allocate(irecT(0:LLmaxT,0:NNmaxT),omtabT(0:LLmaxT,0:NNmaxT))
      allocate(qtabT(0:LLmaxT,0:NNmaxT),gctabT(0:LLmaxT,0:NNmaxT))
      irecT(:,:)=0 
      rewind(unitinfo)
!relecture
 111   read(unitinfo,*,END=198) n,l,ww,qq,gc,indicepointeur,wdiff,flag
         irecT (l,n)=indicepointeur
         omtabT(l,n)=ww
         qtabT (l,n)=qq
         gctabT(l,n)=gc
!         if ((abs(wdiff).gt.1.d-9.and.flag).or.abs(wdiff).gt.1d-5) irecT(l,n)=0
      goto 111
 198   continue     
      close(unitinfo)
      fichierdirect(:)=' '
      fichierinfo(:)=' '
      files_openT=.true.
!
      else
         stop 'open_fctp1: type must be S or T'
      endif
      files_open=.true.
!
!----------------------------------------------------------------------------
      end subroutine open_fctp1
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
      subroutine close_fctp
!----------------------------------------------------------------------------
        implicit none
        close(unitfctpS)
        close(unitfctpT)
        deallocate(rayS,discS_index,discS_rad,irecS,omtabS,qtabS,gctabS)
        deallocate(rayT,discT_index,discT_rad,irecT,omtabT,qtabT,gctabT)
        files_open=.false.
!----------------------------------------------------------------------------
      end subroutine close_fctp
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
      subroutine get_Lmax(fmax,LmaxS,LmaxT)
!----------------------------------------------------------------------------
        use def_gparam
        implicit none
        real(DP), intent(in) :: fmax
        integer, intent(out) :: LmaxS,LmaxT
!
        integer :: l
        real(DP) :: wmax
!
        wmax=2._DP*PI*fmax
        LmaxS=-1;LmaxT=-1
!
        if (.not.files_open) stop 'get_lmax: call open_fctp first!'
!on ne cherche que sur le mode fondamentale
        do l=0,LLmaxS
           if (irecS(l,0)/=0.and.omtabS(l,0)<=wmax.and.l>LmaxS) LmaxS=l
        enddo
        do l=1,LLmaxT
           if (irecT(l,0)/=0.and.omtabT(l,0)<=wmax.and.l>LmaxT) LmaxT=l
        enddo
!----------------------------------------------------------------------------
      end subroutine get_Lmax
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
      subroutine get_fmaxc(Lmax,fmax)
!----------------------------------------------------------------------------
        use def_gparam
        implicit none
        integer, intent(in)  :: Lmax
        real(DP), intent(out):: fmax
        integer :: l
!on ne cherche que sur le mode fondamentale spheroidale        
        if (l>LLmaxS) then
           print*,'get_fmaxc: Attention, l> LLmaxS:',l,LLmaxS
           fmax=omtabS(LLmaxS,0)/2._DP/PI
        else
           l=Lmax+1
           do while(irecS(l,0)==0) 
              l=l-1
              if (l<0) stop 'get_fmaxc: Ca va pas ....!'
           enddo
           fmax=omtabS(l,0)/2._DP/PI
        endif        
!----------------------------------------------------------------------------
      end subroutine get_fmaxc
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
      subroutine get_eigenw(eigw,qw,LmaxS,NmaxS,LmaxT,NmaxT)
!----------------------------------------------------------------------------
        use def_gparam
        implicit none
        integer, intent(in) :: LmaxS,LmaxT
        integer, dimension(0:), intent(in) :: NmaxS,NmaxT        
        real(DP), dimension(0:,0:,:), intent(out) :: eigw,qw
!
        integer :: l,n
!
        do l=0,LmaxS
           do n=0,NmaxS(l)
              eigw(n,l,1)=omtabS(l,n)
              if (qtabS(l,n)/=0._DP) then
                 qw(n,l,1)  =omtabS(l,n)/(2._DP*qtabS(l,n))
              else
                 qw(n,l,1)  =0.0_DP
              endif
           enddo
        enddo
        do l=1,LmaxT
           do n=0,NmaxT(l)
              eigw(n,l,2)=omtabT(l,n)
              if (qtabT(l,n)/=0._DP) then
                 qw(n,l,2)  =omtabT(l,n)/(2._DP*qtabT(l,n))
              else
                 qw(n,l,2)  =0.0_DP
              endif
           enddo
        enddo
!----------------------------------------------------------------------------
      end subroutine get_eigenw
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
      subroutine get_Nmax(fmax,NmaxS,NmaxT)
!no
!----------------------------------------------------------------------------
        use def_gparam
        implicit none
        real(DP), intent(in) :: fmax
        integer, dimension(0:), intent(out) :: NmaxS,NmaxT
!
        integer :: l,Lmax,LmaxT,LmaxS,n
        real(DP) :: wmax
!
        wmax=2._DP*PI*fmax
!on recalcul Lmax. Normalement cette operation a deja ete effectue a 
!l'exterieur pour allouer NmaxS et NmaxT
        call get_lmax(fmax,LmaxS,LmaxT)
!
        NmaxS(:)=-1
        NmaxT(:)=-1
        do l=0,LmaxS
           NmaxS(l)=-1
           do n=0,NNmaxS
              if (irecS(l,n)/=0.and.omtabS(l,n)<=wmax.and.n>NmaxS(l)) NmaxS(l)=n
           enddo
        enddo        
        do l=1,LmaxT
           NmaxT(l)=-1
           do n=0,NNmaxT
              if (irecT(l,n)/=0.and.omtabT(l,n)<=wmax.and.n>NmaxT(l)) NmaxT(l)=n
           enddo
        enddo
        NmaxT(0)=0
!----------------------------------------------------------------------------
      end subroutine get_Nmax
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
      SUBROUTINE get_len(file,LEN)
!retourne la longueur des records des fichiers de fctps a acces direct
!produit par classeT ou classeS.
!INPUT:
!    file: character*90, nom du fichier a acces direct
!OUTPUT
!    LEN: entier, longueur du record
!----------------------------------------------------------------------------
      IMPLICIT NONE
!
      character*90 file
      integer LEN
      
!
      open(22,file=file,status='old',access='direct',RECL=4)
      read(22,rec=1)len
      close(22)
!
!-----------------------------------------------------------------
      end SUBROUTINE get_len
!-----------------------------------------------------------------
!--------------------------------------------------------------------------
  real function get_tref()
!--------------------------------------------------------------------------
    implicit none
    get_tref=modele1D%tref
!--------------------------------------------------------------------------
  end function get_tref
!--------------------------------------------------------------------------
!-----------------------------------------------------------------
      subroutine get_1Dmodel(rho,Vs,Vp,Qs,NBV,rad)
!-----------------------------------------------------------------
        use def_gparam
        use module_spline
        implicit none
        integer, intent(in) :: NBV
        real(DP), dimension(NBV), intent(in ) :: rad
        real(DP), dimension(NBV), intent(out) :: rho,Vs,Vp,Qs
!
        integer :: i,k,ka,kb,nd,si,kas,kbs
        real(DP), dimension(nkT,nb_discT+1,4)   :: y2
        real(DP)        :: yp1,ypn
!d'abord on cherche l'index de de depart: (si=start_index)
        do i=1,modele1D%nbcou
           if (abs(modele1D%r(i)-rayT(1))/rayT(1)< 1.E-7_DP) si=i
        enddo
        si=si-1
!        print*,'start_index=',si
!
        if (modele1D_read) then
           do nd=1,nb_discT+1
              ka=discT_index(nd-1)+1 
              kb=discT_index(nd)     
              kbs=kb+si
              kas=ka+si
!
              yp1=(modele1D%rho(kas+1)-modele1D%rho(kas  ))/(rayT(ka+1)-rayT(ka  ))
              ypn=(modele1D%rho(kbs  )-modele1D%rho(kbs-1))/(rayT(kb  )-rayT(kb-1))
              call spline(rayT(ka:kb),modele1D%rho(kas:kbs),yp1,ypn,y2(ka:kb,nd,1))                 
!
              yp1=(modele1D%vpv(kas+1)-modele1D%vpv(kas  ))/(rayT(ka+1)-rayT(ka  ))
              ypn=(modele1D%vpv(kbs  )-modele1D%vpv(kbs-1))/(rayT(kb  )-rayT(kb-1))
              call spline(rayT(ka:kb),modele1D%vpv(kas:kbs),yp1,ypn,y2(ka:kb,nd,2))
!                 
              yp1=(modele1D%vsv(kas+1)-modele1D%vsv(kas  ))/(rayT(ka+1)-rayT(ka  ))
              ypn=(modele1D%vsv(kbs  )-modele1D%vsv(kbs-1))/(rayT(kb  )-rayT(kb-1))
              call spline(rayT(ka:kb),modele1D%vsv(kas:kbs),yp1,ypn,y2(ka:kb,nd,3))                 
!
              yp1=(modele1D%qshear(kas+1)-modele1D%qshear(kas  ))/(rayT(ka+1)-rayT(ka  ))
              ypn=(modele1D%qshear(kbs  )-modele1D%qshear(kbs-1))/(rayT(kb  )-rayT(kb-1))
              call spline(rayT(ka:kb),modele1D%qshear(kas:kbs),yp1,ypn,y2(ka:kb,nd,4))                 
!
           enddo
           do i=1,NBV
              nd=locate(discT_rad,nb_discT,rad(i))
              ka=discT_index(nd-1)+1 
              kb=discT_index(nd)     
              kbs=kb+si
              kas=ka+si              
              rho(i)=splint(rayT(ka:kb),modele1D%rho(kas:kbs)   ,y2(ka:kb,nd,1),rad(i))
              vp (i)=splint(rayT(ka:kb),modele1D%vpv(kas:kbs)   ,y2(ka:kb,nd,2),rad(i))
              vs (i)=splint(rayT(ka:kb),modele1D%vsv(kas:kbs)   ,y2(ka:kb,nd,3),rad(i))
              qs (i)=splint(rayT(ka:kb),modele1D%qshear(kas:kbs),y2(ka:kb,nd,4),rad(i))
           enddo
        else
           stop 'get_1Dmodel: le model n''a pas ete lu!'
        endif        
!-----------------------------------------------------------------
      end subroutine get_1Dmodel
!-----------------------------------------------------------------

!-----------------------------------------------------------------
      subroutine get_1Dmodelb(modele_name,rho,vph,vpv,vsh,vsv,eta,qs,NBV,rad)
!-----------------------------------------------------------------
        use def_gparam
        use module_spline
        implicit none
        integer, intent(in) :: NBV
        character(len=*), intent(in) :: modele_name
        real(DP), dimension(NBV), intent(in ) :: rad
        real(DP), dimension(NBV), intent(out) :: rho,vph,vpv,vsh,vsv,eta,qs
!
        integer :: i,k,ka,kb,nd,si,kas,kbs
        real(DP), dimension(nkT,nb_discT+1,7)   :: y2
        real(DP)        :: yp1,ypn
        type(modele) :: mod1D
!d'abord on cherche l'index de de depart: (si=start_index)
         call read_modele(modele_name,mod1D)
        do i=1,mod1D%nbcou
           if (abs(mod1D%r(i)-rayT(1))/rayT(1)< 1.E-7_DP) si=i
        enddo
        si=si-1
!        print*,'start_index=',si
!
        do nd=1,nb_discT+1
           ka=discT_index(nd-1)+1 
           kb=discT_index(nd)     
           kbs=kb+si
           kas=ka+si
!
           yp1=(mod1D%rho(kas+1)-mod1D%rho(kas  ))/(rayT(ka+1)-rayT(ka  ))
           ypn=(mod1D%rho(kbs  )-mod1D%rho(kbs-1))/(rayT(kb  )-rayT(kb-1))
           call spline(rayT(ka:kb),mod1D%rho(kas:kbs),yp1,ypn,y2(ka:kb,nd,1))                 
!
           yp1=(mod1D%vpv(kas+1)-mod1D%vpv(kas  ))/(rayT(ka+1)-rayT(ka  ))
           ypn=(mod1D%vpv(kbs  )-mod1D%vpv(kbs-1))/(rayT(kb  )-rayT(kb-1))
           call spline(rayT(ka:kb),mod1D%vpv(kas:kbs),yp1,ypn,y2(ka:kb,nd,2))
!                 
           yp1=(mod1D%vsv(kas+1)-mod1D%vsv(kas  ))/(rayT(ka+1)-rayT(ka  ))
           ypn=(mod1D%vsv(kbs  )-mod1D%vsv(kbs-1))/(rayT(kb  )-rayT(kb-1))
           call spline(rayT(ka:kb),mod1D%vsv(kas:kbs),yp1,ypn,y2(ka:kb,nd,3))                 
!
           yp1=(mod1D%qshear(kas+1)-mod1D%qshear(kas  ))/(rayT(ka+1)-rayT(ka  ))
           ypn=(mod1D%qshear(kbs  )-mod1D%qshear(kbs-1))/(rayT(kb  )-rayT(kb-1))
           call spline(rayT(ka:kb),mod1D%qshear(kas:kbs),yp1,ypn,y2(ka:kb,nd,4))                 
!
           yp1=(mod1D%vph(kas+1)-mod1D%vpv(kas  ))/(rayT(ka+1)-rayT(ka  ))
           ypn=(mod1D%vph(kbs  )-mod1D%vpv(kbs-1))/(rayT(kb  )-rayT(kb-1))
           call spline(rayT(ka:kb),mod1D%vph(kas:kbs),yp1,ypn,y2(ka:kb,nd,5))
!                 
           yp1=(mod1D%vsh(kas+1)-mod1D%vsv(kas  ))/(rayT(ka+1)-rayT(ka  ))
           ypn=(mod1D%vsh(kbs  )-mod1D%vsv(kbs-1))/(rayT(kb  )-rayT(kb-1))
           call spline(rayT(ka:kb),mod1D%vsh(kas:kbs),yp1,ypn,y2(ka:kb,nd,6))                 
!
           yp1=(mod1D%eta(kas+1)-mod1D%eta(kas  ))/(rayT(ka+1)-rayT(ka  ))
           ypn=(mod1D%eta(kbs  )-mod1D%eta(kbs-1))/(rayT(kb  )-rayT(kb-1))
           call spline(rayT(ka:kb),mod1D%eta(kas:kbs),yp1,ypn,y2(ka:kb,nd,7))                 
!
        enddo
        do i=1,NBV
           nd=locate(discT_rad,nb_discT,rad(i))
           ka=discT_index(nd-1)+1 
           kb=discT_index(nd)     
           kbs=kb+si
           kas=ka+si              
           rho(i)=splint(rayT(ka:kb),mod1D%rho(kas:kbs)   ,y2(ka:kb,nd,1),rad(i))
           vph (i)=splint(rayT(ka:kb),mod1D%vpv(kas:kbs)   ,y2(ka:kb,nd,2),rad(i))
           vsh (i)=splint(rayT(ka:kb),mod1D%vsv(kas:kbs)   ,y2(ka:kb,nd,3),rad(i))
           qs (i)=splint(rayT(ka:kb),mod1D%qshear(kas:kbs),y2(ka:kb,nd,4),rad(i))
           vpv (i)=splint(rayT(ka:kb),mod1D%vpv(kas:kbs)   ,y2(ka:kb,nd,5),rad(i))
           vsv (i)=splint(rayT(ka:kb),mod1D%vsv(kas:kbs)   ,y2(ka:kb,nd,6),rad(i))
           eta(i)=splint(rayT(ka:kb),mod1D%qshear(kas:kbs),y2(ka:kb,nd,7),rad(i))
        enddo
!-----------------------------------------------------------------
      end subroutine get_1Dmodelb
!-----------------------------------------------------------------
!-----------------------------------------------------------------
      subroutine write_modeT_ascii(unit,s,n,l)
!-----------------------------------------------------------------
        implicit none
        integer, intent(in) :: unit,n,l,s
        real(DP), dimension(2,nkT) :: xdum
        real(DP), parameter :: Rt=6371000.
        integer :: i,iflag
        if (.not.files_openT) stop 'write_modeT_ascii: use open_fctp1 with T type first'
        call readT(n,l,xdum,iflag)
        if (s.eq.3) then
           do i=1,nkT
              write(unit,*)(Rt-rayT(i)),sqrt(float(l*(l+1)))*xdum(1,i)
           enddo
        else if (s.eq.33) then
           do i=1,nkT
              write(unit,*)(Rt-rayT(i)),sqrt(float(l*(l+1)))*xdum(2,i)
           enddo
        endif
!-----------------------------------------------------------------
      end subroutine write_modeT_ascii
!-----------------------------------------------------------------

!-----------------------------------------------------------------
      subroutine write_modeS_ascii(unit,s,n,l)
!-----------------------------------------------------------------
        implicit none
        integer, intent(in) :: unit,n,l,s
        real(DP), dimension(4,nkS) :: xdum
        real(DP), parameter :: Rt=6371000.
        integer :: i,iflag
        if (.not.files_openS) stop 'write_modeS_ascii: use open_fctp1 with S type first'
        call readS(n,l,xdum,iflag)
        if (s.eq.1) then
           do i=1,nkS
              write(unit,*)(Rt-rayS(i)),xdum(1,i)
           enddo
        else if (s.eq.11) then
           do i=1,nkS
              write(unit,*)(Rt-rayS(i)),xdum(2,i)
           enddo
        else if (s.eq.2) then
           do i=1,nkS
              write(unit,*)(Rt-rayS(i)),sqrt(float(l*(l+1)))*xdum(3,i)
           enddo
        else if (s.eq.22) then
           do i=1,nkS
              write(unit,*)(Rt-rayS(i)),sqrt(float(l*(l+1)))*xdum(4,i)
           enddo
        endif
!-----------------------------------------------------------------
      end subroutine write_modeS_ascii
!-----------------------------------------------------------------
!-----------------------------------------------------------------
      subroutine read_interpT(n,l,NBV,rad,fctp)
!-----------------------------------------------------------------
        use def_gparam
        use module_spline
        implicit none
        integer, intent(in) :: n,l,NBV
        real(DP), dimension(NBV), intent(in)  :: rad
        real(DP), dimension(2,NBV),intent(out):: fctp
!
        integer :: iflag,i,k,ka,kb,nd
        real(DP), dimension(2,nkT) :: xdum
        real(DP), dimension(nkT,nb_discT+1)   :: y2
        real(DP)  :: yp1,ypn
!
        call readT(n,l,xdum,iflag)
        if (iflag==0) then
           do k=1,2
              do nd=1,nb_discT+1
                 ka=discT_index(nd-1)+1
                 kb=discT_index(nd)
                 yp1=(xdum(k,ka+1)-xdum(k,ka  ))/(rayT(ka+1)-rayT(ka  ))
                 ypn=(xdum(k,kb  )-xdum(k,kb-1))/(rayT(kb  )-rayT(kb-1))
                 call spline(rayT(ka:kb),xdum(k,ka:kb),yp1,ypn,y2(ka:kb,nd))                 
              enddo 
              do i=1,NBV
                 nd=locate(discT_rad,nb_discT,rad(i))
                 ka=discT_index(nd-1)+1
                 kb=discT_index(nd)                
                 fctp(k,i)=splint(rayT(ka:kb),xdum(k,ka:kb),y2(ka:kb,nd),rad(i))
              enddo
           enddo
        else
           fctp(:,:)=0.0_DP
        endif
!-----------------------------------------------------------------
      end subroutine read_interpT
!-----------------------------------------------------------------

!-----------------------------------------------------------------
      subroutine read_interpS(n,l,NBV,rad,fctp)
!-----------------------------------------------------------------
        use def_gparam
        use module_spline
        implicit none
        integer, intent(in) :: n,l,NBV
        real(DP), dimension(NBV)  , intent(in):: rad
        real(DP), dimension(4,NBV),intent(out):: fctp
!
        integer :: iflag,i,k,ka,kb,nd
        real(DP), dimension(4,nkS) :: xdum
        real(DP), dimension(nkS,nb_discS+1) :: y2
        real(DP) :: yp1,ypn
!
        call readS(n,l,xdum,iflag)
        if (iflag==0) then
           do k=1,4
              do nd=1,nb_discS+1
                 ka=discS_index(nd-1)+1
                 kb=discS_index(nd)
                 yp1=(xdum(k,ka+1)-xdum(k,ka  ))/(rayS(ka+1)-rayS(ka  ))
                 ypn=(xdum(k,kb  )-xdum(k,kb-1))/(rayS(kb  )-rayS(kb-1))
                 call spline(rayS(ka:kb),xdum(k,ka:kb),yp1,ypn,y2(ka:kb,nd))
              enddo
              do i=1,NBV
                 nd=locate(discS_rad,nb_discS,rad(i))
                 ka=discS_index(nd-1)+1
                 kb=discS_index(nd)                
                 fctp(k,i)=splint(rayS(ka:kb),xdum(k,ka:kb),y2(ka:kb,nd),rad(i))
              enddo
           enddo
        else
           fctp(:,:)=0.0_DP
        endif
!-----------------------------------------------------------------
      end subroutine read_interpS
!-----------------------------------------------------------------

!-----------------------------------------------------------------
      integer function locate(inter,n,x)
!-----------------------------------------------------------------
         use def_gparam
         implicit none      
         integer, intent(in) :: n
         real(DP), dimension(0:n+1), intent(in) :: inter
         real(DP), intent(in) :: x         
!
         integer :: i
         i=0
         if (abs (x-inter(0))/max(x,inter(0)) <1.E-10_DP) then
            i=1
         else if ( abs (x-inter(n+1))/max(x,inter(n+1))  <1.E-10_DP) then
            i=n
         else
            do while (  x < inter(i) .or. x > inter(i+1) )
               i=i+1
               if (i > n) then
                  print*,'i=',i,'n=',n
                  print*,'x=',x
                  print*,'inter=',inter
                  stop 'locate: failed to locate x in inter!'
               endif
            enddo
         endif
         locate=i+1
!-----------------------------------------------------------------
      end function locate
!-----------------------------------------------------------------

!-----------------------------------------------------------------
      subroutine readT(n,l,xdum,iflag)
!lit la fonction propre (n,l) apres ouverture du fichier a acces direct
!par open_fctp
!    
!----------------------------------------------------------------
      use def_gparam
      implicit none
      integer, intent(in) :: l,n
      integer, intent(out):: iflag
      real(DP), dimension(2,nkT), intent(out) ::  xdum
!
      real(DP) :: WW,QQ,GC
      real(SP), dimension(nkT*6) :: tt
      integer ::  i,nn,ll,nvec,nb
!    
!nve!est la longueur necessaire pour atteindre la profondeur
!desire pour w wt w':
      if (.not.files_open) stop 'readT: call open_fctp first!'
      nb=nkT
      nvec=int(nb*2)
      if (l<=LLmaxT.and.n<=NNmaxT) then
         if (irecT(l,n)/=0) then
            read(unitfctpT,rec=irecT(l,n)) nn,ll,WW,QQ,GC,(tt(I),I=1,nvec)
            do i=1,nkT
               xdum(1,i)=dble(tt(2*(nkT-i)+1))
               xdum(2,i)=dble(tt(2*(nkT-i)+2))
            enddo
            iflag=0
         else
            iflag=2
            xdum(:,:)=0.0_DP
         endif
      else
         iflag=1
      endif
!-----------------------------------------------------------------
      end subroutine readT
!-----------------------------------------------------------------

!-----------------------------------------------------------------
      subroutine readS(n,l,xdum,iflag)
!-----------------------------------------------------------------
!lit la fonction propre (n,l) apres ouverture du fichier a acces direct
!par open_fctp
!INPUT:
!    ifd: entier, inite logique du fichier a acces direct
!    irec(0:LLmax,0:NNmax),tableau des numeros de records sortie
!                          de open_fctp
!    n et l, entier
!    nk entier:nombre de couches du model
!    np entier:profondeur jusqu'a laquelle on veut lire ls fctps 
!(de np a nk=surface)
!OUPUT:
!    xdum1..6(NCM),fonction propres et derivees
!Attention les fctps 5 et6 ne sont pas lu, je laisse xdum5 et 6 pour 
!des raisons  de conpatibilites ave!d'autres prog.
!    iflag: 1 pas de pb,
!           0 fin du fichier ( ne devrait jamais arriver)
!          -1 erreur a la lecture
!          -2 l ou n depasse LLmax ou NNmax 
!----------------------------------------------------------------
      use def_gparam
      implicit none
      integer, intent(in) :: l,n
      integer, intent(out):: iflag
      real(DP), dimension(4,nkS), intent(out) ::  xdum
!
      real(DP) :: WW,QQ,GC
      real(SP), dimension(nkS*6) :: tt
      integer :: nvec,i,nn,ll,nb
!
      if (.not.files_open) stop 'readS: call open_fctp first!'
      nb=nkS
!nve!est la longueur necessaire pour atteindre la profondeur
!desire pour v, v' u, u' ( on ne lit pas p et p', si on en a besoin il
!ecrire une autre routine de lecture.)      
      nvec=int(nb*4) 
      if (l<=LLmaxS.and.n<=NNmaxS) then 
         if (irecS(l,n)/=0) then
!reading and killing discontinuities
            read(unitfctpS,rec=irecS(l,n)) nn,ll,WW,QQ,GC,(tt(I),I=1,nvec)
            do i=1,nkS
               xdum(1,i)=dble(tt(4*(nkS-i)+1))
               xdum(2,i)=dble(tt(4*(nkS-i)+2))
               xdum(3,i)=dble(tt(4*(nkS-i)+3))
               xdum(4,i)=dble(tt(4*(nkS-i)+4))
            enddo
            iflag=0
         else
            iflag=2
            xdum(:,:)=0.0_DP
         endif
      else
         iflag=1
      endif
!---------------------------------------------------------------------
      end subroutine readS
!---------------------------------------------------------------------
!---------------------------------------------------------------------
      SUBROUTINE  extension(fichier1,fichier2,ext) 
!---------------------------------------------------------------------
      IMPLICIT NONE
!
      character(len=*), intent(in)  :: fichier1,ext 
      character(len=*), intent(out) :: fichier2
!
      integer :: lenfichier1,lenext 
!
      lenfichier1=INDEX (fichier1,' ') -1   
      lenext=INDEX (ext,' ') -1   
      if ((lenfichier1+lenext).gt.80) then 
         print*,'il y a trop de caracter a  ',fichier1 
      endif 
!
      fichier2(1:lenfichier1)=fichier1(1:lenfichier1) 
      fichier2(lenfichier1+1:(lenfichier1+lenext))=ext(1:lenext) 
      fichier2(lenfichier1+lenext+1:)=' '
!---------------------------------------------------------------------
      end SUBROUTINE  extension
!---------------------------------------------------------------------
!---------------------------------------------------------------------
      logical function mode_existe(n,l,q)
!---------------------------------------------------------------------
        implicit none
        integer :: n,l,q
!
        if (q==1) then
           if (irecS(l,n)/=0) then
              mode_existe=.true.
           else
              mode_existe=.false.
           endif
        else
           if (irecT(l,n)/=0) then
              mode_existe=.true.
           else
              mode_existe=.false.
           endif
        endif
!        
!---------------------------------------------------------------------
      end function mode_existe
!---------------------------------------------------------------------
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
end module util_fctp
