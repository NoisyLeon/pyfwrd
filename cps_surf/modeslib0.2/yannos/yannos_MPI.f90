!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
program yannos_MPI
!version a moi de minos avec MPI
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  use minos
  use neominos
  use bits, only: l,eps,cond_limite,wgrav,nord,epsint,epsint_base
  use rindex, only: n
  use def_gparam
  use earth_modele
  use yannos_flag
  use layer
  implicit none
  include 'mpif.h'
  character(len=100) :: modele_name,fichier_per,prefix,fichierdirect &
                      ,fichierinfo,fichierlost
  integer  :: type_fctp,lmin,lmax,nvec,ieigen,len,iinfo,iper,nbran,lmax_p &
              ,dl,ldeb,i,nmin,ilost,ier,rang,nbproc,j,indicerec,unit
  real(DP) :: epsin,wgravin,fmin,fmax,fdeb,ra,bid
  type(modele):: modele_de_terre
  logical :: flag,flag_deb,first,termine
!
  real(DP), dimension(:,:), allocatable :: f
  integer , dimension(:)  , allocatable :: n_zero
!
  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nbproc,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rang  ,ier)  
  if (rang==0)  print*,'Yannos_MPI est une version de minos faites par Y.Capdeville (140999)'
!lecture des entrees:
  unit=117
  open(unit,file='yannos.dat',status='old')
!lecture des entrees:
!  print*,'Entrer le nom du modele:'
  read(unit,*)
  read(unit,100) modele_name
!  print*,'Entrer le nom du fichier de sortie:'
  read(unit,*)
  read(unit,100) fichier_per
!  print*,'Entrer le prefix des fichiers de fonctions propres'
  read(unit,*)
  read(unit,100) prefix
!  print*,'entrer le type de fctp (2 pour t 3 pour s):'
  read(unit,*)
  read(unit,*) type_fctp
  if (type_fctp/=2.and.type_fctp/=3) stop '2 OU 3 qu''on t''dit!'
!  print*,'enter eps,eps integration and wgrav'
  read(unit,*)
  read(unit,*) epsin,epsint_base,wgravin
!  print*,'enter lmin,lmax,fmin*1000,fmax*1000,nbran'
  read(unit,*)
  read(unit,*) lmin,lmax,fmin,fmax,nbran
  read(unit,*) 
  call read_layers(unit)
  read(unit,*) 
  call read_flags(unit)
  if (wgravin>0.d0 .and. cancel_gravity.and.rang==0) then
     print*,' cancel_gravity requested, I''m setting wgrav to 0!'
     wgravin=0.
  endif
  close(unit)
  if (fmin>fmax) stop 'fmin>fmax !'
  if (lmin>lmax) stop 'lmin>lmax !'
  if (nbran<=0) stop 'Nombre de branche demandee <= 0 !'
!
!bord libre
  cond_limite=1
!
  fmin=fmin/1000._DP
  fmax=fmax/1000._DP
!
! lecture du modele de terre
!
  call read_modele(modele_name,modele_de_terre)
  call model(modele_de_terre) 
  call init_layer(modele_de_terre)
  call MPI_BARRIER(MPI_COMM_WORLD,ier)  
  if (rang==0) then
     print*,'++++++++++++++++++++++++LAYER++++++++++++++++++++++++'
     do i=1,nb_lay
        print*,'Couche:',i
        print*,'Indice debut:',i_la1(i)
        print*,'Indice fin  :',i_la2(i)        
     enddo
     print*,'+++++++++++++++++++++++++++++++++++++++++++++++++++++'
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ier)
!
!ouverture des fichiers sorties
!
  ieigen=31
  iinfo =32
  iper  =33
  ilost =34
  call extension(fichier_per,fichierlost,'.lost ')
  call extension(prefix,fichierdirect,'.direct ')
  call extension(prefix,fichierinfo,'.info ')
  if (type_fctp.eq.3) then
!     nvec=6*n
     nvec=6*nbcou_lay
  else if (type_fctp.eq.2) then
!     nvec=2*n
     nvec=2*nbcou_lay
  endif
  len=8+3*8+nvec*4
  if (rang==0) then
     open(iper,file=fichier_per,status='replace')
     open(ieigen,file=fichierdirect,access='direct',recl=len,status='replace')
     ra=modele_de_terre%r(modele_de_terre%nbcou)
!     write(ieigen,REC=1) len,-1,-1,modele_de_terre%nbceau,modele_de_terre%nbcou &
!          ,(real(modele_de_terre%r(i)/ra,SP),i=1,modele_de_terre%nbcou)        
     write(ieigen,REC=1) len,-1,-1,modele_de_terre%nbceau,nbcou_lay &
          ,((real(modele_de_terre%r(i)/ra,SP),i=i_la1(j),i_la2(j)),j=1,nb_lay)
     open(iinfo,file=fichierinfo,status='replace')
     open(ilost,file=fichierlost,status='replace')
     WRITE(iper,111) real(EPSint_base),real(EPS),real(WGRAV)
111  FORMAT(/,'INTEGRATION PRECISION =',G12.4,'  ROOT PRECISION =',  &
          G12.4,'  GRAVITY CUT OFF =',G12.4,' RAD/S',///,6x,'MODE',      &
          8x,'W(RAD/S)',7x,'W(MHZ)',10x,'T(SECS)',6x,'GRP VEL(KM/S)',    &
          8x,'Q',13x,'RAYLQUO',/)
    write(iper,*)'==============================================================================='
    write(iper,*)'Flag selectionne pour ce run (yannos_flag.f90):'
    if (force_fmin) then
       write(iper,*)'La recherche commence toujours a fmin (pas standard)'
    else
       write(iper,*)'La recherche commence a un fmin flottant (standard)'
    endif
    if (cancel_gravity) then
       write(iper,*)'La gravite a ete desactivee totalement (pas standard)'
    else
       write(iper,*)'Les termes de gravite (pas de redistribution) sont actives (standard) '
    endif
    if (never_use_startlevel) then
       write(iper,*)'Startlevel n''est jamais utilise (pas standard, mais conseille pour les test homogenes)'
    else
       write(iper,*)'Startlevel est utilise a partir de l=',l_startlevel,' (standard)'
    endif
    if (use_tref) then
       write(iper,*)'Utilisation de la periode de reference pour calculer les params elastics'
    else
       write(iper,*)'La periode de reference pour calculer les params elastics est descativee'
    endif
    if (check_modes) then
       write(iper,*)'check_modes active. (Pas standart: vire les modes de la graine)'
    else
       write(iper,*)'check_modes desactive. (standard)'
    endif
    if (use_remedy) then
       write(iper,*)'Utilisation de remdy (standard)'
    else
       write(iper,*)'Desactivation de remdy (standard)'
    endif
    if (rescue) then
       write(iper,*)'Rescue active: recherche systematique des modes perdus'
    else
       write(iper,*)'Rescue desactive'
    endif
    if (keep_bad_modes) then
       print*,'On garde les modes mal calcules (pas une bonne idee ca!)'
    else 
       print*,'On jete les modes mal calcule!'
    endif
    if (force_systemic_search) then
       write(iper,*)' ATTENTION, La recherche systematique a ete activee pour les toroidaux!!!!!!'
       print*,' ATTENTION, La recherche systematique a ete activee pour les toroidaux!!!!!!'
    endif
    write(iper,*)'Format de sortie: ',modout_format
    write(iper,*)'==============================================================================='
    close(iper);close(ieigen);close(iinfo);close(ilost)
  endif  
!
!c'est parti...
!
  if (modout_format/='ipg') stop 'Le format de sortie doit etre ''ipg'' &
                                  & pour la version MPI'
  eps=epsin
  epsint=epsint_base
  wgrav=wgravin
!dans minos:
  call steps(epsint_base)
  if (type_fctp==2) then
     lmin=min(lmin,1)     
  endif
  ldeb=lmin+rang
  dl  =nbproc
  do while(mod(Lmax-lmin+1,dl)/=0)
     Lmax=Lmax+1
  enddo
  termine=.false.
  allocate(f(nbran,lmin:lmax))
  allocate(n_zero(lmin:lmax))
  f(:,:)=0.0_DP
  n_zero(:)=0
  flag_deb=.true.
  lmin=ldeb
  indicerec=1
!=============================
  do l=lmin,lmax,dl
!=============================
     if (never_use_startlevel) then
        use_startlevel=.false.
     else
        if (l>l_startlevel) then
           use_startlevel=.true.
        else
           use_startlevel=.false.
        endif
     endif
     if (l.gt.ldeb) then
        if (flag_deb.or.n_zero(l-dl).ne.0     &
             .or.n_zero(l-dl).ne.0) then
           flag=.true.
           lmax_p=l
        else
           flag=.false.
        endif
     else
        flag=.true.
     endif
     call MPI_BARRIER(MPI_COMM_WORLD,ier)
!"""""""""""""""""""""""""""""""""
     if (.not.termine) then
!"""""""""""""""""""""""""""""""""
     if (flag) then
        if ( l > max(ldeb+dl,1)) then
           if (f(1,l-dl) > 1.0E-12_DP.and.nmin==0.and..not.force_fmin) &
                     fdeb=max(fmin,f(1,l-dl)-2.E-5_DP)
           first  =.false.
        else
           fdeb=fmin 
           first  =.true.
        endif
        print*,' fdeb=',fdeb,' l=',l
     endif
     if (type_fctp==3) then
        call get_fp_modele_R(l,nbran,fdeb,fmax,f(:,l),n_zero(l),first,nmin)
     else
        call get_fp_modele_L(l,nbran,fdeb,fmax,f(:,l),n_zero(l),first,nmin)
     endif
!"""""""""""""""""""""""""""""""""
     else
!"""""""""""""""""""""""""""""""""
        n_zero(l)=0
!"""""""""""""""""""""""""""""""""
     endif
!"""""""""""""""""""""""""""""""""
     do j=0,nbproc-1 
        call MPI_BARRIER(MPI_COMM_WORLD,ier)    
        if (j==rang) then
           if (n_zero(l)>0) then
              open(iper,file=fichier_per,status='old',position='append')
              open(ieigen,file=fichierdirect,status='old',              &
                                             access='direct',recl=len)
!on positionne sur le dernier rec: (bidouille pour que nexrect donne bien indicerec+1)
              read (ieigen,rec=indicerec+1,iostat=ier) bid
!
              open(iinfo,file=fichierinfo,status='old',position='append')
              open(ilost,file=fichierlost,status='old',position='append')           
              print*,'ecriture des fctps pour l=',l,' par le rang ',rang
              do i=1,n_zero(l)
                 nord=i-1+nmin
                 if ((f(i,l)+2._DP)>1.E-15_DP.and.f(i,l)>1.E-15_DP) then
!                    print*,'pour n=',nord
                    if (type_fctp==3) then
                       call write_fctpR(l,f(i,l),iper,iinfo,ieigen,ilost)
                    else
                       call write_fctpL(l,f(i,l),iper,iinfo,ieigen)
                    endif
                 else
                    if (type_fctp==3) then
                       write(ilost,101) nord,'S',l
                    else
                       write(ilost,101) nord,'T',l
                    endif
                    print*,'Pour n=',nord,', la freq p n''a pu etre determinee'
                 endif
              enddo
!on cherche le numero du dernier record:
              inquire(ieigen,nextrec=indicerec)
              indicerec=indicerec-1
!
              close(iper);close(ieigen);close(iinfo);close(ilost)
           else if (.not.termine) then
              if (l>ldeb) then
                 termine=.true.
                 write(*,257) rang
              endif  
           endif
        endif
!on donne indicerec  a tout le monde
        call MPI_BCAST(indicerec,1,MPI_INTEGER,max(j,0) &
             ,MPI_COMM_WORLD,ier)
     enddo
!
257 format('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!',/, &
           'Le rang ',i3,' a terminé son travail.',/,   &
           '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' )
!
!=============================
  enddo
!=============================
  deallocate(f,n_zero)
  call MPI_FINALIZE(ier)
!
100 format(a100)
101 format(i5,1x,a1,1x,i4,' indeterminee, hors precision')
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
end program yannos_MPI
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
