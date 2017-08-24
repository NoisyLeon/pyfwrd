!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
program yannos
!This program comes form minos (form G. Master, J.Woodhouse and berfore 
! them Gilbert)
! I have only rewritten some parts in f90 in order to be able
! to interstand what's going on and to be able to modifie it.
! I haven't change anything to the principle of minos.
! I have add some little things that are suitable for my use, hopping
! I haven't add to much mistakes.
! Yann Capdeville 22/04/99
!version a moi de minos
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  use minos
  use neominos
  use bits, only: l,eps,epsint,epsint_base,cond_limite,wgrav,nord
  use rindex, only: n
  use def_gparam
  use earth_modele
  use yannos_flag
  use layer
  implicit none
  character(len=100) :: modele_name,fichier_per,prefix_keep,fichierdirect &
                      ,fichierinfo,fichierlost,fichier_per_keep,prefix
  integer  :: type_fctp,type_fctp_in,lmin,lmin_in,lmax,nvec,ieigen,len,iinfo,iper,nbran,lmax_p &
              ,dl,ldeb,i,nmin,ilost,j,lf,ll,nn,unit,tdeb,tfin
  real(DP) :: epsin,wgravin,fmin,fmax,fdeb,ra,ww
  type(modele):: modele_de_terre
  logical :: flag,flag_deb,first,termine
!
  real(DP), dimension(:,:), allocatable :: f
  integer , dimension(:)  , allocatable :: n_zero  
   modele_name(:)=' '	
   fichier_per(:)=' ' 
   prefix     (:)=' '
   fichierdirect(:)=' '
   fichier_per_keep(:)=' ' 
   prefix_keep(:)=' '
   fichierinfo(:)=' '
   fichierlost(:)=' '	
  print*,'Yannos est une version de minos faites par Y.Capdeville (a partir de 140999)'
  unit=117
  open(unit,file='yannos.dat',status='old')
!lecture des entrees:
!  print*,'Entrer le nom du modele:'
  read(unit,*)
  read(unit,100) modele_name
!  print*,'Entrer le nom du fichier de sortie:'
  read(unit,*)
  read(unit,100) fichier_per_keep
!  print*,'Entrer le prefix des fichiers de fonctions propres'
  read(unit,*)
  read(unit,100) prefix_keep
!  print*,'entrer le type de fctp (2 pour t 3 pour s):'
  read(unit,*)
  read(unit,*) type_fctp_in
  if (type_fctp_in/=2.and.type_fctp_in/=3.and.type_fctp_in/=0) stop '2 OU 3 qu''on t''dit!'
!  print*,'enter eps,eps integration and wgrav'
  read(unit,*)
  read(unit,*) epsin,epsint_base,wgravin
!  print*,'enter lmin,lmax,fmin*1000,fmax*1000,nbran'
  read(unit,*)
  read(unit,*) lmin_in,lmax,fmin,fmax,nbran
  read(unit,*) 
  call read_layers(unit)
  read(unit,*) 
  call read_flags(unit)
  if (wgravin>0.d0 .and. cancel_gravity) then
     print*,' cancel_gravity requested, I''m setting wgrav to 0!'
     wgravin=0.
  endif
  close(unit)
  if (fmin>fmax) stop 'fmin>fmax !'
  if (lmin_in>lmax) stop 'lmin>lmax !'
  if (nbran<=0) stop 'Nombre de branche demandee <= 0 !'
!===============================================================================
  print*,'==============================================================================='
  print*,'Flag selectionne pour ce run (yannos_flag.f90):'
  if (force_fmin) then
     print*,'La recherche commence toujours a fmin (pas standard)'
  else
     print*,'La recherche commence a un fmin flottant (standard)'
  endif
  if (cancel_gravity) then
     print*,'La gravite a ete desactivee totalement (pas standard)'
  else
     print*,'Les termes de gravite (pas de redistribution) sont actives (standard) '
  endif
  if (never_use_startlevel) then
     print*,'Startlevel n''est jamais utilise (pas standard, mais conseille pour les test homogenes)'
  else
     print*,'Startlevel est utilise a partir de l=',l_startlevel,' (standard)'
  endif
  if (use_tref) then
     print*,'Utilisation de la periode de reference pour calculer les params elastics'
  else
     print*,'La periode de reference pour calculer les params elastics est descativee'
  endif
  if (check_modes) then
     print*,'check_modes active. (Pas standart: vire les modes de la graine)'
  else
     print*,'check_modes desactive. (standard)'
  endif
  if (use_remedy) then
     print*,'Utilisation de remdy (standard)'
  else
     print*,'Desactivation de remdy (standard)'
  endif
  if (rescue) then
     print*,'Rescue active: recherche systematique des modes perdus'
  else
     print*,'Rescue desactive'
  endif
  if (force_systemic_search) print*,' ATTENTION, La recherche systematique a ete activee pour les toroidaux!!!!!!'
  print*,'Format de sortie: ',modout_format
  print*,'==============================================================================='
!===============================================================================

!bord libre ( cond_limite=1)
  cond_limite=1
!  cond_limite=2
!
  fmin=fmin/1000._DP
  fmax=fmax/1000._DP
  if (type_fctp_in/=0) then
     tdeb=type_fctp_in; tfin=type_fctp_in
  else
     tdeb=2; tfin=3
  endif
  do type_fctp=tdeb, tfin
     allocate(f(nbran,lmin_in:lmax))
     allocate(n_zero(lmin_in:lmax))
     f(:,:)=0.0_DP
     n_zero(:)=0
     print*,'*********************************************************************'
     if (type_fctp==2) then
        print*,'             Working on Toroidal modes            '
     else
        print*,'             Working on Spheroidal modes            '
     endif
     print*,'*********************************************************************'
!
! lecture du modele de terre
!
     call read_modele(modele_name,modele_de_terre)
     call model(modele_de_terre)
     call init_layer(modele_de_terre)
     print*,'++++++++++++++++++++++++LAYER++++++++++++++++++++++++'
     do i=1,nb_lay
        print*,'Couche:',i
        print*,'Indice debut:',i_la1(i)
        print*,'Indice fin  :',i_la2(i)        
     enddo
!
!ouverture des fichiers sorties
!
     ieigen=31
     iinfo =32
     iper  =33
     ilost =34

     fichier_per=fichier_per_keep
     prefix     =prefix_keep
     if (type_fctp_in/=0) then
        call extension(fichier_per,fichierlost,'.lost ')
        call extension(prefix,fichierdirect,'.direct ')
        call extension(prefix,fichierinfo,'.info ')
     else
        if (type_fctp==3) then
           call extension(fichier_per_keep,fichier_per,'S ')
           call extension(fichier_per,fichierlost,'.lost ')
           call extension(prefix,fichierdirect,'S.direct ')
           call extension(prefix,fichierinfo,'S.info ')        
        else   if (type_fctp==2) then
           call extension(fichier_per_keep,fichier_per,'T ')
           call extension(fichier_per,fichierlost,'.lost ')
           call extension(prefix,fichierdirect,'T.direct ')
           call extension(prefix,fichierinfo,'T.info ')        
        endif
     endif

     open(iper,file=fichier_per)
!
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
     write(iper,*)'Format de sortie: ',modout_format
     if (force_systemic_search) then
        write(iper,*)' ATTENTION, La recherche systematique a ete activee pour les toroidaux!!!!!!'
     endif
     write(iper,*)'==============================================================================='
!
     if (type_fctp.eq.3) then
        nvec=6*nbcou_lay
     else if (type_fctp.eq.2) then
        nvec=2*nbcou_lay
     endif
     if (modout_format=='ucb'.or.modout_format=='olm') then
!sortie au format Berkeley ou "old minos"
        open(ieigen,file=prefix,form='unformatted')
     else
!sortie au format yann (ipg)
        len=8+3*8+nvec*4
        open(ieigen,file=fichierdirect,access='direct',recl=len)
     endif
!
     ra=modele_de_terre%r(modele_de_terre%nbcou)
     if (modout_format/='ucb'.and.modout_format/='olm') &
          write(ieigen,REC=1) len,-1,-1,modele_de_terre%nbceau,nbcou_lay &
          ,((real(modele_de_terre%r(i)/ra,SP),i=i_la1(j),i_la2(j)),j=1,nb_lay)
     open(iinfo,file=fichierinfo)
     open(ilost,file=fichierlost)
!pour ressembler a minos:
     WRITE(iper,111) real(EPSINT_base),real(EPS),real(WGRAVIN)
111  FORMAT(/,'INTEGRATION PRECISION =',G12.4,'  ROOT PRECISION =',  &
          G12.4,'  GRAVITY CUT OFF =',G12.4,' RAD/S',///,6x,'MODE',      &
          8x,'W(RAD/S)',7x,'W(MHZ)',10x,'T(SECS)',6x,'GRP VEL(KM/S)',    &
          8x,'Q',13x,'RAYLQUO',/)
!
!c'est parti...
!
     termine=.false.
     eps=epsin
     epsint=epsint_base
     wgrav=wgravin
!dans minos:
     call steps(epsint_base)
     if (type_fctp==2) then
        lmin=max(lmin_in,1)
     else
        lmin=lmin_in
     endif
     dl=1
     flag_deb=.true.
     ldeb=lmin
!=============================
     if (restart) then
        n_zero(:)=-1
        f(:,:)   =-2.d0
10      read(iinfo,*,end=75) nn,ll,ww
        n_zero(ll)=max(n_zero(ll),nn+1)
        f(nn+1,ll)=ww/2._DP/PI
        goto 10
75      continue     
     endif
     rewind iinfo
!=============================
     fdeb=fmin
!=============================
     do l=lmin,lmax
!  do l=140,140
!=============================
        if (l>l_startlevel) then
           use_startlevel=.true.
        else
           use_startlevel=.false.
        endif
        if (never_use_startlevel) use_startlevel=.false.
        
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
!"""""""""""""""""""""""""""""""""
        if (.not.termine) then
!"""""""""""""""""""""""""""""""""
           if (flag) then
              if (type_fctp==2) then
                 lf=10
              else
                 lf=1
              endif
              if ( l > max(ldeb,lf)) then
                 
                 if (f(1,l-dl) > 1.0E-12_DP.and.nmin==0.and..not.force_fmin) &
                      fdeb=max(fmin,f(1,l-dl)-4.E-5_DP)
                 first  =.false.
              else
                 fdeb=fmin 
                 first  =.true.
              endif
              print*,' fdeb=',fdeb,' l=',l
           endif
           if (.not.(restart.and.n_zero(l)/=-1)) then
              if (type_fctp==3) then
                 call get_fp_modele_R(l,nbran,fdeb,fmax,f(:,l),n_zero(l),first,nmin)
              else
                 call get_fp_modele_L(l,nbran,fdeb,fmax,f(:,l),n_zero(l),first,nmin)
              endif
           endif
!"""""""""""""""""""""""""""""""""
        else
!"""""""""""""""""""""""""""""""""
           n_zero(l)=0
!"""""""""""""""""""""""""""""""""
        endif
!"""""""""""""""""""""""""""""""""
        if (n_zero(l)>0) then
!"""""""""""""""""""""""""""""""""
!test
!print*,'f(1,30)=',(f(1,30)-4.933379116006998E-003)/f(1,30)
!f(1,30)=4.933379116006998E-003
           print*,'ecriture des fctps pour l=',l,'...'
           do i=1,n_zero(l)
              nord=i-1+nmin
              if ((f(i,l)+2._DP)>1.E-15_DP) then
                 print '(a7,i4)','pour n=',nord
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
!"""""""""""""""""""""""""""""""""
        else if (.not.termine) then
!"""""""""""""""""""""""""""""""""
           if (l>ldeb) then
              termine=.true.
              print*,'J''ai fini ... '
           endif
!"""""""""""""""""""""""""""""""""
        endif
!"""""""""""""""""""""""""""""""""
        !     call flush(iper);call flush(iinfo); call flush(ieigen); call flush(ilost)
!=============================
     enddo
!=============================
     deallocate(f,n_zero)
     close(iper);close(ieigen);close(iinfo);close(ilost)
!-------------------------------
     call clean_minos_memory
  enddo
!-------------------------------
!
100 format(a100)
101 format(i3,1x,a1,1x,i4,' indeterminee, hors precision')
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
end program yannos
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
