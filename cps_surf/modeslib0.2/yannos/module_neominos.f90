!--------------------------------------------------------------------------
module neominos
! module utilisant certaines routine de minos pour calculer l'operatuer 
! de couplage
!--------------------------------------------------------------------------
use def_gparam
private
public :: get_fp_modele_R,get_fp_modele_L,get_ANNp_R,get_ANNp_L &
         ,ANN_for_splineR,ANN_for_splineL,sum_apsom,ANNcR,ANNcL
!
contains
!#########################################################################

!------------------------------------------------------------------
  subroutine get_fp_modele_R(l_in,nmax,fmin,fmax,f_zero,n_zero,first,nmin)
! recherche de Frequences  Propres.
!------------------------------------------------------------------
    use bits, only: jcom,l,epsint,epsint_base
    use def_gparam
    use yannos_flag
    implicit none
!
    integer , intent(in) :: l_in,nmax
    integer , intent(out) :: nmin
    real(DP), intent(in) :: fmin,fmax
    logical, intent(in) :: first
    integer , intent(out):: n_zero
    real(DP), dimension(:), intent(out):: f_zero
!
    real(DP), dimension(:), pointer :: table_freq,table_det
    real(DP):: wmin,wmax,wp,tmp,bid,fc
    integer :: nbmodes_min,i,j,nbi
    logical :: erreur,flag
!
    logical,  dimension(:), allocatable :: rescue_flag
    real(DP), dimension(:), allocatable :: f_zero_tmp
!
!on stipule au module minos qu'on veux les modes de rayleigh
!
    l=l_in
    if (l==0) then
       jcom=1 !modes radiaux
    else
       jcom=3
    endif
!
    wmin=2.d0*PI*fmin
    wmax=2.d0*PI*fmax 
!on etablit la precision d'integration a la presision de base   
    epsint=epsint_base
!
!
! on etablit la table de recherche en freqeunce
  call get_table_freq(l,nmax,wmin,wmax,table_freq,table_det &
                     ,nbmodes_min,n_zero) 


  nbi=0
!do i=nbmodes_min,nbmodes_min+n_zero
!print*,i,table_freq(i)
!enddo
!stop
  do while (any(table_freq(:)<=-2._DP).and.l/=1.and.nbi<5) 
     print*,'Je suis oblige d''incrementer la precision d''integration pour R l=',l
     epsint=epsint*1.d-2
     call get_table_freq(l,nmax,wmin,wmax,table_freq,table_det &
          ,nbmodes_min,n_zero) 
     nbi=nbi+1
  enddo
  epsint=epsint_base
!
  allocate(rescue_flag(1:n_zero),f_zero_tmp(1:n_zero+100))
  f_zero_tmp(:)=-10._DP
  f_zero    (:)=-10._DP
  rescue_flag(:)=.false.
!
  nmin=max(nbmodes_min,0)
  n_zero=n_zero-nmin
  flag=.false.
  if (n_zero>UBOUND(f_zero,dim=1)) then
     print*,'n_zero=',n_zero,UBOUND(f_zero,dim=1)
     stop 'neominos: probleme de taille: increase Nmax'   
  endif
! recherche des frequences:
  do i=1,n_zero
!print*,i,' over ',n_zero
     j=i+nmin
     if (table_freq(j-1)+2._DP<=1.E-15_DP.or.table_freq(j)+2._DP<=1.E-15_DP .or. &
        table_freq(j)<table_freq(j-1)) then
!
        rescue_flag(i)=.true.
        f_zero_tmp(i)=-2._DP
     else
        if (i<=1.and.first) then
           call get_freqp (table_freq(j-1),table_freq(j),wp,erreur)
        else
           call get_freqp2(table_freq(j-1),table_freq(j),wp,erreur)
        endif
!        print*,'Frequence trouvee numero:',i,sngl(wp/2._DP/PI)
        f_zero_tmp(i)=wp/2._DP/PI
     endif
  enddo
  call check_and_rescue(l,table_freq,rescue_flag,f_zero_tmp,f_zero &
                           ,nmin,n_zero,wmin,wmax)
!
  deallocate(table_freq,table_det,f_zero_tmp,rescue_flag)
!------------------------------------------------------------------
  end subroutine get_fp_modele_R
!------------------------------------------------------------------
!------------------------------------------------------------------
  subroutine get_fp_modele_L(l_in,nmax,fmin,fmax,f_zero,n_zero,first,nmin)
!------------------------------------------------------------------
    use bits, only: jcom,l,epsint,epsint_base
    use def_gparam
    use yannos_flag
    implicit none
!
    integer , intent(in) :: l_in,nmax
    integer , intent(out) :: nmin
    real(DP), intent(in) :: fmin,fmax
    logical, intent(in) :: first
    integer , intent(out):: n_zero
    real(DP), dimension(:), intent(out):: f_zero
!
    real(DP), dimension(:), pointer :: table_freq,table_det
    real(DP):: wmin,wmax,wp
    integer :: i,nbmodes_min,j,nbi
    logical :: erreur,flag
    logical,  dimension(:), allocatable :: rescue_flag
    real(DP), dimension(:), allocatable :: f_zero_tmp
!
!on stipule au module minos qu'on veux les modes de rayleigh
!
    l=l_in
    if (l==0) then
       n_zero=0
    else
       jcom=2
!
       wmin=2.d0*PI*fmin
       wmax=2.d0*PI*fmax    
!on etablit la precision d'integration a la presision de base   
    epsint=epsint_base
!
! on etablit la table de recherche en freqeunce
!
       call get_table_freq(l,nmax,wmin,wmax,table_freq,table_det &
            ,nbmodes_min,n_zero)
!-------------------------------------
       if (force_systemic_search) then
! si on veut forcer la recherch systematique 
! ca peut servir dans certain model a D" tres lent
!-------------------------------------
          if (n_zero/=0) call rech_system(l,wmin,wmax,f_zero,n_zero)
!-------------------------------------
       else
!-------------------------------------
       nbi=0
       do while (any(table_freq(:)<=-2._DP).and.l/=1.and.nbi<2) 
          print*,'Je suis oblige d''incrementer la precision d''integration pour L l=',l
          epsint=epsint*1.d-2
          call get_table_freq(l,nmax,wmin,wmax,table_freq,table_det &
               ,nbmodes_min,n_zero) 
          nbi=nbi+1
       enddo
       epsint=epsint_base
       nmin=max(nbmodes_min,0)
       n_zero=n_zero-nmin
       allocate(rescue_flag(1:n_zero),f_zero_tmp(1:n_zero+100))
       f_zero_tmp(:)=-10._DP
       f_zero    (:)=-10._DP
       rescue_flag(:)=.false.
       if (n_zero>size(f_zero)) stop 'neominos: probleme de taille: increase Nmax'
! recherche des frequences:
       do i=1,n_zero
          j=i+nbmodes_min
          if (table_freq(j-1)+2._DP<=1.E-15_DP.or.table_freq(j)+2._DP<=1.E-15_DP .or. &
               table_freq(j)<table_freq(j-1)) then
!
             rescue_flag(i)=.true. 
             f_zero_tmp(i)=-2._DP
          else
             if (i<=1.and.first) then
                call get_freqp (table_freq(j-1),table_freq(j),wp,erreur)
             else
                call get_freqp2(table_freq(j-1),table_freq(j),wp,erreur)
             endif
             if (erreur) stop 'There is a  problem for this frequency'
             f_zero_tmp(i)=wp/2._DP/PI
          endif
       enddo
       call check_and_rescue(l,table_freq,rescue_flag,f_zero_tmp,f_zero &
            ,nmin,n_zero,wmin,wmax)
  !
!-------------------------------------
       endif
!-------------------------------------
       deallocate(table_freq,table_det)
    endif
!------------------------------------------------------------------
  end subroutine get_fp_modele_L
!------------------------------------------------------------------

!------------------------------------------------------------------
  subroutine baisse_wmax(l,nmodemin,nmodemax,nmax,wmax)
!------------------------------------------------------------------
    use def_gparam
    implicit none
    integer, intent(in) :: l,nmodemin,nmodemax,nmax
    real(DP), intent(inout) :: wmax
!
    real(DP) :: bid,wt,wt1
    integer  :: n
!
    n=nmodemax-nmodemin
    if (n>nmax) then
       wt1=wmax
       wt =wmax/2._DP
       do while(n/=nmax)
          call compte_nb_modes(l,wt,n,bid)          
          n=n-nmodemin
          if (n>nmax) then
             wt1=wt
             wt=wt/2._DP
          else if (n<nmax) then
             wt=(wt+wt1)/2._DP
          endif
          if (abs(wt)<1.E-6_DP .or. abs(wt-wmax)/wmax<1.E-6_DP) stop 'baisse_wmax: &
                      & erreur fatale, ca va pas du tout'
       enddo
    else
       print*,'baisse_wmax: appelle inutile!'
    endif
!    print*,'done',wmax,wt
    wmax=wt    
!------------------------------------------------------------------
  end subroutine baisse_wmax
!------------------------------------------------------------------

!------------------------------------------------------------------
  subroutine monte_wmin(l,nmodemin,wmin,wmax)
!------------------------------------------------------------------
    use def_gparam
    implicit none
    integer, intent(in) :: l
    integer, intent(inout) :: nmodemin
    real(DP), intent(in) :: wmax
    real(DP), intent(inout) :: wmin
!
    real(DP) :: bid,wt,wt1
    integer  :: n
!
    n=nmodemin
    if (n<0) then
       wt1=wmin
       wt =(wmax+wmin)/2._DP
       do while(n/=0)
          call compte_nb_modes(l,wt,n,bid)          
          if (n<0) then
             wt1=wt
             wt=2._DP*wt
          else if (n>0) then
             wt=(wt+wt1)/2._DP
          endif
          if (abs(wt)<1.E-6_DP .or. abs(wt-wmax)/wmax<1.E-6_DP) stop 'baisse_wmax: &
                      & erreur fatale, ca va pas du tout'
       enddo
    else
       print*,'monte_wmin: appelle inutile!'
    endif
    print*,'done',wmin,wt
    wmin=wt 
    nmodemin=0
!------------------------------------------------------------------
  end subroutine monte_wmin
!------------------------------------------------------------------

!------------------------------------------------------------------
  subroutine get_freqp2(wmin,wmax,wp,erreur,det1_in,det2_in)
! cherche l'unique frequence wp propre presente entre wmin et wmax
! avec une precision de epsi
!NOTE: Il faut etre passe par get_table_freq AVANT!
!Utilise un iterpolation par spline
!------------------------------------------------------------------
    use bits
    use minos
    use shanks, only: maxo,stepf
    use module_spline
    implicit none
    doubleprecision, intent(in) :: wmin,wmax    
    doubleprecision, intent(out):: wp
    logical, intent(out) :: erreur
    doubleprecision, optional, intent(in) :: det1_in,det2_in
!
    integer, parameter :: sens=0,NBIMAX=5000
    integer :: ibid,nb_iter,nbv,dv
    doubleprecision :: w1,w2,detp,det1,a,b,det2,dw,ypn,yp1,wnorme,wp1,drv
    doubleprecision, dimension(NBIMAX) :: xtab,ytab,ypptab
!
    knsw=0 !kount swich off
    maxo=8 ! precision rks ?
    if (l<4) then
       stepf=0.5d0
    else
       stepf=1.0d0
    endif
!
!
    if ((wmax-wmin)<eps) then
       print*,'wmin, wmax,l:',wmin,wmax,l
       stop 'get_freqp2: wmin doit etre < wmax'
    endif
    w1=wmin
    w2=wmax
    wnorme=max((wmin+wmax)/2.d0,5.d-3)
    nb_iter=0
    erreur=.false.
    if (present(det1_in)) then
       det1=det1_in
    else
       call detqn(w1,ibid,det1,sens)
    endif
    if (present(det2_in)) then
       det2=det2_in
    else
       call detqn(w2,ibid,det2,sens)
    endif
!
!
    wp=(w1+w2)/2.d0
    call detqn(wp,ibid,detp,sens)
    xtab(1)=w1
    xtab(2)=wp
    wp1    =wp
    xtab(3)=w2
    ytab(1)=det1
    ytab(2)=detp
    ytab(3)=det2
    nbv=3
!
    nb_iter=0   
!initialisation pour que ca passe au premier passage: 
    dw=wnorme*eps*1.d10
!
    do while (abs(dw/wnorme)>eps) 
       nb_iter=nb_iter+1
       call derive_bords(xtab,ytab,yp1,ypn,nbv)
       call spline(xtab(1:nbv),ytab(1:nbv),yp1,ypn,ypptab(1:nbv))
       call dicho_spline(xtab,ytab,ypptab,nbv,eps,wp,erreur)
       if (erreur) stop 'erreur apres dicho_spline!'
       call detqn(wp,ibid,detp,sens)
       nbv=nbv+1
       if (nbv>NBIMAX) stop 'get_freqp2: Oooppps nbv>NBIMAX!'
       dw=minval(abs(wp-xtab(:)))
       call add_value(xtab,ytab,wp,detp,nbv)
       wp1=wp
    enddo
!    print*,'get_freqp2: nb_iter:',nb_iter
!
!------------------------------------------------------------------
  end subroutine get_freqp2
!------------------------------------------------------------------

!------------------------------------------------------------------
  subroutine get_freqp3(wmin,wmax,wp,erreur,det1_in,det2_in)
! cherche l'unique frequence wp propre presente entre wmin et wmax
! avec une precision de epsi
!NOTE: Il faut etre passe par get_table_freq AVANT!
!Utilise un iterpolation par spline
!------------------------------------------------------------------
    use bits
    use minos
    use shanks, only: maxo,stepf
    use module_spline
    implicit none
    doubleprecision, intent(in) :: wmin,wmax    
    doubleprecision, intent(out):: wp
    logical, intent(out) :: erreur
    doubleprecision, optional, intent(in) :: det1_in,det2_in
!
    integer, parameter :: sens=0,NBIMAX=5000
    integer :: ibid,nb_iter,nbv,dv,nbv_max
    doubleprecision :: w1,w2,detp,det1,a,b,det2,dw,ypn,yp1,wnorme,wp1,drv
    doubleprecision, dimension(NBIMAX) :: xtab,ytab,ypptab
!
    knsw=0 !kount swich off
    maxo=8 ! precision rks ?
    nbv_max=5
    if (l<4) then
       stepf=0.5d0
    else
       stepf=1.0d0
    endif
!
!
    if ((wmax-wmin)<eps) then
       print*,'wmin, wmax,l:',wmin,wmax,l
       stop 'get_freqp2: wmin doit etre < wmax'
    endif
    w1=wmin
    w2=wmax
    wnorme=max((wmin+wmax)/2.d0,5.d-3)
    nb_iter=0
    erreur=.false.
    if (present(det1_in)) then
       det1=det1_in
    else
       call detqn(w1,ibid,det1,sens)
    endif
    if (present(det2_in)) then
       det2=det2_in
    else
       call detqn(w2,ibid,det2,sens)
    endif
!
!
    wp=(w1+w2)/2.d0
    call detqn(wp,ibid,detp,sens)
    xtab(1)=w1
    xtab(2)=wp
    wp1    =wp
    xtab(3)=w2
    ytab(1)=det1
    ytab(2)=detp
    ytab(3)=det2
    nbv=3
!
    nb_iter=0   
!initialisation pour que ca passe au premier passage: 
    dw=wnorme*eps*1.d10
!
    do while (abs(dw/wnorme)>eps) 
       nb_iter=nb_iter+1
       call derive_bords(xtab,ytab,yp1,ypn,nbv)
       call spline(xtab(1:nbv),ytab(1:nbv),yp1,ypn,ypptab(1:nbv))
       call dicho_spline(xtab,ytab,ypptab,nbv,eps,wp,erreur)
       if (erreur) stop 'erreur apres dicho_spline!'
       call detqn(wp,ibid,detp,sens)
       if (nbv<nbv_max) then
          nbv=nbv+1
          if (nbv>NBIMAX) stop 'get_freqp2: Oooppps nbv>NBIMAX!'
          dw=minval(abs(wp-xtab(:)))
          call add_value(xtab,ytab,wp,detp,nbv)
       else
          dw=minval(abs(wp-xtab(:)))
          call add_value2(xtab,ytab,wp,detp,nbv)
       endif
       wp1=wp
    enddo
    print*,'get_freqp2: nb_iter:',nb_iter
!
!------------------------------------------------------------------
  end subroutine get_freqp3
!------------------------------------------------------------------

!------------------------------------------------------------------
  subroutine ANN_for_splineR(l,nz,fp,nbv,fc,annc,anns,annp,fmax,eps)
!prepare pour calculer ANN en une frequence qq qvec une interpolation 
!par spline
! l: ordre angulaire
! nz: nombre de frequences propres
! fp: frequence propre
! nbv: nombre de frequences stockees entre chaques frequences propres
!     nbv > 5!
! annc(nbv,n) (n de 1 a 4 pour A00,A01,A10 et A11) valeur au nbv points
! anns(nbv,n) derive aux points pour les splines
! fc(nz) : frequences
! fmax: freq max. elle doit correspondre a celle utilise pour caluler annp
!
! Note il faut ne pas avoir rate de frequence propre (proche de 0)
!------------------------------------------------------------------
    use def_gparam
    use module_spline
!
    integer,                       intent(in) :: l,nz,nbv
    real(DP),                      intent(in) :: eps,fmax
    real(DP), dimension(nz,4)    , intent(in) :: annp
    real(DP), dimension(nz)      , intent(in) :: fp
    real(DP), dimension(nbv)     , intent(out):: fc
    real(DP), dimension(nbv,4)   , intent(out):: annc,anns
!
    integer :: nvc,max0
    real(DP) :: A00,A01,A10,A11,df
    real(DP), dimension(4) :: y1p,ynp
    if (nbv<5) stop 'ANN_for_splineR: nbv<5!'
    df=fmax/real(nbv-1,DP)
!==================
    do nvc=1,nbv
!==================
       fc(nvc)=(nvc-1)*df
       if (nvc==1) then
          call ANNcR0(l,annc(1,:),y1p)
          do i=1,4
             annc(1,i)=annc(1,i)-sum_apsom (fc(nvc),nz,fp,annp(:,i),eps)
!derivee au bord [:
             y1p(   i)=y1p(   i)-sum_apsomp(fc(nvc),nz,fp,annp(:,i),eps)
          enddo
       else
!
          if (trop_pres(fc(nvc),fp,nz,eps)) then
!fc tombe sur une frequence propre
             A00=0.0_DP
             A01=0.0_DP
             A10=0.0_DP
             A11=0.0_DP
          else
             call ANNcR(l,fc(nvc),A00,A01,A10,A11)
          endif
          annc(nvc,1)=A00-sum_apsom (fc(nvc),nz,fp,annp(:,1),eps)
          annc(nvc,2)=A01-sum_apsom (fc(nvc),nz,fp,annp(:,2),eps)
          annc(nvc,3)=A10-sum_apsom (fc(nvc),nz,fp,annp(:,3),eps)
          annc(nvc,4)=A11-sum_apsom (fc(nvc),nz,fp,annp(:,4),eps)
       endif
!==================
    enddo
!==================
!derivee au bord ]:
    do i=1,4
!       ynp(   i)=  -sum_apsomp(fc(nbv),nz,fp,annp(:,i),eps)
       ynp(   i)= (annc(nbv,i)-annc(nbv-1,i))/df
    enddo
    if (l==0) then
       max0=1
    else
       max0=4
    endif
    do i=1,max0
       call spline(fc(:),annc(:,i),y1p(i),ynp(i),anns(:,i))
    enddo
!------------------------------------------------------------------
  end subroutine ANN_for_splineR
!------------------------------------------------------------------

!------------------------------------------------------------------
  real(DP) function sum_apsom(f,n,fp,ap,eps)
!------------------------------------------------------------------
    use def_gparam
    implicit none
!
    real(DP) :: f,eps
    integer  :: n
    real(DP), dimension(n) :: fp,ap
!
    real(DP) :: a
    integer  :: i
!
    a=0.0_DP
    do i=1,n
       if (abs(f-fp(i))>eps*fp(i)) a=a+ap(i)/2._DP/PI/(f-fp(i))
    enddo
    sum_apsom=a
!------------------------------------------------------------------
  end function sum_apsom
!------------------------------------------------------------------

!------------------------------------------------------------------
  logical function trop_pres(f,fp,n,eps)
!------------------------------------------------------------------
    use def_gparam
    implicit none
    real(DP) :: f,eps
    integer :: n
    real(DP), dimension(n) :: fp
!
    integer :: i
!
    trop_pres=.false.
    do i=1,n
       if (abs(fp(i))<1.E-15_DP) stop 'function trop_pres: &
            &normalement 0 n''est pas vp!'
         if (abs(f-fp(i))< fp(i)*eps) trop_pres=.true.
    enddo
!    
!------------------------------------------------------------------
  end function trop_pres
!------------------------------------------------------------------

!------------------------------------------------------------------
  real(DP) function sum_apsomp(f,n,fp,ap,eps)
!------------------------------------------------------------------
    use def_gparam
    implicit none
!
    real(DP) :: f,eps
    integer  :: n
    real(DP), dimension(n) :: fp,ap
!
    real(DP) :: a   
    integer  :: i
!
    a=0.0_DP
    do i=1,n
       if (abs(f-fp(i))>eps*fp(i)) a=a-ap(i)/2._DP/PI/(f-fp(i))**2
    enddo
    sum_apsomp=a
!------------------------------------------------------------------
  end function sum_apsomp
!------------------------------------------------------------------

!------------------------------------------------------------------
  subroutine ANN_for_splineL(l,nz,fp,nbv,fc,a5c,a5s,a5p,fmax,eps)
!------------------------------------------------------------------
!prepare pour calculer ANN en une frequence qq qvec une interpolation 
!par spline
! l: ordre angulaire
! nz: nombre de frequences propres
! fp: frequence propre
! nbv: nombre de frequences stockees entre chaques frequences propres
!     nbv > 5!
! a5c(nbv)  valeur au nbv points
! a5s(nbv) derive aux points pour les splines
! fc(nz) : frequences
! fmax: freq max. elle doit correspondre a celle utilise pour caluler a5p
!
! Note il faut ne pas avoir rate de frequence propre (proche de 0)
!------------------------------------------------------------------
    use def_gparam
    use module_spline
!
    integer,                       intent(in) :: l,nz,nbv
    real(DP),                      intent(in) :: eps,fmax
    real(DP), dimension(nz)      , intent(in) :: a5p
    real(DP), dimension(nz)      , intent(in) :: fp
    real(DP), dimension(nbv)     , intent(out):: fc
    real(DP), dimension(nbv)     , intent(out):: a5c,a5s
!
    integer :: nvc
    real(DP) :: y1p,ynp,A5,df
    if (nbv<5) stop 'ANN_for_splineL: nbv<5!'
    df=fmax/(nbv-1)
!+++++++++++++++++++
    if (l/=0) then
!+++++++++++++++++++
!==================
    do nvc=1,nbv
!==================
       fc(nvc)=(nvc-1)*df
       if (nvc==1) then
          call ANNcL0(l,a5c(1),y1p)
          a5c(1)=a5c(1)-sum_apsom (fc(nvc),nz,fp,a5p,eps)

!          y1p   =y1p   -sum_apsomp(fc(nvc),nz,fp,a5p,eps)          
       else
!
          if (minval(abs(fp(:)-fc(nvc))/fc(nvc))<eps) then
!fc tombe sur une frequence propre
             A5=0.0_DP
          else
             call ANNcL(l,fc(nvc),A5)
          endif
          a5c(nvc)=A5-sum_apsom (fc(nvc),nz,fp,a5p,eps)
       endif
!==================
    enddo
!==================
!derivee au bord [:
    y1p   =(a5c(2)-a5c(1))/df
!derivee au bord ]:
!    ynp=-sum_apsomp(fc(nbv),nz,fp,a5p,eps)
    ynp=(a5c(nbv)-a5c(nbv-1))/df
    call spline(fc,a5c,y1p,ynp,a5s)
!++++++++++++++++++
    else
!cas l=0
!++++++++++++++++++
      a5c(:)=0.0_DP
      a5s(:)=0.0_DP
!+++++++++++++++++++
    endif
!+++++++++++++++++++
!------------------------------------------------------------------
  end subroutine ANN_for_splineL
!------------------------------------------------------------------

!------------------------------------------------------------------
  subroutine derive_bords(x,y,yp1,ypn,n)
!------------------------------------------------------------------
    implicit none
    integer, intent(in) :: n
    doubleprecision, dimension(n), intent(in ):: x,y
    doubleprecision,intent(out):: yp1,ypn
!
!
    if (n<2) stop 'derive_bords : nbv doit etre >= a 2'
!
    yp1=(y(2)-y(1))/(x(2)-x(1))
    ypn=(y(n)-y(n-1))/(x(n)-x(n-1))
!------------------------------------------------------------------
  end subroutine derive_bords
!------------------------------------------------------------------
!------------------------------------------------------------------
  subroutine derive(x,y,yp,n)
!------------------------------------------------------------------
    implicit none
    integer, intent(in) :: n
    doubleprecision, dimension(n), intent(in ):: x,y
    doubleprecision, dimension(n), intent(out):: yp
!
    integer :: i
!
    if (n<2) stop 'derive : nbv doit etre >= a 2'
!
    yp(1)=(y(2)-y(1))/(x(2)-x(1))
    yp(n)=(y(n)-y(n-1))/(x(n)-x(n-1))
!
    if (n>2) then
       do i=2,n-1
          yp(i)=(y(i+1)-y(i-1))/(x(i+1)-x(i-1))
       enddo
    endif
!------------------------------------------------------------------
  end subroutine derive
!------------------------------------------------------------------

!------------------------------------------------------------------
  subroutine add_value(x,y,xp,yp,n)
!incere (xp,yp) dans le tableau (x,y)
!------------------------------------------------------------------
    implicit none
    integer, intent(in) :: n
    doubleprecision, dimension(n) , intent(inout) :: x,y
    doubleprecision :: xp,yp
!
    doubleprecision, dimension(n) :: x_tmp, y_tmp
    integer :: i
!
    x_tmp(:)=x(:) 
    y_tmp(:)=y(:)
    i=1
    do while(x_tmp(i) < xp)
       i=i+1
       if (i>n) stop 'add_value: xp n''est pas dans les bornes de x'
    enddo
    i=i-1
    x(1:i)=x_tmp(1:i)
    y(1:i)=y_tmp(1:i)
    x(i+1)=xp
    y(i+1)=yp
    x(i+2:n)=x_tmp(i+1:n-1)
    y(i+2:n)=y_tmp(i+1:n-1)
!------------------------------------------------------------------
  end subroutine add_value
!------------------------------------------------------------------

!------------------------------------------------------------------
  subroutine add_value2(x,y,xp,yp,n)
!incere (xp,yp) dans le tableau (x,y)
!------------------------------------------------------------------
    implicit none
    integer, intent(in) :: n
    doubleprecision, dimension(n) , intent(inout) :: x,y
    doubleprecision :: xp,yp
!
    doubleprecision, dimension(n) :: x_tmp, y_tmp
    integer :: i
!
    x_tmp(:)=x(:) 
    y_tmp(:)=y(:)
    i=1
    do while(x_tmp(i) < xp)
       i=i+1
       if (i>n) stop 'add_value2: xp n''est pas dans les bornes de x'
    enddo
    i=i-1
    if (i>n/2) then
!on pert le premier point
       x(1:i-1)=x_tmp(2:i)
       y(1:i-1)=y_tmp(2:i)
       x(i)=xp
       y(i)=yp
       x(i+1:n)=x_tmp(i+1:n)
       y(i+1:n)=y_tmp(i+1:n)
    else
!on pert le dernier point
       x(1:i)=x_tmp(1:i)
       y(1:i)=y_tmp(1:i)
       x(i+1)=xp
       y(i+1)=yp
       x(i+2:n)=x_tmp(i+1:n-1)
       y(i+2:n)=y_tmp(i+1:n-1)
    endif
!------------------------------------------------------------------
  end subroutine add_value2
!------------------------------------------------------------------

!------------------------------------------------------------------
  subroutine get_driv(x,y,n,yp)
!------------------------------------------------------------------
    implicit none
    integer, intent(in) :: n
    doubleprecision, dimension(n) , intent(in) :: x,y
    doubleprecision, intent(out) :: yp
!
    integer :: i
!
    i=1
    do while(y(i)*y(i+1) > 0.0_DP.and.i<n-1)
       i=i+1
    enddo
    i=i-1
    yp=(y(i+1)-y(i))/max((x(i+1)-x(i)),1.E-16_DP)
!------------------------------------------------------------------
  end subroutine get_driv
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine dicho_spline(x,y,ypp,n,eps,xp,erreur)
!dichotomie sur le spline:
!------------------------------------------------------------------
  use module_spline
  implicit none
!
  integer, intent(in) :: n
  doubleprecision, dimension(n) , intent(in) :: x,y,ypp
  doubleprecision,intent(in)  :: eps
  doubleprecision,intent(out) :: xp
  logical, intent(out) :: erreur
!
  doubleprecision :: x1,x2,xx,yp,y1
  integer :: nb_iter
!
  x1=x(1)
  x2=x(n)
  xx=(x1+x2)
  nb_iter=0
  y1=y(1)
!
  do while(abs((x1-x2)/xx)>eps)
     nb_iter=nb_iter+1 
     xp=(x1+x2)/2.d0
     yp=splint(x,y,ypp,xp)
     if (y1*yp<0.d0) then
        x2=xp
     else
        x1=xp
        y1=yp
     endif
  enddo
  if (abs(x(1)-xp)/xx<eps.or.abs(x(n)-xp)/xx<eps) then
     erreur=.true.
  else
     erreur=.false.
  endif
  xp=(x1+x2)/2.d0
!------------------------------------------------------------------
end subroutine dicho_spline
!------------------------------------------------------------------

!------------------------------------------------------------------
  subroutine get_freqp(wmin,wmax,wp,erreur)
! cherche l'unique frequence wp propre presente entre wmin et wmax
! avec une precision de epsi
!NOTE: Il faut etre passe par get_table_freq AVANT!
! dichotomie simple
!------------------------------------------------------------------
    use bits
    use minos
    use shanks, only: maxo
!
    implicit none
    doubleprecision, intent(in) :: wmin,wmax    
    doubleprecision, intent(out):: wp
    logical, intent(out) :: erreur
!
    integer, parameter :: sens=0
    integer :: ibid,nb_iter,i
    doubleprecision :: w1,w2,detp,det1,det2
!
    maxo=8 ! precision rks ?
    knsw=0 !kount swich off
!
    if ((wmax-wmin)<eps) then
       print*,'wmin, wmax,l:',wmin,wmax,l
       stop 'get_freqp: wmin doit etre < wmax'
    endif
    w1=wmin
    w2=wmax
!
    nb_iter=0
    erreur=.false.
    call detqn(w1,ibid,det1,sens)
    do while(abs((w1-w2)/(wmax+wmin))>eps.and..not.erreur)
       nb_iter=nb_iter+1
       wp=(w1+w2)/2.d0
       if (abs(wmin-wp)<eps.or.abs(wmax-wp)<eps) then
          print*,'get_freqp: il n''y a pas de frequence propre dans [fmin,fmax]=',wmin/2./pi,wmax/2./pi
          erreur=.true.
       endif
       if (.not.erreur) then
          call detqn(wp,ibid,detp,sens)
          if (detp*det1<0.d0) then
             w2=wp
             det2=detp
          else
             w1=wp
             det1=detp
          endif
       endif
    enddo
    wp=(w1+w2)/2.0d0

!------------------------------------------------------------------
  end subroutine get_freqp
!------------------------------------------------------------------
!------------------------------------------------------------------
  subroutine get_table_freq(l_in,nmax,wmin,wmax_,table_freq,table_det &
                             ,nbmodes_min,nbmodes_max)
!
!nmax: nombre de branches maxi a chercher
!
! pour l_in cherches la table frequences
! nbmodes_min: nbmodes_min que l'on va rater
! nbmodes_max-nbmodes_min: nombres de modes present
! pour tout i, il y a une et une seule frequence propre entre
! table_freq(i) et table_freq(i+1)
!------------------------------------------------------------------
    use bits , only: jcom,eps
    implicit none    
    integer, intent(in)  :: l_in,nmax
    doubleprecision, intent(inout) :: wmin
    doubleprecision, intent(in) ::wmax_
    doubleprecision, dimension(:), pointer :: table_freq,table_det
    integer, intent(out)  :: nbmodes_min,nbmodes_max
!
    integer :: i1,i2,n,nmodes,nbmodes_found,num,l,i,j,nb_not_found  &
                    ,nbmodes_min_resc,nbmodes_max_resc,ic
    doubleprecision :: w1,w2,det1,det2,wt,wmax,ww,wmin_resc,wmax_resc ,wmin1
    doubleprecision, dimension(:), pointer :: table_freq_resc,table_det_resc
    logical :: notfound,il_en_reste
    
!
    l=l_in
    wmax=wmax_
!test
!   call get_det_curve(l,wmin,wmax)
!stop 
!

    call compte_nb_modes(l,wmin,nbmodes_min,det1)
    if (nbmodes_min<0.and.l>2) call monte_wmin(l,nbmodes_min,wmin,wmax)
    if (nbmodes_min>0.and.l/=1) print*,'fmin est trop grandes pour determiner&
         & tous les modes'
    call compte_nb_modes(l,wmax,nbmodes_max,det2)

    if (nbmodes_max-max(nbmodes_min,0)>nmax) then
!on cherches wmax tel que nbmodes_found<=nmax:
       call baisse_wmax(l,max(nbmodes_min,0),nbmodes_max,nmax,wmax)
       call compte_nb_modes(l,wmax,nbmodes_max,det2)
    endif


!
! nombres de frequences propres a trouver:
!
    nbmodes_found=max(nbmodes_max-nbmodes_min,0)
    if (jcom==2) then
       write(*,144) nbmodes_found,l
    else
       write(*,143) nbmodes_found,l
    endif
143 format('There is/are ',i4,' frequencies to search for l=',i4,' (S))')
144 format('There is/are ',i4,' frequencies to search for l=',i4,' (T))')
!
    allocate(table_freq(-2:nbmodes_max+2),table_det(-2:nbmodes_max+2))
    table_freq(:)=-1.d0
    table_freq(nbmodes_min)=wmin
    table_freq(nbmodes_max)=wmax
    table_det (nbmodes_min)=det1
    table_det (nbmodes_max)=det2
!
! on remplit la table de frequence:
! il faut chercher nbmodes_found-1 frequences qui separenent les
! nbmodes_found freq propres a trouver
!
    n=0
    if (nbmodes_found>0) then
       il_en_reste=.true.
    else
       il_en_reste=.false.
    endif
!--------------
    do while(il_en_reste)
!--------------
!
! cherche la premiere case vide; la valeur -2 est un code pour une vp 
! non reslue et non resolvable
!
       i=nbmodes_min
       do while ((table_freq(i)>0.0d0.or.(abs(table_freq(i)+2.d0)<1.d-10)) &
                 .and.i<nbmodes_max )
          i=i+1
          if (i>nbmodes_max-1) then
             il_en_reste=.false.
          else
             il_en_reste=.true.
          endif
       enddo 
       i1=i-1
       do while (abs(table_freq(i1)+2.d0)<1.d-10)
          i1=i1-1
       enddo

       if (i1<0) then
!          stop 'erreur dans la table de freq: fatale 2'
       endif
!borne inferieur de recherche:
       w1=table_freq(i1)
!++++++++++++++++++++++++++++
       if (il_en_reste) then
!++++++++++++++++++++++++++++
!
! cherche la prochaine case pleine
!
       i=i1+1
       do while (table_freq(i)<0.0d0)
          i=i+1
          if (i>nbmodes_max) stop 'erreur dans la table de freq: fatale 1'
       enddo
       i2=i
!       if (i<=1) stop 'erreur dans la table de freq: fatale 2'
!borne superieur de recherche:
       w2=table_freq(i2)
!on boucle jusqu'a ce qu'on est trouver la frequnce:
       notfound=.true.
!
!+++++++++++++++++++++
       do while(notfound)
!+++++++++++++++++++++
!on essaie la frequence suivante:
          wt=(w1+w2)/2
          ww=max(w2,wt)
          if (abs(w2/ww-wt/ww)<eps) then
!il y a un probleme
!             print*,'Je ne peux pas separer ces deux frequnces: '
!             print*,'Il faut augmenter la precision (eps)'
!             print*,'Ca peut  venir du mode kount et dans ce cas &
!                   & augmenter le nombre de couche du modele peut etre la solution.'
             do i=i1+1,i2-1
                table_freq(i)=-2.d0
                n=n+1
             enddo
             n=n-1
             notfound=.false.
          else
             if (wt<=0.0_DP) then
                stop 'bordel a queue ... dans get_table_freq' 
             endif
             call compte_nb_modes(l,wt,num,det1)
if (num<nbmodes_min.or.num>i2) then
wt=w1+(w2-w1)/3.; call compte_nb_modes(l,wt,num,det1)
endif
             if (num<i1.or.num>i2) then
!Le mode count a deconne,
!                print*,'pb du mode kount, il faut probablement &
!                     &augmenter le nombre de couches du modele (l=',l,')',num
!                print*,'On essaie un sauvetage...'
                if (num>i2) then
                   table_freq(i2)=-1.d0
                   i2=num
                   w2=wt
                else
                   table_freq(i1)=-1.d0
                   i1=num
                   w1=wt
                endif
                if (num<-2) stop 'get_table_freq: num<-2'
                if (num<nbmodes_min) nbmodes_min=num
                if (num>nbmodes_max+2) then
                   stop 'get_table_freq: num>nbnodes_max+2'
                endif
                if (num>nbmodes_max) nbmodes_max=num
                table_freq(num)=wt	
                table_det (num)=det1
             else if (num==i2) then
                w2=wt
             else if (num==i1) then
                w1=wt
             else
!freq p trouvee:
!test
!                print*,num,wt
!fin test
                notfound=.false.
                table_freq(num)=wt
                table_det (num)=det1
             endif
          endif
!+++++++++++++++++++++
       enddo
!+++++++++++++++++++++           
       n=n+1
!++++++++++++++++++++++++++++
    endif
!++++++++++++++++++++++++++++
!--------------
    enddo
!--------------
!on enleve  les trous
    nb_not_found=0
    do i=nbmodes_min,nbmodes_max
       if (abs(table_freq(i)+2.d0)<1.d-15) then  
!!$          wmin_resc=(table_freq(max(i-2,nbmodes_min)) &
!!$                    +table_freq(max(i-1,nbmodes_min)))/2.d0
!!$          wmax_resc=(table_freq(min(i+3,nbmodes_max)) &
!!$                    +table_freq(min(i+2,nbmodes_max)))/2.d0
!!$          call rescue_table_freq(l_in,nmax,wmin_resc,wmax_resc &
!!$                             ,table_freq_resc,table_det_resc &
!!$                             ,nbmodes_min_resc,nbmodes_max_resc)
!!$          do j=nbmodes_min_resc,nbmodes_max_resc
!!$             if ((table_freq(j)+2.d0)<1.d-15) then
!!$                if ((table_freq_resc(j)<1.d-15)) then
!!$                   print*,'rescue failed pour l=',l_in
!!$                   nb_not_found=nb_not_found+1
!!$                else
!!$                   table_freq(j)=table_freq_resc(j)
!!$                endif
!!$             endif
!!$          enddo
!!$          deallocate(table_freq_resc,table_det_resc)
       endif
    enddo
    nbmodes_max=nbmodes_max- nb_not_found
    if (nb_not_found/=0) then
       write(*,100) nb_not_found,l
       print*,'Une solution a ce probleme peut etre d''augmenter &
              & le nombre de couche du modele'
       print*,'ou de diminuer stepf pour ce l!' 
    endif
100 format(i2,' freq propres oubliees pour l=',i4)
!------------------------------------------------------------------
end subroutine get_table_freq
!------------------------------------------------------------------

!------------------------------------------------------------------
  subroutine rescue_table_freq(l_in,nmax,wmin,wmax_,table_freq,table_det &
                             ,nbmodes_min,nbmodes_max)
!
!nmax: nombre de branches maxi a chercher
!
! pour l_in cherches la table frequences
! nbmodes_min: nbmodes_min que l'on va rater
! nbmodes_max-nbmodes_min: nombres de modes present
! pour tout i, il y a une et une seule frequence propre entre
! table_freq(i) et table_freq(i+1)
!------------------------------------------------------------------
    use bits , only: jcom,eps
    implicit none    
    integer, intent(in)  :: l_in,nmax
    doubleprecision, intent(inout) :: wmin
    doubleprecision, intent(in) ::wmax_
    doubleprecision, dimension(:), pointer :: table_freq,table_det
    integer, intent(out)  :: nbmodes_min,nbmodes_max
!
    integer :: i1,i2,n,nmodes,nbmodes_found,num,l,i,j,nb_not_found
    doubleprecision :: w1,w2,det1,det2,wt,wmax,ww
    logical :: notfound,il_en_reste
    
!
    l=l_in
    wmax=wmax_
!
! determine combient il y a de frequences propres a chercer:
!  
    call compte_nb_modes(l,wmin,nbmodes_min,det1)
    call compte_nb_modes(l,wmax,nbmodes_max,det2)
    print*,'appelle dr rescue_table_freq pour ',l_in,nbmodes_min,nbmodes_max
!
! nombres de frequences propres a trouver:
!
    nbmodes_found=max(nbmodes_max-nbmodes_min,0)
    if (jcom==2) then
       write(*,144) nbmodes_found,l
    else
       write(*,143) nbmodes_found,l
    endif
143 format('There is/are ',i4,' eigenfrequencies to search for rescue l=',i4,' (S))')
144 format('There is/are ',i4,' eigenfrequencies to search for rescue l=',i4,' (T))')
!
    allocate(table_freq(-2:nbmodes_max+2),table_det(-2:nbmodes_max+2))
    table_freq(:)=-1.d0
    table_freq(nbmodes_min)=wmin
    table_freq(nbmodes_max)=wmax
    table_det (nbmodes_min)=det1
    table_det (nbmodes_max)=det2
!
! on remplit la table de frequence:
! il faut chercher nbmodes_found-1 frequences qui separenent les
! nbmodes_found freq propres a trouver
!
    n=0
    if (nbmodes_found>0) then
       il_en_reste=.true.
    else
       il_en_reste=.false.
    endif
!--------------
    do while(il_en_reste)
!--------------
!
! cherche la premiere case vide; la valeur -2 est un code pour une vp 
! non reslue et non resolvable
!
       i=nbmodes_min
       do while ((table_freq(i)>0.0d0.or.(abs(table_freq(i)+2.d0)<1.d-10)) &
                 .and.i<nbmodes_max )
          i=i+1
          if (i>nbmodes_max-1) then
             il_en_reste=.false.
          else
             il_en_reste=.true.
          endif
       enddo
       i1=i-1
       if (i1<0) then
!          stop 'erreur dans la table de freq: fatale 2'
       endif
!borne inferieur de recherche:
       w1=table_freq(i1)
!++++++++++++++++++++++++++++
       if (il_en_reste) then
!++++++++++++++++++++++++++++
!
! cherche la prochaine case pleine
!
       i=i1+1
       do while (table_freq(i)<0.0d0)
          i=i+1
          if (i>nbmodes_max) stop 'erreur dans la rescue table de freq: fatale 1'
       enddo
       i2=i
!       if (i<=1) stop 'erreur dans la table de freq: fatale 2'
!borne superieur de recherche:
       w2=table_freq(i2)
!on boucle jusqu'a ce qu'on est trouver la frequnce:
       notfound=.true.
!
!+++++++++++++++++++++
       do while(notfound)
!+++++++++++++++++++++
!on essaie la frequence suivante:
          wt=(w1+w2)/2
          ww=max(w2,wt)
          if (abs(w2/ww-wt/ww)<eps) then
!il y a un probleme
             print*,'Je ne peux pas separer ces deux frequnces: '
             print*,'Il faut augmenter la precision (eps)'
             print*,'Ca peut  venir du mode kount et dans ce cas &
                   & augmenter le nombre de couche du modele peut etre la solution.'
             do i=i1+1,i2-1
                table_freq(i)=-2.d0
                n=n+1
             enddo
             n=n-1
             notfound=.false.
          else
             call compte_nb_modes(l,wt,num,det1)
             if (num<i1.or.num>i2) then
!Le mode count a deconne,
!                print*,'pb du mode kount, il faut probablement &
!                     &augmenter le nombre de couches du modele (l=',l,')'
!                print*,'On essaie un sauvetage...'
                if (num>i2) then
                   table_freq(i2)=-1.d0
                   i2=num
                   w2=wt
                else
                   table_freq(i1)=-1.d0
                   i1=num
                   w1=wt
                endif
                if (num<-2) stop 'rescue_table_freq: num<-2'
                if (num<nbmodes_min) nbmodes_min=num
                if (num>nbmodes_max+2) stop 'rescue_table_freq: num>nbnodes_max+2'
                if (num>nbmodes_max) nbmodes_max=num
                table_freq(num)=wt	
                table_det (num)=det1
             else if (num==i2) then
                w2=wt
             else if (num==i1) then
                w1=wt
             else
!freq p trouvee:
                notfound=.false.
                table_freq(num)=wt
                table_det (num)=det1
             endif
          endif
!+++++++++++++++++++++
       enddo
!+++++++++++++++++++++           
       n=n+1
!++++++++++++++++++++++++++++
    endif
!++++++++++++++++++++++++++++
!--------------
    enddo
!--------------
!on enleve  les trous
    nb_not_found=0
    do i=nbmodes_min,nbmodes_max
       if (abs(table_freq(i)+2.d0)<1.d-15) then
          nb_not_found=nb_not_found+1
!          do j=i,nbmodes_found-nb_not_found
!             table_freq(j)=table_freq(j+1)
!             table_det (j)=table_det (j+1)
!          enddo
       endif
    enddo
    nbmodes_max=nbmodes_max- nb_not_found
    if (nb_not_found/=0) then
       write(*,100) nb_not_found,l
       print*,'Une solution a ce probleme peut etre d''augmenter &
              & le nombre de couche du modele'
       print*,'ou de diminuer stepf pour ce l!' 
    endif
100 format(i2,' freq propres oubliees, rescue failed pour l=',i4 &
           ,': pas asser de couches pour integrer')
!------------------------------------------------------------------
end subroutine rescue_table_freq
!------------------------------------------------------------------


!------------------------------------------------------------------
  subroutine rescue_table_freq2(l_in,nmax,wmin,wmax_,table_freq,table_det &
                             ,nbmodes_min,nbmodes_max,delta_modes)
!recherche sustematique de modes entre wmin et wmax
!nmax: nombre de branches maxi a chercher
!
! pour l_in cherches la table frequences
! nbmodes_min: nbmodes_min que l'on va rater
! nbmodes_max-nbmodes_min: nombres de modes present
! pour tout i, il y a une et une seule frequence propre entre
! table_freq(i) et table_freq(i+1)
!------------------------------------------------------------------
    use bits , only: jcom,eps
    implicit none    
    integer, intent(in)  :: l_in,nmax
    doubleprecision, intent(inout) :: wmin
    doubleprecision, intent(in) ::wmax_
    doubleprecision, dimension(:), pointer :: table_freq,table_det
    integer, intent(out)  :: nbmodes_min,nbmodes_max,delta_modes
!
    integer, parameter :: NBF=200
    integer :: nmodes,l,i,nbmodes_found
    doubleprecision :: w1,w2,det1,det2,wt,wmax,ww,df,fc,d1,d2
    doubleprecision, dimension(:), allocatable :: ftmp1,ftmp2
    
!
    l=l_in
    wmax=wmax_
!
! determine combient il y a de frequences propres a chercer:
!  
    call compte_nb_modes(l,wmin,nbmodes_min,det1)
    call compte_nb_modes(l,wmax,nbmodes_max,det2)
    print*,'appelle dr rescue_table_freq 2 pour ',l_in,nbmodes_min,nbmodes_max
!
! nombres de frequences propres a trouver:
!
    nbmodes_found=max(nbmodes_max-nbmodes_min,0)
    if (jcom==2) then
       write(*,144) nbmodes_found,l
    else
       write(*,143) nbmodes_found,l
    endif
143 format('There is/are ',i4,' eigenfrequencies to search for rescue2 l=',i4,' (S))')
144 format('There is/are ',i4,' eigenfrequencies to search for rescue2 l=',i4,' (T))')
!
    allocate(table_freq(-2:nbmodes_max+10),table_det(-2:nbmodes_max+10) &
             ,ftmp1(nbmodes_found+10),ftmp2(nbmodes_found+10))
    table_freq(:)=-1.d0
    table_freq(nbmodes_min)=wmin
    table_freq(nbmodes_max)=wmax
    table_det (nbmodes_min)=det1
    table_det (nbmodes_max)=det2
!
! on remplit la table de frequence:
! il faut chercher nbmodes_found-1 frequences qui separenent les
! nbmodes_found freq propres a trouver
!
    df=(wmax-wmin)/(NBF)/2._DP/PI
    nmodes=0
    do i=1,NBF
       fc=wmin/2._DP/PI+(i-1)*df
       if (i==1) then
          call get_det(l,fc,d1)
       else
          d2=d1
          call get_det(l,fc,d1)
          if (d2*d1<0.0_DP) then
!mode found
             nmodes=nmodes+1
             ftmp1(nmodes)=fc-df
             ftmp2(nmodes)=fc
          endif
       endif
    enddo
    do i=1,nmodes
       if (i/=nmodes) then
          table_freq(i+nbmodes_min)=(ftmp1(i+1)+ftmp2(i))*PI
       else
          table_freq(i+nbmodes_min)=ftmp2(nmodes)*PI*2._DP
       endif
    enddo
    delta_modes=nmodes-nbmodes_max+nbmodes_min
    print*,'rescue 2 en a trouve:',nmodes,delta_modes
    nbmodes_max=nmodes+nbmodes_min
    deallocate(ftmp1,ftmp2)
!------------------------------------------------------------------
end subroutine rescue_table_freq2
!------------------------------------------------------------------

!------------------------------------------------------------------
subroutine rech_system(l,wmin,wmax,f_zero,n_zero)
!subroutine de recherche bestiale de freq propre
!------------------------------------------------------------------
  use def_gparam
  implicit none 
!
  integer, intent(in) :: l
  real(DP), intent(in):: wmin,wmax
  real(DP), dimension(:), intent(inout):: f_zero
  integer, intent(out) :: n_zero
!
  real(DP) :: df,det1,det2,bid,fmin,fmax,fc,wp,w1,w2
  integer :: nbmodes_min,nbmodes_max,nbmodes,imodes,i,nbpas
  integer, parameter :: cst_pas=50
  logical :: erreur
!
  call compte_nb_modes(l,wmin,nbmodes_min,det1)
  call compte_nb_modes(l,wmax,nbmodes_max,bid)
  nbmodes=nbmodes_max-nbmodes_min
  nbpas=nbmodes*cst_pas
  nbmodes_max=UBOUND(f_zero(:),dim=1)
!  if (nbmodes_max>UBOUND(f_zero(:),dim=1)) stop 'rech_system: nbmodes_max>ubound(f_zero)'
!  if (nbmodes>UBOUND(f_zero(:),dim=1)) stop 'rech_system: nbmodes_max>ubound(f_zero)'
!
  fmin=wmin/2._DP/PI
  fmax=wmax/2._DP/PI
  df=(fmax-fmin)/nbpas
  imodes=0
  do i=1,nbpas
     fc=i*df+fmin
     if (imodes<nbmodes_max) then 
        call get_det(l,fc,det2)
        if (det1*det2<0._DP) then
           imodes=imodes+1
           w1=(fc-df)*2._DP*PI
           w2=(fc   )*2._DP*PI
           call get_freqp2(w1,w2,wp,erreur)
           f_zero(imodes)=wp/2._DP/PI
        endif
        det1=det2
     endif
  enddo
!
  n_zero=imodes
!------------------------------------------------------------------
end subroutine rech_system 
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine check_and_rescue(l_in,table_freq,rescue_flag,f_zero_tmp,f_zero &
                           ,nmin,n_zero,wmin,wmax)
!------------------------------------------------------------------
  use def_gparam
  use yannos_flag
  implicit none 
!
  
  real(DP), dimension(:), intent(in) :: table_freq
  logical,  dimension(:), intent(inout) :: rescue_flag
  real(DP), dimension(:), intent(inout):: f_zero_tmp
  real(DP), dimension(:), intent(inout):: f_zero
  real(DP), intent(in):: wmin,wmax
  integer, intent(in) :: nmin,l_in
  integer, intent(inout) :: n_zero
!
  real(DP), dimension(n_zero) :: f_zero_tmp2,wmin_resc,wmax_resc &
                                ,wmin_resc1,wmax_resc1
!
  real(DP), dimension(:), pointer :: table_freq_resc,table_det_resc
  integer :: nmax,i,l,j,n_zero_resc,iresc,nresc
  real(DP) :: wp,fc,tmp
  logical :: erreur,flag
!
  l=l_in
  f_zero_tmp2(:)=-10._DP
  iresc=0  
!on determine les bandes de freqeunce dans lequelle on va faire une
!recherche systematique:
  wmin_resc1(:)=0.0_DP
  wmax_resc1(:)=0.0_DP

  do i=1,n_zero
     if (rescue_flag(i)) then 
        iresc=iresc+1
        call get_wmin_resc(wmin,wmin_resc1(iresc),f_zero_tmp,n_zero,i)
        call get_wmax_resc(wmax,wmax_resc1(iresc),f_zero_tmp,n_zero,i)
     endif
  enddo
!on vire les recouverment de bande:
  nresc=iresc
  n_zero_resc=0
  j=0
  if (nresc/=0) then
     wmin_resc(1)=wmin_resc1(1)
     iresc=0
     j=0
     do i=1,nresc-1
        if (wmax_resc1(i)>wmin_resc1(i+1)) then
           j=j+1
        else
           iresc=iresc+1
           wmax_resc(iresc)=wmax_resc1(iresc+j)
           wmin_resc(iresc+1)=wmin_resc1(iresc+j+1)
        endif
     enddo
     wmax_resc(iresc+1)=wmax_resc1(nresc)
     nresc=iresc+1
     if (nresc/=0.and.rescue) then
!recherche systematique:
        n_zero_resc=0
        do i=1,nresc
           print*,'Recherche systematique pour l=',l,' entre ',real(wmin_resc(i)),real(wmax_resc(i))
           call rech_system(l,wmin_resc(i),wmax_resc(i),f_zero_tmp2(n_zero_resc+1:),j)
           n_zero_resc=n_zero_resc+j
        enddo
     endif
  endif  
!
!on inclu f_zero_tmp2 dans f_zero et vire les modes de la graines si demande
!
  call mix(l,f_zero_tmp,n_zero,f_zero_tmp2,n_zero_resc,f_zero,n_zero)

!------------------------------------------------------------------
  end subroutine check_and_rescue
!------------------------------------------------------------------

!------------------------------------------------------------------
  subroutine mix(l,t1,n1,t2,n2,t3,n3)
!------------------------------------------------------------------
    use def_gparam
    implicit none
    integer, intent(in) :: l
    integer, intent(in) :: n1,n2
    integer, intent(out) :: n3
    real(DP), dimension(:), intent(in)   :: t1,t2
    real(DP), dimension(:), intent(out)  :: t3
!
    integer :: i,i1,i2
!
    i1=1
    i2=1
    t3(:)=-1.0_DP
    do while (i1<=n1 .or. i2<=n2 )
       if (i1>n1+1) stop 'mix: i1>n1+1 ???'
       if (i2>n1+1) stop 'mix: i2>n1+1 ???'
       if (t2(i2)>0.0_DP.or.i2>n2) then
          if (t1(i1)>0.0_DP.or.i1>n1) then
             if ((t1(i1)<t2(i2).or.i2>n2).and.i1<=n1) then
                call add_v(l,t1(i1),t3)
                i1=i1+1
             else
                call add_v(l,t2(i2),t3)
                i2=i2+1
             endif 
          else
             i1=i1+1
          endif
       else
          i2=i2+1
       endif
    enddo
    n3=1
    do while (t3(n3)>0.0_DP)
       n3=n3+1
    enddo
    n3=n3-1
!
!------------------------------------------------------------------
  end subroutine mix
!------------------------------------------------------------------

!------------------------------------------------------------------
  subroutine add_v(l,a,t)
!------------------------------------------------------------------
    use def_gparam
    use minos
    use yannos_flag
    implicit none
    integer, intent(in) :: l
    real(DP), intent(in) :: a
    real(DP), dimension(:), intent(inout) :: t
!
    real(DP), dimension(:),allocatable :: tmp
    integer :: i,nr,n,nb_nomodes
    logical :: deja_la,flag,plein
!
    n=UBOUND(t,dim=1)
    nb_nomodes=0
    plein=.false.
    if (a>0.0_DP) then
       nr=1
       do while (t(nr)>0.0_DP)
          nr=nr+1
          if (nr>n) then
             stop 'add_v: le tableau est plein!'
             plein=.true.
          endif
       enddo
       nr=nr-1
       if (nr==0) then
          !le tableau est vide
          t(1)=a
       else
          !la valeur est elle deja presente:      
          deja_la=.false.
          do i=1,nr
             if (abs((t(i)-a)/a)<1.E-5_DP) deja_la=.true.
          enddo
          if (.not.deja_la) then
             if (check_modes) then
                call check_eigenfreq2(l,a*2._DP*PI ,flag)       
             else
                flag=.false.
             endif
             if (.not.flag) then
                if (a>t(nr)) then
                   if (.not.plein) t(nr+1)=a
                else
                   allocate(tmp(n))
                   if (a<t(1)) then
                      tmp(:)=t(:)
                      t(1)=a
                      t(2:n)=tmp(1:n-1)
                   else
                      i=1
                      do while(.not.(a>t(i).and.a<t(i+1)))
                         i=i+1
                         if (i>=n) stop 'add_v: pfff ca chie'                      
                      enddo
                      tmp(:)=t(:)
                      t(i+1)=a
                      t(i+2:n)=tmp(i+1:n-1)
                   endif
                   deallocate(tmp)
                endif
             else
                nb_nomodes=nb_nomodes+1
             endif
          end if
       end if
    end if
    if (nb_nomodes/=0) print*,nb_nomodes,' faux modes pour l=',l
!
!------------------------------------------------------------------
  end subroutine add_v
!------------------------------------------------------------------

!------------------------------------------------------------------
  subroutine get_wmax_resc(wmax,wmax_resc,f_zero_tmp,n_zero,ic)
!------------------------------------------------------------------
    use def_gparam
    implicit none
    real(DP), intent(in) :: wmax
    real(DP), intent(out) :: wmax_resc
    real(DP), dimension(:) :: f_zero_tmp
    integer, intent(in)  :: n_zero,ic
    integer :: ii,ni
          
    ni=2
    ii= ic+ni
    if (ii+ni+1>n_zero) then
       wmax_resc=wmax
    else
       do while ( (f_zero_tmp(ii)<0.0_DP .or. f_zero_tmp(ii+1)<0.0_DP  &
            .or. f_zero_tmp(ii+3)<0.0_DP) &
            .and.ii<n_zero)
          ii=ii+1
       enddo
       if (ii>=n_zero) then
          wmax_resc=wmax
       else       
          wmax_resc=(f_zero_tmp(min(ii+1,n_zero)) &
               +f_zero_tmp(min(ii,n_zero)))*PI
       endif
    endif
!------------------------------------------------------------------
      end subroutine get_wmax_resc
!------------------------------------------------------------------
!------------------------------------------------------------------
  subroutine get_wmin_resc(wmin,wmin_resc,f_zero_tmp,n_zero,ic)
!------------------------------------------------------------------
    use def_gparam
    implicit none
    real(DP), intent(in) :: wmin
    real(DP), intent(out) :: wmin_resc
    real(DP), dimension(:) :: f_zero_tmp
    integer, intent(in)  :: n_zero,ic
    integer :: ii,ni
    ni=3      
    ii= ic-ni
    if (ii<1) then 
       wmin_resc=wmin
    else
       do while ( (f_zero_tmp(ii)<0.0_DP .or. f_zero_tmp(ii-1)<0.0_DP) &
            .and.ii<n_zero)
          ii=ii-1
       enddo
       if (ii<=1) then
          wmin_resc=wmin
       else       
          wmin_resc=(f_zero_tmp(min(ii-1,n_zero)) &
               +f_zero_tmp(min(ii,n_zero)))*PI
       endif
    endif
!------------------------------------------------------------------
      end subroutine get_wmin_resc
!------------------------------------------------------------------
!!$!------------------------------------------------------------------
!!$subroutine resque_table(l_in,wmin,wmax,n,table,det1,det2)
!!$!recherche de facon bestiale n freq propres entre wmin et wmax
!!$!------------------------------------------------------------------
!!$  use def_gaparam
!!$  use bits, only: l,knsw
!!$  use shanks, only: maxo
!!$!
!!$  integer, intent(in) :: l_in,n
!!$  real(DP), intent(in):: wmin,wmax,det1,det2
!!$  real(DP), dimension(n), intent(out):: table
!!$!
!!$  integer, parameter :: sens=0
!!$  real(DP), det2
!!$!  
!!$  knsw=0 !kount swich off
!!$  maxo=8 ! precision rks ?
!!$  l=l_in
!!$!
!!$  if (n>1) stop 'resque_table: je n''ai pas encore prevu n>1!'
!!$!
!!$  do while(not_found)
!!$     wt=(w1+w2)/2
!!$  enddo
!!$!  
!!$!------------------------------------------------------------------
!!$end subroutine resque_table
!!$!------------------------------------------------------------------

!------------------------------------------------------------------
subroutine compte_nb_modes(l_in,wt,k,det)
!------------------------------------------------------------------
  use bits, only: l,fl,fl1,fl2,fl3,sfl3,knsw,jcom,cond_limite
  use shanks, only:maxo,stepf
  use minos
  implicit none
!
  integer, intent(in) :: l_in
  doubleprecision, intent(in) :: wt
  integer, intent(out) :: k
  doubleprecision, optional, intent(out) :: det
!
  integer, parameter :: sens=0,inss=5
  doubleprecision :: bid,w1,w2,det1,det2
  integer :: i
!
  l=l_in  
!
  fl=l
  fl1=fl+1.d0
  fl2=fl+fl1
  fl3=fl*fl1 
  sfl3=dsqrt(fl3) 
  knsw=1 !kount swich on
  if (l>4) then
     maxo=inss
     stepf=1.d0
  else
     maxo=inss
     stepf=0.25d0
  endif
  call detqn(wt,k,bid,sens)
  select case(jcom)
  case(3)
     if (l==1) then
        if (cond_limite==1) then
           k=k+2
        else
           k=k+1
        endif
     else
        k=k+1
     endif
  case(2,1)
     k=k+1
  end select
  if (present(det)) det=bid
!------------------------------------------------------------------
end subroutine compte_nb_modes
!------------------------------------------------------------------

!------------------------------------------------------------------
subroutine get_ANNp_R(l_in,fp,A00,A11,A01,A10)
!------------------------------------------------------------------
  use def_gparam
  use bits, only:l,jcom
  implicit none
!
  integer, intent(in) :: l_in
  real(DP), intent(in) :: fp
  real(DP), intent(out) :: A00,A11,A01,A10
!
  real(DP) :: wp
  real(DP), dimension(5) :: Ann
!pour minos
  l=l_in
  if (l==0) then
     jcom=1
  else
     jcom=3
  endif
!
  wp=2.0_DP*PI*fp
  call get_cop_fctp(wp,Ann)
  A00=Ann(1)
  A01=Ann(2)
  A10=Ann(3)
  A11=Ann(4)  
!------------------------------------------------------------------
end subroutine get_ANNp_R
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine get_ANNp_L(l_in,fp,A5)
!------------------------------------------------------------------
  use def_gparam
  use bits, only:l,jcom
  implicit none
!
  integer, intent(in) :: l_in
  real(DP), intent(in) :: fp
  real(DP), intent(out) :: A5
!
  real(DP) :: wp
  real(DP), dimension(5) :: Ann
!pour minos
  l=l_in
  if (l==0) then
     A5=0.0_DP
  else
     jcom=2
!
     wp=2.0_DP*PI*fp
     call get_cop_fctp(wp,Ann)
     A5=Ann(1)
  endif
!------------------------------------------------------------------
end subroutine get_ANNp_L
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine ANNcR(l_in,fc,A00,A01,A10,A11)
!------------------------------------------------------------------
  use def_gparam
  use bits, only:  jcom, knsw,l,fl,fl1,fl2,fl3,sfl3
  use shanks, only: maxo,stepf
  use minos
  implicit none
!
  integer        , intent(in) :: l_in
  doubleprecision, intent(in) :: fc
  doubleprecision, intent(out):: A00,A01,A10,A11
!
  doubleprecision :: wp,bid
  integer         :: ibid
  integer, parameter :: continu=0,sens=0
  doubleprecision, dimension(5) :: ann
!
!
  if (fc<1.E-30_DP) stop 'ANNcR: freq d''entree nulle!'
  l=l_in  
!
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
  stepf=1.d0
  knsw=0 !kount swich off
  maxo=8 ! precision rks ?
!
  wp=2._DP*fc*PI
  call detqn(wp,ibid,bid,sens,Ann,continu)
!
  A00=Ann(1)
  A01=Ann(2)
  A10=Ann(3)
  A11=Ann(4)
!------------------------------------------------------------------
end subroutine ANNcR
!------------------------------------------------------------------

!------------------------------------------------------------------
subroutine get_det(l_in,fc,det)
!------------------------------------------------------------------
  use def_gparam
  use bits, only:  jcom, knsw,l,fl,fl1,fl2,fl3,sfl3
  use shanks, only: maxo,stepf
  use minos
  implicit none
!
  integer        , intent(in) :: l_in
  doubleprecision, intent(in) :: fc
  doubleprecision, intent(out):: det
!
  doubleprecision :: wp,bid
  integer         :: ibid
  integer, parameter :: sens=0
!
!
  if (fc<1.E-30_DP) stop 'ANNcR: freq d''entree nulle!'
  l=l_in  
!
  if (jcom/=2) then
     if (l==0) then
        jcom=1
     else
        jcom=3
     endif
  endif
!
  fl=l
  fl1=fl+1.d0
  fl2=fl+fl1
  fl3=fl*fl1 
  sfl3=dsqrt(fl3) 
  stepf=1.d0
  knsw=0 !kount swich off
  maxo=8 ! precision rks ?
!
  wp=2._DP*fc*PI
  call detqn(wp,ibid,det,sens)
!
!------------------------------------------------------------------
end subroutine get_det
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine ANNcR0(l,Ann,Annp)
!------------------------------------------------------------------
  use def_gparam
  integer, intent(in) :: l
  real(DP), dimension(4), intent(out) :: Ann,Annp
!
  real(DP) :: A00,A01,A10,A11,a,b,c,d,f1,f2
!
  f1=1.E-6_DP
  f2=2.E-6_DP
  df=f2-f1
!
  call ANNcR(l,f1,A00,A01,A10,A11)
  call ANNcR(l,f2,a  ,b  ,c  ,d  )
!  
  Ann (1)=A00
  Annp(1)=(a-A00)/df
  Ann (2)=A01
  Annp(2)=(b-A01)/df
  Ann (3)=A10
  Annp(3)=(c-A10)/df
  Ann (4)=A11
  Annp(4)=(d-A11)/df
  
!------------------------------------------------------------------
end subroutine ANNcR0
!------------------------------------------------------------------


!------------------------------------------------------------------
subroutine ANNcL0(l,A5,A5p)
!------------------------------------------------------------------
  use def_gparam
  integer, intent(in) :: l
  real(DP),intent(out) :: A5,A5p
!
  real(DP) :: a,f1,f2
!
  f1=1.E-10_DP
  f2=2.E-10_DP
  df=f2-f1
!
  call ANNcL(l,f1,A5)
  call ANNcL(l,f2,a)
!  
  A5p=(a-A5)/df  
!------------------------------------------------------------------
end subroutine ANNcL0
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine ANNcL(l_in,fc,A5)
!------------------------------------------------------------------
  use def_gparam
  use bits, only:  jcom, knsw,l,fl,fl1,fl2,fl3,sfl3
  use shanks, only: maxo,stepf
  use minos
  implicit none
!
  integer        , intent(in) :: l_in
  doubleprecision, intent(in) :: fc
  doubleprecision, intent(out):: A5
!
  doubleprecision :: wp,bid
  integer         :: ibid
  integer, parameter :: continu=0,sens=0
  doubleprecision, dimension(5) :: ann
!
!
  l=l_in  
!
  if (l==0) then
     A5=0.0d0
  else
     jcom=2
!
     fl=l
     fl1=fl+1.d0
     fl2=fl+fl1
     fl3=fl*fl1 
     sfl3=dsqrt(fl3) 
     stepf=1.d0
     knsw=0 !kount swich off
     maxo=8 ! precision rks ?
!
     wp=2._DP*fc*PI
     call detqn(wp,ibid,bid,sens,Ann,continu)
!
     A5=Ann(1)
  endif
!------------------------------------------------------------------
end subroutine ANNcL
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine get_cop_fctp(wp,Ann,fctp)
!------------------------------------------------------------------
  use modes
  use bits, only:  jcom, knsw,sfl3,wn,rn,cg,qinv
  use rindex, only :n
  use eifx, only: ar
  use param_modele, only: ray
  use shanks, only: maxo
  use minos
  implicit none
!
  doubleprecision, intent(in) :: wp
  doubleprecision, dimension(:) :: Ann
  type(mode), optional, intent(inout) :: fctp
!
  integer :: ibid,i,j,sens
  doubleprecision :: bid,scale
  doubleprecision, parameter :: rhobar=5515.d0    
!
  knsw=0 !kount swich off
  maxo=8 ! precision rks ?
!
! on va utiliser detqn pour retrouver les modes a partir des mineurs
! et en integrant de la surface vers le centre
  sens=1
!utilisation de detqn pour calculer les fctps
  call detqn(wp,ibid,bid,sens,Ann)
!
  if (present(fctp)) then
     call allocate_mode(fctp,n,jcom)
     scale=1.d0/dsqrt(rhobar*rn**3)
     fctp%n    =n
     fctp%w    =wp
     fctp%q    =1.d0/qinv
     fctp%g    =cg
     select case(jcom)
     case(3)
!on normalise: 
        ar(1:2,:)=ar(1:2,:)*wp/wn*scale
        ar(3:4,:)=ar(3:4,:)*wp/wn/sfl3*scale
!
        fctp%u (:,1)=ar(1,:)
        fctp%up(:,1)=ar(2,:)
        fctp%u (:,2)=ar(3,:)
        fctp%up(:,2)=ar(4,:)
     case(2)
        ar(1:2,:)=ar(1:2,:)*wp/wn/sfl3*scale
!
        fctp%u (:,1)=ar(1,:)
        fctp%up(:,1)=ar(2,:)
     case(1)
        ar(1:2,:)=ar(1:2,:)*wp/wn*scale
!
        fctp%u (:,1)=ar(1,:)
        fctp%up(:,1)=ar(2,:)
     end select
  endif

   
!------------------------------------------------------------------ 
end subroutine get_cop_fctp
!------------------------------------------------------------------  

!------------------------------------------------------------------
subroutine get_det_curve(l_in,wmi,wma)
!------------------------------------------------------------------
  use bits, only: l,fl,fl1,fl2,fl3,sfl3,knsw,jcom
  use shanks, only:maxo,stepf
  use minos
  implicit none
!
  integer, intent(in) :: l_in
  doubleprecision, intent(in) :: wmi,wma
!
  integer, parameter :: sens=0,inss=5
  doubleprecision :: w,dw,bid
  integer :: i,k
!
  print*,'Using get_det_curve for check for l=',l
  l=l_in  
!
  fl=l
  fl1=fl+1.d0
  fl2=fl+fl1
  fl3=fl*fl1 
  sfl3=dsqrt(fl3) 
  knsw=1 !kount swich on
!  if (l>4) then
     maxo=inss
     stepf=1.d0
!  else
!     maxo=inss
!     stepf=0.25d0
!  endif
!  dw=(wma-wmi)/1000.d0
!  do i=1,1000
  dw=(wma-wmi)/100.d0
  do i=1,100
     w=i*dw+wmi
     print*,i,w
     call detqn(w,k,bid,sens)
     write(111,*)w,bid
     write(1111,*)w,k
  enddo
!------------------------------------------------------------------
end subroutine get_det_curve
!------------------------------------------------------------------
!--------------------------------------------------------------------------
end module neominos
!--------------------------------------------------------------------------
