!--------------------------------------------------------------------------------
module modes
!--------------------------------------------------------------------------------
  use def_gparam
  implicit none
  public :: init_modes,get_and_sum
  private
  real(SP), dimension(:,:,:), allocatable :: source_g
  complex(DP), dimension(:,:,:), allocatable ::source_gw
  real(DP) :: fmax,fmaxt,t0_delay,wmax1,wmax2,duree,deltaf,dtf,dtf2
  integer :: logunit,coef_nbt,NBF,LmaxS,LmaxT,NBF2
  integer :: Lmax_all,Nmax_allS,Nmax_allT,Nmax_all
  real(DP), parameter :: coef_damp1=1.0_DP, coef_damp2=1.2_DP
  real(DP), dimension(:,:), allocatable  :: src_coord,Msrc
  real(DP), dimension(:,:), allocatable  :: rec_coord
  real(DP), dimension(:), allocatable  :: weight
  real(DP), dimension(1) :: rho,Vs,Vp,Qshear,radtmp
  integer, dimension(:)  ,allocatable :: LLLmax 
  integer, dimension(:,:),allocatable :: NNNmax 
  integer, dimension(:), allocatable :: NmaxS,NmaxT
!freq propre et attenuation
  real(DP), dimension(:,:,:), allocatable :: eigw,qw
!
  real(DP), dimension(  :,:), allocatable :: UU,UUp,VV,VVp,WW,WWp
  real(DP), dimension(:,:,:), allocatable :: UUs,UUps,VVs,VVps,WWs,WWps
  real(DP), dimension(:) , allocatable :: al,b0,b1,b2
!
  integer :: length_qnl
  complex*16, dimension(:,:), allocatable :: time_serie
  integer, dimension(2), parameter :: ldeb=(/0,1/)
  integer, parameter :: DEP=1,VEL=2,ACC=3
  integer  :: chanel=VEL
!
  complex(DP), parameter :: czero=(0._DP,0._DP),II=(0._DP,1._DP)
  complex(SP), parameter :: czero_sp=(0._SP,0._SP)
!
  complex(DP), dimension(:,:,:,:,:  ), allocatable :: SkN
  complex(DP), dimension(:,:,:,:,:  ), allocatable :: RkN
!
  real(DP), parameter :: rhobar=5515._DP,rn=6371000._DP 
  real(DP), parameter :: scale2=1.D0/(rhobar*rn*rn*rn),scale1=scale2/rn/rn
  real(DP)            :: scale,dscale,receiver_rad
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------------------------------------------------
  subroutine init_modes(logunit_)
!-------------------------------------------------
    use def_gparam
    use global_main_param, only:source_type,NBT,dt,f0,f1,f2,f3,f4,t0,get_source_dat &
                               ,get_receivers,eigenfileS,eigenfileT,RA  &
                               ,NBE,NBR,nbcomp,nmax_in
    use time_function 
    use util_fctp
!
    implicit none
    integer, intent(in) :: logunit_
!
    real(DP), dimension(:,:), allocatable :: xdum,tmp
    real(DP), dimension(:)  , allocatable :: rad
    integer :: NBFtmp,l,n,ir,ic,unit,i
    real(DP):: t0delay
    character(len=30) :: name
    logical :: flag
    real :: bid
!
    t0delay=0.d0
    unit=121
    logunit=logunit_
    chanel=VEL
    allocate(source_g(NBT,1,1))
    t0_delay=0.0d0
    call def_timefunc(source_type,NBT,dt,f0,f1,f2,f3,f4,t0,source_g(:,1,1),fmax)  
    fmaxt=coef_damp1*fmax
    fmax =coef_damp2*fmax
    wmax1=fmaxt*2._DP*PI
    wmax2=fmax *2._DP*PI
    duree=NBT*dt+t0delay
!pas en frequence pour la FFT
    deltaf=1./duree
!nombre de frequence  , on va a 1.1 fois Nyquist
!mais aussi nombre de pas de temps effectif
    NBFtmp   = int(2*fmax*((NBT-1)*dt+t0delay))+1
    call get_power(int(NBFtmp*1.1),NBF)
    NBF2=2*NBF
!on utilise deux fois plus long pour la convolution (zero padding sur la mointie
!du signal)
!on recupere le plus grand coef_nbt tel que NBF*coef_nbt < NBT
    coef_nbt=1
    do while (coef_nbt*2*NBF < NBT)
       coef_nbt=2*coef_nbt
    enddo
    dtf   = ((NBT-1)*dt)/real(NBF-1)
!pas de temps apres coef_nbt
    dtf2  = dtf/coef_nbt
!lecture des sources:
    allocate(src_coord(3,NBE),Msrc(6,NBE),rec_coord(2,NBR),weight(NBE))
    call get_source_dat(src_coord,Msrc)
    call get_receivers(rec_coord)
    write(logunit,*) 'There are ',NBE,' sources'
    write(logunit,*) 'There are ',NBR,' receivers'
!
    print*,'Getting sources and receivers informations ...'
!
!ouverture de fichiers de fonctions propres
    print*,'Opening eigenfunction files ...'
    call open_fctp(eigenfileS,eigenfileT,rad_recepteur=receiver_rad)    
!calcule des LmaxS, LmaxT, NmaxS,NmaxT
    call get_Lmax(fmax,LmaxS,LmaxT)
    print*,'LmaxS, LmaxT:',LmaxS, LmaxT
    write(logunit,*) 'LmaxS, LmaxT:',LmaxS, LmaxT
    allocate(NmaxS(0:LmaxS),NmaxT(0:LmaxT))
    call get_Nmax(fmax,NmaxS,NmaxT)
    Nmax_all=max(maxval(NmaxS(:)),maxval(NmaxT(:)))
    Lmax_all=max(LmaxS,LmaxT)
    Nmax_allS=maxval(NmaxS(:))
    Nmax_allT=maxval(NmaxT(:))
    allocate(LLLmax(2))
    allocate(NNNmax(0:lmax_all,2))
    NNNmax(:,:)=0
    LLLmax(1)        =LmaxS
    LLLmax(2)        =LmaxT
    NNNmax(0:LmaxS,1)=NmaxS(:)
    NNNmax(0:LmaxT,2)=NmaxT(:)
    if (nmax_in/=-1) NNNmax(:,:)=min(nmax_in,NNNmax(:,:))
! 
!init b0 etc ...
!
    allocate(al(0:Lmax_all),b0(0:Lmax_all),b1(0:Lmax_all),b2(0:Lmax_all))
    do l=0,lmax_all
       al(l)   =dble(l*(l+1))
       b0(l)   =dsqrt(dble(2*l+1)/(4.0d0*PI))
       b1(l)   =dsqrt(al(l))*b0(l)/2.0d0
       b2(l)   =dsqrt(dble((dble(l)-1)*dble(l)*(dble(l)+1) &
            *(dble(l)+2)))*b0(l)/4.0d0
    enddo
    allocate(eigw(0:Nmax_all,0:Lmax_all,2),qw(0:Nmax_all,0:Lmax_all,2))
    call get_eigenw(eigw,qw,LmaxS,NmaxS,LmaxT,NmaxT)
    print*,'Retreving eigenfrequencies for ',1+NBE,' sources and receiver depth'
    print*,'allocating eigenfunctions (',            &
         ((Nmax_all+1)*(LmaxS+1)*4+         &
         (maxval(NmaxT(:))+1)*(LmaxT+1)*2)*8/1024./1024. &
         ,' Mb/proc)'
    allocate(UU(0:Nmax_allS,0:LmaxS),UUp(0:Nmax_allS,0:LmaxS) &
         ,VV(0:Nmax_allS,0:LmaxS),VVp(0:Nmax_allS,0:LmaxS) )
    print*,'allocating temporary eigenfunctions (',            &
         real(NBE)*((Nmax_all+1)*(LmaxS+1)*4+         &
         (maxval(NmaxT(:))+1)*(LmaxT+1)*2)*8/1024./1024. &
         ,' Mb/proc)'
    allocate(UUs(NBE,0:Nmax_allS,0:LmaxS),UUps(NBE,0:Nmax_allS,0:LmaxS) &
         ,VVs(NBE,0:Nmax_allS,0:LmaxS),VVps(NBE,0:Nmax_allS,0:LmaxS) )
    allocate(WW(0:Nmax_allT,0:LmaxT),WWp(0:Nmax_allT,0:LmaxT))
    allocate(WWs(NBE,0:Nmax_allT,0:LmaxT),WWps(NBE,0:Nmax_allT,0:LmaxT))
    UU=0.0_DP;VV=0.0_DP;UUp=0.0_DP;VVp=0.0_DP;WW=0.0_DP;WWp=0.0_DP
!
    print*,'reading and interpolating eigenfunction ...'
!
    allocate(xdum(4,NBE+1),rad(NBE+1))
!
    print*,'receiver_rad=',receiver_rad
    rad(1)=receiver_rad
    rad(2:NBE+1)=src_coord(1,:)
    !

    print*,'->for S modes ...'
    do l=0,LmaxS
       do n=0,NmaxS(l)
          call read_interpS(n,l,1+NBE,rad,xdum)
          UU( n,l) =xdum(1,1)
          UUp(n,l) =xdum(2,1)
          VV( n,l) =xdum(3,1)
          VVp(n,l) =xdum(4,1)
!
          UUs( :,n,l)=xdum(1,2:)
          UUps(:,n,l)=xdum(2,2:)
          VVs( :,n,l)=xdum(3,2:)
          VVps(:,n,l)=xdum(4,2:)
          
       enddo
    enddo
    deallocate(xdum)
    print*,'->for T modes ...'
    allocate(xdum(2,NBE+1))
    do l=0,LmaxT
       do n=0,NmaxT(l)
          call read_interpT(n,l,1+NBE,rad,xdum)
          WW( n,l) =xdum(1,1)
          WWp(n,l) =xdum(2,1)
!  
          WWs( :,n,l)=xdum(1,2:)
          WWps(:,n,l)=xdum(2,2:)
       enddo
    enddo
    deallocate(xdum,rad)    
!
    call init_time_serie()
!
    print*,'init_modes done!'
!-------------------------------------------------
  end subroutine init_modes
!-------------------------------------------------

!-------------------------------------------------
  subroutine get_and_sum
!-------------------------------------------------
    use def_gparam
    use global_main_param, only: NBR,nbcomp,NBT,dt,comp,compr,stations_name,rotation,NBE,stacking_sources
    implicit none
    integer :: ir,ic,i,itp,itt,is
    complex(DP), dimension(length_qnl,nbcomp) :: RSkN
    complex(DP), dimension(NBF2,nbcomp,NBR,NBE)    :: traces
    real   (SP), dimension(NBT,nbcomp) :: trace_time
    real   (SP), dimension(NBT) :: tmptr
    complex*16 :: alpha=(1.d0,0.d0),beta=(0.d0,0.d0)
    character(len=30) :: name
    real(DP) :: t,p
!
    print*,'Computing RkN and SkN terms ...'
    traces(:,:,:,:)=(0.0d0,0.0d0)
    call init_RSkN()
    do is=1,NBE
       do ir=1,NBR
          RSkN(:,:)=czero
          call get_RSkN(ir,is,RSkN)
          do ic=1,nbcomp
             call ZGEMV('T',length_qnl,NBF2,alpha,time_serie ,length_qnl  &
                  ,RSkN(1,ic),1,beta,traces(1,ic,ir,is),1)  
          enddo
       enddo
    enddo
    print*,'converting in time and writing ascii output ...'
    if (NBE==1.or.stacking_sources) then
       do ir=1,NBR
          if (stacking_sources) then
             trace_time(:,:)=0.
             do ic=1,nbcomp
                do i=1,NBE
                   call convert_f2t(traces(:,ic,ir,is),tmptr(:),is)
                   trace_time(:,ic)=trace_time(:,ic)+tmptr(:)
                enddo
             enddo
          else
             do ic=1,nbcomp
                call convert_f2t(traces(:,ic,ir,1),trace_time(:,ic),1)
             enddo
          endif
          call apply_source_filter(trace_time)
          if (rotation) call rotate(trace_time,ir,1)
          do ic=1,nbcomp
             if (rotation) then
                write(name,'("U",a1,"_",a4)') compr(ic),stations_name(ir)
             else
                write(name,'("U",a1,"_",a4)') comp(ic) ,stations_name(ir)
             endif
             open(17,file=name)
             write(17,'(f8.1,1x,e13.7)') ((i-1)*dt,trace_time(i,ic),i=1,NBT)
             close(17)
          enddo
       enddo
    else
       do is=1,NBE
          do ir=1,NBR
             do ic=1,nbcomp
                call convert_f2t(traces(:,ic,ir,is),trace_time(:,ic),is)
             enddo
             call apply_source_filter(trace_time)
             if (rotation) call rotate(trace_time,ir,is)
             do ic=1,nbcomp
                if (rotation) then
                   write(name,'("U",a1,"_",a4,"_",i3.3)') compr(ic),stations_name(ir),is

                else
                   write(name,'("U",a1,"_",a4,"_",i3.3)') comp(ic) ,stations_name(ir),is
                endif
                open(17,file=name)
                write(17,'(f8.1,1x,e13.7)') ((i-1)*dt,trace_time(i,ic),i=1,NBT)
                close(17)
             enddo
          enddo
       enddo
    endif
!-------------------------------------------------
  end subroutine get_and_sum
!-------------------------------------------------
!-------------------------------------------------
  subroutine rotate(trace,ir,is)
!-------------------------------------------------
    use def_gparam
    use util_fctp
    use global_main_param, only: NBT,stations_name
    implicit none
    real(SP), dimension(:,:), intent(out) :: trace
    integer, intent(in) :: ir,is
!
    real(DP) :: tr,pr,ts,ps,asr,bsr,gsr,bazm
    real(SP), dimension(NBT,2) :: tmp
!    

    tr=rec_coord(1,ir)
    pr=rec_coord(2,ir)
    ts=src_coord(2,is)
    ps=src_coord(3,is)
    call euler(ts,ps,tr,pr,asr,bsr,gsr)
    bazm=-asr
    print*,stations_name(ir),', bazm= ',bazm*180./3.1415926,', delta=',bsr*180./3.1415926
    tmp(:,:)=trace(:,2:3)
    trace(:,2)=-tmp(:,1)*cos(bazm)+tmp(:,2)*sin(bazm)
    trace(:,3)=-tmp(:,1)*sin(bazm)-tmp(:,2)*cos(bazm)
!-------------------------------------------------
  end subroutine rotate
!-------------------------------------------------
!-------------------------------------------------
  subroutine convert_f2t(tracef,tracet,is)
!-------------------------------------------------
    use def_gparam
    use util_fctp
    use global_main_param, only: NBT,dt,nbcomp ,t_delay_sources
    implicit none
    complex(DP), dimension(NBF2), intent(in) :: tracef
    real, dimension(NBT), intent(out) :: tracet
    integer, intent(in) :: is
!
!
    complex(DP), dimension(NBT) :: trace_test
    real(DP), dimension(NBF*coef_nbt) :: trace_testr
    real, dimension(NBT,1,1) :: trace_r
    complex(DP), dimension(NBF2*coef_nbt) :: traces_testf
    integer :: ic,ir,i,num
!      
    traces_testf(:)=(0.d0,0.d0)
    traces_testf(1:NBF+1)=tracef(1:NBF+1)
    traces_testf((2*coef_nbt*NBF)-NBF+1:coef_nbt*NBF2) &
               =tracef(NBF+1:NBF2)
    call dfour1(traces_testf,NBF2*coef_nbt,-1)
!SP + * frequency step
    trace_testr(1:NBF*coef_nbt)=real(traces_testf(1:NBF*coef_nbt)/NBF2/dtf,DP)
    call resample(trace_testr,t_delay_sources(is),dtf2,NBF*coef_nbt &
         ,tracet,dt,NBT)
!-------------------------------------------------
  end subroutine convert_f2t
!-------------------------------------------------

!-------------------------------------------------
  subroutine init_RSkN
!-------------------------------------------------
    use def_gparam
    use global_main_param, only: nbcomp,NBE
    implicit none
    integer :: l
    real :: mega
!
    scale =1._DP/sqrt(rhobar*rn**3)
    dscale=1._DP/sqrt(rhobar*rn**3)/rn
!
    mega=(3*nbcomp*(real(Nmax_allS+1)*(LmaxS+1)+real((Nmax_allT+1)*(LmaxT+1))) + &
          5*NBE*(real(Nmax_allS+1)*(LmaxS+1)+real((Nmax_allT+1)*(LmaxT+1))) )/1024./1204.
    allocate(RkN(-1:1,nbcomp,0:Nmax_all,0:Lmax_all,2), &
         SkN(-2:2       ,0:Nmax_all,0:Lmax_all,2,NBE))
    RkN=0.0_DP
    SkN=0.0_DP
    print*,'->SkN ...'
    call SkNfun(SkN)
    print*,'->RkN ...'
    call RkNfun(RkN)
! Cleaning unecessary arrays:
    deallocate(UUs,UUps,VVs,VVps,WWs,WWps,Msrc)   
!-------------------------------------------------
  end subroutine init_RSkN
!-------------------------------------------------
!-------------------------------------------------
  subroutine SkNfun(S)
!-------------------------------------------------
    use def_gparam
    use global_main_param, only:RA,NBE
    IMPLICIT NONE 
    complex(DP), dimension(-2:,0:,0:,:,:), intent(out) :: S
    real(DP) :: ff,xx,vor,dmtp,egu,egv,degu,degv,rr
    real(DP) :: zz,wor,fw,dfw
    real(DP), dimension(6) :: am
    integer :: is,n,l
!
    do l=0,LmaxS
       do n=0,NmaxS(l)
          do is=1,NBE
             rr=src_coord(1,is)/RA
             egu  =UUs (is,n,l)
             egv  =VVs (is,n,l)
             degu =UUps(is,n,l)
             degv =VVps(is,n,l) 
             am(:)=Msrc(:,is)
!    
             ff  =(2.0d0*egu-al(l)*egv)/rr
             xx  =-b1(l)*(degv+(egu-egv)/rr)
             vor =b2(l)*egv/rr
             dmtp=am(3)-am(2)
             S( 0,n,l,1,is)=-b0(l)*(degu*am(1)+ff*0.5*(am(2)+am(3)))*dscale
             S( 1,n,l,1,is)=xx*dcmplx(-am(4),am(5))*dscale
             S(-1,n,l,1,is)=xx*dcmplx( am(4),am(5))*dscale
             S( 2,n,l,1,is)=vor*dcmplx(dmtp, 2.0d0*am(6))*dscale
             S(-2,n,l,1,is)=vor*dcmplx(dmtp,-2.0d0*am(6))*dscale
          enddo
       enddo
    enddo
    do l=1,LmaxT
       do n=0,NmaxT(l)
          do is=1,NBE
             rr=src_coord(1,is)/RA
             fw    =WWs (is,n,l)
             dfw   =WWps(is,n,l)
             am(:) =Msrc(:,is)
!
             zz    =-b1(l)*(dfw-fw/rr)
             wor   = b2(l)*fw/rr
             dmtp  =am(2)-am(3)
             S( 1,n,l,2,is)=zz *dcmplx(  am(5),am(4))*dscale
             S(-1,n,l,2,is)=zz *dcmplx( -am(5),am(4))*dscale
             S( 2,n,l,2,is)=wor*dcmplx(2.0d0*am(6), dmtp)*dscale
             S(-2,n,l,2,is)=wor*dcmplx(2.0d0*am(6),-dmtp)*dscale
          enddo
       enddo
    enddo
!--------------------------------------------------
  end subroutine SkNfun
!--------------------------------------------------

!--------------------------------------------------
  subroutine RkNfun(R)
!--------------------------------------------------
    use def_gparam
    use global_main_param, only: NBE,nbcomp
    IMPLICIT NONE
    complex(DP), dimension(-1:,:,0:,0:,:), intent(out) :: R
    real(DP) :: xu,xv,u,v,vir,vip,vit,x
    integer :: ic,n,l
!
    do l=0,LmaxS
       do n=0,NmaxS(l)
          do ic=1,nbcomp !boucle sur les composantes
             select case(ic)
             case(1) ! verticale
                vir=1.0_DP
                vit=0.0_DP
                vip=0.0_DP
             case(2) !Nord=-theta
                vir=0.0_DP
                vit=-1.0_DP
                vip=0.0_DP
             case(3) !Est=phi
                vir=0.0_DP
                vit=0.0_DP
                vip=1.0_DP
             end select
             u =UU (n,l)
             v =VV (n,l)    
             xu= b0(l)*u
             xv=-b1(l)*v
             R( 0,ic,n,l,1)=xu*vir*scale
             R( 1,ic,n,l,1)=xv*dcmplx( vit,vip)*scale
             R(-1,ic,n,l,1)=xv*dcmplx(-vit,vip)*scale
          enddo
       enddo
    enddo
    do l=1,LmaxT
       do n=0,NmaxT(l)
          do ic=1,nbcomp !boucle sur les composantes
             select case(ic)
             case(1) ! verticale
                vir=1.0_DP
                vit=0.0_DP
                vip=0.0_DP
             case(2) !theta
                vir=0.0_DP
                vit=-1.0_DP
                vip=0.0_DP
             case(3) !phi
                vir=0.0_DP
                vit=0.0_DP
                vip=1.0_DP
             end select
             x=-b1(l)*WW (n,l)             
!
             R ( 1,ic,n,l,2)=x*dcmplx(-vip,vit)*scale
             R (-1,ic,n,l,2)=x*dcmplx( vip,vit)*scale
          enddo
       enddo
    enddo
!--------------------------------------------------
  end subroutine RkNfun
!--------------------------------------------------
!-------------------------------------------------
  subroutine get_RSkN(ir,is,RSkN)
!-------------------------------------------------
    use def_gparam 
    use util_fctp
    use global_main_param, only: nbcomp,source_delay,t_delay_sources,nmin
    implicit none
    integer, intent(in) :: ir,is
    complex(DP), dimension(:,:), intent(out) :: RSkN
!
    integer     :: iqnl,bNp,qdeb,qfin,l,n,q,bN,ird,ic
    complex(DP) :: sd,wq,eiwtd,rs
    real(DP) :: ts,ps,tr,pr,asr,bsr,gsr,x
    complex(DP), dimension(-2:2) :: csd
    complex(DP), dimension(-1:1) :: csr
    real(DP), dimension(-2:2,-2:2,0:Lmax_all) :: dsr
!
    tr=rec_coord(1,ir)
    pr=rec_coord(2,ir)
    ts=src_coord(2,is)
    ps=src_coord(3,is)
    call euler(ts,ps,tr,pr,asr,bsr,gsr)
    do bNp=-2,2
       x=bNp*gsr
       csd(bNp)=dcmplx(dcos(x),dsin(x))
    enddo
    do bNp=-1,1
       x=bNp*asr
       csr(bNp)=dcmplx(dcos(x),dsin(x))
    enddo
    call legendre(bsr,dsr,Lmax_all)   
    do ic=1,nbcomp
       iqnl=0
       qdeb=1;qfin=2
!
       do q=qdeb,qfin ! 1:Spheroidal,2 toroidal
          do l=ldeb(q),LLLmax(q)
             do n=nmin,NNNmax(l,q)
                if (mode_existe(n,l,q)) then
                   iqnl=iqnl+1
                   rs=czero
                   do bN=-1,1
                      sd = czero
                      do bNp=-2,2
                         sd=sd+dsr(bN,bNp,l)*SkN(bNp,n,l,q,is)*csd(bNp)
                      enddo
                      rs=rs+sd*RkN(bN,ic,n,l,q)*csr(bN)
                   enddo
                   RSkN(iqnl,ic)=rs
                endif
             enddo
          enddo
       enddo
    enddo
!-------------------------------------------------
  end subroutine get_RSkN
!-------------------------------------------------

!-------------------------------------------------
  subroutine init_time_serie()
!-------------------------------------------------        
    use  util_fctp
    use global_main_param, only: NBT,dt,nmin
    implicit none
    integer :: q,l,n,it,iqnl,qdeb,qfin,i
    real(DP):: t,wt,tmax,wt1,tshift
    complex(DP) :: eiomt,om,cc
    complex(DP),dimension(NBF2) :: serie
    tmax=(NBF+1)*dtf
!
    length_qnl=0  
    do q=1,2
       do l=ldeb(q),LLLmax(q)
          do n=nmin,NNNmax(l,q)
            length_qnl= length_qnl+1
          enddo
       enddo
    enddo
    print*,'Allocating time serie ' &
       ,real(length_qnl)*real(NBF2)*16./1024./1024.,'Mbytes'
    allocate (time_serie(length_qnl,NBF2))
    time_serie(:,:)=(0.0d0,0.0d0)
    serie(:)=(0.0d0,0.0d0)
    iqnl=0
    qdeb=1;qfin=2
!
    do q=qdeb,qfin
       do l=ldeb(q),LLLmax(q)
          do n=nmin,NNNmax(l,q)
             if (mode_existe(n,l,q)) then
                iqnl=iqnl+1
                serie(:)=(0.d0,0.d0)
                om =cmplx(eigw(n,l,q),qw(n,l,q))
                call wtcoef(eigw(n,l,q),0._DP,0._DP,wmax1,wmax2,wt1)
                do it=1,NBF
                   t=(it-1)*dtf
                   select case(chanel)
                   case(DEP)
                      cc=1._DP/om**2
                   case(VEL)
                      cc=II   /om
                   case(ACC)
                      cc=-1._DP
                   end select
                   cc=cc*wt1
                   eiomt=cc*exp(II*om*t)
                   serie(it)=eiomt
                enddo
                call dfour1(serie,NBF2,1)
                time_serie(iqnl,:)=serie(:)*dtf
             endif             
          enddo
       enddo
    enddo
!-------------------------------------------------
  end subroutine init_time_serie
!-------------------------------------------------

!----------------------------------------------------------------------
    subroutine get_power(n_in,n_out)
! subroutine qui retourne le nombre entier superieur le plus proche de n_in
! qui soit une decomposable en puissance de
! 2 de 3 et de 5
!----------------------------------------------------------------------
      implicit none
      integer, intent(in) :: n_in
      integer, intent(out)::n_out
!
      integer:: n_div,n_tmp,i
!
      n_tmp=n_in-1
      n_div=n_in-1
      do while(n_div/=1) ! on veux un nb divisible par 4
         n_tmp=n_tmp+1
         n_div=n_tmp
         do i=2,5
            if (i/=4) then
               do while (mod(n_div,i)==0.and.(n_div/=1))
                  n_div=n_div/i
               enddo
            endif
         enddo
      enddo
      n_out=n_tmp
!----------------------------------------------------------------------
end subroutine get_power
!----------------------------------------------------------------------
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
    if (f.le.f3.and.f.ge.f2) then
       wt=1.0_DP
    else if (f.gt.f4.or.f.lt.f1 ) then
       wt=0.0_DP
    else if (f.gt.f3.and.f.le.f4) then
       wt=0.5_DP*(1.0+cos(pi*(f-f3)/(f4-f3)))
    else if (f.ge.f1.and.f.lt.f2) then
       wt=0.5_DP*(1.0+cos(pi*(f-f2)/(f2-f1)))
    endif
!----------------------------------------------------------------------
  end subroutine wtcoef
!----------------------------------------------------------------------
!**********************************************************************
      subroutine legendre(thetad,PlmN,ll)
!**********************************************************************
!    Calcule la matrice TLN(theta,phi) de rotation sur la sphere 
!    pour le calcul a la source
!***
!*** la base choisie est celle des fonctions beta+ et eta+ associee
!***  aux modes reels                  ( Lognonne et al)
!***
!   les polynomes de Legendre sont normalises (i.e. * par (2l+1)/2 )
! Attention, version modifiee par yann capdeville en fervier 96
! puis remodifie en 2001
!***********************************************************************
      implicit none
      doubleprecision, intent(in) :: thetad
      integer        , intent(in) :: ll
      doubleprecision, dimension(-2:2,-2:2,0:ll), intent(out) :: PlmN
!
      integer, parameter :: Nmax=2,Mmax=3
      integer :: L,m,Mm,n,N1,Mbeg,in,im
      doubleprecision :: theta,tsq,gaus,aa,tsqi,a,b,c,d,e,tsq1 &
                        ,tsq2,val,coef
      doubleprecision, dimension(0:ll,0:ll+1,0:Nmax)  ::  &
                      Ppleg,Pmleg,Qpleg,Qmleg
!
      theta=thetad
!***
      gaus=dcos(theta)
      tsq=dsqrt(1.d0-gaus*gaus)
      if(tsq.le.1.d-10) then
!boucle sur L
         do  L=0,ll
            aa=dsqrt(dfloat(L)+.5d0)
            Mm=min(L,Mmax)
            do  m=0,Mm
               Ppleg(L,m,0)=0.d0
               Ppleg(L,m,1)=0.d0
               Pmleg(L,m,1)=0.d0
               Ppleg(L,m,2)=0.d0
               Pmleg(L,m,2)=0.d0
            enddo
            if (gaus.gt.0.0d0) then
               Ppleg(L,0,0)=aa
               Ppleg(L,1,1)=aa
               Ppleg(L,2,2)=aa
            else
               Ppleg(L,0,0)=(-1)**L*aa
               Ppleg(L,1,1)=0.0d0
               Ppleg(L,2,2)=0.0d0 
            endif
!fin boucle sur L
         enddo
         Ppleg(0,1,1)=0.d0
         Ppleg(0,2,2)=0.d0
         Ppleg(1,2,2)=0.d0
      else
         tsqi=1.d0/tsq
!------------------> Calcul des polynomes de Legendre N=0
         Ppleg(0,0,0)=dsqrt(.5d0)
         Ppleg(1,0,0)=gaus*dsqrt(1.5d0)
         Ppleg(1,1,0)=-tsq*dsqrt(.75d0)
         Ppleg(2,1,0)=-3.d0*tsq*gaus*dsqrt(5.d0/12.d0)
         Ppleg(2,2,0)=3.d0*tsq**2*dsqrt(5.d0/48.d0)
! Calcul des polynomes m=0,m=1
      do  m=0,1
           do l=m+2,ll
           a=dsqrt((dble(l*l)-.25d0)/dble(l*l-m*m))
           b=dsqrt(dble((2*l+1)*(l-m-1)*(l+m-1))/dble((2*l-3) &
                 *dble(l-m)*dble(l+m)))
           Ppleg(l,m,0)=2.d0*a*gaus*Ppleg(l-1,m,0)-b*Ppleg(l-2,m,0)
        enddo
     enddo
! Calcul des polynomes de Legendre pour m>=0
      do  m=2,Mmax
           do  l=m,ll
              c=dsqrt(dble((2*l+1)*(m+l-1)*(l+m-3))/ &
                   dble((2*l-3)*(l+m)*(m+l-2)))
              d=dsqrt(dble((2*l+1)*(l+m-1)*(l-m+1))/ &
                   dble((2*l-1)*(l+m)*(l+m-2)))
              e=dsqrt(dble((2*l+1)*(l-m))/           &
                   dble((2*l-1)*(l+m)))
              Ppleg(l,m,0)=Ppleg(l-2,m-2,0)*c        &
                  -d*gaus*Ppleg(l-1,m-2,0)+e*gaus*Ppleg(l-1,m,0)     
           enddo
        enddo
! Calcul des pseud0-derivees des polynomes de Legendre
      tsq1=gaus*tsqi
      do  m=0,Mmax
         tsq2=tsq1*dble(m)
         Qpleg(m,m,0)=-tsq2*Ppleg(m,m,0)
         do  l=m+1,ll
            a=dsqrt(dble((l-m)*(l+m+1)))
          Qpleg(l,m,0)=-a*Ppleg(l,m+1,0)-tsq2*Ppleg(l,m,0)
       enddo
    enddo
!--------------------------------------------------------------------
!------------------>  duplication pour N=0
    do  m=0,Mmax
       do l=m,ll
          Pmleg(l,m,0)=Ppleg(l,m,0)
          Qmleg(l,m,0)=Qpleg(l,m,0)
       enddo
    enddo
!------------------> Calcul des polynomes de legendre generalises N
    DO  N=0,Nmax-1
       N1=N+1
       do  m=0,Mmax
          tsq=tsqi*(dble(N)*gaus-dble(m))
          tsq1=tsqi*(dble(N)*gaus+dble(m))
          tsq2=tsqi*gaus
          Mbeg=max0(m,N1)
!*** relation de recurence sur N 
          do  l=Mbeg,ll 
             a=dsqrt(dble((l-N)*(l+N+1)))
             aa=1.d0/a    
             Ppleg(l,m,N1)=aa*(Qpleg(l,m,N)+tsq*Ppleg(l,m,N))
             Pmleg(l,m,N1)=-aa*(Qmleg(l,m,N)+tsq1*Pmleg(l,m,N))
             Qpleg(l,m,N1)=-a*Ppleg(l,m,N)+(tsq+tsq2)*Ppleg(l,m,N1)
             Qmleg(l,m,N1)=a*Pmleg(l,m,N)+(tsq1+tsq2)*Pmleg(l,m,N1)
          enddo
       enddo
    enddo
!-----------------------------------------------------------------------
 endif
!----------------------------------------------------------------------
 do l=0,ll
    coef=sqrt(float((2*l+1))/2.0)
    do n=-2,2
       do m=-2,2
          in=abs(n)
          im=abs(m)
          if (m.gt.0) then
             if (n.gt.0) then
                Val=Ppleg(l,im,in)
             else
                val=Pmleg(l,im,in)
             endif
          else
             if (n.gt.0) then
                val=(-1)**(n+m)*Pmleg(l,in,im)
             else
                val=(-1)**(m+n)*Ppleg(l,im,in)
             endif
          endif
          if (dsin(theta).lt.0.0d0) val=(-1)**(n+m)*val
          PlmN(n,m,l)=val/coef
       enddo
    enddo
 enddo
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
end subroutine legendre
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
!-----------------------------------------------------------  
  subroutine resample(ya,t0,dta,na,yb,dtb,nb)
!-----------------------------------------------------------  
    use module_spline
    implicit none
    integer            , intent(in) :: na,nb
    doubleprecision, dimension(na), intent(in) :: ya
    doubleprecision               , intent(in) :: dta,dtb,t0
    real, dimension(nb), intent(out):: yb
!
    doubleprecision, dimension(na) :: xa,y2
    integer             :: i
    doubleprecision                :: yp1,ypn,x
!
    do i=1,na
       xa(i)=(i-1)*dta+t0
    enddo
    yp1=(ya(2)-ya(1)    )/dta
    ypn=(ya(na)-ya(na-1))/dta
    call spline(xa,ya,yp1,ypn,y2)
    do i=1,nb
       x=(i-1)*dtb
       if (x<=xa(na) .and. x>=xa(1)) then
          yb(i)=real(splint(xa,ya,y2,x))
       else
          yb(i)=0.0
       endif
    enddo
!-----------------------------------------------------------  
  end subroutine resample
!-----------------------------------------------------------  
!--------------------------------------------------------
      subroutine apply_source_filter(data)
!--------------------------------------------------------
        use def_gparam
        use global_main_param, only: NBT,dt,nbcomp
        implicit none
        real(SP), dimension(:,:), intent(inout) :: data
!
        real(SP), dimension(NBT) :: tmp
        integer :: ic,ir,i,jj
!
        do ic=1,nbcomp
           tmp(:)=0.d0
           do i=1,NBT
              do jj=1,i
                 tmp(i)=tmp(i)+data(jj,ic)*source_g(i-jj+1,1,1)
              enddo
           enddo
           data(:,ic)=tmp(:)*dt
        enddo
!--------------------------------------------------------
      end subroutine apply_source_filter
!--------------------------------------------------------
!--------------------------------------------------------------------------------
end module modes
!--------------------------------------------------------------------------------
