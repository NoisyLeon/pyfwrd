!======================
 module time_function
!======================
! use data_main
! use data_temp
!
 implicit none
 contains
!
 double precision function ricker(t,f0,t0)
!-----------------------------------
  use global_main_param, only: amp_source
!
 doubleprecision        :: t,aa,PI,f0,t0
!
 PI = 4.d0*datan(1.d0)
 aa  = (PI*f0)**2

 ricker = - amp_source*(1.d0 - 2.d0*aa*(t - t0)**2)*dexp(-aa*(t - t0)**2)
!
 end function ricker
!
!
 subroutine def_timefunc(source_type,nstep,dt,f0 &
                       ,f1h,f2h,f3h,f4h,t0,g,fmax)
!-----------------------
 use def_gparam
 use global_main_param, only: amp_source
!
  character(len=6), intent(in) :: source_type
  integer, intent(in) :: nstep
  doubleprecision, intent(in)  :: dt,f0,f1h,f2h,f3h,f4h,t0
  real, dimension(:), intent(out) :: g
  doubleprecision, intent(out)  :: fmax
!
 integer		:: it0,istep,j,nstep2,i,ic,istat,itsource_max
 doubleprecision        :: maxg,wt,freq,t1,t2,tmax,seuil
 complex*16             :: dphi
 complex*16, dimension(:), allocatable :: spectre
!
 g(:)   = 0.0d0

 select case (source_type)
 case('ricker')
    do istep = 1,nstep
       g(istep) = ricker(dt * istep,f0,t0)
    enddo
    fmax=2.6d0*f0
 case('heavis')
    fmax=f4h
    nstep2=int(2.d0**(int(log(dble(nstep))/log(2.d0))+1))
!    print*,'nstep,nstep2:',nstep,nstep2
    allocate(spectre(nstep2))
    spectre(:)=cmplx(0.d0,0.d0)
    do j=1,nstep2
         if (j<=nstep2/2) then
            freq=(j-1)/(dt*nstep2)
         else if (j==nstep2/2+1) then
            freq=1/(2.d0*dt)
         else
            freq=-(nstep2-j+1)/(dt*nstep2)        
         endif
         dphi=exp(-2.d0*pi*freq*t0*cmplx(0.d0,1.d0))
         call wtcoef(abs(freq),f1h,f2h,f3h,f4h,wt)
!         if (j/=1) spectre(j)=wt*dphi/cmplx(0._DP,freq*2._DP*PI)
         if (j/=1) spectre(j)=wt*dphi
      enddo   
      call dfour1(spectre,nstep2,1)
      g(:)=amp_source*real(spectre(1:nstep))/nstep2/dt 
!on met les premier pas de temps a zero:
      tmax=nstep2*dt
      t1=0.d0
      t2=t0/5.d0
      do i=1,nstep
         call wtcoef((i-1)*dt,t1,t2,tmax,tmax,wt)
         g(i)=g(i)*wt
      enddo      
!
      deallocate(spectre,stat=istat)
      if (istat/=0) stop 'time_function deallocate error'
 case default
    stop 'time_function ce source_type n''est pas prevu!'
 end select
!
! calcul du temps de fin de la source
!
 seuil=1.d-3*maxval(abs(g(:)))
 itsource_max=-1
 ic=50 !50 pas de temps
 it0=t0/dt+1
 do i=it0,nstep-ic
    if (itsource_max<0) then
       if (maxval(abs(g(i:i+ic)))<=seuil) itsource_max=i
    endif
 enddo 
 if (itsource_max<=0) itsource_max=nstep-ic

 open(10,file='source.gnu',status='unknown')
 do istep = 1,nstep
        write(10,*) dt*istep,real(g(istep))
 enddo
 close(10)
!
 contains
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
 end subroutine def_timefunc
!
!


!=========================
 end module time_function
!=========================
