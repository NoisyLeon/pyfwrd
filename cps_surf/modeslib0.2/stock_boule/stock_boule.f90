!-----------------------------------------------------------------------
program stock_freq
!-----------------------------------------------------------------------
  use def_gparam
  use sol_ana
!-----------------------------------------------------------------------
  implicit none
!
  integer :: Lmax,NCM,NMAX
  integer :: l,ip,i,nz

  real(DP), dimension(2) :: fdeb
  real(DP) :: df,r,vp,vs,rho,lamb,mu,fmax,fmin
!
  character(len=20)  :: junk
  character(len=100) :: filer,filel
  logical  ::  flag_deb=.true.,flag,creat=.true.
!
  
!!$  rho=3000.0_DP
!!$  vp=8500.0_DP
!!$  vs=5000.0_DP
!!$  r=6371000.0_DP
!!$  mu=rho*vs**2
!!$  lamb=rho*vp**2-2.00_DP*mu
!!$!
!!$  filer='totoR'
!!$  filel='totoL'
!!$  fmin=0.0000010_DP
!!$  fmax=0.02_DP             !filtre jusqu'a 20s
  open(12,file='stock_boule.dat',status='old')
  read(12,'(a20,1x,i4)')junk,Lmax
  print*,junk,Lmax
  read(12,'(a20,1x,i4)')junk,NCM
  print*,junk,NCM
  read(12,'(a20,1x,i4)')junk,NMAX
  print*,junk,NMAX
  read(12,'(a20,1x,f7.1)')junk,rho
  print*,junk,rho
  read(12,'(a20,1x,f7.1)')junk,vp
  print*,junk,vp
  read(12,'(a20,1x,f7.1)')junk,vs
  print*,junk,vs
  read(12,'(a20,1x,E9.7)')junk,r
  print*,junk,r
  read(12,'(a20,1x,E9.7)')junk,fmin
  print*,junk,fmin
  read(12,'(a20,1x,E9.7)')junk,fmax
  print*,junk,fmax
  read(12,'(a20,1x,a100)')junk,filer
  print*,junk,filer
  read(12,'(a20,1x,a100)')junk,filel
  print*,junk,filel
  close(12)
!
  mu=rho*vs**2
  lamb=rho*vp**2-2.00_DP*mu

  call set_boule_param(r,vp,vs,rho,lamb,mu,fmin,fmax,NCM,NMAX,LMAX &
                      ,filer,filel)
  call get_freq_param(df_out=df)
  print*,'df=',df
!
  b_l: do l=0,Lmax
!test
!  b_l: do l=16,16
!fin test
     if (l.gt.0) then
        if (flag_deb.or.n_zero(l-1,1).ne.0     &
             .or.n_zero(l-1,2).ne.0) then
           flag=.true.
        else
           flag=.false.
        endif
     else
        flag=.true.
     endif
     if (flag) then
        do i=1,2
           if (l >1 ) then
              if (f(1,l-1,i) > 1.0E-12 ) then
                 fdeb(i)=max(fmin,f(1,l-1,i)-1000.0_DP*df)
              else
                 fdeb(i)=fmin 
              endif
           else
              fdeb(i)=fmin 
           endif
           if (i==1) then
              print*,'----------------R fdeb=',fdeb(i),' l=',l
           else
              print*,'----------------L fdeb=',fdeb(i),' l=',l
           endif
        enddo
        call get_0_sl(l,det1,fdeb(1),f(:,l,1),n_zero(l,1),1)
        call get_0_sl(l,det2,fdeb(2),f(:,l,2),n_zero(l,2),2)
        if (n_zero(l,1).ne.0.or.n_zero(l,2).ne.0) flag_deb=.false.
        do nz=1,n_zero(l,1)
           call write_fctp(f(nz,l,1),l,nz,1)
        enddo
        do nz=1,n_zero(l,2)
           call write_fctp(f(nz,l,2),l,nz,2)
        enddo
     endif
  enddo b_l
  call close_boule
  print*,'fin!'
!-----------------------------------------------------------------------
end program stock_freq
!-----------------------------------------------------------------------
