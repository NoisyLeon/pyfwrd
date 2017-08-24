!-----------------------------------------------------------------
module module_spline
!du Numerical recipes
!-----------------------------------------------------------------
  implicit none
  public :: spline,splint
  private
contains
!------------------------------------------------------------------
SUBROUTINE spline(x,y,yp1,ypn,y2)
!------------------------------------------------------------------
  USE nrtype; USE nrutil, ONLY : assert_eq
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: x,y
  REAL(DP), INTENT(IN) :: yp1,ypn
  REAL(DP), DIMENSION(:), INTENT(OUT) :: y2
  INTEGER(I4B) :: n
  REAL(DP), DIMENSION(size(x)) :: a,b,c,r
  n=assert_eq(size(x),size(y),size(y2),'spline')
  c(1:n-1)=x(2:n)-x(1:n-1)
  r(1:n-1)=6.0_DP*((y(2:n)-y(1:n-1))/c(1:n-1))
  r(2:n-1)=r(2:n-1)-r(1:n-2)
  a(2:n-1)=c(1:n-2)
  b(2:n-1)=2.0_DP*(c(2:n-1)+a(2:n-1))
  b(1)=1.0_DP
  b(n)=1.0_DP
  if (yp1 > 0.99e30_DP) then
     r(1)=0.0_DP
     c(1)=0.0_DP
  else
     r(1)=(3.0_DP/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
     c(1)=0.5_DP
  end if
  if (ypn > 0.99e30_DP) then
     r(n)=0.0_DP
     a(n)=0.0_DP
  else
     r(n)=(-3.0_DP/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
     a(n)=0.5_DP
  end if
  call tridag_ser(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
!------------------------------------------------------------------
END SUBROUTINE spline
!------------------------------------------------------------------
!------------------------------------------------------------------
FUNCTION splint(xa,ya,y2a,x)
!------------------------------------------------------------------
  USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
  REAL(DP), INTENT(IN) :: x
  REAL(DP) :: splint
  INTEGER(I4B) :: khi,klo,n,i
  REAL(DP) :: a,b,h
  n=assert_eq(size(xa),size(ya),size(y2a),'splint')
  klo=max(min(locate(xa,x),n-1),1)
  khi=klo+1
  h=xa(khi)-xa(klo)
  if (h == 0.0_DP) then
        do i=1,n
           print*,i,xa(i),khi,klo
        enddo
	call nrerror('bad xa input in splint')
  endif
  a=(xa(khi)-x)/h
  b=(x-xa(klo))/h
  splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_DP
!------------------------------------------------------------------
END FUNCTION splint
!------------------------------------------------------------------

!------------------------------------------------------------------
SUBROUTINE tridag_ser(a,b,c,r,u)
!------------------------------------------------------------------
  USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: a,b,c,r
  REAL(DP), DIMENSION(:), INTENT(OUT) :: u
  REAL(DP), DIMENSION(size(b)) :: gam
  INTEGER(I4B) :: n,j
  REAL(DP) :: bet
  n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_ser')
  bet=b(1)
  if (bet == 0.0_DP) call nrerror('tridag_ser: Error at code stage 1')
  u(1)=r(1)/bet
  do j=2,n
     gam(j)=c(j-1)/bet
     bet=b(j)-a(j-1)*gam(j)
     if (bet == 0.0_DP) &
          call nrerror('tridag_ser: Error at code stage 2')
     u(j)=(r(j)-a(j-1)*u(j-1))/bet
  end do
  do j=n-1,1,-1
     u(j)=u(j)-gam(j+1)*u(j+1)
  end do
!------------------------------------------------------------------
END SUBROUTINE tridag_ser
!------------------------------------------------------------------
!-----------------------------------------------------------------
FUNCTION locate(xx,x)
!-----------------------------------------------------------------
  USE nrtype
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: xx
  REAL(DP), INTENT(IN) :: x
  INTEGER(I4B) :: locate
  INTEGER(I4B) :: n,jl,jm,ju
  LOGICAL :: ascnd
  n=size(xx)
  ascnd = (xx(n) >= xx(1))
  jl=0
  ju=n+1
  do
     if (ju-jl <= 1) exit
     jm=(ju+jl)/2
     if (ascnd .eqv. (x >= xx(jm))) then
        jl=jm
     else
        ju=jm
     end if
  end do
  if (x == xx(1)) then
     locate=1
  else if (x == xx(n)) then
     locate=n-1
  else
     locate=jl
  end if
!-----------------------------------------------------------------
END FUNCTION locate
!-----------------------------------------------------------------
!-----------------------------------------------------------------
end module module_spline
!-----------------------------------------------------------------
