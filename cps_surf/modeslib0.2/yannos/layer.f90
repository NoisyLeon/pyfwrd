!---------------------------------------------------------------------------
module layer
!module qui permet de limiter le stockage des fonctions propres a qq regions
!---------------------------------------------------------------------------
  use def_gparam
  implicit none
!tout est public
  public
!nombre de couches a selectionne
  integer :: nb_lay
!nombre de couches qui vont etre stockees:
  integer :: nbcou_lay
!indice encadrant les couches selectionnes
  integer, dimension(:), allocatable :: i_la1,i_la2
!rayon normaliser des couches:
  real(DP), dimension(:), allocatable :: rn_la1,rn_la2
!rayon  des couches:
  real(SP), dimension(:), allocatable  :: r_la1,r_la2                                        
!=============================
contains
!=============================
!---------------------------------------------------------------------------
  subroutine read_layers(unit)
!---------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: unit
!
    integer :: i
!
    read(unit,*)
    read(unit,*) nb_lay
    allocate(i_la1(nb_lay),i_la2(nb_lay),rn_la1(nb_lay),rn_la2(nb_lay),r_la1(nb_lay),r_la2(nb_lay))
    read(unit,*)
    do i=1,nb_lay
       read(unit,*) r_la1(i),r_la2(i)
    enddo
!---------------------------------------------------------------------------
  end subroutine read_layers
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  subroutine init_layer(mod)
!---------------------------------------------------------------------------
    use earth_modele
    implicit none
!
    type(modele), intent(in) :: mod
!
    real(DP) :: ra
    integer  :: i,j
!
    ra=mod%r(mod%nbcou)
    do i=1,nb_lay
       rn_la1(i)=r_la1(i)/ra
       rn_la2(i)=r_la2(i)/ra
    enddo
    do i=2,mod%nbcou
       do j=1,nb_lay
          if (r_la1(j) >= mod%r(i-1) .and. r_la1(j) < mod%r(i) ) &
                                                     i_la1(j)=i-1
          if (r_la2(j) >  mod%r(i-1) .and. r_la2(j) <=mod%r(i) ) &
                                                     i_la2(j)=i
       enddo
    enddo
!
    nbcou_lay=0
    do i=1,nb_lay
       nbcou_lay=nbcou_lay+i_la2(i)-i_la1(i)+1
    enddo
    if (nbcou_lay<=0.or.nbcou_lay>mod%nbcou) then
       print*,'nbcou_lay=',nbcou_lay
       stop 'Erreur fatale dans init_layer du module layer!'
    endif
!---------------------------------------------------------------------------
  end subroutine init_layer
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
end module layer
!---------------------------------------------------------------------------
