  program  lovephase
  ! looking for the phase velocity at periods
  ! 05/19/2009 chukren Virginia Tech
  
  implicit none
  
  integer :: max_layer, mlayer = 164      ! the number of layer
  integer :: istat, i, j, cont

  character*10 :: firstline1, firstline2
 
  character*1  :: mode            ! dummy paras
  real :: xvphs, vphs             ! phase velocity 
  real :: xvgrp, vgrp             ! group velocity
  real :: omega                   ! frequency
  real :: depth                   ! tempo para
  real :: Q = 0.                  ! Q value for Love/Rayleigh
  real :: xped, period, ped       ! period of surface wave
  real :: compare, mini = 20.     ! control para
  real :: const, k 
  real :: pi = 3.14159265
  real :: radii, vph, eta, ord, f                 ! dummy papras

  real, dimension(80) :: rho, gama, vp, vsv, vsh  ! -- 80 is set for redundancy 
  real, dimension(80) :: W, dW                    ! eigen function
  real, dimension(80) :: omegaa, omegab           ! scattering coefficent
  real, dimension(80) :: Qkappa, Qmu              ! qulity factor
  real, dimension(80) :: radius                   ! radius of Earth
  real, dimension(80) :: coef_alph                ! scattering coefficient
  
  
  !-------------- reading data  --------------------
  write(*,*) 'Input the period:'
  read(*,*) period

  open(9,file='temp.phase')



  ! -------------- reading the model ------------------------ 
  open(unit=30,file='love.out',iostat=istat)

  do i = 1, 5
    read(30,*)
  end do

  do i = 1, mlayer
    read(30,*) 
  end do

  do i = 1, 10
    read(30,*) 
  end do
  
  !---------------------------------------------------------------------
  !  searching the (frequency) period, phase velocity and group velocity
  !---------------------------------------------------------------------
  do
    read(30,*,iostat=istat) cont, mode, ord, f, xped, xvphs, xvgrp
    if ( istat /= 0 ) exit 
    compare = abs(xped - period)

    if ( cont == 0 .and. compare < mini ) then
        mini = compare
        ped  = xped
        vphs = xvphs * 1.0E+03
        vgrp = xvgrp * 1.0E+03
    end if

    if (mini < 1.E-08) exit

  end do

  if (mini >= 20.) then
      write(*,*) "sorry, can not locate the closest frequency you requested"
      call exit  
  end if
  close (30)

  ! ---------------- calculate local Q   ---------------------------------
  
  vphs = vphs / 1000.                        ! rad/s
  vgrp = vgrp / 1000.                        ! rad/s
 
  write(9,1010)vphs
   1010 format ('vphs = 'F10.4)
  close(9)
  end program lovephase


