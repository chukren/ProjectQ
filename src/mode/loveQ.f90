  program  loveq
  ! Integration on depth to calculate the local Q for Love
  ! waves. According to eqn (12) in the 2009 GJI paper by Zhou.
  ! The "-" sign in first term should be "+" in this paper. 
  !
  ! 05/19/2009 chukren Virginia Tech
  
  implicit none
  
  integer :: max_layer
  integer, parameter :: mlayer = 140  ! the number of PREM layers
  integer :: istat, i, cont

  character*1  :: mode            ! dummy paras
  real :: xvphs, vphs             ! phase velocity 
  real :: xvgrp, vgrp=0.0         ! group velocity
  real :: omega                   ! frequency
  real :: depth                   ! tempo para
  real :: Q = 0.                  ! Q value for Love/Rayleigh
  real :: xped, period, ped=0.0   ! period of surface wave
  real :: compare, mini = 20.     ! control para
  real :: const, k 
  real :: pi = 3.1415926535898
  real :: radii, vph, eta, ord, f                 ! dummy papras

  real, dimension(80) :: rho, vp, vsv, vsh  ! -- 80 is set for redundancy 
  real, dimension(80) :: W, dW                    ! eigen function
  real, dimension(80) :: Qkappa, Qmu              ! qulity factor
  real, dimension(80) :: radius                   ! radius of Earth
  real, dimension(80) :: coef_alph                ! scattering coefficient
  
  
  !-------------- reading data  --------------------
  write(*,*) 'Input the period:'
  read(*,*) period

  open(9,file='temp.Q')

  !--------- counting layer of model --------
  open(unit=20,file='function.dW.0',iostat=istat)
  read(20,*) 
  i = 0
  do
    read(20,*,iostat=istat) 
    if (istat /= 0) exit
    i = i + 1 
  end do
  max_layer = i
  close(20)

  ! -------------- reading eigen functions ---------------
  open(unit=10,file='function.W.0',iostat=istat)
  read(10,*)  ! skip the head line 
  do i = 1, max_layer
    read(10,*)  W(i), depth
    radius(i) = 6.371E+06 - depth * 1.0E+03   
  end do
  close (10)

  open(unit=20,file='function.dW.0',iostat=istat)
  read(20,*) 

  do i = 1, max_layer
    read(20,*) dW(i), depth
  end do

  close (20)

  ! -------------- reading the model ------------------------ 
  open(unit=30,file='love.out',iostat=istat)

  do i = 1, 5
    read(30,*)
  end do

  do 
    read(30,*) cont
    if (cont >= mlayer - max_layer) then
      do i = 1, max_layer
        read(30,*) cont, radii, rho(i), vp(i), vph, vsv(i), vsh(i), eta, Qmu(i), Qkappa(i) 
      end do
      exit
    end if
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
  
  omega = 2.0 * pi / ped                        ! angular freq
  vphs = vphs / 6371000.                        ! rad/s
  vgrp = vgrp / 6371000.                        ! rad/s
  k = omega / vphs    

  ! intergarl alone the depth from 0 to 1 

 
  do i = 1, max_layer-1             
      call LL_coef_alph(radius(i), rho(i), vsh(i), W(i), dW(i), k, coef_alph(i))
      radii = (radius(i+1) + radius(i)) / 2.0 
      Q = Q + coef_alph(i) * radii**2 * (radius(i+1) - radius(i)) / Qmu(i)  
  end do  

  const = vphs * vgrp / ( -2. * omega**2 )
  Q = 1. / (Q * const)
  
  !write(*,*) 'period = ', ped, 'const = ', const, 'Q = ', Q     
  vphs = vphs * 6371.   ! convert to km/s
 
  write(9,1010)Q,vphs
   1010 format ('Q = 'F10.2'  vphs = 'F10.4)
  close(9)
  end program loveq




  subroutine LL_coef_alph(radii, rho, vsh, W, dW, k, coef)

    ! calculate L-L scattering coefficient 
    ! WARNING: eigenfunctions have been scaled to the correct 
    ! normalization and have been multiplied by a factor of 1E09 
    ! for plotting purposes, divide the eigenfunction and
    ! derivatives by 1E09 in calculating the Love or Rayleigh Q factors.

    real :: radii, rho, vsh, W, dW 
    real :: k, coef, para
    real :: term1, term2
   
    para = -2.0 * rho * vsh**2            ! vsh in rad/s
    W = W * 1.0E-09
    dW = dW * 1.0E-09
    term1 = para * ( dW - W / radii ) * ( dW - W / radii )
    term2 = para * k**2 * W**2 / radii**2
    coef = para * ( dW - W / radii ) * ( dW - W / radii )  + para * k**2 * W**2 / radii**2    
    !write(*,*) 'r=',radii,'V=',vsh,'W=',W,'dW=',dW,'k=',k,'coef = ', coef
    !write(*,*) 'term1 = ',term1,'term2 = ',term2
  end subroutine
