  program  rayleighq
  ! Integration on depth to calculate the local Q for Rayleigh 
  ! waves. According to the eq. 12 by Zhou 2009 GJI
  !
  ! 05/19/2009 chukren Virginia Tech
  
  implicit none
  
  integer :: max_layer, mlayer = 140      ! the number of layer
  integer :: istat, i, cont

  character*1  :: mode            ! dummy paras
  real :: xvphs, vphs             ! phase velocity 
  real :: xvgrp, vgrp=0.0         ! group velocity
  real :: omega                   ! frequency
  real :: depth                   ! tempo para
  real :: Q = 0.                  ! Q value for Love/Rayleigh
  real :: xped, period, ped=0.0   ! period of surface wave
  real :: compare, mini = 20.     ! control para
  real :: gamma                   ! control para
  real :: const, k 
  real :: pi = 3.14159265
  real :: radii, vph, eta, ord, f                 ! dummy papras
  real :: tmp1, tmp2                              ! temporary papras

  real, dimension(80) :: rho, vp, vsv, vsh  ! -- 80 is set for redundancy 
  real, dimension(80) :: U, dU, V, dV             ! eigen function
  real, dimension(80) :: Qkappa, Qmu              ! qulity factor
  real, dimension(80) :: radius                   ! radius of Earth
  real, dimension(80) :: coef_alph, coef_beta     ! scattering coefficient
  
  
  !-------------- reading data  --------------------
  write(*,*) 'Input the period:'
  read(*,*) period

  open(9,file='temp.Q')

  !--------- counting layer of eigen functions --------
  open(unit=20,file='function.U.0',iostat=istat)
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
  open(unit=10,file='function.U.0',iostat=istat)
  read(10,*)  ! skip the head line 
  do i = 1, max_layer
    read(10,*)  U(i), depth
    radius(i) = 6.371E+06 - depth * 1.0E+03   
  end do
  close (10)

  open(unit=20,file='function.dU.0',iostat=istat)
  read(20,*) 
  do i = 1, max_layer
    read(20,*) dU(i), depth
  end do
  close (20)

  open(unit=30,file='function.V.0',iostat=istat)
  read(30,*) 
  do i = 1, max_layer
    read(30,*) V(i), depth
  end do
  close (30)

  open(unit=40,file='function.dV.0',iostat=istat)
  read(40,*) 
  do i = 1, max_layer
    read(40,*) dV(i), depth
  end do
  close (40)

  ! -------------- reading the model ------------------------ 
  open(unit=50,file='rayleigh.out',iostat=istat)

  do i = 1, 5
    read(50,*)
  end do

  do 
    read(50,*) cont
    if (cont >= mlayer - max_layer) then
      do i = 1, max_layer
        read(50,*) cont, radii, rho(i), vp(i), vph, vsv(i), vsh(i), eta, Qmu(i), Qkappa(i) 
      end do
    exit
    end if
  end do

  do i = 1, 10
    read(50,*) 
  end do
  
  !---------------------------------------------------------------------
  !  searching the (frequency) period, phase velocity and group velocity
  !---------------------------------------------------------------------
  do
    read(50,*,iostat=istat) cont, mode, ord, f, xped, xvphs, xvgrp
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
  close (50)

  ! ---------------- calculate local Q   ---------------------------------
  
  omega = 2.0 * pi / ped                        ! angular freq
  vphs = vphs / 6371000.                        ! rad/s
  vgrp = vgrp / 6371000.                        ! rad/s
  k = omega / vphs    

  ! intergarl alone the depth from 0 to 1 

  do i = 1, max_layer-1             
      call RR_coef(radius(i), rho(i), vp(i), vsv(i), U(i), dU(i), V(i), dV(i), k, coef_alph(i), coef_beta(i))

      radii = (radius(i+1) + radius(i)) / 2.0 
      gamma = 4.0 * vsv(i)**2 / 3.0 / vp(i)**2 
      tmp1 = (1.0 - gamma)*coef_alph(i) / Qkappa(i) 
      tmp2 = (coef_beta(i) + gamma*coef_alph(i)) / Qmu(i)

      Q = Q + (tmp1 + tmp2) * radii**2 * (radius(i+1) - radius(i))   
  end do  

  const = vphs * vgrp / ( -2. * omega**2 )
  Q = 1. / (Q * const)
  
  vphs = vphs * 6371.   ! convert to km/s
 
  write(9,1010)Q,vphs
   1010 format ('Q = 'F10.2' vphs = 'F10.4)
  close(9)
  end program rayleighq


  subroutine RR_coef(radii, rho, vp, vsv, U, dU, V, dV, k, coefa, coefb)
    ! calculate R-R scattering coefficient 

    ! WARNING: eigenfunctions have been scaled to the correct 
    ! normalization and have been multiplied by a factor of 1E09 
    ! for plotting purposes, divide the eigenfunction and
    ! derivatives by 1E09 in calculating the Love or Rayleigh Q factors.
    implicit none
    real :: radii, rho, vp, vsv, U, dU, V, dV 
    real :: k, coefa, coefb, paraa, parab
    real :: term
   
    paraa = -2.0 * rho * vp**2            
    parab = -2.0 * rho * vsv**2
    U = U * 1.0E-09
    dU = dU * 1.0E-09
    V = V * 1.0E-09
    dV = dV * 1.0E-09

    coefa = paraa * ( dU  +  2.0*U/radii - k * V / radii )**2 
    
    term = -2.0 * ( dU  + 2.0*U/radii - k * V / radii )**2
    term = term + 2 * dU**2 + (2*U - k*V)**2 / radii**2
    term = term + (dV - V / radii + k*U / radii)**2
    term = term + k**2 * V**2 / radii**2 
    coefb = parab * term
  end subroutine
