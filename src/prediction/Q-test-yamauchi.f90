! Generate shear wavespeed profile from a thermal model
! by Y.Y. Ruan 01/16／2017
!=======================================================

  program Qtest
  use thermodel
  implicit none
  
  real, parameter :: waterdepth = 3.163

  real :: tmax, age, seisfreq, temp, pressure, vs, vp, Q
  real :: depth, thick, depth_real
  real :: vs_melt, Q_melt, dvv, Q_gbs
  real :: vs_yamauchi, Q_yamauchi, T_melt, viscosity_yamauchi
  real :: viscosity, shear_modulus
  real :: J1, J2, G1, G2, Q_b, vs_b

  integer :: i

  character (len=120) :: fout_all, fout_solidus, fout_modulus
  character (len=120) :: fout_background, fout_viscosity
      
  call setupzeroage(zeroaget,zerocoeff,pottemp,adbatgrad,srate)

  tmax = zeroaget(201)
  write(*,*)'Maximum temperature:', tmax

  write(*,*)'Input seafloor age (Ma):'
  read(*,*)age

  write(*,*)'Input seismic wave frequency (Hz):'
  read(*,*)seisfreq

  write(*,*)'Output all prediction file name:'
  read(*,*)fout_all

  write(*,*)'Output solidus file name:'
  read(*,*)fout_solidus

  write(*,*)'Output modulus file name:'
  read(*,*)fout_modulus 

  write(*,*)'Output background file name:'
  read(*,*)fout_background

  open(22, file=fout_all)
  open(23, file=fout_solidus)
  open(24, file=fout_modulus)
  open(25, file=fout_background)

  do i = 1,200
      depth = (i-1)    

      ! Calculate temperature profile
      call tempest(zeroaget, zerocoeff, age, depth, tmax, temp)   

      call hydrostatic_pressure(depth, pressure)
      call unrelaxed_shear_modulus(temp, pressure, shear_modulus)

      ! background viscosity
      call shear_viscosity(temp, pressure, 1.0, viscosity)

      ! consider the pre-melting effect on relaxation peak
      call yamauchi_takei_2016(depth, temp, seisfreq, vs_yamauchi, Q_yamauchi, T_melt, &
                               viscosity_yamauchi, J1, J2, Q_b, vs_b)

      depth_real = depth + waterdepth

      write(23, *) T_melt, depth_real
      write(24, *) shear_modulus/1.0E+09, depth_real

      write(22, 1000) depth_real, zeroaget(i), temp, vs_yamauchi, Q_yamauchi, viscosity_yamauchi
      ! Note the viscosity is the one whitout melting effects, but backgroundVs_b and Q_b 
      ! are calculated considering the melting effects
      write(25, 1000) depth_real, zeroaget(i), temp, vs_b, Q_b, viscosity
  enddo    
 
1000 format(5(f12.3,1x), E16.7) 

  close(11)
  close(22)
  close(23)
  close(24)
  close(25)
  end
