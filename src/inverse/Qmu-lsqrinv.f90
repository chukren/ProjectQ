!=======================================================================
!  Inverting 1D attenuation (Q) as function of depth using surface wave 
!  attenuation coefficients inverted using tele-seismic Rayleigh waves 
!
!  linear least squares method.
!  Nov 21 2013 chukren Brown  
!     Rewrite F77 code from Donald W. Forsyth
!
!  Mar 01 2014 chukren Brown  
!     Changed code layout 
!
!  Mar 09 2014 chukren Brown 
!     Add damping Cmm
!
!  Mar 12 2014 chukren Brown
!     Add term to stablize the model variations 1/Cmm * (m-m0)
!
!  May 15 2016 chukren Princeton 
!     Modified for the 1D Q inversion for the Juan de Fuca plate
!     Add the second-order Tichonov regularization  
!=======================================================================
module  const
  implicit none
  real, parameter :: pi  = 3.14159265359
  real, parameter :: eps = 1.0E-05
  real, parameter :: ZERO = 1.0E-30
  real, parameter :: Qmu_min = 30.
  real, parameter :: Qmuinv_max = 0.1
  real, parameter :: Qmuinv_min = 0.002
  integer, parameter :: ndat = 16       ! number of data (16 period band )
  integer, parameter :: nvar = 21       ! number of model variables (21 layers)
  real, parameter :: smooth_mu = 0.5    ! 0.5,  0.0 no application; 1.0 optimum damping
end module

!===
!   Linear inversion:
!     G*m = d  
!       m = (G^T*G)^(-1)*G^T*d = gtginv * gtd
!===
program lininv
  use const
  implicit none
  integer :: i, j, k, icount, n, irow, icol

  real :: dummy1
  real :: rank
  integer :: eof
  integer :: ndata    ! number of used data 
  real, dimension(1:nvar) :: md_00, md, change, mdo
  real, dimension(1:nvar) :: varm
  real, dimension(1:ndat) :: obsd, stdv, pred

  double precision :: b, dpx
  double precision :: d_variance, dstdv
  double precision :: sumsq, rms, resid_rms
  double precision :: dt1, dt2
  double precision :: sum_smth_c
  double precision :: sum_smth_g
  double precision :: smooth_coef
  double precision :: tempcov

  integer, dimension(1:nvar) :: indx
  double precision, dimension(1:ndat) :: delta_d
  double precision, dimension(1:ndat) :: misfit
  double precision, dimension(1:ndat) :: dataimp
  double precision, dimension(1:nvar) :: tempimp
  double precision, dimension(1:nvar) :: stdv_md
  double precision, dimension(1:nvar,1:nvar) :: covinv  ! Cmm, damping, cov of md stdv.
  double precision, dimension(1:ndat,1:nvar) :: g
  double precision, dimension(1:nvar,1:nvar) :: rough
  double precision, dimension(1:nvar) :: gtd, gtdcmm, cmminv_dm          ! G^T*d
  double precision, dimension(1:nvar,1:nvar) :: gtg, gtginv, gtg_reserve ! G^T*G
  double precision, dimension(1:nvar,1:nvar) :: fcov, corr, resolution

  character(len=10) :: dummy2, dummy3
  character(len=120) :: f_obs
  character(len=80) :: f_d
  character(len=80) :: f_m
  character(len=120) :: f_m00
  character(len=80) :: f_g
  character(len=80) :: foutput, f_covariance, f_correlation, f_resolution
  character(len=80) :: fmisfit

  !===
  !   Read input files
  !===
  write(*,*)'Input file of obs and stdv:'
  read(*,*)f_obs
  f_obs = trim(adjustl(f_obs))
  write(*,'(A10,A15)')'read ',f_obs

  write(*,*)'Input file of data:'
  read(*,*)f_d
  write(*,'(A10,A15)')'read ',f_d

  write(*,*)'Input file of model:'
  read(*,*)f_m
  write(*,'(A10,A15)')'read ',f_m

  write(*,*)'Input file of model:'
  read(*,*)f_m00
  write(*,'(A10,A15)')'read ',f_m00

  write(*,*)'Input file of kernel:'
  read(*,*)f_g
  write(*,'(A10,A15)')'read ',f_g

  write(*,*)'file name for updated model output:'
  read(*,*)foutput
  write(*,*)'file name is ',foutput
  open(10,file=foutput,status='new')

  write(*,*)'file name for best fit prediction and misfit:'
  read(*,*)fmisfit
  write(*,*)'file name is ',fmisfit
  open(12,file=fmisfit,status='new')

  write(*,*)'file name for covariance matrix:'
  read(*,*)f_covariance
  write(*,*)'file name is ',f_covariance
  open(14,file=f_covariance,status='new')

  write(*,*)'file name for correlation matrix:'
  read(*,*)f_correlation
  write(*,*)'file name is ',f_correlation
  open(15,file=f_correlation,status='new')

  write(*,*)'file name for resolution matrix:'
  read(*,*)f_resolution
  write(*,*)'file name is ',f_resolution
  open(16,file=f_resolution,status='new')

  !===
  !   Readin g data matrix:
  !         delta_d = syn - obs
  !   and the G(kernel) matrix
  !         g(i,j) = d(syn1(i)-syn2(i))/d(md1(j)-md2(j))
  !   Both are normalized by dividing the standard deviation
  !===

  open(13,file=f_obs,status='old')
  i = 0
  do
    read(13,*,iostat=eof) n, obsd(n), stdv(n)
    if(eof /= 0) exit
    i = i + 1
  enddo
  if(i /= ndat) then
    write(*,*)'number of data is not correct in',f_obs
    stop
  endif
  close(13)

  open(13,file=f_d,status='old')
  i = 0
  do
    read(13,*,iostat=eof) n, delta_d(n)
    !write(*,*) n, delta_d(n)
    if(eof /= 0) exit
    i = i + 1
  enddo
  if(i /= ndat) then
    write(*,*)'number of data is not correct in', f_d
    stop
  endif
  close(13)


  open(13,file=f_m,status='old')
  i = 0
  do
    read(13,*,iostat=eof) n, md(n), varm(n) ! label(n), label2(n)
    !write(*,*) n, md(n), varm(n)
    if(eof /= 0) exit
    i = i + 1
  enddo
  if(i /= nvar) then
    write(*,*)'number of model para is not correct!'
    stop
  endif
  close(13)

  open(13,file=f_m00,status='old')
  i = 0
  do
    read(13, *, iostat=eof)n, md_00(n), dummy1 ! dummy2, dummy3
    if(eof /= 0) exit
    i = i + 1
  enddo
  if(i /= nvar) then
    write(*,*)'number of initial model para is not correct!'
    stop
  endif
  close(13)

  open(13,file=f_g, status='old')
  do
    read(13,*,iostat=eof)i, j, g(i,j)
    !write(*,*) i, j, g(i,j)
    if(eof /= 0) exit
  enddo
  close(13)

  !== test ==
  !stop
  !== end ==

  ndata = ndat

  !===
  !   apriori damping of model parameters, i.e., covariance matrix of model
  !   parameter. In practice use diagonals only, Cmm.
  !   Inverse for non-zero diagonals is 1/(stdv**2)
  !===
  do i = 1, nvar
    do j = 1, nvar
      covinv(i,j) = 0.0d0
    enddo
    covinv(i,i) = 1/(varm(i)**2) 
  enddo


  !===
  !   Rought matrix for minimum curvature 
  !   Note:
  !     This definition of smoothing matrix has minimum first order
  !     derivative constraint on edges
  !===
  do irow = 1, nvar
    do icol = 1, nvar
      rough(irow, icol) = 0.d0
    enddo 
  enddo
    
  ! upper left corner 
  rough(2,1) = -3
  rough(2,2) =  6
  rough(2,3) = -4
  rough(2,4) =  1

  ! lower right corner
  rough(nvar-1, nvar-3) =  1
  rough(nvar-1, nvar-2) = -4
  rough(nvar-1, nvar-1) =  6
  rough(nvar-1, nvar)   = -3

  ! middle 
  do irow = 3, nvar-2 
    icol = irow
    rough(irow, icol-2) =  1
    rough(irow, icol-1) = -4
    rough(irow, icol)   =  6
    rough(irow, icol+1) = -4
    rough(irow, icol+2) =  1
  enddo


  !===
  !   Calculate GT*G and GT*d
  !===
  do j = 1, nvar
    gtd(j) = 0.0

    do i = 1, ndata
      gtd(j) = gtd(j) + g(i,j) * delta_d(i)
    enddo

    cmminv_dm(j) = 0.d0

    do k = 1, j
      gtg(k,j) = 0.0

      do i = 1, ndata
        gtg(k,j) = gtg(k,j) + g(i,k) * g(i,j)
      enddo

      gtg(j,k) = gtg(k,j)

      ! reserve G^T*G for resolution matrix
      gtg_reserve(k,j) = gtg(k,j)
      gtg_reserve(j,k) = gtg(j,k)
    enddo

  enddo

  !===
  !   Combination of minimum length (Cmm) and curvature (rough)
  !   smooth_coef: the smooth coefficient that minimizes, in the least squares
  !   sense, the off-diagonal terms of GT*G, which should lead to minimal
  !   off-diagonal terms of covariance matrix  
  !===
  sum_smth_c = 0.0d0
  sum_smth_g = 0.0d0
  do i = 2, nvar
    do j = 1, i-1
      sum_smth_c = sum_smth_c + rough(i,j)**2
      sum_smth_g = sum_smth_g + rough(i,j)*gtg(i,j)
    enddo
  enddo
  smooth_coef = -sum_smth_g / sum_smth_c
  smooth_coef = smooth_coef * smooth_mu


  do j = 1, nvar
    do i = 1, nvar
      gtg(i,j) = gtg(i,j) + covinv(i,j) + smooth_coef*rough(i,j)  
    enddo
    
    !=== 
    !   Add to GT*d Tarantola term penalizing misfit to original starting model
    !   For on-side correction, do not penalizz changes from original starting
    !   model in the iteration
    !===
    do i = 1, nvar
      tempcov = (covinv(j,i) + smooth_coef * rough(j,i)) * (md(i) - md_00(i))
      cmminv_dm(j) = cmminv_dm(j) + tempcov
    enddo

    !gtdcmm(j) = gtd(j) - covinv(j) * (md(j) - md_00(j))
    gtdcmm(j) = gtd(j) - cmminv_dm(j)
  enddo


  !===
  !   Invert GT*G
  !===
  do i = 1, nvar
    do j = 1, nvar
      gtginv(i,j) = 0.d0
    enddo

    gtginv(i,i) = 1.d0
  enddo

  call dludcmp(gtg, nvar, nvar, indx,b)

  do j = 1, nvar
    call dlubksb(gtg, nvar, nvar, indx, gtginv(1,j))
  enddo

  !===
  !   Find change to starting model
  !===
  do i = 1, nvar
    change(i) = 0.0
    do j = 1, nvar
      change(i) = change(i) + gtdcmm(j) * gtginv(i,j)
    enddo
    mdo(i) = md(i) + change(i)
    
    !make sure Qmu will not be too small or too large, Qmu = [10, 500] 
    if (mdo(i) > Qmuinv_max) mdo(i) = Qmuinv_max
    if (mdo(i) < Qmuinv_min) mdo(i) = Qmuinv_min

    !=== resolution test 
    ! set Qmu to a fixed value to see if the low Q is required by 
    ! data
    !=== Youyi Ruan 09/22/2016 
    !if (i == 10 .or. i == 11 .or. i == 12) then
    !    mdo(i) = 35
    !endif

  enddo

  !===
  !   Find data importance: diagonal of G*(G^TG)inv*G^T
  !===
  do i = 1, ndata
    dataimp(i) = 0.0
    do j = 1, nvar
      tempimp(j) = 0.0

      do k = 1, nvar
        tempimp(j) = tempimp(j) + gtginv(j,k) * g(i,k)
      enddo

      dataimp(i) = dataimp(i) + g(i,j) * tempimp(j)
    enddo
  enddo

  !===
  !   Find resolution matrix (G^T*G + Cmm^{-1})^{-1} * G^T*G
  !===
  rank = 0.0

  do i = 1, nvar
    do j = 1, nvar
      resolution(i,j) = 0.0

      do k = 1, nvar
        resolution(i,j) = resolution(i,j) + gtginv(i,k) * gtg_reserve(k,j)
      enddo
    enddo

    rank = rank + resolution(i,i)
  enddo

  !===
  !   Find normalized residuals (linear problem with zero start
  !   so just use partial derivatives) and sum of squares of errors
  !   The effective model dimension, nvar, is unknown considering the 
  !   correlations between each other, so for simplicity, we use 
  !   (ndata - rank) as degree of freedom assuming rank is a good 
  !   approximation of number of independent information
  !===
  sumsq = 0.0
  resid_rms = 0.0
  do i = 1, ndata
    misfit(i) = 0.0

    do j = 1, nvar
      misfit(i) = misfit(i) + g(i,j) * change(j)
    enddo

    misfit(i) = delta_d(i) - misfit(i)
    sumsq = sumsq + misfit(i)**2

    resid_rms = resid_rms + delta_d(i)**2
  enddo

  rms = sqrt(sumsq/ndata)
  resid_rms = sqrt(resid_rms/ndata)

  d_variance = sumsq/(ndata - rank)
  dstdv = sqrt(d_variance)

  !===
  !   Find standard error of model parameters
  !===
  do j = 1, nvar
    stdv_md(j) = sqrt(gtginv(j,j)*d_variance)
  enddo

  !===
  !   Output results
  !===
  do i = 1, nvar
    write(10,998) i,mdo(i),stdv_md(i),md(i)
  enddo
998 format(i3,1x,e12.5,1x,e12.5,1x,e12.5)


  write(12,*)'Description'
  write(12,*)
  write(12,*)'rank ', rank
  write(12,*)'rms resid_rms', rms, resid_rms
  write(12,*)'Normalized standard deviation data', dstdv
  write(12,*)
  write(12,*)
  write(12,*)'    change      std_err'
  do i = 1, nvar
    write(12,999)change(i),stdv_md(i)
  enddo
999   format(2(f12.5,1x))
  write(12,*)
  write(12,*)

  !===
  !   Output full covariance matrix and correlation matrix
  !===
  do i = 1, nvar
    do j = 1, nvar
      fcov(i,j) = gtginv(i,j)*d_variance
      corr(i,j) = gtginv(i,j)/sqrt(gtginv(i,i)*gtginv(j,j))
    enddo
  enddo

  do i = 1, nvar
    do j = 1, nvar 
      write(14,*)i, j, fcov(i,j)
    enddo
  enddo

  do i = 1, nvar
    do j = 1, nvar
      write(15,*)i, j, corr(i,j)
    enddo
  enddo

  do i = 1, nvar
    do j = 1, nvar
      ! only output diagonal elements
      !if (i == j) write(16,*) i, j, resolution(i,j)
      write(16,*) i, j, resolution(i,j)
    enddo
  enddo

  write(12,*)' pred         obsd        misfit      delta_d      imptnc'
  do i = 1, ndata
    pred(i) = obsd(i) - misfit(i)
    write(12,1001)pred(i), obsd(i), misfit(i), delta_d(i), dataimp(i)
1001 format(4(e11.4,1x),f10.6)
  enddo

  close(10)
  close(12)

end program 


!=======================================================================
! Forward and back substitution algorithm:
! Ax = b -> LUx = Pb and can be solve in two steps:
! (1), Ly = Pb and (2), Ux = y.
!
! Lx = b or Ux = b is very easy to solve by (1) forward substitution
! for lower triangular matrice and, (2) back substitution for upper
! triangular matrices.
!
! Forward substitution (Lx = b),
! l(11)*x(1)                                  = b(1)
! l(21)*x(1)  + l(22)*x(2)                    = b(2)
!   .         .                               .
! l(n1)*x(1)  + l(n2)*x(2) + ... + l(nn)*x(n) = b(n)
!
! soltion,
! x(1) = b(1)/l(11)
! x(2) = ( b(2)-l(21)*x(1) ) / l(22)
!      .
! x(n) = ( b(n)-(sum(l(ni)*x(i)),i=1,...,n-1)) / l(nn)
!=======================================================================
SUBROUTINE dlubksb(a,n,np,indx,b)
      implicit none
      integer :: n,np
      integer, dimension(1:n) :: indx
      double precision :: a(np,np),b(n)
      integer :: i,ii,j,ll
      double precision :: sum

      ii=0

      ! Froward solve Ly = b
      do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do j=ii,i-1
            sum=sum-a(i,j)*b(j)
          enddo
        else if (sum.ne.0.) then
          ii=i
        endif

        b(i)=sum
      enddo

      ! Backward solve Ux = y
      do i=n,1,-1
        sum=b(i)
        do j=i+1,n
          sum=sum-a(i,j)*b(j)
        enddo
        b(i)=sum/a(i,i)
      enddo

      return
END

!=======================================================================
! LU decomposition:
! A -> LU  lower/upper triangular maxtrices when A is not singular matrix
!
!     |l(11)  0 ......  0   |        |  0  u(12)... u(1n) |
!     |l(21) l(22) ...  0   |        |  0   0   ... u(2n) |
! L = | :     :      .  0   |    U = |  :   :    0    :   |
!     | :     :         0   |        |  :   :    .    :   |
!     |l(n1) l(n2) ... l(nn)|        |  0   0 ... ... 0   |
!
! normally, l(nn) = 1, n = 1,...,N
!
!=======================================================================
SUBROUTINE dludcmp(a,n,np,indx,d)
      implicit none
      integer, parameter :: NMAX=500
      double precision, parameter :: TINY=1.0e-20
      integer :: i,imax,j,k
      integer :: n,np
      integer, dimension(1:n) :: indx, vv
      double precision :: d
      double precision, dimension(1:np,1:np) :: a
      double precision :: aamax,dum,sum
      double precision, dimension(1:NMAX) :: evv

      d=1.

      do i=1,n
        aamax=0.
        do j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
        enddo
        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
      enddo

      do j=1,n
        do i=1,j-1
          sum=a(i,j)
          do k=1,i-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
        enddo

        aamax=0.

        do i=j,n
          sum=a(i,j)
          do k=1,j-1
            sum=sum-a(i,k)*a(k,j)
          enddo

          a(i,j)=sum
          dum=vv(i)*abs(sum)

          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
        enddo

        if (j.ne.imax)then
          do k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          enddo

          d=-d
          vv(imax)=vv(j)
        endif

        indx(j)=imax

        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do i=j+1,n
            a(i,j)=a(i,j)*dum
          enddo
        endif

      enddo
      return
END
