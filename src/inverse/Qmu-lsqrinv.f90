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
!=======================================================================
module  const
  implicit none
  real, parameter :: pi  = 3.14159265359
  real, parameter :: eps = 1.0E-05
  real, parameter :: ZERO = 1.0E-30
  integer, parameter :: ndat = 14       ! number of data (14 period band )
  integer, parameter :: nvar = 21       ! number of model variables (21 layers)
end module

!===
!   Linear inversion:
!     G*m = d  
!       m = (G^T*G)^(-1)*G^T*d = gtginv * gtd
!===
program lininv
  use const
  implicit none
  integer :: i, j, k, icount,n

  real :: dummy1
  integer :: eof
  integer :: ndata    ! number of used data 
  real, dimension(1:nvar) :: md_00, md, change, mdo
  real, dimension(1:nvar) :: varm
  real, dimension(1:ndat) :: obsd, stdv, pred

  double precision :: b, dpx
  double precision :: d_variance, dstdv
  double precision :: sumsq
  double precision :: dt1, dt2

  integer, dimension(1:nvar) :: indx
  double precision, dimension(1:ndat) :: delta_d
  double precision, dimension(1:ndat) :: misfit
  double precision, dimension(1:ndat) :: dataimp
  double precision, dimension(1:nvar) :: tempimp
  double precision, dimension(1:nvar) :: stdv_md
  double precision, dimension(1:nvar) :: covinv  ! Cmm, damping, cov of md stdv.
  double precision, dimension(1:ndat,1:nvar) :: g
  double precision, dimension(1:nvar) :: gtd, gtdcmm           ! G^T*d
  double precision, dimension(1:nvar,1:nvar) :: gtg, gtginv   ! G^T*G
  double precision, dimension(1:nvar,1:nvar) :: fcov, corr

  character(len=10) :: dummy2, dummy3
  character(len=10), dimension(1:nvar) :: label, label2
  character(len=80) :: f_obs
  character(len=80) :: f_d
  character(len=80) :: f_m
  character(len=80) :: f_m00
  character(len=80) :: f_g
  character(len=80) :: foutput
  character(len=80) :: ftsfpred

  !===
  !   Read input files
  !===
  write(*,*)'Input file of obs and stdv:'
  read(*,*)f_obs
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

  write(*,*)'file name for parameter output:'
  read(*,*)foutput
  write(*,*)'file name is ',foutput
  open(10,file=foutput,status='new')

  write(*,*)'file name for best fit transf:'
  read(*,*)ftsfpred
  write(*,*)'file name is ',ftsfpred
  open(12,file=ftsfpred,status='new')

  !===
  !   Readin g data matrix:
  !         delta_d = log10(syn) - log10(obs)
  !   and the G(kernel) matrix
  !         g(i,j) = d(syn1(i)-syn2(i))/d(md1(j)-md2(j))
  !   Both are normalized by dividing the standard deviation
  !===

  open(13,file=f_obs,status='old')
  i = 0
  do
    read(13,*,iostat=eof)n, obsd(n), stdv(n)
    if(eof /= 0) exit
    i = i + 1
  enddo
  if(i /= ndat) then
    write(*,*)'number of data is not correct!'
    stop
  endif
  close(13)

  open(13,file=f_d,status='old')
  i = 0
  do
    read(13,*,iostat=eof)n, delta_d(n)
    if(eof /= 0) exit
    i = i + 1
  enddo
  if(i /= ndat) then
    write(*,*)'number of data is not correct!'
    stop
  endif
  close(13)

  open(13,file=f_m,status='old')
  i = 0
  do
    read(13,*,iostat=eof)n, md(n), varm(n), label(n), label2(n)
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
    read(13, *, iostat=eof)n, md_00(n), dummy1, dummy2, dummy3
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
    if(eof /= 0) exit
  enddo
  close(13)

  !== test ==
  !stop
  !== end ==

  ndata = ndat

  do i = 1, nvar
    covinv(i) = 1/(varm(i)**2) 
  enddo

  !===
  !   Calculate GT*G and GT*d
  !===
  do j = 1, nvar
    gtd(j) = 0.0

    do i = 1, ndata
      gtd(j) = gtd(j) + g(i,j) * delta_d(i)
    enddo

    gtdcmm(j) = gtd(j) - covinv(j) * (md(j) - md_00(j))

    do k = 1, j
      gtg(k,j) = 0.0
      do i = 1, ndata
        gtg(k,j) = gtg(k,j) + g(i,k) * g(i,j)
      enddo
      gtg(j,k) = gtg(k,j)
      !write(*,*)gtg(j,k),j,k
    enddo
    gtg(j,j) = gtg(j,j) + covinv(j)
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
  !   Find normalized residuals (linear problem with zero start
  !   so just use partial derivatives) and sum of squares of errors
  !===
  sumsq = 0.0
  do i = 1, ndata
    misfit(i) = 0.0

    do j = 1, nvar
      misfit(i) = misfit(i) + g(i,j) * change(j)
    enddo

    misfit(i) = delta_d(i) - misfit(i)
    sumsq = sumsq + misfit(i)**2
  enddo
  d_variance = sumsq/(ndata-nvar) ! nfree = vs, coef, thickness
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
    write(10,998)label(i),label2(i),mdo(i),stdv_md(i),md(i)
  enddo
998 format(A10,A10,1x,e12.5,1x,e12.5,1x,e12.5)


  write(12,*)'Description'
  write(12,*)
  write(12,*)'Normalized standard deviation data',dstdv
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

  write(12,*)'covariance matrix'
  do i = 1, nvar
    do j = 1, nvar 
      write(12,*)fcov(i,j)
    enddo
  enddo

  write(12,*)
  write(12,*)'correlation matrix'
  do i = 1, nvar
    do j = 1, nvar
      write(12,*)corr(i,j)
    enddo
  enddo

  write(12,*)' pred        obsd        misfit    delta_d    imptnc'
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
