! ======================================================================== !
! ALOGES PROJECT                                                           !
! Quim Ballabrera, April 2022                                              !
! Institut de Ciencies del Mar, CSIC                                       !
! Last Modified: 2022-04-14                                                !
!                                                                          !
! Copyright (C) 2022, Joaquim Ballabrera                                   !
!                                                                          !
! This program is free software: you can redistribute it and/or modify     !
! it under the terms of the GNU Lesser General Public License as published !
! by the Free Software Foundation, either version 3 of the License, or     !
! (at your option) any later version.                                      !
!                                                                          !
! This program is distributed in the hope that it will be useful,          !
! but WITHOUT ANY WARRANTY; without even the implied warranty of           !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                     !
! See the Lesser GNU General Public License for more details.              !
!                                                                          !
! You should have received a copy of the GNU Lesser General                !
! Public License along with this program.                                  !
! If not, see <http://www.gnu.org/licenses/>.                              !
! - randn                                                                  !
! -------------------------------------------------------------------------!

module module_math

use, intrinsic :: IEEE_ARITHMETIC, ONLY : IEEE_VALUE, IEEE_QUIET_NAN
use module_types
use module_constants

implicit none

interface randn
  module procedure randn_r,randn_v,randn_a
end interface randn

interface arange
  module procedure arange_i
  module procedure arange_dp
end interface arange

contains
! ...
! =====================================================================
! =====================================================================
! ...
  function randn_r () result(ff)

    real(dp)                                     :: ff

    ! ... Local variables
    ! ...
    real(dp), parameter                          :: s  =  0.449871D0
    real(dp), parameter                          :: t  = -0.386595D0
    real(dp), parameter                          :: a  =  0.19600D0
    real(dp), parameter                          :: b  =  0.25472D0
    real(dp), parameter                          :: r1 =  0.27597D0
    real(dp), parameter                          :: r2 =  0.27846D0
    real(dp) u,v,x,y,q

    do
      call RANDOM_NUMBER(u)  ! GNU RANDOM GENERATOR
      call RANDOM_NUMBER(v)  ! GNU RANDOM GENERATOR
      v = 1.7156D0 * (v - 0.5D0)

      ! ... Evaluate the quadratic form
      ! ...
      x = u - s
      y = ABS(v) - t
      q = x*x + y*(a*y - b*x)

      if (q .lt. r1) exit
      if (q .gt. r2) cycle
      if (v**2 .LT. -4D0*LOG(u)*u*u) exit
    enddo
    ff = v/u

  end function randn_r
  ! ...
  ! ===================================================================
  ! ...
  function randn_v (m) result(ff)

    integer, intent(in)                          :: m
    real(dp), dimension(m)                       :: ff

    ! ... Local variables
    ! ...
    integer i
    real(dp), parameter                          :: s  =  0.449871D0
    real(dp), parameter                          :: t  = -0.386595D0
    real(dp), parameter                          :: a  =  0.19600D0
    real(dp), parameter                          :: b  =  0.25472D0
    real(dp), parameter                          :: r1 =  0.27597D0
    real(dp), parameter                          :: r2 =  0.27846D0
    real(dp) u,v,x,y,q
  
    do i=1,m
      do
        call RANDOM_NUMBER(u)  ! GNU RANDOM GENERATOR
        call RANDOM_NUMBER(v)  ! GNU RANDOM GENERATOR
        v = 1.7156D0 * (v - 0.5D0)
  
        ! ... Evaluate the quadratic form
        ! ...
        x = u - s
        y = ABS(v) - t
        q = x*x + y*(a*y - b*x)

        if (q .lt. r1) exit
        if (q .gt. r2) cycle
        if (v**2 .LT. -4D0*LOG(u)*u*u) exit
      enddo
      ff(i) = v/u
    enddo

  end function randn_v
  ! ...
  ! ===================================================================
  ! ...
  function randn_a (m,n) result(ff)

    integer, intent(in)                          :: m
    integer, intent(in)                          :: n
    real(dp), dimension(:,:), pointer            :: ff

    ! ... Local variables
    ! ...
    integer i,j
    real(dp), parameter                          :: s  =  0.449871D0
    real(dp), parameter                          :: t  = -0.386595D0
    real(dp), parameter                          :: a  =  0.19600D0
    real(dp), parameter                          :: b  =  0.25472D0
    real(dp), parameter                          :: r1 =  0.27597D0
    real(dp), parameter                          :: r2 =  0.27846D0
    real(dp) u,v,x,y,q

    if (.not.associated(ff)) allocate(ff(m,n))

    do j=1,n
    do i=1,m
      do
        call RANDOM_NUMBER(u)  ! GNU RANDOM GENERATOR
        call RANDOM_NUMBER(v)  ! GNU RANDOM GENERATOR
        v = 1.7156D0 * (v - 0.5D0)

        ! ... Evaluate the quadratic form
        ! ...
        x = u - s
        y = ABS(v) - t
        q = x*x + y*(a*y - b*x)

        if (q .lt. r1) exit
        if (q .gt. r2) cycle
        if (v**2 .LT. -4D0*LOG(u)*u*u) exit
      enddo
      ff(i,j) = v/u
    enddo
    enddo

  end function randn_a
  ! ...
  ! ===================================================================
  ! ...
  function rndname(len,iseed) result(name)

  integer, intent(in)                            :: len
  integer, optional                              :: iseed
  character(len=len)                             :: name
  
  ! ... Local variables
  ! ...
  integer i,io,il,j,n
  integer, dimension(:), allocatable             :: seed
  real(dp) r

  if (present(iseed)) then
    call random_seed(size=n)
    allocate(seed(n))
    seed(:) = iseed
    call random_seed(put=seed)
  endif

  io = ichar('A')
  il = ichar('Z') - io

  do i=1,len
    call random_number(r)
    j = int(io + il*r)
    name(i:i) = char(j)
  enddo

  return
  end function rndname
  ! ...
  ! ===================================================================
  ! ...
  function randseries(n,dlag,periodic) result(g)

  integer, intent(in)                            :: n
  integer, intent(in), optional                  :: dlag
  logical, intent(in), optional                  :: periodic
  real(dp), dimension(n)                         :: g
  
  ! ... Local variables
  ! ...
  logical lper
  integer nn,i,iter,im,ip
  real(dp) xsum
  real(dp), dimension(:), allocatable            :: v,f
  
  lper = .False.
  if (present(periodic)) then
    if (periodic) lper = .True.
  endif
  
  if (lper) then
  
    if (present(dlag)) then
      allocate(f(n))
  
      g = randn(n)
      g(n) = g(1)

      print*, dlag, g(n), g(1)
      do iter=1,6*abs(dlag)
        do i=1,n
          im = i - 1
          if (im.eq.0) im = n
          ip = i + 1
          if (ip.gt.n) ip = 1
          f(i) = 0.25D0*(g(im)+2*g(i)+g(ip))
        enddo
        g(:) = f(:)
      enddo
      deallocate(f)
      xsum = SUM(g)/n
      g(:) = g(:) - xsum
      xsum = sqrt(DOT_PRODUCT(g,g)/n)
      g(:) = g(:) / xsum
    else
      g = randn(n)
      g(n) = g(1)
    endif
  else
    if (present(dlag)) then
  
      ! ... Buffering time series
      ! ...
      nn = n + 20*abs(dlag)
      allocate(v(nn))
      allocate(f(nn))

      v = randn(nn)
      do iter=1,6*abs(dlag)
        f(1) = v(1)
        f(nn) = v(nn)
        do i=2,nn-1
          f(i) = 0.25D0*(v(i-1)+2*v(i)+v(i+1))
        enddo
        v = f
      enddo

      ! ... Extracting far from boundary
      ! ...
      i = 10*abs(dlag)
      g(:) = f(i+1:i+n)

      xsum = SUM(g)/n
      g(:) = g(:) - xsum
      xsum = sqrt(DOT_PRODUCT(g,g)/n)
      g(:) = g(:) / xsum

      deallocate(f)
      deallocate(v)
  
    else
  
      g = randn(n)
  
    endif
  endif

  end function randseries
  ! ...
  ! ===================================================================
  ! ...
  real(dp) function haversine (lon1,lat1,lon2,lat2)
    ! ...
    ! ... Function Haversine
    ! ... Determines the great-circle distance between two points in a
    ! ... sphere. The input lngitudes and latitudes are given in radians
    ! ... The retruned distance is in meters.
    ! ... REarth = 6371315.0_dp      ! m
    ! ...
    real(dp), intent(in)                  :: lon1,lat1
    real(dp), intent(in)                  :: lon2,lat2

    ! ... Local variables
    ! ...
    real(dp) sindx,sindy,dang

    sindx = sin(0.5D0*(lon2-lon1))
    sindy = sin(0.5D0*(lat2-lat1))

    dang = 2.0d0*asin(sqrt(sindy*sindy + cos(lat2)*cos(lat1)*sindx*sindx))
    haversine = REarth * dang

    return
  end function haversine
  ! ...
  ! ===================================================================
  ! ...
  pure function arange_i(amin, amax, adelta) result(arange)
    ! ...
    ! ... Returns a double precision vector from amin to amax with 
    ! ... an optional step of adelta.
    ! ...

    integer, intent(in)                  :: amin 
    integer, intent(in)                  :: amax 
    integer, intent(in), optional        :: adelta 
    integer, dimension(:), allocatable   :: arange

    ! ... Loal variables
    ! ...
    integer  i, nlen
    integer delta

    if (present(adelta))then
      delta = adelta
    else
      delta = 1
    endif

    nlen = (amax-amin+0.5*delta)/delta + 1
    allocate(arange(nlen))

    do i=1,nlen
      arange(i) = amin + (i-1)*delta
    enddo

    return
  end function arange_i
  ! ...
  ! ===================================================================
  ! ...
  pure function arange_dp(amin, amax, adelta) result(arange)
    ! ...
    ! ... Returns a double precision vector from amin to amax with 
    ! ... an optional step of adelta.
    ! ...

    real(dp), intent(in)                  :: amin 
    real(dp), intent(in)                  :: amax 
    real(dp), intent(in), optional        :: adelta 
    real(dp), dimension(:), allocatable   :: arange

    ! ... Loal variables
    ! ...
    integer  i, nlen
    real(dp) delta

    if (present(adelta))then
      delta = adelta
    else
      delta = 1
    endif

    nlen = (amax-amin+0.5*delta)/delta + 1
    allocate(arange(nlen))

    do i=1,nlen
      arange(i) = amin + (i-1.0D0)*delta
    enddo

    return
  end function arange_dp
  ! ...
  ! ===================================================================
  ! ...
  real(dp) function mean(A,W)
    ! ... Calculates the weighted mean = 1/N * Sum W(i)*A(i)
    ! ... Weights are optional.

    real(dp), dimension(:), intent(in)     :: A
    real(dp), dimension(:), optional       :: W

    integer n
    real(dp) nan,Sw

    ! ... nan:
    ! ...
    nan = ieee_value(1.0_dp,ieee_quiet_nan)

    mean = nan
    n = size(A)
    if (n.eq.0) return

    if (present(W)) then
      Sw   = sum(W)
      mean = dot_product(W,A)/Sw
    else
      mean = sum(A)/N
    endif

  end function mean
  ! ...
  ! ===================================================================
  ! ...
  real(dp) function variance (A,W)

    ! ... Calculates the variance = FACTOR * Sum (A(i)-mean(A))**2
    ! ... The value of FACTOR depends of W. By default, W=0, FACTOR = 1/(N-1)
    ! ... However, if W=1, the biased variance is calculated, i.e. FACTOR=1/N
    ! ... If W has to be used as a weight it must be a vector with the same
    ! ... size than A.
    ! ...
    ! ... Weights are optional.
    ! ...

    real(dp), dimension(:), intent(in)     :: A
    real(dp), dimension(:), optional       :: W

    logicaL                                    :: weight = .false.
    logical                                    :: biased = .false.
    integer N,i
    real(dp) xsum1,xsum2,ai,wi,Sw
    real(dp) nan


    ! ... nan:
    ! ...
    nan = ieee_value(1.0_dp,ieee_quiet_nan)

    variance = nan
    N = SIZE(A)
    if (N.EQ.0) return

    if (PRESENT(W)) then
      if (SIZE(W).EQ.1) then
        if (W(1).EQ.0) then
          biased = .false.
        else if (W(1).EQ.1) then
          biased = .true.
        else
          STOP 'ERROR in variance: Invalid normalization. Use 0 (unbiased) or 1 (biased)'
        endif
      else if (SIZE(W).EQ.N) then
        weight = .true.
      else
        STOP 'ERROR in variance: Size W must be 1 or N'
      endif
    endif

    if (weight) then
      Sw    = 0D0
      xsum1 = 0D0
      xsum2 = 0D0
      do i=1,N
        wi = w(i)
        ai = A(i)
        Sw    = Sw    + wi
        xsum1 = xsum1 + wi*ai
        xsum2 = xsum2 + wi*ai*ai
      enddo
      xsum1 = xsum1 / Sw
      xsum2 = xsum2/Sw-xsum1*xsum1
      variance = Sw * xsum2 /(Sw - 1D0)
    else
      xsum1 = 0D0
      xsum2 = 0D0
      do i=1,N
        ai = A(i)
        xsum1 = xsum1 + ai
        xsum2 = xsum2 + ai*ai
      enddo
      xsum1 = xsum1 / N
      xsum2 = xsum2 / N - xsum1*xsum1
      IF (biased) then
        variance = xsum2
      else
        variance = N*xsum2/(N-1D0)
      endif
    endif
    IF (variance.LT.0) variance = 0D0

  end function variance
  ! ...
  ! ===================================================================
  ! ...
end module module_math
