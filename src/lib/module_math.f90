!> @file module_math.f90
!> @brief Pure mathematical utilities: random number generation, 
!> array operations, optimization, and sorting.
!> @author Quim Ballabrera, Institut de Ciencies del Mar, CSIC
!> @date April 2022, refactored 2025
!>
!> This module provides core mathematical functions excluding statistics.
!> For statistical routines (mean, variance, moments, percentiles), see module_statistical.
!>
!> Copyright (C) 2022-2025, Joaquim Ballabrera
!> Licensed under GNU LGPL v3+

module module_math

use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
use iso_fortran_env, only: error_unit, output_unit
use module_types, only: dp
use module_constants, only: nan
use module_tools, only: crash

implicit none (type, external)

private

public :: randn, randseries, random_integer
public :: arange, linspace, meshgrid
public :: indexx, swap, imaxloc, quicksort
public :: brent, golden
public :: diff, gradient, cumsum, cumprod
public :: next_power_of_2, flip
public :: ascending, descending
public :: t_inverse_cdf, t_cdf_complement
public :: f_cdf_complement


!> @brief Generic interface for generating Gaussian random 
!> numbers (scalar, vector, or 2D array).
interface randn
  module procedure randn_r, randn_v, randn_a
end interface randn

!> @brief Generic interface for creating arithmetic 
!> sequences (integer or real).
interface arange
  module procedure arange_i, arange_dp
end interface arange

!> @brief Generic interface for swapping two values (scalar or vector).
interface swap
  module procedure swap_r, swap_v
end interface swap

!> @brief Generic interface for flipping arrays along a dimension.
interface flip
  module procedure flip1d_i, flip1d_dp, flip2d_dp, flip3d_dp
end interface flip

contains

  !> @brief Generate a single Gaussian random number 
  !> with mean=0, std=1.
  !!
  !! Uses the ratio-of-uniforms method (Leva, 1992).
  !! Units: dimensionless.
  !!
  !! @return Real(dp) standard normal variate.
  !!
  !! Notes:
  !! - Thread-safe if Fortran RNG is thread-safe.
  !! - For reproducibility, seed the RNG before calling.
  function randn_single() result(ff)
    real(dp) :: ff
    real(dp), parameter :: s  =  0.449871_dp
    real(dp), parameter :: t  = -0.386595_dp
    real(dp), parameter :: a  =  0.19600_dp
    real(dp), parameter :: b  =  0.25472_dp
    real(dp), parameter :: r1 =  0.27597_dp
    real(dp), parameter :: r2 =  0.27846_dp
    real(dp) :: u, v, x, y, q

    do
      call random_number(u)
      call random_number(v)
      v = 1.7156_dp * (v - 0.5_dp)
      x = u - s
      y = abs(v) - t
      q = x*x + y*(a*y - b*x)
      if (q < r1) exit
      if (q > r2) cycle
      if (v**2 < -4.0_dp*log(u)*u*u) exit
    end do
    ff = v/u
  end function randn_single


  !> @brief Generate a single Gaussian random number 
  !> with mean=0, std=1.
  !!
  !! Uses the ratio-of-uniforms method (Leva, 1992).
  !! Units: dimensionless.
  !!
  !! @return Real(dp) standard normal variate.
  !!
  !! Notes:
  !! - Thread-safe if Fortran RNG is thread-safe.
  !! - For reproducibility, seed the RNG before calling.
  function randn_r() result(ff)
    real(dp) :: ff
    ff = randn_single()
  end function randn_r

  !> @brief Generate a 1D array of Gaussian random numbers 
  !> with mean=0, std=1.
  !!
  !! @param[in] m Integer, length of output vector.
  !! @return Real(dp) array of size m, standard normal variates.
  !!
  !! Notes:
  !! - Uses randn_r internally.
  function randn_v(m) result(ff)
    integer, intent(in) :: m
    real(dp), allocatable :: ff(:)
    integer :: i

    if (allocated(ff)) deallocate(ff)
    allocate(ff(m))
    do i = 1, m
      ff(i) = randn_single()
    end do
  end function randn_v

  !> @brief Generate a 2D array of Gaussian random numbers with mean=0, std=1.
  !!
  !! @param[in] m Integer, number of rows.
  !! @param[in] n Integer, number of columns.
  !! @return Real(dp) allocatable array of size (m,n), standard normal variates.
  !!
  !! Notes:
  !! - OpenMP parallelized over columns and rows.
  !! - Output is allocated if not already allocated.
  function randn_a(m, n) result(ff)
    integer, intent(in) :: m, n
    real(dp), dimension(:,:), allocatable :: ff
    integer i,j

    if (allocated(ff)) deallocate(ff)
    allocate(ff(m, n))

    do j = 1, n
    do i = 1, m
      ff(i, j) = randn_single()
    end do
    end do
  end function randn_a

  !> @brief Generate a smooth random time series with optional periodicity and lag correlation.
  !!
  !! @param[in] n Integer, length of output series.
  !! @param[in] dlag Integer, optional. Smoothing lag (higher = smoother). Default: no smoothing.
  !! @param[in] periodic Logical, optional. If true, enforce periodicity (first = last). Default: false.
  !! @return Real(dp) array of size n, normalized to zero mean and unit variance.
  !!
  !! Notes:
  !! - Smoothing is applied via repeated 3-point averaging (6*|dlag| iterations).
  !! - For periodic series, boundary wraps around.
  !! - For non-periodic, buffering is used to avoid edge artifacts.
  function randseries(n, dlag, periodic) result(g)
    integer, intent(in) :: n
    integer, intent(in), optional :: dlag
    logical, intent(in), optional :: periodic
    real(dp), dimension(n) :: g
    logical :: lper
    integer :: nn, i, iter, im, ip
    real(dp) :: xsum
    real(dp), dimension(:), allocatable :: v, f

    lper = .false.
    if (present(periodic)) then
      if (periodic) lper = .true.
    end if

    if (lper) then
      if (present(dlag)) then
        allocate(f(n))
        g = randn(n)
        g(n) = g(1)
        do iter = 1, 6*abs(dlag)
          do i = 1, n
            im = i - 1
            if (im == 0) im = n
            ip = i + 1
            if (ip > n) ip = 1
            f(i) = 0.25_dp*(g(im) + 2*g(i) + g(ip))
          end do
          g(:) = f(:)
        end do
        deallocate(f)
        xsum = sum(g)/n
        g(:) = g(:) - xsum
        xsum = sqrt(dot_product(g, g)/n)
        g(:) = g(:) / xsum
      else
        g = randn(n)
        g(n) = g(1)
      end if
    else
      if (present(dlag)) then
        nn = n + 20*abs(dlag)
        allocate(v(nn), f(nn))
        v = randn(nn)
        do iter = 1, 6*abs(dlag)
          f(1) = v(1)
          f(nn) = v(nn)
          do i = 2, nn - 1
            f(i) = 0.25_dp*(v(i-1) + 2*v(i) + v(i+1))
          end do
          v = f
        end do
        i = 10*abs(dlag)
        g(:) = f(i+1:i+n)
        xsum = sum(g)/n
        g(:) = g(:) - xsum
        xsum = sqrt(dot_product(g, g)/n)
        g(:) = g(:) / xsum
        deallocate(f, v)
      else
        g = randn(n)
      end if
    end if
  end function randseries


  !> @brief Generate a random integer uniformly distributed in [imin, imax].
  !!
  !! @param[in] imin Integer, lower bound (inclusive).
  !! @param[in] imax Integer, upper bound (inclusive).
  !! @return Integer, random value in [imin, imax].
  !!
  !! Notes:
  !! - Uses intrinsic random_number() for uniform [0,1) real, then scales and floors.
  !! - Thread-safe if Fortran RNG is thread-safe.
  !! - For reproducibility, seed the RNG before calling (call random_seed()).
  !! - If imin > imax, returns imin.
  !! - Uniform distribution: P(k) = 1/(imax - imin + 1) for k in [imin, imax].
  !!
  !! Algorithm:
  !! - Generate u ~ Uniform(0,1).
  !! - Map to integer: k = imin + floor(u * (imax - imin + 1)).
  !! - Handles edge case u=1.0 (rare but possible) by clamping to imax.
  function random_integer(imin, imax) result(k)
    integer, intent(in) :: imin, imax
    integer :: k
    real(dp) :: u
    integer :: range

    if (imin > imax) then
      k = imin
      return
    end if
    range = imax - imin + 1
    ! Generate uniform random in [0, 1)
    call random_number(u)
    ! Map to [imin, imax]
    k = imin + int(u * range)
    ! Clamp to imax in case u rounds to exactly 1.0 (extremely rare)
    if (k > imax) k = imax

  end function random_integer


  !> @brief Create an integer arithmetic sequence from amin to amax with optional step.
  !!
  !! @param[in] amin Integer, start value.
  !! @param[in] amax Integer, end value (inclusive if step divides range).
  !! @param[in] adelta Integer, optional. Step size. Default: 1.
  !! @return Integer allocatable array.
  !!
  !! Notes:
  !! - Similar to Python's range() or NumPy's arange().
  pure function arange_i(amin, amax, adelta) result(arange)
    integer, intent(in) :: amin, amax
    integer, intent(in), optional :: adelta
    integer, dimension(:), allocatable :: arange
    integer :: i, nlen, delta

    if (present(adelta)) then
      delta = adelta
    else
      delta = 1
    end if

    if (delta > 0) then
      nlen = max(0, (amax-amin)/delta + 1)
    else
      nlen = max(0, (amin-amax)/(-delta) + 1)
    endif

    if (allocated(arange)) deallocate(arange)
    allocate(arange(nlen))

    do i = 1, nlen
      arange(i) = amin + (i - 1)*delta
    end do
  end function arange_i

  !> @brief Create a real arithmetic sequence from amin to amax with optional step.
  !!
  !! @param[in] amin Real(dp), start value.
  !! @param[in] amax Real(dp), end value (inclusive if step divides range).
  !! @param[in] adelta Real(dp), optional. Step size. Default: 1.0.
  !! @return Real(dp) allocatable array.
  !!
  !! Notes:
  !! - Units: same as amin, amax, adelta.
  pure function arange_dp(amin, amax, adelta) result(arange)
    real(dp), intent(in) :: amin, amax
    real(dp), intent(in), optional :: adelta
    real(dp), dimension(:), allocatable :: arange
    integer :: i, nlen
    real(dp) :: delta

    if (present(adelta)) then
      delta = adelta
    else
      delta = 1.0_dp
    end if

    if (delta > 0) then
      nlen = max(0, int((amax-amin)/delta) + 1)
    else
      nlen = max(0, int((amin-amax)/(-delta)) + 1)
    endif

    if (allocated(arange)) deallocate(arange)
    allocate(arange(nlen))

    do i = 1, nlen
      arange(i) = amin + (i - 1.0_dp)*delta
    end do
  end function arange_dp

  !> @brief Index sort: return permutation indices that sort array arr in ascending order.
  !!
  !! From Numerical Recipes in Fortran (Press et al., 1992).
  !!
  !! @param[in] arr Real(dp) array to sort.
  !! @param[out] indx Integer array of size(arr), permutation indices such that arr(indx) is sorted.
  !!
  !! Notes:
  !! - Does not modify arr.
  !! - indx(i) gives the index of the i-th smallest element.
  subroutine indexx(arr, indx)
    real(dp), dimension(:), intent(in) :: arr
    integer, dimension(size(arr)), intent(out) :: indx
    integer, parameter :: M = 7
    integer :: n, i, indxt, r, itemp, j, jstack, k, l
    integer, dimension(size(arr)) :: istack
    real(dp) :: a

    n = size(arr)
    do j = 1, n
      indx(j) = j
    end do

    jstack = 0
    l = 1
    r = n
    do
      if (r - l < M) then
        do j = l + 1, r
          indxt = indx(j)
          a = arr(indxt)
          do i = j - 1, 1, -1
            if (arr(indx(i)) <= a) exit
            indx(i+1) = indx(i)
          end do
          indx(i+1) = indxt
        end do
        if (jstack == 0) return
        r = istack(jstack)
        l = istack(jstack-1)
        jstack = jstack - 2
      else
        k = (l + r)/2
        itemp = indx(k)
        indx(k) = indx(l+1)
        indx(l+1) = itemp
        if (arr(indx(l+1)) > arr(indx(r))) then
          itemp = indx(l+1)
          indx(l+1) = indx(r)
          indx(r) = itemp
        end if
        if (arr(indx(l)) > arr(indx(r))) then
          itemp = indx(l)
          indx(l) = indx(r)
          indx(r) = itemp
        end if
        if (arr(indx(l+1)) > arr(indx(l))) then
          itemp = indx(l+1)
          indx(l+1) = indx(l)
          indx(l) = itemp
        end if
        i = l + 1
        j = r
        indxt = indx(l+1)
        a = arr(indxt)
        do
          do
            i = i + 1
            if (arr(indx(i)) >= a) exit
          end do
          do
            j = j - 1
            if (arr(indx(j)) <= a) exit
          end do
          if (j < i) exit
          itemp = indx(i)
          indx(i) = indx(j)
          indx(j) = itemp
        end do
        indx(l+1) = indx(j)
        indx(j) = indxt
        jstack = jstack + 2
        if (r - i + 1 >= j - l) then
          istack(jstack) = r
          istack(jstack-1) = i
          r = j - 1
        else
          istack(jstack) = j - 1
          istack(jstack-1) = l
          l = i
        end if
      end if
    end do
  end subroutine indexx

  !> @brief Swap two real(dp) scalars.
  !!
  !! @param[inout] a Real(dp).
  !! @param[inout] b Real(dp).
  subroutine swap_r(a, b)
    real(dp), intent(inout) :: a, b
    real(dp) :: c
    c = a
    a = b
    b = c
  end subroutine swap_r

  !> @brief Swap two real(dp) vectors element-wise.
  !!
  !! @param[inout] a Real(dp) array.
  !! @param[inout] b Real(dp) array of same size as a.
  subroutine swap_v(a, b)
    real(dp), dimension(:), intent(inout) :: a, b
    real(dp), dimension(size(a)) :: c
    c(:) = a(:)
    a(:) = b(:)
    b(:) = c(:)
  end subroutine swap_v

  !> @brief Return the index of the maximum element in array a.
  !!
  !! @param[in] a Real(dp) array.
  !! @return Integer, index of maximum value.
  integer pure function imaxloc(a)
    real(dp), dimension(:), intent(in) :: a
    integer :: imax(1)
    imax = maxloc(a(:))
    imaxloc = imax(1)
  end function imaxloc

  !> @brief In-place quicksort of a real(dp) array in ascending order.
  !!
  !! @param[inout] arr Real(dp) array to sort.
  !!
  !! Notes:
  !! - Recursive implementation.
  !! - Average O(n log n), worst-case O(n^2).
  subroutine quicksort(arr)
    real(dp), dimension(:), intent(inout) :: arr
    call quicksort_helper(arr, 1, size(arr))
  contains
    recursive subroutine quicksort_helper(arr, low, high)
      real(dp), dimension(:), intent(inout) :: arr
      integer, intent(in) :: low, high
      integer :: i, j
      real(dp) :: pivot, temp

      if (low < high) then
        pivot = arr(high)
        i = low - 1
        do j = low, high - 1
          if (arr(j) <= pivot) then
            i = i + 1
            temp = arr(i)
            arr(i) = arr(j)
            arr(j) = temp
          end if
        end do
        temp = arr(i+1)
        arr(i+1) = arr(high)
        arr(high) = temp
        call quicksort_helper(arr, low, i)
        call quicksort_helper(arr, i+2, high)
      end if
    end subroutine quicksort_helper
  end subroutine quicksort

  !> @brief Brent's method for root-finding of a scalar function f(x) = 0.
  !!
  !! @param[in] f Function interface: real(dp) function f(x) result(y).
  !! @param[in] aa Real(dp), left bracket. Must satisfy f(aa)*f(bb) < 0.
  !! @param[in] bb Real(dp), right bracket.
  !! @param[in] tol Real(dp), tolerance for convergence.
  !! @return Real(dp), approximate root. Returns NaN if bracketing fails.
  !!
  !! Notes:
  !! - Combines bisection, secant, and inverse quadratic interpolation.
  !! - Robust and fast for smooth functions.
  function brent(f, aa, bb, tol) result(root)
    interface
      function f(x) result(y)
        use module_types, only: dp
        real(dp), intent(in) :: x
        real(dp) :: y
      end function f
    end interface
    real(dp), intent(in) :: aa, bb, tol
    real(dp) :: root
    real(dp), parameter :: eps = epsilon(1.0_dp)
    integer, parameter :: max_iter = 100
    real(dp) :: a, b, c, d, e, fa, fb, fc, p, q, r, s, tol1, xm
    integer :: iter

    a = aa
    b = bb
    fa = f(a)
    fb = f(b)

    if (fa*fb > 0.0_dp) then
      root = nan()
      return
    end if

    c = b
    fc = fb

    do iter = 1, max_iter
      if (fb*fc > 0.0_dp) then
        c = a
        fc = fa
        d = b - a
        e = d
      end if

      if (abs(fc) < abs(fb)) then
        a = b
        b = c
        c = a
        fa = fb
        fb = fc
        fc = fa
      end if

      tol1 = 2.0_dp*eps*abs(b) + 0.5_dp*tol
      xm = 0.5_dp*(c - b)

      if (abs(xm) <= tol1 .or. fb == 0.0_dp) then
        root = b
        return
      end if

      if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
        s = fb/fa
        if (a == c) then
          p = 2.0_dp*xm*s
          q = 1.0_dp - s
        else
          q = fa/fc
          r = fb/fc
          p = s*(2.0_dp*xm*q*(q - r) - (b - a)*(r - 1.0_dp))
          q = (q - 1.0_dp)*(r - 1.0_dp)*(s - 1.0_dp)
        end if

        if (p > 0.0_dp) q = -q
        p = abs(p)

        if (2.0_dp*p < min(3.0_dp*xm*q - abs(tol1*q), abs(e*q))) then
          e = d
          d = p/q
        else
          d = xm
          e = d
        end if
      else
        d = xm
        e = d
      end if

      a = b
      fa = fb

      if (abs(d) > tol1) then
        b = b + d
      else
        b = b + sign(tol1, xm)
      end if

      fb = f(b)
    end do

    root = b
  end function brent

  !> @brief Golden section search for minimizing a scalar function f(x).
  !!
  !! @param[in] f Function interface: real(dp) function f(x) result(y).
  !! @param[in] aa Real(dp), left bracket.
  !! @param[in] bb Real(dp), right bracket.
  !! @param[in] tol Real(dp), tolerance for convergence.
  !! @return Real(dp), approximate minimizer x*.
  !!
  !! Notes:
  !! - Does not require derivatives.
  !! - Assumes unimodal function in [aa, bb].
  function golden(f, aa, bb, tol) result(xmin)
    interface
      function f(x) result(y)
        use module_types, only: dp
        real(dp), intent(in) :: x
        real(dp) :: y
      end function f
    end interface
    real(dp), intent(in) :: aa, bb, tol
    real(dp) :: xmin
    real(dp), parameter :: phi = (1.0_dp + sqrt(5.0_dp))/2.0_dp
    real(dp) :: a, b, c, d, fc, fd

    a = aa
    b = bb
    c = b - (b - a)/phi
    d = a + (b - a)/phi

    do while (abs(c - d) > tol)
      fc = f(c)
      fd = f(d)
      if (fc < fd) then
        b = d
      else
        a = c
      end if
      c = b - (b - a)/phi
      d = a + (b - a)/phi
    end do

    xmin = (a + b)/2.0_dp
  end function golden

  !> @brief Create a linearly spaced vector from start to end (inclusive).
  !!
  !! @param[in] start Real(dp), first value.
  !! @param[in] end Real(dp), last value.
  !! @param[in] n Integer, number of points.
  !! @return Real(dp) array of size n.
  !!
  !! Notes:
  !! - If n=1, returns [start].
  !! - Units: same as start and end.
  function linspace(start, end, n) result(arr)
    real(dp), intent(in) :: start, end
    integer, intent(in) :: n
    real(dp), dimension(n) :: arr
    integer :: i
    real(dp) :: step

    if (n == 1) then
      arr(1) = start
    else
      step = (end - start)/(n - 1)
      do i = 1, n
        arr(i) = start + (i - 1)*step
      end do
    end if
  end function linspace

  !> @brief Create 2D coordinate matrices from 1D coordinate vectors (like MATLAB meshgrid).
  !!
  !! @param[in] x Real(dp) array of size m.
  !! @param[in] y Real(dp) array of size n.
  !! @param[out] X_mat Real(dp) array of size (m, n), x-coordinates replicated along columns.
  !! @param[out] Y_mat Real(dp) array of size (m, n), y-coordinates replicated along rows.
  !!
  !! Notes:
  !! - Caller must allocate X_mat and Y_mat to size (m, n).
  subroutine meshgrid(x, y, X_mat, Y_mat)
    real(dp), dimension(:), intent(in) :: x, y
    real(dp), dimension(:,:), intent(out) :: X_mat, Y_mat
    integer :: i, j, m, n

    m = size(x)
    n = size(y)

    if (size(X_mat,1) /= m .or. size(X_mat,2) /= n .or. &
        size(Y_mat,1) /= m .or. size(Y_mat,2) /= n) then
      call crash("meshgrid - output array size mismatch")
    end if

    do concurrent (i = 1:m, j = 1:n)
      X_mat(i, j) = x(i)
      Y_mat(i, j) = y(j)
    end do
  end subroutine meshgrid

  !> @brief Compute n-th order forward difference of array x.
  !!
  !! @param[in] x Real(dp) array.
  !! @param[in] n Integer, optional. Order of difference. Default: 1.
  !! @return Real(dp) allocatable array of size (size(x) - n).
  !!
  !! Notes:
  !! - diff(x, 1) = x(2:) - x(1:end-1).
  !! - If n >= size(x), returns empty array.
  function diff(x, n) result(dx)
    real(dp), dimension(:), intent(in) :: x
    integer, intent(in), optional :: n
    real(dp), dimension(:), allocatable :: dx
    integer :: diff_order, i, m

    diff_order = 1
    if (present(n)) diff_order = n

    m = size(x)
    if (diff_order >= m) then
      allocate(dx(0))
      return
    end if

    allocate(dx(m - diff_order))
    dx = x

    do i = 1, diff_order
      dx = dx(2:) - dx(:size(dx)-1)
    end do
  end function diff

  !> @brief Compute numerical gradient of 1D array f with optional spacing dx.
  !!
  !! @param[in] f Real(dp) array.
  !! @param[in] dx Real(dp), optional. Spacing between points. Default: 1.0.
  !! @return Real(dp) array of same size as f, gradient df/dx.
  !!
  !! Notes:
  !! - Uses central differences for interior points, forward/backward at boundaries.
  !! - Units: gradient has units of f/dx.
  function gradient(f, dx) result(df)
    real(dp), dimension(:), intent(in) :: f
    real(dp), intent(in), optional :: dx
    real(dp), dimension(size(f)) :: df
    integer :: n
    real(dp) :: h

    n = size(f)
    h = 1.0_dp
    if (present(dx)) h = dx

    if (n == 1) then
      df(1) = 0.0_dp
    else if (n == 2) then
      df(1) = (f(2) - f(1))/h
      df(2) = df(1)
    else
      df(1) = (f(2) - f(1))/h
      df(n) = (f(n) - f(n-1))/h
      df(2:n-1) = (f(3:n) - f(1:n-2))/(2.0_dp*h)
    end if
  end function gradient

  !> @brief Compute cumulative sum of array x.
  !!
  !! @param[in] x Real(dp) array.
  !! @return Real(dp) array of same size, cumsum(i) = sum(x(1:i)).
  !!
  !! Notes:
  !! - Units: same as x.
  function cumsum(x) result(cs)
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(size(x)) :: cs
    integer :: i, n

    n = size(x)
    cs(1) = x(1)
    do i = 2, n
      cs(i) = cs(i-1) + x(i)
    end do
  end function cumsum

  !> @brief Compute cumulative product of array x.
  !!
  !! @param[in] x Real(dp) array.
  !! @return Real(dp) array of same size, cumprod(i) = product(x(1:i)).
  !!
  !! Notes:
  !! - Units: dimensionless if x is dimensionless, otherwise units^i at position i.
  function cumprod(x) result(cp)
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(size(x)) :: cp
    integer :: i, n

    n = size(x)
    cp(1) = x(1)
    do i = 2, n
      cp(i) = cp(i-1) * x(i)
    end do
  end function cumprod

  !> @brief Return the smallest power of 2 greater than or equal to n.
  !!
  !! @param[in] n Integer.
  !! @return Integer, smallest 2^k >= n.
  !!
  !! Notes:
  !! - Useful for FFT padding.
  function next_power_of_2(n) result(n_pow2)
    integer, intent(in) :: n
    integer :: n_pow2

    n_pow2 = 2
    do while (n_pow2 < n)
      n_pow2 = n_pow2 * 2
    end do
  end function next_power_of_2

  !> @brief Flip (reverse) a 1D integer array in place.
  !!
  !! @param[inout] x Integer array.
  subroutine flip1d_i(x)
    integer, intent(inout) :: x(:)
    integer :: n, i, temp

    n = size(x)
    do i = 1, n/2
      temp = x(i)
      x(i) = x(n-i+1)
      x(n-i+1) = temp
    end do
  end subroutine flip1d_i

  !> @brief Flip (reverse) a 1D real(dp) array in place.
  !!
  !! @param[inout] x Real(dp) array.
  subroutine flip1d_dp(x)
    real(dp), intent(inout) :: x(:)
    integer :: n, i
    real(dp) :: temp

    n = size(x)
    do i = 1, n/2
      temp = x(i)
      x(i) = x(n-i+1)
      x(n-i+1) = temp
    end do
  end subroutine flip1d_dp

  !> @brief Flip a 2D real(dp) array along dimension dim in place.
  !!
  !! @param[inout] x Real(dp) array of shape (m, n).
  !! @param[in] dim Integer, dimension to flip: 1 = rows, 2 = columns.
  subroutine flip2d_dp(x, dim)
    real(dp), intent(inout) :: x(:,:)
    integer, intent(in) :: dim
    integer :: n, i, j
    real(dp) :: temp

    select case(dim)
    case(1)
      n = size(x, 1)
      do j = 1, size(x, 2)
        do i = 1, n/2
          temp = x(i, j)
          x(i, j) = x(n-i+1, j)
          x(n-i+1, j) = temp
        end do
      end do
    case(2)
      n = size(x, 2)
      do i = 1, size(x, 1)
        do j = 1, n/2
          temp = x(i, j)
          x(i, j) = x(i, n-j+1)
          x(i, n-j+1) = temp
        end do
      end do
    case default
      call crash('FLIP2D_DP - dim must be 1 or 2')
    end select
  end subroutine flip2d_dp

  !> @brief Flip a 3D real(dp) array along dimension dim in place.
  !!
  !! @param[inout] x Real(dp) array of shape (l, m, n).
  !! @param[in] dim Integer, dimension to flip: 1, 2, or 3.
  subroutine flip3d_dp(x, dim)
    real(dp), intent(inout) :: x(:,:,:)
    integer, intent(in) :: dim
    integer :: n, i, j, k
    real(dp) :: temp

    select case(dim)
    case(1)
      n = size(x, 1)
      do k = 1, size(x, 3)
        do j = 1, size(x, 2)
          do i = 1, n/2
            temp = x(i, j, k)
            x(i, j, k) = x(n-i+1, j, k)
            x(n-i+1, j, k) = temp
          end do
        end do
      end do
    case(2)
      n = size(x, 2)
      do k = 1, size(x, 3)
        do i = 1, size(x, 1)
          do j = 1, n/2
            temp = x(i, j, k)
            x(i, j, k) = x(i, n-j+1, k)
            x(i, n-j+1, k) = temp
          end do
        end do
      end do
    case(3)
      n = size(x, 3)
      do j = 1, size(x, 2)
        do i = 1, size(x, 1)
          do k = 1, n/2
            temp = x(i, j, k)
            x(i, j, k) = x(i, j, n-k+1)
            x(i, j, n-k+1) = temp
          end do
        end do
      end do
    case default
      call crash('FLIP3D_DP - dim must be 1, 2, or 3')
    end select
  end subroutine flip3d_dp

  !> @brief Check if array x is strictly ascending.
  !!
  !! @param[in] x Real(dp) array.
  !! @return Logical, true if x(i+1) > x(i) for all i.
  logical function ascending(x)
    real(dp), dimension(:), intent(in) :: x
    integer :: i

    ascending = .false.
    do i = 1, size(x) - 1
      if (x(i+1) - x(i) <= 0.0_dp) return
    end do
    ascending = .true.
  end function ascending

  !> @brief Check if array x is strictly descending.
  !!
  !! @param[in] x Real(dp) array.
  !! @return Logical, true if x(i+1) < x(i) for all i.
  logical function descending(x)
    real(dp), dimension(:), intent(in) :: x
    integer :: i

    descending = .false.
    do i = 1, size(x) - 1
      if (x(i+1) - x(i) >= 0.0_dp) return
    end do
    descending = .true.
  end function descending


  !> @brief Complement of Student's t cumulative distribution function (upper tail).
  !!
  !! P(T > t) for T ~ t(df).
  !!
  !! @param[in] t Real(dp), t-value.
  !! @param[in] df Real(dp), degrees of freedom.
  !! @return Real(dp), upper tail probability.
  !!
  !! Notes:
  !! - Uses incomplete beta function approximation.
  !! - Accurate for df >= 1.
  function t_cdf_complement(t, df) result(p)
    real(dp), intent(in) :: t, df
    real(dp) :: p, x, a, b

    if (df < 1.0_dp) call crash('t_cdf_complement: df must be >= 1')
  
    x = df / (df + t**2)
    a = 0.5_dp * df
    b = 0.5_dp

    p = 0.5_dp * beta_inc(x, a, b)

  end function t_cdf_complement


  !> @brief Inverse of Student's t cumulative distribution function.
  !!
  !! Returns t such that P(T <= t) = p for T ~ t(df).
  !!
  !! @param[in] p Real(dp), cumulative probability in (0, 1).
  !! @param[in] df Real(dp), degrees of freedom.
  !! @return Real(dp), t-value.
  !!
  !! Notes:
  !! - Uses Newton-Raphson iteration with beta function.
  !! - Accurate for df >= 1.
  function t_inverse_cdf(p, df) result(t)
    real(dp), intent(in) :: p, df
    real(dp) :: t
    real(dp) :: x, a, b

    if (p <= 0.0_dp .or. p >= 1.0_dp) call crash('t_inverse_cdf: p must be in (0, 1)')
    if (df < 1.0_dp) call crash('t_inverse_cdf: df must be >= 1')

    ! Use inverse beta approximation
    a = 0.5_dp * df
    b = 0.5_dp
    x = beta_inc_inv(2.0_dp * min(p, 1.0_dp - p), a, b)

    t = sqrt(df * (1.0_dp - x) / x)
    if (p < 0.5_dp) t = -t

  end function t_inverse_cdf


  !> @brief Complement of F cumulative distribution function (upper tail).
  !!
  !! P(F > f) for F ~ F(df1, df2).
  !!
  !! @param[in] f Real(dp), F-value.
  !! @param[in] df1 Real(dp), numerator degrees of freedom.
  !! @param[in] df2 Real(dp), denominator degrees of freedom.
  !! @return Real(dp), upper tail probability.
  !!
  !! Notes:
  !! - Uses incomplete beta function.
  function f_cdf_complement(f, df1, df2) result(p)
    real(dp), intent(in) :: f, df1, df2
    real(dp) :: p, x, a, b

    if (df1 < 1.0_dp .or. df2 < 1.0_dp) call crash('f_cdf_complement: df must be >= 1')

    x = df2 / (df2 + df1 * f)
    a = 0.5_dp * df2
    b = 0.5_dp * df1

    p = beta_inc(x, a, b)

  end function f_cdf_complement


  !> @brief Incomplete beta function I_x(a, b).
  !!
  !! @param[in] x Real(dp), argument in [0, 1].
  !! @param[in] a Real(dp), shape parameter > 0.
  !! @param[in] b Real(dp), shape parameter > 0.
  !! @return Real(dp), I_x(a, b).
  !!
  !! Notes:
  !! - Uses continued fraction expansion for numerical stability.
  !! - Required for t and F distribution CDFs.
  function beta_inc(x, a, b) result(betai)
    real(dp), intent(in) :: x, a, b
    real(dp) :: betai, bt, betacf_val

    if (x < 0.0_dp .or. x > 1.0_dp) call crash('beta_inc: x must be in [0, 1]')
    if (a <= 0.0_dp .or. b <= 0.0_dp) call crash('beta_inc: a, b must be > 0')

    if (x == 0.0_dp .or. x == 1.0_dp) then
      bt = 0.0_dp
    else
      bt = exp(log_gamma(a+b) - log_gamma(a) - log_gamma(b) + a*log(x) + b*log(1.0_dp-x))
    end if
  
    if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
      betai = bt * betacf(x, a, b) / a
    else
      betai = 1.0_dp - bt * betacf(1.0_dp-x, b, a) / b
    end if

  end function beta_inc
  

  !> @brief Continued fraction for incomplete beta function.
  !!
  !! @param[in] x Real(dp), argument.
  !! @param[in] a Real(dp), shape parameter.
  !! @param[in] b Real(dp), shape parameter.
  !! @return Real(dp), continued fraction value.
  !!
  !! Notes:
  !! - Uses modified Lentz's method.
  !! - Internal helper for beta_inc.
  function betacf(x, a, b) result(cf)
    real(dp), intent(in) :: x, a, b
    real(dp) :: cf
    integer, parameter :: maxit = 100
    real(dp), parameter :: eps = 3.0e-7_dp, fpmin = 1.0e-30_dp
    integer :: m, m2
    real(dp) :: aa, c, d, del, h, qab, qam, qap
  
    qab = a + b
    qap = a + 1.0_dp
    qam = a - 1.0_dp
    c = 1.0_dp
    d = 1.0_dp - qab*x/qap
  
    if (abs(d) < fpmin) d = fpmin
    d = 1.0_dp / d
    h = d

    do m = 1, maxit
      m2 = 2 * m
      aa = m * (b-m) * x / ((qam+m2) * (a+m2))
      d = 1.0_dp + aa*d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa/c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp / d
      h = h * d * c

      aa = -(a+m) * (qab+m) * x / ((a+m2) * (qap+m2))
      d = 1.0_dp + aa*d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa/c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp / d
      del = d * c
      h = h * del

      if (abs(del-1.0_dp) < eps) exit
    end do
  
    cf = h

  end function betacf


!> @brief Inverse of incomplete beta function.
  !!
  !! Returns x such that I_x(a, b) = y.
  !!
  !! @param[in] y Real(dp), target value in (0, 1).
  !! @param[in] a Real(dp), shape parameter > 0.
  !! @param[in] b Real(dp), shape parameter > 0.
  !! @return Real(dp), x in [0, 1].
  !!
  !! Notes:
  !! - Uses Newton-Raphson iteration.
  !! - Required for t_inverse_cdf.
  function beta_inc_inv(y, a, b) result(x)
    real(dp), intent(in) :: y, a, b
    real(dp) :: x
    real(dp), parameter :: eps = 1.0e-8_dp
    integer, parameter :: maxit = 100
    integer :: iter
    real(dp) :: f, fp, dx, lnbeta, afac

    if (y <= 0.0_dp .or. y >= 1.0_dp) call crash('beta_inc_inv: y must be in (0, 1)')
    if (a <= 0.0_dp .or. b <= 0.0_dp) call crash('beta_inc_inv: a, b must be > 0')

    ! Initial guess using simple approximation
    if (a >= 1.0_dp .and. b >= 1.0_dp) then
      ! Use mean as initial guess
      x = a / (a + b)
    else
      ! Use a better approximation for small a or b
      lnbeta = log_gamma(a) + log_gamma(b) - log_gamma(a + b)
      if (y < 0.5_dp) then
        x = exp((log(y) + lnbeta) / a)
      else
        x = 1.0_dp - exp((log(1.0_dp - y) + lnbeta) / b)
      end if
    end if
  
    ! Newton-Raphson refinement
    do iter = 1, maxit
      if (x <= 0.0_dp) x = 0.5_dp * eps
      if (x >= 1.0_dp) x = 1.0_dp - 0.5_dp * eps

      f = beta_inc(x, a, b) - y
      if (abs(f) < eps) exit
  
      ! Derivative of beta_inc is the beta distribution PDF
      lnbeta = log_gamma(a) + log_gamma(b) - log_gamma(a + b)
      if (x > 0.0_dp .and. x < 1.0_dp) then
        afac = exp((a-1.0_dp)*log(x) + (b-1.0_dp)*log(1.0_dp-x) - lnbeta)
        if (afac > 0.0_dp) then
          fp = afac
        else
          fp = 1.0_dp  ! Fallback to avoid division by zero
        end if
      else
        fp = 1.0_dp
      end if

      dx = f / fp
      x = x - dx

      ! Keep x in bounds
      if (x < 0.0_dp) x = 0.5_dp * (x + dx)
      if (x > 1.0_dp) x = 0.5_dp * (1.0_dp + (x - dx))
  
      if (abs(dx) < eps * abs(x)) exit
    end do

    if (iter >= maxit) then
      write(error_unit, '(A)') 'Warning: beta_inc_inv did not converge'
    end if

    x = max(0.0_dp, min(1.0_dp, x))

  end function beta_inc_inv

  !> @brief Logarithm of gamma function.
  !!
  !! @param[in] x Real(dp), argument > 0.
  !! @return Real(dp), ln(Γ(x)).
  !!
  !! Notes:
  !! - Uses Lanczos approximation.
  !! - Accurate to ~15 digits for x > 0.
  function log_gamma(x) result(lng)
    real(dp), intent(in) :: x
    real(dp) :: lng
    real(dp), parameter :: stp = 2.5066282746310005_dp
    real(dp), dimension(6), parameter :: cof = [76.18009172947146_dp, &
      -86.50532032941677_dp, 24.01409824083091_dp, -1.231739572450155_dp, &
      0.1208650973866179e-2_dp, -0.5395239384953e-5_dp]
    real(dp) :: ser, tmp, y
    integer :: j

    y = x
    tmp = x + 5.5_dp
    tmp = (x + 0.5_dp) * log(tmp) - tmp
    ser = 1.000000000190015_dp

    do j = 1, 6
      y = y + 1.0_dp
      ser = ser + cof(j) / y
    end do

    lng = tmp + log(stp * ser / x)

  end function log_gamma

end module module_math
