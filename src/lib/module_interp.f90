module module_interp

  use module_types, only : dp
  use module_constants, only: nan

  implicit none (type, external)
  private

  ! Kind selection: adjust if you already define rk elsewhere
  integer, parameter :: rk = kind(1.0d0)

  public :: interp1d, interp1d_vec, interp2d
  public :: spline_coeffs, spline_eval

contains

  !========================
  ! Binary search helper
  ! Returns i such that x(i) <= xv < x(i+1), with 1 <= i <= n-1
  ! If xv < x(1) returns 0; if xv >= x(n) returns n
  !========================
  pure function bracket_index(x, xv) result(i)
    real(rk), intent(in) :: x(:), xv
    integer :: i
    integer :: lo, hi, mid, n
    n = size(x)
    if (n < 2) then
      i = -1
      return
    end if
    if (xv < x(1)) then
      i = 0
      return
    else if (xv >= x(n)) then
      i = n
      return
    end if
    lo = 1
    hi = n
    do
      if (hi - lo <= 1) exit
      mid = (lo + hi)/2
      if (x(mid) <= xv) then
        lo = mid
      else
        hi = mid
      end if
    end do
    i = lo
  end function bracket_index

  !==========================================================
  ! Natural cubic spline coefficients for strictly increasing x
  !
  ! Given x(1:n), y(1:n), returns second derivatives y2(1:n)
  ! such that a cubic spline with natural boundary conditions
  ! is defined. Algorithm based on Numerical Recipes.
  !==========================================================
  pure subroutine spline_coeffs(x, y, y2)
    real(rk), intent(in)  :: x(:), y(:)
    real(rk), intent(out) :: y2(:)
    integer :: n, i, k
    real(rk), allocatable :: u(:)
    real(rk) :: sig, p

    n = size(x)
    if (size(y) /= n .or. size(y2) /= n) then
      y2 = nan()
      return
    end if
    if (n < 2) then
      y2 = 0.0_rk
      return
    else if (n == 2) then
      y2 = 0.0_rk
      return
    end if

    ! Validate strictly increasing x
    do i = 2, n
      if (x(i) <= x(i-1)) then
        y2 = nan()
        return
      end if
    end do

    allocate(u(n-1))
    y2(1) = 0.0_rk
    u(1)  = 0.0_rk

    do i = 2, n-1
      sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
      p   = sig*y2(i-1) + 2.0_rk
      if (p == 0.0_rk) then
        y2 = nan()
        deallocate(u)
        return
      end if
      y2(i) = (sig - 1.0_rk)/p
      u(i)  = (6.0_rk*((y(i+1)-y(i))/(x(i+1)-x(i)) - (y(i)-y(i-1))/(x(i)-x(i-1))) / (x(i+1)-x(i-1)) - sig*u(i-1)) / p
    end do

    y2(n) = 0.0_rk
    do k = n-1, 1, -1
      y2(k) = y2(k)*y2(k+1) + u(k)
    end do

    deallocate(u)
  end subroutine spline_coeffs

  !==========================================================
  ! Evaluate natural cubic spline at xv.
  ! Inputs: x(1:n), y(1:n), y2(1:n) from spline_coeffs
  !==========================================================
  pure function spline_eval(x, y, y2, xv) result(fv)
    real(rk), intent(in) :: x(:), y(:), y2(:), xv
    real(rk) :: fv
    integer :: n, klo, khi
    real(rk) :: h, a, b

    n = size(x)
    if (size(y) /= n .or. size(y2) /= n) then
      fv = nan()
      return
    end if
    if (n < 2) then
      fv = nan()
      return
    end if

    ! Handle outside range by clamping to end intervals
    if (xv <= x(1)) then
      klo = 1;  khi = 2
    else if (xv >= x(n)) then
      klo = n-1; khi = n
    else
      klo = bracket_index(x, xv)
      if (klo < 1) klo = 1
      khi = klo + 1
    end if

    h = x(khi) - x(klo)
    if (h <= 0.0_rk) then
      fv = nan()
      return
    end if

    a = (x(khi) - xv)/h
    b = (xv - x(klo))/h
    fv = a*y(klo) + b*y(khi) + ((a**3 - a)*y2(klo) + (b**3 - b)*y2(khi))*(h*h)/6.0_rk
  end function spline_eval

  !==========================================================
  ! 1D interpolation for a single point
  !
  ! method options: 'nearest' (default), 'linear', 'spline'
  ! fill_value: used when xv is out-of-bounds (default NaN)
  !
  ! Requirements:
  ! - x strictly increasing; y same length as x
  !==========================================================
  pure function interp1d(x, y, xv, method, fill_value) result(fv)
    real(rk), intent(in)           :: x(:), y(:), xv
    character(len=*), intent(in), optional :: method
    real(rk), intent(in), optional :: fill_value
    real(rk) :: fv
    character(len=:), allocatable :: m
    integer :: n, i
    real(rk) :: t, yv
    real(rk), allocatable :: y2(:)

    n = size(x)
    if (size(y) /= n .or. n < 1) then
      fv = nan(); return
    end if

    ! Default method
    if (present(method)) then
      m = to_lower(trim(method))
    else
      m = 'nearest'
    end if

    ! Default fill
    if (present(fill_value)) then
      yv = fill_value
    else
      yv = nan()
    end if

    ! Validate increasing x
    if (n >= 2) then
      do i = 2, n
        if (x(i) <= x(i-1)) then
          fv = nan(); return
        end if
      end do
    end if

    ! Trivial n=1
    if (n == 1) then
      fv = y(1)
      return
    end if

    ! Out-of-bounds
    if (xv < x(1) .or. xv > x(n)) then
      fv = yv
      return
    end if

    select case (m)
    case ('nearest')
      i = bracket_index(x, xv)
      if (i == 0) then
        fv = y(1)
      else if (i >= n) then
        fv = y(n)
      else
        ! Choose nearest between i and i+1
        if (abs(xv - x(i)) <= abs(x(i+1) - xv)) then
          fv = y(i)
        else
          fv = y(i+1)
        end if
      end if

    case ('linear')
      i = bracket_index(x, xv)
      if (i <= 0) then
        fv = y(1); return
      else if (i >= n) then
        fv = y(n); return
      end if
      t = (xv - x(i)) / (x(i+1) - x(i))
      fv = (1.0_rk - t)*y(i) + t*y(i+1)

    case ('spline')
      allocate(y2(n))
      call spline_coeffs(x, y, y2)
      fv = spline_eval(x, y, y2, xv)
      deallocate(y2)

    case default
      ! Fallback to linear for unknown method
      i = bracket_index(x, xv)
      if (i <= 0) then
        fv = y(1); return
      else if (i >= n) then
        fv = y(n); return
      end if
      t = (xv - x(i)) / (x(i+1) - x(i))
      fv = (1.0_rk - t)*y(i) + t*y(i+1)
    end select
  end function interp1d

  !==========================================================
  ! Vectorized 1D interpolation for many query points
  ! method: 'nearest' (default), 'linear', 'spline'
  ! fill_value applied for out-of-bounds
  !==========================================================
  pure subroutine interp1d_vec(x, y, xq, yq, method, fill_value)
    real(rk), intent(in)           :: x(:), y(:)
    real(rk), intent(in)           :: xq(:)
    real(rk), intent(out)          :: yq(:)
    character(len=*), intent(in), optional :: method
    real(rk), intent(in), optional :: fill_value

    character(len=:), allocatable :: m
    integer :: n, nq, i, j
    real(rk) :: t, fv
    real(rk), allocatable :: y2(:)

    n  = size(x)
    nq = size(xq)
    if (size(y) /= n .or. size(yq) /= nq) then
      if (nq > 0) yq = nan()
      return
    end if

    ! Default method and fill
    if (present(method)) then
      m = to_lower(trim(method))
    else
      m = 'nearest'
    end if
    if (present(fill_value)) then
      fv = fill_value
    else
      fv = nan()
    end if

    ! Validate x strictly increasing
    if (n >= 2) then
      do i = 2, n
        if (x(i) <= x(i-1)) then
          yq = nan()
          return
        end if
      end do
    end if

    ! Trivial n cases
    if (n == 0) then
      if (nq > 0) yq = nan()
      return
    else if (n == 1) then
      yq = y(1)
      return
    end if

    select case (m)
    case ('nearest')
      do j = 1, nq
        i = bracket_index(x, xq(j))
        if (i == 0) then
          yq(j) = fv
        else if (i >= n) then
          yq(j) = fv
        else
          if (xq(j) < x(1) .or. xq(j) > x(n)) then
            yq(j) = fv
          else if (abs(xq(j) - x(i)) <= abs(x(i+1) - xq(j))) then
            yq(j) = y(i)
          else
            yq(j) = y(i+1)
          end if
        end if
      end do

    case ('linear')
      do j = 1, nq
        if (xq(j) < x(1) .or. xq(j) > x(n)) then
          yq(j) = fv
        else
          i = bracket_index(x, xq(j))
          if (i <= 0) then
            yq(j) = y(1)
          else if (i >= n) then
            yq(j) = y(n)
          else
            t = (xq(j) - x(i)) / (x(i+1) - x(i))
            yq(j) = (1.0_rk - t)*y(i) + t*y(i+1)
          end if
        end if
      end do

    case ('spline')
      allocate(y2(n))
      call spline_coeffs(x, y, y2)
      do j = 1, nq
        if (xq(j) < x(1) .or. xq(j) > x(n)) then
          yq(j) = fv
        else
          yq(j) = spline_eval(x, y, y2, xq(j))
        end if
      end do
      deallocate(y2)

    case default
      ! Fallback to linear
      do j = 1, nq
        if (xq(j) < x(1) .or. xq(j) > x(n)) then
          yq(j) = fv
        else
          i = bracket_index(x, xq(j))
          if (i <= 0) then
            yq(j) = y(1)
          else if (i >= n) then
            yq(j) = y(n)
          else
            t = (xq(j) - x(i)) / (x(i+1) - x(i))
            yq(j) = (1.0_rk - t)*y(i) + t*y(i+1)
          end if
        end if
      end do
    end select
  end subroutine interp1d_vec

  !==========================================================
  ! 2D interpolation on a rectilinear grid
  !
  ! Inputs:
  ! - x(1:nx), y(1:ny): strictly increasing axes
  ! - z(nx, ny): function values at grid points (Fortran order)
  ! - xq, yq: query points (same length)
  !
  ! method: 'bilinear' (default), 'nearest'
  ! fill_value for out-of-bounds (default NaN)
  !==========================================================
  pure subroutine interp2d(x, y, z, xq, yq, zq, method, fill_value)
    real(rk), intent(in)           :: x(:), y(:)
    real(rk), intent(in)           :: z(:,:)
    real(rk), intent(in)           :: xq(:), yq(:)
    real(rk), intent(out)          :: zq(:)
    character(len=*), intent(in), optional :: method
    real(rk), intent(in), optional :: fill_value

    integer :: nx, ny, nq
    integer :: i, j, ix, iy
    real(rk) :: fv, tx, ty
    character(len=:), allocatable :: m

    nx = size(x); ny = size(y)
    if (size(z,1) /= nx .or. size(z,2) /= ny) then
      if (size(zq) > 0) zq = nan()
      return
    end if
    nq = size(xq)
    if (size(yq) /= nq .or. size(zq) /= nq) then
      if (nq > 0) zq = nan()
      return
    end if

    ! Defaults
    if (present(method)) then
      m = to_lower(trim(method))
    else
      m = 'bilinear'
    end if
    if (present(fill_value)) then
      fv = fill_value
    else
      fv = nan()
    end if

    ! Validate monotonicity
    if (nx < 2 .or. ny < 2) then
      if (nq > 0) zq = nan()
      return
    end if
    do i = 2, nx
      if (x(i) <= x(i-1)) then
        zq = nan(); return
      end if
    end do
    do j = 2, ny
      if (y(j) <= y(j-1)) then
        zq = nan(); return
      end if
    end do

    select case (m)
    case ('nearest')
      do i = 1, nq
        ix = bracket_index(x, xq(i))
        iy = bracket_index(y, yq(i))
        if (xq(i) < x(1) .or. xq(i) > x(nx) .or. yq(i) < y(1) .or. yq(i) > y(ny)) then
          zq(i) = fv
        else
          ! resolve nearest in x
          if (ix <= 0) then
            ix = 1
          else if (ix >= nx) then
            ix = nx
          else
            if (abs(xq(i)-x(ix)) > abs(x(ix+1)-xq(i))) ix = ix+1
          end if
          ! resolve nearest in y
          if (iy <= 0) then
            iy = 1
          else if (iy >= ny) then
            iy = ny
          else
            if (abs(yq(i)-y(iy)) > abs(y(iy+1)-yq(i))) iy = iy+1
          end if
          zq(i) = z(ix, iy)
        end if
      end do

    case default   ! 'bilinear'
      do i = 1, nq
        if (xq(i) < x(1) .or. xq(i) > x(nx) .or. yq(i) < y(1) .or. yq(i) > y(ny)) then
          zq(i) = fv
        else
          ix = bracket_index(x, xq(i))
          iy = bracket_index(y, yq(i))
          if (ix <= 0) then
            ix = 1
          else if (ix >= nx) then
            ix = nx-1
          end if
          if (iy <= 0) then
            iy = 1
          else if (iy >= ny) then
            iy = ny-1
          end if
          tx = (xq(i) - x(ix)) / (x(ix+1) - x(ix))
          ty = (yq(i) - y(iy)) / (y(iy+1) - y(iy))
          zq(i) = (1.0_rk-tx)*(1.0_rk-ty)*z(ix  ,iy  ) + &
                   tx      *(1.0_rk-ty)*z(ix+1,iy  ) + &
                  (1.0_rk-tx)*ty      *z(ix  ,iy+1) + &
                   tx      *ty      *z(ix+1,iy+1)
        end if
      end do
    end select
  end subroutine interp2d

  !========================
  ! to_lower helper
  !========================
  pure function to_lower(s) result(t)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: t
    integer :: i, ia, iz
    ia = ichar('A'); iz = ichar('Z')
    do i = 1, len(s)
      if (ichar(s(i:i)) >= ia .and. ichar(s(i:i)) <= iz) then
        t(i:i) = achar(ichar(s(i:i)) + 32)
      else
        t(i:i) = s(i:i)
      end if
    end do
  end function to_lower

end module module_interp
