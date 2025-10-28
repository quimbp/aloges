module module_interp

use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
use module_types, only : dp
use module_constants, only: nan
use module_tools, only: lowercase, ltrim, crash

implicit none (type, external)
private

  public :: interp1d, interp2d

!> @brief Generic interface for creating arithmetic 
!> sequences (integer or real).
interface interp1d
  module procedure interp1d_r, interp1d_vec
end interface interp1d

contains

  !========================
  ! Binary search helper
  ! Returns i such that x(i) <= xv < x(i+1), with 1 <= i <= n-1
  ! If xv < x(1) returns 0; if xv >= x(n) returns n
  !========================
  pure function bracket_index(x, xv) result(i)
    real(dp), intent(in) :: x(:), xv
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
  ! Natural cubic splines
  !==========================================================

  !==========================================================
  ! Natural cubic spline coefficients for strictly increasing x
  !
  ! Given x(1:n), y(1:n), returns second derivatives y2(1:n)
  ! such that a cubic spline with natural boundary conditions
  ! is defined. Algorithm based on Numerical Recipes.
  !==========================================================
  pure subroutine spline_coeffs(x, y, y2)
    real(dp), intent(in)  :: x(:), y(:)
    real(dp), intent(out) :: y2(:)
    integer :: n, i, k
    real(dp), allocatable :: u(:)
    real(dp) :: sig, p

    n = size(x)
    if (size(y) /= n .or. size(y2) /= n) then
      y2 = nan()
      return
    end if
    if (n < 2) then
      y2 = 0.0_dp
      return
    else if (n == 2) then
      y2 = 0.0_dp
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
    y2(1) = 0.0_dp
    u(1)  = 0.0_dp

    do i = 2, n-1
      sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
      p   = sig*y2(i-1) + 2.0_dp
      if (p == 0.0_dp) then
        y2 = nan()
        deallocate(u)
        return
      end if
      y2(i) = (sig - 1.0_dp)/p
      u(i)  = (6.0_dp*((y(i+1)-y(i))/(x(i+1)-x(i)) - (y(i)-y(i-1))/(x(i)-x(i-1))) / (x(i+1)-x(i-1)) - sig*u(i-1)) / p
    end do

    y2(n) = 0.0_dp
    do k = n-1, 1, -1
      y2(k) = y2(k)*y2(k+1) + u(k)
    end do

    deallocate(u)
  end subroutine spline_coeffs


  !====
  ! General cubic spline coefficients (second derivatives) for strictly increasing x
  ! Supports boundary conditions:
  !   - 'natural' : y''(x1)=0, y''(xn)=0
  !   - 'clamped' : y'(x1)=yp_left, y'(xn)=yp_right
  !   - 'second'  : y''(x1)=yp_left, y''(xn)=yp_right
  !
  ! Inputs:
  !   x(n), y(n)      : data points, x strictly increasing
  !   bc_left         : 'natural' | 'clamped' | 'second'
  !   bc_right        : 'natural' | 'clamped' | 'second'
  !   yp_left         : (optional) first/second derivative at x(1) if bc_left='clamped'
  !   yp_right        : (optional) first/second derivative at x(n) if bc_right='clamped'
  !
  ! Output:
  !   y2(n) : second derivatives at the data points
  !
  ! Notes:
  !   - Algorithm based on Numerical Recipes (tridiagonal system).
  !   - On any input error, returns y2=NaN.
  !====
  pure subroutine spline_coeffs_bc(x, y, y2, bc_left, bc_right, &
                                             yp_left, yp_right)
    real(dp), intent(in)           :: x(:), y(:)
    real(dp), intent(out)          :: y2(:)
    character(len=*), intent(in)   :: bc_left, bc_right
    real(dp), intent(in), optional :: yp_left, yp_right

    integer :: n, i, k
    real(dp), allocatable :: u(:)
    real(dp) :: sig, p
    real(dp) :: dx1, dxn
    logical :: left_natural, right_natural
    logical :: left_clamped, right_clamped
    logical :: left_second,  right_second

    n = size(x)
    if (size(y) /= n .or. size(y2) /= n) then
      y2 = nan()
      return
    end if
    if (n < 2) then
      y2 = 0.0_dp
      return
    else if (n == 2) then
      ! With two points, the spline is linear; set second derivatives to zero
      y2 = 0.0_dp
      return
    end if

    ! Validate strictly increasing x
    do i = 2, n
      if (x(i) <= x(i-1)) then
        y2 = nan()
        return
      end if
    end do

    left_natural  = (bc_left == 'natural')
    right_natural = (bc_right == 'natural')
    left_clamped  = (bc_left == 'clamped' .or.  bc_left == 'first')
    right_clamped = (bc_right == 'clamped' .or. bc_right == 'first')
    left_second   = (bc_left == 'second')
    right_second  = (bc_right == 'second')

    ! Validate needed parameters for boundary conditions
    if (left_clamped .and. .not. present(yp_left)) then
      y2 = nan(); return
    end if
    if (right_clamped .and. .not. present(yp_right)) then
      y2 = nan(); return
    end if
    if (left_second .and. .not. present(yp_left)) then
      y2 = nan(); return
    end if
    if (right_second .and. .not. present(yp_right)) then
      y2 = nan(); return
    end if

    allocate(u(n-1))

    ! Left boundary
    if (left_natural) then
      y2(1) = 0.0_dp
      u(1)  = 0.0_dp
    else if (left_clamped) then
      dx1 = x(2) - x(1)
      if (dx1 <= 0.0_dp) then
        y2 = nan(); deallocate(u); return
      end if
      y2(1) = -0.5_dp
      u(1)  = (3.0_dp / dx1) * ((y(2) - y(1))/dx1 - yp_left)
    else if (left_second) then
      y2(1) = 0.0_dp
      u(1)  = 0.5_dp * yp_left
    else
      ! Unknown bc
      y2 = nan(); deallocate(u); return
    end if

    ! Decomposition of the tridiagonal matrix
    do i = 2, n-1
      sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
      p   = sig*y2(i-1) + 2.0_dp
      if (p == 0.0_dp) then
        y2 = nan(); deallocate(u); return
      end if
      y2(i) = (sig - 1.0_dp) / p
      u(i)  = ( 6.0_dp * ( ( (y(i+1)-y(i))/(x(i+1)-x(i)) ) - ( (y(i)-y(i-1))/(x(i)-x(i-1)) ) ) &
               / (x(i+1)-x(i-1)) - sig*u(i-1) ) / p
    end do

    ! Right boundary
    if (right_natural) then
      y2(n) = 0.0_dp
    else if (right_clamped) then
      dxn = x(n) - x(n-1)
      if (dxn <= 0.0_dp) then
        y2 = nan(); deallocate(u); return
      end if
      y2(n) = 0.5_dp
      u(n-1) = ( 3.0_dp / dxn ) * ( yp_right - (y(n) - y(n-1))/dxn )
    else if (right_second) then
      y2(n) = 0.0_dp
      ! fold right boundary into the back-substitution via u(n-1):
      u(n-1) = u(n-1) - 0.5_dp * yp_right
    else
      y2 = nan(); deallocate(u); return
    end if

    ! Back substitution
    do k = n-1, 1, -1
      y2(k) = y2(k)*y2(k+1) + u(k)
    end do

    deallocate(u)
  end subroutine spline_coeffs_bc


  !==========================================================
  ! Evaluate natural cubic spline at xv.
  ! Inputs: x(1:n), y(1:n), y2(1:n) from spline_coeffs
  !==========================================================
  pure function spline_eval(x, y, y2, xv) result(fv)
    real(dp), intent(in) :: x(:), y(:), y2(:), xv
    real(dp) :: fv
    integer :: n, klo, khi
    real(dp) :: h, a, b

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
    if (h <= 0.0_dp) then
      fv = nan()
      return
    end if

    a = (x(khi) - xv)/h
    b = (xv - x(klo))/h
    fv = a*y(klo) + b*y(khi) + ((a**3 - a)*y2(klo) + (b**3 - b)*y2(khi))*(h*h)/6.0_dp
  end function spline_eval

  !==========================================================
  ! Akima splines
  !==========================================================

  !====  
  ! Akima spline derivative calculation  
  !   
  ! Computes first derivatives at each data point using Akima's algorithm  
  ! (Lancaster and Salkauskas, Curve and Surface Fitting, 1986, Academic Press, p.82)  
  !  
  ! Inputs:  
  !   x(n) : strictly increasing array of sample points  
  !   y(n) : function values at sample points  
  !  
  ! Output:  
  !   dy(n) : first derivatives at each sample point  
  !  
  ! Notes:  
  !   - Akima's method produces a smooth interpolant with good local behavior  
  !   - Less prone to oscillations than cubic splines  
  !   - Returns dy=NaN on error (non-monotonic x, n<2, etc.)  
  !====  

  pure subroutine akima_derivs (x, y, dy) 

    real(dp), intent(in)         :: x(:), y(:)
    real(dp), intent(out)        :: dy(:)

    integer  :: n,i
    real(dp) :: eps, w1, w2, denom
    real(dp), allocatable :: slopes(:)

    n = size(x)

    ! Validate inputs
    if (size(y).ne.n .or. size(dy).ne.n) then
      dy = nan()
      return
    endif
    if (n.lt.2) then
      dy = nan()
      return
    endif

    ! Validate increaing ordinates
    do i=2,n
      if (x(i).le.x(i-1)) then
        dy = nan()
        return
      endif
    enddo

    ! Special case: two points (linear)
    if (n.eq.2) then
      dy(1) = (y(2)-y(1)) / (x(2)-x(1))
      dy(2) = dy(1)
      return
    endif

    ! Compute slopes between consecutive points  
    ! We need slopes from i=-1 to i=n+1 (conceptually)  
    ! Store as slopes(1:n+3) where slopes(i+2) corresponds to slope at segment i  
    allocate(slopes(n+3))  
    
    ! Interior slopes: segments 1 to n-1  
    do i = 1, n-1  
      slopes(i+2) = (y(i+1) - y(i)) / (x(i+1) - x(i))  
    end do  

    ! Extrapolate slopes at boundaries  
    ! Repetition
    !slopes(2) = slopes(3)  
    !slopes(1) = slopes(3)  
    !slopes(n+2) = slopes(n+1)  
    !slopes(n+3) = slopes(n+1)  
    ! Reflection for better boundary behavior
    slopes(2) = 2.0_dp*slopes(3) - slopes(4)  
    slopes(1) = 2.0_dp*slopes(2) - slopes(3)  
    slopes(n+2) = 2.0_dp*slopes(n+1) - slopes(n)  
    slopes(n+3) = 2.0_dp*slopes(n+2) - slopes(n+1)  
    

    ! Small epsilon to avoid division by zero  
    eps = 1.0e-10_dp * abs(x(n) - x(1))  
   
    ! Compute weighted average of slopes for each point  
    ! At point i, we use slopes m_{i-2}, m_{i-1}, m_i, m_{i+1}  
    do i = 1, n  
      w1 = abs(slopes(i+3) - slopes(i+2))  ! |m_{i+1} - m_i|  
      w2 = abs(slopes(i+1) - slopes(i))    ! |m_{i-1} - m_{i-2}|  
      denom = w1 + w2 + eps  
      dy(i) = (w1 * slopes(i+1) + w2 * slopes(i+2)) / denom  
    end do   
    
!    ! Special treatment for endpoints if n >= 3  
!    ! Use parabolic fit for better accuracy  
!    if (n >= 3) then  
!      dy(1) = parabolic_deriv(y(2)-y(1), x(2)-x(1), y(3)-y(1), x(3)-x(1))  
!      dy(n) = parabolic_deriv(y(n-1)-y(n), x(n-1)-x(n), y(n-2)-y(n), x(n-2)-x(n))  
!    end if  
    
    deallocate(slopes)  
    
!  contains  
!    
!    !====  
!    ! Compute derivative from parabolic fit through origin  
!    ! Given two points (x1, u1) and (x2, u2), fits parabola u = a*x + b*x^2  
!    !====  
!    pure function parabolic_deriv(u1, x1, u2, x2) result(deriv)  
!      real(dp), intent(in) :: u1, x1, u2, x2  
!      real(dp) :: deriv  
!      real(dp) :: denom  
!      
!      denom = 1.0_dp/x1 - 1.0_dp/x2  
!      if (abs(denom) < 1.0e-14_dp) then  
!        deriv = u1/x1  ! Fallback to simple slope  
!      else  
!        deriv = (u1/x1**2 - u2/x2**2) / denom  
!      end if  
!    end function parabolic_deriv  
!    
  end subroutine akima_derivs  
  
  
  !====  
  ! Evaluate Akima spline at a single point  
  !  
  ! Inputs:  
  !   x(n)  : strictly increasing array of sample points  
  !   y(n)  : function values at sample points  
  !   dy(n) : first derivatives from akima_derivs  
  !   xv    : evaluation point  
  !  
  ! Output:  
  !   yv : interpolated value at xv  
  !  
  ! Notes:  
  !   - Uses Hermite cubic interpolation  
  !   - Extrapolates using nearest interval if xv is outside [x(1), x(n)]  
  !====  
  pure function akima_eval(x, y, dy, xv) result(yv)  
    real(dp), intent(in) :: x(:), y(:), dy(:), xv  
    real(dp) :: yv  
    integer :: n, i  
    real(dp) :: dx, t, s, h00, h10, h01, h11  
    
    n = size(x)  
    
    ! Validate inputs  
    if (size(y) /= n .or. size(dy) /= n .or. n < 2) then  
      yv = nan()  
      return  
    end if  
    
    ! Find interval containing xv  
    i = bracket_index(x, xv)  
    
    ! Handle boundary cases  
    if (i <= 0) then  
      i = 1  
    else if (i >= n) then  
      i = n - 1  
    end if  
    
    ! Compute normalized position in interval  
    dx = x(i+1) - x(i)  
    t  = (xv - x(i)) / dx  
    s  = 1.0_dp - t  
    
    ! Hermite basis functions  
    h00 = s * s * (1.0_dp + 2.0_dp * t)  
    h10 = t * s * s  
    h01 = t * t * (3.0_dp - 2.0_dp * t)  
    h11 = -t * t * s  
    
    ! Evaluate cubic Hermite polynomial  
    yv = h00 * y(i) + h10 * dx * dy(i) + h01 * y(i+1) + h11 * dx * dy(i+1)  
    
  end function akima_eval  


!====
  ! PCHIP (Piecewise Cubic Hermite Interpolating Polynomial)
  !====

  !====
  ! PCHIP derivative calculation
  !
  ! Computes first derivatives at each data point using PCHIP algorithm
  ! (Fritsch & Carlson, 1980, SIAM J. Numer. Anal.)
  !
  ! Inputs:
  !   x(n) : strictly increasing array of sample points
  !   y(n) : function values at sample points
  !
  ! Output:
  !   dy(n) : first derivatives at each sample point
  !
  ! Notes:
  !   - Monotone-preserving: no overshoots/undershoots
  !   - Shape-preserving: respects local trends in data
  !   - C1 continuous
  !   - Ideal for physical/scientific data
  !====
  pure subroutine pchip_derivs(x, y, dy)
    real(dp), intent(in)  :: x(:), y(:)
    real(dp), intent(out) :: dy(:)
    
    integer :: n, i
    real(dp), allocatable :: h(:), delta(:)
    real(dp) :: w1, w2
    
    n = size(x)
    
    ! Validate inputs
    if (size(y) /= n .or. size(dy) /= n) then
      dy = nan()
      return
    end if
    if (n < 2) then
      dy = nan()
      return
    endif
    
    do i = 2, n
      if (x(i) <= x(i-1)) then
        dy = nan()
        return
      endif
    end do
    
    ! Special case: two points (linear)
    if (n == 2) then
      dy(1) = (y(2) - y(1)) / (x(2) - x(1))
      dy(2) = dy(1)
      return
    end if
    
    ! Compute intervals and slopes
    allocate(h(n-1), delta(n-1))
    
    do i = 1, n-1
      h(i) = x(i+1) - x(i)
      delta(i) = (y(i+1) - y(i)) / h(i)
    end do
    
    ! Compute derivatives at interior points
    do i = 2, n-1
      ! Check for monotonicity
      if (delta(i-1) * delta(i) <= 0.0_dp) then
        ! Local extremum - set derivative to zero
        dy(i) = 0.0_dp
      else
        ! Weighted harmonic mean (preserves monotonicity)
        w1 = 2.0_dp * h(i) + h(i-1)
        w2 = h(i) + 2.0_dp * h(i-1)
        dy(i) = (w1 + w2) / (w1/delta(i-1) + w2/delta(i))
      end if
    end do
    
    ! Endpoint derivatives using non-centered scheme
    ! Left endpoint
    dy(1) = pchip_end_deriv(h(1), h(2), delta(1), delta(2))
    if (dy(1) * delta(1) < 0.0_dp) then
      dy(1) = 0.0_dp
    else if (delta(1) * delta(2) < 0.0_dp .and. abs(dy(1)) > abs(3.0_dp*delta(1))) then
      dy(1) = 3.0_dp * delta(1)
    end if
    
    ! Right endpoint
    dy(n) = pchip_end_deriv(h(n-1), h(n-2), delta(n-1), delta(n-2))
    if (dy(n) * delta(n-1) < 0.0_dp) then
      dy(n) = 0.0_dp
    else if (delta(n-1) * delta(n-2) < 0.0_dp .and. abs(dy(n)) > abs(3.0_dp*delta(n-1))) then
      dy(n) = 3.0_dp * delta(n-1)
    end if
    
    deallocate(h, delta)
    
  contains
    
    !====
    ! Compute endpoint derivative using one-sided formula
    !====
    pure function pchip_end_deriv(h1, h2, d1, d2) result(deriv)
      real(dp), intent(in) :: h1, h2, d1, d2
      real(dp) :: deriv
      
      ! One-sided three-point formula
      deriv = ((2.0_dp*h1 + h2)*d1 - h1*d2) / (h1 + h2)
      
    end function pchip_end_deriv
    
  end subroutine pchip_derivs
  
  !====
  ! Evaluate PCHIP interpolation at a single point
  !
  ! Inputs:
  !   x(n)  : strictly increasing array of sample points
  !   y(n)  : function values at sample points
  !   dy(n) : first derivatives from pchip_derivs
  !   xv    : evaluation point
  !
  ! Output:
  !   yv : interpolated value at xv
  !
  ! Notes:
  !   - Uses Hermite cubic interpolation (same as Akima evaluation)
  !   - Extrapolates using nearest interval if xv is outside [x(1), x(n)]
  !====
  function pchip_eval(x, y, dy, xv) result(yv)
    real(dp), intent(in) :: x(:), y(:), dy(:), xv
    real(dp) :: yv
    integer :: n, i
    real(dp) :: dx, t, s, h00, h10, h01, h11
    
    n = size(x)
    
    ! Validate inputs
    if (size(y) /= n .or. size(dy) /= n .or. n < 2) then
      yv = nan()
      return
    end if
    
    ! Find interval containing xv
    i = bracket_index(x, xv)
    
    ! Handle boundary cases
    if (i <= 0) then
      i = 1
    else if (i >= n) then
      i = n - 1
    end if
    
    ! Compute normalized position in interval
    dx = x(i+1) - x(i)
    t  = (xv - x(i)) / dx
    s  = 1.0_dp - t
    
    ! Hermite basis functions
    h00 = s * s * (1.0_dp + 2.0_dp * t)
    h10 = t * s * s
    h01 = t * t * (3.0_dp - 2.0_dp * t)
    h11 = -t * t * s
    
    ! Evaluate cubic Hermite polynomial
    yv = h00 * y(i) + h10 * dx * dy(i) + h01 * y(i+1) + h11 * dx * dy(i+1)
    
  end function pchip_eval



!====
  ! Compute cubic interpolation coefficients
  !
  ! Uses natural cubic splines in the interior, switches to linear at edges.
  ! Pre-computes all coefficients for efficient vectorized evaluation.
  !
  ! Inputs:
  !   x(n) : strictly increasing array of sample points
  !   y(n) : function values at sample points
  !
  ! Output:
  !   coeffs(4, n-1) : cubic coefficients for each interval
  !    coeffs(:,i) = [a, b, c, d] for interval [x(i), x(i+1)]
  !    where f(t) = a + b*t + c*t^2 + d*t^3
  !    with t = (xv - x(i)) / (x(i+1) - x(i))
  !
  ! Notes:
  !   - First and last intervals use linear interpolation
  !   - Interior intervals use natural cubic splines
  !====
  subroutine cubic_coeffs(x, y, coeffs)
    real(dp), intent(in)  :: x(:), y(:)
    real(dp), intent(out) :: coeffs(:,:)
    
    integer :: n, i
    real(dp), allocatable :: y2(:)
    real(dp) :: h
    
    n = size(x)
    
    ! Validate inputs
    if (size(y) /= n .or. size(coeffs, 2) /= n-1 .or. size(coeffs, 1) /= 4) then
    coeffs = nan()
    return
    end if
    
    if (n < 2) then
    coeffs = nan()
    return
    end if
    
    ! Validate strictly increasing x
    do i = 2, n
    if (x(i) <= x(i-1)) then
    coeffs = nan()
    return
    end if
    end do
    
    ! Special case: only 2 points (linear)
    if (n == 2) then
    coeffs(1, 1) = y(1)
    coeffs(2, 1) = y(2) - y(1)
    coeffs(3, 1) = 0.0_dp
    coeffs(4, 1) = 0.0_dp
    return
    end if
    
    ! Compute natural cubic spline second derivatives
    allocate(y2(n))
    call spline_coeffs(x, y, y2)
    
    ! Check for errors in spline computation
    if (any(y2 /= y2)) then  ! NaN check
      coeffs = nan()
      deallocate(y2)
      return
    end if
    
    ! Compute coefficients for each interval
    do i = 1, n-1

      h = x(i+1) - x(i)
    
      ! Use linear interpolation at edges
      if (i == 1 .or. i == n-1) then
        coeffs(1, i) = y(i)
        coeffs(2, i) = y(i+1) - y(i)
        coeffs(3, i) = 0.0_dp
        coeffs(4, i) = 0.0_dp
      else
    ! Cubic spline in interior
    ! Standard spline formula: f(x) = a*(x-xi) + b*(x-xi)^2 + c*(x-xi)^3 
    ! With normalized t = (x-xi)/h, we have: f(t) = y(i) + b*t + c*t^2 + d*t^3  
    ! From spline formula: f(t) = (1-t)*y(i) + t*y(i+1) + 
    !                             h^2/6*[(1-t)^3-(1-t)]*y2(i) + h^2/6*[t^3-t]*y2(i+1)  
    ! Expanding: f(t) = y(i) + [y(i+1)-y(i)]*t - h^2/6*y2(i)*[(1-t)^3-(1-t)] + 
    !                           h^2/6*y2(i+1)*[t^3-t]  
    ! Collect terms in powers of t:
        coeffs(1, i) = y(i)
        coeffs(2, i) = (y(i+1)-y(i)) - h*h*(2.0_dp*y2(i) + y2(i+1)) / 6.0_dp
        coeffs(3, i) = 0.5_dp * h*h * y2(i)
        coeffs(4, i) = h*h * (y2(i+1) - y2(i)) / 6.0_dp
      end if
    end do
    
    deallocate(y2)
  end subroutine cubic_coeffs


  !====
  ! Evaluate cubic interpolation at a single point using pre-computed coefficients
  !
  ! Inputs:
  !   x(n)           : strictly increasing array of sample points
  !   coeffs(4, n-1) : pre-computed coefficients from cubic_coeffs
  !   xv             : evaluation point
  !
  ! Output:
  !   fv : interpolated value at xv
  !====
  pure function cubic_eval(x, coeffs, xv) result(fv)
    real(dp), intent(in) :: x(:), coeffs(:,:), xv
    real(dp) :: fv
    
    integer :: n, i
    real(dp) :: t
    
    n = size(x)
    
    ! Validate inputs
    if (size(coeffs, 2) /= n-1 .or. size(coeffs, 1) /= 4 .or. n < 2) then
      fv = nan()
      return
    end if
    
    ! Find interval containing xv
    i = bracket_index(x, xv)
    
    ! Handle boundary cases
    if (i <= 0) then
      i = 1
    else if (i >= n) then
      i = n - 1
    end if
    
    ! Compute normalized position in interval
    t = (xv - x(i)) / (x(i+1) - x(i))
    
    ! Evaluate polynomial: a + b*t + c*t^2 + d*t^3
    fv = coeffs(1, i) + t * (coeffs(2, i) + t * (coeffs(3, i) + t * coeffs(4, i)))
    
  end function cubic_eval

  !====
  ! Add a weighted term to a sum when valid
  !
  ! Inputs:
  !   val            : Value to check and add if conditions met
  !   w              : Weight
  !
  ! Input/Outputs:
  !   fsum           : weighted sum
  !   wsum           : sum of weights
  !====
  subroutine add_w(val, w, fsum, wsum)
 
    real(dp), intent(in) :: val, w
    real(dp), intent(inout) :: fsum, wsum

    if (.not. ieee_is_nan(val) .and. w > 0.0_dp) then
      fsum = fsum + w*val 
      wsum = wsum + w
    end if  

  end subroutine add_w

  !=========================================================
  ! Cubic interpolation using Catmull-Rom spline
  ! Given 4 points f0, f1, f2, f3 and parameter t in [0,1]
  ! Interpolates between f1 and f2
  !=========================================================
  pure function cubic_interp(f0, f1, f2, f3, t, fv) result(f)
    real(dp), intent(in) :: f0, f1, f2, f3, t
    real(dp), intent(in) :: fv
    real(dp) :: f
    real(dp) :: t2, t3
    real(dp) :: a0, a1, a2, a3
  
    ! Check for NaN in any input
    if (ieee_is_nan(f0) .or. &
        ieee_is_nan(f1) .or. &
        ieee_is_nan(f2) .or. &
        ieee_is_nan(f3)) then
      f = fv
      return
    endif
  
    t2 = t * t
    t3 = t2 * t
  
    ! Catmull-Rom coefficients
    a0 = -0.5_dp*f0 + 1.5_dp*f1 - 1.5_dp*f2 + 0.5_dp*f3
    a1 =        f0 - 2.5_dp*f1 + 2.0_dp*f2 - 0.5_dp*f3
    a2 = -0.5_dp*f0            + 0.5_dp*f2
    a3 =                   f1
  
    f = a0*t3 + a1*t2 + a2*t + a3
  end function cubic_interp

  !==========================================================
  ! 1D interpolation for a single point
  !
  ! method options: 'nearest', 'linear' (default), 'spline',
  !                 'akima', 'cubic', 'pchip'
  ! fill_value: used when xv is out-of-bounds (default NaN)
  !
  ! Requirements:
  ! - x strictly increasing; y same length as x
  !==========================================================
  function interp1d_r(x, y, xv, method, fill_value, bounds_error, &
                     bc_ini, bc_end, yp_ini, yp_end) result(fv)
    real(dp), intent(in)                   :: x(:), y(:), xv
    character(len=*), intent(in), optional :: method
    real(dp), intent(in), optional         :: fill_value
    logical, intent(in), optional          :: bounds_error
    character(len=*), intent(in), optional :: bc_ini, bc_end
    real(dp), intent(in), optional         :: yp_ini, yp_end
    real(dp)                               :: fv
    logical error_out
    character(len=:), allocatable :: mthd
    character(len=:), allocatable :: bc_left, bc_right
    real(dp) :: yp_left, yp_right
    integer :: n, i
    integer :: i0, i1, i2, i3
    real(dp) :: t, yv
    real(dp), allocatable :: y2(:), dy(:), cub_coeffs(:,:)

    ! Default fill
    if (present(fill_value)) then
      yv = fill_value
    else
      yv = nan()
    end if

    n = size(x)
    if (size(y) /= n .or. n < 1) then
      fv = yv
      return
    end if

    ! Default method
    if (present(method)) then
      mthd = lowercase(ltrim(method))
    else
      mthd = 'linear'
    end if

    ! Boundary_conditions
    ! Default natural splines (ini=left, end=right)
    bc_left  = 'natural'; yp_left  = 0.0_dp
    bc_right = 'natural'; yp_right = 0.0_dp
    
    if (present(bc_ini)) then
      bc_left = ltrim(bc_ini)
      if (present(yp_ini)) yp_left = yp_ini
    endif

    if (present(bc_end)) then
      bc_right = ltrim(bc_end)
      if (present(yp_end)) yp_right = yp_end
    endif

    ! Validate increasing x
    if (n >= 2) then
      do i = 2, n
        if (x(i) <= x(i-1)) then
          fv = yv
          return
        end if
      end do
    end if

    ! Trivial n=1
    if (n == 1) then
      fv = y(1)
      return
    end if

    ! Out-of-bounds
    error_out = .false.
    if (present(bounds_error)) error_out = bounds_error

    if (xv < x(1) .or. xv > x(n)) then
      if (error_out) then
        call crash('interp1d_r - requested point out-of-bounds')
      else
        fv = xv
        return
      endif
    endif

    select case (mthd)
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
      fv = (1.0_dp - t)*y(i) + t*y(i+1)

    case ('spline')
      allocate(y2(n))
      call spline_coeffs_bc(x, y, y2, bc_left, bc_right, yp_left, yp_right)
      fv = spline_eval(x, y, y2, xv)
      deallocate(y2)

    case ('akima')
      allocate(dy(n))
      call akima_derivs(x, y, dy)
      fv = akima_eval(x, y, dy, xv) 
      deallocate(dy)

    case ('pchip')  
      allocate(dy(n))  
      call pchip_derivs(x, y, dy)  
      fv = pchip_eval(x, y, dy, xv)  
      deallocate(dy)  
  
    case('cubic')
!      Catmull-Rom interpolation
!      i = bracket_index(x, xv)
!      ! Clamp to valid range
!      if (i <= 0) then
!        i = 1
!      else if (i >= n) then
!        i = n-1
!      end if
!      ! Get 4-point stencil indices (with clamping at boundaries)
!      i0 = max(i-1, 1)
!      i1 = i
!      i2 = i+1
!      i3 = min(i+2, n)
!      ! Normalized coordinate within interval [0,1]
!      t = (xv - x(i)) / (x(i+1) - x(i))
!      ! Catmull-Rom interpolation
!      fv = cubic_interp(y(i0), y(i1), y(i2), y(i3), t, yv)

      ! ... Use splines
      allocate(cub_coeffs(4, n-1))   
      call cubic_coeffs(x, y, cub_coeffs)  
      fv = cubic_eval(x, cub_coeffs, xv)  
      deallocate(cub_coeffs) 

    case default
      call crash('interp1d_r - unknown interpolation method '//trim(mthd))

    end select
  end function interp1d_r

  !==========================================================
  ! Vectorized 1D interpolation for many query points
  ! method: 'nearest', 'linear' (default), 'spline', 'akima'
  ! fill_value applied for out-of-bounds
  ! Supports same boundary condition as interp1d_r
  !==========================================================
  function interp1d_vec(x, y, xq, method, fill_value, bounds_error, &
                        bc_ini, bc_end, yp_ini, yp_end) result(yq)
    real(dp), intent(in)                   :: x(:), y(:)
    real(dp), intent(in)                   :: xq(:)
    character(len=*), intent(in), optional :: method
    real(dp), intent(in), optional         :: fill_value
    logical, intent(in), optional          :: bounds_error
    character(len=*), intent(in), optional :: bc_ini, bc_end
    real(dp), intent(in), optional         :: yp_ini, yp_end
    real(dp), allocatable                  :: yq(:)

    logical error_out
    character(len=:), allocatable :: mthd
    character(len=:), allocatable :: bc_left, bc_right
    real(dp) :: yp_left, yp_right
    integer :: n, nq, i, j
    integer :: i0, i1, i2, i3
    real(dp) :: t, fv
    real(dp), allocatable :: y2(:), dy(:), cub_coeffs(:,:)

    n  = size(x)
    nq = size(xq)

    if (size(y).ne.n) call crash('interp1d_vec - incompatible dimensions')

    ! Validate x strictly increasing
    if (n >= 2) then
      do i = 2, n
        if (x(i) <= x(i-1)) call crash('interp1d_vec - no strictly increasing x')
      end do
    end if

    if (allocated(yq)) deallocate(yq)
    allocate(yq(nq))

    error_out = .false.
    if (present(bounds_error)) error_out = bounds_error

    ! Default method and fill
    if (present(method)) then
      mthd = lowercase(ltrim(method))
    else
      mthd = 'linear'
    end if

    if (present(fill_value)) then
      fv = fill_value
    else
      fv = nan()
    end if

    ! Boundary conditions (for spline method)  
    bc_left  = 'natural'; yp_left  = 0.0_dp  
    bc_right = 'natural'; yp_right = 0.0_dp  
    
    if (present(bc_ini)) then  
      bc_left = ltrim(bc_ini)  
      if (present(yp_ini)) yp_left = yp_ini  
    end if  
  
    if (present(bc_end)) then  
      bc_right = ltrim(bc_end)  
      if (present(yp_end)) yp_right = yp_end  
    end if  

    ! Trivial n cases
    if (n == 0) then
      if (nq > 0) yq = nan()
      return
    else if (n == 1) then
      yq = y(1)
      return
    end if

    ! Out-of-bounds
    if (error_out) then
      if (any(xq < x(1)) .or. any(xq > x(n))) then
        call crash('interp1d_vec - requested point out-of-bounds')
      endif
    end if

    select case (mthd)
    case ('nearest')
      do j = 1, nq
        if (xq(j).lt.x(1) .or. xq(j).gt.x(n)) then
          yq(j) = fv
        else
          i = bracket_index(x, xq(j))
          if (i <= 0) then
            yq(j) = y(1)
          else if (i >= n) then
            yq(j) = y(n)
          else
            ! Choose nearest between i and i+1
            if (abs(xq(j) - x(i)) <= abs(x(i+1) - xq(j))) then
              yq(j) = y(i)
            else
              yq(j) = y(i+1)
            endif
          endif
        endif
      enddo

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
            yq(j) = (1.0_dp - t)*y(i) + t*y(i+1)
          end if
        end if
      end do

    case ('spline')
      allocate(y2(n))
      call spline_coeffs_bc(x, y, y2, bc_left, bc_right, yp_left, yp_right)
      do j = 1, nq
        if (xq(j) < x(1) .or. xq(j) > x(n)) then
          yq(j) = fv
        else
          yq(j) = spline_eval(x, y, y2, xq(j))
        end if
      end do
      deallocate(y2)

    case ('akima')
      allocate(dy(n))  
      call akima_derivs(x, y, dy)  
      do j = 1, nq  
        if (xq(j) < x(1) .or. xq(j) > x(n)) then  
          yq(j) = fv  
        else  
          yq(j) = akima_eval(x, y, dy, xq(j))  
        end if  
      end do  
      deallocate(dy)  

    case ('pchip')  
      allocate(dy(n))  
      call pchip_derivs(x, y, dy)  
      do j = 1, nq  
        if (xq(j) < x(1) .or. xq(j) > x(n)) then  
          yq(j) = fv  
        else  
          yq(j) = pchip_eval(x, y, dy, xq(j))  
        end if  
      end do  
      deallocate(dy)  

    !case('cubic', 'catmull-rom')
      ! Catmull-rom
      !do j = 1, nq
      !  if (xq(j) < x(1) .or. xq(j) > x(n)) then
      !    yq(j) = fv
      !  else
      !    i = bracket_index(x, xq(j))
      !
      !    ! Clamp to valid range
      !    if (i <= 0) then
      !      i = 1
      !    else if (i >= n) then
      !      i = n-1
      !    end if
      !
      !    ! Get 4-point stencil
      !    i0 = max(i-1, 1)
      !    i1 = i
      !    i2 = i+1
      !    i3 = min(i+2, n)
      !
      !    ! Normalized coordinate
      !    t = (xq(j) - x(i)) / (x(i+1) - x(i))
      !
      !    ! Catmull-Rom interpolation
      !    yq(j) = cubic_interp(y(i0), y(i1), y(i2), y(i3), t, fv)
      !  endif
      !enddo

     case ('cubic')
       ! Pre-compute coefficients once for all evaluation points  
       allocate(cub_coeffs(4, n-1))  
       call cubic_coeffs(x, y, cub_coeffs)  
       do j = 1, nq  
         if (xq(j) < x(1) .or. xq(j) > x(n)) then  
           yq(j) = fv  
         else  
           yq(j) = cubic_eval(x, cub_coeffs, xq(j))  
         end if  
       end do  
       deallocate(cub_coeffs)  

    case default  
      call crash('interp1d_vec - unknown interpolation method '//trim(mthd))  
  
    end select  
  end function interp1d_vec


  !==========================================================
  ! 2D interpolation on a rectilinear grid
  !
  ! Inputs:
  ! - x(1:nx), y(1:ny): strictly increasing axes
  ! - z(nx, ny): function values at grid points (Fortran order)
  ! - xo, yo: query points (same length)
  !
  ! method: 'bilinear' (default), 'nearest', 'inverse_distance' = 'idw'
  !         'bicubic' = 'cubic'
  ! fill_value for out-of-bounds (default NaN)
  !==========================================================
  function interp2d(x, y, z, xo, yo, method, fill_value, bounds_error) result(zo)
    real(dp), intent(in)           :: x(:), y(:)
    real(dp), intent(in)           :: z(:,:)
    real(dp), intent(in)           :: xo, yo
    character(len=*), intent(in), optional :: method
    real(dp), intent(in), optional :: fill_value
    logical, intent(in), optional  :: bounds_error
    real(dp), allocatable          :: zo

    integer :: nx, ny
    integer :: i, j, ix, iy
    integer :: i0, i1, i2, i3, j0, j1, j2, j3  
    real(dp) :: fv, tx, ty
    real(dp) :: z00,z01,z10,z11
    real(dp) :: w00,w01,w10,w11
    real(dp) :: fy0, fy1, fy2, fy3  
    real(dp) :: dx,dy,dx2,dy2,d2,w2,eps2
    real(dp) :: fsum,wsum
    logical error_out
    character(len=:), allocatable :: mth

    nx = size(x); ny = size(y)

    if (size(z,1) /= nx .or. size(z,2) /= ny) call crash('interp2d - incompatible input dims')

    ! Defaults
    if (present(method)) then
      mth = lowercase(ltrim(method))
    else
      mth = 'bilinear'
    end if
    if (present(fill_value)) then
      fv = fill_value
    else
      fv = nan()
    end if

    ! Validate monotonicity
    if (nx < 2 .or. ny < 2) then
      zo = fv
      return
    end if
    do i = 2, nx
      if (x(i) <= x(i-1)) then
        zo = fv
        return
      end if
    end do
    do j = 2, ny
      if (y(j) <= y(j-1)) then
        zo = fv
        return
      end if
    end do

    ! Out-of-bounds
    error_out = .false.
    if (present(bounds_error)) error_out = bounds_error

    if (xo < x(1) .or. xo > x(nx) .or. yo < y(1) .or. yo > y(ny) ) then
      if (error_out) then
        call crash('interp2d - requested point out-of-bounds')
      else
        zo = fv
        return
      endif
    endif

    select case (mth)
    case ('nearest')
      ix = bracket_index(x, xo)
      iy = bracket_index(y, yo)
      ! resolve nearest in x
      if (ix <= 0) then
        ix = 1
      else if (ix >= nx) then
        ix = nx
      else
        if (abs(xo-x(ix)) > abs(x(ix+1)-xo)) ix = ix+1
      end if
      ! resolve nearest in y
      if (iy <= 0) then
        iy = 1
      else if (iy >= ny) then
        iy = ny
      else
        if (abs(yo-y(iy)) > abs(y(iy+1)-yo)) iy = iy+1
      end if
      zo = z(ix, iy)

    case('bilinear','linear')
      ix = bracket_index(x, xo)
      iy = bracket_index(y, yo)
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
      tx = (xo - x(ix)) / (x(ix+1) - x(ix))
      ty = (yo - y(iy)) / (y(iy+1) - y(iy))

      ! Function values
      z00 = z(ix  , iy  )
      z10 = z(ix+1, iy  )
      z01 = z(ix  , iy+1)
      z11 = z(ix+1, iy+1)

      ! Weigths
      w00 = (1.0_dp-tx) * (1.0_dp-ty)
      w10 =      tx     * (1.0_dp-ty)
      w01 = (1.0_dp-tx) *      ty      
      w11 =      tx     *      ty      

      ! Add valid values
      fsum = 0.0_dp; wsum = 0.0_dp
      call add_w(z00,w00,fsum,wsum)
      call add_w(z10,w10,fsum,wsum)
      call add_w(z01,w01,fsum,wsum)
      call add_w(z11,w11,fsum,wsum)

      if (wsum .gt. 0.0_dp) then
        zo = fsum / wsum
      else
        zo = fv
      endif 

    case('inverse_distance','idw')
      dx = (x(nx)-x(1))/real(nx-1,dp)
      dy = (y(ny)-y(1))/real(ny-1,dp)
      d2 = dx*dx + dy*dy
      eps2 = 0.0001d0*d2
      fsum = 0.0_dp; wsum = 0.0_dp
      do j=1,ny
        dy2 = (yo-y(j))**2
        do i=1,nx
          dx2 = (xo-x(i))**2
          d2  = dx2 + dy2
          if (d2 .le. eps2) then
            zo = z(i,j)
            return
          else
            w2 = 1.0_dp/d2
            call add_w(z(i,j),w2,fsum,wsum)
          endif
        enddo
      enddo
      if (wsum .gt. 0.0_dp) then
        zo = fsum / wsum
      else
        zo = fv
      endif

    case('bicubic','cubic')  
      ix = bracket_index(x, xo)  
      iy = bracket_index(y, yo)  
    
      ! Clamp to valid range for 4x4 stencil  
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
    
      ! Normalized coordinates within cell [0,1]  
      tx = (xo - x(ix)) / (x(ix+1) - x(ix))  
      ty = (yo - y(iy)) / (y(iy+1) - y(iy))  
    
      ! Get 4x4 stencil indices (with clamping at boundaries)  
      i0 = max(ix-1, 1)  
      i1 = ix  
      i2 = ix+1  
      i3 = min(ix+2, nx)  
      j0 = max(iy-1, 1)  
      j1 = iy  
      j2 = iy+1  
      j3 = min(iy+2, ny)

      ! Interpolate in x-direction for each of 4 y-rows using cubic  
      fy0 = cubic_interp(z(i0,j0), z(i1,j0), z(i2,j0), z(i3,j0), tx, fv)  
      fy1 = cubic_interp(z(i0,j1), z(i1,j1), z(i2,j1), z(i3,j1), tx, fv)  
      fy2 = cubic_interp(z(i0,j2), z(i1,j2), z(i2,j2), z(i3,j2), tx, fv)  
      fy3 = cubic_interp(z(i0,j3), z(i1,j3), z(i2,j3), z(i3,j3), tx, fv)  
    
      ! Interpolate in y-direction  
      zo = cubic_interp(fy0, fy1, fy2, fy3, ty, fv)  
    
      ! Handle NaN values - fallback to bilinear if any NaN in stencil  
      if (ieee_is_nan(zo)) then
        ! Fallback to bilinear with just the 4 corner points  
        z00 = z(i1, j1)  
        z10 = z(i2, j1)  
        z01 = z(i1, j2)  
        z11 = z(i2, j2)  
      
        w00 = (1.0_dp-tx) * (1.0_dp-ty)  
        w10 =      tx     * (1.0_dp-ty)  
        w01 = (1.0_dp-tx) *      ty        
        w11 =      tx     *      ty        
      
        fsum = 0.0_dp; wsum = 0.0_dp  
        call add_w(z00,w00,fsum,wsum)  
        call add_w(z10,w10,fsum,wsum)  
        call add_w(z01,w01,fsum,wsum)  
        call add_w(z11,w11,fsum,wsum)  
      
        if (wsum .gt. 0.0_dp) then  
          zo = fsum / wsum  
        else  
          zo = fv  
        endif  
      endif  

    case default
      call crash('interp2d - Invalid method '//trim(mth))

    end select
  end function interp2d

end module module_interp
