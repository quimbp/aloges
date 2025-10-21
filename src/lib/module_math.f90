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
! - randstr                                                                !
! - randseries                                                             !
! -------------------------------------------------------------------------!

module module_math

use, intrinsic :: IEEE_ARITHMETIC, ONLY : IEEE_VALUE, IEEE_QUIET_NAN
use module_types
use module_constants

implicit none (type, external)

private
public :: randn, randstr, randseries, arange, mean, variance, &
          indexx, swap, imaxloc, median, corrcoef, &
          brent, golden, spline, quicksort, linspace, meshgrid, &
          diff, gradient, cumsum, cumprod, percentile, interp1

interface randn
  module procedure randn_r,randn_v,randn_a
end interface randn

interface arange
  module procedure arange_i
  module procedure arange_dp
end interface arange

interface swap
  module procedure swap_r
  module procedure swap_v
end interface swap

interface interp1
  module procedure interp1_r, interp1_v
end interface interp1

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
    real(dp), dimension(:,:), allocatable        :: ff

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

    if (.not.allocated(ff)) allocate(ff(m,n))

    !$omp parallel do default(shared) private(i,j,u,v,x,y,q) collapse(2)
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
  function randstr(len,iseed) result(name)

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
  end function randstr
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
  pure function arange_i(amin, amax, adelta) result(arange)
    ! ...
    ! ... Returns an integer-valued vector from amin to amax with 
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
    real(dp) Sw

    mean = nan()          ! Defined in module_constants.f90 
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
  real(dp) function variance (A,W,biased)

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
    logical, intent(in), optional          :: biased

    logicaL                                    :: weight = .false.
    logical                                    :: isbiased = .false.
    integer N,i
    real(dp) xsum1,xsum2,ai,wi,Sw


    variance = nan()          ! Defined in module_constants.f90
    N = SIZE(A)
    if (N.EQ.0) return

    if (present(W)) then
      if (size(W).ne.N) stop 'ERROR in variance: Incompatible weight size'
      weight = .true.
    else
      weight = .false.
    endif

    if (present(biased)) then
      isbiased = biased
    else
      isbiased = .false.
    endif
 
    !if (PRESENT(W)) then
    !  if (SIZE(W).EQ.1) then
    !    if (W(1).EQ.0) then
    !      isbiased = .false.
    !    else if (W(1).EQ.1) then
    !      isbiased = .true.
    !    else
    !      STOP 'ERROR in variance: Invalid normalization. Use 0 (unbiased) or 1 (biased)'
    !    endif
    !  else if (SIZE(W).EQ.N) then
    !    weight = .true.
    !  else
    !    STOP 'ERROR in variance: Size W must be 1 or N'
    !  endif
    !endif

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
      variance = Sw * xsum2 /(Sw - 1.0_dp)
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
      IF (isbiased) then
        variance = xsum2
      else
        variance = N*xsum2/(N-1.0_dp)
      endif
    endif
    IF (variance.LT.0) variance = 0.0_dp

  end function variance
  ! ...
  ! ===================================================================
  ! ...
  subroutine indexx(arr,indx)

  ! ... Subroutine indexx from Numerical Recipes in Fortran.
  ! ... The Art of Scientific Computing. Press et al., 1992.
  ! ...

    real(dp), dimension(:), intent(in)          :: arr
    integer, dimension(size(arr)), intent(out)  :: indx

    ! ... Local varibles:
    ! ...
    integer, parameter                          :: M = 7
    integer                                     :: n,i,indxt,r,itemp
    integer                                     :: j,jstack,k,l
    integer, dimension(size(arr))               :: istack
    real(dp)                                    :: a

    n  =  size(arr)

    do j = 1,n
      indx(j) = j
    enddo

    jstack = 0
    l = 1
    r = n
    do 
      if (r-l.lt.M) then
        do j = l+1,r
          indxt = indx(j)
          a = arr(indxt)
          do i = j-1,1,-1
            if (arr(indx(i)).le.a) exit
            indx(i+1) = indx(i)
          enddo
          indx(i+1) = indxt
        enddo
        if (jstack.eq.0) return
        r = istack(jstack)
        l = istack(jstack-1)
        jstack = jstack-2
      else
        k = (l+r)/2
        itemp = indx(k)
        indx(k) = indx(l+1)
        indx(l+1) = itemp
        if (arr(indx(l+1)).gt.arr(indx(r))) then
          itemp     = indx(l+1)
          indx(l+1) = indx(r)
          indx(r)   = itemp
        endif
        if (arr(indx(l)).gt.arr(indx(r))) then
          itemp   = indx(l)
          indx(l) = indx(r)
          indx(r) = itemp
        endif
        if (arr(indx(l+1)).gt.arr(indx(l))) then
          itemp     = indx(l+1)
          indx(l+1) = indx(l)
          indx(l)   = itemp
        endif
        i = l + 1
        j = r
        indxt = indx(l+1)
        a = arr(indxt)
        do
          do
            i = i + 1
            if (arr(indx(i)).ge.a) exit
          enddo
          do 
            j = j - 1
            if (arr(indx(j)).le.a) exit
          enddo
          if (j.lt.i) exit
          itemp   = indx(i)
          indx(i) = indx(j)
          indx(j) = itemp
        enddo
        indx(l+1) = indx(j)
        indx(j)   = indxt
        jstack    = jstack + 2
        if (r-i+1.ge.j-l) then
          istack(jstack) = r
          istack(jstack-1) = i
          r = j - 1
        else
          istack(jstack) = j-1
          istack(jstack-1) = l
          l = i
        endif
      endif
    enddo

    end subroutine indexx
  ! ...
  ! ===================================================================
  ! ...
  subroutine swap_r(a,b)
    real(dp), intent(inout)  :: a,b
    real(dp)                 :: c
    c = a
    a = b
    b = c
  end subroutine swap_r
  ! ...
  ! ===================================================================
  ! ...
  subroutine swap_v(a,b)
    real(dp), dimension(:), intent(inout)  :: a,b
    real(dp), dimension(size(a)) :: c
    c(:) = a(:)
    a(:) = b(:)
    b(:) = c(:)
  end subroutine swap_v
  ! ...
  ! ===================================================================
  ! ...
  integer pure function imaxloc(a)

    real(8), dimension(:), intent(in)    :: a

    integer imax(1)

    imax = maxloc(a(:))
    imaxloc = imax(1)
  end function imaxloc
  ! ...
  ! ===================================================================
  ! ...
  function median(x) result(med)

    real(dp), dimension(:), intent(in) :: x
    real(dp) :: med

    ! ... Local variables
    ! ...
    real(dp), dimension(size(x)) :: x_sorted
    integer :: n

    n = size(x)
    x_sorted = x
    call quicksort(x_sorted)

    if (mod(n,2) == 0) then
      med = (x_sorted(n/2) + x_sorted(n/2+1)) / 2.0_dp
    else
      med = x_sorted((n+1)/2)
    endif
  end function median
  ! ...
  ! ===================================================================
  ! ...
! Replace the quicksort subroutine in module_math.f90 with this version
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
        ! Choose the rightmost element as pivot
        pivot = arr(high)
        i = low - 1
      
        do j = low, high - 1
          if (arr(j) <= pivot) then
            i = i + 1
            ! Swap arr(i) and arr(j)
            temp = arr(i)
            arr(i) = arr(j)
            arr(j) = temp
          end if
        end do
      
        ! Swap arr(i+1) and arr(high) (pivot)
        temp = arr(i+1)
        arr(i+1) = arr(high)
        arr(high) = temp
      
        ! Recursively sort elements before and after partition
        call quicksort_helper(arr, low, i)
        call quicksort_helper(arr, i+2, high)
    end if
  end subroutine quicksort_helper

end subroutine quicksort
  ! ...
  ! ===================================================================
  ! ...
  function corrcoef(x, y) result(r)
    real(dp), dimension(:), intent(in) :: x, y
    real(dp)                           :: r

    ! ... Local variables
    ! ...
    real(dp) :: mx, my, sx, sy, sxy
    integer :: n

    n = size(x)
    if (size(y) /= n) then
      r = ieee_value(1.0_dp, ieee_quiet_nan)
      return
    endif

    mx = sum(x)/n
    my = sum(y)/n
    sx = sqrt(sum((x - mx)**2)/(n - 1))
    sy = sqrt(sum((y - my)**2)/(n - 1))
    sxy = sum((x - mx)*(y - my))/(n - 1)

    if (sx > 0.0_dp .and. sy > 0.0_dp) then
      r = sxy/(sx*sy)
    else
      r = 0.0_dp
    endif

  end function corrcoef
  ! ...
  ! ===================================================================
  ! ...
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
    
    a  = aa
    b  = bb

    fa = f(a)
    fb = f(b)
    
    if (fa*fb > 0.0_dp) then
      root = ieee_value(1.0_dp, ieee_quiet_nan)
      return
    endif
    
    c = b
    fc = fb
    
    do iter = 1, max_iter
      if (fb*fc > 0.0_dp) then
        c = a
        fc = fa
        d = b - a
        e = d
      endif
      
      if (abs(fc) < abs(fb)) then
        a = b
        b = c
        c = a
        fa = fb
        fb = fc
        fc = fa
      endif
      
      tol1 = 2.0_dp*eps*abs(b) + 0.5_dp*tol
      xm = 0.5_dp*(c - b)
      
      if (abs(xm) <= tol1 .or. fb == 0.0_dp) then
        root = b
        return
      endif
      
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
        endif
        
        if (p > 0.0_dp) q = -q
        p = abs(p)
        
        if (2.0_dp*p < min(3.0_dp*xm*q - abs(tol1*q), abs(e*q))) then
          e = d
          d = p/q
        else
          d = xm
          e = d
        endif
      else
        d = xm
        e = d
      endif
      
      a = b
      fa = fb
      
      if (abs(d) > tol1) then
        b = b + d
      else
        b = b + sign(tol1, xm)
      endif
      
      fb = f(b)
    end do
    
    root = b
  end function brent
  ! ...
  ! ===================================================================
  ! ...
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
      endif
      
      c = b - (b - a)/phi
      d = a + (b - a)/phi
    end do
    
    xmin = (a + b)/2.0_dp

  end function golden
  ! ...
  ! ===================================================================
  ! ...
  function interp1_r(x, y, xi) result(yi)
    real(dp), dimension(:), intent(in) :: x, y
    real(dp), intent(in) :: xi
    real(dp) :: yi
    
    integer :: i, n
    real(dp) :: dx
    
    n = size(x)
    if (size(y) /= n) then
      yi = ieee_value(1.0_dp, ieee_quiet_nan)
      return
    endif
    
    if (xi <= x(1)) then
      yi = y(1)
    else if (xi >= x(n)) then
      yi = y(n)
    else
      do i = 1, n - 1
        if (xi >= x(i) .and. xi <= x(i + 1)) then
          dx = x(i + 1) - x(i)
          if (dx > 0.0_dp) then
            yi = y(i) + (y(i + 1) - y(i))*(xi - x(i))/dx
          else
            yi = (y(i) + y(i + 1))/2.0_dp
          endif
          exit
        endif
      end do
    endif
  end function interp1_r
  ! ...
  ! ===================================================================
  ! ...
  function interp1_v(x, y, xi) result(yi)
    real(dp), dimension(:), intent(in) :: x, y
    real(dp), dimension(:), intent(in) :: xi
    real(dp), dimension(size(xi)) :: yi
    
    integer :: i
    
    do i = 1, size(xi)
      yi(i) = interp1_r(x, y, xi(i))
    end do
  end function interp1_v
  ! ...
  ! ===================================================================
  ! ...
  subroutine spline(x, y, y2)
    real(dp), dimension(:), intent(in) :: x, y
    real(dp), dimension(:), intent(out) :: y2
    
    integer :: i, n
    real(dp) :: p, sig
    real(dp), dimension(:), allocatable :: u
    
    n = size(x)
    if (size(y) /= n .or. size(y2) /= n) then
      error stop "spline: array size mismatch"
    endif
    
    allocate(u(n))
    
    y2(1) = 0.0_dp
    u(1) = 0.0_dp
    
    do i = 2, n - 1
      sig = (x(i) - x(i - 1))/(x(i + 1) - x(i - 1))
      p = sig*y2(i - 1) + 2.0_dp
      y2(i) = (sig - 1.0_dp)/p
      u(i) = (6.0_dp*((y(i + 1) - y(i))/(x(i + 1) - x(i)) - &
              (y(i) - y(i - 1))/(x(i) - x(i - 1)))/(x(i + 1) - x(i - 1)) - &
              sig*u(i - 1))/p
    end do
    
    y2(n) = 0.0_dp
    
    do i = n - 1, 1, -1
      y2(i) = y2(i)*y2(i + 1) + u(i)
    end do
    
    deallocate(u)
  end subroutine spline
  ! ...
  ! ===================================================================
  ! ...
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
    endif
  end function linspace
  ! ...
  ! ===================================================================
  ! ...
  subroutine meshgrid(x, y, X_mat, Y_mat)
    real(dp), dimension(:), intent(in) :: x, y
    real(dp), dimension(:,:), intent(out) :: X_mat, Y_mat
    
    integer :: i, j, m, n
    
    m = size(x)
    n = size(y)
    
    if (size(X_mat,1) /= m .or. size(X_mat,2) /= n .or. &
        size(Y_mat,1) /= m .or. size(Y_mat,2) /= n) then
      error stop "meshgrid: output array size mismatch"
    endif
    
    do concurrent (i = 1:m, j = 1:n)
      X_mat(i,j) = x(i)
      Y_mat(i,j) = y(j)
    end do
  end subroutine meshgrid
  ! ...
  ! ===================================================================
  ! ...
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
    endif
    
    allocate(dx(m - diff_order))
    dx = x
    
    do i = 1, diff_order
      dx = dx(2:) - dx(:size(dx)-1)
    end do
  end function diff
  ! ...
  ! ===================================================================
  ! ...
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
    endif
  end function gradient
  ! ...
  ! ===================================================================
  ! ...
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
  ! ...
  ! ===================================================================
  ! ...
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
  ! ...
  ! ===================================================================
  ! ...
  function percentile(x, p) result(perc)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(in) :: p
    real(dp) :: perc
    
    real(dp), dimension(:), allocatable :: x_sorted
    integer :: n, index
    real(dp) :: frac
    
    if (p < 0.0_dp .or. p > 100.0_dp) then
      perc = ieee_value(1.0_dp, ieee_quiet_nan)
      return
    endif
    
    n = size(x)
    allocate(x_sorted(n))
    x_sorted = x
    call quicksort(x_sorted)
    
    frac = p/100.0_dp * (n - 1) + 1
    index = floor(frac)
    
    if (index < 1) then
      perc = x_sorted(1)
    else if (index >= n) then
      perc = x_sorted(n)
    else
      perc = x_sorted(index) + (frac - index)*(x_sorted(index+1) - x_sorted(index))
    endif
    
    deallocate(x_sorted)

  end function percentile
  ! ...
  ! ===================================================================
  ! ...
end module module_math
