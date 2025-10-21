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
! - interp1d
! - interp2d
! -------------------------------------------------------------------------!

module module_interp

use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
use module_types, only : dp
use module_constants, only : nan
use module_tools, only : lowercase, crash

implicit none (type, external)

private
public interp1d, interp2d

contains
! ...
! =====================================================================
! =====================================================================
! ...
  pure integer function lower_bound(arr, val) result(idx)  
    ! ... Returns largest i | arr(i) <= val for ascending arr  

    real(dp), intent(in) :: arr(:)  
    real(dp), intent(in) :: val  

    ! ... Local variables
    ! ...
    integer :: lo, hi, mid, n  

    n = size(arr)  
    lo = 1; hi = n  
    if (val < arr(1)) then  
      idx = 0  
      return  
    else if (val >= arr(n)) then  
      idx = n-1  
      return  
    end if  
    do  
      if (hi - lo <= 1) exit  
      mid = (lo + hi)/2  
      if (arr(mid) <= val) then  
        lo = mid  
      else  
        hi = mid  
      end if  
    end do  
    idx = lo  

  end function lower_bound  
  ! ...
  ! ====================================================================
  ! ...
  logical function is_monotonic_ascending(arr)  

    real(dp), intent(in) :: arr(:)  

    ! ... Local variables
    ! ...
    integer :: i  

    is_monotonic_ascending = .true.  
    do i = 1, size(arr)-1  
      if (arr(i+1) <= arr(i)) then  
        is_monotonic_ascending = .false.  
        return  
      end if  
    end do  

  end function is_monotonic_ascending  
  ! ...
  ! ====================================================================
  ! ...
  subroutine add_w(val, w, fsum, wsum)
 
    real(dp), intent(in) :: val, w
    real(dp), intent(inout) :: fsum, wsum

    if (.not. ieee_is_nan(val) .and. w > 0.0_dp) then
      fsum = fsum + w*val 
      wsum = wsum + w
    end if  

  end subroutine add_w
  ! ...
  ! ===================================================================
  ! ...
  function interp1d(f, x, xo, method, out_of_bounds) result(f_interp)

    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: f(size(x))
    real(dp), intent(in) :: xo
    character(len=*), intent(in) :: method
    character(len=*), intent(in), optional :: out_of_bounds
    real(dp) :: f_interp

    ! ... Local variables
    ! ...
    character(len=:), allocatable :: mth, oob
    integer :: n, i0
    real(dp) :: dx, wsum, fsum, d2, eps2

    n = size(x)
    if (present(out_of_bounds)) then
      oob = lowercase(adjustl(out_of_bounds))
    else
      oob = 'nan'
    end if
    mth = lowercase(adjustl(method))

    if (n < 1) then
      f_interp = nan(); return
    end if
    if (.not. is_monotonic_ascending(x)) then
      f_interp = nan(); return
    end if

    select case (mth)
    case ('nearest')
      if (xo <= x(1)) then
        if (oob == 'clamp' .or. oob == 'nearest') then
          f_interp = f(1)
        else
          f_interp = nan()
        end if
        return
      else if (xo >= x(n)) then
        if (oob == 'clamp' .or. oob == 'nearest') then
          f_interp = f(n)
        else
          f_interp = nan()
        end if
        return
      end if
      i0 = lower_bound(x, xo)
      if (i0 < n .and. (xo - x(i0)) > (x(i0+1) - xo)) i0 = i0 + 1
      if (i0 < 1) i0 = 1
      f_interp = f(i0)
      return

    case ('linear')
      i0 = lower_bound(x, xo)
      if (i0 < 1 .or. i0 >= n) then
        select case (oob)
        case ('clamp')
          if (xo <= x(1)) then
            f_interp = f(1); return
          else if (xo >= x(n)) then
            f_interp = f(n); return
          else
            f_interp = nan(); return
          end if
        case ('nearest')
          if (xo <= x(1)) then
            f_interp = f(1); return
          else if (xo >= x(n)) then
            f_interp = f(n); return
          else
            f_interp = nan(); return
          end if
        case default
          f_interp = nan(); return
        end select
      end if
      if (x(i0+1) == x(i0)) then
        f_interp = nan(); return
      end if
      dx = (xo - x(i0)) / (x(i0+1) - x(i0))
      if (ieee_is_nan(f(i0)) .and. ieee_is_nan(f(i0+1))) then
        f_interp = nan()
      else if (ieee_is_nan(f(i0))) then
        f_interp = f(i0+1)
      else if (ieee_is_nan(f(i0+1))) then
        f_interp = f(i0)
      else
        f_interp = (1.0_dp - dx)*f(i0) + dx*f(i0+1)
      end if
      return

    case ('inverse_distance','idw')
      fsum = 0.0_dp; wsum = 0.0_dp
      eps2 = 1.0e-12_dp
      do i0 = 1, n
        d2 = (x(i0) - xo)**2
        if (d2 <= eps2) then
          f_interp = f(i0); return
        end if
        if (.not. ieee_is_nan(f(i0))) then
          wsum = wsum + 1.0_dp/d2
          fsum = fsum + f(i0)/d2
        end if
      end do
      if (wsum > 0.0_dp) then
        f_interp = fsum/wsum
      else
        f_interp = nan()
      end if
      return

    case default
      call crash('In interp1d - invalid method '//trim(mth))

    end select

  end function interp1d
  ! ...
  ! ===================================================================
  ! ...
  function interp2d(f,x, y, xo, yo, method, out_of_bounds) result(f_interp)

    ! ... f(x(:), y(:)) : input function
    ! ... x(:), y(:)    : input grid
    ! ... xo, yo        : interpolation point
    ! ... method        : 'nearest', 'bilinear', 'inverse_distance'='idw'
    ! ... out_of_bounds : 'clamp', 'nearest', 'nan'

    real(dp), intent(in)                   :: x(:), y(:)
    real(dp), intent(in)                   :: f(size(x), size(y))
    real(dp), intent(in)                   :: xo, yo
    character(len=*), intent(in)           :: method
    character(len=*), intent(in), optional :: out_of_bounds
    real(dp)                               :: f_interp

    ! ... Local variables
    ! ... 
    logical ascx, ascy
    integer i, j, nx, ny, ix, iy
    real(dp) dx, dy, wsum, fsum, d2, eps2
    real(dp) f00, f10, f01, f11
    real(dp) w00, w10, w01, w11
    character(len=:), allocatable  :: mth, oob

    nx = size(x); ny = size(y)

    ! ... Interpolation method and out of bounds
    ! ...
    mth = lowercase(adjustl(method))
    if (present(out_of_bounds)) then
      oob = lowercase(adjustl(out_of_bounds))
    else
      oob = 'nan'
    endif

    ! ... Basic input validation
    ! ...
    if (nx.lt.1.or.ny.lt.1) then
      f_interp = nan()
      return
    endif
    ascx = is_monotonic_ascending(x)
    ascy = is_monotonic_ascending(y)
    if (.not.ascx.or..not.ascy) then
      write(*,*) 'Warning in interp2d: not ascending unique grid'
      f_interp = nan()
      return
    endif

    select case (mth)
    case ('nearest')
      ! Use lower_bound and clamp based on oob policy
      if (xo < x(1)) then
        if (oob == 'clamp' .or. oob == 'nearest') then
          ix = 1
        else
          f_interp = nan(); return
        end if
      else if (xo > x(nx)) then
        if (oob == 'clamp' .or. oob == 'nearest') then
          ix = nx
        else
          f_interp = nan(); return
        end if
      else
        ix = lower_bound(x, xo)
        if (ix < nx .and. (xo - x(ix)) > (x(ix+1) - xo)) ix = ix + 1
        if (ix == 0) ix = 1
      end if

      if (yo < y(1)) then
        if (oob == 'clamp' .or. oob == 'nearest') then
          iy = 1
        else
          f_interp = nan(); return
        end if
      else if (yo > y(ny)) then
        if (oob == 'clamp' .or. oob == 'nearest') then
          iy = ny
        else
          f_interp = nan(); return
        end if
      else
        iy = lower_bound(y, yo)
        if (iy < ny .and. (yo - y(iy)) > (y(iy+1) - yo)) iy = iy + 1
        if (iy == 0) iy = 1
      end if
      f_interp = f(ix, iy)
      return

    case ('bilinear','linear')
      ! ... Bilinear interpolation between four grid corners:
      ! ...  f(x,y) = (1−dx)*(1−dy)*f00 + dx*(1−dy)*f10 + (1−dx)*dy*f01 + dx*dy*f11
      ! ...
      ix = lower_bound(x, xo)
      iy = lower_bound(y, yo)
      if (ix < 1 .or. iy < 1 .or. ix >= nx .or. iy >= ny) then
        select case (oob)
        case ('clamp')
          ix = max(1, min(nx-1, ix)); iy = max(1, min(ny-1, iy))
        case ('nearest')
          ! Project to nearest valid cell:
          if (xo <= x(1)) ix = 1
          if (xo >= x(nx)) ix = nx-1
          if (yo <= y(1)) iy = 1
          if (yo >= y(ny)) iy = ny-1
          if (ix < 1 .or. iy < 1 .or. ix >= nx .or. iy >= ny) then
            f_interp = nan(); return
          end if
        case default
          f_interp = nan(); return
        end select
      end if
      if (x(ix+1) == x(ix) .or. y(iy+1) == y(iy)) then
        f_interp = nan(); return
      end if

      dx = (xo - x(ix)) / (x(ix+1) - x(ix))
      dy = (yo - y(iy)) / (y(iy+1) - y(iy))

      ! ... Function values
      ! ...
      f00 = f(ix  , iy  )
      f10 = f(ix+1, iy  )
      f01 = f(ix  , iy+1)
      f11 = f(ix+1, iy+1)

      ! ... Weights:
      ! ...
      w00 = (1.0_dp - dx) * (1.0_dp - dy)
      w10 =       dx      * (1.0_dp - dy)
      w01 = (1.0_dp - dx) *       dy
      w11 =       dx      *       dy

      if (all(.not. ieee_is_nan([f00, f10, f01, f11]))) then
        ! ... All finite values
        ! ...
        f_interp = w00*f00 + w10*f10 + w01*f01 + w11*f11
      else
        ! ... Graceful degradation
        ! ...
        fsum = 0.0_dp; wsum = 0.0_dp
        call add_w(f00, w00, fsum, wsum)
        call add_w(f10, w10, fsum, wsum)
        call add_w(f01, w01, fsum, wsum)
        call add_w(f11, w11, fsum, wsum)

        if (wsum > 0.0_dp) then
          f_interp = fsum / wsum
        else
          f_interp = nan()
        end if
      end if
      return

    case ('inverse_distance','idw')
      fsum = 0.0_dp; wsum = 0.0_dp
      eps2 = (1.0e-12_dp)
      do i = 1, nx
      do j = 1, ny
        d2 = (x(i) - xo)**2 + (y(j) - yo)**2
        if (d2 <= eps2) then
          f_interp = f(i, j); return
        end if
        if (.not. ieee_is_nan(f(i,j))) then
          wsum = wsum + 1.0_dp/d2
          fsum = fsum + f(i,j)/d2
        end if
      end do
      end do
      if (wsum > 0.0_dp) then
        f_interp = fsum/wsum
      else
        f_interp = nan()
      end if
      return

    case default
      call crash('In interp2d - invalid method '//trim(mth))

    end select

  end function interp2d
  ! ...
  ! ===================================================================
  ! ...

end module module_interp
