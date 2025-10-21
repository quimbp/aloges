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
! - interpol2d
! -------------------------------------------------------------------------!

module module_interp

use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
use module_types, only : dp
use module_constants

implicit none (type, external)

contains
! ...
! =====================================================================
! =====================================================================
! ...
  function interpol2d(f,x, y,xo, yo, method) result(f_interp)

    real(dp), intent(in) :: x(:), y(:)
    real(dp), intent(in) :: f(size(x), size(y))
    real(dp), intent(in) :: xo, yo
    character(len=*), intent(in) :: method
    real(dp) :: f_interp

    ! ... Local variables
    ! ... 
    logical found
    integer i, j, nx, ny, ix, iy
    real(dp) dx, dy, d, wsum, fsum, dist2, dx2

    nx = size(x)
    ny = size(y)

    ! --- Nearest Neighbor ---
    if (trim(method) == 'nearest') then
        f_interp = 1.0D30
        dist2 = 1.0D30
        do i = 1, nx
            dx2 = (x(i) - xo)**2
            do j = 1, ny
                d = dx2 + (y(j) - yo)**2
                if (d < dist2) then
                    dist2 = d
                    f_interp = f(i, j)
                end if
            end do
        end do
        return
    end if

    ! --- Bilinear Interpolation ---
    if (trim(method) == 'bilinear') then
        found = .false.
        do i = 1, nx-1
            if (xo >= x(i) .and. xo <= x(i+1)) then
                ix = i
                do j = 1, ny-1
                    if (yo >= y(j) .and. yo <= y(j+1)) then
                        iy = j
                        found = .true.
                        exit
                    end if
                end do
                if (found) exit
            end if
        end do

        if (.not. found) then
            f_interp = 1.0D30  
            return
        end if

        dx = (xo - x(ix)) / (x(ix+1) - x(ix))
        dy = (yo - y(iy)) / (y(iy+1) - y(iy))

        f_interp = (1-dx)*(1-dy)*f(ix,iy)   + dx*(1-dy)*f(ix+1,iy) + &
                   (1-dx)*dy*f(ix,iy+1)     + dx*dy*f(ix+1,iy+1)
        return
    end if

    ! --- Inverse Distance Weighting (quadratic decay) ---
    if (trim(method) == 'inverse_distance') then
        fsum = 0.0D0
        wsum = 0.0D0
        do i = 1, nx
            do j = 1, ny
                d = sqrt((x(i) - xo)**2 + (y(j) - yo)**2)
                if (d < 1.0e-6) then
                    f_interp = f(i, j)
                    return
                end if
                wsum = wsum + 1.0 / (d*d)
                fsum = fsum + f(i,j) / (d*d)
            end do
        end do
        if (wsum > 0.0) then
            f_interp = fsum / wsum
        else
            f_interp = 1.0D30
        end if
        return
    end if

    ! Unknown method
    f_interp = 1.0D30

  end function interpol2d
  ! ...
  ! ===================================================================
  ! ...
end module module_interp
