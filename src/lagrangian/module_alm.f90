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
! -------------------------------------------------------------------------!

module module_alm

use module_types

implicit none

! ... Dimension of the Lagrangian model
! ... Three-dimensional model: x, y, z
! ...
integer, parameter                               :: ndims = 3      ! Lagrangian model dimensions

! ... ALM beahvior: Model, fitting, assimilation
! ... By default we run the model
! ...
logical                                          :: option_model = .True.
logical                                          :: option_fitting = .False.
logical                                          :: option_assim = .False.

! ... Forward-Backward logical switch
! ...
logical                                          :: reverse = .False.

! ... Advection terms to activate:
! ...
logical                                          :: Uadv     = .False.
logical                                          :: Vadv     = .False.
logical                                          :: Wadv     = .False.
logical                                          :: Buoyancy = .False.

! ... Wind action activation:
! ...
logical                                          :: Winds = .False.

logical                                          :: Stationary  = .False.
logical                                          :: SingleLayer = .False.

! ... Verbose flag:
! ... 0 - Quiet mode
! ... 1 - Basic output (default) 
! ... 2 - Detailed output
! ... 3 - Maximum verbose level
! ...
integer                                          :: verb = 1       ! Verbose level

! ... The AML time units and calendar:
! ... They will coincide with the ones specified in the 
! ... ocean zonal velocity file or specified by the
! ... user.
! ...
character(len=180)                               :: alm_time_units = ''
character(len=20)                                :: alm_time_calendar = ''

real(dp)                                         :: alm_tini
real(dp)                                         :: alm_tfin
real(dp)                                         :: alm_tlen

! ... The AML 1D-grid axes and mask
! ...
integer                                          :: alm_nx = 1
integer                                          :: alm_ny = 1
integer                                          :: alm_nz = 1
real(dp)                                         :: alm_xmin
real(dp)                                         :: alm_xmax
real(dp)                                         :: alm_ymin
real(dp)                                         :: alm_ymax
real(dp)                                         :: alm_zmin
real(dp)                                         :: alm_zmax
real(dp)                                         :: alm_tmin
real(dp)                                         :: alm_tmax
real(dp)                                         :: alm_dx
real(dp)                                         :: alm_dy
real(dp), dimension(:), allocatable              :: alm_x
real(dp), dimension(:), allocatable              :: alm_y
real(dp), dimension(:), allocatable              :: alm_z
real(dp), dimension(:,:,:), allocatable          :: alm_mask

contains
  ! ...
  ! ==================================================================
  ! ==================================================================
  ! ...
  integer function point_type(xo,yo,zo) result(ptype)

    real(dp), intent(in)                         :: xo,yo,zo

    ! ... Local variables
    ! ...
    integer i,j,k
    real(dp) t,u,f

    ptype = -1
!    if (xo.le.alm_xmin) return
!    if (xo.ge.alm_xmax) return
!    if (yo.le.alm_ymin) return
!    if (yo.ge.alm_ymax) return

    if (xo.lt.alm_x(2)) return
    if (xo.gt.alm_x(alm_nx-1)) return
    if (yo.lt.alm_y(2)) return
    if (yo.gt.alm_y(alm_ny-1)) return

    ! ... On the vertical, the float is allowed to float AT the surface
    ! ... and stranded at the bottom.
    ! ... alm_zmax = 0
    ! ... alm_zmin = -DEPTH
    ! ...
    ptype = -1
    if (zo.gt.0.0D0) return   ! Allowed to float at the surface
    ptype = 0
    if (zo.le.alm_zmin) return   ! Stranded at the bottom

    i = (xo-alm_xmin)/alm_dx + 1
    j = (yo-alm_ymin)/alm_dy + 1

    if (alm_nz.gt.1) then
      do k=1,alm_nz
        if (alm_z(k).le.zo) exit
      enddo
    else
      k = 1
    endif

    t = (xo-alm_x(i))/alm_dx
    u = (yo-alm_y(j))/alm_dy
    f = (1.0D0-t)*(1.0D0-u)*alm_mask(i,j,k) + &
                t*(1.0D0-u)*alm_mask(i+1,j,k) + &
                        t*u*alm_mask(i+1,j+1,k) + &
                (1.0D0-t)*u*alm_mask(i,j+1,k)

    if (f.ge.0.5) then
      ptype = 1
    else
      ptype = 0
    endif

  end function point_type
  ! ...
  ! ==================================================================
  ! ...
end module module_alm

