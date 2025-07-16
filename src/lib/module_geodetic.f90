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
! - haversine
! - bearing
! - compass_to_polar
! - polar_to_compass
! -------------------------------------------------------------------------!

module module_geodetic

use, intrinsic :: IEEE_ARITHMETIC, ONLY : IEEE_VALUE, IEEE_QUIET_NAN
use module_types
use module_constants

implicit none

type type_WGS84_constants
  real(dp)             :: Earth_Radius     = 6371000.0_dp                 ! Mean Earth radius (spherical equivalent) [m]
  !real(dp)             :: Earth_Radius     = 6371008.8_dp
  real(dp)             :: Earth_SemiMajor  = 6378137.0_dp                 ! a, WGS84 semi-major axis [m]
  real(dp)             :: Earth_InvFlatten = 298.257223562997_dp          ! 1/f, WGS84 inverse flattening
 end type type_WGS84_constants
 
 type(type_WGS84_constants)                          :: WGS84

contains
! ...
! =====================================================================
! =====================================================================
! ...
  real(dp) function haversine (lon1,lat1,lon2,lat2)
    ! ...
    ! ... Function Haversine
    ! ... Determines the great-circle distance between two points in a
    ! ... sphere. The input lngitudes and latitudes are given in radians
    ! ... The retruned distance is in meters.
    ! ... Earth_Radius = 6371315.0_dp      ! m
    ! ...
    real(dp), intent(in)                  :: lon1,lat1
    real(dp), intent(in)                  :: lon2,lat2

    ! ... Local variables
    ! ...
    real(dp) sindx,sindy,dang

    sindx = sin(0.5D0*(lon2-lon1))
    sindy = sin(0.5D0*(lat2-lat1))

    dang = 2.0d0*asin(sqrt(sindy*sindy + cos(lat2)*cos(lat1)*sindx*sindx))
    haversine = constants%Earth_Radius * dang

    return
  end function haversine
  ! ...
  ! ===================================================================
  ! ...
  subroutine bearing (lon1,lat1,lon2,lat2,azimuth,rev_azimuth)
    ! ...
    ! ... Calculation of the azimuth and reverse azimuth from (lon1,lat1) to (lon2,lat2)
    ! ... The input lngitudes and latitudes are given in radians

    real(dp), intent(in)                  :: lon1,lat1
    real(dp), intent(in)                  :: lon2,lat2
    real(dp), intent(out)                 :: azimuth,rev_azimuth

    ! ... Local variables
    ! ...
    real(dp) dlon,cdlon,sdlon,x,y,x_rev,y_rev

    dlon = lon2 - lon1
    cdlon = cos(dlon)
    sdlon = sin(dlon)
 
    ! Azimuth
    y = sdlon * cos(lat2)
    x = cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cdlon
    azimuth = mod(atan2(y, x)*rad2deg+360.0d0, 360.0d0)
    
    ! Reverse azimuth
    y_rev = -sdlon * cos(lat1)
    x_rev = cos(lat2)*sin(lat1) - sin(lat2)*cos(lat1)*cdlon
    rev_azimuth = mod(atan2(y_rev, x_rev)*rad2deg+360.0d0, 360.0d0)

    return
  end subroutine bearing
  ! ...
  ! ===================================================================
  ! ...
  real(dp) pure function compass_to_polar(compass)
    ! ... Transform from compass angle (clockwise starting at North) to 
    ! ... polar (counter clockwise starting at east).
    ! ...
    real(dp), intent(in)                  :: compass

    compass_to_polar = mod(450.0 - compass, 360.0)

  end function compass_to_polar
  ! ...
  ! ===================================================================
  ! ...
  real(dp) pure function polar_to_compass(polar)
    !=======================================================================
    ! ... Transform from polar angle (counter clockwise starting at east) to 
    ! ... compass (clockwise starting at north).
    !=======================================================================

    real(dp), intent(in)                  :: polar

    polar_to_compass = mod(450.0 - polar, 360.0)

  end function polar_to_compass
  ! ...
  ! ===================================================================
  ! ...
  subroutine gc_destination_point(lat0_deg, lon0_deg, bearing_deg, dist, lat1_deg, lon1_deg)
  !=======================================================================
  ! ... Compute destination point given start point, bearing, distance
  !=======================================================================

    real(dp), intent(in)  :: lat0_deg, lon0_deg, bearing_deg, dist
    real(dp), intent(out) :: lat1_deg, lon1_deg

    ! ... Local variables
    ! ...
    real(dp) :: lat0, lon0, bearing
    real(dp) :: lat1, lon1

    ! ... Convert input to radians
    ! ...
    lat0 = lat0_deg * deg2rad
    lon0 = lon0_deg * deg2rad
    bearing = bearing_deg * deg2rad

    ! ... Compute new latitude
    ! ...
    lat1 = asin(sin(lat0) * cos(dist / WGS84%Earth_Radius) + cos(lat0) * sin(dist / WGS84%Earth_Radius) * cos(bearing))

    ! ... Compute new longitude
    ! ...
    lon1 = lon0 + atan2(sin(bearing) * sin(dist / WGS84%Earth_Radius) * cos(lat0), &
                        cos(dist / WGS84%Earth_Radius) - sin(lat0) * sin(lat1))

    ! ... Normalize longitude to [-π, π]
    ! ...
    lon1 = modulo(lon1 + pi, 2.0d0 * pi) - pi

    ! ... Convert to degrees
    ! ...
    lat1_deg = lat1 * rad2deg
    lon1_deg = lon1 * rad2deg

  end subroutine gc_destination_point
  ! ...
  ! ===================================================================
  ! ...
  subroutine gc_inverse(lat1_deg, lon1_deg, lat2_deg, lon2_deg, dist, bearing_deg)
  !=======================================================================
  ! ... Compute initial bearing and great-circle distance between two points
  !=======================================================================

    real(dp), intent(in)  :: lat1_deg, lon1_deg, lat2_deg, lon2_deg
    real(dp), intent(out) :: dist, bearing_deg

    real(dp) :: lat1, lon1, lat2, lon2
    real(dp) :: delta_lon, x, y, c, theta

    ! ... Convert to radians
    ! ...
    lat1 = lat1_deg * deg2rad
    lon1 = lon1_deg * deg2rad
    lat2 = lat2_deg * deg2rad
    lon2 = lon2_deg * deg2rad

    delta_lon = lon2 - lon1

    ! ... Compute central angle using haversine
    ! ...
    c = acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(delta_lon))
    dist = WGS84%Earth_Radius * c

    ! Compute initial bearing
    y = sin(delta_lon) * cos(lat2)
    x = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(delta_lon)
    theta = atan2(y, x)

    ! Convert to degrees and normalize to [0, 360)
    bearing_deg = modulo(theta * rad2deg + 360.0d0, 360.0d0)

  end subroutine gc_inverse
  ! ...
  ! ===================================================================
  ! ...
end module module_geodetic
