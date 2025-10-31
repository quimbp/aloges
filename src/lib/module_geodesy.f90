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

module module_geodesy

use, intrinsic :: IEEE_ARITHMETIC, ONLY : IEEE_VALUE, IEEE_QUIET_NAN
use module_types
use module_constants

implicit none (type, external)

type type_WGS84_constants
  !real(dp)             :: Earth_Radius     = 6371000.0_dp                 ! Mean Earth radius (spherical equivalent) [m]
  real(dp)             :: Earth_Radius     = 6371008.8_dp                 ! Mean Earth radius (spherical equivalent) [m]
  real(dp)             :: Earth_SemiMajor  = 6378137.0_dp                 ! a, WGS84 semi-major axis [m]
  real(dp)             :: Earth_InvFlatten = 298.257223562997_dp          ! 1/f, WGS84 inverse flattening
 end type type_WGS84_constants

 type(type_WGS84_constants)                          :: WGS84

interface spdir2uv 
  module procedure spdir2uv_r, spdir2uv_v
end interface spdir2uv

interface uv2spdir
  module procedure uv2spdir_r, uv2spdir_v
end interface uv2spdir

contains
! ...
! =====================================================================
! =====================================================================
! ...
  real(dp) function haversine (lon1,lat1,lon2,lat2)
    !==================================================================
    ! ... Function Haversine
    ! ... Determines the great-circle distance between two points in a
    ! ... sphere. The input lngitudes and latitudes are given in radians.
    ! ... The retruned distance is in meters.
    ! ... Earth_Radius = 6371315.0_dp      ! m
    !==================================================================
    
    real(dp), intent(in)                  :: lon1,lat1         ! Radians
    real(dp), intent(in)                  :: lon2,lat2         ! Radians

    ! ... Local variables
    ! ...
    real(dp) sindx,sindy,dang

    sindx = sin(0.5D0*(lon2-lon1))
    sindy = sin(0.5D0*(lat2-lat1))

    dang = 2.0d0*asin(sqrt(sindy*sindy + cos(lat2)*cos(lat1)*sindx*sindx))
    haversine = constants%Earth_Radius * dang                  ! meters

    return
  end function haversine
  ! ...
  ! ===================================================================
  ! ...
  subroutine Bearing (lon1,lat1,lon2,lat2,azimuth,rev_azimuth)
    !==================================================================
    ! ... Subroutine Bearing
    ! ... Calculation of the azimuth and reverse azimuth from (lon1,lat1) to (lon2,lat2)
    ! ... The input lngitudes and latitudes are given in radians
    ! ... The outputs are given in degrees.
    !==================================================================

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
  end subroutine Bearing
  ! ...
  ! ===================================================================
  ! ...
  real(dp) pure function compass_to_polar(compass)
    !==================================================================
    ! ... Function Compass_to_Polar
    ! ... Transform from compass angle (clockwise starting at North) to 
    ! ... polar (counter clockwise starting at east).
    !==================================================================
    
    real(dp), intent(in)                  :: compass

    compass_to_polar = mod(450.0 - compass, 360.0)

  end function compass_to_polar
  ! ...
  ! ===================================================================
  ! ...
  real(dp) pure function polar_to_compass(polar)
    !==================================================================
    ! ... Transform from polar angle (counter clockwise starting at east) to 
    ! ... compass (clockwise starting at north).
    !==================================================================

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
  subroutine vincenty_direct(lat1_deg, lon1_deg, bearing_deg, dist, lat2_deg, lon2_deg)

    real(dp), intent(in) :: lat1_deg, lon1_deg, bearing_deg, dist
    real(dp), intent(out) :: lat2_deg, lon2_deg

    real(dp) :: a, f, b
    real(dp) :: alpha1, sinAlpha1, cosAlpha1
    real(dp) :: tanU1, cosU1, sinU1, sigma1, sinSigma, cosSigma, sigma
    real(dp) :: sinAlpha, cosSqAlpha, uSq, A_coeff, B_coeff, deltaSigma
    real(dp) :: sigmaP, cos2SigmaM
    real(dp) :: lat1, lon1, lat2, lon2

    a = WGS84%Earth_SemiMajor
    f = 1.0d0 / WGS84%Earth_InvFlatten
    b = a * (1.0d0 - f)

    lat1 = lat1_deg * deg2rad
    lon1 = lon1_deg * deg2rad
    alpha1 = bearing_deg * deg2rad

    sinAlpha1 = sin(alpha1)
    cosAlpha1 = cos(alpha1)

    tanU1 = (1.0d0 - f) * tan(lat1)
    cosU1 = 1.0d0 / sqrt(1.0d0 + tanU1**2)
    sinU1 = tanU1 * cosU1

    sigma1 = atan2(tanU1, cosAlpha1)
    sinAlpha = cosU1 * sinAlpha1
    cosSqAlpha = 1.0d0 - sinAlpha**2
    uSq = cosSqAlpha * (a**2 - b**2) / (b**2)

    A_coeff = 1.0d0 + uSq / 16384.0d0 * (4096.0d0 + uSq * (-768.0d0 + uSq * (320.0d0 - 175.0d0 * uSq)))
    B_coeff = uSq / 1024.0d0 * (256.0d0 + uSq * (-128.0d0 + uSq * (74.0d0 - 47.0d0 * uSq)))

    sigma = dist / (b * A_coeff)
    sigmaP = 2.0d0 * 3.141592653589793d0
    do while (abs(sigma - sigmaP) > 1.0e-12)
      cos2SigmaM = cos(2.0d0 * sigma1 + sigma)
      sinSigma = sin(sigma)
      cosSigma = cos(sigma)
      deltaSigma = B_coeff * sinSigma * (cos2SigmaM + B_coeff / 4.0d0 * (cosSigma * (-1.0d0 + 2.0d0 * cos2SigmaM**2) - &
                    B_coeff / 6.0d0 * cos2SigmaM * (-3.0d0 + 4.0d0 * sinSigma**2) * (-3.0d0 + 4.0d0 * cos2SigmaM**2)))
      sigmaP = sigma
      sigma = dist / (b * A_coeff) + deltaSigma
    end do

    lat2 = atan2(sinU1 * cosSigma + cosU1 * sinSigma * cosAlpha1, &
              (1.0d0 - f) * sqrt(sinAlpha**2 + (sinU1 * sinSigma - cosU1 * cosSigma * cosAlpha1)**2))

    lon2 = lon1 + atan2(sinSigma * sinAlpha1, cosU1 * cosSigma - sinU1 * sinSigma * cosAlpha1)

    lat2_deg = lat2 * rad2deg
    lon2_deg = modulo(lon2 * rad2deg + 540.0d0, 360.0d0) - 180.0d0

  end subroutine vincenty_direct
  ! ...
  ! ===================================================================
  ! ...
  subroutine vincenty_inverse(lat1_deg, lon1_deg, lat2_deg, lon2_deg, dist, azimuth1, azimuth2)

    real(dp), intent(in) :: lat1_deg, lon1_deg, lat2_deg, lon2_deg
    real(dp), intent(out) :: dist, azimuth1, azimuth2

    real(dp) :: a, f, b, L, U1, U2, sinU1, cosU1, sinU2, cosU2
    real(dp) :: lambda, lambdaP, sinLambda, cosLambda, sinSigma, cosSigma, sigma
    real(dp) :: sinAlpha, cosSqAlpha, cos2SigmaM, C, uSq, A_coeff, B_coeff, deltaSigma
    integer :: iterLimit
    real(dp) :: lat1, lon1, lat2, lon2

    a = WGS84%Earth_SemiMajor
    f = 1.0d0 / WGS84%Earth_InvFlatten
    b = a * (1.0d0 - f)

    lat1 = lat1_deg * deg2rad
    lon1 = lon1_deg * deg2rad
    lat2 = lat2_deg * deg2rad
    lon2 = lon2_deg * deg2rad

    L = lon2 - lon1
    U1 = atan((1.0d0 - f) * tan(lat1))
    U2 = atan((1.0d0 - f) * tan(lat2))

    sinU1 = sin(U1)
    cosU1 = cos(U1)
    sinU2 = sin(U2)
    cosU2 = cos(U2)

    lambda = L
    iterLimit = 100
    do
      sinLambda = sin(lambda)
      cosLambda = cos(lambda)
      sinSigma = sqrt((cosU2*sinLambda)**2 + (cosU1*sinU2 - sinU1*cosU2*cosLambda)**2)
      if (sinSigma == 0.0d0) then
        dist = 0.0d0
        azimuth1 = 0.0d0
        azimuth2 = 0.0d0
        return
      end if
      cosSigma = sinU1*sinU2 + cosU1*cosU2*cosLambda
      sigma = atan2(sinSigma, cosSigma)
      sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma
      cosSqAlpha = 1.0d0 - sinAlpha**2
      if (cosSqAlpha == 0.0d0) then
        cos2SigmaM = 0.0d0
      else
        cos2SigmaM = cosSigma - 2.0d0*sinU1*sinU2/cosSqAlpha
      end if
      C = f/16.0d0*cosSqAlpha*(4.0d0 + f*(4.0d0 - 3.0d0*cosSqAlpha))
      lambdaP = lambda
      lambda = L + (1.0d0 - C)*f*sinAlpha*(sigma + C*sinSigma*(cos2SigmaM + C*cosSigma*(-1.0d0 + 2.0d0*cos2SigmaM**2)))
      iterLimit = iterLimit - 1
      if (abs(lambda - lambdaP) < 1.0e-12 .or. iterLimit == 0) exit
    end do

    uSq = cosSqAlpha*(a**2 - b**2)/(b**2)
    A_coeff = 1.0d0 + uSq/16384.0d0*(4096.0d0 + uSq*(-768.0d0 + uSq*(320.0d0 - 175.0d0*uSq)))
    B_coeff = uSq/1024.0d0*(256.0d0 + uSq*(-128.0d0 + uSq*(74.0d0 - 47.0d0*uSq)))
    deltaSigma = B_coeff*sinSigma*(cos2SigmaM + B_coeff/4.0d0*(cosSigma*(-1.0d0 + 2.0d0*cos2SigmaM**2) - &
                 B_coeff/6.0d0*cos2SigmaM*(-3.0d0 + 4.0d0*sinSigma**2)*(-3.0d0 + 4.0d0*cos2SigmaM**2)))

    dist = b*A_coeff*(sigma - deltaSigma)
    azimuth1 = modulo(atan2(cosU2*sinLambda, cosU1*sinU2 - sinU1*cosU2*cosLambda)*rad2deg + 360.0d0, 360.0d0)
    azimuth2 = modulo(atan2(cosU1*sinLambda, -sinU1*cosU2 + cosU1*sinU2*cosLambda)*rad2deg + 360.0d0, 360.0d0)

  end subroutine vincenty_inverse
  ! ...
  ! ===================================================================
  ! ...
  subroutine spdir2uv_r(spd,ang,deg,u,v) 
    ! ... Calculates (u,v) coordinates from velocity and angles measured
    ! ... from north (compass/bearing), increasing clockwise.
    ! ...      ang = 0°   → wind toward north (u=0, v=spd)
    ! ...      ang = 90°  → wind toward east  (u=spd, v=0)
    ! ...      ang = 180° → wind toward south (u=0, v=-spd)
    ! ...      ang = 270° → wind toward west  (u=-spd, v=0)
    ! ...
    real(dp), intent(in)                         :: spd
    real(dp), intent(in)                         :: ang  ! Compass angle
    logical, intent(in), optional                :: deg
    real(dp), intent(out)                        :: u,v

    ! ... Local variables
    ! ...
    real(dp) alpha

    alpha = ang
    if (present(deg)) then
      if (deg) alpha = alpha * deg2rad
    endif
     
    ! ... Calculate u (E-W) and v (N-S) components
    ! ...
    u = spd * sin(alpha)
    v = spd * cos(alpha)

  end subroutine spdir2uv_r
  ! ...
  ! ===================================================================
  ! ...
  subroutine spdir2uv_v(spd,ang,deg,u,v)
    ! ... Calculates (u,v) coordinates from velocity and angles measured
    ! ... from north (compass/bearing), increasing clockwise.
    ! ...      ang = 0°   → wind toward north (u=0, v=spd)
    ! ...      ang = 90°  → wind toward east  (u=spd, v=0)
    ! ...      ang = 180° → wind toward south (u=0, v=-spd)
    ! ...      ang = 270° → wind toward west  (u=-spd, v=0)
    ! ...
    real(dp), dimension(:), intent(in)           :: spd
    real(dp), dimension(:), intent(in)           :: ang  ! Compass angle
    logical, intent(in), optional                :: deg
    real(dp), dimension(size(spd)), intent(out)  :: u,v

    ! ... Local variables
    ! ...
    real(dp), dimension(size(spd))               ::  alpha

    alpha(:) = ang(:)
    if (present(deg)) then
      if (deg) alpha(:) = alpha(:) * deg2rad
    endif
     
    ! ... Calculate u (E-W) and v (N-S) components
    ! ...
    u(:) = spd(:) * sin(alpha(:))
    v(:) = spd(:) * cos(alpha(:))

  end subroutine spdir2uv_v
  ! ...
  ! ===================================================================
  ! ...
  subroutine uv2spdir_r(u,v,spd,ang) 
    ! ... Computes direction (ang) and speed (spd) from U/V wind components
    ! ... Geographic convention: 0° = North, angles increase clockwise.
    ! ...
    real(dp), intent(in)                         :: u
    real(dp), intent(in)                         :: v 
    real(dp), intent(out)                        :: spd
    real(dp), intent(out)                        :: ang  ! Compass angle

    ! ... Speed: magnitude of the velocity vector
    ! ...
    spd = sqrt(u*u + v*v)

    ! ... Compute angle in radians from atan2
    ! ... atan2(y, x) = atan2(v, u) gives math convention from +x axis (east, CCW)
    ! ... To get geographic convention (0°=north, clockwise), use:
    ! ... angle_geo = mod(atan2(u, v) * rad2deg, 360.0)
    ! ...
    ang = atan2(u,v) * rad2deg
    if (ang .lt. 0.0) ang = ang + 360.0

  end subroutine uv2spdir_r
  ! ...
  ! ===================================================================
  ! ...
  subroutine uv2spdir_v(u,v,spd,ang)
    ! ... Computes direction (ang) and speed (spd) from U/V wind components
    ! ... Geographic convention: 0° = North, angles increase clockwise.
    ! ...
    real(dp), dimension(:), intent(in)           :: u
    real(dp), dimension(:), intent(in)           :: v 
    real(dp), dimension(size(u)), intent(out)    :: spd
    real(dp), dimension(size(u)), intent(out)    :: ang  ! Compass angle

    ! ... Local variables
    ! ...
    integer i

    ! ... Speed: magnitude of the velocity vector
    ! ...
    spd(:) = sqrt(u(:)**2 + v(:)**2)

    ! ... Compute angle in radians from atan2
    ! ... atan2(y, x) = atan2(v, u) gives math convention from +x axis (east, CCW)
    ! ... To get geographic convention (0°=north, clockwise), use:
    ! ... angle_geo = mod(atan2(u, v) * rad2deg, 360.0)
    ! ...
    do i = 1, size(u)
       ang(i) = atan2(u(i), v(i)) * rad2deg
       if (ang(i) < 0.0) ang(i) = ang(i) + 360.0
    end do

  end subroutine uv2spdir_v
  ! ...
  ! ===================================================================
  ! ...
  subroutine dms2dec(gg,mm,ss,deg) 
    ! ... Converts degrees, minutes and seconds (integers) to decimal degrees (double precission).

    integer, intent(in)                            :: gg,mm,ss
    real(dp), intent(out)                          :: deg

    deg = abs(gg) + (mm + ss/60.0D0)/60.0D0
    if (gg.lt.0.0D0) deg = -deg
   
  end subroutine dms2dec 
  ! ...
  ! =====================================================================
  ! ...
  subroutine dec2dms(deg,gg,mm,ss) 
    ! ... Converts decimal degrees (double precission) to degrees, minutes and seconds (integers)

    real(dp), intent(in)                           :: deg
    integer, intent(out)                           :: gg,mm,ss

    real(dp) sec
    integer res

    sec = anint(abs(deg) * 3600.0D0)
    ss  = mod(nint(sec),60)
    res = nint((sec - ss) / 60.0D0) 
    mm  = mod(res,60)
    gg  = (res - mm)/60
    if (deg.lt.0) gg = -gg
    
  end subroutine dec2dms 
  ! ...
  ! =====================================================================
  ! ...
end module module_geodesy
