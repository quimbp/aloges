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
! - haversine                                                              !
! - gc_bearing_deg                                                         !
! - compass_to_polar                                                       !
! - polar_to_compass                                                       !
! - gc_destination_point                                                   !
! - gc_inverse                                                             !
! - vincent_direct                                                         !
! - vincent_inverse                                                        !
! - spdir2uv                                                               !
! - uv2spdir                                                               !
! - gms2dec                                                                !
! - dec2gms                                                                !
! - ll2m                                                                   !
! -------------------------------------------------------------------------!

module module_geodesy

  use, intrinsic :: IEEE_ARITHMETIC, only : IEEE_VALUE, IEEE_QUIET_NAN
  use module_types,     only: dp
  use module_constants, only: pi, two_pi, deg2rad, rad2deg

  implicit none (type, external)
  private
  public :: haversine, gc_bearing_deg, compass_to_polar, polar_to_compass
  public :: gc_destination_point, gc_inverse
  public :: vincenty_direct, vincenty_inverse
  public :: spdir2uv, uv2spdir
  public :: dms2dec, dm2dec, dec2dms, dec2dm
  public :: ll2m

  !> @brief Canonical WGS84 constant set (single source of truth).
  !! 
  !! Stores spherical mean radius plus ellipsoidal parameters. Dependent
  !! values (semi-minor axis, eccentricity squared) are derived at runtime.
  !! Units: meters (length), dimensionless (ratios).
  type :: type_WGS84_constants
    real(dp) :: Earth_Radius    = 6371008.8_dp      !< Mean spherical Earth radius [m]
    real(dp) :: Earth_SemiMajor = 6378137.0_dp      !< WGS84 semi-major axis a [m]
    real(dp) :: Earth_InvFlatten= 298.257223563_dp  !< Inverse flattening 1/f
    real(dp) :: Earth_SemiMinor = 0.0_dp            !< WGS84 semi-minor axis b [m] (derived)
    real(dp) :: Earth_Ecc2      = 0.0_dp            !< First eccentricity squared e^2 (derived)
  end type type_WGS84_constants

  type(type_WGS84_constants) :: WGS84

  interface spdir2uv
    module procedure spdir2uv_r, spdir2uv_v
  end interface spdir2uv

  interface uv2spdir
    module procedure uv2spdir_r, uv2spdir_v
  end interface uv2spdir

contains

  !> @brief Initialize WGS84 derived constants once.
  !!
  !! Derives semi-minor axis b and eccentricity squared e^2 from a and 1/f.
  !! Called lazily by public routines; no need to call directly.
  !! Units: meters (a,b), dimensionless (e^2).
  subroutine init_wgs84_once()
    real(dp) :: a, f, b
    if (WGS84%Earth_SemiMinor == 0.0_dp .and. WGS84%Earth_Ecc2 == 0.0_dp) then
      a = WGS84%Earth_SemiMajor
      f = 1.0_dp / WGS84%Earth_InvFlatten
      b = a * (1.0_dp - f)
      WGS84%Earth_SemiMinor = b
      WGS84%Earth_Ecc2      = 1.0_dp - (b*b)/(a*a)
    end if
  end subroutine init_wgs84_once

  !> @brief Great-circle distance using the haversine formula (sphere).
  !!
  !! Inputs are in radians; output distance is in meters using mean Earth radius.
  !!
  !! @param[in] lon1 Longitude of point 1 [rad]
  !! @param[in] lat1 Latitude  of point 1 [rad]
  !! @param[in] lon2 Longitude of point 2 [rad]
  !! @param[in] lat2 Latitude  of point 2 [rad]
  !!
  !! @return Real(dp) distance along great circle [m].
  !!
  !! @note Robust for all separations; numerically stable.
  real(dp) function haversine(lon1,lat1,lon2,lat2)
    real(dp), intent(in) :: lon1, lat1, lon2, lat2
    real(dp) :: sindx, sindy, dang, ccl1, ccl2
    call init_wgs84_once()
    sindx = sin(0.5_dp*(lon2-lon1))
    sindy = sin(0.5_dp*(lat2-lat1))
    ccl1  = cos(lat1); ccl2 = cos(lat2)
    dang  = 2.0_dp*asin( sqrt(sindy*sindy + ccl2*ccl1*sindx*sindx) )
    haversine = WGS84%Earth_Radius * dang
  end function haversine

  !> @brief Forward and reverse azimuth between two points on sphere.
  !!
  !! Inputs are in radians; azimuth outputs in degrees [0, 360).
  !!
  !! @param[in]  lon1 Longitude of start point [rad]
  !! @param[in]  lat1 Latitude  of start point [rad]
  !! @param[in]  lon2 Longitude of end   point [rad]
  !! @param[in]  lat2 Latitude  of end   point [rad]
  !! @param[out] azimuth     Initial bearing from point 1 to 2 [deg]
  !! @param[out] rev_azimuth Reverse bearing from point 2 to 1 [deg]
  !!
  !! @note Spherical geometry; for ellipsoidal bearings use Vincenty inverse.
  subroutine gc_bearing_deg(lon1,lat1,lon2,lat2,azimuth,rev_azimuth)
    real(dp), intent(in)  :: lon1, lat1, lon2, lat2
    real(dp), intent(out) :: azimuth, rev_azimuth
    real(dp) :: dlon, cdlon, sdlon, x, y, xr, yr
    dlon  = lon2 - lon1
    cdlon = cos(dlon); sdlon = sin(dlon)
    y = sdlon * cos(lat2)
    x = cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cdlon
    azimuth = modulo(atan2(y, x)*rad2deg + 360.0_dp, 360.0_dp)
    yr = -sdlon * cos(lat1)
    xr =  cos(lat2)*sin(lat1) - sin(lat2)*cos(lat1)*cdlon
    rev_azimuth = modulo(atan2(yr, xr)*rad2deg + 360.0_dp, 360.0_dp)
  end subroutine gc_bearing_deg

  !> @brief Convert compass angle (0°=North, CW) to math-polar (0°=East, CCW).
  !!
  !! @param[in] compass Compass/bearing angle [deg]
  !! @return Real(dp) polar angle [deg].
  pure real(dp) function compass_to_polar(compass)
    real(dp), intent(in) :: compass
    compass_to_polar = modulo(450.0_dp - compass, 360.0_dp)
  end function compass_to_polar

  !> @brief Convert math-polar angle (0°=East, CCW) to compass (0°=North, CW).
  !!
  !! @param[in] polar Polar angle [deg]
  !! @return Real(dp) compass angle [deg].
  pure real(dp) function polar_to_compass(polar)
    real(dp), intent(in) :: polar
    polar_to_compass = modulo(450.0_dp - polar, 360.0_dp)
  end function polar_to_compass

  !> @brief Destination point given start, bearing, and distance (sphere).
  !!
  !! Great-circle forward problem with mean spherical radius.
  !!
  !! @param[in]  lat0_deg Start latitude  [deg]
  !! @param[in]  lon0_deg Start longitude [deg]
  !! @param[in]  bearing_deg Initial bearing from start [deg]
  !! @param[in]  dist Distance along great circle [m]
  !! @param[out] lat1_deg Destination latitude  [deg]
  !! @param[out] lon1_deg Destination longitude [deg]
  !!
  !! @note Longitudes normalized to [-180, 180).
  subroutine gc_destination_point(lat0_deg, lon0_deg, bearing_deg, dist, lat1_deg, lon1_deg)
    real(dp), intent(in)  :: lat0_deg, lon0_deg, bearing_deg, dist
    real(dp), intent(out) :: lat1_deg, lon1_deg
    real(dp) :: lat0, lon0, bearing, lat1, lon1, R
    call init_wgs84_once()
    R = WGS84%Earth_Radius
    lat0    = lat0_deg * deg2rad
    lon0    = lon0_deg * deg2rad
    bearing = bearing_deg * deg2rad
    lat1 = asin( sin(lat0)*cos(dist/R) + cos(lat0)*sin(dist/R)*cos(bearing) )
    lon1 = lon0 + atan2( sin(bearing)*sin(dist/R)*cos(lat0), cos(dist/R) - sin(lat0)*sin(lat1) )
    lon1 = modulo(lon1 + pi, two_pi) - pi
    lat1_deg = lat1 * rad2deg
    lon1_deg = lon1 * rad2deg
  end subroutine gc_destination_point

  !> @brief Initial bearing and spherical great-circle distance between two points.
  !!
  !! Uses spherical law of cosines with clamped acos for robustness.
  !!
  !! @param[in]  lat1_deg Latitude  of point 1 [deg]
  !! @param[in]  lon1_deg Longitude of point 1 [deg]
  !! @param[in]  lat2_deg Latitude  of point 2 [deg]
  !! @param[in]  lon2_deg Longitude of point 2 [deg]
  !! @param[out] dist Great-circle distance [m]
  !! @param[out] bearing_deg Initial bearing from point 1 to 2 [deg]
  subroutine gc_inverse(lat1_deg, lon1_deg, lat2_deg, lon2_deg, dist, bearing_deg)
    real(dp), intent(in)  :: lat1_deg, lon1_deg, lat2_deg, lon2_deg
    real(dp), intent(out) :: dist, bearing_deg
    real(dp) :: lat1, lon1, lat2, lon2, dlon, cosc, c, y, x, R
    call init_wgs84_once()
    R   = WGS84%Earth_Radius
    lat1 = lat1_deg * deg2rad; lon1 = lon1_deg * deg2rad
    lat2 = lat2_deg * deg2rad; lon2 = lon2_deg * deg2rad
    dlon = lon2 - lon1
    cosc = sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(dlon)
    cosc = max(-1.0_dp, min(1.0_dp, cosc))
    c    = acos(cosc)
    dist = R * c
    y = sin(dlon) * cos(lat2)
    x = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dlon)
    bearing_deg = modulo(atan2(y, x) * rad2deg + 360.0_dp, 360.0_dp)
  end subroutine gc_inverse

  !> @brief Vincenty forward problem on WGS84 ellipsoid.
  !!
  !! Computes destination point given start, initial bearing, and distance.
  !!
  !! @param[in]  lat1_deg Start latitude  [deg]
  !! @param[in]  lon1_deg Start longitude [deg]
  !! @param[in]  bearing_deg Initial azimuth from start [deg]
  !! @param[in]  dist Geodesic distance on ellipsoid [m]
  !! @param[out] lat2_deg Destination latitude  [deg]
  !! @param[out] lon2_deg Destination longitude [deg]
  !!
  !! @note Iterative; typically converges rapidly except near poles with
  !! near-antipodal cases (rare for forward problem).
  subroutine vincenty_direct(lat1_deg, lon1_deg, bearing_deg, dist, lat2_deg, lon2_deg)
    real(dp), intent(in)  :: lat1_deg, lon1_deg, bearing_deg, dist
    real(dp), intent(out) :: lat2_deg, lon2_deg
    real(dp) :: a, f, b, alpha1, sinAlpha1, cosAlpha1
    real(dp) :: tanU1, cosU1, sinU1, sigma1, sinSigma, cosSigma, sigma
    real(dp) :: sinAlpha, cosSqAlpha, uSq, A_coeff, B_coeff, deltaSigma
    real(dp) :: sigmaP, cos2SigmaM
    real(dp) :: lat1, lon1, lat2, lon2
    call init_wgs84_once()
    a = WGS84%Earth_SemiMajor
    f = 1.0_dp / WGS84%Earth_InvFlatten
    b = WGS84%Earth_SemiMinor
    lat1   = lat1_deg   * deg2rad
    lon1   = lon1_deg   * deg2rad
    alpha1 = bearing_deg* deg2rad
    sinAlpha1 = sin(alpha1); cosAlpha1 = cos(alpha1)
    tanU1 = (1.0_dp - f) * tan(lat1)
    cosU1 = 1.0_dp / sqrt(1.0_dp + tanU1*tanU1)
    sinU1 = tanU1 * cosU1
    sigma1   = atan2(tanU1, cosAlpha1)
    sinAlpha = cosU1 * sinAlpha1
    cosSqAlpha = 1.0_dp - sinAlpha*sinAlpha
    uSq = cosSqAlpha * (a*a - b*b) / (b*b)
    A_coeff = 1.0_dp + uSq/16384.0_dp*(4096.0_dp + uSq*(-768.0_dp + uSq*(320.0_dp - 175.0_dp*uSq)))
    B_coeff = uSq/1024.0_dp*(256.0_dp  + uSq*(-128.0_dp + uSq*(74.0_dp  -  47.0_dp*uSq)))
    sigma = dist / (b * A_coeff)
    sigmaP = two_pi
    do while (abs(sigma - sigmaP) > 1.0e-12_dp)
      cos2SigmaM = cos(2.0_dp*sigma1 + sigma)
      sinSigma   = sin(sigma)
      cosSigma   = cos(sigma)
      deltaSigma = B_coeff*sinSigma * ( cos2SigmaM + (B_coeff/4.0_dp)*( cosSigma*(-1.0_dp + 2.0_dp*cos2SigmaM*cos2SigmaM) &
                    - (B_coeff/6.0_dp)*cos2SigmaM*(-3.0_dp + 4.0_dp*sinSigma*sinSigma)*(-3.0_dp + 4.0_dp*cos2SigmaM*cos2SigmaM) ) )
      sigmaP = sigma
      sigma  = dist/(b*A_coeff) + deltaSigma
    end do
    lat2 = atan2( sinU1*cosSigma + cosU1*sinSigma*cosAlpha1, (1.0_dp - f)*sqrt( sinAlpha*sinAlpha + (sinU1*sinSigma - cosU1*cosSigma*cosAlpha1)**2 ) )
    lon2 = lon1 + atan2( sinSigma*sinAlpha1, cosU1*cosSigma - sinU1*sinSigma*cosAlpha1 )
    lat2_deg = lat2 * rad2deg
    lon2_deg = modulo(lon2 * rad2deg + 540.0_dp, 360.0_dp) - 180.0_dp
  end subroutine vincenty_direct

  !> @brief Vincenty inverse problem on WGS84 ellipsoid.
  !!
  !! Computes geodesic distance and forward/back azimuths between two points.
  !!
  !! @param[in]  lat1_deg Latitude  of point 1 [deg]
  !! @param[in]  lon1_deg Longitude of point 1 [deg]
  !! @param[in]  lat2_deg Latitude  of point 2 [deg]
  !! @param[in]  lon2_deg Longitude of point 2 [deg]
  !! @param[out] dist     Geodesic distance on ellipsoid [m]
  !! @param[out] azimuth1 Initial azimuth at point 1 [deg]
  !! @param[out] azimuth2 Reverse azimuth at point 2 [deg]
  !!
  !! @note Iterative; may fail to converge for nearly antipodal points in
  !! some edge cases. A maximum of 100 iterations is used.
  subroutine vincenty_inverse(lat1_deg, lon1_deg, lat2_deg, lon2_deg, dist, azimuth1, azimuth2)
    real(dp), intent(in)  :: lat1_deg, lon1_deg, lat2_deg, lon2_deg
    real(dp), intent(out) :: dist, azimuth1, azimuth2
    real(dp) :: a, f, b, L, U1, U2, sinU1, cosU1, sinU2, cosU2
    real(dp) :: lambda, lambdaP, sinLambda, cosLambda, sinSigma, cosSigma, sigma
    real(dp) :: sinAlpha, cosSqAlpha, cos2SigmaM, C, uSq, A_coeff, B_coeff, deltaSigma
    integer  :: iterLimit
    real(dp) :: lat1, lon1, lat2, lon2
    call init_wgs84_once()
    a = WGS84%Earth_SemiMajor
    f = 1.0_dp / WGS84%Earth_InvFlatten
    b = WGS84%Earth_SemiMinor
    lat1 = lat1_deg * deg2rad; lon1 = lon1_deg * deg2rad
    lat2 = lat2_deg * deg2rad; lon2 = lon2_deg * deg2rad
    L  = lon2 - lon1
    U1 = atan( (1.0_dp - f) * tan(lat1) )
    U2 = atan( (1.0_dp - f) * tan(lat2) )
    sinU1 = sin(U1); cosU1 = cos(U1)
    sinU2 = sin(U2); cosU2 = cos(U2)
    lambda    = L
    iterLimit = 100
    do
      sinLambda = sin(lambda); cosLambda = cos(lambda)
      sinSigma = sqrt( (cosU2*sinLambda)**2 + (cosU1*sinU2 - sinU1*cosU2*cosLambda)**2 )
      if (sinSigma == 0.0_dp) then
        dist = 0.0_dp; azimuth1 = 0.0_dp; azimuth2 = 0.0_dp
        return
      end if
      cosSigma = sinU1*sinU2 + cosU1*cosU2*cosLambda
      sigma    = atan2(sinSigma, cosSigma)
      sinAlpha = (cosU1*cosU2*sinLambda) / sinSigma
      cosSqAlpha = 1.0_dp - sinAlpha*sinAlpha
      if (cosSqAlpha == 0.0_dp) then
        cos2SigmaM = 0.0_dp
      else
        cos2SigmaM = cosSigma - (2.0_dp*sinU1*sinU2)/cosSqAlpha
      end if
      C = (f/16.0_dp) * cosSqAlpha * (4.0_dp + f*(4.0_dp - 3.0_dp*cosSqAlpha))
      lambdaP = lambda
      lambda  = L + (1.0_dp - C)*f*sinAlpha * ( sigma + C*sinSigma*( cos2SigmaM + C*cosSigma*(-1.0_dp + 2.0_dp*cos2SigmaM*cos2SigmaM) ) )
      iterLimit = iterLimit - 1
      if (abs(lambda - lambdaP) < 1.0e-12_dp .or. iterLimit == 0) exit
    end do
    uSq = cosSqAlpha*(a*a - b*b)/(b*b)
    A_coeff = 1.0_dp + uSq/16384.0_dp*(4096.0_dp + uSq*(-768.0_dp + uSq*(320.0_dp - 175.0_dp*uSq)))
    B_coeff = uSq/1024.0_dp *(256.0_dp + uSq*(-128.0_dp + uSq*(74.0_dp  - 47.0_dp*uSq)))
    deltaSigma = B_coeff*sinSigma*( cos2SigmaM + (B_coeff/4.0_dp)*( cosSigma*(-1.0_dp + 2.0_dp*cos2SigmaM*cos2SigmaM) &
                 - (B_coeff/6.0_dp)*cos2SigmaM*(-3.0_dp + 4.0_dp*sinSigma*sinSigma)*(-3.0_dp + 4.0_dp*cos2SigmaM*cos2SigmaM) ) )
    dist     = b*A_coeff*(sigma - deltaSigma)
    azimuth1 = modulo(atan2( cosU2*sinLambda,  cosU1*sinU2 - sinU1*cosU2*cosLambda )*rad2deg + 360.0_dp, 360.0_dp)
    azimuth2 = modulo(atan2( cosU1*sinLambda, -sinU1*cosU2 + cosU1*sinU2*cosLambda )*rad2deg + 360.0_dp, 360.0_dp)
  end subroutine vincenty_inverse

  !> @brief Convert speed and compass direction to (u,v) components (scalar).
  !!
  !! Compass convention: 0°=North, angles increase clockwise.
  !!
  !! @param[in]  spd Speed magnitude [units of choice]
  !! @param[in]  ang Compass/bearing angle [deg if deg=.true., else rad]
  !! @param[in]  deg Logical flag: input angle in degrees if true
  !! @param[out] u Eastward component
  !! @param[out] v Northward component
  subroutine spdir2uv_r(spd, ang, deg, u, v)
    real(dp), intent(in)            :: spd, ang
    logical, intent(in), optional   :: deg
    real(dp), intent(out)           :: u, v
    real(dp) :: alpha
    alpha = ang
    if (present(deg)) then
      if (deg) alpha = alpha * deg2rad
    end if
    u = spd * sin(alpha)
    v = spd * cos(alpha)
  end subroutine spdir2uv_r

  !> @brief Convert speed and compass direction to (u,v) components (vectorized).
  !!
  !! See spdir2uv_r for conventions.
  !!
  !! @param[in]  spd(:) Speeds
  !! @param[in]  ang(:) Angles [deg if deg=.true., else rad]
  !! @param[in]  deg Logical flag: input angle in degrees if true
  !! @param[out] u(:) Eastward components
  !! @param[out] v(:) Northward components
  subroutine spdir2uv_v(spd, ang, deg, u, v)
    real(dp), intent(in)            :: spd(:), ang(:)
    logical, intent(in), optional   :: deg
    real(dp), intent(out)           :: u(size(spd)), v(size(spd))
    real(dp) :: alpha(size(spd))
    alpha = ang
    if (present(deg)) then
      if (deg) alpha = alpha * deg2rad
    end if
    u = spd * sin(alpha)
    v = spd * cos(alpha)
  end subroutine spdir2uv_v

  !> @brief Convert (u,v) velocity components to speed and compass direction (scalar).
  !!
  !! Compass convention: 0°=North, clockwise positive.
  !!
  !! @param[in]  u Eastward component
  !! @param[in]  v Northward component
  !! @param[out] spd Speed magnitude
  !! @param[out] ang Compass/bearing angle [deg]
  subroutine uv2spdir_r(u, v, spd, ang)
    real(dp), intent(in)  :: u, v
    real(dp), intent(out) :: spd, ang
    spd = sqrt(u*u + v*v)
    ang = atan2(u, v) * rad2deg
    if (ang < 0.0_dp) ang = ang + 360.0_dp
  end subroutine uv2spdir_r

  !> @brief Convert (u,v) arrays to speed and compass direction (vectorized).
  !!
  !! @param[in]  u(:) Eastward components
  !! @param[in]  v(:) Northward components
  !! @param[out] spd(:) Speed magnitudes
  !! @param[out] ang(:) Compass/bearing angles [deg]
  subroutine uv2spdir_v(u, v, spd, ang)
    real(dp), intent(in)  :: u(:), v(:)
    real(dp), intent(out) :: spd(size(u)), ang(size(u))
    integer :: i
    spd = sqrt(u*u + v*v)
    do i = 1, size(u)
      ang(i) = atan2(u(i), v(i)) * rad2deg
      if (ang(i) < 0.0_dp) ang(i) = ang(i) + 360.0_dp
    end do
  end subroutine uv2spdir_v

  !> @brief Convert degrees-minutes-seconds to decimal degrees.
  !!
  !! @param[in]  gg Degrees (integer, sign carries hemisphere)
  !! @param[in]  mm Minutes  (0..59)
  !! @param[in]  ss Seconds  (0..59)
  !! @param[out] deg Decimal degrees [deg]
  subroutine dms2dec(gg, mm, ss, deg)
    integer, intent(in) :: gg, mm, ss
    real(dp), intent(out) :: deg
    deg = abs(gg) + (mm + ss/60.0_dp)/60.0_dp
    if (gg < 0) deg = -deg
  end subroutine dms2dec

  !> @brief Convert degrees-minutes to decimal degrees.
  !!
  !! @param[in]  gg Degrees (integer, sign carries hemisphere)
  !! @param[in]  mm Minutes (real, includes seconds as decimal part)
  !! @param[out] deg Decimal degrees [deg]
  subroutine dm2dec(gg, mm, deg)
     integer,  intent(in)           :: gg  
     real(dp), intent(in)           :: mm      ! Fractional minutes  
     real(dp), intent(out)          :: deg  
     real(dp) :: sign_, mm_abs, ss_use  
     sign_  = merge(-1.0_dp, 1.0_dp, gg < 0)  
     mm_abs = abs(mm)  
     ! minutes arefractional already  
     deg = abs(gg) + mm_abs/60.0_dp  
     deg = sign_ * deg  
  end subroutine dm2dec

  !> @brief Convert decimal degrees to degrees-minutes-seconds.
  !!
  !! @param[in]  deg Decimal degrees [deg]
  !! @param[out] gg Degrees (integer, sign carries hemisphere)
  !! @param[out] mm Minutes  (0..59)
  !! @param[out] ss Seconds  (0..59)
  subroutine dec2dms(deg, gg, mm, ss)
    real(dp), intent(in) :: deg
    integer, intent(out) :: gg, mm, ss
    real(dp) :: sec
    integer :: res
    sec = anint(abs(deg) * 3600.0_dp)
    ss  = mod(nint(sec), 60)
    res = nint((sec - ss) / 60.0_dp)
    mm  = mod(res, 60)
    gg  = (res - mm)/60
    if (deg < 0.0_dp) gg = -gg
  end subroutine dec2dms

  !> @brief Convert decimal degrees to degrees-minutes.
  !!
  !! @param[in]  deg Decimal degrees [deg]
  !! @param[out] gg Degrees (integer, sign carries hemisphere)
  !! @param[out] mm Minutes  (real, seconds in the decimal part)
  subroutine dec2dm(deg, gg, mm)
    real(dp), intent(in) :: deg
    integer, intent(out) :: gg
    real(dp), intent(out) :: mm
    real(dp) :: adeg, frac_deg

    adeg = abs(deg)  
    gg   = int(floor(adeg))       ! degrees magnitude  
    frac_deg = adeg - real(gg, dp)  
    mm = frac_deg * 60.0_dp  ! return fractional minutes; avoid integer seconds  
    if (deg < 0.0_dp) gg = -gg  
  end subroutine dec2dm

  !> @brief Project two lon/lat grids to a common local ENU in meters.
  !!
  !! Creates planar east/north coordinates (meters) for both input grids
  !! around a common reference (lon0, lat0). If not provided, the reference
  !! is the midpoint of the combined domain; if provided as -999, it is
  !! overwritten with the computed midpoint.
  !!
  !! @param[in]  xi_deg(:) Longitudes of source grid [deg]
  !! @param[in]  yi_deg(:) Latitudes  of source grid [deg]
  !! @param[in]  xc_deg(:) Longitudes of target grid [deg]
  !! @param[in]  yc_deg(:) Latitudes  of target grid [deg]
  !! @param[out] xim(:) Eastings of source grid [m]
  !! @param[out] yim(:) Northings of source grid [m]
  !! @param[out] xcm(:) Eastings of target grid [m]
  !! @param[out] ycm(:) Northings of target grid [m]
  !! @param[inout] lon0_deg Optional reference longitude [deg]
  !! @param[inout] lat0_deg Optional reference latitude  [deg]
  !!
  !! @note Valid for regional domains (rule of thumb: < ~1000 km).
  !! Uses WGS84 ellipsoidal radii of curvature at reference latitude.
  subroutine ll2m(xi_deg, yi_deg, xc_deg, yc_deg, xim, yim, xcm, ycm, lon0_deg, lat0_deg)
    real(dp), intent(in)  :: xi_deg(:), yi_deg(:), xc_deg(:), yc_deg(:)
    real(dp), intent(out) :: xim(size(xi_deg)), yim(size(yi_deg)), xcm(size(xc_deg)), ycm(size(yc_deg))
    real(dp), intent(inout), optional :: lon0_deg, lat0_deg
    real(dp) :: lon0, lat0, lon0_in_deg, lat0_in_deg
    real(dp) :: lonc_deg, latc_deg, phi0, sinphi0, cosphi0
    real(dp) :: N_phi0, M_phi0
    integer :: n, m

    call init_wgs84_once()

    n = size(xi_deg); m = size(xc_deg)
    if (size(yi_deg) /= n) stop 'll2m: xi/yi size mismatch'
    if (size(yc_deg) /= m) stop 'll2m: xc/yc size mismatch'

    call bbox_midpoint(xi_deg, yi_deg, xc_deg, yc_deg, lonc_deg, latc_deg)

    if (present(lon0_deg)) then
      if (lon0_deg /= -999.0_dp) then
        lon0_in_deg = lon0_deg
      else
        lon0_in_deg = lonc_deg
        lon0_deg    = lonc_deg
      end if
    else
      lon0_in_deg = lonc_deg
    end if

    if (present(lat0_deg)) then
      if (lat0_deg /= -999.0_dp) then
        lat0_in_deg = lat0_deg
      else
        lat0_in_deg = latc_deg
        lat0_deg    = latc_deg
      end if
    else
      lat0_in_deg = latc_deg
    end if

    lat0 = deg2rad * lat0_in_deg
    lon0 = deg2rad * lon0_in_deg
    phi0 = lat0
    sinphi0 = sin(phi0); cosphi0 = cos(phi0)

    N_phi0 = WGS84%Earth_SemiMajor / sqrt(1.0_dp - WGS84%Earth_Ecc2 * sinphi0*sinphi0)
    M_phi0 = WGS84%Earth_SemiMajor * (1.0_dp - WGS84%Earth_Ecc2) / ( (1.0_dp - WGS84%Earth_Ecc2*sinphi0*sinphi0)**1.5_dp )

    call project_set(xi_deg, yi_deg, lon0, lat0, N_phi0, M_phi0, cosphi0, xim, yim)
    call project_set(xc_deg, yc_deg, lon0, lat0, N_phi0, M_phi0, cosphi0, xcm, ycm)
  end subroutine ll2m

  !> @brief Internal helper: project a set of lon/lat to local ENU meters.
  !!
  !! Converts to radians, wraps longitudes around reference, and applies
  !! linearized ENU scaling using radii of curvature at reference latitude.
  !!
  !! @param[in]  lon_deg(:) Longitudes [deg]
  !! @param[in]  lat_deg(:) Latitudes  [deg]
  !! @param[in]  lon0 Reference longitude [rad]
  !! @param[in]  lat0 Reference latitude  [rad]
  !! @param[in]  N0 Prime vertical radius (N) at reference [m]
  !! @param[in]  M0 Meridional radius (M) at reference [m]
  !! @param[in]  cosphi0 cos(latitude0), for efficiency
  !! @param[out] x_m(:) Eastings [m]
  !! @param[out] y_m(:) Northings [m]
  subroutine project_set(lon_deg, lat_deg, lon0, lat0, N0, M0, cosphi0, x_m, y_m)
    real(dp), intent(in)  :: lon_deg(:), lat_deg(:)
    real(dp), intent(in)  :: lon0, lat0, N0, M0, cosphi0
    real(dp), intent(out) :: x_m(size(lon_deg)), y_m(size(lat_deg))
    real(dp), allocatable :: lam(:), phi(:), dlam(:), dphi(:)
    integer :: k, n
    n = size(lon_deg)
    allocate(lam(n), phi(n), dlam(n), dphi(n))
    do k = 1, n
      lam(k) = deg2rad * lon_deg(k)
      phi(k) = deg2rad * lat_deg(k)
    end do
    do k = 1, n
      dlam(k) = wrap_pi(lam(k) - lon0)
      dphi(k) = phi(k) - lat0
      x_m(k)  = dlam(k) * N0 * cosphi0
      y_m(k)  = dphi(k) * M0
    end do
    deallocate(lam, phi, dlam, dphi)
  end subroutine project_set

  !> @brief Compute reference lon/lat as midpoint of combined domain.
  !!
  !! Latitude midpoint is arithmetic between extrema; longitude midpoint
  !! uses circular mean to avoid dateline artifacts.
  !!
  !! @param[in]  xi_deg(:) Source longitudes [deg]
  !! @param[in]  yi_deg(:) Source latitudes  [deg]
  !! @param[in]  xc_deg(:) Target longitudes [deg]
  !! @param[in]  yc_deg(:) Target latitudes  [deg]
  !! @param[out] lon_mid_deg Midpoint longitude [deg]
  !! @param[out] lat_mid_deg Midpoint latitude  [deg]
  subroutine bbox_midpoint(xi_deg, yi_deg, xc_deg, yc_deg, lon_mid_deg, lat_mid_deg)
    real(dp), intent(in)  :: xi_deg(:), yi_deg(:), xc_deg(:), yc_deg(:)
    real(dp), intent(out) :: lon_mid_deg, lat_mid_deg
    real(dp) :: lat_all_min, lat_all_max
    lat_all_min = min(minval(yi_deg), minval(yc_deg))
    lat_all_max = max(maxval(yi_deg), maxval(yc_deg))
    lat_mid_deg = 0.5_dp * (lat_all_min + lat_all_max)
    call circular_midpoint([xi_deg, xc_deg], lon_mid_deg)
  end subroutine bbox_midpoint

  !> @brief Circular mean of longitudes in degrees.
  !!
  !! Computes mean direction on circle via sum of sines/cosines and
  !! returns angle in degrees in [-180, 180].
  !!
  !! @param[in]  lon_deg_all(:) Longitudes [deg]
  !! @param[out] lon_mid_deg Circular mean longitude [deg]
  subroutine circular_midpoint(lon_deg_all, lon_mid_deg)
    real(dp), intent(in)  :: lon_deg_all(:)
    real(dp), intent(out) :: lon_mid_deg
    real(dp) :: th, Cx, Cy
    integer :: k, n
    n = size(lon_deg_all)
    Cx = 0.0_dp; Cy = 0.0_dp
    do k = 1, n
      th = deg2rad * lon_deg_all(k)
      Cx = Cx + cos(th)
      Cy = Cy + sin(th)
    end do
    if (Cx == 0.0_dp .and. Cy == 0.0_dp) then
      lon_mid_deg = 0.0_dp
    else
      lon_mid_deg = rad2deg * atan2(Cy, Cx)
    end if
  end subroutine circular_midpoint

  !> @brief Wrap angle to [-pi, pi] in radians.
  !!
  !! @param[in] theta Angle [rad]
  !! @return Real(dp) wrapped angle in [-pi, pi] [rad].
  pure real(dp) function wrap_pi(theta) result(t)
    real(dp), intent(in) :: theta
    t = modulo(theta + pi, two_pi) - pi
  end function wrap_pi

end module module_geodesy
