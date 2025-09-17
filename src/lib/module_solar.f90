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
!                                                                          !
! List of routines                                                         !
! - Sunrise_and_Sunset                                                     !
! -------------------------------------------------------------------------!

module module_solar

  use module_types
  use module_tools
  use module_time

  implicit none

  contains
  ! ...
  ! =====================================================================
  ! =====================================================================
  ! ...
  subroutine SunRise_and_SunSet (lon,lat,date,hSunrise,hSunset)
    ! ... Fortran 90 subroutines based on the C code distributed by Mike Chirico:
 
    ! ... http://souptonuts.sourceforge.net/code/sunrise.c.html
    ! ... Copyright (GPL) 2004   Mike Chirico mchirico@comcast.net
    ! ... Updated: Sun Nov 28 15:15:05 EST 2004
    ! ... Program adapted by Mike Chirico mchirico@comcast.net
    ! ... Reference:
    ! ... http://prdownloads.sourceforge.net/souptonuts/working_with_time.tar.gz?download
    ! ... http://www.srrb.noaa.gov/highlights/sunrise/sunrise.html

    real(dp), intent(in)                         :: lon
    real(dp), intent(in)                         :: lat
    type(type_date), intent(in)                  :: date
    real(dp), intent(out)                        :: hSunrise
    real(dp), intent(out)                        :: hSunset

    ! ... Local variables:
    ! ...
    integer JD

    JD = julday(date%year,date%month,date%day) 
    hSunrise = anint(SunriseUTC(JD,-lon,lat) * 60.0D0)       ! seconds
    hSunset  = anint(SunsetUTC(JD,-lon,lat) * 60.0D0)        ! seconds

    return
  end subroutine Sunrise_and_Sunset
  ! ...
  ! ===================================================================
  ! ...
  real(dp) function SunriseUTC(JD,lon,lat) result(TSunrise)

    integer, intent(in)                            :: JD
    real(dp), intent(in)                           :: lon
    real(dp), intent(in)                           :: lat

    ! ... Local variables
    ! ...
    real(dp) t,EqTime,SolarDec,HourAngle,Delta,TimeDiff,TimeUTC,newt

    t         = TimeJulianCent(dble(JD))
    EqTime    = EquationOfTime(t)
    SolarDec  = SunDeclination(t)
    HourAngle = HourAngleSunrise(lat,SolarDec)
    Delta     = lon - rad2deg*HourAngle
    TimeDiff  = 4.0D0*Delta
    TimeUTC   = 720.0D0 + TimeDiff - EqTime     ! Minutes

    newt      = TimeJulianCent(JDFromJulianCent(t) + timeUTC/1440.0D0)
    EqTime    = EquationOfTime(newt)
    SolarDec  = SunDeclination(newt)
    HourAngle = HourAngleSunrise(lat,SolarDec);
    Delta     = lon - rad2deg*HourAngle
    TimeDiff  = 4*Delta
    TSunrise  = 720.0D0 + timeDiff - EqTime

  end function SunriseUTC
  ! ...
  ! ===================================================================
  ! ...
  real(dp) function SunsetUTC(JD,lon,lat) result(TSunset)

    integer, intent(in)                            :: JD
    real(dp), intent(in)                           :: lon
    real(dp), intent(in)                           :: lat

    ! ... Local variables
    ! ...
    real(dp) t,EqTime,SolarDec,HourAngle,Delta,TimeDiff,TimeUTC,newt

    t         = TimeJulianCent(dble(JD))
    EqTime    = EquationOfTime(t)
    SolarDec  = SunDeclination(t)
    HourAngle = HourAngleSunset(lat,SolarDec)
    Delta     = lon - rad2deg*HourAngle
    timeDiff  = 4.0D0 * delta
    timeUTC   = 720.0D0 + timeDiff - eqTime     ! Minutes	

    newt      = TimeJulianCent(JDFromJulianCent(t) + timeUTC/1440.0D0)
    eqTime    = EquationOfTime(newt)
    SolarDec  = SunDeclination(newt)
    HourAngle = HourAngleSunset(lat,SolarDec)
    Delta     = lon - rad2deg*HourAngle
    TimeDiff  = 4.0D0 * delta;
    TSunset   = 720.0D0 + timeDiff - EqTime;    ! Minutes

  end function SunsetUTC
  ! ...
  ! ===================================================================
  ! ...
  real(dp) function TimeJulianCent(JD) result(t)
  
    real(dp), intent(in)                         :: JD

    t = (JD - 2451545.0D0)/36525.0D0

  end function TimeJulianCent
  ! ...
  ! ===================================================================
  ! ...
  real(dp) function JDFromJulianCent(t) result(JD)

    real(dp)                                     :: t

    JD = t*36525.0D0 + 2451545.0D0

  end function JDFromJulianCent
  ! ...
  ! ===================================================================
  ! ...
  real(dp) function EquationOfTime(t) result(Etime)

    real(dp), intent(in)                         :: t

    ! ... Local variables
    ! ...
    real(dp) eps,l0,e,m,y,sin2l0,sinm,cos2l0,sin4l0,sin2m

    e   = OrbitEccentricity(t)                 ! Unitless
    eps = deg2rad*ObliquityCorrection(t)       ! Radians
    l0  = deg2rad*GeomMeanLongSun(t)           ! Radians
    m   = deg2rad*GeomMeanAnomalySun(t)        ! Radians

    y      = tan(0.5D0*eps)**2
    sin2l0 = sin(2.0D0*l0)
    sinm   = sin(m)
    cos2l0 = cos(2.0D0*l0)
    sin4l0 = sin(4.0D0*l0)
    sin2m  = sin(2.0D0*m)

    Etime =   y*sin2l0               &
            - 2.0D0*e*sinm           &
            + 4.0D0*e*y*sinm*cos2l0  &
            - 0.5D0*y*y*sin4l0       &
            - 1.25D0*e*e*sin2m

    Etime = 4.0D0*rad2deg*Etime   ! Minutes of time

  end function EquationOfTime
  ! ...
  ! ===================================================================
  ! ...
  real(dp) function ObliquityCorrection(t) result(e)  ! Result in degrees

    real(dp), intent(in)                         :: t

    ! ... Local variables
    ! ...
    real(dp) e0,omega

    omega = 125.04D0 - 1934.136D0*t;

    e0 = MeanObliquity(t);
    e  = e0 + 0.00256D0*cos(deg2rad*omega);

  end function ObliquityCorrection
  ! ...
  ! ===================================================================
  ! ...
  real(dp) function MeanObliquity (t) result(e0)   ! Result in degrees

    real(dp), intent(in)                         :: t

    ! ... Local variables
    ! ...
    real(dp) seconds

    seconds = 21.448D0 - t*(46.8150D0 + t*(0.00059D0 - t*0.001813D0));
    e0 = 23.0D0 + (26.0D0 + seconds/60.0D0)/60.0D0;

  end function MeanObliquity
  ! ...
  ! ===================================================================
  ! ...
  real(dp) function GeomMeanLongSun(t) result(L)   ! in degrees

    real(dp), intent(in)                         :: t

    L = 280.46646D0 + t*(36000.76983D0 + 0.0003032D0*t)
    do while(L.gt.360.0D0)
      L = L - 360.0D0
    enddo
    do while(L.lt.0)
      L = L + 360.0D0
    enddo

  end function GeomMeanLongSun
  ! ...
  ! ===================================================================
  ! ...
  real(dp) function OrbitEccentricity(t) result(e)   ! Unitless

    real(dp), intent(in)                         :: t

    e = 0.016708634D0 - t*(0.000042037D0 + 0.0000001267D0*t)

  end function OrbitEccentricity
  ! ...
  ! ===================================================================
  ! ...
  real(dp) function GeomMeanAnomalySun(t) result(M) ! Degrees

    real(dp), intent(in)                         :: t

    M = 357.52911D0 + t*(35999.05029D0 - 0.0001537*t)

  end function GeomMeanAnomalySun
  ! ...
  ! ===================================================================
  ! ...
  real(dp) function SunDeclination(t) result(theta) ! Degrees

    real(dp), intent(in)                         :: t

    ! ... Local variables
    ! ...
    real(dp) e,lambda,sint

    e      = deg2rad*ObliquityCorrection(t)  ! Radians
    lambda = deg2rad*SunApparentLong(t)      ! Radians

    sint  = sin(e) * sin(lambda)
    theta = rad2deg*asin(sint)

  end function SunDeclination
  ! ...
  ! ===================================================================
  ! ...
  real(dp) function SunApparentLong(t) result(lambda)   ! Degrees

    real(dp), intent(in)                         :: t

    ! ... Local variables
    ! ...
    real(dp) o,omega

    o      = SunTrueLong(t)
    omega  = deg2rad*(125.04D0 - 1934.136D0*t)
    lambda = o - 0.00569D0 - 0.00478D0*sin(omega)

  end function SunApparentLong
  ! ...
  ! ===================================================================
  ! ...
  real(dp) function SunTrueLong(t)  result(o)   ! Degrees

    real(dp), intent(in)                         :: t

    ! ... Local variables
    ! ...
    real(dp) l0,c

    l0 = GeomMeanLongSun(t)
    c  = SunEqOfCenter(t)

    O = l0 + c

  end function SunTrueLong
  ! ...
  ! ===================================================================
  ! ...
  real(dp) function SunEqOfCenter(t) result(C)   ! Degrees

    real(dp), intent(in)                         :: t

    ! ... Local variables
    ! ...
    real(dp) mrad,sinm,sin2m,sin3m

    mrad  = deg2rad*GeomMeanAnomalySun(t)
    sinm  = sin(mrad);
    sin2m = sin(mrad+mrad);
    sin3m = sin(mrad+mrad+mrad);

    C =   sinm*(1.914602D0 - t*(0.004817D0 + 0.000014D0*t)) &
        + sin2m*(0.019993D0 - 0.000101D0*t) &
        + sin3m*0.000289D0

  end function SunEqOfCenter
  ! ...
  ! ===================================================================
  ! ...
  real(dp) function HourAngleSunrise(lat,SolarDec) result(HA)  ! In radians

    real(dp), intent(in)                         :: lat
    real(dp), intent(in)                         :: SolarDec

    ! ... Local variables
    ! ...
    real(dp) latRad,sdRad

    latRad = deg2rad*lat;
    sdRad  = deg2rad*SolarDec;

    HA = acos( cos(deg2rad*90.833)/(cos(latRad)*cos(sdRad)) - tan(latRad)*tan(sdRad) )

  end function HourAngleSunrise
  ! ...
  ! ===================================================================
  ! ...
  real(dp) function HourAngleSunset(lat,SolarDec) result(HA)  ! In radians

    real(dp), intent(in)                         :: lat
    real(dp), intent(in)                         :: SolarDec

    ! ... Local variables
    ! ...
    real(dp) latRad,sdRad

    latRad = deg2rad*lat;
    sdRad  = deg2rad*SolarDec;

    HA = -acos( cos(deg2rad*90.833)/(cos(latRad)*cos(sdRad)) - tan(latRad)*tan(sdRad) )

  end function HourAngleSunset
  ! ...
  ! =====================================================================
  ! ...
  real(dp) function daily_insolation(lat,day) result(Qsw)

  ! ... Estimate solar longitude from calendar day
  ! ... 
    real(dp), intent(in)                         :: lat   ! Latitude in radians
    real(dp), intent(in)                         :: day   ! Day of year: 1 - 365.25

    ! ... Local variables
    ! ...
    real(dp), parameter                          :: S0   = 1365.2D0    ! W/m2
    real(dp), parameter                          :: Ecc  = 0.017236D0  ! Eccentricity
    real(dp), parameter                          :: Peri = deg2rad*281.37D0  ! Long. perihelion (rad)
    real(dp), parameter                          :: Obli = deg2rad*23.446D0    ! Obliquity
    real(dp), parameter                          :: Year = 365.2422D0  ! Days per year

    real(dp) delta_lambda,beta,wrk,lambda_lon,delta,Ho
    real(dp) coszen

    ! ... Get solar longitude (in radians) from calendar day
    ! ... Calendar is referenced to the vernal equinox (21 March), day 81
    ! ...
    delta_lambda = 2.0D0*pi*(day - 81.0D0)/Year
    beta = sqrt(1-Ecc*Ecc)
    wrk = -2.0D0*((0.500D0*Ecc + 0.125D0*Ecc**3)*(1.0D0+beta)*sin(-Peri) - &
                   0.250D0*Ecc*Ecc*(0.50D0+beta)*sin(-2.0D0*Peri) + &
                   0.125D0*Ecc**3*(1.0D0/3.0D0+beta)*sin(-3.0D0*Peri)) &
          + delta_lambda
    lambda_lon = wrk + Ecc*((2.0D0-0.25*Ecc*Ecc)*sin(wrk-Peri) + &
                            1.25D0*Ecc*sin(2.0D0*(wrk-Peri)) + &
                            (13.0D0/12.0D0)*Ecc*Ecc*sin(3.0D0*(wrk-Peri)))

    ! ... Declination angle of the sun:
    ! ...
    delta = asin(sin(Obli)*sin(lambda_lon))

    ! ... Ho, hour angle at sunrise / sunset
    ! ...
    if (abs(delta)-half_pi+abs(lat).lt.0.0D0) then
      ! ... There is sunset/sunrise
      ! ...
      Ho = acos(-tan(lat)*tan(delta))
    else
      ! ... Check if all day or night
      ! ...
      if (lat*delta.gt.0.0D0) then
        Ho = pi
      else
        Ho = 0.0D0
      endif
    endif

    ! ... Integral from sunrise to sunset:
    ! ...
    coszen = Ho*sin(lat)*sin(delta) + cos(lat)*cos(delta)*sin(Ho)   

    ! ... Compute insolation:
    ! ...
    Qsw = S0/pi*coszen*(1.0D0+Ecc*cos(lambda_lon-Peri))**2/(1.0D0-Ecc**2)**2
    

  end function daily_insolation
  ! ...
  ! =====================================================================
  ! ...
end module module_solar
