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

module module_constants

use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan,  &
                  ieee_is_nan, ieee_positive_inf
use module_types

implicit none

! ------------------------------------------------------------
! Logical constants (optional; Fortran already has .True./.False.)
! ------------------------------------------------------------
logical, parameter     :: True    = .True.
logical, parameter     :: False   = .False.

! ------------------------------------------------------------
! Mathematical constants
! ------------------------------------------------------------
real(dp), parameter    :: minus   =-1.0_dp
real(dp), parameter    :: zero    = 0.0_dp
real(dp), parameter    :: one     = 1.0_dp
real(dp), parameter    :: two     = 2.0_dp
real(dp), parameter    :: three   = 3.0_dp
real(dp), parameter    :: four    = 4.0_dp
real(dp), parameter    :: five    = 5.0_dp
real(dp), parameter    :: six     = 6.0_dp
real(dp), parameter    :: seven   = 7.0_dp
real(dp), parameter    :: eight   = 8.0_dp
real(dp), parameter    :: nine    = 9.0_dp
real(dp), parameter    :: ten     = 10.0_dp
real(dp), parameter    :: half    = 0.5_dp
real(dp), parameter    :: quarter = 0.25_dp
real(dp), parameter    :: hundred = 100.0_dp

real(dp), parameter    :: pi      = 3.1415926535897932384626433832795_dp
real(dp), parameter    :: two_pi  = 2.0_dp*pi
real(dp), parameter    :: half_pi = 0.5_dp*pi
real(dp), parameter    :: inv_pi  = 1.0D0/pi
real(dp), parameter    :: Euler   = 2.7182818284590452353602874713527_dp

real(dp), parameter    :: deg2rad = pi/180.0_dp
real(dp), parameter    :: rad2deg = 180.0_dp/pi

complex(dp), parameter :: zj      = (0.0_dp, 1.0_dp)

! ------------------------------------------------------------
! NaN and Infinity (cannot be parameters)
! ------------------------------------------------------------
real(sp)               :: nan4    ! = 0.0_sp/0.0_sp
real(sp)               :: inf4    ! = 1.0_sp/0.0_sp
real(dp)               :: nan     ! = 0.0_dp/0.0_dp
real(dp)               :: inf     ! = 1.0_dp/0.0_dp

! ... Physical constants
! ...
type type_constants
  real(dp)             :: Earth_Gravity     = 9.80665_dp     ! Standard gravity [m/s^2]
  real(dp)             :: Earth_Omega       = 7.2921159D-5   ! Earth's angular velocity [rad/s]
  real(dp)             :: Earth_Radius      = 6371315.0_dp   ! Earth radius used in models [m]
  real(dp)             :: Earth_Mass        = 5.9722D24      ! Earth mass [kg]
  real(dp)             :: Earth_Declination = 23.446_dp      ! Earth orb. decl. deg
  real(dp)             :: Solar_Constant    = 1365.2_dp      ! W/m2
  real(dp)             :: Gravity_constant  = 6.67430e-11_dp ! Universal grav const [m^3/kg/s^2]
  real(dp)             :: sigma_SB          = 5.670374419D-8 ! Stefan-Boltzmann [W/m^2/K^4]
  real(dp)             :: k_Boltzmann       = 1.380649D-23   ! Boltzmann constant [J/K]
  real(dp)             :: R_gas             = 8.314462618_dp ! Universal gas constant [J/mol/K]
  real(dp)             :: Avogadro_Number   = 6.02214076D23  ! Avogadro number [1/mol]
end type type_constants


type(type_Constants)                          :: constants

contains
  ! ...
  ! ===================================================================
  ! ===================================================================
  ! ...
  subroutine set_nan()
    ! ... Initialize NaN and values
    nan4 = ieee_value(0.0_sp, ieee_quiet_nan)
    nan   = ieee_value(0.0_dp, ieee_quiet_nan)
  end subroutine set_nan
  ! ...
  ! ===================================================================
  ! ...
  subroutine set_inf()
    ! ... Initialize Infinity values
    inf4 = ieee_value(0.0_sp, ieee_positive_inf)
    inf   = ieee_value(0.0_dp, ieee_positive_inf)
  end subroutine set_inf

end module module_constants
