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

! -------------------------------------------------------------------
! All control characters in the ASCII table (see www.asciitable.com).
! -------------------------------------------------------------------
character(len=1), public, parameter :: ASCII_NUL = achar(int(z'00')) !! Null
character(len=1), public, parameter :: ASCII_SOH = achar(int(z'01')) !! Start of heading
character(len=1), public, parameter :: ASCII_STX = achar(int(z'02')) !! Start of text
character(len=1), public, parameter :: ASCII_ETX = achar(int(z'03')) !! End of text
character(len=1), public, parameter :: ASCII_EOT = achar(int(z'04')) !! End of transmission
character(len=1), public, parameter :: ASCII_ENQ = achar(int(z'05')) !! Enquiry
character(len=1), public, parameter :: ASCII_ACK = achar(int(z'06')) !! Acknowledge
character(len=1), public, parameter :: ASCII_BEL = achar(int(z'07')) !! Bell
character(len=1), public, parameter :: ASCII_BS  = achar(int(z'08')) !! Backspace
character(len=1), public, parameter :: ASCII_TAB = achar(int(z'09')) !! Horizontal tab
character(len=1), public, parameter :: ASCII_LF  = achar(int(z'0A')) !! NL line feed, new line
character(len=1), public, parameter :: ASCII_VT  = achar(int(z'0B')) !! Vertical tab
character(len=1), public, parameter :: ASCII_FF  = achar(int(z'0C')) !! NP form feed, new page
character(len=1), public, parameter :: ASCII_CR  = achar(int(z'0D')) !! Carriage return
character(len=1), public, parameter :: ASCII_SO  = achar(int(z'0E')) !! Shift out
character(len=1), public, parameter :: ASCII_SI  = achar(int(z'0F')) !! Shift in
character(len=1), public, parameter :: ASCII_DLE = achar(int(z'10')) !! Data link escape
character(len=1), public, parameter :: ASCII_DC1 = achar(int(z'11')) !! Device control 1
character(len=1), public, parameter :: ASCII_DC2 = achar(int(z'12')) !! Device control 2
character(len=1), public, parameter :: ASCII_DC3 = achar(int(z'13')) !! Device control 3
character(len=1), public, parameter :: ASCII_DC4 = achar(int(z'14')) !! Device control 4
character(len=1), public, parameter :: ASCII_NAK = achar(int(z'15')) !! Negative acknowledge
character(len=1), public, parameter :: ASCII_SYN = achar(int(z'16')) !! Synchronous idle
character(len=1), public, parameter :: ASCII_ETB = achar(int(z'17')) !! End of transmission block
character(len=1), public, parameter :: ASCII_CAN = achar(int(z'18')) !! Cancel
character(len=1), public, parameter :: ASCII_EM  = achar(int(z'19')) !! End of medium
character(len=1), public, parameter :: ASCII_SUB = achar(int(z'1A')) !! Substitute
character(len=1), public, parameter :: ASCII_ESC = achar(int(z'1B')) !! Escape
character(len=1), public, parameter :: ASCII_FS  = achar(int(z'1C')) !! File separator
character(len=1), public, parameter :: ASCII_GS  = achar(int(z'1D')) !! Group separator
character(len=1), public, parameter :: ASCII_RS  = achar(int(z'1E')) !! Record separator
character(len=1), public, parameter :: ASCII_US  = achar(int(z'1F')) !! Unit separator
character(len=1), public, parameter :: ASCII_DEL = achar(int(z'7F')) !! Delete

! -------------------------------------------------------------------
! Constant character sequences
! -------------------------------------------------------------------
character(len=*), public, parameter :: digits = '0123456789'
character(len=*), public, parameter :: letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz" 
character(len=*), public, parameter :: uppercase_letters = letters(1:26) !! A .. Z
character(len=*), public, parameter :: lowercase_letters = letters(27:) !! a .. z

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
