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
use module_types, only : dp

implicit none (type, external)


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
! Constant end-status symbols
! -------------------------------------------------------------------
character(len=*), parameter :: ok_check_symbol = '✔'  
character(len=*), parameter :: error_x_symbol  = '❌' 

! ----  
! Symbols units and notation (UTF-8)  
! ----  
type, public :: type_symbol
  character(len=1)         :: minute   = achar(39)                          !! '  
  character(len=1)         :: second   = achar(34)                          !! "  
  character(len=2)         :: degree   = achar(194)//achar(176)             !! °  
  character(len=2)         :: square   = achar(194)//achar(178)             !! ² 
  character(len=2)         :: cube     = achar(194)//achar(179)             !! ³
  character(len=3)         :: forall   = achar(226)//achar(136)//achar(128) !! ∀  U+2200
  character(len=3)         :: exists   = achar(226)//achar(136)//achar(131) !! ∃  U+2203
  character(len=3)         :: partial  = achar(226)//achar(136)//achar(130) !! ∂  U+2202
  character(len=3)         :: emptyset = achar(226)//achar(136)//achar(133) !! ∅  U+2205
  character(len=3)         :: nabla    = achar(226)//achar(136)//achar(135) !! ∇  U+2207
  character(len=3)         :: isin     = achar(226)//achar(136)//achar(136) !! ∈  U+2208
  character(len=3)         :: notin    = achar(226)//achar(136)//achar(137) !! ∉  U+2209
  character(len=3)         :: prod     = achar(226)//achar(136)//achar(143) !! ∏  U+220F
  character(len=3)         :: sum      = achar(226)//achar(136)//achar(145) !! ∑  U+2211
  character(len=3)         :: minus    = achar(226)//achar(136)//achar(146) !! −  U+2212
  character(len=3)         :: radic    = achar(226)//achar(136)//achar(154) !! √  U+221A
  character(len=3)         :: infin    = achar(226)//achar(136)//achar(158) !! ∞  U+221E
  character(len=3)         :: and      = achar(226)//achar(136)//achar(167) !! ∧  U+2227
  character(len=3)         :: or       = achar(226)//achar(136)//achar(168) !! ∨  U+2228
  character(len=3)         :: cap      = achar(226)//achar(136)//achar(169) !! ∩  U+2229
  character(len=3)         :: cup      = achar(226)//achar(136)//achar(170) !! ∪  U+222A
  character(len=3)         :: integral = achar(226)//achar(136)//achar(171) !! ∫  U+222B
  character(len=3)         :: therefore= achar(226)//achar(136)//achar(180) !! ∴  U+2234
  character(len=3)         :: sim      = achar(226)//achar(136)//achar(188) !! ∼  U+223C
  character(len=3)         :: approx   = achar(226)//achar(137)//achar(136) !! ≈  U+2248
  character(len=3)         :: ne       = achar(226)//achar(137)//achar(160) !! ≠  U+2260
  character(len=3)         :: equiv    = achar(226)//achar(137)//achar(161) !! ≡  U+2261
  character(len=3)         :: le       = achar(226)//achar(137)//achar(164) !! ≤  U+2264
  character(len=3)         :: ge       = achar(226)//achar(137)//achar(165) !! ≥  U+2265
  character(len=3)         :: subset   = achar(226)//achar(138)//achar(130) !! ⊂  U+2282
  character(len=3)         :: supset   = achar(226)//achar(138)//achar(131) !! ⊃  U+2283
  character(len=3)         :: sube     = achar(226)//achar(138)//achar(134) !! ⊆  U+2286
  character(len=3)         :: supe     = achar(226)//achar(138)//achar(135) !! ⊇  U+2287
  character(len=3)         :: oplus    = achar(226)//achar(138)//achar(149) !! ⊕  U+2295
  character(len=3)         :: otimes   = achar(226)//achar(138)//achar(151) !! ⊗  U+2297
  character(len=3)         :: perp     = achar(226)//achar(138)//achar(165) !! ⊥  U+22A5
  character(len=3)         :: sdot     = achar(226)//achar(139)//achar(133) !! ⋅  U+22C5
  character(len=2)         :: micro    = achar(194)//achar(181)             !! µ  
  character(len=2)         :: pm       = achar(194)//achar(177)             !! ±  
  character(len=2)         :: div      = achar(195)//achar(183)             !! ÷ (División)  
  character(len=2)         :: mul      = achar(195)//achar(151)             !! × (Multiplicación)  
  character(len=3)         :: arr_right= achar(226)//achar(134)//achar(146) !! →  
  character(len=3)         :: arr_up   = achar(226)//achar(134)//achar(145) !! ↑  
  character(len=3)         :: arr_down = achar(226)//achar(134)//achar(147) !! ↓  
  character(len=3)         :: bullet   = achar(226)//achar(128)//achar(162) !! •  
end type type_symbol

! -------------------------------------------------------------------
! Constant character sequences
! -------------------------------------------------------------------
character(len=*), public, parameter :: digits = '0123456789'
character(len=*), public, parameter :: letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz" 
character(len=*), public, parameter :: uppercase_letters = letters(1:26) !! A .. Z
character(len=*), public, parameter :: lowercase_letters = letters(27:) !! a .. z

! ... Physical constants
! ...
type, public :: type_constants
  real(dp)             :: Earth_Gravity     = 9.80665_dp     ! Standard gravity [m/s^2]
  real(dp)             :: Earth_Omega       = 7.2921159D-5   ! Earth's angular velocity [rad/s]
  real(dp)             :: Earth_Radius      = 6371008.8_dp   ! Earth radius Wikipedia
  real(dp)             :: Earth_Mass        = 5.9722D24      ! Earth mass [kg]
  real(dp)             :: Earth_Declination = 23.446_dp      ! Earth orb. decl. deg
  real(dp)             :: Solar_Constant    = 1365.2_dp      ! W/m2
  real(dp)             :: Gravity_constant  = 6.67430e-11_dp ! Universal grav const [m^3/kg/s^2]
  real(dp)             :: sigma_SB          = 5.670374419D-8 ! Stefan-Boltzmann [W/m^2/K^4]
  real(dp)             :: k_Boltzmann       = 1.380649D-23   ! Boltzmann constant [J/K]
  real(dp)             :: R_gas             = 8.314462618_dp ! Universal gas constant [J/mol/K]
  real(dp)             :: Avogadro_Number   = 6.02214076D23  ! Avogadro number [1/mol]
end type type_constants

type(type_Constants)                        :: constants

! ----  
! Greek alphabet (UTF-8, lowercase stored)  
! ----  
type, public :: type_greek  
  character(len=2) :: alpha   = achar(206)//achar(177) !! α  
  character(len=2) :: beta    = achar(206)//achar(178) !! β  
  character(len=2) :: gamma   = achar(206)//achar(179) !! γ  
  character(len=2) :: delta   = achar(206)//achar(180) !! δ  
  character(len=2) :: epsilon = achar(206)//achar(181) !! ε  
  character(len=2) :: theta   = achar(206)//achar(184) !! θ  
  character(len=2) :: lambda  = achar(206)//achar(187) !! λ  
  character(len=2) :: mu      = achar(206)//achar(188) !! μ  
  character(len=2) :: pi      = achar(207)//achar(128) !! π  
  character(len=2) :: rho     = achar(207)//achar(129) !! ρ  
  character(len=2) :: sigma   = achar(207)//achar(131) !! σ  
  character(len=2) :: phi     = achar(207)//achar(134) !! φ  
  character(len=2) :: omega   = achar(207)//achar(137) !! ω  
contains  
  procedure        :: uppercase => greek_uppercase
end type type_greek  


contains
  ! ...
  ! ===================================================================
  ! ===================================================================
  ! ...
  pure function nan() result(v)  
    real(dp) :: v  
    v = ieee_value(0.0_dp, ieee_quiet_nan)  
  end function nan
  ! ...
  ! ===================================================================
  ! ...
  pure function inf() result(v)  
    real(dp) :: v  
    v = ieee_value(0.0_dp, ieee_positive_inf)  
  end function inf
  ! ...
  ! ===================================================================
  ! ...
  pure function greek_uppercase(this, lower) result(upper_char)  
    class(type_greek), intent(in) :: this  
    character(len=2), intent(in)  :: lower  
    character(len=2)              :: upper_char  
  
    integer :: b1, b2  
  
    b1 = iachar(lower(1:1))  
    b2 = iachar(lower(2:2))  
  
    ! Most Greek lowercase letters → uppercase  
    if (b1 == 206 .or. b1 == 207) then  
       upper_char = achar(206)//achar(b2 - 32)  
    else  
       upper_char = lower  
    end if  
  end function greek_uppercase

end module module_constants
