! ======================================================================== !
! ALOGES PROJECT                                                           !
! Quim Ballabrera, April 2022                                              !
! Institut de Ciencies del Mar, CSIC                                       !
! Last Modified: 2025-10-19                                                !
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
! - print_colored                                                          !
! -------------------------------------------------------------------------!

module module_color

implicit none
character(len=*), parameter :: esc = char(27)

! ... Reset and text attributes
! ...
character(len=*), parameter :: reset       = esc//'[0m'
character(len=*), parameter :: bold        = esc//'[1m'
character(len=*), parameter :: dim         = esc//'[2m'
character(len=*), parameter :: italic      = esc//'[3m'
character(len=*), parameter :: underline   = esc//'[4m'
character(len=*), parameter :: blink       = esc//'[5m'
character(len=*), parameter :: reverse     = esc//'[7m'
character(len=*), parameter :: hidden      = esc//'[8m'
character(len=*), parameter :: strike      = esc//'[9m'

! ... Standard colors
! ...
character(len=*), parameter :: black       = esc//'[30m'
character(len=*), parameter :: red         = esc//'[31m'
character(len=*), parameter :: green       = esc//'[32m'
character(len=*), parameter :: yellow      = esc//'[33m'
character(len=*), parameter :: blue        = esc//'[34m'
character(len=*), parameter :: magenta     = esc//'[35m'
character(len=*), parameter :: cyan        = esc//'[36m'
character(len=*), parameter :: white       = esc//'[37m'

! ... Bright colors
! ...
character(len=*), parameter :: bright_black   = esc//'[90m'
character(len=*), parameter :: bright_red     = esc//'[91m'
character(len=*), parameter :: bright_green   = esc//'[92m'
character(len=*), parameter :: bright_yellow  = esc//'[93m'
character(len=*), parameter :: bright_blue    = esc//'[94m'
character(len=*), parameter :: bright_magenta = esc//'[95m'
character(len=*), parameter :: bright_cyan    = esc//'[96m'
character(len=*), parameter :: bright_white   = esc//'[97m'

! ... Backgrounds
! ...
character(len=*), parameter :: bg_black    = esc//'[40m'
character(len=*), parameter :: bg_red      = esc//'[41m'
character(len=*), parameter :: bg_green    = esc//'[42m'
character(len=*), parameter :: bg_yellow   = esc//'[43m'
character(len=*), parameter :: bg_blue     = esc//'[44m'
character(len=*), parameter :: bg_magenta  = esc//'[45m'
character(len=*), parameter :: bg_cyan     = esc//'[46m'
character(len=*), parameter :: bg_white    = esc//'[47m'

! ... Bright backgrounds
! ...
character(len=*), parameter :: bg_bright_black   = esc//'[100m'
character(len=*), parameter :: bg_bright_red     = esc//'[101m'
character(len=*), parameter :: bg_bright_green   = esc//'[102m'
character(len=*), parameter :: bg_bright_yellow  = esc//'[103m'
character(len=*), parameter :: bg_bright_blue    = esc//'[104m'
character(len=*), parameter :: bg_bright_magenta = esc//'[105m'
character(len=*), parameter :: bg_bright_cyan    = esc//'[106m'
character(len=*), parameter :: bg_bright_white   = esc//'[107m'

contains
! ...
! =====================================================================
! =====================================================================
! ...
  subroutine print_colored(msg, color)
    character(len=*), intent(in) :: msg, color

    print '(A)', trim(color)//msg//reset

  end subroutine print_colored

end module module_color
