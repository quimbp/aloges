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

module aloges

use module_types
use module_constants
use module_nc
use module_math
use module_tools
use module_time
use module_grid
use module_lineargs
use module_help
use module_trajectories
use module_minpack

implicit none

character(len=4)  :: ALOGES_VERSION = "1.0"
character(len=20) :: ALOGES_AUTHOR  = "Joaquim Ballabrera"
character(len=20) :: ALOGES_DATE    = "April, 2022"

end module aloges

