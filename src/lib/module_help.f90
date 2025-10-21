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
! - help_progname                                                          !
! - help_version                                                           !
! - help_command                                                           !
! - help_example                                                           !
! - help_option                                                            !
! - help_write                                                             !
! -------------------------------------------------------------------------!

module module_help

use module_tools
use iso_fortran_env

implicit none (type, external)
integer, parameter                                 :: MAXOPTIONS = 50

type type_help
  integer                                          :: numoptions = 0

  character(len=80)                                :: progname = ''
  character(len=10)                                :: version = 'v0.0'
  character(len=80)                                :: author = 'anonymous'
  character(len=10000)                             :: summary = ''
  character(len=1000)                              :: command = ''
  character(len=1000)                              :: example = ''

  character(len=80), dimension(MAXOPTIONS)         :: option
  character(len=180), dimension(MAXOPTIONS)        :: description
  character(len=80), dimension(MAXOPTIONS)         :: default

  contains

    procedure        :: set_progname   => help_progname
    procedure        :: set_version    => help_version
    procedure        :: set_command    => help_command
    procedure        :: set_example    => help_example
    procedure        :: set_summary    => help_example
    procedure        :: add_option     => help_option 
    procedure        :: write          => help_write

end type type_help

type(type_help) help

contains
! ...
! =============================================================================
! ...
subroutine help_progname(HLP,word)

class(type_help), intent(inout)        :: HLP
character(len=*), intent(in)           :: word

HLP%progname = compress(word)

return
end subroutine help_progname
! ...
! =============================================================================
! ...
subroutine help_version(HLP,word)

class(type_help), intent(inout)        :: HLP
character(len=*), intent(in)           :: word

HLP%version = compress(word)

return
end subroutine help_version
! ...
! =============================================================================
! ...
subroutine help_command(HLP,word)

class(type_help), intent(inout)        :: HLP
character(len=*), intent(in)           :: word

HLP%command = trim(word)
if (HLP%command(1:1).EQ.'>') HLP%command(1:1) = ''
HLP%command = compress(HLP%command)

return
end subroutine help_command
! ...
! =============================================================================
! ...
subroutine help_example(HLP,word)

class(type_help), intent(inout)        :: HLP
character(len=*), intent(in)           :: word

HLP%example = trim(word)
if (HLP%example(1:1).EQ.'>') HLP%example(1:1) = ''
HLP%example = compress(HLP%example)

return
end subroutine help_example
! ...
! =============================================================================
! ...
subroutine help_summary(HLP,word)

class(type_help), intent(inout)        :: HLP
character(len=*), intent(in)           :: word

HLP%summary = compress(word)

return
end subroutine help_summary
! ...
! =============================================================================
! ...
subroutine help_author(HLP,word)

class(type_help), intent(inout)        :: HLP
character(len=*), intent(in)           :: word

HLP%author = compress(word)

return
end subroutine help_author
! ...
! =============================================================================
! ...
subroutine help_option(HLP,option,description,default)

class(type_help), intent(inout)        :: HLP
character(len=*), intent(in)           :: option,description,default

if (len_trim(option).EQ.0) return
if (HLP%numoptions.eq.MAXOPTIONS) then
  call crash('Increase the value of MAXOPTIONS in COSMO_ROOT/src/lib/help.f90')
endif

! ... New option
! ...
HLP%numoptions = HLP%numoptions + 1

HLP%option(HLP%numoptions)      = trim(option)

if (len_trim(option).GT.0) then
  HLP%description(HLP%numoptions) = trim(description)
ELSE
  HLP%description(HLP%numoptions) = ''
endif

if (len_trim(option).GT.0) then
  HLP%default(HLP%numoptions)     = trim(default)
ELSE
  HLP%default(HLP%numoptions)     = ''
endif

return
end subroutine help_option
! ...
! =============================================================================
! ...
subroutine help_write(HLP)

class(type_help), intent(in)           :: HLP

! ... Local variables
! ...
integer i

write(*,*) '======================================================================='
write(*,*) 'Program: ' // trim(HLP%progname)
write(*,*) 'Version: ' // trim(HLP%version)
write(*,*) 'Author : ' // trim(HLP%author)

if (len_trim(HLP%summary).GT.0) then
  write(*,*)
  call say('Summary: '//trim(HLP%summary))
endif
if (len_trim(HLP%command).GT.0) then
  write(*,*)
  write(*,*) 'Command: '
  write(*,*) '> '//trim(HLP%command)
endif
if (HLP%numoptions.GT.0) then
  write(*,*)
  write(*,*) 'Options:'
  do i=1,HLP%numoptions
    write(*,'(T5,A)') trim(HLP%option(i))
    call say(HLP%description(i),30)
    if (len_trim(HLP%default(i)).GT.0) then
      write(*,'(T30,"Default: ",A)') trim(HLP%default(i))
    endif
  enddo
endif
if (len_trim(HLP%example).GT.0) then
  write(*,*)
  write(*,*) 'Example: '
  write(*,*) '> '//trim(HLP%example)
endif
write(*,*) '======================================================================='
  
STOP
end subroutine help_write
! ...
! =============================================================================
! ...
subroutine header(HLP)

class(type_help), intent(in)           :: HLP

write(*,*) '======================================================================='
write(*,*) 'Program: ' // trim(HLP%progname)
write(*,*) 'Version: ' // trim(HLP%version)
write(*,*) 'By ' // trim(HLP%author)
write(*,*) 'Compiler version: ', compiler_version()
write(*,*) 'Compiler options: ', compiler_options()
write(*,*) '======================================================================='

end subroutine header
! ...
! =============================================================================
! ...

end module module_help
