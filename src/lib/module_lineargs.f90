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

module module_lineargs

use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan, ieee_positive_inf
use module_types, only: sp,dp,maxlen
use module_constants, only: nan,nan4
use module_tools

implicit none

integer                                          :: lineargs_nargs = 0
integer, dimension(:), allocatable               :: lineargs_used
character(len=maxlen), dimension(:), allocatable :: lineargs_val

interface linearg
  module procedure argflg,argflt,argint,arglst,argdbl
end interface linearg


contains
! ...
! ========================================================================
! ...
integer function lineargs () result(na)

integer i,j,jj,nn,iu
logical read_from_file
character(len=maxlen) word,line

na = iargc()
if (na.eq.0) return

read_from_file = .False.
call getarg(1,word)
if (word(1:3).eq.'--o') read_from_file = .True.
if (word(1:3).eq.'--O') read_from_file = .True.

if (read_from_file) then
! ... Read execution options from a file
! ...
  call getarg(2,word)
  iu = unitfree()
  open(iu,file=word,status='old')
  nn = numlines(iu,'ascii')
  na = 0
  do i=1,nn
    read(iu,'(A)') line
    if (line(1:1).ne.'#') na = na + numwords(line)
  enddo
  lineargs_nargs = na
  allocate (lineargs_val(na))
  allocate (lineargs_used(na))
  lineargs_used(:)  = 0

  rewind(iu)
  jj = 0
  do i=1,nn
    read(iu,'(A)') line
    if (line(1:1).ne.'#') then
      do j=1,numwords(line)
        jj = jj + 1
        call line_word(line,j,lineargs_val(jj))
      enddo
    endif
  enddo
  if (jj.NE.na) call crash ('something wrong in lineargs_ini')
  close(iu)
else
! ... Read execution options from command line
! ...
  lineargs_nargs = na
  allocate (lineargs_val(na))
  allocate (lineargs_used(na))
  lineargs_used(:)  = 0
  lineargs_val(:) = ''
  do i=1,lineargs_nargs
    call getarg(i,lineargs_val(i))
  enddo
endif

return
end function lineargs
! ...
! ========================================================================
! ...
subroutine argflg (label,flag)

character(len=*), intent(in)     :: label
logical, intent(inout)           :: flag

integer n,k

n = len_trim(label) 

do k=1,lineargs_nargs
 if (lineargs_used(k).EQ.0) then
 if (uppercase(lineargs_val(k)(:n)).eq.uppercase(label(:n))) then 
  lineargs_used(k) = 1
  flag = .True.
  return
 endif
 endif
enddo

return
end subroutine argflg
! ...
! ====================================================================
! ...
subroutine argflt (label,flag,val)
! ... Returns a float value from user arguments

character(len=*), intent(in)        :: label
logical, intent(inout)              :: flag
real(sp), intent(out)               :: val

! ... Local variables
! ...
integer n,k
character lett

if (flag) return

!nan4 = ieee_value(1.0_sp,ieee_quiet_nan)

n = len_trim(label)
DO k=1,lineargs_nargs
 if (lineargs_used(k).EQ.0) then
 if (uppercase(lineargs_val(k)(:n)).eq.uppercase(label(:n))) then 
  if (k+1.GT.lineargs_nargs) call crash('Invalid option')
  if (lineargs_used(k+1).EQ.1) then
    write(6,*) 'Option ',TRIM(lineargs_val(k+1)),' already used.'
    call crash ('Invalid option')
  endif
  flag      = .True.
  lett      = lineargs_val(k+1)(1:1)
!  if (lett.EQ.'N'.OR.lett.EQ.'n') then
!  val     = nan4
!  else 
    read(lineargs_val(k+1),*) val
!  endif
  lineargs_used(k)   = 1
  lineargs_used(k+1) = 1
  return
 endif
 endif
enddo

return
end subroutine argflt
! ...
! ====================================================================
! ...
subroutine argint (label,flag,val)
! ... Returns an integer value from user arguments

character(len=*), intent(in)        :: label
logical, intent(inout)              :: flag
integer, intent(out)                :: val

integer n,k

if (flag) return

n = len_trim(label)
DO k=1,lineargs_nargs
 if (lineargs_used(k).EQ.0) then
 if (uppercase(lineargs_val(k)(:n)).eq.uppercase(label(:n))) then 
  if (k+1.GT.lineargs_nargs) call crash('Invalid option')
  if (lineargs_used(k+1).EQ.1) then
    write(6,*) 'Option ',TRIM(lineargs_val(k+1)),' already used.'
    call crash ('Invalid option')
  endif
  flag      = .True.
  read(lineargs_val(k+1),*) val
  lineargs_used(k)   = 1
  lineargs_used(k+1) = 1
  return
 endif
 endif
enddo

return
end subroutine argint
! ...
! ====================================================================
! ...
subroutine argstr (label,flag,val)

character(len=*), intent(in)        :: label
logical, intent(inout)              :: flag
character(LEN=*), intent(inout)     :: val

integer n,k

if (flag) return
n = len_trim(label)

DO k=1,lineargs_nargs
 if (lineargs_used(k).EQ.0) then
 if (uppercase(lineargs_val(k)(:n)).eq.uppercase(label(:n))) then 
  if (k+1.GT.lineargs_nargs) call crash ('Invalid option')
  if (lineargs_used(k+1).EQ.1) then
    write(6,*) 'Option ',TRIM(lineargs_val(k+1)),' already used.'
    call crash ('Invalid option')
  endif
  flag      = .True.
  val       = TRIM(lineargs_val(k+1))
  lineargs_used(k)   = 1
  lineargs_used(k+1) = 1
  return
 endif
 endif
enddo

return
end subroutine argstr
! ...
! ====================================================================
! ...
subroutine arglst (label,flag,list)

character(LEN=*), intent(in)        :: label
logical, intent(inout)              :: flag
character(LEN=*), intent(inout)     :: list

integer n,j,k
character(len=len(list)) mylist
character(len=maxlen) word

if (flag) return

mylist = list
list = ''
n = len_trim(label)

do k=1,lineargs_nargs
 if (lineargs_used(k).EQ.0) then
 if (uppercase(lineargs_val(k)(:n)).eq.uppercase(label(:n))) then 
  if (k+1.GT.lineargs_nargs) call crash('Invalid option')
  flag = .True.
  lineargs_used(k) = 1

  ! ... Read arguments until the end or till the first -XXXX option
  ! ...
  do j=k+1,lineargs_nargs
    if (lineargs_val(j)(1:1).EQ.'-') EXIT
    if (len_trim(list).eq.0) then
      list = compress(lineargs_val(j))
    else
      word = compress(lineargs_val(j))
      list = trim(list)//' '//trim(word)
    endif
    lineargs_used(j) = 1
  enddo
  return
 endif
 endif
enddo

if (.not.flag) list = mylist

return
end subroutine arglst
! ...
! ====================================================================
! ...
subroutine argdbl (label,flag,val)

character(LEN=*), intent(in)        :: label
logical, intent(inout)              :: flag
REAL(dp), intent(out)               :: val

character(LEN=1) lett
integer n,k

if (flag) return

!nan = ieee_value(1.0_dp,ieee_quiet_nan)

n = len_trim(label)
DO k=1,lineargs_nargs
 if (lineargs_used(k).EQ.0) then
 if (uppercase(lineargs_val(k)(:n)).eq.uppercase(label(:n))) then
  if (k+1.GT.lineargs_nargs) call crash('Invalid option')
  if (lineargs_used(k+1).EQ.1) then
    write(6,*) 'Option ',TRIM(lineargs_val(k+1)),' already used.'
    call crash ('Invalid option')
  endif
  flag      = .True.
  lett      = lineargs_val(k+1)(1:1)
!  if (lett.EQ.'N'.OR.lett.EQ.'n') then
!    val     = nan
!  else 
    read(lineargs_val(k+1),*) val
!  endif
  lineargs_used(k)   = 1
  lineargs_used(k+1) = 1
  return
 endif
 endif
enddo

return
end subroutine argdbl 
! ...
! ====================================================================
! ...
subroutine checkopts 

integer i

if (lineargs_nargs.EQ.0) return

if (SUM(lineargs_used).NE.lineargs_nargs) then
  write(*,'(A)',ADVANCE="no") 'Unknown options : '
  DO i=1,lineargs_nargs
    if (lineargs_used(i).EQ.0) write(*,'(A)',ADVANCE="no") TRIM(lineargs_val(i))//' '
  enddo
  write(*,*)
  call crash('Invalid options')
endif

end subroutine checkopts
! ...
! =======================================================================
! ...
subroutine arglast (file)

character(LEN=*), intent(out)  :: file

if (lineargs_used(lineargs_nargs).EQ.0) then
  file = lineargs_val(lineargs_nargs)
  lineargs_used(lineargs_nargs) = 1
else
  call crash('Last argument already used in arglast')
endif

return
end subroutine arglast
! ...
! ========================================================================
! ...
end module module_lineargs
