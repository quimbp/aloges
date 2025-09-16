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
! List of routines:                                                        !
! - cat                                                                    !
! - compress                                                               !
! - crash                                                                  !
! - filename_split                                                         !
! - get_commandline                                                        !
! - is_numeric                                                             !
! - line_replace                                                           !
! - line_word                                                              !
! - locate                                                                 !
! - lowercase                                                              !
! - ls                                                                     !
! - numlines                                                               !
! - numwords                                                               !
! - unique_elements (integer, real(dp) interface)                          !
! - unitfree                                                               !
! - uppercase                                                              !
! - say                                                                    !
! - strcat                                                                 !
! - token_read                                                             !
! - vprint                                                                 !
! -------------------------------------------------------------------------!

module module_tools

use module_types
use module_constants
use module_math, only: randstr

implicit none

interface unique_elements
  module procedure unique_elements_i
  module procedure unique_elements_r
end interface unique_elements

contains
! ...
! =====================================================================
! =====================================================================
! ...
  subroutine cat(filename)

    character(len=*), intent(in)   :: filename

    ! ... Local variables
    ! ...
    integer unit, ios
    character(len=2048) :: line

    ! ... Use a free unit number
    inquire (iolength=unit) line
    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)

    if (ios /= 0) then
        write(*,*) "Error opening file: ", filename
        return
    end if

    ! ... Read and print line by line
    do
        read(unit,'(A)',iostat=ios) line
        if (ios /= 0) exit
        print '(A)', trim(line)
    end do

    close(unit)
  end subroutine cat
  ! ...
  ! ===================================================================
  ! ...
  function compress(A) result(t)
  ! ... Removes all double whites and spaces before a comma or dot

    character(len=*), intent(in)   :: A
    character(len=len(A))          :: t

    ! ... Local variables
    ! ...
    integer i,n

    ! ... Remove leading whitespaces
    ! ...
    t = adjustl(A)
    n = len_trim(t)

    ! ... Replace tabs by whitespaces
    ! ...
    do i=1,n
      if (iachar(t(i:i)).eq.9) t(i:i) = ' '
    enddo

    ! ... Remove all double whites
    ! ...
    10 i = index(t,'  ')
       if ((i.gt.0).and.(i.lt.n)) then
         t(i+1:) = t(i+2:)
         n = n - 1
         goto 10
       endif

    ! ... Remove space-commas
    ! ...
    20 i = index(t,' ,')
       if ((i.gt.0).and.(i.lt.n)) then
         t(i:) = t(i+1:)
         n = n - 1
         goto 20
       endif

    ! ... Remove space-dots
    ! ...
    30 i = index(t,' .')
       if ((i.gt.0).and.(i.lt.n)) then
         t(i:) = t(i+1:)
         n = n - 1
         goto 30
       endif

  end function compress
  ! ...
  ! ===================================================================
  ! ...
  subroutine crash(message)

    character(len=*)                               :: message

    write(0,*) '> ERROR: ', trim(message)
    stop 1

  end subroutine crash
  ! ...
  ! ===================================================================
  ! ...
  subroutine filename_split (name,path,base,type)

  ! ... Evaluates the path of a given filename. If local folder path is null.
  ! ...
    character*(*), intent(in)                    :: name
    character*(*), intent(out)                   :: path
    character*(*), intent(out)                   :: base
    character*(*), intent(out)                   :: type

    integer io,if
    logical cdot

    io = 1
 01 if = index(name(io:),'/')
    if (if.gt.0) then
       io = io + if
       goto 01
    endif

    path = name(1:io-1)
    base = name(io:)

    cdot = .False.
    io = 1
 02 if = index(base(io:),'.')
    if (if.GT.0) then
       cdot = .True.
       io = io + if
       goto 02
    endif

    if (cdot) then
      type = base(io:)
      base = base(1:io-2)
    else 
      type = ''
    endif

    return
  end subroutine filename_split
  ! ...
  ! ===================================================================
  ! ...
  subroutine get_commandline (commandline)

    character(len=*), intent(out)   :: commandline

    ! ... Local variables
    ! ...
    integer i
    character(len=280) word

    call getarg(0,word)
    commandline = trim(word)
    do i=1,iargc()
      word = ''
      call getarg(i,word)
      commandline = trim(commandline)//' '//trim(adjustl(word))
    enddo

  end subroutine get_commandline
  ! ...
  ! ===================================================================
  ! ...
  logical function is_numeric(string)
  ! ... From Rosetta Code: https://rosettacode.org
  ! ...
    character(len=*), intent(in)                 :: string

    ! ... Local variables
    ! ...
    real(dp) x
    integer e

    read(string,*,iostat=e) x
    is_numeric = (e == 0)

  end function is_numeric
  ! ...
  ! ===================================================================
  ! ...
  function line_replace(A,pattern1,pattern2,ntimes) result(t)
  ! ... Takes pattern1 in A and replaces it by pattern2. Does it ntimes
  ! ... If ntimes < 0 it does it for all appearences

    character(len=*), intent(in)   :: A
    character(len=*), intent(in)   :: pattern1
    character(len=*), intent(in)   :: pattern2
    integer, intent(in), optional  :: ntimes
    character(len=len(A))          :: t

    ! ... Local variables
    ! ...
    integer j,n,n1,n2,nc,ncm

    n = len_trim(A)
    n1 = len(pattern1)
    n2 = len(pattern2)
    if (.not.present(ntimes)) then
      ncm = -1
    else
      ncm = ntimes
    endif

    if (n.eq.0) return
    if (n1.eq.0) return

    t = A
    nc = 0
    10  j = index (t(:n),pattern1)
        nc = nc + 1
        if (j.gt.0) then
          t(j+n2:) = t(j+n1:n)
          t(j:j+n2-1) = pattern2(:n2)
          n = len_trim(t)
          if ((ncm.lt.0).or.(nc.lt.ncm)) goto 10
        endif

  end function line_replace
  ! ...
  ! ===================================================================
  ! ...
  subroutine line_word(line,nw,word)
  ! ... Returns the nw-th word in line

    character(len=*), intent(in)   :: line
    integer, intent(in)            :: nw
    character(len=*)               :: word

    ! ... Local variables
    ! ...
    integer i,j,n,nm,nn,jmax
    character ai
    character(len=len(line))       :: t

    jmax = len(word)
    word = ''

    t    = compress(line)
    nm   = numwords(t)
    if (nw.GT.nm) return

    n    = len_trim(t)
    j = 0
    nn = 0
    do i=1,n-1
      ai = t(i:i)
      if ((ai.ne.' ').and.(ai.ne.',')) then
        j = j + 1
        if (j.le.jmax) word(j:j) = ai
      else
        if (len_trim(word).gt.0) then
          nn = nn + 1
          if (nn.eq.nw) return
          j = 0
          word = ''
        endif
      endif
    enddo

    ai = t(n:n)
    if ((ai.ne.' ').and.(ai.ne.',')) then
      j = j + 1
      if (j.le.jmax) word(j:j) = ai
    endif

    return
  end subroutine line_word
  ! ...
  ! ===================================================================
  ! ...
! ...
! =====================================================================
! ...
  integer function locate(x,xo) result(j)
  ! ... Returns the location of the array x(1:n), such that x(j) < xo < x(j+1)

    real(dp), dimension(:), intent(in)   :: x
    real(dp), intent(in)                 :: xo

    ! ... Local variables
    ! ...
    logical slope
    integer n,jl,jm,ju

    n  = size(x)

    jl = 0
    ju = n+1
    slope = (x(n).ge.x(1))

    do
      if (ju-jl.le.1) exit
      jm = (ju+jl)/2
      if (slope.eqv.(xo.ge.x(jm))) then
        jl = jm
      else
        ju = jm
      endif
    enddo

    if (xo.eq.x(1)) then
      j = 1
    else if (xo.eq.x(n)) then
      j = n - 1
    else
      j = jl
    endif

  end function locate
  ! ...
  ! ===================================================================
  ! ...
  function lowercase(A) result(t)
  ! ... Returns string in lowercase

    character(len=*), intent(in)   :: A
    character(len=len(A))          :: t

    ! ... Local variables
    ! ...
    integer i,nal,nau,nzl,ii,l

    t   = A
    nal = ICHAR('A')
    nau = ICHAR('a')
    nzl = ICHAR('Z')
    l   = nau - nal

    do i=1,len(A)
      ii = ichar(A(i:i))
      if ((ii.ge.nal).and.(ii.le.nzl)) then
        ii = ii + l
        t(i:i) = char(ii)
      endif
    enddo

    return
    end function lowercase
  ! ...
  ! ===================================================================
  ! ...
  function ls(dirname) result(list)

    character(len=*), intent(in)                           :: dirname
    character(maxlen), dimension(:), allocatable           :: list

    integer nl,i,iu
    character(len=8) tmpname
    character(len=maxlen) aa

    ! ... Random filename for temporal storage
    ! ... Send the contents of the selected folder to thar filename
    tmpname = '/tmp/'//randstr(8)
    call system('echo ls '//compress(dirname)//' -1 > '//tmpname)
    call system('ls '//compress(dirname)//' -1 > '//tmpname)

    iu = unitfree()
    open(iu,file=tmpname,status='old')
    nl = numlines(iu)

    allocate(list(nl))

    do i=1,nl
      read(iu,'(A)') aa
      list(i) = trim(aa)
    enddo

    close(iu,status='delete')

  end function ls
  ! ===================================================================
  ! ...
  integer function numlines (iu,type)
  ! ... Returns the number of records in an ASCII or in an UNFORMATTED file
  ! ... For an unformatted file, type must start by 'b' or 'B'.

    integer, intent(in)                      :: iu
    character(len=*), intent(in), optional   :: type

    ! ... Local variables:
    ! ...
    logical ascii
    integer ii

    if (.not.present(type)) then
      ascii = .True.
    else
      if ((type(1:1).eq.'b').or.(type(1:1).eq.'B')) then
        ascii = .False.
      else
        ascii = .True.
      endif
    endif

    rewind(iu)

    ii = 0
    if (ascii) then
  10  read(iu,*,end=20)
         ii = ii + 1
         goto 10
  20  continue
    else
  30  read(iu,end=40)
        ii = ii + 1
        goto 30
  40  continue
    endif

    numlines = ii
    rewind(iu)

    return
  end function numlines
  ! ...
  ! ===================================================================
  ! ...
  integer function numwords(A)
  ! ... Counts the number of words in a string

    character(len=*), intent(in)   :: A

    ! ... Local variables
    ! ...
    integer i,n
    character(len=len(A))          :: t
    character                      :: ai,an,ap

    numwords = 0

    t = compress(A)
    n = len_trim(t)
    if (n.eq.0) return

    numwords = 1
    do i=2,n-1
      ai = t(i:i)
      an = t(i+1:i+1)
      ap = t(i-1:i-1)
      if ((ai.eq.',').and.(an.eq.' ')) then
        numwords = numwords + 1
      else if (ai.eq.',') then
        numwords = numwords + 1
      else if ((ai.eq.' ').and.(ap.ne.',')) then
        numwords = numwords + 1
      endif
    enddo

    return
  end function numwords
  ! ...
  ! ===================================================================
  ! ...
  integer function unitfree()
  ! ... Returns a unit not yet assigned

    ! ... Local variables
    ! ...
    LOGICAL flag

    unitfree = 10
    flag = .TRUE.
    do while (flag)
      unitfree = unitfree + 1
      inquire(unitfree,opened=flag)
    enddo
  
    return
  end function unitfree
  ! ...
  ! ===================================================================
  ! ...
  function uppercase(A) result(t)
  ! ... Returns string in uppercase

    character(len=*), intent(in)   :: A
    character(len=len(A))          :: t

    ! ... Local variables
    ! ...
    integer i,nal,nau,nzl,ii,l

    t   = A
    nal = ICHAR('a')
    nau = ICHAR('A')
    nzl = ICHAR('z')
    l   = nau - nal

    do i=1,len(A)
      ii = ichar(A(i:i))
      if ((ii.ge.nal).and.(ii.le.nzl)) then
        ii = ii + l
        t(i:i) = char(ii)
      endif
    enddo

    return
  end function uppercase
  ! ...
  ! ===================================================================
  ! ...
  subroutine say(text,pos)
  ! ... Writes a text in the screen, using as many lines as required
  ! ... without splitting words.
  ! ... The value fo pos is used for indentation.

    character(len=*), intent(in)         :: text
    integer, intent(in), optional        :: pos

    character(len=180)  :: ww
    character(len=7)  :: fmt
    integer nlen,ipos
    integer i,io,il,n

    if (present(pos)) then
      ipos = pos
    else
      ipos = 1
    endif

    write(fmt,'("(T",i2.2,",A)")') ipos

    ww = ''
    nlen = 80 - ipos
    n = len_trim(text)

    io = 1

10  continue
      il = io + nlen
      if (il.ge.n) then
        il = n
      else
        do i=il,io,-1
          if (text(i-1:i-1).ne.' '.and.text(i:i).EQ.' ') exit
        enddo
        il = i - 1
      endif
      write(6,fmt) text(io:il)
      do i=il,n
        if (text(i-1:i-1).EQ.' '.AND.text(i:i).NE.' ') EXIT
      enddo
      io = i
    if (io.LT.n) goto 10
    
    return
  end subroutine say
  ! ...
  ! ===================================================================
  ! ...
  subroutine strcat(s1,s2)
  ! ... Concatenates s2 into s1. A white character is placed in between.

    character(len=*), intent(inout) :: s1
    character(len=*), intent(in)    :: s2

    s1 = trim(s1)//' '//trim(adjustl(s2))

  end subroutine strcat
  ! ...
  ! ===================================================================
  ! ...
  function token_read(line,prompt) result(ans)
  ! ... Routine to read information of the style key=value from a string.

    character(len=*), intent(in)           :: line
    character(len=*), intent(in)           :: prompt
    character(len=maxlen)                  :: ans

    ! ... Local variables
    ! ...
    integer i,j,i1,i2,nl
    character(len=maxlen) lprompt,lans

    ans = ''

    nl = len_trim(line)
    lprompt = adjustl(prompt)

    do i=2,len_trim(line)-1

      if (line(i:i).eq.'=') then
        ! ... get the index i2
        ! ...
        do j=i-1,1,-1
          if (line(j:j).ne.' ') then
            i2 = j
            exit
          endif
        enddo
        ! ... get the index i1
        ! ...
        do j=i2,1,-1
          if ((line(j:j).eq.' ').or.(line(j:j).eq.',')) then
            i1 = j+1
            exit
          endif
          if (j.eq.1) i1 = 1
        enddo

        if (line(i1:i2).eq.trim(lprompt)) then
          ! ... Ok, we got the left hand side of the token.
          ! ... Now, read the right hand side. There are two
          ! ... options: regular entry or a vector (between square
          ! ... brackets).
          ! ... The equal sign is at position "i"
          
          ! ... Get index i1
          ! ... 
          do j=i+1,nl
            if (line(j:j).ne.' ') then
              i1 = j
              exit
            endif
          enddo
          if (line(i1:i1).eq.'[') then
            ! ... We return a vector
            ! ...
            i2 = -1
            do j=i1+1,nl
              if (line(j:j).eq.']') then
                i2 = j
                exit
              endif
            enddo
            if (i2.eq.-1) call crash('TOKEN_READ: Bad formulated vector')
          else
            do j=i1,nl
              if ((line(j:j).eq.' ').or.(line(j:j).eq.',')) then
                i2 = j-1
                exit
              endif
              if (j.eq.nl) i2 = nl
            enddo
          endif

          ! ... This is the answer. Before returning we check if there was a None
          ! ...
          ans = line(i1:i2)
          lans = lowercase(ans)
          if (trim(lans).eq.'none') ans = ''
          return

        endif
      endif
    enddo

  end function token_read
  ! ...
  ! ===================================================================
  ! ...
  function ReadVector(string) result(A)

    character(len=*), intent(in)                 :: string
    real(dp), dimension(:), allocatable          :: A

    ! ... Local variables
    ! ...
    integer i,n
    character(len=len(string)) att
    character(len=maxlen) word

    if (allocated(A)) deallocate(A)

    att = trim(string)
    att = line_replace(att,',',' ')
    att = line_replace(att,'[',' ')
    att = line_replace(att,']',' ')

    n = numwords(att)
    allocate(A(n))

    do i=1,n
      call line_word(att,i,word)
      read(word,*) A(i)
    enddo

  end function ReadVector
  ! ...
  ! ===================================================================
  ! ...
  subroutine vprint(label,A,mode,fmt)

    character(len=*), intent(in)                 :: label
    real(dp), dimension(:), intent(in)           :: A
    character(len=1), intent(in), optional       :: mode
    character(len=*), intent(in), optional       :: fmt 

    ! ... Local variables:
    ! ...
    integer N,nlines,line,col,rest,ii,i1,i2,i3,i4,i5,ll
    character(len=maxlen) sll
    character(len=maxlen) lfmt

    ll = len(label)
    write(sll,*) ll+2
    sll = compress(sll)

    lfmt = '(T2,A,T'//trim(sll)//',5'//trim(fmt)//')'

    N = size(A)
    rest = mod(N,5)
    nlines = (N - rest)/5

    if ((mode.eq.'H').or.(mode.eq.'h')) then
      ii = 1
      do line=1,nlines
        write(*,lfmt) A(ii:ii+4)
        ii = ii + 5
      enddo
      if (rest.gt.0) write(*,lfmt) A(ii:N)
    else
      do line=1,nlines
        i1 = line
        i2 = i1 + nlines + H(rest.ge.1)
        i3 = i2 + nlines + H(rest.ge.2)
        i4 = i3 + nlines + H(rest.ge.3)
        i5 = i4 + nlines + H(rest.ge.4)
        if (line.eq.1) then
          write(*,lfmt) label, A(i1), A(i2), A(i3), A(i4), A(i5)
        else
          write(*,lfmt) ' ', A(i1), A(i2), A(i3), A(i4), A(i5)
        endif
      enddo
      i1 = line
      i2 = i1 + nlines + 1
      i3 = i2 + nlines + 1
      i4 = i3 + nlines + 1
      if (line.eq.1) then
        sll = label
      else
        sll = ' '
      endif
      if (rest.eq.4) then
        write(*,lfmt) sll(1:ll), A(i1), A(i2), A(i3), A(i4)
      else if (rest.eq.3) then
        write(*,lfmt) sll(1:ll), A(i1), A(i2), A(i3)
      else if (rest.eq.2) then
        write(*,lfmt) sll(1:ll), A(i1), A(i2)
      else 
        write(*,lfmt) sll(1:ll), A(i1)
      endif
    endif

    contains
      integer function H(flag)
  
        logical, intent(in)     :: flag

        if (flag) then
          H = 1
        else
          H = 0
        endif

     end function H

  end subroutine vprint 
  ! ...
  ! ===================================================================
  ! ...
  function unique_elements_i(arr) result(unique_arr)

    integer, intent(in)                          :: arr(:)
    integer, allocatable                         :: unique_arr(:)

    ! ... Local variables
    ! ...
    logical, dimension(size(arr)) :: mask
    integer i

    mask(:) = .true.
    do i = 2, size(arr)
        mask(i) = .not. any(arr(i) == arr(1:i-1))
    end do

    unique_arr = pack(arr, mask)
  end function unique_elements_i
  ! ...
  ! ===================================================================
  ! ...
  function unique_elements_r(arr) result(unique_arr)

    real(dp), intent(in)                          :: arr(:)
    real(dp), allocatable                         :: unique_arr(:)

    ! ... Local variables
    ! ...
    logical, dimension(size(arr)) :: mask
    integer i

    mask(:) = .true.
    do i = 2, size(arr)
        mask(i) = .not. any(arr(i) == arr(1:i-1))
    end do

    unique_arr = pack(arr, mask)
  end function unique_elements_r
  ! ...
  ! ===================================================================
  ! ...
end module module_tools
