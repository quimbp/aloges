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
! - codar_ruv_read
! - codar_att_get
! -------------------------------------------------------------------------!

module module_codar

use module_types
use module_tools
use module_time

implicit none

type type_codar_att
  character(len=maxlen)                   :: key = ""
  character(len=maxlen)                   :: value = ""
end type type_codar_att

type type_codar_table
  character(len=20)                                      :: Type
  integer                                                :: Rows
  integer                                                :: Columns
  character(len=maxlen)                                  :: ColumnTypes
  character(len=maxlen)                                  :: Names,Units
  character(len=10), dimension(:), allocatable           :: Name
  character(len=10), dimension(:), allocatable           :: Unit
  real(sp), dimension(:,:), allocatable                  :: Value
end type type_codar_table

type type_codar
  integer                                                :: nlines = 0
  integer                                                :: natts    = 0
  integer                                                :: ntables  = 0
  real(dp)                                               :: xo,yo
  character(len=20)                                      :: TimeStamp = ""
  type(type_codar_att), dimension(:), allocatable        :: Att
  type(type_codar_table), dimension(:), allocatable      :: Table

  contains
    procedure :: ruv_read => codar_ruv_read
    procedure :: att_get => codar_att_get
end type type_codar

contains 
! ...
! ==========================================================================
! ==========================================================================
! ...
subroutine codar_ruv_read(CODAR,filename)

  class(type_codar), intent(inout)         :: CODAR
  character(len=*), intent(in)             :: filename

  ! ... Local variables
  ! ...
  logical intable
  integer i1,i2,i3,i4,jcol
  integer iu,nlines,natts,ntables,iline,iatt,itable
  integer ncols,nrows,iheader,irow
  character(len=400) line
  character(len=40) word
  character(len=maxlen) key,value

  iu = unitfree()
  write(*,*) 'Opening RUV file: ', trim(filename)
  open(iu,file=filename,status='old')

  ! ... Read the number of attributes
  ! ...
  rewind(iu)
  nlines = 0
  ntables = 0
  do
    read(iu,*,end=10) line
    if (line(1:11).eq.'%TableStart') ntables = ntables + 1
    nlines = nlines + 1
  enddo
10 continue

  if (line(1:5).ne.'%End:') stop 'Invalid RUV file'

  CODAR%nlines = nlines 
  CODAR%ntables = ntables

  nlines  = nlines - 1      ! Skip the last line
  natts   = 0
  itable  = 0
  intable = .False.

  rewind(iu)
  do iline=1,nlines
    read(iu,'(A)',err=20) line
    if (.not.intable.and.line(1:6).eq.'%Table') then
      intable = .true.
      itable = itable + 1
    endif
    if (intable.and.line(1:9).eq.'%TableEnd') intable = .false.
   
    if (.not.intable) then
      if (line(1:2).eq.'%%'.or.line(1:9).eq.'%TableEnd') then
        ! not an attribute
      else
        natts = natts + 1
      endif
    endif

  enddo  

  CODAR%natts = natts
  if (allocated(CODAR%Att)) deallocate(CODAR%Att)
  if (allocated(CODAR%Table)) deallocate(CODAR%Table)

  allocate(CODAR%Att(natts))
  allocate (CODAR%Table(CODAR%ntables))

  intable = .false.
  iatt    = 0
  itable  = 0

  rewind(iu)
  do iline=1,nlines
    read(iu,'(A)',err=20) line
    if (.not.intable.and.line(1:6).eq.'%Table') then
      intable = .true.
      itable = itable + 1
      iheader = 0
      irow = 0
    endif
    if (intable.and.line(1:9).eq.'%TableEnd') intable = .false.
   
    ! ... Read attribute
    ! ...
    if (.not.intable) then
      if (line(1:2).eq.'%%'.or.line(1:9).eq.'%TableEnd') then
        ! not an attribute
      else
        iatt = iatt + 1
        i1 = 2
        i2 = index(line,':') - 1
        i3 = i2 + 2
        i4 = len_trim(line)
        CODAR%Att(iatt)%key = ""
        CODAR%Att(iatt)%value = ""
        CODAR%Att(iatt)%key = line(i1:i2)
        CODAR%Att(iatt)%value = trim(line(i3:i4))
        !print*, iatt, i1, i2, i3, i4, '      ', trim(CODAR%Att(iatt)%key), trim(CODAR%Att(iatt)%value)
      endif
    endif

    ! ... Read table
    ! ...
    if (intable) then
      if (line(1:9).eq.'%TableEnd') then
        ! Nothing to be done
      else if (line(1:10).eq.'%TableType') then
        i3 = index(line,':') + 1
        i4 = len_trim(line)
        CODAR%Table(itable)%Type = line(i3:i4)
      else if (line(1:13).eq.'%TableColumns') then
        i3 = index(line,':') + 1
        i4 = len_trim(line)
        read(line(i3:i4),*) ncols
        CODAR%Table(itable)%Columns = ncols
        allocate(CODAR%Table(itable)%Name(ncols))
        allocate(CODAR%Table(itable)%Unit(ncols))
      else if (line(1:10).eq.'%TableRows') then
        i3 = index(line,':') + 1
        i4 = len_trim(line)
        read(line(i3:i4),*) nrows
        CODAR%Table(itable)%Rows = nrows
        allocate(CODAR%Table(itable)%Value(nrows,ncols))
      else if (line(1:12).eq.'%ColumnTypes') then
        i3 = index(line,':') + 1
        i4 = len_trim(line)
        CODAR%Table(itable)%ColumnTypes = line(i3:i4)
      else if (line(1:2).eq.'%%') then
        iheader = iheader + 1
        if (iheader.eq.1) then
          CODAR%Table(itable)%names = trim(line(3:))
        else if (iheader.eq.2) then
          CODAR%Table(itable)%units = trim(line(3:))
        endif
      else if (line(1:1).ne.'%') then
        irow = irow + 1
        do jcol=1,ncols
          call line_word(line,jcol,word)
          read(word,*) CODAR%Table(itable)%Value(irow,jcol)
        enddo
      endif
    endif
 
  enddo

  ! ... Antenna origin
  ! ...
  key = 'Origin' 
  call codar_att_get(CODAR,key,value)
  call line_word(value,1,word)
  read(word,*) CODAR%yo
  call line_word(value,2,word)
  read(word,*) CODAR%xo

  ! ... Antenna TimeStamp
  ! ...
  key = 'TimeStamp' 
  call codar_att_get(CODAR,key,value)
  CODAR%TimeStamp = trim(value)

  return
  

20 stop 'error' 

end subroutine codar_ruv_read
! ...
! ==========================================================================
! ...
subroutine codar_att_get(CODAR,key,value)

  class(type_codar), intent(inout)         :: CODAR
  character(len=maxlen), intent(inout)     :: key
  character(len=maxlen), intent(out)       :: value

  ! ... Local variables
  ! ...
  integer i

  value = ""
  do i=1,CODAR%natts
    if (trim(CODAR%Att(i)%key).eq.trim(key)) then
      value = trim(CODAR%Att(i)%value)
      return
    endif
  enddo

end subroutine codar_att_get
! ...
! ==========================================================================
! ...
real(dp) function compass_to_polar(compass)
  real(dp), intent(in) :: compass
  compass_to_polar = mod(450.0 - compass, 360.0)
end function compass_to_polar

real(dp) function polar_to_compass(polar)
  real(dp), intent(in) :: polar
  polar_to_compass = mod(450.0 - polar, 360.0)
end function polar_to_compass
! ...
! ==========================================================================
! ...
subroutine combine_radials_wls(vr, theta, w, u, v, cov_uv, gdop, angle_ref)
! ...
! ... Weighted LS to combine radials
! ...
  real(dp), intent(in)            :: vr(:), theta(:), w(:)
  real(dp), intent(out)           :: u, v, cov_uv(2,2), gdop
  character(len=*), optional      :: angle_ref

  ! ... Local variables
  ! ...
  logical polar
  integer  :: n,i
  real(dp) :: A(2,2), b(2), ai, wi, c, s

  n = size(vr)

  ! ... By default, we consider Polar angle reference
  ! ...
  polar = .true.
  if (present(angle_ref)) then
    if ((angle_ref(1:1).eq.'b').or.(angle_ref(1:1).eq.'B')) polar = .false.
  else
    write(*,*) 'WARNING: Reference angle system not provided. Using polar'
  endif
  !write(*,*) 'Polar reference: ', polar
  
  A(:,:) = 0.0; b(:) = 0.0
  do i=1,n
    c = cos(theta(i)); s = sin(theta(i))
    wi = w(i)
    if (polar) then
      A(1,1) = A(1,1) + wi*c*c
      A(1,2) = A(1,2) + wi*c*s
      A(2,2) = A(2,2) + wi*s*s
      b(1) = b(1) + wi*c*vr(i)
      b(2) = b(2) + wi*s*vr(i)
    else
      A(1,1) = A(1,1) + wi*s*s
      A(1,2) = A(1,2) + wi*c*s
      A(2,2) = A(2,2) + wi*c*c
      b(1) = b(1) + wi*s*vr(i)
      b(2) = b(2) + wi*c*vr(i)
    endif
  end do
  A(2,1) = A(1,2)

  call solve22(A, b, u, v, cov_uv)
  gdop = sqrt(cov_uv(1,1) + cov_uv(2,2))

  contains

    subroutine solve22(A, b, x1, x2, cov)
      real(dp), intent(in):: A(2,2), b(2)
      real(dp), intent(out):: x1, x2, cov(2,2)

      ! ... Local variables
      ! ...
      real(dp) :: det

      det = A(1,1)*A(2,2) - A(1,2)*A(2,1)  ! Determinant
      if (abs(det) .le. 1e-10) then
        x1 = 0.; x2 = 0.; cov = huge(1.0_dp)
        return
      end if
      cov(1,1) =  A(2,2)/det
      cov(2,2) =  A(1,1)/det
      cov(1,2) = -A(1,2)/det
      cov(2,1) =  cov(1,2)
      x1 = cov(1,1)*b(1) + cov(1,2)*b(2)
      x2 = cov(2,1)*b(1) + cov(2,2)*b(2)
    end subroutine solve22
  end subroutine combine_radials_wls
! ...
! ==========================================================================
! ...
end module module_codar
