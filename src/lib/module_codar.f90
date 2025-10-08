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

use module_tools

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
  character(len=4), dimension(40)                        :: Heading
  character(len=10), dimension(:), allocatable           :: Name
  character(len=10), dimension(:), allocatable           :: Unit
  real(sp), dimension(:,:), allocatable                  :: Value
  contains
      procedure :: id  => codar_table_typeid
end type type_codar_table

! ... Codar RUV structure Table 1 (Spec LLUVSpec 1.51):
! ...
! ... Columns
! ...  1 - Longitude  (deg)
! ...  2 - Latitude   (deg)
! ...  3 - U comp     (cm/s)
! ...  4 - V comp     (cm/s)
! ...  5 - Vectorflag (GridCode)
! ...  6 - Spatial quality
! ...  7 - Temporal quality
! ...  8 - Pattern distance
! ...  9 - SNR        (dB)
! ... 10 - Velocity Maximum
! ... 11 - Velocity Minimum
! ... 12 - Spatial Count
! ... 13 - Temporal Count
! ... 14 - X Distance (km)
! ... 15 - Y Distance (km)
! ... 16 - Range      (km)
! ... 17 - Bearing    (True)
! ... 18 - Velocity   (cm/s) 
! ... 19 - Direction  (True)
! ... 20 - Spectra    (RngCell)
! ... 21 - Q201       (flag)
! ...  .    
! ...  .    Quality flags
! ...  .
type type_codar
  integer                                                :: nlines = 0
  integer                                                :: natts    = 0
  integer                                                :: ntables  = 0
  integer                                                :: Year,Month,Day
  integer                                                :: Hour,Minute,Second
  integer                                                :: RangeStart
  integer                                                :: RangeEnd
  real(dp)                                               :: TimeCoverage
  real(dp)                                               :: xo,yo
  real(dp)                                               :: RangeResolution
  real(dp)                                               :: AngularResolution
  character(len=4)                                       :: Site = ""
  character(len=20)                                      :: TimeStamp = ""
  character(len=40)                                      :: TimeZone = ""
  type(type_codar_att), dimension(:), allocatable        :: Att
  type(type_codar_table), dimension(:), allocatable      :: Table

  contains
    procedure :: ruv_read  => codar_ruv_read
    procedure :: att_get   => codar_att_get
    procedure :: ruv_write => codar_ruv_write
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
  integer i,i1,i2,i3,i4,jcol
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
      else if (line(1:17).eq.'%TableColumnTypes') then
        i3 = index(line,':') + 1
        i4 = len_trim(line)
        CODAR%Table(itable)%ColumnTypes = line(i3:i4)
        do i=1,CODAR%Table(1)%Columns
          word = ""
          call line_word(CODAR%Table(1)%ColumnTypes,i,word)
          CODAR%Table(itable)%Heading(i) = trim(word)
        enddo
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
  call line_word(value,1,word)
  read(word,*) CODAR%Year
  call line_word(value,2,word)
  read(word,*) CODAR%Month
  call line_word(value,3,word)
  read(word,*) CODAR%Day
  call line_word(value,4,word)
  read(word,*) CODAR%Hour
  call line_word(value,5,word)
  read(word,*) CODAR%Minute
  call line_word(value,6,word)
  read(word,*) CODAR%Second

  ! ... Antenna TimeZone
  ! ...
  key = 'TimeZone' 
  call codar_att_get(CODAR,key,value)
  CODAR%TimeZone = trim(value)

  ! ... Time coverage
  ! ...
  key = 'TimeCoverage' 
  call codar_att_get(CODAR,key,value)
  call line_word(value,1,word)
  read(word,*) CODAR%TimeCoverage

  ! ... Antenna Site
  ! ...
  key = 'Site' 
  call codar_att_get(CODAR,key,value)
  call line_word(value,1,word)
  CODAR%Site = trim(word)

  ! ... RangeStart
  ! ...
  key = 'RangeStart' 
  call codar_att_get(CODAR,key,value)
  call line_word(value,1,word)
  read(word,*) CODAR%RangeStart

  ! ... RangeEnd
  ! ...
  key = 'RangeEnd' 
  call codar_att_get(CODAR,key,value)
  call line_word(value,1,word)
  read(word,*) CODAR%RangeEnd

  ! ... RangeResolution
  ! ...
  key = 'RangeResolutionKMeters'
  call codar_att_get(CODAR,key,value)
  call line_word(value,1,word)
  read(word,*) CODAR%RangeResolution
  CODAR%RangeResolution = CODAR%RangeResolution * 1000.0D0 ! meters

  ! ... AngularResolution
  ! ...
  key = 'AngularResolution' 
  call codar_att_get(CODAR,key,value)
  call line_word(value,1,word)
  read(word,*) CODAR%AngularResolution

  ! ... Get the Ids of the variables:
  ! ... THis information is stored in theColumnTypes of Table 1.
  ! ... LOND LATD VELU VELV VFLG ESPC ETMP EDTP EASN MAXV MINV ERSC ERTC XDST YDST RNGE BEAR VELO HEAD SPRC

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
  real(dp) :: A(2,2), b(2), wi, c, s

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
subroutine codar_ruv_write(CODAR,filename)

  class(type_codar), intent(inout)         :: CODAR
  character(len=*), intent(in)             :: filename

  ! ... Local variables
  ! ...
  integer iu

  iu = unitfree()
  write(*,*) 'Creating RUV file: ', trim(filename)
  open(iu,file=filename,status='unknown')

!      1 %CTF: 1.00
!      2 %FileType: LLUV rdls "RadialMap"
!      3 %LLUVSpec: 1.51  2021 07 01
!      4 %UUID: A23C9091-4F03-4447-BAFA-26D9C32CE744
!      5 %Manufacturer: CODAR Ocean Sensors. SeaSonde
!      6 %Site: AREN ""
!      7 %TimeStamp: 2025 02 05  10 00 00
!      8 %TimeZone: "UTC" +0.000 0 "Atlantic/Reykjavik"
!      9 %TimeCoverage: 75.000 Minutes
!     10 %Origin: 41.5775833    2.5577333
!     11 %GreatCircle: "WGS84" 6378137.000  298.257223562997
!     12 %GeodVersion: "CGEO" 2.00  2019 03 05

  write(iu,'(T1,A)') '%CTF: 1.00'
  write(iu,'(T1,A)') '%FileType: LLUV rdls "RadialMap"'
  write(iu,'(T1,A)') '%LLUVSpec: 1.51  2021 07 01'
  write(iu,'(T1,A)') '%UUID: '
  write(iu,'(T1,A)') '%Manufacturer: CODAR Ocean Sensors. SeaSonde'
  write(iu,'(T1,A)') '%Site: ' // trim(CODAR%Site)//' ""'

  close(iu)

end subroutine codar_ruv_write
! ...
! ==========================================================================
! ...
function codar_table_typeid(TABLE,name) result(id)

  class(type_codar_table), intent(in)          :: TABLE
  character(len=4), intent(in)                 :: name
  integer                                      :: id

  ! ... Local variables
  ! ... 
  integer i

  id = -1
  do i=1,TABLE%Columns
    if (TABLE%Heading(i).eq.name) then
      id = i
      return
    endif
  enddo
  stop 'ERROR in CODAR_TABLE_TYPEID: Column name not found'
 
  end function codar_table_typeid 
! ...
! ==========================================================================
! ...
end module module_codar
