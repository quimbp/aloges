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
! - codar_ruv_read                                                         !
! - codar_att_get                                                          !
! -------------------------------------------------------------------------!

module module_codar

use module_types, only: dp, maxlen
use module_constants, only: deg2rad
use module_tools, only: crash
implicit none (type, external)

type type_codar_table
  character(len=20)                                :: Type = ""
  integer                                          :: nrows = 0
  integer                                          :: ncols = 0
  real(dp)                                         :: missing = -999.0_dp
  character(len=4), allocatable                    :: ColumnNames(:)
  character(len=maxlen)                            :: Names,Units
  real(dp), allocatable                            :: data(:,:) ! (nrows, ncols)
  contains
    procedure                                      :: id  => codar_table_typeid
    procedure                                      :: get => codar_table_get
end type type_codar_table

type type_codar
  character(len=maxlen)                             :: filename = ""
  character(len=4)                                  :: Site = ""
  character(len=30)                                 :: CTF = ""
  character(len=30)                                 :: FileType = ""
  character(len=30)                                 :: LLUVSpec = ""
  integer                                           :: nlines   = 0
  integer                                           :: natts    = 0
  integer                                           :: ntables  = 0
  integer                                           :: Year,Month,Day
  integer                                           :: Hour,Minute,Second
  ! --- Location ---
  real(dp)                                          :: lon0
  real(dp)                                          :: lat0
  ! --- Geometry ---
  ! ... For CODAR-style HFR systems:
  ! ... Bearing is a geometric property of the radar cell
  ! ... It is the look direction of the radar beam
  ! ... It is fixed by antenna geometry and does not depend on the flow
  ! ... Formally, the bearing is the azimuth of the outward radial unit vector
  ! ... (from antenna → ocean cell).
  ! ... bearing is not negotiable and should never be changed per observation.
  ! ...
  integer                                           :: nbearings   ! number of bearings      
  integer                                           :: nranges  ! number of range cells
  integer                                           :: rmin    ! RangeStart
  integer                                           :: rmax    ! RangeEnd
  real(dp)                                          :: bearing_res
  real(dp)                                          :: range_res    ! km

  real(dp), allocatable                             :: bearing(:)   ! [nbearings] degrees
  real(dp), allocatable                             :: range(:)     ! [nranges] meters

  ! --- Radar grid ---
  real(dp), allocatable                             :: cell_lon(:,:)  ! [nbearings,nranges]
  real(dp), allocatable                             :: cell_lat(:,:)  ! [nbearings,nranges]
  integer, allocatable                              :: cell_mask(:,:) ! 1=water,0=land

  ! --- Radial measurements ---
  real(dp), allocatable                             :: cell_vel(:,:)  ! [nbearings,nranges] (m)

  ! --- Measurement characteristics ---
  real(dp), allocatable                             :: sigma_vr(:,:) ! radial velocity std dev
  real(dp)                                          :: freq_hz       ! operating frequency
  real(dp)                                          :: max_range     ! max usable range (m)

  ! --- Precomputed projection ---
  real(dp), allocatable                             :: cos_theta(:,:) ! cos(bearing)
  real(dp), allocatable                             :: sin_theta(:,:) ! sin(bearing)


  real(dp)                                          :: TimeCoverage
  real(dp)                                          :: RangeResolution
  real(dp)                                          :: AngularResolution
  character(len=20)                                 :: TimeStamp = ""
  character(len=40)                                 :: TimeZone = ""
  character(len=maxlen),  allocatable               :: line(:)
  type(type_codar_table), dimension(:), allocatable :: Table

  contains
    procedure :: ruv_read  => codar_ruv_read
    procedure :: att_get   => codar_att_get
    procedure :: row2grid  => codar_row2grid
    procedure :: ruv_write => codar_ruv_write
end type type_codar

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
contains 
! ...
! ==========================================================================
! ==========================================================================
! ...
  subroutine codar_ruv_read(CODAR,filename,verbose_level,status)
    ! ... The radial velocity, vr, is: vr = u ⋅ r^
    ! ... where: u=(u,v) is the horizontal current
    ! ... r^ is the unit vector defined by the bearing
    ! ...
    ! ... Consequences:
    ! ...  vr > 0: flow away from the antenna
    ! ...  vr < 0: flow toward the antenna
    ! ...
    ! ... This sign convention is the one used by:
    ! ...    CODAR
    ! ...    WERA
    ! ...    Most published HFR DA and inversion work
    ! ...
    ! ... Speed + direction in RUV files
    ! ... Some RUV products report:
    ! ...  * a speed (nominally non-negative)
    ! ...  * a direction
    ! ...
    ! ... This direction is usually:
    ! ...  * the direction of flow, not the radar bearing
    ! ...    given in oceanographic convention (direction
    ! ...    toward which the water moves)
    ! ... This direction can indeed differ from the radar bearing by ~180°.
    ! ...
    ! ... ✅ Correct and usual practice
    ! ...  * Keep the bearing fixed
    ! ...  * Convert speed + direction into a signed radial velocity
    ! ...  * Store only the signed radial velocity
    ! ... That is:
    ! ...                vr = ∣u∣ cos(θ_flow−θ_bearing)
    ! ... Once converted:
    ! ... You never touch the bearing
    ! ... The sign of vr carries all the directional information
    ! ... This is what most operational HFR toolchains do internally, even if the 
    ! ... files expose speed/direction.
    ! ...

    class(type_codar), intent(inout)         :: CODAR
    character(len=*), intent(in)             :: filename
    integer, intent(in), optional            :: verbose_level
    integer, intent(out), optional           :: status

    ! ... Local variables
    ! ...
    logical in_table
    integer verbosity
    integer i,rmin,rmax,nbearings,nranges
    integer iu,ios,nlines,natts,ntables,iline,iatt,itable
    integer ncols,nrows,iheader,irow,id_bear
    real(dp) freq_hz, bearing_res, range_res
    real(dp) xbear, min_bearing, closest_theoretical, offset
    character(len=400) line

    verbosity = 1
    if (present(verbose_level)) verbosity = verbose_level
    if (present(status)) status = 0

    if (verbosity.ge.1) write(*,*) 'Opening RUV file: ', trim(filename)
    open(newunit=iu,file=filename,status='old', iostat=ios)
    if (ios.ne.0) then
      if (present(status)) then
        status = -1
        return
      else
        call crash('Cannot open ruv file')
      endif
    endif

    ! ... Read the number of lines
    ! ...
    nlines = 0
    ntables = 0
    rewind(iu)
    do
      read(iu,'(A)',iostat=ios) line
      if (ios.ne.0) exit             ! End-of-file
      if (index(line,'%CTF') > 0) then
        read(line(index(line,':')+1:),'(A)') CODAR%CTF
      else if (index(line,'%FileType') > 0) then
        read(line(index(line,':')+1:),'(A)') CODAR%FileType
      else if (index(line,'%LLUVSpec') > 0) then
        read(line(index(line,':')+1:),'(A)') CODAR%LLUVSpec
      endif
      if (line(1:11).eq.'%TableStart') ntables = ntables + 1
      nlines = nlines + 1
    enddo

    ! ... Check data format and specs
    if (verbosity.ge.2) then
      write(*,*) 'CTF: ', CODAR%CTF
      write(*,*) 'FileType: ', CODAR%FileType
      write(*,*) 'LLUVSpec: ', CODAR%LLUVSpec
    endif
    if (len_trim(CODAR%FileType).eq.0) call crash('Invalid RUV file')
    if (len_trim(CODAR%LLUVSpec).eq.0) call crash('Invalid RUV file')

    ! ... Last line
    ! ...
    if (line(1:5).ne.'%End:') stop 'Invalid RUV file'

    CODAR%nlines  = nlines 
    CODAR%ntables = ntables
    if (verbosity.ge.3) then
      write(*,*) ' > CODAR%nlines: ', CODAR%nlines
      write(*,*) ' > CODAR%ntables: ', CODAR%ntables
    endif

    if (allocated(CODAR%line)) deallocate(CODAR%line)
    if (allocated(CODAR%Table)) deallocate(CODAR%Table)
    allocate (CODAR%line(CODAR%nlines))
    allocate (CODAR%Table(CODAR%ntables))

    natts   = 0

    rewind(iu)
    do iline=1,nlines
      read(iu,'(A)') line
      CODAR%line(iline) = trim(line)

      if (line(1:2).eq.'%%'.or.index(line,'%Table').gt.0) cycle

      natts = natts + 1

      if (index(line,'%Site') > 0) then
        read(line(index(line,':')+1:),*) CODAR%site
      else if (index(line,'%Origin') > 0) then
        read(line(index(line,':')+1:),*) CODAR%lat0, CODAR%lon0 
      else if (index(line,'%TimeStamp') > 0) then
        read(line(index(line,':')+1:),*) CODAR%year, CODAR%month, CODAR%day, &
                                         CODAR%hour, CODAR%minute, CODAR%Second
      else if (index(line,'%TimeZone') > 0) then
        read(line(index(line,':')+1:),*) CODAR%TimeZone
      else if (index(line,'%TimeCoverage') > 0) then
        read(line(index(line,':')+1:),*) CODAR%TimeCoverage
      else if (index(line,'%AngularResolution') > 0) then
        read(line(index(line,':')+1:),*) bearing_res
      else if (index(line,'%RangeResolution') > 0) then
        read(line(index(line,':')+1:),*) range_res
      else if (index(line,'%RangeStart') > 0) then
        read(line(index(line,':')+1:),*) rmin
      else if (index(line,'%RangeEnd') > 0) then
        read(line(index(line,':')+1:),*) rmax
      else if (index(line,'%TransmitCenterFreqMHz') > 0) then
        read(line(index(line,':')+1:),*) freq_hz
      end if

    end do
    close(iu)

    CODAR%natts = natts

    ! --- Geometry ---
    nranges = rmax
    nbearings  = int(360.0_dp / bearing_res)

    CODAR%filename    = trim(filename)
    CODAR%rmin        = rmin
    CODAR%rmax        = rmax
    CODAR%nranges     = nranges
    CODAR%nbearings   = nbearings
    CODAR%freq_hz     = freq_hz
    CODAR%bearing_res = bearing_res
    CODAR%range_res   = range_res  
    CODAR%max_range   = (rmax-0.5_dp)*range_res    ! km

    if (allocated(CODAR%bearing)) deallocate(CODAR%bearing)
    if (allocated(CODAR%range)) deallocate(CODAR%range)
    if (allocated(CODAR%cell_lon)) deallocate(CODAR%cell_lon)
    if (allocated(CODAR%cell_lat)) deallocate(CODAR%cell_lat)
    if (allocated(CODAR%cell_vel)) deallocate(CODAR%cell_vel)
    if (allocated(CODAR%cell_mask)) deallocate(CODAR%cell_mask)
    if (allocated(CODAR%sigma_vr)) deallocate(CODAR%sigma_vr)
    if (allocated(CODAR%cos_theta)) deallocate(CODAR%cos_theta)
    if (allocated(CODAR%sin_theta)) deallocate(CODAR%sin_theta)

    allocate(CODAR%bearing(nbearings))
    allocate(CODAR%range(nranges))
    allocate(CODAR%cell_lon(nbearings,nranges))
    allocate(CODAR%cell_lat(nbearings,nranges))
    allocate(CODAR%cell_vel(nbearings,nranges))
    allocate(CODAR%cell_mask(nbearings,nranges))
    allocate(CODAR%sigma_vr(nbearings,nranges))
    allocate(CODAR%cos_theta(nbearings,nranges))
    allocate(CODAR%sin_theta(nbearings,nranges))
    CODAR%cell_mask = 0
    CODAR%sigma_vr  = 0.0_dp

    ! Bearings in degrees : First cell, by default, is 0.0. 
    ! Latter, we will add the antenna's starting point.
    CODAR%bearing = [( (i-1)*bearing_res, i=1,nbearings )]

    ! Range centers (meters)
    CODAR%range = [ (real(i,dp)*range_res*1000.0_dp, i=1,nranges) ]         ! meters

    in_table = .false.
    iatt    = 0
    itable  = 0

    do iline=1,nlines
      line = CODAR%line(iline)
      if (.not.in_table.and.line(1:6).eq.'%Table') then
        in_table = .true.
        itable = itable + 1
        iheader = 0
        irow = 0
      endif
      if (in_table.and.line(1:9).eq.'%TableEnd') in_table = .false.
   
      ! ... Read table
      ! ...
      if (in_table) then


        if (line(1:9).eq.'%TableEnd') then
          ! Nothing to be done
        else if (index(line,'%TableType') > 0) then
          read(line(index(line,':')+1:),*) CODAR%Table(itable)%Type
          CODAR%Table(itable)%ncols = 0
          CODAR%Table(itable)%nrows = 0
        else if (index(line,'%TableColumns') > 0) then
          read(line(index(line,':')+1:),*,iostat=ios) ncols
          if (ios.eq.0) CODAR%Table(itable)%ncols = ncols
          if (verbosity.ge.3) write(*,*) ' > CODAR%Table [itable, ncols]: ', itable, ncols
        else if (index(line,'%TableRows') > 0) then
          read(line(index(line,':')+1:),*,iostat=ios) nrows
          if (ios.eq.0) CODAR%Table(itable)%nrows = nrows
          if (verbosity.ge.3) write(*,*) ' > CODAR%Table [itable, nrows]: ', itable, nrows
          allocate(CODAR%Table(itable)%data(nrows,ncols))
        else if (line(1:17).eq.'%TableColumnTypes') then
          if (allocated(CODAR%Table(itable)%ColumnNames)) &
                    deallocate(CODAR%Table(itable)%ColumnNames)
          allocate(CODAR%Table(itable)%ColumnNames(ncols))
          read(line(index(line,':')+1:),*) (CODAR%Table(itable)%ColumnNames(i), i=1,ncols)
        else if (line(1:2).eq.'%%') then
          iheader = iheader + 1
          if (iheader.eq.1) then
            CODAR%Table(itable)%names = trim(line(3:))
          else if (iheader.eq.2) then
            CODAR%Table(itable)%units = trim(line(3:))
          endif
        else if (line(1:1).ne.'%') then
          irow = irow + 1
          read(line,*) CODAR%Table(itable)%data(irow,:)
        endif

      endif
 
    enddo

    ! ... Check that all tables are full:
    ! ...
    !do itable=1,CODAR%ntables
    !  print*, itable, CODAR%Table(itable)%nrows, CODAR%Table(itable)%ncols
    !enddo
   
    ! ... Get the ID of the table with the LLUV data (usually the first one).
    ! ... 
    itable = find_table_by_type(CODAR, 'LLUV')
    if (itable.le.0) then
      if (present(status)) then
        status = -1
        return
      else
        call crash('Cannot find LLUV table')
      endif
    endif

    ! ... Check that the table is full
    ! ...
    !print*, itable, CODAR%Table(itable)%nrows, CODAR%Table(itable)%ncols
    if (CODAR%Table(itable)%nrows.le.0.or.CODAR%Table(itable)%ncols.le.0) then
      if (present(status)) then
        status = -1
        return
      else
        call crash('Cannot open ruv file')
      endif
    endif

    ! ... Check units (km) in table header
    ! ...
    if (index(CODAR%Table(itable)%units, 'km') == 0) then
        write(*,*) 'WARNING: Range units may not be km'
    endif 

    id_bear  = CODAR%Table(itable)%id('BEAR')  ! Bearing (True degrees)
    if (id_bear.le.0) then
      if (present(status)) then
        status = -1
        return
      else
        call crash('Cannot find BEAR column')
      endif
    endif

    min_bearing = minval(CODAR%Table(itable)%data(:,id_bear))
    closest_theoretical = floor(min_bearing / bearing_res) * bearing_res 
    offset = min_bearing - closest_theoretical
    if (verbosity.ge.2) write(*,*) ' > Bearing offset: ', offset
    CODAR%bearing(:) = offset + CODAR%bearing(:)
  
    if (any(CODAR%bearing.gt.360)) stop 'Bearing > 360'

    ! ... Populate grid arrays from table data
    call codar_populate_grid(CODAR,verbosity)


  end subroutine codar_ruv_read
! ...
! ==========================================================================
! ...
  function codar_att_get(CODAR,key) result(value)

    class(type_codar), intent(in)            :: CODAR
    character(len=*), intent(in)             :: key
    character(len=:), allocatable            :: value

    ! ... Local variables
    ! ...
    integer i,ilen
    character(len=maxlen) line, answer

    do i=1,CODAR%nlines
      line = CODAR%line(i)
      if (index(line,key) > 0) then
        read(line(index(line,':')+1:),'(A)') answer
        ilen = len_trim(answer)
        if (allocated(value)) deallocate(value)
        allocate(character(len=ilen) :: value)
        value = trim(answer)
        return
      endif
    enddo

  end function codar_att_get
  ! ...
  ! ==========================================================================
  ! ...
  subroutine codar_ruv_write(CODAR,filename)

    class(type_codar), intent(inout)         :: CODAR
    character(len=*), intent(in)             :: filename

    ! ... Local variables
    ! ...
    integer iu

    write(*,*) 'Creating RUV file: ', trim(filename)
    open(newunit=iu,file=filename,status='unknown')

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
    do i=1,TABLE%ncols
      if (TABLE%ColumnNames(i).eq.name) then
        id = i
        return
      endif
    enddo
    ! Returns -1 if not found !
    !stop 'ERROR in CODAR_TABLE_TYPEID: Column name not found'
 
    end function codar_table_typeid 
  ! ...
  ! ==========================================================================
  ! ...
  function codar_table_get(TABLE,name) result(A)

    class(type_codar_table), intent(in)          :: TABLE
    character(len=4), intent(in)                 :: name
    real(dp), dimension(:), allocatable          :: A

    ! ... Local variables
    ! ...
    integer id

    if (TABLE%nrows.le.0) return
    if (allocated(A)) deallocate(A)

    allocate(A(TABLE%nrows))
    id = TABLE%id(name)
    A(:) = TABLE%data(:,id)

  end function codar_table_get
  ! ...
  ! ==========================================================================
  ! ...
  ! ========================================================================
  ! Subroutine to populate CODAR grid arrays from table data
  ! ========================================================================
  ! 
  ! This routine extracts lon, lat, and velocity data from the CODAR
  ! data table and populates the cell_lon, cell_lat, and cell_mask arrays.
  !
  ! Key points:
  ! - The table contains irregularly gridded data (not all bearing/range
  !   combinations have observations)
  ! - cell_mask(ibear, irange) = 1 if data exists, 0 otherwise
  ! - Grid positions are stored as (ibear, iranger) pairs
  ! - Missing velocity flags (VFLG) are checked: 0=good, 128=suspect
  !
  ! ========================================================================

  subroutine codar_populate_grid(CODAR,verbosity)
    
    class(type_codar), intent(inout)  :: CODAR
    integer, intent(in)               :: verbosity
    
    ! Local variables
    integer :: itable, irow, nrows, ncols
    integer :: i, irange, ibear
    real(dp) :: lon, lat, velo, range_km, bearing_deg
    integer :: vflg
    real(dp) :: range_center, bearing_center
    integer :: id_lon, id_lat, id_range, id_bear, id_vflg
    integer :: id_u, id_v, id_velo
    
    ! Find the LLUV radial map table (usually table 1)
    itable = find_table_by_type(CODAR, 'LLUV')
    
    if (itable <= 0) then
        write(*,'(A)') 'ERROR: No LLUV table found in CODAR data'
        return
    endif
    
    ! Initialize mask and arrays
    CODAR%cell_mask = 0
    CODAR%cell_lon  = 0.0_dp
    CODAR%cell_lat  = 0.0_dp
    CODAR%cell_vel  = -999.0_dp
    CODAR%cos_theta = 0.0_dp
    CODAR%sin_theta = 0.0_dp
    if (allocated(CODAR%sigma_vr)) CODAR%sigma_vr = 0.0_dp
    
    nrows = CODAR%Table(itable)%nrows
    ncols = CODAR%Table(itable)%ncols
    
    ! Get column indices for required fields
    ! Column mapping per LLUVSPEC 1.51:
    ! 1=Longitude, 2=Latitude, 16=Range(km), 17=Bearing(deg)
    ! 3=U(cm/s), 4=V(cm/s), 18=Velocity(cm/s), 5=VectorFlag
    
    id_lon   = CODAR%Table(itable)%id('LOND')  ! Longitude
    id_lat   = CODAR%Table(itable)%id('LATD')  ! Latitude
    id_range = CODAR%Table(itable)%id('RNGE')  ! Range (km)
    id_bear  = CODAR%Table(itable)%id('BEAR')  ! Bearing (True degrees)
    id_vflg  = CODAR%Table(itable)%id('VFLG')  ! Vector flag
    id_u     = CODAR%Table(itable)%id('VELU')  ! U component (cm/s)
    id_v     = CODAR%Table(itable)%id('VELV')  ! V component (cm/s)
    id_velo  = CODAR%Table(itable)%id('VELO')  ! Speed (cm/s)
    
    ! Loop through all data rows
    do irow = 1, nrows
        
        ! Get location data from table
        lon  = CODAR%Table(itable)%data(irow, id_lon)
        lat  = CODAR%Table(itable)%data(irow, id_lat)
        velo = CODAR%Table(itable)%data(irow, id_velo) 
        bearing_deg = CODAR%Table(itable)%data(irow, id_bear)
        range_km    = CODAR%Table(itable)%data(irow, id_range)
        vflg        = nint(CODAR%Table(itable)%data(irow, id_vflg))
        
        ! Check for valid data (vector flag should be 0 for good data)
        ! Vector flag = 128 indicates suspect or secondary radial data
        if (vflg /= 0) then
            ! Skip this measurement (marked as suspect or invalid)
            cycle
        endif
        
        ! Find the grid indices (ibear, irange) closest to this measurement
        ! Since the table provides actual lon/lat, use those for grid position
        ! But also validate against range/bearing metadata
        
        call find_grid_indices(CODAR, bearing_deg, range_km, ibear, irange)
        
        ! Check if indices are within valid range
        if (irange < 1 .or. irange > CODAR%nranges .or. &
            ibear < 1 .or. ibear > CODAR%nbearings) then
            write(*,'(A,I5,A,I5,A,I5,A,I5)') &
                'WARNING: Grid indices out of bounds: irange=', irange, &
                ' ibear=', ibear, ' (nranges=', CODAR%nranges, &
                ' nbearings=', CODAR%nbearings, ')'
            cycle
        endif
        
        ! Store grid position
        CODAR%cell_lon(ibear,irange)  = lon
        CODAR%cell_lat(ibear,irange)  = lat
        CODAR%cell_vel(ibear,irange)  = velo / 100.0_dp ! in meter/seconds
        CODAR%cos_theta(ibear,irange) = cos(bearing_deg*deg2rad)
        CODAR%sin_theta(ibear,irange) = sin(bearing_deg*deg2rad)
        CODAR%cell_mask(ibear,irange) = 1  ! Mark as having data
        
        ! Optional: store velocity uncertainty if available
        ! (column 6 = Spatial Quality could be used as proxy for uncertainty)
        
    enddo
    
    if (verbosity.ge.2) then
      write(*,'(A,I5,A)') 'Populated CODAR grid with', & 
                          count(CODAR%cell_mask == 1), &
                          ' valid cells'
    endif
    
  end subroutine codar_populate_grid

  ! ========================================================================
  ! Helper function: Find table by type string
  ! ========================================================================

  function find_table_by_type(CODAR, ttype) result(idx)
    
    class(type_codar), intent(in)     :: CODAR
    character(len=*), intent(in)      :: ttype
    integer                           :: idx
    integer                           :: i
    
    idx = 0
    do i = 1, CODAR%ntables
        if (trim(CODAR%Table(i)%Type) == trim(ttype)) then
            idx = i
            return
        endif
    enddo
    
  end function find_table_by_type

  ! ========================================================================
  ! Helper subroutine: Find grid indices from lon/lat/range/bearing
  ! ========================================================================
  !
  ! Strategy:
  ! 1. Use the range and bearing from the table as primary indices
  ! 2. Fall back to lon/lat matching if range/bearing fails
  ! 3. Implement nearest-neighbor or interpolation as needed
  !
  ! ========================================================================

  subroutine find_grid_indices(CODAR, bearing_deg, range_km, &
                               ibear, irange, status)
    
    class(type_codar), intent(in)     :: CODAR
    real(dp), intent(in)              :: bearing_deg, range_km
    integer, intent(out)              :: ibear, irange
    integer, intent(out), optional    :: status
    integer                           :: i, j
    real(dp)                          :: range_m, dist_min, dist
   
    if (present(status)) status = 0
 
    ! Convert range in km to meters
    range_m = range_km * 1000.0_dp
    
    ! Find range index: closest range cell
    irange = 0
    dist_min = huge(1.0_dp)
    do i = 1, CODAR%nranges
        dist = abs(CODAR%range(i) - range_m)
        if (dist < dist_min) then
            dist_min = dist
            irange = i
        endif
    enddo
    
    ! Find bearing index: closest bearing cell
    ibear = 0
    dist_min = huge(1.0_dp)
    do j = 1, CODAR%nbearings
        ! Handle bearing wrap-around (0° = 360°)
        dist = abs(CODAR%bearing(j) - bearing_deg)
        if (dist > 180.0_dp) then
            dist = 360.0_dp - dist
        endif
        if (dist < dist_min) then
            dist_min = dist
            ibear = j
        endif
    enddo
   
    if (irange.le.0.or.ibear.le.0) then
      if (present(status)) then
        status = -1
      else
        write(*,*) 'WARNING: Could not find valid grid indices'
      endif
    endif 
    
end subroutine find_grid_indices

  subroutine codar_row2grid(CODAR, irow, ibear, irange, status)
    class(type_codar), intent(in)     :: CODAR
    integer, intent(in)               :: irow
    integer, intent(out)              :: ibear,irange
    integer, intent(out), optional    :: status

    integer                           :: itable,id_bear,id_range
    real(dp)                          :: bearing_deg,range_km
    real(dp)                          :: dist_min, dist
   
    if (present(status)) status = 0

    itable      = find_table_by_type(CODAR, 'LLUV')
    id_bear     = CODAR%Table(itable)%id('BEAR')  ! Bearing (True degrees)
    id_range    = CODAR%Table(itable)%id('RNGE')  ! Range (km)
    bearing_deg = CODAR%Table(itable)%data(irow, id_bear)
    range_km    = CODAR%Table(itable)%data(irow, id_range)
    call find_grid_indices(CODAR, bearing_deg, range_km, ibear, irange, status)

  end subroutine codar_row2grid


end module module_codar
