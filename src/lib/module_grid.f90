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

module module_grid

use netcdf
use module_types
use module_constants
use module_tools
use module_time

implicit none

character(len=*), parameter         :: standard_time_units = 'seconds since 2000-01-01 00:00:00'
character(len=*), parameter         :: standard_calendar   = 'gregorian'

! ... Structure for variables
! ...
type type_ncvar
  character(len=maxlen)                          :: name=''
  character(len=maxlen)                          :: units=''
  integer                                        :: ndims=0
  integer                                        :: natts=0
  integer                                        :: type=0
  logical                                        :: dimx = .False.
  logical                                        :: dimy = .False.
  logical                                        :: dimz = .False.
  logical                                        :: dimt = .False.
  logical                                        :: missing=.False.
  integer, dimension(:), pointer                 :: dimids
  character(len=1)                               :: axis=''
  real(dp)                                       :: missing_value=-999.0D0
  real(dp)                                       :: add_offset=0.0D0
  real(dp)                                       :: scale_factor=1.0D0
  real(dp), dimension(:,:,:), pointer            :: mask
end type type_ncvar

! ... Structure for grids
! ...
type type_ncgrid
  logical                                        :: grid2d = .False.
  logical                                        :: Cartesian = .False.
  logical                                        :: Stationary = .False.
  logical                                        :: Climatology = .False.
  character(len=maxlen)                          :: filename=''
  character(len=maxlen)                          :: xname='None'
  character(len=maxlen)                          :: yname='None'
  character(len=maxlen)                          :: zname='None'
  character(len=maxlen)                          :: tname='None'
  character(len=maxlen)                          :: varname='None'
  integer                                        :: fid=0
  integer                                        :: ndims=0
  integer                                        :: nvars=0
  integer                                        :: natts=0
  integer                                        :: idi,idj,idk,idl
  integer                                        :: idx,idy,idz,idt
  integer                                        :: unlimid=-1
  integer                                        :: varid=-1
  integer                                        :: file_nx,ia,ib,nx
  integer                                        :: file_ny,ja,jb,ny
  integer                                        :: file_nz,ka,kb,nz
  integer                                        :: file_nt,la,lb,nt
  integer, dimension(:), pointer                 :: dlen
  real(dp)                                       :: file_xmin,xmin
  real(dp)                                       :: file_xmax,xmax
  real(dp)                                       :: file_ymin,ymin
  real(dp)                                       :: file_ymax,ymax
  real(dp)                                       :: file_zmin,zmin
  real(dp)                                       :: file_zmax,zmax
  real(dp)                                       :: file_tmin,tmin
  real(dp)                                       :: file_tmax,tmax
  real(dp)                                       :: Dconversion
  real(dp)                                       :: Rconversion
  real(dp)                                       :: Xsign = 1.0D0
  real(dp)                                       :: Ysign = 1.0D0
  real(dp)                                       :: Zsign = 1.0D0
  real(dp), dimension(:), pointer                :: file_lon1,lon1,x1
  real(dp), dimension(:), pointer                :: file_lat1,lat1,y1
  real(dp), dimension(:,:), pointer              :: file_lon2,lon2,x2
  real(dp), dimension(:,:), pointer              :: file_lat2,lat2,y2
  real(dp), dimension(:), pointer                :: file_z,z
  real(dp), dimension(:), pointer                :: file_t,t
  type(type_date), dimension(:), pointer         :: file_date,date

  character(len=maxlen), dimension(:), pointer   :: dname
  character(len=maxlen), dimension(:), pointer   :: vname
  type(type_ncvar), dimension(:), pointer        :: var
  logical                                        :: time_arbitrary = .True.
  character(len=maxlen)                          :: time_units = ''
  character(len=maxlen)                          :: calendar = ''

  contains
    procedure                   :: open          => grid_open
    procedure                   :: scan          => grid_scan
    procedure                   :: crop          => grid_crop
    procedure                   :: show          => grid_show
    procedure                   :: interpol      => grid_interpol
    procedure                   :: locate        => grid_locate  
    procedure                   :: read3D        => grid_read3D
    procedure                   :: read2D        => grid_read2D
    procedure                   :: point_type    => grid_point_type

end type type_ncgrid

contains
  ! ...
  ! ==================================================================
  ! ==================================================================
  ! ...
  subroutine grid_open(GRD,filename)

    class(type_ncgrid), intent(inout)            :: GRD
    character(len=*), intent(in)                 :: filename

    ! ... Local variables
    ! ...
    integer err,fid,ndims,nvars,natts,unlimid,iwork
    integer mydim,myvar,myatt
    integer vtype,vndims,vnatts,dimids(10)
    real(dp) rwork
    character(len=maxlen) word
    character(1) vaxis


    !write(*,*) 'Opening grid file: ', trim(filename)
    err = NF90_OPEN(filename,NF90_NOWRITE,fid)
    call cdf_error(err,'opening '//trim(filename))

    GRD%filename = trim(filename)
    GRD%fid      = fid

    ! ... Read axes:
    ! ...
    err = NF90_INQUIRE(fid,ndims,nvars,natts,unlimid)

    GRD%ndims = ndims
    GRD%nvars = nvars
    GRD%natts = natts

    allocate(GRD%dlen(ndims))
    allocate(GRD%dname(ndims))
    do mydim=1,ndims
      word = ''
      err = NF90_INQUIRE_DIMENSION(fid,mydim,word,iwork)
      GRD%dname(mydim) = trim(word)
      GRD%dlen(mydim)  = iwork
    enddo

    allocate(GRD%var(nvars))
    allocate(GRD%vname(nvars))
    do myvar=1,nvars
      word = ''
      err = NF90_INQUIRE_VARIABLE (fid,myvar,word,vtype,vndims,dimids,vnatts)
      GRD%vname(myvar) = trim(word)
      GRD%var(myvar)%name = trim(word)
      GRD%var(myvar)%ndims = vndims
      GRD%var(myvar)%natts = vnatts
      GRD%var(myvar)%type  = vtype
      allocate(GRD%var(myvar)%dimids(vndims))
      GRD%var(myvar)%dimids = dimids(1:vndims)

      ! ... Check for axis
      ! ...
      vaxis = ''
      err  = NF90_GET_ATT(fid,myvar,'axis',vaxis)
      if (err.EQ.NF90_NOERR) GRD%var(myvar)%axis = uppercase(vaxis)

      ! ... Check for missing value
      ! ... Priority: _FillValue, missing_value
      ! ...
      err  = NF90_GET_ATT(fid,myvar,'_FillValue',rwork)
      if (err.eq.NF90_NOERR) then
        GRD%var(myvar)%missing = .True.
        GRD%var(myvar)%missing_value   = rwork
      else
        err = NF90_GET_ATT (fid,myvar,'missing_value',rwork)
        if (err.eq.NF90_NOERR) then
          GRD%var(myvar)%missing = .True.
          GRD%var(myvar)%missing_value   = rwork
        else
          GRD%var(myvar)%missing = .False.
        endif
      endif

      ! ... Check for units, offset and scale factors
      ! ...
      word = ''
      err = NF90_GET_ATT (fid,myvar,'units',word)
      if (err.eq.NF90_NOERR) GRD%var(myvar)%units = trim(word)
      err = NF90_GET_ATT (fid,myvar,'add_offset',rwork)
      if (err.eq.NF90_NOERR) GRD%var(myvar)%add_offset = rwork
      err = NF90_GET_ATT (fid,myvar,'scale_factor',rwork)
      if (err.eq.NF90_NOERR) GRD%var(myvar)%scale_factor = rwork

    enddo


  end subroutine grid_open
  ! ...
  ! ==================================================================
  ! ...
  subroutine grid_scan(GRD,xname,yname,zname,tname,myunits,mycalendar)

    class(type_ncgrid), intent(inout)            :: GRD
    character(len=*), intent(in)                 :: xname
    character(len=*), intent(in)                 :: yname
    character(len=*), intent(in)                 :: zname
    character(len=*), intent(in)                 :: tname
    character(len=*), intent(in)                 :: myunits
    character(len=*), intent(in)                 :: mycalendar

    ! ... Local variables
    ! ...
    integer i,j,err,ndims,var
    integer nx,ny,nz,nt
    real(dp) Coordinate_transform
    real(dp), dimension(:), allocatable          :: wrk1d
    real(dp), dimension(:,:), allocatable        :: wrk2d
    real(dp), dimension(:,:,:), allocatable      :: wrk3d
    character(len=maxlen) time_units,calendar

    if (GRD%Cartesian) then
      Coordinate_transform = 1.0D0
      GRD%Dconversion      = 1.0D0
      GRD%Rconversion      = 1.0D0
    else
      Coordinate_transform = deg2rad
      GRD%Dconversion      = deg2rad
      GRD%Rconversion      = rad2deg
    endif

    ! ... X variable:
    ! ...
    GRD%idx = -1
    GRD%xname = 'NONE'
    if (len_trim(xname).gt.0) then
      ! ... The variable name has been given by the user
      if (xname.ne.'NONE') then
        err = NF90_INQ_VARID(GRD%fid,xname,i)
        if (err.eq.NF90_NOERR) then
          GRD%idx = i
          GRD%xname = xname
          GRD%var(i)%axis = 'X'
        else
          call crash('Variable '//trim(xname)//' not found in file '//trim(GRD%filename))
        endif
      endif
    else
      ! ... The variable name has not been given and we check for Axis att.
      do i=1,GRD%nvars
        if (GRD%var(i)%axis.eq.'X') then
          GRD%idx = i
          GRD%xname = trim(GRD%vname(i))
          exit
        endif
      enddo
    endif
          
    ! ... Y variable:
    ! ...
    GRD%idy = -1
    GRD%yname = 'NONE'
    if (len_trim(yname).gt.0) then
      ! ... The variable name has been given by the user
      if (yname.ne.'NONE') then
        err = NF90_INQ_VARID(GRD%fid,yname,i)
        if (err.eq.NF90_NOERR) then
          GRD%idy = i
          GRD%yname = yname
          GRD%var(i)%axis = 'Y'
        else
          call crash('Variable '//trim(yname)//' not found in file '//trim(GRD%filename))
        endif
      endif
    else
      ! ... The variable name has not been given and we check for Axis att.
      do i=1,GRD%nvars
        if (GRD%var(i)%axis.eq.'Y') then
          GRD%idy = i
          GRD%yname = trim(GRD%vname(i))
          exit
        endif
      enddo
    endif
          
    ! ... Z variable:
    ! ...
    GRD%idz = -1
    GRD%zname = 'NONE'
    if (len_trim(zname).gt.0) then
      ! ... The variable name has been given by the user
      if (zname.ne.'NONE') then
        err = NF90_INQ_VARID(GRD%fid,zname,i)
        if (err.eq.NF90_NOERR) then
          GRD%idz = i
          GRD%zname = zname
          GRD%var(i)%axis = 'Z'
        else
          call crash('Variable '//trim(zname)//' not found in file '//trim(GRD%filename))
        endif
      endif
    else
      ! ... The variable name has not been given and we check for Axis att.
      do i=1,GRD%nvars
        if (GRD%var(i)%axis.eq.'Z') then
          GRD%idz = i
          GRD%zname = trim(GRD%vname(i))
          exit
        endif
      enddo
    endif
          
    ! ... T variable:
    ! ...
    GRD%idt = -1
    GRD%tname = 'NONE'
    if (len_trim(tname).gt.0) then
      ! ... The variable name has been given by the user
      if (tname.ne.'NONE') then
        err = NF90_INQ_VARID(GRD%fid,tname,i)
        if (err.eq.NF90_NOERR) then
          GRD%idt = i
          GRD%tname = tname
          GRD%var(i)%axis = 'T'
        else
          call crash('Variable '//trim(tname)//' not found in file '//trim(GRD%filename))
        endif
      endif
    else
      ! ... The variable name has not been given and we check for Axis att.
      do i=1,GRD%nvars
        if (GRD%var(i)%axis.eq.'T') then
          GRD%idt = i
          GRD%tname = trim(GRD%vname(i))
          exit
        endif
      enddo
    endif

    ! ... Now that we have identified the axis variables, get the
    ! ... information about the system size
    ! ... X
    if (GRD%idx.eq.-1) then
      GRD%idi = -1
      GRD%file_nx  =  1
      GRD%ia  =  1
      GRD%ib  =  1
      GRD%nx  =  1
      GRD%grid2d = .False.
      allocate(GRD%file_lon1(1))
      allocate(GRD%lon1(1))
      allocate(GRD%x1(1))
      GRD%file_lon1(1) = 0.0D0
      GRD%lon1(1) = 0.0D0
      GRD%x1(1) = 0.0D0
      GRD%file_xmin = 0.0D0
      GRD%file_xmax = 0.0D0
    else
      ndims = GRD%var(GRD%idx)%ndims
      select case (ndims)
        case (1)
          GRD%grid2d = .False.
          GRD%idi = GRD%var(GRD%idx)%dimids(1)
          nx = GRD%dlen(GRD%idi)
          GRD%file_nx = nx
          GRD%ia = 1
          GRD%ib = nx
          GRD%nx = nx
        case (2)
          GRD%grid2d = .True.
          GRD%idi = GRD%var(GRD%idx)%dimids(1)
          GRD%idj = GRD%var(GRD%idx)%dimids(2)
          nx = GRD%dlen(GRD%idi)
          ny = GRD%dlen(GRD%idj)
          GRD%file_nx = nx
          GRD%file_ny = ny
          GRD%ia = 1
          GRD%ib = nx
          GRD%nx = nx
          GRD%ja = 1
          GRD%jb = ny
          GRD%ny = ny
        case default
          call crash('Invalid dimensions for '//trim(xname)//' in file '//trim(GRD%filename))
      end select
      if (GRD%grid2d) then
        allocate(GRD%file_lon2(nx,ny))
        allocate(GRD%lon2(nx,ny))
        allocate(GRD%x2(nx,ny))
        err = NF90_GET_VAR(GRD%fid,GRD%idx,GRD%file_lon2)
        call cdf_error(err,'Reading two-dimensional X grid')
        GRD%lon2(:,:) = GRD%file_lon2(:,:)
        GRD%x2(:,:) = Coordinate_transform*GRD%lon2(:,:)
        GRD%file_xmin = minval(GRD%x2)
        GRD%file_xmax = maxval(GRD%x2)
        if (GRD%nx.gt.1) GRD%Xsign = sign(1.0D0,GRD%x2(2,1)-GRD%x2(1,1)) 
      else
        allocate(GRD%file_lon1(nx))
        allocate(GRD%lon1(nx))
        allocate(GRD%x1(nx))
        err = NF90_GET_VAR(GRD%fid,GRD%idx,GRD%file_lon1)
        call cdf_error(err,'Reading one-dimensional X grid')
        GRD%lon1(:) = GRD%file_lon1(:)
        GRD%x1(:) = Coordinate_transform*GRD%lon1(:)
        GRD%file_xmin = minval(GRD%x1)
        GRD%file_xmax = maxval(GRD%x1)
        if (GRD%nx.gt.1) GRD%Xsign = sign(1.0D0,GRD%x1(2)-GRD%x1(1)) 
      endif
    endif

    ! ... Y
    if (GRD%idy.eq.-1) then
      GRD%idj = -1
      GRD%file_ny  =  1
      GRD%ja  =  1
      GRD%jb  =  1
      GRD%ny  =  1
      GRD%grid2d = .False.
      allocate(GRD%file_lat1(1))
      allocate(GRD%lat1(1))
      allocate(GRD%y1(1))
      GRD%file_lat1(1) = 0.0D0
      GRD%lat1(1) = 0.0D0
      GRD%y1(1) = 0.0D0
      GRD%file_ymin = 0.0D0
      GRD%file_ymax = 0.0D0
    else
      ndims = GRD%var(GRD%idy)%ndims
      if (ndims.ne.GRD%var(GRD%idy)%ndims) then
        call crash('Incompatible grid dimensions in file '//trim(GRD%filename))
      endif
      select case (ndims)
        case (1)
          GRD%grid2d = .False.
          GRD%idj = GRD%var(GRD%idy)%dimids(1)
          ny = GRD%dlen(GRD%idj)
          GRD%file_ny = ny
          GRD%ja = 1
          GRD%jb = ny
          GRD%ny = ny
        case (2)
          GRD%grid2d = .True.
          if (GRD%idi.ne.GRD%var(GRD%idy)%dimids(1).or.       &
              GRD%idj.ne.GRD%var(GRD%idy)%dimids(2)) then
            call crash('Incompatible 2D grid dimensions in file '//trim(GRD%filename))
          endif
        case default
          call crash('Invalid dimensions for '//trim(yname)//' in file '//trim(GRD%filename))
      end select
      if (GRD%grid2d) then
        allocate(GRD%file_lat2(nx,ny))
        allocate(GRD%lat2(nx,ny))
        allocate(GRD%y2(nx,ny))
        err = NF90_GET_VAR(GRD%fid,GRD%idy,GRD%file_lat2)
        call cdf_error(err,'Reading two-dimensional Y grid')
        GRD%lat2(:,:) = GRD%file_lat2(:,:)
        GRD%y2(:,:) = Coordinate_transform*GRD%lat2(:,:)
        GRD%file_ymin = minval(GRD%y2)
        GRD%file_ymax = maxval(GRD%y2)
        if (GRD%ny.gt.1) GRD%Ysign = sign(1.0D0,GRD%y2(1,2)-GRD%y2(1,1)) 
      else
        allocate(GRD%file_lat1(ny))
        allocate(GRD%lat1(ny))
        allocate(GRD%y1(ny))
        err = NF90_GET_VAR(GRD%fid,GRD%idy,GRD%file_lat1)
        call cdf_error(err,'Reading one-dimensional Y grid')
        GRD%lat1(:) = GRD%file_lat1(:)
        GRD%y1(:) = Coordinate_transform*GRD%lat1(:)
        GRD%file_ymin = minval(GRD%y1)
        GRD%file_ymax = maxval(GRD%y1)
        if (GRD%ny.gt.1) GRD%Ysign = sign(1.0D0,GRD%y1(2)-GRD%y1(1)) 
      endif
    endif

    ! ... Z
    if (GRD%idz.eq.-1) then
      GRD%idk = -1
      GRD%file_nz  =  1
      GRD%ka  =  1
      GRD%kb  =  1
      GRD%nz  =  1
      allocate(GRD%file_z(1))
      allocate(GRD%z(1))
      GRD%file_z(1) = 0.0D0
      GRD%z(1) = 0.0D0
      GRD%file_zmin = 0.0D0
      GRD%file_zmax = 0.0D0
    else
      ndims = GRD%var(GRD%idz)%ndims
      if (ndims.ne.1) then
        call crash('Invalid dimensions for '//trim(zname)//' in file '//trim(GRD%filename))
      endif
      GRD%idk = GRD%var(GRD%idz)%dimids(1)
      nz = GRD%dlen(GRD%idk)
      GRD%file_nz = nz
      GRD%ka = 1
      GRD%kb = nz
      GRD%nz = nz
      allocate(GRD%file_z(nz))
      allocate(GRD%z(nz))
      err = NF90_GET_VAR(GRD%fid,GRD%idz,GRD%file_z)
      call cdf_error(err,'Reading Z levels')
      GRD%z(:) = -abs(GRD%file_z(:))
      GRD%file_zmin = minval(GRD%z)
      GRD%file_zmax = maxval(GRD%z)
      if (GRD%nz.gt.1) GRD%Zsign = sign(1.0D0,GRD%z(2)-GRD%z(1)) 
    endif
       
    ! ... T
    if (GRD%idt.eq.-1) then
      GRD%idl = -1
      GRD%file_nt  =  1
      GRD%la  =  1
      GRD%lb  =  1
      GRD%nt  =  1
      allocate(GRD%file_t(1))
      allocate(GRD%file_date(1))
      allocate(GRD%t(1))
      allocate(GRD%date(1))
      GRD%file_t(1) = 0.0D0
      GRD%t(1) = GRD%file_t(1)
      GRD%file_tmin = 0.0D0
      GRD%file_tmax = 0.0D0
      if (len_trim(myunits).gt.0) then
        GRD%time_units = trim(myunits)
      else
        GRD%time_units = trim(standard_time_units)
      endif
      if (len_trim(mycalendar).gt.0) then
        GRD%calendar = trim(mycalendar)
      else
        GRD%calendar = 'gregorian'
      endif
      GRD%time_arbitrary = .True.
      GRD%file_date(1) = num2date(GRD%file_t(1),units=GRD%time_units,calendar=GRD%calendar)
      GRD%date(1) = GRD%file_date(1)
    else
      ndims = GRD%var(GRD%idt)%ndims
      if (ndims.ne.1) then
        call crash('Invalid dimensions for '//trim(tname)//' in file '//trim(GRD%filename))
      endif
      GRD%idl = GRD%var(GRD%idt)%dimids(1)
      nt = GRD%dlen(GRD%idl)
      GRD%file_nt = nt
      GRD%la = 1
      GRD%lb = nt
      GRD%nt = nt
      allocate(GRD%file_t(nt))
      allocate(GRD%file_date(nt))
      allocate(GRD%t(nt))
      allocate(GRD%date(nt))
      err = NF90_GET_VAR(GRD%fid,GRD%idt,GRD%file_t)
      call cdf_error(err,'Reading T levels')

      ! ... Now, convert time to dates:
      ! ... 
      time_units = ''
      calendar = ''
      if (len_trim(myunits).eq.0) then
        err = NF90_GET_ATT (GRD%fid,GRD%idt,'units',time_units)
        if (err.ne.NF90_NOERR) then
          time_units = 'seconds'
          GRD%time_arbitrary = .True.
        else
          GRD%time_arbitrary = .False.
        endif
      else
        time_units = trim(myunits)
        GRD%time_arbitrary = .False.
      endif

      if (len_trim(mycalendar).eq.0) then
        err = NF90_GET_ATT (GRD%fid,GRD%idt,'calendar',calendar)
        if (err.ne.NF90_NOERR) calendar = 'gregorian'
      else
        calendar = trim(mycalendar)
      endif
      
      do i=1,GRD%nt
        GRD%file_date(i) = num2date(GRD%file_t(i),units=time_units,calendar=calendar)
        GRD%date(i) = num2date(GRD%file_t(i),units=time_units,calendar=calendar)
        GRD%t(i)    = date2num(GRD%date(i),units=standard_time_units)
      enddo
      GRD%time_units = trim(time_units)
      GRD%calendar   = trim(calendar)
      GRD%file_tmin = minval(GRD%t)
      GRD%file_tmax = maxval(GRD%t)
    endif
       
    GRD%xmin = GRD%file_xmin
    GRD%xmax = GRD%file_xmax
    GRD%ymin = GRD%file_ymin
    GRD%ymax = GRD%file_ymax
    GRD%zmin = GRD%file_zmin
    GRD%zmax = GRD%file_zmax
    GRD%tmin = GRD%file_tmin
    GRD%tmax = GRD%file_tmax

    !... We alreay have almost everything
    ! ...
    do i=1,GRD%nvars
      do j=1,GRD%var(i)%ndims
        if (GRD%var(i)%dimids(j).eq.GRD%idi) GRD%var(i)%dimx = .True.
        if (GRD%var(i)%dimids(j).eq.GRD%idj) GRD%var(i)%dimy = .True.
        if (GRD%var(i)%dimids(j).eq.GRD%idk) GRD%var(i)%dimz = .True.
        if (GRD%var(i)%dimids(j).eq.GRD%idl) GRD%var(i)%dimt = .True.
      enddo 
    enddo

    ! ... Finally, we get the mask for each variable
    ! ...
    do var=1,GRD%nvars
      nx = -1; ny = -1; nz = -1
      if (len_trim(GRD%var(var)%axis).eq.0) then 
        ! The variable is not an axis
        if (GRD%var(var)%dimx) then
          nx = GRD%nx
        else
          nx = 1
        endif
        if (GRD%var(var)%dimy) then
          ny = GRD%ny
        else
          ny = 1
        endif
        if (GRD%var(var)%dimz) then
          nz = GRD%nz
        else
          nz = 1
        endif

        allocate (GRD%var(var)%mask(nx,ny,nz))
        GRD%var(var)%mask(:,:,:) = 1.0_dp          ! By default, all valid points

        ! ... Read an initial field:
        ! ...
        if (GRD%var(var)%ndims.eq.1) then

          if (GRD%var(var)%dimx) then
            ! -------------------------------- X
            allocate (wrk1d(nx))
            err = NF90_GET_VAR(GRD%fid,var,wrk1d)
            call cdf_error(err,'ERROR in getting X mask')
            where(wrk1d.eq.GRD%var(var)%missing_value) GRD%var(var)%mask(:,1,1) = 0.0_dp
            deallocate(wrk1d)
          else if (GRD%var(var)%dimy) then
            ! -------------------------------- Y
            allocate (wrk1d(ny))
            err = NF90_GET_VAR(GRD%fid,var,wrk1d)
            call cdf_error(err,'ERROR in getting Y mask')
            where(wrk1d.eq.GRD%var(var)%missing_value) GRD%var(var)%mask(1,:,1) = 0.0_dp
            deallocate(wrk1d)
          else if (GRD%var(var)%dimz) then
            ! -------------------------------- Z
            allocate (wrk1d(nz))
            err = NF90_GET_VAR(GRD%fid,var,wrk1d)
            call cdf_error(err,'ERROR in getting Z mask')
            where(wrk1d.eq.GRD%var(var)%missing_value) GRD%var(var)%mask(1,1,:) = 0.0_dp
            deallocate(wrk1d)
          endif

        else if (GRD%var(var)%ndims.eq.2) then
        
          if (GRD%var(var)%dimx) then
            if (GRD%var(var)%dimy) then
              ! -------------------------------- XY
              allocate (wrk2d(nx,ny))
              err = NF90_GET_VAR(GRD%fid,var,wrk2d)
              call cdf_error(err,'ERROR in getting XY mask')
              where(wrk2d.eq.GRD%var(var)%missing_value) GRD%var(var)%mask(:,:,1) = 0.0_dp
              deallocate(wrk2d)
            else if (GRD%var(var)%dimz) then
              ! -------------------------------- XZ
              allocate (wrk2d(nx,nz))
              err = NF90_GET_VAR(GRD%fid,var,wrk2d)
              call cdf_error(err,'ERROR in getting XZ mask')
              where(wrk2d.eq.GRD%var(var)%missing_value) GRD%var(var)%mask(:,1,:) = 0.0_dp
              deallocate(wrk2d)
            else
              ! -------------------------------- XT
              allocate (wrk1d(nx))
              err = NF90_GET_VAR(GRD%fid,var,wrk1d,[1,1],[nx,1])
              call cdf_error(err,'ERROR in getting X mask')
              where(wrk1d.eq.GRD%var(var)%missing_value) GRD%var(var)%mask(:,1,1) = 0.0_dp
              deallocate(wrk1d)
            endif
          else if (GRD%var(var)%dimy) then
            if (GRD%var(var)%dimz) then
              ! -------------------------------- YZ
              allocate (wrk2d(ny,nz))
              err = NF90_GET_VAR(GRD%fid,var,wrk2d)
              call cdf_error(err,'ERROR in getting YZ mask')
              where(wrk2d.eq.GRD%var(var)%missing_value) GRD%var(var)%mask(1,:,:) = 0.0_dp
              deallocate(wrk2d)
            else
              ! -------------------------------- YT
              allocate (wrk1d(ny))
              err = NF90_GET_VAR(GRD%fid,var,wrk1d,[1,1],[ny,1])
              call cdf_error(err,'ERROR in getting Y mask')
              where(wrk1d.eq.GRD%var(var)%missing_value) GRD%var(var)%mask(1,:,1) = 0.0_dp
              deallocate(wrk1d)
            endif
          else
            ! -------------------------------- ZT
            allocate (wrk1d(nz))
            err = NF90_GET_VAR(GRD%fid,var,wrk1d,[1,1],[nz,1])
            call cdf_error(err,'ERROR in getting Z mask')
            where(wrk1d.eq.GRD%var(var)%missing_value) GRD%var(var)%mask(1,1,:) = 0.0_dp
            deallocate(wrk1d)
          endif    

        else if (GRD%var(var)%ndims.eq.3) then

          if (.not.GRD%var(var)%dimx) then
            ! -------------------------------- YZT
            allocate (wrk2d(ny,nz))
            err = NF90_GET_VAR(GRD%fid,var,wrk2d,[1,1,1],[ny,nz,1])
            call cdf_error(err,'ERROR in getting YZ mask')
            where(wrk2d.eq.GRD%var(var)%missing_value) GRD%var(var)%mask(1,:,:) = 0.0_dp
            deallocate(wrk2d)
          else if (.not.GRD%var(var)%dimy) then
            ! -------------------------------- XZT
            allocate (wrk2d(nx,nz))
            err = NF90_GET_VAR(GRD%fid,var,wrk2d,[1,1,1],[nx,nz,1])
            call cdf_error(err,'ERROR in getting XZ mask')
            where(wrk2d.eq.GRD%var(var)%missing_value) GRD%var(var)%mask(:,1,:) = 0.0_dp
            deallocate(wrk2d)
          else if (.not.GRD%var(var)%dimz) then
            ! -------------------------------- XYT
            allocate (wrk2d(nx,ny))
            err = NF90_GET_VAR(GRD%fid,var,wrk2d,[1,1,1],[nx,ny,1])
            call cdf_error(err,'ERROR in getting XY mask')
            where(wrk2d.eq.GRD%var(var)%missing_value) GRD%var(var)%mask(:,:,1) = 0.0_dp
            deallocate(wrk2d)
          else 
            ! -------------------------------- XYZ
            allocate (wrk3d(nx,ny,nz))
            err = NF90_GET_VAR(GRD%fid,var,wrk3d)
            call cdf_error(err,'ERROR in getting XYZ mask')
            where(wrk3d.eq.GRD%var(var)%missing_value) GRD%var(var)%mask = 0.0_dp
            deallocate(wrk3d)
          endif
        
        else if (GRD%var(var)%ndims.eq.4) then

          ! -------------------------------- XYZ
          allocate (wrk3d(nx,ny,nz))
          err = NF90_GET_VAR(GRD%fid,var,wrk3d,[1,1,1,1],[nx,ny,nz,1])
          call cdf_error(err,'ERROR in getting XYZ mask')
          where(wrk3d.eq.GRD%var(var)%missing_value) GRD%var(var)%mask = 0.0_dp
          deallocate(wrk3d)

        else

          call crash('Error while getting mask. Invalid number of dimensions')

        endif

      endif
    enddo

  end subroutine grid_scan
  ! ...
  ! ==================================================================
  ! ...
  subroutine grid_crop(GRD,xmin,xmax,ymin,ymax,tmin,tmax)

    class(type_ncgrid), intent(inout)            :: GRD
    real(dp), intent(in)                         :: xmin,xmax
    real(dp), intent(in)                         :: ymin,ymax
    real(dp), intent(in)                         :: tmin,tmax

    ! ... Local variables
    ! ...
    integer i,j,io,jo,il,jl,ni,nj
    real(dp), dimension(:,:), allocatable        :: tmp2
    real(dp), dimension(:), allocatable          :: tmp1
    type(type_date), dimension(:), allocatable   :: tmpd

    if (xmin.lt.GRD%xmin) call crash('Domain xmin < Grid xmin')
    if (xmax.gt.GRD%xmax) call crash('Domain xmax < Grid xmax')
    if (ymin.lt.GRD%ymin) call crash('Domain ymin < Grid ymin')
    if (ymax.gt.GRD%ymax) call crash('Domain ymax < Grid ymax')
    if (tmin.lt.GRD%tmin) call crash('Domain tmin < Grid tmin')
    if (tmax.gt.GRD%tmax) call crash('Domain tmax < Grid tmax')


    ! ... X - Y cropping
    ! ...
    if (GRD%grid2d) then
      ! ... 2D-grid
      ! ...
      io = -1; jo = -1
      do j=1,GRD%ny-1
      do i=1,GRD%nx-1
        if (GRD%x2(i,j).lt.xmin.and.GRD%x2(i+1,j).ge.xmin.and. &
            GRD%y2(i,j).lt.ymin.and.GRD%y2(i,j+1).ge.ymin) then
              io = i
              jo = j
              goto 1
        endif
      enddo
      enddo
      1 continue
      io = max(io,1); jo = max(jo,1)

      il = GRD%nx+1; jl = GRD%ny+1
      do j=1,GRD%ny-1
      do i=1,GRD%nx-1
        if (GRD%x2(i,j).lt.xmax.and.GRD%x2(i+1,j).ge.xmax.and. &
            GRD%y2(i,j).lt.ymax.and.GRD%y2(i,j+1).ge.ymax) then
              il = i + 1
              jl = j + 1
              goto 2
        endif
      enddo
      enddo
      2 continue
      il = min(il,GRD%nx); jl = min(jl,GRD%ny)

      GRD%ia = io; GRD%ib = il; ni = il-io+1; GRD%nx = ni
      GRD%ja = jo; GRD%jb = jl; nj = jl-jo+1; GRD%ny = nj

      allocate(tmp2(ni,nj))

      do j=jo,jo+nj-1
      do i=io,io+ni-1
        tmp2(i-io+1,j-jo+1) = GRD%x2(i,j)
      enddo
      enddo
      deallocate(GRD%x2)
      allocate(GRD%x2(ni,nj))
      do j=1,nj
      do i=1,ni
        GRD%x2(i,j) = tmp2(i,j)
      enddo
      enddo

      do j=jo,jo+nj-1
      do i=io,io+ni-1
        tmp2(i-io+1,j-jo+1) = GRD%y2(i,j)
      enddo
      enddo
      deallocate(GRD%y2)
      allocate(GRD%y2(ni,nj))
      do j=1,nj
      do i=1,ni
        GRD%y2(i,j) = tmp2(i,j)
      enddo
      enddo

      deallocate(tmp2)

      GRD%xmin = minval(GRD%x2)
      GRD%xmax = maxval(GRD%x2)
      GRD%ymin = minval(GRD%y2)
      GRD%ymax = maxval(GRD%y2)

    else
      ! ... 1D-grid
      ! ...
      io = -1; jo = -1
      do i=1,GRD%nx-1
        if (GRD%x1(i).lt.xmin.and.GRD%x1(i+1).ge.xmin) then
          io = i
          exit
        endif
      enddo
      do j=1,GRD%ny-1
        if (GRD%y1(j).lt.ymin.and.GRD%y1(j+1).ge.ymin) then
          jo = j
          exit
        endif
      enddo
      io = max(io,1); jo = max(jo,1)

      il = GRD%nx+1; jl = GRD%ny+1
      do i=1,GRD%nx-1
        if (GRD%x1(i).lt.xmax.and.GRD%x1(i+1).ge.xmax) then
          il = i+1
          exit
        endif
      enddo
      do j=1,GRD%ny-1
        if (GRD%y1(j).lt.ymax.and.GRD%y1(j+1).ge.ymax) then
          jl = j+1
          exit
        endif
      enddo
      il = min(il,GRD%nx); jl = min(jl,GRD%ny)

      GRD%ia = io; GRD%ib = il; ni = il-io+1; GRD%nx = ni
      GRD%ja = jo; GRD%jb = jl; nj = jl-jo+1; GRD%ny = nj

      allocate(tmp1(max(ni,nj)))

      do i=io,io+ni-1
        tmp1(i-io+1) = GRD%x1(i)
      enddo
      deallocate(GRD%x1)
      allocate(GRD%x1(ni))
      do i=1,ni
        GRD%x1(i) = tmp1(i)
      enddo

      do j=jo,jo+nj-1
        tmp1(j-jo+1) = GRD%y1(j)
      enddo
      deallocate(GRD%y1)
      allocate(GRD%y1(nj))
      do j=1,nj
        GRD%y1(j) = tmp1(j)
      enddo

      deallocate(tmp1)

      GRD%xmin = minval(GRD%x1)
      GRD%xmax = maxval(GRD%x1)
      GRD%ymin = minval(GRD%y1)
      GRD%ymax = maxval(GRD%y1)
    endif

    ! ... T cropping
    ! ...
    io = -1
    do i=1,GRD%nt-1
      if (GRD%t(i).lt.tmin.and.GRD%t(i+1).ge.tmin) then
        io = i
        exit
      endif
    enddo
    io = max(io,1)

    il = GRD%nt+1
    do i=1,GRD%nt-1
      if (GRD%t(i).lt.tmax.and.GRD%t(i+1).ge.tmax) then
        il = i+1
        exit
      endif
    enddo
    il = min(il,GRD%nt)

    GRD%la = io; GRD%lb = il; ni = il-io+1; GRD%nt = ni

    allocate(tmp1(ni))
    allocate(tmpd(ni))

    do i=io,io+ni-1
      tmp1(i-io+1) = GRD%t(i)
      tmpd(i-io+1) = GRD%date(i)
    enddo
    deallocate(GRD%t)
    deallocate(GRD%date)
    allocate(GRD%t(ni))
    allocate(GRD%date(ni))
    do i=1,ni
      GRD%t(i) = tmp1(i)
      GRD%date(i) = tmpd(i)
    enddo

    deallocate (tmp1)
    deallocate (tmpd)

    ! ... Update selected bounds
    ! ...
    GRD%zmin = minval(GRD%z)
    GRD%zmax = maxval(GRD%z)
    GRD%tmin = minval(GRD%t)
    GRD%tmax = maxval(GRD%t)

  end subroutine grid_crop
  ! ...
  ! ==================================================================
  ! ...
  subroutine grid_show(GRD,Label)

    class(type_ncgrid), intent(in)               :: GRD
    character(len=*), intent(in)                 :: Label

    ! ... Local variables
    ! ...
    logical x,y,z,t
    type(type_date) di,df

    x = .not.(GRD%idx.eq.-1)
    y = .not.(GRD%idy.eq.-1)
    z = .not.(GRD%idz.eq.-1)
    t = .not.(GRD%idt.eq.-1)

    write(*,*)
    write(*,*) '-------------------------------'
    write(*,*) trim(Label)
    write(*,*) 'Filename         : ', trim(GRD%filename)
    write(*,*) 'Axes X, Y, Z, T  : ',x,y,z,t
    if (GRD%Cartesian) then
      write(*,*) 'Coordinates      : Cartesian '
    else
      write(*,*) 'Coordinates      : Spherical '
    endif
    if (GRD%grid2d) then
      write(*,*) 'Horizontal grid  : 2D '
    else
      write(*,*) 'Horizontal grid  : 1D '
    endif
    if (GRD%time_arbitrary) then
      write(*,*) 'Time units       : seconds '
    else
      write(*,*) 'Time units       : '//trim(GRD%time_units)
      write(*,*) 'Time calendar    : '//trim(GRD%calendar)
    endif


    write(*,*) 'Variable         : ', trim(GRD%varname)
    write(*,*) 
    write(*,*) '== File == '
    write(*,*) 'Nx, Ny, Nz, Nt   : ', GRD%file_Nx, GRD%file_Ny, GRD%file_nz, GRD%file_nt
    if (GRD%Cartesian) then
      write(*,*) 'X(1)  - X(Nx)    : ', GRD%Rconversion*GRD%file_xmin,GRD%Rconversion*GRD%file_xmax
      write(*,*) 'Y(1)  - Y(Ny)    : ', GRD%Rconversion*GRD%file_ymin,GRD%Rconversion*GRD%file_ymax
    else
      write(*,*) 'West  - East     : ', GRD%Rconversion*GRD%file_xmin,GRD%Rconversion*GRD%file_xmax
      write(*,*) 'South - North    : ', GRD%Rconversion*GRD%file_ymin,GRD%Rconversion*GRD%file_ymax
    endif
    write(*,*) 'Z(1)  - Z(Nz)    : ', GRD%file_z(1),GRD%file_z(GRD%file_nz)
!    write(*,*) 'West  (deg, rad) : ', GRD%Rconversion*GRD%file_xmin,GRD%file_xmin
!    write(*,*) 'East  (deg, rad) : ', GRD%Rconversion*GRD%file_xmax,GRD%file_xmax
!    write(*,*) 'South (deg, rad) : ', GRD%Rconversion*GRD%file_ymin,GRD%file_ymin
!    write(*,*) 'North (deg, rad) : ', GRD%Rconversion*GRD%file_ymax,GRD%file_ymax
    if (GRD%time_arbitrary) then
      write(*,*) 'Initial time     : ', GRD%file_tmin
      write(*,*) 'Final   time     : ', GRD%file_tmax
    else
      di = GRD%file_date(1)
      df = GRD%file_date(GRD%nt)
      write(*,*) 'Initial date     :    ', trim(di%iso())
      write(*,*) 'Final   date     :    ', trim(df%iso())
    endif
    write(*,*) 
    write(*,*) '== Selected region == '
    write(*,*) 'Nx, Ny, Nz, Nt   : ', GRD%Nx, GRD%Ny, GRD%nz, GRD%nt
    if (GRD%Cartesian) then
      write(*,*) 'X(1)  - X(Nx)    : ', GRD%Rconversion*GRD%xmin,GRD%Rconversion*GRD%xmax
      write(*,*) 'Y(1)  - Y(Ny)    : ', GRD%Rconversion*GRD%ymin,GRD%Rconversion*GRD%ymax
    else
      write(*,*) 'West  - East     : ', GRD%Rconversion*GRD%xmin,GRD%Rconversion*GRD%xmax
      write(*,*) 'South - North    : ', GRD%Rconversion*GRD%ymin,GRD%Rconversion*GRD%ymax
    endif
      write(*,*) 'Z(1)  - Z(Nz)    : ', GRD%z(1),GRD%z(GRD%nz)
!    write(*,*) 'West  (deg, rad) : ', GRD%Rconversion*GRD%xmin,GRD%xmin
!    write(*,*) 'East  (deg, rad) : ', GRD%Rconversion*GRD%xmax,GRD%xmax
!    write(*,*) 'South (deg, rad) : ', GRD%Rconversion*GRD%ymin,GRD%ymin
!    write(*,*) 'North (deg, rad) : ', GRD%Rconversion*GRD%ymax,GRD%ymax
    if (GRD%time_arbitrary) then
      write(*,*) 'Initial time     : ', GRD%tmin
      write(*,*) 'Final   time     : ', GRD%tmax
    else
      di = GRD%date(1)
      df = GRD%date(GRD%nt)
      write(*,*) 'Initial date     :    ', trim(di%iso())
      write(*,*) 'Final   date     :    ', trim(df%iso())
    endif
    write(*,*)

  end subroutine grid_show
  ! ...
  ! ==================================================================
  ! ...
  function grid_locate(GRD,xo,IND0) result(IND)

    class(type_ncgrid), intent(in)                       :: GRD
    real(dp), dimension(:), intent(in)                   :: xo
    integer, dimension(size(xo)), intent(in), optional   :: IND0
    integer, dimension(size(xo))                         :: IND

    ! ... Local variables
    ! ...
    integer ndim,i,j,k
    logical apply_update

    ndim = size(xo)
    if (ndim.lt.2.or.ndim.gt.4) call crash('Invalid number of dimensions in grid_locate')

    if (present(IND0)) then
      apply_update = .True.
    else
      apply_update = .False.
    endif

    if (apply_update) then

      i = IND0(1)
      j = IND0(2)
      if (ndim.eq.3) k = IND0(3)
     
      if (GRD%grid2d) then

      else
        if (i.gt.1) then
          do while(GRD%Xsign*xo(1).lt.GRD%Xsign*GRD%x1(i))
            i = i - 1
          enddo
        endif
        if (j.gt.1) then
          do while(GRD%Ysign*xo(2).lt.GRD%Ysign*GRD%y1(j))
            j = j - 1
          enddo
        endif
        if (i.lt.GRD%nx-1) then
          do while(GRD%Xsign*xo(1).ge.GRD%Xsign*GRD%x1(i+1))
            i = i + 1
          enddo
        endif
        if (j.lt.GRD%ny-1) then
          do while(GRD%Ysign*xo(2).ge.GRD%Ysign*GRD%y1(j+1))
            j = j + 1
          enddo
        endif

      endif

      IND(1) = i
      IND(2) = j

      if (ndim.eq.3) then
        if (k.gt.1) then
          do while(GRD%Zsign*xo(3).lt.GRD%Zsign*GRD%z(k))
            k = k - 1
          enddo
        endif
        if (k.lt.GRD%nz-1) then
          do while(GRD%Zsign*xo(3).ge.GRD%Zsign*GRD%z(k+1))
            k = k + 1
          enddo
        endif
        IND(3) = k
      endif


    else
      IND(:) = -1

      if (GRD%grid2d) then
      
        do j=1,GRD%ny
        do i=1,GRD%nx
          if ((GRD%Xsign*GRD%x2(i,j).gt.GRD%Xsign*xo(1)).and.(GRD%Ysign*GRD%y2(i,j).gt.GRD%Ysign*xo(2))) then
            IND(1) = i-1
            IND(2) = j-1
            goto 1
          endif
        enddo
        enddo
 1      continue
      else
        if (GRD%Xsign*xo(1).lt.GRD%Xsign*GRD%x1(GRD%nx)) then
          do i=1,GRD%nx
            if (GRD%Xsign*GRD%x1(i).gt.GRD%Xsign*xo(1)) then
              IND(1) = i-1
              exit
            endif
          enddo
        else if (GRD%Xsign*xo(1).gt.GRD%Xsign*GRD%x1(GRD%nx)) then
          IND(1) = GRD%nx
        else
          IND(1) = GRD%nx - 1
        endif
        
        if (GRD%Ysign*xo(2).lt.GRD%Ysign*GRD%y1(GRD%ny)) then
          do j=1,GRD%ny
            if (GRD%Ysign*GRD%y1(j).gt.GRD%Ysign*xo(2)) then
              IND(2) = j-1
              exit
            endif
          enddo
        else if (GRD%Ysign*xo(2).gt.GRD%Ysign*GRD%y1(GRD%ny)) then
          IND(2) = GRD%ny
        else
          IND(2) = GRD%ny - 1
        endif

      endif
 

      if (ndim.eq.3) then
        if (GRD%Zsign*xo(3).gt.GRD%Zsign*GRD%z(GRD%nz)) then
          IND(3) = GRD%nz
        else
          do k=1,GRD%nz
            !if (GRD%Zsign*GRD%z(k).le.GRD%Zsign*xo(3)) then
            if (GRD%z(k).le.xo(3)) then
              IND(3) = k
              exit
            endif
          enddo
        endif
      endif

    endif

  end function grid_locate
  ! ...
  ! ==================================================================
  ! ...
  function grid_interpol(GRD,F,xo,LOC,SingleLayer) result(Fo)

    ! ... 2D/3D interpolation:
    ! ... 2D interpolation uses bilinear interpolation
    ! ... Vertical interpolation uses linear interpolation
    class(type_ncgrid), intent(in)               :: GRD
    real(dp), dimension(:,:,:), intent(in)       :: F
    real(dp), dimension(3), intent(in)           :: xo
    integer, dimension(3), intent(in)            :: LOC
    logical, intent(in)                          :: SingleLayer
    real(dp)                                     :: Fo

    ! ... Local variables:
    ! ...
    logical Good1,Good2,Good3,Good4
    integer ndim,i1,j1,k1,i2,j2,k2,i3,j3,i4,j4,Ngoods
    real(dp) x1,y1,x2,y2,t,u,F1,F2

    Fo = 0.0D0

    i1 = LOC(1); i2 = i1+1; i3 = i2;   i4 = i1
    j1 = LOC(2); j2 = j1;   j3 = j2+1; j4 = j3
    k1 = LOC(3); k2 = k1-1                       ! In Z k1 is the level below ZO

    Good1 = i1.GE.1      .and. j1.GE.1
    Good2 = i2.LE.GRD%nx .and. j2.GE.1
    Good3 = i3.LE.GRD%nx .and. j3.LE.GRD%ny
    Good4 = i4.GE.1      .and. j4.LE.GRD%ny
    Ngoods = count([Good1,Good2,Good3,Good4])

    !print*, i4, j4, '     -     ', i3,j3
    !print*, i1, j1, '     -     ', i2,j2
    !print*, 'k1, k2 = ', k1, k2
    !print*, 'Goods: ', Good1, Good2, Good3, Good4
    !print*, 'Ngoods = ', Ngoods
    !print*, 'SingleLayer: ', SingleLayer
    
    if (Ngoods.eq.4) then
      if (GRD%grid2d) then
          x1 = GRD%x2(i1,j1); x2 = GRD%x2(i2,j2)
          y1 = GRD%y2(i1,j1); y2 = GRD%y2(i4,j4)
      else
          x1 = GRD%x1(i1); x2 = GRD%x1(i2)
          y1 = GRD%y1(j1); y2 = GRD%y1(j4)
      endif
      t = (xo(1)-x1)/(x2-x1)
      u = (xo(2)-y1)/(y2-y1)
      !print*, 't, u = ', t, u
      F1 = (1.0D0-t)*(1.0D0-u)*F(i1,j1,k1) + &
                   t*(1.0D0-u)*F(i2,j2,k1) + &
                           t*u*F(i3,j3,k1) + &
                   (1.0D0-t)*u*F(i4,j4,k1)
      
    else
      ! ... The interpolation point is no longer and
      ! ... interior point. 
      stop 'Ngoods'
      F1 = 0.0D0
      if (Good1) F1 = F1 + F(i1,j1,k1)
      if (Good2) F1 = F1 + F(i2,j2,k1)
      if (Good3) F1 = F1 + F(i3,j3,k1)
      if (Good4) F1 = F1 + F(i4,j4,k1)
      F1 = F1 / Ngoods
    endif
    if (SingleLayer) then
      Fo = F1
    else
      if (k2.gt.0) then
        if (Ngoods.eq.4) then
          F2 = (1.0D0-t)*(1.0D0-u)*F(i1,j1,k2) + &
                       t*(1.0D0-u)*F(i2,j2,k2) + &
                               t*u*F(i3,j3,k2) + &
                       (1.0D0-t)*u*F(i4,j4,k2)
        else
          ! ... The interpolation point is no longer and
          ! ... interior point. 
          F2 = 0.0D0
          if (Good1) F2 = F2 + F(i1,j1,k2)
          if (Good2) F2 = F2 + F(i2,j2,k2)
          if (Good3) F2 = F2 + F(i3,j3,k2)
          if (Good4) F2 = F2 + F(i4,j4,k2)
          F2 = F2 / Ngoods
        endif
        Fo = F1 + (F2-F1)*(xo(3)-GRD%z(k1))/(GRD%z(k2)-GRD%z(k1))
      else
        Fo = F1
      endif
      !print*, 'Vertical interpolation F1, F2, Fo : ', F1, F2, Fo
    endif


  end function grid_interpol
  ! ...
  ! ==================================================================
  ! ...
  function grid_read3D(GRD,Varid,step,missing_value,verbose) result(F)


  ! ... Subroutine to read a 3D grid from a NetCDF
  ! ... Optionally it can read only one unique layer, specified by
  ! ... the optional argument layer.
  ! ...
  ! ... To read a 3D snapshot:
  ! ...               F = GRD%read(Varid,GRD%nz,step=STEP)
  ! ...
  ! ... To read the first layer snapshot:
  ! ...               F = GRD%read(Varid,1,step=STEP,layer=1)
  ! ...

  class(type_ncgrid), intent(in)             :: GRD
  integer, intent(in)                        :: Varid
  integer, intent(in), optional              :: step
  real(dp), intent(in), optional             :: missing_value
  real(dp), dimension(GRD%nx,GRD%ny,GRD%Nz)  :: F
  logical, optional                          :: verbose

  ! ... Local variables
  ! ...
  logical verb
  integer l,err
  real(dp) spval

  if (present(verbose)) then
    verb = verbose
  else
    verb = .False.
  endif

  if (present(missing_value)) then
    spval = missing_value
  else
    spval = GRD%var(Varid)%missing_value
  endif


  if (present(step)) then
    l  = GRD%la + step - 1
  else
    l = GRD%la
  endif

  if (l.GT.GRD%nt) call crash('in GRID_READ3D step > Nt')

  if (verb) write(*,'(T2, "READ_GRID3D: varid = ",i2," step = ", i3)') varid, step

  if (GRD%var(Varid)%ndims.EQ.1) then
    ! ... Read 1-Dim : ERROR
    ! ...
    call crash('in GRID_READ3D. Only one dimension')

  else if (GRD%var(Varid)%ndims.EQ.2) then
    ! ... Read 2-Dim :
    ! ...
    call crash('in GRID_READ3D. Only two dimensions')

  else if (GRD%var(Varid)%ndims.EQ.3) then
    ! ... Read 3-Dim :
    ! ...
    err = NF90_GET_VAR(GRD%fid,Varid,F,(/GRD%ia,GRD%ja,GRD%ka/),(/GRD%nx,GRD%ny,GRD%nz/))
    call cdf_error(err,'Reading three dimensional field in GRID_READ3D')

  else if (GRD%var(Varid)%ndims.EQ.4) then
    ! ... Read 4-Dim :
    ! ...
    err = NF90_GET_VAR(GRD%fid,Varid,F,(/GRD%ia,GRD%ja,GRD%ka,l/),(/GRD%nx,GRD%ny,GRD%nz,1/))
    call cdf_error(err,'Reading four dimensional field in GRID_READ3D')

  else

    call crash('Too many dimensions in GRID_READ3D')

  endif

  ! ... Scale the field:
  ! ...
  if (GRD%var(Varid)%missing) then
    where((F.eq.GRD%var(Varid)%missing_value))
      F = spval
    elsewhere
      F = GRD%var(Varid)%add_offset + GRD%var(Varid)%scale_factor*F
    endwhere
  else
    F = GRD%var(Varid)%add_offset + GRD%var(Varid)%scale_factor*F
  endif

  return
  end function grid_read3D
  ! ...
  ! ==================================================================
  ! ...
  function grid_read2D(GRD,Varid,step,layer,missing_value) result(F)


  ! ... Subroutine to read a 2D grid from a NetCDF
  ! ... Optionally it can read only one unique layer, specified by
  ! ... the optional argument layer.
  ! ...
  ! ... To read a 3D snapshot:
  ! ...               F = GRD%read(Varid,GRD%nz,step=STEP)
  ! ...
  ! ... To read the first layer snapshot:
  ! ...               F = GRD%read(Varid,1,step=STEP,layer=1)
  ! ...

  class(type_ncgrid), intent(in)         :: GRD
  integer, intent(in)                    :: Varid
  integer, intent(in), optional          :: step
  integer, intent(in), optional          :: layer
  real(dp), intent(in), optional         :: missing_value
  real(dp), dimension(GRD%nx,GRD%ny)     :: F

  ! ... Local variables
  ! ...
  integer k,l,err
  real(dp) spval

  if (present(missing_value)) then
    spval = missing_value
  else
    spval = GRD%var(Varid)%missing_value
  endif

  if (present(layer)) then
    k = GRD%ka + layer - 1
  else
    k = GRD%ka
  endif

  if (present(step)) then
    l = GRD%la + step - 1
  else
    l = GRD%la 
  endif

  if (k.GT.GRD%nz) call crash('in GRID_READ2D layer > Nz')
  if (l.GT.GRD%nt) call crash('in GRID_READ2D step > Nt')

  write(*,'(T2, "READ_GRID2D: varid = ",i2," layer = ",i2," step = ", i3)') varid, layer, step

  if (GRD%var(Varid)%ndims.EQ.1) then
    ! ... Read 1-Dim : ERROR
    ! ...
    call crash('in GRID_READ2D. Only one dimension')

  else if (GRD%var(Varid)%ndims.EQ.2) then
    ! ... Read 2-Dim :
    ! ...
    err = NF90_GET_VAR(GRD%fid,Varid,F,(/GRD%ia,GRD%ja/),(/GRD%nx,GRD%ny/))
    call cdf_error(err,'Reading two dimensional firld in GRID_READ2D')

  else if (GRD%var(Varid)%ndims.EQ.3) then
    ! ... Read 3-Dim :
    ! ...
    if (GRD%var(Varid)%dimz) then
      err = NF90_GET_VAR(GRD%fid,Varid,F,(/GRD%ia,GRD%ja,k/),(/GRD%nx,GRD%ny,1/))
      call cdf_error(err,'Reading three dimensional field in GRID_READ2D')
    else if (GRD%var(Varid)%dimt) then
      err = NF90_GET_VAR(GRD%fid,Varid,F,(/GRD%ia,GRD%ja,l/),(/GRD%nx,GRD%ny,1/))
      call cdf_error(err,'Reading three dimensional field in GRID_READ2D')
    else
      call crash('Reading three dimensional field in GRID_READ2D: No Z or T')
    endif

  else if (GRD%var(Varid)%ndims.EQ.4) then
    ! ... Read 4-Dim :
    ! ...
    err = NF90_GET_VAR(GRD%fid,Varid,F,(/GRD%ia,GRD%ja,k,l/),(/GRD%nx,GRD%ny,1,1/))
    call cdf_error(err,'Reading four dimensional field in GRID_READ2D')

  else

    call crash('Too many dimensions in GRID_READ2D')

  endif

  ! ... Scale the field:
  ! ...
  if (GRD%var(Varid)%missing) then
    where((F.eq.GRD%var(Varid)%missing_value))
      F = spval
    elsewhere
      F = GRD%var(Varid)%add_offset + GRD%var(Varid)%scale_factor*F
    endwhere
  else
    F = GRD%var(Varid)%add_offset + GRD%var(Varid)%scale_factor*F
  endif

  return
  end function grid_read2D
  ! ...
  ! ==================================================================
  ! ...
  integer function grid_point_type(GRD,xo,yo,ozo) result(ptype)
  ! ... Bilinear interpolation of the mask
  ! ... It returns 1 if point is in water
  ! ... It returns 0 if point is on land
  ! ... It retusn -1 if point "out of bounds" (no interior point)
  ! ...
    class(type_ncgrid), intent(in)               :: GRD
    real(dp), intent(in)                         :: xo,yo
    real(dp), optional, intent(in)               :: ozo

    ! ... Local variables
    ! ...
    integer i,j,k,ndims
    integer, dimension(3)                        :: IP
    real(dp), dimension(3)                       :: XX
    real(dp) zo
    real(dp) y1,y2,y3,y4
    real(dp) d1,d2,d3,d4
    real(dp) w1,w2,w3,w4
    real(dp) fo

    ndims = 0
    if (GRD%var(GRD%varid)%dimx) ndims = ndims + 1
    if (GRD%var(GRD%varid)%dimy) ndims = ndims + 1
    if (GRD%var(GRD%varid)%dimz) ndims = ndims + 1

    if (ndims.eq.3) then
      if (present(ozo)) then
        zo = ozo
      else
        zo = 0.0D0
      endif
    endif

    ptype = -1
    if (GRD%grid2d) then
      if (xo.le.GRD%x2(1,1)) return
      if (xo.ge.GRD%x2(GRD%nx,GRD%ny)) return
      if (yo.le.GRD%y2(1,1)) return
      if (yo.ge.GRD%y2(GRD%nx,GRD%ny)) return
    else
      if (xo.le.GRD%x1(1)) return
      if (xo.ge.GRD%x1(GRD%nx)) return
      if (yo.le.GRD%y1(1)) return
      if (yo.ge.GRD%y1(GRD%ny)) return
    endif

    ptype = 0
    if (zo.le.GRD%z(GRD%nz)) return

    if (ndims.eq.2) then
      print*, '2D mask interpolation. CAREFULL !!!!!'
      XX = [xo,yo,0.0D0]
      IP(1:2) = GRD%locate([xo,yo])
      IP(3)   = 1
      print*, 'XX: ', XX
      print*, 'IP: ', IP
    else
      XX = [xo,yo,zo]
      IP = GRD%locate(XX)
    endif
    
    fo = GRD%interpol(GRD%var(GRD%varid)%mask(GRD%ia:GRD%ib,GRD%ja:GRD%jb,GRD%ka:GRD%kb), &
                      XX,IP,SingleLayer=.True.)
    if (fo.gt.0.5D0) then
      ptype = 1
      return
    else
      ptype = 0
      return
    endif


  end function grid_point_type
  ! ...
  ! ==================================================================
  ! ...
end module module_grid
