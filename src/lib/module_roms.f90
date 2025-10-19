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
! - roms_open                                                              !
! - roms_state_read                                                        !
! - roms_zselect                                                           !
! - roms_vprofile                                                          !
! - roms_zfields_save                                                      !
! -------------------------------------------------------------------------!

module module_roms

use module_types
use module_time
use module_nc
use module_tools

implicit none

private
public type_roms

type type_roms_grid
  integer                           :: L         = -1
  integer                           :: M         = -1
  integer                           :: N         = -1
  integer                           :: nx_r      = -1
  integer                           :: nx_u      = -1
  integer                           :: nx_v      = -1
  integer                           :: nx_psi    = -1
  integer                           :: ny_r      = -1
  integer                           :: ny_u      = -1
  integer                           :: ny_v      = -1
  integer                           :: ny_psi    = -1
  integer                           :: nz_r      = -1
  integer                           :: nz_w      = -1
  integer                           :: Spherical = -1
  real(dp)                          :: xlength   = 0.0_dp
  real(dp)                          :: ylength   = 0.0_dp
  logical                           :: has_mask = .false.  
  integer                           :: Vtransform, Vstretching  
  real(dp)                          :: theta_s, theta_b, Tcline, hc  
  real(dp), allocatable             :: sc_r(:), sc_w(:)  
  real(dp), allocatable             :: Cs_r(:), Cs_w(:)  
  real(dp), allocatable             :: h(:,:)           ! bathymetry  
  real(dp), allocatable             :: pm(:,:), pn(:,:) ! metric factors (1/dx, 1/dy)  
  real(dp), allocatable             :: x_r(:,:), y_r(:,:) ! grid x, y rho-locations (cell center)
  real(dp), allocatable             :: x_u(:,:), y_u(:,:) ! grid x, y u-locations (cell center)
  real(dp), allocatable             :: x_v(:,:), y_v(:,:) ! grid x, y v-locations (cell center)
  real(dp), allocatable             :: x_psi(:,:), y_psi(:,:) ! grid x, y psi-locations (cell center)
  real(dp), allocatable             :: mask_r(:,:)    ! 0/1  
  real(dp), allocatable             :: Zo_r(:,:,:)
  real(dp), allocatable             :: Zo_w(:,:,:)
  real(dp), allocatable             :: lon_r(:,:), lat_r(:,:)  ! lon/lat values (rho-locations)
  contains
    procedure                       :: show => roms_grid_show
end type type_roms_grid

type type_roms_time
  integer                           :: dimid = -1
  integer                           :: varid = -1
  integer                           :: n     = 0
  character(len=maxlen)             :: units = ""
  character(len=maxlen)             :: calendar = ""
  real(dp), allocatable             :: time(:)
  type(type_date), allocatable      :: date(:)
  contains
    procedure                       :: show => roms_time_show
end type type_roms_time

type type_roms_state
  integer                           :: step
  integer                           :: nx, ny, nz   ! Rho-points
  integer                           :: idz, idt, ids, idu, idv, idw
  real(dp)                          :: time
  type(type_date)                   :: date
  real(dp), allocatable             :: zeta(:,:)
  real(dp), allocatable             :: temp(:,:,:)
  real(dp), allocatable             :: salt(:,:,:)
  real(dp), allocatable             :: u(:,:,:)
  real(dp), allocatable             :: v(:,:,:)
  real(dp), allocatable             :: w(:,:,:)
  real(dp), allocatable             :: z_w(:,:,:)
  real(dp), allocatable             :: z_r(:,:,:)
  real(dp), allocatable             :: Hz(:,:,:)
end type type_roms_state

type type_roms_zfields
  integer                           :: nx, ny, nz   ! Rho-points
  real(dp)                          :: time
  type(type_date)                   :: date
  real(dp), allocatable             :: depths(:)
  real(dp), allocatable             :: zeta(:,:)
  real(dp), allocatable             :: temp(:,:,:)
  real(dp), allocatable             :: salt(:,:,:)
  real(dp), allocatable             :: u(:,:,:)
  real(dp), allocatable             :: v(:,:,:)
  real(dp), allocatable             :: w(:,:,:)
end type type_roms_zfields

type type_roms
  character(len=maxlen)             :: filename   = ""
  character(len=maxlen)             :: grid_file  = ""
  character(len=maxlen)             :: GRID_PATH  = './'
  integer                           :: ncid       = -1
  integer                           :: ntimes     = 0
  integer                           :: ndims      = 0
  integer                           :: nvars      = 0
  integer                           :: natts      = 0
  integer                           :: unlimid    = -1
  integer                           :: time_varid = -1
  real(dp)                          :: missing    = -999.0_dp
  type(type_roms_grid)              :: grid
  type(type_roms_time)              :: time
  type(type_roms_state)             :: state
  type(type_roms_zfields)           :: zfields

  contains
    procedure                       :: open => roms_open
    procedure                       :: read => roms_state_read
    procedure                       :: zselect => roms_zselect
    procedure                       :: vprofile => roms_vprofile
    procedure                       :: zfields_save => roms_zfields_save
end type type_roms

contains
! ...
! ======================================================================
! ======================================================================
! ...
  subroutine roms_open(ROMS,filename,grid_file)
 
  class(type_roms), intent(inout)            :: ROMS
  character(len=*), intent(in)               :: filename
  character(len=*), intent(in), optional     :: grid_file

  ! ... Local variables
  ! ...
  integer ncid,err,ndims,nvars,natts,unlimid
  integer mydim,ilen,myvar,vtype,vndims,dimids(10),vnatts
  integer gid,gidx,gidy
  integer i,j,step
  integer iw1(1)
  real(dp) hwater,fact
  real(dp) rw1(1)
  character(len=maxlen) word,path,base,type

  ! ... Open input file
  ! ...
  err = nf90_open(filename, NF90_NOWRITE, ncid)
  call nc_error(err,'Unable to open input ROMS file')

  ROMS%ncid = ncid

  ! ... Contents
  ! ...
  err = NF90_INQUIRE(ncid,ndims,nvars,natts,unlimid)
  ROMS%ndims = ndims
  ROMS%nvars = nvars
  ROMS%natts = natts
  ROMS%unlimid = unlimid

  ! ... Dimensions
  ! ...
  do mydim=1,ndims
    word = ''
    err = NF90_INQUIRE_DIMENSION(ncid,mydim,word,ilen)
    if (word.eq.'xi_rho')  ROMS%grid%nx_r   = ilen
    if (word.eq.'xi_u')    ROMS%grid%nx_u   = ilen
    if (word.eq.'xi_v')    ROMS%grid%nx_v   = ilen
    if (word.eq.'xi_psi')  ROMS%grid%nx_psi = ilen
    if (word.eq.'eta_rho') ROMS%grid%ny_r   = ilen
    if (word.eq.'eta_u')   ROMS%grid%ny_u   = ilen
    if (word.eq.'eta_v')   ROMS%grid%ny_v   = ilen
    if (word.eq.'eta_psi') ROMS%grid%ny_psi = ilen
    if (word.eq.'s_rho')   ROMS%grid%nz_r   = ilen
    if (word.eq.'s_w')     ROMS%grid%nz_w   = ilen
    if (mydim.eq.unlimid)  ROMS%ntimes      = ilen
  enddo
  ROMS%time%dimid = unlimid
  ROMS%time%n     = ROMS%ntimes
  allocate(ROMS%time%time(ROMS%time%n))
  allocate(ROMS%time%date(ROMS%time%n))

  ROMS%grid%L = ROMS%grid%nx_r - 1
  ROMS%grid%M = ROMS%grid%ny_r - 1
  ROMS%grid%N = ROMS%grid%nz_r

  allocate(ROMS%grid%sc_r(ROMS%grid%nz_r))
  allocate(ROMS%grid%Cs_r(ROMS%grid%nz_r))
  allocate(ROMS%grid%sc_w(0:ROMS%grid%nz_r))
  allocate(ROMS%grid%Cs_w(0:ROMS%grid%nz_r))

  allocate(   ROMS%grid%h(0:ROMS%grid%L,0:ROMS%grid%M))
  allocate(  ROMS%grid%pm(0:ROMS%grid%L,0:ROMS%grid%M))
  allocate(  ROMS%grid%pn(0:ROMS%grid%L,0:ROMS%grid%M))
  allocate( ROMS%grid%x_r(0:ROMS%grid%L,0:ROMS%grid%M))
  allocate( ROMS%grid%y_r(0:ROMS%grid%L,0:ROMS%grid%M))
  allocate( ROMS%grid%x_u(1:ROMS%grid%L,0:ROMS%grid%M))
  allocate( ROMS%grid%y_u(1:ROMS%grid%L,0:ROMS%grid%M))
  allocate( ROMS%grid%x_v(0:ROMS%grid%L,1:ROMS%grid%M))
  allocate( ROMS%grid%y_v(0:ROMS%grid%L,1:ROMS%grid%M))
  allocate(ROMS%grid%Zo_r(0:ROMS%grid%L,0:ROMS%grid%M,ROMS%grid%nz_r))
  allocate(ROMS%grid%Zo_w(0:ROMS%grid%L,0:ROMS%grid%M,0:ROMS%grid%nz_r))
  allocate(ROMS%state%z_w(0:ROMS%grid%L,0:ROMS%grid%M,0:ROMS%grid%nz_r))
  allocate(ROMS%state%z_r(0:ROMS%grid%L,0:ROMS%grid%M,ROMS%grid%nz_r))
  allocate( ROMS%state%Hz(0:ROMS%grid%L,0:ROMS%grid%M,ROMS%grid%nz_r))

  ! ... Grid file
  ! ...
  if (present(grid_file)) then
    ROMS%grid_file = trim(grid_file)
    write(*,*) 'grid_file name '//trim(ROMS%grid_file)//' provided by user'
    call filename_split (ROMS%grid_file,path,base,type)
    ROMS%grid_file = trim(ROMS%GRID_PATH) // '/' //trim(base) // '.' // trim(type)
  else
    err = NF90_GET_ATT(ncid,0,'grd_file',ROMS%grid_file)
    call filename_split (ROMS%grid_file,path,base,type)
    ROMS%grid_file = trim(ROMS%GRID_PATH) // '/' // trim(base) // '.' // trim(type)
    write(*,*) 'grid_file name '//trim(ROMS%grid_file)//' read from ROMS attribute grd_file'
  endif

  allocate(ROMS%grid%lon_r(0:ROMS%grid%L,0:ROMS%grid%M))
  allocate(ROMS%grid%lat_r(0:ROMS%grid%L,0:ROMS%grid%M))

  err = NF90_OPEN(ROMS%grid_file,0,gid)
  call nc_error(err,'Unable to open grid file')
  err = NF90_INQ_VARID(gid,'lon_rho',gidx)
  call nc_error(err,'Variable lon_rho not found in grid file')
  err = NF90_INQ_VARID(gid,'lat_rho',gidy)
  call nc_error(err,'Variable lat_rho not found in grid file')
  err = NF90_GET_VAR(gid,gidx,ROMS%grid%lon_r)
  call nc_error(err,'Unable to read variable lon_rho from grid file')
  err = NF90_GET_VAR(gid,gidy,ROMS%grid%lat_r)
  call nc_error(err,'Unable to read variable lat_rho from grid file')
  err = NF90_CLOSE(gid)


  ! ... Allocate state variables (single time step)
  ! ... nx, ny, nz <= RHO-points
  ! ...
  ROMS%state%nx = ROMS%grid%nx_r
  ROMS%state%ny = ROMS%grid%ny_r
  ROMS%state%nz = ROMS%grid%nz_r
  allocate(ROMS%state%zeta(0:ROMS%grid%L,0:ROMS%grid%M))
  allocate(ROMS%state%temp(0:ROMS%grid%L,0:ROMS%grid%M,ROMS%grid%nz_r))
  allocate(ROMS%state%salt(0:ROMS%grid%L,0:ROMS%grid%M,ROMS%grid%nz_r))
  allocate(   ROMS%state%u(1:ROMS%grid%L,0:ROMS%grid%M,ROMS%grid%nz_r))
  allocate(   ROMS%state%v(0:ROMS%grid%L,1:ROMS%grid%M,ROMS%grid%nz_r))

  ! ... Check for mask:
  ! ...
  if (nc_var_exists(ncid, "mask_rho")) then  
    ROMS%grid%has_mask = .true.  
    allocate(ROMS%grid%mask_r(ROMS%grid%nx_r,ROMS%grid%ny_r))  
  else  
    ROMS%grid%has_mask = .false.  
    stop 'Missing mask_rho'
  end if  

  ! ... Variables
  ! ... Only selected variables
  ! ...
  do myvar=1,nvars
    word = ''
    err = NF90_INQUIRE_VARIABLE (ncid,myvar,word,vtype,vndims,dimids,vnatts)
    if (word.eq.'spherical') then
      err = NF90_GET_VAR(ncid,myvar,iw1)
      ROMS%Grid%Spherical = iw1(1)
    endif
    if (word.eq.'xl') then
      err = NF90_GET_VAR(ncid,myvar,rw1)
      ROMS%Grid%xlength = rw1(1)
    endif
    if (word.eq.'el') then
      err = NF90_GET_VAR(ncid,myvar,rw1)
      ROMS%Grid%ylength = rw1(1)
    endif
    if (word.eq.'Vtransform') then
      err = NF90_GET_VAR(ncid,myvar,iw1)
      ROMS%Grid%Vtransform = iw1(1)
    endif
    if (word.eq.'Vstretching') then
      err = NF90_GET_VAR(ncid,myvar,iw1)
      ROMS%Grid%Vstretching = iw1(1)
    endif
    if (word.eq.'theta_s') then
      err = NF90_GET_VAR(ncid,myvar,rw1)
      ROMS%Grid%theta_s = rw1(1)
    endif
    if (word.eq.'theta_b') then
      err = NF90_GET_VAR(ncid,myvar,rw1)
      ROMS%Grid%theta_b = rw1(1)
    endif
    if (word.eq.'Tcline') then
      err = NF90_GET_VAR(ncid,myvar,rw1)
      ROMS%Grid%Tcline = rw1(1)
    endif
    if (word.eq.'hc') then
      err = NF90_GET_VAR(ncid,myvar,rw1)
      ROMS%Grid%hc = rw1(1)
    endif
    if (word.eq.'s_rho') then
      err = NF90_GET_VAR(ncid,myvar,ROMS%grid%sc_r)
    endif
    if (word.eq.'s_w') then
      err = NF90_GET_VAR(ncid,myvar,ROMS%grid%sc_w)
    endif
    if (word.eq.'Cs_r') then
      err = NF90_GET_VAR(ncid,myvar,ROMS%grid%Cs_r)
    endif
    if (word.eq.'Cs_w') then
      err = NF90_GET_VAR(ncid,myvar,ROMS%grid%Cs_w)
    endif
    if (word.eq.'h') then
      err = NF90_GET_VAR(ncid,myvar,ROMS%grid%h)
    endif
    if (word.eq.'pm') then
      err = NF90_GET_VAR(ncid,myvar,ROMS%grid%pm)
      call nc_error(err,'Unable to read PM')
    endif
    if (word.eq.'pn') then
      err = NF90_GET_VAR(ncid,myvar,ROMS%grid%pn)
      call nc_error(err,'Unable to read PN')
    endif
    if (word.eq.'x_rho') then
      err = NF90_GET_VAR(ncid,myvar,ROMS%grid%x_r)
      call nc_error(err,'Unable to read X_RHO')
    endif
    if (word.eq.'y_rho') then
      err = NF90_GET_VAR(ncid,myvar,ROMS%grid%y_r)
      call nc_error(err,'Unable to read Y_RHO')
    endif
    if (word.eq.'x_u') then
      err = NF90_GET_VAR(ncid,myvar,ROMS%grid%x_u)
      call nc_error(err,'Unable to read X_U')
    endif
    if (word.eq.'y_u') then
      err = NF90_GET_VAR(ncid,myvar,ROMS%grid%y_u)
      call nc_error(err,'Unable to read Y_U')
    endif
    if (word.eq.'x_v') then
      err = NF90_GET_VAR(ncid,myvar,ROMS%grid%x_v)
      call nc_error(err,'Unable to read X_V')
    endif
    if (word.eq.'y_v') then
      err = NF90_GET_VAR(ncid,myvar,ROMS%grid%y_v)
      call nc_error(err,'Unable to read Y_V')
    endif
    if (word.eq.'mask_rho') then
      err = NF90_GET_VAR(ncid,myvar,ROMS%grid%mask_r)
      call nc_error(err,'Unable to read MASK_RHO')
    endif
    if (word.eq.'ocean_time') then
      ROMS%time%varid = myvar
      err = NF90_GET_VAR(ncid,myvar,ROMS%time%time)
      err = NF90_GET_ATT(ncid,myvar,'units',ROMS%time%units)
      err = NF90_GET_ATT(ncid,myvar,'calendar',ROMS%time%calendar)
      do step=1,ROMS%time%n
        ROMS%time%date(step) = num2date(ROMS%time%time(step), &
                                        ROMS%time%units,      &
                                        ROMS%time%calendar)
      enddo
    endif
    if (word.eq.'zeta') then
      ROMS%state%idz = myvar
    endif
    if (word.eq.'temp') then
      ROMS%state%idt = myvar
    endif
    if (word.eq.'salt') then
      ROMS%state%ids = myvar
    endif
    if (word.eq.'u') then
      ROMS%state%idu = myvar
    endif
    if (word.eq.'v') then
      ROMS%state%idv = myvar
    endif
  enddo

  !-----------------------------------------------------------------------
  !  Original formulation: Compute vertical depths (meters, negative) at
  !                        RHO- and W-points, and vertical grid
  !  thicknesses. Various stretching functions are possible.
  !
  !         z_w(x,y,s,t) = Zo_w + zeta(x,y,t) * [1.0 + Zo_w / h(x,y)]
  !                 Zo_w = hc * [s(k) - C(k)] + C(k) * h(x,y)
  !
  !-----------------------------------------------------------------------
  !
  if (ROMS%Grid%Vtransform.eq.1) then

    do j=0,ROMS%Grid%M
    do i=0,ROMS%Grid%L
      hwater = ROMS%Grid%h(i,j)
      ROMS%Grid%Zo_r(i,j,:) = ROMS%Grid%hc*(ROMS%Grid%sc_r(:)-ROMS%Grid%Cs_r(:)) + ROMS%Grid%Cs_r(:)*hwater
      ROMS%Grid%Zo_w(i,j,:) = ROMS%Grid%hc*(ROMS%Grid%sc_w(:)-ROMS%Grid%Cs_w(:)) + ROMS%Grid%Cs_w(:)*hwater
    enddo
    enddo

  !
  !-----------------------------------------------------------------------
  !  New formulation: Compute vertical depths (meters, negative) at
  !                   RHO- and W-points, and vertical grid thicknesses.
  !  Various stretching functions are possible.
  !
  !         z_w(x,y,s,t) = zeta(x,y,t) + [zeta(x,y,t) + h(x,y)] * Zo_w
  !                 Zo_w = [hc * s(k) + C(k) * h(x,y)] / [hc + h(x,y)]
  !
  !-----------------------------------------------------------------------
  !
  else if (ROMS%Grid%Vtransform.eq.2) then

    do j=0,ROMS%Grid%L
    do i=0,ROMS%Grid%M
      hwater = ROMS%Grid%h(i,j)
      fact   = 1.0_dp / ( ROMS%Grid%hc + hwater)
      ROMS%Grid%Zo_r(i,j,:) = fact * (ROMS%Grid%hc * ROMS%Grid%sc_r(:) + ROMS%Grid%Cs_r(:)*hwater )
      ROMS%Grid%Zo_w(i,j,:) = fact * (ROMS%Grid%hc * ROMS%Grid%sc_w(:) + ROMS%Grid%Cs_w(:)*hwater )
    enddo
    enddo

  endif

  end subroutine roms_open
  ! ...
  ! ======================================================================
  ! ...
  subroutine roms_grid_show(GRID)

    class(type_roms_grid), intent(in)          :: GRID 
 
    write(*,*)
    write(*,*) 'Horizontal and Vertical grids'
    write(*,*) 'T dims      : ', GRID%nx_r, GRID%ny_r, GRID%nz_r
    write(*,*) 'U dims      : ', GRID%nx_u,   GRID%ny_u,   GRID%nz_r
    write(*,*) 'V dims      : ', GRID%nx_v,   GRID%ny_v,   GRID%nz_r
    write(*,*) 'Spherical   : ', GRID%Spherical, '   [0-Cartesian, 1-Spherical]'
    write(*,*) 'Vtransform  : ', GRID%Vtransform
    write(*,*) 'Vstretching : ', GRID%Vstretching
    write(*,*) 'Xlen, Ylen  : ', GRID%xlength, GRID%ylength
    write(*,*) 
    write(*,*) 'Hmin, Hmax  : ', minval(GRID%h), maxval(GRID%h)
    write(*,*) 
    write(*,*) 'X   min,max : ', minval(GRID%x_r), maxval(GRID%x_r), ' at RHO-points'
    write(*,*) '            : ', minval(GRID%x_u), maxval(GRID%x_u), ' at U-points'
    write(*,*) '            : ', minval(GRID%x_v), maxval(GRID%x_v), ' at V-points'
    write(*,*) 'Y   min,max : ', minval(GRID%y_r), maxval(GRID%y_r), ' at RHO-points'
    write(*,*) '            : ', minval(GRID%y_u), maxval(GRID%y_u), ' at U-points'
    write(*,*) '            : ', minval(GRID%y_v), maxval(GRID%y_v), ' at V-points'
    write(*,*) 
    write(*,*) 'Lon min,max : ', minval(GRID%lon_r), maxval(GRID%lon_r), ' at RHO-points'
    write(*,*) 'Lat min,max : ', minval(GRID%lat_r), maxval(GRID%lat_r), ' at RHO-points'


  end subroutine roms_grid_show
  ! ...
  ! ======================================================================
  ! ...
  subroutine roms_time_show(TIME)

    class(type_roms_time), intent(in)          :: TIME 
 
    write(*,*)
    write(*,*) 'Number time steps: ', TIME%n
    write(*,*) 'Units            : ', trim(TIME%units)
    write(*,*) 'Calendar         : ', trim(TIME%calendar)
    write(*,*) 'Initial time     : ', TIME%time(1), TIME%date(1)%iso()
    write(*,*) 'Final time       : ', TIME%time(TIME%n), TIME%date(TIME%n)%iso()

  end subroutine roms_time_show
  ! ...
  ! ======================================================================
  ! ...
  subroutine roms_state_read(ROMS,step)

    class(type_roms), intent(inout)             :: ROMS
    integer, intent(in)                         :: step

    ! ... Local variables
    ! ...
    integer i,j,k,L,M
    integer err,nxr,nyr,nzr,nxu,nyu,nxv,nyv
    real(dp) zow,zor,zetaij,hwater,hinv
    !real(kind=8) :: start_time, end_time

    nxr = ROMS%grid%nx_r; nyr = ROMS%grid%ny_r
    nxu = ROMS%grid%nx_u; nyu = ROMS%grid%ny_u
    nxv = ROMS%grid%nx_v; nyv = ROMS%grid%ny_v
    nzr = ROMS%grid%nz_r

    L   = nxr - 1 
    M   = nyr - 1

    ROMS%state%step = step
    ROMS%state%time = ROMS%time%time(step)
    ROMS%state%date = ROMS%time%date(step)

    write(*,*) 'Reading step: ', step, &
               '  Time and date: ', ROMS%state%time, ROMS%state%date%iso()

    err = NF90_GET_VAR(ROMS%ncid,ROMS%state%idz,ROMS%state%zeta,[1,1,step],[nxr,nyr,1])
    call nc_error(err,'Unable to read zeta')

    err = NF90_GET_VAR(ROMS%ncid,ROMS%state%idt,ROMS%state%temp,[1,1,1,step],[nxr,nyr,nzr,1])
    call nc_error(err,'Unable to read temp')

    err = NF90_GET_VAR(ROMS%ncid,ROMS%state%ids,ROMS%state%salt,[1,1,1,step],[nxr,nyr,nzr,1])
    call nc_error(err,'Unable to read salt')

    err = NF90_GET_VAR(ROMS%ncid,ROMS%state%idu,ROMS%state%u,[1,1,1,step],[nxu,nyu,nzr,1])
    call nc_error(err,'Unable to read u')

    err = NF90_GET_VAR(ROMS%ncid,ROMS%state%idv,ROMS%state%v,[1,1,1,step],[nxv,nyv,nzr,1])
    call nc_error(err,'Unable to read v')

    where(ROMS%state%zeta.gt.1000) ROMS%state%zeta = ROMS%missing
    where(ROMS%state%temp.gt.1000) ROMS%state%temp = ROMS%missing
    where(ROMS%state%salt.gt.1000) ROMS%state%salt = ROMS%missing
    where(ROMS%state%u.gt.1000)    ROMS%state%u    = ROMS%missing
    where(ROMS%state%v.gt.1000)    ROMS%state%v    = ROMS%missing

    ! ... Z-depths of S-levels
    ! ...
    if (ROMS%Grid%Vtransform.eq.1) then
      !         z_w(x,y,s,t) = Zo_w + zeta(x,y,t) * [1.0 + Zo_w / h(x,y)]
      !                 Zo_w = hc * [s(k) - C(k)] + C(k) * h(x,y)
      do j=0,M
      do i=0,L
        ROMS%state%z_w(i,j,0) = -ROMS%Grid%h(i,j)
        hwater = ROMS%Grid%h(i,j)
        hinv   = 1.0_dp/hwater
        zetaij = ROMS%state%zeta(i,j)
        do k=1,nzr
          zow = ROMS%grid%Zo_w(i,j,k)
          zor = ROMS%grid%Zo_r(i,j,k)
          ROMS%state%z_w(i,j,k) = zow + zetaij*(1.0_dp+zow*hinv)
          ROMS%state%z_r(i,j,k) = zor + zetaij*(1.0_dp+zor*hinv)
          ROMS%state%Hz(i,j,k)  = ROMS%state%z_w(i,j,k) - ROMS%state%z_w(i,j,k-1)
        enddo
      enddo   
      enddo   
    else if (ROMS%Grid%Vtransform.eq.2) then
      !         z_w(x,y,s,t) = zeta(x,y,t) + [zeta(x,y,t) + h(x,y)] * Zo_w
      !                 Zo_w = [hc * s(k) + C(k) * h(x,y)] / [hc + h(x,y)]

!      call cpu_time(start_time)    ! Start timing
      do j=1,M
      do i=1,L
        hwater = ROMS%Grid%h(i,j)
        zetaij = ROMS%state%zeta(i,j)
        ROMS%state%z_w(i,j,:) = zetaij + (zetaij + hwater)*ROMS%grid%Zo_w(i,j,:)
        ROMS%state%z_r(i,j,:) = zetaij + (zetaij + hwater)*ROMS%grid%Zo_r(i,j,:)
        ROMS%state%Hz(i,j,:)  = ROMS%state%z_w(i,j,1:nzr) - ROMS%state%z_w(i,j,0:nzr-1)
      enddo   
      enddo   
!      call cpu_time(end_time)    ! Start timing
!      print*, end_time - start_time

    endif

  end subroutine roms_state_read
  ! ...
  ! ======================================================================
  ! ...
  subroutine roms_zselect(ROMS,Zout_pos)

    ! ... Transform from original state (C-grid, s-levels) to a A-grid, z-depths
    ! ... The z-values are only defined at rho points, so the velocities U and V
    ! ... must be first set at the rho points. The z(i,j,k) have already been
    ! ... calculated and stored in ROMS%state%z_r(i,j,k) at rho-points.

    class(type_roms), intent(inout)             :: ROMS
    real(dp), dimension(:), intent(in)          :: Zout_pos

    ! ... Local variables
    ! ...
    integer ia,ja,ka,kk,kl,kh,L,M,Lm,Mm
    !integer kl2d(ROMS%grid%nx_r,ROMS%grid%ny_r,size(Zout_pos))  
    integer nxr,nyr,nzr,nzout
    real(dp) Zout(size(Zout_pos))
    real(dp) u1,u2,v1,v2,dz,hh,frac
    real(dp) tl,th,sl,sh,ul,uh,vl,vh

    nxr = ROMS%grid%nx_r; nyr = ROMS%grid%ny_r; nzr = ROMS%grid%nz_r
    L   = ROMS%grid%L; Lm = L - 1
    M   = ROMS%grid%M; Mm = M - 1


    ! ... First, we change the sign of the Zout
    ! ...
    Zout(:) = -abs(Zout_pos(:))     ! Makes sure it is negative
    nzout = size(Zout)

    if (allocated(ROMS%zfields%depths)) deallocate(ROMS%zfields%depths)
    if (allocated(ROMS%zfields%zeta)) deallocate(ROMS%zfields%zeta)
    if (allocated(ROMS%zfields%temp)) deallocate(ROMS%zfields%temp)
    if (allocated(ROMS%zfields%salt)) deallocate(ROMS%zfields%salt)
    if (allocated(ROMS%zfields%u)) deallocate(ROMS%zfields%u)
    if (allocated(ROMS%zfields%v)) deallocate(ROMS%zfields%v)

    ! ... Fiels in A grid
    ! ...
    allocate(ROMS%zfields%depths(nzout))
    allocate(ROMS%zfields%zeta(1:Lm,1:Mm))
    allocate(ROMS%zfields%temp(1:Lm,1:Mm,nzout))
    allocate(ROMS%zfields%salt(1:Lm,1:Mm,nzout))
    allocate(   ROMS%zfields%u(1:Lm,1:Mm,nzout))
    allocate(   ROMS%zfields%v(1:Lm,1:Mm,nzout))

    ROMS%zfields%depths(:) = Zout_pos(:)
    ROMS%zfields%zeta(1:Lm,1:Mm) = ROMS%state%zeta(1:Lm,1:Mm)

    do ja=1,Mm
      do ia=1,Lm

        if (ROMS%grid%mask_r(ia,ja).gt.0.5D0) then

          do ka=1,nzout

            if (Zout(ka).ge.ROMS%state%z_r(ia,ja,nzr)) then

              ! ............................... Level too shallow (assuming fully mixed layer)
              ! ...
              ROMS%zfields%temp(ia,ja,ka) = ROMS%state%temp(ia,ja,nzr)
              ROMS%zfields%salt(ia,ja,ka) = ROMS%state%salt(ia,ja,nzr)

              ! ... U
              ! ...
              u1 = ROMS%state%u(ia  ,ja,nzr)
              u2 = ROMS%state%u(ia+1,ja,nzr)
              if (u1.eq.ROMS%missing.or.u2.eq.ROMS%missing) then
                ROMS%zfields%u(ia,ja,ka) = ROMS%missing
              else
                ROMS%zfields%u(ia,ja,ka) = 0.5_dp * (u1 + u2)
              endif
           
              ! ... V
              ! ...
              v1 = ROMS%state%v(ia,ja  ,nzr)
              v2 = ROMS%state%v(ia,ja+1,nzr)
              if (v1.eq.ROMS%missing.or.v2.eq.ROMS%missing) then
                ROMS%zfields%v(ia,ja,ka) = ROMS%missing
              else
                ROMS%zfields%v(ia,ja,ka) = 0.5_dp * (v1 + v2)
              endif

            else
              kl = 0; kh = 0
              do kk=1,nzr
                if (ROMS%state%z_r(ia,ja,kk).gt.Zout(ka)) then
                  kl = kk - 1
                  kh = kk
                  exit
                endif
              enddo
              if (kl.eq.0) then
                ! ............................... Level too deep
                ! ...
                ROMS%zfields%temp(ia,ja,ka) = ROMS%missing
                ROMS%zfields%salt(ia,ja,ka) = ROMS%missing
                ROMS%zfields%u(ia,ja,ka)    = ROMS%missing
                ROMS%zfields%v(ia,ja,ka)    = ROMS%missing
              else
                ! ............................... Vertical interpolation
                ! ...
                
                ! ... First, temperature and salinity. The depth is given
                ! ... by :
                dz = Zout(ka) - ROMS%state%z_r(ia,ja,kl)
                hh = ROMS%state%z_r(ia,ja,kh) - ROMS%state%z_r(ia,ja,kl)
                frac = dz/hh

                tl = ROMS%state%temp(ia,ja,kl)
                th = ROMS%state%temp(ia,ja,kh)
                ROMS%zfields%temp(ia,ja,ka) = tl + frac*(th-tl)

                sl = ROMS%state%salt(ia,ja,kl)
                sh = ROMS%state%salt(ia,ja,kh)
                ROMS%zfields%salt(ia,ja,ka) = sl + frac*(sh-sl)
                
                u1 = ROMS%state%u(ia  ,ja,kl)
                u2 = ROMS%state%u(ia+1,ja,kl)
                if (u1.eq.ROMS%missing.or.u2.eq.ROMS%missing) then
                  ROMS%zfields%u(ia,ja,ka) = ROMS%missing
                else
                  ul = 0.5_dp * (u1 + u2)
                  u1 = ROMS%state%u(ia  ,ja,kh)
                  u2 = ROMS%state%u(ia+1,ja,kh)
                  uh = 0.5_dp * (u1 + u2)
                  ROMS%zfields%u(ia,ja,ka) = ul + frac*(uh-ul)
                endif
           
                v1 = ROMS%state%v(ia,ja  ,kl)
                v2 = ROMS%state%v(ia,ja+1,kl)
                if (v1.eq.ROMS%missing.or.v2.eq.ROMS%missing) then
                  ROMS%zfields%v(ia,ja,ka) = ROMS%missing
                else
                  vl = 0.5_dp * (v1 + v2)
                  v1 = ROMS%state%v(ia,ja,kh)
                  v2 = ROMS%state%v(ia,ja+1,kh)
                  vh = 0.5_dp * (v1 + v2)
                  ROMS%zfields%v(ia,ja,ka) = vl + frac*(vh-vl)
                endif
              endif 

            endif  ! Shallowness control

          enddo  ! ka-loop

        else   ! Masked rho point

          ROMS%zfields%temp(ia,ja,:) = ROMS%missing
          ROMS%zfields%salt(ia,ja,:) = ROMS%missing
          ROMS%zfields%u(ia,ja,:)    = ROMS%missing
          ROMS%zfields%v(ia,ja,:)    = ROMS%missing

        endif

      enddo   ! ia-loop
    enddo   ! ja-loop
    
  end subroutine roms_zselect
  ! ...
  ! ======================================================================
  ! ...
  subroutine roms_vprofile(ROMS, lon, lat, Zprof_pos, &
                           Tprof, Sprof, Uprof, Vprof, i_rho, j_rho)

    class(type_roms), intent(in)                      :: ROMS
    real(dp), intent(in)                              :: lon, lat
    real(dp), dimension(:), intent(in)                :: Zprof_pos       ! positive meters downward
    real(dp), dimension(:), allocatable, intent(out)  :: Tprof, Sprof, Uprof, Vprof
    integer, intent(out), optional                    :: i_rho, j_rho   ! nearest rho indices (1-based)

    ! ... Local variables
    ! ... 
    integer :: nzr, nzo, ia, ja, ka, kk, kl, kh
    integer :: L, M, Lm, Mm, i0, j0
    real(dp) :: dmin, d, dlon, dlat, lonij, latij
    real(dp), allocatable :: zr(:)              ! z_r vertical for column (nzr)
    real(dp), allocatable :: Zprof(:)
    real(dp) :: tl, th, sl, sh, ul, uh, vl, vh, dz, hh, frac
    real(dp) :: u1, u2, v1, v2
    logical :: masked

    nzr = ROMS%grid%nz_r
    L = ROMS%grid%L; Lm = L - 1
    M = ROMS%grid%M; Mm = M - 1

    if (allocated(Tprof)) deallocate(Tprof)
    if (allocated(Sprof)) deallocate(Sprof)
    if (allocated(Uprof)) deallocate(Uprof)
    if (allocated(Vprof)) deallocate(Vprof)

    nzo = size(Zprof_pos)
    allocate(Tprof(nzo))
    allocate(Sprof(nzo))
    allocate(Uprof(nzo))
    allocate(Vprof(nzo))
    allocate(Zprof(nzo))
    Zprof(:) = -abs(Zprof_pos(:))

    ! ... Set output defaults
    ! ...
    Tprof(:) = ROMS%missing
    Sprof(:) = ROMS%missing
    Uprof(:) = ROMS%missing
    Vprof(:) = ROMS%missing

    ! 1) Find nearest rho point (great-circle is overkill; use Euclidean in lon/lat)
    !    Interior points onlu
    dmin = huge(1.0_dp)
    i0 = -1; j0 = -1
    do ja = 1, Mm
    do ia = 1, Lm
      lonij = ROMS%grid%lon_r(ia, ja)
      latij = ROMS%grid%lat_r(ia, ja)
      dlon = lonij - lon
      dlat = latij - lat
      d = dlon*dlon + dlat*dlat
      if (d < dmin) then
        dmin = d; i0 = ia; j0 = ja
      end if
    enddo
    enddo
    if (present(i_rho)) i_rho = i0
    if (present(j_rho)) j_rho = j0

    ! 2) Check mask if available
    !
    masked = .false.
    if (ROMS%grid%has_mask) then
      if (ROMS%grid%mask_r(i0, j0) <= 0.5_dp) masked = .true.  ! Land
    endif

    if (masked) return

    ! 3) Column z_r and requested depths (negative)
    allocate(zr(nzr))
    zr(:) = ROMS%state%z_r(i0, j0, :)


    ! 4) For each requested depth, do vertical interpolation at rho point

    do ka = 1, nzo

      ! ... Too shallow (above top rho level): assume surface mixed to top rho (nzr)
      ! ...
      if (Zprof(ka) >= zr(nzr)) then
        Tprof(ka) = ROMS%state%temp(i0, j0, nzr)
        Sprof(ka) = ROMS%state%salt(i0, j0, nzr)
        ! U at rho from C-grid U: average adjacent u(ia/ia+1) at same j0
        u1 = ROMS%state%u(i0  , j0, nzr)
        u2 = ROMS%state%u(i0+1, j0, nzr)
        if (u1 == ROMS%missing .or. u2 == ROMS%missing) then
          Uprof(ka) = ROMS%missing
        else
          Uprof(ka) = 0.5_dp * (u1 + u2)
        end if
        ! V at rho from C-grid V: average adjacent v(j0/j0+1) at same i0
        v1 = ROMS%state%v(i0, j0  , nzr)   
        v2 = ROMS%state%v(i0, j0+1, nzr)
        if (v1 == ROMS%missing .or. v2 == ROMS%missing) then
          Vprof(ka) = ROMS%missing
        else
          Vprof(ka) = 0.5_dp * (v1 + v2)
        end if

      else
        ! ... Find bounding levels kl (below) and kh (above) for Zprof(ka)
        kl = 0; kh = 0
        do kk = 1, nzr
          if (zr(kk) > Zprof(ka)) then
            kl = kk - 1
            kh = kk
            exit
          end if
        end do

        if (kl == 0) then
          ! ... Too deep (below bottom): keep missing
          Tprof(ka) = ROMS%missing
          Sprof(ka) = ROMS%missing
          Uprof(ka) = ROMS%missing
          Vprof(ka) = ROMS%missing
        else
          ! ... Vertical interpolation
          dz = Zprof(ka) - zr(kl)
          hh = zr(kh) - zr(kl)
          if (hh == 0.0_dp) then
            frac = 0.0_dp
          else
            frac = dz / hh
          end if

          ! T/S vertical interp at rho
          tl = ROMS%state%temp(i0, j0, kl); th = ROMS%state%temp(i0, j0, kh)
          sl = ROMS%state%salt(i0, j0, kl); sh = ROMS%state%salt(i0, j0, kh)

          if (tl == ROMS%missing .or. th == ROMS%missing) then
            Tprof(ka) = ROMS%missing
            Sprof(ka) = ROMS%missing
          else
            Tprof(ka) = tl + frac * (th - tl)
            Sprof(ka) = sl + frac * (sh - sl)
          end if

          ! ... U at rho: horizontal averaging first at kl,kh then vertical interp
          ! ...
          u1 = ROMS%state%u(i0  , j0, kl)
          u2 = ROMS%state%u(i0+1, j0, kl)
          if (u1 == ROMS%missing .or. u2 == ROMS%missing) then
            Uprof(ka) = ROMS%missing
          else
            ul = 0.5_dp * (u1 + u2)
            u1 = ROMS%state%u(i0  , j0, kh)
            u2 = ROMS%state%u(i0+1, j0, kh)
            if (u1 == ROMS%missing .or. u2 == ROMS%missing) then
              Uprof(ka) = ROMS%missing
            else
              uh = 0.5_dp * (u1 + u2)
              Uprof(ka) = ul + frac * (uh - ul)
            end if
          end if

          ! ... V at rho: horizontal averaging first at kl,kh then vertical interp
          ! ...
          v1 = ROMS%state%v(i0, j0  , kl)
          v2 = ROMS%state%v(i0, j0+1, kl)
          if (v1 == ROMS%missing .or. v2 == ROMS%missing) then
            Vprof(ka) = ROMS%missing
          else
            vl = 0.5_dp * (v1 + v2)
            v1 = ROMS%state%v(i0, j0  , kh)
            v2 = ROMS%state%v(i0, j0+1, kh)
            if (v1 == ROMS%missing .or. v2 == ROMS%missing) then
              Vprof(ka) = ROMS%missing
            else
              vh = 0.5_dp * (v1 + v2)
              Vprof(ka) = vl + frac * (vh - vl)
            end if
          end if

        end if
      end if
    end do

    deallocate(zr, Zprof)
  end subroutine roms_vprofile
  ! ...
  ! ======================================================================
  ! ...
  subroutine roms_zfields_save(ROMS, output_file, record_idx)

    ! ... Write NetCDF file with z-level fields from ROMS%zfields
    ! ... Uses A-grid coordinates from ROMS%grid%lon_r and ROMS%grid%lat_r
    ! ... Time value from ROMS%state%time
    ! ...

    class(type_roms), intent(in)    :: ROMS
    character(len=*), intent(in)    :: output_file
    integer, intent(in)             :: record_idx  ! Time record index to write

    ! ... Local variables
    ! ...
    integer :: ncid_out, err
    integer :: dim_x, dim_y, dim_z, dim_time
    integer :: var_lon, var_lat, var_depth, var_time
    integer :: var_temp, var_salt, var_u, var_v
    integer :: nxa, nya, nzout
    real(dp), allocatable :: lon_a(:,:), lat_a(:,:)
    logical :: file_exists

    ! ... Check if zfields are allocated
    ! ...
    if (.not.allocated(ROMS%zfields%temp)) then
      write(*,*) 'ERROR: ROMS%zfields not allocated. Call zselect first.'
      return
    endif

    nxa = size(ROMS%zfields%temp, 1)
    nya = size(ROMS%zfields%temp, 2)
    nzout = size(ROMS%zfields%depths)

    ! ... Extract A-grid lon/lat (interior points from rho grid)
    ! ...
    allocate(lon_a(nxa, nya))
    allocate(lat_a(nxa, nya))
    lon_a(:,:) = ROMS%grid%lon_r(1:nxa,1:nya)
    lat_a(:,:) = ROMS%grid%lat_r(1:nxa,1:nya)

    ! ... Check if file exists
    ! ...
    inquire(file=trim(output_file), exist=file_exists)

    if (file_exists .and. record_idx > 1) then
      ! ... Open existing file for appending
      ! ...
      err = NF90_OPEN(trim(output_file), NF90_WRITE, ncid_out)
      call nc_error(err, 'Unable to open output file for appending')

      ! ... Get variable IDs
      ! ...
      err = NF90_INQ_VARID(ncid_out, 'time', var_time)
      err = NF90_INQ_VARID(ncid_out, 'temp', var_temp)
      err = NF90_INQ_VARID(ncid_out, 'salt', var_salt)
      err = NF90_INQ_VARID(ncid_out, 'u', var_u)
      err = NF90_INQ_VARID(ncid_out, 'v', var_v)

    else
      ! ... Create new file
      ! ...
      err = NF90_CREATE(trim(output_file), NF90_CLOBBER, ncid_out)
      call nc_error(err, 'Unable to create output NetCDF file')

      ! ... Define dimensions
      ! ...
      err = NF90_DEF_DIM(ncid_out, 'x', nxa, dim_x)
      err = NF90_DEF_DIM(ncid_out, 'y', nya, dim_y)
      err = NF90_DEF_DIM(ncid_out, 'depth', nzout, dim_z)
      err = NF90_DEF_DIM(ncid_out, 'time', NF90_UNLIMITED, dim_time)

      ! ... Define coordinate variables
      ! ...
      err = NF90_DEF_VAR(ncid_out, 'lon', NF90_DOUBLE, [dim_x, dim_y], var_lon)
      err = NF90_PUT_ATT(ncid_out, var_lon, 'long_name', 'longitude')
      err = NF90_PUT_ATT(ncid_out, var_lon, 'units', 'degrees_east')
      err = NF90_PUT_ATT(ncid_out, var_lon, 'axis', 'X')

      err = NF90_DEF_VAR(ncid_out, 'lat', NF90_DOUBLE, [dim_x, dim_y], var_lat)
      err = NF90_PUT_ATT(ncid_out, var_lat, 'long_name', 'latitude')
      err = NF90_PUT_ATT(ncid_out, var_lat, 'units', 'degrees_north')
      err = NF90_PUT_ATT(ncid_out, var_lat, 'axis', 'Y')

      err = NF90_DEF_VAR(ncid_out, 'depth', NF90_DOUBLE, [dim_z], var_depth)
      err = NF90_PUT_ATT(ncid_out, var_depth, 'long_name', 'depth')
      err = NF90_PUT_ATT(ncid_out, var_depth, 'units', 'm')
      err = NF90_PUT_ATT(ncid_out, var_depth, 'positive', 'down')
      err = NF90_PUT_ATT(ncid_out, var_depth, 'axis', 'Z')

      err = NF90_DEF_VAR(ncid_out, 'time', NF90_DOUBLE, [dim_time], var_time)
      err = NF90_PUT_ATT(ncid_out, var_time, 'long_name', 'time')
      err = NF90_PUT_ATT(ncid_out, var_time, 'units', trim(ROMS%time%units))
      err = NF90_PUT_ATT(ncid_out, var_time, 'calendar', trim(ROMS%time%calendar))
      err = NF90_PUT_ATT(ncid_out, var_time, 'axis', 'T')

      ! ... Define data variables
      ! ...
      err = NF90_DEF_VAR(ncid_out, 'temp', NF90_DOUBLE, &
                         [dim_x, dim_y, dim_z, dim_time], var_temp)
      err = NF90_PUT_ATT(ncid_out, var_temp, 'long_name', 'temperature')
      err = NF90_PUT_ATT(ncid_out, var_temp, 'units', 'Celsius')
      err = NF90_PUT_ATT(ncid_out, var_temp, '_FillValue', ROMS%missing)

      err = NF90_DEF_VAR(ncid_out, 'salt', NF90_DOUBLE, &
                         [dim_x, dim_y, dim_z, dim_time], var_salt)
      err = NF90_PUT_ATT(ncid_out, var_salt, 'long_name', 'salinity')
      err = NF90_PUT_ATT(ncid_out, var_salt, 'units', 'PSU')
      err = NF90_PUT_ATT(ncid_out, var_salt, '_FillValue', ROMS%missing)

      err = NF90_DEF_VAR(ncid_out, 'u', NF90_DOUBLE, &
                         [dim_x, dim_y, dim_z, dim_time], var_u)
      err = NF90_PUT_ATT(ncid_out, var_u, 'long_name', 'u-velocity')
      err = NF90_PUT_ATT(ncid_out, var_u, 'units', 'm/s')
      err = NF90_PUT_ATT(ncid_out, var_u, '_FillValue', ROMS%missing)

      err = NF90_DEF_VAR(ncid_out, 'v', NF90_DOUBLE, &
                         [dim_x, dim_y, dim_z, dim_time], var_v)
      err = NF90_PUT_ATT(ncid_out, var_v, 'long_name', 'v-velocity')
      err = NF90_PUT_ATT(ncid_out, var_v, 'units', 'm/s')
      err = NF90_PUT_ATT(ncid_out, var_v, '_FillValue', ROMS%missing)

      ! ... Global attributes
      ! ...
      err = NF90_PUT_ATT(ncid_out, NF90_GLOBAL, 'title', &
                         'ROMS output on z-levels (A-grid)')
      err = NF90_PUT_ATT(ncid_out, NF90_GLOBAL, 'source', &
                         'Converted from ROMS s-coordinates')

      ! ... End define mode
      ! ...
      err = NF90_ENDDEF(ncid_out)
      call nc_error(err, 'Unable to end define mode')

      ! ... Write coordinate variables (only once)
      ! ...
      err = NF90_PUT_VAR(ncid_out, var_lon, lon_a)
      call nc_error(err, 'Unable to write longitude')

      err = NF90_PUT_VAR(ncid_out, var_lat, lat_a)
      call nc_error(err, 'Unable to write latitude')

      err = NF90_PUT_VAR(ncid_out, var_depth, ROMS%zfields%depths)
      call nc_error(err, 'Unable to write depth')

    endif

    ! ... Write time value for this record
    ! ...
    err = NF90_PUT_VAR(ncid_out, var_time, [ROMS%state%time], &
                       start=[record_idx], count=[1])
    call nc_error(err, 'Unable to write time')

    ! ... Write data variables for this record
    ! ...
    err = NF90_PUT_VAR(ncid_out, var_temp, ROMS%zfields%temp, &
                       start=[1,1,1,record_idx], count=[nxa,nya,nzout,1])
    call nc_error(err, 'Unable to write temperature')

    err = NF90_PUT_VAR(ncid_out, var_salt, ROMS%zfields%salt, &
                       start=[1,1,1,record_idx], count=[nxa,nya,nzout,1])
    call nc_error(err, 'Unable to write salinity')

    err = NF90_PUT_VAR(ncid_out, var_u, ROMS%zfields%u, &
                       start=[1,1,1,record_idx], count=[nxa,nya,nzout,1])
    call nc_error(err, 'Unable to write u-velocity')

    err = NF90_PUT_VAR(ncid_out, var_v, ROMS%zfields%v, &
                       start=[1,1,1,record_idx], count=[nxa,nya,nzout,1])
    call nc_error(err, 'Unable to write v-velocity')

    ! ... Close file
    ! ...
    err = NF90_CLOSE(ncid_out)
    call nc_error(err, 'Unable to close output file')

    write(*,*) 'Written record ', record_idx, ' to ', trim(output_file)

    deallocate(lon_a, lat_a)

  end subroutine roms_zfields_save
  ! ...
  ! ======================================================================
  ! ...
end module module_roms

