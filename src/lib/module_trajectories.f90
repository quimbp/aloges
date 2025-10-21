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
! - randn                                                                  !
! -------------------------------------------------------------------------!

module module_trajectories

use netcdf
use module_types
use module_constants
use module_tools
use module_time
use module_nc

implicit none (type, external)

type type_trajectory
  integer                                        :: N         ! Number of trajectories
  integer                                        :: nsteps    ! Number of time steps in trajectories
  character(len=maxlen)                          :: filename
  real(dp), dimension(:,:), pointer              :: x,y,z
  real(dp), dimension(:), pointer                :: t
  type(type_date), dimension(:), pointer         :: date
  real(dp), dimension(:,:), pointer              :: temp,psal,dens

  contains
    procedure                   :: read          => traj_read
    procedure                   :: flip          => traj_flip

end type type_trajectory

contains
! ...
! =====================================================================
! =====================================================================
! ...
  subroutine traj_read(TR,filename,verb)

    class(type_trajectory), intent(inout)        :: TR
    character(len=*), intent(in)                 :: filename
    logical, intent(in)                          :: verb

    ! ... Local variables
    ! ...
    character(len=maxlen)                        :: fpath,fbase,ftype

    call filename_split (filename,fpath,fbase,ftype)
    ftype = lowercase(ftype)

    if (index(ftype,'dat').gt.0) then
      if (verb) write(*,*) 'Reading ASCII file : ', trim(filename)
      !call trajectory_read_ascii (TR,filename)
      stop 'Not coded'

    else if ((ftype.eq.'nc').or.(ftype.eq.'cdf')) then
      if (verb) write(*,*) 'Reading NetCDF file : ', trim(filename)
      call traj_read_nc (TR,filename)

    else if (index(ftype,'json').gt.0) then
      if (verb) write(*,*) 'Reading GEOJSON file : ', trim(filename)
      stop 'Not coded'

    else
      call crash('Input trajectory file: unknown format')
    endif

    if (verb) then                       
     write(*,*) '------------------- '
     write(*,*) 'Number particles : ', TR%N
     write(*,*) 'Number steps     : ', TR%Nsteps
    endif

  end subroutine traj_read
  ! ...
  ! ===================================================================
  ! ...
  subroutine traj_read_nc(TR,filename)

    class(type_trajectory), intent(inout)        :: TR
    character(len=*), intent(in)                 :: filename

    ! ... Local variables
    ! ...
    integer fid,idx,idy,idz,idt,err
    integer np,nt
    character(len=80) dname

    err = NF90_OPEN(filename,0,fid)
    call nc_error(err,'TRAJ_READ: error opening '//trim(filename))

    err = NF90_INQUIRE_DIMENSION(fid,1,name=dname,len=np)
    call nc_error(err,'TRAJ_READ: error reading first dimension')
    err = NF90_INQUIRE_DIMENSION(fid,2,name=dname,len=nt)
    call nc_error(err,'TRAJ_READ: error reading second dimension')

    err = NF90_INQ_VARID(fid,'lon',idx)
    call nc_error(err,'TRAJ_READ: error inquiring about lon')
    err = NF90_INQ_VARID(fid,'lat',idy)
    call nc_error(err,'TRAJ_READ: error inquiring about lat')
    err = NF90_INQ_VARID(fid,'depth',idz)
    call nc_error(err,'TRAJ_READ: error inquiring about depth')
    err = NF90_INQ_VARID(fid,'time',idt)
    call nc_error(err,'TRAJ_READ: error inquiring about time')

    TR%N = np
    TR%nsteps = nt
    TR%filename = trim(filename)

    allocate(TR%t(nt))
    allocate(TR%x(np,nt))
    allocate(TR%y(np,nt))
    allocate(TR%z(np,nt))

    err = NF90_GET_VAR(fid,idt,TR%t)
    call nc_error(err,'TRAJ_READ: error reading time')
    err = NF90_GET_VAR(fid,idx,TR%x)
    call nc_error(err,'TRAJ_READ: error reading lon')
    err = NF90_GET_VAR(fid,idy,TR%y)
    call nc_error(err,'TRAJ_READ: error reading lat')
    err = NF90_GET_VAR(fid,idz,TR%z)
    call nc_error(err,'TRAJ_READ: error reading depth')

    err = NF90_CLOSE(fid)

  end subroutine traj_read_nc
  ! ...
  ! ===================================================================
  ! ...
  subroutine traj_flip(TR)

    class(type_trajectory), intent(inout)        :: TR

    ! ... Local variables
    ! ...
    integer j,n
    real(dp), dimension(TR%Nsteps)       :: tmp1
    real(dp), dimension(TR%N,TR%Nsteps)  :: tmp2

    n = TR%Nsteps

    ! ... Revert longitudes
    ! ...
    tmp2 = TR%x
    do j=1,n
      TR%x(:,j) = tmp2(:,n-j+1)
    enddo

    ! ... Revert longitudes
    ! ...
    tmp2 = TR%y
    do j=1,n
      TR%y(:,j) = tmp2(:,n-j+1)
    enddo

    ! ... Revert depths
    ! ...
    tmp2 = TR%z
    do j=1,n
      TR%z(:,j) = tmp2(:,n-j+1)
    enddo

    ! ... Revert time
    ! ...
    tmp1 = TR%t
    do j=1,n
      TR%t(j) = tmp1(n-j+1)
    enddo


end subroutine traj_flip
  ! ...
  ! ===================================================================
  ! ...
  ! ...
  ! ===================================================================
  ! ...
end module module_trajectories
