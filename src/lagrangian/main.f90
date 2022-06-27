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

program alm

use module_options
use module_forcing
use module_float
use module_model

implicit none

! ... Initial messaje
! ...
call hello()


! ... Managing the users's command line options
! ...
call options()


! ... Opening forcing fields
! ...
call open_forcing()


! ... Set the initial and final time
! ...
alm_time_units = trim(forcing_time_units)
alm_time_calendar = trim(forcing_time_calendar)

model_time_units = trim(forcing_time_units)
model_time_calendar = trim(forcing_time_calendar)


! ... Set a universal frid and land-sea mask to expedite 
! ... the calculation if a floater is in the free-moving
! ... water or in land
! ...
call define_alm_domain()


! ... Check the user-specified initial and final dates
! ... We calculate alm_tini and alm_tfin based on:
! ...       - The common times in forcing files
! ...       - The optional values provided by the user
! ...       - The forward or reverse integration
! ...
if (WithTini) then
  UserDini = strptime(StrgTini)
  UserDini%calendar = trim(model_time_calendar)
  UserTini = date2num(UserDini,units=model_time_units)
else
  UserTini = 0.0D0
endif

if (reverse) then
  ! ... Backward model
  ! ...
  if (WithTini) then
    alm_tini = min(forcing_tmax,UserTini)
  else
    write(*,*) 'WARNING: No user-provided initial date'
    write(*,*) '         Backward Model initial time = Forcing final time.'
    alm_tini = forcing_tmax
  endif
  if (WithTlen) then
    alm_tfin = max(alm_tini-UserTlen,forcing_tmin)
  else
    write(*,*) 'WARNING: No user-provided simulation length'
    write(*,*) '         Backward Model final time = Forcing initial time.'
    alm_tfin = forcing_tmin
  endif
else
  ! ... Forward model
  ! ...
  if (WithTini) then
    alm_tini = max(forcing_tmin,UserTini)
  else
    write(*,*) 'WARNING: No user-provided initial date'
    write(*,*) '         Forward Model initial time = Forcing initial time.'
    alm_tini = forcing_tmin
  endif
  if (WithTlen) then
    alm_tfin = min(alm_tini+UserTlen,forcing_tmax)
  else
    write(*,*) 'WARNING: No user-provided simulation length'
    write(*,*) '         Forward Model final time = Forcing final time.'
    alm_tfin = forcing_tmax
  endif
endif
alm_tlen = nint(alm_tfin - alm_tini)

model_dini = num2date(model_tini,units=model_time_units, &
                      calendar=model_time_calendar)
model_dfin = num2date(model_tfin,units=model_time_units, &
                      calendar=model_time_calendar)


! ... Initialize model
! ...
call model_ini()

! ... Initialize floats
! ....
call float_ini()
if (Nfloats.eq.0) call crash('No valid floats')


! ... Create output trajectory file
! ...
call trajectory_create(trajectory_name,Nfloats,output_missing)


! ... Run the model 
! ...
call model_run()


! ... Close output files:
! ...
call trajectory_close()

contains
! ...
! ====================================================================
! ...
  subroutine hello

    write(*,*) '========================================================='
    write(*,*) '=                 Aloges Lagrangian Model               ='
    write(*,*) '=                     Version: ',version,'                      ='
    write(*,*) '========================================================='

  end subroutine hello
  ! ...
  ! ====================================================================
  ! ...
  subroutine define_alm_domain()

    ! ... Local variables
    ! ...
    logical withU,withV,withW
    integer i,j,k,nx,ny,nz
    integer IP(3)
    real(dp) fu,fv,fw
    real(dp) XX(3)

    withU = len_trim(GOU%filename).gt.0
    withV = len_trim(GOV%filename).gt.0
    withW = len_trim(GOW%filename).gt.0

    nx = 1
    if (withU) nx = max(nx,GOU%nx)
    if (withV) nx = max(nx,GOV%nx)
    if (withW) nx = max(nx,GOW%nx)

    ny = 1
    if (withU) ny = max(ny,GOU%ny)
    if (withV) ny = max(ny,GOV%ny)
    if (withW) ny = max(ny,GOW%ny)

    nz = 1
    if (withU) nz = max(nz,GOU%nz)
    if (withV) nz = max(nz,GOV%nz)
    if (withW) nz = max(nz,GOW%nz)

    alm_xmin = forcing_xmin
    alm_xmax = forcing_xmax
    alm_ymin = forcing_ymin
    alm_ymax = forcing_ymax
    alm_zmin = forcing_zmin
    alm_zmax = 0.0D0                            ! The topmost value is the free surface
    alm_dx   = (alm_xmax-alm_xmin)/(nx-1.0D0)
    alm_dy   = (alm_ymax-alm_ymin)/(ny-1.0D0)
    alm_nx   = nx
    alm_ny   = ny
    alm_nz   = nz

    allocate(alm_x(nx))
    do i=1,nx
      alm_x(i) = alm_xmin + (i-1)*alm_dx
    enddo

    allocate(alm_y(ny))
    do j=1,ny
      alm_y(j) = alm_ymin + (j-1)*alm_dy
    enddo

    allocate(alm_z(nz))
    if (withU) alm_z(:) = GOU%z(:)
    if (withV) alm_z(:) = GOV%z(:)
    if (withW) alm_z(:) = GOW%z(:)

    allocate(alm_mask(nx,ny,nz))

    do k=1,alm_nz
    do j=1,alm_ny
    do i=1,alm_nx

      XX = [alm_x(i), alm_y(j), alm_z(k)]
  
      if (withU) then
        IP = GOU%locate(XX)
        fu = GOU%interpol(GOU%var(GOU%varid)%mask(GOU%ia:GOU%ib,GOU%ja:GOU%jb,GOU%ka:GOU%kb), &
                          XX,IP,SingleLayer=.True.)
      else
        fu = 1.0D0
      endif

      if (withV) then
        IP = GOV%locate(XX)
        fv = GOV%interpol(GOV%var(GOV%varid)%mask(GOV%ia:GOV%ib,GOV%ja:GOV%jb,GOV%ka:GOV%kb), &
                          XX,IP,SingleLayer=.True.)
      else
        fv = 1.0D0
      endif

      if (withW) then
        IP = GOW%locate(XX)
        fw = GOW%interpol(GOW%var(GOW%varid)%mask(GOW%ia:GOW%ib,GOW%ja:GOW%jb,GOW%ka:GOW%kb), &
                          XX,IP,SingleLayer=.True.)
      else
        fw = 1.0D0
      endif

      if (fu.lt.0.5.or.fv.lt.0.5.or.fw.lt.0.5) then
        alm_mask(i,j,k) = 0.0D0
      else
        alm_mask(i,j,k) = 1.0D0
      endif

    enddo
    enddo
    enddo

    if (verb.ge.4) then
      write(*,*) 'Land-sea mask: 1-Water, 0-Land'
      do k=1,alm_nz
        write(*,*) 'K and depth: ', k, alm_z(k)
        do j=alm_ny,1,-1
          write(*,'(200I1)') nint(alm_mask(:,j,k))
        enddo
        write(*,*)
      enddo
    endif

  end subroutine define_alm_domain
  ! ...
  ! ==================================================================
  ! ...
end program alm
