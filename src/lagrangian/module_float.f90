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

module module_float

use module_types
use module_constants
use module_tools
use module_time

use module_alm

implicit none

type type_float
  logical                                        :: released = .False.
  logical                                        :: floating = .False.
  logical                                        :: indomain = .False.
  logical                                        :: surface  = .False.
  real(dp)                                       :: xo,yo,zo,to   ! position and time initial
  real(dp)                                       :: x,y,z,t       ! position and time
  real(dp)                                       :: u = 0.0D0     ! x-velocity (m/s)
  real(dp)                                       :: v = 0.0D0     ! y-velocity (m/s)
  real(dp)                                       :: w = 0.0D0     ! z-velocity (m/s)
  real(dp)                                       :: speed = 0.0D0 ! velocity module (m/s)
  real(dp)                                       :: tlast = 0.0D0 ! Last floating time  (s)
  real(dp)                                       :: dist  = 0.0D0  ! travelled distance (km)
  real(dp)                                       :: Pdens = 0.0D0  ! Particle density 
  real(dp)                                       :: Psize = 0.0D0  ! Particle size
  real(dp)                                       :: Pbuoy = 0.0D0  ! Particle buoyancy
  real(dp)                                       :: Wtemp = 0.0D0  ! Water temperature
  real(dp)                                       :: Wpsal = 0.0D0  ! Water salinity
  real(dp)                                       :: Wdens  = 0.0D0 ! Water density
  real(dp)                                       :: Wtrac = 0.0D0  ! Water tracer conc.
end type type_float

integer                                          :: Nfloats = 0
type(type_float), dimension(:), pointer          :: FLT      ! Array of buoys
type(type_float)                                 :: FLTk     ! A buoy

logical                                          :: Release_by_file = .False.
logical                                          :: Release_by_pos  = .False.
logical                                          :: WithReleaseXo   = .False.
logical                                          :: WithReleaseYo   = .False.
logical                                          :: WithReleaseZo   = .False.
logical                                          :: WithReleaseTime = .False.
logical                                          :: WithReleaseRho  = .False.
logical                                          :: WithReleaseSize = .False.
real(dp)                                         :: Release_xo
real(dp)                                         :: Release_yo
real(dp)                                         :: Release_zo      = 0.0D0
real(dp)                                         :: Release_to      = 0.0D0
real(dp)                                         :: Release_rho     = 1027.0D0
real(dp)                                         :: Release_size    = 1D-4 ! 2r
character(len=maxlen)                            :: Release_file    = ''
character(len=maxlen)                            :: Release_time    = ''

! ... Release of a random cloud of particles:
! ...
logical                                          :: WithRandom = .False.
integer                                          :: Nrandom    = 1
real(dp)                                         :: Radius_x   = 0.02D0 ! deg
real(dp)                                         :: Radius_y   = 0.02D0 ! deg
real(dp)                                         :: Radius_z   = 0.0    ! m
real(dp)                                         :: Radius_t   = 0.0    ! s
character(len=maxlen)                            :: Release_SaveFile= '' ! To save the resulting data

contains

  ! ...
  ! ==================================================================
  ! ...
  subroutine float_ini

    type(type_date) release_date
    logical withdate,valid,is_land,wrngF,wrngR
    integer i,j,k,l,flo,iu
    real(dp) xmin,xmax,ymin,ymax,zmin,zmax,tmin,tmax
    real(dp) xxr,yyr,zzr,ttr
    real(dp) rnd(4)
    character(len=3) asterix

    if (WithReleaseZo) then
      Release_zo = -abs(Release_zo)
    else
      Release_zo = 0.0D0
    endif

    if (WithReleaseTime) then
      ! ... Check how the user has defined the release time
      ! ... Either it has entered a value in seconds
      ! ... or a date. If a date, we process the date and 
      ! ... substract the initial date.
      ! ...
      if (is_numeric(Release_time)) then
        ! ... The value entered is a number:
        ! ...
        read(Release_time,*) Release_to
      else
        ! ... The value entered is a string with a date
        ! ...
        release_date = strptime(Release_time)
        release_date%calendar = alm_time_calendar
        Release_to = anint(date2num(release_date,units=alm_time_units) - alm_tini)
      endif
    else
      Release_to = 0.0D0
    endif
 

    if (Release_by_pos) then
      ! ... We have provided xo and yo coordinates

      Nfloats = Nrandom
      allocate (FLT(Nfloats))

      FLT(1)%xo = deg2rad*Release_xo
      FLT(1)%yo = deg2rad*Release_yo
      FLT(1)%zo = Release_zo
      FLT(1)%to = Release_to

      if (Nfloats.gt.1) then
        do flo=2,Nfloats
          is_land = .True.
          do while (is_land) 
            call random_number(rnd)
            rnd = 2.0D0*(rnd-0.5D0)
            xxr = deg2rad*(Release_xo + Radius_x*rnd(1))
            yyr = deg2rad*(Release_yo + Radius_y*rnd(2))
            zzr = Release_zo + Radius_z*rnd(3) 
            ttr = Release_to + Radius_t*rnd(4)
            i = point_type(xxr,yyr,zzr)
            if (i.eq.1) then 
              is_land = .False.
            else
              !print*, 'unsuitable ', xxr, yyr, zzr
            endif
          enddo
          FLT(flo)%xo = xxr
          FLT(flo)%yo = yyr
          FLT(flo)%zo = zzr
          FLT(flo)%to = ttr
        enddo
      endif

    else if (Release_by_file) then

      Nfloats = release_read(Release_file)
      if (Nfloats.eq.0) return

    else

      ! ... If here, it means that the user just wants a number of floats
      ! ... randomly placed inside the domain.
      ! ...      
      Nfloats = Nrandom
      allocate (FLT(Nfloats))
      xmin = alm_xmin; xmax = alm_xmax
      ymin = alm_ymin; ymax = alm_ymax
      zmin = alm_zmin; zmax = alm_zmax
      tmin = alm_tmin; tmax = alm_tmax

      do flo=1,Nfloats
        is_land = .True.
        do while (is_land) 
          call random_number(rnd)
          xxr = xmin + (xmax-xmin)*rnd(1)
          yyr = ymin + (ymax-ymin)*rnd(2)
          if (Radius_z.gt.0.0D0) then
            zzr = zmin + (zmax-zmin)*rnd(3)
          else
            zzr = release_zo
          endif
          if (Radius_t.gt.0.0D0) then
            ttr = tmin + (tmax-tmin)*rnd(4)
          else
            ttr = release_to
          endif
          i = point_type(xxr,yyr,zzr)
          if (i.eq.1) then 
            is_land = .False.
          else
            !print*, 'unsuitable ', xxr, yyr, zzr
          endif
        enddo
        FLT(flo)%xo = xxr
        FLT(flo)%yo = yyr
        FLT(flo)%zo = zzr
        FLT(flo)%to = ttr
      enddo

    endif 

    if (verb.ge.1) then
      wrngF = .False.
      wrngR = .False.
      if (reverse) then 
        if (any(FLT(:)%to.gt.0)) wrngR = .True.
      else
        if (any(FLT(:)%to.lt.0)) wrngF = .True.
      endif
      write(*,*)
      write(*,*) 'Particles to be released'
      if (wrngF) write(*,*) '(*) Floats with negative release time will be released at initial time'
      if (wrngR) write(*,*) '(*) Floats with positive release time will be released at initial time'
      write(*,*) '   lon      lat     depth         date             secs since initial'
      write(*,*) '====================================================================='
      do i=1,Nfloats
        asterix = ' '
        if (reverse) then
          if (FLT(i)%to.gt.0.0D0) asterix = '(*)'
        else
          if (FLT(i)%to.lt.0.0D0) asterix = '(*)'
        endif
        release_date = num2date(FLT(i)%to+alm_tini,alm_time_units,alm_time_calendar)
        write(*,'(F9.3,F9.3,F7.1,3X,A,F9.0,X,A3)') rad2deg*FLT(i)%xo, rad2deg*FLT(i)%yo, &
                                              FLT(i)%zo, release_date%iso(), FLT(i)%to,  &
                                              asterix
      enddo
      write(*,*)
      write(*,*) 'Number of particles: ', Nfloats
      write(*,*)
    endif

    ! ... Ensure that floats will be released at initial time if release to < 0:
    ! ... Also include the buoyancy-related values:
    ! ...
    do i=1,Nfloats
      if (reverse) then
        FLT(i)%to = min(FLT(i)%to,0.0D0)
      else
        FLT(i)%to = max(FLT(i)%to,0.0D0)
      endif
      FLT(i)%Pdens = Release_rho
      FLT(i)%Psize = Release_size
    enddo

    if (len_trim(Release_SaveFile).gt.0) then
      iu = unitfree()
      open(iu,file=Release_SaveFile,status='unknown')
      do i=1,Nfloats
        write(iu,'(F9.4,2X,F9.4,2X,F6.1,2X,A)') rad2deg*FLT(i)%xo, &
                                                rad2deg*FLT(i)%yo, &
                                                FLT(i)%zo, &
                                                release_date%iso()
      enddo
      close(iu)
    endif
    !stop '99999'

 
  end subroutine float_ini
  ! ...
  ! ==================================================================
  ! ...
  integer function release_read(filename) result(n)

   character(len=*), intent(in)                  :: filename

   ! ... Local variables
   ! ... 
   character(maxlen) path,basename,extension

   call filename_split(filename,path,basename,extension)
   extension = uppercase(extension)

   if (index(extension,'SON').gt.0) then
     ! ... Read a JSON file
     ! ...
     print*, 'to read a JSON file'

   else if (index(extension,'NC')+index(extension,'CDF').gt.0) then
     ! ... Read a NetCDF file
     ! ...
     print*, 'to read a NetCDF file'
   else
     ! ... Read an ASCII file
     ! ...
     n = release_read_ascii(filename)
   endif

  end function release_read
  ! ...
  ! ==================================================================
  ! ...
  integer function release_read_ascii (filename) result(n)

    character(len=*), intent(in)                  :: filename

    ! ... Local variables
    ! ...
    logical withdate,valid,is_land
    integer iu,i,ii,j,nheader,nlines,nmax,flo,err
    real(dp) rnd(4)
    real(dp) xx,yy,zz,tt
    real(dp) xo,yo,zo,to
    real(dp) xxr,yyr,zzr,ttr
    character(len=maxlen) line,str

    real(dp), dimension(:), allocatable         :: rlon
    real(dp), dimension(:), allocatable         :: rlat
    real(dp), dimension(:), allocatable         :: rdepth
    real(dp), dimension(:), allocatable         :: rtime
    type(type_date)                             :: rdate

    if (verb.ge.2) write(*,*) 'Reading release ASCII file: ', trim(filename)

    iu = unitfree()
    open(iu,file=filename,status='old',iostat=err)
    if (err.ne.0) call crash('Input release file not found')

    nlines = numlines(iu)

    ! ... Check for header lines
    ! ...
    nheader = 0
    do i=1,nlines
      read(iu,'(A)') line
      if (line(1:1).eq.'#') nheader = nheader + 1
    enddo
    nlines  = nlines - nheader
    nmax = nlines * Nrandom

    allocate(rlon(nmax))
    allocate(rlat(nmax))
    allocate(rdepth(nmax))
    allocate(rtime(nmax))

    !allocate(FLT(nmax))

    ! ... Check with format has been used
    ! ...
    rewind(iu)
    do i=1,nheader
      read(iu,*) 
    enddo

    read(iu,'(A)') line   ! First release position check if there is a date
    i = index(line,'T')
    if (i.gt.0) then
      withdate = .True.
    else
      withdate = .False.
    endif

    ! ... Now read the file
    ! ...
    rewind(iu)
    do i=1,nheader
      read(iu,*)
    enddo

    ii = 0
    do i=1,nlines
      if (withdate) then
        read(iu,*) xx, yy, zz, str
      else
        read(iu,*) xx, yy, zz, tt
      endif

      xx = deg2rad*xx; yy = deg2rad*yy; zz = -abs(zz)
      valid = point_type(xx,yy,zz).eq.1
      if (valid) then
        ii = ii + 1
        rlon(ii) = xx
        rlat(ii) = yy
        rdepth(ii) = zz

        if (withdate) then
          rdate = strptime(str)
          rdate%calendar = alm_time_calendar
          rtime(ii) = anint(date2num(rdate,units=alm_time_units) - alm_tini)
        else
          rtime(ii) = tt
        endif

        ! ... Now check the random points:
        ! ...
        xo = rlon(ii)   ! Already in radians
        yo = rlat(ii)   ! Already in radians
        zo = rdepth(ii)
        to = rtime(ii)

        do flo=2,NRandom
          is_land = .True.
          do while (is_land) 
            call random_number(rnd)
            rnd = 2.0D0*(rnd-0.5D0)
            xxr = xo + deg2rad*Radius_x*rnd(1)
            yyr = yo + deg2rad*Radius_y*rnd(2)
            zzr = zo + Radius_z*rnd(3) 
            ttr = to + Radius_t*rnd(4)
            is_land = point_type(xxr,yyr,zzr).ne.1
          enddo
          ii = ii + 1
          rlon(ii)   = xxr
          rlat(ii)   = yyr
          rdepth(ii) = zzr
          rtime(ii)  = ttr
        enddo
      else
        if (verb.ge.2) write(*,*) 'WARNING: release position not retained ', xx,yy,zz
      endif
    enddo

    n = ii
    if (n.gt.0) then
      allocate(FLT(n))
      do flo=1,n
        FLT(flo)%xo = rlon(flo)
        FLT(flo)%yo = rlat(flo)
        FLT(flo)%zo = rdepth(flo)
        FLT(flo)%to = rtime(flo)
      enddo
    endif

    deallocate(rlon)
    deallocate(rlat)
    deallocate(rdepth)
    deallocate(rtime)

  end function release_read_ascii
  ! ...
  ! ==================================================================
  ! ...
end module module_float

