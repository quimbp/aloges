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

module module_forcing

use netcdf
use module_types
use module_constants
use module_tools
use module_grid 

use module_alm

implicit none

! ... Input files
! ...
character(len=maxlen)                            :: tifilename = ''  ! Traj initial
character(len=maxlen)                            :: OUfilename = ''  ! Oce u vel
character(len=maxlen)                            :: OVfilename = ''  ! Oce v vel
character(len=maxlen)                            :: OWfilename = ''  ! Oce w vel
character(len=maxlen)                            :: OTfilename = ''  ! Oce temp 
character(len=maxlen)                            :: OSfilename = ''  ! Oce sal 
character(len=maxlen)                            :: ORfilename = ''  ! Oce rho
character(len=maxlen)                            :: OCfilename = ''  ! Oce tracer, c
character(len=maxlen)                            :: AUfilename = ''  ! Atm u vel
character(len=maxlen)                            :: AVfilename = ''  ! Atm v vel

! ... Ocean climatological currents
! ...
logical                                          :: OceClim    = .False.

! ... Input file flags:
! ...
logical                                          :: WithOU     = .False.
logical                                          :: WithOV     = .False.
logical                                          :: WithOW     = .False.
logical                                          :: WithOT     = .False.
logical                                          :: WithOS     = .False.
logical                                          :: WithOR     = .False.
logical                                          :: WithOC     = .False.
logical                                          :: WithAU     = .False.
logical                                          :: WithAV     = .False.

! ... Grids of the input physical parameters:
! ...
type(type_ncgrid)                                :: GOU       ! Water U advection
type(type_ncgrid)                                :: GOV       ! Water V advection
type(type_ncgrid)                                :: GOW       ! Water W advection
type(type_ncgrid)                                :: GOT       ! Water Temperature
type(type_ncgrid)                                :: GOS       ! Water Salinity
type(type_ncgrid)                                :: GOR       ! Water Density
type(type_ncgrid)                                :: GOC       ! Water tracer, C
type(type_ncgrid)                                :: GAU       ! Wind U vomponent
type(type_ncgrid)                                :: GAV       ! Wind V vomponent

! ... User specified initial and final time
! ... This information is used to trim the input
! ... forcing files
! ...
logical                                          :: WithTini = .False.
logical                                          :: WithTlen = .False.
character(maxlen)                                :: StrgTini = ''
type(type_date)                                  :: UserDini 
real(dp)                                         :: UserTini
character(maxlen)                                :: StrgTlen = ''
real(dp)                                         :: UserTlen

! ... Time calendar
! ... The model will assume the time units and calendar of the zonal ocean
! ... current's file.
character(len=180)                               :: forcing_time_units = 'seconds since 2000-01-01 00:00:00'
character(len=20)                                :: forcing_time_calendar = 'gregorian'
real(dp)                                         :: Reference_time = 0.0D0
type(type_date)                                  :: Reference_date

! ... Input variables
! ...
character(len=maxlen)                            :: OUxname = '' ! 'longitude'  
character(len=maxlen)                            :: OUyname = '' ! 'latitude'  
character(len=maxlen)                            :: OUzname = '' ! 'depth'  
character(len=maxlen)                            :: OUtname = '' ! 'time'  
character(len=maxlen)                            :: OUvname = 'uo' ! ''  
character(len=maxlen)                            :: OUunits = ''  
character(len=maxlen)                            :: OUcalen = ''  

character(len=maxlen)                            :: OVxname = '' ! 'longitude'  
character(len=maxlen)                            :: OVyname = '' ! 'latitude'  
character(len=maxlen)                            :: OVzname = '' ! 'depth'  
character(len=maxlen)                            :: OVtname = '' ! 'time'  
character(len=maxlen)                            :: OVvname = 'vo' ! 'None'  
character(len=maxlen)                            :: OVunits = ''  
character(len=maxlen)                            :: OVcalen = ''  

character(len=maxlen)                            :: OWxname = '' ! 'longitude'  
character(len=maxlen)                            :: OWyname = '' ! 'latitude'  
character(len=maxlen)                            :: OWzname = '' ! 'depth'  
character(len=maxlen)                            :: OWtname = '' ! 'time'  
character(len=maxlen)                            :: OWvname = '' ! 'None'  
character(len=maxlen)                            :: OWunits = '' ! ''  
character(len=maxlen)                            :: OWcalen = '' ! ''  

character(len=maxlen)                            :: OTxname = ''  
character(len=maxlen)                            :: OTyname = ''  
character(len=maxlen)                            :: OTzname = ''  
character(len=maxlen)                            :: OTtname = ''  
character(len=maxlen)                            :: OTvname = 'thetao'  
character(len=maxlen)                            :: OTunits = ''  
character(len=maxlen)                            :: OTcalen = ''  

character(len=maxlen)                            :: OSxname = ''  
character(len=maxlen)                            :: OSyname = ''  
character(len=maxlen)                            :: OSzname = ''  
character(len=maxlen)                            :: OStname = ''  
character(len=maxlen)                            :: OSvname = 'so'  
character(len=maxlen)                            :: OSunits = ''  
character(len=maxlen)                            :: OScalen = ''  

character(len=maxlen)                            :: ORxname = ''  
character(len=maxlen)                            :: ORyname = ''  
character(len=maxlen)                            :: ORzname = ''  
character(len=maxlen)                            :: ORtname = ''  
character(len=maxlen)                            :: ORvname = ''  
character(len=maxlen)                            :: ORunits = ''  
character(len=maxlen)                            :: ORcalen = ''  

character(len=maxlen)                            :: OCxname = ''  
character(len=maxlen)                            :: OCyname = ''  
character(len=maxlen)                            :: OCzname = ''  
character(len=maxlen)                            :: OCtname = ''  
character(len=maxlen)                            :: OCvname = ''  
character(len=maxlen)                            :: OCunits = ''  
character(len=maxlen)                            :: OCcalen = ''  

character(len=maxlen)                            :: AUxname = ''  
character(len=maxlen)                            :: AUyname = ''  
character(len=maxlen)                            :: AUzname = '-'  
character(len=maxlen)                            :: AUtname = ''  
character(len=maxlen)                            :: AUvname = 'u10'  
character(len=maxlen)                            :: AUunits = ''  
character(len=maxlen)                            :: AUcalen = ''  

character(len=maxlen)                            :: AVxname = ''  
character(len=maxlen)                            :: AVyname = ''  
character(len=maxlen)                            :: AVzname = '-'  
character(len=maxlen)                            :: AVtname = ''  
character(len=maxlen)                            :: AVvname = 'v10'  
character(len=maxlen)                            :: AVunits = ''  
character(len=maxlen)                            :: AVcalen = ''  

integer                                          :: input_id = -1
integer                                          :: input_timeid = 1

! ... Grid bounds
! ...
real(dp)                                         :: forcing_xmin = -1D10
real(dp)                                         :: forcing_ymin = -1D10
real(dp)                                         :: forcing_zmin = -1D10
real(dp)                                         :: forcing_tmin = -1D10
real(dp)                                         :: forcing_xmax =  1D10
real(dp)                                         :: forcing_ymax =  1D10
real(dp)                                         :: forcing_zmax =  1D10
real(dp)                                         :: forcing_tmax =  1D10

contains
! ...
! ====================================================================
! ====================================================================
! ...
  subroutine open_forcing()

    integer i,j,err,nz
    real(dp) xmin,xmax,ymin,ymax,zmin,zmax,tmin,tmax

    GOU%cartesian = Cartesian
    GOV%cartesian = Cartesian
    GOW%cartesian = Cartesian
    GOT%cartesian = Cartesian
    GOS%cartesian = Cartesian
    GOR%cartesian = Cartesian
    GOC%cartesian = Cartesian
    GAU%cartesian = Cartesian
    GAV%cartesian = Cartesian
 
    if (WithOU) then
      if (verb.ge.2) write(*,*) 'Opening ocean U file: ', trim(OUfilename)
      call GOU%open(OUfilename)
    endif
    if (WithOV) then
      if (verb.ge.2) write(*,*) 'Opening ocean V file: ', trim(OVfilename)
      call GOV%open(OVfilename)
    endif
    if (WithOW) then
      if (verb.ge.2) write(*,*) 'Opening ocean W file: ', trim(OWfilename)
      call GOW%open(OWfilename)
    endif
    if (WithOT) then
      if (verb.ge.2) write(*,*) 'Opening ocean T file: ', trim(OTfilename)
      call GOT%open(OTfilename)
    endif
    if (WithOS) then
      if (verb.ge.2) write(*,*) 'Opening ocean S file: ', trim(OSfilename)
      call GOS%open(OSfilename)
    endif
    if (WithOR) then
      if (verb.ge.2) write(*,*) 'Opening ocean Rho file: ', trim(ORfilename)
      call GOR%open(ORfilename)
    endif
    if (WithOC) then
      if (verb.ge.2) write(*,*) 'Opening ocean C file: ', trim(OCfilename)
      call GOC%open(OCfilename)
    endif
    if (WithAU) then
      if (verb.ge.2) write(*,*) 'Opening atmosphere U: ', trim(AUfilename)
      call GAU%open(AUfilename)
    endif
    if (WithAV) then
      if (verb.ge.2) write(*,*) 'Opening atmosphere V: ', trim(AVfilename)
      call GAV%open(AVfilename)
    endif

    if (WithOU) call GOU%scan(OUxname,OUyname,OUzname,OUtname,OUunits,OUcalen)
    if (WithOV) call GOV%scan(OVxname,OVyname,OVzname,OVtname,OVunits,OVcalen)
    if (WithOW) call GOW%scan(OWxname,OWyname,OWzname,OWtname,OWunits,OWcalen)
    if (WithOT) call GOT%scan(OTxname,OTyname,OTzname,OTtname,OTunits,OTcalen)
    if (WithOS) call GOS%scan(OSxname,OSyname,OSzname,OStname,OSunits,OScalen)
    if (WithOR) call GOR%scan(ORxname,ORyname,ORzname,ORtname,ORunits,ORcalen)
    if (WithOC) call GOC%scan(OCxname,OCyname,OCzname,OCtname,OCunits,OCcalen)
    if (WithAU) call GAU%scan(AUxname,AUyname,AUzname,AUtname,AUunits,AUcalen)
    if (WithAV) call GAV%scan(AVxname,AVyname,AVzname,AVtname,AVunits,AVcalen)

    if(GOU%Climatology.and.(GOU%nt.ne.12)) call crash ('OU Monthly climatologies only')
    if(GOV%Climatology.and.(GOV%nt.ne.12)) call crash ('OV Monthly climatologies only')
    if(GOW%Climatology.and.(GOW%nt.ne.12)) call crash ('OW Monthly climatologies only')
    if(GAU%Climatology.and.(GAU%nt.ne.12)) call crash ('AU Monthly climatologies only')
    if(GAV%Climatology.and.(GAV%nt.ne.12)) call crash ('AV Monthly climatologies only')

    ! ... The ALM can use different time units but cannot mix GREGORIAN and NOLEAP
    ! ... calendars for obvious reasons. Then, we check that the calendars are the
    ! ... the same.
    ! ...
    if (WithOU) then
      if (len_trim(forcing_time_calendar).eq.0) then
        forcing_time_calendar = trim(GOU%calendar)
      else
        if (forcing_time_calendar.ne.GOU%calendar) call crash('Incompatible calendars')
      endif
    endif
    if (WithOV) then
      if (len_trim(forcing_time_calendar).eq.0) then
        forcing_time_calendar = trim(GOV%calendar)
      else
        if (forcing_time_calendar.ne.GOV%calendar) call crash('Incompatible calendars')
      endif
    endif
    if (WithOW) then
      if (len_trim(forcing_time_calendar).eq.0) then
        forcing_time_calendar = trim(GOW%calendar)
      else
        if (forcing_time_calendar.ne.GOW%calendar) call crash('Incompatible calendars')
      endif
    endif
    if (WithOT) then
      if (len_trim(forcing_time_calendar).eq.0) then
        forcing_time_calendar = trim(GOT%calendar)
      else
        if (forcing_time_calendar.ne.GOT%calendar) call crash('Incompatible calendars')
      endif
    endif
    if (WithOS) then
      if (len_trim(forcing_time_calendar).eq.0) then
        forcing_time_calendar = trim(GOS%calendar)
      else
        if (forcing_time_calendar.ne.GOS%calendar) call crash('Incompatible calendars')
      endif
    endif
    if (WithOR) then
      if (len_trim(forcing_time_calendar).eq.0) then
        forcing_time_calendar = trim(GOR%calendar)
      else
        if (forcing_time_calendar.ne.GOR%calendar) call crash('Incompatible calendars')
      endif
    endif
    if (WithOC) then
      if (len_trim(forcing_time_calendar).eq.0) then
        forcing_time_calendar = trim(GOC%calendar)
      else
        if (forcing_time_calendar.ne.GOC%calendar) call crash('Incompatible calendars')
      endif
    endif
    if (WithAU) then
      if (len_trim(forcing_time_calendar).eq.0) then
        forcing_time_calendar = trim(GAU%calendar)
      else
        if (forcing_time_calendar.ne.GAU%calendar) call crash('Incompatible calendars')
      endif
    endif
    if (WithAV) then
      if (len_trim(forcing_time_calendar).eq.0) then
        forcing_time_calendar = trim(GAV%calendar)
      else
        if (forcing_time_calendar.ne.GAV%calendar) call crash('Incompatible calendars')
      endif
    endif
    if (len_trim(forcing_time_calendar).eq.0) forcing_time_calendar = standard_calendar
    forcing_time_units = trim(standard_time_units)

    ! ... Check that U,V,T,S,RHO,C have the same vertical structure
    ! ...
    nz = 1
    if (WithOU) nz = max(nz,GOU%nz)
    if (WithOV) nz = max(nz,GOV%nz)
    !if (WithOT) nz = max(nz,GOT%nz)
    !if (WithOS) nz = max(nz,GOS%nz)
    !if (WithOR) nz = max(nz,GOR%nz)
    !if (WithOC) nz = max(nz,GOC%nz)
    if (verb.ge.3) write(*,*) 'Value of nz = ', nz
    if (WithOU.and.GOU%nz.ne.nz) call crash('GOU: Incompatible vertical grid')
    if (WithOV.and.GOV%nz.ne.nz) call crash('GOV: Incompatible vertical grid')
    !if (WithOT.and.GOT%nz.ne.nz) call crash('GOT: Incompatible vertical grid')
    !if (WithOS.and.GOS%nz.ne.nz) call crash('GOS: Incompatible vertical grid')
    !if (WithOR.and.GOR%nz.ne.nz) call crash('GOR: Incompatible vertical grid')
    !if (WithOC.and.GOC%nz.ne.nz) call crash('GOC: Incompatible vertical grid')
  
    ! ... Common geometrical domain:
    ! ...
    xmin = forcing_xmin
    if (WithOU) xmin = max(xmin,GOU%xmin)
    if (WithOV) xmin = max(xmin,GOV%xmin)
    if (WithOW) xmin = max(xmin,GOW%xmin)
    !if (WithOT) xmin = max(xmin,GOT%xmin)
    !if (WithOS) xmin = max(xmin,GOS%xmin)
    !if (WithOR) xmin = max(xmin,GOR%xmin)
    !if (WithOC) xmin = max(xmin,GOC%xmin)
    if (WithAU) xmin = max(xmin,GAU%xmin)
    if (WithAV) xmin = max(xmin,GAV%xmin)
      
    xmax = forcing_xmax
    if (WithOU) xmax = min(xmax,GOU%xmax)
    if (WithOV) xmax = min(xmax,GOV%xmax)
    if (WithOW) xmax = min(xmax,GOW%xmax)
    !if (WithOT) xmax = min(xmax,GOT%xmax)
    !if (WithOS) xmax = min(xmax,GOS%xmax)
    !if (WithOR) xmax = min(xmax,GOR%xmax)
    !if (WithOC) xmax = min(xmax,GOC%xmax)
    if (WithAU) xmax = min(xmax,GAU%xmax)
    if (WithAV) xmax = min(xmax,GAV%xmax)
    
    ymin = forcing_ymin
    if (WithOU) ymin = max(ymin,GOU%ymin)
    if (WithOV) ymin = max(ymin,GOV%ymin)
    if (WithOW) ymin = max(ymin,GOW%ymin)
    !if (WithOT) ymin = max(ymin,GOT%ymin)
    !if (WithOS) ymin = max(ymin,GOS%ymin)
    !if (WithOR) ymin = max(ymin,GOR%ymin)
    !if (WithOC) ymin = max(ymin,GOC%ymin)
    if (WithAU) ymin = max(ymin,GAU%ymin)
    if (WithAV) ymin = max(ymin,GAV%ymin)
    
    ymax = forcing_ymax
    if (WithOU) ymax = min(ymax,GOU%ymax)
    if (WithOV) ymax = min(ymax,GOV%ymax)
    if (WithOW) ymax = min(ymax,GOW%ymax)
    !if (WithOT) ymax = min(ymax,GOT%ymax)
    !if (WithOS) ymax = min(ymax,GOS%ymax)
    !if (WithOR) ymax = min(ymax,GOR%ymax)
    !if (WithOC) ymax = min(ymax,GOC%ymax)
    if (WithAU) ymax = min(ymax,GAU%ymax)
    if (WithAV) ymax = min(ymax,GAV%ymax)
    
    zmin = forcing_zmin
    if (WithOU) zmin = max(zmin,GOU%zmin)
    if (WithOV) zmin = max(zmin,GOV%zmin)
    if (WithOW) zmin = max(zmin,GOW%zmin)
    !if (WithOT) zmin = max(zmin,GOT%zmin)
    !if (WithOS) zmin = max(zmin,GOS%zmin)
    !if (WithOR) zmin = max(zmin,GOR%zmin)
    !if (WithOC) zmin = max(zmin,GOC%zmin)
    
    zmax = 0.0D0
!    if (WithOU) zmax = min(zmax,GOU%zmax)
!    if (WithOV) zmax = min(zmax,GOV%zmax)
!    if (WithOW) zmax = min(zmax,GOW%zmax)
!    !if (WithOT) zmax = min(zmax,GOT%zmax)
!    !if (WithOS) zmax = min(zmax,GOS%zmax)
!    !if (WithOR) zmax = min(zmax,GOR%zmax)
!    !if (WithOC) zmax = min(zmax,GOC%zmax)
   
 
    ! ... Common time domain:
    ! ...
    tmin = forcing_tmin
    if (WithOU.and..not.GOU%Climatology) tmin = max(tmin,GOU%tmin)
    if (WithOV.and..not.GOV%Climatology) tmin = max(tmin,GOV%tmin)
    if (WithOW.and..not.GOW%Climatology) tmin = max(tmin,GOW%tmin)
    if (WithAU.and..not.GAU%Climatology) tmin = max(tmin,GAU%tmin)
    if (WithAV.and..not.GAV%Climatology) tmin = max(tmin,GAV%tmin)

    tmax = forcing_tmax
    if (WithOU.and..not.GOU%Climatology) tmax = min(tmax,GOU%tmax)
    if (WithOV.and..not.GOV%Climatology) tmax = min(tmax,GOV%tmax)
    if (WithOW.and..not.GOW%Climatology) tmax = min(tmax,GOW%tmax)
    if (WithAU.and..not.GAU%Climatology) tmax = min(tmax,GAU%tmax)
    if (WithAV.and..not.GAV%Climatology) tmax = min(tmax,GAV%tmax)
  
    ! ... Model crop
    ! ... T,S,Rho,C do not modify the valid domain, but they are
    ! ... cropped to reduce the memory required by the model
    ! ...
    if (WithOU) call GOU%crop(xmin,xmax,ymin,ymax,tmin,tmax)
    if (WithOV) call GOV%crop(xmin,xmax,ymin,ymax,tmin,tmax)
    if (WithOW) call GOW%crop(xmin,xmax,ymin,ymax,tmin,tmax)
    if (WithOT) call GOT%crop(xmin,xmax,ymin,ymax,tmin,tmax)
    if (WithOS) call GOS%crop(xmin,xmax,ymin,ymax,tmin,tmax)
    if (WithOR) call GOR%crop(xmin,xmax,ymin,ymax,tmin,tmax)
    if (WithOC) call GOC%crop(xmin,xmax,ymin,ymax,tmin,tmax)
    if (WithAU) call GAU%crop(xmin,xmax,ymin,ymax,tmin,tmax)
    if (WithAV) call GAV%crop(xmin,xmax,ymin,ymax,tmin,tmax)

    ! ... Stationarity conditions
    ! ...
    if (WithOU.and.GOU%nt.eq.1) GOU%Stationary = .True.
    if (WithOV.and.GOV%nt.eq.1) GOV%Stationary = .True.
    if (WithOT.and.GOT%nt.eq.1) GOT%Stationary = .True.
    if (WithOS.and.GOS%nt.eq.1) GOS%Stationary = .True.
    if (WithOR.and.GOR%nt.eq.1) GOR%Stationary = .True.
    if (WithOC.and.GOC%nt.eq.1) GOC%Stationary = .True.
    if (WithAU.and.GAU%nt.eq.1) GAU%Stationary = .True.
    if (WithAV.and.GAV%nt.eq.1) GAV%Stationary = .True.

    forcing_xmin = xmin
    forcing_xmax = xmax
    forcing_ymin = ymin
    forcing_ymax = ymax
    forcing_zmin = zmin
    forcing_zmax = zmax
    forcing_tmin = tmin
    forcing_tmax = tmax

    ! ... Now check for the variable of interest:
    ! ...
    if (withOU) then
      err = NF90_INQ_VARID(GOU%fid,trim(OUvname),GOU%varid)
      call cdf_error(err,'Variable '//trim(OUvname)//' in file '//trim(OUfilename))
      GOU%varname = trim(OUvname)
    endif
    if (withOV) then
      err = NF90_INQ_VARID(GOV%fid,trim(OVvname),GOV%varid)
      call cdf_error(err,'Variable '//trim(OVvname)//' in file '//trim(OVfilename))
      GOV%varname = trim(OVvname)
    endif
    if (withOW) then
      err = NF90_INQ_VARID(GOW%fid,trim(OWvname),GOW%varid)
      call cdf_error(err,'Variable '//trim(OWvname)//' in file '//trim(OWfilename))
      GOW%varname = trim(OWvname)
    endif
    if (withOT) then
      err = NF90_INQ_VARID(GOT%fid,trim(OTvname),GOT%varid)
      call cdf_error(err,'Variable '//trim(OTvname)//' in file '//trim(OTfilename))
      GOT%varname = trim(OTvname)
    endif
    if (withOS) then
      err = NF90_INQ_VARID(GOS%fid,trim(OSvname),GOS%varid)
      call cdf_error(err,'Variable '//trim(OSvname)//' in file '//trim(OSfilename))
      GOS%varname = trim(OSvname)
    endif
    if (withOR) then
      err = NF90_INQ_VARID(GOR%fid,trim(ORvname),GOR%varid)
      call cdf_error(err,'Variable '//trim(ORvname)//' in file '//trim(ORfilename))
      GOR%varname = trim(ORvname)
    endif
    if (withOC) then
      err = NF90_INQ_VARID(GOC%fid,trim(OCvname),GOC%varid)
      call cdf_error(err,'Variable '//trim(OCvname)//' in file '//trim(OCfilename))
      GOC%varname = trim(OCvname)
    endif
    if (withAU) then
      err = NF90_INQ_VARID(GAU%fid,trim(AUvname),GAU%varid)
      call cdf_error(err,'Variable '//trim(AUvname)//' in file '//trim(AUfilename))
      GAU%varname = trim(AUvname)
    endif
    if (withAV) then
      err = NF90_INQ_VARID(GAV%fid,trim(AVvname),GAV%varid)
      call cdf_error(err,'Variable '//trim(AVvname)//' in file '//trim(AVfilename))
      GAV%varname = trim(AVvname)
    endif

    if (verb.ge.1) then
      if (withOU) call GOU%show('Ocean zonal currents')
      if (withOV) call GOV%show('Ocean meridional currents')
      if (withOW) call GOW%show('Ocean vertical currents')
      if (withOT) call GOT%show('Ocean temperature')
      if (withOS) call GOS%show('Ocean practical salinity')
      if (withOR) call GOR%show('Ocean water density')
      if (withOC) call GOC%show('Ocean tracer, C')
      if (withAU) call GAU%show('Atmosphere zonal wind')
      if (withAV) call GAV%show('Atmosphere meridional wind')
    endif
    !stop '7777'


  end subroutine open_forcing
  ! ...
  ! ==================================================================
  ! ...
end module module_forcing
