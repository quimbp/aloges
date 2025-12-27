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
! **************************************************************************
! ... plocan_merge.f90
! ...
! ... Merge a series of PLOCAN, CODAR netcdf files with total velocity components.
! ... Author: Quim Ballabrera
! ... Initial code: 2024, sept 17
! ... Version 0.1: September 2024
! ...  * change run options
! ...  * Add version
! ***************************************************************************

use aloges
use netcdf

implicit none

character(len=5), parameter                           :: Version = 'V0.1'
character(len=20), parameter                          :: Date = 'September, 2024'

integer err,fid,out,dim,var,id,idi,idj,idl,idx,idy
integer ndims,nvars,ngatts,natts,unlimid,vtype,dimids(3)
integer n,nrecs,rec,nx,ny,hours
integer base_year,base_month,base_day,base_hour,year,month,day,hour,base_ord,ord
integer na
character(len=maxlen), dimension(:), allocatable       :: filelist
character(len=maxlen)                                  :: PATH,outfile,word,base_date
real(kind=4), dimension(:), allocatable                :: lon,lat
real(kind=8), dimension(:,:), allocatable              :: f8


! ... Run options
! ...
na = iargc()
if (na.eq.0) then
  write(0,'(T1,A)') 'plocan_merge ERROR: No input files'
  write(*,*)
  stop
endif

! ... If here, at least one option:
! ...
call getarg(1,word)

if (word(1:2).eq.'-v'.or.word(1:3).eq.'--v') then
  write(*,'(T1,4A)') 'PLOCAN_MERGE: ', trim(Version), ' ', trim(Date)
  write(*,'(T1,A)') 'This is free software; see the source for copying conditions.'
  write(*,'(T1,A)') 'There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.'
  write(*,*)
  stop
endif

if (word(1:1).eq.'-') then
  write(*,*)
  write(*,*) 'prompt> plocan_merge list-of-input-files outfile'
  write(*,*) '* Example: plocan_merge CODAR_PLOC_*.nc codar.nc'
  write(*,*)
  stop
endif

! ... If here no, version or help options used.
! ... PATH = './CODAR_PLOC_*.nc'
! ... outfile = 'codar.nc'
! ...
na = iargc()
nrecs = na - 1
if (nrecs.eq.0) then
  write(0,'(T1,A)') 'plocan_merge ERROR: No input files'
  write(*,*)
  stop
endif

allocate(filelist(nrecs))
do rec=1,nrecs
  call getarg(rec,filelist(rec))
enddo
call getarg(na,outfile)

!filelist = ls(PATH)
!nrecs = size(filelist)
!if (nrecs.eq.0) stop 'No valid input files'


err = NF90_CREATE(outfile,NF90_CLOBBER,out)
call nc_error(err,'Unable to create merged file')

do rec=1,nrecs
  print*, trim(filelist(rec))

  err = NF90_OPEN(filelist(rec),NF90_NOWRITE,fid)
  call nc_error(err,'Unable to open input file')

  if (rec.eq.1) then

    err = NF90_INQUIRE(fid,ndims,nvars,ngatts,unlimid)

    ! ... Copy dimensions: Time dimension to be set UNDEFINED
    ! ...
    do dim=1,ndims
      word = ''
      err = NF90_INQUIRE_DIMENSION(fid,dim,name=word,len=n)
      if (word.eq.'lon') nx = n
      if (word.eq.'lat') ny = n
      if (word.eq.'time') n = NF90_UNLIMITED
      err = NF90_DEF_DIM(out,trim(word),n,id)
      call nc_error(err,'Unable to define dimension')
    enddo
    
    ! ... Copy variable properties and attributes
    ! ...
    do var=1,nvars
      word = ''
      err = NF90_INQUIRE_VARIABLE(fid,var,word,vtype,ndims,dimids,natts)
      err = NF90_DEF_VAR(out,word,vtype,dimids(1:ndims),id)
      call nc_error(err,'Unable to define variable')
      call nc_copyatts (.false.,fid,var,out,id,natts)

      if (word.eq.'lon') idx = var
      if (word.eq.'lat') idy = var

      if (word.eq.'time') then
        err = NF90_GET_ATT(fid,var,'base_date',base_date)
        call get_word(base_date,1,word)
        read(word,*) base_year
        call get_word(base_date,2,word)
        read(word,*) base_month
        call get_word(base_date,3,word)
        read(word,*) base_day
        call get_word(base_date,4,word)
        read(word,*) base_hour
      endif
    enddo
    write(*,*) 'BASE_DATE: ', base_year, base_month, base_day, base_hour
    base_ord = ymd2ord(base_year,base_month,base_day)

    ! ... Copy global attributes
    ! ...
    call nc_copyatts (.false.,fid,0,out,0,natts)
    
    ! ... Leave define mode
    ! ...
    err = NF90_ENDDEF(out)
    call nc_error(err,'Unable to leave definition mode')

    allocate(lon(nx)) 
    allocate(lat(ny)) 
    allocate(f8(nx,ny))     ! Double precission field

    err = NF90_GET_VAR(fid,idx,lon)
    err = NF90_PUT_VAR(out,idx,lon)
    call nc_error(err,'Unable to write lon')

    err = NF90_GET_VAR(fid,idy,lat)
    err = NF90_PUT_VAR(out,idy,lat)
    call nc_error(err,'Unable to write lat')

  endif

  do var=1,nvars
    word = ''
    err = NF90_INQUIRE_VARIABLE(fid,var,word,vtype,ndims,dimids,natts)

    if (word.eq.'time') then
      err = NF90_GET_ATT(fid,var,'base_date',base_date)
      call get_word(base_date,1,word)
      read(word,*) year
      call get_word(base_date,2,word)
      read(word,*) month
      call get_word(base_date,3,word)
      read(word,*) day
      call get_word(base_date,4,word)
      read(word,*) hour
      ord = ymd2ord(year,month,day)
      hours = (ord - base_ord)*24 + hour - base_hour
      err = NF90_PUT_VAR(out,var,[hours],[rec],[1])
      call nc_error(err,'Unable to write time')
    else
      ! ... If variable is not time, lon or lat, copy it:
      ! ...
      if (var.ne.idx.and.var.ne.idy) then
        err = NF90_GET_VAR(fid,var,f8)
        err = NF90_PUT_VAR(out,var,f8,[1,1,rec],[nx,ny,1])
        call nc_error(err,'Unable to write double precision variable')
      endif
    endif
  enddo
enddo

err = NF90_CLOSE(out)
call nc_error(err,'Unable to properly close merged file')

stop 'Ok'
end

!dimensions:
!	time = 1 ;
!	lat = 74 ;
!	lon = 70 ;
!variables:
!	int time(time) ;
!		time:long_name = "Valid Time (GMT)" ;
!		time:base_date = "2024, 09, 08, 20" ;
!		time:units = "hours since 2024-09-08 20:00:00" ;
!		time:standard_name = "time" ;
!	float lat(lat) ;
!		lat:long_name = "Latitude" ;
!		lat:units = "degrees north" ;
!		lat:standard_name = "latitude" ;
!	float lon(lon) ;
!		lon:long_name = "Longitude" ;
!		lon:units = "degrees east" ;
!		lon:standard_name = "longitude" ;
!	double u(time, lat, lon) ;
!		u:long_name = "Eastward Water Velocity" ;
!		u:units = "m/s" ;
!		u:standard_name = "eastward_sea_water_velocity" ;
!		u:_FillValue = -9999. ;
!	double v(time, lat, lon) ;
!		v:long_name = "Northward Water Velocity" ;
!		v:units = "m/s" ;
!		v:standard_name = "northward_sea_water_velocity" ;
!		v:_FillValue = -9999. ;
!	double stdu(time, lat, lon) ;
!		stdu:long_name = "standard deviation U" ;
!		stdu:units = "m/s" ;
!		stdu:standard_name = "standard_deviation_u" ;
!		stdu:_FillValue = -9999. ;
!	double stdv(time, lat, lon) ;
!		stdv:long_name = "standard deviation V" ;
!		stdv:units = "m/s" ;
!		stdv:standard_name = "standard_deviation_v" ;
!		stdv:_FillValue = -9999. ;
!	double cov(time, lat, lon) ;
!		cov:long_name = "covariance" ;
!		cov:units = "m^2/s^2" ;
!		cov:standard_name = "covariance" ;
!		cov:_FillValue = -9999. ;
!
!// global attributes:
!		:NC_GLOBAL.netcdf_library_version = "v2" ;
!		:NC_GLOBAL.Conventions = "CF-1.4" ;
!		:NC_GLOBAL.Title = "Near-Real Time Surface Ocean Velocity" ;
!		:NC_GLOBAL.institution = "PLOCAN (Oceanic Platform of the Canary Islands)" ;
!		:NC_GLOBAL.source = "HF Radar Derived Surface Currents obtained from CODAR combine method" ;
!		:NC_GLOBAL.origin = "MAND (measured);TALI (measured);" ;
!		:NC_GLOBAL.history = "09-Sep-2024 00:05:12" ;
!		:NC_GLOBAL.creator_url = "http://www.qualitasinstruments.com/" ;
!		:NC_GLOBAL.institution_references = "http://www.plocan.eu" ;
!		:NC_GLOBAL.contact = "observatory@plocan.eu" ;
!		:NC_GLOBAL.principal_investigator = "Eric Delory" ;
!		:NC_GLOBAL.principal_investigator_mail = "eric.delory@plocan.eu" ;
!		:NC_GLOBAL.sdn_edmo_code = "3497" ;
!		:NC_GLOBAL.distribution_statement = "Data available free of charge. User assumes all risk for use of data. User must display citation in any publication or product using data. User must contact PI prior to any commercial use of data" ;
!		:NC_GLOBAL.citation = "These data were collected and made freely available by PLOCAN (Oceanic Platform of the Canary Island)" ;
!		:NC_GLOBAL.geospatial_lat_min = "27.8505" ;
!		:NC_GLOBAL.geospatial_lat_max = "28.3446" ;
!		:NC_GLOBAL.geospatial_lon_min = "-15.4203" ;
!		:NC_GLOBAL.geospatial_lon_max = "-14.8937" ;
!		:NC_GLOBAL.grid_resolution = "0.75km" ;
!		:NC_GLOBAL.grid_projection = "equidistal cylindrical" ;
!		:NC_GLOBAL.grid_type = "REGULAR" ;
!		:NC_GLOBAL.UUID = "04848e6b-2be9-41d9-839c-bff297accb58" ;
!}
!
