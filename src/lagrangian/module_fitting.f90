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

module module_fitting

use module_types
use module_constants
use module_math
use module_time
use module_grid
use module_minpack

use module_alm
use module_forcing
use module_float
use module_model

implicit none

! ... The fitting option is defined in module_alm.f90:  option_fitting = .False. (default)
! ...
type type_trajectory
  character(len=maxlen)                          :: filename = ''
  integer                                        :: Np = 0
  integer                                        :: Nt = 0
  real(dp)                                       :: dt = 0.0D0
  logical, dimension(:,:),  pointer              :: valid
  type(type_date), dimension(:,:), pointer       :: date
  real(dp), dimension(:,:), pointer              :: x
  real(dp), dimension(:,:), pointer              :: y
  real(dp), dimension(:,:), pointer              :: z
  real(dp), dimension(:,:), pointer              :: t
  real(dp), dimension(:,:), pointer              :: uf       ! File (if stored) u value
  real(dp), dimension(:,:), pointer              :: vf       ! File (if stored) v value
  real(dp), dimension(:,:), pointer              :: wf       ! File (if stored) w value
  real(dp), dimension(:,:), pointer              :: uo       ! Derived value
  real(dp), dimension(:,:), pointer              :: vo       ! Derived value
  real(dp), dimension(:,:), pointer              :: wo       ! Derived value
  real(dp), dimension(:,:), pointer              :: speed    ! Derived value
  real(dp), dimension(:,:), pointer              :: OU       ! Interpolated OU value
  real(dp), dimension(:,:), pointer              :: OV       ! Interpolated OV value
  real(dp), dimension(:,:), pointer              :: OW       ! Interpolated OW value
  real(dp), dimension(:,:), pointer              :: AU       ! Interpolated AU value
  real(dp), dimension(:,:), pointer              :: AV       ! Interpolated AV value

  contains
    procedure                   :: read          => trajectory_read
    procedure                   :: allocate      => trajectory_allocate
    procedure                   :: velocity      => trajectory_velocity
    procedure                   :: compress      => trajectory_compress

end type type_trajectory

real(dp)                                         :: fit_dop1  = 1.0D0
real(dp)                                         :: fit_dop2  = 1.0D0
real(dp)                                         :: fit_dop3  = 1.0D0
real(dp)                                         :: fit_dop4  = 1.0D0
real(dp)                                         :: fit_dop5  = 1.0D0

real(dp)                                         :: fit_vmin  =  0.0D0
real(dp)                                         :: fit_vmax  = 10.0D0
real(dp)                                         :: fit_FGp1  = 1.0D0
real(dp)                                         :: fit_FGp2  = 1.0D0
real(dp)                                         :: fit_FGp3  = 1.0D0
real(dp)                                         :: fit_FGp4  = 1.0D0
real(dp)                                         :: fit_FGp5  = 1.0D0

integer                                          :: fit_Maxfev = 200
integer                                          :: fit_Mode   = 1      ! Self-scaling of variables
integer                                          :: fit_Nprint = 0
real(dp)                                         :: fit_Ftol = sqrt(epsilon(1.0_dp))
real(dp)                                         :: fit_Xtol = sqrt(epsilon(1.0_dp))
real(dp)                                         :: fit_Gtol = 0.0D0
real(dp)                                         :: fit_Factor = 100.0D0
character(len=maxlen)                            :: fit_fout  = 'fitting.dat'
character(len=maxlen)                            :: fit_diag  = 'fitting.diag'


type(type_trajectory), dimension(:), allocatable :: To

contains
! ...
! =====================================================================
! =====================================================================
! ...
  subroutine model_fitting()

    ! ... Local variables
    ! ...
    logical                                          :: firstmin,firstmax
    integer                                          :: i,j,Nf,icount,nvalid,iu
    integer                                          :: pp,ll
    real(dp)                                         :: dtmin,dtmax,dtmean
    real(dp)                                         :: um,vm
    real(dp), dimension(ndims)                       :: xo
    character(maxlen), dimension(:), allocatable     :: filelist

    ! ... Minimization variables
    ! ...
    integer, parameter                               :: Np = 5
    integer                                          :: fit_nobs
    integer                                          :: Ninfo
    integer                                          :: NFev
    integer                                          :: NJev
    integer, dimension(Np)                           :: ipvt
    real(dp), dimension(:), allocatable              :: Fvec
    real(dp), dimension(:,:), allocatable            :: FJac
    real(dp)                                         :: Jcosto,Jcost
    real(dp), dimension(Np)                          :: p,po,pn,diag,qtf
    real(dp), dimension(Np)                          :: dJcost

  
    ! ... First of all, allocate space for forcing fields
    ! ...
    allocate (OUrhs(GOU%nx,GOU%ny,GOU%nz,2))
    allocate (OVrhs(GOV%nx,GOV%ny,GOV%nz,2))
    allocate (OWrhs(GOW%nx,GOW%ny,GOW%nz,2))
    allocate (AUrhs(GAU%nx,GAU%ny,GAU%nz,2))
    allocate (AVrhs(GAV%nx,GAV%ny,GAV%nz,2))
    OUrhs = 0.0D0
    OVrhs = 0.0D0
    OWrhs = 0.0D0
    AUrhs = 0.0D0
    AVrhs = 0.0D0

    ! ... Check out how many trajectories:
    ! ... 
    Nf = numwords(trajectory_name)
    allocate(filelist(Nf))
    do i=1,Nf
      call line_word(trajectory_name,i,filelist(i))
    enddo
    allocate (To(Nf))
    
    ! ... Read and process trajectories
    ! ...
    do i=1,Nf

      if (verb.ge.1) write(*,*) 'Reading file: ', trim(filelist(i))

      ! ... Read (one file at a time):
      ! ...
      call To(i)%read(filelist(i))

      ! ... Pass from deg to radians
      ! ...
      do pp=1,To(i)%Np
      do ll=1,To(i)%Nt
        if (To(i)%valid(pp,ll)) then
          To(i)%x(pp,ll) = deg2rad*To(i)%x(pp,ll)
          To(i)%y(pp,ll) = deg2rad*To(i)%y(pp,ll)
          To(i)%z(pp,ll) =    -abs(To(i)%z(pp,ll))
        endif
      enddo
      enddo

      ! ... Calculate the velocities
      ! ... It calculates the dt between successive valid positions
      ! ...
      call To(i)%velocity()

      if (verb.ge.5) then
        write(*,*)
        write(*,*) 'Estimated velocities: '
        do pp=1,To(i)%Np
        do ll=1,To(i)%Nt-1
          write(*,'(2I6,5F9.4,X,L)') pp, ll, To(i)%x(pp,ll), To(i)%y(pp,ll), To(i)%uo(pp,ll), &
                                            To(i)%vo(pp,ll), To(i)%wo(pp,ll), To(i)%valid(pp,ll)
        enddo
        enddo
      endif

      ! ... Reset the forcing field counters
      ! ...
      GOU%rec1 = -1; GOU%rec2 = -1
      GOV%rec1 = -1; GOV%rec2 = -1
      GOW%rec1 = -1; GOW%rec2 = -1
      GAU%rec1 = -1; GAU%rec2 = -1
      GAV%rec1 = -1; GAV%rec2 = -1

      do ll=1,To(i)%Nt-1
      do pp=1,To(i)%Np
        if (To(i)%valid(pp,ll)) then
          xo = [To(i)%x(pp,ll), To(i)%y(pp,ll), To(i)%z(pp,ll) ]

          To(i)%OU(pp,ll) = 0.0D0; To(i)%OV(pp,ll) = 0.0D0; To(i)%OW(pp,ll) = 0.0D0
          To(i)%AU(pp,ll) = 0.0D0; To(i)%AV(pp,ll) = 0.0D0
          call get_forcing('OU',GOU,To(i)%t(pp,ll),To(i)%date(pp,ll),0.0D0,OUrhs)
          To(i)%OU(pp,ll) = time_interpol(GOU,OUrhs,To(i)%t(pp,ll),xo)
          call get_forcing('OV',GOV,To(i)%t(pp,ll),To(i)%date(pp,ll),0.0D0,OVrhs)
          To(i)%OV(pp,ll) = time_interpol(GOV,OVrhs,To(i)%t(pp,ll),xo)

          if (WithOW) then
            call get_forcing('OW',GOW,To(i)%t(pp,ll),To(i)%date(pp,ll),0.0D0,OWrhs)
            To(i)%OW(pp,ll) = time_interpol(GOW,OWrhs,To(i)%t(pp,ll),xo)
          else
            To(i)%OW(pp,ll) = 0.0D0
          endif
          if (winds) then
            call get_forcing('AU',GAU,To(i)%t(pp,ll),To(i)%date(pp,ll),0.0D0,AUrhs)
            To(i)%AU(pp,ll) = time_interpol(GAU,AUrhs,To(i)%t(pp,ll),xo)
            call get_forcing('AV',GAV,To(i)%t(pp,ll),To(i)%date(pp,ll),0.0D0,AVrhs)
            To(i)%AV(pp,ll) = time_interpol(GAV,AVrhs,To(i)%t(pp,ll),xo)
          else
            To(i)%AU(pp,ll) = 0.0D0
            To(i)%AV(pp,ll) = 0.0D0
          endif
          if (verb.ge.4) then
            write(*,*) 'Particle and step  : ', pp, ll
            write(*,*) 'MODEL_FIT Time, xo : ', To(i)%t(pp,ll), xo
            write(*,*) '      OU, OV, OW   : ', To(i)%OU(pp,ll), To(i)%OV(pp,ll), To(i)%OW(pp,ll)
            write(*,*) '      AU, AV       : ', To(i)%AU(pp,ll), To(i)%AV(pp,ll)
          endif
        endif
      enddo
      enddo

    enddo

    ! ... At this moment we have everything for the evaluation of the cost function
    ! ... The information is on the To structure.

    
    nvalid = 0
    do i=1,Nf
    do pp=1,To(i)%Np
      nvalid = nvalid + count(To(i)%valid(pp,1:To(i)%Nt-1))
    enddo
    enddo

    if (verb.ge.1) then
      write(*,*) 
      write(*,*) '--------------'
      write(*,*) 'Fitting module'
      write(*,*) '--------------'
      write(*,*) 'Number of trajectory files  : ', Nf
      write(*,*) 'Minimum allowed speed       : ', fit_vmin
      write(*,*) 'Maximum allowed speed       : ', fit_vmax
      write(*,*) 'Number valid observations   : ', nvalid
      do i=1,Nf
        write(*,*) 'File name                 : ', trim(filelist(i))
        write(*,*) 'Number particles in file  : ', To(i)%Np
        write(*,*) 'Number of steps in file   : ', To(i)%Nt
        write(*,*) 'Time step                 : ', To(i)%dt
        do pp=1,To(i)%Np
          write(*,*) '  Particle, valid entries : ', pp, count(To(i)%valid(pp,1:To(i)%Nt-1))
        enddo
      enddo
    endif

    ! ... Check if wind data has been loaded
    ! ...
    if (.not.WithAU) then
      fit_FGp2 = 0.0D0
      fit_FGp3 = 0.0D0
      fit_dop2 = 0.0D0
      fit_dop3 = 0.0D0
    endif
    if (.not.WithAV) then
      fit_FGp4 = 0.0D0
      fit_FGp5 = 0.0D0
      fit_dop4 = 0.0D0
      fit_dop5 = 0.0D0
    endif


    ! ... The number of Valid matchings is nvalid. We have to multiply by two to take
    ! ... into account the U ad the V matchings
    ! ...
    fit_nobs = 0
    do i=1,size(To)
    do pp=1,To(i)%Np
    do ll=1,To(i)%Nt-1
      if (To(i)%valid(pp,ll)) fit_nobs = fit_nobs + 2
    enddo
    enddo
    enddo

    allocate(Fvec(fit_nobs))       ! Number of terms   yo(i) - ym(i,p)
    allocate(Fjac(fit_nobs,np))    ! Jacobian of each term by respect p

    p(:)   = [fit_FGp1, fit_FGp2, fit_FGp3, fit_FGp4, fit_FGp5 ]

    Jcosto =  func(p)
    dJcost = dfunc(p)

    iu = unitfree()
    open(iu,file=fit_diag,status='unknown')
    write(iu,'(T2,A,I5)')     'Number obs       : ', fit_nobs
    write(iu,'(T2,A,5F15.5)') 'First guess      : ', p
    write(iu,'(T2,A,G15.5)')  'Initial cost     : ', Jcosto
    write(iu,'(T2,A,5G15.5)') 'Initial gradient : ', dJcost
    write(iu,'(T2,A,5F15.0)') 'Fitting flag     : ', fit_dop1, fit_dop2, fit_dop3, fit_dop4, fit_dop5
    write(iu,'(T2,A,G12.3)')  'Ftol             : ', fit_Ftol
    write(iu,'(T2,A,G12.3)')  'Xtol             : ', fit_Xtol
    write(iu,'(T2,A,G12.3)')  'Gtol             : ', fit_Gtol
    write(iu,'(T2,A,I3)')     'Maxfev           : ', fit_Maxfev
    write(iu,'(T2,A,I3)')     'Mode             : ', fit_Mode   
    write(iu,'(T2,A,G12.3)')  'Factor           : ', fit_Factor
    write(iu,'(T2,A,I3)')     'Nprint           : ', fit_Nprint 

    if (verb.ge.1) then
      write(*,*) 
      write(*,'(T2,A,I5)')     'Number obs       : ', fit_nobs
      write(*,'(T2,A,5F15.5)') 'First guess      : ', p
      write(*,'(T2,A,G15.5)')  'Initial cost     : ', Jcosto
      write(*,'(T2,A,5G15.5)') 'Initial gradient : ', dJcost
      write(*,'(T2,A,5F15.0)') 'Fitting flag     : ', fit_dop1, fit_dop2, fit_dop3, fit_dop4, fit_dop5
      write(*,'(T2,A,G15.5)')  'Ftol             : ', fit_Ftol
      write(*,'(T2,A,G15.5)')  'Xtol             : ', fit_Xtol
      write(*,'(T2,A,G15.5)')  'Gtol             : ', fit_Gtol
      write(*,'(T2,A,I3)')     'Maxfev           : ', fit_Maxfev
      write(*,'(T2,A,I3)')     'Mode             : ', fit_Mode   
      write(*,'(T2,A,G9.3)')   'Factor           : ', fit_Factor
      write(*,'(T2,A,I3)')     'Nprint           : ', fit_Nprint 
    endif

    call lmder (FCN,fit_nobs,Np,p,Fvec,Fjac,fit_nobs, &
                fit_Ftol,fit_Xtol,fit_Gtol,fit_Maxfev, &
                diag,fit_Mode,fit_Factor,fit_Nprint, &
                Ninfo,Nfev,Njev,ipvt,qtf)

    Jcost = func(p)
    dJcost = dfunc(p)
    write(iu,'(T2,A,5F15.5)') 'Solution         : ', p
    write(iu,'(T2,A,F6.4)')   '||A||_2          : ', sqrt(p(2)**2+p(3)**2+p(4)**2+p(5)**2)
    write(iu,'(T2,A,G15.5)')  'Final cost       : ', Jcost
    write(iu,'(T2,A,5G15.5)') 'Final gradient   : ', dJcost
    write(iu,'(T2,A,I3)')     'Exit parameter   : ', Ninfo  
    write(iu,'(T2,A,I3)')     'NFev             : ', Nfev
    write(iu,'(T2,A,I3)')     'NJev             : ', Njev
    close(iu)

    if (verb.ge.1) then
      write(*,*)
      write(*,'(T2,A,5F15.5)') 'Solution         : ', p
      write(*,'(T2,A,F6.4)')   '||A||_2          : ', sqrt(p(2)**2+p(3)**2+p(4)**2+p(5)**2)
      write(*,'(T2,A,G15.5)')  'Final cost       : ', Jcost
      write(*,'(T2,A,5G15.5)') 'Final gradient   : ', dJcost
      write(*,'(T2,A,I3)')     'Exit parameter   : ', Ninfo  
      write(*,'(T2,A,I3)')     'NFev             : ', Nfev
      write(*,'(T2,A,I3)')     'NJev             : ', Njev
    endif

    if (Ninfo.eq.0) call crash('MODEL_FITTING: Inproper input parameters in LMDER')

    open(iu,file=fit_fout,status='unknown')
    write(iu,'(T1,A)') '# File  Part  Stp   UObs     VObs    Umodel   Vmodel   OceU     OceV     AtmU     AtmV'
    write(iu,'(T1,A)') '# --------------------------------------------------------------------------------------'
    do i=1,size(To)
    do pp=1,To(i)%Np
    do ll=1,To(i)%Nt-1
      if (To(i)%valid(pp,ll)) then
        um = p(1)*To(i)%OU(pp,ll) + p(2)*To(i)%AU(pp,ll) + p(3)*To(i)%AV(pp,ll)
        vm = p(1)*To(i)%OV(pp,ll) + p(4)*To(i)%AU(pp,ll) + p(5)*To(i)%AV(pp,ll)
        write(iu,'(I4,I5,I6,8F9.3)') i, pp, ll, To(i)%uo(pp,ll), To(i)%vo(pp,ll), um, vm,  &
                            To(i)%OU(pp,ll), To(i)%OV(pp,ll), To(i)%AV(pp,ll), To(i)%AV(pp,ll)              
      endif
    enddo
    enddo
    enddo
    close(iu)
    
    return


  end subroutine model_fitting
  ! ...
  ! =====================================================================
  ! ...
  subroutine trajectory_read(T,filename)

    class(type_trajectory), intent(inout)        :: T
    character(len=*), intent(in)                 :: filename  

    ! ... Local variables
    ! ...
    character(len=maxlen)                        :: fpath,fbase,ftype
    integer                                      :: pp,ll
    real(dp)                                     :: dt

    call filename_split (filename,fpath,fbase,ftype)
    ftype = lowercase(ftype)

    if (index(ftype,'dat').gt.0) then
      call trajectory_read_ascii (filename,T)
 
    else if ((ftype.eq.'nc').or.(ftype.eq.'cdf')) then
      call trajectory_read_nc (filename,T)

    else if (index(ftype,'json').gt.0) then
      call crash('MODULE_FITTING: JSON files not accepted.')

    else
      call crash('Input trajectory file: unknown format')
    endif

    ! ... Check if time step is homogeneous
    ! ...
    T%dt = anint(T%t(1,2) - T%t(1,1))
    if (T%dt.eq.0) call crash('Time step cannot be zero')

    do pp=1,T%Np
    do ll=2,T%Nt-1
      dt = anint(T%t(pp,ll+1)-T%t(pp,ll))
      if (dt.ne.T%dt) call crash('Time step not uniform')
    enddo
    enddo

  end subroutine trajectory_read
  ! ...
  ! =====================================================================
  ! ...
  subroutine trajectory_read_ascii (filename,T)

    character(len=*), intent(in)                 :: filename  
    type(type_trajectory), intent(inout)         :: T

    ! ... Local variables
    ! ...
    integer iu,i,nl,nh,N,err
    real(dp) x, y, z
    character(len=maxlen) line
    character(len=20) sdate

    iu = unitfree()
    if (verb.ge.2) write(*,*) 'Opening ASCII trajectory ', trim(filename)
    open(iu,file=filename,status='old',iostat=err)
    if (err.ne.0) call crash('File '//trim(filename)//' does not exist')

    nl = numlines(iu,'ascii')
    if (verb.ge.3) write(*,*) 'Number of lines: ', nl

    ! ... Header lines:
    ! ...
    nh = 0 
    do i=1,nl
      read(iu,'(A)') line
      if (line(1:1).ne.'#') exit
      nh = nh + 1
    enddo
    if (verb.ge.3) write(*,*) 'Number of header lines: ', nh

    N = nl - nh

    ! ... Allocate all the tables of the trajectory structure
    ! ...
    call T%allocate(1,N)
    T%valid(1,:) = .True.

    rewind(iu) 
    do i=1,nh
      read(iu,*)
    enddo
    do i=1,N
      read(iu,*) T%x(1,i), T%y(1,i), T%z(1,i), sdate
      T%date(1,i) = strptime(sdate)
      if ((T%x(1,i).lt.-180.001).or.(T%x(1,i).gt.360.001)) T%valid(1,i) = .False.
      if ((T%y(1,i).lt.-90.001).or.(T%y(1,i).gt.90.001))   T%valid(1,i) = .False.
    enddo
    
    ! ... Transform dates into seconds:
    ! ...
    do i=1,N
      T%t(1,i) = anint(date2num(T%date(1,i),units=alm_time_units))
      if ((T%t(1,i).lt.alm_tini).or.(T%t(1,i).gt.alm_tfin)) then
        T%valid(1,i) = .False.
      endif
    enddo

  end subroutine trajectory_read_ascii
  ! ...
  ! =====================================================================
  ! ...
  subroutine trajectory_read_nc (filename,T)

    character(len=*), intent(in)                 :: filename  
    type(type_trajectory), intent(inout)         :: T

    ! ... Local variables
    ! ...
    integer fid,np,nt,idx,idy,idz,idt,err,pp,ll
    character(len=maxlen) dname,time_units,calendar

    if (verb.ge.2) write(*,*) 'Opening NetCDF trajectory ', trim(filename)
    err = NF90_OPEN(filename,0,fid)
    call cdf_error(err,'TRAJECTORY_READ_NC: error opening '//trim(filename))

    err = NF90_INQUIRE_DIMENSION(fid,1,name=dname,len=np)
    call cdf_error(err,'TRAJECTORY_READ_NC: error reading first dimension')
    err = NF90_INQUIRE_DIMENSION(fid,2,name=dname,len=nt)
    call cdf_error(err,'TRAJECTORY_READ_NC: error reading second dimension')
    if (verb.ge.2) write(*,*) 'Np, Nt : ', Np, Nt

    err = NF90_INQ_VARID(fid,'lon',idx)
    call cdf_error(err,'TRAJECTORY_READ_NC: error inquiring about lon')
    err = NF90_INQ_VARID(fid,'lat',idy)
    call cdf_error(err,'TRAJECTORY_READ_NC: error inquiring about lat')
    err = NF90_INQ_VARID(fid,'depth',idz)
    call cdf_error(err,'TRAJECTORY_READ_NC: error inquiring about depth')
    err = NF90_INQ_VARID(fid,'time',idt)
    call cdf_error(err,'TRAJECTORY_READ_NC: error inquiring about time')

    T%filename = trim(filename)
    T%Np = np
    T%Nt = nt

    ! ... Allocate all the tables of the trajectory structure
    ! ...
    call T%allocate(Np,Nt)
    T%valid(:,:) = .True.

    err = NF90_GET_VAR(fid,idx,T%x)
    call cdf_error(err,'TRAJECTORY_READ_NC: error reading lon')
    err = NF90_GET_VAR(fid,idy,T%y)
    call cdf_error(err,'TRAJECTORY_READ_NC: error reading lat')
    err = NF90_GET_VAR(fid,idz,T%z)
    call cdf_error(err,'TRAJECTORY_READ_NC: error reading depth')
    err = NF90_GET_VAR(fid,idt,T%t(1,:))
    call cdf_error(err,'TRAJECTORY_READ_NC: error reading time')


    ! ... Convert time to dates:
    ! ...
    time_units = ''
    calendar = ''
    err = NF90_GET_ATT (fid,idt,'units',time_units)
    call cdf_error(err,'Time variable without UNITS attribute')
    err = NF90_GET_ATT (fid,idt,'calendar',calendar)
    call cdf_error(err,'Time variable without CALENDAR attribute')
    call check_calendar(calendar)

    do ll=1,T%nt
      T%date(1,ll) = num2date(T%t(1,ll),units=time_units,calendar=calendar)
    enddo

    ! ... replicate time and dates:
    ! ...
    do pp=2,np
      T%t(pp,:) = T%t(1,:)
      T%date(pp,:) = T%date(1,:)
    enddo

    ! ... Close file
    ! ...
    err = NF90_CLOSE(fid)

    ! ... Remove not valid points:
    ! ...
    do pp=1,np
    do ll=1,nt
      if ((T%x(pp,ll).lt.-180.001).or.(T%x(pp,ll).gt.360.001))  T%valid(pp,ll) = .False.
      if ((T%y(pp,ll).lt.-90.001).or.(T%y(pp,ll).gt.90.001))    T%valid(pp,ll) = .False.
      if ((T%t(pp,ll).lt.alm_tini).or.(T%t(pp,ll).gt.alm_tfin)) T%valid(pp,ll) = .False.
    enddo
    enddo
    
  end subroutine trajectory_read_nc
  ! ...
  ! =====================================================================
  ! ...
  subroutine trajectory_allocate(T,Np,N)

  class(type_trajectory), intent(inout)          :: T
  integer, intent(in)                            :: Np
  integer, intent(in)                            :: N

  T%valid => NULL()
  T%date  => NULL()
  T%x     => NULL()
  T%y     => NULL()
  T%z     => NULL()
  T%t     => NULL()
  T%uf    => NULL()
  T%vf    => NULL()
  T%wf    => NULL()
  T%uo    => NULL()
  T%vo    => NULL()
  T%wo    => NULL()
  T%speed => NULL()
  T%OU    => NULL()
  T%OV    => NULL()
  T%OW    => NULL()
  T%AU    => NULL()
  T%AV    => NULL()

  T%Np = Np
  T%Nt = N

  allocate(T%valid(Np,N))
  allocate(T%date(Np,N))
  allocate(T%x(Np,N))
  allocate(T%y(Np,N))
  allocate(T%z(Np,N))
  allocate(T%t(Np,N))
  allocate(T%uf(Np,N))
  allocate(T%vf(Np,N))
  allocate(T%wf(Np,N))
  allocate(T%uo(Np,N))
  allocate(T%vo(Np,N))
  allocate(T%wo(Np,N))
  allocate(T%speed(Np,N))
  allocate(T%OU(Np,N))
  allocate(T%OV(Np,N))
  allocate(T%OW(Np,N))
  allocate(T%AU(Np,N))
  allocate(T%AV(Np,N))

  end subroutine trajectory_allocate
  ! ...
  ! =====================================================================
  ! ...
  subroutine trajectory_compress(T)

    class(type_trajectory), intent(inout)          :: T

    ! ... Local variables
    ! ... 
    integer i,iv,N,Nv
    type(type_trajectory) Tv

!    N = T%N
!    Nv = count(T%valid)
!
!    call Tv%allocate(Nv)
!
!    iv = 0
!    do i=1,N
!      if (T%valid(i)) then
!        iv = iv + 1
!        Tv%valid(iv) = T%valid(i)
!        Tv%date(iv)  = T%date(i)
!        Tv%x(iv)     = T%x(i)
!        Tv%y(iv)     = T%y(i)
!        Tv%z(iv)     = T%z(i)
!        Tv%t(iv)     = T%t(i)
!        Tv%dt(iv)    = T%dt(i)
!        Tv%uf(iv)    = T%uf(i)
!        Tv%vf(iv)    = T%vf(i)
!        Tv%wf(iv)    = T%wf(i)
!        Tv%uo(iv)    = T%uo(i)
!        Tv%vo(iv)    = T%vo(i)
!        Tv%wo(iv)    = T%wo(i)
!        Tv%OU(iv)    = T%OU(i)
!        Tv%OV(iv)    = T%OV(i)
!        Tv%OW(iv)    = T%OW(i)
!        Tv%AU(iv)    = T%AU(i)
!        Tv%AV(iv)    = T%AV(i)
!      endif
!    enddo
!    if (iv.ne.Nv) stop 'Incompatible iv and Nv in trajectory_compress'

!    call T%allocate(Nv)
!
!    T%valid(:) = Tv%valid(:)
!    T%date(:)  = Tv%date(:)
!    T%x(:)     = Tv%x(:)
!    T%y(:)     = Tv%y(:)
!    T%z(:)     = Tv%z(:)
!    T%t(:)     = Tv%t(:)
!    T%dt(:)    = Tv%dt(:)
!    T%uf(:)    = Tv%uf(:)
!    T%vf(:)    = Tv%vf(:)
!    T%wf(:)    = Tv%wf(:)
!    T%uo(:)    = Tv%uo(:)
!    T%vo(:)    = Tv%vo(:)
!    T%wo(:)    = Tv%wo(:)
!    T%OU(:)    = Tv%OU(:)
!    T%OV(:)    = Tv%OV(:)
!    T%OW(:)    = Tv%OW(:)
!    T%AU(:)    = Tv%AU(:)
!    T%AV(:)    = Tv%AV(:)
! 
!    Tv%valid => NULL()
!    Tv%date  => NULL()
!    Tv%x     => NULL()
!    Tv%y     => NULL()
!    Tv%z     => NULL()
!    Tv%t     => NULL()
!    Tv%dt    => NULL()
!    Tv%uf    => NULL()
!    Tv%vf    => NULL()
!    Tv%wf    => NULL()
!    Tv%uo    => NULL()
!    Tv%vo    => NULL()
!    Tv%wo    => NULL()
!    Tv%speed => NULL()
!    Tv%OU    => NULL()
!    Tv%OV    => NULL()
!    Tv%OW    => NULL()
!    Tv%AU    => NULL()
!    Tv%AV    => NULL()

  end subroutine trajectory_compress
  ! ...
  ! =====================================================================
  ! ...
  subroutine trajectory_velocity(T)

    class(type_trajectory), intent(inout)          :: T

    ! ... Local variables
    ! ... 
    integer i,iv,N,Nv,pp,ll
    real(dp) dt,ym,dx,dy,dz,dtmean
    real(dp) u,v,w

    T%uo(:,:) = 0.0D0
    T%vo(:,:) = 0.0D0
    T%wo(:,:) = 0.0D0
    dt  = T%dt

    nv      = 0
    do pp=1,T%Np
    do ll=1,T%Nt-1

      dx = T%x(pp,ll+1) - T%x(pp,ll)
      dy = T%y(pp,ll+1) - T%y(pp,ll)
      dz = T%z(pp,ll+1) - T%z(pp,ll)

      if (abs(dx).gt.0.and.abs(dy).gt.0) then

        nv = nv + 1
        if (Spherical) then
          ym =  0.5D0*(T%y(pp,ll+1)+T%y(pp,ll))
          dx = REarth*dx*cos(ym)
          dy = REarth*dy
        endif
        u = dx/dt
        v = dy/dt
        w = dz/dt
        !print*, pp, ll
        !print*, T%x(pp,ll), T%y(pp,ll) 
        !print*, T%x(pp,ll+1), T%y(pp,ll+1) 
        !print*, dx, dy, dz, dt, u, v, w
        !stop '33333'
        T%uo(pp,ll)     = u
        T%vo(pp,ll)     = v
        T%wo(pp,ll)     = w
        T%speed(pp,ll) = sqrt(u*u + v*v + w*w)
        if ((T%speed(pp,ll).lt.fit_vmin).or.(T%speed(pp,ll).gt.fit_vmax)) T%valid(pp,ll) = .False.

      else

        T%valid(pp,ll) = .False.

      endif

    enddo
    enddo
    if (nv.eq.0) call crash ('No valid trajectory data')


  end subroutine trajectory_velocity
  ! ...
  ! =====================================================================
  ! ...
  real(dp) function func(p)

    real(dp), dimension(:), intent(in)       :: p

    ! ... Local variables
    ! ...
    integer f,ll,pp,n
    real(dp) u,v,du,dv

    ! ... p(5)
    ! ... alpha = p(1)
    ! ... a11   = p(2)
    ! ... a12   = p(3)
    ! ... a21   = p(4)
    ! ... a22   = p(5)
    ! ...
    func = 0.0D0
    n    = 0.0D0
    do f=1,size(To)
      do pp=1,To(f)%Np
      do ll=1,To(f)%Nt-1
        if (To(f)%valid(pp,ll)) then
          u = p(1)*To(f)%OU(pp,ll) + p(2)*To(f)%AU(pp,ll) + p(3)*To(f)%AV(pp,ll)
          v = p(1)*To(f)%OV(pp,ll) + p(4)*To(f)%AU(pp,ll) + p(5)*To(f)%AV(pp,ll)
          du = To(f)%uo(pp,ll) - u
          dv = To(f)%vo(pp,ll) - v
          func = func + du*du + dv*dv
          !print*, f, pp, ll,  To(f)%OU(pp,ll), To(f)%OV(pp,ll)           ! DEBUG
          n = n + 1
        endif
      enddo
      enddo
    enddo
    func = 0.5D0 * func / n
    !stop 'DEBUG'

  end function func
  ! ...
  ! ===================================================================
  ! ...
  function dfunc(p) result(xi)

    real(dp), dimension(:), intent(in)       :: p
    real(dp), dimension(size(p))             :: xi

    ! ... Local variables
    ! ...
    integer f,pp,ll,n
    real(dp) u,v,du,dv

    ! ... p(5)
    ! ... alpha = p(1)
    ! ... a11   = p(2)
    ! ... a12   = p(3)
    ! ... a21   = p(4)
    ! ... a22   = p(5)
    ! ...
    xi(:) = 0.0D0
    n     = 0.0D0
    do f=1,size(To)
      do pp=1,To(f)%Np
      do ll=1,To(f)%Nt-1
        if (To(f)%valid(pp,ll)) then
          u = p(1)*To(f)%OU(pp,ll) + p(2)*To(f)%AU(pp,ll) + p(3)*To(f)%AV(pp,ll)
          v = p(1)*To(f)%OV(pp,ll) + p(4)*To(f)%AU(pp,ll) + p(5)*To(f)%AV(pp,ll)
          du = To(f)%uo(pp,ll) - u
          dv = To(f)%vo(pp,ll) - v
          xi(1) = xi(1) + du*To(f)%OU(pp,ll) + dv*To(f)%OV(pp,ll)
          xi(2) = xi(2) + du*To(f)%AU(pp,ll)
          xi(3) = xi(3) + du*To(f)%AV(pp,ll)
          xi(4) = xi(4) + dv*To(f)%AU(pp,ll)
          xi(5) = xi(5) + dv*To(f)%AV(pp,ll)
          n = n + 1
        endif
      enddo
      enddo
    enddo
    xi(:) = -xi(:)/n

  end function dfunc
  ! ...
  ! ===================================================================
  ! ...
  subroutine fcn (M,Np,P,FVEC,FJAC,LDFJAC,IFLAG)

  ! -----------------------------------------------------
  ! ...  IF IFLAG = 1 CALCULATE THE FUNCTIONS AT X AND
  ! ...  RETURN THIS VECTOR IN FVEC.  DO NOT ALTER FJAC.
  ! ...  IF IFLAG = 2 CALCULATE THE JACOBIAN AT X AND
  ! ...  RETURN THIS MATRIX IN FJAC.  DO NOT ALTER FVEC.
  ! -----------------------------------------------------

    integer, intent(in)            :: M,Np,LDFJAC,IFLAG
    real(dp), intent(in)           :: P(Np)
    real(dp), intent(inout)        :: FVEC(M)
    real(dp), intent(inout)        :: FJAC(LDFJAC,Np)

    ! ... Local variables
    ! ...
    integer f,pp,ll,ii
    real(dp) um,vm

    ! ... p(5)
    ! ... alpha = p(1)
    ! ... a11   = p(2)
    ! ... a12   = p(3)
    ! ... a21   = p(4)
    ! ... a22   = p(5)
    ! ...

    if (iflag.eq.1) then
      ii = 0
      do f=1,size(To)
        do pp=1,To(f)%Np
        do ll=1,To(f)%Nt-1
          if (To(f)%valid(pp,ll)) then
            um = p(1)*To(f)%OU(pp,ll) + p(2)*To(f)%AU(pp,ll) + p(3)*To(f)%AV(pp,ll)
            vm = p(1)*To(f)%OV(pp,ll) + p(4)*To(f)%AU(pp,ll) + p(5)*To(f)%AV(pp,ll)
            ii = ii + 1
            Fvec(ii) = To(f)%uo(pp,ll) - um
            ii = ii + 1
            Fvec(ii) = To(f)%vo(pp,ll) - vm
          endif
        enddo
        enddo
      enddo
      if  (ii.ne.M) call crash('Error in number of matchings')
      return
    endif

    if (iflag.eq.2) then
      ii = 0
      do f=1,size(To)
        do pp=1,To(f)%Np
        do ll=1,To(f)%Nt-1
          if (To(f)%valid(pp,ll)) then
            ii = ii + 1
            FJAC(ii,1) = -To(f)%OU(pp,ll) * fit_dop1
            FJAC(ii,2) = -To(f)%AU(pp,ll) * fit_dop2
            FJAC(ii,3) = -To(f)%AV(pp,ll) * fit_dop3
            FJAC(ii,4) =  0.0D0
            FJAC(ii,5) =  0.0D0
            ii = ii + 1
            FJAC(ii,1) = -To(f)%OV(pp,ll) * fit_dop1
            FJAC(ii,2) =  0.0D0
            FJAC(ii,3) =  0.0D0
            FJAC(ii,4) = -To(f)%AU(pp,ll) * fit_dop4
            FJAC(ii,5) = -To(f)%AV(pp,ll) * fit_dop5
          endif
        enddo
        enddo
      enddo
      if  (ii.ne.M) call crash('Error in number of matchings')
      return
    endif

  end subroutine fcn
  ! ...
  ! ===================================================================
  ! ...
end module module_fitting
