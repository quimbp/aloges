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
  integer                                        :: No = 0
  integer                                        :: N  = 0
  integer                                        :: Nm = 0
  character(len=maxlen)                          :: filename = ''
  logical, dimension(:),  pointer                :: valid
  type(type_date), dimension(:), pointer         :: date
  real(dp)                                       :: dtmean = 0.0D0
  real(dp), dimension(:), pointer                :: x
  real(dp), dimension(:), pointer                :: y
  real(dp), dimension(:), pointer                :: z
  real(dp), dimension(:), pointer                :: t
  real(dp), dimension(:), pointer                :: dt
  real(dp), dimension(:), pointer                :: uf       ! File (if stored) u value
  real(dp), dimension(:), pointer                :: vf       ! File (if stored) v value
  real(dp), dimension(:), pointer                :: wf       ! File (if stored) w value
  real(dp), dimension(:), pointer                :: uo       ! Derived value
  real(dp), dimension(:), pointer                :: vo       ! Derived value
  real(dp), dimension(:), pointer                :: wo       ! Derived value
  real(dp), dimension(:), pointer                :: speed    ! Derived value
  real(dp), dimension(:), pointer                :: OU       ! Interpolated OU value
  real(dp), dimension(:), pointer                :: OV       ! Interpolated OV value
  real(dp), dimension(:), pointer                :: OW       ! Interpolated OW value
  real(dp), dimension(:), pointer                :: AU       ! Interpolated AU value
  real(dp), dimension(:), pointer                :: AV       ! Interpolated AV value

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
character(len=maxlen)                            :: fit_fout  = 'fitting.dat'


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
    real(dp)                                         :: dtmin,dtmax,dtmean
    real(dp)                                         :: um,vm
    real(dp), dimension(ndims)                       :: xo
    character(maxlen), dimension(:), allocatable     :: filelist

    ! ... Minimization variables
    ! ...
    integer, parameter                               :: Np = 5
    integer                                          :: fit_nobs
    integer                                          :: Maxfev
    integer                                          :: Mode
    integer                                          :: Nprint
    integer                                          :: Ninfo
    integer                                          :: NFev
    integer                                          :: NJev
    integer, dimension(Np)                           :: ipvt
    real(dp), dimension(:), allocatable              :: Fvec
    real(dp), dimension(:,:), allocatable            :: FJac
    real(dp)                                         :: Ftol,Xtol,Gtol
    real(dp)                                         :: Factor
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
    if (.not.WithAU) AUrhs = 0.0D0
    if (.not.WithAV) AVrhs = 0.0D0

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

      ! ... Read (one file at a time):
      ! ...
      call To(i)%read(filelist(i))

      ! ... Pass from deg to radians
      ! ...
      To(i)%x(:) = deg2rad*To(i)%x(:)
      To(i)%y(:) = deg2rad*To(i)%y(:)
      To(i)%z(:) =    -abs(To(i)%z(:))

      ! ... Calculate the velocities
      ! ... It calculates the dt between successive valid positions
      ! ...
      call To(i)%velocity()

      ! ... Minimum and maximum time intervals between "valid" positions
      ! ...
      dtmin = 1D40
      dtmax = -dtmin
      do j=1,To(i)%Nm
        if (To(i)%valid(j)) then
          dtmin = min(dtmin,To(i)%dt(j))
          dtmax = max(dtmax,To(i)%dt(j))
        endif
      enddo

      ! ... Trimmed mean:
      ! ...
      firstmin = .True.
      firstmax = .True.
      dtmean = 0.0D0
      icount = 0
      do j=1,To(i)%Nm
        if (To(i)%valid(j)) then
          if (firstmin.and.To(i)%dt(j).eq.dtmin) then
            firstmin = .False.
          else if  (firstmax.and.To(i)%dt(j).eq.dtmax) then
            firstmax = .False.
          else
            dtmean = dtmean + To(i)%dt(j)
            icount = icount + 1
          endif
        endif
      enddo
      if (icount.eq.0) call crash('No valid points in trajectory')
      dtmean = dtmean / icount
      To(i)%dtmean = dtmean

      ! ... Large and small dt intervals are not considered:
      ! ... This eliminates the first velocity estimation after a big data gap
      ! ... in the time series. "Big" data gaps correspond to more than 3 times
      ! ... the averaged time interval
      ! ...
      do j=1,To(i)%Nm
        if (To(i)%valid(j)) then
           if (To(i)%dt(j).lt.0.333D0*dtmean) To(i)%valid(j) = .False.
           if (To(i)%dt(j).gt.3.000D0*dtmean) To(i)%valid(j) = .False.
        endif
      enddo
   
      if (verb.ge.5) then
        write(*,*)
        write(*,*) 'Estimated velocities: '
        do j=1,To(i)%Nm
          write(*,'(I4,5F9.4,F7.0,X,L)') j, To(i)%x(j), To(i)%y(j), To(i)%uo(j), To(i)%vo(j), &
                                            To(i)%wo(j), To(i)%dt(j), To(i)%valid(j)
        enddo
      endif

      ! ... Reset the forcing field counters
      ! ...
      GOU%rec1 = -1; GOU%rec2 = -1
      GOV%rec1 = -1; GOV%rec2 = -1
      GOW%rec1 = -1; GOW%rec2 = -1
      GAU%rec1 = -1; GAU%rec2 = -1
      GAV%rec1 = -1; GAV%rec2 = -1

      do j=1,To(i)%Nm
        if (To(i)%valid(j)) then
          xo = [To(i)%x(j), To(i)%y(j), To(i)%z(j) ]

          To(i)%OU(j) = 0.0D0; To(i)%OV(j) = 0.0D0; To(i)%OW(j) = 0.0D0; To(i)%AU(j) = 0.0D0; To(i)%AV(j) = 0.0D0
          call get_forcing('OU',GOU,To(i)%t(j),To(i)%date(j),0.0D0,OUrhs)
          To(i)%OU(j) = time_interpol(GOU,OUrhs,To(i)%t(j),xo)
          call get_forcing('OV',GOV,To(i)%t(j),To(i)%date(j),0.0D0,OVrhs)
          To(i)%OV(j) = time_interpol(GOV,OVrhs,To(i)%t(j),xo)

          if (WithOW) then
            call get_forcing('OW',GOW,To(i)%t(j),To(i)%date(j),0.0D0,OWrhs)
            To(i)%OW(j) = time_interpol(GOW,OWrhs,To(i)%t(j),xo)
          else
            To(i)%OW(j) = 0.0D0
          endif
          if (winds) then
            call get_forcing('AU',GAU,To(i)%t(j),To(i)%date(j),0.0D0,AUrhs)
            To(i)%AU(j) = time_interpol(GAU,AUrhs,To(i)%t(j),xo)
            call get_forcing('AV',GAV,To(i)%t(j),To(i)%date(j),0.0D0,AVrhs)
            To(i)%AV(j) = time_interpol(GAV,AVrhs,To(i)%t(j),xo)
          else
            To(i)%AU(j) = 0.0D0
            To(i)%AV(j) = 0.0D0
          endif
          if (verb.ge.4) then
            write(*,*) 'Position number    : ', j
            write(*,*) 'MODEL_FIT Time, xo : ', To(i)%t(j), xo
            write(*,*) '      OU, OV, OW   : ', To(i)%OU(j), To(i)%OV(j), To(i)%OW(j)
            write(*,*) '      AU, AV       : ', To(i)%AU(j), To(i)%AV(j)
          endif
        endif
      enddo

    enddo

    ! ... At this moment we have everything for the evaluation of the cost function
    ! ... The information is on the To structure.

    nvalid = 0
    do i=1,Nf
      nvalid = nvalid + count(To(i)%valid)
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
        write(*,*) 'Number of entries in file : ', To(i)%No
        write(*,*) 'Number valid entries      : ', count(To(i)%valid(1:To(i)%Nm))
        write(*,*) 'Average time step         : ', To(i)%dtmean
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
    do j=1,To(i)%Nm
      if (To(i)%valid(j)) fit_nobs = fit_nobs + 2
    enddo
    enddo

    allocate(Fvec(fit_nobs))       ! Number of terms   yo(i) - ym(i,p)
    allocate(Fjac(fit_nobs,np))    ! Jacobian of each term by respect p

    p(:)   = [fit_FGp1, fit_FGp2, fit_FGp3, fit_FGp4, fit_FGp5 ]
    Ftol   = 1.0D-6
    Xtol   = 1.0D-7
    Gtol   = 1.0D-18
    Maxfev = 100
    Mode   = 1
    factor = 10.0D0
    Nprint = 0

    if (verb.ge.1) then
        Jcosto =  func(p)
        dJcost = dfunc(p)
        write(*,*) 
        write(*,*) 'Number obs       : ', fit_nobs
        write(*,'(T2,A,5F15.5)') 'First guess      : ', p
        write(*,'(T2,A,G15.5)') 'Initial cost     : ', Jcosto
        write(*,'(T2,A,5G15.5)') 'Initial gradient : ', dJcost
        write(*,'(T2,A,5F15.0)') 'Fitting flag     : ', fit_dop1, fit_dop2, fit_dop3, fit_dop4, fit_dop5
    endif

    call lmder (FCN,fit_nobs,Np,p,Fvec,Fjac,fit_nobs, &
                Ftol,Xtol,Gtol,Maxfev, &
                diag,mode,factor,nprint, &
                Ninfo,Nfev,Njev,ipvt,qtf)

    if (verb.ge.1) then
        Jcost = func(p)
        dJcost = dfunc(p)
        write(*,*)
        write(*,'(T2,A,5F15.5)') 'Solution         : ', p
        write(*,'(T2,A,G15.5)') 'Final cost       : ', Jcost
        write(*,'(T2,A,5G15.5)') 'Final gradient   : ', dJcost
        write(*,*) 'Info             : ', ninfo
        write(*,*) 'Nfev             : ', nfev
        write(*,*) 'Njev             : ', njev
    endif

    open(iu,file=fit_fout,status='unknown')
    do i=1,size(To)
    do j=1,To(i)%Nm
      if (To(i)%valid(j)) then
        um = p(1)*To(i)%OU(j) + p(2)*To(i)%AU(j) + p(3)*To(i)%AV(j)
        vm = p(1)*To(i)%OV(j) + p(4)*To(i)%AU(j) + p(5)*To(i)%AV(j)
        write(iu,'(2I4,8F9.3)') i, j, To(i)%uo(j), To(i)%vo(j), um, vm,  &
                            To(i)%OU(j), To(i)%OV(j), To(i)%AV(j), To(i)%AV(j)              
      endif
    enddo
    enddo
    close(iu)
    

!    p(:) = 1.0D0
!   Jcosto = func(p)
!    print*, 'Cost: ', Jcosto
!    dJcost = dfunc(p)
!    print*, 'd Cost/ dp = ', dJcost
!
!    pn = p; pn(2) = 1.01D0
!    Jcost = func(pn)
!    print*, (Jcost - Jcosto)/0.01D0


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
    integer                                      :: i

    call filename_split (filename,fpath,fbase,ftype)
    ftype = lowercase(ftype)

    if (index(ftype,'dat').gt.0) then
      call trajectory_read_ascii (filename,T)
 
    else if ((ftype.eq.'nc').or.(ftype.eq.'cdf')) then
      print*, 'netcdf'

    else if (index(ftype,'json').gt.0) then
      print*, 'geojson'

    else
      call crash('Input trajectory file: unknown format')
    endif

    ! ... Transform dates into seconds:
    ! ...
    do i=1,T%N
      T%t(i) = anint(date2num(T%date(i),units=alm_time_units))
      if ((T%t(i).lt.alm_tini).or.(T%t(i).gt.alm_tfin)) then
        T%valid(i) = .False.
      endif
    enddo

    T%No = T%N
    call T%compress()
  

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
    call T%allocate(N)
    T%valid(:) = .True.

    rewind(iu) 
    do i=1,nh
      read(iu,*)
    enddo
    do i=1,N
      read(iu,*) T%x(i), T%y(i), T%z(i), sdate
      T%date(i) = strptime(sdate)
      if ((T%x(i).lt.-180.001).or.(T%x(i).gt.360.001)) T%valid(i) = .False.
      if ((T%y(i).lt.-90.001).or.(T%y(i).gt.90.001))   T%valid(i) = .False.
    enddo
    
  end subroutine trajectory_read_ascii
  ! ...
  ! =====================================================================
  ! ...
  subroutine trajectory_allocate(T,N)

  class(type_trajectory), intent(inout)          :: T
  integer, intent(in)                            :: N

  T%valid => NULL()
  T%date  => NULL()
  T%x     => NULL()
  T%y     => NULL()
  T%z     => NULL()
  T%t     => NULL()
  T%dt    => NULL()
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

  T%N  = N
  T%Nm = N - 1

  allocate(T%valid(N))
  allocate(T%date(N))
  allocate(T%x(N))
  allocate(T%y(N))
  allocate(T%z(N))
  allocate(T%t(N))
  allocate(T%dt(N))
  allocate(T%uf(N))
  allocate(T%vf(N))
  allocate(T%wf(N))
  allocate(T%uo(N))
  allocate(T%vo(N))
  allocate(T%wo(N))
  allocate(T%speed(N))
  allocate(T%OU(N))
  allocate(T%OV(N))
  allocate(T%OW(N))
  allocate(T%AU(N))
  allocate(T%AV(N))

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

    N = T%N
    Nv = count(T%valid)

    call Tv%allocate(Nv)

    iv = 0
    do i=1,N
      if (T%valid(i)) then
        iv = iv + 1
        Tv%valid(iv) = T%valid(i)
        Tv%date(iv)  = T%date(i)
        Tv%x(iv)     = T%x(i)
        Tv%y(iv)     = T%y(i)
        Tv%z(iv)     = T%z(i)
        Tv%t(iv)     = T%t(i)
        Tv%dt(iv)    = T%dt(i)
        Tv%uf(iv)    = T%uf(i)
        Tv%vf(iv)    = T%vf(i)
        Tv%wf(iv)    = T%wf(i)
        Tv%uo(iv)    = T%uo(i)
        Tv%vo(iv)    = T%vo(i)
        Tv%wo(iv)    = T%wo(i)
        Tv%OU(iv)    = T%OU(i)
        Tv%OV(iv)    = T%OV(i)
        Tv%OW(iv)    = T%OW(i)
        Tv%AU(iv)    = T%AU(i)
        Tv%AV(iv)    = T%AV(i)
      endif
    enddo
    if (iv.ne.Nv) stop 'Incompatible iv and Nv in trajectory_compress'

    call T%allocate(Nv)

    T%valid(:) = Tv%valid(:)
    T%date(:)  = Tv%date(:)
    T%x(:)     = Tv%x(:)
    T%y(:)     = Tv%y(:)
    T%z(:)     = Tv%z(:)
    T%t(:)     = Tv%t(:)
    T%dt(:)    = Tv%dt(:)
    T%uf(:)    = Tv%uf(:)
    T%vf(:)    = Tv%vf(:)
    T%wf(:)    = Tv%wf(:)
    T%uo(:)    = Tv%uo(:)
    T%vo(:)    = Tv%vo(:)
    T%wo(:)    = Tv%wo(:)
    T%OU(:)    = Tv%OU(:)
    T%OV(:)    = Tv%OV(:)
    T%OW(:)    = Tv%OW(:)
    T%AU(:)    = Tv%AU(:)
    T%AV(:)    = Tv%AV(:)
 
    Tv%valid => NULL()
    Tv%date  => NULL()
    Tv%x     => NULL()
    Tv%y     => NULL()
    Tv%z     => NULL()
    Tv%t     => NULL()
    Tv%dt    => NULL()
    Tv%uf    => NULL()
    Tv%vf    => NULL()
    Tv%wf    => NULL()
    Tv%uo    => NULL()
    Tv%vo    => NULL()
    Tv%wo    => NULL()
    Tv%speed => NULL()
    Tv%OU    => NULL()
    Tv%OV    => NULL()
    Tv%OW    => NULL()
    Tv%AU    => NULL()
    Tv%AV    => NULL()

  end subroutine trajectory_compress
  ! ...
  ! =====================================================================
  ! ...
  subroutine trajectory_velocity(T)

    class(type_trajectory), intent(inout)          :: T

    ! ... Local variables
    ! ... 
    integer i,iv,N,Nv
    real(dp) dt,ym,dx,dy,dz,dtmean
    real(dp) u,v,w

    T%dt(:) = 0.0D0
    T%uo(:) = 0.0D0
    T%vo(:) = 0.0D0
    T%wo(:) = 0.0D0
    dtmean  = 0.0D0
    nv      = 0

    do i=1,T%N-1

      dx = T%x(i+1) - T%x(i)
      dy = T%y(i+1) - T%y(i)
      dz = T%z(i+1) - T%z(i)
      dt = T%t(i+1) - T%t(i)

      if (abs(dt).gt.0.and.abs(dx).gt.0.and.abs(dy).gt.0) then

        T%dt(i) = dt
        dtmean = dtmean + dt
        nv = nv + 1

        if (Spherical) then
          ym =  0.5D0*(T%y(i+1)+T%y(i))
          dx = REarth*dx*cos(ym)
          dy = REarth*dy
        endif
        u = dx/dt
        v = dy/dt
        w = dz/dt
        T%uo(i)     = u
        T%vo(i)     = v
        T%wo(i)     = w
        T%speed(i) = sqrt(u*u + v*v + w*w)
        if ((T%speed(i).lt.fit_vmin).or.(T%speed(i).gt.fit_vmax)) T%valid(i) = .False.

      else

        T%valid(i) = .False.

      endif

    enddo
    if (nv.eq.0) call crash ('No valid trajectory data')
    dtmean   = dtmean / nv
    T%dtmean = dtmean


  end subroutine trajectory_velocity
  ! ...
  ! =====================================================================
  ! ...
  real(dp) function func(p)

    real(dp), dimension(:), intent(in)       :: p

    ! ... Local variables
    ! ...
    integer f,l,n
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
      do l=1,To(f)%Nm
        if (To(f)%valid(l)) then
          u = p(1)*To(f)%OU(l) + p(2)*To(f)%AU(l) + p(3)*To(f)%AV(l)
          v = p(1)*To(f)%OV(l) + p(4)*To(f)%AU(l) + p(5)*To(f)%AV(l)
          du = To(f)%uo(l) - u
          dv = To(f)%vo(l) - v
          func = func + du*du + dv*dv
          n = n + 1
        endif
      enddo
    enddo
    func = 0.5D0 * func / n

  end function func
  ! ...
  ! ===================================================================
  ! ...
  function dfunc(p) result(xi)

    real(dp), dimension(:), intent(in)       :: p
    real(dp), dimension(size(p))             :: xi

    ! ... Local variables
    ! ...
    integer f,l,n
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
      do l=1,To(f)%Nm
        if (To(f)%valid(l)) then
          u = p(1)*To(f)%OU(l) + p(2)*To(f)%AU(l) + p(3)*To(f)%AV(l)
          v = p(1)*To(f)%OV(l) + p(4)*To(f)%AU(l) + p(5)*To(f)%AV(l)
          du = To(f)%uo(l) - u
          dv = To(f)%vo(l) - v
          xi(1) = xi(1) + du*To(f)%OU(l) + dv*To(f)%OV(l)
          xi(2) = xi(2) + du*To(f)%AU(l)
          xi(3) = xi(3) + du*To(f)%AV(l)
          xi(4) = xi(4) + dv*To(f)%AU(l)
          xi(5) = xi(5) + dv*To(f)%AV(l)
          n = n + 1
        endif
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
    integer f,l,ii
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
        do l=1,To(f)%Nm
          if (To(f)%valid(l)) then
            um = p(1)*To(f)%OU(l) + p(2)*To(f)%AU(l) + p(3)*To(f)%AV(l)
            vm = p(1)*To(f)%OV(l) + p(4)*To(f)%AU(l) + p(5)*To(f)%AV(l)
            ii = ii + 1
            Fvec(ii) = To(f)%uo(l) - um
            ii = ii + 1
            Fvec(ii) = To(f)%vo(l) - vm
          endif
        enddo
      enddo
      if  (ii.ne.M) call crash('Error in number of matchings')
      return
    endif

    if (iflag.eq.2) then
      ii = 0
      do f=1,size(To)
        do l=1,To(f)%Nm
          if (To(f)%valid(l)) then
            ii = ii + 1
            FJAC(ii,1) = -To(f)%OU(l) * fit_dop1
            FJAC(ii,2) = -To(f)%AU(l) * fit_dop2
            FJAC(ii,3) = -To(f)%AV(l) * fit_dop3
            FJAC(ii,4) =  0.0D0
            FJAC(ii,5) =  0.0D0
            ii = ii + 1
            FJAC(ii,1) = -To(f)%OV(l) * fit_dop1
            FJAC(ii,2) =  0.0D0
            FJAC(ii,3) =  0.0D0
            FJAC(ii,4) = -To(f)%AU(l) * fit_dop4
            FJAC(ii,5) = -To(f)%AV(l) * fit_dop5
          endif
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
