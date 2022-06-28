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

module module_model

use module_types
use module_constants
use module_math
use module_time
use module_grid

use module_alm
use module_forcing
use module_float

implicit none

character(len=*), parameter                      :: VERSION = '1.0'
character(len=*), parameter                      :: AUTHOR  = 'Joaquim Ballabrera'

! ... Runge-Kutta order
! ...
integer                                          :: rk_order = 5

integer                                          :: model_step  = 0
integer                                          :: model_nsteps= 0
real(dp)                                         :: model_tini  = 0.0D0
real(dp)                                         :: model_tfin  = 0.0D0
real(dp)                                         :: model_tlen  = 0.0D0
real(dp)                                         :: model_time  = 0.0D0
real(dp)                                         :: model_dt    = 3600.0D0
real(dp)                                         :: model_sign  = 1.0D0
real(dp)                                         :: model_velmin = 1D-5  ! m/s ! Pereiro, 2019
type(type_date)                                  :: model_dini, model_dfin
type(type_date)                                  :: model_date 

real(dp), dimension(:,:,:,:), pointer            :: OUrhs        ! (X,Y,Z,2)
real(dp), dimension(:,:,:,:), pointer            :: OVrhs        ! (X,Y,Z,2)
real(dp), dimension(:,:,:,:), pointer            :: OWrhs        ! (X,Y,Z,2)
real(dp), dimension(:,:,:,:), pointer            :: OTrhs        ! (X,Y,Z,2)
real(dp), dimension(:,:,:,:), pointer            :: OSrhs        ! (X,Y,Z,2)
real(dp), dimension(:,:,:,:), pointer            :: ORrhs        ! (X,Y,Z,2)
real(dp), dimension(:,:,:,:), pointer            :: OCrhs        ! (X,Y,Z,2)
real(dp), dimension(:,:,:), pointer              :: AUrhs        ! (X,Y,2)
real(dp), dimension(:,:,:), pointer              :: AVrhs        ! (X,Y,2)

! ... Time records to be read
! ...
integer                                          :: OUr1=-1, OUr2=-1
integer                                          :: OVr1=-1, OVr2=-1
integer                                          :: OWr1=-1, OWr2=-1
integer                                          :: OTr1=-1, OTr2=-1
integer                                          :: OSr1=-1, OSr2=-1
integer                                          :: ORr1=-1, ORr2=-1
integer                                          :: OCr1=-1, OCr2=-1
integer                                          :: AUr1=-1, AUr2=-1
integer                                          :: AVr1=-1, AVr2=-1

! ... Vertical records to be read
! ...
integer                                          :: OUk1=-1, OUk2=-1
integer                                          :: OVk1=-1, OVk2=-1
integer                                          :: OWk1=-1, OWk2=-1
integer                                          :: OTk1=-1, OTk2=-1
integer                                          :: OSk1=-1, OSk2=-1
integer                                          :: ORk1=-1, ORk2=-1
integer                                          :: OCk1=-1, OCk2=-1
integer                                          :: AUk1=-1, AUk2=-1
integer                                          :: AVk1=-1, AVk2=-1

integer                                          :: model_nlayers = 0
integer, dimension(:), pointer                   :: model_layer

integer                                          :: ADVECTION_LAYER = 0

! ... Response Matrix, accounting for direct atmosphere
! ... drag and Stokes drift.
! ... Pereiro et al. (2019).
! ... The atmospheric term only applies near the surface, thus
! ... the use of the surface flag.
! ...
logical                                          :: surface = .True.
real(dp)                                         :: A11     = 0.0_dp
real(dp)                                         :: A12     = 0.0_dp
real(dp)                                         :: A21     = 0.0_dp
real(dp)                                         :: A22     = 0.0_dp

! ... Water speed fraction (alpha)
! ... Dragg-related term. It also allows simulating when
! ... drifter is driven only by the wind (alpha = 0.0).
! ...
real(dp), dimension(ndims)                       :: alpha   = 1.0_dp  ! No units

! ... Time calendar
! ... The model will assume the time units and calendar of the zonal ocean
! ... current's file.
character(len=180)                               :: model_time_units = ''
character(len=20)                                :: model_time_calendar = ''

! ... Middle of the month day (for climatological forcing)
! ...
real(dp), dimension(12) :: HalfMonth = [16.5D0,15.0D0,16.5D0,16.0D0,16.5D0,16.0D0, &
                                        16.5D0,16.5D0,16.0D0,16.5D0,16.0D0,16.5D0]

! ... Random Walk subgrid parameterization
! ...
logical                                          :: noise_model_0 = .False.
logical                                          :: noise_model_1 = .False.
real(dp)                                         :: noise_mu
real(dp)                                         :: noise_KH0 = 0.0D0
real(dp)                                         :: noise_KV0 = 0.0D0
real(dp)                                         :: noise_V0  = 0.0D0
real(dp)                                         :: noise_W0  = 0.0D0
real(dp)                                         :: noise_KH1 = 0.0D0
real(dp)                                         :: noise_KV1 = 0.0D0
real(dp)                                         :: noise_V1  = 0.0D0
real(dp)                                         :: noise_W1  = 0.0D0

! ... Runge-Kutta interpolation point brackets:
! ... [ west, south, bottom ]
! ...
integer, dimension(3)                            :: IOU,IOV,IOW,IAU,IAV
integer, dimension(3)                            :: IOT,IOS,IOR

! ... Saving frequency
! ...
integer                                          :: save_period = 3600  ! seconds
integer                                          :: save_frequency

! ... Output files
! ...
character(len=maxlen)                            :: trajectory_name = 'alm.nc'  ! Traj
character(len=maxlen)                            :: trajectory_final = 'end.dat'  ! Traj final
character(len=maxlen)                            :: output_units=''
character(len=maxlen)                            :: output_calendar=''
integer                                          :: output_nfloats
integer                                          :: output_record = 0
integer                                          :: output_format = 0 !0=nc, 1=ascii
integer                                          :: output_id
integer                                          :: output_nid
integer                                          :: output_idl
integer                                          :: output_timeid
integer                                          :: output_lonid
integer                                          :: output_latid
integer                                          :: output_zid
integer                                          :: output_uid
integer                                          :: output_vid
integer                                          :: output_wid
integer                                          :: output_did
integer                                          :: output_rid
integer                                          :: output_xoid
integer                                          :: output_yoid
integer                                          :: output_zoid
integer                                          :: output_toid
integer                                          :: output_status
integer                                          :: output_kid
real(dp)                                         :: output_missing = -999.0D0


! ... Analytical density
! ... water_visc: Kinematic viscosity [m2/s] = Dynamic viscosity / rho
! ... water_rho: Water density kg/m3
! ...
logical                                          :: model_buoyancy = .False.
integer                                          :: water_density_method=-1
real(dp)                                         :: water_visc = 1D-6 ! m2/s
real(dp)                                         :: water_rho  = 1020.0D0 !kg/m3

contains
! ...
! =====================================================================
! =====================================================================
! ...
  function RHS(xo,t) result(dxdt)
  ! ...
  ! ... Spatio-temporal interpolation of the various terms of the 
  ! ... right-hand-side member of the equation
  ! ...               dx/dt = rhs(x,t)
  ! ...
  ! ... First thing in the Runge-Kutta subroutine, the 3D bracketing of
  ! ... the position is being done. Therefore, the following indices
  ! ... are supposed to be set: 
  ! ...            LOU, LOV, LOW, LAU, LAV
  ! ... where,
  ! ... IOU = [ ii, jj, kk] are the "left, lower, and top" indexes
  ! ... bracketing the integration point.
  ! ... etc.
  ! ...

    real(dp), dimension(ndims), intent(in)         :: xo
    real(dp), intent(in)                           :: t
    real(dp), dimension(ndims)                     :: dxdt

    ! ... Local variables
    ! ...
    real(dp) f1,f2,ff,a1,a2
    real(dp) t1,t2
    real(dp) veps,vx,vy,vz,rhok,wok

    if (noise_model_1) then
      ! ... Horizontal noise
      veps = randn()
      vx = noise_V1*veps*sin(dpi*veps)
      vy = noise_V1*veps*cos(dpi*veps)
      !veps = randn()
      !vx = noise_V1*veps
      !veps = randn()
      !vy = noise_V1*veps
      ! ... Vertical noise
      veps = randn()
      vz = noise_W1*veps
    else
      vx = 0.0D0
      vy = 0.0D0
      vz = 0.0D0
    endif

    dxdt(:) = 0.0D0

    ! ... Zonal advection:
    ! ...
    if (Uadv) then
      if (alpha(1).eq.0) then
        dxdt(1) = vx
      else
        IOU = GOU%locate(xo)
        if (GOU%Stationary) then
          ff  = GOU%interpol(OUrhs(:,:,:,1),xo,IOU,SingleLayer)
        else
          t1 = GOU%t(OUr1)
          t2 = GOU%t(OUr2)
          f1 = GOU%interpol(OUrhs(:,:,:,1),xo,IOU,SingleLayer)
          f2 = GOU%interpol(OUrhs(:,:,:,2),xo,IOU,SingleLayer)
          ff = f1 + (f2-f1)*(model_time-t1)/(t2-t1)
        endif
        dxdt(1) = alpha(1)*ff + vx
      endif
    else
      dxdt(1) = vx
    endif

    ! ... Meridional advection
    ! ...
    if (Vadv) then
      if (alpha(2).eq.0) then
        dxdt(2) = vy
      else
        IOV = GOV%locate(xo)
        if (GOV%Stationary) then
          ff  = GOV%interpol(OVrhs(:,:,:,1),xo,IOV,SingleLayer)
        else
          t1 = GOV%t(OVr1)
          t2 = GOV%t(OVr2)
          f1 = GOV%interpol(OVrhs(:,:,:,1),xo,IOV,SingleLayer)
          f2 = GOV%interpol(OVrhs(:,:,:,2),xo,IOV,SingleLayer)
          ff = f1 + (f2-f1)*(model_time-t1)/(t2-t1)
          !print*, 'IOV: ', IOV
          !print*, GOV%x1(IOV(1)), GOV%y1(IOV(2))
          !print*, 'f1, f2, ff : ', f1, f2, ff
        endif
        dxdt(2) = alpha(2)*ff + vy
      endif
    else
      dxdt(2) = vy
    endif

    ! ... Vertical advection
    ! ...
    if (Wadv) then
      if (alpha(3).eq.0) then
        dxdt(3) = vz
      else
        IOW = GOW%locate(xo)
        if (GOW%Stationary) then
          ff  = GOW%interpol(OWrhs(:,:,:,1),xo,IOW,SingleLayer)
        else
          t1 = GOW%t(OWr1)
          t2 = GOW%t(OWr2)
          f1 = GOW%interpol(OWrhs(:,:,:,1),xo,IOW,SingleLayer)
          f2 = GOW%interpol(OWrhs(:,:,:,2),xo,IOW,SingleLayer)
          ff = f1 + (f2-f1)*(model_time-t1)/(t2-t1)
        endif
        dxdt(3) = alpha(3)*ff + vz
      endif
    else
      dxdt(3) = vz
    endif

    ! ... Add, if requested, a buoyancy term:
    ! ...
    if (model_buoyancy) then
      select case (water_density_method)
      case (0)
        ! ... Density is constant
        ! ... 
        rhok = water_rho
      case (1)
        ! ... Density is a function of z and latitude
        ! ... 
        rhok = analytical_rho(xo(3),rad2deg*xo(2))
      case (2)
        ! ... Density comes from a file
        ! ... 
        IOR = GOR%locate(xo)
        if (GOR%Stationary) then
          ff  = GOR%interpol(ORrhs(:,:,:,1),xo,IOR,SingleLayer)
        else
          t1 = GOR%t(ORr1)
          t2 = GOR%t(ORr2)
          f1 = GOR%interpol(ORrhs(:,:,:,1),xo,IOR,SingleLayer)
          f2 = GOR%interpol(ORrhs(:,:,:,2),xo,IOR,SingleLayer)
          rhok = f1 + (f2-f1)*(model_time-t1)/(t2-t1)
        endif
      end select
      wok = emergence(FLTk,rhok)
      dxdt(3) = dxdt(3) + wok
    endif

    if (winds.and.surface) then
      ! Time interpolation of winds required !
      !a1     = GAU%hinterpol(AUrhs(:,:,ll),x(1),x(2))
      !a2     = GAV%hinterpol(AVrhs(:,:,ll),x(1),x(2))
      dxdt(1) = dxdt(1) +  A11*a1 + A12*a2
      dxdt(2) = dxdt(2) +  A21*a1 + A22*a2
    endif

  end function RHS
  ! ...
  ! ===================================================================
  ! ...
  subroutine model_ini()

    type(type_date)                         :: datemin,datemax
    integer i,j,k
    real(dp) critical_size,yave

    if (reverse) model_sign = -1.0D0

    model_tini = alm_tini
    model_tfin = alm_tfin
    model_tlen = nint(model_tfin - model_tini)
    if (mod(model_tlen,model_dt).ne.0) call crash('Simulation length not a multiple of time step.')
    model_nsteps = nint(model_tlen/model_dt)

    model_dini = num2date(model_tini,units=model_time_units,calendar=model_time_calendar)
    model_dfin = num2date(model_tfin,units=model_time_units,calendar=model_time_calendar)

    !if (minval([GOU%nt,GOV%nt,GOW%nt]).eq.1) Stationary = .True.

    if (WithOU) then
      Uadv = .True.
      allocate (OUrhs(GOU%nx,GOU%ny,GOU%nz,2))
      OUr1 = locate(GOU%t,model_tini) 
      if ((OUr1.lt.1.or.OUr1.ge.GOU%nt)) call crash ('U time not bracketed')
      if (GOU%nt.eq.1.or.Stationary) then
        OUr2 = OUr1
      else
        OUr2 = min(OUr1+1,GOU%nt)
      endif
      OUrhs(:,:,:,1) = GOU%read3D(GOU%varid,step=Our1,missing_value=0.0D0,verbose=verb.ge.3)
      OUrhs(:,:,:,2) = GOU%read3D(GOU%varid,step=Our2,missing_value=0.0D0,verbose=verb.ge.3)
    endif

    if (WithOV) then
      Vadv = .True.
      allocate (OVrhs(GOV%nx,GOV%ny,GOV%nz,2))
      OVr1 = locate(GOV%t,model_tini) 
      if ((OVr1.lt.1.or.OVr1.ge.GOV%nt)) call crash ('V time not bracketed')
      if (GOV%nt.eq.1.or.Stationary) then
        OVr2 = OVr1
      else
        OVr2 = min(OVr1+1,GOV%nt)
      endif
      OVrhs(:,:,:,1) = GOV%read3D(GOV%varid,step=Ovr1,missing_value=0.0D0,verbose=verb.ge.3)
      OVrhs(:,:,:,2) = GOV%read3D(GOV%varid,step=Ovr2,missing_value=0.0D0,verbose=verb.ge.3)
    endif

    if (WithOW) then
      Wadv = .True.
      allocate (OWrhs(GOW%nx,GOW%ny,GOW%nz,2))
      OWr1 = locate(GOW%t,model_tini) 
      if ((OWr1.lt.1.or.OWr1.ge.GOW%nt)) call crash ('W time not bracketed')
      if (GOW%nt.eq.1.or.Stationary) then
        OWr2 = OWr1
      else
        OWr2 = min(OWr1+1,GOW%nt)
      endif
      OWrhs(:,:,:,1) = GOW%read3D(GOW%varid,step=Owr1,missing_value=0.0D0,verbose=verb.ge.3)
      OWrhs(:,:,:,2) = GOW%read3D(GOW%varid,step=Owr2,missing_value=0.0D0,verbose=verb.ge.3)
    endif

    if (WithOT) then
      allocate (OTrhs(GOT%nx,GOT%ny,GOT%nz,2))
      OTr1 = locate(GOT%t,model_tini) 
      if ((OTr1.lt.1.or.OTr1.ge.GOT%nt)) call crash ('T time not bracketed')
      if (GOT%nt.eq.1.or.Stationary) then
        OTr2 = OTr1
      else
        OTr2 = min(OTr1+1,GOT%nt)
      endif
      OTrhs(:,:,:,1) = GOT%read3D(GOT%varid,step=Otr1,verbose=verb.ge.3)
      OTrhs(:,:,:,2) = GOT%read3D(GOT%varid,step=Otr2,verbose=verb.ge.3)
    endif

    if (WithOS) then
      allocate (OSrhs(GOS%nx,GOS%ny,GOS%nz,2))
      OSr1 = locate(GOS%t,model_tini) 
      if ((OSr1.lt.1.or.OSr1.ge.GOS%nt)) call crash ('S time not bracketed')
      if (GOS%nt.eq.1.or.Stationary) then
        OSr2 = OSr1
      else
        OSr2 = min(OSr1+1,GOS%nt)
      endif
      OSrhs(:,:,:,1) = GOS%read3D(GOS%varid,step=Osr1,verbose=verb.ge.3)
      OSrhs(:,:,:,2) = GOS%read3D(GOS%varid,step=Osr2,verbose=verb.ge.3)
    endif

    if (WithOR) then
      allocate (ORrhs(GOR%nx,GOR%ny,GOR%nz,2))
      ORr1 = locate(GOR%t,model_tini) 
      if ((ORr1.lt.1.or.ORr1.ge.GOR%nt)) call crash ('R time not bracketed')
      if (GOR%nt.eq.1.or.Stationary) then
        ORr2 = ORr1
      else
        ORr2 = min(ORr1+1,GOR%nt)
      endif
      ORrhs(:,:,:,1) = GOR%read3D(GOR%varid,step=Orr1,verbose=verb.ge.3)
      ORrhs(:,:,:,2) = GOR%read3D(GOR%varid,step=Orr2,verbose=verb.ge.3)
    endif

    if (WithOC) then
      allocate (OCrhs(GOC%nx,GOC%ny,GOC%nz,2))
      OCr1 = locate(GOC%t,model_tini) 
      if ((OCr1.lt.1.or.OCr1.ge.GOC%nt)) call crash ('C time not bracketed')
      if (GOC%nt.eq.1.or.Stationary) then
        OCr2 = OCr1
      else
        OCr2 = min(OCr1+1,GOC%nt)
      endif
      OCrhs(:,:,:,1) = GOC%read3D(GOC%varid,step=Ocr1,verbose=verb.ge.3)
      OCrhs(:,:,:,2) = GOC%read3D(GOC%varid,step=Ocr2,verbose=verb.ge.3)
    endif

    if (WithAU) then
      allocate (AUrhs(GAU%nx,GAU%ny,2))
      AUr1 = locate(GAU%t,model_tini) 
      if ((AUr1.lt.1.or.AUr1.ge.GAU%nt)) call crash ('AU time not bracketed')
      if (GAU%nt.eq.1.or.Stationary) then
        AUr2 = AUr1
      else
        AUr2 = min(AUr1+1,GAU%nt)
      endif
      AUrhs(:,:,1) = GAU%read2D(GAU%varid,step=Aur1,missing_value=0.0D0)
      AUrhs(:,:,2) = GAU%read2D(GAU%varid,step=Aur2,missing_value=0.0D0)
    endif

    if (WithAV) then
      allocate (AVrhs(GAV%nx,GAV%ny,2))
      AVr1 = locate(GAV%t,model_tini) 
      if ((AVr1.lt.1.or.AVr1.ge.GAV%nt)) call crash ('AV time not bracketed')
      if (GAV%nt.eq.1.or.Stationary) then
        AVr2 = AVr1
      else
        AVr2 = min(AVr1+1,GAV%nt)
      endif
      AVrhs(:,:,1) = GAV%read2D(GAV%varid,step=Avr1,missing_value=0.0D0)
      AVrhs(:,:,2) = GAV%read2D(GAV%varid,step=Avr2,missing_value=0.0D0)
    endif

    if (SingleLayer) Wadv = .False.

    save_frequency = nint(save_period/model_dt)

    noise_KH0 = abs(noise_KH0)
    noise_KV0 = abs(noise_KV0)
    noise_KH1 = abs(noise_KH1)
    noise_KV1 = abs(noise_KV1)

    noise_V0 = sqrt(2.0D0*noise_KH0/abs(model_dt))   ! 
    noise_V1 = sqrt(2.0D0*noise_KH1/abs(model_dt))

!   if (Wadv) then
      noise_W0 = sqrt(2.0D0*noise_KV0/abs(model_dt))   ! 
      noise_W1 = sqrt(2.0D0*noise_KV1/abs(model_dt))
!   else
!     noise_W0 = 0.0D0
!     noise_W1 = 0.0D0
!   endif

    if (noise_V0+noise_W0.gt.0.0D0) noise_model_0 = .True.
    if (noise_V1+noise_W1.gt.0.0D0) noise_model_1 = .True.
   
    ! ... Get the "typical" water density.
    ! ... if water_density_method .le. 0, then it is the
    ! ... value in water_rho.
    ! ..
    if (model_buoyancy) then
      if (water_density_method.eq.1) then
        yave = 0.5D0*(alm_ymin+alm_ymax)*rad2deg
        water_rho = 0.0D0
        do k=1,GOU%Nz
          water_rho = water_rho + analytical_rho(GOU%z(k),yave)
        enddo
        water_rho = water_rho / GOU%Nz
      else if (water_density_method.eq.2) then
        print*, 'WARNING: This needs to be calculated !!!!'
        yave = 0.5D0*(alm_ymin+alm_ymax)*rad2deg
        water_rho = 0.0D0
        do k=1,GOU%Nz
          water_rho = water_rho + analytical_rho(GOU%z(k),yave)
        enddo
        water_rho = water_rho / GOU%Nz
      endif
    endif

    if (verb.ge.1) then
      write(*,*)
      write(*,*) 'Simulation period: '
      if (model_sign.gt.0) then
        write(*,*) 'Forward model'
      else
        write(*,*) 'Backward model'
      endif
      write(*,*) 'Initial time  : ', model_tini, model_dini%iso()
      write(*,*) 'Final time    : ', model_tfin, model_dfin%iso()
      write(*,*) 'Time step     : ', model_dt
      write(*,*) 'Number steps  : ', model_nsteps
      write(*,*) 'Saving period : ', save_period
      write(*,*) 'Advection XYZ : ', Uadv, Vadv, Wadv

      if (model_buoyancy) then
        write(*,*) 'Model buoyancy activated'
        write(*,*) 'Water kinematic viscosity : ', water_visc
        write(*,*) 'Water density method      : ', water_density_method
        write(*,*) 'Reference Water density   : ', water_rho
        write(*,*) 'Particle density          : ', Release_rho
        write(*,*) 'Particle diameter         : ', Release_size
        if (Release_rho.lt.water_rho) then
          ! ... Maximum drop sizes are claculated by Aravamudan et al. (1982)
          ! ... Break up for oil on rough seas - simplified models and step-by
          ! ... step calculations. US Coast Guard Report CG-D-28-82
          ! ... US Department of Transportation, Washington DC.
          ! ...
          critical_size = 9.52D0*(water_visc/gravity)**(2.D0/3.D0)/(1.0D0-Release_rho/water_rho)**(1.0D0/3.0D0)
          write(*,*) 'Maximum Particle diameter : ', critical_size
          !if (Release_size.gt.critical_size) call crash('Particle size > Minumum allowed size')
        endif
      else
        write(*,*) 'Model buoyancy not activated'
      endif

      write(*,*) 'Wind forcing  : ', Winds
      write(*,*) 'Noise model 0 : ', noise_model_0
      write(*,*) 'Noise model 1 : ', noise_model_1
      if (Uadv) write(*,*) 'Zonal velocity record pointers     : ', OUr1, OUr2
      if (Vadv) write(*,*) 'Meridional velocity record pointers: ', OVr1, OVr2
      if (Wadv) write(*,*) 'Vertical velocity record pointers  : ', OWr1, OWr2
      if (WithAU) write(*,*) 'Zonal wind record pointers         : ', AUr1, AUr2
      if (WithAV) write(*,*) 'Meridional wind record pointers    : ', AVr1, AVr2
    endif

    return

  end subroutine model_ini
  ! ...
  ! ===================================================================
  ! ...
  subroutine model_run()

    ! ... Local variables
    ! ...
    integer i,j,k,ii,jj,kk,ll,ifloat,float_status
    real(dp) veps,dx,dy,dz
    real(dp), dimension(ndims)                   :: xo,xn,un


    ! ... First record:
    ! ...
    model_time = model_tini 
    do ifloat=1,Nfloats
      FLT(ifloat)%to = model_tini + FLT(ifloat)%to    ! FLT%to seconds since model_tini
      if (model_time.ge.FLT(ifloat)%to) then
        if (verb.ge.2) write(*,*) 'Initial step releasing for float ', ifloat
        FLT(ifloat)%released = .True.
        FLT(ifloat)%floating = .True.
        FLT(ifloat)%indomain = .True.
        FLT(ifloat)%x        = FLT(ifloat)%xo
        FLT(ifloat)%y        = FLT(ifloat)%yo
        FLT(ifloat)%z        = FLT(ifloat)%zo
        FLT(ifloat)%dist     = 0.0D0
        FLT(ifloat)%tlast    = model_tini
      endif
    enddo
    call trajectory_write(model_time,verb.ge.2) 


    do model_step=1,model_nsteps

      model_time = model_tini + (model_step-1)*model_dt
      model_date = num2date(Reference_time+model_time,model_time_units,model_time_calendar)

      if (verb.ge.1) write(*,*) 'step, time, date :: ', model_step, model_time, model_date%iso()

      do ifloat=1,Nfloats

        ! ... First chek if released or going to be released
        ! ...
        if (.NOT.FLT(ifloat)%released) then
          if (model_time.ge.FLT(ifloat)%to) then
            if (verb.ge.2) write(*,*) 'Releasing float ', ifloat
            FLT(ifloat)%released = .True.
            FLT(ifloat)%floating = .True.
            FLT(ifloat)%indomain = .True.
            FLT(ifloat)%x        = FLT(ifloat)%xo
            FLT(ifloat)%y        = FLT(ifloat)%yo
            FLT(ifloat)%z        = FLT(ifloat)%zo
            FLT(ifloat)%dist     = 0.0D0
            FLT(ifloat)%tlast    = model_time
          endif
        endif

        if (FLT(ifloat)%floating) then

          FLTk = FLT(ifloat)

          ! ... Runge - Kutta
          ! ....................................................
          xo = [FLT(ifloat)%x, FLT(ifloat)%y, FLT(ifloat)%z ]
          call rk(xo,model_time,model_dt,un,xn)
          ! ....................................................

          if (verb.ge.3) then
            write(*,*) 'Float: ', ifloat
            write(*,*) '   xo: ', xo
            write(*,*) '   un: ', un
            write(*,*) '   xn: ', xn
          endif

          FLT(ifloat)%x = xn(1)
          FLT(ifloat)%y = xn(2)
          FLT(ifloat)%z = xn(3)
          FLT(ifloat)%u = un(1)
          FLT(ifloat)%v = un(2)
          FLT(ifloat)%w = un(3)
          ! ...................

          dx = 0.001D0*Rearth*cos(xo(2))*(xn(1)-xo(1))
          dy = 0.001D0*Rearth*(xn(2)-xo(2))
          FLT(ifloat)%dist = FLT(ifloat)%dist + sqrt(dx*dx+dy*dy)

          !float_status =  check_float_status(xn)
          float_status = point_type(xn(1),xn(2),xn(3))  ! = 1, floating inside
                                                        ! = 0, stranded on land
                                                        ! = -1, outside the domain

          select case (float_status)
            case (1)
              FLT(ifloat)%tlast = model_time
            case (0) 
              if (verb.ge.2) write(*,*) 'Particle stranded on land: ', ifloat
              FLT(ifloat)%floating = .False.
            case (-1)
              if (verb.ge.2) write(*,*) 'Particle leaving the domain: ', ifloat
              FLT(ifloat)%indomain = .False.
              FLT(ifloat)%floating = .False.
          end select

        endif

        !print '(I5,F10.0,3F9.4)', model_step, model_time, xn
        !xo = xn
      enddo

      model_time = model_time + model_sign*model_dt

      if (mod(model_step,save_frequency).eq.0) then
        if (verb.ge.1) write(*,*) "Saving trajectory position"
        call trajectory_write(model_time,verb.ge.2)
      endif

      call forcing_update(model_time)

    enddo

  end subroutine model_run
  ! ...
  ! ===================================================================
  ! ...
  subroutine forcing_update(time)
  ! ... Check if a new time bracketed is required
  ! ...
    real(dp), intent(in)                         :: time

    if (WithOU) then
      ! ..................................................
      ! .................................................. OCE U
      ! ..................................................
      if (Our1.lt.0) then
        ! ... First call !
        ! ...
        OUr1 = locate(GOU%t,time) 
        if ((OUr1.lt.1.or.OUr1.ge.GOU%nt)) call crash ('U time not bracketed')
        if (GOU%nt.eq.1.or.Stationary) then
          OUr2 = OUr1
        else
          OUr2 = OUr1 + 1
        endif
        OUrhs(:,:,:,1) = GOU%read3D(GOU%varid,step=Our1,missing_value=0.0D0)
        OUrhs(:,:,:,2) = GOU%read3D(GOU%varid,step=Our2,missing_value=0.0D0)
      else
        ! ... Update fields
        ! ...
        if (reverse) then
          if (time.lt.GOU%t(OUr1)) then
            Our2 = Our1
            Our1 = Our1 - 1
            if (verb.ge.1) print '(T2,"Backward bracket U update:",F9.0," - ",F9.0)', &
                           GOU%t(Our1), GOU%t(Our2)
            OUrhs(:,:,:,2) = OUrhs(:,:,:,1)
            OUrhs(:,:,:,1) = GOU%read3D(GOU%varid,step=Our1,missing_value=0.0D0)
          endif
        else
          if (time.ge.GOU%t(OUr2)) then
            Our1 = Our2
            Our2 = Our2 + 1
            if (verb.ge.1) print '(T2,"Forward bracket U update:",F9.0," - ",F9.0)', &
                           GOU%t(Our1), GOU%t(Our2)
            OUrhs(:,:,:,1) = OUrhs(:,:,:,2)
            OUrhs(:,:,:,2) = GOU%read3D(GOU%varid,step=Our2,missing_value=0.0D0)
          endif
        endif

      endif
    endif


  end subroutine forcing_update
  ! ...
  ! ===================================================================
  ! ...
  function check_float_status(x) result(istatus)
  ! ... 
  ! ... function that returns:
  ! ...    0 if particle in water inside the domain
  ! ...    1 if particle in land  inside the domain
  ! ...    2 if particle outside the domain
  ! ...
    real(dp), dimension(ndims), intent(inout)    :: x
    integer                                      :: istatus

    istatus = 2
    if (x(1).le.forcing_xmin) return
    if (x(1).ge.forcing_xmax) return
    if (x(2).le.forcing_ymin) return
    if (x(2).ge.forcing_ymax) return
   
    istatus = 0 

    return

  end function check_float_status
  ! ...
  ! ===================================================================
  !                     TRAJECTORY FILE SUBROUTINES
  ! ===================================================================
  ! ...
  subroutine trajectory_create(filename,Nfloats,missing)

    character(len=*), intent(in)               :: filename
    integer, intent(in)                        :: Nfloats
    real(dp), intent(in)                       :: missing

    ! ... Local variables
    ! ...
    integer err,i,natts
    character(len=10) ext
    character(len=1000) lcom


    trajectory_name = trim(filename)
    output_nfloats = Nfloats
    output_missing = missing
    output_record = 0

    ! ... Check for extension
    ! ...
    do i=len(filename),1,-1
      if (filename(i:i).eq.'.') then
        ext = filename(i+1:)
        exit
      endif
    enddo

    if (lowercase(ext).eq.'nc') then
      output_format = 0                   ! Netcdf
    else
      output_format = 1                   ! Ascii
    endif

    if (output_format.eq.0) then
      ! ... Netcdf output format
      ! ...
      err = NF90_CREATE(filename,NF90_CLOBBER,output_id);             call check()
      err = NF90_DEF_DIM(output_id,'floats',Nfloats,output_nid);      call check()
      err = NF90_DEF_DIM(output_id,'time',NF90_UNLIMITED,output_idl); call check()

      err = NF90_DEF_VAR(output_id,'release_lon',NF90_DOUBLE,(/output_nid/),output_xoid)
      call check()
      err = NF90_PUT_ATT(output_id,output_xoid,'units','degrees_east')
      call check()

      err = NF90_DEF_VAR(output_id,'release_lat',NF90_DOUBLE,(/output_nid/),output_yoid)
      call check()
      err = NF90_PUT_ATT(output_id,output_yoid,'units','degrees_north')
      call check()

      err = NF90_DEF_VAR(output_id,'release_depth',NF90_DOUBLE,(/output_nid/),output_zoid)
      call check()
      err = NF90_PUT_ATT(output_id,output_zoid,'units','meters')
      call check()

      err = NF90_DEF_VAR(output_id,'release_time',NF90_DOUBLE,(/output_nid/),output_toid)
      call check()
      err = NF90_PUT_ATT(output_id,output_toid,'units','seconds since initial time')
      call check()

      err = NF90_DEF_VAR(output_id,'time',NF90_DOUBLE,(/output_idl/),output_timeid)
      call check()

      if (input_timeid.le.0) then
        call cdf_copyatts(.False.,input_id,input_timeid,output_id,output_timeid,natts)
      else
        err = NF90_PUT_ATT(output_id,output_timeid,'units',trim(model_time_units))
        call check()
        err = NF90_PUT_ATT(output_id,output_timeid,'calendar',trim(model_time_calendar))
        call check()
      endif

      err = NF90_DEF_VAR(output_id,'lon',NF90_REAL,(/output_nid,output_idl/),output_lonid)
      call check()
      err = NF90_PUT_ATT(output_id,output_lonid,'units','degrees_east')
      call check()
      err = NF90_PUT_ATT(output_id,output_lonid,'_FillValue',sngl(missing))
      call check()

      err = NF90_DEF_VAR(output_id,'lat',NF90_REAL,(/output_nid,output_idl/),output_latid)
      call check()
      err = NF90_PUT_ATT(output_id,output_latid,'units','degrees_north')
      call check()
      err = NF90_PUT_ATT(output_id,output_latid,'_FillValue',sngl(missing))
      call check()

      err = NF90_DEF_VAR(output_id,'depth',NF90_REAL,(/output_nid,output_idl/),output_zid)
      call check()
      err = NF90_PUT_ATT(output_id,output_zid,'units','meters')
      call check()
      err = NF90_PUT_ATT(output_id,output_zid,'_FillValue',sngl(missing))
      call check()

      err = NF90_DEF_VAR(output_id,'u',NF90_REAL,(/output_nid,output_idl/),output_uid)
      call check()
      err = NF90_PUT_ATT(output_id,output_uid,'units','m/s')
      call check()
      err = NF90_PUT_ATT(output_id,output_uid,'_FillValue',sngl(missing))
      call check()

      err = NF90_DEF_VAR(output_id,'v',NF90_REAL,(/output_nid,output_idl/),output_vid)
      call check()
      err = NF90_PUT_ATT(output_id,output_vid,'units','m/s')
      call check()
      err = NF90_PUT_ATT(output_id,output_vid,'_FillValue',sngl(missing))
      call check()

      err = NF90_DEF_VAR(output_id,'w',NF90_REAL,(/output_nid,output_idl/),output_wid)
      call check()
      err = NF90_PUT_ATT(output_id,output_wid,'units','m/s')
      call check()
      err = NF90_PUT_ATT(output_id,output_wid,'_FillValue',sngl(missing))
      call check()

      err = NF90_DEF_VAR(output_id,'distance',NF90_REAL,(/output_nid,output_idl/),output_did)
      call check()
      err = NF90_PUT_ATT(output_id,output_did,'units','km')
      call check()
      err = NF90_PUT_ATT(output_id,output_did,'_FillValue',sngl(missing))
      call check()

      err = NF90_DEF_VAR(output_id,'exitcode',NF90_INT,(/output_nid/),output_status)
      call check()
      err = NF90_PUT_ATT(output_id,output_status,'long_name','Float status at the end of the simulation')
      err = NF90_PUT_ATT(output_id,output_status,"is-1",'float has not been released')
      err = NF90_PUT_ATT(output_id,output_status,"is0",'float was moving')
      err = NF90_PUT_ATT(output_id,output_status,"is1",'float left the model area')
      err = NF90_PUT_ATT(output_id,output_status,"is2",'float was stranded')
      call check()

      if (model_buoyancy) then
        err = NF90_PUT_ATT(output_id,0,'Buoyancy','Activated')
        err = NF90_PUT_ATT(output_id,0,'Water_viscosity',Water_visc)
        err = NF90_PUT_ATT(output_id,0,'Water_density_method',water_density_method)
        err = NF90_PUT_ATT(output_id,0,'Reference_Water_density',water_rho)
        err = NF90_PUT_ATT(output_id,0,'Particle_size',Release_size)
        err = NF90_PUT_ATT(output_id,0,'Particle_density',Release_rho)
      else
        err = NF90_PUT_ATT(output_id,0,'Buoyancy','No activated')
      endif

      call get_commandline(lcom)
      err = NF90_PUT_ATT(output_id,0,'Command_line',TRIM(lcom))
      call check()

      err = NF90_ENDDEF(output_id)
      call check()

      ! ... Save release points
      ! ...
      err = NF90_PUT_VAR(output_id,output_xoid,rad2deg*FLT(:)%xo); call check()
      err = NF90_PUT_VAR(output_id,output_yoid,rad2deg*FLT(:)%yo); call check()
      err = NF90_PUT_VAR(output_id,output_zoid,FLT(:)%zo); call check()
      err = NF90_PUT_VAR(output_id,output_toid,FLT(:)%to); call check()

    else
      ! ... ASCII output format
      ! ...
      output_id = unitfree()
      open(output_id,file=filename,status='unknown')
      rewind(output_id)
      write(output_id,*) '# ALM ASCII trajectory'
      write(output_id,*) '# Nfloats '
      write(output_id,*) '# time '
      write(output_id,*) '# x(1), x(2),  ... x(Nfloats)'
      write(output_id,*) '# y(1), y(2),  ... y(Nfloats)'
      write(output_id,*) '# z(1), z(2),  ... z(Nfloats)'
      write(output_id,*) output_nfloats


    endif

    contains
      subroutine check()
        call cdf_error(err,'while creating output NetCDF trajectory file.')
      end subroutine check

  end subroutine trajectory_create
  ! ...
  ! ==================================================================
  ! ...
  subroutine trajectory_write(time,verbose)

    real(dp), intent(in)                         :: time
    logical, intent(in), optional                :: verbose

    ! ... Local variables
    ! ...
    logical verb
    integer i,err
    real(dp) time_out(1), xx(Nfloats)

    if (present(verbose)) then
      verb = verbose
    else
      verb = .False.
    endif

    output_record = output_record + 1

    if (verb) write(*,*) 'Saving record ', output_record, ' at time ', time

    ! ... Time
    ! ...
    time_out(1) = time
    err = NF90_PUT_VAR(output_id,output_timeid,time_out,[output_record],[1])
    call cdf_error(err,'Error writing trajectory time')


    ! ... Longitude
    ! ...
    xx(:) = output_missing
    do i=1,Nfloats
      if (FLT(i)%released) xx(i) = FLT(i)%x*rad2deg
    enddo

    if (output_format.eq.0) then
      err = NF90_PUT_VAR(output_id,output_lonid,xx,[1,output_record],[Nfloats,1])
      call cdf_error(err,'Error writing trajectory longitude')
    else if (output_format.eq.1) then
      write(output_id,*) xx
    endif

    ! ... Latitude
    ! ...
    xx(:) = output_missing
    do i=1,Nfloats
      if (FLT(i)%released) xx(i) = FLT(i)%y*rad2deg
    enddo

    if (output_format.eq.0) then
      err = NF90_PUT_VAR(output_id,output_latid,xx,[1,output_record],[Nfloats,1])
      call cdf_error(err,'Error writing trajectory latitude')
    else if (output_format.eq.1) then
      write(output_id,*) xx
    endif

    ! ... Depth
    ! ...
    xx(:) = output_missing
    do i=1,Nfloats
      if (FLT(i)%released) xx(i) = abs(FLT(i)%z)
    enddo

    if (output_format.eq.0) then
      err = NF90_PUT_VAR(output_id,output_zid,xx,[1,output_record],[Nfloats,1])
      call cdf_error(err,'Error writing trajectory depth')
    else if (output_format.eq.1) then
      write(output_id,*) xx
    endif

    if (output_record.gt.1) then

      ! ... U
      ! ...
      xx(:) = output_missing
      do i=1,Nfloats
        if (FLT(i)%released) xx(i) = FLT(i)%u
      enddo
  
      if (output_format.eq.0) then
        err = NF90_PUT_VAR(output_id,output_uid,xx,[1,output_record-1],[Nfloats,1])
        call cdf_error(err,'Error writing trajectory zonal velocity')
      else if (output_format.eq.1) then
        write(output_id,*) xx
      endif

      ! ... V
      ! ...
      xx(:) = output_missing
      do i=1,Nfloats
        if (FLT(i)%released) xx(i) = FLT(i)%v
      enddo
  
      if (output_format.eq.0) then
        err = NF90_PUT_VAR(output_id,output_vid,xx,[1,output_record-1],[Nfloats,1])
        call cdf_error(err,'Error writing trajectory meridional velocity')
      else if (output_format.eq.1) then
        write(output_id,*) xx
      endif

      ! ... W
      ! ...
      xx(:) = output_missing
      do i=1,Nfloats
        if (FLT(i)%released) xx(i) = FLT(i)%w
      enddo
  
      if (output_format.eq.0) then
        err = NF90_PUT_VAR(output_id,output_wid,xx,[1,output_record-1],[Nfloats,1])
        call cdf_error(err,'Error writing trajectory vertical velocity')
      else if (output_format.eq.1) then
        write(output_id,*) xx
      endif


    endif

    ! ... Distance
    ! ...
    xx(:) = output_missing
    do i=1,Nfloats
      if (FLT(i)%released) xx(i) = abs(FLT(i)%dist)
    enddo

    if (output_format.eq.0) then
      err = NF90_PUT_VAR(output_id,output_did,xx,[1,output_record],[Nfloats,1])
      call cdf_error(err,'Error writing trajectory distance')
    else if (output_format.eq.1) then
      write(output_id,*) xx
    endif


  end subroutine trajectory_write
  ! ...
  ! ==================================================================
  ! ...
  subroutine trajectory_close()

    ! ... Local variables
    ! ...
    integer err,i,code(Nfloats)
   
    do i=1,Nfloats
      if (FLT(i)%released) then
        if (FLT(i)%indomain) then
          if (FLT(i)%floating) then
            code(i) = 0
          else
            code(i) = 2
          endif
        else
          code(i) = 1
        endif
      else
        code(i) = -1
      endif
    enddo 

    if (output_format.eq.0) then
      err = NF90_PUT_VAR(output_id,output_status,code)
      call cdf_error(err,'Error writing output code')
      err = NF90_CLOSE(output_id)
      call cdf_error(err,'while closing output NetCDF trajectory file.')
      if (verb.ge.1) write(*,*) 'Successfully closed NetCDF file ',trim(trajectory_name)
    else
      close(output_id)
      if (verb.ge.1) write(*,*) 'Successfully closed ASCII file ',trim(trajectory_name)
    endif

  end subroutine trajectory_close
  ! ...
  ! ===================================================================
  !                     RUNGE - KUTTA SUBROUTINES
  ! ===================================================================
  ! ...
  subroutine rk(xo,to,dt,un,xn)
  ! ------------------------------------------------------------------
  ! ... Runge-Kutta algorithm to advance the position of the float
  ! ------------------------------------------------------------------
 
    real(dp), dimension(ndims), intent(in)     :: xo
    real(dp), intent(in)                       :: to
    real(dp), intent(in)                       :: dt
    real(dp), dimension(ndims), intent(out)    :: un
    real(dp), dimension(ndims), intent(out)    :: xn

    ! ... Local variables
    ! ...
    real(dp) dt12,dt13,dt14,dt16,dt23,dt34
    real(dp) speed,veps,vx,vy,vz
    real(dp), dimension(ndims)                 :: x2,x3,x4,x5,x6
    real(dp), dimension(ndims)                 :: u1,u2,u3,u4,u5,u6

    dt12 = 0.50D0*dt       ! 1*dt/2
    dt14 = 0.25D0*dt       ! 1*dt/4
    dt16 = dt/6.0D0        ! 1*dt/6
    dt34 = 0.75D0*dt       ! 3*dt/4
  
    if (rk_order.eq.3) then
  
      dt13 =       dt/3.0D0        ! 1*dt/3
      dt23 = 2.0D0*dt/3.0D0        ! 2*dt/3
  
      u1 = RHS(xo,to)
  
      x2 = displacement(xo,u1,dt13)
      u2 = RHS(x2,to+dt13)
  
      x3 = displacement(xo,u2,dt23)
      u3 = RHS(x3,to+dt23)
  
      un = 0.25D0*(u1 + 3*u3)
  
    else if (rk_order.eq.4) then

      u1 = RHS(xo,to)

      x2 = displacement(xo,u1,dt12)
      u2 = RHS(x2,to+dt12)

      x3 = displacement(xo,u2,dt12)
      u3 = RHS(x3,to+dt12)

      x4 = displacement(xo,u3,dt)
      u4 = RHS(x4,to+dt)

      un = (u1 + 2.0_dp*(u2+u3) + u4)/6.0_dp

    else if (rk_order.eq.5) then

      u1 = RHS(xo,to)
      !print*, 'u1 : ', u1

      x2 = displacement(xo,u1,dt14)
      u2 = RHS(x2,to+dt14)
      !print*, 'u2 : ', u2


      x3 = displacement(xo,0.5D0*(u1+u2),dt14)
      u3 = RHS(x3,to+dt14)
      !print*, 'u3 : ', u3

      x4 = displacement(xo,2.0D0*u3-u2,dt12)
      u4 = RHS(x4,to+dt12)
      !print*, 'u4 : ', u4

      x5 = displacement(xo,0.25D0*(u1+3.0_dp*u4),dt34)
      u5 = RHS(x5,to+dt34)
      !print*, 'u5 : ', u5

      x6 = displacement(xo,(-3.0_dp*u1+2.0_dp*u2+12.0_dp*u3-12.0_dp*u4+8.0_dp*u5)/7.0_dp,dt)
      u6 = RHS(x6,to+dt)
      !print*, 'u6 : ', u6

      un = (7.0_dp*u1 + 32.0_dp*u3 + 12.0_dp*u4 + 32.0_dp*u5 + 7.0_dp*u6)/90.0_dp
      
    else
      write(0,*) 'Invalid RK_ORDER in rk function'
      stop 1

    endif

    ! ... Minimum Speed Threshold
    ! ... Pereiro et al. (2019). J. Sea Res.,
    ! ... For computational reasons, the velocity of a beached drifter is small
    ! ... but always different from zero. A minimum threshold to prevent interpolation
    ! ... velocitites near the coast to drive a beached drifter.
    ! ...
    speed = sqrt(dot_product(un,un))
    if (speed.lt.model_velmin) then
      un(:) = 0.0D0
      xn(:) = xo(:)
    else

      if (noise_model_0) then
        ! ... Horizontal noise
        ! ... Sergei A. Lonin: Lagrangian Model for Oil Spill Diffusion at Sea
        ! ... Spill Science and Technology Bulletin, 5 (5/6), 331-336, 1999.
        ! ... The velocity fluctuations can be calculated based on a random
        ! ... walk technique 
        veps = randn()
        vx = noise_V0*veps*sin(dpi*veps)
        vy = noise_V0*veps*cos(dpi*veps)
        ! ... Vertical noise
        veps = randn()
        vz = noise_W0*veps
        un(1) = un(1) + vx
        un(2) = un(2) + vy
        un(3) = un(3) + vz
      endif
      xn = displacement(xo,un,dt)
    endif

  end subroutine rk
  ! ...
  ! =====================================================================
  ! ...
  function displacement(xo,u,dt) result (xn)

    real(dp), parameter                          :: Iearth = 1.0D0/Rearth

    real(dp), dimension(ndims), intent(in)       :: xo                 ! Input in radians
    real(dp), dimension(ndims), intent(in)       :: u
    real(dp), intent(in)                         :: dt
    real(dp), dimension(ndims)                   :: xn                 ! Output in radians

    real(dp) dx,dy,dz,ryo,coslat,ryn,rdx

    dx   = dt*u(1)     !   [seconds] x [meters/seconds]
    dy   = dt*u(2)     !   [seconds] x [meters/seconds]
    dz   = dt*u(3)     !   [seconds] x [meters/seconds]

    if (Spherical) then
      ! ... Spherical metrics
      ! ...
      ryo  = xo(2)
      coslat = cos(ryo)
      ryn = asin(sin(ryo+dy*IEarth)*cos(dx*IEarth))
      rdx = atan2(sin(dx*Iearth)*coslat,cos(dx*Iearth)-sin(ryo)*sin(ryn))
      xn(1)  = xo(1) + rdx             ! Radians
      xn(2)  = ryn                     ! Radians
      xn(3)  = xo(3) + dz              ! meters
    else
      ! ... Cartesian metrics
      ! ...
      xn(1)  = xo(1) + dx              ! meters
      xn(2)  = xo(2) + dy              ! meters
      xn(3)  = xo(3) + dz              ! meters
    endif

    if (xn(3).gt.0) xn(3) = 0.0D0

  end function displacement
  ! ...
  ! =====================================================================
  ! ...
  function analytical_rho(z,lat) result(rho)

    ! ... Analytical model for ocean density
    ! ... V. Gladkikh and R. Tenzer
    ! ... A mathematical model of the Global Ocean Saltwater Density
    ! ... distribution. Pure and Applied Geophysics, 169 (February) 2012
    ! ... 249-257
    ! ... DOI: 10.1007/s00024-011-0275-5
    ! ...
    ! ... rho(D,phi) = 1000 + A(phi) + B(phi)*D^C(phi)
    ! ...

    real(dp), intent(in)                         :: z     ! meters
    real(dp), intent(in)                         :: lat   ! degrees
    real(dp)                                     :: rho

    real(dp) phi,D,A,B,C
    real(dp) mu,Pycno

    phi = abs(lat)
    D   = abs(z)

    mu = 0.928D0 - 0.079D0*cos(0.053*lat)
    Pycno = mu + 0.5D0*(1.0D0-mu)*(1.0D0+tanh(0.00988D0*D-1.01613))

    A = 27.91D0 - 2.06D0*exp(-(0.0161D0*phi)**5.0D0)
    B = 0.00637D0 + 0.00828D0*exp(-(0.017D0*phi)**4.66D0)
    C = 0.964D0 - 0.091D0*exp(-(0.016*phi)**5.0D0)

    rho = 1000.0D0 + A*Pycno + B*D**C

    return

  end function analytical_rho
  ! ...
  ! =====================================================================
  ! ...
  function emergence(F,rw) result (wo)

    type(type_float), intent(in)                 :: F
    real(dp), intent(in)                         :: rw
    real(dp)                                     :: wo

    ! ... Stokes terminal velocity:
    ! ... Sphere of radius Radius_s and density rho_s
    ! ... In a fluid with density rho_m and viscosity eta_m
    ! ... Positive emergence, upward emrgence.
    ! ...
    ! ... wo = - 2*g*Radius_s**2*(rho_s-rho_m)/(9*eta_m)
    ! ...    = - g*Diameter_s**2*(rho_s-rho_m)/(18*eta_m)
    ! ...    = g*Diameter_s**2*(1-rho_s/rho_m)/(18*mu_m)
    ! ... with
    ! ... Diameter_s = Sphere diameter (size)
    ! ... mu_m/rho_m = Kinematic viscosity      ! m2 s-1
    ! ...
    wo = gravity*F%size*F%size*(1.0D0-F%rho/rw)/18.0D0/water_visc

  end function emergence
  ! ...
  ! =====================================================================
  ! ...
end module module_model
