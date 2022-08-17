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

! ... Runge-Kutta order
! ...
integer                                          :: rk_order = 5

integer                                          :: model_step  = 0
integer                                          :: model_nsteps= 0
real(dp)                                         :: model_tini  = 0.0D0
real(dp)                                         :: model_tfin  = 0.0D0
real(dp)                                         :: model_tlen  = 0.0D0
real(dp)                                         :: model_time  = 0.0D0
real(dp)                                         :: model_daysec= 0.0D0
real(dp)                                         :: model_dt    = 3600.0D0
real(dp)                                         :: model_sign  = 1.0D0
real(dp)                                         :: model_velmin = 1D-5  ! m/s ! Pereiro, 2019
type(type_date)                                  :: model_dini, model_dfin
type(type_date)                                  :: model_date 

real(dp), dimension(:,:,:,:), allocatable        :: OUrhs        ! (X,Y,Z,2)
real(dp), dimension(:,:,:,:), allocatable        :: OVrhs        ! (X,Y,Z,2)
real(dp), dimension(:,:,:,:), allocatable        :: OWrhs        ! (X,Y,Z,2)
real(dp), dimension(:,:,:,:), allocatable        :: OTrhs        ! (X,Y,Z,2)
real(dp), dimension(:,:,:,:), allocatable        :: OSrhs        ! (X,Y,Z,2)
real(dp), dimension(:,:,:,:), allocatable        :: ORrhs        ! (X,Y,Z,2)
real(dp), dimension(:,:,:,:), allocatable        :: OCrhs        ! (X,Y,Z,2)
real(dp), dimension(:,:,:,:), allocatable        :: AUrhs        ! (X,Y,1,2)
real(dp), dimension(:,:,:,:), allocatable        :: AVrhs        ! (X,Y,1,2)

! ... Constant values
! ...
logical                                          :: model_fixed_ou = .False.
logical                                          :: model_fixed_ov = .False.
logical                                          :: model_fixed_ow = .False.
logical                                          :: model_fixed_au = .False.
logical                                          :: model_fixed_av = .False.
real(dp)                                         :: model_value_ou = 0.0D0
real(dp)                                         :: model_value_ov = 0.0D0
real(dp)                                         :: model_value_ow = 0.0D0
real(dp)                                         :: model_value_au = 0.0D0
real(dp)                                         :: model_value_av = 0.0D0

! ... Climatology flags
! ...
logical                                          :: OUClim = .False.
logical                                          :: OVClim = .False.
logical                                          :: OWClim = .False.
logical                                          :: AUClim = .False.
logical                                          :: AVClim = .False.

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
! ... The atmospheric term only applies at depths shallower than WindDepth.
! ...
logical                                          :: WindResponse = .False.
real(dp)                                         :: A11       = 0.0_dp
real(dp)                                         :: A12       = 0.0_dp
real(dp)                                         :: A21       = 0.0_dp
real(dp)                                         :: A22       = 0.0_dp

! ... Wind-driven current regression to local winds
! ...       U_wd = beta * exp (alpha*j) * (Windx + j Windy)
! ... It has sense when superposed to a geostrophic ocean current.
! ... Poulain et al. 2009, JAOT, 26, 1144 - 1156
! ... DOI: 10.1175/2008JTECHO618.1
! ...
logical                                          :: WindDriven = .False.
real(dp)                                         :: WDriven_beta  = 0.0D0
real(dp)                                         :: WDriven_theta = 0.0D0

real(dp)                                         :: WindDepth = 0.0D0

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
logical                                          :: noise_mult    = .False.
logical                                          :: noise_model_0 = .False.
logical                                          :: noise_model_1 = .False.
real(dp)                                         :: noise_mu  = 0.0D0
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
character(len=2000)                              :: trajectory_name = 'alm.nc'  ! Traj
character(len=maxlen)                            :: trajectory_final = 'end.dat'  ! Traj final
character(len=maxlen)                            :: output_units=''
character(len=maxlen)                            :: output_calendar=''
character(len=maxlen), dimension(:), allocatable :: Table_filename
integer, dimension(:), allocatable               :: Table_iu
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
integer                                          :: output_tempid
integer                                          :: output_psalid
integer                                          :: output_rid
integer                                          :: output_xoid
integer                                          :: output_yoid
integer                                          :: output_zoid
integer                                          :: output_toid
integer                                          :: output_status
integer                                          :: output_kid


! ... Analytical density
! ... water_visc: Kinematic viscosity [m2/s] = Dynamic viscosity / rho
! ... model_value_rho: Water density kg/m3
! ...
logical                                          :: Particle_buoyant = .False.
integer                                          :: water_density_method=-1
real(dp)                                         :: water_visc = 1D-6 ! m2/s
real(dp)                                         :: model_value_rho  = 1026.0D0 !kg/m3

! ... Linear Equation of State
! ... Simplified EOS (S-EOS) from NEMO ocean model
! ...
logical                                          :: model_density = .False.
logical                                          :: WithEOS  = .False.
real(dp)                                         :: EOS_rho0 = 1026.0D0   ! kg/m3
real(dp)                                         :: EOS_a0   = 1.6550D-1  ! therm. expan.
real(dp)                                         :: EOS_b0   = 7.6554D-1  ! halin. expan.
real(dp)                                         :: EOS_lam1 = 5.9520D-2  ! T2 cabbeling
real(dp)                                         :: EOS_lam2 = 5.4914D-4  ! S2 cabbeling
real(dp)                                         :: EOS_nu   = 2.4341D-3  ! TS cabbeling
real(dp)                                         :: EOS_mu1  = 1.4970D-4  ! T thermobaric
real(dp)                                         :: EOS_mu2  = 1.1090D-5  ! S thermobaric

! ... Diel Vertical Migration
! ...
logical                                          :: Particle_dvm = .False.
real(dp)                                         :: dvm_zday     = 0.0D0
real(dp)                                         :: dvm_znight   = 0.0D0
real(dp)                                         :: dvm_tvm      =  6.0D0*3600.0D0 ! Six hours
real(dp)                                         :: dvm_tdawn    =  7.0D0*3600.0D0 ! 7 am
real(dp)                                         :: dvm_tdusk    = 19.0D0*3600.0D0 ! 7 pm
real(dp)                                         :: dvm_w        = 0.0D0           ! Speed

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
    real(dp) f1,f2,ff,a1,a2,wx,wy,tt,uo,vo,wo
    real(dp) t1,t2
    real(dp) veps,vx,vy,vz,tk,sk,rhok,wok

    if (verb.ge.4) write(*,*) 'RHS time, xo  : ', t, xo

    if (noise_model_1) then
      ! ... Horizontal noise
      veps = randn()
      vx = noise_V1*veps*sin(dpi*veps)
      vy = noise_V1*veps*cos(dpi*veps)
      veps = randn()
      vz = noise_W1*veps
    else
      vx = 0.0D0
      vy = 0.0D0
      vz = 0.0D0
    endif
    if (verb.ge.5) write(*,*) 'RHS Noise model 1 : ', vx, vy, vz

    uo = 0.0D0; vo = 0.0D0; wo = 0.0D0
    dxdt(:) = 0.0D0

    ! ... Zonal advection:
    ! ...
    if (Uadv) then
      if (model_fixed_ou) then
        uo = alpha(1)*model_value_ou
      else if (alpha(1).eq.0) then
        uo = 0.0D0
      else
        uo = time_interpol(GOU,OUrhs,t,xo)
        uo = alpha(1)*uo
      endif
      dxdt(1) = uo + vx
    endif

    ! ... Meridional advection
    ! ...
    if (Vadv) then
      if (model_fixed_ov) then
        vo = alpha(2)*model_value_ov
      else if (alpha(2).eq.0) then
        vo = 0.0D0
      else
        vo = time_interpol(GOV,OVrhs,t,xo)
        vo = alpha(2)*vo
      endif
      dxdt(2) = vo + vy
    endif

    ! ... Vertical advection
    ! ...
    if (Wadv) then
      if (model_fixed_ow) then
        wo = alpha(3)*model_value_ow 
      else if (alpha(3).eq.0) then
        wo = 0.0D0
      else
        wo = time_interpol(GOW,OWrhs,t,xo)
        wo = alpha(3)*wo
      endif
      dxdt(3) = wo + vz
    endif

    if (verb.ge.4) write(*,*) 'RHS u : ', uo, vo, wo

    ! ... Add, if requested, a buoyancy term:
    ! ...
    if (Particle_buoyant) then
      select case (water_density_method)
      case (0)
        ! ... Density is constant
        ! ... 
        rhok = model_value_rho
      case (1)
        ! ... Density is a function of z and latitude
        ! ... 
        rhok = analytical_rho(xo(3),rad2deg*xo(2))
      case (2)
        ! ... Density comes from a RHO file
        ! ... 
        IOR = GOR%locate(xo)
        if (GOR%Stationary) then
          ff  = GOR%interpol(ORrhs(:,:,:,1),xo,IOR,SingleLayer)
        else
          t1 = GOR%t(GOR%rec1)
          t2 = GOR%t(GOR%rec2)
          f1 = GOR%interpol(ORrhs(:,:,:,1),xo,IOR,SingleLayer)
          f2 = GOR%interpol(ORrhs(:,:,:,2),xo,IOR,SingleLayer)
          rhok = f1 + (f2-f1)*(t-t1)/(t2-t1)
        endif
      case (3)
        ! ... Density comes from a TEMP file
        ! ... 
        tk = time_interpol(GOT,OTrhs,t,xo)
        sk = 0.D0
        rhok = eos(tk,sk,abs(xo(3)))
        if (verb.ge.4) write(*,*) 'RHS Temperature, density: ', tk, rhok
      case (4)
        ! ... Density comes from a TEMP and PSAL files
        ! ... 
        tk = time_interpol(GOT,OTrhs,t,xo)
        sk = time_interpol(GOS,OSrhs,t,xo)
        rhok = eos(tk,sk,abs(xo(3)))
        if (verb.ge.4) write(*,*) 'RHS ot, os, or: ', tk, sk, rhok
      end select

      ! ... Emergence velocity
      ! ...
      wok = emergence_model(FLTk,rhok)
      if (verb.ge.4) write(*,*) 'RHS pr, or, wok: ', FLT%Pdens, rhok, wok
      dxdt(3) = dxdt(3) + wok
    endif

    if (Particle_dvm) then
      wok =  dvm_model(xo,t)
      if (verb.ge.4) write(*,*) 'RHS particle motility: ', wok
      dxdt(3) = dxdt(3) + wok
    endif

    ! ... Z is 0 at the surface and negative downwards
    ! ...
    if (winds.and.xo(3).ge.WindDepth) then
      if (model_fixed_au) then
        a1 = model_value_au
      else
        a1 = time_interpol(GAU,AUrhs,t,xo)
      endif
      if (model_fixed_av) then
        a2 = model_value_av
      else
        a2 = time_interpol(GAV,AVrhs,t,xo)
      endif
      wx = A11*a1 + A12*a2
      wy = A21*a1 + A22*a2
      if (verb.ge.4) then
        write(*,*) 'RHS Wind components,   a1, a2 : ', a1, a2
        write(*,*) 'RHS Wind contribution, Wx, Wy : ', wx, wy
      endif
      dxdt(1) = dxdt(1) +  wx
      dxdt(2) = dxdt(2) +  wy
    endif
    if (verb.ge.5) write(*,*) 'RHS dxdt : ', dxdt

  end function RHS
  ! ...
  ! ===================================================================
  ! ...
  subroutine model_ini()

    type(type_date)                         :: datemin,datemax
    integer i,j,k
    real(dp) critical_size,yave,rhok,dx,vmax
    real(dp) ww(ndims)

    if (reverse) model_sign = -1.0D0

    Uadv = model_fixed_ou.or.WithOU
    Vadv = model_fixed_ov.or.WithOV
    Wadv = model_fixed_ow.or.WithOW
    Udrv = model_fixed_au.or.WithAU
    Vdrv = model_fixed_av.or.WithAV
    if (Udrv.neqv.Vdrv) call crash('Incompatible wind forcing options')

    ! ... Model time integration parameters 
    ! ...
    model_tini = alm_tini
    model_tfin = alm_tfin
    model_tlen = nint(model_tfin - model_tini)
    model_dt   = model_sign*model_dt
    if (mod(model_tlen,model_dt).ne.0) call crash('Simulation length not a multiple of time step.')
    model_nsteps = nint(model_tlen/model_dt)

    model_dini = num2date(model_tini,units=model_time_units,calendar=model_time_calendar)
    model_dfin = num2date(model_tfin,units=model_time_units,calendar=model_time_calendar)

    ! ... Propagate the info about the forcing field being a climatology
    ! ...
    GOU%Climatology = OUClim
    GOV%Climatology = OVClim
    GOW%Climatology = OWClim
    GAU%Climatology = AUClim
    GAV%Climatology = AVClim

    !if (model_fixed_ou) GOU%Stationay = .True.
    !if (model_fixed_ov) GOV%Stationay = .True.
    !if (model_fixed_ow) GOW%Stationay = .True.
    !if (model_fixed_au) GAU%Stationay = .True.
    !if (model_fixed_av) GAV%Stationay = .True.

    if (WithOU.and..not.model_fixed_ou) then
      allocate (OUrhs(GOU%nx,GOU%ny,GOU%nz,2))
      call get_forcing('OU',GOU,model_tini,model_dini,0.0D0,OUrhs)
    endif

    if (WithOV.and..not.model_fixed_ov) then
      allocate (OVrhs(GOV%nx,GOV%ny,GOV%nz,2))
      call get_forcing('OV',GOV,model_tini,model_dini,0.0D0,OVrhs)
    endif

    if (WithOW.and..not.model_fixed_ow) then
      allocate (OWrhs(GOW%nx,GOW%ny,GOW%nz,2))
      call get_forcing('OW',GOW,model_tini,model_dini,0.0D0,OWrhs)
    endif

    if (WithOT) then
      allocate (OTrhs(GOT%nx,GOT%ny,GOT%nz,2))
      call get_forcing('OT',GOT,model_tini,model_dini,0.0D0,OTrhs)
    endif

    if (WithOS) then
      allocate (OSrhs(GOS%nx,GOS%ny,GOS%nz,2))
      call get_forcing('OS',GOS,model_tini,model_dini,0.0D0,OSrhs)
    endif

    if (WithOR) then
      stop 'recode'
      allocate (ORrhs(GOR%nx,GOR%ny,GOR%nz,2))
      GOR%rec1 = locate(GOR%t,model_tini) 
      if ((GOR%rec1.lt.1.or.GOR%rec1.ge.GOR%nt)) call crash ('R time not bracketed')
      if (GOR%nt.eq.1.or.GOR%Stationary) then
        GOR%rec2 = GOR%rec1
      else
        GOR%rec2 = min(GOR%rec1+1,GOR%nt)
      endif
      ORrhs(:,:,:,1) = GOR%read3D(GOR%varid,step=GOR%rec1,verbose=verb.ge.3)
      ORrhs(:,:,:,2) = GOR%read3D(GOR%varid,step=GOR%rec2,verbose=verb.ge.3)
    endif

    if (WithOC) then
      stop 'recode'
      allocate (OCrhs(GOC%nx,GOC%ny,GOC%nz,2))
      GOC%rec1 = locate(GOC%t,model_tini) 
      if ((GOC%rec1.lt.1.or.GOC%rec1.ge.GOC%nt)) call crash ('C time not bracketed')
      if (GOC%nt.eq.1.or.GOC%Stationary) then
        GOC%rec2 = GOC%rec1
      else
        GOC%rec2 = min(GOC%rec1+1,GOC%nt)
      endif
      OCrhs(:,:,:,1) = GOC%read3D(GOC%varid,step=GOC%rec1,verbose=verb.ge.3)
      OCrhs(:,:,:,2) = GOC%read3D(GOC%varid,step=GOC%rec2,verbose=verb.ge.3)
    endif

    if (WithAU.and..not.model_fixed_au) then
      allocate (AUrhs(GAU%nx,GAU%ny,1,2))
      call get_forcing('AU',GAU,model_tini,model_dini,0.0D0,AUrhs)
    endif

    if (WithAV.and..not.model_fixed_av) then
      allocate (AVrhs(GAV%nx,GAV%ny,1,2))
      call get_forcing('AV',GAV,model_tini,model_dini,0.0D0,AVrhs)
    endif

    ! ... Stability criteria
    ! ...
    if (GOU%Cartesian) then
      dx = (GOU%xmax-GOU%xmin)/(GOU%nx-1)          ! meters
    else
      dx = Rearth*(GOU%xmax-GOU%xmin)/(GOU%nx-1)   ! meters
    endif
    if (model_fixed_ou) then
      vmax = model_value_ou
    else
      vmax = maxval(OUrhs(:,:,1,1))                ! meters / seconds
    endif
    if (dx.lt.vmax*model_dt) call crash('ERROR: U*dt > dx  -  Reduce dt')

    if (GOV%Cartesian) then
      dx = (GOV%xmax-GOV%xmin)/(GOV%nx-1)          ! meters
    else
      dx = Rearth*(GOV%xmax-GOV%xmin)/(GOV%nx-1)   ! meters
    endif
    if (model_fixed_ov) then
      vmax = model_value_ov
    else
      vmax = maxval(OVrhs(:,:,1,1))                ! meters / seconds
    endif
    if (dx.lt.vmax*model_dt) call crash('ERROR: V*dt > dx  -  Reduce dt')


    if (SingleLayer) Wadv = .False.

    save_frequency = nint(save_period/model_dt)

    noise_KH0 = abs(noise_KH0)
    noise_KV0 = abs(noise_KV0)
    noise_KH1 = abs(noise_KH1)
    noise_KV1 = abs(noise_KV1)

    noise_V0 = sqrt(2.0D0*noise_KH0/abs(model_dt))   ! 
    noise_V1 = sqrt(2.0D0*noise_KH1/abs(model_dt))

    noise_W0 = sqrt(2.0D0*noise_KV0/abs(model_dt))   ! 
    noise_W1 = sqrt(2.0D0*noise_KV1/abs(model_dt))

    if (noise_V0+noise_W0.gt.0.0D0) noise_model_0 = .True.
    if (noise_V1+noise_W1.gt.0.0D0) noise_model_1 = .True.
   
    ! ... Get the "typical" water density.
    ! ... if water_density_method .le. 0, then it is the
    ! ... value in model_value_rho.
    ! ..
    if (model_density) then
      if (water_density_method.eq.1) then
        yave = 0.5D0*(alm_ymin+alm_ymax)*rad2deg
        model_value_rho = 0.0D0
        if (verb.ge.3) then
          write(*,*)
          write(*,*) '   Z         rho(z) '
          write(*,*) '--------------------'
        endif
        do k=1,GOU%Nz
          rhok = analytical_rho(GOU%z(k),yave)
          model_value_rho = model_value_rho + rhok
          if (verb.ge.3) write(*,*) GOU%z(k), rhok
        enddo
        model_value_rho = model_value_rho / GOU%Nz
      else if (water_density_method.eq.2) then
        print*, 'WARNING: This needs to be calculated !!!!'
        yave = 0.5D0*(alm_ymin+alm_ymax)*rad2deg
        model_value_rho = 0.0D0
        do k=1,GOU%Nz
          model_value_rho = model_value_rho + analytical_rho(GOU%z(k),yave)
        enddo
        model_value_rho = model_value_rho / GOU%Nz
      endif
    endif

    ! ... Diel Vertical Movement
    ! ...
    if (Particle_dvm) then
      dvm_zday   = -abs(dvm_zday)
      dvm_znight = -abs(dvm_znight)
      dvm_w      = (dvm_znight - dvm_zday)/dvm_tvm
    endif

    ! ... Wind driven model  U = beta*exp(alpha*j) * W
    ! ... ou = beta* (cos(alpha)*au - sin(alpha)*av)
    ! ... ov = beta* (sin(alpha)*au - cos(alpha)*av)
    ! ...
    if (WindDriven) then
      A11 = WDriven_beta*cos(deg2rad*WDriven_theta); A22 = A11
      A21 = WDriven_beta*sin(deg2rad*WDriven_theta); A12 = -A21
    endif


    if (verb.ge.1) then
      write(*,*)
      write(*,*) 'Verbosity level  : ', verb
      write(*,*) 
      write(*,*) 'Simulation period: '
      if (model_sign.gt.0) then
        write(*,*) 'Forward model'
      else
        write(*,*) 'Backward model'
      endif
      write(*,*) 'Runge Kutta   : ', rk_order
      write(*,*) 'Initial time  : ', model_tini, model_dini%iso()
      write(*,*) 'Final time    : ', model_tfin, model_dfin%iso()
      write(*,*) 'Time step     : ', model_dt
      write(*,*) 'Number steps  : ', model_nsteps
      write(*,*) 'Saving period : ', save_period
      write(*,*) 'Advection XYZ : ', Uadv, Vadv, Wadv
      if (OUClim) then
        write(*,*) 'Advection X   : Climatological'
      else if (model_fixed_ou) then
        write(*,*) 'Advection X   : Constant'
      else if (Uadv) then
        write(*,*) 'Advection X   : Interannual'
      endif
      if (OVClim) then
        write(*,*) 'Advection Y   : Climatological'
      else if (model_fixed_ov) then
        write(*,*) 'Advection Y   : Constant'
      else if (Vadv) then
        write(*,*) 'Advection Y   : Interannual'
      endif
      if (OWClim) then
        write(*,*) 'Advection Z   : Climatological'
      else if (model_fixed_ow) then
        write(*,*) 'Advection Z   : Constant'
      else if (Wadv) then
        write(*,*) 'Advection Z   : Interannual'
      endif

      if (model_density) then
        write(*,*) 'Tracking water density activated'
        write(*,*) 'Using Equation of State   : ', WithEOS
        write(*,*) 'Water density method      : ', water_density_method
      else
        write(*,*) 'Tracking water density not activated'
      endif

      if (Particle_buoyant) then
        write(*,*) 'Particle buoyant activated'
        write(*,*) 'Water kinematic viscosity : ', water_visc
        if (water_density_method.eq.0.or.water_density_method.eq.1) &
          write(*,*) 'Reference Water density   : ', model_value_rho
        write(*,*) 'Particle density          : ', Release_rho
        write(*,*) 'Particle diameter         : ', Release_size
        if (Release_rho.lt.model_value_rho) then
          ! ... Maximum drop sizes are claculated by Aravamudan et al. (1982)
          ! ... Break up for oil on rough seas - simplified models and step-by
          ! ... step calculations. US Coast Guard Report CG-D-28-82
          ! ... US Department of Transportation, Washington DC.
          ! ...
          critical_size = 9.52D0*(water_visc/gravity)**(2.D0/3.D0)/(1.0D0-Release_rho/model_value_rho)**(1.0D0/3.0D0)
          write(*,*) 'Maximum Particle diameter : ', critical_size
          !if (Release_size.gt.critical_size) call crash('Particle size > Minumum allowed size')
        endif
      else
        write(*,*) 'Particle not buoyant'
      endif
      if (Particle_dvm) then
        write(*,*) 'Particle DVM activated'
        write(*,*) 'Particle Day depth        : ', dvm_zday
        write(*,*) 'Particle night depth      : ', dvm_znight
        write(*,*) 'Particle travelling time  : ', dvm_tvm
        write(*,*) 'Particle vertical speed   : ', dvm_w
      else
        write(*,*) 'Particle DVM not activated'
      endif

      write(*,*) 'Wind forcing           : ', Winds

      write(*,*) 'Wind response method   : ', WindResponse
      write(*,*) 'Wind driven method     : ', WindDriven
      if (WindDriven) then
        write(*,*) 'Wind driven beta       : ', WDriven_beta
        write(*,*) 'Wind driven theta      : ', WDriven_theta
      endif
      write(*,*) 'A11, A12               : ', A11, A12
      write(*,*) 'A21, A22               : ', A21, A22
      write(*,*) 'Wind action depth      : ', WindDepth

      write(*,*) 'Muliplicative noise    : ', noise_mu
      write(*,*) 'Additive noise model 0 : ', noise_model_0
      write(*,*) 'Additive noise model 1 : ', noise_model_1
      if (Uadv.and..not.model_fixed_ou) &
         write(*,*) 'Zonal velocity record pointers     : ', GOU%rec1, GOU%rec2
      if (Vadv.and..not.model_fixed_ov) &
         write(*,*) 'Meridional velocity record pointers: ', GOV%rec1, GOV%rec2
      if (Wadv.and..not.model_fixed_ow) &
         write(*,*) 'Vertical velocity record pointers  : ', GOW%rec1, GOW%rec2
      if (WithOT) &
         write(*,*) 'Ocean temperature record pointers  : ', GOT%rec1, GOT%rec2
      if (WithOS) &
         write(*,*) 'Ocean salinity record pointers     : ', GOS%rec1, GOS%rec2
      if (Winds) then
        if (.not.model_fixed_au) &
           write(*,*) 'Zonal wind record pointers         : ', GAU%rec1, GAU%rec2
        if (.not.model_fixed_av) &
           write(*,*) 'Meridional wind record pointers    : ', GAV%rec1, GAV%rec2
      endif
    endif

    !stop '8888'
    return

  end subroutine model_ini
  ! ...
  ! ===================================================================
  ! ...
  subroutine model_run()

    ! ... Local variables
    ! ...
    logical OUupdate,OVupdate,OWupdate,AUupdate,AVupdate
    logical OTupdate,OSupdate
    integer i,j,k,ii,jj,kk,ll,ifloat,float_status,iu
    real(dp) veps,dx,dy,dz
    real(dp), dimension(ndims)                   :: xo,xn,un

    
    ! ... Flags to keep updating forcing fields
    ! ...
    OUupdate = WithOU.and..not.GOU%Stationary.and..not.model_fixed_ou
    OVupdate = WithOV.and..not.GOV%Stationary.and..not.model_fixed_ov
    OWupdate = WithOW.and..not.GOW%Stationary.and..not.model_fixed_ow
    AUupdate = WithAU.and..not.GAU%Stationary.and..not.model_fixed_au
    AVupdate = WithAV.and..not.GAV%Stationary.and..not.model_fixed_av
    OTupdate = WithOT.and..not.GOT%Stationary
    OSupdate = WithOS.and..not.GOS%Stationary
    if (verb.ge.3) then
      write(*,*) 'OU update: ', OUupdate
      write(*,*) 'OV update: ', OVupdate
      write(*,*) 'OW update: ', OWupdate
      if (WithOT) write(*,*) 'OT update: ', OTupdate
      if (WithOS)  write(*,*) 'OS update: ', OSupdate
      write(*,*) 'AU update: ', AUupdate
      write(*,*) 'AV update: ', AVupdate
    endif

    ! ... First record:
    ! ...
    model_time = model_tini 
    do ifloat=1,Nfloats
      FLT(ifloat)%to = model_tini + FLT(ifloat)%to    ! FLT%to seconds since model_tini
      if (model_sign*model_time.ge.model_sign*FLT(ifloat)%to) then
        if (verb.ge.2) write(*,*) 'Initial step releasing for float ', ifloat
        FLT(ifloat)%released = .True.
        FLT(ifloat)%floating = .True.
        FLT(ifloat)%indomain = .True.
        FLT(ifloat)%x        = FLT(ifloat)%xo
        FLT(ifloat)%y        = FLT(ifloat)%yo
        FLT(ifloat)%z        = FLT(ifloat)%zo
        FLT(ifloat)%dist     = 0.0D0
        FLT(ifloat)%tlast    = model_tini
        xn = [FLT(ifloat)%x,FLT(ifloat)%y,FLT(ifloat)%z]
        if (WithOT) then
          FLT(ifloat)%Wtemp  = time_interpol(GOT,OTrhs,model_time,xn)
        endif
        if (WithOS) then
          FLT(ifloat)%Wpsal  = time_interpol(GOS,OSrhs,model_time,xn)
        endif
        if (model_density) then
          select case (water_density_method)
          case (0)
            ! ... Density is constant
            ! ... 
            FLT(ifloat)%Wdens = model_value_rho
          case (1)
            ! ... Density is a function of z and latitude
            ! ... 
            FLT(ifloat)%Wdens = analytical_rho(abs(xn(3)),rad2deg*xn(2))
          case (2)
            ! ... Density comes from a RHO file
            ! ... 
            FLT(ifloat)%Wdens  = time_interpol(GOR,ORrhs,model_time,xn)
          case (3)
            ! ... Density comes from TEMP
            ! ... 
            FLT(ifloat)%Wdens = eos(FLT(ifloat)%Wtemp,35.0D0,abs(xn(3)))
          case (4)
            ! ... Density comes from TEMP and PSAL
            ! ... 
            FLT(ifloat)%Wdens = eos(FLT(ifloat)%Wtemp,FLT(ifloat)%Wpsal,abs(xn(3)))
          end select
        endif
      endif
    enddo
    call trajectory_write(model_time,verb.ge.2) 


    if (verb.ge.1) write(*,*) 'Running model ...'

    ! --------------------------  steps
    do model_step=1,model_nsteps
    ! --------------------------

      model_time = anint(model_tini + (model_step-1)*model_dt)
      model_date = num2date(Reference_time+model_time,model_time_units,model_time_calendar)
      model_daysec = anint(model_date%hour*3600.0D0 + model_date%minute*60.0D0 + model_date%second)

      if (verb.ge.2) write(*,*) 'step, time, date, secs :: ', model_step, model_time, &
                                                              model_date%iso(), model_daysec

      if (OUupdate) call get_forcing('OU',GOU,model_time,model_date,0.0D0,OUrhs)
      if (OVupdate) call get_forcing('OV',GOV,model_time,model_date,0.0D0,OVrhs)
      if (OWupdate) call get_forcing('OW',GOW,model_time,model_date,0.0D0,OWrhs)
      if (AUupdate) call get_forcing('AU',GAU,model_time,model_date,0.0D0,AUrhs)
      if (AVupdate) call get_forcing('AV',GAV,model_time,model_date,0.0D0,AVrhs)
      if (OTupdate) call get_forcing('OT',GOT,model_time,model_date,0.0D0,OTrhs)
      if (OSupdate) call get_forcing('OS',GOS,model_time,model_date,0.0D0,OSrhs)

      do ifloat=1,Nfloats
      ! ------------------ floats

        ! ... First chek if released or going to be released
        ! ...
        if (.NOT.FLT(ifloat)%released) then
          if (model_sign*model_time.ge.model_sign*FLT(ifloat)%to) then
            if (verb.ge.2) write(*,*) 'Releasing float ', ifloat
            FLT(ifloat)%released = .True.
            FLT(ifloat)%floating = .True.
            FLT(ifloat)%indomain = .True.
            FLT(ifloat)%x        = FLT(ifloat)%xo
            FLT(ifloat)%y        = FLT(ifloat)%yo
            FLT(ifloat)%z        = FLT(ifloat)%zo
            FLT(ifloat)%dist     = 0.0D0
            FLT(ifloat)%tlast    = model_time
            xn = [FLT(ifloat)%x,FLT(ifloat)%y,FLT(ifloat)%z]
            if (WithOT) then
              FLT(ifloat)%Wtemp     = time_interpol(GOT,OTrhs,model_time,xn) 
            endif
            if (WithOS) then
              FLT(ifloat)%Wpsal     = time_interpol(GOS,OSrhs,model_time,xn)
            endif
            if (model_density) then
              select case (water_density_method)
              case (0)
                ! ... Density is constant
                ! ... 
                FLT(ifloat)%Wdens = model_value_rho
              case (1)
                ! ... Density is a function of z and latitude
                ! ... 
                FLT(ifloat)%Wdens = analytical_rho(abs(xn(3)),rad2deg*xn(2))
              case (2)
                ! ... Density comes from a RHO file
                ! ... 
                FLT(ifloat)%Wdens  = time_interpol(GOR,ORrhs,model_time,xn)
              case (3)
                ! ... Density comes from TEMP
                ! ... 
                FLT(ifloat)%Wdens = eos(FLT(ifloat)%Wtemp,35.0D0,abs(xn(3)))
              case (4)
                ! ... Density comes from TEMP and PSAL
                ! ... 
                FLT(ifloat)%Wdens = eos(FLT(ifloat)%Wtemp,FLT(ifloat)%Wpsal,abs(xn(3)))
              end select
            endif
          endif
        endif

        if (FLT(ifloat)%floating) then

          FLTk = FLT(ifloat)
          xo = [FLT(ifloat)%x, FLT(ifloat)%y, FLT(ifloat)%z ]

          ! ... Check the dawn and dusk time. It should be expected that these
          ! ... values remain constant during a single time step.
          ! ...
          if (Particle_dvm) then
            call SunRise_and_SunSet (rad2deg*xo(1),rad2deg*xo(2),model_date,dvm_tdawn,dvm_tdusk)
            if (verb.ge.4) write(*,*) 'DVM: tdawn, tdusk (hours): ', dvm_tdawn/3600.0D0, dvm_tdusk/3600.0D0
          endif


          ! ... Runge - Kutta
          ! ....................................................
          call rk(xo,model_time,model_dt,un,xn)
          ! ....................................................

          ! ... Check float status:
          ! ...
          float_status = point_type(xn(1),xn(2),xn(3))  ! = 1, floating inside
                                                        ! = 0, stranded on land
                                                        ! = -1, outside the domain
          if (verb.ge.3) then
            write(*,*) 'Float: ', ifloat
            write(*,*) '   xo: ', xo
            write(*,*) '   un: ', un
            write(*,*) '   xn: ', xn
            write(*,*) ' stat: ', float_status
          endif

          select case (float_status)
            case (1)
            ! ................... The particle is floating inside the domain
            ! ................... Values are updated
              FLT(ifloat)%x = xn(1)
              FLT(ifloat)%y = xn(2)
              FLT(ifloat)%z = xn(3)
              FLT(ifloat)%u = un(1)
              FLT(ifloat)%v = un(2)
              FLT(ifloat)%w = un(3)
              if (Spherical) then
                dx = 0.001D0*Rearth*cos(xo(2))*(xn(1)-xo(1))
                dy = 0.001D0*Rearth*(xn(2)-xo(2))
              else
                dx = xn(1)-xo(1)
                dy = xn(2)-xo(2)
              endif
              FLT(ifloat)%dist  = FLT(ifloat)%dist + sqrt(dx*dx+dy*dy)
              FLT(ifloat)%tlast = model_time + model_dt
            ! ...................
            case (0) 
              if (verb.ge.2) write(*,*) 'Particle stranded on land: ', ifloat
              FLT(ifloat)%floating = .False.
            case (-1)
              if (verb.ge.2) write(*,*) 'Particle leaving the domain: ', ifloat
              FLT(ifloat)%indomain = .False.
              FLT(ifloat)%floating = .False.
          end select

        endif

      ! ------------------ floats
      enddo

      model_time = model_time + model_dt

      if (mod(model_step,save_frequency).eq.0) then
        ! ...
        ! ... If requested, calculate the temp, psal and density of the 
        ! ... water at the particle position 
        ! ...
        do ifloat=1,Nfloats
          xn = [FLT(ifloat)%x, FLT(ifloat)%y, FLT(ifloat)%z ]
          if (WithOT) then
            FLT(ifloat)%Wtemp  = time_interpol(GOT,OTrhs,model_time,xn)
          endif
          if (WithOS) then
            FLT(ifloat)%Wpsal  = time_interpol(GOS,OSrhs,model_time,xn)
          endif
          if (model_density) then
            select case (water_density_method)
            case (0)
              ! ... Density is constant
              ! ... 
              FLT(ifloat)%Wdens = model_value_rho
            case (1)
              ! ... Density is a function of z and latitude
              ! ... 
              FLT(ifloat)%Wdens = analytical_rho(abs(xn(3)),rad2deg*xn(2))
            case (2)
              ! ... Density comes from a RHO file
              ! ... 
              FLT(ifloat)%Wdens  = time_interpol(GOR,ORrhs,model_time,xn)
            case (3)
              ! ... Density comes from TEMP
              ! ... 
              FLT(ifloat)%Wdens = eos(FLT(ifloat)%Wtemp,35.0D0,abs(xn(3)))
            case (4)
              ! ... Density comes from TEMP and PSAL
              ! ... 
              FLT(ifloat)%Wdens = eos(FLT(ifloat)%Wtemp,FLT(ifloat)%Wpsal,abs(xn(3)))
            end select
          endif
        enddo
        call trajectory_write(model_time,verb.ge.2)
      endif


    ! --------------------------  steps
    enddo                            
    ! --------------------------

    ! ... Saving the las position
    ! ...
    if (len_trim(trajectory_final).gt.0) then
      iu = unitfree()
      open(iu,file=trajectory_final,status='unknown')
      do ifloat=1,Nfloats
        model_time = FLT(ifloat)%tlast
        model_date = num2date(Reference_time+model_time,model_time_units,model_time_calendar)
        write(iu,'(F11.6,2X,F11.6,2X,F8.3,2X,A)') FLT(ifloat)%x*rad2deg, &
                                                FLT(ifloat)%y*rad2deg, &
                                                abs(FLT(ifloat)%z), &
                                                model_date%iso()
      enddo
      close(iu)
    endif

    if (verb.ge.1) write(*,*) 'Done !'

  end subroutine model_run
  ! ...
  ! ===================================================================
  ! ...
  subroutine get_forcing(label,G,time,date,missing,f)

    character(len=2), intent(in)                         :: label
    type(type_ncgrid), intent(inout)                     :: G
    real(dp), intent(in)                                 :: time
    type(type_date), intent(in)                          :: date
    real(dp), intent(in)                                 :: missing
    real(dp), dimension(G%nx,G%ny,G%nz,2), intent(inout) :: f

    ! ... Local variables
    ! ...
    integer d1,m1,y1,d2,m2,y2
    real(dp) rd

    ! ... Check if initial call
    ! ...
    if (G%rec1.lt.0) then
      if (verb.ge.5) write(*,*) 'get_forcing initial call for ', label, ' time: ', time
      ! .....................................................
      ! ...                  Initial call                 ...
      ! .....................................................
      if (G%Climatology) then
        rd = date%day + (date%hour*3600.0D0+date%minute*60.0D0+date%second)/86400.0D0
        y1 = date%year
        m1 = date%month
        if (rd.lt.HalfMonth(date%month)) then
          y2 = y1
          m2 = m1
          m1 = m1 - 1
          if (m1.eq.0) then
            y1 = y1 - 1
            m1 = 12
          endif
        else
          y2 = y1
          m2 = m1 + 1
          if (m2.eq.13) then
            y2 = y2 + 1
            m2 = 1
          endif
        endif
        G%rec1 = m1
        call G%drec1%is(y1,m1,1,0,0,0,model_time_calendar)
        G%trec1 = anint(date2num(G%drec1,units=model_time_units) + (HalfMonth(m1)-1)*86400.0D0)
        G%rec2 = m2
        call G%drec2%is(y2,m2,1,0,0,0,model_time_calendar)
        G%trec2 = anint(date2num(G%drec2,units=model_time_units) + (HalfMonth(m2)-1)*86400.0D0)
      else
        if (G%nt.eq.1) then
          G%rec1 = 1
          G%rec2 = 1
        else
          G%rec1 = locate(G%t,time)
          if ((G%rec1.lt.1.or.G%rec1.ge.G%nt)) call crash (label//' time not bracketed')
          if (G%Stationary) then
            G%rec2 = G%rec1
          else
            G%rec2 = min(G%rec1+1,G%nt)
          endif
        endif
        G%trec1 = G%t(G%rec1)
        G%trec2 = G%t(G%rec2)
      endif
      if (verb.ge.5) write(*,*) 'rec1 and rec 2: ', G%rec1, G%rec2
      f(:,:,:,1) = G%read3D(G%varid,step=G%rec1,missing_value=missing,verbose=verb.ge.3)
      f(:,:,:,2) = G%read3D(G%varid,step=G%rec2,missing_value=missing,verbose=verb.ge.3)


    else if (.not.G%Stationary) then
      ! .....................................................
      ! ...                  Update fields                ...
      ! .....................................................
      if (G%Climatology) then
        if (reverse) then                 ! Backward model
          if (time.lt.G%trec1) then
            m1 = G%drec1%month
            y1 = G%drec1%year
            y2 = y1
            m2 = m1
            m1 = m1 - 1
            if (m1.eq.0) then
              y1 = y1 - 1
              m1 = 12
            endif
            G%rec2  = m2
            G%trec2 = G%trec1
            G%drec2 = G%drec1
            G%rec1  = m1
            call G%drec1%is(y1,m1,1,0,0,0,model_time_calendar)
            G%trec1 = anint(date2num(G%drec1,units=model_time_units)+(HalfMonth(m1)-1)*86400.0D0)
            f(:,:,:,2) = f(:,:,:,1)
            f(:,:,:,1) = G%read3D(G%varid,step=G%rec1,missing_value=missing,verbose=verb.ge.3)
            if (verb.ge.2) write(*,'(T2,"Backward bracket ",A2," update:",I2," - ",I2)') &
                           label, G%rec1, G%rec2
          endif
        else                                   ! Forward model
          if (time.gt.G%trec2) then
            m2 = G%drec2%month
            y2 = G%drec2%year
            y1 = y2
            m1 = m2
            m2 = m2 + 1
            if (m2.eq.13) then
              y2 = y2 + 1
              m2 = 1
            endif
            G%rec1  = m1
            G%trec1 = G%trec2
            G%drec1 = G%drec2
            G%rec2  = m2
            call G%drec2%is(y2,m2,1,0,0,0,model_time_calendar)
            G%trec2 = anint(date2num(G%drec2,units=model_time_units)+(HalfMonth(m2)-1)*86400.0D0)
            f(:,:,:,1) = f(:,:,:,2)
            f(:,:,:,2) = G%read3D(G%varid,step=G%rec2,missing_value=missing,verbose=verb.ge.3)
            if (verb.ge.1) write(*,'(T2,"Forward bracket ",A2," update:",I2," - ",I2)') &
                           label, G%rec1, G%rec2
          endif
        endif
      else
        if (reverse) then
          if (time.le.G%t(G%rec1)) then
            G%rec2 = G%rec1
            G%rec1 = G%rec1 - 1
            if (verb.ge.2) write(*,'(T2,"Backward bracket ",A2," update:",F12.0," - ",F12.0)') &
                           label, G%t(G%rec1), G%t(G%rec2)
            f(:,:,:,2) = f(:,:,:,1)
            f(:,:,:,1) = G%read3D(G%varid,step=G%rec1,missing_value=missing,verbose=verb.ge.3)
          endif
        else
          if (time.ge.G%t(G%rec2)) then
            G%rec1 = G%rec2
            G%rec2 = G%rec2 + 1
            if (verb.ge.2) write(*,'(T2,"Forward bracket ",A2," update:",F12.0," - ",F12.0)') &
                           label, G%t(G%rec1), G%t(G%rec2)
            f(:,:,:,1) = f(:,:,:,2)
            f(:,:,:,2) = G%read3D(G%varid,step=G%rec2,missing_value=missing,verbose=verb.ge.3)
          endif
        endif
        G%trec1 = G%t(G%rec1)
        G%trec2 = G%t(G%rec2)
      endif
    endif

    return 
        
  end subroutine get_forcing
  ! ...
  ! ===================================================================
  ! ...
  subroutine forcing_update(label,G,time,f)

    character(len=2), intent(in)                         :: label
    type(type_ncgrid), intent(inout)                     :: G
    real(dp), intent(in)                                 :: time
    real(dp), dimension(G%nx,G%ny,G%nz,2), intent(inout) :: f
        
    if (G%rec1.lt.0) then
      ! ... First call !
      ! ...
      G%rec1 = locate(G%t,time) 
      if ((G%rec1.lt.1.or.G%rec1.ge.G%nt)) call crash (label//' time not bracketed')
      if (G%Stationary) then
        G%rec2 = G%rec1
      else
        G%rec2 = G%rec1 + 1
      endif
      f(:,:,:,1) = G%read3D(G%varid,step=G%rec1,missing_value=0.0D0)
      f(:,:,:,2) = G%read3D(G%varid,step=G%rec2,missing_value=0.0D0)
    else if (.not.G%Stationary) then
      ! ... Update fields
      ! ...
      if (reverse) then
        if (time.le.G%t(G%rec1)) then
          G%rec2 = G%rec1
          G%rec1 = G%rec1 - 1
          if (verb.ge.2) write(*,'(T2,"Backward bracket ",A2," update:",F12.0," - ",F12.0)') &
                         label, G%t(G%rec1), G%t(G%rec2)
          f(:,:,:,2) = f(:,:,:,1)
          f(:,:,:,1) = G%read3D(G%varid,step=G%rec1,missing_value=0.0D0)
        endif
      else
        if (time.ge.G%t(G%rec2)) then
          G%rec1 = G%rec2
          G%rec2 = G%rec2 + 1
          if (verb.ge.2) write(*,'(T2,"Forward bracket ",A2," update:",F12.0," - ",F12.0)') &
                         label, G%t(G%rec1), G%t(G%rec2)
          f(:,:,:,1) = f(:,:,:,2)
          f(:,:,:,2) = G%read3D(G%varid,step=G%rec2,missing_value=0.0D0)
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
  subroutine trajectory_create(filename,Nfloats)

    character(len=*), intent(in)               :: filename
    integer, intent(in)                        :: Nfloats

    ! ... Local variables
    ! ...
    integer err,i,natts,iu
    character(len=1000) lcom
    character(len=maxlen) path,base,ext,afile


    trajectory_name = trim(filename)
    output_nfloats = Nfloats
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

    if (verb.ge.3) write(*,*) 'TRAJECTORY_CREATE, Nfloats, output_format: ', Nfloats, output_format

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

      if (WithOT) then
        err = NF90_DEF_VAR(output_id,'water_temp',NF90_REAL,(/output_nid,output_idl/),output_tempid)
        call check()
        err = NF90_PUT_ATT(output_id,output_tempid,'units','degrees_C')
        call check()
        err = NF90_PUT_ATT(output_id,output_tempid,'_FillValue',sngl(missing))
        call check()
      endif

      if (WithOS) then
        err = NF90_DEF_VAR(output_id,'water_psal',NF90_REAL,(/output_nid,output_idl/),output_psalid)
        call check()
        err = NF90_PUT_ATT(output_id,output_psalid,'units','psu')
        call check()
        err = NF90_PUT_ATT(output_id,output_psalid,'_FillValue',sngl(missing))
        call check()
      endif

      if (model_density) then
        err = NF90_DEF_VAR(output_id,'water_dens',NF90_REAL,(/output_nid,output_idl/),output_rid)
        call check()
        err = NF90_PUT_ATT(output_id,output_rid,'units','kg m-3')
        call check()
        err = NF90_PUT_ATT(output_id,output_rid,'_FillValue',sngl(missing))
        call check()
      endif

      err = NF90_DEF_VAR(output_id,'exitcode',NF90_INT,(/output_nid/),output_status)
      call check()
      err = NF90_PUT_ATT(output_id,output_status,'long_name','Float status at the end of the simulation')
      err = NF90_PUT_ATT(output_id,output_status,"is-1",'float has not been released')
      err = NF90_PUT_ATT(output_id,output_status,"is0",'float was moving')
      err = NF90_PUT_ATT(output_id,output_status,"is1",'float left the model area')
      err = NF90_PUT_ATT(output_id,output_status,"is2",'float was stranded')
      call check()


      err = NF90_PUT_ATT(output_id,0,'Version',VERSION)
      err = NF90_PUT_ATT(output_id,0,'Date',now%iso())
      err = NF90_PUT_ATT(output_id,0,'rk_order',rk_order)
      err = NF90_PUT_ATT(output_id,0,'model_dt',model_dt)

      if (Winds) then
        err = NF90_PUT_ATT(output_id,0,'Wind_forcing','True')
        err = NF90_PUT_ATT(output_id,0,'WindDepth',WindDepth)
        err = NF90_PUT_ATT(output_id,0,'A11',A11)
        err = NF90_PUT_ATT(output_id,0,'A12',A12)
        err = NF90_PUT_ATT(output_id,0,'A21',A21)
        err = NF90_PUT_ATT(output_id,0,'A22',A22)
      else
        err = NF90_PUT_ATT(output_id,0,'Wind_forcing','False')
      endif

      if (WithOT) then
        err = NF90_PUT_ATT(output_id,0,'water_temp','True')
      else
        err = NF90_PUT_ATT(output_id,0,'water_temp','False')
      endif
      if (WithOS) then
        err = NF90_PUT_ATT(output_id,0,'water_psal','True')
      else
        err = NF90_PUT_ATT(output_id,0,'water_psal','False')
      endif
 
      if (model_density) then
        err = NF90_PUT_ATT(output_id,0,'water_dens','True')
        if (WithEOS) then
          err = NF90_PUT_ATT(output_id,0,'WithEOS','True')
        else
          err = NF90_PUT_ATT(output_id,0,'WithEOS','False')
        endif
        err = NF90_PUT_ATT(output_id,0,'Water_density_method',water_density_method)
        if (water_density_method.eq.0.or.water_density_method.eq.1) &
          err = NF90_PUT_ATT(output_id,0,'Reference_Water_density',model_value_rho)
        err = NF90_PUT_ATT(output_id,0,'EOS_rho0',EOS_rho0)
        err = NF90_PUT_ATT(output_id,0,'EOS_a0',EOS_a0)
        err = NF90_PUT_ATT(output_id,0,'EOS_b0',EOS_b0)
        err = NF90_PUT_ATT(output_id,0,'EOS_lam1',EOS_lam1)
        err = NF90_PUT_ATT(output_id,0,'EOS_lam2',EOS_lam2)
        err = NF90_PUT_ATT(output_id,0,'EOS_nu',EOS_nu)
        err = NF90_PUT_ATT(output_id,0,'EOS_mu1',EOS_mu1)
        err = NF90_PUT_ATT(output_id,0,'EOS_mu2',EOS_mu2)
      else
        err = NF90_PUT_ATT(output_id,0,'water_dens','False')
      endif

      if (Particle_buoyant) then
        err = NF90_PUT_ATT(output_id,0,'Particle_buoyant','True')
        err = NF90_PUT_ATT(output_id,0,'Water_viscosity',Water_visc)
        err = NF90_PUT_ATT(output_id,0,'Particle_size',Release_size)
        err = NF90_PUT_ATT(output_id,0,'Particle_density',Release_rho)
      else
        err = NF90_PUT_ATT(output_id,0,'Particle_buoyant','False')
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
      call get_commandline(lcom)
      allocate(table_filename(Nfloats))
      allocate(table_iu(Nfloats))
      call filename_split(filename,path,base,ext)
      do i=1,Nfloats
        write(afile,'(T1,A,A,I4.4,A)') trim(path),trim(base)//'.',i,'.'//trim(ext)
        iu = unitfree()
        Table_filename(i) = trim(afile)
        Table_iu(i) = iu
        open(iu,file=afile,status='unknown')
        rewind(iu)
        write(iu,'(T1,A)') '# ALOGES LAGRANGIAN MODEL trajectory'
        write(iu,'(T1,A)') '# -------------------------------------------------&
                        &--------------------------------------------------------'
        write(iu,'(T1,A)') '# Version: ' // trim(VERSION)
        write(iu,'(T1,A)') '# Runtime: ' // now%iso()
        write(iu,'(T1,A)') '# Command line: ' // trim(lcom)
        write(iu,'(T1,A)') '# -------------------------------------------------&
                        &--------------------------------------------------------'
        write(iu,'(T1,A)') '#  lon      lat      depth         date            &
    &u        v         w        temp     psal     dens     dist  S'
      enddo

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
    integer i,err,iu,code
    type(type_Date) date
    real(dp) time_out(1), xx(Nfloats),axx(10)
    character(len=maxlen) fmt

    if (present(verbose)) then
      verb = verbose
    else
      verb = .False.
    endif

    output_record = output_record + 1

    if (verb) write(*,*) 'Saving record ', output_record, ' at time ', time

    if (output_format.eq.0) then

      ! ... Time
      ! ...
      time_out(1) = time
      err = NF90_PUT_VAR(output_id,output_timeid,time_out,[output_record],[1])
      call cdf_error(err,'Error writing trajectory time')


      ! ... Longitude
      ! ...
      xx(:) = missing
      do i=1,Nfloats
        if (FLT(i)%released) xx(i) = FLT(i)%x*rad2deg
      enddo
      err = NF90_PUT_VAR(output_id,output_lonid,xx,[1,output_record],[Nfloats,1])
      call cdf_error(err,'Error writing trajectory longitude')

      ! ... Latitude
      ! ...
      xx(:) = missing
      do i=1,Nfloats
        if (FLT(i)%released) xx(i) = FLT(i)%y*rad2deg
      enddo
      err = NF90_PUT_VAR(output_id,output_latid,xx,[1,output_record],[Nfloats,1])
      call cdf_error(err,'Error writing trajectory latitude')

      ! ... Depth
      ! ...
      xx(:) = missing
      do i=1,Nfloats
        if (FLT(i)%released) xx(i) = abs(FLT(i)%z)
      enddo
      err = NF90_PUT_VAR(output_id,output_zid,xx,[1,output_record],[Nfloats,1])
      call cdf_error(err,'Error writing trajectory depth')

      ! ... Velocity values are specified at the previus step:
      ! ...
      if (output_record.gt.1) then

        ! ... U
        ! ...
        xx(:) = missing
        do i=1,Nfloats
          if (FLT(i)%released) xx(i) = FLT(i)%u
        enddo
        err = NF90_PUT_VAR(output_id,output_uid,xx,[1,output_record-1],[Nfloats,1])
        call cdf_error(err,'Error writing trajectory zonal velocity')

        ! ... V
        ! ...
        xx(:) = missing
        do i=1,Nfloats
          if (FLT(i)%released) xx(i) = FLT(i)%v
        enddo
        err = NF90_PUT_VAR(output_id,output_vid,xx,[1,output_record-1],[Nfloats,1])
        call cdf_error(err,'Error writing trajectory meridional velocity')

        ! ... W
        ! ...
        xx(:) = missing
        do i=1,Nfloats
          if (FLT(i)%released) xx(i) = FLT(i)%w
        enddo
        err = NF90_PUT_VAR(output_id,output_wid,xx,[1,output_record-1],[Nfloats,1])
        call cdf_error(err,'Error writing trajectory vertical velocity')

      endif

      ! ... Distance
      ! ...
      xx(:) = missing
      do i=1,Nfloats
        if (FLT(i)%released) xx(i) = abs(FLT(i)%dist)
      enddo
      err = NF90_PUT_VAR(output_id,output_did,xx,[1,output_record],[Nfloats,1])
      call cdf_error(err,'Error writing trajectory distance')

      ! ... Temperature
      ! ...
      if (WithOT) then
        xx(:) = missing
        do i=1,Nfloats
          if (FLT(i)%released) xx(i) = FLT(i)%Wtemp
        enddo
        err = NF90_PUT_VAR(output_id,output_tempid,xx,[1,output_record],[Nfloats,1])
        call cdf_error(err,'Error writing water temperature')
      endif

      ! ... Salinity
      ! ...
      if (WithOS) then
        xx(:) = missing
        do i=1,Nfloats
          if (FLT(i)%released) xx(i) = FLT(i)%Wpsal
        enddo
        err = NF90_PUT_VAR(output_id,output_psalid,xx,[1,output_record],[Nfloats,1])
        call cdf_error(err,'Error writing water salinity')
      endif

      ! ... Water density
      ! ...
      if (model_density) then
        xx(:) = missing
        do i=1,Nfloats
          if (FLT(i)%released) xx(i) = FLT(i)%Wdens
        enddo
        err = NF90_PUT_VAR(output_id,output_rid,xx,[1,output_record],[Nfloats,1])
        call cdf_error(err,'Error writing water density')
      endif

  else

    ! ... Write ASCII files
    ! ...
    date = num2date(Reference_time+time,model_time_units,model_time_calendar)

    do i=1,Nfloats
      iu = table_iu(i)
      axx(:) = missing
      if (FLT(i)%released) then
        axx( 1) = rad2deg*FLT(i)%x
        axx( 2) = rad2deg*FLT(i)%y
        axx( 3) = abs(FLT(i)%z)
        axx( 4) = FLT(i)%u
        axx( 5) = FLT(i)%v
        axx( 6) = FLT(i)%w
        axx( 7) = FLT(i)%Wtemp
        axx( 8) = FLT(i)%Wpsal
        axx( 9) = FLT(i)%Wdens
        axx(10) = FLT(i)%dist
      endif

      if (FLT(i)%released) then
        if (FLT(i)%indomain) then
          if (FLT(i)%floating) then
            code = 0
          else
            code = 2
          endif
        else
          code = 1
        endif
      else
        code = -1
      endif

      !fmt = '(T1,A20,2(1X,F9.4),1X,F6.1,X,2(1X,F7.3),1X,G9.2,3(F8.3,X),F6.1)'
      !write(iu,fmt) trim(date%iso()), axx
      fmt = '(T1,2(F9.4,1X),F6.1,X,A20,2(F9.3),1X,G10.3,3(F8.3,X),F7.1,I3)'
      write(iu,fmt) axx(1:3), trim(date%iso()), axx(4:), code
 
    enddo

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
      do i=1,Nfloats
        close(Table_iu(i))
        if (verb.ge.1) write(*,*) 'Successfully closed ASCII file ',trim(Table_filename(i))
      enddo
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

      x2 = displacement(xo,u1,dt14)
      u2 = RHS(x2,to+dt14)

      x3 = displacement(xo,0.5D0*(u1+u2),dt14)
      u3 = RHS(x3,to+dt14)

      x4 = displacement(xo,2.0D0*u3-u2,dt12)
      u4 = RHS(x4,to+dt12)

      x5 = displacement(xo,0.25D0*(u1+3.0_dp*u4),dt34)
      u5 = RHS(x5,to+dt34)

      x6 = displacement(xo,(-3.0_dp*u1+2.0_dp*u2+12.0_dp*u3-12.0_dp*u4+8.0_dp*u5)/7.0_dp,dt)
      u6 = RHS(x6,to+dt)

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
      if (noise_mult) then
        veps = (1.0 + noise_mu*randn())
        if (verb.ge.4) write(*,*) 'RK Multiplicative noise: ', veps
        un(:) = veps*un(:)
      endif
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
        if (verb.ge.4) write(*,*) 'RK Additive noise: ', vx, vy, vz
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

    real(dp) udt,vdt,wdt,yo,yn,cosyo,dx

    udt  = dt*u(1)     !   [seconds] x [meters/seconds]
    vdt  = dt*u(2)     !   [seconds] x [meters/seconds]
    wdt  = dt*u(3)     !   [seconds] x [meters/seconds]

    if (Spherical) then
      ! ... Spherical metrics
      ! ...
      udt = udt*IEarth    ! u*dt / REarth
      vdt = vdt*IEarth    ! v*dt / REarth
      yo  = xo(2)
      cosyo = cos(yo)
      yn = asin(sin(yo+vdt)*cos(udt))
      dx = atan2(sin(udt)*cosyo,cos(udt)-sin(yo)*sin(yn))
      xn(1)  = xo(1) + dx             ! Radians
      xn(2)  = yn                     ! Radians
      xn(3)  = xo(3) + wdt            ! meters
    else
      ! ... Cartesian metrics
      ! ...
      xn(1)  = xo(1) + udt             ! meters
      xn(2)  = xo(2) + vdt             ! meters
      xn(3)  = xo(3) + wdt             ! meters
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
  function emergence_model(F,rw) result (wo)

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
    wo = gravity*F%Psize*F%Psize*(1.0D0-F%Pdens/rw)/18.0D0/water_visc

  end function emergence_model
  ! ...
  ! =====================================================================
  ! ...
  real(dp) function Qswh(lon,lat,date) 

    ! ... Calcultates the clear-sky short wave flux W/m2
    ! ... as a function of the longitue, latitude and date.
    ! ... UTC time (not local time).
    ! ... From Bernie et al. (2007), Clim. Dyn. DOI: 10.1007/s00382-007-0249-6
    ! ...
    real(dp), intent(in)                         :: lon  ! radians
    real(dp), intent(in)                         :: lat  ! radians
    type(type_date), intent(in)                  :: date

    ! ... Local variables
    ! ...
    real(dp) t,Ha,dss,delta

    t = (date%hour*3600 + date%minute*60 + date%second)/86400.0D0
    Ha    = lat + pi*(t+t-1.0D0)
    dss   = date%yearday + 11       ! Days since winter solstice
    delta = -deg2rad*DEarth*cos(dpi*dss/365.25D0)
    Qswh  = max(0.0D0,Solar0*(sin(lat)*sin(delta)+cos(lat)*cos(delta)*cos(Ha)))

  end function Qswh
  ! ...
  ! =====================================================================
  ! ...
  real(dp) function eos(T,S,Z)
  ! ... Simplified equation of state
  ! ... linear case function of T only: rn_alpha<>0, other coefficients = 0
  ! ... linear eos function of T and S: rn_alpha and rn_beta<>0, other coefficients=0
  ! ... Vallis like equation: use default values of coefficients
  ! ... From NEMO 4.2.0

    real(dp), intent(in)                         :: T
    real(dp), intent(in)                         :: S
    real(dp), intent(in)                         :: Z

    ! ... Local variables
    ! ...
    real(dp) Ta,Sa,da

    Ta = T - 10.0D0
    Sa = S - 35.0D0
    da = -EOS_a0*(1.0D0 + 0.5D0*EOS_lam1*Ta + EOS_mu1*z)*Ta + &
          EOS_b0*(1.0D0 - 0.5D0*EOS_lam2*Sa - EOS_mu2*z)*Sa - &
          EOS_nu*Ta*Sa
    eos = EOS_rho0 + da

    return
  end function eos
  ! ...
  ! =====================================================================
  ! ...
  function time_interpol (G,F,t,xo) result(ff)

    type(type_ncgrid), intent(in)                        :: G
    real(dp), dimension(G%nx,G%ny,G%nz,2), intent(in)    :: F
    real(dp), intent(in)                                 :: t
    real(dp), dimension(ndims), intent(in)               :: xo
    real(dp)                                             :: ff

    ! ... Local variables
    ! ...
    integer, dimension(ndims)                            :: P
    real(dp) f1,f2,t1,t2


    P = G%locate(xo)
    if (G%Stationary) then
      ff  = G%interpol(F(:,:,:,1),xo,P,SingleLayer)
    else
      t1 = G%trec1
      t2 = G%trec2
      f1 = G%interpol(F(:,:,:,1),xo,P,SingleLayer)
      f2 = G%interpol(F(:,:,:,2),xo,P,SingleLayer)
      ff = f1 + (f2-f1)*(t-t1)/(t2-t1)
    endif

  end function time_interpol
  ! ...
  ! =====================================================================
  ! ...
  function dvm_model(xo,t) result(wo)

    real(dp), dimension(ndims), intent(in)       :: xo
    real(dp), intent(in)                         :: t
    real(dp)                                     :: wo

    ! ... Local variables
    ! ...
    !type(type_date) date
    !real(dp) lon,lat,z,tday,tdawn,tdusk
    !
    !lon = rad2deg*xo(1); lat = rad2deg*xo(1); z = xo(3)
    !date = num2date(Reference_time+t,model_time_units,model_time_calendar)
    !tday = anint(date%hour*3600.0D0 + date%minute*60.0D0 + date%second)
    !call SunRise_and_SunSet (lon,lat,date,tdawn,tdusk)
    !dvm_tdawn = tdawn
    !dvm_tdusk = tdusk

    !   ...... The previous calculations are done in the main time loop,
    !   ...... neglecting the beginning/ending of the vertical velocity
    !   ...... during a single time step. The smallest the time step, the
    !   ...... lower the error of this approximation.

    if ((model_daysec.ge.dvm_tdawn).and.(model_daysec.le.dvm_tdawn+dvm_tvm)) then
      wo = -dvm_w
    else if ((model_daysec.ge.dvm_tdusk).and.(model_daysec.le.dvm_tdusk+dvm_tvm)) then 
      wo = dvm_w
    else
      wo = 0.0D0
    endif
    
  end function dvm_model
  ! ...
  ! =====================================================================
  ! ...
end module module_model
