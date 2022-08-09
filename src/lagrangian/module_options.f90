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

module module_options

use module_types
use module_constants
use module_tools
use module_lineargs
use module_help

use module_alm
use module_forcing
use module_float
use module_model
use module_fitting

implicit none

logical                                          :: fhlp = .False.
type(type_help)                                  :: HLP

contains

subroutine options

  ! ... Local variables
  ! ...
  logical                                        :: WithVerb  = .False.
  logical                                        :: WithXmin  = .False.
  logical                                        :: WithXmax  = .False.
  logical                                        :: WithYmin  = .False.
  logical                                        :: WithYmax  = .False.
  logical                                        :: WithTname = .False.
  logical                                        :: WithTfnal = .False.
  logical                                        :: WithDt    = .False.
  logical                                        :: WithKH0   = .False.
  logical                                        :: WithKV0   = .False.
  logical                                        :: WithKH1   = .False.
  logical                                        :: WithKV1   = .False.
  logical                                        :: WithMVmin = .False.
  logical                                        :: WithAlfa  = .False.
!  logical                                        :: WithAlfa1 = .False.
!  logical                                        :: WithAlfa2 = .False.
!  logical                                        :: WithAlfa3 = .False.
  logical                                        :: WithRx    = .False.
  logical                                        :: WithRy    = .False.
  logical                                        :: WithRz    = .False.
  logical                                        :: WithRt    = .False.
  logical                                        :: WithA11   = .False.
  logical                                        :: WithA12   = .False.
  logical                                        :: WithA21   = .False.
  logical                                        :: WithA22   = .False.
  logical                                        :: WithVisc  = .False.
  logical                                        :: WithSaveP = .False.
  logical                                        :: WithWDepth= .False.
  logical                                        :: WithUterm = .False.
  logical                                        :: WithVterm = .False.
  logical                                        :: WithWterm = .False.
  logical                                        :: WithXterm = .False.
  logical                                        :: WithYterm = .False.
  logical                                        :: WithRterm = .False.
  logical                                        :: WithDVM   = .False.
  logical                                        :: WithRfn   = .False.
  logical                                        :: WithWind  = .False.
  logical                                        :: WithBeta  = .False.
  logical                                        :: WithTheta = .False.

  integer                                        :: na,i
  real(dp), dimension(:), allocatable            :: AAA
  character(len=maxlen)                          :: word
  character(len=400)                             :: OUlist=''
  character(len=400)                             :: OVlist=''
  character(len=400)                             :: OWlist=''
  character(len=400)                             :: OTlist=''
  character(len=400)                             :: OSlist=''
  character(len=400)                             :: ORlist=''
  character(len=400)                             :: OClist=''
  character(len=400)                             :: AUlist=''
  character(len=400)                             :: AVlist=''
  character(len=400)                             :: GRlist=''
  character(len=400)                             :: uplist=''
  character(len=400)                             :: WDlist=''
  character(len=400)                             :: DVMlist=''
  character(len=400)                             :: FITlist=''
  character(len=400)                             :: Alphalist=''


    ! ... Fill in the help information
    ! ...
    call program_help(HLP,VERSION,AUTHOR)

    na = lineargs()
    if (na.eq.0) call crash('Use option --help for options')

    ! ... Check for help
    ! ...
    call linearg('-he',fhlp)
    call linearg('--he',fhlp)

    if (fhlp) then
      call HLP%write()
      stop -1
    endif

    ! ... Verbose level
    ! ...
    call linearg('-verbose',WithVerb,verb)
    call linearg('-verb',WithVerb,verb)

    call linearg('-OU',WithUterm,OUlist)
    call linearg('-OV',WithVterm,OVlist)
    call linearg('-OW',WithWterm,OWlist)
    call linearg('-OR',WithRterm,ORlist)
    call linearg('-AU',WithXterm,AUlist)
    call linearg('-AV',WithYterm,AVlist)

    call linearg('-OT',WithOT,OTlist)
    call linearg('-OS',WithOS,OSlist)
    call linearg('-OC',WithOC,OClist)
    call linearg('-makegrid',WithMakeGrid,GRlist)

    call linearg('-cart',Cartesian)
    Spherical = .not.Cartesian

    call linearg('-eos',WithEOS)
    if (WithEOS.and.WithRterm) call crash('Incompatible options -eos and -OR')
    if (WithEOS.and..not.WithOT) call crash('Option -eos requires at least -OT option')

    model_density = WithRterm.or.WithEOS

    if (count([WithUterm,WithVterm]).ne.2) call crash('Horizontal currents not specified')
    if (count([WithXterm,WithYterm]).eq.1) call crash('Both wind components are required')
    if (WithXterm) Winds = .True.

    ! ... OU, OV, OW, OT, OS, OR, OC options
    ! ...
    if (WithUterm) then
      OUfilename = token_read(OUlist,'file')
      if (len_trim(OUfilename).ne.0) WithOU = .True.
      word = token_read(OUlist,'var'); if (len_trim(word).gt.0) OUvname = trim(word)
      word = token_read(OUlist,'x');   if (len_trim(word).gt.0) OUxname = trim(word)
      word = token_read(OUlist,'y');   if (len_trim(word).gt.0) OUyname = trim(word)
      word = token_read(OUlist,'z');   if (len_trim(word).gt.0) OUzname = trim(word)
      word = token_read(OUlist,'t');   if (len_trim(word).gt.0) OUtname = trim(word)
      word = token_read(OUlist,'val')
      if (len_trim(word).gt.0) then
        model_fixed_ou = .True.
        read(word,*) model_value_ou
      endif
      uplist = uppercase(OUlist)
      if (index(uplist,'CLIM').gt.0) OUClim = .True.
      if (index(uplist,'CLIM').gt.0) GOU%Climatology = .True.
      if (count([model_fixed_ou,WithOU]).eq.0) call crash('For -OU, missing file or val tokens')
    endif

    if (WithVterm) then
      OVfilename = token_read(OVlist,'file')
      if (len_trim(OVfilename).ne.0) WithOV = .True.
      word = token_read(OVlist,'var'); if (len_trim(word).gt.0) OVvname = trim(word)
      word = token_read(OVlist,'x');   if (len_trim(word).gt.0) OVxname = trim(word)
      word = token_read(OVlist,'y');   if (len_trim(word).gt.0) OVyname = trim(word)
      word = token_read(OVlist,'z');   if (len_trim(word).gt.0) OVzname = trim(word)
      word = token_read(OVlist,'t');   if (len_trim(word).gt.0) OVtname = trim(word)
      word = token_read(OVlist,'val')
      if (len_trim(word).gt.0) then
        model_fixed_ov = .True.
        read(word,*) model_value_ov
      endif
      uplist = uppercase(OVlist)
      if (index(uplist,'CLIM').gt.0) OVClim = .True.
      if (index(uplist,'CLIM').gt.0) GOV%Climatology = .True.
      if (count([model_fixed_ov,WithOV]).eq.0) call crash('For -OV, missing file or val tokens')
    endif
    

    if (WithWterm) then
      OWfilename = token_read(OWlist,'file')
      if (len_trim(OWfilename).ne.0) WithOW = .True.
      word = token_read(OWlist,'var'); if (len_trim(word).gt.0) OWvname = trim(word)
      word = token_read(OWlist,'x');   if (len_trim(word).gt.0) OWxname = trim(word)
      word = token_read(OWlist,'y');   if (len_trim(word).gt.0) OWyname = trim(word)
      word = token_read(OWlist,'z');   if (len_trim(word).gt.0) OWzname = trim(word)
      word = token_read(OWlist,'t');   if (len_trim(word).gt.0) OWtname = trim(word)
      if (len_trim(word).gt.0) then
        model_fixed_ow = .True.
        read(word,*) model_value_ow
      endif
      uplist = uppercase(OWlist)
      if (index(uplist,'CLIM').gt.0) OWClim = .True.
      if (index(uplist,'CLIM').gt.0) GOW%Climatology = .True.
      if (count([model_fixed_ow,WithOW]).eq.0) call crash('For -OW, missing file or val tokens')
    endif

    if (WithXterm) then
      AUfilename = token_read(AUlist,'file')
      if (len_trim(AUfilename).ne.0) WithAU = .True.
      word = token_read(AUlist,'var'); if (len_trim(word).gt.0) AUvname = trim(word)
      word = token_read(AUlist,'x');   if (len_trim(word).gt.0) AUxname = trim(word)
      word = token_read(AUlist,'y');   if (len_trim(word).gt.0) AUyname = trim(word)
      word = token_read(AUlist,'z');   if (len_trim(word).gt.0) AUzname = trim(word)
      word = token_read(AUlist,'t');   if (len_trim(word).gt.0) AUtname = trim(word)
      word = token_read(AUlist,'val')
      if (len_trim(word).gt.0) then
        model_fixed_au = .True.
        read(word,*) model_value_au
      endif
      uplist = uppercase(AUlist)
      if (index(uplist,'CLIM').gt.0) AUClim = .True.
      if (index(uplist,'CLIM').gt.0) GAU%Climatology = .True.
      if (count([model_fixed_au,WithAU]).eq.0) call crash('For -AU, missing file or val tokens')
    endif

    if (WithYterm) then
      AVfilename = token_read(AVlist,'file')
      if (len_trim(AVfilename).ne.0) WithAV = .True.
      word = token_read(AVlist,'var'); if (len_trim(word).gt.0) AVvname = trim(word)
      word = token_read(AVlist,'x');   if (len_trim(word).gt.0) AVxname = trim(word)
      word = token_read(AVlist,'y');   if (len_trim(word).gt.0) AVyname = trim(word)
      word = token_read(AVlist,'z');   if (len_trim(word).gt.0) AVzname = trim(word)
      word = token_read(AVlist,'t');   if (len_trim(word).gt.0) AVtname = trim(word)
      word = token_read(AVlist,'val')
      if (len_trim(word).gt.0) then
        model_fixed_av = .True.
        read(word,*) model_value_av
      endif
      uplist = uppercase(AVlist)
      if (index(uplist,'CLIM').gt.0) AVClim = .True.
      if (index(uplist,'CLIM').gt.0) GAV%Climatology = .True.
      if (count([model_fixed_av,WithAV]).eq.0) call crash('For -AV, missing file or val tokens')
    endif

    if (WithRterm) then
      ORfilename = token_read(ORlist,'file')
      if (len_trim(ORfilename).ne.0) WithOR = .True.
      word = token_read(ORlist,'var'); if (len_trim(word).gt.0) ORvname = trim(word)
      word = token_read(ORlist,'x');   if (len_trim(word).gt.0) ORxname = trim(word)
      word = token_read(ORlist,'y');   if (len_trim(word).gt.0) ORyname = trim(word)
      word = token_read(ORlist,'z');   if (len_trim(word).gt.0) ORzname = trim(word)
      word = token_read(ORlist,'t');   if (len_trim(word).gt.0) ORtname = trim(word)
      word = token_read(ORlist,'val')
      if (len_trim(word).gt.0) then
        word = uppercase(word)
        if (is_numeric(word)) then
          read(word,*) model_value_rho
          if (verb.ge.3) write(*,*) 'Constant water density: ', model_value_rho
          water_density_method = 0
        else
          if (verb.ge.3) write(*,*) 'Water density from analytical model'
          !WithOR = .False.
          water_density_method = 1
        endif
      else if (WithOR) then
        if (verb.ge.3) write(*,*) 'Water density from file'
        water_density_method = 2
      else
        call crash('Invalid density options')
      endif
    endif


    if (WithOT) then
      OTfilename = token_read(OTlist,'file')
      word = token_read(OTlist,'var'); if (len_trim(word).gt.0) OTvname = trim(word)
      word = token_read(OTlist,'x');   if (len_trim(word).gt.0) OTxname = trim(word)
      word = token_read(OTlist,'y');   if (len_trim(word).gt.0) OTyname = trim(word)
      word = token_read(OTlist,'z');   if (len_trim(word).gt.0) OTzname = trim(word)
      word = token_read(OTlist,'t');   if (len_trim(word).gt.0) OTtname = trim(word)
    endif
    if (WithOS) then
      OSfilename = token_read(OSlist,'file')
      word = token_read(OSlist,'var'); if (len_trim(word).gt.0) OSvname = trim(word)
      word = token_read(OSlist,'x');   if (len_trim(word).gt.0) OSxname = trim(word)
      word = token_read(OSlist,'y');   if (len_trim(word).gt.0) OSyname = trim(word)
      word = token_read(OSlist,'z');   if (len_trim(word).gt.0) OSzname = trim(word)
      word = token_read(OSlist,'t');   if (len_trim(word).gt.0) OStname = trim(word)
    endif

    if (WithOC) then
      OCfilename = token_read(OClist,'file')
      word = token_read(OClist,'var'); if (len_trim(word).gt.0) OCvname = trim(word)
      word = token_read(OClist,'x');   if (len_trim(word).gt.0) OCxname = trim(word)
      word = token_read(OClist,'y');   if (len_trim(word).gt.0) OCyname = trim(word)
      word = token_read(OClist,'z');   if (len_trim(word).gt.0) OCzname = trim(word)
      word = token_read(OClist,'t');   if (len_trim(word).gt.0) OCtname = trim(word)
    endif

    if (WithMakeGrid) then
      word = token_read(GRlist,'west'); if (len_trim(word).gt.0) read(word,*) alm_xmin
      word = token_read(GRlist,'east'); if (len_trim(word).gt.0) read(word,*) alm_xmax
      word = token_read(GRlist,'south'); if (len_trim(word).gt.0) read(word,*) alm_ymin
      word = token_read(GRlist,'north'); if (len_trim(word).gt.0) read(word,*) alm_ymax
      word = token_read(GRlist,'depth'); if (len_trim(word).gt.0) read(word,*) alm_depth
      word = token_read(GRlist,'dx'); if (len_trim(word).gt.0) read(word,*) alm_dx
      word = token_read(GRlist,'dy'); if (len_trim(word).gt.0) read(word,*) alm_dy
      word = token_read(GRlist,'dz'); if (len_trim(word).gt.0) read(word,*) alm_dz
      alm_zmin = -abs(alm_depth)
      alm_zmax = 0.0D0
    endif
       
    WithClim = .False.
    if (any([OUClim,OVClim,OWClim,AUClim,AVClim])) WithClim = .True.

    ! ... Model crop:
    ! ...
    call linearg('-xmin',WithXmin,forcing_xmin)
    call linearg('-xmax',WithXmax,forcing_xmax)
    call linearg('-ymin',WithYmin,forcing_ymin)
    call linearg('-ymax',WithYmax,forcing_ymax)
    if (WithXmin) forcing_xmin = deg2rad*forcing_xmin
    if (WithXmax) forcing_xmax = deg2rad*forcing_xmax
    if (WithYmin) forcing_ymin = deg2rad*forcing_ymin
    if (WithYmax) forcing_ymax = deg2rad*forcing_ymax

    call linearg('-mu',noise_mult,noise_mu)
    call linearg('-kh0',WithKH0,noise_KH0)
    call linearg('-kv0',WithKV0,noise_KV0)
    call linearg('-kh1',WithKH1,noise_KH1)
    call linearg('-kv1',WithKV1,noise_KV1)

    ! ... Velocity factor
    ! ...
    !call linearg('-alpha1',WithAlfa1,alpha(1))
    !call linearg('-alpha2',WithAlfa2,alpha(2))
    !call linearg('-alpha2',WithAlfa3,alpha(3))
    call linearg('-alpha',WithAlfa,Alphalist)
    if (WithAlfa) then
      if (Alphalist(1:1).eq.'[') then  
        ! ... A vector has been entered
        AAA = ReadVector(trim(Alphalist))
        na = size(AAA)
        if (na.eq.1) then
          alpha(:) = AAA(1)
        else if (na.eq.2) then
          alpha(1) = AAA(1)
          alpha(2) = AAA(2)
        else if (na.eq.3) then
          alpha(1) = AAA(1)
          alpha(2) = AAA(2)
          alpha(3) = AAA(3)
        endif
      else
        ! ... A single value has been entered
        read(Alphalist,*) alpha(1)
        read(Alphalist,*) alpha(2)
        read(Alphalist,*) alpha(3)
      endif
    endif
    
    ! ... Trajectory name and final point
    ! ...
    call linearg('-traj',WithTname,trajectory_name)
    call linearg('-end',WithTfnal,trajectory_final)
    call linearg('-saveper',WithSaveP,save_period)
    if ((trajectory_final.eq.'NONE').or.(trajectory_final.eq.'NULL')) trajectory_final = ''

    ! ... Pereiro, 2019, minimum thrsthold velovity parameter
    ! ...
    call linearg('-velmin',WithMVmin,model_velmin)
    model_velmin = max(0.0D0,model_velmin)

    ! ... Float release
    ! ...
    call linearg('-rel',Release_by_file,Release_file)
    call linearg('-xo',WithReleaseXo,Release_xo)
    call linearg('-yo',WithReleaseYo,Release_yo)
    call linearg('-zo',WithReleaseZo,Release_zo)
    call linearg('-to',WithReleaseTime,Release_time)
    call linearg('-ro',WithReleaseRho,Release_rho)
    call linearg('-so',WithReleaseSize,Release_size)
    call linearg('-save_release',WithRfn,Release_SaveFile)

    if ( WithReleaseRho.or.WithReleaseSize) Particle_buoyant = .True.
    if (Particle_buoyant) then
      if (.not.model_density) call crash('Buoyant particle needs density specification')
    endif

    call linearg('-vis',WithVisc,water_visc)

    ! ... Random floats
    ! ...
    call linearg('-random',WithRandom,Nrandom)
    call linearg('-cloud',WithRandom,Nrandom)
    call linearg('-ensemble',WithRandom,Nrandom)
    call linearg('-rx',WithRx,Radius_x)
    call linearg('-ry',WithRy,Radius_y)
    call linearg('-rz',WithRz,Radius_z)
    call linearg('-rt',WithRt,Radius_t)


    Release_by_pos = WithReleaseXo .or. WithReleaseYo 

    if (Release_by_file.and.Release_by_pos) then
      call crash('Incompatible release options')
    endif
    !SingleLayer = .True.

    ! ... Reverse simulation
    ! ...
    call linearg('-reverse',reverse)
    call linearg('-backward',reverse)

    ! ... Model time step
    ! ...
    call linearg('-dt',WithDt,model_dt)
    model_dt = anint(model_dt)

    call linearg('-from',WithTini,StrgTini)
    call linearg('-for',WithTlen,StrgTlen)
    call linearg('-during',WithTlen,StrgTlen)


    ! ... Retrieve user-defined model simulation legth
    ! ...
    if (WithTlen) then
      if (is_numeric(StrgTlen)) then
        ! ... no units provided. Default units: "days"
        ! ...
        read(StrgTlen,*) UserTlen
        USerTlen = anint(UserTlen * 86400)              ! Convert to seconds
      else
        StrgTlen = uppercase(StrgTlen)
        i = index(StrgTlen,'D')
        if (i.gt.0) then
          read(StrgTlen(1:i-1),*) UserTlen
          USerTlen = anint(UserTlen * 86400)           ! Convert to seconds
        else 
          i = index(StrgTlen,'H')
          if (i.gt.0) then
            read(StrgTlen(1:i-1),*) UserTlen
            USerTlen = anint(UserTlen * 3600)          ! Convert to seconds
          else
            i = index(StrgTlen,'M')
            if (i.gt.0) then
              read(StrgTlen(1:i-1),*) UserTlen
              USerTlen = anint(UserTlen * 60)          ! Convert to seconds
            else
              i = index(StrgTlen,'S')
              if (i.gt.0) then
                read(StrgTlen(1:i-1),*) UserTlen ! Already in seconds
                UserTlen = anint(UserTlen)
              else
                call crash('Invalid units in simulation length option')
              endif
            endif
          endif
        endif
      endif
    endif

    ! ... Wind method and parameters
    ! ...
    call linearg('-Wind',WithWind,WDlist)
    if (WithWind) then
      WDlist = lowercase(WDlist)
      word = token_read(WDlist,'a')
      if (len_trim(word).gt.0) then
        WindResponse = .True.
        AAA = ReadVector(trim(word))
        if (size(AAA).ne.4) call crash('Incompatible number of A values')
        A11 = AAA(1)
        A12 = AAA(2)
        A21 = AAA(3)
        A22 = AAA(4)
      endif
      word = token_read(WDlist,'beta')
      if (len_trim(word).gt.0) then
        WindDriven = .True.
        Withbeta   = .True.
        read(word,*) WDriven_beta
      endif
      word = token_read(WDlist,'theta')
      if (len_trim(word).gt.0) then
        WindDriven  = .True.
        Withtheta   = .True.
        read(word,*) WDriven_theta
      endif
      word = token_read(WDlist,'depth'); if (len_trim(word).gt.0) read(word,*) WindDepth
    endif
    if (WindResponse.and.WindDriven) call crash('Incompatible use of A, beta and theta parameters')
    if (count([Withbeta,WithTheta]).eq.1) call crash('Both beta and theta values required')
    if (Winds) then
      if (count([WindResponse,WindDriven]).eq.0) call crash('-Wind forcing requires specification of wind paramteres')
    endif

    ! ... Wind Response matrix A11, A12, A21, A22:
    ! ...
    !call linearg('-Winddriven',WindDriven,WDlist)
    !if (WindDriven) then
    !  word = token_read(WDlist,'alpha='); if (len_trim(word).gt.0) read(word,*) WDriven_alpha
    !  word = token_read(WDlist,'beta='); if (len_trim(word).gt.0) read(word,*) WDriven_beta
    !  WindDepth = 0.0D0
    !  if (verb.ge.4) write(*,*) 'WDriven_alpha, WDriven_beta : ', WDriven_alpha, WDriven_beta
    !  if (verb.ge.4) write(*,*) 'WindDepth : ', WindDepth
    !else
    !  call linearg('-winddepth',WithWDepth,WindDepth)
    !  call linearg('-a11',WithA11,A11)
    !  call linearg('-a11',WithA11,A11)
    !  call linearg('-a12',WithA12,A12)
    !  call linearg('-a21',WithA21,A21)
    !  call linearg('-a22',WithA22,A22)
    !endif

    WindDepth = -abs(WindDepth)

    ! ... Options for the Equation of state
    ! ... The model can calculate the density of the water, independently of the 
    ! ... application, or not, of the buoyancy term. The buoyancy term is activated
    ! ... if the option -ro (for setting the density of the particle). This activates
    ! ... the flag: Particle_buoyant
    ! ... The options to track the density are activated via the option -OR or via 
    ! ... the option -EOS, for using temperature and salinity. If any of these two
    ! ... are used, it activates the flag model_density.
    ! ... 
    if (model_density) then
      if (WithRterm) then
        !if (verb.ge.1) write(*,*) 'Density information provided'
      else 
        if (WithOS) then
          water_density_method = 4
          !if (verb.ge.1) write(*,*) 'Simplified Equation of State'
        else
          water_density_method = 3
          !if (verb.ge.1) write(*,*) 'Linear Temperature Equation of State'
          EOS_b0   = 0.0D0  ! halin. expan.
          EOS_lam1 = 0.0D0  ! T2 cabbeling
          EOS_lam2 = 0.0D0  ! S2 cabbeling
          EOS_nu   = 0.0D0  ! TS cabbeling
          EOS_mu1  = 0.0D0  ! T thermobaric
          EOS_mu2  = 0.0D0  ! S thermobaric
        endif
      endif
    endif

    ! ... Options for Diel Vertical Motion (DVM)
    ! ...
    call linearg('-DVM',WithDVM,DVMlist)
    if (WithDVM) then
      Particle_dvm = .True.
      word = token_read(DVMlist,'zday'); if (len_trim(word).gt.0) read(word,*) dvm_zday 
      word = token_read(DVMlist,'znight'); if (len_trim(word).gt.0) read(word,*) dvm_znight 
      word = token_read(DVMlist,'tvm'); if (len_trim(word).gt.0) read(word,*) dvm_tvm
    endif


    ! ... Fitting options
    ! ...
    call linearg('-fitting',option_fitting,FITlist)

    if (option_fitting) then
      option_model = .False.

      ! ... Error messages
      ! ...
      if (.not.WithTname) call crash('Model fitting requires option -trajectory')

      ! ... Warning messages
      ! ...
      if (reverse) then
        if (verb.ge.2) write(*,*) 'WARNING: Overriding option -reverse'
        reverse = .False.
      endif

      word = token_read(FITlist,'vmin'); if (len_trim(word).gt.0) read(word,*) fit_vmin
      word = token_read(FITlist,'vmax'); if (len_trim(word).gt.0) read(word,*) fit_vmax
      word = token_read(FITlist,'out'); if (len_trim(word).gt.0)  read(word,*) fit_fout
      word = token_read(FITlist,'first')
      if (len_trim(word).gt.0)  then
        AAA = ReadVector(trim(word))
        if (size(AAA).ne.5) call crash('Incompatible number of initial parameters')
        fit_FGp1 = AAA(1)
        fit_FGp2 = AAA(2)
        fit_FGp3 = AAA(3)
        fit_FGp4 = AAA(4)
        fit_FGp5 = AAA(5)
      endif
      word = token_read(FITlist,'adjust')
      if (len_trim(word).gt.0)  then
        AAA = ReadVector(trim(word))
        if (size(AAA).ne.5) call crash('Incompatible number of initial parameters')
        if (AAA(1).ge.0.5) then
          fit_dop1 = 1.0D0
        else
          fit_dop1 = 0.0D0
        endif
        if (AAA(2).ge.0.5) then
          fit_dop2 = 1.0D0
        else
          fit_dop2 = 0.0D0
        endif
        if (AAA(3).ge.0.5) then
          fit_dop3 = 1.0D0
        else
          fit_dop3 = 0.0D0
        endif
        if (AAA(4).ge.0.5) then
          fit_dop4 = 1.0D0
        else
          fit_dop4 = 0.0D0
        endif
        if (AAA(5).ge.0.5) then
          fit_dop5 = 1.0D0
        else
          fit_dop5 = 0.0D0
        endif
      endif
    
    endif

    if (option_model) then
      ! ... Check that release information has been provided
      ! ...
      if (.not.Release_by_file.and. &
          .not.WithRandom.and. &
          .not.Release_by_pos) call crash('Missing release information')

      ! ... If a climatology field, we require the user to use the options -from and -for
      ! ...
      if (WithClim) then
        if (.not.WithTini.or..not.WithTlen) call crash('Climatology requires options -from and -for')
      endif
    endif


    ! ... Check options
    ! ...
    call checkopts()

  end subroutine options
  ! ...
  ! ====================================================================
  ! ...
  subroutine program_help(HLP,version,author)

  type(type_help), intent(inout)                         :: HLP
  character(len=*), intent(in)                           :: version
  character(len=*), intent(in)                           :: author

  ! ... The help
  ! ...
  HLP%version = version
  HLP%progname = 'ALM'
  HLP%author = author
  call HLP%set_summary('Reads the zonal and meridional velocity components &
    &from a NetCDF file and calculates trajectories from the time-evolving &
    &velocity field or streamlines for stationary cases. If a list of &
    &floats is not provided, NFLOATS positions will be randomly generated. The &
    &number and position of the floats may be read from an ASCII file &
    &(LON, LAT, DEPTH, RELEASE_TIME [,...]) or passed though command line &
    &using options -xo and -yo (and optionally -to or -do). In the later case, &
    &a random cloud of NFLOATS may also be generated using the option -rand. &
    &The number of internal time steps can be modified using the &
    &option -idt, that specifies the time step of the internal loop. &
    &The program writes a trajectory file (NetCDF) and the final position of &
    &the floats (ASCII). The names can be specified using options -trajectory &
    &and -end, respctively.')
  call HLP%add_option ('-OU token=value [token=value...]','Input ocean U field (required)','')
  call HLP%add_option ('-OV token=value [token=value...]','Input ocean V field (optional)','')
  call HLP%add_option ('-OT token=value [token=value...]','Input ocean temperature field (optional)','')
  call HLP%add_option ('-OS token=value [token=value...]','Input ocean salinity field (optional)','')
  call HLP%add_option ('-OR token=value [token=value...]','Input ocean density field (optional)','')
  call HLP%add_option ('-AU token=value [token=value...]','Input atmosphere U field (optional)','')
  call HLP%add_option ('-AV token=value [token=value...]','Input atmosphere V field (optional)','')
  call HLP%add_option ('-A11        value','Component a11 of the atmosphere response matrix','0.0')
  call HLP%add_option ('-A12        value','Component a12 of the atmosphere response matrix','0.0')
  call HLP%add_option ('-A21        value','Component a21 of the atmosphere response matrix','0.0')
  call HLP%add_option ('-A22        value','Component a22 of the atmosphere response matrix','0.0')
  call HLP%add_option ('-release    filename ','Initial position release file name. &
   &It must exist if no initial coordinates are specified (options -xo and -yo)','')
  call HLP%add_option ('-from       INITIAL_DATE','Date at which the Lagrangian simulation will start','')
  call HLP%add_option ('-for        TIME_PERIOD','Length of the Lagrangian simulation','')
  call HLP%add_option ('-dt         DT (in seconds)','Runge-Kutta time step value','600')
  call HLP%add_option ('-reverse           ','Perform backward integration','')
  call HLP%add_option ('-trajectory filename ','Input/Output trajectory file','out.nc')
  call HLP%add_option ('-endpos     filename ','Output final position file','release.out')
  call HLP%add_option ('-xmin       MIN_LONGITUDE','Option to crop the computational domain','')
  call HLP%add_option ('-xmax       MAX_LONGITUDE','Option to crop the computational domain','')
  call HLP%add_option ('-ymin       MIN_LATITUDE','Option to crop the computational domain','')
  call HLP%add_option ('-ymax       MAX_LATITUDE','Option to crop the computational domain','')
  call HLP%add_option ('-xo         XO','Optional value of the float initial position','')
  call HLP%add_option ('-yo         YO','Optional value of the float initial position','')
  call HLP%add_option ('-zo         ZO','Optional value of the float initial position','0.0')
  call HLP%add_option ('-to         TO/DATE','Optional value of the float initial &
   &release time (seconds after initial simulation time)','0')
  call HLP%add_option ('-ro         RHO','Optional value of the float density. It requires &
    &the use of option -OR to specify water density','')
  call HLP%add_option ('-so         SIZE','Optional value of the float diameter','')
  call HLP%add_option ('-velmin     MIN_THRESHOLD','Minimum velocity threshold. &
    &Set it to zero for pure random walk motion','1D-5')
  call HLP%add_option ('-mu         value','Non-dimensioanl amplitude of the gaussian multiplicative noise','0')
  call HLP%add_option ('-va         value','Amplitude of the gaussian istropic velocity fluctuation','0')
  call HLP%add_option ('-alpha      value','Non-dimensional ocean velocity multiplicator','1.0')
  call HLP%add_option ('-rand       NFLOATS','Option to request a simulation with NFLOATS randomly &
    &generated floats per each release position.','1')
  call HLP%add_option ('-Rx         RADIUS_X','Longitude radius for releasing random floats &
   &around specified release location','0.02')
  call HLP%add_option ('-Ry         RADIUS_Y','Latitude radius for releasing random floats &
   &around specified release location','0.02')
  call HLP%add_option ('-Rz         RADIUS_Z','Depth radius for releasing random floats &
   &around specified release location','0.0')
  call HLP%add_option ('-Rt         RADIUS_T','Time radius for releasing random floats &
   &around specified release location','0.0')
  call HLP%add_option ('-verbose    VERBOSE_LEVEL','To increase output verbosity. (0=Quiet,...3=Maximum verbose)','1')
  call HLP%add_option ('--options   filename','To read the commandline options from a file','')
  call HLP%add_option ('--help','To show this help','')
  
  call HLP%set_example('alm -OU file=roms.nc u=u x=lon y=lat t=time &
   &-V file=roms.nc v=v x=lon y=lat t=time -release release.inp &
   &-trajectory float.nc -end release.out')
  
  end subroutine program_help
  ! ...
  ! ====================================================================
  ! ...
end module module_options
