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
  logical                                        :: WithAu    = .False.
  logical                                        :: WithAv    = .False.
  logical                                        :: WithAw    = .False.
  logical                                        :: WithRx    = .False.
  logical                                        :: WithRy    = .False.
  logical                                        :: WithRz    = .False.
  logical                                        :: WithRt    = .False.
  logical                                        :: withA11   = .False.
  logical                                        :: withA12   = .False.
  logical                                        :: withA21   = .False.
  logical                                        :: withA22   = .False.

  integer                                        :: na,i
  character(len=maxlen)                          :: word
  character(len=maxlen)                          :: rhovalue
  character(len=400)                             :: OUlist=''
  character(len=400)                             :: OVlist=''
  character(len=400)                             :: OWlist=''
  character(len=400)                             :: OTlist=''
  character(len=400)                             :: OSlist=''
  character(len=400)                             :: ORlist=''
  character(len=400)                             :: OClist=''
  character(len=400)                             :: AUlist=''
  character(len=400)                             :: AVlist=''


    ! ... Fill in the help information
    ! ...
    call program_help(HLP,VERSION,AUTHOR)

    na = lineargs()
    if (na.eq.0) call crash('Use option --help for options')

    ! ... Check for help
    ! ...
    call argflg('-he',fhlp)
    call argflg('--he',fhlp)

    if (fhlp) then
      call HLP%write()
      stop -1
    endif

    ! ... Verbose level
    ! ...
    call argint('-verbose',WithVerb,verb)
    call argint('-verb',WithVerb,verb)

    call arglst('-OU',WithOU,OUlist)
    call arglst('-OV',WithOV,OVlist)
    call arglst('-OW',WithOW,OWlist)
    call arglst('-OT',WithOT,OTlist)
    call arglst('-OS',WithOS,OSlist)
    call arglst('-OR',WithOR,ORlist)
    call arglst('-OC',WithOC,OClist)
    call arglst('-AU',WithAU,AUlist)
    call arglst('-AV',WithAV,AVlist)

    call argflg('-cart',Cartesian)
    call argflg('-Cart',Cartesian)
    call argflg('-CART',Cartesian)
    Spherical = .not.Cartesian

    OUfilename = 'fakebasin.nc'
    OUxname    = 'lon'
    OUyname    = 'lat'
    OUzname    = 'depth'
    OUtname    = 'time'
    OUvname    = 'u'

    OVfilename = 'fakebasin.nc'
    OVxname    = 'lon'
    OVyname    = 'lat'
    OVzname    = 'depth'
    OVtname    = 'time'
    OVvname    = 'v'

!  OWfilename = 'fakebasin.nc'
    OWxname    = 'lon'
    OWyname    = 'lat'
    OWzname    = 'depth'
    OWtname    = 'time'
    OWvname    = 'w'

    ! ... OU, OV, OW, OT, OS, OR, OC options
    ! ...
    if (WithOU) then
      OUfilename = token_read(OUlist,'file=')
      word = token_read(OUlist,'var='); if (len_trim(word).gt.0) OUvname = trim(word)
      word = token_read(OUlist,'x=');   if (len_trim(word).gt.0) OUxname = trim(word)
      word = token_read(OUlist,'y=');   if (len_trim(word).gt.0) OUyname = trim(word)
      word = token_read(OUlist,'z=');   if (len_trim(word).gt.0) OUzname = trim(word)
      word = token_read(OUlist,'t=');   if (len_trim(word).gt.0) OUtname = trim(word)
    endif
    if (WithOV) then
      OVfilename = token_read(OVlist,'file=')
      word = token_read(OVlist,'var='); if (len_trim(word).gt.0) OVvname = trim(word)
      word = token_read(OVlist,'x=');   if (len_trim(word).gt.0) OVxname = trim(word)
      word = token_read(OVlist,'y=');   if (len_trim(word).gt.0) OVyname = trim(word)
      word = token_read(OVlist,'z=');   if (len_trim(word).gt.0) OVzname = trim(word)
      word = token_read(OVlist,'t=');   if (len_trim(word).gt.0) OVtname = trim(word)
    endif
    if (WithOW) then
      OWfilename = token_read(OWlist,'file=')
      word = token_read(OWlist,'var='); if (len_trim(word).gt.0) OWvname = trim(word)
      word = token_read(OWlist,'x=');   if (len_trim(word).gt.0) OWxname = trim(word)
      word = token_read(OWlist,'y=');   if (len_trim(word).gt.0) OWyname = trim(word)
      word = token_read(OWlist,'z=');   if (len_trim(word).gt.0) OWzname = trim(word)
      word = token_read(OWlist,'t=');   if (len_trim(word).gt.0) OWtname = trim(word)
    endif

    if (WithOT) then
      OTfilename = token_read(OTlist,'file=')
      word = token_read(OTlist,'var='); if (len_trim(word).gt.0) OTvname = trim(word)
      word = token_read(OTlist,'x=');   if (len_trim(word).gt.0) OTxname = trim(word)
      word = token_read(OTlist,'y=');   if (len_trim(word).gt.0) OTyname = trim(word)
      word = token_read(OTlist,'z=');   if (len_trim(word).gt.0) OTzname = trim(word)
      word = token_read(OTlist,'t=');   if (len_trim(word).gt.0) OTtname = trim(word)
    endif
    if (WithOS) then
      OSfilename = token_read(OSlist,'file=')
      word = token_read(OSlist,'var='); if (len_trim(word).gt.0) OSvname = trim(word)
      word = token_read(OSlist,'x=');   if (len_trim(word).gt.0) OSxname = trim(word)
      word = token_read(OSlist,'y=');   if (len_trim(word).gt.0) OSyname = trim(word)
      word = token_read(OSlist,'z=');   if (len_trim(word).gt.0) OSzname = trim(word)
      word = token_read(OSlist,'t=');   if (len_trim(word).gt.0) OStname = trim(word)
    endif
    if (WithOR) then
      ORfilename = token_read(ORlist,'file=')
      rhovalue = token_read(ORlist,'value=')
      if (len_trim(rhovalue).gt.0) then
        rhovalue = uppercase(rhovalue)
        if(is_numeric(rhovalue)) then
          if (verb.ge.3) write(*,*) 'Constant water density: ', trim(rhovalue)
          WithOR = .False.
          water_density_method = 0
          read(rhovalue,*) water_rho
        else
          if (verb.ge.3) write(*,*) 'Water density from analytical model'
          WithOR = .False.
          water_density_method = 1
        endif
      else
        if (verb.ge.3) write(*,*) 'Water density from file'
        water_density_method = 2
        word = token_read(ORlist,'var='); if (len_trim(word).gt.0) ORvname = trim(word)
        word = token_read(ORlist,'x=');   if (len_trim(word).gt.0) ORxname = trim(word)
        word = token_read(ORlist,'y=');   if (len_trim(word).gt.0) ORyname = trim(word)
        word = token_read(ORlist,'z=');   if (len_trim(word).gt.0) ORzname = trim(word)
        word = token_read(ORlist,'t=');   if (len_trim(word).gt.0) ORtname = trim(word)
      endif
    endif

    if (WithOC) then
      OCfilename = token_read(OClist,'file=')
      word = token_read(OClist,'var='); if (len_trim(word).gt.0) OCvname = trim(word)
      word = token_read(OClist,'x=');   if (len_trim(word).gt.0) OCxname = trim(word)
      word = token_read(OClist,'y=');   if (len_trim(word).gt.0) OCyname = trim(word)
      word = token_read(OClist,'z=');   if (len_trim(word).gt.0) OCzname = trim(word)
      word = token_read(OClist,'t=');   if (len_trim(word).gt.0) OCtname = trim(word)
    endif
    if (WithAU) then
      AUfilename = token_read(AUlist,'file=')
      word = token_read(AUlist,'var='); if (len_trim(word).gt.0) AUvname = trim(word)
      word = token_read(AUlist,'x=');   if (len_trim(word).gt.0) AUxname = trim(word)
      word = token_read(AUlist,'y=');   if (len_trim(word).gt.0) AUyname = trim(word)
      word = token_read(AUlist,'z=');   if (len_trim(word).gt.0) AUzname = trim(word)
      word = token_read(AUlist,'t=');   if (len_trim(word).gt.0) AUtname = trim(word)
    endif
    if (WithAV) then
      AVfilename = token_read(AVlist,'file=')
      word = token_read(AVlist,'var='); if (len_trim(word).gt.0) AVvname = trim(word)
      word = token_read(AVlist,'x=');   if (len_trim(word).gt.0) AVxname = trim(word)
      word = token_read(AVlist,'y=');   if (len_trim(word).gt.0) AVyname = trim(word)
      word = token_read(AVlist,'z=');   if (len_trim(word).gt.0) AVzname = trim(word)
      word = token_read(AVlist,'t=');   if (len_trim(word).gt.0) AVtname = trim(word)
    endif

    ! ... Model crop:
    ! ...
    call argdbl('-xmin',WithXmin,forcing_xmin)
    call argdbl('-xmax',WithXmax,forcing_xmax)
    call argdbl('-ymin',WithYmin,forcing_ymin)
    call argdbl('-ymax',WithYmax,forcing_ymax)
    if (WithXmin) forcing_xmin = deg2rad*forcing_xmin
    if (WithXmax) forcing_xmax = deg2rad*forcing_xmax
    if (WithYmin) forcing_ymin = deg2rad*forcing_ymin
    if (WithYmax) forcing_ymax = deg2rad*forcing_ymax


    call argdbl('-KH0',WithKH0,noise_KH0)
    call argdbl('-Kh0',WithKH0,noise_KH0)
    call argdbl('-kh0',WithKH0,noise_KH0)
    call argdbl('-KV0',WithKV0,noise_KV0)
    call argdbl('-Kv0',WithKV0,noise_KV0)
    call argdbl('-kv0',WithKV0,noise_KV0)
    call argdbl('-KH1',WithKH1,noise_KH1)
    call argdbl('-Kh1',WithKH1,noise_KH1)
    call argdbl('-kh1',WithKH1,noise_KH1)
    call argdbl('-KV1',WithKV1,noise_KV1)
    call argdbl('-Kv1',WithKV1,noise_KV1)
    call argdbl('-kv1',WithKV1,noise_KV1)

    ! ... Velocity factor
    ! ...
    call argdbl('-au',WithAu,alpha(1))
    call argdbl('-av',WithAv,alpha(2))
    call argdbl('-aw',WithAw,alpha(3))
    
    ! ... Trajectory name and final point
    ! ...
    call argstr('-trajectory',WithTname,trajectory_name)
    call argstr('-endpos',WithTfnal,trajectory_final)

    ! ... Pereiro, 2019, minimum thrsthold velovity parameter
    ! ...
    call argdbl('-velmin',WithMVmin,model_velmin)
    model_velmin = max(0.0D0,model_velmin)

    ! ... Float release
    ! ...
    call argstr('-rel',Release_by_file,Release_file)
    call argdbl('-xo',WithReleaseXo,Release_xo)
    call argdbl('-yo',WithReleaseYo,Release_yo)
    call argdbl('-zo',WithReleaseZo,Release_zo)
    call argstr('-to',WithReleaseTime,Release_time)
    call argdbl('-ro',WithReleaseRho,Release_rho)
    call argdbl('-so',WithReleaseSize,Release_size)
    if (WithReleaseRho.and.water_density_method.lt.0) then
      call crash('No water density method has been specified.')
    endif

    if (water_density_method.ge.0.or. &
        WithReleaseRho.or. &
        WithReleaseSize) model_buoyancy = .True.

    ! ... Random floats
    ! ...
    call argint('-random',WithRandom,Nrandom)
    call argint('-cloud',WithRandom,Nrandom)
    call argint('-ensemble',WithRandom,Nrandom)
    call argdbl('-Rx',WithRx,Radius_x)
    call argdbl('-RX',WithRx,Radius_x)
    call argdbl('-rx',WithRx,Radius_x)
    call argdbl('-rX',WithRx,Radius_x)
    call argdbl('-Ry',WithRy,Radius_y)
    call argdbl('-RY',WithRy,Radius_y)
    call argdbl('-ry',WithRy,Radius_y)
    call argdbl('-rY',WithRy,Radius_y)
    call argdbl('-Rz',WithRz,Radius_z)
    call argdbl('-RZ',WithRz,Radius_z)
    call argdbl('-rz',WithRz,Radius_z)
    call argdbl('-rZ',WithRz,Radius_z)
    call argdbl('-Rt',WithRt,Radius_t)
    call argdbl('-RT',WithRt,Radius_t)
    call argdbl('-rt',WithRt,Radius_t)
    call argdbl('-rT',WithRt,Radius_t)


    Release_by_pos = WithReleaseXo .or. WithReleaseYo 

    if (Release_by_file.and.Release_by_pos) then
      call crash('Incompatible release options')
    endif
    if (Release_by_pos) then
      !if (.not.WithReleaseXo) call crash('Option -xo required')
      !if (.not.WithReleaseYo) call crash('Option -yo required')
    endif
    if (option_model) then
      if (.not.Release_by_file.and. &
          .not.WithRandom.and. &
          .not.Release_by_pos) call crash('Missing release information')
    else 
      if (.not.Release_by_file) call crash('Missing release information')
    endif

    !SingleLayer = .True.

    ! ... Reverse simulation
    ! ...
    call argflg('-reverse',reverse)
    call argflg('-backward',reverse)

    ! ... Model time step
    ! ...
    call argdbl('-dt',WithDt,model_dt)

    call argstr('-from',WithTini,StrgTini)
    call argstr('-for',WithTlen,StrgTlen)
    call argstr('-during',WithTlen,StrgTlen)


    ! ... Retrieve user-defined model simulation legth
    ! ...
    if (WithTlen) then
      if (is_numeric(StrgTlen)) then
        ! ... no units provided. Default units: "days"
        ! ...
        read(StrgTlen,*) UserTlen
        USerTlen = UserTlen * 86400              ! Convert to seconds
      else
        StrgTlen = uppercase(StrgTlen)
        i = index(StrgTlen,'D')
        if (i.gt.0) then
          read(StrgTlen(1:i-1),*) UserTlen
          USerTlen = UserTlen * 86400            ! Convert to seconds
        else 
          i = index(StrgTlen,'H')
          if (i.gt.0) then
            read(StrgTlen(1:i-1),*) UserTlen
            USerTlen = UserTlen * 3600           ! Convert to seconds
          else
            i = index(StrgTlen,'M')
            if (i.gt.0) then
              read(StrgTlen(1:i-1),*) UserTlen
              USerTlen = UserTlen * 60           ! Convert to seconds
            else
              i = index(StrgTlen,'S')
              if (i.gt.0) then
                read(StrgTlen(1:i-1),*) UserTlen ! Already in seconds
              else
                call crash('Invalid units in simulation length option')
              endif
            endif
          endif
        endif
      endif
    endif

    ! ... Wind Response matrix A11, A12, A21, A22:
    ! ...
    call argdbl('-a11',WithA11,A11)
    call argdbl('-A11',WithA11,A11)
    call argdbl('-a12',WithA12,A12)
    call argdbl('-A12',WithA12,A12)
    call argdbl('-a21',WithA21,A21)
    call argdbl('-A21',WithA21,A21)
    call argdbl('-a22',WithA22,A22)
    call argdbl('-A22',WithA22,A22)

   

    ! ... Forcing specification flags
    ! ...
    if (len_trim(OUfilename).gt.0) WithOU = .True.
    if (len_trim(OVfilename).gt.0) WithOV = .True.
    if (len_trim(OWfilename).gt.0) WithOW = .True.
    if (len_trim(OTfilename).gt.0) WithOT = .True.
    if (len_trim(OSfilename).gt.0) WithOS = .True.
    if (len_trim(ORfilename).gt.0) WithOR = .True.
    if (len_trim(OCfilename).gt.0) WithOW = .True.
    if (len_trim(AUfilename).gt.0) WithAU = .True.
    if (len_trim(AVfilename).gt.0) WithAV = .True.

!    ! ... Analytical density function
!    ! ...
!    call argflg('-rho_ana',rho_ana)
!    call argflg('-Rho_ana',rho_ana)
!    call argflg('-Rho_Ana',rho_ana)


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
