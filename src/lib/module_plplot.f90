! ======================================================================== !
! ALOGES PROJECT                                                           !
! Quim Ballabrera, April 2022                                              !
! Institut de Ciencies del Mar, CSIC                                       !
! Last Modified: 2025-08-28                                                !
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
!                                                                          !
! To correctly compile this module, the user needs to have installed       !
! the PLPLOT library. In ubuntu, this can be done by executing the command !
! > install sudo apt update                                                !
! > sudo apt install plplot-doc plplot-dev                                 !
!                                                                          !
! Functions defined here:                                                  !
! - PLT%init                                                               !
! - PLT%xlim                                                               !
! - PLT%ylim                                                               !
! - PLT%labels                                                             !
! - PLT%plot                                                               !
! - PLT%show                                                               !
! - AX%new                                                                 !
! - AX%xlim                                                                !
! - AX%ylim                                                                !
! - AX%labels                                                              !
! -------------------------------------------------------------------------!

module module_plplot

use module_types
use plplot

implicit none
private
public type_plplot
public init_custom_colormap0

! ... Custom color indices
! ...
integer, public, parameter :: &
                  COLOR_WHITE      = 0, &
                  COLOR_BLACK      = 1, &
                  COLOR_BLUE       = 2, &
                  COLOR_GREEN      = 3, &
                  COLOR_RED        = 4, &
                  COLOR_CYAN       = 5, &
                  COLOR_MAGENTA    = 6, &
                  COLOR_YELLOW     = 7, &
                  COLOR_ORANGE     = 8, &
                  COLOR_PURPLE     = 9, &
                  COLOR_PINK       = 10, &
                  COLOR_BROWN      = 11, &
                  COLOR_GRAY       = 12, &
                  COLOR_LIGHTBLUE  = 13, &
                  COLOR_LIGHTGREEN = 14, &
                  COLOR_LIGHTRED   = 15

integer, private, parameter                :: nmax = 25

type type_legend
  integer                                  :: nlegend = 0
  character(len=20), dimension(nmax)       :: text_labels = ''
  character(len=20), dimension(nmax)       :: symbols = ''
  integer, dimension(nmax)                 :: opt_array
  integer, dimension(nmax)                 :: text_colors
  integer, dimension(nmax)                 :: line_colors
  integer, dimension(nmax)                 :: box_colors
  integer, dimension(nmax)                 :: box_patterns
  integer, dimension(nmax)                 :: symbol_colors
  integer, dimension(nmax)                 :: line_styles
  integer, dimension(nmax)                 :: symbol_numbers
  real(dp), dimension(nmax)                :: line_widths
  real(dp), dimension(nmax)                :: box_line_widths
  real(dp), dimension(nmax)                :: box_scales
  real(dp), dimension(nmax)                :: symbol_scales
  integer                                  :: opt = 0
  integer                                  :: position = 0
  integer                                  :: bg_color = COLOR_WHITE
  integer                                  :: bb_color = COLOR_BLACK
  integer                                  :: bb_style = 1
  integer                                  :: nrow = 0
  integer                                  :: ncolumn = 0
  real(dp)                                 :: legend_width
  real(dp)                                 :: legend_height
  real(dp)                                 :: xoffset = 0.0D0
  real(dp)                                 :: yoffset = 0.0D0
  real(dp)                                 :: plot_width = 0.1D0
  real(dp)                                 :: text_offset = 1.0D0
  real(dp)                                 :: text_scale = 1.0D0
  real(dp)                                 :: text_spacing = 2.0D0
  real(dp)                                 :: text_justification = 0.D0
end type type_legend

type type_plaxis
  logical                                  :: new_axis = .True.
  logical                                  :: env_done = .False.
  logical                                  :: xlimits_set = .False.
  logical                                  :: ylimits_set = .False.
  logical                                  :: show_grid = .False.
  logical                                  :: axes_equal = .False.
  logical                                  :: xlog = .False.
  logical                                  :: ylog = .False.
  integer                                  :: color_index = 0
  real(dp)                                 :: xmin,xmax,ymin,ymax
  real(dp)                                 :: fontsize = 1.0D0
  character(len=50)                        :: xlabel = ''
  character(len=50)                        :: ylabel = ''
  character(len=50)                        :: title = ''
  character(len=50)                        :: linecolor  = 'black'
  character(len=50)                        :: fontcolor  = 'black'
  character(len=50)                        :: foreground = 'black'
  type(type_legend)                        :: axis_legend
  contains
    procedure                   :: new           => plaxis_new
    procedure                   :: ylim          => plaxis_ylim
    procedure                   :: xlim          => plaxis_xlim
    procedure                   :: labels        => plaxis_labels
    procedure                   :: legend        => plaxis_legend
end type type_plaxis


type type_plplot
  type(type_plaxis)                        :: ax
  logical                                  :: multiplot = .False.
  logical                                  :: init_done = .False.
  logical                                  :: env_done = .False.
  integer                                  :: nrows = 1
  integer                                  :: ncols = 1
  real(dp)                                 :: linewidth = 1.0D0
  character(len=20)                        :: device
  character(len=80)                        :: version
  character(len=maxlen)                    :: output_filename
  contains
    procedure                   :: init          => plplot_init
    procedure                   :: show          => plplot_show
    procedure                   :: ylim          => plplot_ylim
    procedure                   :: xlim          => plplot_xlim
    procedure                   :: labels        => plplot_labels
    procedure                   :: plot_y
    procedure                   :: plot_xy
    procedure                   :: plot_xyy
    generic, public             :: plot => plot_y, plot_xy, plot_xyy

end type type_plplot

! ... The name of the XWINDOW may change from version to version.
! ...
character(len=20), parameter               :: XWINDOW = 'qtwidget'

contains
  ! ...
  ! ===================================================================
  ! ===================================================================
  ! ...
  subroutine init_custom_colormap0()
    ! ... Redefine the default colormap0 with more intuitive colors
        
    ! ... Background and foreground colors
    call plscol0(COLOR_WHITE, 255, 255, 255)   ! White
    call plscol0(COLOR_BLACK, 0, 0, 0)         ! Black
        
    ! ... Basic colors
    call plscol0(COLOR_BLUE, 0, 0, 255)        ! Blue
    call plscol0(COLOR_GREEN, 0, 255, 0)       ! Green
    call plscol0(COLOR_RED, 255, 0, 0)         ! Red
    call plscol0(COLOR_CYAN, 0, 255, 255)      ! Cyan
    call plscol0(COLOR_MAGENTA, 255, 0, 255)   ! Magenta
    call plscol0(COLOR_YELLOW, 255, 255, 0)    ! Yellow
        
    ! ... Additional colors
    call plscol0(COLOR_ORANGE, 255, 165, 0)    ! Orange
    call plscol0(COLOR_PURPLE, 128, 0, 128)    ! Purple
    call plscol0(COLOR_PINK, 255, 192, 203)    ! Pink
    call plscol0(COLOR_BROWN, 165, 42, 42)     ! Brown
    call plscol0(COLOR_GRAY, 128, 128, 128)    ! Gray
    
    ! ... Light colors
    call plscol0(COLOR_LIGHTBLUE, 173, 216, 230)   ! LightBlue
    call plscol0(COLOR_LIGHTGREEN, 144, 238, 144)  ! LightGreen
    call plscol0(COLOR_LIGHTRED, 255, 182, 193)    ! LightRed
        
  end subroutine init_custom_colormap0
    
  subroutine print_color_table()
    ! ... Print a table of available colors
    integer :: i
    character(len=20) :: color_names(0:15)
        
    color_names(0) = "WHITE"
    color_names(1) = "BLACK"
    color_names(2) = "BLUE"
    color_names(3) = "GREEN"
    color_names(4) = "RED"
    color_names(5) = "CYAN"
    color_names(6) = "MAGENTA"
    color_names(7) = "YELLOW"
    color_names(8) = "ORANGE"
    color_names(9) = "PURPLE"
    color_names(10) = "PINK"
    color_names(11) = "BROWN"
    color_names(12) = "GRAY"
    color_names(13) = "LIGHTBLUE"
    color_names(14) = "LIGHTGREEN"
    color_names(15) = "LIGHTRED"
    
    write(*, *) "Custom Color Map:"
    write(*, *) "Index  Name"
    write(*, *) "-----  --------------------"
    do i = 0, 15
      write(*, '(I5, 4X, A)') i, trim(color_names(i))
    end do

  end subroutine print_color_table
  ! ...
  ! ===================================================================
  ! ...
  subroutine plplot_init(PLT,device,filename,background,foreground, &
                         multiplot)
  ! ... Subroutine to initialize the PLPLOT graphical routines.
  ! ... Initialization requires setting the device. Optionally, we can
  ! ... specify output plot filename, and the background and foreground
  ! ... colors.

    class(type_plplot), intent(inout)              :: PLT
    character(len=*), intent(in), optional         :: device
    character(len=*), intent(in), optional         :: filename
    character(len=*), intent(in), optional         :: background
    character(len=*), intent(in), optional         :: foreground
    integer, dimension(2), intent(in), optional    :: multiplot

    ! ... Local variables
    ! ...
    integer RGB(3)

    ! ... Retrieve the PLPLOT version:
    ! ...
    call plgver(PLT%version)

    ! ... Check if version has been specified by the user
    ! ...
    if (present(device)) then
      PLT%device = trim(device)
    else
      PLT%device = XWINDOW
    endif

    if (present(filename)) PLT%output_filename = trim(filename)

    if (PLT%device.eq.XWINDOW) then
      call plsdev(PLT%device)
    else
      call plsdev(PLT%device)
      call plsfnam(PLT%output_filename) 
    endif

    call init_custom_colormap0()

    if (present(background)) then
      RGB = color2rgb(background)
      !call plscol0(0, RGB(1), RGB(2), RGB(3))  ! Background
      call plscolbg(RGB(1), RGB(2), RGB(3))  ! Background
    else
      ! ...            R    G    B
      !call plscol0(0, 255, 255, 255)  ! Background: white
      call plscolbg(255, 255, 255)  ! Background
    endif

    if (present(foreground)) then
      RGB = color2rgb(foreground)
      call plscol0(1, RGB(1), RGB(2), RGB(3))  ! Foreground
    else
      ! ...            R    G    B
      call plscol0(1,   0,   0,   0)        ! Foreground: black, axis, text
    endif

    if (present(multiplot)) then
      PLT%multiplot = .True.
      PLT%nrows = multiplot(1)
      PLT%ncols = multiplot(2)
      call plstar(PLT%ncols,PLT%nrows)
    else
      PLT%multiplot = .False.
      call plinit()
    endif

    PLT%init_done = .True.
    PLT%env_done  = .False.
    PLT%ax%color_index = 0
    PLT%ax%axis_legend%nlegend = 0
    PLT%ax%axis_legend%text_labels(:) = ''

  end subroutine plplot_init
  ! ...
  ! ===================================================================
  ! ...
  subroutine plaxis_new(AX,color,grid,xlog,ylog,xlim,ylim)

    class(type_plaxis), intent(inout)            :: AX
    character(len=*), intent(in), optional       :: color
    logical, intent(in), optional                :: grid
    logical, intent(in), optional                :: xlog
    logical, intent(in), optional                :: ylog
    real(dp), dimension(2), intent(in), optional :: xlim
    real(dp), dimension(2), intent(in), optional :: ylim
    
    ! ... Default colors:
    ! ...
    AX%linecolor  = 'black'
    AX%fontcolor  = 'black'
    AX%foreground = 'black'

    if (present(grid))  AX%show_grid  = grid
    if (present(xlog))  AX%xlog       = xlog
    if (present(ylog))  AX%ylog       = ylog
    if (present(ylog))  AX%ylog       = ylog
    if (present(color)) AX%foreground = trim(color)
    if (present(xlim)) then
      AX%xmin = xlim(1)
      AX%xmax = xlim(2)
      AX%xlimits_set = .True.
    endif
    if (present(ylim)) then
      AX%ymin = ylim(1)
      AX%ymax = ylim(2)
      AX%ylimits_set = .True.
    endif
    
    AX%new_axis = .True.
    AX%env_done = .False.
    AX%color_index = 0
    AX%axis_legend%nlegend = 0
    AX%axis_legend%text_labels(:) = ''
 
  end subroutine plaxis_new
  ! ...
  ! ===================================================================
  ! ...
  function color2rgb(name) result(RGB)

    character(len=*), intent(in)              :: name
    integer, dimension(3)                     :: RGB

    select case(name)
    case('w','white','WHITE','White')
      RGB = [255, 255, 255 ]
    case('k','black','BLACK','Black')
      RGB = [0, 0, 0 ]
    case('r','red','RED','Red')
      RGB = [255, 0, 0 ]
    case('g','green','GREEN','Green')
      RGB = [0, 255, 0 ]
    case('b','blue','BLUE','Blue')
      RGB = [0, 0, 255 ]
    case default
      stop 'invalid color name'
    end select 

  end function color2rgb
  ! ...
  ! ===================================================================
  ! ...
  function get_color_index(name) result(icolor)

    character(len=*), intent(in)                :: name
    integer                                     :: icolor

    select case(name)
      case('w','white','WHITE','White'); icolor = COLOR_WHITE
      case('k','black','BLACK','Black'); icolor = COLOR_BLACK
      case('b','blue','BLUE','Blue'); icolor = COLOR_BLUE
      case('g','green','GREEN','Green'); icolor = COLOR_GREEN
      case('r','red','RED','Red'); icolor = COLOR_RED
      case('c','cyan','CYAN','Cyan'); icolor = COLOR_CYAN
      case('m','magenta','MAGENTA','Magenta'); icolor = COLOR_MAGENTA
      case('y','yellow','YELLOW','Yellow'); icolor = COLOR_YELLOW
      case('o','orange','ORANGE','Orange'); icolor = COLOR_ORANGE
      case('purple','PURPLE','Purple'); icolor = COLOR_PURPLE
      case('pink','PINK','Pink'); icolor = COLOR_PINK
      case('brown','BROWN','Brown'); icolor = COLOR_BROWN
      case('gray','GRAY','Gray'); icolor = COLOR_GRAY
      case('lightblue','LIGHTBLUE','Lightblue','LightBlue'); icolor = COLOR_LIGHTBLUE
      case('lightgreen','LIGHTGREEN','Lightgreen','LightGreen'); icolor = COLOR_LIGHTGREEN
      case('lightred','LIGHTRED','Lightred','LightRed'); icolor = COLOR_LIGHTRED
      case default; stop 'Unknown color name'
    end select

  end function get_color_index
  ! ...
  ! ===================================================================
  ! ...
  integer pure function get_marker_code(c)

    character(len=1), intent(in) :: c

    select case(c)
       case ('s'); get_marker_code = 0
       case ('.'); get_marker_code = 1
       case ('+'); get_marker_code = 2
       case ('*'); get_marker_code = 3
       case ('o'); get_marker_code = 4
       case ('x'); get_marker_code = 5
       case default; get_marker_code = 2
    end select

  end function get_marker_code
  ! ...
  ! ===================================================================
  ! ...
  subroutine plplot_show(PLT)

    class(type_plplot), intent(inout)              :: PLT

    ! Close PLplot library
    call plend()

    PLT%init_done = .False.
    PLT%ax%env_done = .False.
    PLT%ax%new_axis = .True.

  end subroutine plplot_show
  ! ...
  ! ===================================================================
  ! ...
  subroutine plplot_xlim(PLT,xmin,xmax)

    class(type_plplot), intent(inout)              :: PLT
    real(dp), intent(in)                           :: xmin,xmax

    PLT%ax%xmin = xmin
    PLT%ax%xmax = xmax
    PLT%ax%xlimits_set = .True.

  end subroutine plplot_xlim
  ! ...
  ! ===================================================================
  ! ...
  subroutine plplot_ylim(PLT,ymin,ymax)

    class(type_plplot), intent(inout)              :: PLT
    real(dp), intent(in)                           :: ymin,ymax

    PLT%ax%ymin = ymin
    PLT%ax%ymax = ymax
    PLT%ax%ylimits_set = .True.

  end subroutine plplot_ylim
  ! ...
  ! ===================================================================
  ! ...
  subroutine plaxis_xlim(AX,xmin,xmax)

    class(type_plaxis), intent(inout)              :: AX
    real(dp), intent(in)                           :: xmin,xmax

    AX%xmin = xmin
    AX%xmax = xmax
    AX%xlimits_set = .True.

  end subroutine plaxis_xlim
  ! ...
  ! ===================================================================
  ! ...
  subroutine plaxis_ylim(AX,ymin,ymax)

    class(type_plaxis), intent(inout)              :: AX
    real(dp), intent(in)                           :: ymin,ymax

    AX%ymin = ymin
    AX%ymax = ymax
    AX%ylimits_set = .True.

  end subroutine plaxis_ylim
  ! ...
  ! ===================================================================
  ! ...
  subroutine plplot_labels(PLT,xlabel,ylabel,title,color,fontsize)

    class(type_plplot), intent(inout)              :: PLT
    character(len=*), intent(in), optional         :: xlabel
    character(len=*), intent(in), optional         :: ylabel
    character(len=*), intent(in), optional         :: title
    character(len=*), intent(in), optional         :: color
    real(dp), intent(in), optional                 :: fontsize

    if (present(xlabel))   PLT%ax%xlabel = trim(xlabel)
    if (present(ylabel))   PLT%ax%ylabel = trim(ylabel)
    if (present(title))    PLT%ax%title  = trim(title)
    if (present(color))    PLT%ax%fontcolor  = trim(color)
    if (present(fontsize)) PLT%ax%fontsize  = fontsize

  end subroutine plplot_labels
  ! ...
  ! ===================================================================
  ! ...
  subroutine plaxis_labels(AX,xlabel,ylabel,title,color,fontsize)

    class(type_plaxis), intent(inout)              :: AX
    character(len=*), intent(in), optional         :: xlabel
    character(len=*), intent(in), optional         :: ylabel
    character(len=*), intent(in), optional         :: title
    character(len=*), intent(in), optional         :: color
    real(dp), intent(in), optional                 :: fontsize

    if (present(xlabel))   AX%xlabel = trim(xlabel)
    if (present(ylabel))   AX%ylabel = trim(ylabel)
    if (present(title))    AX%title  = trim(title)
    if (present(color))    AX%fontcolor  = trim(color)
    if (present(fontsize)) AX%fontsize  = fontsize

  end subroutine plaxis_labels
  ! ...
  ! ===================================================================
  ! ...
  subroutine plot_y(PLT,y,color,linewidth,linestyle,label)
  ! Case 1: only Y(:)

    class(type_plplot), intent(inout)         :: PLT
    real(dp), dimension(:), intent(in)        :: y
    character(len=*), intent(in), optional    :: color
    real(dp), intent(in), optional            :: linewidth
    character(len=*), intent(in), optional    :: linestyle
    character(len=*), intent(in), optional    :: label

    ! ... Local variables
    ! ...
    integer i
    real(dp) x(size(y))

    do i=1,size(y)
      x(i) = 1.0D0*i
    enddo

    call plot_xy(PLT,x,y,color,linewidth,label)

  end subroutine plot_y
  ! ...
  ! ===================================================================
  ! ...
  subroutine plot_xy(PLT,x,y,color,linewidth,linestyle,label)
  ! Case 2: X(:), Y(:)

    class(type_plplot), intent(inout)         :: PLT
    real(dp), dimension(:), intent(in)        :: x
    real(dp), dimension(:), intent(in)        :: y
    character(len=*), intent(in), optional    :: color
    real(dp), intent(in), optional            :: linewidth
    character(len=*), intent(in), optional    :: linestyle
    character(len=*), intent(in), optional    :: label

    ! ... Local variables
    ! ...
    logical new_environment
    integer ibox,ncode,k,linecolor,RGB(3)
    real(dp) xmin,xmax,ymin,ymax

    new_environment = .False.

    ! ... Check if the plot has been initialized
    ! ...
    if (.not.PLT%init_done) call PLT%init()

    if (PLT%ax%new_axis.and..not.PLT%ax%env_done) then
      new_environment = .True.
      PLT%ax%env_done = .True.
      PLT%ax%new_axis = .False.
    endif
  
    if (new_environment) then
      if (.not.PLT%ax%xlimits_set) then
        xmin = minval(x)
        xmax = maxval(x)
      else
        xmin = PLT%ax%xmin
        xmax = PLT%ax%xmax
      endif

      if (.not.PLT%ax%ylimits_set) then
        ymin = minval(y)
        ymax = maxval(y)
      else
        ymin = PLT%ax%ymin
        ymax = PLT%ax%ymax
      endif

      ! ... Grid
      ! ...
      ibox = 10*merge(1,0,PLT%ax%xlog) + 20*merge(1,0,PLT%ax%ylog) + 2*merge(1,0,PLT%ax%show_grid)

      ! ... Set the new environment
      ! ... The color from the axis foreground
      ! ...
      !RGB = color2rgb(PLT%ax%foreground)
      !call plscol0(10, RGB(1), RGB(2), RGB(3))  
      call plcol0(1)                             ! 1 - Color foreground
      call plenv(xmin,xmax,ymin,ymax,0,ibox)
    endif

    ! ... Labels (black)
    ! ...
    !RGB = color2rgb(PLT%ax%fontcolor)
    !call plscol0(10, RGB(1), RGB(2), RGB(3))  ! Background
    !call plcol0(10)
    call plcol0(get_color_index(PLT%ax%fontcolor))
    call plschr(0.0D0,PLT%ax%fontsize) 
    call pllab( trim(PLT%ax%xlabel), trim(PLT%ax%ylabel), trim(PLT%ax%title) )

    if (present(color)) then
      !RGB = color2rgb(color)
      !call plscol0(10, RGB(1), RGB(2), RGB(3))  ! Background
      !call plcol0(10)
      linecolor = get_color_index(color)
    else
      PLT%ax%color_index = PLT%ax%color_index + 1
      if (PLT%ax%color_index.GT.15) PLT%ax%color_index = 1
      linecolor = PLT%ax%color_index
    endif
    PLT%ax%linecolor = trim(color)
    call plcol0(linecolor)

    ! ... Line width
    if (present(linewidth)) then
      call plwidth(linewidth)
    else
      call plwidth(PLT%linewidth)
    endif
    
    ! ... Line style / markers
    if (present(linestyle)) then
       select case(trim(linestyle))
       CASE ('-')
          CALL plline(x, y)
       CASE ('--',':','-.')
          !CALL plot_dashed(x, y, linestyle)
          stop 'Dashed still not coded'
       case default
          ncode = get_marker_code(linestyle(1:1))
          CALL plpoin(x, y, ncode)
       end select
    else
       call plline(x, y)
    endif

    if (present(label)) then
      PLT%ax%axis_legend%nlegend = PLT%ax%axis_legend%nlegend + 1
      k = PLT%ax%axis_legend%nlegend
      PLT%ax%axis_legend%opt_array(k) = PL_LEGEND_LINE
      PLT%ax%axis_legend%text_colors(k) = COLOR_BLACK
      PLT%ax%axis_legend%text_labels(k) = trim(label)
      PLT%ax%axis_legend%line_colors(k) = linecolor
      PLT%ax%axis_legend%line_styles(k) = 1
      PLT%ax%axis_legend%line_widths(k) = PLT%linewidth
    endif
    
  end subroutine plot_xy
  ! ...
  ! ===================================================================
  ! ...
  subroutine plot_xyy(PLT, x, y)
  ! Case 2: X(:), Y(:,:)

    class(type_plplot), intent(inout)         :: PLT
    real(dp), dimension(:), intent(in)        :: x
    real(dp), dimension(:,:), intent(in)      :: y

    ! your plotting logic here

  end subroutine plot_xyy
  ! ...
  ! ===================================================================
  ! ...
  subroutine plaxis_legend(AX)

    class(type_plaxis), intent(inout)         :: AX

    call plt_legend_show(AX%axis_legend)

  end subroutine plaxis_legend
  ! ...
  ! ===================================================================
  ! ...
  subroutine plt_legend_show(L)
  ! 
  ! ... Subroutine to add a legend in the plot.
  ! ... The properties of each legend entry are defined when the ploty,
  ! ... plotxy are called, it a "label" argument is passed.
  ! ...
  !  !   First legend entry.
  !  opt_array(1)   = PL_LEGEND_LINE
  !  text_colors(1) = COLOR_BLACK
  !  text_labels(1) = 'sin(x)'
  !  line_colors(1) = COLOR_RED
  !  line_styles(1) = 1
  !  line_widths(1) = 1
  !  !   note from the above opt_array the first symbol (and box) indices
  !  !   do not have to be specified, at least in C. For Fortran we need
  !  !   to set the symbols to be something, since the string is always
  !  !   copied as part of the bindings.
  !  symbols(1) = ''
  !
  !  !   Second legend entry.
  !  opt_array(2)      = PL_LEGEND_LINE
  !  text_colors(2)    = COLOR_BLACK
  !  text_labels(2)    = 'my cos(x)'
  !  line_colors(2)    = COLOR_BLUE
  !  line_styles(2)    = 1
  !  line_widths(2)    = 1
  !
  !  ! opt: contains bits controlling the overall legend. If the PL_LEGEND_TEXT_LEFT bit 
  !  ! is set, put the text area on the left of the legend and the plotted area on the 
  !  ! right. Otherwise, put the text area on the right of the legend and the plotted 
  !  ! area on the left. If the PL_LEGEND_BACKGROUND bit is set, plot a (semitransparent) 
  !  ! background for the legend. If the PL_LEGEND_BOUNDING_BOX bit is set, plot a bounding 
  !  ! box for the legend. If the PL_LEGEND_ROW_MAJOR bit is set and (both of the possibly 
  !  ! internally transformed) nrow > 1 and ncolumn > 1, then plot the resulting array of 
  !  ! legend entries in row-major order. Otherwise, plot the legend entries in column-major order.
  !  !
  !  opt = PL_LEGEND_BACKGROUND 
  !  opt = PL_LEGEND_BACKGROUND + PL_LEGEND_BOUNDING_BOX
  !
  !  ! position: contains bits which control the overall position of the legend and the 
  !  ! definition of the adopted coordinates used for positions just like what is done 
  !  ! for the position argument for plcolorbar. However, note that the defaults for the 
  !  ! position bits (see below) are different than the plcolorbar case. The combination 
  !  ! of the PL_POSITION_LEFT, PL_POSITION_RIGHT, PL_POSITION_TOP, PL_POSITION_BOTTOM, 
  !  ! PL_POSITION_INSIDE, and PL_POSITION_OUTSIDE bits specifies one of the 16 possible 
  !  ! standard positions (the 4 corners and centers of the 4 sides for both the inside 
  !  ! and outside cases) of the legend relative to the adopted coordinate system. The 
  !  ! corner positions are specified by the appropriate combination of two of the 
  !  ! PL_POSITION_LEFT, PL_POSITION_RIGHT, PL_POSITION_TOP, and PL_POSITION_BOTTOM bits 
  !  ! while the sides are specified by a single value of one of those bits. The adopted 
  !  ! coordinates are normalized viewport coordinates if the PL_POSITION_VIEWPORT bit is 
  !  ! set or normalized subpage coordinates if the PL_POSITION_SUBPAGE bit is set. 
  !  ! Default position bits: If none of PL_POSITION_LEFT, PL_POSITION_RIGHT, 
  !  ! PL_POSITION_TOP, or PL_POSITION_BOTTOM are set, then use the combination of 
  !  ! PL_POSITION_RIGHT and PL_POSITION_TOP. If neither of PL_POSITION_INSIDE or 
  !  ! PL_POSITION_OUTSIDE is set, use PL_POSITION_INSIDE. If neither of PL_POSITION_VIEWPORT 
  !  ! or PL_POSITION_SUBPAGE is set, use PL_POSITION_VIEWPORT. 
  !  ! 
  !  position = 0
  !  position = PL_POSITION_RIGHT + PL_POSITION_TOP
  !
  !  xoffset = 0.0D0   ! Negative to the right, positive to the left
  !  yoffset = 0.0D0   ! Negative upward, positive downward
  !
  !  plot_width = 0.06D0
  !  bg_color   = COLOR_WHITE   ! legend background's color
  !  bb_color   = COLOR_GREEN   ! legend box's color
  !  bb_style   = 2             ! legend box's line style: 1-8
  !                             ! Integer value between 1 and 8. Line style 1 is a continuous line, 
  !                             ! line style 2 is a line with short dashes and gaps, line style 3 
  !                             ! is a line with long dashes and gaps, line style 4 has long dashes 
  !                             ! and short gaps and so on. 
  !  nrow       = 0
  !  ncolumn    = 0
  !
  !  ! opt_array: A vector of nlegend values of options to control each individual plotted area 
  !  ! corresponding to a legend entry. If the PL_LEGEND_NONE bit is set, then nothing is plotted 
  !  ! in the plotted area. If the PL_LEGEND_COLOR_BOX, PL_LEGEND_LINE, and/or PL_LEGEND_SYMBOL 
  !  ! bits are set, the area corresponding to a legend entry is plotted with a colored box; a 
  !  ! line; and/or a line of symbols. 
  !
  !  ! text_offset: Offset of the text area from the plot area in units of character width. 
  !  ! text_scale: Character height scale for text annotations. 
  !  ! text_spacing: Vertical spacing in units of the character height from one legend entry to the next. 
  !  ! text_justification: Justification parameter used for text justification. The most common values of 
  !  !                     text_justification are 0., 0.5, or 1. corresponding to a text that is left 
  !  !                     justified, centred, or right justified within the text area, but other values 
  !  !                     are allowed as well. 
  !  text_offset  = 1.0D0
  !  text_scale   = 1.0D0
  !  text_spacing = 2.0D0
  !  text_justification = 0.0D0
  !
  !  call pllegend( legend_width, legend_height, &
  !             opt, position, &
  !             xoffset, yoffset, plot_width, bg_color, &
  !             bb_color, bb_style, &
  !             nrow, ncolumn, &
  !             opt_array(1:nlegend), &
  !             text_offset,text_scale,text_spacing,text_justification, &
  !             text_colors(1:nlegend), text(1:nlegend), &
  !             box_colors(1:nlegend), box_patterns(1:nlegend), box_scales(1:nlegend), box_line_widths(1:nlegend), &
  !             line_colors(1:nlegend), line_styles(1:nlegend), line_widths(1:nlegend) , &
  !             symbol_colors(1:nlegend), symbol_scales(1:nlegend), symbol_numbers(1:nlegend), symbols(1:nlegend) )
  !
  !
    class(type_legend), intent(inout)         :: L

    ! ... Local variables
    ! ...
    integer nlegend

    nlegend = L%nlegend

    if (nlegend.LE.0) return

    call pllegend( L%legend_width, L%legend_height, &
                   L%opt,                               &
                   L%position, &
                   L%xoffset, L%yoffset, L%plot_width,  &
                   L%bg_color, L%bb_color, L%bb_style, &
                   L%nrow, L%ncolumn, &
                   L%opt_array(1:nlegend), &
                   L%text_offset, L%text_scale, L%text_spacing, L%text_justification, &
                   L%text_colors(1:nlegend), L%text_labels(1:nlegend), &
                   L%box_colors(1:nlegend), L%box_patterns(1:nlegend), L%box_scales(1:nlegend), L%box_line_widths(1:nlegend), &
                   L%line_colors(1:nlegend), L%line_styles(1:nlegend), L%line_widths(1:nlegend) , &
                   L%symbol_colors(1:nlegend), L%symbol_scales(1:nlegend), L%symbol_numbers(1:nlegend), L%symbols(1:nlegend) )
  
  end subroutine plt_legend_show
  ! ...
  ! ===================================================================
  ! ...

end module module_plplot

