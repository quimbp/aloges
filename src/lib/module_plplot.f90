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
! - plot                                                                   !
! - quiver                                                                 !
! - print_color_table                                                      !
! -------------------------------------------------------------------------!

module module_plplot

use module_types
use plplot

implicit none
private
public type_plplot
public init_custom_colormap0
public line_options, line_opts, quiver_options, quiver_opts
public PIVOT_TAIL, PIVOT_MID, PIVOT_TIP

! Custom color indices
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

character(len=10), dimension(16)  :: &
   color_names = [ "WHITE     ", "BLACK     ", "BLUE      ", "GREEN     ", &
                   "RED       ", "CYAN      ", "MAGENTA   ", "YELLOW    ", &
                   "ORANGE    ", "PURPLE    ", "PINK      ", "BROWN     ", &
                   "GRAY      ", "LIGHTBLUE ", "LIGHTGREEN", "LIGHTRED  " ]
    
! Pivot options for quiver plots
integer, parameter :: &
    PIVOT_TAIL = 1, &
    PIVOT_MID  = 2, &
    PIVOT_TIP  = 3

integer, private, parameter :: nmax = 25

! Line plot options type
type line_options
  real(dp)              :: width = 1.0_dp
  character(len=10)     :: color = 'black'
  character(len=20)     :: label = ''
  character(len=2)      :: style = '-'
end type line_options

! Quiver options type
type quiver_options
    real(dp)             :: scale = 1.0_dp
    real(dp)             :: head_size = 0.02_dp
    character(len=10)    :: color = 'black'
    logical              :: normalize = .true.
    real(dp)             :: min_magnitude = 0.001_dp
    integer              :: pivot = PIVOT_TAIL
end type quiver_options

type type_legend
  integer                              :: nlegend = 0
  character(len=20), dimension(nmax)   :: text_labels = ''
  character(len=20), dimension(nmax)   :: symbols = ''
  integer, dimension(nmax)             :: opt_array
  integer, dimension(nmax)             :: text_colors
  integer, dimension(nmax)             :: line_colors
  integer, dimension(nmax)             :: box_colors
  integer, dimension(nmax)             :: box_patterns
  integer, dimension(nmax)             :: symbol_colors
  integer, dimension(nmax)             :: line_styles
  integer, dimension(nmax)             :: symbol_numbers
  real(dp), dimension(nmax)            :: line_widths
  real(dp), dimension(nmax)            :: box_line_widths
  real(dp), dimension(nmax)            :: box_scales
  real(dp), dimension(nmax)            :: symbol_scales
  integer                              :: opt = 0
  integer                              :: position = 0
  integer                              :: bg_color = COLOR_WHITE
  integer                              :: bb_color = COLOR_BLACK
  integer                              :: bb_style = 1
  integer                              :: nrow = 0
  integer                              :: ncolumn = 0
  real(dp)                             :: legend_width
  real(dp)                             :: legend_height
  real(dp)                             :: xoffset = 0.0_dp
  real(dp)                             :: yoffset = 0.0_dp
  real(dp)                             :: plot_width = 0.1_dp
  real(dp)                             :: text_offset = 1.0_dp
  real(dp)                             :: text_scale = 1.0_dp
  real(dp)                             :: text_spacing = 2.0_dp
  real(dp)                             :: text_justification = 0.0_dp
end type type_legend

type type_plaxis
  logical                              :: new_axis = .true.
  logical                              :: env_done = .false.
  logical                              :: xlimits_set = .false.
  logical                              :: ylimits_set = .false.
  logical                              :: show_grid = .false.
  logical                              :: axes_equal = .false.
  logical                              :: xlog = .false.
  logical                              :: ylog = .false.
  integer                              :: color_index = 3
  real(dp)                             :: xmin, xmax, ymin, ymax
  real(dp)                             :: fontsize = 1.0_dp
  character(len=50)                    :: xlabel = ''
  character(len=50)                    :: ylabel = ''
  character(len=50)                    :: title = ''
  character(len=50)                    :: fontcolor = 'black'
  character(len=50)                    :: foreground = 'black'
  type(type_legend)                    :: axis_legend
  contains
    procedure :: new    => plaxis_new
    procedure :: ylim   => plaxis_ylim
    procedure :: xlim   => plaxis_xlim
    procedure :: labels => plaxis_labels
    procedure :: legend => plaxis_legend
end type type_plaxis

type type_plplot
  type(type_plaxis)                    :: ax
  logical                              :: multiplot = .false.
  logical                              :: init_done = .false.
  integer                              :: nrows = 1
  integer                              :: ncols = 1
  character(len=20)                    :: device
  character(len=80)                    :: version
  character(len=maxlen)                :: output_filename
  contains
    procedure :: init   => plplot_init
    procedure :: show   => plplot_show
    procedure :: ylim   => plplot_ylim
    procedure :: xlim   => plplot_xlim
    procedure :: labels => plplot_labels
    procedure :: legend => plplot_legend
    procedure :: plot_y
    procedure :: plot_xy
    procedure :: plot_yy
    procedure :: plot_xyy
    generic, public :: plot => plot_y, plot_xy, plot_yy, plot_xyy
    procedure :: quiver_2d
    procedure :: quiver_1d
    generic, public :: quiver => quiver_2d, quiver_1d
end type type_plplot

character(len=20), parameter :: XWINDOW = 'qtwidget'

contains

  ! ==== Initialization and utility functions ====
  
  subroutine init_custom_colormap0()
    call plscol0(COLOR_WHITE, 255, 255, 255)
    call plscol0(COLOR_BLACK, 0, 0, 0)
    call plscol0(COLOR_BLUE, 0, 0, 255)
    call plscol0(COLOR_GREEN, 0, 255, 0)
    call plscol0(COLOR_RED, 255, 0, 0)
    call plscol0(COLOR_CYAN, 0, 255, 255)
    call plscol0(COLOR_MAGENTA, 255, 0, 255)
    call plscol0(COLOR_YELLOW, 255, 255, 0)
    call plscol0(COLOR_ORANGE, 255, 165, 0)
    call plscol0(COLOR_PURPLE, 128, 0, 128)
    call plscol0(COLOR_PINK, 255, 192, 203)
    call plscol0(COLOR_BROWN, 165, 42, 42)
    call plscol0(COLOR_GRAY, 128, 128, 128)
    call plscol0(COLOR_LIGHTBLUE, 173, 216, 230)
    call plscol0(COLOR_LIGHTGREEN, 144, 238, 144)
    call plscol0(COLOR_LIGHTRED, 255, 182, 193)
  end subroutine init_custom_colormap0
  ! ...
  ! ===================================================================
  ! ...
  function color2rgb(name) result(RGB)
    character(len=*), intent(in) :: name 
    integer, dimension(3) :: RGB  

    select case(trim(name))
    case('w','white','WHITE','White'); RGB = [255, 255, 255]
    case('k','black','BLACK','Black'); RGB = [0, 0, 0]
    case('b','blue','BLUE','Blue'); RGB = [0, 0, 255]
    case('g','green','GREEN','Green'); RGB = [0, 255, 0]
    case('r','red','RED','Red'); RGB = [255, 0, 0]
    case('c','cyan','Cyan','CYAN'); RGB = [0, 255, 255]
    case('m','magenta','Magenta','MAGENTA'); RGB = [255, 0, 255]
    case('y','yellow','Yellow','YELLOW'); RGB = [255, 255, 0]
    case('o','orange','Orange','ORANGE'); RGB = [255, 165, 0]
    case('purple','Purple','PURPLE'); RGB = [128, 0, 128]
    case('pink','Pink','PINK'); RGB = [255, 192, 203]
    case('brown','Brown','BROWN'); RGB = [165, 42, 42]
    case('gray','Gray','GRAY','grey','Grey','GREY'); RGB = [128, 128, 128]
    case default 
      write(*,*) 'ERROR: Invalid color name: ', trim(name)
      stop 'Invalid color name'
    end select  
  end function color2rgb
  ! ...
  ! ===================================================================
  ! ...
  function get_color_index(name) result(icolor)
    character(len=*), intent(in) :: name
    integer :: icolor

    select case(trim(name))
    case('w','white','WHITE','White')
      icolor = COLOR_WHITE
    case('k','black','BLACK','Black')
      icolor = COLOR_BLACK
    case('b','blue','BLUE','Blue')
      icolor = COLOR_BLUE
    case('g','green','GREEN','Green')
      icolor = COLOR_GREEN
    case('r','red','RED','Red')
      icolor = COLOR_RED
    case('c','cyan','CYAN','Cyan')
      icolor = COLOR_CYAN
    case('m','magenta','MAGENTA','Magenta')
      icolor = COLOR_MAGENTA
    case('y','yellow','YELLOW','Yellow')
      icolor = COLOR_YELLOW
    case('o','orange','ORANGE','Orange')
      icolor = COLOR_ORANGE
    case('purple','PURPLE','Purple')
      icolor = COLOR_PURPLE
    case('pink','PINK','Pink')
      icolor = COLOR_PINK
    case('brown','BROWN','Brown')
      icolor = COLOR_BROWN
    case('gray','GRAY','Gray','grey','Grey','GREY')
      icolor = COLOR_GRAY
    case('lightblue','LIGHTBLUE','Lightblue','LightBlue')
      icolor = COLOR_LIGHTBLUE
    case('lightgreen','LIGHTGREEN','Lightgreen','LightGreen')
      icolor = COLOR_LIGHTGREEN
    case('lightred','LIGHTRED','Lightred','LightRed')
      icolor = COLOR_LIGHTRED
    case default
      write(*,*) 'ERROR: Unknown color name: ', trim(name)
      stop 'Unknown color name'
    end select
  end function get_color_index
  ! ...
  ! ===================================================================
  ! ...
  integer pure function get_line_style_code(style)
    character(len=*), intent(in) :: style

    select case(trim(style))
    case ('-')
      get_line_style_code = 1  ! Solid line
    case ('--')
      get_line_style_code = 2  ! Dashed line
    case (':')
      get_line_style_code = 3  ! Dotted line
    case ('-.')
      get_line_style_code = 4  ! Dash-dot line
    case default
      get_line_style_code = 1  ! Default to solid
    end select
  end function get_line_style_code
  ! ...
  ! ===================================================================
  ! ...
  integer pure function get_marker_code(c)
    character(len=1), intent(in) :: c

    select case(c)
    case ('s')
      get_marker_code = 0  ! Square
    case ('.')
      get_marker_code = 1  ! Dot
    case ('+')
      get_marker_code = 2  ! Plus
    case ('*')
      get_marker_code = 3  ! Asterisk
    case ('o')
      get_marker_code = 4  ! Circle
    case ('x')
      get_marker_code = 5  ! Cross
    case default
      get_marker_code = 2  ! Default to plus
    end select
  end function get_marker_code

  subroutine print_color_table()
    integer :: i
    !character(len=20) :: color_names(0:15)
    !color_names = ["WHITE      ", "BLACK      ", "BLUE       ", "GREEN      ", &
    !               "RED        ", "CYAN       ", "MAGENTA    ", "YELLOW     ", &
    !               "ORANGE     ", "PURPLE     ", "PINK       ", "BROWN      ", &
    !               "GRAY       ", "LIGHTBLUE  ", "LIGHTGREEN ", "LIGHTRED   "]

    write(*, *) "ALOGES PLPLT Color Map:"
    write(*, *) "Index  Name"
    write(*, *) "-----  ----"
    do i = 1, 16
      write(*, '(I5, 4X, A)') i-1, trim(color_names(i))
    end do
  end subroutine print_color_table

  subroutine plplot_init(PLT, device, filename, background, foreground, multiplot)
    class(type_plplot), intent(inout) :: PLT
    character(len=*), intent(in), optional :: device
    character(len=*), intent(in), optional :: filename
    character(len=*), intent(in), optional :: background
    character(len=*), intent(in), optional :: foreground
    integer, dimension(2), intent(in), optional :: multiplot
    integer :: RGB(3)

    call plgver(PLT%version)

    PLT%device = XWINDOW
    if (present(device)) PLT%device = trim(device)
    if (present(filename)) PLT%output_filename = trim(filename)

    if (PLT%device == XWINDOW) then
      call plsdev(PLT%device)
    else
      call plsdev(PLT%device)
      call plsfnam(PLT%output_filename) 
    endif

    call init_custom_colormap0()

    if (present(background)) then
      RGB = color2rgb(background)
      call plscolbg(RGB(1), RGB(2), RGB(3))
    else
      call plscolbg(255, 255, 255)
    endif

    if (present(foreground)) then
      RGB = color2rgb(foreground)
      call plscol0(1, RGB(1), RGB(2), RGB(3))
    else
      call plscol0(1, 0, 0, 0)
    endif

    if (present(multiplot)) then
      PLT%multiplot = .true.
      PLT%nrows = multiplot(1)
      PLT%ncols = multiplot(2)
      call plstar(PLT%ncols, PLT%nrows)
    else
      PLT%multiplot = .false.
      call plinit()
    endif

    PLT%init_done = .true.
    PLT%ax%color_index = 3
    PLT%ax%axis_legend%nlegend = 0
    PLT%ax%axis_legend%text_labels(:) = ''
  end subroutine plplot_init

  subroutine plaxis_new(AX, color, grid, xlog, ylog, xlim, ylim)
    class(type_plaxis), intent(inout) :: AX
    character(len=*), intent(in), optional :: color
    logical, intent(in), optional :: grid
    logical, intent(in), optional :: xlog
    logical, intent(in), optional :: ylog
    real(dp), dimension(2), intent(in), optional :: xlim
    real(dp), dimension(2), intent(in), optional :: ylim
    
    AX%fontcolor = 'black'
    AX%foreground = 'black'

    if (present(grid)) AX%show_grid = grid
    if (present(xlog)) AX%xlog = xlog
    if (present(ylog)) AX%ylog = ylog
    if (present(color)) AX%foreground = trim(color)
    
    if (present(xlim)) then
      AX%xmin = xlim(1)
      AX%xmax = xlim(2)
      AX%xlimits_set = .true.
    endif
    
    if (present(ylim)) then
      AX%ymin = ylim(1)
      AX%ymax = ylim(2)
      AX%ylimits_set = .true.
    endif
    
    AX%new_axis = .true.
    AX%env_done = .false.
    AX%color_index = 0
    AX%axis_legend%nlegend = 0
    AX%axis_legend%text_labels(:) = ''
  end subroutine plaxis_new

  subroutine plplot_show(PLT)
    class(type_plplot), intent(inout) :: PLT

    if (.not.PLT%init_done) then
      write(*,*) 'Warning: in plplot_show, plot not initialized'
      return
    endif

    call plend()
    PLT%init_done = .false.
    PLT%ax%env_done = .false.
    PLT%ax%new_axis = .true.
  end subroutine plplot_show

  ! ==== Axis limit functions ====
  
  subroutine plplot_xlim(PLT, xmin, xmax)
    class(type_plplot), intent(inout) :: PLT
    real(dp), intent(in) :: xmin, xmax
    PLT%ax%xmin = xmin
    PLT%ax%xmax = xmax
    PLT%ax%xlimits_set = .true.
  end subroutine plplot_xlim

  subroutine plplot_ylim(PLT, ymin, ymax)
    class(type_plplot), intent(inout) :: PLT
    real(dp), intent(in) :: ymin, ymax
    PLT%ax%ymin = ymin
    PLT%ax%ymax = ymax
    PLT%ax%ylimits_set = .true.
  end subroutine plplot_ylim

  subroutine plaxis_xlim(AX, xmin, xmax)
    class(type_plaxis), intent(inout) :: AX
    real(dp), intent(in) :: xmin, xmax
    AX%xmin = xmin
    AX%xmax = xmax
    AX%xlimits_set = .true.
  end subroutine plaxis_xlim

  subroutine plaxis_ylim(AX, ymin, ymax)
    class(type_plaxis), intent(inout) :: AX
    real(dp), intent(in) :: ymin, ymax
    AX%ymin = ymin
    AX%ymax = ymax
    AX%ylimits_set = .true.
  end subroutine plaxis_ylim

  ! ==== Label functions ====
  
  subroutine plplot_labels(PLT, xlabel, ylabel, title, color, fontsize)
    class(type_plplot), intent(inout) :: PLT
    character(len=*), intent(in), optional :: xlabel
    character(len=*), intent(in), optional :: ylabel
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: color
    real(dp), intent(in), optional :: fontsize

    if (present(xlabel)) PLT%ax%xlabel = trim(xlabel)
    if (present(ylabel)) PLT%ax%ylabel = trim(ylabel)
    if (present(title)) PLT%ax%title = trim(title)
    if (present(color)) PLT%ax%fontcolor = trim(color)
    if (present(fontsize)) PLT%ax%fontsize = fontsize
  end subroutine plplot_labels

  subroutine plaxis_labels(AX, xlabel, ylabel, title, color, fontsize)
    class(type_plaxis), intent(inout) :: AX
    character(len=*), intent(in), optional :: xlabel
    character(len=*), intent(in), optional :: ylabel
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: color
    real(dp), intent(in), optional :: fontsize

    if (present(xlabel)) AX%xlabel = trim(xlabel)
    if (present(ylabel)) AX%ylabel = trim(ylabel)
    if (present(title)) AX%title = trim(title)
    if (present(color)) AX%fontcolor = trim(color)
    if (present(fontsize)) AX%fontsize = fontsize
  end subroutine plaxis_labels

  ! ==== Shared environment setup ====
  
  subroutine setup_plot_environment(PLT, x, y)
    class(type_plplot), intent(inout) :: PLT
    real(dp), dimension(:), intent(in), optional :: x, y
    integer :: ibox
    real(dp) :: xmin, xmax, ymin, ymax

    if (.not. PLT%init_done) call PLT%init()

    if (PLT%ax%new_axis .and. .not. PLT%ax%env_done) then
      if (.not. PLT%ax%xlimits_set) then
        if (present(x)) then
          xmin = minval(x)
          xmax = maxval(x)
        else
          xmin = 0.0_dp
          xmax = 1.0_dp
        endif
      else
        xmin = PLT%ax%xmin
        xmax = PLT%ax%xmax
      endif

      if (.not. PLT%ax%ylimits_set) then
        if (present(y)) then
          ymin = minval(y)
          ymax = maxval(y)
        else
          ymin = 0.0_dp
          ymax = 1.0_dp
        endif
      else
        ymin = PLT%ax%ymin
        ymax = PLT%ax%ymax
      endif

      ibox = 10 * merge(1, 0, PLT%ax%xlog) + &
             20 * merge(1, 0, PLT%ax%ylog) + &
              2 * merge(1, 0, PLT%ax%show_grid)

      call plcol0(1)
      call plenv(xmin, xmax, ymin, ymax, 0, ibox)
      
      call plcol0(get_color_index(PLT%ax%fontcolor))
      call plschr(0.0_dp, PLT%ax%fontsize) 
      call pllab(trim(PLT%ax%xlabel), trim(PLT%ax%ylabel), trim(PLT%ax%title))

      PLT%ax%env_done = .true.
      PLT%ax%new_axis = .false.
    endif
  end subroutine setup_plot_environment

  ! ==== Line plotting functions ====
  
  subroutine plot_y(PLT, y, options)
    class(type_plplot), intent(inout) :: PLT
    real(dp), dimension(:), intent(in) :: y
    type(line_options), intent(in), optional :: options
    integer :: i
    real(dp) :: x(size(y))

    do i = 1, size(y)
      x(i) = real(i, dp)
    enddo
    
    call plot_xy(PLT, x, y, options)
  end subroutine plot_y
  ! ...
  ! ===================================================================
  ! ...
  subroutine plot_yy(PLT, y, options)

    class(type_plplot), intent(inout) :: PLT
    real(dp), dimension(:,:), intent(in) :: y
    type(line_options), intent(in), optional :: options

    ! ... Local variables
    ! ...
    integer :: i
    real(dp), dimension(size(y,1)) :: x
    type(line_options) :: opts

    opts = line_options()  ! Default values
    if (present(options)) opts = options

    forall(i=1:size(y,1)) x(i) = real(i,dp)

    ! Plot each column of y
    ! cycling color
    do i = 1, size(y, 2)
      opts%color = color_names(PLT%ax%color_index)
      call plot_xy(PLT, x, y(:,i), opts)
      PLT%ax%color_index = PLT%ax%color_index + 1
      if (PLT%ax%color_index.eq.17) PLT%ax%color_index = 1
    end do

  end subroutine plot_yy
  ! ...
  ! ===================================================================
  ! ...
  subroutine plot_xyy(PLT, x, y, options)

    class(type_plplot), intent(inout) :: PLT
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(:,:), intent(in) :: y
    type(line_options), intent(in), optional :: options

    ! ... Local variables
    ! ...
    integer :: i
    type(line_options) :: opts

    opts = line_options()  ! Default values
    if (present(options)) opts = options

    ! Plot each column of y
    ! cycling color
    do i = 1, size(y, 2)
      opts%color = color_names(PLT%ax%color_index)
      call plot_xy(PLT, x, y(:,i), opts)
      PLT%ax%color_index = PLT%ax%color_index + 1
      if (PLT%ax%color_index.eq.17) PLT%ax%color_index = 1
    end do

  end subroutine plot_xyy
  ! ...
  ! ===================================================================
  ! ...
  subroutine plot_xy(PLT, x, y, options)

    class(type_plplot), intent(inout) :: PLT
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(:), intent(in) :: y
    type(line_options), intent(in), optional :: options
    integer :: ncode, k, color_index, line_style
    type(line_options) :: opts

    opts = line_options()  ! Default values
    if (present(options)) opts = options
    
    call setup_plot_environment(PLT, x, y)

    color_index = get_color_index(opts%color)
    call plcol0(color_index)
    call plwidth(opts%width)
    
    ! Check if it's a marker style or line style
    if (len_trim(opts%style) == 1 .and. &
        any(opts%style(1:1) == ['s', '.', '+', '*', 'o', 'x'])) then
      ! It's a marker
      ncode = get_marker_code(opts%style(1:1))
      call plpoin(x, y, ncode)
    else
      ! It's a line style
      line_style = get_line_style_code(opts%style)
      call pllsty(line_style)
      call plline(x, y)
      call pllsty(1)  ! Reset to solid line
    endif

    ! Add to legend if label is provided
    if (len_trim(opts%label) > 0) then
      PLT%ax%axis_legend%nlegend = PLT%ax%axis_legend%nlegend + 1
      k = PLT%ax%axis_legend%nlegend
      PLT%ax%axis_legend%opt_array(k) = PL_LEGEND_LINE
      PLT%ax%axis_legend%text_colors(k) = COLOR_BLACK
      PLT%ax%axis_legend%text_labels(k) = trim(opts%label)
      PLT%ax%axis_legend%line_colors(k) = color_index
      PLT%ax%axis_legend%line_styles(k) = line_style
      PLT%ax%axis_legend%line_widths(k) = opts%width
    endif
  end subroutine plot_xy

  ! ==== Quiver plotting functions ====
  
  subroutine quiver_2d(PLT, x, y, u, v, options)
    class(type_plplot), intent(inout) :: PLT
    real(dp), intent(in) :: x(:), y(:), u(:,:), v(:,:)
    type(quiver_options), intent(in), optional :: options
    type(quiver_options) :: opts
    real(dp) :: u_norm(size(u,1), size(u,2))
    real(dp) :: v_norm(size(v,1), size(v,2))
    integer :: nx, ny, i, j, color_index
    real(dp) :: dx, dy
    
    opts = quiver_options()  ! Default values
    if (present(options)) opts = options
    
    nx = size(x)
    ny = size(y)
    
    if (size(u,1) /= nx .or. size(u,2) /= ny .or. &
        size(v,1) /= nx .or. size(v,2) /= ny) then
      write(*,*) 'ERROR: Dimension mismatch in quiver_2d'
      write(*,*) 'Expected:', nx, 'x', ny
      write(*,*) 'Got u:', size(u,1), 'x', size(u,2)
      write(*,*) 'Got v:', size(v,1), 'x', size(v,2)
      return
    end if
    
    call setup_plot_environment(PLT)
    
    u_norm = u
    v_norm = v
    
    if (opts%normalize) then
      call normalize_vectors(u_norm, v_norm, opts%min_magnitude)
    end if
    
    color_index = get_color_index(opts%color)
    call plcol0(color_index)
    
    do i = 1, nx
      do j = 1, ny
        if (sqrt(u(i,j)**2 + v(i,j)**2) > opts%min_magnitude) then
          dx = u_norm(i,j) * opts%scale
          dy = v_norm(i,j) * opts%scale
          call draw_arrow(x(i), y(j), dx, dy, opts%head_size, opts%pivot)
        end if
      end do
    end do
    
  end subroutine quiver_2d

  subroutine quiver_1d(PLT, x, y, u, v, options)
    class(type_plplot), intent(inout) :: PLT
    real(dp), intent(in) :: x(:), y(:), u(:), v(:)
    type(quiver_options), intent(in), optional :: options
    type(quiver_options) :: opts
    real(dp) :: u_norm(size(u)), v_norm(size(v))
    integer :: n, i, color_index
    real(dp) :: dx, dy
    
    opts = quiver_options()  ! Default values
    if (present(options)) opts = options
    
    n = size(x)
    
    if (size(y) /= n .or. size(u) /= n .or. size(v) /= n) then
      write(*,*) 'ERROR: Dimension mismatch in quiver_1d'
      write(*,*) 'Expected size:', n
      write(*,*) 'Got y:', size(y), 'u:', size(u), 'v:', size(v)
      return
    end if
    
    call setup_plot_environment(PLT, x, y)
    
    u_norm = u
    v_norm = v
    
    if (opts%normalize) then
      call normalize_vectors_1d(u_norm, v_norm, opts%min_magnitude)
    end if
    
    color_index = get_color_index(opts%color)
    call plcol0(color_index)
    
    do i = 1, n
      if (sqrt(u(i)**2 + v(i)**2) > opts%min_magnitude) then
        dx = u_norm(i) * opts%scale
        dy = v_norm(i) * opts%scale
        call draw_arrow(x(i), y(i), dx, dy, opts%head_size, opts%pivot)
      end if
    end do
    
  end subroutine quiver_1d

  ! ==== Quiver helper functions ====
  
  subroutine normalize_vectors(u, v, min_magnitude)
    real(dp), intent(inout) :: u(:,:), v(:,:)
    real(dp), intent(in) :: min_magnitude
    real(dp) :: magnitude
    integer :: i, j
    
    do i = 1, size(u, 1)
      do j = 1, size(u, 2)
        magnitude = sqrt(u(i,j)**2 + v(i,j)**2)
        if (magnitude > min_magnitude) then
          u(i,j) = u(i,j) / magnitude
          v(i,j) = v(i,j) / magnitude
        else
          u(i,j) = 0.0_dp
          v(i,j) = 0.0_dp
        end if
      end do
    end do
  end subroutine normalize_vectors

  subroutine normalize_vectors_1d(u, v, min_magnitude)
    real(dp), intent(inout) :: u(:), v(:)
    real(dp), intent(in) :: min_magnitude
    real(dp) :: magnitude
    integer :: i
    
    do i = 1, size(u)
      magnitude = sqrt(u(i)**2 + v(i)**2)
      if (magnitude > min_magnitude) then
        u(i) = u(i) / magnitude
        v(i) = v(i) / magnitude
      else
        u(i) = 0.0_dp
        v(i) = 0.0_dp
      end if
    end do
  end subroutine normalize_vectors_1d

  subroutine draw_arrow(x, y, dx, dy, head_size, pivot)
    real(dp), intent(in) :: x, y, dx, dy, head_size
    integer, intent(in) :: pivot
    real(dp) :: x_start, y_start, x_end, y_end
    real(dp) :: angle, head_dx1, head_dy1, head_dx2, head_dy2
    real(dp) :: shift_x, shift_y
    
    ! Calculate pivot shift
    select case(pivot)
    case(PIVOT_TAIL)
      shift_x = 0.0_dp
      shift_y = 0.0_dp
    case(PIVOT_MID)
      shift_x = -dx / 2.0_dp
      shift_y = -dy / 2.0_dp
    case(PIVOT_TIP)
      shift_x = -dx
      shift_y = -dy
    case default
      shift_x = 0.0_dp
      shift_y = 0.0_dp
    end select
    
    x_start = x + shift_x
    y_start = y + shift_y
    x_end = x_start + dx
    y_end = y_start + dy
    
    ! Draw arrow shaft
    call pljoin(x_start, y_start, x_end, y_end)
    
    ! Draw arrow head
    if (abs(dx) > 1e-10_dp .or. abs(dy) > 1e-10_dp) then
      angle = atan2(dy, dx)
      
      head_dx1 = head_size * cos(angle + 2.5_dp)
      head_dy1 = head_size * sin(angle + 2.5_dp)
      head_dx2 = head_size * cos(angle - 2.5_dp)
      head_dy2 = head_size * sin(angle - 2.5_dp)
      
      call pljoin(x_end, y_end, x_end + head_dx1, y_end + head_dy1)
      call pljoin(x_end, y_end, x_end + head_dx2, y_end + head_dy2)
    end if
  end subroutine draw_arrow

  ! ==== Option constructor functions ====
  
  function line_opts(width, color, style, label) result(opts)
    real(dp), optional :: width
    character(*), intent(in), optional :: color
    character(*), intent(in), optional :: label
    character(*), intent(in), optional :: style
    type(line_options) :: opts

    if (present(width)) opts%width = width
    if (present(color)) opts%color = trim(color)
    if (present(style)) opts%style = trim(style)
    if (present(label)) opts%label = trim(label)
  end function line_opts

  function quiver_opts(scale, head_size, color, normalize, min_magnitude, &
                       pivot) result(opts)
    real(dp), optional :: scale, head_size, min_magnitude
    integer, optional :: pivot
    character(len=*), optional :: color
    logical, optional :: normalize
    type(quiver_options) :: opts
    
    if (present(scale)) opts%scale = scale
    if (present(head_size)) opts%head_size = head_size
    if (present(color)) opts%color = trim(color)
    if (present(normalize)) opts%normalize = normalize
    if (present(min_magnitude)) opts%min_magnitude = min_magnitude
    if (present(pivot)) opts%pivot = pivot
  end function quiver_opts

  ! ==== Legend functions ====
  
  subroutine plplot_legend(PLT)
    class(type_plplot), intent(inout) :: PLT
    call plt_legend_show(PLT%ax%axis_legend)
  end subroutine plplot_legend

  subroutine plaxis_legend(AX)
    class(type_plaxis), intent(inout) :: AX
    call plt_legend_show(AX%axis_legend)
  end subroutine plaxis_legend

  subroutine plt_legend_show(L)
    class(type_legend), intent(inout) :: L
    integer :: nlegend

    nlegend = L%nlegend
    if (nlegend <= 0) return

    call pllegend(L%legend_width, L%legend_height, &
      L%opt, L%position, &
      L%xoffset, L%yoffset, L%plot_width, &
      L%bg_color, L%bb_color, L%bb_style, &
      L%nrow, L%ncolumn, &
      L%opt_array(1:nlegend), &
      L%text_offset, L%text_scale, L%text_spacing, L%text_justification, &
      L%text_colors(1:nlegend), L%text_labels(1:nlegend), &
      L%box_colors(1:nlegend), L%box_patterns(1:nlegend), &
      L%box_scales(1:nlegend), L%box_line_widths(1:nlegend), &
      L%line_colors(1:nlegend), L%line_styles(1:nlegend), &
      L%line_widths(1:nlegend), &
      L%symbol_colors(1:nlegend), L%symbol_scales(1:nlegend), &
      L%symbol_numbers(1:nlegend), L%symbols(1:nlegend))
  end subroutine plt_legend_show

end module module_plplot
