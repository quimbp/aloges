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

type type_plaxis
  logical                                  :: new_axis = .True.
  logical                                  :: env_done = .False.
  logical                                  :: xlimits_set = .False.
  logical                                  :: ylimits_set = .False.
  logical                                  :: show_grid = .False.
  logical                                  :: axes_equal = .False.
  logical                                  :: xlog = .False.
  logical                                  :: ylog = .False.
  real(dp)                                 :: xmin,xmax,ymin,ymax
  real(dp)                                 :: fontsize = 1.0D0
  character(len=50)                        :: xlabel = ''
  character(len=50)                        :: ylabel = ''
  character(len=50)                        :: title = ''
  character(len=50)                        :: fontcolor  = 'black'
  character(len=50)                        :: foreground = 'black'
  contains
    procedure                   :: new           => plaxis_new
    procedure                   :: ylim          => plaxis_ylim
    procedure                   :: xlim          => plaxis_xlim
    procedure                   :: labels        => plaxis_labels
end type type_plaxis


type type_plplot
  type(type_plaxis)                        :: ax
  logical                                  :: multiplot = .False.
  logical                                  :: init_done = .False.
  logical                                  :: env_done = .False.
  integer                                  :: color_index = 0
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

    if (present(background)) then
      RGB = color2rgb(background)
      print*, 'RGB = ', RGB
      call plscol0(0, RGB(1), RGB(2), RGB(3))  ! Background
    else
      ! ...            R    G    B
      call plscol0(0, 255, 255, 255)  ! Background: white
    endif

    if (present(foreground)) then
      RGB = color2rgb(foreground)
      print*, 'RGB = ', RGB
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
 
  end subroutine plaxis_new
  ! ...
  ! ===================================================================
  ! ...
  function color2rgb(name) result(RGB)

    character(len=*), intent(in)              :: name
    integer, dimension(3)                     :: RGB

    select case(name)
    case('white','WHITE','White')
      RGB = [255, 255, 255 ]
    case('black','BLACK','Black')
      RGB = [0, 0, 0 ]
    case('red','RED','Red')
      RGB = [255, 0, 0 ]
    case('green','GREEN','Green')
      RGB = [0, 255, 0 ]
    case('blue','BLUE','Blue')
      RGB = [0, 0, 255 ]
    case default
      stop 'invalid color name'
    end select 

  end function color2rgb
  ! ...
  ! ===================================================================
  ! ...
  subroutine plplot_show(PLT)

    class(type_plplot), intent(inout)              :: PLT

    ! Close PLplot library
    call plend()

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
  subroutine plot_y(PLT,y,color,linewidth)
  ! Case 1: only Y(:)

    class(type_plplot), intent(inout)         :: PLT
    real(dp), dimension(:), intent(in)        :: y
    character(len=*), intent(in), optional    :: color
    real(dp), intent(in), optional            :: linewidth

    ! ... Local variables
    ! ...
    integer i
    real(dp) x(size(y))

    do i=1,size(y)
      x(i) = 1.0D0*i
    enddo

    call plot_xy(PLT,x,y,color,linewidth)

  end subroutine plot_y
  ! ...
  ! ===================================================================
  ! ...
  subroutine plot_xy(PLT,x,y,color,linewidth)
  ! Case 2: X(:), Y(:)

    class(type_plplot), intent(inout)         :: PLT
    real(dp), dimension(:), intent(in)        :: x
    real(dp), dimension(:), intent(in)        :: y
    character(len=*), intent(in), optional    :: color
    real(dp), intent(in), optional            :: linewidth

    ! ... Local variables
    ! ...
    logical new_environment
    integer ibox,RGB(3)
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
      RGB = color2rgb(PLT%ax%foreground)
      call plscol0(10, RGB(1), RGB(2), RGB(3))  ! Background
      call plcol0(10)
      call plenv(xmin,xmax,ymin,ymax,0,ibox)
    endif

    ! ... Labels (black)
    ! ...
    RGB = color2rgb(PLT%ax%fontcolor)
    call plscol0(10, RGB(1), RGB(2), RGB(3))  ! Background
    call plcol0(10)
    call plschr(0.0D0,PLT%ax%fontsize) 
    call pllab( trim(PLT%ax%xlabel), trim(PLT%ax%ylabel), trim(PLT%ax%title) )

    if (present(color)) then
      RGB = color2rgb(color)
      call plscol0(10, RGB(1), RGB(2), RGB(3))  ! Background
      call plcol0(10)
    else
      PLT%color_index = PLT%color_index + 1
      if (PLT%color_index.GT.10) PLT%color_index = 1
      call plcol0(PLT%color_index)
    endif

    ! ... Line width
    if (present(linewidth)) then
      call plwidth(linewidth)
    else
      call plwidth(PLT%linewidth)
    endif
    
    call plline(x,y)

    
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

    

    

end module module_plplot

