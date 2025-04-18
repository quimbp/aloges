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
! ... rk4                                                                  !
! ...                                                                      !
! -------------------------------------------------------------------------!

module module_sysdyn

use, intrinsic :: IEEE_ARITHMETIC, ONLY : IEEE_VALUE, IEEE_QUIET_NAN
use module_types
use module_constants
use module_tools, only: unitfree

implicit none

type model_lorenz63
  ! Lorenz 1963 model
  ! p[1] = sigma = 10.0D0
  ! p[2] = rho = 28.0D0
  ! p[3] = beta = 8.0D0/3.0D0
  character(len=20)                           :: name = 'LORENZ 1963 MODEL'
  integer                                     :: nsys = 3
  integer                                     :: npar = 3
  real(dp), dimension(3)                      :: p = [10.0D0, 28.0D0, 8.0D0/3.0D0]
  character(len=6), dimension(3)              :: pname = ['sigma ','rho   ','beta  ']
  real(dp)                                    :: dt = 0.01D0
  real(dp)                                    :: to = 0.0D0
  real(dp), dimension(3)                      :: xo = [1.508870D0,-1.531271D0,25.460910D0]
  real(dp), dimension(:), allocatable         :: t
  real(dp), dimension(:,:), allocatable       :: x
  contains
    procedure                   :: config        => l63_config
    procedure                   :: set           => l63_set   
    procedure                   :: run           => l63_run
    procedure                   :: save          => l63_save
    procedure                   :: plot          => l63_plot
    procedure                   :: show          => l63_show
end type model_lorenz63


type model_LotkaVolterra
  character(len=20)                           :: name = 'LOTKA-VOLTERRA PZ'
  integer                                     :: nsys = 2
  integer                                     :: npar = 4
  real(dp), dimension(4)                      :: p = [1.0D0, 3.0D0, 2.0D0, 0.30D0 ]
  character(len=4), dimension(4)              :: pname = ['p(1)','p(2)','p(3)','p(4)']
  real(dp)                                    :: dt = 0.1D0
  real(dp)                                    :: to = 0.0D0
  real(dp), dimension(2)                      :: xo = [0.1D0, 0.1D0 ]
  real(dp), dimension(:), allocatable         :: t
  real(dp), dimension(:,:), allocatable       :: x
  contains
    procedure                   :: config        => LV_config
    procedure                   :: set           => LV_set   
    procedure                   :: run           => LV_run
    procedure                   :: save          => LV_save
    procedure                   :: plot          => LV_plot
    procedure                   :: show          => LV_show
end type model_LotkaVolterra


contains
! ...
! =====================================================================
! =====================================================================
! ...
  function rk1(param,t,x,dt,F) result(xn)
  ! ... Runge-Kutta order 1 schema for time integration of the ODE: dx/dt = F
  ! ... Euler's method

    real(dp), dimension(:), intent(in)     :: param
    real(dp), intent(in)                   :: t
    real(dp), dimension(:), intent(in)     :: x
    real(dp), intent(in)                   :: dt
    real(dp), dimension(size(x))           :: xn

    interface
      pure function F(p, t, x) result(y)
        implicit none
        real(8), dimension(:), intent(in) :: p
        real(8), intent(in)               :: t
        real(8), dimension(:), intent(in) :: x
        real(8), dimension(size(x))       :: y
      end function f
    end interface

    ! ... Local variables:
    ! ...
    real(dp), dimension(size(x))           :: k1

    k1 = dt * F(param, t, x)

    xn(:) = x(:) + k1(:)

  end function rk1
  ! ...
  ! ...
  ! ===================================================================
  ! ...
  function rk2(param,t,x,dt,F) result(xn)
  ! ... Runge-Kutta order 2 schema for time integration of the ODE: dx/dt = F
  ! ... Mid-point's method

    real(dp), dimension(:), intent(in)     :: param
    real(dp), intent(in)                   :: t
    real(dp), dimension(:), intent(in)     :: x
    real(dp), intent(in)                   :: dt
    real(dp), dimension(size(x))           :: xn

    interface
      pure function F(p, t, x) result(y)
        implicit none
        real(8), dimension(:), intent(in) :: p
        real(8), intent(in)               :: t
        real(8), dimension(:), intent(in) :: x
        real(8), dimension(size(x))       :: y
      end function f
    end interface

    ! ... Local variables:
    ! ...
    real(dp), dimension(size(x))           :: k1,k2

    k1 = dt * F(param, t, x)
    k2 = dt * F(param, t+0.5D0*dt, x+0.5D0*k1)

    xn(:) = x(:) + k2(:)

  end function rk2
  ! ...
  ! ===================================================================
  ! ...
  function rk3(param,t,x,dt,F) result(xn)
  ! ... Runge-Kutta order 3 schema for time integration of the ODE: dx/dt = F
  ! ... Heun's or Nystrom method

    real(dp), dimension(:), intent(in)     :: param
    real(dp), intent(in)                   :: t
    real(dp), dimension(:), intent(in)     :: x
    real(dp), intent(in)                   :: dt
    real(dp), dimension(size(x))           :: xn

    interface
      pure function F(p, t, x) result(y)
        implicit none
        real(8), dimension(:), intent(in) :: p
        real(8), intent(in)               :: t
        real(8), dimension(:), intent(in) :: x
        real(8), dimension(size(x))       :: y
      end function f
    end interface

    ! ... Local variables:
    ! ...
    real(dp), dimension(size(x))           :: k1,k2,k3

    k1 = dt * F(param, t, x)
    k2 = dt * F(param, t+0.5D0*dt, x+0.5*k1)
    k3 = dt * F(param, t+dt, x-k1+2.0D0*k2)

    xn(:) = x(:) + (k1(:) + 4.0D0*k2(:) + k3(:)) / 6.0D0

  end function rk3
  ! ...
  ! ===================================================================
  ! ...
  function rk4(param,t,x,dt,F) result(xn)
  ! ... Runge-Kutta order 4 schema for time integration of the ODE: dx/dt = F

    real(dp), dimension(:), intent(in)     :: param
    real(dp), intent(in)                   :: t
    real(dp), dimension(:), intent(in)     :: x
    real(dp), intent(in)                   :: dt
    real(dp), dimension(size(x))           :: xn

    interface
      pure function F(p, t, x) result(y)
        implicit none
        real(8), dimension(:), intent(in) :: p
        real(8), intent(in)               :: t
        real(8), dimension(:), intent(in) :: x
        real(8), dimension(size(x))       :: y
      end function f
    end interface

    ! ... Local variables:
    ! ...
    real(dp), dimension(size(x))           :: k1,k2,k3,k4

    k1 = dt * F(param, t, x)
    k2 = dt * F(param, t+0.5D0*dt, x+0.5*k1)
    k3 = dt * F(param, t+0.5D0*dt, x+0.5*k2)
    k4 = dt * F(param, t+dt, x+k3)

    xn(:) = x(:) + (k1(:) + 2.0D0*k2(:) + 2.0D0*k3(:) + k4(:)) / 6.0D0

  end function rk4
  ! ...
  ! ===================================================================
  ! ...
  function rk5(param,t,x,dt,F) result(xn)
  ! ... Runge-Kutta order 5 schema for time integration of the ODE: dx/dt = F
  ! ... Cash-Karp method

    real(dp), dimension(:), intent(in)     :: param
    real(dp), intent(in)                   :: t
    real(dp), dimension(:), intent(in)     :: x
    real(dp), intent(in)                   :: dt
    real(dp), dimension(size(x))           :: xn

    interface
      pure function F(p, t, x) result(y)
        implicit none
        real(8), dimension(:), intent(in) :: p
        real(8), intent(in)               :: t
        real(8), dimension(:), intent(in) :: x
        real(8), dimension(size(x))       :: y
      end function f
    end interface

    ! ... Local variables:
    ! ...
    real(dp), dimension(size(x))           :: k1,k2,k3,k4,k5,k6

    k1 = dt * F(param, t, x)
    k2 = dt * F(param, t + dt * 0.2d0, x + k1 * 0.2d0)
    k3 = dt * F(param, t + dt * 0.3d0, x + k1 * 3.0d0 / 40.0d0 + k2 * 9.0d0 / 40.0d0)
    k4 = dt * F(param, t + dt * 0.6d0, x + k1 * 0.3d0 - k2 * 0.9d0 + k3 * 1.2d0)
    k5 = dt * F(param, t + dt, x - k1 * 11.0d0 / 54.0d0 + k2 * 2.5d0 - k3 * 70.0d0 / 27.0d0 + k4 * 35.0d0 / 27.0d0)
    k6 = dt * F(param, t + dt * 0.875d0, x + k1 * 1631.0d0 / 55296.0d0 + k2 * 175.0d0 / 512.0d0 + &
              k3 * 575.0d0 / 13824.0d0 + k4 * 44275.0d0 / 110592.0d0 + k5 * 253.0d0 / 4096.0d0)

    xn = x + (k1 * 37.0d0 / 378.0d0 + k3 * 250.0d0 / 621.0d0 + k4 * 125.0d0 / 594.0d0 + &
            k6 * 512.0d0 / 1771.0d0)

  end function rk5
  ! ...
  ! ===================================================================
  ! ...
  subroutine ModelSave(filename,t,X)

    character(len=*), intent(in)          :: filename
    real(dp), dimension(:), intent(in)    :: t
    real(dp), dimension(:,:), intent(in)  :: X

    ! ... Local variables
    ! ...
    integer iu,step,nsteps

    iu = unitfree()
    open(iu,file=filename,status='unknown')
    rewind(iu)
    do step=1,size(t)
      write(iu,'(4F10.4)') t(step), X(:,step)
    enddo
    close(iu)

  end subroutine ModelSave
  ! ...
  ! ===================================================================
  ! ===================================================================
  ! ...
  pure function l63_rhs(p,t,x) result(dxdt)
  ! ...
  ! ... The right hand side equation has to be a vector function f(:) = f(p,t,x)
  ! ... Parameters are passed via 

    real(dp), dimension(:), intent(in)           :: p        ! Model parameters
    real(dp), intent(in)                         :: t        ! time or independent var
    real(dp), dimension(:), intent(in)           :: x        ! state vector
    real(dp), dimension(size(x))                 :: dxdt     ! rate of change

    dxdt(1) = p(1)*(x(2)-x(1))
    dxdt(2) = x(1)*(p(2)-x(3)) - x(2)
    dxdt(3) = x(1)*x(2) - p(3)*x(3)

  end function l63_rhs
  ! ...
  ! ===================================================================
  ! ...
  subroutine l63_config(MODEL,dt,p)

    class(model_lorenz63), intent(inout)                   :: MODEL
    real(dp), intent(in), optional                         :: dt
    real(dp), dimension(MODEL%npar), intent(in), optional  :: p

    if (present(dt))    MODEL%dt = dt
    if (present(p))     MODEL%p(:) = p(:)

  end subroutine l63_config
  ! ...
  ! ===================================================================
  ! ...
  subroutine l63_set(MODEL,to,xo)

    class(model_lorenz63), intent(inout)          :: MODEL
    real(dp), intent(in), optional               :: to
    real(dp), dimension(3), intent(in), optional :: xo

    if (present(to))    MODEL%to = to
    if (present(xo))    MODEL%xo(:) = xo(:)

  end subroutine l63_set
  ! ...
  ! ===================================================================
  ! ...
  subroutine l63_run(MODEL,niter)

    class(model_lorenz63), intent(inout)         :: MODEL
    integer, intent(in), optional                :: niter

    ! ... Local variables
    ! ...
    integer step,nsteps
    real(dp)                                     :: t
    real(dp), dimension(MODEL%nsys)              :: x,xn

    if (present(niter)) then
      nsteps = niter
    else
      nsteps = 1
    endif

    t = MODEL%to
    x = MODEL%xo(:)
    
    if (allocated(MODEL%t)) deallocate(MODEL%t)
    if (allocated(MODEL%x)) deallocate(MODEL%x)
    allocate(MODEL%t(0:nsteps))
    allocate(MODEL%x(MODEL%nsys,0:nsteps))

    ! ... Store initial conditions
    ! ...
    MODEL%t(0)   = t
    MODEL%x(:,0) = x(:)
    MODEL_LOOP: do step=1,nsteps

      ! ... Update the system
      ! ...
      print*, 'rk5'
      xn = rk5(MODEL%p,t,x,MODEL%dt,l63_rhs)
      t = t + MODEL%dt
      x(:) = xn(:)

      MODEL%t(step)   = t
      MODEL%x(:,step) = x(:)

    enddo MODEL_LOOP

  end subroutine l63_run
  ! ...
  ! ===================================================================
  ! ...
  subroutine l63_show(MODEL,label)

    class(model_lorenz63), intent(in)      :: MODEL
    character(len=*), intent(in), optional :: label

    ! ... Local variables
    ! ...
    integer i,n


    WRITE(*,'(A)') "==========================================="
    WRITE(*,'(A)') "         LORENZ 1963 DYNAMICAL MODEL"
    WRITE(*,'(A)') "==========================================="

    if (present(label)) then
      WRITE(*,'(A)') ""
      WRITE(*,'(A)') ">> " // trim(label)
    endif
    WRITE(*,'(A)') ""
    WRITE(*,'(A)') ">> Model Parameters:"
    WRITE(*,'(A,F10.6)') "    Sigma (σ)       = ", MODEL%p(1)
    WRITE(*,'(A,F10.6)') "    Rho (ρ)         = ", MODEL%p(2)
    WRITE(*,'(A,F10.6)') "    Beta (β)        = ", MODEL%p(3)
    WRITE(*,'(A,F10.6)') "    Time Step (dt)  = ", MODEL%dt

    WRITE(*,'(A)') ""
    WRITE(*,'(A)') ">> Initial Conditions:"
    WRITE(*,'(A,F10.6)') "    t₀              = ", MODEL%to
    WRITE(*,'(A,F10.6)') "    x₀              = ", MODEL%xo(1)
    WRITE(*,'(A,F10.6)') "    y₀              = ", MODEL%xo(2)
    WRITE(*,'(A,F10.6)') "    z₀              = ", MODEL%xo(3)

    if (allocated(MODEL%t)) then
      n = size(MODEL%t) - 1

      WRITE(*,'(A)') ""
      WRITE(*,'(A)') ">> Simulation Status:"
      WRITE(*,'(A,I10)') "    Time Steps      = ",   n+1
      WRITE(*,'(A,F10.6)') "    Start Time      = ", MODEL%t(0)
      WRITE(*,'(A,F10.6)') "    End Time        = ", MODEL%t(n)

      WRITE(*,'(A)') ""
      WRITE(*,'(A)') ">> Final Values:"
      WRITE(*,'(A,F10.6)') "    x               = ", MODEL%x(1,n)
      WRITE(*,'(A,F10.6)') "    y               = ", MODEL%x(2,n)
      WRITE(*,'(A,F10.6)') "    z               = ", MODEL%x(3,n)

    else

      WRITE(*,'(A)') ""
      WRITE(*,'(A)') ">> Simulation Status:"
      WRITE(*,'(A)') "    [ No simulation has been run yet ]"

    endif

    WRITE(*,'(A)') "==========================================="


!    if (present(label)) write(*,'(A)') trim(label)
!    write(*,'(A)') trim(MODEL%name)
!    write(*,'(A)') 'Model parameters: ' 
!    do i=1,MODEL%npar
!      write(*,*) trim(MODEL%pname(i))//' = ', MODEL%p(i)
!    enddo
!    write(*,*) 'dt = ', MODEL%dt
!    write(*,'(A)') 'Initial conditions:'
!    write(*,*) 'to: ', MODEL%to
!    write(*,*) 'xo: ', MODEL%xo
!    if (allocated(MODEL%t)) then
!      write(*,'(A)') 'Latest simulation:'
!      write(*,*) 'Time steps: ', size(MODEL%t)
!      write(*,*) 'From: ', MODEL%t(0)
!      write(*,*) 'to:   ', MODEL%t(size(MODEL%t)-1)
!    else
!      write(*,'(A)') 'No simulation yet'
!    endif
! 
      
    
  end subroutine l63_show
  ! ...
  ! ===================================================================
  ! ...
  subroutine l63_save(MODEL,filename)

    class(model_lorenz63), intent(in)    :: MODEL
    character(len=*), intent(in)        :: filename

    call ModelSave(filename,MODEL%t,MODEL%x)

  end subroutine l63_save
  ! ...
  ! ===================================================================
  ! ...
  subroutine l63_plot(MODEL,filename)

    class(model_lorenz63), intent(in)           :: MODEL
    character(len=*), intent(in), optional     :: filename

    ! ... Local variables
    ! ...
    logical saveit
    integer iu

    call MODEL%save('mylorenzdat.dat')

    saveit = .false.
    if (present(filename)) saveit = .true.

    iu = unitfree()
    open(iu,file='myl63plot.py',status='unknown')
    write(iu,'(A)') "import matplotlib"
    write(iu,'(A)') "import matplotlib.pyplot as plt"

    write(iu,'(A)') "font = {'family' : 'DejaVu Sans',"
    write(iu,'(A)') "        'weight'   : 'bold',"
    write(iu,'(A)') "        'size'     : 14}"

    write(iu,'(A)') "fig = {'figsize': (12,8),"
    write(iu,'(A)') "       'dpi'    : 100}"

    write(iu,'(A)') "savefig = {'dpi'       : 300,"
    write(iu,'(A)') "           'directory' : './'}"
    !write(iu,'(A)') "           'format'    : 'pdf'}"

    write(iu,'(A)') "matplotlib.rc('font',**font)"
    write(iu,'(A)') "matplotlib.rc('figure',**fig)"
    write(iu,'(A)') "matplotlib.rc('savefig',**savefig)"

    write(iu,'(A)') "F = open('mylorenzdat.dat','r')"
    write(iu,'(A)') "time = []"
    write(iu,'(A)') "x = []"
    write(iu,'(A)') "y = []"
    write(iu,'(A)') "z = []"
    write(iu,'(A)') "for line in F:"
    write(iu,'(A)') "  tmp = line.split()"
    write(iu,'(A)') "  time.append(float(tmp[0]))"
    write(iu,'(A)') "  x.append(float(tmp[1]))"
    write(iu,'(A)') "  y.append(float(tmp[2]))"
    write(iu,'(A)') "  z.append(float(tmp[3]))"
    write(iu,'(A)') "F.close()"

    write(iu,'(A)') "fig = plt.figure()"

    write(iu,'(A)') "fig.add_subplot(311)"
    write(iu,'(A)') "plt.plot(time,x,'-r',linewidth=2)"
    write(iu,'(A)') "plt.ylabel('x',fontweight='bold',fontsize=14)"
    write(iu,'(A)') "plt.grid(True)"
    write(iu,'(A)') "plt.gca().axes.get_xaxis().set_ticklabels([])"

    write(iu,'(A)') "fig.add_subplot(312)"
    write(iu,'(A)') "plt.plot(time,y,'-b',linewidth=2)"
    write(iu,'(A)') "plt.ylabel('y',fontweight='bold',fontsize=14)"
    write(iu,'(A)') "plt.grid(True)"
    write(iu,'(A)') "plt.gca().axes.get_xaxis().set_ticklabels([])"

    write(iu,'(A)') "fig.add_subplot(313)"
    write(iu,'(A)') "plt.plot(time,z,'-g',linewidth=2)"
    write(iu,'(A)') "plt.xlabel('time',fontweight='bold',fontsize=14)"
    write(iu,'(A)') "plt.ylabel('z',fontweight='bold',fontsize=14)"
    write(iu,'(A)') "plt.grid(True)"

    if (saveit) then
      write(iu,'(A)') "fig.savefig('"//trim(filename)//"',bbox_inches='tight')"
    else
      write(iu,'(A)') "plt.show()"
    endif
    close(iu)
    call system('python myl63plot.py')



  end subroutine l63_plot
  ! ...
  ! ===================================================================
  ! ===================================================================
  ! ...
  pure function LV_rhs(p,t,x) result(dxdt)
  ! ...
  ! ... The right hand side equation has to be a vector function f(:) = f(p,t,x)
  ! ... Parameters are passed via 

    real(dp), dimension(:), intent(in)           :: p        ! Model parameters
    real(dp), intent(in)                         :: t        ! time or independent var
    real(dp), dimension(:), intent(in)           :: x        ! state vector
    real(dp), dimension(size(x))                 :: dxdt     ! rate of change

    dxdt(1) = p(1)*(1.0D0-x(1)/p(2))*x(1) - p(3)*x(1)*x(2)
    dxdt(2) = p(3)*x(1)*x(2) - p(4)*x(2)

  end function LV_rhs
  ! ...
  ! ===================================================================
  ! ...
  subroutine LV_config(MODEL,dt,p)

    class(model_LotkaVolterra), intent(inout)              :: MODEL
    real(dp), intent(in), optional                        :: dt
    real(dp), dimension(MODEL%npar), intent(in), optional :: p

    if (present(dt))    MODEL%dt = dt
    if (present(p))     MODEL%p(:) = p(:)

  end subroutine LV_config
  ! ...
  ! ===================================================================
  ! ...
  subroutine LV_set(MODEL,to,xo)

    class(model_LotkaVolterra), intent(inout)              :: MODEL
    real(dp), intent(in), optional                        :: to
    real(dp), dimension(MODEL%nsys), intent(in), optional :: xo

    if (present(to))    MODEL%to = to
    if (present(xo))    MODEL%xo(:) = xo(:)

  end subroutine LV_set
  ! ...
  ! ===================================================================
  ! ...
  subroutine LV_run(MODEL,niter)

    class(model_LotkaVolterra), intent(inout)     :: MODEL
    integer, intent(in), optional                 :: niter

    ! ... Local variables
    ! ...
    integer step,nsteps
    real(dp)                                     :: t
    real(dp), dimension(MODEL%nsys)              :: x,xn
    
    if (present(niter)) then
      nsteps = niter
    else
      nsteps = 1
    endif

    t = MODEL%to
    x = MODEL%xo(:)
    
    if (allocated(MODEL%t)) deallocate(MODEL%t)
    if (allocated(MODEL%x)) deallocate(MODEL%x)
    allocate(MODEL%t(0:nsteps))
    allocate(MODEL%x(MODEL%nsys,0:nsteps))

    ! ... Store initial conditions
    ! ...
    MODEL%t(0)   = t
    MODEL%x(:,0) = x(:)
    MODEL_LOOP: do step=1,nsteps

      ! ... Update the system
      ! ...
      xn = rk4(MODEL%p,t,x,MODEL%dt,LV_rhs)
      t = t + MODEL%dt
      x(:) = xn(:)

      MODEL%t(step)   = t
      MODEL%x(:,step) = x(:)

    enddo MODEL_LOOP

  end subroutine LV_run
  ! ...
  ! ===================================================================
  ! ...
  subroutine LV_show(MODEL,label)

    class(model_LotkaVolterra), intent(in)   :: MODEL
    character(len=*), intent(in), optional   :: label

    ! ... Local variables
    ! ...
    integer i

    if (present(label)) write(*,'(A)') trim(label)
    write(*,'(A)') trim(MODEL%name)
    write(*,'(A)') 'Model parameters: ' 
    do i=1,MODEL%npar
      write(*,*) trim(MODEL%pname(i))//' = ', MODEL%p(i)
    enddo
    write(*,*) 'dt = ', MODEL%dt
    
  end subroutine LV_show
  ! ...
  ! ===================================================================
  ! ...
  subroutine LV_save(MODEL,filename)

    class(model_LotkaVolterra), intent(in)    :: MODEL
    character(len=*), intent(in)        :: filename

    call ModelSave(filename,MODEL%t,MODEL%x)

  end subroutine LV_save
  ! ...
  ! ===================================================================
  ! ...
  subroutine LV_plot(MODEL,filename)

    class(model_LotkaVolterra), intent(in)      :: MODEL
    character(len=*), intent(in), optional     :: filename

    ! ... Local variables
    ! ...
    logical saveit
    integer iu

    call MODEL%save('mylvdat.dat')

    saveit = .false.
    if (present(filename)) saveit = .true.

    iu = unitfree()
    open(iu,file='myLVplot.py',status='unknown')
    write(iu,'(A)') "import matplotlib"
    write(iu,'(A)') "import matplotlib.pyplot as plt"

    write(iu,'(A)') "font = {'family' : 'DejaVu Sans',"
    write(iu,'(A)') "        'weight'   : 'bold',"
    write(iu,'(A)') "        'size'     : 14}"

    write(iu,'(A)') "fig = {'figsize': (12,8),"
    write(iu,'(A)') "       'dpi'    : 100}"

    write(iu,'(A)') "savefig = {'dpi'       : 300,"
    write(iu,'(A)') "           'directory' : './'}"
    !write(iu,'(A)') "           'format'    : 'pdf'}"

    write(iu,'(A)') "matplotlib.rc('font',**font)"
    write(iu,'(A)') "matplotlib.rc('figure',**fig)"
    write(iu,'(A)') "matplotlib.rc('savefig',**savefig)"

    write(iu,'(A)') "F = open('mylvdat.dat','r')"
    write(iu,'(A)') "time = []"
    write(iu,'(A)') "x = []"
    write(iu,'(A)') "y = []"
    write(iu,'(A)') "for line in F:"
    write(iu,'(A)') "  tmp = line.split()"
    write(iu,'(A)') "  time.append(float(tmp[0]))"
    write(iu,'(A)') "  x.append(float(tmp[1]))"
    write(iu,'(A)') "  y.append(float(tmp[2]))"
    write(iu,'(A)') "F.close()"

    write(iu,'(A)') "fig = plt.figure()"

    write(iu,'(A)') "fig.add_subplot(111)"
    write(iu,'(A)') "plt.plot(time,x,'-k',label='Phytoplankton',linewidth=2)"
    write(iu,'(A)') "plt.plot(time,y,'--k',label='Zooplankton',linewidth=2)"
    write(iu,'(A)') "plt.ylim(0,1.4)"
    write(iu,'(A)') "plt.xlabel('Time (days)',fontweight='bold',fontsize=14)"
    write(iu,'(A)') "plt.ylabel('Plankton (mmol/m3)',fontweight='bold',fontsize=14)"
    write(iu,'(A)') "plt.legend()"
    write(iu,'(A)') "plt.grid(True)"

    if (saveit) then
      write(iu,'(A)') "fig.savefig('"//trim(filename)//"',bbox_inches='tight')"
    else
      write(iu,'(A)') "plt.show()"
    endif
    close(iu)
    call system('python myLVplot.py')

  end subroutine LV_plot
  ! ...
  ! ===================================================================
  ! ===================================================================
  ! ...
end module module_sysdyn
