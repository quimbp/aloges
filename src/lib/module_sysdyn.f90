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
! ...                                                                      !
! ... Model Lorenz 1963                                                    !
! ... Model Lotka-Volterra PZ                                              !
! ... rk1, rk2, rk3, rk4, rk5                                              !
! ... srk4 (stochastic RK4)                                                !
! ... ab2 (Adams-Bashforth 2)                                              !
! ... leapfrog (Symplectic)                                                !
! ... rkf45 (Runge-Kutta-Fehlberg)                                         !
! ... bs (Bulirsch-Stoer)                                                  !
! -------------------------------------------------------------------------!

module module_sysdyn

use module_types, only : dp
use module_math, only : wiener_increment
implicit none (type, external)
private
public :: ode_model, rhs_interface
public :: lorenz_model, L63_rhs
public :: lotka_volterra, LV_rhs
public :: euler_step
public :: rk2_step, rk3_step, rk4_step, rk5_step
public :: srk4_step
public :: ab2_step, leapfrog_step, rkf45_step, bs_step

! ... We use polymorphism to allow multiple models and one interface
! ... We will create different dynamical systems that all work with the
! ... same integrator.
! ...
type, abstract :: ode_model
  character(len=20)                    :: name = 'ODE MODEL'
  integer                              :: ndim
  real(dp)                             :: noise_amplitude = 0.0D0 ! Default: deterministic
  character(len=20)                    :: noise_type = 'additive' ! 'additive' or 'multiplicative'
  real(dp), allocatable                :: f_prev(:)               ! Storage previoius derivative (AB2)
  contains
    procedure(rhs_interface), deferred :: rhs
    procedure                          :: run => ode_run
end type ode_model

! Interface for RHS function
abstract interface
  pure function rhs_interface(this,t,x) result(dxdt)
    import                       :: ode_model, dp
    class(ode_model), intent(in) :: this
    real(dp), intent(in)         :: t
    real(dp), intent(in)         :: x(:)
    real(dp)                     :: dxdt(size(x)) ! Returns array
  end function rhs_interface
end interface


! ================== Lorenz model ==================
type, extends(ode_model) :: lorenz_model
  real(dp)                       :: sigma = 10.0D0
  real(dp)                       :: rho   = 28.0D0
  real(dp)                       :: beta  = 8.0D0/3.0D0
  contains
    procedure                    :: rhs  => L63_rhs
end type lorenz_model

! ================= Lotka-Volterra =================
type, extends(ode_model) :: lotka_volterra
  real(dp)                       :: alpha = 1.0D0
  real(dp)                       :: beta  = 3.0D0
  real(dp)                       :: gamma = 2.0D0
  real(dp)                       :: delta = 0.3D0
  contains
    procedure                    :: rhs => LV_rhs
end type lotka_volterra

! ============= INTERFACES FOR CONSTRUCTORS ========  
interface lorenz_model  
  module procedure :: lorenz_constructor  
end interface  
  
interface lotka_volterra  
  module procedure :: lv_constructor  
end interface  

contains

! ============= CONSTRUCTORS =======================  
  
  function lorenz_constructor(noise_amplitude, noise_type) result(model)  
    real(dp), intent(in), optional         :: noise_amplitude
    character(len=*), intent(in), optional :: noise_type
    type(lorenz_model)                     :: model  
    model%name = 'LORENZ 1963 MODEL'
    model%ndim = 3  ! Set dimension  
    if (present(noise_amplitude)) model%noise_amplitude = noise_amplitude
    if (present(noise_type)) model%noise_type = noise_type
    ! Allocate storage for AB2 method  
    allocate(model%f_prev(model%ndim))  
    model%f_prev = 0.0D0  
  end function lorenz_constructor  
  
  function lv_constructor(noise_amplitude, noise_type) result(model)  
    real(dp), intent(in), optional         :: noise_amplitude
    character(len=*), intent(in), optional :: noise_type
    type(lotka_volterra)                   :: model  
    model%name = 'LOTKA-VOLTERRA PZ'
    model%ndim = 2  ! Set dimension  
    if (present(noise_amplitude)) model%noise_amplitude = noise_amplitude
    if (present(noise_type)) model%noise_type = noise_type
    ! Allocate storage for AB2 method  
    allocate(model%f_prev(model%ndim))  
    model%f_prev = 0.0D0  
  end function lv_constructor  

  ! ============= RHS FUNCTIONS ====================== 

  pure function L63_rhs(this,t,x) result(dxdt)
    class(lorenz_model), intent(in)     :: this
    real(dp), intent(in)                :: t
    real(dp), intent(in)                :: x(:)
    real(dp)                            :: dxdt(size(x))

    dxdt(1) = this%sigma*(x(2) - x(1))
    dxdt(2) = x(1)*(this%rho-x(3)) - x(2)
    dxdt(3) = x(1)*x(2) - this%beta*x(3)
  end function L63_rhs


  pure function LV_rhs(this,t,x) result(dxdt)
    class(lotka_volterra), intent(in)   :: this
    real(dp), intent(in)                :: t
    real(dp), intent(in)                :: x(:)
    real(dp)                            :: dxdt(size(x))

    dxdt(1) = this%alpha*(1.0D0-x(1)/this%beta)*x(1) - this%gamma*x(1)*x(2)
    dxdt(2) = this%gamma*x(1)*x(2) - this%delta*x(2)
  end function LV_rhs

  ! =============== EULER METHOD (for comparison) ====  
  
  function euler_step(model, t, x, dt) result(xn)  
    class(ode_model), intent(in) :: model  
    real(dp), intent(in)          :: t  
    real(dp), intent(in)          :: x(:)  
    real(dp), intent(in)          :: dt  
    real(dp)                      :: xn(size(x))  
      
    xn = x + dt * model%rhs(t, x)  
      
  end function euler_step 

  ! =============== Runge-Kutta Order 2 ==============

  function rk2_step(model,t,x,dt) result(xn)
  ! ... Runge-Kutta order 2 schema for time integration of the ODE: dx/dt = F
  ! ... Mid-point's method

    class(ode_model), intent(in)        :: model
    real(dp), intent(in)                :: t
    real(dp), intent(in)                :: x(:)
    real(dp), intent(in)                :: dt
    real(dp)                            :: xn(size(x))

    ! ... Local variables:
    real(dp), dimension(size(x))        :: k1,k2

    k1 = dt * model%rhs(t, x)
    k2 = dt * model%rhs(t+0.5D0*dt, x+0.5D0*k1)

    xn = x + k2 

  end function rk2_step

  ! =============== Runge-Kutta Order 3 ==============

  function rk3_step(model,t,x,dt) result(xn)
  ! ... Runge-Kutta order 3 schema for time integration of the ODE: dx/dt = F
  ! ... Heun's or Nystrom method

    class(ode_model), intent(in)        :: model
    real(dp), intent(in)                :: t
    real(dp), intent(in)                :: x(:)
    real(dp), intent(in)                :: dt
    real(dp)                            :: xn(size(x))

    ! ... Local variables:
    real(dp), dimension(size(x))        :: k1,k2,k3

    k1 = dt * model%rhs(t, x)
    k2 = dt * model%rhs(t+0.5D0*dt, x+0.5D0*k1)
    k3 = dt * model%rhs(t+dt, x-k1+2.0D0*k2)

    xn = x + (k1 + 4.0D0*k2 + k3) / 6.0D0 

  end function rk3_step

  ! =============== Runge-Kutta Order 4 ==============

  function rk4_step(model,t,x,dt) result(xn)
    class(ode_model), intent(in)        :: model
    real(dp), intent(in)                :: t
    real(dp), intent(in)                :: x(:)
    real(dp), intent(in)                :: dt
    real(dp)                            :: xn(size(x))

    ! ... Local variables:
    real(dp), dimension(size(x))        :: k1,k2,k3,k4

    k1 = dt * model%rhs(t, x)
    k2 = dt * model%rhs(t+0.5D0*dt, x+0.5D0*k1)
    k3 = dt * model%rhs(t+0.5D0*dt, x+0.5D0*k2)
    k4 = dt * model%rhs(t+dt, x+k3)

    xn = x + (k1 + 2.0D0*k2 + 2.0D0*k3 + k4) / 6.0D0 

  end function rk4_step

  ! =============== Runge-Kutta Order 5 ==============

  function rk5_step(model,t,x,dt) result(xn)
    class(ode_model), intent(in)        :: model
    real(dp), intent(in)                :: t
    real(dp), intent(in)                :: x(:)
    real(dp), intent(in)                :: dt
    real(dp)                            :: xn(size(x))

    ! ... Local variables:
    real(dp), dimension(size(x))        :: k1,k2,k3,k4,k5,k6

    k1 = dt * model%rhs(t, x)
    k2 = dt * model%rhs(t + dt * 0.2d0, x + k1 * 0.2d0)
    k3 = dt * model%rhs(t + dt * 0.3d0, x + k1 * 3.0d0 / 40.0d0 + k2 * 9.0d0 / 40.0d0)
    k4 = dt * model%rhs(t + dt * 0.6d0, x + k1 * 0.3d0 - k2 * 0.9d0 + k3 * 1.2d0)
    k5 = dt * model%rhs(t + dt, x - k1 * 11.0d0 / 54.0d0 + k2 * 2.5d0 - k3 * 70.0d0 / 27.0d0 + k4 * 35.0d0 / 27.0d0)
    k6 = dt * model%rhs(t + dt * 0.875d0, x + k1 * 1631.0d0 / 55296.0d0 + k2 * 175.0d0 / 512.0d0 + &
                        k3 * 575.0d0 / 13824.0d0 + k4 * 44275.0d0 / 110592.0d0 + k5 * 253.0d0 / 4096.0d0)

    xn = x + (k1 * 37.0d0 / 378.0d0 + k3 * 250.0d0 / 621.0d0 + k4 * 125.0d0 / 594.0d0 + &
            k6 * 512.0d0 / 1771.0d0)

  end function rk5_step

  ! ============ Stochastic Runge-Kutta ==============

  ! ==== INTERNAL DIFFUSION FUNCTION ====  
  ! Computes diffusion coefficients based on noise type  
  
  pure function compute_diffusion(model, x) result(gx)  
    class(ode_model), intent(in) :: model  
    real(dp), intent(in)          :: x(:)  
    real(dp)                      :: gx(size(x))  
    integer :: i  
    
    select case (trim(model%noise_type))  
      
      case ('additive', 'ADDITIVE')  
        ! Additive noise: g(x) = constant  
        gx = model%noise_amplitude  
      
      case ('multiplicative', 'MULTIPLICATIVE')  
        ! Multiplicative noise: g(x) = amplitude * |x|  
        do i = 1, size(x)  
          gx(i) = model%noise_amplitude * abs(x(i))  
        end do  
      
      case ('state-dependent', 'STATE-DEPENDENT', 'state_dependent')  
        ! State-dependent: g(x) = amplitude * sqrt(|x|)  
        do i = 1, size(x)  
          gx(i) = model%noise_amplitude * sqrt(abs(x(i)))  
        end do  
      
      case default  
        ! Default to additive  
        gx = model%noise_amplitude  
        
    end select  
    
  end function compute_diffusion

  function srk4_step(model, t, x, dt) result(xn)  
  ! For SDEs: dX = f(t,X)dt + g(X)dW  
  ! Handles additive and multiplicative noise internally  

    class(ode_model), intent(in) :: model  
    real(dp), intent(in)         :: t  
    real(dp), intent(in)         :: x(:)  
    real(dp), intent(in)         :: dt  
    real(dp)                     :: xn(size(x))  
    
    ! Local variables  
    real(dp), dimension(size(x)) :: k1, k2, k3, k4  
    real(dp), dimension(size(x)) :: g1, g2  
    real(dp), dimension(size(x)) :: x_temp  
    real(dp) :: dW1, dW2, sqrt_dt  
    
    ! Check if deterministic (no noise)  
    if (model%noise_amplitude == 0.0_dp) then  
      ! Fall back to deterministic RK4  
      xn = rk4_step(model, t, x, dt)  
      return  
    end if  
    
    sqrt_dt = sqrt(dt)  
    
    ! Generate two independent Wiener increments  
    dW1 = wiener_increment(dt)  
    dW2 = wiener_increment(dt)  
    
    ! Stage 1: Evaluate at current point  
    k1 = model%rhs(t, x)  
    g1 = compute_diffusion(model, x)  
    
    ! Stage 2: Evaluate at midpoint  
    x_temp = x + 0.5_dp*dt*k1 + 0.5_dp*sqrt_dt*g1  
    k2 = model%rhs(t + 0.5_dp*dt, x_temp)  
    g2 = compute_diffusion(model, x_temp)  
    
    ! Stage 3: Evaluate at midpoint (alternative)  
    x_temp = x + 0.5_dp*dt*k2 + 0.5_dp*sqrt_dt*g2  
    k3 = model%rhs(t + 0.5_dp*dt, x_temp)  
    
    ! Stage 4: Evaluate at endpoint  
    x_temp = x + dt*k3 + sqrt_dt*g2  
    k4 = model%rhs(t + dt, x_temp)  
    
    ! Combine stages  
    xn = x + (dt/6.0_dp) * (k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4) &  
           + sqrt_dt * g1 * dW1 &  
           + (g2 - g1) * (dW2 / sqrt(2.0_dp))  
    
  end function srk4_step  


 ! ==== Adams-Bashforth Order 2 ====  
  
  function ab2_step(model, t, x, dt, f_prev) result(xn)  
  ! ... Adams-Bashforth 2nd order multistep method  
  ! ... Requires previous derivative f_prev = f(t-dt, x(t-dt))  
  ! ... If f_prev not provided, falls back to Euler for first step  
  
    class(ode_model), intent(in)    :: model  
    real(dp), intent(in)    :: t  
    real(dp), intent(in)    :: x(:)  
    real(dp), intent(in)    :: dt  
    real(dp), intent(in), optional :: f_prev(:)  
    real(dp)    :: xn(size(x))  
  
    ! ... Local variables:  
    real(dp), dimension(size(x))    :: f_curr  
  
    f_curr = model%rhs(t, x)  
  
    if (present(f_prev)) then  
      ! AB2: x(n+1) = x(n) + dt/2 * (3*f(n) - f(n-1))  
      xn = x + dt * (1.5D0 * f_curr - 0.5D0 * f_prev)  
    else  
      ! First step: use Euler  
      xn = x + dt * f_curr  
    end if  
  
  end function ab2_step  


  ! ==== Leapfrog (Verlet) Method ====  
  
  function leapfrog_step(model, t, x, dt, x_prev) result(xn)  
  ! ... Leapfrog/Verlet symplectic integrator  
  ! ... Requires previous position x_prev = x(t-dt)  
  ! ... If x_prev not provided, falls back to RK2 for first step  
  ! ... Best for Hamiltonian systems (energy-conserving)  
  
    class(ode_model), intent(in)    :: model  
    real(dp), intent(in)    :: t  
    real(dp), intent(in)    :: x(:)  
    real(dp), intent(in)    :: dt  
    real(dp), intent(in), optional :: x_prev(:)  
    real(dp)    :: xn(size(x))  
  
    if (present(x_prev)) then  
      ! Leapfrog: x(n+1) = 2*x(n) - x(n-1) + dt^2 * f(n)  
      xn = 2.0D0 * x - x_prev + dt * dt * model%rhs(t, x)  
    else  
      ! First step: use RK2  
      xn = rk2_step(model, t, x, dt)  
    end if  
  
  end function leapfrog_step  
 

  ! ==== Runge-Kutta-Fehlberg 4(5) ====  
  
  function rkf45_step(model, t, x, dt, error_est) result(xn)  
  ! ... Runge-Kutta-Fehlberg embedded method  
  ! ... Provides 4th order solution with 5th order error estimate  
  ! ... Can be used for adaptive step size control  
  
    class(ode_model), intent(in)    :: model  
    real(dp), intent(in)    :: t  
    real(dp), intent(in)    :: x(:)  
    real(dp), intent(in)    :: dt  
    real(dp), intent(out), optional :: error_est(:)  
    real(dp)    :: xn(size(x))  
  
    ! ... Local variables:  
    real(dp), dimension(size(x))    :: k1, k2, k3, k4, k5, k6  
    real(dp), dimension(size(x))    :: x4, x5  
  
    ! Butcher tableau for RKF45  
    k1 = dt * model%rhs(t, x)  
    k2 = dt * model%rhs(t + dt/4.0D0, x + k1/4.0D0)  
    k3 = dt * model%rhs(t + dt*3.0D0/8.0D0, x + k1*3.0D0/32.0D0 + k2*9.0D0/32.0D0)  
    k4 = dt * model%rhs(t + dt*12.0D0/13.0D0, &  
         x + k1*1932.0D0/2197.0D0 - k2*7200.0D0/2197.0D0 + k3*7296.0D0/2197.0D0)  
    k5 = dt * model%rhs(t + dt, &  
         x + k1*439.0D0/216.0D0 - 8.0D0*k2 + k3*3680.0D0/513.0D0 - k4*845.0D0/4104.0D0)  
    k6 = dt * model%rhs(t + dt/2.0D0, &  
         x - k1*8.0D0/27.0D0 + 2.0D0*k2 - k3*3544.0D0/2565.0D0 + k4*1859.0D0/4104.0D0 - k5*11.0D0/40.0D0)  
  
    ! 4th order solution  
    x4 = x + k1*25.0D0/216.0D0 + k3*1408.0D0/2565.0D0 + k4*2197.0D0/4104.0D0 - k5/5.0D0  
  
    ! 5th order solution  
    x5 = x + k1*16.0D0/135.0D0 + k3*6656.0D0/12825.0D0 + k4*28561.0D0/56430.0D0 - k5*9.0D0/50.0D0 + k6*2.0D0/55.0D0  
  
    ! Return 4th order solution  
    xn = x4  
  
    ! Optionally return error estimate  
    if (present(error_est)) then  
      error_est = abs(x5 - x4)  
    end if  
  
  end function rkf45_step  
  
  ! ==== Bulirsch-Stoer Method ====  
  
  function bs_step(model, t, x, dt, max_order) result(xn)  
  ! ... Bulirsch-Stoer extrapolation method  
  ! ... Uses Richardson extrapolation with modified midpoint method  
  ! ... High accuracy for smooth problems  
  
    class(ode_model), intent(in)    :: model  
    real(dp), intent(in)    :: t  
    real(dp), intent(in)    :: x(:)  
    real(dp), intent(in)    :: dt  
    integer, intent(in), optional   :: max_order  
    real(dp)    :: xn(size(x))  
  
    ! ... Local variables:  
    integer :: n, k, max_k  
    integer, parameter :: nseq(8) = [2, 4, 6, 8, 12, 16, 24, 32]  
    real(dp), dimension(size(x), 8) :: W  
    real(dp), dimension(size(x))    :: x_temp, f_temp  
  
    ! Set maximum order (default: 4)  
    max_k = 4  
    if (present(max_order)) max_k = min(max_order, 8)  
  
    ! Compute sequence of approximations with different substeps  
    do k = 1, max_k  
      n = nseq(k)  
      W(:, k) = modified_midpoint(model, t, x, dt, n)  
    end do  
  
    ! Richardson extrapolation  
    do k = 2, max_k  
      do n = k, max_k  
        W(:, n) = W(:, n) + (W(:, n) - W(:, n-1)) / &  
                  (real(nseq(n), dp)**2 / real(nseq(n-k+1), dp)**2 - 1.0D0)  
      end do  
    end do  
  
    xn = W(:, max_k)  
  
  contains  
  
    function modified_midpoint(model, t, x, h, n) result(xn)  
    ! Modified midpoint method for Bulirsch-Stoer  
      class(ode_model), intent(in) :: model  
      real(dp), intent(in)    :: t, x(:), h  
      integer, intent(in)     :: n  
      real(dp)    :: xn(size(x))  
      real(dp)    :: dt_sub, t_sub  
      real(dp), dimension(size(x)) :: x0, x1, x2  
      integer :: i  
  
      dt_sub = h / real(n, dp)  
        
      ! First step  
      x0 = x  
      x1 = x0 + dt_sub * model%rhs(t, x0)  
        
      ! Midpoint steps  
      do i = 1, n-1  
        t_sub = t + real(i, dp) * dt_sub  
        x2 = x0 + 2.0D0 * dt_sub * model%rhs(t_sub, x1)  
        x0 = x1  
        x1 = x2  
      end do  
        
      ! Final smoothing step  
      xn = 0.5D0 * (x1 + x0 + dt_sub * model%rhs(t + h, x1))  
        
    end function modified_midpoint  
  
  end function bs_step  
   


  ! =============== MODEL RUN METHOD =================  
  
  subroutine ode_run(this, t0, x0, dt, nsteps, t, X, method, output_file, form, time_step_output, verbose)  
    class(ode_model), intent(inout)        :: this  
    real(dp), intent(in)                   :: t0      ! Initial time  
    real(dp), intent(in)                   :: x0(:)   ! Initial state  
    real(dp), intent(in)                   :: dt      ! Time step  
    integer, intent(in)                    :: nsteps  ! Number of steps  
    real(dp), intent(out)                  :: t(nsteps+1)      ! Time array  
    real(dp), intent(out)                  :: X(size(x0), nsteps+1)  ! State array  
    character(len=*), intent(in), optional :: method  ! Integration method  
    character(len=*), intent(in), optional :: output_file  ! Optional output file name
    character(len=*), intent(in), optional :: form  ! Optional output file format
    integer, intent(in), optional          :: time_step_output ! Write evey time_step_output steps
    logical, intent(in), optional          :: verbose
      
    ! Local variables  
    logical :: write_output, write_ascii, show_progress
    character(len=20) :: integrator  
    integer :: i,iu,output_every
    real(dp) :: tcur  
    real(dp), dimension(size(x0)) :: xcur, xnew  
      
    ! ... Validate array sizes  
    ! ...
    if (size(x0) /= this%ndim) then  
      error stop 'ode_run: Initial state size does not match model dimension'  
    end if  
      
    ! ... Set integration method (default: rk4)  
    ! ...
    if (present(method)) then  
      integrator = method  
    else  
      integrator = 'rk4'  
    end if  
     
    ! ... Check if output file is requested  
    ! ...
    write_output = present(output_file)   

    ! ... Set output format (default: ASCII)
    write_ascii = .true.  
    if (present(form)) then  
      select case (trim(form))  
        case ('bin', 'BIN', 'binary', 'BINARY', 'unformatted', 'UNFORMATTED')  
          write_ascii = .false.  
        case default
          write_ascii = .true.  
      end select  
    end if  

    ! ... Open output file if requested  
    ! ...
    if (write_output) then  
      if (write_ascii) then
        open(newunit=iu, file=output_file, status='replace', action='write')  
        ! Write header 
        write(iu, '(A)', advance='no') '# time'  
        do i = 1, this%ndim  
          write(iu, '(A,I0)', advance='no') '  x', i  
        end do  
        write(iu, *)  ! New line  
      else
        open(newunit=iu, file=output_file, status='replace', &
             form='unformatted', action='write')  
      endif
    end if  

    ! ... Writes to file every N steps instead of every step  
    ! ... Reduces file size for long simulations
    ! ...
    output_every = 1
    if (present(time_step_output)) then
      output_every = time_step_output
      if (mod(nsteps,time_step_output).ne.0) write(*,*) 'Warning: last point not saved in file'
    endif

    show_progress = .False.
    if (present(verbose)) show_progress = verbose

    ! ... Initialize  
    tcur = t0  
    xcur = x0  
      
    ! ... Store initial condition  
    t(1) = tcur  
    X(:, 1) = xcur  
    
    ! ... Write initial condition to file  
    if (write_output) call write_step()
 
    ! Time integration loop  
    do i = 1, nsteps  
        
      ! Select integration method  
      select case (trim(integrator))  

        case ('euler', 'EULER')  
          xnew = euler_step(this, tcur, xcur, dt)  

        case ('rk2', 'RK2')  
          xnew = rk2_step(this, tcur, xcur, dt)  
          
        case ('rk3', 'RK3')  
          xnew = rk3_step(this, tcur, xcur, dt)  
          
        case ('rk4', 'RK4')  
          xnew = rk4_step(this, tcur, xcur, dt)  
          
        case ('rk5', 'RK5')  
          xnew = rk5_step(this, tcur, xcur, dt)  
          
        case ('srk4', 'SRK4')
          xnew = srk4_step(this, tcur, xcur, dt)

        case ('ab2', 'AB2')  
          if (i > 1) then  
            xnew = ab2_step(this, tcur, xcur, dt, this%f_prev)  
          else  
            xnew = ab2_step(this, tcur, xcur, dt)  
          end if  
      
        case ('leapfrog', 'LEAPFROG')  
          if (i > 1) then  
            xnew = leapfrog_step(this, tcur, xcur, dt, X(:, i-1))  
          else  
            xnew = leapfrog_step(this, tcur, xcur, dt)  
          end if  
      
        case ('rkf45', 'RKF45')  
          xnew = rkf45_step(this, tcur, xcur, dt)  
      
        case ('bs', 'BS')  
          xnew = bs_step(this, tcur, xcur, dt)            

        case default    
          stop 'ode_run: Unknown integration method: ' // trim(integrator)    

      end select  
        
      ! ... Update state  
      ! ...
      tcur = tcur + dt  
      xcur = xnew  
        
      ! ... Store results  
      ! ...
      t(i+1)   = tcur  
      X(:,i+1) = xcur  
      this%f_prev = this%rhs(tcur, xcur)  ! Used by AB2


      ! ... Write to file if requested  
      ! ...
      if (write_output.and.mod(i,output_every).eq.0) call write_step()

      ! ... Progress report
      ! ...
      if (show_progress.and.nsteps.ge.10) then
        if (mod(i,nsteps/10).eq.0) write(*,'(T2,A,F6.1,A)') 'Progress: ', 100.0*real(i)/real(nsteps), '%'
      endif
    end do  

    ! ... Close output file  
    ! ...
    if (write_output) then  
      close(iu)  
      print '(T2,A,A)', 'Output written to: ', trim(output_file)  
    end if  
     
    contains 
      subroutine write_step()     ! It inherits all the variables
        if (write_ascii) then
          write(iu, '(*(ES16.7))') tcur, xcur  
        else
          write(iu) tcur, xcur  
        endif
      end subroutine write_step 

  end subroutine ode_run  

end module module_sysdyn


