module module_ar
  !! Autoregressive models AR(p) for p >= 0.
  !!
  !! Primary estimator: Yule–Walker solved by Levinson–Durbin recursion (Toeplitz).
  !! This is a standard, fast classical approach; for many practical baselines it is
  !! considered “state of the art” in terms of speed/stability for AR-only models.
  !!
  !! Model:
  !!   y_t = c + sum_{i=1..p} phi_i * y_{t-i} + e_t
  !!
  !! Notes on intercept:
  !! - If include_intercept=.true., we estimate mean mu and fit AR to centered data:
  !!     x_t = y_t - mu
  !!     x_t = sum phi_i x_{t-i} + e_t
  !!   then convert to intercept form using: c = mu * (1 - sum phi_i)
  !!
  use module_types, only : dp
  implicit none
  private

  integer, parameter, public :: AR_METHOD_YW_LD = 1

  type, public :: type_ar_model
     integer :: p = 0
     logical        :: include_intercept = .true.
     integer :: method = AR_METHOD_YW_LD

     real(dp)   :: mu = 0.0_dp           !! sample mean (if intercept used)
     real(dp)   :: c  = 0.0_dp           !! intercept in original y-space
     real(dp), allocatable :: phi(:)         !! AR coefficients, size p

     real(dp)   :: sigma2 = 0.0_dp       !! innovation variance estimate
     real(dp)   :: rmse   = 0.0_dp
     real(dp)   :: aic    = 0.0_dp
     real(dp)   :: bic    = 0.0_dp
     integer :: nobs   = 0                !! effective obs used for residuals
  contains
     procedure :: fit => ar_fit
     procedure :: one_step_pred => ar_one_step_pred
     procedure :: forecast => ar_forecast
  end type type_ar_model

contains

  subroutine ar_fit(this, y, p, include_intercept, method)
    class(type_ar_model), intent(inout) :: this
    real(dp),          intent(in)   :: y(:)
    integer,        intent(in)   :: p
    logical,               intent(in), optional :: include_intercept
    integer,        intent(in), optional :: method

    integer :: n, n_eff, npar
    real(dp), allocatable :: x(:), r(:), resid(:)
    real(dp) :: sse, ll, sumphi

    n = size(y)
    if (p < 0) error stop "ar_fit: p must be >= 0"
    if (n < max(3, p+2)) error stop "ar_fit: not enough data"

    this%p = p
    if (present(include_intercept)) this%include_intercept = include_intercept
    if (present(method)) this%method = method

    if (allocated(this%phi)) deallocate(this%phi)
    allocate(this%phi(p))
    this%phi = 0.0_dp

    ! Centering if intercept requested
    allocate(x(n))
    if (this%include_intercept) then
       this%mu = sum(y) / real(n, dp)
       x = y - this%mu
    else
       this%mu = 0.0_dp
       x = y
    end if

    if (p == 0) then
       ! AR(0): y_t = c + e_t. With intercept, c=mu; without, c=0.
       this%phi = 0.0_dp
       if (this%include_intercept) then
          this%c = this%mu
       else
          this%c = 0.0_dp
       end if

       ! residuals are y - c
       n_eff = n
       allocate(resid(n_eff))
       resid = y - this%c
       sse = sum(resid*resid)
       this%rmse = sqrt(sse / real(n_eff, dp))
       this%sigma2 = sse / real(n_eff, dp)

       ll = -0.5_dp * real(n_eff, dp) * ( log(2.0_dp*acos(-1.0_dp)) + log(this%sigma2) + 1.0_dp )
       npar = 1 ! sigma2
       if (this%include_intercept) npar = npar + 1 ! c
       this%aic = -2.0_dp*ll + 2.0_dp*real(npar, dp)
       this%bic = -2.0_dp*ll + log(real(n_eff, dp))*real(npar, dp)
       this%nobs = n_eff

       deallocate(resid, x)
       return
    end if

    select case (this%method)
    case (AR_METHOD_YW_LD)
       ! Autocovariances r(0..p)
       allocate(r(0:p))
       call autocovariances(x, p, r)

       ! Solve Toeplitz system for phi and get innovation variance sigma2
       call levinson_durbin(r, this%phi, this%sigma2)
       deallocate(r)
    case default
       error stop "ar_fit: unknown method"
    end select

    ! Convert centered AR to intercept form if needed:
    ! y_t = c + sum phi_i y_{t-i} + e_t, where c = mu*(1 - sum phi_i)
    sumphi = sum(this%phi)
    if (this%include_intercept) then
       this%c = this%mu * (1.0_dp - sumphi)
    else
       this%c = 0.0_dp
    end if

    ! Compute conditional residuals on original y-space (fair for benchmarking)
    n_eff = n - p
    this%nobs = n_eff
    allocate(resid(n_eff))
    call compute_residuals(y, this%c, this%phi, resid)

    sse = sum(resid*resid)
    this%rmse = sqrt(sse / real(n_eff, dp))

    ! Keep sigma2 consistent with likelihood; you can also overwrite with sse/n_eff.
    ! For YW, sigma2 from LD recursion corresponds to innovation variance on centered series.
    ! We’ll use sse/n_eff as conditional variance for scoring/AIC/BIC.
    this%sigma2 = sse / real(n_eff, dp)

    ll = -0.5_dp * real(n_eff, dp) * ( log(2.0_dp*acos(-1.0_dp)) + log(this%sigma2) + 1.0_dp )

    npar = 1          ! sigma2
    npar = npar + p   ! phi_1..phi_p
    if (this%include_intercept) npar = npar + 1 ! c
    this%aic = -2.0_dp*ll + 2.0_dp*real(npar, dp)
    this%bic = -2.0_dp*ll + log(real(n_eff, dp))*real(npar, dp)

    deallocate(resid, x)
  end subroutine ar_fit

  function ar_one_step_pred(this, y) result(yhat)
    class(type_ar_model), intent(in) :: this
    real(dp),          intent(in) :: y(:)
    real(dp) :: yhat(size(y))
    integer :: n, t, p, i

    n = size(y)
    p = this%p
    yhat = 0.0_dp

    if (p == 0) then
       yhat = this%c
       return
    end if

    do t = p+1, n
       yhat(t) = this%c
       do i = 1, p
          yhat(t) = yhat(t) + this%phi(i) * y(t-i)
       end do
    end do
  end function ar_one_step_pred

  function ar_forecast(this, y_hist, h) result(f)
    class(type_ar_model), intent(in) :: this
    real(dp),          intent(in) :: y_hist(:)  ! need last p values
    integer,        intent(in) :: h
    real(dp) :: f(h)

    integer :: p, i, j
    real(dp), allocatable :: buf(:)  ! holds last p values, updated with forecasts
    real(dp) :: next

    if (h < 1) error stop "ar_forecast: h must be >= 1"
    p = this%p

    if (p == 0) then
       f = this%c
       return
    end if
    if (size(y_hist) < p) error stop "ar_forecast: need at least p history points"

    allocate(buf(p))
    buf = y_hist(size(y_hist)-p+1 : size(y_hist))

    do i = 1, h
       next = this%c
       do j = 1, p
          next = next + this%phi(j) * buf(p-j+1)  ! buf end is most recent
       end do
       f(i) = next

       ! shift left and append next
       if (p > 1) buf(1:p-1) = buf(2:p)
       buf(p) = next
    end do

    deallocate(buf)
  end function ar_forecast

  subroutine autocovariances(x, p, r)
    real(dp), intent(in)  :: x(:)
    integer, intent(in) :: p
    real(dp), intent(out) :: r(0:p)

    integer :: n, k
    real(dp) :: denom

    n = size(x)
    ! Biased autocovariance estimate (divide by n). Common in YW.
    denom = real(n, dp)

    do k = 0, p
       r(k) = sum( x(1:n-k) * x(1+k:n) ) / denom
    end do
  end subroutine autocovariances

  subroutine levinson_durbin(r, a, sigma2)
    !! Levinson–Durbin recursion to solve Toeplitz system for AR coefficients.
    !!
    !! Inputs:
    !!   r(0:p): autocovariances
    !! Outputs:
    !!   a(1:p): AR coefficients
    !!   sigma2: innovation variance
    real(dp), intent(in)  :: r(0:)
    real(dp), intent(out) :: a(:)
    real(dp), intent(out) :: sigma2

    integer :: p, m, i
    real(dp), allocatable :: a_prev(:)
    real(dp) :: kappa, e, acc

    p = size(a)
    if (size(r) < p+1) error stop "levinson_durbin: r too short"

    if (r(0) <= 0.0_dp) error stop "levinson_durbin: r(0) must be > 0"

    allocate(a_prev(p))
    a = 0.0_dp
    a_prev = 0.0_dp

    e = r(0)

    do m = 1, p
       ! reflection coefficient kappa
       acc = r(m)
       do i = 1, m-1
          acc = acc - a_prev(i) * r(m-i)
       end do
       kappa = acc / e

       ! update coefficients
       a(1:m-1) = a_prev(1:m-1) - kappa * a_prev(m-1:1:-1)
       a(m) = kappa

       ! update error
       e = e * (1.0_dp - kappa*kappa)
       if (e <= 0.0_dp) then
          ! numerical issues / non-PD Toeplitz; clamp
          e = max(e, 1.0e-30_dp)
       end if

       a_prev = a
    end do

    sigma2 = e
    deallocate(a_prev)
  end subroutine levinson_durbin

  subroutine compute_residuals(y, c, phi, resid)
    real(dp), intent(in)  :: y(:)
    real(dp), intent(in)  :: c
    real(dp), intent(in)  :: phi(:)
    real(dp), intent(out) :: resid(:)

    integer :: n, p, t, i, idx
    real(dp) :: yhat

    n = size(y)
    p = size(phi)
    if (size(resid) /= n - p) error stop "compute_residuals: bad resid size"

    do t = p+1, n
       yhat = c
       do i = 1, p
          yhat = yhat + phi(i) * y(t-i)
       end do
       idx = t - p
       resid(idx) = y(t) - yhat
    end do
  end subroutine compute_residuals

end module module_ar
