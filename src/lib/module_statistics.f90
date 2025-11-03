!> @file module_statistics.f90
!> @brief Statistical analysis routines: descriptive statistics, moments, percentiles, correlation.
!> @author Quim Ballabrera, Institut de Ciencies del Mar, CSIC
!> @date April 2022, refactored 2025
!>
!> This module provides comprehensive statistical functions for data analysis.
!> For pure math utilities (RNG, sorting, optimization), see module_math.
!>
!> Copyright (C) 2022-2025, Joaquim Ballabrera
!> Licensed under GNU LGPL v3+

module module_statistics

use iso_fortran_env, only: output_unit
use module_types, only: dp
use module_constants, only: nan
use module_tools, only: crash, indexx, quicksort
use module_math, only: t_inverse_cdf, t_cdf_complement, f_cdf_complement

implicit none (type, external)

private

public :: mean, var, std
public :: median, mad, percentile, quantile
public :: moment, skewness, kurtosis
public :: corr, cov
public :: zscore, trimmed_mean
public :: stats, get_stats
public :: spearman, linreg
public :: rms

!> @brief Derived type for comprehensive statistical summary.
!!
!! All fields have units matching the input data, except:
!! - var: units^2
!! - skew, kurt: dimensionless
!! - Percentiles and fences: same units as data
type stats
  integer :: n                      !< Number of data points
  real(dp) :: min                   !< Minimum value
  real(dp) :: max                   !< Maximum value
  real(dp) :: median                !< Median (50th percentile)
  real(dp) :: mean                  !< Arithmetic mean
  real(dp) :: trimmed_mean          !< Trimmed mean (5% each tail)
  real(dp) :: adev                  !< Average absolute deviation from mean
  real(dp) :: iqr                   !< Interquartile range (p75 - p25)
  real(dp) :: std                   !< Standard deviation (unbiased)
  real(dp) :: var                   !< Variance (unbiased)
  real(dp) :: rms                   !< Root mean square
  real(dp) :: skew                  !< Skewness (dimensionless)
  real(dp) :: kurt                  !< Excess kurtosis (dimensionless)
  real(dp) :: p01, p05, p25, p75, p95, p99  !< Percentiles
  real(dp) :: lower_inner_fence     !< p25 - 1.5*IQR
  real(dp) :: lower_outer_fence     !< p25 - 3.0*IQR
  real(dp) :: upper_inner_fence     !< p75 + 1.5*IQR
  real(dp) :: upper_outer_fence     !< p75 + 3.0*IQR
  contains
    procedure :: show => stats_show
end type stats

!> @brief Derived type for t-test results.
!!
!! Contains test statistic, p-value, degrees of freedom, and confidence interval.
type, public :: ttest_result
  real(dp) :: t_stat        !< t-statistic
  real(dp) :: p_value       !< two-tailed p-value
  real(dp) :: df            !< degrees of freedom
  real(dp) :: ci_lower      !< lower bound of 95% CI for mean difference
  real(dp) :: ci_upper      !< upper bound of 95% CI for mean difference
  real(dp) :: mean_diff     !< mean difference (sample mean - mu0, or mean1 - mean2)
  real(dp) :: std_err       !< standard error of the mean difference
end type ttest_result

!> @brief Derived type to hold a group of data for ANOVA.  
type :: data_group  
  real(dp), dimension(:), allocatable :: data  
end type data_group  

!> @brief Derived type for ANOVA results.
!!
!! Contains F-statistic, p-value, degrees of freedom, and sum of squares.
type, public :: anova_result
  real(dp) :: f_stat        !< F-statistic
  real(dp) :: p_value       !< p-value
  integer :: df_between     !< degrees of freedom between groups
  integer :: df_within      !< degrees of freedom within groups
  real(dp) :: ss_between    !< sum of squares between groups
  real(dp) :: ss_within     !< sum of squares within groups
  real(dp) :: ms_between    !< mean square between groups
  real(dp) :: ms_within     !< mean square within groups
end type anova_result

contains

  !> @brief Compute arithmetic mean of array A with optional weights W.
  !!
  !! Units: same as A.
  !!
  !! @param[in] A Real(dp) array, data values.
  !! @param[in] W Real(dp) array, optional. Weights (same size as A). If absent, uniform weights assumed.
  !! @return Real(dp), weighted mean = sum(W*A)/sum(W). Returns NaN if size(A)=0.
  !!
  !! Notes:
  !! - If W is present, must have same size as A.
  !! - Negative weights are allowed but may produce unexpected results.
  real(dp) function mean(A, W)
    real(dp), dimension(:), intent(in) :: A
    real(dp), dimension(:), optional :: W
    integer :: n
    real(dp) :: Sw

    mean = nan()
    n = size(A)
    if (n == 0) return

    if (present(W)) then
      Sw = sum(W)
      mean = dot_product(W, A)/Sw
    else
      mean = sum(A)/n
    end if
  end function mean


  !> @brief Compute variance of array A with optional weights and bias correction.
  !!
  !! Units: (units of A)^2.
  !!
  !! @param[in] A Real(dp) array, data values.
  !! @param[in] W Real(dp) array, optional. Weights (same size as A). If absent, uniform weights.
  !! @param[in] biased Logical, optional. If true, use biased estimator (divide by N). Default: false (unbiased, divide by N-1).
  !! @return Real(dp), variance. Returns NaN if size(A)=0.
  !!
  !! Notes:
  !! - Unbiased variance: sum((A - mean(A))^2)/(N-1).
  !! - Biased variance: sum((A - mean(A))^2)/N.
  !! - For weighted variance, Bessel correction is applied: Sw/(Sw - 1).
  real(dp) function var(A, W, biased)
    real(dp), dimension(:), intent(in) :: A
    real(dp), dimension(:), optional :: W
    logical, intent(in), optional :: biased
    logical :: weight, isbiased
    integer :: N, i
    real(dp) :: xsum1, xsum2, ai, wi, Sw

    var = nan()
    N = size(A)
    if (N == 0) return

    weight = .false.
    if (present(W)) then
      if (size(W) /= N) call crash('In var - Incompatible weight size')
      weight = .true.
    end if

    isbiased = .false.
    if (present(biased)) isbiased = biased

    if (weight) then
      Sw = 0.0_dp
      xsum1 = 0.0_dp
      xsum2 = 0.0_dp
      do i = 1, N
        wi = W(i)
        ai = A(i)
        Sw = Sw + wi
        xsum1 = xsum1 + wi*ai
        xsum2 = xsum2 + wi*ai*ai
      end do
      xsum1 = xsum1 / Sw
      xsum2 = xsum2/Sw - xsum1*xsum1
      var = Sw * xsum2 / (Sw - 1.0_dp)
    else
      xsum1 = 0.0_dp
      xsum2 = 0.0_dp
      do i = 1, N
        ai = A(i)
        xsum1 = xsum1 + ai
        xsum2 = xsum2 + ai*ai
      end do
      xsum1 = xsum1 / N
      xsum2 = xsum2 / N - xsum1*xsum1
      if (isbiased) then
        var = xsum2
      else
        var = N*xsum2/(N - 1.0_dp)
      end if
    end if
    if (var < 0.0_dp) var = 0.0_dp
  end function var


  !> @brief Compute standard deviation (square root of variance).
  !!
  !! Units: same as A.
  !!
  !! @param[in] A Real(dp) array, data values.
  !! @param[in] W Real(dp) array, optional. Weights.
  !! @param[in] biased Logical, optional. If true, use biased variance. Default: false.
  !! @return Real(dp), standard deviation = sqrt(variance(A, W, biased)).
  real(dp) function std(A, W, biased)
    real(dp), dimension(:), intent(in) :: A
    real(dp), dimension(:), optional :: W
    logical, intent(in), optional :: biased
    std = sqrt(var(A, W, biased))
  end function std


  !> @brief Compute root mean square (rms).
  !!
  !! Units: same as A.
  !!
  !! @param[in] A Real(dp) array, data values.
  !! @return Real(dp), rms = sqrt(mean(a*a)).
  real(dp) function rms(A)
    real(dp), dimension(:), intent(in) :: A
    rms = sqrt(dot_product(A,A)/size(A))
  end function rms


  !> @brief Compute median (50th percentile) of array x.
  !!
  !! Units: same as x.
  !!
  !! @param[in] x Real(dp) array, data values.
  !! @return Real(dp), median value. For even N, returns average of two middle values.
  !!
  !! Notes:
  !! - Sorts a copy of x internally; does not modify x.
  function median(x) result(med)
    real(dp), dimension(:), intent(in) :: x
    real(dp) :: med
    real(dp), dimension(size(x)) :: x_sorted
    integer :: n

    n = size(x)
    x_sorted = x
    call quicksort(x_sorted)

    if (mod(n, 2) == 0) then
      med = (x_sorted(n/2) + x_sorted(n/2+1)) / 2.0_dp
    else
      med = x_sorted((n+1)/2)
    end if
  end function median

  !> @brief Local MAD: median(|x - median(x)|).
  real(dp) function mad(v)
    real(dp), intent(in) :: v(:)
    real(dp), allocatable :: w(:)
    real(dp) :: medx
    integer :: n, i

    n = size(v)
    if (n <= 0) call crash('mad - empty vector')

    medx = median(v)
    allocate(w(n))
    do i = 1, n
      w(i) = abs(v(i) - medx)
    end do
    mad = median(w)
    deallocate(w)
  end function mad


  !> @brief Compute the pth percentile of array x using linear interpolation (Type-7, Excel/NumPy default).
  !!
  !! Units: same as x.
  !!
  !! @param[in] x Real(dp) array, data values.
  !! @param[in] p Real(dp), percentile in [0, 100].
  !! @param[in] method Character, optional. 'excel' (default) or 'nist'. Case-insensitive.
  !! @return Real(dp), percentile value. Returns NaN if p out of range.
  !!
  !! Notes:
  !! - Excel method (Type-7): rank = p*(N-1)/100 + 1, linear interpolation.
  !! - NIST method: rank = p*(N+1)/100, linear interpolation.
  !! - For p=0, returns min(x); for p=100, returns max(x).
  !! - Sorts a copy of x internally; does not modify x.
  function percentile(x, p, method) result(perc)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(in) :: p
    character(len=*), optional :: method
    real(dp) :: perc
    logical :: excel, nist
    integer :: n, kk
    integer, dimension(size(x)) :: indx
    real(dp) :: rr, dd

    if (p < 0.0_dp .or. p > 100.0_dp) then
      perc = nan()
      return
    end if

    excel = .true.
    if (present(method)) then
      if ((method(1:1) == 'N') .or. (method(1:1) == 'n')) then
        nist = .true.
        excel = .false.
      end if
    end if

    n = size(x)
    call indexx(x, indx)

    if (excel) then
      rr = p*(n - 1.0_dp)/100.0_dp + 1.0_dp
    else
      rr = p*(n + 1.0_dp)/100.0_dp
    end if

    kk = floor(rr)
    dd = rr - kk

    if (kk == 0) then
      perc = x(indx(1))
    else if (kk == n) then
      perc = x(indx(n))
    else
      perc = x(indx(kk)) + dd*(x(indx(kk+1)) - x(indx(kk)))
    end if
  end function percentile

  !> @brief Local quantile wrapper using percentile(p in [0,100]).
  real(dp) function quantile(v, q)
    real(dp), intent(in) :: v(:), q
    if (q < 0.0_dp .or. q > 1.0_dp) call crash('quantile - q must be in [0,1]')
    quantile = percentile(v, 100.0_dp*q)
  end function quantile

  !> @brief Compute central moments up to 4th order, plus mean and average absolute deviation.
  !!
  !! Based on Numerical Recipes `moment` subroutine.
  !!
  !! @param[in] data Real(dp) array, data values.
  !! @param[out] ave Real(dp), mean. Units: same as data.
  !! @param[out] adev Real(dp), average absolute deviation from mean. Units: same as data.
  !! @param[out] sdev Real(dp), standard deviation (unbiased). Units: same as data.
  !! @param[out] var Real(dp), variance (unbiased). Units: (data units)^2.
  !! @param[out] skew Real(dp), skewness (dimensionless).
  !! @param[out] curt Real(dp), excess kurtosis (dimensionless, normal = 0).
  !!
  !! Notes:
  !! - Requires n >= 2.
  !! - If variance = 0, skew and curt are set to 0.
  subroutine moment(data, ave, adev, sdev, var, skew, curt)
    real(dp), dimension(:), intent(in) :: data
    real(dp), intent(out) :: ave, adev, sdev, var, skew, curt
    integer :: n, j
    real(dp) :: p, s, ep

    n = size(data)
    if (n <= 1) call crash('in moment - n must be at least 2 in moment')

    s = 0.0_dp
    do j = 1, n
      s = s + data(j)
    end do
    ave = s/n

    adev = 0.0_dp
    var = 0.0_dp
    skew = 0.0_dp
    curt = 0.0_dp
    ep = 0.0_dp

    do j = 1, n
      s = data(j) - ave
      ep = ep + s
      adev = adev + abs(s)
      p = s*s
      var = var + p
      p = p*s
      skew = skew + p
      p = p*s
      curt = curt + p
    end do

    adev = adev/n
    var = (var - ep**2/n)/(n - 1)
    sdev = sqrt(var)

    if (var /= 0.0_dp) then
      skew = skew/(n*sdev**3)
      curt = curt/(n*var**2) - 3.0_dp
    else
      skew = 0.0_dp
      curt = 0.0_dp
    end if
  end subroutine moment


  !> @brief Compute skewness (3rd standardized moment) of array x.
  !!
  !! Units: dimensionless.
  !!
  !! @param[in] x Real(dp) array, data values.
  !! @return Real(dp), skewness. Returns 0 if variance = 0.
  !!
  !! Notes:
  !! - Positive skew: right tail longer.
  !! - Negative skew: left tail longer.
  function skewness(x) result(skew)
    real(dp), dimension(:), intent(in) :: x
    real(dp) :: skew, ave, adev, sdev, var, curt
    call moment(x, ave, adev, sdev, var, skew, curt)
  end function skewness


  !> @brief Compute excess kurtosis (4th standardized moment - 3) of array x.
  !!
  !! Units: dimensionless.
  !!
  !! @param[in] x Real(dp) array, data values.
  !! @return Real(dp), excess kurtosis. Normal distribution has kurtosis = 0. Returns 0 if variance = 0.
  !!
  !! Notes:
  !! - Positive kurtosis: heavy tails (leptokurtic).
  !! - Negative kurtosis: light tails (platykurtic).
  function kurtosis(x) result(curt)
    real(dp), dimension(:), intent(in) :: x
    real(dp) :: curt, ave, adev, sdev, var, skew
    call moment(x, ave, adev, sdev, var, skew, curt)
  end function kurtosis


  !> @brief Compute Pearson correlation coefficient between arrays x and y.
  !!
  !! Units: dimensionless, range [-1, 1].
  !!
  !! @param[in] x Real(dp) array.
  !! @param[in] y Real(dp) array of same size as x.
  !! @return Real(dp), correlation coefficient. Returns NaN if size mismatch, 0 if either std = 0.
  !!
  !! Notes:
  !! - r = cov(x, y) / (std(x) * std(y)).
  !! - r = 1: perfect positive linear relationship.
  !! - r = -1: perfect negative linear relationship.
  !! - r = 0: no linear relationship.
  function corr(x, y) result(r)
    real(dp), dimension(:), intent(in) :: x, y
    real(dp) :: r
    real(dp) :: mx, my, sx, sy, sxy
    integer :: n

    n = size(x)
    if (size(y) /= n) then
      r = nan()
      return
    end if

    mx = sum(x)/n
    my = sum(y)/n
    sx = sqrt(sum((x - mx)**2)/(n - 1))
    sy = sqrt(sum((y - my)**2)/(n - 1))
    sxy = sum((x - mx)*(y - my))/(n - 1)

    if (sx > 0.0_dp .and. sy > 0.0_dp) then
      r = sxy/(sx*sy)
    else
      r = 0.0_dp
    end if
  end function corr


  !> @brief Compute covariance between arrays x and y.
  !!
  !! Units: (units of x) * (units of y).
  !!
  !! @param[in] x Real(dp) array.
  !! @param[in] y Real(dp) array of same size as x.
  !! @param[in] biased Logical, optional. If true, divide by N. Default: false (divide by N-1).
  !! @return Real(dp), covariance. Returns NaN if size mismatch.
  !!
  !! Notes:
  !! - cov(x, y) = sum((x - mean(x))*(y - mean(y)))/(N-1).
  real(dp) function cov(x, y, biased) 
    real(dp), dimension(:), intent(in) :: x, y
    logical, intent(in), optional :: biased
    real(dp) :: mx, my
    integer :: n
    logical :: isbiased

    n = size(x)
    if (size(y) /= n) then
      cov = nan()
      return
    end if

    isbiased = .false.
    if (present(biased)) isbiased = biased

    mx = sum(x)/n
    my = sum(y)/n

    if (isbiased) then
      cov = sum((x - mx)*(y - my))/n
    else
      cov = sum((x - mx)*(y - my))/(n - 1)
    end if
  end function cov


  !> @brief Compute z-scores (standardized values) of array x.
  !!
  !! Units: dimensionless.
  !!
  !! @param[in] x Real(dp) array, data values.
  !! @return Real(dp) array of same size, z = (x - mean(x))/std(x).
  !!
  !! Notes:
  !! - If std(x) = 0, returns array of zeros.
  function zscore(x) result(z)
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(size(x)) :: z
    real(dp) :: mx, sx
    integer :: n

    n = size(x)
    mx = sum(x)/n
    sx = sqrt(sum((x - mx)**2)/(n - 1))

    if (sx > 0.0_dp) then
      z = (x - mx)/sx
    else
      z = 0.0_dp
    end if
  end function zscore


  !> @brief Compute trimmed mean (mean after removing extreme values from both tails).
  !!
  !! Units: same as x.
  !!
  !! @param[in] x Real(dp) array, data values.
  !! @param[in] percent Real(dp), percentage to trim from each tail (0 to 50). Default: 5.
  !! @return Real(dp), trimmed mean.
  !!
  !! Notes:
  !! - For percent=5, removes lowest 5% and highest 5%, then computes mean of remaining 90%.
  !! - Robust to outliers.
  function trimmed_mean(x, percent) result(tm)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(in), optional :: percent
    real(dp) :: tm
    real(dp) :: pct
    integer :: n, n_trim
    real(dp), dimension(size(x)) :: x_sorted

    pct = 5.0_dp
    if (present(percent)) pct = percent

    n = size(x)
    x_sorted = x
    call quicksort(x_sorted)

    n_trim = int(n * pct / 100.0_dp)
    if (n_trim >= n/2) n_trim = 0

    if (n_trim == 0) then
      tm = sum(x_sorted)/n
    else
      tm = sum(x_sorted(n_trim+1:n-n_trim))/(n - 2*n_trim)
    end if
  end function trimmed_mean


  !> @brief Compute comprehensive statistical summary of array x.
  !!
  !! @param[in] x Real(dp) array, data values.
  !! @param[out] xstats Type(stats), statistical summary structure.
  !! @param[in] unit Integer, optional. If present, print formatted summary to this unit.
  !!
  !! Notes:
  !! - Computes all fields of the stats type.
  !! - Outlier fences based on Tukey's method (1.5*IQR and 3.0*IQR).
  !! - If unit is provided, prints human-readable summary with outlier flags.
  subroutine get_stats(x, xstats)
    real(dp), dimension(:), intent(in) :: x
    type(stats), intent(out) :: xstats
    integer :: N

    N = size(x)
    if (N <= 0) return

    xstats%n = N
    xstats%min = minval(x)
    xstats%max = maxval(x)
    xstats%median = median(x)
    xstats%trimmed_mean = trimmed_mean(x, 5.0_dp)

    call moment(x, xstats%mean, xstats%adev, xstats%std, xstats%var, &
                xstats%skew, xstats%kurt)

    xstats%rms = sqrt(dot_product(x(1:N), x(1:N))/N)

    xstats%p01 = percentile(x, 1.0_dp)
    xstats%p05 = percentile(x, 5.0_dp)
    xstats%p95 = percentile(x, 95.0_dp)
    xstats%p99 = percentile(x, 99.0_dp)

    xstats%p25 = percentile(x, 25.0_dp)
    xstats%p75 = percentile(x, 75.0_dp)
    xstats%IQR = xstats%p75 - xstats%p25

    xstats%lower_outer_fence = xstats%p25 - 3.0_dp*xstats%IQR
    xstats%lower_inner_fence = xstats%p25 - 1.5_dp*xstats%IQR
    xstats%upper_inner_fence = xstats%p75 + 1.5_dp*xstats%IQR
    xstats%upper_outer_fence = xstats%p75 + 3.0_dp*xstats%IQR

  end subroutine get_stats


  !> @brief Show comprehensive statistical summary of array x.
  !!
  !! @param[in] xstats Type(stats), statistical summary structure.
  !! @param[in] unit Integer, optional. If present, print formatted summary to this unit.
  !!
  !! Notes:
  !! - Prints human-readable summary with outlier flags.
  subroutine stats_show(xstats, unit)

    class(stats), intent(in) :: xstats
    integer, intent(in), optional :: unit
    integer :: iu
    character(len=1) :: minflag, maxflag

    if (present(unit)) then
      iu = unit
    else 
      iu = output_unit
    end if

    minflag = ''
    maxflag = ''
    if (xstats%min < xstats%lower_outer_fence) minflag = '*'
    if (xstats%max > xstats%upper_outer_fence) maxflag = '*'
    write(iu, 30) xstats%n, &
        xstats%min, minflag, xstats%p01, xstats%p05, &
        xstats%p95, xstats%p99, xstats%max, maxflag, &
        xstats%mean, xstats%trimmed_mean, xstats%std, &
        xstats%var, xstats%rms, &
        xstats%median, xstats%adev, xstats%p25, xstats%p75, &
        xstats%iqr, xstats%skew, xstats%kurt, &
        xstats%lower_inner_fence, xstats%upper_inner_fence, &
        xstats%lower_outer_fence, xstats%upper_outer_fence

30  format(10x, '    Number of points = ', i9, /, &
           10x, '    Minimum = ', F9.3, 5X, A1, /, &
           10x, '    Percentil 01% = ', F9.3, /, &
           10x, '    Percentil 05% = ', F9.3, /, &
           10x, '    Percentil 95% = ', F9.3, /, &
           10x, '    Percentil 99% = ', F9.3, /, &
           10x, '    Maximum = ', F9.3, 5X, A1, /, &
           10x, '    Mean = ', F9.3, /, &
           10x, '    Trimmed Mean = ', F9.3, /, &
           10x, '   Standard Deviation = ', F9.3, /, &
           10x, '    Variance = ', F9.3, /, &
           10x, '    RMS = ', F9.3, /, &
           10x, '    Median = ', F9.3, /, &
           10x, '    Average Deviation = ', F9.3, /, &
           10x, '    Percentil 25% = ', F9.3, /, &
           10x, '    Percentil 75% = ', F9.3, /, &
           10x, '    IQR = ', F9.3, /, &
           10x, '    Skewness = ', F9.3, /, &
           10x, '    Kurtosis = ', F9.3, /, &
           10x, '   Inner Fences ', /, &
           10x, '    Lower = ', F9.3/, &
           10x, '    Upper = ', F9.3/, &
           10x, '   Outer Fences ', /, &
           10x, '    Lower = ', F9.3/, &
           10x, '    Upper = ', F9.3/)

  end subroutine stats_show


  !> @brief Compute Spearman rank correlation coefficient between arrays x and y.
  !!
  !! Units: dimensionless, range [-1, 1].
  !!
  !! @param[in] x Real(dp) array.
  !! @param[in] y Real(dp) array of same size as x.
  !! @return Real(dp), Spearman's rho. Returns NaN if size mismatch.
  !!
  !! Notes:
  !! - Measures monotonic (not just linear) association.
  !! - Computed as Pearson correlation of ranks.
  !! - Ties are handled by averaging ranks.
  !! - rho = 1: perfect monotonic increasing relationship.
  !! - rho = -1: perfect monotonic decreasing relationship.
  !! - rho = 0: no monotonic relationship.
  function spearman(x, y) result(rho)
    real(dp), dimension(:), intent(in) :: x, y
    real(dp) :: rho
    real(dp), dimension(size(x)) :: rank_x, rank_y
    integer :: n

    n = size(x)
    if (size(y) /= n) then
      rho = nan()
      return
    end if

    ! Compute ranks for both arrays
    call compute_ranks(x, rank_x)
    call compute_ranks(y, rank_y)

    ! Pearson correlation of ranks
    rho = corr(rank_x, rank_y)

  end function spearman


  !> @brief Compute ranks of array x with tie handling (average rank for ties).
  !!
  !! @param[in] x Real(dp) array, data values.
  !! @param[out] ranks Real(dp) array of same size, ranks from 1 to n.
  !!
  !! Notes:
  !! - Uses indexx from module_math for sorting indices.
  !! - Ties receive the average of their ranks.
  !! - Example: [10, 20, 20, 30] -> [1.0, 2.5, 2.5, 4.0]
  subroutine compute_ranks(x, ranks)
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(:), intent(out) :: ranks
    integer, dimension(size(x)) :: idx
    integer :: n, i, j, tie_start
    real(dp) :: tie_sum

    n = size(x)
  
    ! Get sorted indices
    call indexx(x, idx)
  
    ! Assign ranks with tie handling
    i = 1
    do while (i <= n)
      tie_start = i
      tie_sum = real(i, dp)
    
      ! Find extent of ties
      do while (i < n)
        if (abs(x(idx(i+1)) - x(idx(i))) < epsilon(1.0_dp) * abs(x(idx(i)))) then
          i = i + 1
          tie_sum = tie_sum + real(i, dp)
        else
          exit
        end if
      end do
    
      ! Assign average rank to all tied values
      do j = tie_start, i
        ranks(idx(j)) = tie_sum / real(i - tie_start + 1, dp)
      end do
    
      i = i + 1
    end do

  end subroutine compute_ranks


  !> @brief Fit simple linear regression y = a + b*x using ordinary least squares.
  !!
  !! @param[in] x Real(dp) array, independent variable.
  !! @param[in] y Real(dp) array, dependent variable (same size as x).
  !! @param[out] a Real(dp), intercept. Units: same as y.
  !! @param[out] b Real(dp), slope. Units: (units of y)/(units of x).
  !! @param[out] r2 Real(dp), coefficient of determination (R²), dimensionless [0, 1].
  !!
  !! Notes:
  !! - Uses numerically stable formulas to avoid catastrophic cancellation.
  !! - R² = 1 - SS_res/SS_tot, where:
  !!   - SS_res = sum((y - y_pred)²) (residual sum of squares)
  !!   - SS_tot = sum((y - mean(y))²) (total sum of squares)
  !! - R² = 1: perfect fit.
  !! - R² = 0: model no better than horizontal line at mean(y).
  !! - If var(x) = 0, sets b=0, a=mean(y), r2=NaN.
  !! - Requires n >= 2.
  subroutine linreg(x, y, a, b, r2)
    real(dp), dimension(:), intent(in) :: x, y
    real(dp), intent(out) :: a, b, r2
    integer :: n, i
    real(dp) :: mx, my, sxx, sxy, syy, ss_res, ss_tot, y_pred

    n = size(x)
    if (size(y) /= n) call crash('linreg: x and y must have same size')
    if (n < 2) call crash('linreg: need at least 2 points')

    ! Compute means
    mx = sum(x) / n
    my = sum(y) / n

    ! Compute sums of squares and cross-products using stable formulas
    sxx = 0.0_dp
    sxy = 0.0_dp
    syy = 0.0_dp
  
    do i = 1, n
      sxx = sxx + (x(i) - mx) * (x(i) - mx)
      sxy = sxy + (x(i) - mx) * (y(i) - my)
      syy = syy + (y(i) - my) * (y(i) - my)
    end do

    ! Compute slope and intercept
    if (sxx > 0.0_dp) then
      b = sxy / sxx
      a = my - b * mx
    else
      ! All x values identical: no slope defined
      b = 0.0_dp
      a = my
      r2 = nan()
      return
    end if

    ! Compute R² = 1 - SS_res/SS_tot
    ss_res = 0.0_dp
    do i = 1, n
      y_pred = a + b * x(i)
      ss_res = ss_res + (y(i) - y_pred)**2
    end do
  
    ss_tot = syy
  
    if (ss_tot > 0.0_dp) then
      r2 = 1.0_dp - ss_res / ss_tot
    else
      ! All y values identical: perfect fit (no variance to explain)
      r2 = 1.0_dp
    end if

  end subroutine linreg


  !> @brief One-sample t-test: test if sample mean differs from population mean mu0.
  !!
  !! H0: mean(x) = mu0
  !! H1: mean(x) ≠ mu0 (two-tailed)
  !!
  !! @param[in] x Real(dp) array, sample data.
  !! @param[in] mu0 Real(dp), hypothesized population mean.
  !! @param[in] alpha Real(dp), optional. Significance level for CI (default: 0.05 for 95% CI).
  !! @return Type(ttest_result), test results.
  !!
  !! Notes:
  !! - Assumes x is normally distributed or n is large (CLT).
  !! - Uses Student's t-distribution with n-1 degrees of freedom.
  !! - p-value computed using approximate t-distribution CDF.
  function ttest_1samp(x, mu0, alpha) result(res)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(in) :: mu0
    real(dp), intent(in), optional :: alpha
    type(ttest_result) :: res
    integer :: n
    real(dp) :: xbar, s, se, t_crit, alph

    n = size(x)
    if (n < 2) call crash('ttest_1samp: need at least 2 samples')

    alph = 0.05_dp
    if (present(alpha)) alph = alpha

    xbar = mean(x)
    s = std(x)
    se = s / sqrt(real(n, dp))

    res%mean_diff = xbar - mu0
    res%std_err = se
    res%df = real(n - 1, dp)
    res%t_stat = res%mean_diff / se

    ! Approximate two-tailed p-value using t-distribution
    res%p_value = 2.0_dp * t_cdf_complement(abs(res%t_stat), res%df)

    ! Confidence interval: mean_diff ± t_crit * se
    t_crit = t_inverse_cdf(1.0_dp - alph/2.0_dp, res%df)
    res%ci_lower = res%mean_diff - t_crit * se
    res%ci_upper = res%mean_diff + t_crit * se

  end function ttest_1samp


  !> @brief Two-sample independent t-test: test if two sample means differ.
  !!
  !! H0: mean(x) = mean(y)
  !! H1: mean(x) ≠ mean(y) (two-tailed)
  !!
  !! @param[in] x Real(dp) array, first sample.
  !! @param[in] y Real(dp) array, second sample.
  !! @param[in] equal_var Logical, optional. If true, assume equal variances (pooled t-test). Default: true.
  !! @param[in] alpha Real(dp), optional. Significance level for CI (default: 0.05).
  !! @return Type(ttest_result), test results.
  !!
  !! Notes:
  !! - If equal_var=true: uses pooled variance and df = n1 + n2 - 2.
  !! - If equal_var=false: uses Welch's t-test with Satterthwaite approximation for df.
  function ttest_2samp(x, y, equal_var, alpha) result(res)
    real(dp), dimension(:), intent(in) :: x, y
    logical, intent(in), optional :: equal_var
    real(dp), intent(in), optional :: alpha
    type(ttest_result) :: res
    integer :: n1, n2
    real(dp) :: xbar, ybar, s1, s2, sp, se, t_crit, alph, v1, v2
    logical :: pooled
  
    n1 = size(x)
    n2 = size(y)
    if (n1 < 2 .or. n2 < 2) call crash('ttest_2samp: need at least 2 samples in each group')

    pooled = .true.
    if (present(equal_var)) pooled = equal_var

    alph = 0.05_dp
    if (present(alpha)) alph = alpha

    xbar = mean(x)
    ybar = mean(y)
    s1 = std(x)
    s2 = std(y)

    res%mean_diff = xbar - ybar

    if (pooled) then
      ! Pooled variance t-test
      sp = sqrt(((n1-1)*s1**2 + (n2-1)*s2**2) / (n1 + n2 - 2))
      se = sp * sqrt(1.0_dp/n1 + 1.0_dp/n2)
      res%df = real(n1 + n2 - 2, dp)
    else
      ! Welch's t-test (unequal variances)
      v1 = s1**2 / n1
      v2 = s2**2 / n2
      se = sqrt(v1 + v2)
      ! Welch-Satterthwaite degrees of freedom
      res%df = (v1 + v2)**2 / (v1**2/(n1-1) + v2**2/(n2-1))
    end if

    res%std_err = se
    res%t_stat = res%mean_diff / se
    res%p_value = 2.0_dp * t_cdf_complement(abs(res%t_stat), res%df)

    t_crit = t_inverse_cdf(1.0_dp - alph/2.0_dp, res%df)
    res%ci_lower = res%mean_diff - t_crit * se
    res%ci_upper = res%mean_diff + t_crit * se

  end function ttest_2samp


  !> @brief Paired t-test: test if mean difference between paired samples is zero.
  !!
  !! H0: mean(x - y) = 0
  !! H1: mean(x - y) ≠ 0 (two-tailed)
  !!
  !! @param[in] x Real(dp) array, first sample (before).
  !! @param[in] y Real(dp) array, second sample (after), same size as x.
  !! @param[in] alpha Real(dp), optional. Significance level for CI (default: 0.05).
  !! @return Type(ttest_result), test results.
  !!
  !! Notes:
  !! - Equivalent to one-sample t-test on differences d = x - y.
  !! - Use when samples are naturally paired (e.g., before/after measurements).
  function ttest_paired(x, y, alpha) result(res)
    real(dp), dimension(:), intent(in) :: x, y
    real(dp), intent(in), optional :: alpha
    type(ttest_result) :: res
    real(dp), dimension(size(x)) :: d
    real(dp) :: alph

    if (size(x) /= size(y)) call crash('ttest_paired: x and y must have same size')

    alph = 0.05_dp
    if (present(alpha)) alph = alpha

    d = x - y
    res = ttest_1samp(d, 0.0_dp, alph)

  end function ttest_paired


  !> @brief One-way ANOVA: test if means of multiple groups are equal.
  !!
  !! H0: μ₁ = μ₂ = ... = μₖ
  !! H1: at least one mean differs
  !!
  !! @param[in] groups Array of real(dp) arrays, each containing one group's data.
  !! @return Type(anova_result), ANOVA table results.
  !!
  !! Notes:
  !! - Assumes normality and equal variances across groups.
  !! - F-statistic = MS_between / MS_within.
  !! - p-value computed using approximate F-distribution CDF.
  function anova_oneway(groups) result(res)
    type(data_group), dimension(:), intent(in)  :: groups
    type(anova_result) :: res
    integer :: k, n_total, i, n_i
    real(dp) :: grand_mean, ss_total, group_mean
    real(dp), allocatable :: all_data(:)
    integer :: offset

    k = size(groups)
    if (k < 2) call crash('anova_oneway - need at least 2 groups')

    ! Compute total sample size
    n_total = 0
    do i = 1, k
      n_total = n_total + size(groups(i)%data)
    end do

    ! Concatenate all data for grand mean
    allocate(all_data(n_total))
    offset = 0
    do i = 1, k
      n_i = size(groups(i)%data)
      all_data(offset+1:offset+n_i) = groups(i)%data
      offset = offset + n_i
    end do

    grand_mean = mean(all_data)

    ! Compute sum of squares
    res%ss_between = 0.0_dp
    res%ss_within = 0.0_dp

    do i = 1, k
      n_i = size(groups(i)%data)
      group_mean = mean(groups(i)%data)
      res%ss_between = res%ss_between + n_i * (group_mean - grand_mean)**2
      res%ss_within = res%ss_within + sum((groups(i)%data - group_mean)**2)
    end do

    res%df_between = k - 1
    res%df_within = n_total - k

    res%ms_between = res%ss_between / res%df_between
    res%ms_within = res%ss_within / res%df_within

    res%f_stat = res%ms_between / res%ms_within
    res%p_value = f_cdf_complement(res%f_stat, real(res%df_between, dp), real(res%df_within, dp))

    deallocate(all_data)

  end function anova_oneway


  !> @brief Compute confidence interval for the mean.
  !!
  !! @param[in] x Real(dp) array, sample data.
  !! @param[in] alpha Real(dp), optional. Significance level (default: 0.05 for 95% CI).
  !! @param[out] ci_lower Real(dp), lower bound of CI.
  !! @param[out] ci_upper Real(dp), upper bound of CI.
  !!
  !! Notes:
  !! - Uses Student's t-distribution with n-1 degrees of freedom.
  !! - CI = mean ± t_crit * (std / sqrt(n)).
  subroutine confidence_interval_mean(x, ci_lower, ci_upper, alpha)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(out) :: ci_lower, ci_upper
    real(dp), intent(in), optional :: alpha
    integer :: n
    real(dp) :: xbar, s, se, t_crit, alph, df

    n = size(x)
    if (n < 2) call crash('confidence_interval_mean: need at least 2 samples')

    alph = 0.05_dp
    if (present(alpha)) alph = alpha

    xbar = mean(x)
    s = std(x)
    se = s / sqrt(real(n, dp))
    df = real(n - 1, dp)

    t_crit = t_inverse_cdf(1.0_dp - alph/2.0_dp, df)

    ci_lower = xbar - t_crit * se
    ci_upper = xbar + t_crit * se

  end subroutine confidence_interval_mean


  !> @brief Compute confidence intervals for linear regression coefficients.
  !!
  !! @param[in] x Real(dp) array, independent variable.
  !! @param[in] y Real(dp) array, dependent variable.
  !! @param[in] alpha Real(dp), optional. Significance level (default: 0.05).
  !! @param[out] a Real(dp), intercept estimate.
  !! @param[out] b Real(dp), slope estimate.
  !! @param[out] a_ci_lower Real(dp), lower bound of CI for intercept.
  !! @param[out] a_ci_upper Real(dp), upper bound of CI for intercept.
  !! @param[out] b_ci_lower Real(dp), lower bound of CI for slope.
  !! @param[out] b_ci_upper Real(dp), upper bound of CI for slope.
  !!
  !! Notes:
  !! - Uses standard errors from OLS theory.
  !! - SE(b) = sqrt(MSE / sum((x - mean(x))²)).
  !! - SE(a) = SE(b) * sqrt(sum(x²) / n).
  subroutine confidence_interval_linreg(x, y, a, b, a_ci_lower, a_ci_upper, &
                                       b_ci_lower, b_ci_upper, alpha)
    real(dp), dimension(:), intent(in) :: x, y
    real(dp), intent(out) :: a, b, a_ci_lower, a_ci_upper, b_ci_lower, b_ci_upper
    real(dp), intent(in), optional :: alpha
    integer :: n, i
    real(dp) :: mx, my, sxx, mse, se_b, se_a, t_crit, alph, df, r2_dummy, y_pred

    n = size(x)
    if (size(y) /= n) call crash('confidence_interval_linreg: x and y must have same size')
    if (n < 3) call crash('confidence_interval_linreg: need at least 3 points')

    alph = 0.05_dp
    if (present(alpha)) alph = alpha

    ! Fit regression
    call linreg(x, y, a, b, r2_dummy)

    ! Compute mean squared error (MSE)
    mx = mean(x)
    my = mean(y)
    sxx = sum((x - mx)**2)

    mse = 0.0_dp
    do i = 1, n
      y_pred = a + b * x(i)
      mse = mse + (y(i) - y_pred)**2
    end do
    mse = mse / (n - 2)

    ! Standard errors
    se_b = sqrt(mse / sxx)
    se_a = se_b * sqrt(sum(x**2) / n)

    df = real(n - 2, dp)
    t_crit = t_inverse_cdf(1.0_dp - alph/2.0_dp, df)

    b_ci_lower = b - t_crit * se_b
    b_ci_upper = b + t_crit * se_b
    a_ci_lower = a - t_crit * se_a
    a_ci_upper = a + t_crit * se_a

  end subroutine confidence_interval_linreg


  !> @brief Bootstrap confidence interval for the mean.
  !!
  !! @param[in] x Real(dp) array, sample data.
  !! @param[in] n_boot Integer, number of bootstrap samples (default: 10000).
  !! @param[in] alpha Real(dp), optional. Significance level (default: 0.05).
  !! @param[out] ci_lower Real(dp), lower bound of bootstrap CI.
  !! @param[out] ci_upper Real(dp), upper bound of bootstrap CI.
  !! @param[out] boot_mean Real(dp), mean of bootstrap distribution.
  !!
  !! Notes:
  !! - Uses percentile method: CI = [percentile(alpha/2), percentile(1-alpha/2)] of bootstrap means.
  !! - Requires module_math random number generator (assumed available).
  subroutine bootstrap_mean(x, ci_lower, ci_upper, boot_mean, n_boot, alpha)
    use module_math, only: random_integer
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(out) :: ci_lower, ci_upper, boot_mean
    integer, intent(in), optional :: n_boot
    real(dp), intent(in), optional :: alpha
    integer :: n, n_bootstrap, i, j, idx
    real(dp), allocatable :: boot_means(:), x_boot(:)
    real(dp) :: alph
  
    n = size(x)
    if (n < 2) call crash('bootstrap_mean: need at least 2 samples')

    n_bootstrap = 10000
    if (present(n_boot)) n_bootstrap = n_boot

    alph = 0.05_dp
    if (present(alpha)) alph = alpha

    allocate(boot_means(n_bootstrap), x_boot(n))

    ! Generate bootstrap samples
    do i = 1, n_bootstrap
      do j = 1, n
        idx = random_integer(1, n)
        x_boot(j) = x(idx)
      end do
      boot_means(i) = mean(x_boot)
    end do

    boot_mean = mean(boot_means)
    ci_lower = percentile(boot_means, 100.0_dp * alph/2.0_dp)
    ci_upper = percentile(boot_means, 100.0_dp * (1.0_dp - alph/2.0_dp))

    deallocate(boot_means, x_boot)

  end subroutine bootstrap_mean


  !> @brief Generic bootstrap framework: apply user function to bootstrap samples.
  !!
  !! @param[in] x Real(dp) array, sample data.
  !! @param[in] func Function pointer, user-defined statistic function.
  !! @param[in] n_boot Integer, number of bootstrap samples.
  !! @param[out] boot_stats Real(dp) array, bootstrap distribution of the statistic.
  !!
  !! Notes:
  !! - User provides a function with interface: real(dp) function func(x).
  !! - Returns array of bootstrap statistics for further analysis.
  subroutine bootstrap_generic(x, func, n_boot, boot_stats)
    use module_math, only: random_integer
    real(dp), dimension(:), intent(in) :: x
    integer, intent(in) :: n_boot
    real(dp), dimension(n_boot), intent(out) :: boot_stats
    interface
      function func(x) result(stat)
        import :: dp
        real(dp), dimension(:), intent(in) :: x
        real(dp) :: stat
      end function func
    end interface
    integer :: n, i, j, idx
    real(dp), allocatable :: x_boot(:)
  
    n = size(x)
    allocate(x_boot(n))

    do i = 1, n_boot
      do j = 1, n
        idx = random_integer(1, n)
        x_boot(j) = x(idx)
      end do
      boot_stats(i) = func(x_boot)
    end do

    deallocate(x_boot)

  end subroutine bootstrap_generic


  !> @brief Compute autocorrelation at lag k.
  !!
  !! @param[in] x Real(dp) array, time series data.
  !! @param[in] k Integer, lag (must be >= 0 and < size(x)).
  !! @return Real(dp), autocorrelation coefficient at lag k, range [-1, 1].
  !!
  !! Notes:
  !! - ACF(k) = cov(x(1:n-k), x(1+k:n)) / var(x).
  !! - ACF(0) = 1 by definition.
  !! - For white noise, ACF(k) ≈ 0 for k > 0.
  function autocorr(x, k) result(acf_k)
    real(dp), dimension(:), intent(in) :: x
    integer, intent(in) :: k
    real(dp) :: acf_k
    integer :: n
    real(dp) :: mx, var_x, cov_k

    n = size(x)
    if (k < 0 .or. k >= n) call crash('autocorr: lag k must be in [0, n-1]')
  
    if (k == 0) then
      acf_k = 1.0_dp
      return
    end if

    mx = mean(x)
    var_x = sum((x - mx)**2) / n

    if (var_x == 0.0_dp) then
      acf_k = 0.0_dp
      return
    end if

    cov_k = sum((x(1:n-k) - mx) * (x(1+k:n) - mx)) / n
    acf_k = cov_k / var_x

  end function autocorr


  !> @brief Compute cross-correlation between x and y at lag k.
  !!
  !! @param[in] x Real(dp) array, first time series.
  !! @param[in] y Real(dp) array, second time series (same size as x).
  !! @param[in] k Integer, lag. If k > 0, y leads x. If k < 0, x leads y.
  !! @return Real(dp), cross-correlation coefficient at lag k.
  !!
  !! Notes:
  !! - CCF(k) = cov(x(1:n-|k|), y(1+k:n)) / (std(x) * std(y)) for k >= 0.
  !! - For k < 0, shifts x instead of y.
  function crosscorr(x, y, k) result(ccf_k)
    real(dp), dimension(:), intent(in) :: x, y
    integer, intent(in) :: k
    real(dp) :: ccf_k
    integer :: n, abs_k
    real(dp) :: mx, my, sx, sy, cov_k

    n = size(x)
    if (size(y) /= n) call crash('crosscorr: x and y must have same size')

    abs_k = abs(k)
    if (abs_k >= n) call crash('crosscorr: |lag| must be < n')

    mx = mean(x)
    my = mean(y)
    sx = std(x)
    sy = std(y)

    if (sx == 0.0_dp .or. sy == 0.0_dp) then
      ccf_k = 0.0_dp
      return
    end if

    if (k >= 0) then
      cov_k = sum((x(1:n-k) - mx) * (y(1+k:n) - my)) / (n - k)
    else
      cov_k = sum((x(1+abs_k:n) - mx) * (y(1:n-abs_k) - my)) / (n - abs_k)
    end if

    ccf_k = cov_k / (sx * sy)

  end function crosscorr


  !> @brief Compute autocorrelation function (ACF) for multiple lags.
  !!
  !! @param[in] x Real(dp) array, time series data.
  !! @param[in] max_lag Integer, maximum lag to compute (default: min(n/4, 40)).
  !! @param[out] acf_vals Real(dp) array, ACF values from lag 0 to max_lag.
  !!
  !! Notes:
  !! - acf_vals(1) = ACF(0) = 1.0.
  !! - acf_vals(k+1) = ACF(k) for k = 0, 1, ..., max_lag.
  subroutine acf(x, acf_vals, max_lag)
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(:), allocatable, intent(out) :: acf_vals
    integer, intent(in), optional :: max_lag
    integer :: n, k, kmax

    n = size(x)
    kmax = min(n/4, 40)
    if (present(max_lag)) kmax = max_lag
    if (kmax >= n) kmax = n - 1

    allocate(acf_vals(0:kmax))

    do k = 0, kmax
      acf_vals(k) = autocorr(x, k)
    end do

  end subroutine acf


  !> @brief Compute partial autocorrelation function (PACF) using Durbin-Levinson recursion.
  !!
  !! @param[in] x Real(dp) array, time series data.
  !! @param[in] max_lag Integer, maximum lag to compute (default: min(n/4, 40)).
  !! @param[out] pacf_vals Real(dp) array, PACF values from lag 1 to max_lag.
  !!
  !! Notes:
  !! - PACF(k) = correlation between x(t) and x(t-k) after removing linear dependence on x(t-1), ..., x(t-k+1).
  !! - Uses Durbin-Levinson algorithm for efficiency.
  !! - pacf_vals(1) = ACF(1), pacf_vals(k) computed recursively.
  subroutine pacf(x, pacf_vals, max_lag)
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(:), allocatable, intent(out) :: pacf_vals
    integer, intent(in), optional :: max_lag
    integer :: n, k, kmax, j
    real(dp), allocatable :: acf_vals(:), phi(:,:)
    real(dp) :: num, den

    n = size(x)
    kmax = min(n/4, 40)
    if (present(max_lag)) kmax = max_lag
    if (kmax >= n) kmax = n - 1

    ! Compute ACF first
    call acf(x, acf_vals, kmax)

    allocate(pacf_vals(kmax), phi(kmax, kmax))

    ! Durbin-Levinson recursion
    pacf_vals(1) = acf_vals(1)
    phi(1,1) = acf_vals(1)

    do k = 2, kmax
      num = acf_vals(k)
      do j = 1, k-1
        num = num - phi(k-1, j) * acf_vals(k-j)
      end do

      den = 1.0_dp
      do j = 1, k-1
        den = den - phi(k-1, j) * acf_vals(j)
      end do

      pacf_vals(k) = num / den
      phi(k, k) = pacf_vals(k)

      do j = 1, k-1
        phi(k, j) = phi(k-1, j) - pacf_vals(k) * phi(k-1, k-j)
      end do
    end do

    deallocate(acf_vals, phi)

  end subroutine pacf

end module module_statistics
