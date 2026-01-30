program example_randox_normal

use aloges       ! Defines dp, mean(), var(), autcorr()
use aloges, only: head_aloges
implicit none

integer, parameter :: n = 1000000
real(dp) :: x(n)
real(dp) :: x_mean, x_var, x_skew, x_kurt
real(dp) :: se_mean, se_var, se_skew, se_kurt
real(dp) :: z_mean, z_var, z_skew, z_kurt
real(dp) :: xcorr
    
call head_aloges(os_version=.True.,fortran_version=.True.)

! 1. Generate Gaussian random numbers
! mu = 0.0000; sigma = 1.0000
x = randn(n)
    
! 2. Calculate the four moments
! mean, variance, skewness and excess kurtosis (kurt-3):
call calculate_moments_twopass(x, x_mean, x_var, x_skew, x_kurt)

! 3. Theoretical Standard Errors (SE) for N(0,1)
se_mean = 1.0_dp / sqrt(real(n, dp))
se_var  = sqrt(2.0_dp / (real(n, dp) - 1.0_dp))
se_skew = sqrt(6.0_dp / real(n, dp))
se_kurt = sqrt(24.0_dp / real(n, dp))

! 4. Calculate Z-scores
z_mean = (x_mean - 0.0_dp) / se_mean 
z_var  = (x_var - 1.0_dp) / se_var
z_skew = (x_skew - 0.0_dp) / se_skew 
z_kurt = (x_kurt - 0.0_dp) / se_kurt 

! Test 1: Mean should be ≈ mu = 0.0000
! Sample mean x_mean = mean(r)
! Mean, E(xmean) = 0.0000, Var(xmean) = sigma²/n, Std(xmean) = sigma/sqrt(n)

print *
print '(T2,A25,F8.5,A)', 'Mean:', x_mean, ' (expected: 0.0)'
print '(T2,A25,F8.5)', 'Absolute error: ', abs(x_mean-0.0_dp)
print '(T2,A25,F8.5)', 'Mean error Z-score: ', z_mean
print '(T35,A)', 'Z-score < 1: Great'
print '(T35,A)', 'Z-score < 2: Acceptable'
print '(T35,A)', 'Z-score > 3: RNG biased or broken'
    
! Test 2: Variance should be ≈ 1.000
! Sample variance, xvar = var(r)
! Mean, E(xvar) = sigma² = 1.000, Var(xvar) = 4 sigma⁴/(n-1), Std(xvar) = sqrt(2/(n-1))

print *
print '(T2,A25,F8.5,A)', 'Variance:', x_var, ' (expected: 1.0)'
print '(T2,A25,F8.5)', 'Absolute error: ', abs(x_var-1.0_dp)
print '(T2,A25,F8.5)', 'Variance error Z-score: ', z_var
print '(T35,A)', 'Z-score < 1: Great'
print '(T35,A)', 'Z-score < 2: Acceptable'
print '(T35,A)', 'Z-score > 3: RNG biased or broken'

! Test 3: Skewness should be ≈ 0.000
! Mean, E(skew) = 0.0, Var(skew) ≈ 6 /n, Std(skew) ≈ sqrt(6/n)

print *
print '(T2,A25,F8.5,A)', 'Skewness:', x_skew, ' (expected: 0.0)'
print '(T2,A25,F8.5)', 'Skewness Z-score: ', z_skew
print '(T35,A)', 'Z-score < 1: Great'
print '(T35,A)', 'Z-score < 2: Acceptable'
print '(T35,A)', 'Z-score > 3: RNG biased or broken'

! Test 4: Excess Kurtosis be ≈ 0.000
! Mean, E(kurt) = 0.0, Var(skew) ≈ 24 /n, Std(skew) ≈ sqrt(24/n)

print *
print '(T2,A25,F8.5,A)', 'Kurtosis (Excess):', x_kurt, ' (expected: 0.0)'
print '(T2,A25,F8.5)', 'Kurtosis Z-score: ', z_kurt
print '(T35,A)', 'Z-score < 1: Great'
print '(T35,A)', 'Z-score < 2: Acceptable'
print '(T35,A)', 'Z-score > 3: RNG biased or broken'


! Test 5: Lagged correlation:
xcorr = autocorr(x,1)
print *
print '(T2,A25,F8.5)', 'Lag-1 correlation : ', xcorr
print '(T2,A25,F8.5)', 'Acceptance criterium : ', 3.0/sqrt(real(n,dp))
print *, '     |ρ₁| < accepance criterium : ', xcorr < 3.0/sqrt(real(n,dp))
    
contains

  subroutine calculate_moments_twopass(x, xmean, xvar, xskew, xkurt)

    ! Two-pass shifted algorithm: Usually he fastest method when the data
    ! fits in memory. By shifting the sata by an estimate of the mean (usually
    ! the value of the first element or a quick mean of a small subsample),
    ! floating point errors can be reduced.
    ! First pass: calculate the mean using sasum
    use, intrinsic :: iso_fortran_env
    implicit none
    
    ! Arguments
    real(real64), intent(in) :: x(:)
    real(real64), intent(out) :: xmean, xvar, xskew, xkurt
    
    ! Local variables
    integer :: n, i
    real(real64) :: shift, sum_diff, sum2, sum3, sum4
    real(real64) :: diff, diff2, s_std
    real(real64), allocatable :: x_shifted(:)
    
    ! External BLAS functions
    real(real64), external :: dasum, dnrm2


    xmean = 0.0; xvar = 0.0; xskew = 0.0; xkurt = 0.0

    n = size(x)
    if (n.eq.0) return
    
    allocate(x_shifted(n))

    ! Pass 1: Use the first element as a shift to maintain numerical stability
    shift = x(1)
    
    ! Calculate shifted sum for the mean: sum(x_i - shift)
    ! We use a simple loop here, or could use dasum if we handle signs
    sum_diff = 0.0_real64
    do i = 1, n
        x_shifted(i) = x(i) - shift
        sum_diff = sum_diff + x_shifted(i)
    end do
    
    xmean = shift + (sum_diff / real(n, real64))

    ! Now re-center x_shifted so it is exactly relative to the true sample mean
    ! x_shifted = x - xmean
    do i = 1, n
        x_shifted(i) = x(i) - xmean
    end do

    ! Pass 2: Calculate higher moments
    ! Variance using BLAS dnrm2: sqrt(sum(diff^2))
    ! sum2 = (dnrm2(n, x_shifted, 1))**2
    sum2 = 0.0_real64
    sum3 = 0.0_real64
    sum4 = 0.0_real64
    
    do i = 1, n
        diff  = x_shifted(i)
        diff2 = diff * diff
        sum2  = sum2 + diff2
        sum3  = sum3 + diff2 * diff
        sum4  = sum4 + diff2 * diff2
    end do

    ! Finalize Statistics
    xvar  = sum2 / real(n - 1, real64)
    s_std = sqrt(xvar)
    
    xskew = (sum3 / real(n, real64)) / (s_std**3)
    xkurt = (sum4 / real(n, real64)) / (xvar**2) - 3.0_real64

    deallocate(x_shifted)
  end subroutine calculate_moments_twopass

end program example_randox_normal
