program example_random_uniform

use aloges       ! Defines dp, mean(), var(), autcorr()
use aloges, only: head_aloges
implicit none

integer, parameter :: n = 1000000
real(dp) :: r(n)
real(dp) :: rmean, rvar, rcorr
integer :: i, bins(10)
    
call head_aloges(fortran_version=.True.)

! Generate random numbers
call random_number(r)
    
! Test 1: Mean should be ~0.5

rmean = mean(r)
print *
print '(T2,A20,F8.5,A)', "Mean:", rmean, " (expected: 0.5)"
print '(T2,A20,F8.5)', 'Absolute error: ', rmean-0.5
print '(T2,A20,F8.5)', 'Relative error (%): ', (rmean-0.5)/0.5 * 100.0
    
! Test 2: Variance should be ~1/12 ≈ 0.0833
rvar = var(r)
print *
print '(T2,A20,F8.5,A)', "Variance:", rvar, " (expected: 0.0833)"
print '(T2,A20,F8.5)', 'Absolute error: ', rvar-0.0833
print '(T2,A20,F8.5)', 'Relative error (%): ', (rvar-0.0833)/0.0833 * 100.0
    
! Test 3: Chi-square uniformity
bins = 0
do i = 1, n
    bins(int(r(i) * 10) + 1) = bins(int(r(i) * 10) + 1) + 1
end do
    
print *
print *, "Bin counts (should be ~100000 each):"
print *, bins
    
! Chi-square statistic
print *, "Chi-square:", sum((bins - n/10.0)**2 / (n/10.0))
print *, "(critical value at 95%: 16.919 for 9 df)"

! Test 4: Correlations
rcorr = autocorr(r,1)
print *
print '(T2,A25,F8.5)', 'Lag-1 correlation : ', rcorr
print '(T2,A25,F8.5)', 'Acceptance criterium : ', 3.0/sqrt(real(n,dp))
print *, '     |ρ₁| < accepance criterium : ', rcorr < 3.0/sqrt(real(n,dp))
    
end program example_random_uniform
