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
! - power_spectrum_1d                                                      !
! - power_spectrum_ensemble                                                !
! -------------------------------------------------------------------------!

module module_spectra

use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
use module_types, only : dp
use module_constants, only : pi, two_pi, nan
use module_tools, only: crash
use module_math, only: next_power_of_2
use module_fft, only: fft1d

implicit none
private
    
public :: power_spectrum_1d, power_spectrum_ensemble
public :: check_uniform_spacing, interpolate_to_uniform

contains

  ! ...
  ! ====================================================================
  ! ====================================================================
  ! ...
  !> Calculate ensemble-averaged power spectrum with confidence intervals
  subroutine power_spectrum_ensemble(x, data, wavenumber, wavelength, &
                                     power_mean, power_std, power_ci_lower, power_ci_upper, &
                                     n_spec, detrend, confidence_level)
    
    real(dp), dimension(:), intent(in)     :: x              ! distance (km)
    real(dp), dimension(:,:), intent(in)   :: data           ! data(1:n, 1:nseries)
    integer, intent(out)                   :: n_spec
    real(dp), allocatable, intent(out)     :: wavenumber(:), wavelength(:)
    real(dp), allocatable, intent(out)     :: power_mean(:)  ! Mean power
    real(dp), allocatable, intent(out)     :: power_std(:)   ! Standard deviation
    real(dp), allocatable, intent(out)     :: power_ci_lower(:), power_ci_upper(:)  ! Confidence intervals
    logical, intent(in), optional          :: detrend
    real(dp), intent(in), optional         :: confidence_level  ! Default 0.95 (95%)
    
    ! ... Local variables
    ! ...
    integer :: n, nseries, i, j
    real(dp) :: conf_level, t_value, dof
    real(dp), allocatable :: power_temp(:)
    real(dp), allocatable :: wavenumber_temp(:), wavelength_temp(:)
    real(dp), allocatable :: all_spectra(:,:)  ! (n_spec, nseries)
    integer :: n_spec_temp
    logical :: do_detrend
    
    n = size(x)
    nseries = size(data, 2)
    
    if (n /= size(data, 1)) call crash('In POWER_SPECTRUM_ENSEMBLE - Incompatible input data size')
    if (nseries < 2) call crash('In POWER_SPECTRUM_ENSEMBLE - Need at least 2 series for statistics')
    
    do_detrend = .false.
    if (present(detrend)) do_detrend = detrend
    
    conf_level = 0.95_dp
    if (present(confidence_level)) conf_level = confidence_level
    
    write(*,'(A,I0,A)') ' Computing ensemble-averaged spectrum from ', nseries, ' realizations...'
    
    ! ... Calculate spectrum for first series to get dimensions
    ! ...
    call power_spectrum_1d(x, data(:,1), wavenumber_temp, wavelength_temp, &
                          power_temp, n_spec_temp, detrend=do_detrend)
    
    n_spec = n_spec_temp
    
    ! ... Allocate arrays
    ! ...
    allocate(wavenumber(n_spec), wavelength(n_spec))
    allocate(power_mean(n_spec), power_std(n_spec))
    allocate(power_ci_lower(n_spec), power_ci_upper(n_spec))
    allocate(all_spectra(n_spec, nseries))
    
    wavenumber = wavenumber_temp
    wavelength = wavelength_temp
    all_spectra(:,1) = power_temp
    
    deallocate(wavenumber_temp, wavelength_temp, power_temp)
    
    ! ... Calculate spectra for remaining series
    ! ...
    do i = 2, nseries
      call power_spectrum_1d(x, data(:,i), wavenumber_temp, wavelength_temp, &
                            power_temp, n_spec_temp, detrend=do_detrend)
      if (n_spec_temp /= n_spec) call crash('In POWER_SPECTRUM_ENSEMBLE - Inconsistent spectrum size')
      all_spectra(:,i) = power_temp
      deallocate(wavenumber_temp, wavelength_temp, power_temp)
    end do
    
    ! ... Calculate statistics at each wavenumber
    ! ...
    do j = 1, n_spec
      ! Mean
      power_mean(j) = sum(all_spectra(j,:)) / real(nseries, dp)
      
      ! Standard deviation
      power_std(j) = sqrt(sum((all_spectra(j,:) - power_mean(j))**2) / real(nseries - 1, dp))
      
      ! Confidence intervals using t-distribution
      dof = real(nseries - 1, dp)
      t_value = t_inverse_cdf(conf_level, dof)
      
      power_ci_lower(j) = power_mean(j) - t_value * power_std(j) / sqrt(real(nseries, dp))
      power_ci_upper(j) = power_mean(j) + t_value * power_std(j) / sqrt(real(nseries, dp))
      
      ! Ensure non-negative power
      if (power_ci_lower(j) < 0.0_dp) power_ci_lower(j) = 0.0_dp
    end do
    
    !write(*,'(A)') ' Ensemble-averaged spectrum computed successfully.'
    
    deallocate(all_spectra)
    
  end subroutine power_spectrum_ensemble
  ! ...
  ! ====================================================================
  ! ...
  function t_inverse_cdf(confidence, dof) result(t_val)
    ! ... Approximate t-distribution inverse CDF for confidence intervals

    real(dp), intent(in) :: confidence  ! e.g., 0.95 for 95% CI
    real(dp), intent(in) :: dof         ! degrees of freedom

    ! ... Local variables
    ! ...
    real(dp) :: t_val
    real(dp) :: alpha, z
    
    ! For large dof (>30), use normal approximation
    if (dof > 30.0_dp) then
      ! Standard normal quantile approximation
      alpha = 1.0_dp - confidence
      z = sqrt(2.0_dp) * erfinv(1.0_dp - alpha)
      t_val = z
    else
      ! Simple lookup table for common cases (two-tailed)
      alpha = (1.0_dp - confidence) / 2.0_dp
      
      ! Approximate using polynomial (good enough for most cases)
      ! This is a simplified version - for production, use a proper t-table
      if (confidence >= 0.95_dp) then
        if (dof <= 5.0_dp) then
          t_val = 2.571_dp  ! Approximate for dof=5, 95% CI
        else if (dof <= 10.0_dp) then
          t_val = 2.228_dp  ! Approximate for dof=10, 95% CI
        else if (dof <= 20.0_dp) then
          t_val = 2.086_dp  ! Approximate for dof=20, 95% CI
        else
          t_val = 1.96_dp   ! Normal approximation
        end if
      else if (confidence >= 0.90_dp) then
        if (dof <= 5.0_dp) then
          t_val = 2.015_dp
        else if (dof <= 10.0_dp) then
          t_val = 1.812_dp
        else
          t_val = 1.645_dp
        end if
      else
        t_val = 1.96_dp  ! Default to 95% normal
      end if
    end if
    
  end function t_inverse_cdf
  ! ...
  ! ====================================================================
  ! ...
  function erfinv(x) result(y)
    ! ... Inverse error function (for normal quantiles)

    real(dp), intent(in) :: x
    real(dp) :: y, w, p
    
    ! Approximation for inverse error function
    ! Valid for |x| < 1
    w = -log((1.0_dp - x) * (1.0_dp + x))
    
    if (w < 5.0_dp) then
      w = w - 2.5_dp
      p = 2.81022636e-08_dp
      p = 3.43273939e-07_dp + p * w
      p = -3.5233877e-06_dp + p * w
      p = -4.39150654e-06_dp + p * w
      p = 0.00021858087_dp + p * w
      p = -0.00125372503_dp + p * w
      p = -0.00417768164_dp + p * w
      p = 0.246640727_dp + p * w
      p = 1.50140941_dp + p * w
    else
      w = sqrt(w) - 3.0_dp
      p = -0.000200214257_dp
      p = 0.000100950558_dp + p * w
      p = 0.00134934322_dp + p * w
      p = -0.00367342844_dp + p * w
      p = 0.00573950773_dp + p * w
      p = -0.0076224613_dp + p * w
      p = 0.00943887047_dp + p * w
      p = 1.00167406_dp + p * w
      p = 2.83297682_dp + p * w
    end if
    
    y = p * x
    
  end function erfinv
  ! ...
  ! ====================================================================
  ! ...
  subroutine power_spectrum_1d(x, data, wavenumber, wavelength, power, n_spec, detrend)
  ! ... Calculate 1D power spectrum using FFT (original routine - unchanged)

    real(dp), dimension(:), intent(in)   :: x    ! distance (km)
    real(dp), dimension(:), intent(in)   :: data    ! data
    integer, intent(out)                 :: n_spec
    real(dp), allocatable, intent(out)   :: wavenumber(:), wavelength(:), power(:)
    logical, intent(in), optional        :: detrend
    
    ! ... Local variables
    ! ...
    integer  :: i, n, n_fft, n_pad
    logical  :: is_uniform, do_detrend
    real(dp) :: dx
    real(dp), allocatable :: x_uniform(:), data_uniform(:), data_work(:)
    real(dp), allocatable :: fft_real(:), fft_imag(:)

    n = size(x)
    if (n.ne.size(data)) call crash('In POWER_SPECTRA_1D - Incompatible input data size')
    
    do_detrend = .false.
    if (present(detrend)) do_detrend = detrend
    
    ! ... Check if data is uniformly spaced
    ! ...
    call check_uniform_spacing(x, is_uniform, dx)
    if (is_uniform) then
      allocate(data_work(n))
      data_work = data
      n_fft = n
    else
      write(*,*) "Data not uniformly spaced. Interpolating..."
      allocate(x_uniform(n), data_uniform(n))
      call interpolate_to_uniform(x, data, n, x_uniform, data_uniform, n)
      dx = x_uniform(2) - x_uniform(1)
      allocate(data_work(n))
      data_work = data_uniform
      n_fft = n
      deallocate(x_uniform, data_uniform)
    end if

    ! ... Pad to next power of 2 for FFT efficiency
    ! ...
    n_pad = next_power_of_2(n_fft)
    
    ! ... Detrend if requested
    ! ...
    if (do_detrend) call detrend_linear(data_work, n_fft)
    
    ! ... Apply window function (Hanning)
    ! ...
    call apply_hanning_window(data_work, n_fft)
    
    ! ... Prepare FFT arrays (zero-padded)
    ! ...
    allocate(fft_real(0:n_pad-1), fft_imag(0:n_pad-1))
    fft_real = 0.0_dp
    fft_imag = 0.0_dp
    fft_real(0:n_fft-1) = data_work(1:n_fft)
    
    ! ... Perform FFT 
    ! ...
    call fft1d (n_pad, fft_real, fft_imag, -1)
    
    ! ... Calculate power spectrum (one-sided)
    ! ...
    n_spec = n_pad / 2 + 1
    allocate(wavenumber(n_spec), wavelength(n_spec), power(n_spec))
    
    do i = 1, n_spec
      ! ... Wavenumber (cycles per unit distance)
      wavenumber(i) = real(i - 1, dp) / (real(n_pad, dp) * dx)
    
      ! ... Wavelength: 1/k
      if (i == 1) then
        wavelength(i) = nan() 
      else
        wavelength(i) = 1.0_dp / wavenumber(i)
      endif
    
      ! Power spectral density
      if (i == 1 .or. i == n_spec) then
        power(i) = (fft_real(i-1)**2 + fft_imag(i-1)**2) / real(n_pad, dp)**2
      else
        power(i) = 2.0_dp * (fft_real(i-1)**2 + fft_imag(i-1)**2) / real(n_pad, dp)**2
      endif
    enddo
    
    deallocate(data_work, fft_real, fft_imag)
    
  end subroutine power_spectrum_1d
  ! ...
  ! ====================================================================
  ! ...
  subroutine check_uniform_spacing(x, is_uniform, dx, tolerance)
    ! ... Given an array X, check if it is, up to a tolerance limit, uniformly
    ! ... spaced. Returns the average spacing (dx).
    ! ... On input:
    ! ...   real(dp)  :: x(:)
    ! ...   real(dp)  :: tolerance (Optional, default = 1.0e-6)
    ! ... On output:
    ! ...   logical   :: is_uniform (.True. or .False.)
    ! ...   dx        :: average spacing
    ! ...
    real(dp), intent(in)                            :: x(:)
    logical, intent(out)                            :: is_uniform
    real(dp), intent(out)                           :: dx
    real(dp), intent(in), optional                  :: tolerance

    ! ... Local variables
    ! ...
    integer :: n,i
    real(dp) :: tol, dx_mean, dx_std, dx_i

    n = size(x)
    tol = 1.0e-6_dp
    if (present(tolerance)) tol = tolerance

    dx_mean = (x(n) - x(1)) / real(n - 1, dp)
    dx = dx_mean
    dx_std = 0.0_dp

    do i = 1, n - 1
      dx_i = x(i+1) - x(i)
      dx_std = dx_std + (dx_i - dx_mean)**2
    end do

    ! ... Perform the check, using the deviation of dx over the mean value
    ! ...
    dx_std = sqrt(dx_std / real(n - 1, dp))
    is_uniform = (dx_std / dx_mean) < tol

  end subroutine check_uniform_spacing
  ! ...
  ! ====================================================================
  ! ...
  subroutine interpolate_to_uniform(x_in, y_in, n_in, x_out, y_out, n_out)

    integer, intent(in) :: n_in, n_out
    real(dp), intent(in) :: x_in(n_in), y_in(n_in)
    real(dp), intent(out) :: x_out(n_out), y_out(n_out)

    ! ... Local variables
    integer :: i, j
    real(dp) :: t

    do i = 1, n_out
      x_out(i) = x_in(1) + (x_in(n_in) - x_in(1)) * real(i - 1, dp) / real(n_out - 1, dp)
    end do

    j = 1
    do i = 1, n_out
      do while (j < n_in .and. x_in(j+1) < x_out(i))
        j = j + 1
      end do
      if (j >= n_in) then
        y_out(i) = y_in(n_in)
      else if (x_out(i) <= x_in(1)) then
        y_out(i) = y_in(1)
      else
        t = (x_out(i) - x_in(j)) / (x_in(j+1) - x_in(j))
        y_out(i) = y_in(j) + t * (y_in(j+1) - y_in(j))
      end if
    end do

  end subroutine interpolate_to_uniform
  ! ...
  ! =============================================================================
  ! ...
  subroutine detrend_linear(data, n)

    integer, intent(in) :: n
    real(dp), intent(inout) :: data(n)

    ! ... Local variables
    ! ...
    integer  :: i
    real(dp) :: sum_x, sum_y, sum_xy, sum_xx
    real(dp) :: slope, intercept, xi, di

    sum_x  = 0.0_dp
    sum_y  = 0.0_dp
    sum_xx = 0.0_dp
    sum_xy = 0.0_dp
    do i = 1, n
    xi = real(i, dp)
    di = data(i)
    sum_x  = sum_x  + xi
    sum_y  = sum_y  + di
    sum_xx = sum_xx + xi * xi
    sum_xy = sum_xy + xi * di
    end do
    slope = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x**2)
    intercept = (sum_y - slope * sum_x) / n
    do i = 1, n
    xi = real(i, dp)
    data(i) = data(i) - (slope * xi + intercept)
    enddo
  end subroutine detrend_linear
  ! ...
  ! =============================================================================
  ! ...
  subroutine apply_hanning_window(data,n)
    integer, intent(in)    :: n
    real(dp), intent(inout) :: data(n)
    integer :: i
    real(dp) :: window
    do i = 1, n
    window = 0.5_dp * (1.0_dp - cos(two_pi*real(i-1,dp) / real(n-1,dp)))
    data(i) = data(i) * window
    end do
  end subroutine apply_hanning_window
  ! ...
  ! =============================================================================
  ! ...

end module module_spectra
