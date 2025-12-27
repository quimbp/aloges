!> @file module_wavelet.f90
!> @brief Correct 2D Discrete Wavelet Transform (Daubechies-4) with periodic boundaries
module module_wavelet
  use module_constants, only : dp
  implicit none (type, external)
  private
  public :: dwt_1d, idwt_1d, dwt_2d, idwt_2d, wavelet_denoise_2d

  ! Daubechies-4 filter coefficients (orthonormal)
  ! Low-pass decomposition filter
  real(dp), parameter :: DB4_L0 =  0.4829629131445341_dp
  real(dp), parameter :: DB4_L1 =  0.8365163037378077_dp
  real(dp), parameter :: DB4_L2 =  0.2241438680420134_dp
  real(dp), parameter :: DB4_L3 = -0.1294095225512603_dp
  
  ! High-pass decomposition filter from low-pass using quadrature mirror:
  ! g_k = (-1)^k h_{N-1-k}
  real(dp), parameter :: DB4_H0 = -DB4_L3  ! = (-1)^0 * h_{3} =  0.1294095225512603
  real(dp), parameter :: DB4_H1 =  DB4_L2  ! = (-1)^1 * h_{2} = -0.2241438680420134
  real(dp), parameter :: DB4_H2 = -DB4_L1  ! = (-1)^2 * h_{1} =  0.8365163037378077
  real(dp), parameter :: DB4_H3 =  DB4_L0  ! = (-1)^3 * h_{0} = -0.4829629131445341

contains

  !> Helper function for periodic indexing
  pure function periodic_index(i, n) result(j)
    integer, intent(in) :: i, n
    integer :: j
    ! Convert to 0-based, mod, then back to 1-based
    j = mod(i-1 + n, n) + 1
  end function periodic_index

  !> 1D forward DWT with periodic extension
  subroutine dwt_1d(input, approx, detail)
    real(dp), intent(in)  :: input(:)
    real(dp), intent(out) :: approx(:), detail(:)
    integer :: i, n, n2
    
    n = size(input)
    n2 = size(approx)
    
    ! Input validation
    if (2*n2 /= n) then
      error stop "dwt_1d: size mismatch - output arrays must be half the size of input"
    end if
    
    do i = 1, n2
      ! For Daubechies-4: convolution with stride 2, 4-tap filter
      ! Indexing: x[2i + k - 1] for k = 0..3 (but Fortran is 1-indexed)
      approx(i) = DB4_L0 * input(periodic_index(2*i-1, n)) + &
                  DB4_L1 * input(periodic_index(2*i,   n)) + &
                  DB4_L2 * input(periodic_index(2*i+1, n)) + &
                  DB4_L3 * input(periodic_index(2*i+2, n))
      
      detail(i) = DB4_H0 * input(periodic_index(2*i-1, n)) + &
                  DB4_H1 * input(periodic_index(2*i,   n)) + &
                  DB4_H2 * input(periodic_index(2*i+1, n)) + &
                  DB4_H3 * input(periodic_index(2*i+2, n))
    end do
  end subroutine dwt_1d

  !> 1D inverse DWT with periodic extension
  subroutine idwt_1d(approx, detail, output)
    real(dp), intent(in)  :: approx(:), detail(:)
    real(dp), intent(out) :: output(:)
    integer :: i, n, n2
    
    n2 = size(approx)
    n = size(output)
    
    ! Input validation
    if (2*n2 /= n) then
      error stop "idwt_1d: size mismatch - input arrays must be half the size of output"
    end if
    
    output = 0.0_dp
    
    ! Inverse transform: upsample by 2 and apply filters
    do i = 1, n2
      ! Add contribution to positions
      output(periodic_index(2*i-1, n)) = output(periodic_index(2*i-1, n)) + &
                                         DB4_L0 * approx(i) + DB4_H0 * detail(i)
      output(periodic_index(2*i,   n)) = output(periodic_index(2*i,   n)) + &
                                         DB4_L1 * approx(i) + DB4_H1 * detail(i)
      output(periodic_index(2*i+1, n)) = output(periodic_index(2*i+1, n)) + &
                                         DB4_L2 * approx(i) + DB4_H2 * detail(i)
      output(periodic_index(2*i+2, n)) = output(periodic_index(2*i+2, n)) + &
                                         DB4_L3 * approx(i) + DB4_H3 * detail(i)
    end do
  end subroutine idwt_1d

  !> 2D forward DWT (separable)
  subroutine dwt_2d(input, nx, ny, LL, LH, HL, HH)
    real(dp), intent(in)  :: input(nx, ny)
    integer, intent(in)   :: nx, ny
    real(dp), intent(out) :: LL(nx/2, ny/2), LH(nx/2, ny/2)
    real(dp), intent(out) :: HL(nx/2, ny/2), HH(nx/2, ny/2)
    real(dp), allocatable :: temp_approx_rows(:,:), temp_detail_rows(:,:)
    real(dp), allocatable :: row_in(:), row_approx(:), row_detail(:)
    integer :: i, j, nx2, ny2
    
    ! Input validation
    if (mod(nx, 2) /= 0 .or. mod(ny, 2) /= 0) then
      error stop "dwt_2d: input dimensions must be even"
    end if
    
    nx2 = nx/2
    ny2 = ny/2
    
    allocate(temp_approx_rows(nx2, ny), temp_detail_rows(nx2, ny))
    allocate(row_in(nx), row_approx(nx2), row_detail(nx2))
    
    ! Transform rows
    do j = 1, ny
      row_in = input(:, j)
      call dwt_1d(row_in, row_approx, row_detail)
      temp_approx_rows(:, j) = row_approx
      temp_detail_rows(:, j) = row_detail
    end do
    
    ! Transform columns of approximation rows
    do i = 1, nx2
      call dwt_1d(temp_approx_rows(i, :), LL(i, :), LH(i, :))
    end do
    
    ! Transform columns of detail rows  
    do i = 1, nx2
      call dwt_1d(temp_detail_rows(i, :), HL(i, :), HH(i, :))
    end do
    
    deallocate(temp_approx_rows, temp_detail_rows, row_in, row_approx, row_detail)
  end subroutine dwt_2d

  !> 2D inverse DWT (separable)
  subroutine idwt_2d(LL, LH, HL, HH, nx2, ny2, output)
    real(dp), intent(in)  :: LL(nx2, ny2), LH(nx2, ny2)
    real(dp), intent(in)  :: HL(nx2, ny2), HH(nx2, ny2)
    integer, intent(in)   :: nx2, ny2
    real(dp), intent(out) :: output(2*nx2, 2*ny2)
    real(dp), allocatable :: temp_low_cols(:,:), temp_high_cols(:,:)
    real(dp), allocatable :: col_approx(:), col_detail(:), col_out(:)
    real(dp), allocatable :: row_approx(:), row_detail(:), row_out(:)
    integer :: i, j, nx, ny
    
    nx = 2*nx2
    ny = 2*ny2
    
    allocate(temp_low_cols(nx2, ny), temp_high_cols(nx2, ny))
    allocate(col_approx(ny2), col_detail(ny2), col_out(ny))
    allocate(row_approx(nx2), row_detail(nx2), row_out(nx))
    
    ! Inverse column transforms
    do i = 1, nx2
      col_approx = LL(i, :)
      col_detail = LH(i, :)
      call idwt_1d(col_approx, col_detail, col_out)
      temp_low_cols(i, :) = col_out
      
      col_approx = HL(i, :)
      col_detail = HH(i, :)
      call idwt_1d(col_approx, col_detail, col_out)
      temp_high_cols(i, :) = col_out
    end do
    
    ! Inverse row transforms
    do j = 1, ny
      row_approx = temp_low_cols(:, j)
      row_detail = temp_high_cols(:, j)
      call idwt_1d(row_approx, row_detail, row_out)
      output(:, j) = row_out
    end do
    
    deallocate(temp_low_cols, temp_high_cols, col_approx, col_detail, col_out, &
               row_approx, row_detail, row_out)
  end subroutine idwt_2d

  !> Simple wavelet denoising
  subroutine wavelet_denoise_2d(noisy, nx, ny, threshold, denoised)
    real(dp), intent(in)  :: noisy(nx, ny)
    integer, intent(in)   :: nx, ny
    real(dp), intent(in)  :: threshold
    real(dp), intent(out) :: denoised(nx, ny)
    real(dp), allocatable :: LL(:,:), LH(:,:), HL(:,:), HH(:,:)
    integer :: nx2, ny2
    
    ! Input validation
    if (mod(nx, 2) /= 0 .or. mod(ny, 2) /= 0) then
      error stop "wavelet_denoise_2d: input dimensions must be even"
    end if
    if (threshold < 0.0_dp) then
      error stop "wavelet_denoise_2d: threshold must be non-negative"
    end if
    
    nx2 = nx/2
    ny2 = ny/2
    
    allocate(LL(nx2, ny2), LH(nx2, ny2), HL(nx2, ny2), HH(nx2, ny2))
    
    call dwt_2d(noisy, nx, ny, LL, LH, HL, HH)
    
    ! Apply soft thresholding to detail coefficients
    LH = soft_threshold(LH, threshold)
    HL = soft_threshold(HL, threshold)
    HH = soft_threshold(HH, threshold)
    
    call idwt_2d(LL, LH, HL, HH, nx2, ny2, denoised)
    
    deallocate(LL, LH, HL, HH)
  end subroutine wavelet_denoise_2d

  !> Soft thresholding function
  pure function soft_threshold(x, lambda) result(y)
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(in) :: lambda
    real(dp) :: y(size(x,1), size(x,2))
    real(dp) :: abs_val
    
    integer :: i, j
    
    ! Manual implementation to avoid temporary array creation
    do j = 1, size(x, 2)
      do i = 1, size(x, 1)
        abs_val = abs(x(i,j))
        if (abs_val > lambda) then
          y(i,j) = sign(abs_val - lambda, x(i,j))
        else
          y(i,j) = 0.0_dp
        end if
      end do
    end do
  end function soft_threshold

end module module_wavelet
