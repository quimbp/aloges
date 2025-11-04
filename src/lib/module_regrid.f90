!> @file module_regrid.f90
!> @brief Regrid: regrid2d with multiple methods
!> @author Quim Ballabrera, Institut de Ciencies del Mar, CSIC
!> @date October 2025
!>
!> This module provides regridding functions for coarsening high-resolution
!> fields to lower resolution grids using various conservative methods.
!>
!> Copyright (C) 2022-2025, Joaquim Ballabrera
!> Licensed under GNU LGPL v3+

module module_regrid

use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
use iso_fortran_env, only: error_unit, output_unit
use module_types, only: dp
use module_tools, only: crash, locate, lowercase, ltrim, cell_bounds2d, quicksort

implicit none (type, external)

private

! Minimal CSR type (0-based col indices expected; we’ll handle 1-based conversion)
type :: csr_matrix
  integer :: n_rows = 0
  integer :: n_cols = 0
  integer, allocatable :: row_ptr(:)   ! length n_rows+1, 1-based in Fortran, stores 0-based cumulative nnz
  integer, allocatable :: col_ind(:)   ! length nnz, 1-based in Fortran, stores 1-based col indices
  real(dp), allocatable :: values(:)   ! length nnz
end type csr_matrix

public :: regrid2d
public :: regrid2d_mean, regrid2d_min, regrid2d_max, regrid2d_median
public :: regrid2d_area_weighted
public :: regrid2d_sl_bilinear
!public :: regrid2d_learned  

contains

  pure subroutine estimate_index_window(x_line, valL, valR, i_start, i_end)
    real(dp), intent(in) :: x_line(:)
    real(dp), intent(in) :: valL, valR
    integer, intent(out) :: i_start, i_end
    integer :: n
    n = size(x_line)
    i_start = max(1, locate(x_line, valL))
    i_end   = min(n, locate(x_line, valR) + 1)
  end subroutine estimate_index_window


  !> @brief Generic regridding function with method selection  
  !!  
  !! Dispatches to specific regridding method based on string argument  
  !! Available methods: 'mean', 'min', 'max', 'median', 'area_weighted'  
  !!  
  !! @param xi X-coordinates of fine grid  
  !! @param yi Y-coordinates of fine grid  
  !! @param Fi Field values on fine grid  
  !! @param xc X-coordinates of coarse grid  
  !! @param yc Y-coordinates of coarse grid  
  !! @param method Regridding method ('mean', 'min', 'max', 'median', 'area_weighted')  
  !! @return Fc Regridded field on coarse grid  
  !!  
  !! Example:  
  !!   Fc = regrid2d(xi, yi, Fi, xc, yc, 'min')  
  !!   Fc = regrid2d(xi, yi, Fi, xc, yc, 'area_weighted')  
  !!  
  function regrid2d(xi, yi, Fi, xc, yc, method) result(Fc)
    real(dp), intent(in) :: xi(:,:), yi(:,:), Fi(:,:)
    real(dp), intent(in) :: xc(:,:), yc(:,:)
    character(len=*), intent(in) :: method
    real(dp), allocatable :: Fc(:,:)
    character(len=:), allocatable :: method_lower

    if (any(shape(xi) /= shape(yi))) call crash('regrid2d - xi/yi shape mismatch')
    if (any(shape(xc) /= shape(yc))) call crash('regrid2d - xc/yc shape mismatch')
    if (any(shape(Fi) /= shape(xi))) call crash('regrid2d - Fi shape mismatch with xi/yi')

    method_lower = lowercase(ltrim(method))
    select case (method_lower)
    case ('mean','average','avg');                Fc = regrid2d_mean(xi,yi,Fi,xc,yc)
    case ('min','minimum');                       Fc = regrid2d_min(xi,yi,Fi,xc,yc)
    case ('max','maximum');                       Fc = regrid2d_max(xi,yi,Fi,xc,yc)
    case ('median','med');                        Fc = regrid2d_median(xi,yi,Fi,xc,yc)
    case ('area_weighted','area','conservative'); Fc = regrid2d_area_weighted(xi,yi,Fi,xc,yc)
    case ('slope_limited','sl_bilinear','2nd_order','second_order'); Fc = regrid2d_sl_bilinear(xi,yi,Fi,xc,yc)
!    case ('learned','nn_weights','ml'); call crash('regrid2d(method="learned") requires CSR weights; use regrid2d_learned')
    case default
      call crash('regrid2d - Unknown method: '//trim(method))
    end select
  end function regrid2d
  
  
  !> @brief Mean averaging
  !!
  !! Averages all fine-grid cells within each coarse cell
  !! Preserves volume/mass conservation
  !! Prevents artificial deepening/shallowing
  !! Units: Same as input
  !!
  !! @param xi X-coordinates of fine grid
  !! @param yi Y-coordinates of fine grid
  !! @param Fi Field values on fine grid
  !! @param xc X-coordinates of coarse grid
  !! @param yc Y-coordinates of coarse grid
  !! @return Fc Regridded field on coarse grid
  !!
  function regrid2d_mean(xi,yi,Fi,xc,yc) result(Fc)
    real(dp), intent(in) :: xi(:,:), yi(:,:), Fi(:,:)
    real(dp), intent(in) :: xc(:,:), yc(:,:)
    real(dp), allocatable :: Fc(:,:)
    integer :: nxi, nyi, nxc, nyc
    integer :: ii, ji, ic, jc
    integer :: i_start, i_end, j_start, j_end
    integer :: xcount
    real(dp) :: x_left, x_right, y_bottom, y_top
    real(dp) :: fval, xsum

    nxi = size(xi,1); nyi = size(xi,2)
    if (any(shape(yi) /= [nxi,nyi])) call crash('regrid2d_mean - xi/yi shape mismatch')
    if (any(shape(Fi) /= [nxi,nyi])) call crash('regrid2d_mean - Fi shape mismatch')
    nxc = size(xc,1); nyc = size(xc,2)
    if (any(shape(yc) /= [nxc,nyc])) call crash('regrid2d_mean - xc/yc shape mismatch')
    allocate(Fc(nxc,nyc))

    do jc = 1, nyc
    do ic = 1, nxc
      call cell_bounds2d(xc, yc, ic, jc, x_left, x_right, y_bottom, y_top)
      call estimate_index_window(xi(:,1), x_left,  x_right,  i_start, i_end)
      call estimate_index_window(yi(1,:), y_bottom, y_top,   j_start, j_end)
      xsum = 0.0_dp; xcount = 0
      do ji = j_start, j_end
      do ii = i_start, i_end
        fval = Fi(ii,ji)
        if (.not.isnan(fval)) then
          xsum = xsum + fval
          xcount = xcount + 1
        endif
      enddo
      enddo
      if (xcount > 0) then
        Fc(ic,jc) = xsum / real(xcount, dp)
      else
        Fc(ic,jc) = ieee_value(0.0_dp, ieee_quiet_nan)
      endif
    enddo
    enddo

  end function regrid2d_mean


  !> @brief Minimum value regridding
  !!
  !! Takes the minimum (shallowest for bathymetry) value within each coarse cell
  !! Conservative for safety-critical applications (navigation, obstacle avoidance)
  !! Prevents ships from running aground
  !!
  !! @param xi X-coordinates of fine grid
  !! @param yi Y-coordinates of fine grid
  !! @param Fi Field values on fine grid
  !! @param xc X-coordinates of coarse grid
  !! @param yc Y-coordinates of coarse grid
  !! @return Fc Regridded field on coarse grid
  !!
  function regrid2d_min(xi,yi,Fi,xc,yc) result(Fc)
    real(dp), intent(in) :: xi(:,:), yi(:,:), Fi(:,:)
    real(dp), intent(in) :: xc(:,:), yc(:,:)
    real(dp), allocatable :: Fc(:,:)
    integer :: nxi, nyi, nxc, nyc
    integer :: ii, ji, ic, jc
    integer :: i_start, i_end, j_start, j_end
    integer :: xcount
    real(dp) :: x_left, x_right, y_bottom, y_top
    real(dp) :: fval, xmin

    nxi = size(xi,1); nyi = size(xi,2)
    if (any(shape(yi) /= [nxi,nyi])) call crash('regrid2d_min - xi/yi shape mismatch')
    if (any(shape(Fi) /= [nxi,nyi])) call crash('regrid2d_min - Fi shape mismatch')
    nxc = size(xc,1); nyc = size(xc,2)
    if (any(shape(yc) /= [nxc,nyc])) call crash('regrid2d_min - xc/yc shape mismatch')
    allocate(Fc(nxc,nyc))

    do jc = 1, nyc
    do ic = 1, nxc
      call cell_bounds2d(xc, yc, ic, jc, x_left, x_right, y_bottom, y_top)
      call estimate_index_window(xi(:,1), x_left,  x_right,  i_start, i_end)
      call estimate_index_window(yi(1,:), y_bottom, y_top,   j_start, j_end)
      xmin = huge(1.0_dp); xcount = 0
      do ji = j_start, j_end
      do ii = i_start, i_end
        fval = Fi(ii,ji)
        if (.not.isnan(fval)) then
          xmin = min(xmin, fval)
          xcount = xcount + 1
        endif
      enddo
      enddo
      if (xcount > 0) then
        Fc(ic,jc) = xmin
      else
        Fc(ic,jc) = ieee_value(0.0_dp, ieee_quiet_nan)
      endif
    enddo
    enddo

  end function regrid2d_min


  !> @brief Maximum value regridding
  !!
  !! Takes the maximum (deepest for bathymetry) value within each coarse cell
  !! Used when deep features are critical for ocean models
  !!
  !! @param xi X-coordinates of fine grid
  !! @param yi Y-coordinates of fine grid
  !! @param Fi Field values on fine grid
  !! @param xc X-coordinates of coarse grid
  !! @param yc Y-coordinates of coarse grid
  !! @return Fc Regridded field on coarse grid
  !!
  function regrid2d_max(xi,yi,Fi,xc,yc) result(Fc)
    real(dp), intent(in) :: xi(:,:), yi(:,:), Fi(:,:)
    real(dp), intent(in) :: xc(:,:), yc(:,:)
    real(dp), allocatable :: Fc(:,:)
    integer :: nxi, nyi, nxc, nyc
    integer :: ii, ji, ic, jc
    integer :: i_start, i_end, j_start, j_end
    integer :: xcount
    real(dp) :: x_left, x_right, y_bottom, y_top
    real(dp) :: fval, xmax

    nxi = size(xi,1); nyi = size(xi,2)
    if (any(shape(yi) /= [nxi,nyi])) call crash('regrid2d_max - xi/yi shape mismatch')
    if (any(shape(Fi) /= [nxi,nyi])) call crash('regrid2d_max - Fi shape mismatch')
    nxc = size(xc,1); nyc = size(xc,2)
    if (any(shape(yc) /= [nxc,nyc])) call crash('regrid2d_max - xc/yc shape mismatch')
    allocate(Fc(nxc,nyc))

    do jc = 1, nyc
    do ic = 1, nxc
      call cell_bounds2d(xc, yc, ic, jc, x_left, x_right, y_bottom, y_top)
      call estimate_index_window(xi(:,1), x_left,  x_right,  i_start, i_end)
      call estimate_index_window(yi(1,:), y_bottom, y_top,   j_start, j_end)
      xmax = -huge(1.0_dp); xcount = 0
      do ji = j_start, j_end
      do ii = i_start, i_end
        fval = Fi(ii,ji)
        if (.not.isnan(fval)) then
          xmax = max(xmax, fval)
          xcount = xcount + 1
        endif
      enddo
      enddo
      if (xcount > 0) then
        Fc(ic,jc) = xmax
      else
        Fc(ic,jc) = ieee_value(0.0_dp, ieee_quiet_nan)
      endif
    enddo
    enddo

  end function regrid2d_max


  !> @brief Median value regridding
  !!
  !! Takes the median value within each coarse cell
  !! More robust to outliers than mean
  !! Good compromise between smoothing and feature preservation
  !!
  !! @param xi X-coordinates of fine grid
  !! @param yi Y-coordinates of fine grid
  !! @param Fi Field values on fine grid
  !! @param xc X-coordinates of coarse grid
  !! @param yc Y-coordinates of coarse grid
  !! @return Fc Regridded field on coarse grid
  !!
  function regrid2d_median(xi,yi,Fi,xc,yc) result(Fc)
    real(dp), intent(in) :: xi(:,:), yi(:,:), Fi(:,:)
    real(dp), intent(in) :: xc(:,:), yc(:,:)
    real(dp), allocatable :: Fc(:,:)
    integer :: nxi, nyi, nxc, nyc
    integer :: ii, ji, ic, jc
    integer :: i_start, i_end, j_start, j_end
    integer :: xcount, max_points, mid
    real(dp) :: x_left, x_right, y_bottom, y_top
    real(dp) :: fval
    real(dp), allocatable :: values(:)

    nxi = size(xi,1); nyi = size(xi,2)
    if (any(shape(yi) /= [nxi,nyi])) call crash('regrid2d_median - xi/yi shape mismatch')
    if (any(shape(Fi) /= [nxi,nyi])) call crash('regrid2d_median - Fi shape mismatch')
    nxc = size(xc,1); nyc = size(xc,2)
    if (any(shape(yc) /= [nxc,nyc])) call crash('regrid2d_median - xc/yc shape mismatch')
    allocate(Fc(nxc,nyc))
  
    max_points = ceiling(real(nxi,dp)/real(nxc,dp)) * ceiling(real(nyi,dp)/real(nyc,dp)) * 4
    allocate(values(max_points))

    do jc = 1, nyc
    do ic = 1, nxc
      call cell_bounds2d(xc, yc, ic, jc, x_left, x_right, y_bottom, y_top)
      call estimate_index_window(xi(:,1), x_left,  x_right,  i_start, i_end)
      call estimate_index_window(yi(1,:), y_bottom, y_top,   j_start, j_end)
      xcount = 0
      do ji = j_start, j_end
      do ii = i_start, i_end
        fval = Fi(ii,ji)
        if (.not.isnan(fval)) then
          xcount = xcount + 1
          values(xcount) = fval
        endif
      enddo
      enddo
      if (xcount > 0) then
        call quicksort(values(1:xcount))
        mid = xcount / 2
        if (mod(xcount,2) == 0) then
          Fc(ic,jc) = 0.5_dp*(values(mid) + values(mid+1))
        else
          Fc(ic,jc) = values(mid+1)
        endif
      else
        Fc(ic,jc) = ieee_value(0.0_dp, ieee_quiet_nan)
      endif
    enddo
    enddo
    deallocate(values)

  end function regrid2d_median


  !> @brief Area-weighted averaging regridding
  !!
  !! Weights fine-grid cells by their overlap area with coarse cells
  !! Most accurate for non-uniform grids or when cells don't align
  !! Preserves integral quantities (mass, volume, energy)
  !!
  !! @param xi X-coordinates of fine grid (cell centers)
  !! @param yi Y-coordinates of fine grid (cell centers)
  !! @param Fi Field values on fine grid
  !! @param xc X-coordinates of coarse grid (cell centers)
  !! @param yc Y-coordinates of coarse grid (cell centers)
  !! @return Fc Regridded field on coarse grid
  !!
  function regrid2d_area_weighted(xi,yi,Fi,xc,yc) result(Fc)
    real(dp), intent(in) :: xi(:,:), yi(:,:), Fi(:,:)
    real(dp), intent(in) :: xc(:,:), yc(:,:)
    real(dp), allocatable :: Fc(:,:)
    integer :: nxi, nyi, nxc, nyc
    integer :: ii, ji, ic, jc
    integer :: i_start, i_end, j_start, j_end
    real(dp) :: xc_left, xc_right, yc_bottom, yc_top
    real(dp) :: fval
    real(dp) :: xi_left, xi_right, yi_bottom, yi_top
    real(dp) :: overlap_x, overlap_y, overlap_area
    real(dp) :: weighted_sum, total_area

    nxi = size(xi,1); nyi = size(xi,2)
    if (any(shape(yi) /= [nxi,nyi])) call crash('regrid2d_area_weighted - xi/yi shape mismatch')
    if (any(shape(Fi) /= [nxi,nyi])) call crash('regrid2d_area_weighted - Fi shape mismatch')
    nxc = size(xc,1); nyc = size(xc,2)
    if (any(shape(yc) /= [nxc,nyc])) call crash('regrid2d_area_weighted - xc/yc shape mismatch')
    allocate(Fc(nxc,nyc))

    do jc = 1, nyc
    do ic = 1, nxc
      call cell_bounds2d(xc, yc, ic, jc, xc_left, xc_right, yc_bottom, yc_top)
      call estimate_index_window(xi(:,1), xc_left,  xc_right,  i_start, i_end)
      call estimate_index_window(yi(1,:), yc_bottom, yc_top,   j_start, j_end)
      weighted_sum = 0.0_dp
      total_area   = 0.0_dp
      do ji = j_start, j_end
      do ii = i_start, i_end
        fval = Fi(ii,ji)
        if (isnan(fval)) cycle
        call cell_bounds2d(xi, yi, ii, ji, xi_left, xi_right, yi_bottom, yi_top)
        overlap_x = min(xc_right, xi_right) - max(xc_left, xi_left)
        overlap_y = min(yc_top,  yi_top)    - max(yc_bottom, yi_bottom)
        if (overlap_x > 0.0_dp .and. overlap_y > 0.0_dp) then
          overlap_area = overlap_x * overlap_y
          weighted_sum = weighted_sum + fval * overlap_area
          total_area   = total_area   + overlap_area
        endif
      enddo
      enddo
      if (total_area > 0.0_dp) then
        Fc(ic,jc) = weighted_sum / total_area
      else
        Fc(ic,jc) = ieee_value(0.0_dp, ieee_quiet_nan)
      endif
    enddo
    enddo

  end function regrid2d_area_weighted


  ! Monotonized Central limiter: returns limited slope given left, center, right
  pure elemental real(dp) function mc_limiter(fl, fc, fr) result(s)
    real(dp), intent(in) :: fl, fc, fr
    real(dp) :: dl, dr, dm
    dl = fc - fl
    dr = fr - fc
    if (dl*dr <= 0.0_dp) then
      s = 0.0_dp
    else
      dm = 0.5_dp*(fr - fl)
      s = sign(1.0_dp, dm) * min( abs(dm), 2.0_dp*abs(dl), 2.0_dp*abs(dr) )
    endif
  end function mc_limiter


  ! Integrate bilinear reconstruction over overlap rectangle
  pure function integrate_bilinear_over_rect(f0, x0, y0, sx, sy, xa, xb, ya, yb) result(I)
    real(dp), intent(in) :: f0, x0, y0, sx, sy
    real(dp), intent(in) :: xa, xb, ya, yb  ! rectangle bounds (xa<xb, ya<yb)
    real(dp) :: I, Ax, Ay, cx, cy
    Ax = xb - xa
    Ay = yb - ya
    ! Integral of f(x,y) = f0 + sx*(x-x0) + sy*(y-y0) over [xa,xb]x[ya,yb]
    ! = f0*Ax*Ay + sx*( (xb-x0+xa-x0)/2 * Ay * Ax )? Careful: separate integrals
    ! Do it explicitly:
    ! ∫∫ f0 dA = f0*Ax*Ay
    ! ∫∫ sx*(x-x0) dA = sx * [ ∫(x-x0) dx from xa..xb ] * [ ∫ dy from ya..yb ]
    !  = sx * (0.5*( (xb-x0)^2 - (xa-x0)^2 )) * (Ay)
    ! ∫∫ sy*(y-y0) dA = sy * (0.5*( (yb-y0)^2 - (ya-y0)^2 )) * (Ax)
    cx = 0.5_dp*((xb - x0)**2 - (xa - x0)**2)
    cy = 0.5_dp*((yb - y0)**2 - (ya - y0)**2)
    I = f0*Ax*Ay + sx*cx*Ay + sy*cy*Ax
  end function integrate_bilinear_over_rect


  function regrid2d_sl_bilinear(xi, yi, Fi, xc, yc) result(Fc)
    real(dp), intent(in) :: xi(:,:), yi(:,:), Fi(:,:)
    real(dp), intent(in) :: xc(:,:), yc(:,:)
    real(dp), allocatable :: Fc(:,:)

    integer :: nxi, nyi, nxc, nyc
    integer :: ic, jc, ii, ji
    integer :: i_start, i_end, j_start, j_end
    real(dp) :: xc_left, xc_right, yc_bottom, yc_top
    real(dp) :: xi_left, xi_right, yi_bottom, yi_top
    real(dp) :: xa, xb, ya, yb, Ax, Ay
    real(dp) :: area_total, integral_sum
    real(dp) :: f0, sx_loc, sy_loc

    nxi = size(xi,1); nyi = size(xi,2)
    if (any(shape(yi) /= [nxi,nyi])) call crash('regrid2d_sl_bilinear - xi/yi shape mismatch')
    if (any(shape(Fi) /= [nxi,nyi])) call crash('regrid2d_sl_bilinear - Fi shape mismatch')

    nxc = size(xc,1); nyc = size(xc,2)
    if (any(shape(yc) /= [nxc,nyc])) call crash('regrid2d_sl_bilinear - xc/yc shape mismatch')

    allocate(Fc(nxc,nyc))

    do jc = 1, nyc
    do ic = 1, nxc

      ! Coarse cell bounds via helper (consistent with other methods)
      call cell_bounds2d(xc, yc, ic, jc, xc_left, xc_right, yc_bottom, yc_top)

      ! Restrict fine indices to potential overlaps
      call estimate_index_window(xi(:,1), xc_left,  xc_right,  i_start, i_end)
      call estimate_index_window(yi(1,:), yc_bottom, yc_top,   j_start, j_end)

      area_total   = 0.0_dp
      integral_sum = 0.0_dp

      do ji = j_start, j_end
        ! Fine cell y-bounds
        call cell_bounds2d(xi, yi, max(1,min(nxi, i_start)), ji, xi_left, xi_right, yi_bottom, yi_top)
        ! The above call used i_start only to pass a valid i; recompute y-bounds precisely for this ji:
        call cell_bounds2d(xi, yi, 1, ji, xi_left, xi_right, yi_bottom, yi_top)

        do ii = i_start, i_end
          f0 = Fi(ii,ji)
          if (isnan(f0)) cycle

          ! Fine cell exact bounds
          call cell_bounds2d(xi, yi, ii, ji, xi_left, xi_right, yi_bottom, yi_top)

          ! Overlap rectangle
          xa = max(xc_left,  xi_left)
          xb = min(xc_right, xi_right)
          ya = max(yc_bottom, yi_bottom)
          yb = min(yc_top,    yi_top)

          Ax = xb - xa
          Ay = yb - ya
          if (Ax > 0.0_dp .and. Ay > 0.0_dp) then
            sx_loc = slope_x(ii,ji)
            sy_loc = slope_y(ii,ji)
            integral_sum = integral_sum + integrate_bilinear_over_rect( &
                            f0, xi(ii,ji), yi(ii,ji), sx_loc, sy_loc, xa, xb, ya, yb)
            area_total   = area_total   + Ax*Ay
          endif
        enddo
      enddo
      if (area_total > 0.0_dp) then
        Fc(ic,jc) = integral_sum / area_total
      else
        Fc(ic,jc) = ieee_value(0.0_dp, ieee_quiet_nan)
      end if
    enddo
    enddo

    contains

    ! Compute limited slopes at (i,j) using neighbors if available  
      pure function slope_x(i,j) result(sx)  
        integer, intent(in) :: i, j  
        real(dp) :: sx  
        real(dp) :: fl, fc, fr, s_raw, dxl, dxr  
        sx = 0.0_dp  
        if (i <= 1 .or. i >= nxi) return  
        fl = Fi(i-1,j); fc = Fi(i,j); fr = Fi(i+1,j)  
        if (isnan(fl) .or. isnan(fc) .or. isnan(fr)) return  
        s_raw = mc_limiter(fl, fc, fr)  
        dxl = abs(xi(i,j) - xi(i-1,j))  
        dxr = abs(xi(i+1,j) - xi(i,j))  
        if (dxl>0.0_dp .and. dxr>0.0_dp) sx = s_raw / (0.5_dp*(dxl+dxr))  
      end function slope_x  
  
      pure function slope_y(i,j) result(sy)  
        integer, intent(in) :: i, j  
        real(dp) :: sy  
        real(dp) :: fl, fc, fr, s_raw, dyl, dyr  
        sy = 0.0_dp  
        if (j <= 1 .or. j >= nyi) return  
        fl = Fi(i,j-1); fc = Fi(i,j); fr = Fi(i,j+1)  
        if (isnan(fl) .or. isnan(fc) .or. isnan(fr)) return  
        s_raw = mc_limiter(fl, fc, fr)  
        dyl = abs(yi(i,j) - yi(i,j-1))  
        dyr = abs(yi(i,j+1) - yi(i,j))  
        if (dyl>0.0_dp .and. dyr>0.0_dp) sy = s_raw / (0.5_dp*(dyl+dyr))  
      end function slope_y  

  end function regrid2d_sl_bilinear


!  ! Apply CSR matrix: y = W * x
!  subroutine csr_matvec(W, x, y)
!    type(csr_matrix), intent(in) :: W
!    real(dp), intent(in) :: x(:)
!    real(dp), intent(out) :: y(:)
!    integer :: i, k, k0, k1, col
!
!    if (size(y) /= W%n_rows) call crash('csr_matvec: y size mismatch')  
!    if (size(x) /= W%n_cols) call crash('csr_matvec: x size mismatch')  
!    y = 0.0_dp
!    do i = 1, W%n_rows
!    k0 = W%row_ptr(i)
!    k1 = W%row_ptr(i+1)-1
!    do k = k0, k1
!    col = W%col_ind(k)  ! 1-based
!    y(i) = y(i) + W%values(k) * x(col)
!    enddo
!    enddo
!  end subroutine csr_matvec


!!> @brief Learned regridding using pre-trained neural network weights
!!!
!!! Applies a sparse weight matrix (CSR format) to regrid from fine to coarse grid.
!!! The weights are typically learned from training data using neural networks or
!!! other machine learning methods to preserve specific features or minimize error.
!!!
!!! @param xi X-coordinates of fine grid
!!! @param yi Y-coordinates of fine grid
!!! @param Fi Field values on fine grid
!!! @param xc X-coordinates of coarse grid
!!! @param yc Y-coordinates of coarse grid
!!! @param W CSR sparse weight matrix (n_rows=nxc*nyc, n_cols=nxi*nyi)
!!! @return Fc Regridded field on coarse grid
!!!
!!! HOW TO GENERATE WEIGHTS USING PYTHON:
!!! ----
!!! 1. Train a neural network or compute analytical weights in Python
!!! 2. Export to CSR format and save to NetCDF:
!!!
!!!    import numpy as np
!!!    import scipy.sparse as sp
!!!    from netCDF4 import Dataset
!!!
!!!    # Example: Create weight matrix W (nxc*nyc rows, nxi*nyi cols)
!!!    # W[i_coarse, j_fine] = weight for contribution of fine cell j to coarse cell i
!!!    W_dense = ...  # your trained weights or analytical formula
!!!    W_csr = sp.csr_matrix(W_dense)
!!!
!!!    # Save to NetCDF
!!!    with Dataset('regrid_weights.nc', 'w') as nc:
!!!    nc.createDimension('nnz', W_csr.nnz)
!!!    nc.createDimension('n_rows_plus_1', W_csr.shape[0] + 1)
!!!    
!!!    nc.createVariable('n_rows', 'i4')
!!!    nc.createVariable('n_cols', 'i4')
!!!    nc.createVariable('row_ptr', 'i4', ('n_rows_plus_1',))
!!!    nc.createVariable('col_ind', 'i4', ('nnz',))
!!!    nc.createVariable('values', 'f8', ('nnz',))
!!!    
!!!    nc.variables['n_rows'][:] = W_csr.shape[0]
!!!    nc.variables['n_cols'][:] = W_csr.shape[1]
!!!    nc.variables['row_ptr'][:] = W_csr.indptr  # 0-based cumulative
!!!    nc.variables['col_ind'][:] = W_csr.indices + 1  # Convert to 1-based for Fortran
!!!    nc.variables['values'][:] = W_csr.data
!!!
!!! 3. Load in Fortran using NetCDF library and populate csr_matrix type
!!!
!function regrid2d_learned(xi, yi, Fi, xc, yc, W) result(Fc)
!  real(dp), intent(in) :: xi(:), yi(:), Fi(:,:)
!  real(dp), intent(in) :: xc(:), yc(:)
!  type(csr_matrix), intent(in) :: W
!  real(dp), allocatable :: Fc(:,:)
!
!  integer :: nxi, nyi, nxc, nyc
!  real(dp), allocatable :: x_flat(:), y_flat(:)
!  integer :: i, j, idx
!
!  nxi = size(xi); nyi = size(yi)
!  if (size(Fi,1) /= nxi .or. size(Fi,2) /= nyi) call crash('regrid2d_learned - Incompatible Fi size')
!  nxc = size(xc); nyc = size(yc)
!
!  if (W%n_cols /= nxi*nyi) call crash('regrid2d_learned - W n_cols mismatch with source size')
!  if (W%n_rows /= nxc*nyc) call crash('regrid2d_learned - W n_rows mismatch with target size')
!
!  allocate(Fc(nxc,nyc))
!  allocate(x_flat(nxi*nyi), y_flat(nxc*nyc))
!
!  ! Flatten Fi in Fortran order: (i fastest)
!  idx = 0
!  do j = 1, nyi
!    do i = 1, nxi
!    idx = idx + 1
!    if (isnan(Fi(i,j))) then
!    x_flat(idx) = 0.0_dp
!    else
!    x_flat(idx) = Fi(i,j)
!    endif
!    enddo
!  enddo
!
!  call csr_matvec(W, x_flat, y_flat)
!
!  ! Reshape to Fc; if weights excluded missing properly, y_flat will be valid.
!  ! If some target rows received no contributions, ensure they are set to fill_value.
!  ! We detect it if the row_ptr segment is empty or if the sum of weights is zero.
!  ! Here we assume weights already normalized; otherwise include a sum-of-weights vector.
!  idx = 0
!  do j = 1, nyc
!    do i = 1, nxc
!    idx = idx + 1
!    ! Check bounds and empty row
!    if (idx > W%n_rows .or. idx > size(W%row_ptr)-1) then
!    Fc(i,j) = ieee_value(0.0_dp, ieee_quiet_nan) 
!    else if (W%row_ptr(idx+1) == W%row_ptr(idx)) then
!    Fc(i,j) = ieee_value(0.0_dp, ieee_quiet_nan)
!    else
!    Fc(i,j) = y_flat(idx)
!    endif
!    enddo
!  enddo
!
!  deallocate(x_flat, y_flat)
!end function regrid2d_learned
!
!!

end module module_regrid
