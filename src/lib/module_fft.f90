! ====                                                                     !
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
!>
!> @file module_fft.f90
!> @brief Fast Fourier Transform utilities (1D and 2D wrappers).
!>
!> Provides an in-place radix-2 1D complex-split FFT (real/imag arrays)
!> and simple 2D FFT wrappers built on top of the 1D routine.
!>
!> Notes
!> - The 1D FFT is unnormalized. The inverse requires 1/N scaling by caller.
!> - The 2D inverse here normalizes by 1/(nx*ny).
!> - Inputs must be power-of-two sized for the radix-2 implementation.
!>
!> @author Joaquim Ballabrera
!> @date 2025-10-23

module module_fft

use module_types, only : dp
use module_constants, only : pi, two_pi, nan
use module_tools, only: crash
use module_math, only: next_power_of_2

implicit none (type, external)
private
    
public :: fft1d
public :: fft2d_forward, fft2d_inverse

contains

  !> @brief In-place radix-2 complex-split 1D FFT.
  !>
  !> Computes the discrete Fourier transform (forward if `signo=-1`,
  !> inverse if `signo=+1`) of a complex sequence stored in split form
  !> `funcionR(:)` and `funcionI(:)` with indices [0..dimx-1].
  !>
  !> @param[in]  dimx   Integer length of the transform, must be power of two.
  !> @param[inout] funcionR Real(dp) array(0:dimx-1): real part, in-place.
  !> @param[inout] funcionI Real(dp) array(0:dimx-1): imag part, in-place.
  !> @param[in]  signo  Integer: -1 for forward, +1 for inverse.
  !>
  !> @note No normalization is applied. Caller should divide by dimx for inverse.
  !> @note Uses iterative bit-reversal and Danielson–Lanczos butterfly.
  subroutine fft1d (dimx,funcionR,funcionI,signo)
    integer, intent(in)    :: dimx,signo
    real(dp), dimension(0:dimx-1), intent(inout) :: funcionR,funcionI
    real(dp)    :: tempR,tempI, wpasoR,wpasoI, wwR,wwI, ttR,ttI
    integer    :: ix,je,mm,mmax,istep
    tempR = 0.0D0
    tempI = 0.0D0
    je = 1
    do ix=0,dimx-1
      if (je.GT.ix+1) then
        tempR = funcionR(je-1)
        tempI = funcionI(je-1)
        funcionR(je-1) = funcionR(ix)
        funcionI(je-1) = funcionI(ix)
        funcionR(ix) = tempR
        funcionI(ix) = tempI
      endif
      mm = dimx/2
      do while (mm.GT.1 .AND. je.GT.mm)
        je = je - mm
        mm = mm/2
      enddo
      je = je+mm
    enddo

    mmax = 1
    do while (dimx .GT.mmax)
      istep = 2*mmax
      wpasoR = cos(PI/real(mmax,dp))
      wpasoI = signo*sin(PI/real(mmax,dp))
      wwR = 1.0_dp
      wwI = 0.0_dp
      do mm = 1,mmax
        do ix = mm-1,dimx-1,istep
          je = ix + mmax
          call c_mult(wwR,wwI,funcionR(je),funcionI(je),tempR,tempI)
          funcionR(je) = funcionR(ix) - tempR
          funcionI(je) = funcionI(ix) - tempI
          funcionR(ix) = funcionR(ix) + tempR
          funcionI(ix) = funcionI(ix) + tempI
        enddo
        call c_mult(wwR,wwI,wpasoR,wpasoI,ttR,ttI)
        wwR = ttR
        wwI = ttI
      enddo
      mmax = istep
    enddo

    contains

      subroutine c_mult (ar,ai,br,bi,cr,ci)
        real(dp), INTENT(in)    :: ar,ai,br,bi
        real(dp), INTENT(out)    :: cr,ci

        cr = ar*br - ai*bi
        ci = ar*bi + ai*br

      end subroutine c_mult

  end subroutine fft1d

  !> @brief 2D forward FFT wrapper using fft1d along each dimension.
  !>
  !> Transforms real spatial field u(nx,ny) to complex-split spectrum UfR,UfI.
  !>
  !> @param[in]  u   Real(dp) array(nx,ny)
  !> @param[out] UfR Real(dp) array(nx,ny): real part of spectrum
  !> @param[out] UfI Real(dp) array(nx,ny): imag part of spectrum
  !>
  !> @note Assumes unnormalized 1D FFT; no scaling applied here.
  !> @note Raises crash on size mismatches.
  subroutine fft2d_forward(u, UfR, UfI)
    ! In-place 2D FFT forward: real -> complex(split)
    ! u(nx,ny): real input (spatial)
    ! UfR(nx,ny), UfI(nx,ny): real arrays holding real/imag parts of spectrum

    real(dp), intent(in)  :: u(:,:)
    real(dp), intent(out) :: UfR(:,:), UfI(:,:)
    integer               :: nx, ny, i, j
    real(dp), allocatable :: tmpR(:), tmpI(:)

    nx = size(u,1); ny = size(u,2)
    if (size(UfR,1) /= nx .or. size(UfR,2) /= ny) call crash('fft2d_forward - size mismatch')
    if (size(UfI,1) /= nx .or. size(UfI,2) /= ny) call crash('fft2d_forward - size mismatch')

    ! Copy input
    UfR = u
    UfI = 0.0_dp

    ! FFT along x (first dimension) for each column j
    allocate(tmpR(0:nx-1), tmpI(0:nx-1))
    do j = 1, ny
      do i = 1, nx
        tmpR(i-1) = UfR(i,j)
        tmpI(i-1) = UfI(i,j)
      end do
      call fft1d(nx, tmpR, tmpI, -1)  ! forward
      do i = 1, nx
        UfR(i,j) = tmpR(i-1)
        UfI(i,j) = tmpI(i-1)
      end do
    end do
    deallocate(tmpR, tmpI)

    ! FFT along y (second dimension) for each row i
    allocate(tmpR(0:ny-1), tmpI(0:ny-1))
    do i = 1, nx
      do j = 1, ny
        tmpR(j-1) = UfR(i,j)
        tmpI(j-1) = UfI(i,j)
      end do
      call fft1d(ny, tmpR, tmpI, -1)  ! forward
      do j = 1, ny
        UfR(i,j) = tmpR(j-1)
        UfI(i,j) = tmpI(j-1)
      end do
    end do
    deallocate(tmpR, tmpI)
  end subroutine fft2d_forward


  !> @brief 2D inverse FFT wrapper using fft1d along each dimension.
  !>
  !> Reconstructs real spatial field u(nx,ny) from complex-split spectrum
  !> UfR,UfI. Applies normalization by 1/(nx*ny) at the end.
  !>
  !> @param[inout] UfR Real(dp) array(nx,ny): real part (overwritten in-place)
  !> @param[inout] UfI Real(dp) array(nx,ny): imag part (overwritten in-place)
  !> @param[out]   u   Real(dp) array(nx,ny): reconstructed real field
  !>
  !> @note Raises crash on size mismatches.
  subroutine fft2d_inverse(UfR, UfI, u)
    ! Inverse 2D FFT: complex(split) -> real
    ! Note: Your fft1d uses unnormalized FFT. Normalize by (nx*ny) here.
    real(dp), intent(inout) :: UfR(:,:), UfI(:,:)
    real(dp), intent(out)   :: u(:,:)
    integer :: nx, ny, i, j
    real(dp), allocatable :: tmpR(:), tmpI(:)
    real(dp) :: norm

    nx = size(UfR,1); ny = size(UfR,2)
    if (size(UfI,1) /= nx .or. size(UfI,2) /= ny) call crash('fft2d_inverse - size mismatch')
    if (size(u,1) /= nx .or. size(u,2) /= ny) call crash('fft2d_inverse - size mismatch')

    ! Inverse along y per row
    allocate(tmpR(0:ny-1), tmpI(0:ny-1))
    do i = 1, nx
      do j = 1, ny
        tmpR(j-1) = UfR(i,j)
        tmpI(j-1) = UfI(i,j)
      end do
      call fft1d(ny, tmpR, tmpI, +1)  ! inverse
      do j = 1, ny
        UfR(i,j) = tmpR(j-1)
        UfI(i,j) = tmpI(j-1)
      end do
    end do
    deallocate(tmpR, tmpI)

    ! Inverse along x per column
    allocate(tmpR(0:nx-1), tmpI(0:nx-1))
    do j = 1, ny
      do i = 1, nx
        tmpR(i-1) = UfR(i,j)
        tmpI(i-1) = UfI(i,j)
      end do
      call fft1d(nx, tmpR, tmpI, +1)  ! inverse
      do i = 1, nx
        UfR(i,j) = tmpR(i-1)
        UfI(i,j) = tmpI(i-1)
      end do
    end do
    deallocate(tmpR, tmpI)

    ! Normalize (since your fft1d does no 1/N scaling)
    norm = 1.0_dp / real(nx*ny, dp)
    do j = 1, ny
      do i = 1, nx
        u(i,j) = UfR(i,j) * norm
      end do
    end do
  end subroutine fft2d_inverse


end module module_fft
