! ======================================================================== !
! LINEAR ALGEBRA MODULE                                                    !
! Based on ALOGES PROJECT by Quim Ballabrera                               !
! Last Modified: 2024-01-15                                                !
!                                                                          !
! Copyright (C) 2022, Joaquim Ballabrera                                   !
! This program is free software: you can redistribute it and/or modify     !
! it under the terms of the GNU Lesser General Public License as published !
! by the Free Software Foundation, either version 3 of the License, or     !
! (at your option) any later version.                                      !
! ======================================================================== !

module module_linalg

use, intrinsic :: IEEE_ARITHMETIC, ONLY : IEEE_VALUE, IEEE_QUIET_NAN
use module_types
use module_constants
use module_math

implicit none

private
public :: matrinv, ludcmp, lubksb, det, trace, &
          out_product, norm, cross_product, identity, &
          matrix_multiply_blocked, solve_linear_system, expm, kron

interface matrix_multiply
  module procedure matrix_multiply_blocked
end interface matrix_multiply

contains

! ===================================================================
! ===================================================================
! ...
  pure function out_product(A,B) result(C)
    real(dp), dimension(:), intent(in) :: A, B
    real(dp), dimension(size(A),size(B)) :: C
    integer :: n,m
    
    n = size(A)
    m = size(B)
    C = spread(A,dim=2,ncopies=m) * spread(B,dim=1,ncopies=n)
  end function out_product
! ...
! ===================================================================
! ...
  subroutine ludcmp(a,indx,d)
    real(dp), dimension(:,:), intent(inout) :: a
    integer, dimension(:), intent(out) :: indx
    real(dp), intent(out) :: d

    real(dp), parameter :: TINY=1.0D-20
    integer :: n,i,imax,j,k
    real(dp) :: vv(size(a,1))

    n = size(a,1)
    if (size(a,2) /= n) then
      error stop "ludcmp: Input matrix must be square"
    endif
    
    d = 1.d0
    vv = maxval(abs(a),dim=2)
    if (any(vv.EQ.0.0)) stop 'singular matrix in ludcmp'
    vv(:) = 1.d0/vv(:)

    do j=1,n
      imax = (j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))
      if (j /= imax) then
        call swap(a(imax,:),a(j,:))
        d = -d
        vv(imax) = vv(j)
      endif
      indx(j) = imax
      if (a(j,j).EQ.0.0) a(j,j) = TINY
      a(j+1:n,j) = a(j+1:n,j)/a(j,j)
      a(j+1:n,j+1:n) = a(j+1:n,j+1:n) - out_product(a(j+1:n,j),a(j,j+1:n))
    enddo
  end subroutine ludcmp
! ...
! ===================================================================
! ...
  subroutine lubksb(a,indx,b)
    real(dp), dimension(:,:), intent(in) :: a
    integer, dimension(:), intent(in) :: indx
    real(dp), dimension(:), intent(inout) :: b

    integer :: n,i,ii,ll
    real(dp) :: summ

    n = size(a,1)
    if (size(b) /= n) then
      error stop "lubksb: Vector size doesn't match matrix"
    endif

    ii = 0
    do i=1,n
      ll    = indx(i)
      summ  = b(ll)
      b(ll) = b(i)
      if (ii.ne.0) then
        summ = summ - dot_product(a(i,ii:i-1),b(ii:i-1))
      else if (summ.ne.0.0) then
        ii = i
      endif
      b(i) = summ
    enddo
    do i=n,1,-1
      b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n))) / a(i,i)
    enddo
  end subroutine lubksb
! ...
! ===================================================================
! ...
  function matrinv(A) result(B)
    real(dp), dimension(:,:), intent(in) :: A
    real(dp), dimension(size(A,1),size(A,2)) :: B

    integer :: n,i,j
    integer, dimension(size(A,1)) :: indx
    real(dp) :: d
    real(dp), dimension(size(A,1),size(A,2)) :: C

    n = size(A,1)
    if (size(A,2).ne.n) stop 'ERROR in matrinv: Input array not square'

    C(:,:) = A(:,:)
    B(:,:) = 0.0D0
    do i=1,n
      B(i,i) = 1.0D0
    ENDDO

    CALL ludcmp(C,indx,d)
    DO j=1,N
      CALL lubksb(C,indx,B(:,j))
    ENDDO
  end function matrinv
! ...
! ===================================================================
! ...
  function det(A) result(d)
    real(dp), dimension(:,:), intent(in) :: A
    real(dp) :: d
    integer :: n, i
    integer, dimension(:), allocatable :: indx
    real(dp), dimension(:,:), allocatable :: LU
    real(dp) :: d_sign

    n = size(A, 1)
    if (size(A, 2) /= n) then
      d = ieee_value(1.0_dp, ieee_quiet_nan)
      return
    endif
    
    allocate(LU(n,n), indx(n))
    LU = A
    call ludcmp(LU, indx, d_sign)
    
    d = d_sign
    do i = 1, n
      d = d * LU(i,i)
    end do
    
    deallocate(LU, indx)
  end function det
! ...
! ===================================================================
! ...
  function trace(A) result(tr)
    real(dp), dimension(:,:), intent(in) :: A
    real(dp) :: tr
    integer :: i, n
    
    n = min(size(A,1), size(A,2))
    tr = 0.0_dp
    do i = 1, n
      tr = tr + A(i,i)
    end do
  end function trace
! ...
! ===================================================================
! ...
  function norm(x, p) result(nrm)
    real(dp), dimension(:), intent(in) :: x
    integer, intent(in), optional :: p
    real(dp) :: nrm
    integer :: norm_type
    
    norm_type = 2
    if (present(p)) norm_type = p
    
    select case (norm_type)
    case (1)
      nrm = sum(abs(x))
    case (2)
      nrm = sqrt(sum(x**2))
    case (:0)
      nrm = maxval(abs(x))
    case default
      nrm = (sum(abs(x)**norm_type))**(1.0_dp/norm_type)
    end select
  end function norm
! ...
! ===================================================================
! ...
  function cross_product(a, b) result(c)
    real(dp), dimension(3), intent(in) :: a, b
    real(dp), dimension(3) :: c
    
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross_product
! ...
! ===================================================================
! ...
  function identity(n) result(I_mat)
    integer, intent(in) :: n
    real(dp), dimension(n,n) :: I_mat
    integer :: i
    I_mat = 0.0_dp
    do i = 1, n
      I_mat(i,i) = 1.0_dp
    end do
  end function identity
! ...
! ===================================================================
! ...
  function matrix_multiply_blocked(A, B, block_size) result(C)
    real(dp), dimension(:,:), intent(in) :: A, B
    integer, intent(in), optional :: block_size
    real(dp), dimension(size(A,1), size(B,2)) :: C
    
    integer :: m, n, p, i, j, k, ii, jj, kk, bs
    integer :: i_end, j_end, k_end
    
    m = size(A,1)
    n = size(A,2)
    p = size(B,2)
    
    if (size(B,1) /= n) then
      error stop "matrix_multiply: Dimension mismatch"
    endif
    
    bs = 32  ! Default block size
    if (present(block_size)) bs = block_size
    
    C = 0.0_dp
    
    ! Blocked matrix multiplication for better cache performance
    do ii = 1, m, bs
      i_end = min(ii + bs - 1, m)
      do jj = 1, p, bs
        j_end = min(jj + bs - 1, p)
        do kk = 1, n, bs
          k_end = min(kk + bs - 1, n)
          
          do i = ii, i_end
            do k = kk, k_end
              do j = jj, j_end
                C(i,j) = C(i,j) + A(i,k) * B(k,j)
              end do
            end do
          end do
          
        end do
      end do
    end do
  end function matrix_multiply_blocked
! ...
! ===================================================================
! ...
  function solve_linear_system(A, b) result(x)
    real(dp), dimension(:,:), intent(in) :: A
    real(dp), dimension(:), intent(in) :: b
    real(dp), dimension(size(b)) :: x
    
    integer :: n
    integer, dimension(:), allocatable :: indx
    real(dp), dimension(:,:), allocatable :: LU
    real(dp) :: d
    
    n = size(A,1)
    if (size(A,2) /= n .or. size(b) /= n) then
      error stop "solve_linear_system: Dimension mismatch"
    endif
    
    allocate(LU(n,n), indx(n))
    LU = A
    x = b
    
    call ludcmp(LU, indx, d)
    call lubksb(LU, indx, x)
    
    deallocate(LU, indx)
  end function solve_linear_system
! ...
! ===================================================================
! ...
  function expm(A) result(B)
    ! Matrix exponential using Pade approximation
    real(dp), dimension(:,:), intent(in) :: A
    real(dp), dimension(size(A,1),size(A,2)) :: B
    integer :: n, i, j, q 
    real(dp), dimension(:,:), allocatable :: A2, A4, A6, U, V, N_mat, D_mat
    real(dp) :: normA, c
    
    n = size(A,1)
    if (size(A,2) /= n) error stop "expm: Input matrix must be square" 
    
    normA = maxval(sum(abs(A), dim=2)) 
    q = max(0, ceiling(log(normA)/log(2.0_dp)) + 1)
    
    allocate(A2(n,n), A4(n,n), A6(n,n), U(n,n), V(n,n), N_mat(n,n), D_mat(n,n))
    
    A2 = matmul(A, A) / (2.0_dp**q)
    A4 = matmul(A2, A2)
    A6 = matmul(A2, A4)
    
    ! Pade approximation (6,6)
    U = A6 + 16380.0_dp * A4 + 40840800.0_dp * A2
    U = matmul(A, U + 33522128640.0_dp * identity(n))
    
    V = 136.0_dp * A6 + 24504480.0_dp * A4 + 32590958400.0_dp * A2
    V = V + 10559470521600.0_dp * identity(n)
    
    N_mat = V + U 
    D_mat = V - U 
    
    B = matrinv(D_mat)
    B = matmul(N_mat, B)
    
    ! Repeated squaring
    do i = 1, q
      B = matmul(B, B)
    end do  
    
    deallocate(A2, A4, A6, U, V, N_mat, D_mat)
  end function expm
! ...
! ===================================================================
! ...
  function kron(A, B) result(C)
    ! Kronecker product 
    real(dp), dimension(:,:), intent(in) :: A, B 
    real(dp), dimension(:,:), allocatable :: C
    integer :: i, j, m, n, p, q 
    
    m = size(A,1); n = size(A,2)
    p = size(B,1); q = size(B,2)
    
    allocate(C(m*p, n*q))
    
    do concurrent (i=1:m, j=1:n)
      C((i-1)*p+1:i*p, (j-1)*q+1:j*q) = A(i,j) * B
    end do
  end function kron
! ...
! ===================================================================
! ...
end module module_linalg
