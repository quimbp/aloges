! ==== !
! EMPIRICAL ORTHOGONAL FUNCTIONS (EOF) MODULE                    !
! Based on ALOGES PROJECT by Quim Ballabrera                    !
! Last Modified: 2024-01-15                    !
!                    !
! Copyright (C) 2022, Joaquim Ballabrera                    !
! This program is free software: you can redistribute it and/or modify     !
! it under the terms of the GNU Lesser General Public License as published !
! by the Free Software Foundation, either version 3 of the License, or     !
! (at your option) any later version.                    !
! ==== !

module module_eofs

use, intrinsic :: IEEE_ARITHMETIC, ONLY : IEEE_VALUE, IEEE_QUIET_NAN
use module_types
use module_constants, only: half_pi
use module_tools, only : argsort, lowercase
use module_linalg

implicit none (type, external)

private
public :: svdcmp, svbksb, calc_eofs, eigsort, rotate_eofs, &
          project_eofs, ROTATE, VARIMAX, vari_all, &
          normalize_eofs, check_convergence
public :: calc_cov_eigs

interface  
  subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)  
    import dp
    character :: transa, transb  
    integer :: m, n, k, lda, ldb, ldc  
    real(dp) :: alpha, beta  
    real(dp) :: a(lda,*), b(ldb,*), c(ldc,*)  
  end subroutine dgemm  
    
  subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)  
    import dp
    character :: jobu, jobvt  
    integer :: m, n, lda, ldu, ldvt, lwork, info  
    real(dp) :: a(lda,*), s(*), u(ldu,*), vt(ldvt,*), work(*)  
  end subroutine dgesvd  
    
  subroutine dsyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info)  
    import dp
    character :: jobz, uplo  
    integer :: n, lda, lwork, liwork, info  
    real(dp) :: a(lda,*), w(*), work(*)  
    integer :: iwork(*)  
  end subroutine dsyevd  
end interface  

contains

! ====
! SUBROUTINE: svdcmp
! PURPOSE: Singular Value Decomposition using Numerical Recipes algorithm
! DESCRIPTION: Decomposes matrix A = U * diag(W) * V^T
!              This is the legacy method kept for backward compatibility
! ====
subroutine svdcmp(A,U,W,V,ierr)
  real(dp), dimension(:,:), intent(in)  :: A
  real(dp), dimension(:,:), intent(out) :: U
  real(dp), dimension(:,:), intent(out) :: V
  real(dp), dimension(:), intent(out)   :: W
  integer, intent(out), optional        :: ierr

  integer :: m,n,i,j,k,l,its,nm,error_flag
  real(dp) :: anorm,c,f,g,h,s,scale,x,y,z
  real(dp), dimension(size(A,1)) :: tempm
  real(dp), dimension(size(A,2)) :: rv1,tempn
  integer, parameter :: MAX_ITERATIONS = 30

  error_flag = 0

  ! Bounds checking
  if (size(U,1) /= size(A,1) .or. size(U,2) /= size(A,2)) then
    error_flag = 2
    if (present(ierr)) then
      ierr = error_flag
      return
    else
      error stop "svdcmp: U dimensions don't match A"
    endif
  endif

  if (size(V,1) /= size(A,2) .or. size(V,2) /= size(A,2)) then
    error_flag = 3
    if (present(ierr)) then
      ierr = error_flag
      return
    else
      error stop "svdcmp: V dimensions don't match A"
    endif
  endif

  if (size(W) /= size(A,2)) then
    error_flag = 4
    if (present(ierr)) then
      ierr = error_flag
      return
    else
      error stop "svdcmp: W size doesn't match A columns"
    endif
  endif

  m = size(A,1)
  n = size(A,2)
  U = A
  g = 0.0d0
  scale = 0.0d0

  ! Householder reduction to bidiagonal form
  do i=1,n
    l = i + 1
    rv1(i) = scale*g
    g = 0.0d0
    scale = 0.0d0

    if (i.le.m) then
      scale = sum(ABS(U(i:m,i)))
      if (scale.ne.0.0d0) then
        U(i:m,i) = U(i:m,i)/scale
        s = DOT_PRODUCT(U(i:m,i),U(i:m,i))
        f = U(i,i)
        g = -sign(dsqrt(s),f)
        h = f*g - s
        U(i,i) = f - g
        tempn(l:n) = MATMUL(U(i:m,i),U(i:m,l:n))/h
        U(i:m,l:n) = U(i:m,l:n) + OUT_PRODUCT(U(i:m,i),tempn(l:n))
        U(i:m,i) = scale*U(i:m,i)
      endif
    endif

    W(i) = scale*g
    g = 0.0d0
    scale = 0.0d0

    if ((i.le.m).and.(i.ne.n)) then
      scale = SUM(ABS(U(i,l:n)))
      if (scale.ne.0.0d0) then
        U(i,l:n) = U(i,l:n)/scale
        s = DOT_PRODUCT(U(i,l:n),U(i,l:n))
        f = U(i,l)
        g = -sign(dsqrt(s),f)
        h = f*g - s
        U(i,l) = f - g
        rv1(l:n) = U(i,l:n)/h
        tempm(l:m) = MATMUL(U(l:m,l:n),U(i,l:n))
        U(l:m,l:n) = U(l:m,l:n) + OUT_PRODUCT(tempm(l:m),rv1(l:n))
        U(i,l:n) = scale*U(i,l:n)
      endif
    endif
  enddo

  anorm = maxval(ABS(w)+ABS(rv1))

  ! Accumulation of right-hand transformations
  do i=n,1,-1
    if (i.lt.n) then
      if (g.ne.0.0d0) then
        V(l:n,i) = (U(i,l:n)/U(i,l))/g
        tempn(l:n) = MATMUL(U(i,l:n),V(l:n,l:n))
        V(l:n,l:n) = V(l:n,l:n) + OUT_PRODUCT(V(l:n,i),tempn(l:n))
      endif
      V(i,l:n) = 0.0D0
      V(l:n,i) = 0.0D0
    endif
    V(i,i) = 1.0d0
    g = rv1(i)
    l = i
  enddo

  ! Accumulation of left-hand transformations
  do i=min(m,n),1,-1
    l = i + 1
    g = W(i)
    U(i,l:n) = 0.0D0
    if (g.ne.0.0d0) then
      g = 1.0D0 / g
      tempn(l:n) = g*MATMUL(U(l:m,i),U(l:m,l:n))/U(i,i)
      U(i:m,l:n) = U(i:m,l:n) + OUT_PRODUCT(U(i:m,i),tempn(l:n))
      U(i:m,i) = g*U(i:m,i)
      U(i,i) = U(i,i) + 1.0d0
    else
      U(i:m,i) = 0.0D0
      U(i,i) = 1.0d0
    endif
  enddo

  ! Diagonalization of the bidiagonal form
  do k=n,1,-1
    do its=1,MAX_ITERATIONS
      ! Test for splitting
      do l=k,1,-1
        nm = l - 1
        if ((ABS(rv1(l))+anorm).eq.anorm) exit
        if ((ABS(W(nm))+anorm).eq.anorm) then
          c = 0.0d0
          s = 1.0d0
          do i=l,k
            f = s*rv1(i)
            rv1(i) = c*rv1(i)
            if ((ABS(f)+anorm).eq.anorm) exit
            g = W(i)
            h = hypot(f,g)
            W(i) = h
            h = 1.0D0/h
            c = g*h
            s = -f*h
            tempm(1:m) = U(1:m,nm)
            U(1:m,nm) = c*U(1:m,nm) + s*U(1:m,i)
            U(1:m,i) = -s*tempm(1:m) + c*U(1:m,i)
          enddo
          exit
        endif
      enddo
      z=W(k)

      if (l.eq.k) then
        if (z.lt.0.0d0) then
          W(k) = -z
          V(1:n,k) = -V(1:n,k)
        endif
        exit
      endif
        
      ! Convergence guard
      if (its.eq.MAX_ITERATIONS) then
        if (present(ierr)) then
          ierr = 1
          return
        else
          error stop 'svdcmp: no convergence after maximum iterations'
        endif
      endif
      
      x = W(l)
      nm = k - 1
      y = W(nm)
      g = rv1(nm)
      h = rv1(k)
      f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0d0*h*y)
      g = hypot(f,1.0d0)
      f = ((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
      c = 1.0d0
      s = 1.0d0
      do j=l,nm
        i = j + 1
        g = rv1(i)
        y = W(i)
        h = s*g
        g = c*g
        z = hypot(f,h)
        rv1(j) = z
        c = f/z
        s = h/z
        f = (x*c)+(g*s)
        g = -(x*s)+(g*c)
        h = y*s
        y = y*c
        tempn(1:n) = V(1:n,j)
        V(1:n,j) = c*V(1:n,j) + s*V(1:n,i)
        V(1:n,i) = -s*tempn(1:n) + c*V(1:n,i)
        z = hypot(f,h)
        W(j) = z
        if (z.ne.0.0d0) then
          z = 1.0d0/z
          c = f*z
          s = h*z
        endif
        f = (c*g)+(s*y)
        x = -(s*g)+(c*y)
        tempm(1:m) = U(1:m,j)
        U(1:m,j) = c*U(1:m,j) + s*U(1:m,i)
        U(1:m,i) = -s*tempm(1:m) + c*U(1:m,i)
      enddo
      rv1(l)=0.0d0
      rv1(k)=f
      W(k)=x
    enddo
  enddo

  if (present(ierr)) ierr = error_flag

end subroutine svdcmp

! ====
! FUNCTION: svbksb
! PURPOSE: Back-substitution for SVD solution of linear system
! BUG FIX: tmp array sizing corrected from size(U,1) to size(U,2)
! ====
pure function svbksb(U,W,V,b) result(x)
  real(dp), dimension(:,:), intent(in) :: U,V
  real(dp), dimension(:), intent(in) :: W,b
  real(dp), dimension(size(V,2)) :: x

  integer :: m,n,i
  real(dp) :: threshold
  real(dp), dimension(size(U,2)) :: tmp  ! FIX: was size(U,1), should be size(U,2)

  m = size(U,1)
  n = size(U,2)

  ! Set a threshold for singular values (machine epsilon based)
  threshold = max(m, n) * epsilon(1.0_dp) * maxval(W)

  tmp = 0.0_dp
  do i = 1, n
    if (W(i) > threshold) then
      tmp(i) = dot_product(b, U(:,i)) / W(i)
    endif
  enddo

  x = matmul(V, tmp)

end function svbksb

! ====
! SUBROUTINE: calc_eofs
! PURPOSE: Calculate EOFs using specified method
! METHODS:
!   'legacy'      - Numerical Recipes SVD (works for any nx, nt)
!   'lapack_cov'  - Snapshot method via eigendecomposition (optimal for nx >> nt)
!   'lapack_svd'  - Direct LAPACK SVD (dgesvd, general purpose)
! ARRAY CONVENTION: A(nx,nt), ts(neof,nt), U(nx,nt)
! ====
subroutine calc_eofs(A, U, eigval, ev, cev, rank, method, verbose_level,ierr)
  real(dp), dimension(:,:), intent(in)               :: A
  real(dp), dimension(:,:), allocatable, intent(out) :: U
  real(dp), dimension(:), allocatable, intent(out)   :: eigval, ev, cev
  integer, intent(out)                               :: rank
  character(len=*), intent(in)                       :: method
  integer, intent(in), optional                      :: verbose_level
  integer, intent(out), optional                     :: ierr

  integer :: nx, nt, i, j, error_flag, verbosity
  real(dp) :: aa, total_variance, threshold, xsum
  character(len=50) :: method_lower

  error_flag = 0
  verbosity  = 1
  if (present(verbose_level)) verbosity = verbose_level

  ! Convert method to lowercase for case-insensitive comparison
  method_lower = trim(adjustl(lowercase(method)))

  ! Validate method string
  if (method_lower /= 'legacy' .and. &
      method_lower /= 'lapack_cov' .and. &
      method_lower /= 'lapack_svd') then
    error_flag = 10
    if (present(ierr)) then
      ierr = error_flag
      write(*,'(A)') 'ERROR: Invalid method "'//trim(method)//'"'
      write(*,'(A)') 'Valid methods: legacy, lapack_cov, lapack_svd'
      return
    else
      error stop "calc_eofs: Invalid method specified"
    endif
  endif

  ! System size and memory allocation
  nx = size(A,1)
  nt = size(A,2)

  if (allocated(U))      deallocate(U)
  if (allocated(eigval)) deallocate(eigval)
  if (allocated(ev))     deallocate(ev)
  if (allocated(cev))    deallocate(cev)

  allocate(U(nx,nt), eigval(nt), ev(nt), cev(nt), stat=error_flag)
  if (error_flag /= 0) then
    if (present(ierr)) then
      ierr = error_flag
      return
    else
      error stop "calc_eofs: Memory allocation failed"
    endif
  endif

  ! Calculate total variance
  if (verbosity.ge.1) then
    total_variance = 0.0D0
    do j=1,nt
      do i=1,nx
        aa = A(i,j)
        total_variance = total_variance + aa*aa
      enddo
    enddo
    total_variance = total_variance / real(nt-1, dp)
    write(*,'(A,ES14.6)') 'calc_eofs :: Total variance = ', total_variance
  endif

  ! ====
  ! METHOD SELECTION AND COMPUTATION
  ! ====
  
  select case (method_lower)
  
  case ('legacy')
    ! ----
    ! LEGACY METHOD: Numerical Recipes SVD
    ! Use when: General purpose, small to medium matrices
    ! ----
    if (verbosity.ge.1) then
      write(*,'(A)') 'calc_eofs :: Using LEGACY method (Numerical Recipes SVD)'
    endif
    call calc_eofs_legacy(A, U, eigval, nx, nt, error_flag)
    
  case ('lapack_cov')
    ! ----
    ! SNAPSHOT METHOD: LAPACK eigendecomposition of A^T * A
    ! Use when: nx >> nt (many spatial points, few time steps)
    ! Advantage: Computes eigenvalues of nt x nt matrix instead of nx x nt SVD
    ! ----
    if (verbosity.ge.1) then
      write(*,'(A)') 'calc_eofs :: Using LAPACK_SYEV method (covariance matrix)'
    endif
    call calc_eofs_lapack_cov(A, U, eigval, nx, nt, error_flag)
    
  case ('lapack_svd')
    ! ----
    ! LAPACK SVD METHOD: Direct dgesvd call
    ! Use when: Maximum accuracy needed, LAPACK available
    ! ----
    if (verbosity.ge.1) then
      write(*,'(A)') 'calc_eofs :: Using LAPACK_SVD method (dgesvd)'
    endif
    call calc_eofs_lapack_svd(A, U, eigval, nx, nt, error_flag)
    
  end select

  if (error_flag /= 0) then
    if (present(ierr)) then
      ierr = error_flag
      return
    else
      error stop 'calc_eofs: Computation failed'
    endif
  endif

  ! Sort eigenvalues and eigenvectors in descending order
  call eigsort(eigval, U)

  ! Determine rank based on explained variance threshold
  threshold = 1.0e-6_dp
  rank = 0
  do i=1,nt
    if (eigval(i)/eigval(1) .gt. threshold) then
      rank = rank + 1
    else
      eigval(i) = 0.0D0
    endif
  enddo

  write(*,'(A,I0,A,I0)') 'calc_eofs :: Effective rank = ', rank, ' out of ', nt

  ! Normalize eigenvalues
  eigval(:) = eigval(:) / DBLE(nt-1)
  xsum = sum(eigval)

  ! Calculate explained and cumulative-explained variance
  cev(1) = eigval(1)
  do i=2,nt
    cev(i) = cev(i-1) + eigval(i)
  enddo
  ev(:)  = 100.0_dp * eigval(:) / xsum
  cev(:) = 100.0_dp * cev(:) / xsum

  if (present(ierr)) ierr = 0

end subroutine calc_eofs

! ====
! SUBROUTINE: calc_eofs_legacy
! PURPOSE: Legacy method using Numerical Recipes SVD
! ====
subroutine calc_eofs_legacy(A, U, eigval, nx, nt, ierr)
  real(dp), dimension(:,:), intent(in)    :: A
  real(dp), dimension(:,:), intent(out)   :: U
  real(dp), dimension(:), intent(out)     :: eigval
  integer, intent(in)                    :: nx, nt
  integer, intent(out)                    :: ierr

  real(dp), dimension(:,:), allocatable :: V
  real(dp), dimension(:), allocatable :: D

  allocate(V(nt,nt), D(nt), stat=ierr)
  if (ierr /= 0) return

  call svdcmp(A, U, D, V, ierr=ierr)
  if (ierr /= 0) then
    deallocate(V, D)
    return
  endif

  eigval(:) = D(:)**2

  deallocate(V, D)

end subroutine calc_eofs_legacy

! ====
! SUBROUTINE: calc_eofs_snapshot
! PURPOSE: Snapshot method using LAPACK dsyevd for eigendecomposition
! DESCRIPTION: When nx >> nt, it's more efficient to compute:
!              C = A^T * A (size nt x nt)
!              Then find eigenvalues/eigenvectors of C
!              Finally recover spatial EOFs: U = A * V * Lambda^(-1/2)
! ====
subroutine calc_eofs_lapack_cov(A, U, eigval, nx, nt, ierr)
  real(dp), dimension(:,:), intent(in)    :: A           ! (nx,nt)
  real(dp), dimension(:,:), intent(out)   :: U           ! (nx,nt)
  real(dp), dimension(:), intent(out)     :: eigval      ! (nt)
  integer, intent(in)                     :: nx, nt
  integer, intent(out)                    :: ierr

  real(dp), dimension(:,:), allocatable :: C
  real(dp), dimension(:), allocatable   :: work,eigvalnx
  integer, dimension(:), allocatable    :: iwork
  integer :: lwork, liwork, info, i, j, irank

  ierr = 1  ! default to error
  irank = min(nx,nt)

  if (nx < nt) then
    ! Use C = A AT  (nx x nx)
    allocate(C(nx,nx), eigvalnx(nx), stat=ierr)
    if (ierr /= 0) return

    ! C = A * A^T using BLAS dgemm for efficiency
    call dgemm('N', 'T', nx, nx, nt, 1.0_dp, A, nx, A, nx, 0.0_dp, C, nx)   

    ! Workspace query for dsyevd  
    lwork = -1  
    liwork = -1  
    allocate(work(1), iwork(1))  
    call dsyevd('V', 'U', nx, C, nx, eigvalnx, work, lwork, iwork, liwork, info)  
  
    lwork = int(work(1))  
    liwork = iwork(1)  
    deallocate(work, iwork)  
    allocate(work(lwork), iwork(liwork), stat=ierr)  
    if (ierr /= 0) then  
      deallocate(C)  
      return  
    endif  

    ! Compute eigenvalues and eigenvectors  
    call dsyevd('V', 'U', nx, C, nx, eigvalnx, work, lwork, iwork, liwork, info)  
    print*, 'info = ', info
    if (info /= 0) then  
      ierr = info  
      deallocate(C, work, iwork)  
      return  
    endif  
   
    eigval(:) = 0.0_dp
    U(:,:) = 0.0_dp
    do i=nx,nx-irank+1,-1
      eigval(nx-i+1) = eigvalnx(i)
      U(:,nx-i+1)    = C(:,i)
    enddo
    deallocate(eigvalnx)

  else
    ! Compute covariance matrix C = A^T * A (nt x nt)
    allocate(C(nt,nt), stat=ierr)
    if (ierr /= 0) return

    ! C = A^T * A using BLAS dgemm for efficiency
    call dgemm('T', 'N', nt, nt, nx, 1.0_dp, A, nx, A, nx, 0.0_dp, C, nt)

    ! Workspace query for dsyevd
    lwork = -1
    liwork = -1
    allocate(work(1), iwork(1))
    call dsyevd('V', 'U', nt, C, nt, eigval, work, lwork, iwork, liwork, info)
  
    lwork = int(work(1))
    liwork = iwork(1)
    deallocate(work, iwork)
    allocate(work(lwork), iwork(liwork), stat=ierr)
    if (ierr /= 0) then
      deallocate(C)
      return
    endif

    ! Compute eigenvalues and eigenvectors
    call dsyevd('V', 'U', nt, C, nt, eigval, work, lwork, iwork, liwork, info)
  
    if (info /= 0) then
      ierr = info
      deallocate(C, work, iwork)
      return
    endif

    ! Recover spatial EOFs: U = A * V 
    ! U = A * V_temp, then normalize by singular values to get unit-norm EOFs
    call dgemm('N', 'N', nx, nt, nt, 1.0_dp, A, nx, C, nt, 0.0_dp, U, nx)
  
    do i = 1, nt
      if (eigval(i) > epsilon(1.0d0)) then
        U(:,i) = U(:,i) / sqrt(eigval(i))
      endif
    enddo

  endif

  deallocate(C, work, iwork)
  ierr = 0

end subroutine calc_eofs_lapack_cov

! ====
! SUBROUTINE: calc_eofs_lapack_svd
! PURPOSE: Direct SVD using LAPACK dgesvd
! ====
subroutine calc_eofs_lapack_svd(A, U, eigval, nx, nt, ierr)
  real(dp), dimension(:,:), intent(in)    :: A
  real(dp), dimension(:,:), intent(out)   :: U
  real(dp), dimension(:), intent(out)     :: eigval
  integer, intent(in)                    :: nx, nt
  integer, intent(out)                    :: ierr

  real(dp), dimension(:,:), allocatable :: A_copy, VT
  real(dp), dimension(:), allocatable   :: S, work
  integer :: lwork, info

  allocate(A_copy(nx,nt), VT(nt,nt), S(nt), stat=ierr)
  if (ierr /= 0) return

  A_copy = A

  ! Workspace query
  lwork = -1
  allocate(work(1))
  call dgesvd('S', 'N', nx, nt, A_copy, nx, S, U, nx, VT, nt, work, lwork, info)
  
  lwork = int(work(1))
  deallocate(work)
  allocate(work(lwork), stat=ierr)
  if (ierr /= 0) then
    deallocate(A_copy, VT, S)
    return
  endif

  ! Compute SVD
  call dgesvd('S', 'N', nx, nt, A_copy, nx, S, U, nx, VT, nt, work, lwork, info)
  
  if (info /= 0) then
    ierr = info
    deallocate(A_copy, VT, S, work)
    return
  endif

  ! Convert singular values to eigenvalues
  eigval(:) = S(:)**2

  deallocate(A_copy, VT, S, work)
  ierr = 0

end subroutine calc_eofs_lapack_svd

! ====
! SUBROUTINE: eigsort
! PURPOSE: Sort eigenvalues and eigenvectors in descending order
! ====
subroutine eigsort(d,v)
  real(dp), dimension(:), intent(inout) :: d
  real(dp), dimension(:,:), intent(inout) :: v

  integer n,r,i
  integer, allocatable :: perm(:)
  real(dp), allocatable :: dtmp(:), vtmp(:,:)

  n = size(v,2)
  r = size(v,1)
  allocate(perm(n), dtmp(n), vtmp(r,n))

  call argsort(d,perm,descending=.true.)

  do i = 1, n
    dtmp(i)    = d(perm(i))
    vtmp(:, i) = v(:, perm(i))
  end do

  d = dtmp
  v = vtmp
end subroutine eigsort

! ====
! SUBROUTINE: normalize_eofs
! PURPOSE: Normalize EOF vectors to unit length
! ====
subroutine normalize_eofs(eofs, norms)
  real(dp), dimension(:,:), intent(inout) :: eofs
  real(dp), dimension(:), intent(out), optional :: norms
  
  integer :: i, nmodes
  real(dp) :: norm_val
  
  nmodes = size(eofs, 2)
  
  do i = 1, nmodes
    norm_val = sqrt(sum(eofs(:,i)**2))
    if (norm_val > epsilon(1.0_dp)) then
      eofs(:,i) = eofs(:,i) / norm_val
    endif
    if (present(norms)) then
      if (size(norms) >= i) norms(i) = norm_val
    endif
  enddo

end subroutine normalize_eofs

! ====
! FUNCTION: check_convergence
! PURPOSE: Check convergence based on relative change
! ====
function check_convergence(old_val, new_val, tolerance) result(converged)
  real(dp), intent(in) :: old_val, new_val
  real(dp), intent(in), optional :: tolerance
  logical :: converged

  real(dp) :: rel_change, tol

  tol = 1.0e-6_dp
  if (present(tolerance)) tol = tolerance

  if (abs(old_val) < epsilon(1.0_dp)) then
    converged = .false.
  else
    rel_change = abs(new_val - old_val) / abs(old_val)
    converged = (rel_change <= tol)
  endif
  
end function check_convergence

! ====
! SUBROUTINE: rotate_eofs
! PURPOSE: Rotate EOFs using VARIMAX rotation for better interpretability
! ARRAY CONVENTION: EOF(nx,neof), ts(neof,nt) [CHANGED from ts(nt,neof)]
! BUG FIXES:
!   - Fixed PC array dimensions and transpose logic
!   - Fixed iangle initialization
!   - Added convergence guards
! ====
subroutine rotate_eofs(EOF, ts, irang, method, ierr)
  real(dp), dimension(:,:), intent(inout)        :: EOF, ts
  integer, dimension(:), intent(out)             :: irang
  character(len=*), intent(in)                   :: method
  integer, intent(out), optional                 :: ierr

  integer :: N, NEOF, NT, k, n1, n2, istep, iangle, error_flag
  real(dp) :: vari, glovar, varima, xstart, range, oldvar
  real(dp) :: angle, xdum, change
  real(dp), dimension(:), allocatable :: xvari
  real(dp), dimension(:,:), allocatable :: PAIR, PPAIR, PC
  integer :: iteration_count
  integer, parameter :: MAX_OUTER_ITERATIONS = 100
  character(len=50) :: method_lower

  error_flag = 0
  
  ! Convert method to lowercase
  method_lower = trim(adjustl(lowercase(method)))
  
  ! For now, method parameter is reserved for future use (e.g., different rotation criteria)
  ! Currently only VARIMAX is implemented
  if (method_lower /= 'varimax' .and. method_lower /= 'legacy') then
    write(*,'(A)') 'WARNING: Only VARIMAX rotation implemented, using VARIMAX'
  endif

  N = size(EOF,1)      ! Number of spatial points
  NEOF = size(EOF,2)   ! Number of EOFs
  NT = size(ts,2)      ! Number of time steps (NEW CONVENTION: ts is neof x nt)

  ! Check array sizes
  if (size(ts,1) /= NEOF) then
    error_flag = 20
    if (present(ierr)) then
      ierr = error_flag
      write(*,'(A)') 'ERROR: ts array size mismatch in rotate_eofs'
      write(*,'(A,I0,A,I0)') '       Expected ts(neof,nt) = ts(', NEOF, ',', NT, ')'
      write(*,'(A,I0,A,I0)') '       Got ts(', size(ts,1), ',', size(ts,2), ')'
      return
    else
      error stop 'rotate_eofs: ts array dimension mismatch'
    endif
  endif

  if (size(irang) /= NEOF) then
    error_flag = 21
    if (present(ierr)) then
      ierr = error_flag
      return
    else
      error stop 'rotate_eofs: irang array size mismatch'
    endif
  endif

  ! FIX: Allocate PC with correct dimensions for transpose of ts
  ! ts is (neof, nt), so PC = transpose(ts) should be (nt, neof)
  allocate(PC(NT, NEOF), PAIR(N,2), PPAIR(NT,2), stat=error_flag)
  if (error_flag /= 0) then
    if (present(ierr)) then
      ierr = error_flag
      return
    else
      error stop 'rotate_eofs: Memory allocation failed'
    endif
  endif

  ! Transpose time series for rotation
  PC = transpose(ts)  ! PC(nt, neof)

  vari = vari_all(EOF,N,NEOF)
  write(*,'(A,ES14.6)') 'rotate_eofs :: Initial VARIMAX = ', vari

  ! Main rotation loop with convergence guard
  iteration_count = 0
  
10 continue
  iteration_count = iteration_count + 1
  
  ! Convergence guard for outer loop
  if (iteration_count > MAX_OUTER_ITERATIONS) then
    write(*,'(A,I0,A)') 'WARNING: rotate_eofs reached maximum iterations (', &
                    MAX_OUTER_ITERATIONS, '), stopping rotation'
    goto 999
  endif
  
  glovar = vari
  
  ! Pairwise rotation of all EOF pairs
  do n1=1,NEOF-1
    do n2=n1+1,NEOF
      varima = 0.0D0
      xstart = 0.0D0
      range = half_pi
      istep = 90
      iangle = 0  ! FIX: Initialize iangle to prevent undefined behavior

100   oldvar = varima
      do k=1,istep
        angle = DBLE(k)*range/istep + xstart
        call ROTATE(n1,n2,EOF,N,NEOF,PAIR,angle)
        call VARIMAX(PAIR,N,xdum)
        if (xdum.gt.varima) then
          varima = xdum
          iangle = k
        endif
      enddo

      ! Check for convergence in angle search
      if (oldvar.gt.0.0D0) then
        change = 100.0D0*abs(oldvar-varima)/oldvar
      else
        change = 100.0D0
      endif
      
      if (change.gt.0.1D0) then
        xstart = DBLE(iangle)*range/istep + xstart
        range = 4.0D0*range/istep
        xstart = xstart - 0.5D0*range
        goto 100
      endif

      ! Apply the optimal rotation
      angle = DBLE(iangle)*range/istep + xstart
      call ROTATE(n1,n2,EOF,N,NEOF,PAIR,angle)
      EOF(:,n1) = PAIR(:,1)
      EOF(:,n2) = PAIR(:,2)
      call ROTATE(n1,n2,PC,NT,NEOF,PPAIR,angle)
      PC(:,n1) = PPAIR(:,1)
      PC(:,n2) = PPAIR(:,2)
    enddo
  enddo

  ! Check global convergence
  vari = vari_all(EOF,N,NEOF)
  change = 100.0D0*abs(glovar-vari)/max(glovar, epsilon(1.0_dp))
  write(*,'(A,ES14.6,A,ES12.4,A,I0)') &
    'rotate_eofs :: VARIMAX = ', vari, ' change = ', change, '% (iter ', iteration_count, ')'
  
  if (change.gt.0.1D0) goto 10

999 continue

  ! Compute variance of each rotated EOF for sorting
  allocate(xvari(NEOF), stat=error_flag)
  if (error_flag /= 0) then
    if (present(ierr)) then
      ierr = error_flag
      deallocate(PC, PAIR, PPAIR)
      return
    else
      error stop 'rotate_eofs: Memory allocation failed for xvari'
    endif
  endif

  do k=1,NEOF
    xvari(k) = dot_product(EOF(:,k),EOF(:,k))
  enddo

  ! Transpose PC back to ts convention (neof, nt)
  ts = transpose(PC)
  
  call argsort(xvari,irang,descending=.true.)
  
  deallocate(PC, PAIR, PPAIR, xvari)
  
  if (present(ierr)) ierr = 0

end subroutine rotate_eofs

! ====
! SUBROUTINE: ROTATE
! PURPOSE: Rotate a pair of EOFs by a given angle
! ====
subroutine ROTATE(N1,N2,EOF,N,NEOF,PAIR,angle)
  integer, intent(in) :: N1,N2,N,NEOF
  real(dp), dimension(N,NEOF), intent(in) :: EOF
  real(dp), dimension(N,2), intent(out) :: PAIR
  real(dp), intent(in) :: angle
  integer :: i

  do i=1,N
    PAIR(i,1) = cos(angle)*EOF(i,N1) + sin(angle)*EOF(i,N2)
    PAIR(i,2) = -sin(angle)*EOF(i,N1) + cos(angle)*EOF(i,N2)
  enddo
  
end subroutine ROTATE

! ====
! SUBROUTINE: VARIMAX
! PURPOSE: Compute VARIMAX criterion for a pair of EOFs
! ====
subroutine VARIMAX(PAIR,N,VARIM)
  integer, intent(in) :: N
  real(dp), dimension(N,2), intent(in) :: PAIR
  real(dp), intent(out) :: VARIM
  integer :: i
  real(dp) :: sum11,sum12,sum21,sum22

  sum11 = 0.0D0; sum12 = 0.0D0; sum21 = 0.0D0; sum22 = 0.0D0
  do i=1,N
    sum11 = sum11 + PAIR(i,1)**4
    sum12 = sum12 + PAIR(i,2)**4
    sum21 = sum21 + PAIR(i,1)**2
    sum22 = sum22 + PAIR(i,2)**2
  enddo
  VARIM = N*sum11 - sum21*sum21 + N*sum12 - sum22*sum22
  VARIM = VARIM/(N*N)
  
end subroutine VARIMAX

! ====
! FUNCTION: vari_all
! PURPOSE: Compute total VARIMAX criterion for all EOFs
! ====
function vari_all(EOF,N,NEOF) result(vari)
  integer, intent(in) :: N,NEOF
  real(dp), dimension(N,NEOF), intent(in) :: EOF
  real(dp) :: vari
  integer :: i,k
  real(dp) :: xsum1,xsum2,tmp

  vari = 0.0D0
  do k=1,NEOF
    xsum1 = 0.0D0
    xsum2 = 0.0D0
    do i=1,N
      tmp = EOF(i,k)
      xsum1 = xsum1 + tmp**4
      xsum2 = xsum2 + tmp**2
    enddo
    vari = vari + N*xsum1 - xsum2*xsum2
  enddo
  vari = vari/(N*N)
  
end function vari_all

! ====
! SUBROUTINE: project_eofs
! PURPOSE: Project data onto EOFs to get time series (PCs)
! DESCRIPTION: Compute ts = EOF^T * X where X is the data matrix
! ARRAY CONVENTION: EOF(nx,neof), X(nx,nt), ts(neof,nt) [NEW CONVENTION]
! METHOD: 'direct' uses matrix multiply, 'legacy' uses dot products
! ====
subroutine project_eofs(EOF, X, ts, method, ierr)
  real(dp), dimension(:,:), intent(in)    :: EOF  ! (nx, neof)
  real(dp), dimension(:,:), intent(in)    :: X    ! (nx, nt)
  real(dp), dimension(:,:), intent(out)   :: ts   ! (neof, nt)
  character(len=*), intent(in)            :: method
  integer, intent(out), optional          :: ierr

  integer :: nx, neof, nt, i, j, error_flag
  character(len=50) :: method_lower

  error_flag = 0

  nx = size(EOF, 1)
  neof = size(EOF, 2)
  nt = size(X, 2)

  ! Validate dimensions
  if (size(X,1) /= nx) then
    error_flag = 30
    if (present(ierr)) then
      ierr = error_flag
      write(*,'(A)') 'ERROR: X first dimension does not match EOF'
      return
    else
      error stop 'project_eofs: Dimension mismatch'
    endif
  endif

  if (size(ts,1) /= neof .or. size(ts,2) /= nt) then
    error_flag = 31
    if (present(ierr)) then
      ierr = error_flag
      write(*,'(A)') 'ERROR: ts dimensions incorrect'
      write(*,'(A,I0,A,I0,A)') '       Expected ts(', neof, ',', nt, ')'
      return
    else
      error stop 'project_eofs: ts dimension mismatch'
    endif
  endif

  ! Convert method to lowercase
  method_lower = trim(adjustl(lowercase(method)))

  select case (method_lower)
  
  case ('direct', 'lapack')
    ! Use BLAS matrix multiplication: ts = EOF^T * X
    call dgemm('T', 'N', neof, nt, nx, 1.0_dp, EOF, nx, X, nx, 0.0_dp, ts, neof)
    
  case ('legacy', 'default')
    ! Legacy method using explicit loops
    do j = 1, nt
      do i = 1, neof
        ts(i,j) = dot_product(EOF(:,i), X(:,j))
      enddo
    enddo
    
  case default
    write(*,'(A)') 'WARNING: Unknown method "'//trim(method)//'", using direct'
    call dgemm('T', 'N', neof, nt, nx, 1.0_dp, EOF, nx, X, nx, 0.0_dp, ts, neof)
    
  end select

  if (present(ierr)) ierr = 0

end subroutine project_eofs

! ====
! SUBROUTINE: calc_cov_eigs
! PURPOSE: Eigen-decomposition of a covariance matrix C(nc,nc) using LAPACK
! INPUT:
!   C(nc,nc)  - Symmetric covariance matrix (real(dp))
! OUTPUT:
!   eigval(:) - Allocated on exit, eigenvalues (descending order)
!   eigvec(:,:)
!             - Allocated on exit, eigenvectors (columns), matching eigval
! OPTIONAL:
!   verbose_level - 0: silent, >=1: basic info
!   ierr          - 0 on success, >0 on error (LAPACK info or local code)
! NOTES:
!   - Uses LAPACK dsyevd (already interfaced at module top)
!   - Does not normalize eigenvalues further (C is assumed already a covariance)
!   - Preserves input C (works on an internal copy)
! ====
subroutine calc_cov_eigs(C, eigval, eigvec, verbose_level, ierr)
  real(dp), dimension(:,:), intent(in)               :: C
  real(dp), dimension(:),   allocatable, intent(out) :: eigval
  real(dp), dimension(:,:), allocatable, intent(out) :: eigvec
  integer, intent(in),  optional                     :: verbose_level
  integer, intent(out), optional                     :: ierr

  integer :: nc, lwork, liwork, info
  integer :: error_flag, verbosity
  real(dp), dimension(:,:), allocatable :: A
  real(dp), dimension(:),   allocatable :: work
  integer, dimension(:),    allocatable :: iwork

  error_flag = 0
  verbosity  = 1
  if (present(verbose_level)) verbosity = verbose_level

  ! --- Basic size checks -------------------------------------------------
  if (size(C,1) /= size(C,2)) then
    error_flag = 40     ! Non-square input
    if (present(ierr)) then
      ierr = error_flag
      write(*,'(A)') 'ERROR: calc_cov_eigs :: C must be square'
      return
    else
      error stop 'calc_cov_eigs: C must be square'
    endif
  endif

  nc = size(C,1)

  ! Deallocate outputs if already allocated
  if (allocated(eigval)) deallocate(eigval)
  if (allocated(eigvec)) deallocate(eigvec)

  allocate(eigval(nc), eigvec(nc,nc), stat=error_flag)
  if (error_flag /= 0) then
    if (present(ierr)) then
      ierr = error_flag
      write(*,'(A)') 'ERROR: calc_cov_eigs :: Memory allocation failed (eigval/eigvec)'
      return
    else
      error stop 'calc_cov_eigs: Allocation failed for eigval/eigvec'
    endif
  endif

  allocate(A(nc,nc), stat=error_flag)
  if (error_flag /= 0) then
    if (present(ierr)) then
      ierr = error_flag
      write(*,'(A)') 'ERROR: calc_cov_eigs :: Memory allocation failed (A)'
      return
    else
      error stop 'calc_cov_eigs: Allocation failed for A'
    endif
  endif

  ! Copy input covariance (dsyevd overwrites its matrix argument)
  A = C

  ! --- Workspace query for dsyevd ----------------------------------------
  lwork  = -1
  liwork = -1
  allocate(work(1), iwork(1), stat=error_flag)
  if (error_flag /= 0) then
    if (present(ierr)) then
      ierr = error_flag
      write(*,'(A)') 'ERROR: calc_cov_eigs :: Memory allocation failed (query workspace)'
      deallocate(A)
      return
    else
      error stop 'calc_cov_eigs: Allocation failed for query workspace'
    endif
  endif

  call dsyevd('V', 'U', nc, A, nc, eigval, work, lwork, iwork, liwork, info)

  ! info from this query call is not used; we only need work(1), iwork(1)
  lwork  = int(work(1))
  liwork = iwork(1)

  deallocate(work, iwork)

  allocate(work(lwork), iwork(liwork), stat=error_flag)
  if (error_flag /= 0) then
    if (present(ierr)) then
      ierr = error_flag
      write(*,'(A)') 'ERROR: calc_cov_eigs :: Memory allocation failed (workspace)'
      deallocate(A)
      return
    else
      error stop 'calc_cov_eigs: Allocation failed for workspace'
    endif
  endif

  ! --- Actual eigen-decomposition ----------------------------------------
  call dsyevd('V', 'U', nc, A, nc, eigval, work, lwork, iwork, liwork, info)

  if (info /= 0) then
    error_flag = info
    if (present(ierr)) then
      ierr = error_flag
      write(*,'(A,I0)') 'ERROR: calc_cov_eigs :: LAPACK dsyevd failed, info = ', info
      if (allocated(A))     deallocate(A)
      if (allocated(work))  deallocate(work)
      if (allocated(iwork)) deallocate(iwork)
      return
    else
      error stop 'calc_cov_eigs: LAPACK dsyevd failed'
    endif
  endif

  ! A now contains eigenvectors (columns), eigval contains eigenvalues
  eigvec = A

  ! dsyevd returns eigenvalues in ascending order; sort descending
  call eigsort(eigval, eigvec)

  if (verbosity >= 1) then
    write(*,'(A,I0)') 'calc_cov_eigs :: Matrix dimension nc = ', nc
    write(*,'(A,ES14.6)') 'calc_cov_eigs :: Largest eigenvalue = ', eigval(1)
  endif

  if (allocated(A))     deallocate(A)
  if (allocated(work))  deallocate(work)
  if (allocated(iwork)) deallocate(iwork)

  if (present(ierr)) ierr = 0

end subroutine calc_cov_eigs

! ====
! END OF MODULE
! ====
end module module_eofs
