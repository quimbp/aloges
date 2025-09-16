! ======================================================================== !
! EMPIRICAL ORTHOGONAL FUNCTIONS (EOF) MODULE                              !
! Based on ALOGES PROJECT by Quim Ballabrera                               !
! Last Modified: 2024-01-15                                                !
!                                                                          !
! Copyright (C) 2022, Joaquim Ballabrera                                   !
! This program is free software: you can redistribute it and/or modify     !
! it under the terms of the GNU Lesser General Public License as published !
! by the Free Software Foundation, either version 3 of the License, or     !
! (at your option) any later version.                                      !
! ======================================================================== !

module module_eofs

use, intrinsic :: IEEE_ARITHMETIC, ONLY : IEEE_VALUE, IEEE_QUIET_NAN
use module_types
use module_constants
use module_math, only : indexx
use module_linalg

implicit none

private
public :: svdcmp, pythag, svbksb, get_eof, eigsort, rotate_eofs, &
          ROTATE, VARIMAX, vari_all, &
          normalize_eofs, check_convergence 
    

contains

! ===================================================================
! ===================================================================
! ...
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

    error_flag = 0

    ! ... Bounds checking
    ! ...
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

    do k=n,1,-1
      do its=1,30
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
              h = pythag(f,g)
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
          
        if (its.eq.30) stop 'no convergence in svdcmp'
        x = W(l)
        nm = k - 1
        y = W(nm)
        g = rv1(nm)
        h = rv1(k)
        f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0d0*h*y)
        g = pythag(f,1.0d0)
        f = ((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
        c = 1.0d0
        s = 1.0d0
        do j=l,nm
          i = j + 1
          g = rv1(i)
          y = W(i)
          h = s*g
          g = c*g
          z = pythag(f,h)
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
          z = pythag(f,h)
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
  ! ...
  ! ===================================================================
  ! ...
  pure function pythag(a,b) result(c)

    real(dp), intent(in)                   :: a
    real(dp), intent(in)                   :: b
    real(dp)                               :: c

    ! ... Local variables
    ! ...
    real(dp) aa,ab

    aa = abs(a)
    ab = abs(b)

    if (ab.eq.0) then
      c = 0.0D0
    else if (aa.gt.ab) then
      c = aa*dsqrt(1.0D0 + (ab/aa)**2)
    else
      c = ab*dsqrt(1.0D0 + (aa/ab)**2)
    endif

  end function pythag
  ! ...
  ! ===================================================================
  ! ...
  pure function svbksb(U,W,V,b) result(x)
    real(dp), dimension(:,:), intent(in) :: U,V
    real(dp), dimension(:), intent(in) :: W,b
    real(dp), dimension(size(U,1)) :: x

    ! ... Local variables
    ! ...
    integer :: m,n,i
    real(dp) threshold
    real(dp), dimension(size(U,1)) :: tmp

    m = size(U,1)
    n = size(U,2)

    ! ... Set a threshold for singular values (machine epsilon based)
    ! ...
    threshold = max(m, n) * epsilon(1.0_dp) * maxval(W)
 
    tmp = 0.0_dp
    do i = 1, n
      if (W(i) > threshold) then
        tmp(i) = dot_product(b, U(:,i)) / W(i)
      endif   
    enddo   

    x = matmul(V, tmp)
 
  end function svbksb
! ...
! ===================================================================
! ...
  subroutine get_eof(A,U,eigval,ev,cev,rank,ierr)
    real(dp), dimension(:,:), intent(in)               :: A
    real(dp), dimension(:,:), allocatable, intent(out) :: U
    real(dp), dimension(:), allocatable, intent(out)   :: eigval, ev, cev
    integer, intent(out)                               :: rank
    integer, intent(out), optional                     :: ierr

    integer :: nx,nt,i,j,error_flag
    real(dp) :: aa,total_variance,threshold,xsum
    real(dp), dimension(:,:), allocatable :: V
    real(dp), dimension(:), allocatable :: D

    error_flag = 0

    ! ... System size and memory allocation
    ! ...
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
        error stop "GET_EOF: Memory allocation failed for U, eigval and ev"
      endif
    endif

    ! ... Calculate total variance
    ! ...
    total_variance = 0.0D0
    do j=1,nt
    do i=1,nx
       aa = A(i,j)
       total_variance = total_variance + aa*aa
    enddo
    enddo
    total_variance = total_variance / real(nt-1, dp)
    write(*,*) 'GET_EOF :: Total variance = ', total_variance

    allocate(V(nt,nt), D(nt), stat=error_flag)
    if (error_flag /= 0) then
      if (present(ierr)) then
        ierr = error_flag
        return
      else
        error stop "GET_EOF: Memory allocation failed for V and D"
      endif
    endif

    call svdcmp(A,U,D,V,ierr=error_flag)
    if (error_flag.ne.0) then
      if (present(ierr)) then
        ierr = error_flag
        return
      else
        stop 'GET_EOF: SVD computation failed'
      endif
    endif

    eigval(:) = D(:)**2
    call eigsort(eigval,U)

    ! ... Determine rank based on explained variance threshold
    ! ...
    threshold = 1.0e-6_dp     ! Adjustable threshold
    rank = 0
    do i=1,nt
      if (eigval(i)/eigval(1).gt.threshold) then
        rank = rank + 1
      else
        D(i) = 0.0D0
        eigval(i) = 0.0D0
      endif
    enddo

    ! ... Normalize eigenvalues
    ! ...
    eigval(:) = eigval(:) / DBLE(nt-1)
    xsum = sum(eigval)

    ! ... Calculate explained and cumulative-explained variance
    ! ...
    cev(1) = eigval(1)
    DO i=2,nt
      cev(i) = cev(i-1) + eigval(i)
    ENDDO
    ev(:)  = 100.0_dp * eigval(:) / xsum
    cev(:) = 100.0_dp * cev(:) / xsum

    deallocate(V, D)

    if (present(ierr)) ierr = 0

  end subroutine get_eof
! ...
! ===================================================================
! ...
  subroutine eigsort(d,v)
    real(dp), dimension(:), intent(inout) :: d
    real(dp), dimension(:,:), intent(inout) :: v
    integer :: n,r,i,j,k
    real(dp) :: p

    n = size(v,1)
    r = size(v,2)

    do i=1,r-1
      k=i
      p=d(i)
      do j=i+1,r
        if (d(j).ge.p) then
          k=j
          p=d(j)
        endif
      ENDDO
      IF (k.NE.i) THEN
        d(k)=d(i)
        d(i)=p
        DO j=1,n
          p=v(j,i)
          v(j,i)=v(j,k)
          v(j,k)=p
        enddo
      endif
    enddo
  end subroutine eigsort
  ! ...
  ! ===================================================================
  ! ...     
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
  ! ...
  ! ===================================================================
  ! ...
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
  ! ...
  ! ===================================================================
  ! ...
  subroutine rotate_eofs(EOF,ts,irang)
    real(dp), dimension(:,:), intent(inout) :: EOF, ts
    integer, dimension(size(EOF,2)) :: irang

    integer :: N,NEOF,NT,k,n1,n2,istep,iangle
    real(dp) :: vari,glovar,varima,xstart,range,oldvar
    real(dp) :: angle,xdum,change
    real(dp), dimension(:), allocatable :: xvari
    real(dp), dimension(:,:), allocatable :: PAIR,PPAIR
    real(dp), dimension(size(ts,2),size(ts,1)) :: PC

    PC = TRANSPOSE(ts)
    N = size(EOF,1)
    NEOF = size(EOF,2)
    NT = size(PC,2)

    vari = vari_all(EOF,N,NEOF)
    WRITE(*,*) 'VARIMAX = ', vari

    ALLOCATE(PAIR(N,2), PPAIR(NT,2))

10 CONTINUE
    glovar = vari
    do n1=1,NEOF-1
    do n2=n1+1,NEOF
      varima = 0.0D0
      xstart = 0.0D0
      range = 0.5D0*3.14D0
      istep = 90

100   oldvar = varima
      do k=1,istep
        angle = DBLE(k)*range/istep + xstart
        CALL ROTATE(n1,n2,EOF,N,NEOF,PAIR,angle)
        CALL VARIMAX(PAIR,N,xdum)
        if (xdum.GT.varima) then
          varima = xdum
          iangle = k
        endif
      enddo

      if (oldvar.GT.0.0) then
        change = 100.0D0*DABS(oldvar-varima)/oldvar
      ELSE
        change = 100.0D0
      endif
      if (change.GT.0.1) then
        xstart = DBLE(iangle)*range/istep + xstart
        range = 4.0D0*range/istep
        xstart = xstart - 0.5D0*range
        GOTO 100
      endif

      angle = DBLE(iangle)*range/istep + xstart
      CALL ROTATE(n1,n2,EOF,N,NEOF,PAIR,angle)
      EOF(:,n1) = PAIR(:,1)
      EOF(:,n2) = PAIR(:,2)
      CALL ROTATE(n1,n2,PC,NT,NEOF,PPAIR,angle)
      PC(:,n1) = PPAIR(:,1)
      PC(:,n2) = PPAIR(:,2)
    enddo
    enddo

    vari = vari_all(EOF,N,NEOF)
    change = 100.0D0*DABS(glovar-vari)/glovar
    WRITE(*,*) 'VARIMAX = ', vari, ' change = ', change
    if (change.GT.0.1) GOTO 10

    DEALLOCATE(PAIR,PPAIR)
    ALLOCATE(xvari(NEOF))

    do k=1,NEOF
      xvari(k) = DOT_PRODUCT(EOF(:,k),EOF(:,k))
    enddo

    ts = TRANSPOSE(PC)
    CALL indexx(xvari,irang)
    DEALLOCATE(xvari)
  end subroutine rotate_eofs
! ...
! ===================================================================
! ...
  subroutine ROTATE(N1,N2,EOF,N,NEOF,PAIR,angle)
    integer, intent(in) :: N1,N2,N,NEOF
    real(dp), dimension(N,NEOF), intent(in) :: EOF
    real(dp), dimension(N,2), intent(out) :: PAIR
    real(dp), intent(in) :: angle
    integer :: i

    do i=1,N
      PAIR(i,1)= DCOS(angle)*EOF(i,N1) + DSIN(angle)*EOF(i,N2)
      PAIR(i,2)= -DSIN(angle)*EOF(i,N1) + DCOS(angle)*EOF(i,N2)
    enddo
  end subroutine ROTATE
! ...
! ===================================================================
! ...
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
! ...
! ===================================================================
! ...
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
! ...
! ===================================================================
! ...
end module module_eofs
