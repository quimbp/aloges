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
! - randn                                                                  !
! - randname                                                               !
! - randseries                                                             !
! -------------------------------------------------------------------------!

module module_math

use, intrinsic :: IEEE_ARITHMETIC, ONLY : IEEE_VALUE, IEEE_QUIET_NAN
use module_types
use module_constants

implicit none

interface randn
  module procedure randn_r,randn_v,randn_a
end interface randn

interface arange
  module procedure arange_i
  module procedure arange_dp
end interface arange

contains
! ...
! =====================================================================
! =====================================================================
! ...
  function randn_r () result(ff)

    real(dp)                                     :: ff

    ! ... Local variables
    ! ...
    real(dp), parameter                          :: s  =  0.449871D0
    real(dp), parameter                          :: t  = -0.386595D0
    real(dp), parameter                          :: a  =  0.19600D0
    real(dp), parameter                          :: b  =  0.25472D0
    real(dp), parameter                          :: r1 =  0.27597D0
    real(dp), parameter                          :: r2 =  0.27846D0
    real(dp) u,v,x,y,q

    do
      call RANDOM_NUMBER(u)  ! GNU RANDOM GENERATOR
      call RANDOM_NUMBER(v)  ! GNU RANDOM GENERATOR
      v = 1.7156D0 * (v - 0.5D0)

      ! ... Evaluate the quadratic form
      ! ...
      x = u - s
      y = ABS(v) - t
      q = x*x + y*(a*y - b*x)

      if (q .lt. r1) exit
      if (q .gt. r2) cycle
      if (v**2 .LT. -4D0*LOG(u)*u*u) exit
    enddo
    ff = v/u

  end function randn_r
  ! ...
  ! ===================================================================
  ! ...
  function randn_v (m) result(ff)

    integer, intent(in)                          :: m
    real(dp), dimension(m)                       :: ff

    ! ... Local variables
    ! ...
    integer i
    real(dp), parameter                          :: s  =  0.449871D0
    real(dp), parameter                          :: t  = -0.386595D0
    real(dp), parameter                          :: a  =  0.19600D0
    real(dp), parameter                          :: b  =  0.25472D0
    real(dp), parameter                          :: r1 =  0.27597D0
    real(dp), parameter                          :: r2 =  0.27846D0
    real(dp) u,v,x,y,q
  
    do i=1,m
      do
        call RANDOM_NUMBER(u)  ! GNU RANDOM GENERATOR
        call RANDOM_NUMBER(v)  ! GNU RANDOM GENERATOR
        v = 1.7156D0 * (v - 0.5D0)
  
        ! ... Evaluate the quadratic form
        ! ...
        x = u - s
        y = ABS(v) - t
        q = x*x + y*(a*y - b*x)

        if (q .lt. r1) exit
        if (q .gt. r2) cycle
        if (v**2 .LT. -4D0*LOG(u)*u*u) exit
      enddo
      ff(i) = v/u
    enddo

  end function randn_v
  ! ...
  ! ===================================================================
  ! ...
  function randn_a (m,n) result(ff)

    integer, intent(in)                          :: m
    integer, intent(in)                          :: n
    real(dp), dimension(:,:), pointer            :: ff

    ! ... Local variables
    ! ...
    integer i,j
    real(dp), parameter                          :: s  =  0.449871D0
    real(dp), parameter                          :: t  = -0.386595D0
    real(dp), parameter                          :: a  =  0.19600D0
    real(dp), parameter                          :: b  =  0.25472D0
    real(dp), parameter                          :: r1 =  0.27597D0
    real(dp), parameter                          :: r2 =  0.27846D0
    real(dp) u,v,x,y,q

    if (.not.associated(ff)) allocate(ff(m,n))

    do j=1,n
    do i=1,m
      do
        call RANDOM_NUMBER(u)  ! GNU RANDOM GENERATOR
        call RANDOM_NUMBER(v)  ! GNU RANDOM GENERATOR
        v = 1.7156D0 * (v - 0.5D0)

        ! ... Evaluate the quadratic form
        ! ...
        x = u - s
        y = ABS(v) - t
        q = x*x + y*(a*y - b*x)

        if (q .lt. r1) exit
        if (q .gt. r2) cycle
        if (v**2 .LT. -4D0*LOG(u)*u*u) exit
      enddo
      ff(i,j) = v/u
    enddo
    enddo

  end function randn_a
  ! ...
  ! ===================================================================
  ! ...
  function randname(len,iseed) result(name)

  integer, intent(in)                            :: len
  integer, optional                              :: iseed
  character(len=len)                             :: name
  
  ! ... Local variables
  ! ...
  integer i,io,il,j,n
  integer, dimension(:), allocatable             :: seed
  real(dp) r

  if (present(iseed)) then
    call random_seed(size=n)
    allocate(seed(n))
    seed(:) = iseed
    call random_seed(put=seed)
  endif

  io = ichar('A')
  il = ichar('Z') - io

  do i=1,len
    call random_number(r)
    j = int(io + il*r)
    name(i:i) = char(j)
  enddo

  return
  end function randname
  ! ...
  ! ===================================================================
  ! ...
  function randseries(n,dlag,periodic) result(g)

  integer, intent(in)                            :: n
  integer, intent(in), optional                  :: dlag
  logical, intent(in), optional                  :: periodic
  real(dp), dimension(n)                         :: g
  
  ! ... Local variables
  ! ...
  logical lper
  integer nn,i,iter,im,ip
  real(dp) xsum
  real(dp), dimension(:), allocatable            :: v,f
  
  lper = .False.
  if (present(periodic)) then
    if (periodic) lper = .True.
  endif
  
  if (lper) then
  
    if (present(dlag)) then
      allocate(f(n))
  
      g = randn(n)
      g(n) = g(1)

      print*, dlag, g(n), g(1)
      do iter=1,6*abs(dlag)
        do i=1,n
          im = i - 1
          if (im.eq.0) im = n
          ip = i + 1
          if (ip.gt.n) ip = 1
          f(i) = 0.25D0*(g(im)+2*g(i)+g(ip))
        enddo
        g(:) = f(:)
      enddo
      deallocate(f)
      xsum = SUM(g)/n
      g(:) = g(:) - xsum
      xsum = sqrt(DOT_PRODUCT(g,g)/n)
      g(:) = g(:) / xsum
    else
      g = randn(n)
      g(n) = g(1)
    endif
  else
    if (present(dlag)) then
  
      ! ... Buffering time series
      ! ...
      nn = n + 20*abs(dlag)
      allocate(v(nn))
      allocate(f(nn))

      v = randn(nn)
      do iter=1,6*abs(dlag)
        f(1) = v(1)
        f(nn) = v(nn)
        do i=2,nn-1
          f(i) = 0.25D0*(v(i-1)+2*v(i)+v(i+1))
        enddo
        v = f
      enddo

      ! ... Extracting far from boundary
      ! ...
      i = 10*abs(dlag)
      g(:) = f(i+1:i+n)

      xsum = SUM(g)/n
      g(:) = g(:) - xsum
      xsum = sqrt(DOT_PRODUCT(g,g)/n)
      g(:) = g(:) / xsum

      deallocate(f)
      deallocate(v)
  
    else
  
      g = randn(n)
  
    endif
  endif

  end function randseries
  ! ...
  ! ===================================================================
  ! ...
  real(dp) function haversine (lon1,lat1,lon2,lat2)
    ! ...
    ! ... Function Haversine
    ! ... Determines the great-circle distance between two points in a
    ! ... sphere. The input lngitudes and latitudes are given in radians
    ! ... The retruned distance is in meters.
    ! ... REarth = 6371315.0_dp      ! m
    ! ...
    real(dp), intent(in)                  :: lon1,lat1
    real(dp), intent(in)                  :: lon2,lat2

    ! ... Local variables
    ! ...
    real(dp) sindx,sindy,dang

    sindx = sin(0.5D0*(lon2-lon1))
    sindy = sin(0.5D0*(lat2-lat1))

    dang = 2.0d0*asin(sqrt(sindy*sindy + cos(lat2)*cos(lat1)*sindx*sindx))
    haversine = REarth * dang

    return
  end function haversine
  ! ...
  ! ===================================================================
  ! ...
  pure function arange_i(amin, amax, adelta) result(arange)
    ! ...
    ! ... Returns a double precision vector from amin to amax with 
    ! ... an optional step of adelta.
    ! ...

    integer, intent(in)                  :: amin 
    integer, intent(in)                  :: amax 
    integer, intent(in), optional        :: adelta 
    integer, dimension(:), allocatable   :: arange

    ! ... Loal variables
    ! ...
    integer  i, nlen
    integer delta

    if (present(adelta))then
      delta = adelta
    else
      delta = 1
    endif

    nlen = (amax-amin+0.5*delta)/delta + 1
    allocate(arange(nlen))

    do i=1,nlen
      arange(i) = amin + (i-1)*delta
    enddo

    return
  end function arange_i
  ! ...
  ! ===================================================================
  ! ...
  pure function arange_dp(amin, amax, adelta) result(arange)
    ! ...
    ! ... Returns a double precision vector from amin to amax with 
    ! ... an optional step of adelta.
    ! ...

    real(dp), intent(in)                  :: amin 
    real(dp), intent(in)                  :: amax 
    real(dp), intent(in), optional        :: adelta 
    real(dp), dimension(:), allocatable   :: arange

    ! ... Loal variables
    ! ...
    integer  i, nlen
    real(dp) delta

    if (present(adelta))then
      delta = adelta
    else
      delta = 1
    endif

    nlen = (amax-amin+0.5*delta)/delta + 1
    allocate(arange(nlen))

    do i=1,nlen
      arange(i) = amin + (i-1.0D0)*delta
    enddo

    return
  end function arange_dp
  ! ...
  ! ===================================================================
  ! ...
  real(dp) function mean(A,W)
    ! ... Calculates the weighted mean = 1/N * Sum W(i)*A(i)
    ! ... Weights are optional.

    real(dp), dimension(:), intent(in)     :: A
    real(dp), dimension(:), optional       :: W

    integer n
    real(dp) nan,Sw

    ! ... nan:
    ! ...
    nan = ieee_value(1.0_dp,ieee_quiet_nan)

    mean = nan
    n = size(A)
    if (n.eq.0) return

    if (present(W)) then
      Sw   = sum(W)
      mean = dot_product(W,A)/Sw
    else
      mean = sum(A)/N
    endif

  end function mean
  ! ...
  ! ===================================================================
  ! ...
  real(dp) function variance (A,W)

    ! ... Calculates the variance = FACTOR * Sum (A(i)-mean(A))**2
    ! ... The value of FACTOR depends of W. By default, W=0, FACTOR = 1/(N-1)
    ! ... However, if W=1, the biased variance is calculated, i.e. FACTOR=1/N
    ! ... If W has to be used as a weight it must be a vector with the same
    ! ... size than A.
    ! ...
    ! ... Weights are optional.
    ! ...

    real(dp), dimension(:), intent(in)     :: A
    real(dp), dimension(:), optional       :: W

    logicaL                                    :: weight = .false.
    logical                                    :: biased = .false.
    integer N,i
    real(dp) xsum1,xsum2,ai,wi,Sw
    real(dp) nan


    ! ... nan:
    ! ...
    nan = ieee_value(1.0_dp,ieee_quiet_nan)

    variance = nan
    N = SIZE(A)
    if (N.EQ.0) return

    if (PRESENT(W)) then
      if (SIZE(W).EQ.1) then
        if (W(1).EQ.0) then
          biased = .false.
        else if (W(1).EQ.1) then
          biased = .true.
        else
          STOP 'ERROR in variance: Invalid normalization. Use 0 (unbiased) or 1 (biased)'
        endif
      else if (SIZE(W).EQ.N) then
        weight = .true.
      else
        STOP 'ERROR in variance: Size W must be 1 or N'
      endif
    endif

    if (weight) then
      Sw    = 0D0
      xsum1 = 0D0
      xsum2 = 0D0
      do i=1,N
        wi = w(i)
        ai = A(i)
        Sw    = Sw    + wi
        xsum1 = xsum1 + wi*ai
        xsum2 = xsum2 + wi*ai*ai
      enddo
      xsum1 = xsum1 / Sw
      xsum2 = xsum2/Sw-xsum1*xsum1
      variance = Sw * xsum2 /(Sw - 1D0)
    else
      xsum1 = 0D0
      xsum2 = 0D0
      do i=1,N
        ai = A(i)
        xsum1 = xsum1 + ai
        xsum2 = xsum2 + ai*ai
      enddo
      xsum1 = xsum1 / N
      xsum2 = xsum2 / N - xsum1*xsum1
      IF (biased) then
        variance = xsum2
      else
        variance = N*xsum2/(N-1D0)
      endif
    endif
    IF (variance.LT.0) variance = 0D0

  end function variance
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
  pure function out_product(A,B) result(C)

    real(dp), dimension(:), intent(in)     :: A
    real(dp), dimension(:), intent(in)     :: B
    real(dp), dimension(size(A),size(B))   :: C
    
    ! ... Local variables
    ! ...
    integer n,m
    n = size(A)
    m = size(B)
    C = spread(A,dim=2,ncopies=m) * spread(B,dim=1,ncopies=n)
 
  end function out_product 
  ! ...
  ! ===================================================================
  ! ...
  subroutine svdcmp(A,U,W,V)

    real(dp), dimension(:,:), intent(in)               :: A
    real(dp), dimension(:,:), intent(out), allocatable :: U,V
    real(dp), dimension(:),   intent(out), allocatable :: W

    ! ... local variables
    ! ...
    integer m,n,i,j,k,l,its,nm
    real(dp) anorm,c,f,g,h,s,scale,x,y,z
    real(dp), dimension(size(A,1))              :: tempm
    real(dp), dimension(size(A,2))              :: rv1,tempn

    m = size(A,1)
    n = size(A,2)

    allocate(U(m,n))
    allocate(V(n,n))
    allocate(W(n))

    U = A

    g     = 0.0d0
    scale = 0.0d0

    do i=1,n                         ! ... Loop over columns U(:,i)
      l = i + 1
      rv1(i) = scale*g
      g      = 0.0d0
      scale  = 0.0d0

      if (i.le.m) then
        scale = sum(ABS(U(i:m,i)))
        if (scale.ne.0.0d0) then
          U(i:m,i) = U(i:m,i)/scale
          s =  DOT_PRODUCT(U(i:m,i),U(i:m,i))
          f =  U(i,i)
          g = -sign(dsqrt(s),f)
          h = f*g - s
          U(i,i) = f - g
          tempn(l:n) = MATMUL(U(i:m,i),U(i:m,l:n))/h
          U(i:m,l:n) = U(i:m,l:n) + OUT_PRODUCT(U(i:m,i),tempn(l:n))
          U(i:m,i)   = scale*U(i:m,i)
        endif
      endif

      W(i)  = scale*g
      g     = 0.0d0
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
      else
        U(i:m,i) = 0.0D0
      endif
      U(i,i) = U(i,i) + 1.0d0
    enddo

    do k=n,1,-1
      do its=1,30

        do l=k,1,-1
          nm = l - 1
          if ((ABS(rv1(l))+anorm).eq.anorm) exit
          if ((ABS(W(nm))+anorm).eq.anorm)  then
            c = 0.0d0
            s = 1.0d0
            do i=l,k
              f = s*rv1(i)
              rv1(i) = c*rv1(i)
              if ((ABS(f)+anorm).eq.anorm) exit
              g    =  W(i)
              h    =  pythag(f,g)
              W(i) =  h
              h    =  1.0D0/h
              c    =  g*h
              s    = -f*h
              tempm(1:m) = U(1:m,nm)
              U(1:m,nm) =  c*U(1:m,nm)  + s*U(1:m,i)
              U(1:m,i)  = -s*tempm(1:m) + c*U(1:m,i)
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
          
        if (its.eq.30) stop  'no convergence in svdcmp'
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
          f =  (x*c)+(g*s)
          g = -(x*s)+(g*c)
          h = y*s
          y = y*c
          tempn(1:n) = V(1:n,j)
          V(1:n,j) = c*V(1:n,j)  + s*V(1:n,i)
          V(1:n,i) =-s*tempn(1:n) + c*V(1:n,i)
          z = pythag(f,h)
          W(j) = z
          if (z.ne.0.0d0) then
            z = 1.0d0/z
            c = f*z
            s = h*z
          endif
          f = (c*g)+(s*y)
          x =-(s*g)+(c*y)
          tempm(1:m) = U(1:m,j)
          U(1:m,j) = c*U(1:m,j)  + s*U(1:m,i)
          U(1:m,i) =-s*tempm(1:m) + c*U(1:m,i)
        enddo
        rv1(l)=0.0d0
        rv1(k)=f
        W(k)=x
      enddo
    enddo

  end subroutine svdcmp
  ! ...
  ! ===================================================================
  ! ...
  pure function svbksb(U,W,V,b) result(x)

    ! ... Solves A x = b, wjere A = V W VT as returned by svdcmp.
    ! ... On input, U(m,n), V(n,n), W(n) and b(m)
    ! ... On output, x(m)
    
    real(dp), dimension(:,:), intent(in)                :: U,V
    real(dp), dimension(:), intent(in)                  :: W,b
    real(dp), dimension(size(U,1))                      :: x

    ! ... Local variables
    ! ...
    integer m,n,i,j,jj
    real(dp), dimension(size(U,1))                      :: tmp

    m = size(U,1)
    n = size(U,2)
 
    where (W.ne.0.0D0)
      tmp = matmul(b,U)/W
    elsewhere
      tmp = 0.0D0
    end where

    x = MATMUL(V,tmp)

  end function svbksb
  ! ...
  ! ===================================================================
  ! ...
  function fitting1d (order,x,y,c,trnd) result(err)

    integer, intent(in)                           :: order
    real(dp), dimension(:), intent(in)            :: x,y
    real(dp), dimension(:), intent(out)           :: trnd
    real(dp), dimension(order+1), intent(out)     :: c
    integer                                       :: err

    ! ... Local variables
    ! ...
    integer                                       :: n,i,k,rank,M
    real(dp)                                      :: ym
    real(dp)                                      :: xn
    real(dp)                                      :: ax,bx,ay,by
    real(dp)                                      :: cond,xi,ypol
    real(dp), dimension(size(x))                  :: xx,yy
    real(dp), dimension(:), allocatable           :: RHS,D
    real(dp), dimension(:,:), allocatable         :: A,U,V,S,ST

    err = 0
    n   = size(x)
    M   = order + 1

    ym = mean(y)
    xn = maxval(abs(x))               ! Largets value of x

    if (order.eq.0) then
      trnd(:) = ym
      return
    endif

    ! ... Normalize the data:
    ! ...
    xx = x / xn                       ! Largest value is 1.
    yy = y - ym

    allocate (S(n,M))
    S(:,1) = 1D0
    do k=1,order
    do i=1,n
      S(i,1+k) = xx(i)**k
    enddo
    enddo

    allocate (A(M,M))
    allocate (U(M,M))
    allocate (V(M,M))
    allocate (RHS(M))
    allocate (D(M))
    allocate (ST(M,n))

    ST  = transpose(S)
    A   = matmul(ST,S)
    RHS = matmul(ST,yy)

    ! ... Solve the system
    ! ...
    CALL svdcmp (A,U,D,V)
    WHERE(D.LE.0D0) D = 0D0

    if (count(D.EQ.0).gt.0) then
      write(*,*) 'Singular matrix'
      err = 1
      trnd = ym
      deallocate (A,U,V,RHS,D,S,ST)
      return
    endif

    cond = maxval(D) / minval(D)

    IF (cond.GT.10000) THEN
      write(*,*) 'Condition number > 10000'
      err = 1
      trnd = ym
      deallocate (A,U,V,RHS,D,S,ST)
      return
    ENDIF

    c = svbksb (U,D,V,RHS)       ! Solve A x = b, A = U D VT

    do i=1,n
      xi   = xx(i)
      ypol = c(1)
      do k=1,order
        ypol = ypol + c(1+k)*xi**k
      enddo
      trnd(i) = ym + ypol
    enddo

    ! ... Remember M = order + 1
    ! ... c(1) = c(1) / xn**0 + ym
    ! ...
    c(1) = c(1) + ym
    do k=2,M
      c(k) = c(k)/xn**(k-1)
    enddo

    deallocate (A,U,V,RHS,D,S,ST)

    end function fitting1d
  ! ...
  ! ===================================================================
  ! ...
  subroutine get_eof (A,U,eigval,ev,rank)

    ! ... Subroutine to calculate the EOFs
    ! ... U, eigval and ev must enter as allocatable !!!!!
    ! ...

    real(dp), dimension(:,:), intent(in)                  :: A
    real(dp), dimension(:,:), allocatable, intent(out)    :: U
    real(dp), dimension(:), allocatable, intent(out)      :: eigval
    real(dp), dimension(:), allocatable, intent(out)      :: ev
    integer, intent(out)                                  :: rank

    integer                                     :: nx,nt,i,j
    real(dp)                                    :: C,aa,xsum,Am
    real(dp), dimension(:,:), allocatable       :: V
    real(dp), dimension(:), allocatable         :: D

    nx = size(A,1)
    nt = size(A,2)

    allocate(eigval(nt))
    allocate(ev(nt))

    ! ... Variance: 
    ! ...
    xsum = 0.0D0
    do j=1,nt
    do i=1,nx
       aa = A(i,j)
       xsum = xsum + aa*aa
    enddo
    enddo

    write(*,*) 'GET_EOF :: Total variance     = ', xsum / DBLE(nt-1)

    call svdcmp (A,U,D,V)
    eigval(:) = D(:)**2
    call eigsort(eigval,U)

    rank = nt
    do i=1,nt
      if (eigval(i)/eigval(1).lt.1.0D-6) then
        rank = rank - 1
        D(i) = 0.0D0
        eigval(i) = 0.0D0
      endif
    enddo
    write(*,*) 'rank = ', rank

    eigval(:) = eigval(:) / DBLE(nt-1)

    xsum = 0.0D0
    do i=1,nx
    do j=1,rank
      aa = U(i,j)
      xsum = xsum + eigval(j)*aa*aa
    enddo
    enddo
    write(*,*) 'GET_EOF :: EOF Total variance = ', xsum 

    ! ... Explained variance
    ! ... 
    ev(:) = 0.0D0

    xsum = SUM(eigval)
    ev(1) = eigval(1)
    DO i=2,rank
      ev(i) = ev(i-1) + eigval(i)
    ENDDO
    ev(:) = 100.0D0*ev(:)/xsum

  end subroutine get_eof 
  ! ...
  ! ===================================================================
  ! ...
  subroutine eigsort(d,v)

    real(dp), dimension(:), intent(inout)   :: d
    real(dp), dimension(:,:), intent(inout) :: v

    integer n,r,i,j,k
    real(dp) p

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
  subroutine rotate_eofs (EOF,ts,irang)

    real(dp), dimension(:,:), intent(inout)   :: EOF
    real(dp), dimension(:,:), intent(inout)   :: ts
    integer, dimension(size(EOF,2))           :: irang

    integer N,NEOF,NT,i,j,k,n1,n2,istep,iangle
    real(dp) vari,xsum1,glovar,varima,xstart,range,oldvar
    real(dp) angle,xdum,change,TMP

    real(dp), DIMENSION(:),   ALLOCATABLE            :: xvari
    real(dp), DIMENSION(:,:), ALLOCATABLE            :: PAIR,PPAIR
    real(dp), DIMENSION(size(ts,2),size(ts,1))       :: PC

    PC = TRANSPOSE(ts)

    N    = size(EOF,1)
    NEOF = size(EOF,2)
    NT   = size(PC,2)

    ! ... Calc VARIMAX over all
    ! ...
    vari = vari_all (EOF,N,NEOF)
    WRITE(*,*) 'VARIMAX = ', vari

    ALLOCATE (PAIR(N,2))
    ALLOCATE (PPAIR(NT,2))

10 CONTINUE

       glovar = vari
       do n1=1,NEOF-1
       do n2=n1+1,NEOF

    ! ... Pairwise rotation
    ! ...
          varima = 0.0D0
          xstart = 0.0D0
          range  = 0.5D0*3.14D0
          istep  = 90

 100      oldvar = varima
          do k=1,istep
            angle = DBLE(k)*range/istep + xstart
            CALL rotate (n1,n2,EOF,N,NEOF,PAIR,angle)
            CALL varimax (PAIR,N,xdum)
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
            range  = 4.0D0*range/istep
            xstart = xstart - 0.5D0*range
            GOTO 100
          endif

          angle = DBLE(iangle)*range/istep + xstart
          xstart = 360.0D0/6.28D0*angle
          CALL rotate (n1,n2,EOF,N,NEOF,PAIR,angle)
          EOF(:,n1) = PAIR(:,1)
          EOF(:,n2) = PAIR(:,2)
          CALL rotate (n1,n2,PC,NT,NEOF,PPAIR,angle)
          PC(:,n1) = PPAIR(:,1)
          PC(:,n2) = PPAIR(:,2)

       enddo
       enddo

       vari = vari_all (EOF,N,NEOF)
       change = 100.0D0*DABS(glovar-vari)/glovar
       WRITE(*,*) 'VARIMAX = ', vari, '  change = ', change
       if (change.GT.0.1) GOTO 10


    ! ... Sort EOFs by explained variance
    ! ...
    DEALLOCATE (PAIR)
    DEALLOCATE (PPAIR)

    ALLOCATE (xvari(NEOF))

    do k=1,NEOF
      xvari(k) = DOT_PRODUCT (EOF(:,k),EOF(:,k))
    enddo

    ts = TRANSPOSE(PC)
    CALL indexx (xvari,irang)

  end subroutine rotate_eofs
  ! ...
  ! ===================================================================
  ! ...
  subroutine ROTATE (N1,N2,EOF,N,NEOF,PAIR,angle)

    integer N1,N2,N,NEOF
    real(dp) EOF(N,NEOF)
    real(dp) angle
    real(dp) PAIR(N,2)

    integer i

    ! ... ROTATE

    do i=1,N
      PAIR(i,1)=  DCOS(angle)*EOF(i,N1) + DSIN(angle)*EOF(i,N2)
      PAIR(i,2)= -DSIN(angle)*EOF(i,N1) + DCOS(angle)*EOF(i,N2)
    enddo

  end subroutine ROTATE
  ! ...
  ! ===================================================================
  ! ...
  subroutine VARIMAX (PAIR,N,VARIM)

    integer N
    real(dp) varim
    real(dp) PAIR(N,2)

    integer i
    real(dp) sum11,sum12,sum21,sum22

    ! ... VARIMAX

    sum11 = 0.0D0
    sum12 = 0.0D0
    sum21 = 0.0D0
    sum22 = 0.0D0
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
  function vari_all (EOF,N,NEOF)

    real(dp) vari_all

    integer N,NEOF
    real(dp) EOF(N,NEOF)

    integer i,k
    real(dp) vari,xsum1,xsum2,tmp

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

    vari_all = vari/(N*N)

  end function vari_all
  ! ...
  ! ===================================================================
  ! ...
  subroutine indexx(arr,indx)

  ! ... Subroutine indexx from Numerical Recipes in Fortran.
  ! ... The Art of Scientific Computing. Press et al., 1992.
  ! ...

    real(dp), dimension(:), intent(in)          :: arr
    integer, dimension(size(arr)), intent(out)  :: indx

    ! ... Local varibles:
    ! ...
    integer, parameter                          :: M = 7
    integer                                     :: n,i,indxt,r,itemp
    integer                                     :: j,jstack,k,l
    integer, dimension(size(arr))               :: istack
    real(dp)                                    :: a

    n  =  size(arr)

    do j = 1,n
      indx(j) = j
    enddo

    jstack = 0
    l = 1
    r = n
    do 
      if (r-l.lt.M) then
        do j = l+1,r
          indxt = indx(j)
          a = arr(indxt)
          do i = j-1,1,-1
            if (arr(indx(i)).le.a) exit
            indx(i+1) = indx(i)
          enddo
          indx(i+1) = indxt
        enddo
        if (jstack.eq.0) return
        r = istack(jstack)
        l = istack(jstack-1)
        jstack = jstack-2
      else
        k = (l+r)/2
        itemp = indx(k)
        indx(k) = indx(l+1)
        indx(l+1) = itemp
        if (arr(indx(l+1)).gt.arr(indx(r))) then
          itemp     = indx(l+1)
          indx(l+1) = indx(r)
          indx(r)   = itemp
        endif
        if (arr(indx(l)).gt.arr(indx(r))) then
          itemp   = indx(l)
          indx(l) = indx(r)
          indx(r) = itemp
        endif
        if (arr(indx(l+1)).gt.arr(indx(l))) then
          itemp     = indx(l+1)
          indx(l+1) = indx(l)
          indx(l)   = itemp
        endif
        i = l + 1
        j = r
        indxt = indx(l+1)
        a = arr(indxt)
        do
          do
            i = i + 1
            if (arr(indx(i)).ge.a) exit
          enddo
          do 
            j = j - 1
            if (arr(indx(j)).le.a) exit
          enddo
          if (j.lt.i) exit
          itemp   = indx(i)
          indx(i) = indx(j)
          indx(j) = itemp
        enddo
        indx(l+1) = indx(j)
        indx(j)   = indxt
        jstack    = jstack + 2
        if (r-i+1.ge.j-l) then
          istack(jstack) = r
          istack(jstack-1) = i
          r = j - 1
        else
          istack(jstack) = j-1
          istack(jstack-1) = l
          l = i
        endif
      endif
    enddo

    end subroutine indexx
  ! ...
  ! ===================================================================
  ! ...
end module module_math
