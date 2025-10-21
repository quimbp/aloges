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
! - ensrf
! - randname                                                               !
! - randseries                                                             !
! -------------------------------------------------------------------------!

module module_enkf

use, intrinsic :: IEEE_ARITHMETIC, ONLY : IEEE_VALUE, IEEE_QUIET_NAN
use module_types
use module_constants
use module_math
use module_linalg
use module_eofs

implicit none (type, external)

contains
! ...
! =====================================================================
! =====================================================================
! ...
  function emean (X) result(xm)

    real(dp), dimension(:,:), intent(in)  :: X
    real(dp), dimension(size(X,1))        :: xm

    integer m,n,i

    m = size(X,1)
    n = size(X,2)
    do i=1,m
      xm(i) = sum(X(i,:))/n
    enddo

  end function emean
  ! ...
  ! ===================================================================
  ! ...
  function ecov(X) result(cov)

    real(dp), dimension(:,:), intent(in)       :: X
    real(dp), dimension(size(X,1),size(X,1))   :: cov

    integer m,n,i,j,l
    real(dp) xsum
    real(dp), dimension(size(X,1))             :: xm

    m = size(X,1)
    n = size(X,2)

    xm = emean(X)
    do j=1,m
    do i=1,j
      xsum = 0.0D0
      do l=1,n
        xsum = xsum + (X(i,l)-xm(i))*(X(j,l)-xm(j))
      enddo
      cov(i,j) = xsum/(n-1)
      cov(j,i) = cov(i,j)
    enddo
    enddo

  end function ecov

  ! ... Kalman Filter:
  ! ...              xa = xf + K (yo - H xf) = xf + K d
  ! ...               K = Pf H^T [ H Pf H^T + R]^-1, gain
  ! ...               d = yo - H xf, innovation
  ! ...              Pa = ( I - K H) Pf, analysis error covariance
  ! ... Ensemble Kalman Filter:
  ! ...              Pf = 1/(ne-1) E E^T
  ! ...               K = E (HE)^T ( HE (HE)^T + (ne-1)R )^-1
  ! ...                 = E (HE)^T F^-1
  ! ...               F = HE (HE)^T + (ne-1) R, innovation covariance
  ! ...              Pa = E ( I - (HE)^T F^-1 (HE) ) E^T 
  ! ...                 = E ( I - T T^T ) E^T
  ! ...               E = ensemble perturbation matrix
  ! ...              HE = ensemble perturbation matrix in observational space
  ! ...               T = transform matrix
  ! ...
  ! ... Linear algebra:
  ! ...         If  A = U D U^T, with U^T U = 1, then
  ! ...             A^1/2 = U D^1/2 U^T
  ! ===================================================================

  subroutine seek(E,HE,yo,Robs,alpha,xa,test)
  ! Implementation of the "The SEEK filter method for data assimilation
  ! in oceanography: a synthesis" by Pierre Brasseur and Jacques Verron.
  ! Ocean Dynamics, 2006.
  ! https://doi.org/10.1007/s10236-006-0080-3
  ! 
  ! As input:
  ! E (nsys,ne)  - Ensemble array. nsys: system dimension, ne: ensemble size
  ! HE (nobs,ne) - Sampled ensemble values. nobs: number of observations
  ! yo (nobs)    - Observations 
  ! alpha        - Covariance inflation
  ! Robs (nobs)  - Observation error variances. Diagonal.
  ! test         - Optional flag to print additional diagnostics
  ! 
  ! As output:
  ! E (nsys,ne)  - Analysis ensemble.
  ! xa (nsys)    - Analysis field. 
  !
  ! We will assume that P = E E^T / (ne-1)
  !                     R = diag(Robs)
  !                     xa = mean(Pa)
  !

    real(dp), dimension(:,:), intent(inout)       :: E
    real(dp), dimension(:,:), intent(inout)       :: HE
    real(dp), dimension(:), intent(in)            :: yo
    real(dp), dimension(:), intent(in)            :: Robs
    real(dp), intent(in)                          :: alpha
    real(dp), dimension(size(E,1)), intent(out)   :: xa
    logical,  intent(in), optional                :: test

    ! ... Local variables:
    ! ...
    logical verb
    integer i,j,k,nsys,nobs,ne
    real(dp) xsum
    real(dp), dimension(size(E,1))                :: xf
    real(dp), dimension(size(yo))                 :: Hxf,innov
    real(dp), dimension(size(E,1),size(E,2))      :: Ea      ! Ensemble anom
    real(dp), dimension(size(E,2))                :: HETRi,c
    real(dp), dimension(size(E,2),size(E,2))      :: HETRHE
    real(dp), dimension(size(E,2),size(E,2))      :: UU,VV   ! SVD
    real(dp), dimension(size(E,2))                :: DD      ! SVD
    real(dp), dimension(size(E,2),size(E,2))      :: W

    write(*,*) '================================================='
    write(*,*) '===              SEEK filter                  ==='
    write(*,*) '================================================='
    write(*,*)
    write(*,*) 'System dimension: ', size(E,1)
    write(*,*) 'Ensemble size: ', size(E,2)
    write(*,*)

    ! ... Check verbosity
    ! ...
    verb = .false.
    if (present(test)) verb = test

    nsys = size(E,1)
    ne   = size(E,2)
    nobs = size(yo)
    if (verb) then
      write(*,*) 'TEST: nsys, ne, nobs = ', nsys, ne, nobs
    endif

    ! ... Check dimensions
    ! ...
    if (nsys.le.0) then
      write(*,*) 'Invalid system size'
      stop 1
    endif
    if (ne.le.0) then
      write(*,*) 'Empty ensemble'
      stop 1
    endif

    if (size(HE,1).ne.nobs) stop 'ERROR incompatible dimensions HE/yo'
    if (size(Robs).ne.nobs) stop 'ERROR incompatible dimensions Robs/yo'
    if (size(HE,2).ne.ne) stop 'ERROR incompatible dimensions E/HE'

    ! ... Ensemble mean (xf) and Ensemble anomalies (Ea)
    ! ...     Pf = E E^T / (ne-1) = Ea Ea^T if Ea = E/sqrt(ne-1)
    ! ...

    do i=1,nsys
      xf(i) = sum(E(i,:))/ne
      Ea(i,:) = sqrt(alpha)*(E(i,:) - xf(i))   ! Ensemble anomaly
    enddo

    if (verb) then
      write(*,*) 'TEST: xf = ', xf
      write(*,*) 'TEST: HE'
      do i=1,ne
        write(*,*) HE(:,i)
      enddo
    endif

    ! ... Innovation vector (normalized by Robs):
    ! ...
    do i=1,nobs
      Hxf(i)   = sum(HE(i,:))/ne
      HE(i,:)  = sqrt(alpha)*(HE(i,:) - Hxf(i))
      innov(i) = yo(i) - Hxf(i)
    enddo

    if (verb) then
      write(*,*) 'TEST: yobs   = ', yo
      write(*,*) 'TEST: Hxf    = ', Hxf
      write(*,*) 'TEST: innov  = ', innov
      write(*,*) 'TEST: Robs   = ', Robs
    endif

    ! ... Projections over the reduced rank:
    ! ... HS * Robs^-1 * innov
    ! ...
    do j=1,ne
      xsum = 0.0
      do i=1,nobs
        xsum = xsum + HE(i,j)*innov(i)/Robs(i)
      enddo
      HETRi(j) = xsum 
    enddo

    ! ... Innovation covariance matrix (reduced space):
    ! ... [ (ne-1)I + (HE)T R^-1 HE ]
    ! ...
    do j=1,ne
    do i=1,j
      xsum = 0.0D0
      do k=1,nobs
        xsum = xsum + HE(k,i)*HE(k,j)/Robs(k)
      enddo
      HETRHE(i,j) = xsum 
      HETRHE(j,i) = HETRHE(i,j)
    enddo
    enddo
    do i=1,ne
      HETRHE(i,i) = (ne-1.0D0) + HETRHE(i,i)
    enddo

    ! ... Solving [ (ne-1)I + (HE)T R^-1 HE ]^-1 (HE)^T R^-1 innov
    ! ... HETRHE c = HETRi
    ! ...

    ! ... Option using SVD decomposition: Robust but slow
    ! ...
    call svdcmp(HETRHE,UU,DD,VV)
    c = svbksb(UU,DD,VV,HETRi)

    ! ... Analysis:  xa = xf + Ea [ I + (HE)T R^-1 HE ]^-1 (HE)^T R^-1 innov
    ! ...               = xf + Ea c
    ! ...
    xa(:) = xf(:) + matmul(Ea,c)

    ! ... Posterior analysis
    ! ... W = HRTRHE^-1/2
    ! ...
    do i=1,ne
      if (abs(DD(i)).lt.1E-6*abs(DD(1))) then
        DD(i) = 0.0
      else
        DD(i) = 1.0D0 / sqrt(DD(i))
      endif
    enddo
      
    do j=1,ne
    do i=1,j
      xsum = 0.0D0
      do k=1,ne
        xsum = xsum + UU(i,k)*UU(j,k)*DD(k)   ! DD = 1/sqrt(DD)
      enddo
      W(i,j) = xsum
      W(j,i) = xsum
    enddo
    enddo
    E = matmul(Ea,W)

    do j=1,ne
      E(:,j) = E(:,j) + xa(:)
    enddo

  end subroutine seek
  ! ...
  ! ===================================================================
  ! ...
  subroutine ensrf(E,HE,yo,Robs,alpha,xa,test)
  ! Implementation of the "State-of-the-art stochastic data assimilation
  ! methods for high-dimensional non-Gaussian problems" by Vetar-Carvalho
  ! et al., Tellus A: Dynamic meteorology and oceanography, 2018. 
  ! https://doi.org/10.1080/16000870.2018.1445364
  ! 
  ! As input:
  ! E (nsys,ne)  - Ensemble array. nsys: system dimension, ne: ensemble size
  ! HE (nobs,ne) - Sampled ensemble values. nobs: number of observations
  ! yo (nobs)    - Observations 
  ! Robs (nobs)  - Observation error variances. Diagonal.
  ! alpha        - Covariance inflation
  ! test         - Optional flag to print additional diagnostics
  ! 
  ! As output:
  ! E (nsys,ne)  - Analysis ensemble.
  ! xa (nsys)    - Analysis field. 
  !
  ! We will assume that P = E E^T / (ne-1)
  !                     R = diag(Robs)
  !                     xa = mean(Pa)
  !
  ! This method is based on the equation:
  !
  !                Ea Ea^T = E^T [ I - (HE)^T F^-1 (HE) ] E
  !                      F = (HS) (HS)^T + (ne-1) R
  !

    real(dp), dimension(:,:), intent(inout)       :: E
    real(dp), dimension(:,:), intent(inout)       :: HE
    real(dp), dimension(:), intent(in)            :: yo
    real(dp), dimension(:), intent(in)            :: Robs
    real(dp), intent(in)                          :: alpha
    real(dp), dimension(size(E,1)), intent(out)   :: xa
    logical,  intent(in), optional                :: test

    ! ... Local variables:
    ! ...
    logical verb
    integer i,j,k,nsys,nobs,ne
    real(dp) xsum
    real(dp), dimension(size(E,1))                :: xf
    real(dp), dimension(size(yo))                 :: Hxf,innov
    real(dp), dimension(size(E,1),size(E,2))      :: Ea      ! Ensemble anom
    real(dp), dimension(size(E,2))                :: c
    real(dp), dimension(size(yo),size(yo))        :: F
    real(dp), dimension(size(yo),size(yo))        :: UF,VF   ! EVD F
    real(dp), dimension(size(yo))                 :: DF      ! EVD F
    real(dp), dimension(size(E,2),size(yo))       :: GR
    real(dp), dimension(size(E,2),size(yo))       :: UR      ! SVD GR
    real(dp), dimension(size(yo),size(yo))        :: VR      ! SVD GR
    real(dp), dimension(size(yo))                 :: DR      ! SVD GR
    real(dp), dimension(size(yo))                 :: wrk     
    real(dp), dimension(size(E,2),size(E,2))      :: W       ! Transform weight

    write(*,*) '================================================='
    write(*,*) '===       Ensemble Root Square Filter         ==='
    write(*,*) '================================================='
    write(*,*)
    write(*,*) 'System dimension: ', size(E,1)
    write(*,*) 'Ensemble size: ', size(E,2)
    write(*,*)


    ! ... Check verbosity
    ! ...
    verb = .false.
    if (present(test)) verb = test

    nsys = size(E,1)
    ne   = size(E,2)
    nobs = size(yo)
    if (verb) then
      write(*,*) 'TEST: nsys, ne, nobs = ', nsys, ne, nobs
    endif

    ! ... Check dimensions
    ! ...
    if (nsys.le.0) then
      write(*,*) 'Invalid system size'
      stop 1
    endif
    if (ne.le.0) then
      write(*,*) 'Empty ensemble'
      stop 1
    endif

    if (size(HE,1).ne.nobs) stop 'ERROR incompatible dimensions HE/yo'
    if (size(Robs).ne.nobs) stop 'ERROR incompatible dimensions Robs/yo'
    if (size(HE,2).ne.ne) stop 'ERROR incompatible dimensions E/HE'

    ! ... Ensemble mean (xf) and Ensemble anomalies (Ea)
    ! ...
    do i=1,nsys
      xf(i) = sum(E(i,:))/ne
      Ea(i,:) = sqrt(alpha)*(E(i,:) - xf(i))
    enddo

    if (verb) then
      write(*,*) 'TEST: xf = ', xf
      write(*,*) 'TEST: HE'
      do i=1,ne
        write(*,*) HE(:,i)
      enddo
    endif

    ! ... Innovation vector:
    ! ...
    do i=1,nobs
      Hxf(i) = sum(HE(i,:))/ne
      HE(i,:) = sqrt(alpha)*(HE(i,:) - Hxf(i))
      innov(i) = yo(i) - Hxf(i)
    enddo

    if (verb) then
      write(*,*) 'TEST: yobs  = ', yo
      write(*,*) 'TEST: Hxf   = ', Hxf
      write(*,*) 'TEST: innov = ', innov
      write(*,*) 'TEST: Robs  = ', Robs
    endif
    

    ! ... Innovation covariance matrix (observational space):
    ! ... F = [ (ne-1)*R + (HE) (HE)^T ]
    ! ...
    do j=1,nobs
    do i=1,j
      xsum = 0.0D0
      do k=1,ne
        xsum = xsum + HE(i,k)*HE(j,k)    ! Gram matrix
      enddo
      F(i,j) = xsum
      F(j,i) = xsum
    enddo
    enddo
    do i=1,nobs
      F(i,i) = (ne-1.0D0)*Robs(i) + F(i,i)
    enddo

    ! ... SVD decomposition to calculate the eigenvectors and 
    ! ... eigenvalues of the square, symmetric array F
    ! ... We transform D(i) to 1/D(i) if D(i) != 0
    ! ...
    call svdcmp(F,UF,DF,VF)
    do i=1,nobs
      if (DF(i).lt.1.0D-6*DF(1)) then
        DF(i) = 0.0D0
      else
        DF(i) = 1.0D0 / DF(i)
      endif
    enddo

    ! ... Analysis:  xa = xf + E (HE)^T [ (ne-1)*R + (HE) (HE)^T ]^-1 innov
    ! ...               = xf + E (HE)^T F^-1 innov
    ! ...               = xf + E (HE)^T UF DF^-1 UF^T innov
    ! ...               = xf + E c
    ! ...
    do i=1,nobs
      xsum = 0.0
      do j=1,nobs
        xsum = xsum + UF(j,i)*innov(j)
      enddo
      wrk(i) = xsum*DF(i)                   ! DF^-1 UF^T innov
    enddo
    wrk = matmul(UF,wrk)                    ! UF DF^-1 UF^T innov

    do i=1,ne
      xsum = 0.0D0
      do j=1,nobs
        xsum = xsum + HE(j,i)*wrk(j)
      enddo
      c(i) = xsum                           ! HE^T UF DF^-1 UF^T innov
    enddo

    ! ... Posterior state estimation
    ! ...
    xa(:) = xf(:) + matmul(Ea,c)


    ! ... Posterior ensemble anomalies:
    ! ... GR = (HE)^T UF DF^{-1/2}     (ne,nobs)
    ! ... 
    DO j=1,nobs
    DO i=1,ne
      xsum = 0.0D0
      do k=1,nobs
        xsum = xsum + HE(k,i)*UF(k,j)
      enddo
      GR(i,j) = xsum*sqrt(DF(j))               ! DF = 1/SQRT(DF)
    enddo
    enddo
      
    call svdcmp(GR,UR,DR,VR)

    ! ... W = UR (1-DR^2)^1/2 UR^T      (ne,ne)
    ! ... 
    do i=1,nobs
      DR(i) = 1.0D0 - DR(i)**2
      if (DR(i).GE.0.0D0) then
        DR(i) = sqrt(DR(i))
      else
        stop '1-DR^2 < 0'
      endif
    enddo

    do j=1,ne
    do i=1,j
      xsum = 0.0D0
      do k=1,nobs
        xsum = xsum + UR(i,k)*UR(j,k)*DR(k)
      enddo
      W(i,j) = xsum
      W(j,i) = xsum  
    enddo
    enddo
    E = matmul(Ea,W)

    do j=1,ne
      E(:,j) = E(:,j) + xa(:)
    enddo

  end subroutine ensrf
  ! ...
  ! ===================================================================
  ! ...
end module module_enkf
