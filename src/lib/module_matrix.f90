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
!                                                                          !
! List of routines:                                                        !
! - matrix_close                                                           !
! - matrix_copy                                                            !
! - matrix_open                                                            !
! - matrix_read                                                            !
! - matrix_save                                                            !
! - matrix_show                                                            !
! -------------------------------------------------------------------------!

module module_matrix

use module_types
use module_constants
use module_tools

implicit none (type, external)
private
public type_matrix, matrix_create

type type_matrix
  character(len=maxlen)                 :: filename = ''
  character(len=3)                      :: type = ''
  integer                               :: m=0
  integer                               :: n=0
  real(dp), dimension(:), allocatable   :: t
  real(dp), dimension(:,:), allocatable :: X
  integer                               :: iu = 0
  integer                               :: pos = 0
  character(len=maxlen)                 :: out_filename = ''
  character(len=3)                      :: out_type = ''
  integer                               :: out_iu = 0

  contains
    procedure                   :: open          => matrix_open
    procedure                   :: show          => matrix_show
    procedure                   :: read          => matrix_read
    procedure                   :: save          => matrix_save  
    procedure                   :: close         => matrix_close
    procedure                   :: copy          => matrix_copy

end type type_matrix

contains
! ...
! =====================================================================
! =====================================================================
! ...
  subroutine matrix_open(MAT,filename)

    class(type_matrix), intent(inout)    :: MAT
    character(len=*), intent(in)         :: filename

    ! ... Local variables
    ! ...
    integer m,n
    character(len=3) mtype

    call matrix_type(filename,mtype,m,n)
    if (mtype.eq.'') stop 'ERROR: unknown matrix type in matrix_open'

    MAT%filename = trim(filename)
    MAT%type = mtype
    MAT%m    = m
    MAT%n    = n
    MAT%pos  = 0                   ! File open, but not read
    MAT%iu   = unitfree()
    if (MAT%type.eq.'bin') open(MAT%iu,file=filename,form='unformatted',status='old')
    if (MAT%type.eq.'asc') open(MAT%iu,file=filename,form='formatted',status='old')

  end subroutine matrix_open
  ! ...
  ! ===================================================================
  ! ...
  subroutine matrix_show(MAT,label)

    class(type_matrix), intent(in)         :: MAT
    character(len=*), intent(in), optional :: label

    if (present(label)) write(*,*) trim(label)
    write(*,*) '> Filename: ', trim(MAT%filename)
    write(*,*) '> Type    : ', MAT%type
    write(*,*) '> M, N    : ', MAT%M, MAT%N

  end subroutine matrix_show
  ! ...
  ! ===================================================================
  ! ...
  subroutine matrix_type(filename,type,m,n)
    ! A matrix file can be:
    ! 'asc' - ASCII file
    ! 'bin' - Binary file
    ! 'cdf' - Netcdf file
    ! 
    character(len=*), intent(in)   :: filename
    character(len=*), intent(out)  :: type
    integer, intent(out)           :: m,n

    ! ... Local variables:
    ! ...
    integer iu,ios,i
    real(dp) :: time
    real(dp), dimension(:), allocatable :: A

    type = ''
    iu = unitfree()

    ! ... Try to read a binary file:
    ! ...
    open(iu,file=filename,form='unformatted',status='old',iostat=ios)
    if (ios.eq.0) then
      read(iu,iostat=ios) m,n
      if (ios.eq.0) then
        allocate(A(m))
        do i=1,n
          read(iu,err=100,end=100) time, A
        enddo
        close(iu)
        type = 'bin'
        deallocate(A)
        return
      endif
    endif
100 continue
    if (allocated(A)) deallocate(A)
    close(iu)

    ! ... Try to read an ascii file:
    ! ...
    open(iu,file=filename,form='formatted',status='old',iostat=ios)
    if (ios.eq.0) then
      read(iu,*,iostat=ios) m,n
      if (ios.eq.0) then
        allocate(A(m))
        do i=1,n
          read(iu,*,err=200,end=200) time, A
        enddo
        close(iu)
        type = 'asc'
        deallocate(A)
        return
      endif
    endif
200 continue
    if (allocated(A)) deallocate(A)
    close(iu)

    type = ''
    m = 0
    n = 0
       
  end subroutine matrix_type
  ! ...
  ! ===================================================================
  ! ...
  subroutine matrix_read(MAT,irec,T,X) 

    class(type_matrix), intent(inout)                 :: MAT
    integer, intent(in), optional                     :: irec
    real(dp), intent(out), optional                   :: T
    real(dp), dimension(MAT%m), intent(out), optional :: X

    logical full
    integer i,j
    
    i = count([present(irec),present(T),present(X)])
    if (i.eq.0) then
      full = .true.
    else if (i.eq.3) then
      full = .false.
    else
      stop 'Invalid arguments in matrix_read'
    endif

    if (full) then
      if (allocated(MAT%T)) deallocate(MAT%t)
      if (allocated(MAT%X)) deallocate(MAT%X)
      allocate(MAT%T(MAT%n))
      allocate(MAT%X(MAT%m,MAT%n))
      rewind(MAT%iu)
      if (MAT%type.eq.'bin') then
        read(MAT%iu)
        do j=1,MAT%n
          read(MAT%iu,err=100) MAT%T(j), MAT%X(:,j)
        enddo
      else if (MAT%type.eq.'asc') then
        read(MAT%iu,*)
        do j=1,MAT%n
          read(MAT%iu,*,err=100) MAT%T(j), MAT%X(:,j)
        enddo
      else
        stop 'ERROR: unknown matrix type in matrix_read'
      endif
    else
      rewind(MAT%iu)
      if (MAT%type.eq.'bin') then
        read(MAT%iu)
        do j=1,irec
          read(MAT%iu,err=100) T, X(:)
        enddo
      else if (MAT%type.eq.'asc') then
        read(MAT%iu,*)
        do j=1,irec
          read(MAT%iu,*,err=100) T, X(:)
        enddo
      else
        stop 'ERROR: unknown matrix type in matrix_read'
      endif
    endif
    return

100 continue
    stop 'ERROR in matrix_read'

  end subroutine matrix_read
  ! ...
  ! ===================================================================
  ! ...
  subroutine matrix_close(MAT) 

    class(type_matrix), intent(inout)                 :: MAT

    close(MAT%iu)
    MAT%pos = -1                             ! File closed

  end subroutine matrix_close
  ! ...
  ! ===================================================================
  ! ...
  subroutine matrix_save(MAT,filename,type,irec,T,X) 

    class(type_matrix), intent(inout)                 :: MAT
    character(len=*), intent(in)                      :: filename
    character(len=3), intent(in), optional            :: type
    integer, intent(in), optional                     :: irec
    real(dp), intent(in), optional                    :: T
    real(dp), dimension(MAT%m), intent(in), optional  :: X

    logical full
    integer i,ios
    
    i = count([present(irec),present(T),present(X)])
    if (i.eq.0) then
      full = .true.
    else if (i.eq.3) then
      full = .false.
    else
      stop 'Invalid arguments in matrix_read'
    endif

    MAT%out_filename = trim(filename)
    MAT%out_type = MAT%type
    if (present(type)) MAT%out_type = type

    MAT%out_iu = unitfree()

    if (full) then
      if (MAT%out_type.eq.'bin') then
        open(MAT%out_iu,file=MAT%out_filename,status='unknown',form='unformatted',iostat=ios)
        if (ios.ne.0) stop 'ERROR: opening file in matrix_save'
        write(MAT%out_iu) MAT%m, MAT%n
        do i=1,MAT%n
          write(MAT%out_iu) MAT%T(i), MAT%X(:,i)
        enddo
        close(MAT%out_iu)
      else if (MAT%out_type.eq.'asc') then
        open(MAT%out_iu,file=MAT%out_filename,status='unknown',form='formatted',iostat=ios)
        if (ios.ne.0) stop 'ERROR: opening file in matrix_save'
        write(MAT%out_iu,*) MAT%m, MAT%n
        do i=1,MAT%n
          write(MAT%out_iu,*) MAT%T(i), MAT%X(:,i)
        enddo
        close(MAT%out_iu)
      endif
    endif

  end subroutine matrix_save
  ! ...
  ! ===================================================================
  ! ...
  function matrix_copy (A) result(B)

    class(type_matrix), intent(in)               :: A
    type(type_matrix)                            :: B

    B%filename = ''
    B%type = A%type
    B%m = A%m
    B%n = A%n
    if (allocated(B%T)) deallocate(B%T)
    if (allocated(B%X)) deallocate(B%X)
    allocate(B%T(B%n))
    allocate(B%X(B%m,B%n))
    B%T(:) = A%T(:)
    B%X(:,:) = A%X(:,:)
   
  end function matrix_copy
  ! ...
  ! ===================================================================
  ! ...
  subroutine matrix_create(filename,A,type)

    character(len=*), intent(in)           :: filename
    real(dp), dimension(:,:), intent(in)   :: A
    character(len=*), intent(in), optional :: type
    integer iu,j
    character(len=3) out_type

    out_type = 'asc'  ! Default, ascii type
    if (present(type)) out_type = type

    if (out_type.eq.'asc') then
      open(newunit=iu,file=filename,status='unknown')
      rewind(iu)
      write(iu,*) size(A,1), size(A,2)
      do j=1,size(A,2)
        write(iu,*) real(j,dp), A(:,j)
      enddo   
    else
      open(newunit=iu,file=filename,status='unknown',form='unformatted')
      rewind(iu)
      write(iu) size(A,1), size(A,2)
      do j=1,size(A,2)
        write(iu) real(j,dp), A(:,j)
      enddo   
    endif
    close(iu)

  end subroutine matrix_create

  ! ...
  ! ===================================================================
  ! ...
end module module_matrix
