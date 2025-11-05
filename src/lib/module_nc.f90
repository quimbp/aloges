! ======================================================================== !
! ALOGES PROJECT                                                           !
! Quim Ballabrera, April 2022                                              !
! Institut de Ciencies del Mar, CSIC                                       !
! Last Modified: 2024-04-16                                                !
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
! NF90 core variable types                                                 !
! NetCDF Constant	Description		Fortran Type               !
! NF90_BYTE		Signed 1-byte int	integer(kind=1)     1      !
! NF90_CHAR		ASCII character		character           2      !
! NF90_SHORT		2-byte integer		integer(kind=2)     3      !
! NF90_INT		4-byte integer		integer(kind=4)     4      !
! NF90_FLOAT		4-byte real		real(kind=4)        5      !
! NF90_DOUBLE		8-byte real		real(kind=8)        6      !
!                                                                          !
! NF90 extended variable types                                             !
! NetCDF Constant	Description		Fortran Type               !
! NF90_UBYTE		Unsigned 1-byte int	integer(kind=1)     7      !
! NF90_USHORT		Unsigned 2-byte int	integer(kind=2)     8      !
! NF90_UINT		Unsigned 4-byte int	integer(kind=4)     9      !
! NF90_INT64		Signed 8-byte integer	integer(kind=8)    10      !
! NF90_UINT64		Unsigned 8-byte integer	integer(kind=8)    11      !
! NF90_STRING		Variable-length string	character(len=*)   12      !
!                                                                          !
! List of routines:                                                        !
! - nc_open                                                                !
! - nc_dump                                                                !
! - nc_copyatts                                                            !
! - nc_error                                                               !
! - nc_axis_id                                                             !
! - nc_var_exists                                                          !
! -------------------------------------------------------------------------!

module module_nc

use netcdf
use module_types
use module_tools, only: crash

implicit none (type, external)

type type_nc_dimension
  character(len=maxlen)                              :: name
  integer                                            :: len
end type type_nc_dimension

type type_nc_attribute
  character(len=maxlen)                              :: name
  integer                                            :: type
  integer                                            :: len
  class(*), allocatable                              :: value(:)
!  logical, dimension(:), allocatable                 :: lval
!  integer, dimension(:), allocatable                 :: ival
!  real(dp), dimension(:), allocatable                :: dval
!  character(len=:), allocatable                      :: sval
  contains
    procedure                                        :: get => nc_attribute_get_dp
end type type_nc_attribute

type type_nc_variable
  character(len=maxlen)                              :: name
  integer                                            :: type
  integer                                            :: ndims
  integer                                            :: natts
  logical                                            :: with_missing = .False.
  real(dp)                                           :: missing_value
  real(dp)                                           :: add_offset = 0.0D0
  real(dp)                                           :: scale_factor = 1.0D0
  integer, dimension(:), pointer                     :: dimids
  integer, dimension(:), pointer                     :: shape
  type(type_nc_attribute), dimension(:), pointer     :: attribute
end type type_nc_variable

type type_dataset
  character(len=maxlen)                              :: filename
  integer                                            :: err
  integer                                            :: fid
  integer                                            :: ndims
  integer                                            :: nvars
  integer                                            :: natts
  integer                                            :: unlimid
  character(len=maxlen), dimension(:), allocatable   :: dimlist
  character(len=maxlen), dimension(:), allocatable   :: varlist
  type(type_nc_dimension), dimension(:), allocatable :: dimension
  type(type_nc_variable), dimension(:), allocatable  :: variable
  type(type_nc_attribute), dimension(:), allocatable :: attribute

  contains
    procedure                   :: open          => nc_open
    procedure                   :: dump          => nc_dump
    procedure                   :: read1D        => nc_variable_read1d
    procedure                   :: read2D        => nc_variable_read2d
    procedure                   :: size          => nc_variable_size  
    procedure                   :: varid         => nc_variable_id
end type type_dataset

contains
! ...
! =====================================================================
! =====================================================================
! ...
  subroutine nc_open(SD,filename) 

    class(type_dataset), intent(inout)               :: SD
    character(len=*), intent(in)                :: filename

    ! ... Local variables
    ! ...
    integer err,fid,ndims,nvars,natts,unlimid
    integer i,j,dlen,vtype,dimids(100)
    character(len=maxlen) word,attname


    err = NF90_OPEN(filename,0,fid)
    SD%filename = trim(filename)
    SD%fid      = fid
    SD%err      = err
    if (err.ne.NF90_NOERR) then
      SD%err = err
      return
    endif

    err = NF90_INQUIRE(fid,ndims,nvars,natts,unlimid)
    call nc_error(err,'In NC_OPEN: unable to open file')
    SD%ndims   = ndims
    SD%nvars   = nvars
    SD%natts   = natts
    SD%unlimid = unlimid
    allocate(SD%dimlist(ndims))
    allocate(SD%varlist(nvars))

    if (natts.gt.0) then
      allocate(SD%attribute(natts))
      do i=1,natts
        SD%attribute(i) = nc_read_attribute(fid,0,i)
      enddo
    endif

    if (ndims.gt.0) then
      allocate(SD%dimension(ndims))
      do i=1,ndims
        err = NF90_INQUIRE_DIMENSION(fid,i,word,dlen)
        SD%dimension(i)%name = trim(word)
        SD%dimension(i)%len  = dlen
        SD%dimlist(i) = trim(word)
      enddo
    endif

    if (nvars.gt.0) then
      allocate(SD%variable(nvars))
      do i=1,nvars
        err = NF90_INQUIRE_VARIABLE(fid,i,word,vtype,ndims,dimids,natts)
        SD%variable(i)%name = trim(word)
        SD%variable(i)%type = vtype
        SD%variable(i)%ndims = ndims
        SD%varlist(i) = trim(word)
        allocate(SD%variable(i)%dimids(ndims))
        allocate(SD%variable(i)%shape(ndims))
        SD%variable(i)%dimids(:) = dimids(1:ndims)
        do j=1,ndims
          SD%variable(i)%shape(j) = SD%dimension(dimids(j))%len
        enddo
        SD%variable(i)%natts = natts
        allocate(SD%variable(i)%attribute(natts))
        do j=1,natts
          SD%variable(i)%attribute(j) = nc_read_attribute(fid,i,j)
        enddo
        ! ... Retrieve the missing value, offset and scale factors:
        ! ...
        do j=1,natts
          attname = trim(SD%variable(i)%attribute(j)%name)
          if (attname.eq.'add_offset') SD%variable(i)%add_offset = SD%variable(i)%attribute(j)%get()
          if (attname.eq.'scale_factor') SD%variable(i)%scale_factor = SD%variable(i)%attribute(j)%get()
          !if (attname.eq.'add_offset') SD%variable(i)%add_offset = SD%variable(i)%attribute(j)%dval(1)
          !if (attname.eq.'scale_factor') SD%variable(i)%scale_factor = SD%variable(i)%attribute(j)%dval(1)
          if (attname.eq.'_FillValue'.or.attname.eq.'missing_value') then
            SD%variable(i)%with_missing  = .True.
            SD%variable(i)%missing_value = SD%variable(i)%attribute(j)%get()
            !SD%variable(i)%missing_value = SD%variable(i)%attribute(j)%dval(1)
          endif
        enddo
      enddo
    endif

  end subroutine nc_open
  ! ...
  ! ==================================================================
  ! ...
  subroutine nc_close(SD)

    class(type_dataset), intent(inout)               :: SD
    integer err

    err = NF90_CLOSE(SD%fid)
    call nc_error(err,'nc_close - unable to close file')

  end subroutine nc_close
  ! ...
  ! ==================================================================
  ! ...
  subroutine nc_dump(SD)

    class(type_dataset), intent(inout)               :: SD

    ! ... Local variables
    ! ...
    integer i,j
    character(len=maxlen) word

    ! ... Remove the extension and print basename
    ! ...
    do i=len_trim(SD%filename),1,-1
      if (SD%filename(i:i).eq.'.') exit
    enddo
    i = i - 1
    do j=len_trim(SD%filename),1,-1
      if (SD%filename(j:j).eq.'/') exit
    enddo
    j = j + 1


    write(6,'(T1,A)') 'netcdf '// trim(SD%filename(j:i)) // " {"

    write(6,'(T1,A)') 'dimensions:'
    do i=1,SD%ndims
      if (i.eq.SD%unlimid) then
        write(6,'(T9,A," = UNLIMITED ; // ( ",I4, " currently )" )') trim(SD%dimension(i)%name), &
                 SD%dimension(i)%len
      else
        write(6,'(T9,A," = ",I4, " ;")') trim(SD%dimension(i)%name), &
                 SD%dimension(i)%len
      endif
    enddo


    write(6,'(T1,A)') 'variables:'
    do i=1,SD%nvars
      select case (SD%variable(i)%TYPE)
      case (NF90_CHAR)
        word = 'char'
      case (NF90_SHORT)
        word = 'short'
      case (NF90_INT)
        word = 'int'
      case (NF90_FLOAT)
        word = 'float'
      case (NF90_DOUBLE)
        word = 'double'
      case default
        stop 'TYPE'
      end select

      word = trim(word) // ' ' // trim(SD%variable(i)%name)

      if (SD%variable(i)%ndims.eq.0) then
        word = trim(word) // ' ;' 
      else
        word = trim(word) // '(' 
        do j=SD%variable(i)%ndims,1,-1
          if (j.eq.SD%variable(i)%ndims) then
            word = trim(word) // trim(SD%dimension(SD%variable(i)%dimids(j))%name)
          else
            word = trim(word) // ' ' // trim(SD%dimension(SD%variable(i)%dimids(j))%name)
          endif
          if (j.gt.1) word = trim(word)//','
        enddo
        word = trim(word) // ') ;' 
      endif

      write(6,'(T9,A)') trim(word)
      do j=1,SD%variable(i)%natts
        call nc_print_attribute(SD%variable(i)%name,SD%variable(i)%attribute(j))
      enddo

    enddo

    write(6,*)
    write(6,'(T1,A)') '// global attributes: '
    do i=1, SD%natts
      call nc_print_attribute('',SD%attribute(i))
    enddo

    write(6,'(T1,"}")')

  end subroutine nc_dump
  ! ...
  ! ==================================================================
  ! ...
  function nc_read_attribute(fid,varid,attid) result (ATT)

    integer, intent(in)                      :: fid
    integer, intent(in)                      :: varid
    integer, intent(in)                      :: attid
    type(type_nc_attribute)                  :: ATT

    ! ... Local variables
    ! ...
    integer err,atype,alen
    integer(kind=1), allocatable         :: lval(:)
    integer, allocatable                 :: ival(:)
    real(dp), allocatable                :: dval(:)
    character(len=:), allocatable        :: sval
    character(len=maxlen) word

    err = NF90_INQ_ATTNAME(fid,varid,attid,word)
    call nc_error(err,'Attribute not found')

    err = NF90_INQUIRE_ATTRIBUTE(fid,varid,word,atype,alen)
    ATT%name = trim(word)
    ATT%type = atype
    ATT%len  = alen
    if (allocated(ATT%value)) deallocate(ATT%value)

    if (atype.EQ.NF90_BYTE) then
      allocate(lval(alen))
      err = NF90_GET_ATT(fid,varid,word,lval)
      allocate(ATT%value,source=lval)
      stop 'NF90_BYTE !'
    endif
    if (atype.EQ.NF90_CHAR) then
      allocate(character(len=alen) :: sval)
      err = NF90_GET_ATT(fid,varid,word,sval)
      allocate(ATT%value(1),source=sval)
    endif
    if ((atype.EQ.NF90_SHORT).or.(atype.EQ.NF90_INT)) then
      allocate(ival(alen))
      err = NF90_GET_ATT(fid,varid,word,ival)
      allocate(ATT%value(alen),source=ival)
    endif
    if ((atype.EQ.NF90_FLOAT).or.(atype.EQ.NF90_DOUBLE)) then
      allocate(dval(alen))
      err = NF90_GET_ATT(fid,varid,word,dval)
      allocate(ATT%value(alen),source=dval)
    endif
    if (atype.EQ.NF90_STRING) then
      allocate(character(len=alen) :: sval)
      err = NF90_GET_ATT(fid,varid,word,sval)
      allocate(ATT%value(alen),source=sval)
    endif

  end function nc_read_attribute
  ! ...
  ! ==================================================================
  ! ...
  subroutine nc_print_attribute(Varname,ATT)

    character(len=*), intent(in)             :: Varname
    type(type_nc_attribute), intent(in)      :: ATT

    ! ... Local variables
    ! ...
    character(len=maxlen) word,s

    word = trim(Varname)//':'//trim(ATT%name)// " ="

    select type (v => ATT%value)
      type is (character(*))
        word = trim(word) // ' "' // trim(adjustl(v(1))) // '"'
      type is (integer)
        write(s,*) v
        word = trim(word) // ' ' // trim(adjustl(s)) 
      type is (real)
        write(s,*) v
        word = trim(word) // ' ' // trim(adjustl(s)) // 'f'
      type is (double precision)
        write(s,*) v
        word = trim(word) // ' ' // trim(adjustl(s)) // 'd'
      class default
        print*, ATT%name
        print*, ATT%type
        stop 'uncoded type'
    end select

!    select case (ATT%type)
!    case (NF90_CHAR)
!      word = trim(word) // ' "' // trim(adjustl(ATT%sval)) // '"'
!    case (NF90_FLOAT)
!      if (ATT%len.eq.1) then
!        write(s,*) ATT%dval(1)
!        word = trim(word) // ' ' // trim(adjustl(s)) // 'f'
!      else
!        print*, 'FLOAT more than one dimension'
!        stop 
!      endif
!    case (NF90_STRING)
!      word = trim(word) // ' "' // trim(adjustl(ATT%sval)) // '"'
!    case default
!      print*, 'NF90_CHAR : ', NF90_CHAR
!      print*, 'NF90_SHORT : ', NF90_SHORT
!      print*, 'NF90_FLOAT : ', NF90_FLOAT
!      print*, 'NF90_DOUBLE : ', NF90_DOUBLE
!      print*, ATT%type, NF90_FLOAT
!      stop 'uncoded type'
!    end select

    write(6,'(T17,A," ;")') trim(word)

  end subroutine nc_print_attribute
  ! ...
  ! ==================================================================
  ! ...
  function nc_attribute_get_dp(ATT) result(val)

    class(type_nc_attribute), intent(in)     :: ATT
    real(dp)                                 :: val

    val = 0.0_dp
    select type (v => ATT%value(1))
      type is (integer)
        val = real(v,dp)
      type is (real)
        val = real(v,dp)
      type is (double precision)
        val = v
    end select
      
  end function nc_attribute_get_dp
  ! ...
  ! ==================================================================
  ! ...
  function nc_variable_read1d(SD,varname,po,pf) result (X)

    implicit none

    class(type_dataset), intent(inout)               :: SD
    character(len=*), intent(in)                     :: varname
    integer, dimension(:), intent(in), optional      :: po, pf
    real(dp), dimension(:), allocatable              :: X
   
    ! ... Local variables
    ! ... 
    integer i,varid,dimid,n,err,ndims
    integer, dimension(:), allocatable               :: ppo,ppf
    real(dp) add_offset, scale_factor

    varid = -1
    do i=1,SD%nvars
      if (trim(varname).eq.trim(SD%variable(i)%name)) then
        varid = i
        exit
      endif
    enddo
    if (varid.lt.0) stop 'ERROR nc_variable_read1d: Variable not found'

    ndims = SD%variable(varid)%ndims
    if (ndims.ne.SD%variable(varid)%ndims) stop 'ERROR nc_variable_read1d: incompatible dimensions'

    allocate(ppo(ndims))
    allocate(ppf(ndims))

    do i=1,ndims
      dimid = SD%variable(varid)%dimids(i)
      ppo(i) = 1
      ppf(i) = SD%dimension(dimid)%len
    enddo

    if (present(po)) then
      if (size(po).ne.ndims) stop 'ERROR nc_variable_read1d: incompatible po'
      ppo(:) = po(:)
      if (present(pf)) then
        if (size(pf).ne.ndims) stop 'ERROR nc_variable_read1d: incompatible pf' 
        ppf(:) = pf(:)
      else
        stop 'ERROR nc_variable_read1d: pf required'
      endif
    endif

    n = 1
    do i=1,ndims
      if (ppf(i)-ppo(i)+1.gt.n) n = ppf(i)-ppo(i)+1
    enddo

    allocate(X(n))
  
    err = NF90_GET_VAR(SD%fid,varid,X,ppo,ppf)
    if (err.NE.0) stop 'ERROR reading variable'

    ! ... Process the read data:
    ! ...
    add_offset = SD%variable(varid)%add_offset
    scale_factor = SD%variable(varid)%scale_factor
    if (SD%variable(varid)%with_missing) then
      where(X.ne.SD%variable(varid)%missing_value) X = add_offset + scale_factor*X
    else
      X = add_offset + scale_factor*X
    endif

    deallocate(ppo,ppf)

  end function nc_variable_read1d
  ! ...
  ! ==================================================================
  ! ...
  function nc_variable_read2d(SD,varname,po,pf) result (X)

    implicit none

    class(type_dataset), intent(inout)               :: SD
    character(len=*), intent(in)                     :: varname
    integer, dimension(:), intent(in), optional      :: po, pf
    real(dp), dimension(:,:), allocatable            :: X
   
    ! ... Local variables
    ! ... 
    integer i,varid,dimid,n1,n2,err,ndims,ii,ni(2)
    integer, dimension(:), allocatable               :: ppo,ppf,dpp
    real(dp) add_offset, scale_factor

    varid = -1
    do i=1,SD%nvars
      if (trim(varname).eq.trim(SD%variable(i)%name)) then
        varid = i
        exit
      endif
    enddo
    if (varid.lt.0) stop 'ERROR nc_variable_read2d: Variable not found'

    ndims = SD%variable(varid)%ndims
    !print*, 'varname = ', trim(varname)
    !print*, 'varid = ', varid
    !print*, 'ndims = ', ndims 
    if (ndims.ne.SD%variable(varid)%ndims) stop 'ERROR nc_variable_read2d: incompatible dimensions'

    allocate(ppo(ndims))
    allocate(ppf(ndims))
    allocate(dpp(ndims))

    do i=1,ndims
      dimid = SD%variable(varid)%dimids(i)
      ppo(i) = 1
      ppf(i) = SD%dimension(dimid)%len
    enddo

    if (present(po)) then
      if (size(po).ne.ndims) stop 'ERROR nc_variable_read2d: incompatible po'
      ppo(:) = po(:)
      if (present(pf)) then
        if (size(pf).ne.ndims) stop 'ERROR nc_variable_read2d: incompatible pf' 
        ppf(:) = pf(:)
      else
        stop 'ERROR nc_variable_read2d: pf required'
      endif
    endif

    ! ... Get the size of the 2D grid
    ! ...
    ii = 0
    ni(:) = 0
    do i=1,ndims
      if (ppf(i).gt.1) then
         ii = ii + 1
         if (ii.gt.2) stop 'ERROR nc_variable_read2d: incompatible dpp'
         ni(ii) = ppf(i)
       endif
    enddo
    !print*, 'ni = ', ni
    n1 = ni(1)
    n2 = ni(2)

    allocate(X(n1,n2))
    err = NF90_GET_VAR(SD%fid,varid,X,ppo,ppf)
    if (err.NE.0) stop 'ERROR reading variable'

    ! ... Process the read data:
    ! ...
    add_offset = SD%variable(varid)%add_offset
    scale_factor = SD%variable(varid)%scale_factor
    if (SD%variable(varid)%with_missing) then
      where(X.ne.SD%variable(varid)%missing_value) X = add_offset + scale_factor*X
    else
      X = add_offset + scale_factor*X
    endif

    deallocate(ppo,ppf,dpp)


  end function nc_variable_read2d
  ! ...
  ! ==================================================================
  ! ...
  function nc_variable_size(SD,varname) result (D)

    class(type_dataset), intent(inout)               :: SD
    character(len=*), intent(in)                     :: varname
    integer, dimension(:), allocatable               :: D
   
    ! ... Local variables
    ! ... 
    integer i,varid,dimid,ndims

    varid = -1
    do i=1,SD%nvars
      if (trim(varname).eq.trim(SD%variable(i)%name)) then
        varid = i
        exit
      endif
    enddo
    if (varid.lt.0) stop 'ERROR nc_variable_size: Variable not found'

    ndims = SD%variable(varid)%ndims

    allocate(D(ndims))

    do i=1,ndims
      dimid = SD%variable(varid)%dimids(i)
      D(i) = SD%dimension(dimid)%len
    enddo

  end function nc_variable_size
  ! ...
  ! ===================================================================
  ! ...
  subroutine nc_copyatts (ver,id1,v1,id2,v2,natts,disregard)
  ! ... Copies the attributes from variable v1 in file id1 to the
  ! ... variable v2 in file id2.

    logical, intent(in)                                  :: ver
    integer, intent(in)                                  :: id1,v1,id2,v2
    integer, intent(out)                                 :: natts
    character(len=*), intent(in), optional               :: disregard

    ! ... Local variables:
    ! ...
    logical copy
    integer err,vtype,ndim,j,att_type,att_len
    character(len=120) name,att_name
    integer, dimension(100)  :: dimids
    character(len=maxlen) ldis

    character(len=4000)                     :: tmpt
    integer, dimension(:), allocatable      :: tmpi
    real(sp), dimension(:), allocatable     :: tmp4
    real(dp), dimension(:), allocatable     :: tmp8

    if (present(disregard)) then
      ldis = trim(disregard)
    else
      ldis = ' '
    endif

    ! ... Information from first file:
    ! ..
    if (v1.eq.NF90_GLOBAL) then
      err = NF90_INQUIRE (id1,nAttributes=natts)
      if (ver) write(*,*) 'Number of global attributes ', natts
    else
      err = NF90_INQUIRE_VARIABLE (id1,v1,name,vtype,ndim,dimids,natts)
      call nc_error (err,'CDF_COPYATTS Error: Unable to inquire variable')
      if (ver) write(*,*) 'Variable ', v1, ' has ', natts, ' attributes'
    endif

    do j=1,natts
      err = NF90_INQ_ATTNAME (id1,v1,j,att_name)
      copy = index(ldis,trim(att_name)).le.0
      if (copy) then
        err = NF90_INQUIRE_ATTRIBUTE (id1,v1,att_name,xtype=att_type)
        err = NF90_INQUIRE_ATTRIBUTE (id1,v1,att_name,len=att_len)
        if (att_type.eq.NF90_BYTE) then
          allocate (tmpi(att_len))
          err = NF90_GET_ATT(id1,v1,att_name,tmpi)
          err = NF90_PUT_ATT(id2,v2,TRIM(att_name),tmpi)
          deallocate (tmpi)
        endif
        if (att_type.EQ.NF90_CHAR.or.att_type.EQ.NF90_STRING) then
          if (att_len.gt.len(tmpt)) stop 'ERROR NC_COPYATTS: Increase size tmpt'
          err = NF90_GET_ATT(id1,v1,att_name,tmpt)
          call nc_error (err,'Unable to get text attribute')
          err = NF90_PUT_ATT(id2,v2,TRIM(att_name),tmpt(1:att_len))
          call nc_error (err,'Unable to write text attribute')
        endif
        if (att_type.eq.NF90_SHORT) then
          allocate (tmpi(att_len))
          err = NF90_GET_ATT(id1,v1,att_name,tmpi)
          CALL nc_error (err,'Unable to get short attribute')
          if (TRIM(att_name).NE.'_FillValue') then
            err = NF90_PUT_ATT(id2,v2,TRIM(att_name),tmpi)
          ELSE
            err = NF90_PUT_ATT(id2,v2,TRIM(att_name),tmpi)
          endif
          CALL nc_error (err,'Unable to write short attribute')
          deallocate (tmpi)
        endif
        if (att_type.EQ.NF90_INT) then
          allocate (tmpi(att_len))
          err = NF90_GET_ATT(id1,v1,att_name,tmpi)
          err = NF90_PUT_ATT(id2,v2,TRIM(att_name),tmpi)
          deallocate (tmpi)
        endif
        if (att_type.EQ.NF90_FLOAT) then
          allocate (tmp4(att_len))
          err = NF90_GET_ATT(id1,v1,att_name,tmp4)
          err = NF90_PUT_ATT(id2,v2,TRIM(att_name),tmp4)
          deallocate (tmp4)
        endif
        if (att_type.EQ.NF90_DOUBLE) then
          allocate (tmp8(att_len))
          err = NF90_GET_ATT(id1,v1,att_name,tmp8)
          err = NF90_PUT_ATT(id2,v2,TRIM(att_name),tmp8)
          deallocate (tmp8)
        endif
      endif
    enddo

    return
  end subroutine nc_copyatts
  ! ...
  ! ===================================================================
  ! ...
  subroutine nc_error(err,message)

    integer, intent(in)                            :: err
    character(len=*)                               :: message 

    if (err.eq.0) return 

    ! ... If here, it has been an error:
    ! ...
    call crash(trim(message)//' - '//trim(NF90_STRERROR(err)))

  end subroutine nc_error
  ! ...
  ! ===================================================================
  ! ...
  function nc_axis_id(fid,axis) result(idx)

    integer, intent(in)                                  :: fid
    character(len=1), intent(in)                         :: axis
    integer                                              :: idx

    ! ... Local variables
    ! ...
    integer err,ndims,nvars,ngatts,unlimid
    integer var,vtype,vndims,dimids(10),vnatts
    character(len=40) vname,uname
    character(len=1) vaxis

    idx = -1
    err = NF90_INQUIRE(fid,ndims,nvars,ngatts,unlimid)
    if (err.NE.NF90_NOERR) return

    do var=1,nvars
      vname = ''
      vaxis = ''
      err = NF90_INQUIRE_VARIABLE (fid,var,vname,vtype,vndims,dimids,vnatts)
      err = NF90_GET_ATT(fid,var,'axis',vaxis)
      if (err.EQ.NF90_NOERR) then
        vaxis = to_upper(vaxis)
        if (vaxis.eq.axis) then 
          idx   = var
          return
        endif
      endif
    enddo

    ! ... If we have not found the axis, we check the name
    ! ...
    if (axis.eq.'X') uname = 'LON'
    if (axis.eq.'Y') uname = 'LAT'
    if (axis.eq.'Z') uname = 'DEP'
    if (axis.eq.'T') uname = 'TIM'
 
    if (idx.eq.-1) then
      do var=1,nvars
        vname = ''
        err = NF90_INQUIRE_VARIABLE (fid,var,vname,vtype,vndims,dimids,vnatts)
        vname = to_upper(vname)
        if (vname(1:3).eq.uname) then
          idx = var
          return
        endif
      enddo
    endif

  contains

    function to_upper(A) result(t)
    ! ... Returns string in uppercase

      character(len=*), intent(in)   :: A
      character(len=len(A))          :: t

      ! ... Local variables
      ! ...
      integer i,j

      do i=1,len(A)
        j = iachar(A(i:i))
        if (j>= iachar("a") .and. j<=iachar("z") ) then
          t(i:i) = achar(iachar(A(i:i))-32)
        else
          t(i:i) = A(i:i)
        end if
      enddo

      return
    end function to_upper

  end function nc_axis_id
  ! ...
  ! ===================================================================
  ! ...
  function nc_axis_size(fid,idx) result(nx)

    integer, intent(in)                                  :: fid
    integer, intent(in)                                  :: idx
    integer                                              :: nx

    ! ... Local variables
    ! ...
    integer err
    integer vtype,vndims,dimids(10),vnatts,idi
    character(len=40) vname

    nx = 1
    if (idx.le.0) return

    err = NF90_INQUIRE_VARIABLE (fid,idx,vname,vtype,vndims,dimids,vnatts)
    if (err.NE.NF90_NOERR) return

    idi = dimids(1)
    err = NF90_INQUIRE_DIMENSION(fid,idi,len=nx)

  end function nc_axis_size
  ! ...
  ! ===================================================================
  ! ...
  function nc_axis_dim(fid,idx) result(idi)

    integer, intent(in)                                  :: fid
    integer, intent(in)                                  :: idx
    integer                                              :: idi

    ! ... Local variables
    ! ...
    integer err
    integer vtype,vndims,dimids(10),vnatts
    character(len=40) vname

    idi = -1
    if (idx.le.0) return

    err = NF90_INQUIRE_VARIABLE (fid,idx,vname,vtype,vndims,dimids,vnatts)
    if (err.NE.NF90_NOERR) return

    idi = dimids(1)
    return

  end function nc_axis_dim
  ! ...
  ! ===================================================================
  ! ...
  function nc_variable_id(SD,varname) result(varid)

    class(type_dataset), intent(inout)               :: SD
    character(len=*), intent(in)                     :: varname
    integer                                          :: varid
   
    ! ... Local variables
    ! ... 
    integer i

    varid = -1
    do i=1,SD%nvars
      if (trim(varname).eq.trim(SD%variable(i)%name)) then
        varid = i
        exit
      endif
    enddo
    if (varid.lt.0) stop 'ERROR nc_variable_id: Variable not found'

  end function nc_variable_id
  ! ...
  ! ===================================================================
  ! ...
  logical function nc_var_exists(ncid, varname)  
  ! Check if a variable exists

    integer, intent(in) :: ncid  
    character(len=*), intent(in) :: varname  

    ! ... Local variable
    ! ...
    integer :: varid, stat  
    
    stat = nf90_inq_varid(ncid, trim(varname), varid)  
    nc_var_exists = (stat == nf90_noerr)  

  end function nc_var_exists
  ! ...
  ! ===================================================================
  ! ...

end module module_nc
