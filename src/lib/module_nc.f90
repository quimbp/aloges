module module_nc

use netcdf

implicit none

integer, parameter                                   :: dp = 8
integer, parameter                                   :: maxlen = 360


type type_nc_dimension
  character(len=maxlen)                              :: name
  integer                                            :: len
end type type_nc_dimension

type type_nc_variable
  character(len=maxlen)                              :: name
  integer                                            :: type
  integer                                            :: ndims
  integer                                            :: natts
  integer, dimension(:), pointer                     :: dimids
  type(type_nc_attribute), dimension(:), pointer     :: attribute

end type type_nc_variable

type type_nc_attribute
  character(len=maxlen)                              :: name
  integer                                            :: type
  integer                                            :: len
  logical, dimension(:), allocatable                 :: lval
  integer, dimension(:), allocatable                 :: ival
  real(dp), dimension(:), allocatable                :: dval
  character(len=maxlen)                              :: tval
end type type_nc_attribute

type type_dataset

  character(len=maxlen)                              :: filename
  integer                                            :: err
  integer                                            :: fid
  integer                                            :: ndims
  integer                                            :: nvars
  integer                                            :: natts
  integer                                            :: unlimid
  type(type_nc_dimension), dimension(:), allocatable :: dimension
  type(type_nc_variable), dimension(:), allocatable  :: variable
  type(type_nc_attribute), dimension(:), allocatable :: attribute

  contains
    procedure                   :: open          => nc_open
    procedure                   :: dump          => nc_dump
    procedure                   :: read1D        => nc_variable_read1d
    procedure                   :: read2D        => nc_variable_read2d
    procedure                   :: size          => nc_variable_size  

end type type_dataset

contains

  subroutine nc_open(SD,filename) 

    class(type_dataset), intent(inout)               :: SD
    character(len=maxlen), intent(in)                :: filename

    ! ... Local variables
    ! ...
    integer err,fid,ndims,nvars,natts,unlimid
    integer i,j,dlen,alen,vtype,atype,dimids(100)
    character(len=maxlen) word,word2


    err = NF90_OPEN(filename,0,fid)
    SD%filename = trim(filename)
    SD%fid      = fid
    SD%err      = err
    if (err.ne.NF90_NOERR) then
      SD%err = err
      return
    endif

    err = NF90_INQUIRE(fid,ndims,nvars,natts,unlimid)
    SD%ndims   = ndims
    SD%nvars   = nvars
    SD%natts   = natts
    SD%unlimid = unlimid
    if (err.ne.NF90_NOERR) then
      SD%err = err
      return
    endif

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
      enddo
    endif

    if (nvars.gt.0) then
      allocate(SD%variable(nvars))
      do i=1,nvars
        err = NF90_INQUIRE_VARIABLE(fid,i,word,vtype,ndims,dimids,natts)
        SD%variable(i)%name = trim(word)
        SD%variable(i)%type = vtype
        SD%variable(i)%ndims = ndims
        allocate(SD%variable(i)%dimids(ndims))
        SD%variable(i)%dimids(:) = dimids(1:ndims)
        SD%variable(i)%natts = natts
        allocate(SD%variable(i)%attribute(natts))
        do j=1,natts
          SD%variable(i)%attribute(j) = nc_read_attribute(fid,i,j)
        enddo
      enddo
    endif

  end subroutine nc_open
  ! ...
  ! ==================================================================
  ! ...
  subroutine nc_dump(SD)

    implicit none

    class(type_dataset), intent(inout)               :: SD

    ! ... Local variables
    ! ...
    integer i,j
    character(len=maxlen) word,dname

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
        write(6,'(T9,A," = UNLIMITTED ; // ( ",I4, " currently )" )') trim(SD%dimension(i)%name), &
                 SD%dimension(i)%len
      else
        write(6,'(T9,A," = ",I4, " ;")') trim(SD%dimension(i)%name), &
                 SD%dimension(i)%len
      endif
    enddo


    write(6,'(T1,A)') 'variables:'
    do i=1,SD%nvars
      select case (SD%variable(i)%TYPE)
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
    integer i,atype,alen,err
    character(len=maxlen) word,word2

    err = NF90_INQ_ATTNAME(fid,varid,attid,word)
    err = NF90_INQUIRE_ATTRIBUTE(fid,varid,word,atype,alen)
    ATT%name = trim(word)
    ATT%type = atype
    ATT%len  = alen

    if (atype.EQ.NF90_BYTE) then
      allocate(ATT%lval(alen))
      !err = NF90_GET_ATT(fid,0,word,SD%attribute(i)%lval)
      stop 'NF90_BYTE !'
    endif
    if (atype.EQ.NF90_CHAR) then
      err = NF90_GET_ATT(fid,varid,word,word2)
      ATT%tval = trim(word2)
    endif
    if ((atype.EQ.NF90_SHORT).or.(atype.EQ.NF90_INT)) then
      allocate(ATT%ival(alen))
      err = NF90_GET_ATT(fid,varid,word,ATT%ival)
    endif
    if ((atype.EQ.NF90_FLOAT).or.(atype.EQ.NF90_DOUBLE)) then
      allocate(ATT%dval(alen))
      err = NF90_GET_ATT(fid,varid,word,ATT%dval)
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

    select case (ATT%type)
    
    case (NF90_CHAR)
      word = trim(word) // ' "' // trim(adjustl(ATT%tval)) // '"'

    case (NF90_FLOAT)
      if (ATT%len.eq.1) then
        write(s,*) ATT%dval(1)
        word = trim(word) // ' ' // trim(adjustl(s)) // 'f'
      else
        print*, 'FLOAT more than one dimension'
        stop 
      endif

    case default
      print*, 'NF90_CHAR : ', NF90_CHAR
      print*, 'NF90_SHORT : ', NF90_SHORT
      print*, 'NF90_FLOAT : ', NF90_FLOAT
      print*, 'NF90_DOUBLE : ', NF90_DOUBLE
      print*, ATT%type, NF90_FLOAT
      stop 'uncoded type'
    end select

    write(6,'(T17,A," ;")') trim(word)

  end subroutine nc_print_attribute
  ! ...
  ! ==================================================================
  ! ...
  function nc_variable_read1d(SD,varname,po,pf) result (X)

    implicit none

    class(type_dataset), intent(inout)               :: SD
    character(len=*), intent(in)                     :: varname
    integer, dimension(:), intent(in), optional      :: po, pf
    real(8), dimension(:), allocatable               :: X
   
    ! ... Local variables
    ! ... 
    integer i,varid,dimid,n,err,ndims
    integer, dimension(:), allocatable               :: ppo,ppf

    varid = -1
    do i=1,SD%nvars
      if (trim(varname).eq.trim(SD%variable(i)%name)) then
        varid = i
        exit
      endif
    enddo
    if (varid.lt.0) stop 'ERROR nc_variable_read1d: Variable not found'

    ndims = SD%variable(varid)%ndims
    !print*, 'varname = ', trim(varname)
    !print*, 'varid = ', varid
    !print*, 'ndims = ', ndims 
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

    !print*, 'ppo = ', ppo
    !print*, 'ppf = ', ppf
    n = 1
    do i=1,ndims
      if (ppf(i)-ppo(i)+1.gt.n) n = ppf(i)-ppo(i)+1
    enddo

    allocate(X(n))
  
    err = NF90_GET_VAR(SD%fid,varid,X,ppo,ppf)
    if (err.NE.0) stop 'ERROR reading variable'


  end function nc_variable_read1d
  ! ...
  ! ==================================================================
  ! ...
  function nc_variable_read2d(SD,varname,po,pf) result (X)

    implicit none

    class(type_dataset), intent(inout)               :: SD
    character(len=*), intent(in)                     :: varname
    integer, dimension(:), intent(in), optional      :: po, pf
    real(8), dimension(:,:), allocatable             :: X
   
    ! ... Local variables
    ! ... 
    integer i,varid,dimid,n1,n2,err,ndims,ii,ni(2)
    integer, dimension(:), allocatable               :: ppo,ppf,dpp

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

  end function nc_variable_read2d

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

end module module_nc



program main

  use module_nc

  type(type_dataset)                                 :: SD
  character(len=maxlen)                              :: filename

  integer, dimension(:), allocatable                 :: n
  real(8), dimension(:), allocatable                 :: X
  real(8), dimension(:,:), allocatable               :: U

  filename = './latest.nc'
  call SD%open(filename)
  call SD%dump()

  X = SD%Read1D('lon_rho',[1],[10])
  print*, 'X = ', X
  X = SD%Read1D('lat_rho',[1],[10])
  print*, 'X = ', X

  X = SD%Read1D('temp',[350,200,1,1],[1,1,13,1])
  print*, 'temp = ', X

  n = SD%size('temp')
  print*, 'n = ', n

  U = SD%Read2D('temp',[350,200,1,1],[10,5,1,1])
  print*, 'U = ', U
  print*, SIZE(U,1), SIZE(U,2)
end program main
