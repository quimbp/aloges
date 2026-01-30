program example_color

use module_color ! or: use aloges
use aloges, only: head_aloges
implicit none (type, external)

call head_aloges(os_version=.True.,fortran_version=.True.)

call print_colored('This text is in red', red)      ! red is a parameter constant
call print_colored('This text is in blue', blue)    ! blue is a parameter constant
call print_colored('This text is in green', green)  ! green is a parameter constant
write(*,'(A)') red // 'Red' // reset // ' and ' // blue // ' Blue ' // reset

end program example_color

