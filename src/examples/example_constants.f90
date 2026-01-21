program example_constants

use module_constants ! or:  use aloges
use aloges, only: head_aloges
implicit none (type, external)

call head_aloges(fortran_version=.True.)

write(*,*) 'True, False    = ', True, False
write(*,*) 'One, two, half = ', one, two, half
write(*,*) 'Minus π        = ', minus * pi
write(*,*) 'nan            = ', nan()
write(*,*) '∞              = ', inf()
write(*,*) 'Earth radius   = ', constants%Earth_Radius
write(*,*) 'Constant g     = ', constants%Earth_gravity
write(*,*) 'Uppercase letters = ', uppercase_letters
write(*,*) 'Digits            = ', digits
write(*,*) 'Bell !!', ASCII_BEL

end program example_constants
