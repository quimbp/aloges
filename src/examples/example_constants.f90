use module_constants ! or:  use aloges
implicit none (type, external)

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

end
