use module_color ! or: use aloges
implicit none (type, external)

call print_colored('This text is in red', red)      ! red is a parameter constant
call print_colored('This text is in blue', blue)    ! red is a parameter constant
call print_colored('This text is in green', green)  ! green is a parameter constant

write(*,'(A)') red // 'Red' // reset // ' and ' // blue // ' Blue ' // reset

end

