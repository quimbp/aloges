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
! A generic hash table mapping STRING keys -> ANY Fortran value (class(*)) !
! - Open addressing with linear probing                                    !
! - Automatic resize at load factor ~0.7                                   !
! - String keys are case-sensitive                                         !
!                                                                          !
! Requires: Fortran 2008 (for deferred-length character & allocate)        !
!                                                                          !
! - dict_init                                                              !
! - dict_set                                                               !
! - dict_get                                                               !
! - dict_has_key                                                           !
! - dict_erase                                                             !
! - dict_size                                                              !
! - dict_clear                                                             !
! - dict_keys                                                              !
! -------------------------------------------------------------------------!

module module_dictionary

use module_types

implicit none
private
public :: type_dict
public :: dict_init, dict_set, dict_get, dict_has_key
public :: dict_erase, dict_size, dict_clear, dict_keys

integer, parameter :: EMPTY=0, USED=1, DELETED=-1

! Per-slot boxed string key (deferred-length)
type :: strbox
  character(len=:), allocatable :: s
end type strbox

! Per-slot boxed polymorphic value
type :: box
  class(*), allocatable :: data
end type box

type type_dict
  integer :: capacity = 0
  integer :: count    = 0
  integer, allocatable :: state(:)   ! EMPTY / USED / DELETED
  type(strbox), allocatable :: keys(:)
  type(box),    allocatable :: values(:)
  contains
    procedure :: size   => dict_size
    procedure :: clear  => dict_clear
    procedure :: has    => dict_has_key
    procedure :: set    => dict_set
    procedure :: get    => dict_get
    procedure :: erase  => dict_erase
end type type_dict

contains
    !===============================
    subroutine dict_init(DICT, cap)
        type(type_dict), intent(out) :: DICT
        integer, intent(in), optional :: cap
        integer :: n
        n = merge(cap, 16, present(cap))
        call allocate_storage(DICT, n)
    end subroutine dict_init

    !===============================
    subroutine allocate_storage(DICT, n)
        type(type_dict), intent(inout) :: DICT
        integer, intent(in) :: n
        DICT%capacity = max(4, n)
        DICT%count    = 0
        allocate(DICT%state(DICT%capacity))
        allocate(DICT%keys(DICT%capacity))
        allocate(DICT%values(DICT%capacity))
        DICT%state = EMPTY
        ! keys(:)%s and values(:)%data start unallocated
    end subroutine allocate_storage

    !===============================
    pure integer function hash_string(s) result(h)
        character(len=*), intent(in) :: s
        integer :: i
        h = 0
        do i = 1, len_trim(s)
            h = iachar(s(i:i)) + 31*h   ! simple multiplicative hash
        end do
        if (h < 0) h = -h
    end function hash_string

    !===============================
    subroutine ensure_capacity(DICT)
        type(type_dict), intent(inout) :: DICT
        real(dp) :: load
        if (DICT%capacity == 0) then
            call allocate_storage(DICT, 16)
            return
        end if
        load = real(DICT%count)/real(DICT%capacity)
        if (load > 0.70) call rehash(DICT, max(4, 2*DICT%capacity))
    end subroutine ensure_capacity

    !===============================

    subroutine rehash(DICT, new_cap)
        type(type_dict), intent(inout) :: DICT
        integer, intent(in) :: new_cap     
        type(type_dict) :: newd
        integer :: i

        call allocate_storage(newd, new_cap)
        do i = 1, size(DICT%state)
            if (DICT%state(i) == USED) then
                call DICT_set(newd, DICT%keys(i)%s, DICT%values(i)%data)
            end if
        end do
        DICT = newd
    end subroutine rehash

    !===============================
    subroutine find_slot(DICT, key, found, pos)
        !! Linear probing. If key is found => found=.true., pos=index
        !! Else => found=.false., pos=preferred insertion index (first tombstone or empty)
        type(type_dict), intent(in) :: DICT
        character(len=*),    intent(in) :: key
        logical,             intent(out):: found
        integer,             intent(out):: pos

        integer :: h, start, i, first_deleted
        if (DICT%capacity == 0) then
            found = .false.; pos = 0; return
        end if

        h = mod(hash_string(key), DICT%capacity) + 1  ! 1-based
        start = h
        first_deleted = 0

        do
            select case (DICT%state(h))
            case (EMPTY)
                found = .false.
                if (first_deleted /= 0) then
                    pos = first_deleted
                else
                    pos = h
                end if
                return
            case (USED)
                if (allocated(DICT%keys(h)%s)) then
                    if (DICT%keys(h)%s == key) then
                        found = .true.; pos = h; return
                    end if
                end if
            case (DELETED)
                if (first_deleted == 0) first_deleted = h
            end select
            h = h + 1
            if (h > DICT%capacity) h = 1
            if (h == start) then
                ! Table is full; pick first_deleted or h
                found = .false.
                pos = merge(first_deleted, h, first_deleted /= 0)
                return
            end if
        end do
    end subroutine find_slot

    !===============================
    subroutine dict_set(DICT, key, value)
        class(type_dict), intent(inout) :: DICT
        character(len=*),      intent(in)    :: key
        class(*),              intent(in)    :: value

        logical :: found
        integer :: idx

        call ensure_capacity(DICT)
        call find_slot(DICT, key, found, idx)

        if (found) then
            if (allocated(DICT%values(idx)%data)) deallocate(DICT%values(idx)%data)
            allocate(DICT%values(idx)%data, source=value)
            ! key stays same
        else
            DICT%keys(idx)%s = key  ! auto-allocates to right length
            if (allocated(DICT%values(idx)%data)) deallocate(DICT%values(idx)%data)
            allocate(DICT%values(idx)%data, source=value)
            DICT%state(idx) = USED
            DICT%count = DICT%count + 1
        end if
    end subroutine dict_set

    !===============================
    subroutine dict_get(DICT, key, value, found)
        class(type_dict), intent(in)    :: DICT
        character(len=*),      intent(in)    :: key
        class(*), allocatable, intent(out)   :: value
        logical,      optional, intent(out)  :: found

        logical :: f
        integer :: idx

        if (allocated(value)) deallocate(value)

        if (DICT%capacity == 0) then
            if (present(found)) found = .false.
            return
        end if

        call find_slot(DICT, key, f, idx)
        if (f .and. DICT%state(idx) == USED) then
            if (allocated(DICT%values(idx)%data)) then
                allocate(value, source=DICT%values(idx)%data)
                if (present(found)) found = .true.
                return
            end if
        end if
        if (present(found)) found = .false.
    end subroutine dict_get

    !===============================
    logical function dict_has_key(DICT, key)
        class(type_dict), intent(in) :: DICT
        character(len=*),      intent(in) :: key
        logical :: f
        integer :: idx
        if (DICT%capacity == 0) then
            dict_has_key = .false.; return
        end if
        call find_slot(DICT, key, f, idx)
        dict_has_key = f .and. DICT%state(idx) == USED
    end function dict_has_key

    !===============================
    logical function dict_erase(DICT, key)
        class(type_dict), intent(inout) :: DICT
        character(len=*),      intent(in)    :: key
        logical :: f
        integer :: idx
        if (DICT%capacity == 0) then
            dict_erase = .false.; return
        end if
        call find_slot(DICT, key, f, idx)
        if (f .and. DICT%state(idx) == USED) then
            if (allocated(DICT%values(idx)%data)) deallocate(DICT%values(idx)%data)
            if (allocated(DICT%keys(idx)%s))     deallocate(DICT%keys(idx)%s)
            DICT%state(idx) = DELETED
            DICT%count = DICT%count - 1
            dict_erase = .true.
        else
            dict_erase = .false.
        end if
    end function dict_erase

    !===============================
    integer function dict_size(DICT)
        class(type_dict), intent(in) :: DICT
        dict_size = DICT%count
    end function dict_size

    !===============================
    subroutine dict_clear(DICT, new_cap)
        class(type_dict), intent(inout) :: DICT
        integer, intent(in), optional :: new_cap
        integer :: n
        n = merge(new_cap, max(16, DICT%capacity), present(new_cap))
        ! Deallocate all current content by re-allocating fresh storage
        if (allocated(DICT%state))  deallocate(DICT%state)
        if (allocated(DICT%keys))   deallocate(DICT%keys)
        if (allocated(DICT%values)) deallocate(DICT%values)
        call allocate_storage(DICT, n)
    end subroutine dict_clear

    !===============================
    subroutine dict_keys(DICT, out_keys)
        class(type_dict), intent(in)  :: DICT
        character(len=:), allocatable, intent(out) :: out_keys(:)
        integer :: i, k
        allocate(character(len=0) :: out_keys(0))
        if (DICT%count == 0) return
        allocate(character(len=1) :: out_keys(DICT%count))  ! temp length; will reset below
        k = 0
        do i = 1, DICT%capacity
            if (DICT%state(i) == USED .and. allocated(DICT%keys(i)%s)) then
                k = k + 1
                out_keys(k) = DICT%keys(i)%s  ! reallocation on assignment adjusts length
            end if
        end do
        if (k /= DICT%count) then
            ! Shrink if tombstones existed
            out_keys = out_keys(:k)
        end if
    end subroutine dict_keys

end module module_dictionary
