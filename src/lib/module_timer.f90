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
! - timer_start                                                            !
! - timer_stop                                                             !
! - timer_report                                                           !
! - timer_reset                                                            !
! - timer_cleanup                                                          !
! -------------------------------------------------------------------------!
! Use example:                                                             !
!                                                                          !
! call timer_start('Total Program')                                        !
! call timer_start('Open file')                                            !
! call open_file(filename)                                                 !
! call timer_stop('Open file')                                             !
!                                                                          !
! call timer_stop('Total Program')                                         !
! call timer_report()                                                      !
!                                                                          !
! -------------------------------------------------------------------------!

module module_timer

use module_types
implicit none

private
public :: timer_type, timer_start, timer_stop, timer_report, timer_reset, timer_cleanup

! ... Maximum number of timers
integer, parameter :: MAX_TIMERS = 100

type timer_type
  character(len=64) :: name = ""
  real(dp) :: total_time = 0.0_dp
  real(dp) :: start_time = 0.0_dp
  integer :: call_count = 0
  logical :: is_running = .false.
end type timer_type

! ... Global timer registry
type(timer_type), save :: timers(MAX_TIMERS)
integer, save :: n_timers = 0
logical, save :: timer_enabled = .true.

contains

  ! ====================================================================
  ! Start timing a named section
  ! ====================================================================
  subroutine timer_start(name)
    character(len=*), intent(in) :: name
    integer :: idx
    real(dp) :: current_time

    if (.not. timer_enabled) return

    ! Find or create timer
    idx = find_timer(name)
    if (idx == 0) then
      idx = create_timer(name)
      if (idx == 0) then
        write(*,*) 'WARNING: Maximum number of timers reached'
        return
      endif
    endif

    ! Check if already running
    if (timers(idx)%is_running) then
      write(*,*) 'WARNING: Timer "', trim(name), '" is already running'
      return
    endif

    ! Start timing
    call cpu_time(current_time)
    timers(idx)%start_time = current_time
    timers(idx)%is_running = .true.

  end subroutine timer_start

  ! ====================================================================
  ! Stop timing a named section
  ! ====================================================================
  subroutine timer_stop(name)
    character(len=*), intent(in) :: name
    integer :: idx
    real(dp) :: current_time, elapsed

    if (.not. timer_enabled) return

    ! Find timer
    idx = find_timer(name)
    if (idx == 0) then
      write(*,*) 'WARNING: Timer "', trim(name), '" not found'
      return
    endif

    ! Check if running
    if (.not. timers(idx)%is_running) then
      write(*,*) 'WARNING: Timer "', trim(name), '" is not running'
      return
    endif

    ! Stop timing
    call cpu_time(current_time)
    elapsed = current_time - timers(idx)%start_time
    timers(idx)%total_time = timers(idx)%total_time + elapsed
    timers(idx)%call_count = timers(idx)%call_count + 1
    timers(idx)%is_running = .false.

  end subroutine timer_stop

  ! ====================================================================
  ! Print timing report
  ! ====================================================================
  subroutine timer_report(unit_out, sort_by)
    integer, intent(in), optional :: unit_out  ! Output unit (default: stdout)
    character(len=*), intent(in), optional :: sort_by  ! 'time', 'calls', 'name', 'average'
    
    integer :: i, idx, unit_num, sorted_indices(MAX_TIMERS)
    real(dp) :: total_program_time, percentage, avg_time
    character(len=16) :: sort_method

    if (n_timers == 0) then
      write(*,*) 'No timers recorded'
      return
    endif

    ! Determine output unit
    if (present(unit_out)) then
      unit_num = unit_out
    else
      unit_num = 6  ! stdout
    endif

    ! Determine sort method
    if (present(sort_by)) then
      sort_method = sort_by
    else
      sort_method = 'time'
    endif

    ! Calculate total time
    total_program_time = sum(timers(1:n_timers)%total_time)

    ! Sort timers
    call sort_timers(sorted_indices, sort_method)

    ! Print header
    write(unit_num,*)
    write(unit_num,'(A)') repeat('=', 90)
    write(unit_num,'(A)') '                          TIMING REPORT'
    write(unit_num,'(A)') repeat('=', 90)
    write(unit_num,'(A32,A12,A10,A12,A12,A12)') &
      'Routine', 'Calls', 'Total(s)', 'Average(s)', 'Min/Call(ms)', 'Percent(%)'
    write(unit_num,'(A)') repeat('-', 90)

    ! Print each timer
    do i = 1, n_timers
      idx = sorted_indices(i)
      if (timers(idx)%call_count > 0) then
        avg_time = timers(idx)%total_time / real(timers(idx)%call_count, dp)
        percentage = 100.0_dp * timers(idx)%total_time / max(total_program_time, 1.0e-10_dp)
        
        write(unit_num,'(A32,I12,F10.3,F12.6,F12.3,F12.2)') &
          trim(timers(idx)%name), &
          timers(idx)%call_count, &
          timers(idx)%total_time, &
          avg_time, &
          avg_time * 1000.0_dp, &
          percentage
      endif
    enddo

    write(unit_num,'(A)') repeat('-', 90)
    write(unit_num,'(A32,12X,F10.3,24X,F12.2)') &
      'TOTAL', total_program_time, 100.0_dp
    write(unit_num,'(A)') repeat('=', 90)
    write(unit_num,*)

  end subroutine timer_report

  ! ====================================================================
  ! Reset all timers
  ! ====================================================================
  subroutine timer_reset(name)
    character(len=*), intent(in), optional :: name
    integer :: idx

    if (present(name)) then
      ! Reset specific timer
      idx = find_timer(name)
      if (idx > 0) then
        timers(idx)%total_time = 0.0_dp
        timers(idx)%call_count = 0
        timers(idx)%is_running = .false.
      endif
    else
      ! Reset all timers
      timers(:)%total_time = 0.0_dp
      timers(:)%call_count = 0
      timers(:)%is_running = .false.
    endif

  end subroutine timer_reset

  ! ====================================================================
  ! Cleanup and deallocate
  ! ====================================================================
  subroutine timer_cleanup()
    n_timers = 0
    timers(:)%name = ""
    timers(:)%total_time = 0.0_dp
    timers(:)%call_count = 0
    timers(:)%is_running = .false.
  end subroutine timer_cleanup

  ! ====================================================================
  ! Enable/disable timing
  ! ====================================================================
  subroutine timer_enable(enabled)
    logical, intent(in) :: enabled
    timer_enabled = enabled
  end subroutine timer_enable

  ! ====================================================================
  ! Get timing statistics for a specific timer
  ! ====================================================================
  subroutine timer_get_stats(name, total_time, call_count, avg_time, exists)
    character(len=*), intent(in) :: name
    real(dp), intent(out), optional :: total_time, avg_time
    integer, intent(out), optional :: call_count
    logical, intent(out), optional :: exists
    integer :: idx

    idx = find_timer(name)
    
    if (present(exists)) exists = (idx > 0)
    
    if (idx > 0) then
      if (present(total_time)) total_time = timers(idx)%total_time
      if (present(call_count)) call_count = timers(idx)%call_count
      if (present(avg_time)) then
        if (timers(idx)%call_count > 0) then
          avg_time = timers(idx)%total_time / real(timers(idx)%call_count, dp)
        else
          avg_time = 0.0_dp
        endif
      endif
    else
      if (present(total_time)) total_time = 0.0_dp
      if (present(call_count)) call_count = 0
      if (present(avg_time)) avg_time = 0.0_dp
    endif

  end subroutine timer_get_stats

  ! ====================================================================
  ! PRIVATE HELPER FUNCTIONS
  ! ====================================================================

  ! Find timer by name
  function find_timer(name) result(idx)
    character(len=*), intent(in) :: name
    integer :: idx, i

    idx = 0
    do i = 1, n_timers
      if (trim(timers(i)%name) == trim(name)) then
        idx = i
        return
      endif
    enddo

  end function find_timer

  ! Create new timer
  function create_timer(name) result(idx)
    character(len=*), intent(in) :: name
    integer :: idx

    if (n_timers >= MAX_TIMERS) then
      idx = 0
      return
    endif

    n_timers = n_timers + 1
    idx = n_timers
    timers(idx)%name = trim(name)
    timers(idx)%total_time = 0.0_dp
    timers(idx)%call_count = 0
    timers(idx)%is_running = .false.

  end function create_timer

  ! Sort timers
  subroutine sort_timers(indices, sort_by)
    integer, intent(out) :: indices(MAX_TIMERS)
    character(len=*), intent(in) :: sort_by
    integer :: i, j, temp
    real(dp) :: val_i, val_j

    ! Initialize indices
    do i = 1, n_timers
      indices(i) = i
    enddo

    ! Simple bubble sort (sufficient for small number of timers)
    do i = 1, n_timers - 1
      do j = i + 1, n_timers
        
        ! Get comparison values based on sort method
        select case (trim(sort_by))
        case ('time')
          val_i = timers(indices(i))%total_time
          val_j = timers(indices(j))%total_time
        case ('calls')
          val_i = real(timers(indices(i))%call_count, dp)
          val_j = real(timers(indices(j))%call_count, dp)
        case ('average')
          if (timers(indices(i))%call_count > 0) then
            val_i = timers(indices(i))%total_time / real(timers(indices(i))%call_count, dp)
          else
            val_i = 0.0_dp
          endif
          if (timers(indices(j))%call_count > 0) then
            val_j = timers(indices(j))%total_time / real(timers(indices(j))%call_count, dp)
          else
            val_j = 0.0_dp
          endif
        case ('name')
          ! Alphabetical sort
          if (timers(indices(i))%name > timers(indices(j))%name) then
            temp = indices(i)
            indices(i) = indices(j)
            indices(j) = temp
          endif
          cycle
        case default
          val_i = timers(indices(i))%total_time
          val_j = timers(indices(j))%total_time
        end select

        ! Swap if needed (descending order for numerical values)
        if (val_i < val_j) then
          temp = indices(i)
          indices(i) = indices(j)
          indices(j) = temp
        endif

      enddo
    enddo

  end subroutine sort_timers

end module module_timer
