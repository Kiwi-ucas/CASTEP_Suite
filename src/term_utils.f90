module term_utils
    !! Terminal utilities shared by plotting modules
    !! ANSI color codes, terminal size detection, Bresenham line drawing
    !! Zero dependencies — pure leaf module
    implicit none
    private

    ! ANSI color codes
    character(len=*), parameter, public :: C_RED   = achar(27) // '[1;31m'
    character(len=*), parameter, public :: C_GREEN = achar(27) // '[1;32m'
    character(len=*), parameter, public :: C_YELLOW= achar(27) // '[1;33m'
    character(len=*), parameter, public :: C_CYAN  = achar(27) // '[36m'
    character(len=*), parameter, public :: C_BOLD  = achar(27) // '[1m'
    character(len=*), parameter, public :: C_DIM   = achar(27) // '[2m'
    character(len=*), parameter, public :: C_RESET = achar(27) // '[0m'
    character(len=*), parameter, public :: C_AXIS  = achar(27) // '[1;37m'

    ! Alternate screen buffer (prevents scrollback pollution)
    character(len=*), parameter, public :: C_ALT_ON  = achar(27) // '[?1049h'
    character(len=*), parameter, public :: C_ALT_OFF = achar(27) // '[?1049l'

    public :: get_term_size, draw_line, enter_raw_mode, leave_raw_mode

contains

    subroutine get_term_size(tw, th)
        integer, intent(out) :: tw, th
        character(len=16) :: str
        integer :: ios
        tw = 160; th = 40   ! defaults
        ! query terminal driver directly (works even during stty -icanon)
        call stty_size(tw, th)
        if (tw == 160 .and. th == 40) then
            ! fall back to env vars if stty failed
            call get_environment_variable('COLUMNS', str, status=ios)
            if (ios == 0 .and. len_trim(str) > 0) then
                read(str, *, iostat=ios) tw
                if (ios == 0 .and. tw < 40) tw = 160
            end if
            call get_environment_variable('LINES', str, status=ios)
            if (ios == 0 .and. len_trim(str) > 0) then
                read(str, *, iostat=ios) th
                if (ios == 0 .and. th < 15) th = 40
            end if
        end if
    end subroutine get_term_size

    subroutine stty_size(tw, th)
        integer, intent(out) :: tw, th
        character(len=32) :: line
        integer :: unit, ios, sp
        tw = 160; th = 40
        call execute_command_line('stty size > .stty_tmp 2>/dev/null', &
            wait=.true.)
        open(newunit=unit, file='.stty_tmp', status='old', iostat=ios)
        if (ios /= 0) return
        read(unit, '(a)', iostat=ios) line
        close(unit, status='delete')
        if (ios /= 0 .or. len_trim(line) == 0) return
        sp = index(trim(line), ' ')
        if (sp > 1) then
            read(line(1:sp-1), *, iostat=ios) th
            if (ios == 0) then
                read(line(sp+1:), *, iostat=ios) tw
                if (ios == 0 .and. tw >= 40) return
            end if
        end if
        tw = 160; th = 40
    end subroutine stty_size

    subroutine draw_line(nw, nh, grid, x1, y1, x2, y2)
        integer, intent(in) :: x1, y1, x2, y2, nw, nh
        character(len=1), intent(inout) :: grid(nw, nh)
        integer :: dx, dy, sx, sy, err, e2, x, y
        dx = abs(x2 - x1); dy = -abs(y2 - y1)
        sx = 1; if (x1 > x2) sx = -1
        sy = 1; if (y1 > y2) sy = -1
        err = dx + dy
        x = x1; y = y1
        do
            if (grid(x, y) == ' ') grid(x, y) = '.'
            if (x == x2 .and. y == y2) exit
            e2 = 2 * err
            if (e2 >= dy) then
                if (x == x2) exit
                err = err + dy; x = x + sx
            end if
            if (e2 <= dx) then
                if (y == y2) exit
                err = err + dx; y = y + sy
            end if
        end do
    end subroutine draw_line

    subroutine enter_raw_mode()
        call execute_command_line('stty -icanon -echo min 1', wait=.true.)
        write(*, '(a)', advance='no') C_ALT_ON
    end subroutine enter_raw_mode

    subroutine leave_raw_mode()
        write(*, '(a)', advance='no') C_ALT_OFF
        call execute_command_line('stty sane', wait=.true.)
    end subroutine leave_raw_mode

end module term_utils
