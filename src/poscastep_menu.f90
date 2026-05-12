module poscastep_menu
    !! Interactive CLI menus for PosCASTEP post-processing
    !! Structure: top-level menu -> property-specific sub-menus
    !! Currently implements: Plot Band Structure
    use castep_config, only: dp, HARTREE_TO_EV, bands_data_t, pdos_data_t, &
        MAX_LINE_LEN, IO_INVALID_INPUT, IO_SUCCESS, strip_quotes
    use bands_parser, only: parse_bands_file, free_bands_data
    use bands_plotter, only: BANDS_MODE_ASCII, BANDS_MODE_SVG, &
        plot_bands_ascii, write_bands_svg
    use term_utils, only: get_term_size
    use pdos_parser, only: parse_pdos_file, free_pdos_data
    use dos_compute, only: compute_total_dos, compute_pdos, N_CHANNELS
    use dos_plotter, only: DOS_MODE_ASCII, DOS_MODE_SVG, DOS_MODE_EXPORT, &
        plot_dos_ascii, write_dos_svg, write_dos_csv, plot_pdos_ascii, write_pdos_csv
    implicit none
    private

    public :: run_poscastep_menu

    integer, parameter :: POS_BANDS = 1
    integer, parameter :: POS_DOS   = 2
    integer, parameter :: POS_PDOS  = 3

contains

    subroutine run_poscastep_menu(iostat)
        !! Top-level PosCASTEP post-processing menu.
        !! Displays available analyses, dispatches to sub-menus.
        !! Returns iostat = IO_SUCCESS on Q, non-zero on error.
        integer, intent(out) :: iostat
        integer :: choice, ios
        character(len=MAX_LINE_LEN) :: input

        iostat = IO_SUCCESS

        do
            write(*, '(a)') ''
            write(*, '(a)') '  ================================'
            write(*, '(a)') '            PosCASTEP'
            write(*, '(a)') '  ================================'
            write(*, '(a)') '  1. Plot Band Structure'
            write(*, '(a)') '  2. Plot DOS'
            write(*, '(a)') '  3. Plot pDOS'
            write(*, '(a)') '  Q. Back'
            write(*, '(a)', advance='no') '  Select option: '

            read(*, '(a)', iostat=ios) input
            if (ios /= 0) then
                iostat = IO_INVALID_INPUT
                write(*, '(a)') '  Invalid input.'
                return
            end if

            if (len_trim(input) >= 1) then
                if (input(1:1) == 'q' .or. input(1:1) == 'Q') then
                    iostat = IO_SUCCESS
                    return
                end if
            end if

            read(input, '(I6)', iostat=ios) choice
            if (ios /= 0) then
                write(*, '(a)') '  Invalid input. Enter a number.'
                cycle
            end if

            select case (choice)
            case (POS_BANDS)
                call handle_bands_menu(iostat)
                if (iostat /= 0) return
            case (POS_DOS)
                call handle_dos_menu(iostat)
                if (iostat /= 0) return
            case (POS_PDOS)
                call handle_pdos_menu(iostat)
                if (iostat /= 0) return
            case default
                write(*, '(a)') '  Invalid option. Enter 1, 2, 3, or Q.'
            end select
        end do
    end subroutine run_poscastep_menu


    subroutine handle_bands_menu(iostat)
        integer, intent(out) :: iostat
        type(bands_data_t) :: bands
        character(len=512)  :: bands_path, output_base, output_file
        character(len=256)  :: msg
        integer  :: plot_mode, local_istat
        real(dp) :: fermi_ev

        iostat = 0

        call ask_bands_path('Enter .bands file path: ', bands_path, iostat)
        if (iostat /= 0) return

        call ask_bands_plot_options(plot_mode, output_base, iostat)
        if (iostat /= 0) return

        write(*, '(a)') '  Parsing .bands file...'
        call parse_bands_file(trim(bands_path), bands, iostat, iomsg=msg)
        if (iostat /= 0) then
            write(*, '(a,a)') '  Error: ', trim(msg)
            return
        end if

        fermi_ev = bands%fermi_energy * HARTREE_TO_EV
        write(*, '(a)')        '  ------- Band Structure Summary -------'
        write(*, '(a,i0)')     '  K-points:         ', bands%num_kpoints
        write(*, '(a,i0)')     '  Eigenvalues:      ', bands%num_eigenvalues
        write(*, '(a,i0)')     '  Spin components:  ', bands%num_spin
        write(*, '(a,f10.4)')  '  Num electrons:    ', bands%num_electrons
        write(*, '(a,f10.4)')  '  k-path length:    ', bands%kpath_dist(bands%num_kpoints)
        write(*, '(a,f10.4,a)') '  Fermi energy:     ', fermi_ev, ' eV'

        select case (plot_mode)
        case (BANDS_MODE_ASCII)
            call run_ascii_navigator(bands, fermi_ev)

        case (BANDS_MODE_SVG)
            output_file = trim(output_base) // '.svg'
            call write_bands_svg(bands, output_file, &
                1200, 800, local_istat, iomsg=msg)
            if (local_istat == 0) then
                write(*, '(a,a)') '  SVG saved: ', trim(output_file)
            else
                write(*, '(a,a)') '  Error: ', trim(msg)
            end if
        end select

        call free_bands_data(bands)
        write(*, '(a)') '  ----------------------------------------'
    end subroutine handle_bands_menu


    subroutine run_ascii_navigator(bands, fermi_ev)
        type(bands_data_t), intent(in) :: bands
        real(dp), intent(in) :: fermi_ev
        real(dp) :: e_center, half_range, scroll_step, k_pct, k_width_pct
        character(len=1) :: ch
        integer :: ios, term_w, term_h

        half_range = 10.0_dp
        scroll_step = max(0.5_dp, half_range * 0.5_dp)
        e_center = fermi_ev
        k_pct = 0.5_dp
        k_width_pct = 1.0_dp

        call execute_command_line('stty -icanon -echo min 1', wait=.true.)
        do
            ! re-detect terminal size (handles resize events)
            call get_term_size(term_w, term_h)
            ! clear screen and redraw in-place
            write(*, '(a)', advance='no') achar(27) // '[2J' // achar(27) // '[H'
            call plot_bands_ascii(bands, term_w, term_h, e_center, half_range, k_pct, k_width_pct)
            if (k_width_pct > 0.99_dp) then
                write(*, '(a)') '  [↑↓ scroll  +/- zoom  ←→ scroll k-path  R reset  Q quit]'
            else
                write(*, '(a)') '  [↑↓ scroll  +/- zoom  ←→ scroll k-path  R reset  Q quit]'
            end if

            read(5, '(a1)', advance='no', iostat=ios) ch
            if (ios /= 0) exit
            if (ichar(ch) == 27) then  ! escape sequence
                read(5, '(a1)', advance='no', iostat=ios) ch
                if (ios /= 0) exit
                read(5, '(a1)', advance='no', iostat=ios) ch
                if (ios /= 0) exit
                if (ch == 'A') then
                    e_center = e_center + scroll_step      ! ↑
                else if (ch == 'B') then
                    e_center = e_center - scroll_step      ! ↓
                else if (ch == 'C') then                    ! →
                    if (k_width_pct > 0.99_dp) k_width_pct = 0.5_dp
                    k_pct = min(1.0_dp - k_width_pct/2.0_dp, k_pct + 0.1_dp * k_width_pct)
                else if (ch == 'D') then                    ! ←
                    if (k_width_pct > 0.99_dp) k_width_pct = 0.5_dp
                    k_pct = max(k_width_pct/2.0_dp, k_pct - 0.1_dp * k_width_pct)
                end if
            else if (ch == '+' .or. ch == '=') then
                half_range = max(0.25_dp, half_range * 0.5_dp)
                scroll_step = max(0.25_dp, half_range * 0.5_dp)
            else if (ch == '-') then
                half_range = min(20.0_dp, half_range * 2.0_dp)
                scroll_step = max(0.25_dp, half_range * 0.5_dp)
            else if (ch == 'r' .or. ch == 'R') then
                e_center = fermi_ev
                k_pct = 0.5_dp
                k_width_pct = 1.0_dp
            else if (ch == 'q' .or. ch == 'Q') then
                exit
            end if
        end do

        write(*, '(a)', advance='no') achar(27) // '[2J' // achar(27) // '[H'
        call execute_command_line('stty sane', wait=.true.)
        write(*, '(a)') ''
    end subroutine run_ascii_navigator


    subroutine ask_bands_path(prompt_text, result_path, iostat)
        character(len=*), intent(in)  :: prompt_text
        character(len=*), intent(out) :: result_path
        integer, intent(out)          :: iostat
        integer :: ios
        logical :: exists
        character(len=MAX_LINE_LEN) :: input

        iostat = 0
        do
            write(*, '(a)', advance='no') trim(prompt_text)
            read(*, '(a)', iostat=ios) input
            if (ios /= 0) then
                iostat = IO_INVALID_INPUT
                write(*, '(a)') '  Error reading input.'
                return
            end if
            result_path = adjustl(trim(input))
            call strip_quotes(result_path)
            if (len_trim(result_path) == 0) then
                write(*, '(a)') '  Path cannot be empty. Try again.'
                cycle
            end if
            inquire(file=trim(result_path), exist=exists)
            if (.not. exists) then
                write(*, '(a)') '  File not found: ' // trim(result_path)
                write(*, '(a)') '  Try again or press Ctrl+C to abort.'
                cycle
            end if
            exit
        end do
    end subroutine ask_bands_path


    subroutine ask_bands_plot_options(plot_mode, output_base, iostat)
        integer, intent(out)          :: plot_mode
        character(len=*), intent(out) :: output_base
        integer, intent(out)          :: iostat
        integer :: ios, choice
        character(len=MAX_LINE_LEN) :: input

        iostat = 0
        plot_mode   = BANDS_MODE_ASCII
        output_base = 'bands'

        write(*, '(a)') ''
        write(*, '(a)') '  Output mode:'
        write(*, '(a)') '    1. Terminal ASCII plot (default)'
        write(*, '(a)') '    2. SVG vector graphics file'
        write(*, '(a)', advance='no') '    Enter choice [1]: '
        read(*, '(a)', iostat=ios) input
        if (ios == 0 .and. len_trim(input) > 0) then
            read(input, '(I6)', iostat=ios) choice
            if (ios == 0 .and. choice == 2) plot_mode = BANDS_MODE_SVG
        end if

        if (plot_mode == BANDS_MODE_SVG) then
            write(*, '(a)') ''
            write(*, '(a)', advance='no') '  Output base name [bands]: '
            read(*, '(a)', iostat=ios) input
            if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
                output_base = adjustl(trim(input))
            end if
        end if
    end subroutine ask_bands_plot_options

    ! ----------------------------------------------------------------
    !  DOS handler
    ! ----------------------------------------------------------------
    subroutine handle_dos_menu(iostat)
        integer, intent(out) :: iostat
        type(bands_data_t) :: bands
        character(len=512), save :: last_bands_path = ''
        character(len=512)  :: bands_path, output_file
        character(len=256)  :: msg
        integer  :: plot_mode, local_istat, n_spin, ne, ios
        real(dp) :: fermi_ev, smearing_width
        real(dp), allocatable :: energy_grid(:), dos_result(:,:)
        character(len=MAX_LINE_LEN) :: input

        iostat = 0
        smearing_width = 0.1_dp

        ! --- .bands path ---
        if (len_trim(last_bands_path) > 0) then
            write(*, '(a,a,a)', advance='no') '  Enter .bands file path [', &
                trim(last_bands_path), ']: '
        else
            write(*, '(a)', advance='no') '  Enter .bands file path: '
        end if
        read(*, '(a)', iostat=ios) input
        bands_path = adjustl(trim(input))
        call strip_quotes(bands_path)
        if (len_trim(bands_path) == 0) then
            if (len_trim(last_bands_path) > 0) then
                bands_path = last_bands_path
                write(*, '(a,a)') '  Using: ', trim(bands_path)
            else
                write(*, '(a)') '  No path provided. Aborted.'
                iostat = 1; return
            end if
        else
            block
                logical :: exists
                inquire(file=trim(bands_path), exist=exists)
                if (.not. exists) then
                    write(*, '(a)') '  File not found: ' // trim(bands_path)
                    iostat = 1; return
                end if
            end block
            last_bands_path = bands_path
        end if

        write(*, '(a)') ''
        write(*, '(a)') '  Output mode:'
        write(*, '(a)') '    1. Terminal ASCII plot (default)'
        write(*, '(a)') '    2. SVG vector graphics file'
        write(*, '(a)') '    3. Export data (CSV)'
        write(*, '(a)', advance='no') '    Enter choice [1]: '
        read(*, '(a)', iostat=ios) input
        plot_mode = DOS_MODE_ASCII
        if (ios == 0 .and. len_trim(input) > 0) then
            read(input, *, iostat=ios) plot_mode
            if (ios /= 0) plot_mode = DOS_MODE_ASCII
        end if

        write(*, '(a,f5.2,a)', advance='no') '  Smearing width [', smearing_width, ' eV]: '
        read(*, '(a)', iostat=ios) input
        if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
            read(input, *, iostat=ios) smearing_width
        end if
        smearing_width = max(0.01_dp, min(1.0_dp, smearing_width))

        write(*, '(a)') '  Parsing .bands file...'
        call parse_bands_file(trim(bands_path), bands, iostat, iomsg=msg)
        if (iostat /= 0) then
            write(*, '(a,a)') '  Error: ', trim(msg)
            return
        end if

        fermi_ev = bands%fermi_energy * HARTREE_TO_EV
        n_spin  = bands%num_spin

        write(*, '(a)')        '  ------- DOS Summary -------'
        write(*, '(a,i0)')     '  K-points:         ', bands%num_kpoints
        write(*, '(a,i0)')     '  Bands:            ', bands%num_eigenvalues
        write(*, '(a,i0)')     '  Spin components:  ', n_spin
        write(*, '(a,f10.4)')  '  Fermi energy:     ', fermi_ev, ' eV'
        write(*, '(a,f5.2,a)') '  Smearing:         ', smearing_width, ' eV'

        ne = 4001
        allocate(energy_grid(ne))
        call build_energy_grid(energy_grid, -20.0_dp, 20.0_dp, ne)

        write(*, '(a)') '  Computing total DOS...'
        call compute_total_dos(bands, energy_grid, smearing_width, &
            dos_result, iostat, iomsg=msg)
        if (iostat /= 0) then
            write(*, '(a,a)') '  Error: ', trim(msg)
            deallocate(energy_grid)
            call free_bands_data(bands)
            return
        end if

        select case (plot_mode)
        case (DOS_MODE_ASCII)
            call run_dos_navigator(energy_grid, dos_result, n_spin, &
                fermi_ev, smearing_width)
        case (DOS_MODE_SVG)
            output_file = 'dos.svg'
            write(*, '(a)', advance='no') '  Output file [dos.svg]: '
            read(*, '(a)', iostat=ios) input
            if (ios == 0 .and. len_trim(adjustl(input)) > 0) &
                output_file = adjustl(trim(input))
            call write_dos_svg(energy_grid, dos_result, n_spin, &
                fermi_ev, smearing_width, output_file, 1200, 800, &
                local_istat, iomsg=msg)
            if (local_istat == 0) then
                write(*, '(a,a)') '  SVG saved: ', trim(output_file)
            else
                write(*, '(a,a)') '  Error: ', trim(msg)
            end if
        case (DOS_MODE_EXPORT)
            output_file = 'dos'
            write(*, '(a)', advance='no') '  Output file [dos]: '
            read(*, '(a)', iostat=ios) input
            if (ios == 0 .and. len_trim(adjustl(input)) > 0) &
                output_file = adjustl(trim(input))
            call ensure_ext(output_file, '.csv')
            call write_dos_csv(energy_grid, dos_result, n_spin, &
                output_file, local_istat, iomsg=msg)
            if (local_istat == 0) then
                write(*, '(a,a)') '  Data exported: ', trim(output_file)
            else
                write(*, '(a,a)') '  Error: ', trim(msg)
            end if
        end select

        deallocate(energy_grid, dos_result)
        call free_bands_data(bands)
        write(*, '(a)') '  ------------------------------'
    end subroutine handle_dos_menu

    ! ----------------------------------------------------------------
    !  pDOS handler (prefix-based: <prefix>.bands + <prefix>.pdos_bin)
    ! ----------------------------------------------------------------
    subroutine handle_pdos_menu(iostat)
        integer, intent(out) :: iostat
        type(bands_data_t) :: bands
        type(pdos_data_t)  :: pdos
        character(len=512), save :: last_prefix = ''
        character(len=512)  :: prefix, bands_path, pdos_path, output_file
        character(len=256)  :: msg
        integer  :: plot_mode, local_istat, ne, ios
        real(dp) :: fermi_ev, smearing_width
        real(dp), allocatable :: energy_grid(:), pdos_result(:,:,:)
        character(len=MAX_LINE_LEN) :: input

        iostat = 0
        smearing_width = 0.1_dp

        ! --- prefix prompt ---
        if (len_trim(last_prefix) > 0) then
            write(*, '(a,a,a)', advance='no') '  Enter file prefix (e.g. Cu) [', &
                trim(last_prefix), ']: '
        else
            write(*, '(a)', advance='no') '  Enter file prefix (e.g. Cu): '
        end if
        read(*, '(a)', iostat=ios) input
        prefix = adjustl(trim(input))
        call strip_quotes(prefix)
        if (len_trim(prefix) == 0) then
            if (len_trim(last_prefix) > 0) then
                prefix = last_prefix
                write(*, '(a,a)') '  Using prefix: ', trim(prefix)
            else
                write(*, '(a)') '  No prefix provided. Aborted.'
                iostat = 1; return
            end if
        else
            last_prefix = prefix
        end if

        bands_path = trim(prefix) // '.bands'
        pdos_path  = trim(prefix) // '.pdos_bin'

        ! verify .bands
        block
            logical :: ex_b, ex_p
            inquire(file=trim(bands_path), exist=ex_b)
            inquire(file=trim(pdos_path), exist=ex_p)
            if (.not. ex_b) then
                ! try .pdos_weights
                if (.not. ex_p) then
                    pdos_path = trim(prefix) // '.pdos_weights'
                    inquire(file=trim(pdos_path), exist=ex_p)
                end if
                if (.not. ex_b) then
                    write(*, '(a)') '  File not found: ' // trim(bands_path)
                    iostat = 1; return
                end if
            end if
            if (.not. ex_p) then
                write(*, '(a)') '  File not found: ' // trim(pdos_path)
                write(*, '(a)') '  (need either .pdos_bin or .pdos_weights)'
                iostat = 1; return
            end if
        end block

        write(*, '(a,a)') '  .bands: ', trim(bands_path)
        write(*, '(a,a)') '  .pdos:  ', trim(pdos_path)

        write(*, '(a)') ''
        write(*, '(a)') '  Output mode:'
        write(*, '(a)') '    1. Terminal ASCII plot (default)'
        write(*, '(a)') '    2. Export data (CSV)'
        write(*, '(a)', advance='no') '    Enter choice [1]: '
        read(*, '(a)', iostat=ios) input
        plot_mode = DOS_MODE_ASCII
        if (ios == 0 .and. len_trim(input) > 0) then
            read(input, *, iostat=ios) plot_mode
            if (ios /= 0) plot_mode = DOS_MODE_ASCII
        end if

        write(*, '(a,f5.2,a)', advance='no') '  Smearing width [', smearing_width, ' eV]: '
        read(*, '(a)', iostat=ios) input
        if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
            read(input, *, iostat=ios) smearing_width
        end if
        smearing_width = max(0.01_dp, min(1.0_dp, smearing_width))

        write(*, '(a)') '  Parsing .bands file...'
        call parse_bands_file(trim(bands_path), bands, iostat, iomsg=msg)
        if (iostat /= 0) then
            write(*, '(a,a)') '  Error: ', trim(msg)
            return
        end if

        fermi_ev = bands%fermi_energy * HARTREE_TO_EV

        write(*, '(a)') '  Parsing PDOS weights file...'
        call parse_pdos_file(trim(pdos_path), pdos, local_istat, iomsg=msg)
        if (local_istat /= 0) then
            write(*, '(a,a)') '  Error: ', trim(msg)
            call free_bands_data(bands)
            return
        end if

        write(*, '(a)')        '  ------- pDOS Summary -------'
        write(*, '(a,i0)')     '  K-points:         ', bands%num_kpoints
        write(*, '(a,i0)')     '  Bands:            ', bands%num_eigenvalues
        write(*, '(a,i0)')     '  Spin components:  ', bands%num_spin
        write(*, '(a,f10.4)')  '  Fermi energy:     ', fermi_ev, ' eV'
        write(*, '(a,f5.2,a)') '  Smearing:         ', smearing_width, ' eV'
        write(*, '(a,i0)')     '  PDOS orbitals:    ', pdos%num_orbitals

        ne = 4001
        allocate(energy_grid(ne))
        call build_energy_grid(energy_grid, -20.0_dp, 20.0_dp, ne)

        write(*, '(a)') '  Computing pDOS (s/p/d/f)...'
        call compute_pdos(bands, pdos, energy_grid, smearing_width, &
            pdos_result, local_istat, iomsg=msg)
        if (local_istat /= 0) then
            write(*, '(a,a)') '  Error: ', trim(msg)
            deallocate(energy_grid)
            call free_pdos_data(pdos)
            call free_bands_data(bands)
            return
        end if

        select case (plot_mode)
        case (DOS_MODE_ASCII)
            call run_pdos_navigator(energy_grid, pdos_result(:,:,1), &
                fermi_ev, smearing_width)
        case (DOS_MODE_EXPORT)
            output_file = prefix
            write(*, '(a,a,a)', advance='no') '  Output file [', trim(output_file), '.csv]: '
            read(*, '(a)', iostat=ios) input
            if (ios == 0 .and. len_trim(adjustl(input)) > 0) &
                output_file = adjustl(trim(input))
            call ensure_ext(output_file, '.csv')
            call write_pdos_csv(energy_grid, pdos_result(:,:,1), &
                output_file, local_istat, iomsg=msg)
            if (local_istat == 0) then
                write(*, '(a,a)') '  Data exported: ', trim(output_file)
            else
                write(*, '(a,a)') '  Error: ', trim(msg)
            end if
        end select

        deallocate(energy_grid, pdos_result)
        call free_pdos_data(pdos)
        call free_bands_data(bands)
        write(*, '(a)') '  ------------------------------'
    end subroutine handle_pdos_menu

    subroutine run_dos_navigator(energy_grid, dos_data, nspin, e_fermi, smearing)
        real(dp), intent(in) :: energy_grid(:), dos_data(:,:), e_fermi, smearing
        integer, intent(in) :: nspin
        real(dp) :: e_center, e_half, y_center, y_half, y_half0, x_step, y_step
        character(len=1) :: ch
        integer :: ios, term_w, term_h, is, ie

        e_half  = 10.0_dp
        e_center = 0.0_dp
        ! initial y-axis: [0, y_max0]
        y_half = 1.0_dp
        do is = 1, nspin
            do ie = 1, size(dos_data, 1)
                if (dos_data(ie, is) > y_half) y_half = dos_data(ie, is)
            end do
        end do
        y_half  = y_half * 1.15_dp * 0.5_dp
        if (y_half < 1.0e-12_dp) y_half = 1.0_dp
        y_half0   = y_half
        y_center  = y_half  ! range: [0, 2*y_half]

        call execute_command_line('stty -icanon -echo min 1', wait=.true.)
        do
            x_step = e_half * 0.25_dp
            y_step = y_half * 0.4_dp
            call get_term_size(term_w, term_h)
            write(*, '(a)', advance='no') achar(27) // '[2J' // achar(27) // '[H'
            call plot_dos_ascii(energy_grid, dos_data, nspin, e_fermi, smearing, &
                term_w, term_h, y_center_in=y_center, y_half_in=y_half, &
                e_center_in=e_center, half_range_in=e_half)
            write(*, '(a)') '  [↑↓ y-pan  ←→ x-pan  +/- zoom  R reset  Q quit]'

            read(5, '(a1)', advance='no', iostat=ios) ch
            if (ios /= 0) exit
            if (ichar(ch) == 27) then
                read(5, '(a1)', advance='no', iostat=ios) ch
                if (ios /= 0) exit
                read(5, '(a1)', advance='no', iostat=ios) ch
                if (ios /= 0) exit
                if (ch == 'A') then           ! ↑  vertical pan up
                    y_center = y_center + y_step
                else if (ch == 'B') then      ! ↓  vertical pan down
                    y_center = y_center - y_step
                else if (ch == 'C') then      ! →  horizontal pan right
                    e_center = e_center + x_step
                else if (ch == 'D') then      ! ←  horizontal pan left
                    e_center = e_center - x_step
                end if
            else if (ch == '+' .or. ch == '=') then
                e_half  = max(0.25_dp, e_half * 0.5_dp)
                y_half  = max(0.01_dp, y_half * 0.5_dp)
            else if (ch == '-') then
                e_half  = min(20.0_dp, e_half * 2.0_dp)
                y_half  = y_half * 2.0_dp
            else if (ch == 'r' .or. ch == 'R') then
                e_center = 0.0_dp
                e_half   = 10.0_dp
                y_center = y_half0
                y_half   = y_half0
            else if (ch == 'q' .or. ch == 'Q') then
                exit
            end if
        end do

        write(*, '(a)', advance='no') achar(27) // '[2J' // achar(27) // '[H'
        call execute_command_line('stty sane', wait=.true.)
        write(*, '(a)') ''
    end subroutine run_dos_navigator

    subroutine run_pdos_navigator(energy_grid, pdos_data, e_fermi, smearing)
        real(dp), intent(in) :: energy_grid(:), pdos_data(:,:), e_fermi, smearing
        real(dp) :: e_center, e_half, y_center, y_half, y_half0, x_step, y_step
        character(len=1) :: ch
        integer :: ios, term_w, term_h, ich, ie

        e_half  = 10.0_dp
        e_center = 0.0_dp
        y_half = 1.0_dp
        do ich = 2, size(pdos_data, 2)
            do ie = 1, size(pdos_data, 1)
                if (pdos_data(ie, ich) > y_half) y_half = pdos_data(ie, ich)
            end do
        end do
        y_half  = y_half * 1.15_dp * 0.5_dp
        if (y_half < 1.0e-12_dp) y_half = 1.0_dp
        y_half0   = y_half
        y_center  = y_half  ! range: [0, 2*y_half]

        call execute_command_line('stty -icanon -echo min 1', wait=.true.)
        do
            x_step = e_half * 0.25_dp
            y_step = y_half * 0.4_dp
            call get_term_size(term_w, term_h)
            write(*, '(a)', advance='no') achar(27) // '[2J' // achar(27) // '[H'
            call plot_pdos_ascii(energy_grid, pdos_data, e_fermi, smearing, &
                term_w, term_h, y_center_in=y_center, y_half_in=y_half, &
                e_center_in=e_center, half_range_in=e_half)
            write(*, '(a)') '  [↑↓ y-pan  ←→ x-pan  +/- zoom  R reset  Q quit]'

            read(5, '(a1)', advance='no', iostat=ios) ch
            if (ios /= 0) exit
            if (ichar(ch) == 27) then
                read(5, '(a1)', advance='no', iostat=ios) ch
                if (ios /= 0) exit
                read(5, '(a1)', advance='no', iostat=ios) ch
                if (ios /= 0) exit
                if (ch == 'A') then           ! ↑
                    y_center = y_center + y_step
                else if (ch == 'B') then      ! ↓
                    y_center = y_center - y_step
                else if (ch == 'C') then      ! →
                    e_center = e_center + x_step
                else if (ch == 'D') then      ! ←
                    e_center = e_center - x_step
                end if
            else if (ch == '+' .or. ch == '=') then
                e_half  = max(0.25_dp, e_half * 0.5_dp)
                y_half  = max(0.01_dp, y_half * 0.5_dp)
            else if (ch == '-') then
                e_half  = min(20.0_dp, e_half * 2.0_dp)
                y_half  = y_half * 2.0_dp
            else if (ch == 'r' .or. ch == 'R') then
                e_center = 0.0_dp
                e_half   = 10.0_dp
                y_center = y_half0
                y_half   = y_half0
            else if (ch == 'q' .or. ch == 'Q') then
                exit
            end if
        end do

        write(*, '(a)', advance='no') achar(27) // '[2J' // achar(27) // '[H'
        call execute_command_line('stty sane', wait=.true.)
        write(*, '(a)') ''
    end subroutine run_pdos_navigator

    subroutine build_energy_grid(grid, e_min, e_max, npts)
        real(dp), intent(out) :: grid(:)
        real(dp), intent(in)  :: e_min, e_max
        integer, intent(in)   :: npts
        integer :: i
        do i = 1, npts
            grid(i) = e_min + real(i-1, dp) / real(npts-1, dp) * (e_max - e_min)
        end do
    end subroutine build_energy_grid

    subroutine ensure_ext(filename, ext)
        character(len=*), intent(inout) :: filename
        character(len=*), intent(in)    :: ext
        integer :: nf, ne
        nf = len_trim(filename)
        ne = len_trim(ext)
        if (nf < ne) then
            filename = trim(filename) // ext
        else if (filename(nf-ne+1:nf) /= ext) then
            filename = trim(filename) // ext
        end if
    end subroutine ensure_ext

end module poscastep_menu
