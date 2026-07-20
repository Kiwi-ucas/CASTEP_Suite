module poscastep_menu
    !! Interactive CLI menus for PosCASTEP post-processing
    !! Structure: top-level menu -> property-specific sub-menus
    !! Currently implements: Plot Band Structure
    use castep_config, only: dp, pi, HARTREE_TO_EV, bands_data_t, pdos_data_t, &
        cif_data_t, atom_t, castep_config_t, MAX_LINE_LEN, &
        IO_INVALID_INPUT, IO_SUCCESS, IO_USER_QUIT, &
        IO_FILE_NOT_FOUND, IO_PARSE_ERROR, IO_WRITE_ERROR, IO_PRECASTEP_LAUNCH, &
        TASK_ENERGY, TASK_GEOMETRY_OPT, PSEUDO_C19MK2, KPOINT_GAMMA, &
        default_config, strip_quotes, compute_cartesian_lattice
    use parser, only: parse_cif_inline, parse_pdb_inline, parse_cell_inline, &
        clean_element_symbol
    use phonon_dos, only: phonon_dos_data_t, parse_phonon_file, compute_phonon_dos, &
        free_phonon_dos_data
    use castep_vib, only: vib_data_t, parse_castep_vib, compute_ir_spectrum, &
        compute_raman_spectrum, free_vib_data
    use phonon_modes, only: phonon_modes_data_t, parse_phonon_eigenvectors, &
        parse_castep_born_charges, compute_mode_decomposition, free_phonon_modes_data
    use crystal_json, only: write_crystal_json_modes
    use polarizability, only: pol_data_t, parse_castep_file, &
        parse_cp2k_dipoles, unwrap_dipoles, &
        compute_static_dielectric_windowed, &
        compute_polarizability, free_pol_data
    use bands_parser, only: parse_bands_file, free_bands_data
    use bands_plotter, only: BANDS_MODE_ASCII, plot_bands_ascii
    use term_utils, only: get_term_size, enter_raw_mode, leave_raw_mode
    use pdos_parser, only: parse_pdos_file, free_pdos_data
    use dos_compute, only: compute_total_dos, compute_pdos, N_CHANNELS
    use dos_plotter, only: DOS_MODE_ASCII, DOS_MODE_EXPORT, &
        plot_dos_ascii, write_dos_csv, plot_pdos_ascii, write_pdos_csv
    use crystal_json, only: write_crystal_json_cif, read_crystal_json_to_cif
    use thermodynamics, only: thermo_data_t, compute_thermodynamics, free_thermo_data
    use symmetry, only: expand_cif_symmetry
    use cell_writer, only: write_cell_file
    use param_writer, only: write_param_file
    use cli_menu, only: run_main_menu
    use pes3d, only: pes3d_grid_t, compute_local_grid_bounds, generate_pes3d_grid_points, &
        write_pes3d_cube, collect_pes3d_energies, symmetry_expand_energies
    use pes_scan, only: pes_grid_t, generate_pes_grid_points, write_pes_metadata_json, &
        collect_pes_energies
    implicit none
    private

    public :: run_poscastep_menu, free_cif_data

    integer, parameter :: POS_CONVERTER  = 0
    integer, parameter :: POS_BANDS      = 1
    integer, parameter :: POS_DOS        = 2
    integer, parameter :: POS_PDOS       = 3
    integer, parameter :: POS_PHONON_DOS = 4
    integer, parameter :: POS_IR_SPEC    = 5
    integer, parameter :: POS_RAMAN_SPEC = 6
    integer, parameter :: POS_POLARIZABILITY = 7
    integer, parameter :: POS_PHONON_MODES  = -2
    integer, parameter :: POS_THERMO        = 8
    integer, parameter :: POS_PES_SCAN      = 9
    integer, parameter :: POS_VIEW_STRUCTURE = -1

    ! Module-level storage for PreCASTEP handoff (option 2: no file on disk)
    type(cif_data_t), save, public :: precastep_cif_data
    character(len=MAX_LINE_LEN), save, public :: precastep_source_file = ''
    character(len=MAX_LINE_LEN), save :: precastep_viewer_file = ''
    logical, save, public :: has_precastep_data = .false.

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
            write(*, '(a)') ' -2. Phonon Mode Visualization'
            write(*, '(a)') ' -1. View Crystal Structure (3D)'
            write(*, '(a)') '  0. Format Converter (.cell/.cif/.pdb)'
            write(*, '(a)') '  1. Plot Band Structure'
            write(*, '(a)') '  2. Plot DOS'
            write(*, '(a)') '  3. Plot pDOS'
            write(*, '(a)') '  4. Plot Phonon DOS'
            write(*, '(a)') '  5. Plot IR Spectrum'
            write(*, '(a)') '  6. Plot Raman Spectrum'
            write(*, '(a)') '  7. Static Polarizability'
            write(*, '(a)') '  8. Thermodynamics'
            write(*, '(a)') '  9. PES Scan'
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
            case (POS_CONVERTER)
                call handle_format_converter_menu(iostat)
                if (iostat == IO_USER_QUIT) return
            case (POS_BANDS)
                call handle_bands_menu(iostat)
                if (iostat == IO_USER_QUIT) return
            case (POS_DOS)
                call handle_dos_menu(iostat)
                if (iostat == IO_USER_QUIT) return
            case (POS_PDOS)
                call handle_pdos_menu(iostat)
                if (iostat == IO_USER_QUIT) return
            case (POS_PHONON_DOS)
                call handle_phonon_dos_menu(iostat)
                if (iostat == IO_USER_QUIT) return
            case (POS_IR_SPEC)
                call handle_ir_menu(iostat)
                if (iostat == IO_USER_QUIT) return
            case (POS_RAMAN_SPEC)
                call handle_raman_menu(iostat)
                if (iostat == IO_USER_QUIT) return
            case (POS_POLARIZABILITY)
                call handle_polarizability_menu(iostat)
                if (iostat == IO_USER_QUIT) return
            case (POS_PHONON_MODES)
                call handle_phonon_modes_menu(iostat)
                if (iostat == IO_USER_QUIT) return
            case (POS_THERMO)
                call handle_thermo_menu(iostat)
                if (iostat == IO_USER_QUIT) return
            case (POS_PES_SCAN)
                call handle_pes_scan_menu(iostat)
                if (iostat == IO_USER_QUIT) return
            case (POS_VIEW_STRUCTURE)
                call handle_view_structure(iostat)
                if (iostat == IO_USER_QUIT) return
                if (iostat == IO_PRECASTEP_LAUNCH) return
            case default
                write(*, '(a)') '  Invalid option. Enter 0-9, -1, or Q.'
            end select
        end do
    end subroutine run_poscastep_menu


    subroutine handle_bands_menu(iostat)
        integer, intent(out) :: iostat
        type(bands_data_t) :: bands
        character(len=512)  :: bands_path, output_base
        character(len=256)  :: msg
        integer  :: plot_mode
        real(dp) :: fermi_ev

        iostat = 0

        call ask_bands_path('Enter .bands file path: ', bands_path, iostat)
        if (iostat == IO_USER_QUIT) then
            iostat = 0; return
        end if
        if (iostat /= 0) return

        call ask_bands_plot_options(plot_mode, output_base, iostat)
        if (iostat == IO_USER_QUIT) then
            iostat = 0; return
        end if
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

        call enter_raw_mode
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
                k_width_pct = max(0.05_dp, k_width_pct * 0.5_dp)
            else if (ch == '-') then
                half_range = min(20.0_dp, half_range * 2.0_dp)
                scroll_step = max(0.25_dp, half_range * 0.5_dp)
                k_width_pct = min(1.0_dp, k_width_pct * 2.0_dp)
            else if (ch == 'r' .or. ch == 'R') then
                e_center = fermi_ev
                k_pct = 0.5_dp
                k_width_pct = 1.0_dp
            else if (ch == 'q' .or. ch == 'Q') then
                exit
            end if
        end do

        call leave_raw_mode
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
            if (result_path == 'q' .or. result_path == 'Q') then
                iostat = IO_USER_QUIT
                return
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
        integer :: ios
        character(len=MAX_LINE_LEN) :: input

        iostat = 0
        plot_mode   = BANDS_MODE_ASCII
        output_base = 'bands'

        write(*, '(a)') ''
        write(*, '(a)') '  Output mode:'
        write(*, '(a)') '    1. Terminal ASCII plot (default)'
        write(*, '(a)', advance='no') '    Enter choice [1]: '
        read(*, '(a)', iostat=ios) input
        if (ios == 0 .and. len_trim(input) > 0) then
            if (adjustl(trim(input)) == 'q' .or. adjustl(trim(input)) == 'Q') then
                iostat = IO_USER_QUIT; return
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
                iostat = 0; return
            end if
        else if (bands_path == 'q' .or. bands_path == 'Q') then
            iostat = 0; return
        else
            block
                logical :: exists
                inquire(file=trim(bands_path), exist=exists)
                if (.not. exists) then
                    write(*, '(a)') '  File not found: ' // trim(bands_path)
                    iostat = 0; return
                end if
            end block
            last_bands_path = bands_path
        end if

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
                iostat = 0; return
            end if
        else if (prefix == 'q' .or. prefix == 'Q') then
            iostat = 0; return
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
                    iostat = 0; return
                end if
            end if
            if (.not. ex_p) then
                write(*, '(a)') '  File not found: ' // trim(pdos_path)
                write(*, '(a)') '  (need either .pdos_bin or .pdos_weights)'
                iostat = 0; return
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

    subroutine handle_phonon_dos_menu(iostat)
        integer, intent(out) :: iostat
        type(phonon_dos_data_t) :: phdos
        character(len=MAX_LINE_LEN) :: fname, input, tmp_str
        integer :: ios, plot_mode
        real(dp) :: freq_min_range, freq_max_range, smearing_width
        character(len=MAX_LINE_LEN), save :: last_phonon_path = ''

        iostat = 0

        ! ── File input ──
        write(*, '(a)', advance='no') '  Enter .phonon file path'
        if (len_trim(last_phonon_path) > 0) &
            write(*, '(a)', advance='no') ' [' // trim(last_phonon_path) // ']'
        write(*, '(a)') ': '

        read(*, '(a)', iostat=ios) fname
        if (ios /= 0) return
        fname = adjustl(fname); call strip_quotes(fname)
        if (fname == 'q' .or. fname == 'Q') then
            iostat = 0; return
        end if
        if (len_trim(fname) == 0 .and. len_trim(last_phonon_path) > 0) then
            fname = last_phonon_path
        end if
        if (len_trim(fname) == 0) then
            write(*, '(a)') '  No file specified.'
            return
        end if
        last_phonon_path = trim(fname)

        call parse_phonon_file(trim(fname), phdos, ios)
        if (ios /= 0) then
            write(*, '(a)') '  Error parsing .phonon file.'
            return
        end if

        write(*, '(a,i0,a,i0,a,i0,a)') '  Loaded ', phdos%n_ions, ' ions, ', &
            phdos%n_branches, ' branches, ', phdos%n_qpoints, ' q-points.'

        ! ── Output mode menu (before smearing, matching DOS flow) ──
        write(*, '(a)') '  Output mode:'
        write(*, '(a)') '    1. Terminal ASCII plot (default)'
        write(*, '(a)') '    2. Export data (CSV)'
        write(*, '(a)', advance='no') '    Enter choice [1]: '

        read(*, '(a)', iostat=ios) input
        plot_mode = DOS_MODE_ASCII
        if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
            read(input, *, iostat=ios) plot_mode
            if (ios /= 0) plot_mode = DOS_MODE_ASCII
        end if

        ! ── Smearing ──
        smearing_width = 5.0_dp
        write(tmp_str, '(f5.1)') smearing_width
        write(*, '(a)', advance='no') '  Smearing width [' // trim(adjustl(tmp_str)) // ' cm-1]: '
        read(*, '(a)', iostat=ios) input
        if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
            read(input, *, iostat=ios) smearing_width
            if (ios /= 0) smearing_width = 5.0_dp
        end if
        smearing_width = max(0.5_dp, min(50.0_dp, smearing_width))

        ! ── Compute PHDOS ──
        freq_min_range = max(-200.0_dp, phdos%freq_min - 50.0_dp)
        freq_max_range = phdos%freq_max + 50.0_dp

        call compute_phonon_dos(phdos, freq_min_range, freq_max_range, 4001, smearing_width, ios)
        if (ios /= 0 .or. .not. allocated(phdos%phdos)) then
            write(*, '(a)') '  Error computing phonon DOS.'
            call free_phonon_dos_data(phdos); return
        end if

        ! ── Output dispatch ──
        select case (plot_mode)
        case (DOS_MODE_ASCII)
            call run_phonon_dos_navigator(phdos, freq_min_range, freq_max_range)
        case (DOS_MODE_EXPORT)
            call write_phonon_dos_csv(phdos)
        end select

        call free_phonon_dos_data(phdos)
    end subroutine handle_phonon_dos_menu

    subroutine handle_ir_menu(iostat)
        !! IR spectrum from .castep file (explicit column headers for IR/Raman)
        integer, intent(out) :: iostat
        type(vib_data_t) :: vib
        character(len=MAX_LINE_LEN) :: fname, input, tmp_str
        integer :: ios, plot_mode
        real(dp) :: freq_min, freq_max, freq_min_range, freq_max_range, smearing_width
        character(len=MAX_LINE_LEN), save :: last_castep_path = ''

        iostat = 0

        write(*, '(a)', advance='no') '  Enter .castep file path'
        if (len_trim(last_castep_path) > 0) &
            write(*, '(a)', advance='no') ' [' // trim(last_castep_path) // ']'
        write(*, '(a)') ': '

        read(*, '(a)', iostat=ios) fname
        if (ios /= 0) return
        fname = adjustl(fname); call strip_quotes(fname)
        if (fname == 'q' .or. fname == 'Q') then
            iostat = 0; return
        end if
        if (len_trim(fname) == 0 .and. len_trim(last_castep_path) > 0) then
            fname = last_castep_path
        end if
        if (len_trim(fname) == 0) then
            write(*, '(a)') '  No file specified.'; return
        end if
        last_castep_path = trim(fname)

        call parse_castep_vib(trim(fname), vib, ios)
        if (ios /= 0) then
            write(*, '(a,i0)') '  Error parsing .castep file, code=', ios
            return
        end if

        write(*, '(a,i0,a,i0,a)') '  Loaded ', vib%n_modes, ' Gamma modes, ', &
            count(vib%ir_intensity > 0.0_dp .and. vib%ir_active), ' IR-active.'

        ! ── Output mode menu ──
        write(*, '(a)') '  Output mode:'
        write(*, '(a)') '    1. Terminal ASCII plot (default)'
        write(*, '(a)') '    2. Export data (CSV)'
        write(*, '(a)', advance='no') '    Enter choice [1]: '
        read(*, '(a)', iostat=ios) input
        plot_mode = DOS_MODE_ASCII
        if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
            read(input, *, iostat=ios) plot_mode
            if (ios /= 0) plot_mode = DOS_MODE_ASCII
        end if

        ! ── Smearing ──
        smearing_width = 5.0_dp
        write(tmp_str, '(f5.1)') smearing_width
        write(*, '(a)', advance='no') '  Smearing width [' // trim(adjustl(tmp_str)) // ' cm-1]: '
        read(*, '(a)', iostat=ios) input
        if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
            read(input, *, iostat=ios) smearing_width
            if (ios /= 0) smearing_width = 5.0_dp
        end if
        smearing_width = max(0.5_dp, min(50.0_dp, smearing_width))

        ! ── Compute IR spectrum ──
        freq_min = minval(vib%freq)
        freq_max = maxval(vib%freq)
        freq_min_range = max(0.0_dp, freq_min - 50.0_dp)
        freq_max_range = freq_max + 50.0_dp
        call compute_ir_spectrum(vib, freq_min_range, freq_max_range, 4001, smearing_width, ios)
        if (ios /= 0 .or. .not. allocated(vib%ir_spectrum)) then
            write(*, '(a)') '  Error computing IR spectrum.'
            call free_vib_data(vib); return
        end if

        select case (plot_mode)
        case (DOS_MODE_ASCII)
            call run_ir_navigator(vib, freq_min_range, freq_max_range)
        case (DOS_MODE_EXPORT)
            call write_ir_csv(vib)
        end select

        call free_vib_data(vib)
    end subroutine handle_ir_menu

    subroutine run_ir_navigator(vib, freq_min_range, freq_max_range)
        type(vib_data_t), intent(in) :: vib
        real(dp), intent(in) :: freq_min_range, freq_max_range
        real(dp) :: e_center, half_range, y_center, y_half, y_half0, y_max_val
        integer :: i, tw, th, ios
        character(len=1) :: ch
        real(dp) :: dos_data(size(vib%ir_spectrum), 1)
        character(len=3) :: arrow

        do i = 1, size(vib%ir_spectrum)
            dos_data(i, 1) = vib%ir_spectrum(i)
        end do

        half_range = (freq_max_range - freq_min_range) * 0.5_dp
        e_center = (freq_min_range + freq_max_range) * 0.5_dp
        y_max_val = maxval(vib%ir_spectrum) * 1.15_dp
        y_half = y_max_val * 0.5_dp
        y_center = y_half
        y_half0 = y_half

        call enter_raw_mode
        do
            call get_term_size(tw, th)
            write(*, '(a)', advance='no') achar(27) // '[2J' // achar(27) // '[H'
            call plot_dos_ascii(vib%freq_grid, dos_data, 1, 0.0_dp, &
                vib%smearing, tw, th, &
                y_center_in=y_center, y_half_in=y_half, &
                e_center_in=e_center, half_range_in=half_range, &
                xlabel='Frequency', xunit='cm-1')
            write(*, '(a)') '  [arrows: pan  +/-: zoom  R: reset  Q: quit]'

            read(*, '(a)', advance='no', iostat=ios) ch
            if (ios /= 0) exit
            if (iachar(ch) == 27) then
                read(*, '(a)', advance='no', iostat=ios) arrow(1:1)
                if (ios /= 0) exit
                read(*, '(a)', advance='no', iostat=ios) arrow(2:2)
                if (ios /= 0) exit
                if (arrow(1:2) == '[A') then
                    y_center = y_center + y_half * 0.3_dp
                else if (arrow(1:2) == '[B') then
                    y_center = max(0.0_dp, y_center - y_half * 0.3_dp)
                else if (arrow(1:2) == '[C') then
                    e_center = e_center + half_range * 0.3_dp
                else if (arrow(1:2) == '[D') then
                    e_center = e_center - half_range * 0.3_dp
                end if
                cycle
            end if
            select case (ch)
            case ('+', '=')
                half_range = max(1.0_dp, half_range * 0.5_dp)
                y_half = max(0.001_dp, y_half * 0.5_dp)
            case ('-')
                half_range = min(5000.0_dp, half_range * 2.0_dp)
                y_half = y_half * 2.0_dp
            case ('r', 'R')
                e_center = (freq_min_range + freq_max_range) * 0.5_dp
                half_range = (freq_max_range - freq_min_range) * 0.5_dp
                y_center = y_half0; y_half = y_half0
            case ('q', 'Q')
                exit
            end select
        end do
        call leave_raw_mode
    end subroutine run_ir_navigator

    subroutine write_ir_csv(vib)
        type(vib_data_t), intent(in) :: vib
        character(len=MAX_LINE_LEN) :: csv_file
        integer :: unit, ios, i

        write(*, '(a)', advance='no') '  Enter output CSV file name (without .csv): '
        read(*, '(a)', iostat=ios) csv_file
        if (ios /= 0) return
        csv_file = adjustl(csv_file); call strip_quotes(csv_file)
        if (len_trim(csv_file) == 0) return
        csv_file = trim(csv_file) // '.csv'

        open(newunit=unit, file=trim(csv_file), status='unknown', action='write', iostat=ios)
        if (ios /= 0) then
            write(*, '(a)') '  Cannot write CSV file.'; return
        end if
        write(unit, '(a)') '# Frequency(cm-1),IR_Intensity'
        do i = 1, size(vib%ir_spectrum)
            write(unit, '(f12.4,a,es14.6)') vib%freq_grid(i), ',', vib%ir_spectrum(i)
        end do
        close(unit)
        write(*, '(a)') '  Written ' // trim(csv_file)
    end subroutine write_ir_csv

    subroutine handle_raman_menu(iostat)
        !! Raman spectrum from .castep file (explicit column headers for IR/Raman)
        integer, intent(out) :: iostat
        type(vib_data_t) :: vib
        character(len=MAX_LINE_LEN) :: fname, input, tmp_str
        integer :: ios, plot_mode
        real(dp) :: freq_min, freq_max, freq_min_range, freq_max_range, smearing_width
        character(len=MAX_LINE_LEN), save :: last_castep_path = ''

        iostat = 0

        write(*, '(a)', advance='no') '  Enter .castep file path'
        if (len_trim(last_castep_path) > 0) &
            write(*, '(a)', advance='no') ' [' // trim(last_castep_path) // ']'
        write(*, '(a)') ': '
        read(*, '(a)', iostat=ios) fname
        if (ios /= 0) return
        fname = adjustl(fname); call strip_quotes(fname)
        if (fname == 'q' .or. fname == 'Q') then
            iostat = 0; return
        end if
        if (len_trim(fname) == 0 .and. len_trim(last_castep_path) > 0) then
            fname = last_castep_path
        end if
        if (len_trim(fname) == 0) then
            write(*, '(a)') '  No file specified.'; return
        end if
        last_castep_path = trim(fname)

        call parse_castep_vib(trim(fname), vib, ios)
        if (ios /= 0) then
            write(*, '(a,i0)') '  Error parsing .castep file, code=', ios
            return
        end if

        if (.not. vib%has_raman_numeric) then
            write(*, '(a)') '  Warning: no numeric Raman activity column in .castep.'
            write(*, '(a)') '  This calculation did not compute Raman activities.'
            call free_vib_data(vib); return
        end if

        write(*, '(a,i0,a,i0,a)') '  Loaded ', vib%n_modes, ' Gamma modes, ', &
            count(vib%raman_activity > 0.0_dp .and. vib%raman_active), ' Raman-active.'

        ! ── Output mode menu ──
        write(*, '(a)') '  Output mode:'
        write(*, '(a)') '    1. Terminal ASCII plot (default)'
        write(*, '(a)') '    2. Export data (CSV)'
        write(*, '(a)', advance='no') '    Enter choice [1]: '
        read(*, '(a)', iostat=ios) input
        plot_mode = DOS_MODE_ASCII
        if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
            read(input, *, iostat=ios) plot_mode
            if (ios /= 0) plot_mode = DOS_MODE_ASCII
        end if

        ! ── Smearing ──
        smearing_width = 5.0_dp
        write(tmp_str, '(f5.1)') smearing_width
        write(*, '(a)', advance='no') '  Smearing width [' // trim(adjustl(tmp_str)) // ' cm-1]: '
        read(*, '(a)', iostat=ios) input
        if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
            read(input, *, iostat=ios) smearing_width
            if (ios /= 0) smearing_width = 5.0_dp
        end if
        smearing_width = max(0.5_dp, min(50.0_dp, smearing_width))

        ! ── Compute Raman spectrum ──
        freq_min = minval(vib%freq)
        freq_max = maxval(vib%freq)
        freq_min_range = max(0.0_dp, freq_min - 50.0_dp)
        freq_max_range = freq_max + 50.0_dp
        call compute_raman_spectrum(vib, freq_min_range, freq_max_range, 4001, smearing_width, ios)
        if (ios /= 0 .or. .not. allocated(vib%raman_spectrum)) then
            write(*, '(a)') '  Error computing Raman spectrum.'
            call free_vib_data(vib); return
        end if

        select case (plot_mode)
        case (DOS_MODE_ASCII)
            call run_raman_navigator(vib, freq_min_range, freq_max_range)
        case (DOS_MODE_EXPORT)
            call write_raman_csv(vib)
        end select
        call free_vib_data(vib)
    end subroutine handle_raman_menu

    subroutine run_raman_navigator(vib, freq_min_range, freq_max_range)
        type(vib_data_t), intent(in) :: vib
        real(dp), intent(in) :: freq_min_range, freq_max_range
        real(dp) :: e_center, half_range, y_center, y_half, y_half0, y_max_val
        integer :: i, tw, th, ios
        character(len=1) :: ch
        real(dp) :: dos_data(size(vib%raman_spectrum), 1)
        character(len=3) :: arrow

        do i = 1, size(vib%raman_spectrum)
            dos_data(i, 1) = vib%raman_spectrum(i)
        end do

        half_range = (freq_max_range - freq_min_range) * 0.5_dp
        e_center = (freq_min_range + freq_max_range) * 0.5_dp
        y_max_val = maxval(vib%raman_spectrum) * 1.15_dp
        y_half = y_max_val * 0.5_dp
        y_center = y_half
        y_half0 = y_half

        call enter_raw_mode
        do
            call get_term_size(tw, th)
            write(*, '(a)', advance='no') achar(27) // '[2J' // achar(27) // '[H'
            call plot_dos_ascii(vib%freq_grid, dos_data, 1, 0.0_dp, &
                vib%smearing, tw, th, &
                y_center_in=y_center, y_half_in=y_half, &
                e_center_in=e_center, half_range_in=half_range, &
                xlabel='Frequency', xunit='cm-1')
            write(*, '(a)') '  [arrows: pan  +/-: zoom  R: reset  Q: quit]'

            read(*, '(a)', advance='no', iostat=ios) ch
            if (ios /= 0) exit
            if (iachar(ch) == 27) then
                read(*, '(a)', advance='no', iostat=ios) arrow(1:1)
                if (ios /= 0) exit
                read(*, '(a)', advance='no', iostat=ios) arrow(2:2)
                if (ios /= 0) exit
                if (arrow(1:2) == '[A') then
                    y_center = y_center + y_half * 0.3_dp
                else if (arrow(1:2) == '[B') then
                    y_center = max(0.0_dp, y_center - y_half * 0.3_dp)
                else if (arrow(1:2) == '[C') then
                    e_center = e_center + half_range * 0.3_dp
                else if (arrow(1:2) == '[D') then
                    e_center = e_center - half_range * 0.3_dp
                end if
                cycle
            end if
            select case (ch)
            case ('+', '=')
                half_range = max(1.0_dp, half_range * 0.5_dp)
                y_half = max(0.001_dp, y_half * 0.5_dp)
            case ('-')
                half_range = min(5000.0_dp, half_range * 2.0_dp)
                y_half = y_half * 2.0_dp
            case ('r', 'R')
                e_center = (freq_min_range + freq_max_range) * 0.5_dp
                half_range = (freq_max_range - freq_min_range) * 0.5_dp
                y_center = y_half0; y_half = y_half0
            case ('q', 'Q')
                exit
            end select
        end do
        call leave_raw_mode
    end subroutine run_raman_navigator

    subroutine write_raman_csv(vib)
        type(vib_data_t), intent(in) :: vib
        character(len=MAX_LINE_LEN) :: csv_file
        integer :: unit, ios, i

        write(*, '(a)', advance='no') '  Enter output CSV file name (without .csv): '
        read(*, '(a)', iostat=ios) csv_file
        if (ios /= 0) return
        csv_file = adjustl(csv_file); call strip_quotes(csv_file)
        if (len_trim(csv_file) == 0) return
        csv_file = trim(csv_file) // '.csv'

        open(newunit=unit, file=trim(csv_file), status='unknown', action='write', iostat=ios)
        if (ios /= 0) then
            write(*, '(a)') '  Cannot write CSV file.'; return
        end if
        write(unit, '(a)') '# Frequency(cm-1),Raman_Activity'
        do i = 1, size(vib%raman_spectrum)
            write(unit, '(f12.4,a,es14.6)') vib%freq_grid(i), ',', vib%raman_spectrum(i)
        end do
        close(unit)
        write(*, '(a)') '  Written ' // trim(csv_file)
    end subroutine write_raman_csv

    subroutine run_phonon_dos_navigator(phdos, freq_min_range, freq_max_range)
        type(phonon_dos_data_t), intent(in) :: phdos
        real(dp), intent(in) :: freq_min_range, freq_max_range
        real(dp) :: e_center, half_range, y_center, y_half, y_half0
        real(dp) :: y_max_dos
        integer :: i, tw, th, ios
        character(len=1) :: ch
        real(dp) :: dos_data(size(phdos%phdos), 1)
        character(len=3) :: arrow

        do i = 1, size(phdos%phdos)
            dos_data(i, 1) = phdos%phdos(i)
        end do

        half_range = (freq_max_range - freq_min_range) * 0.5_dp
        e_center = (freq_min_range + freq_max_range) * 0.5_dp
        y_max_dos = maxval(phdos%phdos) * 1.15_dp
        y_half = y_max_dos * 0.5_dp
        y_center = y_half
        y_half0 = y_half

        call enter_raw_mode
        do
            call get_term_size(tw, th)
            write(*, '(a)', advance='no') achar(27) // '[2J' // achar(27) // '[H'
            call plot_dos_ascii(phdos%freq_grid, dos_data, 1, 0.0_dp, &
                phdos%smearing, tw, th, &
                y_center_in=y_center, y_half_in=y_half, &
                e_center_in=e_center, half_range_in=half_range, &
                xlabel='Frequency', xunit='cm-1')
            write(*, '(a)') '  [arrows: pan  +/-: zoom  R: reset  Q: quit]'

            read(*, '(a)', advance='no', iostat=ios) ch
            if (ios /= 0) exit

            ! Decode escape sequences for arrow keys
            if (iachar(ch) == 27) then
                read(*, '(a)', advance='no', iostat=ios) arrow(1:1)
                if (ios /= 0) exit
                read(*, '(a)', advance='no', iostat=ios) arrow(2:2)
                if (ios /= 0) exit
                if (arrow(1:2) == '[A') then       ! ↑
                    y_center = y_center + y_half * 0.3_dp
                else if (arrow(1:2) == '[B') then  ! ↓
                    y_center = max(0.0_dp, y_center - y_half * 0.3_dp)
                else if (arrow(1:2) == '[C') then  ! →
                    e_center = e_center + half_range * 0.3_dp
                else if (arrow(1:2) == '[D') then  ! ←
                    e_center = e_center - half_range * 0.3_dp
                end if
                cycle
            end if

            select case (ch)
            case ('+', '=')
                half_range = max(1.0_dp, half_range * 0.5_dp)
                y_half = max(0.001_dp, y_half * 0.5_dp)
            case ('-')
                half_range = min(5000.0_dp, half_range * 2.0_dp)
                y_half = y_half * 2.0_dp
            case ('r', 'R')
                e_center = (freq_min_range + freq_max_range) * 0.5_dp
                half_range = (freq_max_range - freq_min_range) * 0.5_dp
                y_center = y_half0
                y_half = y_half0
            case ('q', 'Q')
                exit
            end select
        end do
        call leave_raw_mode
    end subroutine run_phonon_dos_navigator

    subroutine write_phonon_dos_csv(phdos)
        type(phonon_dos_data_t), intent(in) :: phdos
        character(len=MAX_LINE_LEN) :: csv_file
        integer :: unit, ios, i

        write(*, '(a)', advance='no') '  Enter output CSV file name (without .csv): '
        read(*, '(a)', iostat=ios) csv_file
        if (ios /= 0) return
        csv_file = adjustl(csv_file); call strip_quotes(csv_file)
        if (len_trim(csv_file) == 0) return
        csv_file = trim(csv_file) // '.csv'

        open(newunit=unit, file=trim(csv_file), status='unknown', action='write', iostat=ios)
        if (ios /= 0) then
            write(*, '(a)') '  Cannot write CSV file.'
            return
        end if
        write(unit, '(a)') '# Frequency(cm-1),PHDOS'
        do i = 1, size(phdos%phdos)
            write(unit, '(f12.4,a,es14.6)') phdos%freq_grid(i), ',', phdos%phdos(i)
        end do
        close(unit)
        write(*, '(a)') '  Written ' // trim(csv_file)
    end subroutine write_phonon_dos_csv

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

        call enter_raw_mode
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

        call leave_raw_mode
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

        call enter_raw_mode
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

        call leave_raw_mode
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

    subroutine handle_polarizability_menu(iostat)
        !! Static polarizability via AIMD polarization fluctuation method
        integer, intent(out) :: iostat
        type(pol_data_t) :: pol
        character(len=MAX_LINE_LEN) :: castep_path, cp2k_dir, input
        character(len=256) :: msg
        integer :: ios
        real(dp) :: time_step_fs
        character(len=MAX_LINE_LEN), save :: last_castep_path = ''
        character(len=MAX_LINE_LEN), save :: last_cp2k_dir = ''

        iostat = 0

        ! --- CASTEP .castep file ---
        write(*, '(a)', advance='no') '  Enter CASTEP .castep file path'
        if (len_trim(last_castep_path) > 0) &
            write(*, '(a)', advance='no') ' [' // trim(last_castep_path) // ']'
        write(*, '(a)') ': '
        read(*, '(a)', iostat=ios) castep_path
        if (ios /= 0) return
        castep_path = adjustl(castep_path); call strip_quotes(castep_path)
        if (castep_path == 'q' .or. castep_path == 'Q') then
            iostat = 0; return
        end if
        if (len_trim(castep_path) == 0 .and. len_trim(last_castep_path) > 0) then
            castep_path = last_castep_path
        end if
        if (len_trim(castep_path) == 0) then
            write(*, '(a)') '  No file specified.'; return
        end if
        last_castep_path = trim(castep_path)

        ! Parse ε_∞ and cell parameters from .castep in one pass
        write(*, '(a)') '  Parsing CASTEP output...'
        call parse_castep_file(trim(castep_path), pol%eps_inf, pol%cell_abc, &
                                pol%volume_ang3, ios, msg)
        if (ios /= 0) then
            write(*, '(a,a)') '  Error: ', trim(msg)
            return
        end if
        write(*, '(a)') '  Optical dielectric tensor ε_∞ (from CASTEP):'
        call print_matrix(pol%eps_inf)
        write(*, '(a,3(f12.6,a))') '  Cell from .castep: a=', pol%cell_abc(1), &
            '  b=', pol%cell_abc(2), '  c=', pol%cell_abc(3), ' Å'
        write(*, '(a,f10.2,a)') '  Cell volume: ', pol%volume_ang3, ' ų'

        ! --- CP2K dipole directory ---
        write(*, '(a)', advance='no') '  Enter CP2K dipole directory path'
        if (len_trim(last_cp2k_dir) > 0) &
            write(*, '(a)', advance='no') ' [' // trim(last_cp2k_dir) // ']'
        write(*, '(a)') ': '
        read(*, '(a)', iostat=ios) cp2k_dir
        if (ios /= 0) return
        cp2k_dir = adjustl(cp2k_dir); call strip_quotes(cp2k_dir)
        if (cp2k_dir == 'q' .or. cp2k_dir == 'Q') then
            iostat = 0; return
        end if
        if (len_trim(cp2k_dir) == 0 .and. len_trim(last_cp2k_dir) > 0) then
            cp2k_dir = last_cp2k_dir
        end if
        if (len_trim(cp2k_dir) == 0) then
            write(*, '(a)') '  No directory specified.'; return
        end if
        last_cp2k_dir = trim(cp2k_dir)

        ! --- Temperature ---
        write(*, '(a)', advance='no') '  Enter temperature (K): '
        read(*, '(a)', iostat=ios) input
        if (ios /= 0) then
            write(*, '(a)') '  Invalid input.'; return
        end if
        call strip_quotes(input)
        if (input == 'q' .or. input == 'Q') then
            iostat = 0; return
        end if
        read(input, *, iostat=ios) pol%temperature
        if (ios /= 0 .or. pol%temperature <= 0.0_dp) then
            write(*, '(a)') '  Invalid temperature.'; return
        end if
        write(*, '(a,f8.1,a)') '  Temperature: ', pol%temperature, ' K'

        ! --- Time step ---
        write(*, '(a)', advance='no') '  Enter MD time step (fs) [1.0]: '
        read(*, '(a)', iostat=ios) input
        time_step_fs = 1.0_dp  ! default 1 fs
        if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
            if (adjustl(trim(input)) == 'q' .or. adjustl(trim(input)) == 'Q') then
                iostat = 0; return
            end if
            read(input, *, iostat=ios) time_step_fs
        end if

        ! --- Parse CP2K dipoles ---
        write(*, '(a)') '  Reading CP2K dipole files...'
        call parse_cp2k_dipoles(trim(cp2k_dir), pol, ios, msg)
        if (ios /= 0) then
            write(*, '(a,a)') '  Error: ', trim(msg)
            call free_pol_data(pol); return
        end if
        write(*, '(a,i0,a)') '  Loaded ', pol%n_frames, ' dipole frames.'

        ! --- Unwrap ---
        write(*, '(a)') '  Unwrapping Berry phase jumps...'
        call unwrap_dipoles(pol)
        write(*, '(a,i0,a)') '  Detected ', pol%n_unwraps, ' jumps.'

        ! --- Window-based dielectric (per-window detrend + W→0 extrapolation) ---
        write(*, '(a)') '  Computing static dielectric (window method)...'
        call compute_static_dielectric_windowed(pol, time_step_fs, ios, msg, verbose=.true.)
        if (ios /= 0) then
            write(*, '(a,a)') '  Error: ', trim(msg)
            call free_pol_data(pol); return
        end if

        ! Compute static dielectric and polarizability from windowed ε_ion
        pol%eps_static = pol%eps_inf + pol%eps_ion
        call compute_polarizability(pol)

        ! --- Print results ---
        write(*, '(a)') ''
        write(*, '(a)') '  ======= Results ======='
        write(*, '(a)') ''
        write(*, '(a)') '  Ionic dielectric tensor ε_ion:'
        call print_matrix(pol%eps_ion)
        write(*, '(a)') ''
        write(*, '(a)') '  Static dielectric tensor ε_static:'
        call print_matrix(pol%eps_static)
        write(*, '(a)') ''
        write(*, '(a,f10.2)') '  ε_∞,iso = ', trace_iso(pol%eps_inf)
        write(*, '(a,f10.2)') '  ε_ion,iso = ', trace_iso(pol%eps_ion)
        write(*, '(a,f10.2)') '  ε_static,iso = ', trace_iso(pol%eps_static)
        write(*, '(a)') ''
        write(*, '(a)') '  Static polarizability α_static (ų):'
        call print_matrix(pol%alpha_static)
        write(*, '(a)') ''
        write(*, '(a,f10.1)') '  α_∞,iso = ', trace_iso(pol%alpha_inf)
        write(*, '(a,f10.1)') '  α_ion,iso = ', trace_iso(pol%alpha_ion)
        write(*, '(a,f10.1)') '  α_static,iso = ', trace_iso(pol%alpha_static)
        write(*, '(a)') ''

        ! Wait for user to press q to return
        do
            write(*, '(a)', advance='no') '  Press q to return to menu: '
            read(*, '(a)', iostat=ios) input
            if (ios /= 0) exit
            if (adjustl(trim(input)) == 'q' .or. adjustl(trim(input)) == 'Q') exit
        end do

        call free_pol_data(pol)
    end subroutine handle_polarizability_menu


    subroutine print_matrix(mat)
        !! Print 3×3 matrix with formatting
        real(dp), intent(in) :: mat(3,3)
        integer :: i
        do i = 1, 3
            write(*, '(a,3f12.4)') '    ', mat(i, 1:3)
        end do
    end subroutine print_matrix


    function trace_iso(mat) result(val)
        !! Compute isotropic average: (mat_xx + mat_yy + mat_zz) / 3
        real(dp), intent(in) :: mat(3,3)
        real(dp) :: val
        val = (mat(1,1) + mat(2,2) + mat(3,3)) / 3.0_dp
    end function trace_iso

    ! ----------------------------------------------------------------
    !  Universal Format Converter (Option 0): CIF ↔ PDB ↔ CELL
    ! ----------------------------------------------------------------

    subroutine handle_format_converter_menu(iostat)
        !! Convert between CIF, PDB, and CASTEP .cell formats
        integer, intent(out) :: iostat
        type(cif_data_t) :: cif
        character(len=512)  :: in_path, out_path, stem
        character(len=256)  :: msg
        character(len=4)    :: in_fmt
        integer :: ios, istat, out_choice
        logical :: exists
        character(len=MAX_LINE_LEN) :: input

        iostat = 0

        write(*, '(a)', advance='no') '  Enter input file (.cif/.pdb/.cell): '
        read(*, '(a)', iostat=ios) input
        if (ios /= 0) return
        in_path = adjustl(trim(input))
        call strip_quotes(in_path)
        if (in_path == 'q' .or. in_path == 'Q') then
            iostat = 0; return
        end if
        if (len_trim(in_path) == 0) then
            write(*, '(a)') '  No file specified.'; return
        end if

        inquire(file=trim(in_path), exist=exists)
        if (.not. exists) then
            write(*, '(a,a)') '  File not found: ', trim(in_path)
            return
        end if

        ! Auto-detect format from extension
        in_fmt = get_ext_lower(in_path)
        select case (trim(in_fmt))
        case ('cif');  continue
        case ('pdb');  continue
        case ('cell'); continue
        case default
            write(*, '(a)') '  Unsupported format. Use .cif, .pdb, or .cell files.'
            return
        end select
        write(*, '(a,a)') '  Input format: ', trim(in_fmt)

        ! Parse
        write(*, '(a)') '  Parsing...'
        select case (trim(in_fmt))
        case ('cif');  call parse_cif_inline(trim(in_path), cif, istat, iomsg=msg)
        case ('pdb');  call parse_pdb_inline(trim(in_path), cif, istat, iomsg=msg)
        case ('cell'); call parse_cell_inline(trim(in_path), cif, istat, iomsg=msg)
        end select
        if (istat /= 0) then
            write(*, '(a,a)') '  Error: ', trim(msg)
            return
        end if

        ! ── Symmetry expansion for format conversion ──
        call expand_cif_symmetry(cif, istat)

        if (cif%n_atoms == 0) then
            write(*, '(a)') '  Warning: no atoms found.'
        end if

        write(*, '(a)')       '  ------- Cell Summary -------'
        write(*, '(a,3f10.4)') '  a,b,c (Ang):       ', cif%a, cif%b, cif%c
        write(*, '(a,3f10.4)') '  alpha,beta,gamma:  ', cif%alpha, cif%beta, cif%gamma
        write(*, '(a,i0)')     '  Atoms:             ', cif%n_atoms
        if (len_trim(cif%space_group) > 0 .and. cif%space_group /= 'P1') &
            write(*, '(a,a)')  '  Space group:       ', trim(cif%space_group)
        if (cif%positions_fractional) then
            write(*, '(a)')    '  Coordinates:       fractional'
        else
            write(*, '(a)')    '  Coordinates:       Cartesian'
        end if

        ! Output format selection
        write(*, '(a)') ''
        write(*, '(a)') '  Select output format:'
        write(*, '(a)') '    1. CASTEP .cell'
        write(*, '(a)') '    2. CIF  (.cif)'
        write(*, '(a)') '    3. PDB  (.pdb)'
        write(*, '(a)') '    Q. Back'
        write(*, '(a)', advance='no') '    Enter choice: '
        read(*, '(a)', iostat=ios) input
        if (ios /= 0) return
        if (input(1:1) == 'q' .or. input(1:1) == 'Q') then
            iostat = 0; return
        end if
        read(input, *, iostat=ios) out_choice
        if (ios /= 0 .or. out_choice < 1 .or. out_choice > 3) then
            write(*, '(a)') '  Invalid choice.'; return
        end if

        ! Derive output filename
        stem = get_file_stem(in_path)

        ! Convert coordinates + write
        select case (out_choice)
        case (1)  ! → CELL (POSITIONS_ABS, Cartesian)
            if (cif%positions_fractional) then
                call convert_to_cartesian(cif, istat)
                if (istat /= 0) then
                    write(*, '(a)') '  Error: zero-volume cell.'; return
                end if
            end if
            out_path = trim(stem) // '.cell'
            call write_cell_simple(cif, trim(out_path), istat, msg)

        case (2)  ! → CIF (fractional)
            if (.not. cif%positions_fractional) then
                call convert_to_fractional(cif, istat)
                if (istat /= 0) then
                    write(*, '(a)') '  Error: zero-volume cell.'; return
                end if
            end if
            out_path = trim(stem) // '.cif'
            call write_cif_file(cif, trim(out_path), istat, msg)

        case (3)  ! → PDB (Cartesian)
            if (cif%positions_fractional) then
                call convert_to_cartesian(cif, istat)
                if (istat /= 0) then
                    write(*, '(a)') '  Error: zero-volume cell.'; return
                end if
            end if
            out_path = trim(stem) // '.pdb'
            call write_pdb_file(cif, trim(out_path), istat, msg)
        end select

        if (istat /= 0) then
            write(*, '(a,a)') '  Error: ', trim(msg)
        else
            write(*, '(a,a)') '  File written: ', trim(out_path)
        end if
        write(*, '(a)') '  ----------------------------'
    end subroutine handle_format_converter_menu


    subroutine convert_to_fractional(cif, iostat)
        !! Convert Cartesian atom coordinates to fractional using inverse lattice
        type(cif_data_t), intent(inout) :: cif
        integer, intent(out) :: iostat
        real(dp) :: lattice(3,3), inv_lattice(3,3), frac(3)
        integer :: i

        iostat = 0
        lattice = compute_cartesian_lattice(cif%a, cif%b, cif%c, &
            cif%alpha, cif%beta, cif%gamma)
        inv_lattice = invert_lattice_3x3(lattice)

        ! Check for zero determinant (degenerate cell)
        if (all(abs(inv_lattice) < 1.0e-30_dp)) then
            iostat = 1; return
        end if

        do i = 1, cif%n_atoms
            frac = cartesian_to_fractional(cif%atoms(i)%x, cif%atoms(i)%y, &
                cif%atoms(i)%z, inv_lattice)
            cif%atoms(i)%x = frac(1)
            cif%atoms(i)%y = frac(2)
            cif%atoms(i)%z = frac(3)
        end do
    end subroutine convert_to_fractional


    subroutine write_cif_file(data, filename, iostat, iomsg)
        !! Write cif_data_t as a standard CIF file
        type(cif_data_t), intent(in) :: data
        character(len=*), intent(in) :: filename
        integer, intent(out) :: iostat
        character(len=*), intent(out) :: iomsg

        integer :: iunit, ios, i
        character(len=10) :: date_str
        character(len=8)  :: date_val
        character(len=128) :: data_name
        integer :: n

        iostat = 0; iomsg = ''

        ! Derive data_ name from filename stem
        data_name = filename
        n = len_trim(data_name)
        do while (n > 0)
            if (data_name(n:n) == '/') then; data_name = data_name(n+1:); exit; end if
            n = n - 1
        end do
        n = len_trim(data_name)
        do while (n > 0)
            if (data_name(n:n) == '.') then; data_name = data_name(1:n-1); exit; end if
            n = n - 1
        end do
        data_name = adjustl(data_name)

        open(newunit=iunit, file=trim(filename), status='replace', &
            action='write', iostat=ios)
        if (ios /= 0) then
            iostat = IO_WRITE_ERROR
            iomsg = 'Cannot open output file: ' // trim(filename)
            return
        end if

        call date_and_time(date=date_val)
        date_str = date_val(1:4) // '-' // date_val(5:6) // '-' // date_val(7:8)

        write(iunit, '(a,a)')       'data_', trim(data_name)
        write(iunit, '(a,a)')       '_audit_creation_date              ', trim(date_str)
        write(iunit, '(a)')         '_audit_creation_method            ''CASTEP Suite PosCASTEP'''
        write(iunit, '(a,a,a)')     "_symmetry_space_group_name_H-M    '", &
            trim(data%space_group), "'"
        write(iunit, '(a,f12.4)')  '_cell_length_a                    ', data%a
        write(iunit, '(a,f12.4)')  '_cell_length_b                    ', data%b
        write(iunit, '(a,f12.4)')  '_cell_length_c                    ', data%c
        write(iunit, '(a,f12.4)')  '_cell_angle_alpha                 ', data%alpha
        write(iunit, '(a,f12.4)')  '_cell_angle_beta                  ', data%beta
        write(iunit, '(a,f12.4)')  '_cell_angle_gamma                 ', data%gamma
        write(iunit, '(a)')        'loop_'
        write(iunit, '(a)')        '_atom_site_label'
        write(iunit, '(a)')        '_atom_site_type_symbol'
        write(iunit, '(a)')        '_atom_site_fract_x'
        write(iunit, '(a)')        '_atom_site_fract_y'
        write(iunit, '(a)')        '_atom_site_fract_z'
        do i = 1, data%n_atoms
            write(iunit, '(a,2x,a,2x,f12.6,2x,f12.6,2x,f12.6)') &
                trim(data%atoms(i)%label), trim(data%atoms(i)%element), &
                data%atoms(i)%x, data%atoms(i)%y, data%atoms(i)%z
        end do

        close(iunit)
    end subroutine write_cif_file




    pure function invert_lattice_3x3(m) result(inv)
        !! Compute the inverse of a 3x3 matrix (adjugate / determinant)
        real(dp), intent(in) :: m(3,3)
        real(dp) :: inv(3,3)
        real(dp) :: det

        det = m(1,1)*(m(2,2)*m(3,3) - m(2,3)*m(3,2)) &
            - m(1,2)*(m(2,1)*m(3,3) - m(2,3)*m(3,1)) &
            + m(1,3)*(m(2,1)*m(3,2) - m(2,2)*m(3,1))

        if (abs(det) < 1.0e-12_dp) then
            inv = 0.0_dp
            return
        end if

        inv(1,1) =  (m(2,2)*m(3,3) - m(2,3)*m(3,2)) / det
        inv(1,2) = -(m(1,2)*m(3,3) - m(1,3)*m(3,2)) / det
        inv(1,3) =  (m(1,2)*m(2,3) - m(1,3)*m(2,2)) / det
        inv(2,1) = -(m(2,1)*m(3,3) - m(2,3)*m(3,1)) / det
        inv(2,2) =  (m(1,1)*m(3,3) - m(1,3)*m(3,1)) / det
        inv(2,3) = -(m(1,1)*m(2,3) - m(1,3)*m(2,1)) / det
        inv(3,1) =  (m(2,1)*m(3,2) - m(2,2)*m(3,1)) / det
        inv(3,2) = -(m(1,1)*m(3,2) - m(1,2)*m(3,1)) / det
        inv(3,3) =  (m(1,1)*m(2,2) - m(1,2)*m(2,1)) / det
    end function invert_lattice_3x3


    pure function cartesian_to_fractional(x, y, z, inv_lattice) result(frac)
        !! Convert Cartesian coordinates to fractional using inverse lattice matrix
        real(dp), intent(in) :: x, y, z, inv_lattice(3,3)
        real(dp) :: frac(3)
        frac(1) = inv_lattice(1,1)*x + inv_lattice(1,2)*y + inv_lattice(1,3)*z
        frac(2) = inv_lattice(2,1)*x + inv_lattice(2,2)*y + inv_lattice(2,3)*z
        frac(3) = inv_lattice(3,1)*x + inv_lattice(3,2)*y + inv_lattice(3,3)*z
    end function cartesian_to_fractional


    subroutine convert_to_cartesian(cif, iostat)
        !! Convert fractional atom coordinates to Cartesian using lattice matrix
        type(cif_data_t), intent(inout) :: cif
        integer, intent(out) :: iostat
        real(dp) :: lattice(3,3), cart(3)
        integer :: i

        iostat = 0
        lattice = compute_cartesian_lattice(cif%a, cif%b, cif%c, &
            cif%alpha, cif%beta, cif%gamma)

        do i = 1, cif%n_atoms
            cart(1) = lattice(1,1)*cif%atoms(i)%x + lattice(1,2)*cif%atoms(i)%y &
                      + lattice(1,3)*cif%atoms(i)%z
            cart(2) = lattice(2,1)*cif%atoms(i)%x + lattice(2,2)*cif%atoms(i)%y &
                      + lattice(2,3)*cif%atoms(i)%z
            cart(3) = lattice(3,1)*cif%atoms(i)%x + lattice(3,2)*cif%atoms(i)%y &
                      + lattice(3,3)*cif%atoms(i)%z
            cif%atoms(i)%x = cart(1)
            cif%atoms(i)%y = cart(2)
            cif%atoms(i)%z = cart(3)
        end do
        cif%positions_fractional = .false.
    end subroutine convert_to_cartesian


    subroutine write_cell_simple(data, filename, iostat, iomsg)
        !! Write a minimal .cell file from cif_data_t (LATTICE_ABC + POSITIONS_ABS)
        type(cif_data_t), intent(in) :: data
        character(len=*), intent(in) :: filename
        integer, intent(out) :: iostat
        character(len=*), intent(out) :: iomsg
        integer :: iunit, ios, i

        iostat = 0; iomsg = ''

        open(newunit=iunit, file=trim(filename), status='replace', &
            action='write', iostat=ios)
        if (ios /= 0) then
            iostat = IO_WRITE_ERROR
            iomsg = 'Cannot open output file: ' // trim(filename)
            return
        end if

        write(iunit, '(a)') '%BLOCK LATTICE_ABC'
        write(iunit, '(3(f12.7,1x))') data%a, data%b, data%c
        write(iunit, '(3(f12.7,1x))') data%alpha, data%beta, data%gamma
        write(iunit, '(a)') '%ENDBLOCK LATTICE_ABC'
        write(iunit, '(a)') ''
        write(iunit, '(a)') '%BLOCK POSITIONS_ABS'
        do i = 1, data%n_atoms
            write(iunit, '(a,2x,3(f12.6,1x))') trim(data%atoms(i)%element), &
                data%atoms(i)%x, data%atoms(i)%y, data%atoms(i)%z
        end do
        write(iunit, '(a)') '%ENDBLOCK POSITIONS_ABS'
        close(iunit)
    end subroutine write_cell_simple


    subroutine write_pdb_file(data, filename, iostat, iomsg)
        !! Write PDB format from cif_data_t (coordinates must be Cartesian)
        type(cif_data_t), intent(in) :: data
        character(len=*), intent(in) :: filename
        integer, intent(out) :: iostat
        character(len=*), intent(out) :: iomsg
        integer :: iunit, ios, i
        character(len=6) :: atom_label
        character(len=80) :: line

        iostat = 0; iomsg = ''

        open(newunit=iunit, file=trim(filename), status='replace', &
            action='write', iostat=ios)
        if (ios /= 0) then
            iostat = IO_WRITE_ERROR
            iomsg = 'Cannot open output file: ' // trim(filename)
            return
        end if

        ! TITLE
        write(iunit, '(a)') 'TITLE     Converted by CASTEP Suite PosCASTEP'

        ! CRYST1: a(F9.3) b(F9.3) c(F9.3) alpha(F7.2) beta(F7.2) gamma(F7.2) spg(11x) z(4x)
        write(iunit, '(a,3f9.3,3f7.2,1x,a,1x,i4)') 'CRYST1', &
            data%a, data%b, data%c, &
            data%alpha, data%beta, data%gamma, &
            'P 1', 1

        ! ATOM records
        do i = 1, data%n_atoms
            atom_label = adjustl(data%atoms(i)%element)
            write(line, '(a,i5,1x,a4,1x,a4,1x,a1,1x,i4,4x,3f8.3,2f6.2,6x,a2)') &
                'ATOM  ', mod(i, 99999), atom_label(1:4), 'MOL ', 'A', &
                mod(i, 9999), data%atoms(i)%x, data%atoms(i)%y, data%atoms(i)%z, &
                1.0, 0.0, atom_label(1:2)
            write(iunit, '(a)') trim(line)
        end do

        write(iunit, '(a)') 'END'
        close(iunit)
    end subroutine write_pdb_file


    pure function get_file_stem(path) result(stem)
        !! Extract basename without directory and extension
        character(len=*), intent(in) :: path
        character(len=512) :: stem
        integer :: n

        stem = trim(path)
        ! Strip directory prefix
        n = len_trim(stem)
        do while (n > 0)
            if (stem(n:n) == '/') then; stem = stem(n+1:); exit; end if
            n = n - 1
        end do
        ! Strip extension
        n = len_trim(stem)
        do while (n > 0)
            if (stem(n:n) == '.') then; stem = stem(1:n-1); exit; end if
            n = n - 1
        end do
    end function get_file_stem


    pure function get_ext_lower(path) result(ext)
        !! Get lowercase file extension (e.g. 'cif', 'pdb', 'cell')
        character(len=*), intent(in) :: path
        character(len=4) :: ext
        integer :: i, n

        ext = ''
        n = len_trim(path)
        ! Find last dot
        i = n
        do while (i > 0)
            if (path(i:i) == '.') then
                ext = path(i+1:n)
                exit
            end if
            i = i - 1
        end do
        ! Convert to lowercase
        do i = 1, len_trim(ext)
            if (ext(i:i) >= 'A' .and. ext(i:i) <= 'Z') &
                ext(i:i) = char(iachar(ext(i:i)) + 32)
        end do
    end function get_ext_lower


    subroutine launch_viewer(json_path)
        character(len=*), intent(in) :: json_path
        character(len=MAX_LINE_LEN) :: viewer_cmd
        integer :: cmdstat
        viewer_cmd = find_viewer()
        if (len_trim(viewer_cmd) == 0) then
            write(*, '(a)') '  crystal-viewer not found.'
            write(*, '(a)') '  Build it: cd crystal-viewer && cargo build --release'
            write(*, '(a)') '  Or download from: https://github.com/Kiwi-ucas/CASTEP_Suite/releases'
        else
            write(*, '(a)') '  Launching Crystal Viewer...'
            call execute_command_line(trim(viewer_cmd) // ' "' // trim(json_path) // '" > /dev/null 2> /dev/null', &
                wait=.true., cmdstat=cmdstat)
        end if
    end subroutine launch_viewer


    function find_viewer() result(cmd)
        !! Auto-detect crystal-viewer relative to CASTEP_Suite executable.
        !! Resolves the real executable path so it works regardless of
        !! which directory the program was launched from.
        character(len=MAX_LINE_LEN) :: cmd
        character(len=MAX_LINE_LEN) :: exe_path, exe_dir
        integer :: slash_pos, unit, ios
        logical :: exists, is_dir

        ! Resolve real absolute path:
        ! 1) Use readlink on /proc/self/exe from a subshell where $PPID is our PID
        call execute_command_line( &
            'readlink -f /proc/$PPID/exe > /tmp/_fv_tmp 2>/dev/null', exitstat=ios)
        if (ios == 0) then
            open(newunit=unit, file='/tmp/_fv_tmp', status='old', action='read', iostat=ios)
            if (ios == 0) then
                read(unit, '(a)', iostat=ios) exe_path
                close(unit, status='delete')
                if (ios == 0 .and. len_trim(exe_path) > 0) goto 10
            end if
        end if
        ! 2) macOS: use lsof to get real executable path
        call execute_command_line( &
            'lsof -p $PPID -a -d txt -Fn 2>/dev/null | grep ''^n'' | head -1 | cut -c2- > /tmp/_fv_tmp', &
            exitstat=ios)
        if (ios == 0) then
            open(newunit=unit, file='/tmp/_fv_tmp', status='old', action='read', iostat=ios)
            if (ios == 0) then
                read(unit, '(a)', iostat=ios) exe_path
                close(unit, status='delete')
                if (ios == 0 .and. len_trim(exe_path) > 0) goto 10
            end if
        end if
        ! 3) Fallback: command argument → make absolute with PWD
        call get_command_argument(0, exe_path)
        if (exe_path(1:1) /= '/') then
            call get_environment_variable('PWD', cmd, ios)
            if (ios > 0) exe_path = trim(cmd) // '/' // adjustl(exe_path)
        end if
10      continue
        slash_pos = index(trim(exe_path), '/', back=.true.)
        if (slash_pos > 0) then
            exe_dir = exe_path(1:slash_pos)
        else
            exe_dir = './'
        end if
        ! 1) Dev layout: crystal-viewer/target/release/crystal-viewer
        cmd = trim(exe_dir) // 'crystal-viewer/target/release/crystal-viewer'
        inquire(file=trim(cmd), exist=exists)
        if (exists) return
        ! 2) Release layout: crystal-viewer in same directory (skip if dir)
        cmd = trim(exe_dir) // 'crystal-viewer'
        inquire(file=trim(cmd), exist=exists)
        if (exists) then
            inquire(file=trim(cmd) // '/.', exist=is_dir)
            if (.not. is_dir) return
        end if
        ! 3) Current directory — dev layout
        cmd = './crystal-viewer/target/release/crystal-viewer'
        inquire(file=trim(cmd), exist=exists)
        if (exists) return
        ! 4) Current directory — bare binary (skip if dir)
        cmd = './crystal-viewer'
        inquire(file=trim(cmd), exist=exists)
        if (exists) then
            inquire(file=trim(cmd) // '/.', exist=is_dir)
            if (.not. is_dir) return
        end if
        cmd = ''
    end function find_viewer


    subroutine handle_view_structure(iostat)
        !! Parse a CIF/PDB/cell file, launch crystal-viewer, detect modifications
        integer, intent(out) :: iostat
        type(cif_data_t) :: cif
        character(len=MAX_LINE_LEN) :: file_path, json_path, ext
        character(len=256) :: msg
        integer :: ios, unit
        logical :: exists, modified

        iostat = 0

        ! ── File input ──
        write(*, '(a)', advance='no') '  Enter structure file (.cif/.pdb/.cell): '
        read(*, '(a)', iostat=ios) file_path
        if (ios /= 0) return
        file_path = adjustl(file_path); call strip_quotes(file_path)
        if (file_path == 'q' .or. file_path == 'Q') then
            iostat = 0; return
        end if
        if (len_trim(file_path) == 0) then
            write(*, '(a)') '  No file specified.'; return
        end if

        ! Check file exists
        inquire(file=trim(file_path), exist=exists)
        if (.not. exists) then
            write(*, '(a)') '  File not found: ' // trim(file_path); return
        end if

        ! ── Detect format ──
        ext = get_ext_lower(file_path)
        select case (trim(ext))
        case ('cif')
            call parse_cif_inline(trim(file_path), cif, ios)
        case ('pdb')
            call parse_pdb_inline(trim(file_path), cif, ios)
        case ('cell')
            call parse_cell_inline(trim(file_path), cif, ios)
        case default
            write(*, '(a)') '  Unsupported format. Use .cif, .pdb, or .cell.'
            return
        end select

        if (ios /= 0) then
            write(*, '(a)') '  Error parsing file.'
            call free_cif_data(cif); return
        end if

        write(*, '(a,i0,a)') '  Parsed ', cif%n_atoms, ' atom(s) in asymmetric unit.'

        ! ── Symmetry expansion for viewer display ──
        if (cif%n_symops > 1) then
            call expand_cif_symmetry(cif, ios)
            if (ios == 0) &
                write(*, '(a,i0,a)') '  Expanded to ', cif%n_atoms, ' atom(s) in full unit cell'
        end if

        ! Remember source file for PreCASTEP handoff and re-edit
        precastep_source_file = get_file_stem(file_path)
        precastep_viewer_file = trim(file_path)

        ! ── Write JSON ──
        json_path = trim(file_path) // '.json'
        call write_crystal_json_cif(cif, json_path, ios)
        call free_cif_data(cif)

        if (ios /= 0) then
            write(*, '(a)') '  Error writing JSON file.'; return
        end if

        ! ── Launch viewer ──
        call launch_viewer(json_path)

        ! ── Read back modified structure ──
        call read_crystal_json_to_cif(json_path, cif, modified, ios, msg)
        if (ios /= 0) then
            write(*, '(a,a)') '  Error reading modified structure: ', trim(msg)
            ! Clean up JSON and return
            open(newunit=unit, file=trim(json_path), status='old', iostat=ios)
            if (ios == 0) close(unit, status='delete')
            return
        end if

        if (modified) then
            write(*, '(a)') '  Structure was modified in viewer.'
            call modified_structure_menu(cif, json_path, iostat)
        else
            ! No modifications — delete JSON and return
            open(newunit=unit, file=trim(json_path), status='old', iostat=ios)
            if (ios == 0) close(unit, status='delete')
            call free_cif_data(cif)
            iostat = 0
        end if

    end subroutine handle_view_structure


    subroutine modified_structure_menu(cif, json_path, iostat)
        !! Sub-menu shown when viewer modified atom positions.
        !! Option 1: save to file (CIF/PDB/cell)
        !! Option 2: hand cif_data_t directly to PreCASTEP (no file on disk)
        !! Option 3: continue editing current structure
        !! Option 4: re-edit from original structure file
        type(cif_data_t), intent(inout) :: cif
        character(len=*), intent(in) :: json_path
        integer, intent(out) :: iostat

        character(len=MAX_LINE_LEN) :: input, out_path, ext_str
        integer :: ios, choice, fmt_choice, menu_choice, unit
        logical :: modified

        iostat = 0

        do
            write(*, '(a)') ''
            write(*, '(a)') '  ================================'
            write(*, '(a)') '    Structure has been modified'
            write(*, '(a)') '  ================================'
            write(*, '(a)') '  1. Save new structure'
            write(*, '(a)') '  2. Use new structure to PreCASTEP'
            write(*, '(a)') '  3. Continue to edit structure'
            write(*, '(a)') '  4. Re-edit from original structure'
            write(*, '(a)') '  Q. Back'
            write(*, '(a)', advance='no') '  Select option: '

            read(*, '(a)', iostat=ios) input
            if (ios /= 0) exit

            if (len_trim(input) >= 1) then
                if (input(1:1) == 'q' .or. input(1:1) == 'Q') then
                    iostat = 0
                    exit
                end if
            end if

            read(input, '(I6)', iostat=ios) choice
            if (ios /= 0) then
                write(*, '(a)') '  Invalid input. Enter 1-4, or Q.'
                cycle
            end if

            if (choice == 1 .or. choice == 2) then
                menu_choice = choice  ! save before format selection overwrites it

                if (menu_choice == 2) then
                    ! Option 2: hand cif_data_t directly to PreCASTEP (no file on disk)
                    call copy_cif_data(cif, precastep_cif_data)
                    has_precastep_data = .true.
                    iostat = IO_PRECASTEP_LAUNCH
                    exit
                end if

                ! ── Option 1: Choose output format ──
                write(*, '(a)') '  Select output format:'
                write(*, '(a)') '    1. CASTEP .cell'
                write(*, '(a)') '    2. CIF  (.cif)'
                write(*, '(a)') '    3. PDB  (.pdb)'
                write(*, '(a)', advance='no') '    Enter choice: '
                read(*, '(a)', iostat=ios) input
                if (ios /= 0) exit
                read(input, *, iostat=ios) fmt_choice
                if (ios /= 0 .or. fmt_choice < 1 .or. fmt_choice > 3) then
                    write(*, '(a)') '  Invalid choice.'
                    cycle
                end if

                ! ── Output filename ──
                write(*, '(a)', advance='no') '  Enter output filename: '
                read(*, '(a)', iostat=ios) out_path
                if (ios /= 0) exit
                out_path = adjustl(out_path)
                call strip_quotes(out_path)
                if (len_trim(out_path) == 0) then
                    write(*, '(a)') '  No filename specified.'
                    cycle
                end if

                ! Convert to appropriate coordinate system and write
                select case (fmt_choice)
                case (1)  ! .cell (Cartesian)
                    if (cif%positions_fractional) call convert_to_cartesian(cif, ios)
                    call ensure_ext(out_path, '.cell')
                    call write_cell_simple(cif, trim(out_path), ios, input)
                case (2)  ! .cif (fractional)
                    if (.not. cif%positions_fractional) call convert_to_fractional(cif, ios)
                    call ensure_ext(out_path, '.cif')
                    call write_cif_file(cif, trim(out_path), ios, input)
                case (3)  ! .pdb (Cartesian)
                    if (cif%positions_fractional) call convert_to_cartesian(cif, ios)
                    call ensure_ext(out_path, '.pdb')
                    call write_pdb_file(cif, trim(out_path), ios, input)
                end select

                if (ios /= 0) then
                    write(*, '(a)') '  Error writing file.'
                    cycle
                end if
                write(*, '(a,a)') '  Saved as ', trim(out_path)

                iostat = 0
                exit
            else if (choice == 3) then
                ! Option 3: Continue editing — re-launch viewer with current cif
                call write_crystal_json_cif(cif, json_path, ios)
                if (ios /= 0) then
                    write(*, '(a)') '  Error writing JSON file.'; cycle
                end if
                call launch_viewer(json_path)
                call read_crystal_json_to_cif(json_path, cif, modified, ios, input)
                if (ios /= 0) then
                    write(*, '(a,a)') '  Error reading modified structure: ', trim(input)
                    cycle
                end if
                if (.not. modified) then
                    write(*, '(a)') '  No changes made in viewer.'
                    cycle
                end if
                write(*, '(a)') '  Structure was modified in viewer.'
                ! Loop back to menu with updated cif
            else if (choice == 4) then
                ! Option 4: Re-edit from original — re-parse original file
                call free_cif_data(cif)
                ext_str = get_ext_lower(precastep_viewer_file)
                select case (trim(ext_str))
                case ('cif');  call parse_cif_inline(trim(precastep_viewer_file), cif, ios)
                case ('pdb');  call parse_pdb_inline(trim(precastep_viewer_file), cif, ios)
                case ('cell'); call parse_cell_inline(trim(precastep_viewer_file), cif, ios)
                case default
                    write(*, '(a)') '  Unknown original file format.'; cycle
                end select
                if (ios /= 0 .or. cif%n_atoms == 0) then
                    write(*, '(a)') '  Error re-reading original file.'; exit
                end if
                ! ── Symmetry expansion for re-edit ──
                call expand_cif_symmetry(cif, ios)
                write(*, '(a,i0,a)') '  Re-loaded ', cif%n_atoms, ' atoms from original.'
                call write_crystal_json_cif(cif, json_path, ios)
                if (ios /= 0) then
                    write(*, '(a)') '  Error writing JSON file.'; cycle
                end if
                call launch_viewer(json_path)
                call read_crystal_json_to_cif(json_path, cif, modified, ios, input)
                if (ios /= 0) then
                    write(*, '(a,a)') '  Error reading modified structure: ', trim(input)
                    cycle
                end if
                if (.not. modified) then
                    write(*, '(a)') '  No changes made in viewer.'
                    cycle
                end if
                write(*, '(a)') '  Structure was modified in viewer.'
                ! Loop back to menu with updated cif
            else
                write(*, '(a)') '  Invalid option. Enter 1-4, or Q.'
            end if
        end do

        ! Clean up JSON (only if we're exiting the loop)
        open(newunit=unit, file=trim(json_path), status='old', iostat=ios)
        if (ios == 0) close(unit, status='delete')

    end subroutine modified_structure_menu


    subroutine copy_cif_data(src, dst)
        !! Deep copy cif_data_t (including allocatable atoms array)
        type(cif_data_t), intent(in) :: src
        type(cif_data_t), intent(out) :: dst
        integer :: n, istat

        dst%a = src%a
        dst%b = src%b
        dst%c = src%c
        dst%alpha = src%alpha
        dst%beta = src%beta
        dst%gamma = src%gamma
        dst%space_group = src%space_group
        dst%n_atoms = src%n_atoms
        dst%positions_fractional = src%positions_fractional

        if (allocated(src%atoms)) then
            n = size(src%atoms)
            allocate(dst%atoms(n), stat=istat)
            if (istat == 0) dst%atoms(1:n) = src%atoms(1:n)
        end if
    end subroutine copy_cif_data


    ! ====================================================================
    ! Option 8: Phonon Mode Visualization
    ! ====================================================================
    subroutine handle_phonon_modes_menu(iostat)
        integer, intent(out) :: iostat
        type(phonon_modes_data_t) :: modes_data
        character(len=MAX_LINE_LEN) :: fname, castep_path, input, json_path, msg
        integer :: ios, mode_idx, best_mode
        real(dp) :: max_ir

        iostat = 0

        ! ── Prompt for .phonon file ──
        write(*, '(a)', advance='no') '  Enter .phonon file path: '
        read(*, '(a)', iostat=ios) fname
        if (ios /= 0) return
        fname = adjustl(fname); call strip_quotes(fname)
        if (fname == 'q' .or. fname == 'Q') then
            iostat = 0; return
        end if
        if (len_trim(fname) == 0) then
            write(*, '(a)') '  No file specified.'; return
        end if

        ! ── Prompt for .castep file ──
        write(*, '(a)', advance='no') '  Enter companion .castep file path: '
        read(*, '(a)', iostat=ios) castep_path
        if (ios /= 0) return
        castep_path = adjustl(castep_path); call strip_quotes(castep_path)
        if (castep_path == 'q' .or. castep_path == 'Q') then
            iostat = 0; return
        end if

        ! ── Parse phonon eigenvectors ──
        call parse_phonon_eigenvectors(trim(fname), modes_data, ios, msg)
        if (ios /= 0) then
            write(*, '(a,a)') '  Error parsing .phonon: ', trim(msg)
            call free_phonon_modes_data(modes_data); return
        end if
        write(*, '(a,i0,a,i0)') '  Parsed ', modes_data%n_ions, ' ions, ', &
            modes_data%n_branches, ' branches.'

        ! ── Parse Born charges from .castep ──
        call parse_castep_born_charges(trim(castep_path), modes_data, ios, msg)
        if (ios /= 0) then
            write(*, '(a,a)') '  Warning: ', trim(msg)
            write(*, '(a)') '  Displacements will be shown without IR decomposition.'
        end if

        ! ── Compute mode decomposition ──
        call compute_mode_decomposition(modes_data, ios, msg)
        if (ios /= 0) then
            write(*, '(a,a)') '  Error in decomposition: ', trim(msg)
            call free_phonon_modes_data(modes_data); return
        end if

        ! ── Find mode with highest IR intensity ──
        best_mode = 1
        max_ir = 0.0_dp
        do ios = 1, modes_data%n_branches
            if (modes_data%modes(ios)%ir_intensity > max_ir) then
                max_ir = modes_data%modes(ios)%ir_intensity
                best_mode = ios
            end if
        end do

        ! ── Print mode list summary ──
        write(*, '(a)') ''
        write(*, '(a)') '  ============================================================'
        write(*, '(a)') '         Phonon Mode Summary'
        write(*, '(a,i0,a,i0,a)') '         ', modes_data%n_ions, ' ions, ', &
            modes_data%n_branches, ' branches'
        write(*, '(a)') '  ============================================================'
        write(*, '(a)') ''
        write(*, '(a)') '  Mode     Freq/cm⁻¹      IR Int.      |p_mode|'
        write(*, '(a)') '  ────     ─────────      ───────      ────────'
        do ios = 1, modes_data%n_branches
            if (ios == best_mode) then
                write(*, '(a,i4,a,f13.2,a,f11.4,a,f11.4,a)') &
                    ' ▶', ios, '   ', modes_data%modes(ios)%frequency, '   ', &
                    modes_data%modes(ios)%ir_intensity, '   ', &
                    modes_data%modes(ios)%mode_charge_norm, &
                    '  ◀ highest IR'
            else
                write(*, '(a,i4,a,f13.2,a,f11.4,a,f11.4)') &
                    '  ', ios, '   ', modes_data%modes(ios)%frequency, '   ', &
                    modes_data%modes(ios)%ir_intensity, '   ', &
                    modes_data%modes(ios)%mode_charge_norm
            end if
        end do
        write(*, '(a)') ''
        write(*, '(a,i0,a,f12.4,a)') '  ▶ Mode ', best_mode, &
            ' has highest IR intensity (', max_ir, ')'
        write(*, '(a)') ''

        ! ── Mode selection loop ──
        json_path = '_crystal_viewer_temp.json'
        do
            write(*, '(a,i0,a)') '  Enter mode number (1-', modes_data%n_branches, &
                ', H=highest, Q=return): '
            write(*, '(a)', advance='no') '  > '
            read(*, '(a)', iostat=ios) input
            if (ios /= 0) exit
            input = adjustl(input)

            if (input(1:1) == 'q' .or. input(1:1) == 'Q') then
                iostat = 0; exit
            end if
            if (input(1:1) == 'h' .or. input(1:1) == 'H') then
                mode_idx = best_mode
            else
                read(input, *, iostat=ios) mode_idx
            end if

            if (ios /= 0 .or. mode_idx < 1 .or. mode_idx > modes_data%n_branches) then
                write(*, '(a,i0,a)') '  Invalid. Enter 1-', modes_data%n_branches, &
                    ', H, or Q.'
                cycle
            end if

            write(*, '(a,i0,a,f10.4,a,f10.4)') '  Mode ', mode_idx, &
                ': freq=', modes_data%modes(mode_idx)%frequency, ' cm-1, IR=', &
                modes_data%modes(mode_idx)%ir_intensity

            ! Write JSON and launch viewer
            call write_crystal_json_modes(modes_data, mode_idx, json_path, 2.0_dp, ios)
            if (ios /= 0) then
                write(*, '(a)') '  Error writing JSON file.'; cycle
            end if
            call launch_viewer(json_path)
        end do

        ! Cleanup
        call free_phonon_modes_data(modes_data)
        block
            integer :: tmp_unit, tmp_ios
            open(newunit=tmp_unit, file=trim(json_path), status='old', iostat=tmp_ios)
            if (tmp_ios == 0) close(tmp_unit, status='delete')
        end block
    end subroutine handle_phonon_modes_menu


    subroutine handle_thermo_menu(iostat)
        !! Thermodynamics: E(T), S(T), F(T), Cv(T) from .phonon file
        integer, intent(out) :: iostat
        type(phonon_dos_data_t) :: phdos
        type(thermo_data_t) :: thermo
        character(len=MAX_LINE_LEN) :: fname, input, tmp_str
        character(len=512) :: csv_file
        integer :: ios, n_pts, unit, i
        real(dp) :: t_min, t_max
        character(len=MAX_LINE_LEN), save :: last_phonon_path = ''

        iostat = 0

        ! ── File input ──
        write(*, '(a)', advance='no') '  Enter .phonon file path'
        if (len_trim(last_phonon_path) > 0) &
            write(*, '(a)', advance='no') ' [' // trim(last_phonon_path) // ']'
        write(*, '(a)') ': '
        read(*, '(a)', iostat=ios) fname
        if (ios /= 0) return
        fname = adjustl(fname); call strip_quotes(fname)
        if (fname == 'q' .or. fname == 'Q') then
            iostat = 0; return
        end if
        if (len_trim(fname) == 0 .and. len_trim(last_phonon_path) > 0) then
            fname = last_phonon_path
        end if
        if (len_trim(fname) == 0) then
            write(*, '(a)') '  No file specified.'; return
        end if
        last_phonon_path = trim(fname)

        call parse_phonon_file(trim(fname), phdos, ios)
        if (ios /= 0) then
            write(*, '(a)') '  Error parsing .phonon file.'; return
        end if
        write(*, '(a,i0,a,i0,a,i0,a)') '  Loaded ', phdos%n_ions, ' ions, ', &
            phdos%n_branches, ' branches, ', phdos%n_qpoints, ' q-points.'

        ! ── Temperature range ──
        t_min = 0.0_dp
        t_max = 1000.0_dp
        n_pts = 200

        write(*, '(a,f8.1,a)', advance='no') '  T_min (K) [', t_min, ']: '
        read(*, '(a)', iostat=ios) input
        if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
            if (adjustl(trim(input)) == 'q' .or. adjustl(trim(input)) == 'Q') then
                call free_phonon_dos_data(phdos); iostat = 0; return
            end if
            read(input, *, iostat=ios) t_min
        end if

        write(*, '(a,f8.1,a)', advance='no') '  T_max (K) [', t_max, ']: '
        read(*, '(a)', iostat=ios) input
        if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
            if (adjustl(trim(input)) == 'q' .or. adjustl(trim(input)) == 'Q') then
                call free_phonon_dos_data(phdos); iostat = 0; return
            end if
            read(input, *, iostat=ios) t_max
        end if
        if (t_max <= t_min) t_max = t_min + 100.0_dp

        write(*, '(a,i0,a)', advance='no') '  Number of points [', n_pts, ']: '
        read(*, '(a)', iostat=ios) input
        if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
            if (adjustl(trim(input)) == 'q' .or. adjustl(trim(input)) == 'Q') then
                call free_phonon_dos_data(phdos); iostat = 0; return
            end if
            read(input, *, iostat=ios) n_pts
        end if
        n_pts = max(2, min(1000, n_pts))

        ! ── Compute ──
        write(*, '(a)') '  Computing thermodynamics...'
        call compute_thermodynamics(phdos, t_min, t_max, n_pts, thermo, ios, tmp_str)
        if (ios /= 0) then
            write(*, '(a,a)') '  Error: ', trim(tmp_str)
            call free_phonon_dos_data(phdos); return
        end if

        write(*, '(a,es14.6,a)') '  Zero-point energy: ', thermo%zpe, ' eV'

        ! ── Output mode ──
        write(*, '(a)') '  Output mode:'
        write(*, '(a)') '    1. Terminal ASCII plot'
        write(*, '(a)') '    2. Export CSV'
        write(*, '(a)', advance='no') '    Enter choice [1]: '
        read(*, '(a)', iostat=ios) input
        if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
            if (adjustl(trim(input)) == 'q' .or. adjustl(trim(input)) == 'Q') then
                call free_phonon_dos_data(phdos); call free_thermo_data(thermo); iostat = 0; return
            end if
            read(input, *, iostat=ios) n_pts
        else
            n_pts = 1
        end if

        if (n_pts == 2) then
            ! CSV export
            write(*, '(a)', advance='no') '  Output CSV file (without .csv): '
            read(*, '(a)', iostat=ios) csv_file
            if (ios /= 0) then
                call free_phonon_dos_data(phdos); call free_thermo_data(thermo); return
            end if
            csv_file = adjustl(csv_file); call strip_quotes(csv_file)
            if (len_trim(csv_file) == 0) then
                call free_phonon_dos_data(phdos); call free_thermo_data(thermo); return
            end if
            call ensure_ext(csv_file, '.csv')
            open(newunit=unit, file=trim(csv_file), status='replace', action='write', iostat=ios)
            if (ios /= 0) then
                write(*, '(a)') '  Cannot write CSV file.'
                call free_phonon_dos_data(phdos); call free_thermo_data(thermo); return
            end if
            write(unit, '(a)') '# T(K),E_vib(eV),F_vib(eV),TS(eV),S(eV/K),Cv(eV/K)'
            do i = 1, thermo%n_temps
                write(unit, '(es12.4,5(a,es14.6))') &
                    thermo%temps(i), ',', &
                    thermo%energy(i) - thermo%zpe, ',', &
                    thermo%free_e(i) - thermo%zpe, ',', &
                    thermo%temps(i) * thermo%entropy(i), ',', &
                    thermo%entropy(i), ',', &
                    thermo%heat_cap(i)
            end do
            close(unit)
            write(*, '(a)') '  Written ' // trim(csv_file)
        else
            ! ASCII plot
            call run_thermo_navigator(thermo)
        end if

        call free_phonon_dos_data(phdos)
        call free_thermo_data(thermo)
    end subroutine handle_thermo_menu


    subroutine run_thermo_navigator(thermo)
        !! Interactive ASCII plot: E, S, F, Cv vs T using unified plot_dos_ascii
        type(thermo_data_t), intent(in) :: thermo
        integer, parameter :: N_CURVES = 2
        integer :: cur_plot, tw, th, ios, nspin
        real(dp) :: x_center, half_range, x_center0, half_range0
        real(dp) :: y_center(N_CURVES), y_half(N_CURVES)
        real(dp) :: y_center0(N_CURVES), y_half0(N_CURVES)
        real(dp), allocatable :: curve_data(:,:), disp_data(:,:)
        character(len=32) :: titles(N_CURVES)
        character(len=1) :: ch
        character(len=3) :: arrow

        cur_plot = 1
        titles(1) = 'E_vib, F_vib & TS vs T'
        titles(2) = 'Heat Capacity vs T'

        ! Build data: E_vib, F_vib, TS, Cv
        allocate(curve_data(thermo%n_temps, 4))
        curve_data(:,1) = thermo%energy - thermo%zpe   ! E_vib (eV)
        curve_data(:,2) = thermo%free_e - thermo%zpe   ! F_vib (eV)
        curve_data(:,3) = thermo%temps * thermo%entropy ! TS (eV)
        curve_data(:,4) = thermo%heat_cap              ! Cv (eV/K)

        ! Display 1 auto-range: covers E_vib, F_vib, TS
        y_half(1) = (max(maxval(curve_data(:,1)), maxval(curve_data(:,2)), &
                      maxval(curve_data(:,3))) &
                  - min(minval(curve_data(:,1)), minval(curve_data(:,2)), &
                      minval(curve_data(:,3)))) * 0.575_dp
        y_center(1) = (max(maxval(curve_data(:,1)), maxval(curve_data(:,2)), &
                        maxval(curve_data(:,3))) &
                    + min(minval(curve_data(:,1)), minval(curve_data(:,2)), &
                        minval(curve_data(:,3)))) * 0.5_dp
        if (y_half(1) < 1.0e-12_dp) y_half(1) = 1.0_dp
        y_center0(1) = y_center(1); y_half0(1) = y_half(1)

        ! Display 2 auto-range: Cv only
        y_half(2) = (maxval(curve_data(:,4)) - minval(curve_data(:,4))) * 0.575_dp
        y_center(2) = (maxval(curve_data(:,4)) + minval(curve_data(:,4))) * 0.5_dp
        if (y_half(2) < 1.0e-12_dp) y_half(2) = 1.0_dp
        y_center0(2) = y_center(2); y_half0(2) = y_half(2)

        x_center = (thermo%temps(1) + thermo%temps(thermo%n_temps)) * 0.5_dp
        half_range = (thermo%temps(thermo%n_temps) - thermo%temps(1)) * 0.5_dp
        x_center0 = x_center
        half_range0 = half_range

        call enter_raw_mode
        do
            call get_term_size(tw, th)
            write(*, '(a)', advance='no') achar(27) // '[2J' // achar(27) // '[H'
            write(*, '(a,es14.6,a)') '  ZPE = ', thermo%zpe, ' eV  (E_vib, F_vib exclude ZPE)'
            write(*, '(a)') ''

            if (cur_plot == 1) then
                ! Overlay E_vib + F_vib + TS (all eV)
                nspin = 3
                allocate(disp_data(thermo%n_temps, 3))
                disp_data(:,1) = curve_data(:,1)  ! E_vib
                disp_data(:,2) = curve_data(:,2)  ! F_vib
                disp_data(:,3) = curve_data(:,3)  ! TS
                call plot_dos_ascii(thermo%temps, disp_data, nspin, &
                    0.0_dp, 0.0_dp, tw, th, &
                    y_center_in=y_center(1), y_half_in=y_half(1), &
                    e_center_in=x_center, half_range_in=half_range, &
                    xlabel='Temperature', xunit='K', hide_fermi=.true., &
                    title=titles(1))
                deallocate(disp_data)
            else
                ! Cv only
                nspin = 1
                allocate(disp_data(thermo%n_temps, 1))
                disp_data(:,1) = curve_data(:,4)
                call plot_dos_ascii(thermo%temps, disp_data, nspin, &
                    0.0_dp, 0.0_dp, tw, th, &
                    y_center_in=y_center(2), y_half_in=y_half(2), &
                    e_center_in=x_center, half_range_in=half_range, &
                    xlabel='Temperature', xunit='K', hide_fermi=.true., &
                    title=titles(2))
                deallocate(disp_data)
            end if

            write(*, '(a)') '  [1=E+F  2=Cv  arrows:pan  +/-:zoom  R:reset  Q:quit]'

            read(*, '(a)', advance='no', iostat=ios) ch
            if (ios /= 0) exit
            if (iachar(ch) == 27) then
                read(*, '(a)', advance='no', iostat=ios) arrow(1:1)
                if (ios /= 0) exit
                read(*, '(a)', advance='no', iostat=ios) arrow(2:2)
                if (ios /= 0) exit
                if (arrow(1:2) == '[A') then
                    y_center(cur_plot) = y_center(cur_plot) + y_half(cur_plot) * 0.3_dp
                else if (arrow(1:2) == '[B') then
                    y_center(cur_plot) = y_center(cur_plot) - y_half(cur_plot) * 0.3_dp
                else if (arrow(1:2) == '[C') then
                    x_center = x_center + half_range * 0.3_dp
                else if (arrow(1:2) == '[D') then
                    x_center = x_center - half_range * 0.3_dp
                end if
                cycle
            end if
            select case (ch)
            case ('1'); cur_plot = 1
            case ('2'); cur_plot = 2
            case ('+', '=')
                half_range = half_range * 0.5_dp
                y_half(cur_plot) = y_half(cur_plot) * 0.5_dp
            case ('-')
                half_range = half_range * 2.0_dp
                y_half(cur_plot) = y_half(cur_plot) * 2.0_dp
            case ('r', 'R')
                x_center = x_center0; half_range = half_range0
                y_center(cur_plot) = y_center0(cur_plot)
                y_half(cur_plot) = y_half0(cur_plot)
            case ('q', 'Q'); exit
            end select
        end do
        call leave_raw_mode
        deallocate(curve_data)
    end subroutine run_thermo_navigator


    subroutine free_cif_data(data)
        !! Deallocate cif_data_t atom arrays
        type(cif_data_t), intent(inout) :: data
        if (allocated(data%atoms)) deallocate(data%atoms)
        data%n_atoms = 0
    end subroutine free_cif_data

    ! ═══════════════════════════════════════════════════════════════
    !  PES Scan sub-menu
    ! ═══════════════════════════════════════════════════════════════

    subroutine handle_pes_scan_menu(iostat)
        integer, intent(out) :: iostat
        integer :: ios
        character(len=256) :: input
        logical :: has_saved  ! tracks if first sub-item has been used

        iostat = 0
        has_saved = .false.

        do
            write(*, '(a)') ''
            write(*, '(a)') '  =================================='
            write(*, '(a)') '            PES Scan'
            write(*, '(a)') '  =================================='
            write(*, '(a)') '  1. Generate 2D Scan Input Files'
            write(*, '(a)') '  2. Collect 2D Results & Visualize'
            write(*, '(a)') '  3. Generate 3D Scan Input Files'
            write(*, '(a)') '  4. Collect 3D Results & Visualize'
            write(*, '(a)') '  Q. Back'
            write(*, '(a)', advance='no') '  Select option: '
            read(*, '(a)', iostat=ios) input
            if (ios /= 0) then; iostat = 0; exit; end if
            input = adjustl(input)

            select case (trim(input))
            case ('1')
                call handle_pes_generate(ios)
                if (ios == 0) has_saved = .true.
            case ('2')
                call handle_pes_collect(ios)
            case ('3')
                call handle_pes3d_generate(ios)
                if (ios == 0) has_saved = .true.
            case ('4')
                call handle_pes3d_collect(ios)
            case ('q', 'Q')
                exit
            case default
                write(*, '(a)') '  Invalid option.'
            end select
        end do
    end subroutine handle_pes_scan_menu


    subroutine handle_pes_generate(iostat)
        !! PES sub-menu 1: Generate scan input files
        integer, intent(out) :: iostat

        type(cif_data_t) :: cif
        type(castep_config_t) :: cfg
        type(pes_grid_t) :: grid
        character(len=MAX_LINE_LEN) :: file_path, scan_dir, sub_dir, stem, input
        character(len=MAX_LINE_LEN) :: cell_path, param_path
        real(dp), allocatable :: frac_points(:,:)
        integer :: ios, i, n_total, mi, n_pts(2)
        real(dp) :: fx, fy, range_vals(4)
        character(len=8) :: plane_choice, mode_choice, mobile_elem, plane_str

        iostat = 0

        ! ── Load structure ──
        write(*, '(a)', advance='no') '  Enter structure file (.cif/.pdb/.cell): '
        read(*, '(a)', iostat=ios) file_path
        if (ios /= 0) return
        file_path = adjustl(file_path); call strip_quotes(file_path)
        if (file_path == 'q' .or. file_path == 'Q') return
        select case (trim(get_ext_lower(file_path)))
        case ('cif')
            call parse_cif_inline(trim(file_path), cif, ios)
        case ('pdb')
            call parse_pdb_inline(trim(file_path), cif, ios)
        case ('cell')
            call parse_cell_inline(trim(file_path), cif, ios)
        case default
            write(*, '(a)') '  Unsupported file format. Use .cif, .pdb, or .cell.'
            iostat = 1; return
        end select
        if (ios /= 0) then
            write(*, '(a)') '  Error parsing file.'; return
        end if
        call expand_cif_symmetry(cif, ios)

        ! ── Populate castep_config_t ──
        call default_config(cfg)
        stem = get_file_stem(file_path)
        cfg%cif_file_path = trim(file_path)
        cfg%cell_length(1) = cif%a
        cfg%cell_length(2) = cif%b
        cfg%cell_length(3) = cif%c
        cfg%cell_angle(1)  = cif%alpha
        cfg%cell_angle(2)  = cif%beta
        cfg%cell_angle(3)  = cif%gamma
        cfg%num_atoms = cif%n_atoms
        cfg%cartesian_coords = .not. cif%positions_fractional
        allocate(cfg%atom_type(cfg%num_atoms))
        allocate(cfg%atom_x(cfg%num_atoms))
        allocate(cfg%atom_y(cfg%num_atoms))
        allocate(cfg%atom_z(cfg%num_atoms))
        do i = 1, cif%n_atoms
            cfg%atom_type(i) = trim(clean_element_symbol(cif%atoms(i)%element))
            cfg%atom_x(i)    = cif%atoms(i)%x
            cfg%atom_y(i)    = cif%atoms(i)%y
            cfg%atom_z(i)    = cif%atoms(i)%z
        end do
        call free_cif_data(cif)

        ! ── Select mobile atom ──
        write(*, '(a)') '  Atom list:'
        do i = 1, cfg%num_atoms
            if (cfg%cartesian_coords) then
                write(*, '(i4, a, a8, a, 3f10.6)') i, '. ', trim(cfg%atom_type(i)), '  cart:', &
                    cfg%atom_x(i), cfg%atom_y(i), cfg%atom_z(i)
            else
                write(*, '(i4, a, a8, a, 3f10.6)') i, '. ', trim(cfg%atom_type(i)), '  frac:', &
                    cfg%atom_x(i), cfg%atom_y(i), cfg%atom_z(i)
            end if
        end do
        write(*, '(a,i0,a)', advance='no') '  Select mobile atom (1-', cfg%num_atoms, '): '
        read(*, *, iostat=ios) mi
        if (ios /= 0 .or. mi < 1 .or. mi > cfg%num_atoms) then
            write(*, '(a)') '  Invalid selection.'; return
        end if
        grid%mobile_atom_idx = mi

        ! ── Select plane ──
        write(*, '(a)') '  Scan plane: 1=xy  2=xz  3=yz'
        write(*, '(a)', advance='no') '  Choice: '
        read(*, '(a)', iostat=ios) plane_choice
        if (plane_choice(1:1) == '1') then; grid%plane_axis = [1, 2]
        else if (plane_choice(1:1) == '2') then; grid%plane_axis = [1, 3]
        else if (plane_choice(1:1) == '3') then; grid%plane_axis = [2, 3]
        else; write(*, '(a)') '  Invalid.'; return; end if

        ! ── Grid parameters ──
        write(*, '(a)', advance='no') '  Grid Nx (default 5): '
        read(*, '(a)', iostat=ios) input
        n_pts(1) = 5; if (len_trim(input) > 0) read(input, *, iostat=ios) n_pts(1)
        write(*, '(a)', advance='no') '  Grid Ny (default 5): '
        read(*, '(a)', iostat=ios) input
        n_pts(2) = 5; if (len_trim(input) > 0) read(input, *, iostat=ios) n_pts(2)
        grid%n_points = n_pts

        write(*, '(a)', advance='no') '  fx range (min max, default 0 1): '
        read(*, '(a)', iostat=ios) input
        range_vals(1) = 0.0_dp; range_vals(2) = 1.0_dp
        if (len_trim(input) > 0) read(input, *, iostat=ios) range_vals(1), range_vals(2)
        grid%frac_range(1,1) = range_vals(1); grid%frac_range(1,2) = range_vals(2)

        write(*, '(a)', advance='no') '  fy range (min max, default 0 1): '
        read(*, '(a)', iostat=ios) input
        range_vals(3) = 0.0_dp; range_vals(4) = 1.0_dp
        if (len_trim(input) > 0) read(input, *, iostat=ios) range_vals(3), range_vals(4)
        grid%frac_range(2,1) = range_vals(3); grid%frac_range(2,2) = range_vals(4)

        ! ── Scan mode ──
        write(*, '(a)') '  Mode: 1=Single-point (SP)  2=Constrained relax (RELAX)'
        write(*, '(a)', advance='no') '  Choice (default SP): '
        read(*, '(a)', iostat=ios) mode_choice
        grid%scan_mode = 'SP'
        if (mode_choice(1:1) == '2') grid%scan_mode = 'RELAX'

        if (grid%scan_mode == 'RELAX') then
            cfg%task_type = TASK_GEOMETRY_OPT
            cfg%cell_opt_mode = 'FIX_CELL'
        else
            cfg%task_type = TASK_ENERGY
        end if

        ! Set PES constraints: mobile atom is fixed in the direction
        ! perpendicular to the scan plane, free to move in the plane.
        ! SP mode: mobile atom fully fixed ([1,1,1]).
        ! RELAX mode: mobile atom fixed only out-of-plane (e.g. xy→[0,0,1]).
        cfg%pes_mobile_idx = mi
        cfg%pes_fix_axes = [1, 1, 1]
        if (grid%scan_mode == 'RELAX') then
            cfg%pes_fix_axes = [0, 0, 0]
            ! The third axis (not in the scan plane) is the fixed direction
            cfg%pes_fix_axes(6 - grid%plane_axis(1) - grid%plane_axis(2)) = 1
        end if

        ! ── PreCASTEP parameter configuration ──
        write(*, '(a)') ''
        write(*, '(a)') '  ── Configure CASTEP parameters ──'
        write(*, '(a)') '  (Select 0. Generate when ready to proceed)'
        cfg%cif_file_path = trim(file_path)
        call run_main_menu(cfg, ios)
        if (ios == IO_USER_QUIT) return

        ! ── Generate grid points ──
        call generate_pes_grid_points(grid, frac_points, n_total, ios)
        if (ios /= 0) then
            write(*, '(a)') '  Error generating grid.'; return
        end if

        ! ── Create output directory ──
        scan_dir = 'pes_scan/' // trim(stem) // '_' // trim(grid%scan_mode)
        sub_dir = 'mkdir -p "' // trim(scan_dir) // '"'
        call execute_command_line(trim(sub_dir), exitstat=ios)

        ! ── Generate files ──
        write(*, '(a,i0,a)') '  Generating ', n_total, ' grid points...'
        do i = 1, n_total
            fx = frac_points(i, 1)
            fy = frac_points(i, 2)

            ! Set mobile atom to grid point fractional coords
            select case (grid%plane_axis(1))
            case (1); cfg%atom_x(mi) = fx
            case (2); cfg%atom_y(mi) = fx
            case (3); cfg%atom_z(mi) = fx
            end select
            select case (grid%plane_axis(2))
            case (1); cfg%atom_x(mi) = fy
            case (2); cfg%atom_y(mi) = fy
            case (3); cfg%atom_z(mi) = fy
            end select

            write(sub_dir, '(a, a, i3.3, a, i3.3)') trim(scan_dir), '/grid_', &
                modulo(i-1, grid%n_points(1)) + 1, '_', &
                (i-1) / grid%n_points(1) + 1

            call execute_command_line('mkdir -p "' // trim(sub_dir) // '"', exitstat=ios)

            cell_path  = trim(sub_dir) // '/scan.cell'
            param_path = trim(sub_dir) // '/scan.param'
            call write_cell_file(cell_path, cfg, ios)
            call write_param_file(param_path, cfg, ios)
        end do

        ! ── Write metadata JSON ──
        call write_pes_metadata_json(trim(scan_dir)//'/pes_metadata.json', grid, cfg, ios)

        ! ── Cleanup ──
        ! Save element and plane label before deallocation
        mobile_elem = trim(cfg%atom_type(mi))
        if (grid%plane_axis(1) == 1 .and. grid%plane_axis(2) == 2) then
            plane_str = 'XY'
        else if (grid%plane_axis(1) == 1 .and. grid%plane_axis(2) == 3) then
            plane_str = 'XZ'
        else
            plane_str = 'YZ'
        end if

        deallocate(frac_points)
        deallocate(cfg%atom_type, cfg%atom_x, cfg%atom_y, cfg%atom_z)

        write(*, '(a)') ''
        write(*, '(a)') '  ========================================'
        write(*, '(a)') '       PES scan files generated!'
        write(*, '(a)') '  ========================================'
        write(*, '(a,a)')   '  Structure:    ', trim(stem)
        write(*, '(a,a)')   '  Scan mode:    ', trim(grid%scan_mode)
        write(*, '(a,i0,a,a,a)') '  Mobile atom:  #', mi, ' (', trim(mobile_elem), ')'
        write(*, '(a,i0,a,i0)') '  Grid:         ', grid%n_points(1), ' x ', grid%n_points(2)
        write(*, '(a,i0)')      '  Total points: ', n_total
        write(*, '(a,a)')   '  Plane:        ', trim(plane_str)
        write(*, '(a,f8.4,a,f8.4,a)') '  fx range:     [', grid%frac_range(1,1), ', ', grid%frac_range(1,2), ']'
        write(*, '(a,f8.4,a,f8.4,a)') '  fy range:     [', grid%frac_range(2,1), ', ', grid%frac_range(2,2), ']'
        write(*, '(a,a)')   '  Directory:    ', trim(scan_dir)
        write(*, '(a)') ''
    end subroutine handle_pes_generate


    subroutine handle_pes_collect(iostat)
        !! PES sub-menu 2: Collect .castep results and launch viewer
        integer, intent(out) :: iostat

        character(len=MAX_LINE_LEN) :: scan_dir, json_path, input
        integer :: ios
        logical :: exists

        iostat = 0

        write(*, '(a)') '  Enter scan directory (e.g. pes_scan/NaCl):'
        write(*, '(a)', advance='no') '  Path: '
        read(*, '(a)', iostat=ios) scan_dir
        if (ios /= 0) return
        scan_dir = adjustl(scan_dir); call strip_quotes(scan_dir)
        if (scan_dir == 'q' .or. scan_dir == 'Q') return

        inquire(file=trim(scan_dir)//'/pes_metadata.json', exist=exists)
        if (.not. exists) then
            write(*, '(a)') '  pes_metadata.json not found in directory.'
            return
        end if

        ! ── Collect energies ──
        call collect_pes_energies(trim(scan_dir), ios)
        if (ios /= 0) then
            write(*, '(a)') '  Error collecting energies.'
            return
        end if

        ! ── Launch viewer? ──
        write(*, '(a)', advance='no') '  Launch 3D Viewer? (y/n): '
        read(*, '(a)', iostat=ios) input
        if (ios == 0 .and. (input(1:1) == 'y' .or. input(1:1) == 'Y')) then
            json_path = trim(scan_dir) // '/pes_metadata.json'
            call launch_viewer(json_path)
        end if
    end subroutine handle_pes_collect


    ! ── 3D PES Scan sub-routines ──

    subroutine handle_pes3d_generate(iostat)
        !! PES 3D sub-menu: Generate 3D voxel grid scan input files.
        !! Supports symmetry-aware scanning: with space-group info, only the
        !! asymmetric sub-volume around a reference atom is scanned; without
        !! symmetry, full [0,1)^3 cell scan.
        integer, intent(out) :: iostat

        type(cif_data_t) :: cif
        type(castep_config_t) :: cfg
        type(pes3d_grid_t) :: grid
        character(len=MAX_LINE_LEN) :: file_path, scan_dir, sub_dir, stem, input
        character(len=MAX_LINE_LEN) :: cell_path, param_path
        real(dp), allocatable :: frac_points(:,:)
        real(dp), allocatable :: dummy_energies(:)
        integer :: ios, i, n_total, mi, n_pts(3)
        real(dp) :: fx, fy, fz, range_vals(6)
        character(len=8) :: mode_choice, sym_choice
        character(len=8) :: mobile_elem

        iostat = 0

        ! ── Load structure ──
        write(*, '(a)', advance='no') '  Enter structure file (.cif/.pdb/.cell): '
        read(*, '(a)', iostat=ios) file_path
        if (ios /= 0) return
        file_path = adjustl(file_path); call strip_quotes(file_path)
        if (file_path == 'q' .or. file_path == 'Q') return
        select case (trim(get_ext_lower(file_path)))
        case ('cif')
            call parse_cif_inline(trim(file_path), cif, ios)
        case ('pdb')
            call parse_pdb_inline(trim(file_path), cif, ios)
        case ('cell')
            call parse_cell_inline(trim(file_path), cif, ios)
        case default
            write(*, '(a)') '  Unsupported file format. Use .cif, .pdb, or .cell.'
            iostat = 1; return
        end select
        if (ios /= 0) then
            write(*, '(a)') '  Error parsing file.'; return
        end if

        ! ── Symmetry choice (only for CIF with sym_ops) ──
        grid%use_symmetry = .false.
        if (cif%n_symops > 1) then
            write(*, '(a)') ''
            write(*, '(a,i0,a)') '  Structure has ', cif%n_symops, ' symmetry operations.'
            write(*, '(a)') '  How to handle symmetry?'
            write(*, '(a)') '    1. Use symmetry — scan only asymmetric sub-volume,'
            write(*, '(a)') '       then expand to full cell (recommended)'
            write(*, '(a)') '    2. Treat as P1 — full-cell scan, no symmetry'
            write(*, '(a)', advance='no') '  Choice (1/2, default 1): '
            read(*, '(a)', iostat=ios) sym_choice
            if (sym_choice(1:1) == '2') then
                grid%use_symmetry = .false.
            else
                grid%use_symmetry = .true.
            end if
        end if

        ! ── Expand symmetry only for atom list display ──
        ! Keep cif%sym_ops intact for later JSON output
        call expand_cif_symmetry(cif, ios)

        ! ── Populate castep_config_t ──
        call default_config(cfg)
        stem = get_file_stem(file_path)
        cfg%cif_file_path = trim(file_path)
        cfg%cell_length(1) = cif%a
        cfg%cell_length(2) = cif%b
        cfg%cell_length(3) = cif%c
        cfg%cell_angle(1)  = cif%alpha
        cfg%cell_angle(2)  = cif%beta
        cfg%cell_angle(3)  = cif%gamma
        cfg%num_atoms = cif%n_atoms
        cfg%cartesian_coords = .not. cif%positions_fractional
        allocate(cfg%atom_type(cfg%num_atoms))
        allocate(cfg%atom_x(cfg%num_atoms))
        allocate(cfg%atom_y(cfg%num_atoms))
        allocate(cfg%atom_z(cfg%num_atoms))
        do i = 1, cif%n_atoms
            cfg%atom_type(i) = trim(clean_element_symbol(cif%atoms(i)%element))
            cfg%atom_x(i)    = cif%atoms(i)%x
            cfg%atom_y(i)    = cif%atoms(i)%y
            cfg%atom_z(i)    = cif%atoms(i)%z
        end do

        ! ── Select reference atom ──
        write(*, '(a)') '  Atom list (expanded cell):'
        do i = 1, cfg%num_atoms
            if (cfg%cartesian_coords) then
                write(*, '(i4, a, a8, a, 3f10.6)') i, '. ', trim(cfg%atom_type(i)), '  cart:', &
                    cfg%atom_x(i), cfg%atom_y(i), cfg%atom_z(i)
            else
                write(*, '(i4, a, a8, a, 3f10.6)') i, '. ', trim(cfg%atom_type(i)), '  frac:', &
                    cfg%atom_x(i), cfg%atom_y(i), cfg%atom_z(i)
            end if
        end do
        write(*, '(a,i0,a)', advance='no') '  Select reference atom (1-', cfg%num_atoms, '): '
        read(*, *, iostat=ios) mi
        if (ios /= 0 .or. mi < 1 .or. mi > cfg%num_atoms) then
            write(*, '(a)') '  Invalid selection.'; return
        end if
        grid%ref_atom_idx = mi
        grid%ref_frac = [cfg%atom_x(mi), cfg%atom_y(mi), cfg%atom_z(mi)]

        ! ── Grid parameters ──
        if (grid%use_symmetry) then
            ! Symmetry mode: auto-compute sub-volume bounds
            call compute_local_grid_bounds(mi, cif%atoms, cif%n_atoms, &
                cif%sym_ops, cif%n_symops, grid%half_dist, ios)
            if (ios /= 0) then
                write(*, '(a)') '  Warning: could not compute symmetry bounds, using defaults.'
                grid%half_dist = 0.3_dp
            end if

            write(*, '(a)') ''
            write(*, '(a)') '  ── Asymmetric Sub-Volume (auto-computed) ──'
            write(*, '(a,3f10.6)') '  Reference atom frac: ', grid%ref_frac
            write(*, '(a,3f10.6)') '  Half-distance:       ', grid%half_dist
            grid%frac_range(1,1) = max(0.0_dp, grid%ref_frac(1) - grid%half_dist(1))
            grid%frac_range(1,2) = min(1.0_dp, grid%ref_frac(1) + grid%half_dist(1))
            grid%frac_range(2,1) = max(0.0_dp, grid%ref_frac(2) - grid%half_dist(2))
            grid%frac_range(2,2) = min(1.0_dp, grid%ref_frac(2) + grid%half_dist(2))
            grid%frac_range(3,1) = max(0.0_dp, grid%ref_frac(3) - grid%half_dist(3))
            grid%frac_range(3,2) = min(1.0_dp, grid%ref_frac(3) + grid%half_dist(3))
            write(*, '(a, 2f10.6)') '  fx range: ', grid%frac_range(1,:)
            write(*, '(a, 2f10.6)') '  fy range: ', grid%frac_range(2,:)
            write(*, '(a, 2f10.6)') '  fz range: ', grid%frac_range(3,:)
        else
            ! No symmetry: user inputs grid range manually
            grid%frac_range(1,1) = 0.0_dp; grid%frac_range(1,2) = 1.0_dp
            grid%frac_range(2,1) = 0.0_dp; grid%frac_range(2,2) = 1.0_dp
            grid%frac_range(3,1) = 0.0_dp; grid%frac_range(3,2) = 1.0_dp
        end if

        write(*, '(a)', advance='no') '  Grid Nx (default 5): '
        read(*, '(a)', iostat=ios) input
        n_pts(1) = 5; if (len_trim(input) > 0) read(input, *, iostat=ios) n_pts(1)
        write(*, '(a)', advance='no') '  Grid Ny (default 5): '
        read(*, '(a)', iostat=ios) input
        n_pts(2) = 5; if (len_trim(input) > 0) read(input, *, iostat=ios) n_pts(2)
        write(*, '(a)', advance='no') '  Grid Nz (default 5): '
        read(*, '(a)', iostat=ios) input
        n_pts(3) = 5; if (len_trim(input) > 0) read(input, *, iostat=ios) n_pts(3)
        grid%n_points = n_pts

        if (.not. grid%use_symmetry) then
            write(*, '(a)', advance='no') '  fx range (min max, default 0 1): '
            read(*, '(a)', iostat=ios) input
            range_vals(1) = 0.0_dp; range_vals(2) = 1.0_dp
            if (len_trim(input) > 0) read(input, *, iostat=ios) range_vals(1), range_vals(2)
            grid%frac_range(1,1) = range_vals(1); grid%frac_range(1,2) = range_vals(2)

            write(*, '(a)', advance='no') '  fy range (min max, default 0 1): '
            read(*, '(a)', iostat=ios) input
            range_vals(3) = 0.0_dp; range_vals(4) = 1.0_dp
            if (len_trim(input) > 0) read(input, *, iostat=ios) range_vals(3), range_vals(4)
            grid%frac_range(2,1) = range_vals(3); grid%frac_range(2,2) = range_vals(4)

            write(*, '(a)', advance='no') '  fz range (min max, default 0 1): '
            read(*, '(a)', iostat=ios) input
            range_vals(5) = 0.0_dp; range_vals(6) = 1.0_dp
            if (len_trim(input) > 0) read(input, *, iostat=ios) range_vals(5), range_vals(6)
            grid%frac_range(3,1) = range_vals(5); grid%frac_range(3,2) = range_vals(6)
        end if

        ! ── Scan mode ──
        write(*, '(a)') '  Mode: 1=Single-point (SP)  2=Constrained relax (RELAX)'
        write(*, '(a)', advance='no') '  Choice (default SP): '
        read(*, '(a)', iostat=ios) mode_choice
        grid%scan_mode = 'SP'
        if (mode_choice(1:1) == '2') grid%scan_mode = 'RELAX'

        if (grid%scan_mode == 'RELAX') then
            cfg%task_type = TASK_GEOMETRY_OPT
            cfg%cell_opt_mode = 'FIX_CELL'
        else
            cfg%task_type = TASK_ENERGY
        end if

        ! PES 3D constraints: mobile atom fixed in SP mode, fully free in RELAX
        cfg%pes_mobile_idx = mi
        cfg%pes_fix_axes = [1, 1, 1]  ! SP: fully fixed
        if (grid%scan_mode == 'RELAX') then
            cfg%pes_fix_axes = [0, 0, 0]  ! RELAX: fully free in 3D volume
        end if

        ! ── PreCASTEP parameter configuration ──
        write(*, '(a)') ''
        write(*, '(a)') '  ── Configure CASTEP parameters ──'
        write(*, '(a)') '  (Select 0. Generate when ready to proceed)'
        cfg%cif_file_path = trim(file_path)
        call run_main_menu(cfg, ios)
        if (ios == IO_USER_QUIT) return

        ! ── Generate grid points ──
        call generate_pes3d_grid_points(grid, frac_points, n_total, ios)
        if (ios /= 0) then
            write(*, '(a)') '  Error generating grid.'; return
        end if

        ! ── Create output directory ──
        scan_dir = 'pes_scan/' // trim(stem) // '_' // trim(grid%scan_mode)
        sub_dir = 'mkdir -p "' // trim(scan_dir) // '"'
        call execute_command_line(trim(sub_dir), exitstat=ios)

        ! ── Generate files ──
        write(*, '(a,i0,a)') '  Generating ', n_total, ' grid points...'
        do i = 1, n_total
            fx = frac_points(i, 1)
            fy = frac_points(i, 2)
            fz = frac_points(i, 3)

            ! Set mobile atom to grid point fractional coords
            cfg%atom_x(mi) = fx
            cfg%atom_y(mi) = fy
            cfg%atom_z(mi) = fz

            write(sub_dir, '(a, a, i3.3, a, i3.3, a, i3.3)') trim(scan_dir), '/grid_', &
                modulo(i-1, grid%n_points(1)) + 1, '_', &
                modulo((i-1)/grid%n_points(1), grid%n_points(2)) + 1, '_', &
                (i-1) / (grid%n_points(1)*grid%n_points(2)) + 1

            call execute_command_line('mkdir -p "' // trim(sub_dir) // '"', exitstat=ios)

            cell_path  = trim(sub_dir) // '/scan.cell'
            param_path = trim(sub_dir) // '/scan.param'
            call write_cell_file(cell_path, cfg, ios)
            call write_param_file(param_path, cfg, ios)
        end do

        ! ── Write cube file (sym_ops in comment line + placeholder energies) ──
        allocate(dummy_energies(n_total)); dummy_energies = 0.0_dp
        call write_pes3d_cube(trim(scan_dir)//'/pes3d.cube', &
            grid, cfg, cif, dummy_energies, n_total, ios)
        deallocate(dummy_energies)

        ! Save element before deallocation (cif is still valid here)
        mobile_elem = trim(clean_element_symbol(cif%atoms(mi)%element))

        ! ── Cleanup ──
        deallocate(frac_points)
        deallocate(cfg%atom_type, cfg%atom_x, cfg%atom_y, cfg%atom_z)
        call free_cif_data(cif)

        write(*, '(a)') ''
        write(*, '(a)') '  ========================================'
        write(*, '(a)') '    3D PES scan files generated!'
        write(*, '(a)') '  ========================================'
        write(*, '(a,a)')   '  Structure:    ', trim(stem)
        write(*, '(a,a)')   '  Scan mode:    ', trim(grid%scan_mode)
        write(*, '(a,i0,a,a,a)') '  Ref. atom:    #', mi, ' (', trim(mobile_elem), ')'
        write(*, '(a,i0,a,i0,a,i0)') '  Grid:         ', grid%n_points(1), &
            ' x ', grid%n_points(2), ' x ', grid%n_points(3)
        write(*, '(a,i0)')      '  Total points: ', n_total
        write(*, '(a,l1)')      '  Use symmetry: ', grid%use_symmetry
        if (grid%use_symmetry) then
            write(*, '(a,f8.4,a,f8.4,a)') '  fx range:     [', grid%frac_range(1,1), &
                ', ', grid%frac_range(1,2), ']'
            write(*, '(a,f8.4,a,f8.4,a)') '  fy range:     [', grid%frac_range(2,1), &
                ', ', grid%frac_range(2,2), ']'
            write(*, '(a,f8.4,a,f8.4,a)') '  fz range:     [', grid%frac_range(3,1), &
                ', ', grid%frac_range(3,2), ']'
        end if
        write(*, '(a,a)')   '  Directory:    ', trim(scan_dir)
        write(*, '(a)') ''
    end subroutine handle_pes3d_generate


    subroutine handle_pes3d_collect(iostat)
        !! PES 3D sub-menu: Collect .castep results from 3D grid, optionally
        !! symmetry-expand, and launch viewer.
        integer, intent(out) :: iostat

        character(len=MAX_LINE_LEN) :: scan_dir, json_path, input
        integer :: ios
        logical :: exists

        iostat = 0

        write(*, '(a)') '  Enter scan directory (e.g. pes_scan/NaCl_SP):'
        read(*, '(a)', iostat=ios) scan_dir
        if (ios /= 0) return
        scan_dir = adjustl(scan_dir); call strip_quotes(scan_dir)
        if (scan_dir == 'q' .or. scan_dir == 'Q') return

        ! Check for pes3d_metadata.json
        json_path = trim(scan_dir) // '/pes3d_metadata.json'
        inquire(file=trim(json_path), exist=exists)
        if (.not. exists) then
            write(*, '(a)') '  pes3d_metadata.json not found in ' // trim(scan_dir)
            return
        end if

        ! Collect CASTEP energies
        call collect_pes3d_energies(trim(scan_dir), ios)
        if (ios /= 0) then
            write(*, '(a)') '  Error collecting energies.'
            return
        end if

        ! Check if symmetry expansion is needed (from JSON)
        ! If use_symmetry is true, expand
        call symmetry_expand_energies(trim(json_path), ios)

        ! Launch viewer
        write(*, '(a)') ''
        write(*, '(a)', advance='no') '  Launch 3D viewer? (y/n): '
        read(*, '(a)', iostat=ios) input
        if (ios == 0 .and. (input(1:1) == 'y' .or. input(1:1) == 'Y')) then
            call launch_viewer(json_path)
        end if
    end subroutine handle_pes3d_collect


end module poscastep_menu
