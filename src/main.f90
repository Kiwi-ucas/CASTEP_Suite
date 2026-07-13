program CASTEP_Suite
    !! CASTEP Suite: CIF-to-CASTEP converter + post-processing tools
    !! Top-level menu dispatches to PreCASTEP (input generation) or PosCASTEP (post-processing)
    use castep_config, only: dp, castep_config_t, cif_data_t, atom_t, &
        IO_USER_QUIT, IO_PRECASTEP_LAUNCH, MAX_LINE_LEN, default_config, &
        compute_cartesian_lattice
    use cell_writer, only: write_cell_file
    use param_writer, only: write_param_file
    use symmetry, only: expand_cif_symmetry
    use cli_menu, only: run_main_menu
    use poscastep_menu, only: run_poscastep_menu, precastep_cif_data, has_precastep_data, &
        precastep_source_file, free_cif_data
    use parser, only: parse_cif_inline, parse_pdb_inline, parse_cell_inline, &
        clean_element_symbol
    implicit none

    integer :: istat, ios, mode_choice
    logical :: should_exit
    character(len=MAX_LINE_LEN) :: input

    do
        write(*, '(a)') ''
        write(*, '(a)') '  =================================='
        write(*, '(a)') '            CASTEP Suite'
        write(*, '(a)') '  =================================='
        write(*, '(a)') '  1. PreCASTEP  (generate CASTEP input files)'
        write(*, '(a)') '  2. PosCASTEP  (post-process CASTEP output)'
        write(*, '(a)') '  Q. Quit'
        write(*, '(a)', advance='no') '  Select mode: '

        read(*, '(a)', iostat=ios) input
        if (ios /= 0) exit

        if (len_trim(input) >= 1) then
            if (input(1:1) == 'q' .or. input(1:1) == 'Q') then
                write(*, '(a)') '  Goodbye!'
                exit
            end if
        end if

        read(input, '(I6)', iostat=ios) mode_choice
        if (ios /= 0) then
            write(*, '(a)') '  Invalid input. Enter 1, 2, or Q.'
            cycle
        end if

        select case (mode_choice)
        case (1)
            call run_precastep_workflow(should_exit)
            if (should_exit) exit
        case (2)
            call run_poscastep_menu(istat)
            if (istat == IO_PRECASTEP_LAUNCH) then
                if (has_precastep_data) then
                    call run_precastep_with_cif(precastep_cif_data, should_exit)
                    has_precastep_data = .false.
                    call free_cif_data(precastep_cif_data)
                    if (should_exit) exit
                end if
            end if
        case default
            write(*, '(a)') '  Invalid option. Enter 1, 2, or Q.'
        end select
    end do

contains

    subroutine run_precastep_workflow(should_exit)
        !! PreCASTEP workflow: file recognition -> config menu -> generate .cell/.param
        logical, intent(out) :: should_exit
        type(castep_config_t) :: cfg
        type(cif_data_t)     :: cif, prod_cif, interm_cif
        integer               :: istat
        character(len=256)    :: iostat_msg
        character(len=256)    :: formula_chars
        integer               :: i, fidx
        character(len=128)    :: sg_name
        character(len=4)      :: file_ext, prod_ext, interm_ext

        should_exit = .false.

        call default_config(cfg)

        write(*, '(a)') ''
        write(*, '(a)') '  =================================='
        write(*, '(a)') '              PreCASTEP'
        write(*, '(a)') '  =================================='
        call run_main_menu(cfg, istat)
        if (istat == IO_USER_QUIT) then
            return
        end if
        if (istat /= 0) then
            write(*, '(a)') '  Configuration aborted.'
            return
        end if

        if (len_trim(cfg%cif_file_path) == 0) then
            write(*, '(a)') '  Error: No CIF file specified.'
            return
        end if

        write(*, '(a)') ''
        file_ext = get_file_extension(trim(cfg%cif_file_path))
        cif%n_atoms = 0

        if (file_ext == 'pdb') then
            write(*, '(a)') '  Parsing PDB file: ' // trim(cfg%cif_file_path)
            call parse_pdb_inline(trim(cfg%cif_file_path), cif, istat, iomsg=iostat_msg)
        else if (file_ext == 'cell') then
            write(*, '(a)') '  Parsing .cell file: ' // trim(cfg%cif_file_path)
            call parse_cell_inline(trim(cfg%cif_file_path), cif, istat, iomsg=iostat_msg)
        else
            write(*, '(a)') '  Parsing CIF file: ' // trim(cfg%cif_file_path)
            call parse_cif_inline(trim(cfg%cif_file_path), cif, istat, iomsg=iostat_msg)
        end if
        if (istat /= 0) then
            write(*, '(a)') '  Error parsing file: ' // trim(iostat_msg)
            return
        end if

        if (cif%a <= 0.0_dp .or. cif%b <= 0.0_dp .or. cif%c <= 0.0_dp) then
            write(*, '(a)') '  Error: Missing or invalid cell parameters.'
            return
        end if
        cfg%cell_length(1) = cif%a
        cfg%cell_length(2) = cif%b
        cfg%cell_length(3) = cif%c
        cfg%cell_angle(1)  = cif%alpha
        cfg%cell_angle(2)  = cif%beta
        cfg%cell_angle(3)  = cif%gamma
        write(*, '(a, f10.3, a, f10.3, a, f10.3)') &
            '  Cell: a=', cif%a, ' b=', cif%b, ' c=', cif%c

        if (cif%alpha > 0.0_dp .and. cif%beta > 0.0_dp .and. cif%gamma > 0.0_dp) then
            write(*, '(a, f8.3, a, f8.3, a, f8.3)') &
                '    alpha=', cif%alpha, ' beta=', cif%beta, ' gamma=', cif%gamma
        end if

        sg_name = trim(cif%space_group)
        if (len_trim(sg_name) > 0 .and. sg_name /= 'P1') then
            cfg%has_space_group = .true.
            cfg%space_group_name = sg_name
            write(*, '(a, a)') '  Space group: ', trim(sg_name)
        end if

        ! ── Symmetry expansion: if symops parsed, expand asymmetric unit to full cell ──
        if (cif%n_symops > 1) then
            write(*, '(a, i0)') '  Symmetry operations found: ', cif%n_symops
            call expand_cif_symmetry(cif, istat)
            if (istat == 0) then
                write(*, '(a, i0, a)') '  Expanded to ', cif%n_atoms, ' atom(s) in full unit cell'
            else
                write(*, '(a)') '  Warning: Symmetry expansion failed, using asymmetric unit'
            end if
        end if

        if (cif%n_atoms > 0) then
            formula_chars = ''
            do i = 1, cif%n_atoms
                fidx = len_trim(formula_chars) + 1
                if (fidx + 6 <= len(formula_chars)) then
                    formula_chars(fidx:fidx+5) = ' ' // trim(cif%atoms(i)%element)
                end if
            end do
            cfg%formula_sum = trim(adjustl(formula_chars))
        end if

        if (cif%n_atoms == 0) then
           write(*, '(a)') '  Error: No atom data found in file.'
            return
        end if

        cfg%num_atoms = cif%n_atoms
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
        write(*, '(a, i0)') '  Found ', cif%n_atoms, ' atom(s)'

        cfg%cartesian_coords = .not. cif%positions_fractional

        write(*, '(a)') '  Computing lattice vectors...'
        cfg%cell_basis = compute_cartesian_lattice(cfg%cell_length(1), cfg%cell_length(2), &
                       cfg%cell_length(3), cfg%cell_angle(1), &
                       cfg%cell_angle(2), cfg%cell_angle(3))

        ! ── CINEB: Parse product and intermediate structures ──
        if (trim(cfg%task_type) == 'TRANSITIONSTATESEARCH') then
            ! Parse product structure
            prod_cif%n_atoms = 0
            write(*, '(a)') '  Parsing product file: ' // trim(cfg%prod_file_path)
            prod_ext = get_file_extension(trim(cfg%prod_file_path))
            if (prod_ext == 'pdb') then
                call parse_pdb_inline(trim(cfg%prod_file_path), prod_cif, istat, iomsg=iostat_msg)
            else if (prod_ext == 'cell') then
                call parse_cell_inline(trim(cfg%prod_file_path), prod_cif, istat, iomsg=iostat_msg)
            else
                call parse_cif_inline(trim(cfg%prod_file_path), prod_cif, istat, iomsg=iostat_msg)
            end if
            if (istat /= 0) then
                write(*, '(a)') '  Error parsing product file: ' // trim(iostat_msg)
                return
            end if
            if (prod_cif%n_atoms /= cfg%num_atoms) then
                write(*, '(a, i0, a, i0)') '  Error: Reactant has ', cfg%num_atoms, &
                    ' atoms but product has ', prod_cif%n_atoms
                write(*, '(a)') '  Reactant and product must have the same number of atoms.'
                return
            end if
            cfg%prod_num_atoms = prod_cif%n_atoms
            allocate(cfg%prod_atom_type(cfg%prod_num_atoms))
            allocate(cfg%prod_atom_x(cfg%prod_num_atoms))
            allocate(cfg%prod_atom_y(cfg%prod_num_atoms))
            allocate(cfg%prod_atom_z(cfg%prod_num_atoms))
            do i = 1, prod_cif%n_atoms
                cfg%prod_atom_type(i) = trim(clean_element_symbol(prod_cif%atoms(i)%element))
                cfg%prod_atom_x(i)    = prod_cif%atoms(i)%x
                cfg%prod_atom_y(i)    = prod_cif%atoms(i)%y
                cfg%prod_atom_z(i)    = prod_cif%atoms(i)%z
            end do
            cfg%prod_cartesian_coords = .not. prod_cif%positions_fractional

            ! Parse intermediate structure
            interm_cif%n_atoms = 0
            write(*, '(a)') '  Parsing intermediate file: ' // trim(cfg%interm_file_path)
            interm_ext = get_file_extension(trim(cfg%interm_file_path))
            if (interm_ext == 'pdb') then
                call parse_pdb_inline(trim(cfg%interm_file_path), interm_cif, istat, iomsg=iostat_msg)
            else if (interm_ext == 'cell') then
                call parse_cell_inline(trim(cfg%interm_file_path), interm_cif, istat, iomsg=iostat_msg)
            else
                call parse_cif_inline(trim(cfg%interm_file_path), interm_cif, istat, iomsg=iostat_msg)
            end if
            if (istat /= 0) then
                write(*, '(a)') '  Error parsing intermediate file: ' // trim(iostat_msg)
                return
            end if
            if (interm_cif%n_atoms /= cfg%num_atoms) then
                write(*, '(a, i0, a, i0)') '  Error: Reactant has ', cfg%num_atoms, &
                    ' atoms but intermediate has ', interm_cif%n_atoms
                write(*, '(a)') '  Intermediate must have the same number of atoms as reactant.'
                return
            end if
            cfg%interm_num_atoms = interm_cif%n_atoms
            allocate(cfg%interm_atom_type(cfg%interm_num_atoms))
            allocate(cfg%interm_atom_x(cfg%interm_num_atoms))
            allocate(cfg%interm_atom_y(cfg%interm_num_atoms))
            allocate(cfg%interm_atom_z(cfg%interm_num_atoms))
            do i = 1, interm_cif%n_atoms
                cfg%interm_atom_type(i) = trim(clean_element_symbol(interm_cif%atoms(i)%element))
                cfg%interm_atom_x(i)    = interm_cif%atoms(i)%x
                cfg%interm_atom_y(i)    = interm_cif%atoms(i)%y
                cfg%interm_atom_z(i)    = interm_cif%atoms(i)%z
            end do
            cfg%interm_cartesian_coords = .not. interm_cif%positions_fractional

            write(*, '(a)') '  Product and intermediate structures validated.'
        end if

        write(*, '(a)') ''
        write(*, '(a)') '  Writing CASTEP .cell file: ' // trim(cfg%cell_output_path)
        call write_cell_file(trim(cfg%cell_output_path), cfg, istat, iomsg=iostat_msg)
        if (istat /= 0) then
            write(*, '(a)') '  Error writing .cell file: ' // trim(iostat_msg)
            return
        end if

        write(*, '(a)') '  Writing CASTEP .param file: ' // trim(cfg%param_output_path)
        call write_param_file(trim(cfg%param_output_path), cfg, istat, iomsg=iostat_msg)
        if (istat /= 0) then
            write(*, '(a)') '  Error writing .param file: ' // trim(iostat_msg)
            return
        end if

        write(*, '(a)') ''
        write(*, '(a)') '  =================================='
        write(*, '(a)') '    PreCASTEP conversion complete!'
        write(*, '(a)') '  =================================='
        write(*, '(a)') ''
        write(*, '(a)') '  Summary:'
        write(*, '(a)') '    Cell parameters:  a=' // trim(real2str_dp(cfg%cell_length(1))) // ' ' // &
            'b=' // trim(real2str_dp(cfg%cell_length(2))) // ' ' // &
            'c=' // trim(real2str_dp(cfg%cell_length(3))) // ' Angstrom'
        write(*, '(a)') '    Cell angles:      alpha=' // trim(real2str_dp(cfg%cell_angle(1))) // ' ' // &
            'beta=' // trim(real2str_dp(cfg%cell_angle(2))) // ' ' // &
            'gamma=' // trim(real2str_dp(cfg%cell_angle(3))) // ' deg'
        write(*, '(a, i0)') '    Atoms:            ', cfg%num_atoms
        write(*, '(a, a)') '    XC functional:    ', trim(cfg%xc_functional)
        write(*, '(a, a, a)') '    Cutoff:           ', trim(real2str_dp(cfg%cutoff_energy)), ' eV'
        write(*, '(a, a)') '    Task:             ', trim(cfg%task_type)
        write(*, '(a, a)') '    vdW:              ', trim(cfg%vdw_method)
        write(*, '(a, a)') '    Pseudo:           ', trim(cfg%pseudopotential)
        write(*, '(a)') ''
        write(*, '(a)') '  Output files:'
        write(*, '(a)') '    ' // trim(cfg%cell_output_path)
        write(*, '(a)') '    ' // trim(cfg%param_output_path)
        write(*, '(a)') ''
        should_exit = .true.
    end subroutine run_precastep_workflow

    ! ── Helper: populate castep_config_t from cif_data_t ──
    subroutine populate_cfg_from_cif(cif, cfg)
        type(cif_data_t), intent(in) :: cif
        type(castep_config_t), intent(inout) :: cfg
        integer :: i, fidx
        character(len=256) :: formula_chars

        if (cif%a <= 0.0_dp .or. cif%b <= 0.0_dp .or. cif%c <= 0.0_dp) then
            write(*, '(a)') '  Error: Missing or invalid cell parameters.'
            return
        end if
        cfg%cell_length(1) = cif%a
        cfg%cell_length(2) = cif%b
        cfg%cell_length(3) = cif%c
        cfg%cell_angle(1)  = cif%alpha
        cfg%cell_angle(2)  = cif%beta
        cfg%cell_angle(3)  = cif%gamma
        write(*, '(a, f10.3, a, f10.3, a, f10.3)') &
            '  Cell: a=', cif%a, ' b=', cif%b, ' c=', cif%c
        if (cif%alpha > 0.0_dp .and. cif%beta > 0.0_dp .and. cif%gamma > 0.0_dp) then
            write(*, '(a, f8.3, a, f8.3, a, f8.3)') &
                '    alpha=', cif%alpha, ' beta=', cif%beta, ' gamma=', cif%gamma
        end if

        if (cif%n_atoms == 0) then
            write(*, '(a)') '  Error: No atom data in structure.'
            return
        end if

        cfg%num_atoms = cif%n_atoms
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

        cfg%cartesian_coords = .not. cif%positions_fractional

        ! Build formula summary
        formula_chars = ''
        do i = 1, cif%n_atoms
            fidx = len_trim(formula_chars) + 1
            if (fidx + 6 <= len(formula_chars)) then
                formula_chars(fidx:fidx+5) = ' ' // trim(cif%atoms(i)%element)
            end if
        end do
        cfg%formula_sum = trim(adjustl(formula_chars))

        write(*, '(a, i0)') '  Atoms: ', cif%n_atoms

        write(*, '(a)') '  Computing lattice vectors...'
        cfg%cell_basis = compute_cartesian_lattice(cfg%cell_length(1), cfg%cell_length(2), &
                           cfg%cell_length(3), cfg%cell_angle(1), &
                           cfg%cell_angle(2), cfg%cell_angle(3))
    end subroutine populate_cfg_from_cif

    ! ── PreCASTEP workflow from in-memory cif_data_t (PosCASTEP option 2) ──
    subroutine run_precastep_with_cif(cif, should_exit)
        type(cif_data_t), intent(inout) :: cif
        logical, intent(out) :: should_exit
        type(castep_config_t) :: cfg
        character(len=256) :: iostat_msg
        integer :: istat

        should_exit = .false.
        call default_config(cfg)

        write(*, '(a)') ''
        write(*, '(a)') '  =================================='
        write(*, '(a)') '        PreCASTEP (from Viewer)'
        write(*, '(a)') '  =================================='

        ! Populate cfg directly from cif_data_t (no file parsing needed)
        call populate_cfg_from_cif(cif, cfg)

        ! Use original input file name for output
        cfg%cif_file_path = trim(precastep_source_file)
        cfg%cell_output_path = trim(precastep_source_file) // '.cell'
        cfg%param_output_path = trim(precastep_source_file) // '.param'

        ! Interactive parameter configuration
        call run_main_menu(cfg, istat)
        if (istat == IO_USER_QUIT) then
            write(*, '(a)') '  Configuration cancelled.'
            return
        end if
        if (istat /= 0) then
            write(*, '(a)') '  Configuration aborted.'
            return
        end if

        ! Write output files
        write(*, '(a)') ''
        write(*, '(a)') '  Writing CASTEP .cell file: ' // trim(cfg%cell_output_path)
        call write_cell_file(trim(cfg%cell_output_path), cfg, istat, iomsg=iostat_msg)
        if (istat /= 0) then
            write(*, '(a)') '  Error writing .cell file: ' // trim(iostat_msg)
            return
        end if

        write(*, '(a)') '  Writing CASTEP .param file: ' // trim(cfg%param_output_path)
        call write_param_file(trim(cfg%param_output_path), cfg, istat, iomsg=iostat_msg)
        if (istat /= 0) then
            write(*, '(a)') '  Error writing .param file: ' // trim(iostat_msg)
            return
        end if

        write(*, '(a)') ''
        write(*, '(a)') '  =================================='
        write(*, '(a)') '    PreCASTEP conversion complete!'
        write(*, '(a)') '  =================================='
        write(*, '(a)') ''
        write(*, '(a)') '  Summary:'
        write(*, '(a)') '    Cell parameters:  a=' // trim(real2str_dp(cfg%cell_length(1))) // ' ' // &
            'b=' // trim(real2str_dp(cfg%cell_length(2))) // ' ' // &
            'c=' // trim(real2str_dp(cfg%cell_length(3))) // ' Angstrom'
        write(*, '(a)') '    Cell angles:      alpha=' // trim(real2str_dp(cfg%cell_angle(1))) // ' ' // &
            'beta=' // trim(real2str_dp(cfg%cell_angle(2))) // ' ' // &
            'gamma=' // trim(real2str_dp(cfg%cell_angle(3))) // ' deg'
        write(*, '(a, i0)') '    Atoms:            ', cfg%num_atoms
        write(*, '(a, a)') '    XC functional:    ', trim(cfg%xc_functional)
        write(*, '(a, a, a)') '    Cutoff:           ', trim(real2str_dp(cfg%cutoff_energy)), ' eV'
        write(*, '(a, a)') '    Task:             ', trim(cfg%task_type)
        write(*, '(a)') ''
        write(*, '(a)') '  Output files:'
        write(*, '(a)') '    ' // trim(cfg%cell_output_path)
        write(*, '(a)') '    ' // trim(cfg%param_output_path)
        write(*, '(a)') ''
        should_exit = .true.
    end subroutine run_precastep_with_cif

    pure function real2str_dp(val) result(s)
        real(dp), intent(in) :: val
        character(20) :: s
        write(s, '(F12.7)') val
        s = adjustl(s)
    end function real2str_dp

    function get_file_extension(filename) result(ext)
        character(len=*), intent(in) :: filename
        character(len=4) :: ext
        integer :: dot_pos, n

        n = len_trim(filename)
        dot_pos = 0
        do while (dot_pos < n)
            n = n - 1
            if (filename(n:n) == '.') then
                dot_pos = n
                exit
            end if
        end do

        if (dot_pos == 0 .or. dot_pos >= len_trim(filename)) then
            ext = ''
            return
        end if

        ext = adjustl(filename(dot_pos+1:))
        call to_lower_inline(ext)
    end function get_file_extension

    pure subroutine to_lower_inline(s)
        character(len=*), intent(inout) :: s
        integer :: i, ic
        do i = 1, len(s)
            ic = iachar(s(i:i))
            if (ic >= iachar('A') .and. ic <= iachar('Z')) then
                s(i:i) = char(ic + 32)
            end if
        end do
    end subroutine to_lower_inline

end program CASTEP_Suite
