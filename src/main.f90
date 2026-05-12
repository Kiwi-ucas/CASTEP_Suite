program CASTEP_Suite
    !! CASTEP Suite: CIF-to-CASTEP converter + post-processing tools
    !! Top-level menu dispatches to PreCASTEP (input generation) or PosCASTEP (post-processing)
    use castep_config, only: dp, castep_config_t, cif_data_t, atom_t, &
        IO_USER_QUIT, MAX_LINE_LEN, default_config, pi
    use cell_writer, only: write_cell_file
    use param_writer, only: write_param_file
    use cli_menu, only: run_main_menu
    use poscastep_menu, only: run_poscastep_menu
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
        case default
            write(*, '(a)') '  Invalid option. Enter 1, 2, or Q.'
        end select
    end do

contains

    subroutine run_precastep_workflow(should_exit)
        !! PreCASTEP workflow: file recognition -> config menu -> generate .cell/.param
        logical, intent(out) :: should_exit
        type(castep_config_t) :: cfg
        type(cif_data_t)     :: cif
        integer               :: istat
        character(len=256)    :: iostat_msg
        character(len=256)    :: formula_chars
        integer               :: i, fidx
        character(len=128)    :: sg_name
        character(len=4)      :: file_ext

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

        cfg%cartesian_coords = (file_ext == 'pdb' .or. file_ext == 'cell')

        write(*, '(a)') '  Computing lattice vectors...'
        cfg%cell_basis = compute_cartesian_lattice(cfg%cell_length(1), cfg%cell_length(2), &
                       cfg%cell_length(3), cfg%cell_angle(1), &
                       cfg%cell_angle(2), cfg%cell_angle(3))

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

    pure function real2str_dp(val) result(s)
        real(dp), intent(in) :: val
        character(20) :: s
        write(s, '(F12.7)') val
    end function real2str_dp

    pure function compute_cartesian_lattice(a, b, c, alpha_deg, beta_deg, gamma_deg) &
        result(lattice)
        real(dp), intent(in) :: a, b, c, alpha_deg, beta_deg, gamma_deg
        real(dp) :: lattice(3, 3)
        real(dp) :: alpha, beta, gamma
        real(dp) :: cos_alpha, cos_beta, cos_gamma, sin_gamma
        real(dp) :: volume_factor

        alpha  = alpha_deg * pi / 180.0_dp
        beta   = beta_deg  * pi / 180.0_dp
        gamma  = gamma_deg * pi / 180.0_dp

        cos_alpha = dcos(alpha)
        cos_beta  = dcos(beta)
        cos_gamma = dcos(gamma)
        sin_gamma = dsin(gamma)

        lattice(1,1) = a
        lattice(2,1) = 0.0_dp
        lattice(3,1) = 0.0_dp

        lattice(1,2) = b * cos_gamma
        lattice(2,2) = b * sin_gamma
        lattice(3,2) = 0.0_dp

        lattice(1,3) = c * cos_beta
        lattice(2,3) = c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma

        volume_factor = 1.0_dp - cos_alpha**2 - cos_beta**2 - cos_gamma**2 &
                        + 2.0_dp * cos_alpha * cos_beta * cos_gamma

        if (volume_factor > 0.0_dp) then
            lattice(3,3) = c * dsqrt(volume_factor) / sin_gamma
        else
            lattice(3,3) = 0.0_dp
        end if
    end function compute_cartesian_lattice

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
