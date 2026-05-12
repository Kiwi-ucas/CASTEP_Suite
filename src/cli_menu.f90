module cli_menu
    !! Interactive CLI menus for PreCASTEP configuration
    !! Main menu loop with cached state in castep_config_t.
    !! Simplified: no Basis Precision, Symmetry = NONE/AUTO only.
    use castep_config, only: &
         dp, TASK_ENERGY, TASK_GEOMETRY_OPT, &
         TASK_MOLECULAR_DYN, TASK_PHONON, TASK_ELECTRONIC_SPECTRO, &
         TASK_TRANSITION_STATE, TASK_EFIELD, TASK_THERMODYNAMICS, &
         TASK_MAGRES, TASK_SPECTRAL, TASK_EPCOUPLING, TASK_PHONON_EFIELD, TASK_GENETIC_ALGO, &
         TASK_SOCKET_DRIVER, TASK_ELASTIC, TASK_AUTOSOLVATION, &
         FUNC_PBE, FUNC_HSE06, FUNC_PBEsol, &
         FUNC_PBE0, FUNC_r2scan, &
         TOL_SUPERFINE, TOL_FINE, TOL_NORMAL, TOL_COARSE, &
         GEO_COARSE, GEO_MEDIUM, GEO_FINE, GEO_EXTREME, &
         CELL_ALL, CELL_INTE, &
         PHONON_FINE_NONE, PHONON_FINE_SUPERCELL, PHONON_FINE_INTERPOLATE, &
         PHONON_METHOD_DFPT, PHONON_METHOD_FD, &
         PHONON_DFPT_DM, PHONON_DFPT_ALLBANDS, &
         PHONON_CUTOFF_CUMULANT, PHONON_CUTOFF_SPHERICAL, &
         PHONON_QPOINT_MP_GRID, PHONON_QPOINT_PATH, &
         PHONON_SUM_NONE, PHONON_SUM_RECIPROCAL, PHONON_SUM_REALSPACE, &
         PHONON_SUM_REAL_RECIP, PHONON_SUM_MOLECULAR, &
         SYM_AUTO, SYM_NONE, &
         VDW_NONE, VDW_D3, VDW_D3_BJ, VDW_D4, &
         PSEUDO_NCP19, PSEUDO_C19MK2, PSEUDO_SOC19, &
         KPOINT_GAMMA, KPOINT_MONKHORST_PACK, &
         OPT_BFGS, OPT_LBFGS, OPT_CG, &
         MAX_LINE_LEN, IO_INVALID_INPUT, IO_SUCCESS, IO_USER_QUIT, &
         castep_config_t, default_config, compare_tags, int2str, strip_quotes
    implicit none
    private

    public :: run_main_menu
    public :: ask_input_file
    public :: ask_task_type
    public :: ask_xc_functional
    public :: ask_cutoff_energy
    public :: ask_vdw_method
    public :: ask_pseudopotential
    public :: ask_kpoint_scheme
    public :: ask_scf_tolerance
    public :: ask_optimizer
    public :: ask_cell_opt_mode
    public :: ask_symmetry_source
    public :: ask_geom_tolerance
    public :: ask_advanced_options
    public :: ask_spectral_task

    ! Menu option labels
    integer, parameter :: &
         MENU_GENERATE    = 0, &
         MENU_TASK        = 1, &
         MENU_XC          = 2, &
         MENU_CUTOFF      = 3, &
         MENU_VDW         = 4, &
         MENU_PSEUDO      = 5, &
         MENU_KPOINT      = 6, &
         MENU_SCF         = 7, &
         MENU_SYMMETRY    = 8, &
         MENU_OPTIMIZER   = 9, &
         MENU_QUIT        = -1, &
         MENU_ADVANCED    = -2, &
         MENU_NONLINEAR   = -3

contains

    subroutine run_main_menu(cfg, iostat)
        !! Main menu loop with cached state
        type(castep_config_t), intent(inout) :: cfg
        integer, intent(out) :: iostat
        integer :: choice, ios
        character(len=MAX_LINE_LEN) :: input
        logical :: cif_ready, needs_geo_params, needs_phonon_params, show_spectral_menu
        character(len=32) :: spectral_display
        character(len=512) :: base_name

        iostat = 0
        cif_ready = .false.

        ! Ask CIF path first
        call ask_input_file('Input structure file path (.cif/.pdb/.cell): ', cfg%cif_file_path, iostat)
        if (iostat /= 0) return
        cif_ready = .true.

        do
            ! Determine if geo params are needed for display
            show_spectral_menu = (trim(cfg%task_type) == TASK_ELECTRONIC_SPECTRO)
            spectral_display = cfg%spectral_task_type
            needs_geo_params   = (trim(cfg%task_type) == TASK_GEOMETRY_OPT)
            needs_phonon_params = (trim(cfg%task_type) == TASK_PHONON &
                              .or. trim(cfg%task_type) == TASK_PHONON_EFIELD &
                              .or. trim(cfg%task_type) == TASK_THERMODYNAMICS &
                              .or. trim(cfg%task_type) == TASK_EFIELD)

            write(*, '(a)') ''
            write(*, '(a)') '  ================================'
            write(*, '(a)') '             PreCASTEP'
            write(*, '(a)') '  ================================'
            if (cif_ready) then
                write(*, '(a, a)') '  CIF: ', trim(cfg%cif_file_path)
            end if
            if (needs_phonon_params .and. &
                (trim(cfg%task_type) == TASK_EFIELD &
                .or. trim(cfg%task_type) == TASK_PHONON_EFIELD)) then
                write(*, '(a, a, a)') ' -3. Nonlinear optics        (', &
                    trim(cfg%efield_calculate_nonlinear), ')'
            end if
            write(*, '(a)') ' -2. Advanced option'
            write(*, '(a)') ' -1. Spin_polarized : ' //   &
                   trim(sp_label(cfg%spin_polarized))
            write(*, '(a)') '  0. Generate InputFile'
            write(*, '(a)') '  1. Task                    (' // trim(task_label(cfg%task_type))  // ')'
            write(*, '(a)') '  2. XC functional           (' // trim(cfg%xc_functional)           // ')'
            write(*, '(a)') '  3. Cutoff energy (eV)      (' // trim(cutoff_label(cfg%cutoff_energy)) // ')'
            write(*, '(a)') '  4. vdW correction          (' // trim(cfg%vdw_method)              // ')'
            write(*, '(a)') '  5. Pseudopotential         (' // trim(cfg%pseudopotential)           // ')'
            write(*, '(a)') '  6. K-point                 (' // trim(kpoint_label(cfg%kpoint_scheme)) // ')'
            write(*, '(a)') '  7. SCF tolerance           (' // trim(scf_label(cfg%scf_tolerance))  // ')'
            write(*, '(a)') '  8. Symmetry                (' // trim(sym_label(cfg%sym_source))     // ')'
            if (needs_geo_params) then
                write(*, '(a)') '  9. Optimizer              (' // trim(cfg%optimizer)             // ')'
                write(*, '(a)') ' 10. Cell opt mode          (' // trim(cfg%cell_opt_mode)         // ')'
                write(*, '(a)') ' 11. Geo tolerance          (' // trim(geom_tol_label(cfg%geom_tolerance)) // ')'
            end if
            if (needs_phonon_params) then
                write(*, '(a, a, a)')  '  9. Phonon q-point scheme   (', &
                    trim(qpoint_label(cfg%phonon_qpoint_scheme, cfg%phonon_kpoint_mp_grid)), ')'
                write(*, '(a, a, a)')  ' 10. Phonon method           (', trim(cfg%phonon_method), ')'
                write(*, '(a, a, a)')  ' 11. Phonon fine method      (', trim(cfg%phonon_fine_method), ')'
                write(*, '(a, es9.1, a)') ' 12. Phonon energy tol       (', cfg%phonon_energy_tol, ')'
                if (trim(cfg%phonon_method) == PHONON_METHOD_FD) then
                    write(*, '(a)') ' 13. Phonon supercell matrix'
                end if
                if (trim(cfg%phonon_fine_method) /= PHONON_FINE_NONE) then
                    write(*, '(a, a, a)') ' 14. Phonon fine q-point     (', &
                        trim(qpoint_label(cfg%phonon_fine_qpoint_scheme, cfg%phonon_fine_kpoint_mp_grid)), ')'
                end if
            end if
            if (show_spectral_menu) then
                write(*, '(a)') ' 20. Spectral Task          (' // trim(spectral_display) // ')'
            end if
            write(*, '(a)') '  Q. Back'
            write(*, '(a)') ' Select option : '

            read(*, '(a)', iostat=ios) input
            if (ios /= 0) then
                iostat = IO_INVALID_INPUT
                write(*, '(a)') '  Invalid input.'
                cycle
            end if

            ! Check for quit
            if (len_trim(input) >= 1) then
                if (input(1:1) == 'q' .or. input(1:1) == 'Q') then
                    write(*, '(a)') '  Aborted.'
                    iostat = IO_USER_QUIT
                    return
                end if
            end if

            read(input, '(I6)', iostat=ios) choice
            if (ios /= 0) then
                iostat = IO_INVALID_INPUT
                write(*, '(a)') '  Invalid input. Enter a number.'
                cycle
            end if

            select case (choice)
            case (MENU_GENERATE)
                base_name = auto_output_name(cfg%cif_file_path, cfg%task_type)
                cfg%cell_output_path  = trim(base_name) // '.cell'
                cfg%param_output_path = trim(base_name) // '.param'
                write(*, '(a, a)') '  Output: ', trim(base_name)
                return
            case (MENU_TASK)
                call ask_task_type('Select calculation type: ', cfg%task_type, iostat)
                if (iostat /= 0) return
                if (trim(cfg%task_type) == TASK_THERMODYNAMICS) then
                    cfg%phonon_fine_method = PHONON_FINE_SUPERCELL
                    cfg%phonon_method = PHONON_METHOD_FD
                else if (trim(cfg%task_type) == TASK_PHONON &
                    .and. trim(cfg%phonon_fine_method) == PHONON_FINE_SUPERCELL) then
                    cfg%phonon_fine_method = PHONON_FINE_INTERPOLATE
                end if
                if (trim(cfg%task_type) == TASK_EFIELD) then
                    cfg%pseudopotential = PSEUDO_NCP19
                else if ((trim(cfg%task_type) == TASK_PHONON &
                    .or. trim(cfg%task_type) == TASK_PHONON_EFIELD) &
                    .and. trim(cfg%phonon_method) == PHONON_METHOD_DFPT) then
                    cfg%pseudopotential = PSEUDO_NCP19
                end if
                if (trim(cfg%task_type) == TASK_PHONON_EFIELD) then
                    cfg%phonon_calc_lo_to_splitting = .true.
                end if
            case (MENU_XC)
                call ask_xc_functional('Select XC functional: ', cfg%xc_functional, iostat)
                if (iostat /= 0) return
            case (MENU_CUTOFF)
                call ask_cutoff_energy('Plane-wave cutoff energy (eV): ', cfg%cutoff_energy, iostat)
                if (iostat /= 0) return
            case (MENU_VDW)
                call ask_vdw_method('Select vdW correction: ', cfg%vdw_method, iostat)
                if (iostat /= 0) return
            case (MENU_PSEUDO)
                call ask_pseudopotential('Select pseudopotential: ', cfg%pseudopotential, iostat)
                if (iostat /= 0) return
                if (trim(cfg%pseudopotential) == PSEUDO_SOC19) then
                    cfg%spin_polarized = .true.
                end if
            case (MENU_KPOINT)
                call ask_kpoint_scheme('Select K-point scheme: ', cfg%kpoint_scheme, cfg%kpoint_grid, iostat)
                if (iostat /= 0) return
            case (MENU_SCF)
                call ask_scf_tolerance('SCF tolerance: ', cfg%scf_tolerance, iostat)
                if (iostat /= 0) return
            case (MENU_SYMMETRY)
                call ask_symmetry_source('Select symmetry handling: ', cfg%sym_source, iostat)
                if (iostat /= 0) return
            case (MENU_OPTIMIZER)
                if (needs_phonon_params) then
                    call ask_phonon_qpoint_scheme(cfg, iostat)
                else
                    call ask_optimizer('Select optimizer: ', cfg%optimizer, iostat)
                end if
                if (iostat /= 0) return
            case (10)
                if (needs_phonon_params) then
                    call ask_phonon_method('Select phonon method (1=DFPT, 2=FiniteDisp): ', &
                        cfg%phonon_method, iostat)
                    if (trim(cfg%phonon_method) == PHONON_METHOD_DFPT) then
                        cfg%pseudopotential = PSEUDO_NCP19
                        if (trim(cfg%phonon_fine_method) == PHONON_FINE_SUPERCELL) &
                            cfg%phonon_fine_method = PHONON_FINE_NONE
                    end if
                else
                    call ask_cell_opt_mode('Select cell optimization mode: ', cfg%cell_opt_mode, iostat)
                end if
                if (iostat /= 0) return
            case (11)
                if (needs_phonon_params) then
                    call ask_phonon_fine_method('Select phonon fine method: ', &
                        cfg%phonon_fine_method, cfg%phonon_method, iostat)
                    if (trim(cfg%phonon_fine_method) == PHONON_FINE_SUPERCELL) then
                        cfg%phonon_method = PHONON_METHOD_FD
                    end if
                else
                    call ask_geom_tolerance('Select geometry optimization tolerance: ', cfg%geom_tolerance, iostat)
                end if
                if (iostat /= 0) return
            case (12)
                if (needs_phonon_params) then
                    call ask_phonon_energy_tol('Phonon energy tolerance (eV/A^2): ', &
                        cfg%phonon_energy_tol, iostat)
                else
                    call ask_geom_tolerance('Select geometry optimization tolerance: ', cfg%geom_tolerance, iostat)
                end if
                if (iostat /= 0) return
            case (13)
                if (needs_phonon_params) then
                    call ask_phonon_supercell_matrix('Phonon supercell matrix (3x3): ', &
                        cfg%phonon_supercell_matrix, iostat)
                end if
                if (iostat /= 0) return
            case (14)
                if (needs_phonon_params) then
                    call ask_phonon_fine_qpoint_scheme(cfg, iostat)
                end if
                if (iostat /= 0) return
            case (20)
                call ask_spectral_task('Select spectral task type: ', cfg%spectral_task_type, iostat)
                if (iostat /= 0) return
            case (MENU_ADVANCED)
                if (needs_phonon_params) then
                    call ask_advanced_options('Advanced options: ', cfg%smearing, &
                        cfg%max_scf_cycles, cfg%elec_convergence_win, &
                        cfg%calculate_elf, cfg%calculate_edd, iostat, &
                        cfg%phonon_calculate_dos, cfg%phonon_dos_spacing, cfg%phonon_sum_rule_method, &
                        cfg%phonon_finite_disp, cfg%phonon_max_cycles, cfg%phonon_dfpt_method, &
                        cfg%phonon_write_force_constants, cfg%phonon_write_dynamical, &
                        cfg%phonon_calc_lo_to_splitting, cfg%phonon_force_constant_cutoff, &
                        cfg%phonon_fine_cutoff_method, cfg%phonon_dos_limit, &
                        cfg%phonon_max_cg_steps, cfg%phonon_use_kpoint_symmetry, &
                        cfg%calculate_born_charges, cfg%calculate_raman, &
                        cfg%raman_method, cfg%efield_dfpt_method, &
                        cfg%efield_max_cycles, cfg%efield_energy_tol, cfg%efield_convergence_win, &
                        cfg%efield_freq_spacing, cfg%efield_oscillator_q, cfg%efield_calc_ion_permittivity, &
                        cfg%efield_ignore_molec_modes, &
                        (trim(cfg%phonon_method) == PHONON_METHOD_DFPT), &
                        (trim(cfg%task_type) == TASK_EFIELD .or. trim(cfg%task_type) == TASK_PHONON_EFIELD))
                else
                    call ask_advanced_options('Advanced options: ', cfg%smearing, &
                        cfg%max_scf_cycles, cfg%elec_convergence_win, &
                        cfg%calculate_elf, cfg%calculate_edd, iostat)
                end if
                if (iostat /= 0) return
            case (MENU_NONLINEAR)
                if (trim(cfg%efield_calculate_nonlinear) == 'NONE') then
                    cfg%efield_calculate_nonlinear = 'CHI2'
                else
                    cfg%efield_calculate_nonlinear = 'NONE'
                end if
            case (MENU_QUIT)
                cfg%spin_polarized = .not. cfg%spin_polarized
            case default
                write(*, '(a)') '  Invalid option. Enter a number.'
            end select
        end do
    end subroutine run_main_menu

    subroutine ask_input_file(prompt_text, result_path, iostat)
        character(len=*), intent(in)  :: prompt_text
        character(len=*), intent(out) :: result_path
        integer, intent(out)           :: iostat
        integer :: ios
        logical :: exists
        character(len=MAX_LINE_LEN) :: input

        iostat = 0
        do
            write(*, '(a)') trim(prompt_text)
            read(*, '(a)', iostat=ios) input
            if (ios /= 0) then
                iostat = IO_INVALID_INPUT
                write(*, '(a)') '  Error reading input.'
                cycle
            end if

            result_path = adjustl(trim(input))
            call strip_quotes(result_path)
            if (len_trim(result_path) == 0) then
                write(*, '(a)') '  Path cannot be empty. Please try again.'
                cycle
            end if

            inquire(file=trim(result_path), exist=exists)
            if (.not. exists) then
                write(*, '(a)') '  File not found. Please try again.'
                cycle
            end if

            exit
        end do
    end subroutine ask_input_file


    subroutine ask_task_type(prompt_text, result_task, iostat)
        character(len=*), intent(in)  :: prompt_text
        character(len=*), intent(out) :: result_task
        integer, intent(out)           :: iostat
        integer :: choice, ios
        character(len=MAX_LINE_LEN) :: input

        iostat = 0
        result_task = TASK_ENERGY

        do
            write(*, '(a)') ''
            write(*, '(a)') '  Select calculation type:'
            write(*, '(a)') '    1. Energy'
            write(*, '(a)') '    2. GeometryOptimisation'
            write(*, '(a)') '    3. ElectronicSpectroscopy'
            write(*, '(a)') '    4. Phonon'
            write(*, '(a)') '    5. Phonon+Efield'
            write(*, '(a)') '    6. Efield'
            write(*, '(a)') '    7. Thermodynamics'
            write(*, '(a)') '    8. MolecularDynamics     (暂未开发)'
            write(*, '(a)') '    9. TransitionState       (暂未开发)'
            write(*, '(a)') '   10. MagneticResponse      (暂未开发)'
            write(*, '(a)') '   11. Spectral              (暂未开发)'
            write(*, '(a)') '   12. Elastic               (暂未开发)'
            write(*, '(a)') '   13. GeneticAlgorithm      (暂未开发)'
            write(*, '(a)') '   14. SocketDriver          (暂未开发)'
            write(*, '(a)') '   15. Autosolvation         (暂未开发)'
            write(*, '(a)') '   16. EpCoupling            (暂未开发)'
            write(*, '(a)') '    Enter choice : ' // trim(prompt_text)

            read(*, '(a)', iostat=ios) input
            if (ios /= 0 .or. len_trim(input) == 0) then
                ! Empty or invalid input → accept default
                exit
            end if

            read(input, '(I6)', iostat=ios) choice
            if (ios /= 0) then
                write(*, '(a)') '  Invalid input. Please enter a number 1-17.'
                cycle
            end if

            select case (choice)
            case (1)  ; result_task = TASK_ENERGY
            case (2)  ; result_task = TASK_GEOMETRY_OPT
            case (3)  ; result_task = TASK_ELECTRONIC_SPECTRO
            case (4)  ; result_task = TASK_PHONON
            case (5)  ; result_task = TASK_PHONON_EFIELD
            case (6)  ; result_task = TASK_EFIELD
            case (7)  ; result_task = TASK_THERMODYNAMICS
            case (8)  ; result_task = TASK_MOLECULAR_DYN
            case (9)  ; result_task = TASK_TRANSITION_STATE
            case (10) ; result_task = TASK_MAGRES
            case (11) ; result_task = TASK_SPECTRAL
            case (12) ; result_task = TASK_ELASTIC
            case (13) ; result_task = TASK_GENETIC_ALGO
            case (14) ; result_task = TASK_SOCKET_DRIVER
            case (15) ; result_task = TASK_AUTOSOLVATION
            case (16) ; result_task = TASK_EPCOUPLING
            case default
                write(*, '(a)') '  Invalid choice. Please enter 1-16.'
                cycle
            end select

            exit
        end do
    end subroutine ask_task_type

    subroutine ask_xc_functional(prompt_text, result_func, iostat)
        character(len=*), intent(in)  :: prompt_text
        character(len=*), intent(out) :: result_func
        integer, intent(out)           :: iostat
        integer :: choice, ios
        character(len=MAX_LINE_LEN) :: input

        iostat = 0
        result_func = FUNC_PBE

        do
            write(*, '(a)') ''
            write(*, '(a)') '  Select XC functional:'
            write(*, '(a)') '    1. PBE'
            write(*, '(a)') '    2. PBEsol'
            write(*, '(a)') '    3. HSE06'
            write(*, '(a)') '    4. PBE0'
            write(*, '(a)') '    5. r2scan'
            write(*, '(a)') '    Enter choice : ' // trim(prompt_text)

            read(*, '(a)', iostat=ios) input
            if (ios /= 0 .or. len_trim(input) == 0) then
                ! Empty or invalid input → accept default
                exit
            end if

            read(input, '(I6)', iostat=ios) choice
            if (ios /= 0) then
                write(*, '(a)') '  Invalid input. Please enter a number 1-5.'
                cycle
            end if

            select case (choice)
            case (1)  ; result_func = FUNC_PBE
            case (2)  ; result_func = FUNC_PBEsol
            case (3)  ; result_func = FUNC_HSE06
            case (4)  ; result_func = FUNC_PBE0
            case (5)  ; result_func = FUNC_r2scan
            case default
                write(*, '(a)') '  Invalid choice. Please enter 1-5.'
                cycle
            end select

            exit
        end do
    end subroutine ask_xc_functional

    subroutine ask_cutoff_energy(prompt_text, result_cutoff, iostat)
        character(len=*), intent(in)  :: prompt_text
        real(dp), intent(out)         :: result_cutoff
        integer, intent(out)           :: iostat
        integer :: ios, ios2, val_int
        character(len=MAX_LINE_LEN) :: input

        iostat = 0
        result_cutoff = 400.0_dp

        do
            write(*, '(a, f8.1)') '  Plane-wave cutoff energy (eV) [default: 400.0]: ' // trim(prompt_text)
            read(*, '(a)', iostat=ios) input
            if (ios /= 0) exit  ! empty input -> default

            if (len_trim(adjustl(input)) == 0) exit

            ! Handle integer input
            if (index(input, '.') == 0 .and. index(input, 'd') == 0 .and. index(input, 'e') == 0) then
                read(input, '(I12)', iostat=ios2) val_int
                if (ios2 == 0) then
                    result_cutoff = dble(val_int)
                end if
            else
                read(input, '(F20.10)', iostat=ios2) result_cutoff
            end if
            if (ios2 /= 0 .or. result_cutoff <= 0.0_dp) then
                write(*, '(a)') '  Invalid value. Please enter a positive number.'
                cycle
            end if

            exit
        end do
    end subroutine ask_cutoff_energy

    subroutine ask_vdw_method(prompt_text, result_vdw, iostat)
        character(len=*), intent(in)  :: prompt_text
        character(len=*), intent(out) :: result_vdw
        integer, intent(out)           :: iostat
        integer :: choice, ios
        character(len=MAX_LINE_LEN) :: input

        iostat = 0
        result_vdw = VDW_NONE

        do
            write(*, '(a)') ''
            write(*, '(a)') '  Select vdW correction:'
            write(*, '(a)') '    1. NONE'
            write(*, '(a)') '    2. D3'
            write(*, '(a)') '    3. D3-BJ'
            write(*, '(a)') '    4. D4'
            write(*, '(a)') '    Enter choice : ' // trim(prompt_text)

            read(*, '(a)', iostat=ios) input
            if (ios /= 0) then
                iostat = IO_INVALID_INPUT
                write(*, '(a)') '  Invalid input. Please enter a number 1-4.'
                cycle
            end if

            read(input, '(I6)', iostat=ios) choice
            if (ios /= 0) then
                write(*, '(a)') '  Invalid input. Please enter a number 1-4.'
                cycle
            end if

            select case (choice)
            case (1) ; result_vdw = VDW_NONE
            case (2) ; result_vdw = VDW_D3
            case (3) ; result_vdw = VDW_D3_BJ
            case (4) ; result_vdw = VDW_D4
            case default
                write(*, '(a)') '  Invalid choice. Please enter 1-4.'
                cycle
            end select

            exit
        end do
    end subroutine ask_vdw_method

    subroutine ask_pseudopotential(prompt_text, result_pseudo, iostat)
        character(len=*), intent(in)  :: prompt_text
        character(len=*), intent(out) :: result_pseudo
        integer, intent(out)           :: iostat
        integer :: choice, ios
        character(len=MAX_LINE_LEN) :: input

        iostat = 0
        result_pseudo = PSEUDO_C19MK2

        do
            write(*, '(a)') ''
            write(*, '(a)') '  Select pseudopotential:'
            write(*, '(a)') '    1. NCP19     (Norm-conserving)'
            write(*, '(a)') '    2. C19MK2    (Ultrasoft)'
            write(*, '(a)') '    3. SOC19     (Spin-orbit coupling)'
            write(*, '(a)') '    Enter choice : ' // trim(prompt_text)

            read(*, '(a)', iostat=ios) input
            if (ios /= 0) then
                iostat = IO_INVALID_INPUT
                write(*, '(a)') '  Invalid input. Please enter a number 1-3.'
                cycle
            end if

            read(input, '(I6)', iostat=ios) choice
            if (ios /= 0) then
                write(*, '(a)') '  Invalid input. Please enter a number 1-3.'
                cycle
            end if

            select case (choice)
            case (1) ; result_pseudo = PSEUDO_NCP19
            case (2) ; result_pseudo = PSEUDO_C19MK2
            case (3) ; result_pseudo = PSEUDO_SOC19
            case default
                write(*, '(a)') '  Invalid choice. Please enter 1-3.'
                cycle
            end select

            exit
        end do
    end subroutine ask_pseudopotential

    subroutine ask_kpoint_scheme(prompt_text, result_scheme, result_grid, iostat)
        character(len=*), intent(in)  :: prompt_text
        character(len=*), intent(out) :: result_scheme
        integer, intent(out)            :: result_grid(3)
        integer, intent(out)           :: iostat
        integer :: choice, ios
        character(len=MAX_LINE_LEN) :: input

        iostat = 0
        result_scheme = KPOINT_GAMMA
        result_grid = 0

        do
            write(*, '(a)') ''
            write(*, '(a)') '  Select K-point scheme:'
            write(*, '(a)') '    1. GAMMA-only'
            write(*, '(a)') '    2. Monkhorst-Pack grid'
            write(*, '(a)') '    Enter choice : ' // trim(prompt_text)

            read(*, '(a)', iostat=ios) input
            if (ios /= 0) then
                iostat = IO_INVALID_INPUT
                write(*, '(a)') '  Invalid input. Please enter a number 1-2.'
                cycle
            end if

            read(input, '(I6)', iostat=ios) choice
            if (ios /= 0) then
                write(*, '(a)') '  Invalid input. Please enter a number 1-2.'
                cycle
            end if

            select case (choice)
            case (1)
                result_scheme = KPOINT_GAMMA
                result_grid = 0
                exit
            case (2)
                result_scheme = KPOINT_MONKHORST_PACK
                write(*, '(a)') '  Enter MP grid [4,4,4]:'
                ! Read as string and split on commas
                read(*, '(a)', iostat=ios) input
                if (ios /= 0) then
                    result_grid = 4
                    write(*, '(a)') '  Using default: 4, 4, 4'
                    exit
                end if
                ! Parse x,y,z from comma-separated string
                call parse_grid_comma(input, result_grid, ios)
                if (ios /= 0) then
                    result_grid = 4
                    write(*, '(a)') '  Using default: 4, 4, 4'
                end if
                exit
            case default
                write(*, '(a)') '  Invalid choice. Please enter 1-2.'
                cycle
            end select
        end do
    end subroutine ask_kpoint_scheme

    subroutine ask_scf_tolerance(prompt_text, result_tol, iostat)
        character(len=*), intent(in)  :: prompt_text
        character(len=*), intent(out) :: result_tol
        integer, intent(out)           :: iostat
        integer :: ios, ios2
        character(len=MAX_LINE_LEN) :: input
        character(len=MAX_LINE_LEN) :: tmp
        real(dp) :: dummy_val

        iostat = 0
        result_tol = '1e-5'

        do
            write(*, '(a, a)') '  SCF tolerance : ' // trim(prompt_text)
            read(*, '(a)', iostat=ios) input
            if (ios /= 0) exit  ! empty -> default

            if (len_trim(adjustl(input)) == 0) exit

            tmp = adjustl(trim(input))

            ! Handle scientific notation or plain number
            if (index(tmp, 'E') > 0 .or. index(tmp, 'e') > 0) then
                ! Directly save user input for scientific notation
                read(tmp, '(E22.10)', iostat=ios2) dummy_val
                if (ios2 == 0 .and. dummy_val > 0.0_dp) then
                    result_tol = tmp
                else
                    write(*, '(a)') '  Invalid value. Please enter a positive number (e.g., 1e-5).'
                    cycle
                end if
            else
                ! Try plain number: validate then accept directly
                read(tmp, '(F22.15)', iostat=ios2) dummy_val
                if (ios2 == 0 .and. dummy_val > 0.0_dp) then
                    result_tol = tmp
                else
                    write(*, '(a)') '  Invalid value. Please enter a positive number (e.g., 1e-5).'
                    cycle
                end if
            end if

            exit
        end do
    end subroutine ask_scf_tolerance

    subroutine ask_optimizer(prompt_text, result_opt, iostat)
        character(len=*), intent(in)  :: prompt_text
        character(len=*), intent(out) :: result_opt
        integer, intent(out)           :: iostat
        integer :: choice, ios
        character(len=MAX_LINE_LEN) :: input

        iostat = 0
        result_opt = OPT_BFGS

        do
            write(*, '(a)') ''
            write(*, '(a)') '  Select optimizer:'
            write(*, '(a)') '    1. BFGS'
            write(*, '(a)') '    2. LBFGS'
            write(*, '(a)') '    3. CG'
            write(*, '(a)') '    Enter choice : ' // trim(prompt_text)

            read(*, '(a)', iostat=ios) input
            if (ios /= 0) then
                iostat = IO_INVALID_INPUT
                write(*, '(a)') '  Invalid input. Please enter a number 1-3.'
                cycle
            end if

            read(input, '(I6)', iostat=ios) choice
            if (ios /= 0) then
                write(*, '(a)') '  Invalid input. Please enter a number 1-3.'
                cycle
            end if

            select case (choice)
            case (1) ; result_opt = OPT_BFGS
            case (2) ; result_opt = OPT_LBFGS
            case (3) ; result_opt = OPT_CG
            case default
                write(*, '(a)') '  Invalid choice. Please enter 1-3.'
                cycle
            end select

            exit
        end do
    end subroutine ask_optimizer

    subroutine ask_cell_opt_mode(prompt_text, result_mode, iostat)
        character(len=*), intent(in)  :: prompt_text
        character(len=*), intent(out) :: result_mode
        integer, intent(out)           :: iostat
        integer :: choice, ios
        character(len=MAX_LINE_LEN) :: input

        iostat = 0
        result_mode = CELL_INTE

        do
            write(*, '(a)') ''
            write(*, '(a)') '  Select cell optimization mode:'
            write(*, '(a)') '    1. ALL           (relax cell and ions)'
            write(*, '(a)') '    2. FIX_CELL      (relax ions only)'
            write(*, '(a)') '    Enter choice : ' // trim(prompt_text)

            read(*, '(a)', iostat=ios) input
            if (ios /= 0) then
                iostat = IO_INVALID_INPUT
                write(*, '(a)') '  Invalid input. Please enter a number 1-2.'
                cycle
            end if

            read(input, '(I6)', iostat=ios) choice
            if (ios /= 0) then
                write(*, '(a)') '  Invalid input. Please enter a number 1-2.'
                cycle
            end if

            select case (choice)
            case (1) ; result_mode = CELL_ALL
            case (2) ; result_mode = CELL_INTE
            case default
                write(*, '(a)') '  Invalid choice. Please enter 1-2.'
                cycle
            end select

            exit
        end do
    end subroutine ask_cell_opt_mode

    subroutine ask_symmetry_source(prompt_text, result_src, iostat)
        !! Simplified symmetry selection: NONE or AUTO only
        character(len=*), intent(in)  :: prompt_text
        character(len=*), intent(out) :: result_src
        integer, intent(out)           :: iostat
        integer :: choice, ios
        character(len=MAX_LINE_LEN) :: input

        iostat = 0
        result_src = SYM_NONE

        do
            write(*, '(a)') ''
            write(*, '(a)') '  Select symmetry handling:'
            write(*, '(a)') '    1. NONE (P1)'
            write(*, '(a)') '    2. AUTO (Auto detect by CASTEP)'
            write(*, '(a)') '    Enter choice : ' // trim(prompt_text)

            read(*, '(a)', iostat=ios) input
            if (ios /= 0) then
                iostat = IO_INVALID_INPUT
                write(*, '(a)') '  Invalid input. Please enter a number 1-2.'
                cycle
            end if

            read(input, '(I6)', iostat=ios) choice
            if (ios /= 0) then
                write(*, '(a)') '  Invalid input. Please enter a number 1-2.'
                cycle
            end if

            select case (choice)
            case (1) ; result_src = SYM_NONE
            case (2) ; result_src = SYM_AUTO
            case default
                write(*, '(a)') '  Invalid choice. Please enter 1-2.'
                cycle
            end select

            exit
        end do
    end subroutine ask_symmetry_source

    subroutine ask_geom_tolerance(prompt_text, result_tol, iostat)
        character(len=*), intent(in)  :: prompt_text
        character(len=*), intent(out) :: result_tol
        integer, intent(out)           :: iostat
        integer :: choice, ios
        character(len=MAX_LINE_LEN) :: input

        iostat = 0
        result_tol = GEO_MEDIUM

        do
            write(*, '(a)') ''
            write(*, '(a)') '  Select geometry optimization tolerance : '
            write(*, '(a)') '    1. COARSE'
            write(*, '(a)') '    2. MEDIUM'
            write(*, '(a)') '    3. FINE'
            write(*, '(a)') '    4. EXTREME (For DFPT Phonon Calculation)'
            write(*, '(a)') '    Enter choice : ' // trim(prompt_text)

            read(*, '(a)', iostat=ios) input
            if (ios /= 0) then
                iostat = IO_INVALID_INPUT
                write(*, '(a)') '  Invalid input. Please enter a number 1-4.'
                cycle
            end if

            read(input, '(I6)', iostat=ios) choice
            if (ios /= 0) then
                write(*, '(a)') '  Invalid input. Please enter a number 1-4.'
                cycle
            end if

            select case (choice)
            case (1) ; result_tol = GEO_COARSE
            case (2) ; result_tol = GEO_MEDIUM
            case (3) ; result_tol = GEO_FINE
            case (4) ; result_tol = GEO_EXTREME
            case default
                write(*, '(a)') '  Invalid choice. Please enter 1-4.'
                cycle
            end select

            exit
        end do
    end subroutine ask_geom_tolerance


    ! ========== Helper functions for menu display ==========

    pure function task_label(task) result(lbl)
        character(len=*), intent(in) :: task
        character(32) :: lbl

        select case (trim(task))
        case (TASK_ENERGY)              ; lbl = 'Energy'
        case (TASK_ELECTRONIC_SPECTRO)  ; lbl = 'ElectronicSpectroscopy'
        case (TASK_GEOMETRY_OPT)        ; lbl = 'GeometryOptimisation'
        case (TASK_MOLECULAR_DYN)       ; lbl = 'MolecularDynamics'
        case (TASK_TRANSITION_STATE)    ; lbl = 'TransitionState'
        case (TASK_PHONON)              ; lbl = 'Phonon'
        case (TASK_EFIELD)              ; lbl = 'Efield'
        case (TASK_PHONON_EFIELD)       ; lbl = 'Phonon+Efield'
        case (TASK_THERMODYNAMICS)      ; lbl = 'Thermodynamics'
        case (TASK_MAGRES)              ; lbl = 'MagneticResponse'
        case (TASK_SPECTRAL)            ; lbl = 'Spectral'
        case (TASK_EPCOUPLING)          ; lbl = 'EpCoupling'
        case (TASK_GENETIC_ALGO)        ; lbl = 'GeneticAlgorithm'
        case (TASK_SOCKET_DRIVER)       ; lbl = 'SocketDriver'
        case (TASK_ELASTIC)             ; lbl = 'Elastic'
        case (TASK_AUTOSOLVATION)       ; lbl = 'Autosolvation'
        case default                    ; lbl = trim(task)
        end select
    end function task_label

    pure function task_short_label(task) result(lbl)
        character(len=*), intent(in) :: task
        character(len=16) :: lbl
        select case (trim(task))
        case (TASK_ENERGY)             ; lbl = 'Energy'
        case (TASK_GEOMETRY_OPT)       ; lbl = 'GeomOpt'
        case (TASK_ELECTRONIC_SPECTRO) ; lbl = 'ElecSpectro'
        case (TASK_PHONON)             ; lbl = 'Phonon'
        case (TASK_PHONON_EFIELD)      ; lbl = 'PhononEfield'
        case (TASK_EFIELD)             ; lbl = 'Efield'
        case (TASK_THERMODYNAMICS)     ; lbl = 'Thermo'
        case (TASK_MOLECULAR_DYN)      ; lbl = 'MolDyn'
        case (TASK_TRANSITION_STATE)   ; lbl = 'TransState'
        case (TASK_MAGRES)             ; lbl = 'MagRes'
        case (TASK_SPECTRAL)           ; lbl = 'Spectral'
        case (TASK_EPCOUPLING)         ; lbl = 'EpCoupling'
        case (TASK_GENETIC_ALGO)       ; lbl = 'GeneticAlgo'
        case (TASK_SOCKET_DRIVER)      ; lbl = 'SocketDrv'
        case (TASK_ELASTIC)            ; lbl = 'Elastic'
        case (TASK_AUTOSOLVATION)      ; lbl = 'Autosolv'
        case default                   ; lbl = trim(task)
        end select
    end function task_short_label

    pure function cutoff_label(cutoff) result(lbl)
        real(dp), intent(in) :: cutoff
        character(16) :: lbl
        integer :: iv

        iv = int(cutoff)
        if (cutoff - dble(iv) >= -1.0e-10_dp .and. cutoff - dble(iv) <= 1.0e-10_dp) then
            write(lbl, '(I6)') iv
        else
            write(lbl, '(F8.1)') cutoff
        end if
        lbl = adjustl(lbl)
    end function cutoff_label

    pure function kpoint_label(scheme) result(lbl)
        character(len=*), intent(in) :: scheme
        character(32) :: lbl

        select case (trim(scheme))
        case (KPOINT_GAMMA)              ; lbl = 'GAMMA'
        case (KPOINT_MONKHORST_PACK)     ; lbl = 'MONKHORST_PACK'
        case default                     ; lbl = trim(scheme)
        end select
    end function kpoint_label

    pure function qpoint_label(scheme, grid) result(lbl)
        character(len=*), intent(in) :: scheme
        integer, intent(in) :: grid(3)
        character(16) :: lbl

        if (trim(scheme) == PHONON_QPOINT_MP_GRID) then
            write(lbl, '(i0, 1x, i0, 1x, i0)') grid(1), grid(2), grid(3)
        else
            lbl = trim(scheme)
        end if
    end function qpoint_label

    pure function scf_label(tol) result(lbl)
        character(len=*), intent(in) :: tol
        character(32) :: lbl
        lbl = adjustl(tol)
    end function scf_label

    pure function geom_tol_label(tol) result(lbl)
        character(len=*), intent(in) :: tol
        character(16) :: lbl

        select case (trim(tol))
        case (GEO_COARSE)  ; lbl = 'COARSE'
        case (GEO_MEDIUM)  ; lbl = 'MEDIUM'
        case (GEO_FINE)    ; lbl = 'FINE'
        case (GEO_EXTREME) ; lbl = 'EXTREME'
        case default       ; lbl = trim(tol)
        end select
    end function geom_tol_label

    pure function sym_label(src) result(lbl)
        character(len=*), intent(in) :: src
        character(32) :: lbl

        select case (trim(src))
        case (SYM_AUTO) ; lbl = 'AUTO'
        case (SYM_NONE) ; lbl = 'NONE'
        case default     ; lbl = trim(src)
        end select
    end function sym_label

    pure function sp_label(spin) result(lbl)
        logical, intent(in) :: spin
        character(8) :: lbl

        if (spin) then
            lbl = 'true'
        else
            lbl = 'false'
        end if
        lbl = adjustl(lbl)
    end function sp_label

    pure function smearing_label(s) result(lbl)
        logical, intent(in) :: s
        character(8) :: lbl

        if (s) then
            lbl = 'on'
        else
            lbl = 'off'
        end if
        lbl = adjustl(lbl)
    end function smearing_label

    subroutine parse_grid_comma(s, grid, ios)
        !! Parse "x,y,z" or "x, y, z" into grid(3)
        character(len=*), intent(in)  :: s
        integer, intent(out)          :: grid(3)
        integer, intent(out)          :: ios
        character(len=len(s)) :: cleaned, field
        integer :: i, n_comma, count, start_pos

        ios = 1
        grid = 0
        cleaned = trim(adjustl(s))
        n_comma = 0
        do i = 1, len_trim(cleaned)
            if (cleaned(i:i) == ',') n_comma = n_comma + 1
        end do
        if (n_comma /= 2) return

        count = 0
        start_pos = 1
        do i = 1, len_trim(cleaned) + 1
            if (i > len_trim(cleaned) .or. cleaned(i:i) == ',') then
                count = count + 1
                if (count >= 4) exit
                field = adjustl(cleaned(start_pos:i-1))
                read(field, '(I12)', iostat=ios) grid(count)
                if (ios /= 0 .or. grid(count) < 1) return
                start_pos = i + 1
            end if
        end do
    end subroutine parse_grid_comma

    subroutine ask_advanced_options(prompt_text, result_smearing, result_max_scf, result_conv_win, &
         result_calculate_elf, result_calculate_edd, iostat, &
         result_phonon_calc_dos, result_phonon_dos_spacing, result_phonon_sum_rule, &
         result_phonon_finite_disp, result_phonon_max_cycles, result_phonon_dfpt_method, &
         result_phonon_write_fc, result_phonon_write_dyn, result_phonon_lo_to, &
         result_phonon_fc_cutoff, result_phonon_fc_cutoff_method, result_phonon_dos_limit, &
         result_phonon_max_cg, result_phonon_kpt_sym, &
         result_calc_born, result_calc_raman, result_raman_method, result_efield_dfpt, &
         result_efield_max_cycles, result_efield_energy_tol, result_efield_conv_win, &
         result_efield_freq_spacing, result_efield_osc_q, result_efield_ion_perm, &
         result_efield_ignore_molec, result_is_dfpt, result_is_efield)
        character(len=*), intent(in)  :: prompt_text
        logical, intent(inout)        :: result_smearing
        integer, intent(inout)        :: result_max_scf
        integer, intent(inout)        :: result_conv_win
        logical, intent(inout)        :: result_calculate_elf
        logical, intent(inout)        :: result_calculate_edd
        integer, intent(out)          :: iostat
        logical, intent(inout), optional :: result_phonon_calc_dos
        real(dp), intent(inout), optional :: result_phonon_dos_spacing
        character(len=*), intent(inout), optional :: result_phonon_sum_rule
        real(dp), intent(inout), optional :: result_phonon_finite_disp
        integer,  intent(inout), optional :: result_phonon_max_cycles
        character(len=*), intent(inout), optional :: result_phonon_dfpt_method
        logical, intent(inout), optional :: result_phonon_write_fc
        logical, intent(inout), optional :: result_phonon_write_dyn
        logical, intent(inout), optional :: result_phonon_lo_to
        real(dp), intent(inout), optional :: result_phonon_fc_cutoff
        character(len=*), intent(inout), optional :: result_phonon_fc_cutoff_method
        real(dp), intent(inout), optional :: result_phonon_dos_limit
        integer,  intent(inout), optional :: result_phonon_max_cg
        logical,  intent(inout), optional :: result_phonon_kpt_sym
        logical,  intent(inout), optional :: result_calc_born
        logical,  intent(inout), optional :: result_calc_raman
        character(len=*), intent(inout), optional :: result_raman_method
        character(len=*), intent(inout), optional :: result_efield_dfpt
        integer,  intent(inout), optional :: result_efield_max_cycles
        real(dp), intent(inout), optional :: result_efield_energy_tol
        integer,  intent(inout), optional :: result_efield_conv_win
        real(dp), intent(inout), optional :: result_efield_freq_spacing
        real(dp), intent(inout), optional :: result_efield_osc_q
        logical,  intent(inout), optional :: result_efield_ion_perm
        character(len=16), intent(inout), optional :: result_efield_ignore_molec
        logical,  intent(in),    optional :: result_is_dfpt, result_is_efield
        integer :: choice, ios
        real(dp) :: rval
        character(len=MAX_LINE_LEN) :: input
        logical :: has_phonon, is_dfpt, is_efield

        iostat = 0
        has_phonon = present(result_phonon_calc_dos)
        is_dfpt    = .false.; if (present(result_is_dfpt)) is_dfpt = result_is_dfpt
        is_efield  = .false.; if (present(result_is_efield)) is_efield = result_is_efield

        do
            write(*, '(a)') ''
            write(*, '(a)') '  Advanced options:'
            write(*, '(a)') '    1. Smearing                  (' // trim(smearing_label(result_smearing)) // ')'
            write(*, '(a)') '    2. Max SCF cycles            (' // trim(int2str(result_max_scf))       // ')'
            write(*, '(a)') '    3. Conv. window              (' // trim(int2str(result_conv_win))    // ')'
            write(*, '(a)') '    4. Calculate ELF             (' // trim(smearing_label(result_calculate_elf)) // ')'
            write(*, '(a)') '    5. Calculate EDD             (' // trim(smearing_label(result_calculate_edd)) // ')'
            if (has_phonon) then
                write(*, '(a)') '    6. Calculate phonon DOS      (' // trim(smearing_label(result_phonon_calc_dos)) // ')'
                if (result_phonon_calc_dos) then
                    write(*, '(a, f6.1, a)') '    7. Phonon DOS spacing        (' , result_phonon_dos_spacing, ' cm-1)'
                    write(*, '(a, f8.1, a)') '    8. Phonon DOS limit          (' , result_phonon_dos_limit, ' cm-1)'
                end if
                write(*, '(a, a, a)')   '    9. Phonon sum rule method    (' , trim(result_phonon_sum_rule), ')'
                if (.not. is_dfpt) &
                    write(*, '(a, f8.3, a)') '   10. Phonon finite disp       (' , result_phonon_finite_disp, ' Bohr)'
                write(*, '(a, i0, a)')  '   11. Phonon max cycles         (' , result_phonon_max_cycles, ')'
                if (is_dfpt) &
                    write(*, '(a, a, a)')   '   12. Phonon DFPT method        (' , trim(result_phonon_dfpt_method), ')'
                write(*, '(a)') '   13. Write force constants    (' // trim(smearing_label(result_phonon_write_fc)) // ')'
                write(*, '(a)') '   14. Write dynamical matrix   (' // trim(smearing_label(result_phonon_write_dyn)) // ')'
                write(*, '(a)') '   15. Calc LO/TO splitting     (' // trim(smearing_label(result_phonon_lo_to)) // ')'
                if (.not. is_dfpt) then
                    if (trim(result_phonon_fc_cutoff_method) == PHONON_CUTOFF_CUMULANT) then
                        write(*, '(a, f8.3, a)') '   16. Force constant cutoff scale (' , result_phonon_fc_cutoff, ')'
                    else
                        write(*, '(a, f8.3, a)') '   16. Force constant cutoff   (' , result_phonon_fc_cutoff, ' Bohr)'
                    end if
                    write(*, '(a, a, a)')   '   17. Fine cutoff method      (' , trim(result_phonon_fc_cutoff_method), ')'
                end if
                write(*, '(a, i0, a)')   '   18. Phonon max CG steps     (' , result_phonon_max_cg, ')'
                write(*, '(a)') '   19. Use k-point symmetry    (' // trim(smearing_label(result_phonon_kpt_sym)) // ')'
                write(*, '(a)') '   20. Calculate Born charges  (' // trim(smearing_label(result_calc_born)) // ')'
                write(*, '(a)') '   21. Calculate Raman         (' // trim(smearing_label(result_calc_raman)) // ')'
                if (result_calc_raman) &
                    write(*, '(a, a, a)')   '   22. Raman method            (' , trim(result_raman_method), ')'
                if (is_efield) then
                    write(*, '(a, a, a)')   '   23. EFIELD DFPT method      (' , trim(result_efield_dfpt), ')'
                    write(*, '(a, i0, a)')   '   24. EFIELD max cycles        (' , result_efield_max_cycles, ')'
                    write(*, '(a, es9.1, a)') '   25. EFIELD energy tol        (' , result_efield_energy_tol, ')'
                    write(*, '(a, i0, a)')   '   26. EFIELD conv. window      (' , result_efield_conv_win, ')'
                    write(*, '(a, f8.1, a)') '   27. EFIELD freq spacing      (' , result_efield_freq_spacing, ' cm-1)'
                    write(*, '(a, f8.1, a)') '   28. EFIELD oscillator Q      (' , result_efield_osc_q, ')'
                    write(*, '(a)') '   29. Calc ion permittivity   (' // trim(smearing_label(result_efield_ion_perm)) // ')'
                    if (present(result_efield_ignore_molec)) &
                        write(*, '(a, a, a)') '   30. EFIELD ignore molec modes (' , trim(result_efield_ignore_molec), ')'
                end if
            end if
            write(*, '(a)') '    0. Back to main menu'
            write(*, '(a)') '    Enter choice : ' // trim(prompt_text)

            read(*, '(a)', iostat=ios) input
            if (ios /= 0 .or. len_trim(input) == 0) exit

            read(input, '(I6)', iostat=ios) choice
            if (ios /= 0) then
                write(*, '(a)') '  Invalid input.'
                cycle
            end if

            select case (choice)
            case (0); exit
            case (1)
                result_smearing = .not. result_smearing
                write(*, '(a)') '  Smearing set to: ' // trim(smearing_label(result_smearing))
            case (2)
                write(*, '(a, a)') '  Enter max SCF cycles [', trim(int2str(result_max_scf)), ']: '
                read(*, '(a)', iostat=ios) input
                if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
                    read(input, *, iostat=ios) choice
                    if (ios == 0 .and. choice > 0) result_max_scf = choice
                end if
            case (3)
                write(*, '(a, a)') '  Enter convergence window [', trim(int2str(result_conv_win)), ']: '
                read(*, '(a)', iostat=ios) input
                if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
                    read(input, *, iostat=ios) choice
                    if (ios == 0 .and. choice >= 2) result_conv_win = choice
                end if
            case (4)
                result_calculate_elf = .not. result_calculate_elf
                write(*, '(a)') '  Calculate ELF set to: ' // trim(smearing_label(result_calculate_elf))
            case (5)
                result_calculate_edd = .not. result_calculate_edd
                write(*, '(a)') '  Calculate EDD set to: ' // trim(smearing_label(result_calculate_edd))
            case (6)
                if (has_phonon) then
                    result_phonon_calc_dos = .not. result_phonon_calc_dos
                    write(*, '(a)') '  Calculate phonon DOS set to: ' // trim(smearing_label(result_phonon_calc_dos))
                end if
            case (7)
                if (has_phonon .and. result_phonon_calc_dos) then
                    write(*, '(a, f6.1, a)') '  Enter phonon DOS spacing [', result_phonon_dos_spacing, ' cm-1]: '
                    read(*, '(a)', iostat=ios) input
                    if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
                        read(input, *, iostat=ios) rval
                        if (ios == 0 .and. rval > 0.0_dp) result_phonon_dos_spacing = rval
                    end if
                end if
            case (8)
                if (has_phonon .and. result_phonon_calc_dos) then
                    write(*, '(a, f8.1, a)') '  Enter phonon DOS limit [', result_phonon_dos_limit, ' cm-1]: '
                    read(*, '(a)', iostat=ios) input
                    if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
                        read(input, *, iostat=ios) rval
                        if (ios == 0 .and. rval > 0.0_dp) result_phonon_dos_limit = rval
                    end if
                end if
            case (9)
                if (has_phonon) then
                    call ask_phonon_sum_rule(result_phonon_sum_rule, ios)
                end if
            case (10)
                if (has_phonon .and. .not. is_dfpt) then
                    write(*, '(a, f8.3, a)') '  Enter finite displacement [', result_phonon_finite_disp, ' Bohr]: '
                    read(*, '(a)', iostat=ios) input
                    if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
                        read(input, *, iostat=ios) rval
                        if (ios == 0 .and. rval > 0.0_dp) result_phonon_finite_disp = rval
                    end if
                end if
            case (11)
                if (has_phonon) then
                    write(*, '(a, i0, a)') '  Enter max cycles [', result_phonon_max_cycles, ']: '
                    read(*, '(a)', iostat=ios) input
                    if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
                        read(input, *, iostat=ios) choice
                        if (ios == 0 .and. choice > 0) result_phonon_max_cycles = choice
                    end if
                end if
            case (12)
                if (has_phonon .and. is_dfpt) then
                    call ask_phonon_dfpt_method(result_phonon_dfpt_method, ios)
                end if
            case (13)
                if (has_phonon) result_phonon_write_fc = .not. result_phonon_write_fc
            case (14)
                if (has_phonon) result_phonon_write_dyn = .not. result_phonon_write_dyn
            case (15)
                if (has_phonon .and. .not. is_efield) then
                    result_phonon_lo_to = .not. result_phonon_lo_to
                else if (has_phonon .and. is_efield) then
                    write(*, '(a)') '  LO/TO splitting is required for EFIELD tasks and cannot be disabled.'
                end if
            case (16)
                if (has_phonon .and. .not. is_dfpt) then
                    if (trim(result_phonon_fc_cutoff_method) == PHONON_CUTOFF_CUMULANT) then
                        write(*, '(a, f8.3, a)') '  Enter force constant cutoff scale [', result_phonon_fc_cutoff, ']: '
                    else
                        write(*, '(a, f8.3, a)') '  Enter force constant cutoff [', result_phonon_fc_cutoff, ' Bohr]: '
                    end if
                    read(*, '(a)', iostat=ios) input
                    if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
                        read(input, *, iostat=ios) rval
                        if (ios == 0 .and. rval >= 0.0_dp) result_phonon_fc_cutoff = rval
                    end if
                end if
            case (17)
                if (has_phonon .and. .not. is_dfpt) then
                    write(*, '(a)') '  Fine cutoff method:  1. CUMULANT  2. SPHERICAL'
                    write(*, '(a)', advance='no') '  Enter choice: '
                    read(*, '(a)', iostat=ios) input
                    if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
                        read(input, *, iostat=ios) choice
                        if (choice == 1) result_phonon_fc_cutoff_method = PHONON_CUTOFF_CUMULANT
                        if (choice == 2) result_phonon_fc_cutoff_method = PHONON_CUTOFF_SPHERICAL
                    end if
                end if
            case (18)
                if (has_phonon) then
                    write(*, '(a, i0, a)') '  Enter max CG steps [', result_phonon_max_cg, ']: '
                    read(*, '(a)', iostat=ios) input
                    if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
                        read(input, *, iostat=ios) choice
                        if (ios == 0 .and. choice >= 0) result_phonon_max_cg = choice
                    end if
                end if
            case (19); if (has_phonon) result_phonon_kpt_sym = .not. result_phonon_kpt_sym
            case (20); if (has_phonon) result_calc_born = .not. result_calc_born
            case (21); if (has_phonon) result_calc_raman = .not. result_calc_raman
            case (22)
                if (has_phonon .and. result_calc_raman) then
                    write(*, '(a)') '  Raman method:  1. DFPT  2. FINITEDISPLACEMENT'
                    write(*, '(a)', advance='no') '  Enter choice: '
                    read(*, '(a)', iostat=ios) input
                    if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
                        read(input, *, iostat=ios) choice
                        if (choice == 1) result_raman_method = 'DFPT'
                        if (choice == 2) result_raman_method = 'FINITEDISPLACEMENT'
                    end if
                end if
            case (23)
                if (has_phonon .and. is_efield) then
                    write(*, '(a)') '  EFIELD DFPT method:  1. ALLBANDS  2. DM'
                    write(*, '(a)', advance='no') '  Enter choice: '
                    read(*, '(a)', iostat=ios) input
                    if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
                        read(input, *, iostat=ios) choice
                        if (choice == 1) result_efield_dfpt = 'ALLBANDS'
                        if (choice == 2) result_efield_dfpt = 'DM'
                    end if
                end if
            case (24)
                if (has_phonon .and. is_efield) then
                    write(*, '(a, i0, a)') '  Enter EFIELD max cycles [', result_efield_max_cycles, ']: '
                    read(*, '(a)', iostat=ios) input
                    if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
                        read(input, *, iostat=ios) choice
                        if (ios == 0 .and. choice > 0) result_efield_max_cycles = choice
                    end if
                end if
            case (25)
                if (has_phonon .and. is_efield) then
                    write(*, '(a, es9.1, a)') '  Enter EFIELD energy tol [', result_efield_energy_tol, ']: '
                    read(*, '(a)', iostat=ios) input
                    if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
                        read(input, *, iostat=ios) rval
                        if (ios == 0 .and. rval > 0.0_dp) result_efield_energy_tol = rval
                    end if
                end if
            case (26)
                if (has_phonon .and. is_efield) then
                    write(*, '(a, i0, a)') '  Enter EFIELD conv. window [', result_efield_conv_win, ']: '
                    read(*, '(a)', iostat=ios) input
                    if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
                        read(input, *, iostat=ios) choice
                        if (ios == 0 .and. choice >= 2) result_efield_conv_win = choice
                    end if
                end if
            case (27)
                if (has_phonon .and. is_efield) then
                    write(*, '(a, f8.1, a)') '  Enter EFIELD freq spacing [', result_efield_freq_spacing, ' cm-1]: '
                    read(*, '(a)', iostat=ios) input
                    if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
                        read(input, *, iostat=ios) rval
                        if (ios == 0 .and. rval > 0.0_dp) result_efield_freq_spacing = rval
                    end if
                end if
            case (28)
                if (has_phonon .and. is_efield) then
                    write(*, '(a, f8.1, a)') '  Enter EFIELD oscillator Q [', result_efield_osc_q, ']: '
                    read(*, '(a)', iostat=ios) input
                    if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
                        read(input, *, iostat=ios) rval
                        if (ios == 0 .and. rval > 0.0_dp) result_efield_osc_q = rval
                    end if
                end if
            case (29)
                if (has_phonon .and. is_efield) &
                    result_efield_ion_perm = .not. result_efield_ion_perm
            case (30)
                if (has_phonon .and. is_efield) then
                    write(*, '(a)') '  Ignore molec modes:  1. CRYSTAL(3)  2. MOLECULE(6)  3. LINEAR_MOLECULE(5)'
                    write(*, '(a)', advance='no') '  Enter choice: '
                    read(*, '(a)', iostat=ios) input
                    if (ios == 0 .and. len_trim(adjustl(input)) > 0) then
                        read(input, *, iostat=ios) choice
                        if (choice == 1) result_efield_ignore_molec = 'CRYSTAL'
                        if (choice == 2) result_efield_ignore_molec = 'MOLECULE'
                        if (choice == 3) result_efield_ignore_molec = 'LINEAR_MOLECULE'
                    end if
                end if
            case default
                write(*, '(a)') '  Invalid choice.'
            end select
        end do
    end subroutine ask_advanced_options

    subroutine ask_spectral_task(prompt_text, result_spectral, iostat)
        !! Spectral task sub-menu: BandStructure or BandStructure_pDOS
        character(len=*), intent(in)  :: prompt_text
        character(len=*), intent(inout) :: result_spectral
        integer, intent(out)           :: iostat
        integer :: choice, ios
        character(len=MAX_LINE_LEN) :: input

        iostat = 0
        result_spectral = 'BandStructure'

        do
            write(*, '(a)') ''
            write(*, '(a)') '  Select spectral task type:'
            write(*, '(a)') '    1. Band Structure and DOS'
            write(*, '(a)') '    2. Band Structure and pDOS'
            write(*, '(a)') '    Enter choice : ' // trim(prompt_text)

            read(*, '(a)', iostat=ios) input
            if (ios /= 0 .or. len_trim(input) == 0) then
                ! Empty or invalid input → accept default
                exit
            end if

            read(input, '(I6)', iostat=ios) choice
            if (ios /= 0) then
                write(*, '(a)') '  Invalid input. Please enter a number 1-2.'
                cycle
            end if

            select case (choice)
            case (1)
                result_spectral = 'BandStructure'
                exit
            case (2)
                result_spectral = 'BandStructure_pDOS'
                exit
            case default
                write(*, '(a)') '  Invalid choice. Please enter 1-2.'
                cycle
            end select
        end do
    end subroutine ask_spectral_task

    subroutine ask_phonon_qpoint_scheme(cfg, iostat)
        type(castep_config_t), intent(inout) :: cfg
        integer,               intent(out)   :: iostat
        call ask_phonon_qpoint_scheme_impl('Phonon q-point scheme', 'k-point path', &
            cfg%phonon_qpoint_scheme, cfg%phonon_kpoint_mp_grid, &
            cfg%phonon_kpoint_path, cfg%phonon_kpoint_path_spacing, iostat)
    end subroutine ask_phonon_qpoint_scheme

    subroutine ask_phonon_fine_qpoint_scheme(cfg, iostat)
        type(castep_config_t), intent(inout) :: cfg
        integer,               intent(out)   :: iostat
        call ask_phonon_qpoint_scheme_impl('Phonon fine q-point scheme', 'fine k-point path', &
            cfg%phonon_fine_qpoint_scheme, cfg%phonon_fine_kpoint_mp_grid, &
            cfg%phonon_fine_kpoint_path, cfg%phonon_fine_kpoint_path_spacing, iostat)
    end subroutine ask_phonon_fine_qpoint_scheme

    subroutine ask_phonon_qpoint_scheme_impl(title, path_prompt, scheme_name, &
            mp_grid, path_str, path_spacing, iostat)
        character(len=*), intent(in)    :: title, path_prompt
        character(len=*), intent(inout) :: scheme_name, path_str
        integer,          intent(inout) :: mp_grid(3)
        real(dp),         intent(inout) :: path_spacing
        integer,          intent(out)   :: iostat
        integer :: choice, ios
        character(len=MAX_LINE_LEN) :: input

        iostat = 0
        write(*, '(a)') ''
        write(*, '(a)') '  ---------------------------------'
        write(*, '(a)') '    ' // trim(title)
        write(*, '(a)') '  ---------------------------------'
        write(*, '(a)') '    1. MP_GRID'
        write(*, '(a)') '    2. PATH   (advanced -- band structure)'
        write(*, '(a)', advance='no') '    Enter choice: '
        read(*, '(a)', iostat=ios) input
        if (ios /= 0) then; iostat = 1; return; end if
        read(input, *, iostat=ios) choice
        if (ios /= 0) return
        select case (choice)
        case (1)
            scheme_name = PHONON_QPOINT_MP_GRID
            call ask_phonon_kpoint_grid('  MP grid (e.g. 4 4 4): ', mp_grid, iostat)
        case (2)
            scheme_name = PHONON_QPOINT_PATH
            call ask_phonon_path_impl(path_prompt, path_str, path_spacing, iostat)
        end select
    end subroutine ask_phonon_qpoint_scheme_impl

    subroutine ask_phonon_kpoint_grid(prompt_text, grid, iostat)
        character(len=*), intent(in)    :: prompt_text
        integer,          intent(inout) :: grid(3)
        integer,          intent(out)   :: iostat
        character(len=MAX_LINE_LEN) :: input
        integer :: ios

        iostat = 0
        do
            write(*, '(a)', advance='no') trim(prompt_text)
            read(*, '(a)', iostat=ios) input
            if (ios /= 0) then; iostat = 1; return; end if
            if (len_trim(input) == 0) return  ! keep current value
            read(input, *, iostat=ios) grid(1), grid(2), grid(3)
            if (ios == 0 .and. all(grid >= 1)) return
            write(*, '(a)') '  Invalid. Use format: n1 n2 n3 (all >= 1). Try again.'
        end do
    end subroutine ask_phonon_kpoint_grid

    subroutine ask_phonon_path_impl(prompt_prefix, path_str, path_spacing, iostat)
        character(len=*), intent(in)    :: prompt_prefix
        character(len=*), intent(inout) :: path_str
        real(dp),         intent(inout) :: path_spacing
        integer,          intent(out)   :: iostat
        character(len=MAX_LINE_LEN) :: line
        integer :: ios, n
        real(dp) :: spacing

        iostat = 0
        write(*, '(a)') '  Enter ' // trim(prompt_prefix) // ' (empty line to finish):'
        write(*, '(a)') '  Format: kx ky kz [! label]'
        write(*, '(a)') '  Use "BREAK" alone on a line to separate segments.'
        path_str = ''
        n = 0
        do
            write(*, '(a,i0,a)', advance='no') '  [', n+1, '] '
            read(*, '(a)', iostat=ios) line
            if (ios /= 0 .or. len_trim(line) == 0) exit
            if (n > 0) path_str = trim(path_str) // achar(10)
            path_str = trim(path_str) // trim(line)
            n = n + 1
        end do
        if (n == 0) return

        spacing = path_spacing
        write(*, '(a,f6.3,a)', advance='no') '  Path spacing [', spacing, ' 1/A]: '
        read(*, '(a)', iostat=ios) line
        if (ios == 0 .and. len_trim(line) > 0) then
            read(line, *, iostat=ios) spacing
            if (ios == 0 .and. spacing > 0.0_dp) path_spacing = spacing
        end if
    end subroutine ask_phonon_path_impl

    subroutine ask_phonon_method(prompt_text, method, iostat)
        character(len=*), intent(in)    :: prompt_text
        character(len=*), intent(inout) :: method
        integer,          intent(out)   :: iostat
        integer :: choice, ios
        character(len=MAX_LINE_LEN) :: input

        iostat = 0
        write(*, '(a)') trim(prompt_text)
        write(*, '(a)') '  1. DFPT (Linear Response)'
        write(*, '(a)') '  2. FINITEDISPLACEMENT'
        write(*, '(a)', advance='no') '  Enter choice: '
        read(*, '(a)', iostat=ios) input
        if (ios /= 0) then; iostat = 1; return; end if
        read(input, *, iostat=ios) choice
        if (ios == 0) then
            select case (choice)
            case (1); method = PHONON_METHOD_DFPT
            case (2); method = PHONON_METHOD_FD
            end select
        end if
    end subroutine ask_phonon_method

    subroutine ask_phonon_energy_tol(prompt_text, tol, iostat)
        character(len=*), intent(in)    :: prompt_text
        real(dp),         intent(inout) :: tol
        integer,          intent(out)   :: iostat
        real(dp) :: val
        integer :: ios
        character(len=MAX_LINE_LEN) :: input

        iostat = 0
        write(*, '(a, es9.1, a)', advance='no') trim(prompt_text) // ' [', tol, ']: '
        read(*, '(a)', iostat=ios) input
        if (ios /= 0) then; iostat = 1; return; end if
        if (len_trim(input) > 0) then
            read(input, *, iostat=ios) val
            if (ios == 0 .and. val > 0.0_dp) tol = val
        end if
    end subroutine ask_phonon_energy_tol

    subroutine ask_phonon_fine_method(prompt_text, method, current_method, iostat)
        character(len=*), intent(in)    :: prompt_text
        character(len=*), intent(inout) :: method
        character(len=*), intent(in)    :: current_method
        integer,          intent(out)   :: iostat
        integer :: choice, ios
        character(len=MAX_LINE_LEN) :: input
        logical :: is_dfpt

        is_dfpt = (trim(current_method) == PHONON_METHOD_DFPT)
        iostat = 0
        write(*, '(a)') trim(prompt_text)
        write(*, '(a)') '  1. INTERPOLATE'
        if (.not. is_dfpt) then
            write(*, '(a)') '  2. SUPERCELL (incompatible with DFPT)'
        end if
        write(*, '(a)', advance='no') '  Enter choice: '
        read(*, '(a)', iostat=ios) input
        if (ios /= 0) then; iostat = 1; return; end if
        read(input, *, iostat=ios) choice
        if (ios /= 0) return
        if (is_dfpt) then
            select case (choice)
            case (1); method = PHONON_FINE_INTERPOLATE
            end select
        else
            select case (choice)
            case (1); method = PHONON_FINE_INTERPOLATE
            case (2); method = PHONON_FINE_SUPERCELL
            end select
        end if
    end subroutine ask_phonon_fine_method

    subroutine ask_phonon_dfpt_method(method, iostat)
        character(len=*), intent(inout) :: method
        integer,          intent(out)   :: iostat
        integer :: choice, ios
        character(len=MAX_LINE_LEN) :: input

        iostat = 0
        write(*, '(a)') '  DFPT method:'
        write(*, '(a)') '    1. DM (Baroni-Green)'
        write(*, '(a)') '    2. ALLBANDS (Gonze variational)'
        write(*, '(a)', advance='no') '  Enter choice: '
        read(*, '(a)', iostat=ios) input
        if (ios /= 0) then; iostat = 1; return; end if
        read(input, *, iostat=ios) choice
        if (ios /= 0) return
        select case (choice)
        case (1); method = PHONON_DFPT_DM
        case (2); method = PHONON_DFPT_ALLBANDS
        end select
    end subroutine ask_phonon_dfpt_method

    subroutine ask_phonon_sum_rule(method, iostat)
        character(len=*), intent(inout) :: method
        integer,          intent(out)   :: iostat
        integer :: choice, ios
        character(len=MAX_LINE_LEN) :: input

        iostat = 0
        write(*, '(a)') '  Phonon sum rule methods:'
        write(*, '(a)') '   0. NONE'
        write(*, '(a)') '   1. RECIPROCAL'
        write(*, '(a)') '   2. REALSPACE'
        write(*, '(a)') '   3. REAL-RECIP'
        write(*, '(a)') '   4. MOLECULAR'
        write(*, '(a)', advance='no') '  Enter choice: '
        read(*, '(a)', iostat=ios) input
        if (ios /= 0) return
        read(input, *, iostat=ios) choice
        if (ios /= 0) return
        select case (choice)
        case (0); method = PHONON_SUM_NONE
        case (1); method = PHONON_SUM_RECIPROCAL
        case (2); method = PHONON_SUM_REALSPACE
        case (3); method = PHONON_SUM_REAL_RECIP
        case (4); method = PHONON_SUM_MOLECULAR
        end select
    end subroutine ask_phonon_sum_rule

    subroutine ask_phonon_supercell_matrix(prompt_text, matrix, iostat)
        character(len=*), intent(in)    :: prompt_text
        integer,          intent(inout) :: matrix(3,3)
        integer,          intent(out)   :: iostat
        integer :: i, ios
        character(len=MAX_LINE_LEN) :: input

        iostat = 0
        write(*, '(a)') trim(prompt_text)
        write(*, '(a)') '  Enter 3 rows (e.g. "2 0 0"):'
        do i = 1, 3
            write(*, '(a,i0,a)', advance='no') '  Row ', i, ': '
            read(*, '(a)', iostat=ios) input
            if (ios /= 0) then; iostat = 1; return; end if
            if (len_trim(input) > 0) then
                read(input, *, iostat=ios) matrix(i, 1), matrix(i, 2), matrix(i, 3)
                if (ios /= 0) then
                    write(*, '(a)') '  Invalid. Use format: n1 n2 n3'
                    matrix(i,:) = [1,1,1]  ! on error for this row, use 1s
                end if
            end if
        end do
    end subroutine ask_phonon_supercell_matrix

    pure function auto_output_name(path, task) result(name)
        character(len=*), intent(in) :: path, task
        character(len=512) :: name, stem
        integer :: n

        stem = trim(path)
        n = len_trim(stem)
        do while (n > 0)
            if (stem(n:n) == '/') then; stem = stem(n+1:); exit; end if
            n = n - 1
        end do
        n = len_trim(stem)
        do while (n > 0)
            if (stem(n:n) == '.') then; stem = stem(1:n-1); exit; end if
            n = n - 1
        end do
        name = trim(stem) // '_' // trim(task_short_label(task))
    end function auto_output_name

end module cli_menu
