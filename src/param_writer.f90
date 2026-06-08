module param_writer
    !! CASTEP .param file generator using simple key-value format
    !! Structure: task line first, then common keywords, then task-diff keywords.
    !! No %BLOCK or section wrapping.
    use castep_config, only: &
         dp, TASK_ENERGY, TASK_GEOMETRY_OPT, &
         TASK_MOLECULAR_DYN, TASK_PHONON, TASK_ELECTRONIC_SPECTRO, &
         TASK_TRANSITION_STATE, TASK_EFIELD, TASK_PHONON_EFIELD, TASK_THERMODYNAMICS, &
         TASK_MAGRES, TASK_SPECTRAL, TASK_EPCOUPLING, TASK_GENETIC_ALGO, &
         TASK_SOCKET_DRIVER, TASK_ELASTIC, TASK_AUTOSOLVATION, &
         TOL_SUPERFINE, TOL_FINE, TOL_NORMAL, TOL_COARSE, &
         GEO_COARSE, GEO_MEDIUM, GEO_FINE, GEO_EXTREME, &
         CELL_ALL, CELL_INTE, &
         PHONON_METHOD_DFPT, PHONON_METHOD_FD, &
         PHONON_FINE_NONE, PHONON_FINE_SUPERCELL, PHONON_FINE_INTERPOLATE, &
         PHONON_DFPT_DM, PHONON_DFPT_ALLBANDS, &
         PHONON_CUTOFF_CUMULANT, PHONON_CUTOFF_SPHERICAL, &
         PHONON_SUM_NONE, PHONON_SUM_RECIPROCAL, PHONON_SUM_REALSPACE, &
         PHONON_SUM_REAL_RECIP, PHONON_SUM_MOLECULAR, &
         OPT_BFGS, OPT_LBFGS, OPT_CG, &
         VDW_NONE, PSEUDO_SOC19, IO_WRITE_FAIL, &
         castep_config_t, get_castep_task_name, int2str
    implicit none
    private

    public :: write_param_file

contains

    subroutine write_kv(unit, key, value)
        !! Write a single key : value line
        integer, intent(in) :: unit
        character(len=*), intent(in) :: key
        character(len=*), intent(in) :: value
        write(unit, '(a, a, a)') trim(key), ' : ', trim(value)
    end subroutine write_kv

    subroutine write_param_file(filename, cfg, iostat, iomsg)
        !! Write a CASTEP .param file in key-value format
        character(len=*), intent(in) :: filename
        type(castep_config_t), intent(in) :: cfg
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg
        integer :: unit, ios
        character(len=64) :: castep_task

        iostat = 0

        open(newunit=unit, file=trim(filename), status='unknown', action='write', iostat=ios)
        if (ios /= 0) then
            iostat = IO_WRITE_FAIL
            if (present(iomsg)) iomsg = 'Cannot write file: ' // trim(filename)
            return
        end if

        ! Determine CASTEP task name
        castep_task = get_castep_task_name(cfg%task_type)

        ! First line: task
        write(unit, '(a)') 'task : ' // trim(castep_task)

        ! Common keywords
        call write_kv(unit, 'comment', 'CASTEP calculation from CASTEP Suite')
        call write_kv(unit, 'xc_functional', trim(cfg%xc_functional))
        if (cfg%spin_polarized) then
            call write_kv(unit, 'spin_polarized', 'true')
            call write_kv(unit, 'spin', '0          #need to set by user')
            if (trim(cfg%pseudopotential) == PSEUDO_SOC19) then
                call write_kv(unit, 'spin_orbit_coupling', 'true')
            end if
        else
            call write_kv(unit, 'spin_polarized', 'false')
        end if
        call write_kv(unit, 'opt_strategy', 'speed')
        call write_cutoff_energy(unit, cfg%cutoff_energy)
        call write_kv(unit, 'grid_scale', '2')
        call write_kv(unit, 'fine_grid_scale', '3')
        call write_kv(unit, 'finite_basis_corr', '2')
        call write_kv(unit, 'finite_basis_npoints', '3')
        call write_kv(unit, 'elec_energy_tol', energy_tol_str(cfg%scf_tolerance))
        call write_kv(unit, 'max_scf_cycles', trim(int2str(cfg%max_scf_cycles)))
        if (cfg%smearing) then
            call write_kv(unit, 'fix_occupancy', 'false')
            call write_kv(unit, 'perc_extra_bands', '50            #50% of the total valence electrons')
            call write_kv(unit, 'smearing_width', '0.1')
        else
            call write_kv(unit, 'fix_occupancy', 'true')
        end if
        call write_kv(unit, 'metals_method', 'dm          #density mixing, also could be EDFT')
        call write_kv(unit, 'mixing_scheme', 'Pulay           #also could be Broyden')
        call write_kv(unit, 'mix_charge_amp', '0.5')
        call write_kv(unit, 'mix_charge_gmax', '1.5')
        call write_kv(unit, 'mix_history_length', '20')
        call write_kv(unit, 'elec_convergence_win', trim(int2str(cfg%elec_convergence_win)))
        if (cfg%calculate_elf) then
            call write_kv(unit, 'calculate_ELF', 'true')
        end if
        if (cfg%calculate_edd) then
            call write_kv(unit, 'calculate_densdiff', 'true')
        end if

        ! vdW-DED keywords
        if (trim(cfg%vdw_method) /= VDW_NONE) then
            call write_kv(unit, 'sedc_apply', 'true')
            call write_kv(unit, 'sedc_scheme', trim(cfg%vdw_method))
        end if

        ! Task-diff keywords
        select case (trim(cfg%task_type))
        case (TASK_ENERGY)
            call write_kv(unit, 'write_cell_structure', 'false')
        case (TASK_GEOMETRY_OPT)
            call write_kv(unit, 'write_cell_structure', 'true')
            call write_kv(unit, 'geom_method', trim(cfg%optimizer))
            call write_kv(unit, 'geom_max_iter', '500')
            call write_kv(unit, 'calculate_stress', 'true')
            select case (trim(cfg%geom_tolerance))
            case (GEO_COARSE)
                call write_kv(unit, 'geom_energy_tol', '5e-5')
                call write_kv(unit, 'geom_force_tol', '0.1')
                call write_kv(unit, 'geom_stress_tol', '0.2')
                call write_kv(unit, 'geom_disp_tol', '0.005')
             
            case (GEO_MEDIUM)
                call write_kv(unit, 'geom_energy_tol', '2e-5')
                call write_kv(unit, 'geom_force_tol', '0.05')
                call write_kv(unit, 'geom_stress_tol', '0.1')
                call write_kv(unit, 'geom_disp_tol', '0.002')
                
            case (GEO_FINE)
                call write_kv(unit, 'geom_energy_tol', '1e-5')
                call write_kv(unit, 'geom_force_tol', '0.03')
                call write_kv(unit, 'geom_stress_tol', '0.05')
                call write_kv(unit, 'geom_disp_tol', '0.001')
                
            case (GEO_EXTREME)
                call write_kv(unit, 'geom_energy_tol', '5e-6')
                call write_kv(unit, 'geom_force_tol', '0.01')
                call write_kv(unit, 'geom_stress_tol', '0.02')
                call write_kv(unit, 'geom_disp_tol', '5e-4')
                
            case default
                call write_kv(unit, 'geom_energy_tol', '2e-5')
                call write_kv(unit, 'geom_force_tol', '0.05')
                call write_kv(unit, 'geom_stress_tol', '0.1')
                call write_kv(unit, 'geom_disp_tol', '0.002')
                
            end select
            if (trim(cfg%cell_opt_mode) == CELL_ALL) then
                call write_kv(unit, 'geom_modulus_est', '500 GPa')
            end if
        case (TASK_MOLECULAR_DYN)
            call write_kv(unit, 'md_steps', '100')
            call write_kv(unit, 'time_step', '0.5')
            call write_kv(unit, 'thermostat', 'nvt')
        case (TASK_PHONON)
            call write_phonon_params(unit, cfg)
        case (TASK_ELECTRONIC_SPECTRO)
            call write_kv(unit, 'spectral_task', 'BandStructure')
            call write_kv(unit, 'spectral_write_eigenvalues', 'true')
            call write_kv(unit, 'spectral_perc_extra_bands', '30')
            call write_kv(unit, 'spectral_eigenvalue_tol', '1e-5')
            if (trim(cfg%spectral_task_type) == 'BandStructure_pDOS') then
                call write_kv(unit, 'pdos_calculate_weights', 'true')
            else
                call write_kv(unit, 'pdos_calculate_weights', 'false')
            end if
        case (TASK_TRANSITION_STATE)
            call write_cineb_params(unit, cfg)
        case (TASK_EFIELD)
            call write_efield_params(unit, cfg)
        case (TASK_PHONON_EFIELD)
            call write_phonon_params(unit, cfg)
            if (cfg%calculate_born_charges) then
                call write_kv(unit, 'calculate_born_charges', 'true')
                call write_kv(unit, 'born_charge_sum_rule', 'true')
            end if
            call write_efield_params(unit, cfg)
        case (TASK_THERMODYNAMICS)
            ! TASK=THERMODYNAMICS implies thermodynamics; no separate toggle keyword
            call write_kv(unit, 'phonon_method', 'FINITEDISPLACEMENT')
            call write_kv(unit, 'phonon_fine_method', 'SUPERCELL')
            call write_kv(unit, 'phonon_finite_disp', real2str(cfg%phonon_finite_disp))
            ! phonon_finite_disp: unit is implicit (Bohr), no _unit keyword
            if (trim(cfg%phonon_sum_rule_method) /= PHONON_SUM_NONE) &
                call write_kv(unit, 'phonon_sum_rule_method', trim(cfg%phonon_sum_rule_method))
        case (TASK_MAGRES)
            call write_kv(unit, 'magnetic_response', 'true')
            call write_kv(unit, 'field_direction', '0, 0, 1')
        case (TASK_SPECTRAL)
            call write_kv(unit, 'spectral', 'true')
        case (TASK_EPCOUPLING)
            call write_kv(unit, 'ep_coupling', 'true')
        case (TASK_GENETIC_ALGO)
            call write_kv(unit, 'genetic_algo', 'true')
            call write_kv(unit, 'population_size', '10')
        case (TASK_SOCKET_DRIVER)
            call write_kv(unit, 'socket_driver', 'true')
        case (TASK_ELASTIC)
            call write_kv(unit, 'elastic', 'true')
            call write_kv(unit, 'num_deformations', '6')
        case (TASK_AUTOSOLVATION)
            call write_kv(unit, 'autosolvation', 'true')
        case default
            call write_kv(unit, 'geom_opt', 'false')
            call write_kv(unit, 'max_geom_iter', '300')
            call write_kv(unit, 'relax_strain', '.false.')
        end select

        close(unit)

    contains

        subroutine write_cutoff_energy(unit, val)
            integer, intent(in) :: unit
            real(dp), intent(in) :: val
            integer :: iv
            character(20) :: s

            iv = int(val)
            if (val - dble(iv) >= -1.0e-10_dp .and. val - dble(iv) <= 1.0e-10_dp) then
                write(s, '(I6)') iv
            else
                write(s, '(F12.4)') val
            end if
            call write_kv(unit, 'cut_off_energy', adjustl(s))
        end subroutine write_cutoff_energy

    end subroutine write_param_file

    ! ----------------------------------------------------------------
    !  Shared phonon params (used by TASK_PHONON and TASK_PHONON_EFIELD)
    ! ----------------------------------------------------------------
    subroutine write_phonon_params(unit, cfg)
        integer, intent(in) :: unit
        type(castep_config_t), intent(in) :: cfg

        call write_kv(unit, 'phonon_method', trim(cfg%phonon_method))
        if (trim(cfg%phonon_method) == PHONON_METHOD_DFPT) &
            call write_kv(unit, 'phonon_dfpt_method', trim(cfg%phonon_dfpt_method))
        call write_kv(unit, 'phonon_energy_tol', energy_tol_str(cfg%phonon_energy_tol))
        ! phonon_energy_tol: unit is implicit (eV/A^2), no _unit keyword
        call write_kv(unit, 'phonon_max_cycles', int2str(cfg%phonon_max_cycles))
        call write_kv(unit, 'phonon_convergence_win', int2str(cfg%phonon_convergence_win))

        if (trim(cfg%phonon_fine_method) /= PHONON_FINE_NONE) then
            call write_kv(unit, 'phonon_fine_method', trim(cfg%phonon_fine_method))
            if (trim(cfg%phonon_fine_method) == PHONON_FINE_SUPERCELL) &
                call write_kv(unit, 'phonon_finite_disp', real2str(cfg%phonon_finite_disp))
        end if

        if (trim(cfg%phonon_sum_rule_method) /= PHONON_SUM_NONE) &
            call write_kv(unit, 'phonon_sum_rule_method', trim(cfg%phonon_sum_rule_method))

        if (cfg%phonon_calculate_dos) then
            call write_kv(unit, 'phonon_calculate_dos', 'true')
            call write_kv(unit, 'phonon_dos_spacing', real2str(cfg%phonon_dos_spacing))
            call write_kv(unit, 'phonon_dos_spacing_unit', 'cm-1')
            call write_kv(unit, 'phonon_dos_limit', real2str(cfg%phonon_dos_limit))
            call write_kv(unit, 'phonon_dos_limit_unit', 'cm-1')
        end if
        if (cfg%phonon_write_force_constants) &
            call write_kv(unit, 'phonon_write_force_constants', 'true')
        if (cfg%phonon_write_dynamical) &
            call write_kv(unit, 'phonon_write_dynamical', 'true')
        if (.not. cfg%phonon_calc_lo_to_splitting) &
            call write_kv(unit, 'phonon_calc_lo_to_splitting', 'false')
        if (cfg%phonon_force_constant_cutoff > 0.0_dp) then
            if (trim(cfg%phonon_fine_cutoff_method) == PHONON_CUTOFF_CUMULANT) then
                call write_kv(unit, 'phonon_force_constant_cutoff_scale', real2str(cfg%phonon_force_constant_cutoff))
            else
                call write_kv(unit, 'phonon_force_constant_cutoff', real2str(cfg%phonon_force_constant_cutoff))
            end if
        end if
        if (trim(cfg%phonon_fine_cutoff_method) /= PHONON_CUTOFF_CUMULANT) &
            call write_kv(unit, 'phonon_fine_cutoff_method', trim(cfg%phonon_fine_cutoff_method))
        if (cfg%phonon_max_cg_steps > 0) &
            call write_kv(unit, 'phonon_max_cg_steps', int2str(cfg%phonon_max_cg_steps))
        if (.not. cfg%phonon_use_kpoint_symmetry) &
            call write_kv(unit, 'phonon_use_kpoint_symmetry', 'false')
        if (cfg%calculate_raman) then
            call write_kv(unit, 'calculate_raman', 'true')
            call write_kv(unit, 'raman_method', trim(cfg%raman_method))
            call write_kv(unit, 'raman_range_low', '0.00e+00 cm-1')
            call write_kv(unit, 'raman_range_high', '1.00e+04 cm-1')
        end if
    end subroutine write_phonon_params

    subroutine write_efield_params(unit, cfg)
        integer, intent(in) :: unit
        type(castep_config_t), intent(in) :: cfg
        call write_kv(unit, 'efield_dfpt_method', trim(cfg%efield_dfpt_method))
        call write_kv(unit, 'efield_max_cycles', int2str(cfg%efield_max_cycles))
        call write_kv(unit, 'efield_energy_tol', energy_tol_str(cfg%efield_energy_tol))
        call write_kv(unit, 'efield_convergence_win', int2str(cfg%efield_convergence_win))
        if (cfg%efield_calc_ion_permittivity) &
            call write_kv(unit, 'efield_calc_ion_permittivity', 'true')
        call write_kv(unit, 'efield_freq_spacing', real2str(cfg%efield_freq_spacing))
        call write_kv(unit, 'efield_oscillator_q', real2str(cfg%efield_oscillator_q))
        if (trim(cfg%efield_calculate_nonlinear) /= 'NONE') &
            call write_kv(unit, 'efield_calculate_nonlinear', 'CHI2')
        call write_kv(unit, 'efield_ignore_molec_modes', trim(cfg%efield_ignore_molec_modes))
    end subroutine write_efield_params


    subroutine write_cineb_params(unit, cfg)
        !! Write CINEB (NEB + Climbing Image) transition state search parameters
        integer, intent(in) :: unit
        type(castep_config_t), intent(in) :: cfg

        call write_kv(unit, 'tssearch_method', 'NEB')
        call write_kv(unit, 'tssearch_max_path_points', trim(cfg%cineb_max_images))
        call write_kv(unit, 'tssearch_neb_spring_constant', &
            trim(cfg%cineb_spring_constant) // ' eV/ANG**2')
        call write_kv(unit, 'tssearch_neb_tangent_mode', trim(cfg%cineb_tangent_mode))
        call write_kv(unit, 'tssearch_neb_method', trim(cfg%cineb_neb_method))
        call write_kv(unit, 'tssearch_neb_max_iter', trim(cfg%cineb_max_iter))
        call write_kv(unit, 'tssearch_neb_climbing', trim(cfg%cineb_climbing))

        ! TS convergence tolerance (same 4-level system as geom tolerance)
        select case (trim(cfg%ts_geom_tolerance))
        case (GEO_COARSE)
            call write_kv(unit, 'tssearch_energy_tol', '5e-5')
            call write_kv(unit, 'tssearch_force_tol',  '0.1')
            call write_kv(unit, 'tssearch_disp_tol',   '0.005')
        case (GEO_MEDIUM)
            call write_kv(unit, 'tssearch_energy_tol', '2e-5')
            call write_kv(unit, 'tssearch_force_tol',  '0.05')
            call write_kv(unit, 'tssearch_disp_tol',   '0.002')
        case (GEO_FINE)
            call write_kv(unit, 'tssearch_energy_tol', '1e-5')
            call write_kv(unit, 'tssearch_force_tol',  '0.03')
            call write_kv(unit, 'tssearch_disp_tol',   '0.001')
        case (GEO_EXTREME)
            call write_kv(unit, 'tssearch_energy_tol', '5e-6')
            call write_kv(unit, 'tssearch_force_tol',  '0.01')
            call write_kv(unit, 'tssearch_disp_tol',   '5e-4')
        case default
            call write_kv(unit, 'tssearch_energy_tol', '2e-5')
            call write_kv(unit, 'tssearch_force_tol',  '0.05')
            call write_kv(unit, 'tssearch_disp_tol',   '0.002')
        end select
    end subroutine write_cineb_params

    pure function real2str(val) result(s)
        real(dp), intent(in) :: val
        character(20) :: s
        integer :: iv
        iv = int(val)
        if (abs(val - dble(iv)) <= 1.0e-10_dp) then
            write(s, '(I0)') iv
        else
            write(s, '(F14.8)') val
            s = adjustl(s)
        end if
    end function real2str

    pure function energy_tol_str(val) result(s)
        real(dp), intent(in) :: val
        character(10) :: s
        write(s, '(es9.1)') val
        s = adjustl(s)
    end function energy_tol_str

end module param_writer
