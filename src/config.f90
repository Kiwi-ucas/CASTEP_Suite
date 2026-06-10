module castep_config
    !! Central configuration types, constants, and initialization procedures
    !! for the CASTEP Suite CIF-to-CASTEP converter.
    implicit none

    ! Double precision kind
    integer, parameter, public :: dp = selected_real_kind(15, 307)

    ! Maximum sizes
    integer, parameter, public :: &
        MAX_ATOMS      = 10000, &
        MAX_LINE_LEN   = 1024, &
        MAX_TAGS       = 5000, &
        MAX_LOOP_ROWS  = 50000, &
        MAX_LOOP_COLS  = 50, &
        MAX_SYM_OPS    = 400

    ! PI constant
    real(dp), parameter, public :: pi = 3.14159265358979323846_dp
    real(dp), parameter, public :: HARTREE_TO_EV = 27.211386245988_dp

    ! Tag names used throughout
    character(len=*), parameter, public :: &
        TAG_A            = '_cell_length_a', &
        TAG_B            = '_cell_length_b', &
        TAG_C            = '_cell_length_c', &
        TAG_ALPHA        = '_cell_angle_alpha', &
        TAG_BETA         = '_cell_angle_beta', &
        TAG_GAMMA        = '_cell_angle_gamma', &
        TAG_SYM_HM       = '_symmetry_space_group_name_H-M', &
        TAG_SYM_IT       = '_symmetry_Int_Tables_number', &
        TAG_SYM_POS_XYZ  = '_symmetry_equiv_pos_as_xyz', &
        TAG_SPACEGROUP   = '_space_group_IT_number', &
        TAG_FORMULA_SUM  = '_chemical_formula_sum', &
        TAG_ATOM_LABEL   = '_atom_site_label', &
        TAG_ATOM_TYPE    = '_atom_site_type_symbol', &
        TAG_ATOM_FRAC_X  = '_atom_site_fract_x', &
        TAG_ATOM_FRAC_Y  = '_atom_site_fract_y', &
        TAG_ATOM_FRAC_Z  = '_atom_site_fract_z', &
        TAG_ATOM_OCC     = '_atom_site_occupancy', &
        TAG_ATOM_B       = '_atom_site_b_iso_or_equiv', &
        TAG_CELL_VOL     = '_cell_volume'

    ! Task type constants (16 types)
    character(len=*), parameter, public :: &
        TASK_ENERGY               = 'SINGLEPOINT', &
        TASK_ELECTRONIC_SPECTRO   = 'ElectronicSpectroscopy', &
        TASK_GEOMETRY_OPT         = 'GEOMETRYOPTIMISATION', &
        TASK_MOLECULAR_DYN        = 'MOLECULAR_DYNAMICS', &
        TASK_TRANSITION_STATE     = 'TRANSITIONSTATESEARCH', &
        TASK_PHONON               = 'PHONON', &
        TASK_EFIELD               = 'ELECTRIC_FIELD', &
        TASK_PHONON_EFIELD        = 'PHONON_ELECTRIC_FIELD', &
        TASK_THERMODYNAMICS       = 'THERMODYNAMICS', &
        TASK_MAGRES               = 'MAGNETIC_RESPONSE', &
        TASK_SPECTRAL             = 'SPECTRAL', &
        TASK_EPCOUPLING           = 'EP_COUPLING', &
        TASK_GENETIC_ALGO         = 'GENETIC_ALGORITHM', &
        TASK_SOCKET_DRIVER        = 'SOCKET_DRIVER', &
        TASK_ELASTIC              = 'ELASTIC', &
        TASK_AUTOSOLVATION        = 'AUTOSOLVATION'

    ! XC functional constants
    character(len=*), parameter, public :: &
        FUNC_PBE    = 'PBE', &
        FUNC_PBEsol   = 'PBEsol', &
        FUNC_HSE06  = 'HSE06', &
        FUNC_PBE0   = 'PBE0', &
        FUNC_r2scan = 'r2scan'

    ! vdW correction constants
    character(len=*), parameter, public :: &
        VDW_NONE   = 'NONE', &
        VDW_D3     = 'D3', &
        VDW_D3_BJ  = 'D3-BJ', &
        VDW_D4     = 'D4'

    ! Pseudopotential constants
    character(len=*), parameter, public :: &
        PSEUDO_NCP19 = 'NCP19', &
        PSEUDO_C19MK2 = 'C19MK2', &
        PSEUDO_SOC19 = 'SOC19'

    ! K-point scheme constants
    character(len=*), parameter, public :: &
        KPOINT_GAMMA              = 'GAMMA', &
        KPOINT_MONKHORST_PACK     = 'MONKHORST_PACK', &
        KPOINT_EXPLICIT           = 'EXPLICIT'

    ! Optimizer constants
    character(len=*), parameter, public :: &
        OPT_BFGS   = 'BFGS', &
        OPT_LBFGS  = 'LBFGS', &
        OPT_CG     = 'CG'

    ! Tolerance level constants
    character(len=*), parameter, public :: &
        TOL_SUPERFINE  = 'SUPERFINE', &
        TOL_FINE       = 'FINETOLERANT', &
        TOL_NORMAL     = 'NORMAL', &
        TOL_COARSE     = 'COARSE'

    ! Geo tolerance constants
    character(len=*), parameter, public :: &
        GEO_COARSE  = 'COARSE', &
        GEO_MEDIUM  = 'MEDIUM', &
        GEO_FINE    = 'FINE', &
        GEO_EXTREME = 'EXTREME'

    ! Symmetry source constants
    character(len=*), parameter, public :: &
        SYM_NONE  = 'NONE', &
        SYM_AUTO  = 'AUTO'

    ! Phonon method constants
    character(len=*), parameter, public :: &
        PHONON_METHOD_DFPT = 'DFPT', &
        PHONON_METHOD_FD   = 'FINITEDISPLACEMENT'
    character(len=*), parameter, public :: &
        PHONON_FINE_NONE        = 'NONE', &
        PHONON_FINE_SUPERCELL   = 'SUPERCELL', &
        PHONON_FINE_INTERPOLATE = 'INTERPOLATE'
    character(len=*), parameter, public :: &
        PHONON_DFPT_DM       = 'DM', &
        PHONON_DFPT_ALLBANDS = 'ALLBANDS'
    character(len=*), parameter, public :: &
        PHONON_SUM_NONE       = 'NONE', &
        PHONON_SUM_RECIPROCAL = 'RECIPROCAL', &
        PHONON_SUM_REALSPACE  = 'REALSPACE', &
        PHONON_SUM_REAL_RECIP = 'REAL-RECIP', &
        PHONON_SUM_MOLECULAR  = 'MOLECULAR'
    character(len=*), parameter, public :: &
        PHONON_QPOINT_MP_GRID = 'MP_GRID', &
        PHONON_QPOINT_PATH    = 'PATH'
    character(len=*), parameter, public :: &
        PHONON_CUTOFF_CUMULANT  = 'CUMULANT', &
        PHONON_CUTOFF_SPHERICAL = 'SPHERICAL'

    ! CINEB tangent mode constants
    character(len=*), parameter, public :: &
        CINEB_TANGENT_NONE   = 'NONE', &
        CINEB_TANGENT_BISECT = 'BISECT', &
        CINEB_TANGENT_HIGH_E = 'HIGH_E', &
        CINEB_TANGENT_SPLINE = 'SPLINE'

    ! CINEB NEB method (optimizer) constants
    character(len=*), parameter, public :: &
        CINEB_METHOD_TPSD   = 'TPSD', &
        CINEB_METHOD_FIRE   = 'FIRE', &
        CINEB_METHOD_ODE12R = 'ODE12R'

    ! Cell optimization constants
    character(len=*), parameter, public :: &
        CELL_ALL   = 'ALL', &
        CELL_INTE  = 'FIX_CELL'

    ! Iostat error codes
    integer, parameter, public :: &
        IO_SUCCESS       = 0, &
        IO_FILE_NOT_FOUND  = 100, &
        IO_PARSE_ERROR   = 101, &
        IO_MISSING_CELL  = 102, &
        IO_MISSING_ATOMS = 103, &
        IO_BAD_NUMERIC   = 104, &
    IO_WRITE_ERROR   = 105, &
        IO_INVALID_INPUT = 106, &
        IO_WRITE_FAIL    = 107, &
        IO_BANDS_NOT_FOUND   = 108, &
        IO_BANDS_PARSE_ERROR = 109, &
        IO_PDOS_NOT_FOUND    = 110, &
        IO_PDOS_PARSE_ERROR  = 111, &
        IO_USER_QUIT     = -1, &
        IO_PRECASTEP_LAUNCH = -2

    !! Atom data type
    type :: atom_t
        character(len=12) :: label
        character(len=6)  :: element
        real(dp) :: x, y, z
    end type atom_t

    !! CIF parsed data type
    type :: cif_data_t
        real(dp) :: a, b, c, alpha, beta, gamma
        character(len=32) :: space_group
        integer :: n_atoms
        type(atom_t), allocatable :: atoms(:)
        logical :: positions_fractional = .false.
    end type cif_data_t

    !! CASTEP configuration type
    type :: castep_config_t
        ! File paths
        character(len=1024) :: cif_file_path
        character(len=512) :: cell_output_path
        character(len=512) :: param_output_path

        ! Computation settings
        character(len=64)  :: task_type
        character(len=16)  :: xc_functional
        real(dp)           :: cutoff_energy
        character(len=16)  :: geom_tolerance
        character(len=16)  :: cell_opt_mode
        character(len=16)  :: sym_source
        character(len=16)  :: vdw_method
        character(len=8)   :: pseudopotential
        character(len=16)  :: kpoint_scheme
        integer            :: kpoint_grid(3)
        real(dp)           :: scf_tolerance = 1.0e-5_dp
        character(len=16)  :: optimizer

        ! Cell parameters
        real(dp) :: cell_length(3)    ! a, b, c (Angstrom)
        real(dp) :: cell_angle(3)     ! alpha, beta, gamma (degrees)
        real(dp) :: cell_basis(3, 3)  ! Cartesian lattice vectors

        ! Spin
        logical :: spin_polarized

        ! Coordinate system (CIF=fractional, PDB/.cell=Cartesian)
        logical :: cartesian_coords

        ! Atom data
        integer :: num_atoms
        character(len=8), allocatable :: atom_type(:)
        real(dp), allocatable         :: atom_x(:)
        real(dp), allocatable         :: atom_y(:)
        real(dp), allocatable         :: atom_z(:)

        ! CINEB product atom data
        integer :: prod_num_atoms = 0
        character(len=8), allocatable :: prod_atom_type(:)
        real(dp), allocatable         :: prod_atom_x(:)
        real(dp), allocatable         :: prod_atom_y(:)
        real(dp), allocatable         :: prod_atom_z(:)
        logical :: prod_cartesian_coords = .false.

        ! CINEB intermediate atom data
        integer :: interm_num_atoms = 0
        character(len=8), allocatable :: interm_atom_type(:)
        real(dp), allocatable         :: interm_atom_x(:)
        real(dp), allocatable         :: interm_atom_y(:)
        real(dp), allocatable         :: interm_atom_z(:)
        logical :: interm_cartesian_coords = .false.

        ! CINEB file paths
        character(len=1024) :: prod_file_path    = ''
        character(len=1024) :: interm_file_path  = ''

        ! CINEB parameters
        character(len=16) :: cineb_max_images      = '11'
        character(len=16) :: cineb_spring_constant = '0.1'
        character(len=16) :: cineb_max_iter        = '50'
        character(len=16) :: cineb_tangent_mode    = 'SPLINE'
        character(len=16) :: cineb_neb_method      = 'ODE12R'
        character(len=16) :: cineb_climbing        = 'TRUE'
        character(len=16) :: ts_geom_tolerance     = 'MEDIUM'

        ! Symmetry
        logical :: has_space_group
        character(len=128) :: space_group_name

        ! Metadata
        character(len=256) :: formula_sum

        ! Advanced options
        logical :: smearing
        integer  :: max_scf_cycles
        integer  :: elec_convergence_win
        logical :: calculate_elf
        logical :: calculate_edd

        ! Spectral task sub-options
        character(len=32) :: spectral_task_type

        ! Phonon settings
        character(len=20) :: phonon_method          = 'DFPT'
        character(len=20) :: phonon_fine_method     = 'INTERPOLATE'
        character(len=16) :: phonon_dfpt_method     = 'DM'
        character(len=16) :: phonon_sum_rule_method = 'RECIPROCAL'
        real(dp)          :: phonon_energy_tol      = 1.0e-5_dp
        integer           :: phonon_max_cycles      = 50
        integer           :: phonon_convergence_win = 3
        logical           :: phonon_calculate_dos   = .false.
        real(dp)          :: phonon_dos_spacing     = 10.0_dp
        real(dp)          :: phonon_finite_disp     = 0.01_dp
        integer           :: phonon_kpoint_mp_grid(3) = [1, 1, 1]
        character(len=16) :: phonon_qpoint_scheme        = 'MP_GRID'
        character(len=1024) :: phonon_kpoint_path         = ''
        real(dp)          :: phonon_kpoint_path_spacing  = 0.1_dp
        integer           :: phonon_supercell_matrix(3,3) = &
            reshape([1,0,0, 0,1,0, 0,0,1], [3,3])
        character(len=16) :: phonon_fine_qpoint_scheme        = 'MP_GRID'
        integer           :: phonon_fine_kpoint_mp_grid(3)    = [1, 1, 1]
        character(len=1024) :: phonon_fine_kpoint_path         = ''
        real(dp)          :: phonon_fine_kpoint_path_spacing  = 0.1_dp
        logical           :: phonon_write_force_constants = .false.
        logical           :: phonon_write_dynamical       = .false.
        logical           :: phonon_calc_lo_to_splitting  = .true.
        real(dp)          :: phonon_force_constant_cutoff = 0.0_dp
        character(len=16) :: phonon_fine_cutoff_method   = 'CUMULANT'
        real(dp)          :: phonon_dos_limit             = 5000.0_dp
        integer           :: phonon_max_cg_steps          = 0
        logical           :: phonon_use_kpoint_symmetry   = .true.
        logical           :: calculate_born_charges       = .true.
        logical           :: calculate_raman              = .false.
        character(len=20) :: raman_method                = 'DFPT'
        character(len=16) :: efield_dfpt_method          = 'ALLBANDS'
        integer           :: efield_max_cycles           = 50
        integer           :: efield_convergence_win      = 2
        real(dp)          :: efield_energy_tol           = 1.0e-5_dp
        logical           :: efield_calc_ion_permittivity = .true.
        character(len=16) :: efield_ignore_molec_modes   = 'CRYSTAL'
        real(dp)          :: efield_freq_spacing         = 1.0_dp
        real(dp)          :: efield_oscillator_q         = 50.0_dp
        character(len=16) :: efield_calculate_nonlinear  = 'NONE'
    contains
        final :: finalize_castep_config
    end type castep_config_t

    type :: bands_data_t
        integer  :: num_kpoints       = 0
        integer  :: num_spin          = 1
        real(dp) :: num_electrons     = 0.0_dp
        integer  :: num_eigenvalues   = 0
        real(dp) :: fermi_energy      = 0.0_dp
        real(dp) :: cell_vectors(3,3) = 0.0_dp
        integer,  allocatable :: kpoint_indices(:)
        real(dp), allocatable :: kpoint_coords(:,:)
        real(dp), allocatable :: kpath_dist(:)
        real(dp), allocatable :: eigenvalues(:,:,:)
    end type bands_data_t

    !! PDOS data type — parsed .pdos_weights / .pdos_bin content
    type :: pdos_data_t
        integer  :: total_kpoints  = 0
        integer  :: num_spins      = 1
        integer  :: num_orbitals   = 0
        integer  :: max_bands      = 0
        integer,  allocatable :: orbital_species(:)
        integer,  allocatable :: orbital_ion(:)
        integer,  allocatable :: orbital_am(:)         ! 0=S, 1=P, 2=D, 3=F
        real(dp), allocatable :: kpoint_coords(:,:)    ! (3, nk)
        integer,  allocatable :: kpoint_indices(:)      ! (nk)
        real(dp), allocatable :: orbital_weights(:,:,:,:) ! (norbs, nbands, nk, nspin)
    end type pdos_data_t

    ! Public procedures
    public :: default_config
    public :: strip_quotes
    public :: new_castep_config
    public :: normalize_tag
    public :: cif_data_t
    public :: atom_t
    public :: bands_data_t
    public :: pdos_data_t
    public :: compare_tags
    public :: string_to_real
    public :: get_castep_task_name, int2str

contains

    subroutine default_config(cfg)
        !! Set all castep_config_t fields to CASTEP sensible defaults
        type(castep_config_t), intent(out) :: cfg

        cfg%cif_file_path    = ''
        cfg%cell_output_path  = ''
        cfg%param_output_path = ''
        cfg%task_type         = TASK_ENERGY
        cfg%xc_functional     = FUNC_PBE
        cfg%cutoff_energy     = 400.0_dp
        cfg%geom_tolerance    = GEO_MEDIUM
        cfg%cell_opt_mode     = CELL_INTE
        cfg%sym_source        = SYM_NONE
        cfg%vdw_method        = VDW_NONE
        cfg%pseudopotential   = PSEUDO_C19MK2
        cfg%kpoint_scheme     = KPOINT_GAMMA
        cfg%kpoint_grid       = 0
        cfg%scf_tolerance     = 1.0e-5_dp
        cfg%optimizer         = OPT_BFGS
        cfg%cell_length      = 0.0_dp
        cfg%cell_angle       = 90.0_dp
        cfg%cell_basis       = 0.0_dp
        cfg%num_atoms        = 0
        cfg%has_space_group  = .false.
        cfg%space_group_name = ''
        cfg%formula_sum      = ''
        cfg%spin_polarized     = .false.
        cfg%cartesian_coords   = .false.
        cfg%smearing             = .false.
        cfg%max_scf_cycles       = 256
        cfg%elec_convergence_win = 3
        cfg%calculate_elf        = .false.
        cfg%calculate_edd        = .false.
        cfg%spectral_task_type   = 'BandStructure'
        cfg%phonon_method        = PHONON_METHOD_DFPT
        cfg%phonon_fine_method   = PHONON_FINE_INTERPOLATE
        cfg%phonon_dfpt_method   = PHONON_DFPT_DM
        cfg%phonon_sum_rule_method = PHONON_SUM_RECIPROCAL
        cfg%phonon_energy_tol    = 1.0e-5_dp
        cfg%phonon_max_cycles    = 50
        cfg%phonon_convergence_win = 3
        cfg%phonon_calculate_dos = .false.
        cfg%phonon_dos_spacing   = 10.0_dp
        cfg%phonon_finite_disp   = 0.01_dp
        cfg%phonon_kpoint_mp_grid  = [1, 1, 1]
        cfg%phonon_qpoint_scheme   = PHONON_QPOINT_MP_GRID
        cfg%phonon_kpoint_path     = ''
        cfg%phonon_kpoint_path_spacing = 0.1_dp
        cfg%phonon_supercell_matrix = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])
        cfg%phonon_fine_qpoint_scheme  = PHONON_QPOINT_MP_GRID
        cfg%phonon_fine_kpoint_mp_grid = [1, 1, 1]
        cfg%phonon_fine_kpoint_path    = ''
        cfg%phonon_fine_kpoint_path_spacing = 0.1_dp
        cfg%phonon_write_force_constants = .false.
        cfg%phonon_write_dynamical       = .false.
        cfg%phonon_calc_lo_to_splitting  = .true.
        cfg%phonon_force_constant_cutoff = 0.0_dp
        cfg%phonon_fine_cutoff_method    = PHONON_CUTOFF_CUMULANT
        cfg%phonon_dos_limit             = 5000.0_dp
        cfg%phonon_max_cg_steps         = 0
        cfg%phonon_use_kpoint_symmetry  = .true.
        cfg%calculate_born_charges      = .true.
        cfg%calculate_raman             = .false.
        cfg%raman_method                = 'DFPT'
        cfg%efield_dfpt_method          = 'ALLBANDS'
        cfg%efield_max_cycles           = 50
        cfg%efield_convergence_win      = 2
        cfg%efield_energy_tol           = 1.0e-5_dp
        cfg%efield_calc_ion_permittivity = .true.
        cfg%efield_ignore_molec_modes   = 'CRYSTAL'
        cfg%efield_freq_spacing         = 1.0_dp
        cfg%efield_oscillator_q         = 50.0_dp
        cfg%efield_calculate_nonlinear  = 'NONE'

        ! CINEB defaults
        cfg%prod_num_atoms        = 0
        cfg%interm_num_atoms      = 0
        cfg%prod_cartesian_coords  = .false.
        cfg%interm_cartesian_coords = .false.
        cfg%prod_file_path         = ''
        cfg%interm_file_path       = ''
        cfg%cineb_max_images       = '11'
        cfg%cineb_spring_constant  = '0.1'
        cfg%cineb_max_iter         = '50'
        cfg%cineb_tangent_mode     = CINEB_TANGENT_SPLINE
        cfg%cineb_neb_method       = CINEB_METHOD_ODE12R
        cfg%cineb_climbing         = 'TRUE'
        cfg%ts_geom_tolerance      = GEO_MEDIUM

        if (allocated(cfg%atom_type)) deallocate(cfg%atom_type)
        if (allocated(cfg%atom_x)) deallocate(cfg%atom_x)
        if (allocated(cfg%atom_y)) deallocate(cfg%atom_y)
        if (allocated(cfg%atom_z)) deallocate(cfg%atom_z)
        if (allocated(cfg%prod_atom_type)) deallocate(cfg%prod_atom_type)
        if (allocated(cfg%prod_atom_x)) deallocate(cfg%prod_atom_x)
        if (allocated(cfg%prod_atom_y)) deallocate(cfg%prod_atom_y)
        if (allocated(cfg%prod_atom_z)) deallocate(cfg%prod_atom_z)
        if (allocated(cfg%interm_atom_type)) deallocate(cfg%interm_atom_type)
        if (allocated(cfg%interm_atom_x)) deallocate(cfg%interm_atom_x)
        if (allocated(cfg%interm_atom_y)) deallocate(cfg%interm_atom_y)
        if (allocated(cfg%interm_atom_z)) deallocate(cfg%interm_atom_z)
    end subroutine default_config

    type(castep_config_t) function new_castep_config() result(cfg)
        !! Factory constructor with defaults
        call default_config(cfg)
end function new_castep_config

    pure logical function compare_tags(a, b)
        !! Case-insensitive tag comparison
        character(len=*), intent(in) :: a
        character(len=*), intent(in) :: b
        integer :: i, la, lb
        character(1) :: ca, cb

        compare_tags = .false.
        la = len_trim(a)
        lb = len_trim(b)
        if (la /= lb) return

        do i = 1, la
            ca = a(i:i)
            cb = b(i:i)
            ! Convert to lowercase for comparison
            if (ca >= 'A' .and. ca <= 'Z') ca = char(iachar(ca) + 32)
            if (cb >= 'A' .and. cb <= 'Z') cb = char(iachar(cb) + 32)
            if (ca /= cb) return
        end do
        compare_tags = .true.
    end function compare_tags

    subroutine normalize_tag(tag, normalized)
        !! Normalize a CIF tag: strip leading underscore/dots,
        !! convert CamelCase to snake_case
        character(len=*), intent(in) :: tag
        character(len=*), intent(out) :: normalized
        integer :: i, len_tag, idx
        character(1) :: ch
        logical :: prev_upper, needs_underscore

        len_tag = len_trim(tag)
        idx = 1
        normalized = ''

        ! Strip leading underscores and dots
        i = 1
        do while (i <= len_tag)
            ch = tag(i:i)
            if (ch /= '_' .and. ch /= '.') exit
            i = i + 1
        end do

        ! Convert CamelCase and collect
        prev_upper = .false.
        needs_underscore = .false.
        do while (i <= len_tag .and. idx <= len(normalized))
            ch = tag(i:i)

            if (ch >= 'A' .and. ch <= 'Z') then
                if (prev_upper .and. .not. needs_underscore) then
                    ! Check if next char is lower case (CamelCase boundary)
                    if (i + 1 <= len_tag) then
                        if (tag(i+1:i+1) >= 'a' .and. tag(i+1:i+1) <= 'z') then
                            if (idx > 1 .and. normalized(idx-1:idx-1) /= '_') then
                                       normalized(idx:idx) = '_'
                                       idx = idx + 1
                            end if
                        end if
                    end if
                end if
                normalized(idx:idx) = char(iachar(ch) + 32)
                prev_upper = .true.
            else if (ch >= 'a' .and. ch <= 'z') then
                normalized(idx:idx) = ch
                prev_upper = .false.
            else if (ch >= '0' .and. ch <= '9') then
                normalized(idx:idx) = ch
                prev_upper = .false.
            else if (ch == '_') then
                normalized(idx:idx) = ch
                prev_upper = .false.
            else
                ! Other non-alphanumeric: skip
                i = i + 1
                cycle
            end if
            idx = idx + 1
            i = i + 1
        end do

        normalized = adjustl(normalized)
    end subroutine normalize_tag

    function string_to_real(s, val) result(success)
        !! Safely convert a string to real(dp)
        character(len=*), intent(in) :: s
        real(dp), intent(out) :: val
        logical :: success
        character(len=256) :: trimmed
        integer :: ios

        trimmed = adjustl(trim(s))
        if (len_trim(trimmed) == 0) then
            success = .false.
            return
        end if

        read(trimmed, '(F20.10)', iostat=ios) val
        success = (ios == 0)
    end function string_to_real

    pure function get_castep_task_name(task) result(castep_name)
        !! Map internal task constant to CASTEP internal task name
        character(len=*), intent(in) :: task
        character(len=64) :: castep_name

        select case (trim(task))
        case (TASK_ENERGY)
            castep_name = 'SINGLEPOINT'
        case (TASK_GEOMETRY_OPT)
            castep_name = 'GEOMETRYOPTIMISATION'
        case (TASK_MOLECULAR_DYN)
            castep_name = 'MOLECULARDYNAMICS'
        case (TASK_PHONON)
            castep_name = 'PHONON'
        case (TASK_ELECTRONIC_SPECTRO)
            castep_name = 'ElectronicSpectroscopy'
        case (TASK_TRANSITION_STATE)
            castep_name = TASK_TRANSITION_STATE
        case (TASK_EFIELD)
            castep_name = 'EFIELD'
        case (TASK_PHONON_EFIELD)
            castep_name = 'Phonon+Efield'
        case (TASK_THERMODYNAMICS)
            castep_name = 'THERMODYNAMICS'
        case (TASK_MAGRES)
            castep_name = 'MAGNETIC_RESPONSE'
        case (TASK_SPECTRAL)
            castep_name = 'SPECTRAL'
        case (TASK_EPCOUPLING)
            castep_name = 'EP_COUPLING'
        case (TASK_GENETIC_ALGO)
            castep_name = 'GENETIC_ALGORITHM'
        case (TASK_SOCKET_DRIVER)
            castep_name = 'SOCKET_DRIVER'
        case (TASK_ELASTIC)
            castep_name = 'ELASTIC'
        case (TASK_AUTOSOLVATION)
            castep_name = 'AUTOSOLVATION'
        case default
            castep_name = trim(task)
        end select
    end function get_castep_task_name

    pure function int2str(val) result(s)
        !! Convert integer to trimmed string
        integer, intent(in) :: val
        character(16) :: s
        write(s, '(I6)') val
        s = adjustl(s)
    end function int2str

    pure subroutine strip_quotes(s)
        !! Remove leading/trailing ' or " characters in place
        character(len=*), intent(inout) :: s
        integer :: n
        n = len_trim(s)
        if (n < 1) return
        if (s(1:1) == "'" .or. s(1:1) == '"') s = s(2:)
        n = len_trim(s)
        if (n < 1) return
        if (s(n:n) == "'" .or. s(n:n) == '"') s = s(:n-1)
    end subroutine strip_quotes

    subroutine finalize_castep_config(this)
        !! Finalizer: deallocate allocatable components of castep_config_t
        type(castep_config_t), intent(inout) :: this
        if (allocated(this%atom_type)) deallocate(this%atom_type)
        if (allocated(this%atom_x))    deallocate(this%atom_x)
        if (allocated(this%atom_y))    deallocate(this%atom_y)
        if (allocated(this%atom_z))    deallocate(this%atom_z)
        if (allocated(this%prod_atom_type)) deallocate(this%prod_atom_type)
        if (allocated(this%prod_atom_x))    deallocate(this%prod_atom_x)
        if (allocated(this%prod_atom_y))    deallocate(this%prod_atom_y)
        if (allocated(this%prod_atom_z))    deallocate(this%prod_atom_z)
        if (allocated(this%interm_atom_type)) deallocate(this%interm_atom_type)
        if (allocated(this%interm_atom_x))    deallocate(this%interm_atom_x)
        if (allocated(this%interm_atom_y))    deallocate(this%interm_atom_y)
        if (allocated(this%interm_atom_z))    deallocate(this%interm_atom_z)
    end subroutine finalize_castep_config

end module castep_config
