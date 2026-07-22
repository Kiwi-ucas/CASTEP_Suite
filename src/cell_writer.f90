module cell_writer
    !! CASTEP .cell file generator using %BLOCK format
    !! Each %BLOCK type is written by an independent subroutine.
    !! Extensible: add new write_block_xxx subroutines and call them
    !! from write_cell_file.
    use castep_config, only: &
         dp, CELL_ALL, CELL_INTE, TASK_GEOMETRY_OPT, TASK_ELECTRONIC_SPECTRO, SYM_AUTO, &
         TASK_PHONON, TASK_PHONON_EFIELD, TASK_THERMODYNAMICS, TASK_TRANSITION_STATE, &
         PHONON_FINE_NONE, PHONON_QPOINT_MP_GRID, PHONON_QPOINT_PATH, &
         PHONON_METHOD_FD, PHONON_FINE_SUPERCELL, &
         KPOINT_GAMMA, KPOINT_MONKHORST_PACK, &
         compare_tags, IO_WRITE_FAIL, &
         compute_cartesian_lattice, &
         castep_config_t
    implicit none
    private

    public :: write_cell_file
    public :: write_block_lattice_abc, write_block_lattice_cart
    public :: write_block_positions_abs, write_block_positions_frac
    public :: write_block_species_pot, write_block_kpoint_grid
    public :: write_block_ionic_constraints

contains

    subroutine write_block_lattice_abc(unit, cfg)
        !! Write %BLOCK LATTICE_ABC with lattice parameters (degrees)
        integer, intent(in) :: unit
        type(castep_config_t), intent(in) :: cfg

        write(unit, '(a)') '#%BLOCK LATTICE_ABC  (commented out — for reference)'
        write(unit, '(a,3(f12.7,1x))') '#', cfg%cell_length(1), cfg%cell_length(2), cfg%cell_length(3)
        write(unit, '(a,3(f12.7,1x))') '#', cfg%cell_angle(1), cfg%cell_angle(2), cfg%cell_angle(3)
        write(unit, '(a)') '#%ENDBLOCK LATTICE_ABC'
        write(unit, '(a)') ''
    end subroutine write_block_lattice_abc

    subroutine write_block_lattice_cart(unit, cfg)
        !! Write %BLOCK LATTICE_CART from cfg%cell_basis.
        integer, intent(in) :: unit
        type(castep_config_t), intent(in) :: cfg
        integer :: i
        write(unit, '(a)') '%BLOCK LATTICE_CART'
        do i = 1, 3
            write(unit, '(3(f12.7, 1x))') cfg%cell_basis(i,1), cfg%cell_basis(i,2), cfg%cell_basis(i,3)
        end do
        write(unit, '(a)') '%ENDBLOCK LATTICE_CART'
        write(unit, '(a)') ''
    end subroutine write_block_lattice_cart

    subroutine write_block_species_pot(unit, cfg)
        !! Write %BLOCK SPECIES_POT with unique atom kinds
        integer, intent(in) :: unit
        type(castep_config_t), intent(in) :: cfg
        integer :: i, j, n_kinds
        logical, allocatable :: is_kind(:)
        character(len=8), allocatable :: kind_labels(:)

        if (cfg%num_atoms == 0) return

        n_kinds = 0
        allocate(kind_labels(cfg%num_atoms))
        allocate(is_kind(cfg%num_atoms))
        is_kind = .false.

        do i = 1, cfg%num_atoms
            if (.not. is_kind(i)) then
                n_kinds = n_kinds + 1
                kind_labels(n_kinds) = trim(cfg%atom_type(i))
                is_kind(i) = .true.
                do j = i + 1, cfg%num_atoms
                    if (cfg%atom_type(j) == cfg%atom_type(i)) then
                        is_kind(j) = .true.
                    end if
                end do
            end if
        end do

        write(unit, '(a)') '%BLOCK SPECIES_POT'
        do i = 1, n_kinds
            write(unit, '(a, 1x, a)') trim(kind_labels(i)), trim(cfg%pseudopotential)
        end do
        write(unit, '(a)') '%ENDBLOCK SPECIES_POT'
        write(unit, '(a)') ''

        deallocate(kind_labels, is_kind)
    end subroutine write_block_species_pot

    subroutine write_block_positions_abs(unit, cfg)
        !! Write %BLOCK POSITIONS_ABS with Cartesian coordinates.
        !! Converts from fractional via cfg%cell_basis if cartesian_coords is false.
        integer, intent(in) :: unit
        type(castep_config_t), intent(in) :: cfg
        real(dp), allocatable :: cart_coords(:,:)
        integer :: i, j

        if (cfg%num_atoms == 0) return

        allocate(cart_coords(cfg%num_atoms, 3))
        if (cfg%cartesian_coords) then
            do i = 1, cfg%num_atoms
                cart_coords(i,1) = cfg%atom_x(i)
                cart_coords(i,2) = cfg%atom_y(i)
                cart_coords(i,3) = cfg%atom_z(i)
            end do
        else
            do i = 1, cfg%num_atoms
                do j = 1, 3
                    cart_coords(i,j) = cfg%cell_basis(j,1) * cfg%atom_x(i) &
                                     + cfg%cell_basis(j,2) * cfg%atom_y(i) &
                                     + cfg%cell_basis(j,3) * cfg%atom_z(i)
                end do
            end do
        end if

        write(unit, '(a)') '#%BLOCK POSITIONS_ABS  (commented out — for reference)'
        do i = 1, cfg%num_atoms
            write(unit, '(a, a, 3(1x, f10.7))') '#', trim(cfg%atom_type(i)), &
                cart_coords(i,1), cart_coords(i,2), cart_coords(i,3)
        end do
        write(unit, '(a)') '#%ENDBLOCK POSITIONS_ABS'
        write(unit, '(a)') ''

        deallocate(cart_coords)
    end subroutine write_block_positions_abs


    subroutine write_block_positions_abs_product(unit, cfg)
        !! Write %BLOCK POSITIONS_ABS_PRODUCT for CINEB product structure
        integer, intent(in) :: unit
        type(castep_config_t), intent(in) :: cfg
        real(dp), allocatable :: cart_coords(:,:)
        integer :: i, j

        if (cfg%prod_num_atoms == 0) return

        allocate(cart_coords(cfg%prod_num_atoms, 3))
        if (cfg%prod_cartesian_coords) then
            do i = 1, cfg%prod_num_atoms
                cart_coords(i,1) = cfg%prod_atom_x(i)
                cart_coords(i,2) = cfg%prod_atom_y(i)
                cart_coords(i,3) = cfg%prod_atom_z(i)
            end do
        else
            do i = 1, cfg%prod_num_atoms
                do j = 1, 3
                    cart_coords(i,j) = cfg%cell_basis(j,1) * cfg%prod_atom_x(i) &
                                     + cfg%cell_basis(j,2) * cfg%prod_atom_y(i) &
                                     + cfg%cell_basis(j,3) * cfg%prod_atom_z(i)
                end do
            end do
        end if

        write(unit, '(a)') '%BLOCK POSITIONS_ABS_PRODUCT'
        do i = 1, cfg%prod_num_atoms
            write(unit, '(a, 3(1x, f10.7))') trim(cfg%prod_atom_type(i)), &
                cart_coords(i,1), cart_coords(i,2), cart_coords(i,3)
        end do
        write(unit, '(a)') '%ENDBLOCK POSITIONS_ABS_PRODUCT'
        write(unit, '(a)') ''

        deallocate(cart_coords)
    end subroutine write_block_positions_abs_product


    subroutine write_block_positions_abs_intermediate(unit, cfg)
        !! Write %BLOCK POSITIONS_ABS_INTERMEDIATE for CINEB intermediate structure
        integer, intent(in) :: unit
        type(castep_config_t), intent(in) :: cfg
        real(dp), allocatable :: cart_coords(:,:)
        integer :: i, j

        if (cfg%interm_num_atoms == 0) return

        allocate(cart_coords(cfg%interm_num_atoms, 3))
        if (cfg%interm_cartesian_coords) then
            do i = 1, cfg%interm_num_atoms
                cart_coords(i,1) = cfg%interm_atom_x(i)
                cart_coords(i,2) = cfg%interm_atom_y(i)
                cart_coords(i,3) = cfg%interm_atom_z(i)
            end do
        else
            do i = 1, cfg%interm_num_atoms
                do j = 1, 3
                    cart_coords(i,j) = cfg%cell_basis(j,1) * cfg%interm_atom_x(i) &
                                     + cfg%cell_basis(j,2) * cfg%interm_atom_y(i) &
                                     + cfg%cell_basis(j,3) * cfg%interm_atom_z(i)
                end do
            end do
        end if

        write(unit, '(a)') '%BLOCK POSITIONS_ABS_INTERMEDIATE'
        do i = 1, cfg%interm_num_atoms
            write(unit, '(a, 3(1x, f10.7))') trim(cfg%interm_atom_type(i)), &
                cart_coords(i,1), cart_coords(i,2), cart_coords(i,3)
        end do
        write(unit, '(a)') '%ENDBLOCK POSITIONS_ABS_INTERMEDIATE'
        write(unit, '(a)') ''

        deallocate(cart_coords)
    end subroutine write_block_positions_abs_intermediate

    subroutine write_block_cell_constraints(unit, cfg)
        !! Write %BLOCK CELL_CONSTRAINTS for optimization tasks
        integer, intent(in) :: unit
        type(castep_config_t), intent(in) :: cfg

        if (trim(cfg%cell_opt_mode) == CELL_ALL) then
            write(unit, '(a)') '%BLOCK CELL_CONSTRAINTS'
            write(unit, '(a)') '  1        2       3'
            write(unit, '(a)') '  4        5       6'
            write(unit, '(a)') '%ENDBLOCK CELL_CONSTRAINTS'
            write(unit, '(a)') ''
        else if (trim(cfg%cell_opt_mode) == CELL_INTE) then
            write(unit, '(a)') 'FIX_ALL_CELL : true'
            write(unit, '(a)') ''
        end if
    end subroutine write_block_cell_constraints

    subroutine write_block_kpoint_grid(unit, cfg)
        !! Write KPOINTS_MP_GRID based on k-point scheme
        integer, intent(in) :: unit
        type(castep_config_t), intent(in) :: cfg

        if (trim(cfg%kpoint_scheme) == KPOINT_GAMMA) then
            write(unit, '(a, a, a)') 'KPOINTS_MP_GRID : ', '1  1  1'
        else if (trim(cfg%kpoint_scheme) == KPOINT_MONKHORST_PACK) then
            write(unit, '(a, i4, 1x, i4, 1x, i4)') &
                'KPOINTS_MP_GRID : ', cfg%kpoint_grid(1), &
                cfg%kpoint_grid(2), cfg%kpoint_grid(3)
        end if
    end subroutine write_block_kpoint_grid

    pure logical function is_phonon_task(cfg)
        type(castep_config_t), intent(in) :: cfg
        is_phonon_task = (trim(cfg%task_type) == TASK_PHONON &
                     .or. trim(cfg%task_type) == TASK_PHONON_EFIELD &
                     .or. trim(cfg%task_type) == TASK_THERMODYNAMICS)
    end function is_phonon_task

    subroutine write_phonon_kpoint_mp(unit, keyword, grid)
        integer, intent(in) :: unit
        character(len=*), intent(in) :: keyword
        integer, intent(in) :: grid(3)
        write(unit, '(a, i4, 1x, i4, 1x, i4)') &
            trim(keyword) // ' : ', grid(1), grid(2), grid(3)
    end subroutine write_phonon_kpoint_mp

    subroutine write_phonon_path_block(unit, block_name, path_str, spacing_kw, spacing_val)
        integer, intent(in) :: unit
        character(len=*), intent(in) :: block_name, path_str, spacing_kw
        real(dp), intent(in) :: spacing_val
        character(len=1024) :: buf
        integer :: i, n, start
        if (len_trim(path_str) == 0) return
        write(unit, '(a)') '%BLOCK ' // trim(block_name)
        buf = path_str
        n = len_trim(buf)
        start = 1
        do i = 1, n
            if (buf(i:i) == achar(10)) then
                if (i > start) write(unit, '(a)') trim(buf(start:i-1))
                start = i + 1
            end if
        end do
        if (start <= n) write(unit, '(a)') trim(buf(start:n))
        write(unit, '(a)') '%ENDBLOCK ' // trim(block_name)
        write(unit, '(a, f6.3)') trim(spacing_kw) // ' : ', spacing_val
        write(unit, '(a)') ''
    end subroutine write_phonon_path_block

    subroutine write_block_phonon_supercell_matrix(unit, cfg)
        integer, intent(in) :: unit
        type(castep_config_t), intent(in) :: cfg
        integer :: i
        write(unit, '(a)') '%BLOCK PHONON_SUPERCELL_MATRIX'
        do i = 1, 3
            write(unit, '(3(i4, 1x))') cfg%phonon_supercell_matrix(i, 1), &
                cfg%phonon_supercell_matrix(i, 2), cfg%phonon_supercell_matrix(i, 3)
        end do
        write(unit, '(a)') '%ENDBLOCK PHONON_SUPERCELL_MATRIX'
        write(unit, '(a)') ''
    end subroutine write_block_phonon_supercell_matrix

    subroutine write_cell_file(filename, cfg, iostat, iomsg)
        !! Write a CASTEP .cell file in %BLOCK format.
        !! Sets cfg%cell_basis if not already computed.
        character(len=*), intent(in) :: filename
        type(castep_config_t), intent(inout) :: cfg
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg
        integer :: unit, ios

        iostat = 0

        open(newunit=unit, file=trim(filename), status='unknown', action='write', iostat=ios)
        if (ios /= 0) then
            iostat = IO_WRITE_FAIL
            if (present(iomsg)) iomsg = 'Cannot write file: ' // trim(filename)
            return
        end if

        write(unit, '(a)') '! CASTEP Suite generated cell file'
        write(unit, '(a)') ''

        ! Ensure lattice basis is computed
        if (all(abs(cfg%cell_basis) < 1.0e-12_dp)) then
            cfg%cell_basis = compute_cartesian_lattice(cfg%cell_length(1), cfg%cell_length(2), &
                               cfg%cell_length(3), cfg%cell_angle(1), &
                               cfg%cell_angle(2), cfg%cell_angle(3))
        end if

        call write_block_lattice_abc(unit, cfg)
        call write_block_lattice_cart(unit, cfg)

        call write_block_positions_abs(unit, cfg)
        call write_block_positions_frac(unit, cfg)
        if (cfg%pes_mobile_idx > 0) then
            call write_block_ionic_constraints(unit, cfg, cfg%pes_mobile_idx, cfg%pes_fix_axes)
        end if
        if (trim(cfg%task_type) == TASK_GEOMETRY_OPT) then
            call write_block_cell_constraints(unit, cfg)
            write(unit, '(a)') 'FIX_COM : false'
            write(unit, '(a)') ''
        end if
        if (trim(cfg%task_type) == TASK_TRANSITION_STATE) then
            call write_block_positions_abs_product(unit, cfg)
            call write_block_positions_abs_intermediate(unit, cfg)
        end if

        call write_block_species_pot(unit, cfg)
        if (trim(cfg%sym_source) == SYM_AUTO) then
            write(unit, '(a)') 'SYMMETRY_GENERATE'
            write(unit, '(a)') ''
        end if
        
        ! K-point grid settings
        call write_block_kpoint_grid(unit, cfg)

        ! Spectral task k-point settings
        if (trim(cfg%task_type) == TASK_ELECTRONIC_SPECTRO) then
            write(unit, '(a, f6.2)') 'SPECTRAL_KPOINT_MP_SPACING ', 0.05_dp
            write(unit, '(a)') '%BLOCK SPECTRAL_KPOINT_PATH'
            write(unit, '(a)') '%ENDBLOCK SPECTRAL_KPOINT_PATH'
            write(unit, '(a)') ''
        end if

        ! Phonon task settings
        if (is_phonon_task(cfg)) then
            if (trim(cfg%phonon_qpoint_scheme) == PHONON_QPOINT_PATH) then
                call write_phonon_path_block(unit, 'PHONON_KPOINT_PATH', &
                    cfg%phonon_kpoint_path, 'PHONON_KPOINT_PATH_SPACING', &
                    cfg%phonon_kpoint_path_spacing)
            else
                call write_phonon_kpoint_mp(unit, 'PHONON_KPOINT_MP_GRID', cfg%phonon_kpoint_mp_grid)
            end if
            if (trim(cfg%phonon_method) == PHONON_METHOD_FD &
                .or. trim(cfg%phonon_fine_method) == PHONON_FINE_SUPERCELL) then
                call write_block_phonon_supercell_matrix(unit, cfg)
            end if
            if (trim(cfg%phonon_fine_method) /= PHONON_FINE_NONE) then
                if (trim(cfg%phonon_fine_qpoint_scheme) == PHONON_QPOINT_PATH) then
                    call write_phonon_path_block(unit, 'PHONON_FINE_KPOINT_PATH', &
                        cfg%phonon_fine_kpoint_path, 'PHONON_FINE_KPOINT_PATH_SPACING', &
                        cfg%phonon_fine_kpoint_path_spacing)
                else
                    call write_phonon_kpoint_mp(unit, 'PHONON_FINE_KPOINT_MP_GRID', &
                        cfg%phonon_fine_kpoint_mp_grid)
                end if
            end if
        end if

        close(unit)
    end subroutine write_cell_file

    ! ── PES scan support ──

    subroutine write_block_positions_frac(unit, cfg)
        !! Write %BLOCK POSITIONS_FRAC with fractional coordinates.
        !! Used by PES scan to place mobile atom at exact fractional grid points.
        integer, intent(in) :: unit
        type(castep_config_t), intent(in) :: cfg
        integer :: i

        if (cfg%num_atoms == 0) return

        write(unit, '(a)') '%BLOCK POSITIONS_FRAC'
        do i = 1, cfg%num_atoms
            write(unit, '(a, 3(1x, f12.8))') trim(cfg%atom_type(i)), &
                cfg%atom_x(i), cfg%atom_y(i), cfg%atom_z(i)
        end do
        write(unit, '(a)') '%ENDBLOCK POSITIONS_FRAC'
        write(unit, '(a)') ''
    end subroutine write_block_positions_frac


    subroutine write_block_ionic_constraints(unit, cfg, mobile_idx, fix_axes)
        !! Write %BLOCK IONIC_CONSTRAINTS for PES scan.
        !!
        !! CASTEP format: each constraint removes ONE degree of freedom.
        !!   constraint_num  element  elem_atom_idx  x_wt  y_wt  z_wt
        !! where elem_atom_idx counts atoms PER ELEMENT (not global).
        !!
        !! mobile_idx: 1-based GLOBAL index of the mobile atom.
        !! fix_axes(3): 1=fix, 0=free for x,y,z (e.g. [0,0,1] = fix z only).
        !! Non-mobile atoms are always fully fixed (3 constraints each).
        integer, intent(in) :: unit
        type(castep_config_t), intent(in) :: cfg
        integer, intent(in) :: mobile_idx
        integer, intent(in) :: fix_axes(3)
        integer :: i, dir, cnum, elem_count
        character(len=16) :: prev_elem

        if (cfg%num_atoms == 0) return

        write(unit, '(a)') '%BLOCK IONIC_CONSTRAINTS'
        cnum = 0
        prev_elem = ''
        elem_count = 0
        do i = 1, cfg%num_atoms
            ! Per-element atom index (CASTEP counts atoms within each element)
            if (trim(cfg%atom_type(i)) /= trim(prev_elem)) then
                elem_count = 1
                prev_elem = trim(cfg%atom_type(i))
            else
                elem_count = elem_count + 1
            end if

            if (i == mobile_idx) then
                ! Mobile atom: constrain only the directions in fix_axes
                do dir = 1, 3
                    if (fix_axes(dir) == 1) then
                        cnum = cnum + 1
                        call write_one_constraint(unit, cnum, &
                            trim(cfg%atom_type(i)), elem_count, dir)
                    end if
                end do
            else if (cfg%pes_constrain_others) then
                ! Other atoms: fully fixed (3 constraints: x, y, z)
                do dir = 1, 3
                    cnum = cnum + 1
                    call write_one_constraint(unit, cnum, &
                        trim(cfg%atom_type(i)), elem_count, dir)
                end do
            end if
        end do
        write(unit, '(a)') '%ENDBLOCK IONIC_CONSTRAINTS'
        write(unit, '(a)') ''
    end subroutine write_block_ionic_constraints

    subroutine write_one_constraint(unit, cnum, element, atom_idx, dir)
        !! Write a single constraint line fixing one Cartesian direction.
        integer, intent(in) :: unit, cnum, atom_idx, dir
        character(len=*), intent(in) :: element
        real(dp) :: w(3)
        w = [0.0_dp, 0.0_dp, 0.0_dp]
        w(dir) = 1.0_dp
        write(unit, '(i7, 1x, a, i8, 3(1x, f12.10))') &
            cnum, trim(element), atom_idx, w(1), w(2), w(3)
    end subroutine write_one_constraint

end module cell_writer
