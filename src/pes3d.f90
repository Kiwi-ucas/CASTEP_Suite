module pes3d
    !! 3D Potential Energy Surface scan — voxel grid generation, JSON output,
    !! and symmetry-based energy expansion.
    !!
    !! Generates .cell + .param input files for a 3D voxel grid of a mobile atom.
    !! Supports symmetry-aware scanning: only computes the asymmetric sub-volume
    !! around a reference atom, then uses space-group operations to expand to
    !! the full unit cell.
    use castep_config, only: dp, castep_config_t, cif_data_t, atom_t, sym_op_t, &
        IO_WRITE_FAIL, IO_PARSE_ERROR, IO_FILE_NOT_FOUND
    implicit none
    private

    public :: pes3d_grid_t
    public :: compute_local_grid_bounds
    public :: generate_pes3d_grid_points
    public :: write_pes3d_metadata_json
    public :: write_pes3d_cube
    public :: collect_pes3d_energies
    public :: symmetry_expand_energies

    ! Maximum grid size per dimension
    integer, parameter :: MAX_GRID_3D = 50

    type :: pes3d_grid_t
        integer  :: n_points(3)      = [5, 5, 5]     ! Nx, Ny, Nz
        real(dp) :: frac_range(3,2)  = 0.0_dp         ! [fmin,fmax] per axis
        integer  :: ref_atom_idx     = 1               ! 1-based index of reference atom
        character(len=8) :: scan_mode = 'SP'           ! 'SP' or 'RELAX'
        logical  :: use_symmetry     = .false.         ! whether to use space-group symmetry
        real(dp) :: ref_frac(3)      = 0.0_dp          ! reference atom's original fractional coords
        real(dp) :: half_dist(3)     = 0.0_dp          ! half-distance per axis (symmetry mode)
    end type pes3d_grid_t

contains

    ! ── Compute asymmetric sub-volume bounds ──

    subroutine compute_local_grid_bounds(ref_idx, atoms, n_atoms, sym_ops, n_symops, &
                                          half_dist, iostat, iomsg)
        !! For the selected reference atom, compute the half-distance along each
        !! fractional axis to the nearest space-group-equivalent atom.
        !! Uses minimum-image convention (|dx| with periodicity 1.0).
        !! Clamps to [1/MAX_GRID_3D, 0.3] as a safety margin.
        integer, intent(in) :: ref_idx             ! 1-based index of reference atom
        type(atom_t), intent(in) :: atoms(:)
        integer, intent(in) :: n_atoms
        type(sym_op_t), intent(in) :: sym_ops(:)
        integer, intent(in) :: n_symops
        real(dp), intent(out) :: half_dist(3)      ! output half-distance per axis
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        real(dp) :: frac0(3), new_frac(3), delta(3)
        real(dp) :: min_dist(3)
        integer :: iop
        logical :: found_other
        real(dp), parameter :: TOL = 1.0e-5_dp
        real(dp), parameter :: MAX_HALF = 0.3_dp

        iostat = 0
        half_dist = MAX_HALF  ! default large value

        if (ref_idx < 1 .or. ref_idx > n_atoms) then
            iostat = 1
            if (present(iomsg)) iomsg = 'Reference atom index out of range'
            return
        end if

        frac0 = [atoms(ref_idx)%x, atoms(ref_idx)%y, atoms(ref_idx)%z]
        min_dist = 1.0_dp  ! large initial (larger than any minimum-image distance)

        if (n_symops <= 1) return  ! P1: no symmetry → keep default MAX_HALF

        do iop = 1, n_symops
            ! Apply symmetry operation: new_frac = rot * frac0 + trans
            new_frac(1) = sym_ops(iop)%rot(1,1)*frac0(1) &
                        + sym_ops(iop)%rot(1,2)*frac0(2) &
                        + sym_ops(iop)%rot(1,3)*frac0(3) &
                        + sym_ops(iop)%trans(1)
            new_frac(2) = sym_ops(iop)%rot(2,1)*frac0(1) &
                        + sym_ops(iop)%rot(2,2)*frac0(2) &
                        + sym_ops(iop)%rot(2,3)*frac0(3) &
                        + sym_ops(iop)%trans(2)
            new_frac(3) = sym_ops(iop)%rot(3,1)*frac0(1) &
                        + sym_ops(iop)%rot(3,2)*frac0(2) &
                        + sym_ops(iop)%rot(3,3)*frac0(3) &
                        + sym_ops(iop)%trans(3)

            ! Wrap to [0, 1)
            call wrap_to_unit(new_frac(1))
            call wrap_to_unit(new_frac(2))
            call wrap_to_unit(new_frac(3))

            ! Minimum-image displacement (periodic distance)
            delta(1) = new_frac(1) - frac0(1)
            delta(2) = new_frac(2) - frac0(2)
            delta(3) = new_frac(3) - frac0(3)
            ! Minimum-image: map to [-0.5, 0.5)
            delta = delta - anint(delta)

            ! Check if this is a DIFFERENT position (not identity)
            found_other = abs(delta(1)) > TOL .or. abs(delta(2)) > TOL .or. abs(delta(3)) > TOL
            if (found_other) then
                ! Update minimum absolute distance per axis
                if (abs(delta(1)) > TOL .and. abs(delta(1)) < min_dist(1)) &
                    min_dist(1) = abs(delta(1))
                if (abs(delta(2)) > TOL .and. abs(delta(2)) < min_dist(2)) &
                    min_dist(2) = abs(delta(2))
                if (abs(delta(3)) > TOL .and. abs(delta(3)) < min_dist(3)) &
                    min_dist(3) = abs(delta(3))
            end if
        end do

        ! Half-distance = half of minimum separation, clamped
        select case (3)
        case (1)  ! dummy for dimension loop
        end select
        if (min_dist(1) < 1.0_dp) half_dist(1) = max(1.0_dp/MAX_GRID_3D, &
            min(MAX_HALF, min_dist(1) * 0.45_dp))
        if (min_dist(2) < 1.0_dp) half_dist(2) = max(1.0_dp/MAX_GRID_3D, &
            min(MAX_HALF, min_dist(2) * 0.45_dp))
        if (min_dist(3) < 1.0_dp) half_dist(3) = max(1.0_dp/MAX_GRID_3D, &
            min(MAX_HALF, min_dist(3) * 0.45_dp))

    end subroutine compute_local_grid_bounds


    ! ── 3D grid point generation ──

    subroutine generate_pes3d_grid_points(grid, frac_points, n_total, iostat, iomsg)
        !! Generate all fractional coordinate triplets for the 3D voxel grid.
        !! Returns frac_points(n_total, 3) — columns are (fx, fy, fz).
        !! Order: x fastest, then y, then z (matches Fortran inner-loop convention).
        type(pes3d_grid_t), intent(in) :: grid
        real(dp), allocatable, intent(out) :: frac_points(:,:)
        integer, intent(out) :: n_total
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        integer :: i, j, k, idx
        real(dp) :: dx, dy, dz, fx, fy, fz

        iostat = 0
        n_total = 0

        if (any(grid%n_points < 2)) then
            iostat = 1
            if (present(iomsg)) iomsg = 'Grid must have at least 2 points per dimension'
            return
        end if
        if (any(grid%n_points > MAX_GRID_3D)) then
            iostat = 2
            if (present(iomsg)) iomsg = 'Grid size exceeds maximum'
            return
        end if

        n_total = grid%n_points(1) * grid%n_points(2) * grid%n_points(3)
        allocate(frac_points(n_total, 3), stat=iostat)
        if (iostat /= 0) then
            if (present(iomsg)) iomsg = 'Memory allocation failed for grid points'
            return
        end if

        ! Compute step sizes
        dx = 0.0_dp; dy = 0.0_dp; dz = 0.0_dp
        if (grid%n_points(1) > 1) dx = (grid%frac_range(1,2) - grid%frac_range(1,1)) &
                                       / real(grid%n_points(1) - 1, dp)
        if (grid%n_points(2) > 1) dy = (grid%frac_range(2,2) - grid%frac_range(2,1)) &
                                       / real(grid%n_points(2) - 1, dp)
        if (grid%n_points(3) > 1) dz = (grid%frac_range(3,2) - grid%frac_range(3,1)) &
                                       / real(grid%n_points(3) - 1, dp)

        idx = 0
        do k = 0, grid%n_points(3) - 1
            fz = grid%frac_range(3,1) + k * dz
            do j = 0, grid%n_points(2) - 1
                fy = grid%frac_range(2,1) + j * dy
                do i = 0, grid%n_points(1) - 1
                    fx = grid%frac_range(1,1) + i * dx
                    idx = idx + 1
                    frac_points(idx, 1) = fx
                    frac_points(idx, 2) = fy
                    frac_points(idx, 3) = fz
                end do
            end do
        end do
    end subroutine generate_pes3d_grid_points


    ! ── JSON metadata output ──

    subroutine write_pes3d_metadata_json(json_path, grid, cfg, cif, iostat, iomsg)
        !! Write pes3d_metadata.json describing the 3D scan grid, lattice,
        !! atom info, and (if use_symmetry) symmetry operations.
        !! Energies are written as null placeholders.
        !! Symmetry operations are written at generation time so the collect
        !! step can work independently of the original CIF file.
        character(len=*), intent(in) :: json_path
        type(pes3d_grid_t), intent(in) :: grid
        type(castep_config_t), intent(in), target :: cfg
        type(cif_data_t), intent(in) :: cif
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        integer :: unit, ios, i, n_total, iop
        character(len=256) :: rot_str(3)

        iostat = 0

        open(newunit=unit, file=trim(json_path), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            iostat = IO_WRITE_FAIL
            if (present(iomsg)) iomsg = 'Cannot write: ' // trim(json_path)
            return
        end if

        write(unit, '(a)') '{'
        write(unit, '(a)') '  "type": "pes_scan_3d",'
        write(unit, '(a,i0,a,i0,a,i0,a)') '  "nx": ', grid%n_points(1), &
            ', "ny": ', grid%n_points(2), ', "nz": ', grid%n_points(3), ','

        write(unit, '(a, f12.8, a, f12.8, a)') &
            '  "fx_range": [', grid%frac_range(1,1), ', ', grid%frac_range(1,2), '],'
        write(unit, '(a, f12.8, a, f12.8, a)') &
            '  "fy_range": [', grid%frac_range(2,1), ', ', grid%frac_range(2,2), '],'
        write(unit, '(a, f12.8, a, f12.8, a)') &
            '  "fz_range": [', grid%frac_range(3,1), ', ', grid%frac_range(3,2), '],'

        write(unit, '(a,i0,a,a,a)') '  "mobile_atom": { "index": ', grid%ref_atom_idx - 1, &
            ', "element": "', trim(cfg%atom_type(grid%ref_atom_idx)), '" },'
        write(unit, '(a, f10.6, a, f10.6, a, f10.6, a)') &
            '  "ref_atom_fractional": [', grid%ref_frac(1), ', ', &
            grid%ref_frac(2), ', ', grid%ref_frac(3), '],'
        write(unit, '(a,a,a)') '  "scan_mode": "', trim(grid%scan_mode), '",'

        write(unit, '(a, f10.6, a, f10.6, a, f10.6, a)') &
            '  "lattice": { "a": ', cfg%cell_length(1), &
            ', "b": ', cfg%cell_length(2), ', "c": ', cfg%cell_length(3), ','
        write(unit, '(a, f10.6, a, f10.6, a, f10.6, a)') &
            '    "alpha": ', cfg%cell_angle(1), &
            ', "beta": ', cfg%cell_angle(2), ', "gamma": ', cfg%cell_angle(3), ' },'

        ! Structure atoms (fractional coords)
        write(unit, '(a)') '  "structure_atoms": ['
        do i = 1, cfg%num_atoms
            write(unit, '(a,a,a, f10.6, a, f10.6, a, f10.6, a)', advance='no') &
                '    { "element": "', trim(cfg%atom_type(i)), &
                '", "fx": ', cfg%atom_x(i), ', "fy": ', cfg%atom_y(i), &
                ', "fz": ', cfg%atom_z(i), ' }'
            if (i < cfg%num_atoms) then
                write(unit, '(a)') ','
            else
                write(unit, '(a)') ''
            end if
        end do
        write(unit, '(a)') '  ],'

        ! Symmetry info (written at generation time — self-contained JSON)
        write(unit, '(a,l1,a)') '  "use_symmetry": ', grid%use_symmetry, ','
        if (grid%use_symmetry .and. cif%n_symops > 0) then
            write(unit, '(a,a,a)') '  "space_group": "', &
                trim(cif%space_group_name), '",'
            write(unit, '(a,i0,a)') '  "n_symops": ', cif%n_symops, ','
            write(unit, '(a)') '  "sym_ops": ['
            do iop = 1, cif%n_symops
                write(rot_str(1), '(a,i2,a,i2,a,i2,a)') '[', &
                    cif%sym_ops(iop)%rot(1,1), ',', &
                    cif%sym_ops(iop)%rot(1,2), ',', &
                    cif%sym_ops(iop)%rot(1,3), ']'
                write(rot_str(2), '(a,i2,a,i2,a,i2,a)') '[', &
                    cif%sym_ops(iop)%rot(2,1), ',', &
                    cif%sym_ops(iop)%rot(2,2), ',', &
                    cif%sym_ops(iop)%rot(2,3), ']'
                write(rot_str(3), '(a,i2,a,i2,a,i2,a)') '[', &
                    cif%sym_ops(iop)%rot(3,1), ',', &
                    cif%sym_ops(iop)%rot(3,2), ',', &
                    cif%sym_ops(iop)%rot(3,3), ']'
                write(unit, '(a,a,a,a,a,a,a,f14.10,a,f14.10,a,f14.10,a)', advance='no') &
                    '    { "rot": [', &
                    trim(adjustl(rot_str(1))), ', ', &
                    trim(adjustl(rot_str(2))), ', ', &
                    trim(adjustl(rot_str(3))), &
                    '], "trans": [', &
                    cif%sym_ops(iop)%trans(1), ', ', &
                    cif%sym_ops(iop)%trans(2), ', ', &
                    cif%sym_ops(iop)%trans(3), '] }'
                if (iop < cif%n_symops) then
                    write(unit, '(a)') ','
                else
                    write(unit, '(a)') ''
                end if
            end do
            write(unit, '(a)') '  ],'
        else
            write(unit, '(a)') '  "space_group": "P1",'
            write(unit, '(a)') '  "n_symops": 0,'
            write(unit, '(a)') '  "sym_ops": [],'
        end if

        ! Energy grid: write null entries
        write(unit, '(a)') '  "energies": ['
        n_total = grid%n_points(1) * grid%n_points(2) * grid%n_points(3)
        do i = 1, n_total
            if (i < n_total) then
                write(unit, '(a)') '    null,'
            else
                write(unit, '(a)') '    null'
            end if
        end do
        write(unit, '(a)') '  ],'
        write(unit, '(a)') '  "energies_expanded": [],'
        write(unit, '(a)') '  "has_energies": false,'
        write(unit, '(a)') '  "has_expanded": false'
        write(unit, '(a)') '}'

        close(unit)
    end subroutine write_pes3d_metadata_json


    ! ── Cube file output (Gaussian Cube format) ──

    subroutine write_pes3d_cube(cube_path, grid, cfg, cif, energies, n_energies, iostat, iomsg)
        !! Write a Gaussian Cube file for the 3D PES scalar field.
        !! Line 1: description, Line 2: compact JSON metadata.
        !! Compatible with VMD, VESTA, Avogadro.
        character(len=*), intent(in) :: cube_path
        type(pes3d_grid_t), intent(in) :: grid
        type(castep_config_t), intent(in), target :: cfg
        type(cif_data_t), intent(in) :: cif
        real(dp), intent(in) :: energies(:)
        integer, intent(in) :: n_energies
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        real(dp) :: lattice_vecs(3,3), sp_vecs(3,3), cart(3)
        real(dp) :: a, b, c, alpha, beta, gamma
        integer :: unit, ios, i, nx, ny, nz
        real(dp) :: e_val
        character(len=256) :: meta_json
        integer :: z_num
        character(len=8) :: el

        iostat = 0
        nx = grid%n_points(1); ny = grid%n_points(2); nz = grid%n_points(3)

        a = cfg%cell_length(1); b = cfg%cell_length(2); c = cfg%cell_length(3)
        alpha = cfg%cell_angle(1); beta = cfg%cell_angle(2); gamma = cfg%cell_angle(3)

        ! Compute lattice vectors (a along x, b in xy-plane, c general)
        alpha = alpha * 3.141592653589793_dp / 180.0_dp
        beta  = beta  * 3.141592653589793_dp / 180.0_dp
        gamma = gamma * 3.141592653589793_dp / 180.0_dp
        lattice_vecs(1,1) = a
        lattice_vecs(2,1) = 0.0_dp
        lattice_vecs(3,1) = 0.0_dp
        lattice_vecs(1,2) = b * cos(gamma)
        lattice_vecs(2,2) = b * sin(gamma)
        lattice_vecs(3,2) = 0.0_dp
        lattice_vecs(1,3) = c * cos(beta)
        lattice_vecs(2,3) = c * (cos(alpha) - cos(beta)*cos(gamma)) / max(sin(gamma), 1.0e-10_dp)
        lattice_vecs(3,3) = sqrt(max(c*c - lattice_vecs(1,3)**2 - lattice_vecs(2,3)**2, 0.0_dp))

        ! Voxel spacing vectors = lattice / (N-1)
        sp_vecs(:,1) = lattice_vecs(:,1) / max(1, nx-1)
        sp_vecs(:,2) = lattice_vecs(:,2) / max(1, ny-1)
        sp_vecs(:,3) = lattice_vecs(:,3) / max(1, nz-1)

        ! Compact metadata JSON
        write(meta_json, '(a,l1,a,a,a,i0,a)') &
            '{"use_symmetry":', grid%use_symmetry, &
            ',"scan_mode":"', trim(grid%scan_mode), &
            '","n_symops":', cif%n_symops, '}'

        open(newunit=unit, file=trim(cube_path), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            iostat = IO_WRITE_FAIL
            if (present(iomsg)) iomsg = 'Cannot write: ' // trim(cube_path)
            return
        end if

        ! Line 1: description
        write(unit, '(a,i0,a,i0,a,i0)') '3D PES: ', nx, 'x', ny, 'x', nz
        ! Line 2: metadata
        write(unit, '(a)') trim(meta_json)
        ! Lines 3-6: header
        write(unit, '(i5,3f12.6)') cfg%num_atoms, 0.0_dp, 0.0_dp, 0.0_dp  ! natom, origin=0
        write(unit, '(i5,3f12.6)') nx, sp_vecs(1,1), sp_vecs(2,1), sp_vecs(3,1)
        write(unit, '(i5,3f12.6)') ny, sp_vecs(1,2), sp_vecs(2,2), sp_vecs(3,2)
        write(unit, '(i5,3f12.6)') nz, sp_vecs(1,3), sp_vecs(2,3), sp_vecs(3,3)
        ! Atom lines: Z, charge, x, y, z (Cartesian)
        do i = 1, cfg%num_atoms
            el = trim(cfg%atom_type(i))
            z_num = element_to_z(el)
            ! Fractional → Cartesian
            cart(:) = cfg%atom_x(i)*lattice_vecs(:,1) &
                    + cfg%atom_y(i)*lattice_vecs(:,2) &
                    + cfg%atom_z(i)*lattice_vecs(:,3)
            write(unit, '(i5,4f12.6)') z_num, 0.0_dp, cart(1), cart(2), cart(3)
        end do
        ! Volumetric data (6 per line)
        do i = 1, n_energies
            e_val = energies(i)
            write(unit, '(es13.5)', advance='no') e_val
            if (mod(i, 6) == 0 .or. i == n_energies) write(unit, '(a)') ''
        end do

        close(unit)
        write(*, '(a)') '  Cube file written: ' // trim(cube_path)
    end subroutine write_pes3d_cube


    ! Element symbol → atomic number (common elements for battery materials)
    integer function element_to_z(el) result(z)
        character(len=*), intent(in) :: el
        character(len=8) :: uel
        uel = adjustl(el)
        ! Uppercase first char
        if (uel(1:1) >= 'a' .and. uel(1:1) <= 'z') uel(1:1) = achar(iachar(uel(1:1)) - 32)
        select case (trim(uel))
        case ('H');  z = 1
        case ('LI'); z = 3
        case ('C');  z = 6
        case ('N');  z = 7
        case ('O');  z = 8
        case ('F');  z = 9
        case ('NA'); z = 11
        case ('MG'); z = 12
        case ('AL'); z = 13
        case ('SI'); z = 14
        case ('P');  z = 15
        case ('S');  z = 16
        case ('CL'); z = 17
        case ('K');  z = 19
        case ('CA'); z = 20
        case ('TI'); z = 22
        case ('V');  z = 23
        case ('CR'); z = 24
        case ('MN'); z = 25
        case ('FE'); z = 26
        case ('CO'); z = 27
        case ('NI'); z = 28
        case ('CU'); z = 29
        case ('ZN'); z = 30
        case ('GE'); z = 32
        case ('AS'); z = 33
        case ('SE'); z = 34
        case ('BR'); z = 35
        case ('ZR'); z = 40
        case ('NB'); z = 41
        case ('MO'); z = 42
        case ('AG'); z = 47
        case ('SN'); z = 50
        case ('SB'); z = 51
        case ('TE'); z = 52
        case ('I');  z = 53
        case ('BA'); z = 56
        case ('W');  z = 74
        case ('PT'); z = 78
        case ('AU'); z = 79
        case ('PB'); z = 82
        case ('BI'); z = 83
        case default; z = 0
        end select
    end function element_to_z


    ! ── Result collection ──

    subroutine collect_pes3d_energies(scan_dir, iostat, iomsg)
        !! Scan all grid_i_j_k/ subdirectories, parse .castep files, fill
        !! energies into pes3d_metadata.json. If use_symmetry is true,
        !! also run symmetry expansion.
        character(len=*), intent(in) :: scan_dir
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        character(len=1024) :: json_path, grid_dir, castep_file, line
        integer :: nx, ny, nz, n_total, collected, missing, ios, json_unit
        integer :: i, j, k, idx
        real(dp), allocatable :: energies(:)
        logical, allocatable :: has_energy(:)
        logical :: exists, use_sym
        real(dp) :: e_val, e_min, e_max

        iostat = 0
        collected = 0; missing = 0
        e_min = huge(1.0_dp); e_max = -huge(1.0_dp)
        use_sym = .false.

        json_path = trim(scan_dir) // '/pes3d_metadata.json'
        inquire(file=trim(json_path), exist=exists)
        if (.not. exists) then
            iostat = IO_FILE_NOT_FOUND
            if (present(iomsg)) iomsg = 'pes3d_metadata.json not found in: ' // trim(scan_dir)
            return
        end if

        ! Read nx, ny, nz, use_symmetry from JSON (simple line-based parsing)
        nx = 0; ny = 0; nz = 0
        open(newunit=json_unit, file=trim(json_path), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            iostat = IO_PARSE_ERROR
            if (present(iomsg)) iomsg = 'Cannot open: ' // trim(json_path)
            return
        end if
        do
            read(json_unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, '"nx"') > 0 .and. nx <= 0) &
                call extract_json_int(line, '"nx"', nx, ios)
            if (index(line, '"ny"') > 0 .and. ny <= 0) &
                call extract_json_int(line, '"ny"', ny, ios)
            if (index(line, '"nz"') > 0 .and. nz <= 0) &
                call extract_json_int(line, '"nz"', nz, ios)
            if (index(line, '"use_symmetry"') > 0) &
                use_sym = index(line, 'true') > 0
        end do
        close(json_unit)

        if (nx <= 0 .or. ny <= 0 .or. nz <= 0) then
            iostat = IO_PARSE_ERROR
            if (present(iomsg)) iomsg = 'Failed to read grid dimensions from JSON'
            return
        end if

        n_total = nx * ny * nz
        allocate(energies(n_total), stat=ios)
        allocate(has_energy(n_total), stat=iostat)
        if (ios /= 0 .or. iostat /= 0) then
            iostat = 1
            if (present(iomsg)) iomsg = 'Memory allocation failed'
            return
        end if
        energies = 0.0_dp
        has_energy = .false.

        ! Scan grid directories (x fastest, then y, then z)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    write(grid_dir, '(a, a, i3.3, a, i3.3, a, i3.3)') &
                        trim(scan_dir), '/grid_', i, '_', j, '_', k
                    castep_file = find_castep_in_dir(grid_dir)
                    if (len_trim(castep_file) == 0) then
                        missing = missing + 1
                        cycle
                    end if

                    idx = (k-1)*nx*ny + (j-1)*nx + i
                    e_val = parse_castep_energy(castep_file, ios)
                    if (ios == 0) then
                        energies(idx) = e_val
                        has_energy(idx) = .true.
                        collected = collected + 1
                        if (e_val < e_min) e_min = e_val
                        if (e_val > e_max) e_max = e_val
                    else
                        missing = missing + 1
                    end if
                end do
            end do
        end do

        write(*, '(a)') '  ── Collection Summary ──'
        write(*, '(a, i0, a, i0)') '  Collected: ', collected, ' / ', n_total
        if (missing > 0) write(*, '(a, i0)') '  Missing:   ', missing
        if (collected > 0) then
            write(*, '(a, f18.8)') '  E min (eV): ', e_min
            write(*, '(a, f18.8)') '  E max (eV): ', e_max
        end if

        ! Write updated JSON with energies
        call rewrite_json_with_energies_3d(json_path, energies, has_energy, &
            nx, ny, nz, use_sym, ios)
        if (ios /= 0) then
            iostat = IO_WRITE_FAIL
            if (present(iomsg)) iomsg = 'Failed to rewrite JSON with energies'
            deallocate(energies); deallocate(has_energy)
            return
        end if

        deallocate(energies)
        deallocate(has_energy)
    end subroutine collect_pes3d_energies


    ! ── Symmetry expansion ──

    subroutine symmetry_expand_energies(json_path, iostat, iomsg)
        !! Read energies and sym_ops from pes3d_metadata.json, apply symmetry
        !! operations to expand the local energy grid to the full [0,1)^3 cell.
        !! Uses FORWARD symmetry ops: each local grid point is mapped by each
        !! sym op to the expanded cell. Overlapping cells take minimum energy.
        character(len=*), intent(in) :: json_path
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        integer, parameter :: MAX_SYM_OPS = 256
        integer, parameter :: MAX_EXP = 100   ! max expanded grid points per dimension

        ! Local grid data
        integer :: nx, ny, nz
        real(dp) :: fx_range(2), fy_range(2), fz_range(2)
        real(dp), allocatable :: local_energies(:)
        integer :: n_local

        ! Symmetry operations
        integer :: n_symops
        integer :: rot(3,3,MAX_SYM_OPS)
        real(dp) :: trans(3,MAX_SYM_OPS)
        logical :: has_sym

        ! Expanded grid
        integer :: exp_nx, exp_ny, exp_nz, n_exp
        real(dp), allocatable :: exp_energies(:)
        logical, allocatable :: exp_filled(:)
        real(dp) :: sp_x, sp_y, sp_z
        real(dp) :: local_frac(3), exp_frac(3)
        real(dp) :: e_min_final, e_max_final

        character(len=1024) :: line, tmp_path
        integer :: unit_in, unit_out, ios, json_ios
        integer :: i, j, k, idx, iop, iexp
        logical :: exists, in_energies, wrote_exp, use_sym
        integer :: ei, ej, ek
        real(dp) :: el

        iostat = 0
        has_sym = .false.
        n_symops = 0
        rot = 0
        trans = 0.0_dp
        nx = 0; ny = 0; nz = 0
        fx_range = [0.0_dp, 1.0_dp]
        fy_range = [0.0_dp, 1.0_dp]
        fz_range = [0.0_dp, 1.0_dp]
        use_sym = .false.

        ! ── Parse JSON header: nx, ny, nz, frac_range, use_symmetry ──
        inquire(file=trim(json_path), exist=exists)
        if (.not. exists) then
            iostat = IO_FILE_NOT_FOUND
            if (present(iomsg)) iomsg = 'JSON not found'
            return
        end if

        open(newunit=unit_in, file=trim(json_path), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            iostat = IO_PARSE_ERROR
            if (present(iomsg)) iomsg = 'Cannot open JSON'
            return
        end if

        do
            read(unit_in, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, '"nx"') > 0 .and. index(line, '"n_symops"') == 0) &
                call extract_json_int(line, '"nx"', nx, ios)
            if (index(line, '"ny"') > 0) call extract_json_int(line, '"ny"', ny, ios)
            if (index(line, '"nz"') > 0) call extract_json_int(line, '"nz"', nz, ios)
            if (index(line, '"use_symmetry"') > 0) &
                use_sym = index(line, 'true') > 0
            if (index(line, '"fx_range"') > 0) &
                call extract_json_two_reals(line, fx_range(1), fx_range(2))
            if (index(line, '"fy_range"') > 0) &
                call extract_json_two_reals(line, fy_range(1), fy_range(2))
            if (index(line, '"fz_range"') > 0) &
                call extract_json_two_reals(line, fz_range(1), fz_range(2))
        end do
        close(unit_in)

        if (.not. use_sym) then
            if (present(iomsg)) iomsg = 'No symmetry — skipping expansion'
            return  ! nothing to do
        end if

        if (nx <= 0 .or. ny <= 0 .or. nz <= 0) then
            iostat = IO_PARSE_ERROR
            if (present(iomsg)) iomsg = 'Invalid grid dimensions'
            return
        end if

        ! ── Parse sym_ops from JSON ──
        n_symops = 0
        open(newunit=unit_in, file=trim(json_path), status='old', action='read', iostat=ios)
        do
            read(unit_in, '(a)', iostat=ios) line
            if (ios /= 0) exit
            ! Count sym_ops: look for "rot" lines
            if (index(line, '"rot"') > 0) n_symops = n_symops + 1
        end do
        close(unit_in)

        if (n_symops < 1 .or. n_symops > MAX_SYM_OPS) then
            if (present(iomsg)) iomsg = 'No symmetry ops found in JSON'
            return
        end if

        ! Parse each sym_op
        open(newunit=unit_in, file=trim(json_path), status='old', action='read', iostat=ios)
        iop = 0
        do
            read(unit_in, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, '"rot"') > 0) then
                iop = iop + 1
                if (iop > MAX_SYM_OPS) exit
                ! Parse the rot matrix from this line and following 2 lines
                call parse_symop_line(line, rot(:,:,iop))
                ! read next 2 lines for rows 2 and 3
                read(unit_in, '(a)', iostat=ios) line  ! row 2
                if (index(line, '[') > 0) &
                    call parse_symop_matrix_row(line, rot(2,:,iop))
                read(unit_in, '(a)', iostat=ios) line  ! row 3
                if (index(line, '[') > 0) &
                    call parse_symop_matrix_row(line, rot(3,:,iop))
            end if
            if (index(line, '"trans"') > 0 .and. iop > 0) then
                ! Parse the trans vector from the same line as rot
                call parse_symop_trans_from_line(line, trans(:,iop))
                ! If trans not on rot line, it might be on the rot row 3 line
                if (all(abs(trans(:,iop)) < 1.0e-10_dp)) then
                    ! Try to find trans on the current line
                    call parse_trans_from_json_line(line, trans(:,iop))
                end if
            end if
        end do
        close(unit_in)

        ! Fallback: if trans parsing failed, parse from lines containing "trans"
        if (all(abs(trans(:,1)) < 1.0e-10_dp)) then
            open(newunit=unit_in, file=trim(json_path), status='old', action='read', iostat=ios)
            iop = 0
            do
                read(unit_in, '(a)', iostat=ios) line
                if (ios /= 0) exit
                if (index(line, '"rot"') > 0) iop = iop + 1
                if (index(line, '"trans"') > 0 .and. iop > 0 .and. iop <= n_symops) then
                    call parse_trans_from_json_line(line, trans(:,iop))
                end if
            end do
            close(unit_in)
        end if

        has_sym = (n_symops > 0)

        if (.not. has_sym) then
            if (present(iomsg)) iomsg = 'Failed to parse sym_ops from JSON'
            iostat = IO_PARSE_ERROR
            return
        end if

        ! ── Read local energies from JSON ──
        n_local = nx * ny * nz
        allocate(local_energies(n_local), stat=ios)
        if (ios /= 0) then
            iostat = 1; return
        end if
        local_energies = 0.0_dp

        open(newunit=unit_in, file=trim(json_path), status='old', action='read', iostat=ios)
        in_energies = .false.
        i = 0
        do
            read(unit_in, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, '"energies"') > 0 .and. &
                index(line, '"energies_expanded"') == 0) then
                in_energies = .true.
                cycle
            end if
            if (in_energies) then
                if (index(line, ']') > 0) exit
                if (index(line, 'null') > 0) then
                    i = i + 1
                    local_energies(i) = 0.0_dp  ! null -> 0 (will be skipped)
                    cycle
                end if
                read(line, *, iostat=json_ios) el
                if (json_ios == 0) then
                    i = i + 1
                    local_energies(i) = el
                end if
            end if
        end do
        close(unit_in)

        if (i < n_local) then
            iostat = IO_PARSE_ERROR
            if (present(iomsg)) iomsg = 'Too few energy values in JSON'
            deallocate(local_energies); return
        end if

        ! ── Build expanded grid ──
        sp_x = (fx_range(2) - fx_range(1)) / max(1, nx - 1)
        sp_y = (fy_range(2) - fy_range(1)) / max(1, ny - 1)
        sp_z = (fz_range(2) - fz_range(1)) / max(1, nz - 1)

        exp_nx = min(MAX_EXP, nint(1.0_dp / sp_x) + 1)
        exp_ny = min(MAX_EXP, nint(1.0_dp / sp_y) + 1)
        exp_nz = min(MAX_EXP, nint(1.0_dp / sp_z) + 1)
        n_exp = exp_nx * exp_ny * exp_nz

        allocate(exp_energies(n_exp), stat=ios)
        allocate(exp_filled(n_exp), stat=iostat)
        if (ios /= 0 .or. iostat /= 0) then
            iostat = 1
            deallocate(local_energies); return
        end if
        exp_energies = huge(1.0_dp)
        exp_filled = .false.

        ! ── Apply FORWARD symmetry ops to each local grid point ──
        do iop = 1, n_symops
            do k = 0, nz - 1
                local_frac(3) = fz_range(1) + k * sp_z
                do j = 0, ny - 1
                    local_frac(2) = fy_range(1) + j * sp_y
                    do i = 0, nx - 1
                        local_frac(1) = fx_range(1) + i * sp_x

                        ! Forward operation: exp_frac = R * local_frac + trans
                        exp_frac(1) = rot(1,1,iop)*local_frac(1) &
                                    + rot(1,2,iop)*local_frac(2) &
                                    + rot(1,3,iop)*local_frac(3) + trans(1,iop)
                        exp_frac(2) = rot(2,1,iop)*local_frac(1) &
                                    + rot(2,2,iop)*local_frac(2) &
                                    + rot(2,3,iop)*local_frac(3) + trans(2,iop)
                        exp_frac(3) = rot(3,1,iop)*local_frac(1) &
                                    + rot(3,2,iop)*local_frac(2) &
                                    + rot(3,3,iop)*local_frac(3) + trans(3,iop)

                        call wrap_to_unit(exp_frac(1))
                        call wrap_to_unit(exp_frac(2))
                        call wrap_to_unit(exp_frac(3))

                        ! Map to nearest expanded grid cell
                        ei = nint(exp_frac(1) / sp_x)
                        ej = nint(exp_frac(2) / sp_y)
                        ek = nint(exp_frac(3) / sp_z)

                        ! Wrap indices to [0, exp_n-1]
                        ei = modulo(ei, exp_nx)
                        ej = modulo(ej, exp_ny)
                        ek = modulo(ek, exp_nz)

                        idx = k * nx * ny + j * nx + i + 1  ! 1-based index
                        iexp = ek * exp_nx * exp_ny + ej * exp_nx + ei + 1

                        if (iexp >= 1 .and. iexp <= n_exp) then
                            if (.not. exp_filled(iexp) .or. &
                                local_energies(idx) < exp_energies(iexp)) then
                                exp_energies(iexp) = local_energies(idx)
                                exp_filled(iexp) = .true.
                            end if
                        end if
                    end do
                end do
            end do
        end do

        ! ── Write expanded energies to JSON ──
        tmp_path = trim(json_path) // '.tmp'

        open(newunit=unit_in, file=trim(json_path), status='old', action='read', iostat=ios)
        open(newunit=unit_out, file=trim(tmp_path), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            close(unit_in); deallocate(local_energies)
            deallocate(exp_energies); deallocate(exp_filled); return
        end if

        wrote_exp = .false.
        in_energies = .false.
        e_min_final = huge(1.0_dp)
        e_max_final = -huge(1.0_dp)
        do iexp = 1, n_exp
            if (exp_filled(iexp)) then
                if (exp_energies(iexp) < e_min_final) e_min_final = exp_energies(iexp)
                if (exp_energies(iexp) > e_max_final) e_max_final = exp_energies(iexp)
            end if
        end do

        do
            read(unit_in, '(a)', iostat=ios) line
            if (ios /= 0) exit

            if (index(line, '"energies_expanded"') > 0 .and. .not. wrote_exp) then
                write(unit_out, '(a)') '  "energies_expanded": ['
                do iexp = 1, n_exp
                    if (exp_filled(iexp)) then
                        if (iexp < n_exp) then
                            write(unit_out, '(a, f20.8, a)') '    ', exp_energies(iexp), ','
                        else
                            write(unit_out, '(a, f20.8)') '    ', exp_energies(iexp)
                        end if
                    else
                        if (iexp < n_exp) then
                            write(unit_out, '(a)') '    null,'
                        else
                            write(unit_out, '(a)') '    null'
                        end if
                    end if
                end do
                write(unit_out, '(a)') '  ],'
                wrote_exp = .true.
                cycle
            end if

            if (index(line, '"has_expanded"') > 0) then
                write(unit_out, '(a)') '  "has_expanded": true'
                cycle
            end if

            ! Write expanded grid dimensions
            if (index(line, '"exp_nx"') > 0) then
                write(unit_out, '(a,i0,a)') '  "exp_nx": ', exp_nx, ','
                cycle
            end if
            if (index(line, '"exp_ny"') > 0) then
                write(unit_out, '(a,i0,a)') '  "exp_ny": ', exp_ny, ','
                cycle
            end if
            if (index(line, '"exp_nz"') > 0) then
                write(unit_out, '(a,i0,a)') '  "exp_nz": ', exp_nz, ','
                cycle
            end if

            ! Update energy range
            if (index(line, '"e_min"') > 0) then
                write(unit_out, '(a, f20.8, a)') '  "e_min": ', e_min_final, ','
                cycle
            end if
            if (index(line, '"e_max"') > 0) then
                write(unit_out, '(a, f20.8, a)') '  "e_max": ', e_max_final, ','
                cycle
            end if

            write(unit_out, '(a)') trim(line)
        end do

        close(unit_in)
        close(unit_out)

        call execute_command_line('mv "' // trim(tmp_path) // '" "' // trim(json_path) // '"', &
            exitstat=ios)

        deallocate(local_energies)
        deallocate(exp_energies)
        deallocate(exp_filled)

        write(*, '(a)') '  ── Symmetry Expansion ──'
        write(*, '(a, i0, a, i0, a, i0)') '  Expanded grid: ', &
            exp_nx, ' x ', exp_ny, ' x ', exp_nz, ' = ', n_exp
        if (e_min_final < huge(1.0_dp)) then
            write(*, '(a, f18.8)') '  E min (expanded): ', e_min_final
            write(*, '(a, f18.8)') '  E max (expanded): ', e_max_final
        end if
    end subroutine symmetry_expand_energies


    ! ── Private helpers ──

    function find_castep_in_dir(dir_path) result(castep_path)
        !! Find a .castep file: hardcoded scan.castep first, then wildcard ls.
        character(len=*), intent(in) :: dir_path
        character(len=1024) :: castep_path, cmd, tmp_file
        logical :: exists
        integer :: unit, ios

        castep_path = trim(dir_path) // '/scan.castep'
        inquire(file=trim(castep_path), exist=exists)
        if (exists) return

        tmp_file = '/tmp/pes3d_find_castep.tmp'
        cmd = 'ls "' // trim(dir_path) // '"/*.castep 2>/dev/null | head -1 > ' // trim(tmp_file)
        call execute_command_line(trim(cmd), exitstat=ios)
        if (ios /= 0) then; castep_path = ''; return; end if
        open(newunit=unit, file=trim(tmp_file), status='old', action='read', iostat=ios)
        if (ios /= 0) then; castep_path = ''; return; end if
        read(unit, '(a)', iostat=ios) castep_path
        close(unit, status='delete')
        if (ios /= 0 .or. len_trim(castep_path) == 0) castep_path = ''
    end function find_castep_in_dir


    function parse_castep_energy(castep_file, ios) result(energy)
        !! Parse "Final energy =  -XXXX.XXXX eV" from a .castep file.
        character(len=*), intent(in) :: castep_file
        integer, intent(out) :: ios
        real(dp) :: energy

        character(len=256) :: line
        integer :: unit, eq_pos

        energy = 0.0_dp
        ios = 0

        open(newunit=unit, file=trim(castep_file), status='old', action='read', iostat=ios)
        if (ios /= 0) return

        do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, 'Final energy') > 0 .or. &
                index(line, 'final energy') > 0 .or. &
                index(line, 'Final Energy') > 0) then
                eq_pos = index(line, '=')
                if (eq_pos > 0) then
                    read(line(eq_pos+1:), *, iostat=ios) energy
                end if
                ios = 0
                close(unit)
                return
            end if
        end do
        ios = 1
        close(unit)
    end function parse_castep_energy


    subroutine rewrite_json_with_energies_3d(json_path, energies, has_energy, &
                                              nx, ny, nz, use_sym, ios)
        !! Rewrite pes3d_metadata.json with collected energy values.
        character(len=*), intent(in) :: json_path
        real(dp), intent(in) :: energies(:)
        logical, intent(in) :: has_energy(:)
        integer, intent(in) :: nx, ny, nz
        logical, intent(in) :: use_sym
        integer, intent(out) :: ios

        character(len=1024) :: tmp_path
        integer :: unit_in, unit_out, j, n_total
        character(len=4096) :: line
        logical :: in_energies, wrote_energy

        ! Drop unused warning
        if (use_sym) continue

        tmp_path = trim(json_path) // '.tmp'

        open(newunit=unit_in, file=trim(json_path), status='old', action='read', iostat=ios)
        if (ios /= 0) return
        open(newunit=unit_out, file=trim(tmp_path), status='replace', action='write', iostat=ios)
        if (ios /= 0) then; close(unit_in); return; end if

        n_total = nx * ny * nz
        in_energies = .false.
        wrote_energy = .false.

        do
            read(unit_in, '(a)', iostat=ios) line
            if (ios /= 0) exit

            if (index(line, '"energies"') > 0 .and. index(line, '"energies_expanded"') == 0) then
                in_energies = .true.
                write(unit_out, '(a)') trim(line)
                cycle
            end if

            if (in_energies .and. .not. wrote_energy) then
                if (index(line, ']') > 0) then
                    do j = 1, n_total
                        if (has_energy(j)) then
                            if (j < n_total) then
                                write(unit_out, '(a, f20.8, a)') '    ', energies(j), ','
                            else
                                write(unit_out, '(a, f20.8)') '    ', energies(j)
                            end if
                        else
                            if (j < n_total) then
                                write(unit_out, '(a)') '    null,'
                            else
                                write(unit_out, '(a)') '    null'
                            end if
                        end if
                    end do
                    write(unit_out, '(a)') '  ],'
                    wrote_energy = .true.
                    in_energies = .false.
                    cycle
                else
                    cycle
                end if
            end if

            if (index(line, '"has_energies"') > 0) then
                write(unit_out, '(a)') '  "has_energies": true'
                cycle
            end if

            write(unit_out, '(a)') trim(line)
        end do

        close(unit_in)
        close(unit_out)
        call execute_command_line('mv "' // trim(tmp_path) // '" "' // trim(json_path) // '"', &
            exitstat=ios)
    end subroutine rewrite_json_with_energies_3d


    subroutine extract_json_int(line, key, val, ios)
        !! Extract integer value for a JSON key from a line.
        character(len=*), intent(in) :: line, key
        integer, intent(out) :: val, ios
        integer :: kp, cp
        kp = index(line, trim(key))
        if (kp == 0) then; ios = 1; return; end if
        cp = index(line(kp:), ':')
        if (cp == 0) then; ios = 1; return; end if
        read(line(kp+cp:), *, iostat=ios) val
    end subroutine extract_json_int


    subroutine extract_json_two_reals(line, r1, r2)
        !! Extract two real values from a JSON array like "fx_range": [0.2, 0.4]
        character(len=*), intent(in) :: line
        real(dp), intent(out) :: r1, r2
        integer :: bp, ep, ios
        character(len=256) :: sub
        bp = index(line, '[')
        ep = index(line, ']')
        r1 = 0.0_dp; r2 = 1.0_dp
        if (bp > 0 .and. ep > bp) then
            sub = line(bp+1:ep-1)
            read(sub, *, iostat=ios) r1, r2
        end if
    end subroutine extract_json_two_reals


    subroutine parse_symop_line(line, rot)
        !! Parse the first row of a sym_op rot matrix from a JSON line.
        !! Example: "rot": [[1,0,0],[0,1,0],[0,0,1]]
        character(len=*), intent(in) :: line
        integer, intent(out) :: rot(3,3)
        integer :: p1, p2, p3, ios
        character(len=256) :: sub
        ! Extract row 1 from the first [...] after "rot"
        p1 = index(line, '[[', .false.) + 1  ! after the first [[
        if (p1 <= 1) p1 = index(line, '[') + 1
        p2 = index(line(p1:), ']') + p1 - 1
        p3 = index(line(p2+1:), '[') + p2
        if (p1 > 1 .and. p2 > p1) then
            sub = line(p1:p2-1)
            read(sub, *, iostat=ios) rot(1,1), rot(1,2), rot(1,3)
        end if
    end subroutine parse_symop_line


    subroutine parse_symop_matrix_row(line, row)
        !! Parse one row [a,b,c] of a 3x3 matrix.
        character(len=*), intent(in) :: line
        integer, intent(out) :: row(3)
        integer :: p1, p2, ios
        p1 = index(line, '[')
        p2 = index(line, ']')
        if (p1 > 0 .and. p2 > p1) then
            read(line(p1+1:p2-1), *, iostat=ios) row(1), row(2), row(3)
        end if
        if (ios /= 0) row = 0
    end subroutine parse_symop_matrix_row


    subroutine parse_symop_trans_from_line(line, trans)
        !! Parse trans vector from the third rot row line.
        character(len=*), intent(in) :: line
        real(dp), intent(out) :: trans(3)
        integer :: p1, p2, ios
        trans = 0.0_dp
        p1 = index(line, '"trans":')
        if (p1 > 0) then
            p1 = index(line(p1:), '[') + p1
            p2 = index(line(p1:), ']') + p1 - 1
            if (p1 > 0 .and. p2 > p1) then
                read(line(p1:p2-1), *, iostat=ios) trans(1), trans(2), trans(3)
            end if
        end if
    end subroutine parse_symop_trans_from_line


    subroutine parse_trans_from_json_line(line, trans)
        !! Parse translation vector from a JSON line containing "trans": [...]
        character(len=*), intent(in) :: line
        real(dp), intent(out) :: trans(3)
        integer :: p1, p2, ios
        trans = 0.0_dp
        p1 = index(line, '"trans":')
        if (p1 > 0) then
            p1 = index(line(p1:), '[') + p1
            p2 = index(line(p1:), ']') + p1 - 1
            if (p1 > 0 .and. p2 > p1) then
                read(line(p1:p2-1), *, iostat=ios) trans(1), trans(2), trans(3)
            end if
        end if
    end subroutine parse_trans_from_json_line


    pure subroutine wrap_to_unit(x)
        !! Wrap fractional coordinate into [0, 1).
        real(dp), intent(inout) :: x
        real(dp), parameter :: ONE = 1.0_dp, ZERO = 0.0_dp
        x = x - floor(x)
        if (x < ZERO) x = x + ONE
        if (x >= ONE .or. x < ZERO) x = x - aint(x)
    end subroutine wrap_to_unit

end module pes3d
