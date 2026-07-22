module pes
    !! Unified Potential Energy Surface scan module — 2D and 3D grid generation,
    !! Gaussian Cube file output, result collection, and symmetry expansion.
    !!
    !! Replaces the formerly split pes_scan.f90 (2D) and pes3d.f90 (3D) modules.
    !! A 2D PES is represented as a 3D cube with nz=1; nz>1 is a full 3D scan.
    !! Line 2 of each cube file carries a compact JSON metadata block that
    !! describes the scan plane, mobile atom, lattice, and symmetry info.
    use castep_config, only: dp, castep_config_t, cif_data_t, atom_t, sym_op_t, &
        IO_WRITE_FAIL, IO_PARSE_ERROR, IO_FILE_NOT_FOUND
    implicit none
    private

    ! ── Public types and constants ──
    public :: pes_grid_t
    public :: compute_local_grid_bounds
    public :: generate_pes_grid_points
    public :: write_pes_cube
    public :: collect_pes_energies
    public :: symmetry_expand_energies

    ! Legacy wrappers — kept for backward compatibility with poscastep_menu.f90
    public :: write_pes_metadata_json   ! deprecated: use write_pes_cube
    public :: write_pes3d_cube          ! deprecated: use write_pes_cube

    ! Maximum grid sizes
    integer, parameter :: MAX_GRID_2D = 200
    integer, parameter :: MAX_GRID_3D = 50

    ! ── Unified grid type ──

    type :: pes_grid_t
        integer  :: ndim = 2                       ! 2 or 3
        integer  :: plane_axis(2) = [1, 2]         ! 2D: scan plane axis indices (1-based)
        integer  :: n_points(3) = [5, 5, 1]        ! Nx, Ny[, Nz] (Nz=1 for 2D)
        real(dp) :: frac_range(3,2) = 0.0_dp       ! [fmin, fmax] per axis
        integer  :: mobile_atom_idx = 1             ! 1-based index of mobile/ref atom
        character(len=8) :: scan_mode = 'SP'        ! 'SP' or 'RELAX'
        ! 3D-specific
        logical  :: use_symmetry = .false.          ! whether space-group symmetry is used
        real(dp) :: ref_frac(3) = 0.0_dp            ! reference atom original fractional coords
        real(dp) :: half_dist(3) = 0.0_dp           ! symmetry half-distance per axis
    end type pes_grid_t

contains

    ! ═══════════════════════════════════════════════════════════════════════════
    !  Grid point generation
    ! ═══════════════════════════════════════════════════════════════════════════

    subroutine generate_pes_grid_points(grid, frac_points, n_total, iostat, iomsg)
        !! Generate fractional coordinate points for the scan grid.
        !! 2D: returns frac_points(n_total, 2) with columns (fx, fy).
        !! 3D: returns frac_points(n_total, 3) with columns (fx, fy, fz).
        !! Order: x fastest, then y, then z (Fortran inner-loop convention).
        type(pes_grid_t), intent(in) :: grid
        real(dp), allocatable, intent(out) :: frac_points(:,:)
        integer, intent(out) :: n_total
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        if (grid%ndim == 2) then
            call generate_2d_points(grid, frac_points, n_total, iostat, iomsg)
        else
            call generate_3d_points(grid, frac_points, n_total, iostat, iomsg)
        end if
    end subroutine generate_pes_grid_points


    subroutine generate_2d_points(grid, frac_points, n_total, iostat, iomsg)
        type(pes_grid_t), intent(in) :: grid
        real(dp), allocatable, intent(out) :: frac_points(:,:)
        integer, intent(out) :: n_total, iostat
        character(len=*), optional, intent(out) :: iomsg

        integer :: i, j, idx, nx, ny
        real(dp) :: dx, dy, fx, fy

        iostat = 0; n_total = 0
        nx = grid%n_points(1); ny = grid%n_points(2)

        if (nx < 2 .or. ny < 2) then
            iostat = 1
            if (present(iomsg)) iomsg = 'Grid must have at least 2 points per dimension'
            return
        end if
        if (nx > MAX_GRID_2D .or. ny > MAX_GRID_2D) then
            iostat = 2
            if (present(iomsg)) iomsg = '2D grid size exceeds maximum'
            return
        end if

        n_total = nx * ny
        allocate(frac_points(n_total, 2), stat=iostat)
        if (iostat /= 0) then
            if (present(iomsg)) iomsg = 'Memory allocation failed for grid points'
            return
        end if

        dx = 0.0_dp; dy = 0.0_dp
        if (nx > 1) dx = (grid%frac_range(1,2) - grid%frac_range(1,1)) / real(nx - 1, dp)
        if (ny > 1) dy = (grid%frac_range(2,2) - grid%frac_range(2,1)) / real(ny - 1, dp)

        idx = 0
        do j = 0, ny - 1
            fy = grid%frac_range(2,1) + j * dy
            do i = 0, nx - 1
                fx = grid%frac_range(1,1) + i * dx
                idx = idx + 1
                frac_points(idx, 1) = fx
                frac_points(idx, 2) = fy
            end do
        end do
    end subroutine generate_2d_points


    subroutine generate_3d_points(grid, frac_points, n_total, iostat, iomsg)
        type(pes_grid_t), intent(in) :: grid
        real(dp), allocatable, intent(out) :: frac_points(:,:)
        integer, intent(out) :: n_total, iostat
        character(len=*), optional, intent(out) :: iomsg

        integer :: i, j, k, idx, nx, ny, nz
        real(dp) :: dx, dy, dz, fx, fy, fz

        iostat = 0; n_total = 0
        nx = grid%n_points(1); ny = grid%n_points(2); nz = grid%n_points(3)

        if (any(grid%n_points < 2)) then
            iostat = 1
            if (present(iomsg)) iomsg = 'Grid must have at least 2 points per dimension'
            return
        end if
        if (any(grid%n_points > MAX_GRID_3D)) then
            iostat = 2
            if (present(iomsg)) iomsg = '3D grid size exceeds maximum'
            return
        end if

        n_total = nx * ny * nz
        allocate(frac_points(n_total, 3), stat=iostat)
        if (iostat /= 0) then
            if (present(iomsg)) iomsg = 'Memory allocation failed for grid points'
            return
        end if

        dx = 0.0_dp; dy = 0.0_dp; dz = 0.0_dp
        if (nx > 1) dx = (grid%frac_range(1,2) - grid%frac_range(1,1)) / real(nx - 1, dp)
        if (ny > 1) dy = (grid%frac_range(2,2) - grid%frac_range(2,1)) / real(ny - 1, dp)
        if (nz > 1) dz = (grid%frac_range(3,2) - grid%frac_range(3,1)) / real(nz - 1, dp)

        idx = 0
        do k = 0, nz - 1
            fz = grid%frac_range(3,1) + k * dz
            do j = 0, ny - 1
                fy = grid%frac_range(2,1) + j * dy
                do i = 0, nx - 1
                    fx = grid%frac_range(1,1) + i * dx
                    idx = idx + 1
                    frac_points(idx, 1) = fx
                    frac_points(idx, 2) = fy
                    frac_points(idx, 3) = fz
                end do
            end do
        end do
    end subroutine generate_3d_points


    ! ═══════════════════════════════════════════════════════════════════════════
    !  Symmetry-aware local grid bounds (3D only)
    ! ═══════════════════════════════════════════════════════════════════════════

    subroutine compute_local_grid_bounds(ref_idx, atoms, n_atoms, sym_ops, n_symops, &
                                          half_dist, iostat, iomsg)
        !! For the selected reference atom, compute the half-distance along each
        !! fractional axis to the nearest space-group-equivalent atom.
        !! Uses minimum-image convention. Clamps to [1/MAX_GRID_3D, 0.3].
        integer, intent(in) :: ref_idx
        type(atom_t), intent(in) :: atoms(:)
        integer, intent(in) :: n_atoms
        type(sym_op_t), intent(in) :: sym_ops(:)
        integer, intent(in) :: n_symops
        real(dp), intent(out) :: half_dist(3)
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        real(dp) :: frac0(3), new_frac(3), delta(3), min_dist(3)
        integer :: iop
        logical :: found_other
        real(dp), parameter :: TOL = 1.0e-5_dp, MAX_HALF = 0.3_dp

        iostat = 0
        half_dist = MAX_HALF

        if (ref_idx < 1 .or. ref_idx > n_atoms) then
            iostat = 1
            if (present(iomsg)) iomsg = 'Reference atom index out of range'
            return
        end if

        frac0 = [atoms(ref_idx)%x, atoms(ref_idx)%y, atoms(ref_idx)%z]
        min_dist = 1.0_dp

        if (n_symops <= 1) return  ! P1: no symmetry → keep default MAX_HALF

        do iop = 1, n_symops
            new_frac(1) = sym_ops(iop)%rot(1,1)*frac0(1) &
                        + sym_ops(iop)%rot(1,2)*frac0(2) &
                        + sym_ops(iop)%rot(1,3)*frac0(3) + sym_ops(iop)%trans(1)
            new_frac(2) = sym_ops(iop)%rot(2,1)*frac0(1) &
                        + sym_ops(iop)%rot(2,2)*frac0(2) &
                        + sym_ops(iop)%rot(2,3)*frac0(3) + sym_ops(iop)%trans(2)
            new_frac(3) = sym_ops(iop)%rot(3,1)*frac0(1) &
                        + sym_ops(iop)%rot(3,2)*frac0(2) &
                        + sym_ops(iop)%rot(3,3)*frac0(3) + sym_ops(iop)%trans(3)

            call wrap_to_unit(new_frac(1))
            call wrap_to_unit(new_frac(2))
            call wrap_to_unit(new_frac(3))

            delta = new_frac - frac0
            delta = delta - anint(delta)   ! minimum image

            found_other = abs(delta(1)) > TOL .or. abs(delta(2)) > TOL .or. abs(delta(3)) > TOL
            if (found_other) then
                if (abs(delta(1)) > TOL .and. abs(delta(1)) < min_dist(1)) min_dist(1) = abs(delta(1))
                if (abs(delta(2)) > TOL .and. abs(delta(2)) < min_dist(2)) min_dist(2) = abs(delta(2))
                if (abs(delta(3)) > TOL .and. abs(delta(3)) < min_dist(3)) min_dist(3) = abs(delta(3))
            end if
        end do

        if (min_dist(1) < 1.0_dp) half_dist(1) = max(1.0_dp/MAX_GRID_3D, min(MAX_HALF, min_dist(1)*0.45_dp))
        if (min_dist(2) < 1.0_dp) half_dist(2) = max(1.0_dp/MAX_GRID_3D, min(MAX_HALF, min_dist(2)*0.45_dp))
        if (min_dist(3) < 1.0_dp) half_dist(3) = max(1.0_dp/MAX_GRID_3D, min(MAX_HALF, min_dist(3)*0.45_dp))
    end subroutine compute_local_grid_bounds


    ! ═══════════════════════════════════════════════════════════════════════════
    !  Gaussian Cube file output (unified 2D / 3D)
    ! ═══════════════════════════════════════════════════════════════════════════

    subroutine write_pes_cube(cube_path, grid, cfg, energies, has_energy, iostat, iomsg)
        !! Write a Gaussian Cube file for 2D or 3D PES scalar field.
        !!
        !! Line 1: description string
        !! Line 2: compact JSON metadata (plane, mobile atom, lattice, ranges)
        !! Lines 3-6: standard cube header (natom, origin, nx/dv_x, ny/dv_y, nz/dv_z)
        !! Atom lines: Z 0.0  cart_x  cart_y  cart_z
        !! Volumetric data: 6 values per line, NaN for missing points
        !!
        !! 2D PES: nz=1, dv_z = full out-of-plane lattice vector.
        character(len=*), intent(in) :: cube_path
        type(pes_grid_t), intent(in) :: grid
        type(castep_config_t), intent(in), target :: cfg
        real(dp), intent(in) :: energies(:)
        logical,  intent(in) :: has_energy(:)
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        real(dp) :: lattice_vecs(3,3), dv(3,3), origin(3), cart(3)
        real(dp) :: a, b, c, alpha, beta, gamma
        real(dp) :: scan_domain_1, scan_domain_2
        integer :: unit, ios, i, nx, ny, nz, n_total, pa0, pa1, pa2
        real(dp) :: e_val
        character(len=512) :: meta_json
        character(len=64) :: desc, plane_label
        integer :: z_num
        character(len=8) :: el
        logical :: is_2d, has_nan_val

        iostat = 0
        is_2d = (grid%ndim == 2)
        nx = grid%n_points(1); ny = grid%n_points(2)
        nz = grid%n_points(3)
        if (is_2d) nz = 1

        a = cfg%cell_length(1); b = cfg%cell_length(2); c = cfg%cell_length(3)
        alpha = cfg%cell_angle(1); beta = cfg%cell_angle(2); gamma = cfg%cell_angle(3)

        ! ── Compute Cartesian lattice vectors (a ‖ x, b in xy, c general) ──
        call compute_lattice_vectors(a, b, c, alpha, beta, gamma, lattice_vecs)

        ! ── Voxel spacing vectors ──
        pa0 = grid%plane_axis(1); pa1 = grid%plane_axis(2)
        pa2 = 6 - pa0 - pa1  ! the third axis (1+2+3=6)

        if (is_2d) then
            scan_domain_1 = grid%frac_range(1,2) - grid%frac_range(1,1)
            scan_domain_2 = grid%frac_range(2,2) - grid%frac_range(2,1)
            dv(:,1) = lattice_vecs(:,pa0) * scan_domain_1 / max(1, nx - 1)
            dv(:,2) = lattice_vecs(:,pa1) * scan_domain_2 / max(1, ny - 1)
            dv(:,3) = lattice_vecs(:,pa2)   ! full out-of-plane vector (nz=1)
            ! origin = first scan point Cartesian coords + out-of-plane pos
            origin(:) = grid%frac_range(1,1) * lattice_vecs(:,pa0) &
                      + grid%frac_range(2,1) * lattice_vecs(:,pa1) &
                      + cfg_atom_frac(grid, cfg, pa2) * lattice_vecs(:,pa2)
        else
            dv(:,1) = lattice_vecs(:,1) / max(1, nx - 1)
            dv(:,2) = lattice_vecs(:,2) / max(1, ny - 1)
            dv(:,3) = lattice_vecs(:,3) / max(1, nz - 1)
            origin = [0.0_dp, 0.0_dp, 0.0_dp]
        end if

        ! ── Description + metadata ──
        plane_label = plane_name(grid)
        if (is_2d) then
            write(desc, '(a,i0,a,i0,a,a)') '2D PES: ', nx, 'x', ny, ' plane=', trim(plane_label)
        else
            write(desc, '(a,i0,a,i0,a,i0)') '3D PES: ', nx, 'x', ny, 'x', nz
        end if

        meta_json = build_metadata_json(grid, cfg)

        ! ── Write file ──
        open(newunit=unit, file=trim(cube_path), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            iostat = IO_WRITE_FAIL
            if (present(iomsg)) iomsg = 'Cannot write: ' // trim(cube_path)
            return
        end if

        write(unit, '(a)') trim(desc)
        write(unit, '(a)') trim(meta_json)
        write(unit, '(i5,3f12.6)') cfg%num_atoms, origin(1), origin(2), origin(3)
        write(unit, '(i5,3f12.6)') nx, dv(1,1), dv(2,1), dv(3,1)
        write(unit, '(i5,3f12.6)') ny, dv(1,2), dv(2,2), dv(3,2)
        write(unit, '(i5,3f12.6)') nz, dv(1,3), dv(2,3), dv(3,3)

        ! Atom lines: Z, charge=0.0, x_cart, y_cart, z_cart
        do i = 1, cfg%num_atoms
            el = trim(cfg%atom_type(i))
            z_num = element_to_z(el)
            cart(:) = cfg%atom_x(i)*lattice_vecs(:,1) &
                    + cfg%atom_y(i)*lattice_vecs(:,2) &
                    + cfg%atom_z(i)*lattice_vecs(:,3)
            write(unit, '(i5,4f12.6)') z_num, 0.0_dp, cart(1), cart(2), cart(3)
        end do

        ! Volumetric data: 6 per line, NaN for missing
        n_total = nx * ny * nz
        has_nan_val = .false.
        do i = 1, n_total
            if (has_energy(i)) then
                e_val = energies(i)
            else
                e_val = ieee_nan()
                has_nan_val = .true.
            end if
            write(unit, '(es13.5)', advance='no') e_val
            if (mod(i, 6) == 0 .or. i == n_total) write(unit, '(a)') ''
        end do

        close(unit)
        write(*, '(a)') '  Cube file written: ' // trim(cube_path)
        if (has_nan_val) write(*, '(a)') '  (NaN entries for failed/missing grid points)'
    end subroutine write_pes_cube


    ! ── Legacy: write_pes_metadata_json (deprecated, kept for backward compat) ──

    subroutine write_pes_metadata_json(json_path, grid, cfg, iostat, iomsg)
        !! Deprecated: use write_pes_cube() instead.
        !! Redirects to cube writing at the same directory with name scan.cube.
        character(len=*), intent(in) :: json_path
        type(pes_grid_t), intent(in) :: grid
        type(castep_config_t), intent(in) :: cfg
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        real(dp), allocatable :: energies(:)
        logical, allocatable :: has_en(:)
        integer :: n_total
        character(len=1024) :: cube_path, dir_path
        integer :: slash_pos

        iostat = 0
        n_total = grid%n_points(1) * grid%n_points(2)
        if (grid%ndim == 3) n_total = n_total * grid%n_points(3)

        allocate(energies(n_total), has_en(n_total), stat=iostat)
        if (iostat /= 0) return
        energies = 0.0_dp
        has_en = .false.

        ! Derive cube path from json path
        slash_pos = index(json_path, '/', back=.true.)
        if (slash_pos > 0) then
            dir_path = json_path(1:slash_pos)
            cube_path = trim(dir_path) // 'scan.cube'
        else
            cube_path = 'scan.cube'
        end if

        call write_pes_cube(cube_path, grid, cfg, energies, has_en, iostat, iomsg)
        deallocate(energies, has_en)
    end subroutine write_pes_metadata_json


    ! ── Legacy: write_pes3d_cube (deprecated, kept for backward compat) ──

    subroutine write_pes3d_cube(cube_path, grid, cfg, cif, energies, n_energies, iostat, iomsg)
        !! Deprecated: use write_pes_cube() instead.
        character(len=*), intent(in) :: cube_path
        type(pes_grid_t), intent(in) :: grid
        type(castep_config_t), intent(in), target :: cfg
        type(cif_data_t), intent(in) :: cif
        real(dp), intent(in) :: energies(:)
        integer, intent(in) :: n_energies
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        integer :: i
        logical, allocatable :: has_en(:)

        ! suppress unused
        if (cif%n_symops >= 0) continue

        allocate(has_en(n_energies), stat=iostat)
        if (iostat /= 0) return
        has_en = .true.

        call write_pes_cube(cube_path, grid, cfg, energies, has_en, iostat, iomsg)
        deallocate(has_en)
    end subroutine write_pes3d_cube


    ! ═══════════════════════════════════════════════════════════════════════════
    !  Result collection (unified 2D / 3D)
    ! ═══════════════════════════════════════════════════════════════════════════

    subroutine collect_pes_energies(scan_dir, iostat, iomsg)
        !! Scan grid_*/ subdirectories, parse .castep files, rewrite the cube
        !! file with collected energies (NaN for missing).
        character(len=*), intent(in) :: scan_dir
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        character(len=1024) :: cube_path, grid_dir, castep_file
        integer :: i, j, k, idx, ios, unit_cube
        integer :: natom, nx, ny, nz, n_total, collected, missing
        real(dp), allocatable :: energies(:)
        logical, allocatable :: has_energy(:)
        logical :: exists, is_2d
        real(dp) :: e_val, e_min, e_max
        ! Buffer for header + atom lines
        character(len=4096), allocatable :: header_buf(:)
        integer :: n_header_lines

        iostat = 0
        collected = 0; missing = 0
        e_min = huge(1.0_dp); e_max = -huge(1.0_dp)

        ! Find cube file
        cube_path = trim(scan_dir) // '/scan.cube'
        inquire(file=trim(cube_path), exist=exists)
        if (.not. exists) then
            ! backward compat: try pes_metadata.json path
            cube_path = trim(scan_dir) // '/pes3d.cube'
            inquire(file=trim(cube_path), exist=exists)
        end if
        if (.not. exists) then
            iostat = IO_FILE_NOT_FOUND
            if (present(iomsg)) iomsg = 'scan.cube not found in: ' // trim(scan_dir)
            return
        end if

        ! Parse cube header for dimensions
        call parse_cube_header(cube_path, natom, nx, ny, nz, ios)
        if (ios /= 0) then
            iostat = IO_PARSE_ERROR
            if (present(iomsg)) iomsg = 'Failed to read cube header from: ' // trim(cube_path)
            return
        end if

        is_2d = (nz == 1)
        n_total = nx * ny * nz
        n_header_lines = 6 + natom  ! lines to preserve (header + atoms)

        ! Read and buffer header + atom lines
        allocate(header_buf(n_header_lines), stat=ios)
        if (ios /= 0) then
            iostat = 1; return
        end if
        open(newunit=unit_cube, file=trim(cube_path), status='old', action='read', iostat=ios)
        do i = 1, n_header_lines
            read(unit_cube, '(a)', iostat=ios) header_buf(i)
            if (ios /= 0) then
                close(unit_cube); deallocate(header_buf); iostat = IO_PARSE_ERROR; return
            end if
        end do
        close(unit_cube)

        ! Allocate energy tracking arrays
        allocate(energies(n_total), has_energy(n_total), stat=ios)
        if (ios /= 0) then
            deallocate(header_buf); iostat = 1; return
        end if
        energies = 0.0_dp
        has_energy = .false.

        ! Scan grid directories
        if (is_2d) then
            do j = 1, ny
                do i = 1, nx
                    write(grid_dir, '(a,a,i3.3,a,i3.3)') trim(scan_dir), '/grid_', i, '_', j
                    castep_file = find_castep_in_dir(grid_dir)
                    if (len_trim(castep_file) == 0) then
                        missing = missing + 1; cycle
                    end if
                    e_val = parse_castep_energy(castep_file, ios)
                    idx = (j-1)*nx + i
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
        else
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        write(grid_dir, '(a,a,i3.3,a,i3.3,a,i3.3)') &
                            trim(scan_dir), '/grid_', i, '_', j, '_', k
                        castep_file = find_castep_in_dir(grid_dir)
                        if (len_trim(castep_file) == 0) then
                            missing = missing + 1; cycle
                        end if
                        e_val = parse_castep_energy(castep_file, ios)
                        idx = (k-1)*nx*ny + (j-1)*nx + i
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
        end if

        ! Summary
        write(*, '(a)') '  ── Collection Summary ──'
        write(*, '(a, i0, a, i0)') '  Collected: ', collected, ' / ', n_total
        if (missing > 0) write(*, '(a, i0)') '  Missing:   ', missing
        if (collected > 0) then
            write(*, '(a, f18.8)') '  E min (eV): ', e_min
            write(*, '(a, f18.8)') '  E max (eV): ', e_max
        end if

        ! Rewrite cube with collected energies
        call rewrite_cube_with_energies(cube_path, header_buf, n_header_lines, &
            energies, has_energy, nx, ny, nz, ios)
        if (ios /= 0) then
            iostat = IO_WRITE_FAIL
            if (present(iomsg)) iomsg = 'Failed to rewrite cube with energies'
        end if

        deallocate(energies, has_energy, header_buf)
    end subroutine collect_pes_energies


    ! ═══════════════════════════════════════════════════════════════════════════
    !  Symmetry expansion (3D only)
    ! ═══════════════════════════════════════════════════════════════════════════

    subroutine symmetry_expand_energies(scan_dir, iostat, iomsg)
        !! Read energies from cube, apply symmetry ops to expand local grid
        !! to full [0,1)^3 cell. Writes pes3d_expanded.cube.
        !!
        !! Currently works with the older JSON metadata path for sym_ops.
        !! The cube line 2 carries n_symops for flagging, but sym_ops matrices
        !! are parsed from pes3d_metadata.json if present.
        character(len=*), intent(in) :: scan_dir
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        integer, parameter :: MAX_SYM_OPS = 256, MAX_EXP = 100

        integer :: nx, ny, nz, n_local, natom
        real(dp) :: fx_range(2), fy_range(2), fz_range(2)
        real(dp), allocatable :: local_energies(:)
        logical, allocatable :: local_has(:)

        integer :: n_symops
        integer :: rot(3,3,MAX_SYM_OPS)
        real(dp) :: trans(3,MAX_SYM_OPS)

        integer :: exp_nx, exp_ny, exp_nz, n_exp
        real(dp), allocatable :: exp_energies(:)
        logical, allocatable :: exp_filled(:)
        real(dp) :: sp_x, sp_y, sp_z
        real(dp) :: local_frac(3), exp_frac(3)
        real(dp) :: e_min_final, e_max_final
        real(dp) :: el

        character(len=1024) :: cube_path, exp_path, json_path, line
        integer :: unit_in, ios, i, j, k, idx, iop, iexp, ei, ej, ek
        logical :: exists, use_sym
        character(len=4096), allocatable :: header_buf(:)
        integer :: n_header_lines

        iostat = 0
        n_symops = 0
        rot = 0; trans = 0.0_dp
        fx_range = [0.0_dp, 1.0_dp]
        fy_range = [0.0_dp, 1.0_dp]
        fz_range = [0.0_dp, 1.0_dp]
        use_sym = .false.

        cube_path = trim(scan_dir) // '/scan.cube'
        inquire(file=trim(cube_path), exist=exists)
        if (.not. exists) then
            cube_path = trim(scan_dir) // '/pes3d.cube'
            inquire(file=trim(cube_path), exist=exists)
        end if
        if (.not. exists) then
            iostat = IO_FILE_NOT_FOUND
            if (present(iomsg)) iomsg = 'Cube file not found in: ' // trim(scan_dir)
            return
        end if

        ! Parse cube header for dims + line 2 JSON for use_symmetry
        call parse_cube_header(cube_path, natom, nx, ny, nz, ios)
        if (ios /= 0 .or. nz == 1) then
            if (present(iomsg)) iomsg = 'Symmetry expansion requires 3D cube (nz>1)'
            if (ios == 0) iostat = 1
            if (ios /= 0) iostat = ios
            return
        end if

        ! Read line 2 to check use_symmetry
        open(newunit=unit_in, file=trim(cube_path), status='old', action='read', iostat=ios)
        read(unit_in, '(a)') line   ! line 1: skip
        read(unit_in, '(a)') line   ! line 2: metadata JSON
        close(unit_in)
        use_sym = index(line, '"use_symmetry":true') > 0
        if (.not. use_sym) then
            if (present(iomsg)) iomsg = 'No symmetry — skipping expansion'
            return
        end if

        ! Parse fractional ranges from JSON
        call extract_json_two_reals_by_key(line, '"fx_range"', fx_range(1), fx_range(2))
        call extract_json_two_reals_by_key(line, '"fy_range"', fy_range(1), fy_range(2))
        call extract_json_two_reals_by_key(line, '"fz_range"', fz_range(1), fz_range(2))

        ! Parse sym_ops from sidecar JSON if available
        json_path = trim(scan_dir) // '/pes3d_metadata.json'
        inquire(file=trim(json_path), exist=exists)
        if (exists) then
            call parse_symops_from_json(json_path, rot, trans, n_symops, MAX_SYM_OPS, ios)
            if (ios /= 0 .or. n_symops < 1) then
                if (present(iomsg)) iomsg = 'Failed to parse sym_ops from JSON'
                iostat = IO_PARSE_ERROR; return
            end if
        else
            if (present(iomsg)) iomsg = 'Symmetry metadata JSON not found — cannot expand'
            iostat = IO_FILE_NOT_FOUND; return
        end if

        ! Read energies from cube
        n_local = nx * ny * nz
        n_header_lines = 6 + natom
        allocate(local_energies(n_local), local_has(n_local), header_buf(n_header_lines), stat=ios)
        if (ios /= 0) then; iostat = 1; return; end if
        local_has = .true.

        open(newunit=unit_in, file=trim(cube_path), status='old', action='read', iostat=ios)
        do i = 1, n_header_lines
            read(unit_in, '(a)') header_buf(i)
        end do
        do i = 1, n_local
            read(unit_in, *, iostat=ios) el
            if (ios /= 0) then
                local_energies(i) = 0.0_dp
                local_has(i) = .false.
                ios = 0
            else
                local_energies(i) = el
                local_has(i) = .not. ieee_is_nan(el)
            end if
        end do
        close(unit_in)

        ! Build expanded grid
        sp_x = (fx_range(2) - fx_range(1)) / max(1, nx - 1)
        sp_y = (fy_range(2) - fy_range(1)) / max(1, ny - 1)
        sp_z = (fz_range(2) - fz_range(1)) / max(1, nz - 1)
        exp_nx = min(MAX_EXP, nint(1.0_dp / sp_x) + 1)
        exp_ny = min(MAX_EXP, nint(1.0_dp / sp_y) + 1)
        exp_nz = min(MAX_EXP, nint(1.0_dp / sp_z) + 1)
        n_exp = exp_nx * exp_ny * exp_nz

        allocate(exp_energies(n_exp), exp_filled(n_exp), stat=ios)
        if (ios /= 0) then
            deallocate(local_energies, local_has, header_buf); iostat = 1; return
        end if
        exp_energies = huge(1.0_dp)
        exp_filled = .false.

        ! Apply FORWARD symmetry ops to each local grid point
        do iop = 1, n_symops
            do k = 0, nz - 1
                local_frac(3) = fz_range(1) + k * sp_z
                do j = 0, ny - 1
                    local_frac(2) = fy_range(1) + j * sp_y
                    do i = 0, nx - 1
                        local_frac(1) = fx_range(1) + i * sp_x
                        exp_frac(1) = rot(1,1,iop)*local_frac(1) + rot(1,2,iop)*local_frac(2) &
                                    + rot(1,3,iop)*local_frac(3) + trans(1,iop)
                        exp_frac(2) = rot(2,1,iop)*local_frac(1) + rot(2,2,iop)*local_frac(2) &
                                    + rot(2,3,iop)*local_frac(3) + trans(2,iop)
                        exp_frac(3) = rot(3,1,iop)*local_frac(1) + rot(3,2,iop)*local_frac(2) &
                                    + rot(3,3,iop)*local_frac(3) + trans(3,iop)
                        call wrap_to_unit(exp_frac(1))
                        call wrap_to_unit(exp_frac(2))
                        call wrap_to_unit(exp_frac(3))
                        ei = nint(exp_frac(1) / sp_x)
                        ej = nint(exp_frac(2) / sp_y)
                        ek = nint(exp_frac(3) / sp_z)
                        ei = modulo(ei, exp_nx)
                        ej = modulo(ej, exp_ny)
                        ek = modulo(ek, exp_nz)
                        idx = k * nx * ny + j * nx + i + 1
                        iexp = ek * exp_nx * exp_ny + ej * exp_nx + ei + 1
                        if (iexp >= 1 .and. iexp <= n_exp) then
                            if (local_has(idx) .and. &
                                (.not. exp_filled(iexp) .or. local_energies(idx) < exp_energies(iexp))) then
                                exp_energies(iexp) = local_energies(idx)
                                exp_filled(iexp) = .true.
                            end if
                        end if
                    end do
                end do
            end do
        end do

        ! Compute final energy range
        e_min_final = huge(1.0_dp); e_max_final = -huge(1.0_dp)
        do iexp = 1, n_exp
            if (exp_filled(iexp)) then
                if (exp_energies(iexp) < e_min_final) e_min_final = exp_energies(iexp)
                if (exp_energies(iexp) > e_max_final) e_max_final = exp_energies(iexp)
            end if
        end do

        ! Write expanded cube (header from original, with updated dims)
        exp_path = trim(scan_dir) // '/pes3d_expanded.cube'
        call write_expanded_cube(exp_path, header_buf, n_header_lines, &
            exp_energies, exp_filled, exp_nx, exp_ny, exp_nz, ios)

        write(*, '(a)') '  ── Symmetry Expansion ──'
        write(*, '(a, i0, a, i0, a, i0, a, i0)') '  Expanded grid: ', &
            exp_nx, ' x ', exp_ny, ' x ', exp_nz, ' = ', n_exp
        if (e_min_final < huge(1.0_dp)) then
            write(*, '(a, f18.8)') '  E min (expanded): ', e_min_final
            write(*, '(a, f18.8)') '  E max (expanded): ', e_max_final
        end if

        deallocate(local_energies, local_has, exp_energies, exp_filled, header_buf)
    end subroutine symmetry_expand_energies


    ! ═══════════════════════════════════════════════════════════════════════════
    !  Private helpers
    ! ═══════════════════════════════════════════════════════════════════════════

    ! ── Lattice vector computation ──

    subroutine compute_lattice_vectors(a, b, c, alpha_deg, beta_deg, gamma_deg, vecs)
        !! Compute Cartesian lattice vectors: a along x, b in xy-plane, c general.
        real(dp), intent(in) :: a, b, c, alpha_deg, beta_deg, gamma_deg
        real(dp), intent(out) :: vecs(3,3)
        real(dp) :: al, be, ga

        al = alpha_deg * 3.141592653589793_dp / 180.0_dp
        be = beta_deg  * 3.141592653589793_dp / 180.0_dp
        ga = gamma_deg * 3.141592653589793_dp / 180.0_dp

        vecs(1,1) = a
        vecs(2,1) = 0.0_dp
        vecs(3,1) = 0.0_dp
        vecs(1,2) = b * cos(ga)
        vecs(2,2) = b * sin(ga)
        vecs(3,2) = 0.0_dp
        vecs(1,3) = c * cos(be)
        vecs(2,3) = c * (cos(al) - cos(be)*cos(ga)) / max(sin(ga), 1.0e-10_dp)
        vecs(3,3) = sqrt(max(c*c - vecs(1,3)**2 - vecs(2,3)**2, 0.0_dp))
    end subroutine compute_lattice_vectors


    ! ── Extract mobile atom's fractional coord on a given axis ──

    function cfg_atom_frac(grid, cfg, axis) result(val)
        type(pes_grid_t), intent(in) :: grid
        type(castep_config_t), intent(in), target :: cfg
        integer, intent(in) :: axis
        real(dp) :: val
        integer :: mi
        mi = grid%mobile_atom_idx
        select case (axis)
        case (1); val = cfg%atom_x(mi)
        case (2); val = cfg%atom_y(mi)
        case (3); val = cfg%atom_z(mi)
        case default; val = 0.0_dp
        end select
    end function cfg_atom_frac


    ! ── Plane name helper ──

    function plane_name(grid) result(name)
        type(pes_grid_t), intent(in) :: grid
        character(len=2) :: name
        if (grid%plane_axis(1) == 1 .and. grid%plane_axis(2) == 2) then
            name = 'xy'
        else if (grid%plane_axis(1) == 1 .and. grid%plane_axis(2) == 3) then
            name = 'xz'
        else
            name = 'yz'
        end if
    end function plane_name


    ! ── Build line-2 metadata JSON ──

    function build_metadata_json(grid, cfg) result(json)
        type(pes_grid_t), intent(in) :: grid
        type(castep_config_t), intent(in), target :: cfg
        character(len=512) :: json
        character(len=128) :: lat_str, plane_l, mob_el
        integer :: mi

        mi = grid%mobile_atom_idx
        mob_el = trim(cfg%atom_type(mi))

        write(lat_str, '(a,f12.6,a,f12.6,a,f12.6,a,f10.4,a,f10.4,a,f10.4)') &
            '"a":', cfg%cell_length(1), ',"b":', cfg%cell_length(2), &
            ',"c":', cfg%cell_length(3), ',"alpha":', cfg%cell_angle(1), &
            ',"beta":', cfg%cell_angle(2), ',"gamma":', cfg%cell_angle(3)

        plane_l = plane_name(grid)

        if (grid%ndim == 2) then
            write(json, '(a,a,a,a,a,i0,a,a,a,a,f12.8,a,f12.8,a,a,f12.8,a,f12.8,a,a,a,a)') &
                '{"type":"pes_2d","plane":"', trim(plane_l), &
                '","scan_mode":"', trim(grid%scan_mode), &
                '","mobile_idx":', mi - 1, &
                ',"mobile_el":"', trim(mob_el), '",', &
                '"fx_range":[', grid%frac_range(1,1), ',', grid%frac_range(1,2), '],', &
                '"fy_range":[', grid%frac_range(2,1), ',', grid%frac_range(2,2), '],', &
                '"lattice":{', trim(lat_str), '}}'
        else
            write(json, '(a,a,a,a,i0,a,a,a,a,f12.8,a,f12.8,a,a,f12.8,a,f12.8,a,a,f12.8,a,f12.8,a,a,a,l1,a)') &
                '{"type":"pes_3d","scan_mode":"', trim(grid%scan_mode), &
                '","mobile_idx":', mi - 1, &
                ',"mobile_el":"', trim(mob_el), '",', &
                '"fx_range":[', grid%frac_range(1,1), ',', grid%frac_range(1,2), '],', &
                '"fy_range":[', grid%frac_range(2,1), ',', grid%frac_range(2,2), '],', &
                '"fz_range":[', grid%frac_range(3,1), ',', grid%frac_range(3,2), '],', &
                '"lattice":{', trim(lat_str), '},', &
                '"use_symmetry":', grid%use_symmetry, '}'
        end if
    end function build_metadata_json


    ! ── IEEE NaN helper ──

    function ieee_nan() result(val)
        real(dp) :: val
        ! Portable way to generate quiet NaN without ieee_arithmetic module
        val = 0.0_dp
        val = val / val   ! 0/0 → NaN on virtually all IEEE 754 systems
    end function ieee_nan

    function ieee_is_nan(x) result(is_nan)
        real(dp), intent(in) :: x
        logical :: is_nan
        is_nan = .not. (x == x)   ! NaN is the only value not equal to itself
    end function ieee_is_nan


    ! ── Cube header parser ──

    subroutine parse_cube_header(cube_path, natom, nx, ny, nz, ios)
        !! Parse natom, nx, ny, nz from a cube file header (lines 3-6).
        character(len=*), intent(in) :: cube_path
        integer, intent(out) :: natom, nx, ny, nz, ios
        integer :: unit
        character(len=256) :: line
        real(dp) :: dummy(3)

        natom = 0; nx = 0; ny = 0; nz = 0; ios = 0

        open(newunit=unit, file=trim(cube_path), status='old', action='read', iostat=ios)
        if (ios /= 0) return
        read(unit, '(a)', iostat=ios) line  ! line 1: skip
        read(unit, '(a)', iostat=ios) line  ! line 2: skip
        read(unit, *, iostat=ios) natom, dummy(1), dummy(2), dummy(3)
        if (ios /= 0) then; close(unit); return; end if
        read(unit, *, iostat=ios) nx, dummy(1), dummy(2), dummy(3)
        read(unit, *, iostat=ios) ny, dummy(1), dummy(2), dummy(3)
        read(unit, *, iostat=ios) nz, dummy(1), dummy(2), dummy(3)
        close(unit)
    end subroutine parse_cube_header


    ! ── Cube rewrite with collected energies ──

    subroutine rewrite_cube_with_energies(cube_path, header_buf, n_header_lines, &
                                           energies, has_energy, nx, ny, nz, ios)
        character(len=*), intent(in) :: cube_path
        character(len=*), intent(in) :: header_buf(:)
        integer, intent(in) :: n_header_lines
        real(dp), intent(in) :: energies(:)
        logical, intent(in) :: has_energy(:)
        integer, intent(in) :: nx, ny, nz
        integer, intent(out) :: ios

        character(len=1024) :: tmp_path
        integer :: unit_out, i, n_total
        real(dp) :: e_val

        tmp_path = trim(cube_path) // '.tmp'
        n_total = nx * ny * nz

        open(newunit=unit_out, file=trim(tmp_path), status='replace', action='write', iostat=ios)
        if (ios /= 0) return

        ! Write header + atom lines unchanged
        do i = 1, n_header_lines
            write(unit_out, '(a)') trim(header_buf(i))
        end do

        ! Write energy data (6 per line)
        do i = 1, n_total
            if (has_energy(i)) then
                e_val = energies(i)
            else
                e_val = ieee_nan()
            end if
            write(unit_out, '(es13.5)', advance='no') e_val
            if (mod(i, 6) == 0 .or. i == n_total) write(unit_out, '(a)') ''
        end do

        close(unit_out)
        call execute_command_line('mv "' // trim(tmp_path) // '" "' // trim(cube_path) // '"', exitstat=ios)
    end subroutine rewrite_cube_with_energies


    ! ── Write expanded cube (symmetry expansion output) ──

    subroutine write_expanded_cube(exp_path, header_buf, n_header_lines, &
                                    energies, has_energy, nx, ny, nz, ios)
        character(len=*), intent(in) :: exp_path
        character(len=*), intent(in) :: header_buf(:)
        integer, intent(in) :: n_header_lines
        real(dp), intent(in) :: energies(:)
        logical, intent(in) :: has_energy(:)
        integer, intent(in) :: nx, ny, nz
        integer, intent(out) :: ios

        integer :: unit_out, i, n_total, natom
        real(dp) :: e_val, origin(3), dv(3,3)
        character(len=4096) :: line

        n_total = nx * ny * nz

        open(newunit=unit_out, file=trim(exp_path), status='replace', action='write', iostat=ios)
        if (ios /= 0) return

        ! Rewrite header with updated dimensions (preserving origin, atoms)
        write(unit_out, '(a,i0,a,i0,a,i0)') '3D PES (expanded): ', nx, 'x', ny, 'x', nz
        write(unit_out, '(a)') trim(header_buf(2))  ! preserve line-2 metadata

        ! Re-parse natom and origin from original header
        read(header_buf(3), *) natom, origin(1), origin(2), origin(3)

        ! Recompute voxel vectors for expanded grid
        ! Read original dv from header_buf(4-6), scale by (old_n-1)/(new_n-1)
        ! ...complex. For now, write the header with placeholder voxel info
        ! and let the metadata JSON carry the lattice.
        read(header_buf(3), *) natom, origin(1), origin(2), origin(3)
        write(unit_out, '(i5,3f12.6)') natom, origin(1), origin(2), origin(3)

        ! Read old voxel vectors and scale
        read(header_buf(4), *) i, dv(1,1), dv(2,1), dv(3,1)
        read(header_buf(5), *) i, dv(1,2), dv(2,2), dv(3,2)
        read(header_buf(6), *) i, dv(1,3), dv(2,3), dv(3,3)
        ! Scale: old_n-1 new voxels cover same domain
        ! dv_new = dv_old * (old_n-1) / (new_n-1)
        if (i > 1) dv(:,1) = dv(:,1) * real(i-1, dp) / real(max(1, nx-1), dp)
        read(header_buf(5), *) i
        if (i > 1) dv(:,2) = dv(:,2) * real(i-1, dp) / real(max(1, ny-1), dp)
        read(header_buf(6), *) i
        if (i > 1) dv(:,3) = dv(:,3) * real(i-1, dp) / real(max(1, nz-1), dp)

        write(unit_out, '(i5,3f12.6)') nx, dv(1,1), dv(2,1), dv(3,1)
        write(unit_out, '(i5,3f12.6)') ny, dv(1,2), dv(2,2), dv(3,2)
        write(unit_out, '(i5,3f12.6)') nz, dv(1,3), dv(2,3), dv(3,3)

        ! Atom lines unchanged (lines 7..6+natom)
        do i = 7, n_header_lines
            write(unit_out, '(a)') trim(header_buf(i))
        end do

        ! Energy data
        do i = 1, n_total
            if (has_energy(i)) then
                e_val = energies(i)
            else
                e_val = ieee_nan()
            end if
            write(unit_out, '(es13.5)', advance='no') e_val
            if (mod(i, 6) == 0 .or. i == n_total) write(unit_out, '(a)') ''
        end do

        close(unit_out)
        write(*, '(a)') '  Expanded cube written: ' // trim(exp_path)
    end subroutine write_expanded_cube


    ! ── Parse sym_ops from JSON (for symmetry_expand) ──

    subroutine parse_symops_from_json(json_path, rot, trans, n_symops, max_ops, ios)
        character(len=*), intent(in) :: json_path
        integer, intent(in) :: max_ops
        integer, intent(out) :: rot(3,3,max_ops), n_symops, ios
        real(dp), intent(out) :: trans(3,max_ops)
        character(len=4096) :: line
        integer :: unit, iop

        n_symops = 0; ios = 0
        rot = 0; trans = 0.0_dp

        open(newunit=unit, file=trim(json_path), status='old', action='read', iostat=ios)
        if (ios /= 0) return

        iop = 0
        do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, '"rot"') > 0) then
                iop = iop + 1
                if (iop > max_ops) exit
                call parse_symop_block(unit, line, rot(:,:,iop), trans(:,iop), ios)
            end if
        end do
        close(unit)
        n_symops = iop
    end subroutine parse_symops_from_json


    subroutine parse_symop_block(unit, first_line, rot_mat, trans_vec, ios)
        integer, intent(in) :: unit
        character(len=*), intent(in) :: first_line
        integer, intent(out) :: rot_mat(3,3), ios
        real(dp), intent(out) :: trans_vec(3)
        character(len=4096) :: line

        ios = 0
        call parse_rot_row(first_line, rot_mat(1,:))
        read(unit, '(a)', iostat=ios) line
        if (ios == 0) call parse_rot_row(line, rot_mat(2,:))
        read(unit, '(a)', iostat=ios) line
        if (ios == 0) call parse_rot_row(line, rot_mat(3,:))
        call parse_trans_from_line(line, trans_vec)
        if (all(abs(trans_vec) < 1.0e-10_dp)) call parse_trans_from_line(first_line, trans_vec)
    end subroutine parse_symop_block


    subroutine parse_rot_row(line, row)
        character(len=*), intent(in) :: line
        integer, intent(out) :: row(3)
        integer :: p1, p2, ios
        row = 0
        p1 = index(line, '[')
        p2 = index(line, ']')
        if (p1 > 0 .and. p2 > p1) read(line(p1+1:p2-1), *, iostat=ios) row(1), row(2), row(3)
    end subroutine parse_rot_row


    subroutine parse_trans_from_line(line, trans)
        character(len=*), intent(in) :: line
        real(dp), intent(out) :: trans(3)
        integer :: p1, p2, ios
        trans = 0.0_dp
        p1 = index(line, '"trans":')
        if (p1 > 0) then
            p1 = index(line(p1:), '[') + p1
            p2 = index(line(p1:), ']') + p1 - 1
            if (p1 > 0 .and. p2 > p1) read(line(p1:p2-1), *, iostat=ios) trans(1), trans(2), trans(3)
        end if
    end subroutine parse_trans_from_line


    ! ── JSON utility helpers ──

    subroutine extract_json_int(line, key, val, ios)
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


    subroutine extract_json_two_reals_by_key(line, key, r1, r2)
        !! Like extract_json_two_reals but searches for a specific key first.
        character(len=*), intent(in) :: line, key
        real(dp), intent(out) :: r1, r2
        integer :: kp
        character(len=512) :: sub
        kp = index(line, trim(key))
        if (kp > 0) then
            sub = line(kp:)
            call extract_json_two_reals(sub, r1, r2)
        else
            r1 = 0.0_dp; r2 = 1.0_dp
        end if
    end subroutine extract_json_two_reals_by_key


    ! ── Element symbol → atomic number ──

    integer function element_to_z(el) result(z)
        character(len=*), intent(in) :: el
        character(len=8) :: uel
        integer :: ic
        uel = adjustl(el)
        do ic = 1, len_trim(uel)
            if (uel(ic:ic) >= 'a' .and. uel(ic:ic) <= 'z') &
                uel(ic:ic) = achar(iachar(uel(ic:ic)) - 32)
        end do
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


    ! ── Directory + Castep helpers ──

    function find_castep_in_dir(dir_path) result(castep_path)
        character(len=*), intent(in) :: dir_path
        character(len=1024) :: castep_path, cmd, tmp_file
        logical :: exists
        integer :: unit, ios

        castep_path = trim(dir_path) // '/scan.castep'
        inquire(file=trim(castep_path), exist=exists)
        if (exists) return

        tmp_file = '/tmp/pes_find_castep.tmp'
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
        character(len=*), intent(in) :: castep_file
        integer, intent(out) :: ios
        real(dp) :: energy
        character(len=256) :: line
        integer :: unit, eq_pos

        energy = 0.0_dp; ios = 0
        open(newunit=unit, file=trim(castep_file), status='old', action='read', iostat=ios)
        if (ios /= 0) return
        do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, 'Final energy') > 0 .or. &
                index(line, 'final energy') > 0 .or. &
                index(line, 'Final Energy') > 0) then
                eq_pos = index(line, '=')
                if (eq_pos > 0) read(line(eq_pos+1:), *, iostat=ios) energy
                ios = 0
                close(unit); return
            end if
        end do
        ios = 1
        close(unit)
    end function parse_castep_energy


    ! ── Fractional coordinate wrapping ──

    pure subroutine wrap_to_unit(x)
        real(dp), intent(inout) :: x
        x = x - floor(x)
        if (x < 0.0_dp) x = x + 1.0_dp
        if (x >= 1.0_dp .or. x < 0.0_dp) x = x - aint(x)
    end subroutine wrap_to_unit

end module pes
