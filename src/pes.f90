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
    public :: get_irreducible_grid
    public :: generate_pes_grid_points
    public :: write_pes_cube
    public :: collect_pes_energies
    public :: symops_translation_lcm

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
        integer  :: n_irred = 0                     ! number of irreducible grid points
        real(dp), allocatable :: irred_coords(:,:)   ! (3, n_irred) irreducible coords
        ! Per-point flag: 0 = normal (symmetry-expanded), 1 = special (mobile
        ! atom's own site — fill only its own grid point, do NOT expand; its
        ! energy is the equilibrium value, different from the rest of its
        ! symmetry orbit), 2 = rejected (too close to a fixed atom, no energy).
        integer, allocatable :: irred_flags(:)
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

    subroutine get_irreducible_grid(Na, Nb, Nc, sym_ops, n_symops, &
                                     n_irred, irred_coords, iostat)
        !! Grid symmetry reduction via orbit mapping.
        !! Lays a uniform Na×Nb×Nc grid over the full cell [0,1)^3,
        !! then uses symmetry operations to group equivalent points
        !! into orbits.  Each orbit contributes exactly one irreducible
        !! representative — the minimal set of CASTEP scan points.
        !!
        !! IMPORTANT: Axis-swapping symmetry operations require Na == Nb == Nc.
        !! The caller must enforce cubic grids when the space group has such operations.
        integer, intent(in) :: Na, Nb, Nc, n_symops
        type(sym_op_t), intent(in) :: sym_ops(:)
        integer, intent(out) :: n_irred
        real(dp), allocatable, intent(out) :: irred_coords(:,:)
        integer, intent(out) :: iostat

        logical, allocatable :: visited(:,:,:)
        real(dp) :: x(3), xp(3)
        integer :: i, j, k, s, ip, jp, kp, max_pts

        iostat = 0; n_irred = 0
        max_pts = Na * Nb * Nc

        ! Defensive check: axis-swapping operations require cubic grids
        if (Na /= Nb .or. Nb /= Nc) then
            write(*, '(a)') '  ERROR: orbit mapping requires Na == Nb == Nc for axis-swapping operations.'
            iostat = IO_PARSE_ERROR
            return
        end if

        allocate(visited(0:Na-1, 0:Nb-1, 0:Nc-1), stat=iostat)
        if (iostat /= 0) return
        visited = .false.

        allocate(irred_coords(3, max_pts), stat=iostat)
        if (iostat /= 0) then
            deallocate(visited); return
        end if

        do k = 0, Nc - 1
            do j = 0, Nb - 1
                do i = 0, Na - 1
                    if (visited(i, j, k)) cycle

                    ! Record irreducible representative
                    n_irred = n_irred + 1
                    irred_coords(1, n_irred) = real(i, dp) / real(Na, dp)
                    irred_coords(2, n_irred) = real(j, dp) / real(Nb, dp)
                    irred_coords(3, n_irred) = real(k, dp) / real(Nc, dp)
                    x = irred_coords(:, n_irred)

                    ! Mark entire symmetry orbit as visited
                    do s = 1, n_symops
                        xp(1) = sym_ops(s)%rot(1,1)*x(1) &
                              + sym_ops(s)%rot(1,2)*x(2) &
                              + sym_ops(s)%rot(1,3)*x(3) + sym_ops(s)%trans(1)
                        xp(2) = sym_ops(s)%rot(2,1)*x(1) &
                              + sym_ops(s)%rot(2,2)*x(2) &
                              + sym_ops(s)%rot(2,3)*x(3) + sym_ops(s)%trans(2)
                        xp(3) = sym_ops(s)%rot(3,1)*x(1) &
                              + sym_ops(s)%rot(3,2)*x(2) &
                              + sym_ops(s)%rot(3,3)*x(3) + sym_ops(s)%trans(3)

                        ! Wrap to [0, 1) — use floor(x+eps) for boundary safety:
                        ! when xp ≈ 0.9999999999999, floor(xp)=0 (wrong), floor(xp+eps)=1 (correct)
                        xp(1) = xp(1) - floor(xp(1) + 1.0d-5)
                        xp(2) = xp(2) - floor(xp(2) + 1.0d-5)
                        xp(3) = xp(3) - floor(xp(3) + 1.0d-5)

                        ! Map to nearest grid index with periodic wrap
                        ip = nint(xp(1) * Na)
                        jp = nint(xp(2) * Nb)
                        kp = nint(xp(3) * Nc)
                        ip = modulo(ip, Na)
                        jp = modulo(jp, Nb)
                        kp = modulo(kp, Nc)

                        visited(ip, jp, kp) = .true.
                    end do
                end do
            end do
        end do

        deallocate(visited)
        write(*, '(a, i0, a, i0, a)') '  ── Orbit Mapping ──'
        write(*, '(a, i0, a, i0)') '  Full grid: ', max_pts, &
            '  →  Irreducible: ', n_irred
    end subroutine get_irreducible_grid


    ! ═══════════════════════════════════════════════════════════════════════════
    !  Symmetry translation commensurability check
    ! ═══════════════════════════════════════════════════════════════════════════

    integer function symops_translation_lcm(sym_ops, n_symops) result(m)
        !! Computes the smallest integer m such that every symmetry translation
        !! component t satisfies t*m == integer (within 1e-6 tolerance).
        !!
        !! A grid with N divisible by m is commensurate with the space group:
        !! every symmetry image of a grid point lands exactly on another grid point.
        !!
        !! Returns 0 if any translation component cannot be represented by a
        !! crystallographic denominator (1, 2, 3, 4, 6), indicating the space
        !! group is incompatible with a uniform grid.
        type(sym_op_t), intent(in) :: sym_ops(:)
        integer, intent(in) :: n_symops
        integer :: s, d, denom
        real(dp) :: t
        integer, parameter :: CAND(5) = [1, 2, 3, 4, 6]  ! crystallographic denominators
        logical :: found

        m = 1
        do s = 1, n_symops
            do d = 1, 3
                t = sym_ops(s)%trans(d)
                ! Find smallest candidate denominator that makes t*denom ~ integer
                found = .false.
                do denom = 1, size(CAND)
                    if (abs(t * CAND(denom) - nint(t * CAND(denom))) < 1.0d-6) then
                        m = lcm(m, CAND(denom))
                        found = .true.
                        exit
                    end if
                end do
                if (.not. found) then
                    m = 0  ! non-crystallographic translation
                    return
                end if
            end do
        end do
    end function symops_translation_lcm


    integer function lcm(a, b) result(res)
        !! Least common multiple of two positive integers.
        integer, intent(in) :: a, b
        res = (a * b) / gcd(a, b)
    end function lcm


    integer function gcd(a, b) result(res)
        !! Greatest common divisor via Euclidean algorithm.
        integer, intent(in) :: a, b
        integer :: x, y, r
        x = a; y = b
        do while (y /= 0)
            r = mod(x, y)
            x = y
            y = r
        end do
        res = x
    end function gcd


    ! ═══════════════════════════════════════════════════════════════════════════
    !  Gaussian Cube file output (unified 2D / 3D)
    ! ═══════════════════════════════════════════════════════════════════════════

    subroutine write_pes_cube(cube_path, grid, cfg, energies, has_energy, iostat, iomsg, &
                              sym_ops, n_symops_stored)
        !! Write a Gaussian Cube file for 2D or 3D PES scalar field.
        !!
        !! Line 1: description string
        !! Line 2: compact JSON metadata (plane, mobile atom, lattice, ranges, sym_ops)
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
        type(sym_op_t), intent(in), optional :: sym_ops(:)
        integer, intent(in), optional :: n_symops_stored
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        real(dp) :: lattice_vecs(3,3), dv(3,3), origin(3), cart(3)
        real(dp) :: a, b, c, alpha, beta, gamma
        real(dp) :: scan_domain_1, scan_domain_2
        integer :: unit, ios, i, nx, ny, nz, n_total, pa0, pa1, pa2
        real(dp) :: e_val
        character(len=16384) :: meta_json
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

        meta_json = build_metadata_json(grid, cfg, sym_ops, n_symops_stored)

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
        character(len=16384), allocatable :: header_buf(:)
        integer :: n_header_lines

        ! Symmetry expansion locals
        integer, parameter :: MAX_SYM_OPS = 256
        integer :: n_irr, unit_irr, ii
        logical :: was_expanded, coords_are_exact
        real(dp), allocatable :: irred_energies(:), irred_coords(:,:)
        integer, allocatable :: irred_idx(:,:)
        integer, allocatable :: irred_flags(:)   ! 0=normal, 1=mobile-site special, 2=rejected
        integer :: na_f, nb_f, nc_f
        integer :: n_symops, rot(3,3,MAX_SYM_OPS)
        real(dp) :: trans(3,MAX_SYM_OPS), fx, fy, fz
        integer :: exp_nx, exp_ny, exp_nz, n_exp, iexp, ei, ej, ek, n_filled, n_holes
        integer :: n_holes_extra
        real(dp), allocatable :: exp_energies(:)
        logical, allocatable :: exp_filled(:)
        real(dp) :: exp_frac(3)
        character(len=16384) :: json_line
        character(len=128) :: first_line
        integer :: ix, iy, iz, isrc, idst, old_nx, old_ny, old_nz, ti, tj, tk, n_face_bad
        real(dp) :: dv(3,3)

        iostat = 0
        collected = 0; missing = 0
        e_min = huge(1.0_dp); e_max = -huge(1.0_dp)

        ! Find cube file
        cube_path = trim(scan_dir) // '/scan.cube'
        inquire(file=trim(cube_path), exist=exists)
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
            ! Detect format: irred_coords.dat (orbit mapping) or grid_III_JJJ_KKK
            inquire(file=trim(scan_dir)//'/irred_coords.dat', exist=exists)
            if (exists) then
                ! ═══════════════════════════════════════════════════════════
                !  Symmetry mode: collect → expand in memory → full-cell cube
                ! ═══════════════════════════════════════════════════════════
                ! Step 1: Read irreducible coordinates (dual format support)

                open(newunit=unit_irr, file=trim(scan_dir)//'/irred_coords.dat', &
                     status='old', action='read', iostat=ios)
                if (ios /= 0) then
                    deallocate(energies, has_energy, header_buf)
                    iostat = IO_PARSE_ERROR; return
                end if

                read(unit_irr, '(a)') first_line
                if (index(first_line, '# irred_index_v3') > 0) then
                    ! v3: integer indices + per-point flag (ix iy iz flag)
                    !   flag 0 = normal, 1 = mobile-site special (own grid
                    !   point only), 2 = rejected (overlap, no energy)
                    read(unit_irr, *) n_irr, na_f, nb_f, nc_f
                    allocate(irred_idx(3, n_irr), irred_energies(n_irr), &
                             irred_flags(n_irr), stat=ios)
                    if (ios /= 0) then
                        close(unit_irr); deallocate(energies, has_energy, header_buf)
                        iostat = 1; return
                    end if
                    irred_flags = 0
                    do ii = 1, n_irr
                        read(unit_irr, *, iostat=ios) irred_idx(:, ii), irred_flags(ii)
                        if (ios /= 0) exit
                    end do
                    coords_are_exact = .true.
                else if (index(first_line, '# irred_index_v2') > 0) then
                    ! v2: integer indices only (all points normal)
                    read(unit_irr, *) n_irr, na_f, nb_f, nc_f
                    allocate(irred_idx(3, n_irr), irred_energies(n_irr), &
                             irred_flags(n_irr), stat=ios)
                    if (ios /= 0) then
                        close(unit_irr); deallocate(energies, has_energy, header_buf)
                        iostat = 1; return
                    end if
                    irred_flags = 0
                    do ii = 1, n_irr
                        read(unit_irr, *, iostat=ios) irred_idx(:, ii)
                        if (ios /= 0) exit
                    end do
                    coords_are_exact = .true.
                else
                    ! Legacy format: floating-point fractional coordinates
                    read(first_line, *) n_irr
                    allocate(irred_coords(3, n_irr), irred_energies(n_irr), stat=ios)
                    if (ios /= 0) then
                        close(unit_irr); deallocate(energies, has_energy, header_buf)
                        iostat = 1; return
                    end if
                    do ii = 1, n_irr
                        read(unit_irr, *, iostat=ios) irred_coords(:, ii)
                        if (ios /= 0) exit
                    end do
                    coords_are_exact = .false.
                    write(*, '(a)') '  NOTE: legacy irred_coords.dat — regenerate scan for exact expansion'
                end if
                close(unit_irr)

                ! Step 2: Collect CASTEP energies per irreducible point
                collected = 0; missing = 0
                do ii = 1, n_irr
                    ! Rejected points (flag=2) are NEVER collected — even if a
                    ! stale directory exists (e.g. copied from an old scan),
                    ! their energies are overlap artifacts and would corrupt
                    ! the E_min reference (relative-energy zero).
                    if (allocated(irred_flags)) then
                        if (irred_flags(ii) == 2) cycle
                    end if
                    write(grid_dir, '(a,a,i5.5)') trim(scan_dir), '/irred_', ii
                    castep_file = find_castep_in_dir(grid_dir)
                    if (len_trim(castep_file) > 0) then
                        e_val = parse_castep_energy(castep_file, ios)
                        if (ios == 0) then
                            irred_energies(ii) = e_val
                            collected = collected + 1
                            if (e_val < e_min) e_min = e_val
                            if (e_val > e_max) e_max = e_val
                        else
                            irred_energies(ii) = huge(1.0_dp)
                            missing = missing + 1
                        end if
                    else
                        irred_energies(ii) = huge(1.0_dp)
                        missing = missing + 1
                    end if
                end do

                ! Step 3: Parse sym_ops from cube line-2 JSON
                json_line = header_buf(2)
                if (index(json_line, '"expanded":true') > 0) then
                    write(*, '(a)') '  ── Already expanded — skipping'
                    deallocate(energies, has_energy, header_buf, irred_energies)
                    if (allocated(irred_coords)) deallocate(irred_coords)
                    if (allocated(irred_idx)) deallocate(irred_idx)
                    iostat = 0
                    return  ! Direct return to protect already-expanded data
                end if
                n_symops = 0; rot = 0; trans = 0.0_dp
                call parse_sym_ops_from_json_str(json_line, rot, trans, n_symops, MAX_SYM_OPS, ios)
                if (ios /= 0 .or. n_symops < 1) then
                    write(*, '(a)') '  WARNING: Cannot parse sym_ops from cube JSON — writing partial energies only'
                else
                    ! Step 4: Forward orbit expansion
                    exp_nx = nx; exp_ny = ny; exp_nz = nz
                    n_exp = (exp_nx + 1) * (exp_ny + 1) * (exp_nz + 1)
                    allocate(exp_energies(n_exp), exp_filled(n_exp), stat=ios)
                    if (ios /= 0) then
                        deallocate(energies, has_energy, irred_energies, header_buf)
                        if (allocated(irred_coords)) deallocate(irred_coords)
                        if (allocated(irred_idx)) deallocate(irred_idx)
                        iostat = 1; return
                    end if
                    exp_energies = huge(1.0_dp)
                    exp_filled = .false.

                    ! Forward orbit expansion
                    if (coords_are_exact) then
                        ! ── Pass 1: normal points (flag=0) expanded by all
                        ! symmetry ops. Rejected points (flag=2) are skipped
                        ! unconditionally — their grid points stay NaN (holes)
                        ! even if a stale directory happens to exist.
                        do ii = 1, n_irr
                            if (allocated(irred_flags)) then
                                if (irred_flags(ii) == 2) cycle  ! rejected
                                if (irred_flags(ii) == 1) cycle  ! special — pass 2
                            end if
                            if (irred_energies(ii) > huge(1.0_dp) * 0.5_dp) cycle  ! missing
                            do i = 1, n_symops
                                ! Integer rotation + translation (trans*N is guaranteed integer)
                                ti = nint(trans(1,i) * exp_nx)
                                tj = nint(trans(2,i) * exp_ny)
                                tk = nint(trans(3,i) * exp_nz)
                                ei = rot(1,1,i)*irred_idx(1,ii) + rot(1,2,i)*irred_idx(2,ii) &
                                   + rot(1,3,i)*irred_idx(3,ii) + ti
                                ej = rot(2,1,i)*irred_idx(1,ii) + rot(2,2,i)*irred_idx(2,ii) &
                                   + rot(2,3,i)*irred_idx(3,ii) + tj
                                ek = rot(3,1,i)*irred_idx(1,ii) + rot(3,2,i)*irred_idx(2,ii) &
                                   + rot(3,3,i)*irred_idx(3,ii) + tk
                                ei = modulo(ei, exp_nx)
                                ej = modulo(ej, exp_ny)
                                ek = modulo(ek, exp_nz)
                                iexp = ek * (exp_nx+1) * (exp_ny+1) + ej * (exp_nx+1) + ei + 1
                                if (.not. exp_filled(iexp) .or. irred_energies(ii) < exp_energies(iexp)) then
                                    exp_energies(iexp) = irred_energies(ii)
                                    exp_filled(iexp) = .true.
                                end if
                            end do
                        end do

                        ! ── Pass 2: special points (flag=1) — the mobile atom's
                        ! own site. Its equilibrium energy is NOT symmetry-
                        ! equivalent to the other sites of its orbit (those are
                        ! doubly-occupied → high energy), so it is written ONLY
                        ! at its own grid point, overriding whatever the orbit
                        ! expansion wrote there.
                        if (allocated(irred_flags)) then
                            do ii = 1, n_irr
                                if (irred_flags(ii) /= 1) cycle
                                if (irred_energies(ii) > huge(1.0_dp) * 0.5_dp) cycle
                                ei = irred_idx(1, ii)
                                ej = irred_idx(2, ii)
                                ek = irred_idx(3, ii)
                                iexp = ek * (exp_nx+1) * (exp_ny+1) + ej * (exp_nx+1) + ei + 1
                                if (iexp >= 1 .and. iexp <= n_exp) then
                                    exp_energies(iexp) = irred_energies(ii)  ! override
                                    exp_filled(iexp) = .true.
                                end if
                            end do
                        end if
                    else
                        ! Legacy format: floating-point with rounding errors
                        do ii = 1, n_irr
                            if (irred_energies(ii) > huge(1.0_dp) * 0.5_dp) cycle  ! missing
                            fx = irred_coords(1, ii)
                            fy = irred_coords(2, ii)
                            fz = irred_coords(3, ii)
                            do i = 1, n_symops
                                exp_frac(1) = rot(1,1,i)*fx + rot(1,2,i)*fy + rot(1,3,i)*fz + trans(1,i)
                                exp_frac(2) = rot(2,1,i)*fx + rot(2,2,i)*fy + rot(2,3,i)*fz + trans(2,i)
                                exp_frac(3) = rot(3,1,i)*fx + rot(3,2,i)*fy + rot(3,3,i)*fz + trans(3,i)
                                call wrap_to_unit(exp_frac(1))
                                call wrap_to_unit(exp_frac(2))
                                call wrap_to_unit(exp_frac(3))
                                ei = nint(exp_frac(1) * exp_nx)
                                ej = nint(exp_frac(2) * exp_ny)
                                ek = nint(exp_frac(3) * exp_nz)
                                ei = modulo(ei, exp_nx)
                                ej = modulo(ej, exp_ny)
                                ek = modulo(ek, exp_nz)
                                iexp = ek * (exp_nx+1) * (exp_ny+1) + ej * (exp_nx+1) + ei + 1
                                if (iexp >= 1 .and. iexp <= n_exp) then
                                    if (.not. exp_filled(iexp) .or. irred_energies(ii) < exp_energies(iexp)) then
                                        exp_energies(iexp) = irred_energies(ii)
                                        exp_filled(iexp) = .true.
                                    end if
                                end if
                            end do
                        end do
                    end if

                    ! Periodic boundary: explicit copy layer N from layer 0
                    do iz = 0, exp_nz
                        do iy = 0, exp_ny
                            do ix = 0, exp_nx
                                if (ix < exp_nx .and. iy < exp_ny .and. iz < exp_nz) cycle
                                isrc = modulo(iz, exp_nz) * (exp_nx+1) * (exp_ny+1) &
                                     + modulo(iy, exp_ny) * (exp_nx+1) + modulo(ix, exp_nx) + 1
                                idst = iz * (exp_nx+1) * (exp_ny+1) + iy * (exp_nx+1) + ix + 1
                                exp_energies(idst) = exp_energies(isrc)
                                exp_filled(idst) = exp_filled(isrc)
                            end do
                        end do
                    end do

                    ! Count fill and check for holes
                    n_filled = 0; n_holes = 0
                    do iexp = 1, n_exp
                        if (exp_energies(iexp) < huge(1.0_dp) * 0.5_dp) then
                            n_filled = n_filled + 1
                        else
                            n_holes = n_holes + 1
                        end if
                    end do

                    write(*, '(a)') '  ── Symmetry Expansion (in-line) ──'
                    write(*, '(a, i0, a, i0)') '  Irreducible: ', collected, ' / ', n_irr
                    write(*, '(a, i0, a, i0, a, i0, a, i0)') '  Orbit grid: ', &
                        exp_nx, 'x', exp_ny, 'x', exp_nz, ' = ', exp_nx*exp_ny*exp_nz
                    write(*, '(a, i0, a, i0, a, i0, a, i0)') '  N+1 cube:   ', &
                        exp_nx+1, 'x', exp_ny+1, 'x', exp_nz+1, ' = ', n_exp
                    write(*, '(a, i0, a, i0)') '  Filled: ', n_filled, ' / ', n_exp

                    if (n_holes > 0) then
                        write(*, '(a)') '  ── Expansion Diagnostics ──'
                        write(*, '(a,i0,a,i0)') '  INFO: ', n_holes, &
                            ' grid points left empty (NaN) out of ', n_exp
                        write(*, '(a)') '  Expected: points rejected by the atom-overlap filter ' // &
                            '(mobile atom closer than 1.0 Å to a fixed atom) have no energy.'
                        write(*, '(a)') '  If holes are far from atoms, CASTEP results may be missing.'
                    end if

                    ! Periodic boundary self-check (x-face only as representative)
                    n_face_bad = 0
                    do iz = 0, exp_nz - 1
                        do iy = 0, exp_ny - 1
                            isrc = iz * (exp_nx+1) * (exp_ny+1) + iy * (exp_nx+1) + 1
                            idst = isrc + exp_nx
                            if (abs(exp_energies(idst) - exp_energies(isrc)) > 1.0d-9) then
                                n_face_bad = n_face_bad + 1
                            end if
                        end do
                    end do
                    write(*, '(a,i0)') '  Periodic boundary check (x-face mismatches): ', n_face_bad

                    ! ── Per-grid-point distance filter ──
                    ! Orbit representatives may be safe while some of their
                    ! symmetry images overlap fixed atoms (relocated reps
                    ! expand into overlap points). NaN those grid points.
                    block
                        real(dp) :: atom_cart(3, natom), frac_p(3), cart_p(3), dfc(3)
                        real(dp) :: lvec(3,3), dminp, dtmp
                        integer :: mobile_idx0, kp, ix, iy, iz, ia, iidx
                        real(dp), parameter :: MIND = 1.0_dp   ! Å
                        integer :: zz
                        real(dp) :: chg

                        ! Lattice vectors: dv (N+1 header, = lattice/N) * N
                        lvec(:,1) = dv(:,1) * real(exp_nx, dp)
                        lvec(:,2) = dv(:,2) * real(exp_ny, dp)
                        lvec(:,3) = dv(:,3) * real(exp_nz, dp)

                        ! Atoms (header lines 7..6+natom: Z charge x y z)
                        do ia = 1, natom
                            read(header_buf(6 + ia), *) zz, chg, &
                                atom_cart(1, ia), atom_cart(2, ia), atom_cart(3, ia)
                        end do

                        ! Mobile atom index from JSON ("mobile_idx":N, 0-based)
                        mobile_idx0 = 0
                        kp = index(json_line, '"mobile_idx"')
                        if (kp > 0) then
                            kp = index(json_line(kp:), ':') + kp - 1
                            read(json_line(kp+1:), *, iostat=ios) mobile_idx0
                            mobile_idx0 = mobile_idx0 + 1
                        end if

                        n_holes_extra = 0
                        do iz = 0, exp_nz
                            do iy = 0, exp_ny
                                do ix = 0, exp_nx
                                    iidx = iz * (exp_nx+1) * (exp_ny+1) + iy * (exp_nx+1) + ix + 1
                                    if (.not. exp_filled(iidx)) cycle
                                    frac_p(1) = real(ix, dp) / real(exp_nx, dp)
                                    frac_p(2) = real(iy, dp) / real(exp_ny, dp)
                                    frac_p(3) = real(iz, dp) / real(exp_nz, dp)
                                    cart_p = frac_p(1)*lvec(:,1) + frac_p(2)*lvec(:,2) &
                                           + frac_p(3)*lvec(:,3)
                                    dminp = huge(1.0_dp)
                                    do ia = 1, natom
                                        if (ia == mobile_idx0) cycle
                                        dfc = cart_p - atom_cart(:, ia)
                                        dtmp = sqrt(dfc(1)**2 + dfc(2)**2 + dfc(3)**2)
                                        dminp = min(dminp, dtmp)
                                    end do
                                    if (dminp < MIND) then
                                        exp_filled(iidx) = .false.
                                        n_holes_extra = n_holes_extra + 1
                                    end if
                                end do
                            end do
                        end do
                        if (n_holes_extra > 0) then
                            write(*, '(a, i0, a)') '  Per-point overlap filter: ', n_holes_extra, &
                                ' grid points NaN (inside fixed-atom exclusion radius)'
                        end if
                    end block

                    ! Replace energies/has_energy with expanded arrays
                    deallocate(energies, has_energy)
                    allocate(energies(n_exp), has_energy(n_exp), stat=ios)
                    if (ios /= 0) then
                        deallocate(exp_energies, exp_filled, irred_energies, header_buf)
                        if (allocated(irred_coords)) deallocate(irred_coords)
                        if (allocated(irred_idx)) deallocate(irred_idx)
                        iostat = 1; return
                    end if
                    ! Set has_energy from actual values
                    do iexp = 1, n_exp
                        has_energy(iexp) = (exp_energies(iexp) < huge(1.0_dp) * 0.5_dp)
                    end do
                    energies = exp_energies
                    deallocate(exp_energies, exp_filled)

                    ! Update header for N+1 format
                    n_total = n_exp
                    ! Patch header_buf line 4-6: N+1, dv scaled to lattice/N
                    read(header_buf(4), *) old_nx, dv(1,1), dv(2,1), dv(3,1)
                    read(header_buf(5), *) old_ny, dv(1,2), dv(2,2), dv(3,2)
                    read(header_buf(6), *) old_nz, dv(1,3), dv(2,3), dv(3,3)
                    if (old_nx > 1) dv(:,1) = dv(:,1) * real(old_nx-1, dp) / real(exp_nx, dp)
                    if (old_ny > 1) dv(:,2) = dv(:,2) * real(old_ny-1, dp) / real(exp_ny, dp)
                    if (old_nz > 1) dv(:,3) = dv(:,3) * real(old_nz-1, dp) / real(exp_nz, dp)
                    write(header_buf(4), '(i5,3f12.6)') exp_nx+1, dv(1,1), dv(2,1), dv(3,1)
                    write(header_buf(5), '(i5,3f12.6)') exp_ny+1, dv(1,2), dv(2,2), dv(3,2)
                    write(header_buf(6), '(i5,3f12.6)') exp_nz+1, dv(1,3), dv(2,3), dv(3,3)
                    nx = exp_nx + 1; ny = exp_ny + 1; nz = exp_nz + 1

                    ! Patch line 1 description
                    write(header_buf(1), '(a,i0,a,i0,a,i0)') '3D PES (expanded): ', exp_nx, 'x', exp_ny, 'x', exp_nz
                    ! Patch line 2 JSON: set frac ranges to [0,1], use_symmetry:false, mark expanded
                    call patch_json_full_cell(header_buf(2), json_line)
                    ! Append expanded flag before the closing }
                    i = len_trim(json_line)
                    if (json_line(i:i) == '}') then
                        json_line(i:i) = ','
                        json_line = trim(json_line) // '"expanded":true}'
                    end if
                    header_buf(2) = trim(json_line)
                    was_expanded = .true.
                end if  ! sym_ops parsed
                deallocate(irred_energies)
                if (allocated(irred_coords)) deallocate(irred_coords)
                if (allocated(irred_idx)) deallocate(irred_idx)
                if (allocated(irred_flags)) deallocate(irred_flags)
            else
                ! Rectangular grid format
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
        end if

        ! Summary (skip if expansion already printed its own summary)
        if (.not. was_expanded) then
            write(*, '(a)') '  ── Collection Summary ──'
            write(*, '(a, i0, a, i0)') '  Collected: ', collected, ' / ', n_total
            if (missing > 0) write(*, '(a, i0)') '  Missing:   ', missing
            if (collected > 0) then
                write(*, '(a, f18.8)') '  E min (eV): ', e_min
                write(*, '(a, f18.8)') '  E max (eV): ', e_max
            end if
        end if

        ! Convert to relative energies (E - E_min) for better isosurface visualization
        if (collected > 0 .and. e_min < huge(1.0_dp)) then
            do i = 1, n_total
                if (has_energy(i)) then
                    energies(i) = energies(i) - e_min
                end if
            end do
            write(*, '(a)') '  ── Energy converted to relative (E - E_min) ──'
            write(*, '(a, f18.8)') '  Reference (E_min): ', e_min
            write(*, '(a, f18.8, a)') '  New range: 0.00 to ', e_max - e_min, ' eV'
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
    ! ── Parse sym_ops from JSON string (for inline expansion) ──

    subroutine parse_sym_ops_from_json_str(json, rot, trans, n_symops, max_ops, ios)
        !! Parse sym_ops from cube line-2 JSON: {"sym_ops":[{"rot":[...],"trans":[...]},...]}
        character(len=*), intent(in) :: json
        integer, intent(in) :: max_ops
        integer, intent(out) :: rot(3,3,max_ops), n_symops, ios
        real(dp), intent(out) :: trans(3,max_ops)
        integer :: p, q, p1, p2, iop

        n_symops = 0; ios = 0
        rot = 0; trans = 0.0_dp

        ! Find "sym_ops"
        p = index(json, '"sym_ops"')
        if (p == 0) then; ios = 1; return; end if
        ! Find first '[' after sym_ops
        p = index(json(p:), '[') + p - 1
        if (p <= 0) then; ios = 1; return; end if

        ! Parse each {"rot":[...],"trans":[...]} object
        iop = 0
        q = p + 1
        do while (iop < max_ops)
            p1 = index(json(q:), '"rot"')
            if (p1 == 0) exit
            p1 = p1 + q - 1
            p1 = index(json(p1:), '[') + p1
            p2 = index(json(p1:), ']') + p1 - 2
            if (p1 <= 0 .or. p2 <= p1) exit
            iop = iop + 1
            read(json(p1:p2), *, iostat=ios) &
                rot(1,1,iop), rot(1,2,iop), rot(1,3,iop), &
                rot(2,1,iop), rot(2,2,iop), rot(2,3,iop), &
                rot(3,1,iop), rot(3,2,iop), rot(3,3,iop)
            ! Find trans array
            p1 = index(json(p2:), '"trans"')
            if (p1 == 0) exit
            p1 = p1 + p2 - 1
            p1 = index(json(p1:), '[') + p1
            p2 = index(json(p1:), ']') + p1 - 2
            if (p1 <= 0 .or. p2 <= p1) exit
            read(json(p1:p2), *, iostat=ios) trans(1,iop), trans(2,iop), trans(3,iop)
            q = p2 + 2
        end do
        n_symops = iop
    end subroutine parse_sym_ops_from_json_str




    ! ═══════════════════════════════════════════════════════════════════════════


    ! ═══════════════════════════════════════════════════════════════════════════


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

    function build_metadata_json(grid, cfg, sym_ops, n_symops_stored) result(json)
        type(pes_grid_t), intent(in) :: grid
        type(castep_config_t), intent(in), target :: cfg
        type(sym_op_t), intent(in), optional :: sym_ops(:)
        integer, intent(in), optional :: n_symops_stored
        character(len=16384) :: json
        character(len=128) :: lat_str, plane_l, mob_el, sym_buf
        character(len=5) :: sym_str
        integer :: mi, iop, n_sym

        mi = grid%mobile_atom_idx
        mob_el = trim(cfg%atom_type(mi))

        write(lat_str, '(a,f12.6,a,f12.6,a,f12.6,a,f10.4,a,f10.4,a,f10.4)') &
            '"a":', cfg%cell_length(1), ',"b":', cfg%cell_length(2), &
            ',"c":', cfg%cell_length(3), ',"alpha":', cfg%cell_angle(1), &
            ',"beta":', cfg%cell_angle(2), ',"gamma":', cfg%cell_angle(3)

        plane_l = plane_name(grid)
        if (grid%use_symmetry) then
            sym_str = 'true '
        else
            sym_str = 'false'
        end if

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
            ! 3D: include sym_ops in JSON for self-contained expansion
            n_sym = 0
            if (present(sym_ops) .and. present(n_symops_stored)) n_sym = n_symops_stored
            if (n_sym > 0) then
                ! Build sym_ops JSON array inline
                sym_buf = ' '
                write(json, '(a,a,a,i0,a,a,a,a,f12.8,a,f12.8,a,a,f12.8,a,f12.8,a,a,f12.8,a,f12.8,a,a,a,a,a,a,a,i0,a)') &
                    '{"type":"pes_3d","scan_mode":"', trim(grid%scan_mode), &
                    '","mobile_idx":', mi - 1, &
                    ',"mobile_el":"', trim(mob_el), '",', &
                    '"fx_range":[', grid%frac_range(1,1), ',', grid%frac_range(1,2), '],', &
                    '"fy_range":[', grid%frac_range(2,1), ',', grid%frac_range(2,2), '],', &
                    '"fz_range":[', grid%frac_range(3,1), ',', grid%frac_range(3,2), '],', &
                    '"lattice":{', trim(lat_str), '},', &
                    '"use_symmetry":', trim(sym_str), ',"n_symops":', n_sym, ',"sym_ops":['
                ! Append each sym_op
                do iop = 1, n_sym
                    if (iop < n_sym) then
                        write(sym_buf, '(a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,f10.8,a,f10.8,a,f10.8,a)') &
                            '{"rot":[', sym_ops(iop)%rot(1,1), ',', sym_ops(iop)%rot(1,2), ',', sym_ops(iop)%rot(1,3), ',', &
                            sym_ops(iop)%rot(2,1), ',', sym_ops(iop)%rot(2,2), ',', sym_ops(iop)%rot(2,3), ',', &
                            sym_ops(iop)%rot(3,1), ',', sym_ops(iop)%rot(3,2), ',', sym_ops(iop)%rot(3,3), &
                            '],"trans":[', sym_ops(iop)%trans(1), ',', sym_ops(iop)%trans(2), ',', sym_ops(iop)%trans(3), ']},'
                        json = trim(json) // trim(sym_buf)
                    else
                        write(sym_buf, '(a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,f10.8,a,f10.8,a,f10.8,a)') &
                            '{"rot":[', sym_ops(iop)%rot(1,1), ',', sym_ops(iop)%rot(1,2), ',', sym_ops(iop)%rot(1,3), ',', &
                            sym_ops(iop)%rot(2,1), ',', sym_ops(iop)%rot(2,2), ',', sym_ops(iop)%rot(2,3), ',', &
                            sym_ops(iop)%rot(3,1), ',', sym_ops(iop)%rot(3,2), ',', sym_ops(iop)%rot(3,3), &
                            '],"trans":[', sym_ops(iop)%trans(1), ',', sym_ops(iop)%trans(2), ',', sym_ops(iop)%trans(3), ']}]}'
                        json = trim(json) // trim(sym_buf)
                    end if
                end do
            else
                write(json, '(a,a,a,i0,a,a,a,a,f12.8,a,f12.8,a,a,f12.8,a,f12.8,a,a,f12.8,a,f12.8,a,a,a,a,a,a,a)') &
                    '{"type":"pes_3d","scan_mode":"', trim(grid%scan_mode), &
                    '","mobile_idx":', mi - 1, &
                    ',"mobile_el":"', trim(mob_el), '",', &
                    '"fx_range":[', grid%frac_range(1,1), ',', grid%frac_range(1,2), '],', &
                    '"fy_range":[', grid%frac_range(2,1), ',', grid%frac_range(2,2), '],', &
                    '"fz_range":[', grid%frac_range(3,1), ',', grid%frac_range(3,2), '],', &
                    '"lattice":{', trim(lat_str), '},', &
                    '"use_symmetry":', trim(sym_str), '}'
            end if
        end if
    end function build_metadata_json


    ! ── IEEE NaN helper ──

    function ieee_nan() result(val)
        real(dp) :: val
        ! Portable way to generate quiet NaN without ieee_arithmetic module
        val = 0.0_dp
        val = val / val   ! 0/0 → NaN on virtually all IEEE 754 systems
    end function ieee_nan



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

    subroutine patch_json_full_cell(json_in, json_out)
        !! Patch JSON line 2 for expanded cube: set all frac ranges to [0,1]
        !! and remove symmetry flag (already expanded).
        character(len=*), intent(in) :: json_in
        character(len=*), intent(out) :: json_out
        integer :: kp

        json_out = json_in

        ! Replace fx_range values
        call replace_json_array(json_out, '"fx_range"', 0.0_dp, 1.0_dp)
        call replace_json_array(json_out, '"fy_range"', 0.0_dp, 1.0_dp)
        call replace_json_array(json_out, '"fz_range"', 0.0_dp, 1.0_dp)

        ! Change use_symmetry:true → use_symmetry:false. 'false' is one
        ! character longer than 'true', so this must SPLICE the string — the
        ! old in-place overwrite clobbered the following '"' and produced
        ! `"use_symmetry":false"n_symops":…`, invalid JSON that made
        ! cube_reader drop the PES metadata (2D/3D detection broken).
        kp = index(json_out, '"use_symmetry":true')
        if (kp > 0) then
            json_out = json_out(:kp-1) // '"use_symmetry":false' // json_out(kp+19:)
        end if
    end subroutine patch_json_full_cell

    subroutine replace_json_array(json, key, v1, v2)
        !! Replace [old1, old2] with [v1, v2] for a given JSON key.
        character(len=*), intent(inout) :: json
        character(len=*), intent(in) :: key
        real(dp), intent(in) :: v1, v2
        integer :: kp, bp, ep
        character(len=48) :: new_val
        kp = index(json, trim(key))
        if (kp == 0) return
        bp = index(json(kp:), '[') + kp - 1
        ep = index(json(kp:), ']') + kp - 1
        if (bp <= 0 .or. ep <= bp) return
        write(new_val, '(f12.8,a,f12.8)') v1, ',', v2
        json(bp+1:ep-1) = trim(new_val)
    end subroutine replace_json_array




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
        integer :: unit, eq_pos, n_found

        energy = 0.0_dp; ios = 0
        n_found = 0
        open(newunit=unit, file=trim(castep_file), status='old', action='read', iostat=ios)
        if (ios /= 0) return
        do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            ! CASTEP appends every run of the same seed to the existing .castep,
            ! so a file can accumulate several SCF histories (an old run plus
            ! one or more reruns). Always take the LAST occurrence:
            !   1. "Total energy corrected for finite basis set" — the final
            !      converged value (finite-basis extrapolation + dispersion)
            !   2. fallback: plain "Final energy" (last basis size), skipping
            !      "Dispersion corrected final energy*" (vdW-inclusive variant)
            if (index(line, 'corrected for finite basis set') > 0) then
                eq_pos = index(line, '=')
                if (eq_pos > 0) then
                    read(line(eq_pos+1:), *, iostat=ios) energy
                    if (ios == 0) n_found = n_found + 1
                end if
            else if ((index(line, 'Final energy') > 0 .or. &
                      index(line, 'final energy') > 0 .or. &
                      index(line, 'Final Energy') > 0) .and. &
                     index(line, 'ispersion') == 0) then
                eq_pos = index(line, '=')
                if (eq_pos > 0) then
                    read(line(eq_pos+1:), *, iostat=ios) energy
                    if (ios == 0) n_found = n_found + 1
                end if
            end if
        end do
        close(unit)
        if (n_found == 0) then
            ios = 1   ! no parseable final energy found
        else
            ios = 0   ! success — the loop's EOF (-1) must not leak out
        end if
    end function parse_castep_energy


    ! ── Fractional coordinate wrapping ──

    pure subroutine wrap_to_unit(x)
        real(dp), intent(inout) :: x
        real(dp), parameter :: EPS = 1.0d-5
        x = x - floor(x + EPS)
        if (x < 0.0_dp) x = x + 1.0_dp
        if (x >= 1.0_dp .or. x < 0.0_dp) x = x - aint(x)
    end subroutine wrap_to_unit

end module pes
