module pes_scan
    !! Potential Energy Surface scan — batch file generation and result collection.
    !!
    !! Generates .cell + .param input files for a 2D grid scan of a mobile atom
    !! on a selected crystallographic plane (xy, xz, or yz) with optional
    !! constrained geometry optimization. Collects CASTEP output energies into
    !! a JSON metadata file for 3D visualization in crystal-viewer.
    use castep_config, only: dp, castep_config_t, IO_WRITE_FAIL, &
        IO_PARSE_ERROR, IO_FILE_NOT_FOUND
    implicit none
    private

    public :: pes_grid_t
    public :: generate_pes_grid_points
    public :: write_pes_metadata_json
    public :: collect_pes_energies

    ! Maximum grid size per dimension
    integer, parameter :: MAX_GRID = 200

    type :: pes_grid_t
        integer  :: plane_axis(2)   = [1, 2]     ! scan plane: (1,2)=xy, (1,3)=xz, (2,3)=yz
        integer  :: n_points(2)     = [5, 5]     ! Nx, Ny
        real(dp) :: frac_range(2,2) = 0.0_dp     ! [fx_min,fx_max], [fy_min,fy_max]
        integer  :: mobile_atom_idx = 1           ! 1-based index of mobile atom
        character(len=8) :: scan_mode = 'SP'      ! 'SP' or 'RELAX'
    end type pes_grid_t

contains

    ! ── Grid point generation ──

    subroutine generate_pes_grid_points(grid, frac_points, n_total, iostat, iomsg)
        !! Generate all fractional coordinate pairs for the 2D scan grid.
        !! Returns frac_points(n_total, 2) — columns are (fx, fy) on the scan plane.
        type(pes_grid_t), intent(in) :: grid
        real(dp), allocatable, intent(out) :: frac_points(:,:)
        integer, intent(out) :: n_total
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        integer :: i, j, idx
        real(dp) :: dx, dy, fx, fy

        iostat = 0
        n_total = 0

        if (grid%n_points(1) < 2 .or. grid%n_points(2) < 2) then
            iostat = 1
            if (present(iomsg)) iomsg = 'Grid must have at least 2 points per dimension'
            return
        end if
        if (grid%n_points(1) > MAX_GRID .or. grid%n_points(2) > MAX_GRID) then
            iostat = 2
            if (present(iomsg)) iomsg = 'Grid size exceeds maximum'
            return
        end if

        n_total = grid%n_points(1) * grid%n_points(2)
        allocate(frac_points(n_total, 2), stat=iostat)
        if (iostat /= 0) then
            if (present(iomsg)) iomsg = 'Memory allocation failed for grid points'
            return
        end if

        dx = 0.0_dp
        dy = 0.0_dp
        if (grid%n_points(1) > 1) dx = (grid%frac_range(1,2) - grid%frac_range(1,1)) &
                                       / real(grid%n_points(1) - 1, dp)
        if (grid%n_points(2) > 1) dy = (grid%frac_range(2,2) - grid%frac_range(2,1)) &
                                       / real(grid%n_points(2) - 1, dp)

        idx = 0
        do j = 0, grid%n_points(2) - 1
            fy = grid%frac_range(2,1) + j * dy
            do i = 0, grid%n_points(1) - 1
                fx = grid%frac_range(1,1) + i * dx
                idx = idx + 1
                frac_points(idx, 1) = fx
                frac_points(idx, 2) = fy
            end do
        end do
    end subroutine generate_pes_grid_points


    ! ── JSON metadata ──

    subroutine write_pes_metadata_json(json_path, grid, cfg, iostat, iomsg)
        !! Write pes_metadata.json describing the scan grid, lattice, and atom info.
        !! Energies are written as null placeholders.
        character(len=*), intent(in) :: json_path
        type(pes_grid_t), intent(in) :: grid
        type(castep_config_t), intent(in) :: cfg
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        integer :: unit, ios, i, n_total
        character(len=8) :: plane_str

        iostat = 0

        ! Determine plane name
        if (grid%plane_axis(1) == 1 .and. grid%plane_axis(2) == 2) then
            plane_str = 'xy'
        else if (grid%plane_axis(1) == 1 .and. grid%plane_axis(2) == 3) then
            plane_str = 'xz'
        else
            plane_str = 'yz'
        end if

        open(newunit=unit, file=trim(json_path), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            iostat = IO_WRITE_FAIL
            if (present(iomsg)) iomsg = 'Cannot write: ' // trim(json_path)
            return
        end if

        write(unit, '(a)') '{'
        write(unit, '(a)') '  "type": "pes_scan",'
        write(unit, '(a,a,a)') '  "plane": "', trim(plane_str), '",'
        write(unit, '(a,i0,a,i0,a)') '  "nx": ', grid%n_points(1), ', "ny": ', grid%n_points(2), ','

        write(unit, '(a, f12.8, a, f12.8, a)') &
            '  "fx_range": [', grid%frac_range(1,1), ', ', grid%frac_range(1,2), '],'
        write(unit, '(a, f12.8, a, f12.8, a)') &
            '  "fy_range": [', grid%frac_range(2,1), ', ', grid%frac_range(2,2), '],'

        write(unit, '(a,i0,a,a,a)') '  "mobile_atom": { "index": ', grid%mobile_atom_idx - 1, &
            ', "element": "', trim(cfg%atom_type(grid%mobile_atom_idx)), '" },'
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

        ! Energy grid: write null entries
        write(unit, '(a)') '  "energies": ['
        n_total = grid%n_points(1) * grid%n_points(2)
        do i = 1, n_total
            if (i < n_total) then
                write(unit, '(a)') '    null,'
            else
                write(unit, '(a)') '    null'
            end if
        end do
        write(unit, '(a)') '  ],'
        write(unit, '(a)') '  "has_energies": false'
        write(unit, '(a)') '}'

        close(unit)
    end subroutine write_pes_metadata_json


    ! ── Result collection ──

    subroutine collect_pes_energies(scan_dir, iostat, iomsg)
        !! Scan all grid_*/ subdirectories, parse .castep files, fill energies
        !! into pes_metadata.json, and set has_energies = true.
        character(len=*), intent(in) :: scan_dir
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        character(len=1024) :: json_path, grid_dir, castep_file, line
        integer :: i, j, nx, ny, n_total, collected, missing, ios
        integer :: json_unit
        real(dp), allocatable :: energies(:)
        logical, allocatable :: has_energy(:)
        logical :: exists
        real(dp) :: e_val, e_min, e_max

        iostat = 0
        collected = 0; missing = 0
        e_min = huge(1.0_dp); e_max = -huge(1.0_dp)

        json_path = trim(scan_dir) // '/pes_metadata.json'
        inquire(file=trim(json_path), exist=exists)
        if (.not. exists) then
            iostat = IO_FILE_NOT_FOUND
            if (present(iomsg)) iomsg = 'pes_metadata.json not found in: ' // trim(scan_dir)
            return
        end if

        ! Read nx, ny from JSON (simple line-based parsing)
        nx = 0; ny = 0
        open(newunit=json_unit, file=trim(json_path), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            iostat = IO_PARSE_ERROR
            if (present(iomsg)) iomsg = 'Cannot open: ' // trim(json_path)
            return
        end if
        do
            read(json_unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            ! "nx": 5, "ny": 5,  — may be on same line or separate lines
            if (index(line, '"nx"') > 0 .and. nx <= 0) &
                call extract_json_int(line, '"nx"', nx, ios)
            if (index(line, '"ny"') > 0 .and. ny <= 0) &
                call extract_json_int(line, '"ny"', ny, ios)
        end do
        close(json_unit)

        if (nx <= 0 .or. ny <= 0) then
            iostat = IO_PARSE_ERROR
            if (present(iomsg)) iomsg = 'Failed to read grid dimensions from JSON'
            return
        end if

        n_total = nx * ny
        allocate(energies(n_total), stat=ios)
        allocate(has_energy(n_total), stat=iostat)
        if (ios /= 0 .or. iostat /= 0) then
            iostat = 1
            if (present(iomsg)) iomsg = 'Memory allocation failed'
            return
        end if
        energies = 0.0_dp
        has_energy = .false.

        ! Scan grid directories
        do j = 1, ny
            do i = 1, nx
                write(grid_dir, '(a, a, i3.3, a, i3.3)') trim(scan_dir), '/grid_', i, '_', j
                castep_file = find_castep_in_dir(grid_dir)
                if (len_trim(castep_file) == 0) then
                    missing = missing + 1
                    cycle
                end if

                e_val = parse_castep_energy(castep_file, ios)
                if (ios == 0) then
                    energies((j-1)*nx + i) = e_val
                    has_energy((j-1)*nx + i) = .true.
                    collected = collected + 1
                    if (e_val < e_min) e_min = e_val
                    if (e_val > e_max) e_max = e_val
                else
                    missing = missing + 1
                end if
            end do
        end do

        ! Write updated JSON with energies (null for missing entries)
        call rewrite_json_with_energies(json_path, energies, has_energy, nx, ny, ios)
        if (ios /= 0) then
            iostat = IO_WRITE_FAIL
            if (present(iomsg)) iomsg = 'Failed to rewrite pes_metadata.json with energies'
            deallocate(energies); deallocate(has_energy)
            return
        end if

        write(*, '(a)') '  ── Collection Summary ──'
        write(*, '(a, i0, a, i0)') '  Collected: ', collected, ' / ', n_total
        if (missing > 0) write(*, '(a, i0)') '  Missing:   ', missing
        if (collected > 0) then
            write(*, '(a, f18.8)') '  E min (eV): ', e_min
            write(*, '(a, f18.8)') '  E max (eV): ', e_max
        end if

        deallocate(energies)
        deallocate(has_energy)
    end subroutine collect_pes_energies


    ! ── Private helpers ──

    function find_castep_in_dir(dir_path) result(castep_path)
        !! Find a .castep file: hardcoded scan.castep first, then wildcard ls.
        character(len=*), intent(in) :: dir_path
        character(len=1024) :: castep_path, cmd, tmp_file
        logical :: exists
        integer :: unit, ios

        ! 1) Hardcoded: scan.castep
        castep_path = trim(dir_path) // '/scan.castep'
        inquire(file=trim(castep_path), exist=exists)
        if (exists) return

        ! 2) Wildcard fallback via ls > tmp file
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
        ios = 1  ! not found
        close(unit)
    end function parse_castep_energy


    subroutine rewrite_json_with_energies(json_path, energies, has_energy, nx, ny, ios)
        !! Rewrite pes_metadata.json with collected energy values.
        !! Missing entries (has_energy == .false.) are written as null.
        character(len=*), intent(in) :: json_path
        real(dp), intent(in) :: energies(:)
        logical, intent(in) :: has_energy(:)
        integer, intent(in) :: nx, ny
        integer, intent(out) :: ios

        character(len=1024) :: tmp_path
        integer :: unit_in, unit_out, j, n_total
        character(len=4096) :: line
        logical :: in_energies, wrote_energy

        tmp_path = trim(json_path) // '.tmp'

        open(newunit=unit_in, file=trim(json_path), status='old', action='read', iostat=ios)
        if (ios /= 0) return

        open(newunit=unit_out, file=trim(tmp_path), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            close(unit_in); return
        end if

        n_total = nx * ny
        in_energies = .false.
        wrote_energy = .false.

        do
            read(unit_in, '(a)', iostat=ios) line
            if (ios /= 0) exit

            if (index(line, '"energies"') > 0) then
                in_energies = .true.
                write(unit_out, '(a)') trim(line)
                cycle
            end if

            if (in_energies .and. .not. wrote_energy) then
                ! Skip old energy lines until we find ']'
                if (index(line, ']') > 0 .and. index(line, '"has_energies"') == 0) then
                    ! Write new energy values (null for missing)
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
                    cycle  ! skip old energy values / empty entries
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

        ! Replace original with tmp
        call execute_command_line('mv "' // trim(tmp_path) // '" "' // trim(json_path) // '"', exitstat=ios)
    end subroutine rewrite_json_with_energies

    subroutine extract_json_int(line, key, val, ios)
        !! Extract integer value for a JSON key from a line.
        !! Example: line = '"nx": 5, "ny": 3,'  key='"nx"' → val=5
        character(len=*), intent(in) :: line, key
        integer, intent(out) :: val, ios
        integer :: kp, cp
        kp = index(line, trim(key))
        if (kp == 0) then; ios = 1; return; end if
        cp = index(line(kp:), ':')
        if (cp == 0) then; ios = 1; return; end if
        read(line(kp+cp:), *, iostat=ios) val
    end subroutine extract_json_int

end module pes_scan
