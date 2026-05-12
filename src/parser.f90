module parser
    !! File format parsers for CIF, PDB, and CASTEP .cell files
    use castep_config, only: &
         dp, MAX_ATOMS, IO_FILE_NOT_FOUND, IO_PARSE_ERROR, &
         cif_data_t, atom_t, pi, compare_tags
    implicit none
    private

    public :: parse_cif_inline
    public :: parse_pdb_inline
    public :: parse_cell_inline
    public :: clean_element_symbol

contains

    subroutine grow_atoms(atoms, current, new_size)
        !! Grow the atoms array to new_size, preserving all existing entries.
        type(atom_t), allocatable, intent(inout) :: atoms(:)
        integer, intent(in) :: current, new_size
        type(atom_t), allocatable :: tmp(:)

        if (.not. allocated(atoms)) then
            allocate(atoms(new_size))
        else
            allocate(tmp(new_size))
            tmp(1:current) = atoms(1:current)
            call move_alloc(tmp, atoms)
        end if
    end subroutine grow_atoms

    subroutine parse_pdb_inline(filename, data, iostat, iomsg)
        !! Parse PDB format: CRYST1 record for cell, ATOM/HETATM records for atoms
        !! All coordinates are Cartesian (Angstrom).
        character(len=*), intent(in) :: filename
        type(cif_data_t), intent(out) :: data
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        integer :: ios, iunit, n_atoms
        character(len=256) :: line
        logical :: file_found

        iostat = 0
        data%a = 0.0_dp; data%b = 0.0_dp; data%c = 0.0_dp
        data%alpha = 90.0_dp; data%beta = 90.0_dp; data%gamma = 90.0_dp
        data%space_group = 'P1'
        data%n_atoms = 0
        n_atoms = 0
        if (allocated(data%atoms)) deallocate(data%atoms)

        inquire(file=trim(filename), exist=file_found)
        if (.not. file_found) then
            iostat = IO_FILE_NOT_FOUND
            if (present(iomsg)) iomsg = 'File not found: ' // trim(filename)
            return
        end if

        open(newunit=iunit, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            iostat = IO_PARSE_ERROR
            if (present(iomsg)) iomsg = 'Cannot open file: ' // trim(filename)
            return
        end if

        do
            read(iunit, '(A)', iostat=ios) line
            if (ios /= 0) exit

            ! Handle CRYST1 record (fixed-width columns per PDB spec)
            if (line(1:6) == 'CRYST1') then
                read(line(7:15), '(F9.3)') data%a
                read(line(16:24), '(F9.3)') data%b
                read(line(25:33), '(F9.3)') data%c
                read(line(34:40), '(F7.2)') data%alpha
                read(line(41:47), '(F7.2)') data%beta
                read(line(48:54), '(F7.2)') data%gamma
                ! Space group: columns 56-66
                if (len_trim(line) >= 66) then
                    data%space_group = line(56:66)
                end if
                cycle
            end if

            ! Handle UNITCELL record (orthorhombic: a12=a13=a23=0)
            if (line(1:8) == 'UNITCELL') then
                read(line(9:17), '(F9.3)') data%a
                read(line(18:26), '(F9.3)') data%b
                read(line(27:35), '(F9.3)') data%c
                ! UNITCELL implies orthorhombic (90-degree angles)
                data%alpha = 90.0_dp; data%beta = 90.0_dp; data%gamma = 90.0_dp
                cycle
            end if

            ! Handle ATOM and HETATM records (fixed-width columns per PDB spec)
            if (line(1:4) == 'ATOM' .or. line(1:7) == 'HETATM') then
                n_atoms = n_atoms + 1
                if (n_atoms > MAX_ATOMS) exit

                if (.not. allocated(data%atoms)) then
                    call grow_atoms(data%atoms, 0, 64)
                else if (n_atoms > size(data%atoms)) then
                    call grow_atoms(data%atoms, n_atoms - 1, max(2 * size(data%atoms), n_atoms))
                end if

                ! Element symbol: columns 77-78 (right-justified)
                data%atoms(n_atoms)%element = line(77:78)
                data%atoms(n_atoms)%element = adjustl(data%atoms(n_atoms)%element)

                ! Atom name: columns 13-16
                data%atoms(n_atoms)%label = adjustl(line(13:16))

                ! Cartesian coordinates: columns 31-38, 39-46, 47-54
                read(line(31:38), '(F8.3)', iostat=ios) data%atoms(n_atoms)%x
                read(line(39:46), '(F8.3)', iostat=ios) data%atoms(n_atoms)%y
                read(line(47:54), '(F8.3)', iostat=ios) data%atoms(n_atoms)%z
            end if
        end do

        close(iunit)
        data%n_atoms = n_atoms
    end subroutine parse_pdb_inline

    pure function compute_abc_from_cartesian(cart) result(abc)
        !! Convert Cartesian lattice vectors (column-major 3x3) back to
        !! lattice lengths (Angstrom) and angles (degrees).
        real(dp), intent(in) :: cart(3, 3)
        real(dp) :: abc(6)  ! a, b, c, alpha, beta, gamma
        real(dp) :: dot_ab, dot_ac, dot_bc
        real(dp) :: cos_gamma, cos_beta, cos_alpha

        ! Magnitudes
        abc(1) = dsqrt(cart(1,1)**2 + cart(2,1)**2 + cart(3,1)**2)
        abc(2) = dsqrt(cart(1,2)**2 + cart(2,2)**2 + cart(3,2)**2)
        abc(3) = dsqrt(cart(1,3)**2 + cart(2,3)**2 + cart(3,3)**2)

        ! Dot products
        dot_ab = cart(1,1)*cart(1,2) + cart(2,1)*cart(2,2) + cart(3,1)*cart(3,2)
        dot_ac = cart(1,1)*cart(1,3) + cart(2,1)*cart(2,3) + cart(3,1)*cart(3,3)
        dot_bc = cart(1,2)*cart(1,3) + cart(2,2)*cart(2,3) + cart(3,2)*cart(3,3)

        ! cos(gamma) = dot(a_vec, b_vec) / (|a|*|b|)
        if (abc(1) > 0.0_dp .and. abc(2) > 0.0_dp) then
            cos_gamma = dmax1(-1.0_dp, dmin1(1.0_dp, dot_ab / (abc(1)*abc(2))))
            abc(6) = dacos(cos_gamma) * 180.0_dp / pi
        else
            abc(6) = 90.0_dp
        end if
        ! cos(beta) = dot(a_vec, c_vec) / (|a|*|c|)
        if (abc(1) > 0.0_dp .and. abc(3) > 0.0_dp) then
            cos_beta = dmax1(-1.0_dp, dmin1(1.0_dp, dot_ac / (abc(1)*abc(3))))
            abc(5) = dacos(cos_beta) * 180.0_dp / pi
        else
            abc(5) = 90.0_dp
        end if
        ! cos(alpha) = dot(b_vec, c_vec) / (|b|*|c|)
        if (abc(2) > 0.0_dp .and. abc(3) > 0.0_dp) then
            cos_alpha = dmax1(-1.0_dp, dmin1(1.0_dp, dot_bc / (abc(2)*abc(3))))
            abc(4) = dacos(cos_alpha) * 180.0_dp / pi
        else
            abc(4) = 90.0_dp
        end if
    end function compute_abc_from_cartesian

    subroutine parse_cell_inline(filename, data, iostat, iomsg)
        !! Parse CASTEP .cell format: %BLOCK LATTICE_ABC for cell,
        !! %BLOCK POSITIONS_ABS for atoms (Cartesian coords).
        character(len=*), intent(in) :: filename
        type(cif_data_t), intent(out) :: data
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        integer :: ios, iunit, n_atoms, ios_local, n_cell_cols
        character(len=256) :: line, trimmed
        logical :: file_found, in_lattice, in_positions, in_lattice_cart
        logical :: read_a_line, read_angle_line, cart_seen, skip_units = .false.
        integer :: cart_line
        character(len=6) :: elem
        real(dp) :: x, y, z
        real(dp) :: lattice_cart(3, 3), abc(6)
        character(len=256) :: row_vals_cell(50)

        iostat = 0
        data%a = 0.0_dp; data%b = 0.0_dp; data%c = 0.0_dp
        data%alpha = 90.0_dp; data%beta = 90.0_dp; data%gamma = 90.0_dp
        data%space_group = 'P1'
        data%n_atoms = 0
        in_lattice = .false.
        in_lattice_cart = .false.
        in_positions = .false.
        read_a_line = .false.
        read_angle_line = .false.
        cart_seen = .false.
        cart_line = 0
        n_atoms = 0
        if (allocated(data%atoms)) deallocate(data%atoms)

        inquire(file=trim(filename), exist=file_found)
        if (.not. file_found) then
            iostat = IO_FILE_NOT_FOUND
            if (present(iomsg)) iomsg = 'File not found: ' // trim(filename)
            return
        end if

        open(newunit=iunit, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            iostat = IO_PARSE_ERROR
            if (present(iomsg)) iomsg = 'Cannot open file: ' // trim(filename)
            return
        end if

        do
            read(iunit, '(A)', iostat=ios) line
            if (ios /= 0) exit

            trimmed = adjustl(trim(line))

            ! Skip empty lines
            if (len_trim(trimmed) == 0) cycle
            ! Skip comment lines
            if (trimmed(1:1) == '!' .or. trimmed(1:1) == '#') cycle

            ! Detect LATTICE_ABC block (case-insensitive)
            if (compare_tags(trimmed, '%BLOCK LATTICE_ABC')) then
                in_lattice = .true.
                read_a_line = .false.
                read_angle_line = .false.
                cycle
            end if
            if (compare_tags(trimmed, '%ENDBLOCK LATTICE_ABC')) then
                in_lattice = .false.
                cycle
            end if

            ! Detect LATTICE_CART block (case-insensitive)
            if (compare_tags(trimmed, '%BLOCK LATTICE_CART')) then
                in_lattice_cart = .true.
                cart_line = 0
                skip_units = .true.  ! next line may be ANG/BOHR
                lattice_cart = 0.0_dp
                cycle
            end if
            if (compare_tags(trimmed, '%ENDBLOCK LATTICE_CART')) then
                abc = compute_abc_from_cartesian(lattice_cart)
                data%a = abc(1); data%b = abc(2); data%c = abc(3)
                data%alpha = abc(4); data%beta = abc(5); data%gamma = abc(6)
                in_lattice_cart = .false.
                cart_seen = .true.
                cycle
            end if

            ! Detect POSITIONS_ABS or POSITIONS_FRAC blocks (case-insensitive)
            if (compare_tags(trimmed, '%BLOCK POSITIONS_ABS') .or. &
                compare_tags(trimmed, '%BLOCK POSITIONS_FRAC')) then
                in_positions = .true.
                cycle
            end if
            if (compare_tags(trimmed, '%ENDBLOCK POSITIONS_ABS') .or. &
                compare_tags(trimmed, '%ENDBLOCK POSITIONS_FRAC')) then
                in_positions = .false.
                cycle
            end if

            ! Skip other % directives (ENDBLOCK, etc.)
            if (trimmed(1:1) == '%') cycle

            ! Parse LATTICE_ABC content (skip if LATTICE_CART was already seen)
            if (in_lattice .and. .not. cart_seen) then
                if (.not. read_a_line) then
                    read(trimmed, '(3F12.7)', iostat=ios_local) data%a, data%b, data%c
                    if (ios_local == 0) read_a_line = .true.
                else if (.not. read_angle_line) then
                    read(trimmed, '(3F12.7)', iostat=ios_local) data%alpha, data%beta, data%gamma
                    if (ios_local == 0) read_angle_line = .true.
                end if
                cycle
            end if

             ! Parse LATTICE_CART: 3 lines of 3 values each
            ! Each file row = one lattice vector (row-major in file).
            ! Need to transpose to column-major for compute_abc_from_cartesian:
            ! row 1 (a_x,a_y,a_z) -> cart(:,1), row 2 (b_x,b_y,b_z) -> cart(:,2), etc.
            if (in_lattice_cart) then
                if (skip_units) then
                    skip_units = .false.
                    if (trimmed == 'ANG' .or. trimmed == 'BOHR' .or. &
                        trimmed == 'ang' .or. trimmed == 'bohr') cycle
                end if
                cart_line = cart_line + 1
                ios_local = 0
                read(trimmed, *, iostat=ios_local) x, y, z
                if (ios_local == 0) then
                    ! Transpose: file row i -> cart(:,i)
                    lattice_cart(cart_line, 1) = x
                    lattice_cart(cart_line, 2) = y
                    lattice_cart(cart_line, 3) = z
                end if
                if (cart_line >= 3) then
                    abc = compute_abc_from_cartesian(lattice_cart)
                    data%a = abc(1); data%b = abc(2); data%c = abc(3)
                    data%alpha = abc(4); data%beta = abc(5); data%gamma = abc(6)
                    in_lattice_cart = .false.
                    cart_seen = .true.
                end if
                cycle
            end if

            ! Parse POSITIONS_ABS content: Element X Y Z (token-based)
            if (in_positions) then
                ios_local = 0
                call tokenize_inline(trimmed, row_vals_cell, n_cell_cols)
                if (n_cell_cols >= 4) then
                    call copy_str_no_quotes(row_vals_cell(1), elem)
                    read(row_vals_cell(2), '(F12.7)', iostat=ios_local) x
                    read(row_vals_cell(3), '(F12.7)', iostat=ios_local) y
                    read(row_vals_cell(4), '(F12.7)', iostat=ios_local) z
                end if
                if (ios_local == 0 .and. n_cell_cols >= 4) then
                    n_atoms = n_atoms + 1
                    if (n_atoms <= MAX_ATOMS) then
                        if (.not. allocated(data%atoms)) then
                            call grow_atoms(data%atoms, 0, 64)
                        else if (n_atoms > size(data%atoms)) then
                            call grow_atoms(data%atoms, n_atoms - 1, max(2 * size(data%atoms), n_atoms))
                        end if
                        data%atoms(n_atoms)%label = adjustl(elem)
                        data%atoms(n_atoms)%element = adjustl(elem)
                        data%atoms(n_atoms)%x = x
                        data%atoms(n_atoms)%y = y
                        data%atoms(n_atoms)%z = z
                    end if
                end if
                cycle
            end if
        end do

        close(iunit)
        data%n_atoms = n_atoms
    end subroutine parse_cell_inline

    subroutine parse_cif_inline(filename, data, iostat, iomsg)
        !! CIF parser: tag-value pairs and loop_ blocks with atom_site data
        character(len=*), intent(in) :: filename
        type(cif_data_t), intent(out) :: data
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        integer :: ios, iunit, i, n, n_cols, line_len, start_pos, end_pos
        integer :: idx_col(5)
        character(len=256) :: line, tag, val
        logical :: in_atom_loop, in_loop, file_found
        character(len=256) :: row_vals(50)

        iostat = 0
        data%a = 0.0_dp; data%b = 0.0_dp; data%c = 0.0_dp
        data%alpha = 90.0_dp; data%beta = 90.0_dp; data%gamma = 90.0_dp
        data%space_group = 'P1'
        data%n_atoms = 0
        idx_col(1:5) = 0
        if (allocated(data%atoms)) deallocate(data%atoms)

        inquire(file=trim(filename), exist=file_found)
        if (.not. file_found) then
            iostat = IO_FILE_NOT_FOUND
            if (present(iomsg)) iomsg = 'File not found: ' // trim(filename)
            return
        end if

        open(newunit=iunit, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            iostat = IO_PARSE_ERROR
            if (present(iomsg)) iomsg = 'Cannot open file: ' // trim(filename)
            return
        end if

        in_atom_loop = .false.
        in_loop = .false.
        n = 0
        n_cols = 0

        do
            read(iunit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = trim(adjustl(line))
            line_len = len_trim(line)
            if (line_len == 0 .or. line(1:1) == '#') cycle

            if (trim(line) == 'loop_') then
                in_atom_loop = .false.
                in_loop = .true.
                n_cols = 0
                idx_col(1:5) = 0
                cycle
            end if

            ! Inside a loop block
            if (in_loop) then
                ! Underscore lines are column headers
                if (line_len >= 1 .and. line(1:1) == '_') then
                    if (index(line, '_atom_site_label') > 0)       idx_col(1) = n_cols + 1
                    if (index(line, '_atom_site_type_symbol') > 0) idx_col(2) = n_cols + 1
                    if (index(line, '_atom_site_fract_x') > 0)     idx_col(3) = n_cols + 1
                    if (index(line, '_atom_site_fract_y') > 0)     idx_col(4) = n_cols + 1
                    if (index(line, '_atom_site_fract_z') > 0)     idx_col(5) = n_cols + 1
                    if (idx_col(3) > 0 .and. idx_col(4) > 0 .and. idx_col(5) > 0) then
                        in_atom_loop = .true.
                    end if
                    n_cols = n_cols + 1
                    cycle
                end if
                ! Non-underscore lines are data rows
                if (in_atom_loop) then
                    row_vals = ''
                    n_cols = 0
                    i = 1
                    do while (i <= line_len)
                        do while (i <= line_len .and. (line(i:i) == ' ' .or. line(i:i) == char(9)))
                            i = i + 1
                        end do
                        if (i > line_len) exit
                        n_cols = n_cols + 1
                        start_pos = i
                        do while (i <= line_len .and. line(i:i) /= ' ' .and. line(i:i) /= char(9))
                            i = i + 1
                        end do
                        end_pos = i - 1
                        if (n_cols <= 50) row_vals(n_cols) = adjustl(line(start_pos:end_pos))
                    end do
                    if (n_cols >= 4 .and. idx_col(3) > 0 .and. idx_col(4) > 0 .and. idx_col(5) > 0) then
                        n = n + 1
                        if (.not. allocated(data%atoms)) then
                            call grow_atoms(data%atoms, 0, 64)
                        else if (n > size(data%atoms)) then
                            call grow_atoms(data%atoms, n - 1, max(2 * size(data%atoms), n))
                        end if
                        call copy_str_no_quotes(row_vals(idx_col(1)), data%atoms(n)%label)
                        call copy_str_no_quotes(row_vals(idx_col(2)), data%atoms(n)%element)
                        read(row_vals(idx_col(3)), *, iostat=ios) data%atoms(n)%x
                        read(row_vals(idx_col(4)), *, iostat=ios) data%atoms(n)%y
                        read(row_vals(idx_col(5)), *, iostat=ios) data%atoms(n)%z
                    end if
                    cycle
                end if
                ! Non-atom data row - parse as tag-value pair
                tag = ''; val = ''
                i = index(line, ' ')
                if (i > 0) then
                    tag = trim(adjustl(line(1:i-1)))
                    val = clean_str_inline(line(i+1:))
                end if
                cycle
            end if

            ! Outside loop block: all lines are tag-value pairs
            tag = ''; val = ''
            i = index(line, ' ')
            if (i > 0) then
                tag = trim(adjustl(line(1:i-1)))
                val = clean_str_inline(line(i+1:))
            end if

            select case (trim(tag))
            case ('_cell_length_a')    ; read(val, *, iostat=ios) data%a
            case ('_cell_length_b')    ; read(val, *, iostat=ios) data%b
            case ('_cell_length_c')    ; read(val, *, iostat=ios) data%c
            case ('_cell_angle_alpha') ; read(val, *, iostat=ios) data%alpha
            case ('_cell_angle_beta')  ; read(val, *, iostat=ios) data%beta
            case ('_cell_angle_gamma') ; read(val, *, iostat=ios) data%gamma
            case ('_symmetry_space_group_name_H-M', '_space_group_name_H-M_alt')
                data%space_group = adjustl(val(1:min(len_trim(val), len(data%space_group))))
            case ('_symmetry_Int_Tables_number')
                ! Only set from IT number if not already set by an H-M name
                if (trim(data%space_group) == 'P1' .or. trim(data%space_group) == '') then
                    data%space_group = adjustl(val(1:min(len_trim(val), len(data%space_group))))
                end if
            case default
                ! Handle atom_site tags as tag-value pairs (fallback)
                if (index(tag, '_atom_site_label') > 0)       idx_col(1) = n_cols + 1
                if (index(tag, '_atom_site_type_symbol') > 0) idx_col(2) = n_cols + 1
                if (index(tag, '_atom_site_fract_x') > 0)     idx_col(3) = n_cols + 1
                if (index(tag, '_atom_site_fract_y') > 0)     idx_col(4) = n_cols + 1
                if (index(tag, '_atom_site_fract_z') > 0)     idx_col(5) = n_cols + 1
                if (idx_col(3) > 0 .and. idx_col(4) > 0 .and. idx_col(5) > 0) then
                    in_atom_loop = .true.
                end if
                n_cols = n_cols + 1
            end select
        end do

        close(iunit)
        data%n_atoms = n
    end subroutine parse_cif_inline

    subroutine tokenize_inline(line, tokens, n_tokens)
        !! Split whitespace-delimited tokens into array
        character(len=*), intent(in) :: line
        character(len=*), intent(out) :: tokens(:)
        integer, intent(out) :: n_tokens
        integer :: i, n, n_cols, start_pos, end_pos

        n = len_trim(line)
        n_tokens = 0
        n_cols = 0
        i = 1

        do while (i <= n)
            ! Skip whitespace
            do while (i <= n .and. (line(i:i) == ' ' .or. line(i:i) == char(9)))
                i = i + 1
            end do
            if (i > n) exit

            n_cols = n_cols + 1
            start_pos = i
            do while (i <= n .and. line(i:i) /= ' ' .and. line(i:i) /= char(9))
                i = i + 1
            end do
            end_pos = i - 1
            if (n_cols <= size(tokens)) then
                tokens(n_cols) = adjustl(line(start_pos:end_pos))
            end if
        end do
        n_tokens = n_cols
    end subroutine tokenize_inline

    function clean_str_inline(s) result(res)
        character(len=*), intent(in) :: s
        character(len=64) :: res, tmp
        tmp = trim(adjustl(s))
        if (len_trim(tmp) == 0) then; res = ''; return; end if
        if (tmp(1:1) == "'" .or. tmp(1:1) == '"') tmp = tmp(2:)
        if (len_trim(tmp) > 0 .and. tmp(len_trim(tmp):len_trim(tmp)) == "'") &
            tmp = tmp(:len_trim(tmp)-1)
        if (len_trim(tmp) > 0 .and. tmp(len_trim(tmp):len_trim(tmp)) == '"') &
            tmp = tmp(:len_trim(tmp)-1)
        res = trim(adjustl(tmp))
    end function clean_str_inline

    function clean_element_symbol(sym) result(elem)
        !! Strip oxidation state suffix from element symbol (e.g., 'Cu0+' -> 'Cu')
        character(len=*), intent(in) :: sym
        character(len=6) :: elem
        integer :: i, n

        elem = sym
        n = len_trim(sym)
        if (n == 0) return

        ! Element symbol is at most 3 chars: 1 uppercase + (1 or 2 lowercase)
        i = 2
        if (n >= 2 .and. sym(2:2) >= 'a' .and. sym(2:2) <= 'z') i = 3
        if (n >= 3 .and. sym(3:3) >= 'a' .and. sym(3:3) <= 'z') i = 4
        elem = adjustl(sym(1:i-1))
    end function clean_element_symbol

    subroutine copy_str_no_quotes(src, tgt)
        !! Copy src to tgt, stripping leading/trailing quotes and whitespace
        !! tgt must be character(len=*) - copy only meaningful characters
        character(len=*), intent(in)  :: src
        character(len=*), intent(out) :: tgt
        character(len=len(src)) :: tmp
        integer :: i, j, n

        n = len_trim(src)
        if (n == 0) then
            tgt = repeat(' ', len(tgt))
            return
        end if
        tmp = adjustl(src)
        j = 1
        if (tmp(1:1) == "'" .or. tmp(1:1) == '"') tmp = tmp(2:)
        if (n > 0 .and. tmp(len_trim(tmp):len_trim(tmp)) == "'") tmp = tmp(:len_trim(tmp)-1)
        if (n > 0 .and. tmp(len_trim(tmp):len_trim(tmp)) == '"') tmp = tmp(:len_trim(tmp)-1)
        do i = 1, min(len_trim(tmp), len(tgt))
            tgt(i:i) = tmp(i:i)
        end do
    end subroutine copy_str_no_quotes

end module parser
