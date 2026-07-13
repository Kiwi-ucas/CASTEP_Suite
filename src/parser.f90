module parser
    !! File format parsers for CIF, PDB, and CASTEP .cell files
    use castep_config, only: &
         dp, MAX_ATOMS, IO_FILE_NOT_FOUND, IO_PARSE_ERROR, IO_MISSING_ATOMS, &
         cif_data_t, atom_t, pi, compare_tags, sym_op_t
    use symmetry, only: parse_symop_xyz_string
    implicit none
    private

    public :: parse_cif_inline
    public :: parse_pdb_inline
    public :: parse_cell_inline
    public :: clean_element_symbol

contains

    subroutine grow_atoms(atoms, current, new_size, iostat)
        !! Grow the atoms array to new_size, preserving all existing entries.
        !! Returns iostat=0 on success, iostat=1 on allocation failure.
        type(atom_t), allocatable, intent(inout) :: atoms(:)
        integer, intent(in) :: current, new_size
        integer, intent(out) :: iostat
        type(atom_t), allocatable :: tmp(:)
        integer :: ios

        if (.not. allocated(atoms)) then
            allocate(atoms(new_size), stat=ios)
        else
            allocate(tmp(new_size), stat=ios)
            if (ios == 0) then
                tmp(1:current) = atoms(1:current)
                call move_alloc(tmp, atoms)
            end if
        end if
        if (ios /= 0) then
            iostat = 1
        else
            iostat = 0
        end if
    end subroutine grow_atoms

    subroutine grow_symops(symops, current, new_size, iostat)
        !! Grow the sym_ops array to new_size, preserving all existing entries.
        !! Returns iostat=0 on success, iostat=1 on allocation failure.
        type(sym_op_t), allocatable, intent(inout) :: symops(:)
        integer, intent(in) :: current, new_size
        integer, intent(out) :: iostat
        type(sym_op_t), allocatable :: tmp(:)
        integer :: ios

        if (.not. allocated(symops)) then
            allocate(symops(new_size), stat=ios)
        else
            allocate(tmp(new_size), stat=ios)
            if (ios == 0) then
                tmp(1:current) = symops(1:current)
                call move_alloc(tmp, symops)
            end if
        end if
        if (ios /= 0) then
            iostat = 1
        else
            iostat = 0
        end if
    end subroutine grow_symops

    ! ── Private CIF helpers ──

    subroutine strip_uncertainty(val_str, clean_str)
        !! Strip parenthetical standard uncertainty from numeric string.
        !! "8.6559(9)" → "8.6559", "0.71073" → "0.71073"
        character(len=*), intent(in)  :: val_str
        character(len=*), intent(out) :: clean_str
        integer :: paren_pos
        clean_str = trim(adjustl(val_str))
        paren_pos = index(clean_str, '(')
        if (paren_pos > 0) then
            clean_str = trim(clean_str(1:paren_pos-1))
        end if
    end subroutine strip_uncertainty


    logical function is_textfield_start(line)
        !! Return .true. if the trimmed line starts with ';' (semicolon),
        !! indicating the start of a CIF TextField.
        character(len=*), intent(in) :: line
        character(len=256) :: t
        t = adjustl(line)
        is_textfield_start = (len_trim(t) > 0 .and. t(1:1) == ';')
    end function is_textfield_start


    subroutine parse_textfield_body(iunit, body, ios)
        !! Read TextField body from current position until closing ';'.
        !! The opening ';' has already been consumed from the first line.
        !! Reads lines until a line with only ';' (or ';' + whitespace) is encountered.
        !! Returns the concatenated text (lines joined by newline).
        integer,           intent(in)  :: iunit
        character(len=:),  allocatable, intent(out) :: body
        integer,           intent(out) :: ios

        character(len=256) :: line, trimmed
        character(len=4096) :: accumulator
        integer :: body_len, line_len

        body_len = 0
        accumulator = ''
        ios = 0

        do
            read(iunit, '(A)', iostat=ios) line
            if (ios /= 0) return  ! EOF or error

            trimmed = adjustl(line)
            line_len = len_trim(trimmed)

            ! Check for closing semicolon (a line with only ';', possibly with whitespace)
            if (line_len >= 1 .and. trimmed(1:1) == ';') then
                ! If the rest of the line is whitespace-only, this is the closing delimiter
                if (line_len == 1) exit
                if (len_trim(trimmed(2:)) == 0) exit
                ! Otherwise it contains content — likely a line starting with ';'
                ! In the CIF standard, a line starting with ';' inside a textfield
                ! is supposed to be the closing delimiter, so we exit
                exit
            end if

            ! Append this line to accumulator
            if (body_len + line_len + 1 <= len(accumulator)) then
                if (body_len > 0) then
                    accumulator(body_len+1:body_len+1) = char(10)
                    body_len = body_len + 1
                end if
                if (line_len > 0) then
                    accumulator(body_len+1:body_len+line_len) = trim(line)
                else
                    accumulator(body_len+1:body_len+1) = ' '
                end if
                body_len = body_len + max(line_len, 1)
            end if
        end do

        allocate(character(len=body_len) :: body)
        if (body_len > 0) then
            body = accumulator(1:body_len)
        else
            body = ''
        end if
    end subroutine parse_textfield_body

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
            if (line(1:4) == 'ATOM' .or. line(1:6) == 'HETATM') then
                n_atoms = n_atoms + 1
                if (n_atoms > MAX_ATOMS) exit

                if (.not. allocated(data%atoms)) then
                    call grow_atoms(data%atoms, 0, 64, ios)
                    if (ios /= 0) exit
                else if (n_atoms > size(data%atoms)) then
                    call grow_atoms(data%atoms, n_atoms - 1, max(2 * size(data%atoms), n_atoms), ios)
                    if (ios /= 0) exit
                end if

                ! Element symbol: columns 77-78 (right-justified)
                data%atoms(n_atoms)%element = line(77:78)
                data%atoms(n_atoms)%element = adjustl(data%atoms(n_atoms)%element)

                ! Atom name: columns 13-16
                data%atoms(n_atoms)%label = adjustl(line(13:16))

                ! Cartesian coordinates: columns 31-38, 39-46, 47-54
                read(line(31:38), '(F8.3)', iostat=ios) data%atoms(n_atoms)%x
                if (ios /= 0) then; n_atoms = n_atoms - 1; cycle; end if
                read(line(39:46), '(F8.3)', iostat=ios) data%atoms(n_atoms)%y
                if (ios /= 0) then; n_atoms = n_atoms - 1; cycle; end if
                read(line(47:54), '(F8.3)', iostat=ios) data%atoms(n_atoms)%z
                if (ios /= 0) then; n_atoms = n_atoms - 1; cycle; end if
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
                if (compare_tags(trimmed, '%BLOCK POSITIONS_FRAC')) &
                    data%positions_fractional = .true.
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
                    read(trimmed, *, iostat=ios_local) data%a, data%b, data%c
                    if (ios_local == 0) read_a_line = .true.
                else if (.not. read_angle_line) then
                    read(trimmed, *, iostat=ios_local) data%alpha, data%beta, data%gamma
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
                    lattice_cart(1, cart_line) = x
                    lattice_cart(2, cart_line) = y
                    lattice_cart(3, cart_line) = z
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
                            call grow_atoms(data%atoms, 0, 64, ios_local)
                            if (ios_local /= 0) exit
                        else if (n_atoms > size(data%atoms)) then
                            call grow_atoms(data%atoms, n_atoms - 1, max(2 * size(data%atoms), n_atoms), ios_local)
                            if (ios_local /= 0) exit
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

        integer :: ios, iunit, i, n, n_cols, line_len, start_pos, end_pos, j
        integer :: idx_col(5), idx_symop_col
        character(len=256) :: line, tag, val, clean_val
        logical :: in_atom_loop, in_symop_loop, in_loop, file_found
        character(len=256) :: row_vals(50)
        character(len=:), allocatable :: text_body
        type(sym_op_t) :: sop
        ! Warning counters for parse failures
        integer :: n_warn_cell, n_warn_atom, n_warn_symop, n_warn_text
        character(len=256) :: accum_warnings

        n_warn_cell = 0; n_warn_atom = 0; n_warn_symop = 0; n_warn_text = 0
        accum_warnings = ''
        iostat = 0
        data%a = 0.0_dp; data%b = 0.0_dp; data%c = 0.0_dp
        data%alpha = 90.0_dp; data%beta = 90.0_dp; data%gamma = 90.0_dp
        data%space_group = 'P1'
        data%space_group_name = ''
        data%space_group_it = 0
        data%n_atoms = 0
        data%n_symops = 0
        data%positions_fractional = .true.   ! CIF coordinates are always fractional
        idx_col(1:5) = 0
        idx_symop_col = 0
        if (allocated(data%atoms)) deallocate(data%atoms)
        if (allocated(data%sym_ops)) deallocate(data%sym_ops)

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
        in_symop_loop = .false.
        in_loop = .false.
        n = 0
        n_cols = 0
        idx_symop_col = 0

        do
            read(iunit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = trim(adjustl(line))
            line_len = len_trim(line)
            if (line_len == 0 .or. line(1:1) == '#') cycle

            if (trim(line) == 'loop_') then
                in_atom_loop = .false.
                in_symop_loop = .false.
                in_loop = .true.
                n_cols = 0
                idx_col(1:5) = 0
                idx_symop_col = 0
                cycle
            end if

            ! Inside a loop block
            if (in_loop) then
                ! Underscore lines may be loop headers or tag-value pairs (end of loop)
                if (line_len >= 1 .and. line(1:1) == '_') then
                    ! Check if this is a tag-value pair (has a value after the tag)
                    ! If so, exit loop mode and process as a regular tag-value
                    j = index(line, ' ')
                    if (j > 0) then
                        val = clean_str_inline(line(j+1:))
                        if (len_trim(val) > 0) then
                            in_loop = .false.
                            in_atom_loop = .false.
                            goto 100
                        end if
                    end if
                    ! Detect atom_site columns using exact tag matching
                    j = index(line, ' ')
                    if (j > 0) then
                        tag = adjustl(line(1:j-1))
                    else
                        tag = adjustl(line)
                    end if
                    if (compare_tags(tag, '_atom_site_label') &
                        .or. compare_tags(tag(1:18), '_atom_site_label')) &
                        idx_col(1) = n_cols + 1
                    if (compare_tags(tag, '_atom_site_type_symbol')) &
                        idx_col(2) = n_cols + 1
                    if (compare_tags(tag, '_atom_site_fract_x')) &
                        idx_col(3) = n_cols + 1
                    if (compare_tags(tag, '_atom_site_fract_y')) &
                        idx_col(4) = n_cols + 1
                    if (compare_tags(tag, '_atom_site_fract_z')) &
                        idx_col(5) = n_cols + 1
                    if (idx_col(3) > 0 .and. idx_col(4) > 0 .and. idx_col(5) > 0) then
                        in_atom_loop = .true.
                    end if
                    ! Detect symmetry operation column
                    if (compare_tags(tag, '_symmetry_equiv_pos_as_xyz') &
                        .or. compare_tags(tag, '_space_group_symop_operation_xyz')) then
                        idx_symop_col = n_cols + 1
                        in_symop_loop = .true.
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
                        if (n >= MAX_ATOMS) cycle  ! hard limit
                        n = n + 1
                        if (.not. allocated(data%atoms)) then
                            call grow_atoms(data%atoms, 0, 64, ios)
                            if (ios /= 0) exit
                        else if (n > size(data%atoms)) then
                            call grow_atoms(data%atoms, n - 1, max(2 * size(data%atoms), n), ios)
                            if (ios /= 0) exit
                        end if
                        call copy_str_no_quotes(row_vals(idx_col(1)), data%atoms(n)%label)
                        call copy_str_no_quotes(row_vals(idx_col(2)), data%atoms(n)%element)
                        call strip_uncertainty(row_vals(idx_col(3)), clean_val)
                        read(clean_val, *, iostat=ios) data%atoms(n)%x
                        if (ios /= 0) then; n = n - 1; n_warn_atom = n_warn_atom + 1; cycle; end if
                        call strip_uncertainty(row_vals(idx_col(4)), clean_val)
                        read(clean_val, *, iostat=ios) data%atoms(n)%y
                        if (ios /= 0) then; n = n - 1; n_warn_atom = n_warn_atom + 1; cycle; end if
                        call strip_uncertainty(row_vals(idx_col(5)), clean_val)
                        read(clean_val, *, iostat=ios) data%atoms(n)%z
                        if (ios /= 0) then; n = n - 1; n_warn_atom = n_warn_atom + 1; cycle; end if
                    end if
                    cycle
                end if
                ! Symmetry operation data row — values may be quoted
                if (in_symop_loop) then
                    row_vals = ''
                    n_cols = 0
                    call tokenize_quoted_inline(line, row_vals, n_cols)
                    if (idx_symop_col > 0 .and. idx_symop_col <= n_cols) then
                        if (data%n_symops >= 256) cycle  ! safety cap
                        if (.not. allocated(data%sym_ops)) then
                            allocate(data%sym_ops(64))
                        else if (data%n_symops >= size(data%sym_ops)) then
                            call grow_symops(data%sym_ops, data%n_symops, &
                                             2 * size(data%sym_ops), ios)
                            if (ios /= 0) cycle
                        end if
                        call copy_str_no_quotes(row_vals(idx_symop_col), val)
                        call parse_symop_xyz_string(val, sop, ios)
                        if (ios == 0) then
                            data%n_symops = data%n_symops + 1
                            data%sym_ops(data%n_symops) = sop
                        else
                            n_warn_symop = n_warn_symop + 1
                        end if
                    end if
                    cycle
                end if

                ! Non-atom data row - skip
                cycle
            end if

            ! Outside loop block: all lines are tag-value pairs
100         continue
            tag = ''; val = ''
            i = index(line, ' ')
            if (i > 0) then
                tag = trim(adjustl(line(1:i-1)))
                val = clean_str_inline(line(i+1:))
            end if

            ! Handle TextField: if value starts with ';', read multiline body
            if (is_textfield_start(line(i+1:))) then
                call parse_textfield_body(iunit, text_body, ios)
                if (ios == 0) then
                    val = text_body
                    deallocate(text_body)
                else
                    n_warn_text = n_warn_text + 1
                end if
            end if

            ! Strip uncertainty before numeric parsing
            call strip_uncertainty(val, clean_val)

            select case (trim(tag))
            case ('_cell_length_a')    ; read(clean_val, *, iostat=ios) data%a
                                         if (ios /= 0) n_warn_cell = n_warn_cell + 1
            case ('_cell_length_b')    ; read(clean_val, *, iostat=ios) data%b
                                         if (ios /= 0) n_warn_cell = n_warn_cell + 1
            case ('_cell_length_c')    ; read(clean_val, *, iostat=ios) data%c
                                         if (ios /= 0) n_warn_cell = n_warn_cell + 1
            case ('_cell_angle_alpha') ; read(clean_val, *, iostat=ios) data%alpha
                                         if (ios /= 0) n_warn_cell = n_warn_cell + 1
            case ('_cell_angle_beta')  ; read(clean_val, *, iostat=ios) data%beta
                                         if (ios /= 0) n_warn_cell = n_warn_cell + 1
            case ('_cell_angle_gamma') ; read(clean_val, *, iostat=ios) data%gamma
                                         if (ios /= 0) n_warn_cell = n_warn_cell + 1
            case ('_symmetry_space_group_name_H-M', '_space_group_name_H-M_alt', &
                  '_symmetry_space_group_name_H-M_alt')
                data%space_group_name = adjustl(val(1:min(len_trim(val), len(data%space_group_name))))
                if (trim(data%space_group_name) /= '') &
                    data%space_group = adjustl(val(1:min(len_trim(val), len(data%space_group))))
            case ('_symmetry_Int_Tables_number', '_space_group_IT_number')
                ! Store IT number for symmetry expansion
                read(clean_val, *, iostat=ios) data%space_group_it
                if (ios /= 0) data%space_group_it = 0  ! reset on parse failure
                ! Only set space_group from IT number if H-M name not already set
                if (trim(data%space_group) == 'P1' .or. trim(data%space_group) == '') then
                    data%space_group = adjustl(clean_val(1:min(len_trim(clean_val), len(data%space_group))))
                end if
            case default
                ! Handle atom_site tags as tag-value pairs (fallback)
                if (compare_tags(tag, '_atom_site_label') &
                    .or. (len_trim(tag) >= 17 .and. compare_tags(tag(1:17), '_atom_site_label'))) &
                    idx_col(1) = n_cols + 1
                if (compare_tags(tag, '_atom_site_type_symbol')) &
                    idx_col(2) = n_cols + 1
                if (compare_tags(tag, '_atom_site_fract_x')) &
                    idx_col(3) = n_cols + 1
                if (compare_tags(tag, '_atom_site_fract_y')) &
                    idx_col(4) = n_cols + 1
                if (compare_tags(tag, '_atom_site_fract_z')) &
                    idx_col(5) = n_cols + 1
                if (idx_col(3) > 0 .and. idx_col(4) > 0 .and. idx_col(5) > 0) then
                    in_atom_loop = .true.
                end if
                n_cols = n_cols + 1
            end select
        end do

        close(iunit)
        data%n_atoms = n

        ! ── Post-parse validation and warnings ──
        ! Hard errors (set iostat)
        if (n_warn_cell > 0) then
            iostat = IO_PARSE_ERROR
            if (present(iomsg)) write(iomsg, &
                '("Failed to read ", i0, " cell parameter(s)")') n_warn_cell
        end if
        if (n == 0 .and. iostat == 0) then
            iostat = IO_MISSING_ATOMS
            if (present(iomsg)) iomsg = 'No atoms found in CIF file'
        end if
        ! Soft warnings (print only, don't block)
        if (n_warn_atom > 0) &
            write(*, '(a,i0,a)') '  Warning: Skipped ', n_warn_atom, &
                ' unreadable atom coordinate(s)'
        if (n_warn_symop > 0) &
            write(*, '(a,i0,a)') '  Warning: Skipped ', n_warn_symop, &
                ' unparseable symmetry operation(s)'
        if (n_warn_text > 0) &
            write(*, '(a,i0,a)') '  Warning: ', n_warn_text, &
                ' TextField(s) could not be read'
    end subroutine parse_cif_inline

    subroutine tokenize_quoted_inline(line, tokens, n_tokens)
        !! Split tokens preserving quoted strings (single and double quotes).
        !! e.g., `1  'x, y, z'` → ["1", "x, y, z"]
        character(len=*), intent(in) :: line
        character(len=*), intent(out) :: tokens(:)
        integer, intent(out) :: n_tokens
        integer :: i, n, n_cols, start_pos, end_pos
        character(len=1) :: quote_char

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
            ! Check for quoted token
            if (line(i:i) == "'" .or. line(i:i) == '"') then
                quote_char = line(i:i)
                i = i + 1  ! skip opening quote
                start_pos = i
                do while (i <= n .and. line(i:i) /= quote_char)
                    i = i + 1
                end do
                end_pos = i - 1
                if (i <= n) i = i + 1  ! skip closing quote
                if (n_cols <= size(tokens)) then
                    if (end_pos >= start_pos) then
                        tokens(n_cols) = adjustl(line(start_pos:end_pos))
                    else
                        tokens(n_cols) = ''
                    end if
                end if
            else
                ! Unquoted token
                start_pos = i
                do while (i <= n .and. line(i:i) /= ' ' .and. line(i:i) /= char(9))
                    i = i + 1
                end do
                end_pos = i - 1
                if (n_cols <= size(tokens)) then
                    tokens(n_cols) = adjustl(line(start_pos:end_pos))
                end if
            end if
        end do
        n_tokens = n_cols
    end subroutine tokenize_quoted_inline

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
        !! Strip leading/trailing quotes and unescape CIF doubled quotes.
        !! 'it''s'   → it's
        !! "say ""hi""" → say "hi"
        character(len=*), intent(in) :: s
        character(len=64) :: res, tmp
        integer :: i, n, j
        character(len=1) :: quote_char
        tmp = trim(adjustl(s))
        n = len_trim(tmp)
        if (n == 0) then; res = ''; return; end if
        ! Detect quote type
        quote_char = ' '
        if (tmp(1:1) == "'" .or. tmp(1:1) == '"') quote_char = tmp(1:1)
        if (quote_char /= ' ') then
            ! Strip outer quotes
            tmp = tmp(2:n)
            n = n - 1
            if (n > 0 .and. tmp(n:n) == quote_char) then
                tmp(n:n) = ' '; n = n - 1
            end if
            ! Unescape doubled inner quotes: '' → '  or  "" → "
            j = 1; i = 1
            do while (i <= n)
                if (i < n .and. tmp(i:i) == quote_char .and. tmp(i+1:i+1) == quote_char) then
                    tmp(j:j) = quote_char
                    j = j + 1
                    i = i + 2  ! skip both quotes
                else
                    tmp(j:j) = tmp(i:i)
                    j = j + 1
                    i = i + 1
                end if
            end do
            n = j - 1
            tmp(n+1:) = ' '
        end if
        res = trim(adjustl(tmp(1:n)))
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
        !! Copy src to tgt, stripping leading/trailing quotes and unescaping
        !! CIF doubled quotes ('' → ', "" → ").
        !! tgt must be character(len=*) - copy only meaningful characters
        character(len=*), intent(in)  :: src
        character(len=*), intent(out) :: tgt
        character(len=len(src)) :: tmp
        character(len=1) :: quote_char
        integer :: i, n, j

        n = len_trim(src)
        if (n == 0) then
            tgt = repeat(' ', len(tgt))
            return
        end if
        tmp = adjustl(src)
        n = len_trim(tmp)

        ! Detect quote type
        quote_char = ' '
        if (tmp(1:1) == "'" .or. tmp(1:1) == '"') quote_char = tmp(1:1)

        if (quote_char /= ' ') then
            ! Strip outer quotes
            tmp = tmp(2:n)
            n = n - 1
            if (n > 0 .and. tmp(n:n) == quote_char) then
                tmp(n:n) = ' '; n = n - 1
            end if
            ! Unescape doubled inner quotes: '' → '  or  "" → "
            j = 1; i = 1
            do while (i <= n)
                if (i < n .and. tmp(i:i) == quote_char .and. tmp(i+1:i+1) == quote_char) then
                    tmp(j:j) = quote_char
                    j = j + 1
                    i = i + 2
                else
                    tmp(j:j) = tmp(i:i)
                    j = j + 1
                    i = i + 1
                end if
            end do
            n = j - 1
        end if

        n = min(n, len(tgt))
        do i = 1, n
            tgt(i:i) = tmp(i:i)
        end do
        do i = n + 1, len(tgt)
            tgt(i:i) = ' '
        end do
    end subroutine copy_str_no_quotes

end module parser
