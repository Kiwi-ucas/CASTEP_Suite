module symmetry
    !! Symmetry operation parsing and asymmetric unit expansion.
    !!
    !! Parses CIF symmetry operation strings (_symmetry_equiv_pos_as_xyz)
    !! of the form "-x, y+1/2, -z+1/2" into 3×3 rotation matrices and
    !! fractional translation vectors, then applies them to generate
    !! the full unit cell from the asymmetric unit.
    use castep_config, only: dp, atom_t, sym_op_t, cif_data_t, pi
    implicit none
    private

    public :: parse_symop_xyz_string, expand_asymmetric_unit, expand_cif_symmetry, free_sym_ops

    ! Maximum number of symmetry operations (space group max = 192 for Fm-3m)
    integer, parameter :: MAX_SYM_OPS = 256

    ! Tolerance for deduplication (fractional coordinates)
    real(dp), parameter :: DUP_TOL = 1.0e-5_dp

contains

    ! ── Public interface ──

    subroutine parse_symop_xyz_string(str, sym_op, iostat, iomsg)
        !! Parse a CIF symmetry operation string like "-x, y+1/2, -z+1/2"
        !! into rotation matrix rot(3,3) and translation vector trans(3).
        !! Rotation entries are 0, ±1, ±2. Translation entries are fractions
        !! like 0, 1/2, 1/3, 2/3, 1/4, 3/4, 1/6, 5/6.
        character(len=*), intent(in)  :: str
        type(sym_op_t),   intent(out) :: sym_op
        integer,          intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        character(len=64) :: parts(3)
        integer :: ios_part(3)

        iostat = 0
        sym_op%rot = 0
        sym_op%trans = 0.0_dp

        ios_part = 0
        call parse_xyz_components(trim(adjustl(str)), parts, iostat, iomsg)
        if (iostat /= 0) return

        call parse_component(parts(1), 'x', 'y', 'z', &
                             sym_op%rot(1,1), sym_op%rot(1,2), sym_op%rot(1,3), &
                             sym_op%trans(1), ios_part(1))
        call parse_component(parts(2), 'x', 'y', 'z', &
                             sym_op%rot(2,1), sym_op%rot(2,2), sym_op%rot(2,3), &
                             sym_op%trans(2), ios_part(2))
        call parse_component(parts(3), 'x', 'y', 'z', &
                             sym_op%rot(3,1), sym_op%rot(3,2), sym_op%rot(3,3), &
                             sym_op%trans(3), ios_part(3))

        if (any(ios_part /= 0)) then
            iostat = 1
            if (present(iomsg)) iomsg = 'Failed to parse symop component: ' // trim(str)
            return
        end if
    end subroutine parse_symop_xyz_string


    subroutine expand_asymmetric_unit(atoms, n_atoms, sym_ops, n_symops, &
                                       expanded_atoms, n_expanded, iostat, iomsg)
        !! Generate the full unit cell by applying each symmetry operation
        !! to each asymmetric unit atom. Deduplicates and sorts the result.
        type(atom_t),     intent(in)  :: atoms(:)
        integer,          intent(in)  :: n_atoms
        type(sym_op_t),   intent(in)  :: sym_ops(:)
        integer,          intent(in)  :: n_symops
        type(atom_t),     allocatable, intent(out) :: expanded_atoms(:)
        integer,          intent(out) :: n_expanded
        integer,          intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        type(atom_t), allocatable :: tmp_atoms(:)
        real(dp) :: frac(3), new_frac(3)
        integer  :: i, j, max_atoms, n_unique
        logical  :: is_dup

        iostat = 0
        n_expanded = 0

        if (n_atoms <= 0) then
            iostat = 1
            if (present(iomsg)) iomsg = 'No atoms to expand'
            return
        end if
        if (n_symops <= 0) then
            iostat = 2
            if (present(iomsg)) iomsg = 'No symmetry operations provided'
            return
        end if

        max_atoms = n_atoms * n_symops
        allocate(tmp_atoms(max_atoms))

        ! Generate all equivalent positions
        n_expanded = 0
        do i = 1, n_atoms
            frac = [atoms(i)%x, atoms(i)%y, atoms(i)%z]
            do j = 1, n_symops
                new_frac(1) = sym_ops(j)%rot(1,1)*frac(1) &
                            + sym_ops(j)%rot(1,2)*frac(2) &
                            + sym_ops(j)%rot(1,3)*frac(3) &
                            + sym_ops(j)%trans(1)
                new_frac(2) = sym_ops(j)%rot(2,1)*frac(1) &
                            + sym_ops(j)%rot(2,2)*frac(2) &
                            + sym_ops(j)%rot(2,3)*frac(3) &
                            + sym_ops(j)%trans(2)
                new_frac(3) = sym_ops(j)%rot(3,1)*frac(1) &
                            + sym_ops(j)%rot(3,2)*frac(2) &
                            + sym_ops(j)%rot(3,3)*frac(3) &
                            + sym_ops(j)%trans(3)

                ! Wrap into [0, 1)
                call wrap_to_unit(new_frac(1))
                call wrap_to_unit(new_frac(2))
                call wrap_to_unit(new_frac(3))

                n_expanded = n_expanded + 1
                tmp_atoms(n_expanded)%label   = atoms(i)%label
                tmp_atoms(n_expanded)%element = atoms(i)%element
                tmp_atoms(n_expanded)%x = new_frac(1)
                tmp_atoms(n_expanded)%y = new_frac(2)
                tmp_atoms(n_expanded)%z = new_frac(3)
            end do
        end do

        ! Deduplicate
        n_unique = 0
        do i = 1, n_expanded
            is_dup = .false.
            do j = 1, n_unique
                if (same_position(tmp_atoms(i), tmp_atoms(j), DUP_TOL)) then
                    is_dup = .true.
                    exit
                end if
            end do
            if (.not. is_dup) then
                n_unique = n_unique + 1
                if (n_unique /= i) tmp_atoms(n_unique) = tmp_atoms(i)
            end if
        end do

        ! Sort by z, then y, then x
        call sort_atoms(tmp_atoms, n_unique)

        allocate(expanded_atoms(n_unique))
        expanded_atoms(1:n_unique) = tmp_atoms(1:n_unique)
        n_expanded = n_unique
        deallocate(tmp_atoms)

    end subroutine expand_asymmetric_unit


    subroutine expand_cif_symmetry(cif, iostat)
        !! Convenience wrapper: expand cif%atoms using cif%sym_ops if n_symops > 1.
        !! On success, cif%atoms is replaced with the expanded full-cell atoms.
        !! If n_symops <= 1, this is a no-op (P1 structure).
        type(cif_data_t), intent(inout) :: cif
        integer,          intent(out)   :: iostat

        type(atom_t), allocatable :: expanded(:)
        integer :: n_expanded

        iostat = 0
        if (cif%n_symops <= 1) return

        call expand_asymmetric_unit(cif%atoms, cif%n_atoms, &
                                     cif%sym_ops, cif%n_symops, &
                                     expanded, n_expanded, iostat)
        if (iostat == 0) then
            if (allocated(cif%atoms)) deallocate(cif%atoms)
            allocate(cif%atoms(n_expanded))
            cif%atoms = expanded
            cif%n_atoms = n_expanded
            deallocate(expanded)
        end if
    end subroutine expand_cif_symmetry


    subroutine free_sym_ops(data_symops, data_n_symops)
        !! Deallocate symmetry operations array and reset count.
        type(sym_op_t), allocatable, intent(inout) :: data_symops(:)
        integer, intent(inout) :: data_n_symops
        if (allocated(data_symops)) deallocate(data_symops)
        data_n_symops = 0
    end subroutine free_sym_ops


    ! ── Private helpers ──

    subroutine parse_xyz_components(str, parts, iostat, iomsg)
        !! Split "x-component, y-component, z-component" into 3 parts.
        !! Handles commas that may appear inside quoted strings (unlikely but robust).
        character(len=*), intent(in)  :: str
        character(len=*), intent(out) :: parts(3)
        integer,          intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        integer :: comma_pos(2), n

        iostat = 0
        n = len_trim(str)

        ! Find the two comma positions that separate the 3 components.
        ! We scan and ignore commas inside single quotes (rare for symop strings).
        comma_pos = 0
        call find_separating_commas(str, n, comma_pos(1), comma_pos(2))

        if (comma_pos(1) <= 1 .or. comma_pos(2) <= comma_pos(1) .or. comma_pos(2) >= n) then
            iostat = 1
            if (present(iomsg)) iomsg = 'Cannot split symop into 3 components: ' // trim(str)
            return
        end if

        ! 1st component: start to comma(1)-1
        parts(1) = trim(adjustl(str(1:comma_pos(1)-1)))
        ! 2nd component: comma(1)+1 to comma(2)-1
        parts(2) = trim(adjustl(str(comma_pos(1)+1:comma_pos(2)-1)))
        ! 3rd component: comma(2)+1 to end
        parts(3) = trim(adjustl(str(comma_pos(2)+1:n)))
    end subroutine parse_xyz_components


    subroutine find_separating_commas(str, n, c1, c2)
        !! Find the two comma positions that separate the 3 xyz components.
        character(len=*), intent(in) :: str
        integer, intent(in) :: n
        integer, intent(out) :: c1, c2
        integer :: i, depth
        logical :: in_squote, in_dquote
        in_squote = .false.
        in_dquote = .false.
        depth = 0
        c1 = 0
        c2 = 0
        do i = 1, n
            if (str(i:i) == "'" .and. .not. in_dquote) in_squote = .not. in_squote
            if (str(i:i) == '"' .and. .not. in_squote) in_dquote = .not. in_dquote
            if (in_squote .or. in_dquote) cycle
            if (str(i:i) == ',') then
                if (c1 == 0) then
                    c1 = i
                else if (c2 == 0) then
                    c2 = i
                    return
                end if
            end if
        end do
    end subroutine find_separating_commas


    subroutine parse_component(comp, xvar, yvar, zvar, rx, ry, rz, t, iostat)
        !! Parse a single component of a symmetry operation like "-x+1/2" or "y".
        !! Sets exactly one of rx, ry, rz to ±1 (or ±2) and sets t to the
        !! translation fraction (0 if none).
        character(len=*), intent(in)  :: comp
        character(len=*), intent(in)  :: xvar, yvar, zvar
        integer,          intent(out) :: rx, ry, rz
        real(dp),         intent(out) :: t
        integer,          intent(out) :: iostat

        character(len=64) :: s
        integer :: i, n, sign, coeff
        logical :: found_var

        iostat = 0
        rx = 0; ry = 0; rz = 0
        t = 0.0_dp
        s = trim(adjustl(comp))
        n = len_trim(s)
        if (n == 0) then
            iostat = 1; return
        end if

        ! Determine overall sign
        sign = 1
        if (s(1:1) == '-') then
            sign = -1
            s = s(2:)
            n = n - 1
        else if (s(1:1) == '+') then
            s = s(2:)
            n = n - 1
        end if

        ! Scan for variable token (x, y, or z) with optional numeric prefix
        found_var = .false.
        i = 1
        do while (i <= n)
            ! Check for variable at current position
            if (s(i:i) == 'x' .or. s(i:i) == 'y' .or. s(i:i) == 'z') then
                if (.not. found_var) then
                    ! Determine coefficient
                    if (i > 1 .and. is_digit(s(i-1:i-1))) then
                        read(s(i-1:i-1), '(i1)') coeff
                    else
                        coeff = 1
                    end if
                    coeff = sign * coeff
                    select case (s(i:i))
                    case ('x'); rx = coeff
                    case ('y'); ry = coeff
                    case ('z'); rz = coeff
                    end select
                    found_var = .true.
                end if
            end if
            i = i + 1
        end do

        ! Parse the translation part (everything not a variable or its coefficient)
        t = parse_translation_fraction(s)
    end subroutine parse_component


    real(dp) function parse_translation_fraction(s) result(t)
        !! Extract fractional translation like "1/2", "1/3", "2/3", "", "1" from
        !! the non-variable part of a symop component. Handles sign.
        character(len=*), intent(in) :: s

        integer :: i, n, num, den, sign, slash_pos, num_end, num_start
        character(len=64) :: num_str

        t = 0.0_dp
        n = len_trim(s)
        if (n == 0) return

        ! Scan for a fraction pattern: digits/digits, possibly signed
        ! Look for '/' in the string (excluding positions adjacent to variables)
        slash_pos = 0
        do i = 1, n
            if (s(i:i) == '/') then
                slash_pos = i
                exit
            end if
        end do

        if (slash_pos > 0) then
            ! Fraction found: extract numerator and denominator
            ! Find start of numerator (first digit before slash)
            num_start = slash_pos - 1
            do while (num_start >= 1)
                if (is_digit(s(num_start:num_start))) then
                    num_start = num_start - 1
                else
                    exit
                end if
            end do
            num_start = num_start + 1

            ! Check for sign before numerator
            sign = 1
            if (num_start > 1) then
                if (s(num_start-1:num_start-1) == '-') sign = -1
            end if

            num_str = s(num_start:slash_pos-1)
            read(num_str, '(i10)', err=99) num

            num_str = s(slash_pos+1:n)
            ! Only take consecutive digits after slash
            num_end = 1
            do while (num_end <= len_trim(num_str))
                if (.not. is_digit(num_str(num_end:num_end))) exit
                num_end = num_end + 1
            end do
            num_end = num_end - 1
            if (num_end < 1) goto 99
            read(num_str(1:num_end), '(i10)', err=99) den

            if (den > 0) then
                t = real(sign * num, dp) / real(den, dp)
            end if
            return
        end if

        ! No fraction — look for a bare integer at the end
        ! A bare "+1" or "-1" without a variable preceding would be a translation
        ! But because variables also consume their sign, a bare translation looks
        ! like "+1" or "-1" where the sign is part of the overall sign consumed above.
        ! However the parse_component already strips the leading sign, so a bare
        ! "1" would remain. Let's look for isolated digit at the very end.
        do i = n, 1, -1
            if (is_digit(s(i:i))) then
                ! Check this digit is not adjacent to a variable
                if (i < n) then  ! digit not at end… skip
                else
                    read(s(i:i), '(i1)') num
                    t = real(num, dp)
                end if
                exit
            end if
        end do
        return

99      continue
        t = 0.0_dp
    end function parse_translation_fraction


    pure logical function is_digit(ch)
        character(len=1), intent(in) :: ch
        is_digit = (ch >= '0' .and. ch <= '9')
    end function is_digit


    pure subroutine wrap_to_unit(x)
        !! Wrap fractional coordinate x into [0, 1).
        real(dp), intent(inout) :: x
        real(dp), parameter :: ONE = 1.0_dp, ZERO = 0.0_dp
        ! Floor to integer → subtract to get [0,1)
        x = x - floor(x)
        if (x < ZERO) x = x + ONE
        if (x >= ONE .or. x < ZERO) x = x - aint(x)  ! safety
    end subroutine wrap_to_unit


    pure logical function same_position(a, b, tol)
        !! Compare two atom fractional positions within tolerance.
        type(atom_t), intent(in) :: a, b
        real(dp),     intent(in) :: tol
        same_position = abs(a%x - b%x) < tol .and. &
                        abs(a%y - b%y) < tol .and. &
                        abs(a%z - b%z) < tol
    end function same_position


    subroutine sort_atoms(atoms, n)
        !! Sort atoms by z, then y, then x (stable-ish insertion sort).
        type(atom_t), intent(inout) :: atoms(:)
        integer,      intent(in)    :: n
        type(atom_t) :: tmp
        integer :: i, j
        do i = 2, n
            tmp = atoms(i)
            j = i - 1
            do while (j >= 1)
                if (atoms(j)%z < tmp%z - 1.0e-10_dp) exit
                if (abs(atoms(j)%z - tmp%z) < 1.0e-10_dp) then
                    if (atoms(j)%y < tmp%y - 1.0e-10_dp) exit
                    if (abs(atoms(j)%y - tmp%y) < 1.0e-10_dp) then
                        if (atoms(j)%x <= tmp%x + 1.0e-10_dp) exit
                    end if
                end if
                atoms(j+1) = atoms(j)
                j = j - 1
            end do
            atoms(j+1) = tmp
        end do
    end subroutine sort_atoms

end module symmetry
