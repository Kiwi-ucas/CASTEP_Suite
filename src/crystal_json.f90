module crystal_json
    !! Generate crystal_data.json for the Rust Crystal Viewer
    use castep_config, only: dp, pi, cif_data_t, IO_PARSE_ERROR
    use parser, only: clean_element_symbol
    implicit none
    private

    public :: write_crystal_json_cif, read_crystal_json_to_cif
    public :: write_crystal_json_modes

contains

    subroutine write_crystal_json_cif(cif, json_path, iostat)
        !! Write crystal structure from cif_data_t (auto-converts fractional→Cartesian)
        type(cif_data_t), intent(in) :: cif
        character(len=*), intent(in) :: json_path
        integer, intent(out) :: iostat
        real(dp) :: va(3), vb(3), vc(3), x, y, z
        integer :: unit, i

        iostat = 0
        open(newunit=unit, file=trim(json_path), status='replace', action='write', iostat=iostat)
        if (iostat /= 0) return

        ! Compute lattice vectors
        call lattice_vectors(cif%a, cif%b, cif%c, cif%alpha, cif%beta, cif%gamma, va, vb, vc)

        write(unit, '(a)') '{'
        write(unit, '(a)') '  "lattice": {'
        write(unit, '(a,f10.4,a)') '    "a": ', cif%a, ','
        write(unit, '(a,f10.4,a)') '    "b": ', cif%b, ','
        write(unit, '(a,f10.4,a)') '    "c": ', cif%c, ','
        write(unit, '(a,f10.4,a)') '    "alpha": ', cif%alpha, ','
        write(unit, '(a,f10.4,a)') '    "beta": ', cif%beta, ','
        write(unit, '(a,f10.4)') '    "gamma": ', cif%gamma
        write(unit, '(a)') '  },'
        write(unit, '(a)') '  "positions_fractional": false,'
        write(unit, '(a)') '  "atoms": ['

        do i = 1, cif%n_atoms
            if (cif%positions_fractional) then
                x = cif%atoms(i)%x * va(1) + cif%atoms(i)%y * vb(1) + cif%atoms(i)%z * vc(1)
                y = cif%atoms(i)%x * va(2) + cif%atoms(i)%y * vb(2) + cif%atoms(i)%z * vc(2)
                z = cif%atoms(i)%x * va(3) + cif%atoms(i)%y * vb(3) + cif%atoms(i)%z * vc(3)
            else
                x = cif%atoms(i)%x
                y = cif%atoms(i)%y
                z = cif%atoms(i)%z
            end if
            write(unit, '(a,a,a,a,f10.6,a,a,f10.6,a,a,f10.6,a)', advance='no') &
                '    { "element": "', trim(clean_element_symbol(cif%atoms(i)%element)), '", ', &
                '"x": ', x, ', ', &
                '"y": ', y, ', ', &
                '"z": ', z, ' }'
            if (i < cif%n_atoms) write(unit, '(a)') ','
            if (i == cif%n_atoms) write(unit, '(a)') ''
        end do

        write(unit, '(a)') '  ]'
        write(unit, '(a)') '}'
        close(unit)
    end subroutine write_crystal_json_cif


    pure subroutine lattice_vectors(a, b, c, alpha, beta, gamma, va, vb, vc)
        !! Compute Cartesian lattice vectors from cell parameters
        real(dp), intent(in)  :: a, b, c, alpha, beta, gamma
        real(dp), intent(out) :: va(3), vb(3), vc(3)
        real(dp) :: ar, br, gr, v3_z

        ar = alpha * pi / 180.0_dp
        br = beta  * pi / 180.0_dp
        gr = gamma * pi / 180.0_dp

        va = [a, 0.0_dp, 0.0_dp]
        vb = [b * cos(gr), b * sin(gr), 0.0_dp]
        vc(1) = c * cos(br)
        vc(2) = c * (cos(ar) - cos(br) * cos(gr)) / sin(gr)
        v3_z = c*c - vc(1)*vc(1) - vc(2)*vc(2)
        if (v3_z > 0.0_dp) then
            vc(3) = sqrt(v3_z)
        else
            vc(3) = 0.0_dp
        end if
    end subroutine lattice_vectors


    subroutine read_crystal_json_to_cif(json_path, cif, modified, iostat, iomsg)
        !! Read crystal_data.json back into cif_data_t.
        !! Sets modified=.true. if the JSON has "modified": true.
        character(len=*), intent(in) :: json_path
        type(cif_data_t), intent(out) :: cif
        logical, intent(out) :: modified
        integer, intent(out) :: iostat
        character(len=*), intent(out) :: iomsg

        integer :: unit, ios, atom_idx
        character(len=512) :: line
        logical :: in_atoms

        iostat = 0
        iomsg = ''
        modified = .false.
        cif%n_atoms = 0
        in_atoms = .false.

        open(newunit=unit, file=trim(json_path), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            iostat = IO_PARSE_ERROR
            iomsg = 'Cannot open JSON file: ' // trim(json_path)
            return
        end if

        ! ── First pass: count atoms ──
        do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (index(line, '"element"') > 0) then
                cif%n_atoms = cif%n_atoms + 1
            end if
        end do

        if (cif%n_atoms == 0) then
            iostat = IO_PARSE_ERROR
            iomsg = 'No atoms found in JSON file'
            close(unit)
            return
        end if

        allocate(cif%atoms(cif%n_atoms), stat=ios)
        if (ios /= 0) then
            iostat = IO_PARSE_ERROR
            iomsg = 'Memory allocation failed'
            close(unit)
            return
        end if

        ! ── Second pass: parse values ──
        rewind(unit)
        atom_idx = 0
        do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)

            ! Lattice parameters (guard against element names like Ca/Ba/Na)
            if (is_json_key(line, '"a"') .and. index(line, '"alpha"') == 0) then
                call extract_json_real(line, cif%a, ios)
            else if (is_json_key(line, '"b"') .and. index(line, '"beta"') == 0) then
                call extract_json_real(line, cif%b, ios)
            else if (is_json_key(line, '"c"')) then
                call extract_json_real(line, cif%c, ios)
            else if (index(line, '"alpha"') > 0) then
                call extract_json_real(line, cif%alpha, ios)
            else if (index(line, '"beta"') > 0) then
                call extract_json_real(line, cif%beta, ios)
            else if (index(line, '"gamma"') > 0) then
                call extract_json_real(line, cif%gamma, ios)

            ! positions_fractional flag
            else if (index(line, '"positions_fractional"') > 0) then
                if (index(line, 'true') > 0) cif%positions_fractional = .true.

            ! modified flag
            else if (index(line, '"modified"') > 0) then
                if (index(line, 'true') > 0) modified = .true.

            ! Atom data
            else if (index(line, '"element"') > 0) then
                atom_idx = atom_idx + 1
                if (atom_idx <= cif%n_atoms) then
                    call extract_json_string(line, cif%atoms(atom_idx)%element)
                end if
            else if (index(line, '"x"') > 0 .and. atom_idx > 0 .and. atom_idx <= cif%n_atoms) then
                call extract_json_real(line, cif%atoms(atom_idx)%x, ios)
            else if (index(line, '"y"') > 0 .and. atom_idx > 0 .and. atom_idx <= cif%n_atoms) then
                call extract_json_real(line, cif%atoms(atom_idx)%y, ios)
            else if (index(line, '"z"') > 0 .and. atom_idx > 0 .and. atom_idx <= cif%n_atoms) then
                call extract_json_real(line, cif%atoms(atom_idx)%z, ios)
            end if
        end do

        close(unit)
        cif%space_group = 'P1'
        cif%positions_fractional = .false.  ! viewer always outputs Cartesian
    end subroutine read_crystal_json_to_cif


    pure logical function is_json_key(line, key)
        !! True if key appears as a JSON key (not inside an element name like Ca/Ba)
        character(len=*), intent(in) :: line, key
        integer :: i
        is_json_key = .false.
        i = index(line, key)
        if (i == 0) return
        ! Character before key must not be a letter (a-z, A-Z)
        if (i > 1) then
            if ((line(i-1:i-1) >= 'A' .and. line(i-1:i-1) <= 'Z') .or. &
                (line(i-1:i-1) >= 'a' .and. line(i-1:i-1) <= 'z')) return
        end if
        is_json_key = .true.
    end function is_json_key

    subroutine extract_json_real(line, val, ios)
        !! Extract a real number from a JSON line like  "key": value,
        character(len=*), intent(in) :: line
        real(dp), intent(out) :: val
        integer, intent(out) :: ios
        integer :: colon_pos, end_pos
        character(len=256) :: tmp

        val = 0.0_dp
        ios = 0

        colon_pos = index(line, ':')
        if (colon_pos == 0) then
            ios = 1; return
        end if

        tmp = adjustl(line(colon_pos+1:))
        ! Remove trailing comma
        end_pos = len_trim(tmp)
        if (end_pos > 0) then
            if (tmp(end_pos:end_pos) == ',' .or. tmp(end_pos:end_pos) == '}') then
                tmp(end_pos:end_pos) = ' '
            end if
        end if

        read(tmp, *, iostat=ios) val
    end subroutine extract_json_real


    subroutine extract_json_string(line, str)
        !! Extract a quoted string from a JSON line like  "key": "value",
        character(len=*), intent(in) :: line
        character(len=*), intent(out) :: str
        integer :: colon_pos, quote_start, quote_end

        str = ''
        colon_pos = index(line, ':')
        if (colon_pos == 0) return
        ! Find opening quote after colon
        quote_start = index(line(colon_pos:), '"')
        if (quote_start == 0) return
        quote_start = colon_pos + quote_start - 1  ! convert to absolute position
        ! Find closing quote
        quote_end = index(line(quote_start+1:), '"')
        if (quote_end == 0) return
        str = line(quote_start+1:quote_start+quote_end-1)
    end subroutine extract_json_string

    ! ====================================================================
    ! Write crystal structure + phonon mode displacement data for viewer
    ! ====================================================================
    subroutine write_crystal_json_modes(modes_data, mode_index, json_path, scale_factor, iostat)
        use phonon_modes, only: phonon_modes_data_t
        type(phonon_modes_data_t), intent(in) :: modes_data
        integer, intent(in) :: mode_index
        character(len=*), intent(in) :: json_path
        real(dp), intent(in) :: scale_factor
        integer, intent(out) :: iostat

        real(dp) :: va(3), vb(3), vc(3), x, y, z, dx, dy, dz, contrib
        integer :: unit, i, m

        iostat = 0
        m = mode_index
        if (m < 1 .or. m > modes_data%n_branches) then
            iostat = 1; return
        end if

        open(newunit=unit, file=trim(json_path), status='replace', action='write', iostat=iostat)
        if (iostat /= 0) return

        ! Compute lattice vectors for Cartesian conversion
        call lattice_vectors(modes_data%cell_a, modes_data%cell_b, modes_data%cell_c, &
            modes_data%cell_alpha, modes_data%cell_beta, modes_data%cell_gamma, va, vb, vc)

        write(unit, '(a)') '{'
        write(unit, '(a)') '  "lattice": {'
        write(unit, '(a,f10.4,a)') '    "a": ', modes_data%cell_a, ','
        write(unit, '(a,f10.4,a)') '    "b": ', modes_data%cell_b, ','
        write(unit, '(a,f10.4,a)') '    "c": ', modes_data%cell_c, ','
        write(unit, '(a,f10.4,a)') '    "alpha": ', modes_data%cell_alpha, ','
        write(unit, '(a,f10.4,a)') '    "beta": ', modes_data%cell_beta, ','
        write(unit, '(a,f10.4)') '    "gamma": ', modes_data%cell_gamma
        write(unit, '(a)') '  },'
        write(unit, '(a)') '  "positions_fractional": false,'
        write(unit, '(a)') '  "atoms": ['

        ! Write atoms in Cartesian coordinates
        do i = 1, modes_data%n_ions
            x = modes_data%ion_positions_frac(1,i) * va(1) &
              + modes_data%ion_positions_frac(2,i) * vb(1) &
              + modes_data%ion_positions_frac(3,i) * vc(1)
            y = modes_data%ion_positions_frac(1,i) * va(2) &
              + modes_data%ion_positions_frac(2,i) * vb(2) &
              + modes_data%ion_positions_frac(3,i) * vc(2)
            z = modes_data%ion_positions_frac(1,i) * va(3) &
              + modes_data%ion_positions_frac(2,i) * vb(3) &
              + modes_data%ion_positions_frac(3,i) * vc(3)
            write(unit, '(a,a,a,a,f10.6,a,f10.6,a,f10.6,a)', advance='no') &
                '    { "element": "', trim(modes_data%ion_species(i)), '", ', &
                '"x": ', x, ', "y": ', y, ', "z": ', z, ' }'
            if (i < modes_data%n_ions) then
                write(unit, '(a)') ','
            else
                write(unit, '(a)') ''
            end if
        end do

        write(unit, '(a)') '  ],'
        write(unit, '(a)') '  "phonon_modes": {'
        write(unit, '(a,i0,a)') '    "mode_index": ', mode_index, ','
        write(unit, '(a,f12.4,a)') '    "frequency": ', modes_data%modes(m)%frequency, ','
        write(unit, '(a,f12.4,a)') '    "ir_intensity": ', modes_data%modes(m)%ir_intensity, ','
        write(unit, '(a,f12.4,a)') '    "raman_activity": ', modes_data%modes(m)%raman_activity, ','
        write(unit, '(a,f12.6,a)') '    "mode_charge_norm": ', modes_data%modes(m)%mode_charge_norm, ','
        write(unit, '(a)') '    "atom_displacements": ['

        do i = 1, modes_data%n_ions
            dx = modes_data%modes(m)%displacements(i,1) * scale_factor
            dy = modes_data%modes(m)%displacements(i,2) * scale_factor
            dz = modes_data%modes(m)%displacements(i,3) * scale_factor
            contrib = modes_data%modes(m)%atom_contributions(i)
            write(unit, '(a,f10.6,a,f10.6,a,f10.6,a,f10.6,a)', advance='no') &
                '      { "dx": ', dx, ', "dy": ', dy, ', "dz": ', dz, &
                ', "contribution": ', contrib, ' }'
            if (i < modes_data%n_ions) then
                write(unit, '(a)') ','
            else
                write(unit, '(a)') ''
            end if
        end do

        write(unit, '(a)') '    ]'
        write(unit, '(a)') '  }'
        write(unit, '(a)') '}'
        close(unit)
    end subroutine write_crystal_json_modes

end module crystal_json
