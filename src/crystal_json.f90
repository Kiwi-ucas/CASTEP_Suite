module crystal_json
    !! Generate crystal_data.json for the Rust Crystal Viewer
    use castep_config, only: dp, pi, castep_config_t, cif_data_t
    use parser, only: clean_element_symbol
    implicit none
    private

    public :: write_crystal_json, write_crystal_json_cif

contains

    subroutine write_crystal_json(cfg, json_path, iostat)
        !! Write crystal structure data from castep_config_t
        type(castep_config_t), intent(in) :: cfg
        character(len=*), intent(in) :: json_path
        integer, intent(out) :: iostat
        integer :: unit, i

        iostat = 0
        open(newunit=unit, file=trim(json_path), status='replace', action='write', iostat=iostat)
        if (iostat /= 0) return

        write(unit, '(a)') '{'
        write(unit, '(a)') '  "lattice": {'
        write(unit, '(a,f10.4,a)') '    "a": ', cfg%cell_length(1), ','
        write(unit, '(a,f10.4,a)') '    "b": ', cfg%cell_length(2), ','
        write(unit, '(a,f10.4,a)') '    "c": ', cfg%cell_length(3), ','
        write(unit, '(a,f10.4,a)') '    "alpha": ', cfg%cell_angle(1), ','
        write(unit, '(a,f10.4,a)') '    "beta": ', cfg%cell_angle(2), ','
        write(unit, '(a,f10.4)') '    "gamma": ', cfg%cell_angle(3)
        write(unit, '(a)') '  },'
        write(unit, '(a)') '  "atoms": ['

        do i = 1, cfg%num_atoms
            write(unit, '(a,a,a,a,f10.6,a,a,f10.6,a,a,f10.6,a)', advance='no') &
                '    { "element": "', trim(clean_element_symbol(cfg%atom_type(i))), '", ', &
                '"x": ', cfg%atom_x(i), ', ', &
                '"y": ', cfg%atom_y(i), ', ', &
                '"z": ', cfg%atom_z(i), ' }'
            if (i < cfg%num_atoms) write(unit, '(a)') ','
            if (i == cfg%num_atoms) write(unit, '(a)') ''
        end do

        write(unit, '(a)') '  ]'
        write(unit, '(a)') '}'
        close(unit)
    end subroutine write_crystal_json


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

end module crystal_json
