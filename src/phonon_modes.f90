! phonon_modes.f90 — Phonon eigenvector parsing, Born charge parsing, mode decomposition
module phonon_modes
    use castep_config, only: dp, MAX_LINE_LEN, IO_FILE_NOT_FOUND, IO_USER_QUIT
    use iso_fortran_env, only: iostat_end
    implicit none
    private

    public :: born_charge_t, phonon_mode_t, phonon_modes_data_t
    public :: parse_phonon_eigenvectors, parse_castep_born_charges
    public :: compute_mode_decomposition, free_phonon_modes_data

    ! ── Error codes ──
    integer, parameter, public :: IO_EIGENVECTORS_NOT_FOUND = 120
    integer, parameter, public :: IO_BORN_MISMATCH = 121
    integer, parameter, public :: IO_BORN_NOT_FOUND = 122

    ! ── Data types ──

    type :: born_charge_t
        real(dp) :: tensor(3,3) = 0.0_dp
    end type born_charge_t

    type :: phonon_mode_t
        integer  :: mode_index = 0
        real(dp) :: frequency = 0.0_dp       ! cm⁻¹
        real(dp) :: ir_intensity = 0.0_dp    ! (D/A)²/amu
        real(dp) :: raman_activity = 0.0_dp  ! A⁴/amu
        real(dp) :: mode_charge(3) = 0.0_dp  ! p_m vector
        real(dp) :: mode_charge_norm = 0.0_dp
        real(dp), allocatable :: eigenvectors(:,:)    ! (n_ions, 6)
        real(dp), allocatable :: displacements(:,:)  ! (n_ions, 3) Cartesian
        real(dp), allocatable :: atom_contributions(:) ! (n_ions) 0..1
    end type phonon_mode_t

    type :: phonon_modes_data_t
        integer :: n_ions = 0
        integer :: n_branches = 0
        integer :: n_qpoints = 0
        ! Lattice (from .phonon header)
        real(dp) :: lattice_vectors(3,3) = 0.0_dp
        real(dp) :: cell_a = 0.0_dp, cell_b = 0.0_dp, cell_c = 0.0_dp
        real(dp) :: cell_alpha = 0.0_dp, cell_beta = 0.0_dp, cell_gamma = 0.0_dp
        ! Structure
        character(len=6), allocatable :: ion_species(:)      ! (n_ions)
        real(dp), allocatable :: ion_positions_frac(:,:)      ! (3, n_ions)
        real(dp), allocatable :: ion_masses(:)                ! (n_ions) AMU
        ! Born charges
        type(born_charge_t), allocatable :: born_charges(:)   ! (n_ions)
        logical :: has_born_charges = .false.
        ! Modes
        type(phonon_mode_t), allocatable :: modes(:)          ! (n_branches)
    end type phonon_modes_data_t

contains

    ! ====================================================================
    ! Parse .phonon file: header (structure) + frequencies + eigenvectors
    ! ====================================================================
    subroutine parse_phonon_eigenvectors(filename, data, iostat_out, iomsg)
        character(len=*), intent(in) :: filename
        type(phonon_modes_data_t), intent(inout) :: data
        integer, intent(out) :: iostat_out
        character(len=*), intent(out) :: iomsg

        integer :: unit, ios, i, j, mode_idx, ion_idx, nb, ni
        real(dp) :: frac(3), amass, freq, ir, raman
        real(dp) :: eig(6)
        character(len=MAX_LINE_LEN) :: line, species
        character(len=6) :: tmp_species

        iostat_out = 0; iomsg = ''

        open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios, iomsg=iomsg)
        if (ios /= 0) then
            iostat_out = IO_FILE_NOT_FOUND; return
        end if

        ! ── Parse header ──
        do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, 'Number of ions') > 0) then
                read(line, *, iostat=ios) tmp_species, tmp_species, tmp_species, data%n_ions
            else if (index(line, 'Number of branches') > 0) then
                read(line, *, iostat=ios) tmp_species, tmp_species, tmp_species, data%n_branches
            else if (index(line, 'Number of wavevectors') > 0) then
                read(line, *, iostat=ios) tmp_species, tmp_species, tmp_species, data%n_qpoints
            else if (index(line, 'Unit cell vectors (A)') > 0) then
                do i = 1, 3
                    read(unit, *, iostat=ios) (data%lattice_vectors(i,j), j=1,3)
                end do
            else if (index(line, 'END header') > 0) then
                exit
            end if
        end do

        ! Compute cell params from lattice vectors
        call compute_cell_params(data)

        ni = data%n_ions
        nb = data%n_branches
        if (ni <= 0 .or. nb <= 0) then
            iostat_out = 1; iomsg = 'Invalid ion/branch count in .phonon header'; close(unit); return
        end if

        ! Allocate
        allocate(data%ion_species(ni), data%ion_positions_frac(3,ni), data%ion_masses(ni))
        allocate(data%modes(nb))
        do i = 1, nb
            allocate(data%modes(i)%eigenvectors(ni, 6))
            allocate(data%modes(i)%displacements(ni, 3))
            allocate(data%modes(i)%atom_contributions(ni))
            data%modes(i)%mode_index = i
            data%modes(i)%eigenvectors = 0.0_dp
            data%modes(i)%displacements = 0.0_dp
            data%modes(i)%atom_contributions = 0.0_dp
        end do

        ! ── Parse fractional coordinates + species + masses from header ──
        ! Re-open to re-read the fractional coords block
        rewind(unit)
        do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, 'Fractional Co-ordinates') > 0) exit
        end do

        if (ios == 0) then
            do i = 1, ni
                read(unit, *, iostat=ios) ion_idx, frac(1), frac(2), frac(3), species, amass
                if (ios /= 0) then
                    ! Try reading without mass (some formats lack it)
                    ! The current format has mass, let's just handle errors
                    iostat_out = 2; iomsg = 'Error reading fractional coords in .phonon'; close(unit); return
                end if
                data%ion_species(i) = trim(adjustl(species))
                data%ion_positions_frac(1,i) = frac(1)
                data%ion_positions_frac(2,i) = frac(2)
                data%ion_positions_frac(3,i) = frac(3)
                data%ion_masses(i) = amass
            end do
        end if

        ! ── Skip to frequency block ──
        rewind(unit)
        do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (line(1:5) == '     q-pt=' .or. index(line, 'q-pt=') > 0) exit
        end do

        if (ios /= 0) then
            iostat_out = 3; iomsg = 'No frequency data found in .phonon'; close(unit); return
        end if

        ! ── Read frequencies, IR, Raman for each branch ──
        ! Format: mode_idx  freq  ir_intensity  [raman_activity]
        ! Raman column is optional (only present if calculate_raman=true in CASTEP)
        do i = 1, nb
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) then
                iostat_out = 4; iomsg = 'Error reading frequency table'; close(unit); return
            end if
            ! Try 4 values first (with Raman), fall back to 3 (no Raman)
            read(line, *, iostat=ios) mode_idx, freq, ir, raman
            if (ios /= 0) then
                read(line, *, iostat=ios) mode_idx, freq, ir
                raman = 0.0_dp
            end if
            if (ios /= 0) then
                iostat_out = 4; iomsg = 'Error reading frequency table'; close(unit); return
            end if
            data%modes(i)%frequency = freq
            data%modes(i)%ir_intensity = ir
            data%modes(i)%raman_activity = raman
        end do

        ! ── Skip to "Phonon Eigenvectors" block ──
        do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, 'Phonon Eigenvectors') > 0) exit
        end do

        if (ios /= 0) then
            iostat_out = IO_EIGENVECTORS_NOT_FOUND
            iomsg = 'No eigenvector data found in .phonon — ensure PHONON_WRITE_EIGENVECTORS=TRUE'
            close(unit); return
        end if

        ! Skip the header line: "Mode Ion                X ..."
        read(unit, '(a)', iostat=ios) line

        ! ── Read eigenvectors (only first q-point = Gamma) ──
        do i = 1, nb
            do j = 1, ni
                read(unit, *, iostat=ios) mode_idx, ion_idx, &
                    eig(1), eig(2), eig(3), eig(4), eig(5), eig(6)
                if (ios /= 0) then
                    iostat_out = 5; iomsg = 'Error reading eigenvectors'; close(unit); return
                end if
                data%modes(i)%eigenvectors(ion_idx, 1) = eig(1)  ! X_real
                data%modes(i)%eigenvectors(ion_idx, 2) = eig(2)  ! X_imag
                data%modes(i)%eigenvectors(ion_idx, 3) = eig(3)  ! Y_real
                data%modes(i)%eigenvectors(ion_idx, 4) = eig(4)  ! Y_imag
                data%modes(i)%eigenvectors(ion_idx, 5) = eig(5)  ! Z_real
                data%modes(i)%eigenvectors(ion_idx, 6) = eig(6)  ! Z_imag
            end do
        end do

        close(unit)

        ! Pre-compute Cartesian displacements from eigenvectors (Gamma: real parts only)
        do i = 1, nb
            do j = 1, ni
                if (data%ion_masses(j) > 0.0_dp) then
                    data%modes(i)%displacements(j, 1) = &
                        data%modes(i)%eigenvectors(j, 1) / sqrt(data%ion_masses(j))
                    data%modes(i)%displacements(j, 2) = &
                        data%modes(i)%eigenvectors(j, 3) / sqrt(data%ion_masses(j))
                    data%modes(i)%displacements(j, 3) = &
                        data%modes(i)%eigenvectors(j, 5) / sqrt(data%ion_masses(j))
                end if
            end do
        end do

        iostat_out = 0
    end subroutine parse_phonon_eigenvectors

    ! ====================================================================
    ! Parse Born effective charges from .castep file
    ! ====================================================================
    subroutine parse_castep_born_charges(filename, data, iostat_out, iomsg)
        character(len=*), intent(in) :: filename
        type(phonon_modes_data_t), intent(inout) :: data
        integer, intent(out) :: iostat_out
        character(len=*), intent(out) :: iomsg

        integer :: unit, ios, i, n_born
        character(len=MAX_LINE_LEN) :: line
        character(len=6) :: element
        integer :: ion_idx

        iostat_out = 0; iomsg = ''

        open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios, iomsg=iomsg)
        if (ios /= 0) then
            iostat_out = IO_FILE_NOT_FOUND; return
        end if

        ! ── Find "Born Effective Charges" section header ──
        do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) then
                iostat_out = IO_BORN_NOT_FOUND
                iomsg = 'Born Effective Charges section not found in .castep'
                close(unit); return
            end if
            if (index(line, 'Born Effective Charges') > 0) exit
        end do

        ! Skip the dashed line
        read(unit, '(a)', iostat=ios) line

        ! Allocate Born charges
        if (.not. allocated(data%born_charges)) then
            allocate(data%born_charges(data%n_ions))
        end if
        do i = 1, data%n_ions
            data%born_charges(i)%tensor = 0.0_dp
        end do

        ! ── Read Born charges ──
        n_born = 0
        do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit

            ! Check for end of section (line of = signs)
            if (len_trim(line) > 0) then
                if (line(1:3) == '===') exit
            end if

            ! Parse: element, ion_idx, tensor row 1 (xx, xy, xz)
            read(line, *, iostat=ios) element, ion_idx, &
                data%born_charges(ion_idx)%tensor(1,1), &
                data%born_charges(ion_idx)%tensor(1,2), &
                data%born_charges(ion_idx)%tensor(1,3)
            if (ios /= 0) cycle  ! skip non-data lines

            ! Read row 2
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            read(line, *, iostat=ios) &
                data%born_charges(ion_idx)%tensor(2,1), &
                data%born_charges(ion_idx)%tensor(2,2), &
                data%born_charges(ion_idx)%tensor(2,3)

            ! Read row 3
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            read(line, *, iostat=ios) &
                data%born_charges(ion_idx)%tensor(3,1), &
                data%born_charges(ion_idx)%tensor(3,2), &
                data%born_charges(ion_idx)%tensor(3,3)

            n_born = n_born + 1
        end do

        close(unit)

        if (n_born /= data%n_ions) then
            iostat_out = IO_BORN_MISMATCH
            write(iomsg, '(a,i0,a,i0)') 'Born charge count (', n_born, &
                ') does not match ion count (', data%n_ions, ')'
            data%has_born_charges = .false.
            return
        end if

        data%has_born_charges = .true.
        iostat_out = 0
    end subroutine parse_castep_born_charges

    ! ====================================================================
    ! Compute mode decomposition using Born effective charges
    ! ====================================================================
    subroutine compute_mode_decomposition(data, iostat_out, iomsg)
        type(phonon_modes_data_t), intent(inout) :: data
        integer, intent(out) :: iostat_out
        character(len=*), intent(out) :: iomsg

        integer :: m, k, i, j
        real(dp) :: u(3), zu(3), p_m(3), pm_norm, ck, max_disp
        real(dp) :: eps

        iostat_out = 0; iomsg = ''
        eps = 1.0e-12_dp

        do m = 1, data%n_branches
            p_m = 0.0_dp
            zu = 0.0_dp

            if (data%has_born_charges) then
                do k = 1, data%n_ions
                    ! Cartesian displacement: already pre-computed in displacements(:,:)
                    u(1) = data%modes(m)%displacements(k, 1)
                    u(2) = data%modes(m)%displacements(k, 2)
                    u(3) = data%modes(m)%displacements(k, 3)

                    ! Z*(kappa) · u(kappa,m)
                    do i = 1, 3
                        zu(i) = 0.0_dp
                        do j = 1, 3
                            zu(i) = zu(i) + data%born_charges(k)%tensor(i,j) * u(j)
                        end do
                    end do

                    p_m = p_m + zu
                end do
            end if

            pm_norm = sqrt(dot_product(p_m, p_m))
            data%modes(m)%mode_charge = p_m
            data%modes(m)%mode_charge_norm = pm_norm

            ! Per-atom contributions
            if (data%has_born_charges .and. pm_norm > eps) then
                do k = 1, data%n_ions
                    u(1) = data%modes(m)%displacements(k, 1)
                    u(2) = data%modes(m)%displacements(k, 2)
                    u(3) = data%modes(m)%displacements(k, 3)
                    do i = 1, 3
                        zu(i) = 0.0_dp
                        do j = 1, 3
                            zu(i) = zu(i) + data%born_charges(k)%tensor(i,j) * u(j)
                        end do
                    end do
                    ck = sqrt(dot_product(zu, zu)) / pm_norm
                    data%modes(m)%atom_contributions(k) = min(ck, 1.0_dp)
                end do
            else
                data%modes(m)%atom_contributions = 0.0_dp
            end if
        end do

        ! Scale displacements to visual magnitude
        ! Normalize: max displacement across all modes = 1.0 Å
        max_disp = 0.0_dp
        do m = 1, data%n_branches
            do k = 1, data%n_ions
                do i = 1, 3
                    max_disp = max(max_disp, abs(data%modes(m)%displacements(k,i)))
                end do
            end do
        end do

        if (max_disp > eps) then
            do m = 1, data%n_branches
                data%modes(m)%displacements = data%modes(m)%displacements / max_disp
            end do
        end if

        iostat_out = 0
    end subroutine compute_mode_decomposition

    ! ====================================================================
    ! Compute cell parameters (a,b,c,alpha,beta,gamma) from lattice vectors
    ! ====================================================================
    subroutine compute_cell_params(data)
        type(phonon_modes_data_t), intent(inout) :: data
        real(dp) :: norm1, norm2, norm3, dot_prod, pi
        pi = acos(-1.0_dp)

        norm1 = sqrt(sum(data%lattice_vectors(1,:)**2))
        norm2 = sqrt(sum(data%lattice_vectors(2,:)**2))
        norm3 = sqrt(sum(data%lattice_vectors(3,:)**2))

        data%cell_a = norm1
        data%cell_b = norm2
        data%cell_c = norm3

        if (norm2 > 0.0_dp .and. norm3 > 0.0_dp) then
            dot_prod = dot_product(data%lattice_vectors(2,:), data%lattice_vectors(3,:))
            data%cell_alpha = acos(dot_prod / (norm2 * norm3)) * 180.0_dp / pi
        else
            data%cell_alpha = 90.0_dp
        end if

        if (norm1 > 0.0_dp .and. norm3 > 0.0_dp) then
            dot_prod = dot_product(data%lattice_vectors(1,:), data%lattice_vectors(3,:))
            data%cell_beta = acos(dot_prod / (norm1 * norm3)) * 180.0_dp / pi
        else
            data%cell_beta = 90.0_dp
        end if

        if (norm1 > 0.0_dp .and. norm2 > 0.0_dp) then
            dot_prod = dot_product(data%lattice_vectors(1,:), data%lattice_vectors(2,:))
            data%cell_gamma = acos(dot_prod / (norm1 * norm2)) * 180.0_dp / pi
        else
            data%cell_gamma = 90.0_dp
        end if
    end subroutine compute_cell_params

    ! ====================================================================
    ! Free allocated memory
    ! ====================================================================
    subroutine free_phonon_modes_data(data)
        type(phonon_modes_data_t), intent(inout) :: data
        integer :: i

        if (allocated(data%modes)) then
            do i = 1, size(data%modes)
                if (allocated(data%modes(i)%eigenvectors)) deallocate(data%modes(i)%eigenvectors)
                if (allocated(data%modes(i)%displacements)) deallocate(data%modes(i)%displacements)
                if (allocated(data%modes(i)%atom_contributions)) deallocate(data%modes(i)%atom_contributions)
            end do
            deallocate(data%modes)
        end if

        if (allocated(data%ion_species)) deallocate(data%ion_species)
        if (allocated(data%ion_positions_frac)) deallocate(data%ion_positions_frac)
        if (allocated(data%ion_masses)) deallocate(data%ion_masses)
        if (allocated(data%born_charges)) deallocate(data%born_charges)

        data%n_ions = 0
        data%n_branches = 0
        data%n_qpoints = 0
        data%has_born_charges = .false.
    end subroutine free_phonon_modes_data

end module phonon_modes
