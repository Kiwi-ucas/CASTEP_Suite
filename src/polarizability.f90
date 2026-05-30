module polarizability
    !! Static polarizability via AIMD polarization fluctuation method
    !! Combines CASTEP DFPT optical dielectric tensor with CP2K dipole trajectory
    use castep_config, only: dp, pi, EPSILON_0, KBOLTZMANN, DEBYE_TO_CM, &
        ANG3_TO_M3, DEBYE_PER_ANG, MAX_LINE_LEN, &
        IO_EPS_NOT_FOUND, IO_EPS_PARSE_ERROR, IO_DIPOLE_ERROR
    implicit none
    private

    public :: pol_data_t, parse_castep_epsilon, parse_cp2k_dipoles, &
              unwrap_dipoles, detrend_dipoles, &
              compute_static_dielectric, compute_polarizability, &
              free_pol_data

    type :: pol_data_t
        ! CASTEP optical dielectric tensor ε_∞ (3×3)
        real(dp) :: eps_inf(3,3) = 0.0_dp
        ! CP2K metadata
        integer  :: n_frames = 0
        real(dp) :: temperature = 0.0_dp
        real(dp) :: volume_ang3 = 0.0_dp
        real(dp) :: cell_abc(3) = 0.0_dp  ! a, b, c in Angstrom
        ! Raw dipole trajectory (n_frames, 3) in Debye
        real(dp), allocatable :: dipoles(:,:)
        ! Unwrapped + detrended dipoles in Debye
        real(dp), allocatable :: dipoles_clean(:,:)
        ! Results (3×3 tensors)
        real(dp) :: eps_static(3,3) = 0.0_dp
        real(dp) :: eps_ion(3,3) = 0.0_dp
        real(dp) :: alpha_static(3,3) = 0.0_dp
        real(dp) :: alpha_ion(3,3) = 0.0_dp
        real(dp) :: alpha_inf(3,3) = 0.0_dp
        ! Unwrap statistics
        integer :: n_unwraps = 0
        ! Detrend statistics
        real(dp) :: drift_rate(3) = 0.0_dp  ! Debye/ps
    end type

    integer, parameter :: MAX_FILES = 200000
    character(len=*), parameter :: EOM = 'No optical dielectric tensor found'

contains

    subroutine parse_castep_epsilon(filename, eps_inf, iostat, iomsg)
        !! Extract optical dielectric tensor ε_∞ from CASTEP .castep file
        character(len=*), intent(in)  :: filename
        real(dp), intent(out)         :: eps_inf(3,3)
        integer, intent(out)          :: iostat
        character(len=*), intent(out), optional :: iomsg
        character(len=MAX_LINE_LEN) :: line
        integer :: unit, ios, idx, i
        logical :: found

        iostat = 0
        eps_inf = 0.0_dp

        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            iostat = IO_EPS_NOT_FOUND
            if (present(iomsg)) iomsg = 'Cannot open file: ' // trim(filename)
            return
        end if

        found = .false.
        do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            idx = index(line, 'Optical Permittivity')
            if (idx > 0) then
                ! Skip separator line
                read(unit, '(a)', iostat=ios) line
                if (ios /= 0) exit
                ! Read 3×3 matrix
                do i = 1, 3
                    read(unit, '(a)', iostat=ios) line
                    if (ios /= 0) exit
                    read(line, *, iostat=ios) eps_inf(i, 1:3)
                    if (ios /= 0) exit
                end do
                if (ios == 0) found = .true.
                exit
            end if
        end do
        close(unit)

        if (.not. found) then
            iostat = IO_EPS_PARSE_ERROR
            if (present(iomsg)) iomsg = trim(EOM)
        end if
    end subroutine parse_castep_epsilon


    subroutine parse_cp2k_dipoles(dir_path, data, n_frames_exp, iostat, iomsg)
        !! Read CP2K dipole files from directory
        !! Each file has format:
        !!   Dipole moment [Debye]
        !!     X= xxxxx Y= yyyyy Z= zzzzz  Total= ttttt
        character(len=*), intent(in)  :: dir_path
        type(pol_data_t), intent(inout) :: data
        integer, intent(in)           :: n_frames_exp  ! expected number of frames
        integer, intent(out)          :: iostat
        character(len=*), intent(out), optional :: iomsg
        character(len=MAX_LINE_LEN) :: line, fname, cmd
        character(len=32) :: file_list
        integer :: unit, ios, iframe, n_files
        real(dp) :: dx, dy, dz

        iostat = 0

        ! Generate file list
        file_list = '.dipole_list_tmp'
        cmd = 'ls ' // trim(dir_path) // '/dipole* > ' // trim(file_list) // ' 2>/dev/null'
        call execute_command_line(trim(cmd), wait=.true.)

        open(newunit=unit, file=trim(file_list), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            iostat = IO_DIPOLE_ERROR
            if (present(iomsg)) iomsg = 'Cannot list dipole files in: ' // trim(dir_path)
            return
        end if

        ! Count files
        n_files = 0
        do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) > 0) n_files = n_files + 1
        end do

        if (n_files == 0) then
            close(unit, status='delete')
            iostat = IO_DIPOLE_ERROR
            if (present(iomsg)) iomsg = 'No dipole files found in: ' // trim(dir_path)
            return
        end if

        ! Allocate arrays
        data%n_frames = n_files
        allocate(data%dipoles(n_files, 3))
        allocate(data%dipoles_clean(n_files, 3))
        data%dipoles = 0.0_dp
        data%dipoles_clean = 0.0_dp

        ! Read all files
        rewind(unit)
        do iframe = 1, n_files
            read(unit, '(a)', iostat=ios) fname
            if (ios /= 0) then
                iostat = IO_DIPOLE_ERROR
                if (present(iomsg)) iomsg = 'Error reading file list at frame'
                close(unit, status='delete')
                return
            end if

            open(newunit=unit, file=trim(fname), status='old', action='read', iostat=ios)
            if (ios /= 0) then
                iostat = IO_DIPOLE_ERROR
                if (present(iomsg)) iomsg = 'Cannot open: ' // trim(fname)
                close(unit, status='delete')
                return
            end if

            ! Read header line
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) then
                close(unit)
                iostat = IO_DIPOLE_ERROR
                if (present(iomsg)) iomsg = 'Error reading header in: ' // trim(fname)
                return
            end if

            ! Read dipole values: "  X= xxxxx Y= yyyyy Z= zzzzz  Total= ttttt"
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) then
                close(unit)
                iostat = IO_DIPOLE_ERROR
                if (present(iomsg)) iomsg = 'Error reading dipole in: ' // trim(fname)
                return
            end if
            close(unit)

            ! Parse X= Y= Z= values
            call parse_dipole_line(line, dx, dy, dz, ios)
            if (ios /= 0) then
                iostat = IO_DIPOLE_ERROR
                if (present(iomsg)) iomsg = 'Cannot parse dipole values in: ' // trim(fname)
                return
            end if
            data%dipoles(iframe, 1) = dx
            data%dipoles(iframe, 2) = dy
            data%dipoles(iframe, 3) = dz
        end do
        close(unit, status='delete')

        ! Warn if file count differs from expected
        if (n_frames_exp > 0 .and. n_files /= n_frames_exp) then
            ! Not an error, just note it
        end if
    end subroutine parse_cp2k_dipoles


    subroutine parse_dipole_line(line, dx, dy, dz, iostat)
        !! Parse "  X= xxxxx Y= yyyyy Z= zzzzz  Total= ttttt" line
        character(len=*), intent(in) :: line
        real(dp), intent(out) :: dx, dy, dz
        integer, intent(out)  :: iostat
        character(len=MAX_LINE_LEN) :: tmp
        integer :: ix, iy, iz

        iostat = 0
        dx = 0.0_dp; dy = 0.0_dp; dz = 0.0_dp

        ! Find positions of X= Y= Z=
        ix = index(line, 'X=')
        iy = index(line, 'Y=')
        iz = index(line, 'Z=')

        if (ix == 0 .or. iy == 0 .or. iz == 0) then
            iostat = -1
            return
        end if

        ! Extract values
        read(line(ix+2:iy-1), *, iostat=iostat) dx
        if (iostat /= 0) return
        read(line(iy+2:iz-1), *, iostat=iostat) dy
        if (iostat /= 0) return
        ! Z extends to end of line or "Total="
        tmp = line(iz+2:)
        iz = index(tmp, 'Total')
        if (iz > 0) then
            read(tmp(1:iz-1), *, iostat=iostat) dz
        else
            read(tmp, *, iostat=iostat) dz
        end if
    end subroutine parse_dipole_line


    subroutine unwrap_dipoles(data)
        !! Unwrap Berry phase polarization quantum jumps
        !! Pq,α = DEBYE_PER_ANG × a_α [Debye]
        !! If |ΔD| > Pq/2, apply correction ±Pq
        type(pol_data_t), intent(inout) :: data
        real(dp) :: pq(3), delta
        integer :: i, j

        ! Compute polarization quantum for each direction
        do j = 1, 3
            pq(j) = DEBYE_PER_ANG * data%cell_abc(j)
        end do

        data%n_unwraps = 0

        do i = 2, data%n_frames
            do j = 1, 3
                delta = data%dipoles(i, j) - data%dipoles(i-1, j)
                if (delta > pq(j) * 0.5_dp) then
                    data%dipoles(i, j) = data%dipoles(i, j) - pq(j)
                    data%n_unwraps = data%n_unwraps + 1
                else if (delta < -pq(j) * 0.5_dp) then
                    data%dipoles(i, j) = data%dipoles(i, j) + pq(j)
                    data%n_unwraps = data%n_unwraps + 1
                end if
            end do
        end do

        ! Copy to clean array
        data%dipoles_clean = data%dipoles
    end subroutine unwrap_dipoles


    subroutine detrend_dipoles(data, time_step_ps)
        !! Remove linear drift from Li+ diffusion
        !! D_detrend(t) = D(t) - (s*t + b) via least squares
        type(pol_data_t), intent(inout) :: data
        real(dp), intent(in) :: time_step_ps  ! MD time step in ps
        real(dp) :: sum_t, sum_t2, sum_d(3), sum_td(3)
        real(dp) :: denom, s, b, n
        integer :: i, j

        n = real(data%n_frames, dp)
        sum_t = 0.0_dp; sum_t2 = 0.0_dp
        sum_d = 0.0_dp; sum_td = 0.0_dp

        ! Accumulate sums
        do i = 1, data%n_frames
            sum_t = sum_t + real(i, dp)
            sum_t2 = sum_t2 + real(i, dp)**2
            do j = 1, 3
                sum_d(j) = sum_d(j) + data%dipoles_clean(i, j)
                sum_td(j) = sum_td(j) + real(i, dp) * data%dipoles_clean(i, j)
            end do
        end do

        denom = n * sum_t2 - sum_t * sum_t
        if (abs(denom) < 1.0e-30_dp) return

        ! Linear fit and detrend for each direction
        do j = 1, 3
            s = (n * sum_td(j) - sum_t * sum_d(j)) / denom
            b = (sum_d(j) - s * sum_t) / n
            ! Store drift rate (Debye/ps): s is per frame, convert with time_step_ps
            data%drift_rate(j) = s / time_step_ps
            do i = 1, data%n_frames
                data%dipoles_clean(i, j) = data%dipoles_clean(i, j) - (s * real(i, dp) + b)
            end do
        end do
    end subroutine detrend_dipoles


    subroutine compute_static_dielectric(data, iostat, iomsg)
        !! Compute static dielectric tensor
        !! ε_αβ = ε_∞,αβ + Cov(M_α, M_β) / (ε₀ k_B T Ω)
        type(pol_data_t), intent(inout) :: data
        integer, intent(out)  :: iostat
        character(len=*), intent(out), optional :: iomsg
        real(dp) :: m_sum(3), m2_sum(3,3), m_mean(3), cov(3,3)
        real(dp) :: vol_m3, kt, denom
        integer :: i, j, k, n

        iostat = 0
        n = data%n_frames

        if (n < 2) then
            iostat = IO_DIPOLE_ERROR
            if (present(iomsg)) iomsg = 'Need at least 2 frames for statistics'
            return
        end if

        if (data%temperature <= 0.0_dp) then
            iostat = IO_DIPOLE_ERROR
            if (present(iomsg)) iomsg = 'Temperature must be positive'
            return
        end if

        if (data%volume_ang3 <= 0.0_dp) then
            iostat = IO_DIPOLE_ERROR
            if (present(iomsg)) iomsg = 'Cell volume must be positive'
            return
        end if

        ! Convert volume to m³
        vol_m3 = data%volume_ang3 * ANG3_TO_M3

        ! k_B * T in Joules
        kt = KBOLTZMANN * data%temperature

        ! Precompute denominator: ε₀ * k_B * T * Ω
        denom = EPSILON_0 * kt * vol_m3

        ! Compute mean of each component
        m_sum = 0.0_dp
        do i = 1, n
            do j = 1, 3
                m_sum(j) = m_sum(j) + data%dipoles_clean(i, j)
            end do
        end do
        m_mean = m_sum / real(n, dp)

        ! Compute covariance matrix
        m2_sum = 0.0_dp
        do i = 1, n
            do j = 1, 3
                do k = 1, 3
                    m2_sum(j, k) = m2_sum(j, k) + &
                        (data%dipoles_clean(i, j) - m_mean(j)) * &
                        (data%dipoles_clean(i, k) - m_mean(k))
                end do
            end do
        end do
        cov = m2_sum / real(n, dp)

        ! Convert covariance from Debye² to (C·m)²
        ! M_α = D_α × DEBYE_TO_CM
        cov = cov * DEBYE_TO_CM**2

        ! ε_ion = Cov / (ε₀ k_B T Ω)
        data%eps_ion = cov / denom

        ! ε_static = ε_∞ + ε_ion
        data%eps_static = data%eps_inf + data%eps_ion
    end subroutine compute_static_dielectric


    subroutine compute_polarizability(data)
        !! Convert dielectric tensor to polarizability (Å³)
        !! α_αβ = (ε_αβ − δ_αβ) × Ω / (4π)
        type(pol_data_t), intent(inout) :: data
        real(dp) :: factor, kronecker
        integer :: i, j

        factor = data%volume_ang3 / (4.0_dp * pi)

        do i = 1, 3
            do j = 1, 3
                kronecker = 0.0_dp
                if (i == j) kronecker = 1.0_dp
                data%alpha_static(i, j) = (data%eps_static(i, j) - kronecker) * factor
                data%alpha_ion(i, j) = data%eps_ion(i, j) * factor
                data%alpha_inf(i, j) = (data%eps_inf(i, j) - kronecker) * factor
            end do
        end do
    end subroutine compute_polarizability


    subroutine free_pol_data(data)
        !! Deallocate all allocatable arrays
        type(pol_data_t), intent(inout) :: data
        if (allocated(data%dipoles)) deallocate(data%dipoles)
        if (allocated(data%dipoles_clean)) deallocate(data%dipoles_clean)
        data%n_frames = 0
        data%n_unwraps = 0
        data%drift_rate = 0.0_dp
        data%eps_static = 0.0_dp
        data%eps_ion = 0.0_dp
        data%alpha_static = 0.0_dp
        data%alpha_ion = 0.0_dp
        data%alpha_inf = 0.0_dp
    end subroutine free_pol_data

end module polarizability
