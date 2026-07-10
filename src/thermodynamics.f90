module thermodynamics
    !! Compute thermodynamic properties from phonon frequencies.
    !!
    !! Follows CASTEP official formalism (Eqs. 80-84, Baroni et al. 2001):
    !!   E(T) = ZPE + Σ_q w_q Σ_i ℏω_i / (eˣ−1)       x = ℏω/kT
    !!   S(T) = k_B Σ_q w_q Σ_i [x/(eˣ−1) − ln(1−e⁻ˣ)]
    !!   Cv(T)= k_B Σ_q w_q Σ_i x² eˣ / (eˣ−1)²
    !!   F(T) = E(T) − T·S(T)
    !!
    !! This is the exact integral ∫F(ω)g(ω)dω with F(ω)=Σ w_q δ(ω−ω_i),
    !! i.e. the δ-function limit of the phonon DOS. Gaussian-smeared DOS
    !! introduces systematic error at low T (artificial ω≈0 weight from
    !! finite-σ tails), so we use direct summation which is exact.
    !!
    !! Acoustic modes at Gamma (3 translational d.o.f.) are excluded.
    use castep_config, only: dp
    implicit none
    private

    public :: thermo_data_t, compute_thermodynamics, free_thermo_data

    real(dp), parameter :: CM1_TO_MEV = 0.123984_dp
    real(dp), parameter :: KB_MEV     = 0.0861733_dp
    real(dp), parameter :: MEV_TO_J   = 1.602176634e-22_dp
    real(dp), parameter :: NA         = 6.02214076e23_dp

    type :: thermo_data_t
        integer :: n_temps = 0
        real(dp) :: zpe = 0.0_dp
        real(dp), allocatable :: temps(:)
        real(dp), allocatable :: energy(:)
        real(dp), allocatable :: free_e(:)
        real(dp), allocatable :: entropy(:)
        real(dp), allocatable :: heat_cap(:)
    end type thermo_data_t

contains

    subroutine compute_thermodynamics(phdos, t_min, t_max, n_pts, thermo, &
                                       iostat, iomsg, smearing)
        !! Direct summation over phonon modes with q-point weights.
        !! Equivalent to ∫F(ω)g(ω)dω with F(ω)=Σw_q δ(ω-ω_i).
        !! The `smearing` argument is accepted for interface compatibility
        !! but not used (δ-function limit).
        use phonon_dos, only: phonon_dos_data_t
        type(phonon_dos_data_t), intent(in) :: phdos
        real(dp), intent(in) :: t_min, t_max
        integer,  intent(in) :: n_pts
        type(thermo_data_t), intent(out) :: thermo
        integer,  intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg
        real(dp), intent(in), optional :: smearing   ! accepted, not used

        integer  :: i, j, n_modes, nb, i_mode, idx
        real(dp) :: t, dt, freq_mev, x, expx, expx_m1, ln_term
        real(dp) :: e_sum, s_sum, cv_sum, zpe_accum, w

        iostat = 0

        if (.not. allocated(phdos%freqs) .or. .not. allocated(phdos%weights)) then
            iostat = 1
            if (present(iomsg)) iomsg = 'Phonon data not loaded'
            return
        end if

        n_modes = size(phdos%freqs)
        nb      = phdos%n_branches
        if (n_modes == 0 .or. nb <= 0) then
            iostat = 2
            if (present(iomsg)) iomsg = 'Empty or invalid phonon data'
            return
        end if

        ! ── Zero-point energy: ½ Σ w_q Σ_i ℏω_i  (skip acoustic at Gamma) ──
        zpe_accum = 0.0_dp
        do i = 1, n_modes
            j = (i - 1) / nb + 1          ! q-point index (1..n_q)
            i_mode = mod(i - 1, nb) + 1    ! mode index within q-point
            if (j == 1 .and. i_mode <= 3) cycle  ! skip Gamma acoustic
            if (phdos%freqs(i) <= 0.0_dp) cycle   ! skip negative (ASR artifact)
            freq_mev = phdos%freqs(i) * CM1_TO_MEV
            zpe_accum = zpe_accum + 0.5_dp * freq_mev * phdos%weights(j)
        end do
        thermo%zpe = zpe_accum

        ! ── Temperature loop ──
        thermo%n_temps = n_pts
        allocate(thermo%temps(n_pts), thermo%energy(n_pts), thermo%free_e(n_pts), &
                 thermo%entropy(n_pts), thermo%heat_cap(n_pts))

        dt = (t_max - t_min) / max(n_pts - 1, 1)

        do idx = 1, n_pts
            t = t_min + (idx - 1) * dt
            if (t < 1.0e-6_dp) t = 1.0e-6_dp

            thermo%temps(idx) = t
            e_sum  = 0.0_dp
            s_sum  = 0.0_dp
            cv_sum = 0.0_dp

            do i = 1, n_modes
                j = (i - 1) / nb + 1
                i_mode = mod(i - 1, nb) + 1
                if (j == 1 .and. i_mode <= 3) cycle
                if (phdos%freqs(i) <= 0.0_dp) cycle
                w = phdos%weights(j)

                freq_mev = phdos%freqs(i) * CM1_TO_MEV
                x = freq_mev / (KB_MEV * t)
                if (x > 50.0_dp) cycle

                expx    = exp(x)
                expx_m1 = expx - 1.0_dp

                ! E: ℏω / (eˣ−1)
                e_sum = e_sum + freq_mev / expx_m1 * w

                ! S: k [ x/(eˣ−1) − ln(1−e⁻ˣ) ]
                if (expx_m1 > 1.0e-30_dp) then
                    ln_term = log(1.0_dp - exp(-x))
                    s_sum = s_sum + (x / expx_m1 - ln_term) * KB_MEV * w
                end if

                ! Cv: k x² eˣ / (eˣ−1)²
                cv_sum = cv_sum + KB_MEV * x * x * expx / (expx_m1 * expx_m1) * w
            end do

            thermo%energy(idx)  = thermo%zpe + e_sum
            thermo%free_e(idx)  = thermo%zpe + e_sum - t * s_sum
            ! Convert meV/K → J/mol/K: × MEV_TO_J × NA ≈ × 96.5
            thermo%entropy(idx) = s_sum * MEV_TO_J * NA
            thermo%heat_cap(idx)= cv_sum * MEV_TO_J * NA
        end do

    end subroutine compute_thermodynamics


    subroutine free_thermo_data(thermo)
        type(thermo_data_t), intent(inout) :: thermo
        if (allocated(thermo%temps))    deallocate(thermo%temps)
        if (allocated(thermo%energy))   deallocate(thermo%energy)
        if (allocated(thermo%free_e))   deallocate(thermo%free_e)
        if (allocated(thermo%entropy))  deallocate(thermo%entropy)
        if (allocated(thermo%heat_cap)) deallocate(thermo%heat_cap)
        thermo%n_temps = 0
    end subroutine free_thermo_data

end module thermodynamics
