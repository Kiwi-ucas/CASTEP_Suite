module dos_compute
    !! Density of states via Gaussian smearing from .bands eigenvalues
    !!   Total DOS: DOS(E) = Σ_k Σ_i w_k · G_σ(E − ε_i)
    !!   PDOS: PDOS_ch(E) = Σ_k Σ_i w_k · G_σ(E−ε_i) · W_ch(i,k)
    !!     where W_ch sums orbital weights per angular momentum channel
    !!   G_σ(dE) = (S / (σ √2π)) · exp(−½(dE/σ)²)
    !! S = 2 for non-spin-polarized (spin degeneracy), 1 for spin-polarized
    use castep_config, only: dp, pi, bands_data_t, pdos_data_t, HARTREE_TO_EV
    implicit none
    private

    public :: compute_total_dos, compute_pdos

    real(dp), parameter :: GAUSS_CUTOFF = 5.0_dp   ! |dE|/sigma beyond which exp() is negligible
    integer, parameter, public :: CH_TOT = 1, CH_S = 2, CH_P = 3, CH_D = 4, CH_F = 5
    integer, parameter, public :: N_CHANNELS = 5

contains

    subroutine compute_total_dos(bands, energy_grid, smearing, dos_result, iostat, iomsg)
        type(bands_data_t), intent(in) :: bands
        real(dp), intent(in) :: energy_grid(:)
        real(dp), intent(in) :: smearing
        real(dp), allocatable, intent(out) :: dos_result(:,:)
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        integer  :: ne, nk, nspin, nbands, ie, ik, is, ien
        real(dp) :: spin_coeff, norm, inv_two_s2, w_k, dE, fermi_e
        real(dp), allocatable :: eig_ev(:,:,:)

        iostat = 0
        if (present(iomsg)) iomsg = ''

        nk     = bands%num_kpoints
        nbands = bands%num_eigenvalues
        nspin  = bands%num_spin
        ne     = size(energy_grid)

        if (nk < 1 .or. nbands < 1) then
            iostat = 1
            if (present(iomsg)) iomsg = 'No band data for DOS'; return
        end if

        spin_coeff = 1.0_dp
        if (nspin == 1) spin_coeff = 2.0_dp
        norm = spin_coeff / (smearing * sqrt(2.0_dp * pi))
        inv_two_s2 = 0.5_dp / (smearing * smearing)

        fermi_e = bands%fermi_energy * HARTREE_TO_EV

        allocate(eig_ev(nbands, nk, nspin))
        do is = 1, nspin
            do ik = 1, nk
                do ie = 1, nbands
                    eig_ev(ie, ik, is) = bands%eigenvalues(ie, ik, is) &
                                         * HARTREE_TO_EV - fermi_e
                end do
            end do
        end do

        allocate(dos_result(ne, nspin), stat=iostat)
        if (iostat /= 0) then
            if (present(iomsg)) iomsg = 'Allocation failed'; deallocate(eig_ev); return
        end if
        dos_result = 0.0_dp

        do is = 1, nspin
            do ik = 1, nk
                w_k = bands%kpoint_coords(4, ik)
                do ie = 1, nbands
                    do ien = 1, ne
                        dE = energy_grid(ien) - eig_ev(ie, ik, is)
                        if (abs(dE) > GAUSS_CUTOFF * smearing) cycle
                        dos_result(ien, is) = dos_result(ien, is) &
                            + w_k * norm * exp(-dE*dE * inv_two_s2)
                    end do
                end do
            end do
        end do

        deallocate(eig_ev)
    end subroutine compute_total_dos

    ! ----------------------------------------------------------------
    !  Projected DOS: s, p, d, f channels from orbital weights
    ! ----------------------------------------------------------------
    subroutine compute_pdos(bands, pdos, energy_grid, smearing, pdos_result, iostat, iomsg)
        type(bands_data_t), intent(in) :: bands
        type(pdos_data_t),  intent(in) :: pdos
        real(dp), intent(in) :: energy_grid(:)
        real(dp), intent(in) :: smearing
        real(dp), allocatable, intent(out) :: pdos_result(:,:,:)  ! (ne, N_CHANNELS, nspin)
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        integer  :: ne, nk, nspin, nbands, norbs, nocc, ie, ik, is, jo, ich, ien
        real(dp) :: spin_coeff, norm, inv_two_s2, w_k, dE, fermi_e, gauss, ch_w
        real(dp), allocatable :: eig_ev(:,:,:)   ! (nbands, nk, nspin)
        integer,  allocatable :: ch_idx(:)        ! which channel for each orbital

        iostat = 0
        if (present(iomsg)) iomsg = ''

        nk     = bands%num_kpoints
        nbands = bands%num_eigenvalues
        nspin  = bands%num_spin
        ne     = size(energy_grid)
        norbs  = pdos%num_orbitals

        if (nk < 1 .or. nbands < 1) then
            iostat = 1; if (present(iomsg)) iomsg = 'No band data'; return
        end if
        if (norbs < 1) then
            iostat = 1; if (present(iomsg)) iomsg = 'No PDOS orbital data'; return
        end if

        ! map each orbital to its angular momentum channel (0=S,1=P,2=D,3=F → 1=s,2=p,3=d,4=f)
        allocate(ch_idx(norbs))
        do jo = 1, norbs
            ch_idx(jo) = pdos%orbital_am(jo) + 2   ! 0→2(s), 1→3(p), 2→4(d), 3→5(f)
        end do

        spin_coeff = 1.0_dp
        if (nspin == 1) spin_coeff = 2.0_dp
        norm = spin_coeff / (smearing * sqrt(2.0_dp * pi))
        inv_two_s2 = 0.5_dp / (smearing * smearing)

        fermi_e = bands%fermi_energy * HARTREE_TO_EV

        allocate(eig_ev(nbands, nk, nspin))
        do is = 1, nspin
            do ik = 1, nk
                do ie = 1, nbands
                    eig_ev(ie, ik, is) = bands%eigenvalues(ie, ik, is) &
                                         * HARTREE_TO_EV - fermi_e
                end do
            end do
        end do

        allocate(pdos_result(ne, N_CHANNELS, nspin), stat=iostat)
        if (iostat /= 0) then
            if (present(iomsg)) iomsg = 'Allocation failed for pdos_result'
            deallocate(eig_ev, ch_idx); return
        end if
        pdos_result = 0.0_dp

        nocc = min(nbands, pdos%max_bands)

        do is = 1, nspin
            do ik = 1, nk
                w_k = bands%kpoint_coords(4, ik)
                do ie = 1, nocc
                    do ien = 1, ne
                        dE = energy_grid(ien) - eig_ev(ie, ik, is)
                        if (abs(dE) > GAUSS_CUTOFF * smearing) cycle
                        gauss = w_k * norm * exp(-dE*dE * inv_two_s2)
                        pdos_result(ien, CH_TOT, is) = pdos_result(ien, CH_TOT, is) + gauss
                        do jo = 1, norbs
                            ich = ch_idx(jo)
                            ch_w = pdos%orbital_weights(jo, ie, ik, is)
                            pdos_result(ien, ich, is) = pdos_result(ien, ich, is) + gauss * ch_w
                        end do
                    end do
                end do
            end do
        end do

        deallocate(eig_ev, ch_idx)
    end subroutine compute_pdos

end module dos_compute
