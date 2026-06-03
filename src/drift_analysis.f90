module drift_analysis
    !! Drift rate diagnostics for anisotropic diffusion analysis
    !! Extracted from polarizability module for future development.
    !! NOT compiled into the main program — saved as standalone analysis tool.
    !!
    !! Potential application: analyze directional Li+ diffusion rates in
    !! superionic conductors to characterize anisotropic transport.
    !!
    !! Usage (when integrated):
    !!   use drift_analysis, only: compute_drift_rates, compute_global_dielectric
    !!   call compute_drift_rates(pol, time_step_ps)
    !!   call compute_global_dielectric(pol, iostat, iomsg)
    use castep_config, only: dp
    use polarizability, only: pol_data_t, EPSILON_0, KBOLTZMANN, DEBYE_TO_CM, ANG3_TO_M3
    implicit none
    private

    public :: compute_drift_rates, compute_global_dielectric

contains

    subroutine compute_drift_rates(data, time_step_ps)
        !! Compute per-direction drift rates from linear fit to unwrapped dipoles.
        !! Drift rate = slope of linear fit (D/frame) / time_step_ps → D/ps.
        !! Large drift rates indicate preferred diffusion directions.
        type(pol_data_t), intent(inout) :: data
        real(dp), intent(in) :: time_step_ps  ! MD time step in ps
        real(dp) :: sum_t, sum_t2, sum_d(3), sum_td(3)
        real(dp) :: denom, s, b, n
        integer :: i, j

        n = real(data%n_frames, dp)
        sum_t = 0.0_dp; sum_t2 = 0.0_dp
        sum_d = 0.0_dp; sum_td = 0.0_dp

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

        do j = 1, 3
            s = (n * sum_td(j) - sum_t * sum_d(j)) / denom
            data%drift_rate(j) = s / time_step_ps
        end do
    end subroutine compute_drift_rates


    subroutine compute_global_dielectric(data, iostat, iomsg)
        !! Compute ionic dielectric tensor from global dipole covariance.
        !! ε_ion,αβ = Cov(M_α, M_β) / (ε₀ k_B T Ω)
        !!
        !! NOTE: This uses the full trajectory and includes residual diffusion
        !! contamination. For the physically meaningful vibrational contribution,
        !! use compute_static_dielectric_windowed from the polarizability module.
        !!
        !! This subroutine is useful for anisotropy analysis:
        !! compare ε_ion,xx / ε_ion,yy / ε_ion,zz to identify preferred
        !! ionic polarization directions.
        type(pol_data_t), intent(inout) :: data
        integer, intent(out)  :: iostat
        character(len=*), intent(out), optional :: iomsg
        real(dp) :: m_sum(3), m2_sum(3,3), m_mean(3), cov(3,3)
        real(dp) :: vol_m3, kt, denom
        integer :: i, j, k, n

        iostat = 0
        n = data%n_frames

        if (n < 2) then
            iostat = -1
            if (present(iomsg)) iomsg = 'Need at least 2 frames for statistics'
            return
        end if

        if (data%temperature <= 0.0_dp) then
            iostat = -1
            if (present(iomsg)) iomsg = 'Temperature must be positive'
            return
        end if

        if (data%volume_ang3 <= 0.0_dp) then
            iostat = -1
            if (present(iomsg)) iomsg = 'Cell volume must be positive'
            return
        end if

        vol_m3 = data%volume_ang3 * ANG3_TO_M3
        kt = KBOLTZMANN * data%temperature
        denom = EPSILON_0 * kt * vol_m3

        m_sum = 0.0_dp
        do i = 1, n
            do j = 1, 3
                m_sum(j) = m_sum(j) + data%dipoles_clean(i, j)
            end do
        end do
        m_mean = m_sum / real(n, dp)

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
        cov = cov * DEBYE_TO_CM**2
        data%eps_ion = cov / denom
    end subroutine compute_global_dielectric

end module drift_analysis
