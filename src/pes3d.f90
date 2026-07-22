module pes3d
    !! Deprecated wrapper — use the unified `pes` module directly.
    !! Re-exports the 3D-relevant symbols from `pes` for backward compatibility.
    use pes, only: pes_grid_t, compute_local_grid_bounds, &
        generate_pes_grid_points, write_pes_cube, write_pes3d_cube, &
        collect_pes_energies, symmetry_expand_energies
    implicit none
    private
    ! Re-export with old names for compatibility
    public :: pes3d_grid_t, compute_local_grid_bounds, &
        generate_pes3d_grid_points, write_pes3d_cube, &
        collect_pes3d_energies, symmetry_expand_energies

    ! Type alias: pes3d_grid_t is the same as pes_grid_t
    type, extends(pes_grid_t) :: pes3d_grid_t
    end type pes3d_grid_t

contains

    subroutine generate_pes3d_grid_points(grid, frac_points, n_total, iostat, iomsg)
        !! Wrapper: calls unified generate_pes_grid_points with ndim=3.
        use pes, only: pes_grid_t
        class(pes3d_grid_t), intent(in) :: grid
        real(dp), allocatable, intent(out) :: frac_points(:,:)
        integer, intent(out) :: n_total, iostat
        character(len=*), optional, intent(out) :: iomsg
        type(pes_grid_t) :: ugrid
        ugrid%ndim = 3
        ugrid%n_points = grid%n_points
        ugrid%frac_range = grid%frac_range
        ugrid%mobile_atom_idx = grid%mobile_atom_idx
        ugrid%scan_mode = grid%scan_mode
        ugrid%use_symmetry = grid%use_symmetry
        ugrid%ref_frac = grid%ref_frac
        ugrid%half_dist = grid%half_dist
        call generate_pes_grid_points(ugrid, frac_points, n_total, iostat, iomsg)
    end subroutine generate_pes3d_grid_points

    subroutine collect_pes3d_energies(scan_dir, iostat, iomsg)
        !! Wrapper: calls unified collect_pes_energies.
        character(len=*), intent(in) :: scan_dir
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg
        call collect_pes_energies(scan_dir, iostat, iomsg)
    end subroutine collect_pes3d_energies

end module pes3d
