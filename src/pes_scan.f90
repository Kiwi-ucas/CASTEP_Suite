module pes_scan
    !! Deprecated wrapper — use the unified `pes` module directly.
    !! Re-exports the 2D-relevant symbols from `pes` for backward compatibility.
    use pes, only: pes_grid_t, generate_pes_grid_points, &
        write_pes_metadata_json, collect_pes_energies
    implicit none
    public :: pes_grid_t, generate_pes_grid_points, &
        write_pes_metadata_json, collect_pes_energies
end module pes_scan
