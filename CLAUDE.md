# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CASTEP Suite is a Fortran 2018 CLI suite with a Rust/Bevy 3D crystal structure viewer:

1. **PreCASTEP** — converts crystallographic structure files (CIF, PDB, .cell) into CASTEP DFT input files (`.cell` and `.param`)
2. **PosCASTEP** — post-processes CASTEP output: band structure with gap analysis, total DOS, projected DOS (s/p/d/f), phonon DOS, IR spectrum, Raman spectrum, static polarizability, PES (potential energy surface) scan with 3D visualization, and 3D crystal structure viewing. All plots ASCII terminal or CSV export
3. **crystal-viewer** — standalone Rust/Bevy 3D viewer for interactive crystal structure visualization, atom picking, and editing (launched from PosCASTEP -1 or standalone)

Launching `./CASTEP_Suite` shows a top-level suite menu. The PreCASTEP mode exits after successful file generation. PosCASTEP returns to its sub-menu (Q. Back) and then to the suite menu (Q. Quit). No `stop` statements remain — all error/recovery paths use `return`.

## Build

```bash
make          # Release build (gfortran -O2 + cargo --release)
make debug    # Debug build (-O0, -fcheck=all)
make run      # Build and run
make clean    # Remove obj/, CASTEP_Suite, and crystal-viewer/target/
```

Requirements: gfortran >= 7.0, cargo/rustc >= 1.75, Make.

If cargo is not installed, `make` builds Fortran only and skips the viewer.

## Running

```bash
./CASTEP_Suite
```

### Suite top-level menu (main.f90)

```
  ==================================
         CASTEP Suite
  ==================================
  1. PreCASTEP  (generate CASTEP input files)
  2. PosCASTEP  (post-process CASTEP output)
  Q. Quit
```

### PreCASTEP (option 1)

Interactive menu prompts for CIF/PDB/.cell file path via `ask_input_file`, auto-generates output base name from input stem + task type, then loops with numbered options for CASTEP parameters:
- `-2. Advanced option` — Smearing (off/on), Max SCF cycles, Convergence window (min 2), Calculate ELF, Calculate EDD
- `-1. Spin_polarized : true/false` — Toggle spin polarization; auto-set to true when pseudopotential is SOC19
- `0` — Generate `.cell` and `.param` files and **exit the program** (does not return to suite menu)
- `1` — Configure task type (17 options, some marked "暂未开发")
- `2` — Configure XC functional (5 options)
- `3` — Configure cutoff energy (eV)
- `4` — Configure vdW correction
- `5` — Configure pseudopotential
- `6` — Configure K-point scheme
- `7` — Configure SCF tolerance
- `8` — Configure symmetry handling (NONE/AUTO)
- `9` — Spectral Task (only shown when task is ElectronicSpectroscopy; sub-menu: BandStructure / BandStructure_pDOS)
- `9/10/11` — Optimizer, Cell opt mode, Geo tolerance (only shown when task is GEOMETRYOPTIMISATION)
- `9-14` — CINEB parameters: max images, spring constant, tangent mode, NEB method, max iterations, TS tolerance (only shown when task is CINEB)
- `9-14` — Phonon q-point scheme, method, fine method, energy tol, supercell matrix, fine q-point (only shown for phonon/EFIELD/thermo tasks)
- `Q. Back` — Return to suite menu (returns `IO_USER_QUIT`)

### PosCASTEP (option 2)

Post-processing menu for CASTEP output analysis and structure visualization:
- `-1. View Crystal Structure (3D)` — prompts for CIF/PDB/.cell file → auto-generates JSON → launches crystal-viewer (3D interactive rendering, atom picking, editing with 6-direction movement). If the structure was modified in the viewer, shows a sub-menu:
  - `1. Save new structure` — save as CIF/PDB/cell to file
  - `2. Save new structure and PreCASTEP` — pass modified structure directly to PreCASTEP (no intermediate file), then configure and generate .cell/.param
- `1. Plot Band Structure` — prompts for `.bands` file path, interactive ASCII terminal plot
- `2. Plot DOS` — total density of states from `.bands` only; ASCII (interactive) or CSV export
- `3. Plot pDOS` — projected DOS: enter file prefix (e.g. `Cu`), auto-loads `<prefix>.bands` + `<prefix>.pdos_bin`; ASCII (interactive) or CSV export with s/p/d/f columns
- `4. Plot Phonon DOS` — phonon density of states from `.phonon` file; ASCII (interactive) or CSV export
- `5. Plot IR Spectrum` — infrared absorption spectrum from `.phonon` file; ASCII (interactive) or CSV export
- `6. Plot Raman Spectrum` — Raman scattering spectrum from `.phonon` file; ASCII (interactive) or CSV export
- `7. Static Polarizability` — AIMD polarization fluctuation method: prompts for CASTEP `.castep` file (ε_∞), CP2K dipole directory, cell parameters, temperature, and MD time step; computes ionic dielectric constant via window-based extrapolation to W→0
- `8. Thermodynamics` — compute E(T), S(T), F(T), Cv(T) from `.phonon` file via direct summation over phonon modes (CASTEP/Baroni formalism, δ-function limit). User-configurable temperature range + number of points; interactive ASCII plot (4 merged curves: E_vib, F_vib, TS, Cv) with ↑↓/←→/+−/R controls; CSV export
- `9. PES Scan` — 2D/3D potential energy surface scan: sub-menu (1. generate 2D grid with IONIC_CONSTRAINTS, 2. collect 2D .castep energies → scan.cube → viewer, 3. generate 3D voxel grid, 4. collect 3D results → symmetry expansion). Output is Gaussian Cube format (`.cube`) for both 2D (`nz=1`) and 3D; line-2 JSON metadata carries plane, mobile atom, lattice, and scan ranges. **3D symmetry mode**: user N is internally mapped to N−1 so that the expanded cube header displays exactly the user's input (orbit grid with N−1 points → cube header `(N-1)+1 = N`).
- `10. Frozen Phonon Scan` — `.phonon` → Γ-only validation → mode list (index/frequency) → user selects mode → raw `e/√m` normalized to max|u|=1 → ±0.1/±0.2/±0.5 Å Cartesian displacements converted to fractional and wrapped → preview → PreCASTEP with task locked to `SINGLEPOINT` (`run_main_menu(..., lock_task=.true.)`) → writes six `frozen_phonon/<stem>_mode<M>/d±A/` directories each with `.cell`+`.param`
- `-2. Phonon Mode Visualization` — parse eigenvectors + Born charges from `.phonon`/`.castep`, decompose modes, launch viewer with displacement arrows

At any file path prompt, type `q` to cancel and return to the PosCASTEP menu. All interactive ASCII plots use the alternate screen buffer (`ESC[?1049h/l`) to avoid polluting terminal scrollback history.

**Band structure analysis** (bands_plotter.f90):
- Auto-detects terminal size via `stty size` (queries terminal driver directly, works even in `-icanon` mode), falls back to `COLUMNS`/`LINES` env vars, default 160×40
- Re-detects terminal size on every redraw — adapts to window resize events
- Detects VBM (highest occupied) and CBM (lowest unoccupied) across all k-points
- Marks VBM with green ◆ and CBM with yellow ◈ in the plot
- Calculates band gap size and determines direct vs indirect type
- Gap info (E_F, Band Gap, VBM, CBM, Dir gap) displayed **above** the plot in 2 compact lines (Dir gap merged into VBM/CBM line when indirect)
- Grid height = terminal rows − 9 (reserved for headers + footer), clamped to min 15
- Energy range: Fermi level ±10 eV, scrollable with ↑↓ keys, zoomable with +/-
- K-path: horizontal scrolling with ← → keys (first press auto-zooms to 50% window, subsequent presses pan within the path)
- R key resets energy range, k-path window, and k-path zoom to defaults
- Grid fills terminal: `label_width + 1(gap) + nw`, label-plot spacing = 1 column
- Raw terminal mode uses `stty -icanon -echo min 1` (preserves `opost` for correct newline translation, avoiding staircase alignment bug)
- All energies in eV only (Hartree units removed)
- SVG output similarly uses eV and Fermi ±10 eV range
- ANSI color codes are module-level parameters (`C_RED`, `C_GREEN`, `C_YELLOW`, `C_CYAN`, `C_BOLD`, `C_DIM`, `C_RESET`, `C_AXIS`)

**DOS/PDOS analysis** (dos_compute.f90, dos_plotter.f90, pdos_parser.f90):
- Total DOS via Gaussian smearing: `DOS(E) = Σ_k Σ_i w_k · G_σ(E − ε_i)`, smearing width configurable (default 0.1 eV)
- PDOS via angular-momentum-resolved orbital weights from `.pdos_bin`/`.pdos_weights` binary files (big-endian record-delimited format)
- PDOS channels: s (● cyan), p (○ yellow), d (△ green), f (▽ red) with legend line below plot
- Interactive ASCII controls: ↑↓ y-axis pan, ←→ x-axis pan, +/- overall zoom (both axes simultaneously), R reset all, Q quit
- y-axis uses y_center + y_half model symmetric with x-axis (e_center + e_half); initial range [0, y_max0]
- y=0 horizontal reference line (├───┤) visible when y-axis panned below zero
- Grid height = terminal rows − 9 (DOS) or − 10 (pDOS, legend)
- CSV export: Energy,Total,s,p,d,f columns for pDOS; Energy,DOS[_up,_down] for total DOS
- Path memory: `.bands` path and pDOS prefix remembered across repeated calls via SAVE
- Auto `.csv` extension: user enters stem name, `.csv` appended automatically (no double-append)
- Energy grid: Fermi ±20 eV precomputed at 4001 points, interactive viewport shows subset

## Architecture

Twenty-two Fortran source files in `src/` (21 compiled, 1 uncompiled development artifact), plus a Rust subproject:

1. **config.f90** — `castep_config` module. Defines kinds (`dp`), physical constants (`HARTREE_TO_EV`), constants (16 task types including `TASK_TRANSITION_STATE='TRANSITIONSTATESEARCH'` for CINEB, 7 CINEB constants for tangent mode + NEB method, 5 XC functionals: PBE/PBEsol/HSE06/PBE0/r2scan, 4 vdW corrections, 3 pseudopotentials, 3 K-point schemes, 3 optimizers, 4 geo tolerance levels, 3 cell opt modes, 2 symmetry modes, phonon constant groups, I/O error codes), max sizes (MAX_ATOMS=10000, MAX_TAGS=5000, MAX_LOOP_ROWS=50000, MAX_SYM_OPS=400), CIF/CASTEP tag names, and types (`atom_t`, `cif_data_t`, `castep_config_t`, `bands_data_t`). Allocatable fields: `cif_data_t%atoms`, `bands_data_t%kpoint_indices`, `bands_data_t%kpoint_coords`, `bands_data_t%kpath_dist`, `bands_data_t%eigenvalues`. CINEB fields in `castep_config_t`: product/intermediate atom arrays (`prod_atom_type/x/y/z`, `interm_atom_type/x/y/z`), CINEB parameters (`cineb_max_images`, `cineb_spring_constant`, `cineb_max_iter`, `cineb_tangent_mode`, `cineb_neb_method`, `cineb_climbing` — all `character(16)` strings; `ts_geom_tolerance`). `scf_tolerance` is `real(dp)` (was `character(32)`). `castep_config_t` has a `final :: finalize_castep_config` procedure that auto-deallocates reactant, product, and intermediate atom arrays on scope exit. Error codes: `IO_FILE_NOT_FOUND=100`…`IO_USER_QUIT=-1`, `IO_PRECASTEP_LAUNCH=-2`. Physical constants and polarizability error codes moved to `polarizability` module. Provides `default_config`, `new_castep_config`, tag normalization (`normalize_tag`, `compare_tags`), `string_to_real`, `int2str`, `get_castep_task_name`, `strip_quotes`, `compute_cartesian_lattice` (shared lattice computation merged from main.f90/poscastep_menu.f90).

2. **term_utils.f90** — `term_utils` module. Leaf module (zero dependencies). Provides shared ANSI color constants (`C_RED`, `C_GREEN`, `C_YELLOW`, `C_CYAN`, `C_BOLD`, `C_DIM`, `C_RESET`, `C_AXIS`), alternate screen buffer constants (`C_ALT_ON`/`C_ALT_OFF`), terminal size detection (`get_term_size` + private `stty_size`), Bresenham line-drawing (`draw_line`), and raw terminal mode management (`enter_raw_mode`/`leave_raw_mode` — encapsulate `stty -icanon -echo min 1` + alternate screen buffer on/off). Used by `bands_plotter`, `dos_plotter`, and `poscastep_menu`.

3. **parser.f90** — `parser` module. File format parsers: `parse_cif_inline`, `parse_pdb_inline`, `parse_cell_inline`. Private helpers: `tokenize_inline`, `clean_str_inline`, `copy_str_no_quotes`. Exports `clean_element_symbol` for oxidation state stripping (e.g., "Cu0+" -> "Cu"). Handles CIF tag-value pairs and `loop_` blocks, PDB CRYST1/UNITCELL/ATOM records, CASTEP .cell `%BLOCK LATTICE_ABC`, `%BLOCK LATTICE_CART`, `%BLOCK POSITIONS_ABS`, and `%BLOCK POSITIONS_FRAC`. `parse_cell_inline` includes `compute_abc_from_cartesian` for Cartesian-to-lattice conversion. `cif_data_t` has a `positions_fractional` field set by the parser — CIF and `POSITIONS_FRAC` set it `.true.`, PDB and `POSITIONS_ABS` leave it `.false.`.

4. **cell_writer.f90** — `cell_writer` module. Generates CASTEP `.cell` file in `%BLOCK` format. Each `%BLOCK` type written by independent subroutine: `write_block_lattice_abc`, `write_block_species_pot`, `write_block_positions_abs`, `write_block_cell_constraints`, `write_block_kpoint_grid`. CINEB blocks: `write_block_positions_abs_product`, `write_block_positions_abs_intermediate` (written when task is TRANSITIONSTATESEARCH). Phonon-specific blocks: `write_block_phonon_kpoint_mp`, `write_block_phonon_kpoint_path`, `write_block_phonon_fine_kpoint_path`, `write_block_phonon_supercell_matrix`. `write_block_positions_abs` handles both fractional (CIF) and Cartesian (PDB/.cell) input. `write_block_cell_constraints` only written when task is GEOMETRYOPTIMISATION and cell_opt_mode is ALL. `SYMMETRY_GENERATE` written when sym_source is AUTO. `KPOINTS_MP_GRID` written for GAMMA and MONKHORST_PACK schemes. `PHONON_KPOINT_MP_GRID` and `PHONON_FINE_KPOINT_MP_GRID` written for phonon tasks. Private helper `is_phonon_task` gates phonon blocks.

5. **param_writer.f90** — `param_writer` module. Generates CASTEP `.param` file in `key : value` format. Task line first, then common keywords, then task-diff keywords (16 task-specific blocks including full CINEB/PHONON/PHONON+EFIELD/THERMODYNAMICS support). CINEB params written by `write_cineb_params` (tssearch_method, max_path_points, spring_constant, tangent_mode, neb_method, max_iter, climbing, plus TS tolerance select-case). SCF tolerance written via `energy_tol_str` helper. `elec_energy_tol` unit is implicit (eV). Phonon keywords written by shared `write_phonon_params` helper (method, energy_tol, max_cycles, dfpt_method, DOS, sum_rule, output control, Raman, etc.). `calculate_raman` also writes `raman_range_low : 0.00e+00 cm-1` and `raman_range_high : 1.00e+04 cm-1`. `calculate_born_charges` also writes `born_charge_sum_rule : true`. `write_efield_params` writes `efield_ignore_molec_modes` (CRYSTAL/MOLECULE/LINEAR_MOLECULE). vdW-DED keywords when vdW is not NONE. Force constant cutoff branches on method: SPHERICAL→`phonon_force_constant_cutoff`, CUMULANT→`phonon_force_constant_cutoff_scale`. Each keyword written by `write_kv(unit, key, value)`. CASTEP 25.12 compatibility: no `_unit` suffix keywords, no explicit toggle keywords (`thermo`, `electric_field`) — units are implicit and tasks self-enable.

6. **bands_parser.f90** — `bands_parser` module. Parses CASTEP `.bands` output files into `bands_data_t`. Public: `parse_bands_file(filename, bands, iostat, iomsg)`, `free_bands_data(bands)`. Reads header metadata (num_kpoints, num_spin, num_electrons, num_eigenvalues, fermi_energy), k-point coordinates with cumulative path distances, and eigenvalue data. Eigenvalues stored in Hartree as `eigenvalues(ie, ik, is)`.

7. **bands_plotter.f90** — `bands_plotter` module. Band structure visualization with **gap-focused analysis**. Public constants: `BANDS_MODE_ASCII=1`. Public subroutines: `plot_bands_ascii(bands, term_w_in, term_h_in, e_center, half_range, k_pct_in, k_width_pct_in)` — interactive scrollable ANSI-color terminal plot with k-path windowing and gap info above plot. Uses `term_utils` for ANSI colors, terminal detection, and line drawing. Private: `print_grid_row`, `char_type`, `type_color`, `type_char`. Features: symmetric band display (first k-point copied to end), Unicode box frame, multi-symbol bands, k-point path labels (auto-detected K1..Kn), right-aligned energy labels with dynamic width, gap info in 2 lines above plot, dynamic grid height (nh−9) to fit terminal. **Gap analysis**: scans all eigenvalues in eV, finds VBM (highest below E_F) and CBM (lowest above E_F), computes indirect gap (VBM→CBM) and direct gap at VBM k-point, classifies as direct/indirect/metallic (gap < 0.005 eV). VBM marker: green ◆, CBM marker: yellow ◈. Optional `k_pct_in`/`k_width_pct_in` parameters define a k-path window for horizontal scrolling.

8. **pdos_parser.f90** — `pdos_parser` module. Parses CASTEP `.pdos_weights` / `.pdos_bin` binary files into `pdos_data_t`. Binary format: big-endian, record-delimited (4B u32 size prefix + data + 4B u32 suffix). `.pdos_bin` has two extra prefix records (f64 version + version string). Uses `access='stream', form='unformatted'` for byte-level I/O with explicit big-endian→little-endian conversion. Public: `parse_pdos_file(filename, pdos, iostat, iomsg)`, `free_pdos_data(pdos)`. Private: `read_be_u32`, `read_be_f64`, `read_payload`, `skip_record`, `be_u32_at`, `be_f64_at`, `read_u32_record`, `read_int_vec`. Output: `pdos_data_t` with header fields (total_kpoints, num_spins, num_orbitals, max_bands), orbital metadata arrays (species, ion, am[0=S,1=P,2=D,3=F]), k-point coordinates, and 4D orbital_weights(norbs, nbands, nk, nspin).

9. **dos_compute.f90** — `dos_compute` module. Gaussian smearing DOS computation. Total DOS: `DOS(E) = Σ_k Σ_i w_k · G_σ(E − ε_i)`. PDOS: `PDOS_ch(E) = Σ_k Σ_i w_k · G_σ(E−ε_i) · W_ch(i,k)` where W_ch sums orbital weights by angular momentum channel. Spin coefficient: 2 for non-polarized, 1 for spin-polarized. Gaussian cutoff at 5σ for performance. Uses `pi` from `castep_config`. Public: `compute_total_dos(bands, energy_grid, smearing, dos_result, iostat, iomsg)` — allocates `dos_result(ne, nspin)`, `compute_pdos(bands, pdos, energy_grid, smearing, pdos_result, iostat, iomsg)` — allocates `pdos_result(ne, N_CHANNELS, nspin)`. Channel constants: `CH_TOT=1, CH_S=2, CH_P=3, CH_D=4, CH_F=5`, `N_CHANNELS=5`.

10. **dos_plotter.f90** — `dos_plotter` module. DOS/PDOS visualization: interactive ASCII terminal plots + CSV export. Public constants: `DOS_MODE_ASCII=1`, `DOS_MODE_EXPORT=3`. Public: `plot_dos_ascii(energy_grid, dos_data, nspin, e_fermi, smearing, term_w_in, term_h_in, y_center_in, y_half_in, e_center_in, half_range_in)` — interactive total DOS ASCII plot with y_center/y_half axis model, `plot_pdos_ascii` — s/p/d/f multi-channel PDOS plot with legend (● s, ○ p, △ d, ▽ f), `write_dos_csv`, `write_pdos_csv`. Uses `term_utils` for ANSI colors, terminal detection, and line drawing. y=0 horizontal reference line with ├┤ junctions. PDOS character types extended: 'S'/'P'/'L'/'F' for s/p/d/f orbitals, 'Y'/'Z' for y=0 junctions.

11. **cli_menu.f90** — `cli_menu` module. PreCASTEP main configuration menu loop with cached `castep_config_t` state. Q returns `IO_USER_QUIT`. Optional `lock_task` dummy locks menu item 1 (used by frozen-phonon scan; displays `[locked]` and refuses task changes). Task switch auto-configures DFPT→NCP19, THERMO→FD+SUPERCELL, PHONON+EFIELD→LO/TO ON. Uses `strip_quotes` from config. Helper functions: `task_label`, `cutoff_label`, `kpoint_label`, `qpoint_label`, `geom_tol_label`, `sym_label`, `sp_label`, `smearing_label`, `sci_str` (scientific notation), `real_str` (formatted real), `spectral_label` (maps internal constant to display label). **Unified input**: `ask_positive_real(prompt, default, result, iostat)` — single subroutine for all real-number inputs with do-loop retry and error messages, replacing 4 separate ask subroutines (`ask_cutoff_energy`, `ask_scf_tolerance`, `ask_phonon_energy_tol` deleted) and 7 inline patterns in advanced options. CINEB-specific: `ask_cineb_max_images` (odd-enforcing), `ask_cineb_spring_constant` (wraps `ask_positive_real`), `ask_cineb_tangent_mode`, `ask_cineb_neb_method`, `ask_cineb_max_iter`.

12. **poscastep_menu.f90** — `poscastep_menu` module. PosCASTEP post-processing menu loop. Public: `run_poscastep_menu(iostat)`. Menu options: -2. Phonon Mode Visualization, -1. View Crystal Structure (3D), 0. Format Converter, 1. Plot Band Structure, 2. Plot DOS, 3. Plot pDOS, 4. Plot Phonon DOS, 5. Plot IR Spectrum, 6. Plot Raman Spectrum, 7. Static Polarizability, 8. Thermodynamics, 9. PES Scan, 10. Frozen Phonon Scan, Q. Back. Frozen phonon helpers: `handle_frozen_phonon_menu` (prompt .phonon → Γ-only check → mode/frequency list → raw `e/√m` normalization to max|u|=1 → ±0.1/±0.2/±0.5 Å preview → PreCASTEP with `lock_task=.true.` → six `frozen_phonon/<stem>_mode<M>/d±A/` .cell/.param sets), `frozen_phonon_preview`, `frozen_phonon_populate_cfg`, `frozen_phonon_write_inputs`, `frozen_phonon_amp_label`. Uses `IO_USER_QUIT` for consistent 'q' handling at file input prompts — returns to PosCASTEP menu. Private: `handle_bands_menu` → `run_ascii_navigator` (↑↓ energy scroll, ← → k-path, +/- zoom both axes, R reset, Q quit), `handle_dos_menu` — total DOS (prompts .bands path with SAVE memory; output: ASCII/CSV; 4001-point Fermi ±20 eV grid), `handle_pdos_menu` — projected DOS (prompts file prefix, auto-derives .bands + .pdos_bin paths with SAVE memory; output: ASCII/CSV), `handle_phonon_dos_menu` — phonon DOS from `.phonon` file (ASCII/CSV), `handle_ir_menu` — IR spectrum (ASCII/CSV), `handle_raman_menu` — Raman spectrum (ASCII/CSV), `handle_phonon_modes_menu` — phonon mode visualization (parse eigenvectors + Born charges, decompose modes, write displacement JSON, launch viewer), `run_dos_navigator` (↑↓ y-pan, ←→ x-pan, +/- both-axes zoom, R reset, Q quit), `run_pdos_navigator` (same controls), `run_phonon_dos_navigator`, `run_ir_navigator`, `run_raman_navigator`, `build_energy_grid`, `ensure_ext` (auto-append file extension). Uses `enter_raw_mode`/`leave_raw_mode` from `term_utils` for terminal mode + alternate screen buffer.

13. **phonon_dos.f90** — `phonon_dos` module. Parses CASTEP `.phonon` files (frequencies, IR intensities, Raman activities). Computes phonon DOS via Gaussian smearing. Computes IR absorption and Raman scattering spectra. Public: `parse_phonon_file`, `compute_phonon_dos`, `compute_ir_spectrum`, `compute_raman_spectrum`, `free_phonon_dos_data`. Output types: `phonon_dos_data_t` with frequency grid, phonon DOS, IR spectrum, and Raman spectrum arrays.

14. **polarizability.f90** — `polarizability` module. Static polarizability via AIMD polarization fluctuation method. Defines `pol_data_t` type (ε_∞ tensor, dipole trajectory, result tensors, unwrap statistics). Public: `parse_castep_epsilon` — extract optical dielectric tensor from `.castep` file, `parse_cp2k_dipoles` — read CP2K Berry phase dipole files via `find | sort` (handles 10k+ files without ARG_MAX overflow), `unwrap_dipoles` — unwrap polarization quantum jumps (raw-to-raw comparison with cumulative offset, nint for multi-quantum), `compute_static_dielectric_windowed` — window-based ε_ion with per-window detrend + median + W→0 extrapolation, `compute_polarizability` — convert ε tensor to α (ų). Private helpers: `detrend_window` (linear detrend), `compute_covariance_single` (3×3 covariance), `median` (selection sort). Module-level constants: `EPSILON_0`, `KBOLTZMANN`, `DEBYE_TO_CM`, `ANG3_TO_M3`, `DEBYE_PER_ANG`. Error codes: `IO_EPS_NOT_FOUND=112`, `IO_EPS_PARSE_ERROR=113`, `IO_DIPOLE_ERROR=114`.

15. **pes.f90** — `pes` module. Unified PES (Potential Energy Surface) scan. Defines `pes_grid_t` type with `ndim` (2/3), `n_irred`, `irred_coords` (allocatable, for orbit mapping). Public: `generate_pes_grid_points` — dispatches to 2D/3D grid generation; `get_irreducible_grid` — **orbit mapping**: lays uniform Na×Nb×Nc grid over [0,1)³, uses symmetry ops to group equivalent points, returns irreducible representatives (minimal CASTEP scan set); `write_pes_cube` — unified Gaussian Cube output (line-2 JSON metadata with plane, mobile atom, lattice, ranges); `collect_pes_energies` — auto-detects `irred_NNNNN` vs `grid_III_JJJ_KKK` format, reads `irred_coords.dat` sidecar for index mapping, rewrites cube with energies; `symmetry_expand_energies` — **forward orbit expansion**: reads `irred_coords.dat`, maps each irreducible point to cube energy via nearest-neighbour, applies ALL symmetry ops (min energy for duplicates) → 100% fill of expanded grid; `write_expanded_cube` — N+1 format with boundary wrapping (e.g. 9-point orbit grid → header 10, data 10³, periodic wrap at boundary). Private: `find_castep_in_dir`, `parse_castep_energy`, `element_to_z`, `compute_lattice_vectors`, `wrap_to_unit` (uses `floor(x+eps)` for boundary safety). Depends on `castep_config` + `parser`.

16. **phonon_modes.f90** — `phonon_modes` module. Phonon eigenvector parsing, Born effective charge parsing, and mode decomposition for 3D visualization. Defines `phonon_modes_data_t` (now also stores the first q-point as `qpoint`, `qpoint_weight`, `qpoint_parsed`), `phonon_mode_t`, `born_charge_t` types. Public: `parse_phonon_eigenvectors` — reads .phonon file header (structure + masses + frequencies + eigenvectors), `parse_castep_born_charges` — extracts 3×3 Born effective charge tensors per atom from .castep, `compute_mode_decomposition` — decomposes each phonon mode to per-atom contributions using Z*·u formalism (mode effective charge vector p_m, atom contribution fractions 0..1), `free_phonon_modes_data`. Error codes: `IO_EIGENVECTORS_NOT_FOUND=120`, `IO_BORN_MISMATCH=121`, `IO_BORN_NOT_FOUND=122`. Depends on `castep_config`.

17. **crystal_json.f90** — `crystal_json` module. JSON bridge between Fortran and Rust crystal-viewer. Public: `write_crystal_json` (from `castep_config_t`), `write_crystal_json_cif` (from `cif_data_t` — auto-converts fractional→Cartesian), `write_crystal_json_modes` (from `phonon_modes_data_t` — writes crystal structure + per-mode per-atom displacement vectors), `read_crystal_json_to_cif` (parse modified JSON back into `cif_data_t`). Private: `lattice_vectors` (Cartesian lattice from cell params), `extract_json_real`, `extract_json_string`. Depends on `castep_config` + `parser` + `phonon_modes`.

18. **drift_analysis.f90** — `drift_analysis` module (NOT compiled — development artifact for future anisotropic diffusion analysis). Contains `compute_drift_rates` (per-direction linear drift from unwrapped dipoles) and `compute_global_dielectric` (global covariance ε_ion tensor). Depends on `config` + `polarizability`.

19. **main.f90** — `CASTEP_Suite` program. Suite top-level `do` loop (1. PreCASTEP, 2. PosCASTEP, Q. Quit). PreCASTEP logic in `run_precastep_workflow(should_exit)`: init defaults → file recognition → `run_main_menu` → parse reactant → **if CINEB: parse product + intermediate structures with atom-count validation** → compute lattice → write .cell/.param → exit. Also contains `run_precastep_with_cif(cif, should_exit)` for viewer→PreCASTEP handoff (no file parsing) and `populate_cfg_from_cif` helper. All `stop` replaced with `return`. Coordinate system data-driven from parser. Helpers: `real2str_dp`, `get_file_extension`, `to_lower_inline`.

### Rust crystal-viewer subproject (`crystal-viewer/`)

Standalone Bevy 0.15 3D application (13 source files):

- **`crystal.rs`** — `CrystalData`, `Lattice`, `AtomData` types with serde JSON. `Lattice::to_vectors()` (Cartesian lattice), `Lattice::from_cartesian_vectors()` (norms+angles from arbitrary Cartesian vectors, for slab/vacuum reconstruction), `Lattice::inverse_vectors()` (3×3 inverse), `Lattice::apply_inverse()` (M⁻¹ × vector), `CrystalData::expand_to_cell()` (asymmetric→full unit cell via ±1 fractional translations) + `display_positions()` (two-mode display gate: (1) plain structures tile via `expand_to_cell` (±1 images on all 3 axes); (2) ANY slab/vacuum/supercell structure — a provenance field set — is rendered by `display_boundary_complete()`: the stored atoms plus equivalent face/edge/corner copies on every PERIODIC axis (`periodic_display_axes()`: plain and plain-confirmed-superlaced structures are periodic on all 3 axes; a slab/vacuum cut region is periodic on its two in-plane axes (`display_inplane_axes()` helper: slab → a,b; vacuum-only → the two axes other than the vacuum axis) and on the out-of-plane axis (slab normal c, or the vacuum axis) ONLY if `supercell[axis] > 1`). An atom ON a periodic cell face (frac ≈ 0 or ≈ 1 within 1e-3) is shared with the neighbouring cell, so copies land on BOTH faces of the drawn box (corner atom → all corners, edge atom → both midpoints, interior atom → single in-cell copy), wrapped to [0,1] so no copy lands outside the box — i.e. the edge/face equivalent atoms of the drawn box are always rendered, matching the plain ±1-tiling look (Cu (001) T=2.5 primitive slab: 4 corner copies of the bottom-layer atom + 1 top-layer atom at ½,½ = 5 spheres, all inside the drawn cell; a 2×2×2 plain supercell shows the 16 interior in-plane columns PLUS the 9 face/edge/corner equivalents on the merged-box boundary = 25 top-down blobs)) + `display_box_edges()` (plain structures and confirmed supercells → full 12-edge box; slab-only (no vacuum yet) → only the 4 in-plane bottom-face edges (in-plane rectangle); slab+vacuum or vacuum-only → full 12-edge box of the c=T+V slab+vacuum cell) + `to_json()`/`write_to_file()`. `InPlaneBasis` enum — `primitive` (true 2D Bravais lattice of the layer), `orthogonal` (right-angle in-plane cell → 90/90/90), `conventional` (conventional-cell translations only); serde `rename_all = "snake_case"`, `Default` = Primitive, `as_u8()`/`from_u8()` for the E2E hook. `CrystalData` carries optional provenance fields `slab: Option<SlabMeta>` (h,k,l + start/thickness in Å + `u`/`v` in-plane expansion + `basis` + MS-style explicit in-plane vectors `u_vec`/`v_vec: Option<[i32;3]>` (conventional-basis `(i j k)` rows that override `u`/`v`/`basis` when set)) and `vacuum: Option<VacuumMeta>` (axis + thickness + legacy `both_sides` + `position` 0..=1 vertical placement: fraction of the vacuum below the structure, 0 = structure at the bottom, 1 = at the top) and `supercell: Option<[i32; 3]>` (last confirmed x/y/z supercell multipliers); all `#[serde(default)]` — absent in Fortran-written JSON, and the Fortran line-based JSON reader ignores the unknown keys when the viewer writes them back.
- **`slab.rs`** — slab cross-section + vacuum layer builders (pure math, unit-tested). `build_slab(data, SlabParams{h,k,l,start_ang,thickness_ang,u,v,basis,orth_idx,u_vec,v_vec})`: optional MS-style explicit in-plane vectors `u_vec`/`v_vec: Option<[i32;3]>` (conventional-basis `(i j k)` rows, Materials Studio supercell-matrix format) override the basis × U/V expansion when set — validated in-plane (n̂·ap ≈ 0) and non-collinear via `explicit_inplane_vectors()` (out-of-plane or collinear → `Err` with an example hint); `slab_layer_info()` reports the explicit cell's a′/b′/γ and area-scaled atoms-per-layer. General (hkl) Miller-index termination — surface normal n̂ from reciprocal vectors, slab region [s, s+T) along n̂ (s/T in Å, T=0 → 3·d_hkl), atom images enumerated over the FULL conventional translation group in a z-window (`collect_slab_images`, not the 1-D Bézout w-strip, so fine in-plane cells get complete content), in-plane cell = U×V expansion of the chosen 2D basis (`InPlaneBasis`: `primitive` = true 2D Bravais lattice from the in-plane translation group — half-integer structure-preserving candidates, shortest-pair Gauss reduction, `find_layer_bravais`; `orthogonal` = smallest-area right-angle candidate pair of that lattice, `ortho_candidates`, → 90/90/90 cell; `conventional` = exactly-in-plane conventional combinations, `find_inplane_basis`), coincident wrapped points deduped (with a near-1.0→0.0 fractional snap), c' = T·n̂, atoms re-expressed in the new standard-setting frame via `Lattice::from_cartesian_vectors`. `build_vacuum(data, VacuumParams{axis,thickness_ang,both_sides,position})`: extends the chosen a/b/c axis by V Å; vertical placement by `vacuum_position()` — `position` 0..=1 = fraction of the vacuum below the structure (0 = structure at the bottom of the vacuum, 1 = at the top; 0/0.5/1 presets in the UI), legacy `both_sides` → 0.5, absent → 0. `build_supercell(data, SupercellParams{x,y,z})`: multiplies the lattice vectors by (x, y, z) (multipliers < 1 clamp to 1) and replicates every atom into the x·y·z copies at the integer translations (i, j, k), re-expressed in the enlarged cell; keeps the slab/vacuum provenance and records `supercell: Some([x,y,z])` — a c-superlaced slab+vacuum stack keeps its vacuum gaps (the catalysis use). All three builders return `modified: true` + provenance metadata (incl. `u`/`v`/`basis` on `SlabMeta`, `position` on `VacuumMeta`, the multipliers on `supercell`); callers rebuild entities + auto-save. `slab_d_hkl()` / `slab_period()` (gidx·d_hkl, the layer period the UI's fractional s/T are fractions of) / `slab_layer_snap(data,h,k,l,s,up)` (nearest atomic plane above/below s — the panel's ↑/↓ snap buttons) / `inplane_basis_params()` (chosen U/V/basis in-plane vectors for preview + UI) / `slab_layer_info()` (live a′/b′/γ + atoms-per-layer + orthogonal candidates for the panel).
- **`pes.rs`** — `PesData` type (plane, nx/ny grid, fractional ranges, lattice, structure atoms, energies, has_energies). `generate_surface()` — builds 3D triangle mesh from grid energy data: fractional coords → Cartesian, energy → height along plane normal, jet colormap 256×1 texture (blue→red) with UV mapping by energy value, front+back faces for double-sided rendering. `crystal_data_from_pes()` — converts PesData structure atoms to CrystalData for shared atom/cell rendering pipeline.
- **`render_export.rs`** — Bevy-native offscreen render export (WYSIWYG). `start_offscreen_render` (consumes `RenderSettings.request`) creates a large `Image` (Rgba8UnormSrgb, usage RENDER_ATTACHMENT|TEXTURE_BINDING|COPY_DST) + a temporary camera copying the main camera's Transform/Projection/Tonemapping/Exposure/Msaa with `Camera.target = RenderTarget::Image`, plus `Screenshot::image(handle)`; after the next frame `save_to_disk` (bevy built-in, PNG/TIFF by extension) writes `render.<ext>` and `finish_offscreen` despawns the temp camera/frees the texture/publishes `last_status`. Any vis mode renders (no mesh rebuild — the GPU mesh is reused). `TonemapChoice` (TonyMcMapface/ACES/AgX/Reinhard — no None/linear) + `ImgFormat` (Png/Tiff). Headless E2E hook: `CRYSTAL_VIEWER_AUTORENDER=<WxH>` fires one request 4 s after startup, pair with `CRYSTAL_VIEWER_AUTOEXIT` and check `render.png`.
- **`main.rs`** — App setup, `orbit_camera` (yaw/pitch quaternion via right-drag, ortho/perspective toggle via P), `rotate_camera_keys` (arrow keys ←→↑↓ with configurable angle in right panel), `move_selected_atom` (HKUMIN 6-direction movement with step via DragValue in right panel, default 0.5 Å), `display_mode_system` (1/2/3 ball-stick/space-filling/wireframe, A axes, B bonds, C cell), `spawn_cell_axes` (red X/green Y/blue Z arrows from origin, 1.5× lattice vectors with cone tips), `toggle_projection` (P key), `ortho_projection()` helper. `apply_render_params` (runs on `RenderSettings` change) pushes the dialog's live scene parameters into the running viewer: `KeyLight`/`FillLight` `DirectionalLight.illuminance` (defaults 4000/2000 lx) + key-light `shadows_enabled` (default off), `AmbientLight.brightness` (80), the camera `Msaa` (Off/2x/4x — 8x is clamped to 4x because Metal/WebGPU only guarantee 1-4 samples for Rgba8UnormSrgb and requesting 8 aborts the viewer), the camera clear color (RGB sliders + Black/White/Light-gray 170/Dark-gray/Viewer-default 43,44,47 presets), atom `perceptual_roughness` (0.5) + `metallic` (0.2, via `PickingState.atom_material_handles`), and the camera `Tonemapping`. SSAO was removed on purpose: Bevy 0.15's SSAO strength is hardcoded in the shader (no intensity parameter) and was visually negligible in these scenes. NOTE: the two light queries MUST stay merged in one `Query<(&mut DirectionalLight, Option<&KeyLight>, Option<&FillLight>)>` — separate `With<KeyLight>`/`With<FillLight>` queries conflict (B0001). **PES 2D mode**: auto-detects `.cube` file with `nz==1` + `"type":"pes_2d"` metadata → `PesData::from_cube()` → surface mesh with colormap (S toggles visibility, +/- adjusts color clip). **PES 3D mode**: detects `.cube` with `nz>1` + `"type":"pes_3d"` → MC isosurface (4), volume render (5), slice planes (6), fixed-radius sphere (7), radial-stationary migration surface (8). Default for PES 2D: **orthographic projection**, bonds hidden, axes hidden, camera perpendicular to XY plane. **Slab/vacuum systems** (see `slab.rs`): `SlabState` resource (h/k/l text fields, dual s/T inputs — fractional-of-layer-period and absolute Å, last-edited-wins sync, ↑/↓ atomic-plane snap buttons — plus in-plane basis choice `InPlaneBasis` + U/V expansion + orthogonal index + MS-style explicit in-plane vector fields `u_vec_str`/`v_vec_str` ("0 1 0" / "0,1,0" / "(0 1 0)" parsed by `parse_iv_vec`, empty = off, override basis × U/V) + vacuum axis/thickness (coupled with the total cell c) + placement preset/position 0–1 + popup visibility flags `slab_open`/`vacuum_open`/`supercell_open` (the right panel's three buttons open the settings windows, × closes), supercell multipliers `sc_x`/`sc_y`/`sc_z` (default 1,1,1 = the current cell), request + preview + snapshot flags; `last_hkl` re-derives the fractional displays when the termination changes), `apply_slab_system` / `apply_vacuum_system` / `apply_supercell_system` (MAX_SLAB_ATOMS=50000 guard; snapshot original for Reset; rebuild all structure entities via `rebuild_structure_entities()`; **`refit_camera_to_cell()` re-centres + re-scales the camera to the new cell after slab/vacuum/supercell/reset so the taller/narrower slab stays in view** — also refreshes `CameraInit` so the R-key reset follows the new framing; auto-save modified JSON), `rebuild_structure_system` (re-spawn atoms/bonds/cell/axes, refill PickingState/ImageOffsets/AtomInfo/LatticeData), `slab_preview_system` (ghost cut plane / slab region / vacuum boxes aligned to the U/V-aware in-plane basis via `Quat::from_mat3`, despawn+respawn on change), `sphere_section.rs`) — PES 3D pipeline, `apply_reset_system`. E2E headless hooks: `CRYSTAL_VIEWER_AUTO_SLAB=h,k,l,s,T[,u,v,basis,orth]` (basis 0=primitive/1=orthogonal/2=conventional), `CRYSTAL_VIEWER_AUTO_SLAB_UV=ui,uj,uk,vi,vj,vk` (MS-style explicit in-plane vectors, paired with AUTO_SLAB), and `CRYSTAL_VIEWER_AUTO_VACUUM=axis,V,both[,pos]`, and `CRYSTAL_VIEWER_AUTO_POPUP=slab|vacuum|supercell|off` (open the settings popup window at the same 4 s mark, for offscreen/UI verification) fire once 4 s after startup via `auto_slab_system`; `CRYSTAL_VIEWER_AUTO_SUPERCELL=x,y,z[,preview]` sets the multipliers and applies the supercell (4th token `preview` = stay in the live preview state, structure untouched), pair with `CRYSTAL_VIEWER_AUTOEXIT` and assert on the auto-saved JSON; `CRYSTAL_VIEWER_AUTO_RENDER_DELAY=<s>` overrides the 4 s offscreen-capture delay of `CRYSTAL_VIEWER_AUTORENDER` (use > 4 s to capture AFTER an auto slab/vacuum apply); `CRYSTAL_VIEWER_AUTORENDER_UI=[file]` (`auto_ui_screenshot_system`) writes a `Screenshot::primary_window()` framebuffer capture of the running viewer 5 s after startup (3D scene included; the egui overlay is NOT in that capture — use a macOS `screencapture` of the real window for UI-layout checks); `CRYSTAL_VIEWER_DEBUG_UI=1` logs each frame the slab/vacuum popup branch executes (popup open-state + egui screen rect).
- **`picking.rs`** — MVP-projection-based atom picking. `PickingState` with parent indices for symmetry-aware selection. `click_pick`/`hover_pick`/`highlight_atoms` (highlights all symmetry-equivalent atoms). Pick radius = 6% of viewport height. Custom `over_egui_panel` check prevents picks/hover over UI panels.
- **`ui.rs`** — egui panels: left 180px (asymmetric-unit atom list + periodic table popup), right 220px (cell params + three settings buttons, **Slab cut / Vacuum / Supercell**, that open the slab/vacuum/supercell settings in floating popup windows + the compact applied-state provenance lines + selected atom + Rotation Angle / Move Step DragValue + the "Render" button), bottom toolbar (concise key hints). The slab/vacuum/supercell settings live in three popup windows — `egui::Area` + `Frame::window`, same pattern as the Render dialog (draggable from anywhere, position remembered, × button closes; open flags `SlabState.slab_open`/`vacuum_open`/`supercell_open`, content built by `slab_section_ui`/`vacuum_section_ui`/`supercell_section_ui` at the end of ui.rs) — keeping the right panel compact. The slab popup hosts: typed h/k/l + slab position s (dual inputs: fraction of the layer period gidx·d_hkl + absolute Å, last-edited-wins sync) with "up"/"down" snap-to-nearest-atomic-plane buttons (disabled when no plane lies in the scan window) + thickness T (dual frac/Å, 0 = auto 3·period, live d_hkl/period display) + a **U/V definition mode selector** ("U/V def" ComboBox: "integer × basis" [default] or "vector (i j k)") — in **integer mode** the in-plane basis ComboBox (primitive = true 2D Bravais lattice / 90° = smallest right-angle in-plane cell with candidate index / conventional = conventional-cell translations) + U/V expansion DragValues (1..=8, in-plane supercell) are shown and the explicit vec text fields are hidden; in **vector mode** the MS-style explicit in-plane vector text fields "U vec (i j k)" / "V vec (i j k)" (conventional-basis rows, Materials Studio supercell-matrix format) are shown and the basis combo + DragValues are hidden (empty vec fields fall back to basis × U/V internally) + live `slab_layer_info` (in-plane a′/b′/γ, atoms per layer, orthogonal candidates — reports explicit-cell dimensions in vector mode with active vectors, basis × U/V expansion otherwise), Preview/Apply slab buttons + the error label. The `SlabState.uv_mode` field (0 = integer, 1 = vector) gates which fields are displayed; `slab_params_from_state` enforces the gate: integer mode forces `u_vec`/`v_vec` to `None` even if stale vec text remains, vector mode passes the parsed vectors through. The vacuum popup hosts: thickness V DragValue (default 15 Å) + total cell c TEXT field coupled to it (empty field shows the hint "Auto" = c follows V as base + V; type a value to pin c, then V = c − base; base = axis length of the current slab cell; type "Auto"/clear to go back to auto; last-edited-wins; the apply itself only uses V) + direction a/b/c + placement bottom/center/top presets (center = legacy both-sides split) + position 0–1 slider (fraction of the vacuum below the structure) + Preview/Apply + Reset. The supercell popup hosts: x/y/z integer multipliers (DragValues 1..=8; 1,1,1 = the current cell, 2,1,1 = one extra cell along a) + a live cells/atoms count + the applied-multipliers line when `supercell` provenance is set + **Preview** button (toggles the live ghost preview: every original cell keeps its own box edges — boxes side by side — plus all atom copies, so the tiling is visible; the preview is capped at 512 cells) + **Apply supercell** button (confirms: lattice ×(x,y,z), atoms replicated, the internal box edges merge into one box, auto-save with `supercell` provenance). Applying rebuilds all structure entities and auto-saves the modified JSON (slab/vacuum provenance keys, incl. U/V/basis and vacuum position), so the Fortran modified-structure menu picks up the result. The Render dialog (built on `egui::Area` with `default_pos` — draggable from anywhere and the position is remembered; never use `anchor`, it re-applies every frame and snaps the window back to center): resolution W×H, format (PNG/TIFF), MSAA dropdown (Off/2x/4x), and the live scene parameters (key/fill/ambient lux sliders, Shadows checkbox, atom roughness + metallic sliders, tonemapping dropdown, background RGB sliders + 5 presets) — all edits apply to the viewer immediately (WYSIWYG), then the Render button triggers the offscreen export; bottom row: Render | Reset all | Close (natural button widths). The right panel hosts "Render" + "Export PLY" side by side (PLY only shown in mode 8 with a cube, E key handled by `ply_export_key_system` so `ui_system` stays under the 16-system-parameter limit of Bevy 0.15 fn systems). Panels fixed-width, non-resizable.
- **`resources.rs`** — Periodic table data: covalent radii and CPK/Jmol colors for element→color mapping.
- **`cube_reader.rs` / `marching_cubes.rs` / `volume_render.rs` / `slice_plane.rs` / `sphere_section.rs`** — PES 3D pipeline: cube parsing, MC isosurface, volume rendering, slice planes, radial migration surface (mode 8, with cage-center detection and shell welding). The right-panel **"Export PLY"** button (next to Render, mode 8 with a loaded cube, or E key) writes `migration_surface.ply` — the shell mesh with baked jet vertex colors — for external rendering (e.g. Blender: import PLY; the vertex colors become a vertex-color layer). Publication images come from the Bevy-native render dialog (offscreen WYSIWYG) or a viewer screenshot. The sphere_section tests read a local synthetic cube from `/tmp/sphere_test.cube` (a real PES-3D export snapshot kept outside the repo; tests fail if it is absent).

Key features:
- **Arrow key rotation**: ←→ yaw, ↑↓ pitch, angle configurable via DragValue (1°–90°, default 45°)
- **egui context access rule**: ALWAYS `EguiContexts::try_ctx_mut()` (never `ctx_mut()`) — on the closing frame Bevy despawns the window before our Update systems run, bevy_egui then skips context init and `ctx_mut()` panics (`uninitialized context`) on a Compute Task Pool thread; with `panic = "abort"` (release profile) that is the SIGABRT crash macOS reported on every viewer close. Pattern: `let Some(ctx) = contexts.try_ctx_mut() else { return; };` for UI entry points, `contexts.try_ctx_mut().is_some_and(|c| c.wants_pointer_input())` for guards. Reproduce/verify with `CRYSTAL_VIEWER_AUTOEXIT=<seconds>` env (sends AppExit after N s; same teardown as closing the window) — exit code must be 0.
- **A-key axes**: toggle red/green/blue axes arrows from cell origin along lattice vectors (1.5× length)
- **Cell expansion**: asymmetric unit → full unit cell images via ±1 fractional translations
- **Symmetry-equivalent sync**: moving one atom moves all its equivalent images
- **Zoom coupling**: perspective radius ↔ orthographic scale (1.09× factor), no jump on P toggle
- **Non-orthogonal cell support**: correct M⁻¹ × vector (row-based, not column-dot) for triclinic/hexagonal cells
- **Gimbal lock prevention**: `looking_at` instead of `look_at` for camera light
- **Keyboard isolation**: egui field focus blocks Bevy shortcuts (no accidental triggers while typing)
- **Auto-save**: modified positions written to JSON on each move; Fortran reads back after viewer exits
- **Slab cross-section + vacuum layer** (catalysis workflow, viewer-only, zero Fortran changes): Cell Parameters panel → general (hkl) termination with position s / thickness T as DUAL inputs (fraction of the layer period gidx·d_hkl + absolute Å, last-edited-wins sync; "up"/"down" snap s to the nearest atomic plane above/below; T=0 → 3·d_hkl). The in-plane cell is the U×V expansion of a user-chosen 2D basis — primitive (the layer's true 2D Bravais lattice, structure-preserving translation search), 90° (smallest right-angle in-plane cell, giving a 90/90/90 slab cell), or conventional (conventional-cell translations); the default is the primitive lattice. **MS-style explicit in-plane vectors**: optional `(i j k)` triple input (Materials Studio supercell-matrix rows in the conventional basis, e.g. U=`(0 1 0)`, V=`(0 0 1)` for a (100) slab) overrides basis × U/V — vectors must be in-plane and non-collinear, otherwise the apply shows an error with an example; the single-number U/V DragValues remain as the uniform-expansion shortcut. The settings live in two popup windows opened from the right panel's **Slab cut / Vacuum** buttons (`egui::Area` + window Frame, draggable, position remembered, × closes; both windows hug their content — no forced min width, and the title row must NOT use a `right_to_left` sub-layout, which claims the whole available width on the sizing pass and stretches the frame to the screen edge) — the panel itself keeps only the buttons, a compact applied-state line, and the error label. After a slab cut (no vacuum yet) the structure is **not** a periodic unit cell: the viewer renders the cut in-plane cell with boundary sharing — each slab atom at its in-cell position, plus equivalent copies at the cell corners/edges only when the atom sits on the in-plane cell boundary (a corner atom appears at all 4 corners of the in-plane cell as one shared Cu; an interior atom, e.g. the Cu (001) top-layer atom at ½,½, is drawn ONCE — Cu (001) T=2.5 primitive → 4+1 = 5 spheres, all inside the drawn cell) — with NO copies along the out-of-plane/vacuum axis. Slab-only (no vacuum): draws the 4 in-plane (bottom-face) boundary edges — a flat in-plane rectangle. Once a vacuum layer is built, the full 12-edge box of the slab+vacuum cell (c=T+V) is drawn — the complete 3D frame the user expects. Plain (non-slab) structures keep the full-cell ±1 tiling and the 12-edge box. Applying slab/vacuum/reset re-fits the camera (`refit_camera_to_cell`) to the displayed region (cell plus the 2×2 in-plane block for slab/vacuum) so the new cell stays centred and in view. Vacuum builder: axis a/b/c, thickness V default 15 Å, coupled with the total cell c (c = base + V), vertical placement bottom/center/top or a 0–1 position slider (0 = structure at the bottom of the vacuum, 1 = at the top; center = legacy both-sides split). Applying rebuilds the structure and auto-saves the modified JSON; the existing modified-structure menu (save .cell / hand off to PreCASTEP) picks up the result. `slab`/`vacuum` JSON keys are provenance only — the Fortran JSON reader ignores them. E2E: `bash test/slab_e2e.sh` (auto slab (001, 6 layers × 1 atom, L2 primitive) + 15 Å vacuum on `test/Cu.cif` → c = 25.8 Å through both .cell paths, plus (111) primitive/90°/conventional cases → 3/6/12 atoms, c ≈ 21.26 Å).
- **Supercell expansion** (viewer-only, zero Fortran changes): the right panel's **Supercell** button opens a popup with x/y/z integer multipliers (1..=8; 1,1,1 = the current cell, 2,1,1 = one extra cell along a). **Preview** shows every ORIGINAL cell keeping its own box edges (boxes side by side — for a slab/vacuum structure each copy shows the 4 in-plane edges) plus the full post-apply display set: tiled atom copies AND face/edge/corner equivalents on the merged-box boundary (computed via `slab::build_supercell` + `display_positions()`, so the preview is a true WYSIWYG of what Apply will produce); the structure is NOT modified until **Apply** (grid capped at 512 cells). After the apply the internal edges merge into one box, the lattice vectors scale by (x,y,z), the atoms replicate, the camera re-fits, and the JSON auto-saves with `supercell` provenance. A c-superlaced slab+vacuum stack keeps its vacuum gaps (catalysis workflow; the c copies on the top face are also rendered because the c axis becomes periodic when `supercell[c] > 1`); an in-plane supercell on a slab disables the 2×2 in-plane block display (`display_inplane_axes()` → None, the cell already spans the pattern) and the enlarged in-plane cell is then boundary-completed on its own faces/edges like any periodic cell. After any supercell apply the merged box renders its face/edge/corner equivalent atoms (the same ±1 completion as the plain tiling, restricted to the periodic axes — `display_boundary_complete` in `crystal.rs`). E2E: `test/ui_screenshot_check.sh` cases `sc_preview` (AUTO_SUPERCELL=2,2,2,preview → the saved JSON stays identical to the input) and `sc_apply` (2,2,2 → a=b=c = 2× the original, 32 atoms, `supercell` [2,2,2], 25 top-down atom blobs = 16 interior 4×4 in-plane columns + 9 face/edge/corner equivalents on the merged-box boundary). Headless: `CRYSTAL_VIEWER_AUTO_SUPERCELL=x,y,z[,preview]` (+ `AUTO_POPUP=supercell` for the popup window).
- **Profile**: LTO, strip, opt-level="z" → 11 MB binary; auto-cleans intermediate build files

### Module dependency chain

```
castep_config (leaf)
term_utils    (leaf)
  ├── symmetry        (config)
  ├── parser          (config)
  ├── cell_writer     (config)
  ├── param_writer    (config)
  ├── bands_parser    (config)
  ├── pdos_parser     (config)
  ├── phonon_dos      (config)
  ├── phonon_modes    (config)
  ├── dos_compute     (config)
  ├── polarizability  (config)
  ├── thermodynamics  (config)
  ├── castep_vib      (config)
  ├── pes             (config)
  ├── crystal_json    (config + parser + phonon_modes)
  ├── bands_plotter   (config + term_utils)
  ├── dos_plotter     (config + term_utils)
  ├── cli_menu        (config)          ─┐
  └── poscastep_menu  (config +          ├─ siblings
                        term_utils +      │  (no cross-dep)
                        bands_parser +    │
                        bands_plotter +   │
                        pdos_parser +     │
                        phonon_dos +      │
                        phonon_modes +    │
                        dos_compute +     │
                        dos_plotter +     │
                        polarizability +  │
                        crystal_json +   ─┘
                        parser)
  └── main            (config + parser + cell_writer + param_writer
                        + cli_menu + poscastep_menu)
  drift_analysis      (config + polarizability)  ← NOT compiled

crystal-viewer/       (Rust — standalone 3D viewer subproject)
```

`cli_menu` and `poscastep_menu` are sibling modules — they share types/utilities from `castep_config` but do not depend on each other. Post-processing code lives in `poscastep_menu` + `pes` + `bands_parser` + `bands_plotter` + `pdos_parser` + `dos_compute` + `dos_plotter` + `polarizability` + `phonon_modes` + `term_utils`. `drift_analysis` depends on `config` + `polarizability` but is not compiled into the main program — it is a standalone development artifact for future anisotropic diffusion analysis.

## Key Design Notes

- **No external dependencies** — pure Fortran 2018, no libraries beyond stdlib.
- **Suite dual-mode architecture** — `main.f90` top-level loop dispatches to `run_precastep_workflow` (PreCASTEP input generation) or `run_poscastep_menu` (PosCASTEP post-processing). The PreCASTEP mode exits the program after successful generation; PosCASTEP returns to suite.
- **No `stop` statements** — all error/recovery paths use `return`. `IO_USER_QUIT` (-1) signals user-requested exit from sub-menus.
- **Q label convention** — sub-menus show "Q. Back" (PreCASTEP config, PosCASTEP); only the suite top-level menu shows "Q. Quit".
- **PreCASTEP exit after generation** — `run_precastep_workflow` sets `should_exit=.true.` on successful completion, causing the suite loop to exit. Errors/quits return to suite via `return`.
- **`strip_quotes` in config.f90** — single canonical implementation shared by cli_menu and poscastep_menu.
- **Bands gap analysis** — `bands_plotter` detects VBM (highest eigenvalue below E_F) and CBM (lowest above E_F), computes indirect gap (VBM→CBM regardless of k-point) and direct gap at VBM k-point, classifies as direct/indirect/metallic (gap < 0.005 eV).
- **VBM/CBM markers** — green ◆ for VBM, yellow ◈ for CBM in ASCII plot. Rendered via `char_type` types 6 and 7, `type_color` cases 6→C_GREEN and 7→C_YELLOW.
- **eV-only display** — Hartree energy units removed entirely from CLI output and plot code. All energies displayed in eV using `HARTREE_TO_EV = 27.211386245988_dp`.
- **Energy range** — default `e_center ± half_range` (10 eV), scrollable with ↑↓ keys, zoomable with +/-. `e_center` defaults to Fermi energy.
- **K-path windowing** — optional `k_pct_in`/`k_width_pct_in` parameters (percentages 0–1) define a k-path display window. Defaults to full path (center=0.5, width=1.0). Left/right arrow keys auto-zoom to 50% window on first press, then pan in 10% increments. `plot_bands_ascii` signature: `(bands, term_w_in, term_h_in, e_center, half_range, k_pct_in, k_width_pct_in)` where the last two are optional.
- **Gap info above plot** — band gap information displayed in 2 compact lines above the plot box. Line 1: E_F + Band Gap + type. Line 2: VBM + CBM + Dir gap (if indirect). Grid height dynamically reduced by 9 rows to keep all content visible without scrolling.
- **Resize-adaptive** — terminal size re-detected via `stty size` on every redraw. Resize the terminal and press any key to adapt.
- **Terminal mode** — `enter_raw_mode`/`leave_raw_mode` in `term_utils` encapsulate `stty -icanon -echo min 1` + alternate screen buffer (`ESC[?1049h`/`ESC[?1049l`). The alternate screen buffer prevents interactive plot redraws from polluting terminal scrollback history. `stty -icanon` (not `stty raw`) preserves `opost` for correct newline translation, avoiding staircase alignment bug.
- **Compact layout** — label-plot gap = 1 column, no "k-path" text on x-axis, maximizing grid space.
- **Consistent 'q' handling** — typing `q` at any file path prompt in PosCASTEP returns to the PosCASTEP menu (not the suite menu). Uses `IO_USER_QUIT` internally, intercepted by each handler.
- **Phonon DOS/IR/Raman** — `phonon_dos` module parses `.phonon` files, computes phonon DOS (Gaussian smearing), IR absorption spectrum, and Raman scattering spectrum. All three use the same `plot_dos_ascii` renderer with configurable `xlabel`/`xunit`.
- **Multi-format parsing in parser.f90** — CIF, PDB, and .cell parsing extracted to `parser` module. File type determined by extension.
- **PDB parsing** — CRYST1 record (columns), UNITCELL record, and ATOM/HETATM records. All PDB coordinates are Cartesian.
- **LATTICE_CART support** — `parse_cell_inline` handles both `%BLOCK LATTICE_ABC` and `%BLOCK LATTICE_CART`; `compute_abc_from_cartesian` converts back. CASTEP .cell files are column-major.
- **Output always LATTICE_ABC and POSITIONS_ABS** — regardless of input format. Fractional→Cartesian conversion uses `cell_basis` matrix; LATTICE_CART→ABC uses `compute_abc_from_cartesian`.
- **Shared types in castep_config** — `cif_data_t`, `atom_t`, `bands_data_t` all defined in config module.
- **Allocatable working arrays** — avoids stack overflow for large atom counts.
- **Character truncation safety** — `copy_str_no_quotes` strips quotes by character copy.
- **Element symbol cleaning** — `clean_element_symbol` strips oxidation state suffixes.
- **Space group priority** — H-M name tags take priority over IT number.
- **Conditional CELL_CONSTRAINTS** — only for GEOMETRYOPTIMISATION with cell_opt_mode ALL.
- **Conditional geo params in menu** — only displayed for geometry/structure tasks.
- **Symmetry handling simplified** — NONE/AUTO only.
- **Geo tolerance thresholds** — COARSE (5e-5/0.1/0.2/0.005), MEDIUM (2e-5/0.05/0.1/0.002), FINE (1e-5/0.03/0.05/0.001), EXTREME (5e-6/0.01/0.02/5e-4).
- **SOC19 ↔ spin_polarized linkage** — selecting SOC19 auto-sets `spin_polarized = .true.`; `spin_orbit_coupling : true` written when both are active.
- **Smearing logic** — off → `fix_occupancy : true`; on → `fix_occupancy : false` + `nextra_bands` + `smearing_width`.
- **vdW sedc** — `sedc_apply : true` + `sedc_scheme : <method>` when vdW ≠ NONE.
- **Cutoff energy formatting** — integer when whole-number value, otherwise F12.4.
- **Crystal viewer integration** — PosCASTEP -1 parses CIF/PDB/cell → `write_crystal_json_cif` generates JSON (Cartesian coords) → `launch_viewer` spawns `crystal-viewer` process → `wait=.true.` blocks until viewer closes → `read_crystal_json_to_cif` reads modified JSON → if modified, shows sub-menu → JSON auto-deleted after handling.
- **Viewer→PreCASTEP handoff** — Option 2 in modified-structure menu copies `cif_data_t` to module-level `precastep_cif_data` in `poscastep_menu`, returns `IO_PRECASTEP_LAUNCH` to main, which calls `run_precastep_with_cif` skipping file parsing. `run_main_menu` skips `ask_input_file` when `cfg%num_atoms > 0`.
- **Viewer auto-detection** — `find_viewer()` in poscastep_menu searches relative to executable: dev layout (`crystal-viewer/target/release/crystal-viewer`) then release layout (`crystal-viewer` in same dir).
- **Viewer defaults**: orthographic projection, camera perpendicular to XY (yaw=0, pitch=0), bonds hidden, zoom auto-fit to cell.
- **Fractional↔Cartesian conversion** — Rust uses row-based M⁻¹×v (not column-dot) for correct non-orthogonal cells. Both Fortran and Rust compute lattice vectors identically.

## Extending PosCASTEP

### Adding a new analysis type to PosCASTEP:

```fortran
! In poscastep_menu.f90, add constant:
integer, parameter :: POS_DOS = 2

! In run_poscastep_menu, add menu line and case:
write(*, '(a)') '  2. Plot Density of States'
case (POS_DOS)
    call handle_dos_menu(iostat)

! Add private subroutine:
subroutine handle_dos_menu(iostat)
    ...
end subroutine
```

## Extending the Writer Modules

### Adding new %BLOCK to .cell:
```fortran
subroutine write_block_newtype(unit, cfg)
    ...
end subroutine

! Then in write_cell_file:
call write_block_newtype(unit, cfg)
```

### Adding new keyword to .param:
```fortran
call write_kv(unit, 'new_key', 'new_value')
```
Place it in the appropriate section of `write_param_file` (common or task-diff).

## Testing

```bash
bash test/run_tests.sh <test_name>
```

Quick tests (200 eV, -np 1): `singlepoint`, `geomopt`, `phonon_fd`, `phonon_dfpt`, `efield`, `phonon_efield`, `all`.

Phonon suite (800 eV, NCP19, smearing ON, -np 44, GeomOpt→4 phonon variants): `phonon_suite`.

Test structure: `test/Cu.cif` (4-atom FCC Cu). Output in `test/output/<test_name>/`. Each test builds PreCASTEP, generates input via printf, runs CASTEP, checks `.castep` for errors, reports PASS/FAIL.

Viewer slab/vacuum E2E (no CASTEP needed): `bash test/slab_e2e.sh` — builds release, drives PosCASTEP -1 on `test/Cu.cif` with the viewer's headless hooks (`CRYSTAL_VIEWER_AUTO_SLAB=h,k,l,s,T[,u,v,basis,orth]`, `CRYSTAL_VIEWER_AUTO_VACUUM=axis,V,both[,pos]`, `CRYSTAL_VIEWER_AUTOEXIT=12`), then verifies both modified-structure menu paths: option 1 (.cell save, (001) case expects c = 25.8 Å + 6 atoms — 6 (001) layers × 1 atom of the L2 primitive square in-plane cell) and option 2 (PreCASTEP handoff → `Cu_Energy.cell`/`.param`), plus (111) cases with the new in-plane basis hook: primitive → 3 atoms, 90° orthogonal → 6 atoms, conventional → 12 atoms (all c ≈ 21.26 Å = 3·d₁₁₁ + 15 Å vacuum). Viewer unit tests: `cd crystal-viewer && cargo test` (43 pass — incl. the 27 slab tests with the explicit MS-style `u_vec`/`v_vec` cases: (100) U=(0 2 0), V=(0 0 1) → 2b×c cell, 24 atoms, and out-of-plane / collinear rejection, the 4 cut-region display tests (`slab_structure_inplane_display_rule` (corner atom → 4 in-cell corner copies at [0,1]², interior atom → 1 in-cell copy), `slab_structure_shows_inplane_rectangle_not_box`, `vacuum_only_axes_follow_vacuum_axis`, `plain_structure_still_tiles_the_cell`: slab/vacuum structures render the cut in-plane cell with boundary sharing and draw the 4-edge in-plane rectangle, vacuum-only structures tile on the two axes other than the vacuum axis, plain structures keep the ±1 tiling + 12-edge box), the 5 supercell display-gate tests (`confirmed_plain_supercell_renders_face_edge_equivalents`: confirmed 2×2×2 / 1×1×1 cells render stored atoms + face/edge/corner copies, single merged 12-edge box; `inplane_supercell_still_completes_inplane_faces`: slab + [2,1,1] → in-plane boundary completion with no z copies; slab + [1,1,2] → c-axis copies enabled; `slab_vacuum_supercell_adds_c_face_copies`: slab+vacuum + c-superlaced → all 3 axes periodic; `plain_supercell_display_completes_merged_box_boundary`: build_supercell 2×2×2 on plain fcc → 35 display positions incl. far corner (2,2,2) copy; `slab_c_supercell_display_adds_top_face_copies`: slab [1,1,2] → 14 display positions, top-face z copies present, 4-edge in-plane box), and the 5 `build_supercell` tests (identity no-op, 2×2×2 on fcc → 32 atoms + doubled lattice + corner/face-copy spot checks, [2,1,1] doubles only a, slab provenance kept by a c-supercell, clamping < 1 → 1); 3 pre-existing `sphere_section` failures are expected without `/tmp/sphere_test.cube`). **Display/UI screenshot-verification hook**: `bash test/ui_screenshot_check.sh` — runs the release viewer headlessly (auto slab / slab+vacuum cases + offscreen WYSIWYG render + `CRYSTAL_VIEWER_DEBUG_UI` log) and asserts the visual state: atom-blob count in the 3D render (thresholded connected components — the offscreen render is the 3D scene only, egui NOT included) and the slab/vacuum/supercell popup window width from the `[debug-ui] … popup rect=` log lines (popup layout cannot come from any render capture). The two supercell cases assert on the auto-saved JSON instead: `sc_preview` (CRYSTAL_VIEWER_AUTO_SUPERCELL=2,2,2,preview → JSON byte-identical to the input fixture) and `sc_apply` (2,2,2 → lattice doubled, 32 atoms, `supercell` [2,2,2], 25 top-down atom blobs — 16 interior 4×4 in-plane columns plus the 9 face/edge/corner equivalents rendered on the merged-box boundary by `display_boundary_complete`). Standing workflow: after every completed display/UI modification, run this script AND launch a one-shot subagent to independently screenshot-verify the modified part before it is considered done. Popup UI can also be exercised headlessly with `CRYSTAL_VIEWER_AUTO_POPUP=slab|vacuum|supercell` (the auto-saved JSON / auto-render still reflect the applied state).

The `precastep-tester` agent wraps the script — invoke it after modifying any source file.

## Default Parameter Summary

| Parameter | Options | Default |
|-----------|---------|---------|
| Task | 17 types | SINGLEPOINT |
| XC functional | PBE, PBEsol, HSE06, PBE0, r2scan | PBE |
| Cutoff (eV) | Any positive number | 400.0 |
| vdW | NONE, D3, D3-BJ, D4 | NONE |
| Pseudopotential | NCP19, C19MK2, SOC19 | C19MK2 |
| K-point | GAMMA, MONKHORST_PACK | GAMMA |
| SCF tolerance | Any positive number (scientific notation supported) | 1e-5 |
| Optimizer | BFGS, LBFGS, CG | BFGS |
| Cell opt mode | ALL, FIX_CELL | FIX_CELL |
| Symmetry | NONE, AUTO | NONE |
| Geo tolerance | COARSE, MEDIUM, FINE, EXTREME | MEDIUM |
| Advanced — Smearing | off, on | off |
| Advanced — Max SCF | Any positive integer | 256 |
| Advanced — Conv Window | Any integer >= 2 | 3 |
| Advanced — ELF | on, off | off |
| Advanced — EDD | on, off | off |
| Spin polarized | true, false (toggled via -1) | false |
| CINEB max images | >= 3, odd | 11 |
| CINEB spring constant | Any positive number | 0.1 eV/Å² |
| CINEB tangent mode | NONE, BISECT, HIGH_E, SPLINE | SPLINE |
| CINEB NEB method | TPSD, FIRE, ODE12R | ODE12R |
| CINEB max iterations | Any positive integer | 50 |
| CINEB climbing image | TRUE/FALSE (always on, not user-configurable) | TRUE |
| TS tolerance | COARSE, MEDIUM, FINE, EXTREME | MEDIUM |
