# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CASTEP Suite is a Fortran 2008 CLI suite with two modes:

1. **PreCASTEP** — converts crystallographic structure files (CIF, PDB, .cell) into CASTEP DFT input files (`.cell` and `.param`)
2. **PosCASTEP** — post-processes CASTEP output: band structure with gap analysis, total DOS, and projected DOS (s/p/d/f)

Launching `./CASTEP_Suite` shows a top-level suite menu. The PreCASTEP mode exits after successful file generation. PosCASTEP returns to its sub-menu (Q. Back) and then to the suite menu (Q. Quit). No `stop` statements remain — all error/recovery paths use `return`.

## Build

```bash
make          # Release build (gfortran, -O2)
make debug    # Debug build (-O0, -fcheck=all)
make run      # Build and run
make clean    # Remove obj/ and CASTEP_Suite binary
```

Requirements: gfortran >= 7.0, Make.

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
- `9/10/11` — Optimizer, Cell opt mode, Geo tolerance (only shown when task is GEOMETRYOPTIMISATION or ElectronicSpectroscopy)
- `20` — Spectral Task (only shown when task is ElectronicSpectroscopy; sub-menu: BandStructure / BandStructure_pDOS)
- `Q. Back` — Return to suite menu (returns `IO_USER_QUIT`)

### PosCASTEP (option 2)

Post-processing menu for CASTEP output analysis:
- `1. Plot Band Structure` — prompts for `.bands` file path, then ASCII terminal plot or SVG output
- `2. Plot DOS` — total density of states from `.bands` only; ASCII (interactive), SVG, or CSV export
- `3. Plot pDOS` — projected DOS: enter file prefix (e.g. `Cu`), auto-loads `<prefix>.bands` + `<prefix>.pdos_bin` (or `.pdos_weights`); ASCII (interactive) or CSV export with s/p/d/f columns
- `Q. Back` — Return to suite menu

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
- SVG output available for total DOS only
- Energy grid: Fermi ±20 eV precomputed at 4001 points, interactive viewport shows subset

## Architecture

Thirteen source files in `src/`:

1. **config.f90** — `castep_config` module. Defines kinds (`dp`), physical constants (`HARTREE_TO_EV`), constants (16 task types, 5 XC functionals: PBE/PBEsol/HSE06/PBE0/r2scan, 4 vdW corrections, 3 pseudopotentials, 3 K-point schemes, 3 optimizers, 4 geo tolerance levels, 3 cell opt modes, 2 symmetry modes, phonon constant groups (method, fine, DFPT, sum_rule, cutoff, q-point), I/O error codes), max sizes (MAX_ATOMS=10000, MAX_TAGS=5000, MAX_LOOP_ROWS=50000, MAX_SYM_OPS=400), CIF/CASTEP tag names, and types (`atom_t`, `cif_data_t`, `castep_config_t`, `bands_data_t`). Allocatable fields: `cif_data_t%atoms`, `bands_data_t%kpoint_indices`, `bands_data_t%kpoint_coords`, `bands_data_t%kpath_dist`, `bands_data_t%eigenvalues`. Advanced option fields in `castep_config_t`: `smearing`, `max_scf_cycles`, `elec_convergence_win`, `calculate_elf`, `calculate_edd`. `castep_config_t` has a `final :: finalize_castep_config` procedure that auto-deallocates `atom_type`, `atom_x`, `atom_y`, `atom_z` on scope exit. Types include `atom_t`, `cif_data_t`, `castep_config_t`, `bands_data_t`, `pdos_data_t` (for `.pdos_bin`/`.pdos_weights` binary data). Error codes: `IO_FILE_NOT_FOUND=100`, `IO_PARSE_ERROR=101`, `IO_MISSING_CELL=102`, `IO_MISSING_ATOMS=103`, `IO_BAD_NUMERIC=104`, `IO_WRITE_ERROR=105`, `IO_INVALID_INPUT=106`, `IO_WRITE_FAIL=107`, `IO_BANDS_NOT_FOUND=108`, `IO_BANDS_PARSE_ERROR=109`, `IO_PDOS_NOT_FOUND=110`, `IO_PDOS_PARSE_ERROR=111`, `IO_USER_QUIT=-1`. Provides `default_config`, `new_castep_config`, tag normalization/case-insensitive comparison (`normalize_tag`, `compare_tags`), `string_to_real`, `int2str`, `get_castep_task_name` (maps user-facing names to CASTEP internal names), and `strip_quotes` (removes leading/trailing quote characters — shared by cli_menu and poscastep_menu).

2. **term_utils.f90** — `term_utils` module. Leaf module (zero dependencies). Provides shared ANSI color constants (`C_RED`, `C_GREEN`, `C_YELLOW`, `C_CYAN`, `C_BOLD`, `C_DIM`, `C_RESET`, `C_AXIS`), terminal size detection (`get_term_size` + private `stty_size`), and Bresenham line-drawing (`draw_line`). Used by `bands_plotter`, `dos_plotter`, and `poscastep_menu`.

3. **parser.f90** — `parser` module. File format parsers: `parse_cif_inline`, `parse_pdb_inline`, `parse_cell_inline`. Private helpers: `tokenize_inline`, `clean_str_inline`, `copy_str_no_quotes`. Exports `clean_element_symbol` for oxidation state stripping (e.g., "Cu0+" -> "Cu"). Handles CIF tag-value pairs and `loop_` blocks, PDB CRYST1/UNITCELL/ATOM records, CASTEP .cell `%BLOCK LATTICE_ABC`, `%BLOCK LATTICE_CART`, `%BLOCK POSITIONS_ABS`, and `%BLOCK POSITIONS_FRAC`. `parse_cell_inline` includes `compute_abc_from_cartesian` for Cartesian-to-lattice conversion. `cif_data_t` has a `positions_fractional` field set by the parser — CIF and `POSITIONS_FRAC` set it `.true.`, PDB and `POSITIONS_ABS` leave it `.false.`.

4. **cell_writer.f90** — `cell_writer` module. Generates CASTEP `.cell` file in `%BLOCK` format. Each `%BLOCK` type written by independent subroutine: `write_block_lattice_abc`, `write_block_species_pot`, `write_block_positions_abs`, `write_block_cell_constraints`, `write_block_kpoint_grid`. Phonon-specific blocks: `write_block_phonon_kpoint_mp`, `write_block_phonon_kpoint_path`, `write_block_phonon_fine_kpoint_path`, `write_block_phonon_supercell_matrix`. `write_block_positions_abs` handles both fractional (CIF) and Cartesian (PDB/.cell) input. `write_block_cell_constraints` only written when task is GEOMETRYOPTIMISATION and cell_opt_mode is ALL. `SYMMETRY_GENERATE` written when sym_source is AUTO. `KPOINTS_MP_GRID` written for GAMMA and MONKHORST_PACK schemes. `PHONON_KPOINT_MP_GRID` and `PHONON_FINE_KPOINT_MP_GRID` written for phonon tasks. Private helper `is_phonon_task` gates phonon blocks.

5. **param_writer.f90** — `param_writer` module. Generates CASTEP `.param` file in `key : value` format. Task line first, then common keywords, then task-diff keywords (16 task-specific blocks including full PHONON/PHONON+EFIELD/THERMODYNAMICS support). Phonon keywords written by shared `write_phonon_params` helper (method, energy_tol, max_cycles, dfpt_method, DOS, sum_rule, output control, Raman, etc.). `calculate_raman` also writes `raman_range_low : 0.00e+00 cm-1` and `raman_range_high : 1.00e+04 cm-1`. `calculate_born_charges` also writes `born_charge_sum_rule : true`. `write_efield_params` writes `efield_ignore_molec_modes` (CRYSTAL/MOLECULE/LINEAR_MOLECULE). vdW-DED keywords when vdW is not NONE. Force constant cutoff branches on method: SPHERICAL→`phonon_force_constant_cutoff`, CUMULANT→`phonon_force_constant_cutoff_scale`. Each keyword written by `write_kv(unit, key, value)`. CASTEP 25.12 compatibility: no `_unit` suffix keywords, no explicit toggle keywords (`thermo`, `electric_field`) — units are implicit and tasks self-enable.

6. **bands_parser.f90** — `bands_parser` module. Parses CASTEP `.bands` output files into `bands_data_t`. Public: `parse_bands_file(filename, bands, iostat, iomsg)`, `free_bands_data(bands)`. Reads header metadata (num_kpoints, num_spin, num_electrons, num_eigenvalues, fermi_energy), k-point coordinates with cumulative path distances, and eigenvalue data. Eigenvalues stored in Hartree as `eigenvalues(ie, ik, is)`.

7. **bands_plotter.f90** — `bands_plotter` module. Band structure visualization with **gap-focused analysis**. Public constants: `BANDS_MODE_ASCII=1`, `BANDS_MODE_SVG=2`. Public subroutines: `plot_bands_ascii(bands, term_w_in, term_h_in, e_center, half_range, k_pct_in, k_width_pct_in)` — interactive scrollable ANSI-color terminal plot with k-path windowing and gap info above plot, `write_bands_svg(bands, svg_file, w_px, h_px, iostat, iomsg)` — SVG vector output. Uses `term_utils` for ANSI colors, terminal detection, and line drawing. Private: `print_grid_row`, `char_type`, `type_color`, `type_char`, `band_color`, `hsl_to_rgb`, `write_svg_text`, `real2str_short`. Features: symmetric band display (first k-point copied to end), Unicode box frame, multi-symbol bands, k-point path labels (auto-detected K1..Kn), right-aligned energy labels with dynamic width, gap info in 2 lines above plot, dynamic grid height (nh−9) to fit terminal. **Gap analysis**: scans all eigenvalues in eV, finds VBM (highest below E_F) and CBM (lowest above E_F), computes indirect gap (VBM→CBM) and direct gap at VBM k-point, classifies as direct/indirect/metallic (gap < 0.005 eV). VBM marker: green ◆, CBM marker: yellow ◈. Optional `k_pct_in`/`k_width_pct_in` parameters define a k-path window for horizontal scrolling.

8. **pdos_parser.f90** — `pdos_parser` module. Parses CASTEP `.pdos_weights` / `.pdos_bin` binary files into `pdos_data_t`. Binary format: big-endian, record-delimited (4B u32 size prefix + data + 4B u32 suffix). `.pdos_bin` has two extra prefix records (f64 version + version string). Uses `access='stream', form='unformatted'` for byte-level I/O with explicit big-endian→little-endian conversion. Public: `parse_pdos_file(filename, pdos, iostat, iomsg)`, `free_pdos_data(pdos)`. Private: `read_be_u32`, `read_be_f64`, `read_payload`, `skip_record`, `be_u32_at`, `be_f64_at`, `read_u32_record`, `read_int_vec`. Output: `pdos_data_t` with header fields (total_kpoints, num_spins, num_orbitals, max_bands), orbital metadata arrays (species, ion, am[0=S,1=P,2=D,3=F]), k-point coordinates, and 4D orbital_weights(norbs, nbands, nk, nspin).

9. **dos_compute.f90** — `dos_compute` module. Gaussian smearing DOS computation. Total DOS: `DOS(E) = Σ_k Σ_i w_k · G_σ(E − ε_i)`. PDOS: `PDOS_ch(E) = Σ_k Σ_i w_k · G_σ(E−ε_i) · W_ch(i,k)` where W_ch sums orbital weights by angular momentum channel. Spin coefficient: 2 for non-polarized, 1 for spin-polarized. Gaussian cutoff at 5σ for performance. Uses `pi` from `castep_config`. Public: `compute_total_dos(bands, energy_grid, smearing, dos_result, iostat, iomsg)` — allocates `dos_result(ne, nspin)`, `compute_pdos(bands, pdos, energy_grid, smearing, pdos_result, iostat, iomsg)` — allocates `pdos_result(ne, N_CHANNELS, nspin)`. Channel constants: `CH_TOT=1, CH_S=2, CH_P=3, CH_D=4, CH_F=5`, `N_CHANNELS=5`.

10. **dos_plotter.f90** — `dos_plotter` module. DOS/PDOS visualization: interactive ASCII terminal plots + SVG + CSV export. Public constants: `DOS_MODE_ASCII=1`, `DOS_MODE_SVG=2`, `DOS_MODE_EXPORT=3`. Public: `plot_dos_ascii(energy_grid, dos_data, nspin, e_fermi, smearing, term_w_in, term_h_in, y_center_in, y_half_in, e_center_in, half_range_in)` — interactive total DOS ASCII plot with y_center/y_half axis model, `plot_pdos_ascii` — s/p/d/f multi-channel PDOS plot with legend (● s, ○ p, △ d, ▽ f), `write_dos_svg`, `write_dos_csv`, `write_pdos_csv`. Uses `term_utils` for ANSI colors, terminal detection, and line drawing. y=0 horizontal reference line with ├┤ junctions. PDOS character types extended: 'S'/'P'/'L'/'F' for s/p/d/f orbitals, 'Y'/'Z' for y=0 junctions.

11. **cli_menu.f90** — `cli_menu` module. PreCASTEP main configuration menu loop with cached `castep_config_t` state. Q returns `IO_USER_QUIT`. Task switch auto-configures DFPT→NCP19, THERMO→FD+SUPERCELL, and PHONON+EFIELD→LO/TO ON. Uses `strip_quotes` from config. Helper functions: `task_label`, `cutoff_label`, `kpoint_label`, `qpoint_label` (shows MP_GRID as "i j k" values), `scf_label`, `geom_tol_label`, `sym_label`, `sp_label`, `smearing_label`. Public: `ask_input_file` (generic path reader with `inquire` existence check). Conditionally shown items: optimizer/cell_opt_mode/geo_tolerance for GEOMETRY_OPT; phonon q-point scheme, method, fine method, energy_tol, supercell matrix for phonon tasks (items 9-14); EFIELD ignore molec modes (item 30). Extended `ask_advanced_options` accepts optional phonon/EFIELD parameters for the Advanced sub-menu (items 6-30, with DOS spacing+limit grouped under item 6, FD-only keywords gated by `.not. is_dfpt`, LO/TO toggle locked for EFIELD tasks). Phonon-specific subroutines: `ask_phonon_qpoint_scheme`, `ask_phonon_kpoint_grid`, `ask_phonon_path`, `ask_phonon_method`, `ask_phonon_fine_method`, `ask_phonon_energy_tol`, `ask_phonon_supercell_matrix`, `ask_phonon_fine_qpoint_scheme`, `ask_phonon_fine_path`, `ask_phonon_sum_rule`, `ask_phonon_dfpt_method`. K-point MONKHORST_PACK accepts comma-separated grid. EFIELD: `efield_ignore_molec_modes` (CRYSTAL/MOLECULE/LINEAR_MOLECULE).

12. **poscastep_menu.f90** — `poscastep_menu` module. PosCASTEP post-processing menu loop. Public: `run_poscastep_menu(iostat)`. Menu options: 1. Plot Band Structure, 2. Plot DOS, 3. Plot pDOS, Q. Back. Private: `handle_bands_menu` → `run_ascii_navigator` (↑↓ energy scroll, ← → k-path, +/- zoom both axes, R reset, Q quit), `handle_dos_menu` — total DOS (prompts .bands path with SAVE memory; output: ASCII/SVG/CSV; 4001-point Fermi ±20 eV grid), `handle_pdos_menu` — projected DOS (prompts file prefix, auto-derives .bands + .pdos_bin paths with SAVE memory; output: ASCII/CSV), `run_dos_navigator` (↑↓ y-pan, ←→ x-pan, +/- both-axes zoom, R reset, Q quit), `run_pdos_navigator` (same controls), `build_energy_grid`, `ensure_ext` (auto-append file extension). Gets `get_term_size` from `term_utils`. Uses `stty -icanon -echo min 1` for character input.

13. **main.f90** — `CASTEP_Suite` program. Suite top-level `do` loop with three options (1. PreCASTEP, 2. PosCASTEP, Q. Quit). PreCASTEP logic extracted into internal subroutine `run_precastep_workflow(should_exit)`: init defaults → file recognition → `run_main_menu` → parse input → compute lattice → write .cell/.param → set `should_exit=.true.` to exit program. All `stop` replaced with `return`. Coordinate system determined by `castep_config_t%cartesian_coords = .not. cif%positions_fractional` (data-driven from parser, not file extension). Helper functions: `real2str_dp`, `compute_cartesian_lattice`, `get_file_extension`, `to_lower_inline`.

### Module dependency chain

```
castep_config (leaf)
term_utils    (leaf)
  ├── parser          (config)
  ├── cell_writer     (config)
  ├── param_writer    (config)
  ├── bands_parser    (config)
  ├── pdos_parser     (config)
  ├── dos_compute     (config)
  ├── bands_plotter   (config + term_utils)
  ├── dos_plotter     (config + term_utils)
  ├── cli_menu        (config)          ─┐
  └── poscastep_menu  (config +          ├─ siblings
                        term_utils +      │  (no cross-dep)
                        bands_parser +    │
                        bands_plotter +   │
                        pdos_parser +     │
                        dos_compute +    ─┘
                        dos_plotter)
  └── main            (config + parser + cell_writer + param_writer
                        + cli_menu + poscastep_menu)
```

`cli_menu` and `poscastep_menu` are sibling modules — they share types/utilities from `castep_config` but do not depend on each other. Post-processing code lives in `poscastep_menu` + `bands_parser` + `bands_plotter` + `pdos_parser` + `dos_compute` + `dos_plotter` + `term_utils`.

## Key Design Notes

- **No external dependencies** — pure Fortran 2008, no libraries beyond stdlib.
- **Suite dual-mode architecture** — `main.f90` top-level loop dispatches to `run_precastep_workflow` (PreCASTEP input generation) or `run_poscastep_menu` (PosCASTEP post-processing). The PreCASTEP mode exits the program after successful generation; PosCASTEP returns to suite.
- **No `stop` statements** — all error/recovery paths use `return`. `IO_USER_QUIT` (-1) signals user-requested exit from sub-menus.
- **Q label convention** — sub-menus show "Q. Back" (PreCASTEP config, PosCASTEP); only the suite top-level menu shows "Q. Quit".
- **PreCASTEP exit after generation** — `run_precastep_workflow` sets `should_exit=.true.` on successful completion, causing the suite loop to exit. Errors/quits return to suite via `return`.
- **`strip_quotes` in config.f90** — single canonical implementation shared by cli_menu and poscastep_menu.
- **Bands gap analysis** — `bands_plotter` detects VBM (highest eigenvalue below E_F) and CBM (lowest above E_F), computes indirect gap (VBM→CBM regardless of k-point) and direct gap at VBM k-point, classifies as direct/indirect/metallic (gap < 0.005 eV).
- **VBM/CBM markers** — green ◆ for VBM, yellow ◈ for CBM in ASCII plot. Rendered via `char_type` types 6 and 7, `type_color` cases 6→C_GREEN and 7→C_YELLOW.
- **eV-only display** — Hartree energy units removed entirely from CLI output, SVG output, and plot code. All energies displayed in eV using `HARTREE_TO_EV = 27.211386245988_dp`.
- **Energy range** — default `e_center ± half_range` (10 eV), scrollable with ↑↓ keys, zoomable with +/-. `e_center` defaults to Fermi energy.
- **K-path windowing** — optional `k_pct_in`/`k_width_pct_in` parameters (percentages 0–1) define a k-path display window. Defaults to full path (center=0.5, width=1.0). Left/right arrow keys auto-zoom to 50% window on first press, then pan in 10% increments. `plot_bands_ascii` signature: `(bands, term_w_in, term_h_in, e_center, half_range, k_pct_in, k_width_pct_in)` where the last two are optional.
- **Gap info above plot** — band gap information displayed in 2 compact lines above the plot box. Line 1: E_F + Band Gap + type. Line 2: VBM + CBM + Dir gap (if indirect). Grid height dynamically reduced by 9 rows to keep all content visible without scrolling.
- **Resize-adaptive** — terminal size re-detected via `stty size` on every redraw. Resize the terminal and press any key to adapt.
- **Terminal mode** — uses `stty -icanon -echo min 1` (not `stty raw`). This disables line buffering and echo while preserving `opost` output processing, which is required for `\n`→`\r\n` translation. Without this, the plot displays with a staircase alignment bug.
- **Compact layout** — label-plot gap = 1 column, no "k-path" text on x-axis, maximizing grid space.
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
