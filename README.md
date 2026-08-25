# CASTEP Suite

A Fortran 2018 CLI toolkit for [CASTEP] DFT calculations — generate input files from crystal structures, and post-process data.

## What it does

| Mode | What it does |
|------|-------------|
| **PreCASTEP** | Converts `.cif` / `.pdb` / `.cell` structure files into CASTEP `.cell` + `.param` input files via an interactive menu |
| **PosCASTEP** | Post-processes CASTEP output: band structure plots with gap analysis, total DOS, projected DOS (s/p/d/f), phonon DOS, IR spectrum, Raman spectrum, and static polarizability — all as interactive ASCII terminal plots or CSV export |

## Requirements

- **gfortran** >= 7.0 (Fortran 2018)
- **cargo/rustc** >= 1.75 (for crystal-viewer; Fortran-only builds work without it)
- **Make**
- No external libraries

## Quick start

```bash
make
./CASTEP_Suite
```

```
  ==================================
            CASTEP Suite
  ==================================
  1. PreCASTEP  (generate CASTEP input files)
  2. PosCASTEP  (post-process CASTEP output)
  Q. Quit
  Select mode:
```

## PreCASTEP — input file generation

Converts crystal structure files into CASTEP input. Supports three input formats:

| Format | Extension | Source |
|--------|-----------|--------|
| CIF / mmCIF | `.cif` | Crystallography databases |
| PDB | `.pdb` | Protein Data Bank |
| CASTEP cell | `.cell` | Existing CASTEP files |

Output files are auto-named from the input stem and task type (e.g. `Cu_Phonon.cell`, `Cu_Phonon.param`).

### Menu

```
  ================================
             PreCASTEP
  ================================
  CIF: Cu.cif
 -3. Nonlinear optics       (NONE)       ← EFIELD / Phonon+Efield only
 -2. Advanced option
 -1. Spin_polarized : false
  0. Generate InputFile
  1. Task                   (Energy)
  2. XC functional          (PBE)
  3. Cutoff energy (eV)     (400)
  4. vdW correction         (NONE)
  5. Pseudopotential        (C19MK2)
  6. K-point                (GAMMA)
  7. SCF tolerance          (1e-5)
  8. Symmetry               (NONE)
  9. Phonon q-point scheme   (1 1 1)     ← phonon tasks only
 10. Phonon method           (DFPT)
 11. Phonon fine method      (INTERPOLATE)
 12. Phonon energy tol       (  1.0E-05)
 14. Phonon fine q-point     (1 1 1)
  Q. Back
```

### Task types (8 developed, 8 stubs)

| # | Task | Status |
|---|------|--------|
| 1 | Energy | Done |
| 2 | GeometryOptimisation | Done |
| 3 | ElectronicSpectroscopy | Done |
| 4 | Phonon | Done |
| 5 | Phonon+Efield | Done |
| 6 | Efield | Done |
| 7 | Thermodynamics | Done |
| 8 | CINEB (NEB+Climbing Image) | Done |
| 9-16 | MolecularDynamics … EpCoupling | Stub |

### Key options

- **XC functionals**: PBE, PBEsol, HSE06, PBE0, r2scan
- **vdW corrections**: NONE, D3, D3-BJ, D4
- **Pseudopotentials**: NCP19, C19MK2, SOC19 (auto-enables spin polarization)
- **K-point schemes**: Gamma, Monkhorst-Pack
- **Geometry optimizers**: BFGS, LBFGS, CG
- **Geo tolerance**: COARSE, MEDIUM, FINE, EXTREME
- **Symmetry**: NONE, AUTO

### Phonon support

Full phonon calculation support with Finite Displacement (FD) and DFPT methods:

- **Q-point sampling**: MP_GRID or custom k-point PATH
- **Phonon DOS**: Gaussian smearing with configurable spacing and limit
- **Sum rule**: NONE or RECIPROCAL
- **Force constants**: write to disk, cutoff with SPHERICAL or CUMULANT method
- **Born charges** and **Raman** intensities
- **LO/TO splitting** (auto-enabled for Phonon+Efield)
- **Fine method**: INTERPOLATE or SUPERCELL with independent q-point sampling

### EFIELD support

- DFPT polarizability with configurable max cycles, energy tolerance, convergence window
- Frequency-dependent response: spacing, oscillator Q, ion permittivity
- **Nonlinear optics**: CHI2 calculation
- Molecular mode exclusion: CRYSTAL(3), MOLECULE(6), LINEAR_MOLECULE(5)
- Phonon+Efield combined task with full phonon + EFIELD parameter sets

### CINEB transition state search

NEB + Climbing Image transition state search (`TASK : TRANSITIONSTATESEARCH`). Requires 3 structure files:

- **Reaction path**: reactant + product + intermediate guess (separate `.cif`/`.pdb`/`.cell` files)
- **Climbing Image**: always ON (hardcoded), drives the highest-energy image to the exact saddle point
- **Path images**: odd number enforced (auto-corrected from even)
- **Convergence tolerance**: reuses the same 4-level system as geometry optimization (COARSE/MEDIUM/FINE/EXTREME)
- **Product/intermediate validation**: atom counts must match reactant (hard CASTEP requirement)

| Parameter | Options | Default |
|-----------|---------|---------|
| CINEB max images | >= 3, odd | 11 |
| CINEB spring constant | Any positive number | 0.1 eV/Å² |
| CINEB tangent mode | NONE, BISECT, HIGH_E, SPLINE | SPLINE |
| CINEB NEB method | TPSD, FIRE, ODE12R | ODE12R |
| CINEB max iterations | Any positive integer | 50 |
| TS tolerance | COARSE, MEDIUM, FINE, EXTREME | MEDIUM |

### Advanced options

Additional parameters accessible via `-2. Advanced option`:

| Group | Items |
|-------|-------|
| SCF | Smearing, Max SCF cycles, Convergence window |
| Electronic | Calculate ELF, Calculate EDD |
| Phonon | DOS, sum rule, finite displacement, max cycles, DFPT method, force constants, dynamical matrix, LO/TO, cutoff, max CG steps, k-point symmetry, Born charges, Raman |
| EFIELD | DFPT method, max cycles, energy tol, conv window, freq spacing, oscillator Q, ion permittivity, ignore molec modes |

### Output

Generates a `.cell` file in `%BLOCK` format and a `.param` file in `key : value` format. Always writes `LATTICE_ABC` regardless of input format. CASTEP 25.12 compatible.

## PosCASTEP — post-processing

All interactive plots use the alternate screen buffer — no scrollback pollution. Type `q` at any file path prompt to cancel and return to the PosCASTEP menu.

```
  ================================
            PosCASTEP
  ================================
 -2. Phonon Mode Visualization
 -1. View Crystal Structure (3D)
  0. Format Converter (.cell/.cif/.pdb)
  1. Plot Band Structure
  2. Plot DOS
  3. Plot pDOS
  4. Plot Phonon DOS
  5. Plot IR Spectrum
  6. Plot Raman Spectrum
  7. Static Polarizability
  8. Thermodynamics
  9. PES Scan
 10. Frozen Phonon Scan
  Q. Back
```

### 1. Plot Band Structure

- Parses CASTEP `.bands` files
- Interactive ASCII terminal plot with ANSI colors
- **Gap analysis**: auto-detects VBM (green ◆) and CBM (yellow ◈), computes direct/indirect band gap
- **Terminal-adaptive**: re-detects size on every redraw
- Unicode box frame, multi-symbol bands (● ○ □ △ ▽)
- k-point path labels auto-detected from direction changes
- Controls: `↑↓` scroll energy, `← →` scroll k-path, `+/-` zoom both axes, `R` reset, `Q` quit

### 2. Plot DOS (total density of states)

- Total DOS from `.bands` eigenvalues using Gaussian smearing: `DOS(E) = Σ_k Σ_i w_k · G_σ(E − ε_i)`
- Smearing width configurable (default 0.1 eV)
- Interactive ASCII plot with y=0 reference line
- Controls: `↑↓` y-axis pan, `← →` x-axis pan, `+/-` overall zoom (both axes), `R` reset
- Output modes: **ASCII** (terminal), **CSV** (for Origin etc.)
- Path memory: Enter to reuse last `.bands` path

### 3. Plot pDOS (partial density of states)

- Requires `.bands` + `.pdos_bin` (or `.pdos_weights`) from the same CASTEP calculation
- Enter file prefix once — auto-finds `<prefix>.bands` + `<prefix>.pdos_bin`
- Parses CASTEP's binary PDOS weights format (big-endian, record-delimited)
- Decomposes DOS into **s**, **p**, **d**, **f** angular momentum channels
- Multi-curve plot with legend: `● s  ○ p  △ d  ▽ f`
- Controls: same as total DOS (`↑↓ ← → +/- R Q`)
- Output modes: **ASCII**, **CSV** (6 columns: Energy, Total, s, p, d, f)

### 4. Plot Phonon DOS (phonon density of states)

- Parses CASTEP `.phonon` files for vibrational frequencies
- Gaussian smearing with configurable width (default 5 cm⁻¹)
- Interactive ASCII plot with frequency axis (cm⁻¹)
- Controls: same as total DOS (`↑↓ ← → +/- R Q`)
- Output modes: **ASCII**, **CSV**

### 5. Plot IR Spectrum (infrared absorption)

- IR absorption intensities from `.phonon` file (Gamma point)
- Gaussian broadening with configurable width
- Interactive ASCII plot, CSV export
- Controls: same as total DOS

### 6. Plot Raman Spectrum (Raman scattering)

- Raman scattering activities from `.phonon` file (Gamma point)
- Gaussian broadening with configurable width
- Interactive ASCII plot, CSV export
- Controls: same as total DOS

### -1. View Crystal Structure (3D) & Format Converter (0)

**View Crystal Structure** launches a standalone 3D crystal viewer (Rust/Bevy):

- Auto-parses CIF, PDB, or CASTEP .cell files
- **Full unit cell display**: automatic expansion from asymmetric unit (applies ±1 fractional translations)
- **Interactive 3D**: right-drag or arrow keys (←↓↑→) to rotate, scroll to zoom, click/hover to select atoms
- **Atom editing**: select an atom → use HKUMIN keys (±X/±Y/±Z) to move; all symmetry-equivalent copies move in sync
- **Display modes**: 1=ball-stick, 2=space-filling, 3=wireframe; A=toggle cell axes, B=toggle bonds, C=toggle cell frame, P=toggle ortho/perspective, R=reset camera
- **Right panel controls**: Rotation Angle DragValue (1°–90°, default 45°), Move Step DragValue (0.01–10 Å, default 0.5 Å)
- **Cell axes**: A toggles red X/green Y/blue Z arrows from cell origin along lattice vectors (1.5× length)
- **Modified structure**: on close, if atoms were moved, prompts to (1) save as CIF/PDB/cell or (2) pass directly to PreCASTEP for input generation
- JSON auto-cleaned after viewing; output files follow original input name

**Format Converter** converts between CIF, PDB, and CASTEP .cell formats.

### 7. Static Polarizability

- Computes static dielectric constant and polarizability via AIMD polarization fluctuation method
- Combines CASTEP DFPT optical dielectric tensor (ε_∞) with CP2K Berry phase dipole trajectory
- Automatic cell parameter extraction from `.castep` file
- Unwraps Berry phase polarization quantum jumps
- Window-based analysis: per-window detrend + median + linear extrapolation to W→0 for vibrational limit
- Outputs: ε_ion, ε_static, α_static (tensors + isotropic scalars)
- Input: CASTEP `.castep` file, CP2K dipole directory, temperature, MD time step

### -2. Phonon Mode Visualization

- Parses phonon eigenvectors from `.phonon` files and Born effective charges from `.castep` files
- Decomposes each phonon mode to per-atom contributions using the Born charge formalism: p_m = Σ Z*(κ)·u(κ,m)
- Launches the 3D crystal viewer with colored displacement arrows showing vibration patterns
- **Displacement arrows**: cone-tipped cylinders colored by atom contribution (green=low → red=high)
- **Mode selector**: choose mode by index or jump to the highest-IR-intensity mode
- Supports `.phonon` files with or without Raman activity column
- Born charge decomposition is optional — arrows still render without it (in white)

### 9. PES Scan (Potential Energy Surface)

2D/3D potential energy surface scan with symmetry reduction and isosurface visualization:

**2D PES**: single-atom constrained scan on a plane → colored energy surface mesh in viewer.

**3D PES** with **orbit mapping** symmetry reduction:
- Automatically detects space group symmetry from CIF files
- Lays full uniform grid over the unit cell, groups symmetry-equivalent points into orbits
- Each orbit contributes exactly one irreducible representative — the minimal CASTEP scan set
- e.g., Li6PS5Cl (F-43m) with 5×5×5 grid: 125 → 9 irreducible points (93% reduction)
- Generates `irred_NNNNN` subdirectories (or `grid_III_JJJ_KKK` for non-symmetry scans)
- Forward orbit expansion fills the full-cell output grid at 100% coverage
- **Collection**: auto-detects directory format, reads `irred_coords.dat` sidecar for mapping
- **3D viewer**: MC isosurface (4), volume render (5), slice planes (6), fixed-radius sphere (7), radial-stationary migration surface (8); `-/+` adjust isovalue; semi-transparent rendering
- **Export**: "Export PLY" button (next to Render, mode 8 with a loaded cube, or E key) writes the migration shell as `migration_surface.ply` (jet-colored mesh for Blender/external rendering)
- **Render (WYSIWYG)**: bottom-right "Render" button → draggable dialog with resolution (up to 8192²), MSAA (Off/2x/4x), format (PNG/TIFF), and live scene parameters: key/fill/ambient light intensity, key-light shadows, atom roughness + metallic, tonemapping (TonyMcMapface/ACES/AgX/Reinhard), background color (RGB sliders + Black/White/Light-gray/Dark-gray/Viewer-default presets). Every edit applies to the viewer in real time — the export is rendered by the SAME Bevy PBR pipeline you see on screen (any display mode), saved as `render.png`/`render.tiff`; "Reset all" restores the defaults

**2D PES**:
- Select mobile atom, scan plane (XY/XZ/YZ), fractional coordinate range, and grid size
- Generates `grid_001_001/` … `grid_Nx_Ny/` subdirectories with CASTEP input + constraints
- Result collection parses `.castep` files, writes `scan.cube` with energy grid
- 3D viewer: colored surface mesh (jet colormap), S toggles visibility, -/+ adjusts color clip

### 10. Frozen Phonon Scan

Generates CASTEP single-point inputs for frozen-phonon total-energy scans directly from a `.phonon` file:

- **Input**: `.phonon` file with `Phonon Eigenvectors` block (run the phonon job with `PHONON_WRITE_EIGENVECTORS : true`)
- **Gamma only**: validates and displays the first q-point; only `q = (0,0,0)` is accepted
- **Mode selector**: lists mode index + frequency; user chooses a mode
- **Displacement convention**: raw mass-weighted eigenvectors `e/√m` are normalized so that `max_i |u_i| = 1`, then scaled to `±0.1`, `±0.2`, `±0.5 Å` (6 structures, no 0-point)
- **Preview**: prints every atom's actual `|u|` and new fractional coordinates for all six displacements
- **PreCASTEP handoff**: structure/lattice/species are taken from the `.phonon` header; task is locked to `SINGLEPOINT` (`ENERGY`), all other parameters are user-selectable
- **Output**:
  ```text
  frozen_phonon/<stem>_mode<M>/
  ├── d-0.50/<stem>_mode<M>_d-0.50.cell + .param
  ├── d-0.20/...
  ├── d-0.10/...
  ├── d0.10/...
  ├── d0.20/...
  └── d0.50/...
  ```


## Project structure

```
CASTEP_Suite/
├── Makefile
├── README.md
├── CLAUDE.md
├── src/
│   ├── config.f90           # Types, constants, physical parameters
│   ├── term_utils.f90       # ANSI colors, terminal size, Bresenham, alt screen
│   ├── parser.f90           # CIF / PDB / .cell file parsers
│   ├── cell_writer.f90      # CASTEP .cell file generator (%BLOCK format)
│   ├── param_writer.f90     # CASTEP .param file generator (key-value)
│   ├── bands_parser.f90     # CASTEP .bands file parser
│   ├── bands_plotter.f90    # Band structure ASCII plotter
│   ├── pdos_parser.f90      # Binary .pdos_bin / .pdos_weights parser
│   ├── phonon_dos.f90       # .phonon parser, phonon DOS, IR & Raman spectra
│   ├── phonon_modes.f90     # Eigenvector/Born charge parser, mode decomposition
│   ├── dos_compute.f90      # Gaussian smearing DOS + PDOS computation
│   ├── dos_plotter.f90      # DOS / PDOS ASCII + CSV plotter
│   ├── cli_menu.f90         # PreCASTEP configuration menu
│   ├── poscastep_menu.f90   # PosCASTEP post-processing + viewer integration
│   ├── crystal_json.f90     # JSON bridge for Rust crystal-viewer
│   ├── polarizability.f90   # Static polarizability (AIMD fluctuation method)
│   ├── pes.f90              # Unified PES scan (2D/3D), orbit mapping, symmetry expansion
│   ├── drift_analysis.f90   # Drift rate diagnostics (not compiled, dev artifact)
│   └── main.f90             # Entry point, suite menu dispatcher
└── crystal-viewer/          # Rust/Bevy 3D viewer subproject
    ├── Cargo.toml
    └── src/
        ├── main.rs           # App setup, camera, movement, PES surface, display modes
        ├── crystal.rs        # Data types, lattice math, cell expansion
        ├── pes.rs            # PES data types, surface mesh generation, colormap
        ├── picking.rs        # MVP-projection atom picking/highlighting
        ├── ui.rs             # egui panels (atom list, info, toolbar)
        └── resources.rs      # Periodic table (radii, CPK/Jmol colors)
```

## Build options

| Command | Description |
|---------|-------------|
| `make` | Release build (gfortran, `-O2`) |
| `make debug` | Debug build (`-O0 -fcheck=all -fbacktrace`) |
| `make run` | Build and run |
| `make clean` | Remove `obj/` and binary |

## Compilation flags

| Flag | Purpose |
|------|---------|
| `-std=f2018` | Fortran 2018 standard |
| `-fimplicit-none` | Require explicit declarations |
| `-Wall -Wextra` | Comprehensive warnings |
| `-O2` | Optimization |
| `-g` | Debug symbols |

## License

MIT
