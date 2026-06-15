# CASTEP Suite

A Fortran 2018 CLI toolkit for [CASTEP] DFT calculations — generate input files from crystal structures, and post-process data.

## What it does

| Mode | What it does |
|------|-------------|
| **PreCASTEP** | Converts `.cif` / `.pdb` / `.cell` structure files into CASTEP `.cell` + `.param` input files via an interactive menu |
| **PosCASTEP** | Post-processes CASTEP output: band structure plots with gap analysis, total DOS, projected DOS (s/p/d/f), phonon DOS, IR spectrum, Raman spectrum, and static polarizability — all as interactive ASCII terminal plots, SVG, or CSV |

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
 -1. View Crystal Structure (3D)
  0. Format Converter (.cell/.cif/.pdb)
  1. Plot Band Structure
  2. Plot DOS
  3. Plot pDOS
  4. Plot Phonon DOS
  5. Plot IR Spectrum
  6. Plot Raman Spectrum
  7. Static Polarizability
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
- SVG vector output available

### 2. Plot DOS (total density of states)

- Total DOS from `.bands` eigenvalues using Gaussian smearing: `DOS(E) = Σ_k Σ_i w_k · G_σ(E − ε_i)`
- Smearing width configurable (default 0.1 eV)
- Interactive ASCII plot with y=0 reference line
- Controls: `↑↓` y-axis pan, `← →` x-axis pan, `+/-` overall zoom (both axes), `R` reset
- Output modes: **ASCII** (terminal), **SVG**, **CSV** (for Origin etc.)
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
- **Interactive 3D**: right-drag to rotate, scroll to zoom, click/hover to select atoms
- **Atom editing**: select an atom → use IJKLUO keys (±X/±Y/±Z) to move; all symmetry-equivalent copies move in sync
- **Display modes**: 1=ball-stick, 2=space-filling, 3=wireframe; B=toggle bonds, C=toggle cell frame, P=toggle ortho/perspective
- **Step size**: `[`/`]` cycles through 0.01/0.05/0.1/0.5/1.0 Å
- **Modified structure**: on close, if atoms were moved, prompts to (1) save as CIF/PDB/cell or (2) pass directly to PreCASTEP for input generation
- JSON auto-cleaned after viewing; output files follow original input name

**Format Converter** converts between CIF, PDB, and CASTEP .cell formats.

### 8. Static Polarizability

- Computes static dielectric constant and polarizability via AIMD polarization fluctuation method
- Combines CASTEP DFPT optical dielectric tensor (ε_∞) with CP2K Berry phase dipole trajectory
- Automatic cell parameter extraction from `.castep` file
- Unwraps Berry phase polarization quantum jumps
- Window-based analysis: per-window detrend + median + linear extrapolation to W→0 for vibrational limit
- Outputs: ε_ion, ε_static, α_static (tensors + isotropic scalars)
- Input: CASTEP `.castep` file, CP2K dipole directory, temperature, MD time step

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
│   ├── bands_plotter.f90    # Band structure ASCII + SVG plotter
│   ├── pdos_parser.f90      # Binary .pdos_bin / .pdos_weights parser
│   ├── phonon_dos.f90       # .phonon parser, phonon DOS, IR & Raman spectra
│   ├── dos_compute.f90      # Gaussian smearing DOS + PDOS computation
│   ├── dos_plotter.f90      # DOS / PDOS ASCII + SVG + CSV plotter
│   ├── cli_menu.f90         # PreCASTEP configuration menu
│   ├── poscastep_menu.f90   # PosCASTEP post-processing + viewer integration
│   ├── crystal_json.f90     # JSON bridge for Rust crystal-viewer
│   ├── polarizability.f90   # Static polarizability (AIMD fluctuation method)
│   ├── drift_analysis.f90   # Drift rate diagnostics (not compiled, dev artifact)
│   └── main.f90             # Entry point, suite menu dispatcher
└── crystal-viewer/          # Rust/Bevy 3D viewer subproject
    ├── Cargo.toml
    └── src/
        ├── main.rs           # App setup, camera, movement, display modes
        ├── crystal.rs        # Data types, lattice math, cell expansion
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
