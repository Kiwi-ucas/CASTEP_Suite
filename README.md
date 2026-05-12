# CASTEP Suite

A Fortran 2008 CLI toolkit for [CASTEP] DFT calculations — generate input files from crystal structures, and post-process band structures, DOS, and projected DOS.

## What it does

| Mode | What it does |
|------|-------------|
| **PreCASTEP** | Converts `.cif` / `.pdb` / `.cell` structure files into CASTEP `.cell` + `.param` input files via an interactive menu |
| **PosCASTEP** | Post-processes CASTEP output: band structure plots with gap analysis, total DOS, and projected DOS (s/p/d/f) — all as interactive ASCII terminal plots, SVG, or CSV |

## Requirements

- **gfortran** >= 7.0 (Fortran 2008)
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

### Menu

```
  ==================================
              PreCASTEP
  ==================================
  -2. Advanced option
  -1. Spin_polarized : false
   0. Generate .cell and .param files
   1. Task : SINGLEPOINT
   2. XC_functional : PBE
   3. Cutoff energy : 400.0 eV
   4. vdW correction : NONE
   5. Pseudopotential : C19MK2
   6. K-point scheme : GAMMA
   7. SCF tolerance : 1e-5
   8. Symmetry : NONE
   Q. Back
```

Key options:

- **Task types**: 17 options including SINGLEPOINT, GEOMETRYOPTIMISATION, CELLOPTIMISATION, MOLECULAR_DYNAMICS, PHONON, ElectronicSpectroscopy, and more
- **XC functionals**: PBE, PBEsol, HSE06, PBE0, r2scan
- **vdW corrections**: D3, D3-BJ, D4
- **Pseudopotentials**: NCP19, C19MK2, SOC19 (auto-enables spin polarization)
- **K-point schemes**: Gamma, Monkhorst-Pack
- **Optimizer** (geometry tasks): BFGS, LBFGS, CG

### Output

Generates a `.cell` file in `%BLOCK` format and a `.param` file in `key : value` format. Always writes `LATTICE_ABC` regardless of input format.

## PosCASTEP — post-processing

```
  ================================
            PosCASTEP
  ================================
  1. Plot Band Structure
  2. Plot DOS
  3. Plot pDOS
  Q. Back
```

### 1. Plot Band Structure

- Parses CASTEP `.bands` files
- Interactive ASCII terminal plot with ANSI colors
- **Gap analysis**: auto-detects VBM (green ◆) and CBM (yellow ◈), computes direct/indirect band gap
- **Terminal-adaptive**: re-detects size on every redraw
- Unicode box frame, multi-symbol bands (● ○ □ △ ▽)
- k-point path labels auto-detected from direction changes
- Controls: `↑↓` scroll energy, `← →` scroll k-path, `+/-` zoom, `R` reset, `Q` quit
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

## Project structure

```
CASTEP_Suite/
├── Makefile
├── README.md
├── CLAUDE.md
└── src/
    ├── config.f90           # Types, constants, physical parameters
    ├── parser.f90           # CIF / PDB / .cell file parsers
    ├── cell_writer.f90      # CASTEP .cell file generator (%BLOCK format)
    ├── param_writer.f90     # CASTEP .param file generator (key-value)
    ├── bands_parser.f90     # CASTEP .bands file parser
    ├── bands_plotter.f90    # Band structure ASCII + SVG plotter
    ├── pdos_parser.f90      # Binary .pdos_bin / .pdos_weights parser
    ├── dos_compute.f90      # Gaussian smearing DOS + PDOS computation
    ├── dos_plotter.f90      # DOS / PDOS ASCII + SVG + CSV plotter
    ├── cli_menu.f90         # PreCASTEP configuration menu
    ├── poscastep_menu.f90   # PosCASTEP post-processing menu
    └── main.f90             # Entry point, suite menu dispatcher
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
| `-std=f2008` | Fortran 2008 standard |
| `-fimplicit-none` | Require explicit declarations |
| `-Wall -Wextra` | Comprehensive warnings |
| `-O2` | Optimization |
| `-g` | Debug symbols |
