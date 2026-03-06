# SCT27: SCFLY + ZELDA — NLTE Plasma Code with Non-Maxwellian EEDF

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Reference paper:** M. S. Cho et al., *Phys. Rev. E* **109**, 045207 (2024)  
DOI: [10.1103/PhysRevE.109.045207](https://doi.org/10.1103/PhysRevE.109.045207)

---

## What is SCT27?

SCT27 is a self-consistent NLTE (non-local thermodynamic equilibrium) plasma simulation code that couples:

- **SCFLY** — Super-Configuration atomic kinetics and rate equation solver (NLTE populations)
- **ZELDA** — Time-dependent Boltzmann equation solver for the Electron Energy Distribution Function (EEDF)

It was developed to study femtosecond-timescale plasma dynamics in XFEL (X-ray Free Electron Laser) experiments, where the EEDF is strongly non-Maxwellian due to:
- Photoionization producing ~keV photoelectrons
- Auger cascade producing secondary electrons
- Collisional ionization avalanche before thermalization

The code is an extension of CT27 (FLYCHK + ZELDA), replacing FLYCHK with SCFLY for improved atomic physics.

---

## Key Result (PRE 2024)

For neutral Neon (500 torr) irradiated by a 2 keV XFEL pulse:

| Model | Required density for best CSD match |
|-------|--------------------------------------|
| SCFLY alone (Maxwellian) | Ni = 2×10¹⁹ cm⁻³ |
| SCT27 (non-Maxwellian EEDF) | Ni = 10¹⁸ cm⁻³ |

The non-Maxwellian EEDF increases collisional ionization rates by ~10×, resolving an order-of-magnitude discrepancy in required density.

---

## Code Structure

```
sct27_github/
├── src/                    # All source and include files
│   ├── scmn_v3.f           # SCFLY main: LSODE rate equation driver
│   ├── scdrv.f             # Rate matrix driver (dndt)
│   ├── sclsd.f             # LSODE interface
│   ├── scslv.f             # Matrix solver
│   ├── scstr.f             # Atomic structure
│   ├── sctrn.f             # Transition rates
│   ├── zchkxc.f            # SCFLY↔ZELDA interface (runzelda, getzqs)
│   ├── zchkfly.f           # Atomic cross-section functions
│   ├── zmodules.f90        # F90 module definitions (dimensions, flags)
│   ├── zchk_stuff.f90      # F90 modules (FofE_STUFF, TIMING_STUFF, etc.)
│   ├── zcoeff_elastic.f90  # DCOEFF: collision matrix (elastic + inelastic)
│   ├── zelda.f90           # ZELDA: Boltzmann solver (FDIFF, EECOLL1, etc.)
│   ├── mainstuf            # SCFLY COMMON include: plasma parameters
│   ├── popstuf             # Population arrays
│   ├── runstuf             # Run parameters
│   ├── timestuf            # Time-stepping parameters
│   ├── xc_stuf             # Cross-section arrays
│   ├── flystuf             # FLY atomic data
│   └── kinstuf             # Kinetics parameters
├── input/
│   └── neon_2keV/          # Input files for PRE 2024 Neon case
│       ├── runfile_fig1b   # Run control file
│       ├── fe_200g_2000ev  # Initial EEDF (200 bins, 5–1995 eV, ~0 electrons)
│       ├── hvfile          # XFEL photon flux vs time (Gaussian pulse)
│       ├── his_ne_29_200step # Plasma history (density, Te fixed)
│       ├── initial_ne.dat  # Initial population (100% neutral Ne)
│       ├── atomic.data     # Atomic structure data
│       └── atomic_inp.10   # Atomic input for Ne (Z=10 index)
├── scripts/
│   ├── run_fig1b.sh        # Reproduce Fig. 1(b): XFEL pulse dynamics
│   └── run_fig1c.sh        # Reproduce Fig. 1(c): EEDF thermalization
├── docs/
│   └── physics_notes.md    # Physics details and known limitations
├── Makefile
├── LICENSE
└── README.md
```

---

## Data Flow (Per SCFLY Time Step)

```
SCFLY populations (pop[])
        │
        ▼
getzqs()  ─────────────────────────────────────┐
  Computes electron source/sink vectors:        │
  zqs[j]: electrons created (photoionization,   │
          Auger, bremsstrahlung absorption)      │
  zqr[j]: electrons lost (3-body recomb,        │
          stimulated recombination)              │
        │                                       │
        ▼                                       │
ZELDA (isfly=20 sub-steps per SCFLY step)       │
  EECOLL1() ─── builds aee[][] (e-e, Rockwood)  │
  DCOEFF()  ─── builds cf[][] matrix:           │
    elastic (Rockwood FD, e-ion)                │
    bound-bound (excitation/de-excitation)       │
    bound-free (ionization energy transfer)      │
    3-body recombination                         │
  ESOURCE() ─── applies zqs[] as source term    │
  BKSOLVE() ─── implicit solve (I - dt*cf)b = b │
  ee_Collisions() ─── explicit RK45 for e-e     │
  → fe1[i]: updated EEDF                        │
        │                                       │
        ▼                                       │
LSODE(dndt) ─── updates populations with        │
  new EEDF-weighted rate coefficients ──────────┘
```

---

## Building

### Requirements
- `gfortran` >= 9.0
- Linux or macOS (tested on Linux)

### Compile
```bash
cd src/
make
```

This produces the executable `sct27` in `src/`.

### Build options
```bash
make           # debug build (bounds checking, no optimization)
make release   # optimized build (O2)
make clean     # remove objects and executable
```

---

## Running: Reproduce PRE 2024 Results

### Fig. 1(b) — XFEL pulse dynamics

```bash
cd scripts/
bash run_fig1b.sh
```

**Physics:** Neutral Ne (Nt=2×10¹⁹ cm⁻³, 500 torr) + 2 keV XFEL (FWHM~230 fs).  
**Output:** EEDF evolving from zero electrons → photoelectron peak (~1130 eV) + Auger peak (~770 eV) → collisional cascade.

### Fig. 1(c) — EEDF thermalization to Maxwellian

To accelerate the thermalization study, set `ee_scale ~ 100` in `src/scmn_v3.f` (line ~88) and recompile:

```bash
# Edit src/scmn_v3.f: ee_scale = 100.0d0
make clean && make
cd scripts/
bash run_fig1c.sh
```

The time axis should then be divided by `ee_scale` to recover physical time. The EEDF converges to Maxwell-Boltzmann at Te ~ 100 eV by ~20 ps.

---

## Bug Fixes vs Original Code

The following corrections were applied relative to the internal development version:

| # | File | Issue | Fix |
|---|------|-------|-----|
| 1 | `zcoeff_elastic.f90` L.109 | `tgev` hardcoded to 1.0 eV (ion temperature), ignoring actual `tiev`. For Ne gas at 300 K (0.026 eV), this overestimated elastic thermalization rate by ~40×. | Removed overwrite line; now uses `tiev` correctly. |
| 2 | `zelda.f90` L.936 | `duby2=0` instead of `0.5*du` in Rockwood finite-difference scheme. Deviates from original Rockwood (1973) derivation at low energy boundary. | Restored to `duby2 = 0.5*du`. |
| 3 | `zelda.f90` L.389 | `open(..., status='old')` for `source_monitor` fails on first run when file does not yet exist. | Changed to `status='unknown'`. |

---

## Known Limitations

These are recognized physics gaps documented in the original thesis, left for future work:

1. **IPD (Ionization Potential Depression):** SCFLY computes IPD corrections but these are not fed back into ZELDA's energy thresholds. Relevant at Ni > 10²⁰ cm⁻³.

2. **Secondary electrons from ionization (`esflag`):** Currently disabled. The energy of secondary (low-energy) electrons from impact ionization is not separately tracked in the EEDF. These are handled implicitly through `getzqs` for the population rate equations.

3. **Time-splitting accuracy:** Source/sink terms (photoionization, Auger) are computed once per SCFLY step and held constant through all `isfly=20` ZELDA sub-steps. This is first-order operator splitting. Strang splitting (second-order) would improve accuracy for rapidly changing sources.

4. **Fixed `isfly=20`:** No adaptive sub-stepping in ZELDA. The stability criterion requires `dt < usqrt(i)*denergy(i)/zqr(i)` for all bins. For fine energy grids (> 300 bins) this can become marginal.

5. **Dense collision matrix:** `cf(NPTSF, NPTSF)` stores a 500×500 dense matrix. Many off-diagonal elements are structurally zero; sparse storage (e.g., band-diagonal for elastic terms) would reduce memory from ~2 MB to ~100 KB.

6. **Double Auger:** Ne⁴⁺ formation via double Auger decay is not included.

---

## Citation

If you use this code, please cite:

```bibtex
@article{Cho2024,
  author  = {Cho, M. S. and others},
  title   = {Non-equilibrium electron energy distribution in XFEL-irradiated Ne plasma},
  journal = {Phys. Rev. E},
  volume  = {109},
  pages   = {045207},
  year    = {2024},
  doi     = {10.1103/PhysRevE.109.045207}
}
```

---

## License

MIT License. See [LICENSE](LICENSE).

Original SCFLY/FLY code: Copyright Hyun-Kyung Chung and R. W. Lee.  
ZELDA Boltzmann solver: derived from ELENDIF (W. L. Morgan, Kinema Software, 1996).  
SCT27 integration and modifications: M. S. Cho, 2020–2024.
