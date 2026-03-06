# SCT27 Physics Notes

## Key Variables

| Variable | Location | Meaning |
|----------|----------|---------|
| `fe(i)` | `FofE_STUFF` | EEDF: Ne × f(E), normalized so ∫ fe(E)√E dE = Ne |
| `fe0(i)` | `FofE_STUFF` | EEDF at start of SCFLY time step (input to ZELDA) |
| `fe1(i)` | `FofE_STUFF` | EEDF at end of ZELDA sub-stepping (output) |
| `b(i)` | `zelda.f90` FDIFF | Working variable: b = fe0 × √E × Ne during ZELDA solve |
| `zqs(j)` | `xc_stuf` | Electron source rate [/cc/eV/s] per energy bin j (photoion, Auger) |
| `zqr(j)` | `xc_stuf` | Electron sink rate coefficient [1/s] per bin j (3-body recomb) |
| `ee_scale` | `TIMING_STUFF` | Multiplier on e-e collision rate (default 1.0; set ~100 for Fig.1c) |
| `tgev` | `zcoeff_elastic.f90` | Ion temperature [eV] used in elastic scattering (= tiev, NOT hardcoded) |
| `isfly` | `scmn_v3.f` | Number of ZELDA sub-steps per SCFLY step (default 20) |

---

## b(i) Definition and Normalization

ZELDA internally works with the quantity:

```
b(i) = f(E_i) × √(E_i) × Ne     [units: /cc/eV]
```

This is initialized from `fe0` and at the end of FDIFF converted back:

```fortran
b(1:n) = b(1:n) / usqrt(1:n)    ! → f(E) × Ne
```

INTGRL integrates:
```
Ne = ∫ b(E) dE   (when b = f√E × Ne)
<E> = (2/3) × ∫ b(E) × E^1.5 dE / Ne
```

---

## Fig. 1(c) Reproduction Method

The ~20 ps thermalization in Fig. 1(c) cannot be directly simulated because:
- Required time step: ~0.01 fs (stability: `dt < u(i)^(1/2) × dE / zqr(i)`)
- Physical duration: ~20 ps = 20,000 fs
- Required steps: ~2,000,000 (impractical run time)

**Workaround used in the paper:**
1. Run the XFEL pulse simulation (<600 fs) → obtain non-Maxwellian EEDF
2. Set `ee_scale = 100` → artificially accelerates e-e collisions by 100×
3. Simulate for ~200 fs (100× less, same number of steps)
4. Rescale time axis linearly: `t_plot = t_sim × ee_scale`

**Physical justification:** e-e collisions (Coulomb) are the dominant thermalization mechanism after the XFEL pulse ends. Scaling them uniformly by a constant preserves the shape of the relaxation but compresses the timescale.

**Variable in code:** `ee_scale` in `src/scmn_v3.f` (line ~88) and `TIMING_STUFF` module.  
Edit and recompile to change. The factor is applied in `derivs()` (zelda.f90) as `alf = alf * ee_scale`.

---

## Elastic Scattering: Rockwood Finite-Difference Scheme

The electron-ion elastic collision term follows Rockwood (1973):

```
cf(i,i)   = -(a_i + b_i)
cf(i,i+1) =  b_{i+1}
cf(i,i-1) =  a_{i-1}
```

where `a_k`, `b_k+1` are computed from the momentum-transfer cross section Q(E):

```fortran
ak(z,q)   = const1 * alpha * (d2 / (const2 * sqrt(z) * q)) * (z + duby4) + ...
bkp1(z,q) = const1 * alpha * (d2 / (const2 * sqrt(z) * q)) * (z - duby4) + ...
```

`tgev` (ion temperature) appears in the thermal correction terms. Using `tgev = tiev` (actual ion temperature) is essential for correct thermalization at low density.

---

## Coulomb Logarithm (Zollweg & Liebermann 1987)

Both e-i elastic (DCOEFF/QMOMNTM) and e-e (derivs) use the same Coulomb log:

```
λ = √[ (743 √(Te/Ne))² + (4π Ne/3)^{-2/3} ] / (Ze × 4.8e-8 / Te)
ln Λ = log( √(1 + 1.4 λ²) )
```

This interpolates between the quantum limit (low Te) and classical limit (high Te/low Ne).

---

## Time-Splitting Structure

Per SCFLY step (tlast → tnext):

```
1. getzqs()            ← compute zqs, zqr from current populations
2. ZELDA loop (isfly sub-steps):
     a. b += dt × Source    (explicit, Source fixed = zqs/dE)
     b. b -= dt × Sink      (explicit, Sink ∝ zqr × b(i), updated each step)
     c. BKSOLVE(b)          (implicit solve for elastic + inelastic cf)
     d. ee_Collisions(b)    (explicit RK45 for e-e)
3. fe = fe1               ← pass new EEDF back to SCFLY
4. LSODE(dndt)            ← update populations with new EEDF rates
```

This is first-order Lie-Trotter splitting. The source term uses the population at the *start* of the SCFLY step (held constant through all isfly sub-steps), which introduces an O(dt_SCFLY) splitting error.
