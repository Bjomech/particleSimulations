# oilWaterAirParticle — complete runnable case (small / thesis case)

This archive contains two directories:

```
oilWaterAirParticle/        <- the complete case (this is the "small" case)
oilWaterAirParticleFoam/    <- the solver source (same as the last bundle, included so
                               the archive is self-contained)
```

## Which case is the small one

`oilWaterAirParticle` is the small case, and it is the one included here:

| | small (this one) | large |
|---|---|---|
| `scale` in blockMeshDict | 0.07 | 1.0 |
| domain | 5.6 × 7 × 0.07 mm | 80 × 100 × 1 mm |
| cells | 400 × 500 = 200 000 | 150 × 500 = 75 000 |
| dx = dy | 14 µm | 200 µm / 140 µm |
| air bubble diameter | 0.7 mm | 10 mm |
| oil droplets | 3, d = 0.28–0.42 mm | 3, d = 4–6 mm |
| parcels | 3 | 3 |

The `controlDict` in the previous bundle was tuned for the **large** case. The one
here is tuned for the small case — that was the missing piece.

## Build and run

```bash
# 1) build the solver (once)
cd oilWaterAirParticleFoam
wclean ; wmake            # -> $FOAM_USER_APPBIN/oilWaterAirParticleFoam

# 2) run the case
cd ../oilWaterAirParticle
./Allrun                  # blockMesh -> setFields -> decomposePar -> mpirun -> reconstructPar
```

Keep the path free of spaces and em dashes. `wmake` tolerated
`.../oilWaterAirParticleFoam — complete buildable solver/`, but several OpenFOAM
utilities do not.

`./Allclean` resets the case (removes time dirs, `0/`, logs, `*.foam`).

## What changed in the scripts you sent

Your `Allrun-2` is used as `Allrun`, with four changes:

1. **`restore0Dir`** at the top. `setFields` overwrites `0/alpha.*` in place, so a
   case that has run once no longer has pristine initial fields. Initial fields now
   live in `0.orig/` and `0/` is regenerated on every run — reproducible restarts.
2. **Error checks** on `blockMesh` / `setFields` / `decomposePar`. Previously
   `blockMesh > log 2>&1 && echo done` silently continued to the next step on
   failure; now the script stops and tells you which log to read. This is exactly
   the class of failure that produced your empty-`$application` `log.solver`.
3. **The subnormal sanitizer now checks the FoamFile header for `format ascii`**
   before rewriting a file. A byte-level regex over binary field data would corrupt
   it, so this makes the switch to `writeFormat binary` safe. The sanitizer should
   now be a no-op anyway: the subnormals came from the global
   `alpha *= sumBefore/sumAfter` rescale in the old `alphaEqn.H`, which the fixed
   solver does not do.
4. **`bubbleRise.py` runs at the end**, printing the numbers that tell you within
   seconds whether the run is physical.

`Allclean` additionally removes `0/` (regenerated from `0.orig/`) and `*.foam`.

## bubbleRise.py

Reads the reconstructed ASCII time directories and prints per output time: air and
oil volume (mm³), air and oil centroid height, max |U|, live parcel count and mean
parcel height; then a linear fit of the air centroid to give the mean rise velocity,
and the volume drift over the run. Mesh size is read from `system/blockMeshDict`, so
it keeps working if you change the resolution.

Run on **your previous small-case results** it reports:

```
mesh: 400 x 500 x 1 = 200000 cells, dx = 14.0 um, dy = 14.0 um, dz = 70.0 um

  time[s]     V_air    y_air     V_oil    y_oil     maxU  n_p   y_parc
              [mm3]     [mm]     [mm3]     [mm]    [m/s]          [mm]
    0.004     0.027   1.1788     0.021   3.6855   0.7893    3   2.6389
     0.02     0.027   2.3015     0.021   3.9734   0.5530    3   2.4986
    0.308     0.027   3.3309     0.020   4.0542   1.7642    3   1.4480
    0.328     0.027   3.3296     0.020   4.0557   4.8501    3   1.4638
    0.596     0.027   3.4627     0.019   4.1596   0.9864    3   1.8868

mean air-centroid rise velocity (linear fit): 2.84 mm/s
air volume drift: -0.00 %   oil volume drift: -8.75 %
```

Three failures visible in five lines: the bubble rises 1.18 → 2.30 mm in the first
20 ms (≈ 56 mm/s, already low) and then stalls at ≈ 3.33 mm for the remaining
0.58 s; `maxU` reaches 4.85 m/s, which for a bubble whose terminal velocity is
0.13 m/s is numerical noise from the unresolved capillary waves; and the oil loses
8.8 % of its volume. That is the signature of the unbalanced CSF force plus the 6×
oversized time step, both addressed by the fixed solver and this `controlDict`.

## Time step and end time

- Capillary (Brackbill) limit at dx = 14 µm:
  \(\Delta t_\sigma=\sqrt{(\rho_w+\rho_a)\Delta x^3/(4\pi\sigma)} = 1.77\ \mu s\).
  Your run used ≈ 11 µs. `maxDeltaT` is now 1.5 µs.
- `endTime 0.08` instead of `0.6`. A 0.7 mm bubble reaches terminal velocity within
  ~2 diameters and crosses the remaining 5.5 mm in 0.04–0.07 s.
- Cost: ≈ 53 000 steps at the ≈ 0.27 s/step your last run averaged on 5 ranks ≈ 4 h.

## Particle diameter: read this before running

`kinematicCloudProperties` now sets d_p = **14 µm** (was 35 µm), because the
Lagrangian cloud is a *point-particle* formulation and needs d_p ≤ dx = 14 µm. The
consequence, which the dictionary documents in full:

| d_p | settling velocity | Stokes number St |
|---|---|---|
| 14 µm | 0.16 mm/s | 0.005 |
| 35 µm (your value) | 1.0 mm/s | 0.032 |
| 250 µm (large case) | 33 mm/s | 0.17 |
| 500 µm (large case) | 73 mm/s | 0.69 |

At St ≈ 5 × 10⁻³ the particles are passive tracers: they cannot be captured by
inertia. Any upward transport you see comes from the bubble wake, not from capture.
**This is a scale problem, not a bug** — shrinking the geometry by 0.07 kept the
bubble physics reasonable (Eo = 0.069, Mo = 2.9 × 10⁻¹¹) but destroyed the particle
inertia. True bubble–particle attachment needs a custom `ParticleForce`
(F ∝ −c·σ·d_p·∇α_air); a sketch is in the solver bundle's README §3. Also note
`solution { coupled no; }` and the absence of `particles.SU(U)` in `UEqn.H` mean the
coupling is one-way regardless: particles feel the flow, the flow does not feel them.

## Contents

```
oilWaterAirParticle/
├── Allrun                          your Allrun-2, with the four changes above
├── Allclean                        your Allclean, + removes 0/ and *.foam
├── bubbleRise.py                   post-run quantitative check
├── 0.orig/                         pristine initial fields (U, p_rgh, alpha.air/oil/water)
├── constant/
│   ├── g, hRef, dynamicMeshDict, turbulenceProperties
│   ├── transportProperties         air 1.2/1.8e-5, oil 850/0.0425, water 1000/0.001
│   ├── kinematicCloudPositions     3 parcels at y = 1.6, 2.8, 3.5 mm
│   └── kinematicCloudProperties    gravity added, d_p = 14 µm, documented
└── system/
    ├── blockMeshDict               scale 0.07, 400 x 500 x 1
    ├── setFieldsDict               bubble r = 0.35 mm at (2.8, 1.05) mm + 3 oil drops
    ├── controlDict                 small-case time control + function objects
    ├── fvSchemes                   fixed (interfaceCompression, no steady terms)
    ├── fvSolution                  fixed (transient, no relaxation)
    └── decomposeParDict            scotch; Allrun sets numberOfSubdomains
```

Note: OpenFOAM is not installed in the environment where these files were prepared,
so nothing here has been compiled or test-run. `bubbleRise.py` was executed against
your existing result files, and the mesh/dictionary parsing is verified; the solver
build and the run itself still need to happen on your machine.
