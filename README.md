# OpenFOAM Multiphase Simulations — Oil / Water / Air / Particles

Custom OpenFOAM solvers and cases developed for a mechanical engineering thesis on turbulent, multiphase flows: oil droplets, air bubbles, and solid particles interacting in a water column. Everything here builds on the stock `interFoam` family (VOF multiphase flow) extended with Lagrangian particle tracking (`basicKinematicCloud`) and, in the most recent solver, particle–bubble attachment modelling.

Recommended platform: **OpenFOAM v2512** (the newest solver, `fourPhaseDNSFoam`, is built and tested against it; the earlier three solvers use a similar `interFoam`-family API and should also compile on v2512, though they haven't been re-verified against every intermediate version).

## Repository layout

This repo intentionally exposes the same 4 solvers and 5 cases through **two parallel views**:

- **`solver/` + `case/`** — flat, browse-by-type view. All solvers live under `solver/`, all cases under `case/`.
- **Per-simulation folders** (`FourPhaseSimulation/`, `InterFoamSimulation/`, `InterIsoFoamSimulation/`, `oilWaterAirParticleSimulation/`) — browse-by-simulation view, each folder bundling its solver + case(s) together.

Both views point at the same underlying work; pick whichever you prefer to browse. Note that folder names aren't always identical between the two views (see mapping tables below) — worth aligning at some point to avoid confusion.

```
.
├── case/                              # all cases, flat
│   ├── DNSCase01
│   ├── DNSCase02
│   ├── oilWaterAirParticle
│   ├── oilWaterParticlesInterFoam
│   └── oilWaterParticlesInterIsoFoam
├── solver/                            # all solvers, flat
│   ├── fourPhaseDNSFoam
│   ├── oilWaterAirParticleFoam
│   ├── turbulentThreePhaseInterFoam
│   └── turbulentThreePhaseInterIsoFoam
├── FourPhaseSimulation/                # = solver/fourPhaseDNSFoam + case/DNSCase0{1,2}
├── InterFoamSimulation/                # = solver/turbulentThreePhaseInterFoam + case/oilWaterParticlesInterFoam
├── InterIsoFoamSimulation/             # = solver/turbulentThreePhaseInterIsoFoam + case/oilWaterParticlesInterIsoFoam
├── oilWaterAirParticleSimulation/       # = solver/oilWaterAirParticleFoam + case/oilWaterAirParticle
└── README.md
```

### Solver name mapping

| Compiled executable | Physics | `solver/` folder | Per-simulation folder |
|---|---|---|---|
| `fourPhaseDNSFoam` | 4-phase DNS: water/oil/air VOF + falling particle cloud + rising bubble cloud + attachment model | `solver/fourPhaseDNSFoam` | `FourPhaseSimulation/fourPhaseDNSFoam` |
| `turbulentThreePhaseFoam` | 2-phase (water/oil) VOF + particle cloud | `solver/turbulentThreePhaseInterFoam` | `InterFoamSimulation/turbulentThreePhaseFoam` |
| `interisofoamSolver` | Same physics as above, isoAdvector interface capturing instead of MULES | `solver/turbulentThreePhaseInterIsoFoam` | `InterIsoFoamSimulation/turbulentThreePhaseFoam` |
| `oilWaterAirParticleFoam` | 3-phase (water/oil/air) VOF + particle cloud | `solver/oilWaterAirParticleFoam` | `oilWaterAirParticleSimulation/oilWaterAirParticleFoam` |

Note: `InterFoamSimulation/turbulentThreePhaseFoam` and `InterIsoFoamSimulation/turbulentThreePhaseFoam` are two **different** solvers (`turbulentThreePhaseFoam.C` vs. `interisofoamSolver.C` inside) that happen to share the same folder name in the per-simulation view — check `Make/files` inside each if in doubt.

### Case name mapping

| Case (physics) | `case/` folder | Per-simulation folder | Solver it runs with |
|---|---|---|---|
| Gravity-driven, pre-baked mesh, `particleInlet`/`outlet` patches | `case/DNSCase01` | `FourPhaseSimulation/DNSCase01` | `fourPhaseDNSFoam` |
| No baked mesh (needs `blockMesh`), combined `top`/`bottom` patches | `case/DNSCase02` | `FourPhaseSimulation/DNSCase02` | `fourPhaseDNSFoam` |
| Oil droplet + air bubble + particles in water, dynamic mesh | `case/oilWaterAirParticle` | `oilWaterAirParticleSimulation/oilWaterAirParticle_case` | `oilWaterAirParticleFoam` |
| Water/oil + particles, dynamic mesh | `case/oilWaterParticlesInterFoam` | `InterFoamSimulation/oilWaterParticles` | `turbulentThreePhaseFoam` |
| Same as above, isoAdvector comparison | `case/oilWaterParticlesInterIsoFoam` | `InterIsoFoamSimulation/oilWaterParticles` | `interisofoamSolver` |

Note: the two "InterFoam"/"InterIso" case folders are both named `oilWaterParticles` in the per-simulation view — only their parent folder (`InterFoamSimulation/` vs `InterIsoFoamSimulation/`) tells them apart.

## Solvers, in detail

### `fourPhaseDNSFoam`

Four-phase DNS-style solver (water/oil/air VOF + two `basicKinematicCloud`s: falling solid `particleCloud` and rising `bubbleCloud`), with a custom `attachmentModel` intended to let particles and bubbles stick together on collision. Forked from the stock `interMixingFoam` solver (its `Make/files` links directly against `.../interFoam/interMixingFoam/{immiscibleIncompressibleThreePhaseMixture,incompressibleThreePhaseMixture,threePhaseInterfaceProperties}` — no separate library build step needed, just `wmake` the solver).

Cases:
- **DNSCase01** — pre-generated mesh, correct gravity `(0 -9.81 0)`, dedicated `particleInlet` patch.
- **DNSCase02** — mesh must be generated with `blockMesh` first, combined `top`/`bottom` patches instead of separate inlet/outlet/particleInlet.

### `turbulentThreePhaseFoam` (InterFoam-based)

Two-phase (water/oil) `interFoam`-derived solver with a Lagrangian particle cloud added. Case fluids: water (ρ = 1000 kg/m³, ν = 1 cSt) and oil (ρ = 850 kg/m³, ν = 50 cSt), σ = 0.03 N/m. Case uses a `dynamicMeshDict` (adaptive/moving mesh) alongside the particle cloud. Includes post-run analysis scripts (`sim_dashboard.py`, `graph.py`).

### `interisofoamSolver` (InterIso-based)

Same physical setup as the InterFoam case above (fluid properties kept in sync deliberately for a like-for-like comparison), but built on the isoAdvector interface-capturing scheme instead of standard MULES/VOF advection. Also uses `dynamicMeshDict`. Includes turbulence post-processing (`kolmogorov_from_log.py` / `kolmogorov_scales_from_log.csv`) alongside the shared dashboard script.

### `oilWaterAirParticleFoam`

Three-phase (water/oil/air) solver with particle tracking — a small air bubble rising through an oil droplet suspended in water, with solid particles also present. Case includes `dynamicMeshDict`, `hRef`, and writes an extra `DUcDt` field (carrier-phase acceleration, used by particle force models). Comes with the most extensive post-processing setup: `oilWaterAirParticle_postprocess.py` and variants, CSV output, dashboard PNGs, and a rendered `p_rgh.avi` animation.

## Building

From inside each solver directory:

```bash
wclean
wmake
```

## Running a case

Where an `Allrun`/`Allrun_1`/`Allrun_2` script exists, use it:

```bash
cd <case_dir>
./Allrun
```

Where none exists (currently `DNSCase01` and `DNSCase02`), run manually:

```bash
cd <case_dir>
blockMesh      # only if constant/polyMesh isn't already present (e.g. DNSCase02)
decomposePar   # if running in parallel
<solverName>   # e.g. fourPhaseDNSFoam
```

## Known issues / things to check before trusting results

- **`fourPhaseDNSFoam` — buoyancy sign bug (fixed):** `pEqn.H`'s pressure-correction flux (`phig`) previously used `+ ghf*fvc::snGrad(rho)` instead of the correct `- ghf*fvc::snGrad(rho)`, which is what `UEqn.H`'s momentum predictor and the stock OpenFOAM `interFoam`/`interMixingFoam` convention both use. This reversed the buoyancy-driven flow of the carrier fluid itself. Verify your copy has `-`, not `+`.
- **DNSCase02 — zero gravity (fixed):** `constant/g` was `(0 0 0)`, which disables the `basicKinematicCloud` gravity force submodel entirely. Set to `(0 -9.81 0)` to match DNSCase01.
- **DNSCase01 — `blockMeshDict` patch conflict (not fixed, dormant):** `outlet` and `particleInlet` both claim the same top face set. Masked today because the case ships a pre-baked `constant/polyMesh` with the two patches already correctly split — running `blockMesh` from scratch will fail or silently produce a wrong mesh. Don't regenerate the mesh for this case until `blockMeshDict` is reworked.
- **`attachmentModel` is a stub:** particle/bubble velocities are synced and logged on "attachment," but particles are never actually merged or removed from the cloud (calls commented out as "read-only in OF2512" / "use cloud remove API instead"). Rise/fall physics work independently of this; the attachment *feature* doesn't yet do anything physical.
- Several headers in `fourPhaseDNSFoam` (`createPorosity.H`, `alphaEqns.H`, `createUf.H`, `porousCourantNo.H`, etc.) are leftovers from the `interMixingFoam` template and aren't `#include`d by `fourPhaseDNSFoam.C` — harmless, can be deleted for clarity.
- **Duplicate `kinematicCloudPositions`:** the InterFoam/InterIso cases have both `constant/kinematicCloudPositions` and `constant/kinematicCloud/kinematicCloudPositions`. Confirm which one the case actually reads before editing seed positions — likely only one is live and the other is a leftover copy.

## Post-processing

Each case with Python scripts (`sim_dashboard.py`, `graph.py`, `*_postprocess*.py`, `kolmogorov_from_log.py`) expects to be run from inside the case directory after a completed run, reading OpenFOAM's `postProcessing/` output and/or solver log files directly.

