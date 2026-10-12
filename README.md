# OpenFOAM Multiphase Particle Simulations

This repository contains custom OpenFOAM solvers, simulation cases, and post-processing tools developed for multiphase-flow research. The projects model combinations of water, oil, air, and Lagrangian solid particles, with a focus on bubble/droplet/particle interactions.

> **Project status:** This is an active research/development repository. The two-phase and three-phase cases are useful development baselines. The four-phase solver and its cases are experimental and **not working reliably yet**; do not use four-phase results for physical conclusions or validation.

## Requirements

- macOS or Linux
- OpenFOAM v2512, with the environment sourced before building or running
- A C++ compiler configured for OpenFOAM (`wmake`)
- Python 3 for optional post-processing scripts

Check that OpenFOAM is available:

```bash
which wmake
which interFoam
```

## Repository layout

```text
.
├── FourPhaseSimulation/
│   ├── DNSCase/                    # Experimental four-phase case
│   └── fourPhaseDNSFoam/           # Experimental four-phase solver
├── InterFoamSimulation/
│   ├── oilWaterParticles/          # Water/oil VOF + particle case
│   ├── oilWaterParticlesTurbulances/
│   └── turbulentThreePhaseFoam/    # InterFoam-derived solver
├── InterIsoFoamSimulation/
│   ├── oilWaterParticles/          # isoAdvector comparison case
│   ├── oilWaterParticlesTurbulances/
│   └── turbulentThreePhaseFoam/    # InterIsoFoam-derived solver
├── oilWaterAirParticleSimulation/
│   ├── oilWaterAirParticle*/       # Water/oil/air + particle cases
│   └── oilWaterAirParticleFoam/    # Three-phase custom solver
├── oilAirWaterParticleSimulation_coarse/
├── oilWaterAirParticleSimulation_large/
├── collect_diagnostics.sh
├── collect_trajectory.sh
├── comparison_dashboard.py
└── README.md
```

The repository may contain archived cases, previous attempts, compiled artifacts, logs, and macOS metadata files. These are useful for development history but are not necessarily clean or reproducible reference cases.

## Simulation overview

| Folder | Solver / approach | Modelled system | Status |
|---|---|---|---|
| `InterFoamSimulation/` | `turbulentThreePhaseFoam`, derived from `interFoam` | Water/oil VOF flow with Lagrangian particles | Development baseline |
| `InterIsoFoamSimulation/` | InterIsoFoam / isoAdvector-derived implementation | Water/oil VOF flow with particles; comparison with InterFoam formulation | Development baseline |
| `oilWaterAirParticleSimulation/` | `oilWaterAirParticleFoam` | Water/oil/air VOF flow with Lagrangian particles | Experimental, partially developed |
| `FourPhaseSimulation/` | `fourPhaseDNSFoam` | Water/oil/air VOF with particle and bubble clouds, intended particle-bubble attachment | **Not working yet** |

## Four-phase simulation status

`FourPhaseSimulation/fourPhaseDNSFoam` is intended to model a four-component system:

- Water as the carrier phase
- Oil droplets represented through VOF
- Air bubbles represented through VOF and/or a Lagrangian bubble cloud
- Solid particles represented by a Lagrangian particle cloud

The solver includes an experimental `attachmentModel` intended to describe particle-bubble attachment. However, this development is incomplete.

### Important warning

The four-phase simulation is currently **not operational as a validated simulation workflow**. It may compile or start depending on the local OpenFOAM installation and case state, but this does not mean that the implementation is physically or numerically correct.

Known limitations include:

- The particle-bubble attachment model is incomplete and should be regarded as a placeholder.
- Particle and bubble clouds are not yet coupled into a validated attachment, merging, or removal mechanism.
- Case setup, boundary conditions, and mesh-generation workflow still require verification.
- The four-phase pressure/buoyancy treatment and multiphase coupling require further testing.
- Numerical stability, mass conservation, interface behaviour, and particle/bubble trajectories have not been validated.
- Output from this solver must not be used as a research result without systematic verification and validation.

Use the working lower-complexity cases as references while developing the four-phase implementation.

## Building solvers

Source the OpenFOAM v2512 environment first. The exact command depends on the installation method. For example:

```bash
source /path/to/OpenFOAM-v2512/etc/bashrc
```

Build a solver from its directory:

```bash
cd InterFoamSimulation/turbulentThreePhaseFoam
wclean
wmake
```

Other solver directories include:

```bash
InterIsoFoamSimulation/turbulentThreePhaseFoam
oilWaterAirParticleSimulation/oilWaterAirParticleFoam
FourPhaseSimulation/fourPhaseDNSFoam
```

After a successful build, OpenFOAM normally places the executable in `$FOAM_USER_APPBIN` or the configured application-binary location.

## Running cases

Run commands from inside the intended case directory.

### Cases with an `Allrun` script

```bash
cd InterFoamSimulation/oilWaterParticles
chmod +x Allrun
./Allrun
```

Some folders use alternate script names such as `Allrun_1` or `Allrun_2`:

```bash
chmod +x Allrun_1
./Allrun_1
```

Always inspect an `Allrun` script before using it, especially in archived or experimental folders.

### Manual serial workflow

A typical OpenFOAM sequence is:

```bash
blockMesh
setFields
<solverName> | tee log.solver
```

For example, after building the InterFoam-derived solver:

```bash
cd InterFoamSimulation/oilWaterParticles
blockMesh
setFields
turbulentThreePhaseFoam | tee log.solver
```

Only run `blockMesh` if the case is designed to create its mesh from `system/blockMeshDict`. Do not overwrite a pre-existing `constant/polyMesh` unless that is intentional.

### Parallel workflow

If `system/decomposeParDict` is present:

```bash
decomposePar
mpirun -np 4 <solverName> -parallel | tee log.solver
reconstructPar
```

Replace `4` with the desired number of processes and `<solverName>` with the actual compiled executable.

## Post-processing

Several cases contain Python scripts and collection utilities, for example:

- `sim_dashboard.py`
- `graph.py`
- `comparison_dashboard.py`
- `oilWaterAirParticle_postprocess.py`
- `kolmogorov_from_log.py`
- `collect_diagnostics.sh`
- `collect_trajectory.sh`

Run post-processing from the relevant case directory where possible, because scripts commonly assume relative paths to `postProcessing/`, solver logs, or OpenFOAM time directories.

Example:

```bash
cd InterFoamSimulation/oilWaterParticles
python3 sim_dashboard.py
```

Review each script’s paths and expected input files before execution.



## Development notes

- Treat each case as version-specific: solver APIs and dictionary syntax can differ across OpenFOAM releases.
- Keep a known-working lower-complexity case before adding a new physical mechanism.
- Test solver changes incrementally: compile first, then test mesh creation, initialization, serial execution, and only then parallel execution.
- Archive or remove duplicated historical cases once a reproducible reference case has been identified.
- Before interpreting results, check Courant numbers, phase-volume conservation, residuals, interface behaviour, cloud statistics, and force/buoyancy directions.

## Code provenance and original contributions

This repository contains modifications of existing OpenFOAM solvers and uses
existing OpenFOAM libraries. It is not a set of solvers written entirely from
scratch. Development was done mainly with OpenFOAM v2512; the exact upstream
version from which each solver was originally copied has not been verified.

| Solver | Upstream basis (OpenFOAM) |
|---|---|
| `turbulentThreePhaseFoam` (InterFoamSimulation) | `interFoam` |
| `interisofoamSolver` (InterIsoFoamSimulation) | `interIsoFoam` |
| `oilWaterAirParticleFoam` | `interFoam` / `interIsoFoam`, three-phase mixture classes from `interMixingFoam` |
| `fourPhaseDNSFoam` | `interMixingFoam`; porosity terms follow Jensen et al. (2014), see `UEqnAddPorosity.H` |

Upstream copyright holders and the upstream reference are documented in the
header of each source file. Changes are documented in the file headers and the
git history.

Original contributions and modifications by Bjoern Klepsch include:
- Coupling of a Lagrangian `basicKinematicCloud` to the interFoam-, interIsoFoam-
  and interMixingFoam-based VOF solvers.
- A side-by-side setup comparing MULES/interface compression with geometric
  isoAdvector VOF on identical oil-water-particle cases.
- A three-phase (air/oil/water) particle solver, including the rewritten
  conservative alpha advection and the modified momentum/pressure equations.
- The custom `InterfaceAttachmentForce` particle force model and its
  run-time registration.
- Case setups, diagnostic scripts (`collect_diagnostics.sh`,
  `collect_trajectory.sh`) and the post-processing dashboard
  (`comparison_dashboard.py`).

Models, algorithms, or code adopted from other sources are identified in the
relevant source files and in the accompanying thesis.

## License

See [`LICENSE`](LICENSE).