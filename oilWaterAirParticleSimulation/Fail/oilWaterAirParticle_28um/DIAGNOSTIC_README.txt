28um INTERMEDIATE-RESOLUTION DIAGNOSTIC VARIANT
=================================================

Purpose: pinpoint where along the resolution spectrum the bubble's behavior
switches from "rises correctly" to "frozen in place."

Changed from the fine (14um) baseline case:
  - system/blockMeshDict: cell counts halved in both x and y
    (400 -> 200, 500 -> 250), giving dx = dy = 28um instead of 14um.
    Same physical domain size (80mm x 100mm x 1mm), same geometry.

Everything else is UNCHANGED from the original fine-mesh baseline case:
  - constant/transportProperties: real sigma (sigma12=0.03, sigma13=0.07,
    sigma23=0.05), real rho/mu values.
  - system/fvSolution: original relaxation factors (p_rgh=0.3, U=0.7),
    original nOuterCorrectors=2.
  - system/setFieldsDict: unchanged (geometric regions, resolution-independent).
  - system/controlDict: endTime=0.6 (same as the coarse 56um run, for direct
    comparability), writeInterval=0.004.
  - system/decomposeParDict: numberOfSubdomains 5, method scotch (unchanged;
    scotch decomposition works for any cell count).

Comparison points so far (all with the current canonical fixed solver):
  - 14um (fine):   bubble frozen at y~0.7-0.9mm across every test variant
                   (real sigma, 20% sigma, zero sigma, relaxation-fix).
  - 56um (coarse): bubble rises cleanly, y: 0.81mm -> 6.24mm over t=0-0.6s.
  - 28um (this case): unknown -- this is the test.

To revert to the fine-mesh baseline: change blockMeshDict blocks back to
(400 500 1).

FIX applied (v2 of this case, replaces the earlier broken zip):
  The copied 0/alpha.air, 0/alpha.oil, 0/alpha.water still held a fully
  resolved "nonuniform List<scalar> 200000(...)" field -- a leftover
  per-cell list sized for the FINE mesh's 200000 cells (400x500), not this
  mesh's 50000 cells (200x250). setFields/decomposePar choke trying to read
  a mismatched list length ("Size 200000 is not equal to the expected
  length 50000"). Reset all three to plain "uniform" values (air=0, oil=0,
  water=1) -- setFieldsDict's defaultFieldValues + sphereToCell regions use
  absolute physical coordinates, so this is fully resolution-independent
  and setFields will correctly re-populate the bubble/droplet regions on
  any mesh. Also removed the stale 0.orig/C, Cx, Cy, Cz cell-center files
  (leftover diagnostic output from the fine mesh, not read by setFields or
  the solver, just unused clutter at the wrong size).

Run with the usual Allrun script; upload log.solver when done.
