DIAGNOSTIC CASE: oilWaterAirParticle_28um_noDroplet1
======================================================

PURPOSE
-------
A true grid-refinement test: same shrunk geometry as your coarse case
(5.6mm x 7mm domain, bubble R=0.35mm, real sigma12/13/23), just with the
mesh doubled in both directions (dx=dy=28um instead of 56um -- 200x250 =
50,000 cells instead of 12,500). Oil droplet 1 removed, same as the
coarse-mesh no-droplet1 test, so nothing but resolution differs from that
earlier test.

BACKGROUND
----------
On the 56um coarse mesh, with droplet 1 removed, the v3 alpha-smoothing fix
still let the bubble drift sideways into a wall and freeze there -- LEFT in
a 5-core parallel run, RIGHT in a serial rerun of the identical case. That
direction-flip (while the drift itself persisted either way) ruled out the
droplet and parallel decomposition as root causes and pointed to curvature
noise on an under-resolved bubble (~12.5 cells across).

A separate test on the ORIGINAL, pre-shrink full-scale case (bubble R=5mm,
~50 cells across, 30x100mm domain -- also with its own "droplet 1" removed)
showed no such wall-drift-and-freeze at all: the bubble stayed mobile
within the domain interior the entire run (x oscillated ~5-27mm within a
30mm-wide domain, never pinned against a wall). It did lose ~52% of its oil
mass and plateau in height once it rose far enough to contact Oil droplet
2 -- a separate, already-understood collision effect, not the runaway
sideways drift.

However, that full-scale case differs from the coarse case in much more
than cell count (absolute bubble size, domain size, relative wall
confinement -- only 3 bubble-diameters wide there vs. the coarse case's
deliberately widened 8 bubble-diameters). So it couldn't cleanly prove
resolution alone was the fix. This case closes that gap: identical
geometry to the coarse case, only the mesh is finer.

THE ONE CHANGE IN THIS CASE
----------------------------
system/setFieldsDict: the sphereToCell block that seeds Oil droplet 1 has
been commented out (not deleted -- easy to restore), exactly as in the
coarse no-droplet1 test. Nothing else (mesh resolution, transportProperties,
controlDict, droplets 2/3, bubble) has changed from your existing
oilWaterAirParticle_28um case.

WHAT YOU NEED
-------------
Your local solver must already be compiled with
threePhaseInterfaceProperties_FIXED_v3.C (the alpha-smoothing fix) -- no
new solver changes needed.

HOW TO RUN
----------
    cd oilWaterAirParticle_28um_noDroplet1
    ./Allrun

Regenerates mesh + initial fields from scratch (blockMesh + setFields), so
the setFieldsDict edit is all that's needed. This has 4x the cells of the
coarse case, so expect a longer run than the coarse tests but shorter than
the 75,000-cell full-scale case.

Upload log.solver when it finishes.

WHAT THE RESULT WOULD MEAN
---------------------------
- Bubble rises cleanly with no sideways drift/freeze: confirms resolution
  alone (not absolute physical scale) was the cause on the coarse mesh --
  the v3 fix is validated, just needs at least ~25 cells across the bubble
  to behave reliably. The coarse (56um) mesh would then be understood as
  simply too coarse for this curvature scheme, independent of anything
  else.

- Bubble still drifts into a wall and freezes: resolution alone doesn't
  fix it even at 2x, and we'd want to try 14um (4x, oilWaterAirParticle
  case already in your workspace) before concluding the smoothing
  formulation itself needs another look (e.g. the corner-curvature
  artifact, or the smoothing stencil's treatment near boundaries).
