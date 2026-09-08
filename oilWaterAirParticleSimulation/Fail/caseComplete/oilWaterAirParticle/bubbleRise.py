#!/usr/bin/env python3
"""
bubbleRise.py - quantitative post-check for the oilWaterAirParticle case.

Reads the reconstructed ASCII time directories of this case and prints, per
written time:

    V_air      air volume            [mm^3]  -> must stay constant (VOF mass
                                               conservation; this is the single
                                               most useful sanity check)
    y_air      air centroid height   [mm]    -> its slope is the rise velocity
    V_oil      oil volume            [mm^3]  -> must also stay constant; the old
                                               solver lost 27% of it
    y_oil      oil centroid height   [mm]
    maxU       max velocity magnitude [m/s]  -> if this blows up, dt is too big
    n_parc     number of live parcels
    y_parc     mean parcel height    [mm]    -> rising with the bubble = capture

and finally the mean rise velocity from a linear fit of y_air(t).

Usage:
    python3 bubbleRise.py            # auto-detects mesh + case dir
    python3 bubbleRise.py <caseDir>

Requires numpy only. Mesh dimensions are read from system/blockMeshDict, so it
keeps working if you change the resolution or the domain size.
"""

import os
import re
import sys
import glob

try:
    import numpy as np
except ImportError:
    sys.exit("bubbleRise.py needs numpy:  pip3 install numpy")


# --------------------------------------------------------------------------- #
# blockMeshDict -> domain size and cell counts
# --------------------------------------------------------------------------- #
def read_mesh(case):
    txt = open(os.path.join(case, "system", "blockMeshDict")).read()
    txt = re.sub(r"/\*.*?\*/", " ", txt, flags=re.S)
    txt = re.sub(r"//[^\n]*", " ", txt)

    m = re.search(r"(?:convertToMeters|scale)\s+([-\d.eE+]+)\s*;", txt)
    scale = float(m.group(1)) if m else 1.0

    verts = re.search(r"vertices\s*\((.*?)\)\s*;", txt, flags=re.S).group(1)
    pts = np.array(
        [[float(v) for v in t.split()] for t in re.findall(r"\(([^()]*)\)", verts)]
    ) * scale

    m = re.search(r"hex\s*\([^)]*\)\s*\(\s*(\d+)\s+(\d+)\s+(\d+)\s*\)", txt)
    if m is None:
        m = re.search(r"hex\s*\([^)]*\)\s*\w*\s*\(\s*(\d+)\s+(\d+)\s+(\d+)\s*\)", txt)
    nx, ny, nz = (int(g) for g in m.groups())

    lo, hi = pts.min(axis=0), pts.max(axis=0)
    return nx, ny, nz, lo, hi


# --------------------------------------------------------------------------- #
# ASCII volScalarField / volVectorField internalField reader
# --------------------------------------------------------------------------- #
def read_scalar(path, ncells):
    txt = open(path, errors="replace").read()
    m = re.search(r"internalField\s+nonuniform[^(]*?(\d+)\s*\(", txt, flags=re.S)
    if m:
        start = txt.index("(", m.end() - 1) + 1
        end = txt.index(")", start)
        vals = np.fromstring(txt[start:end], sep=" ")
        if len(vals) != ncells:
            raise ValueError(f"{path}: got {len(vals)} values, mesh has {ncells}")
        return vals
    m = re.search(r"internalField\s+uniform\s+([-\d.eE+]+)\s*;", txt)
    if m:
        return np.full(ncells, float(m.group(1)))
    raise ValueError(f"{path}: cannot parse internalField")


def read_vector_mag(path, ncells):
    txt = open(path, errors="replace").read()
    m = re.search(r"internalField\s+nonuniform[^(]*?(\d+)\s*\(", txt, flags=re.S)
    if m:
        start = txt.index("(", m.end() - 1) + 1
        body = txt[start : txt.index("\n)", start) if "\n)" in txt[start:] else len(txt)]
        trip = re.findall(r"\(([^()]*)\)", body)
        arr = np.array([[float(v) for v in t.split()] for t in trip])
        return np.linalg.norm(arr, axis=1)
    m = re.search(r"internalField\s+uniform\s+\(([^)]*)\)", txt)
    if m:
        v = np.linalg.norm([float(x) for x in m.group(1).split()])
        return np.full(ncells, v)
    raise ValueError(f"{path}: cannot parse internalField")


def read_positions(path):
    """lagrangian/kinematicCloud/positions -> array of (x y z)."""
    txt = open(path, errors="replace").read()
    body = txt[txt.index("(", txt.index("FoamFile")) :]
    pts = []
    for t in re.findall(r"\(\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\)", body):
        pts.append([float(x) for x in t])
    return np.array(pts) if pts else np.zeros((0, 3))


# --------------------------------------------------------------------------- #
def main():
    case = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
    nx, ny, nz, lo, hi = read_mesh(case)
    ncells = nx * ny * nz
    dx = (hi[0] - lo[0]) / nx
    dy = (hi[1] - lo[1]) / ny
    dz = (hi[2] - lo[2]) / nz
    Vcell = dx * dy * dz

    idx = np.arange(ncells)
    yc = lo[1] + ((idx // nx) % ny + 0.5) * dy

    print(f"mesh: {nx} x {ny} x {nz} = {ncells} cells, "
          f"dx = {dx*1e6:.1f} um, dy = {dy*1e6:.1f} um, dz = {dz*1e6:.1f} um")

    times = []
    for d in os.listdir(case):
        if os.path.isdir(os.path.join(case, d)):
            try:
                times.append((float(d), d))
            except ValueError:
                pass
    times.sort()
    if not times:
        sys.exit("no time directories found - run the case first, "
                 "or reconstructPar has not produced full-domain data")

    print(f"\n{'time[s]':>9} {'V_air':>9} {'y_air':>8} {'V_oil':>9} {'y_oil':>8} "
          f"{'maxU':>8} {'n_p':>4} {'y_parc':>8}")
    print(f"{'':>9} {'[mm3]':>9} {'[mm]':>8} {'[mm3]':>9} {'[mm]':>8} "
          f"{'[m/s]':>8} {'':>4} {'[mm]':>8}")

    ts, ys = [], []
    for tval, tdir in times:
        d = os.path.join(case, tdir)
        try:
            a = read_scalar(os.path.join(d, "alpha.air"), ncells)
            o = read_scalar(os.path.join(d, "alpha.oil"), ncells)
        except (OSError, ValueError) as e:
            print(f"{tval:9g}  skipped ({e})")
            continue

        Va, Vo = a.sum() * Vcell, o.sum() * Vcell
        ya = (a * yc).sum() / max(a.sum(), 1e-30)
        yo = (o * yc).sum() / max(o.sum(), 1e-30)

        try:
            maxU = read_vector_mag(os.path.join(d, "U"), ncells).max()
        except (OSError, ValueError):
            maxU = float("nan")

        pfile = os.path.join(d, "lagrangian", "kinematicCloud", "positions")
        if os.path.isfile(pfile):
            P = read_positions(pfile)
            npar, ypar = len(P), (P[:, 1].mean() if len(P) else float("nan"))
        else:
            npar, ypar = 0, float("nan")

        print(f"{tval:9g} {Va*1e9:9.3f} {ya*1e3:8.4f} {Vo*1e9:9.3f} {yo*1e3:8.4f} "
              f"{maxU:8.4f} {npar:4d} {ypar*1e3:8.4f}")

        ts.append(tval)
        ys.append(ya)

    if len(ts) > 2:
        ts, ys = np.array(ts), np.array(ys)
        slope = np.polyfit(ts, ys, 1)[0]
        print(f"\nmean air-centroid rise velocity (linear fit): {slope*1e3:.2f} mm/s "
              f"= {slope:.4f} m/s")
        print("expected terminal velocity for a 0.7 mm air bubble in water:")
        print("   ~0.13 m/s  clean/mobile interface     ~0.079 m/s  contaminated/rigid")
        # volume drift
        try:
            a0 = read_scalar(os.path.join(case, times[0][1], "alpha.air"), ncells).sum() * Vcell
            aN = read_scalar(os.path.join(case, times[-1][1], "alpha.air"), ncells).sum() * Vcell
            o0 = read_scalar(os.path.join(case, times[0][1], "alpha.oil"), ncells).sum() * Vcell
            oN = read_scalar(os.path.join(case, times[-1][1], "alpha.oil"), ncells).sum() * Vcell
            print(f"air volume drift: {(aN-a0)/max(a0,1e-30)*100:+.2f} %   "
                  f"oil volume drift: {(oN-o0)/max(o0,1e-30)*100:+.2f} %")
            print("   both should be within a fraction of a percent.")
        except (OSError, ValueError):
            pass


if __name__ == "__main__":
    main()
