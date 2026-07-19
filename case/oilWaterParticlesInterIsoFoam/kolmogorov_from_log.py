#!/usr/bin/env python3
"""
Kolmogorov & Batchelor scale checker for OpenFOAM logs.

- Parses time step info from an OpenFOAM log:
    deltaT = 0.00012
    Time = 0.00012

- Computes Kolmogorov microscales (eta, tau_eta, u_eta)
- Computes Batchelor scalar scale (eta_B) for given Schmidt number
- Reports resolution ratios: dx/eta, dx/eta_B, dt/tau_eta
"""

import re
import sys
import numpy as np
from pathlib import Path

# ============================================================
# USER INPUTS — EDIT THESE FOR YOUR CASE
# ============================================================
logfile       = "log.interisofoamSolver"  # path to OpenFOAM log
nu            = 1.0e-6    # kinematic viscosity [m^2/s]
Sc            = 1.0       # Schmidt number (nu / D); set >1 to get Batchelor scale
dx            = 1.0e-3    # representative cell size [m] (use min dx if you want)
epsilon_const = 1.0e-3    # dissipation rate [m^2/s^3]; set from U^3/L or field average

# ============================================================
# PARSE deltaT / Time FROM LOG
# ============================================================

# Examples in your log:
#   Time   : 10:29:02
#   deltaT = 0.00012
#   Time = 0.00012
#
# We want to:
#   - ignore "Time   : 10:29:02"
#   - pair "deltaT = ..." with the next "Time = ..."

DT_RE   = re.compile(r"deltaT\s*=\s*([0-9.eE+\-]+)")
TIME_RE = re.compile(r"\bTime\s*=\s*([0-9.eE+\-]+)")  # only 'Time =', not 'Time   :'

times = []
dts   = []
pending_dt = None  # holds the last seen deltaT until we see Time = ...

try:
    with open(logfile, "r", errors="replace") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue

            # match deltaT line
            m_dt = DT_RE.search(line)
            if m_dt:
                pending_dt = float(m_dt.group(1))
                continue

            # match simulation time line
            m_t = TIME_RE.search(line)
            if m_t and pending_dt is not None:
                t = float(m_t.group(1))
                times.append(t)
                dts.append(pending_dt)
                pending_dt = None
                continue

except FileNotFoundError:
    sys.exit(f"ERROR: Log file '{logfile}' not found in {Path.cwd()}")

if not times:
    raise RuntimeError("No deltaT/Time pairs found. Check logfile name or regex patterns.")

times = np.array(times, dtype=float)
dts   = np.array(dts,   dtype=float)

print(
    f"Parsed {len(times)} time steps from {logfile}  |  "
    f"t_start = {times[0]:.4g} s,  t_end = {times[-1]:.4g} s,  "
    f"last dt = {dts[-1]:.3e} s"
)

# ============================================================
# KOLMOGOROV & BATCHELOR SCALES
# ============================================================

# Use constant epsilon for now (same for all times)
eps = np.full_like(times, epsilon_const, dtype=float)

# Kolmogorov microscales
eta     = (nu**3 / eps)**0.25      # length scale [m]
tau_eta = (nu / eps)**0.5          # time scale [s]
u_eta   = (nu * eps)**0.25         # velocity scale [m/s]

# Batchelor (scalar) scale for given Sc
eta_B = eta / np.sqrt(Sc)          # [m]

# Resolution ratios
dx_over_eta    = dx / eta
dx_over_etaB   = dx / eta_B
dt_over_tauEta = dts / tau_eta

# ============================================================
# SUMMARY AT LAST TIME STEP
# ============================================================

i = -1  # index for last time step
sep = "-" * 60

print(f"\n{sep}")
print("Kolmogorov / Batchelor Scale Check")
print(sep)
print(f"nu             = {nu:.3e}  m^2/s")
print(f"Sc             = {Sc:.3g}")
print(f"dx             = {dx:.3e}  m")
print(f"epsilon (input)= {epsilon_const:.3e}  m^2/s^3")
print(sep)
print(f"At t = {times[i]:.6g} s:")
print(f"  eta          = {eta[i]:.3e}  m      (Kolmogorov length)")
print(f"  tau_eta      = {tau_eta[i]:.3e}  s      (Kolmogorov time)")
print(f"  u_eta        = {u_eta[i]:.3e}  m/s    (Kolmogorov velocity)")
print(f"  eta_B        = {eta_B[i]:.3e}  m      (Batchelor scale)")
print(sep)
print("Resolution ratios (target: < 1 for fully resolved):")
print(f"  dx / eta     = {dx_over_eta[i]:.3f}")
print(f"  dx / eta_B   = {dx_over_etaB[i]:.3f}")
print(f"  dt / tau_eta = {dt_over_tauEta[i]:.3f}")
print(sep)

# ============================================================
# SAVE FULL TIME HISTORY TO CSV (OPTIONAL)
# ============================================================

out = Path("kolmogorov_scales_from_log.csv")
header = "time,dt,epsilon,eta,tau_eta,u_eta,eta_B,dx_over_eta,dx_over_etaB,dt_over_tauEta"
data = np.column_stack(
    (times, dts, eps, eta, tau_eta, u_eta, eta_B, dx_over_eta, dx_over_etaB, dt_over_tauEta)
)
np.savetxt(out, data, delimiter=",", header=header, comments="")
print(f"\nSaved full time history to: {out}")
