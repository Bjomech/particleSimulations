#!/usr/bin/env python3
"""
oilWaterAirParticle_postprocess.py
====================================
Enhanced visual dashboard for oilWaterAirParticle_case_final.
Reads OpenFOAM binary/ASCII fields and produces a polished PNG.
"""

import os, re, sys
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

# ----------------------------------------------------------------------
# Configuration
CASE = sys.argv[1] if len(sys.argv) > 1 else \
    os.path.expanduser("~/openfoam/cases/oilWaterAirParticle_case_final")
OUTPUT_PNG = os.path.join(CASE, "postprocess_dashboard.png")
OUTPUT_CSV = os.path.join(CASE, "postprocess_data.csv")

# Physical parameters (same as original)
rho_water = 1000.0;  rho_oil = 800.0;  rho_air = 1.2;  rho_p = 1500.0
mu_water  = 1e-3;    sigma_aw = 0.07;  sigma_ow = 0.03
d_p_nom   = 3e-4;    R_bubble = 0.05;  d_bubble = 0.1
g_acc     = 9.81
NX, NY    = 50, 150
DX, DY    = 0.5/NX, 1.0/NY

# ----------------------------------------------------------------------
# OpenFOAM binary / ASCII reader functions (unchanged from original)
# ----------------------------------------------------------------------
def read_foam_header(data_bytes):
    try:
        head = data_bytes[:2048].decode("latin-1")
    except Exception:
        return {}
    info = {}
    for key in ["format", "class", "arch", "object"]:
        m = re.search(rf'{key}\s+([^;]+);', head)
        if m:
            info[key] = m.group(1).strip().strip('"')
    arch = info.get("arch", "")
    lm = re.search(r'label=(\d+)', arch)
    sm = re.search(r'scalar=(\d+)', arch)
    info["label_bytes"]  = int(lm.group(1))//8 if lm else 4
    info["scalar_bytes"] = int(sm.group(1))//8 if sm else 8
    info["endian"]       = "<" if "LSB" in arch else ">"
    return info

def find_data_start(data_bytes):
    head = data_bytes[:4096].decode("latin-1", errors="replace")
    m = re.search(r'(\d+)\s*\n\s*\(', head)
    if not m:
        return None, None
    n         = int(m.group(1))
    paren_pos = data_bytes.find(b'(', m.start())
    return n, paren_pos + 1

def read_binary_scalar_field(filepath):
    try:
        with open(filepath, "rb") as f:
            data = f.read()
        info = read_foam_header(data)
        if info.get("format", "ascii") != "binary":
            return read_ascii_scalar_field(filepath)
        n, start = find_data_start(data)
        if n is None:
            return None
        sb     = info.get("scalar_bytes", 8)
        endian = info.get("endian", "<")
        dtype  = np.dtype(f"{endian}f{sb}")
        raw    = data[start:start + n*sb]
        if len(raw) < n*sb:
            return None
        return np.frombuffer(raw, dtype=dtype).copy().astype(np.float64)
    except Exception:
        return None

def read_ascii_scalar_field(filepath):
    try:
        with open(filepath) as f:
            content = f.read()
        u = re.search(r'internalField\s+uniform\s+([-\d.eE+]+)', content)
        if u:
            return np.array([float(u.group(1))])
        m = re.search(r'(\d+)\s*\n\s*\(', content)
        if not m:
            return None
        n     = int(m.group(1))
        start = content.find('(', m.start()) + 1
        end   = content.find('\n)', start)
        if end == -1:
            end = content.find(')', start)
        return np.array([float(v) for v in content[start:end].split()])
    except Exception:
        return None

def read_binary_vector_field(filepath):
    try:
        with open(filepath, "rb") as f:
            data = f.read()
        info = read_foam_header(data)
        if info.get("format", "ascii") != "binary":
            return read_ascii_vector_field(filepath)
        n, start = find_data_start(data)
        if n is None:
            return None
        sb     = info.get("scalar_bytes", 8)
        endian = info.get("endian", "<")
        dtype  = np.dtype(f"{endian}f{sb}")
        raw    = data[start:start + n*3*sb]
        if len(raw) < n*3*sb:
            return None
        arr = np.frombuffer(raw, dtype=dtype).copy().reshape(n, 3).astype(np.float64)
        return arr
    except Exception:
        return None

def read_ascii_vector_field(filepath):
    try:
        with open(filepath) as f:
            content = f.read()
        tuples = re.findall(
            r'\(\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\)',
            content)
        return np.array([[float(v) for v in t] for t in tuples]) if tuples else None
    except Exception:
        return None

def read_particle_positions(lag_dir):
    """
    Try multiple file names / formats used across OpenFOAM versions:
      1. lagrangian/kinematicCloud/positions  (classic, binary vector field)
      2. lagrangian/kinematicCloud/coordinates (OF >=2112, same layout)
      3. ASCII tuple scan on both files
    Returns array shape (N,3) or None.
    """
    for fname in ["positions", "coordinates"]:
        fp = os.path.join(lag_dir, fname)
        if not os.path.isfile(fp):
            continue
        # Try binary first
        arr = read_binary_vector_field(fp)
        if arr is not None and len(arr) > 0:
            # Sanity check: coordinates should be in a reasonable domain
            if np.all(np.abs(arr) < 1e6):
                return arr
        # Try ASCII tuple scan (handles both formats)
        arr2 = read_ascii_vector_field(fp)
        if arr2 is not None and len(arr2) > 0:
            return arr2
    return None

def bubble_centroid_y(alpha_vals, ny=NY, nx=NX):
    try:
        arr  = alpha_vals[:ny*nx].reshape((ny, nx))
        mask = arr > 0.5
        if not mask.any():
            return np.nan, 0
        rows = np.argwhere(mask)[:, 0]
        return (rows.mean() + 0.5) * DY, int(mask.sum())
    except Exception:
        return np.nan, 0

# ----------------------------------------------------------------------
# Collect data – with Co-only time-dir detection
# ----------------------------------------------------------------------
def get_times(case):
    dirs = []
    for d in os.listdir(case):
        try: dirs.append((float(d), d))
        except ValueError: pass
    return sorted(dirs)

def has_field_data(t_dir_path):
    """Return True if this time directory contains actual field files
    (not just Co / functionObject output)."""
    for fname in ["alpha.water", "alpha.air", "alpha.oil", "U", "p", "p_rgh"]:
        if os.path.isfile(os.path.join(t_dir_path, fname)):
            return True
    return False

print(f"Reading: {CASE}")
times = get_times(CASE)
print(f"  {len(times)} time directories found")

# Warn if time dirs contain only Co (functionObject) output
co_only = [d for _, d in times
           if os.path.isfile(os.path.join(CASE, d, "Co"))
           and not has_field_data(os.path.join(CASE, d))]
field_dirs = [(tv, td) for tv, td in times
              if has_field_data(os.path.join(CASE, td))]

if co_only:
    print(f"\n  ⚠  WARNING: {len(co_only)} time director(y/ies) contain only 'Co'")
    print("     (Courant-number functionObject output) — NO alpha.*, U, p found.")
    print("     These steps are SKIPPED for field post-processing.")
    print()
    print("  ► FIX in system/controlDict:")
    print("       writeControl    timeStep;   // or adjustableRunTime")
    print("       writeInterval   1;          // write fields every step")
    print("     The 'Co' file comes from a separate functions{} entry and")
    print("     writes independently of the main field writeInterval.")
    print("     Make sure 'functions > CourantNo > writeControl' does NOT")
    print("     override your main writeControl.\n")

if not field_dirs:
    print("  No time directories with field data found — using all dirs for Lagrangian only.")
    field_dirs = times   # fall back so particle data is still attempted

print(f"  {len(field_dirs)} time director(y/ies) with field data")

records = []
for t_val, t_dir in field_dirs:
    t_path = os.path.join(CASE, t_dir)
    lag    = os.path.join(t_path, "lagrangian", "kinematicCloud")
    rec    = {"time": t_val}

    # Particle position (robust multi-format reader)
    pos_arr = read_particle_positions(lag)
    if pos_arr is not None and len(pos_arr) > 0:
        rec["px"], rec["py"], rec["pz"] = pos_arr[0]
    else:
        rec["px"] = rec["py"] = rec["pz"] = np.nan

    # Particle velocity
    vel_arr = read_binary_vector_field(os.path.join(lag, "U"))
    if vel_arr is None:
        vel_arr = read_ascii_vector_field(os.path.join(lag, "U"))
    if vel_arr is not None and len(vel_arr) > 0:
        rec["pvx"], rec["pvy"], rec["pvz"] = vel_arr[0]
        rec["pv_mag"] = float(np.linalg.norm(vel_arr[0]))
    else:
        rec["pvx"] = rec["pvy"] = rec["pvz"] = rec["pv_mag"] = np.nan

    # Particle diameter
    d_arr = read_binary_scalar_field(os.path.join(lag, "d"))
    if d_arr is None:
        d_arr = read_ascii_scalar_field(os.path.join(lag, "d"))
    rec["d_p"] = float(d_arr[0]) if d_arr is not None and len(d_arr) > 0 else d_p_nom

    # Particle age
    age_arr = read_binary_scalar_field(os.path.join(lag, "age"))
    rec["age"] = float(age_arr[0]) if age_arr is not None and len(age_arr) > 0 else np.nan

    # Alpha fields
    for phase in ["air", "oil", "water"]:
        fpath = os.path.join(t_path, f"alpha.{phase}")
        arr   = read_binary_scalar_field(fpath)
        if arr is None:
            arr = read_ascii_scalar_field(fpath)
        if arr is not None and len(arr) >= NX*NY:
            rec[f"alpha_{phase}_mean"] = float(arr[:NX*NY].mean())
            rec[f"alpha_{phase}_max"]  = float(arr[:NX*NY].max())
            if phase == "air":
                cy, nc = bubble_centroid_y(arr)
                rec["bubble_cy"]    = cy
                rec["bubble_cells"] = nc
        else:
            rec[f"alpha_{phase}_mean"] = np.nan
            rec[f"alpha_{phase}_max"]  = np.nan
            if phase == "air":
                rec["bubble_cy"]    = np.nan
                rec["bubble_cells"] = 0

    records.append(rec)

df = pd.DataFrame(records)
df = df[df["time"] > 0].reset_index(drop=True)

# Clamp physically impossible particle positions
df.loc[(df["px"].abs() > 10) | (df["py"].abs() > 10), ["px","py","pz"]] = np.nan
df.loc[df["pv_mag"] > 100, "pv_mag"] = np.nan

if "bubble_cy" not in df.columns:
    df["bubble_cy"] = np.nan
df["bubble_cy"] = df["bubble_cy"].ffill()
df["bubble_vy"] = np.gradient(df["bubble_cy"].ffill().values,
                               df["time"].values)
df["v_rel"]     = (df["pvy"] - df["bubble_vy"]).abs()
df["in_bubble"] = (np.abs(df["py"] - df["bubble_cy"]) < R_bubble).astype(int)

df.to_csv(OUTPUT_CSV, index=False)
print(f"  CSV → {OUTPUT_CSV}")
print(df[["time","px","py","pv_mag","bubble_cy","bubble_cells"]].head(8).to_string(index=False))

# ----------------------------------------------------------------------
# Flotation numbers
# ----------------------------------------------------------------------
def flotation_numbers(row):
    dp    = row["d_p"]    if pd.notna(row["d_p"])    else d_p_nom
    U_p   = row["pv_mag"] if pd.notna(row["pv_mag"]) else 0.0
    U_b   = abs(row["bubble_vy"]) if pd.notna(row["bubble_vy"]) else 1e-4
    U_rel = max(row["v_rel"] if pd.notna(row["v_rel"]) else abs(U_p - U_b), 1e-12)

    Re_p   = rho_water * U_rel * dp / mu_water
    St     = rho_p * dp**2 * U_b / (9 * mu_water * R_bubble)
    We_p   = rho_water * U_rel**2 * dp / sigma_aw
    Bo     = (rho_water - rho_air) * g_acc * d_bubble**2 / sigma_aw
    Ar     = rho_water * abs(rho_water - rho_air) * g_acc * d_bubble**3 / mu_water**2
    Ca     = mu_water * U_rel / sigma_aw
    Fr     = U_b / max(np.sqrt(g_acc * R_bubble), 1e-12)
    P_coll = St / (St + 0.25)
    t_star = (row["age"] * U_rel / dp) if pd.notna(row["age"]) else np.nan
    return pd.Series({"Re_p":Re_p, "St":St, "We_p":We_p, "Bo":Bo,
                      "Ar":Ar, "Ca":Ca, "Fr":Fr, "P_coll":P_coll, "t_star":t_star})

df = pd.concat([df, df.apply(flotation_numbers, axis=1)], axis=1)

# ----------------------------------------------------------------------
# ENHANCED DASHBOARD with improved styling
# ----------------------------------------------------------------------
pio.templates.default = "plotly_white"

# Modern color palette
colors = {
    'air': '#4A90E2',       # bright blue
    'oil': '#F5A623',       # warm orange
    'water': '#50E3C2',     # teal
    'particle': '#D0021B',  # strong red
    'bubble': '#8B9DC3',    # muted periwinkle
    'Re': '#9013FE',        # purple
    'St': '#417505',        # dark green
    'We': '#F8E71C',        # bright yellow
    'Ca': '#9B9B9B',        # grey
    'Fr': '#E94F6F',        # pinkish red
    'P': '#4A4A4A',         # dark grey
    'Bo': '#B8E986',        # light green
    'Ar': '#FFB347'         # light orange
}

# Subplot layout: 6 rows × 2 cols, top row for parameters
ROWS, COLS   = 6, 2
row_heights  = [0.05] + [0.19]*5   # first row (param) small, then 5 equal
vertical_spacing = 0.07
horizontal_spacing = 0.12

fig = make_subplots(
    rows=ROWS, cols=COLS,
    specs=[[{"colspan": 2}, None]] + [[{}, {}] for _ in range(5)],
    subplot_titles=[
        "",  # param row – title added manually
        "① Particle trajectory (x–y)",
        "② Particle y(t) vs. bubble centroid",
        "③ Particle speed |Uₚ|(t)",
        "④ Bubble rise height cy(t)",
        "⑤ Phase volume fractions α(t)",
        "⑥ Relative velocity |ΔU|(t)  (shaded = inside bubble)",
        "⑦ Particle Reynolds number Reₚ(t)",
        "⑧ Stokes, Weber, Capillary, Froude numbers",
        "⑨ Sutherland collision probability P_coll(t)",
        "⑩ Bond & Archimedes numbers (bubble, constant)"
    ],
    row_heights=row_heights,
    vertical_spacing=vertical_spacing,
    horizontal_spacing=horizontal_spacing
)

# Style subplot titles
for ann in fig.layout.annotations:
    if ann.text and ann.text[0] in "①②③④⑤⑥⑦⑧⑨⑩":
        ann.font = dict(size=15, family="Arial Black", color="#1F2A3A")
        ann.bgcolor = "rgba(230, 242, 255, 0.8)"
        ann.borderpad = 4
        ann.borderwidth = 1
        ann.bordercolor = "#A3C6FF"

# Hide axes for parameter row
fig.update_xaxes(visible=False, row=1, col=1)
fig.update_yaxes(visible=False, row=1, col=1)

# Parameter panel (styled as a bordered box with line breaks)
param_text = (
    "<b>Simulation parameters</b><br>"
    f"ρ<sub>water</sub>={rho_water:.0f} kg/m³ ρ<sub>oil</sub>={rho_oil:.0f} kg/m³ "
    f"ρ<sub>air</sub>={rho_air:.1f} kg/m³ ρ<sub>p</sub>={rho_p:.0f} kg/m³ "
    f"μ<sub>water</sub>={mu_water:.0e} Pa·s<br>"
    f"σ<sub>aw</sub>={sigma_aw:.3f} N/m σ<sub>ow</sub>={sigma_ow:.3f} N/m "
    f"d<sub>p</sub>={d_p_nom*1e3:.2f} mm R<sub>bubble</sub>={R_bubble*100:.0f} cm "
    f"d<sub>bubble</sub>={d_bubble*100:.0f} cm g={g_acc:.2f} m/s²<br>"
    f"Mesh {NX}×{NY} Δx={DX*1e3:.1f} mm Δy={DY*1e3:.2f} mm"
)
fig.add_annotation(
    text=param_text,
    x=0.5, y=1.0 - row_heights[0]/2,
    xref="paper", yref="paper",
    showarrow=False,
    font=dict(size=12, family="monospace", color="#1F3A4B"),
    align="center", xanchor="center", yanchor="middle",
    bgcolor="rgba(240, 248, 255, 0.98)",
    bordercolor="#7FA8C9",
    borderwidth=1.5,
    borderpad=8,
    opacity=0.95
)

# Helper for axis styling
def style_axes(fig, row, col, x_title="", y_title=""):
    if x_title:
        fig.update_xaxes(title_text=x_title, title_font=dict(size=12, family="Arial"),
                         gridcolor="lightgrey", linecolor="black", mirror=True,
                         row=row, col=col)
    if y_title:
        fig.update_yaxes(title_text=y_title, title_font=dict(size=12, family="Arial"),
                         gridcolor="lightgrey", linecolor="black", mirror=True,
                         row=row, col=col)
    fig.update_xaxes(minor=dict(gridcolor="whitesmoke"), row=row, col=col)
    fig.update_yaxes(minor=dict(gridcolor="whitesmoke"), row=row, col=col)

t = df["time"]

# 1) Particle trajectory (2,1)
fig.add_trace(go.Scatter(
    x=df["px"], y=df["py"], mode="lines+markers",
    marker=dict(color=t, colorscale="Plasma", size=5,
                colorbar=dict(title="Time [s]", x=-0.12, len=0.35)),
    line=dict(color="rgba(100,100,100,0.5)", width=1),
    name="trajectory",
    hovertemplate="x=%{x:.4f} m<br>y=%{y:.4f} m"
), row=2, col=1)
style_axes(fig, 2, 1, x_title="x [m]", y_title="y [m]")

# 2) py vs bubble centroid (2,2)
fig.add_trace(go.Scatter(
    x=t, y=df["py"], mode="lines",
    line=dict(color=colors["particle"], width=2.5), name="Particle y",
    fill="tozeroy", fillcolor="rgba(208,2,27,0.1)"
), row=2, col=2)
fig.add_trace(go.Scatter(
    x=t, y=df["bubble_cy"], mode="lines",
    line=dict(color=colors["air"], width=2, dash="dash"),
    name="Bubble centroid y"
), row=2, col=2)
style_axes(fig, 2, 2, x_title="t [s]", y_title="y [m]")

# 3) |U_p| (3,1)
fig.add_trace(go.Scatter(
    x=t, y=df["pv_mag"], mode="lines",
    line=dict(color=colors["particle"], width=2.5),
    fill="tozeroy", fillcolor="rgba(208,2,27,0.05)",
    name="|Uₚ|"
), row=3, col=1)
style_axes(fig, 3, 1, x_title="t [s]", y_title="|Uₚ| [m/s]")

# 4) Bubble cy (3,2)
fig.add_trace(go.Scatter(
    x=t, y=df["bubble_cy"], mode="lines",
    line=dict(color=colors["air"], width=2.5),
    fill="tozeroy", fillcolor="rgba(74,144,226,0.15)",
    name="Bubble centroid y"
), row=3, col=2)
style_axes(fig, 3, 2, x_title="t [s]", y_title="y [m]")
fig.update_yaxes(range=[0, 1.05], row=3, col=2)

# 5) Phase fractions (4,1)
for phase, col, label in [("air", colors["air"], "α_air"),
                           ("oil", colors["oil"], "α_oil"),
                           ("water", colors["water"], "α_water")]:
    fig.add_trace(go.Scatter(
        x=t, y=df[f"alpha_{phase}_mean"], mode="lines",
        line=dict(color=col, width=2),
        fill="tozeroy", fillcolor=f"rgba{tuple(int(col[i:i+2],16) for i in (1,3,5)) + (0.08,)}",
        name=label
    ), row=4, col=1)
style_axes(fig, 4, 1, x_title="t [s]", y_title="α [-]")
fig.add_hline(y=0.5, line_dash="dot", line_color="darkgrey", row=4, col=1)

# 6) Relative velocity + bubble shading (4,2)
fig.add_trace(go.Scatter(
    x=t, y=df["v_rel"], mode="lines",
    line=dict(color=colors["Fr"], width=2.5),
    fill="tozeroy", fillcolor="rgba(233,79,111,0.1)",
    name="|ΔU|"
), row=4, col=2)
# Shade when particle inside bubble
in_bub = df["in_bubble"].values
for i in range(1, len(t)):
    if in_bub[i-1] == 1 and in_bub[i] == 1:
        fig.add_vrect(x0=float(t.iloc[i-1]), x1=float(t.iloc[i]),
                      fillcolor="rgba(74,144,226,0.2)", line_width=0,
                      row=4, col=2)
style_axes(fig, 4, 2, x_title="t [s]", y_title="|ΔU| [m/s]")

# 7) Re_p (5,1)
fig.add_trace(go.Scatter(
    x=t, y=df["Re_p"], mode="lines",
    line=dict(color=colors["Re"], width=2.5),
    fill="tozeroy", fillcolor="rgba(144,19,254,0.08)",
    name="Reₚ"
), row=5, col=1)
style_axes(fig, 5, 1, x_title="t [s]", y_title="Reₚ [-]")
fig.add_hline(y=1, line_dash="dot", line_color="grey", annotation_text="Stokes (Re=1)",
              annotation_position="top left", row=5, col=1)
fig.add_hline(y=1000, line_dash="dot", line_color="grey", annotation_text="Inertial (Re=1000)",
              annotation_position="top left", row=5, col=1)

# 8) St, We, Ca, Fr (log) (5,2)
for key, col, label in [("St", colors["St"], "St"), ("We_p", colors["We"], "We"),
                         ("Ca", colors["Ca"], "Ca"), ("Fr", colors["Fr"], "Fr")]:
    yvals = df[key].clip(lower=1e-12)
    fig.add_trace(go.Scatter(
        x=t, y=yvals, mode="lines",
        line=dict(color=col, width=2),
        name=label
    ), row=5, col=2)
style_axes(fig, 5, 2, x_title="t [s]", y_title="[-]")
fig.update_yaxes(type="log", title_text="[-]", row=5, col=2)
fig.add_hline(y=1, line_dash="dot", line_color="grey", annotation_text="= 1",
              row=5, col=2)

# 9) P_coll (6,1)
fig.add_trace(go.Scatter(
    x=t, y=df["P_coll"], mode="lines",
    line=dict(color=colors["P"], width=2.5),
    fill="tozeroy", fillcolor="rgba(74,74,74,0.1)",
    name="P_coll"
), row=6, col=1)
style_axes(fig, 6, 1, x_title="t [s]", y_title="P_coll [-]")
fig.add_hline(y=0.5, line_dash="dot", line_color="grey", annotation_text="P=0.5",
              row=6, col=1)

# 10) Bo & Ar (constant) (6,2)
Bo_val = float(df["Bo"].mean())
Ar_val = float(df["Ar"].mean())
fig.add_trace(go.Scatter(
    x=t, y=np.full(len(t), Bo_val), mode="lines",
    line=dict(color=colors["Bo"], width=2.5, dash="dash"),
    name=f"Bo = {Bo_val:.1f}"
), row=6, col=2)
fig.add_trace(go.Scatter(
    x=t, y=np.full(len(t), Ar_val), mode="lines",
    line=dict(color=colors["Ar"], width=2.5, dash="dot"),
    name=f"Ar = {Ar_val:.3g}"
), row=6, col=2)
style_axes(fig, 6, 2, x_title="t [s]", y_title="[-]")
fig.update_yaxes(type="log", row=6, col=2)

# ------------------------------------------------------------------
# Add formulae as small footnotes (bottom-left of each subplot)
# ------------------------------------------------------------------
formulas = {
    (2,1): "x(t), y(t) from Lagrangian positions",
    (2,2): "bubble cy = mean row(α_air > 0.5)·Δy",
    (3,1): "|Uₚ| = √(Ux²+Uy²+Uz²)",
    (3,2): "cy(t) = Σ(y_i·α_air,i)/Σα_air,i",
    (4,1): "α = domain mean of volScalarField",
    (4,2): "|ΔU| = |Uₚ,y − dcy/dt|; blue = inside bubble",
    (5,1): "Reₚ = ρ_w·|ΔU|·dₚ/μ_w",
    (5,2): "St = ρₚ dₚ² U_b/(9μ_w R_b) We = ρ_w|ΔU|² dₚ/σ_aw Ca = μ_w|ΔU|/σ_aw Fr = U_b/√(gR_b)",
    (6,1): "P_coll = St / (St + 0.25)  [Sutherland 1948]",
    (6,2): "Bo = (ρ_w−ρ_a)g d_b²/σ_aw Ar = ρ_w|Δρ| g d_b³/μ_w²"
}

# Compute y positions for the bottom of each subplot (paper coordinates)
gap = vertical_spacing
cum_top = {}
y_cursor = 1.0
for ri, rh in enumerate(row_heights):
    cum_top[ri+1] = y_cursor
    y_cursor -= rh
    y_cursor -= gap

for (r, c), txt in formulas.items():
    y_pos = cum_top[r] - row_heights[r-1] + 0.01  # 0.01 above bottom edge
    x_pos = 0.02 if c == 1 else 0.52  # left side of each column
    fig.add_annotation(
        text=f"<i>{txt}</i>",
        x=x_pos, y=y_pos,
        xref="paper", yref="paper",
        showarrow=False,
        font=dict(size=9, color="#6B6B6B"),
        align="left", xanchor="left", yanchor="bottom"
    )

# ------------------------------------------------------------------
# Final layout adjustments
# ------------------------------------------------------------------
fig.update_layout(
    title=dict(
        text="<b>Oil‑Water‑Air Flotation & Particle‑Bubble Interaction</b><br>"
             "<sup>Reₚ · St · We · Ca · Fr · P_coll (Sutherland) | Blue shading = particle inside bubble</sup>",
        font=dict(size=18, family="Arial Black"),
        x=0.5,
        xanchor="center"
    ),
    width=1900,
    height=2300,
    legend=dict(
        orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5,
        font=dict(size=10), bgcolor="rgba(255,255,255,0.8)", bordercolor="lightgrey", borderwidth=1
    ),
    paper_bgcolor="white",
    plot_bgcolor="white",
    margin=dict(l=60, r=40, t=120, b=40)
)

# Save high-resolution PNG
fig.write_image(OUTPUT_PNG, scale=2.5)
print(f"Dashboard saved → {OUTPUT_PNG}")

# ------------------------------------------------------------------
# Console summary
# ------------------------------------------------------------------
print("\n" + "="*64)
print("TIME-AVERAGED FLOTATION NUMBERS")
print("="*64)
for c in ["Re_p", "St", "We_p", "Bo", "Ar", "Ca", "Fr", "P_coll"]:
    print(f"  {c:<10s} = {df[c].mean():.4g}")

st = df["St"].mean(); re = df["Re_p"].mean()
print("\nRegime interpretation:")
print(f"  Reₚ = {re:.3g} → " + (
    "Stokes (Re<1)" if re < 1 else
    "Transitional (1<Re<1000)" if re < 1000 else "Inertial (Re>1000)"))
print(f"  St  = {st:.3g} → " + (
    "Follows streamlines (St<<1)" if st < 0.1 else
    "Intermediate inertia" if st < 1 else "High inertia → efficient collision"))
