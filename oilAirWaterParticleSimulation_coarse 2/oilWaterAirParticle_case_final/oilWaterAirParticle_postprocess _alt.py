#!/usr/bin/env python3
"""
oilWaterAirParticle_postprocess.py
====================================
Post-processing PNG dashboard for oilWaterAirParticle_case_final.
Handles OpenFOAM binary format (format binary; arch LSB;label=32;scalar=64)

Panels:
  1  Particle Trajectory (x-y)
  2  Particle py(t) vs Bubble Centroid
  3  Particle |U_p|(t)
  4  Bubble Rise cy(t)
  5  Phase Volume Fractions
  6  Particle-Bubble |ΔU|(t)
  7  Re_p(t)                        ← dedicated Re panel
  8  St / We_p / Ca / Fr (log)
  9  P_coll(t)  Sutherland
  10 Bo & Ar (constant reference)

Usage:  python3 oilWaterAirParticle_postprocess.py [CASE_PATH]
"""

import os, re, sys, struct
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

CASE = sys.argv[1] if len(sys.argv) > 1 else \
       os.path.expanduser("~/openfoam/cases/oilWaterAirParticle_case_final")
OUTPUT_PNG = os.path.join(CASE, "postprocess_dashboard.png")
OUTPUT_CSV = os.path.join(CASE, "postprocess_data.csv")

# ── Physical parameters ──────────────────────────────────────────
rho_water = 1000.0; rho_oil = 800.0; rho_air = 1.2; rho_p = 1500.0
mu_water  = 1e-3;   sigma_aw = 0.07; sigma_ow = 0.03
d_p_nom   = 3e-4;   R_bubble = 0.05; d_bubble = 0.1
g_acc     = 9.81
NX, NY    = 50, 150
DX, DY    = 0.5/NX, 1.0/NY


# ═══════════════════════════════════════════════════════════════════
# BINARY / ASCII OPENFOAM READER
# ═══════════════════════════════════════════════════════════════════

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
    info["endian"] = "<" if "LSB" in arch else ">"
    return info


def find_data_start(data_bytes):
    head = data_bytes[:4096].decode("latin-1", errors="replace")
    m = re.search(r'(\d+)\s*\n\s*\(', head)
    if not m:
        return None, None
    n = int(m.group(1))
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
        sb = info.get("scalar_bytes", 8)
        endian = info.get("endian", "<")
        dtype = np.dtype(f"{endian}f{sb}")
        raw = data[start:start + n*sb]
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
        n = int(m.group(1))
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
        sb = info.get("scalar_bytes", 8)
        endian = info.get("endian", "<")
        dtype = np.dtype(f"{endian}f{sb}")
        raw = data[start:start + n*3*sb]
        if len(raw) < n*3*sb:
            return None
        return np.frombuffer(raw, dtype=dtype).copy().reshape(n, 3).astype(np.float64)
    except Exception:
        return None


def read_ascii_vector_field(filepath):
    try:
        with open(filepath) as f:
            content = f.read()
        m = re.search(r'(\d+)\s*\n\s*\(', content)
        if not m:
            return None
        start = content.find('(', m.start()) + 1
        end   = content.rfind(')')
        tuples = re.findall(
            r'\(\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\)',
            content[start:end])
        return np.array([[float(v) for v in t] for t in tuples]) if tuples else None
    except Exception:
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


# ═══════════════════════════════════════════════════════════════════
# COLLECT DATA
# ═══════════════════════════════════════════════════════════════════

def get_times(case):
    dirs = []
    for d in os.listdir(case):
        try: dirs.append((float(d), d))
        except ValueError: pass
    return sorted(dirs)

print(f"Reading: {CASE}")
times = get_times(CASE)
print(f"  {len(times)} time directories")

records = []
for t_val, t_dir in times:
    lag = os.path.join(CASE, t_dir, "lagrangian", "kinematicCloud")
    rec = {"time": t_val}

    # Particle position
    pos_arr = read_binary_vector_field(os.path.join(lag, "positions"))
    if pos_arr is not None and len(pos_arr) > 0:
        rec["px"], rec["py"], rec["pz"] = pos_arr[0]
    else:
        rec["px"] = rec["py"] = rec["pz"] = np.nan

    # Particle velocity
    vel_arr = read_binary_vector_field(os.path.join(lag, "U"))
    if vel_arr is not None and len(vel_arr) > 0:
        rec["pvx"], rec["pvy"], rec["pvz"] = vel_arr[0]
        rec["pv_mag"] = float(np.linalg.norm(vel_arr[0]))
    else:
        rec["pvx"] = rec["pvy"] = rec["pvz"] = rec["pv_mag"] = np.nan

    # Particle diameter
    d_arr = read_binary_scalar_field(os.path.join(lag, "d"))
    rec["d_p"] = float(d_arr[0]) if d_arr is not None and len(d_arr) > 0 else d_p_nom

    # Particle age
    age_arr = read_binary_scalar_field(os.path.join(lag, "age"))
    rec["age"] = float(age_arr[0]) if age_arr is not None and len(age_arr) > 0 else np.nan

    # Alpha fields
    for phase in ["air", "oil", "water"]:
        fpath = os.path.join(CASE, t_dir, f"alpha.{phase}")
        arr = read_binary_scalar_field(fpath)
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

# Clamp garbage particle positions from binary parse artifacts
df.loc[(df["px"].abs() > 10) | (df["py"].abs() > 10), ["px","py","pz"]] = np.nan
df.loc[df["pv_mag"] > 100, "pv_mag"] = np.nan

if "bubble_cy" not in df.columns:
    df["bubble_cy"] = np.nan

df["bubble_cy"] = df["bubble_cy"].ffill()
df["bubble_vy"] = np.gradient(df["bubble_cy"].values, df["time"].values)
df["v_rel"]     = (df["pvy"] - df["bubble_vy"]).abs()
df["in_bubble"] = (np.abs(df["py"] - df["bubble_cy"]) < R_bubble).astype(int)

df.to_csv(OUTPUT_CSV, index=False)
print(f"  CSV → {OUTPUT_CSV}")
print(df[["time","px","py","pv_mag","bubble_cy","bubble_cells"]].head(8).to_string(index=False))


# ═══════════════════════════════════════════════════════════════════
# FLOTATION NUMBERS
# ═══════════════════════════════════════════════════════════════════

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
    return pd.Series({"Re_p":Re_p,"St":St,"We_p":We_p,"Bo":Bo,
                      "Ar":Ar,"Ca":Ca,"Fr":Fr,"P_coll":P_coll,"t_star":t_star})

df = pd.concat([df, df.apply(flotation_numbers, axis=1)], axis=1)


# ═══════════════════════════════════════════════════════════════════
# FORMULA ANNOTATIONS
# each entry: (panel_row, panel_col, formula_string)
# ═══════════════════════════════════════════════════════════════════

FORMULAS = {
    # (row, col): formula text shown below the subplot
    (1,1): "trajectory: (x(t), y(t)) from Lagrangian positions field",
    (1,2): "py(t) from Lagrangian positions  |  bubble cy = mean row of α_air > 0.5",
    (2,1): "|U_p| = √(Ux² + Uy² + Uz²)  from Lagrangian U field",
    (2,2): "cy(t) = Σ(row_i · DY) / N_cells,  cells where α_air > 0.5",
    (3,1): "α_phase = domain mean of internalField  (volScalarField)",
    (3,2): "|ΔU| = |U_p,y − dcy/dt|  |  blue = particle inside bubble",
    (4,1): "Re_p = ρ_w · |ΔU| · d_p / μ_w",
    (4,2): "St = ρ_p·d_p²·U_b / (9·μ_w·R_b)   We = ρ_w·|ΔU|²·d_p / σ_aw\n"
           "Ca = μ_w·|ΔU| / σ_aw               Fr = U_b / √(g·R_b)",
    (5,1): "P_coll = St / (St + 0.25)  [Sutherland 1948 flotation collision model]",
    (5,2): "Bo = (ρ_w−ρ_air)·g·d_b² / σ_aw     Ar = ρ_w·|Δρ|·g·d_b³ / μ_w²",
}


# ═══════════════════════════════════════════════════════════════════
# DASHBOARD  (5 rows × 2 cols)
# ═══════════════════════════════════════════════════════════════════

pio.templates.default = "plotly_white"

C = dict(
    air="#64B5F6", oil="#FFD54F", water="#4DB6AC",
    particle="#EF5350", bubble="#90CAF9",
    Re="#AB47BC", St="#26A69A", We="#FFA726",
    Ca="#78909C", Fr="#EC407A", P="#66BB6A",
    Bo="#5C6BC0", Ar="#FF7043",
)

ROWS, COLS = 5, 2
fig = make_subplots(
    rows=ROWS, cols=COLS,
    subplot_titles=[
        "① Particle Trajectory  (x–y, colour = time)",
        "② Particle py(t)  vs  Bubble Centroid",
        "③ Particle  |U_p|(t)",
        "④ Bubble Rise  cy(t)",
        "⑤ Phase Volume Fractions  α(t)",
        "⑥ Particle–Bubble  |ΔU|(t)",
        "⑦ Particle Reynolds Number  Re_p(t)",
        "⑧ St · We_p · Ca · Fr  [log scale]",
        "⑨ Sutherland Collision Probability  P_coll(t)",
        "⑩ Bond (Bo) & Archimedes (Ar)  [bubble, constant]",
    ],
    vertical_spacing=0.07,
    horizontal_spacing=0.10,
)

t = df["time"]

# ① Trajectory
fig.add_trace(go.Scatter(
    x=df["px"], y=df["py"], mode="lines+markers",
    marker=dict(color=t, colorscale="Viridis", size=5,
                colorbar=dict(title="t [s]", x=-0.07, len=0.36, y=0.90)),
    line=dict(color="rgba(150,150,150,0.4)", width=1),
    name="Trajectory",
    hovertemplate="x=%{x:.4f} m<br>y=%{y:.4f} m",
), row=1, col=1)

# ② py + bubble centroid
fig.add_trace(go.Scatter(
    x=t, y=df["py"], mode="lines",
    line=dict(color=C["particle"], width=2.5), name="Particle py",
    fill="tozeroy", fillcolor="rgba(239,83,80,0.08)",
), row=1, col=2)
fig.add_trace(go.Scatter(
    x=t, y=df["bubble_cy"], mode="lines",
    line=dict(color=C["air"], width=2, dash="dash"),
    name="Bubble cy",
), row=1, col=2)

# ③ |U_p|
fig.add_trace(go.Scatter(
    x=t, y=df["pv_mag"], mode="lines",
    line=dict(color=C["particle"], width=2), name="|U_p|",
), row=2, col=1)

# ④ Bubble cy
fig.add_trace(go.Scatter(
    x=t, y=df["bubble_cy"], mode="lines",
    line=dict(color=C["air"], width=2),
    fill="tozeroy", fillcolor="rgba(100,181,246,0.12)",
    name="Bubble cy",
), row=2, col=2)

# ⑤ Phase fractions
for phase, col, label in [("air",C["air"],"α_air"),
                           ("oil",C["oil"],"α_oil"),
                           ("water",C["water"],"α_water")]:
    fig.add_trace(go.Scatter(
        x=t, y=df[f"alpha_{phase}_mean"], mode="lines",
        line=dict(color=col, width=2), name=label,
    ), row=3, col=1)

# ⑥ Relative velocity + in-bubble shading
fig.add_trace(go.Scatter(
    x=t, y=df["v_rel"], mode="lines",
    line=dict(color=C["Fr"], width=2),
    fill="tozeroy", fillcolor="rgba(236,64,122,0.08)",
    name="|ΔU|",
), row=3, col=2)
in_bub = df["in_bubble"].values
for i in range(1, len(t)):
    if in_bub[i-1] == 1 and in_bub[i] == 1:
        fig.add_vrect(x0=float(t.iloc[i-1]), x1=float(t.iloc[i]),
                      fillcolor="rgba(100,181,246,0.25)", line_width=0,
                      row=3, col=2)

# ⑦ Re_p  — dedicated panel
fig.add_trace(go.Scatter(
    x=t, y=df["Re_p"], mode="lines",
    line=dict(color=C["Re"], width=2.5),
    fill="tozeroy", fillcolor="rgba(171,71,188,0.10)",
    name="Re_p",
), row=4, col=1)
fig.add_hline(y=1.0,  line_dash="dot", line_color="gray",
              annotation_text="Re=1 (Stokes limit)", row=4, col=1)
fig.add_hline(y=1000, line_dash="dot", line_color="gray",
              annotation_text="Re=1000 (inertial)", row=4, col=1)

# ⑧ St / We_p / Ca / Fr  (log)
for key, col, label in [("St",C["St"],"St"),("We_p",C["We"],"We_p"),
                          ("Ca",C["Ca"],"Ca"),("Fr",C["Fr"],"Fr")]:
    fig.add_trace(go.Scatter(
        x=t, y=df[key].clip(lower=1e-12), mode="lines",
        line=dict(color=col, width=2), name=label,
    ), row=4, col=2)
fig.add_hline(y=1.0, line_dash="dot", line_color="gray",
              annotation_text="= 1", row=4, col=2)
fig.update_yaxes(type="log", row=4, col=2)

# ⑨ P_coll
fig.add_trace(go.Scatter(
    x=t, y=df["P_coll"], mode="lines",
    line=dict(color=C["P"], width=2.5),
    fill="tozeroy", fillcolor="rgba(102,187,106,0.12)",
    name="P_coll",
), row=5, col=1)
fig.add_hline(y=0.5, line_dash="dot", line_color="gray",
              annotation_text="P=0.5", row=5, col=1)

# ⑩ Bo & Ar  (constant values shown as horizontal lines)
Bo_val = float(df["Bo"].mean())
Ar_val = float(df["Ar"].mean())
fig.add_trace(go.Scatter(
    x=t, y=np.full(len(t), Bo_val), mode="lines",
    line=dict(color=C["Bo"], width=2.5, dash="dash"),
    name=f"Bo = {Bo_val:.1f}",
), row=5, col=2)
fig.add_trace(go.Scatter(
    x=t, y=np.full(len(t), Ar_val), mode="lines",
    line=dict(color=C["Ar"], width=2.5, dash="dot"),
    name=f"Ar = {Ar_val:.3g}",
), row=5, col=2)
fig.update_yaxes(type="log", row=5, col=2)

# ── Axis labels ──────────────────────────────────────────────────
for axis, title in {
    "xaxis":"x [m]",     "yaxis":"y [m]",
    "xaxis2":"t [s]",    "yaxis2":"y [m]",
    "xaxis3":"t [s]",    "yaxis3":"|U| [m/s]",
    "xaxis4":"t [s]",    "yaxis4":"y [m]",
    "xaxis5":"t [s]",    "yaxis5":"α [-]",
    "xaxis6":"t [s]",    "yaxis6":"|ΔU| [m/s]",
    "xaxis7":"t [s]",    "yaxis7":"Re_p [-]",
    "xaxis8":"t [s]",    "yaxis8":"[-]",
    "xaxis9":"t [s]",    "yaxis9":"P_coll [-]",
    "xaxis10":"t [s]",   "yaxis10":"[-]",
}.items():
    fig.update_layout(**{axis: {"title": {"text": title}}})

# ── Formula annotations below each subplot ───────────────────────
# Map (row,col) to normalised paper coordinates
# Each subplot occupies approx 1/ROWS height and 1/COLS width
subplot_x_centres = [0.225, 0.775]  # col 1, col 2
row_heights = 1.0 / ROWS
formula_y_offsets = []
for r in range(ROWS):
    # y of bottom of subplot in paper coords (plotly counts from bottom)
    y_bottom = 1.0 - (r + 1) * row_heights + 0.005
    formula_y_offsets.append(y_bottom)

annotations = []
for (r, c), formula in FORMULAS.items():
    annotations.append(dict(
        text=f"<i>{formula}</i>",
        x=subplot_x_centres[c-1],
        y=formula_y_offsets[r-1],
        xref="paper", yref="paper",
        showarrow=False,
        font=dict(size=9, color="#666666"),
        align="center",
        xanchor="center",
        yanchor="top",
    ))

fig.update_layout(
    title=dict(
        text=(
            "oilWaterAirParticle — Flotation & Surface Interaction Dashboard"
            "<br><span style='font-size:13px;font-weight:normal;color:#888'>"
            "Re_p · St · We · Bo · Ar · Ca · Fr · P_coll (Sutherland) | "
            "blue shading = particle inside bubble</span>"
        ),
        font=dict(size=18),
    ),
    width=1600, height=2000,
    legend=dict(orientation="h", x=0.5, xanchor="center", y=-0.02),
    paper_bgcolor="white", plot_bgcolor="white",
    annotations=annotations,
)

fig.write_image(OUTPUT_PNG, scale=2)
print(f"\nDashboard PNG → {OUTPUT_PNG}")

# ── Console summary ───────────────────────────────────────────────
print("\n" + "="*62)
print("TIME-AVERAGED FLOTATION NUMBERS")
print("="*62)
for c, formula in [
    ("Re_p",   "ρ_w·|ΔU|·d_p / μ_w"),
    ("St",     "ρ_p·d_p²·U_b / (9·μ_w·R_b)"),
    ("We_p",   "ρ_w·|ΔU|²·d_p / σ_aw"),
    ("Bo",     "(ρ_w−ρ_air)·g·d_b² / σ_aw"),
    ("Ar",     "ρ_w·|Δρ|·g·d_b³ / μ_w²"),
    ("Ca",     "μ_w·|ΔU| / σ_aw"),
    ("Fr",     "U_b / √(g·R_b)"),
    ("P_coll", "St / (St + 0.25)"),
]:
    print(f"  {c:<10s} = {df[c].mean():.4g:<12}  [{formula}]")

st=df["St"].mean(); we=df["We_p"].mean(); re=df["Re_p"].mean()
print("\nRegime:")
print(f"  Re_p = {re:.3g}  → " + (
    "Stokes viscous drag (Re<1)" if re<1 else
    "Transitional (1<Re<1000)"   if re<1000 else
    "Inertial drag (Re>1000)"))
print(f"  St   = {st:.3g}  → " + (
    "Particle follows streamlines (St<<1)" if st<0.1 else
    "Intermediate inertia"                  if st<1 else
    "High inertia — efficient bubble collision"))
print(f"  We_p = {we:.3g}  → " + (
    "Surface tension dominant (We<1)" if we<1 else
    "Interface deformation possible (We>1)"))
Bo_v = df["Bo"].mean()
print(f"  Bo   = {Bo_v:.1f}  → " + (
    "Surface tension dominates gravity (Bo<1)" if Bo_v<1 else
    f"Gravity dominates surface tension (Bo={Bo_v:.0f} >> 1)"))
