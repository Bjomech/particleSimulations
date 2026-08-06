import re
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ── Case configuration ─────────────────────────────────────────────────
CASES = {
    'InterFoam':    'oilWaterParticlesInterFoam',
    'InterIsoFoam': 'oilWaterParticlesInterIsoFoam',
}
LOG_NAMES = {
    'InterFoam':    'log.turbulentThreePhaseFoam',
    'InterIsoFoam': 'log.interisofoamSolver',
}
CLOUD   = 'kinematicCloud'
COLORS  = {'InterFoam': '#2196f3', 'InterIsoFoam': '#ff6f00'}
MARKERS = {'InterFoam': 'o',       'InterIsoFoam': 's'}

# ── Physical constants ─────────────────────────────────────────────────
RHO_W  = 998.2
RHO_P  = 1000.0
MU_W   = 1.002e-3
NU_W   = MU_W / RHO_W
SIGMA  = 0.072
D_P    = 1e-3
N_P    = 40
M_ONE  = 5.23599e-08 / N_P
DX     = 3e-3
EPS_EST = 0.01**3 / D_P
ETA_K  = (NU_W**3 / EPS_EST)**0.25
R_RE   = NU_W / 0.01
TAU_P  = (RHO_P * D_P**2) / (18.0 * MU_W)

plt.rcParams.update({
    'text.usetex': False,
    'font.family': 'monospace',
    'figure.facecolor': 'white',
    'axes.facecolor': '#f6f8fa',
    'axes.edgecolor': '#444444',
    'axes.labelcolor': '#222222',
    'xtick.color': '#222222',
    'ytick.color': '#222222',
    'grid.color': '#cccccc',
    'text.color': '#111111',
    'axes.titlesize': 11,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'axes.titlepad': 10,
    'axes.labelpad': 8,
    'xtick.major.pad': 5,
    'ytick.major.pad': 5,
    'legend.fontsize': 9,
    'legend.borderpad': 0.6,
    'legend.labelspacing': 0.5,
    'legend.handlelength': 2.0,
})

# ── Log parser ─────────────────────────────────────────────────────────
LOG_PAT = {
    'time': re.compile(r'^Time\s*=\s*([\d.eE+\-]+)'),
    'co':   re.compile(r'Courant Number mean:\s+([\d.eE+\-]+)'),
    'prgh': re.compile(r'Solving for p_rgh, Initial residual\s*=\s*([\d.eE+\-]+)'),
    'ux':   re.compile(r'Solving for Ux, Initial residual\s*=\s*([\d.eE+\-]+)'),
    'dt':   re.compile(r'^deltaT\s*=\s*([\d.eE+\-]+)'),
    'we_log': re.compile(r'Average Weber Number:\s+([\d.eE+\-]+)'),
    'ke':   re.compile(r'Linear kinetic energy\s+=\s+([\d.eE+\-]+)'),
    'mom':  re.compile(r'\|Linear momentum\|\s+=\s+([\d.eE+\-]+)'),
}

def parse_log(path):
    if not os.path.exists(path):
        print('  WARNING: log not found:', path)
        return []
    rows, cur = [], {}
    with open(path) as f:
        for line in f:
            ls = line.strip()
            for key, pat in LOG_PAT.items():
                m = pat.search(ls)
                if m:
                    if key == 'time':
                        if cur and 'time' in cur:
                            rows.append(cur)
                        cur = {'time': float(m.group(1))}
                    elif key not in cur:
                        cur[key] = float(m.group(1))
                    break
    if cur and 'time' in cur:
        rows.append(cur)
    return rows

# ── Lagrangian field readers ───────────────────────────────────────────
def read_scalar_field(path):
    vals = []
    with open(path) as f:
        in_list = False
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s == '(':
                in_list = True
                continue
            if s == ')':
                break
            if in_list:
                try:
                    vals.append(float(s))
                except ValueError:
                    pass
    return np.array(vals)

def read_vector_field(path):
    vecs = []
    with open(path) as f:
        in_list = False
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s == '(':
                in_list = True
                continue
            if s == ')':
                break
            if in_list and s.startswith('('):
                parts = s.strip('()').split()
                try:
                    vecs.append([float(p) for p in parts])
                except ValueError:
                    pass
    return np.array(vecs) if vecs else np.empty((0, 3))

def get_time_dirs(case_dir):
    dirs = []
    for name in os.listdir(case_dir):
        p = os.path.join(case_dir, name)
        if os.path.isdir(p):
            try:
                float(name)
                dirs.append(name)
            except ValueError:
                pass
    return sorted(dirs, key=float)

def read_lagrangian(case_dir):
    time_dirs = get_time_dirs(case_dir)
    times, we, re_p, up_mag, mom = [], [], [], [], []
    for td in time_dirs:
        lag = os.path.join(case_dir, td, 'lagrangian', CLOUD)
        if not os.path.isdir(lag):
            continue
        t = float(td)
        we_file = os.path.join(lag, 'We')
        re_file = os.path.join(lag, 'Re')
        u_file  = os.path.join(lag, 'U')
        we_val = re_val = up_val = None
        if os.path.exists(we_file):
            v = read_scalar_field(we_file)
            if len(v):
                we_val = float(np.mean(v))
        if os.path.exists(re_file):
            v = read_scalar_field(re_file)
            if len(v):
                re_val = float(np.mean(v))
        if os.path.exists(u_file):
            v = read_vector_field(u_file)
            if len(v):
                mag = np.linalg.norm(v, axis=1)
                up_val = float(np.mean(mag))
                mom.append(M_ONE * N_P * up_val)
        if we_val is None and up_val is not None:
            we_val = RHO_W * up_val**2 * D_P / SIGMA
        if re_val is None and up_val is not None:
            re_val = up_val * D_P / NU_W
        if we_val is None:
            continue
        times.append(t)
        we.append(we_val)
        re_p.append(re_val if re_val is not None else float('nan'))
        up_mag.append(up_val if up_val is not None else float('nan'))
        if up_val is None:
            mom.append(float('nan'))
    return (np.array(times), np.array(we), np.array(re_p),
            np.array(up_mag), np.array(mom))

def build_dataset(case_dir, log_path):
    nan = float('nan')
    times_lag, we, re_p, up_mag, mom_lag = read_lagrangian(case_dir)
    has_lag = len(times_lag) > 0
    rows    = parse_log(log_path)
    t_log   = np.array([r['time']            for r in rows]) if rows else np.array([])
    co      = np.array([r.get('co',   nan)   for r in rows]) if rows else np.array([])
    pr      = np.array([r.get('prgh', nan)   for r in rows]) if rows else np.array([])
    ux      = np.array([r.get('ux',   nan)   for r in rows]) if rows else np.array([])
    dt      = np.array([r.get('dt',   nan)   for r in rows]) if rows else np.array([])
    we_log  = np.array([r.get('we_log', nan) for r in rows]) if rows else np.array([])
    mom_log = np.array([r.get('mom',  nan)   for r in rows]) if rows else np.array([])
    ke_log  = np.array([r.get('ke',   nan)   for r in rows]) if rows else np.array([])
    if not has_lag:
        print('  No lagrangian data -- falling back to log')
        up_log  = np.sqrt(np.clip(2.0 * ke_log / (M_ONE * N_P), 0, None))
        we      = np.where(np.isfinite(we_log), we_log,
                           RHO_W * up_log**2 * D_P / SIGMA)
        re_p    = up_log * D_P / NU_W
        up_mag  = up_log
        mom_lag = mom_log
        times_lag = t_log
    if len(t_log) > 1 and len(times_lag) > 0:
        dt_interp = np.interp(times_lag, t_log, dt)
    else:
        dt_interp = np.full_like(times_lag, nan)
    dt_s   = np.where(dt_interp > 0, dt_interp, nan)
    pe_num = up_mag * D_P / (DX**2 / dt_s)
    a_p    = (np.gradient(up_mag, times_lag)
              if len(times_lag) > 1 else np.zeros_like(up_mag))
    f_net  = M_ONE * N_P * a_p
    n = len(times_lag)
    t_norm = np.arange(n) / float(max(n - 1, 1))
    return dict(
        source='lagrangian' if has_lag else 'log',
        t=times_lag, t_norm=t_norm,
        t_log=t_log, co=co, pr=pr, ux=ux, dt=dt,
        we=we, re_p=re_p, up=up_mag, mom=mom_lag,
        pe_num=pe_num, a_p=a_p, f_net=f_net,
    )

# ── Helpers ────────────────────────────────────────────────────────────
def style_ax(ax):
    ax.grid(True, linestyle='--', alpha=0.45, which='major')
    ax.minorticks_on()
    ax.grid(True, which='minor', linestyle=':', alpha=0.2)
    for sp in ax.spines.values():
        sp.set_edgecolor('#444444')

def add_formula(ax, text, dy=-0.30):
    ax.text(0.5, dy, text, transform=ax.transAxes,
            ha='center', va='top', fontsize=8,
            color='#555555', style='italic', clip_on=False)

def lag_panel(ax, datasets, key, ylabel, title, formula,
              use_tnorm=True, ref=None, ref_label=None, logy=False):
    for name, d in datasets.items():
        x = d['t_norm'] if use_tnorm else d['t']
        y = d[key]
        if not len(x):
            continue
        ax.plot(x, y, color=COLORS[name], lw=1.5,
                marker=MARKERS[name], markersize=1.2,
                label='%s (%s)' % (name, d['source']), alpha=0.88)
    if ref is not None:
        ax.axhline(ref, color='#e05252', lw=1.0, ls='--',
                   alpha=0.75, label=ref_label or '')
    if logy:
        ax.set_yscale('log')
    ax.set_xlabel('t/t_end' if use_tnorm else 't [s]', fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title(title, pad=10)
    if use_tnorm:
        ax.set_xlim([0.0, 1.0])
    ax.legend(facecolor='white', edgecolor='#aaaaaa',
              labelcolor='#111111', loc='best')
    style_ax(ax)
    add_formula(ax, formula)

def log_panel(ax, datasets, key, ylabel, title, formula, logy=False):
    for name, d in datasets.items():
        x = d['t_log']
        y = d[key]
        if not len(x):
            continue
        ax.plot(x, y, color=COLORS[name], lw=1.5,
                marker=MARKERS[name], markersize=1.2,
                label=name, alpha=0.88)
    if logy:
        ax.set_yscale('log')
    ax.set_xlabel('t [s]', fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title(title, pad=10)
    ax.legend(facecolor='white', edgecolor='#aaaaaa', loc='best')
    style_ax(ax)
    add_formula(ax, formula)

def info_panel(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    ax.set_facecolor('white')
    lines = [
        '=== Case Comparison ===', '',
        '-- Solvers --',
        '  A (blue):',
        '  turbulentThreePhaseFoam',
        '  B (orange):',
        '  interisofoamSolver',
        '', '-- Particle data --',
        '  Source: Lagrangian fields',
        '  (fallback: log file)',
        '', '-- Solver diagnostics --',
        '  Source: log file only',
        '', '-- Parameters --',
        '  d_p   = 1.00e-3 m',
        '  rho_p = 1000  kg/m^3',
        '  rho_w = 998.2 kg/m^3',
        '  nu_w  = 1.004e-6 m^2/s',
        '  sigma = 0.072 N/m',
        '  N_p   = 40 parcels',
        '  dx    = 3.0e-3 m',
        '', '-- Kolmogorov scales --',
        '  eta_K = %.2e m' % ETA_K,
        '  r_Re  = %.2e m' % R_RE,
        '  tau_p = %.2e s' % TAU_P,
    ]
    y0, dy = 0.97, 0.038
    for i, ln in enumerate(lines):
        col = '#b45000' if ln.startswith('===') else (
              '#005fa3' if ln.startswith('--') else '#111111')
        bold = ln.startswith('===') or ln.startswith('--')
        ax.text(0.04, y0 - i * dy, ln,
                transform=ax.transAxes,
                ha='left', va='top',
                fontsize=8.0, color=col,
                fontfamily='monospace',
                fontweight='bold' if bold else 'normal')

def cd_panel(ax, datasets):
    re_ref = np.logspace(-2, 4, 600)
    ax.loglog(re_ref, (24./re_ref)*(1+0.15*re_ref**0.687),
              color='#333333', lw=1.8, label='Schiller-Naumann')
    ax.loglog(re_ref, 24./re_ref,
              color='#888888', lw=1.1, ls='--', label='Stokes')
    for name, d in datasets.items():
        r = d['re_p']
        r = r[np.isfinite(r) & (r > 1e-9)]
        if not len(r):
            continue
        cd = (24./r)*(1+0.15*r**0.687)
        ax.scatter(r, cd, c=COLORS[name], s=8, zorder=5,
                   label=name, alpha=0.75, marker=MARKERS[name])
    ax.set_xlabel('Re_p  [-]', fontsize=10)
    ax.set_ylabel('C_D  [-]', fontsize=10)
    ax.set_title('Drag Coefficient vs. Particle Reynolds Number', pad=10)
    ax.legend(facecolor='white', edgecolor='#aaaaaa', fontsize=9)
    ax.grid(True, which='both', linestyle='--', alpha=0.4)
    ax.set_xlim([1e-2, 1e4])
    ax.set_ylim([1e-1, 1e5])
    style_ax(ax)
    add_formula(ax, 'C_D = 24/Re_p*(1+0.15*Re_p^0.687)  [Schiller-Naumann]', dy=-0.07)

def main():
    datasets = {}
    for name, case_dir in CASES.items():
        print('Reading:', name)
        log_path = os.path.join(case_dir, LOG_NAMES[name])
        datasets[name] = build_dataset(case_dir, log_path)
        print('  source:', datasets[name]['source'],
              '| timesteps:', len(datasets[name]['t']))

    fig = plt.figure(figsize=(26, 34))
    fig.patch.set_facecolor('white')
    gs = gridspec.GridSpec(
        6, 4, figure=fig,
        left=0.07, right=0.97,
        top=0.95, bottom=0.04,
        hspace=1.30, wspace=0.40,
    )

    # Row 0: Info | We (2 cols) | Re_p
    info_panel(fig.add_subplot(gs[0, 0]))
    ax_we = fig.add_subplot(gs[0, 1:3])
    lag_panel(ax_we, datasets, 'we',
              'We  [-]', 'Weber Number  [Lagrangian]',
              'We = rho_f * |Up - Uf|^2 * dp / sigma',
              ref=1.0, ref_label='We = 1')
    ax_we.set_xlim([0.001, 1.0])
    lag_panel(fig.add_subplot(gs[0, 3]), datasets, 're_p',
              'Re_p  [-]', 'Particle Reynolds Number  [Lagrangian]',
              'Re_p = |Up| * dp / nu_f')

    # Row 1: |Up| | a_p | F_net | Momentum
    lag_panel(fig.add_subplot(gs[1, 0]), datasets, 'up',
              '|Up|  [m/s]', 'Particle Speed  [Lagrangian]',
              '|Up| from lagrangian U field')
    lag_panel(fig.add_subplot(gs[1, 1]), datasets, 'a_p',
              'a_p  [m/s^2]', 'Particle Acceleration',
              'a_p = dUp/dt', ref=0.0, ref_label='a_p = 0')
    lag_panel(fig.add_subplot(gs[1, 2]), datasets, 'f_net',
              'F_net  [N]', 'Net Force',
              'F_net = m_total * a_p', ref=0.0, ref_label='F_net = 0')
    lag_panel(fig.add_subplot(gs[1, 3]), datasets, 'mom',
              'p  [kg m/s]', '|Linear Momentum|  [Lagrangian]',
              'p = m_total * |Up|')

    # Row 2: Pe_num | Co | dT | Stokes Number
    lag_panel(fig.add_subplot(gs[2, 0]), datasets, 'pe_num',
              'Pe_num  [-]', 'Numerical Peclet Number',
              'Pe_num = Up*dp*dt / dx^2', ref=1.0, ref_label='Pe_num = 1')
    log_panel(fig.add_subplot(gs[2, 1]), datasets, 'co',
              'Co  [-]', 'Courant Number  [Log]',
              'Co = |U| * dt / dx')
    log_panel(fig.add_subplot(gs[2, 2]), datasets, 'dt',
              'dT  [s]', 'Adaptive Time Step  [Log]',
              'dT^(n+1) = dT^n * Co_max / Co^n')
    ax_st = fig.add_subplot(gs[2, 3])
    tau_k = np.sqrt(NU_W / max(EPS_EST, 1e-30))
    for name, d in datasets.items():
        re = d['re_p']
        cd = np.where(re > 1e-9, (24./re)*(1+0.15*re**0.687), np.nan)
        tau_p_loc = (4./3.) * (RHO_P/RHO_W) * (D_P / np.maximum(d['up'], 1e-12)) / cd
        st = tau_p_loc / tau_k
        ax_st.plot(d['t_norm'], st, color=COLORS[name], lw=1.5,
                   marker=MARKERS[name], markersize=1.2,
                   label=name, alpha=0.88)
    ax_st.axhline(1.0, color='#e05252', lw=1.0, ls='--', label='St = 1')
    ax_st.set_xlabel('t/t_end', fontsize=10)
    ax_st.set_ylabel('St  [-]', fontsize=10)
    ax_st.set_title('Stokes Number  [derived]', pad=10)
    ax_st.set_xlim([0, 1])
    ax_st.legend(facecolor='white', edgecolor='#aaaaaa', loc='best')
    style_ax(ax_st)
    add_formula(ax_st, 'St = tau_p / tau_Kolmogorov')

    # Row 3: p_rgh (2 cols) | Ux (2 cols)
    log_panel(fig.add_subplot(gs[3, 0:2]), datasets, 'pr',
              'residual', 'p_rgh Initial Residual  [Log]',
              'r0 = ||b - A*x0|| / ||b||  (DICPCG)', logy=True)
    log_panel(fig.add_subplot(gs[3, 2:4]), datasets, 'ux',
              'residual', 'Ux Initial Residual  [Log]',
              'r0 = ||b - A*x0|| / ||b||  (smoothSolver)', logy=True)

    # Row 4: C_D vs Re_p (full width)
    cd_panel(fig.add_subplot(gs[4, :]), datasets)

    # Row 5: Phase diagram We vs Re_p | F_net vs |Up|
    ax_ph = fig.add_subplot(gs[5, 0:2])
    for name, d in datasets.items():
        re = d['re_p']
        we = d['we']
        mask = np.isfinite(re) & np.isfinite(we)
        ax_ph.scatter(re[mask], we[mask], c=COLORS[name],
                      s=10, alpha=0.6, label=name, marker=MARKERS[name])
    ax_ph.axhline(1.0, color='#e05252', lw=1.0, ls='--', label='We = 1')
    ax_ph.axvline(1.0, color='#888888', lw=1.0, ls=':', label='Re_p = 1')
    ax_ph.set_xlabel('Re_p  [-]', fontsize=10)
    ax_ph.set_ylabel('We  [-]', fontsize=10)
    ax_ph.set_title('Phase Diagram: We vs. Re_p', pad=10)
    ax_ph.legend(facecolor='white', edgecolor='#aaaaaa')
    style_ax(ax_ph)
    add_formula(ax_ph, 'We = rho_f*|dU|^2*dp/sigma   Re_p = |Up|*dp/nu')

    ax_fv = fig.add_subplot(gs[5, 2:4])
    for name, d in datasets.items():
        up = d['up']
        fn = d['f_net']
        mask = np.isfinite(up) & np.isfinite(fn)
        ax_fv.scatter(up[mask], fn[mask], c=COLORS[name],
                      s=10, alpha=0.6, label=name, marker=MARKERS[name])
    ax_fv.axhline(0, color='#e05252', lw=1.0, ls='--', label='F_net = 0')
    ax_fv.set_xlabel('|Up|  [m/s]', fontsize=10)
    ax_fv.set_ylabel('F_net  [N]', fontsize=10)
    ax_fv.set_title('Net Force vs. Particle Speed', pad=10)
    ax_fv.legend(facecolor='white', edgecolor='#aaaaaa')
    style_ax(ax_fv)
    add_formula(ax_fv, 'F_net = m_total * dUp/dt')

    fig.suptitle(
        'Solver Comparison:  turbulentThreePhaseFoam  vs.  interisofoamSolver\n'
        'Particle data from Lagrangian fields  |  Diagnostics from Log',
        fontsize=13, color='#111111', y=0.975, linespacing=1.5,
    )

    out = 'comparison_dashboard_2.png'
    plt.savefig(out, dpi=150, bbox_inches='tight', facecolor='white')
    print('Saved', out)

main()
