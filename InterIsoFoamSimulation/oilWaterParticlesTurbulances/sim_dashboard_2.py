import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# -- Physical parameters --
RHO_W    = 998.2
RHO_OIL  = 870.0
RHO_P    = 1000.0
MU_W     = 1.002e-3
MU_OIL   = 30e-3
NU_W     = MU_W / RHO_W
SIGMA    = 0.072
D_P      = 1e-3
R_P      = D_P / 2.0
M_TOTAL  = 5.23599e-08
N_P      = 40
ALPHA_W  = 0.8912
LOG_FILE = 'log.turbulentThreePhaseFoam'

U_REF   = 0.01
L_REF   = D_P
EPS_EST = U_REF**3 / L_REF
ETA_K   = (NU_W**3 / EPS_EST)**0.25
R_RE    = NU_W / U_REF

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
    'axes.titlesize': 10,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
})

PAT = {
    'time': re.compile(r'^Time\s*=\s*([\d.eE+\-]+)'),
    'co':   re.compile(r'Courant Number mean:\s+([\d.eE+\-]+)'),
    'ke':   re.compile(r'Linear kinetic energy\s+=\s+([\d.eE+\-]+)'),
    'mom':  re.compile(r'\|Linear momentum\|\s+=\s+([\d.eE+\-]+)'),
    'prgh': re.compile(r'Solving for p_rgh, Initial residual\s*=\s*([\d.eE+\-]+)'),
    'ux':   re.compile(r'Solving for Ux, Initial residual\s*=\s*([\d.eE+\-]+)'),
    'dt':   re.compile(r'^deltaT\s*=\s*([\d.eE+\-]+)'),
}

def parse_log(path):
    rows, cur = [], {}
    with open(path) as f:
        for line in f:
            ls = line.strip()
            for key, pat in PAT.items():
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

def ax_plot(ax, x, y, color, xlabel, ylabel, title, formula):
    ax.plot(x, y, color=color, linewidth=1.3, marker='o', markersize=1.2)
    ax.fill_between(x, y, alpha=0.12, color=color)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.set_title(title, color='#111111', pad=4, fontsize=10)
    ax.grid(True, linestyle='--', alpha=0.4, which='both')
    ax.minorticks_on()
    ax.grid(True, which='minor', linestyle=':', alpha=0.2)
    ax.set_xlim(x[0], x[-1])
    for sp in ax.spines.values(): sp.set_edgecolor('#444444')
    ax.text(0.5, -0.28, formula, transform=ax.transAxes,
            ha='center', va='top', fontsize=7, color='#555555', style='italic')

def info_box(ax, eta_k, r_re):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    ax.set_facecolor('white')
    lines_text = [
        '=== Simulation Parameters ===',
        '',
        '-- Geometry --',
        'Particle diameter   d_p = 1.00e-3 m',
        'Particle radius     r_p = 5.00e-4 m',
        'Particle mass    m_tot = 5.24e-8 kg',
        'Parcels               N = 40',
        '',
        '-- Densities --',
        'Water          rho_w = 998.2 kg/m^3',
        'Oil            rho_o = 870.0 kg/m^3',
        'Particle       rho_p = 1000  kg/m^3',
        '',
        '-- Fluid Properties --',
        'nu_water = 1.004e-6 m^2/s',
        'sigma    = 0.072 N/m',
        'alpha_w  = 0.8912',
        '',
        '-- Turbulence Scales --',
        'eta_K  = %.3e m' % eta_k,
        '       (nu^3/eps)^(1/4)',
        'r_Re   = %.3e m' % r_re,
        '       (nu / U_ref)',
    ]
    y0, dy = 0.97, 0.055
    for i, ln in enumerate(lines_text):
        col = '#b45000' if ln.startswith('===') else (
              '#005fa3' if ln.startswith('--') else '#111111')
        ax.text(0.04, y0 - i * dy, ln, transform=ax.transAxes,
                ha='left', va='top', fontsize=7.2, color=col,
                fontfamily='monospace')

def cd_panel(ax, re_sim):
    # Schiller-Naumann reference curve
    re_ref = np.logspace(-2, 3, 400)
    cd_sn  = (24.0 / re_ref) * (1.0 + 0.15 * re_ref**0.687)
    cd_stokes = 24.0 / re_ref

    ax.loglog(re_ref, cd_sn,     color='#e3b341', lw=1.8,
              label='Schiller-Naumann')
    ax.loglog(re_ref, cd_stokes, color='#58a6ff', lw=1.2,
              ls='--', label='Stokes  $C_D=24/Re$')

    # Simulation trajectory as scatter (Re vs CD from simulation)
    re_sim_f = re_sim[np.isfinite(re_sim) & (re_sim > 1e-6)]
    cd_sim   = (24.0 / re_sim_f) * (1.0 + 0.15 * re_sim_f**0.687)
    sc = ax.scatter(re_sim_f, cd_sim, c=np.linspace(0, 1, len(re_sim_f)),
                    cmap='plasma', s=4, zorder=5, label='Simulation (colored by time)')
    plt.colorbar(sc, ax=ax, label='$t/t_{end}$', pad=0.02)

    ax.set_xlabel('$Re_p$', fontsize=9)
    ax.set_ylabel('$C_D$', fontsize=9)
    ax.set_title('Drag Coefficient vs. Particle Reynolds Number', color='#111111',
                 pad=4, fontsize=10)
    ax.legend(fontsize=7, facecolor='white', edgecolor='#aaaaaa')
    ax.grid(True, which='both', linestyle='--', alpha=0.4)
    ax.set_xlim([1e-3, 1e3])
    ax.set_ylim([1e0, 1e5])
    for sp in ax.spines.values(): sp.set_edgecolor('#444444')
    ax.text(0.5, -0.18,
            r'$C_D = \dfrac{24}{Re_p}(1 + 0.15\,Re_p^{0.687})$  [Schiller-Naumann]',
            transform=ax.transAxes, ha='center', va='top',
            fontsize=7, color='#555555', style='italic')

def main():
    rows = parse_log(LOG_FILE)
    if not rows:
        print('No data found.')
        return
    nan = float('nan')
    t    = np.array([r['time']          for r in rows])
    ke   = np.array([r.get('ke',   nan) for r in rows])
    mom  = np.array([r.get('mom',  nan) for r in rows])
    co   = np.array([r.get('co',   nan) for r in rows])
    pr   = np.array([r.get('prgh', nan) for r in rows])
    ux   = np.array([r.get('ux',   nan) for r in rows])
    dt   = np.array([r.get('dt',   nan) for r in rows])
    Up   = np.sqrt(2.0 * ke / M_TOTAL)
    we   = RHO_W * Up**2 * D_P / SIGMA
    re   = Up * D_P / NU_W
    pe   = Up * D_P / NU_W
    n    = len(t)
    t_norm = np.arange(n) / float(n - 1)

    fig = plt.figure(figsize=(22, 18), constrained_layout=False)
    fig.patch.set_facecolor('white')
    fig.subplots_adjust(left=0.05, right=0.97, top=0.94, bottom=0.08,
                        hspace=0.70, wspace=0.38)
    gs = gridspec.GridSpec(4, 4, figure=fig)

    # Row 0
    ax_info = fig.add_subplot(gs[0, 0])
    info_box(ax_info, ETA_K, R_RE)

    a = fig.add_subplot(gs[0, 1:3])
    ax_plot(a, t_norm, we, '#e3b341', '$t/t_{end}$', '$We$',
            'Average Weber Number',
            r'$We = \dfrac{\rho_f \, U_p^2 \, d_p}{\sigma}$')
    a.set_xlim([0.001, 1])
    a.set_ylim([0, max(0.1, float(np.nanmax(we)) * 1.1)])
    a.axhline(1.0, color='#e3b341', lw=0.8, ls='--', alpha=0.6, label='We=1')
    a.legend(fontsize=7, facecolor='white', edgecolor='#aaaaaa', labelcolor='#111111')

    ax_plot(fig.add_subplot(gs[0, 3]), t_norm, ke, '#ffa657',
            '$t/t_{end}$', 'KE [J]', 'Particle Kinetic Energy',
            r'$KE = \dfrac{1}{2} \, m \, |U_p|^2$')

    # Row 1
    ax_plot(fig.add_subplot(gs[1, 0]), t_norm, re, '#f78166',
            '$t/t_{end}$', '$Re_p$', 'Particle Reynolds Number',
            r'$Re_p = \dfrac{U_p \, d_p}{\nu}$')
    ax_plot(fig.add_subplot(gs[1, 1]), t_norm, pe, '#d2a8ff',
            '$t/t_{end}$', '$Pe_p$', 'Particle Peclet Number',
            r'$Pe_p = \dfrac{U_p \, d_p}{D}$')
    ax_plot(fig.add_subplot(gs[1, 2]), t_norm, Up, '#79c0ff',
            '$t/t_{end}$', '|Up| [m/s]', 'Particle Speed',
            r'$|U_p| = \sqrt{\dfrac{2\,KE}{m_{total}}}$')
    ax_plot(fig.add_subplot(gs[1, 3]), t_norm, mom, '#3fb950',
            '$t/t_{end}$', 'p [kg m/s]', '|Linear Momentum|',
            r'$\mathbf{p} = m_{total} \cdot |U_p|$')

    # Row 2
    ax_plot(fig.add_subplot(gs[2, 0]), t, co, '#58a6ff',
            't [s]', 'Co [-]', 'Courant Number',
            r'$Co = \dfrac{|U| \, \Delta t}{\Delta x}$')
    ax_plot(fig.add_subplot(gs[2, 1]), t, pr, '#f85149',
            't [s]', 'residual', 'p_rgh Initial Residual',
            r'$r_0 = \|b - A x_0\| / \|b\|$  (DICPCG)')
    ax_plot(fig.add_subplot(gs[2, 2]), t, ux, '#bc8cff',
            't [s]', 'residual', 'Ux Initial Residual',
            r'$r_0 = \|b - A x_0\| / \|b\|$  (smoothSolver)')
    ax_plot(fig.add_subplot(gs[2, 3]), t, dt, '#39d353',
            't [s]', 'dT [s]', 'Adaptive Time Step',
            r'$\Delta t^{n+1} = \Delta t^n \cdot Co_{max} / Co^n$')

    # Row 3 -- CD vs Re_p  (wide panel spanning all 4 cols)
    ax_cd = fig.add_subplot(gs[3, :])
    cd_panel(ax_cd, re)

    fig.suptitle('turbulentThreePhaseFoam  --  oilWaterParticles',
                 fontsize=13, color='#111111', y=0.97)
    plt.savefig('simulation_dashboard.png', dpi=150, bbox_inches='tight',
                facecolor='white')
    print('Saved simulation_dashboard.png')
    plt.show()

main()
