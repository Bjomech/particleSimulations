import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# -- Physical parameters --
RHO_W    = 998.2
RHO_OIL  = 870.0
RHO_P    = 1000.0
MU_W     = 1.002e-3
NU_W     = MU_W / RHO_W
SIGMA    = 0.072
D_P      = 1e-3
M_TOTAL  = 5.23599e-08
N_P      = 40
ALPHA_W  = 0.8912
L_REF    = D_P
DX       = 3e-3        # grid cell size estimate [m]  (from Co~0.003, U~0.01)
D_OIL_WATER = 7e-10    # molecular D_AB oil/water [m^2/s]
SC_OW    = NU_W / D_OIL_WATER
LOG_FILE = 'log.interisofoamSolver'

# Turbulence / relaxation scales
U_REF   = 0.01
EPS_EST = U_REF**3 / L_REF
ETA_K   = (NU_W**3 / EPS_EST)**0.25
R_RE    = NU_W / U_REF
TAU_P   = (RHO_P * D_P**2) / (18.0 * MU_W)

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
                        if cur and 'time' in cur: rows.append(cur)
                        cur = {'time': float(m.group(1))}
                    elif key not in cur:
                        cur[key] = float(m.group(1))
                    break
    if cur and 'time' in cur: rows.append(cur)
    return rows

def ax_plot(ax, x, y, color, xlabel, ylabel, title, formula, ref=None, ref_label=None):
    ax.plot(x, y, color=color, linewidth=1.3, marker='o', markersize=1.2)
    ax.fill_between(x, y, alpha=0.12, color=color)
    if ref is not None:
        ax.axhline(ref, color='#e05252', lw=1.1, ls='--', alpha=0.8,
                   label=ref_label or '')
        ax.legend(fontsize=7, facecolor='white', edgecolor='#aaaaaa')
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.set_title(title, color='#111111', pad=4, fontsize=10)
    ax.grid(True, linestyle='--', alpha=0.4, which='both')
    ax.minorticks_on()
    ax.grid(True, which='minor', linestyle=':', alpha=0.2)
    if len(x) > 1:
        ax.set_xlim(x[0], x[-1])
    for sp in ax.spines.values(): sp.set_edgecolor('#444444')
    ax.text(0.5, -0.30, formula, transform=ax.transAxes,
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
        'd_p  = 1.00e-3 m',
        'r_p  = 5.00e-4 m',
        'm    = 5.24e-8 kg',
        'N    = 40  parcels',
        'dx   = 3.0e-3 m (grid)',
        '',
        '-- Densities --',
        'rho_w = 998.2 kg/m^3',
        'rho_o = 870.0 kg/m^3',
        'rho_p = 1000  kg/m^3',
        '',
        '-- Fluid / Diffusion --',
        'nu_w  = 1.004e-6 m^2/s',
        'sigma = 0.072 N/m',
        'alpha_w = 0.8912',
        'D_AB  = 7.0e-10 m^2/s',
        'Sc_ow = %.0f' % SC_OW,
        '',
        '-- Scales (constant) --',
        'eta_K  = %.2e m' % eta_k,
        'r_Re   = %.2e m' % r_re,
        'tau_p  = %.2e s' % TAU_P,
        '',
        '-- Pe (molecular, ref.) --',
        'Pe_mol = Re * Sc_ow',
        '       = Re * %.0f' % SC_OW,
    ]
    y0, dy = 0.98, 0.046
    for i, ln in enumerate(lines_text):
        col = '#b45000' if ln.startswith('===') else (
              '#005fa3' if ln.startswith('--') else '#111111')
        ax.text(0.04, y0 - i * dy, ln, transform=ax.transAxes,
                ha='left', va='top', fontsize=6.9, color=col,
                fontfamily='monospace')

def cd_panel(ax, re_sim):
    re_ref    = np.logspace(-2, 3, 400)
    cd_sn     = (24.0 / re_ref) * (1.0 + 0.15 * re_ref**0.687)
    cd_stokes = 24.0 / re_ref
    ax.loglog(re_ref, cd_sn,     color='#e3b341', lw=1.8, label='Schiller-Naumann')
    ax.loglog(re_ref, cd_stokes, color='#58a6ff', lw=1.2, ls='--',
              label='Stokes  $C_D=24/Re$')
    re_f = re_sim[np.isfinite(re_sim) & (re_sim > 1e-9)]
    cd_f = (24.0 / re_f) * (1.0 + 0.15 * re_f**0.687)
    sc   = ax.scatter(re_f, cd_f, c=np.linspace(0, 1, len(re_f)),
                      cmap='plasma', s=4, zorder=5, label='Simulation')
    plt.colorbar(sc, ax=ax, label='$t/t_{end}$', pad=0.02)
    ax.set_xlabel('$Re_p$', fontsize=9)
    ax.set_ylabel('$C_D$', fontsize=9)
    ax.set_title('Drag Coefficient vs. Particle Reynolds Number',
                 color='#111111', pad=4)
    ax.legend(fontsize=7, facecolor='white', edgecolor='#aaaaaa')
    ax.grid(True, which='both', linestyle='--', alpha=0.4)
    ax.set_xlim([1e-3, 1e3])
    ax.set_ylim([1e0,  1e5])
    for sp in ax.spines.values(): sp.set_edgecolor('#444444')
    ax.text(0.5, -0.10,
            r'$C_D = \frac{24}{Re_p}(1+0.15\,Re_p^{0.687})$  [Schiller-Naumann]',
            transform=ax.transAxes, ha='center', va='top',
            fontsize=7, color='#555555', style='italic')

def main():
    rows = parse_log(LOG_FILE)
    if not rows:
        print('No data found.')
        return
    nan = float('nan')
    t   = np.array([r['time']          for r in rows])
    ke  = np.array([r.get('ke',   nan) for r in rows])
    mom = np.array([r.get('mom',  nan) for r in rows])
    co  = np.array([r.get('co',   nan) for r in rows])
    pr  = np.array([r.get('prgh', nan) for r in rows])
    ux  = np.array([r.get('ux',   nan) for r in rows])
    dt  = np.array([r.get('dt',   nan) for r in rows])

    Up  = np.sqrt(2.0 * ke / M_TOTAL)
    we  = RHO_W * Up**2 * D_P / SIGMA
    re  = Up * D_P / NU_W

    # -- New independent panels --
    # 1) Particle acceleration  a_p = dUp/dt
    dUp = np.gradient(Up, t)
    a_p = dUp

    # 2) Net force  F_net = m * a_p
    f_net = M_TOTAL * a_p

    # 3) Pe_num = Up * d_p / D_num  with  D_num = dx^2 / dt  (numerical diffusion)
    dt_safe = np.where(dt > 0, dt, np.nan)
    D_num   = DX**2 / dt_safe
    pe_num  = Up * D_P / D_num   # = Up * d_p * dt / dx^2

    n      = len(t)
    t_norm = np.arange(n) / float(n - 1)

    # -- Layout 5 x 4 --
    fig = plt.figure(figsize=(22, 22), constrained_layout=False)
    fig.patch.set_facecolor('white')
    fig.subplots_adjust(left=0.06, right=0.97, top=0.95, bottom=0.06,
                        hspace=0.72, wspace=0.38)
    gs = gridspec.GridSpec(5, 4, figure=fig)

    # Row 0: info | We (x2) | KE
    ax_info = fig.add_subplot(gs[0, 0])
    info_box(ax_info, ETA_K, R_RE)

    a = fig.add_subplot(gs[0, 1:3])
    ax_plot(a, t_norm, we, '#e3b341', '$t/t_{end}$', '$We$',
            'Average Weber Number',
            r'$We = \dfrac{\rho_f U_p^2 d_p}{\sigma}$')
    a.set_xlim([0.001, 1])
    a.set_ylim([0, max(0.1, float(np.nanmax(we)) * 1.1)])
    a.axhline(1.0, color='#e3b341', lw=0.8, ls='--', alpha=0.6, label='We=1')
    a.legend(fontsize=7, facecolor='white', edgecolor='#aaaaaa',
             labelcolor='#111111')

    ax_plot(fig.add_subplot(gs[0, 3]), t_norm, ke, '#ffa657',
            '$t/t_{end}$', 'KE [J]', 'Particle Kinetic Energy',
            r'$KE = \dfrac{1}{2} m |U_p|^2$')

    # Row 1: Re | a_p | F_net | |Up|
    ax_plot(fig.add_subplot(gs[1, 0]), t_norm, re, '#f78166',
            '$t/t_{end}$', '$Re_p$', 'Particle Reynolds Number',
            r'$Re_p = \dfrac{U_p d_p}{\nu}$')

    ax_plot(fig.add_subplot(gs[1, 1]), t_norm, a_p, '#ff9f43',
            '$t/t_{end}$', '$a_p$ [m/s$^2$]', 'Particle Acceleration',
            r'$a_p = \dfrac{\Delta U_p}{\Delta t}$  (numerical gradient)',
            ref=0.0, ref_label='$a_p = 0$')

    ax_plot(fig.add_subplot(gs[1, 2]), t_norm, f_net, '#e84393',
            '$t/t_{end}$', '$F_{net}$ [N]', 'Net Force on Particle',
            r'$F_{net} = m \cdot a_p$  (gravity - drag)',
            ref=0.0, ref_label='$F_{net}=0$  (equilibrium)')

    ax_plot(fig.add_subplot(gs[1, 3]), t_norm, Up, '#79c0ff',
            '$t/t_{end}$', '|Up| [m/s]', 'Particle Speed',
            r'$|U_p| = \sqrt{2\,KE / m_{total}}$')

    # Row 2: Pe_num | Drag Force | mom | dt
    ax_plot(fig.add_subplot(gs[2, 0]), t_norm, pe_num, '#80deea',
            '$t/t_{end}$', '$Pe_{num}$', 'Numerical Peclet Number',
            r'$Pe_{num} = \dfrac{U_p d_p}{D_{num}},\quad D_{num}=\dfrac{\Delta x^2}{\Delta t}$',
            ref=1.0, ref_label='$Pe_{num}=1$')

    cd_arr = np.where(re > 1e-10,
                      (24.0 / re) * (1.0 + 0.15 * re**0.687),
                      24.0 / 1e-10)
    A_p    = np.pi * (D_P / 2.0)**2
    fd     = 0.5 * cd_arr * RHO_W * A_p * Up**2
    ax_plot(fig.add_subplot(gs[2, 1]), t_norm, fd, '#e84370',
            '$t/t_{end}$', '$F_D$ [N]', 'Drag Force on Particle',
            r'$F_D = \dfrac{1}{2} C_D \rho_f A_p U_p^2$')

    ax_plot(fig.add_subplot(gs[2, 2]), t_norm, mom, '#3fb950',
            '$t/t_{end}$', 'p [kg m/s]', '|Linear Momentum|',
            r'$\mathbf{p} = m_{total} \cdot |U_p|$')

    ax_plot(fig.add_subplot(gs[2, 3]), t, dt, '#39d353',
            't [s]', 'dT [s]', 'Adaptive Time Step',
            r'$\Delta t^{n+1} = \Delta t^n \cdot Co_{max}/Co^n$')

    # Row 3: Co | p_rgh | Ux | tau_p textbox
    ax_plot(fig.add_subplot(gs[3, 0]), t, co, '#58a6ff',
            't [s]', 'Co [-]', 'Courant Number',
            r'$Co = |U|\,\Delta t / \Delta x$')
    ax_plot(fig.add_subplot(gs[3, 1]), t, pr, '#f85149',
            't [s]', 'residual', 'p_rgh Initial Residual',
            r'$r_0 = \|b - Ax_0\|/\|b\|$  (DICPCG)')
    ax_plot(fig.add_subplot(gs[3, 2]), t, ux, '#bc8cff',
            't [s]', 'residual', 'Ux Initial Residual',
            r'$r_0 = \|b - Ax_0\|/\|b\|$  (smoothSolver)')

    ax_tau = fig.add_subplot(gs[3, 3])
    ax_tau.axis('off')
    ax_tau.set_facecolor('white')
    tau_txt = [
        '=== Particle Relaxation Time ===',
        '',
        'tau_p = rho_p * d_p^2 / (18 * mu_f)',
        '',
        'tau_p = %.4e s' % TAU_P,
        '',
        'CONSTANT -- independent of U_p.',
        '',
        'tau_p << t_flow:',
        '  particle follows fluid instantly',
        'tau_p >> t_flow:',
        '  particle barely responds',
        '',
        '--- Numerical D_AB (reference) ---',
        'D_num = dx^2 / dt',
        'dx = %.1e m' % DX,
        'D_mol = %.1e m^2/s' % D_OIL_WATER,
        '',
        'Pe_num < 1:  diffusion dominates',
        'Pe_num > 1:  convection dominates',
    ]
    y0, dy = 0.97, 0.060
    for i, ln in enumerate(tau_txt):
        col = '#b45000' if ln.startswith('===') else (
              '#005fa3' if ln.startswith('---') else '#111111')
        ax_tau.text(0.04, y0 - i * dy, ln, transform=ax_tau.transAxes,
                    ha='left', va='top', fontsize=6.7, color=col,
                    fontfamily='monospace')

    # Row 4: CD vs Re_p (full width)
    ax_cd = fig.add_subplot(gs[4, :])
    cd_panel(ax_cd, re)

    fig.suptitle('interisofoamSolver  --  oilWaterParticles',
                 fontsize=14, color='#111111', y=0.97)
    plt.savefig('simulation_dashboard.png', dpi=150, bbox_inches='tight',
                facecolor='white')
    print('Saved simulation_dashboard.png')
    plt.show()

main()
