"""
Plotting functions

Sources to help create manifold plotting function:
Parker, J., & Anderson, R. (2013). LOW-ENERGY LUNAR TRAJECTORY DESIGN. https://descanso.jpl.nasa.gov/monograph/series12/LunarTraj--Overall.pdf

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from lagrange_point import lib_gen
from systems import get_system_dim
from dynamics import main_dynamics

# --------------------------------
# Hex colors for plotting
# --------------------------------
_BGC = "#000000" # black background
_GRID = "#313138" # dark grey grid
_PANE = "#000000" # pane fill
_MOON = "#c9c9c9" # moon color
_EARTH = "#00d2ff" # earth color
_L_COL = "#0eebf5" # Collinear Lagrange points (L1, L2, L3)
_L_TRI = "#ffcf56" # Triangular Lagrange points (L4, L5)
_CMAP = "plasma" # Color gradient map for halo family
_L_LABELS = ["$L_1$", "$L_2$", "$L_3$", "$L_4$", "$L_5$"] # List of strings for Lagrangian Points plotting

def plot_family(
        family,
        system,
        Lpt,
        color_by="index",
        cmap=_CMAP,
        show_primaries=True,
        show_Lpoint=True,
        show_shadow=True,
        show_stability=True,
        figsize=(11,9),
        save_path=None,
        dpi=200,
        elev=22,
        azim=-55
    ):
    """
    Plots the halo family with options

    Parameters:
    -------------
    family : dict
        Dictionary that holds all necessary information about the halo orbit
    system : str
        Name of system in which the family is propagated in
    Lpt : int
        Integer used to indicate which Lagrangian point the halo family is propagated from
    title : str (set to None)
        Title of halo family plot
    color_by : str
        The method in which the color gradient is made by
    cmap : str
        Color gradient that is used in halo family plot
    show_primaries : bool
        Boolean used to show primaries in system or not
    show_Lpoint : bool
        Boolean used to show Lagrangian points or not
    show_shadow : bool
        Boolean used to show xy projection of halo family on xy plane or not
    show_stability : bool
        Boolean used to plot stability index of halo family or not
    figsize : tuple
        Figure size
    save_path : str
        File path to save halo family plot to
    dpi : int
        Resolution of the saved image
    elev : float
        Initial viewing elevation in the halo family plot
    axim:  float
        Initial viewing azimuth in the halo family plot
    
    """
    if not family:
        raise ValueError("family list is empty")

    # Plotting parameters initialization
    cr3bp_sys = get_system_dim(system)
    mu = cr3bp_sys.mu  
    orbits = [m["states"] for m in family]
    jc_vals = np.array([m["jc"] for m in family])
    tf_vals = np.array([m["tf"] for m in family])
    periods = 2 * tf_vals
    az_vals = np.array([np.max(np.abs(m["states"][:,2])) for m in family])
    n = len(family)
    primary_1 = cr3bp_sys.primary
    primary_2 = cr3bp_sys.secondary

    scalar_map = {
        "index": np.arange(1, n+1, dtype=float),
        "jc": jc_vals,
        "period": periods,
        "az": az_vals
    }
    if color_by not in scalar_map:
        raise ValueError(f"color_by must be one of {list(scalar_map)}")
    values = scalar_map[color_by]
    norm = mcolors.Normalize(vmin=values.min(), vmax=values.max())
    cmap_obj = cm.get_cmap(cmap)
    colors = [cmap_obj(norm(v)) for v in values]

    label_map = {
        "index": "Family member index",
        "jc": "Jacobi constant $C_J$",
        "period": "Orbital period $T$ [TU]",
        "az": "z-amplitude $A_z$ [DU]",
    }

    fig = plt.figure(figsize=figsize, facecolor=_BGC)
    if show_stability:
        gs = fig.add_gridspec(2, 1, height_ratios=[3,1])
        ax = fig.add_subplot(gs[0], projection="3d", facecolor=_PANE)
        ax_stab = fig.add_subplot(gs[1], facecolor=_PANE)
    else:
        ax = fig.add_subplot(111, projection="3d", facecolor=_PANE)

    for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
        pane.fill = True
        pane.set_facecolor(_PANE)
        pane.set_edgecolor(_GRID)
    
    ax.grid(True, color=_GRID, linewidth=0.4, linestyle="--")
    ax.set_facecolor(_PANE)

    label_kw = dict(color="#d6d6d6", fontsize=9, labelpad=8)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.label.set_color("#d6d6d6")
        plt.setp(axis.get_ticklabels(), color="#b8b8b8", fontsize=8)
    
    ax.set_xlabel("x [DU]", **label_kw)
    ax.set_ylabel("y [DU]", **label_kw)
    ax.set_zlabel("z [DU]", **label_kw)

    # Initial viewing angles
    ax.view_init(elev=elev, azim=azim)

    # Using mirror configuration theorem to get full orbit
    def _full_orbit(half):
        mirror = half.copy()
        mirror[:, 1] *= -1
        mirror[:, 3] *= -1
        mirror[:, 4] *= -1
        return np.vstack([half, mirror[::-1]])

    # Initialize halo orbit plotting and shadow
    scale = 1.05
    full_orbits = [_full_orbit(orbit) for orbit in orbits]
    z_floor = np.min([orbit[:,2].min() for orbit in full_orbits]) * scale
    ax.set_zlim(z_floor, np.max([o[:,2].max() for o in full_orbits]) * scale)

    # Plot orbits and shadow
    for orbit, col in zip(full_orbits, colors):
        xs, ys, zs = orbit[:,0], orbit[:,1], orbit[:,2]

        if show_shadow:
            ax.plot(xs, ys, np.full_like(zs, z_floor),
                    color=col, lw=0.25, alpha=0.18, zorder=1)
        
        ax.plot(xs, ys, zs, color=col, lw=0.9, alpha=0.85, zorder=3)

    # Plot primaries
    if show_primaries: 
        ax.scatter([1 - mu], [0], [0], 
                   color=_MOON, s=10, marker="o",
                   zorder=10, depthshade=False)
        ax.text(1-mu, 0, 0, f"{primary_2}",
                color=_MOON, fontsize=9, va="center", zorder=11)
        if Lpt == 3:
            ax.scatter([-mu], [0], [0], 
                    color=_EARTH, s=10, marker="o",
                    zorder=10, depthshade=False)
            ax.text(-mu, 0, 0, f"{primary_1}",
                    color=_EARTH, fontsize=9, va="center", zorder=12)

    if show_Lpoint:
        xlocs, ylocs = lib_gen(mu)

        if Lpt in (1, 2, 3, 4, 5):
            ids = [Lpt - 1]
        else:
            ids = range(5)

        for i in ids:
            xL, yL = xlocs[i], ylocs[i]
            label = _L_LABELS[i]
            is_triangular = (i >= 3)
            col = _L_TRI if is_triangular else _L_COL
            marker = "^" if is_triangular else "*"
            size = 22 if is_triangular else 30

            ax.scatter([xL], [yL], [0], color=col, s=size, marker=marker,
                    zorder=10, depthshade=False)
            ax.text(xL, yL, 0, f"  {label}", color=col, fontsize=8,
                    va="center", zorder=11)
    
    if show_stability:
        stab_vals = np.array([fam["stability"] for fam in family])
        idxs = np.arange(1,n+1)

        imin = int(np.nanargmin(stab_vals))
        min_idx = idxs[imin]
        min_val = stab_vals[imin]

        box_txt = (
            f"Min stability\n"
            f"Index: {min_idx}\n"
            f"Value: {min_val:.10f}"
        )

        ax_stab.text(
            0.02, 0.98, box_txt,
            transform=ax_stab.transAxes,
            ha="left", va="top",
            fontsize=7.5, color="#cccccc",
            bbox=dict(
                boxstyle="round,pad=0.35",
                facecolor="#1a1a20",
                edgecolor="#333333",
                alpha=0.85
            )
        ) 
        
        pts = np.column_stack([idxs, stab_vals]).reshape(-1, 1, 2)
        segs = np.concatenate([pts[:-1], pts[1:]], axis=1)

        lc = LineCollection(segs, cmap=cmap_obj, norm=norm)
        lc.set_array(values[:-1])
        lc.set_linewidth(1.5)
        ax_stab.add_collection(lc)

        ax_stab.set_xlim(idxs.min(), idxs.max())
        hi = max(1.2, 1.05 * np.nanpercentile(stab_vals, 95))
        ax_stab.set_ylim(0.0, hi)
        ax_stab.scatter([min_idx], [min_val], s=55, facecolors="none",
                edgecolors="#ffffff", linewidths=1.2, zorder=10)
        ax_stab.set_xlabel(label_map["index"], color="#dedddd")
        ax_stab.set_ylabel("Stability Index", color="#dedddd")
        ax_stab.grid(True, color=_GRID, linestyle="--", lw=0.4)
        ax_stab.tick_params(colors="#b8b8b8")
        ax_stab.axhline(1.0, color="white", linestyle="--", alpha=0.4)

    sm = cm.ScalarMappable(cmap=cmap_obj, norm=norm)
    sm.set_array([])

    if show_stability:
        cbar = fig.colorbar(sm, ax=[ax, ax_stab], pad=0.02, shrink=0.55, aspect=22)
    else:
        cbar = fig.colorbar(sm, ax=ax, pad=0.02, shrink=0.55, aspect=22)

    cbar.set_label(label_map[color_by], color="#aaaaaa", fontsize=8.5)
    cbar.ax.yaxis.set_tick_params(color="#868585", labelsize=7.5)
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color="#868585")
    cbar.outline.set_edgecolor(_GRID)

    stats = (
        f"$\\mu$ = {mu:.6f}\n"
        f"Members: {n}\n"
        f"$A_z$: {az_vals.min():.3f} - {az_vals.max():.3f} DU\n"
        f"$T$: {periods.min():.3f} - {periods.max():.3f} TU\n"
        f"$C_J$: {jc_vals.min():.4f} - {jc_vals.max():.4f}"
    )
    ax.text2D(
        0.985, 0.03, stats,
        transform=ax.transAxes,
        fontsize=7, color="#999999",
        ha="right", va="bottom",
        bbox=dict(boxstyle="round,pad=0.4",
                  facecolor="#1a1a20", edgecolor="#333333",
                  alpha=0.7),
    )

    fig.suptitle(rf"CR3BP $L_{Lpt}$ Halo Orbit Family ($\mu$ = {mu:.9f})" + "\n" +
                 rf"System: {primary_1} - {primary_2}",
                 color="#dddddd", fontsize=12,
                 fontweight="semibold", y=0.97)
    
    fig.subplots_adjust(left=0.06, right=0.7, top=0.93, bottom=0.08, hspace=0.18)

    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight",
                    facecolor=_BGC, edgecolor="none")
        print(f"Figure saved -> {save_path}")
    
    plt.show()
    return fig, ax

def plot_manifold(system,
                  family,
                  member = 27,
                  prop_time_scale = 2.5,
                  gr=False,
                  show_earth=False):
    """
    Plots unstable and stable manifolds

    Parameter:
    -----------
    family : dict
        Dictionary that holds all necessary information about the halo orbit
    member : int
        Index of halo orbit within the halo orbit the manifolds are plotting from
    prop_time_scale : float
        A time scale of the half period of the halo orbit used to propagate manifold trajectories
    gr : bool
        Boolean used to determine if the grids and panes will be enabled or not
    """
    # Initializing states and STMs
    cr3bp_sys = get_system_dim(system)
    phi_tf_half = family[member-1]['stm'][-1]
    init_state  = family[member-1]['ic']
    tf_half     = family[member-1]['tf']
    mu          = cr3bp_sys.mu
    x, y, z, vx, vy, vz = init_state[:6]
    L = get_system_dim(system).L
    eps = (100/L)/np.sqrt(vx**2 + vy**2 + vz**2)
    primary_1 = cr3bp_sys.primary
    primary_2 = cr3bp_sys.secondary

    # full orbit with STMs (assumed: STM is Phi(t, t0) stacked in states.y[6:])
    orbit  = main_dynamics(mu, init_state, 2*tf_half, stm_bool=True)
    states = orbit.propagate()
    vectors = states.y[:6, :]                       # x(t_i)
    stms    = states.y[6:, :].T.reshape(-1, 6, 6)   # Phi(t_i, t0)

    # Build full-period monodromy from half-period using symmetry
    A = np.diag([1, -1, 1, -1, 1, -1])
    Phi_tf_full = A @ np.linalg.inv(phi_tf_half) @ A @ phi_tf_half

    # Eigenvectors and indices of (un)stable and center eigenvectors
    eigvals, eigvecs = np.linalg.eig(Phi_tf_full)
    unstable_idxs, stable_idxs, center_idxs = find_idxs(Phi_tf_full)
    stab_idx = np.max(np.abs(eigvals))
    # pick unstable eigenvector (column)
    u_idx = unstable_idxs[0]
    v_u = np.real(eigvecs[:, u_idx])
    v_u = v_u / np.linalg.norm(v_u)

    # pick stable eigenvector (column)
    s_idx = stable_idxs[0]
    v_s = np.real(eigvecs[:, s_idx])
    v_s = v_s / np.linalg.norm(v_s)

    # Plotting
    fig = plt.figure(figsize=[11,6], facecolor=_BGC)
    ax = fig.add_subplot(111, projection="3d", facecolor=_PANE)
    ax.plot(vectors[0,:], vectors[1,:], vectors[2,:], color='white', lw=4)

    for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
        pane.fill = True
        pane.set_facecolor(_PANE)
        pane.set_edgecolor(_GRID)

    ax.grid(True, color=_GRID, linewidth=0.4, linestyle="--")
    ax.set_facecolor(_PANE)

    label_kw = dict(color="#d6d6d6", fontsize=9, labelpad=8)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.label.set_color("#d6d6d6")
        plt.setp(axis.get_ticklabels(), color="#b8b8b8", fontsize=8)

    ax.set_xlabel("x [DU]", **label_kw)
    ax.set_ylabel("y [DU]", **label_kw)
    ax.set_zlabel("z [DU]", **label_kw)

    ax.scatter([1 - mu], [0], [0], color=_MOON, s=10, marker="o", zorder=10, depthshade=False)
    ax.text(1-mu, 0, 0, f"{primary_2}", color=_MOON, fontsize=9, va="center", zorder=11)

    if show_earth:
        ax.scatter([-mu], [0], [0], color=_EARTH, s=10, marker="o", zorder=10, depthshade=False)
        ax.text(-mu, 0, 0, f"{primary_1}", color=_EARTH, fontsize=9, va="center", zorder=12)

    fig.suptitle(
        rf"Stable and Unstable Manifold of Halo Orbit Family Member {member} ($\pm\,\epsilon = {eps:.4f}$)" + "\n" +
        rf"$t_{{prop}} = {prop_time_scale*tf_half:.5f}\,\mathrm{{TU}}$ or ${prop_time_scale/2} \times T_{{halo}}$ Stability Index = ${stab_idx}$",
        color="#dddddd",
        fontsize=12,
        fontweight="semibold",
        y=0.97
    )
    
    # Inject manifold off each point on the orbit correctly
    for i in range(len(states.t)):
        state_i = vectors[:, i]
        Phi_i   = stms[i]
        delta_u = Phi_i @ v_u
        delta_s = Phi_i @ v_s
        delta_u = delta_u / np.linalg.norm(delta_u)
        delta_s = delta_s / np.linalg.norm(delta_s)
        # two branches: +eps and -eps
        for sgn in (+1.0, -1.0):
            new_state_u = state_i + sgn * eps * delta_u
            new_orbit_u = main_dynamics(mu, new_state_u, prop_time_scale*tf_half, stm_bool=False)
            vis_states_u = new_orbit_u.propagate()
            xu, yu, zu = vis_states_u.y[0,:], vis_states_u.y[1,:], vis_states_u.y[2,:]

            new_state_s = state_i + sgn * eps * delta_s
            new_orbit_s = main_dynamics(mu, new_state_s, -prop_time_scale*tf_half, stm_bool=False)
            vis_states_s = new_orbit_s.propagate()
            xs, ys, zs = vis_states_s.y[0,:], vis_states_s.y[1,:], vis_states_s.y[2,:]

            if sgn == -1.0:
                color_u = '#0009ff'
                color_s = '#f300ff'
            else:
                color_u = '#ff0000'
                color_s = '#00db12'
            
            ax.plot(xu, yu, zu, color=color_u, alpha=0.4, lw=1)
            ax.plot(xs, ys, zs, color=color_s, alpha=0.4, lw=1)
    
    # Legend Creation
    legend_elements = [
    Line2D([0], [0], color='white', lw=6, label='Halo Orbit'),
    Line2D([0], [0], color='#ff0000',   lw=2, label=r'Unstable $W_{+}^{U}$ (+ε)'),
    Line2D([0], [0], color='#0009ff',  lw=2, label=r'Unstable $W_{-}^{U}$ (−ε)'),
    Line2D([0], [0], color='#00db12', lw=2, label=r'Stable $W_{+}^{S}$ (+ε)'),
    Line2D([0], [0], color='#f300ff',  lw=2, label=r'Stable $W_{-}^{S}$ (−ε)')
    ]

    ax.legend(
        handles=legend_elements,
        loc='upper left',
        frameon=True,
        facecolor=_PANE,
        edgecolor=_GRID,
        labelcolor='#dddddd',
        fontsize=9
    )

    # Deletes grid and panes if gr == False
    if gr == False:
        ax.grid(gr)
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis._axinfo["grid"]["linewidth"] = 0
        for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
            pane.set_visible(False)

    ax.view_init(azim=-90, elev=90)
    plt.show()

def find_idxs(monodromy):
    """
    Gets the indices of the unstable, stable and center eigvectors
    by searching the eigenvalues for the corresponding indices

    unstable: if the magnitude of the eigenvalue is less than -1 or greater than 1
    stable: if the magnitude of the eigenvalue is within [-1, 1]
    center: if magnitude of the eigenvalue is equal to 1

    Uses magnitudes helps with issues of tight margins in determining stable vs center
    eigenvectors

    Parameters:
    monodromy: ndarray (6x6)
        The state transition matrix of the full orbital period Φ(T,0)
    
    Returns:
    u_idx : ndarray
        The indices corresponding to the unstable eigenvectors
    s_idx : ndarray
        The indices corresponding to the stable eigenvectors
    c_idx : ndarray
        The indices corresponding to the center eigenvectors
    """
    tol = 1e-6
    eigvals, _ = np.linalg.eig(monodromy)

    mags = np.abs(eigvals)

    u_idx = [i for i, m in enumerate(mags) if m > 1 + tol]
    s_idx = [i for i, m in enumerate(mags) if m < 1 - tol]
    c_idx = [i for i, m in enumerate(mags) if abs(m - 1.0) <= tol]

    print("="*70)
    print("Eigenvalues (mag):", mags)
    print("Unstable:", u_idx)
    print("Stable:", s_idx)
    print("Center:", c_idx)
    print("="*70)

    return u_idx, s_idx, c_idx
