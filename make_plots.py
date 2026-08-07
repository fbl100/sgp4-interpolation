#!/usr/bin/env python3
"""Render validation figures (PNG) from data/results.json for the HTML report."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white",
    "font.size": 11, "axes.grid": True, "grid.alpha": 0.3,
    "axes.spines.top": False, "axes.spines.right": False,
})

RES = json.load(open("data/results.json"))
REG = {r["regime"]: r for r in RES["regimes"]}
ORDER = ["LEO", "MEO", "GEO", "HEO"]
COLOR = {"LEO": "#2563eb", "MEO": "#059669", "GEO": "#d97706", "HEO": "#dc2626"}
TARGET = RES["meta"]["target_m"]
NOM = RES["meta"]["nominal_h_s"]


def fig_convergence():
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
    for ax, reg in zip(axes.flat, ORDER):
        r = REG[reg]; sw = r["sweep"]
        h = np.array([s["h"] for s in sw])
        obs = np.array([s["max_m"] for s in sw])
        pred = np.array([s["pred_m"] for s in sw])
        cons = np.array([s["max_consistent_v_m"] for s in sw])
        c = COLOR[reg]
        ax.loglog(h, obs, "o-", color=c, lw=2, ms=6,
                  label="Observed max error (SGP4 reported v)")
        ax.loglog(h, pred, "--", color="#111", lw=1.5,
                  label="A-priori bound  jounce·h⁴/384 + 0.096·δᵥ·h")
        ax.loglog(h, cons, "s:", color="#7c3aed", lw=1.5, ms=4,
                  label="Consistent (differentiated) v → pure O(h⁴)")
        ax.axhline(TARGET, color="#888", lw=1, ls="-.")
        ax.text(h[0], TARGET*1.25, f"{TARGET:.0f} m budget", color="#666", fontsize=8)
        ax.axvline(NOM, color="#bbb", lw=1)
        st = r["slope_truncation"]
        sl = f"O(h⁴) slope (truncation region) ≈ {st:.2f}" if st else ""
        ax.set_title(f"{reg} — {r['name']}\n"
                     f"e={r['ecc']:.4f}, T={r['period_min']:.0f} min, "
                     f"h≤{r['required_h_for_target_s']:.0f}s for <1 m", fontsize=10)
        ax.set_xlabel("Node spacing h (s)")
        ax.set_ylabel("Max 3D position error (m)")
        if sl:
            ax.text(0.04, 0.94, sl, transform=ax.transAxes, fontsize=8.5,
                    va="top", color=c)
        if reg == "LEO":
            ax.legend(fontsize=7.5, loc="lower right", framealpha=0.9)
    fig.suptitle("Cubic Hermite interpolation error vs node spacing — real TLEs, SGP4 truth",
                 fontsize=13, y=0.998)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig("plots/convergence.png", dpi=130)
    plt.close(fig)


def fig_bound_scatter():
    """Observed vs predicted across ALL regimes & spacings — bound validity."""
    fig, ax = plt.subplots(figsize=(6.6, 6.0))
    lo, hi = 1e-2, 1e5
    ax.loglog([lo, hi], [lo, hi], color="#111", lw=1.2, label="observed = predicted")
    ax.fill_between([lo, hi], [lo, hi], [lo*1e-3, hi*1e-3], color="#e5f0ff",
                    alpha=0.5, label="observed ≤ bound (conservative)")
    for reg in ORDER:
        sw = REG[reg]["sweep"]
        obs = [s["max_m"] for s in sw]; pred = [s["pred_m"] for s in sw]
        ax.loglog(pred, obs, "o", color=COLOR[reg], ms=8, label=reg, alpha=0.9)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel("A-priori predicted bound E(h)  (m)")
    ax.set_ylabel("Observed max error  (m)")
    ax.set_title("Every point sits on or below the bound\n(4 regimes × 6 spacings, zero free parameters)")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig("plots/bound_scatter.png", dpi=130)
    plt.close(fig)


def fig_ric():
    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    x = np.arange(len(ORDER)); w = 0.26
    labels = ["Radial", "In-track", "Cross-track"]
    cols = ["#0ea5e9", "#ef4444", "#22c55e"]
    for k in range(3):
        vals = [REG[reg]["nominal"]["ric_max_m"][k] for reg in ORDER]
        vals = [abs(v) for v in vals]
        ax.bar(x + (k-1)*w, vals, w, label=labels[k], color=cols[k])
    ax.set_yscale("log")
    ax.set_xticks(x); ax.set_xticklabels([f"{r}\n{REG[r]['name']}" for r in ORDER], fontsize=9)
    ax.axhline(TARGET, color="#888", ls="-.", lw=1)
    ax.text(len(ORDER)-0.5, TARGET*1.1, "1 m budget", color="#666", fontsize=8, ha="right")
    ax.set_ylabel("Max |error component| at 120 s (m)")
    ax.set_title("Dominant error axis reveals the mechanism (RIC, 120 s nodes)\n"
                 "truncation-limited LEO/HEO → RADIAL · velocity-floor GEO → IN-TRACK",
                 fontsize=11)
    ax.legend()
    fig.tight_layout()
    fig.savefig("plots/ric.png", dpi=130)
    plt.close(fig)


def fig_timeseries():
    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    for ax, reg in zip(axes.flat, ORDER):
        s = RES["series"][reg]
        t = np.array(s["t_frac"]); e = np.array(s["err_m"])
        ax.plot(t, e, color=COLOR[reg], lw=0.9)
        ax.set_title(f"{reg} — {REG[reg]['name']}", fontsize=10)
        ax.set_xlabel("Time (orbital periods)")
        ax.set_ylabel("|error| at 120 s (m)")
        ax.axhline(TARGET, color="#888", ls="-.", lw=1)
        if reg == "HEO":
            ax.annotate("perigee passages\n(fast motion → error spikes)",
                        xy=(t[np.argmax(e)], np.max(e)), xytext=(0.35, 0.7),
                        textcoords="axes fraction", fontsize=8.5, color="#dc2626",
                        arrowprops=dict(arrowstyle="->", color="#dc2626"))
    fig.suptitle("Where the error lives over the orbit (120 s nodes, reported velocity)",
                 fontsize=12, y=0.999)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig("plots/timeseries.png", dpi=130)
    plt.close(fig)


def fig_velocity_compare():
    """Per-orbit: max error at 120 s using SGP4 reported velocity vs a consistent
    (differentiated) velocity for the Hermite node slopes."""
    fig, ax = plt.subplots(figsize=(9.2, 5.2))
    x = np.arange(len(ORDER)); w = 0.38
    rep = [REG[r]["nominal"]["max_m"] for r in ORDER]
    con = [REG[r]["nominal"]["max_consistent_v_m"] for r in ORDER]
    for i, reg in enumerate(ORDER):
        c = COLOR[reg]
        ax.bar(x[i]-w/2, rep[i], w, color=c,
               label="SGP4 reported velocity" if i == 0 else None)
        ax.bar(x[i]+w/2, con[i], w, color=c, alpha=0.4, hatch="////", edgecolor=c,
               label="Consistent (differentiated) velocity" if i == 0 else None)
        factor = rep[i] / con[i]
        ftxt = f"×{factor:,.0f}" if factor >= 10 else f"×{factor:.1f}"
        ax.annotate(ftxt, xy=(x[i], max(rep[i], con[i])), xytext=(0, 6),
                    textcoords="offset points", ha="center", fontsize=10,
                    fontweight="bold", color=c)
    ax.set_yscale("log")
    ax.set_ylim(3e-4, 5e1)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{r}\n{REG[r]['name']}" for r in ORDER], fontsize=9)
    ax.axhline(TARGET, color="#888", ls="-.", lw=1)
    ax.text(len(ORDER)-0.5, TARGET*1.15, "1 m budget", color="#666", fontsize=8, ha="right")
    ax.set_ylabel("Max 3D position error at 120 s (m)")
    ax.set_title("Node-velocity source, per orbit (120 s spacing)\n"
                 "consistent velocity removes the data floor → big win where the floor dominates (GEO, MEO)",
                 fontsize=11)
    ax.legend(fontsize=9, loc="upper left")
    fig.tight_layout()
    fig.savefig("plots/velocity_compare.png", dpi=130)
    plt.close(fig)


def fig_population():
    """Obs vs a-priori bound for ~40 diverse real objects x 3 spacings, colored
    by eccentricity. The robustness claim: the bound holds everywhere."""
    pop = json.load(open("data/population.json"))
    objs = pop["objects"]; m = pop["meta"]
    preds, obss, eccs, ref = [], [], [], []
    for o in objs:
        for s in o["samples"]:
            if s["pred"] > 0 and s["obs"] > 0:
                preds.append(s["pred"]); obss.append(s["obs"])
                eccs.append(o["ecc"]); ref.append(o["reference"])
    preds = np.array(preds); obss = np.array(obss)
    eccs = np.array(eccs); ref = np.array(ref)

    fig, ax = plt.subplots(figsize=(8.2, 6.8))
    lo, hi = 1e-4, 1e5
    ax.fill_between([lo, hi], [lo, hi], [lo*1e-4, hi*1e-4], color="#eef4ff",
                    alpha=0.7, zorder=0, label="observed ≤ bound (valid)")
    ax.plot([lo, hi], [lo, hi], color="#111", lw=1.2, zorder=1,
            label="observed = bound")
    sc = ax.scatter(preds[~ref], obss[~ref], c=eccs[~ref], cmap="viridis",
                    s=42, edgecolor="white", lw=0.4, vmin=0, vmax=0.92, zorder=2)
    ax.scatter(preds[ref], obss[ref], c=eccs[ref], cmap="viridis",
               s=150, marker="*", edgecolor="#111", lw=0.9, vmin=0, vmax=0.92,
               zorder=3, label="reference satellites (§ 01– 07)")
    cb = fig.colorbar(sc, ax=ax, pad=0.02); cb.set_label("orbit eccentricity")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel("A-priori predicted bound  E(h)  (m)")
    ax.set_ylabel("Observed max error  (m)")
    ax.set_title(f"The bound holds across {m['n_analyzed']} diverse real orbits\n"
                 f"{m['n_points']} (object × spacing) points, e = 0 → 0.91, "
                 f"LEO → 3.5-day HEO — {m['n_above_bound']} above the bound",
                 fontsize=11.5)
    ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout()
    fig.savefig("plots/population.png", dpi=130)
    plt.close(fig)


def fig_ecc_ladder():
    """How the velocity floor (left) and the payoff of a consistent velocity
    (right) scale with eccentricity, across all analyzed objects. Scatter, not a
    line: these are different real objects (period varies), so we show the trend,
    not a function. Labels on a few famous ones only."""
    pop = json.load(open("data/population.json"))
    objs = pop["objects"]
    e = np.array([o["ecc"] for o in objs])
    dv = np.array([o["dv_m_s"] for o in objs]) * 100.0        # cm/s
    adv = np.array([o["req_cons_s"] / o["req_sgp4_s"] for o in objs])
    reg = [o["regime"] for o in objs]
    cols = [COLOR[r] for r in reg]
    LABEL = {25544: "ISS", 26407: "GPS BIIR-5", 40129: "Galileo 6 (e≈0.17)",
             25989: "XMM-Newton", 25847: "Molniya 3-50", 27540: "INTEGRAL",
             26464: "Cluster II"}

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.6, 5.5))

    a1.scatter(e, dv, c=cols, s=48, edgecolor="white", lw=0.4, zorder=3)
    for o, x, y in zip(objs, e, dv):
        if o["num"] in LABEL:
            a1.annotate(LABEL[o["num"]], (x, y), fontsize=8, color="#444",
                        xytext=(5, 4), textcoords="offset points", zorder=4)
    a1.set_yscale("log")
    a1.set_xlabel("Orbit eccentricity")
    a1.set_ylabel("Node-velocity gap  δᵥ  (cm/s)")
    a1.set_title("The floor grows ~100× with eccentricity\n"
                 "SGP4 reported velocity vs. d/dt of SGP4 position", fontsize=11)
    a1.grid(alpha=0.3)

    a2.scatter(e, adv, c=cols, s=48, edgecolor="white", lw=0.4, zorder=3)
    for o, x, y in zip(objs, e, adv):
        if o["num"] in LABEL:
            a2.annotate(LABEL[o["num"]], (x, y), fontsize=8, color="#444",
                        xytext=(5, 4), textcoords="offset points", zorder=4)
    a2.axhline(1, color="#999", lw=1, ls="-.")
    a2.set_yscale("log")
    a2.set_xlabel("Orbit eccentricity")
    a2.set_ylabel("Consistent-v advantage  (× coarser spacing for <1 m)")
    a2.set_title("The payoff of a consistent velocity\ngrows with eccentricity",
                 fontsize=11)
    a2.grid(alpha=0.3)
    # regime legend
    from matplotlib.lines import Line2D
    handles = [Line2D([0],[0], marker="o", ls="", mfc=COLOR[r], mec="white",
                      ms=8, label=r) for r in ORDER]
    a1.legend(handles=handles, fontsize=8.5, loc="lower right", title="regime")
    fig.tight_layout()
    fig.savefig("plots/ecc_ladder.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    import os
    os.makedirs("plots", exist_ok=True)
    fig_convergence(); fig_bound_scatter(); fig_ric(); fig_timeseries()
    fig_velocity_compare(); fig_population(); fig_ecc_ladder()
    print("Wrote plots/*.png (incl. population.png, ecc_ladder.png)")
