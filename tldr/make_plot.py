#!/usr/bin/env python3
"""TL;DR calibration chart: max error vs node spacing (1-300 s), two panels
(raw SGP4 velocity | central-difference velocity), from tldr/sweep.json."""
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

SW = json.load(open("tldr/sweep.json"))
REG = {r["regime"]: r for r in SW["regimes"]}
ORDER = ["LEO", "MEO", "GEO", "HEO"]
COLOR = {"LEO": "#2563eb", "MEO": "#059669", "GEO": "#d97706", "HEO": "#dc2626"}
# plain-language descriptor for a non-orbit audience (matches the calibrator artifact)
DESC = {"LEO": "Low Earth, ~90 min · ISS, Starlink",
        "MEO": "Medium Earth, ~12 h · GPS/nav",
        "GEO": "Geostationary, 24 h · TV, weather",
        "HEO": "Highly elliptical · Molniya"}
H = np.array(SW["h_s"])
DT = SW["meta"]["dt_central_s"]
BUDGETS = [0.1, 1.0, 10.0]   # reference "error you can accept" lines


def panel(ax, key, title):
    for reg in ORDER:
        r = REG[reg]
        y = np.array(r[key])
        ax.loglog(H, y, "-", color=COLOR[reg], lw=2,
                  label=f"{reg} — {DESC[reg]}")
    for b in BUDGETS:
        ax.axhline(b, color="#888", lw=0.9, ls="--", alpha=0.7)
        ax.text(H[-1]*0.97, b*1.18, f"accept {b:g} m", color="#666", fontsize=8,
                ha="right")
    ax.axvline(120, color="#bbb", lw=1)
    ax.text(120, ax.get_ylim()[0], " common 120 s", color="#999", fontsize=8,
            rotation=90, va="bottom", ha="left")
    ax.set_xlim(H[0], H[-1])
    ax.set_xlabel("Spacing between stored points (s)")
    ax.set_ylabel("Worst-case reconstruction error (m)")
    ax.set_title(title, fontsize=11)


def main():
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.6), sharey=True)
    panel(axes[0], "err_raw_m",
          "As stored  (raw velocity from the propagator)")
    panel(axes[1], "err_diff_m",
          f"With smoothing  (node velocity re-derived from position, dt={DT:g} s)")
    axes[1].set_ylabel("")
    axes[0].legend(fontsize=8.5, loc="upper left", framealpha=0.95,
                   title="Find your budget on the y-axis → read the spacing off the x-axis",
                   title_fontsize=8.5)
    fig.suptitle("How far apart can stored orbit points be? — worst-case error when you "
                 "rebuild the path between them (cubic Hermite, SGP4 truth)",
                 fontsize=12.5, y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig("tldr/calibration.png", dpi=140)
    plt.close(fig)
    print("Wrote tldr/calibration.png")


if __name__ == "__main__":
    main()
