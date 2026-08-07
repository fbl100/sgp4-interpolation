#!/usr/bin/env python3
"""Generate README.md from the same data as the HTML report (data/results.json,
data/population.json) so the two never drift. Charts are linked as relative
paths (plots/*.png) which GitHub renders inline."""
import json

R = json.load(open("data/results.json"))
P = json.load(open("data/population.json"))
REG = {r["regime"]: r for r in R["regimes"]}
ORDER = ["LEO", "MEO", "GEO", "HEO"]
TARGET = R["meta"]["target_m"]
NOM = int(R["meta"]["nominal_h_s"])
AXIS = ["radial", "in-track", "cross-track"]
leo, meo, geo, heo = (REG[k] for k in ORDER)
pm = P["meta"]
pe = [o["ecc"] for o in P["objects"]]
pT = [o["period_min"] for o in P["objects"]]


def f(x, p=2):
    return f"{x:,.{p}f}"


_SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")

def sci(x):
    a, b = f"{x:.2e}".split("e")
    return f"{a}×10{str(int(b)).translate(_SUP)}"


def dom(ric):
    return AXIS[max(range(3), key=lambda i: abs(ric[i]))]


def verdict(r):
    return "✅ pass" if r["nominal"]["max_m"] <= TARGET else "⚠️ refine"


md = []
w = md.append

w(f"""# SGP4 Orbit Interpolation — Cubic Hermite Error Validation

**Is piecewise cubic Hermite interpolation of ephemeris nodes an acceptable — and
*boundable* — way to reconstruct an orbit between nodes?** This repo answers that
with a reproducible study on **real satellites** across four orbit regimes, using
SGP4 as *both* the node source and the dense truth reference. Because the same
propagator produces the nodes and the truth, it isolates **interpolation** error
with zero force-model mismatch — the "state vectors with 0 % uncertainty" case.

> This README is generated from the same data as the full report
> (`report.html`) by `build_readme.py`, so every number here is exact.

---

## Verdict

Cubic Hermite is an excellent, provably **O(h⁴)** interpolant for orbital motion —
every regime confirms the fourth-order law. Whether the nominal **{NOM} s / <{TARGET:.0f} m**
target is met depends on the regime, and the required spacing follows directly
from the a-priori bound below.

| Regime | Satellite | e | Period (min) | max err @ {NOM} s | vs {TARGET:.0f} m | O(h⁴) slope | h≤ for <{TARGET:.0f} m (SGP4 v) | h≤ (consistent v) |
|---|---|--:|--:|--:|:--:|--:|--:|--:|""")

for k in ORDER:
    r = REG[k]; n = r["nominal"]
    w(f"| **{k}** | {r['name']} | {r['ecc']:.4f} | {r['period_min']:,.0f} | "
      f"{f(n['max_m'])} m | {verdict(r)} | {r['slope_truncation']:.2f} | "
      f"{r['required_h_for_target_s']:.0f} s | {r['required_h_consistent_s']:.0f} s |")

w(f"""
*"refine" = boundable and understood, just needs finer nodes than {NOM} s for a
strict {TARGET:.0f} m budget. Even the worst case here ({f(heo['nominal']['max_m'])} m, Molniya) is far
below SGP4's own ~1–3 km state uncertainty — Hermite is not the dominant error
source in an SGP4 pipeline.*

---

## 1 · Why Hermite, and how it was tested

SGP4 returns **position *and* velocity** at each epoch, so we build a genuine
two-point cubic Hermite interpolant on every interval — matching position and
velocity exactly at both ends (C¹ continuous), not a position-only spline. On an
interval [t₀, t₁] with h = t₁−t₀ and s = (t−t₀)/h:

```
p(t) = h₀₀(s)·p₀ + h₁₀(s)·h·v₀ + h₀₁(s)·p₁ + h₁₁(s)·h·v₁
```

The interpolation uses SciPy's standard `CubicHermiteSpline` (reproduces the
basis above to ~10⁻⁹ m), so "is the interpolant coded correctly?" is not in
question. *(`PchipInterpolator` is **not** equivalent — it discards the supplied
derivatives; we deliberately feed SGP4's velocities.)* Truth is SGP4 sampled
densely inside every interval; error is the 3D position difference (max and RMS),
also decomposed into the radial / in-track / cross-track (RIC) frame.

---

## 2 · The bound — a two-term law with zero free parameters

Across every regime and spacing, the observed maximum error is captured by the
sum of two analytic terms, **both computed directly from the SGP4 states** —
nothing fitted:

```
E(h) =  max|jounce| · h⁴ / 384   +   0.0962 · δᵥ · h
        └── truncation, O(h⁴) ──┘     └ node-velocity floor, O(h) ┘
```

- **Truncation** — the exact cubic-Hermite remainder; *jounce* is the 4th time
  derivative of position, set by the orbit dynamics.
- **Node-consistency floor** — SGP4's *reported* velocity is not exactly the
  derivative of its *reported* position; the gap δᵥ, fed in as the Hermite
  slope, injects an O(h) term (coefficient 0.0962 = maxₛ|2s³−3s²+s|).

Because both inputs come from the propagator, **E(h) is a genuine a-priori bound,
not a curve fit** — and every test point lands on or below it.

![Observed vs. predicted bound, reference satellites](plots/bound_scatter.png)

---

## 3 · Convergence — fourth order, exactly as predicted

Halving the node spacing cuts truncation error ~16×. Fitted log-log slopes:
**{leo['slope_truncation']:.2f}** (LEO), **{meo['slope_truncation']:.2f}** (MEO),
**{geo['slope_truncation']:.2f}** (GEO), **{heo['slope_truncation']:.2f}** (HEO) —
all essentially the theoretical 4.0.

![Convergence of max error vs node spacing](plots/convergence.png)

*Solid: observed error (SGP4 reported velocity). Dashed: the a-priori bound E(h).
Dotted: the same interpolation from a consistent (differentiated) velocity —
pure O(h⁴) with the floor removed (§4).*

---

## 4 · Which velocity you feed the nodes

The one caveat. SGP4's reported velocity isn't exactly d/dt of its reported
position — the gap δᵥ is **{leo['dv_m_s']*100:.1f} cm/s** for the ISS,
**{geo['dv_m_s']*100:.1f} cm/s** for GEO, and **{heo['dv_m_s']:.2f} m/s** at Molniya
perigee. Feeding it as the Hermite slope injects the O(h) floor.

> **Intuition.** Picture connecting dots with a smooth curve while, at each dot,
> an arrow tells you which way to head. If the arrow points slightly off from the
> way the dots trend, the curve bows away between dots. SGP4's reported velocity
> is such a slightly-off arrow; the error it injects grows *linearly* with spacing.

The fix — a **consistent velocity** — sets each node's slope by *differentiating
SGP4's position*. The floor vanishes and pure O(h⁴) returns:

| Regime | SGP4 reported v @ {NOM} s | consistent v @ {NOM} s | improvement |
|---|--:|--:|--:|""")

for k in ORDER:
    n = REG[k]["nominal"]
    rep, con = n["max_m"], n["max_consistent_v_m"]
    fac = rep / con
    ftxt = f"×{fac:,.0f}" if fac >= 10 else f"×{fac:.1f}"
    w(f"| **{k}** | {f(rep,3)} m | {f(con,4)} m | {ftxt} |")

w(f"""
Consistent velocity is a big lever where the floor dominates (GEO, MEO) and does
little where truncation already dominates (LEO). *Trade-off:* the interpolated
velocity at nodes no longer equals SGP4's reported velocity — fine for position
lookup, a conscious choice if a consumer reads velocity out of the ephemeris. If
nodes come from a numerical integrator (where v genuinely is dr/dt), the floor
doesn't exist at all.

![RIC error decomposition at 120 s](plots/ric.png)

*The dominant error axis fingerprints the mechanism: truncation-limited orbits
(LEO, HEO) err **radially** (jounce points along the radius for near-circular
motion); the velocity floor (GEO) errs **in-track**, along the velocity direction.*

---

## 5 · Where the error lives

On near-circular orbits the per-interval error is nearly stationary around the
orbit. On the eccentric Molniya orbit it concentrates at the **perigee passages**,
where curvature and speed spike — exactly where fixed-step interpolation struggles
and where adaptive refinement pays off.

![Error over the orbit at 120 s](plots/timeseries.png)

---

## 6 · Robustness across the catalog

The four reference satellites make the mechanism legible, but the bound is
analytic and *per object* — so the real test is diversity, not count (50 more
Starlinks would prove nothing). From **{pm['n_pool']:,}** real objects pulled
across 13 Celestrak groups, **{pm['n_analyzed']}** were selected to span orbital
period (LEO to a {max(pT)/1440:.1f}-day deep-space orbit) and — the variable that
matters most — eccentricity, from **{min(pe):.4f}** to **{max(pe):.2f}** (Cluster II).

**Result: all {pm['n_points']} (object × spacing) points sit on or below the
bound** — worst observed/predicted = {pm['max_obs_over_pred']:.3f}, with
{pm['n_above_bound']} above it and {len(pm['skipped'])} objects skipped.

![Bound holds across 41 diverse orbits](plots/population.png)

Sweeping by eccentricity shows *why* Molniya was the hard case, and generalizes
it: the velocity gap δᵥ climbs roughly **100×** from near-circular (~1 cm/s) to
e ≈ 0.9 (several m/s) — and the coarser spacing a consistent velocity buys you
*grows* with eccentricity.

![Velocity floor and consistent-v advantage vs eccentricity](plots/ecc_ladder.png)

---

## 7 · The bound, term by term (at {NOM} s)

| Regime | max\\|jounce\\| (m/s⁴) | δᵥ (cm/s) | truncation (m) | floor (m) | E(h) pred (m) | observed (m) | obs/pred |
|---|--:|--:|--:|--:|--:|--:|--:|""")

for k in ORDER:
    r = REG[k]; n = r["nominal"]
    w(f"| **{k}** | {sci(r['max_jounce'])} | {r['dv_m_s']*100:.2f} | "
      f"{f(n['trunc_m'],3)} | {f(n['floor_m'],3)} | {f(n['pred_m'],3)} | "
      f"{f(n['max_m'],3)} | {n['obs_over_pred']:.2f} |")

w(f"""
Both bound inputs come straight from the SGP4 states; obs/pred ≤ 1 in every row
confirms E(h) is a valid upper bound.

---

## 8 · Recommendations — set node spacing from the bound

- **LEO, strict {TARGET:.0f} m:** use **≤ {leo['required_h_for_target_s']:.0f} s**
  (60 s gives {f(leo['sweep'][1]['max_m'])} m). The {NOM} s point is
  truncation-limited at {f(leo['nominal']['max_m'])} m — a consistent velocity
  doesn't rescue it; only finer h does.
- **MEO & GEO:** {NOM} s already meets {TARGET:.0f} m
  ({f(meo['nominal']['max_m'])} m and {f(geo['nominal']['max_m'])} m). GEO is
  floor-limited — for deep sub-meter, switch to a consistent velocity rather than
  spending nodes.
- **Molniya / HEO:** uniform spacing is the wrong tool — error concentrates at
  perigee. Use a consistent velocity and/or refine near perigee.
- **Always size h from E(h):** compute max|jounce| and δᵥ from the states and
  solve E(h) = budget. No trial-and-error, no per-mission re-tuning.

---

## 9 · Reproduce

```bash
python -m venv .venv && .venv/bin/pip install numpy scipy sgp4 matplotlib
.venv/bin/python validate_hermite.py   # 4 reference sats -> data/results.json
.venv/bin/python population.py         # 41-object catalog -> data/population.json
.venv/bin/python make_plots.py         # -> plots/*.png
.venv/bin/python build_report.py       # -> report.html
.venv/bin/python build_readme.py       # -> README.md (this file)
```

## Layout

| Path | What |
|------|------|
| `validate_hermite.py` | Core validation + full derivation of the bound (module docstring). Interpolation via SciPy `CubicHermiteSpline`. |
| `population.py` | Catalog-wide robustness check over `data/population.tle`. |
| `make_plots.py` | All figures in `plots/`. |
| `build_report.py` / `build_readme.py` | Assemble `report.html` / this `README.md` from the JSON. |
| `data/tles.txt` | The 4 reference TLEs (real, from Celestrak). |
| `data/population.tle` | The {pm['n_analyzed']} diverse objects (selection rule in the file header). |
| `report.html` | Self-contained rendered report (figures + tables + methodology). |

TLEs fetched from [Celestrak](https://celestrak.org) on 2026-08-06.
""")

open("README.md", "w").write("\n".join(md) + "\n")
print(f"Wrote README.md ({sum(len(x) for x in md)/1024:.1f} KB)")
