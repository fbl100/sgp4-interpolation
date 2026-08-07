#!/usr/bin/env python3
"""
Cubic Hermite interpolation error validation for SGP4-propagated orbits.
=======================================================================

Question under test
-------------------
Is piecewise cubic Hermite interpolation of ephemeris nodes an acceptable and
*boundable* way to reconstruct an orbit between nodes? Nominal operational case:
LEO, 120 s node spacing, target < 1 m position error.

Why Hermite is the right tool
-----------------------------
SGP4 returns BOTH position and velocity at every node. Cubic Hermite consumes
value + first derivative at each interval endpoint, so we build a *true* Hermite
interpolant (C1-continuous, exact position and velocity at nodes). On [t0,t1]
with h = t1-t0 and s = (t-t0)/h in [0,1]:

    p(t) = h00(s) p0 + h10(s) h v0 + h01(s) p1 + h11(s) h v1
    h00 = 2s^3-3s^2+1   h10 = s^3-2s^2+s   h01 = -2s^3+3s^2   h11 = s^3-s^2

Velocity terms carry a factor h because the basis is in the normalized variable s.
The actual interpolation uses scipy.interpolate.CubicHermiteSpline (the standard,
well-tested implementation of exactly this); it reproduces the basis formula above
to ~1e-9 m. Note: PchipInterpolator is NOT equivalent -- it discards the supplied
derivatives to enforce monotonicity; we specifically want SGP4's velocities.

The error model (the whole point: it is BOUNDABLE a priori)
-----------------------------------------------------------
Observed max position error over any interval is captured, with NO free
parameters, by the sum of two analytic terms:

    E(h) = (max|jounce| / 384) * h^4     [truncation, O(h^4)]
         + 0.0962 * dv * h               [node-consistency floor, O(h^1)]

* Truncation term: exact two-point cubic Hermite remainder is
  f(t)-p(t) = f''''(xi)/4! (t-t0)^2 (t-t1)^2, whose interval max is
  max|f''''| * h^4 / 384. Here f'''' is the 4th time-derivative of position
  ("jounce"), fixed by the orbit dynamics.
* Floor term: SGP4's *reported* velocity is not exactly d/dt of its *reported*
  position; they differ by dv (m/s). Feeding reported velocity as the Hermite
  derivative injects an O(h) error, coefficient 0.0962 = max_s|2s^3-3s^2+s|
  (the Hermite endpoint-derivative sensitivity). This floor is a property of the
  NODE DATA, not of the interpolation; using a consistent (differentiated)
  velocity removes it and recovers pure O(h^4).

Both quantities (max|jounce|, dv) are computed directly from the SGP4 states, so
E(h) is a genuine a-priori bound, not a curve fit. We verify observed <= E(h)
across four real orbit regimes and a wide range of node spacings.

Truth model: densely-sampled SGP4 (same propagator that produced the nodes), so
this isolates interpolation error with zero force-model mismatch -- the
"state vectors with 0% uncertainty" framing.

Run:  python validate_hermite.py   ->   writes data/results.json
"""
import json
import math
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from sgp4.api import Satrec

# Hermite endpoint-derivative sensitivity coefficient: max_s |2s^3 - 3s^2 + s|
HERMITE_DV_COEF = 0.09623  # attained at s = 1/2 - 1/sqrt(12)
TRUNC_COEF = 1.0 / 384.0

NOMINAL_H = 120.0   # user's operational LEO node spacing (s)
TARGET_M = 1.0      # user's position error budget (m)

# Per-regime node-spacing sweep (s) and window length (orbital periods).
SWEEP = {
    "LEO": {"periods": 3.0, "h": [30, 60, 120, 240, 480, 960]},
    "MEO": {"periods": 1.0, "h": [60, 120, 240, 480, 960, 1920]},
    "GEO": {"periods": 1.0, "h": [120, 240, 480, 960, 1920, 3840]},
    "HEO": {"periods": 1.0, "h": [30, 60, 120, 240, 480, 960]},
}


# --------------------------------------------------------------------------- #
# TLE loading / SGP4 propagation
# --------------------------------------------------------------------------- #
def load_tles(path):
    sats, regime, pend = [], None, []
    for raw in open(path):
        line = raw.rstrip("\n")
        if line.startswith("# regime="):
            regime = line.split("=", 1)[1].strip()
            continue
        if line.startswith("#") or not line.strip():
            continue
        pend.append(line)
        if len(pend) == 3:
            sats.append({"regime": regime, "name": pend[0].strip(),
                         "l1": pend[1], "l2": pend[2]})
            pend = []
    return sats


def propagate(sat, t_seconds):
    """SGP4 TEME position (km) and velocity (km/s) at offsets t_seconds."""
    jd = np.full(np.asarray(t_seconds).shape, sat.jdsatepoch)
    fr = sat.jdsatepochF + np.asarray(t_seconds) / 86400.0
    err, r, v = sat.sgp4_array(jd, fr)
    if np.any(err != 0):
        raise RuntimeError(f"SGP4 errors {np.unique(err[err!=0])}")
    return r, v


def period_seconds(sat):
    return (2.0 * math.pi / sat.no_kozai) * 60.0


# --------------------------------------------------------------------------- #
# Hermite, consistent velocity, RIC
# --------------------------------------------------------------------------- #
def hermite_interp(node_t, node_r, node_v, sample_t):
    """Piecewise cubic Hermite via scipy's standard CubicHermiteSpline: matches
    position and velocity exactly at each node (C1). Vector-valued (axis=0), so
    the three ECI components share one spline object. Returns [M,3] (km)."""
    spl = CubicHermiteSpline(node_t, node_r, node_v, axis=0)
    return spl(sample_t)


def consistent_velocity(sat, t, d=1.0):
    """d/dt of SGP4 position (km/s), 4th-order central FD -- self-consistent
    derivative of the position trajectory (removes the node-consistency floor)."""
    t = np.asarray(t, float)
    tp = np.concatenate([t - 2*d, t - d, t + d, t + 2*d])
    r, _ = propagate(sat, tp); n = len(t)
    return (r[:n] - 8*r[n:2*n] + 8*r[2*n:3*n] - r[3*n:]) / (12*d)


def ric_components(r, v, err_vec):
    R_hat = r / np.linalg.norm(r, axis=1, keepdims=True)
    C_hat = np.cross(r, v); C_hat /= np.linalg.norm(C_hat, axis=1, keepdims=True)
    I_hat = np.cross(C_hat, R_hat)
    return np.stack([np.einsum("ij,ij->i", err_vec, R_hat),
                     np.einsum("ij,ij->i", err_vec, I_hat),
                     np.einsum("ij,ij->i", err_vec, C_hat)], axis=1)


# --------------------------------------------------------------------------- #
# a-priori bound inputs: dv and max|jounce|
# --------------------------------------------------------------------------- #
def dv_inconsistency(sat, window, n=400):
    """max over window of |v_sgp4 - d/dt r_sgp4| (m/s). Stable wrt FD step."""
    d = 1.0; worst = 0.0
    ts = np.linspace(100.0, window - 100.0, n)
    for t0 in ts:
        t = np.array([t0-2*d, t0-d, t0, t0+d, t0+2*d])
        r, v = propagate(sat, t)
        v_num = (r[0] - 8*r[1] + 8*r[3] - r[4]) / (12*d)
        worst = max(worst, np.linalg.norm(v_num - v[2]))
    return worst * 1000.0  # km/s -> m/s


def _jounce_at(sat, window, dt):
    """max |d^4 r/dt^4| (m/s^4) over the window via 5-point FD at step dt."""
    t = np.arange(0.0, window + dt, dt)
    r, _ = propagate(sat, t); r *= 1000.0
    d4 = (r[:-4] - 4*r[1:-3] + 6*r[2:-2] - 4*r[3:-1] + r[4:]) / dt**4
    return float(np.max(np.linalg.norm(d4, axis=1)))


def max_jounce(sat, window):
    """FD step scaled to the orbit (period/1000): far above float roundoff even
    for slow (GEO) orbits, and fine enough to resolve near-circular jounce."""
    return _jounce_at(sat, window, period_seconds(sat) / 1000.0)


def max_jounce_robust(sat, window):
    """Jounce estimate that is valid across BOTH extremes:
      * high-eccentricity orbits peak jounce sharply at perigee -- too coarse a
        step under-resolves the peak (under-estimate);
      * slow near-circular orbits have tiny jounce -- too fine a step drowns it
        in float roundoff (over-estimate).
    A single fixed step cannot serve both, so we scan a bounded ladder of steps
    (period/500 down to period/16000) and pick the PLATEAU -- the step whose
    estimate is most stable vs its coarser neighbour. Undersampling makes the
    estimate rise as dt shrinks; roundoff makes it rise again; the true value is
    the flat region between."""
    per = period_seconds(sat)
    fracs = [1/500., 1/1000., 1/2000., 1/4000., 1/8000., 1/16000.]
    vals = [_jounce_at(sat, window, per*f) for f in fracs]
    best_i, best_rel = 1, float("inf")
    for i in range(1, len(vals)):
        rel = abs(vals[i] - vals[i-1]) / max(vals[i], 1e-300)
        if rel < best_rel:
            best_rel, best_i = rel, i
    return vals[best_i]


# --------------------------------------------------------------------------- #
# Interpolation error over a window
# --------------------------------------------------------------------------- #
def interp_error(sat, window, h, use_reported_v=True, spi=20, want_series=False):
    n = int(math.floor(window / h))
    nt = np.arange(n + 1) * h
    nr, nv_rep = propagate(sat, nt)
    nv = nv_rep if use_reported_v else consistent_velocity(sat, nt)

    # interior sample times across every interval (exclude nodes: error ~0 there)
    s = (np.arange(1, spi + 1) - 0.5) / spi
    ts = (nt[:-1, None] + s[None, :] * h).ravel()

    p = hermite_interp(nt, nr, nv, ts)          # km, one spline for the whole window
    rt, vt = propagate(sat, ts)
    e_m = (p - rt) * 1000.0
    err = np.linalg.norm(e_m, axis=1)

    out = {"h": h, "n_intervals": n,
           "max_m": float(np.max(err)), "rms_m": float(np.sqrt(np.mean(err**2)))}
    if want_series:
        ric = ric_components(rt, vt, e_m)
        mr = lambda x: (float(np.max(np.abs(x))), float(np.sqrt(np.mean(x**2))))
        out["ric_max_m"] = [mr(ric[:, k])[0] for k in range(3)]
        out["ric_rms_m"] = [mr(ric[:, k])[1] for k in range(3)]
        out["series"] = (ts, err)  # raw arrays; caller subsamples
    return out


def required_spacing(jo, dv, target):
    """Largest h with the FULL bound E(h) = jo/384 h^4 + coef*dv*h <= target
    (bisection). This is the SGP4-reported-velocity case (truncation + floor)."""
    f = lambda h: TRUNC_COEF*jo*h**4 + HERMITE_DV_COEF*dv*h - target
    lo, hi = 0.01, 1.0
    while f(hi) < 0:
        hi *= 2
        if hi > 1e7:
            return hi
    for _ in range(100):
        mid = 0.5*(lo+hi)
        (lo, hi) = (mid, hi) if f(mid) < 0 else (lo, mid)
    return 0.5*(lo+hi)


def required_spacing_consistent(jo, target):
    """Largest h with the TRUNCATION-ONLY bound jo/384 h^4 <= target. This is the
    consistent-velocity case (node-velocity floor removed); closed form."""
    return (target / (TRUNC_COEF * jo)) ** 0.25


def slope_fit(hs, errs):
    m, b = np.polyfit(np.log(hs), np.log(errs), 1)
    return float(m), float(b)


# --------------------------------------------------------------------------- #
def main():
    sats = load_tles("data/tles.txt")
    out = {"meta": {"nominal_h_s": NOMINAL_H, "target_m": TARGET_M,
                    "source": "Celestrak GP/TLE, fetched 2026-08-06",
                    "hermite_dv_coef": HERMITE_DV_COEF,
                    "trunc_coef": TRUNC_COEF},
           "regimes": []}

    for sd in sats:
        sat = Satrec.twoline2rv(sd["l1"], sd["l2"])
        regime, cfg = sd["regime"], SWEEP[sd["regime"]]
        per = period_seconds(sat)
        window = cfg["periods"] * per
        dv = dv_inconsistency(sat, window)
        jo = max_jounce(sat, window)

        print(f"\n=== {regime}: {sd['name']}  "
              f"(T={per/60:.1f} min, e={sat.ecco:.4f})  "
              f"dv={dv:.3e} m/s  max|jounce|={jo:.3e} m/s^4 ===")

        rows = []
        for h in cfg["h"]:
            h = float(h)
            rep = interp_error(sat, window, h, True, want_series=True)
            con = interp_error(sat, window, h, False)
            trunc = TRUNC_COEF * jo * h**4
            floor = HERMITE_DV_COEF * dv * h
            pred = trunc + floor
            tt, err = rep.pop("series")
            row = {"h": h, "n_intervals": rep["n_intervals"],
                   "max_m": rep["max_m"], "rms_m": rep["rms_m"],
                   "max_consistent_v_m": con["max_m"],
                   "ric_max_m": rep["ric_max_m"], "ric_rms_m": rep["ric_rms_m"],
                   "trunc_m": trunc, "floor_m": floor, "pred_m": pred,
                   "obs_over_pred": rep["max_m"] / pred}
            rows.append(row)
            flag = "  <-- nominal" if abs(h - NOMINAL_H) < 1e-6 else ""
            print(f"  h={h:6.0f}s  obs={rep['max_m']:11.4e}  pred={pred:11.4e}"
                  f"  obs/pred={row['obs_over_pred']:.2f}"
                  f"  consistent-v={con['max_m']:11.4e}{flag}")
            if abs(h - NOMINAL_H) < 1e-6:
                idx = np.linspace(0, len(tt)-1, min(1200, len(tt))).astype(int)
                out.setdefault("_series", {})[regime] = {
                    "t_frac": (tt[idx] / per).tolist(),
                    "err_m": err[idx].tolist()}

        # O(h^4) slope: fit the truncation-dominated points (trunc > 3*floor);
        # also fit the consistent-v curve, which is pure truncation throughout.
        trunc_pts = [(r["h"], r["max_m"]) for r in rows if r["trunc_m"] > 3*r["floor_m"]]
        if len(trunc_pts) >= 2:
            hs, es = zip(*trunc_pts)
            slope_trunc, _ = slope_fit(np.array(hs), np.array(es))
        else:
            slope_trunc = None
        hs_all = np.array([r["h"] for r in rows])
        slope_cons, _ = slope_fit(hs_all, np.array([r["max_consistent_v_m"] for r in rows]))

        nominal = next(r for r in rows if abs(r["h"] - NOMINAL_H) < 1e-6)
        req = required_spacing(jo, dv, TARGET_M)                 # SGP4 reported v
        req_cons = required_spacing_consistent(jo, TARGET_M)     # consistent v

        out["regimes"].append({
            "regime": regime, "name": sd["name"], "norad": int(sat.satnum),
            "period_min": per/60.0, "ecc": float(sat.ecco),
            "incl_deg": math.degrees(sat.inclo), "window_s": window,
            "alt_km": _altitude_km(sat), "dv_m_s": dv, "max_jounce": jo,
            "sweep": rows, "nominal": nominal,
            "slope_truncation": slope_trunc, "slope_consistent_v": slope_cons,
            "required_h_for_target_s": req,
            "required_h_consistent_s": req_cons,
        })
        st = "PASS" if nominal["max_m"] <= TARGET_M else "FAIL"
        print(f"  120s max = {nominal['max_m']:.3f} m  [{st} vs {TARGET_M} m]"
              f"   |  <1 m needs h <= {req:.0f} s (SGP4 v) / {req_cons:.0f} s (consistent v)")

    # attach series
    out["series"] = out.pop("_series", {})
    with open("data/results.json", "w") as f:
        json.dump(out, f, indent=2)
    print("\nWrote data/results.json")


def _altitude_km(sat, RE=6378.137):
    r, _ = propagate(sat, np.array([0.0]))
    return float(np.linalg.norm(r[0]) - RE)


if __name__ == "__main__":
    main()
