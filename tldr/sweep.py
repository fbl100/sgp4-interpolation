#!/usr/bin/env python3
"""
TL;DR calibration sweep: max Hermite position error vs node spacing (1 s -> 300 s).
==================================================================================

Standalone companion to the main study (does NOT modify or import any of the
existing analysis code paths except the well-tested SGP4/Hermite helpers). For
each of the four reference regimes it sweeps node spacing h across [1, 300] s and
records the max 3D position error of a piecewise cubic Hermite interpolant, for
TWO choices of the node-velocity fed to the interpolant:

  * "raw"  -- SGP4's reported velocity at each node (what the propagator hands you)
  * "diff" -- velocity from a simple 2-point central difference of SGP4 POSITION:
                  v(t) = ( r(t+dt) - r(t-dt) ) / (2*dt)
              i.e. a self-consistent derivative of the position trajectory.

Note on dt: a 2-point central difference subtracts two nearly-equal position
vectors, so it is roundoff-sensitive. SGP4's internal Julian-date arithmetic has
an effective time resolution of ~5e-5 s, so dt=1e-3 s occasionally lands two
samples across a rounding boundary and produces a catastrophic-cancellation spike
(seen as isolated 0.04-0.7 m outliers). dt=1e-2 s sits at the roundoff/truncation
sweet spot -- the residual velocity floor bottoms out (~1e-5 m of position error)
and the spikes drop ~10x to <0.005 m, negligible against any practical budget.
We therefore use dt=1e-2 s (same method, just above the precision floor).

The point is CALIBRATION: given an acceptable error budget, read off the largest
node spacing you may use -- per regime, and per velocity choice.

Method (per (regime, h, velocity)):
  Truth is densely-sampled SGP4 (same propagator that produced the nodes), so this
  isolates interpolation error with zero force-model mismatch. To find the max over
  a full orbit at any spacing without cost exploding as h->0, we measure a bounded
  set of real, contiguous-geometry intervals whose starts are spread across one
  orbital period (plus a cluster around perigee so eccentric orbits' error spike is
  never missed), and take the max. Validated against the main report's known points
  (LEO@120s~6.0 m, MEO~0.11 m, GEO~0.63 m, HEO~13.9 m for raw velocity).

Run:  python tldr/sweep.py   ->   writes tldr/sweep.json
"""
import os
import sys
import json
import math
import numpy as np
from sgp4.api import Satrec

# reuse only the vetted, side-effect-free helpers from the main module
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from validate_hermite import load_tles, propagate, period_seconds  # noqa: E402

# ---- sweep configuration --------------------------------------------------- #
H_MIN, H_MAX, N_H = 1.0, 300.0, 60          # node spacings (s), log-spaced
DT_CENTRAL = 1e-2                            # central-difference step (s); see docstring
SPI = 24                                     # interior truth samples per interval
N_STARTS = 400                              # interval-start phases spread over 1 period
N_PERIGEE = 120                             # extra start phases clustered at perigee
PERIGEE_FRAC = 0.06                         # +/- fraction of period around perigee

# reference points from the main report, for a sanity check on the method
REPORT_120 = {"LEO": 6.00, "MEO": 0.11, "GEO": 0.63, "HEO": 13.93}


# ---- Hermite basis (vectorised across intervals, matches CubicHermiteSpline) #
def hermite_eval(p0, v0, p1, v1, h, s):
    """p0..v1: [K,3] node states (km, km/s); h: scalar; s: [S] in [0,1].
    Returns [K,S,3] interpolated positions (km)."""
    s = s[None, :, None]                      # [1,S,1]
    h00 = 2*s**3 - 3*s**2 + 1
    h10 = s**3 - 2*s**2 + s
    h01 = -2*s**3 + 3*s**2
    h11 = s**3 - s**2
    p0 = p0[:, None, :]; v0 = v0[:, None, :]  # [K,1,3]
    p1 = p1[:, None, :]; v1 = v1[:, None, :]
    return h00*p0 + h10*(h*v0) + h01*p1 + h11*(h*v1)


def central_velocity(sat, t):
    """2-point central difference of SGP4 position: (r(t+dt)-r(t-dt))/(2*dt)."""
    t = np.asarray(t, float)
    n = len(t)
    tp = np.concatenate([t + DT_CENTRAL, t - DT_CENTRAL])
    r, _ = propagate(sat, tp)
    return (r[:n] - r[n:]) / (2.0 * DT_CENTRAL)   # km/s


def perigee_time(sat, per, n=4000):
    """Time offset (s) of minimum geocentric radius within one period."""
    t = np.linspace(0.0, per, n)
    r, _ = propagate(sat, t)
    return float(t[int(np.argmin(np.linalg.norm(r, axis=1)))])


def start_times(sat, per):
    """Interval-start phases over one period: uniform coverage + a perigee cluster
    (so an eccentric orbit's sharp perigee error is always sampled)."""
    uni = np.linspace(0.0, per, N_STARTS, endpoint=False)
    tp = perigee_time(sat, per)
    span = PERIGEE_FRAC * per
    peri = np.linspace(tp - span, tp + span, N_PERIGEE)
    ts = np.concatenate([uni, peri]) % per
    return np.unique(ts)


def max_error(sat, per, h, starts, use_reported_v):
    """Max 3D position error (m) over the sampled contiguous intervals at spacing h."""
    t0 = starts
    t1 = t0 + h
    # node states
    r0, v0_rep = propagate(sat, t0)
    r1, v1_rep = propagate(sat, t1)
    if use_reported_v:
        v0, v1 = v0_rep, v1_rep
    else:
        v0 = central_velocity(sat, t0)
        v1 = central_velocity(sat, t1)
    # interior truth samples (midpoint-shifted, exclude the nodes)
    s = (np.arange(1, SPI + 1) - 0.5) / SPI            # [S]
    p = hermite_eval(r0, v0, r1, v1, h, s)             # [K,S,3] km
    ts = (t0[:, None] + s[None, :] * h).ravel()        # [K*S]
    rt, _ = propagate(sat, ts)                          # truth [K*S,3]
    rt = rt.reshape(len(t0), SPI, 3)
    err = np.linalg.norm((p - rt) * 1000.0, axis=2)    # [K,S] meters
    return float(np.max(err))


def main():
    sats = load_tles("data/tles.txt")
    grid = np.geomspace(H_MIN, H_MAX, N_H)
    anchors = [30.0, 60.0, 120.0]           # operational spacings the report cites
    hs = np.unique(np.round(np.concatenate([grid, anchors]), 4))

    out = {"meta": {"h_min_s": H_MIN, "h_max_s": H_MAX, "n_h": N_H,
                    "dt_central_s": DT_CENTRAL, "spi": SPI,
                    "source": "Celestrak GP/TLE, fetched 2026-08-06",
                    "note": "raw = SGP4 reported velocity; diff = central-difference "
                            f"velocity of SGP4 position, dt={DT_CENTRAL:g} s"},
           "h_s": hs.tolist(), "regimes": []}

    for sd in sats:
        sat = Satrec.twoline2rv(sd["l1"], sd["l2"])
        regime = sd["regime"]
        per = period_seconds(sat)
        starts = start_times(sat, per)

        raw = [max_error(sat, per, float(h), starts, True) for h in hs]
        dif = [max_error(sat, per, float(h), starts, False) for h in hs]

        # sanity check vs the main report at exactly h=120 s
        i120 = int(np.argmin(np.abs(hs - 120.0)))
        chk = REPORT_120[regime]
        ratio = raw[i120] / chk if chk else float("nan")

        print(f"{regime:3s} {sd['name']:22s} T={per/60:7.1f}m e={sat.ecco:.4f} "
              f"| raw@120s={raw[i120]:8.3f} m (report {chk:.2f}, x{ratio:.2f}) "
              f"| diff@120s={dif[i120]:8.4f} m")

        out["regimes"].append({
            "regime": regime, "name": sd["name"], "norad": int(sat.satnum),
            "period_min": per / 60.0, "ecc": float(sat.ecco),
            "err_raw_m": raw, "err_diff_m": dif,
            "check_h_s": float(hs[i120]),
            "check_raw_m": raw[i120], "check_report_m": chk,
        })

    with open("tldr/sweep.json", "w") as f:
        json.dump(out, f, indent=2)
    print("\nWrote tldr/sweep.json")


if __name__ == "__main__":
    main()
