#!/usr/bin/env python3
"""Population robustness study.

The per-satellite bound E(h) = max|jounce|*h^4/384 + 0.0962*dv*h is analytic, not
fitted -- so robustness is a question of whether it holds across DIVERSE orbits,
not of raw sample count. This script:

  1. loads the frozen, curated population (data/population.tle -- 41 real objects
     chosen by a documented diversity rule; see that file's header),
  2. for each, checks observed max error <= predicted bound at several
     period-fraction spacings (period-fair across very different orbits),
  3. builds an eccentricity ladder (how the velocity floor and required node
     spacing scale with e),
  4. skips objects SGP4 cannot propagate cleanly (logged, never silently),

writing data/population.json for the plots/report.
"""
import math
import json
import numpy as np
from sgp4.api import Satrec

import validate_hermite as V

REFERENCE = {25544, 26407, 19548, 25847}   # the four in-depth reference sats
# period-fraction spacings: span the floor-dominated (fine) and truncation-
# dominated (coarse) regimes for every object, comparably.
FRACS = [1/500.0, 1/125.0, 1/30.0]
TARGET = V.TARGET_M


# Source-pull size (documented provenance): the frozen population below was
# selected from this many unique objects across 13 Celestrak groups + a few
# named eccentric satellites, pulled 2026-08-06. The raw dump is not retained --
# see data/population.tle's header for the deterministic selection rule.
N_POOL_SOURCE = 985


# --------------------------------------------------------------------------- #
def load_frozen(path="data/population.tle"):
    """Load the frozen, curated population (name/line1/line2 triples). These 41
    objects were chosen by a deterministic diversity rule (documented in the
    file header): every object with e>0.1, plus near-circular objects sampled
    evenly by period band. Frozen so results reproduce without a live re-pull."""
    lines = [l.rstrip("\n") for l in open(path)
             if l.strip() and not l.startswith("#")]
    objs = []
    for i in range(0, len(lines) - 2, 3):
        name, l1, l2 = lines[i].strip(), lines[i+1], lines[i+2]
        if not (l1.startswith("1 ") and l2.startswith("2 ")):
            continue
        sat = Satrec.twoline2rv(l1, l2)
        per = (2*math.pi/sat.no_kozai)*60/60.0 if sat.no_kozai > 0 else 0.0
        objs.append({"num": sat.satnum, "name": name, "sat": sat,
                     "ecc": float(sat.ecco), "per_min": per})
    return objs


def classify(ecc, per_min):
    if ecc > 0.25:
        return "HEO"
    if per_min < 200:
        return "LEO"
    if per_min < 800:
        return "MEO"
    return "GEO"


# --------------------------------------------------------------------------- #
def analyze(o):
    """Return a result dict for one object, or None (with reason) if SGP4 fails."""
    sat = o["sat"]
    per = V.period_seconds(sat)
    window = per  # one orbital period
    try:
        dv = V.dv_inconsistency(sat, window)
        jo = V.max_jounce_robust(sat, window)
        samples = []
        for fr in FRACS:
            h = fr * per
            res = V.interp_error(sat, window, h, use_reported_v=True)
            trunc = V.TRUNC_COEF * jo * h**4
            floor = V.HERMITE_DV_COEF * dv * h
            pred = trunc + floor
            samples.append({"h": h, "obs": res["max_m"], "pred": pred})
        req_sgp4 = V.required_spacing(jo, dv, TARGET)
        req_cons = V.required_spacing_consistent(jo, TARGET)
    except Exception as e:                       # SGP4 error code, decayed, etc.
        return None, f"{o['name']} ({o['num']}): {type(e).__name__} {e}"
    return {
        "num": o["num"], "name": o["name"],
        "regime": classify(o["ecc"], o["per_min"]),
        "ecc": o["ecc"], "period_min": o["per_min"],
        "reference": o["num"] in REFERENCE,
        "dv_m_s": dv, "max_jounce": jo,
        "req_sgp4_s": req_sgp4, "req_cons_s": req_cons,
        "samples": samples,
    }, None


def main():
    sel = load_frozen()
    print(f"frozen population: {len(sel)} objects "
          f"(selected from {N_POOL_SOURCE} source objects)")

    results, skipped = [], []
    for o in sorted(sel, key=lambda o: o["ecc"]):
        r, err = analyze(o)
        if r is None:
            skipped.append(err); print("  SKIP", err); continue
        results.append(r)
        worst = max(s["obs"]/s["pred"] for s in r["samples"] if s["pred"] > 0)
        print(f"  e={r['ecc']:.4f} T={r['period_min']:6.0f}m {r['regime']:3s} "
              f"{r['name'][:22]:22s} worst obs/pred={worst:.2f}")

    # global check
    all_pts = [(s["pred"], s["obs"]) for r in results for s in r["samples"]]
    ratios = [o/p for p, o in all_pts if p > 0]
    n_over = sum(1 for x in ratios if x > 1.0)
    print(f"\n{len(all_pts)} (object x spacing) points; "
          f"max obs/pred = {max(ratios):.3f}; points above bound = {n_over}")

    # eccentricity ladder: nearest analyzed object to each target e
    targets = [0.0007, 0.012, 0.05, 0.10, 0.17, 0.30, 0.47, 0.65, 0.78, 0.87, 0.91]
    ladder, used = [], set()
    for te in targets:
        cand = min((r for r in results if r["num"] not in used),
                   key=lambda r: abs(r["ecc"] - te), default=None)
        if cand and abs(cand["ecc"] - te) < 0.15:
            ladder.append(cand["num"]); used.add(cand["num"])

    out = {
        "meta": {"target_m": TARGET, "n_pool": N_POOL_SOURCE, "n_analyzed": len(results),
                 "n_points": len(all_pts), "n_above_bound": n_over,
                 "max_obs_over_pred": max(ratios), "fracs": FRACS,
                 "skipped": skipped},
        "objects": results,
        "ladder_norads": ladder,
    }
    with open("data/population.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote data/population.json  (ladder: {len(ladder)} rungs, "
          f"{len(skipped)} skipped)")


if __name__ == "__main__":
    main()
