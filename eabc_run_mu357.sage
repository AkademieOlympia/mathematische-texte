# eabc_run_mu357.sage
# SageMath 10.5+
#
# Full first geometric EABC Omega=3 run for (p,q,r)=(3,5,7).
#
# Requires these files in the same directory:
#   eabc_rehm_fiber_generator.sage
#   eabc_s6_9_octavian_backend.sage
# Optionally:
#   eabc_rehm_modp_audit.sage
#
# IMPORTANT:
# The G2 deficit is computed after equivariant projection to Im(O):
#
#   u_pi = Im(pi), u_rho = Im(rho), u_sigma = Im(sigma)
#
# and
#
#   t = ||[u_pi,u_rho,u_sigma]||^2 /
#       (4 ||u_pi||^2 ||u_rho||^2 ||u_sigma||^2).
#
# This is NOT the same as dividing by 4*N(pi)*N(rho)*N(sigma)
# when the Rehm representatives have nonzero real part.

from sage.all import *
from collections import Counter
import os
import json
import time

BASE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else "."

load(os.path.join(BASE, "eabc_rehm_fiber_generator.sage"))
load(os.path.join(BASE, "eabc_s6_9_octavian_backend.sage"))

def initialize_rehm_357(run_modp_audit=False):
    """
    Generate and audit the Rehm fibers P(3), P(5), P(7).

    If run_modp_audit=True, also loads the independent finite-field
    kernel audit.  Pairwise class-intersection checking is enabled for
    p=3 and p=5; for p=7 it is skipped by default because O(|P(7)|^2)
    subspace intersections are a separate cost.
    """
    print("=== Initialisiere Rehm-Fasern P(3), P(5), P(7) ===")

    P3 = generate_rehm_P(3)
    P5 = generate_rehm_P(5)
    P7 = generate_rehm_P(7)

    assert len(P3) == 3**3 + 1 == 28
    assert len(P5) == 5**3 + 1 == 126
    assert len(P7) == 7**3 + 1 == 344

    # exact norm checks
    assert all(norm2(x) == 3 for x in P3)
    assert all(norm2(x) == 5 for x in P5)
    assert all(norm2(x) == 7 for x in P7)

    print("PASS: |P(3)| =", len(P3))
    print("PASS: |P(5)| =", len(P5))
    print("PASS: |P(7)| =", len(P7))

    if run_modp_audit:
        load(os.path.join(BASE, "eabc_rehm_modp_audit.sage"))

        print("\n=== Backend-B Mod-p Audit ===")
        print(audit_rehm_mod_p(3, P3, pairwise=True))
        print("Direktzählung |N_3| =", direct_count_Np_small(3))

        print(audit_rehm_mod_p(5, P5, pairwise=True))
        print("Direktzählung |N_5| =", direct_count_Np_small(5))

        # For p=7: kernel/class audit without all pairwise intersections.
        print(audit_rehm_mod_p(7, P7, pairwise=False))

    return P3, P5, P7


def run_geometric_audit_357(P3, P5, P7, progress_every=100000):
    """
    Exact streaming computation of mu_357.

    Memory use is O(|supp(mu)|), not O(1,213,632):
    no list of all t-values is stored.
    """
    total_expected = 28 * 126 * 344
    total = len(P3) * len(P5) * len(P7)
    assert total == total_expected

    print("\n=== Starte geometrischen Omega=3-Lauf (3,5,7) ===")
    print("Konfigurationen:", total)

    hist = Counter()
    assoc_count = ZZ(0)
    sum_t = QQ(0)
    sum_t2 = QQ(0)

    # Optional audit data for the richer (detG,phi^2,t) spectrum
    joint_hist = Counter()

    done = 0
    t0 = time.time()

    for pi in P3:
        for rho in P5:
            for sigma in P7:
                # Uses exact rational G2 geometry after imaginary projection.
                g = triple_geometry(pi, rho, sigma)
                t = QQ(g["t"])
                d = QQ(g["detG_norm"])
                ph2 = QQ(g["phi2_norm"])

                # identity already checked inside triple_geometry, repeat as
                # a hard guard at the aggregate layer
                assert t == d - ph2

                hist[t] += 1
                joint_hist[(d, ph2, t)] += 1

                sum_t += t
                sum_t2 += t*t

                if t == 0:
                    assoc_count += 1

                done += 1
                if progress_every and done % progress_every == 0:
                    elapsed = time.time() - t0
                    print(f"  {done}/{total} verarbeitet; "
                          f"|supp(mu)|={len(hist)}; elapsed={elapsed:.1f}s")

    assert done == total
    assert sum(hist.values()) == total

    M1 = sum_t / total
    M2 = sum_t2 / total
    variance = M2 - M1*M1
    P_assoc = QQ(assoc_count) / total

    mu = {t: QQ(n)/total for t,n in sorted(hist.items(), key=lambda kv: kv[0])}

    # exact probability normalization
    assert sum(mu.values()) == 1

    result = {
        "pqr": (3,5,7),
        "color_type": color_multiset(3,5,7),
        "v4_total": v4_total(3,5,7),
        "total": total,
        "support_size": len(hist),
        "joint_support_size": len(joint_hist),
        "assoc_count": assoc_count,
        "P_assoc": P_assoc,
        "M1": M1,
        "M2": M2,
        "variance": variance,
        "histogram": hist,
        "mu": mu,
        "joint_histogram": joint_hist,
    }

    print("\n=== Ergebnisse mu_357 ===")
    print("V4-Farbtyp             :", result["color_type"])
    print("V4-Gesamtgrad          :", result["v4_total"])
    print("Total                  :", total)
    print("|supp(mu_357)|         :", len(hist))
    print("|supp(detG,phi^2,t)|   :", len(joint_hist))
    print("assoziative Tripel     :", assoc_count)
    print("P_assoc exakt          :", P_assoc)
    print("P_assoc numerisch      :", RDF(P_assoc))
    print("M1 exakt               :", M1)
    print("M1 numerisch           :", RDF(M1))
    print("M2 exakt               :", M2)
    print("Var exakt              :", variance)
    print("Var numerisch          :", RDF(variance))

    return result


def save_mu357(result, basename="mu_357_exact"):
    """
    Save compact exact summaries.
    Rational values are serialized as strings so no precision is lost.
    """
    txt_path = os.path.join(BASE, basename + ".txt")
    json_path = os.path.join(BASE, basename + ".json")

    with open(txt_path, "w") as f:
        f.write("EABC exact mu_357 audit\n")
        f.write("=======================\n")
        f.write(f"pqr = {result['pqr']}\n")
        f.write(f"color_type = {result['color_type']}\n")
        f.write(f"v4_total = {result['v4_total']}\n")
        f.write(f"total = {result['total']}\n")
        f.write(f"support_size = {result['support_size']}\n")
        f.write(f"joint_support_size = {result['joint_support_size']}\n")
        f.write(f"assoc_count = {result['assoc_count']}\n")
        f.write(f"P_assoc = {result['P_assoc']}\n")
        f.write(f"M1 = {result['M1']}\n")
        f.write(f"M2 = {result['M2']}\n")
        f.write(f"variance = {result['variance']}\n")
        f.write("\nmu_357:\n")
        for t, prob in result["mu"].items():
            f.write(f"{t}\t{prob}\tcount={result['histogram'][t]}\n")

    payload = {
        "pqr": list(result["pqr"]),
        "color_type": list(result["color_type"]),
        "v4_total": list(result["v4_total"]),
        "total": int(result["total"]),
        "support_size": int(result["support_size"]),
        "joint_support_size": int(result["joint_support_size"]),
        "assoc_count": int(result["assoc_count"]),
        "P_assoc": str(result["P_assoc"]),
        "M1": str(result["M1"]),
        "M2": str(result["M2"]),
        "variance": str(result["variance"]),
        "mu": [
            {
                "t": str(t),
                "probability": str(prob),
                "count": int(result["histogram"][t]),
            }
            for t,prob in result["mu"].items()
        ],
    }

    with open(json_path, "w") as f:
        json.dump(payload, f, indent=2)

    print("Gespeichert:", txt_path)
    print("Gespeichert:", json_path)
    return txt_path, json_path


def main(run_modp_audit=False):
    # Geometry regression test first
    print("=== Fano-Geometrie-Selbsttest ===")
    print(self_test())

    P3, P5, P7 = initialize_rehm_357(run_modp_audit=run_modp_audit)
    result = run_geometric_audit_357(P3, P5, P7)
    save_mu357(result)
    return result


# Typical first production run:
#
#   sage eabc_run_mu357.sage
#
# If you want the expensive independent p=3,p=5 mod-p closeout first:
#
#   sage -c 'load("eabc_run_mu357.sage"); main(run_modp_audit=True)'
#
if __name__ == "__main__":
    main(run_modp_audit=False)
