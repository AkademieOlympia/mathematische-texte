import math
from collections import defaultdict
import numpy as np

# --------------------------------------------------
# 1. Grundfunktionen: Primfaktorzerlegung / Signatur
# --------------------------------------------------

def factorize(n: int) -> dict[int, int]:
    d = 2
    out = {}
    x = n
    while d * d <= x:
        while x % d == 0:
            out[d] = out.get(d, 0) + 1
            x //= d
        d += 1
    if x > 1:
        out[x] = out.get(x, 0) + 1
    return out

def smooth_part(n: int, basis=(2, 3, 5)) -> int:
    fac = factorize(n)
    s = 1
    for p in basis:
        if p in fac:
            s *= p ** fac[p]
    return s

def rest_kernel(n: int, basis=(2, 3, 5)) -> int:
    return n // smooth_part(n, basis)

def smooth_vector(n: int) -> tuple[int, int, int]:
    fac = factorize(n)
    return (fac.get(2, 0), fac.get(3, 0), fac.get(5, 0))

def family_of_prime(p: int) -> str | None:
    if p in (2, 3):
        return None
    r = p % 12
    if r == 1:
        return "E"
    if r == 5:
        return "A"
    if r == 7:
        return "B"
    if r == 11:
        return "C"
    return None

def rho_profile(n: int, basis=(2, 3, 5)) -> dict[str, float]:
    k = rest_kernel(n, basis)
    fac = factorize(k)
    weights = {"E": 0.0, "A": 0.0, "B": 0.0, "C": 0.0}
    total = 0.0
    for p, v in fac.items():
        fam = family_of_prime(p)
        if fam is None:
            continue
        w = v * math.log(p)
        weights[fam] += w
        total += w
    if total == 0:
        return weights
    return {k: v / total for k, v in weights.items()}

def c_fam(n: int, basis=(2, 3, 5)) -> float:
    rho = rho_profile(n, basis)
    s = 0.0
    for val in rho.values():
        if val > 0:
            s -= val * math.log(val)
    return s

def epsilon_B(n: int, basis=(2, 3, 5)) -> float:
    if n <= 1:
        return 0.0
    m = rest_kernel(n, basis)
    return math.log(m) / math.log(n)

def omega_sig(n: int, basis=(2, 3, 5)) -> float:
    # erste einfache Signaturgewichtung
    return epsilon_B(n, basis) * c_fam(n, basis)


# ============================================================
# 1. 5-Regime / 23-235-Umschaltung
# ============================================================

def v_p(n: int, p: int) -> int:
    c = 0
    x = n
    while x % p == 0 and x > 0:
        x //= p
        c += 1
    return c


def v5(n: int) -> int:
    return v_p(n, 5)


def m23(n: int) -> int:
    fac = factorize(n)
    return n // (2 ** fac.get(2, 0) * 3 ** fac.get(3, 0))


def m235(n: int) -> int:
    fac = factorize(n)
    return n // (2 ** fac.get(2, 0) * 3 ** fac.get(3, 0) * 5 ** fac.get(5, 0))


def delta_35(n: int) -> float:
    """
    Misst die 3/5-Umschaltung über die Differenz
    log m23 - log m235 = v5(n) * log 5.
    """
    return math.log(m23(n)) - math.log(m235(n))


def v5_bucket(n: int) -> str:
    k = v5(n)
    if k == 0:
        return "v5=0"
    if k == 1:
        return "v5=1"
    return "v5>=2"


# ============================================================
# 1a. Hilfsgrößen für Signaturdistanz
# ============================================================

def eps_sig(n: int, basis=(2, 3, 5)) -> float:
    return epsilon_B(n, basis)


def rho_vec(n: int, basis=(2, 3, 5)) -> np.ndarray:
    rho = rho_profile(n, basis)
    return np.array([rho["E"], rho["A"], rho["B"], rho["C"]], dtype=float)


def kappa_vec(n: int, basis=(2, 3, 5)) -> np.ndarray:
    """
    Vorläufige Minimalversion:
    normierte Loggewichte der Kleinträger 5, 7, 11, 13 im Restkern.
    """
    k = rest_kernel(n, basis)
    fac = factorize(k)
    ks = []
    total = 0.0
    for p in [5, 7, 11, 13]:
        w = fac.get(p, 0) * math.log(p)
        ks.append(w)
        total += w
    if total == 0:
        return np.zeros(4, dtype=float)
    return np.array([x / total for x in ks], dtype=float)


def d_smooth(m: int, z: int) -> float:
    sv_m = np.array(smooth_vector(m), dtype=float)
    sv_z = np.array(smooth_vector(z), dtype=float)
    return float(np.sum(np.abs(sv_m - sv_z)))


def d_rho(m: int, z: int, basis=(2, 3, 5)) -> float:
    return float(np.linalg.norm(rho_vec(m, basis) - rho_vec(z, basis), ord=2))


def d_kappa(m: int, z: int, basis=(2, 3, 5)) -> float:
    return float(np.linalg.norm(kappa_vec(m, basis) - kappa_vec(z, basis), ord=2))


def d_eps(m: int, z: int, basis=(2, 3, 5)) -> float:
    return abs(eps_sig(m, basis) - eps_sig(z, basis))


def d_log(m: int, z: int) -> float:
    return abs(math.log(m) - math.log(z))


# ============================================================
# 2. Kostenfunktion für signaturgeometrische Flankenwahl
# ============================================================

def flank_cost(
    m: int,
    z: int,
    basis=(2, 3, 5),
    alpha=1.0,
    beta=1.0,
    gamma=1.0,
    delta=1.0,
    eta=0.5,
    require_same_smooth=False,
) -> float:
    """
    Kostenfunktion C(m,z).
    Kleinere Werte = bessere signaturgeometrische Flanke.
    """
    if require_same_smooth and smooth_vector(m) != smooth_vector(z):
        return float("inf")

    return (
        alpha * d_log(m, z)
        + beta * d_smooth(m, z)
        + gamma * d_rho(m, z, basis)
        + delta * d_kappa(m, z, basis)
        + eta * d_eps(m, z, basis)
    )


# --------------------------------------------------
# 2. Flankenfinder
# --------------------------------------------------

def family_of_number(n: int, basis=(2, 3, 5)) -> set[str]:
    """Familien, die im Restkern von n aktiv sind."""
    k = rest_kernel(n, basis)
    fac = factorize(k)
    fams = set()
    for p in fac:
        fam = family_of_prime(p)
        if fam is not None:
            fams.add(fam)
    return fams


# ============================================================
# 3. Familienfilter
# ============================================================

def family_support(n: int, basis=(2, 3, 5)) -> set[str]:
    return family_of_number(n, basis)


def matches_family_exact(n: int, fam: str, basis=(2, 3, 5)) -> bool:
    return family_support(n, basis) == {fam}


def matches_family_contains(n: int, fam: str, basis=(2, 3, 5)) -> bool:
    return fam in family_support(n, basis)


# ============================================================
# 4. Signaturgeometrische Flankenwahl
# ============================================================

def best_family_below(
    m: int,
    fam: str,
    basis=(2, 3, 5),
    max_search=500,
    match_mode="exact",
    require_same_smooth=False,
    alpha=1.0,
    beta=1.0,
    gamma=1.0,
    delta=1.0,
    eta=0.5,
):
    best_x = None
    best_c = float("inf")

    for x in range(m - 1, max(1, m - max_search), -1):
        ok = (
            matches_family_exact(x, fam, basis)
            if match_mode == "exact"
            else matches_family_contains(x, fam, basis)
        )
        if not ok:
            continue

        c = flank_cost(
            m,
            x,
            basis=basis,
            alpha=alpha,
            beta=beta,
            gamma=gamma,
            delta=delta,
            eta=eta,
            require_same_smooth=require_same_smooth,
        )
        if c < best_c:
            best_c = c
            best_x = x

    return best_x, best_c


def best_family_above(
    m: int,
    fam: str,
    basis=(2, 3, 5),
    max_search=500,
    match_mode="exact",
    require_same_smooth=False,
    alpha=1.0,
    beta=1.0,
    gamma=1.0,
    delta=1.0,
    eta=0.5,
):
    best_x = None
    best_c = float("inf")

    for x in range(m + 1, m + max_search):
        ok = (
            matches_family_exact(x, fam, basis)
            if match_mode == "exact"
            else matches_family_contains(x, fam, basis)
        )
        if not ok:
            continue

        c = flank_cost(
            m,
            x,
            basis=basis,
            alpha=alpha,
            beta=beta,
            gamma=gamma,
            delta=delta,
            eta=eta,
            require_same_smooth=require_same_smooth,
        )
        if c < best_c:
            best_c = c
            best_x = x

    return best_x, best_c


def flank_pair_signature(
    m: int,
    left_fam: str,
    right_fam: str,
    basis=(2, 3, 5),
    max_search=500,
    match_mode="exact",
    require_same_smooth=False,
    alpha=1.0,
    beta=1.0,
    gamma=1.0,
    delta=1.0,
    eta=0.5,
):
    left, c_left = best_family_below(
        m,
        left_fam,
        basis=basis,
        max_search=max_search,
        match_mode=match_mode,
        require_same_smooth=require_same_smooth,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        delta=delta,
        eta=eta,
    )
    right, c_right = best_family_above(
        m,
        right_fam,
        basis=basis,
        max_search=max_search,
        match_mode=match_mode,
        require_same_smooth=require_same_smooth,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        delta=delta,
        eta=eta,
    )

    if left is None or right is None:
        return None

    return {
        "left": left,
        "right": right,
        "cost_left": c_left,
        "cost_right": c_right,
        "cost_total": c_left + c_right,
    }


def nearest_family_below(
    m: int, fam: str, basis=(2, 3, 5), max_search: int = 500
) -> int | None:
    for x in range(m - 1, max(1, m - max_search), -1):
        fams = family_of_number(x, basis)
        if fams == {fam}:
            return x
    return None


def nearest_family_above(
    m: int, fam: str, basis=(2, 3, 5), max_search: int = 500
) -> int | None:
    for x in range(m + 1, m + max_search):
        fams = family_of_number(x, basis)
        if fams == {fam}:
            return x
    return None


def flank_pair(
    m: int, left_fam: str, right_fam: str, basis=(2, 3, 5), max_search: int = 500
):
    left = nearest_family_below(m, left_fam, basis, max_search)
    right = nearest_family_above(m, right_fam, basis, max_search)
    if left is None or right is None:
        return None
    return left, right


# --------------------------------------------------
# 3. Produktdaten erzeugen
# --------------------------------------------------

def wallis_factor(m: int, left: int, right: int) -> float:
    return (m * m) / (left * right)


# ============================================================
# Hilfsfilter
# ============================================================

def center_passes_filter(
    m: int,
    basis=(2, 3, 5),
    omega_min: float = 0.0,
    require_s000: bool = False,
) -> bool:
    """
    Filter für Zentren m:
    - omega_sig(m) > omega_min
    - optional glatte Bühne S=(0,0,0)
    """
    om = omega_sig(m, basis)
    if om <= omega_min:
        return False
    if require_s000 and smooth_vector(m) != (0, 0, 0):
        return False
    return True


def product_data_filtered(
    N: int,
    mode="BC",
    basis=(2, 3, 5),
    weighted=False,
    omega_min: float = 0.0,
    require_s000: bool = False,
    max_search: int = 500,
):
    """
    mode: 'BC' oder 'AE'
    weighted=False: logw = log(m^2/(left*right))
    weighted=True : logw += log(1+omega_sig(m))
    """
    if mode == "BC":
        left_fam, right_fam = "B", "C"
    elif mode == "AE":
        left_fam, right_fam = "A", "E"
    else:
        raise ValueError("mode muss 'BC' oder 'AE' sein")

    rows = []
    for m in range(10, N + 1):
        if not center_passes_filter(
            m,
            basis=basis,
            omega_min=omega_min,
            require_s000=require_s000,
        ):
            continue

        pair = flank_pair(
            m,
            left_fam,
            right_fam,
            basis=basis,
            max_search=max_search,
        )
        if pair is None:
            continue

        left, right = pair
        w = wallis_factor(m, left, right)
        logw = math.log(w)
        om = omega_sig(m, basis)

        if weighted:
            logw += math.log(1.0 + om)

        rows.append(
            {
                "m": m,
                "left": left,
                "right": right,
                "logm": math.log(m),
                "w": w,
                "logw": logw,
                "omega": om,
                "smooth_vec": smooth_vector(m),
                "families_m": sorted(list(family_of_number(m, basis))),
            }
        )

    return rows


# ============================================================
# 2. Produktdaten mit 5-Regime-Markierung
# ============================================================

def product_data_regime(
    N: int,
    mode="BC",
    basis=(2, 3, 5),
    weighted=False,
    omega_min: float = 0.0,
    require_s000: bool = False,
    max_search: int = 500,
):
    """
    Verwendet bewusst die alte lokale familienbasierte Flankenwahl.
    """
    rows = product_data_filtered(
        N=N,
        mode=mode,
        basis=basis,
        weighted=weighted,
        omega_min=omega_min,
        require_s000=require_s000,
        max_search=max_search,
    )

    out = []
    for r in rows:
        m = r["m"]
        left = r["left"]
        right = r["right"]

        r2 = dict(r)
        r2["v5_center"] = v5(m)
        r2["v5_left"] = v5(left)
        r2["v5_right"] = v5(right)
        r2["bucket_center"] = v5_bucket(m)
        r2["bucket_triplet"] = f"{v5_bucket(left)}|{v5_bucket(m)}|{v5_bucket(right)}"
        r2["delta35_center"] = delta_35(m)
        r2["delta35_left"] = delta_35(left)
        r2["delta35_right"] = delta_35(right)
        out.append(r2)

    return out


# ============================================================
# 1. Vereinfachte 5-Profilklassen
# ============================================================

def bucket_v5_small(n: int) -> str:
    k = v5(n)
    if k == 0:
        return "0"
    if k == 1:
        return "1"
    return "2+"


def triplet_profile_small(row) -> str:
    """
    Kompaktes Profil der Form left|center|right mit 0/1/2+.
    """
    return f"{bucket_v5_small(row['left'])}|{bucket_v5_small(row['m'])}|{bucket_v5_small(row['right'])}"


def flank_5_pattern(row) -> str:
    """
    Noch groebere Klassifikation:
    - none: keine Flanke traegt 5
    - left: nur linke Flanke traegt 5
    - right: nur rechte Flanke traegt 5
    - both: beide Flanken tragen 5
    """
    l = v5(row["left"]) > 0
    r = v5(row["right"]) > 0
    if not l and not r:
        return "none"
    if l and not r:
        return "left"
    if not l and r:
        return "right"
    return "both"


def center_flank_regime(row) -> str:
    """
    Kombiniert Zentrum und Flanken grob.
    Fuer den aktuellen Lauf ist das Zentrum meist 0,
    aber wir lassen es allgemein.
    """
    c = bucket_v5_small(row["m"])
    f = flank_5_pattern(row)
    return f"center={c},flanks={f}"


# ============================================================
# 3. Gruppierung nach 5-Regime
# ============================================================

def split_by_center_bucket(rows):
    groups = defaultdict(list)
    for r in rows:
        groups[r["bucket_center"]].append(r)
    return groups


def split_by_triplet_bucket(rows):
    groups = defaultdict(list)
    for r in rows:
        groups[r["bucket_triplet"]].append(r)
    return groups


def group_by_full_triplet_profile(rows):
    groups = defaultdict(list)
    for r in rows:
        key = f"{bucket_v5_small(r['left'])}|{bucket_v5_small(r['m'])}|{bucket_v5_small(r['right'])}"
        groups[key].append(r)
    return groups


def group_by_flank_pattern(rows):
    groups = defaultdict(list)
    for r in rows:
        key = center_flank_regime(r)
        groups[key].append(r)
    return groups


# ============================================================
# 5. Neue Wallis-Daten mit signaturgeometrischen Flanken
# ============================================================

def product_data_signature(
    N: int,
    mode="BC",
    basis=(2, 3, 5),
    weighted=False,
    omega_min=0.0,
    require_center_s000=False,
    flank_same_smooth=False,
    match_mode="exact",
    max_search=500,
    alpha=1.0,
    beta=1.0,
    gamma=1.0,
    delta=1.0,
    eta=0.5,
):
    """
    mode: BC oder AE
    weighted: logw += log(1+omega_sig(m))
    require_center_s000: Zentrum m muss S=(0,0,0) haben
    flank_same_smooth: Flanken müssen gleiche glatte Bühne wie Zentrum haben
    match_mode: 'exact' oder 'contains'
    """
    if mode == "BC":
        left_fam, right_fam = "B", "C"
    elif mode == "AE":
        left_fam, right_fam = "A", "E"
    else:
        raise ValueError("mode muss 'BC' oder 'AE' sein")

    rows = []

    for m in range(10, N + 1):
        om = omega_sig(m, basis)
        if om <= omega_min:
            continue
        if require_center_s000 and smooth_vector(m) != (0, 0, 0):
            continue

        pair = flank_pair_signature(
            m,
            left_fam,
            right_fam,
            basis=basis,
            max_search=max_search,
            match_mode=match_mode,
            require_same_smooth=flank_same_smooth,
            alpha=alpha,
            beta=beta,
            gamma=gamma,
            delta=delta,
            eta=eta,
        )
        if pair is None:
            continue

        left = pair["left"]
        right = pair["right"]

        w = wallis_factor(m, left, right)
        logw = math.log(w)
        if weighted:
            logw += math.log(1.0 + om)

        rows.append(
            {
                "m": m,
                "left": left,
                "right": right,
                "logm": math.log(m),
                "w": w,
                "logw": logw,
                "omega": om,
                "smooth_vec_m": smooth_vector(m),
                "smooth_vec_left": smooth_vector(left),
                "smooth_vec_right": smooth_vector(right),
                "cost_total": pair["cost_total"],
                "rho_m": rho_profile(m, basis),
                "rho_left": rho_profile(left, basis),
                "rho_right": rho_profile(right, basis),
            }
        )

    return rows


def cumulative_log_product(rows):
    out = []
    s = 0.0
    for r in rows:
        s += r["logw"]
        out.append((r["m"], r["logm"], s))
    return out


# --------------------------------------------------
# 5. Residuenanalyse
# --------------------------------------------------

def sign_change_rate(resid):
    signs = np.sign(resid)
    signs = signs[signs != 0]
    if len(signs) < 2:
        return 0.0
    return np.mean(signs[1:] != signs[:-1])


def lag1_corr(resid):
    if len(resid) < 3:
        return float("nan")
    a = resid[:-1]
    b = resid[1:]
    if np.std(a) == 0 or np.std(b) == 0:
        return 0.0
    return np.corrcoef(a, b)[0, 1]


def trig_scan(x, resid, omegas):
    x = np.asarray(x, dtype=float)
    resid = np.asarray(resid, dtype=float)
    if len(x) < 5:
        return None

    best = None
    ss_tot = np.sum((resid - resid.mean()) ** 2)

    for w in omegas:
        X = np.column_stack([np.cos(w * x), np.sin(w * x)])
        beta, *_ = np.linalg.lstsq(X, resid, rcond=None)
        fit = X @ beta
        ss_res = np.sum((resid - fit) ** 2)
        r2 = 0.0 if ss_tot == 0 else 1 - ss_res / ss_tot
        corr = 0.0
        if np.std(fit) > 0 and np.std(resid) > 0:
            corr = np.corrcoef(resid, fit)[0, 1]
        cand = {"omega": w, "r2": r2, "corr": corr}
        if best is None or cand["r2"] > best["r2"]:
            best = cand
    return best


# ============================================================
# Fits
# ============================================================

def linear_fit(cum_rows):
    x = np.array([r[1] for r in cum_rows], dtype=float)
    y = np.array([r[2] for r in cum_rows], dtype=float)
    X = np.column_stack([np.ones_like(x), x])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    fit = X @ beta
    resid = y - fit
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1.0 if ss_tot == 0 else 1 - np.sum((y - fit) ** 2) / ss_tot
    rmse = math.sqrt(np.mean((y - fit) ** 2))
    return {
        "x": x,
        "y": y,
        "fit": fit,
        "resid": resid,
        "beta": beta,
        "r2": r2,
        "rmse": rmse,
    }


def semismooth_fit(cum_rows):
    x = np.array([r[1] for r in cum_rows], dtype=float)
    y = np.array([r[2] for r in cum_rows], dtype=float)
    xp = x + 1.0
    X = np.column_stack([
        np.ones_like(x),
        x,
        np.log(xp),
        1.0 / xp,
        1.0 / (xp ** 3),
    ])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    fit = X @ beta
    resid = y - fit
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1.0 if ss_tot == 0 else 1 - np.sum((y - fit) ** 2) / ss_tot
    rmse = math.sqrt(np.mean((y - fit) ** 2))
    return {
        "x": x,
        "y": y,
        "fit": fit,
        "resid": resid,
        "beta": beta,
        "r2": r2,
        "rmse": rmse,
    }


def structure_record(n: int, basis=(2, 3, 5)) -> dict:
    return {
        "n": n,
        "smooth": smooth_part(n, basis),
        "kernel": rest_kernel(n, basis),
        "smooth_vec": smooth_vector(n),
        "C_fam": c_fam(n, basis),
        "epsilon": epsilon_B(n, basis),
        "omega_sig": omega_sig(n, basis),
    }


def records_up_to(limit: int, basis=(2, 3, 5), start: int = 2) -> list[dict]:
    return [structure_record(n, basis) for n in range(start, limit + 1)]


def top_by(records: list[dict], key: str, k: int = 10) -> list[dict]:
    return sorted(records, key=lambda r: r[key], reverse=True)[:k]


def print_records(records: list[dict]) -> None:
    print(
        f"{'n':>6} | {'smooth':>8} | {'kernel':>8} | {'v_(2,3,5)':<12} | "
        f"{'C_fam':>7} | {'epsilon':>7} | {'omega_sig':>9}"
    )
    print("-" * 76)
    for r in records:
        print(
            f"{r['n']:6d} | "
            f"{r['smooth']:8d} | "
            f"{r['kernel']:8d} | "
            f"{str(r['smooth_vec']):<12} | "
            f"{r['C_fam']:7.3f} | "
            f"{r['epsilon']:7.3f} | "
            f"{r['omega_sig']:9.3f}"
        )


# ============================================================
# Zusammenfassung eines Laufs
# ============================================================

def summarize_run(
    N: int,
    mode: str,
    weighted: bool,
    basis=(2, 3, 5),
    omega_min: float = 0.0,
    require_s000: bool = False,
    max_search: int = 500,
):
    rows = product_data_filtered(
        N=N,
        mode=mode,
        basis=basis,
        weighted=weighted,
        omega_min=omega_min,
        require_s000=require_s000,
        max_search=max_search,
    )

    if len(rows) < 3:
        return {
            "N": N,
            "mode": mode,
            "weighted": weighted,
            "omega_min": omega_min,
            "require_s000": require_s000,
            "count": len(rows),
            "m_min": None,
            "m_max": None,
            "mean_logw": None,
            "lin_r2": None,
            "semi_r2": None,
            "semi_rmse": None,
            "sign_change_rate": None,
            "lag1_corr": None,
            "best_omega": None,
            "best_r2": None,
            "best_corr": None,
        }

    cum = cumulative_log_product(rows)
    lin = linear_fit(cum)
    semi = semismooth_fit(cum)

    omegas = np.linspace(0.5, 30.0, 500)
    best = trig_scan(semi["x"], semi["resid"], omegas)

    return {
        "N": N,
        "mode": mode,
        "weighted": weighted,
        "omega_min": omega_min,
        "require_s000": require_s000,
        "count": len(rows),
        "m_min": rows[0]["m"],
        "m_max": rows[-1]["m"],
        "mean_logw": float(np.mean([r["logw"] for r in rows])),
        "lin_r2": lin["r2"],
        "semi_r2": semi["r2"],
        "semi_rmse": semi["rmse"],
        "sign_change_rate": sign_change_rate(semi["resid"]),
        "lag1_corr": lag1_corr(semi["resid"]),
        "best_omega": None if best is None else best["omega"],
        "best_r2": None if best is None else best["r2"],
        "best_corr": None if best is None else best["corr"],
    }


# ============================================================
# 4. Analyse einer Regime-Teilgruppe
# ============================================================

def analyze_rows(rows):
    if len(rows) < 6:
        return None

    rows = sorted(rows, key=lambda r: r["m"])
    cum = cumulative_log_product(rows)
    lin = linear_fit(cum)
    semi = semismooth_fit(cum)
    omegas = np.linspace(0.5, 30.0, 500)
    best = trig_scan(semi["x"], semi["resid"], omegas)

    return {
        "count": len(rows),
        "m_min": rows[0]["m"],
        "m_max": rows[-1]["m"],
        "mean_logw": float(np.mean([r["logw"] for r in rows])),
        "mean_omega": float(np.mean([r["omega"] for r in rows])),
        "mean_delta35_center": float(np.mean([delta_35(r["m"]) for r in rows])),
        "mean_v5_left": float(np.mean([v5(r["left"]) for r in rows])),
        "mean_v5_center": float(np.mean([v5(r["m"]) for r in rows])),
        "mean_v5_right": float(np.mean([v5(r["right"]) for r in rows])),
        "lin_r2": lin["r2"],
        "semi_r2": semi["r2"],
        "semi_rmse": semi["rmse"],
        "signchg": sign_change_rate(semi["resid"]),
        "lag1": lag1_corr(semi["resid"]),
        "omega_star": None if best is None else best["omega"],
        "r2_star": None if best is None else best["r2"],
        "corr_star": None if best is None else best["corr"],
    }


# ============================================================
# 6. Vergleichslauf alt vs. signaturgeometrisch
# ============================================================

def summarize_rows(rows):
    if len(rows) < 3:
        return None

    cum = cumulative_log_product(rows)
    lin = linear_fit(cum)
    semi = semismooth_fit(cum)
    omegas = np.linspace(0.5, 30.0, 500)
    best = trig_scan(semi["x"], semi["resid"], omegas)

    return {
        "count": len(rows),
        "m_min": rows[0]["m"],
        "m_max": rows[-1]["m"],
        "mean_logw": float(np.mean([r["logw"] for r in rows])),
        "mean_cost": float(np.mean([r.get("cost_total", 0.0) for r in rows])),
        "lin_r2": lin["r2"],
        "semi_r2": semi["r2"],
        "semi_rmse": semi["rmse"],
        "signchg": sign_change_rate(semi["resid"]),
        "lag1": lag1_corr(semi["resid"]),
        "omega_star": None if best is None else best["omega"],
        "r2_star": None if best is None else best["r2"],
        "corr_star": None if best is None else best["corr"],
    }


def print_compare(label, stats):
    if stats is None:
        print(f"{label}: zu wenige Daten")
        return

    omega_star = "nan" if stats["omega_star"] is None else f"{stats['omega_star']:.4f}"
    r2_star = "nan" if stats["r2_star"] is None else f"{stats['r2_star']:.6f}"

    print(
        f"{label:28s} | "
        f"count={stats['count']:4d} | "
        f"m=[{stats['m_min']},{stats['m_max']}] | "
        f"mean_logw={stats['mean_logw']:+.6f} | "
        f"mean_cost={stats['mean_cost']:.4f} | "
        f"lin_R2={stats['lin_r2']:.6f} | "
        f"semi_R2={stats['semi_r2']:.6f} | "
        f"RMSE={stats['semi_rmse']:.6f} | "
        f"signchg={stats['signchg']:.4f} | "
        f"lag1={stats['lag1']:.4f} | "
        f"ω*={omega_star} | "
        f"R2*={r2_star}"
    )


def print_regime_table(title, analyses):
    print("\n" + "=" * len(title))
    print(title)
    print("=" * len(title))
    header = (
        f"{'bucket':<28} | {'count':>5} | {'m_min':>5} | {'m_max':>5} | "
        f"{'mean_logw':>10} | {'mean_omega':>10} | {'mean_d35':>10} | "
        f"{'lin_R2':>8} | {'semi_R2':>8} | {'RMSE':>8} | "
        f"{'signchg':>8} | {'lag1':>8} | {'ω*':>8} | {'R2*':>8}"
    )
    print(header)
    print("-" * len(header))

    for bucket in sorted(analyses):
        stats = analyses[bucket]
        if stats is None:
            print(
                f"{bucket:<28} | {'<6':>5} | {'-':>5} | {'-':>5} | "
                f"{'-':>10} | {'-':>10} | {'-':>10} | "
                f"{'-':>8} | {'-':>8} | {'-':>8} | "
                f"{'-':>8} | {'-':>8} | {'-':>8} | {'-':>8}"
            )
            continue

        print(
            f"{bucket:<28} | "
            f"{stats['count']:>5d} | "
            f"{stats['m_min']:>5d} | "
            f"{stats['m_max']:>5d} | "
            f"{stats['mean_logw']:>10.6f} | "
            f"{stats['mean_omega']:>10.6f} | "
            f"{stats['mean_delta35_center']:>10.6f} | "
            f"{stats['lin_r2']:>8.6f} | "
            f"{stats['semi_r2']:>8.6f} | "
            f"{stats['semi_rmse']:>8.6f} | "
            f"{stats['signchg']:>8.4f} | "
            f"{stats['lag1']:>8.4f} | "
            f"{fmt(stats['omega_star'], 4):>8} | "
            f"{fmt(stats['r2_star'], 6):>8}"
        )


def print_group_table(title, groups, min_count=6):
    print("\n" + "=" * len(title))
    print(title)
    print("=" * len(title))

    header = (
        f"{'group':>20} | {'count':>5} | {'m_min':>5} | {'m_max':>5} | "
        f"{'mean_logw':>10} | {'v5L':>5} | {'v5C':>5} | {'v5R':>5} | "
        f"{'lin_R2':>8} | {'semi_R2':>8} | {'RMSE':>8} | "
        f"{'signchg':>8} | {'lag1':>8} | {'ω*':>8} | {'R2*':>8}"
    )
    print(header)
    print("-" * len(header))

    for key, rows in sorted(groups.items(), key=lambda kv: (len(kv[1]), kv[0]), reverse=True):
        stats = analyze_rows(rows)
        if stats is None or stats["count"] < min_count:
            continue

        print(
            f"{key:>20} | "
            f"{stats['count']:>5d} | "
            f"{stats['m_min']:>5d} | "
            f"{stats['m_max']:>5d} | "
            f"{stats['mean_logw']:>10.6f} | "
            f"{stats['mean_v5_left']:>5.2f} | "
            f"{stats['mean_v5_center']:>5.2f} | "
            f"{stats['mean_v5_right']:>5.2f} | "
            f"{stats['lin_r2']:>8.6f} | "
            f"{stats['semi_r2']:>8.6f} | "
            f"{stats['semi_rmse']:>8.6f} | "
            f"{stats['signchg']:>8.4f} | "
            f"{stats['lag1']:>8.4f} | "
            f"{stats['omega_star']:>8.4f} | "
            f"{stats['r2_star']:>8.6f}"
        )


def build_group_analyses(groups, min_count=6):
    analyses = {}
    for key, rows in groups.items():
        stats = analyze_rows(rows)
        if stats is None or stats["count"] < min_count:
            continue
        analyses[key] = stats
    return analyses


def print_mode_diff_table(title, analyses_bc, analyses_ae):
    print("\n" + "=" * len(title))
    print(title)
    print("=" * len(title))

    header = (
        f"{'group':>20} | {'n_BC':>5} | {'n_AE':>5} | "
        f"{'d_logw':>10} | {'d_linR2':>9} | {'d_semiR2':>10} | "
        f"{'d_RMSE':>9} | {'d_sign':>8} | {'d_lag1':>8} | {'d_ω*':>8}"
    )
    print(header)
    print("-" * len(header))

    common_keys = sorted(
        set(analyses_bc) & set(analyses_ae),
        key=lambda k: (analyses_bc[k]["count"] + analyses_ae[k]["count"], k),
        reverse=True,
    )

    for key in common_keys:
        bc = analyses_bc[key]
        ae = analyses_ae[key]
        d_omega = (
            None
            if bc["omega_star"] is None or ae["omega_star"] is None
            else bc["omega_star"] - ae["omega_star"]
        )
        print(
            f"{key:>20} | "
            f"{bc['count']:>5d} | "
            f"{ae['count']:>5d} | "
            f"{(bc['mean_logw'] - ae['mean_logw']):>10.6f} | "
            f"{(bc['lin_r2'] - ae['lin_r2']):>9.6f} | "
            f"{(bc['semi_r2'] - ae['semi_r2']):>10.6f} | "
            f"{(bc['semi_rmse'] - ae['semi_rmse']):>9.6f} | "
            f"{(bc['signchg'] - ae['signchg']):>8.4f} | "
            f"{(bc['lag1'] - ae['lag1']):>8.4f} | "
            f"{fmt(d_omega, 4):>8}"
        )


# ============================================================
# Tabellenausgabe
# ============================================================

def fmt(x, nd=6):
    if x is None:
        return "-"
    if isinstance(x, bool):
        return str(x)
    if isinstance(x, int):
        return str(x)
    return f"{x:.{nd}f}"


def print_summary_table(results, title):
    print("\n" + "=" * len(title))
    print(title)
    print("=" * len(title))
    header = (
        f"{'N':>6} | {'mode':>3} | {'wt':>2} | {'ω>':>5} | {'S000':>4} | "
        f"{'count':>5} | {'m_min':>5} | {'m_max':>5} | "
        f"{'mean_logw':>10} | {'lin_R2':>8} | {'semi_R2':>8} | "
        f"{'RMSE':>8} | {'signchg':>8} | {'lag1':>8} | "
        f"{'ω*':>8} | {'R2*':>8} | {'corr*':>8}"
    )
    print(header)
    print("-" * len(header))

    for r in results:
        print(
            f"{r['N']:>6} | "
            f"{r['mode']:>3} | "
            f"{('Y' if r['weighted'] else 'N'):>2} | "
            f"{fmt(r['omega_min'], 2):>5} | "
            f"{('Y' if r['require_s000'] else 'N'):>4} | "
            f"{fmt(r['count'], 0):>5} | "
            f"{fmt(r['m_min'], 0):>5} | "
            f"{fmt(r['m_max'], 0):>5} | "
            f"{fmt(r['mean_logw']):>10} | "
            f"{fmt(r['lin_r2']):>8} | "
            f"{fmt(r['semi_r2']):>8} | "
            f"{fmt(r['semi_rmse']):>8} | "
            f"{fmt(r['sign_change_rate']):>8} | "
            f"{fmt(r['lag1_corr']):>8} | "
            f"{fmt(r['best_omega']):>8} | "
            f"{fmt(r['best_r2']):>8} | "
            f"{fmt(r['best_corr']):>8}"
        )


# ============================================================
# Hauptlauf
# ============================================================

if __name__ == "__main__":
    basis = (2, 3, 5)
    N = 2000
    omega_min = 0.1
    require_s000 = True

    for weighted in [False, True]:
        mode_rows = {}
        full_analyses_by_mode = {}
        coarse_analyses_by_mode = {}

        for mode in ["BC", "AE"]:
            rows = product_data_regime(
                N=N,
                mode=mode,
                basis=basis,
                weighted=weighted,
                omega_min=omega_min,
                require_s000=require_s000,
                max_search=500,
            )
            mode_rows[mode] = rows

            print("\n" + "#" * 80)
            print(f"MODE={mode}, weighted={weighted}, N={N}, S=(0,0,0), omega>{omega_min}")
            print("#" * 80)

            total_stats = analyze_rows(rows)
            if total_stats is not None:
                print(
                    f"\nGesamt: count={total_stats['count']}, "
                    f"m=[{total_stats['m_min']},{total_stats['m_max']}], "
                    f"mean_logw={total_stats['mean_logw']:.6f}, "
                    f"lin_R2={total_stats['lin_r2']:.6f}, "
                    f"semi_R2={total_stats['semi_r2']:.6f}, "
                    f"RMSE={total_stats['semi_rmse']:.6f}, "
                    f"signchg={total_stats['signchg']:.4f}, "
                    f"lag1={total_stats['lag1']:.4f}, "
                    f"ω*={total_stats['omega_star']:.4f}, "
                    f"R2*={total_stats['r2_star']:.6f}"
                )
            else:
                print("\nGesamt: zu wenige Datenpunkte")
                continue

            full_groups = group_by_full_triplet_profile(rows)
            full_analyses = build_group_analyses(full_groups, min_count=6)
            full_analyses_by_mode[mode] = full_analyses
            print_group_table(
                f"Volle Tripelprofile | mode={mode}, weighted={weighted}",
                full_groups,
                min_count=6,
            )

            coarse_groups = group_by_flank_pattern(rows)
            coarse_analyses = build_group_analyses(coarse_groups, min_count=6)
            coarse_analyses_by_mode[mode] = coarse_analyses
            print_group_table(
                f"Grobe Flankenmuster | mode={mode}, weighted={weighted}",
                coarse_groups,
                min_count=6,
            )

        if "BC" in full_analyses_by_mode and "AE" in full_analyses_by_mode:
            print_mode_diff_table(
                f"Differenz BC - AE | volle Tripelprofile | weighted={weighted}",
                full_analyses_by_mode["BC"],
                full_analyses_by_mode["AE"],
            )

        if "BC" in coarse_analyses_by_mode and "AE" in coarse_analyses_by_mode:
            print_mode_diff_table(
                f"Differenz BC - AE | grobe Flankenmuster | weighted={weighted}",
                coarse_analyses_by_mode["BC"],
                coarse_analyses_by_mode["AE"],
            )