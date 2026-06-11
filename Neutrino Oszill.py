# Übergangsanalyse: Phase 3×3 (mod 36 // 12) und Restklassen 12×12 (mod 36)
import math
import random
from math import gcd
from pathlib import Path

import numpy as np
import pandas as pd
from collections import Counter

try:
    from scipy.stats import chisquare
except ImportError:
    chisquare = None

# erlaubte Klassen mod 36 (teilerfremd zu 6)
classes36 = [1, 5, 7, 11, 13, 17, 19, 23, 25, 29, 31, 35]
index36 = {r: i for i, r in enumerate(classes36)}

# Einheiten mod 210 (= 2·3·5·7), φ(210) = 48
mod = 210
classes210 = [r for r in range(mod) if gcd(r, mod) == 1]
index210 = {r: i for i, r in enumerate(classes210)}


def split_class(r: int) -> tuple[int, int]:
    r = int(r)
    return (r % 30, r % 7)


def class36(p: int) -> int:
    return int(p) % 36


def transition36(quadruplets: np.ndarray, max_gap: int) -> np.ndarray:
    M = np.zeros((12, 12), dtype=int)

    for i in range(len(quadruplets) - 1):
        p = int(quadruplets[i])
        q = int(quadruplets[i + 1])
        gap = q - p

        if gap <= max_gap:
            r1 = class36(p)
            r2 = class36(q)

            if r1 in index36 and r2 in index36:
                M[index36[r1], index36[r2]] += 1

    return M


def transition210(quadruplets: np.ndarray, max_gap: int) -> np.ndarray:
    """Übergänge zwischen Restklassen p mod 210 (nur Einheiten, index210)."""
    n = len(classes210)
    M = np.zeros((n, n), dtype=int)

    for i in range(len(quadruplets) - 1):
        p = int(quadruplets[i])
        q = int(quadruplets[i + 1])
        gap = q - p

        if gap <= max_gap:
            r1 = p % mod
            r2 = q % mod

            if r1 in index210 and r2 in index210:
                M[index210[r1], index210[r2]] += 1

    return M


def phase(p: int) -> int:
    """Phase des Vierling-Starts p: (p % 36) // 12 ∈ {0,1,2}."""
    return int(p) % 36 // 12


def transition_matrix(phases: list[int], centers: np.ndarray, max_gap: int) -> np.ndarray:
    M = np.zeros((3, 3), dtype=int)

    for i in range(len(centers) - 1):
        gap = int(centers[i + 1] - centers[i])
        if gap <= max_gap:
            a = phases[i]
            b = phases[i + 1]
            M[a, b] += 1

    return M


def v_p(n: int, p: int) -> int:
    """p-adische Bewertung v_p(n). Für n==0 wird 0 zurückgegeben."""
    if n == 0:
        return 0
    c = 0
    while n % p == 0:
        n //= p
        c += 1
    return c


def print_gap_valuation_summary(
    quadruplets: np.ndarray,
    centers: np.ndarray,
    *,
    small_thr: int = 100,
    large_thr: int = 1000,
) -> None:
    """Vergleicht kleine vs. große Abstände aufeinanderfolgender Vierling-Starts (p) bzgl. v_p(m), m=p+4."""
    q = np.asarray(quadruplets, dtype=np.int64)
    c = np.asarray(centers, dtype=np.int64)
    gaps = np.diff(q)

    data: list[dict] = []
    for i in range(len(q) - 1):
        m = int(c[i])
        gap = int(gaps[i])
        data.append({
            "gap": gap,
            "v3": v_p(m, 3),
            "v5": v_p(m, 5),
            "v7": v_p(m, 7),
        })

    small = [d for d in data if d["gap"] <= small_thr]
    large = [d for d in data if d["gap"] > large_thr]

    def avg(lst: list[dict], key: str) -> float:
        return float(np.mean([x[key] for x in lst])) if lst else 0.0

    def fraction_ge(lst: list[dict], key: str, k: int) -> float:
        if not lst:
            return 0.0
        return float(np.mean([x[key] >= k for x in lst]))

    print("\n=== Gaps vs. p-adische Bewertung der Zentren m = p+4 ===")
    print(
        f"(small: gap <= {small_thr}, large: gap > {large_thr}; "
        f"n_small={len(small)}, n_large={len(large)})"
    )
    print("\n--- Mittelwerte ---")
    for key in ("v3", "v5", "v7"):
        print(f"{key}: small={avg(small, key):.3f}, large={avg(large, key):.3f}")
    print("\n--- Anteil >= 2 ---")
    for key in ("v3", "v5", "v7"):
        print(f"{key}: small={fraction_ge(small, key, 2):.3f}, large={fraction_ge(large, key, 2):.3f}")


def gap_class(d: int, *, s_max: int = 100, m_max: int = 500) -> str:
    """S: d <= s_max, M: d <= m_max, sonst L."""
    if d <= s_max:
        return "S"
    if d <= m_max:
        return "M"
    return "L"


def shuffled_sequences(
    gaps: np.ndarray,
    *,
    trials: int = 1000,
    s_max: int = 100,
    m_max: int = 500,
    seed: int | None = 42,
) -> Counter[tuple[str, str, str]]:
    """
    Zählt S/M/L-Tripel auf zufällig permutierten Gap-Folgen.
    Über ``trials`` Läufe akkumuliert (wie im Original-Snippet).
    """
    gaps = np.asarray(gaps, dtype=np.int64).ravel()
    if len(gaps) < 3:
        return Counter()

    rng = random.Random(seed) if seed is not None else random.Random()
    counts: Counter[tuple[str, str, str]] = Counter()
    g_list = [int(x) for x in gaps]

    for _ in range(trials):
        g = g_list.copy()
        rng.shuffle(g)
        for i in range(len(g) - 2):
            seq = (
                gap_class(g[i], s_max=s_max, m_max=m_max),
                gap_class(g[i + 1], s_max=s_max, m_max=m_max),
                gap_class(g[i + 2], s_max=s_max, m_max=m_max),
            )
            counts[seq] += 1

    return counts


def print_top_gap_sequences(
    quadruplets: np.ndarray,
    *,
    top_n: int = 15,
    s_max: int = 100,
    m_max: int = 500,
    shuffle_trials: int = 1000,
    shuffle_seed: int | None = 42,
) -> None:
    """Zählt aufeinanderfolgende Tripel von Gap-Klassen (S/M/L) zwischen Vierling-Starts."""
    q = np.asarray(quadruplets, dtype=np.int64)
    gaps = np.diff(q)
    seq_counter: Counter[tuple[str, str, str]] = Counter()

    for i in range(len(gaps) - 2):
        g1 = gap_class(int(gaps[i]), s_max=s_max, m_max=m_max)
        g2 = gap_class(int(gaps[i + 1]), s_max=s_max, m_max=m_max)
        g3 = gap_class(int(gaps[i + 2]), s_max=s_max, m_max=m_max)
        seq_counter[(g1, g2, g3)] += 1

    print("\n=== Top Gap-Sequenzen (S/M/L) ===")
    print(f"(S: gap <= {s_max}, M: gap <= {m_max}, L: größer)")
    for seq, count in seq_counter.most_common(top_n):
        print(seq, count)

    if shuffle_trials > 0 and len(gaps) >= 3:
        baseline = shuffled_sequences(
            gaps,
            trials=shuffle_trials,
            s_max=s_max,
            m_max=m_max,
            seed=shuffle_seed,
        )
        print("\n--- vs. Permutations-Baseline (Zufallsreihenfolge der Gaps) ---")
        print(
            f"trials={shuffle_trials}, seed={shuffle_seed!r}; "
            "Spalten: beobachtet | Mittel/Trial | Summe baseline"
        )
        for seq, obs in seq_counter.most_common(top_n):
            bl_sum = baseline[seq]
            bl_mean = bl_sum / shuffle_trials
            print(f"  {seq}  {obs}  {bl_mean:.2f}  {bl_sum}")


def print_gap_class_pair_transitions(
    quadruplets: np.ndarray,
    *,
    s_max: int = 100,
    m_max: int = 500,
) -> None:
    """Zählt Übergänge (S/M/L) zwischen zwei aufeinanderfolgenden Gaps."""
    gaps = np.diff(np.asarray(quadruplets, dtype=np.int64))
    trans: Counter[tuple[str, str]] = Counter()
    for i in range(len(gaps) - 1):
        a = gap_class(int(gaps[i]), s_max=s_max, m_max=m_max)
        b = gap_class(int(gaps[i + 1]), s_max=s_max, m_max=m_max)
        trans[(a, b)] += 1

    print("\n=== Übergänge ===")
    print(f"(S: gap <= {s_max}, M: gap <= {m_max}, L: größer)")
    for k, v in sorted(trans.items(), key=lambda t: (-t[1], t[0])):
        print(k, v)


def _default_csv() -> Path:
    here = Path(__file__).resolve().parent
    for base in (here, Path.cwd()):
        cand = base / "primzahlvierlinge_1000.csv"
        if cand.is_file():
            return cand
    return here / "primzahlvierlinge_1000.csv"


def load_p_start_and_centers(csv_path: Path) -> tuple[np.ndarray, np.ndarray]:
    df = pd.read_csv(csv_path)
    keep = [c for c in df.columns if not str(c).startswith("Unnamed")]
    df = df[keep].copy()
    if "Mittel" not in df.columns:
        raise ValueError("CSV braucht Spalte 'Mittel' (Zentrum p+4).")

    df["Mittel"] = pd.to_numeric(df["Mittel"], errors="raise")
    df = df.sort_values("Mittel").reset_index(drop=True)

    p_start = (df["Mittel"] - 4).astype(np.int64).to_numpy()
    centers = df["Mittel"].to_numpy(dtype=np.int64)
    return p_start, centers


def print_rowwise_chisquare_uniform(M: np.ndarray, title: str) -> None:
    """H0: Bedingte Sprungverteilung in jeder Zeile ist gleichmäßig über k Spalten."""
    if chisquare is None:
        print(f"\n{title}: scipy nicht installiert, kein Chi².")
        return
    M = np.asarray(M, dtype=float)
    k = M.shape[1]
    print(f"\n{title}")
    print(f"(H0: pro Zeile gleichmäßig über {k} Spalten, Erwartung n/{k})")
    for i, row in enumerate(M):
        n = float(np.sum(row))
        if n <= 0:
            print(f"Row {i}: n=0 (übersprungen)")
            continue
        expected = np.full(k, n / k)
        chi, p = chisquare(row, expected)
        print(f"Row {i}: chi={chi:.2f}, p={p:.3f}, n={int(n)}")


def print_chisquare_m36_nonempty_targets(M36: np.ndarray) -> None:
    """H0: unter den tatsächlich getroffenen Zielklassen ist die Verteilung gleichmäßig."""
    if chisquare is None:
        print("\nM36 χ² (nur erlaubte Targets): scipy fehlt.")
        return
    M36 = np.asarray(M36, dtype=float)
    print("\nKORRIGIERTER χ²-TEST M36 (nur Spalten mit Treffer > 0)\n")
    for i, row in enumerate(M36):
        n = float(np.sum(row))
        if n == 0:
            continue
        nonzero_idx = np.where(row > 0)[0]
        observed = row[nonzero_idx]
        k = len(observed)
        if k < 2:
            print(f"Startklasse {classes36[i]}: nur ein Ziel (k=1), kein χ².\n")
            continue
        expected = np.full(k, n / k)
        chi, p = chisquare(observed, expected)
        tgs = [classes36[int(j)] for j in nonzero_idx]
        print(f"Startklasse {classes36[i]}:")
        print(f"  Targets: {tgs}")
        print(f"  counts : {observed.astype(int).tolist()}")
        print(f"  chi={chi:.2f}, p={p:.3f}, k={k}, n={int(n)}")
        print()


def print_chisquare_m36_robust(M36: np.ndarray) -> None:
    """Nur Ziele mit >=2 Treffern; Erwartung gleichmäßig über diese Ziele."""
    if chisquare is None:
        print("\nM36 χ² (robust): scipy fehlt.")
        return
    M36 = np.asarray(M36, dtype=float)
    print("\nROBUSTER χ²-TEST M36 (nur Targets mit count >= 2, Zeile n >= 10)\n")
    for i, row in enumerate(M36):
        n = float(np.sum(row))
        if n < 10:
            continue
        idx = np.where(row >= 2)[0]
        if len(idx) < 2:
            continue
        observed = row[idx]
        k = len(observed)
        n_sub = float(np.sum(observed))
        expected = np.full(k, n_sub / k)
        chi, p = chisquare(observed, expected)
        tgs = [classes36[int(j)] for j in idx]
        print(f"Startklasse {classes36[i]}:")
        print(f"  Targets>=2: {tgs}")
        print(f"  counts    : {observed.astype(int).tolist()}")
        print(f"  chi={chi:.2f}, p={p:.3f}")
        print()


def plot_m36_transition(M36: np.ndarray, save_path: Path | None = None) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nmatplotlib nicht installiert, Heatmap übersprungen.")
        return

    M36 = np.asarray(M36, dtype=float)
    if save_path is None:
        save_path = Path(__file__).resolve().parent / "Neutrino_M36_transition.png"

    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(M36, cmap="viridis")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.set_title("Transition Matrix mod 36")
    ax.set_xticks(range(12))
    ax.set_xticklabels([str(c) for c in classes36], rotation=45, ha="right")
    ax.set_yticks(range(12))
    ax.set_yticklabels([str(c) for c in classes36])
    ax.set_xlabel("to")
    ax.set_ylabel("from")
    fig.tight_layout()
    fig.savefig(save_path, dpi=150)
    print(f"\nHeatmap gespeichert: {save_path}")
    plt.close(fig)


def print_triplet_neutralitaet_mod210(p_start: np.ndarray, *, n_sim: int = 10000, seed: int = 12345) -> None:
    """
    Tripelneutralität der drei mod-210-Hauptfamilien {11, 101, 191}:
    gleitende Dreierfenster, Monte-Carlo, Permutation, Monochromie, Zyklus-Orientierung.
    Ankerklasse 5 wird vor der Familienliste entfernt (wie im Snippet).
    """
    p_start = np.asarray(p_start, dtype=np.int64)
    mask = (p_start % mod) != 5
    fam_residuals = (p_start[mask] % mod).astype(np.int64)
    allowed = {11, 101, 191}
    families = np.array([int(x) for x in fam_residuals if x in allowed], dtype=int)

    rng = np.random.default_rng(seed)
    fam_values = np.array([11, 101, 191])

    print("\n------------------------------------------------------------")
    print("8. Tripelneutralität der drei mod-210-Hauptfamilien")
    print("------------------------------------------------------------")

    print(f"Anzahl nach Entfernung der Ankerklasse: {len(families)}")

    fam_counts = Counter(families)
    print("\nFamilienhäufigkeiten:")
    for f in [11, 101, 191]:
        print(f"{f:3d}: {fam_counts.get(f, 0)}")

    triples: list[tuple[int, int, int]] = []
    neutral_flags: list[bool] = []

    for i in range(len(families) - 2):
        tri = (int(families[i]), int(families[i + 1]), int(families[i + 2]))
        triples.append(tri)
        neutral_flags.append(set(tri) == allowed)

    neutral_arr = np.array(neutral_flags, dtype=bool)
    n_triples = len(neutral_arr)
    n_neutral = int(neutral_arr.sum())
    neutral_rate = n_neutral / n_triples if n_triples else float("nan")

    print(f"\nGleitende Dreierfenster: {n_triples}")
    print(f"neutrale Tripel:         {n_neutral}")
    print(f"Neutralitätsrate:        {neutral_rate:.6f}")

    expected_rate_equal = 2 / 9
    print(f"Erwartung bei 1:1:1 unabhängig: {expected_rate_equal:.6f}")

    neutral_order_counts = Counter(
        tri for tri, flag in zip(triples, neutral_flags) if flag
    )

    print("\nNeutrale Ordnungen:")
    for tri, c in neutral_order_counts.most_common():
        print(f"{tri}: {c}")

    all_triple_counts = Counter(triples)

    print("\nHäufigste Tripel insgesamt:")
    for tri, c in all_triple_counts.most_common(12):
        print(f"{tri}: {c}")

    fam_probs = np.array(
        [
            float(fam_counts.get(11, 0)),
            float(fam_counts.get(101, 0)),
            float(fam_counts.get(191, 0)),
        ],
        dtype=float,
    )
    fam_probs = fam_probs / fam_probs.sum()

    sim_neutral_counts: list[int] = []

    for _ in range(n_sim):
        sim = rng.choice(fam_values, size=len(families), replace=True, p=fam_probs)
        sim_flags = [set(sim[i : i + 3]) == allowed for i in range(len(sim) - 2)]
        sim_arr = np.array(sim_flags, dtype=bool)
        sim_neutral_counts.append(int(sim_arr.sum()))

    sim_neutral_counts_arr = np.array(sim_neutral_counts, dtype=np.int64)

    p_greater = (np.sum(sim_neutral_counts_arr >= n_neutral) + 1) / (n_sim + 1)
    p_less = (np.sum(sim_neutral_counts_arr <= n_neutral) + 1) / (n_sim + 1)

    print("\nMonte-Carlo gegen unabhängige empirische Familienverteilung:")
    print(f"sim_mean_count = {sim_neutral_counts_arr.mean():.3f}")
    print(f"sim_sd_count   = {sim_neutral_counts_arr.std(ddof=1):.3f}")
    print(f"p_greater      = {p_greater:.5f}")
    print(f"p_less         = {p_less:.5f}")

    perm_neutral_counts: list[int] = []

    for _ in range(n_sim):
        perm = rng.permutation(families)
        perm_flags = [set(perm[i : i + 3]) == allowed for i in range(len(perm) - 2)]
        perm_neutral_counts.append(int(np.sum(perm_flags)))

    perm_neutral_counts_arr = np.array(perm_neutral_counts, dtype=np.int64)

    p_perm_greater = (np.sum(perm_neutral_counts_arr >= n_neutral) + 1) / (n_sim + 1)
    p_perm_less = (np.sum(perm_neutral_counts_arr <= n_neutral) + 1) / (n_sim + 1)

    print("\nPermutationstest, Familienzahlen exakt erhalten:")
    print(f"perm_mean_count = {perm_neutral_counts_arr.mean():.3f}")
    print(f"perm_sd_count   = {perm_neutral_counts_arr.std(ddof=1):.3f}")
    print(f"p_greater       = {p_perm_greater:.5f}")
    print(f"p_less          = {p_perm_less:.5f}")

    mono_flags = np.array([len(set(tri)) == 1 for tri in triples], dtype=bool)
    n_mono = int(mono_flags.sum())
    mono_rate = n_mono / n_triples if n_triples else float("nan")

    print("\nMonochrome Tripel:")
    print(f"monochrome Tripel: {n_mono}")
    print(f"Monochrom-Rate:    {mono_rate:.6f}")

    perm_mono_counts: list[int] = []

    for _ in range(n_sim):
        perm = rng.permutation(families)
        flags = [len(set(perm[i : i + 3])) == 1 for i in range(len(perm) - 2)]
        perm_mono_counts.append(int(np.sum(flags)))

    perm_mono_counts_arr = np.array(perm_mono_counts, dtype=np.int64)

    p_mono_greater = (np.sum(perm_mono_counts_arr >= n_mono) + 1) / (n_sim + 1)
    p_mono_less = (np.sum(perm_mono_counts_arr <= n_mono) + 1) / (n_sim + 1)

    print("\nMonochromie gegen Permutationstest:")
    print(f"perm_mean_mono = {perm_mono_counts_arr.mean():.3f}")
    print(f"perm_sd_mono   = {perm_mono_counts_arr.std(ddof=1):.3f}")
    print(f"p_greater      = {p_mono_greater:.5f}")
    print(f"p_less         = {p_mono_less:.5f}")

    cycle_plus = {
        (11, 101, 191),
        (101, 191, 11),
        (191, 11, 101),
    }
    cycle_minus = {
        (11, 191, 101),
        (191, 101, 11),
        (101, 11, 191),
    }

    n_plus = sum(1 for tri in triples if tri in cycle_plus)
    n_minus = sum(1 for tri in triples if tri in cycle_minus)

    print("\nZyklische Orientierung neutraler Tripel:")
    print(f"plus-Zyklen:  {n_plus}")
    print(f"minus-Zyklen: {n_minus}")
    print(f"Differenz:    {n_plus - n_minus}")

    if n_plus + n_minus > 0:
        z_orient = (n_plus - 0.5 * (n_plus + n_minus)) / math.sqrt(
            0.25 * (n_plus + n_minus)
        )
        print(f"z Orientierung gegen 50:50 ≈ {z_orient:.3f}")


if __name__ == "__main__":
    path = _default_csv()
    if not path.is_file():
        raise SystemExit(
            f"Keine primzahlvierlinge_1000.csv unter {path.parent} oder CWD."
        )

    quadruplets, centers = load_p_start_and_centers(path)
    phases = [phase(int(p)) for p in quadruplets]
    residues = [int(p) % mod for p in quadruplets]
    counts = Counter(residues)
    print("Unique residues:", sorted(set(int(p) % mod for p in quadruplets)))

    print("\n=== mod 210 Häufigkeiten ===")
    for r in classes210:
        print(r, counts[r])

    if chisquare is not None:
        obs = np.array([counts[r] for r in classes210], dtype=float)
        expected = np.full(len(classes210), np.sum(obs) / len(classes210))
        chi, p = chisquare(obs, expected)
        print(f"\nChi² mod210: {chi:.2f}, p={p:.3f}")
    else:
        print("\nChi² mod210: scipy nicht verfügbar.")

    valid = [5, 11, 101, 191]
    if chisquare is not None:
        obs = np.array([counts[r] for r in valid], dtype=float)
        expected = np.full(len(valid), np.sum(obs) / len(valid))
        chi, p = chisquare(obs, expected)
        print(chi, p)

    pairs = [split_class(int(p) % mod) for p in quadruplets]
    counts_pairs = Counter(pairs)
    print("\n=== Paare (p mod 210 → mod 30, mod 7) ===")
    for pair, c in counts_pairs.most_common():
        print(pair, c)

    gaps = np.diff(np.asarray(quadruplets, dtype=np.int64))
    small = [int(p) % mod for i, p in enumerate(quadruplets[:-1]) if int(gaps[i]) <= 100]
    large = [int(p) % mod for i, p in enumerate(quadruplets[:-1]) if int(gaps[i]) > 1000]
    c_small = Counter(small)
    c_large = Counter(large)
    print("\n=== small vs large ===")
    for r in classes210[:20]:
        print(r, c_small[r], c_large[r])

    max_gap = 1000

    M210 = transition210(quadruplets, max_gap=max_gap)
    print("\nÜbergangsmatrix 48×48 (mod 210, Zeilen/Spalten: classes210), gap <=", max_gap)
    print(M210)

    M = transition_matrix(phases, centers, max_gap=max_gap)
    print("Übergangsmatrix 3×3 (Phase i → j), gap <=", max_gap)
    print(M)
    print("Zeilensummen:", M.sum(axis=1))
    print_rowwise_chisquare_uniform(M, "χ² je Zeile (3×3-M, Zielphasen gleichverteilt)")

    print_gap_valuation_summary(quadruplets, centers)
    print_top_gap_sequences(quadruplets)
    print_gap_class_pair_transitions(quadruplets)

    M36 = transition36(quadruplets, max_gap=max_gap)
    print("\nÜbergangsmatrix 12×12 (Rest mod 36), Zeile/Spalte: classes36")
    print("classes36 =", classes36)
    print(M36)
    print_chisquare_m36_nonempty_targets(M36)
    print_chisquare_m36_robust(M36)
    plot_m36_transition(M36)

    print_triplet_neutralitaet_mod210(quadruplets)
