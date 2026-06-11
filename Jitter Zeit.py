# ============================================================
# Bamberger Quantum-Dice-Modul / Zeitwürfel
# Entanglement-Tempered Timestamp Spectroscopy
# ============================================================

import math
import cmath
import random
import numpy as np
import pandas as pd


# ------------------------------------------------------------
# 1. Kleine Quaternionen-Klasse für EABC-Zustände
# ------------------------------------------------------------

class Quaternion:
    """
    Einfache Quaternion q = e + a i + b j + c k.
    Für das Bamberger EABC-Modul:
        E -> 1
        A -> i
        B -> j
        C -> k
    """

    def __init__(self, e=0.0, a=0.0, b=0.0, c=0.0):
        self.e = float(e)
        self.a = float(a)
        self.b = float(b)
        self.c = float(c)

    def __add__(self, other):
        return Quaternion(
            self.e + other.e,
            self.a + other.a,
            self.b + other.b,
            self.c + other.c
        )

    def __sub__(self, other):
        return Quaternion(
            self.e - other.e,
            self.a - other.a,
            self.b - other.b,
            self.c - other.c
        )

    def __mul__(self, other):
        """
        Hamilton-Produkt:
        i^2 = j^2 = k^2 = ijk = -1
        ij = k, jk = i, ki = j
        ji = -k, kj = -i, ik = -j
        """
        e1, a1, b1, c1 = self.e, self.a, self.b, self.c
        e2, a2, b2, c2 = other.e, other.a, other.b, other.c

        return Quaternion(
            e1*e2 - a1*a2 - b1*b2 - c1*c2,
            e1*a2 + a1*e2 + b1*c2 - c1*b2,
            e1*b2 - a1*c2 + b1*e2 + c1*a2,
            e1*c2 + a1*b2 - b1*a2 + c1*e2
        )

    def conj(self):
        return Quaternion(self.e, -self.a, -self.b, -self.c)

    def norm2(self):
        return self.e**2 + self.a**2 + self.b**2 + self.c**2

    def norm(self):
        return math.sqrt(self.norm2())

    def as_tuple(self):
        return (self.e, self.a, self.b, self.c)

    def __repr__(self):
        return f"({self.e:+.3f} {self.a:+.3f}i {self.b:+.3f}j {self.c:+.3f}k)"


# EABC-Basis
Q_E = Quaternion(1, 0, 0, 0)
Q_A = Quaternion(0, 1, 0, 0)
Q_B = Quaternion(0, 0, 1, 0)
Q_C = Quaternion(0, 0, 0, 1)

EABC_TO_Q = {
    "E": Q_E,
    "A": Q_A,
    "B": Q_B,
    "C": Q_C
}


# ------------------------------------------------------------
# 2. EABC-Klassifikation modulo 12
# ------------------------------------------------------------

def eabc_label_mod12(n):
    """
    Klassifiziert eine ganze Zahl n nach dem EABC-Schema modulo 12.

    E ≡ 1 mod 12
    A ≡ 5 mod 12
    B ≡ 7 mod 12
    C ≡ 11 mod 12

    Andere Reste gelten als glatte oder ausgeschlossene Zwischenzustände.
    """
    r = n % 12

    if r == 1:
        return "E"
    elif r == 5:
        return "A"
    elif r == 7:
        return "B"
    elif r == 11:
        return "C"
    else:
        return "X"


def eabc_quaternion(label):
    """
    Gibt den zugehörigen Quaternionenbasiszustand zurück.
    Für X wird der Nullzustand verwendet.
    """
    return EABC_TO_Q.get(label, Quaternion(0, 0, 0, 0))


# ------------------------------------------------------------
# 3. Entanglement-Jitter und effektive Breite
# ------------------------------------------------------------

def sigma_t_entanglement(S, sigma_min=1e-3, sigma_0=1.0, gamma=1.0):
    """
    Entanglement-abhängiger Zeitjitter:

        sigma_t(S) = sigma_min + sigma_0 * exp(-gamma S)

    S ist die zugängliche Verschränkungsentropie.
    """
    return sigma_min + sigma_0 * math.exp(-gamma * S)


def delta_eff(T, S, lam=1.0):
    """
    Temperatur- und verschränkungsabhängige effektive Vierlingsbreite:

        Delta_eff(T,S)
        = 8 + 2(1 - exp(-T/(2pi)))(1 - exp(-lambda S))

    Grenzfälle:
        T -> 0      : Delta_eff -> 8
        T ~ 2pi, S groß : Delta_eff nähert sich 8 + 2(1-e^-1) ≈ 9.264
        T groß, S groß  : Delta_eff -> 10

    Hinweis:
    Wenn man bei T = 2pi bereits näher an 10 liegen möchte,
    kann man den Temperaturfaktor schärfer wählen.
    """
    temp_factor = 1.0 - math.exp(-T / (2.0 * math.pi))
    ent_factor = 1.0 - math.exp(-lam * S)
    return 8.0 + 2.0 * temp_factor * ent_factor


def delta_eff_sharp(T, S, lam=1.0, beta=3.0):
    """
    Alternative schärfere Temperierung:

        Delta_eff(T,S)
        = 8 + 2(1 - exp(-beta T/(2pi)))(1 - exp(-lambda S))

    Für beta > 1 öffnet sich die Hülle bereits bei T ≈ 2pi stärker.
    """
    temp_factor = 1.0 - math.exp(-beta * T / (2.0 * math.pi))
    ent_factor = 1.0 - math.exp(-lam * S)
    return 8.0 + 2.0 * temp_factor * ent_factor


# ------------------------------------------------------------
# 4. Zeitwürfel-Simulation
# ------------------------------------------------------------

def simulate_zeitwuerfel(
    N=1000,
    T=2.0 * math.pi,
    S=2.0,
    omega_0=1.0,
    sigma_min=1e-3,
    sigma_0=0.05,
    gamma=1.0,
    lam=1.0,
    alpha_scale=1.0,
    seed=42,
    use_sharp_delta=True,
    beta=3.0
):
    """
    Simuliert den Bamberger Zeitwürfel.

    Parameter:
        N             Anzahl der Zeitintervalle
        T             Temperaturparameter
        S             Verschränkungsparameter
        omega_0       Grundfrequenz
        sigma_min     minimaler Jitterboden
        sigma_0       anfänglicher Jitter
        gamma         Entanglement-Dämpfung des Jitters
        lam           Entanglement-Faktor der effektiven Breite
        alpha_scale   Skalenfaktor zur EABC-Klassifikation
        seed          Zufallsseed
        use_sharp_delta  verwendet schärfere Delta_eff-Form
        beta          Temperatur-Schärfeparameter

    Rückgabe:
        DataFrame mit t_k, Delta t_k, EABC-Label, Quaternion, Schwung.
    """

    random.seed(seed)
    np.random.seed(seed)

    sig = sigma_t_entanglement(
        S,
        sigma_min=sigma_min,
        sigma_0=sigma_0,
        gamma=gamma
    )

    if use_sharp_delta:
        Deff = delta_eff_sharp(T, S, lam=lam, beta=beta)
    else:
        Deff = delta_eff(T, S, lam=lam)

    base_dt = Deff / omega_0

    intervals = []
    timestamps = []

    t = 0.0

    for k in range(N):
        eps = np.random.normal(loc=0.0, scale=sig)
        dt = base_dt + eps

        # Sicherheit: keine negativen Intervalle
        if dt <= 0:
            dt = base_dt

        t += dt

        intervals.append(dt)
        timestamps.append(t)

    rows = []

    q_list = []

    for k in range(N):
        n_k = int(round(alpha_scale * timestamps[k]))
        label = eabc_label_mod12(n_k)
        q = eabc_quaternion(label)
        q_list.append(q)

        rows.append({
            "k": k,
            "t_k": timestamps[k],
            "Delta_t": intervals[k],
            "n_k": n_k,
            "mod12": n_k % 12,
            "label": label,
            "q_e": q.e,
            "q_a": q.a,
            "q_b": q.b,
            "q_c": q.c,
            "sigma_t": sig,
            "Delta_eff": Deff,
            "T": T,
            "S_ent": S
        })

    # Schwung u_k = conjugate(q_k) * q_{k+1}
    for k in range(N - 1):
        u = q_list[k].conj() * q_list[k + 1]
        rows[k]["u_e"] = u.e
        rows[k]["u_a"] = u.a
        rows[k]["u_b"] = u.b
        rows[k]["u_c"] = u.c
        rows[k]["u_norm"] = u.norm()

    # letzter Eintrag hat keinen Nachfolger
    rows[-1]["u_e"] = np.nan
    rows[-1]["u_a"] = np.nan
    rows[-1]["u_b"] = np.nan
    rows[-1]["u_c"] = np.nan
    rows[-1]["u_norm"] = np.nan

    return pd.DataFrame(rows)


# ------------------------------------------------------------
# 5. Quaternionisches Spektrum
# ------------------------------------------------------------

def quaternionic_spectrum(df, Omega_values):
    """
    Berechnet das quaternionische Spektrum:

        P_H(Omega) =
        | sum_k q_k exp(-i Omega Delta_t_k) |^2

    Praktisch wird q_k komponentenweise komplex gewichtet:
        Q_E(Omega), Q_A(Omega), Q_B(Omega), Q_C(Omega)

    und dann
        P = |Q_E|^2 + |Q_A|^2 + |Q_B|^2 + |Q_C|^2

    berechnet.
    """

    Delta_t = df["Delta_t"].values
    q_e = df["q_e"].values
    q_a = df["q_a"].values
    q_b = df["q_b"].values
    q_c = df["q_c"].values

    spec_rows = []

    for Om in Omega_values:
        phase = np.exp(-1j * Om * Delta_t)

        QE = np.sum(q_e * phase)
        QA = np.sum(q_a * phase)
        QB = np.sum(q_b * phase)
        QC = np.sum(q_c * phase)

        P_E = abs(QE)**2
        P_A = abs(QA)**2
        P_B = abs(QB)**2
        P_C = abs(QC)**2
        P_total = P_E + P_A + P_B + P_C

        spec_rows.append({
            "Omega": Om,
            "P_E": P_E,
            "P_A": P_A,
            "P_B": P_B,
            "P_C": P_C,
            "P_total": P_total
        })

    return pd.DataFrame(spec_rows)


# ------------------------------------------------------------
# 5b. Aktive EABC-Folge, Schwung nur zwischen EABC, zentriertes Spektrum
# ------------------------------------------------------------

def active_eabc_sequence(df):
    """
    Filtert die aktive EABC-Folge.
    Entfernt alle X-Zustände und berechnet aktive Intervalle
    zwischen aufeinanderfolgenden aktiven Ereignissen (Zeitabstände in t_k).
    """
    active = df[df["label"].isin(["E", "A", "B", "C"])].copy()
    active = active.reset_index(drop=True)

    active["Delta_t_active"] = active["t_k"].diff()
    active.loc[0, "Delta_t_active"] = np.nan

    return active


def schwung_label(row):
    """
    Liest u_e, u_a, u_b, u_c aus einer Zeile (Series oder dict-kompatibel).
    Dominante reelle Komponente von u (nach Betrag); Vorzeichen als +/- vor dem Kanal.
    """
    comps = {
        "E": row["u_e"],
        "A": row["u_a"],
        "B": row["u_b"],
        "C": row["u_c"],
    }

    key = max(comps, key=lambda k: abs(float(comps[k])))
    val = float(comps[key])

    if abs(val) < 1e-12:
        return "0"

    sign = "+" if val > 0 else "-"
    return sign + key


def active_transition_schwung(active_df):
    """
    Berechnet den relativen Schwung u = conj(q_k) * q_{k+1} nur zwischen
    aufeinanderfolgenden Einträgen der bereits gefilterten aktiven Folge.
    """
    tuples = active_df.itertuples(index=False)
    q_list = [
        Quaternion(row.q_e, row.q_a, row.q_b, row.q_c)
        for row in tuples
    ]

    rows = []
    n = len(q_list)
    for k in range(n - 1):
        q0 = q_list[k]
        q1 = q_list[k + 1]
        u = q0.conj() * q1

        row_dict = {
            "k_active": k,
            "label_from": active_df.loc[k, "label"],
            "label_to": active_df.loc[k + 1, "label"],
            "t_from": active_df.loc[k, "t_k"],
            "t_to": active_df.loc[k + 1, "t_k"],
            "Delta_t_active": active_df.loc[k + 1, "t_k"] - active_df.loc[k, "t_k"],
            "u_e": u.e,
            "u_a": u.a,
            "u_b": u.b,
            "u_c": u.c,
            "u_norm": u.norm(),
        }
        row_dict["schwung_label"] = schwung_label(pd.Series(row_dict))
        rows.append(row_dict)

    return pd.DataFrame(rows)


def transition_matrix(active_transitions):
    """
    Erzeugt eine Übergangsmatrix label_from -> label_to
    für die aktive EABC-Folge.
    """
    labels = ["E", "A", "B", "C"]

    mat = pd.crosstab(
        active_transitions["label_from"],
        active_transitions["label_to"],
    )

    mat = mat.reindex(index=labels, columns=labels, fill_value=0)

    return mat


def schwung_distribution(active_transitions):
    """
    Zählt die Schwunglabels, z.B. -A, +C, +A, +B.
    """
    return active_transitions["schwung_label"].value_counts().sort_index()


def print_transition_report(active_transitions):
    """
    Kompakter Bericht zur aktiven EABC-Dynamik.
    """
    mat = transition_matrix(active_transitions)
    schwung = schwung_distribution(active_transitions)

    print("\n=== Aktive EABC-Übergangsmatrix ===")
    print(mat)

    print("\n=== Schwung-Verteilung ===")
    print(schwung)

    total = len(active_transitions)

    ideal_patterns = {
        ("A", "E"),
        ("E", "C"),
        ("C", "B"),
        ("B", "A"),
    }

    ideal_count = 0

    for _, row in active_transitions.iterrows():
        pair = (row["label_from"], row["label_to"])
        if pair in ideal_patterns:
            ideal_count += 1

    print("\n=== Idealrotations-Anteil ===")
    print(f"Idealübergänge: {ideal_count} / {total}")
    if total == 0:
        print("Anteil: — (keine Übergänge)")
    else:
        print(f"Anteil: {ideal_count / total:.6f}")

    return mat, schwung


# Funktionen
def rotation_coherence_8(active_transitions):
    """
    Kohärenz der nackten 8-Breitenstruktur.

    Delta = 8 entspricht 8 ≡ -4 mod 12.

    Aktiv entstehen zwei Zweierzyklen:
        A <-> E
        C <-> B
    """
    ideal = {
        ("A", "E"),
        ("E", "A"),
        ("C", "B"),
        ("B", "C"),
    }

    total = len(active_transitions)
    if total == 0:
        return float("nan")

    hits = 0
    for _, row in active_transitions.iterrows():
        pair = (row["label_from"], row["label_to"])
        if pair in ideal:
            hits += 1

    return hits / total


def rotation_coherence_10(active_transitions):
    """
    Kohärenz der solenoidischen 10-Hüllrotation.

    Delta = 10 entspricht 10 ≡ -2 mod 12.

    Aktiv entsteht:
        A -> E -> C -> B -> A
    """
    ideal = {
        ("A", "E"),
        ("E", "C"),
        ("C", "B"),
        ("B", "A"),
    }

    total = len(active_transitions)
    if total == 0:
        return float("nan")

    hits = 0
    for _, row in active_transitions.iterrows():
        pair = (row["label_from"], row["label_to"])
        if pair in ideal:
            hits += 1

    return hits / total


def run_simple_critical_scan(
    T_over_min=0.500,
    T_over_max=0.580,
    steps=161,
    S=3.0,
    N=3000,
    seed=123,
    beta=3.0,
    csv_path="bamberger_zeitwuerfel_simple_critical_scan.csv",
):
    rows = []

    T_over_values = np.linspace(T_over_min, T_over_max, steps)

    for T_over in T_over_values:
        T = T_over * 2.0 * math.pi

        df_tmp = simulate_zeitwuerfel(
            N=N,
            T=T,
            S=S,
            seed=seed,
            use_sharp_delta=True,
            beta=beta,
        )

        active_tmp = active_eabc_sequence(df_tmp)
        trans_tmp = active_transition_schwung(active_tmp)

        R8 = rotation_coherence_8(trans_tmp)
        R10 = rotation_coherence_10(trans_tmp)

        D = R10 - R8
        O = (R10 - R8) / (R10 + R8) if (R10 + R8) != 0 else float("nan")

        rows.append({
            "T_over_2pi": T_over,
            "T": T,
            "Delta_eff": df_tmp["Delta_eff"].iloc[0],
            "sigma_t": df_tmp["sigma_t"].iloc[0],
            "N_active": len(active_tmp),
            "N_trans": len(trans_tmp),
            "R8": R8,
            "R10": R10,
            "D_order": D,
            "O_order_norm": O,
        })

    scan = pd.DataFrame(rows)
    scan.to_csv(csv_path, index=False)

    idx_crit = scan["D_order"].abs().idxmin()
    crit = scan.loc[idx_crit]

    print("\n=== Einfacher kritischer Scan ===")
    print("Datei:", csv_path)

    print("\n=== Kandidat T_c ===")
    print("T/(2pi)  :", crit["T_over_2pi"])
    print("T        :", crit["T"])
    print("Delta_eff:", crit["Delta_eff"])
    print("sigma_t  :", crit["sigma_t"])
    print("N_active :", crit["N_active"])
    print("N_trans  :", crit["N_trans"])
    print("R8       :", crit["R8"])
    print("R10      :", crit["R10"])
    print("D_order  :", crit["D_order"])
    print("O_order  :", crit["O_order_norm"])

    return scan


def run_temperature_contrast():
    """
    Vergleicht den Low-T-8-Modus mit dem High-T-10-Modus (gleicher Seed S etc.).
    Schreibt CSV-Dateien für Übergangsmatrizen, Schwungverteilungen und Kontrasttabelle.
    """
    configs = [
        ("Low-T", 0.01),
        ("High-T", 2.0 * math.pi),
    ]

    rows = []

    for name, T in configs:
        df_tmp = simulate_zeitwuerfel(
            N=2000,
            T=T,
            S=3.0,
            seed=123,
            use_sharp_delta=True,
            beta=3.0,
        )

        active_tmp = active_eabc_sequence(df_tmp)
        trans_tmp = active_transition_schwung(active_tmp)

        R8 = rotation_coherence_8(trans_tmp)
        R10 = rotation_coherence_10(trans_tmp)

        mat_tmp = transition_matrix(trans_tmp)
        schwung_tmp = schwung_distribution(trans_tmp)

        rows.append({
            "mode": name,
            "T": T,
            "T_over_2pi": T / (2.0 * math.pi),
            "Delta_eff": df_tmp["Delta_eff"].iloc[0],
            "sigma_t": df_tmp["sigma_t"].iloc[0],
            "N_active": len(active_tmp),
            "N_trans": len(trans_tmp),
            "R8": R8,
            "R10": R10,
            "D_order": R10 - R8,
        })

        print(f"\n=== {name} ===")
        print("Delta_eff:", df_tmp["Delta_eff"].iloc[0])
        print("sigma_t  :", df_tmp["sigma_t"].iloc[0])
        print("Aktive Ereignisse:", len(active_tmp))
        print("Aktive Übergänge :", len(trans_tmp))
        print("R8 :", R8)
        print("R10:", R10)
        print("D_order (R10 - R8):", R10 - R8)

        print("\nÜbergangsmatrix:")
        print(mat_tmp)

        print("\nSchwung-Verteilung:")
        print(schwung_tmp)

        mat_tmp.to_csv(f"bamberger_zeitwuerfel_transition_matrix_{name}.csv")
        schwung_tmp.to_csv(f"bamberger_zeitwuerfel_schwung_distribution_{name}.csv")

    contrast = pd.DataFrame(rows)
    contrast.to_csv("bamberger_zeitwuerfel_temperature_contrast.csv", index=False)

    print("\n=== Temperatur-Kontrast ===")
    print(contrast)

    return contrast


def run_temperature_order_scan(
    T_min=0.0,
    T_max=2.0 * math.pi,
    steps=41,
    S=3.0,
    N=2000,
    seed=123,
    beta=3.0,
):
    """
    Scan über T: Kohärenzen R8, R10, Differenz D=R10-R8 und
    normierte Größe O=(R10-R8)/(R10+R8), falls R10+R8 ungleich 0.

    Schreibt die Datei bamberger_zeitwuerfel_temperature_order_scan.csv.
    """
    T_values = np.linspace(T_min, T_max, steps)

    rows = []

    for T in T_values:
        df_tmp = simulate_zeitwuerfel(
            N=N,
            T=T,
            S=S,
            seed=seed,
            use_sharp_delta=True,
            beta=beta,
        )

        active_tmp = active_eabc_sequence(df_tmp)
        trans_tmp = active_transition_schwung(active_tmp)

        R8 = rotation_coherence_8(trans_tmp)
        R10 = rotation_coherence_10(trans_tmp)

        D = R10 - R8
        denom = R10 + R8
        O = (R10 - R8) / denom if denom != 0 else np.nan

        rows.append({
            "T": T,
            "T_over_2pi": T / (2.0 * math.pi),
            "Delta_eff": df_tmp["Delta_eff"].iloc[0],
            "sigma_t": df_tmp["sigma_t"].iloc[0],
            "N_active": len(active_tmp),
            "N_trans": len(trans_tmp),
            "R8": R8,
            "R10": R10,
            "D_order": D,
            "O": O,
        })

    scan = pd.DataFrame(rows)
    scan.to_csv("bamberger_zeitwuerfel_temperature_order_scan.csv", index=False)

    print("\n=== Temperatur-Order-Scan ===")
    print(scan.head(8))
    print("...")
    print(scan.tail(5))

    return scan


def shannon_entropy_from_counts(counts):
    """
    Shannon-Entropie einer Zählverteilung (natürlicher Logarithmus).
    counts kann eine Liste, Series oder ndarray sein.
    """
    arr = np.asarray(counts, dtype=float)
    total = arr.sum()

    if total <= 0:
        return 0.0

    p = arr[arr > 0] / total
    return float(-np.sum(p * np.log(p)))


def transition_entropy(active_transitions):
    """
    Entropie der aktiven Übergangsmatrix.
    Maximal wären 16 Übergänge E/A/B/C -> E/A/B/C.
    """
    mat = transition_matrix(active_transitions)
    return shannon_entropy_from_counts(mat.values.flatten())


def schwung_entropy(active_transitions):
    """
    Entropie der Schwung-Verteilung.
    """
    dist = schwung_distribution(active_transitions)
    return shannon_entropy_from_counts(dist.values)


def active_interval_stats(active_df):
    """
    Statistiken der aktiven Zeitabstände.
    Nutzt Delta_t_active, also Abstände zwischen aktiven Ereignissen.
    """
    if "Delta_t_active" not in active_df.columns:
        active_df = active_df.copy()
        active_df["Delta_t_active"] = active_df["t_k"].diff()

    x = active_df["Delta_t_active"].dropna().values

    if len(x) == 0:
        return {
            "active_dt_mean": np.nan,
            "active_dt_std": np.nan,
            "active_dt_cv": np.nan,
            "active_dt_min": np.nan,
            "active_dt_max": np.nan,
        }

    mean = float(np.mean(x))
    std = float(np.std(x))

    return {
        "active_dt_mean": mean,
        "active_dt_std": std,
        "active_dt_cv": std / mean if mean != 0 else np.nan,
        "active_dt_min": float(np.min(x)),
        "active_dt_max": float(np.max(x)),
    }


def estimate_hurst_rs(x, min_chunk=16, max_chunk=None):
    """
    Einfache R/S-Hurst-Schätzung.
    Kein Beweis für Fraktalität, aber ein nützlicher Indikator.

    H ~ 0.5: weißes Rauschen
    H > 0.5: persistente Korrelation
    H < 0.5: antipersistente Struktur
    """
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]

    N = len(x)

    if N < 4 * min_chunk:
        return np.nan

    if max_chunk is None:
        max_chunk = N // 4

    sizes = []
    rs_values = []

    chunk = min_chunk
    while chunk <= max_chunk:
        n_chunks = N // chunk
        if n_chunks < 2:
            break

        rs_local = []

        for i in range(n_chunks):
            seg = x[i * chunk : (i + 1) * chunk]
            seg = seg - np.mean(seg)

            y = np.cumsum(seg)
            R = np.max(y) - np.min(y)
            S = np.std(seg)

            if S > 0:
                rs_local.append(R / S)

        if len(rs_local) > 0:
            sizes.append(chunk)
            rs_values.append(np.mean(rs_local))

        chunk *= 2

    if len(sizes) < 2:
        return np.nan

    log_sizes = np.log(np.asarray(sizes))
    log_rs = np.log(np.asarray(rs_values))

    slope, _intercept = np.polyfit(log_sizes, log_rs, 1)

    return float(slope)


def run_critical_temperature_scan(
    T_over_2pi_min=0.500,
    T_over_2pi_max=0.580,
    steps=401,
    S=3.0,
    N=5000,
    seed=123,
    beta=3.0,
    save_prefix="bamberger_zeitwuerfel_critical_scan",
):
    """
    Hochauflösender Scan um den vermuteten kritischen Bereich.

    Misst:
        Delta_eff
        R8
        R10
        D_order = R10 - R8
        O_order_norm
        Übergangsentropie
        Schwungentropie
        aktive Intervallstatistik
        Hurst-R/S-Schätzung
    """

    rows = []

    T_over_values = np.linspace(T_over_2pi_min, T_over_2pi_max, steps)

    for idx, T_over in enumerate(T_over_values):
        T = T_over * 2.0 * math.pi

        df_tmp = simulate_zeitwuerfel(
            N=N,
            T=T,
            S=S,
            seed=seed,
            use_sharp_delta=True,
            beta=beta,
        )

        active_tmp = active_eabc_sequence(df_tmp)
        trans_tmp = active_transition_schwung(active_tmp)

        R8 = rotation_coherence_8(trans_tmp)
        R10 = rotation_coherence_10(trans_tmp)
        D = R10 - R8
        denom = R10 + R8
        O = (R10 - R8) / denom if denom != 0 else np.nan

        H_trans = transition_entropy(trans_tmp)
        H_schwung = schwung_entropy(trans_tmp)

        stats = active_interval_stats(active_tmp)

        active_dt = active_tmp["Delta_t_active"].dropna().values
        if len(active_dt) > 0:
            active_fluct = active_dt / np.mean(active_dt) - 1.0
            H_hurst = estimate_hurst_rs(active_fluct)
        else:
            H_hurst = np.nan

        rows.append({
            "idx": idx,
            "T": T,
            "T_over_2pi": T_over,
            "Delta_eff": df_tmp["Delta_eff"].iloc[0],
            "sigma_t": df_tmp["sigma_t"].iloc[0],
            "N_active": len(active_tmp),
            "N_trans": len(trans_tmp),
            "R8": R8,
            "R10": R10,
            "D_order": D,
            "O_order_norm": O,
            "H_transition": H_trans,
            "H_schwung": H_schwung,
            "H_hurst_rs": H_hurst,
            **stats,
        })

    scan = pd.DataFrame(rows)

    out_csv = f"{save_prefix}.csv"
    scan.to_csv(out_csv, index=False)

    print("\n=== Kritischer Hochauflösungs-Scan geschrieben ===")
    print(out_csv)

    idx_crit = scan["D_order"].abs().idxmin()
    crit = scan.loc[idx_crit]

    print("\n=== Kandidat kritischer Punkt ===")
    print("T/(2pi)    :", crit["T_over_2pi"])
    print("T          :", crit["T"])
    print("Delta_eff  :", crit["Delta_eff"])
    print("R8         :", crit["R8"])
    print("R10        :", crit["R10"])
    print("D_order    :", crit["D_order"])
    print("O_order    :", crit["O_order_norm"])
    print("H_trans    :", crit["H_transition"])
    print("H_schwung  :", crit["H_schwung"])
    print("H_hurst_rs :", crit["H_hurst_rs"])

    return scan


def centered_quaternionic_spectrum(
    df,
    Omega_values,
    active_only=True,
    dt_column="Delta_t",
):
    """
    Zentriertes quaternionisches Spektrum (DC der Kanäle q_e…q_c entfernt).

    active_only=True: nur Zeilen mit Label E, A, B, C.

    dt_column:
        \"Delta_t\" — ursprüngliches Simulationsintervall pro Zeile (wie bei quaternionic_spectrum).
        \"Delta_t_active\" — Abstände zwischen aufeinanderfolgenden aktiven t_k
        (Spalte aus active_eabc_sequence); erste Zeile mit NaN wird verworfen.

    Phasengewichtung: exp(-i * Omega * Delta_t) mit der gewählten Intervallspalte.
    """
    if active_only:
        data = df[df["label"].isin(["E", "A", "B", "C"])].copy()
    else:
        data = df.copy()

    if dt_column not in data.columns:
        raise ValueError(f'Spalte "{dt_column}" fehlt im DataFrame.')

    data = data.dropna(subset=[dt_column])
    if data.empty:
        raise ValueError(
            "centered_quaternionic_spectrum: keine Zeilen nach Filter/Dropna — "
            "prüfen Sie active_only und dt_column."
        )

    Delta_t = data[dt_column].values

    q_e = data["q_e"].values - np.mean(data["q_e"].values)
    q_a = data["q_a"].values - np.mean(data["q_a"].values)
    q_b = data["q_b"].values - np.mean(data["q_b"].values)
    q_c = data["q_c"].values - np.mean(data["q_c"].values)

    rows = []

    for Om in Omega_values:
        phase = np.exp(-1j * Om * Delta_t)

        QE = np.sum(q_e * phase)
        QA = np.sum(q_a * phase)
        QB = np.sum(q_b * phase)
        QC = np.sum(q_c * phase)

        P_E = abs(QE) ** 2
        P_A = abs(QA) ** 2
        P_B = abs(QB) ** 2
        P_C = abs(QC) ** 2

        rows.append({
            "Omega": Om,
            "P_E_centered": P_E,
            "P_A_centered": P_A,
            "P_B_centered": P_B,
            "P_C_centered": P_C,
            "P_total_centered": P_E + P_A + P_B + P_C,
        })

    return pd.DataFrame(rows)


# ------------------------------------------------------------
# 6. Riemann-Vergleichsmodul
# ------------------------------------------------------------

def normalize_sequence(x):
    """
    Normiert eine Folge auf Mittelwert 1.
    Nützlich für Abstandsvergleiche.
    """
    x = np.asarray(x, dtype=float)
    m = np.mean(x)
    if m == 0:
        return x
    return x / m


def spectral_rigidity_proxy(intervals, max_lag=50):
    """
    Einfacher Proxy für spektrale Steifigkeit:
    Autokorrelation der normierten Intervallfluktuationen.

        x_k = Delta_t_k / mean(Delta_t) - 1

    Ausgabe:
        DataFrame mit Lag und Autokorrelation.
    """
    intervals = normalize_sequence(intervals)
    x = intervals - 1.0

    rows = []
    N = len(x)

    for lag in range(1, max_lag + 1):
        if lag >= N:
            break

        a = x[:-lag]
        b = x[lag:]

        denom = np.std(a) * np.std(b)

        if denom == 0:
            corr = 0.0
        else:
            corr = np.mean((a - np.mean(a)) * (b - np.mean(b))) / denom

        rows.append({
            "lag": lag,
            "corr": corr
        })

    return pd.DataFrame(rows)


def compare_interval_statistics(intervals_a, intervals_b):
    """
    Vergleicht zwei Intervallfolgen statistisch.
    Gedacht für:
        Zeitwürfel-Intervalle vs. Riemann-Nullstellenabstände
        Zeitwürfel-Intervalle vs. Primlücken
        Zeitwürfel-Intervalle vs. Primvierlingsmuster
    """

    a = normalize_sequence(intervals_a)
    b = normalize_sequence(intervals_b)

    n = min(len(a), len(b))
    a = a[:n]
    b = b[:n]

    diff = a - b

    return {
        "n": n,
        "mean_a": float(np.mean(a)),
        "mean_b": float(np.mean(b)),
        "std_a": float(np.std(a)),
        "std_b": float(np.std(b)),
        "mean_abs_diff": float(np.mean(np.abs(diff))),
        "rmse": float(np.sqrt(np.mean(diff**2))),
        "corr": float(np.corrcoef(a, b)[0, 1]) if n > 2 else np.nan
    }


# ------------------------------------------------------------
# 7. Hauptlauf
# ------------------------------------------------------------

if __name__ == "__main__":

    # bisheriger Zeitwürfel-Lauf
    df = simulate_zeitwuerfel(
        N=2000,
        T=2.0 * math.pi,
        S=3.0,
        seed=123,
        use_sharp_delta=True,
        beta=3.0,
    )

    print("\n=== Bamberger Zeitwürfel: Kopf ===")
    print(df.head(20))

    print("\n=== Parameter ===")
    print("Delta_eff =", df["Delta_eff"].iloc[0])
    print("sigma_t   =", df["sigma_t"].iloc[0])

    print("\n=== EABC-Verteilung ===")
    print(df["label"].value_counts())

    # Spektrum (Omega=0 vermeiden → gleiches Raster auch für zentriertes Spektrum)
    Omega_values = np.linspace(0.001, 2.0 * math.pi, 512)
    spec = quaternionic_spectrum(df, Omega_values)

    print("\n=== Spektrum: Kopf ===")
    print(spec.head())

    active = active_eabc_sequence(df)
    trans = active_transition_schwung(active)

    spec_centered = centered_quaternionic_spectrum(df, Omega_values, active_only=True)

    active.to_csv("bamberger_zeitwuerfel_active_events.csv", index=False)
    trans.to_csv("bamberger_zeitwuerfel_active_transitions.csv", index=False)
    spec_centered.to_csv("bamberger_zeitwuerfel_centered_spectrum.csv", index=False)

    print("\nAktive EABC-Ereignisse:", len(active))
    print("\nAktive Übergänge:")
    print(trans.head(20))

    print_transition_report(trans)

    print("\nZentriertes Spektrum:")
    print(spec_centered.head())

    # Steifigkeitsproxy
    rigidity = spectral_rigidity_proxy(df["Delta_t"].values, max_lag=50)

    print("\n=== Spektrale Steifigkeit Proxy ===")
    print(rigidity.head(20))

    # CSV-Ausgabe
    df.to_csv("bamberger_zeitwuerfel_events.csv", index=False)
    spec.to_csv("bamberger_zeitwuerfel_spectrum.csv", index=False)
    rigidity.to_csv("bamberger_zeitwuerfel_rigidity.csv", index=False)

    print("\nDateien geschrieben:")
    print("  bamberger_zeitwuerfel_events.csv")
    print("  bamberger_zeitwuerfel_spectrum.csv")
    print("  bamberger_zeitwuerfel_rigidity.csv")
    print("  bamberger_zeitwuerfel_active_events.csv")
    print("  bamberger_zeitwuerfel_active_transitions.csv")
    print("  bamberger_zeitwuerfel_centered_spectrum.csv")

    # Temperatur-Kontrast
    contrast = run_temperature_contrast()
    print("\nZusätzlich (Temperatur-Kontrast):")
    print("  bamberger_zeitwuerfel_temperature_contrast.csv")
    print("  bamberger_zeitwuerfel_transition_matrix_Low-T.csv")
    print("  bamberger_zeitwuerfel_transition_matrix_High-T.csv")
    print("  bamberger_zeitwuerfel_schwung_distribution_Low-T.csv")
    print("  bamberger_zeitwuerfel_schwung_distribution_High-T.csv")

    # Kritischer Scan (einfach)
    critical_scan = run_simple_critical_scan(
        T_over_min=0.500,
        T_over_max=0.580,
        steps=161,
        S=3.0,
        N=3000,
        seed=123,
        beta=3.0,
    )
    print("\nZusätzlich (einfacher kritischer Scan):")
    print("  bamberger_zeitwuerfel_simple_critical_scan.csv")

    critical_scan_fine = run_simple_critical_scan(
        T_over_min=0.490,
        T_over_max=0.520,
        steps=601,
        S=3.0,
        N=8000,
        seed=123,
        beta=3.0,
        csv_path="bamberger_zeitwuerfel_simple_critical_scan_fine.csv",
    )
    print("\nZusätzlich (feiner kritischer Scan):")
    print("  bamberger_zeitwuerfel_simple_critical_scan_fine.csv")