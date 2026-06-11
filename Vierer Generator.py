import csv
from math import isqrt
from pathlib import Path

import numpy as np

try:
    import tkinter as tk
    from tkinter import messagebox, simpledialog
except Exception:
    tk = None
    messagebox = None
    simpledialog = None

# -----------------------------------------
# Schnelles Sieb (bool array)
# -----------------------------------------
def sieve(n):
    sieve = np.ones(n+1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(n**0.5)+1):
        if sieve[i]:
            sieve[i*i::i] = False
    return sieve

# -----------------------------------------
# Vierlinge effizient
# -----------------------------------------
def find_quadruplets(mask):
    idx = np.where(mask)[0]
    s = set(idx)
    return [p for p in idx if p+2 in s and p+6 in s and p+8 in s]


def format_quadruplet(p):
    return f"{p:>12} {p + 2:>12} {p + 6:>12} {p + 8:>12}"


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    limit = isqrt(n)
    divisor = 3
    while divisor <= limit:
        if n % divisor == 0:
            return False
        divisor += 2
    return True


def validate_quadruplet_csv(path):
    valid_rows = []
    errors = []

    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for line_number, row in enumerate(reader, start=2):
            try:
                p = int(row["p"])
                p2 = int(row["p+2"])
                p6 = int(row["p+6"])
                p8 = int(row["p+8"])
            except Exception as exc:
                errors.append((line_number, f"Ungueltige Zahlenwerte: {exc}"))
                continue

            actual = (p, p2, p6, p8)
            expected = (p, p + 2, p + 6, p + 8)
            if actual != expected:
                errors.append(
                    (
                        line_number,
                        f"Musterfehler: erwartet {expected}, gefunden {actual}",
                    )
                )
                continue

            non_primes = [value for value in actual if not is_prime(value)]
            if non_primes:
                errors.append((line_number, f"Nicht prim: {non_primes}"))
                continue

            valid_rows.append(actual)

    return valid_rows, errors


def choose_mode():
    if tk is not None and messagebox is not None:
        root = tk.Tk()
        root.withdraw()
        try:
            validate_csv = messagebox.askyesno(
                "Vierer Generator",
                "Soll eine vorhandene CSV-Datei auf korrekte Vierlinge geprueft werden?",
            )
        finally:
            root.destroy()
        return "validate_csv" if validate_csv else "generate"

    raw = input("Vorhandene CSV-Datei pruefen? [j/N]: ").strip().lower()
    return "validate_csv" if raw in {"j", "ja", "y", "yes"} else "generate"


def resolve_csv_path():
    base_dir = Path.cwd()
    candidates = sorted(
        base_dir.glob("vierlinge_bis_*.csv"),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
    if not candidates:
        raise SystemExit(
            "Im Startordner wurde keine Datei 'vierlinge_bis_*.csv' gefunden."
        )
    return candidates[0]


def ask_run_config(default_n=1_000_000):
    if tk is not None and simpledialog is not None:
        root = tk.Tk()
        root.withdraw()
        try:
            n = simpledialog.askinteger(
                "Vierer Generator",
                "Bitte Obergrenze N eingeben:",
                initialvalue=default_n,
                minvalue=11,
            )
            if n is None:
                raise SystemExit("Abgebrochen.")

            create_csv = False
            if messagebox is not None:
                create_csv = messagebox.askyesno(
                    "CSV-Ausgabe",
                    "Soll eine CSV-Datei mit allen gefundenen Vierlingen erstellt werden?",
                )
        finally:
            root.destroy()
        return n, create_csv

    raw = input(f"Bitte Obergrenze N eingeben [{default_n}]: ").strip()
    n = int(raw) if raw else default_n
    create_csv_raw = input("CSV-Datei erstellen? [j/N]: ").strip().lower()
    return n, create_csv_raw in {"j", "ja", "y", "yes"}


def write_quadruplets_csv(quads, path):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["p", "p+2", "p+6", "p+8"])
        for p in quads:
            writer.writerow([p, p + 2, p + 6, p + 8])


def show_validation_result(path, valid_rows, errors):
    print(f"CSV-Datei: {path}")
    print(f"Gepruefte Datenzeilen: {len(valid_rows) + len(errors)}")
    print(f"Fehlerhafte Zeilen: {len(errors)}")
    if valid_rows:
        print("Letzte vier gueltige Vierlinge:")
        print(f"{'p':>12} {'p+2':>12} {'p+6':>12} {'p+8':>12}")
        for quadruplet in valid_rows[-4:]:
            print(
                f"{quadruplet[0]:>12} {quadruplet[1]:>12} "
                f"{quadruplet[2]:>12} {quadruplet[3]:>12}"
            )
    if errors:
        print("Beispielhafte Fehler:")
        for line_number, message in errors[:10]:
            print(f"Zeile {line_number}: {message}")

    if tk is not None and messagebox is not None:
        root = tk.Tk()
        root.withdraw()
        try:
            summary = (
                f"Datei: {path}\n"
                f"Gepruefte Datenzeilen: {len(valid_rows) + len(errors)}\n"
                f"Fehlerhafte Zeilen: {len(errors)}"
            )
            if valid_rows:
                last_four = "\n".join(str(item) for item in valid_rows[-4:])
                summary += f"\n\nLetzte vier gueltige Vierlinge:\n{last_four}"
            if errors:
                sample_errors = "\n".join(
                    f"Zeile {line_number}: {message}"
                    for line_number, message in errors[:5]
                )
                summary += f"\n\nBeispielhafte Fehler:\n{sample_errors}"

            messagebox.showinfo("CSV-Pruefung", summary)
        finally:
            root.destroy()

# -----------------------------------------
# EABC
# -----------------------------------------
def classify(p):
    return {1:0,5:1,7:2,11:3}.get(p % 12, -1)

weights = np.array([0.0, 1.0, -1.0, 0.5])

# -----------------------------------------
# Potential (log dominant!)
# -----------------------------------------
def potential(p):
    return np.sin(np.log(p)) * weights[classify(p)]

# -----------------------------------------
# Matrix (fix oder aus Bestwerten einsetzen!)
# -----------------------------------------
def internal_matrix():
    return np.array([
        [0,1.0,0.8,0.5],
        [1.0,0,1.2,0.7],
        [0.8,1.2,0,0.9],
        [0.5,0.7,0.9,0]
    ])

# -----------------------------------------
# Operator
# -----------------------------------------
def build_operator(quads, eps):
    Nq = len(quads)
    dim = 4*Nq
    H = np.zeros((dim, dim))
    M = internal_matrix()
    
    for i,p in enumerate(quads):
        H[4*i:4*i+4, 4*i:4*i+4] = M
        
        for j,shift in enumerate([0,2,6,8]):
            H[4*i+j,4*i+j] += potential(p+shift)
    
    # Kopplung
    for i in range(Nq-1):
        for k in range(4):
            H[4*i+k,4*(i+1)+k] = 1
            H[4*(i+1)+k,4*i+k] = 1
    
    # Hybrid
    R = np.random.normal(size=(dim,dim))
    R = (R+R.T)/2
    H += eps * R
    
    return H

# -----------------------------------------
# KS-Distanz
# -----------------------------------------
def ks(a,b):
    a = np.sort(a)
    b = np.sort(b)
    x = np.sort(np.concatenate([a,b]))
    
    c1 = np.searchsorted(a,x,side='right')/len(a)
    c2 = np.searchsorted(b,x,side='right')/len(b)
    
    return np.max(np.abs(c1-c2))

# -----------------------------------------
# Zeta laden
# -----------------------------------------
z = np.load("zeros6.npz")
zeros = z[list(z.keys())[0]][:30000]
z_sp = np.diff(zeros)
z_sp /= np.mean(z_sp)

# -----------------------------------------
# HAUPTLAUF
# -----------------------------------------
mode = choose_mode()

if mode == "validate_csv":
    csv_path = resolve_csv_path()
    valid_rows, errors = validate_quadruplet_csv(csv_path)
    show_validation_result(csv_path, valid_rows, errors)
else:
    N, create_csv = ask_run_config()
    mask = sieve(N)
    quads = find_quadruplets(mask)

    print("Vierlinge:", len(quads))
    if quads:
        print("Letzte vier gefundene Vierlinge:")
        print(f"{'p':>12} {'p+2':>12} {'p+6':>12} {'p+8':>12}")
        for quad in quads[-4:]:
            print(format_quadruplet(quad))

    if create_csv:
        csv_path = f"vierlinge_bis_{N}.csv"
        write_quadruplets_csv(quads, csv_path)
        print(f"CSV gespeichert: {csv_path}")

    H = build_operator(quads, eps=0.6)

    ev = np.linalg.eigvalsh(H)
    sp = np.diff(np.sort(ev))
    sp /= np.mean(sp)

    print("Var:", np.var(sp))
    print("KS:", ks(sp, z_sp))