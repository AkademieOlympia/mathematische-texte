#!/usr/bin/env python3
"""EABC Hybrid-Shor Simulator (1000-Qubit Gottesman–Knill Tableau).

Claim-Wall:
  - STABILIZER TABLEAU PHASE   : SIMULATED
  - NON-CLIFFORD / ORACLE STEP : CLASSICAL / CLAIM-WALL
  - SHOR HARDWARE CLAIM        : NOT MADE

Release-Ledger-Referenz (N=1701):
  a=1502, r=54, a^{r/2} ≡ 1457 (mod 1701) — echte nichttriviale 1-Wurzel.
  Der Postprocessing-Wert 1216 war falsch (1216² ≡ 487 ≢ 1 mod 1701).
"""

from __future__ import annotations

import argparse
import math
import random
import sys
import warnings

import numpy as np

try:
    # SymPy ≥1.12: n_order; ältere Alias-Namen optional
    from sympy.ntheory import n_order as _sympy_n_order
except ImportError:  # pragma: no cover
    try:
        from sympy import multiplicative_order as _sympy_n_order  # type: ignore
    except ImportError:
        _sympy_n_order = None  # type: ignore[assignment]

# --- Release-Ledger (exakt, Lean: EabcRenorm.QEC.HybridShorCRTExample) ---
LEDGER_N = 1701
LEDGER_A = 1502
LEDGER_R = 54
LEDGER_HALF = 27
LEDGER_MODULAR_POWER = 1457  # korrigiert; nicht 1216
LEDGER_P = 7
LEDGER_Q = 243
REJECTED_POSTPROCESS = 1216

# --- BIG-SMOKE Ledger (allgemeiner Shor-Schnitt, nicht Semiprim) ---
BIG_N = 12_345_678_907
BIG_A = 11_036_832_490
BIG_R = 3_065_905_626
BIG_HALF_POW = 760_221_407
BIG_SHOR_D1 = 1559  # prim
BIG_SHOR_D2 = 7_918_973  # zusammengesetzt: 167 × 47419
BIG_PRIMES = (167, 1559, 47_419)

# --- 47-Bit Smoke: beide Shor-Teiler zusammengesetzt ---
BIT47_N = 123_456_789_012_345
BIT47_A = 51_935_289_148_391
BIT47_R = 58_565_992_650
BIT47_HALF_POW = 76_109_610_082_991
BIT47_SHOR_D1 = 5_449_165  # 5 × 283 × 3851
BIT47_SHOR_D2 = 22_656_093  # 3 × 7_552_031
BIT47_PRIMES = (3, 5, 283, 3851, 7_552_031)

# --- 134-Bit Smoke: großer zusammengesetzter Teiler mit 43² ---
BIT134_N = 12_345_678_901_234_567_890_123_456_789_012_345_678_901
BIT134_A = 5_414_613_428_567_653_958_871_976_718_907_210_683_804
BIT134_R = 4_331_382_856_009_368_439_019_253_117_910_730_368
BIT134_HALF_POW = 2_265_680_664_137_665_268_250_570_897_783_108_255_044
BIT134_SHOR_D1 = 55_202_177  # prim
BIT134_SHOR_D2 = 223_644_783_089_525_036_125_358_910_917_813
# 29 × 43² × 2928109491677 × 1424415039563189
BIT134_PRIMES = {
    29: 1,
    43: 2,
    2_928_109_491_677: 1,
    55_202_177: 1,
    1_424_415_039_563_189: 1,
}


class GottesmanKnill1000QC:
    """Symplektischer Stabilizer-Simulator für 1000 Qubits über 𝔽₂."""

    def __init__(self, num_qubits: int = 1000):
        self.n = num_qubits
        self.tableau = np.zeros((2 * self.n, 2 * self.n + 1), dtype=np.uint8)
        for i in range(self.n):
            self.tableau[i, i] = 1  # Destabilisator X_i
            self.tableau[self.n + i, self.n + i] = 1  # Stabilisator Z_i

    def reset_registers(self):
        self.tableau.fill(0)
        for i in range(self.n):
            self.tableau[i, i] = 1
            self.tableau[self.n + i, self.n + i] = 1

    def apply_h(self, q: int):
        x = self.tableau[:, q].copy()
        z = self.tableau[:, self.n + q].copy()
        self.tableau[:, 2 * self.n] ^= x & z
        self.tableau[:, q] = z
        self.tableau[:, self.n + q] = x

    def apply_cnot(self, control: int, target: int):
        x_c = self.tableau[:, control]
        z_c = self.tableau[:, self.n + control]
        x_t = self.tableau[:, target]
        z_t = self.tableau[:, self.n + target]
        self.tableau[:, 2 * self.n] ^= x_c & z_t & (x_t ^ z_c ^ 1)
        self.tableau[:, target] ^= x_c
        self.tableau[:, self.n + control] ^= z_t


class DialogShorConsole:
    """
    Dialoggesteuerte Konsole für Shors Faktorisierungsalgorithmus.
    Verbindet klassische Arithmetik mit dem 1000-Qubit Clifford-Simulator.
    """

    def __init__(self):
        self.qc = GottesmanKnill1000QC(num_qubits=1000)

    def print_header(self):
        print("=" * 65)
        print("  EABC HYBRID-SHOR SIMULATOR (1000 QUBITS)")
        print("  Claim-Wall: PhysicalForce := False | Non-Clifford Gate Protection")
        print("  Ledger-Korrektur: a^{r/2} ≡ 1457 (mod 1701), nicht 1216")
        print("=" * 65)

    @staticmethod
    def register_layout(N: int) -> tuple[int, int, int]:
        """Steuer- / Ziel- / Frei-Qubits wie im Release-Ledger."""
        target_bits = math.ceil(math.log2(N)) if N > 1 else 1
        control_bits = 2 * target_bits
        free_bits = 1000 - control_bits - target_bits
        return control_bits, target_bits, free_bits

    @staticmethod
    def crt_sign_pair(x: int, p: int, q: int) -> tuple[int, int]:
        """Vorzeichenpaar (x mod p, x mod q) im CRT-Produktring."""
        return x % p, x % q

    def factor_from_period(self, a: int, N: int, r: int) -> dict | None:
        """Klassisches Postprocessing; verlangt echte nichttriviale 1-Wurzel."""
        if r % 2 != 0:
            return None
        half_pow = pow(a, r // 2, N)
        if half_pow in (1, N - 1):
            return None
        if pow(half_pow, 2, N) != 1:
            raise RuntimeError(
                f"Postprocessing-Fehler: {half_pow}^2 ≢ 1 (mod {N}). "
                f"Bekannter Fehlwert: {REJECTED_POSTPROCESS}."
            )
        p = math.gcd(half_pow - 1, N)
        q = math.gcd(half_pow + 1, N)
        return {
            "half_pow": half_pow,
            "p": p,
            "q": q,
            "ok": p * q == N and 1 < p < N and 1 < q < N,
        }

    def print_crt_audit(self, half_pow: int, N: int, p: int, q: int) -> None:
        s_p, s_q = self.crt_sign_pair(half_pow, p, q)
        # Normiere auf ±1-Darstellung, falls möglich
        label_p = "+1" if s_p == 1 else ("-1" if s_p == p - 1 else str(s_p))
        label_q = "+1" if s_q == 1 else ("-1" if s_q == q - 1 else str(s_q))
        print("\n  [CRT-Audit]")
        print(f"  {half_pow} ≡ {label_p} (mod {p}),  {half_pow} ≡ {label_q} (mod {q})")
        print(f"  => ({label_p})^2 × ({label_q})^2 = (+1,+1) ≡ 1 (mod {N})")
        print(f"  ggT({half_pow}-1, {N}) = {math.gcd(half_pow - 1, N)}")
        print(f"  ggT({half_pow}+1, {N}) = {math.gcd(half_pow + 1, N)}")
        print(
            f"  Verworfener Fehlwert {REJECTED_POSTPROCESS}: "
            f"{REJECTED_POSTPROCESS}^2 ≡ {pow(REJECTED_POSTPROCESS, 2, N)} ≢ 1 (mod {N})"
        )

    @staticmethod
    def _is_prime(n: int) -> bool:
        try:
            from sympy import isprime

            return bool(isprime(n))
        except ImportError:
            if n < 2:
                return False
            if n % 2 == 0:
                return n == 2
            d = 3
            while d * d <= n:
                if n % d == 0:
                    return False
                d += 2
            return True

    @staticmethod
    def _complete_prime_factors(n: int) -> dict[int, int]:
        """Klassische Vervollständigung nach dem ersten Shor-Schnitt (kein Hardware-Claim)."""
        try:
            from sympy import factorint

            return {int(p): int(e) for p, e in factorint(n).items()}
        except ImportError as exc:
            raise RuntimeError("SymPy factorint benötigt für Primfaktor-Vervollständigung") from exc

    @staticmethod
    def mainstream_factorint(n: int) -> dict[int, int]:
        """Mainstream-Referenz: SymPy `factorint` (Pollard/ECM/… — rein klassisch)."""
        return DialogShorConsole._complete_prime_factors(n)

    @staticmethod
    def _format_factors(fac: dict[int, int]) -> str:
        return " × ".join(str(p) if e == 1 else f"{p}^{e}" for p, e in sorted(fac.items()))

    @staticmethod
    def resolve_full_prime_factorization(d1: int, d2: int) -> tuple[str, dict[int, int], dict[int, int]]:
        """Löst beide Shor-Schnitt-Teiler rekursiv in echte Primfaktoren auf."""
        factors_d1 = DialogShorConsole._complete_prime_factors(d1)
        factors_d2 = DialogShorConsole._complete_prime_factors(d2)
        full_factors = DialogShorConsole._complete_prime_factors(d1 * d2)
        prime_list: list[str] = []
        for p, exp in sorted(full_factors.items()):
            prime_list.extend([str(p)] * exp)
        return " × ".join(prime_list), factors_d1, factors_d2

    def hybrid_shor_factorize(
        self, N: int, a: int | None = None, *, max_a_trials: int = 32
    ) -> dict:
        """Hybrid-Shor-Pipeline: Tableau-Ritual + klassisches Order-Finding + Postprocessing.

        Liefert ersten Schnitt und vollständige Primfaktoren (Vervollständigung klassisch).
        Kein Quantum-/Hardware-Claim.
        """
        import time

        if N < 3 or N % 2 == 0:
            raise ValueError(f"N={N} ungeeignet (gerade oder zu klein)")

        rng = random.Random(1701 ^ N)
        candidates: list[int] = []
        if a is not None:
            candidates.append(a)
        # Ledger-Basen zuerst, dann Zufall
        if N == LEDGER_N:
            candidates.append(LEDGER_A)
        if N == BIG_N:
            candidates.append(BIG_A)
        if N == BIT47_N:
            candidates.append(BIT47_A)
        if N == BIT134_N:
            candidates.append(BIT134_A)
        while len(candidates) < max_a_trials:
            candidates.append(rng.randint(2, N - 1))

        control_bits, target_bits, _ = self.register_layout(N)
        self.qc.reset_registers()
        for q in range(control_bits):
            self.qc.apply_h(q)

        last_err = "keine geeignete Basis gefunden"
        for trial_a in candidates:
            if math.gcd(trial_a, N) != 1:
                g = math.gcd(trial_a, N)
                # Glücksfall: klassischer Teiler ohne Periodensuche
                fac = self._complete_prime_factors(N)
                return {
                    "N": N,
                    "a": trial_a,
                    "method": "gcd_lucky",
                    "r": None,
                    "order_ms": 0.0,
                    "first_cut": tuple(sorted((g, N // g))),
                    "primes": fac,
                    "control_qubits": control_bits,
                    "target_qubits": target_bits,
                }
            t0 = time.perf_counter()
            try:
                r = self.find_period_classical_oracle(trial_a, N)
            except Exception as exc:  # noqa: BLE001
                last_err = str(exc)
                continue
            order_ms = (time.perf_counter() - t0) * 1000
            try:
                cut = self.factor_from_period(trial_a, N, r)
            except RuntimeError as exc:
                last_err = str(exc)
                continue
            if cut is None:
                last_err = f"trivialer Fall für a={trial_a}, r={r}"
                continue
            d1, d2 = sorted((cut["p"], cut["q"]))
            full_str, fac_d1, fac_d2 = self.resolve_full_prime_factorization(d1, d2)
            primes = self._complete_prime_factors(d1 * d2)
            return {
                "N": N,
                "a": trial_a,
                "method": "hybrid_shor",
                "r": r,
                "order_ms": order_ms,
                "half_pow": cut["half_pow"],
                "first_cut": (d1, d2),
                "factors_d1": fac_d1,
                "factors_d2": fac_d2,
                "full_prime_str": full_str,
                "primes": primes,
                "control_qubits": control_bits,
                "target_qubits": target_bits,
            }
        raise RuntimeError(f"Hybrid-Shor-Faktorisierung fehlgeschlagen für N={N}: {last_err}")

    def compare_with_mainstream(
        self, N: int, a: int | None = None
    ) -> dict:
        """Vergleich Hybrid-Shor vs. Mainstream-`factorint` (beide klassisch)."""
        import time

        hybrid = self.hybrid_shor_factorize(N, a)
        t0 = time.perf_counter()
        mainstream = self.mainstream_factorint(N)
        mainstream_ms = (time.perf_counter() - t0) * 1000
        match = hybrid["primes"] == mainstream
        prod_h = 1
        for p, e in hybrid["primes"].items():
            prod_h *= p**e
        prod_m = 1
        for p, e in mainstream.items():
            prod_m *= p**e
        return {
            "N": N,
            "hybrid": hybrid,
            "mainstream": mainstream,
            "mainstream_ms": mainstream_ms,
            "primes_match": match,
            "product_ok": prod_h == N and prod_m == N,
            "claim_wall": {
                "stabilizer": "SIMULATED",
                "order_finding": "CLASSICAL (n_order)",
                "mainstream": "CLASSICAL (sympy.factorint)",
                "hardware": "NOT CLAIMED",
                "quantum_speedup": "NOT CLAIMED",
            },
        }

    def run_compare(self, numbers: list[int] | None = None) -> int:
        """CLI-Vergleichstest Hybrid-Shor vs. Mainstream-Faktorisierer."""
        self.print_header()
        print("\n[COMPARE] Hybrid-Shor-Postprocessing vs. SymPy factorint")
        print("  Beide Pfade sind klassisch — kein Quanten-Speedup-Claim.")
        if numbers is None:
            numbers = [LEDGER_N, BIG_N, BIT47_N, BIT134_N]
        all_ok = True
        for N in numbers:
            print("\n" + "-" * 65)
            print(f"  N = {N}  ({N.bit_length()} Bit)")
            try:
                report = self.compare_with_mainstream(N)
            except Exception as exc:  # noqa: BLE001
                print(f"  FEHLER: {exc}")
                all_ok = False
                continue
            h = report["hybrid"]
            print(f"  Hybrid-Methode:     {h['method']}")
            print(f"  Basis a:            {h['a']}")
            if h["r"] is not None:
                print(f"  Ordnung r:          {h['r']}  (n_order: {h['order_ms']:.3f} ms)")
                print(f"  Erster Shor-Schnitt:{h['first_cut'][0]} × {h['first_cut'][1]}")
            print(f"  Hybrid-Primfaktoren:{self._format_factors(h['primes'])}")
            print(
                f"  Mainstream factorint:{self._format_factors(report['mainstream'])}  "
                f"({report['mainstream_ms']:.3f} ms)"
            )
            print(f"  Primfaktoren gleich: {report['primes_match']}")
            print(f"  Produkt = N:         {report['product_ok']}")
            print(
                "  Claim-Wall: "
                f"{report['claim_wall']['order_finding']} | "
                f"{report['claim_wall']['mainstream']} | "
                f"{report['claim_wall']['hardware']}"
            )
            if not (report["primes_match"] and report["product_ok"]):
                all_ok = False
        print("\n" + "=" * 65)
        print(f"  COMPARE RESULT: {'PASS' if all_ok else 'FAIL'}")
        print("=" * 65)
        return 0 if all_ok else 1

    def run_big_smoke(self) -> int:
        """34-Bit-Smoke: allgemeiner Shor-Schnitt (Teiler ≠ automatisch Primzahlen)."""
        import time

        self.print_header()
        N, a = BIG_N, BIG_A
        print(f"\n[BIG-SMOKE] N = {N} ({N.bit_length()} Bit), a = {a}")
        print("  Didaktik: Shor liefert echte Teiler d1·d2=N — nicht zwingend Primfaktoren.")
        if _sympy_n_order is None:
            print("FEHLER: SymPy nicht verfügbar — Orakel würde O(r) hängen.")
            return 1
        if math.gcd(a, N) != 1:
            print(f"gcd(a,N)={math.gcd(a, N)} > 1 — wähle andere Basis.")
            return 1

        control_bits, target_bits, free_bits = self.register_layout(N)
        active = control_bits + target_bits
        print(
            f"Register: Steuer={control_bits}, Ziel={target_bits}, "
            f"aktiv={active}/1000, frei={free_bits}"
        )
        print(
            f"Tableau: Kapazität {2 * self.qc.n} × {2 * self.qc.n + 1} über 𝔽₂; "
            f"aktive Skala {2 * active} × {2 * active + 1}"
        )
        self.qc.reset_registers()
        for q in range(control_bits):
            self.qc.apply_h(q)
        print("Tableau: Hadamard auf Steuerregister (SIMULATED).")
        print("ORDER FINDING: CLASSICALLY EVALUATED (n_order) — kein Quantenorakel-Claim.")

        t0 = time.perf_counter()
        r = self.find_period_classical_oracle(a, N)
        dt_ms = (time.perf_counter() - t0) * 1000
        print(f"-> r = {r}  (lokale n_order-Wandzeit: {dt_ms:.3f} ms)")
        print(
            "   Hinweis: ms-Angabe = klassische SymPy-CPU-Zeit auf dieser Instanz, "
            "kein Komplexitätsmaß eines Quantenorakels."
        )
        if r != BIG_R:
            print(f"WARNUNG: r={r} ≠ Ledger-Soll {BIG_R}")

        try:
            result = self.factor_from_period(a, N, r)
        except RuntimeError as exc:
            print(f"  [!] {exc}")
            return 1
        if result is None:
            print("  [-] Trivialer Fall (ungerade r oder ±1). Kein Faktor aus diesem a.")
            return 1

        half_pow = result["half_pow"]
        d1, d2 = result["p"], result["q"]
        # kanonische Reihenfolge wie im Ledger
        if d1 > d2:
            d1, d2 = d2, d1

        print(f"\n  [ERSTER SHOR-SCHNITT]")
        print(f"  a^(r/2) ≡ {half_pow} (mod {N})")
        print(f"  gcd(x-1, N) = {d1}   prim? {self._is_prime(d1)}")
        print(f"  gcd(x+1, N) = {d2}   prim? {self._is_prime(d2)}")
        print(f"  Zerlegung: {d1} × {d2} = {d1 * d2}  (ok: {d1 * d2 == N})")
        if half_pow != BIG_HALF_POW or {d1, d2} != {BIG_SHOR_D1, BIG_SHOR_D2}:
            print("  WARNUNG: weicht vom BIG-SMOKE-Ledger ab.")

        print(f"\n  [ZWEITER SCHRITT — rekursive Primfaktor-Auflösung]")
        print("  Klassisch: beide Shor-Teiler via factorint auflösen (kein Hardware-Claim).")
        full, fac_d1, fac_d2 = self.resolve_full_prime_factorization(d1, d2)
        print(f"  d1={d1}: {self._format_factors(fac_d1)}  (prim? {self._is_prime(d1)})")
        print(f"  d2={d2}: {self._format_factors(fac_d2)}  (prim? {self._is_prime(d2)})")
        primes = self._complete_prime_factors(d1 * d2)
        prod = 1
        for p, e in primes.items():
            prod *= p**e
        expected = set(BIG_PRIMES)
        got = set(primes)
        print(f"\n  VOLLSTÄNDIGE PRIMFAKTORZERLEGUNG: {full}")
        print(f"  Verifikation Produkt = N: {prod == N}")
        print(f"  Ledger-Soll {{167,1559,47419}}: {got == expected}")
        print("  STABILIZER PHASE: SIMULATED")
        print("  PHYSICAL HARDWARE CLAIM: NOT MADE")
        return 0 if result["ok"] and prod == N and got == expected else 1

    def run_bit47_smoke(self) -> int:
        """47-Bit-Smoke: beide euklidischen Shor-Teiler sind zusammengesetzt."""
        import time

        self.print_header()
        N, a = BIT47_N, BIT47_A
        print(f"\n[47-BIT SMOKE] N = {N} ({N.bit_length()} Bit), a = {a}")
        print("  Didaktik: erster Schnitt spaltet auf BEIDEN Seiten in zusammengesetzte Teiler.")
        if _sympy_n_order is None:
            print("FEHLER: SymPy nicht verfügbar — Orakel würde O(r) hängen.")
            return 1
        if math.gcd(a, N) != 1:
            print(f"gcd(a,N)={math.gcd(a, N)} > 1 — wähle andere Basis.")
            return 1

        control_bits, target_bits, free_bits = self.register_layout(N)
        active = control_bits + target_bits
        print(
            f"Register: Steuer={control_bits}, Ziel={target_bits}, "
            f"aktiv={active}/1000, frei={free_bits}"
        )
        print(
            f"Tableau: Kapazität {2 * self.qc.n} × {2 * self.qc.n + 1} über 𝔽₂; "
            f"aktive Skala {2 * active} × {2 * active + 1}"
        )
        self.qc.reset_registers()
        for q in range(control_bits):
            self.qc.apply_h(q)
        print("Tableau: Hadamard auf Steuerregister (SIMULATED).")
        print("ORDER FINDING: CLASSICAL ORACLE / n_order — kein Quantenorakel-Claim.")

        t0 = time.perf_counter()
        r = self.find_period_classical_oracle(a, N)
        dt_ms = (time.perf_counter() - t0) * 1000
        print(f"-> r = {r}  (lokale n_order-Wandzeit: {dt_ms:.3f} ms)")
        print(
            "   Hinweis: ms-Angabe = klassische SymPy-CPU-Zeit auf dieser Instanz, "
            "kein Komplexitätsmaß eines Quantenorakels."
        )
        if r != BIT47_R:
            print(f"WARNUNG: r={r} ≠ Ledger-Soll {BIT47_R}")

        try:
            result = self.factor_from_period(a, N, r)
        except RuntimeError as exc:
            print(f"  [!] {exc}")
            return 1
        if result is None:
            print("  [-] Trivialer Fall (ungerade r oder ±1). Kein Faktor aus diesem a.")
            return 1

        half_pow = result["half_pow"]
        d1, d2 = sorted((result["p"], result["q"]))

        print(f"\n  [ERSTER SHOR-SCHNITT]")
        print(f"  a^(r/2) ≡ {half_pow} (mod {N})")
        print(f"  gcd(x-1, N) = {d1}   prim? {self._is_prime(d1)}")
        print(f"  gcd(x+1, N) = {d2}   prim? {self._is_prime(d2)}")
        print(f"  Zerlegung: {d1} × {d2} = {d1 * d2}  (ok: {d1 * d2 == N})")
        if half_pow != BIT47_HALF_POW or {d1, d2} != {BIT47_SHOR_D1, BIT47_SHOR_D2}:
            print("  WARNUNG: weicht vom 47-Bit-Ledger ab.")
        if self._is_prime(d1) or self._is_prime(d2):
            print("  WARNUNG: erwartet waren ZWEI zusammengesetzte Teiler.")

        print(f"\n  [REKURSIVE AUFLÖSUNG BEIDER TEILER]")
        full, fac_d1, fac_d2 = self.resolve_full_prime_factorization(d1, d2)
        print(f"  {d1} = {self._format_factors(fac_d1)}")
        print(f"  {d2} = {self._format_factors(fac_d2)}")
        primes = self._complete_prime_factors(d1 * d2)
        prod = 1
        for p, e in primes.items():
            prod *= p**e
        expected = set(BIT47_PRIMES)
        got = set(primes)
        print(f"\n  VOLLSTÄNDIGE PRIMFAKTORZERLEGUNG: {full}")
        print(f"  Verifikation Produkt = N: {prod == N}")
        print(f"  Ledger-Soll {sorted(expected)}: {got == expected}")
        print("  STABILIZER TABLEAU: SIMULATED")
        print("  HARDWARE CLAIM: NOT MADE")
        print("  QUANTUM SPEEDUP CLAIM: NOT MADE")
        return 0 if result["ok"] and prod == N and got == expected else 1

    def run_bit134_smoke(self) -> int:
        """134-Bit-Smoke: großer zusammengesetzter Teiler inkl. 43²-Korrektur."""
        import time

        self.print_header()
        N, a = BIT134_N, BIT134_A
        print(f"\n[134-BIT SMOKE] N = {N} ({N.bit_length()} Bit), a = {a}")
        print("  Korrektur: großer Teiler = 29 × 43² × 2928109491677 × 1424415039563189")
        if _sympy_n_order is None:
            print("FEHLER: SymPy nicht verfügbar — Orakel würde O(r) hängen.")
            return 1
        if math.gcd(a, N) != 1:
            print(f"gcd(a,N)={math.gcd(a, N)} > 1 — wähle andere Basis.")
            return 1

        control_bits, target_bits, free_bits = self.register_layout(N)
        active = control_bits + target_bits
        print(
            f"Register: Steuer={control_bits}, Ziel={target_bits}, "
            f"aktiv={active}/1000, frei={free_bits}"
        )
        print(
            f"Tableau: Kapazität {2 * self.qc.n} × {2 * self.qc.n + 1} über 𝔽₂; "
            f"aktive Skala {2 * active} × {2 * active + 1}"
        )
        self.qc.reset_registers()
        for q in range(control_bits):
            self.qc.apply_h(q)
        print("Tableau: Hadamard auf Steuerregister (SIMULATED / Gottesman-Knill 𝔽₂-Scaling).")
        print("ORDER FINDING: CLASSICAL ORACLE / n_order — kein Quantenorakel-Claim.")

        t0 = time.perf_counter()
        r = self.find_period_classical_oracle(a, N)
        dt_ms = (time.perf_counter() - t0) * 1000
        print(f"-> r = {r}  (lokale n_order-Wandzeit: {dt_ms:.3f} ms)")
        print(
            "   Hinweis: ms-Angabe = klassische SymPy-CPU-Zeit auf dieser Instanz, "
            "kein Komplexitätsmaß eines Quantenorakels."
        )
        if r != BIT134_R:
            print(f"WARNUNG: r={r} ≠ Ledger-Soll {BIT134_R}")
        if pow(a, r, N) != 1:
            print("FEHLER: a^r ≢ 1 (mod N)")
            return 1

        try:
            result = self.factor_from_period(a, N, r)
        except RuntimeError as exc:
            print(f"  [!] {exc}")
            return 1
        if result is None:
            print("  [-] Trivialer Fall (ungerade r oder ±1). Kein Faktor aus diesem a.")
            return 1

        half_pow = result["half_pow"]
        d1, d2 = sorted((result["p"], result["q"]))

        print(f"\n  [ERSTER SHOR-SCHNITT]")
        print(f"  a^(r/2) ≡ {half_pow} (mod {N})")
        print(f"  gcd(x-1, N) / gcd(x+1, N) sortiert: {d1} , {d2}")
        print(f"  {d1} prim? {self._is_prime(d1)}")
        print(f"  {d2} prim? {self._is_prime(d2)}")
        print(f"  PRODUCT CHECK: {d1} × {d2} = {d1 * d2}  (ok: {d1 * d2 == N})")
        if half_pow != BIT134_HALF_POW or {d1, d2} != {BIT134_SHOR_D1, BIT134_SHOR_D2}:
            print("  WARNUNG: weicht vom 134-Bit-Ledger ab.")

        print(f"\n  [REKURSIVE AUFLÖSUNG]")
        full, fac_d1, fac_d2 = self.resolve_full_prime_factorization(d1, d2)
        print(f"  {d1} = {self._format_factors(fac_d1)}")
        print(f"  {d2} = {self._format_factors(fac_d2)}")
        primes = self._complete_prime_factors(d1 * d2)
        prod = 1
        for p, e in primes.items():
            prod *= p**e
        print(f"\n  VOLLSTÄNDIGE PRIMFAKTORZERLEGUNG: {full}")
        print(f"  Verifikation Produkt = N: {prod == N}")
        print(f"  Ledger-Map Match: {primes == BIT134_PRIMES}")
        print("  STABILIZER TABLEAU: SIMULATED (Gottesman-Knill 𝔽₂-Scaling)")
        print("  ORDER FINDING: CLASSICAL ORACLE / n_order")
        print("  PHYSICAL HARDWARE CLAIM: NOT MADE")
        print("  QUANTUM SPEEDUP CLAIM: NOT MADE")
        print("  CLAIM-WALL: PhysicalForce := False")
        return 0 if result["ok"] and prod == N and primes == BIT134_PRIMES else 1

    def run_ledger_demo(self, interactive: bool = False) -> int:
        """Fester Release-Lauf N=1701, a=1502 → 1457 → 7×243."""
        self.print_header()
        N, a = LEDGER_N, LEDGER_A
        print(f"\n[DEMO / RELEASE-LEDGER] N = {N} = {LEDGER_P} × {LEDGER_Q}, a = {a}")

        if math.gcd(a, N) != 1:
            print("Ledger-Inkonsistenz: gcd(a,N) ≠ 1")
            return 1

        control_bits, target_bits, free_bits = self.register_layout(N)
        active = control_bits + target_bits
        print(f"\n[Schritt 3] Initialisiere 1000-Qubit Quanten-Register...")
        print(f"  - Steuerregister (Superposition/QFT): {control_bits} Qubits")
        print(f"  - Zielregister (Modul-Exp):           {target_bits} Qubits")
        print(f"  - Aktive Tableau-Qubits:               {active} / 1000")
        print(f"  - Verbleibendes Freiregister:          {free_bits} Qubits")

        self.qc.reset_registers()
        print("\n-> Wende Walsh-Hadamard-Transformation auf Steuerregister an...")
        for q in range(control_bits):
            self.qc.apply_h(q)
        print(f"   [OK] {control_bits} Steuer-Qubits in Superposition (SIMULATED).")

        print("\n" + "-" * 65)
        print("  [EPISTEMISCHE SCHRANKE / EABC CLAIM-WALL]")
        print("  Tableau: SIMULATED | Oracle/QFT: CLASSICAL | Hardware: NOT MADE")
        print("-" * 65)

        if interactive:
            proceed = input(
                "\nSoll die Periodensuche r = ord_N(a) klassisch evaluiert werden? (j/n): "
            ).strip().lower()
            if proceed != "j":
                print("Dialog durch Anwender abgebrochen.")
                return 0

        print(f"\n[Schritt 4] Klassische Periodensuche (Oracle / CLAIM-WALL)...")
        r = self.find_period_classical_oracle(a, N)
        print(f"-> Gemessene Periode: r = {r} (Ledger-Soll: {LEDGER_R})")
        if r != LEDGER_R:
            print(f"FEHLER: Ordnung {r} ≠ Ledger {LEDGER_R}")
            return 1

        print(f"\n[Schritt 5] Klassische Post-Processing Auswertung:")
        result = self.factor_from_period(a, N, r)
        if result is None:
            print("  [-] Trivialer Fall oder ungerade Periode.")
            return 1

        half_pow = result["half_pow"]
        if half_pow != LEDGER_MODULAR_POWER:
            print(
                f"FEHLER: a^(r/2) ≡ {half_pow}, erwartet Ledger-Wert "
                f"{LEDGER_MODULAR_POWER} (nicht {REJECTED_POSTPROCESS})."
            )
            return 1

        print(f"  [+] a^(r/2) ≡ {half_pow} (mod {N})  ← echte nichttriviale 1-Wurzel")
        self.print_crt_audit(half_pow, N, LEDGER_P, LEDGER_Q)

        p, q = result["p"], result["q"]
        print("\n" + "=" * 65)
        print("  ERFOLG: Faktoren erfolgreich berechnet!")
        print(f"  Gegebene Zahl:    N = {N}")
        print(f"  Basis:            a = {a}")
        print(f"  Gemessene Ord:    r = {r}")
        print(f"  Halbperiode:      r/2 = {r // 2}")
        print(f"  Modulare Potenz:  {half_pow}")
        print(f"  Berechnete Kerne: p = {p}, q = {q}")
        print(f"  Verifikation:     {p} * {q} = {p * q} (Korrekt: {result['ok']})")
        print("  FACTORIZATION:    EXACT")
        print("  HARDWARE CLAIM:   NOT MADE")
        print("=" * 65)
        return 0 if result["ok"] else 1

    def run_dialog(self):
        self.print_header()

        raw_input = input(
            "\n[Schritt 1] Gib eine ungerade zusammengesetzte Zahl N ein "
            "(z. B. 15, 21, 35, 77, 1701) oder 'demo': "
        ).strip()
        if raw_input.lower() == "demo":
            raise SystemExit(self.run_ledger_demo(interactive=True))

        if not raw_input.isdigit():
            print("Fehler: Bitte eine gültige Ganzzahl eingeben.")
            return
        N = int(raw_input)

        if N % 2 == 0 or N < 3:
            print(f"Abbruch: N = {N} ist gerade oder zu klein. Klassischer Faktor: 2.")
            return

        print(f"\n-> Ziel: Finde Faktoren von N = {N} im 1000-Qubit-Stabilizerraum.")

        a_raw = input(
            f"\n[Schritt 2] Basis a ∈ [2, N-1] (Enter = zufällig"
            f"{', Ledger a=1502' if N == LEDGER_N else ''}): "
        ).strip()
        if a_raw.isdigit():
            a = int(a_raw)
        elif N == LEDGER_N and a_raw == "":
            a = LEDGER_A
            print(f"  -> Ledger-Basis a = {a}")
        else:
            a = random.randint(2, N - 1)
            print(f"  -> Zufallsbasis a = {a}")

        if not (2 <= a < N):
            print(f"Fehler: a = {a} liegt nicht in [2, {N - 1}].")
            return

        gcd_val = math.gcd(a, N)
        if gcd_val > 1:
            print(f"\n[Glücksfall!] ggT({a}, {N}) = {gcd_val} > 1.")
            print(
                f"-> Echter Faktor klassisch gefunden ohne Quantenregister: "
                f"p = {gcd_val}, q = {N // gcd_val}"
            )
            return

        control_bits, target_bits, free_bits = self.register_layout(N)
        active = control_bits + target_bits
        print(f"\n[Schritt 3] Initialisiere 1000-Qubit Quanten-Register...")
        print(f"  - Steuerregister (Superposition/QFT): {control_bits} Qubits")
        print(f"  - Zielregister (Modul-Exp):           {target_bits} Qubits")
        print(f"  - Aktive Tableau-Qubits:               {active} / 1000")
        print(f"  - Verbleibendes Freiregister:          {free_bits} Qubits")

        self.qc.reset_registers()

        print("\n-> Wende Walsh-Hadamard-Transformation auf Steuerregister an...")
        for q in range(control_bits):
            self.qc.apply_h(q)
        print(
            f"   [OK] {control_bits} Qubits befinden sich in maximaler Superposition "
            f"(Tableau SIMULATED)."
        )

        print("\n" + "-" * 65)
        print("  [EPISTEMISCHE SCHRANKE / EABC CLAIM-WALL]")
        print("  Das 2000 × 2001 Gottesman-Knill Tableau verarbeitet Clifford-Gates.")
        print("  Die modulare Exponentiation f(x) = a^x mod N und die universelle QFT")
        print("  benötigen Nicht-Clifford T-Gates (π/4 Phase).")
        print("  Oracle-Schritt: CLASSICAL / CLAIM-WALL | Hardware: NOT MADE")
        print("-" * 65)

        proceed = input(
            "\nSoll die Periodensuche r = ord_N(a) im Simulator evaluiert werden? (j/n): "
        ).strip().lower()
        if proceed != "j":
            print("Dialog durch Anwender abgebrochen.")
            return

        print(f"\n[Schritt 4] Führe Phasenschätzung & mod-N Exponentiation aus...")
        print("  (klassischer Oracle-Ersatz für Non-Clifford-Schritt)")
        r = self.find_period_classical_oracle(a, N)
        print(f"-> Gemessene Periode der Funktion f(x) = {a}^x mod {N}: r = {r}")

        print(f"\n[Schritt 5] Klassische Post-Processing Auswertung:")
        try:
            result = self.factor_from_period(a, N, r)
        except RuntimeError as exc:
            print(f"  [!] {exc}")
            return

        if result is None:
            if r % 2 != 0:
                print(f"  [-] Periode r = {r} ist ungerade. Versuche andere Basis a.")
            else:
                half_pow = pow(a, r // 2, N)
                print(
                    f"  [-] a^(r/2) ≡ {half_pow} mod N (trivial: ±1). Versuche andere Basis."
                )
            return

        half_pow = result["half_pow"]
        p, q = result["p"], result["q"]
        print(f"  [+] a^(r/2) ≡ {half_pow} (mod {N})")
        if half_pow == REJECTED_POSTPROCESS:
            print("  [!] Warnung: verworfener Fehlwert 1216 — sollte nicht auftreten.")
        if N == LEDGER_N and a == LEDGER_A:
            self.print_crt_audit(half_pow, N, LEDGER_P, LEDGER_Q)

        print("\n" + "=" * 65)
        print("  ERFOLG: Faktoren erfolgreich berechnet!")
        print(f"  Gegebene Zahl:    N = {N}")
        print(f"  Basis:            a = {a}")
        print(f"  Gemessene Ord:    r = {r}")
        print(f"  Halbperiode:      r/2 = {r // 2}")
        print(f"  Modulare Potenz:  {half_pow}")
        print(f"  Berechnete Kerne: p = {p}, q = {q}")
        print(f"  Verifikation:     {p} * {q} = {p * q} (Korrekt: {result['ok']})")
        print("  HARDWARE CLAIM:   NOT MADE")
        print("=" * 65)

    @staticmethod
    def find_period_classical_oracle(a: int, N: int) -> int:
        """Schnelles klassisches Orakel (Lagrange / Faktorisierung) statt O(r).

        Claim-Wall: ersetzt den Nicht-Clifford-QFT-Schritt; kein Hardware-Claim.
        Bei fehlendem SymPy: begrenzter O(r)-Fallback (warnt bei großen N).
        """
        if math.gcd(a, N) != 1:
            raise ValueError(f"gcd({a}, {N}) ≠ 1 — Ordnung in (ℤ/Nℤ)× undefiniert")
        if _sympy_n_order is not None:
            try:
                return int(_sympy_n_order(a, N))
            except Exception as exc:  # noqa: BLE001 — Fallback behält Demo lauffähig
                warnings.warn(
                    f"SymPy-Orakel fehlgeschlagen ({exc!r}); nutze O(r)-Fallback.",
                    RuntimeWarning,
                    stacklevel=2,
                )
        return DialogShorConsole._ord_fallback(a, N)

    @staticmethod
    def _ord_fallback(a: int, N: int, max_steps: int | None = None) -> int:
        """Lineare Suche — nur Notfall. Für ≳30-Bit-N praktisch unbrauchbar."""
        bitlen = N.bit_length()
        if max_steps is None:
            # 34-Bit-N: r kann ~10^9 sein → harte Kappe verhindert „hängendes“ Terminal
            max_steps = 2_000_000 if bitlen >= 28 else N
        if bitlen >= 28:
            warnings.warn(
                f"O(r)-Fallback für N≈2^{bitlen} (max_steps={max_steps}). "
                "Installiere SymPy für schnelles n_order-Orakel.",
                RuntimeWarning,
                stacklevel=3,
            )
        val = a % N
        r = 1
        while val != 1:
            val = (val * a) % N
            r += 1
            if r > max_steps:
                raise RuntimeError(
                    f"Periodensuche abgebrochen nach {max_steps} Schritten "
                    f"(a={a}, N={N}, bitlen={bitlen}). "
                    "Ursache: klassischer Flaschenhals O(r) ohne SymPy-Orakel — "
                    "genau der Nicht-Clifford-Schritt, den Shor quantenbeschleunigt."
                )
        return r


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="EABC Hybrid-Shor 1000-Qubit Simulator")
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Nicht-interaktiver Release-Lauf N=1701, a=1502 → 1457 → 7×243",
    )
    parser.add_argument(
        "--interactive-demo",
        action="store_true",
        help="Ledger-Demo mit Bestätigungsdialog",
    )
    parser.add_argument(
        "--big",
        action="store_true",
        help="34-Bit-Smoke: N=12345678907 mit schnellem SymPy-Orakel",
    )
    parser.add_argument(
        "--bit47",
        action="store_true",
        help="47-Bit-Smoke: N=123456789012345, beide Shor-Teiler zusammengesetzt",
    )
    parser.add_argument(
        "--bit134",
        action="store_true",
        help="134-Bit-Smoke: korrigierte Primstruktur inkl. 43²",
    )
    parser.add_argument(
        "--compare",
        action="store_true",
        help="Vergleichstest Hybrid-Shor vs. Mainstream sympy.factorint",
    )
    parser.add_argument(
        "--N",
        type=int,
        action="append",
        dest="compare_ns",
        help="Mit --compare: zu testende N (mehrfach möglich; Default: 1701 und BIG-SMOKE)",
    )
    args = parser.parse_args(argv)

    console = DialogShorConsole()
    if args.demo:
        return console.run_ledger_demo(interactive=False)
    if args.interactive_demo:
        return console.run_ledger_demo(interactive=True)
    if args.big:
        return console.run_big_smoke()
    if args.bit47:
        return console.run_bit47_smoke()
    if args.bit134:
        return console.run_bit134_smoke()
    if args.compare:
        return console.run_compare(args.compare_ns)
    console.run_dialog()
    return 0


if __name__ == "__main__":
    sys.exit(main())
1701
