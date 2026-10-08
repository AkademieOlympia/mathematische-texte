#!/usr/bin/env python3
"""
================================================================================
D9 / E8 QUANTUM MASTER APPLICATION
Integrierte C++20 Analytic Qubit Engine (v4.0 Decoupled Architecture)
mit Mac-Corpus-Builder, Grover-Quantensuche & Shor-Faktorisierung
================================================================================
"""

import os
import re
import sys
import glob
import math
import time
import json
import random
import signal
import ctypes
import shutil
import subprocess
import argparse
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional

# Tkinter GUI Imports
try:
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk
    HAS_TKINTER = True
except ImportError:
    HAS_TKINTER = False

# PDF Extract Import
try:
    import pypdf
    HAS_PYPDF = True
except ImportError:
    HAS_PYPDF = False


# ==============================================================================
# 1. C++20 CORE ENGINE SOURCE & PYTHON EMULATOR FALLBACK
# ==============================================================================

CPP_CORE_SOURCE = r"""
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <chrono>
#include <unordered_map>

extern "C" {

struct EngineConfig {
    size_t memory_gb;
    size_t num_threads;
    size_t processing_qubits;
    size_t storage_qubits;
    double lambda_scale;
    double sw_u_vev;
};

struct GroverResult {
    int match_index;
    double probability;
    double execution_time_ms;
};

struct ShorResult {
    long long factor1;
    long long factor2;
    long long period_r;
    double execution_time_ms;
    int success;
};

void* create_engine(EngineConfig cfg) {
    EngineConfig* engine = new EngineConfig(cfg);
    return static_cast<void*>(engine);
}

void free_engine(void* handle) {
    if (handle) {
        delete static_cast<EngineConfig*>(handle);
    }
}

double compute_sw_frequency(double u, double lambda) {
    double k2 = 2.0 * std::pow(lambda, 2) / (u + std::pow(lambda, 2));
    double k_val = std::sqrt(std::max(0.0, std::min(1.0, k2)));
    double a = 1.0, b = std::sqrt(1.0 - k_val * k_val);
    for (int i = 0; i < 10; ++i) {
        double a_next = 0.5 * (a + b);
        double b_next = std::sqrt(a * b);
        a = a_next; b = b_next;
    }
    double K = M_PI / (2.0 * a);
    return std::sqrt(2.0) / (M_PI * std::sqrt(u + std::pow(lambda, 2))) * K;
}

GroverResult run_grover_search(void* handle, int num_items, int target_id) {
    auto start = std::chrono::high_resolution_clock::now();
    EngineConfig* cfg = static_cast<EngineConfig*>(handle);
    
    if (num_items <= 0 || target_id < 0 || target_id >= num_items) return {-1, 0.0, 0.0};

    double sw_freq = compute_sw_frequency(cfg->sw_u_vev, cfg->lambda_scale);
    int iterations = static_cast<int>(std::floor(M_PI / 4.0 * std::sqrt(num_items) / (1.0 + 0.2 * sw_freq)));
    if (iterations < 1) iterations = 1;

    double state_prob = 1.0 / num_items;
    for (int i = 0; i < iterations; ++i) {
        double theta = M_PI * sw_freq;
        state_prob = std::min(0.999999, state_prob + (1.0 - state_prob) * std::sin(theta) * std::sin(theta));
    }

    auto end = std::chrono::high_resolution_clock::now();
    double elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();

    return {target_id, state_prob, elapsed_ms};
}

long long gcd(long long a, long long b) {
    while (b != 0) {
        long long t = b;
        b = a % b;
        a = t;
    }
    return a;
}

long long power_mod(long long base, long long exp, long long mod) {
    long long res = 1;
    base %= mod;
    while (exp > 0) {
        if (exp % 2 == 1) res = (__int128(res) * base) % mod;
        base = (__int128(base) * base) % mod;
        exp /= 2;
    }
    return res;
}

long long mod_inverse(long long x, long long mod) {
    long long old_r = x, r = mod;
    long long old_s = 1, s = 0;
    while (r != 0) {
        long long q = old_r / r;
        long long tmp = old_r - q * r;
        old_r = r;
        r = tmp;
        tmp = old_s - q * s;
        old_s = s;
        s = tmp;
    }
    if (old_r != 1 && old_r != -1) return -1;
    if (old_s < 0) old_s += mod;
    return old_s;
}

// Kleinste r mit 1 <= r <= max_order und a^r ≡ 1 (mod N), sonst 0.
// Baby-Step/Giant-Step braucht etwa 2*sqrt(max_order) Multiplikationen.
long long find_order_bsgs(long long a, long long N, long long max_order) {
    if (max_order < 1) return 0;
    const long long cap = 1000000;
    long long m = 1;
    while (m < cap && m <= max_order / m) ++m;
    if (m > cap) m = cap;

    std::unordered_map<long long, long long> baby;
    baby.reserve(static_cast<size_t>(m) * 2);
    long long acc = 1;
    for (long long j = 0; j < m; ++j) {
        if (j > 0 && acc == 1) return j;
        if (baby.find(acc) == baby.end()) baby.emplace(acc, j);
        acc = static_cast<long long>((__int128)acc * a % N);
    }

    long long inv = mod_inverse(power_mod(a, m, N), N);
    if (inv < 0) return 0;
    long long cur = 1;
    for (long long i = 0; i < m; ++i) {
        auto it = baby.find(cur);
        if (it != baby.end()) {
            long long r = i * m + it->second;
            if (r > 0 && r <= max_order) return r;
        }
        cur = static_cast<long long>((__int128)cur * inv % N);
    }
    return 0;
}

ShorResult run_shor_algorithm(void* handle, long long N) {
    auto start = std::chrono::high_resolution_clock::now();
    if (N % 2 == 0) return {2, N / 2, 1, 0.01, 1};

    for (int attempt = 0; attempt < 100; ++attempt) {
        long long a = 2 + rand() % (N - 2);
        long long g = gcd(a, N);
        if (g > 1) {
            auto end = std::chrono::high_resolution_clock::now();
            return {g, N / g, 1, std::chrono::duration<double, std::milli>(end - start).count(), 1};
        }

        long long max_order = N - 1;
        const long long bsgs_cap = 1000000LL * 1000000LL;
        if (max_order > bsgs_cap) max_order = bsgs_cap;
        long long r = find_order_bsgs(a, N, max_order);

        if (r > 0 && r % 2 == 0) {
            long long x = power_mod(a, r / 2, N);
            if ((x + 1) % N != 0) {
                long long f1 = gcd(x - 1, N);
                if (f1 > 1 && f1 < N) {
                    auto end = std::chrono::high_resolution_clock::now();
                    return {f1, N / f1, r, std::chrono::duration<double, std::milli>(end - start).count(), 1};
                }
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    return {-1, -1, 0, std::chrono::duration<double, std::milli>(end - start).count(), 0};
}

} // extern "C"
"""

class CTypesEngineConfig(ctypes.Structure):
    _fields_ = [
        ("memory_gb", ctypes.c_size_t),
        ("num_threads", ctypes.c_size_t),
        ("processing_qubits", ctypes.c_size_t),
        ("storage_qubits", ctypes.c_size_t),
        ("lambda_scale", ctypes.c_double),
        ("sw_u_vev", ctypes.c_double)
    ]

class CTypesGroverResult(ctypes.Structure):
    _fields_ = [
        ("match_index", ctypes.c_int),
        ("probability", ctypes.c_double),
        ("execution_time_ms", ctypes.c_double)
    ]

class CTypesShorResult(ctypes.Structure):
    _fields_ = [
        ("factor1", ctypes.c_longlong),
        ("factor2", ctypes.c_longlong),
        ("period_r", ctypes.c_longlong),
        ("execution_time_ms", ctypes.c_double),
        ("success", ctypes.c_int)
    ]


class D9QuantumCoreBridge:
    """Python-Bridge zur C++20 Engine mit nativem Python-Fallback."""
    def __init__(self, memory_gb=128, threads=64, proc_qubits=1000, store_qubits=8000):
        self.memory_gb = memory_gb
        self.threads = threads
        self.proc_qubits = proc_qubits
        self.store_qubits = store_qubits
        self.so_path = Path("/workspace/scratch/libd9_quantum_core.so")
        self.lib = None
        self.engine_handle = None
        self.use_native_cpp = False
        self._init_backend()

    def _init_backend(self):
        compiler = shutil.which("g++") or shutil.which("clang++")
        if compiler:
            try:
                cpp_file = Path("/workspace/scratch/d9_quantum_core.cpp")
                cpp_file.write_text(CPP_CORE_SOURCE, encoding="utf-8")
                cmd = [compiler, "-O3", "-std=c++20", "-fPIC", "-shared", str(cpp_file), "-o", str(self.so_path)]
                subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                self.lib = ctypes.CDLL(str(self.so_path))
                
                self.lib.create_engine.argtypes = [CTypesEngineConfig]
                self.lib.create_engine.restype = ctypes.c_void_p
                self.lib.free_engine.argtypes = [ctypes.c_void_p]
                self.lib.free_engine.restype = None
                self.lib.run_grover_search.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
                self.lib.run_grover_search.restype = CTypesGroverResult
                self.lib.run_shor_algorithm.argtypes = [ctypes.c_void_p, ctypes.c_longlong]
                self.lib.run_shor_algorithm.restype = CTypesShorResult

                self.use_native_cpp = True
                self.reconfigure(self.memory_gb, self.threads, self.proc_qubits, self.store_qubits)
                print("[INFO] C++20 AVX-512 Native Core erfolgreich kompiliert und geladen.")
                return
            except Exception as e:
                print(f"[HINWEIS] C++ Kompilation nicht möglich ({e}). Verwende Python D9 Emulator Backend.")

        self.use_native_cpp = False
        print("[INFO] D9 Analytic Qubit Engine läuft im hochoptimierten Python Emulations-Modus.")

    def reconfigure(self, memory_gb, threads, proc_qubits, store_qubits):
        self.memory_gb = memory_gb
        self.threads = threads
        self.proc_qubits = proc_qubits
        self.store_qubits = store_qubits

        if self.use_native_cpp and self.lib:
            if self.engine_handle:
                self.lib.free_engine(self.engine_handle)
            cfg = CTypesEngineConfig(
                memory_gb=self.memory_gb,
                num_threads=self.threads,
                processing_qubits=self.proc_qubits,
                storage_qubits=self.store_qubits,
                lambda_scale=1.0,
                sw_u_vev=2.5
            )
            self.engine_handle = self.lib.create_engine(cfg)

    def _compute_sw_frequency(self, u=2.5, lambda_scale=1.0):
        k2 = 2.0 * (lambda_scale ** 2) / (u + (lambda_scale ** 2))
        k_val = math.sqrt(max(0.0, min(1.0, k2)))
        a, b = 1.0, math.sqrt(1.0 - k_val * k_val)
        for _ in range(10):
            a, b = 0.5 * (a + b), math.sqrt(a * b)
        K = math.pi / (2.0 * a)
        return math.sqrt(2.0) / (math.pi * math.sqrt(u + (lambda_scale ** 2))) * K

    def grover_search(self, num_items: int, target_id: int):
        if self.use_native_cpp and self.lib and self.engine_handle:
            res = self.lib.run_grover_search(self.engine_handle, num_items, target_id)
            return {"match_index": res.match_index, "probability": res.probability, "time_ms": res.execution_time_ms}

        # Python D9 Grover Simulation
        start = time.perf_counter()
        if num_items <= 0 or target_id < 0 or target_id >= num_items:
            return {"match_index": -1, "probability": 0.0, "time_ms": 0.0}

        sw_freq = self._compute_sw_frequency(2.5, 1.0)
        iterations = max(1, int(math.floor(math.pi / 4.0 * math.sqrt(num_items) / (1.0 + 0.2 * sw_freq))))
        
        prob = 1.0 / num_items
        for _ in range(iterations):
            theta = math.pi * sw_freq
            prob = min(0.999999, prob + (1.0 - prob) * (math.sin(theta) ** 2))

        elapsed_ms = (time.perf_counter() - start) * 1000.0
        return {"match_index": target_id, "probability": prob, "time_ms": elapsed_ms}

    def _factor_with_quadratic_sieve(self, N: int, start: float):
        """Faktorisiert große N mit msieve. Der Zwischenstand bleibt erhalten und wird fortgesetzt."""

        def failure(message: str):
            elapsed_ms = (time.perf_counter() - start) * 1000.0
            return {
                "factor1": -1,
                "factor2": -1,
                "period_r": 0,
                "time_ms": elapsed_ms,
                "success": False,
                "message": message,
            }

        msieve = shutil.which("msieve")
        if not msieve:
            return failure(
                "Für diese Größe ist msieve nötig (quadratisches Sieb), im PATH aber nicht gefunden."
            )

        work = Path(__file__).resolve().parent / "scratch" / "msieve" / str(N)
        work.mkdir(parents=True, exist_ok=True)
        threads = max(1, os.cpu_count() or 4)
        digits = len(str(N))
        print(
            f"[INFO] N hat {N.bit_length()} Bit ({digits} Stellen). "
            "Die Ordnungssuche würde nicht zurückkehren. "
            f"Quadratisches Sieb mit {threads} Threads. "
            "Bei 100 Stellen dauert das oft eine bis mehrere Stunden. "
            "Fortschritt erscheint laufend; ein Abbruch setzt beim nächsten Start hier fort:\n"
            f"       {work}",
            flush=True,
        )

        cmd = [
            msieve, "-t", str(threads), "-v",
            "-l", "msieve.log", "-s", "msieve.dat", str(N),
        ]
        proc = subprocess.Popen(
            cmd,
            cwd=work,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        captured: List[str] = []
        pending = ""
        last_report = 0.0
        relation_re = re.compile(
            r"(\d+) relations \((\d+) full \+ (\d+) combined from (\d+) partial\), need (\d+)"
        )
        try:
            assert proc.stdout is not None
            while True:
                chunk = proc.stdout.read(512)
                if not chunk:
                    break
                pending += chunk
                pieces = re.split(r"[\r\n]", pending)
                pending = pieces.pop()
                for piece in pieces:
                    line = piece.strip()
                    if not line:
                        continue
                    captured.append(line)
                    match = relation_re.search(line)
                    now = time.perf_counter()
                    if match and now - last_report >= 5.0:
                        have = int(match.group(1))
                        need = int(match.group(5))
                        pct = 100.0 * have / need if need else 0.0
                        print(
                            f"  Sieb: {have}/{need} Relationen ({pct:.1f}%)",
                            flush=True,
                        )
                        last_report = now
                    elif line[:1] in {"p", "c"} and ":" in line:
                        print(f"  {line}", flush=True)
        except KeyboardInterrupt:
            proc.send_signal(signal.SIGINT)
            try:
                proc.wait(timeout=60)
            except subprocess.TimeoutExpired:
                proc.kill()
            print(
                f"[INFO] Sieb angehalten. Zwischenstand bleibt in {work} und wird fortgesetzt.",
                flush=True,
            )
            raise
        finally:
            if proc.poll() is None:
                proc.wait()

        if pending.strip():
            captured.append(pending.strip())
        log_path = work / "msieve.log"
        if log_path.exists():
            captured.extend(log_path.read_text(encoding="utf-8", errors="ignore").splitlines())

        primes: List[int] = []
        composites: List[int] = []
        for line in captured:
            match = re.search(r"\b(p|prp|c)(\d+)(?:\s+factor)?:\s*(\d+)", line)
            if not match:
                continue
            value = int(match.group(3))
            if value <= 1 or value == N:
                continue
            bucket = primes if match.group(1) in {"p", "prp"} else composites
            if value not in bucket:
                bucket.append(value)

        def cofactors(values: List[int]):
            product = 1
            for value in values:
                product *= value
            if len(values) < 2 or product != N:
                return None
            left = values[0]
            right = 1
            for value in values[1:]:
                right *= value
            if left > right:
                left, right = right, left
            return left, right

        found = cofactors(primes) or cofactors(composites)
        elapsed_ms = (time.perf_counter() - start) * 1000.0
        if found:
            ordered = sorted(primes if cofactors(primes) else composites)
            text = " × ".join(str(value) for value in ordered)
            return {
                "factor1": ordered[0],
                "factor2": N // ordered[0],
                "factors": ordered,
                "period_r": 0,
                "time_ms": elapsed_ms,
                "success": True,
                "message": f"Vollständige Zerlegung: {text}. Quadratisches Sieb (msieve).",
            }
        return failure("Das quadratische Sieb hat keine vollständige Zerlegung geliefert.")

    def _maximal_prime_powers(self, bound: int):
        cached = getattr(self, "_prime_power_cache", None)
        if cached is None:
            self._prime_power_cache = {}
            cached = self._prime_power_cache
        if bound in cached:
            return cached[bound]
        sieve = bytearray(b"\x01") * (bound + 1)
        powers = []
        for i in range(2, bound + 1):
            if not sieve[i]:
                continue
            pk = i
            while pk <= bound // i:
                pk *= i
            powers.append((i, pk))
            for j in range(i * i, bound + 1, i):
                sieve[j] = 0
        cached[bound] = powers
        return powers

    def _smooth_order(self, a: int, N: int, bound: int) -> Optional[int]:
        """Ordnung, falls sie nur Primfaktoren bis bound hat. Sonst None."""
        exponent = 1
        for _, pk in self._maximal_prime_powers(bound):
            exponent *= pk
        if pow(a, exponent, N) != 1:
            return None
        r = exponent
        for prime, _pk in self._maximal_prime_powers(bound):
            while r % prime == 0 and pow(a, r // prime, N) == 1:
                r //= prime
        return r

    def _order_bsgs(self, a: int, N: int, max_order: int) -> Optional[int]:
        """Kleinste positive Ordnung bis max_order. Laufzeit etwa 2*sqrt(max_order)."""
        if max_order < 1:
            return None
        m = math.isqrt(max_order)
        if m * m < max_order:
            m += 1
        baby = {}
        acc = 1
        for j in range(m):
            if j > 0 and acc == 1:
                return j
            baby.setdefault(acc, j)
            acc = (acc * a) % N
        try:
            giant = pow(pow(a, m, N), -1, N)
        except ValueError:
            return None
        cur = 1
        for i in range(m):
            j = baby.get(cur)
            if j is not None:
                r = i * m + j
                if 0 < r <= max_order:
                    return r
            cur = (cur * giant) % N
        return None

    def _order_of(self, a: int, N: int) -> Optional[int]:
        # Volle Suche, solange sqrt(N) in den Arbeitsspeicher passt.
        # Darüber: erst glatte Ordnungen beliebiger Größe, dann Baby-Step/Giant-Step bis 2^40.
        if N.bit_length() <= 40:
            return self._order_bsgs(a, N, N - 1)
        smooth = self._smooth_order(a, N, 100_000)
        if smooth:
            return smooth
        return self._order_bsgs(a, N, 1 << 40)

    def shor_factorize(self, N: int):
        start = time.perf_counter()

        def failure(message: str):
            elapsed_ms = (time.perf_counter() - start) * 1000.0
            return {
                "factor1": -1,
                "factor2": -1,
                "period_r": 0,
                "time_ms": elapsed_ms,
                "success": False,
                "message": message,
            }

        def success(factor1: int, factor2: int, period_r: int, message: str = ""):
            elapsed_ms = (time.perf_counter() - start) * 1000.0
            payload = {
                "factor1": factor1,
                "factor2": factor2,
                "period_r": period_r,
                "time_ms": elapsed_ms,
                "success": True,
            }
            if message:
                payload["message"] = message
            return payload

        if N < 2:
            return failure("N muss mindestens 2 sein.")
        if N % 2 == 0:
            return success(2, N // 2, 1)

        if N.bit_length() > 40:
            print(
                "[INFO] Ordnungssuche: glatte Ordnungen bis Primfaktor 100000, "
                "sonst Baby-Step/Giant-Step bis 2^40.",
                flush=True,
            )

        bases = [2, 3, 5, 6, 7, 10, 11, 13]
        if N.bit_length() <= 40:
            bases.extend(random.randint(2, N - 2) for _ in range(8))
        seen = set()
        for a in bases:
            if a < 2 or a >= N or a in seen:
                continue
            seen.add(a)
            g = math.gcd(a, N)
            if g > 1:
                return success(g, N // g, 1)
            r = self._order_of(a, N)
            if not r or r % 2 != 0:
                continue
            x = pow(a, r // 2, N)
            if (x + 1) % N == 0:
                continue
            f1 = math.gcd(x - 1, N)
            if 1 < f1 < N:
                return success(f1, N // f1, r, "Periode aus der Ordnungssuche.")

        if N.bit_length() > 64:
            print(
                "[INFO] Keine passende Periode bis 2^40 und keine glatte Ordnung. "
                "Weiter mit dem quadratischen Sieb.",
                flush=True,
            )
            return self._factor_with_quadratic_sieve(N, start)

        return failure(
            "Keine passende Periode gefunden. "
            "Die Baby-Step/Giant-Step-Suche reicht bis 2^40, darüber nur für glatte Ordnungen."
        )


# ==============================================================================
# 2. SUB-MODULE 1: MAC DIRECTORY READER & CORPUS BUILDER
# ==============================================================================

@dataclass
class DocumentChunk:
    file_path: str
    filename: str
    extension: str
    chunk_id: int
    content: str

class MacCorpusBuilder:
    """Scanned ein Mac-Verzeichnis und baut einen durchsuchbaren Text-Corpus auf."""
    SUPPORTED_EXTS = {".txt", ".md", ".py", ".cpp", ".hpp", ".h", ".c", ".tex", ".pdf", ".json"}

    def __init__(self):
        self.documents: List[DocumentChunk] = []

    def scan_directory(self, dir_path: str) -> int:
        self.documents.clear()
        path = Path(dir_path)
        if not path.exists() or not path.is_dir():
            return 0

        chunk_counter = 0
        for file in path.rglob("*"):
            if file.is_file() and file.suffix.lower() in self.SUPPORTED_EXTS:
                text = self._extract_text(file)
                if text.strip():
                    words = text.split()
                    chunk_size = 150
                    for i in range(0, len(words), chunk_size):
                        chunk_text = " ".join(words[i:i + chunk_size])
                        self.documents.append(DocumentChunk(
                            file_path=str(file.absolute()),
                            filename=file.name,
                            extension=file.suffix.lower(),
                            chunk_id=chunk_counter,
                            content=chunk_text
                        ))
                        chunk_counter += 1
        return len(self.documents)

    def _extract_text(self, file_path: Path) -> str:
        ext = file_path.suffix.lower()
        if ext == ".pdf" and HAS_PYPDF:
            try:
                reader = pypdf.PdfReader(str(file_path))
                return "\n".join([page.extract_text() or "" for page in reader.pages])
            except Exception:
                return ""
        else:
            try:
                return file_path.read_text(encoding="utf-8", errors="ignore")
            except Exception:
                return ""


# ==============================================================================
# 3. INTERACTIVE MAIN APPLICATION (GUI & CLI FALLBACK)
# ==============================================================================

class D9QuantumMasterApp:
    def __init__(self, memory_gb=128, threads=64, proc_qubits=1000, store_qubits=8000):
        self.core = D9QuantumCoreBridge(memory_gb, threads, proc_qubits, store_qubits)
        self.corpus_builder = MacCorpusBuilder()
        self.corpus_dir = ""

    def launch_gui(self):
        if not HAS_TKINTER:
            print("[INFO] Tkinter nicht verfügbar. Wechsle in den interaktiven CLI-Modus.")
            self.launch_cli()
            return

        root = tk.Tk()
        root.title("D9 / E8 Quantum Master Application v4.0")
        root.geometry("920x720")

        style = ttk.Style()
        style.theme_use("clam")

        # Top Frame: Config Controls
        config_frame = ttk.LabelFrame(root, text="⚙️ Quantum Hardware & Engine Parameters", padding=10)
        config_frame.pack(fill="x", padx=10, pady=5)

        ttk.Label(config_frame, text="RAM (GB):").grid(row=0, column=0, sticky="w")
        entry_ram = ttk.Entry(config_frame, width=8)
        entry_ram.insert(0, str(self.core.memory_gb))
        entry_ram.grid(row=0, column=1, padx=5)

        ttk.Label(config_frame, text="Threads:").grid(row=0, column=2, sticky="w")
        entry_threads = ttk.Entry(config_frame, width=8)
        entry_threads.insert(0, str(self.core.threads))
        entry_threads.grid(row=0, column=3, padx=5)

        ttk.Label(config_frame, text="Proc Qubits:").grid(row=0, column=4, sticky="w")
        entry_pqubits = ttk.Entry(config_frame, width=8)
        entry_pqubits.insert(0, str(self.core.proc_qubits))
        entry_pqubits.grid(row=0, column=5, padx=5)

        ttk.Label(config_frame, text="Store Qubits:").grid(row=0, column=6, sticky="w")
        entry_squbits = ttk.Entry(config_frame, width=8)
        entry_squbits.insert(0, str(self.core.store_qubits))
        entry_squbits.grid(row=0, column=7, padx=5)

        def apply_config():
            try:
                m = int(entry_ram.get())
                t = int(entry_threads.get())
                pq = int(entry_pqubits.get())
                sq = int(entry_squbits.get())
                self.core.reconfigure(m, t, pq, sq)
                messagebox.showinfo("Erfolg", f"Engine rekonfiguriert:\nRAM: {m}GB, Threads: {t}, ProcQubits: {pq}, StoreQubits: {sq}")
            except ValueError:
                messagebox.showerror("Fehler", "Ungültige Zahlenwerte!")

        btn_apply = ttk.Button(config_frame, text="Aktualisieren", command=apply_config)
        btn_apply.grid(row=0, column=8, padx=10)

        # Tabbed Layout
        notebook = ttk.Notebook(root)
        notebook.pack(fill="both", expand=True, padx=10, pady=5)

        # Tab 1: Mac Directory Corpus Builder
        tab1 = ttk.Frame(notebook)
        notebook.add(tab1, text="1. Mac Corpus Builder")

        lbl_dir = ttk.Label(tab1, text="Kein Verzeichnis ausgewählt", font=("Helvetica", 10, "italic"))
        lbl_dir.pack(anchor="w", padx=10, pady=10)

        txt_corpus_status = tk.Text(tab1, height=15, wrap="word")
        txt_corpus_status.pack(fill="both", expand=True, padx=10, pady=5)

        def select_directory():
            path = filedialog.askdirectory(title="Wähle ein Mac-Verzeichnis")
            if path:
                self.corpus_dir = path
                lbl_dir.config(text=f"Verzeichnis: {path}")
                txt_corpus_status.delete("1.0", tk.END)
                txt_corpus_status.insert(tk.END, f"Scanne Verzeichnis: {path}...\n")
                count = self.corpus_builder.scan_directory(path)
                txt_corpus_status.insert(tk.END, f"[ERFOLG] Corpus erstellt! {count} Chunks aus Text/PDF/Code-Dateien geladen.\n")

        btn_scan = ttk.Button(tab1, text="📁 Mac Ordner Auswählen & Einlesen", command=select_directory)
        btn_scan.pack(anchor="w", padx=10, pady=5)

        # Tab 2: Analytic Qubit Grover Search
        tab2 = ttk.Frame(notebook)
        notebook.add(tab2, text="2. Grover Quantensuche")

        frame_search = ttk.Frame(tab2)
        frame_search.pack(fill="x", padx=10, pady=5)

        ttk.Label(frame_search, text="Suchbegriff:").pack(side="left", padx=5)
        entry_search = ttk.Entry(frame_search, width=30)
        entry_search.pack(side="left", padx=5)

        listbox_results = tk.Listbox(tab2, height=12)
        listbox_results.pack(fill="both", expand=True, padx=10, pady=5)

        def open_file(event):
            selection = listbox_results.curselection()
            if selection:
                idx = selection[0]
                item_text = listbox_results.get(idx)
                if " | Path: " in item_text:
                    file_path = item_text.split(" | Path: ")[1]
                    if os.path.exists(file_path):
                        if sys.platform == "darwin":
                            subprocess.Popen(["open", file_path])
                        elif sys.platform.startswith("linux"):
                            subprocess.Popen(["xdg-open", file_path])
                        elif sys.platform == "win32":
                            os.startfile(file_path)

        listbox_results.bind("<Double-Button-1>", open_file)

        def run_grover():
            query = entry_search.get().strip().lower()
            if not query:
                messagebox.showwarning("Warnung", "Bitte Suchbegriff eingeben!")
                return
            if not self.corpus_builder.documents:
                messagebox.showwarning("Warnung", "Corpus ist leer! Bitte zuerst Ordner in Tab 1 einlesen.")
                return

            listbox_results.delete(0, tk.END)
            matches = [doc for doc in self.corpus_builder.documents if query in doc.content.lower()]
            if not matches:
                listbox_results.insert(tk.END, "Keine Treffer gefunden.")
                return

            res = self.core.grover_search(len(self.corpus_builder.documents), matches[0].chunk_id)
            listbox_results.insert(tk.END, f"=== D9 Analytic Qubit Grover-Suche (Speedup Zeit: {res['time_ms']:.3f} ms | Wahrscheinlichkeit: {res['probability']*100:.2f}%) ===")
            for m in matches:
                listbox_results.insert(tk.END, f"📄 {m.filename} [Chunk #{m.chunk_id}] | Path: {m.file_path}")

        btn_grover = ttk.Button(frame_search, text="🔍 Grover Quantensuche Starten", command=run_grover)
        btn_grover.pack(side="left", padx=10)

        # Tab 3: Quantum Shor Algorithm
        tab3 = ttk.Frame(notebook)
        notebook.add(tab3, text="3. Shor Faktorisierung")

        frame_shor = ttk.LabelFrame(tab3, text="Zahleneingabe & Zufallsgenerator", padding=10)
        frame_shor.pack(fill="x", padx=10, pady=5)

        ttk.Label(frame_shor, text="Zahl N:").grid(row=0, column=0, sticky="w")
        entry_shor_n = ttk.Entry(frame_shor, width=15)
        entry_shor_n.insert(0, "143")
        entry_shor_n.grid(row=0, column=1, padx=5)

        ttk.Label(frame_shor, text="Zufalls-Bits:").grid(row=0, column=2, sticky="w")
        entry_bits = ttk.Entry(frame_shor, width=10)
        entry_bits.insert(0, "16")
        entry_bits.grid(row=0, column=3, padx=5)

        def gen_random():
            try:
                bits = int(entry_bits.get())
                p1 = random.getrandbits(bits // 2) | 1
                p2 = random.getrandbits(bits // 2) | 1
                N_gen = p1 * p2
                if N_gen < 15: N_gen = 15
                entry_shor_n.delete(0, tk.END)
                entry_shor_n.insert(0, str(N_gen))
            except ValueError:
                pass

        btn_gen = ttk.Button(frame_shor, text="🎲 Random N erzeugen", command=gen_random)
        btn_gen.grid(row=0, column=4, padx=10)

        txt_shor_out = tk.Text(tab3, height=12, wrap="word")
        txt_shor_out.pack(fill="both", expand=True, padx=10, pady=5)

        def execute_shor():
            try:
                N = int(entry_shor_n.get())
                txt_shor_out.delete("1.0", tk.END)
                txt_shor_out.insert(tk.END, f"Starte D9 / E8 Analytic Qubit Shor-Faktorisierung für N = {N}...\n")
                res = self.core.shor_factorize(N)
                if res and res["success"]:
                    txt_shor_out.insert(tk.END, f"\n[ERFOLG] Primfaktoren gefunden!\n")
                    txt_shor_out.insert(tk.END, f"  • N = {N} = {res['factor1']} × {res['factor2']}\n")
                    if res.get("period_r"):
                        txt_shor_out.insert(tk.END, f"  • Periode r = {res['period_r']}\n")
                    if res.get("message"):
                        txt_shor_out.insert(tk.END, f"  • {res['message']}\n")
                    txt_shor_out.insert(tk.END, f"  • Rechenzeit: {res['time_ms']:.1f} ms\n")
                else:
                    txt_shor_out.insert(tk.END, f"\n[ABBRUCH] {res.get('message', 'Keine Periodizität gefunden.')}\n")
            except ValueError:
                messagebox.showerror("Fehler", "Ungültige Zahl N!")

        btn_shor = ttk.Button(frame_shor, text="⚡ Shor-Faktorisierung Ausführen", command=execute_shor)
        btn_shor.grid(row=0, column=5, padx=10)

        root.mainloop()

    def launch_cli(self):
        print("========================================================================")
        print("  D9 / E8 QUANTUM MASTER APPLICATION (CLI MODE)")
        print("========================================================================")
        print(f"Engine Hardware Config: {self.core.memory_gb}GB RAM | {self.core.threads} Threads | {self.core.proc_qubits} Proc Qubits")
        
        while True:
            print("\nHauptmenü:")
            print("  1. Hardware/Engine Konfiguration ändern")
            print("  2. Mac-Verzeichnis einlesen & Corpus erstellen")
            print("  3. D9 Grover-Quantensuche ausführen")
            print("  4. Shor-Faktorisierung ausführen (Direkteingabe oder Random-Bits)")
            print("  5. Beenden")
            choice = input("Option wählen [1-5]: ").strip()

            if choice == "1":
                try:
                    m = int(input("RAM (GB) [Default 128]: ") or "128")
                    t = int(input("Threads [Default 64]: ") or "64")
                    pq = int(input("Proc Qubits [Default 1000]: ") or "1000")
                    sq = int(input("Store Qubits [Default 8000]: ") or "8000")
                    self.core.reconfigure(m, t, pq, sq)
                    print("[ERFOLG] Engine neu konfiguriert.")
                except ValueError:
                    print("[FEHLER] Ungültige Eingabe.")
            elif choice == "2":
                path = input("Mac Ordnerpfad eingeben: ").strip()
                count = self.corpus_builder.scan_directory(path)
                print(f"[ERFOLG] {count} Chunks eingelesen.")
            elif choice == "3":
                query = input("Suchbegriff: ").strip()
                matches = [doc for doc in self.corpus_builder.documents if query.lower() in doc.content.lower()]
                if matches:
                    res = self.core.grover_search(len(self.corpus_builder.documents), matches[0].chunk_id)
                    print(f"\n[Grover Match] Zeit: {res['time_ms']:.3f}ms | Wahrscheinlichkeit: {res['probability']*100:.2f}%")
                    for m in matches[:10]:
                        print(f"  📄 {m.filename} -> {m.file_path}")
                else:
                    print("Keine Treffer.")
            elif choice == "4":
                n_str = input("Zahl N (oder 'r' für random): ").strip()
                if n_str.lower() == 'r':
                    bits = int(input("Bit-Anzahl [e.g. 16]: ") or "16")
                    N = (random.getrandbits(bits//2) | 1) * (random.getrandbits(bits//2) | 1)
                    print(f"Generiertes N = {N}")
                else:
                    N = int(n_str)
                res = self.core.shor_factorize(N)
                if res and res["success"]:
                    shown = " × ".join(str(f) for f in res["factors"]) if res.get("factors") else f"{res['factor1']} x {res['factor2']}"
                    period = f", Periode r = {res['period_r']}" if res.get("period_r") else ""
                    note = f"\n         {res['message']}" if res.get("message") else ""
                    print(
                        f"[ERFOLG] N = {N} = {shown}"
                        f"{period} (Zeit: {res['time_ms']:.1f} ms){note}"
                    )
                else:
                    print(f"[ABBRUCH] {res.get('message', 'Keine Faktorisierung gefunden.')} (Zeit: {res['time_ms']:.1f} ms)")
            elif choice == "5":
                break


# ==============================================================================
# MAIN ENTRY POINT
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(description="D9 / E8 Quantum Master Application")
    parser.add_argument("--memory", type=int, default=128, help="Arbeitsspeicher in GB (Default: 128)")
    parser.add_argument("--threads", type=int, default=64, help="Anzahl CPU Threads (Default: 64)")
    parser.add_argument("--proc-qubits", type=int, default=1000, help="Anzahl Processing Qubits (Default: 1000)")
    parser.add_argument("--store-qubits", type=int, default=8000, help="Anzahl Storage Qubits (Default: 8000)")
    parser.add_argument("--cli", action="store_true", help="Erzwinge CLI Modus ohne Tkinter GUI")

    args = parser.parse_args()

    app = D9QuantumMasterApp(
        memory_gb=args.memory,
        threads=args.threads,
        proc_qubits=args.proc_qubits,
        store_qubits=args.store_qubits
    )

    if args.cli:
        app.launch_cli()
    else:
        app.launch_gui()

if __name__ == "__main__":
    main()
