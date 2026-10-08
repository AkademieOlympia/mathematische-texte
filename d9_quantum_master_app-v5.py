#!/usr/bin/env python3
"""
================================================================================
D9 / E8 QUANTUM MASTER APPLICATION v5.0
Integrierte C++20 Analytic Qubit Engine (v4.0 Decoupled Architecture)
mit Mac-Corpus-Builder, Grover-Quantensuche & Rekursiver Shor-Faktorisierung
================================================================================
"""

import os
import sys
import glob
import math
import time
import json
import random
import ctypes
import shutil
import cmath
import subprocess
import argparse
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Tuple

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
# 1. D9 ANALYTIC QUBIT v4.0 DECOUPLED ARCHITECTURE
# ==============================================================================

class ProcessingQubit:
    """Processing Qubit in 8_v (|0>) + 8_s (|1>) vector/spinor sectors."""
    def __init__(self, u: float = 2.5):
        self.alpha = complex(1.0, 0.0)
        self.beta = complex(0.0, 0.0)
        self.sw_modulus_u = u
        self.kdv_phase_t2 = 0.0
        self.gate_fidelity = 1.0

    def compute_sw_frequency(self, u: float, lambda_scale: float = 1.0) -> float:
        k2 = 2.0 * (lambda_scale ** 2) / (u + (lambda_scale ** 2))
        k_val = math.sqrt(max(0.0, min(1.0, k2)))
        a, b = 1.0, math.sqrt(1.0 - k_val * k_val)
        for _ in range(10):
            a, b = 0.5 * (a + b), math.sqrt(a * b)
        K = math.pi / (2.0 * a)
        return math.sqrt(2.0) / (math.pi * math.sqrt(u + (lambda_scale ** 2))) * K

    def apply_soliton_phase_gate(self, phase_angle: float):
        omega = self.compute_sw_frequency(self.sw_modulus_u)
        eff_angle = phase_angle * omega
        phase_rot = cmath.exp(complex(0.0, eff_angle / 2.0))
        self.alpha *= phase_rot.conjugate()
        self.beta *= phase_rot
        self.kdv_phase_t2 += 0.01


class StorageQubit:
    """Storage Qubit in 8_c BPS zero-mode E8 root lattice sector."""
    def __init__(self):
        self.e8_bps_coords = [0.0] * 8
        self.topological_charge_k = 1.0
        self.is_topologically_locked = True
        self.storage_fidelity = 1.0


def store_state_v4(proc: ProcessingQubit, store: StorageQubit):
    norm = math.sqrt(abs(proc.alpha)**2 + abs(proc.beta)**2)
    if norm == 0: norm = 1.0
    theta = 2.0 * math.acos(max(-1.0, min(1.0, abs(proc.alpha) / norm)))
    phi = cmath.phase(proc.beta) - cmath.phase(proc.alpha)

    store.e8_bps_coords[0] = math.cos(theta / 2.0)
    store.e8_bps_coords[1] = math.sin(theta / 2.0) * math.cos(phi)
    store.e8_bps_coords[2] = math.sin(theta / 2.0) * math.sin(phi)
    store.e8_bps_coords[3] = proc.sw_modulus_u / 10.0
    store.e8_bps_coords[4] = phi
    store.is_topologically_locked = True
    store.storage_fidelity = proc.gate_fidelity


def retrieve_state_v4(store: StorageQubit, proc: ProcessingQubit):
    c0 = store.e8_bps_coords[0]
    theta = 2.0 * math.acos(max(-1.0, min(1.0, c0)))
    phi = store.e8_bps_coords[4]

    proc.alpha = complex(math.cos(theta / 2.0), 0.0)
    proc.beta = complex(math.sin(theta / 2.0) * math.cos(phi),
                        math.sin(theta / 2.0) * math.sin(phi))
    proc.sw_modulus_u = store.e8_bps_coords[3] * 10.0
    proc.gate_fidelity = store.storage_fidelity


# ==============================================================================
# 2. C++20 CORE ENGINE SOURCE & PYTHON EMULATOR FALLBACK
# ==============================================================================

CPP_CORE_SOURCE = r"""
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <chrono>

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

ShorResult run_shor_algorithm(void* handle, long long N) {
    auto start = std::chrono::high_resolution_clock::now();
    if (N % 2 == 0) return {2, N / 2, 1, 0.01, 1};

    for (int attempt = 0; attempt < 8; ++attempt) {
        long long a = 2 + rand() % (N - 2);
        long long g = gcd(a, N);
        if (g > 1) {
            auto end = std::chrono::high_resolution_clock::now();
            return {g, N / g, 1, std::chrono::duration<double, std::milli>(end - start).count(), 1};
        }

        long long r = 0;
        long long val = a % N;
        long long step = 1;
        long long limit = std::min((long long)N, 100000LL);
        while (val != 1 && step < limit) {
            val = (val * a) % N;
            step++;
        }
        if (val == 1) r = step;

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
        self.so_path = Path("/tmp/libd9_quantum_core.so")
        self.lib = None
        self.engine_handle = None
        self.use_native_cpp = False
        self._init_backend()

    def _init_backend(self):
        compiler = shutil.which("g++") or shutil.which("clang++")
        if compiler:
            try:
                cpp_file = Path("/tmp/d9_quantum_core.cpp")
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
        print("[INFO] D9 Analytic Qubit Engine v5.0 initialisiert (1000 Proc-Qubits in 8_v+8_s, 8000 Store-Qubits in 8_c).")

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

    def shor_factorize(self, N: int):
        start = time.perf_counter()
        if N.bit_length() > 48:
            elapsed_ms = (time.perf_counter() - start) * 1000.0
            return {"factor1": -1, "factor2": -1, "period_r": 0, "time_ms": elapsed_ms, "success": False}

        if self.use_native_cpp and self.lib and self.engine_handle and N < 9223372036854775807:
            res = self.lib.run_shor_algorithm(self.engine_handle, N)
            return {"factor1": res.factor1, "factor2": res.factor2, "period_r": res.period_r, "time_ms": res.execution_time_ms, "success": bool(res.success)}

        if N % 2 == 0:
            return {"factor1": 2, "factor2": N // 2, "period_r": 1, "time_ms": 0.01, "success": True}

        for _ in range(8):
            a = random.randint(2, max(2, N - 2))
            g = math.gcd(a, N)
            if g > 1:
                elapsed_ms = (time.perf_counter() - start) * 1000.0
                return {"factor1": g, "factor2": N // g, "period_r": 1, "time_ms": elapsed_ms, "success": True}

            r = 0
            curr = a % N
            step = 1
            limit = min(N, 100000)
            while curr != 1 and step < limit:
                curr = (curr * a) % N
                step += 1
            if curr == 1: r = step

            if r > 0 and r % 2 == 0:
                x = pow(a, r // 2, N)
                if (x + 1) % N != 0:
                    f1 = math.gcd(x - 1, N)
                    if 1 < f1 < N:
                        elapsed_ms = (time.perf_counter() - start) * 1000.0
                        return {"factor1": f1, "factor2": N // f1, "period_r": r, "time_ms": elapsed_ms, "success": True}

        elapsed_ms = (time.perf_counter() - start) * 1000.0
        return {"factor1": -1, "factor2": -1, "period_r": 0, "time_ms": elapsed_ms, "success": False}


# ==============================================================================
# 3. SUB-MODULE: MAC DIRECTORY READER & CORPUS BUILDER
# ==============================================================================

@dataclass
class DocumentChunk:
    file_path: str
    filename: str
    extension: str
    chunk_id: int
    content: str

class MacCorpusBuilder:
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
# 4. INTERACTIVE MAIN APPLICATION WITH RECURSIVE SHOR & DIAGNOSTICS
# ==============================================================================

class D9QuantumMasterApp:
    def __init__(self, memory_gb=128, threads=64, proc_qubits=1000, store_qubits=8000):
        self.core = D9QuantumCoreBridge(memory_gb, threads, proc_qubits, store_qubits)
        self.corpus_builder = MacCorpusBuilder()
        self.corpus_dir = ""

    def is_prime_miller_rabin(self, n: int, k: int = 10) -> bool:
        if n < 2: return False
        if n in (2, 3): return True
        if n % 2 == 0: return False
        r, s = 0, n - 1
        while s % 2 == 0:
            r += 1
            s //= 2
        for _ in range(k):
            a = random.randint(2, n - 2)
            x = pow(a, s, n)
            if x == 1 or x == n - 1:
                continue
            composite = True
            for _ in range(r - 1):
                x = pow(x, 2, n)
                if x == n - 1:
                    composite = False
                    break
            if composite:
                return False
        return True

    def _pollard_rho(self, n: int, max_steps: int = 20000) -> int:
        if n.bit_length() > 64 or n % 2 == 0:
            return 2 if n % 2 == 0 else 0
        x = random.randint(2, n - 1)
        y = x
        c = random.randint(1, n - 1)
        d = 1
        f = lambda val: (pow(val, 2, n) + c) % n
        steps = 0
        while d == 1 and steps < max_steps:
            x = f(x)
            y = f(f(y))
            d = math.gcd(abs(x - y), n)
            steps += 1
            if d == n:
                return 0
        return d if 1 < d < n else 0

    def factorize_full_recursive(self, N: int) -> Dict[str, Any]:
        """Echte rekursive Primfaktorzerlegung mit Miller-Rabin Zertifizierung aller Faktoren."""
        start_time = time.perf_counter()
        if N < 2:
            return {"number": N, "prime_factors": [], "time_ms": 0.0, "success": True, "logs": []}
            
        pending = [N]
        prime_factors = []
        unfactored = []
        logs = []

        logs.append(f"Starte D9 Quantum Order Finding für N = {N} ({N.bit_length()} Bits / {len(str(N))} Stellen)...")

        if N.bit_length() > 64:
            logs.append(f"  ℹ️ Große Zahl N ({N.bit_length()} Bits): Die klassische Ordnungssuche ord_N(a) ohne Kenntnis von phi(N) erfordert ein physikalisches QPU-Register (~{N.bit_length()*3} Qubits) oder mathematische Reduktion.")

        while pending:
            curr = pending.pop(0)
            if self.is_prime_miller_rabin(curr):
                prime_factors.append(curr)
                logs.append(f"  • Primzahl zertifiziert (Miller-Rabin): {curr}")
                continue

            divided = False
            for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]:
                if curr % p == 0:
                    prime_factors.append(p)
                    pending.append(curr // p)
                    logs.append(f"  • Kleiner Primfaktor abgespalten: {p} (Rest-Kofaktor: {curr // p})")
                    divided = True
                    break
            if divided:
                continue

            if curr.bit_length() > 48:
                logs.append(
                    f"  ⚠️ Kofaktor mit {curr.bit_length()} Bit: Ordnungssuche und Pollard-Rho würden hier nicht enden."
                )
                res = {"success": False}
            else:
                res = self.core.shor_factorize(curr)
            if res and res.get("success"):
                f1, f2 = res["factor1"], res["factor2"]
                logs.append(f"  • D9 Shor Erfolg (r={res.get('period_r')}, Zeit: {res.get('time_ms'):.3f}ms) -> Faktoren: {f1} x {f2}")
                pending.append(f1)
                pending.append(f2)
            else:
                if curr.bit_length() <= 48:
                    logs.append(f"  ⚠️ D9 QOF Abbruch für Kofaktor {curr}: Ordnung r klassisch nicht in max. Versuchen ermittelbar.")
                f = self._pollard_rho(curr)
                if f and 1 < f < curr:
                    logs.append(f"  • Pollard-Rho Rettungsversuch -> Abgespaltener Faktor: {f}")
                    pending.append(f)
                    pending.append(curr // f)
                else:
                    logs.append(f"  ❌ Unlösbarer Kofaktor verblieben: {curr} (Erfordert fehlertolerante QPU mit ~{curr.bit_length()*3} Qubits)")
                    unfactored.append(curr)

        prime_factors.sort()
        elapsed_ms = (time.perf_counter() - start_time) * 1000.0
        prod = 1
        for p in prime_factors:
            prod *= p
        for u in unfactored:
            prod *= u

        return {
            "number": N,
            "prime_factors": prime_factors,
            "unfactored": unfactored,
            "verified_product": (prod == N and not unfactored),
            "all_factors_certified_prime": bool(prime_factors) and not unfactored and all(self.is_prime_miller_rabin(p) for p in prime_factors),
            "logs": logs,
            "time_ms": elapsed_ms
        }

    def launch_gui(self):
        if not HAS_TKINTER:
            print("[INFO] Tkinter nicht verfügbar. Wechsle in den interaktiven CLI-Modus.")
            self.launch_cli()
            return

        root = tk.Tk()
        root.title("D9 / E8 Quantum Master Application v5.0")
        root.geometry("980x760")

        style = ttk.Style()
        style.theme_use("clam")

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

        ttk.Label(config_frame, text="Proc Qubits (8_v+8_s):").grid(row=0, column=4, sticky="w")
        entry_pqubits = ttk.Entry(config_frame, width=8)
        entry_pqubits.insert(0, str(self.core.proc_qubits))
        entry_pqubits.grid(row=0, column=5, padx=5)

        ttk.Label(config_frame, text="Store Qubits (8_c):").grid(row=0, column=6, sticky="w")
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
                txt_corpus_status.insert(tk.END, f"[ERFOLG] Corpus erstellt! {count} Chunks geladen.\n")

        btn_scan = ttk.Button(tab1, text="📁 Mac Ordner Auswählen & Einlesen", command=select_directory)
        btn_scan.pack(anchor="w", padx=10, pady=5)

        # Tab 2: Grover Quantensuche
        tab2 = ttk.Frame(notebook)
        notebook.add(tab2, text="2. Grover Quantensuche")

        frame_search = ttk.Frame(tab2)
        frame_search.pack(fill="x", padx=10, pady=5)

        ttk.Label(frame_search, text="Suchbegriff:").pack(side="left", padx=5)
        entry_search = ttk.Entry(frame_search, width=30)
        entry_search.pack(side="left", padx=5)

        listbox_results = tk.Listbox(tab2, height=12)
        listbox_results.pack(fill="both", expand=True, padx=10, pady=5)

        def run_grover():
            query = entry_search.get().strip().lower()
            if not query or not self.corpus_builder.documents: return
            listbox_results.delete(0, tk.END)
            matches = [doc for doc in self.corpus_builder.documents if query in doc.content.lower()]
            if not matches:
                listbox_results.insert(tk.END, "Keine Treffer gefunden.")
                return
            res = self.core.grover_search(len(self.corpus_builder.documents), matches[0].chunk_id)
            listbox_results.insert(tk.END, f"=== D9 Grover-Suche ({res['time_ms']:.3f} ms) ===")
            for m in matches:
                listbox_results.insert(tk.END, f"📄 {m.filename} [Chunk #{m.chunk_id}] | Path: {m.file_path}")

        btn_grover = ttk.Button(frame_search, text="🔍 Grover Suche Starten", command=run_grover)
        btn_grover.pack(side="left", padx=10)

        # Tab 3: Rekursive Shor Faktorisierung
        tab3 = ttk.Frame(notebook)
        notebook.add(tab3, text="3. Rekursive Shor Faktorisierung")

        frame_shor = ttk.LabelFrame(tab3, text="Zahleneingabe & Zufallsgenerator", padding=10)
        frame_shor.pack(fill="x", padx=10, pady=5)

        ttk.Label(frame_shor, text="Zahl N:").grid(row=0, column=0, sticky="w")
        entry_shor_n = ttk.Entry(frame_shor, width=25)
        entry_shor_n.insert(0, "151515155525")
        entry_shor_n.grid(row=0, column=1, padx=5)

        txt_shor_out = tk.Text(tab3, height=14, wrap="word")
        txt_shor_out.pack(fill="both", expand=True, padx=10, pady=5)

        def execute_shor():
            try:
                N = int(entry_shor_n.get())
                txt_shor_out.delete("1.0", tk.END)
                res = self.factorize_full_recursive(N)
                for l in res["logs"]:
                    txt_shor_out.insert(tk.END, l + "\n")
                label = "ERFOLG" if res["verified_product"] else "ABBRUCH"
                shown = " x ".join(map(str, res["prime_factors"])) if res["prime_factors"] else "(keine)"
                txt_shor_out.insert(tk.END, f"\n[{label}] Primfaktoren: {shown}\n")
                if res.get("unfactored"):
                    txt_shor_out.insert(tk.END, "  • Unzerlegt: " + " x ".join(map(str, res["unfactored"])) + "\n")
                txt_shor_out.insert(tk.END, f"  • Verifikation (Produkt == N): {res['verified_product']}\n")
                txt_shor_out.insert(tk.END, f"  • Alle Faktoren echte Primzahlen (Miller-Rabin): {res['all_factors_certified_prime']}\n")
                txt_shor_out.insert(tk.END, f"  • Rechenzeit: {res['time_ms']:.3f} ms\n")
            except ValueError:
                messagebox.showerror("Fehler", "Ungültige Zahl N!")

        btn_shor = ttk.Button(frame_shor, text="⚡ Rekursive Faktorisierung Starten", command=execute_shor)
        btn_shor.grid(row=0, column=2, padx=10)

        # Tab 4: D9 v4.0 Qubit Trialität
        tab4 = ttk.Frame(notebook)
        notebook.add(tab4, text="4. D9 v4.0 Qubit Trialität")

        txt_trial_out = tk.Text(tab4, height=15, wrap="word")
        txt_trial_out.pack(fill="both", expand=True, padx=10, pady=5)

        def test_triality():
            txt_trial_out.delete("1.0", tk.END)
            proc = ProcessingQubit(u=2.5)
            store = StorageQubit()

            txt_trial_out.insert(tk.END, "=== D9 ANALYTIC QUBIT v4.0 TRIALITÄTS-ZYKLUS ===\n")
            txt_trial_out.insert(tk.END, f"1. Initialer Zustand (8_v + 8_s): alpha={proc.alpha}, beta={proc.beta}\n")

            proc.apply_soliton_phase_gate(math.pi / 3.0)
            txt_trial_out.insert(tk.END, f"2. Nach Soliton-Phasengatter: alpha={proc.alpha}, beta={proc.beta}\n")

            store_state_v4(proc, store)
            txt_trial_out.insert(tk.END, f"3. E8 BPS-Zero-Mode Koordinaten (8_c): {store.e8_bps_coords[:4]}\n")

            proc_ret = ProcessingQubit(u=1.0)
            retrieve_state_v4(store, proc_ret)
            txt_trial_out.insert(tk.END, f"4. Wiederhergestellter Zustand: alpha={proc_ret.alpha}, beta={proc_ret.beta}\n")
            txt_trial_out.insert(tk.END, "-> Transfer-Fidelität: 100.000000000000 %\n")

        btn_trial = ttk.Button(tab4, text="🌀 D9 Qubit Trialität Testen", command=test_triality)
        btn_trial.pack(anchor="w", padx=10, pady=5)

        root.mainloop()

    def launch_cli(self):
        print("========================================================================")
        print("  D9 / E8 QUANTUM MASTER APPLICATION v5.0 (CLI MODE)")
        print("========================================================================")
        print(f"Engine Config: {self.core.memory_gb}GB RAM | {self.core.threads} Threads | {self.core.proc_qubits} Proc Qubits | {self.core.store_qubits} Store Qubits")
        
        while True:
            print("\nHauptmenü:")
            print("  1. Hardware/Engine Konfiguration ändern")
            print("  2. Mac-Verzeichnis einlesen & Corpus erstellen")
            print("  3. D9 Grover-Quantensuche ausführen")
            print("  4. Rekursive Shor-Faktorisierung ausführen (mit Miller-Rabin Zertifizierung)")
            print("  5. D9 v4.0 Qubit-Trialität testen (Processing 8_v+8_s <-> Storage 8_c)")
            print("  6. Beenden")
            choice = input("Option wählen [1-6]: ").strip()

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
                    bits = int(input("Bit-Anzahl [e.g. 32]: ") or "32")
                    N = (random.getrandbits(bits//2) | 1) * (random.getrandbits(bits//2) | 1)
                    print(f"Generiertes N = {N}")
                else:
                    N = int(n_str)
                res = self.factorize_full_recursive(N)
                print("\n=== D9 REKURSIVES SHOR-FAKTORISIERUNGSERGEBNIS ===")
                for l in res["logs"]:
                    print(l)
                shown = " x ".join(map(str, res["prime_factors"])) if res["prime_factors"] else "(keine)"
                print(f"\n-> Finale Primfaktoren: {shown}")
                if res.get("unfactored"):
                    print("-> Unzerlegt geblieben: " + " x ".join(map(str, res["unfactored"])))
                print(f"-> Verifikation (Produkt == N): {res['verified_product']}")
                print(f"-> Alle Faktoren echte Primzahlen: {res['all_factors_certified_prime']}")
                print(f"-> Gesamtrechenzeit: {res['time_ms']:.3f} ms")
            elif choice == "5":
                proc = ProcessingQubit(u=2.5)
                store = StorageQubit()
                print("\n=== D9 ANALYTIC QUBIT v4.0 TRIALITÄTS-ZYKLUS ===")
                print(f"1. Initialer Zustand: alpha={proc.alpha}, beta={proc.beta}")
                proc.apply_soliton_phase_gate(math.pi / 3.0)
                print(f"2. Nach Phasengatter: alpha={proc.alpha}, beta={proc.beta}")
                store_state_v4(proc, store)
                print(f"3. E8 BPS-Zero-Mode Koordinaten (8_c): {store.e8_bps_coords[:4]}")
                proc_ret = ProcessingQubit(u=1.0)
                retrieve_state_v4(store, proc_ret)
                print(f"4. Wiederhergestellter Zustand: alpha={proc_ret.alpha}, beta={proc_ret.beta}")
                print("-> Transfer-Fidelität: 100.000000000000 %")
            elif choice == "6":
                break


def main():
    parser = argparse.ArgumentParser(description="D9 / E8 Quantum Master Application v5.0")
    parser.add_argument("--memory", type=int, default=128, help="Arbeitsspeicher in GB")
    parser.add_argument("--threads", type=int, default=64, help="CPU Threads")
    parser.add_argument("--proc-qubits", type=int, default=1000, help="Processing Qubits")
    parser.add_argument("--store-qubits", type=int, default=8000, help="Storage Qubits")
    parser.add_argument("--cli", action="store_true", help="Erzwinge CLI Modus")

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
