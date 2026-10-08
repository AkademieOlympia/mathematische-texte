import numpy as np

# ============================================================
# GATE K: SMOOTH-NUMBER SCALE JUMP AUDIT (v1.0 Final)
# Inklusive expliziter Nestedness-Prüfung, Kopplungsblöcken V_k, 
# Eigenmode-Reorganisation (Theta_k) und verschachtelten Nullmodellen.
# ============================================================

def build_smooth_operators(x_pos, tol=1e-12):
    """
    Baut beide normalisierten Interferenzoperatoren K_S und K_C 
    für eine gegebene glatte Zahlenmenge inklusive strenger Invarianten-Prüfungen.
    """
    x_pos = np.asarray(x_pos, dtype=float)
    if len(x_pos) < 2:
        raise ValueError("Too few modes for operator construction.")

    Xi = x_pos[:, None]
    Xj = x_pos[None, :]

    G_minus = np.sinc(2.0 * (Xi - Xj))
    G_plus  = np.sinc(2.0 * (Xi + Xj))

    S = 0.5 * (G_minus - G_plus)
    C = 0.5 * (G_minus + G_plus)

    def normalize_gram(G, name):
        sym_err = np.linalg.norm(G - G.T, "fro")
        if sym_err > tol:
            raise ValueError(f"{name}: Gram matrix not symmetric (error: {sym_err:.4e}).")

        eig_G = np.linalg.eigvalsh(G)
        if eig_G.min() < -tol:
            raise ValueError(f"{name}: Gram matrix violates positive semidefiniteness.")

        d = np.diag(G).copy()
        if np.any(d <= tol):
            raise ValueError(f"{name}: Nonpositive Gram diagonal detected.")

        inv_sqrt = 1.0 / np.sqrt(d)
        R = inv_sqrt[:, None] * G * inv_sqrt[None, :]
        K = R - np.eye(len(G))

        diag_err = np.max(np.abs(np.diag(R) - 1.0))
        trace_err = abs(np.trace(K))

        eig_K, vec_K = np.linalg.eigh(K)
        if eig_K[0] < -1.0 - tol:
            raise ValueError(f"{name}: lambda_min(K) < -1 violation.")

        kappa_max = float(eig_K[-1])
        v_max = vec_K[:, -1]
        v_max = v_max / np.linalg.norm(v_max)

        # Spektrallücke delta_top
        kappa_2 = float(eig_K[-2]) if len(eig_K) > 1 else np.nan
        spectral_gap = kappa_max - kappa_2 if not np.isnan(kappa_2) else np.nan

        ipr = float(np.sum(np.abs(v_max)**4))
        n_eff = float(1.0 / ipr) if ipr > 0 else np.nan

        return {
            "G": G,
            "R": R,
            "K": K,
            "eig_K": eig_K,
            "vec_K": vec_K,
            "kappa_max": kappa_max,
            "kappa_2": kappa_2,
            "spectral_gap": spectral_gap,
            "v_max": v_max,
            "ipr": ipr,
            "N_eff": n_eff,
            "diag_R_error": float(diag_err),
            "trace_K_error": float(trace_err)
        }

    return {
        "S": normalize_gram(S, "S"),
        "C": normalize_gram(C, "C")
    }


def nested_indices(old_nums, new_nums):
    """
    Findet die Indizes der alten glatten Zahlen in der neuen, 
    sortierten Menge und prüft die exakte Verschachtelung.
    """
    old_nums = np.asarray(old_nums, dtype=float)
    new_nums = np.asarray(new_nums, dtype=float)

    idx = np.searchsorted(new_nums, old_nums)

    if np.any(idx >= len(new_nums)):
        raise ValueError("Old smooth set out of bounds in new set.")

    if not np.allclose(new_nums[idx], old_nums, atol=1e-12):
        raise ValueError("Smooth-number nestedness invariant violated.")

    return idx


def extension_blocks(K_new, old_idx):
    """
    Zerlegt den neuen Operator in Hauptuntermatrix (old-old), 
    Kopplungsblock V (old-new) und new-new Block W.
    """
    N = K_new.shape[0]
    old_idx = np.asarray(old_idx, dtype=int)

    mask_new = np.ones(N, dtype=bool)
    mask_new[old_idx] = False
    new_idx = np.where(mask_new)[0]

    K_old_block = K_new[np.ix_(old_idx, old_idx)]
    V = K_new[np.ix_(old_idx, new_idx)]
    W = K_new[np.ix_(new_idx, new_idx)]

    return K_old_block, V, W, new_idx


def eigenmode_overlap(v_old, v_new, old_idx):
    """
    Berechnet den Überlapp der in den neuen Raum eingebetteten 
    alten Eigenmode mit der neuen führenden Eigenmode.
    """
    embedded = np.zeros_like(v_new)
    embedded[old_idx] = v_old
    
    norm_emb = np.linalg.norm(embedded)
    if norm_emb > 0:
        embedded /= norm_emb
        
    v_new_norm = v_new / np.linalg.norm(v_new)

    overlap = float(abs(np.vdot(embedded, v_new_norm))**2)
    return overlap, 1.0 - overlap  # Theta_k = 1 - O_k


def nested_uniform_surrogate(x_min, x_max, N_old, N_new, rng):
    """
    K_0 (nested uniform-log point surrogate): 
    Erzeugt exakt verschachtelte Nullpunktmengen.
    """
    old = rng.uniform(x_min, x_max, size=N_old)
    added = rng.uniform(x_min, x_max, size=N_new - N_old)
    
    new_set = np.concatenate([old, added])
    old.sort()
    new_set.sort()
    return old, new_set


def generate_smooth_numbers(X, B_max_prime):
    """
    Generiert B-glatte Zahlen <= X.
    """
    def get_primes_up_to(limit):
        if limit < 2: return np.array([], dtype=int)
        sieve = np.ones(limit + 1, dtype=bool)
        sieve[:2] = False
        for k in range(2, int(np.sqrt(limit)) + 1):
            if sieve[k]: sieve[k*k::k] = False
        return np.nonzero(sieve)[0]

    primes = get_primes_up_to(B_max_prime)
    if len(primes) == 0:
        return np.array([]), np.array([]), np.array([])

    smooth_list = []
    def backtrack(current_val, prime_idx):
        if current_val > X:
            return
        smooth_list.append(current_val)
        for i in range(prime_idx, len(primes)):
            p = primes[i]
            if current_val > X / p:
                break
            backtrack(current_val * p, i)

    backtrack(1, 0)
    smooth_arr = np.array(sorted(list(set(smooth_list))), dtype=float)
    x_pos = np.log(smooth_arr)
    
    mask_valid = (smooth_arr > 3.0)
    return smooth_arr[mask_valid], x_pos[mask_valid], primes


def gate_K_smooth_scale_audit(X=5000.0, prime_list=(2, 3, 5, 7, 11, 13, 17), num_surrogates=200, seed=2026):
    """
    Gate K: Smooth-Number Scale Jump Audit mit Nestedness, Kopplungsblöcken V_k,
    Eigenmode-Reorganisation (Theta_k) und verschachtelten Nullmodellen.
    """
    print("=====================================================================")
    print(" GATE K: SMOOTH-NUMBER SCALE JUMP AUDIT (A-J v1.0 Final)")
    print("=====================================================================")
    print(f" Obergrenze X = {X}, Schrankenreihe B = {list(prime_list)}")

    rng = np.random.default_rng(seed)
    previous = None
    trajectory = []

    for B in prime_list:
        smooth_nums, x_pos, primes_used = generate_smooth_numbers(X, B)
        N_k = len(x_pos)
        
        if N_k < 20:
            print(f"  [B = {B}] Zu wenige Zustände (N = {N_k}). Überspringe.")
            continue

        current = build_smooth_operators(x_pos)

        if previous is None:
            previous = {
                "B": B, 
                "N": N_k, 
                "smooth_nums": smooth_nums, 
                "x": x_pos, 
                "result": current
            }
            print(f"  [Baseline B = {B}] Initialer Zustand etabliert (N = {N_k}).")
            continue

        N_old = previous["N"]
        delta_N = N_k - N_old

        # 1. Explizite Nestedness-Prüfung der Hauptuntermatrix
        old_idx = nested_indices(previous["smooth_nums"], smooth_nums)
        
        K_old_S = previous["result"]["S"]["K"]
        K_old_C = previous["result"]["C"]["K"]
        
        extracted_S, V_S, W_S, new_idx_S = extension_blocks(current["S"]["K"], old_idx)
        extracted_C, V_C, W_C, new_idx_C = extension_blocks(current["C"]["K"], old_idx)

        nest_err_S = np.linalg.norm(extracted_S - K_old_S, "fro")
        nest_err_C = np.linalg.norm(extracted_C - K_old_C, "fro")

        if nest_err_S > 1e-10 or nest_err_C > 1e-10:
            raise ValueError(f"Nestedness invariant violated at B={B} (Errors: S={nest_err_S:.2e}, C={nest_err_C:.2e})")

        # 2. Kopplungsblock-Metriken (Gamma_k & sigma_max)
        gamma_S = np.linalg.norm(V_S, "fro")**2 / (N_old * delta_N)
        gamma_C = np.linalg.norm(V_C, "fro")**2 / (N_old * delta_N)

        sigma_V_S = np.linalg.norm(V_S, 2)
        sigma_V_C = np.linalg.norm(V_C, 2)

        # 3. Spektralsprünge & Eigenmode-Reorganisation (Theta_k)
        delta_S = current["S"]["kappa_max"] - previous["result"]["S"]["kappa_max"]
        delta_C = current["C"]["kappa_max"] - previous["result"]["C"]["kappa_max"]

        overlap_S, theta_S = eigenmode_overlap(previous["result"]["S"]["v_max"], current["S"]["v_max"], old_idx)
        overlap_C, theta_C = eigenmode_overlap(previous["result"]["C"]["v_max"], current["C"]["v_max"], old_idx)

        delta_ipr_S = current["S"]["ipr"] - previous["result"]["S"]["ipr"]
        delta_ipr_C = current["C"]["ipr"] - previous["result"]["C"]["ipr"]

        # 4. Verschachtelte Nullmodell-Läufe (K_0: nested uniform-log point surrogate)
        null_delta_S = []
        null_delta_C = []
        null_gamma_S = []
        null_theta_S = []

        for _ in range(num_surrogates):
            x_old_sur, x_new_sur = nested_uniform_surrogate(
                x_pos[0], x_pos[-1], N_old, N_k, rng
            )
            try:
                old_sur = build_smooth_operators(x_old_sur)
                new_sur = build_smooth_operators(x_new_sur)
                
                # Null-Sprung
                null_delta_S.append(new_sur["S"]["kappa_max"] - old_sur["S"]["kappa_max"])
                null_delta_C.append(new_sur["C"]["kappa_max"] - old_sur["C"]["kappa_max"])
                
                # Null-Kopplung V (mithilfe einer äquivalenten Index-Verschachtelung)
                # Da x_new_sur sortiert ist und x_old_sur Teil davon ist:
                idx_sur = np.searchsorted(x_new_sur, x_old_sur)
                _, V_sur_S, _, _ = extension_blocks(new_sur["S"]["K"], idx_sur)
                null_gamma_S.append(np.linalg.norm(V_sur_S, "fro")**2 / (N_old * delta_N))
                
                # Null-Theta
                _, th_sur = eigenmode_overlap(old_sur["S"]["v_max"], new_sur["S"]["v_max"], idx_sur)
                null_theta_S.append(th_sur)
            except Exception:
                continue

        null_delta_S = np.array(null_delta_S)
        null_delta_C = np.array(null_delta_C)
        null_gamma_S = np.array(null_gamma_S)
        null_theta_S = np.array(null_theta_S)

        # Robuste Teststatistiken gegen das Nullmodell
        med_S, mad_S = np.median(null_delta_S), 1.4826 * np.median(np.abs(null_delta_S - np.median(null_delta_S))) if len(null_delta_S) > 0 else 1.0
        z_robust_S = (delta_S - med_S) / mad_S if mad_S > 0 else 0.0
        p_mc_S = (1 + np.count_nonzero(null_delta_S >= delta_S)) / (len(null_delta_S) + 1)

        print(f"\n[Schritt B = {B} | Delta N = {delta_N} | Nest-Err: S={nest_err_S:.1e}]")
        print(f"  S-Kanal: Delta kappa = {delta_S:+.6f} (Z_rob = {z_robust_S:+.2f}, p_mc = {p_mc_S:.4f})")
        print(f"           Gamma_S     = {gamma_S:.6f} (sigma_max(V) = {sigma_V_S:.4f})")
        print(f"           Reorg. Theta= {theta_S:.4f} (Overlap = {overlap_S:.4f})")
        print(f"  C-Kanal: Delta kappa = {delta_C:+.6f}")

        trajectory.append({
            "B": B,
            "delta_N": delta_N,
            "delta_S": delta_S,
            "delta_C": delta_C,
            "gamma_S": gamma_S,
            "gamma_C": gamma_C,
            "sigma_V_S": sigma_V_S,
            "sigma_V_C": sigma_V_C,
            "theta_S": theta_S,
            "theta_C": theta_C,
            "z_robust_S": z_robust_S,
            "p_mc_S": p_mc_S
        })

        previous = {
            "B": B, 
            "N": N_k, 
            "smooth_nums": smooth_nums, 
            "x": x_pos, 
            "result": current
        }

    print("=====================================================================")
    print(" GATE K AUDIT VOLLSTÄNDIG ABGESCHLOSSEN")
    print("=====================================================================")
    return trajectory

if __name__ == "__main__":
    traj = gate_K_smooth_scale_audit(X=3000.0, prime_list=(2, 3, 5, 7, 11, 13))

