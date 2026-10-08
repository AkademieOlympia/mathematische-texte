# ordinal_rank_explorer_1M.sage
# Großflächiger Ordinal-Rang-Explorer bis o = 1,000,000
# Testet lexikographische Ränge unter Einbeziehung von akkumulierten Exponentensummen.

from collections import Counter
import itertools

def v2(n):
    n = ZZ(n)
    return n.valuation(2)

def collatz_step_with_acc(o, s_acc):
    o = ZZ(o)
    val = 3 * o + 1
    k = v2(val)
    nxt_o = val // (2**k)
    nxt_s_acc = s_acc + k
    return nxt_o, nxt_s_acc

def evaluate_rank_state(o, s_acc, cand):
    o = ZZ(o)
    if o <= 1:
        return (ZZ(0), ZZ(0), ZZ(0))
    fA, fB, fC = cand
    return (fA(o), fB(o, s_acc), fC(o))

# --- Bausteine ---
def scale_bit_length(o):
    return ZZ(o.nbits() - 1)

def ast_current_k(o, s_acc):
    return v2(3 * o + 1)

def ast_accumulated_s(o, s_acc):
    return ZZ(s_acc)

def rest_mod_12(o):
    return o % 12

def rest_mod_8(o):
    return o % 8

def failure_type(R0, R1):
    A0, B0, C0 = R0
    A1, B1, C1 = R1
    if A1 > A0:
        return "A grows"
    if A1 == A0 and B1 > B0:
        return "A same, B grows"
    if A1 == A0 and B1 == B0 and C1 >= C0:
        return "A,B same, C not smaller"
    return "other"

def test_orbit_with_buffer(cand_name, cand, max_val=1000000, max_steps=150):
    total_orbits = 0
    failures = 0
    first_failure_at = None
    kinds = Counter()
    
    for start_o in range(3, max_val + 1, 2):
        total_orbits += 1
        o = ZZ(start_o)
        s_acc = 0
        
        for _ in range(max_steps):
            if o <= 1:
                break
                
            R0 = evaluate_rank_state(o, s_acc, cand)
            nxt_o, nxt_s_acc = collatz_step_with_acc(o, s_acc)
            R1 = evaluate_rank_state(nxt_o, nxt_s_acc, cand)
            
            if not (R1 < R0):
                failures += 1
                kind = failure_type(R0, R1)
                kinds[kind] += 1
                if first_failure_at is None:
                    first_failure_at = start_o
                break
                
            o, s_acc = nxt_o, nxt_s_acc
            
    success_rate = 1.0 - (RR(failures) / RR(total_orbits))
    return {
        "name": cand_name,
        "orbits": total_orbits,
        "failures": failures,
        "success_rate": success_rate,
        "first_failure": first_failure_at,
        "kinds": kinds
    }

def run_large_explorer():
    max_val = 1000000
    print(f"=== Großflächiger Ordinal-Rang-Explorer bis o = {max_val} ===")
    print("Starte Auswertung...\n")
    
    A_list = [("nbits-1", scale_bit_length)]
    B_list = [("k = v2(3o+1)", ast_current_k), ("S_acc (kumuliert)", ast_accumulated_s)]
    C_list = [("mod 12", rest_mod_12), ("mod 8", rest_mod_8)]
    
    results = []
    
    for (a_name, a_fn), (b_name, b_fn), (c_name, c_fn) in itertools.product(A_list, B_list, C_list):
        cand_name = f"A:{a_name} | B:{b_name} | C:{c_name}"
        cand_tuple = (a_fn, b_fn, c_fn)
        
        print(f"Teste Kandidat: {cand_name} ...", flush=True)
        res = test_orbit_with_buffer(cand_name, cand_tuple, max_val=max_val)
        results.append(res)
        print(f"  -> Erster Fehler bei o = {res['first_failure']} | Fehlerausfälle: {res['failures']} / {res['orbits']}\n")
        
    results.sort(key=lambda x: (x['first_failure'] is not None, x['first_failure'] or 0), reverse=True)
    
    print("="*65)
    print(f"FINALES RANKING BIS o = {max_val}")
    print("="*65)
    for i, r in enumerate(results, 1):
        ff = r['first_failure'] if r['first_failure'] is not None else "Keine!"
        print(f"{i:2d}. {r['name']}")
        print(f"    Erster Fehler bei o = {ff} | Ausfälle: {r['failures']} / {r['orbits']}")
        print(f"    Fehlertypen: {dict(r['kinds'])}\n")

if __name__ == "__main__":
    run_large_explorer()