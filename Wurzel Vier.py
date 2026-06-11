def run_experiment(N, eps, runs=3):
    results = []
    
    for i in range(runs):
        print(f"\nRun {i+1} bei N={N}")
        
        mask = sieve(N)
        quads = find_quadruplets(mask)
        
        print("Vierlinge:", len(quads))
        
        H = build_operator(quads, eps)
        
        ev = np.linalg.eigvalsh(H)
        sp = np.diff(np.sort(ev))
        sp /= np.mean(sp)
        
        v = np.var(sp)
        k = ks(sp, z_sp)
        
        print("Var:", v)
        print("KS :", k)
        
        results.append((v, k))
    
    print("\n--- Mittelwerte ---")
    print("Var mean:", np.mean([r[0] for r in results]))
    print("KS  mean:", np.mean([r[1] for r in results]))
    
    return results