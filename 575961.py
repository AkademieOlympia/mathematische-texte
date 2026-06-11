import numpy as np

def is_prime(n):
    if n < 2: return False
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0: return False
    return True

def metric_seal_eabc(a, b, c):
    """
    Operatives Protokoll zur Identifikation der Packungs-Grenzwerte.
    Berechnet die Kondensation der Primzahl 'e' aus dem triadischen Feld (a,b,c).
    """
    # I. KONSTANTEN DER ORDNUNG (Die 380.000er Schranke)
    T_CRIT = 380000
    C_EABC = 2000 * np.pi  # Die Druck-Konstante (6283.18)
    
    # II. DIE OKTONISCHE DEFEKT-FREQUENZ
    # nu_defekt entspricht der 'Taktfrequenz' der topologischen Korrektur.
    nu_defekt = T_CRIT / C_EABC  # ~60.48 Hz (relativ zur Spektraldichte)
    
    # III. DAS TRIADISCHE FENSTER (Packungs-Shift 57-59-61)
    # Prüfung, ob die Frequenz in der Resonanz-Zone der 'jüdischen Kusszahl' liegt.
    TRIAD = [57, 59, 61]
    is_stable = 57 <= nu_defekt <= 61 # Resonanz-Check
    
    # IV. DIE 1/3 INVERSION (Topologische Linse)
    # Transformation des Raumes auf die Innenkugel zur Bertrand-Fokussierung.
    target_e_approx = (a * b * c) / 3.0
    
    # V. DER 180°-UMSCHLAG (Kondensations-Suche)
    # Wir suchen die stabilisierende Primzahl 'e' (Familie E), 
    # die das G2-Invariante Quadrat n = a*b*c*e = K^2 vervollständigt.
    results = []
    search_window = int(target_e_approx * 0.5)
    
    for e_cand in range(max(2, int(target_e_approx - search_window)), int(target_e_approx + search_window)):
        if is_prime(e_cand):
            n_total = a * b * c * e_cand
            k_val = np.sqrt(n_total)
            
            if k_val == int(k_val): # Quadrat-Bedingung erfüllt
                results.append({
                    'a': a, 'b': b, 'c': c, 'e': e_cand,
                    'n': n_total, 'K': int(k_val),
                    'nu': round(nu_defekt, 4),
                    'shift': "HEXAGONAL -> QUATERNIONAL",
                    'status': "CRYSTALLIZED" if is_stable else "DECOHERENT"
                })
    return results

# Beispiel: Initialisierung der Pulsation mit einem stabilen Triplett
# eabc_crystal = metric_seal_eabc(a_prim, b_prim, c_prim)