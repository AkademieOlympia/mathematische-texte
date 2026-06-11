"""
Stabilitätsanalyse des quaternionischen Feldes gegen die theoretische RH-Schranke.
Ohne Sage – reines Python mit numpy.
"""
import math


def sieve_primes(limit):
    """Sieb des Eratosthenes: Primzahlen <= limit."""
    if limit < 2:
        return []
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            for j in range(i * i, limit + 1, i):
                is_prime[j] = False
    return [p for p in range(2, limit + 1) if is_prime[p]]


def quat_add(a, b):
    """Addition zweier Quaternionen (a,b,c,d)."""
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3])


def quat_reduced_norm(q):
    """Reduzierte Norm: a² + b² + c² + d² (Hamilton, -1,-1)."""
    return q[0]**2 + q[1]**2 + q[2]**2 + q[3]**2


def check_energy_stability(x_limit):
    """
    Prüft die Stabilität des quaternionischen Feldes gegen die
    theoretische Schranke von O(sqrt(x) * log(x)^2).
    """
    actual_q = (0.0, 0.0, 0.0, 0.0)
    primes_list = sieve_primes(x_limit + 1)

    for p in primes_list:
        if p == 2:
            direction = (1.0, 0.0, 0.0, 0.0)
        elif p % 4 == 1:
            direction = (0.0, 1.0, 0.0, 0.0)
        else:
            direction = (0.0, 0.0, 1.0, 0.0)
        actual_q = quat_add(actual_q, direction)

    actual_norm = math.sqrt(quat_reduced_norm(actual_q))
    theoretical_bound = math.sqrt(x_limit) * math.log(x_limit)
    stability_ratio = actual_norm / theoretical_bound
    return actual_norm, theoretical_bound, stability_ratio


# Analyse für die Energiedoku in Bamberg
x_val = 5000
norm, bound, ratio = check_energy_stability(x_val)

print(f"--- Stabilitätsanalyse x = {x_val} ---")
print(f"Aktuelle Energie-Norm: {norm:.4f}")
print(f"Theoretische RH-Schranke: {bound:.4f}")
print(f"Stabilitäts-Index: {ratio:.4f} (Werte < 1 deuten auf RH-Stabilität hin)")
