import numpy as np
import matplotlib.pyplot as plt

zeros = np.load("zeros6.npy")

# Comparison: Around Carmichael 561 vs Prime 541 (nearby prime)
c_num = 561
p_num = 541

x_range_c = np.linspace(c_num - 10, c_num + 10, 400)
x_range_p = np.linspace(p_num - 10, p_num + 10, 400)

def zeta_signal_vec(x_array, gamma_list):
    signals = []
    for x in x_array:
        signals.append(np.sum(np.cos(gamma_list * np.log(x))))
    return np.array(signals)

# Use 50k zeros for faster signal processing
gamma_fast = zeros[:50000]

sig_c = zeta_signal_vec(x_range_c, gamma_fast)
sig_p = zeta_signal_vec(x_range_p, gamma_fast)

plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
plt.plot(x_range_p, sig_p, color='red', label='Signal around Prime 541')
plt.axvline(x=p_num, color='black', linestyle='--', label='Prime Position')
plt.title('Topological Singularity (True Prime)')
plt.ylabel('Zeta Resonance')
plt.grid(True, alpha=0.3)
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(x_range_c, sig_c, color='blue', label='Signal around Carmichael 561')
plt.axvline(x=c_num, color='black', linestyle='--', label='Carmichael Position')
plt.title('Topological Glideway (Almost Prime)')
plt.grid(True, alpha=0.3)
plt.legend()

plt.tight_layout()
plt.savefig('prime_vs_carmichael_landscape.png')
print("Visual comparison generated.")