# Python Verification Script for Fibonacci Anyon Hadamard Braid Word
import numpy as np

phi_inv = (np.sqrt(5) - 1) / 2
F = np.array([[phi_inv, np.sqrt(phi_inv)], [np.sqrt(phi_inv), -phi_inv]], dtype=np.complex128)
R = np.array([[np.exp(-1j * 4 * np.pi / 5), 0], [0, np.exp(1j * 3 * np.pi / 5)]], dtype=np.complex128)

s1 = R
s2_inv = np.linalg.inv(F @ R @ F)

word_ops = [s1, s2_inv, s1, s2_inv, s1, s2_inv @ s2_inv @ s2_inv, s1, s2_inv, s1, s2_inv, s1]
U_w = np.eye(2, dtype=np.complex128)
for op in word_ops:
    U_w = U_w @ op

H = (1 / np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=np.complex128)
fidelity = 0.5 * np.abs(np.trace(U_w.conj().T @ H))

print(f"Braid Word: s1 s2^-1 s1 s2^-1 s1 s2^-3 s1 s2^-1 s1 s2^-1 s1")
print(f"Verified Fidelity F = {fidelity:.8f} ({fidelity*100:.4f}%)")
