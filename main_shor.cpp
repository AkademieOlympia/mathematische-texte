#include <iostream>
#include <chrono>
#include "e8_qram_shor_avx512_kernel.hpp"

int main() {
    // Beispiel: Faktorisierung von N = 10403 (103 x 101)
    uint64_t N = 10403;
    uint64_t a = 2;

    std::cout << "=========================================================\n";
    std::cout << "  Nativer C++20 E8-QRAM Shor-Kernel (Apple Silicon NEON)\n";
    std::cout << "=========================================================\n";
    std::cout << "Faktorisiere N = " << N << " mit Generator a = " << a << "...\n";

    auto start = std::chrono::high_resolution_clock::now();
    
    // Initialisierung und Ausführung des C++20 Kernels
    d9_shor::E8QRAMShorAVX512Kernel<double> kernel(N, a);
    auto [p, q] = kernel.run_factorization();

    auto end = std::chrono::high_resolution_clock::now();
    double elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();

    if (p > 0 && q > 0) {
        std::cout << "--> ERFOLG: N = " << N << " -> p = " << p << ", q = " << q << "\n";
        std::cout << "--> Laufzeit: " << elapsed_ms << " ms\n";
    } else {
        std::cout << "--> Kein voller Periodentreffer für a = " << a << " (Anderen Generator wählen)\n";
    }

    return 0;
}