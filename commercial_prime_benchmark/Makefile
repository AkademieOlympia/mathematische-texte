# Commercial next-10 primes benchmark: GMP / OpenSSL vs Prime-Peano
# Same binary, same seeds, same hardware.

CXX      ?= clang++
GMP_PREFIX ?= /opt/homebrew
SSL_PREFIX ?= /opt/homebrew/opt/openssl@3

CXXFLAGS += -std=c++17 -O3 -march=native
CXXFLAGS += -Wall -Wextra -Wno-deprecated-literal-operator
CXXFLAGS += -I$(GMP_PREFIX)/include -I$(SSL_PREFIX)/include
LDFLAGS  += -L$(GMP_PREFIX)/lib -L$(SSL_PREFIX)/lib
LDFLAGS  += -lgmpxx -lgmp -lcrypto
LDFLAGS  += -Wl,-rpath,$(GMP_PREFIX)/lib -Wl,-rpath,$(SSL_PREFIX)/lib

.PHONY: all clean run run-quick run-eabc-quick

all: bench_next10 Primes_EABC

bench_next10: bench_next10.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

Primes_EABC: Primes_EABC.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

run: bench_next10
	./bench_next10 --trials 5 --bits 512,1024,1536,2048,3072,4096

run-quick: bench_next10
	./bench_next10 --trials 1 --bits 512,1024

run-eabc-quick: Primes_EABC
	./Primes_EABC --bench --trials 1 --bits 512 --out out_eabc_quick

clean:
	rm -f bench_next10 Primes_EABC
	rm -rf out out_eabc out_eabc_quick
