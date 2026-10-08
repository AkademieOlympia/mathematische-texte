from sympy import isprime

F = {1: [], 5: [], 7: [], 11: []}
p = 5
while min(map(len, F.values())) < 100:
    if isprime(p) and p % 12 in F and len(F[p % 12]) < 100:
        F[p % 12].append(p)
    p += 1

for r, values in F.items():
    assert len(values) == 100
    assert all(isprime(q) and q % 12 == r for q in values)
    print(r, values)
