import math

def eabc_class(p):
    return {5: "A", 7: "B", 11: "C", 1: "E"}.get(p % 12, "—")

primes = [5, 7, 11, 13, 17]

print("p  Klasse  a_p             E_main  E_side_each  E_side_total  E_total")
for p in primes:
    a = p / math.sqrt(3)
    print(f"{p:2d}   {eabc_class(p):>1}     {a:10.6f}      {p**2:5d}    {p**2/3:10.6f}      {2*p**2:5d}       {3*p**2:5d}")
