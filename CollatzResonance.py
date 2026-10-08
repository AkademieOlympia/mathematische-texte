import math
from collections import defaultdict


def collatz(n, steps):
    orbit = []
    for _ in range(steps):
        orbit.append(n)
        if n % 2 == 0:
            n //= 2
        else:
            n = 3*n + 1
    orbit.append(n)
    return orbit


def parity_signature(orbit):
    return [x % 2 for x in orbit[:-1]]


def gap_signature(bits):
    """
    Abstand zwischen ungeraden Ereignissen
    """
    pos = [i for i,b in enumerate(bits) if b==1]

    gaps=[]
    for i in range(len(pos)-1):
        gaps.append(pos[i+1]-pos[i])

    return gaps


def resonance(gaps):
    """
    Analogie:
    rho = |M 1|/sqrt(k)

    Hier M als Gap-Differenzmatrix
    """

    k=len(gaps)

    if k<2:
        return 0

    s=sum(gaps)

    # Zentrum entfernen
    mean=s/k

    M=[]

    for g in gaps:
        M.append(g-mean)

    norm=sum(x*x for x in M)

    return math.sqrt(norm/k)


def collatz_class(n,steps=100):

    orbit=collatz(n,steps)

    bits=parity_signature(orbit)

    gaps=gap_signature(bits)

    rho=resonance(gaps)

    return rho, tuple(gaps)


def test(N=10000):

    atlas=defaultdict(list)

    for n in range(2,N):

        rho,gap=collatz_class(n)

        key=round(rho,10)

        atlas[key].append(n)


    collisions=0

    for r,values in atlas.items():

        if len(values)>1:
            collisions+=1


    print("Collatz Resonanzatlas")
    print("---------------------")
    print("N =",N)
    print("Klassen =",len(atlas))
    print("Kollisionen =",collisions)


    print("\nBeispiele:")
    for r,v in list(atlas.items())[:10]:
        print(r,v[:10])



if __name__=="__main__":
    test(10000)