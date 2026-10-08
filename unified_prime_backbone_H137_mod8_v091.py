from math import isqrt

def qmul(x,y):
    a,b,c,d=x; e,f,g,h=y
    return (
        a*e-b*f-c*g-d*h,
        a*f+b*e+c*h-d*g,
        a*g-b*h+c*e+d*f,
        a*h+b*g-c*f+d*e,
    )

def qconj(x):
    a,b,c,d=x
    return (a,-b,-c,-d)

def qnorm(x):
    return sum(v*v for v in x)

E1=(0,1,0,0)

def is_prime(n):
    if n<2: return False
    if n%2==0: return n==2
    d=3
    while d*d<=n:
        if n%d==0: return False
        d+=2
    return True

def canonical_two_squares(N, require_odd_pair=False):
    lim=isqrt(N)
    best=None
    for s in range(1,lim+1):
        if require_odd_pair and s%2==0:
            continue
        t2=N-s*s
        if t2<s*s:
            break
        t=isqrt(t2)
        if t*t!=t2:
            continue
        if require_odd_pair and t%2==0:
            continue
        cand=(s,t)
        if best is None or cand<best:
            best=cand
    return best

def gauss_anchor(p):
    if p%4!=1:
        raise ValueError("p must be 1 mod 4")
    candidates=[]
    for a in range(0,isqrt(p)+1,2):
        b2=p-a*a
        if b2<0: break
        b=isqrt(b2)
        if b*b==b2 and b%2==1:
            candidates.append((b,a,0,0))
    if not candidates:
        raise RuntimeError(f"No Gauss anchor for p={p}")
    return min(candidates)

def inert_anchor_mod8(p):
    if p%4!=3:
        raise ValueError("p must be 3 mod 4")

    r = 0 if p%8==3 else 2
    m = p-r*r
    if m%8!=3:
        raise AssertionError("mod-8 reduction failed")

    for u in range(1,isqrt(m)+1,2):
        N=m-u*u
        if N<=0:
            break
        pair=canonical_two_squares(N, require_odd_pair=True)
        if pair is not None:
            s,t=pair
            x=(r,s,t,u)   # r + s e1 + t e3 + u e7
            if qnorm(x)!=p:
                raise AssertionError("norm failure")
            return x

    raise RuntimeError(f"No inert H137 anchor found for p={p}")

def canonical_state(p):
    if not is_prime(p) or p==2:
        raise ValueError("odd prime required")
    return gauss_anchor(p) if p%4==1 else inert_anchor_mod8(p)

def alpha(x):
    # alpha_p = -e1*x
    return tuple(-v for v in qmul(E1,x))

def bridge(p,q):
    xp=canonical_state(p)
    xq=canonical_state(q)

    ap=alpha(xp)
    aq=alpha(xq)

    pivot_num=qmul(xp,qconj(ap))
    if any(v%p for v in pivot_num):
        raise ArithmeticError("pivot divisibility failure")

    pivot=tuple(v//p for v in pivot_num)
    if pivot!=E1:
        raise AssertionError(f"pivot failure: {pivot}")

    numerator=qmul(pivot_num,aq)
    if any(v%p for v in numerator):
        raise ArithmeticError("bridge integrality failure")

    y=tuple(v//p for v in numerator)
    if y!=xq:
        raise AssertionError("target mismatch")
    if qnorm(y)!=q:
        raise AssertionError("target norm mismatch")

    return {
        "p":p,
        "q":q,
        "x_p":xp,
        "x_q":xq,
        "pivot":pivot,
        "output":y,
    }

if __name__=="__main__":
    tests=[
        (17,19),
        (19,23),
        (23,29),
        (1009,1013),
        (1013,1019),
        (1019,1021),
        (100003,100019),
        (1000003,1000033),
    ]
    for p,q in tests:
        r=bridge(p,q)
        print(p,"->",q,"PASS",r)
