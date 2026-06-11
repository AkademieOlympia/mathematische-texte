import numpy as np, math

def angle(u, v):
    # Winkel zwischen u,v
    cuv = np.dot(u, v) / (np.linalg.norm(u)*np.linalg.norm(v))
    cuv = max(-1.0, min(1.0, cuv))
    return math.acos(cuv)

def tet_volume(A, B, C, D):
    return abs(np.linalg.det(np.vstack([B-A, C-A, D-A]).T)) / 6.0

def tet_deficits(A, B, C, D):
    # Defizite an den vier Ecken des Tetraeders
    P = [A,B,C,D]
    deficits = []
    for i in range(4):
        a = P[i]
        others = [P[j] for j in range(4) if j != i]  # b,c,d
        # drei Face-Winkel an a:
        # (a, o0, o1), (a, o0, o2), (a, o1, o2)
        v0, v1, v2 = [o - a for o in others]
        th01 = angle(v0, v1)
        th02 = angle(v0, v2)
        th12 = angle(v1, v2)
        delta = 2*math.pi - (th01 + th02 + th12)
        deficits.append(delta)
    return np.array(deficits)

def D_curvature(A,B,C,D, mode="sq"):
    d = tet_deficits(A,B,C,D)
    if mode == "sq":
        return float(np.sum(d*d))
    elif mode == "center_pi":
        return float(np.sum((d-math.pi)**2))
    else:
        raise ValueError("mode unknown")

def P_press(A,B,C,D, p, Vref, mode="lin"):
    V = tet_volume(A,B,C,D)
    if mode == "lin":
        return p*(Vref - V)
    elif mode == "sq":
        return p*(Vref - V)**2
    else:
        raise ValueError("mode unknown")

def E_geom(A,B,C,D, p, Vref, alpha=1.0, beta=1.0,
           Dmode="sq", Pmode="lin"):
    return alpha*D_curvature(A,B,C,D, Dmode) + beta*P_press(A,B,C,D, p, Vref, Pmode)
