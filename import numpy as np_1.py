import numpy as np

zeros = np.load("zeros6.npy")          # ~sehr schnell
# zeros = np.load("zeros6.npz")["zeros"]  # Alternative

def nearest_zero(t):
    i = np.searchsorted(zeros, t)
    if i == 0: 
        return zeros[0]
    if i == len(zeros): 
        return zeros[-1]
    left, right = zeros[i-1], zeros[i]
    return left if (t-left) <= (right-t) else right

def nearest_dist(t):
    z = nearest_zero(t)
    return abs(t - z), z