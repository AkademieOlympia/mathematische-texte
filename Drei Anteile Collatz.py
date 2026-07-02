def reaches_three_smooth(m, max_iter=10000000000):
    """Test if Collatz orbit from m hits a 3-smooth number."""
    n = m
    for k in range(max_iter):
        # Check if n = 2^a * 3^b
        temp = n
        while temp % 2 == 0:
            temp //= 2
        while temp % 3 == 0:
            temp //= 3
        if temp == 1:
            return True, k  # Hit 3-smooth number at step k
        
        # Collatz step
        n = 3*n + 1 if n % 2 == 1 else n // 2
    
    return False, max_iter