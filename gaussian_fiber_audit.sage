from sage.all import I, Primes, ZZ, factor


class GaussianFiberAudit:
    def __init__(self):
        self.R = ZZ[I]
        self.primes = Primes()

    def projection(self, element):
        """Scalar projection pi(z) = N(z)."""
        return ZZ(element.norm())

    def compute_fiber(self, n, search_bound=5):
        """Compute all Gaussian integers of norm n in a finite box."""
        fiber = []
        for real_part in range(-search_bound, search_bound + 1):
            for imag_part in range(-search_bound, search_bound + 1):
                z = self.R(real_part + imag_part * I)
                if self.projection(z) == n:
                    fiber.append(z)
        return fiber

    def are_associated(self, x, y):
        """Test association using the four units of ZZ[i]."""
        units = [self.R(1), self.R(-1), self.R(I), self.R(-I)]
        return any(x == unit * y for unit in units)

    def verify(self, a, b, x, y):
        if a not in self.primes or b not in self.primes:
            raise ValueError("a and b must be natural primes")

        scalar = ZZ(a * b)
        fiber = self.compute_fiber(scalar)
        result = {
            "scalar": scalar,
            "scalar_factorization": factor(scalar),
            "pi_x": self.projection(x),
            "pi_y": self.projection(y),
            "same_fiber": self.projection(x) == scalar
                          and self.projection(y) == scalar,
            "identical": x == y,
            "associated": self.are_associated(x, y),
            "x_factorization": x.factor(),
            "y_factorization": y.factor(),
            "fiber": fiber,
        }
        return result


if __name__ == "__main__":
    audit = GaussianFiberAudit()

    # Scalar target: 25 = 5^2.
    a = 5
    b = 5

    # Two distinct representatives in the norm fiber over 25.
    x = audit.R(5)
    y = audit.R(3 + 4 * I)

    result = audit.verify(a, b, x, y)

    assert result["same_fiber"]
    assert not result["identical"]
    assert not result["associated"]
    assert len(result["fiber"]) == 12

    for key, value in result.items():
        print(f"{key}: {value}")
