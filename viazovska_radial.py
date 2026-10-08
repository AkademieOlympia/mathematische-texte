from dataclasses import dataclass
import mpmath as mp

mp.mp.dps = 50

@dataclass(frozen=True)
class FourierInterpolationKernel:
    n: int

    def a(self, r: mp.mpf) -> mp.mpf:
        """Basisfunktion a_n(r) zur Isolation von f(sqrt(2n))"""
        raise NotImplementedError("G3b.2b.3: Kernel a_n noch nicht aus CKMRV transkribiert.")

    def da(self, r: mp.mpf) -> mp.mpf:
        """Radiale Ableitung a'_n(r)"""
        raise NotImplementedError("G3b.2b.3: Kernel da_n noch nicht aus CKMRV transkribiert.")

    def b(self, r: mp.mpf) -> mp.mpf:
        """Basisfunktion b_n(r) zur Isolation von f'(sqrt(2n))"""
        raise NotImplementedError("G3b.2b.3: Kernel b_n noch nicht aus CKMRV transkribiert.")

    def db(self, r: mp.mpf) -> mp.mpf:
        """Radiale Ableitung b'_n(r)"""
        raise NotImplementedError("G3b.2b.3: Kernel db_n noch nicht aus CKMRV transkribiert.")

    def Fa(self, r: mp.mpf) -> mp.mpf:
        """Fouriertransformierte Korrespondenz zu a_n"""
        raise NotImplementedError("G3b.2b.3: Fourier-Kern Fa_n noch nicht implementiert.")

    def dFa(self, r: mp.mpf) -> mp.mpf:
        """Radiale Ableitung der Fouriertransformierten zu a_n"""
        raise NotImplementedError("G3b.2b.3: Fourier-Kern dFa_n noch nicht implementiert.")

    def Fb(self, r: mp.mpf) -> mp.mpf:
        """Fouriertransformierte Korrespondenz zu b_n"""
        raise NotImplementedError("G3b.2b.3: Fourier-Kern Fb_n noch nicht implementiert.")

    def dFb(self, r: mp.mpf) -> mp.mpf:
        """Radiale Ableitung der Fouriertransformierten zu b_n"""
        raise NotImplementedError("G3b.2b.3: Fourier-Kern dFb_n noch nicht implementiert.")