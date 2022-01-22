
# -*- coding: utf-8 -*-

# Test base classes.

__all__ = ('Tests',)
__version__ = '22.01.20'

from base import startswith, TestsBase

from pygeodesy import Fsum, fsum, fsum_, fsums

from math import ceil, floor
from random import random, gauss, shuffle


class Tests(TestsBase):

    def testFsums(self):  # MCCABE 14

        # U{Neumaier<https://WikiPedia.org/wiki/Kahan_summation_algorithm>}
        t = 1, 1e101, 1, -1e101
        for _ in range(10):
            s = float(len(t) / 2)  # number of ones
            self.test('sum', sum(t), s, known=True)
            self.test('fsum', fsum(t), s)
            self.test('Fsum', Fsum().fsum(t), s)
            t += t

        # <https://GitHub.com/ActiveState/code/blob/master/recipes/Python/
        #        393090_Binary_floating_point_summatiaccurate_full/recipe-393090.py>
        t = 1.0, 0.0050, 0.0000000000010
        s = 1.0050000000010
        self.test('sum', sum(t), s, known=True)
        self.test('fsum', fsum(t), s)
        self.test('Fsum', Fsum().fsum(t), s)

        # <https://GitHub.com/python/cpython/blob/master/Modules/mathmodule.c>
        t = 1e-16, 1, 1e16
        s = '1.0000000000000002e+16'  # half-even rounded
        self.test('fsum', fsum(t), s, prec=-16)
        self.test('Fsum', Fsum().fsum(t), s, prec=-16)

        # <https://GitHub.com/ActiveState/code/blob/master/recipes/Python/
        #        393090_Binary_floating_point_summatiaccurate_full/recipe-393090.py>
        for _ in range(100):
            t = [7, 1e100, -9e-20, -7, -1e100, 8e-20] * 10
            s = 0
            for _ in range(20):
                v = gauss(0, random())**7 - s
                s += v
                t.append(v)
            shuffle(t)
            s = fsum(t)
            self.test('sum', sum(t), s, known=True)
            self.test('fsum', s, s)
            n = len(t) // 2
            f = Fsum()
            f.fsum(t[:n])  # test ps
            self.test('Fsum', f.fsum(t[n:]), s)
            f = Fsum()
            f.fsum(t[n:])  # test ps
            self.test('Fsum', f.fsum(t[:n]), s)

        p = f * (f * 1e10)  # coverage Fsum.__imul__
        f *= f * 1e10
        self.test('fmul', p.fsum(), f.fsum(), prec=8)

        t = 1, 1e101, 1, -1e101
        for _ in range(10):
            a = Fsum(*t)
            b = a.fcopy()
            c = a + b
            self.test('FSum+', c.fsum(), a.fsum() + b.fsum())
            c -= a
            self.test('FSum-', c.fsum(), b.fsum())
            c -= b
            self.test('FSum-', c.fsum(), 0.0)
            b = a * 2
            a += a
            self.test('FSum*', a.fsum(), b.fsum())
            t += t
            self.testCopy(a, '_fsum2_', '_n', '_ps', deep=True)  # _n

        c = a - b
        self.test('FSum0', c.fsum(), 0.0)
        c -= 0
        self.test('FSum0', c.fsum(), 0.0)
        c -= c
        self.test('FSum0', c.fsum(), 0.0)
        c *= Fsum(1.0)
        self.test('FSum0', c.fsum(), 0.0)
        t = a + a.fcopy().fmul(-1)
        self.test('FSum0', t.fsum(), 0.0)
        a.fsub_(*a._ps)
        self.test('FSum0', a.fsum(), 0.0)
        self.test('Fsum#', len(a), 2049)
        self.test('Fsum#', len(a._ps), 1)
        self.test('FSum.', a, 'fsums.Fsum[2049] (0.0, 0.0)')

        self.test('FsumI', c.imag, 0.0)
        self.test('FsumR', c.real, float(c))

        try:
            self.test('_2sum', fsums._2sum(1e308, 1e803), OverflowError.__name__)
        except OverflowError as x:
            self.test('_2sum', repr(x), repr(x))

        # Knuth/Kulisch, TAOCP, vol 2, p 245, sec 4.2.2, item 31, see also .testKarney.py
        # <https://SeriousComputerist.Atariverse.com/media/pdf/book/
        #  Art%20of%20Computer%20Programming%20-%20Volume%202%20(Seminumerical%20Algorithms).pdf>
        x = 408855776
        y = 708158977
        self.test('ints', 2*y**2 +  9*x**4 - y**4, 1)
        self.test('ints', 2*y**2 + (3*x**2 - y**2) * (3*x**2 + y**2), 1)
        t = 2*float(y)**2, 9*float(x)**4, -(float(y)**4)
        self.test('fsum ', fsum(t),          '1.0', prec=-8, known=True)  # -3.589050987400773e+19
        self.test('fsum_', fsum_(*t),        '1.0', prec=-8, known=True)
        self.test('Fsum ', Fsum().fsum_(*t), '1.0', prec=-8, known=True)
        self.test('sum  ',  sum(t),          '1.0', prec=-8, known=True)

        t = Fsum(1, 1e101, 1, -1e101)
        d = t * 2
        self.test('F * 2', float(d), 4.0, prec=4, known=True)
        s = d / 2  # /= chockes PyChecker
        self.test('F / 2', float(s), float(t), prec=4, known=True)
        self.test('F / F', s / s == Fsum(1.0), True, known=True)  # PYCHOK "s / s is always 1 or ZeroDivisionError"
        self.test('F / F', float(s / t), '1.0', known=True)
        self.test('F / F', float(d / t), '2.0', known=True)
        self.test('abs  ', abs(t), t, known=True)
        self.test('int  ', int(t), 2, known=True)

        self.test('ceil ', ceil(t.fcopy().fadd_(1e-15)), '3', known=True)
        self.test('floor', floor(t), '2', known=True)
        self.test('divmod ', divmod(d, 2), "(2.0, <fsums.Fsum '__divmod__'[2]", known=startswith)
        self.test('mod ', d % 2, "fsums.Fsum '__divmod__'[2] (0.0, 0.0)")
        self.test('pos ', +t, t, known=True)  # PYCHOK "Unary positive (+) usually has no effect"
        self.test('is_int', t.is_integer(), True, known=True)
        x = Fsum(1, 1e-101, -1, -1e-102)
        self.test('float', float(x), '9e-102', known=True)
        self.test('is_int', x.is_integer(), False, known=True)
        try:
            self.test('pow(F)', pow(x, 2.1), ValueError.__name__)
        except Exception as X:
            self.test('pow(F)', repr(X), ValueError.__name__, known=startswith)
        try:
            self.test('pow(F, m)', pow(x, 2.1, 2), NotImplementedError.__name__, known=True)
        except Exception as X:
            self.test('pow(F, m)', repr(X), NotImplementedError.__name__, known=True)
        try:
            self.test('Z**-2', Fsum(0.0)**-2, ValueError.__name__)
        except Exception as X:
            self.test('Z**-2', repr(X), ValueError.__name__, known=startswith)
        self.test('pow(0)', x**0, '1.000', prec=3)
        self.test('pow(1)', x**1, '0.000', prec=3)
        self.test('pow(2)', x**2, '0.000', prec=3)
        x = Fsum(2, 3)
        self.test('pow(-1)',   x**-1,   '0.200', prec=3)
        self.test('pow(-2)',   x**-2,   '0.040', prec=3)
        self.test('pow(-2.5)', x**-2.5, '0.018', prec=3)

        self.test('eq F', t == s,  True)
        self.test('ge F', t >= s,  True)
        self.test('gt F', t > s,  False)
        self.test('le F', t <= s,  True)
        self.test('lt F', t <  s, False)
        self.test('ne F', t != s, False)
        self.test('if F', bool(s), True)

        self.test('gt 0', t > 0, True)
        self.test('lt 0', t < 0, False)
        self.test('eq 0', t == 0, False)
        s = -t  # t is positive
        self.test('lt 0', s < 0, True)
        self.test('gt 0', t > 0, True)


if __name__ == '__main__':

    t = Tests(__file__, __version__, fsums)
    t.testFsums()
    t.results()
    t.exit()
