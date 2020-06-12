
# -*- coding: utf-8 -*-

# Test base classes.

__all__ = ('Tests',)
__version__ = '20.06.02'

from base import coverage, TestsBase

from pygeodesy import cbrt, cbrt2, Ellipsoids, fhorner, fmath, \
                      fpolynomial, fpowers, Fsum, fsum, fsum_, \
                      hypot3, hypot_, sqrt3

from math import sqrt
from random import random, gauss, shuffle


class Tests(TestsBase):

    def testFmath(self):

        # see function _p2 in ellipsoidalVincenty.py
        for i in range(-16, 17):
            x = pow(1.0, i)

            p = fpolynomial(x, 16384, 4096, -768, 320, -175) / 16384.0
            a = x / 16384.0 * (4096 + x * (-768 + x * (320 - 175 * x))) + 1
            self.test('fpolynomialA', p, a)
            h = fhorner(x, 16384, 4096, -768, 320, -175) / 16384.0
            self.test('fhornerA', h, p)

            p = fpolynomial(x, 0, 256, -128, 74, -47) / 1024.0
            b = x / 1024.0 * (256 + x * (-128 + x * (74 - 47 * x)))
            self.test('fpolynomialB', p, b)
            h = fhorner(x, 0, 256, -128, 74, -47) / 1024.0
            self.test('fhornerB', h, p)

        # U{Neumaier<https://WikiPedia.org/wiki/Kahan_summation_algorithm>}
        t = 1, 1e101, 1, -1e101
        for _ in range(10):
            s = float(len(t) / 2)  # number of ones
            self.test('sum', sum(t), s, known=True)
            self.test('fsum', fsum(t), s)
            self.test('Fsum', Fsum().fsum(t), s)
            t += t

        # <https://code.ActiveState.com/recipes/393090>
        t = 1.0, 0.0050, 0.0000000000010
        s = 1.0050000000010
        self.test('sum', sum(t), s, known=True)
        self.test('fsum', fsum(t), s)
        self.test('Fsum', Fsum().fsum(t), s)

        # <https://GitHub.com/python/cpython/blob/master/Modules/mathmodule.c>
        t = 1e-16, 1, 1e16
        s = '1.0000000000000002e+16'  # half-even rounded
        self.test('fsum', fsum(t), s, fmt='%.16e')
        self.test('Fsum', Fsum().fsum(t), s, fmt='%.16e')

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
        self.test('fmul', p.fsum(), f.fsum(), fmt='%0.8f')

        p = fpowers(2, 10)  # PYCHOK false!
        self.test('fpowers', len(p), 10)
        self.test('fpowers', p[0], 2)
        self.test('fpowers', p[9], 2**10)
        p = fpowers(2, 10, 4)  # PYCHOK false!
        self.test('fpowers', len(p), 4)
        self.test('fpowers', p[0], 2**4)
        self.test('fpowers', p[3], 2**10)
        p = fpowers(2, 10, 3)  # PYCHOK false!
        self.test('fpowers', len(p), 4)
        self.test('fpowers', p[0], 2**3)
        self.test('fpowers', p[3], 2**9)

        for _, E in sorted(Ellipsoids.items()):
            Ah = E.a / (1 + E.n) * fhorner(E.n**2, 1., 1./4, 1./64, 1./256, 25./16384)
            self.test(E.name, Ah, E.A, fmt='%.8f')
            Ah = E.a / (1 + E.n) * (fhorner(E.n**2, 16384, 4096, 256, 64, 25) / 16384)
            self.test(E.name, Ah, E.A, fmt='%.8f')

            Ah = E.a / (1 + E.n) * fhorner(E.n**2, 1., 1./4, 1./64, 1./256, 25./16384, 49./65536)
            self.test(E.name, Ah, E.A, fmt='%.10f')
            Ah = E.a / (1 + E.n) * (fhorner(E.n**2, 65536, 16384, 1024, 256, 100, 49) / 65536)
            self.test(E.name, Ah, E.A, fmt='%.10f')

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
            self.testCopy(a, '_fsum2_', '_n', '_ps')

        if coverage:  # for test coverage
            c = a - b
            self.test('FSum0', c.fsum(), 0.0)
            c -= 0
            self.test('FSum0', c.fsum(), 0.0)
            c -= c
            self.test('FSum0', c.fsum(), 0.0)
            c *= Fsum(1.0)
            self.test('FSum0', c.fsum(), 0.0)
            a.fsub_(*a._ps)
            self.test('FSum0', a.fsum(), 0.0)
            self.test('Fsum#', len(a), 2049)
            self.test('Fsum#', len(a._ps), 1)
            self.test('FSum.', a, 'fmath.Fsum()')

        try:
            self.test('_2sum', fmath._2sum(1e308, 1e803), OverflowError.__name__)
        except OverflowError as x:
            self.test('_2sum', repr(x), repr(x))

        h = hypot_(1.0, 0.0050, 0.0000000000010)
        self.test('hypot_', h, '1.0000124999219', fmt='%.13f')
        s = hypot3(1.0, 0.0050, 0.0000000000010)  # DEPRECATED
        self.test('hypot3', s, h, fmt='%.13f')

        h = hypot_(3000, 2000, 100)  # note, all C{int}
        self.test('hypot_', h, '3606.937759', fmt='%.6f')
        s = hypot3(3000, 2000, 100)  # DEPRECATED
        self.test('hypot3', s, h, fmt='%.6f')

        h = hypot_(40000, 3000, 200, 10.0)
        s = sqrt(fsum_(40000**2, 3000**2, 200**2, 100))
        self.test('hypot_', h, s, fmt='%.3f')

        self.test('cbrt', cbrt(27),   '3.0')
        self.test('cbrt', cbrt(-27), '-3.0')
        self.test('cbrt2', cbrt2(27),  '9.0')
        self.test('cbrt2', cbrt2(-27), '9.0')
        self.test('sqrt3', sqrt3(9),  '27.0')

        # Knuth/Kulisch, TAOCP, vol 2, p 245, sec 4.2.2, item 31, see also .testKarney.py
        # <https://SeriousComputerist.Atariverse.com/media/pdf/book/
        #          Art%20of%20Computer%20Programming%20-%20Volume%202%20(Seminumerical%20Algorithms).pdf>
        x = 408855776
        y = 708158977
        self.test('ints', 2*y**2 +  9*x**4 - y**4, 1)
        self.test('ints', 2*y**2 + (3*x**2 - y**2) * (3*x**2 + y**2), 1)
        t = 2*float(y)**2, 9*float(x)**4, -(float(y)**4)
        self.test('fsum ', fsum(t),          '1.0', fmt='%.12e', known=True)  # -3.589050987400773e+19
        self.test('fsum_', fsum_(*t),        '1.0', fmt='%.12e', known=True)
        self.test('Fsum ', Fsum().fsum_(*t), '1.0', fmt='%.12e', known=True)
        self.test('sum  ',  sum(t),          '1.0', fmt='%.12e', known=True)


if __name__ == '__main__':

    t = Tests(__file__, __version__, fmath)
    t.testFmath()
    t.results()
    t.exit()
