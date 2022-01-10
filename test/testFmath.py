
# -*- coding: utf-8 -*-

# Test base classes.

__all__ = ('Tests',)
__version__ = '22.01.06'

from base import coverage, TestsBase

from pygeodesy import cbrt, cbrt2, euclid_, Ellipsoids, facos1, fasin1, \
                      fatan, fatan2, fhorner, fmath, fpolynomial, fpowers, \
                      Fsum, fsum, fsum_, hypot, hypot_, hypot2_, norm_, \
                      signOf, sqrt3  # hypot3  # DEPRECATED

from math import acos, asin, atan, atan2, ceil, floor, sqrt
from random import random, gauss, shuffle


class Tests(TestsBase):

    def testFmath(self):  # MCCABE 14

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
            self.test(E.name, Ah, E.A, prec=10, known=abs(Ah - E.A) < 1e-5)  # b_None, f_None on iPhone
            Ah = E.a / (1 + E.n) * (fhorner(E.n**2, 16384, 4096, 256, 64, 25) / 16384)
            self.test(E.name, Ah, E.A, prec=10, known=abs(Ah - E.A) < 1e-5)  # b_None, f_None on iPhone

            Ah = E.a / (1 + E.n) * fhorner(E.n**2, 1., 1./4, 1./64, 1./256, 25./16384, 49./65536)
            self.test(E.name, Ah, E.A, prec=10, known=abs(Ah - E.A) < 1e-9)
            Ah = E.a / (1 + E.n) * (fhorner(E.n**2, 65536, 16384, 1024, 256, 100, 49) / 65536)
            self.test(E.name, Ah, E.A, prec=10, known=abs(Ah - E.A) < 1e-9)

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
            self.test('FSum.', a, 'fmath.Fsum[2049]')

        try:
            self.test('_2sum', fmath._2sum(1e308, 1e803), OverflowError.__name__)
        except OverflowError as x:
            self.test('_2sum', repr(x), repr(x))

        # cffk <https://Bugs.Python.org/issue43088> and
        # <https://Bugs.Python.org/file49783/hypot_bug.py>
        x  = 0.6102683302836215
        y1 = 0.7906090004346522
        y2 = y1 + 1e-16
        z1 = hypot(x, y1)
        z2 = hypot(x, y2)
        self.test('hypot', signOf(y2 - y1) == signOf(z2 - z1), True, known=True)  # (3, 7) < sys.version_info[:2] < (3, 10))

        h = hypot_(1.0, 0.0050, 0.0000000000010)
        self.test('hypot_ ', h, '1.00001250', prec=8)
        e = euclid_(1.0, 0.0050, 0.0000000000010)
        self.test('euclid_', e, h, prec=8, known=abs(e - h) < h * 0.01)
        t = hypot2_(1.0, 0.0050, 0.0000000000010)
        self.test('hypot2_', t, '1.00002500', prec=8)
        s = hypot_(*norm_(1.0, 0.0050, 0.0000000000010))  # PYCHOK norm_
        self.test('norm_  ', s, 1.0, prec=8)

        h = hypot_(3000, 2000, 100)  # note, all C{int}
        self.test('hypot_ ', h, '3606.937759', prec=6)
        e = euclid_(3000, 2000, 100)
        self.test('euclid_', e, h * 1.07, prec=6, known=abs(e - h) < h * 0.07)
        t = hypot2_(3000, 2000, 100)  # note, all C{int}
        s = fsum_(3000**2, 2000**2, 100**2)
        self.test('hypot2_', t, s, prec=1)
        s = hypot_(*norm_(3000, 2000, 100))  # PYCHOK norm_
        self.test('norm_  ', s, 1.0, prec=1)

        h = hypot_(40000, 3000, 200, 10.0)
        s = fsum_(40000**2, 3000**2, 200**2, 100)
        self.test('hypot_ ', h, sqrt(s), prec=3)
        t = hypot2_(40000, 3000, 200, 10.0)
        self.test('hypot2_', t, s, prec=1)
        e = euclid_(40000, 3000, 200, 10.0)
        self.test('euclid_', e, h * 1.03, prec=3, known=abs(e - h) < h * 0.03)

        self.test('cbrt',  cbrt(27),   '3.00', prec=2)
        self.test('cbrt',  cbrt(-27), '-3.00', prec=2)
        self.test('cbrt2', cbrt2(27),  '9.00', prec=2)
        self.test('cbrt2', cbrt2(-27), '9.00', prec=2)
        self.test('sqrt3', sqrt3(9),  '27.00', prec=2)

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
        self.test('F * 2', d.fsum(), '2.0', prec=-8, known=True)
        s = d / 2  # /= chockes PyChecker
        self.test('F / 2', float(s), float(t), prec=-8, known=True)
        self.test('F / F', s / s == Fsum(1.0), True, known=True)  # PYCHOK "s / s is always 1 or ZeroDivisionError"
        self.test('F / F', float(s / t), '1.0', known=True)
        self.test('F / F', float(d / t), '2.0', known=True)
        self.test('abs  ', abs(t), t, known=True)
        self.test('int  ', int(t), 2, known=True)
        if coverage:
            self.test('ceil ', ceil(t.fcopy().fadd_(1e-15)), '3.0', known=True)
            self.test('floor', floor(t), '2.0', known=True)
            self.test('divmod ', divmod(d, 2), (2, d), known=True)
            self.test('is_int', t.is_integer(), True, known=True)
            x = Fsum(1, 1e-101, -1, -1e-102)
            self.test('is_int', float(x), '9e-102', known=True)
            self.test('is_int', x.is_integer(), False, known=True)
            self.test('mod ', d % 2, t, known=True)
            self.test('pos ', +t, t, known=True)  # PYCHOK "Unary positive (+) usually has no effect"
            try:
                self.test('pow(F)', pow(x, 2), NotImplementedError.__name__)
            except Exception as x:
                self.test('pow(F)', repr(x), NotImplementedError.__name__, known=True)

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

        def _percent(a, x, f):
            return max(a, abs((f - x) * 100 / x) if x else 0)

        c = s = a = t = 0
        for f in range(100):
            f =  float(f)
            a = _percent(a, atan(f), fatan(f))
            t = _percent(t, atan2(f, 2.0), fatan2(f, 2.0))
            f *= 0.01
            s = _percent(s, asin(f), fasin1(f))
            c = _percent(c, acos(f), facos1(f))
        self.test('facos1', c, '0.005%', fmt='%.3f%%', known=True)
        self.test('fasin1', s, '0.439%', fmt='%.3f%%', known=True)
        self.test('fatan ', a, '0.508%', fmt='%.3f%%', known=True)
        self.test('fatan2', t, '1.214%', fmt='%.3f%%', known=True)


if __name__ == '__main__':

    t = Tests(__file__, __version__, fmath)
    t.testFmath()
    t.results()
    t.exit()
