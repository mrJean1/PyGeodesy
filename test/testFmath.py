
# -*- coding: utf-8 -*-

# Test base classes.

__all__ = ('Tests',)
__version__ = '22.02.27'

from base import TestsBase

from pygeodesy import cbrt, cbrt2, euclid_, Ellipsoids, facos1, fasin1, \
                      fatan, fatan1, fatan2, fhorner, fmath, fpolynomial, fpowers, \
                      fsum_, hypot, hypot_, hypot2_, norm_, signOf, sqrt3, sqrt_a

from math import acos, asin, atan, atan2, sqrt


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

        # cffk <https://Bugs.Python.org/issue43088> and <https://Bugs.Python.org/file49783/hypot_bug.py>
        x  = 0.6102683302836215
        y1 = 0.7906090004346522
        y2 = y1 + 1e-16
        z1 = hypot(x, y1)
        z2 = hypot(x, y2)
        self.test('hypot', signOf(y2 - y1) == signOf(z2 - z1), True, known=True)  # (3, 7) < sys.version_info[:2] < (3, 10))
        self.test('sqrt_a', sqrt_a(z1, y1), x, prec=9)

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

        def _percent(a, x, f):
            return max(a, abs((f - x) * 100 / x) if x else 0)

        c = s = t = t1 = t2 = 0
        for f in range(100):
            f  =  float(f)
            t  = _percent(t,  atan(f), fatan(f))
            t2 = _percent(t2, atan2(f, 2.0), fatan2(f, 2.0))
            f *= 0.01
            c  = _percent(c,  acos(f), facos1(f))
            s  = _percent(s,  asin(f), fasin1(f))
            t1 = _percent(t1, atan(f), fatan1(f))
        self.test('facos1', c,  '0.005%', fmt='%.3f%%', known=True)
        self.test('fasin1', s,  '0.439%', fmt='%.3f%%', known=True)
        self.test('fatan ', t,  '0.134%', fmt='%.3f%%', known=True)
        self.test('fatan1', t1, '2.834%', fmt='%.3f%%', known=True)
        self.test('fatan2', t2, '0.321%', fmt='%.3f%%', known=True)


if __name__ == '__main__':

    t = Tests(__file__, __version__, fmath)
    t.testFmath()
    t.results()
    t.exit()
