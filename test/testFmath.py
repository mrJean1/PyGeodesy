
# -*- coding: utf-8 -*-

# Test L{fmath} module.

__all__ = ('Tests',)
__version__ = '23.08.30'

from bases import endswith, isWindows, randoms, startswith, TestsBase

from pygeodesy import EPS, Fcbrt, Fhypot, INF, Fn_rt, Fpowers, Fsqrt, Fsum, \
                      bqrt, cbrt, cbrt2, euclid_, Ellipsoids, facos1, fasin1, \
                      fatan, fatan1, fatan2, fhorner, fmath, fpolynomial, \
                      fpowers, fsum_, hypot, hypot_, hypot2_, norm_, signOf, \
                      sqrt3, sqrt_a, zcrt, zqrt

from math import acos, asin, atan, atan2, sqrt
from sys import version_info as _version_info
# hypot issue in Python 3.8 and 3.9 plus 2.7 32-bit on Windows
_hypot_known = (3, 8) <= _version_info[:2] <= (3, 9) or (isWindows and
               (2, 7) == _version_info[:2])


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
            self.test('fhornerB', h, p, nt=1)

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
        self.test('fpowers', p[3], 2**9, nt=1)

        for _, E in Ellipsoids.items(all=True, asorted=True):
            a_n = E.a / (1 + E.n)
            n2  = E.n**2
            Ah = a_n * fhorner(n2, 1, 1./4, 1./64, 1./256, 25./16384)
            self.test(E.name, Ah, E.A, prec=9, known=abs(Ah - E.A) < 1e-5)  # b_None, f_None on iPhone
            Ah = a_n * (fhorner(n2, 16384, 4096, 256, 64, 25) / 16384)
            self.test(E.name, Ah, E.A, prec=9, known=abs(Ah - E.A) < 1e-5)  # b_None, f_None on iPhone

            Ah = a_n * fhorner(n2, 1, 1./4, 1./64, 1./256, 25./16384, 49./65536)
            self.test(E.name, Ah, E.A, prec=9, known=abs(Ah - E.A) < 1e-8)
            Ah = a_n * (fhorner(n2, 65536, 16384, 1024, 256, 100, 49) / 65536)
            self.test(E.name, Ah, E.A, prec=9, known=abs(Ah - E.A) < 1e-8, nt=1)

        # cffk <https://Bugs.Python.org/issue43088> and <https://Bugs.Python.org/file49783/hypot_bug.py>
        x  = 0.6102683302836215
        y1 = 0.7906090004346522
        y2 = y1 + 1e-16
        h1 = hypot(x, y1)
        h2 = hypot(x, y2)
        self.test('hypot', signOf(y2 - y1), signOf(h2 - h1), known=_hypot_known)
        self.test('sqrt_a', sqrt_a(h1, y1), x, prec=13)
        self.test('sqrt_a', sqrt_a(h2, y2), x, prec=13)

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

        self.test('bqrt',  bqrt(16),   '2.00', prec=2)
        self.test('cbrt',  cbrt(27),   '3.00', prec=2)
        self.test('cbrt',  cbrt(-27), '-3.00', prec=2)
        self.test('cbrt2', cbrt2(27),  '9.00', prec=2)
        self.test('cbrt2', cbrt2(-27), '9.00', prec=2)
        self.test('sqrt3', sqrt3(9),  '27.00', prec=2)
        self.test('zcrt',  zcrt(64),   '2.00', prec=2)
        self.test('zqrt',  zqrt(256),  '2.00', prec=2)

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

        f = Fsum(3)
        n = Fhypot.__name__
        self.test(n, Fhypot(f, 4.0), '(5.0, 0)', known=endswith, nl=1)
        f = -f
        self.test(n, Fhypot(f, -4, 8), '(9.43398, 0)', known=endswith)
        self.test(n, Fhypot(f, Fsum(-4, 8)), '(5.0, 0)', known=endswith)
        self.test(n, Fhypot(f, -4, 8, power=-1), '(-2.18182, 0)', known=endswith)
        self.test(n, Fhypot(f, Fsum(-4, 8), power=-1), '(-12, 0)', known=endswith)
        f = Fsum(-1)
        self.test(n, Fhypot(f, -1), '(1.41421, 0)', known=endswith)
        self.test(n, Fhypot(f, -1, power=-1), '(-0.5, 0)', known=endswith)
        try:  # pow(INF, -1) == 0.0
            self.test(n, Fhypot(f, INF, f, f, power=-1), ValueError.__name__, known=True)
        except Exception as x:
            self.test(n, str(x), ' not finite', known=endswith)
        try:
            self.test(n, Fhypot(f, -1, power=0), ZeroDivisionError.__name__)
        except Exception as x:
            self.test(n, str(x), ' float division by zero', known=endswith)

        f = Fsum(3)
        n = Fsqrt.__name__
        r = Fn_rt.__name__
        self.test(n, Fsqrt(f, 6), '(3.0, 0)', known=endswith, nl=1)
        self.test(n, Fcbrt(20, 4, f), '(3.0, 0)', known=endswith)
        f = -f
        self.test(n, Fsqrt(f, -4, 9),       '(1.41421, 0)', known=endswith)
        self.test(n, Fsqrt(f, Fsum(-4, 9)), '(1.41421, 0)', known=endswith)
        self.test(r, Fn_rt(-1, f, -4, 9), '(0.5, 0)', known=endswith)
        self.test(r, Fn_rt(-1, f, Fsum(-4, 9)), '(0.5, 0)', known=endswith)
        f = Fsum(-1)
        try:
            self.test(n, Fsqrt(f, -1), ValueError.__name__)
        except Exception as x:
            self.test(n, str(x), "fmath.Fsqrt(<fsums.Fsum[1] (-1, 0) ", known=startswith)
        try:
            self.test(r, Fn_rt(0, f, -1), ValueError.__name__)
        except Exception as x:
            self.test(n, str(x), ' float division by zero', known=endswith)
        self.test(r, Fn_rt(-1, f, -3), '(-0.25, 0)', known=endswith)

        i = 0
        while i < 9:
            t = randoms(i)
            F = Fsum(*t)
            s = fsum_(F, *t)
            if s > 0:
                i += 1
                n  = str(i)
                h  = hypot_(F, *t)
                c  = fsum_(F**3, *(_**3 for _ in t))
                self.test(Fhypot.__name__  + n, float(Fhypot( F, *t, RESIDUAL=EPS)), h,       fmt='%.9g', nl=1)
                self.test(Fpowers.__name__ + n, float(Fpowers(3, F, *t)),            c,       fmt='%.9g')
                self.test(Fsqrt.__name__   + n, float(Fsqrt(  F, *t, RESIDUAL=EPS)), sqrt(s), fmt='%.9g')
                self.test(Fcbrt.__name__   + n, float(Fcbrt(  F, *t, RESIDUAL=EPS)), cbrt(s), fmt='%.9g')


if __name__ == '__main__':

    t = Tests(__file__, __version__, fmath)
    t.testFmath()
    t.results()
    t.exit()
