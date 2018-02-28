
# -*- coding: utf-8 -*-

# Test base classes.

__all__ = ('Tests',)
__version__ = '18.02.09'

from base import TestsBase

from pygeodesy import fpolynomial, fpowers, Fsum, fsum, isfinite


class Tests(TestsBase):

    def testFmath(self):

        # see function _p2 in ellipsoidalVincenty.py
        for x in (0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000):
            p = fpolynomial(x, 16384, 4096, -768, 320, -175) / 16384.0
            a = x / 16384.0 * (4096 + x * (-768 + x * (320 - 175 * x))) + 1
            self.test('fpolynomialA', p, a)

            p = fpolynomial(x, 0, 256, -128, 74, -47) / 1024.0
            b = x / 1024.0 * (256 + x * (-128 + x * (74 - 47 * x)))
            self.test('fpolynomialB', p, b)

        # U{Neumaier<http://wikipedia.org/wiki/Kahan_summation_algorithm>}
        t = 1, 1e101, 1, -1e101
        for _ in range(10):
            s = float(len(t) / 2)  # number of ones
            self.test('sum', sum(t), s, known=True)
            self.test('fsum', fsum(t), s)
            f = Fsum()
            f.fadd2(t)
            self.test('Fsum', f.fsum(), s)
            t += t

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

        self.test('isfinite', isfinite(0), 'True')
        self.test('isfinite', isfinite(1e300), 'True')
        self.test('isfinite', isfinite(-1e300), 'True')
        self.test('isfinite', isfinite(1e1234), 'False')
        self.test('isfinite', isfinite(float('inf')), 'False')
        self.test('isfinite', isfinite(float('nan')), 'False')


if __name__ == '__main__':

    from pygeodesy import fmath

    t = Tests(__file__, __version__, fmath)
    t.testFmath()
    t.results(nl=0)
    t.exit()
