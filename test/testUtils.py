
# -*- coding: utf-8 -*-

# Test base classes.

__all__ = ('Tests',)
__version__ = '17.09.18'

from base import TestsBase

from pygeodesy import fpolynomial, fsum, isfinite


class Tests(TestsBase):

    def testUtils(self):

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
            t += t

        self.test('isfinite', isfinite(0), 'True')
        self.test('isfinite', isfinite(1e300), 'True')
        self.test('isfinite', isfinite(-1e300), 'True')
        self.test('isfinite', isfinite(1e1234), 'False')
        self.test('isfinite', isfinite(float('inf')), 'False')
        self.test('isfinite', isfinite(float('nan')), 'False')


if __name__ == '__main__':

    from pygeodesy import utils  # private

    t = Tests(__file__, __version__, utils)
    t.testUtils()
    t.results(nl=0)
    t.exit()
