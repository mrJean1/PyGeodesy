
# -*- coding: utf-8 -*-

# Test base classes.

__all__ = ('Tests',)
__version__ = '18.07.21'

from base import TestsBase
from random import random, gauss, shuffle

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

        # <http://code.activestate.com/recipes/393090>
        t = 1.00, 0.00500, 0.00000000000100
        s = 1.00500000000100
        self.test('sum', sum(t), s, known=True)
        self.test('fsum', fsum(t), s)
        f = Fsum()
        f.fadd2(t)
        self.test('Fsum', f.fsum(), s)

        # <http://github.com/ActiveState/code/blob/master/recipes/Python/
        #       393090_Binary_floating_point_summatiaccurate_full/recipe-393090.py>
        for _ in range(100):
            t = [7, 1e100, -7, -1e100, -9e-20, 8e-20] * 10
            s = 0
            for _ in range(20):
                v = gauss(0, random())**7 - s
                s += v
                t.append(v)
            shuffle(t)
            s = fsum(t)
            # self.test('fsum', s, s)
            f = Fsum()
            f.fadd2(t)
            self.test('Fsum', f.fsum(), s)

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
