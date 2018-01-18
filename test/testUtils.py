
# -*- coding: utf-8 -*-

# Test base classes.

__all__ = ('Tests',)
__version__ = '18.01.16'

from base import TestsBase

from pygeodesy import R_M, fpolynomial, fStr, fsum, heightOf, horizon, \
                           isfinite, unroll180


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

        self.test('heightof0', heightOf(0, R_M), 2638958.23912, fmt='%.5f')
        self.test('heightof45', heightOf(45, R_M), 5401080.43931, fmt='%.5f')
        self.test('heightof90', heightOf(90, R_M), R_M)
        self.test('heightof135', heightOf(135, R_M), 5401080.43931, fmt='%.5f')

        self.test('horizon0', horizon(0), 0.0)
        self.test('horizon10Km',  horizon(10000), '357099.672', fmt='%.3f')
        self.test('horizon30Kft', horizon(10000, refraction=True), '392310.704', fmt='%.3f')
        self.test('horizon10Kft', horizon( 3000, refraction=True), '214877.422', fmt='%.3f')

        self.test('isfinite', isfinite(0), 'True')
        self.test('isfinite', isfinite(1e300), 'True')
        self.test('isfinite', isfinite(-1e300), 'True')
        self.test('isfinite', isfinite(1e1234), 'False')
        self.test('isfinite', isfinite(float('inf')), 'False')
        self.test('isfinite', isfinite(float('nan')), 'False')

        self.test('unroll180', fStr(unroll180(-90, 110, wrap=True)), '-160.0, -250.0')
        self.test('unroll180', fStr(unroll180(-90, 110, wrap=False)), '200.0, 110.0')

        self.test('unroll180', fStr(unroll180(-90, 830, wrap=True)), '-160.0, -250.0')
        self.test('unroll180', fStr(unroll180(-90, 830, wrap=False)), '920.0, 830.0')

        self.test('unroll180', fStr(unroll180(-110, 90, wrap=True)), '-160.0, -270.0')
        self.test('unroll180', fStr(unroll180(-110, 90, wrap=False)), '200.0, 90.0')

        self.test('unroll180', fStr(unroll180(-830, 90, wrap=True)), '-160.0, -990.0')
        self.test('unroll180', fStr(unroll180(-830, 90, wrap=False)), '920.0, 90.0')


if __name__ == '__main__':

    from pygeodesy import utils  # private

    t = Tests(__file__, __version__, utils)
    t.testUtils()
    t.results(nl=0)
    t.exit()
