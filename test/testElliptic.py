
# -*- coding: utf-8 -*-

u'''Test Elliptic Python implementation.
'''

__all__ = ('Tests',)
__version__ = '19.05.19'

from base import TestsBase

from pygeodesy import elliptic, EPS


class Tests(TestsBase):

    def testElliptic(self):

        _RC = elliptic._RC
        _RD = elliptic._RD
        _RF = elliptic._RF
        _RJ = elliptic._RJ

        eps = EPS * 4
        self.test('eps', eps, eps, fmt='%.12e')

        y = 0.1
        for x in range(2, 100):
            x *= 0.01
            rc = _RC(x, y)
            rf = _RF(x, y, y)
            k = abs(rc - rf) < eps
            t = '_RC, _RF(%.3f, ...)' % (x,)
            self.test(t, rc, rf, fmt='%.12f', known=k)
            y = x

        for kp2 in range(1, 100):
            kp2 *= 0.01
            rd = _RD(0, kp2, 1)
            rj = _RJ(0, kp2, 1, 1)
            k = abs(rd - rj) < eps
            t = '_RD, _RJ(%.3f, ...)' % (kp2,)
            self.test(t, rd, rj, fmt='%.12f', known=k)

        self.test('eps', eps, eps, fmt='%.12e')


if __name__ == '__main__':

    t = Tests(__file__, __version__)
    t.testElliptic()
    t.results()
    t.exit()
