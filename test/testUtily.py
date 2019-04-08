
# -*- coding: utf-8 -*-

# Test base classes.

__all__ = ('Tests',)
__version__ = '19.04.05'

from base import TestsBase, geographiclib

from pygeodesy import EPS, fStr, map1, sincos2, sincos2d, splice, unroll180

from math import cos, radians, sin

if geographiclib:
    from geographiclib.geomath import Math
    sincosd = Math.sincosd
    del Math
else:
    sincosd = None


class Tests(TestsBase):

    def testUtily(self):

        self.test('unroll180', fStr(unroll180(-90, 110, wrap=True)), '-160.0, -250.0')
        self.test('unroll180', fStr(unroll180(-90, 110, wrap=False)), '200.0, 110.0')

        self.test('unroll180', fStr(unroll180(-90, 830, wrap=True)), '-160.0, -250.0')
        self.test('unroll180', fStr(unroll180(-90, 830, wrap=False)), '920.0, 830.0')

        self.test('unroll180', fStr(unroll180(-110, 90, wrap=True)), '-160.0, -270.0')
        self.test('unroll180', fStr(unroll180(-110, 90, wrap=False)), '200.0, 90.0')

        self.test('unroll180', fStr(unroll180(-830, 90, wrap=True)), '-160.0, -990.0')
        self.test('unroll180', fStr(unroll180(-830, 90, wrap=False)), '920.0, 90.0')

        e = d = g = f = 0
        for a in range(-1000, 1000):
            a *= 0.47
            r = radians(a)
            sr, cr = sin(r), cos(r)

            s, c = sincos2(r)
            e = max(e, abs(sr - s), abs(cr - c))

            sd, cd = sincos2d(a)
            d = max(d, abs(sr - sd), abs(cr - cd))
            if sincosd:  # compare with geographiclib
                s, c = sincosd(a)
                g = max(g, abs(sr - s), abs(cr - c))
                f = max(f, abs(sd - s), abs(cd - c))

        EPS_ = EPS * 8
        self.test('sincos2',  e, EPS, known=e < EPS_)
        self.test('sincos2d', d, EPS, known=d < EPS_)
        if sincosd:
            self.test('sincosd ', g, EPS, known=g < EPS_)
            self.test('sincos*d', f, EPS, known=f < EPS_)

        a, b = splice(range(10))  # PYCHOK false
        self.test('splice', (a, b), map1(type(a), (0, 2, 4, 6, 8), (1, 3, 5, 7, 9)))
        a, b, c = splice(range(10), n=3)  # PYCHOK false
        self.test('splice', (a, b, c), map1(type(a), (0, 3, 6, 9), (1, 4, 7), (2, 5, 8)))
        a, b, c = splice(range(10), n=3, fill=-1)  # PYCHOK false
        self.test('splice', (a, b, c), map1(type(a), (0, 3, 6, 9), (1, 4, 7, -1), (2, 5, 8, -1)))
        t = tuple(splice(range(12), n=5))  # PYCHOK false
        self.test('splice', t, map1(type(t[0]), (0, 5, 10), (1, 6, 11), (2, 7), (3, 8), (4, 9)))


if __name__ == '__main__':

    from pygeodesy import utily  # private

    t = Tests(__file__, __version__, utily)
    t.testUtily()
    t.results()
    t.exit()
