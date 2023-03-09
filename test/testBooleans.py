
# -*- coding: utf-8 -*-

# Test base classes.

__all__ = ('Tests',)
__version__ = '23.03.09'

from base import TestsBase

from pygeodesy import BooleanFHP, BooleanGH


class Tests(TestsBase):

    def testBooleans(self, module):  # MCCABE 15

        self.subtitle(module)
        LatLon = module.LatLon

        p = LatLon(0, 0, height=1.),  LatLon(7, 5, height=2.), LatLon(0, 10, height=3.)  # (0, 0)
        q = LatLon(10, 0, height=1.), LatLon(3, 5, height=2.), LatLon(10, 10, height=3.)  # (5, 0)

        s = BooleanGH(p, name='subject')
        c = BooleanGH(q, name='clipper')
        for n, r, x in (('and',   s & c, 'BooleanGH[4]((lat=5, lon=3.5714286, height=1.7142857), (lat=7, lon=5, height=2), (lat=5, lon=6.4285714, height=2.2857143), (lat=3, lon=5, height=2))'),
                        ('or',    s | c, 'BooleanGH[6]((lat=5, lon=3.5714286, height=1.7142857), (lat=0, lon=0, height=1), (lat=0, lon=10, height=3), (lat=5, lon=6.4285714, height=2.2857143), (lat=10, lon=10, height=3), (lat=10, lon=0, height=1))'),
                        ('minus', s - c, 'BooleanGH[5]((lat=5, lon=3.5714286, height=1.7142857), (lat=0, lon=0, height=1), (lat=0, lon=10, height=3), (lat=5, lon=6.4285714, height=2.2857143), (lat=3, lon=5, height=2))'),
                        ('rev_d', c - s, 'BooleanGH[5]((lat=5, lon=3.5714286, height=1.7142857), (lat=10, lon=0, height=1), (lat=10, lon=10, height=3), (lat=5, lon=6.4285714, height=2.2857143), (lat=7, lon=5, height=2))')):
            self.test(n, repr(r), x)

        s &= c
        self.test('iand', repr(s), 'BooleanGH[4]((lat=5, lon=3.5714286, height=1.7142857), (lat=7, lon=5, height=2), (lat=5, lon=6.4285714, height=2.2857143), (lat=3, lon=5, height=2))')
        s = BooleanGH(p)
        s |= c
        self.test('ior', repr(s), 'BooleanGH[6]((lat=5, lon=3.5714286, height=1.7142857), (lat=0, lon=0, height=1), (lat=0, lon=10, height=3), (lat=5, lon=6.4285714, height=2.2857143), (lat=10, lon=10, height=3), (lat=10, lon=0, height=1))')

        # s -= s
        # self.test('isub', repr(s), '')

        s = BooleanFHP(p, name='subject')
        c = BooleanFHP(q, name='clipper')
        for n, r, x in (('and',   s & c, 'BooleanFHP[4]((lat=7, lon=5, height=2), (lat=5, lon=6.4285714, height=2.2857143), (lat=3, lon=5, height=2), (lat=5, lon=3.5714286, height=1.7142857))'),
                        ('or',    s | c, 'BooleanFHP[6]((lat=0, lon=0, height=1), (lat=0, lon=10, height=3), (lat=5, lon=6.4285714, height=2.2857143), (lat=10, lon=10, height=3), (lat=10, lon=0, height=1), (lat=5, lon=3.5714286, height=1.7142857))')):
            self.test(n, repr(r), x)

        s &= c
        self.test('iand', repr(s), 'BooleanFHP[4]((lat=7, lon=5, height=2), (lat=5, lon=6.4285714, height=2.2857143), (lat=3, lon=5, height=2), (lat=5, lon=3.5714286, height=1.7142857))')
        s = BooleanFHP(p)
        s |= c
        self.test('ior', repr(s), 'BooleanFHP[6]((lat=0, lon=0, height=1), (lat=0, lon=10, height=3), (lat=5, lon=6.4285714, height=2.2857143), (lat=10, lon=10, height=3), (lat=10, lon=0, height=1), (lat=5, lon=3.5714286, height=1.7142857))')

        # s -= s
        # self.test('isub', repr(s), '')


if __name__ == '__main__':

    from pygeodesy import ellipsoidalNvector, ellipsoidalVincenty, \
                          sphericalNvector, sphericalTrigonometry

    t = Tests(__file__, __version__)
    t.testBooleans(ellipsoidalNvector)
    t.testBooleans(ellipsoidalVincenty)
    t.testBooleans(sphericalNvector)
    t.testBooleans(sphericalTrigonometry)
    t.results()
    t.exit()
