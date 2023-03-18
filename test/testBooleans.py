
# -*- coding: utf-8 -*-

# Test base classes.

__all__ = ('Tests',)
__version__ = '23.03.18'

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

        b = BooleanFHP(p) + BooleanFHP(q)
        self.test('sum', repr(b), 'BooleanFHP[2][6]((lat=0, lon=0, height=1), (lat=7, lon=5, height=2), (lat=0, lon=10, height=3), (lat=10, lon=0, height=1, clipid=1), (lat=3, lon=5, height=2, clipid=1), (lat=10, lon=10, height=3, clipid=1))')
        self.test('==', b == (BooleanFHP(q) + BooleanFHP(p)), True)

        t = tuple(b.toLatLon(LatLon))[:3]
        self.test('toLatLon[0:3]', repr(t), '(LatLon(00°00′00.0″N, 000°00′00.0″E, +1.00m), LatLon(07°00′00.0″N, 005°00′00.0″E, +2.00m), LatLon(00°00′00.0″N, 010°00′00.0″E, +3.00m))', nl=1)
        t = tuple(b.toLatLon())[-3:]
        self.test('toLatLon[-3:]', repr(t), '((lat=10, lon=0, height=1, clipid=1), (lat=3, lon=5, height=2, clipid=1), (lat=10, lon=10, height=3, clipid=1))')


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
