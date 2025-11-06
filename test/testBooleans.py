
# -*- coding: utf-8 -*-

# Test L{booleans} module.

__all__ = ('Tests',)
__version__ = '23.03.31'

from bases import TestsBase

from pygeodesy import BooleanFHP, BooleanGH, isenclosedBy


class Tests(TestsBase):

    def testBooleans(self, module):  # MCCABE 15

        self.subtitle(module)
        LatLon = module.LatLon
        areaOf = module.areaOf
        periOf = module.perimeterOf

        p = LatLon(0, 0, height=1.),  LatLon(7, 5, height=2.), LatLon(0, 10, height=3.)  # (0, 0)
        q = LatLon(10, 0, height=1.), LatLon(3, 5, height=2.), LatLon(10, 10, height=3.)  # (5, 0)

        # BooleanGH
        s = BooleanGH(p, name='subject')
        c = BooleanGH(q, name='clipper')
        for n, r, x in (('and',   s & c, 'BooleanGH[4]((lat=5.0, lon=3.5714286, height=1.7142857), (lat=7.0, lon=5.0, height=2.0), (lat=5.0, lon=6.4285714, height=2.2857143), (lat=3.0, lon=5.0, height=2.0))'),
                        ('or',    s | c, 'BooleanGH[6]((lat=5.0, lon=3.5714286, height=1.7142857), (lat=0.0, lon=0.0, height=1.0), (lat=0.0, lon=10.0, height=3.0), (lat=5.0, lon=6.4285714, height=2.2857143), (lat=10.0, lon=10.0, height=3.0), (lat=10.0, lon=0.0, height=1.0))'),
                        ('minus', s - c, 'BooleanGH[5]((lat=5.0, lon=3.5714286, height=1.7142857), (lat=0.0, lon=0.0, height=1.0), (lat=0.0, lon=10.0, height=3.0), (lat=5.0, lon=6.4285714, height=2.2857143), (lat=3.0, lon=5.0, height=2.0))'),
                        ('rev_d', c - s, 'BooleanGH[5]((lat=5.0, lon=3.5714286, height=1.7142857), (lat=10.0, lon=0.0, height=1.0), (lat=10.0, lon=10.0, height=3.0), (lat=5.0, lon=6.4285714, height=2.2857143), (lat=7.0, lon=5.0, height=2.0))')):
            self.test(n, repr(r), x)

        s &= c
        self.test('iand', repr(s), 'BooleanGH[4]((lat=5.0, lon=3.5714286, height=1.7142857), (lat=7.0, lon=5.0, height=2.0), (lat=5.0, lon=6.4285714, height=2.2857143), (lat=3.0, lon=5.0, height=2.0))')
        s = BooleanGH(p)
        s |= c
        self.test('ior', repr(s), 'BooleanGH[6]((lat=5.0, lon=3.5714286, height=1.7142857), (lat=0.0, lon=0.0, height=1.0), (lat=0.0, lon=10.0, height=3.0), (lat=5.0, lon=6.4285714, height=2.2857143), (lat=10.0, lon=10.0, height=3.0), (lat=10.0, lon=0.0, height=1.0))')

        # s -= s
        # self.test('isub', repr(s), '')

        b = BooleanGH(p) + BooleanGH(q)
        self.test('sum', repr(b), 'BooleanGH[2][6]((lat=0.0, lon=0.0, height=1.0), (lat=7.0, lon=5.0, height=2.0), (lat=0.0, lon=10.0, height=3.0), (lat=10.0, lon=0.0, height=1.0, clipid=1), (lat=3.0, lon=5.0, height=2.0, clipid=1), (lat=10.0, lon=10.0, height=3.0, clipid=1))')
        t = BooleanGH(q) + BooleanGH(p)
        self.test('GH ==', b == t, True, nl=1)
        self.test('equalTo', b.isequalTo(t, eps=1e-9), True)

        self.test(areaOf.__name__, areaOf(b) == areaOf(t), True)
        self.test('enclosed', isenclosedBy(LatLon(9, 5), b), True)
        self.test('enclosed', isenclosedBy(LatLon(5, 5), t), False)
        self.test(periOf.__name__, periOf(b, closed=True) == periOf(t, closed=True), True)

        t = tuple(b.toLatLon(LatLon))[:3]
        self.test('toLatLon[0:3]', repr(t), '(LatLon(00°00′00.0″N, 000°00′00.0″E, +1.00m), LatLon(07°00′00.0″N, 005°00′00.0″E, +2.00m), LatLon(00°00′00.0″N, 010°00′00.0″E, +3.00m))', nl=1)
        t = tuple(b.toLatLon())[-3:]
        self.test('toLatLon[-3:]', repr(t), '((lat=10.0, lon=0.0, height=1.0, clipid=1), (lat=3.0, lon=5.0, height=2.0, clipid=1), (lat=10.0, lon=10.0, height=3.0, clipid=1))', nt=1)

        # Boolean FHP
        s = BooleanFHP(p, name='subject')
        c = BooleanFHP(q, name='clipper')
        for n, r, x in (('and',   s & c, 'BooleanFHP[4]((lat=7.0, lon=5.0, height=2.0), (lat=5.0, lon=6.4285714, height=2.2857143), (lat=3.0, lon=5.0, height=2.0), (lat=5.0, lon=3.5714286, height=1.7142857))'),
                        ('or',    s | c, 'BooleanFHP[6]((lat=0.0, lon=0.0, height=1.0), (lat=0.0, lon=10.0, height=3.0), (lat=5.0, lon=6.4285714, height=2.2857143), (lat=10.0, lon=10.0, height=3.0), (lat=10.0, lon=0.0, height=1.0), (lat=5.0, lon=3.5714286, height=1.7142857))')):
            self.test(n, repr(r), x)

        s &= c
        self.test('iand', repr(s), 'BooleanFHP[4]((lat=7.0, lon=5.0, height=2.0), (lat=5.0, lon=6.4285714, height=2.2857143), (lat=3.0, lon=5.0, height=2.0), (lat=5.0, lon=3.5714286, height=1.7142857))')
        s = BooleanFHP(p)
        s |= c
        self.test('ior', repr(s), 'BooleanFHP[6]((lat=0.0, lon=0.0, height=1.0), (lat=0.0, lon=10.0, height=3.0), (lat=5.0, lon=6.4285714, height=2.2857143), (lat=10.0, lon=10.0, height=3.0), (lat=10.0, lon=0.0, height=1.0), (lat=5.0, lon=3.5714286, height=1.7142857))')

        # s -= s
        # self.test('isub', repr(s), '')

        b = BooleanFHP(p) + BooleanFHP(q)
        self.test('sum', repr(b), 'BooleanFHP[2][6]((lat=0.0, lon=0.0, height=1.0), (lat=7.0, lon=5.0, height=2.0), (lat=0.0, lon=10.0, height=3.0), (lat=10.0, lon=0.0, height=1.0, clipid=1), (lat=3.0, lon=5.0, height=2.0, clipid=1), (lat=10.0, lon=10.0, height=3.0, clipid=1))')
        t = BooleanFHP(q) + BooleanFHP(p)
        self.test('FHP ==', b == t, True, nl=1)
        self.test('equalTo', b.isequalTo(t, eps=1e-9), True)

        self.test(areaOf.__name__, areaOf(b) == areaOf(t), True)
        self.test('enclosed', isenclosedBy(LatLon(9, 5), b), True)
        self.test('enclosed', isenclosedBy(LatLon(5, 5), t), False)
        self.test(periOf.__name__, periOf(b, closed=True) == periOf(t, closed=True), True)

        t = tuple(b.toLatLon(LatLon))[:3]
        self.test('toLatLon[0:3]', repr(t), '(LatLon(00°00′00.0″N, 000°00′00.0″E, +1.00m), LatLon(07°00′00.0″N, 005°00′00.0″E, +2.00m), LatLon(00°00′00.0″N, 010°00′00.0″E, +3.00m))', nl=1)
        t = tuple(b.toLatLon())[-3:]
        self.test('toLatLon[-3:]', repr(t), '((lat=10.0, lon=0.0, height=1.0, clipid=1), (lat=3.0, lon=5.0, height=2.0, clipid=1), (lat=10.0, lon=10.0, height=3.0, clipid=1))')


if __name__ == '__main__':

    from pygeodesy import ellipsoidalExact, ellipsoidalVincenty, \
                          sphericalNvector, sphericalTrigonometry

    t = Tests(__file__, __version__)
    t.testBooleans(ellipsoidalExact)
#   t.testBooleans(ellipsoidalNvector)
    t.testBooleans(ellipsoidalVincenty)
    t.testBooleans(sphericalNvector)
    t.testBooleans(sphericalTrigonometry)
    t.results()
    t.exit()
