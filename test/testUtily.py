
# -*- coding: utf-8 -*-

# Test base classes.

__all__ = ('Tests',)
__version__ = '20.10.20'

from base import TestsBase, geographiclib

from pygeodesy import EPS, PI, PI2, PI_2, \
                      acre2ha, acre2m2, atan2d, chain2m, \
                      degrees90, degrees180, degrees360, degrees2m, \
                      fathom2m, ft2m, furlong2m, isPoints2, \
                      m2degrees, m2ft, m2yard, \
                      radiansPI, radiansPI2, radiansPI_2, \
                      sincos2, sincos2d, unroll180, \
                      wrap90, wrap180, wrap360, \
                      wrapPI, wrapPI2, wrapPI_2, \
                      yard2m, fstr  # DEPRECATED, use fstr

from math import cos, radians, sin

if geographiclib:
    from geographiclib.geomath import Math
    atand   = Math.atan2d
    sincosd = Math.sincosd
else:
    Math = None


class Tests(TestsBase):

    def testUtily(self):

        # Python 2.6.9 on Travis Ubuntu 14.04 produces -0.0

        self.test('degrees90(PI_2)', degrees90(PI_2), 90.0)
        self.test('degrees90(PI)',   degrees90(PI), -180.0)  # XXX
        self.test('degrees90(PI2)',  degrees90(PI2),   0.0)
        self.test('degrees90(-PI_2)',    degrees90(-PI_2), -90.0)
        self.test('degrees90(-PI)',      degrees90(-PI),  -180.0)  # XXX
        self.test('degrees90(-PI2)', abs(degrees90(-PI2)),   0.0)  # -0.0

        self.test('degrees180(PI_2)', degrees180(PI_2), 90.0)
        self.test('degrees180(PI)',   degrees180(PI),  180.0)  # XXX
        self.test('degrees180(PI2)',  degrees180(PI2),   0.0)
        self.test('degrees180(-PI_2)',    degrees180(-PI_2), -90.0)
        self.test('degrees180(-PI)',      degrees180(-PI),  -180.0)  # XXX
        self.test('degrees180(-PI2)', abs(degrees180(-PI2)),   0.0)  # -0.0

        self.test('degrees360(PI_2)', degrees360(PI_2), 90.0)
        self.test('degrees360(PI)',   degrees360(PI),  180.0)  # XXX
        self.test('degrees360(PI2)',  degrees360(PI2),   0.0)
        self.test('degrees360(-PI_2)',    degrees360(-PI_2), 270.0)
        self.test('degrees360(-PI)',      degrees360(-PI),   180.0)  # XXX
        self.test('degrees360(-PI2)', abs(degrees360(-PI2)),   0.0)  # -0.0

        self.test('radiansPI_2(90)',  radiansPI_2(90), PI_2)
        self.test('radiansPI_2(180)', radiansPI_2(180), -PI)
        self.test('radiansPI_2(360)', radiansPI_2(360), 0.0)
        self.test('radiansPI_2(-90)',      radiansPI_2(-90), -PI_2)
        self.test('radiansPI_2(-180)',     radiansPI_2(-180),  -PI)
        self.test('radiansPI_2(-360)', abs(radiansPI_2(-360)), 0.0)  # -0.0

        self.test('radiansPI(90)',  radiansPI(90), PI_2)
        self.test('radiansPI(180)', radiansPI(180),  PI)
        self.test('radiansPI(360)', radiansPI(360), 0.0)
        self.test('radiansPI(-90)',      radiansPI(-90), -PI_2)
        self.test('radiansPI(-180)',     radiansPI(-180),  -PI)
        self.test('radiansPI(-360)', abs(radiansPI(-360)), 0.0)  # -0.0

        self.test('radiansPI2(90)',  radiansPI2(90), PI_2)
        self.test('radiansPI2(180)', radiansPI2(180),  PI)
        self.test('radiansPI2(360)', radiansPI2(360), 0.0)
        self.test('radiansPI2(-90)',      radiansPI2(-90), PI_2+PI)
        self.test('radiansPI2(-180)',     radiansPI2(-180),     PI)
        self.test('radiansPI2(-360)', abs(radiansPI2(-360)),   0.0)  # -0.0

        self.test('wrap90(90)',   wrap90(90),     90.0)
        self.test('wrap90(180)',  wrap90(180),  -180.0)
        self.test('wrap90(360)',  wrap90(360),     0.0)
        self.test('wrap90(-90)',  wrap90(-90),   -90.0)
        self.test('wrap90(-180)', wrap90(-180), -180.0)
        self.test('wrap90(-360)', wrap90(-360),    0.0)

        self.test('wrap180(90)',   wrap180(90),     90.0)
        self.test('wrap180(180)',  wrap180(180),   180.0)
        self.test('wrap180(360)',  wrap180(360),     0.0)
        self.test('wrap180(-90)',  wrap180(-90),   -90.0)
        self.test('wrap180(-180)', wrap180(-180), -180.0)
        self.test('wrap180(-360)', wrap180(-360),    0.0)

        self.test('wrap360(90)',   wrap360(90),    90.0)
        self.test('wrap360(180)',  wrap360(180),  180.0)
        self.test('wrap360(360)',  wrap360(360),    0.0)
        self.test('wrap360(-90)',  wrap360(-90),  270.0)
        self.test('wrap360(-180)', wrap360(-180), 180.0)
        self.test('wrap360(-360)', wrap360(-360),   0.0)

        self.test('wrapPI_2(PI_2)', wrapPI_2(PI_2), PI_2)
        self.test('wrapPI_2(PI)',   wrapPI_2(PI),    -PI)  # XXX
        self.test('wrapPI_2(PI2)',  wrapPI_2(PI2),   0.0)
        self.test('wrapPI_2(-PI_2)',    wrapPI_2(-PI_2), -PI_2)
        self.test('wrapPI_2(-PI)',      wrapPI_2(-PI),     -PI)  # XXX
        self.test('wrapPI_2(-PI2)', abs(wrapPI_2(-PI2)),   0.0)

        self.test('wrapPI(PI_2)', wrapPI(PI_2), PI_2)
        self.test('wrapPI(PI)',   wrapPI(PI),     PI)  # XXX
        self.test('wrapPI(PI2)',  wrapPI(PI2),   0.0)
        self.test('wrapPI(-PI_2)',    wrapPI(-PI_2), -PI_2)
        self.test('wrapPI(-PI)',      wrapPI(-PI),     -PI)  # XXX
        self.test('wrapPI(-PI2)', abs(wrapPI(-PI2)),   0.0)  # -0.0

        self.test('wrapPI2(PI_2)', wrapPI2(PI_2), PI_2)
        self.test('wrapPI2(PI)',   wrapPI2(PI),     PI)  # XXX
        self.test('wrapPI2(PI2)',  wrapPI2(PI2),   0.0)
        self.test('wrapPI2(-PI_2)',    wrapPI2(-PI_2), PI_2+PI)
        self.test('wrapPI2(-PI)',      wrapPI2(-PI),        PI)  # XXX
        self.test('wrapPI2(-PI2)', abs(wrapPI2(-PI2)),     0.0)  # -0.0

        self.test('unroll180', fstr(unroll180(-90, 110, wrap=True)), '-160.0, -250.0')
        self.test('unroll180', fstr(unroll180(-90, 110, wrap=False)), '200.0, 110.0')

        self.test('unroll180', fstr(unroll180(-90, 830, wrap=True)), '-160.0, -250.0')
        self.test('unroll180', fstr(unroll180(-90, 830, wrap=False)), '920.0, 830.0')

        self.test('unroll180', fstr(unroll180(-110, 90, wrap=True)), '-160.0, -270.0')
        self.test('unroll180', fstr(unroll180(-110, 90, wrap=False)), '200.0, 90.0')

        self.test('unroll180', fstr(unroll180(-830, 90, wrap=True)), '-160.0, -990.0')
        self.test('unroll180', fstr(unroll180(-830, 90, wrap=False)), '920.0, 90.0')

        e = d = g = f = t = 0
        for a in range(-1000, 1000, 2):
            a *= 0.47
            r = radians(a)
            sr, cr = sin(r), cos(r)

            s, c = sincos2(r)
            e = max(e, abs(sr - s), abs(cr - c))

            sd, cd = sincos2d(a)
            d = max(d, abs(sr - sd), abs(cr - cd))
            if Math:  # compare with geographiclib
                t = max(t, abs(atan2d(sr, cr) - atand(sr, cr)))
                s, c = sincosd(a)
                g = max(g, abs(sr - s), abs(cr - c))
                f = max(f, abs(sd - s), abs(cd - c))

        EPS_ = EPS * 8
        self.test('sincos2',  e, EPS_, known=e < EPS_)
        self.test('sincos2d', d, EPS_, known=d < EPS_)
        if Math:
            self.test('atand',    t, EPS,  known=t < EPS)
            self.test('sincosd ', g, EPS_, known=g < EPS_)
            self.test('sincos*d', f, EPS_, known=f < EPS_)

        # <https://www.CivilGeo.com/when-a-foot-isnt-really-a-foot/>
        self.test('iFt2m', ft2m( 614963.91), 187441, fmt='%.0f')
        self.test('iFt2m', ft2m(2483759.84), 757050, fmt='%.0f')
        self.test('sFt2m', ft2m( 614962.68, usurvey=True), 187441, fmt='%.0f')
        self.test('sFt2m', ft2m(2483754.87, usurvey=True), 757050, fmt='%.0f')

        self.test('m2iFt', m2ft(187441),  614963.91, prec=2)
        self.test('m2iFt', m2ft(757050), 2483759.84, prec=2)
        self.test('m2sFt', m2ft(187441, usurvey=True),  614962.68, prec=2)
        self.test('m2sFt', m2ft(757050, usurvey=True), 2483754.87, prec=2)

        for f, m in ((m2yard,      '1.093613'),
                     (acre2ha,     '0.404686'), (acre2m2, '4046.856422'),
                     (chain2m,    '20.116800'), (fathom2m,   '1.828800'),
                     (furlong2m, '201.168000'), (yard2m,     '0.914400')):
            self.test(f.__name__, f(1), m, prec=6)

        self.test('degrees2m', fstr(degrees2m(90), prec=4),        '10007557.1761')
        self.test('degrees2m', fstr(degrees2m(90, lat=30), prec=4), '8666798.7443')
        self.test('m2degrees', fstr(m2degrees(degrees2m(90)), prec=1),   '90.0')

        self.test('degrees2m', fstr(degrees2m(180), prec=4),          '20015114.3522')
        self.test('degrees2m', fstr(degrees2m(180, lat=3-0), prec=4), '19987684.3336')
        self.test('m2degrees', fstr(m2degrees(degrees2m(180)), prec=1),    '180.0')

        t = 'm2degrees2m(%s, lat=%s)'
        for a in range(0, 90, 7):
            d = m2degrees(degrees2m(45, lat=a), lat=a)
            self.test(t % (45, a), d, '45.00', prec=2)

        self.test('isPoints2', isPoints2(None), False)


if __name__ == '__main__':

    from pygeodesy import utily  # private

    t = Tests(__file__, __version__, utily)
    t.testUtily()
    t.results()
    t.exit()
