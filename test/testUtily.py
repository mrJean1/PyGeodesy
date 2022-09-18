
# -*- coding: utf-8 -*-

# Test base classes.

__all__ = ('Tests',)
__version__ = '22.09.15'

from base import TestsBase, geographiclib

from pygeodesy import EPS, INF, NEG0, NINF, PI, PI2, PI_2, PI3_2, \
                      acre2ha, acre2m2, atan2d, chain2m, cot_, cotd_, \
                      degrees90, degrees180, degrees360, degrees2m, \
                      fathom2m, ft2m, furlong2m, \
                      grades400, degrees2grades, grades2degrees, grades2radians, \
                      isPoints2, map1, \
                      m2chain, m2degrees, m2fathom, m2ft, m2furlong, m2toise, m2yard, \
                      radiansPI, radiansPI2, radiansPI_2, \
                      sincos2, sincos2d, sincostan3, tan_2, unroll180, \
                      wrap90, wrap180, wrap360, wrapPI, wrapPI2, wrapPI_2, \
                      toise2m, yard2m, fstr  # DEPRECATED, use fstr

from math import cos, radians, sin, tan

if geographiclib:
    from geographiclib.geomath import Math
    atand   = Math.atan2d
    sincosd = Math.sincosd
else:
    Math = None


class Tests(TestsBase):

    def testUtily(self):

        # Python 2.6.9 on Travis Ubuntu 14.04 produces -0.0

        self.test('degrees90(PI_2)', degrees90(PI_2), 90.0, nl=1)
        self.test('degrees90(PI)',   degrees90(PI), -180.0)  # XXX
        self.test('degrees90(PI2)',  degrees90(PI2),   0.0)
        self.test('degrees90(-PI_2)',    degrees90(-PI_2), -90.0)
        self.test('degrees90(-PI)',      degrees90(-PI),  -180.0)  # XXX
        self.test('degrees90(-PI2)', abs(degrees90(-PI2)),   0.0)  # -0.0

        self.test('degrees180(PI_2)', degrees180(PI_2), 90.0, nl=1)
        self.test('degrees180(PI)',   degrees180(PI),  180.0)  # XXX
        self.test('degrees180(PI2)',  degrees180(PI2),   0.0)
        self.test('degrees180(-PI_2)',    degrees180(-PI_2), -90.0)
        self.test('degrees180(-PI)',      degrees180(-PI),  -180.0)  # XXX
        self.test('degrees180(-PI2)', abs(degrees180(-PI2)),   0.0)  # -0.0

        self.test('degrees360(PI_2)', degrees360(PI_2), 90.0, nl=1)
        self.test('degrees360(PI)',   degrees360(PI),  180.0)  # XXX
        self.test('degrees360(PI2)',  degrees360(PI2),   0.0)
        self.test('degrees360(-PI_2)',    degrees360(-PI_2), 270.0)
        self.test('degrees360(-PI)',      degrees360(-PI),   180.0)  # XXX
        self.test('degrees360(-PI2)', abs(degrees360(-PI2)),   0.0)  # -0.0

        self.test('degrees2grades(90)',  degrees2grades(90),  100.0, nl=1)
        self.test('degrees2grades(180)', degrees2grades(180), 200.0)
        self.test('degrees2grades(360)', degrees2grades(360), 400.0)
        self.test('degrees2grades(-90)',  degrees2grades(-90),  -100.0)
        self.test('degrees2grades(-180)', degrees2grades(-180), -200.0)
        self.test('degrees2grades(-360)', degrees2grades(-360), -400.0)

        self.test('grades400(PI_2)', grades400(PI_2), 100.0, nl=1)
        self.test('grades400(PI)',   grades400(PI),   200.0)  # XXX
        self.test('grades400(PI2)',  grades400(PI2),    0.0)
        self.test('grades400(-PI_2)',    grades400(-PI_2), 300.0)
        self.test('grades400(-PI)',      grades400(-PI),   200.0)  # XXX
        self.test('grades400(-PI2)', abs(grades400(-PI2)),   0.0)  # -0.0

        self.test('grades2degrees(100)', grades2degrees(100),  90.0, nl=1)
        self.test('grades2degrees(200)', grades2degrees(200), 180.0)
        self.test('grades2degrees(400)', grades2degrees(400), 360.0)
        self.test('grades2degrees(-100)', grades2degrees(-100),  -90.0)
        self.test('grades2degrees(-200)', grades2degrees(-200), -180.0)
        self.test('grades2degrees(-400)', grades2degrees(-400), -360.0)

        self.test('grades2radians(100)', grades2radians(100), PI_2, nl=1)
        self.test('grades2radians(200)', grades2radians(200), PI)
        self.test('grades2radians(400)', grades2radians(400), PI2)
        self.test('grades2radians(-100)', grades2radians(-100), -PI_2)
        self.test('grades2radians(-200)', grades2radians(-200), -PI)
        self.test('grades2radians(-400)', grades2radians(-400), -PI2)

        self.test('radiansPI_2(90)',  radiansPI_2(90), PI_2, nl=1)
        self.test('radiansPI_2(180)', radiansPI_2(180), -PI)
        self.test('radiansPI_2(360)', radiansPI_2(360), 0.0)
        self.test('radiansPI_2(-90)',      radiansPI_2(-90), -PI_2)
        self.test('radiansPI_2(-180)',     radiansPI_2(-180), -PI)
        self.test('radiansPI_2(-360)', abs(radiansPI_2(-360)), 0.0)  # -0.0

        self.test('radiansPI(90)',  radiansPI(90), PI_2, nl=1)
        self.test('radiansPI(180)', radiansPI(180),  PI)
        self.test('radiansPI(360)', radiansPI(360), 0.0)
        self.test('radiansPI(-90)',      radiansPI(-90), -PI_2)
        self.test('radiansPI(-180)',     radiansPI(-180),  -PI)
        self.test('radiansPI(-360)', abs(radiansPI(-360)), 0.0)  # -0.0

        self.test('radiansPI2(90)',  radiansPI2(90), PI_2, nl=1)
        self.test('radiansPI2(180)', radiansPI2(180),  PI)
        self.test('radiansPI2(360)', radiansPI2(360), 0.0)
        self.test('radiansPI2(-90)',      radiansPI2(-90), PI_2+PI)
        self.test('radiansPI2(-180)',     radiansPI2(-180),     PI)
        self.test('radiansPI2(-360)', abs(radiansPI2(-360)),   0.0)  # -0.0

        self.test('wrap90(90)',   wrap90(90),     90.0, nl=1)
        self.test('wrap90(180)',  wrap90(180),  -180.0)
        self.test('wrap90(360)',  wrap90(360),     0.0)
        self.test('wrap90(-90)',  wrap90(-90),   -90.0)
        self.test('wrap90(-180)', wrap90(-180), -180.0)
        self.test('wrap90(-360)', wrap90(-360),    0.0)

        self.test('wrap180(90)',   wrap180(90),     90.0, nl=1)
        self.test('wrap180(180)',  wrap180(180),   180.0)
        self.test('wrap180(360)',  wrap180(360),     0.0)
        self.test('wrap180(-90)',  wrap180(-90),   -90.0)
        self.test('wrap180(-180)', wrap180(-180), -180.0)
        self.test('wrap180(-360)', wrap180(-360),    0.0)

        self.test('wrap360(90)',   wrap360(90),    90.0, nl=1)
        self.test('wrap360(180)',  wrap360(180),  180.0)
        self.test('wrap360(360)',  wrap360(360),    0.0)
        self.test('wrap360(-90)',  wrap360(-90),  270.0)
        self.test('wrap360(-180)', wrap360(-180), 180.0)
        self.test('wrap360(-360)', wrap360(-360),   0.0)

        self.test('wrapPI_2(PI_2)', wrapPI_2(PI_2), PI_2, nl=1)
        self.test('wrapPI_2(PI)',   wrapPI_2(PI),    -PI)  # XXX
        self.test('wrapPI_2(PI2)',  wrapPI_2(PI2),   0.0)
        self.test('wrapPI_2(-PI_2)',    wrapPI_2(-PI_2), -PI_2)
        self.test('wrapPI_2(-PI)',      wrapPI_2(-PI),     -PI)  # XXX
        self.test('wrapPI_2(-PI2)', abs(wrapPI_2(-PI2)),   0.0)

        self.test('wrapPI(PI_2)', wrapPI(PI_2), PI_2, nl=1)
        self.test('wrapPI(PI)',   wrapPI(PI),     PI)  # XXX
        self.test('wrapPI(PI2)',  wrapPI(PI2),   0.0)
        self.test('wrapPI(-PI_2)',    wrapPI(-PI_2), -PI_2)
        self.test('wrapPI(-PI)',      wrapPI(-PI),     -PI)  # XXX
        self.test('wrapPI(-PI2)', abs(wrapPI(-PI2)),   0.0)  # -0.0

        self.test('wrapPI2(PI_2)', wrapPI2(PI_2), PI_2, nl=1)
        self.test('wrapPI2(PI)',   wrapPI2(PI),     PI)  # XXX
        self.test('wrapPI2(PI2)',  wrapPI2(PI2),   0.0)
        self.test('wrapPI2(-PI_2)',    wrapPI2(-PI_2), PI_2+PI)
        self.test('wrapPI2(-PI)',      wrapPI2(-PI),        PI)  # XXX
        self.test('wrapPI2(-PI2)', abs(wrapPI2(-PI2)),     0.0)  # -0.0

        self.test('unroll180', fstr(unroll180(-90, 110, wrap=True)), '-160.0, -250.0', nl=1)
        self.test('unroll180', fstr(unroll180(-90, 110, wrap=False)), '200.0, 110.0')

        self.test('unroll180', fstr(unroll180(-90, 830, wrap=True)), '-160.0, -250.0', nl=1)
        self.test('unroll180', fstr(unroll180(-90, 830, wrap=False)), '920.0, 830.0')

        self.test('unroll180', fstr(unroll180(-110, 90, wrap=True)), '-160.0, -270.0', nl=1)
        self.test('unroll180', fstr(unroll180(-110, 90, wrap=False)), '200.0, 90.0')

        self.test('unroll180', fstr(unroll180(-830, 90, wrap=True)), '-160.0, -990.0', nl=1)
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

        if sr:  # coverage
            c, _ = cot_(r, r)  # PYCHOK .next() or __next__()
            self.test('cot_ ', c, cr / sr, prec=12, nl=1)
        if sd:  # coverage
            c, _ = cotd_(a, a)  # PYCHOK .next() or __next__()
            self.test('cotd_', c, cd / sd, prec=12, nl=1)
        EPS_ = EPS * 8
        self.test('sincos2',  e, EPS_, known=e < EPS_)
        self.test('sincos2d', d, EPS_, known=d < EPS_)
        if Math:
            self.test('atand',    t, EPS,  known=t < EPS)
            self.test('sincosd ', g, EPS_, known=g < EPS_)
            self.test('sincos*d', f, EPS_, known=f < EPS_)

        # <https://www.CivilGeo.com/when-a-foot-isnt-really-a-foot/>
        self.test('iFt2m', ft2m( 614963.91), 187441, fmt='%.0f', nl=1)
        self.test('iFt2m', ft2m(2483759.84), 757050, fmt='%.0f')
        self.test('sFt2m', ft2m( 614962.68, usurvey=True), 187441, fmt='%.0f')
        self.test('sFt2m', ft2m(2483754.87, usurvey=True), 757050, fmt='%.0f')

        self.test('m2iFt', m2ft(187441),  614963.91, prec=2, nl=1)
        self.test('m2iFt', m2ft(757050), 2483759.84, prec=2)
        self.test('m2sFt', m2ft(187441, usurvey=True),  614962.68, prec=2)
        self.test('m2sFt', m2ft(757050, usurvey=True), 2483754.88, prec=2, nt=1)

        for f, m in ((acre2ha,     '0.404686'), (acre2m2, '4046.856422'),
                     (chain2m,    '20.116800'), (fathom2m,   '1.828800'),
                     (furlong2m, '201.168000'), (toise2m,    '1.949044'),
                     (yard2m,      '0.914400'), (m2chain,    '0.049710'),
                     (m2fathom,    '0.546807'), (m2furlong,  '0.004971'),
                     (m2toise,     '0.513072'), (m2yard,     '1.093613')):
            self.test(f.__name__, f(1), m, prec=6)

        self.test('degrees2m', fstr(degrees2m(90), prec=4),        '10007557.1761', nl=1)
        self.test('degrees2m', fstr(degrees2m(90, lat=30), prec=4), '8666798.7443')
        self.test('m2degrees', fstr(m2degrees(degrees2m(90)), prec=1),   '90.0')

        self.test('degrees2m', fstr(degrees2m(180), prec=4),          '20015114.3522', nl=1)
        self.test('degrees2m', fstr(degrees2m(180, lat=3-0), prec=4), '19987684.3336')
        self.test('m2degrees', fstr(m2degrees(degrees2m(180)), prec=1),    '180.0', nt=1)

        t = 'm2degrees2m(%s, lat=%s)'
        for a in range(0, 90, 7):
            d = m2degrees(degrees2m(45, lat=a), lat=a)
            self.test(t % (45, a), d, '45.00', prec=2)

        self.test('isPoints2', isPoints2(None), False, nl=1)

        try:  # coverage
            self.test('tan_2_semi', tan_2(PI, PI=1), ValueError.__name__, nl=1)
        except ValueError as x:
            t = str(x)
            t = t[:t.find('(')+9] + t[t.find(')'):]
            self.test('tan_2_semi', t, 'PI[1] edge (3.141592): semi-circular', nl=1)

        def _fin(x):
            return (NEG0 if x < 0 else 0.0) if abs(x) < 3e-16 else (
                   (NINF if x < 0 else INF) if abs(x) > 5e+15 else x)

        f = sincostan3.__name__ + '(%+.4f)'
        for a in (0, NEG0, PI_2, -PI_2, PI, -PI, PI3_2, -PI_2, PI2, -PI2):
            t = map1(_fin, sin(a), cos(a), tan(a))
            self.test(f % (a,), sincostan3(a), t, known=True)  # =a in (PI3_2, PI2)


if __name__ == '__main__':

    from pygeodesy import utily  # private

    t = Tests(__file__, __version__, utily)
    t.testUtily()
    t.results()
    t.exit()
