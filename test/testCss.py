
# -*- coding: utf-8 -*-

# Test LCC functions and methods.

__all__ = ('Tests',)
__version__ = '19.04.12'

from base import TestsBase, geographiclib

from pygeodesy import CassiniSoldner, fStr, haversine, hypot, toCss


class Tests(TestsBase):

    def testCss(self, *LatLons):

        P = CassiniSoldner(48 + 50/60.0, 2 + 20/60.0, name='Paris')
        self.test(repr(P), P, P)

        f = P.forward(50.9, 1.8)  # Calais
        self.test('forward', fStr(f, prec=6), '-37518.854545, 230003.561828')
        r = P.reverse(*f)
        self.test('reverse', fStr(r, prec=6), '50.9, 1.8')
        f = P.forward4(50.9, 1.8)  # Calais
        self.test('forward4', fStr(f, prec=6), '-37518.854545, 230003.561828, 89.586104, 0.999983')

        r = P.reverse(-38e3, 230e3)
        self.test('reverse', fStr(r, prec=6), '50.899937, 1.793161')
        f = P.forward(*r)
        self.test('forward', fStr(f, prec=6), '-38000.0, 230000.0')
        r = P.reverse4(-38e3, 230e3)
        self.test('reverse4', fStr(r, prec=6), '50.899937, 1.793161, 89.580797, 0.999982')

        for LL in LatLons:
            r = P.reverse(-38e3, 230e3, LatLon=LL)
            self.test('reverse', repr(r), 'LatLon(50°53′59.77″N, 001°47′35.38″E)')

        g = CassiniSoldner(51.4934, 0.0098, name='Greenwich')
        self.test(repr(g), g, g)
        f = g.forward(48 + 50/60.0, 2 + 20/60.0)  # Paris
        self.test('forward', fStr(f, prec=6), '170557.151692, -293280.6051')
        r = g.reverse(*f)
        self.test('reverse', fStr(r, prec=6), '48.833333, 2.333333')

        h = hypot(*f)  # easting + norting ~= distance
        d = haversine(*(g.latlon0 + r))
        self.test('hypot', h, d, fmt='%.3f', known=abs(d - h) < 1000)

        c = toCss(LL(50.9, 1.8, height=1), cs0=P, name='Calais')
        self.test('toCss', c, '-37518.854545 230003.561828 +1.00m')
        self.test('toCss', c.toStr2(C=True), "[E:-37518.854545, N:230003.561828, H:+1.00m, name:'Calais', C:CassiniSoldner(48.833333, 2.333333, name='Paris')]")
        for a, f, x in (('easting',  '%.6f', '-37518.854545'),
                        ('northing', '%.6f', '230003.561828'),
                      # ('latlon',   '%r',   '(50.9, 1.8)'),  # Python 2.6 Ubuntu (50.899999999999999, 1.8)
                        ('height',   '%.1f', '1.0'),
                        ('azi',      '%.9f', '89.586103815'),
                        ('rk',       '%.9f', '0.999982722'),
                        ('name',     '%s',   'Calais'),
                        ('cs0',      '%s',   '48.833333 2.333333')):
            v = getattr(c, a)
            self.test('Css.'+a, v, x, fmt=f)
        r = c.toLatLon(LatLon=LL)
        self.test('Css.'+'toLatLon', repr(r), 'LatLon(50°54′00.0″N, 001°48′00.0″E, +1.00m)')
        self.test('Css.'+'toLatLon.height', r.height, '1.0')
        self.test('Css.'+'toLatLon.name', r.name, 'Calais')
        self.test('Css.'+'toLatLon.datum.name', r.datum.name, 'WGS84')


if __name__ == '__main__':

    from pygeodesy import css, ellipsoidalKarney, \
                               ellipsoidalNvector, \
                               ellipsoidalVincenty

    t = Tests(__file__, __version__, css)

    if geographiclib:
        t.testCss(ellipsoidalKarney.LatLon,
                  ellipsoidalNvector.LatLon,
                  ellipsoidalVincenty.LatLon)
    else:
        t.skip('no geographiclib', n=14)

    t.results()
    t.exit()
