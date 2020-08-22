
# -*- coding: utf-8 -*-

# Test Albers Equal-Area projection.

__all__ = ('Tests',)
__version__ = '20.08.22'

from base import TestsBase  # RandomLatLon

from pygeodesy import albers, AlbersError, AlbersEqualArea, \
                      AlbersEqualArea2, AlbersEqualArea4, \
                      AlbersEqualAreaCylindrical, \
                      AlbersEqualAreaNorth, AlbersEqualAreaSouth, \
                      Datums, fstr, sincos2d

_WGS84 = Datums.WGS84
_NAD27 = Datums.NAD27  # Clarke1866


class Tests(TestsBase):

    def testAlbers2(self):

        A = AlbersEqualArea2(40 + 58/60.0, 39 + 56/60.0, name='Karney_example')  # WGS84
        self.test('name',      A.named,                A.named)
        self.test('datum',     A.datum.name,           _WGS84.name)
        self.test('ellipsoid', A.datum.ellipsoid.name, _WGS84.ellipsoid.name)
        for n, x in (('lat0',        '40.451991337063'),
                     ('scale0',      '0.999959500363'),
                     ('equatoradius', A.datum.ellipsoid.a),
                     ('flattening',   A.datum.ellipsoid.f),
                     ('_sign',       '1.000000000000'),
                     ('_m02',        '0.580681094922'),
                     ('_n0',         '0.648810669236'),
                     ('_txi0',       '0.848822476849')):
            self.test(n, getattr(A, n), x, fmt='%.12f')
        self.test('iteration', A.iteration, 3, known=A.iteration > 0)
        self.test('ispolar', A.ispolar, False)

        for lon0, x in ((0,     '-5675721.76113534, 2516917.91242155, 39.95, -75.17, 311.23285234, 0.99999745'),
                        (-77.5, '199089.12574012, -53115.52801838, 39.95, 2.33, 1.51160641, 0.99999745')):
            f = A.forward(39.95, -75.17, lon0=lon0)
            self.test('forward', fstr(f[:6], prec=8), x, known=round(f.x, 7) == -5675721.7611353)
            r = A.reverse(f.x, f.y, lon0=lon0)
            self.test('reverse', fstr(r[:6], prec=8), x, known=round(r.lon, 2) == -75.17)

        for lon0, x in ((0,     '220000.0, -53000.0, 39.94581132, 2.57463362, 1.67031446, 0.99999808'),
                        (-77.5, '220000.0, -53000.0, 39.94581132, -74.92536638, 1.67031446, 0.99999808')):
            r = A.reverse(220e3, -53e3,  lon0=lon0)
            self.test('reverse', fstr(r[:6], prec=8), x)
            f = A.forward(r.lat, r.lon, lon0=lon0)
            self.test('forward', fstr(f[:6], prec=8), x, known=round(f.lon, 8) == 2.57463362)

        self.subtitle(albers, 'Page292')
        A = AlbersEqualArea2(29.5, 45.5, datum=_NAD27, name='Snyder_p292')
        self.test('name',      A.named,                A.named)
        self.test('datum',     A.datum.name,           _NAD27.name)
        self.test('ellipsoid', A.datum.ellipsoid.name, _NAD27.ellipsoid.name)
        for n, x in (('lat0',        '37.934543880726'),
                     ('scale0',       '0.990309187872'),
                     ('equatoradius',  A.datum.ellipsoid.a),
                     ('flattening',    A.datum.ellipsoid.f),
                     ('_sign',        '1.000000000000'),
                     ('_m02',         '0.623664507732'),
                     ('_n0',          '0.614760830736'),
#                    ('_nrho0', '5037024.736824393272'),
                     ('_txi0',        '0.775925617021')):
            self.test(n, getattr(A, n), x, fmt='%.12f')
        self.test('iteration', A.iteration, 4, known=A.iteration > 0)
        self.test('ispolar', A.ispolar, False)

        for lon0, x in ((0,   '-6105839.22928148, 2214046.74930274, 35.0, -75.0, 314.78223745, 0.99155461'),
                        (-96, '1885472.72581347, -119505.66687766, 35.0, 21.0, 12.66097351, 0.99155461')):
            f = A.forward(35, -75, lon0=lon0)
            self.test('forward', fstr(f[:6], prec=8), x, known=round(f.x, 7) == -6105839.2292815 or
                                                               round(f.y, 7) ==  -119505.6668777)
            r = A.reverse(f.x, f.y, lon0=lon0)
            self.test('reverse', fstr(r[:6], prec=8), x, known=round(r.lon, 1) == -75.0)

        for lon0, x in ((0,   '1885427.7, 1535925.0, 49.40436665, 25.93001383, 15.63329611, 1.01436109'),
                        (-96, '1885427.7, 1535925.0, 49.40436665, -70.06998617, 15.63329611, 1.01436109')):
            r = A.reverse(1885427.7, 1535925, lon0=lon0)
            self.test('reverse', fstr(r[:6], prec=8), x)
            f = A.forward(r.lat, r.lon, lon0=lon0)
            self.test('forward', fstr(f[:6], prec=8), x, known=round(f.lon, 8) == 25.93001383)

        self.subtitle(albers, 'Table15')
        A = AlbersEqualArea2(29.5, 45.5, datum=Datums.NAD27, name='Snyder_p103')
        for lat, k in ((52,   1.02863),
                       (50,   1.01727),
                       (45.5, 1.00000),
                       (45,   0.99869),
                       (40,   0.99097),
                       (35,   0.99155),
                       (30,   0.99893),
                       (29.5, 1.00000),
                       (25,   1.01222),
                       (22,   1.02283)):
            t = A.forward(lat, 0)
            self.test(str(lat) + ' k', t.scale, k, fmt='%.5f')

    def testLats(self):
        self.subtitle(albers, 'Lats')
        s, c = sincos2d(30)
        for lat, A in (( 45, AlbersEqualArea(45)),
                       ( 40, AlbersEqualArea2(40, 40)),
                       ( 30, AlbersEqualArea4( s, c,  s, c)),
                       (-30, AlbersEqualArea4(-s, c, -s, c)),
                       (  0, AlbersEqualAreaCylindrical()),
                       ( 90, AlbersEqualAreaNorth()),
                       (-90, AlbersEqualAreaSouth())):
            self.test(A.named + '.lat0', A.lat0, lat, fmt='%.1f')
            self.test(A.named + '.lat1', A.lat1, lat, fmt='%.1f')
            self.test(A.named + '.lat2', A.lat2, lat, fmt='%.1f')

        try:
            self.test('error', AlbersEqualArea4(s, -c, 0, c), AlbersError.__name__)
        except Exception as x:
            self.test('error', str(x), 'lat1 (150.0): above 90 limit')
        try:
            self.test('error', AlbersEqualArea4(-0.5, c, 0.5, c), AlbersError.__name__)
        except Exception as x:
            self.test('error', str(x), 'slat1 (-0.5) or slat2 (0.5): invalid')


if __name__ == '__main__':

    t = Tests(__file__, __version__, albers)
    t.testAlbers2()
    t.testLats()
    t.results()
    t.exit()
