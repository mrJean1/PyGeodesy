
# -*- coding: utf-8 -*-

# Test L{etm} projection L{ExactTransverseMercator}.

__all__ = ('Tests',)
__version__ = '23.08.30'

from bases import TestsBase

from pygeodesy import etm, ExactTransverseMercator
from pygeodesy.karney import _K_2_4


class Tests(TestsBase):

    def testEtm(self, LatLon):
        self.subtitle(etm, LatLon.__name__)

        e = etm.toEtm8(ellipsoidalVincenty.LatLon(-2, 88), name='test')  # coverage
        t = '45 S -19368242 3386457' if _K_2_4 else '45 S -20297797 5336899'
        self.test('toEtm8', e, t)
        self.test('name', e.name, 'test')
        u = e.toUtm()
        self.test('toUtm', u, t)
        self.test('name', u.name, 'test')
        self.test('toETM5', etm.parseETM5(t), e)

        self.testCopy(e, 'name')

        e = e.parse('31 N 448251, 5411932', name='parse')  # coverage
        t = '31 N 448251 5411932'
        self.test('parse', e, t)
        self.test('name', e.name, 'parse')
        e = e.parse(t.replace(' ', ','))
        self.test('parse', e, t)
        self.test('name', e.name, 'parse')

    def testExactTM(self, extendp):
        self.subtitle(etm, ExactTransverseMercator.__name__)

        xtm = ExactTransverseMercator(extendp=extendp, name='test')
        self.test('name', xtm.name, 'test')  # coverage
        t = xtm.toStr()
        self.test('toStr', t, t)  # coverage

        e, n, g, k = xtm.forward(40.4, -3.7, lon0=-3)  # Madrid , UTM zone 30?
        self.test('easting',  e,  '-59401.799843' if _K_2_4 else  '-59401.921148', prec=6, nl=1)
        self.test('northing', n, '4470036.484128' if _K_2_4 else '4472390.031129', prec=6)
        self.test('gamma',    g,      '-0.453309' if _K_2_4 else      '-0.453697', prec=6)
        self.test('scale',    k,       '0.999643', prec=6)
        lat, lon, g, k = xtm.reverse(e, n, lon0=-3)
        self.test('lat', lat, '40.357615' if _K_2_4 else '40.400000', prec=6)
        self.test('lon', lon, '-3.699485' if _K_2_4 else '-3.700000', prec=6)
        self.test('gamma', g, '-0.452969' if _K_2_4 else '-0.453697', prec=6)
        self.test('scale', k,  '0.999643', prec=6)

        e, n, g, k = xtm.forward(lat, lon, lon0=-3)
        self.test('easting',  e,  '-59395.300630' if _K_2_4 else  '-59401.921148', prec=6)
        self.test('northing', n, '4465335.557639' if _K_2_4 else '4472390.031129', prec=6)
        self.test('gamma',    g,      '-0.452582' if _K_2_4 else      '-0.453697', prec=6)
        self.test('scale',    k,       '0.999643', prec=6)

        e, n, g, k = xtm.forward(40.3, -74.7, lon0=-75)
        self.test('easting',  e,   '25495.451559' if _K_2_4 else   '25495.511523', prec=6, nl=1)
        self.test('northing', n, '4458754.643158' if _K_2_4 else '4461098.320889', prec=6)
        self.test('gamma',    g,       '0.193888' if _K_2_4 else       '0.194038', prec=6)
        self.test('scale',    k,       '0.999608', prec=6)

        lat, lon, g, k = xtm.reverse(e, n, lon0=-75)
        self.test('lat', lat,  '40.257789' if _K_2_4 else  '40.300000', prec=6)
        self.test('lon', lon, '-74.700195' if _K_2_4 else '-74.700000', prec=6)
        self.test('gamma', g,   '0.193744' if _K_2_4 else   '0.194038', prec=6)
        self.test('scale', k,      '0.999608', prec=6)

        e, n, g, k = xtm.forward(lat, lon, lon0=-75)
        self.test('easting',  e,   '25494.753036' if _K_2_4 else   '25495.511523', prec=6)
        self.test('northing', n, '4454073.433614' if _K_2_4 else '4461098.320889', prec=6)
        self.test('gamma',    g,       '0.193594' if _K_2_4 else       '0.194038', prec=6)
        self.test('scale',    k,       '0.999608', prec=6)

        # <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1TransverseMercatorExact.html>
        # lat, lon, _, _ = xtm.reverse(22051449.037349, -7131237.022729)
        lat, lon, g, k = xtm.reverse(29735142.378357, 4235043.607933)
        self.test('lat', lat,  '6.40143705' if _K_2_4 else '-2.00000000', prec=8, nl=1)
        self.test('lon', lon, '88.26853332' if _K_2_4 else '88.00000000', prec=8)
        self.test('gamma', g, '79.02693716' if _K_2_4 else '67.63332900', prec=8)
        self.test('scale', k,  '7.30720576' if _K_2_4 else '26.33699547', prec=8)

        if extendp:
            e, n, g, k = xtm.forward(lat, lon)
            self.test('easting',  e, '16883115.869330' if _K_2_4 else '29735142.378357', prec=6)
            self.test('northing', n,  '8853236.377995' if _K_2_4 else  '4235043.607933', prec=6)
            self.test('gamma',    g,       '79.822275' if _K_2_4 else       '67.633329', prec=6)
            self.test('scale',    k,        '3.873352' if _K_2_4 else       '26.336995', prec=6, nt=1)

        else:
            e, n, g, k = xtm.forward(-90, -120)  # coverage
            self.test('easting', abs(e),        '0.000' if _K_2_4 else        '0.000', prec=3, nl=1)
            self.test('northing',    n, '-10003578.604' if _K_2_4 else '-9997964.943', prec=3, known=True)
            self.test('gamma',       g,       '120.000' if _K_2_4 else      '120.000', prec=3)
            self.test('scale',       k,         '1.000' if _K_2_4 else        '1.000', prec=3)

            lat, lon, g, k = xtm.reverse(e, n)
            self.test('lat', lat,  '-89.899' if _K_2_4 else  '-90.000', prec=3)
            self.test('lon', lon,  '180.000' if _K_2_4 else    '0.000', prec=3, known=True)
            self.test('gamma', g, '-180.000' if _K_2_4 else '-180.000', prec=3)
            self.test('scale', k,    '1.000' if _K_2_4 else    '1.000', prec=3, nt=1)

        self.testCopy(xtm)


if __name__ == '__main__':

    from pygeodesy import ellipsoidalNvector, ellipsoidalVincenty

    t = Tests(__file__, __version__, etm)
    t.testExactTM(True)
    t.testExactTM(False)
    t.testEtm(ellipsoidalNvector.LatLon)
    t.testEtm(ellipsoidalVincenty.LatLon)
    t.results()
    t.exit()
