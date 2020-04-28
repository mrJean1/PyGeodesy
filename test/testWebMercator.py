
# -*- coding: utf-8 -*-

# Test Web Mercator classes functions and methods.

__all__ = ('Tests',)
__version__ = '20.04.22'

from math import log, radians, tan
from base import TestsBase

from pygeodesy import F_D, F_DMS, R_M, R_MA, Datums, LatLon_, \
                      fstr, toWm, webmercator, Wm


class Tests(TestsBase):

    def testWebMercator(self, LatLon, LatLonE):
        # <https://AlastairA.WordPress.com/2011/01/23/
        #        the-google-maps-bing-maps-spherical-mercator-projection>
        lat = 52.4827802220782
        w = toWm(lat, -5.625)
        self.test('toWm1', w.toStr(prec=8), '-626172.13571216 6887893.4928338')
        y = R_MA * log(tan(radians((90 + lat) * 0.5)))
        self.test('Wm1.y', y, '6887893.49283380', fmt='%.8f')
        self.testCopy(w)

        w = Wm(448251.795, 5411932.678)
        self.test('Wm2', w, '448251.795 5411932.678')
        self.test('Wm2', w.toStr(prec=0), '448252 5411933')
        self.test('Wm2', w.toStr(prec=1), '448251.8 5411932.7')
        self.testCopy(w)

        ll = w.to2ll(None)  # 2-tuple
        self.test('Wm2.to2ll', fstr(ll, prec=8), '43.65321741, 4.02671439')

        ll = w.toLatLon(LatLon)
        self.test('Wm2.toLatLon', ll, '43.653217°N, 004.026714°E')
        self.test('Wm2.toLatLon', ll.toStr(form=F_DMS), '43°39′11.58″N, 004°01′36.17″E')

        w = ll.toWm()  # 448251.795205746 5411932.67761691
        self.test('toWm1', w, '448251.795 5411932.678')
        self.test('toWm2', w.toStr(prec=0), '448252 5411933')
        self.test('toWm3', w.toRepr(prec=0, radius=True), '[x:448252, y:5411933, radius:6378137]')
        self.test('toWm3', w.toStr2(prec=0, radius=True), '[x:448252, y:5411933, radius:6378137]')

        t = w.copy()
        self.test('copy', t, w)
        t = w.parseWM(w.toStr(prec=3))
        self.test('parseWM', t, w)
        t = w.parseWM(w.toStr(prec=3, radius=True))
        self.test('parseWM', t, w)
        t = w.parseWM(w.toStr(prec=3, radius=0))
        self.test('parseWM', t.toStr2(radius=True), w.toStr2(radius=True))

        ll = LatLon(13.4125, 103.8667)
        w = toWm(ll)
        self.test('toWm4', w.toStr(prec=0), '11562388 1506899')
        self.test('toWm4', w.toStr(prec=6), '11562388.154378 1506899.04498')

        ll = LatLonE(13.4125, 103.8667)
        w = toWm(ll)
        self.test('toWm4E', w.toStr(prec=0), '11562388 1496994')
        self.test('toWm4E', w.toStr(prec=6), '11562388.154378 1496993.698095')

        # <https://www.EPSG.org/Portals/0/373-07-02.pdf> page 42
        ll = LatLon('''24°22'54.433"N''', '''100°20'0"W''')
        w = ll.toWm()
        self.test('toWm5', w.toStr(prec=0), '-11169056 2800000')
        self.test('toWm5', w.toStr(prec=6), '-11169055.576258 2800000.003136')

        ll = LatLonE('''24°22'54.433"N''', '''100°20'00.0"W''')
        w = ll.toWm()
        self.test('toWm5E', w.toStr(prec=0), '-11169056 2782367')
        self.test('toWm5E', w.toStr(prec=6), '-11169055.576258 2782367.05923')

        w = Wm(-11169055.58, 2810000)
        ll = w.toLatLon(LatLon)
        self.test('Wm6.toLatLon', ll, '24.46358°N, 100.333333°W')
        self.test('Wm6.toLatLon', ll.toStr(form=F_DMS), '24°27′48.89″N, 100°20′00.0″W')
        ll = w.toLatLon(LatLonE, datum=Datums.WGS84)
        self.test('Wm6.toLatLonE', ll, '24.299812°N, 100.333333°W')
        self.test('Wm6.toLatLonE', ll.toStr(form=F_DMS), '24°17′59.32″N, 100°20′00.0″W')

        w = Wm(-11169055.58, 2800000)
        ll = w.toLatLon(LatLon)
        self.test('Wm7.toLatLon', ll, '24.381787°N, 100.333333°W')
        self.test('Wm7.toLatLon', ll.toStr(form=F_DMS), '24°22′54.43″N, 100°20′00.0″W')
        ll = w.toLatLon(LatLonE, datum=Datums.WGS84)
        self.test('Wm7.toLatLonE', ll, '24.218566°N, 100.333333°W')
        self.test('Wm7.toLatLonE', ll.toStr(form=F_DMS), '24°13′06.84″N, 100°20′00.0″W')

        # <https://AlastairA.WordPress.com/2011/01/23/
        #        the-google-maps-bing-maps-spherical-mercator-projection>
        w = toWm(51.4085960537841, -0.304339270784791, Wm=None)
        self.test('Wm8.toWm', fstr(w, prec=3), '-33878.893, 6693890.382, 6378137.0')
        w = toWm(51.4085960537841, -0.304339270784791)
        self.test('Wm8.toWm', w.toRepr(), '[x:-33878.893, y:6693890.382]')
        self.test('Wm8.toWm', w.toStr2(), '[x:-33878.893, y:6693890.382]')
        self.test('Wm8.toWm', w.toStr(radius=R_M), '-33878.893 6693890.382 6371008.771')
        self.test('Wm8.toWm.x', w.x, '-33878.893',  fmt='%.3f')
        self.test('Wm8.toWm.y', w.y, '6693890.382', fmt='%.3f')
        self.test('Wm8.toWm.latlon', fstr(w.latlon, prec=6), '51.408596, -0.304339')
        self.test('Wm8.toWm.philam', fstr(w.philam, prec=6),  '0.897249, -0.005312')
        ll = w.toLatLon(LatLon)
        self.test('Wm8.toLatLon', ll.toStr(form=F_D, prec=12),  '51.408596053784°N, 000.304339270785°W')
        self.test('Wm8.toLatLon', ll.toStr(form=F_DMS, prec=6), '51°24′30.945794″N, 000°18′15.621375″W')

        for LL, datum in ((LatLon , Datums.WGS84),
                          (LatLon_, Datums.WGS84),
                          (LatLonE, None),
                          (None,    Datums.WGS84),
                          (None,    None)):
            try:
                ll = w.toLatLon(LL, datum=datum)
            except TypeError:
                ll = 'TypeError'
            self.test('Wm9.toLatLon', ll, 'TypeError')

        # <https://Earth-Info.NGA.mil/GandG/wgs84/web_mercator/
        #         %28U%29%20NGA_SIG_0011_1.0.0_WEBMERC.pdf>
        self._TableA1(LatLon, '', 0,
                      '1118889.97', '2273030.93', '3503549.84',
                      '4865942.28', '6446275.84', '8399737.89')
        self._TableA1(LatLon, '', 1.0 / 3600,
                      '1118921.37', '2273063.83', '3503585.55',
                      '4865982.65', '6446323.95', '8399799.73')
        self._TableA1(LatLonE, 'E', 0,
                      '1111475.10', '2258423.65', '3482189.09',
                      '4838471.40', '6413524.59', '8362698.55')
        self._TableA1(LatLonE, 'E', 1.0 / 3600,
                      '1111506.30', '2258456.36', '3482224.61',
                      '4838511.61', '6413572.57', '8362760.29')

    def _TableA1(self, LL, E, secs, *ns):
        lat = secs
        for n in ns:
            lat += 10
            t = 'toWm(LatLon%s(%.4f, 0)).y' % (E, lat)
            self.test(t, toWm(LL(lat, 0)).y, n, fmt='%.2f')


if __name__ == '__main__':

    from pygeodesy import ellipsoidalVincenty, sphericalTrigonometry

    t = Tests(__file__, __version__, webmercator)
    t.testWebMercator(sphericalTrigonometry.LatLon,
                      ellipsoidalVincenty.LatLon)
    t.results()
    t.exit()
