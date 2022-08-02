
# -*- coding: utf-8 -*-

# Test UTM functions and methods.

__all__ = ('Tests',)
__version__ = '22.07.27'

from base import endswith, TestsBase

from pygeodesy import degDMS, F_DMS, parseUTMUPS5, RangeError, strs, \
                      toUps8, toUtmUps8, ups, Ups, UtmUps


class Tests(TestsBase):

    def testUps(self, LL):
        u = Ups(0, 'N', 448251, 5411932.0001)  # NOT Ups
        self.test('Ups', u.toStr(4), '00 N 448251.0 5411932.0001')

        u = Ups(0, 'N', 448251.795, 5411932.678, falsed=False)  # NOT Ups
        self.test('Ups', u, '00 N 448252 5411933')
        self.test('Ups', u.toStr(prec=3), '00 N 448251.795 5411932.678')
        self.test('Ups', u.toStr(prec=1, B=True, cs=True), '00Z N 448251.8 5411932.7 n/a n/a')
        self.test('Ups2', u.toRepr(), '[Z:00, H:N, E:448252, N:5411933]')
#       self.test('Ups2', u.toStr2(), '[Z:00, H:N, E:448252, N:5411933]')

        # u = ll.toUps(falsed=False)  # UTM N 448251.795205746 5411932.67761691
        # self.test('LL.toUps', u, 'N 448252 5411933')
        # self.test('LL.toUps', u.toStr(prec=3), 'N 448251.795 5411932.678')
        # self.test('LL.toUps', u.toStr2(B=True, cs=True), '[Z:00Z P:N E:448252 N:5411933 C:+175.26519494° S:1.17547892]')

        u  = UtmUps(60, 'N', 360176.69112, 4838249.4217, falsed=True)
        ll = u.toLatLon(LL, unfalse=False)  # UTM 48.85820000°N, 002.29450000°E
        self.test('UtmUps.toLatLon', ll, '43.610051°N, 004.46308°E')
        self.test('UtmUps.toLatLon', ll.toStr(form=F_DMS), '43°36′36.18″N, 004°27′47.09″E')

        m = u.toMgrs()
        self.test('UtmUps.toMgrs', m, '60T UP 60176 38249')
        try:
            t = u.toUps()
            self.test('toUps', str(u), RangeError.__name__)
        except Exception as x:
            self.test('toUps', str(x), 'inside UTM range [-79.5, 83.5]', known=endswith)
        t = u.toUtm(u.zone)
        self.test('UtmUps.toUtm', t, '60 N 360177 4838249')

        # TM8358-2 pg 3-7 ID 1
        u = toUps8('84 17 14.042N', '132 14 52.761W')  # -132.247988889
        self.test('toUpsID1', u.toStr(prec=2, cs=True), '00 N 1530125.78 2426773.6 -132.24798917° 0.99647445')

        # TM8358-2 pg 3-7 ID 2
        u = toUtmUps8('73N', '44E')  # Karney 38n 467368 8100752 (73°00'00.0"N 044°00'00.0"E) MGRS 38XMG6736700752
        self.test('toUtmUps8ID2', u.toStr(prec=2, cs=True), '38 N 3320416.75 632668.43 +44.0° 1.01619505', known=True)
        u = toUps8('73N', '44E', strict=False)  # allow lat outside UPS range
        self.test('toUtmUps8ID2', u.toStr(prec=2, cs=True), '00 N 3320416.75 632668.43 +44.0° 1.01619505')

        # TM8358-2 pg 3-7 ID 3
        u = toUps8('87 17 14.4S', '132 14 52.303E', pole='S')  # -132.247861111
        self.test('toUpsID3', u.toStr(prec=2, cs=True), '00 S 2222979.47 1797474.9 -132.24786194° 0.99455723')

        # TM8358-2 pg 3-7 ID 4
        u = Ups(0, 'N', '1530125.78', '2426773.6')
        ll = u.toLatLon(LL)
        self.test('Ups.toLatLonID4', ll.toStr(form=F_DMS), '84°17′14.04″N, 132°14′52.76″W')
        self.test('Ups.toLatLonID4', ll, '84.287234°N, 132.247989°W')

        # TM8358-2 pg 3-7 ID 5
        u = Ups(0, 'N', '3320416.75', '632668.43')
        ll = u.toLatLon(LL)
        self.test('Ups.toLatLonID5', ll.toStr(form=F_DMS), '73°00′00.0″N, 044°00′00.0″E')  # '72°59′60.0″N, ...
        self.test('Ups.toLatLonID5', ll, '73.0°N, 044.0°E')

        # TM8358-2 pg 3-7 ID 6
        u = Ups(0, 'S', '2222979.47', '1797474.9')
        ll = u.toLatLon(LL)
        self.test('Ups.toLatLonID6', ll.toStr(form=F_DMS), '87°17′14.4″S, 132°14′52.3″E')
        self.test('Ups.toLatLonID6', ll, '87.287333°S, 132.247861°E')

        # <https://GeographicLib.SourceForge.io/cgi-bin/GeoConvert>
        ll = LL(84, 84)
        self.test('latlon', ll, ll)
        u = toUps8(ll)  # n 2663075 1930308 (84°00'00.0"N 084°00'00.0"E) MGRS ZJG6307530307
        self.test('toUps', u, '00 N 2663075 1930308')
        self.test('toUps', u.toStr(prec=6, cs=True), '00 N 2663075.299562 1930307.977716 +84.0° 0.99673')
        # self.test('toMgrs5', u.toMgrs(), 'Z JG 63075 30307')

        t = ' '.join(toUps8(ll, Ups=None).toStr(prec=6).split()[:5] + ['...)'])
        self.test('toUps(None)', t, "(0, 'N', 2663075.299562, 1930307.977716, 'Z', ...)", known=True)  # coverage
        self.test('.scale0', u.scale0, '0.994000', fmt='%.6f')
        u.rescale0(84, 1.0)
        self.test('rescale0', u.scale0, '0.997261', fmt='%.6f')

        # <https://Earth-Info.NGA.mil/GandG/publications/NGA_SIG_0012_2_0_0_UTMUPS/NGA.SIG.0012_2.0.0_UTMUPS.pdf>
        # 10.2 Examples of computng {x, y, sigma, gamma}, given {lambda, phi, Z} page 41
        # replaced with Karney's results from <https://GeographicLib.SourceForge.io/cgi-bin/GeoConvert>
        #  8: lat lon 83   90  UTM/UPS 46n 459200.256323 9217519.441609  MGRS 46X DT 5920025632317519441609
        #  9: lat lon 82   91  UTM/UPS 46n 468930.934996 9105366.008486  MGRS 46X DS 6893093499605366008486
        # 10: lat lon 81  179  UTM/UPS 60n 534921.971582 8993806.415149  MGRS 60X WQ 3492197158193806415148
        # 11: lat lon 80  180  UTM/UPS 01n 441867.784867 8883084.955948  MGRS 01X DJ 4186778486783084955948
        # 12: lat lon 40    0  UTM/UPS 31n 243900.352030 4432069.056899  MGRS 31T BE 4390035202932069056898
        # 13: lat lon  3 -179  UTM/UPS 01n 277707.830749  331796.291679  MGRS 01N BD 7770783074931796291678
        # 14: lat lon  2  -90  UTM/UPS 16n 166223.907623  221366.166030  MGRS 16N AH 6622390762321366166030
        # 15: lat lon  1   -1  UTM/UPS 30n 722561.736479  110597.972524  MGRS 30N YG 2256173647810597972523
        # 16: lat lon  0    0  UTM/UPS 31n 166021.443081       0.000000  MGRS 31N AA 6602144308000000000000
        # 17: lat lon -1    1  UTM/UPS 31s 277438.263521 9889402.027476  MGRS 31M BU 7743826352189402027476
        # 18: lat lon -2   90  UTM/UPS 46s 166223.907623 9778633.833970  MGRS 46M AC 6622390762378633833969
        # 19: lat lon -3  179  UTM/UPS 60s 722292.169251 9668203.708321  MGRS 60M YB 2229216925068203708321
        # 20: lat lon -4  180  UTM/UPS 01s 166831.065275 9557263.747314  MGRS 01M AR 6683106527557263747313
        for t in ('1   -0 90 N   2000000.0        2000000.0      0.994      -0',
                  '2 -179 89 N   1998062.320046   2111009.610243 0.994076 -179',
                  '3  -90 88 N   1777930.731071   2000000.0      0.994303  -90',
                  '4   -1 87 N   1994185.827038   1666906.254073 0.994682   -1',
                  '5    0 86 N   2000000.0        1555731.570643 0.995212    0',
                  '6    1 85 N   2009694.068153   1444627.207468 0.995895    1',
                  '7   89 84 N   2666626.157825   1988363.997132 0.996730   89',
               #  '8   90 83 N   2778095.750322   2000000.0      0.997718   90',
                  '8   90 83 N    459200.256323   9217519.441609 0.997718   -2.97767886',
               #  '9   91 82 N   2889442.490749   2015525.276426 0.998860   91',
                  '9   91 82 N    468930.934996   9105366.008486 0.998860   -1.98055172',
               # '10  179 81 N   2017473.190606   3001038.419357 1.000156  179',
                 '10  179 81 N    534921.971582   8993806.415149 1.000156   +1.97539632',
               # '11  180 80 N   2000000.0        3112951.136955 1.001608  180',
                 '11  180 80 N    441867.784867   8883084.955948 1.001608   -2.95450468',
               # '12    0 40 N   2000000.0       -3918313.984953 1.209619    0',
                 '12    0 40 N   243900.35203     4432069.056899 1.0004075  -1.92940969',
               # '13 -179  3 N   1790630.987261  13994742.706481 1.883453 -179',
                 '13 -179  3 N    277707.830749    331796.291679 1.00021172 -0.1047151895',
               # '14  -90  2 N -10206568.118587   2000000.0      1.914973  -90',
                 '14  -90  2 N    166223.907623    221366.16603  1.00097936 -0.104796101',
               # '15   -1  1 N   1783239.204558 -10418217.653909 1.947589   -1',
                 '15   -1  1 N    722561.736479    110597.972524 1.00021322  0.03491928033333334',
               # '16    0  0 N   2000000.0      -10637318.498257 1.981349    0',
                 '16    0  0 N    166021.443081         0.0      1.00098106  0',
               # '17    1 -1 N   2224408.737826 -10856367.979638 2.016305    1',
                 '17    1 -1 S    277438.263521   9889402.027476 1.00021322  0.03491928033333334',
               # '18   90 -2 N  15083269.373905   2000000.0      2.052510   90',
                 '18   90 -2 S    166223.907623   9778633.83397  1.00097936  0.104796101',
               # '19  179 -3 N   2232331.498720  15310262.647286 2.090020  179',
                 '19  179 -3 S    722292.169251   9668203.708321 1.00021172 -0.1047151895',
               # '20  180 -4 N   2000000.0       15545537.944524 2.128897  180',
                 '20  180 -4 S    166831.065275   9557263.747314 1.00097428  0.209463796167',):
            i, lon, lat, p, e, n, s, c = t.split()
            u = toUtmUps8(lat, lon)
            z = '%02d' % (u.zone,)
            x = ' '.join((z, p, e, n, degDMS(float(c), prec=8, pos='+'), s))
            t = u.toStr(prec=6, cs=True)
            if abs(float(s) - u.scale) < 1e-2:
                t = ' '.join(t.split()[:-1] + [s])
            self.test('NGA-10.2-' + i, t, x)

        # <https://Earth-Info.NGA.mil/GandG/publications/NGA_SIG_0012_2_0_0_UTMUPS/NGA.SIG.0012_2.0.0_UTMUPS.pdf>
        # 10.3 Examples of computing {lambda, phi}, given {Z, x, y} page 41
        for t in ('1 S       0       0 -135.0          -64.9164123332',
                  '2 S 1000000       0 -153.4349488229 -70.0552944014',
                  '3 S 2000000       0 -180.0          -72.1263610163',
                  '4 S 3000000       0  153.4349488229 -70.0552944014',
                  '5 S 4000000       0  135.0          -64.9164123332',
                  '6 S       0 1000000 -116.5650511771 -70.0552944014',
                  '7 S 1000000 1000000 -135.0          -77.3120791908',
                  '8 S 2000000 1000000  180.0          -81.0106632645',
                  '9 S 3000000 1000000  135.0          -77.3120791908',
                 '10 S 4000000 1000000  116.5650511771 -70.0552944014',
                 '11 S       0 2000000  -90.0          -72.1263610163',
                 '12 S 1000000 2000000  -90.0          -81.0106632645',
                 '13 S 2000000 2000000    0.0          -90.0',
                 '14 S 3000000 2000000   90.0          -81.0106632645',
                 '15 S 4000000 2000000   90.0          -72.1263610163',
                 '16 S       0 3000000  -63.4349488229 -70.0552944014',
                 '17 S 1000000 3000000  -45.0          -77.3120791908',
                 '18 S 2000000 3000000    0.0          -81.0106632645',
                 '19 S 3000000 3000000   45.0          -77.3120791908',
                 '20 S 4000000 3000000   63.4349488229 -70.0552944014',
                 '21 S       0 4000000  -45.0          -64.9164123332',
                 '22 S 1000000 4000000  -26.5650511771 -70.0552944014',
                 '23 S 2000000 4000000    0.0          -72.1263610163',
                 '24 S 3000000 4000000   26.5650511771 -70.0552944014',
                 '25 S 4000000 4000000   45.0          -64.9164123332'):
            i, p, e, n, lon, lat = t.split()
            x = lat + ' ' + lon
            u = parseUTMUPS5(' '.join(('00', p, e, n)))
            ll = u.toLatLon(LL)
            t = ' '.join(strs(ll.latlon, prec=10))
            self.test('NGA-10.3-' + i, t, x, known=i == '3')

        u = LL(83.6, 0).toUps()  # coverage
        self.test('toUps', str(u.toUps()), '00 N 2000000 1288738')
        self.test('toUtm', str(u.toUtm(2)), '02 N 611555 10703765')
        u = Ups()  # default kwds
        self.test('toUtm', repr(u), '[Z:00Z, H:N, E:2000000, N:2000000]')


if __name__ == '__main__':

    from pygeodesy import ellipsoidalVincenty

    t = Tests(__file__, __version__, ups)
    t.testUps(ellipsoidalVincenty.LatLon)
    t.results()
    t.exit()
