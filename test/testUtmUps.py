
# -*- coding: utf-8 -*-

# Test UTM functions and methods.

__all__ = ('Tests',)
__version__ = '22.01.03'

from base import TestsBase

from pygeodesy import F_DMS, parseUTMUPS5, toUps8, toUtmUps8, \
                      utmups, UtmUps, utmupsValidateOK


class Tests(TestsBase):

    def testUtmUps(self, LL):
        OK = True  # deprecated, was ok='OK', OK='OK'

        u = UtmUps(0, 'N', 448251, 5411932.0001, falsed=False)
        self.test('UtmUps', u.toStr(4), '00 N 448251.0 5411932.0001')
        self.test('UtmUps.ValidateOK', utmupsValidateOK(u), OK)

        u = UtmUps(0, 'N', 448251.795, 5411932.678, falsed=False)
        self.test('UtmUps', u, '00 N 448252 5411933')
        self.test('UtmUps', u.toStr(prec=3), '00 N 448251.795 5411932.678')
        self.test('UtmUps', u.toStr(prec=1, B=True, cs=True), '00Z N 448251.8 5411932.7 n/a n/a')
        self.test('UtmUps.ValidateOK', utmupsValidateOK(u), OK)

        ll = u.toLatLon(LL, unfalse=False)
        self.test('UtmUps.toLatLon', ll, '43.684097°N, 175.265195°E')
        self.test('UtmUps.toLatLon', ll.toStr(form=F_DMS), '43°41′02.75″N, 175°15′54.7″E')

        u = ll.toUtmUps()
        self.test('LL.toUtmUps', u, '60 N 360177 4838249')
        self.test('LL.toUtmUps', u.toStr(prec=3), '60 N 360176.686 4838249.416')
        self.test('LL.toUtmUps', u.toRepr(B=True, cs=True), '[Z:60T, H:N, E:360177, N:4838249, C:-1.19839167°, S:0.99984048]')
#       self.test('LL.toUtmUps', u.toStr2(B=True, cs=True), '[Z:60T, H:N, E:360177, N:4838249, C:-1.19839167°, S:0.99984048]')
        self.test('LL.toUtmUps.ValidateOK', utmupsValidateOK(u), OK)

        # TM8358-2 pg 3-7 ID 1
        u = toUtmUps8('84 17 14.042N', '132 14 52.761W')  # -132.247988889
        self.test('toUtmUps8ID1', u.toStr(prec=2, B=True, cs=True), '00Y N 1530125.78 2426773.6 -132.24798917° 0.99647445')
        self.test('toUtmUps8ID1.ValidateOK', utmupsValidateOK(u), OK)
        u = toUtmUps8('84 17 14.042N', '132 14 52.761W', Utm=None, Ups=None)  # coverage
        self.test('toUtmUps8ID1.ValidateOK', utmupsValidateOK(u), OK)

        # TM8358-2 pg 3-7 ID 2
        u = toUtmUps8('73N', '44E')  # Karney 38n 467368 8100752 (73°00'00.0"N 044°00'00.0"E) MGRS 38XMG6736700752
        self.test('toUtmUps8ID2', u.toStr(prec=2, cs=True), '38 N 3320416.75 632668.43 +44.0° 1.01619505', known=True)
        self.test('toUtmUps8ID2.ValidateOK', utmupsValidateOK(u), OK)
        u = toUtmUps8('73N', '44E', Utm=None, Ups=None)  # coverage
        self.test('toUtmUps8ID2.ValidateOK', utmupsValidateOK(u), OK)
        u = toUps8('73N', '44E', strict=False)  # allow lat outside UPS range
        self.test('toUtmUps8ID2', u.toStr(prec=2, cs=True), '00 N 3320416.75 632668.43 +44.0° 1.01619505')
        self.test('toUtmUps8ID2.ValidateOK', utmupsValidateOK(u), OK, known=True)

        # TM8358-2 pg 3-7 ID 3
        u = toUtmUps8('87 17 14.4S', '132 14 52.303E')  # -132.247861111
        self.test('toUtmUps8ID3', u.toStr(prec=2, B=True, cs=True), '00B S 2222979.47 1797474.9 -132.24786194° 0.99455723')
        self.test('toUtmUps8ID3.ValidateOK', utmupsValidateOK(u), OK)
        u = toUtmUps8('87 17 14.4S', '132 14 52.303E', Utm=None, Ups=None)  # coverage
        self.test('toUtmUps8ID3.ValidateOK', utmupsValidateOK(u), OK)

        # TM8358-2 pg 3-7 ID 4
        u = UtmUps(0, 'N', '1530125.78', '2426773.6')
        self.test('UtmUps.toLatLonID4.ValidateOK', utmupsValidateOK(u), OK)
        ll = u.toLatLon(LL)
        self.test('UtmUps.toLatLonID4', ll.toStr(form=F_DMS), '84°17′14.04″N, 132°14′52.76″W')
        self.test('UtmUps.toLatLonID4', ll, '84.287234°N, 132.247989°W')

        # TM8358-2 pg 3-7 ID 5
        u = UtmUps(0, 'N', '3320416.75', '632668.43')
        self.test('UtmUps.ValidateOK', utmupsValidateOK(u), OK, known=True)
        ll = u.toLatLon(LL)
        self.test('UtmUps.toLatLonID5', ll.toStr(form=F_DMS), '73°00′00.0″N, 044°00′00.0″E')  # 72°59′60.0″N, ...
        self.test('UtmUps.toLatLonID5', ll, '73.0°N, 044.0°E')

        # TM8358-2 pg 3-7 ID 6
        u = UtmUps(0, 'S', '2222979.47', '1797474.9')
        self.test('UtmUps.ValidateOK', utmupsValidateOK(u), OK)
        ll = u.toLatLon(LL)
        self.test('UtmUps.toLatLonID6', ll.toStr(form=F_DMS), '87°17′14.4″S, 132°14′52.3″E')
        self.test('UtmUps.toLatLonID6', ll, '87.287333°S, 132.247861°E')

        # <https://GeographicLib.SourceForge.io/cgi-bin/GeoConvert>
        ll = LL(61.2, -149.9, name='Anchorage')
        self.test('latlon1', ll, ll)
        u = toUtmUps8(ll)  # 06n 344174 6788521 (61°12'00.0"N 149°54'00.0"W) MGRS 06VUN4417388521
        self.test('toUtmUps8', u, '06 N 344174 6788521')
        self.test('toUtmUps8', u.toStr(prec=6, B=True, cs=True), '06V N 344173.864114 6788521.418164 -2.54179531° 0.99989751')
        self.test('toUtmUps8.ValidateOK', utmupsValidateOK(u), OK)
        self.test('toMgrs', u.toMgrs(), '06V UN 44173 88521')

        ll = LL(83.627, -32.664, name='Cape Morris Jesup, Greenland')
        self.test('latlon2', repr(ll), 'LatLon(83°37′37.2″N, 032°39′50.4″W)')
        u = toUtmUps8(ll)  # 25n 504164 9286466 (83°37'37.2"N 032°39'50.4"W) MGRS 25XEN0416386465
        self.test('toUtmUps8', u, '25 N 504164 9286466')
        self.test('toUtmUps8', repr(u), '[Z:25X, H:N, E:504164, N:9286466]')
        self.test('toUtmUps8', u.toStr(prec=6, B=True, cs=True), '25X N 504163.899383 9286465.664902 +20.03542083′ 0.99960021')
        self.test('toUtmUps8.ValidateOK', utmupsValidateOK(u), OK)
        self.test('toMgrs', u.toMgrs(), '25X EN 04163 86465')

        ll = LL(33.33, 44.44)
        self.test('latlon3', ll, ll)
        u = toUtmUps8(ll)  # 38n 447882 3688012 (33°19'48.0"N 044°26'24.0"E) MGRS 38SMB4788288011
        self.test('toUtmUps8', u, '38 N 447882 3688012')
        self.test('toUtmUps8', u.toStr(prec=6, B=True, cs=True), '38S N 447882.413169 3688011.692733 -18.46228466′ 0.99963349')
        self.test('toUtmUps8.ValidateOK', utmupsValidateOK(u), OK)
        self.test('toMgrs', u.toMgrs(), '38S MB 47882 88011')

        ll = LL(-79, -79)
        self.test('latlon4', ll, ll)
        u = toUtmUps8(ll)  # 17s 542594 1229296 (79°00'00.0"S 079°00'00.0"W) MGRS 17CNN4259429296
        self.test('toUtmUps8', u, '17 S 542594 1229296')
        self.test('toUtmUps8', u.toStr(prec=6, B=True, cs=True), '17C S 542594.134555 1229296.157301 -1.96328341° 0.99962217')
        self.test('toUtmUps8.ValidateOK', utmupsValidateOK(u), OK)
        self.test('toMgrs', u.toMgrs(), '17C NN 42594 29296')

        ll = LL(84, 84)
        self.test('latlon5', ll, ll)
        u = toUtmUps8(ll)  # n 2663075 1930308 (84°00'00.0"N 084°00'00.0"E) MGRS ZJG6307530307
        self.test('toUtmUps8', u, '00 N 2663075 1930308')
        self.test('toUtmUps8', u.toStr(prec=6, B=True, cs=True), '00Z N 2663075.299562 1930307.977716 +84.0° 0.99673')
        self.test('toUtmUps8.ValidateOK', utmupsValidateOK(u), OK)
        # self.test('toMgrs', u.toMgrs(), 'Z JG 63075 30307')

        ll = LL(13.4125, 103.8667)
        self.test('latlon6', ll, ll)
        u = toUtmUps8(ll)  # N 377302.354182663 1483034.77706381 -000.26291348° 0.999786229
        self.test('toUtmUps8', u, '48 N 377302 1483035')
        self.test('toUtmUps8', u.toStr(prec=6, B=True, cs=True), '48P N 377302.354183 1483034.777084 -15.77480856′ 0.99978623')
        self.test('toUtmUps8.ValidateOK', utmupsValidateOK(u), OK)
        self.test('toMgrs', u.toMgrs(), '48P UV 77302 83034')

        ll = LL(-13.4125, -103.8667)
        self.test('latlon7', ll, ll)
        u = ll.toUtmUps()  # S 622697.645817337 8516965.22293619 -000.26291348° 0.999786229
        self.test('LL.toUtmUps', u, '13 S 622698 8516965')
        self.test('LL.toUtmUps', u.toStr(prec=6, B=True, cs=True), '13L S 622697.645817 8516965.222916 -15.77480856′ 0.99978623')
        self.test('LL.toUtmUps.ValidateOK', utmupsValidateOK(u), OK)
        self.test('toMgrs', u.toMgrs(), '13L FF 22697 16965')

        ll = LL(43.684097, 175.265195)
        self.test('latlon8', ll, ll)
        u = ll.toUtmUps()
        self.test('LL.toUtmUps', u, '60 N 360177 4838249')
        self.test('LL.toUtmUps', u.toStr(prec=3), '60 N 360176.691 4838249.422')
        self.test('LL.toUtmUps', u.toRepr(B=True, cs=True), '[Z:60T, H:N, E:360177, N:4838249, C:-1.19839163°, S:0.99984048]')
#       self.test('LL.toUtmUps', u.toStr2(B=True, cs=True), '[Z:60T, H:N, E:360177, N:4838249, C:-1.19839163°, S:0.99984048]')
        self.test('LL.toUtmUps.ValidateOK', utmupsValidateOK(u), OK)
        self.test('toMgrs', u.toMgrs(), '60T UP 60176 38249')
        ll = u.toLatLon(LL)
        self.test('UtmUps.toLatLon', ll, '43.684097°N, 175.265195°E')
        self.test('UtmUps.toLatLon', ll.toStr(form=F_DMS), '43°41′02.75″N, 175°15′54.7″E')

        ll = LL(41.321801, -74.801413)  # Milford, PA
        self.test('latlon9', ll, ll)
        u = ll.toUtmUps()
        self.test('LL.toUtmUps', u, '18 N 516620 4574500')
        self.test('LL.toUtmUps', u.toRepr(B=True, cs=True), '[Z:18T, H:N, E:516620, N:4574500, C:+7.86748851′, S:0.9996034]')
#        self.test('LL.toUtmUps', u.toStr2(B=True, cs=True), '[Z:18T, H:N, E:516620, N:4574500, C:+7.86748851′, S:0.9996034]')
        self.test('LL.toUtmUps.ValidateOK', utmupsValidateOK(u), OK)
        self.test('toMgrs', u.toMgrs(), '18T WL 16619 74500')
        ll = u.toLatLon(LL)
        self.test('UtmUps.toLatLon', ll, '41.321801°N, 074.801413°W')
        self.test('UtmUps.toLatLon', ll.toStr(form=F_DMS), '41°19′18.48″N, 074°48′05.09″W')
        u = parseUTMUPS5(str(u))
        self.test('parseUTMUPS5', u, '18 N 516620 4574500')
        self.test('parseUTMUPS5.ValidateOK', utmupsValidateOK(u), OK)
        ll = u.toLatLon(LL)
        self.test('UtmUps.toLatLon', ll, '41.321801°N, 074.801413°W')

        # courtesy of U{sumnamazu<https://GitHub.com/mrJean1/PyGeodesy/issues/26>}
        u = UtmUps(0, 'S', 321441.0425108216, 5810117.133231169)  # falsed=False
        self.test('UtmUps', u.toStr(B=True), '00A S 321441 5810117')
        self.test('UtmUps.ValidateOK', utmupsValidateOK(u), OK, known=True)
        ll = u.toLatLon(LL)
        self.test('UtmUps.toLatLon', ll, '53.713776°S, 023.77604°W')
        self.test('UtmUps.toLatLon', ll.toStr(form=F_DMS), '53°42′49.59″S, 023°46′33.74″W')
        u = ll.toUtmUps()  # 27s 316807 4044745 (53°42'49.6"S 023°46'33.7"W) MGRS 27FUA1680744744
        self.test('LL.toUtmUps', u, '27 S 316807 4044745')
        self.test('LL.toUtmUps.ValidateOK', utmupsValidateOK(u), OK)
        self.test('LL.toUtmUps', u.toStr(prec=3), '27 S 316807.326 4044744.532')
        self.test('LL.toUtmUps', u.toRepr(B=True, cs=True), '[Z:27F, H:S, E:316807, N:4044745, C:+2.23830171°, S:1.00001184]')
#       self.test('LL.toUtmUps', u.toStr2(B=True, cs=True), '[Z:27F, H:S, E:316807, N:4044745, C:+2.23830171°, S:1.00001184]')

        u = UtmUps(0, 'N', 400000, 5000000, falsed=False)
        self.test('UtmUps', u.toStr(B=True), '00Z N 400000 5000000')
        self.test('UtmUps.ValidateOK', utmupsValidateOK(u), OK)

        u = parseUTMUPS5('31X N 446000,8436100', Utm=None, Ups=None)  # Svalbard
        self.test('parseUTMUPS5', u, "(31, 'N', 446000.0, 8436100.0, 'X')")
        u = parseUTMUPS5('00A S 506346 1057743', Utm=None, Ups=None)
        self.test('parseUTMUPS5', u, "(0, 'S', 506346.0, 1057743.0, 'A')")


if __name__ == '__main__':

    from pygeodesy import ellipsoidalVincenty

    t = Tests(__file__, __version__, utmups)
    t.testUtmUps(ellipsoidalVincenty.LatLon)
    t.results()
    t.exit()
