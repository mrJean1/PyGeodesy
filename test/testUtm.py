
# -*- coding: utf-8 -*-

# Test UTM functions and methods.

__all__ = ('Tests',)
__version__ = '17.06.21'

from base import TestsBase

from pygeodesy import F_DEG, F_DMS, parseUTM, utm


class Tests(TestsBase):

    def testUtm(self, LatLon):
        u = utm.Utm(3, 'N', 448251, 5411932.0001)
        self.test('Utm1', u.toStr(4), '03 N 448251.0 5411932.0001')

        u = utm.Utm(31, 'N', 448251.795, 5411932.678)
        self.test('Utm2', u, '31 N 448252 5411933')
        self.test('Utm2', u.toStr(prec=3), '31 N 448251.795 5411932.678')
        self.test('Utm2', u.toStr(prec=1, cs=True), '31 N 448251.8 5411932.7 n/a n/a')

        ll = u.toLatLon(LatLon)  # 48.85820000°N, 002.29450000°E
        self.test('Utm.toLatLon1', ll, '48.8582°N, 002.2945°E')
        self.test('Utm.toLatLon1', ll.toStr(form=F_DMS),  '48°51′29.52″N, 002°17′40.2″E')

        u = ll.toUtm()  # 31U N 448251.795205746 5411932.67761691
        self.test('toUtm1', u, '31 N 448252 5411933')
        self.test('toUtm1', u.toStr(prec=3), '31 N 448251.795 5411932.678')
        self.test('toUtm2', u.toStr2(cs=True), '[Z:31, H:N, E:448252, N:5411933, C:-000.53131221°, S:0.9996329]')

        ll = LatLon(13.4125, 103.8667)
        u = utm.toUtm(ll)  # 48P N 377302.354182663 1483034.77706381 -000.26291348° 0.999786229
        self.test('toUtm4', u, '48 N 377302 1483035')
        self.test('toUtm5', u.toStr(prec=6, B=True, cs=True), '48P N 377302.354183 1483034.777084 -000.26291348° 0.99978623')

        ll = LatLon(-13.4125, -103.8667)
        u = ll.toUtm()  # 13L S 622697.645817337 8516965.22293619 -000.26291348° 0.999786229
        self.test('toUtm6', u, '13 S 622698 8516965')
        self.test('toUtm7', u.toStr(prec=6, B=True, cs=True), '13L S 622697.645817 8516965.222916 -000.26291348° 0.99978623')

        m = u.toMgrs()
        self.test('toMgrs1', m, '13L FF 22697 16965')

        m = utm.Utm('31U', 'N', 448251, 5411932).toMgrs()
        self.test('toMgrs2', m, '31U DQ 48251 11932')

        u = parseUTM('18 N 516620 4574500')  # Milford, PA
        self.test('Utm8', u, '18 N 516620 4574500')
        ll = u.toLatLon(LatLon)
        self.test('Utm8.toLatLon', ll, '41.321801°N, 074.801413°W')
        self.test('Utm8.toLatLon', ll.toStr(F_DEG), '41.321801N, 074.801413W')

        for lat, lon, x in (( 61.44,      25.4,    '35V N 414668 6812845'),  # 35V N 414668.257431168 6812844.72764648
                            (-47.04,     -73.48,   '18G S 615472 4789270'),  # 18G S 615471.65815765  4789269.76738578
                            ( 40.4,      -74.7,    '18T N 525458 4472198'),  # 18T N 525457.882388688 4472198.04072697
                            ( 44.5,      -88.5,    '16T N 380753 4928503'),  # 16T N 380753.114847639 4928503.38224615
                            ( 50.8694,  -115.6508, '11U N 594937 5636169'),  # 11U N 594936.575444796 5636168.98481247
                            (  0.0,        0.0,    '31N N 166021 0'),        # 31N N 166021.443080537       0
                            (  0.13,      -0.2324, '30N N 808084 14386'),    # 30N N 808084.436750719   14385.7989105346
                            (-45.6456,    23.3545, '34G S 683474 4942631'),  # 34G S 683473.746903862 4942631.26945221
                            (-12.765,    -33.8765, '25L S 404859 8588691'),  # 25L S 404859.139809849 8588691.00770755
                            (-80.5434,  -170.654,  '02A S 506346 1057743'),  # outside 02C? S 506346 1057743
                            (  90.0,     177.0,    '60Z N 500000 9997965'),  # outside
                            ( -90.0,    -177.0,    '01A S 500000 2035'),     # outside
                            (  90.0,       3.0,    '31Z N 500000 9997965'),  # outside
                            (  23.4578, -135.4545, '08Q N 453580 2594273'),  # 08Q N 453580 2594273
                            (  77.345,   156.9876, '57X N 450794 8586116'),  # 57X N 450793.553276976 8586116.22730171
                            ( -89.3454,  -48.9306, '22A S 502639 75073'),    # outside
                            (  60.0,       1.0,    '31V N 388456 6653097'),  # southern Norway
                            (  60.0,       3.0,    '32V N 165640 6666594'),
                            (  60.0,       6.0,    '32V N 332705 6655205'),
                            (  60.0,       9.0,    '32V N 500000 6651411'),
                            (  60.0,      12.0,    '33V N 332705 6655205'),
                            (  76.0,       1.0,    '31X N 446000 8436100'),  # Svalbard
                            (  76.0,       7.0,    '31X N 607943 8438843'),
                            (  76.0,      13.0,    '33X N 446000 8436100'),
                            (  76.0,      19.0,    '33X N 607943 8438843'),
                            (  76.0,      25.0,    '35X N 446000 8436100'),
                            (  76.0,      31.0,    '35X N 607943 8438843'),
                            (  76.0,      37.0,    '37X N 446000 8436100')):
            p = LatLon(lat, lon)
            try:
                u = p.toUtm().toStr(prec=0, B=True)
            except ValueError as e:
                if x[2] in 'ABYZ':
                    x = u = str(e)
            self.test('toUtm(%s)' % (p,), u, x)


if __name__ == '__main__':

    from pygeodesy import ellipsoidalVincenty

    t = Tests(__file__, __version__, utm)
    t.testUtm(ellipsoidalVincenty.LatLon)
    t.results(nl=0)
    t.exit()
