
# -*- coding: utf-8 -*-

# Test MGRS functions and methods.

__all__ = ('Tests',)
__version__ = '22.07.22'

from base import startswith, TestsBase

from pygeodesy import mgrs


class Tests(TestsBase):

    def testMgrs(self, LatLon):

        # courtesy of Richard Wright
        m = mgrs.parseMGRS('42SXD0970538646')
        self.test('Mgrs1', str(m), '42S XD 09705 38646')
        self.test('Mgrs1', repr(m), '[Z:42S, G:XD, E:09705, N:38646]')

        # courtesy of Richard Wright
        m = mgrs.parseMGRS('42SXD1970508646')
        self.test('Mgrs2', str(m), '42S XD 19705 08646')
        self.test('Mgrs2', repr(m), '[Z:42S, G:XD, E:19705, N:08646]')

        m = mgrs.parseMGRS('42SXD1938')  # 2 digits means Km
        self.test('Mgrs3', str(m), '42S XD 19000 38000')  # meter
        self.test('Mgrs3', repr(m), '[Z:42S, G:XD, E:19000, N:38000]')

        s = '31U DQ 48251 11932'
        r = '[Z:31U, G:DQ, E:48251, N:11932]'

        m = mgrs.Mgrs('31U', 'DQ', 48251, 11932)
        self.test('Mgrs4', str(m), s)
        self.test('Mgrs4', repr(m), r)

        m = mgrs.parseMGRS('31U DQ 48251, 11932')
        self.test('Mgrs5', str(m), s)
        self.test('Mgrs5', repr(m), r)

        m = mgrs.parseMGRS('31UDQ4825111932')
        self.test('Mgrs6', str(m), s)
        self.test('Mgrs6', repr(m), r)

        m = mgrs.parseMGRS('31UDQ 4825111932')  # coverage
        self.test('Mgrs7', str(m), s)
        self.test('Mgrs7', repr(m), r)

        # courtesy of U{abubelinha<https://github.com/mrJean1/PyGeodesy/issues/54>}
        m = mgrs.parseMGRS('31TDF3182')
        self.test('Mgrs8', m.toUtm(), '31 N 431000 4582000')  # center=False
        self.test('Mgrs8', m.toLatLon(center=False), "(41.38657, 2.174726, Datum(name='WGS84', ellipsoid=Ellipsoids.WGS84, transform=Transforms.WGS84), -0.54564, 0.999659)")
        u = m.toUtm(center=True)
        self.test('Mgrs8', u, '31 N 431500 4582500')
        self.test('Mgrs8', u.toMgrs(center=True), '31T DF 31000 82000')
        c = m.toLatLon(LatLon=LatLon)  # center=True
        self.test('Mgrs8', c, '41.391116°N, 002.180649°E')
        self.test('Mgrs8', c.toMgrs(center=True), '31T DF 31000 82000')

        m = mgrs.parseMGRS('31UDQ 48251 11932')  # coverage
        self.test('Mgrs9', str(m), s)
        self.test('Mgrs9', repr(m), r)

        u = m.toUtm()
        self.test('toUtm1', str(u), '31 N 448251 5411932')
        self.test('toUtm1', repr(u), '[Z:31U, H:N, E:448251, N:5411932]')

        m = mgrs.toMgrs(u)
        self.test('toMgrs1', str(m), s)
        self.test('toMgrs1', repr(m), r)

        p = m.parse('31UDQ4825111932')  # coverage
        t = p.toUtm(None)
        self.test('toUtm(None)', t, "(31, 'N', 448251.0, 5411932.0)", known=True)  # OBSOLETE UtmUps4Tuple
        self.test('toUtm(None)', t, "(31, 'N', 448251.0, 5411932.0, 'U')")  # UtmUps5Tuple
        for a, x in (('easting',  48251.0),
                     ('northing', 11932.0),
                     ('EN', 'DQ'),
                     ('en100k', 'DQ'),
                     ('digraph', 'DQ'),
                     ('zone', 31),
                     ('band', 'U'),
                     ('bandLatitude', 48)):
            self.test(a, getattr(p, a), x)

        m = u.toMgrs()
        self.test('toMgrs', str(m), s)
        m = mgrs.toMgrs(u, Mgrs=None)
        self.test('toMgrs(None)', m.classname, mgrs.Mgrs6Tuple.__name__)

        for lat, lon, x in ((60.0,  1.0, '31V CG 88455 53097'),  # southern Norway
                            (60.0,  3.0, '32V JM 65640 66593'),
                            (60.0,  9.0, '32V NM 00000 51411'),
                            (76.0,  1.0, '31X DE 45999 36099'),  # Svalbard
                            (76.0, 13.0, '33X VE 45999 36099'),
                            (76.0, 25.0, '35X ME 45999 36099'),
                            (76.0, 37.0, '37X DE 45999 36099'),
                            (84.0, 42.0, '00Z GC 46127 04524')):  # polar
            p = LatLon(lat, lon)
            m = p.toUtmUps().toMgrs()
            self.test('toUtm(%s).toMgrs' % (p,), m, x)

        m = mgrs.parseMGRS('YUB17770380')  # GeoConvert
        self.test('MgrsY', str(m), '00Y UB 17770 03800', nl=1)
        self.test('MgrsY', repr(m), '[Z:00Y, G:UB, E:17770, N:03800]')
        t = m.toUtmUps(center=True)
        self.test('toUtmUps', str(t), '00 N 1617775 1403805')
        self.test('toUtmUps', repr(t), '[Z:00Y, H:N, E:1617775, N:1403805]')
        t = m.toLatLon()  # always centered
        self.test('toLatLon', str(t), "(83.627518, -32.664231, Datum(name='WGS84', ", known=startswith)
        self.test('toLatLon', repr(t), "LatLonDatum5Tuple(lat=83.627518, lon=-32.664231, ", known=startswith)
        m = LatLon(83.627518, -32.664231).toMgrs()
        self.test('toMgrsX!', str(m), '25X EN 04160 86523')  # not '00Y UB 17770 03800'

        # courtesy of U{CF-FHB-X<https://GitHub.com/mrJean1/PyGeodesy/issues/70>}
        m = mgrs.parseMGRS('BFS7751499182')  # GeoConvert
        self.test('MgrsB', str(m), '00B FS 77514 99182', nl=1)
        self.test('MgrsB', repr(m), '[Z:00B, G:FS, E:77514, N:99182]')
        t = m.toUtmUps(center=False)
        self.test('toUtmUps', str(t), '00 S 2377514 2499182')
        self.test('toUtmUps', repr(t), '[Z:00B, H:S, E:2377514, N:2499182]')
        t = m.toLatLon()  # always centered
        self.test('toLatLon', str(t), "(-84.367192, 37.098959, Datum(name='WGS84', ", known=startswith)
        self.test('toLatLon', repr(t), "LatLonDatum5Tuple(lat=-84.367192, lon=37.098959, ", known=startswith)
        m = LatLon(-84.367192, 37.098959).toMgrs()
        self.test('toMgrsB', str(m), '00B FS 77514 99182')

        # DMA TM 8358.1 Appendix B <https://Apps.DTIC.mil/sti/pdfs/ADA247651.pdf>
        m = mgrs.parseMGRS('45SXT4791')  # GeoConvert
        self.test('MgrsS', str(m), '45S XT 47000 91000', nl=1)
        self.test('MgrsS', repr(m), '[Z:45S, G:XT, E:47000, N:91000]')
        t = m.toUtmUps(center=True)
        self.test('toUtmUps', str(t), '45 N 647500 3791500')
        self.test('toUtmUps', repr(t), '[Z:45S, H:N, E:647500, N:3791500]')
        t = m.toLatLon()  # always centered
        self.test('toLatLon', str(t), "(34.254177, 88.601932, Datum(name='WGS84', ", known=startswith)
        self.test('toLatLon', repr(t), "LatLonDatum5Tuple(lat=34.254177, lon=88.601932, ", known=startswith)
        m = LatLon(34.254177, 88.601932).toMgrs()
        self.test('toMgrsS', str(m), '45S XT 47499 91499')

        # DMA TM 8358.1 Appendix B <https://Apps.DTIC.mil/sti/pdfs/ADA247651.pdf>
        m = mgrs.parseMGRS('YXK3543')  # GeoConvert
        self.test('MgrsY', str(m), '00Y XK 35000 43000', nl=1)
        self.test('MgrsY', repr(m), '[Z:00Y, G:XK, E:35000, N:43000]')
        t = m.toUtmUps(center=True)
        self.test('toUtmUps', str(t), '00 N 1735500 2243500')
        self.test('toUtmUps', repr(t), '[Z:00Y, H:N, E:1735500, N:2243500]')
        t = m.toLatLon()  # always centered
        self.test('toLatLon', str(t), "(86.762629, -132.632821, Datum(name='WGS84', ", known=startswith)
        self.test('toLatLon', repr(t), "LatLonDatum5Tuple(lat=86.762629, lon=-132.632821, ", known=startswith)
        m = LatLon(86.762629, -132.632821).toMgrs()
        self.test('toMgrsY', str(m), '00Y XK 35499 43500', nt=1)


if __name__ == '__main__':

    from pygeodesy import ellipsoidalVincenty

    t = Tests(__file__, __version__, mgrs)
    t.testMgrs(ellipsoidalVincenty.LatLon)
    t.results()
    t.exit()
