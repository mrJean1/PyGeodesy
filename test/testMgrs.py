
# -*- coding: utf-8 -*-

# Test MGRS functions and methods.

__all__ = ('Tests',)
__version__ = '20.05.01'

from base import TestsBase

from pygeodesy import mgrs


class Tests(TestsBase):

    def testMgrs(self, LatLon):

        # courtesy Richard Wright
        m = mgrs.parseMGRS('42SXD0970538646')
        self.test('Mgrs1', str(m), '42S XD 09705 38646')
        self.test('Mgrs1', repr(m), '[Z:42S, G:XD, E:09705, N:38646]')

        # courtesy Richard Wright
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

        m = mgrs.parseMGRS('31UDQ 48251 11932')  # coverage
        self.test('Mgrs8', str(m), s)
        self.test('Mgrs8', repr(m), r)

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
                            (76.0, 37.0, '37X DE 45999 36099')):
            p = LatLon(lat, lon)
            m = p.toUtm().toMgrs()
            self.test('toUtm(%s).toMgrs' % (p,), m, x)


if __name__ == '__main__':

    from pygeodesy import ellipsoidalVincenty

    t = Tests(__file__, __version__, mgrs)
    t.testMgrs(ellipsoidalVincenty.LatLon)
    t.results()
    t.exit()
