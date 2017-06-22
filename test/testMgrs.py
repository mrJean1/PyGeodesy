
# -*- coding: utf-8 -*-

# Test MGRS functions and methods.

__all__ = ('Tests',)
__version__ = '17.06.21'

from base import TestsBase

from pygeodesy import mgrs


class Tests(TestsBase):

    def testMgrs(self, LatLon):

        m = mgrs.Mgrs('31U', 'DQ', 48251, 11932)
        self.test('Mgrs1', str(m), '31U DQ 48251 11932')
        self.test('Mgrs1', repr(m), '[Z:31U, G:DQ, E:48251, N:11932]')

        m = mgrs.parseMGRS('31U DQ 48251, 11932')
        self.test('Mgrs2', str(m), '31U DQ 48251 11932')
        self.test('Mgrs2', repr(m), '[Z:31U, G:DQ, E:48251, N:11932]')

        m = mgrs.parseMGRS('31UDQ4825111932')
        self.test('Mgrs3', str(m), '31U DQ 48251 11932')
        self.test('Mgrs3', repr(m), '[Z:31U, G:DQ, E:48251, N:11932]')

        u = m.toUtm()
        self.test('toUtm1', str(u), '31 N 448251 5411932')
        self.test('toUtm1', repr(u), '[Z:31, H:N, E:448251, N:5411932]')

        m = u.toMgrs()
        self.test('toMgrs', str(m), '31U DQ 48251 11932')

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
    t.results(nl=0)
    t.exit()
