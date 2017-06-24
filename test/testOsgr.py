
# -*- coding: utf-8 -*-

# Test OSGR functions and methods.

__all__ = ('Tests',)
__version__ = '17.06.23'

from base import TestsBase

from pygeodesy import F_DMS, Datums, osgr


class Tests(TestsBase):

    def testOSgr(self, LatLon):

        # check convertDatum and back
        p = LatLon(51.4778, -0.0016, datum=Datums.WGS84)
        self.test('WGS84', p, '51.4778°N, 000.0016°W')
        r = p.convertDatum(Datums.OSGB36)
        r.height = 0
        self.test('OSGB36', r, '51.477284°N, 000.00002°E')
        r = r.convertDatum(Datums.WGS84)
        r.height = 0
        self.test('WGS84', r, '51.4778°N, 000.0016°W')

        g = osgr.Osgr(651409.903, 313177.270)
        self.test('OSgr1', g, 'TG 51409 13177')
        self.test('OSgr1', repr(g), '[G:TG, E:51409, N:13177]')

        p = g.toLatLon(LatLon)
        p.height = 0
        self.test('toLatLon1', p.toStr(F_DMS), '52°39′28.72″N, 001°42′57.74″E', known=True)
        self.test('toLatLon1', p, '52.657977°N, 001.716038°E', known=True)
        r = p.toOsgr()
        self.test('toOsgr1', r.toStr(0), '651409.903, 313177.270', known=True)

        p = g.toLatLon(LatLon, datum=Datums.OSGB36)
        self.test('toLatLon2', p.toStr(F_DMS), '52°39′27.25″N, 001°43′04.47″E', known=True)
        self.test('toLatLon2', p, '52.657568°N, 001.717908°E', known=True)
        r = osgr.toOsgr(p)
        self.test('toOsgr2', r.toStr(0), '651409,313177', known=True)

        p = LatLon(52.65798, 1.71605)
        r = osgr.toOsgr(p)  # TG 51409 13177
        self.test('toOsgr3', r, 'TG 51409 13177')

        r = osgr.toOsgr(52.65757, 1.71791, datum=Datums.OSGB36)
        self.test('toOsgr4', r, 'TG 51409 13177')

        g = osgr.parseOSGR('TG 48251 11932')
        self.test('OSGR1', g, 'TG 48251 11932')
        self.test('OSGR1', repr(g), '[G:TG, E:48251, N:11932]')

        g = osgr.parseOSGR('TG51409 13177')
        self.test('OSGR2', g, 'TG 51409 13177')
        self.test('OSGR2', repr(g), '[G:TG, E:51409, N:13177]')

        g = osgr.parseOSGR('TG5140913177')
        self.test('OSGR3', g, 'TG 51409 13177')
        self.test('OSGR3', repr(g), '[G:TG, E:51409, N:13177]')

        g = osgr.parseOSGR('651409,313177')
        self.test('OSGR4', g, 'TG 51409 13177')
        self.test('OSGR4', repr(g), '[G:TG, E:51409, N:13177]')

        self.test('OSGR5', g.toStr(prec=0), '651409,313177')
        self.test('OSGR5', g.toStr2(prec=-3), '[OSGR:651409.000,313177.000]')

        r = osgr.parseOSGR(g.toStr(prec=-3))
        self.test('OSGR6', r.toStr(prec=0), '651409,313177')


if __name__ == '__main__':

    from pygeodesy import ellipsoidalNvector

    t = Tests(__file__, __version__, osgr)
    t.testOSgr(ellipsoidalNvector.LatLon)
    t.results(nl=0)
    t.exit()
