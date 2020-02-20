
# -*- coding: utf-8 -*-

# Test OSGR functions and methods.

__all__ = ('Tests',)
__version__ = '20.02.19'

from base import TestsBase

from pygeodesy import F_DEG, F_DMS, fStr, Datums, osgr


class Tests(TestsBase):

    def testOSgr(self, LatLon):

        m = LatLon.__module__
        self.test('LatLon', m, m)

        # check convertDatum and back
        p = LatLon(51.4778, -0.0016, datum=Datums.WGS84)
        self.test('WGS84', p, '51.4778°N, 000.0016°W')
        r = p.convertDatum(Datums.OSGB36)
        r.height = 0
        self.test('OSGB36', r, '51.477284°N, 000.00002°E')
        r = r.convertDatum(Datums.WGS84)
        r.height = 0
        self.test('WGS84', r, '51.4778°N, 000.0016°W', known=True)

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

        self.test('issue', 38, 38)
        # courtesy jaluebbe <https://GitHub.com/mrJean1/PyGeodesy/issues/38>, with expected
        # results from <https://www.Movable-Type.co.UK/scripts/latlong-os-gridref.html>
        p = LatLon(52, -0.12, datum=Datums.WGS84)
        g = osgr.toOsgr(p)
        self.test('toOsgr', g.toStr2(), '[G:TL, E:29158, N:35174]')
        self.test('toOsgr', fStr((g.easting, g.northing), prec=3), '529158.072, 235174.785')
        self.test('toOsgr', g.datum.name, 'OSGB36')
        r = g.toLatLon(LatLon, datum=Datums.OSGB36)
        self.test('toLatLonOSGB36', r.toStr(m=None), '51°59′58.37″N, 000°07′06.14″W')
        self.test('toLatLonOSGB36', r.datum.name, 'OSGB36')
        p = g.toLatLon(LatLon, datum=Datums.WGS84)
        self.test('toLatLonWGS84 ', p.toStr(m=None), '52°00′00.0″N, 000°07′12.0″W')
        self.test('toLatLonWGS84 ', p.datum.name, 'WGS84')

        # courtesy jaluebbe <https://GitHub.com/mrJean1/PyGeodesy/issues/38>, with expected
        # results from <https://www.Movable-Type.co.UK/scripts/latlong-os-gridref.html>
        g = osgr.Osgr(532014, 123971)
        self.test('Osgr', g.toStr2(), '[G:TQ, E:32014, N:23971]')
        self.test('Osgr', fStr((g.easting, g.northing), prec=1), '532014.0, 123971.0')
        self.test('Osgr', g.datum.name, 'OSGB36')
        r = g.toLatLon(LatLon, datum=Datums.OSGB36)
        self.test('toLatLonOSGB36', r.toStr(form=F_DEG, m=None), '50.999425N, 000.118417W', known=True)
        self.test('toLatLonOSGB36', r.datum.name, 'OSGB36')
        p = g.toLatLon(LatLon, datum=Datums.WGS84)
        self.test('toLatLonWGS84 ', p.toStr(form=F_DEG, m=None), '50.999995N, 000.120004W', known=True)
        self.test('toLatLonWGS84 ', p.datum.name, 'WGS84')

        g = osgr.parseOSGR('TQ3201423971')
        self.test('parseOSGR', g.toStr2(), '[G:TQ, E:32014, N:23971]')
        self.test('parseOSGR', fStr((g.easting, g.northing), prec=1), '532014.0, 123971.0')
        self.test('parseOSGR', g.datum.name, 'OSGB36')
        r = g.toLatLon(LatLon, datum=Datums.OSGB36)
        self.test('toLatLonOSGB36', r.toStr(form=F_DEG, m=None), '50.999425N, 000.118417W', known=True)
        self.test('toLatLonOSGB36', r.datum.name, 'OSGB36')
        p = g.toLatLon(LatLon, datum=Datums.WGS84)
        self.test('toLatLonWGS84', p.toStr(form=F_DEG, m=None), '50.999995N, 000.120004W', known=True)
        self.test('toLatLonWGS84', p.datum.name, 'WGS84')

        # courtesy jaluebbe <https://GitHub.com/mrJean1/PyGeodesy/issues/38>, with expected
        # results from <https://www.Movable-Type.co.UK/scripts/latlong-os-gridref.html>
        g = osgr.toOsgr(LatLon(50.999995, -0.120004, datum=Datums.WGS84))
        self.test('toOsgr', g.toStr2(), '[G:TQ, E:32013, N:23971]')
        self.test('toOsgr', fStr((g.easting, g.northing), prec=3), '532013.969, 123971.046')
        self.test('toOsgr', g.datum.name, 'OSGB36')
        r = g.toLatLon(LatLon, datum=Datums.OSGB36)
        self.test('toLatLonOSGB36', r.toStr(form=F_DEG, m=None), '50.999426N, 000.118417W')
        self.test('toLatLonOSGB36', r.datum.name, 'OSGB36')
        p = g.toLatLon(LatLon, datum=Datums.WGS84)
        self.test('toLatLonWGS84 ', p.toStr(form=F_DEG, m=None), '50.999995N, 000.120004W')
        self.test('toLatLonWGS84 ', p.datum.name, 'WGS84')

        # courtesy jaluebbe <https://GitHub.com/mrJean1/PyGeodesy/issues/38>, with expected
        # results from <https://www.Movable-Type.co.UK/scripts/latlong-os-gridref.html>
        g = osgr.toOsgr(LatLon(50.999995, +0.120004, datum=Datums.WGS84))
        self.test('toOsgr', g.toStr2(), '[G:TQ, E:48853, N:24427]')
        self.test('toOsgr', fStr((g.easting, g.northing), prec=3), '548853.602, 124427.985')
        self.test('toOsgr', g.datum.name, 'OSGB36')
        r = g.toLatLon(LatLon, datum=Datums.OSGB36)
        self.test('toLatLonOSGB36', r.toStr(form=F_DEG, m=None), '50.999422N, 000.121618E')
        self.test('toLatLonOSGB36', r.datum.name, 'OSGB36')
        p = g.toLatLon(LatLon, datum=Datums.WGS84)
        self.test('toLatLonWGS84 ', p.toStr(form=F_DEG, m=None), '50.999995N, 000.120004E')
        self.test('toLatLonWGS84 ', p.datum.name, 'WGS84')

        for d in (Datums.WGS84, Datums.OSGB36):
            p = LatLon(52, -0.12, datum=d)
            g = osgr.toOsgr(p)
            r = g.toLatLon(LatLon, datum=d)
            t = r.toStr(form=F_DEG, m=None)  # '52.0N, 000.12W'
            self.test('toLatLon', r.toStr(form=F_DEG, m=None), t)
            for _ in range(3):
                g = g.copy()
                r = g.toLatLon(LatLon, datum=d)
                self.test('toLatLon', r.toStr(form=F_DEG, m=None), t)


if __name__ == '__main__':

    from pygeodesy import ellipsoidalNvector, ellipsoidalVincenty  # sphericalNvector

    t = Tests(__file__, __version__, osgr)
    t.testOSgr(ellipsoidalNvector.LatLon)
#   t.testOSgr(sphericalNvector.LatLon)
    t.testOSgr(ellipsoidalVincenty.LatLon)
    t.results()
    t.exit()
