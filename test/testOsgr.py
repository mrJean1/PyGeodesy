
# -*- coding: utf-8 -*-

# Test L{osgr} module.

__all__ = ('Tests',)
__version__ = '23.03.27'

from bases import GeodSolve, geographiclib, startswith, TestsBase

from pygeodesy import F_D, F_DEG, F_DMS, fstr, Datums, \
                      Osgr, parseOSGR, toOsgr


class Tests(TestsBase):

    def testOSgr(self, module):

        self.subtitle(module, 'OSgr')
        LatLon = module.LatLon

        # check convertDatum and back
        p = LatLon(51.4778, -0.0016)  # datum=Datums.WGS84)
        self.test('WGS84', p, '51.4778°N, 000.0016°W')
        r = p.convertDatum(Datums.OSGB36)
        self.test('OSGB36', r.toStr(form=F_D, m=None), '51.477284°N, 000.00002°E')
        r = r.convertDatum(Datums.WGS84)
        self.test('WGS84', r.toStr(form=F_D, m=None), '51.4778°N, 000.0016°W', known=True)

        g = Osgr(651409.903, 313177.270)
        self.test('OSgr1', g, 'TG 51409 13177')
        self.test('OSgr1', repr(g), '[G:TG, E:51409, N:13177]')
        self.test('iteration', g.iteration, g.iteration)  # coverage

        p = g.toLatLon(LatLon)
        self.test('toLatLon1', p.toStr(F_DMS, m=None), '52°39′28.72″N, 001°42′57.79″E', known=True)
        self.test('toLatLon1', p, '52.657979°N, 001.716052°E', known=True)
        self.test('iteration', p.iteration, g.iteration)  # coverage
        r = p.toOsgr()
        self.test('toOsgr1', r.toStr(prec=-3), '651409.903,313177.270', known=True)  # OLD
        self.test('toOsgr1', r.toStr(prec=3, GD=False),  '651409.903,313177.270', known=True)
        self.test('toOsgr1', r.toStr(prec=3, GD=True), 'TG5140990313177270', known=True)

        p = g.toLatLon(LatLon, datum=Datums.OSGB36)
        self.test('toLatLon2', p.toStr(F_DMS), '52°39′27.25″N, 001°43′04.52″E', known=True)
        self.test('toLatLon2', p, '52.657570°N, 001.717922°E', known=True)
        self.test('iteration', p.iteration, g.iteration)  # coverage
        r = toOsgr(p)
        self.test('toOsgr2', r.toStr(prec=0), '651409,313177', known=True)  # OLD
        self.test('toOsgr2', r.toStr(GD=False), '651409,313177', known=True)

        p = LatLon(52.65798, 1.71605)
        r = p.toOsgr()  # TG 51409 13177
        self.test('toOsgr3', r, 'TG 51409 13177')
        p = r.toLatLon()
        self.test('toLatLon3', p, '(52.65798, 1.71605, ', known=startswith)

        p = LatLon(52.65757, 1.71791, datum=Datums.OSGB36)
        r = toOsgr(p, kTM=True)
        self.test('toOsgr4', r, 'TG 51409 13177')
        t = r.toLatLon(kTM=True, datum=Datums.OSGB36).toStr(prec=-9)
        self.test('toLatLon4', t, '(52.657570000, 1.717910000, ', known=startswith)
        t = r.toLatLon(kTM=False, datum=Datums.OSGB36).toStr(prec=-9)
        self.test('toLatLon4', t, '(52.657569999, 1.717910045, ', known=startswith)
        t = r.toLatLon(kTM=True).toStr(prec=-9)
        self.test('toLatLon4', t, '(52.657978296, 1.716040366, ', known=startswith)
        t = r.toLatLon(kTM=False).toStr(prec=-9)
        self.test('toLatLon4', t, '(52.657978295, 1.716040411, ', known=startswith)

        r = parseOSGR('TG5140900013177000')
        self.test('toOsgr5', r.resolution, 0.001)
        self.test('toOsgr5', r.toStr(GD=True, prec=3), 'TG5140900013177000')
        self.test('toOsgr5', r.toStr(GD=False, prec=3), '651409.000,313177.000')
        t = r.toLatLon(kTM=True).toStr(prec=-9)
        self.test('toLatLon5', t, '(52.657976595, 1.716038422, ', known=startswith)

        g = parseOSGR('TG 48251 11932')
        self.test('OSGR1', g, 'TG 48251 11932', nl=1)
        self.test('OSGR1', repr(g), '[G:TG, E:48251, N:11932]')

        g = parseOSGR('TG51409 13177')
        self.test('OSGR2', g, 'TG 51409 13177')
        self.test('OSGR2', repr(g), '[G:TG, E:51409, N:13177]')

        g = parseOSGR('TG5140913177')
        self.test('OSGR3', g, 'TG 51409 13177')
        self.test('OSGR3', repr(g), '[G:TG, E:51409, N:13177]')

        g = parseOSGR('651409,313177')
        self.test('OSGR4', g, 'TG 51409 13177')
        self.test('OSGR4', repr(g), '[G:TG, E:51409, N:13177]')

        self.test('OSGR5', g.toStr(prec=0), '651409,313177')  # OLD
        self.test('OSGR5', g.toStr(prec=7), 'TG51409001317700')  # OLD

        self.test('OSGR5', g.toStr(GD=False), '651409,313177')
        self.test('OSGR5', g.toStr(GD=False, prec=7), '651409.000000,313177.000000')  # um max
        self.test('OSGR5', g.toStr(GD=False, prec=2), '651409.00,313177.00')

        self.test('OSGR5', g.toRepr(prec=-3), '[OSGR:651409.000,313177.000]')  # OLD

        self.test('OSGR5', g.toRepr(GD=False), '[OSGR:651409,313177]')
        self.test('OSGR5', g.toRepr(GD=False, prec=3), '[OSGR:651409.000,313177.000]')
        self.test('OSGR5', g.toRepr(GD=False, prec=-3), '[OSGR:651,313]')

        r = g.toStr(prec=-3)  # OLD
        self.test('OSGR6', r, '651409.000,313177.000')
        r = parseOSGR(r)  # OLD
        self.test('OSGR6', r.toStr(prec=0), '651409,313177')  # OLD
        r = parseOSGR('651409,313177', Osgr=None)  # coverage  # OLD
        self.test('OSGR6', r.toStr(prec=0), '(651409, 313177)')  # OLD
        r = g.parse('651409, 313177')  # coverage  # OLD
        self.test('OSGR6', r.toStr(prec=0), '651409,313177')  # OLD

        r = parseOSGR(g.toStr(GD=False))
        self.test('OSGR6', r.toStr(GD=False), '651409,313177')
        r = parseOSGR('651409,313177', Osgr=None)  # coverage
        self.test('OSGR6', r.toStr(prec=0), '(651409, 313177)')
        r = g.parse('651409, 313177')  # coverage
        self.test('OSGR6', r.toStr(GD=False, sep=' '), '651409 313177')

        self.test('issue', 38, 38, nl=1)
        # courtesy of U{jaluebbe<https://GitHub.com/mrJean1/PyGeodesy/issues/38>}, with expected
        # results from <https://www.Movable-Type.co.UK/scripts/latlong-os-gridref.html>
        p = LatLon(52, -0.12, datum=Datums.WGS84)
        g = toOsgr(p)
        self.test('toOsgr', g.toRepr(), '[G:TL, E:29158, N:35174]')
        self.test('toOsgr', fstr((g.easting, g.northing), prec=3), '529158.072, 235174.785')
        self.test('toOsgr', g.datum.name, 'OSGB36')
        r = g.toLatLon(LatLon, datum=Datums.OSGB36)
        self.test('toLatLonOSGB36', r.toStr(m=None), '51°59′58.37″N, 000°07′06.14″W')
        self.test('toLatLonOSGB36', r.datum.name, 'OSGB36')
        p = g.toLatLon(LatLon, datum=Datums.WGS84)
        self.test('toLatLonWGS84 ', p.toStr(m=None), '52°00′00.0″N, 000°07′12.0″W')
        self.test('toLatLonWGS84 ', p.datum.name, 'WGS84')

        # courtesy of U{jaluebbe<https://GitHub.com/mrJean1/PyGeodesy/issues/38>}, with expected
        # results from <https://www.Movable-Type.co.UK/scripts/latlong-os-gridref.html>
        g = Osgr(532014, 123971)
        self.test('Osgr', g.toRepr(), '[G:TQ, E:32014, N:23971]', nl=1)
        self.test('Osgr', fstr((g.easting, g.northing), prec=1), '532014.0, 123971.0')
        self.test('Osgr', g.datum.name, 'OSGB36')
        r = g.toLatLon(LatLon, datum=Datums.OSGB36)
        self.test('toLatLonOSGB36', r.toStr(form=F_DEG, m=None), '50.999425N, 000.118417W', known=True)
        self.test('toLatLonOSGB36', r.datum.name, 'OSGB36')
        p = g.toLatLon(LatLon, datum=Datums.WGS84)
        self.test('toLatLonWGS84 ', p.toStr(form=F_DEG, m=None), '50.999995N, 000.120004W', known=True)
        self.test('toLatLonWGS84 ', p.datum.name, 'WGS84')

        g = parseOSGR('TQ3201423971')
        self.test('parseOSGR', g.toRepr(), '[G:TQ, E:32014, N:23971]')
        self.test('parseOSGR', fstr((g.easting, g.northing), prec=1), '532014.0, 123971.0')
        self.test('parseOSGR', g.datum.name, 'OSGB36')
        r = g.toLatLon(LatLon, datum=Datums.OSGB36)
        self.test('toLatLonOSGB36', r.toStr(form=F_DEG, m=None), '50.999425N, 000.118417W', known=True)
        self.test('toLatLonOSGB36', r.datum.name, 'OSGB36')
        p = g.toLatLon(LatLon, datum=Datums.WGS84)
        self.test('toLatLonWGS84', p.toStr(form=F_DEG, m=None), '50.999995N, 000.120004W', known=True)
        self.test('toLatLonWGS84', p.datum.name, 'WGS84')

        # courtesy of U{jaluebbe<https://GitHub.com/mrJean1/PyGeodesy/issues/38>}, with expected
        # results from <https://www.Movable-Type.co.UK/scripts/latlong-os-gridref.html>
        g = toOsgr(LatLon(50.999995, -0.120004, datum=Datums.WGS84))
        self.test('toOsgr', g.toRepr(), '[G:TQ, E:32013, N:23971]')
        self.test('toOsgr', fstr((g.easting, g.northing), prec=3), '532013.969, 123971.046')
        self.test('toOsgr', g.datum.name, 'OSGB36')
        r = g.toLatLon(LatLon, datum=Datums.OSGB36)
        self.test('toLatLonOSGB36', r.toStr(form=F_DEG, m=None), '50.999426N, 000.118417W')
        self.test('toLatLonOSGB36', r.datum.name, 'OSGB36')
        p = g.toLatLon(LatLon, datum=Datums.WGS84)
        self.test('toLatLonWGS84 ', p.toStr(form=F_DEG, m=None), '50.999995N, 000.120004W')
        self.test('toLatLonWGS84 ', p.datum.name, 'WGS84')

        # courtesy of U{jaluebbe<https://GitHub.com/mrJean1/PyGeodesy/issues/38>}, with expected
        # results from <https://www.Movable-Type.co.UK/scripts/latlong-os-gridref.html>
        g = toOsgr(LatLon(50.999995, +0.120004, datum=Datums.WGS84))
        self.test('toOsgr', g.toRepr(), '[G:TQ, E:48853, N:24427]')
        self.test('toOsgr', fstr((g.easting, g.northing), prec=3), '548853.602, 124427.985')
        self.test('toOsgr', g.datum.name, 'OSGB36')
        r = g.toLatLon(LatLon, datum=Datums.OSGB36)
        self.test('toLatLonOSGB36', r.toStr(form=F_DEG, m=None), '50.999422N, 000.121618E')
        self.test('toLatLonOSGB36', r.datum.name, 'OSGB36')
        p = g.toLatLon(LatLon, datum=Datums.WGS84)
        self.test('toLatLonWGS84 ', p.toStr(form=F_DEG, m=None), '50.999995N, 000.120004E')
        self.test('toLatLonWGS84 ', p.datum.name, 'WGS84')

        # courtesy of U{jaluebbe<https://GitHub.com/mrJean1/PyGeodesy/issues/38>}, with expected
        # results from <https://www.Movable-Type.co.UK/scripts/latlong-os-gridref.html>
        p = LatLon(49.926244, -6.297934)
        self.test('LatLon', p, '49.926244°N, 006.297934°W', nl=1)
        self.test('datum', p.datum.name, 'WGS84')
        g = p.toOsgr()
        self.test('datum', g.datum.name, 'OSGB36')
        self.test('toOsgr', g.toRepr(), '[G:SV, E:91645, N:11753]')
        g = Osgr(g.easting, g.northing)
        self.test('datum', g.datum.name, 'OSGB36')
        q = g.toLatLon(LatLon=LatLon)
        t = q.distanceTo(p)
        k = abs(t) < 0.015
        self.test('LatLon', q.toStr(form=F_D, m=None), '49.926244°N, 006.297934°W', known=k)
        self.test('datum', q.datum.name, 'WGS84')
        self.test('distanceTo', t, 0.0104, fmt='%.4f', known=k, nt=1)

        g = p.toOsgr(prec=-2)  # coverage
        self.test('prec=-2', g.toRepr(), '[G:SV, E:91600, N:11700]', nt=1)

        for d in (Datums.WGS84, Datums.OSGB36):
            p = LatLon(52, -0.12, datum=d)
            g = toOsgr(p)
            r = g.toLatLon(LatLon, datum=d)
            t = r.toStr(form=F_DEG, m=None)  # '52.0N, 000.12W'
            self.test('toLatLon', t, t)
            for _ in range(3):
                g = g.copy()
                r = g.toLatLon(LatLon, datum=d)
                self.test('toLatLon', r.toStr(form=F_DEG, m=None), t)


if __name__ == '__main__':

    from pygeodesy import ellipsoidalExact, ellipsoidalKarney, osgr, \
                          ellipsoidalNvector, ellipsoidalVincenty  # sphericalNvector

    t = Tests(__file__, __version__, osgr)
    t.testOSgr(ellipsoidalNvector)
#   t.testOSgr(sphericalNvector)
    t.testOSgr(ellipsoidalVincenty)

    if geographiclib:
        t.testOSgr(ellipsoidalKarney)
    t.testOSgr(ellipsoidalExact)
    if GeodSolve:
        from pygeodesy import ellipsoidalGeodSolve
        t.testOSgr(ellipsoidalGeodSolve)

    t.results()
    t.exit()
