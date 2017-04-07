
# -*- coding: utf-8 -*-

# Test spherical earth model functions and methods.

__all__ = ('Tests',)
__version__ = '17.03.21'

from tests import Tests as _Tests

from geodesy import F_D, F_DMS, lonDMS


class Tests(_Tests):

    def testSpherical(self, LatLon):
        p = LatLon(52.205, 0.119)
        q = LatLon(48.857, 2.351)
        i = p.intermediateTo(q, 0.25)
        self.test('intermediateTo', i, '51.372084°N, 000.707337°E')

        if hasattr(LatLon, 'intermediateChordTo'):
            i = p.intermediateChordTo(q, 0.25)
            self.test('intermediateChordTo', i, '51.372294°N, 000.707192°E')

        p = LatLon(51.8853, 0.2545)
        q = LatLon(49.0034, 2.5735)
        i = p.intersection(108.55, q, 32.44)
        self.test('intersection', i.toStr(F_D),  '50.907608°N, 004.508575°E')  # 50.9076°N, 004.5086°E  # Trig
        self.test('intersection', i.toStr(F_DMS), '50°54′27.39″N, 004°30′30.87″E')

        REO = LatLon(42.600, -117.866)
        BKE = LatLon(44.840, -117.806)
        i = REO.intersection(51, BKE, 137)
        self.test('intersection', i.toStr(F_D), '43.5719°N, 116.188757°W')  # 43.572°N, 116.189°W
        self.test('intersection', i.toStr(F_DMS), '43°34′18.84″N, 116°11′19.53″W')

        p = LatLon(0, 0)
        self.test('maxLat0',  p.maxLat( 0), '90.0')
        self.test('maxLat1',  p.maxLat( 1), '89.0')
        self.test('maxLat90', p.maxLat(90),  '0.0')

        if hasattr(LatLon, 'crossingParallels'):
            ps = p.crossingParallels(LatLon(60, 30), 30)
            t = ', '.join(map(lonDMS, ps))
            self.test('crossingParallels', t, '009°35′38.65″E, 170°24′21.35″E')

        p = LatLon(51.127, 1.338)
        q = LatLon(50.964, 1.853)
        b = p.rhumbBearingTo(q)
        self.test('rhumbBearingTo', b, '116.722', '%.3f')  # 116.7

        d = p.rhumbDistanceTo(q)
        self.test('rhumbDistanceTo', d, '40307.8', '%.1f')  # 40310 ?

        m = p.rhumbMidpointTo(q)
        self.test('rhumbMidpointo', m, '51.0455°N, 001.595727°E')  # 51.0455°N, 001.5957°E


if __name__ == '__main__':

    from geodesy import sphericalNvector as N
    t = Tests(__file__, __version__, N)
    t.testLatLon(N.LatLon)
    t.testSpherical(N.LatLon)
    t.testVectorial(N.LatLon, N.Nvector, N.sumOf)
    t.results()

    from geodesy import sphericalTrigonometry as T
    t = Tests(__file__, __version__, T)
    t.testLatLon(T.LatLon)
    t.testSpherical(T.LatLon)
    t.results()
    t.exit()
