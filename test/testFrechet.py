
# -*- coding: utf-8 -*-

# Test the Frechet distances.

__all__ = ('Tests',)
__version__ = '21.05.17'

from base import coverage, GeodSolve, geographiclib, \
                 isPython3, isWindows, TestsBase

from pygeodesy import FrechetError, fstr, LatLon_, randomrangenerator

_rr = randomrangenerator('R')
_ms = [LatLon_(*_ll) for _ll in zip(_rr( -90,  90, 2),  # 90
                                    _rr(-180, 180, 4))]

_ps = [LatLon_(*_ll) for _ll in zip(_rr( -89,  90, 3),  # 60
                                    _rr(-179, 180, 6))]


class Tests(TestsBase):

    def test2(self, Frechet, x, y, **kwds):

        k = kwds.pop('known', False)
        p = kwds.pop('prec',  5)

        def _tstr(t):
            s = list(t[:5])
            s[0] = fstr(t.fd, prec=p)
            return '(%s)' % (', '.join(map(str, s)),)

        f = Frechet(_ms, fraction=1, **kwds)
        n = '%s (%s)' % (f.named, f.units)

        try:
            t = _tstr(f.discrete(_ps))
            self.test(n, t, x, known=k)  # + (t.units,)
        except FrechetError as e:
            self.test(n, str(e), x, known=k)  # + (t.units,)

        try:
            t = _tstr(f.discrete(_ps, fraction=0.5))
            self.test(n, t, y, known=k)  # + (t.units,)
        except FrechetError as e:
            self.test(n, str(e), y, known=k)  # + (t.units,)

        self.testCopy(f)

        f.units = ''  # coverage


if __name__ == '__main__':  # MCCABE 13

    from pygeodesy import fractional, frechet_, FrechetCosineAndoyerLambert, \
                          FrechetCosineForsytheAndoyerLambert, FrechetCosineLaw, \
                          FrechetDegrees, FrechetDistanceTo, FrechetEquirectangular, \
                          FrechetEuclidean, FrechetExact, FrechetFlatLocal, \
                          FrechetFlatPolar, FrechetKarney, FrechetHaversine, \
                          FrechetHubeny, FrechetRadians, FrechetThomas, FrechetVincentys

    def _distance(p1, p2):
        dy, dx = abs(p1.lat - p2.lat), abs(p1.lon - p2.lon)
        if dx < dy:
            dx, dy = dy, dx
        return dx + dy * 0.5

    class FrechetDegrees_(FrechetDegrees):
        '''Custom Frechet.'''
        def distance(self, p1, p2):
            return _distance(p1, p2)

    class FrechetRadians_(FrechetRadians):
        '''Custom Frechet.'''
        def distance(self, p1, p2):
            dy, dx = abs(p1.phi - p2.phi), abs(p1.lam - p2.lam)
            if dx < dy:
                dx, dy = dy, dx
            return dx + dy * 0.5

    t = Tests(__file__, __version__)

    if isPython3:  # XXX different Random?
        t.test2(FrechetDegrees_, (178.5, 74, 56,   19,  5400),
                                 (175.5, 74, 52.5, 29, 10710))

        t.test2(FrechetRadians_, (3.11541, 74, 56,   19,  5400),
                                 (3.06305, 74, 52.5, 29, 10710))

        t.test2(FrechetCosineAndoyerLambert, (2.6319, 0, 0, 149,  5400),
                                             (2.6319, 0, 0, 208, 10710))

        t.test2(FrechetCosineForsytheAndoyerLambert, (2.6319, 0, 0, 149,  5400),
                                                     (2.6319, 0, 0, 208, 10710))

        t.test2(FrechetCosineLaw, (2.63867, 0, 0, 149,  5400),
                                  (2.63867, 0, 0, 208, 10710))

        t.test2(FrechetEquirectangular, (7.1331,  8, 3, 138,  5400),
                                        (7.01295, 0, 0, 208, 10710))

        t.test2(FrechetEuclidean, (2.81124, 8, 3, 138,  5400),
                                  (2.74296, 0, 0, 208, 10710))

        t.test2(FrechetFlatLocal, (7.13778, 8, 3, 138,  5400),
                                  (6.92262, 0, 0, 208, 10710))

        t.test2(FrechetFlatPolar, (2.65039, 0, 0, 149,  5400),
                                  (2.65039, 0, 0, 208, 10710))

        t.test2(FrechetHaversine, (2.63867, 0, 0, 149,  5400),
                                  (2.63867, 0, 0, 208, 10710))

        t.test2(FrechetHubeny, (7.13778, 8, 3, 138,  5400),
                               (6.92262, 0, 0, 208, 10710))

        t.test2(FrechetThomas, (2.63187, 0, 0, 149,  5400),
                               (2.63187, 0, 0, 208, 10710))

        t.test2(FrechetVincentys, (2.63867, 0, 0, 149,  5400),
                                  (2.63867, 0, 0, 208, 10710))

        if geographiclib:
            t.test2(FrechetKarney, (151.09508, 0, 0, 149,  5400),
                                   (151.09508, 0, 0, 208, 10710))

        if t._testX:
            t.test2(FrechetExact, (151.09508, 0, 0, 149,  5400),
                                  (151.09508, 0, 0, 208, 10710))

    elif isWindows:  # Python 2
        t.test2(FrechetDegrees_, (182.5,  83, 45,   21,  5400),
                                 (175.75, 83, 56.5, 12, 10710))

        t.test2(FrechetRadians_, (3.18523, 83, 45,   21,  5400),
                                 (3.06742, 83, 56.5, 12, 10710))

        t.test2(FrechetCosineAndoyerLambert, (1.74692, 49, 27,  73,  5400),
                                             (1.74692, 49, 27, 105, 10710))

        t.test2(FrechetCosineForsytheAndoyerLambert, (1.74692, 49, 27,  73,  5400),
                                                     (1.74692, 49, 27, 105, 10710))

        t.test2(FrechetCosineLaw, (1.75068, 49, 27,  73,  5400),
                                  (1.75068, 49, 27, 105, 10710))

        t.test2(FrechetEquirectangular, (5.88254, 41, 18,    90,  5400),
                                        (5.90078, 40, 15.5, 137, 10710))

        t.test2(FrechetEuclidean, (2.52802, 50, 26,    73,  5400),
                                  (2.48948, 39, 22.5, 124, 10710))

        t.test2(FrechetFlatLocal, (5.85735, 42, 19,    88,  5400),
                                  (5.89746, 40, 15.5, 137, 10710))

        t.test2(FrechetFlatPolar, (2.10357, 84, 40,   25,  5400),
                                  (2.00246, 57, 40.5, 70, 10710))

        t.test2(FrechetHaversine, (1.75068, 49, 27,  73,  5400),
                                  (1.75068, 49, 27, 105, 10710))

        t.test2(FrechetHubeny, (5.85735, 42, 19,    88,  5400),
                               (5.89746, 40, 15.5, 137, 10710))

        t.test2(FrechetThomas, (1.74692, 49, 27,  73,  5400),
                               (1.74692, 49, 27, 105, 10710))

        t.test2(FrechetVincentys, (1.75068, 49, 27,  73,  5400),
                                  (1.75068, 49, 27, 105, 10710))

        if geographiclib:
            t.test2(FrechetKarney, (100.27538, 49, 27,  73,  5400),
                                   (100.27538, 49, 27, 105, 10710))

        if t._testX:
            t.test2(FrechetExact, (100.27538, 49, 27,  73,  5400),
                                  (100.27538, 49, 27, 105, 10710))

    else:  # Python 2, elsewhere
        t.test2(FrechetDegrees_, (288.0, 1, 1, 147,  5400),
                                 (288.0, 1, 1, 205, 10710))

        t.test2(FrechetRadians_, (5.02655, 1, 1, 147,  5400),
                                 (5.02655, 1, 1, 205, 10710))

        t.test2(FrechetCosineAndoyerLambert, (1.81201, 18, 14,   117,  5400),
                                             (1.83135,  3,  4.5, 196, 10710))

        t.test2(FrechetCosineForsytheAndoyerLambert, (1.81201, 18, 14,   117,  5400),
                                                     (1.83135,  3,  4.5, 196, 10710))

        t.test2(FrechetCosineLaw, (1.81341, 18, 14,   117,  5400),
                                  (1.83289,  3,  4.5, 196, 10710))

        t.test2(FrechetEquirectangular, ( 7.53702, 1, 3,   145,  5400),
                                        (12.58507, 0, 2.5, 203, 10710))

        t.test2(FrechetEuclidean, (2.80593, 1, 3,   145,  5400),
                                  (3.83981, 0, 2.5, 203, 10710))

        t.test2(FrechetFlatLocal, ( 7.55994, 1, 3,   145,  5400),
                                  (12.61423, 0, 2.5, 203, 10710))

        t.test2(FrechetFlatPolar, (2.88199, 89, 59, 1,  5400),
                                  (2.88199, 89, 59, 1, 10710))

        t.test2(FrechetHaversine, (1.81341, 18, 14,   117,  5400),
                                  (1.83289,  3,  4.5, 196, 10710))

        t.test2(FrechetHubeny, ( 7.55994, 1, 3,   145,  5400),
                               (12.61423, 0, 2.5, 203, 10710))

        t.test2(FrechetThomas, (1.81201, 18, 14,   117,  5400),
                               (1.83135,  3,  4.5, 196, 10710))

        t.test2(FrechetVincentys, (1.81341, 18, 14,   117,  5400),
                                  (1.83289,  3,  4.5, 196, 10710))

        if geographiclib:
            t.test2(FrechetKarney, (104.00172, 18, 14,   117,  5400),
                                   (105.26515,  3,  4.5, 196, 10710))

        if t._testX:
            t.test2(FrechetExact, (104.00172, 18, 14,   117,  5400),
                                  (105.26515,  3,  4.5, 196, 10710))

    if isPython3:
        from pygeodesy import ellipsoidalNvector, ellipsoidalVincenty, \
                              sphericalNvector, sphericalTrigonometry

        ms = (ellipsoidalVincenty,)
        if geographiclib:
            from pygeodesy import ellipsoidalKarney
            ms += (ellipsoidalKarney,)
        if t._testX:
            from pygeodesy import ellipsoidalExact
            ms += (ellipsoidalExact,)
            if GeodSolve:
                from pygeodesy import ellipsoidalGeodSolve
                ms += (ellipsoidalGeodSolve,)

        for m in ms:
            _ms = [m.LatLon(*ll.latlon) for ll in _ms]
            _ps = [m.LatLon(*ll.latlon) for ll in _ps]
            t.test2(FrechetDistanceTo, (16786640.7064, 0, 0, 149, 5400),
                                       (16786640.7064, 0, 0, 149, 5400), prec=4, known=True)

        for m in (ellipsoidalNvector, sphericalNvector, sphericalTrigonometry):
            _ms = [m.LatLon(*ll.latlon) for ll in _ms]
            _ps = [m.LatLon(*ll.latlon) for ll in _ps]
            t.test2(FrechetDistanceTo, (16810959.0015, 0, 0, 149, 5400),
                                       (16810959.0015, 0, 0, 149, 5400), prec=4, known=True)

    if coverage:  # for test coverage
        f = frechet_(_ms, _ps, distance=_distance, units='test')
        t.test('frechet_', f, "(178.5, 74, 56, 19, 5400, 'test')" if isPython3 else
                              "(288.0, 1, 1, 147, 5400, 'test')", known=True)

        t.test('[fi1]', fractional(_ms, f.fi1), '64.0°S, 096.0°E' if isPython3 else '38.0°S, 116.0°W', known=True)
        t.test('[fi2]', fractional(_ps, f.fi2), '41.0°S, 071.0°W' if isPython3 else '64.0°N, 121.0°E', known=True)

        t.test('[fi1]', fractional(_ms, f.fi1, LatLon=LatLon_).toRepr(), 'LatLon_(64.0°S, 096.0°E)' if isPython3 else
                                                                         'LatLon_(38.0°S, 116.0°W)', known=True)
        t.test('[fi2]', fractional(_ps, f.fi2, LatLon=LatLon_).toRepr(), 'LatLon_(41.0°S, 071.0°W)' if isPython3 else
                                                                         'LatLon_(64.0°N, 121.0°E)', known=True)

    t.results()
    t.exit()
