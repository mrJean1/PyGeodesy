
# -*- coding: utf-8 -*-

# Test the Hausdorff distances.  Very slow with coverage!

__all__ = ('Tests',)
__version__ = '21.08.14'

from base import GeodSolve, geographiclib, isPython3, isWindows, \
                 TestsBase

from pygeodesy import Datums, fstr, hausdorff_, \
                      HausdorffCosineAndoyerLambert, \
                      HausdorffCosineForsytheAndoyerLambert, \
                      HausdorffDegrees, HausdorffDistanceTo, \
                      HausdorffRadians, HausdorffCosineLaw, \
                      HausdorffEquirectangular, HausdorffEuclidean, \
                      HausdorffExact, HausdorffFlatLocal, \
                      HausdorffFlatPolar, HausdorffHaversine, \
                      HausdorffHubeny, HausdorffKarney, \
                      HausdorffThomas, HausdorffVincentys, \
                      LatLon_, randomrangenerator


class HausdorffDegrees_(HausdorffDegrees):
    '''Custom Hausdorff.'''
    def distance(self, p1, p2):
        dy, dx = abs(p1.lat - p2.lat), abs(p1.lon - p2.lon)
        if dx < dy:
            dx, dy = dy, dx
        return dx + dy * 0.5


class HausdorffRadians_(HausdorffRadians):
    '''Custom Hausdorff.'''
    def distance(self, p1, p2):
        dy, dx = abs(p1.phi - p2.phi), abs(p1.lam - p2.lam)
        if dx < dy:
            dx, dy = dy, dx
        return dx + dy * 0.5


def _tstr(t):
    s = list(t[:5])
    s[0] = fstr(t.hd, prec=5)
    s[4] = None if t.md is None else fstr(t.md, prec=5)
    return '(%s)' % (', '.join(map(str, s)),)


_rr = randomrangenerator('R')
# like <https://GitHub.com/mavillan/py-hausdorff> Example
_ms = [LatLon_(*_ll) for _ll in zip(_rr( -90,  90, 2),  # 90
                                    _rr(-180, 180, 4))]

_ps = [LatLon_(*_ll) for _ll in zip(_rr( -89,  90, 3),  # 60
                                    _rr(-179, 180, 6))]


class Tests(TestsBase):

    def test4(self, Hausdorff, *tests, **kwds):

        k = kwds.pop('known', False)
        for s, e, x, y in tests:
            h = Hausdorff(_ms, seed=s, **kwds)
            n = '%s (%s)' % (h.named, h.units)

            t = _tstr(h.directed(_ps, early=e))
            self.test(n, t, x, known=k)  # + (h.units,)

            t = _tstr(h.symmetric(_ps, early=e))
            self.test(n, t, y, known=k)  # + (h.units,)

        self.testCopy(h)

    def test4_(self, H_, *tests):

        h = HausdorffDegrees_(_ms)
        d = h.distance

        for s, e, x, y in tests:

            t = _tstr(H_(_ms, _ps, both=False, distance=d, early=e, seed=s))
            self.test(hausdorff_.__name__, t, x)

            t = _tstr(H_(_ms, _ps, both=True,  distance=d, early=e, seed=s))
            self.test(hausdorff_.__name__, t, y)


if __name__ == '__main__':  # MCCABE 13

    def _4(t5directed, t5symmetric):
        # generate Hausdorff(*args) and results for
        # 4 different seed= and early= test cases
        yield None, False, t5directed, t5symmetric
        yield 'X',  False, t5directed, t5symmetric

        t4directed  = t5directed[:4] + (None,)
        t4symmetric = t5symmetric[:4] + (None,)
        yield None, True, t4directed, t4symmetric
        yield 'X',  True, t4directed, t4symmetric

    t = Tests(__file__, __version__)

    # check repeatable random ranges
    r1 = randomrangenerator('R')
    r2 = randomrangenerator('R')
    for n in (0, 1, 2, 8, 32, 128):
        t.test('randomrange[%d]' % (n,), tuple(r1(n)), tuple(r2(n)))

    if isPython3:  # XXX different Random?
        t.test4(HausdorffDegrees_, *_4((40.0, 22,  6,  90, 18.16111),
                                       (48.0, 38, 36, 150, 17.30667)))

        t.test4(HausdorffRadians_, *_4((0.69813, 22,  6,  90, 0.31697),
                                       (0.83776, 38, 36, 150, 0.30206)))

        t.test4(HausdorffCosineAndoyerLambert, *_4((0.501, 35, 3,  90, 0.21184),
                                                   (0.501, 35, 3, 150, 0.20084)))

        t.test4(HausdorffCosineForsytheAndoyerLambert, *_4((0.501, 35, 3,  90, 0.21184),
                                                           (0.501, 35, 3, 150, 0.20084)))

        t.test4(HausdorffCosineLaw, *_4((0.50097, 35, 3,  90, 0.212),
                                        (0.50097, 35, 3, 150, 0.20099)))

        t.test4(HausdorffEquirectangular, *_4((0.25113, 35,  3,  90, 0.05965),
                                              (0.25113, 35,  3, 150, 0.05532)))

        t.test4(HausdorffEuclidean, *_4((0.52803, 56, 51,  90, 0.22571),  # XXX different i, j?
                                        (0.52803, 56, 51, 150, 0.21579)))

        t.test4(HausdorffFlatLocal, *_4((0.25116, 35, 3,  90, 0.05951),
                                        (0.25116, 35, 3, 150, 0.0552)))

        t.test4(HausdorffFlatPolar, *_4((0.69043, 4, 31,  90, 0.22286),
                                        (0.69043, 4, 31, 150, 0.21318)))

        t.test4(HausdorffHaversine, *_4((0.50097, 35, 3,  90, 0.212),
                                        (0.50097, 35, 3, 150, 0.20099)))

        t.test4(HausdorffHubeny, *_4((0.25116, 35, 3,  90, 0.05951),
                                     (0.25116, 35, 3, 150, 0.0552)))

        t.test4(HausdorffThomas, *_4((0.501, 35, 3,  90, 0.21184),
                                     (0.501, 35, 3, 150, 0.20084)))

        t.test4(HausdorffVincentys, *_4((0.50097, 35, 3,  90, 0.212),
                                        (0.50097, 35, 3, 150, 0.20099)))

        if geographiclib:
            t.test4(HausdorffKarney, *_4((28.79903, 35, 3,  90, 12.16138),
                                         (28.79903, 35, 3, 150, 11.53021)),
                                     datum=Datums.WGS84)

        if t._testX:
            t.test4(HausdorffExact, *_4((28.79903, 35, 3,  90, 12.16138),
                                        (28.79903, 35, 3, 150, 11.53021)),
                                    datum=Datums.WGS84)

        t.test4_(hausdorff_, *_4((40.0, 22,  6,  90, 18.16111),
                                 (48.0, 38, 36, 150, 17.30667)))

    elif isWindows:  # Python 2
        t.test4(HausdorffDegrees_, *_4((51.0, 10, 50,  90, 15.92222),
                                       (51.0, 10, 50, 150, 14.32667)))

        t.test4(HausdorffRadians_, *_4((0.89012, 10, 50,  90, 0.2779),
                                       (0.89012, 10, 50, 150, 0.25005)))

        t.test4(HausdorffCosineAndoyerLambert, *_4((0.52412, 74, 45,  90, 0.18171),
                                                   (0.52412, 74, 45, 150, 0.16541)))

        t.test4(HausdorffCosineForsytheAndoyerLambert, *_4((0.52412, 74, 45,  90, 0.18171),
                                                           (0.52412, 74, 45, 150, 0.16541)))

        t.test4(HausdorffCosineLaw, *_4((0.52674, 74, 45,  90, 0.18192),
                                        (0.52674, 74, 45, 150, 0.16555)))

        t.test4(HausdorffEquirectangular, *_4((0.30982, 40, 0,  90, 0.04843),
                                              (0.30982, 40, 0, 150, 0.04072)))

        t.test4(HausdorffEuclidean, *_4((0.58262, 69, 13,  90, 0.19621),
                                        (0.58262, 69, 13, 150, 0.1787)))

        t.test4(HausdorffFlatLocal, *_4((0.3105, 40, 0,  90, 0.04826),
                                        (0.3105, 40, 0, 150, 0.04061)))

        t.test4(HausdorffFlatPolar, *_4((0.89288, 53, 15,  90, 0.22047),
                                        (0.89288, 53, 15, 150, 0.20447)))

        t.test4(HausdorffHaversine, *_4((0.52674, 74, 45,  90, 0.18192),
                                        (0.52674, 74, 45, 150, 0.16555)))

        t.test4(HausdorffHubeny, *_4((0.3105, 40, 0,  90, 0.04826),
                                     (0.3105, 40, 0, 150, 0.04061)))

        t.test4(HausdorffThomas, *_4((0.52412, 74, 45,  90, 0.18171),
                                     (0.52412, 74, 45, 150, 0.16541)))

        t.test4(HausdorffVincentys, *_4((0.52674, 74, 45,  90, 0.18192),
                                        (0.52674, 74, 45, 150, 0.16555)))

        if geographiclib:
            t.test4(HausdorffKarney, *_4((30.11794, 74, 45,  90, 10.43166),
                                         (30.11794, 74, 45, 150,  9.49554)),
                                     datum=Datums.WGS84)

        if t._testX:
            t.test4(HausdorffExact, *_4((30.11794, 74, 45,  90, 10.43166),
                                        (30.11794, 74, 45, 150,  9.49554)),
                                    datum=Datums.WGS84)

        t.test4_(hausdorff_, *_4((51.0, 10, 50,  90, 15.92222),
                                 (51.0, 10, 50, 150, 14.32667)))

    else:  # Python 2, elsewhere
        t.test4(HausdorffDegrees_, *_4((50.5, 61, 7,  90, 18.45556),
                                       (50.5, 61, 7, 150, 16.05)))

        t.test4(HausdorffRadians_, *_4((0.88139, 61, 7,  90, 0.32211),
                                       (0.88139, 61, 7, 150, 0.28013)))

        t.test4(HausdorffCosineAndoyerLambert, *_4((0.56263, 49, 29,  90, 0.19764),
                                                   (0.56263, 49, 29, 150, 0.17591)))

        t.test4(HausdorffCosineForsytheAndoyerLambert, *_4((0.56263, 49, 29,  90, 0.19764),
                                                           (0.56263, 49, 29, 150, 0.17591)))

        t.test4(HausdorffCosineLaw, *_4((0.56202, 49, 29,  90, 0.19776),
                                        (0.56202, 49, 29, 150, 0.176)))

        t.test4(HausdorffEquirectangular, *_4((0.33016, 49, 29,  90, 0.06075),
                                              (0.33016, 49, 29, 150, 0.04918)))

        t.test4(HausdorffEuclidean, *_4((0.61785, 49, 29,  90, 0.21327),
                                        (0.61785, 49, 29, 150, 0.19061)))

        t.test4(HausdorffFlatLocal, *_4((0.33092, 49, 29,  90, 0.06067),
                                        (0.33092, 49, 29, 150, 0.04912)))

        t.test4(HausdorffFlatPolar, *_4((0.96853, 20, 46,  90, 0.24569),
                                        (0.96853, 20, 46, 150, 0.21335)))

        t.test4(HausdorffHaversine, *_4((0.56202, 49, 29,  90, 0.19776),
                                        (0.56202, 49, 29, 150, 0.176)))

        t.test4(HausdorffHubeny, *_4((0.33092, 49, 29,  90, 0.06067),
                                     (0.33092, 49, 29, 150, 0.04912)))

        t.test4(HausdorffThomas, *_4((0.56263, 49, 29,  90, 0.19764),
                                     (0.56263, 49, 29, 150, 0.17591)))

        t.test4(HausdorffVincentys, *_4((0.56202, 49, 29,  90, 0.19776),
                                        (0.56202, 49, 29, 150, 0.176)))

        if geographiclib:
            t.test4(HausdorffKarney, *_4((32.28234, 49, 29,  90, 11.34653),
                                         (32.28234, 49, 29, 150, 10.09914)),
                                     datum=Datums.WGS84)
        if t._testX:
            t.test4(HausdorffExact, *_4((32.28234, 49, 29,  90, 11.34653),
                                        (32.28234, 49, 29, 150, 10.09914)),
                                    datum=Datums.WGS84)

        t.test4_(hausdorff_, *_4((50.5, 61, 7,  90, 18.45556),
                                 (50.5, 61, 7, 150, 16.05)))

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
            t.test4(HausdorffDistanceTo, *_4((3195418.34044, 35, 3,  90, 1351164.35981),
                                             (3195418.34044, 35, 3, 150, 1280992.80341)), known=True)

        for m in (ellipsoidalNvector, sphericalNvector, sphericalTrigonometry):
            _ms = [m.LatLon(*ll.latlon) for ll in _ms]
            _ps = [m.LatLon(*ll.latlon) for ll in _ps]
            t.test4(HausdorffDistanceTo, *_4((3191685.70841, 35, 3,  90, 1350684.04012),
                                             (3191685.70841, 35, 3, 150, 1280478.58707)), known=True)

    t.results()
    t.exit()
