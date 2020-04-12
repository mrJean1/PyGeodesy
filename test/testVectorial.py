
# -*- coding: utf-8 -*-

# Test module attributes.

__all__ = ('Tests',)
__version__ = '20.04.06'

from base import coverage, TestsBase

from pygeodesy import F_D, fstr, sphericalNvector


class Tests(TestsBase):

    def testNvectorBase(self, module, **kwds):

        try:
            Nvector = module.Nvector
            c = Nvector.__name__
        except AttributeError:
            Nvector = module.NvectorBase
            c = 'Vector4Tuple'
        self.subtitle(module, Nvector.__name__)

        v = Nvector(0.500, 0.500, 0.707, **kwds)
        s = module.sumOf((v, v), h=0, name='sumOf')
        self.test('sumOf', s.__class__.__name__, c)

        p = v.toLatLon(LatLon=None)
        c = v.toCartesian(Cartesian=None)
        self.test('ecef.x, .y, .z', fstr(p[:3],  prec=5), fstr(c[:3],  prec=5))
        self.test('ecef.lat, .lon', fstr(p[3:5], prec=6), fstr(c[3:5], prec=6))
        self.test('ecef.height', fstr(p.height, prec=6), fstr(c.height, prec=6), known=True)
        if c.M is not None:
            self.test('ecef.M', fstr(p.M, prec=9), fstr(c.M, prec=9))

        if coverage:
            from pygeodesy.named import LatLon2Tuple, LatLon3Tuple, \
                                        PhiLam2Tuple, PhiLam3Tuple

            self.test('.isEllipsoidal', v.isEllipsoidal, not v.isSpherical)
            self.test('.isSpherical',   v.isSpherical,   not v.isEllipsoidal)

            self.test('.latlon', v.latlon, LatLon2Tuple(v.lat, v.lon))
            self.test('.philam', v.philam, PhiLam2Tuple(v.phi, v.lam))

            self.test('.latlonheight', v.latlonheight, LatLon3Tuple(v.lat, v.lon, float(v.h)))
            self.test('.philamheight', v.philamheight, PhiLam3Tuple(v.phi, v.lam, float(v.h)))

            t = v.parse('0.5, 0.5, 0.707')
            self.test('parse', t, v)
            self.test('cmp', t.cmp(v), 0)

            self.test('eq', t == v, True)
            self.test('ge', t >= v, True)
            self.test('gt', t >  v, False)
            self.test('le', t <= v, True)
            self.test('lt', t <  v, False)
            self.test('ne', t != v, False)

            m = t * 2
            self.test('*', m, '(1.0, 1.0, 1.414)')
            self.test('+', t + v, m)
            self.test('/', m / 2, t)
            self.test('-', m - t, t)

            m = v.__matmul__(t)
            self.test('@', m, '(0.0, 0.0, 0.0)')
            r = t.__rmatmul__(m)
            self.test('@', r, m)

            r = v.rotate(m, 45)
            self.test('rotate', r, '(0.26268, 0.26268, 0.37143)')

    def testVectorial(self, module):  # MCCABE 14

        self.subtitle(module, 'Vectorial')

        LatLon, Nvector, meanOf, sumOf = module.LatLon, module.Nvector, module.meanOf, module.sumOf

        if hasattr(LatLon, 'crossTrackDistanceTo'):
            p = LatLon(53.2611, -0.7972)
            s = LatLon(53.3206, -1.7297)
            d = p.crossTrackDistanceTo(s, 96.0)
            self.test('crossTrackDistanceTo', d, '-305.67', fmt='%.2f')  # -305.7
            e = LatLon(53.1887, 0.1334)
            d = p.crossTrackDistanceTo(s, e)
            self.test('crossTrackDistanceTo', d, '-307.55', fmt='%.2f')  # -307.5

        if hasattr(LatLon, 'enclosedby'):
            r = LatLon(45,1), LatLon(45,2), LatLon(46,2), LatLon(46,1)
            p = LatLon(45.1, 1.1)
            self.test('enclosedby', p.enclosedby(r), True)
            r = LatLon(45,1), LatLon(46,2), LatLon(45,2), LatLon(46,1)
            self.test('enclosedby', p.enclosedby(r), False)

        v = Nvector(0.500, 0.500, 0.707)
        p = v.toLatLon()
        self.test('toLatLon', p, '44.995674°N, 045.0°E')  # 45.0°N, 45.0°E
        c = p.toNvector()
        self.test('toNvector', c, '(0.50004, 0.50004, 0.70705)')  # 0.500, 0.500, 0.707
        self.test('isequalTo', c.isequalTo(v), False)
        self.test('isequalTo', c.isequalTo(v, units=True), True)
        self.test('length', v.length, '0.99992449715', fmt='%.11f')
        self.test('length', c.length, '1.0')

        s = meanOf((p, p, p, p), height=0, LatLon=LatLon)
        self.test('meanOf', s, '44.995674°N, 045.0°E')
        self.test('meanOf', s.__class__.__name__, LatLon.__name__)

        class Nv(Nvector):
            pass
        v = Nvector(52.205, 0.119, 0.0)
        s = sumOf((v, c), Vector=Nv, h=0, name='sumOf')
        self.test('sumOf', s, '(52.70504, 0.61904, 0.70705)')
        self.test('sumOf', s.__class__.__name__, 'Nv')
        self.test('sumOf', s._name, 'sumOf')
        self.test('length', s.length, '52.7134151513', fmt='%.10f')

        c = v.copy()
        self.test('copy', c.isequalTo(v), True)
        self.test('length', v.length, '52.2051356286', fmt='%.10f')
        self.test('length', c.length, '52.2051356286', fmt='%.10f')

        if module is sphericalNvector:  # coverage
            c = p.toCartesian()
            self.test('toCartesian', c, '[3185744.919, 3185744.919, 4504643.315]')
            self.test('toLatLon',  c.toLatLon(), p, known=True)  # '44.995674°N, 045.0°E, -0.00m'
            self.test('toNvector', c.toNvector(), '(0.50004, 0.50004, 0.70705, -0.00)', known=True)

        if hasattr(LatLon, 'intersection'):
            # <https://GitHub.com/ChrisVeness/geodesy/blob/master/test/latlon-vectors-tests.js>
            p = LatLon(1, 1)
            i = p.intersection(LatLon(2, 2), LatLon(1, 4), LatLon(2, 3))
            self.test('intersection', i, '02.499372°N, 002.5°E')

        if hasattr(LatLon, 'isEnclosedBy'):
            p = LatLon(45.1, 1.1)

            b = LatLon(45, 1), LatLon(45, 2), LatLon(46, 2), LatLon(46, 1)
            for _ in self.testiter():
                self.test('isEnclosedBy', p.isEnclosedBy(b), True)

            b = LatLon(45, 1), LatLon(45, 3), LatLon(46, 2), LatLon(47, 3), LatLon(47, 1)
            for _ in self.testiter():
                try:
                    self.test('isEnclosedBy', p.isEnclosedBy(b), True)  # Nvector
                except ValueError as x:
                    t = ' '.join(str(x).split()[:3] + ['...)'])
                    self.test('isEnclosedBy', t, 'non-convex: (LatLon(45°00′00.0″N, 001°00′00.0″E), ...)')  # Trig

        if hasattr(LatLon, 'isWithin'):
            # courtesy of Paulius Šarka  psarka  Aug 30, 2017
            p = LatLon(1, 1).isWithin(LatLon(2, 2), LatLon(2, 2))
            self.test('isWithin', p, False)
            p = LatLon(2, 2).isWithin(LatLon(2, 2), LatLon(2, 2))
            self.test('isWithin', p, True)

        if hasattr(LatLon, 'nearestOn'):
            s1 = LatLon(51.0, 1.0)
            s2 = LatLon(51.0, 2.0)
            s = LatLon(51.0, 1.9)
            p = s.nearestOn(s1, s2)  # 51.0004°N, 001.9000°E
            self.test('nearestOn', p.toStr(F_D, prec=4), '51.0004°N, 001.9°E')
            self.test('nearestOn', isinstance(p, LatLon), True)

            d = p.distanceTo(s)  # 42.71 m
            self.test('distanceTo', d, 42.712, fmt='%.3f')
            s = LatLon(51.0, 2.1)
            p = s.nearestOn(s1, s2)  # 51.0000°N, 002.0000°E
            self.test('nearestOn', p.toStr(F_D), '51.0°N, 002.0°E')
            self.test('nearestOn', isinstance(p, LatLon), True)

            # courtesy AkimboEG on GitHub
            s1 = LatLon(0, 0)
            s2 = LatLon(0, 1)
            s = LatLon(1, 0)
            p = s.nearestOn(s1, s2)  # 0.0°N, 0.0°E
            self.test('nearestOn', p, '00.0°N, 000.0°E')
            self.test('nearestOn', isinstance(p, LatLon), True)

            p = LatLon(10, -140).nearestOn(LatLon(0, 20), LatLon(0, 40))
            self.test('nearestOn', p, '00.0°N, 020.0°E')
            self.test('nearestOn', isinstance(p, LatLon), True)

            # courtesy of Paulius Šarka  psarka  Aug 30, 2017
            p = LatLon(1, 1).nearestOn(LatLon(2, 2), LatLon(2, 2))
            self.test('nearestOn', p, '02.0°N, 002.0°E')
            p = LatLon(2, 2).nearestOn(LatLon(2, 2), LatLon(2, 2))
            self.test('nearestOn', p, '02.0°N, 002.0°E')

        if hasattr(LatLon, 'triangulate'):
            # courtesy of pvezid  Feb 10, 2017
            p = LatLon("47°18.228'N","002°34.326'W")  # Basse Castouillet
            self.test('BasseC', p, '47.3038°N, 002.5721°W')
            s = LatLon("47°18.664'N","002°31.717'W")  # Basse Hergo
            self.test('BasseH', s, '47.311067°N, 002.528617°W')
            t = p.triangulate(7, s, 295)
            self.test('triangulate', t, '47.323667°N, 002.568501°W')
            self.test('triangulate', isinstance(t, LatLon), True)

        if hasattr(LatLon, 'trilaterate'):
            # <https://GitHub.com/chrisveness/geodesy/blob/master/test/latlon-nvector-spherical-tests.js>
            p = LatLon(37.418436, -121.963477)
            t = p.trilaterate(265.710701754, LatLon(37.417243, -121.961889), 234.592423446,
                                             LatLon(37.418692, -121.960194), 54.8954278262)
            self.test('trilaterate', t, '37.419078°N, 121.960579°W')
            self.test('trilaterate', isinstance(t, LatLon), True)

            # courtesy Carlos Freitas <https://GitHub.com/mrJean1/PyGeodesy/issues/33>
            b1 = LatLon(-8.068361, -34.892722)
            b2 = LatLon(-8.075917, -34.894611)
            b3 = LatLon(-8.076361, -34.908000)
            p  = LatLon(-8.068912, -34.888699)
            d1 = b1.distanceTo(p)
            d2 = b2.distanceTo(p)
            d3 = b3.distanceTo(p)
            t = b1.trilaterate(d1, b2, d2, b3, d3)
            self.test('trilaterate', t, p)
            self.test('trilaterate', isinstance(t, LatLon), True)
            t = b1.trilaterate(d1, b2, d2, b3, d3, useZ=True)
            self.test('trilaterate', t, p)
            self.test('trilaterate', isinstance(t, LatLon), True)

        self.testNvectorBase(module)


if __name__ == '__main__':

    from pygeodesy import Datums, ellipsoidalNvector, nvectorBase

    t = Tests(__file__, __version__)

    t.testVectorial(ellipsoidalNvector)
    t.testVectorial(sphericalNvector)

    t.testNvectorBase(nvectorBase, datum=Datums.Sphere)
    t.testNvectorBase(nvectorBase, datum=Datums.WGS84)

    t.results()
    t.exit()
