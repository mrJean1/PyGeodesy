
# -*- coding: utf-8 -*-

# Test module attributes.

__all__ = ('Tests',)
__version__ = '17.08.31'

from base import TestsBase

from pygeodesy import F_D


class Tests(TestsBase):

    def testVectorial(self, module):

        self.subtitle(module, 'Vectorial')

        LatLon, Nvector, sumOf = module.LatLon, module.Nvector, module.sumOf

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

#       p = meanOf(r)
#       self.test('meanOf', p, '')

        v = Nvector(0.500, 0.500, 0.707)
        p = v.toLatLon()
        self.test('toLatLon', p, '44.995674°N, 045.0°E')  # 45.0°N, 45.0°E
        c = p.toNvector()
        self.test('toNvector', c, '(0.50004, 0.50004, 0.70705)')  # 0.500, 0.500, 0.707
        self.test('equals', c.equals(v), False)
        self.test('equals', c.equals(v, units=True), True)
        self.test('length', v.length, '0.99992449715',  fmt='%.11f')
        self.test('length', c.length, '1.0')

        class Nv(Nvector):
            pass
        v = Nvector(52.205, 0.119, 0.0)
        s = sumOf((v, c), Vector=Nv, h=0)
        self.test('sumOf', s, '(52.70504, 0.61904, 0.70705)')
        self.test('sumOf', s.__class__.__name__, 'Nv')
        self.test('length', s.length, '52.7134151513',  fmt='%.10f')

        c = v.copy()
        self.test('copy', c.equals(v), True)
        self.test('length', v.length, '52.2051356286',  fmt='%.10f')
        self.test('length', c.length, '52.2051356286',  fmt='%.10f')

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


if __name__ == '__main__':

    from pygeodesy import ellipsoidalNvector, sphericalNvector

    t = Tests(__file__, __version__)

    t.testVectorial(ellipsoidalNvector)

    t.testVectorial(sphericalNvector)

    t.results()
    t.exit()
