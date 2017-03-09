
# -*- coding: utf-8 -*-

# Test ellipsoidal earth model functions and methods.

__all__ = ('Tests',)
__version__ = '17.03.08'

from tests import Tests as _Tests

from geodesy import F_D, VincentyError, compassDMS, Datums, fStr


class Tests(_Tests):

    def testEllipsoidal(self, LatLon, Nvector=None, Cartesian=None):
        # ellipsoidal modules tests
        p = LatLon(51.4778, -0.0016, 0, Datums.WGS84)
        d = p.convertDatum(Datums.OSGB36)
        self.test('convertDatum', d, '51.477284°N, 000.00002°E, -45.91m')  # 51.4773°N, 000.0000°E, -45.91m
        self.test('convertDatum', d.toStr(F_D, prec=4), '51.4773°N, 000.0°E, -45.91m')

        if Cartesian:
            c = Cartesian(3980581, 97, 4966825)
            n = c.toNvector()  # {x: 0.6228, y: 0.0000, z: 0.7824, h: 0.0000}  # XXX height
            self.test('toNVector', n.toStr(4), '(0.6228, 0.0, 0.7824, +0.24)')
            c = n.toCartesian()
            self.test('toCartesian', c.toStr(0), '[3980581, 97, 4966825]')

        if Nvector:
            n = Nvector(0.5, 0.5, 0.7071)
            c = n.toCartesian()  # [3194434, 3194434, 4487327]
            self.test('toCartesian', c, '[3194434.411, 3194434.411, 4487326.82]')
            p = c.toLatLon()  # 45.0°N, 45.0°E
            self.test('toLatLon', p.toStr('d', 2), '45.0°N, 045.0°E, +0.00m')  # 45.0°N, 45.0°E

            self.test('Nvector', n, '(0.5, 0.5, 0.7071)')
            n = Nvector(0.5, 0.5, 0.7071, 1).toStr(3)
            self.test('Nvector', n, '(0.5, 0.5, 0.707, +1.00)')

    def testVincenty(self, LatLon, datum):
        d = datum
        n = ' (%s)' % (d.name,)

        Newport_RI = LatLon(41.49008, -71.312796, datum=d)
        Cleveland_OH = LatLon(41.499498, -81.695391, datum=d)
        m = Newport_RI.distanceTo(Cleveland_OH)
        self.test('distanceTo' + n, '%.5f' % m, '866455.43292')

        try:
            t = None
            m = Newport_RI.distanceTo(Newport_RI)
        except VincentyError as x:
            t = x  # Python 3+
        self.test('VincentyError' + n, t, 'LatLon(41°29′24.29″N, 071°18′46.07″W) coincident with LatLon(41°29′24.29″N, 071°18′46.07″W)')

        if hasattr(LatLon, 'toCartesian'):
            try:
                m = Newport_RI.distanceTo(Cleveland_OH.convertDatum(Datums.OSGB36))
                self.test('ValueError' + n, None, 'other Ellipsoid mistmatch: ...' + d.ellipsoid.name)
            except ValueError as x:
                self.test('ValueError' + n, x, 'other Ellipsoid mistmatch: Ellipsoids.Airy1830 vs Ellipsoids.' + d.ellipsoid.name)
            except Exception as x:
                self.test('ValueError' + n, x, 'ValueError ...' + d.ellipsoid.name)

        p = LatLon(50.06632, -5.71475, datum=d)
        q = LatLon(58.64402, -3.07009, datum=d)
        m = p.distanceTo(q)
        self.test('distanceTo' + n, '%.4f' % m, '969954.1663')

        self.test('copy', p.copy().equals(p), 'True')

        t = p.distanceTo3(q)
        t = fStr(t, prec=6)
        self.test('distanceTo3' + n, t, '969954.166314, 9.141877, 11.29722')

        p = LatLon(37.95103, 144.42487, datum=d)
        q = LatLon(37.65280, 143.9265, datum=d)
        m = p.distanceTo(q)
        self.test('distanceTo' + n, '%.3f' % m, '54973.295')

        t = p.distanceTo3(q)
        t = fStr(t, prec=5)
        self.test('distanceTo3' + n, t, '54973.29527, 126.86992, 127.17539')

        p = LatLon(-37.95103, 144.42487, datum=d)
        p, f = p.destination2(54972.271, 306.86816)
        t = p.toStr(F_D) + ', ' + compassDMS(f, prec=4)
        self.test('destination2' + n, t, '37.652818°S, 143.926498°E, 307.1736°NW')


if __name__ == '__main__':

    from geodesy import ellipsoidalNvector as N
    t = Tests(__file__, __version__, N)
    t.testLatLon(N.LatLon)
    t.testVectorial(N.LatLon, N.Nvector, N.sumOf)
    t.testEllipsoidal(N.LatLon, N.Nvector, N.Cartesian)
    t.results()

    from geodesy import ellipsoidalVincenty as V
    t = Tests(__file__, __version__, V)
    t.testLatLon(V.LatLon, Sph=False)
    for d in (Datums.WGS84, Datums.NAD83,):  # Datums.Sphere):
        t.testVincenty(V.LatLon, d)
    t.results()
    t.exit()
