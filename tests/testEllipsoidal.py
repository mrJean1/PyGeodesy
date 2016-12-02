
# -*- coding: utf-8 -*-

# Test ellipsoidal earth model functions and methods.

__version__ = '16.11.28'

if __name__ == '__main__':

    from tests import Tests

    from geodesy import ellipsoidalNvector as N
    t = Tests(__file__, __version__, N)
#   t.testLatLon(N.LatLon)
    t.testVectorial(N.LatLon, N.Nvector)
    t.testEllipsoidal(N.LatLon, N.Nvector, N.Cartesian)
    t.results()

    from geodesy import Datums, ellipsoidalVincenty as V, VincentyError
    t = Tests(__file__, __version__, V)
    for d in (Datums.WGS84, Datums.NAD83,):  # Datums.Sphere):
        t.testVincenty(V.LatLon, d, VincentyError)
    t.results()

    # Typical test results (on MacOS X):

    # testing geodesy.ellipsoidalNvector version 16.11.11
    # test 1 toLatLon: 44.995674°N, 045.0°E
    # test 2 toNvector: (0.50004, 0.50004, 0.70705)
    # test 3 equals: False
    # test 4 equals: True
    # test 5 copy: True
    # test 6 convertDatum: 51.477284°N, 000.00002°E, -45.91m
    # test 7 convertDatum: 51.4773°N, 000.0°E, -45.91m
    # test 8 toNVector: (0.6228, 0.0, 0.7824, +0.24)
    # test 9 toCartesian: [3980581, 97, 4966825]
    # test 10 toCartesian: [3194434.411, 3194434.411, 4487326.82]
    # test 11 toLatLon: 45.0°N, 045.0°E, +0.00m
    # test 12 Nvector: (0.5, 0.5, 0.7071)
    # test 13 Nvector: (0.5, 0.5, 0.707, +1.00)
    # all geodesy.ellipsoidalNvector tests passed (Python 2.7.10)

    # testing geodesy.ellipsoidalVincenty version 16.11.28
    # test 1 distanceTo (WGS84): 866455.43292
    # test 2 VincentyError (WGS84): LatLon(41°29′24.29″N, 071°18′46.07″W) coincident with LatLon(41°29′24.29″N, 071°18′46.07″W)
    # test 3 ValueError (WGS84): other Ellipsoid mistmatch: Ellipsoids.Airy1830 vs Ellipsoids.WGS84
    # test 4 distanceTo (WGS84): 969954.1663
    # test 5 copy: True
    # test 6 distanceTo3 (WGS84): 969954.166314, 9.141877, 11.29722
    # test 7 distanceTo (WGS84): 54973.295
    # test 8 distanceTo3 (WGS84): 54973.29527, 126.86992, 127.17539
    # test 9 destination2 (WGS84): 37.652818°S, 143.926498°E, 307.1736°NW
    # test 10 distanceTo (NAD83): 866455.43292
    # test 11 VincentyError (NAD83): LatLon(41°29′24.29″N, 071°18′46.07″W) coincident with LatLon(41°29′24.29″N, 071°18′46.07″W)
    # test 12 ValueError (NAD83): other Ellipsoid mistmatch: Ellipsoids.Airy1830 vs Ellipsoids.GRS80
    # test 13 distanceTo (NAD83): 969954.1663
    # test 14 copy: True
    # test 15 distanceTo3 (NAD83): 969954.166314, 9.141877, 11.29722
    # test 16 distanceTo (NAD83): 54973.295
    # test 17 distanceTo3 (NAD83): 54973.29527, 126.86992, 127.17539
    # test 18 destination2 (NAD83): 37.652818°S, 143.926498°E, 307.1736°NW
    # all geodesy.ellipsoidalVincenty tests passed (Python 2.7.10)

    # testing ellipsoidalNvector version 16.11.11
    # test 1 toLatLon: 44.995674°N, 045.0°E
    # test 2 toNvector: (0.50004, 0.50004, 0.70705)
    # test 3 equals: False
    # test 4 equals: True
    # test 5 copy: True
    # test 6 convertDatum: 51.477284°N, 000.00002°E, -45.91m
    # test 7 convertDatum: 51.4773°N, 000.0°E, -45.91m
    # test 8 toNVector: (0.6228, 0.0, 0.7824, +0.24)
    # test 9 toCartesian: [3980581, 97, 4966825]
    # test 10 toCartesian: [3194434.411, 3194434.411, 4487326.82]
    # test 11 toLatLon: 45.0°N, 045.0°E, +0.00m
    # test 12 Nvector: (0.5, 0.5, 0.7071)
    # test 13 Nvector: (0.5, 0.5, 0.707, +1.00)
    # all ellipsoidalNvector tests passed (Python 3.5.2)

    # testing ellipsoidalVincenty version 16.11.28
    # test 1 distanceTo (WGS84): 866455.43292
    # test 2 VincentyError (WGS84): LatLon(41°29′24.29″N, 071°18′46.07″W) coincident with LatLon(41°29′24.29″N, 071°18′46.07″W)
    # test 3 ValueError (WGS84): other Ellipsoid mistmatch: Ellipsoids.Airy1830 vs Ellipsoids.WGS84
    # test 4 distanceTo (WGS84): 969954.1663
    # test 5 copy: True
    # test 6 distanceTo3 (WGS84): 969954.166314, 9.141877, 11.29722
    # test 7 distanceTo (WGS84): 54973.295
    # test 8 distanceTo3 (WGS84): 54973.29527, 126.86992, 127.17539
    # test 9 destination2 (WGS84): 37.652818°S, 143.926498°E, 307.1736°NW
    # test 10 distanceTo (NAD83): 866455.43292
    # test 11 VincentyError (NAD83): LatLon(41°29′24.29″N, 071°18′46.07″W) coincident with LatLon(41°29′24.29″N, 071°18′46.07″W)
    # test 12 ValueError (NAD83): other Ellipsoid mistmatch: Ellipsoids.Airy1830 vs Ellipsoids.GRS80
    # test 13 distanceTo (NAD83): 969954.1663
    # test 14 copy: True
    # test 15 distanceTo3 (NAD83): 969954.166314, 9.141877, 11.29722
    # test 16 distanceTo (NAD83): 54973.295
    # test 17 distanceTo3 (NAD83): 54973.29527, 126.86992, 127.17539
    # test 18 destination2 (NAD83): 37.652818°S, 143.926498°E, 307.1736°NW
    # all ellipsoidalVincenty tests passed (Python 3.5.2)
