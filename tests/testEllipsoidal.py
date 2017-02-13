
# -*- coding: utf-8 -*-

# Test ellipsoidal earth model functions and methods.

__version__ = '17.02.09'

if __name__ == '__main__':

    from tests import Tests

    from geodesy import ellipsoidalNvector as N
    t = Tests(__file__, __version__, N)
#   t.testLatLon(N.LatLon)
    t.testVectorial(N.LatLon, N.Nvector, N.sumOf)
    t.testEllipsoidal(N.LatLon, N.Nvector, N.Cartesian)
    t.results()

    from geodesy import Datums, ellipsoidalVincenty as V, VincentyError
    t = Tests(__file__, __version__, V)
    for d in (Datums.WGS84, Datums.NAD83,):  # Datums.Sphere):
        t.testVincenty(V.LatLon, d, VincentyError)
    t.results()
    t.exit()

    # Typical test results (on MacOS 10.12.3):

    # testing geodesy.ellipsoidalNvector version 17.02.13
    # test 1 toLatLon: 44.995674°N, 045.0°E
    # test 2 toNvector: (0.50004, 0.50004, 0.70705)
    # test 3 equals: False
    # test 4 equals: True
    # test 5 sumOf: (52.70504, 0.61904, 0.70705)
    # test 6 sumOf: Nv
    # test 7 copy: True
    # test 8 convertDatum: 51.477284°N, 000.00002°E, -45.91m
    # test 9 convertDatum: 51.4773°N, 000.0°E, -45.91m
    # test 10 toNVector: (0.6228, 0.0, 0.7824, +0.24)
    # test 11 toCartesian: [3980581, 97, 4966825]
    # test 12 toCartesian: [3194434.411, 3194434.411, 4487326.82]
    # test 13 toLatLon: 45.0°N, 045.0°E, +0.00m
    # test 14 Nvector: (0.5, 0.5, 0.7071)
    # test 15 Nvector: (0.5, 0.5, 0.707, +1.00)
    # all geodesy.ellipsoidalNvector tests passed (Python 2.7.13 64bit)

    # testing geodesy.ellipsoidalVincenty version 17.02.13
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
    # all geodesy.ellipsoidalVincenty tests passed (Python 2.7.13 64bit)

    # testing ellipsoidalNvector version 17.02.13
    # test 1 toLatLon: 44.995674°N, 045.0°E
    # test 2 toNvector: (0.50004, 0.50004, 0.70705)
    # test 3 equals: False
    # test 4 equals: True
    # test 5 sumOf: (52.70504, 0.61904, 0.70705)
    # test 6 sumOf: Nv
    # test 7 copy: True
    # test 8 convertDatum: 51.477284°N, 000.00002°E, -45.91m
    # test 9 convertDatum: 51.4773°N, 000.0°E, -45.91m
    # test 10 toNVector: (0.6228, 0.0, 0.7824, +0.24)
    # test 11 toCartesian: [3980581, 97, 4966825]
    # test 12 toCartesian: [3194434.411, 3194434.411, 4487326.82]
    # test 13 toLatLon: 45.0°N, 045.0°E, +0.00m
    # test 14 Nvector: (0.5, 0.5, 0.7071)
    # test 15 Nvector: (0.5, 0.5, 0.707, +1.00)
    # all ellipsoidalNvector tests passed (Python 3.6.0 64bit)

    # testing ellipsoidalVincenty version 17.02.13
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
    # all ellipsoidalVincenty tests passed (Python 3.6.0 64bit)
