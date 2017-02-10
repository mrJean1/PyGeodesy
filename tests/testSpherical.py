
# -*- coding: utf-8 -*-

# Test spherical earth model functions and methods.

__version__ = '17.02.02'

if __name__ == '__main__':

    from tests import Tests

    from geodesy import sphericalNvector as N
    t = Tests(__file__, __version__, N)
    t.testLatLon(N.LatLon)
    t.testSpherical(N.LatLon, N.Nvector)
    t.testVectorial(N.LatLon, N.Nvector, N.sumOf)
    t.results()

    from geodesy import sphericalTrigonometry as T
    t = Tests(__file__, __version__, T)
    t.testLatLon(T.LatLon)
    t.testSpherical(T.LatLon)
    t.results()
    t.exit()

    # Typical test results (on MacOS 10.12.2):

    # testing geodesy.sphericalNvector version 17.02.10
    # test 1 lat/lonDMS: 52.20472°N, 000.14056°E
    # test 2 lat/lonDMS F_DM: 52°12.283′N, 000°08.434′E
    # test 3 lat/lonDMS F_DM: 52°12.2832′N, 000°08.4336′E
    # test 4 lat/lonDMS F_DMS: 52°12′17″N, 000°08′26″E
    # test 5 lat/lonDMS F_DMS: 52°12′17.0″N, 000°08′26.0″E
    # test 6 lat/lonDMS F_RAD: 0.911144N, 0.002453E
    # test 7 equals: True
    # test 8 equals: False
    # test 9 bearingTo: 156.1666
    # test 10 finalBearingTo: 157.8904
    # test 11 bearingTo: 65.8921
    # test 12 copy: True
    # test 13 distanceTo: 404279.720589
    # test 14 distanceTo: 404279.720589
    # test 15 distanceTo: 2145
    # test 16 midpointTo: 50.536327°N, 001.274614°E
    # test 17 destination: 51.513546°N, 000.098345°W
    # test 18 destination: 51°30′49″N, 000°05′54″W
    # test 19 destination: 34°37′N, 116°33′W
    # test 20 destination: 34.613643°N, 116.551171°W
    # test 21 crossTrackDistanceTo: -305.67
    # test 22 crossTrackDistanceTo: -307.55
    # test 23 greatCircle: (-0.79408, 0.12856, 0.59406)
    # test 24 greatCircleTo: (-0.79408, 0.12859, 0.59406)
    # test 25 intermediateTo: 51.372084°N, 000.707337°E
    # test 26 intermediateChordTo: 51.372294°N, 000.707192°E
    # test 27 intersection: 50.907608°N, 004.508575°E
    # test 28 intersection: 50°54′27.39″N, 004°30′30.87″E
    # test 29 intersection: 43.5719°N, 116.188757°W
    # test 30 intersection: 43°34′18.84″N, 116°11′19.53″W
    # test 31 maxLat0: 90.0
    # test 32 maxLat1: 89.0
    # test 33 maxLat90: 0.0
    # test 34 rhumbBearingTo: 116.722
    # test 35 rhumbDistanceTo: 40307.8
    # test 36 rhumbMidpointo: 51.0455°N, 001.595727°E
    # test 37 crossTrackDistanceTo: -305.67
    # test 38 crossTrackDistanceTo: -307.55
    # test 39 toLatLon: 44.995674°N, 045.0°E
    # test 40 toNvector: (0.50004, 0.50004, 0.70705)
    # test 41 equals: False
    # test 42 equals: True
    # test 43 sumOf: (52.70504, 0.61904, 0.70705)
    # test 44 sumOf: Nv
    # test 45 copy: True
    # test 46 nearestOn: 51.0004°N, 001.9°E
    # test 47 distanceTo: 42.712
    # test 48 nearestOn: 51.0°N, 002.0°E
    # test 49 nearestOn: 00.0°N, 000.0°E
    # test 50 nearestOn: 00.0°N, 020.0°E
    # test 51 BasseC: 47.3038°N, 002.5721°W
    # test 52 BasseH: 47.311067°N, 002.528617°W
    # test 53 triangulate: 47.323667°N, 002.568501°W
    # all geodesy.sphericalNvector tests passed (Python 2.7.13 64bit)

    # testing geodesy.sphericalTrigonometry version 17.02.09
    # test 1 lat/lonDMS: 52.20472°N, 000.14056°E
    # test 2 lat/lonDMS F_DM: 52°12.283′N, 000°08.434′E
    # test 3 lat/lonDMS F_DM: 52°12.2832′N, 000°08.4336′E
    # test 4 lat/lonDMS F_DMS: 52°12′17″N, 000°08′26″E
    # test 5 lat/lonDMS F_DMS: 52°12′17.0″N, 000°08′26.0″E
    # test 6 lat/lonDMS F_RAD: 0.911144N, 0.002453E
    # test 7 equals: True
    # test 8 equals: False
    # test 9 bearingTo: 156.1666
    # test 10 finalBearingTo: 157.8904
    # test 11 bearingTo: 65.8921
    # test 12 copy: True
    # test 13 distanceTo: 404279.720589
    # test 14 distanceTo: 404279.720589
    # test 15 distanceTo: 2145
    # test 16 midpointTo: 50.536327°N, 001.274614°E
    # test 17 destination: 51.513546°N, 000.098345°W
    # test 18 destination: 51°30′49″N, 000°05′54″W
    # test 19 destination: 34°37′N, 116°33′W
    # test 20 destination: 34.613643°N, 116.551171°W
    # test 21 crossTrackDistanceTo: LatLon.crossTrackDistanceTo(end=bearing)
    # test 22 crossTrackDistanceTo: -307.55
    # test 23 greatCircle: (-0.79408, 0.12856, 0.59406)
    # test 24 intermediateTo: 51.372084°N, 000.707337°E
    # test 25 intersection: 50.907608°N, 004.508575°E
    # test 26 intersection: 50°54′27.39″N, 004°30′30.87″E
    # test 27 intersection: 43.5719°N, 116.188757°W
    # test 28 intersection: 43°34′18.84″N, 116°11′19.53″W
    # test 29 maxLat0: 90.0
    # test 30 maxLat1: 89.0
    # test 31 maxLat90: 0.0
    # test 32 crossingParallels: 009°35′38.65″E, 170°24′21.35″E
    # test 33 rhumbBearingTo: 116.722
    # test 34 rhumbDistanceTo: 40307.8
    # test 35 rhumbMidpointo: 51.0455°N, 001.595727°E
    # all geodesy.sphericalTrigonometry tests passed (Python 2.7.13 64bit)

    # testing sphericalNvector version 17.02.10
    # test 1 lat/lonDMS: 52.20472°N, 000.14056°E
    # test 2 lat/lonDMS F_DM: 52°12.283′N, 000°08.434′E
    # test 3 lat/lonDMS F_DM: 52°12.2832′N, 000°08.4336′E
    # test 4 lat/lonDMS F_DMS: 52°12′17″N, 000°08′26″E
    # test 5 lat/lonDMS F_DMS: 52°12′17.0″N, 000°08′26.0″E
    # test 6 lat/lonDMS F_RAD: 0.911144N, 0.002453E
    # test 7 equals: True
    # test 8 equals: False
    # test 9 bearingTo: 156.1666
    # test 10 finalBearingTo: 157.8904
    # test 11 bearingTo: 65.8921
    # test 12 copy: True
    # test 13 distanceTo: 404279.720589
    # test 14 distanceTo: 404279.720589
    # test 15 distanceTo: 2145
    # test 16 midpointTo: 50.536327°N, 001.274614°E
    # test 17 destination: 51.513546°N, 000.098345°W
    # test 18 destination: 51°30′49″N, 000°05′54″W
    # test 19 destination: 34°37′N, 116°33′W
    # test 20 destination: 34.613643°N, 116.551171°W
    # test 21 crossTrackDistanceTo: -305.67
    # test 22 crossTrackDistanceTo: -307.55
    # test 23 greatCircle: (-0.79408, 0.12856, 0.59406)
    # test 24 greatCircleTo: (-0.79408, 0.12859, 0.59406)
    # test 25 intermediateTo: 51.372084°N, 000.707337°E
    # test 26 intermediateChordTo: 51.372294°N, 000.707192°E
    # test 27 intersection: 50.907608°N, 004.508575°E
    # test 28 intersection: 50°54′27.39″N, 004°30′30.87″E
    # test 29 intersection: 43.5719°N, 116.188757°W
    # test 30 intersection: 43°34′18.84″N, 116°11′19.53″W
    # test 31 maxLat0: 90.0
    # test 32 maxLat1: 89.0
    # test 33 maxLat90: 0.0
    # test 34 rhumbBearingTo: 116.722
    # test 35 rhumbDistanceTo: 40307.8
    # test 36 rhumbMidpointo: 51.0455°N, 001.595727°E
    # test 37 crossTrackDistanceTo: -305.67
    # test 38 crossTrackDistanceTo: -307.55
    # test 39 toLatLon: 44.995674°N, 045.0°E
    # test 40 toNvector: (0.50004, 0.50004, 0.70705)
    # test 41 equals: False
    # test 42 equals: True
    # test 43 sumOf: (52.70504, 0.61904, 0.70705)
    # test 44 sumOf: Nv
    # test 45 copy: True
    # test 46 nearestOn: 51.0004°N, 001.9°E
    # test 47 distanceTo: 42.712
    # test 48 nearestOn: 51.0°N, 002.0°E
    # test 49 nearestOn: 00.0°N, 000.0°E
    # test 50 nearestOn: 00.0°N, 020.0°E
    # test 51 BasseC: 47.3038°N, 002.5721°W
    # test 52 BasseH: 47.311067°N, 002.528617°W
    # test 53 triangulate: 47.323667°N, 002.568501°W
    # all sphericalNvector tests passed (Python 3.6.0 64bit)

    # testing sphericalTrigonometry version 17.02.09
    # test 1 lat/lonDMS: 52.20472°N, 000.14056°E
    # test 2 lat/lonDMS F_DM: 52°12.283′N, 000°08.434′E
    # test 3 lat/lonDMS F_DM: 52°12.2832′N, 000°08.4336′E
    # test 4 lat/lonDMS F_DMS: 52°12′17″N, 000°08′26″E
    # test 5 lat/lonDMS F_DMS: 52°12′17.0″N, 000°08′26.0″E
    # test 6 lat/lonDMS F_RAD: 0.911144N, 0.002453E
    # test 7 equals: True
    # test 8 equals: False
    # test 9 bearingTo: 156.1666
    # test 10 finalBearingTo: 157.8904
    # test 11 bearingTo: 65.8921
    # test 12 copy: True
    # test 13 distanceTo: 404279.720589
    # test 14 distanceTo: 404279.720589
    # test 15 distanceTo: 2145
    # test 16 midpointTo: 50.536327°N, 001.274614°E
    # test 17 destination: 51.513546°N, 000.098345°W
    # test 18 destination: 51°30′49″N, 000°05′54″W
    # test 19 destination: 34°37′N, 116°33′W
    # test 20 destination: 34.613643°N, 116.551171°W
    # test 21 crossTrackDistanceTo: LatLon.crossTrackDistanceTo(end=bearing)
    # test 22 crossTrackDistanceTo: -307.55
    # test 23 greatCircle: (-0.79408, 0.12856, 0.59406)
    # test 24 intermediateTo: 51.372084°N, 000.707337°E
    # test 25 intersection: 50.907608°N, 004.508575°E
    # test 26 intersection: 50°54′27.39″N, 004°30′30.87″E
    # test 27 intersection: 43.5719°N, 116.188757°W
    # test 28 intersection: 43°34′18.84″N, 116°11′19.53″W
    # test 29 maxLat0: 90.0
    # test 30 maxLat1: 89.0
    # test 31 maxLat90: 0.0
    # test 32 crossingParallels: 009°35′38.65″E, 170°24′21.35″E
    # test 33 rhumbBearingTo: 116.722
    # test 34 rhumbDistanceTo: 40307.8
    # test 35 rhumbMidpointo: 51.0455°N, 001.595727°E
    # all sphericalTrigonometry tests passed (Python 3.6.0 64bit)
