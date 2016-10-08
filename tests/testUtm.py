
# -*- coding: utf-8 -*-

# Test spherical earth model functions and methods.

__version__ = '16.10.07'

if __name__ == '__main__':

    from tests import Tests

    from geodesy import ellipsoidalVincenty as eV, utm

    t = Tests(__file__, __version__, utm)
    t.testUtm(utm, eV.LatLon)
    t.results()

    # Typical test results (on MacOS X):

    # testing geodesy.utm version 16.10.08
    # test 1 Utm1: 03 N 448251.0 5411932.0001
    # test 2 Utm2: 31 N 448252 5411933
    # test 3 Utm2: 31 N 448251.795 5411932.678
    # test 4 Utm2: 31 N 448251.8 5411932.7 n/a n/a
    # test 5 Utm.toLatLon1: 48.8582°N, 002.2945°E
    # test 6 Utm.toLatLon1: 48°51′29.52″N, 002°17′40.2″E
    # test 7 toUtm1: 31U N 448252 5411933
    # test 8 toUtm1: 31U N 448251.795 5411932.678
    # test 9 toUtm2: [Z:31U, H:N, E:448252, N:5411933, C:-000.53131221°, S:0.9996329]
    # test 10 toUtm4: 48P N 377302 1483035
    # test 11 toUtm5: 48P N 377302.354183 1483034.777084 -000.26291348° 0.99978623
    # test 12 toUtm6: 13L S 622698 8516965
    # test 13 toUtm7: 13L S 622697.645817 8516965.222916 -000.26291348° 0.99978623
    # test 14 toUtm(61.44°N, 025.4°E): 35V N 414668 6812845
    # test 15 toUtm(47.04°S, 073.48°W): 18G S 615472 4789270
    # test 16 toUtm(40.4°N, 074.7°W): 18T N 525458 4472198
    # test 17 toUtm(44.5°N, 088.5°W): 16T N 380753 4928503
    # test 18 toUtm(50.8694°N, 115.6508°W): 11U N 594937 5636169
    # test 19 toUtm(00.0°N, 000.0°E): 31N N 166021 0
    # test 20 toUtm(00.13°N, 000.2324°W): 30N N 808084 14386
    # test 21 toUtm(45.6456°S, 023.3545°E): 34G S 683474 4942631
    # test 22 toUtm(12.765°S, 033.8765°W): 25L S 404859 8588691
    # test 23 toUtm(80.5434°S, 170.654°W): 02C S 506346 1057743
    # test 24 toUtm(90.0°N, 177.0°E): 60Z N 500000 9997965
    # test 25 toUtm(90.0°S, 177.0°W): 01A S 500000 2035
    # test 26 toUtm(90.0°N, 003.0°E): 31Z N 500000 9997965
    # test 27 toUtm(23.4578°N, 135.4545°W): 08Q N 453580 2594273
    # test 28 toUtm(77.345°N, 156.9876°E): 57X N 450794 8586116
    # test 29 toUtm(89.3454°S, 048.9306°W): 22A S 502639 75073
    # all geodesy.utm tests passed (Python 2.7.10)

    # testing geodesy.utm version 16.10.08
    # test 1 Utm1: 03 N 448251.0 5411932.0001
    # test 2 Utm2: 31 N 448252 5411933
    # test 3 Utm2: 31 N 448251.795 5411932.678
    # test 4 Utm2: 31 N 448251.8 5411932.7 n/a n/a
    # test 5 Utm.toLatLon1: 48.8582°N, 002.2945°E
    # test 6 Utm.toLatLon1: 48°51′29.52″N, 002°17′40.2″E
    # test 7 toUtm1: 31U N 448252 5411933
    # test 8 toUtm1: 31U N 448251.795 5411932.678
    # test 9 toUtm2: [Z:31U, H:N, E:448252, N:5411933, C:-000.53131221°, S:0.9996329]
    # test 10 toUtm4: 48P N 377302 1483035
    # test 11 toUtm5: 48P N 377302.354183 1483034.777084 -000.26291348° 0.99978623
    # test 12 toUtm6: 13L S 622698 8516965
    # test 13 toUtm7: 13L S 622697.645817 8516965.222916 -000.26291348° 0.99978623
    # test 14 toUtm(61.44°N, 025.4°E): 35V N 414668 6812845
    # test 15 toUtm(47.04°S, 073.48°W): 18G S 615472 4789270
    # test 16 toUtm(40.4°N, 074.7°W): 18T N 525458 4472198
    # test 17 toUtm(44.5°N, 088.5°W): 16T N 380753 4928503
    # test 18 toUtm(50.8694°N, 115.6508°W): 11U N 594937 5636169
    # test 19 toUtm(00.0°N, 000.0°E): 31N N 166021 0
    # test 20 toUtm(00.13°N, 000.2324°W): 30N N 808084 14386
    # test 21 toUtm(45.6456°S, 023.3545°E): 34G S 683474 4942631
    # test 22 toUtm(12.765°S, 033.8765°W): 25L S 404859 8588691
    # test 23 toUtm(80.5434°S, 170.654°W): 02C S 506346 1057743
    # test 24 toUtm(90.0°N, 177.0°E): 60Z N 500000 9997965
    # test 25 toUtm(90.0°S, 177.0°W): 01A S 500000 2035
    # test 26 toUtm(90.0°N, 003.0°E): 31Z N 500000 9997965
    # test 27 toUtm(23.4578°N, 135.4545°W): 08Q N 453580 2594273
    # test 28 toUtm(77.345°N, 156.9876°E): 57X N 450794 8586116
    # test 29 toUtm(89.3454°S, 048.9306°W): 22A S 502639 75073
    # all geodesy.utm tests passed (Python 3.5.1)
