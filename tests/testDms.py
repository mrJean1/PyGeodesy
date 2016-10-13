
# -*- coding: utf-8 -*-

# Test degrees, minutes, seconds functions.

__version__ = '16.10.13'

if __name__ == '__main__':

    from tests import Tests

    from geodesy import dms
    t = Tests(__file__, __version__, dms)
    t.testDMS()
    t.results()

    # Typical test results (on MacOS X):

    # testing geodesy.dms version 16.10.10
    # test 1 parseDMS: 0.0
    # test 2 parseDMS: 0.0
    # test 3 parseDMS: 0.0
    # test 4 parseDMS: 0.0
    # test 5 parseDMS: 0.0
    # test 6 parseDMS: 0.0
    # test 7 parse3llh: 51.477811, -0.001475, 0.000000
    # test 8 toDMS: 45°45′45.36″
    # test 9 toDMS: 45.7626°
    # test 10 toDMS: 45°45.756′
    # test 11 toDMS: 45°45′45.36″
    # test 12 toDMS: 45.7626°
    # test 13 toDMS: 45°45.7560′
    # test 14 toDMS: 45°45′45.36″
    # test 15 compassPoint: N
    # test 16 compassPoint: N
    # test 17 compassPoint: N
    # test 18 compassPoint: N
    # test 19 compassPoint: NNE
    # test 20 compassPoint: N
    # test 21 compassPoint: NE
    # test 22 compassPoint: NNE
    # test 23 compassPoint: SW
    # test 24 compassPoint: W
    # test 25 compassPoint: SW
    # test 26 compassPoint: SW
    # test 27 compassPoint: WSW
    # test 28 compassPoint: W
    # test 29 compassPoint: SW
    # test 30 compassPoint: WSW
    # all geodesy.dms tests passed (Python 2.7.10)

    # testing geodesy.dms version 16.10.10
    # test 1 parseDMS: 0.0
    # test 2 parseDMS: 0.0
    # test 3 parseDMS: 0.0
    # test 4 parseDMS: 0.0
    # test 5 parseDMS: 0.0
    # test 6 parseDMS: 0.0
    # test 7 parse3llh: 51.477811, -0.001475, 0.000000
    # test 8 toDMS: 45°45′45.36″
    # test 9 toDMS: 45.7626°
    # test 10 toDMS: 45°45.756′
    # test 11 toDMS: 45°45′45.36″
    # test 12 toDMS: 45.7626°
    # test 13 toDMS: 45°45.7560′
    # test 14 toDMS: 45°45′45.36″
    # test 15 compassPoint: N
    # test 16 compassPoint: N
    # test 17 compassPoint: N
    # test 18 compassPoint: N
    # test 19 compassPoint: NNE
    # test 20 compassPoint: N
    # test 21 compassPoint: NE
    # test 22 compassPoint: NNE
    # test 23 compassPoint: SW
    # test 24 compassPoint: W
    # test 25 compassPoint: SW
    # test 26 compassPoint: SW
    # test 27 compassPoint: WSW
    # test 28 compassPoint: W
    # test 29 compassPoint: SW
    # test 30 compassPoint: WSW
    # all geodesy.dms tests passed (Python 3.5.1)
