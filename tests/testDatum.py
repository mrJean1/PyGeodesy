
# -*- coding: utf-8 -*-

# Test datums, ellipsoids and transforms.

__version__ = '17.01.14'

if __name__ == '__main__':

    from tests import Tests

    import geodesy
    import geodesy.datum as D

    t = Tests(__file__, __version__, D)
    t.testDatum(geodesy)
    t.results()
    t.exit()

    # Typical test results (on MacOS X):

    # testing geodesy.datum version 16.12.06
    # test 1 ellipsoid: True
    # test 2 transform: True
    # test 3 datum: True
    # test 4 meanR: 6371008.77142
    # test 5 WGS84: a=6378137.0, b=6356752.3142499998, f=0.0033528107, e2=0.00669438, e22=0.0067394967, R=6371008.7714166669, Rm=6367435.6797186071, name='WGS84'
    # test 6 WGS84: A=6367449.1458234154, e=0.0818191908, f=1/298.2572235630, n=0.0016792204(-3.7914875232e-13)
    # test 7 WGS84: Alpha6=(0, 8.377318206245e-04, 7.608527773572e-07, 1.197645503329e-09, 2.429170607201e-12, 5.711757677866e-15, 1.491117731258e-17)
    # test 8 WGS84: Beta6=(0, 8.377321640579e-04, 5.905870152220e-08, 1.673482665284e-1, 2.164798040063e-13, 3.787978046169e-16, 7.248748890694e-19)
    # all geodesy.datum tests passed (Python 2.7.13 64bit)

    # testing geodesy.datum version 16.12.06
    # test 1 ellipsoid: True
    # test 2 transform: True
    # test 3 datum: True
    # test 4 meanR: 6371008.77142
    # test 5 WGS84: a=6378137.0, b=6356752.3142499998, f=0.0033528107, e2=0.00669438, e22=0.0067394967, R=6371008.7714166669, Rm=6367435.6797186071, name='WGS84'
    # test 6 WGS84: A=6367449.1458234154, e=0.0818191908, f=1/298.2572235630, n=0.0016792204(-3.7914875232e-13)
    # test 7 WGS84: Alpha6=(0, 8.377318206245e-04, 7.608527773572e-07, 1.197645503329e-09, 2.429170607201e-12, 5.711757677866e-15, 1.491117731258e-17)
    # test 8 WGS84: Beta6=(0, 8.377321640579e-04, 5.905870152220e-08, 1.673482665284e-1, 2.164798040063e-13, 3.787978046169e-16, 7.248748890694e-19)
    # all geodesy.datum tests passed (Python 3.6.0 64bit)
