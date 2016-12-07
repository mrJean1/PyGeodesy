
# -*- coding: utf-8 -*-

# Test datums, ellipsoids and transforms.

__version__ = '16.12.07'

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
    # test 6 WGS84: A=6367449.145823415, e=0.0818191908, f=1/298.2572235630, n=0.0016792204(-3.7914875232e-13)
    # test 7 WGS84: Alpha6=(0, 0.0008377318206244698, 7.608527773572307e-07, 1.1976455033294527e-09, 2.4291706072013587e-12, 5.711757677865804e-15, 1.4911177312583895e-17)
    # test 8 WGS84: Beta6=(0, 0.0008377321640579486, 5.905870152220203e-08, 1.6734826652839968e-10, 2.1647980400627059e-13, 3.7879780461686053e-16, 7.2487488906941545e-19)
    # all geodesy.datum tests passed (Python 2.7.10)

    # testing geodesy.datum version 16.12.06
    # test 1 ellipsoid: True
    # test 2 transform: True
    # test 3 datum: True
    # test 4 meanR: 6371008.77142
    # test 5 WGS84: a=6378137.0, b=6356752.3142499998, f=0.0033528107, e2=0.00669438, e22=0.0067394967, R=6371008.7714166669, Rm=6367435.6797186071, name='WGS84'
    # test 6 WGS84: A=6367449.145823415, e=0.0818191908, f=1/298.2572235630, n=0.0016792204(-3.7914875232e-13)
    # test 7 WGS84: Alpha6=(0, 0.0008377318206244698, 7.608527773572307e-07, 1.1976455033294527e-09, 2.4291706072013587e-12, 5.711757677865804e-15, 1.4911177312583895e-17)
    # test 8 WGS84: Beta6=(0, 0.0008377321640579486, 5.905870152220203e-08, 1.6734826652839968e-10, 2.1647980400627059e-13, 3.7879780461686053e-16, 7.2487488906941545e-19)
    # all geodesy.datum tests passed (Python 3.5.2)
