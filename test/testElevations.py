
# -*- coding: utf-8 -*-

# Test elevations functions and methods.

__all__ = ('Tests',)
__version__ = '18.08.21'

from base import TestsBase

from pygeodesy import elevation2, Datums, geoidHeight2


class Tests(TestsBase):

    def testElevations(self, LatLon, datum, timeout):

        # <http://WikiPedia.org/wiki/Mount_Diablo>
        m, _ = elevation2(37.8816, -121.9142, timeout=timeout)
        self.test('elevation2', m, 1173.7)
        m, _ = geoidHeight2(37.8816, -121.9142, timeout=timeout)
        self.test('geoidHeight2', m, -31.703)

        Mount_Diablo = LatLon(37.8816, -121.9142)
        NewYork = LatLon(40.7791472, -73.9680804)
        Newport_RI = LatLon(41.49008, -71.312796)
        Cleveland_OH = LatLon(41.499498, -81.695391)
        # <http://github.com/maurycyp/vincenty> Maurycy Pietrzak
        Boston = LatLon(42.3541165, -71.0693514)
        for p, e, h in ((Mount_Diablo, 1173.7, -31.703),
                        (Boston,          2,   -27.765),
                        (Cleveland_OH,  199.1, -34.366),
                        (Newport_RI,      8.5, -30.009),
                        (NewYork,        32.7, -31.668)):
            m, _ = p.elevation2(datum=datum, timeout=timeout)
            self.test('elevation2', m, e, fmt='%s' if m is None else '%.3f')
            m, _ = p.geoidHeight2(datum=datum, timeout=timeout)
            self.test('geodHeight2', m, h, fmt='%s' if m is None else '%.3f')  # PYCHOK invoked

        m, x = elevation2(0, 0, timeout=timeout)
        self.test('elevation2', m, None)
        self.test('elevation2', x, 'elevation2(0.0, 0.0): non-CONUS')
        m, x = geoidHeight2(0, 0, timeout=timeout)
        self.test('geoidHeight2', m, None)
        # self.test('geoidHeight2', x, 'HTTP Error 403: Forbidden ...')

        non_CONUS = LatLon(51.4778, 0.0016)
        m, x = non_CONUS.elevation2(datum=datum, timeout=timeout)
        self.test('elevation2', m, None)
        self.test('elevation2', x, 'elevation2(51.4778, 0.0016): non-CONUS')
        m, _ = non_CONUS.geoidHeight2(datum=datum, timeout=timeout)
        self.test('geodHeight2', m, None)


if __name__ == '__main__':

    import sys

    t = Tests(__file__, __version__)

    if len(sys.argv) > 1:
        # avoid hitting web services too often/hard
        timeout = float(sys.argv[1])

        from pygeodesy import ellipsoidalVincenty
        LL = ellipsoidalVincenty.LatLon

        t.testElevations(LL, Datums.NAD83, timeout)
#       t.testElevations(LL, Datums.Sphere, timeout)
        t.testElevations(LL, Datums.WGS84, timeout)

        t.results()
    else:
        t.results(passed='SKIPPED', nl=0)

    t.exit()


# Typical test results (for 2 second timeout.  Web
# services often fail or just hang when hit too often).
#
# % pypy test/testElevations.py 2
#
#     testing testElevations.py 18.08.21
#     test 1 elevation2: 1173.7
#     test 2 geoidHeight2: -31.703
#     test 3 elevation2: 1173.700
#     test 4 geodHeight2: -31.703
#     test 5 elevation2: 2.000
#     test 6 geodHeight2: -27.765
#     test 7 elevation2: 199.100
#     test 8 geodHeight2: -34.366
#     test 9 elevation2: 8.500
#     test 10 geodHeight2: -30.009
#     test 11 elevation2: 32.700
#     test 12 geodHeight2: -31.668
#     test 13 elevation2: None
#     test 14 elevation2: elevation2(0.0, 0.0): non-CONUS
#     test 15 geoidHeight2: None
#     test 16 elevation2: None
#     test 17 elevation2: elevation2(51.4778, 0.0016): non-CONUS
#     test 18 geodHeight2: None
#     test 19 elevation2: 1173.7
#     test 20 geoidHeight2: -31.703
#     test 21 elevation2: 1173.700
#     test 22 geodHeight2: -31.703
#     test 23 elevation2: 2.000
#     test 24 geodHeight2: -27.765
#     test 25 elevation2: 199.100
#     test 26 geodHeight2: -34.366
#     test 27 elevation2: 8.500
#     test 28 geodHeight2: -30.009
#     test 29 elevation2: 32.700
#     test 30 geodHeight2: -31.668
#     test 31 elevation2: None
#     test 32 elevation2: elevation2(0.0, 0.0): non-CONUS
#     test 33 geoidHeight2: None
#     test 34 elevation2: None
#     test 35 elevation2: elevation2(51.4778, 0.0016): non-CONUS
#     test 36 geodHeight2: None
#
#     all testElevations.py tests passed (PyGeodesy 18.8.21 PyPy-Python 2.7.13 64bit macOS 10.13.6) 6.512 sec
#
# % python2.7 test/testElevations.py 2
#
#     testing testElevations.py 18.08.21
#     test 1 elevation2: 1173.7
#     test 2 geoidHeight2: -31.703
#     test 3 elevation2: 1173.700
#     test 4 geodHeight2: -31.703
#     test 5 elevation2: 2.000
#     test 6 geodHeight2: -27.765
#     test 7 elevation2: 199.100
#     test 8 geodHeight2: -34.366
#     test 9 elevation2: 8.500
#     test 10 geodHeight2: -30.009
#     test 11 elevation2: 32.700
#     test 12 geodHeight2: -31.668
#     test 13 elevation2: None
#     test 14 elevation2: elevation2(0.0, 0.0): non-CONUS
#     test 15 geoidHeight2: None
#     test 16 elevation2: None
#     test 17 elevation2: elevation2(51.4778, 0.0016): non-CONUS
#     test 18 geodHeight2: None
#     test 19 elevation2: 1173.7
#     test 20 geoidHeight2: -31.703
#     test 21 elevation2: 1173.700
#     test 22 geodHeight2: -31.703
#     test 23 elevation2: 2.000
#     test 24 geodHeight2: -27.765
#     test 25 elevation2: 199.100
#     test 26 geodHeight2: -34.366
#     test 27 elevation2: 8.500
#     test 28 geodHeight2: -30.009
#     test 29 elevation2: 32.700
#     test 30 geodHeight2: -31.668
#     test 31 elevation2: None
#     test 32 elevation2: elevation2(0.0, 0.0): non-CONUS
#     test 33 geoidHeight2: None
#     test 34 elevation2: None
#     test 35 elevation2: elevation2(51.4778, 0.0016): non-CONUS
#     test 36 geodHeight2: None
#
#     all testElevations.py tests passed (PyGeodesy 18.8.21 Python 2.7.15 64bit geographiclib 1.49 numpy 1.14.0 macOS 10.13.6) 6.058 sec
#
# % intelpython3 test/testElevations.py 2
#
#     testing testElevations.py 18.08.21
#     test 1 elevation2: 1173.7
#     test 2 geoidHeight2: -31.703
#     test 3 elevation2: 1173.700
#     test 4 geodHeight2: -31.703
#     test 5 elevation2: 2.000
#     test 6 geodHeight2: -27.765
#     test 7 elevation2: 199.100
#     test 8 geodHeight2: -34.366
#     test 9 elevation2: 8.500
#     test 10 geodHeight2: -30.009
#     test 11 elevation2: 32.700
#     test 12 geodHeight2: -31.668
#     test 13 elevation2: None
#     test 14 elevation2: elevation2(0.0, 0.0): non-CONUS
#     test 15 geoidHeight2: None
#     test 16 elevation2: None
#     test 17 elevation2: elevation2(51.4778, 0.0016): non-CONUS
#     test 18 geodHeight2: None
#     test 19 elevation2: 1173.7
#     test 20 geoidHeight2: -31.703
#     test 21 elevation2: 1173.700
#     test 22 geodHeight2: -31.703
#     test 23 elevation2: 2.000
#     test 24 geodHeight2: -27.765
#     test 25 elevation2: 199.100
#     test 26 geodHeight2: -34.366
#     test 27 elevation2: 8.500
#     test 28 geodHeight2: -30.009
#     test 29 elevation2: 32.700
#     test 30 geodHeight2: -31.668
#     test 31 elevation2: None
#     test 32 elevation2: elevation2(0.0, 0.0): non-CONUS
#     test 33 geoidHeight2: None
#     test 34 elevation2: None
#     test 35 elevation2: elevation2(51.4778, 0.0016): non-CONUS
#     test 36 geodHeight2: None
#
#     all testElevations.py tests passed (PyGeodesy 18.8.21 Intel-Python 3.5.3 64bit numpy 1.11.3 macOS 10.13.6) 6.131 sec
#
# % python3 test/testElevations.py 2
#
#     testing testElevations.py 18.08.21
#     test 1 elevation2: 1173.7
#     test 2 geoidHeight2: -31.703
#     test 3 elevation2: 1173.700
#     test 4 geodHeight2: -31.703
#     test 5 elevation2: 2.000
#     test 6 geodHeight2: -27.765
#     test 7 elevation2: 199.100
#     test 8 geodHeight2: -34.366
#     test 9 elevation2: 8.500
#     test 10 geodHeight2: -30.009
#     test 11 elevation2: 32.700
#     test 12 geodHeight2: -31.668
#     test 13 elevation2: None
#     test 14 elevation2: elevation2(0.0, 0.0): non-CONUS
#     test 15 geoidHeight2: None
#     test 16 elevation2: None
#     test 17 elevation2: elevation2(51.4778, 0.0016): non-CONUS
#     test 18 geodHeight2: None
#     test 19 elevation2: 1173.7
#     test 20 geoidHeight2: -31.703
#     test 21 elevation2: 1173.700
#     test 22 geodHeight2: -31.703
#     test 23 elevation2: 2.000
#     test 24 geodHeight2: -27.765
#     test 25 elevation2: 199.100
#     test 26 geodHeight2: -34.366
#     test 27 elevation2: 8.500
#     test 28 geodHeight2: -30.009
#     test 29 elevation2: 32.700
#     test 30 geodHeight2: -31.668
#     test 31 elevation2: None
#     test 32 elevation2: elevation2(0.0, 0.0): non-CONUS
#     test 33 geoidHeight2: None
#     test 34 elevation2: None
#     test 35 elevation2: elevation2(51.4778, 0.0016): non-CONUS
#     test 36 geodHeight2: None
#
#     all testElevations.py tests passed (PyGeodesy 18.8.21 Python 3.7.0 64bit macOS 10.13.6) 8.430 sec
