
# -*- coding: utf-8 -*-

# Test elevations functions and methods.

__all__ = ('Tests',)
__version__ = '18.08.22'

from base import TestsBase

from pygeodesy import elevation2, Datums, geoidHeight2


class Tests(TestsBase):

    def testElevations(self, LatLon, datum, timeout):

        # <http://WikiPedia.org/wiki/Mount_Diablo>
        m, _ = elevation2(37.8816, -121.9142, timeout=timeout)
        self.test('elevation2', m, 1173.79)
        m, _ = geoidHeight2(37.8816, -121.9142, timeout=timeout)
        self.test('geoidHeight2', m, -31.703)

        MtDiablo = LatLon(37.8816, -121.9142)
        NewYork = LatLon(40.7791472, -73.9680804)
        Newport_RI = LatLon(41.49008, -71.312796)
        Cleveland_OH = LatLon(41.499498, -81.695391)
        # <http://github.com/maurycyp/vincenty> Maurycy Pietrzak
        Boston = LatLon(42.3541165, -71.0693514)
        for p, e, h in ((MtDiablo,    1173.79, -31.703),
                        (Boston,         2.03, -27.765),
                        (Cleveland_OH, 199.18, -34.366),
                        (Newport_RI,     8.52, -30.009),
                        (NewYork,       32.79, -31.668)):
            m, _ = p.elevation2(datum=datum, timeout=timeout)
            self.test('elevation2', m, e, fmt='%s' if m is None else '%.3f')
            m, _ = p.geoidHeight2(datum=datum, timeout=timeout)
            self.test('geodHeight2', m, h, fmt='%s' if m is None else '%.3f')  # PYCHOK invoked

        m, x = elevation2(0, 0, timeout=timeout)
        self.test('elevation2', m, None)
        self.test('elevation2', x, "elevation2(0.0, 0.0): ValueError('non-CONUS',)")
        m, x = geoidHeight2(0, 0, timeout=timeout)
        self.test('geoidHeight2', m, None)
        self.test('geoidHeight2', x, x)  # 'HTTP Error 403: Forbidden ...'

        non_CONUS = LatLon(51.4778, 0.0016)
        m, x = non_CONUS.elevation2(datum=datum, timeout=0.01)  # force timeout
        self.test('elevation2', m, None)
        self.test('elevation2', x, x)
        m, x = non_CONUS.geoidHeight2(datum=datum, timeout=0.01)  # force timeout
        self.test('geodHeight2', m, None)
        self.test('geoidHeight2', x, x)  # 'HTTP Error 403: Forbidden ...'


if __name__ == '__main__':

    import sys

    t = Tests(__file__, __version__)

    if len(sys.argv) > 1:
        timeout = float(sys.argv[1])

        from pygeodesy import ellipsoidalVincenty
        LL = ellipsoidalVincenty.LatLon

        t.testElevations(LL, Datums.NAD83, timeout)
#       t.testElevations(LL, Datums.Sphere, timeout)
        t.testElevations(LL, Datums.WGS84, timeout)

        t.results()
    else:
        # sites fail/hang when hit too often
        t.results(passed='SKIPPED', nl=0)

    t.exit()


# Typical test results (for 2 second timeout).
#
# python test/testElevations.py 2
#
#   testing testElevations.py 18.08.22
#   test 1 elevation2: 1173.79
#   test 2 geoidHeight2: -31.703
#   test 3 elevation2: 1173.790
#   test 4 geodHeight2: -31.703
#   test 5 elevation2: 2.030
#   test 6 geodHeight2: -27.765
#   test 7 elevation2: 199.180
#   test 8 geodHeight2: -34.366
#   test 9 elevation2: 8.520
#   test 10 geodHeight2: -30.009
#   test 11 elevation2: 32.790
#   test 12 geodHeight2: -31.668
#   test 13 elevation2: None
#   test 14 elevation2: elevation2(0.0, 0.0): ValueError('non-CONUS',)
#   test 15 geoidHeight2: None
#   test 16 geoidHeight2: geoidHeight2(0.0, 0.0): <HTTPError 403: 'Forbidden'>
#   test 17 elevation2: None
#   test 18 elevation2: elevation2(51.4778, 0.0016): URLError(timeout('timed out',),)
#   test 19 geodHeight2: None
#   test 20 geoidHeight2: geoidHeight2(51.4778, 0.0016): URLError(timeout('timed out',),)
#   test 21 elevation2: 1173.79
#   test 22 geoidHeight2: -31.703
#   test 23 elevation2: 1173.790
#   test 24 geodHeight2: -31.703
#   test 25 elevation2: 2.030
#   test 26 geodHeight2: -27.765
#   test 27 elevation2: 199.180
#   test 28 geodHeight2: -34.366
#   test 29 elevation2: 8.520
#   test 30 geodHeight2: -30.009
#   test 31 elevation2: 32.790
#   test 32 geodHeight2: -31.668
#   test 33 elevation2: None
#   test 34 elevation2: elevation2(0.0, 0.0): ValueError('non-CONUS',)
#   test 35 geoidHeight2: None
#   test 36 geoidHeight2: geoidHeight2(0.0, 0.0): <HTTPError 403: 'Forbidden'>
#   test 37 elevation2: None
#   test 38 elevation2: elevation2(51.4778, 0.0016): URLError(timeout('timed out',),)
#   test 39 geodHeight2: None
#   test 40 geoidHeight2: geoidHeight2(51.4778, 0.0016): URLError(timeout('timed out',),)
#
#   all testElevations.py tests passed (PyGeodesy 18.8.22 Python 3.6.1 64bit numpy 1.8.0 iOS 11.4.1) 9.439 sec
