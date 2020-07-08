
# -*- coding: utf-8 -*-

# Test elevations functions.

__all__ = ('Tests',)
__version__ = '20.07.03'

from base import isPython2, isPython3, TestsBase

from pygeodesy import elevation2, Datums, geoidHeight2


class Tests(TestsBase):

    def testApprox(self, name, m, e, x):
        # allow margin and errors
        if m:
            self.test(name, m, e, fmt='%.3f', known=abs(m - e) < 0.1)
        else:  # m is None
            self.test(name, x, e, fmt='%s', known=True)

    def testElevations(self, LatLon, datum, timeout):

        # <https://WikiPedia.org/wiki/Mount_Diablo>
        m, x = elevation2(37.8816, -121.9142, timeout=timeout)
        self.testApprox('elevation2', m, 1173.79, x)
        m, x = geoidHeight2(37.8816, -121.9142, timeout=timeout)
        self.testApprox('geoidHeight2', m, -31.699, x)

        MtDiablo = LatLon(37.8816, -121.9142)
        NewYork = LatLon(40.7791472, -73.9680804)
        Newport_RI = LatLon(41.49008, -71.312796)
        Cleveland_OH = LatLon(41.499498, -81.695391)
        # <https://GitHub.com/maurycyp/vincenty> Maurycy Pietrzak
        Boston = LatLon(42.3541165, -71.0693514)
        for p, e, h in ((MtDiablo,    1173.79, -31.699),
                        (Boston,         2.03, -27.773),
                        (Cleveland_OH, 199.18, -34.337),
                        (Newport_RI,     8.52, -30.000),
                        (NewYork,       32.79, -31.666)):
            m, x = p.elevation2(datum=datum, timeout=timeout)
            self.testApprox('elevation2', m, e, x)
            m, x = p.geoidHeight2(datum=datum, timeout=timeout)
            self.testApprox('geodHeight2', m, h, x)  # PYCHOK test attr?

        m, x = elevation2(0, 0, timeout=timeout)
        self.testError('elevation2', m, x, 'non-CONUS -1000000.00')
        m, x = geoidHeight2(0, 0, timeout=timeout)
        self.testError('geoidHeight2', m, x, "<HTTPError 403: 'Forbidden'>" if isPython3
                                        else ('HTTPError()' if isPython2 else 'isPython?'))

        non_CONUS = LatLon(51.4, 0.1)
        m, x = non_CONUS.elevation2(datum=datum, timeout=0.01)  # force timeout
        self.testError('elevation2', m, x, "URLError(timeout('timed out'))")
        m, x = non_CONUS.geoidHeight2(datum=datum, timeout=0.01)  # force timeout
        self.testError('geodHeight2', m, x, "URLError(timeout('timed out'))")

    def testError(self, name, m, x, error):
        # check an error test result
        x = x.replace(',)', ')').split('): ')[-1]
        self.test(name, (m, x), (None, error), known=True)


if __name__ == '__main__':

    import os
    import sys

    t = Tests(__file__, __version__)

    if len(sys.argv) > 1:
        timeout = float(sys.argv[1])
    else:
        timeout = float(os.environ.get('PYGEODESY_COVERAGE', '0'))

    if timeout > 4:
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

# python test/testElevations.py  2
#
#   testing testElevations.py 18.08.25
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
#   test 13 elevation2: (None, "ValueError('non-CONUS')")
#   test 14 geoidHeight2: (None, 'HTTPError()')
#   test 15 elevation2: (None, "URLError(timeout('timed out'))")
#   test 16 geodHeight2: (None, "URLError(timeout('timed out'))")
#   test 17 elevation2: 1173.79
#   test 18 geoidHeight2: -31.703
#   test 19 elevation2: 1173.790
#   test 20 geodHeight2: -31.703
#   test 21 elevation2: 2.030
#   test 22 geodHeight2: -27.765
#   test 23 elevation2: 199.180
#   test 24 geodHeight2: -34.366
#   test 25 elevation2: 8.520
#   test 26 geodHeight2: -30.009
#   test 27 elevation2: 32.790
#   test 28 geodHeight2: -31.668
#   test 29 elevation2: (None, "ValueError('non-CONUS')")
#   test 30 geoidHeight2: (None, 'HTTPError()')
#   test 31 elevation2: (None, "URLError(timeout('timed out'))")
#   test 32 geodHeight2: (None, "URLError(timeout('timed out'))")

#   all testElevations.py tests passed (PyGeodesy 18.8.24 Python 2.7.15 64bit geographiclib 1.49 numpy 1.14.0 macOS 10.13.6) 5.757 sec


# pypy test/testElevations.py 2
#
#   testing testElevations.py 18.08.25
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
#   test 13 elevation2: (None, "ValueError('non-CONUS')")
#   test 14 geoidHeight2: (None, 'HTTPError()')
#   test 15 elevation2: (None, "URLError(timeout('timed out'))")
#   test 16 geodHeight2: (None, "URLError(timeout('timed out'))")
#   test 17 elevation2: 1173.79
#   test 18 geoidHeight2: -31.703
#   test 19 elevation2: 1173.790
#   test 20 geodHeight2: -31.703
#   test 21 elevation2: 2.030
#   test 22 geodHeight2: -27.765
#   test 23 elevation2: 199.180
#   test 24 geodHeight2: -34.366
#   test 25 elevation2: 8.520
#   test 26 geodHeight2: -30.009
#   test 27 elevation2: 32.790
#   test 28 geodHeight2: -31.668
#   test 29 elevation2: (None, "ValueError('non-CONUS')")
#   test 30 geoidHeight2: (None, 'HTTPError()')
#   test 31 elevation2: (None, "URLError(timeout('timed out'))")
#   test 32 geodHeight2: (None, "URLError(timeout('timed out'))")
#
#   all testElevations.py tests passed (PyGeodesy 18.8.26 PyPy-Python 2.7.13 64bit macOS 10.13.6) 7.180 sec


# python3 test/testElevations.py 2
#
#   testing testElevations.py 18.08.25
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
#   test 13 elevation2: (None, "ValueError('non-CONUS')")
#   test 14 geoidHeight2: (None, "<HTTPError 403: 'Forbidden'>")
#   test 15 elevation2: (None, "URLError(timeout('timed out'))")
#   test 16 geodHeight2: (None, "URLError(timeout('timed out'))")
#   test 17 elevation2: 1173.79
#   test 18 geoidHeight2: -31.703
#   test 19 elevation2: 1173.790
#   test 20 geodHeight2: -31.703
#   test 21 elevation2: 2.030
#   test 22 geodHeight2: -27.765
#   test 23 elevation2: 199.180
#   test 24 geodHeight2: -34.366
#   test 25 elevation2: 8.520
#   test 26 geodHeight2: -30.009
#   test 27 elevation2: 32.790
#   test 28 geodHeight2: -31.668
#   test 29 elevation2: (None, "ValueError('non-CONUS')")
#   test 30 geoidHeight2: (None, "<HTTPError 403: 'Forbidden'>")
#   test 31 elevation2: (None, "URLError(timeout('timed out'))")
#   test 32 geodHeight2: (None, "URLError(timeout('timed out'))")
#
#   all testElevations.py tests passed (PyGeodesy 18.8.26 Python 3.7.0 64bit macOS 10.13.6) 8.479 sec


# intelpython3 test/testElevations.py 2
#
#   testing testElevations.py 18.08.25
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
#   test 13 elevation2: (None, "ValueError('non-CONUS')")
#   test 14 geoidHeight2: (None, "<HTTPError 403: 'Forbidden'>")
#   test 15 elevation2: (None, "URLError(timeout('timed out'))")
#   test 16 geodHeight2: (None, "URLError(timeout('timed out'))")
#   test 17 elevation2: 1173.79
#   test 18 geoidHeight2: -31.703
#   test 19 elevation2: 1173.790
#   test 20 geodHeight2: -31.703
#   test 21 elevation2: 2.030
#   test 22 geodHeight2: -27.765
#   test 23 elevation2: 199.180
#   test 24 geodHeight2: -34.366
#   test 25 elevation2: 8.520
#   test 26 geodHeight2: -30.009
#   test 27 elevation2: 32.790
#   test 28 geodHeight2: -31.668
#   test 29 elevation2: (None, "ValueError('non-CONUS')")
#   test 30 geoidHeight2: (None, "<HTTPError 403: 'Forbidden'>")
#   test 31 elevation2: (None, "URLError(timeout('timed out'))")
#   test 32 geodHeight2: (None, "URLError(timeout('timed out'))")
#
#   all testElevations.py tests passed (PyGeodesy 18.8.26 Intel-Python 3.5.3 64bit numpy 1.11.3 macOS 10.13.6) 7.794 sec


# Pythonista3: python2.7 test/testElevations.py 2

#   testing testElevations.py 18.08.25
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
#   test 13 elevation2: (None, "ValueError('non-CONUS')")
#   test 14 geoidHeight2: (None, 'HTTPError()')
#   test 15 elevation2: (None, "URLError(timeout('timed out'))")
#   test 16 geodHeight2: (None, "URLError(timeout('timed out'))")
#   test 17 elevation2: 1173.79
#   test 18 geoidHeight2: -31.703
#   test 19 elevation2: 1173.790
#   test 20 geodHeight2: -31.703
#   test 21 elevation2: 2.030
#   test 22 geodHeight2: -27.765
#   test 23 elevation2: 199.180
#   test 24 geodHeight2: -34.366
#   test 25 elevation2: 8.520
#   test 26 geodHeight2: -30.009
#   test 27 elevation2: 32.790
#   test 28 geodHeight2: -31.668
#   test 29 elevation2: (None, "ValueError('non-CONUS')")
#   test 30 geoidHeight2: (None, 'HTTPError()')
#   test 31 elevation2: (None, "URLError(timeout('timed out'))")
#   test 32 geodHeight2: (None, "URLError(timeout('timed out'))")

#   all testElevations.py tests passed (PyGeodesy 18.8.28 Python 2.7.12 64bit numpy 1.8.0 iOS 11.4.1) 11.339 sec


# Pythonista3: python3.6 test/testElevations.py 2

#   testing testElevations.py 18.08.25
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
#   test 13 elevation2: (None, "ValueError('non-CONUS')")
#   test 14 geoidHeight2: (None, "<HTTPError 403: 'Forbidden'>")
#   test 15 elevation2: (None, "URLError(timeout('timed out'))")
#   test 16 geodHeight2: (None, "URLError(timeout('timed out'))")
#   test 17 elevation2: 1173.79
#   test 18 geoidHeight2: -31.703
#   test 19 elevation2: 1173.790
#   test 20 geodHeight2: -31.703
#   test 21 elevation2: 2.030
#   test 22 geodHeight2: -27.765
#   test 23 elevation2: 199.180
#   test 24 geodHeight2: -34.366
#   test 25 elevation2: 8.520
#   test 26 geodHeight2: -30.009
#   test 27 elevation2: 32.790
#   test 28 geodHeight2: -31.668
#   test 29 elevation2: (None, "ValueError('non-CONUS')")
#   test 30 geoidHeight2: (None, "<HTTPError 403: 'Forbidden'>")
#   test 31 elevation2: (None, "URLError(timeout('timed out'))")
#   test 32 geodHeight2: (None, "URLError(timeout('timed out'))")

#   all testElevations.py tests passed (PyGeodesy 18.8.28 Python 3.6.1 64bit numpy 1.8.0 iOS 11.4.1) 7.694 sec
