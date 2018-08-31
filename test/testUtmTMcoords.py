
# -*- coding: utf-8 -*-

u'''Test UTM functions with the C(TMcoords.dat} from
U{Test data for the transverse Mercator projection
<http://geographiclib.sourceforge.io/html/transversemercator.html>},
also available U{here<http://zenodo.org/record/32470>}.
'''

__all__ = ('Tests',)
__version__ = '18.08.28'

from base import TestsBase

from pygeodesy import RangeError, toUtm, utm


class Tests(TestsBase):

    def testUtmTMcoord(self, coord, line):
        # <http://Geographiclib.SourceForge.io/html/transversemercator.html>
        # <http://Zenodo.org/record/32470#.W4LEJS2ZON8>
        # format: lat lon easting northing convergence scale
        lat, lon, e1, n1, c1, s1 = coord.split()
        try:
            _, e2, n2, _, c2, s2 = toUtm(lat, lon, cmoff=False)
            self.test(line + 'easting',  e2, float(e1), fmt='%.2f')
            self.test(line + 'northing', n2, float(n1), fmt='%.2f')
            # self.test(line + 'convergence', c2, float(c1), fmt='%.4f')
            # self.test(line + 'scale',       s2, float(s1), fmt='%.4f')
        except RangeError as x:
            self.test(line + 'RangeError', x, None, known=True)


if __name__ == '__main__':

    import sys

    t = Tests(__file__, __version__, utm, verbose=False)

    if len(sys.argv) > 1:
        coords = sys.argv[1]
        with open(coords) as f:
            n, testUtmTMcoord = 0, t.testUtmTMcoord
            for coord in f.readlines():
                n += 1
                if n > 250000:
                    break
                testUtmTMcoord(coord, 'line %d ' % (n,))

        t.results(nl=0)
    else:
        # XXX TMcoords.dat file is 30+ MB
        t.results(passed='SKIPPED', nl=0)

    t.exit()


# Typical results on macOS, verbose=True:

#   testing testUtmTMcoords.py 18.08.25 (module pygeodesy.utm 18.08.26)
#   test 1 line 1 easting: 1548706.79
#   test 2 line 1 northing: 8451449.20
#   test 3 line 2 easting: 2624150.74
#   test 4 line 2 northing: 1204434.04
#   test 5 line 3 easting: 9855841.23
#   test 6 line 3 northing: 6145496.12
#   test 7 line 4 easting: 3206390.69
#   test 8 line 4 northing: 2650745.40
#   test 9 line 5 easting: 4328154.08
#   ....
#   test 45 line 23 easting: 263004.77
#   test 46 line 23 northing: 4493669.76
#   test 47 line 24 easting: 3217221.74
#   test 48 line 24 northing: 437776.12
#   test 49 line 25 easting: 14661143.10  FAILED, expected 14661142.44
#   test 50 line 25 northing: 7476099.31  FAILED, expected 7476100.82
#   ....
#   test 110 line 55 northing: 8382823.86
#   test 111 line 56 easting: 5090358.13
#   test 112 line 56 northing: 4318294.13
#   test 113 line 57 RangeError: lat outside UTM: 84.9869301372  FAILED, KNOWN, expected None
#   test 114 line 58 easting: 3447670.55
#   test 115 line 58 northing: 3680238.89
#   ...
#   test 155 line 78 northing: 1762592.04
#   test 156 line 79 easting: 4880570.30
#   test 157 line 79 northing: 3325433.99
#   test 158 line 80 easting: -55111265728907.72  FAILED, expected 23930680.08
#   test 159 line 80 northing: 35534561035326.49  FAILED, expected 7491462.10
#   ....
#   test 205 line 103 northing: 6793804.30
#   test 206 line 104 easting: 8385525.15
#   test 207 line 104 northing: 6924932.34
#   test 208 line 105 easting: 2297510.79
#   test 209 line 105 northing: 2805666.96
#   test 210 line 106 easting: 13749544.95  FAILED, expected 13749544.92
#   test 211 line 106 northing: 8288728.50  FAILED, expected 8288728.38
#   ....
#   test 498619 line 249999 easting: 2093722.41
#   test 498620 line 249999 northing: 3819079.20
#   test 498621 line 250000 easting: 2539602.48
#   test 498622 line 250000 northing: 581070.98
#   19290 testUtmTMcoords.py tests (3.9%) FAILED, incl. 1378 KNOWN (PyGeodesy 18.8.26 Python 2.7.15 64bit geographiclib 1.49 numpy 1.14.0 macOS 10.13.6) 128.797 sec


#   testing testUtmTMcoords.py 18.08.25 (module pygeodesy.utm 18.08.26)
#   ....
#   test 498621 line 250000 easting: 2539602.48
#   test 498622 line 250000 northing: 581070.98
#   19290 testUtmTMcoords.py tests (3.9%) FAILED, incl. 1378 KNOWN (PyGeodesy 18.8.26 PyPy-Python 2.7.13 64bit macOS 10.13.6) 24.763 sec


#   testing testUtmTMcoords.py 18.08.25 (module pygeodesy.utm 18.08.26)
#   ....
#   test 498621 line 250000 easting: 2539602.48
#   test 498622 line 250000 northing: 581070.98
#   19290 testUtmTMcoords.py tests (3.9%) FAILED, incl. 1378 KNOWN (PyGeodesy 18.8.26 Python 3.7.0 64bit macOS 10.13.6) 81.966 sec


#   testing testUtmTMcoords.py 18.08.25 (module pygeodesy.utm 18.08.26)
#   ....
#   test 498621 line 250000 easting: 2539602.48
#   test 498622 line 250000 northing: 581070.98
#   19290 testUtmTMcoords.py tests (3.9%) FAILED, incl. 1378 KNOWN (PyGeodesy 18.8.26 Intel-Python 3.5.3 64bit numpy 1.11.3 macOS 10.13.6) 84.423 sec


# Typical results on iOS/Pythonista3, verbose=False:

#   testing testUtmTMcoords.py 18.08.28 (module utm 18.08.28)

#   test 498457 line 249918 easting: 15832598.46  FAILED, expected 15832602.28
#   test 498458 line 249918 northing: 8465626.52  FAILED, expected 8465566.94
#   19290 testUtmTMcoords.py tests (3.9%) FAILED, incl. 1378 KNOWN (PyGeodesy 18.8.28 Python 2.7.12 64bit numpy 1.8.0 iOS 11.4.1) 272.860 sec


#   testing testUtmTMcoords.py 18.08.28 (module utm 18.08.28)

#   test 498457 line 249918 easting: 15832598.46  FAILED, expected 15832602.28
#   test 498458 line 249918 northing: 8465626.52  FAILED, expected 8465566.94
#   19290 testUtmTMcoords.py tests (3.9%) FAILED, incl. 1378 KNOWN (PyGeodesy 18.8.28 Python 3.6.1 64bit numpy 1.8.0 iOS 11.4.1) 197.901 sec
