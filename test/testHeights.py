
# -*- coding: utf-8 -*-

# Test the height interpolators.

__all__ = ('Tests',)
__version__ = '19.02.24'

import warnings  # PYCHOK expected
# RuntimeWarning: numpy.ufunc size changed, may indicate binary
# incompatibility. Expected 192 from C header, got 216 from PyObject
# warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings('ignore')  # or 'error'

from base import scipy, TestsBase

from pygeodesy import fStr, HeightError, \
                      HeightCubic, HeightIDW, HeightLinear, \
                      HeightLSQBiSpline, HeightSmoothBiSpline
from pygeodesy.sphericalTrigonometry import LatLon


def _fStrs(floats):
    fmt = '[%s]'  if isinstance(floats, list)  else \
          '(%s,)' if isinstance(floats, tuple) else '%s'
    return fmt % (fStr(floats, prec=9),)


class Tests(TestsBase):

    def testHeight(self, H, kts, lli, expected, lats, lons):
        interpolator = H(kts)
        self.testHeightError(interpolator)
        h = interpolator(lli)
        self.test(H.__name__, h, expected, fmt='%.9f')
        self.test(H.__name__+'-float', type(h), float)
        h2 = interpolator.height(lli.lat, lli.lon)
        self.test(H.__name__+'-latlon', h2 == h, True)
        hs = interpolator(*kts)  # return tuple
        self.test(H.__name__+'-tuple', type(hs), tuple)
        self.test(H.__name__+'-tuple-float', type(hs[0]), float)
        self.test(H.__name__+'-tuple-float', type(hs[1]), float)
        hs = interpolator(kts)  # return list
        self.test(H.__name__+'-list', type(hs), list)
        self.test(H.__name__+'-list-float', type(hs[0]), float)
        self.test(H.__name__+'-list-float', type(hs[1]), float)
        hs2 = interpolator.height(lats, lons)
        self.test(H.__name__+'-latlon', hs2 == hs, True)

    def testHeightError(self, interpolator):
        try:  # force an error
            h = interpolator(0.0, 1.0)
        except HeightError as x:
            h = str(x)
        self.test('HeightError', h, "'float' object has no attribute 'lon': 0.0")

    def testHeightIDW(self, adjust, kts, lli, expected):
        interpolator = HeightIDW(kts, adjust=adjust)
        self.testHeightError(interpolator)
        kt, adjust = kts[2], str(adjust)
        expected1 = '%.1f' % (kt.height,)
        h = interpolator(lli)  # return scalar
        self.test('HeightIDW-'+adjust, h, expected, fmt='%.9f')
        self.test('HeightIDW-float', type(h), float)
        h2 = interpolator.height(lli.lat, lli.lon)
        self.test('HeightIDW-latlon', h2 == h, True)
        h = interpolator(kt)  # return scalar
        self.test('HeightIDW-'+adjust, h, expected1, fmt='%.1f')
        self.test('HeightIDW-float', type(h), float)
        h2 = interpolator.height(kt.lat, kt.lon)
        self.test('HeightIDW-latlon', h2 == h, True)
        hs = interpolator(lli, kt)  # return tuple
        self.test('HeightIDW-'+adjust, _fStrs(hs), '(%s, %s,)' % (expected, expected1))
        self.test('HeightIDW-tuple', type(hs), tuple)
        self.test('HeightIDW-tuple-float', type(hs[0]), float)
        self.test('HeightIDW-tuple-float', type(hs[1]), float)
        hs = interpolator([lli, kt])  # return list
        self.test('HeightIDW-'+adjust, _fStrs(hs), '[%s, %s]' % (expected, expected1))
        self.test('HeightIDW-list', type(hs), list)
        self.test('HeightIDW-list-float', type(hs[0]), float)
        self.test('HeightIDW-list-float', type(hs[1]), float)

    def testHeights(self):
        kts = LatLon(0.4, 0.9, 1), LatLon(1.5, 1.5, 3), \
              LatLon(1, 0.5, 5), LatLon(0.5, 1.4, 7), LatLon(1.2, 1, 7)
        lli = LatLon(1, 1)
        self.testHeightIDW(None,  kts, lli, '6.108538037')
        self.testHeightIDW(True,  kts, lli, '6.167071011')
        self.testHeightIDW(False, kts, lli, '6.166920194')

        kts = LatLon(1.1, 1, 2), LatLon(2.1, 2, 2), \
              LatLon(1.2, 4, 3), LatLon(2.2, 3, 3)
        lli = kts[0].intersection(*kts[1:])
        self.test('intersection', lli, '02.64932°N, 002.550079°E, +2.50m')  # mean height
        self.testHeightIDW(None,  kts, lli, '2.592742938')
        self.testHeightIDW(True,  kts, lli, '2.592742049')
        self.testHeightIDW(False, kts, lli, '2.592735027')

        if scipy:
            interpolator = HeightLinear(kts)
            h = interpolator(lli)
            self.test('HeightLinear', h, '2.536626441', fmt='%.9f')
            self.test('HeightLinear-float', type(h), float)
            h2 = interpolator.height(lli.lat, lli.lon)
            self.test('HeightLinear-latlon', h2 == h, True)

            kts = tuple(kts)  # XXX grid-, mesh-like
            kts += LatLon(1.1, 2, 2), LatLon(2.1, 3, 2), LatLon(1.2, 5, 3), LatLon(2.2, 4, 3)
            kts += LatLon(1.1, 3, 2), LatLon(2.1, 4, 2), LatLon(1.2, 6, 3), LatLon(2.2, 5, 3)
            kts += LatLon(1.1, 4, 2), LatLon(2.1, 5, 2), LatLon(1.2, 7, 3), LatLon(2.2, 6, 3)
            # list of lats and lons
            lats, lons = zip(*[(ll.lat, ll.lon) for ll in kts])

            self.testHeight(HeightCubic,          kts, lli, '3.000000000', lats, lons)
            self.testHeight(HeightIDW,            kts, lli, '2.407864110', lats, lons)
            self.testHeight(HeightLinear,         kts, lli, '3.000000000', lats, lons)
            self.testHeight(HeightLSQBiSpline,    kts, lli, '6.419251669', lats, lons)
            self.testHeight(HeightSmoothBiSpline, kts, lli, '2.598922541', lats, lons)

        else:
            self.skip('no scipy', n=53)


if __name__ == '__main__':  # PYCHOK internal error?

    t = Tests(__file__, __version__)
    t.testHeights()
    t.results()
    t.exit()

    _ = '''
    import numpy as np
    from scipy.interpolate import SmoothSphereBivariateSpline

    theta = np.linspace(0., np.pi, 7)
    phi = np.linspace(0., 2*np.pi, 9)
    data = np.zeros((theta.shape[0], phi.shape[0]))
    data[1:-1,1], data[1:-1,-1] = 1., 1.
    data[1,1:-1], data[-2,1:-1] = 1., 1.
    data[2:-2,2], data[2:-2,-2] = 2., 2.
    data[2,2:-2], data[-3,2:-2] = 2., 2.
    data[3,3:-2] = 3.
    data = np.roll(data, 4, 1)
    print('_data', data)

    lats, lons = np.meshgrid(theta, phi)
    print('lats', lats)
    print('lons', lons)
    lut = SmoothSphereBivariateSpline(lats.ravel(), lons.ravel(), data.T.ravel(), s=3.5)

    data_orig = lut(theta, phi)
    print('_orig', data_orig)

    fine_lats = np.linspace(0., np.pi, 15)
    fine_lons = np.linspace(0., 2*np.pi, 19)
    data_smth = lut(fine_lats, fine_lons)
    print('_smth', data_smth)
'''
    del _
