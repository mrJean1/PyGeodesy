
# -*- coding: utf-8 -*-

# Test L{heights} interpolators.

__all__ = ('Tests',)
__version__ = '23.05.14'

import warnings  # PYCHOK expected
# RuntimeWarning: numpy.ufunc size changed, may indicate binary
# incompatibility. Expected 192 from C header, got 216 from PyObject
# warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings('ignore')  # or 'error'

from bases import coverage, geographiclib, startswith, TestsBase

from pygeodesy import Datums, fstr, HeightError, SciPyError, \
                      HeightCubic, HeightIDWcosineAndoyerLambert, \
                      HeightIDWcosineForsytheAndoyerLambert, HeightIDWcosineLaw, \
                      HeightIDWdistanceTo, HeightIDWequirectangular, \
                      HeightIDWeuclidean, HeightIDWexact, HeightIDWflatLocal, \
                      HeightIDWflatPolar, HeightIDWhaversine, HeightIDWhubeny, \
                      HeightIDWkarney, HeightIDWthomas, HeightIDWvincentys, \
                      HeightLinear, HeightLSQBiSpline, HeightSmoothBiSpline
#                     HeightIDW, HeightIDW2, HeightIDW3
from pygeodesy.interns import _DOT_
from pygeodesy.sphericalTrigonometry import LatLon


def _fstr(floats):
    fmt = '[%s]'  if isinstance(floats, list)  else \
          '(%s,)' if isinstance(floats, tuple) else '%s'
    return fmt % (fstr(floats, prec=9),)


def _kwdstr(kwds):
    return ', '.join('%s=%r' % t for t in sorted(kwds.items()))


class Tests(TestsBase):

    def testHeight(self, H, kts, lli, expected, lats, lons, **kwds):
        interpolator = H(kts, **kwds)
        self.testCopy(interpolator)
        self.testHeightError(interpolator)
        h = interpolator(lli)
        self.test(H.__name__, h, expected, fmt='%.9f')
        self.test(H.__name__+'(float)', type(h), float)
        h2 = interpolator.height(lli.lat, lli.lon)
        self.test(H.__name__+'(latlon)', h2 == h, True)
        hs = interpolator(*kts)  # return tuple
        self.test(H.__name__+'(tuple)', type(hs), tuple)
        self.test(H.__name__+'(tuple-float)', type(hs[0]), float)
        self.test(H.__name__+'(tuple-float)', type(hs[1]), float)
        hs = interpolator(kts)  # return list
        self.test(H.__name__+'(list)', type(hs), list)
        self.test(H.__name__+'(list-float)', type(hs[0]), float)
        self.test(H.__name__+'(list-float)', type(hs[1]), float)
        hs2 = interpolator.height(lats, lons)
        self.test(H.__name__+'(latlon)', hs2 == hs, True, nt=1)

    def testHeightError(self, interpolator, *attrs):
        try:  # force an error
            h = interpolator(9.0, 18.0)
        except HeightError as x:
            h = str(x)
        self.test('HeightError', h, "llis[0] (9.0): 'float' object has no attribute ", known=startswith)
        if coverage:
            n = interpolator.__class__.__name__
            for a in ('datum', 'kmin', 'wrap') + attrs:  # XXX 'adjust'
                t = getattr(interpolator, a)
                self.test(_DOT_(n, a), t, t)

    def testIDW(self, IDW, kts, lli, expected, **kwds):
        interpolator = IDW(kts, **kwds)
        self.testCopy(interpolator)
        kt = kts[2]
        expected1 = '%.1f' % (kt.height,)
        kwdstr = '(%s)' % (_kwdstr(kwds),)
        h = interpolator(lli)  # return scalar
        k = False  # f = float(expected)
        self.test(IDW.__name__+kwdstr, h, expected, fmt='%.9f', known=k)
        self.test(IDW.__name__+'(float)', type(h), float)
        h2 = interpolator.height(lli.lat, lli.lon)
        self.test(IDW.__name__+'(latlon)', h2 == h, True)
        h = interpolator(kt)  # return scalar
        t = '%s, %s' % (expected, expected1)
        self.test(IDW.__name__+kwdstr, h, expected1, fmt='%.1f', known=k)
        self.test(IDW.__name__+'(float)', type(h), float)
        h2 = interpolator.height(kt.lat, kt.lon)
        self.test(IDW.__name__+'(latlon)', h2 == h, True)
        hs = interpolator(lli, kt)  # return tuple
        self.test(IDW.__name__+kwdstr, _fstr(hs), '(%s,)' % (t,), known=k)
        self.test(IDW.__name__+'(tuple)', type(hs), tuple)
        self.test(IDW.__name__+'(tuple-float)', type(hs[0]), float)
        self.test(IDW.__name__+'(tuple-float)', type(hs[1]), float)
        hs = interpolator([lli, kt])  # return list
        self.test(IDW.__name__+kwdstr, _fstr(hs), '[%s]' % (t,), known=k)
        self.test(IDW.__name__+'(list', type(hs), list)
        self.test(IDW.__name__+'(list-float)', type(hs[0]), float)
        self.test(IDW.__name__+'(list-float)', type(hs[1]), float)
        self.testHeightError(interpolator, 'beta')

        if coverage:
            for a in ('adjust', 'beta', 'kmin', 'wrap'):
                t = getattr(interpolator, a)
                self.test(_DOT_(IDW.__name__, a), t, t)

    def testHeights(self):
        kts = LatLon(0.4, 0.9, 1), LatLon(1.5, 1.5, 3), \
              LatLon(1, 0.5, 5), LatLon(0.5, 1.4, 7), LatLon(1.2, 1, 7)
        lli = LatLon(1, 1)
#       self.testIDW(HeightIDW,                kts, lli, '6.142945781', adjust=True)
        self.testIDW(HeightIDWcosineAndoyerLambert,         kts, lli, '6.111157919', wrap=False)  # XXX 6.108538037
        self.testIDW(HeightIDWcosineForsytheAndoyerLambert, kts, lli, '6.111158354', wrap=False)  # XXX 6.108538037
        self.testIDW(HeightIDWcosineLaw,       kts, lli, '6.108538037', wrap=True)
        self.testIDW(HeightIDWcosineLaw,       kts, lli, '6.108538037', wrap=False)
        self.testIDW(HeightIDWdistanceTo,      kts, lli, '6.108538037')
#       self.testIDW(HeightIDW2,               kts, lli, '6.108538529', adjust=True, wrap=False)
        self.testIDW(HeightIDWequirectangular, kts, lli, '6.860459007', adjust=True, wrap=True)  # XXX 6.108538529
#       self.testIDW(HeightIDW2,               kts, lli, '6.108614369', adjust=False, wrap=False)
        self.testIDW(HeightIDWequirectangular, kts, lli, '6.860498092', adjust=False, wrap=True)  # XXX 6.108614369
        self.testIDW(HeightIDWeuclidean,       kts, lli, '6.142945781')
        self.testIDW(HeightIDWexact,           kts, lli, '6.108538037')
        self.testIDW(HeightIDWflatLocal,       kts, lli, '6.111158784')  # XXX 6.860459007
        self.testIDW(HeightIDWflatPolar,       kts, lli, '6.261469975')
#       self.testIDW(HeightIDW3,               kts, lli, '6.108538037', wrap=True)
        self.testIDW(HeightIDWhaversine,       kts, lli, '6.108538037')
        self.testIDW(HeightIDWhubeny,          kts, lli, '6.111158784')  # XXX 6.860459007
        if geographiclib:
            self.testIDW(HeightIDWkarney,      kts, lli, '6.111158743', wrap=True,  datum=Datums.WGS84)
            self.testIDW(HeightIDWkarney,      kts, lli, '6.111158743', wrap=False, datum=Datums.WGS84)
            self.testIDW(HeightIDWkarney,      kts, lli, '6.108538037', wrap=True,  datum=Datums.Sphere)
            self.testIDW(HeightIDWkarney,      kts, lli, '6.108538037', wrap=False, datum=Datums.Sphere)
        self.testIDW(HeightIDWthomas,          kts, lli, '6.111159454', wrap=True)   # XXX 6.108538037
        self.testIDW(HeightIDWthomas,          kts, lli, '6.111159454', wrap=False)  # XXX 6.108538037
        self.testIDW(HeightIDWvincentys,       kts, lli, '6.108538037', wrap=True)
        self.testIDW(HeightIDWvincentys,       kts, lli, '6.108538037', wrap=False)

        kts = LatLon(1.1, 1, 2), LatLon(2.1, 2, 2), \
              LatLon(1.2, 4, 3), LatLon(2.2, 3, 3)
        lli = kts[0].intersection(*kts[1:])
        self.test('intersection', lli, '02.64932°N, 002.550079°E, +2.50m', nt=1)  # mean height

#       self.testIDW(HeightIDW,                kts, lli, '2.592748835', adjust=True)
        self.testIDW(HeightIDWcosineAndoyerLambert,         kts, lli, '2.592742781', wrap=False)  # XXX 2.592742938
        self.testIDW(HeightIDWcosineForsytheAndoyerLambert, kts, lli, '2.592742781', wrap=False)  # XXX 2.592742938
        self.testIDW(HeightIDWcosineLaw,       kts, lli, '2.592742938', wrap=True)
        self.testIDW(HeightIDWcosineLaw,       kts, lli, '2.592742938', wrap=False)
#       self.testIDW(HeightIDW2,               kts, lli, '2.592743455', adjust=True, wrap=False)
        self.testIDW(HeightIDWequirectangular, kts, lli, '2.689429914', adjust=True, wrap=True)  # XXX 2.592743455
#       self.testIDW(HeightIDW2,               kts, lli, '2.592732915', adjust=False, wrap=False)
        self.testIDW(HeightIDWequirectangular, kts, lli, '2.689413401', adjust=False, wrap=True)  # XXX 2.592732915
        self.testIDW(HeightIDWeuclidean,       kts, lli, '2.592735541', adjust=False)
        self.testIDW(HeightIDWexact,           kts, lli, '2.592742938')
        self.testIDW(HeightIDWflatLocal,       kts, lli, '2.592743288')  # XXX 2.689429914
        self.testIDW(HeightIDWflatPolar,       kts, lli, '2.592973059')
#       self.testIDW(HeightIDW3,               kts, lli, '2.592742938', wrap=True)
        self.testIDW(HeightIDWhaversine,       kts, lli, '2.592742938')
        self.testIDW(HeightIDWhubeny,          kts, lli, '2.592743288')  # XXX 2.689429914
        if geographiclib:
            self.testIDW(HeightIDWkarney,      kts, lli, '2.592742915', wrap=True,  datum=Datums.WGS84)
            self.testIDW(HeightIDWkarney,      kts, lli, '2.592742915', wrap=False, datum=Datums.WGS84)
            self.testIDW(HeightIDWkarney,      kts, lli, '2.592742938', wrap=True,  datum=Datums.Sphere)
            self.testIDW(HeightIDWkarney,      kts, lli, '2.592742938', wrap=False, datum=Datums.Sphere)
        self.testIDW(HeightIDWthomas,          kts, lli, '2.592742781', wrap=True)   # XXX 2.592742938
        self.testIDW(HeightIDWthomas,          kts, lli, '2.592742781', wrap=False)  # XXX 2.592742938
        self.testIDW(HeightIDWvincentys,       kts, lli, '2.592742938', wrap=True)
        self.testIDW(HeightIDWvincentys,       kts, lli, '2.592742938', wrap=False)

        try:
            interpolator = HeightLinear(kts)
            self.testCopy(interpolator)
            h = interpolator(lli)
            self.test(HeightLinear.__name__, h, '2.536626441', fmt='%.9f')
            self.test(HeightLinear.__name__+'(float)', type(h), float)
            h2 = interpolator.height(lli.lat, lli.lon)
            self.test(HeightLinear.__name__+'(latlon)', h2 == h, True)

            kts = tuple(kts)  # XXX grid-, mesh-like
            kts += LatLon(1.1, 2, 2), LatLon(2.1, 3, 2), LatLon(1.2, 5, 3), LatLon(2.2, 4, 3)
            kts += LatLon(1.1, 3, 2), LatLon(2.1, 4, 2), LatLon(1.2, 6, 3), LatLon(2.2, 5, 3)
            kts += LatLon(1.1, 4, 2), LatLon(2.1, 5, 2), LatLon(1.2, 7, 3), LatLon(2.2, 6, 3)
            # unzip into a list of lats and of lons
            lats, lons = zip(*[(ll.lat, ll.lon) for ll in kts])

            self.testHeight(HeightCubic,              kts, lli, '3.000000000', lats, lons)
            self.testHeight(HeightIDWcosineAndoyerLambert,         kts, lli, '2.401661864', lats, lons)  # 2.402157442
            self.testHeight(HeightIDWcosineForsytheAndoyerLambert, kts, lli, '2.401661859', lats, lons)  # 2.402157442
            self.testHeight(HeightIDWcosineLaw,       kts, lli, '2.402157442', lats, lons)
            self.testHeight(HeightIDWequirectangular, kts, lli, '2.469718302', lats, lons)  # 2.402157181
            self.testHeight(HeightIDWeuclidean,       kts, lli, '2.409288552', lats, lons)
            self.testHeight(HeightIDWexact,           kts, lli, '2.402157442', lats, lons)
            self.testHeight(HeightIDWflatLocal,       kts, lli, '2.401661562', lats, lons)  # 2.469718302
            self.testHeight(HeightIDWflatPolar,       kts, lli, '2.370266641', lats, lons)
            self.testHeight(HeightIDWhaversine,       kts, lli, '2.402157442', lats, lons)
            self.testHeight(HeightIDWhubeny,          kts, lli, '2.401661562', lats, lons)  # 2.469718302
            if geographiclib:
                self.testHeight(HeightIDWkarney,      kts, lli, '2.402157442', lats, lons)  # datum=Datums.Sphere
            self.testHeight(HeightIDWthomas,          kts, lli, '2.401661653', lats, lons)  # 2.402157442
            self.testHeight(HeightIDWvincentys,       kts, lli, '2.402157442', lats, lons)
            self.testHeight(HeightLinear,             kts, lli, '3.000000000', lats, lons)
            self.testHeight(HeightLSQBiSpline,        kts, lli, '6.419251669', lats, lons)
            try:  # SciPy 1.9.0 issue
                self.testHeight(HeightLSQBiSpline,    kts, lli, '6.419251669', lats, lons, weight=2)  # coverage
            except SciPyError as x:
                self.test('SciPy 1.9 issue', str(x), str(x))
            try:  # SciPy 1.9.0 issue
                self.testHeight(HeightLSQBiSpline,    kts, lli, '6.419251669', lats, lons, weight=[1] * len(kts))
            except SciPyError as x:
                self.test('SciPy 1.19 issue', str(x), str(x))
            self.testHeight(HeightSmoothBiSpline,     kts, lli, '2.598922541', lats, lons)

        except ImportError as x:
            self.skip(str(x), n=80)


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
