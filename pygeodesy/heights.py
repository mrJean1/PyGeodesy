
# -*- coding: utf-8 -*-

u'''Classes to interpolate the height of C{LatLon} locations or
separate lat-/longitudes.

Except for the L{HeightIDW} interpolator, the height interpolators
in this module require U{numpy<http://PyPI.org/project/numpy>} and
U{scipy<http://SciPy.org>} to be installed.

Typical usage is as follows.  First create an interpolator from a
given set of C{LatLon} points with known heights, called knots.

C{hinterpolator = HeightXyz(knots, **options)}

Get the interpolated height of other C{LatLon} location(s) with

C{h = hinterpolator(ll)}

or

C{h0, h1, h2, ... = hinterpolator(ll0, ll1, ll2, ...)}

or

C{hs = hinterpolator(lls)}  # C{list, tuple, generator, ...}

For separate lat-/longitudes invoke the C{.height} method

C{h = hinterpolator.height(lat, lon)}

or

C{h0, h1, h2, ... = hinterpolator.height(lats, lons)}  # C{list, ...}


The knots do not need to be ordered for any of the height interpolators.

Errors from C{scipy} as raised as L{SciPyError}s.  Warnings issued by
C{scipy} can be thrown as L{SciPyWarning} exceptions, provided Python
C{warnings} are filtered accordingly, see L{SciPyWarning}.

@see: U{SciPy<http://docs.SciPy.org/doc/scipy/reference/interpolate.html>}.
'''

from fmath import EPS, fdot, fsum, isscalar, len2, map1
from formy import euclidean_, haversine_, _scaler
from lazily import _ALL_LAZY
from points import LatLon_
from utily import PI, PI2, PI_2, radiansPI, radiansPI2, unrollPI

__all__ = _ALL_LAZY.heights
__version__ = '19.04.03'


class HeightError(ValueError):  # imported by .geoids
    '''Height interpolator or interpolation error.
    '''
    pass


class SciPyError(HeightError):
    '''Error raised for C{SciPy} errors.
    '''
    pass


class SciPyWarning(HeightError):
    '''Exception thrown for C{SciPy} warnings.

       To raise C{SciPy} warnings as L{SciPyWarning} exceptions, Python
       C{warnings} must be filtered as U{warnings.filterwarnings('error')
       <http://docs.Python.org/3/library/warnings.html#the-warnings-filter>}
       I{prior to} C{import scipy} or by setting environment variable
       U{PYTHONWARNINGS<http://docs.Python.org/3/using/cmdline.html
       #envvar-PYTHONWARNINGS>} or with C{python} command line option
       U{-W<http://docs.Python.org/3/using/cmdline.html#cmdoption-w>}
       as C{error}.
    '''
    pass


def _alist(ais):
    # return tuple of floats, not numpy.float64s
    return list(map(float, ais))


def _allis2(llis, m=1, Error=HeightError):  # imported by .geoids
    # dtermine return type and convert lli C{LatLon}s to list
    if not isinstance(llis, tuple):  # llis are *args
        raise AssertionError('type(%s): %r' % ('*llis', llis))

    n = len(llis)
    if n == 1:  # convert single lli to 1-item list
        llis = llis[0]
        try:
            n, llis = len2(llis)
            _as = _alist  # return list of interpolated heights
        except TypeError:  # single lli
            n, llis = 1, [llis]
            _as = _ascalar  # return single interpolated heights
    else:  # of 0, 2 or more llis
        _as = _atuple  # return tuple of interpolated heights

    if n < m:
        raise Error('insufficient %s: %s, need %s' % ('llis', n, m))
    return _as, llis


def _ascalar(ais):  # imported by .geoids
    # return single float, not numpy.float64
    ais = list(ais)  # np.array, etc. to list
    if len(ais) != 1:
        raise AssertionError('len[%r] = %s == 1' % (ais, len(ais)))
    return float(ais[0])  # remove np.<type>


def _atuple(ais):
    # return tuple of floats, not numpy.float64s
    return tuple(map(float, ais))


def _axyllis4(atype, llis, m=1, off=True):
    # convert lli C{LatLon}s to tuples or C{NumPy} arrays of
    # C{SciPy} sphericals and determine the return type
    _as, llis = _allis2(llis, m=m)
    xis, yis, _ =  zip(*_xyhs(llis, off=off))  # unzip
    return _as, atype(xis), atype(yis), llis


def _ordedup(ts, lo=EPS, hi=PI2-EPS):
    # clip, order and remove duplicates
    p, ks = 0, []
    for k in sorted(max(lo, min(hi, t)) for t in ts):
        if k > p:
            ks.append(k)
            p = k
    return ks


def _SciPyIssue(x, *extras):  # imported by .geoids
    t = ' '.join(str(x).strip().split() + map(str, extras))
    if isinstance(x, (RuntimeWarning, UserWarning)):
        return SciPyWarning(t)
    else:
        return SciPyError(t)  # PYCHOK not really


def _xyhs(lls, off=True):
    # map (lat, lon, h) to (x, y, h) in radians, offset as
    # x: 0 <= lon <= PI2, y: 0 <= lat <= PI if off is True
    # else x: -PI <= lon <= PI, y: -PI_2 <= lat <= PI_2
    if off:
        xf = yf = 0.0
    else:  # undo offset
        xf, yf = PI, PI_2
    try:
        for ll in lls:
            yield (max(0.0, radiansPI2(ll.lon + 180.0)) - xf), \
                  (max(0.0, radiansPI( ll.lat +  90.0)) - yf), ll.height
    except AttributeError as x:
        raise HeightError('%s: %r' % (x, ll))


def _xyhs3(atype, m, knots, off=True):
    # convert knot C{LatLon}s to tuples or C{NumPy} arrays and C{SciPy} sphericals
    xs, ys, hs = zip(*_xyhs(knots, off=off))  # PYCHOK unzip
    n = len(hs)
    if n < m:
        raise HeightError('insufficient %s: %s, need %s' % ('knots', n, m))
    return map1(atype, xs, ys, hs)


class _HeightBase(object):  # imported by .geoids
    '''Interpolator base class.
    '''
    _kmin = 2     # min number of knots
    _np   = None  # numpy
    _np_v = None  # version
    _spi  = None  # scipy.interpolate
    _sp_v = None  # version

    def __call__(self, *unused):
        raise AssertionError('%s.%s not overloaded' % (self.__class__.__name__, '_()'))

    def _axyllis4(self, llis):
        return _axyllis4(self._np.array, llis)

    def _ev(self, *unused):
        raise AssertionError('%s.%s not overloaded' % (self.__class__.__name__, '_ev'))

    def _eval(self, llis):  # XXX single arg, not *args
        _as, xis, yis, _ = self._axyllis4(llis)
        try:  # SciPy .ev signature: y first, then x!
            return _as(self._ev(yis, xis))
        except Exception as x:
            raise _SciPyIssue(x)

    def _height(self, lats, lons, Error=HeightError):
        if isscalar(lats) and isscalar(lons):
            llis = LatLon_(lats, lons)
        else:
            n, lats = len2(lats)
            m, lons = len2(lons)
            if n != m:
                raise Error('non-matching %s: %s vs %s' % ('len', n, m))
            llis = [LatLon_(*ll) for ll in zip(lats, lons)]
        return self(llis)  # __call__(lli) or __call__(llis)

    def _NumSciPy(self, throwarnings=False):
        # import numpy and scipy
        if throwarnings:  # raise SciPyWarnings, but ...
            # ... not if scipy has been imported already
            import sys
            if 'scipy' not in sys.modules:
                import warnings
                warnings.filterwarnings('error')

        import scipy as sp
        import scipy.interpolate as spi
        import numpy as np

        _HeightBase._np   = np
        _HeightBase._np_v = np.__version__
        _HeightBase._spi  = spi
        _HeightBase._sp_v = sp.__version__

        return np, spi

    def _xyhs3(self, knots):
        return _xyhs3(self._np.array, self._kmin, knots)


class HeightCubic(_HeightBase):
    '''Height interpolator based on C{SciPy} U{interp2d<http://docs.SciPy.org/
       doc/scipy/reference/generated/scipy.interpolate.interp2d.html>}
       C{kind='cubic'}.
    '''
    _interp2d = None
    _kind     = 'cubic'
    _kmin     = 16

    def __init__(self, knots):
        '''New L{HeightCubic} interpolator.

           @param knots: The points with known height (C{LatLon}s).

           @raise HeightError: Insufficient number of I{knots} or
                               invalid I{knot}.

           @raise ImportError: Package C{numpy} or C{scipy} not found
                               or not installed.

           @raise SciPyError: A C{scipy.interpolate.interp2d} issue.

           @raise SciPyWarning: A C{scipy.interpolate.interp2d} warning
                                as exception.
        '''
        _, spi = self._NumSciPy()

        xs, ys, hs = self._xyhs3(knots)
        try:  # SciPy.interpolate.interp2d kind 'linear' or 'cubic'
            self._interp2d = spi.interp2d(xs, ys, hs, kind=self._kind)
        except Exception as x:
            raise _SciPyIssue(x)

    def __call__(self, *llis):
        '''Interpolate the height for one or several locations.

           @param llis: The location or locations (C{LatLon}, ... or
                        C{LatLon}s).

           @return: A single interpolated height (C{float}) or a list
                    or tuple of interpolated heights (C{float}s).

           @raise HeightError: Insufficient number of I{llis} or
                               invalid I{lli}.

           @raise SciPyError: A C{scipy.interpolate.interp2d} issue.

           @raise SciPyWarning: A C{scipy.interpolate.interp2d} warning
                                as exception.
        '''
        return _HeightBase._eval(self, llis)

    def _ev(self, yis, xis):  # PYCHOK expected
        # to make SciPy .interp2d signature(x, y), single (x, y)
        # match SciPy .ev signature(ys, xs), flipped multiples
        return map(self._interp2d, xis, yis)

    def height(self, lats, lons):
        '''Interpolate the height for one or several lat-/longitudes.

           @param lats: Latitude or latitudes (C{degrees} or C{degrees}s).
           @param lons: Longitude or longitudes (C{degrees} or C{degrees}s).

           @return: A single interpolated height (C{float}) or a list of
                    interpolated heights (C{float}s).

           @raise HeightError: Insufficient or non-matching number of
                               I{lats} and I{lons}.

           @raise SciPyError: A C{scipy.interpolate.interp2d} issue.

           @raise SciPyWarning: A C{scipy.interpolate.interp2d} warning
                                as exception.
        '''
        return _HeightBase._height(self, lats, lons)


class HeightLinear(HeightCubic):
    '''Height interpolator based on C{SciPy} U{interp2d<http://docs.SciPy.org/
       doc/scipy/reference/generated/scipy.interpolate.interp2d.html>}
       C{kind='linear]}.
    '''
    _kind = 'linear'
    _kmin = 2

    def __init__(self, knots):
        '''New L{HeightLinear} interpolator.

           @param knots: The points with known height (C{LatLon}s).

           @raise HeightError: Insufficient number of I{knots} or
                               invalid I{knot}.

           @raise ImportError: Package C{numpy} or C{scipy} not found
                               or not installed.

           @raise SciPyError: A C{scipy.interpolate.interp2d} issue.

           @raise SciPyWarning: A C{scipy.interpolate.interp2d} warning
                                as exception.
        '''
        HeightCubic.__init__(self, knots)

    __call__ = HeightCubic.__call__
    height   = HeightCubic.height


class _HeightIDW(_HeightBase):
    '''Base class for U{Inverse Distance Weighting
       <http://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW)
       height interpolators.
    '''
    _b  = 0  # negative distance power
    _hs = ()  # known heights
    _xs = ()  # knot lons
    _ys = ()  # knot lats

    def __init__(self, knots, beta=2):
        '''New L{_HeightIDW} interpolator.
        '''
        self._xs, self._ys, self._hs = _xyhs3(tuple, self._kmin, knots, off=False)

        self._b = -int(beta)
        if not (0 < beta < 4 and beta == -self._b):
            raise HeightError('invalid %s=%r' % ('beta', beta))

    def __call__(self, *llis):
        '''Interpolate the height for one or several locations.

           @param llis: The location or locations (C{LatLon}, ... or
                        C{LatLon}s).

           @return: A single interpolated height (C{float}) or a list
                    or tuple of interpolated heights (C{float}s).

           @raise HeightError: Insufficient number of I{llis} or
                               invalid I{lli}.
        '''
        _as, xis, yis, _ = _axyllis4(tuple, llis, off=False)
        return _as(map(self._hIDW, xis, yis))

    def _distances(self, x, y):  # PYCHOK unused (x, y) radians
        '''Must be overloaded.
        '''
        raise NotImplementedError('method: %s' % ('_distances',))

    def _hIDW(self, x, y):
        # interpolate height at (x, y) radians
        ws = tuple(self._distances(x, y))
        w, h = min(zip(ws, self._hs))
        if w > EPS:
            ws = tuple(w**self._b for w in ws)
            h = fdot(ws, *self._hs) / fsum(ws)
        return h

    def height(self, lats, lons):
        '''Interpolate the height for one or several lat-/longitudes.

           @param lats: Latitude or latitudes (C{degrees} or C{degrees}s).
           @param lons: Longitude or longitudes (C{degrees} or C{degrees}s).

           @return: A single interpolated height (C{float}) or a list of
                    interpolated heights (C{float}s).

           @raise HeightError: Insufficient or non-matching number of
                               I{lats} and I{lons}.
        '''
        return _HeightBase._height(self, lats, lons)


class HeightIDW(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <http://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW) and
       the I{angular} C{Euclidean} distance from function L{euclidean_}.

       @see: U{IDW<http://www.Geo.FU-Berlin.De/en/v/soga/Geodata-analysis/
             geostatistics/Inverse-Distance-Weighting/index.html>},
             U{SHEPARD_INTERP_2D<http://People.SC.FSU.edu/~jburkardt/c_src/
             shepard_interp_2d/shepard_interp_2d.html>} and function L{euclidean_}.
    '''
    _adjust = True

    def __init__(self, knots, adjust=True, beta=2):
        '''New L{HeightIDW} interpolator.

           @param knots: The points with known height (C{LatLon}s).
           @keyword adjust: Adjust the longitudinal delta by the cosine
                            of the mean latitude for I{adjust}=C{True}.
           @keyword beta: Inverse distance power (C{int} 1, 2, or 3).

           @raise HeightError: Insufficient number of I{knots} or invalid
                               I{knot}, I{adjust} or I{beta}.
        '''
        if adjust not in (True, False):
            raise HeightError('invalid %s=%r' % ('adjust', adjust))
        self._adjust = adjust
        _HeightIDW.__init__(self, knots, beta=beta)

    def _distances(self, x, y):  # (x, y) radians
        for xk, yk in zip(self._xs, self._ys):
            yield euclidean_(y, yk, x - xk, adjust=self._adjust)

    __call__ = _HeightIDW.__call__
    height   = _HeightIDW.height


class HeightIDW2(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <http://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW)
       and the C{equirectangular} distance (in radians squared) like
       function L{equirectangular_}.

       @see: U{IDW<http://www.Geo.FU-Berlin.De/en/v/soga/Geodata-analysis/
             geostatistics/Inverse-Distance-Weighting/index.html>},
             U{SHEPARD_INTERP_2D<http://People.SC.FSU.edu/~jburkardt/c_src/
             shepard_interp_2d/shepard_interp_2d.html>} and function L{euclidean_}.
    '''
    _adjust = True
    _wrap   = False

    def __init__(self, knots, adjust=True, wrap=False):
        '''New L{HeightIDW2} interpolator.

           @param knots: The points with known height (C{LatLon}s).
           @keyword adjust: Adjust the wrapped, unrolled longitudinal
                            delta by the cosine of the mean latitude (C{bool}).
           @keyword wrap: Wrap and L{unroll180} longitudes (C{bool}).

           @raise HeightError: Insufficient number of I{knots} or invalid
                               I{knot}.
        '''
        if not adjust:
            self._adjust = False
        if wrap:
            self._wrap = True
        _HeightIDW.__init__(self, knots, beta=1)  # distance**2

    def _distances(self, x, y):  # (x, y) radians
        for xk, yk in zip(self._xs, self._ys):
            dx, _ = unrollPI(xk, x, wrap=self._wrap)
            if self._adjust:
                dx *= _scaler(yk, y)
            yield dx**2 + (y - yk)**2  # like equirectangular_ distance2

    __call__ = _HeightIDW.__call__
    height   = _HeightIDW.height


class HeightIDW3(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <http://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW) and
       the I{angular} C{Haversine} distance from function L{haversine_}.

       @see: U{IDW<http://www.Geo.FU-Berlin.De/en/v/soga/Geodata-analysis/
             geostatistics/Inverse-Distance-Weighting/index.html>},
             U{SHEPARD_INTERP_2D<http://People.SC.FSU.edu/~jburkardt/c_src/
             shepard_interp_2d/shepard_interp_2d.html>} and function L{euclidean_}.
    '''
    _wrap = False

    def __init__(self, knots, beta=2, wrap=False):
        '''New L{HeightIDW3} interpolator.

           @param knots: The points with known height (C{LatLon}s).
           @keyword beta: Inverse distance power (C{int} 1, 2, or 3).
           @keyword wrap: Wrap and L{unroll180} longitudes (C{bool}).

           @raise HeightError: Insufficient number of I{knots} or invalid
                               I{knot} or I{beta}.
        '''
        if wrap:
            self._warp = True
        _HeightIDW.__init__(self, knots, beta=beta)

    def _distances(self, x, y):  # (x, y) radians
        for xk, yk in zip(self._xs, self._ys):
            d, _ = unrollPI(xk, x, wrap=self._wrap)
            yield haversine_(y, yk, d)

    __call__ = _HeightIDW.__call__
    height   = _HeightIDW.height


class HeightLSQBiSpline(_HeightBase):
    '''Height interpolator using C{SciPy} U{LSQSphereBivariateSpline
       <http://docs.SciPy.org/doc/scipy/reference/generated/scipy.
       interpolate.LSQSphereBivariateSpline.html>}.
    '''
    _kmin = 16  # k = 3, always

    def __init__(self, knots, weight=None):
        '''New L{HeightLSQBiSpline} interpolator.

           @param knots: The points with known height (C{LatLon}s).
           @keyword weight: Optional weight or weights for each I{knot}
                            (C{scalar} or C{scalar}s).

           @raise HeightError: Insufficient number of I{knots} or
                               I{weight}s or invalid I{knot} or I{weight}.

           @raise ImportError: Package C{numpy} or C{scipy} not found
                               or not installed.

           @raise SciPyError: A C{LSQSphereBivariateSpline} issue.

           @raise SciPyWarning: A C{LSQSphereBivariateSpline} warning
                                as exception..
        '''
        np, spi = self._NumSciPy()

        xs, ys, hs = self._xyhs3(knots)
        m = len(hs)

        if not weight:
            w = None  # default
        elif isscalar(weight):
            w = float(weight)
            if w <= 0:
                raise HeightError('invalid %s: %.6f' % ('weight', w))
        else:
            n, w = len2(weight)
            if n != m:
                raise HeightError('invalid %s: %s, not %s' % (
                                  'number of weights', n, m))
            w = np.array(map(float, w))
            for i in range(m):
                if w[i] <= 0:
                    raise HeightError('invalid %s[%s]: %.6f' % (
                                      'weight', i, w[i]))

        T = 1.0e-4  # like SciPy example
        ps = np.array(_ordedup(xs, T, PI2 - T))
        ts = np.array(_ordedup(ys, T, PI  - T))

        try:
            self._ev = spi.LSQSphereBivariateSpline(ys, xs, hs,
                                                    ts, ps, eps=EPS, w=w).ev
        except Exception as x:
            raise _SciPyIssue(x)

    def __call__(self, *llis):
        '''Interpolate the height for one or several locations.

           @param llis: The location or locations (C{LatLon}, ... or
                        C{LatLon}s).

           @return: A single interpolated height (C{float}) or a list
                    or tuple of interpolated heights (C{float}s).

           @raise HeightError: Insufficient number of I{llis} or
                               invalid I{lli}.

           @raise SciPyError: A C{LSQSphereBivariateSpline} issue.

           @raise SciPyWarning: A C{LSQSphereBivariateSpline} warning
                                as exception.
        '''
        return _HeightBase._eval(self, llis)

    def height(self, lats, lons):
        '''Interpolate the height for one or several lat-/longitudes.

           @param lats: Latitude or latitudes (C{degrees} or C{degrees}s).
           @param lons: Longitude or longitudes (C{degrees} or C{degrees}s).

           @return: A single interpolated height (C{float}) or a list of
                    interpolated heights (C{float}s).

           @raise HeightError: Insufficient or non-matching number of
                               I{lats} and I{lons}.

           @raise SciPyError: A C{LSQSphereBivariateSpline} issue.

           @raise SciPyWarning: A C{LSQSphereBivariateSpline} warning
                                as exception.
        '''
        return _HeightBase._height(self, lats, lons)


class HeightSmoothBiSpline(_HeightBase):
    '''Height interpolator using C{SciPy} U{SmoothSphereBivariateSpline
       <http://docs.SciPy.org/doc/scipy/reference/generated/scipy.
       interpolate.SmoothSphereBivariateSpline.html>}.
    '''
    _kmin = 16  # k = 3, always

    def __init__(self, knots, s=4):
        '''New L{HeightSmoothBiSpline} interpolator.

           @param knots: The points with known height (C{LatLon}s).
           @keyword s: The spline smoothing factor (C{4}).

           @raise HeightError: Insufficient number of I{knots} or
                               invalid I{knot} or I{s}.

           @raise ImportError: Package C{numpy} or C{scipy} not found
                               or not installed.

           @raise SciPyError: A C{SmoothSphereBivariateSpline} issue.

           @raise SciPyWarning: A C{SmoothSphereBivariateSpline} warning
                                as exception.
        '''
        _, spi = self._NumSciPy()

        if s < 4:
            raise HeightError('%s too small: %s' % ('smoothing', s))

        xs, ys, hs = self._xyhs3(knots)
        try:
            self._ev = spi.SmoothSphereBivariateSpline(ys, xs, hs,
                                                       eps=EPS, s=s).ev
        except Exception as x:
            raise _SciPyIssue(x)

    def __call__(self, *llis):
        '''Interpolate the height for one or several locations.

           @param llis: The location or locations (C{LatLon}, ... or
                        C{LatLon}s).

           @return: A single interpolated height (C{float}) or a list
                    or tuple of interpolated heights (C{float}s).

           @raise HeightError: Insufficient number of I{llis} or
                               invalid I{lli}.

           @raise SciPyError: A C{SmoothSphereBivariateSpline} issue.

           @raise SciPyWarning: A C{SmoothSphereBivariateSpline} warning
                                as exception.
        '''
        return _HeightBase._eval(self, llis)

    def height(self, lats, lons):
        '''Interpolate the height for one or several lat-/longitudes.

           @param lats: Latitude or latitudes (C{degrees} or C{degrees}s).
           @param lons: Longitude or longitudes (C{degrees} or C{degrees}s).

           @return: A single interpolated height (C{float}) or a list of
                    interpolated heights (C{float}s).

           @raise HeightError: Insufficient or non-matching number of
                               I{lats} and I{lons}.

           @raise SciPyError: A C{SmoothSphereBivariateSpline} issue.

           @raise SciPyWarning: A C{SmoothSphereBivariateSpline} warning
                                as exception.
        '''
        return _HeightBase._height(self, lats, lons)

# **) MIT License
#
# Copyright (C) 2016-2019 -- mrJean1 at Gmail dot com
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
