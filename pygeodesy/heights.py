
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


The knots do not need to be ordered for any of the interpolators.

Errors from C{scipy} as raised as L{SciPyError}s.  Warnings issued
by C{scipy} can be thrown as L{SciPyWarning} exceptions, provided
Python C{warnings} are filtered accordingly, see L{SciPyWarning}.

@see: U{SciPy<http://docs.SciPy.org/doc/scipy/reference/interpolate.html>}.
'''

from fmath import EPS, fdot, fsum, isscalar, len2, map1
from formy import euclidean_, haversine_
from lazily import _ALL_LAZY
from points import LatLon_
from utily import PI, PI2, PI_2, radiansPI, radiansPI2

__all__ = _ALL_LAZY.heights
__version__ = '19.02.24'


class HeightError(ValueError):
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


def _ascalar(ais):
    # return single float, not numpy.float64
    ais = list(ais)  # np.array, etc. to list
    if len(ais) != 1:
        raise AssertionError('len[%r] = %s == 1' % (ais, len(ais)))
    return float(ais[0])  # remove np.<type>


def _atuple(ais):
    # return tuple of floats, not numpy.float64s
    return tuple(map(float, ais))


def _axyllis4(atype, llis, m=1, off=True):
    # convert lli C{LatLon}s to tuples or C{NumPy} arrays of C{SciPy} sphericals
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
        raise HeightError('insufficient %s: %s, need %s' % ('llis', n, m))
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


def _SciPyIssue(x):
    t = ' '.join(str(x).strip().split())
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
            yield (max(0.0, radiansPI2(ll.lon + 180.0)  - xf),
                   max(0.0, radiansPI( ll.lat +  90.0)) - yf, ll.height)
    except AttributeError as x:
        raise HeightError('%s: %r' % (x, ll))


def _xyhs3(atype, m, knots, off=True):
    # convert knot C{LatLon}s to tuples or C{NumPy} arrays and C{SciPy} sphericals
    xs, ys, hs = zip(*_xyhs(knots, off=off))  # PYCHOK unzip
    n = len(hs)
    if n < m:
        raise HeightError('insufficient %s: %s, need %s' % ('knots', n, m))
    return map1(atype, xs, ys, hs)


class _HeightBase(object):
    '''Interpolator base class.
    '''
    _kmin = 2  # min number of knots
    _np   = None  # numpy
    _sp   = None  # scipy

    def __call__(self, *unused):
        raise AssertionError('%s.%s not overloaded' % (self.__class__.__name__, '_()'))

    def _axyllis4(self, llis):
        return _axyllis4(self._np.array, llis)

    def _ev(self, *unused):
        raise AssertionError('%s.%s not overloaded' % (self.__class__.__name__, '_ev'))

    def _eval(self, llis):  # XXX not *llis
        _as, xis, yis, _ = self._axyllis4(llis)
        try:  # SciPy .ev signature: y first, then x!
            return _as(self._ev(yis, xis))
        except Exception as x:
            raise _SciPyIssue(x)

    def _height(self, lats, lons):
        if isscalar(lats) and isscalar(lons):
            llis = LatLon_(lats, lons)
        else:
            n, lats = len2(lats)
            m, lons = len2(lons)
            if n != m:
                raise HeightError('non-matching %s: %s vs %s' % ('len', n, m))
            llis = [LatLon_(*ll) for ll in zip(lats, lons)]
        return self(llis)

    def _NumSciPy(self, throwarnings=False):
        # import numpy and scipy
        if throwarnings:  # raise SciPyWarnings
            import warnings
            warnings.filterwarnings('error')
        import scipy.interpolate as sp
        import numpy as np

        self._np = np
        self._sp = sp
        return np, sp

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

           @raise ImportError: No C{numpy} or no C{scipy} module.

           @raise SciPyError: A C{scipy.interpolate.interp2d} issue.

           @raise SciPyWarning: A C{scipy.interpolate.interp2d} warning
                                as exception.
        '''
        _, sp = self._NumSciPy()

        xs, ys, hs = self._xyhs3(knots)
        try:  # SciPy.interpolate.interp2d kind 'linear' or 'cubic'
            self._interp2d = sp.interp2d(xs, ys, hs, kind=self._kind)
        except Exception as x:
            raise _SciPyIssue(x)

    def __call__(self, *llis):
        '''Interpolate the height for one or several locations.

           @param llis: The location or locations (C{LatLon}, ... or
                        C{LatLon}s).

           @return: A single interpolated height (C{float}) or a list
                    or tuple of interpolated heights (C{float}s).

           @raise HeightError: Insufficient number of I{llis} or
                               invalide I{lli}.

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

           @raise ImportError: No C{numpy} or no C{scipy} module.

           @raise SciPyError: A C{scipy.interpolate.interp2d} issue.

           @raise SciPyWarning: A C{scipy.interpolate.interp2d} warning
                                as exception.
        '''
        HeightCubic.__init__(self, knots)

    __call__ = HeightCubic.__call__
    height   = HeightCubic.height


class HeightIDW(_HeightBase):
    '''Height interpolator using U{Inverse Distance Weighting
       <http://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW).

       The distance is either the C{Euclidean} or C{Haversine} I{angular}
       distance from function L{euclidean_} respectively L{haversine_}.

       @see: U{IDW<http://www.Geo.FU-Berlin.De/en/v/soga/Geodata-analysis/
             geostatistics/Inverse-Distance-Weighting/index.html>},
             U{SHEPARD_INTERP_2D<http://People.SC.FSU.edu/~jburkardt/c_src/
             shepard_interp_2d/shepard_interp_2d.html>} and function L{euclidean_}.
    '''
    _adjust    = None
    _b         = 0  # negative distance power
    _distances = None
    _hs        = ()  # known heights
    _xs        = ()  # knot lons
    _ys        = ()  # knot lats

    def __init__(self, knots, adjust=True, beta=2):
        '''New L{HeightIDW} interpolator.

           @param knots: The points with known height (C{LatLon}s).
           @keyword adjust: Set I{adjust}=C{None} to use the L{haversine_}
                            distance.  Set I{adjust}=C{True} or C{False}
                            for the L{euclidean_} distance and adjust
                            the longitudinal delta by the cosine of the
                            mean latitude for I{adjust}=C{True}.
           @keyword beta: Inverse distance power (C{int} 1, 2, or 3).

           @raise HeightError: Insufficient number of I{knots} or invalid
                               I{knot}, I{adjust} or I{beta}.
        '''
        self._xs, self._ys, self._hs = _xyhs3(tuple, self._kmin, knots, off=False)

        if adjust in (True, False):
            self._adjust = adjust
            self._distances = self._euclidean
        elif adjust is None:
            self._distances = self._haversine
        else:
            raise HeightError('invalid %s=%r' % ('adjust', adjust))

        self._b = -int(beta)
        if not (0 < beta < 4 and beta == -self._b):
            raise HeightError('invalid %s=%r' % ('beta', beta))

    def _euclidean(self, x, y):
        for xk, yk in zip(self._xs, self._ys):
            yield euclidean_(y, yk, x - xk, adjust=self._adjust)

    def _haversine(self, x, y):
        for xk, yk in zip(self._xs, self._ys):
            yield haversine_(y, yk, x - xk)

    def _hIDW(self, x, y):
        ws = tuple(self._distances(x, y))
        w, h = min(zip(ws, self._hs))
        if w > EPS:
            ws = tuple(w**self._b for w in ws)
            h = fdot(ws, *self._hs) / fsum(ws)
        return h

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

           @raise ImportError: No C{numpy} or no C{scipy} module.

           @raise SciPyError: A C{LSQSphereBivariateSpline} issue.

           @raise SciPyWarning: A C{LSQSphereBivariateSpline} warning
                                as exception..
        '''
        np, sp = self._NumSciPy()

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
            self._ev = sp.LSQSphereBivariateSpline(ys, xs, hs,
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

           @raise ImportError: No C{numpy} or no C{scipy} module.

           @raise SciPyError: A C{SmoothSphereBivariateSpline} issue.

           @raise SciPyWarning: A C{SmoothSphereBivariateSpline} warning
                                as exception.
        '''
        _, sp = self._NumSciPy()

        if s < 4:
            raise HeightError('%s too small: %s' % ('smoothing', s))

        xs, ys, hs = self._xyhs3(knots)
        try:
            self._ev = sp.SmoothSphereBivariateSpline(ys, xs, hs, eps=EPS, s=s).ev
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
