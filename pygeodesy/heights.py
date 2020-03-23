
# -*- coding: utf-8 -*-

u'''Classes L{HeightCubic}, L{HeightIDWequirectangular},
L{HeightIDWeuclidean}, L{HeightIDWhaversine}, L{HeightIDWkarney},
L{HeightIDWvincentys}, L{HeightLinear}, L{HeightLSQBiSpline} and
L{HeightSmoothBiSpline} to interpolate the height of C{LatLon}
locations or separate lat-/longitudes from a set of C{LatLon}
points with known heights.

Except for L{HeightIDWequirectangular}, L{HeightIDWeuclidean},
L{HeightIDWhaversine} and L{HeightIDWvincentys}, the height
interpolators in this module require the packages U{numpy
<https://PyPI.org/project/numpy>} and U{scipy<https://SciPy.org>}
or U{geographiclib<https://PyPI.org/project/geographiclib>} to be
installed.

Typical usage is as follows.  First create an interpolator from a
given set of C{LatLon} points with known heights, called C{knots}.

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


The C{knots} do not need to be ordered for any of the height
interpolators.

Errors from C{scipy} as raised as L{SciPyError}s.  Warnings issued
by C{scipy} can be thrown as L{SciPyWarning} exceptions, provided
Python C{warnings} are filtered accordingly, see L{SciPyWarning}.

@see: U{SciPy<https://docs.SciPy.org/doc/scipy/reference/interpolate.html>}.
'''

from pygeodesy.basics import EPS, _isnotError, isscalar, \
                             len2, map1, map2, property_RO
from pygeodesy.datum import Datum
from pygeodesy.fmath import fidw, hypot2
from pygeodesy.formy import euclidean_, haversine_, _scaler, vincentys_
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY
from pygeodesy.named import _Named, notOverloaded
from pygeodesy.points import LatLon_
from pygeodesy.utily import PI, PI2, PI_2, radiansPI, radiansPI2, \
                            unroll180, unrollPI

__all__ = _ALL_LAZY.heights + _ALL_DOCS('_HeightBase')
__version__ = '20.03.23'


class HeightError(ValueError):  # imported by .geoids
    '''Height interpolator C{Height...} or interpolation issue.
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
       <https://docs.Python.org/3/library/warnings.html#the-warnings-filter>}
       I{prior to} C{import scipy} or by setting environment variable
       U{PYTHONWARNINGS<https://docs.Python.org/3/using/cmdline.html
       #envvar-PYTHONWARNINGS>} or with C{python} command line option
       U{-W<https://docs.Python.org/3/using/cmdline.html#cmdoption-w>}
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


def _SciPyIssue(x, *extras):  # imported by .geoids  # PYCHOK no cover
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


class _HeightBase(_Named):  # imported by .geoids
    '''(INTERNAL) Interpolator base class.
    '''
    _adjust = None  # not applicable
    _kmin   = 2     # min number of knots
    _np     = None  # numpy
    _np_v   = None  # version
    _spi    = None  # scipy.interpolate
    _sp_v   = None  # version
    _wrap   = None  # not applicable

    def __call__(self, *args):
        '''(INTERNAL) I{Must be overloaded}.
        '''
        notOverloaded(self, '__call__', *args)

    @property_RO
    def adjust(self):
        '''Get the adjust setting (C{bool} or C{None} if not applicable).
        '''
        return self._adjust

    def _axyllis4(self, llis):
        return _axyllis4(self._np.array, llis)

    def _ev(self, *args):
        '''(INTERNAL) I{Must be overloaded}.
        '''
        notOverloaded(self, self._ev.__name__, *args)

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

    @property_RO
    def kmin(self):
        '''Get the minimum number of knots (C{int}).
        '''
        return self._kmin

    def _NumSciPy(self, throwarnings=False):
        '''(INTERNAL) Import C{numpy} and C{scipy}.
        '''
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

    @property_RO
    def wrap(self):
        '''Get the wrap setting (C{bool} or C{None} if not applicable).
        '''
        return self._wrap


class HeightCubic(_HeightBase):
    '''Height interpolator based on C{SciPy} U{interp2d<https://docs.SciPy.org/
       doc/scipy/reference/generated/scipy.interpolate.interp2d.html>}
       C{kind='cubic'}.
    '''
    _interp2d = None
    _kind     = 'cubic'
    _kmin     = 16

    def __init__(self, knots, name=''):
        '''New L{HeightCubic} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg name: Optional height interpolator name (C{str}).

           @raise HeightError: Insufficient number of B{C{knots}} or
                               invalid B{C{knot}}.

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

        if name:
            self.name = name

    def __call__(self, *llis):
        '''Interpolate the height for one or several locations.

           @arg llis: The location or locations (C{LatLon}, ... or
                      C{LatLon}s).

           @return: A single interpolated height (C{float}) or a list
                    or tuple of interpolated heights (C{float}s).

           @raise HeightError: Insufficient number of B{C{llis}} or
                               invalid B{C{lli}}.

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

           @arg lats: Latitude or latitudes (C{degrees} or C{degrees}s).
           @arg lons: Longitude or longitudes (C{degrees} or C{degrees}s).

           @return: A single interpolated height (C{float}) or a list of
                    interpolated heights (C{float}s).

           @raise HeightError: Insufficient or non-matching number of
                               B{C{lats}} and B{C{lons}}.

           @raise SciPyError: A C{scipy.interpolate.interp2d} issue.

           @raise SciPyWarning: A C{scipy.interpolate.interp2d} warning
                                as exception.
        '''
        return _HeightBase._height(self, lats, lons)


class HeightLinear(HeightCubic):
    '''Height interpolator based on C{SciPy} U{interp2d<https://docs.SciPy.org/
       doc/scipy/reference/generated/scipy.interpolate.interp2d.html>}
       C{kind='linear}.
    '''
    _kind = 'linear'
    _kmin = 2

    def __init__(self, knots, name=''):
        '''New L{HeightLinear} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg name: Optional height interpolator name (C{str}).

           @raise HeightError: Insufficient number of B{C{knots}} or
                               invalid B{C{knot}}.

           @raise ImportError: Package C{numpy} or C{scipy} not found
                               or not installed.

           @raise SciPyError: A C{scipy.interpolate.interp2d} issue.

           @raise SciPyWarning: A C{scipy.interpolate.interp2d} warning
                                as exception.
        '''
        HeightCubic.__init__(self, knots, name=name)

    __call__ = HeightCubic.__call__  # for __doc__
    height   = HeightCubic.height    # for __doc__


class _HeightIDW(_HeightBase):
    '''(INTERNAL) Base class for U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW)
       height interpolators.
    '''
    _beta = 0   # inverse distance power
    _hs   = ()  # known heights
    _xs   = ()  # knot lons
    _ys   = ()  # knot lats

    def __init__(self, knots, beta=2, name=''):
        '''New L{_HeightIDW} interpolator.
        '''
        self._xs, self._ys, self._hs = _xyhs3(tuple, self._kmin, knots, off=False)
        self.beta = beta
        if name:
            self.name = name

    def __call__(self, *llis):
        '''Interpolate the height for one or several locations.

           @arg llis: The location or locations (C{LatLon}, ... or
                      C{LatLon}s).

           @return: A single interpolated height (C{float}) or a list
                    or tuple of interpolated heights (C{float}s).

           @raise HeightError: Insufficient number of B{C{llis}},
                               invalid B{C{lli}} or an L{fidw}
                               issue.
        '''
        _as, xis, yis, _ = _axyllis4(tuple, llis, off=False)
        return _as(map(self._hIDW, xis, yis))

    def _distances(self, x, y):  # PYCHOK unused (x, y) radians
        '''Must be overloaded.
        '''
        raise NotImplementedError('method: %s' % ('_distances',))

    def _hIDW(self, x, y):
        # interpolate height at (x, y) radians or degrees
        try:
            ds = self._distances(x, y)
            return fidw(self._hs, ds, beta=self._beta)
        except ValueError as x:
            raise HeightError(str(x))

    @property
    def beta(self):
        '''Get the inverse distance power (C{int}).
        '''
        return self._beta

    @beta.setter  # PYCHOK setter!
    def beta(self, beta):
        '''Set the inverse distance power.

           @arg beta: New inverse distance power (C{int} 1, 2, or 3).

           @raise HeightError: Invalid B{C{beta}}.
        '''
        try:
            b = int(beta)
            if b != beta or 1 > b or b > 3:
                raise ValueError
        except (TypeError, ValueError):
            raise HeightError('%s invalid: %r' % ('beta', beta))
        self._beta = b

    def height(self, lats, lons):
        '''Interpolate the height for one or several lat-/longitudes.

           @arg lats: Latitude or latitudes (C{degrees} or C{degrees}s).
           @arg lons: Longitude or longitudes (C{degrees} or C{degrees}s).

           @return: A single interpolated height (C{float}) or a list of
                    interpolated heights (C{float}s).

           @raise HeightError: Insufficient or non-matching number of
                               B{C{lats}} and B{C{lons}} or an L{fidw}
                               issue.
        '''
        return _HeightBase._height(self, lats, lons)


class HeightIDWequirectangular(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW)
       and the C{equirectangular} distance (in radians squared) like
       function L{equirectangular_}.

       @see: L{HeightIDWeuclidean}, L{HeightIDWhaversine},
             L{HeightIDWvincentys}, U{Inverse distance weighting
             <https://WikiPedia.org/wiki/Inverse_distance_weighting>},
             U{IDW<https://www.Geo.FU-Berlin.De/en/v/soga/Geodata-analysis/
             geostatistics/Inverse-Distance-Weighting/index.html>} and
             U{SHEPARD_INTERP_2D<https://People.SC.FSU.edu/~jburkardt/c_src/
             shepard_interp_2d/shepard_interp_2d.html>}.
    '''
    _adjust = True
    _wrap   = False

    def __init__(self, knots, adjust=True, wrap=False, name=''):
        '''New L{HeightIDWequirectangular} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg adjust: Adjust the wrapped, unrolled longitudinal
                          delta by the cosine of the mean latitude (C{bool}).
           @kwarg wrap: Wrap and L{unroll180} longitudes (C{bool}).
           @kwarg name: Optional height interpolator name (C{str}).

           @raise HeightError: Insufficient number of B{C{knots}} or
                               invalid B{C{knot}}.
        '''
        if not adjust:
            self._adjust = False
        if wrap:
            self._wrap = True
        _HeightIDW.__init__(self, knots, beta=1, name=name)  # distance**2

    def _distances(self, x, y):  # (x, y) radians
        for xk, yk in zip(self._xs, self._ys):
            d, _ = unrollPI(xk, x, wrap=self._wrap)
            if self.adjust:
                d *= _scaler(yk, y)
            yield hypot2(d, yk - y)  # like equirectangular_ distance2

    __call__ = _HeightIDW.__call__  # for __doc__
    height   = _HeightIDW.height    # for __doc__


class HeightIDWeuclidean(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW)
       and the C{Euclidean} distance from function L{euclidean_}.

       @see: L{HeightIDWequirectangular}, L{HeightIDWhaversine},
             L{HeightIDWkarney}, L{HeightIDWvincentys}, U{Inverse distance
             weighting<https://WikiPedia.org/wiki/Inverse_distance_weighting>},
             U{IDW<https://www.Geo.FU-Berlin.De/en/v/soga/Geodata-analysis/
             geostatistics/Inverse-Distance-Weighting/index.html>} and
             U{SHEPARD_INTERP_2D<https://People.SC.FSU.edu/~jburkardt/c_src/
             shepard_interp_2d/shepard_interp_2d.html>}.
    '''
    _adjust = True

    def __init__(self, knots, adjust=True, beta=2, name=''):
        '''New L{HeightIDWeuclidean} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg adjust: Adjust the longitudinal delta by the cosine
                          of the mean latitude for B{C{adjust}}=C{True}.
           @kwarg beta: Inverse distance power (C{int} 1, 2, or 3).
           @kwarg name: Optional height interpolator name (C{str}).

           @raise HeightError: Insufficient number of B{C{knots}} or
                               invalid B{C{knot}} or B{C{beta}}.
        '''
        if not adjust:
            self._adjust = False
        _HeightIDW.__init__(self, knots, beta=beta, name=name)

    def _distances(self, x, y):  # (x, y) radians
        for xk, yk in zip(self._xs, self._ys):
            yield euclidean_(yk, y, xk - x, adjust=self.adjust)

    __call__ = _HeightIDW.__call__  # for __doc__
    height   = _HeightIDW.height    # for __doc__


class HeightIDWhaversine(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW) and
       the I{angular} C{Haversine} distance from function L{haversine_}.

       @note: See note under L{HeightIDWvincentys}.

       @see: L{HeightIDWequirectangular}, L{HeightIDWeuclidean},
             L{HeightIDWkarney}, L{HeightIDWvincentys}, U{Inverse distance
             weighting<https://WikiPedia.org/wiki/Inverse_distance_weighting>},
             U{IDW<https://www.Geo.FU-Berlin.De/en/v/soga/Geodata-analysis/
             geostatistics/Inverse-Distance-Weighting/index.html>} and
             U{SHEPARD_INTERP_2D<https://People.SC.FSU.edu/~jburkardt/c_src/
             shepard_interp_2d/shepard_interp_2d.html>}.
    '''
    _wrap = False

    def __init__(self, knots, beta=2, wrap=False, name=''):
        '''New L{HeightIDWhaversine} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg beta: Inverse distance power (C{int} 1, 2, or 3).
           @kwarg wrap: Wrap and L{unroll180} longitudes (C{bool}).
           @kwarg name: Optional height interpolator name (C{str}).

           @raise HeightError: Insufficient number of B{C{knots}} or
                               invalid B{C{knot}} or B{C{beta}}.
        '''
        if wrap:
            self._wrap = True
        _HeightIDW.__init__(self, knots, beta=beta, name=name)

    def _distances(self, x, y):  # (x, y) radians
        for xk, yk in zip(self._xs, self._ys):
            d, _ = unrollPI(xk, x, wrap=self._wrap)
            yield haversine_(yk, y, d)

    __call__ = _HeightIDW.__call__  # for __doc__
    height   = _HeightIDW.height    # for __doc__


class HeightIDWkarney(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW) and
       the I{angular} distance from I{Charles F. F. Karney's}
       U{GeographicLib<https://PyPI.org/project/geographiclib>} U{Geodesic
       <https://geographiclib.sourceforge.io/1.49/python/code.html>}
       Inverse method.

       @see: L{HeightIDWequirectangular}, L{HeightIDWeuclidean},
             L{HeightIDWhaversine}, L{HeightIDWvincentys}, U{Inverse distance
             weighting<https://WikiPedia.org/wiki/Inverse_distance_weighting>},
             U{IDW<https://www.Geo.FU-Berlin.De/en/v/soga/Geodata-analysis/
             geostatistics/Inverse-Distance-Weighting/index.html>} and
             U{SHEPARD_INTERP_2D<https://People.SC.FSU.edu/~jburkardt/c_src/
             shepard_interp_2d/shepard_interp_2d.html>}.
    '''
    _datum    = None
    _geodesic = None
    _wrap     = False

    def __init__(self, knots, datum=None, beta=2, wrap=False, name=''):
        '''New L{HeightIDWhaversine} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg datum: Optional datum (L{Datum} to use, overriding
                         the default B{C{knots[0].datum}}
           @kwarg beta: Inverse distance power (C{int} 1, 2, or 3).
           @kwarg wrap: Wrap and L{unroll180} longitudes (C{bool}).
           @kwarg name: Optional height interpolator name (C{str}).

           @raise HeightError: Insufficient number of B{C{knots}} or
                               invalid B{C{knot}}, B{C{datum}} or
                               B{C{beta}}.

           @raise ImportError: Package U{GeographicLib
                  <https://PyPI.org/project/geographiclib>} missing.

           @raise TypeError: Invalid B{C{datum}}.
        '''
        n, self._lls = len2(knots)
        if n < self._kmin:
            raise HeightError('insufficient %s: %s, need %s' % ('knots', n, self._kmin))
        try:
            self._datum = self._lls[0].datum if datum is None else datum
            if not isinstance(self.datum, Datum):
                raise TypeError
        except (AttributeError, TypeError):
            raise _isnotError('valid', datum=self.datum or datum)
        self._geodesic = self.datum.ellipsoid.geodesic

        self.beta = beta
        if wrap:
            self._wrap = True
        if name:
            self.name = name

    def _distances(self, x, y):  # (x, y) degrees
        g = self._geodesic
        for ll in self._lls:
            # see .ellipsoidalKarney.LatLon._inverse
            _, lon = unroll180(x, ll.lon, wrap=self._wrap)  # g.LONG_UNROLL
            # XXX g.DISTANCE needed for 's12', distance in meters?
            yield abs(g.Inverse(y, x, ll.lat, lon)['a12'])

    @property_RO
    def _hs(self):
        for ll in self._lls:
            yield ll.height

    def __call__(self, *llis):
        '''Interpolate the height for one or several locations.

           @arg llis: The location or locations (C{LatLon}, ... or
                      C{LatLon}s).

           @return: A single interpolated height (C{float}) or a list
                    or tuple of interpolated heights (C{float}s).

           @raise HeightError: Insufficient number of B{C{llis}},
                               invalid B{C{lli}} or an L{fidw}
                               issue.
        '''
        def _xy2(lls):
            try:  # like _xyhs above, but keeping degrees
                for ll in lls:
                    yield ll.lon, ll.lat
            except AttributeError as x:
                raise HeightError('%s: %r' % (x, ll))

        _as, llis = _allis2(llis)
        return _as(map(self._hIDW, *zip(*_xy2(llis))))

    @property_RO
    def datum(self):
        '''Get the datum of this interpolator (L{Datum}).
        '''
        return self._datum

    height = _HeightIDW.height  # for __doc__


class HeightIDWvincentys(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW)
       and the I{angular} C{Vincenty} distance from function L{vincentys_}.

       @note: See note under L{vincentys_}.

       @see: L{HeightIDWequirectangular}, L{HeightIDWeuclidean},
             L{HeightIDWhaversine}, L{HeightIDWkarney}, U{Inverse distance
             weighting<https://WikiPedia.org/wiki/Inverse_distance_weighting>},
             U{IDW<https://www.Geo.FU-Berlin.De/en/v/soga/Geodata-analysis/
             geostatistics/Inverse-Distance-Weighting/index.html>} and
             U{SHEPARD_INTERP_2D<https://People.SC.FSU.edu/~jburkardt/c_src/
             shepard_interp_2d/shepard_interp_2d.html>}.
    '''
    _wrap = False

    def __init__(self, knots, beta=2, wrap=False, name=''):
        '''New L{HeightIDWvincentys} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg beta: Inverse distance power (C{int} 1, 2, or 3).
           @kwarg wrap: Wrap and L{unroll180} longitudes (C{bool}).
           @kwarg name: Optional height interpolator name (C{str}).

           @raise HeightError: Insufficient number of B{C{knots}} or
                               invalid B{C{knot}} or B{C{beta}}.
        '''
        if wrap:
            self._wrap = True
        _HeightIDW.__init__(self, knots, beta=beta, name=name)

    def _distances(self, x, y):  # (x, y) radians
        for xk, yk in zip(self._xs, self._ys):
            d, _ = unrollPI(xk, x, wrap=self._wrap)
            yield vincentys_(yk, y, d)

    __call__ = _HeightIDW.__call__  # for __doc__
    height   = _HeightIDW.height    # for __doc__


class HeightLSQBiSpline(_HeightBase):
    '''Height interpolator using C{SciPy} U{LSQSphereBivariateSpline
       <https://docs.SciPy.org/doc/scipy/reference/generated/scipy.
       interpolate.LSQSphereBivariateSpline.html>}.
    '''
    _kmin = 16  # k = 3, always

    def __init__(self, knots, weight=None, name=''):
        '''New L{HeightLSQBiSpline} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg weight: Optional weight or weights for each B{C{knot}}
                          (C{scalar} or C{scalar}s).
           @kwarg name: Optional height interpolator name (C{str}).

           @raise HeightError: Insufficient number of B{C{knots}} or
                               B{C{weight}}s or invalid B{C{knot}} or B{C{weight}}.

           @raise ImportError: Package C{numpy} or C{scipy} not found
                               or not installed.

           @raise SciPyError: A C{LSQSphereBivariateSpline} issue.

           @raise SciPyWarning: A C{LSQSphereBivariateSpline} warning
                                as exception.
        '''
        np, spi = self._NumSciPy()

        xs, ys, hs = self._xyhs3(knots)
        n = len(hs)

        w = weight
        if isscalar(w):
            w = float(w)
            if w <= 0:
                raise HeightError('%s invalid: %.6F' % ('weight', w))
            w = [w] * n
        elif w is not None:
            m, w = len2(w)
            if m != n:
                raise HeightError('%s invalid: %s, not %s' % (
                                  'number of weights', m, n))
            w = map2(float, w)
            m = min(w)
            if m <= 0:
                raise HeightError('%s[%s] invalid: %.6F' % (
                                  'weight', w.find(m), m))
        try:
            T = 1.0e-4  # like SciPy example
            ps = np.array(_ordedup(xs, T, PI2 - T))
            ts = np.array(_ordedup(ys, T, PI  - T))
            self._ev = spi.LSQSphereBivariateSpline(ys, xs, hs,
                                                    ts, ps, eps=EPS, w=w).ev
        except Exception as x:
            raise _SciPyIssue(x)

        if name:
            self.name = name

    def __call__(self, *llis):
        '''Interpolate the height for one or several locations.

           @arg llis: The location or locations (C{LatLon}, ... or
                      C{LatLon}s).

           @return: A single interpolated height (C{float}) or a list
                    or tuple of interpolated heights (C{float}s).

           @raise HeightError: Insufficient number of B{C{llis}} or
                               invalid B{C{lli}}.

           @raise SciPyError: A C{LSQSphereBivariateSpline} issue.

           @raise SciPyWarning: A C{LSQSphereBivariateSpline} warning
                                as exception.
        '''
        return _HeightBase._eval(self, llis)

    def height(self, lats, lons):
        '''Interpolate the height for one or several lat-/longitudes.

           @arg lats: Latitude or latitudes (C{degrees} or C{degrees}s).
           @arg lons: Longitude or longitudes (C{degrees} or C{degrees}s).

           @return: A single interpolated height (C{float}) or a list of
                    interpolated heights (C{float}s).

           @raise HeightError: Insufficient or non-matching number of
                               B{C{lats}} and B{C{lons}}.

           @raise SciPyError: A C{LSQSphereBivariateSpline} issue.

           @raise SciPyWarning: A C{LSQSphereBivariateSpline} warning
                                as exception.
        '''
        return _HeightBase._height(self, lats, lons)


class HeightSmoothBiSpline(_HeightBase):
    '''Height interpolator using C{SciPy} U{SmoothSphereBivariateSpline
       <https://docs.SciPy.org/doc/scipy/reference/generated/scipy.
       interpolate.SmoothSphereBivariateSpline.html>}.
    '''
    _kmin = 16  # k = 3, always

    def __init__(self, knots, s=4, name=''):
        '''New L{HeightSmoothBiSpline} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg s: The spline smoothing factor (C{4}).
           @kwarg name: Optional height interpolator name (C{str}).

           @raise HeightError: Insufficient number of B{C{knots}} or
                               invalid B{C{knot}} or B{C{s}}.

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

        if name:
            self.name = name

    def __call__(self, *llis):
        '''Interpolate the height for one or several locations.

           @arg llis: The location or locations (C{LatLon}, ... or
                      C{LatLon}s).

           @return: A single interpolated height (C{float}) or a list
                    or tuple of interpolated heights (C{float}s).

           @raise HeightError: Insufficient number of B{C{llis}} or
                               invalid B{C{lli}}.

           @raise SciPyError: A C{SmoothSphereBivariateSpline} issue.

           @raise SciPyWarning: A C{SmoothSphereBivariateSpline} warning
                                as exception.
        '''
        return _HeightBase._eval(self, llis)

    def height(self, lats, lons):
        '''Interpolate the height for one or several lat-/longitudes.

           @arg lats: Latitude or latitudes (C{degrees} or C{degrees}s).
           @arg lons: Longitude or longitudes (C{degrees} or C{degrees}s).

           @return: A single interpolated height (C{float}) or a list of
                    interpolated heights (C{float}s).

           @raise HeightError: Insufficient or non-matching number of
                               B{C{lats}} and B{C{lons}}.

           @raise SciPyError: A C{SmoothSphereBivariateSpline} issue.

           @raise SciPyWarning: A C{SmoothSphereBivariateSpline} warning
                                as exception.
        '''
        return _HeightBase._height(self, lats, lons)

# **) MIT License
#
# Copyright (C) 2016-2020 -- mrJean1 at Gmail -- All Rights Reserved.
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
