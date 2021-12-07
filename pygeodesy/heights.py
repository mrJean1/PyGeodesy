
# -*- coding: utf-8 -*-

u'''Height interpolations of C{LatLon} points.

Classes L{HeightCubic}, L{HeightIDWcosineAndoyerLambert},
L{HeightIDWcosineForsytheAndoyerLambert}, L{HeightIDWcosineLaw},
L{HeightIDWdistanceTo}, L{HeightIDWequirectangular}, L{HeightIDWeuclidean},
L{HeightIDWflatLocal}, L{HeightIDWflatPolar}, L{HeightIDWhaversine},
L{HeightIDWhubeny}, L{HeightIDWkarney}, L{HeightIDWthomas}, L{HeightIDWvincentys},
L{HeightLinear}, L{HeightLSQBiSpline} and L{HeightSmoothBiSpline}
to interpolate the height of C{LatLon} locations or separate
lat-/longitudes from a set of C{LatLon} points with I{known heights}.

Typical usage
=============

1. Get or create a set of C{LatLon} points with I{known heights},
called C{knots}.  The C{knots} do not need to be ordered in any
particular way.

2. Select one of the C{Height} classes for height interpolation

C{>>> from pygeodesy import HeightCubic  # or other Height... as HeightXyz}

3. Instantiate a height interpolator with the C{knots} and use keyword
arguments to select different interpolation options

C{>>> hinterpolator = HeightXyz(knots, **options)}

4. Get the interpolated height of other C{LatLon} location(s) with

C{>>> h = hinterpolator(ll)}

or

C{>>> h0, h1, h2, ... = hinterpolator(ll0, ll1, ll2, ...)}

or a list, tuple, generator, etc. of C{LatLon}s

C{>>> hs = hinterpolator(lls)}

5. For separate lat- and longitudes invoke the C{.height} method

C{>>> h = hinterpolator.height(lat, lon)}

or as 2 lists, tuples, etc.

C{>>> hs = hinterpolator.height(lats, lons)}

@note: Classes L{HeightCubic} and L{HeightLinear} require package U{numpy
       <https://PyPI.org/project/numpy>}, classes L{HeightLSQBiSpline} and
       L{HeightSmoothBiSpline} require package U{scipy<https://SciPy.org>}.
       Classes L{HeightIDWkarney} and L{HeightIDWdistanceTo} -if used with
       L{ellipsoidalKarney.LatLon} points- require I{Karney}'s U{geographiclib
       <https://PyPI.org/project/geographiclib>} to be installed.

@note: Errors from C{scipy} are raised as L{SciPyError}s.  Warnings issued
       by C{scipy} can be thrown as L{SciPyWarning} exceptions, provided
       Python C{warnings} are filtered accordingly, see L{SciPyWarning}.

@see: U{SciPy<https://docs.SciPy.org/doc/scipy/reference/interpolate.html>}
      Interpolation.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import isscalar, len2, map1, map2, _xnumpy, _xscipy
from pygeodesy.datums import _ellipsoidal_datum, _WGS84
from pygeodesy.errors import _AssertionError, LenError, PointsError, _SciPyIssue
from pygeodesy.fmath import fidw, hypot2
from pygeodesy.formy import cosineAndoyerLambert_, cosineForsytheAndoyerLambert_, \
                            cosineLaw_, euclidean_, flatPolar_, haversine_, \
                            _scale_rad, thomas_, vincentys_
from pygeodesy.interns import EPS, NN, PI, PI2, PI_2, _COMMASPACE_, _cubic_, \
                             _datum_, _distanceTo_, _knots_, _len_, _linear_, \
                             _scipy_, _0_0, _90_0, _180_0
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _FOR_DOCS
from pygeodesy.named import _Named, notOverloaded
from pygeodesy.points import LatLon_
from pygeodesy.props import Property_RO
from pygeodesy.streprs import _boolkwds, Fmt
from pygeodesy.units import Float_, Int_
from pygeodesy.utily import radiansPI, radiansPI2, unrollPI

__all__ = _ALL_LAZY.heights
__version__ = '21.11.30'

_error_        = 'error'
_insufficient_ = 'insufficient'
_llis_         = 'llis'


class HeightError(PointsError):
    '''Height interpolator C{Height...} or interpolation issue.
    '''
    pass


def _alist(ais):
    # return list of floats, not numpy.float64s
    return list(map(float, ais))


def _allis2(llis, m=1, Error=HeightError):  # imported by .geoids
    # dtermine return type and convert lli C{LatLon}s to list
    if not isinstance(llis, tuple):  # llis are *args
        raise _AssertionError('%s(%s): %r' % ('type_', '*llis', llis))

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
        raise _insufficientError(m, Error=Error, llis=n)
    return _as, llis


def _ascalar(ais):  # imported by .geoids
    # return single float, not numpy.float64
    ais = list(ais)  # np.array, etc. to list
    if len(ais) != 1:
        raise _AssertionError('%s(%r): %s != 1' % (_len_, ais, len(ais)))
    return float(ais[0])  # remove np.<type>


def _atuple(ais):
    # return tuple of floats, not numpy.float64s
    return tuple(map(float, ais))


def _axyllis4(atype, llis, m=1, off=True):
    # convert lli C{LatLon}s to tuples or C{NumPy} arrays of
    # C{SciPy} sphericals and determine the return type
    _as, llis = _allis2(llis, m=m)
    xis, yis, _ =  zip(*_xyhs(llis, off=off))  # PYCHOK unzip
    return _as, atype(xis), atype(yis), llis


def _insufficientError(need, Error=HeightError, **name_value):
    # create an insufficient Error instance
    t = _COMMASPACE_(_insufficient_, str(need))
    return Error(txt=t, **name_value)


def _ordedup(ts, lo=EPS, hi=PI2-EPS):
    # clip, order and remove duplicates
    # p, ks = 0, []
    # for k in sorted(max(lo, min(hi, t)) for t in ts):
    #     if k > p:
    #         ks.append(k)
    #         p = k
    # return ks
    return sorted(set(max(lo, min(hi, t)) for t in ts))  # list


def _xyhs(lls, off=True, name=_llis_):
    # map (lat, lon, h) to (x, y, h) in radians, offset as
    # x: 0 <= lon <= PI2, y: 0 <= lat <= PI if off is True
    # else x: -PI <= lon <= PI, y: -PI_2 <= lat <= PI_2
    if off:
        xf = yf = _0_0
    else:  # undo offset
        xf, yf = PI, PI_2
    try:
        for i, ll in enumerate(lls):
            yield (max(_0_0, radiansPI2(ll.lon + _180_0)) - xf), \
                  (max(_0_0, radiansPI( ll.lat +  _90_0)) - yf), ll.height
    except AttributeError as x:
        i = Fmt.SQUARE(name, i)
        raise HeightError(i, ll, txt=str(x))


def _xyhs3(atype, m, knots, off=True):
    # convert knot C{LatLon}s to tuples or C{NumPy} arrays and C{SciPy} sphericals
    xs, ys, hs = zip(*_xyhs(knots, off=off, name=_knots_))  # PYCHOK unzip
    n = len(hs)
    if n < m:
        raise _insufficientError(m, knots=n)
    return map1(atype, xs, ys, hs)


class _HeightBase(_Named):  # imported by .geoids
    '''(INTERNAL) Interpolator base class.
    '''
    _adjust = None     # not applicable
    _datum  = None     # ._height datum
    _kmin   = 2        # min number of knots
    _LLis   = LatLon_  # ._height class
    _np_sp  = None     # (numpy, scipy)
    _wrap   = None     # not applicable

    def __call__(self, *args):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}, see function C{notOverloaded}.
        '''
        notOverloaded(self, callername='__call__', *args)

    @Property_RO
    def adjust(self):
        '''Get the adjust setting (C{bool} or C{None} if not applicable).
        '''
        return self._adjust

    def _axyllis4(self, llis):
        return _axyllis4(self.numpy.array, llis)

    @Property_RO
    def datum(self):
        '''Get the datum (L{Datum} or C{None} if not applicable).
        '''
        return self._datum

    def _ev(self, *args):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}, see function C{notOverloaded}.
        '''
        notOverloaded(self, *args)

    def _eval(self, llis):  # XXX single arg, not *args
        _as, xis, yis, _ = self._axyllis4(llis)
        try:  # SciPy .ev signature: y first, then x!
            return _as(self._ev(yis, xis))
        except Exception as x:
            raise _SciPyIssue(x)

    def _height(self, lats, lons, Error=HeightError):
        LLis, d = self._LLis, self.datum
        if isscalar(lats) and isscalar(lons):
            llis = LLis(lats, lons, datum=d)
        else:
            n, lats = len2(lats)
            m, lons = len2(lons)
            if n != m:
                # format a LenError, but raise an Error
                e = LenError(self.__class__, lats=n, lons=m, txt=None)
                raise e if Error is LenError else Error(str(e))
            llis = [LLis(*ll, datum=d) for ll in zip(lats, lons)]
        return self(llis)  # __call__(lli) or __call__(llis)

    @Property_RO
    def kmin(self):
        '''Get the minimum number of knots (C{int}).
        '''
        return self._kmin

    def _np_sp2(self, throwarnings=False):
        '''(INTERNAL) Import C{numpy} and C{scipy}, once.
        '''
        t = _HeightBase._np_sp
        if not t:
            # raise SciPyWarnings, but not if
            # scipy has been imported already
            if throwarnings:
                import sys
                if _scipy_ not in sys.modules:
                    import warnings
                    warnings.filterwarnings(_error_)

            sp = _xscipy(self.__class__, 1, 2)
            np = _xnumpy(self.__class__, 1, 9)

            _HeightBase._np_sp = t = np, sp
        return t

    @Property_RO
    def numpy(self):
        '''Get the C{numpy} module or C{None}.
        '''
        np, _ = self._np_sp2()
        return np

    @Property_RO
    def scipy(self):
        '''Get the C{scipy} module or C{None}.
        '''
        _, sp = self._np_sp2()
        return sp

    @Property_RO
    def scipy_interpolate(self):
        '''Get the C{scipy.interpolate} module or C{None}.
        '''
        _ = self.scipy
        import scipy.interpolate as spi
        return spi

    def _xyhs3(self, knots):
        return _xyhs3(self.numpy.array, self._kmin, knots)

    @Property_RO
    def wrap(self):
        '''Get the wrap setting (C{bool} or C{None} if not applicable).
        '''
        return self._wrap


class HeightCubic(_HeightBase):
    '''Height interpolator based on C{SciPy} U{interp2d<https://docs.SciPy.org/
       doc/scipy/reference/generated/scipy.interpolate.interp2d.html>}
       C{kind='cubic'}.
    '''
    _interp2d =  None
    _kind     = _cubic_
    _kmin     =  16

    def __init__(self, knots, name=NN):
        '''New L{HeightCubic} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg name: Optional name for this height interpolator (C{str}).

           @raise HeightError: Insufficient number of B{C{knots}} or
                               invalid B{C{knot}}.

           @raise ImportError: Package C{numpy} or C{scipy} not found
                               or not installed.

           @raise SciPyError: A C{scipy.interpolate.interp2d} issue.

           @raise SciPyWarning: A C{scipy.interpolate.interp2d} warning
                                as exception.
        '''
        spi = self.scipy_interpolate

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
                               an invalid B{C{lli}}.

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
       C{kind='linear'}.
    '''
    _kind = _linear_
    _kmin =  2

    def __init__(self, knots, name=NN):
        '''New L{HeightLinear} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg name: Optional name for this height interpolator (C{str}).

           @raise HeightError: Insufficient number of B{C{knots}} or
                               an invalid B{C{knot}}.

           @raise ImportError: Package C{numpy} or C{scipy} not found
                               or not installed.

           @raise SciPyError: A C{scipy.interpolate.interp2d} issue.

           @raise SciPyWarning: A C{scipy.interpolate.interp2d} warning
                                as exception.
        '''
        HeightCubic.__init__(self, knots, name=name)

    if _FOR_DOCS:
        __call__ = HeightCubic.__call__
        height   = HeightCubic.height


class _HeightIDW(_HeightBase):
    '''(INTERNAL) Base class for U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW)
       height interpolators.
    '''
    _beta = 0     # inverse distance power
    _hs   = ()    # known heights
    _xs   = ()    # knot lons
    _ys   = ()    # knot lats

    def __init__(self, knots, beta=2, name=NN, **wrap_adjust):
        '''New L{_HeightIDW} interpolator.
        '''
        self._xs, self._ys, self._hs = _xyhs3(tuple, self._kmin, knots, off=False)
        self.beta = beta
        if name:
            self.name = name
        if wrap_adjust:
            _boolkwds(self, **wrap_adjust)

    def __call__(self, *llis):
        '''Interpolate the height for one or several locations.

           @arg llis: The location or locations (C{LatLon}, ... or
                      C{LatLon}s).

           @return: A single interpolated height (C{float}) or a list
                    or tuple of interpolated heights (C{float}s).

           @raise HeightError: Insufficient number of B{C{llis}}, an
                               invalid B{C{lli}} or L{pygeodesy.fidw}
                               issue.
        '''
        _as, xis, yis, _ = _axyllis4(tuple, llis, off=False)
        return _as(map(self._hIDW, xis, yis))

    def _datum_setter(self, datum, knots):
        '''(INTERNAL) Set the datum.

           @raise TypeError: Invalid B{C{datum}}.
        '''
        d = datum or getattr(knots[0], _datum_, datum)
        if d not in (None, self._datum):
            self._datum = _ellipsoidal_datum(d, name=self.name)

    def _distances(self, x, y):  # PYCHOK unused (x, y) radians
        '''(INTERNAL) I{Must be overloaded}, see function C{notOverloaded}.
        '''
        notOverloaded(self, x, y)

    def _distances_angular_(self, func_, x, y):
        # return angular distances from func_
        for xk, yk in zip(self._xs, self._ys):
            r, _ = unrollPI(xk, x, wrap=self._wrap)
            yield func_(yk, y, r)

    def _distances_angular_datum_(self, func_, x, y):
        # return angular distances from func_
        for xk, yk in zip(self._xs, self._ys):
            r, _ = unrollPI(xk, x, wrap=self._wrap)
            yield func_(yk, y, r, datum=self._datum)

    def _hIDW(self, x, y):
        # interpolate height at (x, y) radians or degrees
        try:
            ds = self._distances(x, y)
            return fidw(self._hs, ds, beta=self._beta)
        except (TypeError, ValueError) as x:
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
        self._beta = Int_(beta=beta, Error=HeightError, low=1, high=3)

    def height(self, lats, lons):
        '''Interpolate the height for one or several lat-/longitudes.

           @arg lats: Latitude or latitudes (C{degrees} or C{degrees}s).
           @arg lons: Longitude or longitudes (C{degrees} or C{degrees}s).

           @return: A single interpolated height (C{float}) or a list of
                    interpolated heights (C{float}s).

           @raise HeightError: Insufficient or non-matching number of
                               B{C{lats}} and B{C{lons}} or L{pygeodesy.fidw}
                               issue.
        '''
        return _HeightBase._height(self, lats, lons)


class HeightIDWcosineAndoyerLambert(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW)
       and the I{angular} distance in C{radians} from function
       L{pygeodesy.cosineAndoyerLambert_}.

       @see: L{HeightIDWcosineForsytheAndoyerLambert}, L{HeightIDWdistanceTo},
             L{HeightIDWflatLocal}, L{HeightIDWhubeny}, L{HeightIDWthomas},
             U{IDW<https://www.Geo.FU-Berlin.DE/en/v/soga/Geodata-analysis/
             geostatistics/Inverse-Distance-Weighting/index.html>} and
             U{SHEPARD_INTERP_2D<https://People.SC.FSU.edu/~jburkardt/c_src/
             shepard_interp_2d/shepard_interp_2d.html>}.
    '''
    _datum = _WGS84
    _wrap  =  False

    def __init__(self, knots, datum=None, beta=2, wrap=False, name=NN):
        '''New L{HeightIDWcosineAndoyerLambert} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg datum: Optional datum overriding the default C{Datums.WGS84}
                         and first B{C{knots}}' datum (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg beta: Inverse distance power (C{int} 1, 2, or 3).
           @kwarg wrap: Wrap and L{pygeodesy.unrollPI} longitudes (C{bool}).
           @kwarg name: Optional name for this height interpolator (C{str}).

           @raise HeightError: Insufficient number of B{C{knots}} or
                               an invalid B{C{knot}} or B{C{beta}}.

           @raise TypeError: Invalid B{C{datum}}.
        '''
        _HeightIDW.__init__(self, knots, beta=beta, name=name, wrap=wrap)
        self._datum_setter(datum, knots)

    def _distances(self, x, y):  # (x, y) radians
        return self._distances_angular_datum_(cosineAndoyerLambert_, x, y)

    if _FOR_DOCS:
        __call__ = _HeightIDW.__call__
        height   = _HeightIDW.height


class HeightIDWcosineForsytheAndoyerLambert(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW)
       and the I{angular} distance in C{radians} from function
       L{pygeodesy.cosineForsytheAndoyerLambert_}.

       @see: L{HeightIDWcosineAndoyerLambert}, L{HeightIDWdistanceTo},
             L{HeightIDWflatLocal}, L{HeightIDWhubeny}, L{HeightIDWthomas},
             U{IDW<https://www.Geo.FU-Berlin.DE/en/v/soga/Geodata-analysis/
             geostatistics/Inverse-Distance-Weighting/index.html>} and
             U{SHEPARD_INTERP_2D<https://People.SC.FSU.edu/~jburkardt/c_src/
             shepard_interp_2d/shepard_interp_2d.html>}.
    '''
    _datum = _WGS84
    _wrap  =  False

    def __init__(self, knots, datum=None, beta=2, wrap=False, name=NN):
        '''New L{HeightIDWcosineForsytheAndoyerLambert} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg datum: Optional datum overriding the default C{Datums.WGS84}
                         and first B{C{knots}}' datum (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg beta: Inverse distance power (C{int} 1, 2, or 3).
           @kwarg wrap: Wrap and L{pygeodesy.unrollPI} longitudes (C{bool}).
           @kwarg name: Optional name for this height interpolator (C{str}).

           @raise HeightError: Insufficient number of B{C{knots}} or
                               an invalid B{C{knot}} or B{C{beta}}.

           @raise TypeError: Invalid B{C{datum}}.
        '''
        _HeightIDW.__init__(self, knots, beta=beta, name=name, wrap=wrap)
        self._datum_setter(datum, knots)

    def _distances(self, x, y):  # (x, y) radians
        return self._distances_angular_datum_(cosineForsytheAndoyerLambert_, x, y)

    if _FOR_DOCS:
        __call__ = _HeightIDW.__call__
        height   = _HeightIDW.height


class HeightIDWcosineLaw(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW) and the
       I{angular} distance in C{radians} from function L{pygeodesy.cosineLaw_}.

       @see: L{HeightIDWequirectangular}, L{HeightIDWeuclidean},
             L{HeightIDWflatPolar}, L{HeightIDWhaversine}, L{HeightIDWvincentys},
             U{IDW<https://www.Geo.FU-Berlin.DE/en/v/soga/Geodata-analysis/
             geostatistics/Inverse-Distance-Weighting/index.html>} and
             U{SHEPARD_INTERP_2D<https://People.SC.FSU.edu/~jburkardt/c_src/
             shepard_interp_2d/shepard_interp_2d.html>}.

       @note: See note at function L{pygeodesy.vincentys_}.
    '''
    _wrap = False

    def __init__(self, knots, beta=2, wrap=False, name=NN):
        '''New L{HeightIDWcosineLaw} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg beta: Inverse distance power (C{int} 1, 2, or 3).
           @kwarg wrap: Wrap and L{pygeodesy.unrollPI} longitudes (C{bool}).
           @kwarg name: Optional name for this height interpolator (C{str}).

           @raise HeightError: Insufficient number of B{C{knots}} or
                               an invalid B{C{knot}} or B{C{beta}}.
        '''
        _HeightIDW.__init__(self, knots, beta=beta, name=name, wrap=wrap)

    def _distances(self, x, y):  # (x, y) radians
        return self._distances_angular_(cosineLaw_, x, y)

    if _FOR_DOCS:
        __call__ = _HeightIDW.__call__
        height   = _HeightIDW.height


class HeightIDWdistanceTo(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW)
       and the distance from the points' C{LatLon.distanceTo} method,
       conventionally in C{meter}.

       @see: L{HeightIDWcosineAndoyerLambert}, L{HeightIDWcosineForsytheAndoyerLambert},
             L{HeightIDWflatPolar}, L{HeightIDWkarney}, L{HeightIDWthomas},
             U{IDW<https://www.Geo.FU-Berlin.DE/en/v/soga/Geodata-analysis/
             geostatistics/Inverse-Distance-Weighting/index.html>} and
             U{SHEPARD_INTERP_2D<https://People.SC.FSU.edu/~jburkardt/c_src/
             shepard_interp_2d/shepard_interp_2d.html>}.
    '''
    _distanceTo_kwds = {}
    _ks              = ()

    def __init__(self, knots, beta=2, name=NN, **distanceTo_kwds):
        '''New L{HeightIDWdistanceTo} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg beta: Inverse distance power (C{int} 1, 2, or 3).
           @kwarg name: Optional name for this height interpolator (C{str}).
           @kwarg distanceTo_kwds: Optional keyword arguments for the
                                   B{C{points}}' C{LatLon.distanceTo}
                                   method.

           @raise HeightError: Insufficient number of B{C{knots}} or
                               an invalid B{C{knot}} or B{C{beta}}.

           @raise ImportError: Package U{geographiclib
                  <https://PyPI.org/project/geographiclib>} missing
                  iff B{C{points}} are L{ellipsoidalKarney.LatLon}s.

           @note: All B{C{points}} I{must} be instances of the same
                  ellipsoidal or spherical C{LatLon} class, I{not
                  checked however}.
        '''
        n, self._ks = len2(knots)
        if n < self._kmin:
            raise _insufficientError(self._kmin, knots=n)
        for i, k in enumerate(self._ks):
            if not callable(getattr(k, _distanceTo_, None)):
                i = Fmt.SQUARE(_knots_, i)
                raise HeightError(i, k, txt=_distanceTo_)

        # use knots[0] class and datum to create
        # compatible points in _HeightBase._height
        # instead of class LatLon_ and datum None
        self._datum = self._ks[0].datum
        self._LLis  = self._ks[0].classof

        self.beta = beta
        if name:
            self.name = name
        if distanceTo_kwds:
            self._distanceTo_kwds = distanceTo_kwds

    def __call__(self, *llis):
        '''Interpolate the height for one or several locations.

           @arg llis: The location or locations (C{LatLon}, ... or
                      C{LatLon}s).

           @return: A single interpolated height (C{float}) or a list
                    or tuple of interpolated heights (C{float}s).

           @raise HeightError: Insufficient number of B{C{llis}}, an
                               invalid B{C{lli}} or L{pygeodesy.fidw}
                               issue.
        '''
        _as, llis = _allis2(llis)
        return _as(map(self._hIDW, llis))

    if _FOR_DOCS:
        height = _HeightIDW.height

    def _hIDW(self, lli):  # PYCHOK expected
        # interpolate height at point lli
        try:
            kwds = self._distanceTo_kwds
            ds = (k.distanceTo(lli, **kwds) for k in self._ks)
            return fidw(self._hs, ds, beta=self._beta)
        except (TypeError, ValueError) as x:
            raise HeightError(str(x))

    @property  # NOT _RO, no caching
    def _hs(self):  # see HeightIDWkarney
        for k in self._ks:
            yield k.height


class HeightIDWequirectangular(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW) and
       the I{angular} distance in C{radians squared} like function
       L{pygeodesy.equirectangular_}.

       @see: L{HeightIDWeuclidean}, L{HeightIDWflatPolar},
             L{HeightIDWhaversine}, L{HeightIDWvincentys},
             U{IDW<https://www.Geo.FU-Berlin.DE/en/v/soga/Geodata-analysis/
             geostatistics/Inverse-Distance-Weighting/index.html>} and
             U{SHEPARD_INTERP_2D<https://People.SC.FSU.edu/~jburkardt/c_src/
             shepard_interp_2d/shepard_interp_2d.html>}.
    '''
    _adjust = True
    _wrap   = False

    def __init__(self, knots, adjust=True, wrap=False, name=NN):
        '''New L{HeightIDWequirectangular} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg adjust: Adjust the wrapped, unrolled longitudinal
                          delta by the cosine of the mean latitude (C{bool}).
           @kwarg wrap: Wrap and L{pygeodesy.unrollPI} longitudes (C{bool}).
           @kwarg name: Optional name for this height interpolator (C{str}).

           @raise HeightError: Insufficient number of B{C{knots}} or
                               an invalid B{C{knot}}.
        '''
        _HeightIDW.__init__(self, knots, beta=1, name=name, wrap=wrap,
                                                          adjust=adjust)

    def _distances(self, x, y):  # (x, y) radians**2
        for xk, yk in zip(self._xs, self._ys):
            d, _ = unrollPI(xk, x, wrap=self._wrap)
            if self._adjust:
                d *= _scale_rad(yk, y)
            yield hypot2(d, yk - y)  # like equirectangular_ distance2

    if _FOR_DOCS:
        __call__ = _HeightIDW.__call__
        height   = _HeightIDW.height


class HeightIDWeuclidean(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW) and the
       I{angular} distance in C{radians} from function L{pygeodesy.euclidean_}.

       @see: L{HeightIDWcosineLaw}, L{HeightIDWequirectangular},
             L{HeightIDWhaversine}, L{HeightIDWvincentys},
             U{IDW<https://www.Geo.FU-Berlin.DE/en/v/soga/Geodata-analysis/
             geostatistics/Inverse-Distance-Weighting/index.html>} and
             U{SHEPARD_INTERP_2D<https://People.SC.FSU.edu/~jburkardt/c_src/
             shepard_interp_2d/shepard_interp_2d.html>}.
    '''
    _adjust = True

    def __init__(self, knots, adjust=True, beta=2, name=NN):
        '''New L{HeightIDWeuclidean} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg adjust: Adjust the longitudinal delta by the cosine
                          of the mean latitude for B{C{adjust}}=C{True}.
           @kwarg beta: Inverse distance power (C{int} 1, 2, or 3).
           @kwarg name: Optional name for this height interpolator (C{str}).

           @raise HeightError: Insufficient number of B{C{knots}} or
                               an invalid B{C{knot}} or B{C{beta}}.
        '''
        _HeightIDW.__init__(self, knots, beta=beta, name=name, adjust=adjust)

    def _distances(self, x, y):  # (x, y) radians
        for xk, yk in zip(self._xs, self._ys):
            yield euclidean_(yk, y, xk - x, adjust=self._adjust)

    if _FOR_DOCS:
        __call__ = _HeightIDW.__call__
        height   = _HeightIDW.height


class HeightIDWflatLocal(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW)
       and the I{angular} distance in C{radians squared} like function
       L{pygeodesy.flatLocal_}/L{pygeodesy.hubeny_}.

       @see: L{HeightIDWcosineAndoyerLambert}, L{HeightIDWcosineForsytheAndoyerLambert},
             L{HeightIDWdistanceTo}, L{HeightIDWhubeny}, L{HeightIDWthomas},
             U{IDW<https://www.Geo.FU-Berlin.DE/en/v/soga/Geodata-analysis/
             geostatistics/Inverse-Distance-Weighting/index.html>} and
             U{SHEPARD_INTERP_2D<https://People.SC.FSU.edu/~jburkardt/c_src/
             shepard_interp_2d/shepard_interp_2d.html>}.
    '''
    _datum = _WGS84
    _wrap  =  False

    def __init__(self, knots, datum=None, beta=2, wrap=False, name=NN):
        '''New L{HeightIDWflatLocal}/L{HeightIDWhubeny} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg datum: Optional datum overriding the default C{Datums.WGS84}
                         and first B{C{knots}}' datum (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg beta: Inverse distance power (C{int} 1, 2, or 3).
           @kwarg wrap: Wrap and L{pygeodesy.unrollPI} longitudes (C{bool}).
           @kwarg name: Optional name for this height interpolator (C{str}).

           @raise HeightError: Insufficient number of B{C{knots}} or
                               an invalid B{C{knot}} or B{C{beta}}.

           @raise TypeError: Invalid B{C{datum}}.
        '''
        _HeightIDW.__init__(self, knots, beta=beta, name=name, wrap=wrap)
        self._datum_setter(datum, knots)

    def _distances(self, x, y):  # (x, y) radians
        _r2_ = self._datum.ellipsoid._hubeny2_
        return self._distances_angular_(_r2_, x, y)  # radians**2

    if _FOR_DOCS:
        __call__ = _HeightIDW.__call__
        height   = _HeightIDW.height


class HeightIDWflatPolar(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW)
       and the I{angular} distance in C{radians} from function
       L{pygeodesy.flatPolar_}.

       @see: L{HeightIDWcosineLaw}, L{HeightIDWequirectangular},
             L{HeightIDWeuclidean}, L{HeightIDWhaversine}, L{HeightIDWvincentys},
             U{IDW<https://www.Geo.FU-Berlin.DE/en/v/soga/Geodata-analysis/
             geostatistics/Inverse-Distance-Weighting/index.html>} and
             U{SHEPARD_INTERP_2D<https://People.SC.FSU.edu/~jburkardt/c_src/
             shepard_interp_2d/shepard_interp_2d.html>}.
    '''
    _wrap = False

    def __init__(self, knots, beta=2, wrap=False, name=NN):
        '''New L{HeightIDWflatPolar} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg beta: Inverse distance power (C{int} 1, 2, or 3).
           @kwarg wrap: Wrap and L{pygeodesy.unrollPI} longitudes (C{bool}).
           @kwarg name: Optional name for this height interpolator (C{str}).

           @raise HeightError: Insufficient number of B{C{knots}} or
                               an invalid B{C{knot}} or B{C{beta}}.
        '''
        _HeightIDW.__init__(self, knots, beta=beta, name=name, wrap=wrap)

    def _distances(self, x, y):  # (x, y) radians
        return self._distances_angular_(flatPolar_, x, y)

    if _FOR_DOCS:
        __call__ = _HeightIDW.__call__
        height   = _HeightIDW.height


class HeightIDWhaversine(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW) and the
       I{angular} distance in C{radians} from function L{pygeodesy.haversine_}.

       @see: L{HeightIDWcosineLaw}, L{HeightIDWequirectangular}, L{HeightIDWeuclidean},
             L{HeightIDWflatPolar}, L{HeightIDWvincentys},
             U{IDW<https://www.Geo.FU-Berlin.DE/en/v/soga/Geodata-analysis/
             geostatistics/Inverse-Distance-Weighting/index.html>} and
             U{SHEPARD_INTERP_2D<https://People.SC.FSU.edu/~jburkardt/c_src/
             shepard_interp_2d/shepard_interp_2d.html>}.

       @note: See note at function L{pygeodesy.vincentys_}.
    '''
    _wrap = False

    def __init__(self, knots, beta=2, wrap=False, name=NN):
        '''New L{HeightIDWhaversine} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg beta: Inverse distance power (C{int} 1, 2, or 3).
           @kwarg wrap: Wrap and L{pygeodesy.unrollPI} longitudes (C{bool}).
           @kwarg name: Optional name for this height interpolator (C{str}).

           @raise HeightError: Insufficient number of B{C{knots}} or
                               an  B{C{knot}} or B{C{beta}}.
        '''
        _HeightIDW.__init__(self, knots, beta=beta, name=name, wrap=wrap)

    def _distances(self, x, y):  # (x, y) radians
        return self._distances_angular_(haversine_, x, y)

    if _FOR_DOCS:
        __call__ = _HeightIDW.__call__
        height   = _HeightIDW.height


class HeightIDWhubeny(HeightIDWflatLocal):  # for Karl Hubeny
    if _FOR_DOCS:
        __doc__  = HeightIDWflatLocal.__doc__
        __init__ = HeightIDWflatLocal.__init__
        __call__ = HeightIDWflatLocal.__call__
        height   = HeightIDWflatLocal.height


class HeightIDWkarney(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW) and
       the I{angular} distance in C{degrees} from I{Karney}'s
       U{geographiclib<https://PyPI.org/project/geographiclib>} U{Geodesic
       <https://GeographicLib.SourceForge.io/html/python/code.html>}
       Inverse method.

       @see: L{HeightIDWcosineAndoyerLambert},
             L{HeightIDWcosineForsytheAndoyerLambert}, L{HeightIDWdistanceTo},
             L{HeightIDWflatLocal}, L{HeightIDWhubeny}, L{HeightIDWthomas},
             U{IDW<https://www.Geo.FU-Berlin.DE/en/v/soga/Geodata-analysis/
             geostatistics/Inverse-Distance-Weighting/index.html>} and
             U{SHEPARD_INTERP_2D<https://People.SC.FSU.edu/~jburkardt/c_src/
             shepard_interp_2d/shepard_interp_2d.html>}.
    '''
    _datum    = _WGS84
    _Inverse1 =  None
    _ks       = ()
    _wrap     =  False

    def __init__(self, knots, datum=None, beta=2, wrap=False, name=NN):
        '''New L{HeightIDWkarney} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg datum: Optional datum overriding the default C{Datums.WGS84}
                         and first B{C{knots}}' datum (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg beta: Inverse distance power (C{int} 1, 2, or 3).
           @kwarg wrap: Wrap and L{pygeodesy.unroll180} longitudes (C{bool}).
           @kwarg name: Optional name for this height interpolator (C{str}).

           @raise HeightError: Insufficient number of B{C{knots}} or
                               an invalid B{C{knot}}, B{C{datum}} or
                               B{C{beta}}.

           @raise ImportError: Package U{geographiclib
                  <https://PyPI.org/project/geographiclib>} missing.

           @raise TypeError: Invalid B{C{datum}}.
        '''
        n, self._ks = len2(knots)
        if n < self._kmin:
            raise _insufficientError(self._kmin, knots=n)
        self._datum_setter(datum, self._ks)
        self._Inverse1 = self.datum.ellipsoid.geodesic.Inverse1

        self.beta = beta
        if wrap:
            self._wrap = True
        if name:
            self.name = name

    def _distances(self, x, y):  # (x, y) degrees
        for k in self._ks:
            # non-negative I{angular} distance in C{degrees}
            yield self._Inverse1(y, x, k.lat, k.lon, wrap=self._wrap)

    @property  # NOT _RO, no caching
    def _hs(self):  # see HeightIDWdistanceTo
        for k in self._ks:
            yield k.height

    def __call__(self, *llis):
        '''Interpolate the height for one or several locations.

           @arg llis: The location or locations (C{LatLon}, ... or
                      C{LatLon}s).

           @return: A single interpolated height (C{float}) or a list
                    or tuple of interpolated heights (C{float}s).

           @raise HeightError: Insufficient number of B{C{llis}}, an
                               invalid B{C{lli}} or L{pygeodesy.fidw}
                               issue.
        '''
        def _xy2(lls):
            try:  # like _xyhs above, but keeping degrees
                for i, ll in enumerate(lls):
                    yield ll.lon, ll.lat
            except AttributeError as x:
                i = Fmt.SQUARE(llis=i)
                raise HeightError(i, ll, txt=str(x))

        _as, llis = _allis2(llis)
        return _as(map(self._hIDW, *zip(*_xy2(llis))))

    if _FOR_DOCS:
        height = _HeightIDW.height


class HeightIDWthomas(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW) and the
       I{angular} distance in C{radians} from function L{pygeodesy.thomas_}.

       @see: L{HeightIDWcosineAndoyerLambert}, L{HeightIDWcosineForsytheAndoyerLambert},
             L{HeightIDWdistanceTo}, L{HeightIDWflatLocal}, L{HeightIDWhubeny},
             U{IDW<https://www.Geo.FU-Berlin.DE/en/v/soga/Geodata-analysis/
             geostatistics/Inverse-Distance-Weighting/index.html>} and
             U{SHEPARD_INTERP_2D<https://People.SC.FSU.edu/~jburkardt/c_src/
             shepard_interp_2d/shepard_interp_2d.html>}.
    '''
    _datum = _WGS84
    _wrap  =  False

    def __init__(self, knots, datum=None, beta=2, wrap=False, name=NN):
        '''New L{HeightIDWthomas} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg datum: Optional datum overriding the default C{Datums.WGS84}
                         and first B{C{knots}}' datum (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg beta: Inverse distance power (C{int} 1, 2, or 3).
           @kwarg wrap: Wrap and L{pygeodesy.unrollPI} longitudes (C{bool}).
           @kwarg name: Optional name for this height interpolator (C{str}).

           @raise HeightError: Insufficient number of B{C{knots}} or
                               an invalid B{C{knot}} or B{C{beta}}.

           @raise TypeError: Invalid B{C{datum}}.
        '''
        _HeightIDW.__init__(self, knots, beta=beta, name=name, wrap=wrap)
        self._datum_setter(datum, knots)

    def _distances(self, x, y):  # (x, y) radians
        return self._distances_angular_datum_(thomas_, x, y)

    if _FOR_DOCS:
        __call__ = _HeightIDW.__call__
        height   = _HeightIDW.height


class HeightIDWvincentys(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW) and the
       I{angular} distance in C{radians} from function L{pygeodesy.vincentys_}.

       @see: L{HeightIDWcosineLaw}, L{HeightIDWequirectangular},
             L{HeightIDWeuclidean}, L{HeightIDWflatPolar}, L{HeightIDWhaversine},
             U{IDW<https://www.Geo.FU-Berlin.DE/en/v/soga/Geodata-analysis/
             geostatistics/Inverse-Distance-Weighting/index.html>} and
             U{SHEPARD_INTERP_2D<https://People.SC.FSU.edu/~jburkardt/c_src/
             shepard_interp_2d/shepard_interp_2d.html>}.

       @note: See note at function L{pygeodesy.vincentys_}.
    '''
    _wrap = False

    def __init__(self, knots, beta=2, wrap=False, name=NN):
        '''New L{HeightIDWvincentys} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg beta: Inverse distance power (C{int} 1, 2, or 3).
           @kwarg wrap: Wrap and L{pygeodesy.unrollPI} longitudes (C{bool}).
           @kwarg name: Optional name for this height interpolator (C{str}).

           @raise HeightError: Insufficient number of B{C{knots}} or
                               an invalid B{C{knot}} or B{C{beta}}.
        '''
        _HeightIDW.__init__(self, knots, beta=beta, name=name, wrap=wrap)

    def _distances(self, x, y):  # (x, y) radians
        return self._distances_angular_(vincentys_, x, y)

    if _FOR_DOCS:
        __call__ = _HeightIDW.__call__
        height   = _HeightIDW.height


class HeightLSQBiSpline(_HeightBase):
    '''Height interpolator using C{SciPy} U{LSQSphereBivariateSpline
       <https://docs.SciPy.org/doc/scipy/reference/generated/scipy.
       interpolate.LSQSphereBivariateSpline.html>}.
    '''
    _kmin = 16  # k = 3, always

    def __init__(self, knots, weight=None, name=NN):
        '''New L{HeightLSQBiSpline} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg weight: Optional weight or weights for each B{C{knot}}
                          (C{scalar} or C{scalar}s).
           @kwarg name: Optional name for this height interpolator (C{str}).

           @raise HeightError: Insufficient number of B{C{knots}} or
                               an invalid B{C{knot}} or B{C{weight}}.

           @raise LenError: Unequal number of B{C{knots}} and B{C{weight}}s.

           @raise ImportError: Package C{numpy} or C{scipy} not found
                               or not installed.

           @raise SciPyError: A C{scipy} C{LSQSphereBivariateSpline} issue.

           @raise SciPyWarning: A C{scipy} C{LSQSphereBivariateSpline}
                                warning as exception.
        '''
        np  = self.numpy
        spi = self.scipy_interpolate

        xs, ys, hs = self._xyhs3(knots)
        n = len(hs)

        w = weight
        if isscalar(w):
            w = float(w)
            if w <= 0:
                raise HeightError(weight=w)
            w = [w] * n
        elif w is not None:
            m, w = len2(w)
            if m != n:
                raise LenError(HeightLSQBiSpline, weight=m, knots=n)
            w = map2(float, w)
            m = min(w)
            if m <= 0:
                w = Fmt.SQUARE(weight=w.find(m))
                raise HeightError(w, m)
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
                               an invalid B{C{lli}}.

           @raise SciPyError: A C{scipy} C{LSQSphereBivariateSpline} issue.

           @raise SciPyWarning: A C{scipy} C{LSQSphereBivariateSpline}
                                warning as exception.
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

           @raise SciPyError: A C{scipy} C{LSQSphereBivariateSpline} issue.

           @raise SciPyWarning: A C{scipy} C{LSQSphereBivariateSpline}
                                warning as exception.
        '''
        return _HeightBase._height(self, lats, lons)


class HeightSmoothBiSpline(_HeightBase):
    '''Height interpolator using C{SciPy} U{SmoothSphereBivariateSpline
       <https://docs.SciPy.org/doc/scipy/reference/generated/scipy.
       interpolate.SmoothSphereBivariateSpline.html>}.
    '''
    _kmin = 16  # k = 3, always

    def __init__(self, knots, s=4, name=NN):
        '''New L{HeightSmoothBiSpline} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg s: The spline smoothing factor (C{scalar}), default C{4}.
           @kwarg name: Optional name for this height interpolator (C{str}).

           @raise HeightError: Insufficient number of B{C{knots}} or
                               an invalid B{C{knot}} or B{C{s}}.

           @raise ImportError: Package C{numpy} or C{scipy} not found
                               or not installed.

           @raise SciPyError: A C{scipy} C{SmoothSphereBivariateSpline} issue.

           @raise SciPyWarning: A C{scipy} C{SmoothSphereBivariateSpline}
                                warning as exception.
        '''
        spi = self.scipy_interpolate

        s = Float_(s, name='smoothing', Error=HeightError, low=4)

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
                               an invalid B{C{lli}}.

           @raise SciPyError: A C{scipy} C{SmoothSphereBivariateSpline} issue.

           @raise SciPyWarning: A C{scipy} C{SmoothSphereBivariateSpline}
                                warning as exception.
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

           @raise SciPyError: A C{scipy} C{SmoothSphereBivariateSpline} issue.

           @raise SciPyWarning: A C{scipy} C{SmoothSphereBivariateSpline}
                                warning as exception.
        '''
        return _HeightBase._height(self, lats, lons)


__all__ += _ALL_DOCS(_HeightBase)

# **) MIT License
#
# Copyright (C) 2016-2022 -- mrJean1 at Gmail -- All Rights Reserved.
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
