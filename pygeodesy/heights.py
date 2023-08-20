
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

C{>>> ll = LatLon(1, 2, ...)}
C{>>> h = hinterpolator(ll)}

or

C{>>> h0, h1, h2, ... = hinterpolator(ll0, ll1, ll2, ...)}

or a list, tuple, generator, etc. of C{LatLon}s

C{>>> hs = hinterpolator(lls)}

5. For separate lat- and longitudes invoke the C{.height} method

C{>>> h = hinterpolator.height(lat, lon)}

or as 2 lists, 2 tuples, etc.

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
from pygeodesy.constants import EPS, PI, PI2, _0_0, _90_0, _180_0
from pygeodesy.datums import _ellipsoidal_datum, _WGS84
from pygeodesy.errors import _AssertionError, LenError, PointsError, \
                             _SciPyIssue, _xattr, _xkwds, _xkwds_get
from pygeodesy.fmath import fidw,  Float_, Int_
from pygeodesy.formy import cosineAndoyerLambert, cosineForsytheAndoyerLambert, \
                            cosineLaw, equirectangular_, euclidean, flatLocal, \
                            flatPolar, haversine, thomas, vincentys,  radians
from pygeodesy.interns import NN, _COMMASPACE_, _cubic_, _insufficient_, _knots_, \
                             _linear_, _NOTEQUAL_, _PLUS_, _scipy_, _SPACE_, _STAR_
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _FOR_DOCS
from pygeodesy.named import _Named, notOverloaded
from pygeodesy.points import _distanceTo, LatLon_,  Fmt, _Wrap
from pygeodesy.props import Property_RO, property_RO
# from pygeodesy.streprs import Fmt  # from .points
# from pygeodesy.units import Float_, Int_  # from .fmath
# from pygeodesy.utily import _Wrap  # from .points

# from math import radians  # from .formy

__all__ = _ALL_LAZY.heights
__version__ = '23.08.15'

_error_     = 'error'
_llis_      = 'llis'
_smoothing_ = 'smoothing'


class HeightError(PointsError):
    '''Height interpolator C{Height...} or interpolation issue.
    '''
    pass


def _alist(ais):
    # return list of floats, not numpy.float64s
    return list(map(float, ais))


def _ascalar(ais):  # in .geoids
    # return single float, not numpy.float64
    ais = list(ais)  # np.array, etc. to list
    if len(ais) != 1:
        n =  Fmt.PAREN(len=repr(ais))
        t = _SPACE_(len(ais), _NOTEQUAL_, 1)
        raise _AssertionError(n, txt=t)
    return float(ais[0])  # remove np.<type>


def _atuple(ais):
    # return tuple of floats, not numpy.float64s
    return tuple(map(float, ais))


def _as_llis2(llis, m=1, Error=HeightError):  # in .geoids
    # determine return type and convert lli C{LatLon}s to list
    if not isinstance(llis, tuple):  # llis are *args
        n = Fmt.PAREN(type_=_STAR_(NN, _llis_))
        raise _AssertionError(n, txt=repr(llis))

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


def _height_called(inst, lats, lons, Error=HeightError, **wrap):  # in .geoids
    LLis, d = inst._LLis, inst.datum
    if isscalar(lats) and isscalar(lons):
        llis = LLis(lats, lons, datum=d)
    else:
        n, lats = len2(lats)
        m, lons = len2(lons)
        if n != m:
            # format a LenError, but raise an Error
            e = LenError(inst.__class__, lats=n, lons=m, txt=None)
            raise e if Error is LenError else Error(str(e))
        llis = [LLis(*t, datum=d) for t in zip(lats, lons)]
    return inst(llis, **wrap)  # __call__(lli) or __call__(llis)


def _insufficientError(need, Error=HeightError, **name_value):  # PYCHOK no cover
    # create an insufficient Error instance
    t = _COMMASPACE_(_insufficient_, str(need) + _PLUS_)
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


def _xyhs(lls, wrap=False, name=_llis_):
    # map (lat, lon, h) to (x, y, h) in radians, offset
    # x as 0 <= lon <= PI2 and y as 0 <= lat <= PI
    _0, _90, _180 = _0_0, _90_0,    _180_0
    _m, _r,  _w   =  max,  radians, _Wrap._latlonop(wrap)
    try:
        for i, ll in enumerate(lls):
            y, x = _w(ll.lat, ll.lon)
            yield _m(_0, _r(x + _180)), \
                  _m(_0, _r(y +  _90)), ll.height
    except AttributeError as x:
        i = Fmt.SQUARE(name, i)
        raise HeightError(i, ll, cause=x)


class _HeightBase(_Named):  # in .geoids
    '''(INTERNAL) Interpolator base class.
    '''
    _datum = _WGS84    # default
    _kmin  =  2        # min number of knots
    _LLis  =  LatLon_  # ._height class
    _np_sp =  None     # (numpy, scipy)
    _wrap  =  None     # wrap knots and llis

    @property_RO
    def datum(self):
        '''Get the C{datum} setting or the default (L{Datum}).
        '''
        return self._datum

    @property_RO
    def kmin(self):
        '''Get the minimum number of knots (C{int}).
        '''
        return self._kmin

    @property_RO
    def wrap(self):
        '''Get the C{wrap} setting (C{bool}) or C{None}.
        '''
        return self._wrap


class _HeightsBase(_HeightBase):  # in .geoids
    '''(INTERNAL) Interpolator base class.
    '''
    _np_sp = None  # (numpy, scipy)

    def __call__(self, *llis, **wrap):  # PYCHOK no cover
        '''Interpolate the height for one or several locations.

           @arg llis: One or more locations (C{LatLon}s), all positional.
           @kwarg wrap: If C{True}, wrap or I{normalize} all B{C{llis}}
                        locations (C{bool}), overriding the B{C{knots}}
                        setting.

           @return: A single interpolated height (C{float}) or a list
                    or tuple of interpolated heights (C{float}s).

           @raise HeightError: Insufficient number of B{C{llis}} or
                               an invalid B{C{lli}}.

           @raise SciPyError: A C{scipy} issue.

           @raise SciPyWarning: A C{scipy} warning as exception.
        '''
        notOverloaded(self, callername='__call__', *llis, **wrap)

    def _as_xyllis4(self, llis, **wrap):
        # convert lli C{LatLon}s to tuples or C{NumPy} arrays of
        # C{SciPy} sphericals and determine the return type
        atype = self.numpy.array
        wrap = _xkwds(wrap, wrap=self._wrap)
        _as, llis   = _as_llis2(llis)
        xis, yis, _ =  zip(*_xyhs(llis, **wrap))  # PYCHOK yield
        return _as, atype(xis), atype(yis), llis

    def _ev(self, *args):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}, see function C{notOverloaded}.
        '''
        notOverloaded(self, *args)

    def _eval(self, llis, **wrap):  # XXX single arg, not *args
        _as, xis, yis, _ = self._as_xyllis4(llis, **wrap)
        try:  # SciPy .ev signature: y first, then x!
            return _as(self._ev(yis, xis))
        except Exception as x:
            raise _SciPyIssue(x)

    def height(self, lats, lons, **wrap):  # PYCHOK no cover
        '''Interpolate the height for one or several lat-/longitudes.

           @arg lats: Latitude or latitudes (C{degrees} or C{degrees}s).
           @arg lons: Longitude or longitudes (C{degrees} or C{degrees}s).
           @kwarg wrap: If C{True}, wrap or I{normalize} all B{C{lats}}
                        and B{C{lons}} locationts (C{bool}), overriding
                        the B{C{knots}}  setting.

           @return: A single interpolated height (C{float}) or a list of
                    interpolated heights (C{float}s).

           @raise HeightError: Insufficient or non-matching number of
                               B{C{lats}} and B{C{lons}}.

           @raise SciPyError: A C{scipy} issue.

           @raise SciPyWarning: A C{scipy} warning as exception.
        '''
        notOverloaded(self, lats, lons, **wrap)

    def _np_sp2(self, throwarnings=False):
        '''(INTERNAL) Import C{numpy} and C{scipy}, once.
        '''
        t = _HeightsBase._np_sp
        if not t:
            # raise SciPyWarnings, but not if
            # scipy has been imported already
            if throwarnings:  # PYCHOK no cover
                import sys
                if _scipy_ not in sys.modules:
                    import warnings
                    warnings.filterwarnings(_error_)

            sp = _xscipy(self.__class__, 1, 2)
            np = _xnumpy(self.__class__, 1, 9)

            _HeightsBase._np_sp = t = np, sp
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

    def _xyhs3(self, knots, **wrap):
        # convert knot C{LatLon}s to tuples or C{NumPy} arrays and C{SciPy} sphericals
        xs, ys, hs = zip(*_xyhs(knots, name=_knots_, **wrap))  # PYCHOK yield
        n = len(hs)
        if n < self.kmin:
            raise _insufficientError(self.kmin, knots=n)
        return map1(self.numpy.array, xs, ys, hs)


class HeightCubic(_HeightsBase):
    '''Height interpolator based on C{SciPy} U{interp2d<https://docs.SciPy.org/
       doc/scipy/reference/generated/scipy.interpolate.interp2d.html>}
       C{kind='cubic'}.
    '''
    _interp2d =  None
    _kind     = _cubic_
    _kmin     =  16

    def __init__(self, knots, name=NN, **wrap):
        '''New L{HeightCubic} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg name: Optional name for this height interpolator (C{str}).
           @kwarg wrap: If C{True}, wrap or I{normalize} all B{C{knots}}
                        and B{C{llis}} locations (C{bool}).

           @raise HeightError: Insufficient number of B{C{knots}} or
                               invalid B{C{knot}}.

           @raise ImportError: Package C{numpy} or C{scipy} not found
                               or not installed.

           @raise SciPyError: A C{scipy.interpolate.interp2d} issue.

           @raise SciPyWarning: A C{scipy.interpolate.interp2d} warning
                                as exception.
        '''
        spi = self.scipy_interpolate

        xs, ys, hs = self._xyhs3(knots, **wrap)
        try:  # SciPy.interpolate.interp2d kind 'linear' or 'cubic'
            self._interp2d = spi.interp2d(xs, ys, hs, kind=self._kind)
        except Exception as x:
            raise _SciPyIssue(x)

        if name:
            self.name = name

    def __call__(self, *llis, **wrap):
        '''Interpolate the height for one or several locations.

           @raise SciPyError: A C{scipy.interpolate.interp2d} issue.

           @raise SciPyWarning: A C{scipy.interpolate.interp2d} warning
                                as exception.

           @see: L{_HeightsBase.__call__} for details about B{C{llis}},
                  B{C{wrap}}, return values and other exceptions.
        '''
        return _HeightsBase._eval(self, llis, **wrap)

    def _ev(self, yis, xis):  # PYCHOK expected
        # to make SciPy .interp2d signature(x, y), single (x, y)
        # match SciPy .ev signature(ys, xs), flipped multiples
        return map(self._interp2d, xis, yis)

    def height(self, lats, lons, **wrap):
        '''Interpolate the height for one or several lat-/longitudes.

           @raise SciPyError: A C{scipy.interpolate.interp2d} issue.

           @raise SciPyWarning: A C{scipy.interpolate.interp2d} warning
                                as exception.

           @see: L{_HeightsBase.height} for details about B{C{lats}},
                 B{C{lons}}, B{C{wrap}}, return values and other exceptions.
        '''
        return _height_called(self, lats, lons, **wrap)


class HeightLinear(HeightCubic):
    '''Height interpolator based on C{SciPy} U{interp2d<https://docs.SciPy.org/
       doc/scipy/reference/generated/scipy.interpolate.interp2d.html>}
       C{kind='linear'}.
    '''
    _kind = _linear_
    _kmin =  2

    def __init__(self, knots, name=NN, **wrap):
        '''New L{HeightLinear} interpolator.

           @see: L{HeightCubic.__init__} for details about B{C{knots}},
                  B{C{name}}, B{C{wrap}} and other exceptions.
        '''
        HeightCubic.__init__(self, knots, name=name, **wrap)

    if _FOR_DOCS:
        __call__ = HeightCubic.__call__
        height   = HeightCubic.height


class HeightLSQBiSpline(_HeightsBase):
    '''Height interpolator using C{SciPy} U{LSQSphereBivariateSpline
       <https://docs.SciPy.org/doc/scipy/reference/generated/scipy.
       interpolate.LSQSphereBivariateSpline.html>}.
    '''
    _kmin = 16  # k = 3, always

    def __init__(self, knots, weight=None, name=NN, **wrap):
        '''New L{HeightLSQBiSpline} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg weight: Optional weight or weights for each B{C{knot}}
                          (C{scalar} or C{scalar}s).
           @kwarg name: Optional name for this height interpolator (C{str}).
           @kwarg wrap: If C{True}, wrap or I{normalize} all B{C{knots}}
                        and B{C{llis}} locations (C{bool}).

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

        xs, ys, hs = self._xyhs3(knots, **wrap)
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
            if m <= 0:  # PYCHOK no cover
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

    def __call__(self, *llis, **wrap):
        '''Interpolate the height for one or several locations.

           @raise SciPyError: A C{scipy} C{LSQSphereBivariateSpline} issue.

           @raise SciPyWarning: A C{scipy} C{LSQSphereBivariateSpline}
                                warning as exception.

           @see: L{_HeightsBase.__call__} for details about B{C{llis}},
                  B{C{wrap}}, return values and other exceptions.
        '''
        return _HeightsBase._eval(self, llis, **wrap)

    def height(self, lats, lons, **wrap):
        '''Interpolate the height for one or several lat-/longitudes.

           @raise SciPyError: A C{scipy} C{LSQSphereBivariateSpline} issue.

           @raise SciPyWarning: A C{scipy} C{LSQSphereBivariateSpline}
                                warning as exception.

           @see: L{_HeightsBase.height} for details about B{C{lats}},
                 B{C{lons}}, B{C{wrap}}, return values and other exceptions.
        '''
        return _height_called(self, lats, lons, **wrap)


class HeightSmoothBiSpline(_HeightsBase):
    '''Height interpolator using C{SciPy} U{SmoothSphereBivariateSpline
       <https://docs.SciPy.org/doc/scipy/reference/generated/scipy.
       interpolate.SmoothSphereBivariateSpline.html>}.
    '''
    _kmin = 16  # k = 3, always

    def __init__(self, knots, s=4, name=NN, **wrap):
        '''New L{HeightSmoothBiSpline} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg s: The spline smoothing factor (C{scalar}), default C{4}.
           @kwarg name: Optional name for this height interpolator (C{str}).
           @kwarg wrap: If C{True}, wrap or I{normalize} the B{C{knots}}
                        and any called B{C{llis}} and height B{C{lats}}
                        and B{C{lons}} locations (C{bool}).

           @raise HeightError: Insufficient number of B{C{knots}} or
                               an invalid B{C{knot}} or B{C{s}}.

           @raise ImportError: Package C{numpy} or C{scipy} not found
                               or not installed.

           @raise SciPyError: A C{scipy} C{SmoothSphereBivariateSpline} issue.

           @raise SciPyWarning: A C{scipy} C{SmoothSphereBivariateSpline}
                                warning as exception.
        '''
        spi = self.scipy_interpolate

        s = Float_(s, name=_smoothing_, Error=HeightError, low=4)

        xs, ys, hs = self._xyhs3(knots, **wrap)
        try:
            self._ev = spi.SmoothSphereBivariateSpline(ys, xs, hs,
                                                       eps=EPS, s=s).ev
        except Exception as x:
            raise _SciPyIssue(x)

        if name:
            self.name = name

    def __call__(self, *llis, **wrap):
        '''Interpolate the height for one or several locations.

           @raise SciPyError: A C{scipy} C{SmoothSphereBivariateSpline} issue.

           @raise SciPyWarning: A C{scipy} C{SmoothSphereBivariateSpline}
                                warning as exception.

           @see: L{_HeightsBase.__call__} for details about B{C{llis}},
                  B{C{wrap}}, return values and other exceptions.
        '''
        return _HeightsBase._eval(self, llis, **wrap)

    def height(self, lats, lons, **wrap):
        '''Interpolate the height for one or several lat-/longitudes.

           @raise SciPyError: A C{scipy} C{SmoothSphereBivariateSpline} issue.

           @raise SciPyWarning: A C{scipy} C{SmoothSphereBivariateSpline}
                                warning as exception.

           @see: L{_HeightsBase.height} for details about B{C{lats}},
                 B{C{lons}}, B{C{wrap}}, return values and other exceptions.
        '''
        return _height_called(self, lats, lons, **wrap)


class _HeightIDW(_HeightBase):
    '''(INTERNAL) Base class for U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW) height
       interpolators.

       @see: U{IDW<https://www.Geo.FU-Berlin.DE/en/v/soga/Geodata-analysis/
             geostatistics/Inverse-Distance-Weighting/index.html>},
             U{SHEPARD_INTERP_2D<https://People.SC.FSU.edu/~jburkardt/c_src/
             shepard_interp_2d/shepard_interp_2d.html>} and other C{_HeightIDW*}
             classes.
    '''
    _beta  =  0     # fidw inverse power
    _func  =  None  # formy function
    _ks    = ()     # knots list or tuple
    _kwds  = {}     # func_ options

    def __init__(self, knots, beta=2, name=NN, **kwds):
        '''New C{_HeightIDW*} interpolator.

           @arg knots: The points with known height (C{LatLon}s).
           @kwarg beta: Inverse distance power (C{int} 1, 2, or 3).
           @kwarg name: Optional name for this height interpolator (C{str}).
           @kwarg kwds: Optional keyword argument for distance function,
                        retrievable with property C{kwds}.

           @raise HeightError: Insufficient number of B{C{knots}} or
                               an invalid B{C{knot}} or B{C{beta}}.
        '''
        if name:
            self.name = name
        n, self._ks = len2(knots)
        if n < self.kmin:
            raise _insufficientError(self.kmin, knots=n)
        self.beta = beta
        self._kwds = kwds or {}

    def __call__(self, *llis, **wrap):
        '''Interpolate the height for one or several locations.

           @arg llis: One or more locations (C{LatLon}s), all positional.
           @kwarg wrap: If C{True}, wrap or I{normalize} all B{C{llis}}
                        locations (C{bool}).

           @return: A single interpolated height (C{float}) or a list
                    or tuple of interpolated heights (C{float}s).

           @raise HeightError: Insufficient number of B{C{llis}}, an
                               invalid B{C{lli}} or L{pygeodesy.fidw}
                               issue.
        '''
        def _xy2(lls, wrap=False):
            _w = _Wrap._latlonop(wrap)
            try:  # like _xyhs above, but degrees
                for i, ll in enumerate(lls):
                    yield _w(ll.lon, ll.lat)
            except AttributeError as x:
                i = Fmt.SQUARE(llis=i)
                raise HeightError(i, ll, cause=x)

        _as, llis = _as_llis2(llis)
        return _as(map(self._hIDW, *zip(*_xy2(llis, **wrap))))

    @property_RO
    def adjust(self):
        '''Get the C{adjust} setting (C{bool}) or C{None}.
        '''
        return _xkwds_get(self._kwds, adjust=None)

    @property
    def beta(self):
        '''Get the inverse distance power (C{int}).
        '''
        return self._beta

    @beta.setter  # PYCHOK setter!
    def beta(self, beta):
        '''Set the inverse distance power (C{int} 1, 2, or 3).

           @raise HeightError: Invalid B{C{beta}}.
        '''
        self._beta = Int_(beta=beta, Error=HeightError, low=1, high=3)

    @property_RO
    def datum(self):
        '''Get the C{datum} setting or the default (L{Datum}).
        '''
        return _xkwds_get(self._kwds, datum=self._datum)

    def _datum_setter(self, datum):
        '''(INTERNAL) Set the default C{datum}.
        '''
        d = datum or _xattr(self._ks[0], datum=None)
        if d and d is not self._datum:
            self._datum = _ellipsoidal_datum(d, name=self.name)

    def _distances(self, x, y):
        '''(INTERNAL) Yield distances to C{(x, y)}.
        '''
        _f, kwds = self._func, self._kwds
        if not callable(_f):  # PYCHOK no cover
            notOverloaded(self, distance_function=_f)
        for _, k in enumerate(self._ks):
            yield _f(y, x, k.lat, k.lon, **kwds)

    def _hIDW(self, x, y):
        '''(INTERNAL) Return the IDW-interpolate height at
           location (x, y), both C{degrees} or C{radians}.
        '''
        try:
            ds =  self._distances(x, y)
            hs = (k.height for k in self._ks)
            return fidw(hs, ds, beta=self._beta)
        except (TypeError, ValueError) as x:
            raise HeightError(str(x), cause=x)

    def height(self, lats, lons, **wrap):
        '''Interpolate the height for one or several lat-/longitudes.

           @arg lats: Latitude or latitudes (C{degrees} or C{degrees}s).
           @arg lons: Longitude or longitudes (C{degrees} or C{degrees}s).
           @kwarg wrap: If C{True}, wrap or I{normalize} all B{C{lats}}
                        and B{C{lons}} (C{bool}).

           @return: A single interpolated height (C{float}) or a list of
                    interpolated heights (C{float}s).

           @raise HeightError: Insufficient or non-matching number of
                               B{C{lats}} and B{C{lons}} or L{pygeodesy.fidw}
                               issue.
        '''
        return _height_called(self, lats, lons, **wrap)

    @property_RO
    def hypot(self):
        '''Get the C{hypot} setting (C{callable}) or C{None}.
        '''
        return _xkwds_get(self._kwds, hypot=None)

    @property_RO
    def knots(self):
        '''Get the B{C{knots}} (C{list} or C{tuple}).
        '''
        return self._ks

    @property_RO
    def kwds(self):
        '''Get the optional keyword arguments (C{dict}).
        '''
        return self._kwds

    @property_RO
    def limit(self):
        '''Get the C{limit} setting (C{degrees}) or C{None}.
        '''
        return _xkwds_get(self._kwds, limit=None)

    @property_RO
    def radius(self):
        '''Get the C{radius} setting (C{bool}) or C{None}.
        '''
        return _xkwds_get(self._kwds, radius=None)

    @property_RO
    def scaled(self):
        '''Get the C{scaled} setting (C{bool}) or C{None}.
        '''
        return _xkwds_get(self._kwds, scaled=None)

    @property_RO
    def wrap(self):
        '''Get the C{wrap} setting or the default (C{bool}) or C{None}.
        '''
        return _xkwds_get(self._kwds, wrap=self._wrap)


class HeightIDWcosineAndoyerLambert(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW)
       and the I{angular} distance in C{radians} from function
       L{pygeodesy.cosineAndoyerLambert_}.
    '''
    def __init__(self, knots, beta=2, name=NN, **datum_wrap):
        '''New L{HeightIDWcosineAndoyerLambert} interpolator.

           @kwarg datum_wrap: Optional keyword arguments for function
                              L{pygeodesy.cosineAndoyerLambert}.

           @see: L{_HeightIDW.__init__} for details about B{C{knots}},
                 B{C{beta}}, B{C{name}} and other exceptions.
        '''
        _HeightIDW.__init__(self, knots, beta=beta, name=name, **datum_wrap)
        self._func = cosineAndoyerLambert

    if _FOR_DOCS:
        __call__ = _HeightIDW.__call__
        height   = _HeightIDW.height


class HeightIDWcosineForsytheAndoyerLambert(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW)
       and the I{angular} distance in C{radians} from function
       L{pygeodesy.cosineForsytheAndoyerLambert_}.
    '''
    def __init__(self, knots, beta=2, name=NN, **datum_wrap):
        '''New L{HeightIDWcosineForsytheAndoyerLambert} interpolator.

           @kwarg datum_wrap: Optional keyword arguments for function
                              L{pygeodesy.cosineAndoyerLambert}.

           @see: L{_HeightIDW.__init__} for details about B{C{knots}},
                 B{C{beta}}, B{C{name}} and other exceptions.
        '''
        _HeightIDW.__init__(self, knots, beta=beta, name=name, **datum_wrap)
        self._func = cosineForsytheAndoyerLambert

    if _FOR_DOCS:
        __call__ = _HeightIDW.__call__
        height   = _HeightIDW.height


class HeightIDWcosineLaw(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW) and the
       I{angular} distance in C{radians} from function L{pygeodesy.cosineLaw_}.

       @note: See note at function L{pygeodesy.vincentys_}.
    '''
    def __init__(self, knots, beta=2, name=NN, **radius_wrap):
        '''New L{HeightIDWcosineLaw} interpolator.

           @kwarg radius_wrap: Optional keyword arguments for function
                              L{pygeodesy.cosineLaw}.

           @see: L{_HeightIDW.__init__} for details about B{C{knots}},
                 B{C{beta}}, B{C{name}} and other exceptions.
        '''
        _HeightIDW.__init__(self, knots, beta=beta, name=name, **radius_wrap)
        self._func = cosineLaw

    if _FOR_DOCS:
        __call__ = _HeightIDW.__call__
        height   = _HeightIDW.height


class HeightIDWdistanceTo(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW)
       and the distance from the points' C{LatLon.distanceTo} method,
       conventionally in C{meter}.
    '''
    def __init__(self, knots, beta=2, name=NN, **distanceTo_kwds):
        '''New L{HeightIDWdistanceTo} interpolator.

           @kwarg distanceTo_kwds: Optional keyword arguments for the
                                   B{C{knots}}' C{LatLon.distanceTo}
                                   method.

           @see: L{_HeightIDW.__init__} for details about B{C{knots}},
                 B{C{beta}}, B{C{name}} and other exceptions.

           @note: All B{C{points}} I{must} be instances of the same
                  ellipsoidal or spherical C{LatLon} class, I{not
                  checked however}.
        '''
        _HeightIDW.__init__(self, knots, beta=beta, name=name,
                                       **distanceTo_kwds)
        ks = _distanceTo(HeightError, knots=self._ks)
        # use knots[0] class and datum to create
        # compatible points in _height_called
        # instead of class LatLon_ and datum None
        self._datum = ks[0].datum
        self._LLis  = ks[0].classof

    def _distances(self, x, y):
        '''(INTERNAL) Yield distances to C{(x, y)}.
        '''
        kwds, ll = self._kwds, self._LLis(y, x)
        for _, k in enumerate(self._ks):
            yield k.distanceTo(ll, **kwds)

    if _FOR_DOCS:
        __call__ = _HeightIDW.__call__
        height   = _HeightIDW.height


class HeightIDWequirectangular(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW) and
       the I{angular} distance in C{radians squared} like function
       L{pygeodesy.equirectangular_}.
    '''
    def __init__(self, knots, beta=2, name=NN, **adjust_limit_wrap):  # XXX beta=1
        '''New L{HeightIDWequirectangular} interpolator.

           @kwarg adjust_limit_wrap: Optional keyword arguments for
                               function L{pygeodesy.equirectangular_}.

           @see: L{_HeightIDW.__init__} for details about B{C{knots}},
                 B{C{beta}}, B{C{name}} and other exceptions.
        '''
        _HeightIDW.__init__(self, knots, beta=beta, name=name,
                                       **adjust_limit_wrap)

    def _distances(self, x, y):
        '''(INTERNAL) Yield distances to C{(x, y)}.
        '''
        _f, kwds = equirectangular_, self._kwds
        for _, k in enumerate(self._ks):
            yield _f(y, x, k.lat, k.lon, **kwds).distance2

    if _FOR_DOCS:
        __call__ = _HeightIDW.__call__
        height   = _HeightIDW.height


class HeightIDWeuclidean(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW) and the
       I{angular} distance in C{radians} from function L{pygeodesy.euclidean_}.
    '''
    def __init__(self, knots, beta=2, name=NN, **adjust_radius_wrap):
        '''New L{HeightIDWeuclidean} interpolator.

           @kwarg adjust_radius_wrap: Optional keyword arguments for
                                function L{pygeodesy.euclidean_}.

           @see: L{_HeightIDW.__init__} for details about B{C{knots}},
                 B{C{beta}}, B{C{name}} and other exceptions.
        '''
        _HeightIDW.__init__(self, knots, beta=beta, name=name,
                                       **adjust_radius_wrap)
        self._func = euclidean

    if _FOR_DOCS:
        __call__ = _HeightIDW.__call__
        height   = _HeightIDW.height


class HeightIDWexact(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW)
       and the I{angular} distance in C{degrees} from method
       L{GeodesicExact}C{.Inverse}.
    '''
    def __init__(self, knots, beta=2, name=NN, datum=None, **wrap):
        '''New L{HeightIDWexact} interpolator.

           @kwarg datum: Datum to override the default C{Datums.WGS84} and
                         first B{C{knots}}' datum (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg wrap: Optional keyword argument for method C{Inverse1}
                        of class L{geodesicx.GeodesicExact}.

           @raise TypeError: Invalid B{C{datum}}.

           @see: L{_HeightIDW.__init__} for details about B{C{knots}},
                 B{C{beta}}, B{C{name}} and other exceptions.
        '''
        _HeightIDW.__init__(self, knots, beta=beta, name=name, **wrap)
        self._datum_setter(datum)
        self._func = self.datum.ellipsoid.geodesicx.Inverse1

    if _FOR_DOCS:
        __call__ = _HeightIDW.__call__
        height   = _HeightIDW.height


class HeightIDWflatLocal(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW)
       and the I{angular} distance in C{radians squared} like function
       L{pygeodesy.flatLocal_}/L{pygeodesy.hubeny_}.
    '''
    def __init__(self, knots, beta=2, name=NN, **datum_hypot_scaled_wrap):
        '''New L{HeightIDWflatLocal}/L{HeightIDWhubeny} interpolator.

           @kwarg datum_hypot_scaled_wrap: Optional keyword arguments for
                              function L{pygeodesy.flatLocal_}.

           @see: L{_HeightIDW.__init__} for details about B{C{knots}},
                 B{C{beta}}, B{C{name}} and other exceptions.
        '''
        _HeightIDW.__init__(self, knots, beta=beta, name=name,
                                       **datum_hypot_scaled_wrap)
        self._func = flatLocal

    if _FOR_DOCS:
        __call__ = _HeightIDW.__call__
        height   = _HeightIDW.height


class HeightIDWflatPolar(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW)
       and the I{angular} distance in C{radians} from function
       L{pygeodesy.flatPolar_}.
    '''
    def __init__(self, knots, beta=2, name=NN, **radius_wrap):
        '''New L{HeightIDWflatPolar} interpolator.

           @kwarg radius_wrap: Optional keyword arguments for function
                               L{pygeodesy.flatPolar}.

           @see: L{_HeightIDW.__init__} for details about B{C{knots}},
                 B{C{beta}}, B{C{name}} and other exceptions.
        '''
        _HeightIDW.__init__(self, knots, beta=beta, name=name, **radius_wrap)
        self._func = flatPolar

    if _FOR_DOCS:
        __call__ = _HeightIDW.__call__
        height   = _HeightIDW.height


class HeightIDWhaversine(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW) and the
       I{angular} distance in C{radians} from function L{pygeodesy.haversine_}.

       @note: See note at function L{pygeodesy.vincentys_}.
    '''
    def __init__(self, knots, beta=2, name=NN, **radius_wrap):
        '''New L{HeightIDWhaversine} interpolator.

           @kwarg radius_wrap: Optional keyword arguments for function
                               L{pygeodesy.haversine}.

           @see: L{_HeightIDW.__init__} for details about B{C{knots}},
                 B{C{beta}}, B{C{name}} and other exceptions.
        '''
        _HeightIDW.__init__(self, knots, beta=beta, name=name, **radius_wrap)
        self._func = haversine

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
       <https://GeographicLib.SourceForge.io/Python/doc/code.html>}
       Inverse method.
    '''
    def __init__(self, knots, beta=2, name=NN, datum=None, **wrap):
        '''New L{HeightIDWkarney} interpolator.

           @kwarg datum: Datum to override the default C{Datums.WGS84} and
                         first B{C{knots}}' datum (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg wrap: Optional keyword argument for method C{Inverse1}
                        of class L{geodesicw.Geodesic}.

           @raise ImportError: Package U{geographiclib
                  <https://PyPI.org/project/geographiclib>} missing.

           @raise TypeError: Invalid B{C{datum}}.

           @see: L{_HeightIDW.__init__} for details about B{C{knots}},
                 B{C{beta}}, B{C{name}} and other exceptions.
        '''
        _HeightIDW.__init__(self, knots, beta=beta, name=name, **wrap)
        self._datum_setter(datum)
        self._func = self.datum.ellipsoid.geodesic.Inverse1

    if _FOR_DOCS:
        __call__ = _HeightIDW.__call__
        height   = _HeightIDW.height


class HeightIDWthomas(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW) and the
       I{angular} distance in C{radians} from function L{pygeodesy.thomas_}.
    '''
    def __init__(self, knots, beta=2, name=NN, **datum_wrap):
        '''New L{HeightIDWthomas} interpolator.

           @kwarg datum_wrap: Optional keyword argument for function
                              L{pygeodesy.thomas}.

           @see: L{_HeightIDW.__init__} for details about B{C{knots}},
                 B{C{beta}}, B{C{name}} and other exceptions.
        '''
        _HeightIDW.__init__(self, knots, beta=beta, name=name, **datum_wrap)
        self._func = thomas

    if _FOR_DOCS:
        __call__ = _HeightIDW.__call__
        height   = _HeightIDW.height


class HeightIDWvincentys(_HeightIDW):
    '''Height interpolator using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW) and the
       I{angular} distance in C{radians} from function L{pygeodesy.vincentys_}.

       @note: See note at function L{pygeodesy.vincentys_}.
    '''
    def __init__(self, knots, beta=2, name=NN, **radius_wrap):
        '''New L{HeightIDWvincentys} interpolator.

           @kwarg radius_wrap: Optional keyword arguments for function
                               L{pygeodesy.vincentys}.

           @see: L{_HeightIDW.__init__} for details about B{C{knots}},
                 B{C{beta}}, B{C{name}} and other exceptions.
        '''
        _HeightIDW.__init__(self, knots, beta=beta, name=name, **radius_wrap)
        self._func = vincentys

    if _FOR_DOCS:
        __call__ = _HeightIDW.__call__
        height   = _HeightIDW.height


__all__ += _ALL_DOCS(_HeightBase, _HeightIDW, _HeightsBase)

# **) MIT License
#
# Copyright (C) 2016-2023 -- mrJean1 at Gmail -- All Rights Reserved.
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
