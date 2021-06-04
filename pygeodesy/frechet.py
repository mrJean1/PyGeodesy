
# -*- coding: utf-8 -*-

u'''Fréchet distances.

Classes L{Frechet}, L{FrechetDegrees}, L{FrechetRadians},
L{FrechetCosineAndoyerLambert}, L{FrechetCosineForsytheAndoyerLambert},
L{FrechetCosineLaw}, L{FrechetDistanceTo}< L{FrechetEquirectangular},
L{FrechetEuclidean}, L{FrechetExact}, L{FrechetFlatLocal}, L{FrechetFlatPolar},
L{FrechetHaversine}, L{FrechetHubeny}, L{FrechetKarney}, L{FrechetThomas}
and L{FrechetVincentys} to compute I{discrete} U{Fréchet
<https://WikiPedia.org/wiki/Frechet_distance>} distances between two sets
of C{LatLon}, C{NumPy}, C{tuples} or other types of points.

Only L{FrechetDistanceTo} -iff used with L{ellipsoidalKarney.LatLon}
points- and L{FrechetKarney} requires installation of I{Charles Karney}'s
U{geographiclib<https://PyPI.org/project/geographiclib>}.

Typical usage is as follows.  First, create a C{Frechet} calculator
from one set of C{LatLon} points.

C{f = FrechetXyz(points1, ...)}

Get the I{discrete} Fréchet distance to another set of C{LatLon} points
by

C{t6 = f.discrete(points2)}

Or, use function C{frechet_} with a proper C{distance} function passed
as keyword arguments as follows

C{t6 = frechet_(points1, points2, ..., distance=...)}.

In both cases, the returned result C{t6} is a L{Frechet6Tuple}.

For C{(lat, lon, ...)} points in a C{NumPy} array or plain C{tuples},
wrap the points in a L{Numpy2LatLon} respectively L{Tuple2LatLon}
instance, more details in the documentation thereof.

For other points, create a L{Frechet} sub-class with the appropriate
C{distance} method overloading L{Frechet.distance} as in this example.

    >>> from pygeodesy import Frechet, hypot_
    >>>
    >>> class F3D(Frechet):
    >>>     """Custom Frechet example.
    >>>     """
    >>>     def distance(self, p1, p2):
    >>>         return hypot_(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z)
    >>>
    >>> f3D = F3D(xyz1, ..., units="...")
    >>> t6 = f3D.discrete(xyz2)

Transcribed from the original U{Computing Discrete Fréchet Distance
<https://www.kr.TUWien.ac.AT/staff/eiter/et-archive/cdtr9464.pdf>} by
Eiter, T. and Mannila, H., 1994, April 25, Technical Report CD-TR 94/64,
Information Systems Department/Christian Doppler Laboratory for Expert
Systems, Technical University Vienna, Austria.

This L{Frechet.discrete} implementation optionally generates intermediate
points for each point set separately.  For example, using keyword argument
C{fraction=0.5} adds one additional point halfway between each pair of
points.  Or using C{fraction=0.1} interpolates nine additional points
between each points pair.

The L{Frechet6Tuple} attributes C{fi1} and/or C{fi2} will be I{fractional}
indices of type C{float} if keyword argument C{fraction} is used.  Otherwise,
C{fi1} and/or C{fi2} are simply type C{int} indices into the respective
points set.

For example, C{fractional} index value 2.5 means an intermediate point
halfway between points[2] and points[3].  Use function L{fractional}
to obtain the intermediate point for a I{fractional} index in the
corresponding set of points.

The C{Fréchet} distance was introduced in 1906 by U{Maurice Fréchet
<https://WikiPedia.org/wiki/Maurice_Rene_Frechet>}, see U{reference
[6]<https://www.kr.TUWien.ac.AT/staff/eiter/et-archive/cdtr9464.pdf>}.
It is a measure of similarity between curves that takes into account the
location and ordering of the points.  Therefore, it is often a better metric
than the well-known C{Hausdorff} distance, see the L{hausdorff} module.
'''

from pygeodesy.basics import isscalar, _xinstanceof
from pygeodesy.datums import Datum, _WGS84
from pygeodesy.errors import _IsnotError, PointsError
from pygeodesy.fmath import hypot2
from pygeodesy.formy import cosineAndoyerLambert_, cosineForsytheAndoyerLambert_, \
                            cosineLaw_, euclidean_, flatPolar_, haversine_, \
                            thomas_, vincentys_, _scale_rad
from pygeodesy.interns import EPS, EPS1, INF, NN, _datum_, _distanceTo_, _DOT_, \
                             _n_, _points_, _units_
from pygeodesy.iters import points2 as _points2
from pygeodesy.lazily import _ALL_LAZY, _FOR_DOCS
from pygeodesy.named import _Named, _NamedTuple, notOverloaded, _Pass
from pygeodesy.namedTuples import PhiLam2Tuple
from pygeodesy.points import _fractional
from pygeodesy.props import Property_RO, property_doc_
from pygeodesy.streprs import _boolkwds, Fmt
from pygeodesy.units import FIx, Float, Number_, _Str_degrees, _Str_meter, \
                           _Str_NN, _Str_radians, _Str_radians2, _xUnit, _xUnits
from pygeodesy.utily import unrollPI

from collections import defaultdict as _defaultdict
from math import radians

__all__ = _ALL_LAZY.frechet
__version__ = '21.06.01'


def _fraction(fraction, n):
    f = 1  # int, no fractional indices
    if fraction in (None, 1):
        pass
    elif not (isscalar(fraction) and EPS < fraction < EPS1
                                 and (float(n) - fraction) < n):
        raise FrechetError(fraction=fraction)
    elif fraction < EPS1:
        f = float(fraction)
    return f


class FrechetError(PointsError):
    '''Fréchet issue.
    '''
    pass


class Frechet(_Named):
    '''Frechet base class, requires method L{Frechet.distance} to
       be overloaded.
    '''
    _adjust =  None  # not applicable
    _datum  =  None  # not applicable
    _f1     =  1
    _n1     =  0
    _ps1    =  None
    _units  = _Str_NN  # XXX Str to _Pass and for backward compatibility
    _wrap   =  None  # not applicable

    def __init__(self, points, fraction=None, name=NN, units=NN, **wrap_adjust):
        '''New C{Frechet...} calculator/interpolator.

           @arg points: First set of points (C{LatLon}[], L{Numpy2LatLon}[],
                        L{Tuple2LatLon}[] or C{other}[]).
           @kwarg fraction: Index fraction (C{float} in L{EPS}..L{EPS1}) to
                            interpolate intermediate B{C{points}} or use
                            C{None}, C{0} or C{1} for no intermediate
                            B{C{points}} and no I{fractional} indices.
           @kwarg name: Optional calculator/interpolator name (C{str}).
           @kwarg units: Optional, the distance units (C{Unit} or C{str}).
           @kwarg wrap_adjust: Optionally, C{wrap} and unroll longitudes, iff
                               applicable (C{bool}) and C{adjust} wrapped,
                               unrolled longitudinal delta by the cosine
                               of the mean latitude, iff applicable.

           @raise FrechetError: Insufficient number of B{C{points}} or
                                invalid B{C{fraction}} or B{C{wrap}} or
                                B{C{ajust}} not applicable.

        '''
        self._n1, self._ps1 = self._points2(points)
        if fraction:
            self.fraction = fraction
        if name:
            self.name = name
        if units:  # and not self.units:
            self.units = units
        if wrap_adjust:
            _boolkwds(self, **wrap_adjust)

    @Property_RO
    def adjust(self):
        '''Get the adjust setting (C{bool} or C{None} if not applicable).
        '''
        return self._adjust

    @Property_RO
    def datum(self):
        '''Get the datum (L{Datum} or C{None} if not applicable).
        '''
        return self._datum

    def _datum_setter(self, datum):
        '''(INTERNAL) Set the datum.
        '''
        d = datum or getattr(self._ps1[0], _datum_, datum)
        if d and d != self.datum:  # PYCHOK no cover
            _xinstanceof(Datum, datum=d)
            self._datum = d

    def discrete(self, points, fraction=None):
        '''Compute the C{forward, discrete Fréchet} distance.

           @arg points: Second set of points (C{LatLon}[], L{Numpy2LatLon}[],
                        L{Tuple2LatLon}[] or C{other}[]).
           @kwarg fraction: Index fraction (C{float} in L{EPS}..L{EPS1}) to
                            interpolate intermediate B{C{points}} or use
                            C{None}, C{0} or C{1} for no intermediate
                            B{C{points}} and no I{fractional} indices.

           @return: A L{Frechet6Tuple}C{(fd, fi1, fi2, r, n, units)}.

           @raise FrechetError: Insufficient number of B{C{points}} or an
                                invalid B{C{point}} or B{C{fraction}}.

           @raise RecursionError: Recursion depth exceeded, see U{sys.getrecursionlimit()
                                  <https://docs.Python.org/3/library/sys.html#sys.getrecursionlimit>}.
        '''
        n2, ps2 = self._points2(points)

        f2 = _fraction(fraction, n2)
        p2 = self.points_fraction if f2 < EPS1 else self.points_  # PYCHOK expected

        f1 = self.fraction
        p1 = self.points_fraction if f1 < EPS1 else self.points_  # PYCHOK expected

        def dF(fi1, fi2):
            return self.distance(p1(self._ps1, fi1), p2(ps2, fi2))

        try:
            return _frechet_(self._n1, f1, n2, f2, dF, self.units)
        except TypeError as x:
            t = _DOT_(self.classname, self.discrete.__name__)
            raise FrechetError(t, txt=str(x))

    def distance(self, point1, point2):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}, see function C{notOverloaded}.
        '''
        notOverloaded(self, point1, point2)

    @property_doc_(''' the index fraction (C{float}).''')
    def fraction(self):
        '''Get the index fraction (C{float} or C{1}).
        '''
        return self._f1

    @fraction.setter  # PYCHOK setter!
    def fraction(self, fraction):
        '''Set the index fraction (C{float} or C{1}).

           @arg fraction: Index fraction (C{float} in L{EPS}..L{EPS1}) to
                          interpolate intermediate B{C{points}} or use
                          C{None}, C{0} or C{1} for no intermediate
                          B{C{points}} and no I{fractional} indices.

           @raise FrechetError: Invalid B{C{fraction}}.
        '''
        self._f1 = _fraction(fraction, self._n1)

    def point(self, point):
        '''Convert a point for the C{.distance} method.

           @arg point: The point to convert ((C{LatLon}, L{Numpy2LatLon},
                       L{Tuple2LatLon} or C{other}).

           @return: The converted B{C{point}}.
        '''
        return point  # pass thru

    def points_(self, points, i):
        '''Get and convert a point for the C{.distance} method.

           @arg points: The orignal B{C{points}} to convert ((C{LatLon}[],
                        L{Numpy2LatLon}[], L{Tuple2LatLon}[] or C{other}[]).
           @arg i: The B{C{points}} index (C{int}).

           @return: The converted B{C{points[i]}}.
        '''
        return self.point(points[i])

    def _points2(self, points):
        '''(INTERNAL) Check a set of points.
        '''
        return _points2(points, closed=False, Error=FrechetError)

    def points_fraction(self, points, fi):
        '''Get and convert a I{fractional} point for the C{.distance} method.

           @arg points: The orignal B{C{points}} to convert ((C{LatLon}[],
                        L{Numpy2LatLon}[], L{Tuple2LatLon}[] or C{other}[]).
           @arg fi: The I{fractional} index in B{C{points}} (C{float} or C{int}).

           @return: The interpolated, converted, intermediate B{C{points[fi]}}.
        '''
        return self.point(_fractional(points, fi))

    @property_doc_(''' the distance units (C{Unit} or C{str}).''')
    def units(self):
        '''Get the distance units (C{Unit} or C{str}).
        '''
        return self._units

    @units.setter  # PYCHOK setter!
    def units(self, units):
        '''Set the distance units.

           @arg units: New units (C{Unit} or C{str}).

           @raise TypeError: Invalid B{C{units}}.
        '''
        self._units = _xUnits(units, Base=Float)

    @Property_RO
    def wrap(self):
        '''Get the wrap setting (C{bool} or C{None} if not applicable).
        '''
        return self._wrap


class FrechetDegrees(Frechet):
    '''L{Frechet} base class for distances in C{degrees} from
       C{LatLon} points in C{degrees}.
    '''
    _units = _Str_degrees

    if _FOR_DOCS:
        __init__ = Frechet.__init__
        discrete = Frechet.discrete


class FrechetRadians(Frechet):
    '''L{Frechet} base class for distances in C{radians} or C{radians
       squared} from C{LatLon} points converted from C{degrees} to
       C{radians}.
    '''
    _units = _Str_radians

    if _FOR_DOCS:
        __init__ = Frechet.__init__
        discrete = Frechet.discrete

    def point(self, point):
        '''Convert C{(lat, lon)} point in degrees to C{(a, b)}
           in radians.

           @return: An L{PhiLam2Tuple}C{(phi, lam)}.
        '''
        try:
            return point.philam
        except AttributeError:  # PYCHOK no cover
            return PhiLam2Tuple(radians(point.lat), radians(point.lon))


class FrechetCosineAndoyerLambert(FrechetRadians):
    '''Compute the C{Frechet} distance based on the I{angular} distance
       in C{radians} from function L{cosineAndoyerLambert_}.

       @see: L{FrechetCosineForsytheAndoyerLambert}, L{FrechetDistanceTo},
             L{FrechetExact}, L{FrechetFlatLocal}, L{FrechetHubeny},
             L{FrechetThomas} and L{FrechetKarney}.
    '''
    _datum = _WGS84
    _wrap  =  False

    def __init__(self, points, datum=None, wrap=False, fraction=None, name=NN):
        '''New L{FrechetCosineAndoyerLambert} calculator/interpolator.

           @arg points: First set of points (C{LatLon}[], L{Numpy2LatLon}[],
                        L{Tuple2LatLon}[] or C{other}[]).
           @kwarg datum: Optional datum overriding the default C{Datums.WGS84}
                         and first B{C{points}}' datum (L{Datum}).
           @kwarg wrap: Wrap and L{unrollPI} longitudes (C{bool})
           @kwarg fraction: Index fraction (C{float} in L{EPS}..L{EPS1}) to
                            interpolate intermediate B{C{points}} or use
                            C{None}, C{0} or C{1} for no intermediate
                            B{C{points}} and no L{fractional} indices.
           @kwarg name: Optional calculator/interpolator name (C{str}).

           @raise FrechetError: Insufficient number of B{C{points}} or
                                invalid B{C{fraction}}.

           @raise TypeError: Invalid B{C{datum}}.
        '''
        FrechetRadians.__init__(self, points, fraction=fraction, name=name,
                                                                 wrap=wrap)
        self._datum_setter(datum)

    if _FOR_DOCS:
        discrete = Frechet.discrete

    def distance(self, p1, p2):
        '''Return the L{cosineAndoyerLambert_} distance in C{radians}.
        '''
        r, _ = unrollPI(p1.lam, p2.lam, wrap=self._wrap)
        return cosineAndoyerLambert_(p2.phi, p1.phi, r, datum=self._datum)


class FrechetCosineForsytheAndoyerLambert(FrechetRadians):
    '''Compute the C{Frechet} distance based on the I{angular} distance
       in C{radians} from function L{cosineForsytheAndoyerLambert_}.

       @see: L{FrechetCosineAndoyerLambert}, L{FrechetDistanceTo},
             L{FrechetExact}, L{FrechetFlatLocal}, L{FrechetHubeny},
             L{FrechetThomas} and L{FrechetKarney}.
    '''
    _datum = _WGS84
    _wrap  =  False

    def __init__(self, points, datum=None, wrap=False, fraction=None, name=NN):
        '''New L{FrechetCosineForsytheAndoyerLambert} calculator/interpolator.

           @arg points: First set of points (C{LatLon}[], L{Numpy2LatLon}[],
                        L{Tuple2LatLon}[] or C{other}[]).
           @kwarg datum: Optional datum overriding the default C{Datums.WGS84}
                         and first B{C{points}}' datum (L{Datum}).
           @kwarg wrap: Wrap and L{unrollPI} longitudes (C{bool})
           @kwarg fraction: Index fraction (C{float} in L{EPS}..L{EPS1}) to
                            interpolate intermediate B{C{points}} or use
                            C{None}, C{0} or C{1} for no intermediate
                            B{C{points}} and no L{fractional} indices.
           @kwarg name: Optional calculator/interpolator name (C{str}).

           @raise FrechetError: Insufficient number of B{C{points}} or
                                invalid B{C{fraction}}.

           @raise TypeError: Invalid B{C{datum}}.
        '''
        FrechetRadians.__init__(self, points, fraction=fraction, name=name,
                                                                 wrap=wrap)
        self._datum_setter(datum)

    if _FOR_DOCS:
        discrete = Frechet.discrete

    def distance(self, p1, p2):
        '''Return the L{cosineForsytheAndoyerLambert_} distance in C{radians}.
        '''
        r, _ = unrollPI(p1.lam, p2.lam, wrap=self._wrap)
        return cosineForsytheAndoyerLambert_(p2.phi, p1.phi, r, datum=self._datum)


class FrechetCosineLaw(FrechetRadians):
    '''Compute the C{Frechet} distance based on the I{angular} distance
       in C{radians} from function L{cosineLaw_}.

       @see: L{FrechetEquirectangular}, L{FrechetEuclidean},
             L{FrechetExact}, L{FrechetFlatPolar}, L{FrechetHaversine}
             and L{FrechetVincentys}.

       @note: See note at function L{vincentys_}.
    '''
    _wrap = False

    def __init__(self, points, wrap=False, fraction=None, name=NN):
        '''New L{FrechetCosineLaw} calculator/interpolator.

           @arg points: First set of points (C{LatLon}[], L{Numpy2LatLon}[],
                        L{Tuple2LatLon}[] or C{other}[]).
           @kwarg wrap: Wrap and L{unrollPI} longitudes (C{bool}).
           @kwarg fraction: Index fraction (C{float} in L{EPS}..L{EPS1}) to
                            interpolate intermediate B{C{points}} or use
                            C{None}, C{0} or C{1} for no intermediate
                            B{C{points}} and no I{fractional} indices.
           @kwarg name: Optional calculator/interpolator name (C{str}).

           @raise FrechetError: Insufficient number of B{C{points}} or
                                invalid B{C{fraction}}.
        '''
        FrechetRadians.__init__(self, points, fraction=fraction, name=name,
                                                                 wrap=wrap)

    if _FOR_DOCS:
        discrete = Frechet.discrete

    def distance(self, p1, p2):
        '''Return the L{cosineLaw_} distance in C{radians}.
        '''
        r, _ = unrollPI(p1.lam, p2.lam, wrap=self._wrap)
        return cosineLaw_(p2.phi, p1.phi, r)


class FrechetDistanceTo(Frechet):
    '''Compute the C{Frechet} distance based on the distance from the
       points' C{LatLon.distanceTo} method, conventionally in C{meter}.

       @see: L{FrechetCosineAndoyerLambert}, L{FrechetCosineForsytheAndoyerLambert},
             L{FrechetExact}, L{FrechetFlatLocal}, L{FrechetHubeny}, L{FrechetThomas}
             and L{FrechetKarney}.
    '''
    _distanceTo_kwds =  {}
    _units           = _Str_meter

    def __init__(self, points, fraction=None, name=NN, **distanceTo_kwds):
        '''New L{FrechetDistanceTo} calculator/interpolator.

           @arg points: First set of points (C{LatLon}[]).
           @kwarg fraction: Index fraction (C{float} in L{EPS}..L{EPS1}) to
                            interpolate intermediate B{C{points}} or use
                            C{None}, C{0} or C{1} for no intermediate
                            B{C{points}} and no I{fractional} indices.
           @kwarg name: Optional calculator/interpolator name (C{str}).
           @kwarg distanceTo_kwds: Optional keyword arguments for the
                                   B{C{points}}' C{LatLon.distanceTo}
                                   method.

           @raise FrechetError: Insufficient number of B{C{points}} or an
                                invalid B{C{point}} or B{C{fraction}}.

           @raise ImportError: Package U{geographiclib
                  <https://PyPI.org/project/geographiclib>} missing
                  iff B{C{points}} are L{ellipsoidalKarney.LatLon}s.

           @note: All B{C{points}} I{must} be instances of the same
                  ellipsoidal or spherical C{LatLon} class.
        '''
        FrechetRadians.__init__(self, points, fraction=fraction, name=name)
        if distanceTo_kwds:
            self._distanceTo_kwds = distanceTo_kwds

    if _FOR_DOCS:
        discrete = Frechet.discrete

    def distance(self, p1, p2):
        '''Return the distance in C{meter}.
        '''
        return p1.distanceTo(p2, **self._distanceTo_kwds)

    def _points2(self, points):
        '''(INTERNAL) Check a set of points.
        '''
        np, ps = Frechet._points2(self, points)
        for i, p in enumerate(ps):
            if not callable(getattr(p, _distanceTo_, None)):
                i = Fmt.SQUARE(_points_, i)
                raise FrechetError(i, p, txt=_distanceTo_)
        return np, ps


class FrechetEquirectangular(FrechetRadians):
    '''Compute the C{Frechet} distance based on the C{equirectangular}
       distance in C{radians squared} like function L{equirectangular_}.

       @see: L{FrechetCosineLaw}, L{FrechetEuclidean}, L{FrechetExact},
             L{FrechetFlatPolar}, L{FrechetHaversine} and L{FrechetVincentys}.
    '''
    _adjust =  True
    _units  = _Str_radians2
    _wrap   =  False

    def __init__(self, points, adjust=True, wrap=False, fraction=None, name=NN):
        '''New L{FrechetEquirectangular} calculator/interpolator.

           @arg points: First set of points (C{LatLon}[], L{Numpy2LatLon}[],
                        L{Tuple2LatLon}[] or C{other}[]).
           @kwarg adjust: Adjust the wrapped, unrolled longitudinal
                          delta by the cosine of the mean latitude (C{bool}).
           @kwarg wrap: Wrap and L{unrollPI} longitudes (C{bool}).
           @kwarg fraction: Index fraction (C{float} in L{EPS}..L{EPS1}) to
                            interpolate intermediate B{C{points}} or use
                            C{None}, C{0} or C{1} for no intermediate
                            B{C{points}} and no L{fractional} indices.
           @kwarg name: Optional calculator/interpolator name (C{str}).

           @raise FrechetError: Insufficient number of B{C{points}} or
                                invalid B{C{adjust}} or B{C{seed}}.
        '''
        FrechetRadians.__init__(self, points, fraction=fraction, name=name,
                                                adjust=adjust,   wrap=wrap)

    if _FOR_DOCS:
        discrete = Frechet.discrete

    def distance(self, p1, p2):
        '''Return the L{equirectangular_} distance in C{radians squared}.
        '''
        r, _ = unrollPI(p1.lam, p2.lam, wrap=self._wrap)
        if self._adjust:
            r *= _scale_rad(p1.phi, p2.phi)
        return hypot2(r, p2.phi - p1.phi)  # like equirectangular_ d2


class FrechetEuclidean(FrechetRadians):
    '''Compute the C{Frechet} distance based on the C{Euclidean}
       distance in C{radians} from function L{euclidean_}.

       @see: L{FrechetCosineLaw}, L{FrechetEquirectangular},
             L{FrechetExact}, L{FrechetFlatPolar}, L{FrechetHaversine}
             and L{FrechetVincentys}.
    '''
    _adjust = True
    _wrap   = True  # fixed

    def __init__(self, points, adjust=True, fraction=None, name=NN):
        '''New L{FrechetEuclidean} calculator/interpolator.

           @arg points: First set of points (C{LatLon}[], L{Numpy2LatLon}[],
                        L{Tuple2LatLon}[] or C{other}[]).
           @kwarg adjust: Adjust the wrapped, unrolled longitudinal
                          delta by the cosine of the mean latitude (C{bool}).
           @kwarg fraction: Index fraction (C{float} in L{EPS}..L{EPS1}) to
                            interpolate intermediate B{C{points}} or use
                            C{None}, C{0} or C{1} for no intermediate
                            B{C{points}} and no L{fractional} indices.
           @kwarg name: Optional calculator/interpolator name (C{str}).

           @raise FrechetError: Insufficient number of B{C{points}} or
                                invalid B{C{fraction}}.
        '''
        FrechetRadians.__init__(self, points, fraction=fraction, name=name,
                                                adjust=adjust)

    if _FOR_DOCS:
        discrete = Frechet.discrete

    def distance(self, p1, p2):
        '''Return the L{euclidean_} distance in C{radians}.
        '''
        return euclidean_(p2.phi, p1.phi, p2.lam - p1.lam, adjust=self._adjust)


class FrechetExact(FrechetDegrees):
    '''Compute the C{Frechet} distance based on the I{angular}
       distance in C{degrees} from method L{GeodesicExact}C{.Inverse}.

       @see: L{FrechetCosineAndoyerLambert}, L{FrechetCosineForsytheAndoyerLambert},
             L{FrechetDistanceTo}, L{FrechetFlatLocal}, L{FrechetHubeny},
             L{FrechetKarney} and L{FrechetThomas}.
    '''
    _datum    = _WGS84
    _Inverse1 =  None
    _wrap     =  False

    def __init__(self, points, datum=None, wrap=False, fraction=None, name=NN):
        '''New L{FrechetExact} calculator/interpolator.

           @arg points: First set of points (C{LatLon}[], L{Numpy2LatLon}[],
                        L{Tuple2LatLon}[] or C{other}[]).
           @kwarg datum: Optional datum overriding the default C{Datums.WGS84}
                         and first B{C{points}}' datum (L{Datum}).
           @kwarg wrap: Wrap and L{unroll180} longitudes (C{bool})
           @kwarg fraction: Index fraction (C{float} in L{EPS}..L{EPS1}) to
                            interpolate intermediate B{C{points}} or use
                            C{None}, C{0} or C{1} for no intermediate
                            B{C{points}} and no L{fractional} indices.
           @kwarg name: Optional calculator/interpolator name (C{str}).

           @raise FrechetError: Insufficient number of B{C{points}} or
                                invalid B{C{fraction}}.

           @raise TypeError: Invalid B{C{datum}}.
        '''
        FrechetDegrees.__init__(self, points, fraction=fraction, name=name,
                                                                 wrap=wrap)
        self._datum_setter(datum)
        self._Inverse1 = self.datum.ellipsoid.geodesicx.Inverse1  # note -x

    if _FOR_DOCS:
        discrete = Frechet.discrete

    def distance(self, p1, p2):
        '''Return the non-negative I{angular} distance in C{degrees}.
        '''
        return self._Inverse1(p1.lat, p1.lon, p2.lat, p2.lon, wrap=self._wrap)


class FrechetFlatLocal(FrechetRadians):
    '''Compute the C{Frechet} distance based on the I{angular} distance
       in C{radians squared} like function L{flatLocal_}/L{hubeny_}.

       @see: L{FrechetCosineAndoyerLambert}, L{FrechetCosineForsytheAndoyerLambert},
             L{FrechetDistanceTo}, L{FrechetExact}, L{FrechetHubeny},
             L{FrechetKarney} and L{FrechetThomas}.
    '''
    _datum    = _WGS84
    _hubeny2_ =  None
    _units    = _Str_radians2
    _wrap     =  False

    def __init__(self, points, datum=None, wrap=False, fraction=None, name=NN):
        '''New L{FrechetFlatLocal}/L{FrechetHubeny} calculator/interpolator.

           @arg points: First set of points (C{LatLon}[], L{Numpy2LatLon}[],
                        L{Tuple2LatLon}[] or C{other}[]).
           @kwarg datum: Optional datum overriding the default C{Datums.WGS84}
                         and first B{C{points}}' datum (L{Datum}).
           @kwarg wrap: Wrap and L{unrollPI} longitudes (C{bool})
           @kwarg fraction: Index fraction (C{float} in L{EPS}..L{EPS1}) to
                            interpolate intermediate B{C{points}} or use
                            C{None}, C{0} or C{1} for no intermediate
                            B{C{points}} and no L{fractional} indices.
           @kwarg name: Optional calculator/interpolator name (C{str}).

           @raise FrechetError: Insufficient number of B{C{points}} or
                                invalid B{C{fraction}}.

           @raise TypeError: Invalid B{C{datum}}.
        '''
        FrechetRadians.__init__(self, points, fraction=fraction, name=name,
                                                                 wrap=wrap)
        self._datum_setter(datum)
        self._hubeny2_ = self.datum.ellipsoid._hubeny2_

    if _FOR_DOCS:
        discrete = Frechet.discrete

    def distance(self, p1, p2):
        '''Return the L{flatLocal_}/L{hubeny_} distance in C{radians squared}.
        '''
        d, _ = unrollPI(p1.lam, p2.lam, wrap=self._wrap)
        return self._hubeny2_(p2.phi, p1.phi, d)


class FrechetFlatPolar(FrechetRadians):
    '''Compute the C{Frechet} distance based on the I{angular} distance
       in C{radians} from function L{flatPolar_}.

       @see: L{FrechetCosineLaw}, L{FrechetEquirectangular},
             L{FrechetEuclidean}, L{FrechetExact}, L{FrechetHaversine}
             and L{FrechetVincentys}.
    '''
    _wrap = False

    def __init__(self, points, wrap=False, fraction=None, name=NN):
        '''New L{FrechetFlatPolar} calculator/interpolator.

           @arg points: First set of points (C{LatLon}[], L{Numpy2LatLon}[],
                        L{Tuple2LatLon}[] or C{other}[]).
           @kwarg wrap: Wrap and L{unrollPI} longitudes (C{bool}).
           @kwarg fraction: Index fraction (C{float} in L{EPS}..L{EPS1}) to
                            interpolate intermediate B{C{points}} or use
                            C{None}, C{0} or C{1} for no intermediate
                            B{C{points}} and no I{fractional} indices.
           @kwarg name: Optional calculator/interpolator name (C{str}).

           @raise FrechetError: Insufficient number of B{C{points}} or
                                invalid B{C{fraction}}.
        '''
        FrechetRadians.__init__(self, points, fraction=fraction, name=name,
                                                                 wrap=wrap)

    if _FOR_DOCS:
        discrete = Frechet.discrete

    def distance(self, p1, p2):
        '''Return the L{flatPolar_} distance in C{radians}.
        '''
        d, _ = unrollPI(p1.lam, p2.lam, wrap=self._wrap)
        return flatPolar_(p2.phi, p1.phi, d)


class FrechetHaversine(FrechetRadians):
    '''Compute the C{Frechet} distance based on the I{angular}
       distance in C{radians} from function L{haversine_}.

       @see: L{FrechetCosineLaw}, L{FrechetEquirectangular},
             L{FrechetEuclidean}, L{FrechetExact}, L{FrechetFlatPolar}
             and L{FrechetVincentys}.

       @note: See note at function L{vincentys_}.
    '''
    _wrap = False

    def __init__(self, points, wrap=False, fraction=None, name=NN):
        '''New L{FrechetHaversine} calculator/interpolator.

           @arg points: First set of points (C{LatLon}[], L{Numpy2LatLon}[],
                        L{Tuple2LatLon}[] or C{other}[]).
           @kwarg wrap: Wrap and L{unrollPI} longitudes (C{bool}).
           @kwarg fraction: Index fraction (C{float} in L{EPS}..L{EPS1}) to
                            interpolate intermediate B{C{points}} or use
                            C{None}, C{0} or C{1} for no intermediate
                            B{C{points}} and no I{fractional} indices.
           @kwarg name: Optional calculator/interpolator name (C{str}).

           @raise FrechetError: Insufficient number of B{C{points}} or
                                invalid B{C{fraction}}.
        '''
        FrechetRadians.__init__(self, points, fraction=fraction, name=name,
                                                                 wrap=wrap)

    if _FOR_DOCS:
        discrete = Frechet.discrete

    def distance(self, p1, p2):
        '''Return the L{haversine_} distance in C{radians}.
        '''
        d, _ = unrollPI(p1.lam, p2.lam, wrap=self._wrap)
        return haversine_(p2.phi, p1.phi, d)


class FrechetHubeny(FrechetFlatLocal):  # for Karl Hubeny
    if _FOR_DOCS:
        __doc__  = FrechetFlatLocal.__doc__
        __init__ = FrechetFlatLocal.__init__
        discrete = FrechetFlatLocal.discrete
        distance = FrechetFlatLocal.discrete


class FrechetKarney(FrechetExact):
    '''Compute the C{Frechet} distance based on the I{angular}
       distance in C{degrees} from I{Karney}'s U{geographiclib
       <https://PyPI.org/project/geographiclib>} U{Geodesic
       <https://GeographicLib.SourceForge.io/html/python/code.html>}
       Inverse method.

       @see: L{FrechetCosineAndoyerLambert}, L{FrechetCosineForsytheAndoyerLambert},
             L{FrechetDistanceTo}, L{FrechetExact}, L{FrechetFlatLocal},
             L{FrechetHubeny} and L{FrechetThomas}.
    '''
    def __init__(self, points, datum=None, wrap=False, fraction=None, name=NN):
        '''New L{FrechetKarney} calculator/interpolator.

           @arg points: First set of points (C{LatLon}[], L{Numpy2LatLon}[],
                        L{Tuple2LatLon}[] or C{other}[]).
           @kwarg datum: Optional datum overriding the default C{Datums.WGS84}
                         and first B{C{points}}' datum (L{Datum}).
           @kwarg wrap: Wrap and L{unroll180} longitudes (C{bool})
           @kwarg fraction: Index fraction (C{float} in L{EPS}..L{EPS1}) to
                            interpolate intermediate B{C{points}} or use
                            C{None}, C{0} or C{1} for no intermediate
                            B{C{points}} and no L{fractional} indices.
           @kwarg name: Optional calculator/interpolator name (C{str}).

           @raise FrechetError: Insufficient number of B{C{points}} or
                                invalid B{C{fraction}}.

           @raise ImportError: Package U{geographiclib
                  <https://PyPI.org/project/geographiclib>} missing.

           @raise TypeError: Invalid B{C{datum}}.
        '''
        FrechetDegrees.__init__(self, points, fraction=fraction, name=name,
                                                                 wrap=wrap)
        self._datum_setter(datum)
        self._Inverse1 = self.datum.ellipsoid.geodesic.Inverse1


class FrechetThomas(FrechetRadians):
    '''Compute the C{Frechet} distance based on the I{angular} distance
       in C{radians} from function L{thomas_}.

       @see: L{FrechetCosineAndoyerLambert}, L{FrechetCosineForsytheAndoyerLambert},
             L{FrechetDistanceTo}, L{FrechetExact}, L{FrechetFlatLocal},
             L{FrechetHubeny} and L{FrechetKarney}.
    '''
    _datum = _WGS84
    _wrap  =  False

    def __init__(self, points, datum=None, wrap=False, fraction=None, name=NN):
        '''New L{FrechetThomas} calculator/interpolator.

           @arg points: First set of points (C{LatLon}[], L{Numpy2LatLon}[],
                        L{Tuple2LatLon}[] or C{other}[]).
           @kwarg datum: Optional datum overriding the default C{Datums.WGS84}
                         and first B{C{points}}' datum (L{Datum}).
           @kwarg wrap: Wrap and L{unrollPI} longitudes (C{bool})
           @kwarg fraction: Index fraction (C{float} in L{EPS}..L{EPS1}) to
                            interpolate intermediate B{C{points}} or use
                            C{None}, C{0} or C{1} for no intermediate
                            B{C{points}} and no L{fractional} indices.
           @kwarg name: Optional calculator/interpolator name (C{str}).

           @raise FrechetError: Insufficient number of B{C{points}} or
                                invalid B{C{fraction}}.

           @raise TypeError: Invalid B{C{datum}}.
        '''
        FrechetRadians.__init__(self, points, fraction=fraction, name=name,
                                                                 wrap=wrap)
        self._datum_setter(datum)

    if _FOR_DOCS:
        discrete = Frechet.discrete

    def distance(self, p1, p2):
        '''Return the L{thomas_} distance in C{radians}.
        '''
        r, _ = unrollPI(p1.lam, p2.lam, wrap=self._wrap)
        return thomas_(p2.phi, p1.phi, r, datum=self._datum)


class FrechetVincentys(FrechetRadians):
    '''Compute the C{Frechet} distance based on the I{angular}
       distance in C{radians} from function L{vincentys_}.

       @see: L{FrechetCosineLaw}, L{FrechetEquirectangular},
             L{FrechetEuclidean}, L{FrechetExact}, L{FrechetFlatPolar}
             and L{FrechetHaversine}.

       @note: See note at function L{vincentys_}.
    '''
    _wrap = False

    def __init__(self, points, wrap=False, fraction=None, name=NN):
        '''New L{FrechetVincentys} calculator/interpolator.

           @arg points: First set of points (C{LatLon}[], L{Numpy2LatLon}[],
                        L{Tuple2LatLon}[] or C{other}[]).
           @kwarg wrap: Wrap and L{unrollPI} longitudes (C{bool}).
           @kwarg fraction: Index fraction (C{float} in L{EPS}..L{EPS1}) to
                            interpolate intermediate B{C{points}} or use
                            C{None}, C{0} or C{1} for no intermediate
                            B{C{points}} and no I{fractional} indices.
           @kwarg name: Optional calculator/interpolator name (C{str}).

           @raise FrechetError: Insufficient number of B{C{points}} or
                                invalid B{C{fraction}}.
        '''
        FrechetRadians.__init__(self, points, fraction=fraction, name=name,
                                                                 wrap=wrap)

    if _FOR_DOCS:
        discrete = Frechet.discrete

    def distance(self, p1, p2):
        '''Return the L{vincentys_} distance in C{radians}.
        '''
        d, _ = unrollPI(p1.lam, p2.lam, wrap=self._wrap)
        return vincentys_(p2.phi, p1.phi, d)


def _frechet_(ni, fi, nj, fj, dF, units):  # MCCABE 14
    '''(INTERNAL) Recursive core of function L{frechet_}
       and method C{discrete} of C{Frechet...} classes.
    '''
    iFs = {}

    def iF(i):  # cache index, depth ints and floats
        return iFs.setdefault(i, i)

    cF = _defaultdict(dict)

    def rF(i, j, r):  # recursive Fréchet
        i = iF(i)
        j = iF(j)
        try:
            t = cF[i][j]
        except KeyError:
            r = iF(r + 1)
            try:
                if i > 0:
                    if j > 0:
                        t = min(rF(i - fi, j,      r),
                                rF(i - fi, j - fj, r),
                                rF(i,      j - fj, r))
                    elif j < 0:
                        raise IndexError
                    else:  # j == 0
                        t = rF(i - fi, 0, r)
                elif i < 0:
                    raise IndexError

                elif j > 0:  # i == 0
                    t = rF(0, j - fj, r)
                elif j < 0:  # i == 0
                    raise IndexError
                else:  # i == j == 0
                    t = (-INF, i, j, r)

                d = dF(i, j)
                if d > t[0]:
                    t = (d, i, j, r)
            except IndexError:
                t = (INF, i, j, r)
            cF[i][j] = t
        return t

    t  = rF(ni - 1, nj - 1, 0)
    t += (sum(map(len, cF.values())), units)
#   del cF, iFs
    return Frechet6Tuple(*t)


def frechet_(points1, points2, distance=None, units=NN):
    '''Compute the I{discrete} U{Fréchet<https://WikiPedia.org/wiki/Frechet_distance>}
       distance between two paths given as sets of points.

       @arg points1: First set of points (C{LatLon}[], L{Numpy2LatLon}[],
                     L{Tuple2LatLon}[] or C{other}[]).
       @arg points2: Second set of points (C{LatLon}[], L{Numpy2LatLon}[],
                     L{Tuple2LatLon}[] or C{other}[]).
       @kwarg distance: Callable returning the distance between a B{C{points1}}
                        and a B{C{points2}} point (signature C{(point1, point2)}).
       @kwarg units: Optional, the distance units (C{Unit} or C{str}).

       @return: A L{Frechet6Tuple}C{(fd, fi1, fi2, r, n, units)} where C{fi1}
                and C{fi2} are type C{int} indices into B{C{points1}} respectively
                B{C{points2}}.

       @raise FrechetError: Insufficient number of B{C{points1}} or B{C{points2}}.

       @raise RecursionError: Recursion depth exceeded, see U{sys.getrecursionlimit()
                              <https://docs.Python.org/3/library/sys.html#sys.getrecursionlimit>}.

       @raise TypeError: If B{C{distance}} is not a callable.

       @note: Function L{frechet_} does not support I{fractional} indices for
              intermediate B{C{points1}} and B{C{points2}}.
    '''
    if not callable(distance):
        raise _IsnotError(callable.__name__, distance=distance)

    n1, ps1 = _points2(points1, closed=False, Error=FrechetError)
    n2, ps2 = _points2(points2, closed=False, Error=FrechetError)

    def dF(i1, i2):
        return distance(ps1[i1], ps2[i2])

    return _frechet_(n1, 1, n2, 1, dF, units)


class Frechet6Tuple(_NamedTuple):
    '''6-Tuple C{(fd, fi1, fi2, r, n, units)} with the I{discrete}
       U{Fréchet<https://WikiPedia.org/wiki/Frechet_distance>} distance
       C{fd}, I{fractional} indices C{fi1} and C{fi2} as C{FIx}, the
       recursion depth C{r}, the number of distances computed C{n} and
       the L{units} class or class or name of the distance C{units}.

       If I{fractional} indices C{fi1} and C{fi2} are C{int}, the
       returned C{fd} is the distance between C{points1[fi1]} and
       C{points2[fi2]}.  For C{float} indices, the distance is
       between an intermediate point along C{points1[int(fi1)]} and
       C{points1[int(fi1) + 1]} respectively an intermediate point
       along C{points2[int(fi2)]} and C{points2[int(fi2) + 1]}.

       Use function L{fractional} to compute the point at a
       I{fractional} index.
    '''
    _Names_ = ('fd',  'fi1', 'fi2', 'r',     _n_,      _units_)
    _Units_ = (_Pass,  FIx,   FIx,   Number_, Number_, _Pass)

    def toUnits(self, **Error):  # PYCHOK expected
        '''Overloaded C{_NamedTuple.toUnits} for C{fd} units.
        '''
        U = _xUnit(self.units, Float)  # PYCHOK expected
        self._Units_ = (U,) + Frechet6Tuple._Units_[1:]
        return _NamedTuple.toUnits(self, **Error)

#   def __gt__(self, other):
#       _xinstanceof(Frechet6Tuple, other=other)
#       return self if self.fd > other.fd else other  # PYCHOK .fd=[0]
#
#   def __lt__(self, other):
#       _xinstanceof(Frechet6Tuple, other=other)
#       return self if self.fd < other.fd else other  # PYCHOK .fd=[0]

# **) MIT License
#
# Copyright (C) 2016-2021 -- mrJean1 at Gmail -- All Rights Reserved.
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
