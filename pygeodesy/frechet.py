
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

C{f = FrechetXyz(point1s, ...)}

Get the I{discrete} Fréchet distance to another set of C{LatLon} points
by

C{t6 = f.discrete(point2s)}

Or, use function C{frechet_} with a proper C{distance} function passed
as keyword arguments as follows

C{t6 = frechet_(point1s, point2s, ..., distance=...)}.

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

# from pygeodesy.basics import isscalar  # from .points
from pygeodesy.constants import EPS, EPS1, INF, NINF
from pygeodesy.datums import _ellipsoidal_datum, _WGS84
from pygeodesy.errors import _IsnotError, PointsError, _xattr, \
                             _xkwds, _xkwds_get
import pygeodesy.formy as _formy
from pygeodesy.interns import NN, _DOT_, _n_, _units_
# from pygeodesy.iters import points2 as _points2  # from .points
from pygeodesy.lazily import _ALL_LAZY, _FOR_DOCS
from pygeodesy.named import _Named, _NamedTuple, notOverloaded, _Pass
# from pygeodesy.namedTuples import PhiLam2Tuple  # from .points
from pygeodesy.points import _distanceTo, _fractional,  isscalar, \
                              PhiLam2Tuple, points2 as _points2, radians
from pygeodesy.props import property_doc_, property_RO
from pygeodesy.units import FIx, Float, Number_, _xUnit, _xUnits
from pygeodesy.unitsBase import _Str_degrees, _Str_meter, _Str_NN, \
                                _Str_radians, _Str_radians2

from collections import defaultdict as _defaultdict
# from math import radians  # from .points

__all__ = _ALL_LAZY.frechet
__version__ = '23.08.06'


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
    _datum = _WGS84
    _func  =  None  # formy function
    _f1    =  1
    _kwds  = {}     # func_ options
    _n1    =  0
    _ps1   =  None
    _units = _Str_NN  # XXX Str to _Pass and for backward compatibility

    def __init__(self, point1s, fraction=None, name=NN, units=NN, **kwds):
        '''New C{Frechet...} calculator/interpolator.

           @arg point1s: First set of points (C{LatLon}[], L{Numpy2LatLon}[],
                         L{Tuple2LatLon}[] or C{other}[]).
           @kwarg fraction: Index fraction (C{float} in L{EPS}..L{EPS1}) to
                            interpolate intermediate B{C{point1s}} or use
                            C{None}, C{0} or C{1} for no intermediate
                            B{C{point1s}} and no I{fractional} indices.
           @kwarg name: Optional calculator/interpolator name (C{str}).
           @kwarg units: Optional distance units (C{Unit} or C{str}).
           @kwarg kwds: Optional keyword argument for distance function,
                        retrievable with property C{kwds}.

           @raise FrechetError: Insufficient number of B{C{point1s}} or
                                an invalid B{C{point1}}, B{C{fraction}}
                                or B{C{units}}.
        '''
        self._n1, self._ps1 = self._points2(point1s)
        if fraction:
            self.fraction = fraction
        if name:
            self.name = name
        if units:  # and not self.units:
            self.units = units
        if kwds:
            self._kwds = kwds

    @property_RO
    def adjust(self):
        '''Get the C{adjust} setting (C{bool} or C{None}).
        '''
        return _xkwds_get(self._kwds, adjust=None)

    @property_RO
    def datum(self):
        '''Get the datum (L{Datum} or C{None} if not applicable).
        '''
        return self._datum

    def _datum_setter(self, datum):
        '''(INTERNAL) Set the datum.
        '''
        d = datum or _xattr(self._ps1[0], datum=None)
        if d and d is not self._datum:  # PYCHOK no cover
            self._datum = _ellipsoidal_datum(d, name=self.name)

    def discrete(self, point2s, fraction=None):
        '''Compute the C{forward, discrete Fréchet} distance.

           @arg point2s: Second set of points (C{LatLon}[], L{Numpy2LatLon}[],
                         L{Tuple2LatLon}[] or C{other}[]).
           @kwarg fraction: Index fraction (C{float} in L{EPS}..L{EPS1}) to
                            interpolate intermediate B{C{point2s}} or use
                            C{None}, C{0} or C{1} for no intermediate
                            B{C{point2s}} and no I{fractional} indices.

           @return: A L{Frechet6Tuple}C{(fd, fi1, fi2, r, n, units)}.

           @raise FrechetError: Insufficient number of B{C{point2s}} or
                                an invalid B{C{point2}} or B{C{fraction}}.

           @raise RecursionError: Recursion depth exceeded, see U{sys.getrecursionlimit
                                  <https://docs.Python.org/3/library/sys.html#sys.getrecursionlimit>}.
        '''
        return self._discrete(point2s, fraction, self.distance)

    def _discrete(self, point2s, fraction, distance):
        '''(INTERNAL) Detailed C{discrete} with C{disance}.
        '''
        n2, ps2 = self._points2(point2s)

        f2 = _fraction(fraction, n2)
        p2 = self.points_fraction if f2 < EPS1 else self.points_  # PYCHOK expected

        f1 = self.fraction
        p1 = self.points_fraction if f1 < EPS1 else self.points_  # PYCHOK expected

        def _dF(fi1, fi2):
            return distance(p1(self._ps1, fi1), p2(ps2, fi2))

        try:
            return _frechet_(self._n1, f1, n2, f2, _dF, self.units)
        except TypeError as x:
            t = _DOT_(self.classname, self.discrete.__name__)
            raise FrechetError(t, cause=x)

    def distance(self, point1, point2):
        '''Return the distance between B{C{point1}} and B{C{point2s}},
           subject to the supplied optional keyword arguments, see
           property C{kwds}.
        '''
        return self._func(point1.lat, point1.lon,
                          point2.lat, point2.lon, **self._kwds)

    @property_doc_(''' the index fraction (C{float}).''')
    def fraction(self):
        '''Get the index fraction (C{float} or C{1}).
        '''
        return self._f1

    @fraction.setter  # PYCHOK setter!
    def fraction(self, fraction):
        '''Set the index fraction (C{float} in C{EPS}..L{EPS1}) to interpolate
           intermediate B{C{point1s}} or use C{None}, C{0} or C{1} for no
           intermediate B{C{point1s}} and no I{fractional} indices.

           @raise FrechetError: Invalid B{C{fraction}}.
        '''
        self._f1 = _fraction(fraction, self._n1)

    def _func(self, *args, **kwds):  # PYCHOK no cover
        notOverloaded(self, *args, **kwds)

    @property_RO
    def kwds(self):
        '''Get the supplied, optional keyword arguments (C{dict}).
        '''
        return self._kwds

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

    def points_fraction(self, points, fi):
        '''Get and convert a I{fractional} point for the C{.distance} method.

           @arg points: The orignal B{C{points}} to convert ((C{LatLon}[],
                        L{Numpy2LatLon}[], L{Tuple2LatLon}[] or C{other}[]).
           @arg fi: The I{fractional} index in B{C{points}} (C{float} or C{int}).

           @return: The interpolated, converted, intermediate B{C{points[fi]}}.
        '''
        return self.point(_fractional(points, fi, None, wrap=None))  # was=self.wrap

    def _points2(self, points):
        '''(INTERNAL) Check a set of points, overloaded in L{FrechetDistanceTo}.
        '''
        return _points2(points, closed=False, Error=FrechetError)

    @property_doc_(''' the distance units (C{Unit} or C{str}).''')
    def units(self):
        '''Get the distance units (C{Unit} or C{str}).
        '''
        return self._units

    @units.setter  # PYCHOK setter!
    def units(self, units):
        '''Set the distance units (C{Unit} or C{str}).

           @raise TypeError: Invalid B{C{units}}.
        '''
        self._units = _xUnits(units, Base=Float)

    @property_RO
    def wrap(self):
        '''Get the C{wrap} setting (C{bool} or C{None}).
        '''
        return _xkwds_get(self._kwds, wrap=None)


class FrechetDegrees(Frechet):
    '''DEPRECATED, use an other C{Frechet*} class.
    '''
    _units = _Str_degrees

    if _FOR_DOCS:
        __init__ = Frechet.__init__
        discrete = Frechet.discrete

    def distance(self, point1, point2, *args, **kwds):  # PYCHOK no cover
        '''I{Must be overloaded} to return the distance between
           B{C{point1}} and B{C{point2}} in C{degrees}.
        '''
        notOverloaded(self, point1, point2, *args, **kwds)


class FrechetRadians(Frechet):
    '''DEPRECATED, use an other C{Frechet*} class.
    '''
    _units = _Str_radians

    if _FOR_DOCS:
        __init__ = Frechet.__init__
        discrete = Frechet.discrete

    def distance(self, point1, point2, *args, **kwds):  # PYCHOK no cover
        '''I{Must be overloaded} to return the distance between
           B{C{point1}} and B{C{point2}} in C{radians}.
        '''
        notOverloaded(self, point1, point2, *args, **kwds)

    def point(self, point):
        '''Return B{C{point}} as L{PhiLam2Tuple} to maintain
           I{backward compatibility} of L{FrechetRadians}.

           @return: A L{PhiLam2Tuple}C{(phi, lam)}.
        '''
        try:
            return point.philam
        except AttributeError:
            return PhiLam2Tuple(radians(point.lat), radians(point.lon))


class _FrechetMeterRadians(Frechet):
    '''(INTERNAL) Returning C{meter} or C{radians} depending on
       the optional keyword arguments supplied at instantiation
       of the C{Frechet*} sub-class.
    '''
    _units  = _Str_meter
    _units_ = _Str_radians

    def discrete(self, point2s, fraction=None):
        '''Overloaded method L{Frechet.discrete} to determine
           the distance function and units from the optional
           keyword arguments given at this instantiation, see
           property C{kwds}.

           @see: L{Frechet.discrete} for other details.
        '''
        return self._discrete(point2s, fraction, _formy._radistance(self))

    def _func_(self, *args, **kwds):  # PYCHOK no cover
        notOverloaded(self, *args, **kwds)


class FrechetCosineAndoyerLambert(_FrechetMeterRadians):
    '''Compute the C{Frechet} distance based on the I{angular} distance
       in C{radians} from function L{pygeodesy.cosineAndoyerLambert}.
    '''
    def __init__(self, point1s, fraction=None, name=NN, **datum_wrap):
        '''New L{FrechetCosineAndoyerLambert} calculator/interpolator.

           @kwarg datum_wrap: Optional keyword arguments for function
                              L{pygeodesy.cosineAndoyerLambert}.

           @see: L{Frechet.__init__} for details about B{C{point1s}},
                 B{C{fraction}}, B{C{name}} and other exceptions.
        '''
        Frechet.__init__(self, point1s, fraction=fraction, name=name,
                                      **datum_wrap)
        self._func  = _formy.cosineAndoyerLambert
        self._func_ = _formy.cosineAndoyerLambert_

    if _FOR_DOCS:
        discrete = Frechet.discrete


class FrechetCosineForsytheAndoyerLambert(_FrechetMeterRadians):
    '''Compute the C{Frechet} distance based on the I{angular} distance
       in C{radians} from function L{pygeodesy.cosineForsytheAndoyerLambert}.
    '''
    def __init__(self, point1s, fraction=None, name=NN, **datum_wrap):
        '''New L{FrechetCosineForsytheAndoyerLambert} calculator/interpolator.

           @kwarg datum_wrap: Optional keyword arguments for function
                              L{pygeodesy.cosineAndoyerLambert}.

           @see: L{Frechet.__init__} for details about B{C{point1s}},
                 B{C{fraction}}, B{C{name}} and other exceptions.
        '''
        Frechet.__init__(self, point1s, fraction=fraction, name=name,
                                      **datum_wrap)
        self._func  = _formy.cosineForsytheAndoyerLambert
        self._func_ = _formy.cosineForsytheAndoyerLambert_

    if _FOR_DOCS:
        discrete = Frechet.discrete


class FrechetCosineLaw(_FrechetMeterRadians):
    '''Compute the C{Frechet} distance based on the I{angular} distance
       in C{radians} from function L{pygeodesy.cosineLaw}.

       @note: See note at function L{pygeodesy.vincentys_}.
    '''
    def __init__(self, point1s, fraction=None, name=NN, **radius_wrap):
        '''New L{FrechetCosineLaw} calculator/interpolator.

           @kwarg radius_wrap: Optional keyword arguments for function
                               L{pygeodesy.cosineLaw}.

           @see: L{Frechet.__init__} for details about B{C{point1s}},
                 B{C{fraction}}, B{C{name}} and other exceptions.
        '''
        Frechet.__init__(self, point1s, fraction=fraction, name=name,
                                      **radius_wrap)
        self._func  = _formy.cosineLaw
        self._func_ = _formy.cosineLaw_

    if _FOR_DOCS:
        discrete = Frechet.discrete


class FrechetDistanceTo(Frechet):  # FrechetMeter
    '''Compute the C{Frechet} distance based on the distance from the
       point1s' C{LatLon.distanceTo} method, conventionally in C{meter}.
    '''
    _units = _Str_meter

    def __init__(self, point1s, fraction=None, name=NN, **distanceTo_kwds):
        '''New L{FrechetDistanceTo} calculator/interpolator.

           @kwarg distanceTo_kwds: Optional keyword arguments for each
                                   B{C{point1s}}' C{LatLon.distanceTo}
                                   method.

           @see: L{Frechet.__init__} for details about B{C{point1s}},
                 B{C{fraction}}, B{C{name}} and other exceptions.

           @note: All B{C{point1s}} I{must} be instances of the same
                  ellipsoidal or spherical C{LatLon} class.
        '''
        Frechet.__init__(self, point1s, fraction=fraction, name=name,
                                      **distanceTo_kwds)

    if _FOR_DOCS:
        discrete = Frechet.discrete

    def distance(self, p1, p2):
        '''Return the distance in C{meter}.
        '''
        return p1.distanceTo(p2, **self._kwds)

    def _points2(self, points):
        '''(INTERNAL) Check a set of points.
        '''
        np, ps = Frechet._points2(self, points)
        return np, _distanceTo(FrechetError, points=ps)


class FrechetEquirectangular(Frechet):
    '''Compute the C{Frechet} distance based on the I{equirectangular}
       distance in C{radians squared} like function L{pygeodesy.equirectangular}.
    '''
    _units = _Str_radians2

    def __init__(self, point1s, fraction=None, name=NN, **adjust_limit_wrap):
        '''New L{FrechetEquirectangular} calculator/interpolator.

           @kwarg adjust_limit_wrap: Optional keyword arguments for function
                               L{pygeodesy.equirectangular_} I{with default}
                               C{B{limit}=0} for I{backward compatibility}.

           @see: L{Frechet.__init__} for details about B{C{point1s}},
                 B{C{fraction}}, B{C{name}} and other exceptions.
        '''
        adjust_limit_wrap = _xkwds(adjust_limit_wrap, limit=0)
        Frechet.__init__(self, point1s, fraction=fraction, name=name,
                                      **adjust_limit_wrap)
        self._func = _formy._equirectangular  # helper

    if _FOR_DOCS:
        discrete = Frechet.discrete


class FrechetEuclidean(_FrechetMeterRadians):
    '''Compute the C{Frechet} distance based on the I{Euclidean}
       distance in C{radians} from function L{pygeodesy.euclidean}.
    '''
    def __init__(self, point1s, fraction=None, name=NN, **adjust_radius_wrap):  # was=True
        '''New L{FrechetEuclidean} calculator/interpolator.

           @kwarg adjust_radius_wrap: Optional keyword arguments for
                                function L{pygeodesy.euclidean}.

           @see: L{Frechet.__init__} for details about B{C{point1s}},
                 B{C{fraction}}, B{C{name}} and other exceptions.
        '''
        Frechet.__init__(self, point1s, fraction=fraction, name=name,
                                      **adjust_radius_wrap)
        self._func  = _formy.euclidean
        self._func_ = _formy.euclidean_

    if _FOR_DOCS:
        discrete = Frechet.discrete


class FrechetExact(Frechet):
    '''Compute the C{Frechet} distance based on the I{angular} distance
       in C{degrees} from method L{GeodesicExact}C{.Inverse}.
    '''
    _units = _Str_degrees

    def __init__(self, point1s, fraction=None, name=NN, datum=None, **wrap):
        '''New L{FrechetExact} calculator/interpolator.

           @kwarg datum: Datum to override the default C{Datums.WGS84} and
                         first B{C{point1s}}' datum (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg wrap: Optional keyword argument for method C{Inverse1}
                        of class L{geodesicx.GeodesicExact}.

           @raise TypeError: Invalid B{C{datum}}.

           @see: L{Frechet.__init__} for details about B{C{point1s}},
                 B{C{fraction}}, B{C{name}} and other exceptions.
        '''
        Frechet.__init__(self, point1s, fraction=fraction, name=name,
                                      **wrap)
        self._datum_setter(datum)
        self._func = self.datum.ellipsoid.geodesicx.Inverse1  # note -x

    if _FOR_DOCS:
        discrete = Frechet.discrete


class FrechetFlatLocal(_FrechetMeterRadians):
    '''Compute the C{Frechet} distance based on the I{angular} distance in
       C{radians squared} like function L{pygeodesy.flatLocal_}/L{pygeodesy.hubeny}.
    '''
    _units_ = _Str_radians2  # see L{flatLocal_}

    def __init__(self, point1s, fraction=None, name=NN, **datum_scaled_wrap):
        '''New L{FrechetFlatLocal}/L{FrechetHubeny} calculator/interpolator.

           @kwarg datum_scaled_wrap: Optional keyword arguments for
                                     function L{pygeodesy.flatLocal}.

           @see: L{Frechet.__init__} for details about B{C{point1s}},
                 B{C{fraction}}, B{C{name}} and other exceptions.

           @note: The distance C{units} are C{radians squared}, not C{radians}.
        '''
        Frechet.__init__(self, point1s, fraction=fraction, name=name,
                                      **datum_scaled_wrap)
        self._func  = _formy.flatLocal
        self._func_ =  self.datum.ellipsoid._hubeny_2

    if _FOR_DOCS:
        discrete = Frechet.discrete


class FrechetFlatPolar(_FrechetMeterRadians):
    '''Compute the C{Frechet} distance based on the I{angular} distance
       in C{radians} from function L{flatPolar_}.
    '''
    def __init__(self, point1s, fraction=None, name=NN, **radius_wrap):
        '''New L{FrechetFlatPolar} calculator/interpolator.

           @kwarg radius_wrap: Optional keyword arguments for function
                               L{pygeodesy.flatPolar}.

           @see: L{Frechet.__init__} for details about B{C{point1s}},
                 B{C{fraction}}, B{C{name}} and other exceptions.
       '''
        Frechet.__init__(self, point1s, fraction=fraction, name=name,
                                      **radius_wrap)
        self._func  = _formy.flatPolar
        self._func_ = _formy.flatPolar_

    if _FOR_DOCS:
        discrete = Frechet.discrete


class FrechetHaversine(_FrechetMeterRadians):
    '''Compute the C{Frechet} distance based on the I{angular}
       distance in C{radians} from function L{pygeodesy.haversine_}.

       @note: See note at function L{pygeodesy.vincentys_}.
    '''
    def __init__(self, point1s, fraction=None, name=NN, **radius_wrap):
        '''New L{FrechetHaversine} calculator/interpolator.

           @kwarg radius_wrap: Optional keyword arguments for function
                               L{pygeodesy.haversine}.

           @see: L{Frechet.__init__} for details about B{C{point1s}},
                 B{C{fraction}}, B{C{name}} and other exceptions.
        '''
        Frechet.__init__(self, point1s, fraction=fraction, name=name,
                                      **radius_wrap)
        self._func  = _formy.haversine
        self._func_ = _formy.haversine_

    if _FOR_DOCS:
        discrete = Frechet.discrete


class FrechetHubeny(FrechetFlatLocal):  # for Karl Hubeny
    if _FOR_DOCS:
        __doc__  = FrechetFlatLocal.__doc__
        __init__ = FrechetFlatLocal.__init__
        discrete = FrechetFlatLocal.discrete
        distance = FrechetFlatLocal.discrete


class FrechetKarney(Frechet):
    '''Compute the C{Frechet} distance based on the I{angular}
       distance in C{degrees} from I{Karney}'s U{geographiclib
       <https://PyPI.org/project/geographiclib>} U{geodesic.Geodesic
       <https://GeographicLib.SourceForge.io/Python/doc/code.html>}
       C{Inverse} method.
    '''
    _units = _Str_degrees

    def __init__(self, point1s, fraction=None, name=NN, datum=None, **wrap):
        '''New L{FrechetKarney} calculator/interpolator.

           @kwarg datum: Datum to override the default C{Datums.WGS84} and
                         first B{C{knots}}' datum (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg wrap: Optional keyword argument for method C{Inverse1}
                        of class L{geodesicw.Geodesic}.

           @raise ImportError: Package U{geographiclib
                  <https://PyPI.org/project/geographiclib>} missing.

           @raise TypeError: Invalid B{C{datum}}.

           @see: L{Frechet.__init__} for details about B{C{point1s}},
                 B{C{fraction}}, B{C{name}} and other exceptions.
        '''
        Frechet.__init__(self, point1s, fraction=fraction, name=name,
                                      **wrap)
        self._datum_setter(datum)
        self._func = self.datum.ellipsoid.geodesic.Inverse1

    if _FOR_DOCS:
        discrete = Frechet.discrete


class FrechetThomas(_FrechetMeterRadians):
    '''Compute the C{Frechet} distance based on the I{angular} distance
       in C{radians} from function L{pygeodesy.thomas_}.
    '''
    def __init__(self, point1s, fraction=None, name=NN, **datum_wrap):
        '''New L{FrechetThomas} calculator/interpolator.

           @kwarg datum_wrap: Optional keyword argument for function
                              L{pygeodesy.thomas}.

           @see: L{Frechet.__init__} for details about B{C{point1s}},
                 B{C{fraction}}, B{C{name}} and other exceptions.
        '''
        Frechet.__init__(self, point1s, fraction=fraction, name=name,
                                      **datum_wrap)
        self._func  = _formy.thomas
        self._func_ = _formy.thomas_

    if _FOR_DOCS:
        discrete = Frechet.discrete


class FrechetVincentys(_FrechetMeterRadians):
    '''Compute the C{Frechet} distance based on the I{angular}
       distance in C{radians} from function L{pygeodesy.vincentys_}.

       @note: See note at function L{pygeodesy.vincentys_}.
    '''
    def __init__(self, point1s, fraction=None, name=NN, **radius_wrap):
        '''New L{FrechetVincentys} calculator/interpolator.

           @kwarg radius_wrap: Optional keyword arguments for function
                               L{pygeodesy.vincentys}.

           @see: L{Frechet.__init__} for details about B{C{point1s}},
                 B{C{fraction}}, B{C{name}} and other exceptions.
        '''
        Frechet.__init__(self, point1s, fraction=fraction, name=name,
                                      **radius_wrap)
        self._func  = _formy.vincentys
        self._func_ = _formy.vincentys_

    if _FOR_DOCS:
        discrete = Frechet.discrete


def _frechet_(ni, fi, nj, fj, dF, units):  # MCCABE 14
    '''(INTERNAL) Recursive core of function L{frechet_}
       and method C{discrete} of C{Frechet...} classes.
    '''
    iFs = {}

    def iF(i):  # cache index, depth ints and floats
        return iFs.setdefault(i, i)

    cF = _defaultdict(dict)

    def _rF(i, j, r):  # recursive Fréchet
        i = iF(i)
        j = iF(j)
        try:
            t = cF[i][j]
        except KeyError:
            r = iF(r + 1)
            try:
                if i > 0:
                    if j > 0:
                        t = min(_rF(i - fi, j,      r),
                                _rF(i - fi, j - fj, r),
                                _rF(i,      j - fj, r))
                    elif j < 0:
                        raise IndexError
                    else:  # j == 0
                        t = _rF(i - fi, 0, r)
                elif i < 0:
                    raise IndexError

                elif j > 0:  # i == 0
                    t = _rF(0, j - fj, r)
                elif j < 0:  # i == 0
                    raise IndexError
                else:  # i == j == 0
                    t = (NINF, i, j, r)

                d = dF(i, j)
                if d > t[0]:
                    t = (d, i, j, r)
            except IndexError:
                t = (INF, i, j, r)
            cF[i][j] = t
        return t

    t  = _rF(ni - 1, nj - 1, 0)
    t += (sum(map(len, cF.values())), units)
#   del cF, iFs
    return Frechet6Tuple(t)  # *t


def frechet_(point1s, point2s, distance=None, units=NN):
    '''Compute the I{discrete} U{Fréchet<https://WikiPedia.org/wiki/Frechet_distance>}
       distance between two paths, each given as a set of points.

       @arg point1s: First set of points (C{LatLon}[], L{Numpy2LatLon}[],
                     L{Tuple2LatLon}[] or C{other}[]).
       @arg point2s: Second set of points (C{LatLon}[], L{Numpy2LatLon}[],
                     L{Tuple2LatLon}[] or C{other}[]).
       @kwarg distance: Callable returning the distance between a B{C{point1s}}
                        and a B{C{point2s}} point (signature C{(point1, point2)}).
       @kwarg units: Optional, the distance units (C{Unit} or C{str}).

       @return: A L{Frechet6Tuple}C{(fd, fi1, fi2, r, n, units)} where C{fi1}
                and C{fi2} are type C{int} indices into B{C{point1s}} respectively
                B{C{point2s}}.

       @raise FrechetError: Insufficient number of B{C{point1s}} or B{C{point2s}}.

       @raise RecursionError: Recursion depth exceeded, see U{sys.getrecursionlimit()
                              <https://docs.Python.org/3/library/sys.html#sys.getrecursionlimit>}.

       @raise TypeError: If B{C{distance}} is not a callable.

       @note: Function L{frechet_} does I{not} support I{fractional} indices
              for intermediate B{C{point1s}} and B{C{point2s}}.
    '''
    if not callable(distance):
        raise _IsnotError(callable.__name__, distance=distance)

    n1, ps1 = _points2(point1s, closed=False, Error=FrechetError)
    n2, ps2 = _points2(point2s, closed=False, Error=FrechetError)

    def _dF(i1, i2):
        return distance(ps1[i1], ps2[i2])

    return _frechet_(n1, 1, n2, 1, _dF, units)


class Frechet6Tuple(_NamedTuple):
    '''6-Tuple C{(fd, fi1, fi2, r, n, units)} with the I{discrete}
       U{Fréchet<https://WikiPedia.org/wiki/Frechet_distance>} distance
       C{fd}, I{fractional} indices C{fi1} and C{fi2} as C{FIx}, the
       recursion depth C{r}, the number of distances computed C{n} and
       the L{units} class or class or name of the distance C{units}.

       If I{fractional} indices C{fi1} and C{fi2} are C{int}, the
       returned C{fd} is the distance between C{point1s[fi1]} and
       C{point2s[fi2]}.  For C{float} indices, the distance is
       between an intermediate point along C{point1s[int(fi1)]} and
       C{point1s[int(fi1) + 1]} respectively an intermediate point
       along C{point2s[int(fi2)]} and C{point2s[int(fi2) + 1]}.

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
