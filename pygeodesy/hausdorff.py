
# -*- coding: utf-8 -*-

u'''Hausdorff distances.

Classes L{Hausdorff}, L{HausdorffDegrees}, L{HausdorffRadians},
L{HausdorffCosineAndoyerLambert}, L{HausdorffCosineForsytheAndoyerLambert},
L{HausdorffCosineLaw}, L{HausdorffDistanceTo}, L{HausdorffEquirectangular},
L{HausdorffEuclidean}, L{HausdorffFlatLocal}, L{HausdorffFlatPolar},
L{HausdorffHaversine}, L{HausdorffHubeny}, L{HausdorffKarney},
L{HausdorffThomas} and L{HausdorffVincentys} to compute U{Hausdorff
<https://WikiPedia.org/wiki/Hausdorff_distance>} distances between two
sets of C{LatLon}, C{NumPy}, C{tuples} or other types of points.

Only L{HausdorffDistanceTo} -iff used with L{ellipsoidalKarney.LatLon}
points- and L{HausdorffKarney} requires installation of I{Charles Karney}'s
U{geographiclib<https://PyPI.org/project/geographiclib>}.

Typical usage is as follows.  First, create a C{Hausdorff} calculator
from a given set of C{LatLon} points, called the C{model} or C{template}
points.

C{h = HausdorffXyz(point1s, ...)}

Get the C{directed} or C{symmetric} Hausdorff distance to a second set
of C{LatLon} points, named the C{target} points, by using

C{t6 = h.directed(point2s)}

respectively

C{t6 = h.symmetric(point2s)}.

Or, use function C{hausdorff_} with a proper C{distance} function and
optionally a C{point} function passed as keyword arguments as follows

C{t6 = hausdorff_(point1s, point2s, ..., distance=..., point=...)}.

In all cases, the returned result C{t6} is a L{Hausdorff6Tuple}.

For C{(lat, lon, ...)} points in a C{NumPy} array or plain C{tuples},
wrap the points in a L{Numpy2LatLon} respectively L{Tuple2LatLon}
instance, more details in the documentation thereof.

For other points, create a L{Hausdorff} sub-class with the appropriate
C{distance} method overloading L{Hausdorff.distance} and optionally a
C{point} method overriding L{Hausdorff.point} as the next example.

    >>> from pygeodesy import Hausdorff, hypot_
    >>>
    >>> class H3D(Hausdorff):
    >>>     """Custom Hausdorff example.
    >>>     """
    >>>     def distance(self, p1, p2):
    >>>         return hypot_(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z)
    >>>
    >>> h3D = H3D(xyz1, ..., units="...")
    >>> d6 = h3D.directed(xyz2)

Transcribed from the original SciPy U{Directed Hausdorff Code
<https://GitHub.com/scipy/scipy/blob/master/scipy/spatial/_hausdorff.pyx>}
version 0.19.0, Copyright (C) Tyler Reddy, Richard Gowers, and Max Linke,
2016, distributed under the same BSD license as SciPy, including C{early
breaking} and C{random sampling} as in U{Abdel Aziz Taha, Allan Hanbury
"An Efficient Algorithm for Calculating the Exact Hausdorff Distance"
<https://Publik.TUWien.ac.AT/files/PubDat_247739.pdf>}, IEEE Trans. Pattern
Analysis Machine Intelligence (PAMI), vol 37, no 11, pp 2153-2163, Nov 2015.
'''

from pygeodesy.constants import INF, NINF, _0_0
from pygeodesy.datums import _ellipsoidal_datum, _WGS84
from pygeodesy.errors import _IsnotError, PointsError, _xattr, _xkwds, _xkwds_get
import pygeodesy.formy as _formy
from pygeodesy.interns import NN, _i_, _j_, _units_
# from pygeodesy.iters import points2  # from .points
from pygeodesy.lazily import _ALL_LAZY, _FOR_DOCS
from pygeodesy.named import _Named, _NamedTuple, notOverloaded, _Pass
# from pygeodesy.namedTuples import PhiLam2Tuple  # from .points
from pygeodesy.points import _distanceTo, points2 as _points2,  PhiLam2Tuple, radians
from pygeodesy.props import Property_RO, property_doc_, property_RO
from pygeodesy.units import Float, Number_, _xUnit, _xUnits
from pygeodesy.unitsBase import _Str_degrees, _Str_degrees2, _Str_meter, _Str_NN, \
                                _Str_radians, _Str_radians2

# from math import radians  # from .points
from random import Random

__all__ = _ALL_LAZY.hausdorff
__version__ = '23.08.06'


class HausdorffError(PointsError):
    '''Hausdorff issue.
    '''
    pass


class Hausdorff(_Named):
    '''Hausdorff base class, requires method L{Hausdorff.distance} to
       be overloaded.
    '''
    _datum = _WGS84
    _func  =  None  # formy function
    _kwds  = {}     # func_ options
    _model = ()
    _seed  =  None
    _units = _Str_NN  # XXX Str to _Pass and for backward compatibility

    def __init__(self, point1s, seed=None, name=NN, units=NN, **kwds):
        '''New C{Hausdorff...} calculator.

           @arg point1s: Initial set of points, aka the C{model} or
                         C{template} (C{LatLon}[], C{Numpy2LatLon}[],
                         C{Tuple2LatLon}[] or C{other}[]).
           @kwarg seed: Random sampling seed (C{any}) or C{None}, C{0}
                        or C{False} for no U{random sampling<https://
                        Publik.TUWien.ac.AT/files/PubDat_247739.pdf>}.
           @kwarg name: Optional name for this interpolator (C{str}).
           @kwarg units: Optional, the distance units (C{Unit} or C{str}).
           @kwarg kwds: Optional keyword argument for distance function,
                        retrievable with property C{kwds}.

           @raise HausdorffError: Insufficient number of B{C{point1s}}
                                  or an invalid B{C{point1}}, B{C{seed}}
                                  or B{C{units}}.
        '''
        _, self._model = self._points2(point1s)
        if seed:
            self.seed = seed
        if name:
            self.name = name
        if units:  # and not self.units:
            self.units = units
        if kwds:
            self._kwds = kwds

    @Property_RO
    def adjust(self):
        '''Get the adjust setting (C{bool} or C{None} if not applicable).
        '''
        return _xkwds_get(self._kwds, adjust=None)

    @Property_RO
    def datum(self):
        '''Get the datum of this calculator (L{Datum} or C{None} if not applicable).
        '''
        return self._datum

    def _datum_setter(self, datum):
        '''(INTERNAL) Set the datum.
        '''
        d = datum or _xattr(self._model[0], datum=datum)
        if d not in (None, self._datum):  # PYCHOK no cover
            self._datum = _ellipsoidal_datum(d, name=self.name)

    def directed(self, point2s, early=True):
        '''Compute only the C{forward Hausdorff} distance.

           @arg point2s: Second set of points, aka the C{target} (C{LatLon}[],
                         C{Numpy2LatLon}[], C{Tuple2LatLon}[] or C{other}[]).
           @kwarg early: Enable or disable U{early breaking<https://
                         Publik.TUWien.ac.AT/files/PubDat_247739.pdf>} (C{bool}).

           @return: A L{Hausdorff6Tuple}C{(hd, i, j, mn, md, units)}.

           @raise HausdorffError: Insufficient number of B{C{point2s}} or
                                  an invalid B{C{point2}}.

           @note: See B{C{point2s}} note at L{HausdorffDistanceTo}.
        '''
        return self._hausdorff_(point2s, False, early, self.distance)

    def distance(self, point1, point2):
        '''Return the distance between B{C{point1}} and B{C{point2s}},
           subject to the supplied optional keyword arguments, see
           property C{kwds}.
        '''
        return self._func(point1.lat, point1.lon,
                          point2.lat, point2.lon, **self._kwds)

    def _hausdorff_(self, point2s, both, early, distance):
        _, ps2 = self._points2(point2s)
        return _hausdorff_(self._model, ps2, both, early, self.seed,
                           self.units, distance, self.point)

    @property_RO
    def kwds(self):
        '''Get the supplied, optional keyword arguments (C{dict}).
        '''
        return self._kwds

    def point(self, point):
        '''Convert a C{model} or C{target} point for the C{.distance} method.
        '''
        return point  # pass thru

    def _points2(self, points):
        '''(INTERNAL) Check a set of points.
        '''
        return _points2(points, closed=False, Error=HausdorffError)

    @property_doc_(''' the random sampling seed (C{Random}).''')
    def seed(self):
        '''Get the random sampling seed (C{any} or C{None}).
        '''
        return self._seed

    @seed.setter  # PYCHOK setter!
    def seed(self, seed):
        '''Set the random sampling seed (C{Random(seed)}) or
           C{None}, C{0} or C{False} for no U{random sampling
           <https://Publik.TUWien.ac.AT/files/PubDat_247739.pdf>}.

           @raise HausdorffError: Invalid B{C{seed}}.
        '''
        if seed:
            try:
                Random(seed)
            except (TypeError, ValueError) as x:
                raise HausdorffError(seed=seed, cause=x)
            self._seed = seed
        else:
            self._seed = None

    def symmetric(self, point2s, early=True):
        '''Compute the combined C{forward and reverse Hausdorff} distance.

           @arg point2s: Second set of points, aka the C{target} (C{LatLon}[],
                         C{Numpy2LatLon}[], C{Tuple2LatLon}[] or C{other}[]).
           @kwarg early: Enable or disable U{early breaking<https://
                         Publik.TUWien.ac.AT/files/PubDat_247739.pdf>} (C{bool}).

           @return: A L{Hausdorff6Tuple}C{(hd, i, j, mn, md, units)}.

           @raise HausdorffError: Insufficient number of B{C{point2s}} or
                                  an invalid B{C{point2}}.

           @note: See B{C{point2s}} note at L{HausdorffDistanceTo}.
        '''
        return self._hausdorff_(point2s, True, early, self.distance)

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

    @Property_RO
    def wrap(self):
        '''Get the wrap setting (C{bool} or C{None} if not applicable).
        '''
        return _xkwds_get(self._kwds, adjust=None)


class HausdorffDegrees(Hausdorff):
    '''L{Hausdorff} base class for distances from C{LatLon}
       points in C{degrees}.
    '''
    _units = _Str_degrees

    if _FOR_DOCS:
        __init__  = Hausdorff.__init__
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric

    def distance(self, point1, point2):  # PYCHOK no cover
        '''I{Must be overloaded} to return the distance between
           B{C{point1}} and B{C{point2}} in C{degrees}.
        '''
        notOverloaded(self, point1, point2)


class HausdorffRadians(Hausdorff):
    '''L{Hausdorff} base class for distances from C{LatLon}
       points converted from C{degrees} to C{radians}.
    '''
    _units = _Str_radians

    if _FOR_DOCS:
        __init__  = Hausdorff.__init__
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric

    def distance(self, point1, point2):  # PYCHOK no cover
        '''I{Must be overloaded} to return the distance between
           B{C{point1}} and B{C{point2}} in C{radians}.
        '''
        notOverloaded(self, point1, point2)

    def point(self, point):
        '''Return B{C{point}} as L{PhiLam2Tuple} to maintain
           I{backward compatibility} of L{HausdorffRadians}.

           @return: A L{PhiLam2Tuple}C{(phi, lam)}.
        '''
        try:
            return point.philam
        except AttributeError:
            return PhiLam2Tuple(radians(point.lat), radians(point.lon))


class _HausdorffMeterRadians(Hausdorff):
    '''(INTERNAL) Returning C{meter} or C{radians} depending on
       the optional keyword arguments supplied at instantiation
       of the C{Hausdorff*} sub-class.
    '''
    _units  = _Str_meter
    _units_ = _Str_radians

    def directed(self, point2s, early=True):
        '''Overloaded method L{Hausdorff.directed} to determine
           the distance function and units from the optional
           keyword arguments given at this instantiation, see
           property C{kwds}.

           @see: L{Hausdorff.directed} for other details.
        '''
        return self._hausdorff_(point2s, False, early, _formy._radistance(self))

    def symmetric(self, point2s, early=True):
        '''Overloaded method L{Hausdorff.symmetric} to determine
           the distance function and units from the optional
           keyword arguments given at this instantiation, see
           property C{kwds}.

           @see: L{Hausdorff.symmetric} for other details.
        '''
        return self._hausdorff_(point2s, True, early, _formy._radistance(self))

    def _func_(self, *args, **kwds):  # PYCHOK no cover
        notOverloaded(self, *args, **kwds)


class HausdorffCosineAndoyerLambert(_HausdorffMeterRadians):
    '''Compute the C{Hausdorff} distance based on the I{angular} distance
       in C{radians} from function L{pygeodesy.cosineAndoyerLambert}.
    '''
    def __init__(self, point1s, seed=None, name=NN, **datum_wrap):
        '''New L{HausdorffCosineAndoyerLambert} calculator.

          @see: L{Hausdorff.__init__} for details about B{C{point1s}},
                 B{C{seed}}, B{C{name}} and other exceptions.

           @kwarg datum_wrap: Optional keyword arguments for function
                              L{pygeodesy.cosineAndoyerLambert}.
        '''
        Hausdorff.__init__(self, point1s, seed=seed, name=name,
                                        **datum_wrap)
        self._func  = _formy.cosineAndoyerLambert
        self._func_ = _formy.cosineAndoyerLambert_

    if _FOR_DOCS:
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric


class HausdorffCosineForsytheAndoyerLambert(_HausdorffMeterRadians):
    '''Compute the C{Hausdorff} distance based on the I{angular} distance
       in C{radians} from function L{pygeodesy.cosineForsytheAndoyerLambert}.
    '''
    def __init__(self, point1s, seed=None, name=NN, **datum_wrap):
        '''New L{HausdorffCosineForsytheAndoyerLambert} calculator.

           @see: L{Hausdorff.__init__} for details about B{C{point1s}},
                 B{C{seed}}, B{C{name}} and other exceptions.

           @kwarg datum_wrap: Optional keyword arguments for function
                              L{pygeodesy.cosineAndoyerLambert}.
        '''
        Hausdorff.__init__(self, point1s, seed=seed, name=name,
                                        **datum_wrap)
        self._func  = _formy.cosineForsytheAndoyerLambert
        self._func_ = _formy.cosineForsytheAndoyerLambert_

    if _FOR_DOCS:
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric


class HausdorffCosineLaw(_HausdorffMeterRadians):
    '''Compute the C{Hausdorff} distance based on the I{angular}
       distance in C{radians} from function L{pygeodesy.cosineLaw_}.

       @note: See note at function L{pygeodesy.vincentys_}.
    '''
    def __init__(self, point1s, seed=None, name=NN, **radius_wrap):
        '''New L{HausdorffCosineLaw} calculator.

           @kwarg radius_wrap: Optional keyword arguments for function
                               L{pygeodesy.cosineLaw}.

           @see: L{Hausdorff.__init__} for details about B{C{point1s}},
                 B{C{seed}}, B{C{name}} and other exceptions.
        '''
        Hausdorff.__init__(self, point1s, seed=seed, name=name,
                                        **radius_wrap)
        self._func  = _formy.cosineLaw
        self._func_ = _formy.cosineLaw_

    if _FOR_DOCS:
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric


class HausdorffDistanceTo(Hausdorff):
    '''Compute the C{Hausdorff} distance based on the distance from the
       points' C{LatLon.distanceTo} method, conventionally in C{meter}.
    '''
    _units = _Str_meter

    def __init__(self, point1s, seed=None, name=NN, **distanceTo_kwds):
        '''New L{HausdorffDistanceTo} calculator.

           @kwarg distanceTo_kwds: Optional keyword arguments for each
                                   B{C{point1s}}' C{LatLon.distanceTo}
                                   method.

           @see: L{Hausdorff.__init__} for details about B{C{point1s}},
                 B{C{seed}}, B{C{name}} and other exceptions.

           @note: All C{model}, C{template} and C{target} B{C{points}}
                  I{must} be instances of the same ellipsoidal or
                  spherical C{LatLon} class.
        '''
        Hausdorff.__init__(self, point1s, seed=seed, name=name,
                                        **distanceTo_kwds)

    if _FOR_DOCS:
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric

    def distance(self, p1, p2):
        '''Return the distance in C{meter}.
        '''
        return p1.distanceTo(p2, **self._kwds)

    def _points2(self, points):
        '''(INTERNAL) Check a set of points.
        '''
        np, ps = Hausdorff._points2(self, points)
        return np, _distanceTo(HausdorffError, points=ps)


class HausdorffEquirectangular(Hausdorff):
    '''Compute the C{Hausdorff} distance based on the C{equirectangular} distance
       in C{radians squared} like function L{pygeodesy.equirectangular_}.
    '''
    _units = _Str_degrees2

    def __init__(self, point1s, seed=None, name=NN, **adjust_limit_wrap):
        '''New L{HausdorffEquirectangular} calculator.

           @kwarg adjust_limit_wrap: Optional keyword arguments for function
                               L{pygeodesy.equirectangular_} I{with default}
                               C{B{limit}=0} for I{backward compatibility}.

           @see: L{Hausdorff.__init__} for details about B{C{point1s}},
                 B{C{seed}}, B{C{name}} and other exceptions.
        '''
        adjust_limit_wrap = _xkwds(adjust_limit_wrap, limit=0)
        Hausdorff.__init__(self, point1s, seed=seed, name=name,
                                        **adjust_limit_wrap)
        self._func = _formy._equirectangular  # helper

    if _FOR_DOCS:
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric


class HausdorffEuclidean(_HausdorffMeterRadians):
    '''Compute the C{Hausdorff} distance based on the C{Euclidean}
       distance in C{radians} from function L{pygeodesy.euclidean_}.
    '''
    def __init__(self, point1s, seed=None, name=NN, **adjust_wrap):
        '''New L{HausdorffEuclidean} calculator.

           @kwarg adjust_radius_wrap: Optional keyword arguments for
                                function L{pygeodesy.euclidean}.

           @see: L{Hausdorff.__init__} for details about B{C{point1s}},
                 B{C{seed}}, B{C{name}} and other exceptions.
        '''
        Hausdorff.__init__(self, point1s, seed=seed, name=name,
                                        **adjust_wrap)
        self._func  = _formy.euclidean
        self._func_ = _formy.euclidean_

    if _FOR_DOCS:
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric


class HausdorffExact(Hausdorff):
    '''Compute the C{Hausdorff} distance based on the I{angular}
       distance in C{degrees} from method L{GeodesicExact}C{.Inverse}.
    '''
    _units = _Str_degrees

    def __init__(self, point1s, seed=None, name=NN, datum=None, **wrap):
        '''New L{HausdorffKarney} calculator.

           @kwarg datum: Datum to override the default C{Datums.WGS84} and
                         first B{C{point1s}}' datum (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg wrap: Optional keyword argument for method C{Inverse1}
                        of class L{geodesicx.GeodesicExact}.

           @see: L{Hausdorff.__init__} for details about B{C{point1s}},
                 B{C{seed}}, B{C{name}} and other exceptions.

           @raise TypeError: Invalid B{C{datum}}.
        '''
        Hausdorff.__init__(self, point1s, seed=seed, name=name,
                                                   **wrap)
        self._datum_setter(datum)
        self._func = self.datum.ellipsoid.geodesicx.Inverse1  # note -x

    if _FOR_DOCS:
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric


class HausdorffFlatLocal(_HausdorffMeterRadians):
    '''Compute the C{Hausdorff} distance based on the I{angular} distance in
       C{radians squared} like function L{pygeodesy.flatLocal_}/L{pygeodesy.hubeny_}.
    '''
    _units = _Str_radians2

    def __init__(self, point1s, seed=None, name=NN, **datum_scaled_wrap):
        '''New L{HausdorffFlatLocal}/L{HausdorffHubeny} calculator.

           @kwarg datum_scaled_wrap: Optional keyword arguments for
                                     function L{pygeodesy.flatLocal}.

           @see: L{Hausdorff.__init__} for details about B{C{point1s}},
                 B{C{seed}}, B{C{name}} and other exceptions.

           @note: The distance C{units} are C{radians squared}, not C{radians}.
        '''
        Hausdorff.__init__(self, point1s, seed=seed, name=name,
                                        **datum_scaled_wrap)
        self._func  = _formy.flatLocal
        self._func_ =  self.datum.ellipsoid._hubeny_2

    if _FOR_DOCS:
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric


class HausdorffFlatPolar(_HausdorffMeterRadians):
    '''Compute the C{Hausdorff} distance based on the I{angular}
       distance in C{radians} from function L{pygeodesy.flatPolar_}.
    '''
    _wrap = False

    def __init__(self, points, seed=None, name=NN, **radius_wrap):
        '''New L{HausdorffFlatPolar} calculator.

           @kwarg radius_wrap: Optional keyword arguments for function
                               L{pygeodesy.flatPolar}.

           @see: L{Hausdorff.__init__} for details about B{C{point1s}},
                 B{C{seed}}, B{C{name}} and other exceptions.
        '''
        Hausdorff.__init__(self, points, seed=seed, name=name,
                                       **radius_wrap)
        self._func  = _formy.flatPolar
        self._func_ = _formy.flatPolar_

    if _FOR_DOCS:
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric


class HausdorffHaversine(_HausdorffMeterRadians):
    '''Compute the C{Hausdorff} distance based on the I{angular}
       distance in C{radians} from function L{pygeodesy.haversine_}.

       @note: See note under L{HausdorffVincentys}.
    '''
    _wrap = False

    def __init__(self, points, seed=None, name=NN, **radius_wrap):
        '''New L{HausdorffHaversine} calculator.

           @kwarg radius_wrap: Optional keyword arguments for function
                               L{pygeodesy.haversine}.

           @see: L{Hausdorff.__init__} for details about B{C{point1s}},
                 B{C{seed}}, B{C{name}} and other exceptions.
       '''
        Hausdorff.__init__(self, points, seed=seed, name=name,
                                       **radius_wrap)
        self._func  = _formy.haversine
        self._func_ = _formy.haversine_

    if _FOR_DOCS:
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric


class HausdorffHubeny(HausdorffFlatLocal):  # for Karl Hubeny
    if _FOR_DOCS:
        __doc__   = HausdorffFlatLocal.__doc__
        __init__  = HausdorffFlatLocal.__init__
        directed  = HausdorffFlatLocal.directed
        distance  = HausdorffFlatLocal.distance
        symmetric = HausdorffFlatLocal.symmetric


class HausdorffKarney(Hausdorff):
    '''Compute the C{Hausdorff} distance based on the I{angular}
       distance in C{degrees} from I{Karney}'s U{geographiclib
       <https://PyPI.org/project/geographiclib>} U{Geodesic
       <https://GeographicLib.SourceForge.io/Python/doc/code.html>}
       Inverse method.
    '''
    _units = _Str_degrees

    def __init__(self, point1s, datum=None, seed=None, name=NN, **wrap):
        '''New L{HausdorffKarney} calculator.

           @kwarg datum: Datum to override the default C{Datums.WGS84} and
                         first B{C{knots}}' datum (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg wrap: Optional keyword argument for method C{Inverse1}
                        of class L{geodesicw.Geodesic}.

           @raise ImportError: Package U{geographiclib
                  <https://PyPI.org/project/geographiclib>} missing.

           @raise TypeError: Invalid B{C{datum}}.

           @see: L{Hausdorff.__init__} for details about B{C{point1s}},
                 B{C{seed}}, B{C{name}} and other exceptions.
        '''
        Hausdorff.__init__(self, point1s, seed=seed, name=name,
                                        **wrap)
        self._datum_setter(datum)
        self._func = self.datum.ellipsoid.geodesic.Inverse1


class HausdorffThomas(_HausdorffMeterRadians):
    '''Compute the C{Hausdorff} distance based on the I{angular}
       distance in C{radians} from function L{pygeodesy.thomas_}.
    '''
    def __init__(self, point1s, seed=None, name=NN, **datum_wrap):
        '''New L{HausdorffThomas} calculator.

           @kwarg datum_wrap: Optional keyword argument for function
                              L{pygeodesy.thomas}.

           @see: L{Hausdorff.__init__} for details about B{C{point1s}},
                 B{C{seed}}, B{C{name}} and other exceptions.
        '''
        Hausdorff.__init__(self, point1s, seed=seed, name=name,
                                        **datum_wrap)
        self._func  = _formy.thomas
        self._func_ = _formy.thomas_

    if _FOR_DOCS:
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric


class HausdorffVincentys(_HausdorffMeterRadians):
    '''Compute the C{Hausdorff} distance based on the I{angular}
       distance in C{radians} from function L{pygeodesy.vincentys_}.

       @note: See note at function L{pygeodesy.vincentys_}.
    '''
    _wrap = False

    def __init__(self, point1s, seed=None, name=NN, **radius_wrap):
        '''New L{HausdorffVincentys} calculator.

           @kwarg radius_wrap: Optional keyword arguments for function
                               L{pygeodesy.vincentys}.

           @see: L{Hausdorff.__init__} for details about B{C{point1s}},
                 B{C{seed}}, B{C{name}} and other exceptions.
        '''
        Hausdorff.__init__(self, point1s, seed=seed, name=name,
                                        **radius_wrap)
        self._func  = _formy.vincentys
        self._func_ = _formy.vincentys_

    if _FOR_DOCS:
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric


def _hausdorff_(ps1, ps2, both, early, seed, units, distance, point):
    '''(INTERNAL) Core of function L{hausdorff_} and methods C{directed}
       and C{symmetric} of classes C{hausdorff.Hausdorff...}.
    '''
    # shuffling the points generally increases the
    # chance of an early break in the inner j loop
    rr = randomrangenerator(seed) if seed else range

    hd = NINF
    hi = hj = m = mn = 0
    md = _0_0

    # forward or forward and backward
    for fb in range(2 if both else 1):
        n = len(ps2)
        for i in rr(len(ps1)):
            p1 = point(ps1[i])
            dh, dj = INF, 0
            for j in rr(n):
                p2 = point(ps2[j])
                d = distance(p1, p2)
                if early and d < hd:
                    break  # early
                elif d < dh:
                    dh, dj = d, j
            else:  # no early break
                if hd < dh:
                    hd = dh
                    if fb:
                        hi, hj = dj, i
                    else:
                        hi, hj = i, dj
                md += dh
                mn += 1
            m += 1
        # swap model and target
        ps1, ps2 = ps2, ps1

    md = None if mn < m else (md / float(m))
    return Hausdorff6Tuple(hd, hi, hj, m, md, units)


def _point(p):
    '''Default B{C{point}} callable for function L{hausdorff_}.

       @arg p: The original C{model} or C{target} point (C{any}).

       @return: The point, suitable for the L{hausdorff_}
                B{C{distance}} callable.
    '''
    return p


def hausdorff_(model, target, both=False, early=True, seed=None, units=NN,
                              distance=None, point=_point):
    '''Compute the C{directed} or C{symmetric} U{Hausdorff
       <https://WikiPedia.org/wiki/Hausdorff_distance>} distance between 2 sets of points
       with or without U{early breaking<https://Publik.TUWien.ac.AT/files/PubDat_247739.pdf>}
       and U{random sampling<https://Publik.TUWien.ac.AT/files/PubDat_247739.pdf>}.

       @arg model: First set of points (C{LatLon}[], C{Numpy2LatLon}[],
                   C{Tuple2LatLon}[] or C{other}[]).
       @arg target: Second set of points (C{LatLon}[], C{Numpy2LatLon}[],
                    C{Tuple2LatLon}[] or C{other}[]).
       @kwarg both: Return the C{directed} (forward only) or the C{symmetric}
                    (combined forward and reverse) C{Hausdorff} distance (C{bool}).
       @kwarg early: Enable or disable U{early breaking<https://Publik.TUWien.ac.AT/
                     files/PubDat_247739.pdf>} (C{bool}).
       @kwarg seed: Random sampling seed (C{any}) or C{None}, C{0} or C{False} for no
                    U{random sampling<https://Publik.TUWien.ac.AT/files/PubDat_247739.pdf>}.
       @kwarg units: Optional, the distance units (C{Unit} or C{str}).
       @kwarg distance: Callable returning the distance between a B{C{model}}
                        and B{C{target}} point (signature C{(point1, point2)}).
       @kwarg point: Callable returning the B{C{model}} or B{C{target}} point
                     suitable for B{C{distance}} (signature C{(point)}).

       @return: A L{Hausdorff6Tuple}C{(hd, i, j, mn, md, units)}.

       @raise HausdorffError: Insufficient number of B{C{model}} or B{C{target}} points.

       @raise TypeError: If B{C{distance}} or B{C{point}} is not callable.
    '''
    if not callable(distance):
        raise _IsnotError(callable.__name__, distance=distance)
    if not callable(point):
        raise _IsnotError(callable.__name__, point=point)

    _, ps1 = _points2(model,  closed=False, Error=HausdorffError)  # PYCHOK non-sequence
    _, ps2 = _points2(target, closed=False, Error=HausdorffError)  # PYCHOK non-sequence
    return _hausdorff_(ps1, ps2, both, early, seed, units, distance, point)


class Hausdorff6Tuple(_NamedTuple):
    '''6-Tuple C{(hd, i, j, mn, md, units)} with the U{Hausdorff
       <https://WikiPedia.org/wiki/Hausdorff_distance>} distance C{hd},
       indices C{i} and C{j}, the total count C{mn}, the C{I{mean}
       Hausdorff} distance C{md} and the class or name of both distance
       C{units}.

       For C{directed Hausdorff} distances, count C{mn} is the number
       of model points considered. For C{symmetric Hausdorff} distances
       count C{mn} twice that.

       Indices C{i} and C{j} are the C{model} respectively C{target}
       point with the C{hd} distance.

       Mean distance C{md} is C{None} if an C{early break} occurred and
       U{early breaking<https://Publik.TUWien.ac.AT/files/PubDat_247739.pdf>}
       was enabled by keyword argument C{early=True}.
    '''
    _Names_ = ('hd', _i_,     _j_,     'mn',     'md',  _units_)
    _Units_ = (_Pass, Number_, Number_, Number_, _Pass, _Pass)

    def toUnits(self, **Error):  # PYCHOK expected
        '''Overloaded C{_NamedTuple.toUnits} for C{hd} and C{md} units.
        '''
        U = _xUnit(self.units, Float)  # PYCHOK expected
        M = _Pass if self.md is None else U  # PYCHOK expected
        self._Units_ = (U,) + Hausdorff6Tuple._Units_[1:4] \
                     + (M,) + Hausdorff6Tuple._Units_[5:]
        return _NamedTuple.toUnits(self, **Error)


def randomrangenerator(seed):
    '''Return a C{seed}ed random range function generator.

       @arg seed: Initial, internal L{Random} state (C{hashable}
                  or C{None}).

       @note: L{Random} with C{B{seed} is None} seeds from the
              current time or from a platform-specific randomness
              source, if available.

       @return: A function to generate random ranges.

       @example:

        >>> rrange = randomrangenerator('R')
        >>> for r in rrange(n):
        >>>    ...  # r is random in 0..n-1
    '''
    R = Random(seed)

    def _range(n, *stop_step):
        '''Like standard L{range}C{start, stop=..., step=...)},
           except the returned values are in random order.

           @note: Especially C{range(n)} behaves like standard
                  L{Random.sample}C{(range(n), n)} but avoids
                  creating a tuple with the entire C{population}
                  and a list containing all sample values (for
                  large C{n}).
        '''
        if stop_step:
            s = range(n, *stop_step)

        elif n > 32:
            r = R.randrange  # Random._randbelow
            s = set()
            for _ in range(n - 32):
                i = r(n)
                while i in s:
                    i = r(n)
                s.add(i)
                yield i
            s = set(range(n)) - s  # [i for i in range(n) if i not in s]
        else:
            s = range(n)

        s = list(s)
        R.shuffle(s)
        while s:
            yield s.pop(0)

    return _range

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
