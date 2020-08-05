
# -*- coding: utf-8 -*-

u'''Classes L{Hausdorff}, L{HausdorffDegrees}, L{HausdorffRadians},
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

C{h = HausdorffXyz(points1, ...)}

Get the C{directed} or C{symmetric} Hausdorff distance to a second set
of C{LatLon} points, named the C{target} points, by using

C{t6 = h.directed(points2)}

respectively

C{t6 = h.symmetric(points2)}.

Or, use function C{hausdorff_} with a proper C{distance} function and
optionally a C{point} function passed as keyword arguments as follows

C{t6 = hausdorff_(points1, points2, ..., distance=..., point=...)}.

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

from pygeodesy.basics import _bkwds, INF, property_doc_, property_RO, \
                             _xinstanceof
from pygeodesy.datum import Datums, Datum
from pygeodesy.errors import _IsnotError, PointsError
from pygeodesy.fmath import hypot2
from pygeodesy.formy import cosineAndoyerLambert_, cosineForsytheAndoyerLambert_, \
                            cosineLaw_, euclidean_, flatPolar_, haversine_, \
                            points2, _scale_rad, thomas_, vincentys_
from pygeodesy.interns import _datum_, _degrees_, _distanceTo_, _item_sq, _meter_, \
                               NN, _points_, _radians_, _radians2_, _units_
from pygeodesy.lazily import _ALL_LAZY, _ALL_DOCS, _FOR_DOCS
from pygeodesy.named import _Named, _NamedTuple, notOverloaded, PhiLam2Tuple
from pygeodesy.utily import unrollPI

from math import radians
from random import Random

__all__ = _ALL_LAZY.hausdorff
__version__ = '20.08.04'


class HausdorffError(PointsError):
    '''Hausdorff issue.
    '''
    pass


class Hausdorff6Tuple(_NamedTuple):
    '''6-Tuple C{(hd, i, j, mn, md, units)} with the U{Hausdorff
       <https://WikiPedia.org/wiki/Hausdorff_distance>} distance C{hd},
       indices C{i} and C{j}, the total count C{mn}, the C{I{mean}
       Hausdorff} distance C{md} and the name of the distance C{units}.

       For C{directed Hausdorff} distances, count C{mn} is the number
       of model points considered. For C{symmetric Hausdorff} distances
       count C{mn} twice that.

       Indices C{i} and C{j} are the C{model} respectively C{target}
       point with the C{hd} distance.

       Mean distance C{md} is C{None} if an C{early break} occurred and
       U{early breaking<https://Publik.TUWien.ac.AT/files/PubDat_247739.pdf>}
       was enabled by keyword argument C{early=True}.
    '''
    _Names_ = ('hd', 'i', 'j', 'mn', 'md', _units_)


class Hausdorff(_Named):
    '''Hausdorff base class, requires method L{Hausdorff.distance} to
       be overloaded.
    '''
    _adjust = None  # not applicable
    _datum  = None  # not applicable
    _model  = ()
    _seed   = None
    _units  = NN
    _wrap   = None  # not applicable

    def __init__(self, points, seed=None, name=NN, units=NN, **wrap_adjust):
        '''New C{Hausdorff...} calculator.

           @arg points: Initial set of points, aka the C{model} or
                        C{template} (C{LatLon}[], C{Numpy2LatLon}[],
                        C{Tuple2LatLon}[] or C{other}[]).
           @kwarg seed: Random sampling seed (C{any}) or C{None}, C{0}
                        or C{False} for no U{random sampling<https://
                        Publik.TUWien.ac.AT/files/PubDat_247739.pdf>}.
           @kwarg name: Optional name for this interpolator (C{str}).
           @kwarg units: Optional, distance units (C{str}).
           @kwarg wrap_adjust: Optionally, C{wrap} and unroll longitudes, iff
                               applicable (C{bool}) and C{adjust} wrapped,
                               unrolled longitudinal delta by the cosine
                               of the mean latitude, iff applicable.

           @raise HausdorffError: Insufficient number of B{C{points}} or an
                                  invalid B{C{point}}, B{C{seed}} or B{{wrap}}
                                  or B{C{ajust}} not applicable.
        '''
        _, self._model = self._points2(points)
        if seed:
            self.seed = seed
        if name:
            self.name = name
        if units and not self.units:
            self.units = units
        if wrap_adjust:
            _bkwds(self, Error=HausdorffError, **wrap_adjust)

    @property_RO
    def adjust(self):
        '''Get the adjust setting (C{bool} or C{None} if not applicable).
        '''
        return self._adjust

    @property_RO
    def datum(self):
        '''Get the datum of this calculator (L{Datum} or C{None} if not applicable).
        '''
        return self._datum

    def _datum_setter(self, datum):
        '''(INTERNAL) Set the datum.
        '''
        d = datum or getattr(self._model[0], _datum_, datum)
        if d and d != self.datum:  # PYCHOK no cover
            _xinstanceof(Datum, datum=d)
            self._datum = d

    def directed(self, points, early=True):
        '''Compute only the C{forward Hausdorff} distance.

           @arg points: Second set of points, aka the C{target} (C{LatLon}[],
                        C{Numpy2LatLon}[], C{Tuple2LatLon}[] or C{other}[]).
           @kwarg early: Enable or disable U{early breaking<https://
                         Publik.TUWien.ac.AT/files/PubDat_247739.pdf>} (C{bool}).

           @return: A L{Hausdorff6Tuple}C{(hd, i, j, mn, md, units)}.

           @raise HausdorffError: Insufficient number of B{C{points}} or
                                  an invalid B{C{point}}.

           @note: See B{C{points}} note at L{HausdorffDistanceTo}.
        '''
        _, ps2 = self._points2(points)
        return _hausdorff_(self._model, ps2, False, early, self.seed,
                           self.units, self.distance, self.point)

    def distance(self, point1, point2):
        '''(INTERNAL) I{Must be overloaded}.
        '''
        notOverloaded(self, self.distance, point1, point2)  # PYCHOK no cover

    def point(self, point):
        '''Convert a C{model} or C{target} point for the C{.distance} method.
        '''
        return point  # pass thru

    def _points2(self, points):
        '''(INTERNAL) Check a set of points.
        '''
        return points2(points, closed=False, Error=HausdorffError)

    @property_doc_(''' the random sampling seed (C{Random}).''')
    def seed(self):
        '''Get the random sampling seed (C{any} or C{None}).
        '''
        return self._seed

    @seed.setter  # PYCHOK setter!
    def seed(self, seed):
        '''Set the random sampling seed.

           @arg seed: Valid L{Random(seed)} or C{None}, C{0} or
                      C{False} for no U{random sampling<https://
                      Publik.TUWien.ac.AT/files/PubDat_247739.pdf>}.

           @raise HausdorffError: Invalid B{C{seed}}.
        '''
        if seed:
            try:
                Random(seed)
            except (TypeError, ValueError) as x:
                raise HausdorffError(seed=seed, txt=str(x))
            self._seed = seed
        else:
            self._seed = None

    def symmetric(self, points, early=True):
        '''Compute the combined C{forward and reverse Hausdorff} distance.

           @arg points: Second set of points, aka the C{target} (C{LatLon}[],
                        C{Numpy2LatLon}[], C{Tuple2LatLon}[] or C{other}[]).
           @kwarg early: Enable or disable U{early breaking<https://
                         Publik.TUWien.ac.AT/files/PubDat_247739.pdf>} (C{bool}).

           @return: A L{Hausdorff6Tuple}C{(hd, i, j, mn, md, units)}.

           @raise HausdorffError: Insufficient number of B{C{points}} or
                                  an invalid B{C{point}}.

           @note: See B{C{points}} note at L{HausdorffDistanceTo}.
        '''
        _, ps2 = self._points2(points)
        return _hausdorff_(self._model, ps2, True, early, self.seed,
                           self.units, self.distance, self.point)

    @property_doc_(''' the distance units (C{str}).''')
    def units(self):
        '''Get the distance units (C{str} or C{""}).
        '''
        return self._units

    @units.setter  # PYCHOK setter!
    def units(self, units):
        '''Set the distance units.

           @arg units: New units name (C{str}).
        '''
        self._units = str(units or NN)

    @property_RO
    def wrap(self):
        '''Get the wrap setting (C{bool} or C{None} if not applicable).
        '''
        return self._wrap


class HausdorffDegrees(Hausdorff):
    '''L{Hausdorff} base class for distances from C{LatLon}
       points in C{degrees}.
    '''
    _units = _degrees_

    if _FOR_DOCS:
        __init__  = Hausdorff.__init__
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric


class HausdorffRadians(Hausdorff):
    '''L{Hausdorff} base class for distances from C{LatLon}
       points converted from C{degrees} to C{radians}.
    '''
    _units = _radians_

    if _FOR_DOCS:
        __init__  = Hausdorff.__init__
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric

    def point(self, point):
        '''Convert C{(lat, lon)} point in degrees to C{(a, b)}
           in radians.

           @return: An L{PhiLam2Tuple}C{(phi, lam)}.
        '''
        try:
            return point.philam
        except AttributeError:  # PYCHOK no cover
            return PhiLam2Tuple(radians(point.lat), radians(point.lon))


class HausdorffCosineAndoyerLambert(HausdorffRadians):
    '''Compute the C{Hausdorff} distance based on the I{angular}
       distance in C{radians} from function L{cosineAndoyerLambert_}.

       @see: L{HausdorffCosineForsytheAndoyerLambert}, L{HausdorffDistanceTo},
             L{HausdorffFlatLocal}, L{HausdorffHubeny},
             L{HausdorffThomas} and L{HausdorffKarney}.
    '''
    _datum = Datums.WGS84
    _wrap  = False

    def __init__(self, points, datum=None, wrap=False, seed=None, name=NN):
        '''New L{HausdorffCosineAndoyerLambert} calculator.

           @arg points: Initial set of points, aka the C{model} or
                        C{template} (C{LatLon}[], C{Numpy2LatLon}[],
                        C{Tuple2LatLon}[] or C{other}[]).
           @kwarg datum: Optional datum overriding the default C{Datums.WGS84}
                         and first B{C{points}}' datum (L{Datum}).
           @kwarg wrap: Optionally, wrap and L{unrollPI} longitudes (C{bool}).
           @kwarg seed: Random seed (C{any}) or C{None}, C{0} or
                        C{False} for no U{random sampling<https://
                        Publik.TUWien.ac.AT/files/PubDat_247739.pdf>}.
           @kwarg name: Optional name for this interpolator (C{str}).

           @raise HausdorffError: Insufficient number of B{C{points}} or
                                  invalid B{C{seed}}.
        '''
        HausdorffRadians.__init__(self, points, seed=seed, name=name,
                                                           wrap=wrap)
        self._datum_setter(datum)

    if _FOR_DOCS:
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric

    def distance(self, p1, p2):
        '''Return the L{cosineAndoyerLambert_} distance in C{radians}.
        '''
        r, _ = unrollPI(p1.lam, p2.lam, wrap=self._wrap)
        return cosineAndoyerLambert_(p2.phi, p1.phi, r, datum=self._datum)


class HausdorffCosineForsytheAndoyerLambert(HausdorffRadians):
    '''Compute the C{Hausdorff} distance based on the I{angular}
       distance in C{radians} from function L{cosineForsytheAndoyerLambert_}.

       @see: L{HausdorffCosineAndoyerLambert}, L{HausdorffDistanceTo},
             L{HausdorffFlatLocal}, L{HausdorffHubeny},
             L{HausdorffThomas} and L{HausdorffKarney}.
    '''
    _datum = Datums.WGS84
    _wrap  = False

    def __init__(self, points, datum=None, wrap=False, seed=None, name=NN):
        '''New L{HausdorffCosineForsytheAndoyerLambert} calculator.

           @arg points: Initial set of points, aka the C{model} or
                        C{template} (C{LatLon}[], C{Numpy2LatLon}[],
                        C{Tuple2LatLon}[] or C{other}[]).
           @kwarg datum: Optional datum overriding the default C{Datums.WGS84}
                         and first B{C{points}}' datum (L{Datum}).
           @kwarg wrap: Optionally, wrap and L{unrollPI} longitudes (C{bool}).
           @kwarg seed: Random seed (C{any}) or C{None}, C{0} or
                        C{False} for no U{random sampling<https://
                        Publik.TUWien.ac.AT/files/PubDat_247739.pdf>}.
           @kwarg name: Optional name for this interpolator (C{str}).

           @raise HausdorffError: Insufficient number of B{C{points}} or
                                  invalid B{C{seed}}.
        '''
        HausdorffRadians.__init__(self, points, seed=seed, name=name,
                                                           wrap=wrap)
        self._datum_setter(datum)

    if _FOR_DOCS:
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric

    def distance(self, p1, p2):
        '''Return the L{cosineForsytheAndoyerLambert_} distance in C{radians}.
        '''
        r, _ = unrollPI(p1.lam, p2.lam, wrap=self._wrap)
        return cosineForsytheAndoyerLambert_(p2.phi, p1.phi, r, datum=self._datum)


class HausdorffCosineLaw(HausdorffRadians):
    '''Compute the C{Hausdorff} distance based on the I{angular}
       distance in C{radians} from function L{cosineLaw_}.

       @note: See note at function L{vincentys_}.

       @see: L{HausdorffEquirectangular}, L{HausdorffEuclidean},
             L{HausdorffFlatLocal}, L{HausdorffHubeny}, L{HausdorffFlatPolar},
             L{HausdorffHaversine}, L{HausdorffKarney} and
             L{HausdorffVincentys}.
    '''
    _wrap = False

    def __init__(self, points, wrap=False, seed=None, name=NN):
        '''New L{HausdorffCosineLaw} calculator.

           @arg points: Initial set of points, aka the C{model} or
                        C{template} (C{LatLon}[], C{Numpy2LatLon}[],
                        C{Tuple2LatLon}[] or C{other}[]).
           @kwarg wrap: Optionally, wrap and L{unrollPI} longitudes (C{bool}).
           @kwarg seed: Random seed (C{any}) or C{None}, C{0} or
                        C{False} for no U{random sampling<https://
                        Publik.TUWien.ac.AT/files/PubDat_247739.pdf>}.
           @kwarg name: Optional name for this interpolator (C{str}).

           @raise HausdorffError: Insufficient number of B{C{points}} or
                                  invalid B{C{seed}}.
        '''
        HausdorffRadians.__init__(self, points, seed=seed, name=name,
                                                           wrap=wrap)

    if _FOR_DOCS:
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric

    def distance(self, p1, p2):
        '''Return the L{cosineLaw_} distance in C{radians}.
        '''
        d, _ = unrollPI(p1.lam, p2.lam, wrap=self._wrap)
        return cosineLaw_(p2.phi, p1.phi, d)


class HausdorffDistanceTo(Hausdorff):
    '''Compute the C{Hausdorff} distance based on the distance from the
       points' C{LatLon.distanceTo} method, conventionally in C{meter}.

       @see: L{HausdorffCosineAndoyerLambert},
             L{HausdorffCosineForsytheAndoyerLambert},
             L{HausdorffFlatLocal}, L{HausdorffHubeny},
             L{HausdorffThomas} and L{HausdorffKarney}.
    '''
    _distanceTo_kwds =  {}
    _units           = _meter_

    def __init__(self, points, seed=None, name=NN, **distanceTo_kwds):
        '''New L{HausdorffDistanceTo} calculator.

           @arg points: Initial set of points, aka the C{model} or
                        C{template} (C{LatLon}[]).
           @kwarg seed: Random sampling seed (C{any}) or C{None}, C{0}
                        or C{False} for no U{random sampling<https://
                        Publik.TUWien.ac.AT/files/PubDat_247739.pdf>}.
           @kwarg name: Optional name for this interpolator (C{str}).
           @kwarg distanceTo_kwds: Optional keyword arguments for the
                                   B{C{points}}' C{LatLon.distanceTo}
                                   method.

           @raise HausdorffError: Insufficient number of B{C{points}} or
                                  an invalid B{C{point}} or B{C{seed}}.

           @raise ImportError: Package U{GeographicLib
                  <https://PyPI.org/project/geographiclib>} missing
                  iff B{C{points}} are L{ellipsoidalKarney.LatLon}s.

           @note: All C{model}, C{template} and C{target} B{C{points}}
                  I{must} be instances of the same ellipsoidal or
                  spherical C{LatLon} class.
        '''
        Hausdorff.__init__(self, points, seed=seed, name=name)
        if distanceTo_kwds:
            self._distanceTo_kwds = distanceTo_kwds

    if _FOR_DOCS:
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric

    def distance(self, p1, p2):
        '''Return the distance in C{meter}.
        '''
        return p1.distanceTo(p2, **self._distanceTo_kwds)

    def _points2(self, points):
        '''(INTERNAL) Check a set of points.
        '''
        np, ps = Hausdorff._points2(self, points)
        for i, p in enumerate(ps):
            if not callable(getattr(p, _distanceTo_, None)):
                raise HausdorffError(_item_sq(_points_, i), p, txt=_distanceTo_)
        return np, ps


class HausdorffEquirectangular(HausdorffRadians):
    '''Compute the C{Hausdorff} distance based on the C{equirectangular}
       distance in C{radians squared} like function L{equirectangular_}.

       @see: L{HausdorffCosineLaw}, L{HausdorffEuclidean}
             L{HausdorffFlatPolar}, L{HausdorffHaversine}
             and L{HausdorffVincentys}.
    '''
    _adjust =  True
    _units  = _radians2_
    _wrap   =  False

    def __init__(self, points, adjust=True, wrap=False, seed=None, name=NN):
        '''New L{HausdorffEquirectangular} calculator.

           @arg points: Initial set of points, aka the C{model} or
                        C{template} (C{LatLon}[], C{Numpy2LatLon}[],
                        C{Tuple2LatLon}[] or C{other}[]).
           @kwarg adjust: Adjust the wrapped, unrolled longitudinal
                          delta by the cosine of the mean latitude (C{bool}).
           @kwarg wrap: Optionally, wrap and L{unrollPI} longitudes (C{bool}).
           @kwarg seed: Random seed (C{any}) or C{None}, C{0} or
                        C{False} for no U{random sampling<https://
                        Publik.TUWien.ac.AT/files/PubDat_247739.pdf>}.
           @kwarg name: Optional name for this interpolator (C{str}).

           @raise HausdorffError: Insufficient number of B{C{points}} or
                                  invalid B{C{seed}}.
        '''
        HausdorffRadians.__init__(self, points, seed=seed,   name=name,
                                              adjust=adjust, wrap=wrap)

    if _FOR_DOCS:
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric

    def distance(self, p1, p2):
        '''Return the L{equirectangular_} distance in C{radians squared}.
        '''
        r, _ = unrollPI(p1.lam, p2.lam, wrap=self._wrap)
        if self._adjust:
            r *= _scale_rad(p1.phi, p2.phi)
        return hypot2(r, p2.phi - p1.phi)  # like equirectangular_ d2


class HausdorffEuclidean(HausdorffRadians):
    '''Compute the C{Hausdorff} distance based on the C{Euclidean}
       distance in C{radians} from function L{euclidean_}.

       @see: L{HausdorffCosineLaw}, L{HausdorffEquirectangular},
             L{HausdorffFlatPolar}, L{HausdorffHaversine} and
             L{HausdorffVincentys}.
    '''
    _adjust = True
    _wrap   = True

    def __init__(self, points, adjust=True, seed=None, name=NN):
        '''New L{HausdorffEuclidean} calculator.

           @arg points: Initial set of points, aka the C{model} or
                        C{template} (C{LatLon}[], C{Numpy2LatLon}[],
                        C{Tuple2LatLon}[] or C{other}[]).
           @kwarg adjust: Adjust the wrapped, unrolled longitudinal
                          delta by the cosine of the mean latitude (C{bool}).
           @kwarg seed: Random sampling seed (C{any}) or C{None}, C{0}
                        or C{False} for no U{random sampling<https://
                        Publik.TUWien.ac.AT/files/PubDat_247739.pdf>}.
           @kwarg name: Optional name for this interpolator (C{str}).

           @raise HausdorffError: Insufficient number of B{C{points}} or
                                  invalid B{C{seed}}.
        '''
        HausdorffRadians.__init__(self, points, seed=seed, name=name,
                                                           wrap=True)
        if not adjust:
            self._adjust = False

    if _FOR_DOCS:
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric

    def distance(self, p1, p2):
        '''Return the L{euclidean_} distance in C{radians}.
        '''
        return euclidean_(p2.phi, p1.phi, p2.lam - p1.lam, adjust=self._adjust)


class HausdorffFlatLocal(HausdorffRadians):
    '''Compute the C{Hausdorff} distance based on the I{angular} distance
       in C{radians squared} like function L{flatLocal_}/L{hubeny_}.

       @see: L{HausdorffCosineAndoyerLambert},
             L{HausdorffCosineForsytheAndoyerLambert},
             L{HausdorffDistanceTo}, L{HausdorffHubeny},
             L{HausdorffThomas} and L{HausdorffKarney}.
    '''
    _datum    =  Datums.WGS84
    _hubeny2_ =  None
    _units    = _radians2_
    _wrap     =  False

    def __init__(self, points, datum=None, wrap=False, seed=None, name=NN):
        '''New L{HausdorffFlatLocal}/L{HausdorffHubeny} calculator.

           @arg points: Initial set of points, aka the C{model} or
                        C{template} (C{LatLon}[], C{Numpy2LatLon}[],
                        C{Tuple2LatLon}[] or C{other}[]).
           @kwarg datum: Optional datum overriding the default C{Datums.WGS84}
                         and first B{C{points}}' datum (L{Datum}).
           @kwarg wrap: Optionally, wrap and L{unrollPI} longitudes (C{bool}).
           @kwarg seed: Random sampling seed (C{any}) or C{None}, C{0}
                        or C{False} for no U{random sampling<https://
                        Publik.TUWien.ac.AT/files/PubDat_247739.pdf>}.
           @kwarg name: Optional name for this interpolator (C{str}).

           @raise HausdorffError: Insufficient number of B{C{points}} or
                                  invalid B{C{seed}}.

           @raise TypeError: Invalid B{C{datum}}.
        '''
        HausdorffRadians.__init__(self, points, seed=seed, name=name,
                                                           wrap=wrap)
        self._datum_setter(datum)
        self._hubeny2_ = self.datum.ellipsoid._hubeny2_

    if _FOR_DOCS:
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric

    def distance(self, p1, p2):
        '''Return the L{flatLocal_}/L{hubeny_} distance in C{radians squared}.
        '''
        d, _ = unrollPI(p1.lam, p2.lam, wrap=self._wrap)
        return self._hubeny2_(p2.phi, p1.phi, d)


class HausdorffFlatPolar(HausdorffRadians):
    '''Compute the C{Hausdorff} distance based on the I{angular}
       distance in C{radians} from function L{flatPolar_}.

       @see: L{HausdorffCosineLaw}, L{HausdorffEquirectangular},
             L{HausdorffEuclidean}, L{HausdorffHaversine} and
             L{HausdorffVincentys}.
    '''
    _wrap = False

    def __init__(self, points, wrap=False, seed=None, name=NN):
        '''New L{HausdorffFlatPolar} calculator.

           @arg points: Initial set of points, aka the C{model} or
                        C{template} (C{LatLon}[], C{Numpy2LatLon}[],
                        C{Tuple2LatLon}[] or C{other}[]).
           @kwarg wrap: Optionally, wrap and L{unrollPI} longitudes (C{bool}).
           @kwarg seed: Random seed (C{any}) or C{None}, C{0} or
                        C{False} for no U{random sampling<https://
                        Publik.TUWien.ac.AT/files/PubDat_247739.pdf>}.
           @kwarg name: Optional name for this interpolator (C{str}).

           @raise HausdorffError: Insufficient number of B{C{points}} or
                                  invalid B{C{seed}}.
        '''
        HausdorffRadians.__init__(self, points, seed=seed, name=name,
                                                           wrap=wrap)

    if _FOR_DOCS:
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric

    def distance(self, p1, p2):
        '''Return the L{flatPolar_} distance in C{radians}.
        '''
        r, _ = unrollPI(p1.lam, p2.lam, wrap=self._wrap)
        return flatPolar_(p2.phi, p1.phi, r)


class HausdorffHaversine(HausdorffRadians):
    '''Compute the C{Hausdorff} distance based on the I{angular}
       distance in C{radians} from function L{haversine_}.

       @note: See note under L{HausdorffVincentys}.

       @see: L{HausdorffEquirectangular}, L{HausdorffEuclidean},
             L{HausdorffFlatPolar} and L{HausdorffVincentys}.
    '''
    _wrap = False

    def __init__(self, points, wrap=False, seed=None, name=NN):
        '''New L{HausdorffHaversine} calculator.

           @arg points: Initial set of points, aka the C{model} or
                        C{template} (C{LatLon}[], C{Numpy2LatLon}[],
                        C{Tuple2LatLon}[] or C{other}[]).
           @kwarg wrap: Optionally, wrap and L{unrollPI} longitudes (C{bool}).
           @kwarg seed: Random sampling seed (C{any}) or C{None}, C{0}
                        or C{False} for no U{random sampling<https://
                        Publik.TUWien.ac.AT/files/PubDat_247739.pdf>}.
           @kwarg name: Optional name for this interpolator (C{str}).

           @raise HausdorffError: Insufficient number of B{C{points}} or
                                  invalid B{C{seed}}.
        '''
        HausdorffRadians.__init__(self, points, seed=seed, name=name,
                                                           wrap=wrap)

    if _FOR_DOCS:
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric

    def distance(self, p1, p2):
        '''Return the L{haversine_} distance in C{radians}.
        '''
        r, _ = unrollPI(p1.lam, p2.lam, wrap=self._wrap)
        return haversine_(p2.phi, p1.phi, r)


class HausdorffHubeny(HausdorffFlatLocal):  # for Karl Huebny
    if _FOR_DOCS:
        __doc__   = HausdorffFlatLocal.__doc__
        __init__  = HausdorffFlatLocal.__init__
        directed  = HausdorffFlatLocal.directed
        distance  = HausdorffFlatLocal.distance
        symmetric = HausdorffFlatLocal.symmetric


class HausdorffKarney(HausdorffDegrees):
    '''Compute the C{Hausdorff} distance based on the I{angular}
       distance in C{degrees} from I{Karney}'s U{GeographicLib
       <https://PyPI.org/project/geographiclib>} U{Geodesic
       <https://GeographicLib.SourceForge.io/html/python/code.html>}
       Inverse method.

       @see: L{HausdorffCosineAndoyerLambert},
             L{HausdorffCosineForsytheAndoyerLambert},
             L{HausdorffDistanceTo}, L{HausdorffFlatLocal},
             L{HausdorffHubeny} and L{HausdorffThomas}.
    '''
    _datum    =  Datums.WGS84
    _Inverse1 =  None
    _units    = _degrees_
    _wrap     =  False

    def __init__(self, points, datum=None, wrap=False, seed=None, name=NN):
        '''New L{HausdorffKarney} calculator.

           @arg points: Initial set of points, aka the C{model} or
                        C{template} (C{LatLon}[], C{Numpy2LatLon}[],
                        C{Tuple2LatLon}[] or C{other}[]).
           @kwarg datum: Optional datum overriding the default C{Datums.WGS84}
                         and first B{C{points}}' datum (L{Datum}).
           @kwarg wrap: Optionally, wrap and L{unroll180} longitudes (C{bool}).
           @kwarg seed: Random sampling seed (C{any}) or C{None}, C{0}
                        or C{False} for no U{random sampling<https://
                        Publik.TUWien.ac.AT/files/PubDat_247739.pdf>}.
           @kwarg name: Optional name for this interpolator (C{str}).

           @raise HausdorffError: Insufficient number of B{C{points}} or
                                  invalid B{C{seed}}.

           @raise ImportError: Package U{GeographicLib
                  <https://PyPI.org/project/geographiclib>} missing.

           @raise TypeError: Invalid B{C{datum}}.
        '''
        HausdorffDegrees.__init__(self, points, seed=seed, name=name,
                                                           wrap=wrap)
        self._datum_setter(datum)
        self._Inverse1 = self.datum.ellipsoid.geodesic.Inverse1

    if _FOR_DOCS:
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric

    def distance(self, p1, p2):
        '''Return the non-negative I{angular} distance in C{degrees}.
        '''
        return self._Inverse1(p1.lat, p1.lon, p2.lat, p2.lon, wrap=self._wrap)


class HausdorffThomas(HausdorffRadians):
    '''Compute the C{Hausdorff} distance based on the I{angular}
       distance in C{radians} from function L{thomas_}.

       @see: L{HausdorffCosineAndoyerLambert},
             L{HausdorffCosineForsytheAndoyerLambert},
             L{HausdorffDistanceTo}, L{HausdorffFlatLocal},
             L{HausdorffHubeny} and L{HausdorffKarney}.
    '''
    _datum = Datums.WGS84
    _wrap  = False

    def __init__(self, points, datum=None, wrap=False, seed=None, name=NN):
        '''New L{HausdorffThomas} calculator.

           @arg points: Initial set of points, aka the C{model} or
                        C{template} (C{LatLon}[], C{Numpy2LatLon}[],
                        C{Tuple2LatLon}[] or C{other}[]).
           @kwarg datum: Optional datum overriding the default C{Datums.WGS84}
                         and first B{C{points}}' datum (L{Datum}).
           @kwarg wrap: Optionally, wrap and L{unrollPI} longitudes (C{bool}).
           @kwarg seed: Random seed (C{any}) or C{None}, C{0} or
                        C{False} for no U{random sampling<https://
                        Publik.TUWien.ac.AT/files/PubDat_247739.pdf>}.
           @kwarg name: Optional name for this interpolator (C{str}).

           @raise HausdorffError: Insufficient number of B{C{points}} or
                                  invalid B{C{seed}}.
        '''
        HausdorffRadians.__init__(self, points, seed=seed, name=name,
                                                           wrap=wrap)
        self._datum_setter(datum)

    if _FOR_DOCS:
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric

    def distance(self, p1, p2):
        '''Return the L{thomas_} distance in C{radians}.
        '''
        r, _ = unrollPI(p1.lam, p2.lam, wrap=self._wrap)
        return thomas_(p2.phi, p1.phi, r, datum=self._datum)


class HausdorffVincentys(HausdorffRadians):
    '''Compute the C{Hausdorff} distance based on the I{angular}
       distance in C{radians} from function L{vincentys_}.

       @note: See note at function L{vincentys_}.

       @see: L{HausdorffCosineLaw}, L{HausdorffEquirectangular},
             L{HausdorffEuclidean}, L{HausdorffFlatPolar} and
             L{HausdorffHaversine}.
    '''
    _wrap = False

    def __init__(self, points, wrap=False, seed=None, name=NN):
        '''New L{HausdorffVincentys} calculator.

           @arg points: Initial set of points, aka the C{model} or
                        C{template} (C{LatLon}[], C{Numpy2LatLon}[],
                        C{Tuple2LatLon}[] or C{other}[]).
           @kwarg wrap: Optionally, wrap and L{unrollPI} longitudes (C{bool}).
           @kwarg seed: Random sampling seed (C{any}) or C{None}, C{0}
                        or C{False} for no U{random sampling<https://
                        Publik.TUWien.ac.AT/files/PubDat_247739.pdf>}.
           @kwarg name: Optional name for this interpolator (C{str}).

           @raise HausdorffError: Insufficient number of B{C{points}} or
                                  invalid B{C{seed}}.
        '''
        HausdorffRadians.__init__(self, points, seed=seed, name=name,
                                                           wrap=wrap)

    if _FOR_DOCS:
        directed  = Hausdorff.directed
        symmetric = Hausdorff.symmetric

    def distance(self, p1, p2):
        '''Return the L{vincentys_} distance in C{radians}.
        '''
        r, _ = unrollPI(p1.lam, p2.lam, wrap=self._wrap)
        return vincentys_(p2.phi, p1.phi, r)


def _hausdorff_(ps1, ps2, both, early, seed, units, distance, point):
    '''(INTERNAL) Core of function L{hausdorff.hausdorff} and methods
       C{directed} and C{symmetric} of classes C{hausdorff.Hausdorff...}.
    '''
    # shuffling the points generally increases the
    # chance of an early break in the inner j loop
    rr = randomrangenerator(seed) if seed else range

    hd = -INF
    hi = hj = m = mn = 0
    md = 0.0

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
    '''Compute the C{directed} or C{symmetric} U{Hausdorff distance<https://
       WikiPedia.org/wiki/Hausdorff_distance>} between 2 sets of points with or
       without U{early breaking<https://Publik.TUWien.ac.AT/files/PubDat_247739.pdf>}
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
       @kwarg units: Optional, name of the distance units (C{str}).
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

    _, ps1 = points2(model,  closed=False, Error=HausdorffError)
    _, ps2 = points2(target, closed=False, Error=HausdorffError)
    return _hausdorff_(ps1, ps2, both, early, seed, units, distance, point)


def randomrangenerator(seed):
    '''Return a C{seed}ed random range function generator.

       @arg seed: Initial, internal L{Random} state (C{hashable}
                  or C{None}).

       @note: L{Random} B{C{seed=None}} seeds from the current
              time or from a platform-specific randomness source,
              if available.

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


__all__ += _ALL_DOCS(Hausdorff6Tuple)

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
