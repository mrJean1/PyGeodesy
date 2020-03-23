
# -*- coding: utf-8 -*-

u'''Classes L{Frechet}, L{FrechetDegrees}, L{FrechetRadians},
L{FrechetEquirectangular}, L{FrechetEuclidean}, L{FrechetHaversine},
and L{FrechetVincentys} to compute I{discrete} U{Fréchet
<https://WikiPedia.org/wiki/Frechet_distance>} distances between two
sets of C{LatLon}, C{NumPy}, C{tuples} or other types of points.

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

from pygeodesy.basics import EPS, EPS1, INF, _isnotError, \
                             isscalar, property_doc_, property_RO
from pygeodesy.fmath import favg, hypot2
from pygeodesy.formy import euclidean_, haversine_, points2 as _points2, \
                           _scaler, vincentys_
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY
from pygeodesy.named import LatLon2Tuple, _Named, _NamedTuple, \
                            notOverloaded, PhiLam2Tuple
from pygeodesy.utily import unrollPI

from collections import defaultdict
from math import radians

__all__ = _ALL_LAZY.frechet + _ALL_DOCS('Frechet6Tuple')
__version__ = '20.03.23'

if not 0 < EPS < EPS1 < 1:
    raise AssertionError('%s < %s: 0 < %.6g < %.6g < 1' % ('EPS', 'EPS1', EPS, EPS1))


class FrechetError(ValueError):
    '''Fréchet issue.
    '''
    pass


def _fraction(fraction, n):
    f = 1  # int, no fractional indices
    if fraction in (None, 1):
        pass
    elif not (isscalar(fraction) and EPS < fraction < EPS1
                                 and (float(n) - fraction) < n):
        raise FrechetError('%s invalid: %r' % ('fraction', fraction))
    elif fraction < EPS1:
        f = float(fraction)
    return f


class Frechet6Tuple(_NamedTuple):
    '''6-Tuple C{(fd, fi1, fi2, r, n, units)} with the I{discrete}
       U{Fréchet<https://WikiPedia.org/wiki/Frechet_distance>} distance
       C{fd}, I{fractional} indices C{fi1} and C{fi2}, the recursion
       depth C{r}, the number of distances computed C{n} and the name
       of the distance C{units}.

       If I{fractional} indices C{fi1} and C{fi2} are type C{int}, the
       returned C{fd} is the distance between C{points1[fi1]} and
       C{points2[fi2]}.  For type C{float} indices, the distance is
       between an intermediate point along C{points1[int(fi1)]} and
       C{points1[int(fi1) + 1]} respectively an intermediate point
       along C{points2[int(fi2)]} and C{points2[int(fi2) + 1]}.

       Use function L{fractional} to compute the point at a fractional
       index.
    '''
    _Names_ = ('fd', 'fi1', 'fi2', 'r', 'n', 'units')

#   def __gt__(self, other):
#       _TypeError(Frechet6Tuple, other=other)
#       return self if self.fd > other.fd else other  # PYCHOK .fd=[0]
#
#   def __lt__(self, other):
#       _TypeError(Frechet6Tuple, other=other)
#       return self if self.fd < other.fd else other  # PYCHOK .fd=[0]


class Frechet(_Named):
    '''Frechet base class, requires method L{Frechet.distance} to
       be overloaded.
    '''
    _adjust = None
    _f1     = 1
    _n1     = 0
    _ps1    = None
    _units  = ''
    _wrap   = None

    def __init__(self, points, fraction=None, name='', units=''):
        '''New L{Frechet} calculator/interpolator.

           @arg points: First set of points (C{LatLon}[], L{Numpy2LatLon}[],
                        L{Tuple2LatLon}[] or C{other}[]).
           @kwarg fraction: Index fraction (C{float} in L{EPS}..L{EPS1}) to
                            interpolate intermediate B{C{points}} or use
                            C{None}, C{0} or C{1} for no intermediate
                            B{C{points}} and no I{fractional} indices.
           @kwarg name: Optional calculator/interpolator name (C{str}).
           @kwarg units: Optional, distance units (C{str}).

           @raise FrechetError: Insufficient number of B{C{points}} or
                                invalid B{C{fraction}}.
        '''
        self._n1, self._ps1 = _points2(points, closed=False, Error=FrechetError)
        if fraction:
            self.fraction = fraction
        if name:
            self.name = name
        if units:
            self.units = units

    @property_RO
    def adjust(self):
        '''Get the adjust setting (C{bool} or C{None} if not applicable).
        '''
        return self._adjust

    def discrete(self, points, fraction=None):
        '''Compute the C{forward, discrete Fréchet} distance.

           @arg points: Second set of points (C{LatLon}[], L{Numpy2LatLon}[],
                        L{Tuple2LatLon}[] or C{other}[]).
           @kwarg fraction: Index fraction (C{float} in L{EPS}..L{EPS1}) to
                            interpolate intermediate B{C{points}} or use
                            C{None}, C{0} or C{1} for no intermediate
                            B{C{points}} and no I{fractional} indices.

           @return: A L{Frechet6Tuple}C{(fd, fi1, fi2, r, n, units)}.

           @raise FrechetError: Insufficient number of B{C{points}} or
                                invalid B{C{fraction}}.

           @raise RecursionError: Recursion depth exceeded, see U{sys.getrecursionlimit()
                                  <https://docs.Python.org/3/library/sys.html#sys.getrecursionlimit>}.
        '''
        n2, ps2 = _points2(points, closed=False, Error=FrechetError)

        f2 = _fraction(fraction, n2)
        p2 = self.points_fraction if f2 < EPS1 else self.points_  # PYCHOK expected

        f1 = self.fraction
        p1 = self.points_fraction if f1 < EPS1 else self.points_  # PYCHOK expected

        def dF(fi1, fi2):
            return self.distance(p1(self._ps1, fi1), p2(ps2, fi2))

        return _frechet_(self._n1, f1, n2, f2, dF, self.units)

    def distance(self, point1, point2):
        '''(INTERNAL) I{must be overloaded}.
        '''
        notOverloaded(self, self.distance.__name__, point1, point2)

    @property_doc_(''' the index fraction (C{float}).''')
    def fraction(self):
        '''Get the index fraction (C{float} or C{1}).
        '''
        return self._f1

    @fraction.setter  # PYCHOK setter!
    def fraction(self, fraction):
        '''Set the the index fraction (C{float} or C{1}).

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

    def points_fraction(self, points, fi):
        '''Get and convert I{fractional} index point for the C{.distance} method.

           @arg points: The orignal B{C{points}} to convert ((C{LatLon}[],
                        L{Numpy2LatLon}[], L{Tuple2LatLon}[] or C{other}[]).
           @arg fi: The I{fractional} index in B{C{points}} (C{float} or C{int}).

           @return: The interpolated, converted, intermediate B{C{points[fi]}}.
        '''
        return self.point(_fractional(points, fi))

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
        self._units = str(units or "")

    @property_RO
    def wrap(self):
        '''Get the wrap setting (C{bool} or C{None} if not applicable).
        '''
        return self._wrap


class FrechetDegrees(Frechet):
    '''L{Frechet} base class for distances in C{degrees} from
       C{LatLon} points in C{degrees}.
    '''
    _units = 'degrees'

    discrete = Frechet.discrete  # for __doc__


class FrechetRadians(Frechet):
    '''L{Frechet} base class for distances in C{radians} from
       C{LatLon} points converted from C{degrees} to C{radians}.
    '''
    _units = 'radians'

    discrete = Frechet.discrete  # for __doc__

    def point(self, point):
        '''Convert C{(lat, lon)} point in degrees to C{(a, b)}
           in radians.

           @return: An L{PhiLam2Tuple}C{(phi, lam)}.
        '''
        return PhiLam2Tuple(radians(point.lat), radians(point.lon))


class FrechetEquirectangular(FrechetRadians):
    '''Compute the C{Frechet} distance based on the C{equirectangular}
       distance (in radians squared) like function L{equirectangular_}.

       @see: L{FrechetEuclidean}, L{FrechetHaversine} and L{FrechetVincentys}.
    '''
    _adjust = True
    _wrap   = False

    def __init__(self, points, adjust=True, wrap=False, fraction=None, name=''):
        '''New L{FrechetEquirectangular} calculator/interpolator.

           @arg points: First set of points (C{LatLon}[], L{Numpy2LatLon}[],
                        L{Tuple2LatLon}[] or C{other}[]).
           @kwarg adjust: Adjust the wrapped, unrolled longitudinal
                          delta by the cosine of the mean latitude (C{bool}).
           @kwarg wrap: Wrap and L{unroll180} longitudes (C{bool}).
           @kwarg fraction: Index fraction (C{float} in L{EPS}..L{EPS1}) to
                            interpolate intermediate B{C{points}} or use
                            C{None}, C{0} or C{1} for no intermediate
                            B{C{points}} and no I{fractional} indices.
           @kwarg name: Optional calculator/interpolator name (C{str}).

           @raise FrechetError: Insufficient number of B{C{points}} or
                                invalid B{C{adjust}} or B{C{seed}}.
        '''
        if not adjust:
            self._adjust = False
        if wrap:
            self._wrap = True
        super(FrechetRadians, self).__init__(points, fraction=fraction, name=name)  # distance**2

    def distance(self, p1, p2):
        '''Return the L{equirectangular_} distance in C{radians squared}.
        '''
        d, _ = unrollPI(p1.lam, p2.lam, wrap=self.wrap)
        if self.adjust:
            d *= _scaler(p1.phi, p2.phi)
        return hypot2(d, p2.phi - p1.phi)  # like equirectangular_ d2

    discrete = Frechet.discrete  # for __doc__


class FrechetEuclidean(FrechetRadians):
    '''Compute the C{Frechet} distance based on the C{Euclidean}
       distance (in radians) from function L{euclidean_}.

       @see: L{FrechetEquirectangular}, L{FrechetHaversine} and L{FrechetVincentys}.
    '''
    _adjust = True

    def __init__(self, points, adjust=True, fraction=None, name=''):
        '''New L{FrechetEuclidean} calculator/interpolator.

           @arg points: First set of points (C{LatLon}[], L{Numpy2LatLon}[],
                        L{Tuple2LatLon}[] or C{other}[]).
           @kwarg adjust: Adjust the wrapped, unrolled longitudinal
                          delta by the cosine of the mean latitude (C{bool}).
           @kwarg fraction: Index fraction (C{float} in L{EPS}..L{EPS1}) to
                            interpolate intermediate B{C{points}} or use
                            C{None}, C{0} or C{1} for no intermediate
                            B{C{points}} and no I{fractional} indices.
           @kwarg name: Optional calculator/interpolator name (C{str}).

           @raise FrechetError: Insufficient number of B{C{points}} or
                                invalid B{C{fraction}}.
        '''
        if not adjust:
            self._adjust = False
        super(FrechetRadians, self).__init__(points, fraction=fraction, name=name)

    def distance(self, p1, p2):
        '''Return the L{euclidean_} distance in C{radians}.
        '''
        return euclidean_(p2.phi, p1.phi, p2.lam - p1.lam, adjust=self.adjust)

    discrete = Frechet.discrete  # for __doc__


class FrechetHaversine(FrechetRadians):
    '''Compute the C{Frechet} distance based on the I{angular}
       C{Haversine} distance (in radians) from function L{haversine_}.

       @note: See note under L{FrechetVincentys}.

       @see: L{FrechetEquirectangular}, L{FrechetEuclidean} and L{FrechetVincentys}.
    '''
    _wrap = False

    def __init__(self, points, wrap=False, fraction=None, name=''):
        '''New L{FrechetHaversine} calculator/interpolator.

           @arg points: First set of points (C{LatLon}[], L{Numpy2LatLon}[],
                        L{Tuple2LatLon}[] or C{other}[]).
           @kwarg wrap: Wrap and L{unroll180} longitudes (C{bool}).
           @kwarg fraction: Index fraction (C{float} in L{EPS}..L{EPS1}) to
                            interpolate intermediate B{C{points}} or use
                            C{None}, C{0} or C{1} for no intermediate
                            B{C{points}} and no I{fractional} indices.
           @kwarg name: Optional calculator/interpolator name (C{str}).

           @raise FrechetError: Insufficient number of B{C{points}} or
                                invalid B{C{fraction}}.
        '''
        if wrap:
            self._wrap = True
        super(FrechetRadians, self).__init__(points, fraction=fraction, name=name)

    def distance(self, p1, p2):
        '''Return the L{haversine_} distance in C{radians}.
        '''
        d, _ = unrollPI(p1.lam, p2.lam, wrap=self.wrap)
        return haversine_(p2.phi, p1.phi, d)

    discrete = Frechet.discrete  # for __doc__


class FrechetVincentys(FrechetRadians):
    '''Compute the C{Frechet} distance based on the I{angular}
       C{Vincenty} distance (in radians) from function L{vincentys_}.

       @note: See note under L{vincentys_}.

       @see: L{FrechetEquirectangular}, L{FrechetEuclidean} and L{FrechetHaversine}.
    '''
    _wrap = False

    def __init__(self, points, wrap=False, fraction=None, name=''):
        '''New L{FrechetVincentys} calculator/interpolator.

           @arg points: First set of points (C{LatLon}[], L{Numpy2LatLon}[],
                        L{Tuple2LatLon}[] or C{other}[]).
           @kwarg wrap: Wrap and L{unroll180} longitudes (C{bool}).
           @kwarg fraction: Index fraction (C{float} in L{EPS}..L{EPS1}) to
                            interpolate intermediate B{C{points}} or use
                            C{None}, C{0} or C{1} for no intermediate
                            B{C{points}} and no I{fractional} indices.
           @kwarg name: Optional calculator/interpolator name (C{str}).

           @raise FrechetError: Insufficient number of B{C{points}} or
                                invalid B{C{fraction}}.
        '''
        if wrap:
            self._wrap = True
        super(FrechetRadians, self).__init__(points, fraction=fraction, name=name)

    def distance(self, p1, p2):
        '''Return the L{vincentys_} distance in C{radians}.
        '''
        d, _ = unrollPI(p1.lam, p2.lam, wrap=self.wrap)
        return vincentys_(p2.phi, p1.phi, d)

    discrete = Frechet.discrete  # for __doc__


def _fractional(points, fi):
    '''(INTERNAL) Compute point at I{fractional} index.
    '''
    i = int(fi)
    p = points[i]
    f = fi - float(i)
    if f > EPS:
        if f < EPS1:
            q = points[i + 1]
            p = LatLon2Tuple(favg(p.lat, q.lat, f=f),
                             favg(p.lon, q.lon, f=f))
        else:
            p = points[i + 1]
    return p


def fractional(points, fi, LatLon=None, **LatLon_kwds):
    '''Return the point at a given I{fractional} index.

       @arg points: The points (C{LatLon}[], L{Numpy2LatLon}[],
                    L{Tuple2LatLon}[] or C{other}[]).
       @arg fi: The fractional index (C{float} or C{int}).
       @kwarg LatLon: Optional class to return the I{intermediate}
                      point (C{LatLon}) or C{None}.
       @kwarg LatLon_kwds: Optional B{C{LatLon}} keyword arguments,
                           ignored of B{C{LatLon=None}}.

       @return: A B{C{LatLon}} or if B{C{LatLon}} is C{None}, a
                L{LatLon2Tuple}C{(lat, lon)} for B{C{points[fi]}} if
                I{fractional} index B{C{fi}} is C{int}, otherwise the
                intermediate point between B{C{points[int(fi)]}} and
                B{C{points[int(fi) + 1]}} for C{float} I{fractional}
                index B{C{fi}}.

       @raise IndexError: Fractional index B{C{fi}} invalid or
                          B{C{points}} not subscriptable.
    '''
    try:
        if not (isscalar(fi) and 0 <= fi < len(points)):
            raise IndexError
        p = _fractional(points, fi)
    except (IndexError, TypeError):
        raise IndexError('%s invalid: %r' % (fractional.__name__, fi))

    if LatLon and isinstance(p, LatLon2Tuple):
        p = LatLon(*p, **LatLon_kwds)
    return p


def _frechet_(ni, fi, nj, fj, dF, units):  # MCCABE 14
    '''(INTERNAL) Recursive core of function L{frechet_}
       and method C{discrete} of C{Frechet...} classes.
    '''
    iFs = {}

    def iF(i):  # cache index, depth ints and floats
        return iFs.setdefault(i, i)

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

    cF = defaultdict(dict)
    t  = rF(ni - 1, nj - 1, 0)
    t += (sum(map(len, cF.values())), units)
#   del cF, iFs
    return Frechet6Tuple(*t)


def frechet_(points1, points2, distance=None, units=''):
    '''Compute the I{discrete} U{Fréchet<https://WikiPedia.org/wiki/Frechet_distance>}
       distance between two paths given as sets of points.

       @arg points1: First set of points (C{LatLon}[], L{Numpy2LatLon}[],
                     L{Tuple2LatLon}[] or C{other}[]).
       @arg points2: Second set of points (C{LatLon}[], L{Numpy2LatLon}[],
                     L{Tuple2LatLon}[] or C{other}[]).
       @kwarg distance: Callable returning the distance between a B{C{points1}}
                        and a B{C{points2}} point (signature C{(point1, point2)}).
       @kwarg units: Optional, name of the distance units (C{str}).

       @return: A L{Frechet6Tuple}C{(fd, fi1, fi2, r, n, units)} where C{fi1}
                and C{fi2} are type C{int} indices into B{C{points1}} respectively
                B{C{points2}}.

       @raise FrechetError: Insufficient number of B{C{points1}} or B{C{points2}}.

       @raise RecursionError: Recursion depth exceeded, see U{sys.getrecursionlimit()
                              <https://docs.Python.org/3/library/sys.html#sys.getrecursionlimit>}.

       @raise TypeError: If B{C{distance}} is not a callable.

       @note: Keyword arguments C{fraction}, intermediate B{C{points1}} and
              B{C{points2}} and I{fractional} indices are I{not} supported in
              this L{frechet_} function.
    '''
    if not callable(distance):
        raise _isnotError(callable.__name__, distance=distance)

    n1, ps1 = _points2(points1, closed=False, Error=FrechetError)
    n2, ps2 = _points2(points2, closed=False, Error=FrechetError)

    def dF(i1, i2):
        return distance(ps1[i1], ps2[i2])

    return _frechet_(n1, 1, n2, 1, dF, units)

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
