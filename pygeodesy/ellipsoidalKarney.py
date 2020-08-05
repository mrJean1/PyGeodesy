
# -*- coding: utf-8 -*-

u'''Ellipsoidal classes geodetic (lat-/longitude) L{LatLon} and
geocentric (ECEF) L{Cartesian} and functions L{areaOf}, L{isclockwise}
and L{perimeterOf}, all based on I{Charles Karney}'s Python implementation
of U{GeographicLib<https://PyPI.org/project/geographiclib>}.

Here's an example usage of C{ellipsoidalKarney}:

    >>> from pygeodesy.ellipsoidalKarney import LatLon
    >>> Newport_RI = LatLon(41.49008, -71.312796)
    >>> Cleveland_OH = LatLon(41.499498, -81.695391)
    >>> Newport_RI.distanceTo(Cleveland_OH)
    866,455.4329098687  # meter

You can change the ellipsoid model used by the Karney formulae
as follows:

    >>> from pygeodesy import Datums
    >>> from pygeodesy.ellipsoidalKarney import LatLon
    >>> p = LatLon(0, 0, datum=Datums.OSGB36)

or by converting to anothor datum:

    >>> p = p.convertDatum(Datums.OSGB36)

@newfield example: Example, Examples
'''

from pygeodesy.basics import property_RO, _xkwds
from pygeodesy.datum import Datums
from pygeodesy.ecef import EcefKarney
from pygeodesy.ellipsoidalBase import _intersections2, _TOL_M, \
                                       CartesianEllipsoidalBase, \
                                       LatLonEllipsoidalBase
from pygeodesy.errors import _ValueError
from pygeodesy.formy import points2
from pygeodesy.lazily import _ALL_LAZY, _ALL_OTHER
from pygeodesy.named import Bearing2Tuple, Destination2Tuple
from pygeodesy.points import _areaError, ispolar  # PYCHOK exported
from pygeodesy.utily import unroll180, wrap90, wrap180, wrap360

__all__ = _ALL_LAZY.ellipsoidalKarney
__version__ = '20.08.04'


class Cartesian(CartesianEllipsoidalBase):
    '''Extended to convert C{Karney}-based L{Cartesian} to
       C{Karney}-based L{LatLon} points.
    '''
    _Ecef = EcefKarney  #: (INTERNAL) Preferred C{Ecef...} class.

    def toLatLon(self, **LatLon_datum_kwds):  # PYCHOK LatLon=LatLon, datum=None
        '''Convert this cartesian point to a C{Karney}-based
           geodetic point.

           @kwarg LatLon_datum_kwds: Optional L{LatLon}, B{C{datum}} and
                  other keyword arguments, ignored if B{C{LatLon=None}}.
                  Use B{C{LatLon=...}} to override this L{LatLon} class
                  or specify B{C{LatLon=None}}.

           @return: The geodetic point (L{LatLon}) or if B{C{LatLon}}
                    is C{None}, an L{Ecef9Tuple}C{(x, y, z, lat, lon,
                    height, C, M, datum)} with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{LatLon}}, B{C{datum}} or other
                             B{C{LatLon_datum_kwds}}.
        '''
        kwds = _xkwds(LatLon_datum_kwds, LatLon=LatLon, datum=self.datum)
        return CartesianEllipsoidalBase.toLatLon(self, **kwds)


class LatLon(LatLonEllipsoidalBase):
    '''An ellipsoidal L{LatLon} similar to L{ellipsoidalVincenty.LatLon}
       but using I{Charles F. F. Karney}'s Python U{GeographicLib
       <https://PyPI.org/project/geographiclib>} to compute the geodesic
       distance, initial and final bearing (azimuths) between two given
       points or the destination point given a start point and an initial
       bearing.

       @note: This L{LatLon}'s methods require the U{GeographicLib
              <https://PyPI.org/project/geographiclib>} package.
    '''
    _Ecef = EcefKarney  #: (INTERNAL) Preferred C{Ecef...} class.

    def bearingTo(self, other, wrap=False):  # PYCHOK no cover
        '''DEPRECATED, use method C{initialBearingTo}.
        '''
        return self.initialBearingTo(other, wrap=wrap)

    def bearingTo2(self, other, wrap=False):
        '''Compute the initial and final bearing (forward and reverse
           azimuth) from this to an other point, using I{Karney}'s
           C{Inverse} method.  See methods L{initialBearingTo} and
           L{finalBearingTo} for more details.

           @arg other: The other point (L{LatLon}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: A L{Bearing2Tuple}C{(initial, final)}.

           @raise ImportError: Package U{geographiclib
                               <https://PyPI.org/project/geographiclib>}
                               not installed or not found.

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: If this and the B{C{other}} point's L{Datum}
                              ellipsoids are not compatible.
        '''
        r = self._inverse(other, wrap)
        return self._xnamed(Bearing2Tuple(r.initial, r.final))

    def destination(self, distance, bearing, height=None):
        '''Compute the destination point after having travelled
           for the given distance from this point along a geodesic
           given by an initial bearing, using I{Karney}'s C{Direct}
           method.  See method L{destination2} for more details.

           @arg distance: Distance (C{meter}).
           @arg bearing: Initial bearing in (compass C{degrees360}).
           @kwarg height: Optional height, overriding the default
                          height (C{meter}, same units as C{distance}).

           @return: The destination point (L{LatLon}).

           @raise ImportError: Package U{geographiclib
                               <https://PyPI.org/project/geographiclib>}
                               not installed or not found.

           @example:

           >>> p = LatLon(-37.95103, 144.42487)
           >>> d = p.destination(54972.271, 306.86816)
           >>> d
           LatLon(37°39′10.14″S, 143°55′35.39″E)  # 37.652818°S, 143.926498°E
        '''
        r = self._direct(distance, bearing, self.classof, height)
        return r.destination

    def destination2(self, distance, bearing, height=None):
        '''Compute the destination point and the final bearing (reverse
           azimuth) after having travelled for the given distance from
           this point along a geodesic given by an initial bearing,
           using I{Karney}'s C{Direct} method.

           The distance must be in the same units as this point's datum
           axes, conventionally C{meter}.  The distance is measured on
           the surface of the ellipsoid, ignoring this point's height.

           The initial and final bearing (forward and reverse azimuth)
           are in compass C{degrees360}.

           The destination point's height and datum are set to this
           point's height and datum, unless the former is overridden.

           @arg distance: Distance (C{meter}).
           @arg bearing: Initial bearing (compass C{degrees360}).
           @kwarg height: Optional height, overriding the default
                          height (C{meter}, same units as C{distance}).

           @return: A L{Destination2Tuple}C{(destination, final)}.

           @raise ImportError: Package U{geographiclib
                               <https://PyPI.org/project/geographiclib>}
                               not installed or not found.

           @example:

           >>> p = LatLon(-37.95103, 144.42487)
           >>> d, f = p.destination2(54972.271, 306.86816)
           >>> d
           LatLon(37°39′10.14″S, 143°55′35.39″E)  # 37.652818°S, 143.926498°E
           >>> f
           307.1736313846665
        '''
        r = self._direct(distance, bearing, self.classof, height)
        return self._xnamed(r)

    def distanceTo(self, other, wrap=False, **unused):  # for -DistanceTo
        '''Compute the distance between this and an other point
           along a geodesic, using I{Karney}'s C{Inverse} method.
           See method L{distanceTo3} for more details.

           @arg other: The other point (L{LatLon}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Distance (C{meter}).

           @raise ImportError: Package U{geographiclib
                               <https://PyPI.org/project/geographiclib>}
                               not installed or not found.

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: If this and the B{C{other}} point's L{Datum}
                              ellipsoids are not compatible.

           @example:

           >>> p = LatLon(50.06632, -5.71475)
           >>> q = LatLon(58.64402, -3.07009)
           >>> d = p.distanceTo(q)  # 969,954.1663142084 m
        '''
        return self._inverse(other, wrap).distance

    def distanceTo3(self, other, wrap=False):
        '''Compute the distance, the initial and final bearing along a
           geodesic between this and an other point, using I{Karney}'s
           C{Inverse} method.

           The distance is in the same units as this point's datum axes,
           conventionally meter.  The distance is measured on the surface
           of the ellipsoid, ignoring this point's height.

           The initial and final bearing (forward and reverse azimuth)
           are in compass C{degrees360} from North.

           @arg other: Destination point (L{LatLon}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: A L{Distance3Tuple}C{(distance, initial, final)}.

           @raise ImportError: Package U{geographiclib
                               <https://PyPI.org/project/geographiclib>}
                               not installed or not found.

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: If this and the B{C{other}} point's L{Datum}
                              ellipsoids are not compatible.
        '''
        return self._xnamed(self._inverse(other, wrap))

    def finalBearingOn(self, distance, bearing):
        '''Compute the final bearing (reverse azimuth) after having
           travelled for the given distance along a geodesic given
           by an initial bearing from this point, using I{Karney}'s
           C{Direct} method.  See method L{destination2} for more details.

           @arg distance: Distance (C{meter}).
           @arg bearing: Initial bearing (compass C{degrees360}).

           @return: Final bearing (compass C{degrees360}).

           @raise ImportError: Package U{geographiclib
                               <https://PyPI.org/project/geographiclib>}
                               not installed or not found.

           @example:

           >>> p = LatLon(-37.95103, 144.42487)
           >>> b = 306.86816
           >>> f = p.finalBearingOn(54972.271, b)  # 307.1736313846665°
        '''
        return self._direct(distance, bearing, None, None).final

    def finalBearingTo(self, other, wrap=False):
        '''Compute the final bearing (reverse azimuth) after having
           travelled along a geodesic from this point to an other
           point, using I{Karney}'s C{Inverse} method.  See method
           L{distanceTo3} for more details.

           @arg other: The other point (L{LatLon}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Final bearing (compass C{degrees360}).

           @raise ImportError: Package U{geographiclib
                               <https://PyPI.org/project/geographiclib>}
                               not installed or not found.

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: If this and the B{C{other}} point's L{Datum}
                              ellipsoids are not compatible.

           @example:

           >>> p = new LatLon(50.06632, -5.71475)
           >>> q = new LatLon(58.64402, -3.07009)
           >>> f = p.finalBearingTo(q)  # 11.297220414306684°

           >>> p = LatLon(52.205, 0.119)
           >>> q = LatLon(48.857, 2.351)
           >>> f = p.finalBearingTo(q)  # 157.83449958372714°
        '''
        return self._inverse(other, wrap).final

    @property_RO
    def geodesic(self):
        '''Get this C{LatLon}'s I{wrapped} U{Karney Geodesic
           <https://GeographicLib.SourceForge.io/html/python/code.html>},
           provided package U{geographiclib
           <https://PyPI.org/project/geographiclib>} is installed.
        '''
        return self.datum.ellipsoid.geodesic

    def initialBearingTo(self, other, wrap=False):
        '''Compute the initial bearing (forward azimuth) to travel
           along a geodesic from this point to an other point,
           using I{Karney}'s C{Inverse} method.  See method
           L{distanceTo3} for more details.

           @arg other: The other point (L{LatLon}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Initial bearing (compass C{degrees360}).

           @raise ImportError: Package U{geographiclib
                               <https://PyPI.org/project/geographiclib>}
                               not installed or not found.

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: If this and the B{C{other}} point's L{Datum}
                              ellipsoids are not compatible.

           @example:

           >>> p = LatLon(50.06632, -5.71475)
           >>> q = LatLon(58.64402, -3.07009)
           >>> b = p.initialBearingTo(q)  # 9.141877488906045°

           >>> p = LatLon(52.205, 0.119)
           >>> q = LatLon(48.857, 2.351)
           >>> b = p.initialBearingTo(q)  # 156.1106404059787°

           @JSname: I{bearingTo}.
        '''
        return self._inverse(other, wrap).initial

    def toCartesian(self, **Cartesian_datum_kwds):  # PYCHOK Cartesian=Cartesian, datum=None
        '''Convert this point to C{Karney}-based cartesian (ECEF) coordinates.

           @kwarg Cartesian_datum_kwds: Optional L{Cartesian}, B{C{datum}}
                  and other keyword arguments, ignored if B{C{Cartesian=None}}.
                  Use B{C{Cartesian=...}} to override this L{Cartesian} class
                  or set B{C{Cartesian=None}}.

           @return: The cartesian (ECEF) coordinates (L{Cartesian}) or if
                    B{C{Cartesian}} is C{None}, an L{Ecef9Tuple}C{(x, y, z,
                    lat, lon, height, C, M, datum)} with C{C} and C{M} if
                    available.

           @raise TypeError: Invalid B{C{Cartesian}}, B{C{datum}} or other
                             B{C{Cartesian_datum_kwds}}.
        '''
        kwds = _xkwds(Cartesian_datum_kwds, Cartesian=Cartesian,
                                                datum=self.datum)
        return LatLonEllipsoidalBase.toCartesian(self, **kwds)

    def _direct(self, distance, bearing, LL, height):
        '''(INTERNAL) I{Karney}'s C{Direct} method.

           @return: A L{Destination2Tuple}C{(destination, final)} or
                    a L{Destination3Tuple}C{(lat, lon, final)} if
                    B{C{LL}} is C{None}.
        '''
        g = self.datum.ellipsoid.geodesic
        r = g.Direct3(self.lat, self.lon, bearing, distance)
        if LL:
            h = self.height if height is None else height
            d = LL(wrap90(r.lat), wrap180(r.lon), height=h, datum=self.datum)
            r = Destination2Tuple(self._xnamed(d), wrap360(r.final))
        return r

    def _inverse(self, other, wrap):
        '''(INTERNAL) I{Karney}'s C{Inverse} method.

           @return: A L{Distance3Tuple}C{(distance, initial, final)}.

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: If this and the B{C{other}} point's
                              L{Datum} ellipsoids are not compatible.
        '''
        g = self.ellipsoids(other).geodesic
        _, lon = unroll180(self.lon, other.lon, wrap=wrap)
        return g.Inverse3(self.lat, self.lon, other.lat, lon)


def _geodesic(datum, points, closed, line, wrap):
    # Compute the area or perimeter of a polygon,
    # using the GeographicLib package, iff installed
    g = datum.ellipsoid.geodesic

    if not wrap:  # capability LONG_UNROLL can't be off
        raise _ValueError(wrap=wrap)

    _, points = points2(points, closed=closed)  # base=LatLonEllipsoidalBase(0, 0)

    g = g.Polygon(line)

    # note, lon deltas are unrolled, by default
    for p in points:
        g.AddPoint(p.lat, p.lon)
    if closed and line:
        p = points[0]
        g.AddPoint(p.lat, p.lon)

    # g.Compute returns (number_of_points, perimeter, signed area)
    return g.Compute(False, True)[1 if line else 2]


def areaOf(points, datum=Datums.WGS84, wrap=True):
    '''Compute the area of a (n ellipsoidal) polygon.

       @arg points: The polygon points (L{LatLon}[]).
       @kwarg datum: Optional datum (L{Datum}).
       @kwarg wrap: Wrap and unroll longitudes (C{bool}).

       @return: Area (C{meter}, same as units of the B{C{datum}}
                ellipsoid, squared).

       @raise ImportError: Package U{GeographicLib
              <https://PyPI.org/project/geographiclib>} missing.

       @raise TypeError: Some B{C{points}} are not L{LatLon}.

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise ValueError: Invalid B{C{wrap}}, longitudes not
                          wrapped, unrolled.

       @note: This function requires installation of the U{GeographicLib
              <https://PyPI.org/project/geographiclib>} package.

       @see: L{pygeodesy.areaOf}, L{sphericalNvector.areaOf} and
             L{sphericalTrigonometry.areaOf}.
    '''
    return abs(_geodesic(datum, points, True, False, wrap))


def intersections2(center1, rad1, center2, rad2, height=None, wrap=False,
                   equidistant=None, tol=_TOL_M, LatLon=LatLon, **LatLon_kwds):
    '''Iteratively compute the intersection points of two circles each defined
       by an (ellipsoidal) center point and a radius.

       @arg center1: Center of the first circle (L{LatLon}).
       @arg rad1: Radius of the second circle (C{meter}).
       @arg center2: Center of the second circle (L{LatLon}).
       @arg rad2: Radius of the second circle (C{meter}).
       @kwarg height: Optional height for the intersection points,
                      overriding the "radical height" at the "radical
                      line" between both centers (C{meter}).
       @kwarg wrap: Wrap and unroll longitudes (C{bool}).
       @kwarg equidistant: An azimuthal equidistant projection class
                           (L{Equidistant} or L{equidistant})
                           or C{None} for L{EquidistantKarney}.
       @kwarg tol: Convergence tolerance (C{meter}).
       @kwarg LatLon: Optional class to return the intersection points
                      (L{LatLon}) or C{None}.
       @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                           arguments, ignored if B{C{LatLon=None}}.

       @return: 2-Tuple of the intersection points, each a B{C{LatLon}}
                instance or L{LatLon4Tuple}C{(lat, lon, height, datum)}
                if B{C{LatLon}} is C{None}.  For abutting circles, the
                intersection points are the same instance.

       @raise ImportError: Package U{geographiclib
                           <https://PyPI.org/project/geographiclib>}
                           not installed or not found.

       @raise IntersectionError: Concentric, antipodal, invalid or
                                 non-intersecting circles or no
                                 convergence for the B{C{tol}}.

       @raise TypeError: If B{C{center1}} or B{C{center2}} not ellipsoidal.

       @raise UnitError: Invalid B{C{rad1}}, B{C{rad2}} or B{C{height}}.

       @see: U{The B{ellipsoidal} case<https://GIS.StackExchange.com/questions/48937/
             calculating-intersection-of-two-circles>}, U{Karney's paper
             <https://ArXiv.org/pdf/1102.1215.pdf>}, pp 20-21, section 14 I{Maritime Boundaries},
             U{circle-circle<https://MathWorld.Wolfram.com/Circle-CircleIntersection.html>} and
             U{sphere-sphere<https://MathWorld.Wolfram.com/Sphere-SphereIntersection.html>}
             intersections.
    '''
    from pygeodesy.azimuthal import EquidistantKarney
    E = EquidistantKarney if equidistant is None else equidistant
    return _intersections2(center1, rad1, center2, rad2, height=height, wrap=wrap,
                                 equidistant=E, tol=tol, LatLon=LatLon, **LatLon_kwds)


def isclockwise(points, datum=Datums.WGS84, wrap=True):
    '''Determine the direction of a path or polygon.

       @arg points: The path or polygon points (C{LatLon}[]).
       @kwarg datum: Optional datum (L{Datum}).
       @kwarg wrap: Wrap and unroll longitudes (C{bool}).

       @return: C{True} if B{C{points}} are clockwise, C{False} otherwise.

       @raise TypeError: Some B{C{points}} are not C{LatLon}.

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise ValueError: The B{C{points}} enclose a pole or zero
                          area.

       @note: This function requires installation of the U{GeographicLib
              <https://PyPI.org/project/geographiclib>} package.

       @see: L{pygeodesy.isclockwise}.
    '''
    a = _geodesic(datum, points, True, False, wrap)
    if a > 0:
        return True
    elif a < 0:
        return False
    raise _areaError(points)


def perimeterOf(points, closed=False, datum=Datums.WGS84, wrap=True):
    '''Compute the perimeter of a (n ellipsoidal) polygon.

       @arg points: The polygon points (L{LatLon}[]).
       @kwarg closed: Optionally, close the polygon (C{bool}).
       @kwarg datum: Optional datum (L{Datum}).
       @kwarg wrap: Wrap and unroll longitudes (C{bool}).

       @return: Perimeter (C{meter}, same as units of the B{C{datum}}
                ellipsoid).

       @raise ImportError: Package U{GeographicLib
              <https://PyPI.org/project/geographiclib>} missing.

       @raise TypeError: Some B{C{points}} are not L{LatLon}.

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise ValueError: Invalid B{C{wrap}}, longitudes not
                          wrapped, unrolled.

       @note: This function requires installation of the U{GeographicLib
              <https://PyPI.org/project/geographiclib>} package.

       @see: L{pygeodesy.perimeterOf} and L{sphericalTrigonometry.perimeterOf}.
    '''
    return _geodesic(datum, points, closed, True, wrap)


__all__ += _ALL_OTHER(Cartesian, LatLon,  # classes
                      areaOf, intersections2, isclockwise, ispolar,  # functions
                      perimeterOf)

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
