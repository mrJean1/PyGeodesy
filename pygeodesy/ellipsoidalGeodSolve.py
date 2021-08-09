
# -*- coding: utf-8 -*-

u'''Exact ellipsoidal geodesy, intended I{for testing purposes only}.

Ellipsoidal geodetic (lat-/longitude) L{LatLon} and geocentric
(ECEF) L{Cartesian} classes and functions L{areaOf}, L{intersections2},
L{isclockwise}, L{nearestOn} and L{perimeterOf} based on module
L{geodsolve}, a wrapper invoking I{Karney}'s U{GeodSolve
<https://GeographicLib.SourceForge.io/html/GeodSolve.1.html>} utility.
'''

from pygeodesy.datums import _WGS84
from pygeodesy.ellipsoidalBase import CartesianEllipsoidalBase
from pygeodesy.ellipsoidalBaseDI import LatLonEllipsoidalBaseDI, \
                                       _intersection3, _intersections2, \
                                       _nearestOn
# from pygeodesy.errors import _xkwds  # from .karney
from pygeodesy.karney import _polygon, Property_RO, _TOL_M, _xkwds
from pygeodesy.lazily import _ALL_LAZY, _ALL_OTHER
from pygeodesy.points import _areaError, ispolar  # PYCHOK exported
# from pygeodesy.props import Property_RO  # from .karney
# from pygeodesy.units import _1mm as _TOL_M  # from .karney

__all__ = _ALL_LAZY.ellipsoidalGeodSolve
__version__ = '21.08.07'


class Cartesian(CartesianEllipsoidalBase):
    '''Extended to convert exact L{Cartesian} to exact L{LatLon} points.
    '''

    def toLatLon(self, **LatLon_datum_kwds):  # PYCHOK LatLon=LatLon, datum=None
        '''Convert this cartesian point to an exact geodetic point.

           @kwarg LatLon_datum_kwds: Optional L{LatLon}, B{C{datum}} and
                  other keyword arguments, ignored if C{B{LatLon} is None}.
                  Use C{B{LatLon}=...} to override this L{LatLon} class
                  or specify C{B{LatLon} is None}.

           @return: The geodetic point (L{LatLon}) or if B{C{LatLon}}
                    is C{None}, an L{Ecef9Tuple}C{(x, y, z, lat, lon,
                    height, C, M, datum)} with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{LatLon}}, B{C{datum}} or other
                             B{C{LatLon_datum_kwds}}.
        '''
        kwds = _xkwds(LatLon_datum_kwds, LatLon=LatLon, datum=self.datum)
        return CartesianEllipsoidalBase.toLatLon(self, **kwds)


class LatLon(LatLonEllipsoidalBaseDI):
    '''An ellipsoidal L{LatLon} like L{ellipsoidalKarney.LatLon} but using (exact)
       geodesic I{wrapper} L{GeodesicSolve} to compute the geodesic distance,
       initial and final bearing (azimuths) between two given points or the
       destination point given a start point and an (initial) bearing.
    '''

    @Property_RO
    def Equidistant(self):
        '''Get the prefered azimuthal equidistant projection I{class} (L{EquidistantGeodSolve}).
        '''
        from pygeodesy.azimuthal import EquidistantGeodSolve
        return EquidistantGeodSolve

    @Property_RO
    def geodesicx(self):
        '''Get this C{LatLon}'s (exact) geodesic (L{GeodesicSolve}).
        '''
        return self.datum.ellipsoid.geodsolve

    geodesic = geodesicx  # for C{._Direct} and C{._Inverse}

    def toCartesian(self, **Cartesian_datum_kwds):  # PYCHOK Cartesian=Cartesian, datum=None
        '''Convert this point to exact cartesian (ECEF) coordinates.

           @kwarg Cartesian_datum_kwds: Optional L{Cartesian}, B{C{datum}}
                  and other keyword arguments, ignored if C{B{Cartesian} is None}.
                  Use C{B{Cartesian}=...} to override this L{Cartesian} class
                  or set C{B{Cartesian} is None}.

           @return: The cartesian (ECEF) coordinates (L{Cartesian}) or if
                    B{C{Cartesian}} is C{None}, an L{Ecef9Tuple}C{(x, y, z,
                    lat, lon, height, C, M, datum)} with C{C} and C{M} if
                    available.

           @raise TypeError: Invalid B{C{Cartesian}}, B{C{datum}} or other
                             B{C{Cartesian_datum_kwds}}.
        '''
        kwds = _xkwds(Cartesian_datum_kwds, Cartesian=Cartesian, datum=self.datum)
        return LatLonEllipsoidalBaseDI.toCartesian(self, **kwds)


def areaOf(points, datum=_WGS84, wrap=True):
    '''Compute the area of an (ellipsoidal) polygon.

       @arg points: The polygon points (L{LatLon}[]).
       @kwarg datum: Optional datum (L{Datum}).
       @kwarg wrap: Wrap and unroll longitudes (C{bool}).

       @return: Area (C{meter}, same as units of the
                B{C{datum}}'s ellipsoid axes, I{squared}).

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise TypeError: Some B{C{points}} are not L{LatLon}.

       @raise ValueError: Invalid C{B{wrap}=False}, unwrapped,
                          unrolled longitudes not supported.

       @see: L{pygeodesy.areaOf}, L{ellipsoidalExact.areaOf},
             L{ellipsoidalKarney.areaOf}, L{sphericalNvector.areaOf}
             and L{sphericalTrigonometry.areaOf}.
    '''
    return abs(_polygon(datum.ellipsoid.geodsolve, points, True, False, wrap))


def intersection3(start1, end1, start2, end2, height=None, wrap=True,
                  equidistant=None, tol=_TOL_M, LatLon=LatLon, **LatLon_kwds):
    '''Interatively compute the intersection point of two paths,
       each defined by an (ellipsoidal) start and end point.

       @arg start1: Start point of the first path (L{LatLon}).
       @arg end1: End point of the first path (L{LatLon}).
       @arg start2: Start point of the second path (L{LatLon}).
       @arg end2: End point of the second path (L{LatLon}).
       @kwarg height: Optional height at the intersection (C{meter},
                      conventionally) or C{None} for the mean height.
       @kwarg wrap: Wrap and unroll longitudes (C{bool}).
       @kwarg equidistant: An azimuthal equidistant projection (I{class}
                           or function L{equidistant}) or C{None} for
                           the preferred C{B{start1}.Equidistant}.
       @kwarg tol: Tolerance for skew line distance and length and for
                   convergence (C{meter}, conventionally).
       @kwarg LatLon: Optional class to return the intersection point
                      (L{LatLon}) or C{None}.
       @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                           arguments, ignored if C{B{LatLon} is None}.

       @return: An L{Intersection3Tuple}C{(point, outside1, outside2)}
                with C{point} a B{C{LatLon}} or if C{B{LatLon} is None},
                a L{LatLon4Tuple}C{(lat, lon, height, datum)}.

       @raise IntersectionError: Skew, colinear, parallel or otherwise
                                 non-intersecting paths or no convergence
                                 for the B{C{tol}}.

       @raise TypeError: Invalid or non-ellipsoidal B{C{start1}}, B{C{end1}},
                         B{C{start2}} or B{C{end2}} or invalide B{C{equidistant}}.
    '''
    return _intersection3(start1, end1, start2, end2, height=height, wrap=wrap,
                          equidistant=equidistant, tol=tol, LatLon=LatLon, **LatLon_kwds)


def intersections2(center1, radius1, center2, radius2, height=None, wrap=True,
                   equidistant=None, tol=_TOL_M, LatLon=LatLon, **LatLon_kwds):
    '''Iteratively compute the intersection points of two circles each defined
       by an (ellipsoidal) center point and a radius.

       @arg center1: Center of the first circle (L{LatLon}).
       @arg radius1: Radius of the first circle (C{meter}, conventionally).
       @arg center2: Center of the second circle (L{LatLon}).
       @arg radius2: Radius of the second circle (C{meter}, same units as
                     B{C{radius1}}).
       @kwarg height: Optional height for the intersection points (C{meter},
                      conventionally) or C{None} for the I{"radical height"}
                      at the I{radical line} between both centers.
       @kwarg wrap: Wrap and unroll longitudes (C{bool}).
       @kwarg equidistant: An azimuthal equidistant projection (I{class}
                           or function L{equidistant}) or C{None} for
                           the preferred C{B{center1}.Equidistant}.
       @kwarg tol: Convergence tolerance (C{meter}, same units as
                   B{C{radius1}} and B{C{radius2}}).
       @kwarg LatLon: Optional class to return the intersection points
                      (L{LatLon}) or C{None}.
       @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                           arguments, ignored if C{B{LatLon} is None}.

       @return: 2-Tuple of the intersection points, each a B{C{LatLon}}
                instance or L{LatLon4Tuple}C{(lat, lon, height, datum)}
                if C{B{LatLon} is None}.  For abutting circles, both
                points are the same instance, aka I{radical center}.

       @raise IntersectionError: Concentric, antipodal, invalid or
                                 non-intersecting circles or no
                                 convergence for the B{C{tol}}.

       @raise TypeError: Invalid or non-ellipsoidal B{C{center1}} or B{C{center2}}
                         or invalid B{C{equidistant}}.

       @raise UnitError: Invalid B{C{radius1}}, B{C{radius2}} or B{C{height}}.

       @see: U{The B{ellipsoidal} case<https://GIS.StackExchange.com/questions/48937/
             calculating-intersection-of-two-circles>}, U{Karney's paper
             <https://ArXiv.org/pdf/1102.1215.pdf>}, pp 20-21, section B{14. MARITIME BOUNDARIES},
             U{circle-circle<https://MathWorld.Wolfram.com/Circle-CircleIntersection.html>} and
             U{sphere-sphere<https://MathWorld.Wolfram.com/Sphere-SphereIntersection.html>}
             intersections.
    '''
    return _intersections2(center1, radius1, center2, radius2, height=height, wrap=wrap,
                           equidistant=equidistant, tol=tol, LatLon=LatLon, **LatLon_kwds)


def isclockwise(points, datum=_WGS84, wrap=True):
    '''Determine the direction of a path or polygon.

       @arg points: The path or polygon points (C{LatLon}[]).
       @kwarg datum: Optional datum (L{Datum}).
       @kwarg wrap: Wrap and unroll longitudes (C{bool}).

       @return: C{True} if B{C{points}} are clockwise, C{False} otherwise.

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise TypeError: Some B{C{points}} are not C{LatLon}.

       @raise ValueError: The B{C{points}} enclose a pole or zero
                          area.

       @see: L{pygeodesy.isclockwise}.
    '''
    a = _polygon(datum.ellipsoid.geodsolve, points, True, False, wrap)
    if a > 0:
        return True
    elif a < 0:
        return False
    raise _areaError(points)


def nearestOn(point, point1, point2, within=True, height=None, wrap=False,
              equidistant=None, tol=_TOL_M, LatLon=LatLon, **LatLon_kwds):
    '''Iteratively locate the closest point on the arc between two
       other (ellipsoidal) points.

       @arg point: Reference point (C{LatLon}).
       @arg point1: Start point of the arc (C{LatLon}).
       @arg point2: End point of the arc (C{LatLon}).
       @kwarg within: If C{True} return the closest point I{between}
                      B{C{point1}} and B{C{point2}}, otherwise the
                      closest point elsewhere on the arc (C{bool}).
       @kwarg height: Optional height for the closest point (C{meter},
                      conventionally) or C{None} or C{False} for the
                      interpolated height.  If C{False}, the distance
                      between points takes height into account.
       @kwarg wrap: Wrap and unroll longitudes (C{bool}).
       @kwarg equidistant: An azimuthal equidistant projection (I{class}
                           or function L{equidistant}) or C{None} for
                           the preferred C{B{point}.Equidistant}.
       @kwarg tol: Convergence tolerance (C{meter}).
       @kwarg LatLon: Optional class to return the closest point
                      (L{LatLon}) or C{None}.
       @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                           arguments, ignored if C{B{LatLon} is None}.

       @return: Closest point, a B{C{LatLon}} instance or if C{B{LatLon}
                is None}, a L{LatLon4Tuple}C{(lat, lon, height, datum)}.

       @raise TypeError: Invalid or non-ellipsoidal B{C{point}}, B{C{point1}}
                         or B{C{point2}} or invalid B{C{equidistant}}.

       @raise ValueError: No convergence for the B{C{tol}}.
    '''
    return _nearestOn(point, point1, point2, within=within, height=height, wrap=wrap,
                      equidistant=equidistant, tol=tol, LatLon=LatLon, **LatLon_kwds)


def perimeterOf(points, closed=False, datum=_WGS84, wrap=True):
    '''Compute the perimeter of an (ellipsoidal) polygon.

       @arg points: The polygon points (L{LatLon}[]).
       @kwarg closed: Optionally, close the polygon (C{bool}).
       @kwarg datum: Optional datum (L{Datum}).
       @kwarg wrap: Wrap and unroll longitudes (C{bool}).

       @return: Perimeter (C{meter}, same as units of the
                B{C{datum}}'s ellipsoid axes).

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise TypeError: Some B{C{points}} are not L{LatLon}.

       @raise ValueError: Invalid C{B{wrap}=False}, unwrapped,
                          unrolled longitudes not supported.

       @see: L{pygeodesy.perimeterOf}, L{ellipsoidalExact.perimeterOf},
             L{ellipsoidalKarney.perimeterOf}, L{sphericalNvector.perimeterOf}
             and L{sphericalTrigonometry.perimeterOf}.
    '''
    return _polygon(datum.ellipsoid.geodsolve, points, closed, True, wrap)


__all__ += _ALL_OTHER(Cartesian, LatLon,  # classes
                      areaOf,  # functions
                      intersection3, intersections2, isclockwise, ispolar,
                      nearestOn, perimeterOf)

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
