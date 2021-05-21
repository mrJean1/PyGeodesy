
# -*- coding: utf-8 -*-

u'''Exact ellipsoidal geodesy.

Ellipsoidal geodetic (lat-/longitude) L{LatLon} and geocentric
(ECEF) L{Cartesian} classes and functions L{areaOf}, L{intersections2},
L{isclockwise}, L{nearestOn} and L{perimeterOf} based on classes
L{GeodesicExact}, L{GeodesicAreaExact} and L{GeodesicLineExact}.
'''

from pygeodesy.datums import _WGS84
from pygeodesy.ellipsoidalBase import _intermediateTo, _intersections2, \
                                       CartesianEllipsoidalBase, \
                                       LatLonEllipsoidalBase, _nearestOn
from pygeodesy.errors import _xellipsoidal, _xkwds
from pygeodesy.karney import _Equidistant, _polygon
from pygeodesy.lazily import _ALL_LAZY, _ALL_OTHER
from pygeodesy.namedTuples import Bearing2Tuple, Destination2Tuple
from pygeodesy.points import _areaError, ispolar  # PYCHOK exported
from pygeodesy.props import Property_RO
from pygeodesy.units import _1mm as _TOL_M
from pygeodesy.utily import unroll180, wrap90, wrap180, wrap360

__all__ = _ALL_LAZY.ellipsoidalExact
__version__ = '21.05.18'


class Cartesian(CartesianEllipsoidalBase):
    '''Extended to convert exact L{Cartesian} to exact L{LatLon} points.
    '''

    def toLatLon(self, **LatLon_datum_kwds):  # PYCHOK LatLon=LatLon, datum=None
        '''Convert this cartesian point to an exact geodetic point.

           @kwarg LatLon_datum_kwds: Optional L{LatLon}, B{C{datum}} and
                  other keyword arguments, ignored if C{B{LatLon}=None}.
                  Use C{B{LatLon}=...} to override this L{LatLon} class
                  or specify C{B{LatLon}=None}.

           @return: The geodetic point (L{LatLon}) or if B{C{LatLon}}
                    is C{None}, an L{Ecef9Tuple}C{(x, y, z, lat, lon,
                    height, C, M, datum)} with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{LatLon}}, B{C{datum}} or other
                             B{C{LatLon_datum_kwds}}.
        '''
        kwds = _xkwds(LatLon_datum_kwds, LatLon=LatLon, datum=self.datum)
        return CartesianEllipsoidalBase.toLatLon(self, **kwds)


class LatLon(LatLonEllipsoidalBase):
    '''An ellipsoidal L{LatLon} like L{ellipsoidalKarney.LatLon} but using
       exact geodesic classes L{GeodesicExact} and L{GeodesicLineExact} to
       compute the geodesic distance, initial and final bearing (azimuths)
       between two given points or the destination point given a start point
       and an (initial) bearing.
    '''

    def bearingTo2(self, other, wrap=False):
        '''Compute the initial and final bearing (forward and reverse
           azimuth) from this to an other point, using I{Karney}'s
           C{Inverse} method.  See methods L{initialBearingTo} and
           L{finalBearingTo} for more details.

           @arg other: The other point (L{LatLon}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: A L{Bearing2Tuple}C{(initial, final)}.

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: If this and the B{C{other}} point's L{Datum}
                              ellipsoids are not compatible.
        '''
        r = self._inverse(other, wrap)
        return Bearing2Tuple(r.initial, r.final, name=self.name)

    def destination(self, distance, bearing, height=None):
        '''Compute the destination point after having travelled for
           the given distance from this point along a geodesic given
           by an initial bearing, using the C{GeodesicExact.Direct}
           method.  See method L{destination2} for more details.

           @arg distance: Distance (C{meter}).
           @arg bearing: Initial bearing in (compass C{degrees360}).
           @kwarg height: Optional height, overriding the default
                          height (C{meter}, same units as C{distance}).

           @return: The destination point (L{LatLon}).
        '''
        return self._direct(distance, bearing, self.classof, height).destination

    def destination2(self, distance, bearing, height=None):
        '''Compute the destination point and the final bearing (reverse
           azimuth) after having travelled for the given distance from
           this point along a geodesic given by an initial bearing,
           using the C{GeodesicExact.Direct} method.

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
        '''
        r = self._direct(distance, bearing, self.classof, height)
        return self._xnamed(r)

    def distanceTo(self, other, wrap=False, **unused):  # ignore radius=R_M
        '''Compute the distance between this and an other point
           along a geodesic, using I{Karney}'s C{Inverse} method.
           See method L{distanceTo3} for more details.

           @arg other: The other point (L{LatLon}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Distance (C{meter}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: If this and the B{C{other}} point's L{Datum}
                              ellipsoids are not compatible.
        '''
        return self._inverse(other, wrap).distance

    def distanceTo3(self, other, wrap=False):
        '''Compute the distance, the initial and final bearing along
           a geodesic between this and an other point, using the
           C{GeodesicExact.Inverse} method.

           The distance is in the same units as this point's datum axes,
           conventionally meter.  The distance is measured on the surface
           of the ellipsoid, ignoring this point's height.

           The initial and final bearing (forward and reverse azimuth)
           are in compass C{degrees360} from North.

           @arg other: Destination point (L{LatLon}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: A L{Distance3Tuple}C{(distance, initial, final)}.

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: If this and the B{C{other}} point's L{Datum}
                              ellipsoids are not compatible.
        '''
        return self._xnamed(self._inverse(other, wrap))

    def finalBearingOn(self, distance, bearing):
        '''Compute the final bearing (reverse azimuth) after having
           travelled for the given distance along a geodesic given by
           an initial bearing from this point, using the
           C{GeodesicExact.Direct} method.  See method L{destination2}
           for more details.

           @arg distance: Distance (C{meter}).
           @arg bearing: Initial bearing (compass C{degrees360}).

           @return: Final bearing (compass C{degrees360}).
        '''
        return self._direct(distance, bearing, None, None).final

    def finalBearingTo(self, other, wrap=False):
        '''Compute the final bearing (reverse azimuth) after having
           travelled along a geodesic from this point to an other
           point, using the C{GeodesicExact.Inverse} method.  See
           method L{distanceTo3} for more details.

           @arg other: The other point (L{LatLon}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Final bearing (compass C{degrees360}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: If this and the B{C{other}} point's L{Datum}
                              ellipsoids are not compatible.
        '''
        return self._inverse(other, wrap).final

    @Property_RO
    def geodesicx(self):
        '''Get this C{LatLon}'s exact geodesic (L{GeodesicExact}).
        '''
        return self.datum.ellipsoid.geodesicx

    geodesic = geodesicx  # for convenience

    def initialBearingTo(self, other, wrap=False):
        '''Compute the initial bearing (forward azimuth) to travel
           along a geodesic from this point to an other point,
           using the C{GeodesicExact.Inverse} method.  See method
           L{distanceTo3} for more details.

           @arg other: The other point (L{LatLon}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Initial bearing (compass C{degrees360}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: If this and the B{C{other}} point's L{Datum}
                              ellipsoids are not compatible.
        '''
        return self._inverse(other, wrap).initial

    def intermediateTo(self, other, fraction, height=None, wrap=False):
        '''Return the point at given fraction along the geodesic between
           this and an other point, using the C{GeodesicExact.Direct} and
           C{GeodesicExact.Inverse} methods.

           @arg other: The other point (L{LatLon}).
           @arg fraction: Fraction between both points ranging from
                          0, meaning this to 1, the other point (C{float}).
           @kwarg height: Optional height, overriding the fractional
                          height (C{meter}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Intermediate point (L{LatLon}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise UnitError: Invalid B{C{fraction}} or B{C{height}}.

           @raise ValueError: If this and the B{C{other}} point's L{Datum}
                              ellipsoids are not compatible.

           @see: Methods L{distanceTo3} and L{destination}.
        '''
        return _intermediateTo(self, other, fraction, height, wrap)

    def intersections2(self, radius1, other, radius2, height=None, wrap=True,  # PYCHOK expected
                                                      tol=_TOL_M):
        '''Compute the intersection points of two circles each defined
           by a center point and a radius.

           @arg radius1: Radius of the this circle (C{meter}, conventionally).
           @arg other: Center of the other circle (C{LatLon}).
           @arg radius2: Radius of the other circle (C{meter}, same units as
                         B{C{radius1}}).
           @kwarg height: Optional height for the intersection points,
                          overriding the "radical height" at the "radical
                          line" between both centers (C{meter}, conventionally)
                          or C{None}.
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).
           @kwarg tol: Convergence tolerance (C{meter}, same units as B{C{radius1}}
                       and B{C{radius2}}).

           @return: 2-Tuple of the intersection points, each a C{LatLon}
                    instance.  For abutting circles, both intersection
                    points are the same instance.

           @raise IntersectionError: Concentric, antipodal, invalid or
                                     non-intersecting circles or no
                                     convergence.

           @raise TypeError: Invalid B{C{other}} or B{C{equidistant}}.

           @raise UnitError: Invalid B{C{radius1}}, B{C{radius2}} or B{C{height}}.

           @see: U{The B{ellipsoidal} case<https://GIS.StackExchange.com/questions/48937/
                 calculating-intersection-of-two-circles>}, U{Karney's paper
                 <https://ArXiv.org/pdf/1102.1215.pdf>}, pp 20-21, section B{14. MARITIME BOUNDARIES},
                 U{circle-circle<https://MathWorld.Wolfram.com/Circle-CircleIntersection.html>} and
                 U{sphere-sphere<https://MathWorld.Wolfram.com/Sphere-SphereIntersection.html>}
                 intersections.
        '''
        self.others(other)
        return intersections2(self, radius1, other, radius2, height=height, wrap=wrap,
                                    tol=tol, LatLon=self.classof, datum=self.datum)

    def nearestOn(self, point1, point2, within=True, height=None, wrap=False,  # PYCHOK expected
                                                     tol=_TOL_M):
        '''Locate the closest point on the arc between two other points
           and this point.

           @arg point1: Start point of the arc (C{LatLon}).
           @arg point2: End point of the arc (C{LatLon}).
           @kwarg within: If C{True} return the closest point I{between}
                          B{C{point1}} and B{C{point2}}, otherwise the
                          closest point elsewhere on the arc (C{bool}).
           @kwarg height: Optional height for the closest point (C{meter})
                          or C{None}.
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).
           @kwarg tol: Convergence tolerance (C{meter}).

           @return: Closest point (B{C{LatLon}}).

           @raise TypeError: Invalid or non-ellipsoidal B{C{point1}} or B{C{point2}}.
        '''
        # use .nearestOn to get C{azimuthal.EquidistantKarney} or ImportError
        return nearestOn(self, point1, point2, within=within, height=height, wrap=wrap,
                               tol=tol, LatLon=self.classof, datum=self.datum)

    def toCartesian(self, **Cartesian_datum_kwds):  # PYCHOK Cartesian=Cartesian, datum=None
        '''Convert this point to exact cartesian (ECEF) coordinates.

           @kwarg Cartesian_datum_kwds: Optional L{Cartesian}, B{C{datum}}
                  and other keyword arguments, ignored if C{B{Cartesian}=None}.
                  Use C{B{Cartesian}=...} to override this L{Cartesian} class
                  or set C{B{Cartesian}=None}.

           @return: The cartesian (ECEF) coordinates (L{Cartesian}) or if
                    B{C{Cartesian}} is C{None}, an L{Ecef9Tuple}C{(x, y, z,
                    lat, lon, height, C, M, datum)} with C{C} and C{M} if
                    available.

           @raise TypeError: Invalid B{C{Cartesian}}, B{C{datum}} or other
                             B{C{Cartesian_datum_kwds}}.
        '''
        kwds = _xkwds(Cartesian_datum_kwds, Cartesian=Cartesian, datum=self.datum)
        return LatLonEllipsoidalBase.toCartesian(self, **kwds)

    def _direct(self, distance, bearing, LL, height):
        '''(INTERNAL) C{GeodesicExact.Direct} method.

           @return: A L{Destination2Tuple}C{(destination, final)} or
                    a L{Destination3Tuple}C{(lat, lon, final)} if
                    B{C{LL}} is C{None}.
        '''
        gX = self.geodesicx
        r = gX.Direct3(self.lat, self.lon, bearing, distance)
        if LL:
            h = self.height if height is None else height
            d = LL(wrap90(r.lat), wrap180(r.lon), height=h, datum=self.datum)
            r = Destination2Tuple(self._xnamed(d), wrap360(r.final))
        return r

    def _inverse(self, other, wrap):
        '''(INTERNAL) C{GeodesicExact.Inverse} method.

           @return: A L{Distance3Tuple}C{(distance, initial, final)}.

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: If this and the B{C{other}} point's
                              L{Datum} ellipsoids are not compatible.
        '''
        gX = self.ellipsoids(other).geodesicx
        _, lon = unroll180(self.lon, other.lon, wrap=wrap)
        return gX.Inverse3(self.lat, self.lon, other.lat, lon)


def areaOf(points, datum=_WGS84, wrap=True):
    '''Compute the area of an (ellipsoidal) polygon.

       @arg points: The polygon points (L{LatLon}[]).
       @kwarg datum: Optional datum (L{Datum}).
       @kwarg wrap: Wrap and unroll longitudes (C{bool}).

       @return: Area (C{meter}, same as units of the B{C{datum}}
                ellipsoid, squared).

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise TypeError: Some B{C{points}} are not L{LatLon}.

       @raise ValueError: Invalid C{B{wrap}=False}, unwrapped,
                          unrolled longitudes not supported.

       @see: L{pygeodesy.areaOf}, L{sphericalNvector.areaOf} and
             L{sphericalTrigonometry.areaOf}.
    '''
    return abs(_polygon(datum.ellipsoid.geodesicx, points, True, False, wrap))


def intersections2(center1, radius1, center2, radius2, height=None, wrap=True,
                   equidistant=None, tol=_TOL_M, LatLon=LatLon, **LatLon_kwds):
    '''Iteratively compute the intersection points of two circles each defined
       by an (ellipsoidal) center point and a radius.

       @arg center1: Center of the first circle (L{LatLon}).
       @arg radius1: Radius of the first circle (C{meter}, conventionally).
       @arg center2: Center of the second circle (L{LatLon}).
       @arg radius2: Radius of the second circle (C{meter}, same units as
                     B{C{radius1}}).
       @kwarg height: Optional height for the intersection points,
                      overriding the "radical height" at the "radical
                      line" between both centers (C{meter}, conventionally)
                      or C{None}.
       @kwarg wrap: Wrap and unroll longitudes (C{bool}).
       @kwarg equidistant: An azimuthal equidistant projection (class
                           L{Equidistant} or function L{equidistant}) or
                           C{None} for L{EquidistantKarney}.
       @kwarg tol: Convergence tolerance (C{meter}, same units as
                   B{C{radius1}} and B{C{radius2}}).
       @kwarg LatLon: Optional class to return the intersection points
                      (L{LatLon}) or C{None}.
       @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                           arguments, ignored if C{B{LatLon}=None}.

       @return: 2-Tuple of the intersection points, each a B{C{LatLon}}
                instance or L{LatLon4Tuple}C{(lat, lon, height, datum)}
                if B{C{LatLon}} is C{None}.  For abutting circles, both
                intersection points are the same instance.

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
    E = _Equidistant(equidistant, exact=True)
    return _intersections2(center1, radius1, center2, radius2, height=height, wrap=wrap,
                                    equidistant=E, tol=tol, LatLon=LatLon, **LatLon_kwds)


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
    a = _polygon(datum.ellipsoid.geodesicx, points, True, False, wrap)
    if a > 0:
        return True
    elif a < 0:
        return False
    raise _areaError(points)


def nearestOn(point, point1, point2, within=True, height=None, wrap=False,
              equidistant=None, tol=_TOL_M, LatLon=LatLon, **LatLon_kwds):
    '''Iteratively locate the closest point on the arc between two other points.

       @arg point: Reference point (C{LatLon}).
       @arg point1: Start point of the arc (C{LatLon}).
       @arg point2: End point of the arc (C{LatLon}).
       @kwarg within: If C{True} return the closest point I{between}
                      B{C{point1}} and B{C{point2}}, otherwise the
                      closest point elsewhere on the arc (C{bool}).
       @kwarg height: Optional height for the closest point (C{meter})
                      or C{None}.
       @kwarg wrap: Wrap and unroll longitudes (C{bool}).
       @kwarg equidistant: An azimuthal equidistant projection (class
                           L{Equidistant} or function L{equidistant}) or
                           C{None} for L{EquidistantKarney}.
       @kwarg tol: Convergence tolerance (C{meter}).
       @kwarg LatLon: Optional class to return the closest point
                      (L{LatLon}) or C{None}.
       @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                           arguments, ignored if C{B{LatLon}=None}.

       @return: Closest point (B{C{LatLon}}).

       @raise TypeError: Invalid or non-ellipsoidal B{C{point}}, B{C{point1}}
                         or B{C{point2}} or invalid B{C{equidistant}}.
    '''
    p = _xellipsoidal(point=point)
    p1 = p.others(point1=point1)
    p2 = p.others(point2=point2)
    E = _Equidistant(equidistant, exact=True)
    return _nearestOn(p, p1, p2, within=within, height=height, wrap=wrap,
                      equidistant=E, tol=tol, LatLon=LatLon, **LatLon_kwds)


def perimeterOf(points, closed=False, datum=_WGS84, wrap=True):
    '''Compute the perimeter of an (ellipsoidal) polygon.

       @arg points: The polygon points (L{LatLon}[]).
       @kwarg closed: Optionally, close the polygon (C{bool}).
       @kwarg datum: Optional datum (L{Datum}).
       @kwarg wrap: Wrap and unroll longitudes (C{bool}).

       @return: Perimeter (C{meter}, same as units of the B{C{datum}}
                ellipsoid).

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise TypeError: Some B{C{points}} are not L{LatLon}.

       @raise ValueError: Invalid C{B{wrap}=False}, unwrapped,
                          unrolled longitudes not supported.

       @see: L{pygeodesy.perimeterOf} and L{sphericalTrigonometry.perimeterOf}.
    '''
    return _polygon(datum.ellipsoid.geodesicx, points, closed, True, wrap)


__all__ += _ALL_OTHER(Cartesian, LatLon,  # classes
                      areaOf, intersections2, isclockwise, ispolar,  # functions
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
