
# -*- coding: utf-8 -*-

u'''Trigonometric spherical geodetic (lat-/longitude) class L{LatLon}
and functions L{areaOf}, L{intersection}, L{isPoleEnclosedBy},
L{meanOf}, L{nearestOn2} and L{perimeterOf}.

Pure Python implementation of geodetic (lat-/longitude) methods using
spherical trigonometry, transcribed from JavaScript originals by
I{(C) Chris Veness 2011-2016} published under the same MIT Licence**, see
U{Latitude/Longitude<http://www.movable-type.co.uk/scripts/latlong.html>}.

@newfield example: Example, Examples
'''

from datum import R_M
from fmath import EPS, favg, fmean, fsum, map1
from sphericalBase import LatLonSphericalBase
from utils import PI2, PI_2, degrees90, degrees180, degrees360, \
                  equirectangular_, haversine_, iterNumpy2, \
                  radians, tan_2, unrollPI, wrap180, wrapPI
from vector3d import Vector3d, sumOf

from math import acos, asin, atan2, copysign, cos, hypot, sin, sqrt

# all public contants, classes and functions
__all__ = ('LatLon',  # classes
           'areaOf',  # functions
           'intersection', 'isPoleEnclosedBy',
           'meanOf',
           'nearestOn2',
           'perimeterOf')
__version__ = '18.03.02'


class LatLon(LatLonSphericalBase):
    '''New point on spherical model earth model.

       @example:

       >>> p = LatLon(52.205, 0.119)  # height=0
    '''

    _v3d = None  # cache Vector3d

    def _update(self, updated):
        '''(INTERNAL) Clear caches if updated.
        '''
        if updated:  # reset caches
            self._v3d = None
            LatLonSphericalBase._update(self, updated)

    def _trackDistanceTo3(self, start, end, radius, wrap):
        '''(INTERNAL) Helper for along-/crossTrackDistanceTo.
        '''
        self.others(start, name='start')
        self.others(end, name='end')

        r = float(radius)
        if r < EPS:
            raise ValueError('%s invalid: %r' % ('radius', radius))

        r = start.distanceTo(self, r, wrap=wrap) / r
        b = radians(start.initialBearingTo(self, wrap=wrap))
        e = radians(start.initialBearingTo(end, wrap=wrap))

        x = asin(sin(r) * sin(b - e))
        return r, x, e - b

    def alongTrackDistanceTo(self, start, end, radius=R_M, wrap=False):
        '''Compute the (signed) distance from the start to the closest
           point on the great circle path defined by a start and an
           end point.

           That is, if a perpendicular is drawn from this point to the
           great circle path, the along-track distance is the distance
           from the start point to the point where the perpendicular
           crosses the path.

           @param start: Start point of great circle path (L{LatLon}).
           @param end: End point of great circle path (L{LatLon}).
           @keyword radius: Optional, mean earth radius (meter).
           @keyword wrap: Wrap and unroll longitudes (bool).

           @return: Distance along the great circle path (positive if
                    after the start toward the end point of the path
                    or negative if before the start point).

           @raise TypeError: The I{start} or I{end} point is not L{LatLon}.

           @example:

           >>> p = LatLon(53.2611, -0.7972)

           >>> s = LatLon(53.3206, -1.7297)
           >>> e = LatLon(53.1887, 0.1334)
           >>> d = p.alongTrackDistanceTo(s, e)  # 62331.58
        '''
        r, x, b = self._trackDistanceTo3(start, end, radius, wrap)
        cx = cos(x)
        if abs(cx) > EPS:
            return copysign(acos(cos(r) / cx), cos(b)) * radius
        else:
            return 0.0

    def crossingParallels(self, other, lat, wrap=False):
        '''Return the pair of meridians at which a great circle defined
           by this and an other point crosses the given latitude.

           @param other: The other point defining great circle (L{LatLon}).
           @param lat: Latitude at the crossing (degrees).
           @keyword wrap: Wrap and unroll longitudes (bool).

           @return: 2-Tuple (lon1, lon2) in (degrees180) or None if the
                    great circle doesn't reach the given I{lat}.
        '''
        self.others(other)

        a1, b1 = self.to2ab()
        a2, b2 = other.to2ab()

        a = radians(lat)
        db, b2 = unrollPI(b1, b2, wrap=wrap)

        ca, ca1, ca2, cdb = map1(cos, a, a1, a2, db)
        sa, sa1, sa2, sdb = map1(sin, a, a1, a2, db)

        x = sa1 * ca2 * ca * sdb
        y = sa1 * ca2 * ca * cdb - ca1 * sa2 * ca
        z = ca1 * ca2 * sa * sdb

        h = hypot(x, y)
        if h < EPS or abs(z) > h:
            return None  # great circle doesn't reach latitude

        m = atan2(-y, x) + b1  # longitude at max latitude
        d = acos(z / h)  # delta longitude to intersections

        return degrees180(m - d), degrees180(m + d)

    def crossTrackDistanceTo(self, start, end, radius=R_M, wrap=False):
        '''Compute the (signed) distance from this point to the great
           circle defined by a start and an end point.

           @param start: Start point of great circle path (L{LatLon}).
           @param end: End point of great circle path (L{LatLon}).
           @keyword radius: Optional, mean earth radius (meter).
           @keyword wrap: Wrap and unroll longitudes (bool).

           @return: Distance to great circle (negative if to the
                    left or positive if to the right of the path).

           @raise TypeError: The I{start} or I{end} point is not L{LatLon}.

           @example:

           >>> p = LatLon(53.2611, -0.7972)

           >>> s = LatLon(53.3206, -1.7297)
           >>> e = LatLon(53.1887, 0.1334)
           >>> d = p.crossTrackDistanceTo(s, e)  # -307.5
        '''
        _, x, _ = self._trackDistanceTo3(start, end, radius, wrap)
        return x * radius

    def destination(self, distance, bearing, radius=R_M, height=None):
        '''Locate the destination from this point after having
           travelled the given distance on the given initial bearing.

           @param distance: Distance travelled (same units as radius).
           @param bearing: Bearing from this point (compass degrees).
           @keyword radius: Optional, mean earth radius (meter).
           @keyword height: Optional height at destination (meter).

           @return: Destination point (L{LatLon}).

           @example:

           >>> p1 = LatLon(51.4778, -0.0015)
           >>> p2 = p1.destination(7794, 300.7)
           >>> p2.toStr()  # '51.5135°N, 000.0983°W'

           @JSname: I{destinationPoint}.
        '''
        a, b = self.to2ab()

        r = float(distance) / float(radius)  # angular distance in radians
        t = radians(bearing)

        a, b = _destination2(a, b, r, t)
        h = self.height if height is None else height
        return self.classof(a, b, height=h)

    def distanceTo(self, other, radius=R_M, wrap=False):
        '''Compute the distance from this to an other point.

           @param other: The other point (L{LatLon}).
           @keyword radius: Optional, mean earth radius (meter).
           @keyword wrap: Wrap and unroll longitudes (bool).

           @return: Distance between this and the I{other} point
                    (in the same units as radius).

           @raise TypeError: The I{other} point is not L{LatLon}.

           @example:

           >>> p1 = LatLon(52.205, 0.119)
           >>> p2 = LatLon(48.857, 2.351);
           >>> d = p1.distanceTo(p2)  # 404300
        '''
        self.others(other)

        a1, b1 = self.to2ab()
        a2, b2 = other.to2ab()

        db, b2 = unrollPI(b1, b2, wrap=wrap)
        r = haversine_(a2, a1, db)
        return r * float(radius)

    def greatCircle(self, bearing):
        '''Compute the vector normal to great circle obtained by heading
           on the given initial bearing from this point.

           Direction of vector is such that initial bearing vector
           b = c × n, where n is an n-vector representing this point.

           @param bearing: Bearing from this point (compass degrees).

           @return: Vector representing great circle (L{Vector3d}).

           @example:

           >>> p = LatLon(53.3206, -1.7297)
           >>> g = p.greatCircle(96.0)
           >>> g.toStr()  # (-0.794, 0.129, 0.594)
        '''
        a, b = self.to2ab()
        t = radians(bearing)

        ca, cb, ct = map1(cos, a, b, t)
        sa, sb, st = map1(sin, a, b, t)

        return Vector3d(sb * ct - cb * sa * st,
                       -cb * ct - sb * sa * st,
                        ca * st)  # XXX .unit()?

    def initialBearingTo(self, other, wrap=False):
        '''Compute the initial bearing (aka forward azimuth) from
           this to an other point.

           @param other: The other point (L{LatLon}).
           @keyword wrap: Wrap and unroll longitudes (bool).

           @return: Initial bearing (compass degrees).

           @raise TypeError: The I{other} point is not L{LatLon}.

           @example:

           >>> p1 = LatLon(52.205, 0.119)
           >>> p2 = LatLon(48.857, 2.351)
           >>> b = p1.initialBearingTo(p2)  # 156.2

           @JSname: I{bearingTo}.
        '''
        self.others(other)

        a1, b1 = self.to2ab()
        a2, b2 = other.to2ab()

        db, b2 = unrollPI(b1, b2, wrap=wrap)
        ca1, ca2, cdb = map1(cos, a1, a2, db)
        sa1, sa2, sdb = map1(sin, a1, a2, db)

        # see <http://mathforum.org/library/drmath/view/55417.html>
        x = ca1 * sa2 - sa1 * ca2 * cdb
        y = sdb * ca2

        return degrees360(atan2(y, x))

    bearingTo = initialBearingTo  # for backward compatibility

    def intermediateTo(self, other, fraction, height=None, wrap=False):
        '''Locate the point at given fraction between this and an
           other point.

           @param other: The other point (L{LatLon}).
           @param fraction: Fraction between both points (float, 0.0 =
                            this point, 1.0 = the other point).
           @keyword height: Optional height, overriding the fractional
                            height (meter).
           @keyword wrap: Wrap and unroll longitudes (bool).

           @return: Intermediate point (L{LatLon}).

           @raise TypeError: The I{other} point is not L{LatLon}.

           @example:

           >>> p1 = LatLon(52.205, 0.119)
           >>> p2 = LatLon(48.857, 2.351)
           >>> p = p1.intermediateTo(p2, 0.25)  # 51.3721°N, 000.7073°E

           @JSname: I{intermediatePointTo}.
        '''
        self.others(other)

        a1, b1 = self.to2ab()
        a2, b2 = other.to2ab()

        db, b2 = unrollPI(b1, b2, wrap=wrap)
        r = haversine_(a2, a1, db)
        if r > EPS:
            cb1, cb2, ca1, ca2     = map1(cos, b1, b2, a1, a2)
            sb1, sb2, sa1, sa2, sr = map1(sin, b1, b2, a1, a2, r)

            A = sin((1 - fraction) * r) / sr
            B = sin(     fraction  * r) / sr

            x = A * ca1 * cb1 + B * ca2 * cb2
            y = A * ca1 * sb1 + B * ca2 * sb2
            z = A * sa1       + B * sa2

            a = atan2(z, hypot(x, y))
            b = atan2(y, x)

        else:  # points too close
            a = favg(a1, a2, f=fraction)
            b = favg(b1, b2, f=fraction)

        if height is None:
            h = self._havg(other, f=fraction)
        else:
            h = height
        return self.classof(degrees90(a), degrees180(b), height=h)

    def intersection(self, bearing, start2, bearing2,
                           height=None, wrap=False):
        '''Locate the intersection of two paths each defined by
           a start point and an initial bearing.

           @param bearing: Initial bearing from this point (compass degrees).
           @param start2: Start point of second path (L{LatLon}).
           @param bearing2: Initial bearing from start2 (compass degrees).
           @keyword height: Optional height for intersection point,
                            overriding the mean height (meter).
           @keyword wrap: Wrap and unroll longitudes (bool).

           @return: Intersection point (L{LatLon}).

           @raise TypeError: Point I{start2} is not a L{LatLon}.

           @raise ValueError: Intersection is ambiguous or infinite
                              or the paths are parallel or coincide.

           @example:

           >>> p = LatLon(51.8853, 0.2545)
           >>> s = LatLon(49.0034, 2.5735)
           >>> i = p.intersection(108.547, s, 32.435)  # '50.9078°N, 004.5084°E'
        '''
        return intersection(self, bearing, start2, bearing2,
                                  height=height, wrap=wrap,
                                  LatLon=self.classof)

    def isEnclosedBy(self, points):
        '''Test whether this point is enclosed by the polygon defined
           as a list, sequence, set or tuple of points.

           @param points: The points defining the polygon (L{LatLon}[]).

           @return: True if the polygon encloses this point (bool).

           @raise ValueError: Insufficient number of I{points} or
                              non-convex polygon.

           @raise TypeError: Some I{points} are not L{LatLon}.

           @example:

           >>> b = LatLon(45,1), LatLon(45,2), LatLon(46,2), LatLon(46,1)
           >>> p = LatLon(45,1, 1.1)
           >>> inside = p.isEnclosedBy(b)  # True
        '''
        n, points = self.points(points, closed=True)

        n0 = self.toVector3d()

        if iterNumpy2(points):

            v1 = points[n-1].toVector3d()
            v2 = points[n-2].toVector3d()
            gc1 = v2.cross(v1)
            t0 = gc1.angleTo(n0) > PI_2
            for i in range(n):
                v2 = points[i].toVector3d()
                gc = v1.cross(v2)
                v1 = v2

                ti = gc.angleTo(n0) > PI_2
                if ti != t0:
                    return False  # outside

                if gc1.angleTo(gc, vSign=n0) < 0:
                    raise ValueError('non-convex: %r...' % (points[:2],))
                gc1 = gc

        else:
            # get great-circle vector for each edge
            gc, v1 = [], points[n-1].toVector3d()
            for i in range(n):
                v2 = points[i].toVector3d()
                gc.append(v1.cross(v2))
                v1 = v2

            # check whether this point on same side of all
            # polygon edges (to the left or right depending
            # on anti-/clockwise polygon direction)
            t0 = gc[0].angleTo(n0) > PI_2  # True if on the right
            for i in range(1, n):
                ti = gc[i].angleTo(n0) > PI_2
                if ti != t0:  # different sides of edge i
                    return False  # outside

            # check for convex polygon (otherwise
            # the test above is not reliable)
            gc1 = gc[n-1]
            for gc2 in gc:
                # angle between gc vectors, signed by direction of n0
                if gc1.angleTo(gc2, vSign=n0) < 0:
                    raise ValueError('non-convex: %r...' % (points[:2],))
                gc1 = gc2

        return True  # inside

    def midpointTo(self, other, height=None, wrap=False):
        '''Find the midpoint between this and an other point.

           @param other: The other point (L{LatLon}).
           @keyword height: Optional height for midpoint, overriding
                            the mean height (meter).
           @keyword wrap: Wrap and unroll longitudes (bool).

           @return: Midpoint (L{LatLon}).

           @raise TypeError: The I{other} point is not L{LatLon}.

           @example:

           >>> p1 = LatLon(52.205, 0.119)
           >>> p2 = LatLon(48.857, 2.351)
           >>> m = p1.midpointTo(p2)  # '50.5363°N, 001.2746°E'
        '''
        self.others(other)

        # see <http://mathforum.org/library/drmath/view/51822.html>
        a1, b1 = self.to2ab()
        a2, b2 = other.to2ab()

        db, b2 = unrollPI(b1, b2, wrap=wrap)

        ca1, ca2, cdb = map1(cos, a1, a2, db)
        sa1, sa2, sdb = map1(sin, a1, a2, db)

        x = ca2 * cdb + ca1
        y = ca2 * sdb

        a = atan2(sa1 + sa2, hypot(x, y))
        b = atan2(y, x) + b1

        if height is None:
            h = self._havg(other)
        else:
            h = height
        return self.classof(degrees90(a), degrees180(b), height=h)

    def nearestOn(self, point1, point2, **options):
        '''Locate the closest point between two points and this point.

           If this point is within the extent of the segment between
           both end points, the returned point is on the segment.
           Otherwise the returned point is the closest of the segment
           end points.

           Distances are approximated by function L{equirectangular_},
           subject to the supplied I{options}.

           @param point1: Start point of the segment (L{LatLon}).
           @param point2: End point of the segment (L{LatLon}).
           @keyword options: Optional keyword argument for function
                             L{equirectangular_}.

           @return: Closest point on segment (L{LatLon}).

           @raise LimitError: Lat- and/or longitudinal delta exceeds
                              I{limit}, see function L{equirectangular_}.

           @raise TypeError: If I{point1} or I{point2} is not L{LatLon}.
        '''
        p, _ = nearestOn2(self, [point1, point2], **options)
        return p

    def nearestOn2(self, points, radius=R_M, **options):
        '''Locate the closest point on any segment between two
           consecutive points of a path.

           If this point is within the extent of any segment, the
           closest point is on the segment.  Otherwise the closest
           point is the nearest of the segment end points.

           Distances are approximated by function L{equirectangular_},
           subject to the supplied I{options}.

           @param points: The points of the path (L{LatLon}[]).
           @keyword radius: Optional, mean earth radius (meter).
           @keyword options: Optional keyword argument for function
                             L{equirectangular_}.

           @return: 2-Tuple (closest, distance) of the closest point
           (L{LatLon}) on the path and the distance to that point in
           meter, rather in the units of I{radius}.

           @raise LimitError: Lat- and/or longitudinal delta exceeds
                              I{limit}, see function L{equirectangular_}.

           @raise TypeError: Some I{points} are not I{LatLon}.

           @raise ValueError: Insufficient number of I{points}.
        '''
        return nearestOn2(self, points, radius=radius, **options)

    def toVector3d(self):
        '''Convert this point to a vector normal to earth's surface.

           @return: Vector representing this point (L{Vector3d}).
        '''
        if self._v3d is None:
            x, y, z = self.to3xyz()
            self._v3d = Vector3d(x, y, z)  # XXX .unit()
        return self._v3d


_Trll = LatLon(0, 0)  #: (INTERNAL) Reference instance (L{LatLon}).


def _destination2(a, b, r, t):
    '''(INTERNAL) Computes destination lat-/longitude.

       @param a: Latitude (radians).
       @param b: Longitude (radians).
       @param r: Angular distance (radians).
       @param t: Bearing (compass radians).

       @return: 2-Tuple (lat, lon) of (degrees90, degrees180).
    '''
    # see <http://www.edwilliams.org/avform.htm#LL>
    ca, cr, ct = map1(cos, a, r, t)
    sa, sr, st = map1(sin, a, r, t)

    a = asin(ct * sr * ca + cr * sa)
    d = atan2(st * sr * ca, cr - sa * sin(a))
    # note, use addition, not subtraction of d (since
    # in EdWilliams.org/avform.htm W is + and E is -)
    return degrees90(a), degrees180(b + d)


def areaOf(points, radius=R_M, wrap=True):
    '''Calculate the area of a spherical polygon where the sides
       of the polygon are great circle arcs joining the points.

       @param points: The points defining the polygon (L{LatLon}[]).
       @keyword radius: Optional, mean earth radius (meter).
       @keyword wrap: Wrap and unroll longitudes (bool).

       @return: Polygon area (float, same units as radius squared).

       @raise TypeError: Some I{points} are not L{LatLon}.

       @raise ValueError: Insufficient number of I{points}.

       @note: The area is based on Karney U{'Area of a spherical polygon'
              <http://osgeo-org.1560.x6.nabble.com/Area-of-a-spherical-polygon-td3841625.html>}.

       @see: L{pygeodesy.areaOf}, L{sphericalNvector.areaOf} and
             L{ellipsoidalVincenty.areaOf}.

       @example:

       >>> b = LatLon(45, 1), LatLon(45, 2), LatLon(46, 2), LatLon(46, 1)
       >>> areaOf(b)  # 8666058750.718977

       >>> c = LatLon(0, 0), LatLon(1, 0), LatLon(0, 1)
       >>> areaOf(c)  # 6.18e9
    '''
    n, points = _Trll.points(points, closed=True)

    # Area method due to Karney: for each edge of the polygon,
    #
    #                tan(Δλ/2) · (tan(φ1/2) + tan(φ2/2))
    #     tan(E/2) = ------------------------------------
    #                     1 + tan(φ1/2) · tan(φ2/2)
    #
    # where E is the spherical excess of the trapezium obtained by
    # extending the edge to the equator-circle vector for each edge

    if iterNumpy2(points):

        def _exs(n, points):  # iterate over spherical edge excess
            a1, b1 = points[n-1].to2ab()
            ta1 = tan_2(a1)
            for i in range(n):
                a2, b2 = points[i].to2ab()
                db, b2 = unrollPI(b1, b2, wrap=wrap)
                ta2, tdb = map1(tan_2, a2, db)
                yield atan2(tdb * (ta1 + ta2), 1 + ta1 * ta2)
                ta1, b1 = ta2, b2

        s = fsum(_exs(n, points)) * 2

    else:
        a1, b1 = points[n-1].to2ab()
        s, ta1 = [], tan_2(a1)
        for i in range(n):
            a2, b2 = points[i].to2ab()
            db, b2 = unrollPI(b1, b2, wrap=wrap)
            ta2, tdb = map1(tan_2, a2, db)
            s.append(atan2(tdb * (ta1 + ta2), 1 + ta1 * ta2))
            ta1, b1 = ta2, b2

        s = fsum(s) * 2

    if isPoleEnclosedBy(points):
        s = abs(s) - PI2

    return abs(s * radius**2)


def intersection(start1, bearing1, start2, bearing2,
                 height=None, wrap=False, LatLon=LatLon):
    '''Compute the intersection point of two paths each defined
       by a start point and an initial bearing.

       @param start1: Start point of first path (L{LatLon}).
       @param bearing1: Initial bearing from start1 (compass degrees).
       @param start2: Start point of second path (L{LatLon}).
       @param bearing2: Initial bearing from start2 (compass degrees).
       @keyword height: Optional height for the intersection point,
                        overriding the mean height (meter).
       @keyword wrap: Wrap and unroll longitudes (bool).
       @keyword LatLon: Optional LatLon class for the intersection
                        point (L{LatLon}) or None.

       @return: Intersection point (L{LatLon}) or 3-tuple
                (lat, lon, height) if I{LatLon} is None.

       @raise TypeError: Point I{start1} or I{start2} is not L{LatLon}.

       @raise ValueError: Intersection is ambiguous or infinite
                          or the paths are parallel or coincide.

       @example:

       >>> p = LatLon(51.8853, 0.2545)
       >>> s = LatLon(49.0034, 2.5735)
       >>> i = intersection(p, 108.547, s, 32.435)  # '50.9078°N, 004.5084°E'
    '''
    _Trll.others(start1, name='start1')
    _Trll.others(start2, name='start2')

    # see <http://www.edwilliams.org/avform.htm#Intersection>
    a1, b1 = start1.to2ab()
    a2, b2 = start2.to2ab()

    db, b2 = unrollPI(b1, b2, wrap=wrap)
    r12 = haversine_(a2, a1, db)
    if abs(r12) < EPS:
        raise ValueError('intersection %s: %r vs %r' % ('parallel', start1, start2))

    ca1, ca2       = map1(cos, a1, a2)
    sa1, sa2, sr12 = map1(sin, a1, a2, r12)
    x1, x2 = (sr12 * ca1), (sr12 * ca2)
    if min(abs(x1), abs(x2)) < EPS:
        raise ValueError('intersection %s: %r vs %r' % ('parallel', start1, start2))

    cr12 = cos(r12)
    t1, t2 = map1(acos, (sa2 - sa1 * cr12) / x1,
                        (sa1 - sa2 * cr12) / x2)
    if sin(db) > 0:
        t12, t21 = t1, PI2 - t2
    else:
        t12, t21 = PI2 - t1, t2

    t13, t23 = map1(radians, bearing1, bearing2)
    x1, x2 = map1(wrapPI, t13 - t12,  # angle 2-1-3
                          t21 - t23)  # angle 1-2-3

    sx1, sx2 = map1(sin, x1, x2)
    if sx1 == 0 and sx2 == 0:
        raise ValueError('intersection %s: %r vs %r' % ('infinite', start1, start2))
    sx3 = sx1 * sx2
    if sx3 < 0:
        raise ValueError('intersection %s: %r vs %r' % ('ambiguous', start1, start2))
    cx1, cx2 = map1(cos, x1, x2)

    x3 = acos(cr12 * sx3 - cx2 * cx1)
    r13 = atan2(sr12 * sx3, cx2 + cx1 * cos(x3))

    a, b = _destination2(a1, b1, r13, t13)
    h = start1._havg(start2) if height is None else height
    return (a, b, h) if LatLon is None else LatLon(a, b, height=h)


def isPoleEnclosedBy(points, wrap=False):
    '''Test whether a pole is enclosed by a polygon defined by a list,
       sequence, set or tuple of points.

       @param points: The points defining the polygon (L{LatLon}[]).
       @keyword wrap: Wrap and unroll longitudes (bool).

       @return: True if the polygon encloses a pole (bool).

       @raise ValueError: Insufficient number of I{points}.

       @raise TypeError: Some I{points} are not L{LatLon}.
    '''
    n, points = _Trll.points(points)

    def _cds(n, points):  # iterate over course deltas
        p1 = points[n-1]
        b1 = p1.initialBearingTo(points[0], wrap=wrap)  # XXX p1.finalBearingTo(points[0])?
        for i in range(n):
            p2 = points[i]
            if not p2.equals(p1, EPS):
                b = p1.initialBearingTo(p2, wrap=wrap)
                yield wrap180(b - b1)  # (b - b1 + 540) % 360 - 180
                b2 = p1.finalBearingTo(p2, wrap=wrap)
                yield wrap180(b2 - b)  # (b2 - b + 540) % 360 - 180
                p1, b1 = p2, b2

    # sum of course deltas around pole is 0° rather than normally ±360°
    # <http://blog.element84.com/determining-if-a-spherical-polygon-contains-a-pole.html>
    s = fsum(_cds(n, points))

    # XXX fix (intermittant) edge crossing pole - eg (85,90), (85,0), (85,-90)
    return abs(s) < 90  # "zero-ish"


def meanOf(points, height=None, LatLon=LatLon):
    '''Compute the geographic mean of the supplied points.

       @param points: Points to be averaged (L{LatLon}[]).
       @keyword height: Optional height at mean point overriding
                        the mean height (meter).
       @keyword LatLon: Optional LatLon class to return mean point
                        (L{LatLon}) or None.

       @return: Point at geographic mean and height (L{LatLon}) or
                3-tuple (mean_lat, mean_lon, height) if I{LatLon}
                is None.

       @raise TypeError: Some I{points} are not L{LatLon}.

       @raise ValueError: No I{points}.
    '''
    # geographic mean
    n, points = _Trll.points(points, closed=False)

    m = sumOf(points[i].toVector3d() for i in range(n))
    a, b = m.to2ll()

    if height is None:
        h = fmean(points[i].height for i in range(n))
    else:
        h = height
    return (a, b, h) if LatLon is None else LatLon(a, b, height=h)


def nearestOn2(point, points, radius=R_M, **options):
    '''Locate the closest point on any segment between two consecutive
       points of a path.

       If the given point is within the extent of any segment, the
       closest point is on the segment.  Otherwise the closest point
       is the nearest of the segment end points.

       Distances are approximated by function L{equirectangular_},
       subject to the supplied I{options}.

       @param point: The reference point (L{LatLon}).
       @param points: The points of the path (L{LatLon}[]).
       @keyword radius: Optional, mean earth radius (meter).
       @keyword options: Optional keyword arguments for function
                         L{equirectangular_}.

       @return: 2-Tuple (closest, distance) of the closest point
                (L{LatLon}) on the path and the distance to that
                point.  The distance is the L{equirectangular_}
                distance between the given and the closest point
                in meter, rather the units of I{radius}.

       @raise LimitError: Lat- and/or longitudinal delta exceeds
                          I{limit}, see function L{equirectangular_}.

       @raise TypeError: Some I{points} are not I{LatLon}.

       @raise ValueError: Insufficient number of I{points}.
    '''
    def _d2yx(p2, p1, u):
        # equirectangular_ returns a 4-tuple (distance in
        # degrees squared, delta lat, delta lon, p2.lon
        # unroll/wrap); the previous p2.lon unroll/wrap
        # is also applied to the next edge's p1.lon
        return equirectangular_(p1.lat, p1.lon + u,
                                p2.lat, p2.lon, **options)

    # point (x, y) on axis rotated by angle a ccw:
    #   x' = y * sin(a) + x * cos(a)
    #   y' = y * cos(a) - x * sin(a)
    #
    # distance (w) along and perpendicular (h) to
    # a line thru point (dx, dy) and the origin:
    #   w = (y * dy + x * dx) / hypot(dx, dy)
    #   h = (y * dx - x * dy) / hypot(dx, dy)
    #
    # closest point on that line thru (dx, dy):
    #   xc = w * dx / hypot(dx, dy)
    #   yc = w * dy / hypot(dx, dy)

    n, points = point.points(points, closed=False)

    u = 0
    c = p2 = points[0]
    d, _, _, _ = _d2yx(p2, point, 0)
    for i in range(1, n):
        p1, p2 = p2, points[i]
        # iff wrapped, unroll lon1 (actually previous
        # lon2) like function unroll180/-PI would've
        d21, y21, x21, u = _d2yx(p2, p1, u)
        if d21 > EPS:
            # distance point to p1
            d2, y01, x01, _ = _d2yx(point, p1, 0)
            if d2 > EPS:
                w = y01 * y21 + x01 * x21
                if w > 0:
                    if w < d21:
                        # closest is between p1 and p2, use
                        # original delta's, not y21 and x21
                        f = w / d21
                        p1 = point.classof(favg(p1.lat, p2.lat, f=f),
                                           favg(p1.lon, p2.lon + u, f=f))
                        d2, _, _, _ = _d2yx(point, p1, 0)
                    else:  # p2 is closer
                        d2, _, _, _ = _d2yx(point, p2, u)
                        p1 = p2  # in case d2 < d
            if d2 < d:  # p1 is closer
                c, d = p1, d2

    # distance in degrees squared to meter
    d = radians(sqrt(d)) * float(radius)
    return c, d


def perimeterOf(points, closed=False, radius=R_M, wrap=True):
    '''Compute the perimeter of a polygon/-line defined by an array,
       list, sequence, set or tuple of points.

       @param points: The points defining the polygon (L{LatLon}[]).
       @keyword closed: Optionally, close the polygon/-line (bool).
       @keyword radius: Optional, mean earth radius (meter).
       @keyword wrap: Wrap and unroll longitudes (bool).

       @return: Polygon perimeter (float, same units as radius).

       @raise TypeError: Some I{points} are not L{LatLon}.

       @raise ValueError: Insufficient number of I{points}.

       @note: This perimeter is based on the L{haversine} formula.

       @see: L{pygeodesy.perimeterOf} and L{ellipsoidalVincenty.perimeterOf}.
    '''
    n, points = _Trll.points(points, closed=closed)

    def _rads(n, points, closed):  # angular edge lengths in radians
        if closed:
            m, i = 0, n - 1
        else:
            m, i = 1, 0
        a1, b1 = points[i].to2ab()
        for i in range(m, n):
            a2, b2 = points[i].to2ab()
            db, b2 = unrollPI(b1, b2, wrap=wrap)
            yield haversine_(a2, a1, db)
            a1, b1 = a2, b2

    r = fsum(_rads(n, points, closed))
    return r * float(radius)

# **) MIT License
#
# Copyright (C) 2016-2018 -- mrJean1 at Gmail dot com
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
