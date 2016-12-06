
# -*- coding: utf-8 -*-

# Python implementation of vector-based spherical geodetic (lat-/
# longitude) functions.  Transcribed from JavaScript originals by
# (C) Chris Veness 2011-2016, published under the same MIT Licence**.
# See <http://www.movable-type.co.uk/scripts/latlong-vectors.html> and
# <http://www.movable-type.co.uk/scripts/geodesy/docs/module-latlon-
# nvector-spherical.html>.

# Tools for working with points and paths on (a spherical model of)
# th Earth’s surface using using n-vectors rather than the more common
# spherical trigonometry.  N-vectors make many calculations much simpler,
# and easier to follow, compared with the trigonometric equivalents.

# Based on Kenneth Gade’s ‘Non-singular Horizontal Position Representation’,
# Journal of Navigation (2010), vol 63, 395-417.

# Note that the formulations below take x => 0°N,0°E, y => 0°N,90°E and
# z => 90°N while Gade uses x => 90°N, y => 0°N,90°E, z => 0°N,0°E.

# Also note that on a spherical model earth, an n-vector is equivalent
# to a normalised version of an (ECEF) cartesian coordinate.

from datum import R_M
from nvector import NorthPole, _LatLonNvectorBase, \
                    Nvector as _NvectorBase, sumOf
from sphericalBase import _LatLonSphericalBase
from utils import EPS, EPS1, PI, PI2, PI_2, degrees360, \
                  fsum, isscalar, len2
from math import atan2, cos, radians, sin

# all public contants, classes and functions
__all__ = ('LatLon',  # classes
           'areaOf', 'intersection', 'meanOf',  # functions
           'triangulate', 'trilaterate')
__version__ = '16.11.20'


class LatLon(_LatLonNvectorBase, _LatLonSphericalBase):
    '''Create a LatLon point on spherical Earth model.

       @classdesc Tools for working with points and paths
       on (a spherical model of) the earth's surface using
       spherical trigonometry.

       @constructor
       @param {number} lat - Latitude in degrees.
       @param {number} lon - Longitude in degrees.

       @example
       from sphericalNvector import LatLon
       p = LatLon(52.205, 0.119)
    '''
    _Nv = None  # cache Nvector

    def _update(self, updated):
        if updated:  # reset caches
            self._Nv = None
#           _LatLonNvectorBase._update(self, updated)

    def _gc3(self, start, end, name):
        # private, return great circle, start and end
        s = start.toNvector()
        if isscalar(end):  # bearing
            gc = s.greatCircle(end)
            e = None
        else:
            self.others(end, name=name)
            e = end.toNvector()
            gc = s.cross(e)  # XXX .unit()?
        return gc, s, e

    def bearingTo(self, other):
        '''Return the initial bearing (forward azimuth) from this
           to an other point, in compass degrees from North.

           @param {LatLon} other - The other LatLon point.

           @returns {degrees360} Initial bearing in degrees from North.

           @throws {TypeError} If other is not a LatLon.

           @example
           p1 = LatLon(52.205, 0.119)
           p2 = LatLon(48.857, 2.351)
           b = p1.bearingTo(p2)  # 156.2
        '''
        gc1 = self.greatCircleTo(other)
        gc2 = self.toNvector().cross(NorthPole)
#       gc2 = self.greatCircleTo(NorthPole)

        return degrees360(gc1.angleTo(gc2, vSign=self.toNvector()))

    def crossTrackDistanceTo(self, start, end, radius=R_M):
        '''Return (signed) distance from this point to great circle
           defined by a start and end point.

           @param {LatLon} start - Start point of great circle path.
           @param {LatLon} end - End point of great circle path.
           @param {number} [radius=R_M] - Mean radius of earth,
                                          defaults to meter.

           @returns {number} Distance to great circle (negative if to
                             left or positive if to right of path).

           @example
           p = LatLon(53.2611, -0.7972)

           s = LatLon(53.3206, -1.7297)
           d = p.crossTrackDistanceTo(s, 96)  # -305.7

           e = LatLon(53.1887, 0.1334)
           d = p.crossTrackDistanceTo(s, e)  # -307.5
        '''
        self.others(start, name='start')
        gc, _, _ = self._gc3(start, end, 'end')

        p = self.toNvector()
        return (gc.angleTo(p) - PI_2) * radius

    def destination(self, distance, bearing, radius=R_M):
        '''Return the destination LatLon from this point after having
           travelled the given distance on the given initial bearing.

           @param {number} distance - Distance travelled (in same units
                                      as the given earth radius).
           @param {degrees} bearing - Initial bearing in degrees from North.
           @param {number} [radius=R_M] - Mean radius of earth (default
                                          the WGS84 mean in meter).

           @returns {LatLon} Destination point.

           @example
           p = LatLon(51.4778, -0.0015)
           q = p.destination(7794, 300.7)
           q.toStr()  # 51.513546°N, 000.098345°W
        '''
        p = self.toNvector()
        e = NorthPole.cross(p).unit()  # east vector at p
        n = p.cross(e)  # north vector at p

        t = radians(bearing)
        q = n.times(cos(t)).plus(e.times(sin(t)))  # direction vector @ p

        r = float(distance) / float(radius)  # angular distance in radians
        n = p.times(cos(r)).plus(q.times(sin(r)))
        return Nvector(n.x, n.y, n.z).toLatLon()

    destinationPoint = destination  # XXX original name

    def distanceTo(self, other, radius=R_M):
        '''Returns distance from this to an other point.

           @param {LatLon} other - the other LatLon point.
           @param {number} [radius=R_M] - Mean radius of earth (default
                                          the WGS84 mean in meter).

           @returns {number} Distance between this and the other point
                             (in the same units as radius).

           @example
           p = LatLon(52.205, 0.119)
           q = LatLon(48.857, 2.351);
           d = p.distanceTo(q)  # 404300
        '''
        self.others(other)

        return self.toNvector().angleTo(other.toNvector()) * radius

    def greatCircle(self, bearing):
        '''Vector normal to great circle obtained by heading on the
           given compass bearing from this point.

           Direction of vector is such that initial bearing vector
           b = c x n, where n is an n-vector representing this point.

           @param {degrees360} bearing - Compass bearing.

           @returns {Nvector} Normalised vector representing great circle.

           @example
           p = LatLon(53.3206, -1.7297)
           gc = p.greatCircle(96.0)
           gc.toStr()  # [-0.794, 0.129, 0.594]
        '''
        a, b = self.toradians()
        t = radians(bearing)

        ca, sa = cos(a), sin(a)
        cb, sb = cos(b), sin(b)
        ct, st = cos(t), sin(t)

        return Nvector(sb * ct - sa * cb * st,
                      -cb * ct - sa * sb * st,
                       ca * st)

    def greatCircleTo(self, other):
        '''Vector normal to great circle obtained by heading from
           this point to another point on a given compass bearing
           from this point.

           Direction of vector is such that initial bearing vector
           b = c x n, where n is an n-vector representing this point.

           @param {LatLon|degrees360} other - The other point or the
                                              compass bearing from this
                                              point.

           @returns {Nvector} Normalised vector representing great circle.

           @example
           p = LatLon(53.3206, -1.7297)
           gc = p.greatCircle(96.0)
           gc.toStr()  # (-0.79408, 0.12856, 0.59406)

           q = LatLon(53.1887, 0.1334)
           g = p.greatCircleTo(q)
           g.toStr()  # (-0.79408, 0.12859, 0.59406)
        '''
        gc, _, _ = self._gc3(self, other, 'other')
        return gc.unit()

    def intermediateChordTo(self, other, fraction):
        '''Returns the point projected from the point at given fraction
           on a straight line (chord) between this and an other point.

           @param {LatLon} other - The other point.
           @param {number} fraction - Fraction between both points
                                      0 = this point, 1 = other point.

           @returns {LatLon} Intermediate point.

           @example
           p = LatLon(52.205, 0.119)
           q = LatLon(48.857, 2.351)
           i = p.intermediateChordTo(q, 0.25)  # 51.3723°N, 000.7072°E
        '''
        self.others(other)

        if fraction > EPS1:
            i = other
        elif fraction < EPS:  # EPS2
            i = self
        else:
            i = self.toNvector().times(1 - fraction).plus(
               other.toNvector().times(fraction))
#           i = other.toNvector() * fraction + \
#                self.toNvector() * (1 - fraction))
            i = Nvector(i.x, i.y, i.z).toLatLon()
        return i

    intermediatePointDirectlyTo = intermediateChordTo  # XXX original name

    def intermediateTo(self, other, fraction):
        '''Returns the point at given fraction between this and
           an other point.

           @param {LatLon} other - The other point.
           @param {number} fraction - Fraction between both points
                                      0 = this point, 1 = other point.

           @returns {LatLon} Intermediate point.

           @example
           p = LatLon(52.205, 0.119)
           q = LatLon(48.857, 2.351)
           i = p.intermediateTo(q, 0.25)  # 51.3721°N, 000.7074°E
        '''
        self.others(other)

        if fraction > EPS1:
            i = other
        elif fraction < EPS:  # EPS2
            i = self
        else:
            p = self.toNvector()
            q = other.toNvector()
            x = p.cross(q)
            d = x.unit().cross(p)  # unit(p × q) × p
            # angular distance tan(a) = |p × q| / p ⋅ q
            a = atan2(x.length(), p.dot(q)) * fraction  # interpolated
            i = p.times(cos(a)).plus(d.times(sin(a)))  # p * cosα + d * sinα
            i = Nvector(i.x, i.y, i.z).toLatLon()
        return i

    intermediatePointTo = intermediateTo  # XXX original name

    def intersection(self, end1, start2, end2):
        '''Return the point of intersection of two paths each defined
           by two points or a start point and bearing.

           @param {LatLon|degrees360} end1 - End point of first path
                              or the initial bearing from this point.
           @param {LatLon} start2 - Start point of the second path.
           @param {LatLon|degrees360} end2 - End point of second path
                   or the initial bearing from the second start point.

           @returns {LatLon} Destination point (null if no unique intersection defined)

           @throws {TypeError} Parameter is not LatLon object.

           @example
           s = LatLon(51.8853, 0.2545)
           e = LatLon(49.0034, 2.5735)
           i = s.intersection(108.55, e, 32.44)  # 50.9076°N, 004.5086°E
        '''
        # If gc1 and gc2 are great circles through start and end points
        # (or defined by start point and bearing), then the candidate
        # intersections are simply gc1 × gc2 and gc2 × gc1.  Most of the
        # work is deciding the correct intersection point to select!  If
        # bearing is given, that determines the intersection, but if both
        # paths are defined by start/end points, take closer intersection.
        gc1, s1, e1 = self._gc3(self, end1, 'end1')
        gc2, s2, e2 = self._gc3(start2, end2, 'end2')

        # there are two (antipodal) candidate intersection
        # points ... we have to choose the one to return
        i1 = gc1.cross(gc2)
        i2 = gc2.cross(gc1)

        # selection of intersection point depends on how
        # paths are defined (by bearings or endpoints)
        if e1 and e2:  # endpoint+endpoint
            d = sumOf((s1, s2, e1, e2)).dot(i1)
        elif e1 and not e2:  # endpoint+bearing
            # gc2 x v2 . i1 +ve means v2 bearing points to i1
            d = gc2.cross(s2).dot(i1)
        elif e2 and not e1:  # bearing+endpoint
            # gc1 x v1 . i1 +ve means v1 bearing points to i1
            d = gc1.cross(s1).dot(i1)
        else:  # bearing+bearing
            # if gc x v . i1 is +ve, initial bearing is
            # towards i1, otherwise towards antipodal i2
            d1 = gc1.cross(s1).dot(i1)  # +ve means p1 bearing points to i1
            d2 = gc2.cross(s2).dot(i1)  # +ve means p2 bearing points to i1
            if d1 > 0 and d2 > 0:
                d = 1  # both point to i1
            elif d1 < 0 and d2 < 0:
                d = -1  # both point to i2
            else:  # d1, d2 opposite signs
                # intersection is at further-away intersection
                # point, take opposite intersection from mid-
                # point of v1 and v2 [is this always true?]
                d = -s1.plus(s2).dot(i1)

        i = i1 if d > 0 else i2
        return Nvector(i.x, i.y, i.z).toLatLon()

    def isEnclosedBy(self, points):
        '''Test whether this point is enclosed by the polygon
           defined by a list, set or tuple of points.

           @param {LatLon[]} points - Ordered set of points defining the
                                      vertices of polygon.

           @returns {bool} Whether this point is enclosed by the polygon.

           @throws {ValueError} If polygon is not convex or has too few
                                points.

           @example
           b = LatLon(45,1), LatLon(45,2), LatLon(46,2), LatLon(46,1)
           p = LatLon(45,1, 1.1);
           inside = p.isEnclosedBy(b)  # True
        '''
        # this method uses angle summation test; on a plane, angles for
        # an enclosed point will sum to 360°, angles for an exterior
        # point will sum to 0°. On a sphere, enclosed point angles will
        # sum to less than 360° (due to spherical excess), exterior point
        # angles will be small but non-zero.
        # XXX are winding number optimisations applicable to spherical surface?

        # close the polygon so that the last point equals the first point
        n, points = len2(points)
        if n > 0 and points[0].equals(points[n-1]):
            n -= 1
        if n < 3:
            raise ValueError('too few polygon points: %s' % (n,))

        p = self.toNvector()
        # get vectors from p to each point
        vs = [p.minus(points[i].toNvector()) for i in range(n)]
        vs.append(vs[0])

        # sum subtended angles of each edge (using vector p to determine sign)
        s = fsum(vs[i].angleTo(vs[i + 1], vSign=p) for i in range(n))
        return abs(s) > PI

    enclosedBy = isEnclosedBy  # XXX original name

    def isWithin(self, point1, point2):
        '''Returns whether this point is within the extent of a
           segment joining two other points.

           If this point is not on the great circle defined by both
           points, return whether it is within the area bound by
           perpendiculars to the great circle at each point (in
           the same hemispere).

           @param {LatLon} point1 - First point defining the segment.
           @param {LatLon} point2 - Second point defining the segment.

           @returns {bool} True if this point is within the segment extent.
        '''
        self.others(point1, name='point1')
        self.others(point2, name='point2')

        n0 = self.toNvector()
        n1 = point1.toNvector()
        n2 = point2.toNvector()

        if n0.dot(n1) < 0 or n0.dot(n2) < 0:
            return False  # different hemisphere

        # get vectors representing d0=p0->p1 and d2=p2->p1
        # and dot product d0⋅d2 tells us if p0 is on the
        # p2 side of p1 or on the other side (similarly
        # for d0=p0->p2 and d1=p1->p2 and dot product
        # d0⋅d1 and p0 on the p1 side of p2 or not)
        return n0.minus(n1).dot(n2.minus(n1)) >= 0 and \
               n0.minus(n2).dot(n1.minus(n2)) >= 0

    isWithinExtent = isWithin  # XXX original name

    def midpointTo(self, other):
        '''Return the midpoint between this and an other point.

           @param {LatLon} other - the other LatLon point.

           @returns {LatLon} Midpoint between both points.

           @example
           p1 = LatLon(52.205, 0.119)
           p2 = LatLon(48.857, 2.351)
           m = p1.midpointTo(p2)  # '50.5363°N, 001.2746°E'
        '''
        self.others(other)

        m = self.toNvector().plus(other.toNvector())
        return Nvector(m.x, m.y, m.z).toLatLon()

    def nearestOn(self, point1, point2):
        '''Return the point closest on great circle segment between
           two points and this point.

           If this point is within the extent of the segment, the
           returned point is on the segment between both endpoints.
           Otherwise the returned is the closest of the segment
           endpoints.

           @param {LatLon} point1 - Start point the segment.
           @param {LatLon} point2 - End point the segment.

           @returns {LatLon} Closet point on segment.

           @example
           s1 = LatLon(51.0, 1.0)
           s2 = LatLon(51.0, 2.0)

           s = LatLon(51.0, 1.9)
           p = s.nearestOn(s1, s2)  # 51.0004°N, 001.9000°E

           d = p.distanceTo(s)  # 42.71 m

           s = LatLon(51.0, 2.1)
           p = s.nearestOn(s1, s2)  # 51.0000°N, 002.0000°E
        '''
        if self.isWithin(point1, point2):
            # closer to segment than to its endpoints,
            # find the closest point on the segment
            gc1 = point1.toNvector().cross(point2.toNvector())
            gc2 = self.toNvector().cross(gc1)
            p = gc1.cross(gc2).toLatLon()

            # beyond segment extent, take closer endpoint
        elif self.distanceTo(point1) < self.distanceTo(point2):
            p = point1
        else:
            p = point2

        return p

    nearestPointOnSegment = nearestOn  # XXX original name

    def toNvector(self):
        '''Convert this (geodetic) LatLon point to a (spherical) n-vector
           (normal to the earth's surface).

           @returns {Nvector} N-vector representing this LatLon.

           @example
           p = LatLon(45, 45)
           n = p.toNvector()
           n.toStr()  # [0.50000, 0.50000, 0.70710]
        '''
        if self._Nv is None:
            x, y, z, h = self.to4xyzh()
            self._Nv = Nvector(x, y, z, h)
        return self._Nv

    def triangulate(self, bearing1, other, bearing2):
        '''Locate a LatLon point given two points and bearings
           from this and the other point.

           @param {degrees360} bearing1 - Bearing from this point
                                          (in degrees from North).
           @param {LatLon} other - Second reference point.
           @param {degrees360} bearing2 - Bearing from second point
                                          (in degrees from North).

           @returns {LatLon} Triangulated point.
        '''
        return triangulate(self, bearing1, other, bearing2)

    def trilaterate(self, distance1, other2, distance2, other3, distance3, radius=R_M):
        '''Locates a LatLon point at given distances from three points,
           this and two other ones.

           @param {number} distance1 - Distance to this point (same units as radius).
           @param {LatLon} point2 - Second reference point.
           @param {number} distance2 - Distance to second point (same units as radius).
           @param {LatLon} point3 - Third reference point.
           @param {number} distance3 - Distance to third  point (same units as radius).
           @param {number} [radius=R_M] - Radius of earth (defaults to mean WGS84 radius in meter).

           @returns {LatLon} Triliterated point.

           See <http://wikipedia.org/wiki/Trilateration>
        '''
        return trilaterate(self, distance1, other2, distance2,
                                            other3, distance3, radius=radius)


class Nvector(_NvectorBase):
    '''An n-vector is a position representation using a (unit) vector
       normal to the earth's surface.  Unlike lat-/longitude points,
       n-vectors have no singularities or discontinuities.

       For many applications, n-vectors are more convenient to work
       with than other position representations like lat-/longitude,
       earth-centred earth-fixed (ECEF) vectors, UTM coordinates, etc.

       On a spherical model earth, an n-vector is equivalent to an
       earth-centred earth-fixed (ECEF) vector.
    '''

    def toLatLon(self, height=None):
        '''Convert this n-vector to a (sphericalNvector) LatLon point.

           @param {meter} [height=None] - Height above earth radius.

           @returns {LatLon} Point equivalent to this n-vector.

           @example
           v = Nvector(0.5, 0.5, 0.7071)
           p = v.toLatLon()  # 45.0°N, 45.0°E
        '''
        a, b, h = self.to3latlonheight()
        return LatLon(a, b, height=h if height is None else height)

    def greatCircle(self, bearing):
        '''Return n-vector normal to great circle obtained by heading
           on given bearing from this point given by its n-vector.

           Direction of vector is such that initial bearing vector b =
           c × p.

           @param {degrees360} bearing - Initial bearing in degrees
                                         from North.

           @returns {Nvector} Normalised vector representing great circle.

           @example
           n = LatLon(53.3206, -1.7297).toNvector()
           gc = n.greatCircle(96.0)  # [-0.794, 0.129, 0.594]
        '''
        t = radians(bearing)

        e = NorthPole.cross(self)  # easting
        n = self.cross(e)  # northing

        e = e.times(cos(t) / e.length())
        n = n.times(sin(t) / n.length())
        return n.minus(e)


def areaOf(points, radius=R_M):
    '''Calculate the area of a spherical polygon where the sides
       of the polygon are great circle arcs joining the vertices.

       @param {LatLon[]} points - Ordered set of points defining
                                  the vertices of polygon.
       @param {number} [radius=R_M] - Earth radius (default, mean
                                      WGS84 radius in meter).

       @returns {number} Polygon area in the same units as radius
                         squared.
    '''
    # uses Girard’s theorem: A = [Σθᵢ − (n−2)·π]·R²
    n, points = len2(points)
    if n > 0 and points[0].equals(points[n-1]):
        n -= 1
    if n < 3:
        raise ValueError('too few polygon points: %s' % (n,))

    # get great-circle vector for each edge
    gc, v2 = [], points[n-1].toNvector()
    for i in range(n):
        v1 = v2
        v2 = points[i].toNvector()
        gc.append(v1.cross(v2))
    gc.append(gc[0])  # XXX needed?

    # sum interior angles
    s = fsum(gc[i].angleTo(gc[i + 1]) for i in range(n))
    return (s - PI2 - n * PI) * radius * radius


def intersection(start1, end1, start2, end2):
    '''Return the point of intersection of two paths each defined
       by two points or a start point and bearing.

       @param {LatLon} start1 - Start point of thefirst path.
       @param {LatLon|number} end1 - End point of first path or the
                                    initial bearing from this point.
       @param {LatLon} start2 - Start point of the second path.
       @param {LatLon|number} end2 - End point of second path or the
                                     initial bearing from the second
                                     start point.

       @returns {LatLon} Destination point.

       @throws {TypeError} Parameter is not LatLon object.

       @example
       p = LatLon(51.8853, 0.2545)
       q = LatLon(49.0034, 2.5735)
       i = intersection(p, 108.55, q, 32.44)  # 50.9076°N, 004.5086°E
    '''
    return start1.intersection(end1, start2, end2)


def meanOf(points, height=None):
    '''Return the geographic mean of the supplied points.

       @param {LatLon[]} points - Array of LatLon pionts to be averaged.

       @returns {LatLon} Point at the geographic mean and mean height.
    '''
    # geographic mean
    m = sumOf(p.toNvector() for p in points)
    lat, lon, _ = m.to3latlonheight()
    return LatLon(lat, lon, height=m.h if height is None else height)


def triangulate(point1, bearing1, point2, bearing2):
    '''Locate a LatLon point given two known points and bearings
       from those points.

       @param {LatLon} point1 - First reference point.
       @param {degrees360} bearing1 - Bearing from first point
                                      (in degrees from North).
       @param {LatLon} point2 - Second reference point.
       @param {degrees360} bearing2 - Bearing from second point
                                      (in degrees from North).

       @returns {LatLon} Triangulated point.
    '''
    point1.others(point2, name='point2')

    def _gc(p, b):
        n, t = p.toNvector(), radians(b)
        de = NorthPole.cross(n).unit()  # east vector @ n
        dn = n.cross(de)  # north vector @ n
        dest = de.times(sin(t))
        dnct = dn.times(cos(t))
        d = dnct.plus(dest)  # direction vector @ n
        return n.cross(d)  # great circle point + bearing

    gc1 = _gc(point1, bearing1)  # great circle p1 + b1
    gc2 = _gc(point2, bearing2)  # great circle p2 + b2

    h = point1._alter(point2)
    i = gc1.cross(gc2)  # n-vector of intersection point
    return Nvector(i.x, i.y, i.z).toLatLon(height=h, datum=point1.datum)


def trilaterate(point1, distance1, point2, distance2, point3, distance3, radius=R_M):
    '''Locate a LatLon point at given distances from three other points.

       @param {LatLon} point1 - First reference point.
       @param {number} distance1 - Distance to first point (same units as radius).
       @param {LatLon} point2 - Second reference point.
       @param {number} distance2 - Distance to second point (same units as radius).
       @param {LatLon} point3 - Third reference point.
       @param {number} distance3 - Distance to third  point (same units as radius).
       @param {number} [radius=R_M] - Radius of earth (defaults to mean WGS84 radius in meter).

       @returns {LatLon} Triliterated point.

       See <http://wikipedia.org/wiki/Trilateration>
    '''
    point1.others(point2, name='point2')
    point1.others(point3, name='point3')

    n1, d1 = point1.toNvector(), float(distance1) / radius
    n2, d2 = point2.toNvector(), float(distance2) / radius
    n3, d3 = point3.toNvector(), float(distance3) / radius

    # the following uses x,y coordinate system with origin at n1, x axis n1->n2
    X = n2.minus(n1).unit()  # unit vector in x direction n1->n2
    i = X.dot(n3.minus(n1))  # signed magnitude of x component of n1->n3
    Y = n3.minus(n1).minus(X.times(i)).unit()  # unit vector in y direction
    d = n2.minus(n1).length()  # distance n1->n2
    j = Y.dot(n3.minus(n1))  # signed magnitude of y component of n1->n3

    d12 = d1 * d1
    x = (d12 - d2 * d2 + d * d) / (2 * d)  # x component of n1 -> intersection
    y = (d12 - d3 * d3 + i * i + j * j) / (2 * j) - x * i / j  # y component of n1 -> intersection

    z = x * x + y * y
    if z >= d12:
        raise ValueError('no %s for %r, %r, %r' % ('trilaterate',
                          point1, point2, point3))
#   Z = X.cross(Y)  # unit vec to r perpendicular to plane
    # note don't use z component; assume points at same height
#   z = sqrt(d12 - z)  # z will be NaN for no intersections
    h = sum(p.height for p in (point1, point2, point3)) / 3

    n = n1.plus(X.times(x)).plus(Y.times(y))
    return Nvector(n.x, n.y, n.z).toLatLon(height=h, datum=point1.datum)

# **) MIT License
#
# Copyright (c) 2016-2017 -- mrJean1@Gmail.com
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
