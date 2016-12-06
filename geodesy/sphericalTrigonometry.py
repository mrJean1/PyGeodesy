
# -*- coding: utf-8 -*-

# Python implementation of geodetic (lat-/longitude) functions using
# spherical trigonometry.  Transcribed from JavaScript originals by
# (C) Chris Veness 2011-2016 published under the same MIT Licence**,
# see <http://www.movable-type.co.uk/scripts/latlong.html>

from datum import R_M
from sphericalBase import _LatLonSphericalBase
from utils import EPS, EPS1, EPS2, PI2, PI_2, \
                  degrees90, degrees180, degrees360, \
                  fsum, isscalar, len2, radians, sin_2, wrapPI
from vector3d import Vector3d, sumOf
from math import acos, asin, atan2, cos, hypot, sin, sqrt

# all public contants, classes and functions
__all__ = ('LatLon',  # classes
           'meanOf')  # functions
__version__ = '16.11.11'


class LatLon(_LatLonSphericalBase):
    '''Create a LatLon point on spherical model earth.

       @classdesc Tools for working with points and paths
       on (a spherical model of) the earth's surface using
       spherical trigonometry.

       @constructor
       @param {number} lat - Latitude in degrees.
       @param {number} lon - Longitude in degrees.

       @example
       p = LatLon(52.205, 0.119)  # height=0
    '''

    _v3d = None  # cache Vector3d

    def _update(self, updated):
        if updated:  # reset caches
            self._v3d = None
#           _LatLonHeightBase._update(self, updated)

    def bearingTo(self, other):
        '''Return the initial bearing (forward azimuth) from this
           to an other point, in compass degrees from North.

           @param {LatLon} other - The other LatLon point.

           @returns {degrees360} Initial bearing in degrees from North.

           @example
           p1 = LatLon(52.205, 0.119)
           p2 = LatLon(48.857, 2.351)
           b = p1.bearingTo(p2)  # 156.2
        '''
        self.others(other)

        a1, b1 = self.toradians()
        a2, b2 = other.toradians()
        db = b2 - b1

        ca2 = cos(a2)
        # see http://mathforum.org/library/drmath/view/55417.html
        x = cos(a1) * sin(a2) - sin(a1) * ca2 * cos(db)
        y = sin(db) * ca2

        return degrees360(atan2(y, x))

    def crossingParallels(self, other, lat):
        '''Return the pair of meridians at which a great circle defined
           by this and an other point crosses the given latitude.

           @param {LatLon} point1 - First point defining great circle.
           @param {LatLon} point2 - Second point defining great circle.
           @param {degrees} lat - Latitude at the crossings.

           @returns {(degrees180, degrees180)} 2-Tuple of the meridians.
           @returns {None} If the great circle doesn't reach the given
                           latitude.
        '''
        self.others(other)

        a = radians(lat)
        ca, sa = cos(a), sin(a)

        a1, b1 = self.toradians()
        a2, b2 = other.toradians()

        ca1, sa1 = cos(a1), sin(a1)
        ca2, sa2 = cos(a2), sin(a2)

        db = b2 - b1
        cdb, sdb = cos(db), sin(db)

        x = sa1 * ca2 * ca * sdb
        y = sa1 * ca2 * ca * cdb - ca1 * sa2 * ca
        z = ca1 * ca2 * sa * sdb

        h = hypot(x, y)
        if abs(z) > h:
            return None  # great circle doesn't reach latitude

        bm = atan2(-y, x) + b1  # longitude at max latitude
        bi = acos(z / h)  # delta longitude to intersections

        return degrees180(bm - bi), degrees180(bm + bi)

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
           e = LatLon(53.1887, 0.1334)
           d = p.crossTrackDistanceTo(s, e)  # -307.5
        '''
        if isscalar(end):
            raise self.notImplemented('crossTrackDistanceTo(end=bearing)')
        self.others(start, name='start')
        self.others(end, name='end')

        r = start.distanceTo(self, radius) / float(radius)
        b = radians(start.bearingTo(self))
        e = radians(start.bearingTo(end))

        return asin(sin(r) * sin(b - e)) * radius

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
           p1 = LatLon(51.4778, -0.0015)
           p2 = p1.destination(7794, 300.7)
           p2.toStr()  # '51.5135°N, 000.0983°W'
        '''
        d = float(distance) / float(radius)  # angular distance in radians
        t = radians(bearing)

        # see http://williams.best.vwh.net/avform.htm#LL
        a1, b1 = self.toradians()
        ca1, sa1 = cos(a1), sin(a1)

        cd, sd = cos(d), sin(d)
        a2 = asin(sa1 * cd + ca1 * sd * cos(t))
        b2 = atan2(sin(t) * sd * ca1, cd - sa1 * sin(a2)) + b1

        return LatLon(degrees90(a2), degrees180(b2), height=self.height)

    destinationPoint = destination  # XXX original name

    def distanceTo(self, other, radius=R_M):
        '''Returns distance from this to an other point.

           @param {LatLon} other - the other LatLon point.
           @param {number} [radius=R_M] - Mean radius of earth (default
                                          the WGS84 mean in meter).

           @returns {number} Distance between this and the other point
                             (in the same units as radius).

           @example
           p1 = LatLon(52.205, 0.119)
           p2 = LatLon(48.857, 2.351);
           d = p1.distanceTo(p2)  # 404300
        '''
        self.others(other)

        a1, b1 = self.toradians()
        a2, b2 = other.toradians()

        # see http://williams.best.vwh.net/avform.htm#Dist
        sa = sin_2(a2 - a1)
        sb = sin_2(b2 - b1)

        a = sa * sa + cos(a1) * cos(a2) * sb * sb
        c = 2 * atan2(sqrt(a), sqrt(1 - a))
#       c = 2 * asin(sqrt(a))

        return c * float(radius)

    def greatCircle(self, bearing):
        '''Great circle heading on given bearing from this point.

           Direction of vector is such that initial bearing vector b = c x p.

           @param {degrees} bearing - Compass bearing in degrees from North.

           @returns {Vector3d} Normalised vector representing great circle.

           @example
           p = LatLon(53.3206, -1.7297)
           g = p.greatCircle(96.0)
           g.toStr()  # (-0.794, 0.129, 0.594)
        '''
        a, b = self.toradians()
        t = radians(bearing)

        ca, sa = cos(a), sin(a)
        cb, sb = cos(b), sin(b)
        ct, st = cos(t), sin(t)

        return Vector3d(sb * ct - cb * sa * st,
                       -cb * ct - sb * sa * st,
                        ca * st)  # XXX .unit()?

    def intermediateTo(self, other, fraction):
        '''Returns the point at given fraction between this and
           an other point.

           @param {LatLon} other - The other point.
           @param {number} fraction - Fraction between both points, 0
                                      for this or 1 for other point.

           @returns {LatLon} Intermediate point.

           @example
           p1 = LatLon(52.205, 0.119)
           p2 = LatLon(48.857, 2.351)
           p = p1.intermediateTo(p2, 0.25)  # 51.3721°N, 000.7073°E
        '''
        self.others(other)

        if fraction < EPS:
            i = self
        elif fraction > EPS1:
            i = other
        else:
            a1, b1 = self.toradians()
            a2, b2 = other.toradians()

            ca1, sa1, cb1, sb1 = cos(a1), sin(a1), cos(b1), sin(b1)
            ca2, sa2, cb2, sb2 = cos(a2), sin(a2), cos(b2), sin(b2)

            # distance between points
            da, db = (a2 - a1), (b2 - b1)
            sa, sb = sin_2(da), sin_2(db)

            a = sa * sa + ca1 * ca2 * sb * sb
            d = atan2(sqrt(a), sqrt(1 - a)) * 2
            sd = sin(d)
            if abs(sd) > EPS:
                A = sin((1 - fraction) * d) / sd
                B = sin(     fraction  * d) / sd

                x = A * ca1 * cb1 + B * ca2 * cb2
                y = A * ca1 * sb1 + B * ca2 * sb2
                z = A * sa1       + B * sa2

                a3 = atan2(z, hypot(x, y))
                b3 = atan2(y, x)
            else:  # points too close
                a3 = a1 + fraction * da
                b3 = b1 + fraction * db
            i = LatLon(degrees90(a3), degrees180(b3),
                       height=self._alter(other, f=fraction))
        return i

    def intersection(self, bearing, start2, bearing2):
        '''Return the intersection point of two paths each defined
           by a start point and initial bearing.

           @param {degrees} bearing - Initial bearing from this point.
           @param {LatLon} start2 - Start point of second path.
           @param {degrees} bearing2 - Initial bearing from start2.

           @returns {LatLon} Intersection point.

           @throws {TypeError} Point start2 is not a LatLon.

           @throws {ValueError} Intersection is ambiguous or infinite
                                or the paths are parallel or coincide.

           @example
           p = LatLon(51.8853, 0.2545)
           s = LatLon(49.0034, 2.5735)
           i = p.intersection(108.547, s, 32.435)  # '50.9078°N, 004.5084°E'
        '''
        self.others(start2, name='start2')

        # see http://williams.best.vwh.net/avform.htm#Intersection
        a1, b1 = self.toradians()
        a2, b2 = start2.toradians()

        sa = sin_2(a2 - a1)
        sb = sin_2(b2 - b1)

        # angular distance
        r12 = asin(sqrt(sa * sa + cos(a1) * cos(a2) * sb * sb)) * 2
        if abs(r12) < EPS2:
            raise ValueError('intersection %s: %r vs %r' % ('parallel', self, start2))
        cr12 = cos(r12)
        sr12 = sin(r12)

        sa1 = sin(a1)
        sa2 = sin(a2)
        t1 = acos((sa2 - sa1 * cr12) / (sr12 * cos(a1)))
        t2 = acos((sa1 - sa2 * cr12) / (sr12 * cos(a2)))
        if sin(b2 - b1) > 0:
            t12, t21 = t1, PI2 - t2
        else:
            t12, t21 = PI2 - t1, t2

        t13 = radians(bearing)
        t23 = radians(bearing2)
        x1 = wrapPI(t13 - t12)  # angle 2-1-3
        x2 = wrapPI(t21 - t23)  # angle 1-2-3
        sx1 = sin(x1)
        sx2 = sin(x2)
        if sx1 == 0 and sx2 == 0:
            raise ValueError('intersection %s: %r vs %r' % ('infinite', self, start2))
        sx3 = sx1 * sx2
        if sx3 < 0:
            raise ValueError('intersection %s: %r vs %r' % ('ambiguous', self, start2))
        cx1 = cos(x1)
        cx2 = cos(x2)
        ca1 = cos(a1)

        x3 = acos(-cx1 * cx2 + sx3 * cr12)
        r13 = atan2(sr12 * sx3, cx2 + cx1 * cos(x3))
        cr13 = cos(r13)
        sr13 = sin(r13)

        a3 = asin(sa1 * cr13 + ca1 * sr13 * cos(t13))
        b3 = atan2(sin(t13) * sr13 * ca1, cr13 - sa1 * sin(a3)) + b1

        return LatLon(degrees90(a3), degrees180(b3), height=self._alter(start2))

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
           p = LatLon(45,1, 1.1)
           inside = p.isEnclosedBy(b)  # True
        '''
        n, points = len2(points)
        if n > 0 and points[0].equals(points[n-1]):
            n -= 1
        if n < 3:
            raise ValueError('too few polygon points: %s' % (n,))

        # get great-circle vector for each edge
        gc, v2 = [], points[n-1].toVector3d()
        for i in range(n):
            v1 = v2
            v2 = points[i].toVector3d()
            gc.append(v1.cross(v2))

        v = self.toVector3d()
        # check whether this point on same side of all
        # polygon edges (to the left or right depending
        # on anti-/clockwise polygon direction)
        t0 = gc[0].angleTo(v) > PI_2  # True if on the right
        for i in range(1, n):
            ti = gc[i].angleTo(v) > PI_2
            if ti != t0:  # different sides of edge i
                return False  # outside

        # check for convex polygon (otherwise
        # the test above is not reliable)
        gc2 = gc[n-1]
        for i in range(n):
            gc1 = gc2
            gc2 = gc[i]
            # angle between gc vectors,
            # signed by direction of v
            if gc1.angleTo(gc2, v) < 0:
                raise ValueError('polygon is not convex: %r' % (points[:3],))

        return True  # inside

    def isWithin(self, point1, point2):
        self.others(point1, name='point1')
        self.others(point2, name='point2')
        raise self.notImplemented('isWithin')

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

        # see http://mathforum.org/library/drmath/view/51822.html
        a1, b1 = self.toradians()
        a2, b2 = other.toradians()

        d = b2 - b1
        c = cos(a2)
        x = c * cos(d) + cos(a1)
        y = c * sin(d)

        a3 = atan2(sin(a1) + sin(a2), hypot(x, y))
        b3 = atan2(y, x) + b1

        return LatLon(degrees90(a3), degrees180(b3), height=self._alter(other))

    def nearestOn(self, point1, point2):
        self.others(point1, name='point1')
        self.others(point2, name='point2')
        raise self.notImplemented('nearestOn')

    nearestPointOnSegment = nearestOn  # XXX original name

    def toVector3d(self):
        '''Converts this point to a Vector3d (normal to earth's surface).

           @returns {Vector3d} Normalised n-vector representing LatLon point.
        '''
        if self._v3d is None:
            x, y, z = self.to3xyz()
            self._v3d = Vector3d(x, y, z)  # .unit()
        return self._v3d


def meanOf(points):
    '''Return the geographic mean of the supplied points.

       @param {LatLon[]} points - Array of LatLon pionts to be averaged.

       @returns {LatLon} Point at the geographic mean and mean height.
    '''
    # geographic mean
    n, points = len2(points)
    m = sumOf(p.Vector3d() for p in points)
    lat, lon = m.to2latlon()
    h = fsum(p.height for p in points) / (n or 1)
    return LatLon(lat, lon, height=h)

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
