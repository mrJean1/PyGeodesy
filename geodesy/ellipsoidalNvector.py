
# -*- coding: utf-8 -*-

# Python implementation of vector-based geodetic (lat-/longitude) functions
# by (C) Chris Veness 2011-2015 published under the same MIT Licence**,
# see <http://www.movable-type.co.uk/scripts/latlong-vectors.html>
#
# These functions work with
# a) geodesic (polar) lat-/longitude points on the earth's surface
# b) 3D vectors used as n-vectors representing points on the earth's
#    surface or vectors normal to the plane of a great circle

# See Kenneth Gade, "A Non-singular Horizontal Position Representation",
# The Journal of Navigation (2010), vol 63, nr 3, pp 395-417.  Also at
# <http://www.navlab.net/Publications/A_Nonsingular_Horizontal_Position_Representation.pdf>

from datum import Datums
from dms import F_D, toDMS
from ellipsoidalBase import _CartesianBase, _LatLonHeightDatumBase
from nvector import NorthPole, _LatLonNvectorBase, \
                    Nvector as _NvectorBase, sumOf
from utils import EPS, EPS1, degrees90, degrees360, \
                  cbrt, fdot, hypot3, radians, fStr
from vector3d import Vector3d
from math import asin, atan2, cos, hypot, sin, sqrt

# all public contants, classes and functions
__all__ = ('Cartesian', 'LatLon', 'Ned', 'Nvector',  # classes
           'meanOf', 'toNed')  # functions
__version__ = '16.11.11'


class Cartesian(_CartesianBase):
    '''Extend with method to convert Cartesian to Nvector and
       Nvector-based LatLon.
    '''
    _Nv = None  # cache Nvector

    def toLatLon(self, datum=Datums.WGS84):  # PYCHOK XXX
        '''Converts this (geocentric) Cartesian (x/y/z) point to
           (ellipsoidal geodetic) LatLon point on the specified datum.

           @param {Datum} [datum=Datums.WGS84] - Datum to use.

           @returns {LatLon} The (ellipsoidal) LatLon point.
        '''
        a, b, h = self.to3latlonheight(datum)
        return LatLon(a, b, height=h, datum=datum)

    def toNvector(self, datum=Datums.WGS84):
        '''Converts this cartesian to an (ellipsoidal) n-vector.

           @param {Datum} [datum=WGS84] - Datum to use for conversion.

           @returns {Nvector} N-vector equivalent to this cartesian.

           @example
           c = Cartesian(3980581, 97, 4966825)
           n = c.toNvector()  # (0.6228, 0.0, 0.7824, 0.0)
        '''
        if self._Nv is None or datum != self._Nv.datum:
            E = datum.ellipsoid
            x, y, z = self.to3tuple()

            # Kenneth Gade eqn 23
            p = (x * x + y * y) * E.a2
            q = (z * z * E.e12) * E.a2
            r = (p + q - E.e4) / 6
            s = (p * q * E.e4) / (4 * r * r * r)
            t = cbrt(1 + s + sqrt(s * (2 + s)))

            u = r * (1 + t + 1 / t)
            v = sqrt(u * u + E.e4 * q)
            w = E.e2 * (u + v - q) / (2 * v)

            k = sqrt(u + v + w * w) - w
            e = k / (k + E.e2)
            d = e * hypot(x, y)

            t = hypot(d, z)
            h = (k + E.e2 - 1) / k * t

            s = e / t
            self._Nv = Nvector(x * s, y * s, z / t, h=h, datum=datum)
        return self._Nv


class LatLon(_LatLonNvectorBase, _LatLonHeightDatumBase):
    '''An n-vector-based ellipsoidal LatLon.

       @constructor
       @param {degrees} lat - Latitude in degrees.
       @param {degrees} lon - Longitude in degrees.
       @param {number} [height=0] - Elevation in meter (must be
                       the same units as the datum's half-axes).
       @param {Datum} [datum=Datums.WGS84] - Datum to use.

       @example
       from ellipsoidalNvector import LatLon
       p = LatLon(52.205, 0.119)  # height=0, datum=Datums.WGS84
    '''
    _Nv = None  # cache Nvector
    _r3 = None  # cache _rotation3

    def _rotation3(self):
        # build rotation matrix from n-vector coordinate frame axes
        if self._r3 is None:
            nv = self.toNvector()  # local (n-vector) coordinate frame

            d = nv.negate()  # down (opposite to n-vector)
            e = NorthPole.cross(nv).unit()  # east (pointing perpendicular to the plane)
            n = e.cross(d)  # north (by right hand rule)

            self._r3 = n, e, d  # matrix rows
        return self._r3

    def _update(self, updated):
        if updated:  # reset caches
            self._Nv = self._r3 = None

#     def bearingTo(self, other):
#         '''Return the initial bearing (forward azimuth) from this
#            to an other point, in compass degrees from North.
#
#            @param {LatLon} other - The other LatLon point.
#            @returns {number} Initial bearing in degrees from North.
#
#            @example
#            p1 = LatLon(52.205, 0.119)
#            p2 = LatLon(48.857, 2.351)
#            b = p1.bearingTo(p2)  # 156.2
#         '''
#         self.others(other)
#
#         v1 = self.toVector3d()
#         v2 = other.toVector3d()
#
#         gc1 = v1.cross(v2)  # gc through v1 & v2
#         gc2 = v1.cross(_NP3)  # gc through v1 & North pole
#
#         # bearing is (signed) angle between gc1 & gc2
#         return degrees360(gc1.angleTo(gc2, v1))

#     def crossTrackDistanceTo(self, start, end, radius=R_M):
#         '''Return (signed) distance from this point to great circle
#            defined by a start point and an end point or bearing.
#
#            @param {LatLon} start - Start point of great circle path.
#            @param {LatLon|number} end - End point of great circle path
#                                         or initial bearing (in degrees)
#                                         from great circle start point.
#            @param {number} [radius=R_M] - Mean radius of earth,
#                                           defaults to meter.
#            @returns {number} Distance to great circle (negative if to
#                              left or positive if to right of path).
#
#            @example
#            p = LatLon(53.2611, -0.7972)
#
#            s = LatLon(53.3206, -1.7297)
#            b = 96.0
#            d = p.crossTrackDistanceTo(s, b)  # -305.7
#
#            e = LatLon(53.1887, 0.1334)
#            d = p.crossTrackDistanceTo(s, e)  # -307.5
#         '''
#         self.others(start, name='start')
#
#         if isinstance(end, LatLon):  # gc by two points
#             gc = start.toVector3d().cross(end.toVector3d())
#         else:  # gc from point and bearing
#             gc = start.greatCircle(end)
#
#         # (signed) angle between point and gc normal vector
#         v = self.toVector3d()
#         a = gc.angleTo(v, v.cross(gc))
#         if a < 0:
#             a = -PI_2 - a
#         else:
#             a =  PI_2 - a
#         return a * float(radius)

    def deltaTo(self, other):
        '''Calculates delta from this point to an other LatLon point.

           The delta is given as a north-east-down NED vector.  Note
           that this is a linear delta, unrelated to a geodesic on
           the ellipsoid.

           Points need not be defined on the same datum.

           @param {LatLon} other - Point delta is to be determined to.

           @returns {Ned} Delta from this point to the other point
                          in local tangent plane of this point.

           @example
           a = LatLon(49.66618, 3.45063)
           b = LatLon(48.88667, 2.37472)
           delta = a.deltaTo(b)  # [N:-86126, E:-78900, D:1069]
           d = delta.length  # 116807.681 m
           b = delta.bearing  # 222.493°
           e = delta.elevation  # -0.5245°
        '''
        self.ellipsoids(other)  # throws TypeError and ValueError

        n, e, d = self._rotation3()
        # get delta in cartesian frame
        dc = other.toCartesian().minus(self.toCartesian())
        # rotate dc to get delta in n-vector reference
        # frame using the rotation matrix row vectors
        return Ned(dc.dot(n), dc.dot(e), dc.dot(d))

#     def destination(self, distance, bearing, radius=R_M):
#         '''Return the destination point after traveling from this point
#            the given distance on the given initial bearing.
#
#            @param {number} distance - Distance traveled, in same units
#                                       as the given earth radius.
#            @param {number} bearing - Initial bearing in degrees from North.
#            @param {number} [radius=R_M] - Mean radius of earth, defaults
#                                           to meter.
#            @returns {LatLon} Destination point.
#
#            @example
#            p = LatLon(51.4778, -0.0015)
#            q = p.destination(7794, 300.7)
#            q.toStr()  # '51.5135°N, 000.0983°W' ?
#         '''
#         r = float(distance) / float(radius)  # angular distance in radians
#         # great circle by starting from this point on given bearing
#         gc = self.greatCircle(bearing)
#
#         v1 = self.toVector3d()
#         x = v1.times(cos(r))  # component of v2 parallel to v1
#         y = gc.cross(v1).times(sin(r))  # component of v2 perpendicular to v1
#
#         v2 = x.plus(y).unit()
#         return v2.toLatLon(height=self.height)

    def destinationNed(self, delta):
        '''Calculates destination point using supplied delta from this point.

           @param {Ned} delta - Delta from this to the other point
                                in local tangent plane of this point.

           @returns {Cartesian} Destination point.

           @example
           a = LatLon(49.66618, 3.45063)
           delta = toNed(116807.681, 222.493, -0.5245)  # [N:-86126, E:-78900, D:1069]
           b = a.destinationNed(delta)  # 48.88667°N, 002.37472°E
        '''
        if not isinstance(delta, Ned):
            raise TypeError('%s not a %s.%s' % ('delta', Ned.__module__, Ned.__name__))

        n, e, d = self._rotation3()
        # convert NED delta to standard Vector3d in coordinate frame of n-vector
        dn = delta.toVector3d().to3tuple()
        # rotate dn to get delta in cartesian (ECEF) coordinate
        # reference frame using the rotation matrix column vectors
        dc = Cartesian(fdot(dn, n.x, e.x, d.x),
                       fdot(dn, n.y, e.y, d.y),
                       fdot(dn, n.z, e.z, d.z))

        # apply (cartesian) delta to this Cartesian to
        # obtain destination point as cartesian
        v = self.toCartesian().plus(dc)  # the plus() gives a plain vector

        return Cartesian(v.x, v.y, v.z).toLatLon(datum=self.datum)

    destinationPoint = destinationNed  # XXX original name

#     def distanceTo(self, other):
#         '''Returns distance from this to an other point.
#
#            @param {LatLon} other - the other LatLon point.
#
#            @returns {meter} Distance between this and the other point.
#
#            @example
#            p = LatLon(52.205, 0.119)
#            q = LatLon(48.857, 2.351);
#            d = p.distanceTo(q)  # 404300
#         '''
#         self.others(other)
#
#         v1 = self.toNvector()
#         v2 = other.toNvector()
#         return v1.angleTo(v2) * self.datum.ellipsoid.R
#
#     distanceTo = distanceTo  # XXX original name
#
#     def distanceTo(self, other, radius=R_M):
#         '''Returns distance from this to an other point.
#
#            @param {LatLon} other - the other LatLon point.
#            @param {number} [radius=R_M] - Mean radius of earth,
#                                           defaults to meter.
#
#            @returns {number} Distance between this and the other
#                              point, in the same units as radius.
#
#            @example
#            p = LatLon(52.205, 0.119)
#            q = LatLon(48.857, 2.351);
#            d = p.distanceTo(q)  # 404300
#         '''
#         self.others(other)
#
#         v1 = self.toVector3d()
#         v2 = other.toVector3d()
#         return v1.angleTo(v2) * float(radius)

    def equals(self, other, eps=None):
        '''Check if this point is equal to an other point.

           @param {LatLon} other - The other point.

           @returns {bool} True if points are identical.

           @example
           p = LatLon(52.205, 0.119)
           q = LatLon(52.205, 0.119)
           e = p.equals(q)  # True
        '''
        return _LatLonHeightDatumBase.equals(self, other, eps=eps) and \
               self.height == other.height and self.datum == other.datum

#     def greatCircle(self, bearing):
#         '''Great circle heading on given bearing from this point.
#
#            Direction of vector is such that initial bearing vector b = c x p.
#
#            @param {number} bearing - Compass bearing in degrees.
#            @returns {Vector3d} Normalised vector representing great circle.
#
#            @example
#            p = LatLon(53.3206, -1.7297)
#            g = p.greatCircle(96.0)
#            g.toStr()  # '(-0.794, 0.129, 0.594)'
#         '''
#         b, a = self.toradians()
#         c = radians(bearing)
#
#         ca, sa = cos(a), sin(a)
#         cb, sb = cos(b), sin(b)
#         cc, sc = cos(c), sin(c)
#
#         return Vector3d(sa * cc - ca * sb * sc,
#                        -ca * cc - sa * sb * sc,
#                         cb * sc)

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
           p = p.intermediateTo(q, 0.25)  # 51.3721°N, 000.7073°E
        '''
        self.others(other)

        if fraction > EPS1:
            i = other
        elif fraction < EPS:  # EPS2
            i = self
        else:
            i = other.toNvector().times(fraction).plus(
                 self.toNvector().times(1 - fraction))
#           i = other.toNvector() * fraction + \
#                self.toNvector() * (1 - fraction)
            i = Nvector(i.x, i.y, i.z).toLatLon()
        return i

    intermediatePointTo = intermediateTo  # XXX original name

#     def intersection(self, end, start2, end2):
#         '''Return the point of intersection of two paths, each defined
#            by a start and end point or a start point and bearing.
#
#            @param {LatLon|degrees} end - End point of first path or
#                                          (initial!) bearing from this
#                                          point.
#            @param {LatLon} start2 - Start point of second path.
#            @param {LatLon|number} end2 - End point of second path or
#                                          (initial!) bearing from start2.
#
#            @returns {LatLon} Intersection point (or None if no unique
#                              intersection is defined).
#
#            @example
#            p1 = LatLon(51.8853, 0.2545); b1 = 108.55
#            p2 = LatLon(49.0034, 2.5735); b2 =  32.44
#            i = p1.intersection(b1, p2, b2)
#            i.toStr()  # 50.9076°N, 004.5086°E
#         '''
#         self.others(start2, name='start2')
#
#         # If gc1 & gc2 are great circles through start and end points
#         # (or defined by start point plus bearing), then candidate
#         # intersections are simply gc1 x gc2 and gc2 x gc1.  Most of
#         # the work is deciding the correct intersection point to select!
#
#         # If bearing is given, that determines which intersection, if both
#         # paths are defined by start/end points, take closer intersection
#
#         # gc1 and gc2 are vectors defining great circles through
#         # start and end points; v x gc gives initial bearing vector
#
#         v1 = self.toVector3d()
#         e1 = isinstance(end, LatLon)
#         if e1:  # path1 defined by end point
#             e1 = end.toVector3d()
#             gc1 = v1.cross(e1)
#         else:  # path1 defined by initial bearing
#             gc1 = self.greatCircle(end)
#
#         v2 = start2.toVector3d()
#         e2 = isinstance(end2, LatLon)
#         if e2:  # path2 defined by end point
#             e2 = end2.toVector3d()
#             gc2 = v2.cross(e2)
#         else:  # path2 defined by initial bearing
#             gc2 = start2.greatCircle(end2)
#
#         # there are two antipodal candidate intersection
#         # points ... we have to choose the one to return
#         i1 = gc1.cross(gc2)
#         i2 = gc2.cross(gc1)
#
#         # selection of intersection point depends on how
#         # paths are defined (by bearings or endpoints)
#         if e1 and e2:  # endpoint+endpoint
#             d = sumOf((v1, v2, e1, e2)).dot(i1)
#         elif e1 and not e2:  # endpoint+bearing
#             # gc2 x v2 . i1 +ve means v2 bearing points to i1
#             d = gc2.cross(v2).dot(i1)
#         elif e2 and not e1:  # bearing+endpoint
#             # gc1 x v1 . i1 +ve means v1 bearing points to i1
#             d = gc1.cross(v1).dot(i1)
#         else:  # bearing+bearing
#             # if gc x v . i1 is +ve, initial bearing is
#             # towards i1, otherwise towards antipodal i2
#             d1 = gc1.cross(v1).dot(i1)
#             d2 = gc2.cross(v2).dot(i1)
#             if d1 > 0 and d2 > 0:
#                 d = 1  # both point to i1
#             elif d1 < 0 and d2 < 0:
#                 d = -1  # both point to i2
#             else:  # d1, d2 opposite signs
#                 # intersection is at further-away intersection
#                 # point, take opposite intersection from mid-
#                 # point of v1 and v2 [is this always true?]
#                 d = -v1.plus(v2).dot(i1)
#         i = i1 if d > 0 else i2
#         return i.toLatLon()

#     def isEnclosedBy(self, points):
#         '''Test whether this point is enclosed by the (convex)
#            polygon defined by a set of points.
#
#            @param {LatLon[]} points - Ordered set of points
#                                       defining the polygon.
#
#            @returns {bool} Whether this point is enclosed.
#
#            @throws {ValueError} If polygon is not convex.
#
#            @example
#            r = LatLon(45,1), LatLon(45,2), LatLon(46,2), LatLon(46,1)
#            p = LatLon(45.1, 1.1)
#            i = p.isEnclosedBy(r)  # True
#         '''
#         n, points = len2(points)
#         if n > 0 and points[0].equals(points[n-1]):
#             n -= 1
#         if n < 3:
#             raise ValueError('too few polygon points')
#
#         # get great-circle vector for each edge
#         gc, v2 = [], points[n-1].toVector3d()
#         for i in range(n):
#             v1 = v2
#             v2 = points[i].toVector3d()
#             gc.append(v1.cross(v2))
#
#         v = self.toVector3d()
#         # check whether this point on same side of all
#         # polygon edges (to the left or right depending
#         # on anti-/clockwise polygon direction)
#         t0 = gc[0].angleTo(v) > PI_2  # True if on the right
#         for i in range(1, n):
#             ti = gc[i].angleTo(v) > PI_2
#             if ti != t0:  # different sides of edge i
#                 return False  # outside
#
#         # check for convex polygon (otherwise
#         # the test above is not reliable)
#         gc2 = gc[n-1]
#         for i in range(n):
#             gc1 = gc2
#             gc2 = gc[i]
#             # angle between gc vectors,
#             # signed by direction of v
#             if gc1.angleTo(gc2, v) < 0:
#                 raise ValueError('polygon is not convex')
#
#         return True  # inside

#     def midpointTo(self, other):
#         '''Return the midpoint between this and an other point.
#
#            @param {LatLon} other - the other LatLon point.
#            @returns {LatLon} Midpoint between both points.
#
#            @example
#            p = LatLon(52.205, 0.119)
#            q = LatLon(48.857, 2.351)
#            m = p.midpointTo(q)
#            m.toStr()  # '50.5363°N, 001.2746°E'
#         '''
#         self.others(other)
#
#         v1 = self.toVector3d()
#         v2 = other.toVector3d()
#         m = v1.plus(v2).unit()
#         return m.toLatLon(self._alter(other))

    def toCartesian(self):
        '''Convert this (geodetic) LatLon point to (geocentric) x/y/z
           cartesian coordinates.

           @returns {Cartesian} Cartesian point equivalent, with x,
                                y and z in meter from earth center.
        '''
        x, y, z = self.to3xyz()  # ellipsoidalBase._LatLonHeightDatumBase
        return Cartesian(x, y, z)  # this ellipsoidalNvector Cartesian

    def toNvector(self):  # note: replicated in LatLonNvectorSpherical
        '''Convert this (geodetic) LatLon point to n-vector (normal
           to the earth's surface).

           @returns {Nvector} N-vector representing this LatLon.

           @example
           p = LatLon(45, 45)
           n = p.toNvector()
           n.toStr()  # [0.50000, 0.50000, 0.70710]
        '''
        if self._Nv is None:
            x, y, z, h = self.to4xyzh()  # nvector._LatLonNvectorBase
            self._Nv = Nvector(x, y, z, h=h, datum=self.datum)
        return self._Nv

#     def toVector3d(self):
#         '''Converts this point to a Vector3d (normal to earth's surface).
#
#            @returns {Vector3d} Normalised n-vector representing LatLon point.
#
#            @example
#            p = LatLon(45, 45)
#            v = p.toVector3d()
#            v.toStr()  # '(0.500. 0.500. 0.707)'
#         '''
#         if self._v3d is None:
#             a, b = self.toradians()
#             ca = cos(a)
#             self._v3d = Vector3d(ca * cos(b), ca * sin(b), sin(a))
#         return self._v3d


class Ned(object):
    '''North-Eeast-Down (NED), also known as Local Tangent Plane (LTP),
       is a vector in the local coordinate frame of a body.
    '''
    _bearing   = None  # cache bearing
    _elevation = None  # cache elevation
    _length    = None  # cache length

    def __init__(self, north, east, down):
        '''Create North-East-Down vector.

           @param {number} north - North component in meter.
           @param {number} east - East component in meter.
           @param {number} down - Down component (normal to the surface
                                  of the ellipsoid) in meter.

           @example
           from ellipsiodalNvector import Ned
           delta = Ned(110569, 111297, 1936)
           delta.toStr(prec=0)  #  [N:110569, E:111297, D:1936]
        '''
        self.north = north
        self.east  = east
        self.down  = down

    def __str__(self):
        return self.toStr()

    @property
    def bearing(self):
        '''Bearing of this NED vector.

           @returns {degrees360} NED vector bearing in degrees from North.
        '''
        if self._bearing is None:
            self._bearing = degrees360(atan2(self.east, self.north))
        return self._bearing

    @property
    def elevation(self):
        '''Elevation, tilt of this NED vector.

           @returns {degrees90} NED vector elevation in degrees from
                                horizontal (i.e. tangent to ellipsoid
                                surface).
        '''
        if self._elevation is None:
            self._elevation = -degrees90(asin(self.down / self.length))
        return self._elevation

    @property
    def length(self):
        '''Length of this NED vector.

           @returns {number} Length of NED vector in meter.
        '''
        if self._length is None:
            self._length = hypot3(self.north, self.east, self.down)
        return self._length

    def toStr(self, prec=3, fmt='[%s]', sep=', '):  # PYCHOK expected
        '''Return a string representation of this NED vector.

           @param {number} [prec=3] - Number of decimals, unstripped.
           @param {string} [fmt='[%s]'] - Enclosing backets format.
           @param {string} [sep=', '] - Separator between NEDs.

           @returns {string} This Ned as "[N:f, E:f, D:f]" string.
        '''
        t3 = fStr(self.to3tuple(), prec=prec, sep=' ').split()
        return fmt % (sep.join('%s:%s' % t for t in zip('NED', t3)),)

    def toStr2(self, prec=None, fmt='[%s]', sep=', '):  # PYCHOK expected
        '''Return a string representation of this NED vector as
           length, bearing and elevation.

           @param {number} [prec=3] - Number of decimals, unstripped.
           @param {string} [fmt='[%s]'] - Enclosing backets format.
           @param {string} [sep=', '] - Separator between NEDs.

           @returns {string} This Ned as "[L:f, B:degrees360, E:degrees90]" string.
        '''
        t3 = (fStr(self.length, prec=3 if prec is None else prec),
              toDMS(self.bearing, form=F_D, prec=prec, ddd=0),
              toDMS(self.elevation, form=F_D, prec=prec, ddd=0))
        return fmt % (sep.join('%s:%s' % t for t in zip('LBE', t3)),)

    def to3tuple(self):
        '''Return this NED vector as a 3-tuple.

           @returns {(degrees, degrees, number)} 3-Tuple (north, east, down).
        '''
        return self.north, self.east, self.down

    def toVector3d(self):
        '''Return this NED vector as a Vector3d.

           @returns {Vector3d} From (north, east, down).
        '''
        return Vector3d(*self.to3tuple())


class Nvector(_NvectorBase):
    '''An n-vector is a position representation using a (unit) vector
       normal to the Earth ellipsoid.  Unlike lat-/longitude points,
       n-vectors have no singularities or discontinuities.

       For many applications, n-vectors are more convenient to work
       with than other position representations like lat-/longitude,
       earth-centred earth-fixed (ECEF) vectors, UTM coordinates, etc.

       Note commonality with sphericalNvector.Nvector.
    '''
    datum = Datums.WGS84

    def __init__(self, x, y, z, h=0, datum=None):
        '''Create a 3d n-vector normal to the earth's surface.

           @param {number} x - X component.
           @param {number} y - Y component.
           @param {number} z - Z component.
           @param {number} [h=0] - Height above ellipsoid surface in metres.
           @param {Datum} [datum=WGS84] - Datum this n-vector is defined within.

           @example
           from ellipsoidalNvector import Nvector
           v = Nvector(0.5, 0.5, 0.7071, 1)
           v.toLatLon()  # 45.0000°N, 045.0000°E, +1.00m
        '''
        _NvectorBase.__init__(self, x, y, z, h)
        if datum and datum is not Nvector.datum:
            self.datum = datum

    def copy(self):
        '''Return a copy of this vector.

           @returns {Nvector} Copy of this vector.
        '''
        n = _NvectorBase.copy(self)
        if self.datum != n.datum:
            n.datum = self.datum
        return n

    def toLatLon(self):
        '''Converts this n-vector to an (ellipsoidalNvector) LatLon point.

           @returns {LatLon} Point equivalent to this n-vector.

           @example
           v = Nvector(0.5, 0.5, 0.7071)
           p = v.toLatLon()  # 45.0°N, 45.0°E
        '''
        a, b, h = self.to3latlonheight()
        return LatLon(a, b, height=h, datum=self.datum)

    def toCartesian(self):
        '''Convert this n-vector to an (ellipsoidalNvector) Cartesian.

           @returns {Cartesian} A cartesian vector quivalent to this n-vector.

           @example
           v = Nvector(0.5, 0.5, 0.7071)
           c = v.toCartesian()  # [3194434, 3194434, 4487327]
           p = c.toLatLon()  # 45.0°N, 45.0°E
        '''
        E = self.datum.ellipsoid

        x, y, z, h = self.to4tuple()
        # Kenneth Gade eqn (22)
        n = E.b / sqrt(z * z + (x * x + y * y) * E.a2b2)
        r = E.a2b2 * n + h

        return Cartesian(x * r, y * r, z * (n + h))

    def unit(self):
        '''Normalize this vector to unit length.

           @returns {Nvector} Normalised vector.
        '''
        if self._united is None:
            u = _NvectorBase.unit(self)
            if u.datum != self.datum:
                u.datum = self.datum
            self._united = u._united = u
        return self._united


def meanOf(points, datum=Datums.WGS84):
    '''Return the geographic mean of the supplied points.

       @param {LatLon[]} points - Array of LatLon pionts to be averaged.
       @param {Datum} [datum=Datums.WGS84] - Datum to use.

       @returns {LatLon} Point at the geographic mean and mean height.
    '''
    # geographic mean
    m = sumOf(p.toNvector() for p in points)
    lat, lon, _ = m.to3latlonheight()
    return LatLon(lat, lon, height=m.h, datum=datum)


def toNed(distance, bearing, elevation):
    '''Create an NED vector from distance, bearing and elevation
       (in local coordinate system).

       @param {number} distance - NED vector length in meter.
       @param {degrees} bearing - NED vector bearing (in degrees
                                  from North).
       @param {degrees} elevation - NED vector elevation (in
                                    degrees from local coordinate
                                    frame horizontal).

       @returns {Ned} NED vector equivalent to distance, bearing
                      and elevation.
    '''
    b, e = radians(bearing), radians(elevation)

    d = float(distance)
    dce = d * cos(e)

    return Ned(cos(b) * dce,
               sin(b) * dce,
              -sin(e) * d)

fromDistanceBearingElevation = toNed  # XXX original name

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
