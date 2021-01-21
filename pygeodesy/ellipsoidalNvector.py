
# -*- coding: utf-8 -*-

u'''Ellispdoidal, N-vector-based classes geodetic (lat-/longitude) L{LatLon},
geocentric (ECEF) L{Cartesian}, L{Ned} and L{Nvector} and functions L{meanOf}.
L{sumOf} and L{toNed}.

Pure Python implementation of n-vector-based geodetic (lat-/longitude)
methods by I{(C) Chris Veness 2011-2016} published under the same MIT
Licence**, see U{Vector-based geodesy
<https://www.Movable-Type.co.UK/scripts/latlong-vectors.html>}.

These classes and functions work with: (a) geodesic (polar) lat-/longitude
points on the earth's surface and (b) 3-D vectors used as n-vectors
representing points on the earth's surface or vectors normal to the plane
of a great circle.

See also Kenneth Gade U{'A Non-singular Horizontal Position Representation'
<https://www.NavLab.net/Publications/A_Nonsingular_Horizontal_Position_Representation.pdf>},
The Journal of Navigation (2010), vol 63, nr 3, pp 395-417.

@newfield example: Example, Examples
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division

from pygeodesy.basics import neg, _xinstanceof
from pygeodesy.datums import Datums, _ellipsoidal_datum
from pygeodesy.ellipsoidalBase import CartesianEllipsoidalBase, \
                                      LatLonEllipsoidalBase
from pygeodesy.errors import _xkwds
from pygeodesy.fmath import fdot, hypot_
from pygeodesy.interns import NN, _COMMASPACE_
from pygeodesy.interns import _pole_  # PYCHOK used!
from pygeodesy.lazily import _ALL_LAZY, _ALL_OTHER
from pygeodesy.named import _Named, _NamedTuple, _xnamed
from pygeodesy.namedTuples import LatLon3Tuple
from pygeodesy.nvectorBase import NorthPole, LatLonNvectorBase, \
                                  NvectorBase, sumOf as _sumOf
from pygeodesy.props import Property_RO, property_RO
from pygeodesy.streprs import Fmt, fstr, strs, _xzipairs
from pygeodesy.units import Bearing, Degrees, Distance, Float,\
                            Height, Radius, Scalar
from pygeodesy.utily import asin1, atan2b, degrees90, sincos2d

__all__ = _ALL_LAZY.ellipsoidalNvector
__version__ = '21.01.19'

_down_  = 'down'
_east_  = 'east'
_north_ = 'north'


class Cartesian(CartesianEllipsoidalBase):
    '''Extended to convert geocentric, L{Cartesian} points to
       L{Nvector} and n-vector-based, geodetic L{LatLon}.
    '''

    _Ecef = None  # preferred C{EcefVeness} class

    @property_RO
    def Ecef(self):
        '''Get the ECEF I{class} (L{EcefVeness}).
        '''
        if Cartesian._Ecef is None:
            from pygeodesy.ecef import EcefVeness
            Cartesian._Ecef = EcefVeness
        return Cartesian._Ecef

    def toLatLon(self, **LatLon_datum_kwds):  # PYCHOK LatLon=LatLon, datum=None
        '''Convert this cartesian point to an C{Nvector}-based
           geodetic point.

           @kwarg LatLon_datum_kwds: Optional L{LatLon}, B{C{datum}} and
                                     other keyword arguments, ignored if
                                     C{B{LatLon}=None}.  Use
                                     B{C{LatLon=...}} to override this
                                     L{LatLon} class or specify
                                     C{B{LatLon}=None}.

           @return: The geodetic point (L{LatLon}) or if B{C{LatLon}}
                    is C{None}, an L{Ecef9Tuple}C{(x, y, z, lat, lon,
                    height, C, M, datum)} with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{LatLon}}, B{C{datum}} or other
                             B{C{LatLon_datum_kwds}}.
        '''
        kwds = _xkwds(LatLon_datum_kwds, LatLon=LatLon, datum=self.datum)
        return CartesianEllipsoidalBase.toLatLon(self, **kwds)

    def toNvector(self, **Nvector_datum_kwds):  # PYCHOK Datums.WGS84
        '''Convert this cartesian to L{Nvector} components,
           I{including height}.

           @kwarg Nvector_datum_kwds: Optional L{Nvector}, B{C{datum}} and
                                      other keyword arguments, ignored if
                                      B{C{Nvector=None}}.  Use
                                      B{C{Nvector=...}} to override this
                                      L{Nvector} class or specify
                                      B{C{Nvector=None}}.

           @return: The C{n-vector} components (L{Nvector}) or a
                    L{Vector4Tuple}C{(x, y, z, h)} if B{C{Nvector=None}}.

           @raise TypeError: Invalid B{C{Nvector}}, B{C{datum}} or other
                             B{C{Nvector_datum_kwds}}.

           @example:

           >>> from ellipsoidalNvector import LatLon
           >>> c = Cartesian(3980581, 97, 4966825)
           >>> n = c.toNvector()  # (0.62282, 0.000002, 0.78237, +0.24)
        '''
        kwds = _xkwds(Nvector_datum_kwds, Nvector=Nvector, datum=self.datum)
        return CartesianEllipsoidalBase.toNvector(self, **kwds)


class LatLon(LatLonNvectorBase, LatLonEllipsoidalBase):
    '''An n-vector-based, ellipsoidal L{LatLon} point.

       @example:

       >>> from ellipsoidalNvector import LatLon
       >>> p = LatLon(52.205, 0.119)  # height=0, datum=Datums.WGS84
    '''
    _Ecef = None  # preferred C{EcefVeness} class
    _Nv   = None  # cached toNvector (L{Nvector})

    def _update(self, updated, *attrs):  # PYCHOK args
        '''(INTERNAL) Zap cached attributes if updated.
        '''
        if updated:
            LatLonNvectorBase._update(self, updated, _Nv=self._Nv)  # special case
            LatLonEllipsoidalBase._update(self, updated, *attrs)

#     def crossTrackDistanceTo(self, start, end, radius=R_M):
#         '''Return the (signed) distance from this point to the great
#            circle defined by a start point and an end point or bearing.
#
#            @arg start: Start point of great circle path (L{LatLon}).
#            @arg end: End point of great circle path (L{LatLon}) or
#                      initial bearing (compass C{degrees360}) at the
#                      start point.
#            @kwarg radius: Mean earth radius (C{meter}).
#
#            @return: Distance to great circle, negative if to left or
#                     positive if to right of path (C{meter}, same units
#                     as B{C{radius}}).
#
#            @raise TypeError: If B{C{start}} or B{C{end}} point is not L{LatLon}.
#
#            @example:
#
#            >>> p = LatLon(53.2611, -0.7972)
#
#            >>> s = LatLon(53.3206, -1.7297)
#            >>> b = 96.0
#            >>> d = p.crossTrackDistanceTo(s, b)  # -305.7
#
#            >>> e = LatLon(53.1887, 0.1334)
#            >>> d = p.crossTrackDistanceTo(s, e)  # -307.5
#         '''
#         self.others(start=start)
#
#         if isscalar(end):  # gc from point and bearing
#             gc = start.greatCircle(end)
#         else:  # gc by two points
#             gc = start.toNvector().cross(end.toNvector())
#
#         # (signed) angle between point and gc normal vector
#         v = self.toNvector()
#         a = gc.angleTo(v, vSign=v.cross(gc))
#         a = (-PI_2 - a) if a < 0 else (PI_2 - a)
#         return a * float(radius)

    def deltaTo(self, other):
        '''Calculate the NED delta from this to an other point.

           The delta is returned as a North-East-Down (NED) vector.

           Note, this is a linear delta, unrelated to a geodesic
           on the ellipsoid.  The points need not be defined on
           the same datum.

           @arg other: The other point (L{LatLon}).

           @return: Delta of this point (L{Ned}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: If ellipsoids are incompatible.

           @example:

           >>> a = LatLon(49.66618, 3.45063)
           >>> b = LatLon(48.88667, 2.37472)
           >>> delta = a.deltaTo(b)  # [N:-86126, E:-78900, D:1069]
           >>> d = delta.length  # 116807.681 m
           >>> b = delta.bearing  # 222.493°
           >>> e = delta.elevation  # -0.5245°
        '''
        self.ellipsoids(other)  # throws TypeError and ValueError

        # get delta in cartesian frame
        dc = other.toCartesian().minus(self.toCartesian())
        # rotate dc to get delta in n-vector reference
        # frame using the rotation matrix row vectors
        n, e, d = self._rotation3
        return Ned(dc.dot(n), dc.dot(e), dc.dot(d), name=self.name)

#     def destination(self, distance, bearing, radius=R_M, height=None):
#         '''Return the destination point after traveling from this
#            point the given distance on the given initial bearing.
#
#            @arg distance: Distance traveled (C{meter}, same units as
#                           given earth B{C{radius}}).
#            @arg bearing: Initial bearing (compass C{degrees360}).
#            @kwarg radius: Mean earth radius (C{meter}).
#            @kwarg height: Optional height at destination point,
#                           overriding default (C{meter}, same units
#                           as B{C{radius}}).
#
#            @return: Destination point (L{LatLon}).
#
#            @example:
#
#            >>> p = LatLon(51.4778, -0.0015)
#            >>> q = p.destination(7794, 300.7)
#            >>> q.toStr()  # '51.5135°N, 000.0983°W' ?
#         '''
#         r = _angular(distance, radius)  # angular distance in radians
#         # great circle by starting from this point on given bearing
#         gc = self.greatCircle(bearing)
#
#         v1 = self.toNvector()
#         x = v1.times(cos(r))  # component of v2 parallel to v1
#         y = gc.cross(v1).times(sin(r))  # component of v2 perpendicular to v1
#
#         v2 = x.plus(y).unit()
#         return v2.toLatLon(height=self.height if height is C{None} else height)

    def destinationNed(self, delta):
        '''Calculate the destination point using the supplied NED delta
           from this point.

           @arg delta: Delta from this to the other point in the local
                       tangent plane (LTP) of this point (L{Ned}).

           @return: Destination point (L{Cartesian}).

           @raise TypeError: If B{C{delta}} is not L{Ned}.

           @example:

           >>> a = LatLon(49.66618, 3.45063)
           >>> delta = toNed(116807.681, 222.493, -0.5245)  # [N:-86126, E:-78900, D:1069]
           >>> b = a.destinationNed(delta)  # 48.88667°N, 002.37472°E

           @JSname: I{destinationPoint}.
        '''
        _xinstanceof(Ned, delta=delta)

        n, e, d = self._rotation3
        # convert NED delta to standard coordinate frame of n-vector
        dn = delta.ned
        # rotate dn to get delta in cartesian (ECEF) coordinate
        # reference frame using the rotation matrix column vectors
        dc = Cartesian(fdot(dn, n.x, e.x, d.x),
                       fdot(dn, n.y, e.y, d.y),
                       fdot(dn, n.z, e.z, d.z))

        # apply (cartesian) delta to this Cartesian to
        # obtain destination point as cartesian
        v = self.toCartesian().plus(dc)  # the plus() gives a plain vector

        return v.toLatLon(datum=self.datum, LatLon=self.classof)  # Cartesian(v.x, v.y, v.z).toLatLon(...)

    def distanceTo(self, other, radius=None, wrap=False):
        '''Approximate the distance from this to an other point.

           @arg other: The other point (L{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}).
           @kwarg wrap: Wrap/unroll the angular distance (C{bool}).

           @return: Distance (C{meter}, same units as B{C{radius}}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: Invalid B{C{radius}}.

           @example:

           >>> p = LatLon(52.205, 0.119)
           >>> q = LatLon(48.857, 2.351);
           >>> d = p.distanceTo(q)  # 404300
        '''
        self.others(other)

        a = self._N_vector.angleTo(other._N_vector, wrap=wrap)
        r = Radius(radius) if radius else self.datum.ellipsoid.R1
        return abs(a) * r

    @property_RO
    def Ecef(self):
        '''Get the ECEF I{class} (L{EcefVeness}).
        '''
        if LatLon._Ecef is None:
            from pygeodesy.ecef import EcefVeness
            LatLon._Ecef = EcefVeness
        return LatLon._Ecef

    def equals(self, other, eps=None):  # PYCHOK no cover
        '''DEPRECATED, use method C{isequalTo}.
        '''
        return self.isequalTo(other, eps=eps)

    def isequalTo(self, other, eps=None):
        '''Compare this point with an other point.

           @arg other: The other point (L{LatLon}).
           @kwarg eps: Optional margin (C{float}).

           @return: C{True} if points are identical, including
                    datum, I{ignoring height}, C{False} otherwise.

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: Invalid B{C{eps}}.

           @see: Use method L{isequalTo3} to include I{height}.

           @example:

           >>> p = LatLon(52.205, 0.119)
           >>> q = LatLon(52.205, 0.119)
           >>> e = p.isequalTo(q)  # True
        '''
        return LatLonEllipsoidalBase.isequalTo(self, other, eps=eps) \
               if self.datum == other.datum else False

#     def greatCircle(self, bearing):
#         '''Return the great circle heading on the given bearing
#            from this point.
#
#            Direction of vector is such that initial bearing vector
#            b = c × p, where p is representing this point.
#
#            @arg bearing: Bearing from this point (compass C{degrees360}).
#
#            @return: N-vector representing great circle (L{Nvector}).
#
#            @example:
#
#            >>> p = LatLon(53.3206, -1.7297)
#            >>> g = p.greatCircle(96.0)
#            >>> g.toStr()  # '(-0.794, 0.129, 0.594)'
#         '''
#         a, b, _ = self.philamheight
#         t = radians(bearing)
#
#         sa, ca, sb, cb, st, ct = sincos2(a, b, t)
#         return self._xnamed(Nvector(sb * ct - sa * cb * st,
#                                    -cb * ct - sa * sb * st,
#                                     ca * st)

#     def initialBearingTo(self, other):
#         '''Return the initial bearing (forward azimuth) from this
#            to an other point.
#
#            @arg other: The other point (L{LatLon}).
#
#            @return: Initial bearing (compass C{degrees360}).
#
#            @raise TypeError: The B{C{other}} point is not L{LatLon}.
#
#            @example:
#
#            >>> p1 = LatLon(52.205, 0.119)
#            >>> p2 = LatLon(48.857, 2.351)
#            >>> b = p1.bearingTo(p2)  # 156.2
#
#            @JSname: I{bearingTo}.
#         '''
#         self.others(other)
#
#         v1 = self.toNvector()
#         v2 = other.toNvector()
#
#         gc1 = v1.cross(v2)  # gc through v1 & v2
#         gc2 = v1.cross(_NP3)  # gc through v1 & North pole
#
#         # bearing is (signed) angle between gc1 & gc2
#         return degrees360(gc1.angleTo(gc2, vSign=v1))

    def intermediateTo(self, other, fraction, height=None):
        '''Return the point at given fraction between this and
           an other point.

           @arg other: The other point (L{LatLon}).
           @arg fraction: Fraction between both points ranging from
                          0, meaning this to 1, the other point (C{float}).
           @kwarg height: Optional height, overriding the fractional
                          height (C{meter}).

           @return: Intermediate point (L{LatLon}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @example:

           >>> p = LatLon(52.205, 0.119)
           >>> q = LatLon(48.857, 2.351)
           >>> p = p.intermediateTo(q, 0.25)  # 51.3721°N, 000.7073°E

           @JSname: I{intermediatePointTo}.
        '''
        self.others(other)

        i = other.toNvector().times(fraction).plus(
             self.toNvector().times(1 - fraction))

        if height is None:
            h = self._havg(other, f=fraction)
        else:
            h = height
        return i.toLatLon(height=h, LatLon=self.classof)  # Nvector(i.x, i.y, i.z).toLatLon(...)

    @Property_RO
    def _rotation3(self):
        '''(INTERNAL) Get the rotation matrix from n-vector coordinate frame axes.
        '''
        nv = self.toNvector()  # local (n-vector) coordinate frame

        d = nv.negate()  # down (opposite to n-vector)
        e = NorthPole.cross(nv, raiser=_pole_).unit()  # east (pointing perpendicular to the plane)
        n = e.cross(d)  # north (by right hand rule)
        return n, e, d  # matrix rows

    def toCartesian(self, **Cartesian_datum_kwds):  # PYCHOK Cartesian=Cartesian, datum=None
        '''Convert this point to an C{Nvector}-based geodetic point.

           @kwarg Cartesian_datum_kwds: Optional L{Cartesian}, B{C{datum}}
                                        and other keyword arguments, ignored
                                        if B{C{Cartesian=None}}.  Use
                                        B{C{Cartesian=...}} to override this
                                        L{Cartesian} class or specify
                                        B{C{Cartesian=None}}.

           @return: The geodetic point (L{LatLon}) or if B{C{LatLon}}
                    is C{None}, an L{Ecef9Tuple}C{(x, y, z, lat, lon,
                    height, C, M, datum)} with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{LatLon}} or other
                             B{C{Cartesian_datum_kwds}}.
        '''
        kwds = _xkwds(Cartesian_datum_kwds, Cartesian=Cartesian, datum=self.datum)
        return LatLonEllipsoidalBase.toCartesian(self, **kwds)

    def toNvector(self, **Nvector_datum_kwds):  # PYCHOK signature
        '''Convert this point to L{Nvector} components, I{including
           height}.

           @kwarg Nvector_datum_kwds: Optional L{Nvector}, B{C{datum}} or
                                      other keyword arguments, ignored if
                                      B{C{Nvector=None}}.  Use
                                      B{C{Nvector=...}} to override this
                                      L{Nvector} class or specify
                                      B{C{Nvector=None}}.

           @return: The C{n-vector} components (L{Nvector}) or a
                    L{Vector4Tuple}C{(x, y, z, h)} if B{C{Nvector}}
                    is C{None}.

           @raise TypeError: Invalid B{C{Nvector}}, B{C{datum}} or other
                             B{C{Nvector_datum_kwds}}.

           @example:

           >>> p = LatLon(45, 45)
           >>> n = p.toNvector()
           >>> n.toStr()  # [0.50, 0.50, 0.70710]
        '''
        kwds = _xkwds(Nvector_datum_kwds, Nvector=Nvector, datum=self.datum)
        return LatLonNvectorBase.toNvector(self, **kwds)


class Ned(_Named):
    '''North-Eeast-Down (NED), also known as Local Tangent Plane (LTP),
       is a vector in the local coordinate frame of a body.
    '''
#   _bearing   = None  # cached bearing (compass C{degrees360})
    _down      = None  # Down component (C{meter})
    _east      = None  # East component (C{meter})
#   _elevation = None  # cached elevation (C{degrees})
#   _length    = None  # cached length (C{float})
    _north     = None  # North component (C{meter})

    def __init__(self, north, east, down, name=NN):
        '''New North-East-Down vector.

           @arg north: North component (C{meter}).
           @arg east: East component (C{meter}).
           @arg down: Down component, normal to the surface of
                      the ellipsoid (C{meter}).
           @kwarg name: Optional name (C{str}).

           @raise ValueError: Invalid B{C{north}}, B{C{east}}
                              or B{C{down}}.
           @example:

           >>> from ellipsiodalNvector import Ned
           >>> delta = Ned(110569, 111297, 1936)
           >>> delta.toStr(prec=0)  #  [N:110569, E:111297, D:1936]
        '''
        self._north = Scalar(north=north or 0)
        self._east  = Scalar(east =east  or 0)
        self._down  = Scalar(down =down  or 0)
        if name:
            self.name = name

    def __str__(self):
        return self.toStr()

    @Property_RO
    def bearing(self):
        '''Get the bearing of this NED vector (compass C{degrees360}).
        '''
        return Bearing(atan2b(self.east, self.north))

    @Property_RO
    def down(self):
        '''Gets the Down component of this NED vector (C{meter}).
        '''
        return self._down

    @Property_RO
    def east(self):
        '''Gets the East component of this NED vector (C{meter}).
        '''
        return self._east

    @Property_RO
    def elevation(self):
        '''Get the elevation, tilt of this NED vector in degrees from
           horizontal, i.e. tangent to ellipsoid surface (C{degrees90}).
        '''
        return Float(elevation=neg(degrees90(asin1(self.down / self.length))))

    @Property_RO
    def length(self):
        '''Gets the length of this NED vector (C{meter}).
        '''
        return Float(length=hypot_(self.north, self.east, self.down))

    @Property_RO
    def ned(self):
        '''Get the C{(north, east, down)} components of the NED vector (L{Ned3Tuple}).
        '''
        r = Ned3Tuple(self.north, self.east, self.down)
        return self._xnamed(r)

    @Property_RO
    def north(self):
        '''Gets the North component of this NED vector (C{meter}).
        '''
        return self._north

    def to3ned(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{ned}.

           @return: An L{Ned3Tuple}C{(north, east, down)}.
        '''
        return self.ned

    def toRepr(self, prec=None, fmt=Fmt.SQUARE, sep=_COMMASPACE_, **unused):  # PYCHOK expected
        '''Return a string representation of this NED vector as
           length, bearing and elevation.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).
           @kwarg fmt: Optional enclosing backets format (C{str}).
           @kwarg sep: Optional separator between NEDs (C{str}).

           @return: This Ned as "[L:f, B:degrees360, E:degrees90]" (C{str}).
        '''
        from pygeodesy.dms import F_D, toDMS
        t = (fstr(self.length, prec=3 if prec is None else prec),
             toDMS(self.bearing,   form=F_D, prec=prec, ddd=0),
             toDMS(self.elevation, form=F_D, prec=prec, ddd=0))
        return _xzipairs('LBE', t, sep=sep, fmt=fmt)

    def toStr(self, prec=3, fmt=Fmt.SQUARE, sep=_COMMASPACE_):  # PYCHOK expected
        '''Return a string representation of this NED vector.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).
           @kwarg fmt: Optional enclosing backets format (C{str}).
           @kwarg sep: Optional separator between NEDs (C{str}).

           @return: This Ned as "[N:f, E:f, D:f]" (C{str}).
        '''
        t = strs(self.ned, prec=prec)
        return _xzipairs('NED', t, sep=sep, fmt=fmt)

    def toVector3d(self):
        '''Return this NED vector as a 3-d vector.

           @return: The vector(north, east, down) (L{Vector3d}).
        '''
        from pygeodesy.vector3d import Vector3d
        return Vector3d(*self.ned, name=self.name)


class Ned3Tuple(_NamedTuple):  # .ellipsoidalNvector.py
    '''3-Tuple C{(north, east, down)}, all in C{degrees}.
    '''
    _Names_ = (_north_, _east_,  _down_)
    _Units_ = ( Degrees, Degrees, Degrees)


_Nvll = LatLon(0, 0)  # Reference instance (L{LatLon})


class Nvector(NvectorBase):
    '''An n-vector is a position representation using a (unit) vector
       normal to the earth ellipsoid.  Unlike lat-/longitude points,
       n-vectors have no singularities or discontinuities.

       For many applications, n-vectors are more convenient to work
       with than other position representations like lat-/longitude,
       earth-centred earth-fixed (ECEF) vectors, UTM coordinates, etc.

       Note commonality with L{sphericalNvector.Nvector}.
    '''
    _datum = Datums.WGS84  # default datum (L{Datum})

    def __init__(self, x, y, z, h=0, datum=None, ll=None, name=NN):
        '''New n-vector normal to the earth's surface.

           @arg x: X component (C{meter}).
           @arg y: Y component (C{meter}).
           @arg z: Z component (C{meter}).
           @kwarg h: Optional height above model surface (C{meter}).
           @kwarg datum: Optional datum this n-vector is defined in
                         (L{Datum}, L{Ellipsoid}, L{Ellipsoid2} or
                         L{a_f2Tuple}).
           @kwarg ll: Optional, original latlon (C{LatLon}).
           @kwarg name: Optional name (C{str}).

           @raise TypeError: If B{C{datum}} is not a L{Datum}.

           @example:

           >>> from ellipsoidalNvector import Nvector
           >>> v = Nvector(0.5, 0.5, 0.7071, 1)
           >>> v.toLatLon()  # 45.0°N, 045.0°E, +1.00m
        '''
        NvectorBase.__init__(self, x, y, z, h=h, ll=ll, name=name)
        if datum not in (None, self._datum):
            self._datum = _ellipsoidal_datum(datum, name=name)

    @Property_RO
    def datum(self):
        '''Get this n-vector's datum (L{Datum}).
        '''
        return self._datum

    def toCartesian(self, **Cartesian_h_datum_kwds):  # PYCHOK Cartesian=Cartesian
        '''Convert this n-vector to C{Nvector}-based cartesian
           (ECEF) coordinates.

           @kwarg Cartesian_h_datum_kwds: Optional L{Cartesian}, B{C{h}},
                  B{C{datum}} and other keyword arguments, ignored if
                  B{C{Cartesian=None}}.  Use B{C{Cartesian=...}}
                  to override this L{Cartesian} class or specify
                  B{C{Cartesian=None}}.

           @return: The cartesian point (L{Cartesian}) or if B{C{Cartesian}}
                    is C{None}, an L{Ecef9Tuple}C{(x, y, z, lat, lon, height,
                    C, M, datum)} with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{Cartesian}}, B{C{h}}, B{C{datum}} or
                             other B{C{Cartesian_h_datum_kwds}}.

           @example:

           >>> v = Nvector(0.5, 0.5, 0.7071)
           >>> c = v.toCartesian()  # [3194434, 3194434, 4487327]
           >>> p = c.toLatLon()  # 45.0°N, 45.0°E
        '''
        kwds = _xkwds(Cartesian_h_datum_kwds, h=self.h, Cartesian=Cartesian,
                                                            datum=self.datum)
        return NvectorBase.toCartesian(self, **kwds)  # class or .classof

    def toLatLon(self, **LatLon_height_datum_kwds):  # PYCHOK height=None, LatLon=LatLon
        '''Convert this n-vector to an C{Nvector}-based geodetic point.

           @kwarg LatLon_height_datum_kwds: Optional L{LatLon}, B{C{height}},
                                            B{C{datum}} and other keyword arguments,
                                            ignored if C{B{LatLon}=None}.  Use
                                            B{C{LatLon=...}} to override this
                                            L{LatLon} class or set C{B{LatLon}=None}.

           @return: The geodetic point (L{LatLon}) or a L{LatLon3Tuple}C{(lat,
                    lon, height)} if B{C{LatLon}} is C{None}.

           @raise TypeError: Invalid B{C{LatLon}}, B{C{height}}, B{C{datum}}
                             or other B{C{LatLon_height_datum_kwds}}.

           @example:

           >>> v = Nvector(0.5, 0.5, 0.7071)
           >>> p = v.toLatLon()  # 45.0°N, 45.0°E
        '''
        kwds = _xkwds(LatLon_height_datum_kwds, height=self.h,
                                                LatLon=LatLon,
                                                 datum=self.datum)
        return NvectorBase.toLatLon(self, **kwds)  # class or .classof

    def unit(self, ll=None):
        '''Normalize this vector to unit length.

           @kwarg ll: Optional, original latlon (C{LatLon}).

           @return: Normalised vector (L{Nvector}).
        '''
        if self._united is None:
            u = NvectorBase.unit(self, ll=ll)
            if u.datum != self.datum:
                u._datum = self.datum
#           self._united = u._united = u
        return self._united


def meanOf(points, datum=Datums.WGS84, height=None, LatLon=LatLon,
                                                  **LatLon_kwds):
    '''Compute the geographic mean of several points.

       @arg points: Points to be averaged (L{LatLon}[]).
       @kwarg datum: Optional datum to use (L{Datum}).
       @kwarg height: Optional height at mean point, overriding
                      the mean height (C{meter}).
       @kwarg LatLon: Optional class to return the mean point
                      (L{LatLon}) or C{None}.
       @kwarg LatLon_kwds: Optional, additional B{C{LatLon}}
                           keyword arguments, ignored if
                           C{B{LatLon}=None}.

       @return: Geographic mean point and mean height (B{C{LatLon}})
                or a L{LatLon3Tuple}C{(lat, lon, height)} if
                B{C{LatLon}} is C{None}.

       @raise ValueError: Insufficient number of B{C{points}}.
    '''
    _, points = _Nvll.points2(points, closed=False)
    # geographic mean
    m = sumOf(p._N_vector for p in points)
    lat, lon, h = m._N_vector.latlonheight

    if height is not None:
        h = height
    if LatLon is None:
        r = LatLon3Tuple(lat, lon, h)
    else:
        kwds = _xkwds(LatLon_kwds, height=h, datum=datum)
        r = LatLon(lat, lon, **kwds)
    return _xnamed(r, meanOf.__name__)


def sumOf(nvectors, Vector=Nvector, h=None, **Vector_kwds):
    '''Return the vectorial sum of two or more n-vectors.

       @arg nvectors: Vectors to be added (L{Nvector}[]).
       @kwarg Vector: Optional class for the vectorial sum (L{Nvector}).
       @kwarg h: Optional height, overriding the mean height (C{meter}).
       @kwarg Vector_kwds: Optional, additional B{C{Vector}} keyword
                           arguments, ignored if C{B{Vector}=None}.

       @return: Vectorial sum (B{C{Vector}}).

       @raise VectorError: No B{C{nvectors}}.
    '''
    return _sumOf(nvectors, Vector=Vector, h=h, **Vector_kwds)


def toNed(distance, bearing, elevation, Ned=Ned, name=NN):
    '''Create an NED vector from distance, bearing and elevation
       (in local coordinate system).

       @arg distance: NED vector length (C{meter}).
       @arg bearing: NED vector bearing (compass C{degrees360}).
       @arg elevation: NED vector elevation from local coordinate
                       frame horizontal (C{degrees}).
       @kwarg Ned: Optional class to return the NED (L{Ned}) or
                   C{None}.
       @kwarg name: Optional name (C{str}).

       @return: An NED vector equivalent to this B{C{distance}},
                B{C{bearing}} and B{C{elevation}} (L{Ned}) or
                if B{C{Ned=None}}, an L{Ned3Tuple}C{(north, east,
                down)}.

       @raise ValueError: Invalid B{C{distance}}, B{C{bearing}}
                          or B{C{elevation}}.

       @JSname: I{fromDistanceBearingElevation}.
    '''
    d = Distance(distance)

    sb, cb, se, ce = sincos2d(Bearing(bearing),
                               Height(elevation=elevation))
    n  = cb * d * ce
    e  = sb * d * ce
    d *= se

    r = Ned3Tuple(n, e, -d) if Ned is None else \
              Ned(n, e, -d)
    return _xnamed(r, name)


__all__ += _ALL_OTHER(Cartesian, LatLon, Ned, Nvector,  # classes
                      meanOf, sumOf, toNed)

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
