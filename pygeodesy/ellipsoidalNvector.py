
# -*- coding: utf-8 -*-

u'''Ellipsoidal, C{N-vector}-based geodesy.

Ellipsoidal classes geodetic (lat-/longitude) L{LatLon}, geocentric
(ECEF) L{Cartesian}, DEPRECATED L{Ned} and L{Nvector} and functions
L{meanOf}, L{sumOf} and DEPRECATED L{toNed}.

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
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import issubclassof, map2, _xinstanceof
from pygeodesy.datums import _ellipsoidal_datum, _spherical_datum, _WGS84
# from pygeodesy.dms import toDMS  # _MODS
from pygeodesy.ellipsoidalBase import CartesianEllipsoidalBase, _TOL_M, \
                                      LatLonEllipsoidalBase, _nearestOn
from pygeodesy.errors import _IsnotError, _xkwds
# from pygeodesy.fmath import fdot  # from .nvectorBase
from pygeodesy.interns import NN, _Nv00_, _COMMASPACE_
from pygeodesy.interns import _down_, _east_, _north_, _pole_  # PYCHOK used!
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS, _ALL_OTHER
# from pygeodesy.ltp import Ltp  # _MODS
from pygeodesy.ltpTuples import Aer as _Aer, Ned as _Ned, Ned4Tuple, \
                                sincos2d_, _xnamed
# from pygeodesy.named import _xnamed  # from .ltpTuples
from pygeodesy.nvectorBase import fdot, NorthPole, LatLonNvectorBase, \
                                  NvectorBase, sumOf as _sumOf
from pygeodesy.props import deprecated_class, deprecated_function, \
                            deprecated_method, Property_RO
from pygeodesy.streprs import Fmt, fstr, _xzipairs
from pygeodesy.units import Bearing, Distance, Height, Scalar
# from pygeodesy.utily import sincos2d_  # from .ltpTuples

__all__ = _ALL_LAZY.ellipsoidalNvector
__version__ = '22.09.20'


class Ned(_Ned):
    '''DEPRECATED, use class L{pygeodesy.Ned}.'''

    def __init__(self, north, east, down, name=NN):
        deprecated_class(self.__class__)
        _Ned.__init__(self, north, east, down, name=name)

    @deprecated_method  # PYCHOK expected
    def toRepr(self, prec=None, fmt=Fmt.SQUARE, sep=_COMMASPACE_, **unused):
        '''DEPRECATED, use class L{pygeodesy.Ned}.

           @kwarg prec: Number of (decimal) digits, unstripped (C{int}).
           @kwarg fmt: Enclosing backets format (C{str}).
           @kwarg sep: Separator between NEDs (C{str}).

           @return: This Ned as "[L:f, B:degrees360, E:degrees90]" (C{str})
                    with length or slantrange C{L}, bearing or azimuth C{B}
                    and elevation C{E}.
        '''
        dms = _MODS.dms
        t = (fstr(self.slantrange, prec=3 if prec is None else prec),
             dms.toDMS(self.azimuth,   form=dms.F_D, prec=prec, ddd=0),
             dms.toDMS(self.elevation, form=dms.F_D, prec=prec, ddd=0))
        return _xzipairs('LBE', t, sep=sep, fmt=fmt)


class Cartesian(CartesianEllipsoidalBase):
    '''Extended to convert geocentric, L{Cartesian} points to
       L{Nvector} and n-vector-based, geodetic L{LatLon}.
    '''
    @Property_RO
    def Ecef(self):
        '''Get the ECEF I{class} (L{EcefVeness}), I{lazily}.
        '''
        return _MODS.ecef.EcefVeness

    def toLatLon(self, **LatLon_and_kwds):  # PYCHOK LatLon=LatLon, datum=None
        '''Convert this cartesian to an C{Nvector}-based geodetic point.

           @kwarg LatLon_and_kwds: Optional L{LatLon}, B{C{datum}} and other
                                   keyword arguments.  Use C{B{LatLon}=...} to
                                   override this L{LatLon} class or specify
                                   C{B{LatLon} is None}.

           @return: The geodetic point (L{LatLon}) or if B{C{LatLon}} is set
                    to C{None}, an L{Ecef9Tuple}C{(x, y, z, lat, lon, height,
                    C, M, datum)} with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{LatLon_and_kwds}}.
        '''
        kwds = _xkwds(LatLon_and_kwds, LatLon=LatLon, datum=self.datum)
        return CartesianEllipsoidalBase.toLatLon(self, **kwds)

    def toNvector(self, **Nvector_and_kwds):  # PYCHOK Datums.WGS84
        '''Convert this cartesian to L{Nvector} components, I{including height}.

           @kwarg Nvector_and_kwds: Optional L{Nvector}, B{C{datum}} and other
                                    keyword arguments.  Use C{B{Nvector}=...} to
                                    override this L{Nvector} class or specify
                                    C{B{Nvector} is None}.

           @return: The C{n-vector} components (L{Nvector}) or if B{C{Nvector}}
                    is set to C{None}, a L{Vector4Tuple}C{(x, y, z, h)}

           @raise TypeError: Invalid B{C{Nvector_and_kwds}}.

           @example:

            >>> from ellipsoidalNvector import LatLon
            >>> c = Cartesian(3980581, 97, 4966825)
            >>> n = c.toNvector()  # (0.62282, 0.000002, 0.78237, +0.24)
        '''
        kwds = _xkwds(Nvector_and_kwds, Nvector=Nvector, datum=self.datum)
        return CartesianEllipsoidalBase.toNvector(self, **kwds)


class LatLon(LatLonNvectorBase, LatLonEllipsoidalBase):
    '''An n-vector-based, ellipsoidal L{LatLon} point.

       @example:

        >>> from ellipsoidalNvector import LatLon
        >>> p = LatLon(52.205, 0.119)  # height=0, datum=Datums.WGS84
    '''
    _Nv = None  # cached toNvector (L{Nvector})

    def _update(self, updated, *attrs, **setters):  # PYCHOK args
        '''(INTERNAL) Zap cached attributes if updated.
        '''
        if updated:
            LatLonNvectorBase._update(self, updated, _Nv=self._Nv)  # special case
            LatLonEllipsoidalBase._update(self, updated, *attrs, **setters)

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
#             >>> p = LatLon(53.2611, -0.7972)
#
#             >>> s = LatLon(53.3206, -1.7297)
#             >>> b = 96.0
#             >>> d = p.crossTrackDistanceTo(s, b)  # -305.7
#
#             >>> e = LatLon(53.1887, 0.1334)
#             >>> d = p.crossTrackDistanceTo(s, e)  # -307.5
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

    def deltaTo(self, other, Ned=Ned):
        '''Calculate the local delta from this to an other point.

           @note: This is a linear delta, I{unrelated} to a geodesic
                  on the ellipsoid.

           @arg other: The other point (L{LatLon}).
           @kwarg Ned: Class to use (L{pygeodesy.Ned} or
                       L{pygeodesy.Ned4Tuple}), overriding the
                       default DEPRECATED L{Ned}.

           @return: Delta from this to the other point (B{C{Ned}}).

           @raise TypeError: The B{C{other}} point is not L{LatLon} or
                             B{C{Ned}} is not L{pygeodesy.Ned} nor
                             L{pygeodesy.Ned4Tuple} nor DEPRECATED L{Ned}.

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
        ned_ = map2(dc.dot, self._rotation3)
        if issubclassof(Ned, Ned4Tuple):
            ned_ += (_MODS.ltp.Ltp(self, ecef=self.Ecef(self.datum)),)
        elif not issubclassof(Ned, _Ned):
            raise _IsnotError(Fmt.sub_class(_Ned, Ned4Tuple), Ned=Ned)
        return Ned(*ned_, name=self.name)

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
#             >>> p = LatLon(51.4778, -0.0015)
#             >>> q = p.destination(7794, 300.7)
#             >>> q.toStr()  # '51.5135°N, 000.0983°W' ?
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

           @return: Destination point (L{LatLon}).

           @raise TypeError: If B{C{delta}} is not L{pygeodesy.Ned} or
                             DEPRECATED L{Ned}.

           @example:

            >>> a = LatLon(49.66618, 3.45063)
            >>> delta = Ned(-86126, -78900, 1069)  # from Aer(222.493, -0.5245, 116807.681)
            >>> b = a.destinationNed(delta)  # 48.886669°N, 002.374721°E or 48°53′12.01″N, 002°22′29.0″E   +0.20m
        '''
        _xinstanceof(_Ned, delta=delta)

        nv, ev, dv = self._rotation3
        # convert NED delta to standard coordinate frame of n-vector
        dn = delta.ned
        # rotate dn to get delta in cartesian (ECEF) coordinate
        # reference frame using the rotation matrix column vectors
        dc = Cartesian(fdot(dn, nv.x, ev.x, dv.x),
                       fdot(dn, nv.y, ev.y, dv.y),
                       fdot(dn, nv.z, ev.z, dv.z))

        # apply (cartesian) delta to this Cartesian to obtain destination as cartesian
        v = self.toCartesian().plus(dc)
        return v.toLatLon(datum=self.datum, LatLon=self.classof)  # Cartesian(v.x, v.y, v.z).toLatLon(...)

    def distanceTo(self, other, radius=None, wrap=False):
        '''I{Approximate} the distance from this to an other point.

           @arg other: The other point (L{LatLon}).
           @kwarg radius: Mean earth radius, ellipsoid or datum
                          (C{meter}, L{Ellipsoid}, L{Ellipsoid2},
                          L{Datum} or L{a_f2Tuple}), overriding the
                          mean radius C{R1} of this point's datum..
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
        d = self.datum if radius is None else _spherical_datum(radius)
        return abs(a) * d.ellipsoid.R1  # see .utily.radians2m

    @Property_RO
    def Ecef(self):
        '''Get the ECEF I{class} (L{EcefVeness}), I{lazily}.
        '''
        return _MODS.ecef.EcefVeness

    @deprecated_method
    def equals(self, other, eps=None):  # PYCHOK no cover
        '''DEPRECATED, use method L{isequalTo}.
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

           @see: Method C{isequalTo3} to include I{height}.

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
#             >>> p = LatLon(53.3206, -1.7297)
#             >>> g = p.greatCircle(96.0)
#             >>> g.toStr()  # '(-0.794, 0.129, 0.594)'
#         '''
#         a, b, _ = self.philamheight
#         t = radians(bearing)
#
#         sa, ca, sb, cb, st, ct = sincos2_(a, b, t)
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
#             >>> p1 = LatLon(52.205, 0.119)
#             >>> p2 = LatLon(48.857, 2.351)
#             >>> b = p1.bearingTo(p2)  # 156.2
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

    def intermediateTo(self, other, fraction, height=None, **unused):  # PYCHOK wrap=False
        '''Return the point at given fraction between this and
           an other point.

           @arg other: The other point (L{LatLon}).
           @arg fraction: Fraction between both points (C{scalar},
                          0.0 at this to 1.0 at the other point.
           @kwarg height: Optional height, overriding the fractional
                          height (C{meter}).

           @return: Intermediate point (L{LatLon}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @example:

            >>> p = LatLon(52.205, 0.119)
            >>> q = LatLon(48.857, 2.351)
            >>> p = p.intermediateTo(q, 0.25)  # 51.3721°N, 000.7073°E
        '''
        self.others(other)

        f = Scalar(fraction=fraction)
        i = self.toNvector().intermediateTo(other.toNvector(), f)

        h = self._havg(other, f=f) if height is None else Height(height)
        return i.toLatLon(height=h, LatLon=self.classof)  # Nvector(i.x, i.y, i.z).toLatLon(...)

    @Property_RO
    def _rotation3(self):
        '''(INTERNAL) Get the rotation matrix from n-vector coordinate frame axes.
        '''
        nv = self.toNvector()  # local (n-vector) coordinate frame

        dv = nv.negate()  # down, opposite to n-vector
        ev = NorthPole.cross(nv, raiser=_pole_).unit()  # east, pointing perpendicular to the plane
        nv = ev.cross(dv)  # north, by right hand rule
        return nv, ev, dv  # matrix rows

    def toCartesian(self, **Cartesian_and_kwds):  # PYCHOK Cartesian=Cartesian, datum=None
        '''Convert this point to an C{Nvector}-based geodetic point.

           @kwarg Cartesian_and_kwds: Optional L{Cartesian}, B{C{datum}} and other
                                      keyword arguments.  Use C{B{Cartesian}=...}
                                      to override this L{Cartesian} class or specify
                                      C{B{Cartesian} is None}.

           @return: The geodetic point (L{Cartesian}) or if B{C{Cartesian}} is set
                    to C{None}, an L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M,
                    datum)} with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{Cartesian}} or other B{C{Cartesian_and_kwds}}.
        '''
        kwds = _xkwds(Cartesian_and_kwds, Cartesian=Cartesian, datum=self.datum)
        return LatLonEllipsoidalBase.toCartesian(self, **kwds)

    def toNvector(self, **Nvector_and_kwds):  # PYCHOK signature
        '''Convert this point to L{Nvector} components, I{including height}.

           @kwarg Nvector_and_kwds: Optional L{Nvector}, B{C{datum}} and other
                                    keyword arguments.  Use C{B{Nvector}=...}
                                    to override this L{Nvector} class or specify
                                    C{B{Nvector} is None}.

           @return: The C{n-vector} components (L{Nvector}) or if B{C{Nvector}}
                    is set to C{None}, a L{Vector4Tuple}C{(x, y, z, h)}.

           @raise TypeError: Invalid B{C{Nvector}} or other B{C{Nvector_and_kwds}}.

           @example:

            >>> p = LatLon(45, 45)
            >>> n = p.toNvector()
            >>> n.toStr()  # [0.50, 0.50, 0.70710]
        '''
        kwds = _xkwds(Nvector_and_kwds, Nvector=Nvector, datum=self.datum)
        return LatLonNvectorBase.toNvector(self, **kwds)


_Nvll = LatLon(0, 0, name=_Nv00_)  # reference instance (L{LatLon})


class Nvector(NvectorBase):
    '''An n-vector is a position representation using a (unit) vector
       normal to the earth ellipsoid.  Unlike lat-/longitude points,
       n-vectors have no singularities or discontinuities.

       For many applications, n-vectors are more convenient to work
       with than other position representations like lat-/longitude,
       earth-centred earth-fixed (ECEF) vectors, UTM coordinates, etc.

       Note commonality with L{sphericalNvector.Nvector}.
    '''
    _datum = _WGS84  # default datum (L{Datum})

    def __init__(self, x_xyz, y=None, z=None, h=0, datum=None, ll=None, name=NN):
        '''New n-vector normal to the earth's surface.

           @arg x_xyz: X component of vector (C{scalar}) or (3-D) vector
                       (C{Nvector}, L{Vector3d}, L{Vector3Tuple} or
                       L{Vector4Tuple}).
           @kwarg y: Y component of vector (C{scalar}), ignored if B{C{x_xyz}}
                     is not C{scalar}, otherwise same units as B{C{x_xyz}}.
           @kwarg z: Z component of vector (C{scalar}), ignored if B{C{x_xyz}}
                     is not C{scalar}, otherwise same units as B{C{x_xyz}}.
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
        NvectorBase.__init__(self, x_xyz, y=y, z=z, h=h, ll=ll, name=name)
        if datum not in (None, self._datum):
            self._datum = _ellipsoidal_datum(datum, name=name)

    @Property_RO
    def datum(self):
        '''Get this n-vector's datum (L{Datum}).
        '''
        return self._datum

    def toCartesian(self, **Cartesian_and_kwds):  # PYCHOK Cartesian=Cartesian
        '''Convert this n-vector to C{Nvector}-based cartesian (ECEF) coordinates.

           @kwarg Cartesian_and_kwds: Optional L{Cartesian}, B{C{h}}, B{C{datum}} and
                                      other keyword arguments.  Use C{B{Cartesian}=...}
                                      to override this L{Cartesian} class or specify
                                      C{B{Cartesian} is None}.

           @return: The cartesian point (L{Cartesian}) or if B{C{Cartesian}} is set
                    to C{None}, an L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M,
                    datum)} with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{Cartesian_and_kwds}}.

           @example:

            >>> v = Nvector(0.5, 0.5, 0.7071)
            >>> c = v.toCartesian()  # [3194434, 3194434, 4487327]
            >>> p = c.toLatLon()  # 45.0°N, 45.0°E
        '''
        kwds = _xkwds(Cartesian_and_kwds, h=self.h, Cartesian=Cartesian,
                                                        datum=self.datum)
        return NvectorBase.toCartesian(self, **kwds)  # class or .classof

    def toLatLon(self, **LatLon_and_kwds):  # PYCHOK height=None, LatLon=LatLon
        '''Convert this n-vector to an C{Nvector}-based geodetic point.

           @kwarg LatLon_and_kwds: Optional L{LatLon}, B{C{height}}, B{C{datum}}
                                   and other keyword arguments.  Use C{B{LatLon}=...}
                                   to override this L{LatLon} class or specify
                                   C{B{LatLon} is None}.

           @return: The geodetic point (L{LatLon}) or if B{C{LatLon}} is set
                    to C{None}, an L{Ecef9Tuple}C{(x, y, z, lat, lon, height,
                    C, M, datum)} with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{LatLon_and_kwds}}.

           @example:

            >>> v = Nvector(0.5, 0.5, 0.7071)
            >>> p = v.toLatLon()  # 45.0°N, 45.0°E
        '''
        kwds = _xkwds(LatLon_and_kwds, height=self.h, datum=self.datum, LatLon=LatLon)
        return NvectorBase.toLatLon(self, **kwds)  # class or .classof

    def unit(self, ll=None):
        '''Normalize this vector to unit length.

           @kwarg ll: Optional, original latlon (C{LatLon}).

           @return: Normalised vector (L{Nvector}).
        '''
        u = NvectorBase.unit(self, ll=ll)
        if u.datum != self.datum:
            u._update(False, datum=self.datum)
        return u


def meanOf(points, datum=_WGS84, height=None, LatLon=LatLon,
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
                           C{B{LatLon} is None}.

       @return: Geographic mean point and mean height (B{C{LatLon}})
                or if B{C{LatLon}} is C{None}, an L{Ecef9Tuple}C{(x,
                y, z, lat, lon, height, C, M, datum)} with C{C} and
                C{M} if available.

       @raise ValueError: Insufficient number of B{C{points}}.
    '''
    Ps = _Nvll.PointsIter(points)
    # geographic mean
    m = sumOf(p._N_vector for p in Ps.iterate(closed=False))
    kwds = _xkwds(LatLon_kwds, height=height, datum=datum,
                               LatLon=LatLon, name=meanOf.__name__)
    return m.toLatLon(**kwds)


def nearestOn(point, point1, point2, within=True, height=None, wrap=False,
              equidistant=None, tol=_TOL_M, LatLon=LatLon, **LatLon_kwds):
    '''Iteratively locate the closest point on the geodesic between
       two other (ellipsoidal) points.

       @arg point: Reference point (C{LatLon}).
       @arg point1: Start point of the geodesic (C{LatLon}).
       @arg point2: End point of the geodesic (C{LatLon}).
       @kwarg within: If C{True} return the closest point I{between}
                      B{C{point1}} and B{C{point2}}, otherwise the
                      closest point elsewhere on the geodesic (C{bool}).
       @kwarg height: Optional height for the closest point (C{meter},
                      conventionally) or C{None} or C{False} for the
                      interpolated height.  If C{False}, the closest
                      takes the heights of the points into account.
       @kwarg wrap: Wrap and unroll longitudes (C{bool}).
       @kwarg equidistant: An azimuthal equidistant projection (I{class}
                           or function L{pygeodesy.equidistant}) or C{None}
                           for the preferred C{B{point}.Equidistant}.
       @kwarg tol: Convergence tolerance (C{meter}).
       @kwarg LatLon: Optional class to return the closest point
                      (L{LatLon}) or C{None}.
       @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                           arguments, ignored if C{B{LatLon} is None}.

       @return: Closest point, a B{C{LatLon}} instance or if C{B{LatLon}
                is None}, a L{LatLon4Tuple}C{(lat, lon, height, datum)}.

       @raise ImportError: Package U{geographiclib
                           <https://PyPI.org/project/geographiclib>}
                           not installed or not found.

       @raise TypeError: Invalid or non-ellipsoidal B{C{point}}, B{C{point1}}
                         or B{C{point2}} or invalid B{C{equidistant}}.

       @raise ValueError: No convergence for the B{C{tol}}.

       @see: U{The B{ellipsoidal} case<https://GIS.StackExchange.com/questions/48937/
             calculating-intersection-of-two-circles>} and U{Karney's paper
             <https://ArXiv.org/pdf/1102.1215.pdf>}, pp 20-21, section B{14. MARITIME
             BOUNDARIES} for more details about the iteration algorithm.
    '''
    return _nearestOn(point, point1, point2, within=within, height=height, wrap=wrap,
                      equidistant=equidistant, tol=tol, LatLon=LatLon, **LatLon_kwds)


def sumOf(nvectors, Vector=Nvector, h=None, **Vector_kwds):
    '''Return the vectorial sum of two or more n-vectors.

       @arg nvectors: Vectors to be added (L{Nvector}[]).
       @kwarg Vector: Optional class for the vectorial sum (L{Nvector}).
       @kwarg h: Optional height, overriding the mean height (C{meter}).
       @kwarg Vector_kwds: Optional, additional B{C{Vector}} keyword
                           arguments, ignored if C{B{Vector} is None}.

       @return: Vectorial sum (B{C{Vector}}).

       @raise VectorError: No B{C{nvectors}}.
    '''
    return _sumOf(nvectors, Vector=Vector, h=h, **Vector_kwds)


@deprecated_function
def toNed(distance, bearing, elevation, Ned=Ned, name=NN):
    '''DEPRECATED, use L{pygeodesy.Aer}C{(bearing, elevation,
       distance).xyzLocal.toNed(B{Ned}, name=B{name})} or
       L{XyzLocal}C{(pygeodesy.Aer(bearing, elevation,
       distance)).toNed(B{Ned}, name=B{name})}.

       Create an NED vector from distance, bearing and elevation
       (in local coordinate system).

       @arg distance: NED vector length (C{meter}).
       @arg bearing: NED vector bearing (compass C{degrees360}).
       @arg elevation: NED vector elevation from local coordinate
                       frame horizontal (C{degrees}).
       @kwarg Ned: Optional class to return the NED (C{Ned}) or
                   C{None}.
       @kwarg name: Optional name (C{str}).

       @return: An NED vector equivalent to this B{C{distance}},
                B{C{bearing}} and B{C{elevation}} (DEPRECATED L{Ned})
                or a DEPRECATED L{Ned3Tuple}C{(north, east, down)}
                if C{B{Ned} is None}.

       @raise ValueError: Invalid B{C{distance}}, B{C{bearing}}
                          or B{C{elevation}}.
    '''
    if True:  # use new Aer class
        n, e, d, _ = _Aer(bearing, elevation, distance).xyz4
    else:  # DEPRECATED
        d = Distance(distance)

        sb, cb, se, ce = sincos2d_(Bearing(bearing),
                                   Height(elevation=elevation))
        n  = cb * d * ce
        e  = sb * d * ce
        d *= se

    r = _MODS.deprecated.Ned3Tuple(n, e, -d) if Ned is None else \
                               Ned(n, e, -d)
    return _xnamed(r, name)


__all__ += _ALL_OTHER(Cartesian, LatLon, Ned, Nvector,  # classes
                      meanOf, sumOf, toNed)

# **) MIT License
#
# Copyright (C) 2016-2022 -- mrJean1 at Gmail -- All Rights Reserved.
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
