
# -*- coding: utf-8 -*-

u'''(INTERNAL) Base classes for elliposiodal and spherical C{Nvector}s.

Classes L{LatLonNvectorBase} for C{n-vectorial} ellipsoidal and spherical
C{LatLon}s, class L{NvectorBase} for C{Cartesian}s and function L{sumOf}.

Pure Python implementation of C{n-vector}-based geodesy tools for
ellipsoidal earth models, transcoded from JavaScript originals by
I{(C) Chris Veness 2005-2016} and published under the same MIT Licence**,
see U{Vector-based geodesy
<https://www.Movable-Type.co.UK/scripts/latlong-vectors.html>}.
'''

from pygeodesy.basics import len2, map1
from pygeodesy.datums import _spherical_datum
from pygeodesy.errors import IntersectionError, _ValueError, \
                            _xkwds, _xkwds_pop
from pygeodesy.fmath import fidw, fsum, fsum_, hypot_
from pygeodesy.formy import n_xyz2latlon, n_xyz2philam
from pygeodesy.interns import EPS, EPS1, EPS_2, MISSING, NN, R_M, \
                             _bearing_, _coincident_, _COMMASPACE_, \
                             _distance_, _intersection_, _no_, _NorthPole_, \
                             _points_, _pole_, _SPACE_, _SouthPole_, \
                             _1_, _2_, _3_, _2_0
from pygeodesy.latlonBase import LatLonBase
from pygeodesy.lazily import _ALL_DOCS
from pygeodesy.named import notImplemented, _xother3
from pygeodesy.namedTuples import Trilaterate5Tuple, Vector3Tuple, \
                                  Vector4Tuple
from pygeodesy.props import deprecated_method, property_doc_, Property_RO
from pygeodesy.streprs import Fmt, hstr, unstr, _xattrs
from pygeodesy.units import Bearing, Height, Radius_, Scalar
from pygeodesy.utily import sincos2d
from pygeodesy.vector3d import Vector3d, VectorError, \
                               sumOf as _sumOf, _xyzhdn6

from math import fabs, sqrt  # atan2, cos, sin

__all__ = (_NorthPole_, _SouthPole_)  # constants
__version__ = '21.08.06'


class NvectorBase(Vector3d):  # XXX kept private
    '''Base class for ellipsoidal and spherical C{Nvector}s.
    '''
    _datum = None  # L{Datum}, overriden
    _h     = 0     # height (C{meter})
    _H     = NN    # height prefix (C{str}), '↑' in JS version

    def __init__(self, x, y=None, z=None, h=0, ll=None, datum=None, name=NN):
        '''New n-vector normal to the earth's surface.

           @arg x: An C{Nvector}, L{Vector3Tuple}, L{Vector4Tuple} or
                   the C{X} coordinate (C{scalar}).
           @arg y: The C{Y} coordinate (C{scalar}) if B{C{x}} C{scalar}.
           @arg z: The C{Z} coordinate (C{scalar}) if B{C{x}} C{scalar}.
           @kwarg h: Optional height above surface (C{meter}).
           @kwarg ll: Optional, original latlon (C{LatLon}).
           @kwarg datum: Optional, I{pass-thru} datum (L{Datum}).
           @kwarg name: Optional name (C{str}).

           @raise TypeError: Non-scalar B{C{x}}, B{C{y}} or B{C{z}}
                             coordinate or B{C{x}} not an C{Nvector},
                             L{Vector3Tuple} or L{Vector4Tuple} or
                             invalid B{C{datum}}.

           @example:

            >>> from pygeodesy.sphericalNvector import Nvector
            >>> v = Nvector(0.5, 0.5, 0.7071, 1)
            >>> v.toLatLon()  # 45.0°N, 045.0°E, +1.00m
        '''
        x, y, z, h, d, n = _xyzhdn6(x, y, z, h, datum, ll)
        Vector3d.__init__(self, x, y, z, ll=ll, name=name or n)
        if h:
            self.h = h
        if d:
            self._datum = _spherical_datum(d, name=self.name)  # pass-thru

    @Property_RO
    def datum(self):
        '''Get the I{pass-thru} datum (C{Datum}) or C{None}.
        '''
        return self._datum

    @Property_RO
    def Ecef(self):
        '''Get the ECEF I{class} (L{EcefKarney}), I{lazily}.
        '''
        from pygeodesy.ecef import EcefKarney
        return EcefKarney  # default

    @property_doc_(''' the height above surface (C{meter}).''')
    def h(self):
        '''Get the height above surface (C{meter}).
        '''
        return self._h

    @h.setter  # PYCHOK setter!
    def h(self, h):
        '''Set the height above surface.

           @arg h: New height (C{meter}).

           @raise TypeError: If B{C{h}} invalid.

           @raise VectorError: If B{C{h}} invalid.
        '''
        h = Height(h=h, Error=VectorError)
        self._update(h != self._h)
        self._h = h

    @property_doc_(''' the height prefix (C{str}).''')
    def H(self):
        '''Get the height prefix (C{str}).
        '''
        return self._H

    @H.setter  # PYCHOK setter!
    def H(self, H):
        '''Set the height prefix.

           @arg H: New height prefix (C{str}).
        '''
        self._H = str(H) if H else NN

    def hStr(self, prec=-2, m=NN):
        '''Return a string for the height B{C{h}}.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).
           @kwarg m: Optional unit of the height (C{str}).

           @see: Function L{hstr}.
        '''
        return NN(self.H, hstr(self.h, prec=prec, m=m))

    @Property_RO
    def isEllipsoidal(self):
        '''Check whether this n-vector is ellipsoidal (C{bool} or C{None} if unknown).
        '''
        return self.datum.isEllipsoidal if self.datum else None

    @Property_RO
    def isSpherical(self):
        '''Check whether this n-vector is spherical (C{bool} or C{None} if unknown).
        '''
        return self.datum.isSpherical if self.datum else None

    @Property_RO
    def lam(self):
        '''Get the (geodetic) longitude in C{radians} (C{float}).
        '''
        return self.philam.lam

    @Property_RO
    def lat(self):
        '''Get the (geodetic) latitude in C{degrees} (C{float}).
        '''
        return self.latlon.lat

    @Property_RO
    def latlon(self):
        '''Get the (geodetic) lat-, longitude in C{degrees} (L{LatLon2Tuple}C{(lat, lon)}).
        '''
        return n_xyz2latlon(self.x, self.y, self.z, name=self.name)

    @Property_RO
    def latlonheight(self):
        '''Get the (geodetic) lat-, longitude in C{degrees} and height (L{LatLon3Tuple}C{(lat, lon, height)}).
        '''
        return self.latlon.to3Tuple(self.h)

    @Property_RO
    def latlonheightdatum(self):
        '''Get the lat-, longitude in C{degrees} with height and datum (L{LatLon4Tuple}C{(lat, lon, height, datum)}).
        '''
        return self.latlonheight.to4Tuple(self.datum)

    @Property_RO
    def lon(self):
        '''Get the (geodetic) longitude in C{degrees} (C{float}).
        '''
        return self.latlon.lon

    @Property_RO
    def phi(self):
        '''Get the (geodetic) latitude in C{radians} (C{float}).
        '''
        return self.philam.phi

    @Property_RO
    def philam(self):
        '''Get the (geodetic) lat-, longitude in C{radians} (L{PhiLam2Tuple}C{(phi, lam)}).
        '''
        return n_xyz2philam(self.x, self.y, self.z, name=self.name)

    @Property_RO
    def philamheight(self):
        '''Get the (geodetic) lat-, longitude in C{radians} and height (L{PhiLam3Tuple}C{(phi, lam, height)}).
        '''
        return self.philam.to3Tuple(self.h)

    @Property_RO
    def philamheightdatum(self):
        '''Get the lat-, longitude in C{radians} with height and datum (L{PhiLam4Tuple}C{(phi, lam, height, datum)}).
        '''
        return self.philamheight.to4Tuple(self.datum)

    @deprecated_method
    def to2ab(self):  # PYCHOK no cover
        '''DEPRECATED, use property L{philam}.

           @return: A L{PhiLam2Tuple}C{(phi, lam)}.
        '''
        return self.philam

    @deprecated_method
    def to3abh(self, height=None):  # PYCHOK no cover
        '''DEPRECATED, use method L{philamheight} or C{philam.to3Tuple(B{height})}.

           @kwarg height: Optional height, overriding this
                          n-vector's height (C{meter}).

           @return: A L{PhiLam3Tuple}C{(phi, lam, height)}.

           @raise ValueError: Invalid B{C{height}}.
        '''
        return self.philamheight if height in (None, self.h) else \
               self.philam.to3Tuple(height)

    def toCartesian(self, h=None, Cartesian=None, datum=None, **Cartesian_kwds):
        '''Convert this n-vector to C{Nvector}-based cartesian (ECEF) coordinates.

           @kwarg h: Optional height, overriding this n-vector's height (C{meter}).
           @kwarg Cartesian: Optional class to return the (ECEF) coordinates
                             (C{Cartesian}).
           @kwarg datum: Optional datum (C{Datum}), overriding this datum.
           @kwarg Cartesian_kwds: Optional, additional B{C{Cartesian}}
                                  keyword arguments, ignored if
                                  C{B{Cartesian} is None}.

           @return: The cartesian (ECEF) coordinates (B{C{Cartesian}}) or
                    if C{B{Cartesian} is None}, an L{Ecef9Tuple}C{(x, y, z,
                    lat, lon, height, C, M, datum)} with C{C} and C{M} if
                    available.

           @raise TypeError: Invalid B{C{Cartesian}}.

           @raise ValueError: Invalid B{C{h}}.

           @example:

            >>> v = Nvector(0.5, 0.5, 0.7071)
            >>> c = v.toCartesian()  # [3194434, 3194434, 4487327]
            >>> p = c.toLatLon()  # 45.0°N, 45.0°E
        '''
        d = _spherical_datum(datum or self.datum, name=self.name)
        E =  d.ellipsoid
        h =  self.h if h is None else Height(h=h)

        x, y, z = self.x, self.y, self.z
        # Kenneth Gade eqn (22)
        n = E.b / hypot_(x * E.a_b, y * E.a_b, z)
        r = h + n * E.a2_b2

        x *= r
        y *= r
        z *= n + h

        if Cartesian is None:
            r = self.Ecef(d).reverse(x, y, z, M=True)
        else:
            kwds = _xkwds(Cartesian_kwds, datum=d)  # h=0
            r = Cartesian(x, y, z, **kwds)
        return self._xnamed(r)

    @deprecated_method
    def to2ll(self):  # PYCHOK no cover
        '''DEPRECATED, use property L{latlon}.

           @return: A L{LatLon2Tuple}C{(lat, lon)}.
        '''
        return self.latlon

    @deprecated_method
    def to3llh(self, height=None):  # PYCHOK no cover
        '''DEPRECATED, use property C{latlonheight} or C{latlon.to3Tuple}C{)}B{C{height}}C{)}.

           @kwarg height: Optional height, overriding this
                          n-vector's height (C{meter}).

           @return: A L{LatLon3Tuple}C{(lat, lon, height)}.

           @raise ValueError: Invalid B{C{height}}.
        '''
        return self.latlonheight if height in (None, self.h) else \
               self.latlon.to3Tuple(height)

    def toLatLon(self, height=None, LatLon=None, datum=None, **LatLon_kwds):
        '''Convert this n-vector to an C{Nvector}-based geodetic point.

           @kwarg height: Optional height, overriding this n-vector's
                          height (C{meter}).
           @kwarg LatLon: Optional class to return the geodetic point
                          (C{LatLon}) or C{None}.
           @kwarg datum: Optional, spherical datum (C{Datum}).
           @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                               arguments, ignored if C{B{LatLon} is None}.

           @return: The geodetic point (C{LatLon}) or if C{B{LatLon} is None},
                    an L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M,
                    datum)} with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{LatLon}}.

           @raise ValueError: Invalid B{C{height}}.

           @example:

             >>> v = Nvector(0.5, 0.5, 0.7071)
             >>> p = v.toLatLon()  # 45.0°N, 45.0°E
        '''
        d = _spherical_datum(datum or self.datum, name=self.name)
        h =  self.h if height is None else Height(height)
        # use self.Cartesian(Cartesian=None) for better accuracy of the height
        # than self.Ecef(d).forward(self.lat, self.lon, height=h, M=True)
        if LatLon is None:
            r = self.toCartesian(h=h, Cartesian=None, datum=d)
        else:
            kwds = _xkwds(LatLon_kwds, height=h, datum=d)
            r = self._xnamed(LatLon(self.lat, self.lon, **kwds))
        return r

    def toStr(self, prec=5, fmt=Fmt.PAREN, sep=_COMMASPACE_):  # PYCHOK expected
        '''Return a string representation of this n-vector.

           Height component is only included if non-zero.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).
           @kwarg fmt: Optional enclosing backets format (C{str}).
           @kwarg sep: Optional separator between components (C{str}).

           @return: Comma-separated C{"(x, y, z [, h])"} enclosed in
                    B{C{fmt}} brackets (C{str}).

           @example:

           >>> Nvector(0.5, 0.5, 0.7071).toStr()  # (0.5, 0.5, 0.7071)
           >>> Nvector(0.5, 0.5, 0.7071, 1).toStr(-3)  # (0.500, 0.500, 0.707, +1.00)
        '''
        t = Vector3d.toStr(self, prec=prec, fmt=NN, sep=sep)
        if self.h:
            t = sep.join((t, self.hStr()))
        return (fmt % (t,)) if fmt else t

    def toVector3d(self, norm=True):
        '''Convert this n-vector to a 3-D vector, I{ignoring
           the height}.

           @kwarg norm: Normalize the 3-D vector (C{bool}).

           @return: The (normalized) vector (L{Vector3d}).
        '''
        v = Vector3d.unit(self) if norm else self
        return Vector3d(v.x, v.y, v.z, name=self.name)

    @deprecated_method
    def to4xyzh(self, h=None):  # PYCHOK no cover
        '''DEPRECATED, use property L{xyzh} or C{xyz.to4Tuple(B{h})}.
        '''
        return self.xyzh if h in (None, self.h) else Vector4Tuple(
               self.x, self.y, self.z, h, name=self.name)

    def unit(self, ll=None):
        '''Normalize this n-vector to unit length.

           @kwarg ll: Optional, original latlon (C{LatLon}).

           @return: Normalized vector (C{Nvector}).
        '''
        return _xattrs(Vector3d.unit(self, ll=ll), '_h')

    @Property_RO
    def xyzh(self):
        '''Get this n-vector's components (L{Vector4Tuple}C{(x, y, z, h)})
        '''
        return self.xyz.to4Tuple(self.h)


NorthPole = NvectorBase(0, 0, +1, name=_NorthPole_)  # North pole (C{Nvector})
SouthPole = NvectorBase(0, 0, -1, name=_SouthPole_)  # South pole (C{Nvector})


class _N_vector_(NvectorBase):
    '''(INTERNAL) Minimal, low-overhead C{n-vector}.
    '''
    def __init__(self, x, y, z, h=0, name=NN):
        self._x, self._y, self._z = x, y, z
        if h:
            self._h = h
        if name:
            self.name = name


class LatLonNvectorBase(LatLonBase):
    '''(INTERNAL) Base class for n-vector-based ellipsoidal
       and spherical C{LatLon} classes.
    '''

    def _update(self, updated, *attrs, **kwds):  # PYCHOK _Nv=None
        '''(INTERNAL) Zap cached attributes if updated.

           @see: C{ellipsoidalNvector.LatLon} and C{sphericalNvector.LatLon}
                 for the special case of B{C{_Nv}}.
        '''
        if updated:
            _Nv = _xkwds_pop(kwds, _Nv=None)
            if _Nv is not None:
                if _Nv._fromll is not None:
                    _Nv._fromll = None
                self._Nv = None
            LatLonBase._update(self, updated, *attrs)

#   def distanceTo(self, other, **kwds):  # PYCHOK no cover
#       '''(INTERNAL) I{Must be overloaded}, see function C{notOverloaded}.
#       '''
#       from pygeodesy.named import notOverloaded
#       notOverloaded(self, other, **kwds)

    def intersections2(self, radius1, other, radius2, **kwds):  # PYCHOK expected
        '''B{Not implemented}, throws a C{NotImplementedError} always.
        '''
        notImplemented(self, radius1, other, radius2, **kwds)

    def others(self, *other, **name_other_up):
        '''Refined class comparison.

           @arg other: The other instance (C{LatLonNvectorBase}).
           @kwarg name_other_up: Overriding C{name=other} and C{up=1}
                                 keyword arguments.

           @return: The B{C{other}} if compatible.

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        if other:
            other0 = other[0]
            if isinstance(other0, (self.__class__, LatLonNvectorBase)):  # XXX NvectorBase?
                return other0

        other, name, up = _xother3(self, other, **name_other_up)
        if not isinstance(other, (self.__class__, LatLonNvectorBase)):  # XXX NvectorBase?
            LatLonBase.others(self, other, name=name, up=up + 1)
        return other

    def toNvector(self, Nvector=NvectorBase, **Nvector_kwds):  # PYCHOK signature
        '''Convert this point to C{Nvector} components, I{including
           height}.

           @kwarg Nvector_kwds: Optional, additional B{C{Nvector}} keyword
                                arguments, ignored if C{B{Nvector} is None}.

           @return: An B{C{Nvector}} or a L{Vector4Tuple}C{(x, y, z, h)} if
                    B{C{Nvector}} is C{None}.

           @raise TypeError: Invalid B{C{Nvector}} or B{C{Nvector_kwds}}.
        '''
        return LatLonBase.toNvector(self, Nvector=Nvector, **Nvector_kwds)

    def triangulate(self, bearing1, other, bearing2, height=None):
        '''Locate a point given this and an other point and a bearing
           at this and the other point.

           @arg bearing1: Bearing at this point (compass C{degrees360}).
           @arg other: The other point (C{LatLon}).
           @arg bearing2: Bearing at the other point (compass C{degrees360}).
           @kwarg height: Optional height at the triangulated point,
                          overriding the mean height (C{meter}).

           @return: Triangulated point (C{LatLon}).

           @raise TypeError: Invalid B{C{other}} point.

           @raise Valuerror: Points coincide.

           @example:

            >>> p = LatLon("47°18.228'N","002°34.326'W")  # Basse Castouillet
            >>> q = LatLon("47°18.664'N","002°31.717'W")  # Basse Hergo
            >>> t = p.triangulate(7, q, 295)  # 47.323667°N, 002.568501°W'
        '''
        return _triangulate(self, bearing1, self.others(other), bearing2,
                                  height=height, LatLon=self.classof)

    def trilaterate(self, distance1, point2, distance2, point3, distance3,
                          radius=R_M, height=None, useZ=False):
        '''Locate a point at given distances from this and two other points.

           @arg distance1: Distance to this point (C{meter}, same units
                           as B{C{radius}}).
           @arg point2: Second reference point (C{LatLon}).
           @arg distance2: Distance to point2 (C{meter}, same units as
                           B{C{radius}}).
           @arg point3: Third reference point (C{LatLon}).
           @arg distance3: Distance to point3 (C{meter}, same units as
                           B{C{radius}}).
           @kwarg radius: Mean earth radius (C{meter}).
           @kwarg height: Optional height at trilaterated point, overriding
                          the mean height (C{meter}, same units as B{C{radius}}).
           @kwarg useZ: Include Z component iff non-NaN, non-zero (C{bool}).

           @return: Trilaterated point (C{LatLon}).

           @raise IntersectionError: No intersection, trilateration failed.

           @raise TypeError: Invalid B{C{point2}} or B{C{point3}}.

           @raise ValueError: Some B{C{points}} coincide or invalid B{C{distance1}},
                              B{C{distance2}}, B{C{distance3}} or B{C{radius}}.

           @see: U{Trilateration<https://WikiPedia.org/wiki/Trilateration>},
                 Veness' JavaScript U{Trilateration<https://www.Movable-Type.co.UK/
                 scripts/latlong-vectors.html>} and method C{LatLon.trilaterate2}
                 of other, non-C{Nvector LatLon} classes.
        '''
        return _trilaterate(self, distance1,
                                  self.others(point2=point2), distance2,
                                  self.others(point3=point3), distance3,
                                  radius=radius, height=height, useZ=useZ,
                                  LatLon=self.classof)

    def trilaterate5(self, distance1, point2, distance2, point3, distance3,  # PYCHOK signature
                           area=False, eps=EPS1, radius=R_M, wrap=False):
        '''B{Not implemented} for C{B{area}=True} or C{B{wrap}=True}
           and falls back to method C{trilaterate} otherwise.

           @return: A L{Trilaterate5Tuple}C{(min, minPoint, max, maxPoint, n)}
                    with a single trilaterated intersection C{minPoint I{is}
                    maxPoint}, C{min I{is} max} the nearest intersection
                    margin and count C{n = 1}.

           @raise IntersectionError: No intersection, trilateration failed.

           @raise NotImplementedError: Keyword argument C{B{area}=True} or
                                       C{B{wrap}=True} not (yet) supported.

           @raise TypeError: Invalid B{C{point2}} or B{C{point3}}.

           @raise ValueError: Coincident B{C{points}} or invalid B{C{distance1}},
                              B{C{distance2}}, B{C{distance3}} or B{C{radius}}.
        '''
        if area or wrap:
            notImplemented(self, area=area, wrap=wrap)

        t = _trilaterate(self, distance1, self.others(point2=point2), distance2,
                                          self.others(point3=point3), distance3,
                                          radius=radius, height=None, useZ=True,
                                          LatLon=self.classof)
        # ... and handle B{C{eps}} and C{IntersectionError}
        # like function C{.latlonBase._trilaterate5}
        d = self.distanceTo(t, radius=radius, wrap=wrap)  # PYCHOK distanceTo
        d = min(fabs(distance1 - d), fabs(distance2 - d), fabs(distance3 - d))
        if d < eps:  # min is max, minPoint is maxPoint
            return Trilaterate5Tuple(d, t, d, t, 1)  # n = 1
        t = _SPACE_(_no_(_intersection_), Fmt.PAREN(min.__name__, Fmt.f(d, prec=3)))
        raise IntersectionError(area=area, eps=eps, wrap=wrap, txt=t)


def sumOf(nvectors, Vector=None, h=None, **Vector_kwds):
    '''Return the vectorial sum of two or more n-vectors.

       @arg nvectors: Vectors to be added (C{Nvector}[]).
       @kwarg Vector: Optional class for the vectorial sum (C{Nvector})
                      or C{None}.
       @kwarg h: Optional height, overriding the mean height (C{meter}).
       @kwarg Vector_kwds: Optional, additional B{C{Vector}} keyword
                           arguments, ignored if C{B{Vector} is None}.

       @return: Vectorial sum (B{C{Vector}}) or a L{Vector4Tuple}C{(x, y,
                z, h)} if B{C{Vector}} is C{None}.

       @raise VectorError: No B{C{nvectors}}.
    '''
    n, nvectors = len2(nvectors)
    if n < 1:
        raise VectorError(nvectors=n, txt=MISSING)

    if h is None:
        h = fsum(v.h for v in nvectors) / float(n)

    if Vector is None:
        r = _sumOf(nvectors, Vector=Vector3Tuple).to4Tuple(h)
    else:
        r = _sumOf(nvectors, Vector=Vector, h=h, **Vector_kwds)
    return r


def _triangulate(point1, bearing1, point2, bearing2, height=None,
                                 **LatLon_LatLon_kwds):
    # (INTERNAL)Locate a point given two known points and initial
    # bearings from those points, see LatLon.triangulate above

    def _gc(p, b, _i_):
        n = p.toNvector()
        de = NorthPole.cross(n, raiser=_pole_).unit()  # east vector @ n
        dn = n.cross(de)  # north vector @ n
        s, c = sincos2d(Bearing(b, name=_bearing_ + _i_))
        dest = de.times(s)
        dnct = dn.times(c)
        d = dnct.plus(dest)  # direction vector @ n
        return n.cross(d)  # great circle point + bearing

    if point1.isequalTo(point2, EPS):
        raise _ValueError(points=point2, txt=_coincident_)

    gc1 = _gc(point1, bearing1, _1_)  # great circle p1 + b1
    gc2 = _gc(point2, bearing2, _2_)  # great circle p2 + b2

    n = gc1.cross(gc2, raiser=_points_)  # n-vector of intersection point

    h = point1._havg(point2) if height is None else Height(height)
    kwds = _xkwds(LatLon_LatLon_kwds, height=h)
    return n.toLatLon(**kwds)  # Nvector(n.x, n.y, n.z).toLatLon(...)


def _trilaterate(point1, distance1, point2, distance2, point3, distance3,
                                    radius=R_M, height=None, useZ=False,
                                    **LatLon_LatLon_kwds):
    # (INTERNAL) Locate a point at given distances from
    # three other points, see LatLon.triangulate above

    def _nd2(p, d, r, _i_, *qs):  # .toNvector and angular distance squared
        for q in qs:
            if p.isequalTo(q, EPS):
                raise _ValueError(points=p, txt=_coincident_)
        return p.toNvector(), (Scalar(d, name=_distance_ + _i_) / r)**2

    r = Radius_(radius)

    n1, r12 = _nd2(point1, distance1, r, _1_)
    n2, r22 = _nd2(point2, distance2, r, _2_, point1)
    n3, r32 = _nd2(point3, distance3, r, _3_, point1, point2)

    # the following uses x,y coordinate system with origin at n1, x axis n1->n2
    y = n3.minus(n1)
    x = n2.minus(n1)
    z = None

    d = x.length  # distance n1->n2
    if d > EPS_2:  # and y.length > EPS_2:
        X = x.unit()  # unit vector in x direction n1->n2
        i = X.dot(y)  # signed magnitude of x component of n1->n3
        Y = y.minus(X.times(i)).unit()  # unit vector in y direction
        j = Y.dot(y)  # signed magnitude of y component of n1->n3
        if abs(j) > EPS_2:
            # courtesy Carlos Freitas <https://GitHub.com/mrJean1/PyGeodesy/issues/33>
            x = fsum_(r12, -r22, d**2) / (_2_0 * d)  # n1->intersection x- and ...
            y = fsum_(r12, -r32, i**2, j**2, -2 * x * i) / (_2_0 * j)  # ... y-component
            # courtesy AleixDev <https://GitHub.com/mrJean1/PyGeodesy/issues/43>
            z = fsum_(max(r12, r22, r32), -(x**2), -(y**2))  # XXX not just r12!
            if z > EPS:
                n = n1.plus(X.times(x)).plus(Y.times(y))
                if useZ:  # include Z component
                    Z = X.cross(Y)  # unit vector perpendicular to plane
                    n = n.plus(Z.times(sqrt(z)))
                if height is None:
                    h = fidw((point1.height, point2.height, point3.height),
                             map1(fabs, distance1, distance2, distance3))
                else:
                    h = Height(height)
                kwds = _xkwds(LatLon_LatLon_kwds, height=h)
                return n.toLatLon(**kwds)  # Nvector(n.x, n.y, n.z).toLatLon(...)

    # no intersection, d < EPS_2 or abs(j) < EPS_2 or z < EPS
    t = _SPACE_(_no_, _intersection_, NN)
    raise IntersectionError(point1=point1, distance1=distance1,
                            point2=point2, distance2=distance2,
                            point3=point3, distance3=distance3,
                            txt=unstr(t, z=z, useZ=useZ))


__all__ += _ALL_DOCS(LatLonNvectorBase, NvectorBase, sumOf)  # classes

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
