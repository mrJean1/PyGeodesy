
# -*- coding: utf-8 -*-

u'''(INTERNAL) Private elliposiodal and spherical C{Nvector} base classes
L{LatLonNvectorBase} and L{NvectorBase} and function L{sumOf}.

Pure Python implementation of C{n-vector}-based geodesy tools for ellipsoidal
earth models, transcoded from JavaScript originals by I{(C) Chris Veness 2005-2016}
and published under the same MIT Licence**, see U{Vector-based geodesy
<https://www.Movable-Type.co.UK/scripts/latlong-vectors.html>}.
'''

from pygeodesy.basics import _isin, map1
from pygeodesy.constants import EPS, EPS0, EPS1, EPS_2, R_M, \
                               _0_0, _1_0, _2_0, _N_2_0
# from pygeodesy.datums import _spherical_datum  # from .formy
from pygeodesy.errors import IntersectionError, _ValueError, VectorError, \
                            _xattrs, _xkwds, _xkwds_pop2
from pygeodesy.fmath import fidw, hypot
from pygeodesy.fsums import Fsum, fsumf_
from pygeodesy.formy import _isequalTo,  _spherical_datum
# from pygeodesy.internals import _under  # from .named
from pygeodesy.interns import NN, _1_, _2_, _3_, _bearing_, _coincident_, \
                             _COMMASPACE_, _distance_, _h_, _insufficient_, \
                             _intersection_, _no_, _point_, _pole_, _SPACE_
from pygeodesy.latlonBase import LatLonBase,  _ALL_DOCS, _ALL_LAZY, _MODS
# from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS  # from .latlonBase
from pygeodesy.named import _xother3,  _under
from pygeodesy.namedTuples import LatLon2Tuple, PhiLam2Tuple, Trilaterate5Tuple, \
                                  Vector3Tuple, Vector4Tuple
from pygeodesy.props import deprecated_method, Property_RO, property_doc_, \
                            property_RO, property_ROnce, _update_all
from pygeodesy.streprs import Fmt, hstr, unstr
from pygeodesy.units import Bearing, Height, Radius_, Scalar
from pygeodesy.utily import atan2, sincos2d, _unrollon, _unrollon3
from pygeodesy.vector3d import Vector3d, _xyzhdlln4

from math import degrees, fabs, sqrt

__all__ = _ALL_LAZY.nvectorBase
__version__ = '25.05.12'


class NvectorBase(Vector3d):  # XXX kept private
    '''(INTERNAL) Base class for ellipsoidal and spherical C{Nvector}s.
    '''
    _datum = None         # L{Datum}, overriden
    _h     = Height(h=0)  # height (C{meter})
    _H     = NN           # height prefix (C{str}), '↑' in JS version

    def __init__(self, x_xyz, y=None, z=None, h=0, datum=None, **ll_name):
        '''New n-vector normal to the earth's surface.

           @arg x_xyz: X component of vector (C{scalar}) or (3-D) vector (C{Nvector},
                       L{Vector3d}, L{Vector3Tuple} or L{Vector4Tuple}).
           @kwarg y: Y component of vector (C{scalar}), required if B{C{x_xyz}} is
                     C{scalar} and same units as B{C{x_xyz}}, ignored otherwise.
           @kwarg z: Z component of vector (C{scalar}), like B{C{y}}.
           @kwarg h: Optional height above surface (C{meter}).
           @kwarg datum: Optional, I{pass-thru} datum (L{Datum}).
           @kwarg ll_name: Optional C{B{name}=NN} (C{str}) and optional, original
                           latlon C{B{ll}=None} (C{LatLon}).

           @raise TypeError: Non-scalar B{C{x}}, B{C{y}} or B{C{z}} coordinate or
                             B{C{x_xyz}} not an C{Nvector}, L{Vector3Tuple} or
                             L{Vector4Tuple} or invalid B{C{datum}}.
        '''
        h, d, ll, n = _xyzhdlln4(x_xyz, h, datum, **ll_name)
        Vector3d.__init__(self, x_xyz, y=y, z=z, ll=ll, name=n)
        if h:
            self.h = h
        if d is not None:
            self._datum = _spherical_datum(d, name=n)  # pass-thru

    @Property_RO
    def datum(self):
        '''Get the I{pass-thru} datum (C{Datum}) or C{None}.
        '''
        return self._datum

    @property_ROnce
    def Ecef(self):
        '''Get the ECEF I{class} (L{EcefKarney}), I{once}.
        '''
        return _MODS.ecef.EcefKarney

    @property_RO
    def ellipsoidalNvector(self):
        '''Get the C{Nvector type} iff ellipsoidal, overloaded in L{pygeodesy.ellipsoidalNvector.Nvector}.
        '''
        return False

    @property_doc_(''' the height above surface (C{meter}).''')
    def h(self):
        '''Get the height above surface (C{meter}).
        '''
        return self._h

    @h.setter  # PYCHOK setter!
    def h(self, h):
        '''Set the height above surface (C{meter}).

           @raise TypeError: If B{C{h}} invalid.

           @raise VectorError: If B{C{h}} invalid.
        '''
        h = Height(h=h, Error=VectorError)
        if self._h != h:
            _update_all(self)
            self._h = h

    @property_doc_(''' the height prefix (C{str}).''')
    def H(self):
        '''Get the height prefix (C{str}).
        '''
        return self._H

    @H.setter  # PYCHOK setter!
    def H(self, H):
        '''Set the height prefix (C{str}).
        '''
        self._H = str(H) if H else NN

    def hStr(self, prec=-2, m=NN):
        '''Return a string for the height B{C{h}}.

           @kwarg prec: Number of (decimal) digits, unstripped (C{int}).
           @kwarg m: Optional unit of the height (C{str}).

           @see: Function L{pygeodesy.hstr}.
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
        return n_xyz2latlon(self, name=self.name)

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
        return n_xyz2philam(self, name=self.name)

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

    @property_RO
    def sphericalNvector(self):
        '''Get the C{Nvector type} iff spherical, overloaded in L{pygeodesy.sphericalNvector.Nvector}.
        '''
        return False

    @deprecated_method
    def to2ab(self):  # PYCHOK no cover
        '''DEPRECATED, use property L{philam}.'''
        return self.philam

    @deprecated_method
    def to3abh(self, height=None):  # PYCHOK no cover
        '''DEPRECATED, use property L{philamheight} or C{philam.to3Tuple(B{height})}.'''
        return self.philamheight if _isin(height, None, self.h) else \
               self.philam.to3Tuple(height)

    def toCartesian(self, h=None, Cartesian=None, datum=None, **name_Cartesian_kwds):  # PYCHOK signature
        '''Convert this n-vector to C{Nvector}-based cartesian (ECEF) coordinates.

           @kwarg h: Optional height, overriding this n-vector's height (C{meter}).
           @kwarg Cartesian: Optional class to return the (ECEF) coordinates (C{Cartesian}).
           @kwarg datum: Optional datum (C{Datum}), overriding this datum.
           @kwarg name_Cartesian_kwds: Optional C{B{name}=NN} (C{str}) and optionally, additional
                       B{C{Cartesian}} keyword arguments, ignored if C{B{Cartesian} is None}.

           @return: The (ECEF) coordinates (B{C{Cartesian}}) or if C{B{Cartesian} is None}, an
                    L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)} with C{C} and C{M}
                    if available.

           @raise TypeError: Invalid B{C{Cartesian}} or B{C{name_Cartesian_kwds}} argument.

           @raise ValueError: Invalid B{C{h}}.
        '''
        if h is None:
            h = self.h
        elif not isinstance(h, Height):
            h = Height(h=h, Error=VectorError)
        _, r, v = self._toEcefDrv3(Cartesian, None, datum, h, **name_Cartesian_kwds)
        if r is None:
            r = v.toCartesian(Cartesian, **self._name1__(name_Cartesian_kwds))  # h=0
        return r

    def _toEcefDrv3(self, CC, LL, datum, h, name=NN, **unused):
        '''(INTERNAL) Helper for methods C{toCartesian} and C{toLatLon}.
        '''
        D = self.datum if _isin(datum, None, self.datum) else \
           _spherical_datum(datum, name=self.name)
        if LL is None:
            v = Vector3d(self, name=name or self.name)  # .toVector3d(norm=False)
            E = D.ellipsoid
            r = E.a_b  # Kenneth Gade eqn 22
            n = v.times_(r, r, _1_0).length
            n = (E.b / n) if n > EPS0 else _0_0
            r = E.a2_b2 * n + h  # fma
            v = v.times_(r, r, n + h)
            r = self.Ecef(D).reverse(v, M=True) if CC is None else None
        else:
            r = v = None
        return D, r, v

    @deprecated_method
    def to2ll(self):  # PYCHOK no cover
        '''DEPRECATED, use property L{latlon}.'''
        return self.latlon

    @deprecated_method
    def to3llh(self, height=None):  # PYCHOK no cover
        '''DEPRECATED, use property C{latlonheight} or C{latlon.to3Tuple(B{height})}.'''
        return self.latlonheight if _isin(height, None, self.h) else \
               self.latlon.to3Tuple(height)

    def toLatLon(self, height=None, LatLon=None, datum=None, **name_LatLon_kwds):
        '''Convert this n-vector to an C{Nvector}-based geodetic point.

           @kwarg height: Optional height, overriding this n-vector's
                          height (C{meter}).
           @kwarg LatLon: Optional class to return the geodetic point
                          (C{LatLon}) or C{None}.
           @kwarg datum: Optional, spherical datum (C{Datum}).
           @kwarg name_LatLon_kwds: Optional C{B{name}=NN} (C{str}) and optionally,
                       additional B{C{LatLon}} keyword arguments, ignored if
                       C{B{LatLon} is None}.

           @return: The geodetic point (C{LatLon}) or if C{B{LatLon} is None},
                    an L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)}
                    with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{LatLon}} or B{C{name_LatLon_kwds}} argument.

           @raise ValueError: Invalid B{C{height}}.
        '''
        h = self.h if height is None else (
            height if isinstance(height, Height) else
            Height(height, Error=VectorError))
        # use the .toCartesian() logic for better height accuracy instead of
        # r = self.Ecef(D).forward(self.lat, self.lon, height=h, M=True)
        D, r, _ = self._toEcefDrv3(None, LatLon, datum, h, **name_LatLon_kwds)
        if r is None:
            kwds = _xkwds(name_LatLon_kwds, height=h, datum=D)
            r = LatLon(self.lat, self.lon, **self._name1__(kwds))
        return r

    def toStr(self, prec=5, fmt=Fmt.PAREN, sep=_COMMASPACE_):  # PYCHOK expected
        '''Return a string representation of this n-vector.

           Height component is only included if non-zero.

           @kwarg prec: Number of (decimal) digits, unstripped (C{int}).
           @kwarg fmt: Enclosing backets format (C{str}).
           @kwarg sep: Optional separator between components (C{str}).

           @return: Comma-separated C{"(x, y, z [, h])"} enclosed in
                    B{C{fmt}} brackets (C{str}).
        '''
        t = Vector3d.toStr(self, prec=prec, fmt=NN, sep=sep)
        if self.h:
            t = sep.join((t, self.hStr()))
        return (fmt % (t,)) if fmt else t

    def toVector3d(self, norm=True):
        '''Convert this n-vector to a 3-D vector, I{ignoring height}.

           @kwarg norm: If C{True}, normalize the 3-D vector (C{bool}).

           @return: The (normalized) vector (L{Vector3d}).
        '''
        v = Vector3d.unit(self) if norm else self
        return Vector3d(v, name=self.name)

    @deprecated_method
    def to4xyzh(self, h=None):  # PYCHOK no cover
        '''DEPRECATED, use property L{xyzh} or C{xyz.to4Tuple(B{h})}.'''
        return self.xyzh if _isin(h, None, self.h) else Vector4Tuple(
               self.x, self.y, self.z, h, name=self.name)

    def unit(self, ll=None):
        '''Normalize this n-vector to unit length.

           @kwarg ll: Optional, original latlon (C{LatLon}).

           @return: Normalized vector (C{Nvector}).
        '''
        return _xattrs(Vector3d.unit(self, ll=ll), self, _under(_h_))

    @Property_RO
    def xyzh(self):
        '''Get this n-vector's components (L{Vector4Tuple}C{(x, y, z, h)})
        '''
        return self.xyz.to4Tuple(self.h)


class _N_vector_(NvectorBase):
    '''(INTERNAL) Minimal, low-overhead C{n-vector}.
    '''
    def __init__(self, x, y, z, h=0, **name):
        self._x, self._y, self._z = x, y, z
        if h:
            self._h = h
        if name:
            self.name = name


NorthPole = _N_vector_(0, 0, +1, name='NorthPole')  # North pole
SouthPole = _N_vector_(0, 0, -1, name='SouthPole')  # South pole


class LatLonNvectorBase(LatLonBase):
    '''(INTERNAL) Base class for n-vector-based ellipsoidal and spherical C{LatLon}s.
    '''

    def _update(self, updated, *attrs, **setters):  # PYCHOK _Nv=None
        '''(INTERNAL) Zap cached attributes if updated.

           @see: C{ellipsoidalNvector.LatLon} and C{sphericalNvector.LatLon} for
                 the special case of B{C{_Nv}}.
        '''
        if updated:
            _Nv, setters = _xkwds_pop2(setters, _Nv=None)
            if _Nv is not None:
                if _Nv._fromll is not None:
                    _Nv._fromll = None
                self._Nv = None
            LatLonBase._update(self, updated, *attrs, **setters)

#   def distanceTo(self, other, **kwds):  # PYCHOK no cover
#       '''I{Must be overloaded}.'''
#       self._notOverloaded(other, **kwds)

    def intersections2(self, radius1, other, radius2, **kwds):  # PYCHOK expected
        '''B{Not implemented}, throws a C{NotImplementedError} always.'''
        self._notImplemented(radius1, other, radius2, **kwds)

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

    def toNvector(self, **Nvector_and_kwds):  # PYCHOK signature
        '''Convert this point to C{Nvector} components, I{including height}.

           @kwarg Nvector_and_kwds: Optional C{Nvector} class and C{Nvector} keyword arguments,
                                    Specify C{B{Nvector}=...} to override this C{Nvector} class
                                    or use C{B{Nvector}=None}.

           @return: An C{Nvector} or if C{Nvector is None}, a L{Vector4Tuple}C{(x, y, z, h)}.

           @raise TypeError: Invalid C{Nvector} or other B{C{Nvector_and_kwds}} item.
        '''
        return LatLonBase.toNvector(self, **_xkwds(Nvector_and_kwds, Nvector=NvectorBase))

    def triangulate(self, bearing1, other, bearing2, height=None, wrap=False):  # PYCHOK signature
        '''Locate a point given this, an other point and the (initial) bearing
           from this and the other point.

           @arg bearing1: Bearing at this point (compass C{degrees360}).
           @arg other: The other point (C{LatLon}).
           @arg bearing2: Bearing at the other point (compass C{degrees360}).
           @kwarg height: Optional height at the triangulated point, overriding
                          the mean height (C{meter}).
           @kwarg wrap: If C{True}, use this and the B{C{other}} point
                        I{normalized} (C{bool}).

           @return: Triangulated point (C{LatLon}).

           @raise TypeError: Invalid B{C{other}} point.

           @raise Valuerror: Points coincide.
        '''
        return _triangulate(self, bearing1, self.others(other), bearing2,
                                  height=height, wrap=wrap, LatLon=self.classof)

    def trilaterate(self, distance1, point2, distance2, point3, distance3,
                          radius=R_M, height=None, useZ=False, wrap=False):
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
           @kwarg wrap: If C{True}, use this, B{C{point2}} and B{C{point3}}
                        I{normalized} (C{bool}).

           @return: Trilaterated point (C{LatLon}).

           @raise IntersectionError: No intersection, trilateration failed.

           @raise TypeError: Invalid B{C{point2}} or B{C{point3}}.

           @raise ValueError: Some B{C{points}} coincide or invalid B{C{distance1}},
                              B{C{distance2}}, B{C{distance3}} or B{C{radius}}.

           @see: U{Trilateration<https://WikiPedia.org/wiki/Trilateration>}, I{Veness}'
                 JavaScript U{Trilateration<https://www.Movable-Type.co.UK/scripts/
                 latlong-vectors.html>} and method C{LatLon.trilaterate5} of other,
                 non-C{Nvector LatLon} classes.
        '''
        return _trilaterate(self, distance1, self.others(point2=point2), distance2,
                                             self.others(point3=point3), distance3,
                                             radius=radius, height=height, useZ=useZ,
                                             wrap=wrap, LatLon=self.classof)

    def trilaterate5(self, distance1, point2, distance2, point3, distance3,  # PYCHOK signature
                           area=False, eps=EPS1, radius=R_M, wrap=False):
        '''B{Not implemented} for C{B{area}=True} and falls back to method C{trilaterate}.

           @return: A L{Trilaterate5Tuple}C{(min, minPoint, max, maxPoint, n)} with a
                    single trilaterated intersection C{minPoint I{is} maxPoint}, C{min
                    I{is} max} the nearest intersection margin and count C{n = 1}.

           @raise NotImplementedError: Keyword argument C{B{area}=True} not (yet) supported.

           @see: Method L{trilaterate} for other and more details.
        '''
        if area:
            self._notImplemented(area=area)

        t = _trilaterate(self, distance1, self.others(point2=point2), distance2,
                                          self.others(point3=point3), distance3,
                                          radius=radius, useZ=True, wrap=wrap,
                                          LatLon=self.classof)
        # ... and handle B{C{eps}} and C{IntersectionError}
        # like function C{.latlonBase._trilaterate5}
        d = self.distanceTo(t, radius=radius, wrap=wrap)  # PYCHOK distanceTo
        d = min(fabs(distance1 - d), fabs(distance2 - d), fabs(distance3 - d))
        if d < eps:  # min is max, minPoint is maxPoint
            return Trilaterate5Tuple(d, t, d, t, 1)  # n = 1
        t = _SPACE_(_no_(_intersection_), Fmt.PAREN(min.__name__, Fmt.f(d, prec=3)))
        raise IntersectionError(area=area, eps=eps, radius=radius, wrap=wrap, txt=t)


def n_xyz2latlon(x_xyz, y=0, z=0, **name):
    '''Convert C{n-vector} to (geodetic) lat- and longitude in C{degrees}.

       @arg x_xyz: X component (C{scalar}) or (3-D) vector (C{Nvector},
                   L{Vector3d}, L{Vector3Tuple} or L{Vector4Tuple}).
       @kwarg y: Y component of vector (C{scalar}), required if C{B{x_xyz} is
                 scalar} and same units as B{C{x_xyz}}, ignored otherwise.
       @kwarg z: Z component of vector (C{scalar}), like B{C{y}}.
       @kwarg name: Optional C{B{name}=NN} (C{str}).

       @return: A L{LatLon2Tuple}C{(lat, lon)}.

       @see: Function L{n_xyz2philam}.
    '''
    ll = map(degrees, n_xyz2philam(x_xyz, y, z))
    return LatLon2Tuple(*ll, **name)


def n_xyz2philam(x_xyz, y=0, z=0, **name):
    '''Convert C{n-vector} to (geodetic) lat- and longitude in C{radians}.

       @arg x_xyz: X component (C{scalar}) or (3-D) vector (C{Nvector},
                   L{Vector3d}, L{Vector3Tuple} or L{Vector4Tuple}).
       @kwarg y: Y component of vector (C{scalar}), required if C{B{x_xyz} is
                 scalar} and same units as B{C{x_xyz}}, ignored otherwise.
       @kwarg z: Z component of vector (C{scalar}), like B{C{y}}.
       @kwarg name: Optional C{B{name}=NN} (C{str}).

       @return: A L{PhiLam2Tuple}C{(phi, lam)}.

       @see: Function L{n_xyz2latlon}.
    '''
    try:
        x, y, z = x_xyz.xyz
    except AttributeError:
        x = x_xyz
    return PhiLam2Tuple(atan2(z, hypot(x, y)), atan2(y, x), **name)


def _nsumOf(nvs, h_None, Vector, Vector_kwds):  # .sphericalNvector, .vector3d
    '''(INTERNAL) Separated to allow callers to embellish exceptions.
    '''
    X, Y, Z, n = Fsum(), Fsum(), Fsum(), 0
    H = Fsum() if h_None is None else n
    for n, v in enumerate(nvs or ()):  # one pass
        X += v.x
        Y += v.y
        Z += v.z
        H += v.h
    if n < 1:
        raise ValueError(_SPACE_(Fmt.PARENSPACED(len=n), _insufficient_))

    x, y, z = map1(float, X, Y, Z)
    h = H.fover(n) if h_None is None else h_None
    return Vector3Tuple(x, y, z).to4Tuple(h) if Vector is None else \
                 Vector(x, y, z, **_xkwds(Vector_kwds, h=h))


def sumOf(nvectors, Vector=None, h=None, **Vector_kwds):
    '''Return the I{vectorial} sum of two or more n-vectors.

       @arg nvectors: Vectors to be added (C{Nvector}[]).
       @kwarg Vector: Optional class for the vectorial sum (C{Nvector})
                      or C{None}.
       @kwarg h: Optional height, overriding the mean height (C{meter}).
       @kwarg Vector_kwds: Optional, additional B{C{Vector}} keyword
                           arguments, ignored if C{B{Vector} is None}.

       @return: Vectorial sum (B{C{Vector}}) or a L{Vector4Tuple}C{(x, y,
                z, h)} if C{B{Vector} is None}.

       @raise VectorError: No B{C{nvectors}}.
    '''
    try:
        return _nsumOf(nvectors, h, Vector, Vector_kwds)
    except (TypeError, ValueError) as x:
        raise VectorError(nvectors=nvectors, Vector=Vector, cause=x)


def _triangulate(point1, bearing1, point2, bearing2, height=None,
                                   wrap=False, **LatLon_and_kwds):
    # (INTERNAL) Locate a point given two known points and initial
    # bearings from those points, see C{LatLon.triangulate} above

    def _gc(p, b, _i_):
        n  = p.toNvector()
        de = NorthPole.cross(n, raiser=_pole_).unit()  # east vector @ n
        dn = n.cross(de)  # north vector @ n
        s, c = sincos2d(Bearing(b, name=_bearing_ + _i_))
        dest = de.times(s)
        dnct = dn.times(c)
        d = dnct.plus(dest)  # direction vector @ n
        return n.cross(d)  # great circle point + bearing

    if wrap:
        point2 = _unrollon(point1, point2, wrap=wrap)
    if _isequalTo(point1, point2, eps=EPS):
        raise _ValueError(points=point2, wrap=wrap, txt=_coincident_)

    gc1 = _gc(point1, bearing1, _1_)  # great circle p1 + b1
    gc2 = _gc(point2, bearing2, _2_)  # great circle p2 + b2

    n = gc1.cross(gc2, raiser=_point_)  # n-vector of intersection point
    h = point1._havg(point2, h=height)
    kwds = _xkwds(LatLon_and_kwds, height=h)
    return n.toLatLon(**kwds)  # Nvector(n.x, n.y, n.z).toLatLon(...)


def _trilaterate(point1, distance1, point2, distance2, point3, distance3,
                                    radius=R_M, height=None, useZ=False,
                                    wrap=False, **LatLon_and_kwds):
    # (INTERNAL) Locate a point at given distances from
    # three other points, see LatLon.triangulate above

    def _nr2(p, d, r, _i_, *qs):  # .toNvector and angular distance squared
        for q in qs:
            if _isequalTo(p, q, eps=EPS):
                raise _ValueError(points=p, txt=_coincident_)
        return p.toNvector(), (Scalar(d, name=_distance_ + _i_) / r)**2

    p1, r     =  point1, Radius_(radius)
    p2, p3, _ = _unrollon3(p1, point2, point3, wrap)

    n1, r12 = _nr2(p1, distance1, r, _1_)
    n2, r22 = _nr2(p2, distance2, r, _2_, p1)
    n3, r32 = _nr2(p3, distance3, r, _3_, p1, p2)

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
        if fabs(j) > EPS_2:
            # courtesy of U{Carlos Freitas<https://GitHub.com/mrJean1/PyGeodesy/issues/33>}
            x = fsumf_(r12, -r22, d**2) / (d * _2_0)  # n1->intersection x- and ...
            y = fsumf_(r12, -r32, i**2, j**2, x * i * _N_2_0) / (j * _2_0)  # ... y-component
            # courtesy of U{AleixDev<https://GitHub.com/mrJean1/PyGeodesy/issues/43>}
            z = fsumf_(max(r12, r22, r32), -(x**2), -(y**2))  # XXX not just r12!
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
                kwds = _xkwds(LatLon_and_kwds, height=h)
                return n.toLatLon(**kwds)  # Nvector(n.x, n.y, n.z).toLatLon(...)

    # no intersection, d < EPS_2 or fabs(j) < EPS_2 or z < EPS
    t = _SPACE_(_no_, _intersection_, NN)
    raise IntersectionError(point1=point1, distance1=distance1,
                            point2=point2, distance2=distance2,
                            point3=point3, distance3=distance3,
                            txt=unstr(t, z=z, useZ=useZ, wrap=wrap))


__all__ += _ALL_DOCS(LatLonNvectorBase, NvectorBase, sumOf)  # classes

# **) MIT License
#
# Copyright (C) 2016-2025 -- mrJean1 at Gmail -- All Rights Reserved.
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
