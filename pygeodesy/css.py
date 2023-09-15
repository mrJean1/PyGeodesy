
# -*- coding: utf-8 -*-

u'''Cassini-Soldner (CSS) projection.

Classes L{CassiniSoldner}, L{Css} and L{CSSError} requiring I{Charles Karney}'s
U{geographiclib <https://PyPI.org/project/geographiclib>} Python package to be
installed.
'''

from pygeodesy.basics import islistuple, neg, _xinstanceof, _xsubclassof
from pygeodesy.constants import _umod_360, _0_0, _0_5, _90_0
from pygeodesy.datums import _ellipsoidal_datum, _WGS84
from pygeodesy.ellipsoidalBase import LatLonEllipsoidalBase as _LLEB
from pygeodesy.errors import _ValueError, _xdatum, _xellipsoidal, \
                             _xattr, _xkwds
from pygeodesy.interns import NN, _azimuth_, _COMMASPACE_, _easting_, \
                             _lat_, _lon_, _m_, _name_, _northing_, \
                             _reciprocal_, _SPACE_
from pygeodesy.interns import _C_  # PYCHOK used!
from pygeodesy.karney import _atan2d, _copysign, _diff182, _norm2, \
                             _norm180, _signBit, _sincos2d, fabs
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.named import _NamedBase, _NamedTuple, nameof
from pygeodesy.namedTuples import EasNor2Tuple, EasNor3Tuple, \
                                  LatLon2Tuple, LatLon4Tuple, _LL4Tuple
from pygeodesy.props import deprecated_Property_RO, Property, \
                                       Property_RO, _update_all
from pygeodesy.streprs import Fmt, _fstrENH2, _fstrLL0, _xzipairs
from pygeodesy.units import Bearing, Degrees, Easting, Height, _heigHt, \
                            Lat_, Lon_, Northing, Scalar

# from math import fabs  # from .karney

__all__ = _ALL_LAZY.css
__version__ = '23.09.07'


def _CS0(cs0):
    '''(INTERNAL) Get/set default projection.
    '''
    if cs0 is None:
        cs0 = Css._CS0
        if cs0 is None:
            Css._CS0 = cs0 = CassiniSoldner(_0_0, _0_0, name='Default')
    else:
        _xinstanceof(CassiniSoldner, cs0=cs0)
    return cs0


class CSSError(_ValueError):
    '''Cassini-Soldner (CSS) conversion or other L{Css} issue.
    '''
    pass


class CassiniSoldner(_NamedBase):
    '''Cassini-Soldner projection, a Python version of I{Karney}'s C++ class U{CassiniSoldner
       <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1CassiniSoldner.html>}.
    '''
    _cb0      = _0_0
    _datum    = _WGS84  # L{Datum}
    _geodesic =  None
    _latlon0  = ()
    _meridian =  None
    _sb0      = _0_0

    def __init__(self, lat0, lon0, datum=_WGS84, name=NN):
        '''New L{CassiniSoldner} projection.

           @arg lat0: Latitude of center point (C{degrees90}).
           @arg lon0: Longitude of center point (C{degrees180}).
           @kwarg datum: Optional datum or ellipsoid (L{Datum},
                         L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg name: Optional name (C{str}).

           @raise CSSError: Invalid B{C{lat}} or B{C{lon}}.

           @example:

            >>> p = CassiniSoldner(48 + 50/60.0, 2 + 20/60.0)  # Paris
            >>> p.forward(50.9, 1.8)  # Calais
            (-37518.854545, 230003.561828)

            >>> p.reverse4(-38e3, 230e3)
            (50.899937, 1.793161, 89.580797, 0.999982)
        '''
        if datum not in (None, self._datum):
            self._datum = _xellipsoidal(datum=_ellipsoidal_datum(datum, name=name))
        if name:
            self.name = name

        self.reset(lat0, lon0)

    @Property
    def datum(self):
        '''Get the datum (L{Datum}).
        '''
        return self._datum

    @datum.setter  # PYCHOK setter!
    def datum(self, datum):
        '''Set the datum or ellipsoid (L{Datum}, L{Ellipsoid}, L{Ellipsoid2}
           or L{a_f2Tuple}) or C{None} for the default.
        '''
        d = CassiniSoldner._datum if datum is None else \
              _xellipsoidal(datum=_ellipsoidal_datum(datum, name=self.name))
        if self._datum != d:
            self._datum = d
            self.geodesic = None if self._geodesic is None else self.isExact

    def _datumatch(self, latlon):
        '''Check for matching datum ellipsoids.

           @raise CSSError: Ellipsoid mismatch of B{C{latlon}} and this projection.
        '''
        d = _xattr(latlon, datum=None)
        if d:
            _xdatum(self.datum, d, Error=CSSError)

    @Property_RO
    def equatoradius(self):
        '''Get the ellipsoid's equatorial radius, semi-axis (C{meter}).
        '''
        return self.geodesic.a

    a = equatoradius

    @Property_RO
    def flattening(self):
        '''Get the ellipsoid's flattening (C{scalar}).
        '''
        return self.geodesic.f

    f = flattening

    def forward(self, lat, lon, name=NN):
        '''Convert an (ellipsoidal) geodetic location to Cassini-Soldner
           easting and northing.

           @arg lat: Latitude of the location (C{degrees90}).
           @arg lon: Longitude of the location (C{degrees180}).
           @kwarg name: Name inlieu of this projection's name (C{str}).

           @return: An L{EasNor2Tuple}C{(easting, northing)}.

           @see: Methods L{CassiniSoldner.forward4}, L{CassiniSoldner.reverse}
                 and L{CassiniSoldner.reverse4}.

           @raise CSSError: Invalid B{C{lat}} or B{C{lon}}.
        '''
        t = self.forward6(lat, lon, name=name)
        return EasNor2Tuple(t.easting, t.northing, name=t.name)

    def forward4(self, lat, lon, name=NN):
        '''Convert an (ellipsoidal) geodetic location to Cassini-Soldner
           easting and northing.

           @arg lat: Latitude of the location (C{degrees90}).
           @arg lon: Longitude of the location (C{degrees180}).
           @kwarg name: Name inlieu of this projection's name (C{str}).

           @return: An L{EasNorAziRk4Tuple}C{(easting, northing,
                    azimuth, reciprocal)}.

           @see: Method L{CassiniSoldner.forward}, L{CassiniSoldner.forward6},
                 L{CassiniSoldner.reverse} and L{CassiniSoldner.reverse4}.

           @raise CSSError: Invalid B{C{lat}} or B{C{lon}}.
        '''
        t = self.forward6(lat, lon, name=name)
        return EasNorAziRk4Tuple(t.easting, t.northing,
                                 t.azimuth, t.reciprocal, name=t.name)

    def forward6(self, lat, lon, name=NN):
        '''Convert an (ellipsoidal) geodetic location to Cassini-Soldner
           easting and northing.

           @arg lat: Latitude of the location (C{degrees90}).
           @arg lon: Longitude of the location (C{degrees180}).
           @kwarg name: Name inlieu of this projection's name (C{str}).

           @return: An L{EasNorAziRkEqu6Tuple}C{(easting, northing,
                    azimuth, reciprocal, equatorarc, equatorazimuth)}.

           @see: Method L{CassiniSoldner.forward}, L{CassiniSoldner.forward4},
                 L{CassiniSoldner.reverse} and L{CassiniSoldner.reverse4}.

           @raise CSSError: Invalid B{C{lat}} or B{C{lon}}.
        '''
        g = self.geodesic

        lat  =  Lat_(lat, Error=CSSError)
        d, _ = _diff182(self.lon0, Lon_(lon, Error=CSSError))  # _2sum
        D    =  fabs(d)

        r = g.Inverse(lat, -D, lat, D)
        z1, a = r.azi1, (r.a12 * _0_5)
        z2, e = r.azi2, (r.s12 * _0_5)
        if e == 0:  # PYCHOK no cover
            z = _diff182(z1, z2)[0] * _0_5  # _2sum
            c = _copysign(_90_0, 90 - D)  # -90 if D > 90 else 90
            z1, z2 = c - z, c + z
        if _signBit(d):
            a, e, z2 = neg(a), neg(e), z1

        z = _norm180(z2)  # azimuth of easting direction
        p =  g.Line(lat, d, z, g.DISTANCE | g.GEODESICSCALE | g.LINE_OFF)
        # reciprocal of azimuthal northing scale
        rk = p.ArcPosition(neg(a), g.GEODESICSCALE).M21
        # rk = p._GenPosition(True, -a, g.DISTANCE)[7]

        s, c = _sincos2d(p.azi0)  # aka equatorazimuth
        sb1  = _copysign(c, lat)
        cb1  = _copysign(s, 90 - D)  # -abs(s) if D > 90 else abs(s)
        d    = _atan2d(sb1 * self._cb0 - cb1 * self._sb0,
                       cb1 * self._cb0 + sb1 * self._sb0)
        n = self._meridian.ArcPosition(d, g.DISTANCE).s12
        # n = self._meridian._GenPosition(True, d, g.DISTANCE)[4]
        return EasNorAziRkEqu6Tuple(e, n, z, rk, p.a1, p.azi0,
                                    name=name or self.name)

    @Property
    def geodesic(self):
        '''Get this projection's I{wrapped} U{geodesic.Geodesic
           <https://GeographicLib.SourceForge.io/Python/doc/code.html>}, provided
           I{Karney}'s U{geographiclib<https://PyPI.org/project/geographiclib>}
           package is installed, otherwise an I{exact} L{GeodesicExact} instance.
        '''
        g = self._geodesic
        if g is None:
            E = self.datum.ellipsoid
            try:
                g = E.geodesicw
            except ImportError:
                g = E.geodesicx
            self._geodesic = g
        return g

    @geodesic.setter  # PYCHOK setter!
    def geodesic(self, exact):
        '''Set this projection's geodesic (C{bool}) to L{GeodesicExact}
           or I{wrapped Karney}'s or C{None} for the default.

           @raise ImportError: Package U{geographiclib<https://PyPI.org/
                               project/geographiclib>} not installed or
                               not found and C{B{exact}=False}.
        '''
        E = self.datum.ellipsoid
        self._geodesic = None if exact is None else (
                         E.geodesicx if exact else E.geodesicw)
        self.reset(*self.latlon0)

    @Property_RO
    def isExact(self):
        '''Return C{True} if this projection's geodesic is L{GeodesicExact}.
        '''
        return isinstance(self.geodesic, _MODS.geodesicx.GeodesicExact)

    @Property_RO
    def lat0(self):
        '''Get the center latitude (C{degrees90}).
        '''
        return self.latlon0.lat

    @property
    def latlon0(self):
        '''Get the center lat- and longitude (L{LatLon2Tuple}C{(lat, lon)})
           in (C{degrees90}, (C{degrees180}).
        '''
        return self._latlon0

    @latlon0.setter  # PYCHOK setter!
    def latlon0(self, latlon0):
        '''Set the center lat- and longitude (ellipsoidal C{LatLon},
           L{LatLon2Tuple}, L{LatLon4Tuple} or a C{tuple} or C{list}
           with the C{lat}- and C{lon}gitude in C{degrees}).

           @raise CSSError: Invalid B{C{latlon0}} or ellipsoid mismatch
                            of B{C{latlon0}} and this projection.
        '''
        if islistuple(latlon0, 2):
            lat0, lon0 = latlon0[:2]
        else:
            try:
                lat0, lon0 = latlon0.lat, latlon0.lon
                self._datumatch(latlon0)
            except (AttributeError, TypeError, ValueError) as x:
                raise CSSError(latlon0=latlon0, cause=x)
        self.reset(lat0, lon0)

    @Property_RO
    def lon0(self):
        '''Get the center longitude (C{degrees180}).
        '''
        return self.latlon0.lon

    @deprecated_Property_RO
    def majoradius(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{equatoradius}.'''
        return self.equatoradius

    def reset(self, lat0, lon0):
        '''Set or reset the center point of this Cassini-Soldner projection.

           @arg lat0: Center point latitude (C{degrees90}).
           @arg lon0: Center point longitude (C{degrees180}).

           @raise CSSError: Invalid B{C{lat0}} or B{C{lon0}}.
        '''
        _update_all(self)

        g = self.geodesic
        self._meridian = m = g.Line(Lat_(lat0=lat0, Error=CSSError),
                                    Lon_(lon0=lon0, Error=CSSError), _0_0,
                                    g.STANDARD | g.DISTANCE_IN | g.LINE_OFF)
        self._latlon0 = LatLon2Tuple(m.lat1, m.lon1)
        s, c = _sincos2d(m.lat1)  # == self.lat0 == self.LatitudeOrigin()
        self._sb0, self._cb0 = _norm2(s * g.f1, c)

    def reverse(self, easting, northing, name=NN, LatLon=None, **LatLon_kwds):
        '''Convert a Cassini-Soldner location to (ellipsoidal) geodetic
           lat- and longitude.

           @arg easting: Easting of the location (C{meter}).
           @arg northing: Northing of the location (C{meter}).
           @kwarg name: Name inlieu of this projection's name (C{str}).
           @kwarg LatLon: Optional, ellipsoidal class to return the
                          geodetic location as (C{LatLon}) or C{None}.
           @kwarg LatLon_kwds: Optional (C{LatLon}) keyword arguments,
                               ignored if C{B{LatLon} is None}.

           @return: Geodetic location B{C{LatLon}} or if B{C{LatLon}}
                    is C{None}, a L{LatLon2Tuple}C{(lat, lon)}.

           @raise CSSError: Ellipsoidal mismatch of B{C{LatLon}} and this projection.

           @raise TypeError: Invalid B{C{LatLon}} or B{C{LatLon_kwds}}.

           @see: Method L{CassiniSoldner.reverse4}, L{CassiniSoldner.forward}.
                 L{CassiniSoldner.forward4} and L{CassiniSoldner.forward6}.
        '''
        r = self.reverse4(easting, northing, name=name)
        if LatLon is None:
            r = LatLon2Tuple(r.lat, r.lon, name=r.name)  # PYCHOK expected
        else:
            _xsubclassof(_LLEB, LatLon=LatLon)
            kwds = _xkwds(LatLon_kwds, datum=self.datum, name=r.name)
            r = LatLon(r.lat, r.lon, **kwds)  # PYCHOK expected
            self._datumatch(r)
        return r

    def reverse4(self, easting, northing, name=NN):
        '''Convert a Cassini-Soldner location to (ellipsoidal) geodetic
           lat- and longitude.

           @arg easting: Easting of the location (C{meter}).
           @arg northing: Northing of the location (C{meter}).
           @kwarg name: Name inlieu of this projection's name (C{str}).

           @return: A L{LatLonAziRk4Tuple}C{(lat, lon, azimuth, reciprocal)}.

           @see: Method L{CassiniSoldner.reverse}, L{CassiniSoldner.forward}
                 and L{CassiniSoldner.forward4}.
        '''
        g =  self.geodesic
        n =  self._meridian.Position(northing)
        r =  g.Direct(n.lat2, n.lon2, n.azi2 + _90_0, easting, outmask=g.STANDARD | g.GEODESICSCALE)
        z = _umod_360(r.azi2)  # -180 <= r.azi2 < 180 ... 0 <= z < 360
        # include z azimuth of easting direction and rk reciprocal
        # of azimuthal northing scale (see C++ member Direct() 5/6
        # <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Geodesic.html>)
        return LatLonAziRk4Tuple(r.lat2, r.lon2, z, r.M12, name=name or self.name)

    toLatLon = reverse  # XXX not reverse4

    def toRepr(self, prec=6, **unused):  # PYCHOK expected
        '''Return a string representation of this projection.

           @kwarg prec: Number of (decimal) digits, unstripped (C{int}).

           @return: This projection as C{"<classname>(lat0, lon0, ...)"}
                    (C{str}).
        '''
        return _fstrLL0(self, prec, True)

    def toStr(self, prec=6, sep=_SPACE_, **unused):  # PYCHOK expected
        '''Return a string representation of this projection.

           @kwarg prec: Number of (decimal) digits, unstripped (C{int}).
           @kwarg sep: Separator to join (C{str}).

           @return: This projection as C{"lat0 lon0"} (C{str}).
        '''
        t = _fstrLL0(self, prec, False)
        return t if sep is None else sep.join(t)


class Css(_NamedBase):
    '''Cassini-Soldner East-/Northing location.
    '''
    _CS0      =  None  # default projection (L{CassiniSoldner})
    _cs0      =  None  # projection (L{CassiniSoldner})
    _easting  = _0_0   # easting (C{float})
    _height   =  0     # height (C{meter})
    _northing = _0_0   # northing (C{float})

    def __init__(self, e, n, h=0, cs0=None, name=NN):
        '''New L{Css} Cassini-Soldner position.

           @arg e: Easting (C{meter}).
           @arg n: Northing (C{meter}).
           @kwarg h: Optional height (C{meter}).
           @kwarg cs0: Optional, the Cassini-Soldner projection
                       (L{CassiniSoldner}).
           @kwarg name: Optional name (C{str}).

           @return: The Cassini-Soldner location (L{Css}).

           @raise CSSError: If B{C{e}} or B{C{n}} is invalid.

           @raise TypeError: If B{C{cs0}} is not L{CassiniSoldner}.

           @raise ValueError: Invalid B{C{h}}.

           @example:

            >>> cs = Css(448251, 5411932.0001)
        '''
        self._cs0 = _CS0(cs0)
        self._easting  = Easting(e,  Error=CSSError)
        self._northing = Northing(n, Error=CSSError)
        if h:
            self._height = Height(h=h)
        if name:
            self.name = name

    @Property_RO
    def azi(self):
        '''Get the azimuth of easting direction (C{degrees}).
        '''
        return self.reverse4.azimuth

    azimuth = azi

    @Property
    def cs0(self):
        '''Get the projection (L{CassiniSoldner}).
        '''
        return self._cs0 or Css._CS0

    @cs0.setter  # PYCHOK setter!
    def cs0(self, cs0):
        '''Set the I{Cassini-Soldner} projection (L{CassiniSoldner}).

           @raise TypeError: Invalid B{C{cs0}}.
        '''
        cs0 = _CS0(cs0)
        if cs0 != self._cs0:
            _update_all(self)
            self._cs0 = cs0

#   def dup(self, name=NN, **e_n_h_cs0):  # PYCHOK signature
#       '''Duplicate this position with some attributes modified.
#
#          @kwarg e_n_h_cs0: Use keyword argument C{B{e}=...}, C{B{n}=...},
#                            C{B{h}=...} and/or C{B{cs0}=...} to override
#                            the current C{easting}, C{northing} C{height}
#                            or C{cs0} projectio, respectively.
#       '''
#       def _args_kwds(e=None, n=None, **kwds):
#           return (e, n), kwds
#
#       kwds = _xkwds(e_n_h_cs0, e=self.easting, n=self.northing,
#                                h=self.height, cs0=self.cs0,
#                                name=name or self.name)
#       args, kwds = _args_kwds(**kwds)
#       return self.__class__(*args, **kwds)  # .classof

    @Property_RO
    def easting(self):
        '''Get the easting (C{meter}).
        '''
        return self._easting

    @Property_RO
    def height(self):
        '''Get the height (C{meter}).
        '''
        return self._height

    @Property_RO
    def latlon(self):
        '''Get the lat- and longitude (L{LatLon2Tuple}).
        '''
        r = self.reverse4
        return LatLon2Tuple(r.lat, r.lon, name=self.name)

    @Property_RO
    def northing(self):
        '''Get the northing (C{meter}).
        '''
        return self._northing

    @Property_RO
    def reverse4(self):
        '''Get the lat, lon, azimuth and reciprocal (L{LatLonAziRk4Tuple}).
        '''
        return self.cs0.reverse4(self.easting, self.northing, name=self.name)

    @Property_RO
    def rk(self):
        '''Get the reciprocal of azimuthal northing scale (C{scalar}).
        '''
        return self.reverse4.reciprocal

    reciprocal = rk

    def toLatLon(self, LatLon=None, height=None, **LatLon_kwds):
        '''Convert this L{Css} to an (ellipsoidal) geodetic point.

           @kwarg LatLon: Optional, ellipsoidal class to return the
                          geodetic point (C{LatLon}) or C{None}.
           @kwarg height: Optional height for the point, overriding the
                          default height (C{meter}).
           @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                               arguments, ignored if C{B{LatLon} is None}.

           @return: The geodetic point (B{C{LatLon}}) or if B{C{LatLon}}
                    is C{None}, a L{LatLon4Tuple}C{(lat, lon, height,
                    datum)}.

           @raise TypeError: If B{C{LatLon}} or B{C{datum}} is not
                             ellipsoidal or invalid B{C{height}} or
                             B{C{LatLon_kwds}}.
        '''
        if LatLon:
            _xsubclassof(_LLEB, LatLon=LatLon)

        lat, lon = self.latlon
        h = _heigHt(self, height)
        return _LL4Tuple(lat, lon, h, self.cs0.datum, LatLon, LatLon_kwds,
                                                      inst=self, name=self.name)

    def toRepr(self, prec=6, fmt=Fmt.SQUARE, sep=_COMMASPACE_, m=_m_, C=False):  # PYCHOK expected
        '''Return a string representation of this L{Css} position.

           @kwarg prec: Number of (decimal) digits, unstripped (C{int}).
           @kwarg fmt: Enclosing backets format (C{str}).
           @kwarg sep: Optional separator between name:values (C{str}).
           @kwarg m: Optional unit of the height, default meter (C{str}).
           @kwarg C: Optionally, include name of projection (C{bool}).

           @return: This position as C{"[E:meter, N:meter, H:m, name:'',
                    C:Conic.Datum]"} (C{str}).
        '''
        t, T = _fstrENH2(self, prec, m)
        if self.name:
            t +=  repr(self.name),
            T += _name_,
        if C:
            t +=  self.cs0.toRepr(prec=prec),
            T += _C_,
        return _xzipairs(T, t, sep=sep, fmt=fmt)

    def toStr(self, prec=6, sep=_SPACE_, m=_m_):  # PYCHOK expected
        '''Return a string representation of this L{Css} position.

           @kwarg prec: Number of (decimal) digits, unstripped (C{int}).
           @kwarg sep: Optional separator to join (C{str}) or C{None}
                       to return an unjoined C{tuple} of C{str}s.
           @kwarg m: Height units, default C{meter} (C{str}).

           @return: This position as C{"easting nothing"} C{str} in
                    C{meter} plus C{" height"} and C{'m'} if height
                    is non-zero (C{str}).
        '''
        t, _ = _fstrENH2(self, prec, m)
        return t if sep is None else sep.join(t)


class EasNorAziRk4Tuple(_NamedTuple):
    '''4-Tuple C{(easting, northing, azimuth, reciprocal)} for the
       Cassini-Soldner location with C{easting} and C{northing} in
       C{meters} and the C{azimuth} of easting direction and
       C{reciprocal} of azimuthal northing scale, both in C{degrees}.
    '''
    _Names_ = (_easting_, _northing_, _azimuth_, _reciprocal_)
    _Units_ = ( Easting,   Northing,   Bearing,   Scalar)


class EasNorAziRkEqu6Tuple(_NamedTuple):
    '''6-Tuple C{(easting, northing, azimuth, reciprocal, equatorarc,
       equatorazimuth)} for the Cassini-Soldner location with
       C{easting} and C{northing} in C{meters} and the C{azimuth} of
       easting direction, C{reciprocal} of azimuthal northing scale,
       C{equatorarc} and C{equatorazimuth}, all in C{degrees}.
    '''
    _Names_ = EasNorAziRk4Tuple._Names_ + ('equatorarc', 'equatorazimuth')
    _Units_ = EasNorAziRk4Tuple._Units_ + ( Degrees,      Bearing)


class LatLonAziRk4Tuple(_NamedTuple):
    '''4-Tuple C{(lat, lon, azimuth, reciprocal)}, all in C{degrees}
       where C{azimuth} is the azimuth of easting direction and
       C{reciprocal} the reciprocal of azimuthal northing scale.
    '''
    _Names_ = (_lat_, _lon_, _azimuth_, _reciprocal_)
    _Units_ = ( Lat_,  Lon_,  Bearing,   Scalar)


def toCss(latlon, cs0=None, height=None, Css=Css, name=NN):
    '''Convert an (ellipsoidal) geodetic point to a Cassini-Soldner
       location.

       @arg latlon: Ellipsoidal point (C{LatLon} or L{LatLon4Tuple}).
       @kwarg cs0: Optional, the Cassini-Soldner projection to use
                   (L{CassiniSoldner}).
       @kwarg height: Optional height for the point, overriding the
                      default height (C{meter}).
       @kwarg Css: Optional class to return the location (L{Css}) or C{None}.
       @kwarg name: Optional B{C{Css}} name (C{str}).

       @return: The Cassini-Soldner location (B{C{Css}}) or an
                L{EasNor3Tuple}C{(easting, northing, height)}
                if B{C{Css}} is C{None}.

       @raise CSSError: Ellipsoidal mismatch of B{C{latlon}} and B{C{cs0}}.

       @raise ImportError: Package U{geographiclib<https://PyPI.org/
                           project/geographiclib>} not installed or
                           not found.

       @raise TypeError: If B{C{latlon}} is not ellipsoidal.
    '''
    _xinstanceof(_LLEB, LatLon4Tuple, latlon=latlon)

    cs = _CS0(cs0)
    cs._datumatch(latlon)

    c =  cs.forward4(latlon.lat, latlon.lon)
    h = _heigHt(latlon, height)
    n =  name or nameof(latlon)

    if Css is None:
        r = EasNor3Tuple(c.easting, c.northing, h, name=n)
    else:
        r = Css(c.easting, c.northing, h=h, cs0=cs, name=n)
        r._latlon = LatLon2Tuple(latlon.lat, latlon.lon, name=n)
        r._azi, r._rk = c.azimuth, c.reciprocal
    return r

# **) MIT License
#
# Copyright (C) 2016-2023 -- mrJean1 at Gmail -- All Rights Reserved.
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
