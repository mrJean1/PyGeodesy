
# -*- coding: utf-8 -*-

u'''Equidistant, Equal-Area, and other Azimuthal projections.

Classes L{Equidistant}, L{EquidistantExact}, L{EquidistantGeodSolve},
L{EquidistantKarney}, L{Gnomonic}, L{GnomonicExact}, L{GnomonicKarney},
L{LambertEqualArea}, L{Orthographic} and L{Stereographic}, classes
L{AzimuthalError}, L{Azimuthal7Tuple} and functions L{equidistant}
and L{gnomonic}.

L{EquidistantExact} and L{GnomonicExact} are based on exact geodesic classes
L{GeodesicExact} and L{GeodesicLineExact}, Python versions of I{Charles Karney}'s
C++ original U{GeodesicExact<https://GeographicLib.SourceForge.io/C++/doc/
classGeographicLib_1_1GeodesicExact.html>}, respectively U{GeodesicLineExact
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1GeodesicLineExact.html>}.

Using L{EquidistantGeodSolve} requires I{Karney}'s utility U{GeodSolve
<https://GeographicLib.SourceForge.io/C++/doc/GeodSolve.1.html>} to be
executable and set in env variable C{PYGEODESY_GEODSOLVE}, see module
L{geodsolve} for more details.

L{EquidistantKarney} and L{GnomonicKarney} require I{Karney}'s Python package
U{geographiclib<https://PyPI.org/project/geographiclib>} to be installed.

Other azimuthal classes implement only (***) U{Snyder's FORMULAS FOR THE SPHERE
<https://Pubs.USGS.gov/pp/1395/report.pdf>} and use those for any datum,
spherical and ellipsoidal.  The radius used for the latter is the ellipsoid's
I{mean radius of curvature} at the latitude of the projection center point.  For
further justification, see the first paragraph under U{Snyder's FORMULAS FOR THE
ELLIPSOID, page 197<https://Pubs.USGS.gov/pp/1395/report.pdf>}.

Page numbers in C{Snyder} references apply to U{John P. Snyder, "Map Projections
-- A Working Manual", 1987<https://Pubs.USGS.gov/pp/1395/report.pdf>}.

See also U{here<https://WikiPedia.org/wiki/Azimuthal_equidistant_projection>},
especially the U{Comparison of the Azimuthal equidistant projection and some
azimuthal projections centred on 90° N at the same scale, ordered by projection
altitude in Earth radii<https://WikiPedia.org/wiki/Azimuthal_equidistant_projection
#/media/File:Comparison_azimuthal_projections.svg>}.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import isnon0, _umod_360, _xinstanceof
from pygeodesy.ellipsoidalBase import LatLonEllipsoidalBase as _LLEB
from pygeodesy.datums import _spherical_datum, _WGS84
from pygeodesy.errors import _ValueError, _xdatum, _xkwds
from pygeodesy.fmath import euclid, Fsum, hypot
# from pygeodesy.fsums import Fsum  # from .fmath
# from pygeodesy.formy import antipode  # from latlonBase
from pygeodesy.interns import EPS, EPS0, EPS1, _EPStol, NAN, NN, \
                             _azimuth_, _datum_, _lat_, _lon_, \
                             _no_, _scale_, _SPACE_, _x_, _y_, \
                             _0_0, _0_1, _0_5, _1_0, _N_1_0, _2_0
from pygeodesy.karney import _norm180
from pygeodesy.latlonBase import antipode, LatLonBase as _LLB
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _FOR_DOCS
from pygeodesy.named import _NamedBase, _NamedTuple, notOverloaded, _Pass
from pygeodesy.namedTuples import LatLon2Tuple, LatLon4Tuple
from pygeodesy.props import deprecated_Property_RO, Property_RO, \
                            property_doc_, _update_all
from pygeodesy.streprs import Fmt, _fstrLL0
from pygeodesy.units import Bearing, Easting, Lat_, Lon_, \
                            Northing, Scalar, Scalar_
from pygeodesy.utily import asin1, atan2b, atan2d, sincos2, \
                            sincos2d, sincos2d_

from math import acos, atan, atan2, degrees, sin, sqrt

__all__ = _ALL_LAZY.azimuthal
__version__ = '22.07.25'

_EPS_K         = _EPStol * _0_1  # Karney's eps_ or _EPSmin * _0_1?
_over_horizon_ = 'over horizon'
_TRIPS         =  21  # numit, 4 sufficient


def _enzh4(x, y, *h):
    '''(INTERNAL) Return 4-tuple (easting, northing, azimuth, hypot).
    '''
    e = Easting( x=x)
    n = Northing(y=y)
    z = atan2b(e, n)  # (x, y) for azimuth from true North
    return e, n, z, (h[0] if h else hypot(e, n))


class _AzimuthalBase(_NamedBase):
    '''(INTERNAL) Base class for azimuthal projections.

       @see: I{Karney}'s C++ class U{AzimuthalEquidistant<https://GeographicLib.SourceForge.io/
       html/classGeographicLib_1_1AzimuthalEquidistant.html>} and U{Gnomonic
       <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Gnomonic.html>} or the
       C{PyGeodesy} versions thereof L{EquidistantKarney} respectively L{GnomonicKarney}.
    '''
    _datum   = _WGS84  # L{Datum}
    _latlon0 =  LatLon2Tuple(_0_0, _0_0)  # lat0, lon0 (L{LatLon2Tuple})
    _sc0     = _0_0, _1_0  # 2-Tuple C{sincos2d(lat0)}

    def __init__(self, lat0, lon0, datum=None, name=NN):
        '''New azimuthal projection.

           @arg lat0: Latitude of the center point (C{degrees90}).
           @arg lon0: Longitude of the center point (C{degrees180}).
           @kwarg datum: Optional datum or ellipsoid (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}) or I{scalar} earth
                         radius (C{meter}).
           @kwarg name: Optional name for the projection (C{str}).

           @raise AzimuthalError: Invalid B{C{lat0}} or B{C{lon0}} or (spherical) B{C{datum}}.

           @raise TypeError: Invalid B{C{datum}}.
       '''
        if datum not in (None, self._datum):
            self._datum = _spherical_datum(datum, name=name)
        if name:
            self.name = name

        if lat0 or lon0:  # often both 0
            self._reset(lat0, lon0)

    @Property_RO
    def datum(self):
        '''Get the datum (L{Datum}).
        '''
        return self._datum

    @Property_RO
    def equatoradius(self):
        '''Get the geodesic's equatorial radius, semi-axis (C{meter}).
        '''
        return self.datum.ellipsoid.a

    @Property_RO
    def flattening(self):
        '''Get the geodesic's flattening (C{float}).
        '''
        return self.datum.ellipsoid.f

    def forward(self, lat, lon, name=NN):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}, see function C{notOverloaded}.
        '''
        notOverloaded(self, lat, lon, name=name)

    def _forward(self, lat, lon, name, _k_t_2):
        '''(INTERNAL) Azimuthal (spherical) forward C{lat, lon} to C{x, y}.
        '''
        lat, lon = Lat_(lat), Lon_(lon)
        sa, ca, sb, cb = sincos2d_(lat, lon - self.lon0)
        s0, c0 = self._sc0

        cb  *= ca
        k, t = _k_t_2(s0 * sa + c0 * cb)
        if t:
            r = k * self.radius
            y = r * (c0 * sa - s0 * cb)
            e, n, z, _ = _enzh4(r * sb * ca, y, None)
        else:  # 0 or 180
            e = n = z = _0_0

        t = Azimuthal7Tuple(e, n, lat, lon, z, k, self.datum,
                                  name=name or self.name)
        return t

    @Property_RO
    def lat0(self):
        '''Get the center latitude (C{degrees90}).
        '''
        return self._latlon0.lat

    @property
    def latlon0(self):
        '''Get the center lat- and longitude (L{LatLon2Tuple}C{(lat, lon)}) in (C{degrees90}, C{degrees180}).
        '''
        return self._latlon0

    @latlon0.setter  # PYCHOK setter!
    def latlon0(self, latlon0):
        '''Set the center lat- and longitude (C{LatLon}, L{LatLon2Tuple} or L{LatLon4Tuple}).

           @raise AzimuthalError: Invalid B{C{lat0}} or B{C{lon0}} or ellipsoidal mismatch
                                  of B{C{latlon0}} and this projection.
        '''
        B = _LLEB if self.datum.isEllipsoidal else _LLB
        _xinstanceof(B, LatLon2Tuple, LatLon4Tuple, latlon0=latlon0)
        if hasattr(latlon0, _datum_):
            _xdatum(self.datum, latlon0.datum, Error=AzimuthalError)
        self.reset(latlon0.lat, latlon0.lon)

    @Property_RO
    def lon0(self):
        '''Get the center longitude (C{degrees180}).
        '''
        return self._latlon0.lon

    @deprecated_Property_RO
    def majoradius(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{equatoradius}.'''
        return self.equatoradius

    @Property_RO
    def radius(self):
        '''Get this projection's mean radius of curvature (C{meter}).
        '''
        return self.datum.ellipsoid.rocMean(self.lat0)

    def reset(self, lat0, lon0):
        '''Set or reset the center point of this azimuthal projection.

           @arg lat0: Center point latitude (C{degrees90}).
           @arg lon0: Center point longitude (C{degrees180}).

           @raise AzimuthalError: Invalid B{C{lat0}} or B{C{lon0}}.
        '''
        _update_all(self)  # zap caches
        self._reset(lat0, lon0)

    def _reset(self, lat0, lon0):
        '''(INTERNAL) Update the center point.
        '''
        self._latlon0 = LatLon2Tuple(Lat_(lat0=lat0, Error=AzimuthalError),
                                     Lon_(lon0=lon0, Error=AzimuthalError))
        self._sc0     = sincos2d(self.lat0)

    def reverse(self, x, y, name=NN, **LatLon_and_kwds):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}, see function C{notOverloaded}.
        '''
        notOverloaded(self, x, y, name=name, **LatLon_and_kwds)

    def _reverse(self, x, y, name, _c, lea, LatLon=None, **LatLon_kwds):
        '''(INTERNAL) Azimuthal (spherical) reverse C{x, y} to C{lat, lon}.
        '''
        e, n, z, r = _enzh4(x, y)

        c = _c(r / self.radius)
        if c is None:
            lat, lon = self.latlon0
            k, z = _1_0, _0_0
        else:
            s0, c0 = self._sc0
            sc, cc = sincos2(c)
            k = c  / sc
            s = s0 * cc
            if r > EPS0:
                s += c0 * sc * (n / r)
            lat = degrees(asin1(s))
            if lea or abs(c0) > EPS:
                d = atan2d(e * sc, r * c0 * cc - n * s0 * sc)
            else:
                d = atan2d(e, (n if s0 < 0 else -n))
            lon = _norm180(self.lon0 + d)

        if LatLon is None:
            t = Azimuthal7Tuple(e, n, lat, lon, z, k, self.datum,
                                      name=name or self.name)
        else:
            t = self._toLatLon(lat, lon, LatLon, LatLon_kwds, name)
        return t

    def _reverse2(self, x, y):
        '''(INTERNAL) See iterating functions .ellipsoidalBaseDI._intersect3,
           .ellipsoidalBaseDI._intersects2 and .ellipsoidalBaseDI._nearestOne.
        '''
        t = self.reverse(x, y)  # LatLon=None
        d = euclid(t.lat - self.lat0, t.lon - self.lon0)  # degrees
        return t, d

    def _toLatLon(self, lat, lon, LatLon, LatLon_kwds, name):
        '''(INTERNAL) Check B{C{LatLon}} and return an instance.
        '''
        kwds = _xkwds(LatLon_kwds, datum=self.datum)
        r =  self._xnamed(LatLon(lat, lon, **kwds), name=name)  # handle .classof
        B = _LLEB if self.datum.isEllipsoidal else _LLB
        _xinstanceof(B, LatLon=r)
        return r

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


class AzimuthalError(_ValueError):
    '''An azimuthal L{Equidistant}, L{EquidistantKarney}, L{Gnomonic},
       L{LambertEqualArea}, L{Orthographic}, L{Stereographic} or
       L{Azimuthal7Tuple} issue.
    '''
    pass


class Azimuthal7Tuple(_NamedTuple):
    '''7-Tuple C{(x, y, lat, lon, azimuth, scale, datum)}, in C{meter}, C{meter},
       C{degrees90}, C{degrees180}, compass C{degrees}, C{scalar} and C{Datum}
       where C{(x, y)} is the easting and northing of a projected point, C{(lat,
       lon)} the geodetic location, C{azimuth} the azimuth, clockwise from true
       North and C{scale} is the projection scale, either C{1 / reciprocal} or
       C{1} or C{-1} in the L{Equidistant} case.
    '''
    _Names_ = (_x_,     _y_,      _lat_, _lon_, _azimuth_, _scale_, _datum_)
    _Units_ = ( Easting, Northing, Lat_,  Lon_,  Bearing,   Scalar, _Pass)

    def antipodal(self, azimuth=None):
        '''Return this tuple with the antipodal C{lat} and C{lon}.

           @kwarg azimuth: Optional azimuth, overriding the current azimuth
                           (C{compass degrees360}).
        '''
        a = antipode(self.lat, self.lon)  # PYCHOK named
        z = self.azimuth if azimuth is None else Bearing(azimuth=azimuth)  # PYCHOK named
        return _NamedTuple.dup(self, lat=a.lat, lon=a.lon, azimuth=z)


class Equidistant(_AzimuthalBase):
    '''Azimuthal equidistant projection for the sphere***, see U{Snyder, pp 195-197
       <https://Pubs.USGS.gov/pp/1395/report.pdf>} and U{MathWorld-Wolfram
       <https://MathWorld.Wolfram.com/AzimuthalEquidistantProjection.html>}.

       @note: Results from this L{Equidistant} and an L{EquidistantExact},
              L{EquidistantGeodSolve} or L{EquidistantKarney} projection
              C{may differ} by 10% or more.  For an example, see method
              C{testDiscrepancies} in module C{testAzimuthal.py}.
    '''
    if _FOR_DOCS:
        __init__ = _AzimuthalBase.__init__

    def forward(self, lat, lon, name=NN):
        '''Convert a geodetic location to azimuthal equidistant east- and northing.

           @arg lat: Latitude of the location (C{degrees90}).
           @arg lon: Longitude of the location (C{degrees180}).
           @kwarg name: Optional name for the location (C{str}).

           @return: An L{Azimuthal7Tuple}C{(x, y, lat, lon, azimuth, scale, datum)}
                    with easting C{x} and northing C{y} of point in C{meter} and C{lat}
                    and C{lon} in C{degrees} and C{azimuth} clockwise from true North.
                    The C{scale} of the projection is C{1} in I{radial} direction and
                    is C{1 / reciprocal} in the direction perpendicular to this.

           @raise AzimuthalError: Invalid B{C{lat}} or B{C{lon}}.

           @note: The C{scale} will be C{-1} if B{C{(lat, lon)}} is antipodal to
                  the projection center C{(lat0, lon0)}.
        '''
        def _k_t(c):
            k = _N_1_0 if c < 0 else _1_0
            t =  abs(c) < EPS1
            if t:
                a = acos(c)
                s = sin(a)
                if s:
                    k = a / s
            return k, t

        return self._forward(lat, lon, name, _k_t)

    def reverse(self, x, y, name=NN, **LatLon_and_kwds):
        '''Convert an azimuthal equidistant location to geodetic lat- and longitude.

           @arg x: Easting of the location (C{meter}).
           @arg y: Northing of the location (C{meter}).
           @kwarg name: Optional name for the location (C{str}).
           @kwarg LatLon_and_kwds: Optional, C{B{LatLon}=None} class to use and
                                   additional B{C{LatLon}} keyword arguments,
                                   ignored if C{B{LatLon} is None} or not given.

           @return: The geodetic (C{LatLon}) or if B{C{LatLon}} is C{None} an
                    L{Azimuthal7Tuple}C{(x, y, lat, lon, azimuth, scale, datum)}.

           @note: The C{lat} will be in the range C{[-90..90] degrees} and C{lon}
                  in the range C{[-180..180] degrees}.  The C{scale} of the
                  projection is C{1} in I{radial} direction, C{azimuth} clockwise
                  from true North and is C{1 / reciprocal} in the direction
                  perpendicular to this.
        '''
        def _c(c):
            return c if c > EPS else None

        return self._reverse(x, y, name, _c, False, **LatLon_and_kwds)


def equidistant(lat0, lon0, datum=_WGS84, exact=False, geodsolve=False, name=NN):
    '''Return an L{EquidistantExact}, L{EquidistantGeodSolve} or (if I{Karney}'s
       U{geographiclib<https://PyPI.org/project/geographiclib>} package is
       installed) an L{EquidistantKarney}, otherwise an L{Equidistant} instance.

       @arg lat0: Latitude of center point (C{degrees90}).
       @arg lon0: Longitude of center point (C{degrees180}).
       @kwarg datum: Optional datum or ellipsoid (L{Datum}, L{Ellipsoid},
                     L{Ellipsoid2} or L{a_f2Tuple}) or I{scalar} earth
                     radius (C{meter}).
       @kwarg exact: Return an L{EquidistantExact} instance.
       @kwarg geodsolve: Return an L{EquidistantGeodSolve} instance.
       @kwarg name: Optional name for the projection (C{str}).

       @return: An L{EquidistantExact}, L{EquidistantGeodSolve},
                L{EquidistantKarney} or L{Equidistant} instance.

       @raise AzimuthalError: Invalid B{C{lat0}} or B{C{lon0}} or (spherical) B{C{datum}}.

       @raise GeodesicError: Issue with L{GeodesicExact}, L{GeodesicSolve}
                             or I{Karney}'s wrapped C{Geodesic}.

       @raise TypeError: Invalid B{C{datum}}.
    '''

    E = EquidistantExact if exact else (EquidistantGeodSolve if geodsolve else Equidistant)
    if E is Equidistant:
        try:
            return EquidistantKarney(lat0, lon0, datum=datum, name=name)  # PYCHOK types
        except ImportError:
            pass
    return E(lat0, lon0, datum=datum, name=name)  # PYCHOK types


class _AzimuthalGeodesic(_AzimuthalBase):
    '''(INTERNAL) Base class for azimuthal projections using U{Karney Geodesic
       <https://GeographicLib.SourceForge.io/C++/doc/python/code.html>} or
       exact geodesic classes L{GeodesicExact} and L{GeodesicLineExact}.
    '''
    _mask = 0

    @Property_RO
    def geodesic(self):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}, see function C{notOverloaded}.
        '''
        notOverloaded(self)

    def _7Tuple(self, e, n, r, M=None, name=NN):
        '''(INTERNAL) Return an C{Azimuthal7Tuple}.
        '''
        s = M if M is not None else (  # reciprocal, azimuthal scale
            (r.m12 / r.s12) if r.a12 > _EPS_K else _1_0)
        z = _umod_360(r.azi2)  # -180 <= r.azi2 < 180 ... 0 <= z < 360
        return Azimuthal7Tuple(e, n, r.lat2, r.lon2, z, s, self.datum,
                                     name=name or self.name)


class _EquidistantBase(_AzimuthalGeodesic):
    '''(INTERNAL) Base for classes L{EquidistantExact}, L{EquidistantGeodSolve}
       and L{EquidistantKarney}.
    '''
    def __init__(self, lat0, lon0, datum=_WGS84, name=NN):
        '''New azimuthal L{EquidistantExact}, L{EquidistantGeodSolve} or
           L{EquidistantKarney} projection.
        '''
        _AzimuthalGeodesic.__init__(self, lat0, lon0, datum=datum, name=name)

        g = self.geodesic
        # g.STANDARD = g.AZIMUTH | g.DISTANCE | g.LATITUDE | g.LONGITUDE
        self._mask = g.REDUCEDLENGTH | g.STANDARD  # | g.LONG_UNROLL

    def forward(self, lat, lon, name=NN):
        '''Convert an (ellipsoidal) geodetic location to azimuthal equidistant east- and northing.

           @arg lat: Latitude of the location (C{degrees90}).
           @arg lon: Longitude of the location (C{degrees180}).
           @kwarg name: Optional name for the location (C{str}).

           @return: An L{Azimuthal7Tuple}C{(x, y, lat, lon, azimuth, scale, datum)}
                    with easting C{x} and northing C{y} of point in C{meter} and C{lat}
                    and C{lon} in C{degrees} and C{azimuth} clockwise from true North.
                    The C{scale} of the projection is C{1} in I{radial} direction and
                    is C{1 / reciprocal} in the direction perpendicular to this.

           @raise AzimuthalError: Invalid B{C{lat}} or B{C{lon}}.

           @note: A call to C{.forward} followed by a call to C{.reverse} will return
                  the original C{lat, lon} to within roundoff.
       '''
        r = self.geodesic.Inverse(self.lat0, self.lon0,
                                  Lat_(lat), Lon_(lon), self._mask)
        x, y = sincos2d(r.azi1)
        return self._7Tuple(x * r.s12, y * r.s12, r, name=name)

    def reverse(self, x, y, name=NN, LatLon=None, **LatLon_kwds):  # PYCHOK signature
        '''Convert an azimuthal equidistant location to (ellipsoidal) geodetic lat- and longitude.

           @arg x: Easting of the location (C{meter}).
           @arg y: Northing of the location (C{meter}).
           @kwarg name: Optional name for the location (C{str}).
           @kwarg LatLon: Class to use (C{LatLon}) or C{None}.
           @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                               arguments, ignored if C{B{LatLon} is None}.

           @return: The geodetic (C{LatLon}) or if B{C{LatLon}} is C{None} an
                    L{Azimuthal7Tuple}C{(x, y, lat, lon, azimuth, scale, datum)}.

           @note: The C{lat} will be in the range C{[-90..90] degrees} and C{lon}
                  in the range C{[-180..180] degrees}.  The scale of the projection
                  is C{1} in I{radial} direction, C{azimuth} clockwise from true
                  North and is C{1 / reciprocal} in the direction perpendicular
                  to this.
        '''
        e, n, z, s = _enzh4(x, y)

        r = self.geodesic.Direct(self.lat0, self.lon0, z, s, self._mask)
        return self._7Tuple(e, n, r, name=name) if LatLon is None else \
               self._toLatLon(r.lat2, r.lon2, LatLon, LatLon_kwds, name)


class EquidistantExact(_EquidistantBase):
    '''Azimuthal equidistant projection, a Python version of I{Karney}'s C++ class U{AzimuthalEquidistant
       <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1AzimuthalEquidistant.html>},
       based on exact geodesic classes L{GeodesicExact} and L{GeodesicLineExact}.

       An azimuthal equidistant projection is centered at an arbitrary position on the ellipsoid.
       For a point in projected space C{(x, y)}, the geodesic distance from the center position
       is C{hypot(x, y)} and the C{azimuth} of the geodesic from the center point is C{atan2(x, y)},
       clockwise from true North.

       The C{.forward} and C{.reverse} methods also return the C{azimuth} of the geodesic at C{(x,
       y)} and the C{scale} in the azimuthal direction which, together with the basic properties
       of the projection, serve to specify completely the local affine transformation between
       geographic and projected coordinates.
    '''
    def __init__(self, lat0, lon0, datum=_WGS84, name=NN):
        '''New azimuthal L{EquidistantExact} projection.

           @arg lat0: Latitude of center point (C{degrees90}).
           @arg lon0: Longitude of center point (C{degrees180}).
           @kwarg datum: Optional datum or ellipsoid (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}) or I{scalar} earth
                         radius (C{meter}).
           @kwarg name: Optional name for the projection (C{str}).

           @raise AzimuthalError: Invalid B{C{lat0}} or B{C{lon0}} or B{C{datum}}.
        '''
        _EquidistantBase.__init__(self, lat0, lon0, datum=datum, name=name)

    if _FOR_DOCS:
        forward = _EquidistantBase.forward
        reverse = _EquidistantBase.reverse

    @Property_RO
    def geodesic(self):
        '''Get this projection's exact geodesic (L{GeodesicExact}).
        '''
        return self.datum.ellipsoid.geodesicx


class EquidistantGeodSolve(_EquidistantBase):
    '''Azimuthal equidistant projection, a Python version of I{Karney}'s C++ class U{AzimuthalEquidistant
       <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1AzimuthalEquidistant.html>},
       based on (exact) geodesic I{wrappers} L{GeodesicSolve} and L{GeodesicLineSolve} and intended
       I{for testing purposes only}.

       @see: L{EquidistantExact} and module L{geodsolve}.
    '''
    def __init__(self, lat0, lon0, datum=_WGS84, name=NN):
        '''New azimuthal L{EquidistantGeodSolve} projection.

           @arg lat0: Latitude of center point (C{degrees90}).
           @arg lon0: Longitude of center point (C{degrees180}).
           @kwarg datum: Optional datum or ellipsoid (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}) or I{scalar} earth
                         radius (C{meter}).
           @kwarg name: Optional name for the projection (C{str}).

           @raise AzimuthalError: Invalid B{C{lat0}} or B{C{lon0}} or B{C{datum}}.
        '''
        _EquidistantBase.__init__(self, lat0, lon0, datum=datum, name=name)

    if _FOR_DOCS:
        forward = _EquidistantBase.forward
        reverse = _EquidistantBase.reverse

    @Property_RO
    def geodesic(self):
        '''Get this projection's (exact) geodesic (L{GeodesicSolve}).
        '''
        return self.datum.ellipsoid.geodsolve


class EquidistantKarney(_EquidistantBase):
    '''Azimuthal equidistant projection, a Python version of I{Karney}'s C++ class U{AzimuthalEquidistant
       <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1AzimuthalEquidistant.html>},
       requiring package U{geographiclib<https://PyPI.org/project/geographiclib>} to be installed.

       @see: L{EquidistantExact}.
    '''
    def __init__(self, lat0, lon0, datum=_WGS84, name=NN):
        '''New azimuthal L{EquidistantKarney} projection.

           @arg lat0: Latitude of center point (C{degrees90}).
           @arg lon0: Longitude of center point (C{degrees180}).
           @kwarg datum: Optional datum or ellipsoid (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}) or I{scalar} earth
                         radius (C{meter}).
           @kwarg name: Optional name for the projection (C{str}).

           @raise ImportError: Package U{geographiclib<https://PyPI.org/project/geographiclib>}
                               not installed or not found.

           @raise AzimuthalError: Invalid B{C{lat0}} or B{C{lon0}} or B{C{datum}}.
        '''
        _EquidistantBase.__init__(self, lat0, lon0, datum=datum, name=name)

    if _FOR_DOCS:
        forward = _EquidistantBase.forward
        reverse = _EquidistantBase.reverse

    @Property_RO
    def geodesic(self):
        '''Get this projection's I{wrapped} U{Karney Geodesic
           <https://GeographicLib.SourceForge.io/C++/doc/python/code.html>},
           provided package U{geographiclib
           <https://PyPI.org/project/geographiclib>} is installed.
        '''
        return self.datum.ellipsoid.geodesic


_Equidistants = (Equidistant, EquidistantExact, EquidistantGeodSolve,
                 EquidistantKarney)  # PYCHOK in .ellipsoidalBaseDI


class Gnomonic(_AzimuthalBase):
    '''Azimuthal gnomonic projection for the sphere***, see U{Snyder, pp 164-168
       <https://Pubs.USGS.gov/pp/1395/report.pdf>} and U{MathWorld-Wolfram
       <https://MathWorld.Wolfram.com/GnomonicProjection.html>}.
    '''
    if _FOR_DOCS:
        __init__ = _AzimuthalBase.__init__

    def forward(self, lat, lon, name=NN):
        '''Convert a geodetic location to azimuthal equidistant east- and northing.

           @arg lat: Latitude of the location (C{degrees90}).
           @arg lon: Longitude of the location (C{degrees180}).
           @kwarg name: Optional name for the location (C{str}).

           @return: An L{Azimuthal7Tuple}C{(x, y, lat, lon, azimuth, scale, datum)}
                    with easting C{x} and northing C{y} of point in C{meter} and C{lat}
                    and C{lon} in C{degrees} and C{azimuth} clockwise from true North.
                    The C{scale} of the projection is C{1} in I{radial} direction and
                    is C{1 / reciprocal} in the direction perpendicular to this.

           @raise AzimuthalError: Invalid B{C{lat}} or B{C{lon}}.
        '''
        def _k_t(c):
            t = c > EPS
            k = (_1_0 / c) if t else _1_0
            return k, t

        return self._forward(lat, lon, name, _k_t)

    def reverse(self, x, y, name=NN, **LatLon_and_kwds):
        '''Convert an azimuthal equidistant location to geodetic lat- and longitude.

           @arg x: Easting of the location (C{meter}).
           @arg y: Northing of the location (C{meter}).
           @kwarg name: Optional name for the location (C{str}).
           @kwarg LatLon_and_kwds: Optional, C{B{LatLon}=None} class to use and
                                   additional B{C{LatLon}} keyword arguments,
                                   ignored if C{B{LatLon} is None} or not given.

           @return: The geodetic (C{LatLon}) or if B{C{LatLon}} is C{None} an
                    L{Azimuthal7Tuple}C{(x, y, lat, lon, azimuth, scale, datum)}.

           @note: The C{lat} will be in the range C{[-90..90] degrees} and C{lon}
                  in the range C{[-180..180] degrees}.  The C{scale} of the
                  projection is C{1} in I{radial} direction, C{azimuth} clockwise
                  from true North and C{1 / reciprocal} in the direction
                  perpendicular to this.
        '''
        def _c(c):
            return atan(c) if c > EPS else None

        return self._reverse(x, y, name, _c, False, **LatLon_and_kwds)


def gnomonic(lat0, lon0, datum=_WGS84, exact=False, geodsolve=False, name=NN):
    '''Return a L{GnomonicExact} or (if I{Karney}'s U{geographiclib
       <https://PyPI.org/project/geographiclib>} package is installed)
       a L{GnomonicKarney}, otherwise a L{Gnomonic} instance.

       @arg lat0: Latitude of center point (C{degrees90}).
       @arg lon0: Longitude of center point (C{degrees180}).
       @kwarg datum: Optional datum or ellipsoid (L{Datum}, L{Ellipsoid},
                     L{Ellipsoid2} or L{a_f2Tuple}) or I{scalar} earth
                     radius (C{meter}).
       @kwarg exact: Return a L{GnomonicExact} instance.
       @kwarg geodsolve: Return an L{GnomonicGeodSolve} instance.
       @kwarg name: Optional name for the projection (C{str}).

       @return: A L{GnomonicExact}, L{GnomonicGeodSolve},
                L{GnomonicKarney} or L{Gnomonic} instance.

       @raise AzimuthalError: Invalid B{C{lat0}} or B{C{lon0}} or
                              (spherical) B{C{datum}}.

       @raise GeodesicError: Issue with L{GeodesicExact}, L{GeodesicSolve}
                             or I{Karney}'s wrapped C{Geodesic}.

       @raise TypeError: Invalid B{C{datum}}.
    '''
    G = GnomonicExact if exact else (GnomonicGeodSolve if geodsolve else Gnomonic)
    if G is Gnomonic:
        try:
            return GnomonicKarney(lat0, lon0, datum=datum, name=name)  # PYCHOK types
        except ImportError:
            pass
    return G(lat0, lon0, datum=datum, name=name)  # PYCHOK types


class _GnomonicBase(_AzimuthalGeodesic):
    '''(INTERNAL) Base for classes L{GnomonicExact}, L{GnomonicGeodSolve}
       and L{GnomonicKarney}.
    '''
    def __init__(self, lat0, lon0, datum=_WGS84, name=NN):
        '''New azimuthal L{GnomonicExact} or L{GnomonicKarney} projection.
        '''
        _AzimuthalGeodesic.__init__(self, lat0, lon0, datum=datum, name=name)

        g = self.geodesic
        self._mask = g.ALL  # | g.LONG_UNROLL

    def forward(self, lat, lon, name=NN, raiser=True):  # PYCHOK signature
        '''Convert an (ellipsoidal) geodetic location to azimuthal gnomonic east-
           and northing.

           @arg lat: Latitude of the location (C{degrees90}).
           @arg lon: Longitude of the location (C{degrees180}).
           @kwarg name: Optional name for the location (C{str}).
           @kwarg raiser: Do or don't throw an error (C{bool}) if
                          the location lies over the horizon.

           @return: An L{Azimuthal7Tuple}C{(x, y, lat, lon, azimuth, scale, datum)}
                    with easting C{x} and northing C{y} in C{meter} and C{lat} and
                    C{lon} in C{degrees} and C{azimuth} clockwise from true North.
                    The C{scale} of the projection is C{1 / reciprocal**2} in I{radial}
                    direction and C{1 / reciprocal} in the direction perpendicular to
                    this.  Both C{x} and C{y} will be C{NAN} if the (geodetic) location
                    lies over the horizon and C{B{raiser}=False}.

           @raise AzimuthalError: Invalid B{C{lat}}, B{C{lon}} or the location lies
                                  over the horizon and C{B{raiser}=True}.
        '''
        self._iteration = 0

        r = self.geodesic.Inverse(self.lat0, self.lon0,
                                  Lat_(lat), Lon_(lon), self._mask)
        M = r.M21
        if M > EPS0:
            q = r.m12 / M  # .M12
            e, n = sincos2d(r.azi1)
            e *= q
            n *= q
        elif raiser:  # PYCHOK no cover
            raise AzimuthalError(lat=lat, lon=lon, txt=_over_horizon_)
        else:  # PYCHOK no cover
            e = n = NAN

        t = self._7Tuple(e, n, r, M=M, name=name)
        t._iteraton = 0
        return t

    def reverse(self, x, y, name=NN, LatLon=None, **LatLon_kwds):  # PYCHOK signature
        '''Convert an azimuthal gnomonic location to (ellipsoidal) geodetic lat- and longitude.

           @arg x: Easting of the location (C{meter}).
           @arg y: Northing of the location (C{meter}).
           @kwarg name: Optional name for the location (C{str}).
           @kwarg LatLon: Class to use (C{LatLon}) or C{None}.
           @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                               arguments, ignored if C{B{LatLon} is None}.

           @return: The geodetic (C{LatLon}) or if B{C{LatLon}} is C{None} an
                    L{Azimuthal7Tuple}C{(x, y, lat, lon, azimuth, scale, datum)}.

           @raise AzimuthalError: No convergence.

           @note: The C{lat} will be in the range C{[-90..90] degrees} and C{lon}
                  in the range C{[-180..180] degrees}.  The C{azimuth} is clockwise
                  from true North.  The scale is C{1 / reciprocal**2} in C{radial}
                  direction and C{1 / reciprocal} in the direction perpendicular
                  to this.
        '''
        e, n, z, q = _enzh4(x, y)

        d = a = self.equatoradius
        s = a * atan(q / a)
        if q > a:  # PYCHOK no cover
            def _d(r, q):
                return (r.M12 - q * r.m12) * r.m12  # negated

            q = _1_0 / q
        else:  # little == True
            def _d(r, q):  # PYCHOK _d
                return (q * r.M12 - r.m12) * r.M12  # negated

        a *= _EPS_K
        m  =  self._mask
        g  =  self.geodesic

        P   = g.Line(self.lat0, self.lon0, z, m | g.LINE_OFF).Position
        S2_ = Fsum(s).fsum2_
        for i in range(1, _TRIPS):
            r = P(s, m)
            if abs(d) < a:
                break
            s, d = S2_(_d(r, q))
        else:  # PYCHOK no cover
            self._iteration = _TRIPS
            raise AzimuthalError(x=x, y=y, txt=_no_(Fmt.convergence(a)))

        t = self._7Tuple(e, n, r, M=r.M12, name=name) if LatLon is None else \
            self._toLatLon(r.lat2, r.lon2, LatLon, LatLon_kwds, name)
        t._iteration = self._iteration = i
        return t


class GnomonicExact(_GnomonicBase):
    '''Azimuthal gnomonic projection, a Python version of I{Karney}'s C++ class U{Gnomonic
       <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Gnomonic.html>},
       based on exact geodesic classes L{GeodesicExact} and L{GeodesicLineExact}.

       @see: I{Karney}'s U{Detailed Description<https://GeographicLib.SourceForge.io/C++/doc/
             classGeographicLib_1_1Gnomonic.html>}, especially the B{Warning}.
    '''
    def __init__(self, lat0, lon0, datum=_WGS84, name=NN):
        '''New azimuthal L{GnomonicExact} projection.

           @arg lat0: Latitude of center point (C{degrees90}).
           @arg lon0: Longitude of center point (C{degrees180}).
           @kwarg datum: Optional datum or ellipsoid (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}) or I{scalar} earth
                         radius (C{meter}).
           @kwarg name: Optional name for the projection (C{str}).

           @raise AzimuthalError: Invalid B{C{lat0}} or B{C{lon0}}.
        '''
        _GnomonicBase.__init__(self, lat0, lon0, datum=datum, name=name)

    if _FOR_DOCS:
        forward = _GnomonicBase.forward
        reverse = _GnomonicBase.reverse

    @Property_RO
    def geodesic(self):
        '''Get this projection's exact geodesic (L{GeodesicExact}).
        '''
        return self.datum.ellipsoid.geodesicx


class GnomonicGeodSolve(_GnomonicBase):
    '''Azimuthal gnomonic projection, a Python version of I{Karney}'s C++ class U{Gnomonic
       <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Gnomonic.html>},
       based on (exact) geodesic I{wrappers} L{GeodesicSolve} and L{GeodesicLineSolve} and
       intended I{for testing purposes only}.

       @see: L{GnomonicExact} and module L{geodsolve}.
    '''
    def __init__(self, lat0, lon0, datum=_WGS84, name=NN):
        '''New azimuthal L{GnomonicGeodSolve} projection.

           @arg lat0: Latitude of center point (C{degrees90}).
           @arg lon0: Longitude of center point (C{degrees180}).
           @kwarg datum: Optional datum or ellipsoid (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}) or I{scalar} earth
                         radius (C{meter}).
           @kwarg name: Optional name for the projection (C{str}).

           @raise AzimuthalError: Invalid B{C{lat0}} or B{C{lon0}}.
        '''
        _GnomonicBase.__init__(self, lat0, lon0, datum=datum, name=name)

    if _FOR_DOCS:
        forward = _GnomonicBase.forward
        reverse = _GnomonicBase.reverse

    @Property_RO
    def geodesic(self):
        '''Get this projection's (exact) geodesic (L{GeodesicSolve}).
        '''
        return self.datum.ellipsoid.geodsolve


class GnomonicKarney(_GnomonicBase):
    '''Azimuthal gnomonic projection, a Python version of I{Karney}'s C++ class U{Gnomonic
       <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Gnomonic.html>},
       requiring package U{geographiclib<https://PyPI.org/project/geographiclib>} to be installed.

       @see: L{GnomonicExact}.
    '''
    def __init__(self, lat0, lon0, datum=_WGS84, name=NN):
        '''New azimuthal L{GnomonicKarney} projection.

           @arg lat0: Latitude of center point (C{degrees90}).
           @arg lon0: Longitude of center point (C{degrees180}).
           @kwarg datum: Optional datum or ellipsoid (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}) or I{scalar} earth
                         radius (C{meter}).
           @kwarg name: Optional name for the projection (C{str}).

           @raise ImportError: Package U{geographiclib<https://PyPI.org/project/geographiclib>}
                               not installed or not found.

           @raise AzimuthalError: Invalid B{C{lat0}} or B{C{lon0}}.
        '''
        _GnomonicBase.__init__(self, lat0, lon0, datum=datum, name=name)

    if _FOR_DOCS:
        forward = _GnomonicBase.forward
        reverse = _GnomonicBase.reverse

    @Property_RO
    def geodesic(self):
        '''Get this projection's I{wrapped} U{Karney Geodesic
           <https://GeographicLib.SourceForge.io/C++/doc/python/code.html>}, provided package
           U{geographiclib<https://PyPI.org/project/geographiclib>} is installed.
        '''
        return self.datum.ellipsoid.geodesic


class LambertEqualArea(_AzimuthalBase):
    '''Lambert-equal-area projection for the sphere*** (aka U{Lambert zenithal equal-area
       projection<https://WikiPedia.org/wiki/Lambert_azimuthal_equal-area_projection>}, see
       U{Snyder, pp 185-187<https://Pubs.USGS.gov/pp/1395/report.pdf>} and U{MathWorld-Wolfram
       <https://MathWorld.Wolfram.com/LambertAzimuthalEqual-AreaProjection.html>}.
    '''
    if _FOR_DOCS:
        __init__ = _AzimuthalBase.__init__

    def forward(self, lat, lon, name=NN):
        '''Convert a geodetic location to azimuthal Lambert-equal-area east- and northing.

           @arg lat: Latitude of the location (C{degrees90}).
           @arg lon: Longitude of the location (C{degrees180}).
           @kwarg name: Optional name for the location (C{str}).

           @return: An L{Azimuthal7Tuple}C{(x, y, lat, lon, azimuth, scale, datum)}
                    with easting C{x} and northing C{y} of point in C{meter} and C{lat}
                    and C{lon} in C{degrees} and C{azimuth} clockwise from true North.
                    The C{scale} of the projection is C{1} in I{radial} direction and
                    is C{1 / reciprocal} in the direction perpendicular to this.

           @raise AzimuthalError: Invalid B{C{lat}} or B{C{lon}}.
        '''
        def _k_t(c):
            c += _1_0
            t  =  c > EPS0
            k  =  sqrt(_2_0 / c) if t else _1_0
            return k, t

        return self._forward(lat, lon, name, _k_t)

    def reverse(self, x, y, name=NN, **LatLon_and_kwds):
        '''Convert an azimuthal Lambert-equal-area location to geodetic lat- and longitude.

           @arg x: Easting of the location (C{meter}).
           @arg y: Northing of the location (C{meter}).
           @kwarg name: Optional name for the location (C{str}).
           @kwarg LatLon_and_kwds: Optional, C{B{LatLon}=None} class to use and
                                   additional B{C{LatLon}} keyword arguments,
                                   ignored if C{B{LatLon} is None} or not given.

           @return: The geodetic (C{LatLon}) or if B{C{LatLon}} is C{None} an
                    L{Azimuthal7Tuple}C{(x, y, lat, lon, azimuth, scale, datum)}.

           @note: The C{lat} will be in the range C{[-90..90] degrees} and C{lon}
                  in the range C{[-180..180] degrees}.  The C{scale} of the
                  projection is C{1} in I{radial} direction, C{azimuth} clockwise
                  from true North and is C{1 / reciprocal} in the direction
                  perpendicular to this.
        '''
        def _c(c):
            c *= _0_5
            return (asin1(c) * _2_0) if c > EPS else None

        return self._reverse(x, y, name, _c, True, **LatLon_and_kwds)


class Orthographic(_AzimuthalBase):
    '''Orthographic projection for the sphere***, see U{Snyder, pp 148-153
       <https://Pubs.USGS.gov/pp/1395/report.pdf>} and U{MathWorld-Wolfram
       <https://MathWorld.Wolfram.com/OrthographicProjection.html>}.
    '''
    if _FOR_DOCS:
        __init__ = _AzimuthalBase.__init__

    def forward(self, lat, lon, name=NN):
        '''Convert a geodetic location to azimuthal orthographic east- and northing.

           @arg lat: Latitude of the location (C{degrees90}).
           @arg lon: Longitude of the location (C{degrees180}).
           @kwarg name: Optional name for the location (C{str}).

           @return: An L{Azimuthal7Tuple}C{(x, y, lat, lon, azimuth, scale, datum)}
                    with easting C{x} and northing C{y} of point in C{meter} and C{lat}
                    and C{lon} in C{degrees} and C{azimuth} clockwise from true North.
                    The C{scale} of the projection is C{1} in I{radial} direction and
                    is C{1 / reciprocal} in the direction perpendicular to this.

           @raise AzimuthalError: Invalid B{C{lat}} or B{C{lon}}.
        '''
        def _k_t(c):
            return _1_0, (c >= 0)

        return self._forward(lat, lon, name, _k_t)

    def reverse(self, x, y, name=NN, **LatLon_and_kwds):
        '''Convert an azimuthal orthographic location to geodetic lat- and longitude.

           @arg x: Easting of the location (C{meter}).
           @arg y: Northing of the location (C{meter}).
           @kwarg name: Optional name for the location (C{str}).
           @kwarg LatLon_and_kwds: Optional, C{B{LatLon}=None} class to use and
                                   additional B{C{LatLon}} keyword arguments,
                                   ignored if C{B{LatLon} is None} or not given.

           @return: The geodetic (C{LatLon}) or if B{C{LatLon}} is C{None} an
                    L{Azimuthal7Tuple}C{(x, y, lat, lon, azimuth, scale, datum)}.

           @note: The C{lat} will be in the range C{[-90..90] degrees} and C{lon}
                  in the range C{[-180..180] degrees}.  The C{scale} of the
                  projection is C{1} in I{radial} direction, C{azimuth} clockwise
                  from true North and is C{1 / reciprocal} in the direction
                  perpendicular to this.
        '''
        def _c(c):
            return asin1(c) if c > EPS else None

        return self._reverse(x, y, name, _c, False, **LatLon_and_kwds)


class Stereographic(_AzimuthalBase):
    '''Stereographic projection for the sphere***, see U{Snyder, pp 157-160
       <https://Pubs.USGS.gov/pp/1395/report.pdf>} and U{MathWorld-Wolfram
       <https://MathWorld.Wolfram.com/StereographicProjection.html>}.
    '''
    _k0  = _1_0  # central scale factor (C{scalar})
    _k02 = _2_0  # double ._k0

    if _FOR_DOCS:
        __init__ = _AzimuthalBase.__init__

    def forward(self, lat, lon, name=NN):
        '''Convert a geodetic location to azimuthal stereographic east- and northing.

           @arg lat: Latitude of the location (C{degrees90}).
           @arg lon: Longitude of the location (C{degrees180}).
           @kwarg name: Optional name for the location (C{str}).

           @return: An L{Azimuthal7Tuple}C{(x, y, lat, lon, azimuth, scale, datum)}
                    with easting C{x} and northing C{y} of point in C{meter} and C{lat}
                    and C{lon} in C{degrees} and C{azimuth} clockwise from true North.
                    The C{scale} of the projection is C{1} in I{radial} direction and
                    is C{1 / reciprocal} in the direction perpendicular to this.

           @raise AzimuthalError: Invalid B{C{lat}} or B{C{lon}}.
        '''
        def _k_t(c):
            c += _1_0
            t  =  isnon0(c)
            k  = (self._k02 / c) if t else _1_0
            return k, t

        return self._forward(lat, lon, name, _k_t)

    @property_doc_(''' optional, central scale factor (C{scalar}).''')
    def k0(self):
        '''Get the central scale factor (C{scalar}).
        '''
        return self._k0

    @k0.setter  # PYCHOK setter!
    def k0(self, factor):
        '''Set the central scale factor (C{scalar}).
        '''
        n = Stereographic.k0.fget.__name__
        self._k0  = Scalar_(factor, name=n, low=EPS, high=2)  # XXX high=1, 2, other?
        self._k02 = self._k0 * _2_0

    def reverse(self, x, y, name=NN, **LatLon_and_kwds):
        '''Convert an azimuthal stereographic location to geodetic lat- and longitude.

           @arg x: Easting of the location (C{meter}).
           @arg y: Northing of the location (C{meter}).
           @kwarg name: Optional name for the location (C{str}).
           @kwarg LatLon_and_kwds: Optional, C{B{LatLon}=None} class to use and
                                   additional B{C{LatLon}} keyword arguments,
                                   ignored if C{B{LatLon} is None} or not given.

           @return: The geodetic (C{LatLon}) or if B{C{LatLon}} is C{None} an
                    L{Azimuthal7Tuple}C{(x, y, lat, lon, azimuth, scale, datum)}.

           @note: The C{lat} will be in the range C{[-90..90] degrees} and C{lon}
                  in the range C{[-180..180] degrees}.  The C{scale} of the
                  projection is C{1} in I{radial} direction, C{azimuth} clockwise
                  from true North and is C{1 / reciprocal} in the direction
                  perpendicular to this.
        '''
        def _c(c):
            return (_2_0 * atan2(c, self._k02)) if c > EPS else None

        return self._reverse(x, y, name, _c, False, **LatLon_and_kwds)


__all__ += _ALL_DOCS(_AzimuthalBase, _AzimuthalGeodesic, _EquidistantBase, _GnomonicBase)

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
