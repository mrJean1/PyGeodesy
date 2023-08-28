
# -*- coding: utf-8 -*-

u'''Web Mercator (WM) projection.

Classes L{Wm} and L{WebMercatorError} and functions L{parseWM} and L{toWm}.

Pure Python implementation of a U{Web Mercator<https://WikiPedia.org/wiki/Web_Mercator>}
(aka I{Pseudo-Mercator}) class and conversion functions for spherical and near-spherical
earth models.

@see: U{Google Maps / Bing Maps Spherical Mercator Projection
<https://AlastairA.WordPress.com/2011/01/23/the-google-maps-bing-maps-spherical-mercator-projection>},
U{Geomatics Guidance Note 7, part 2<https://www.IOGP.org/wp-content/uploads/2019/09/373-07-02.pdf>}
and U{Implementation Practice Web Mercator Map Projection<https://Web.Archive.org/web/20141009142830/
http://earth-info.nga.mil/GandG/wgs84/web_mercator/(U)%20NGA_SIG_0011_1.0.0_WEBMERC.pdf>}.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import isscalar, _splituple, _xinstanceof
from pygeodesy.constants import PI_2, R_MA, _2_0
from pygeodesy.datums import Datum, _spherical_datum
from pygeodesy.dms import clipDegrees, parseDMS2
from pygeodesy.errors import _parseX, _ValueError, _xattr, _xkwds
from pygeodesy.interns import NN, _COMMASPACE_, _datum_, _earth_, _easting_, \
                             _northing_, _radius_, _SPACE_, _x_, _y_
# from pygeodesy.lazily import _ALL_LAZY  from .named
from pygeodesy.named import _NamedBase, _NamedTuple,  _ALL_LAZY
from pygeodesy.namedTuples import LatLon2Tuple, LatLonDatum3Tuple, PhiLam2Tuple
from pygeodesy.props import deprecated_method, Property_RO
from pygeodesy.streprs import Fmt, strs, _xzipairs
from pygeodesy.units import Easting, Lat, Northing, Radius
from pygeodesy.utily import degrees90, degrees180

from math import atan, atanh, exp, radians, sin, tanh

__all__ = _ALL_LAZY.webmercator
__version__ = '23.08.24'

# _FalseEasting  = 0   # false Easting (C{meter})
# _FalseNorthing = 0   # false Northing (C{meter})
_LatLimit = Lat(limit=85.051129)  # latitudinal limit (C{degrees})
# _LonOrigin     = 0   # longitude of natural origin (C{degrees})


class EasNorRadius3Tuple(_NamedTuple):
    '''3-Tuple C{(easting, northing, radius)}, all in C{meter}.
    '''
    _Names_ = (_easting_, _northing_, _radius_)
    _Units_ = ( Easting,   Northing,   Radius)


class WebMercatorError(_ValueError):
    '''Web Mercator (WM) parser or L{Wm} issue.
    '''
    pass


class Wm(_NamedBase):
    '''Web Mercator (WM) coordinate.
    '''
    _datum  = None  # set further below
    _earths = ()    # dito
    _radius = R_MA  # earth radius (C{meter})
    _x      = 0     # Easting (C{meter})
    _y      = 0     # Northing (C{meter})

    def __init__(self, x, y, earth=R_MA, name=NN, **radius):
        '''New L{Wm} Web Mercator (WM) coordinate.

           @arg x: Easting from central meridian (C{meter}).
           @arg y: Northing from equator (C{meter}).
           @kwarg earth: Earth radius (C{meter}), datum or
                         ellipsoid (L{Datum}, L{a_f2Tuple},
                         L{Ellipsoid} or L{Ellipsoid2}).
           @kwarg name: Optional name (C{str}).
           @kwarg radius: DEPRECATED, use keyword argument B{C{earth}}.

           @note: WM is strictly defined for spherical and WGS84
                  ellipsoidal earth models only.

           @raise WebMercatorError: Invalid B{C{x}}, B{C{y}} or B{C{radius}}.
        '''
        self._x = Easting( x=x, Error=WebMercatorError)
        self._y = Northing(y=y, Error=WebMercatorError)

        R = radius.get(_radius_, earth)
        if R not in Wm._earths:
            self._datum = _datum(R, _radius_ if radius else _earth_)
            self._radius = self.datum.ellipsoid.a

        if name:
            self.name = name

    @Property_RO
    def datum(self):
        '''Get the datum (C{Datum}).
        '''
        return self._datum

    @Property_RO
    def ellipsoid(self):
        '''Get the ellipsoid (C{Ellipsoid}).
        '''
        return self.datum.ellipsoid

    @Property_RO
    def latlon(self):
        '''Get the lat- and longitude (L{LatLon2Tuple}).
        '''
        return self.latlon2()

    def latlon2(self, datum=None):
        '''Convert this WM coordinate to a lat- and longitude.

           @kwarg datum: Optional datum (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}) or earth
                         radius (C{meter}), overriding this WM's
                         C{radius} and C{datum}.

           @return: A L{LatLon2Tuple}C{(lat, lon)}.

           @note: WM is strictly defined for spherical and WGS84
                  ellipsoidal earth models only.

           @raise TypeError: Invalid or non-ellipsoidal B{C{datum}}.

           @see: Method C{toLatLon} for other return types.
        '''
        d = self.datum if datum in (None, self.datum, self.radius) else _datum(datum)
        E = d.ellipsoid
        R = self.radius
        x = self.x / R
        y = atan(exp(self.y / R)) * _2_0 - PI_2
        if E.es or E.a != R:  # strictly, WGS84 only
            # <https://Web.Archive.org/web/20141009142830/http://earth-info.nga.mil/
            #        GandG/wgs84/web_mercator/(U)%20NGA_SIG_0011_1.0.0_WEBMERC.pdf>
            y  = y / R  # /= chokes PyChecker
            y -= E.es_atanh(tanh(y))
            y *= E.a
            x *= E.a / R

        return LatLon2Tuple(degrees90(y), degrees180(x), name=self.name)

    def parse(self, strWM, name=NN):
        '''Parse a string to a similar L{Wm} instance.

           @arg strWM: The WM coordinate (C{str}), see
                       function L{parseWM}.
           @kwarg name: Optional instance name (C{str}),
                        overriding this name.

           @return: The similar instance (L{Wm}).

           @raise WebMercatorError: Invalid B{C{strWM}}.
        '''
        return parseWM(strWM, radius=self.radius, Wm=self.classof,
                              name=name or self.name)

    @deprecated_method
    def parseWM(self, strWM, name=NN):  # PYCHOK no cover
        '''DEPRECATED, use method L{Wm.parse}.'''
        return self.parse(strWM, name=name)

    @Property_RO
    def philam(self):
        '''Get the lat- and longitude ((L{PhiLam2Tuple}).
        '''
        return PhiLam2Tuple(*map(radians, self.latlon), name=self.name)

    @Property_RO
    def radius(self):
        '''Get the earth radius (C{meter}).
        '''
        return self._radius

    @deprecated_method
    def to2ll(self, datum=None):  # PYCHOK no cover
        '''DEPRECATED, use method C{latlon2}.

           @return: A L{LatLon2Tuple}C{(lat, lon)}.
        '''
        return self.latlon2(datum=datum)

    def toLatLon(self, LatLon=None, datum=None, **LatLon_kwds):
        '''Convert this WM coordinate to a geodetic point.

           @kwarg LatLon: Ellipsoidal or sphperical C{LatLon} class to
                          return the geodetic point (C{LatLon}) or C{None}.
           @kwarg datum: Optional, datum (C{Datum}) overriding this WM's.
           @kwarg LatLon_kwds: Optional, additional B{C{LatLon}}
                               keyword arguments, ignored if
                               C{B{LatLon} is None}.

           @return: This WM coordinate as B{C{LatLon}} or if
                    C{B{LatLon} is None} a L{LatLonDatum3Tuple}.

           @raise TypeError: If B{C{LatLon}} and B{C{datum}} are
                             incompatible or if B{C{datum}} is
                             invalid.
        '''
        d = datum or self.datum
        _xinstanceof(Datum, datum=d)
        r = self.latlon2(datum=d)
        r = LatLonDatum3Tuple(r.lat, r.lon, d, name=r.name) if LatLon is None else \
                       LatLon(r.lat, r.lon, **_xkwds(LatLon_kwds, datum=d, name=r.name))
        return r

    def toRepr(self, prec=3, fmt=Fmt.SQUARE, sep=_COMMASPACE_, radius=False, **unused):  # PYCHOK expected
        '''Return a string representation of this WM coordinate.

           @kwarg prec: Number of (decimal) digits, unstripped (C{int}).
           @kwarg fmt: Enclosing backets format (C{str}).
           @kwarg sep: Optional separator between name:value pairs (C{str}).
           @kwarg radius: If C{True} include the radius (C{bool}) or
                          C{scalar} to override this WM's radius.

           @return: This WM as "[x:meter, y:meter]" (C{str}) or as "[x:meter,
                    y:meter], radius:meter]" if B{C{radius}} is C{True} or
                    C{scalar}.

           @raise WebMercatorError: Invalid B{C{radius}}.
        '''
        t = self.toStr(prec=prec, sep=None, radius=radius)
        n = (_x_, _y_, _radius_)[:len(t)]
        return _xzipairs(n, t, sep=sep, fmt=fmt)

    def toStr(self, prec=3, fmt=Fmt.F, sep=_SPACE_, radius=False, **unused):  # PYCHOK expected
        '''Return a string representation of this WM coordinate.

           @kwarg prec: Number of (decimal) digits, unstripped (C{int}).
           @kwarg fmt: The C{float} format (C{str}).
           @kwarg sep: Optional separator to join (C{str}) or C{None}
                       to return an unjoined C{tuple} of C{str}s.
           @kwarg radius: If C{True} include the radius (C{bool}) or
                          C{scalar} to override this WM's radius.

           @return: This WM as "meter meter" (C{str}) or as "meter meter
                    radius" if B{C{radius}} is C{True} or C{scalar}.

           @raise WebMercatorError: Invalid B{C{radius}}.
        '''
        fs = self.x, self.y
        if isscalar(radius):
            fs += (radius,)
        elif radius:  # is True:
            fs += (self.radius,)
        elif radius not in (None, False):
            raise WebMercatorError(radius=radius)
        t = strs(fs, prec=prec)
        return t if sep is None else sep.join(t)

    @Property_RO
    def x(self):
        '''Get the easting (C{meter}).
        '''
        return self._x

    @Property_RO
    def y(self):
        '''Get the northing (C{meter}).
        '''
        return self._y

Wm._datum  = _spherical_datum(Wm._radius, name=Wm.__name__, raiser=_radius_)  # PYCHOK defaults
Wm._earths = (Wm._radius, Wm._datum, Wm._datum.ellipsoid)


def _datum(earth, name=_datum_):
    '''(INTERNAL) Make a datum from an C{earth} radius, datum or ellipsoid.
    '''
    if earth in Wm._earths:
        return  Wm._datum
    try:
        return _spherical_datum(earth, name=name)
    except Exception as x:
        raise WebMercatorError(name, earth, cause=x)


def parseWM(strWM, radius=R_MA, Wm=Wm, name=NN):
    '''Parse a string C{"e n [r]"} representing a WM coordinate,
       consisting of easting, northing and an optional radius.

       @arg strWM: A WM coordinate (C{str}).
       @kwarg radius: Optional earth radius (C{meter}), needed in
                      case B{C{strWM}} doesn't include C{r}.
       @kwarg Wm: Optional class to return the WM coordinate (L{Wm})
                  or C{None}.
       @kwarg name: Optional name (C{str}).

       @return: The WM coordinate (B{C{Wm}}) or an
                L{EasNorRadius3Tuple}C{(easting, northing, radius)}
                if B{C{Wm}} is C{None}.

       @raise WebMercatorError: Invalid B{C{strWM}}.
    '''
    def _WM(strWM, radius, Wm, name):
        w = _splituple(strWM)

        if len(w) == 2:
            w += (radius,)
        elif len(w) != 3:
            raise ValueError
        x, y, R = map(float, w)

        return EasNorRadius3Tuple(x, y, R, name=name) if Wm is None else \
                               Wm(x, y, earth=R, name=name)

    return _parseX(_WM, strWM, radius, Wm, name,
                        strWM=strWM, Error=WebMercatorError)


def toWm(latlon, lon=None, earth=R_MA, Wm=Wm, name=NN, **Wm_kwds):
    '''Convert a lat-/longitude point to a WM coordinate.

       @arg latlon: Latitude (C{degrees}) or an (ellipsoidal or
                    spherical) geodetic C{LatLon} point.
       @kwarg lon: Optional longitude (C{degrees} or C{None}).
       @kwarg earth: Earth radius (C{meter}), datum or ellipsoid
                     (L{Datum}, L{a_f2Tuple}, L{Ellipsoid} or
                     L{Ellipsoid2}), overridden by B{C{latlon}}'s
                     datum if present.
       @kwarg Wm: Optional class to return the WM coordinate (L{Wm})
                  or C{None}.
       @kwarg name: Optional name (C{str}).
       @kwarg Wm_kwds: Optional, additional B{C{Wm}} keyword arguments,
                       ignored if C{B{Wm} is None}.

       @return: The WM coordinate (B{C{Wm}}) or if B{C{Wm}} is C{None}
                an L{EasNorRadius3Tuple}C{(easting, northing, radius)}.

       @raise ValueError: If B{C{lon}} value is missing, if B{C{latlon}} is not
                          scalar, if B{C{latlon}} is beyond the valid WM range
                          and L{pygeodesy.rangerrors} is set to C{True} or if
                          B{C{earth}} is invalid.
    '''
    if _radius_ in Wm_kwds:  # remove DEPRECATED, radius
        d = _datum(Wm_kwds.pop(_radius_), _radius_)
    else:
        d = _datum(earth, _earth_)
    try:
        y, x = latlon.lat, latlon.lon
        y =  clipDegrees(y, _LatLimit)
        d = _xattr(latlon, datum=d)
        n =  name or _xattr(latlon, name=NN)
    except AttributeError:
        y, x = parseDMS2(latlon, lon, clipLat=_LatLimit)
        n = name

    E = d.ellipsoid
    R = E.a
    s = sin(radians(y))
    y = atanh(s)  # == log(tand((90 + lat) / 2)) == log(tanPI_2_2(radians(lat)))
    if E.es:
        y -= E.es_atanh(s)  # strictly, WGS84 only
    y *= R
    x  = R * radians(x)
    r  = EasNorRadius3Tuple(x, y, R, name=n) if Wm is None else \
                         Wm(x, y, **_xkwds(Wm_kwds, earth=d, name=n))
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
