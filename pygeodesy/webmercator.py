
# -*- coding: utf-8 -*-

u'''Web Mercator (WM) classes L{Wm} and L{WebMercatorError} and functions
L{parseWM} and L{toWm}.

Pure Python implementation of a U{Web Mercator<https://WikiPedia.org/wiki/Web_Mercator>}
(aka I{Pseudo-Mercator}) class and conversion functions for spherical and
near-spherical earth models.

References U{Google Maps / Bing Maps Spherical Mercator Projection
<https://AlastairA.WordPress.com/2011/01/23/the-google-maps-bing-maps-spherical-mercator-projection>},
U{Geomatics Guidance Note 7, part 2<https://www.EPSG.org/Portals/0/373-07-02.pdf>} and
U{Implementation Practice Web Mercator Map Projection
<https://Earth-Info.NGA.mil/GandG/wgs84/web_mercator/%28U%29%20NGA_SIG_0011_1.0.0_WEBMERC.pdf>}.

@newfield example: Example, Examples
'''

from pygeodesy.basics import PI_2, isscalar, issubclassof, property_RO, \
                            _xinstanceof, _xkwds, _xzipairs
from pygeodesy.datum import Datum, R_MA
from pygeodesy.dms import clipDegrees, parseDMS2
from pygeodesy.ellipsoidalBase import LatLonEllipsoidalBase as _LLEB
from pygeodesy.errors import _IsnotError, _parseX, _TypeError, _ValueError
from pygeodesy.interns import _COMMA_, _COMMA_SPACE_, _easting_, \
                              _ellipsoidal_, NN, _northing_, _radius_, \
                              _SPACE_, _SQUARE_, _x_, _y_
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY
from pygeodesy.named import LatLon2Tuple, _NamedBase, _NamedTuple, \
                            nameof, PhiLam2Tuple, _xnamed
from pygeodesy.streprs import strs
from pygeodesy.units import Easting, Lam_, Lat, Lon, Northing, Phi_, \
                            Radius, Radius_
from pygeodesy.utily import degrees90, degrees180

from math import atan, atanh, exp, radians, sin, tanh

__all__ = _ALL_LAZY.webmercator
__version__ = '20.07.08'

# _FalseEasting  = 0   #: (INTERNAL) False Easting (C{meter}).
# _FalseNorthing = 0   #: (INTERNAL) False Northing (C{meter}).
_LatLimit = 85.051129  #: (INTERNAL) Latitudinal limit (C{degrees}).
# _LonOrigin     = 0   #: (INTERNAL) Longitude of natural origin (C{degrees}).


class EasNorRadius3Tuple(_NamedTuple):
    '''3-Tuple C{(easting, northing, radius)}, all in C{meter}.
    '''
    _Names_ = (_easting_, _northing_, _radius_)

    def __new__(cls, e, n, r):
        return _NamedTuple.__new__(cls, Easting( e, Error=WebMercatorError),
                                        Northing(n, Error=WebMercatorError),
                                        Radius(  r, Error=WebMercatorError))


class WebMercatorError(_ValueError):
    '''Web Mercator (WM) parser or L{Wm} issue.
    '''
    pass


class Wm(_NamedBase):
    '''Web Mercator (WM) coordinate.
    '''
    _latlon = None  #: (INTERNAL) cached (L{LatLon2Tuple}).
    _philam = None  #: (INTERNAL) Cached (L{PhiLam2Tuple}).
    _radius = 0     #: (INTERNAL) earth radius (C{meter}).
    _x      = 0     #: (INTERNAL) easting (C{meter}).
    _y      = 0     #: (INTERNAL) northing (C{meter}).

    def __init__(self, x, y, radius=R_MA, name=NN):
        '''New L{Wm} Web Mercator (WM) coordinate.

           @arg x: Easting from central meridian (C{meter}).
           @arg y: Northing from equator (C{meter}).
           @kwarg radius: Optional earth radius (C{meter}).
           @kwarg name: Optional name (C{str}).

           @raise WebMercatorError: Invalid B{C{x}}, B{C{y}} or B{C{radius}}.

           @example:

           >>> import pygeodesy
           >>> w = pygeodesy.Wm(448251, 5411932)
        '''
        self._x = Easting( x, name=_x_, Error=WebMercatorError)
        self._y = Northing(y, name=_y_, Error=WebMercatorError)
        self._radius = Radius_(radius, Error=WebMercatorError)

        if name:
            self.name = name

    @property_RO
    def latlon(self):
        '''Get the lat- and longitude (L{LatLon2Tuple}C{(lat, lon)}).
        '''
        if self._latlon is None:
            self._latlon = self.latlon2()
        return self._xnamed(self._latlon)

    def latlon2(self, datum=None):
        '''Convert this WM coordinate to a lat- and longitude.

           @kwarg datum: Optional ellipsoidal datum (C{Datum}).

           @return: A L{LatLon2Tuple}C{(lat, lon)}.

           @raise TypeError: Non-ellipsoidal B{C{datum}}.

           @see: Method C{toLatLon}.
        '''
        r = self.radius
        x = self._x / r
        y = 2 * atan(exp(self._y / r)) - PI_2
        if datum is not None:
            _xinstanceof(Datum, datum=datum)
            E = datum.ellipsoid
            if not E.isEllipsoidal:
                raise _IsnotError(_ellipsoidal_, datum=datum)
            # <https://Earth-Info.NGA.mil/GandG/wgs84/web_mercator/
            #        %28U%29%20NGA_SIG_0011_1.0.0_WEBMERC.pdf>
            y = y / r
            if E.e:
                y -= E.e * atanh(E.e * tanh(y))  # == E.es_atanh(tanh(y))
            y *= E.a
            x *= E.a / r

        r = LatLon2Tuple(Lat(degrees90(y)), Lon(degrees180(x)))
        return self._xnamed(r)

    def parseWM(self, strWM, name=NN):
        '''Parse a string to a WM coordinate.

           For more details, see function L{parseWM} in
           this module L{webmercator}.
        '''
        return parseWM(strWM, radius=self.radius, Wm=self.classof, name=name)

    @property_RO
    def philam(self):
        '''Get the lat- and longitude ((L{PhiLam2Tuple}C{(phi, lam)}).
        '''
        if self._philam is None:
            r = self.latlon
            self._philam = PhiLam2Tuple(Phi_(r.lat), Lam_(r.lon))
        return self._xnamed(self._philam)

    @property_RO
    def radius(self):
        '''Get the earth radius (C{meter}).
        '''
        return self._radius

    def to2ll(self, datum=None):  # PYCHOK no cover
        '''DEPRECATED, use method C{latlon2}.

           @return: A L{LatLon2Tuple}C{(lat, lon)}.
        '''
        return self.latlon2(datum=datum)

    def toLatLon(self, LatLon, datum=None, **LatLon_kwds):
        '''Convert this WM coordinate to a geodetic point.

           @arg LatLon: Ellipsoidal class to return the geodetic
                        point (C{LatLon}).
           @kwarg datum: Optional datum for ellipsoidal or C{None}
                         for spherical B{C{LatLon}} (C{Datum}).
           @kwarg LatLon_kwds: Optional, additional B{C{LatLon}}
                               keyword arguments.

           @return: Point of this WM coordinate (B{C{LatLon}}).

           @raise TypeError: If B{C{LatLon}} and B{C{datum}} are
                             incompatible or if B{C{datum}} is not
                             ellipsoidal.

           @example:

           >>> w = Wm(448251.795, 5411932.678)
           >>> from pygeodesy import sphericalTrigonometry as sT
           >>> ll = w.toLatLon(sT.LatLon)  # 43°39′11.58″N, 004°01′36.17″E
        '''
        e = issubclassof(LatLon, _LLEB)
        if e and datum:
            kwds = _xkwds(LatLon_kwds, datum=datum)
        elif not (e or datum):  # and LatLon
            kwds = LatLon_kwds
            datum = None
        else:
            raise _TypeError(LatLon=LatLon, datum=datum)

        r = self.latlon2(datum=datum)
        r = LatLon(r.lat, r.lon, **kwds)
        return self._xnamed(r)

    def toRepr(self, prec=3, fmt=_SQUARE_, sep=_COMMA_SPACE_, radius=False, **unused):  # PYCHOK expected
        '''Return a string representation of this WM coordinate.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).
           @kwarg fmt: Optional, enclosing backets format (C{str}).
           @kwarg sep: Optional separator between name:value pairs (C{str}).
           @kwarg radius: Optionally, include radius (C{bool} or C{scalar}).

           @return: This WM as "[x:meter, y:meter]" (C{str}) plus
                    ", radius:meter]" if B{C{radius}} is C{True} or
                    C{scalar}.

           @raise WebMercatorError: Invalid B{C{radius}}.
        '''
        t = self.toStr(prec=prec, sep=None, radius=radius)
        return _xzipairs((_x_, _y_, _radius_), t, sep=sep, fmt=fmt)

    toStr2 = toRepr  # PYCHOK for backward compatibility
    '''DEPRECATED, use method L{Wm.toRepr}.'''

    def toStr(self, prec=3, sep=_SPACE_, radius=False, **unused):  # PYCHOK expected
        '''Return a string representation of this WM coordinate.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).
           @kwarg sep: Optional separator to join (C{str}) or C{None}
                       to return an unjoined C{tuple} of C{str}s.
           @kwarg radius: Optionally, include radius (C{bool} or C{scalar}).

           @return: This WM as "meter meter" (C{str}) plus " radius"
                    if B{C{radius}} is C{True} or C{scalar}.

           @raise WebMercatorError: Invalid B{C{radius}}.

           @example:

           >>> w = Wm(448251, 5411932.0001)
           >>> w.toStr(4)  # 448251.0 5411932.0001
           >>> w.toStr(sep=', ')  # 448251, 5411932
        '''
        fs = self._x, self._y
        if radius in (False, None):
            pass
        elif radius is True:
            fs += (self._radius,)
        elif isscalar(radius):
            fs += (radius,)
        else:
            raise WebMercatorError(radius=radius)
        t = strs(fs, prec=prec)
        return t if sep is None else sep.join(t)

    @property_RO
    def x(self):
        '''Get the easting (C{meter}).
        '''
        return self._x

    @property_RO
    def y(self):
        '''Get the northing (C{meter}).
        '''
        return self._y


def parseWM(strWM, radius=R_MA, Wm=Wm, name=NN):
    '''Parse a string representing a WM coordinate, consisting
       of easting, northing and an optional radius.

       @arg strWM: A WM coordinate (C{str}).
       @kwarg radius: Optional earth radius (C{meter}).
       @kwarg Wm: Optional class to return the WM coordinate (L{Wm})
                  or C{None}.
       @kwarg name: Optional name (C{str}).

       @return: The WM coordinate (B{C{Wm}}) or an
                L{EasNorRadius3Tuple}C{(easting, northing, radius)}
                if B{C{Wm}} is C{None}.

       @raise WebMercatorError: Invalid B{C{strWM}}.

       @example:

       >>> u = parseWM('448251 5411932')
       >>> u.toStr2()  # [E:448251, N:5411932]
    '''
    def _WM_(strWM, radius, Wm, name):
        w = strWM.replace(_COMMA_, _SPACE_).strip().split()

        if len(w) == 2:
            w += [radius]
        elif len(w) != 3:
            raise ValueError
        x, y, r = map(float, w)

        r = EasNorRadius3Tuple(x, y, r) if Wm is None else \
                            Wm(x, y, radius=r)
        return _xnamed(r, name)

    return _parseX(_WM_, strWM, radius, Wm, name,
                         strWM=strWM, Error=WebMercatorError)


def toWm(latlon, lon=None, radius=R_MA, Wm=Wm, name=NN, **Wm_kwds):
    '''Convert a lat-/longitude point to a WM coordinate.

       @arg latlon: Latitude (C{degrees}) or an (ellipsoidal or
                    spherical) geodetic C{LatLon} point.
       @kwarg lon: Optional longitude (C{degrees} or C{None}).
       @kwarg radius: Optional earth radius (C{meter}).
       @kwarg Wm: Optional class to return the WM coordinate
                  (L{Wm}) or C{None}.
       @kwarg name: Optional name (C{str}).
       @kwarg Wm_kwds: Optional, additional B{C{Wm}} keyword
                       arguments, ignored if B{C{Wm=None}}.

       @return: The WM coordinate (B{C{Wm}}) or an
                L{EasNorRadius3Tuple}C{(easting, northing, radius)}
                if B{C{Wm}} is C{None}.

       @raise ValueError: If B{C{lon}} value is missing, if B{C{latlon}}
                          is not scalar, if B{C{latlon}} is beyond the
                          valid WM range and L{rangerrors} is set
                          to C{True} or if B{C{radius}} is invalid.

       @example:

       >>> p = LatLon(48.8582, 2.2945)  # 448251.8 5411932.7
       >>> w = toWm(p)  # 448252 5411933
       >>> p = LatLon(13.4125, 103.8667)  # 377302.4 1483034.8
       >>> w = toWm(p)  # 377302 1483035
    '''
    e, r = None, Radius(radius)
    try:
        lat, lon = latlon.lat, latlon.lon
        if isinstance(latlon, _LLEB):
            r = latlon.datum.ellipsoid.a
            e = latlon.datum.ellipsoid.e
            if not name:  # use latlon.name
                name = nameof(latlon)
        lat = clipDegrees(lat, _LatLimit)
    except AttributeError:
        lat, lon = parseDMS2(latlon, lon, clipLat=_LatLimit)

    s = sin(radians(lat))
    y = atanh(s)  # == log(tan((90 + lat) / 2)) == log(tanPI_2_2(radians(lat)))
    if e:
        y -= e * atanh(e * s)

    e = r * radians(lon)
    n = r * y
    r = EasNorRadius3Tuple(e, n, r) if Wm is None else \
                        Wm(e, n, **_xkwds(Wm_kwds, radius=r))
    return _xnamed(r, name)


__all__ += _ALL_DOCS(EasNorRadius3Tuple)

# **) MIT License
#
# Copyright (C) 2016-2020 -- mrJean1 at Gmail -- All Rights Reserved.
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
