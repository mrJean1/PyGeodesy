
# -*- coding: utf-8 -*-

u'''Web Mercator (WM) projection.

Classes L{Wm} and L{WebMercatorError} and functions L{parseWM} and L{toWm}.

Pure Python implementation of a U{Web Mercator<https://WikiPedia.org/wiki/Web_Mercator>}
(aka I{Pseudo-Mercator}) class and conversion functions for spherical and
near-spherical earth models.

References U{Google Maps / Bing Maps Spherical Mercator Projection
<https://AlastairA.WordPress.com/2011/01/23/the-google-maps-bing-maps-spherical-mercator-projection>},
U{Geomatics Guidance Note 7, part 2<https://www.EPSG.org/Portals/0/373-07-02.pdf>} and
U{Implementation Practice Web Mercator Map Projection
<https://Earth-Info.NGA.mil/GandG/wgs84/web_mercator/%28U%29%20NGA_SIG_0011_1.0.0_WEBMERC.pdf>}.
'''

from pygeodesy.basics import isscalar, issubclassof
from pygeodesy.datums import _ellipsoidal_datum
from pygeodesy.dms import clipDegrees, parseDMS2
from pygeodesy.ellipsoidalBase import LatLonEllipsoidalBase as _LLEB
from pygeodesy.ellipsoids import R_MA
from pygeodesy.errors import _IsnotError, _parseX, _TypeError, \
                             _ValueError, _xkwds
from pygeodesy.interns import NN, PI_2, _COMMA_, _COMMASPACE_, \
                             _easting_, _ellipsoidal_, _northing_, \
                             _radius_, _SPACE_
from pygeodesy.interns import _x_, _y_  # PYCHOK used!
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.named import _NamedBase, _NamedTuple, nameof
from pygeodesy.namedTuples import LatLon2Tuple, PhiLam2Tuple
from pygeodesy.props import deprecated_method, Property_RO
from pygeodesy.streprs import Fmt, strs, _xzipairs
from pygeodesy.units import Easting, Lam_, Lat, Lon, Northing, Phi_, \
                            Radius, Radius_
from pygeodesy.utily import degrees90, degrees180

from math import atan, atanh, exp, radians, sin, tanh

__all__ = _ALL_LAZY.webmercator
__version__ = '21.08.10'

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
    _radius = 0  # earth radius (C{meter})
    _x      = 0  # Easting (C{meter})
    _y      = 0  # Northing (C{meter})

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
        self._x = Easting( x=x, Error=WebMercatorError)
        self._y = Northing(y=y, Error=WebMercatorError)
        self._radius = Radius_(radius, Error=WebMercatorError)

        if name:
            self.name = name

    @Property_RO
    def latlon(self):
        '''Get the lat- and longitude (L{LatLon2Tuple}C{(lat, lon)}).
        '''
        return self.latlon2()

    def latlon2(self, datum=None):
        '''Convert this WM coordinate to a lat- and longitude.

           @kwarg datum: Optional, ellipsoidal datum (L{Datum},
                         L{Ellipsoid}, L{Ellipsoid2} or
                         L{a_f2Tuple}) or C{None}.

           @return: A L{LatLon2Tuple}C{(lat, lon)}.

           @raise TypeError: Invalid or non-ellipsoidal B{C{datum}}.

           @see: Method C{toLatLon}.
        '''
        r = self.radius
        x = self._x / r
        y = 2 * atan(exp(self._y / r)) - PI_2
        if datum is not None:
            E = _ellipsoidal_datum(datum, name=self.name).ellipsoid
            if not E.isEllipsoidal:
                raise _IsnotError(_ellipsoidal_, datum=datum)
            # <https://Earth-Info.NGA.mil/GandG/wgs84/web_mercator/
            #        %28U%29%20NGA_SIG_0011_1.0.0_WEBMERC.pdf>
            y = y / r
            if E.e:
                y -= E.e * atanh(E.e * tanh(y))  # == E.es_atanh(tanh(y))
            y *= E.a
            x *= E.a / r

        return LatLon2Tuple(Lat(degrees90( y)),
                            Lon(degrees180(x)), name=self.name)

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
        '''Get the lat- and longitude ((L{PhiLam2Tuple}C{(phi, lam)}).
        '''
        r = self.latlon
        return PhiLam2Tuple(Phi_(r.lat), Lam_(r.lon), name=r.name)

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
                             incompatible or if B{C{datum}} is
                             invalid or not ellipsoidal.

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
        return LatLon(r.lat, r.lon, **_xkwds(kwds, name=r.name))

    def toRepr(self, prec=3, fmt=Fmt.SQUARE, sep=_COMMASPACE_, radius=False, **unused):  # PYCHOK expected
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


def parseWM(strWM, radius=R_MA, Wm=Wm, name=NN):
    '''Parse a string C{"e n [r]"} representing a WM coordinate,
       consisting of easting, northing and an optional radius.

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
        >>> u.toRepr()  # [E:448251, N:5411932]
    '''
    def _WM_(strWM, radius, Wm, name):
        w = strWM.replace(_COMMA_, _SPACE_).strip().split()

        if len(w) == 2:
            w += [radius]
        elif len(w) != 3:
            raise ValueError
        x, y, r = map(float, w)

        return EasNorRadius3Tuple(x, y, r, name=name) if Wm is None else \
                               Wm(x, y, radius=r, name=name)

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
                       arguments, ignored if C{B{Wm} is None}.

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
    r = EasNorRadius3Tuple(e, n, r, name=name) if Wm is None else \
                        Wm(e, n, **_xkwds(Wm_kwds, radius=r, name=name))
    return r

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
