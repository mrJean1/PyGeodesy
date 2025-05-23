
# -*- coding: utf-8 -*-

u'''World Geographic Reference System (WGRS) en-/decoding, aka GEOREF.

Class L{Georef} and several functions to encode, decode and inspect WGRS
(or GEOREF) references.

Transcoded from I{Charles Karney}'s C++ class U{Georef
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Georef.html>},
but with modified C{precision} and extended with C{height} and C{radius}.

@see: U{World Geographic Reference System
      <https://WikiPedia.org/wiki/World_Geographic_Reference_System>}.
'''
from pygeodesy.basics import isstr,  typename
from pygeodesy.constants import INT0, _float, _off90, _0_001, \
                               _0_5, _1_0, _2_0, _60_0, _1000_0
from pygeodesy.dms import parse3llh
from pygeodesy.errors import _ValueError, _xattr, _xStrError
# from pygeodesy.internals import typename  # from .basics
from pygeodesy.interns import NN, _0to9_, _AtoZnoIO_, _COMMA_, \
                             _height_, _INV_, _radius_, _SPACE_
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY
from pygeodesy.named import _name2__, nameof,  Property_RO
from pygeodesy.namedTuples import LatLon2Tuple, LatLonPrec3Tuple
# from pygeodesy.props import Property_RO  # from .named
from pygeodesy.streprs import Fmt, _0wd
from pygeodesy.units import Height, Int, Lat, Lon, Precision_, \
                            Radius, Scalar_, Str
from pygeodesy.utily import ft2m, m2ft, m2NM

from math import floor

__all__ = _ALL_LAZY.wgrs
__version__ = '25.04.14'

_Base    =  10
_BaseLen =  4
_DegChar = _AtoZnoIO_.tillQ
_Digits  = _0to9_
_LatOrig = -90
_LatTile = _AtoZnoIO_.tillM
_LonOrig = -180
_LonTile = _AtoZnoIO_
_60B     =  60000000000  # == 60_000_000_000 == 60e9
_MaxPrec =  11
_Tile    =  15  # tile size in degrees

_MaxLen  = _BaseLen + 2 * _MaxPrec
_MinLen  = _BaseLen - 2

_LatOrig_60B = _LatOrig * _60B
_LonOrig_60B = _LonOrig * _60B

_float_Tile   = _float(_Tile)
_LatOrig_Tile = _float(_LatOrig) / _Tile
_LonOrig_Tile = _float(_LonOrig) / _Tile


def _divmod3(x, _Orig_60B):
    '''(INTERNAL) Convert B{C{x}} to 3_tuple C{(tile, modulo, fraction)}/
    '''
    i = int(floor(x * _60B))
    i, x = divmod(i - _Orig_60B, _60B)
    xt, xd = divmod(i, _Tile)
    return xt, xd, x


def _2fll(lat, lon):
    '''(INTERNAL) Convert lat, lon.
    '''
    # lat, lon = parseDMS2(lat, lon)
    return (Lat(lat, Error=WGRSError),
            Lon(lon, Error=WGRSError))


def _2geostr2(georef):
    '''(INTERNAL) Check a georef string.
    '''
    try:
        n, g = len(georef), georef.upper()
        p, o = divmod(n, 2)
        if o or n < _MinLen or n > _MaxLen \
             or g.startswith(_INV_) or not g.isalnum():
            raise ValueError()
        return g, _2Precision(p - 1)

    except (AttributeError, TypeError, ValueError) as x:
        raise WGRSError(typename(Georef), georef, cause=x)


def _2Precision(precision):
    '''(INTERNAL) Return a L{Precision_} instance.
    '''
    return Precision_(precision, Error=WGRSError, low=0, high=_MaxPrec)


class WGRSError(_ValueError):
    '''World Geographic Reference System (WGRS) encode, decode or other L{Georef} issue.
    '''
    pass


class Georef(Str):
    '''Georef class, a named C{str}.
    '''
    # no str.__init__ in Python 3
    def __new__(cls, lat_gll, lon=None, height=None, precision=3, name=NN):
        '''New L{Georef} from an other L{Georef} instance or georef
           C{str} or from a C{LatLon} instance or lat-/longitude C{str}.

           @arg lat_gll: Latitude (C{degrees90}), a georef (L{Georef},
                         C{str}) or a location (C{LatLon}, C{LatLon*Tuple}).
           @kwarg lon: Logitude (C{degrees180)}, required if B{C{lat_gll}}
                       is C{degrees90}, ignored otherwise.
           @kwarg height: Optional height in C{meter}, used if B{C{lat_gll}}
                          is a location.
           @kwarg precision: The desired georef resolution and length (C{int}
                             0..11), see L{encode<pygeodesy.wgrs.encode>}.
           @kwarg name: Optional name (C{str}).

           @return: New L{Georef}.

           @raise RangeError: Invalid B{C{lat_gll}} or B{C{lon}}.

           @raise TypeError: Invalid B{C{lat_gll}} or B{C{lon}}.

           @raise WGRSError: INValid B{C{lat_gll}}.
        '''
        if lon is None:
            if isinstance(lat_gll, Georef):
                g, ll, p = str(lat_gll), lat_gll.latlon, lat_gll.precision
            elif isstr(lat_gll):
                if _COMMA_ in lat_gll or _SPACE_ in lat_gll:
                    lat, lon, h = parse3llh(lat_gll, height=height)
                    g, ll, p = _encode3(lat, lon, precision, h=h)
                else:
                    g, ll = lat_gll.upper(), None
                    try:
                        _, p = _2geostr2(g)  # validate
                    except WGRSError:  # R00H00?
                        p = None  # = decode5(g).precision?
            else:  # assume LatLon
                try:
                    g, ll, p = _encode3(lat_gll.lat, lat_gll.lon, precision,
                                        h=_xattr(lat_gll, height=height))
                except AttributeError:
                    raise _xStrError(Georef, gll=lat_gll)  # Error=WGRSError
        else:
            g, ll, p = _encode3(lat_gll, lon, precision, h=height)

        self = Str.__new__(cls, g, name=name or nameof(lat_gll))
        self._latlon    = ll
        self._precision = p
        return self

    @Property_RO
    def decoded3(self):
        '''Get this georef's attributes (L{LatLonPrec3Tuple}).
        '''
        lat, lon = self.latlon
        return LatLonPrec3Tuple(lat, lon, self.precision, name=self.name)

    @Property_RO
    def decoded5(self):
        '''Get this georef's attributes (L{LatLonPrec5Tuple}) with
           height and radius set to C{None} if missing.
        '''
        return self.decoded3.to5Tuple(self.height, self.radius)

    @Property_RO
    def _decoded5(self):
        '''(INTERNAL) Initial L{LatLonPrec5Tuple}.
        '''
        return decode5(self)

    @Property_RO
    def height(self):
        '''Get this georef's height in C{meter} or C{None} if missing.
        '''
        return self._decoded5.height

    @Property_RO
    def latlon(self):
        '''Get this georef's (center) lat- and longitude (L{LatLon2Tuple}).
        '''
        lat, lon = self._latlon or self._decoded5[:2]
        return LatLon2Tuple(lat, lon, name=self.name)

    @Property_RO
    def latlonheight(self):
        '''Get this georef's (center) lat-, longitude and height (L{LatLon3Tuple}),
           with height set to C{INT0} if missing.
        '''
        return self.latlon.to3Tuple(self.height or INT0)

    @Property_RO
    def precision(self):
        '''Get this georef's precision (C{int}).
        '''
        p = self._precision
        return self._decoded5.precision if p is None else p

    @Property_RO
    def radius(self):
        '''Get this georef's radius in C{meter} or C{None} if missing.
        '''
        return self._decoded5.radius

    def toLatLon(self, LatLon=None, height=None, **name_LatLon_kwds):
        '''Return (the center of) this georef cell as a C{LatLon}.

           @kwarg LatLon: Class to use (C{LatLon}) or C{None}.
           @kwarg height: Optional height (C{meter}), overriding this height.
           @kwarg name_LatLon_kwds: Optional C{B{name}=NN} (C{str}) and optionally,
                       additional B{C{LatLon}} keyword arguments, ignored if C{B{LatLon}
                       is None}.

           @return: This georef location (B{C{LatLon}}) or if C{B{LatLon} is None},
                    a L{LatLon3Tuple}C{(lat, lon, height)}.

           @raise TypeError: Invalid B{C{LatLon}} or B{C{name_LatLon_kwds}}.
        '''
        n, kwds = _name2__(name_LatLon_kwds, _or_nameof=self)
        h = (self.height or INT0) if height is None else height  # _heigHt
        r =  self.latlon.to3Tuple(h) if LatLon is None else LatLon(
            *self.latlon, height=h, **kwds)
        return r.renamed(n) if n else r


def decode3(georef, center=True):
    '''Decode a C{georef} to lat-, longitude and precision.

       @arg georef: To be decoded (L{Georef} or C{str}).
       @kwarg center: If C{True}, use the georef's center, otherwise
                      the south-west, lower-left corner (C{bool}).

       @return: A L{LatLonPrec3Tuple}C{(lat, lon, precision)}.

       @raise WGRSError: Invalid B{C{georef}}, INValid, non-alphanumeric
                         or odd length B{C{georef}}.
    '''
    def _digit(ll, g, i, m):
        d = _Digits.find(g[i])
        if d < 0 or d >= m:
            raise _Error(i)
        return ll * m + d

    def _Error(i):
        return WGRSError(Fmt.SQUARE(georef=i), georef)

    def _index(chars, g, i):
        k = chars.find(g[i])
        if k < 0:
            raise _Error(i)
        return k

    g, precision = _2geostr2(georef)
    lon = _index(_LonTile, g, 0) + _LonOrig_Tile
    lat = _index(_LatTile, g, 1) + _LatOrig_Tile

    u = _1_0
    if precision > 0:
        lon = lon * _Tile + _index(_DegChar, g, 2)
        lat = lat * _Tile + _index(_DegChar, g, 3)
        m, p = 6, precision - 1
        for i in range(_BaseLen, _BaseLen + p):
            lon = _digit(lon, g, i,     m)
            lat = _digit(lat, g, i + p, m)
            u *= m
            m  = _Base
        u *= _Tile

    if center:
        lon = lon * _2_0 + _1_0
        lat = lat * _2_0 + _1_0
        u *= _2_0
    u = _Tile / u
    return LatLonPrec3Tuple(Lat(lat * u, Error=WGRSError),
                            Lon(lon * u, Error=WGRSError),
                            precision, name=nameof(georef))


def decode5(georef, center=True):
    '''Decode a C{georef} to lat-, longitude, precision, height and radius.

       @arg georef: To be decoded (L{Georef} or C{str}).
       @kwarg center: If C{True}, use the georef's center, otherwise the
                      south-west, lower-left corner (C{bool}).

       @return: A L{LatLonPrec5Tuple}C{(lat, lon, precision, height, radius)}
                where C{height} and/or C{radius} are C{None} if missing.

       @raise WGRSError: Invalid B{C{georef}}.
    '''
    def _h2m(kft, g_n):
        return Height(ft2m(kft * _1000_0), name=g_n, Error=WGRSError)

    def _r2m(NM, g_n):
        return Radius(NM / m2NM(1), name=g_n, Error=WGRSError)

    def _split2(g, Unit, _2m):
        n = typename(Unit)
        i = max(g.find(n[0]), g.rfind(n[0]))
        if i > _BaseLen:
            return g[:i], _2m(int(g[i+1:]), _SPACE_(georef, n))
        else:
            return g, None

    g = Str(georef, Error=WGRSError)

    g, h = _split2(g, Height, _h2m)  # H is last
    g, r = _split2(g, Radius, _r2m)  # R before H

    return decode3(g, center=center).to5Tuple(h, r)


def encode(lat, lon, precision=3, height=None, radius=None):  # MCCABE 14
    '''Encode a lat-/longitude as a C{georef} of the given precision.

       @arg lat: Latitude (C{degrees}).
       @arg lon: Longitude (C{degrees}).
       @kwarg precision: Optional, the desired C{georef} resolution and length
                         (C{int} 0..11).
       @kwarg height: Optional, height in C{meter}, see U{Designation of area
                      <https://WikiPedia.org/wiki/World_Geographic_Reference_System>}.
       @kwarg radius: Optional, radius in C{meter}, see U{Designation of area
                      <https://WikiPedia.org/wiki/World_Geographic_Reference_System>}.

       @return: The C{georef} (C{str}).

       @raise RangeError: Invalid B{C{lat}} or B{C{lon}}.

       @raise WGRSError: Invalid B{C{precision}}, B{C{height}} or B{C{radius}}.

       @note: The B{C{precision}} value differs from I{Karney}'s U{Georef<https://
              GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Georef.html>}.
              The C{georef} length is M{2 * (precision + 1)} and the C{georef}
              resolution is I{15°} for B{C{precision}} 0:, I{1°} for 1, I{1′} for 2,
              I{0.1′} for 3, I{0.01′} for 4, ... up to I{10**(2 - precision)′}.
    '''
    def _option(name, m, m2_, K):
        f = Scalar_(m, name=name, Error=WGRSError)
        return NN(name[0].upper(), int(m2_(f * K) + _0_5))

    g, _, _ = _encode3(lat, lon, precision)
    if radius is not None:  # R before H
        g += _option(_radius_, radius, m2NM, _1_0)
    if height is not None:  # H is last
        g += _option(_height_, height, m2ft, _0_001)
    return g


def _encode3(lat, lon, precision, h=None):
    '''Return 3-tuple C{(georef, (lat, lon), p)}.
    '''
    p = _2Precision(precision)

    lat, lon = _2fll(lat, lon)

    if h is None:
        xt, xd, x = _divmod3(       lon,  _LonOrig_60B)
        yt, yd, y = _divmod3(_off90(lat), _LatOrig_60B)

        g = _LonTile[xt], _LatTile[yt]
        if p > 0:
            g += _DegChar[xd], _DegChar[yd]
            p -= 1
            if p > 0:
                d =  pow(_Base, _MaxPrec - p)
                x = _0wd(p, x // d)
                y = _0wd(p, y // d)
                g += x, y
        g = NN.join(g)
    else:
        g = encode(lat, lon, precision=p, height=h)

    return g, (lat, lon), p  # XXX Georef(''.join(g))


def precision(res):
    '''Determine the L{Georef} precision to meet a required (geographic)
       resolution.

       @arg res: The required resolution (C{degrees}).

       @return: The L{Georef} precision (C{int} 0..11).

       @raise ValueError: Invalid B{C{res}}.

       @see: Function L{wgrs.encode} for more C{precision} details.
    '''
    r = Scalar_(res=res)
    for p in range(_MaxPrec):
        if resolution(p) <= r:
            return p
    return _MaxPrec


def resolution(prec):
    '''Determine the (geographic) resolution of a given L{Georef} precision.

       @arg prec: The given precision (C{int}).

       @return: The (geographic) resolution (C{degrees}).

       @raise ValueError: Invalid B{C{prec}}.

       @see: Function L{wgrs.encode} for more C{precision} details.
    '''
    p = Int(prec=prec, Error=WGRSError)
    if p > 1:
        p =  min(p, _MaxPrec) - 1
        r = _1_0 / (pow(_Base, p) * _60_0)
    elif p < 1:
        r = _float_Tile
    else:
        r = _1_0
    return r


__all__ += _ALL_DOCS(decode3, decode5,  # functions
                     encode, precision, resolution)

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
