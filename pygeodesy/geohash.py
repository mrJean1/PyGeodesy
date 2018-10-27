
# -*- coding: utf-8 -*-

u'''Class L{Geohash} and several functions to encode, decode and
inspect I{geohashes}.

Transcribed from JavaScript originals by I{(C) Chris Veness 2011-2015}
and published under the same MIT Licence**, see U{Geohashes
<http://www.Movable-Type.co.UK/scripts/geohash.html>}.

See also U{Geohash<http://WikiPedia.org/wiki/Geohash>},
U{Geohash<http://GitHub.com/vinsci/geohash>},
U{PyGeohash<http://PyPI.org/project/pygeohash>} and
U{Geohash-Javascript<http://GitHub.com/DaveTroy/geohash-js>}.

@newfield example: Example, Examples
'''

from bases import Named, _xnamed
from dms import parse3llh, parseDMS2
from fmath import EPS, favg, fStr, map2
from formy import equirectangular, equirectangular_, haversine_
from utily import R_M, property_RO, _Strs, unrollPI

from math import log10, radians

# all public contants, classes and functions
__all__ = ('Geohash',  # classes
           'bounds', 'decode', 'decode_error',  # functions
           'distance1', 'distance2', 'distance3',
           'encode', 'neighbors', 'sizes')
__version__ = '18.10.26'

_Border = dict(
    N=('prxz',     'bcfguvyz'),
    S=('028b',     '0145hjnp'),
    E=('bcfguvyz', 'prxz'),
    W=('0145hjnp', '028b'))

_Neighbor = dict(
    N=('p0r21436x8zb9dcf5h7kjnmqesgutwvy', 'bc01fg45238967deuvhjyznpkmstqrwx'),
    S=('14365h7k9dcfesgujnmqp0r2twvyx8zb', '238967debc01fg45kmstqrwxuvhjyznp'),
    E=('bc01fg45238967deuvhjyznpkmstqrwx', 'p0r21436x8zb9dcf5h7kjnmqesgutwvy'),
    W=('238967debc01fg45kmstqrwxuvhjyznp', '14365h7k9dcfesgujnmqp0r2twvyx8zb'))

# lat-, longitudinal and radial cell size (in meter)
_Sizes = (  # radius = sqrt(latHeight * lonWidth / PI)
    (20032e3, 20000e3, 11292815.096),  # 0
    ( 5003e3,  5000e3,  2821794.075),  # 1
    (  650e3,  1225e3,   503442.397),  # 2
    (  156e3,   156e3,    88013.575),  # 3
    (  19500,   39100,    15578.683),  # 4
    (   4890,    4890,     2758.887),  # 5
    (    610,    1220,      486.710),  # 6
    (    153,     153,       86.321),  # 7
    (     19.1,    38.2,     15.239),  # 8
    (      4.77,    4.77,     2.691),  # 9
    (      0.596,   1.19,     0.475),  # 10
    (      0.149,   0.149,    0.084),  # 11
    (      0.0186,  0.0372,   0.015))  # 12

# Geohash-specific base32 map
_GeohashBase32 = '0123456789bcdefghjkmnpqrstuvwxyz'
# ... and the inverse map
_DecodedBase32 = dict((c, i) for i, c in enumerate(_GeohashBase32))
c = i = None
del c, i


def _2fll(lat, lon, *unused):
    '''(INTERNAL) Convert lat, lon to 2-tuple of floats.
    '''
    return parseDMS2(lat, lon)


def _2Geohash(geohash):
    '''(INTERNAL) Check or create a Geohash instance.
    '''
    if not isinstance(geohash, Geohash):
        try:
            geohash = Geohash(geohash)
        except (TypeError, ValueError):
            raise TypeError('%r not %s, %s or str' % (geohash,
                             Geohash.__name__, 'LatLon'))
    return geohash


def _2geostr(geohash):
    '''(INTERNAL) Check a geohash string.
    '''
    try:
        if not (0 < len(geohash) < 13):
            raise ValueError
        geostr = geohash.lower()
        for c in geostr:
            if c not in _DecodedBase32:
                raise ValueError
    except (AttributeError, ValueError, TypeError):
        raise ValueError('%s: %r[%s]' % (Geohash.__name__,
                          geohash, len(geohash)))
    return geostr


class Geohash(str, Named):
    '''Geohash class, sub-class of str.
    '''
    _bounds = None  # cached bounds property
    _latlon = None  # cached latlon property

    _N = None  # cached neighbors properties
    _E = None
    _S = None
    _W = None

    # no str.__init__ in Python 3
    def __new__(cls, cll, precision=None, name=''):
        '''New Geohash from a L{Geohash} instance or C{str} or from
           a C{LatLon} instance or C{str}.

           @param cll: Cell or location (L{Geohash} or C{str}, C{LatLon}
                       or C{str}).
           @keyword precision: Optional, desired geohash length (C{int}),
                               see function L{geohash.encode} for more
                               details.
           @keyword name: Optional name (C{str}).

           @return: New L{Geohash}.
        '''
        if isinstance(cll, Geohash):
            self = str.__new__(cls, _2geostr('%s' % (cll,)))

        elif isinstance(cll, _Strs):
            if ',' in cll:
                lat, lon = _2fll(*parse3llh(cll))
                cll = encode(lat, lon, precision=precision)
                self = str.__new__(cls, cll)
                self._latlon = lat, lon
            else:
                self = str.__new__(cls, _2geostr(cll))

        else:  # assume LatLon
            try:
                lat, lon = _2fll(cll.lat, cll.lon)
            except AttributeError:
                raise TypeError('%s: %r' % (Geohash.__name__, cll))
            cll = encode(lat, lon, precision=precision)
            self = str.__new__(cls, cll)
            self._latlon = lat, lon

        if name:
            self.name = name
        return self

    def __repr__(self):
        return "%s('%s')" % (Geohash.__name__, self)  # str.__str__(self))

    @property_RO
    def ab(self):
        '''Get the lat- and longitude of (the approximate center of)
           this geohash as a 2-tuple (lat, lon) in radians.
        '''
        return map2(radians, self.latlon)

    def adjacent(self, direction):
        '''Determine the adjacent cell in given compass direction.

           @param direction: Compass direction ('N', 'S', 'E' or 'W').

           @return: Geohash of adjacent cell (L{Geohash}).

           @raise ValueError: If this geohash or I{direction} invalid.
        '''
        # based on <http://GitHub.com/DaveTroy/geohash-js>

        d = direction[:1].upper()
        if d not in _Neighbor:
            raise ValueError('%s invalid: %s' % ('direction', direction))

        e = len(self) & 1  # % 2

        c = self[-1:]  # last hash char
        i = _Neighbor[d][e].find(c)
        if i < 0:
            raise ValueError('%s invalid: %s' % ('geohash', self))

        p = self[:-1]  # hash without last char
        # check for edge-cases which don't share common prefix
        if p and (c in _Border[d][e]):
            p = Geohash(p).adjacent(d)

        # append letter for direction to parent
        return Geohash(p + _GeohashBase32[i])

    def bounds(self, LatLon, **kwds):
        '''Return the SW and NE bounds of this geohash cell.

           @param LatLon: Class to use (C{LatLon}).
           @keyword kwds: Optional keyword arguments for I{LatLon}.

           @return: 2-Tuple (LatLonSW, LatLonNE) of I{LatLon} for
                    the lower-left respectively upper-right corner.
        '''
        if not self._bounds:
            self._bounds = bounds(self)
        s, w, n, e = self._bounds
        return LatLon(s, w, **kwds), LatLon(n, e, **kwds)

    def distance1(self, other):
        '''DEPRECATED, use method C{distance1To}.
        '''
        return self.distance1To(other)

    def distance1To(self, other):
        '''Estimate the distance between this and an other geohash
           (from the cell sizes).

           @param other: The other geohash (L{Geohash}).

           @return: Approximate distance (C{meter}).

           @raise TypeError: The I{other} is not a L{Geohash}, C{LatLon}
                             or C{str}.
        '''
        other = _2Geohash(other)

        n = min(len(self), len(other), len(_Sizes))
        for n in range(n):
            if self[n] != other[n]:
                break
        return float(_Sizes[n][2])

    def distance2(self, other, radius=R_M, adjust=False, wrap=False):
        '''DEPRECATED, use method C{distance2To}.
        '''
        return self.distance2To(other, radius=radius, adjust=adjust, wrap=wrap)

    def distance2To(self, other, radius=R_M, adjust=False, wrap=False):
        '''Compute the distance between this and an other geohash
           using the U{Equirectangular Approximation / Projection
           <http://www.Movable-Type.co.UK/scripts/latlong.html>}.

           @param other: The other geohash (L{Geohash}).
           @keyword radius: Optional, mean earth radius (C{meter}) or
                            C{None}.
           @keyword adjust: Adjust the wrapped, unrolled longitudinal
                            delta by the cosine of the mean latitude
                            (C{bool}).
           @keyword wrap: Wrap and unroll longitudes (C{bool}).

           @return: Approximate distance (C{meter}, same units as I{radius})
                    or distance squared (C{degrees} squared) if I{radius}
                    is C{None} or 0.

           @raise TypeError: The I{other} is not a L{Geohash}, C{LatLon}
                             or C{str}.

           @see: U{Local, flat earth approximation
                 <http://www.EdWilliams.org/avform.htm#flat>}, functions
        '''
        other = _2Geohash(other)

        a1, b1 = self.latlon
        a2, b2 = other.latlon
        if radius:
            return equirectangular(a1, b1, a2, b2, radius=radius,
                                   adjust=adjust, limit=None, wrap=wrap)
        else:
            return equirectangular_(a1, b1, a2, b2,
                                   adjust=adjust, limit=None, wrap=wrap)[0]

    def distance3(self, other, radius=R_M, wrap=False):
        '''DEPRECATED, use method C{distance3To}.
        '''
        return self.distance3To(other, radius=radius, wrap=wrap)

    def distance3To(self, other, radius=R_M, wrap=False):
        '''Compute the great-circle distance between this and an other
           geohash using the U{Haversine
           <http://www.Movable-Type.co.UK/scripts/latlong.html>} formula.

           @param other: The other geohash (L{Geohash}).
           @keyword radius: Optional, mean earth radius (C{meter}).
           @keyword wrap: Wrap and unroll longitudes (C{bool}).

           @return: Great-circle distance (C{meter}, same units as I{radius}).

           @raise TypeError: The I{other} is not a L{Geohash}, C{LatLon}
                             or C{str}.
        '''
        other = _2Geohash(other)

        a1, b1 = self.ab
        a2, b2 = other.ab

        db, b2 = unrollPI(b1, b2, wrap=wrap)
        return haversine_(a2, a1, db) * float(radius)

    @property_RO
    def latlon(self):
        '''Get the lat- and longitude of (the approximate center of)
           this geohash as a 2-tuple (lat, lon) in degrees.

           B{Example:}

           >>> geohash.Geohash('geek').latlon  # 65.478515625, -17.75390625
           >>> geohash.decode('geek')  # '65.48', '-17.75'
        '''
        if not self._latlon:
            if not self._bounds:
                self._bounds = bounds(self)
            s, w, n, e = self._bounds
            self._latlon = favg(n, s), favg(e, w)
        return self._latlon

    @property_RO
    def neighbors(self):
        '''Get all 8 adjacent cells as a dict(N=, NE=, E= ..., SW=)
           of L{Geohash}es.

           B{JSname:} I{neighbours}.
        '''
        return dict((d, getattr(self, d)) for d in ('N', 'NE', 'E', 'SE',
                                                    'S', 'SW', 'W', 'NW'))

    @property_RO
    def sizes(self):
        '''Get the lat- and longitudinal size of this cell as
           a 2-tuple (latHeight, lonWidth) in meter.
        '''
        n = min(len(_Sizes) - 1, len(self) or 1)
        return map2(float, _Sizes[n][:2])

    def toLatLon(self, LatLon, **kwds):
        '''Return (the approximate center of) this geohash cell
           as an instance of the supplied C{LatLon} class.

           @param LatLon: Class to use (C{LatLon}).
           @keyword kwds: Optional keyword arguments for I{LatLon}.

           @return: This geohash location (I{LatLon}).

           @example:

           >>> from sphericalTrigonometry import LatLon
           >>> ll = Geohash('u120fxw').toLatLon(LatLon)
           >>> print(repr(ll))  # LatLon(52°12′17.9″N, 000°07′07.64″E)
           >>> print(ll)  # 52.204971°N, 000.11879°E
        '''
        if LatLon is None:
            raise ValueError('%s invalid: %r' % ('LatLon', LatLon))

        return _xnamed(LatLon(*self.latlon, **kwds), self.name)

    @property_RO
    def N(self):
        '''Get the cell North of this (L{Geohash}).
        '''
        if self._N is None:
            self._N = self.adjacent('N')
        return self._N

    @property_RO
    def S(self):
        '''Get the cell South of this (L{Geohash}).
        '''
        if self._S is None:
            self._S = self.adjacent('S')
        return self._S

    @property_RO
    def E(self):
        '''Get the cell East of this (L{Geohash}).
        '''
        if self._E is None:
            self._E = self.adjacent('E')
        return self._E

    @property_RO
    def W(self):
        '''Get the cell West of this (L{Geohash}).
        '''
        if self._W is None:
            self._W = self.adjacent('W')
        return self._W

    @property_RO
    def NE(self):
        '''Get the cell NorthEast of this (L{Geohash}).
        '''
        return self.N.E

    @property_RO
    def NW(self):
        '''Get the cell NorthWest of this (L{Geohash}).
        '''
        return self.N.W

    @property_RO
    def SE(self):
        '''Get the cell SouthEast of this (L{Geohash}).
        '''
        return self.S.E

    @property_RO
    def SW(self):
        '''Get the cell SouthWest of this (L{Geohash}).
        '''
        return self.S.W


def bounds(geohash):
    '''Returns the SW and NE lat-/longitude bounds of a geohash.

       @return: 4-Tuple (latS, lonW, latN, lonE) of bounds in (C{degrees}).

       @raise TypeError: The I{geohash} is not a L{Geohash}, C{LatLon}
                         or C{str}.

       @raise ValueError: Invalid or null I{geohash}.

       @example:

       >>> geohash.bounds('u120fxw')  #  52.20428467, 0.11810303,
                                      #  52.20565796, 0.11947632
       >>> geohash.decode('u120fxw')  # '52.205',    '0.1188'
    '''
    geohash = _2Geohash(geohash)
    if len(geohash) < 1:
        raise ValueError('%s invalid: %s' % ('geohash', geohash))

    latS, latN =  -90,  90
    lonW, lonE = -180, 180

    e = True
    for c in geohash:  # .lower():
        try:
            i = _DecodedBase32[c]
        except KeyError:
            raise ValueError('%s invalid: %s' % ('geohash', geohash))

        for m in (16, 8, 4, 2, 1):
            if e:  # longitude
                lon = favg(lonW, lonE)
                if i & m:
                    lonW = lon
                else:
                    lonE = lon
            else:  # latitude
                lat = favg(latS, latN)
                if i & m:
                    latS = lat
                else:
                    latN = lat
            e = not e

    return latS, lonW, latN, lonE


def decode(geohash):
    '''Decode a geohash to lat-/longitude of the (approximate
       centre of) geohash cell, to reasonable precision.

       @param geohash: To be decoded (L{Geohash}).

       @return: 2-Tuple (latStr, lonStr) in (C{str}).

       @raise TypeError: The I{geohash} is not a L{Geohash}, C{LatLon}
                         or C{str}.

       @raise ValueError: Invalid or null I{geohash}.

       @example:

       >>> geohash.decode('u120fxw')  # '52.205', '0.1188'
       >>> geohash.decode('sunny')  # '23.708', '42.473'  Saudi Arabia
       >>> geohash.decode('fur')  # '69.6', '-45.7'  Greenland
       >>> geohash.decode('reef')  # '-24.87', '162.95'  Coral Sea
       >>> geohash.decode('geek')  # '65.48', '-17.75'  Iceland
    '''
    s, w, n, e = bounds(geohash)

    # round to near centre without excessive precision
    # ⌊2-log10(Δ°)⌋ decimal places, strip trailing zeros
    lat = fStr(favg(n, s), prec=int(2 - log10(n - s)))
    lon = fStr(favg(e, w), prec=int(2 - log10(e - w)))

    return lat, lon  # strings


def decode_error(geohash):
    '''Return the relative lat-/longitude decoding errors for
       this geohash.

       @param geohash: To be decoded (L{Geohash}).

       @return: 2-Tuple (latErr, lonErr) in (C{degrees}).

       @raise TypeError: The I{geohash} is not a L{Geohash}, C{LatLon}
                         or C{str}.

       @raise ValueError: Invalid or null I{geohash}.

       @example:

       >>> geohash.decode_error('u120fxw')  # 0.00068665, 0.00068665
       >>> geohash.decode_error('fur')  # 0.703125, 0.703125
       >>> geohash.decode_error('fu')  # 2.8125, 5.625
       >>> geohash.decode_error('f')  # 22.5, 22.5
    '''
    s, w, n, e = bounds(geohash)
    return (n - s) * 0.5, (e - w) * 0.5


def distance1(geohash1, geohash2):
    '''Estimate the distance between two geohash (from the cell sizes).

       @param geohash1: First geohash (L{Geohash}).
       @param geohash2: Second geohash (L{Geohash}).

       @return: Approximate distance (C{meter}).

       @raise TypeError: A I{geohash*} is not a L{Geohash}, C{LatLon}
              or C{str}.

       @example:

       >>> geohash.distance1('u120fxwsh', 'u120fxws0')  # 15.239
    '''
    return _2Geohash(geohash1).distance1(geohash2)


def distance2(geohash1, geohash2, radius=R_M):
    '''Approximate the distance between two geohashes (with
       Pythagoras' theorem).

       @param geohash1: First geohash (L{Geohash}).
       @param geohash2: Second geohash (L{Geohash}).
       @keyword radius: Optional, mean earth radius (C{meter}) or C{None}.

       @return: Approximate distance (C{meter}, same units as I{radius}).

       @raise TypeError: A I{geohash*} is not a L{Geohash}, C{LatLon}
              or C{str}.

       @example:

       >>> geohash.distance2('u120fxwsh', 'u120fxws0')  # 19.0879
    '''
    return _2Geohash(geohash1).distance2(geohash2, radius=radius)


def distance3(geohash1, geohash2, radius=R_M):
    '''Compute the great-circle distance between two geohashes
       (using the Haversine formula).

       @param geohash1: First geohash (L{Geohash}).
       @param geohash2: Second geohash (L{Geohash}).
       @keyword radius: Optional, mean earth radius (C{meter}).

       @return: Great-circle distance (C{meter}, same units as I{radius}).

       @raise TypeError: A I{geohash*} is not a L{Geohash}, C{LatLon}
              or C{str}.

       @example:

       >>> geohash.distance3('u120fxwsh', 'u120fxws0')  # 11.6978
    '''
    return _2Geohash(geohash1).distance3(geohash2, radius=radius)


def encode(lat, lon, precision=None):
    '''Encode a lat-/longitude as a geohash, either to the specified
       or if not given, an automatically evaluated precision.

       @param lat: Latitude (C{degrees}).
       @param lon: Longitude (C{degrees}).
       @keyword precision: Optional, desired geohash length (C{int}).

       @return: The geohash (C{str}).

       @raise ValueError: Invalid I{lat}, I{lon} or I{precision}.

       @example:

       >>> geohash.encode(52.205, 0.119,   7)  # 'u120fxw'
       >>> geohash.encode(52.205, 0.119,  12)  # 'u120fxwshvkg'
       >>> geohash.encode(52.205, 0.1188, 12)  # 'u120fxws0jre'
       >>> geohash.encode(52.205, 0.1188)      # 'u120fxw'
       >>> geohash.encode(     0, 0)           # 's00000000000'
    '''
    lat, lon = _2fll(lat, lon)

    if not precision:
        # Infer precision by refining geohash until
        # it matches precision of supplied lat/lon.
        for prec in range(1, 13):
            gh = encode(lat, lon, prec)
            ll = map2(float, decode(gh))
            if abs(lat - ll[0]) < EPS and \
               abs(lon - ll[1]) < EPS:
                return gh
        prec = 12  # maximum
    else:
        try:
            prec = int(precision)
            if not 0 < prec < 13:
                raise ValueError
        except ValueError:
            raise ValueError('%s invalid: %r' % ('precision', precision))

    latS, latN =  -90,  90
    lonW, lonE = -180, 180

    b = i = 0
    e, gh = True, []

    while len(gh) < prec:
        i += i
        if e:  # bisect longitude
            m = favg(lonE, lonW)
            if lon < m:
                lonE = m
            else:
                lonW = m
                i += 1
        else:  # bisect latitude
            m = favg(latN, latS)
            if lat < m:
                latN = m
            else:
                latS = m
                i += 1
        e = not e

        b += 1
        if b == 5:
            # 5 bits gives a character:
            # append it and start over
            gh.append(_GeohashBase32[i])
            b = i = 0

    return ''.join(gh)


def neighbors(geohash):
    '''Return the L{Geohash}es for all 8 adjacent cells.

       @param geohash: Cell for which neighbors are required (L{Geohash} or C{str}).

       @return: Neighbors (dict(N=, NE=, E= ..., SW=) of L{Geohash}es).

       @raise TypeError: The I{geohash} is not a L{Geohash}, C{LatLon}
              or C{str}.

       @JSname: I{neighbours}.
    '''
    return _2Geohash(geohash).neighbors


def sizes(geohash):
    '''Return the lat- and longitudinal size of this L{Geohash} cell.

       @param geohash: Cell for which size are required (L{Geohash} or C{str}).

       @return: 2-Tuple (latHeight, lonWidth) in (C{meter}).

       @raise TypeError: The I{geohash} is not a L{Geohash}, C{LatLon}
              or C{str}.
    '''
    return _2Geohash(geohash).sizes

# **) MIT License
#
# Copyright (C) 2016-2018 -- mrJean1 at Gmail dot com
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
