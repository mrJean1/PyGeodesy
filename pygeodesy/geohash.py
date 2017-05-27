
# -*- coding: utf-8 -*-

'''
Class L{Geohash} and several functions to encode, decode and inspect
geohashes.

Transcribed from JavaScript originals by I{(C) Chris Veness 2011-2015}
and published under the same MIT Licence**.

More details at U{http://www.movable-type.co.uk/scripts/geohash.html}.
See also U{http://github.com/vinsci/geohash},
U{http://github.com/davetroy/geohash-js} and
U{http://pypi.python.org/pypi/pygeohash}.

@newfield example: Example, Examples
'''

from dms import parseDMS
from utils import EPS, R_M, favg, fStr, hsin3, map2

from math import cos, hypot, log10, radians

# all public contants, classes and functions
__all__ = ('Geohash',  # classes
           'bounds', 'decode', 'decode_error',  # functions
           'distance1', 'distance2', 'distance3',
           'encode', 'neighbors', 'sizes')
__version__ = '17.05.27'

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

try:
    _Str = str, basestring
except NameError:
    _Str = str

# Geohash-specific base32 map
_GeohashBase32 = '0123456789bcdefghjkmnpqrstuvwxyz'
# ... and the inverse map
_DecodedBase32 = dict((c, i) for i, c in enumerate(_GeohashBase32))
c = i = None
del c, i


def _2fll(lat, lon):
    '''(INTERNAL) Convert lat, lon to 2-tuple of floats.
    '''
    try:
        lat = parseDMS(lat, 'NS')
        if abs(lat) > 90:
            raise ValueError
    except ValueError:
        raise ValueError('%s invalid: %r' % ('lat', lat))

    try:
        lon = parseDMS(lon, 'EW')
        if abs(lon) > 180:
            raise ValueError
    except ValueError:
        raise ValueError('%s invalid: %r' % ('lon', lon))

    return lat, lon


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
        if not 0 < len(geohash) < 13:
            raise ValueError
        geostr = geohash.lower()
        for c in geostr:
            if c not in _DecodedBase32:
                raise ValueError
    except (AttributeError, ValueError, TypeError):
        raise ValueError('%s: %r[%s]' % (Geohash.__name__,
                          geohash, len(geohash)))
    return geostr


class Geohash(str):
    '''Geohash class, sub-class of str.
    '''
    _bounds = None
    _latlon = None

    _N  = None
    _E  = None
    _S  = None
    _W  = None
    _NE = None
    _NW = None
    _SE = None
    _SW = None

    # no str.__init__ in Python 3
    def __new__(cls, cll, precision=None):
        '''New Geohash.

           @param cll: Cell or location (L{Geohash}, I{LatLon} or string).

           @return: New L{Geohash}.
        '''
        if isinstance(cll, Geohash):
            self = str.__new__(cls, _2geostr('%s' % (cll,)))

        elif isinstance(cll, _Str):
            self = str.__new__(cls, _2geostr(cll))

        else:  # assume LatLon
            try:
                lat, lon = _2fll(cll.lat, cll.lon)
            except AttributeError:
                raise TypeError('%s: %r' % (Geohash.__name__, cll))
            cll = encode(lat, lon, precision)
            self = str.__new__(cls, cll)
            self._latlon = lat, lon

        return self

    def __repr__(self):
        return "%s('%s')" % (Geohash.__name__, self)  # str.__str__(self))

    def adjacent(self, direction):
        '''Determines adjacent cell in given compass direction.

           @param direction: Compass direction ('N', 'S', 'E' or 'W').

           @return: Geohash of adjacent cell (L{Geohash}).

           @raise ValueError: Invalid direction or geohash.
        '''
        # based on <http://github.com/davetroy/geohash-js>

        d = direction.upper()
        if not d or d not in _Neighbor:
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

    def bounds(self, LatLon):
        '''Returns SW/NE lat-/longitude bounds of this geohash cell.

           @param LatLon: LatLon class to use (I{LatLon}).

           @return: 2-Tuple (LatLonSW, LatLonNE) of I{LatLon}s for
                    the lower-left respectively upper-right corner.
        '''
        if not self._bounds:
            self._bounds = bounds(self)
        s, w, n, e = self._bounds
        return LatLon(s, w), LatLon(n, e)

    def distance1(self, other):
        '''Estimates the distance between this and an other geohash
           (from the cell sizes).

           @param other: The other geohash (L{Geohash}).

           @return: Approximate distance (meter).

           @raise TypeError: The other is not a L{Geohash}, I{LatLon} or str.
        '''
        other = _2Geohash(other)

        n = min(len(self), len(other), len(_Sizes))
        for n in range(n):
            if self[n] != other[n]:
                break
        return float(_Sizes[n][2])

    def distance2(self, other, radius=R_M):
        '''Approximates the distance between this and an other geohash
           (with Pythagoras' theorem).

           @param other: The other geohash (L{Geohash}).
           @keyword radius: Optional earth radius (meter) or None.

           @return: Approximate distance (meter, same units as radius).

           @raise TypeError: The other is not a L{Geohash}, I{LatLon} or str.
        '''
        other = _2Geohash(other)

        a1, b1 = self.latlon
        a2, b2 = other.latlon

        x = radians(a2 - a1) * cos(radians(favg(b2, b1)))
        y = radians(b2 - b1)

        if radius:
            d = hypot(x, y) * float(radius)
        else:
            d = x * x + y * y
        return d

    def distance3(self, other, radius=R_M):
        '''Computes the great-circle distance between this and
           an other geohash (using the Haversine formula).

           @param other: The other geohash (L{Geohash}).
           @keyword radius: Optional earth radius (meter).

           @return: Great-circle distance (meter, same units as radius).

           @raise TypeError: The other is not a L{Geohash}, I{LatLon} or str.
        '''
        other = _2Geohash(other)

        a1, b1 = self.latlon
        a2, b2 = other.latlon

        r, _, _ = hsin3(radians(a2), radians(a1), radians(b2 - b1))
        return r * float(radius)

    @property
    def latlon(self):
        '''Gets the lat-/longitude of (the approximate center of)
           this geohash as 2-Tuple (lat, lon) in degrees.

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

    @property
    def neighbors(self):
        '''Gets all 8 adjacent cells as a dict(N=, NE=, E= ..., SW=)
           of L{Geohash}es.

           B{JSname:} I{neighbours}.
        '''
        return dict((d, getattr(self, d)) for d in ('N', 'NE', 'E', 'SE',
                                                    'S', 'SW', 'W', 'NW'))

    @property
    def sizes(self):
        '''Gets the lat- and longitudinal size of this cell as a 2-tuple
           (latHeight, lonWidth) in meter.
        '''
        n = min(len(self) or 1, len(_Sizes)) - 1
        return tuple(map(float, _Sizes[n][:2]))

    def toLatLon(self, LatLon, **kwds):
        '''Returns (the approximate center of) this geohash cell
           as an instance of the supplied I{LatLon} class.

           @param LatLon: Class to use (I{LatLon}).
           @keyword kwds: Optional keyword arguments for I{LatLon}.

           @return: This geohash location (I{LatLon}).

           @example:

           >>> from sphericalTrigonometry import LatLon
           >>> ll = Geohash('u120fxw').toLatLon(LatLon)
           >>> print(repr(ll))  # LatLon(52°12′17.9″N, 000°07′07.64″E)
           >>> print(ll)  # 52.204971°N, 000.11879°E
        '''
        return LatLon(*self.latlon, **kwds)

    @property
    def N(self):
        '''Gets the cell North of this (L{Geohash}).
        '''
        if self._N is None:
            self._N = self.adjacent('N')
        return self._N

    @property
    def S(self):
        '''Gets the cell South of this (L{Geohash}).
        '''
        if self._S is None:
            self._S = self.adjacent('S')
        return self._S

    @property
    def E(self):
        '''Gets the cell East of this (L{Geohash}).
        '''
        if self._E is None:
            self._E = self.adjacent('E')
        return self._E

    @property
    def W(self):
        '''Gets the cell West of this (L{Geohash}).
        '''
        if self._W is None:
            self._W = self.adjacent('W')
        return self._W

    @property
    def NE(self):
        '''Gets the cell NorthEast of this (L{Geohash}).
        '''
        if self._NE is None:
            self._NE = self.N.E
        return self._NE

    @property
    def NW(self):
        '''Gets the cell NorthWest of this (L{Geohash}).
        '''
        if self._NW is None:
            self._NW = self.N.W
        return self._NW

    @property
    def SE(self):
        '''Gets the cell SouthEast of this (L{Geohash}).
        '''
        if self._SE is None:
            self._SE = self.S.E
        return self._SE

    @property
    def SW(self):
        '''Gets the cell SouthWest of this (L{Geohash}).
        '''
        if self._SW is None:
            self._SW = self.S.W
        return self._SW


def bounds(geohash):
    '''Returns SW and NE lat-/longitude bounds of a geohash.

       @return: 4-Tuple (latS, lonW, latN, lonE) in (degrees).

       @raise TypeError: The geohash is not a L{Geohash}, I{LatLon} or str.

       @raise ValueError: Invalid or null geohash.

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
    '''Decodes a geohash to lat-/longitude of the (approximate
       centre of) geohash cell, to reasonable precision.

       @param geohash: To be decoded (L{Geohash}).

       @return: 2-Tuple (latStr, lonStr) in (strings).

       @raise TypeError: The geohash is not a L{Geohash}, I{LatLon} or str.

       @raise ValueError: Invalid or null geohash.

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
    '''Returns relative lat-/longitude decoding errors for
       this geohash.

       @param geohash: To be decoded (L{Geohash}).

       @return: 2-Tuple (latErr, lonErr) in (degrees).

       @raise TypeError: The geohash is not a L{Geohash}, I{LatLon} or str.

       @raise ValueError: Invalid or null geohash.

       @example:

       >>> geohash.decode_error('u120fxw')  # 0.00068665, 0.00068665
       >>> geohash.decode_error('fur')  # 0.703125, 0.703125
       >>> geohash.decode_error('fu')  # 2.8125, 5.625
       >>> geohash.decode_error('f')  # 22.5, 22.5
    '''
    s, w, n, e = bounds(geohash)
    return (n - s) * 0.5, (e - w) * 0.5


def distance1(geohash1, geohash2):
    '''Estimates the distance between two geohash (from the cell sizes).

       @param geohash1: First geohash (L{Geohash}).
       @param geohash2: Second geohash (L{Geohash}).

       @return: Approximate distance (meter).

       @raise TypeError: A geohash is not a L{Geohash}, I{LatLon} or str.

       @example:

       >>> geohash.distance1('u120fxwsh', 'u120fxws0')  # 15.239
    '''
    return _2Geohash(geohash1).distance1(geohash2)


def distance2(geohash1, geohash2, radius=R_M):
    '''Approximates the distance between two geohashes (with
       Pythagoras' theorem).

       @param geohash1: First geohash (L{Geohash}).
       @param geohash2: Second geohash (L{Geohash}).
       @keyword radius: Optional earth radius (meter) or None.

       @return: Approximate distance (meter, same units as radius).

       @raise TypeError: A geohash is not a L{Geohash}, I{LatLon} or str.

       @example:

       >>> geohash.distance2('u120fxwsh', 'u120fxws0')  # 19.0879
    '''
    return _2Geohash(geohash1).distance2(geohash2, radius=radius)


def distance3(geohash1, geohash2, radius=R_M):
    '''Computes the great-circle distance between two geohashes
       (using the Haversine formula).

       @param geohash1: First geohash (L{Geohash}).
       @param geohash2: Second geohash (L{Geohash}).
       @keyword radius: Optional earth radius (meter).

       @return: Great-circle distance (meter, same units as radius).

       @raise TypeError: A geohash is not a L{Geohash}, I{LatLon} or str.

       @example:

       >>> geohash.distance3('u120fxwsh', 'u120fxws0')  # 11.6978
    '''
    return _2Geohash(geohash1).distance3(geohash2, radius=radius)


def encode(lat, lon, precision=None):
    '''Encodes lat-/longitude to geohash, either to the
       specified or an automatically evaluated precision.

       @param lat: Latitude in degrees (scalar).
       @param lon: Longitude in degrees (scalar).
       @keyword precision: Desired geohash length (integer).

       @return: The geohash (string).

       @raise ValueError: Invalid lat, lon or precision.

       @example:

       >>> geohash.encode(52.205, 0.119,   7)  # 'u120fxw'
       >>> geohash.encode(52.205, 0.119,  12)  # 'u120fxwshvkg'
       >>> geohash.encode(52.205, 0.1188, 12)  # 'u120fxws0jre'
       >>> geohash.encode(52.205, 0.1188)      # 'u120fxw'
       >>> geohash.encode(     0, 0)           # 's00000000000'
    '''
    lat, lon = _2fll(lat, lon)

    if precision is None:
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
    '''Returns the L{Geohash}es for all 8 adjacent cells.

       @param geohash: Cell for which neighbors are required (L{Geohash} or str).

       @return: Neighbors (dict(N=, NE=, E= ..., SW=) of L{Geohash}es).

       @raise TypeError: The geohash is not a L{Geohash}, I{LatLon} or str.

       @JSname: I{neighbours}.
    '''
    return _2Geohash(geohash).neighbors


def sizes(geohash):
    '''Returns the lat- and longitudinal size of this L{Geohash} cell.

       @param geohash: Cell for which size are required (L{Geohash} or str).

       @return: 2-Tuple (latHeight, lonWidth) in (meter).

       @raise TypeError: The geohash is not a L{Geohash}, I{LatLon} or str.
    '''
    return _2Geohash(geohash).sizes

# **) MIT License
#
# Copyright (C) 2016-2017 -- mrJean1 at Gmail dot com
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
