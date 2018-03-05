
# -*- coding: utf-8 -*-

u'''Functions to parse and format bearing, compass, lat- and longitudes
in various forms of degrees, minutes and seconds.

After I{(C) Chris Veness 2011-2015} published under the same MIT Licence**, see
U{Latitude/Longitude<http://www.movable-type.co.uk/scripts/latlong.html>} and
U{Vector-based geodesy<http://www.movable-type.co.uk/scripts/latlong-vectors.html>}.

@newfield example: Example, Examples
'''

from fmath import fStr, fStrzs, isint
from utils import degrees360

from math import atan2, copysign, radians
try:
    from string import letters as LETTERS
except ImportError:  # Python 3+
    from string import ascii_letters as LETTERS

# all public contants, classes and functions
__all__ = ('F_D', 'F_DM', 'F_DMS',  # forms
           'F_DEG', 'F_MIN', 'F_SEC', 'F_RAD',
           'S_DEG', 'S_MIN', 'S_SEC', 'S_RAD', 'S_SEP',  # symbols
           'RangeError',  # classes
           'bearingDMS',  # functions
           'clipDMS', 'compassAngle', 'compassDMS', 'compassPoint',
           'latDMS', 'lonDMS', 'normDMS',
           'parseDMS', 'parseDMS2', 'parse3llh', 'precision',
           'rangerrors', 'toDMS')
__version__ = '18.03.04'

F_D   = 'd'    #: Format degrees as deg° (string).
F_DM  = 'dm'   #: Format degrees as deg°min′ (string).
F_DMS = 'dms'  #: Format degrees as deg°min′sec″ (string).
F_DEG = 'deg'  #: Format degrees as [D]DD without symbol (string).
F_MIN = 'min'  #: Format degrees as [D]DDMM without symbols (string).
F_SEC = 'sec'  #: Format degrees as [D]DDMMSS without symbols (string).
F_RAD = 'rad'  #: Convert degrees to radians and format as RR.r (string).

S_DEG = '°'  #: Degrees ° symbol (string).
S_MIN = '′'  #: Minutes ′ symbol (string).
S_SEC = '″'  #: Seconds ″ symbol (string).
S_RAD = ''   #: Radians symbol (string).
S_SEP = ''   #: Separator between deg, min and sec (string).

_F_prec = {F_D:   6, F_DM:  4, F_DMS: 2,  #: (INTERNAL) default precs.
           F_DEG: 6, F_MIN: 4, F_SEC: 2, F_RAD: 5}

_S_norm = {'^': S_DEG, '˚': S_DEG,  #: (INTERNAL) normalized DMS.
           "'": S_MIN, '’': S_MIN, '′': S_MIN,
           '"': S_SEC, '″': S_SEC, '”': S_SEC}
_S_ALL  = (S_DEG, S_MIN, S_SEC) + tuple(_S_norm.keys())  #: (INTERNAL) alternates.

_rangerrors = True


class RangeError(ValueError):
    '''Error raised for lat- or longitude values outside the I{clip},
       I{clipLat}, I{clipLon} or I{limit} range in function L{clipDMS},
       L{parse3llh}, L{parseDMS} or L{parseDMS2}.
    '''
    pass


def _toDMS(deg, form, prec, sep, ddd):
    '''(INTERNAL) Convert degrees to string, without sign or suffix.
    '''
    try:
        d = abs(float(deg))
    except ValueError:
        raise ValueError('%s invalid: %r' % ('deg', deg))

    if prec is None:
        z = p = _F_prec.get(form, 6)
    else:
        z = int(prec)
        p = abs(z)
    w = p + (1 if p else 0)

    f = form.lower()
    if f in (F_DEG, F_MIN, F_SEC):
        s_deg = s_min = s_sec = ''  # no symbols
    else:
        s_deg, s_min, s_sec = S_DEG, S_MIN, S_SEC

    if f in (F_D, F_DEG, 'degrees'):  # deg°, degrees
        t = '%0*.*f' % (ddd+w,p,d)
        s = s_deg

    elif f in (F_RAD, 'radians'):
        t = '%.*f' % (p,radians(d))
        s = S_RAD

    elif f in (F_DM, F_MIN, 'deg+min'):
        d, m = divmod(d * 60, 60)
        t = "%0*d%s%s%0*.*f" % (ddd,int(d),s_deg, sep, w+2,p,m)
        s = s_min

    else:  # F_DMS, F_SEC, 'deg+min+sec'
        d, s = divmod(d * 3600, 3600)
        m, s = divmod(s, 60)
        t = "%0*d%s%s%02d%s%s%0*.*f" % (ddd,int(d),s_deg, sep,
                                            int(m),s_min, sep,
                                          w+2,p,s)
        s = s_sec

    if z > 1:
        t = fStrzs(t)
    return t + s


def bearingDMS(bearing, form=F_D, prec=None, sep=S_SEP):
    '''Convert bearing to a string.

       @param bearing: Bearing from North (compass degrees).
       @keyword form: Optional format, L{F_D}, L{F_DM}, L{F_DMS},
                      L{F_DEG}, L{F_MIN}, L{F_SEC} or L{F_RAD} for
                      deg°, deg°min′, deg°min′sec″, [D]DD, [D]DDMM,
                      [D]DDMMSS or radians (string).
       @keyword prec: Optional number of decimal digits (0..9 or
                      None for default).  Trailing zero decimals
                      are stripped for prec values of 1 and above,
                      but kept for negative prec values.
       @keyword sep: Optional separator (string).

       @return: Compass degrees per the specified form (string).

       @JSname: I{toBrng}.
    '''
    return _toDMS(bearing % 360, form, prec, sep, 1)


def clipDMS(deg, limit):
    '''Clip a lat- or longitude to the given range.

       @param deg: Unclipped lat- or longitude (degrees).
       @param limit: Valid I{-limit..+limit} range (degrees).

       @return: Clipped value (degrees).

       @raise RangeError: If I{deg} beyond I{limit} and L{rangerrors}
                          set to True.
    '''
    if limit > 0:
        c = min(limit, max(-limit, deg))
        if _rangerrors and deg != c:
            raise RangeError('%s beyond %s degrees' % (fStr(deg, prec=6),
                             fStr(copysign(limit, deg), prec=3, ints=True)))
        deg = c
    return deg


def compassAngle(lat0, lon0, lat1, lon1):
    '''Return the angle from North for the direction vector
       M{(lon1 - lon0, lat1 - lat0)} between two points.

       Suitable only for short, non-near-polar vectors up to a few
       hundred Km or Miles.  Use I{LatLon} methods I{initialBearingTo}
       or I{forward azimuth} for larger distances.

       @param lat0: From latitude (degrees).
       @param lon0: From longitude (degrees).
       @param lat1: To latitude (degrees).
       @param lon1: To longitude (degrees).

       @return: Angle from North (degrees360).

       @note: Courtesy Martin Schultz.

       @see: U{Local, Flat Earth<http://www.edwilliams.org/avform.htm#flat>}.
    '''
    return degrees360(atan2(lon1 - lon0, lat1 - lat0))


def compassDMS(bearing, form=F_D, prec=None, sep=S_SEP):
    '''Convert bearing to a string suffixed with compass point.

       @param bearing: Bearing from North (compass degrees).
       @keyword form: Optional format, L{F_D}, L{F_DM}, L{F_DMS},
                      L{F_DEG}, L{F_MIN}, L{F_SEC} or L{F_RAD} for
                      deg°, deg°min′, deg°min′sec″, [D]DD, [D]DDMM,
                      [D]DDMMSS or radians (string).
       @keyword prec: Optional number of decimal digits (0..9 or
                      None for default).  Trailing zero decimals
                      are stripped for prec values of 1 and above,
                      but kept for negative prec values.
       @keyword sep: Optional separator (string).

       @return: Compass degrees and point in the specified form (string).
    '''
    t = bearingDMS(bearing, form, prec, sep), compassPoint(bearing)
    return sep.join(t)


_COMPASS = ('N', 'NbE', 'NNE', 'NEbN', 'NE', 'NEbE', 'ENE', 'EbN',
            'E', 'EbS', 'ESE', 'SEbE', 'SE', 'SEbS', 'SSE', 'SbE',
            'S', 'SbW', 'SSW', 'SWbS', 'SW', 'SWbW', 'WSW', 'WbS',
            'W', 'WbN', 'WNW', 'NWbW', 'NW', 'NWbN', 'NNW', 'NbW')  #: (INTERNAL) points

_M_X = {1: (4, 8), 2: (8, 4), 3: (16, 2), 4: (32, 1)}  #: (INTERNAL) precs


def compassPoint(bearing, prec=3):
    '''Convert bearing to a compass point.

       @param bearing: Bearing from North (compass degrees).
       @keyword prec: Optional precision (1 for cardinal or basic winds,
                      2 for intercardinal or ordinal or principal winds,
                      3 for secondary-intercardinal or half-winds or
                      4 for quarter-winds).

       @return: Compass point (1-, 2-, 3- or 4-letter string).

       @raise ValueError: Invalid I{prec}.

       @see: U{Compass rose<http://wikipedia.org/wiki/Compass_rose>}

       @example:

       >>> p = compassPoint(24)     # 'NNE'
       >>> p = compassPoint(24, 1)  # 'N'
       >>> p = compassPoint(24, 2)  # 'NE'
       >>> p = compassPoint(18, 3)  # 'NNE'
       >>> p = compassPoint(12, 4)  # 'NbE'
       >>> p = compassPoint(30, 4)  # 'NEbN'
    '''
    try:  # m = 2 << prec; q = 32 // m
        m, x = _M_X[prec]
    except KeyError:
        raise ValueError('%s invalid: %r' % ('prec', prec))

    q = int(round((bearing % 360) * m / 360.0)) % m
    return _COMPASS[q * x]


def latDMS(deg, form=F_DMS, prec=2, sep=S_SEP):
    '''Convert latitude to a string suffixed with N or S.

       @param deg: Latitude to be formatted (degrees).
       @keyword form: Optional format, L{F_D}, L{F_DM}, L{F_DMS},
                      L{F_DEG}, L{F_MIN}, L{F_SEC} or L{F_RAD} for
                      deg°, deg°min′, deg°min′sec″, DD, DDMM, DDMMSS
                      or radians (string).
       @keyword prec: Optional number of decimal digits (0..9 or
                      None for default).  Trailing zero decimals
                      are stripped for prec values of 1 and above,
                      but kept for negative prec values.
       @keyword sep: Optional separator (string).

       @return: Degrees in the specified form (string).

       @JSname: I{toLat}.
    '''
    t = _toDMS(deg, form, prec, sep, 2), ('S' if deg < 0 else 'N')
    return sep.join(t)


def lonDMS(deg, form=F_DMS, prec=2, sep=S_SEP):
    '''Convert longitude to a string suffixed with E or W.

       @param deg: Longitude to be formatted (degrees).
       @keyword form: Optional format, L{F_D}, L{F_DM}, L{F_DMS},
                      L{F_DEG}, L{F_MIN}, L{F_SEC} or L{F_RAD} for
                      deg°, deg°min′, deg°min′sec″, DDD, DDDMM,
                      DDDMMSS or radians (string).
       @keyword prec: Optional number of decimal digits (0..9 or
                      None for default).  Trailing zero decimals
                      are stripped for prec values of 1 and above,
                      but kept for negative prec values.
       @keyword sep: Optional separator (string).

       @return: Degrees in the specified form (string).

       @JSname: I{toLon}.
    '''
    t = _toDMS(deg, form, prec, sep, 3), ('W' if deg < 0 else 'E')
    return sep.join(t)


def normDMS(strDMS, norm=''):
    '''Normalize all degree ˚, minute ' and second " symbols in a
       string to the default symbols %s, %s and %s.

       @param strDMS: DMS (string).
       @keyword norm: Optional replacement symbol, default symbol
                      otherwise (string).

       @return: Normalized DMS (string).
    '''
    if norm:
        for s in _S_ALL:
            strDMS = strDMS.replace(s, norm)
        strDMS = strDMS.rstrip(norm)
    else:
        for s, S in _S_norm.items():
            strDMS = strDMS.replace(s, S)
    return strDMS


if __debug__:  # no __doc__ at -O and -OO
    normDMS.__doc__  %= (S_DEG, S_MIN, S_SEC)


def parse3llh(strll, height=0, sep=',', clipLat=90, clipLon=180):
    '''Parse a string representing lat-, longitude and height point.

       The lat- and longitude value must be separated by a separator
       character.  If height is present it must follow, separated by
       another separator.

       The lat- and longitude values may be swapped, provided at least
       one ends with the proper compass point.

       @param strll: Latitude, longitude[, height] (string).
       @keyword height: Optional, default for missing height (meter).
       @keyword sep: Optional separator (string).
       @keyword clipLat: Keep latitude in I{-clipLat..+clipLat} (degrees).
       @keyword clipLon: Keep longitude in I{-clipLon..+clipLon} range (degrees).

       @return: 3-Tuple (lat, lon, height) as (scalars).

       @raise RangeError: Lat- or longitude value of I{strll} outside
                          valid range and I{rangerrrors} set to True.

       @raise ValueError: Invalid I{strll}.

       @see: Functions L{parseDMS} and L{parseDMS2} for more details
             on the forms and symbols accepted.

       @example:

       >>> parse3llh('000°00′05.31″W, 51° 28′ 40.12″ N')
       (51.4778°N, 000.0015°W, 0)
    '''
    try:
        ll = strll.strip().split(sep)
        if len(ll) > 2:  # XXX interpret height unit
            h = float(ll.pop(2).strip().rstrip(LETTERS).rstrip())
        else:
            h = height
        if len(ll) != 2:
            raise ValueError
    except (AttributeError, TypeError, ValueError):
        return ValueError('parsing %r failed' % (strll,))

    a, b = [_.strip() for _ in ll]
    if a[-1:] in 'EW' or b[-1:] in 'NS':
        a, b = b, a
    a, b = parseDMS2(a, b, clipLat=clipLat, clipLon=clipLon)
    return a, b, h


def parseDMS(strDMS, suffix='NSEW', sep=S_SEP, clip=0):
    '''Parse a string representing deg°min′sec″ to degrees.

       This is very flexible on formats, allowing signed decimal
       degrees, degrees and minutes or degrees minutes and seconds
       optionally suffixed by compass direction NSEW.

       A variety of symbols, separators and suffixes are accepted,
       for example 3° 37′ 09″W.  Minutes and seconds may be omitted.

       @param strDMS: Degrees in any of several forms (string or scalar).
       @keyword suffix: Optional, valid compass directions (NSEW).
       @keyword sep: Optional separator between deg°, min′ and sec″ ('').
       @keyword clip: Optionally, limit value to -clip..+clip (degrees).

       @return: Degrees (float).

       @raise RangeError: Value of I{strDMS} outside the valid range
                          and I{rangerrrors} set to True.

       @raise ValueError: Invalid I{strDMS}.

       @see: Function L{parse3llh} to parse a string with lat-,
             longitude and height values.
    '''
    try:  # signed decimal degrees without NSEW
        d = float(strDMS)

    except (TypeError, ValueError):
        try:
            strDMS = strDMS.strip()

            t = strDMS.lstrip('-+').rstrip(suffix.upper())
            if sep:
                t = t.replace(sep, ' ')
                for s in _S_ALL:
                    t = t.replace(s, '')
            else:
                for s in _S_ALL:
                    t = t.replace(s, ' ')
            t = list(map(float, t.strip().split())) + [0, 0]
            d = t[0] + (t[1] + t[2] / 60.0) / 60.0
            if strDMS[:1] == '-' or strDMS[-1:] in 'SW':
                d = -d

        except (IndexError, ValueError):
            raise ValueError('parsing %r failed' % (strDMS,))

    return clipDMS(d, clip)


def parseDMS2(strLat, strLon, sep=S_SEP, clipLat=90, clipLon=180):
    '''Parse lat- and longitude representions.

       @param strLat: Latitude in any of several forms (string or scalar).
       @param strLon: Longitude in any of several forms (string or scalar).
       @keyword sep: Optional separator between deg°, min′ and sec″ ('').
       @keyword clipLat: Keep latitude in I{-clipLat..+clipLat} range (degrees).
       @keyword clipLon: Keep longitude in I{-clipLon..+clipLon} range (degrees).

       @return: 2-Tuple (lat, lon) in (degrees, degrees).

       @raise RangeError: Value of I{strLat} or I{strLon} outside the
                          valid range and I{rangerrrors} set to True.

       @raise ValueError: Invalid I{strLat} or I{strLon}.

       @see: Function L{parse3llh} to parse a string with lat-,
             longitude and height values and function L{parseDMS}
             to parse individual lat- or longitudes.
    '''
    return (parseDMS(strLat, suffix='NS', sep=sep, clip=clipLat),
            parseDMS(strLon, suffix='EW', sep=sep, clip=clipLon))


def precision(form, prec=None):
    '''Set the default precison for a given F_ form.

       @param form: L{F_D}, L{F_DM}, L{F_DMS}, L{F_DEG}, L{F_MIN},
                    L{F_SEC} or L{F_RAD} (string).
       @keyword prec: Optional number of decimal digits (0..9 or
                      None for default).  Trailing zero decimals
                      are stripped for I{prec} values of 1 and
                      above, but kept for negative I{prec} values.

       @return: Previous precision (int).

       @raise ValueError: If I{form} or I{prec} or beyond valid range.
    '''
    try:
        p = _F_prec[form]
    except KeyError:
        raise ValueError('%s invalid: %s' % ('form', form))
    if isint(prec):
        if not -10 < prec < 10:
            raise ValueError('%s invalid: %s' % ('prec', prec))
        _F_prec[form] = prec
    return p


def rangerrors(raiser=None):
    '''Gert/set raising of range errors.

       @keyword raiser: Choose True to raise or False to not raise
                        L{RangeError} exceptions.  Use None to leave
                        the setting unchanged.

       @return: Previous setting (bool).

       @note: Out-of-range lat- and longitude values are always
              clipped to the nearest range limit.
    '''
    global _rangerrors
    t = _rangerrors
    if raiser in (True, False):
        _rangerrors = raiser
    return t


def toDMS(deg, form=F_DMS, prec=2, sep=S_SEP, ddd=2, neg='-', pos=''):
    '''Convert signed degrees to string, without suffix.

       @param deg: Degrees to be formatted (scalar).
       @keyword form: Optional format, L{F_D}, L{F_DM}, L{F_DMS},
                      L{F_DEG}, L{F_MIN}, L{F_SEC} or L{F_RAD} for
                      deg°, deg°min′, deg°min′sec″, [D]DD, [D]DDMM,
                      [D]DDMMSS or radians (string).
       @keyword prec: Optional number of decimal digits (0..9 or
                      None for default).  Trailing zero decimals
                      are stripped for prec values of 1 and above,
                      but kept for negative prec values.
       @keyword sep: Optional separator (string).
       @keyword ddd: Optional number of digits for deg° (2 or 3).
       @keyword neg: Optional sign for negative degrees ('-').
       @keyword pos: Optional sign for positive degrees ('').

       @return: Degrees in the specified form (string).
    '''
    t = _toDMS(deg, form, prec, sep, ddd)
    s = neg if deg < 0 else pos
    return s + t

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
