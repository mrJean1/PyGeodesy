
# -*- coding: utf-8 -*-

'''Functions to parse and format bearing, lat- and longitudes in various
forms of degrees, minutes and seconds.

After I{(C) Chris Veness 2011-2015} published under the same MIT Licence**,
see U{http://www.movable-type.co.uk/scripts/latlong.html}
and U{http://www.movable-type.co.uk/scripts/latlong-vectors.html}.

@newfield example: Example, Examples
'''

from utils import fStrzs

from math import radians
try:
    from string import letters as LETTERS
except ImportError:  # Python 3+
    from string import ascii_letters as LETTERS

# all public contants, classes and functions
__all__ = ('F_D', 'F_DM', 'F_DMS', 'F_DEG', 'F_RAD',  # forms
           'S_DEG', 'S_MIN', 'S_SEC', 'S_RAD', 'S_SEP',  # symbols
           'bearingDMS', 'compassDMS', 'compassPoint',  # functions
           'latDMS', 'lonDMS', 'normDMS',
           'parseDMS', 'parse3llh', 'precision', 'toDMS')
__version__ = '17.05.29'

F_D   = 'd'    #: Format degrees as deg° (string).
F_DM  = 'dm'   #: Format degrees as deg°min′ (string).
F_DMS = 'dms'  #: Format degrees as deg°min′sec″ (string).
F_DEG = 'deg'  #: Format degrees as degrees without symbol (string).
F_RAD = 'rad'  #: Convert degrees to radians and format (string).

S_DEG = '°'  #: Degrees ° symbol (string).
S_MIN = '′'  #: Minutes ′ symbol (string).
S_SEC = '″'  #: Seconds ″ symbol (string).
S_RAD = ''   #: Radians symbol (string).
S_SEP = ''   #: Separator between deg, min and sec (string).

_F_prec = {F_D: 6, F_DM: 4, F_DMS: 2, F_DEG: 6, F_RAD: 5}  #: (INTERNAL) default precs.

_S_norm = {'^': S_DEG, '˚': S_DEG,  #: (INTERNAL) normalized DMS.
           "'": S_MIN, '’': S_MIN, '′': S_MIN,
           '"': S_SEC, '″': S_SEC, '”': S_SEC}
_S_ALL  = (S_DEG, S_MIN, S_SEC) + tuple(_S_norm.keys())  #: (INTERNAL) alternates.


def _toDMS(deg, form, prec, ddd):
    '''(INTERNAL) Converts degrees to string, without sign or suffix.
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
    if f in (F_D, F_DEG, 'degrees'):  # deg°, degrees
        t = '%0*.*f' % (ddd+w,p,d)
        s = S_DEG if f == F_D else ''

    elif f in (F_RAD, 'radians'):
        t = '%.*f' % (p, radians(d))
        s = S_RAD

    elif f in (F_DM, 'deg+min'):
        d, m = divmod(d * 60, 60)
        t = "%0*d%s%s%0*.*f" % (ddd,int(d),S_DEG, S_SEP, w+2,p,m)
        s = S_MIN

    else:  # F_DMS, 'deg+min+sec'
        d, s = divmod(d * 3600, 3600)
        m, s = divmod(s, 60)
        t = "%0*d%s%s%02d%s%s%0*.*f" % (ddd,int(d),S_DEG, S_SEP,
                                            int(m),S_MIN, S_SEP,
                                          w+2,p,s)
        s = S_SEC

    if z > 1:
        t = fStrzs(t)
    return t + s


def bearingDMS(bearing, form=F_D, prec=None):
    '''Converts bearing to string.

       @param bearing: Bearing from North (compass degrees).
       @keyword form: Use F_D, F_DM, F_DMS, F_DEG or F_RAD for deg°,
                      deg°min′, deg°min′sec″ or degrees or radians.
       @keyword prec: Optional, number of decimal digits (0..9 or
                      None for default).  Trailing zero decimals
                      are stripped for prec values of 1 and above,
                      but kept for negative prec values.

       @return: Compass degrees per the specified form (string).

       @JSname: I{toBrng}.
    '''
    return _toDMS(bearing % 360, form, prec, 1)


_COMPASS = ('N', 'NNE', 'NE', 'ENE',
            'E', 'ESE', 'SE', 'SSE',
            'S', 'SSW', 'SW', 'WSW',
            'W', 'WNW', 'NW', 'NNW')  #: (INTERNAL) points

_M_X = {1: (4, 4), 2: (8, 2), 3: (16, 1)}  #: (INTERNAL) precs


def compassDMS(bearing, form=F_D, prec=None):
    '''Converts bearing to string suffixed with compass point.

       @param bearing: Bearing from North (compass degrees).
       @keyword form: Use F_D, F_DM, F_DMS, F_DEG or F_RAD for deg°,
                      deg°min′, deg°min′sec″ or degrees or radians.
       @keyword prec: Optional, number of decimal digits (0..9 or
                      None for default).  Trailing zero decimals
                      are stripped for prec values of 1 and above,
                      but kept for negative prec values.

       @return: Compass degrees and point per the specified form (string).
    '''
    t = _toDMS(bearing % 360, form, prec, 1), compassPoint(bearing)
    return S_SEP.join(t)


def compassPoint(bearing, prec=3):
    '''Converts bearing to compass point.

       @param bearing: Bearing from North (compass degrees).
       @keyword prec: Optional precision (1 for cardinal, 2 for
                      intercardinal or 3 for secondary-intercardinal).

       @return: Compass point (1-, 2- or 3-letter string).

       @raise ValueError: Invalid prec.

       @example:

       >>> p = compass(24)     # 'NNE'
       >>> p = compass(24, 1)  # 'N'
    '''
    try:
        m, x = _M_X[prec]
    except KeyError:
        raise ValueError('%s invalid: %r' % ('prec', prec))

    q = int(round((bearing % 360) * m / 360.0)) % m
    return _COMPASS[q * x]


def latDMS(deg, form=F_DMS, prec=2):
    '''Converts latitude to string suffixed with N or S.

       @param deg: Latitude to be formatted (degrees).
       @keyword form: Use F_D, F_DM, F_DMS, F_DEG or F_RAD for deg°,
                      deg°min′, deg°min′sec″ or degrees or radians.
       @keyword prec: Optional, number of decimal digits (0..9 or
                      None for default).  Trailing zero decimals
                      are stripped for prec values of 1 and above,
                      but kept for negative prec values.

       @return: Degrees per the specified form (string).

       @JSname: I{toLat}.
    '''
    return _toDMS(deg, form, prec, 2) + ('S' if deg < 0 else 'N')


def lonDMS(deg, form=F_DMS, prec=2):
    '''Converts longitude to string suffixed with E or W.

       @param deg: Longitude to be formatted (degrees).
       @keyword form: Use F_D, F_DM, F_DMS, F_DEG or F_RAD for deg°,
                      deg°min′, deg°min′sec″ or degrees or radians.
       @keyword prec: Optional, number of decimal digits (0..9 or
                      None for default).  Trailing zero decimals
                      are stripped for prec values of 1 and above,
                      but kept for negative prec values.

       @return: Degrees per the specified form (string).

       @JSname: I{toLon}.
    '''
    return _toDMS(deg, form, prec, 3) + ('W' if deg < 0 else 'E')


toLon = lonDMS  # XXX original name


def normDMS(strDMS, norm=''):
    '''Normalizes all degrees ˚, minutes ' and seconds " symbols
       in a string to the defaults %s, %s and %s.

       @param strDMS: DMS (string).
       @keyword norm: Replacement, default symbol otherwise (string).

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


def parse3llh(strll, height=0, sep=','):
    '''Parses a string representing lat-, longitude and height point.

       The lat- and longitude value must be separated by a separator
       character.  If height is present it must follow, separated by
       another separator.

       The lat- and longitude values may be swapped, provided at least
       one ends with the proper compass direction.

       See function parseDMS for more details on the forms and
       symbols accepted.

       @param strll: Lat, lon[, height] (string).
       @keyword height: Default for missing height (meter).
       @keyword sep: Optional, separator (string).

       @return: 3-Tuple (lat, lon, height) as (scalars).

       @raise ValueError: Invalid strll.

       @example:

       >>> t = parse3llh('000°00′05.31″W, 51° 28′ 40.12″ N')  # (51.4778°N, 000.0015°W, 0)
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
    return parseDMS(a, suffix='NS'), parseDMS(b, suffix='EW'), h


def parseDMS(strDMS, suffix='NSEW', sep=S_SEP):
    '''Parses a string representing deg°min′sec″ to degrees.

       This is very flexible on formats, allowing signed decimal
       degrees, degrees and minutes or degrees minutes and seconds
       optionally suffixed by compass direction NSEW.

       A variety of symbols, separators and suffixes are accepted,
       for example 3° 37′ 09″W.  Minutes and seconds may be omitted.

       @param strDMS: Degrees in any of several forms (string).
       @keyword suffix: Optional, valid compass directions (NSEW).
       @keyword sep: Optional separator between deg°, min′ and sec″ ('').

       @return: Degrees (float).

       @raise ValueError: Invalid strDMS.
    '''
    try:  # signed decimal degrees without NSEW
        return float(strDMS)
    except ValueError:
        pass

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

    return d


def precision(form, prec=None):
    '''Sets the default precison for a given F_ form.

       @param form: F_D, F_DM, F_DMS, F_DEG or F_RAD (string).
       @keyword prec: Optional, number of decimal digits (0..9 or
                      None for default).  Trailing zero decimals
                      are stripped for prec values of 1 and above,
                      but kept for negative prec values.

       @return: Previous precision (int).

       @raise ValueError: Invalid precision.
    '''
    try:
        p  = _F_prec[form]
    except KeyError:
        raise ValueError('%s invalid: %s' % ('form', form))
    if prec is not None:
        if -10 < prec < 10:
            _F_prec[form] = prec
        else:
            raise ValueError('%s invalid: %s' % ('prec', prec))
    return p


def toDMS(deg, form=F_DMS, prec=2, ddd=2, neg='-', pos=''):
    '''Converts signed degrees to string, without suffix.

       @param deg: Degrees to be formatted (scalar).
       @keyword form: F_D, F_DM, F_DMS, F_DEG or F_RAD for deg°,
                      deg°min′, deg°min′sec″ or degrees or radians.
       @keyword prec: Optional, number of decimal digits (0..9 or
                      None for default).  Trailing zero decimals
                      are stripped for prec values of 1 and above,
                      but kept for negative prec values.
       @keyword ddd: Optional, number of digits for deg° (2 or 3).
       @keyword neg: Optional, sign for negative degrees ('-').
       @keyword pos: Optional, sign for positive degrees ('').

       @return: Degrees per the specified form (string).
    '''
    t = _toDMS(deg, form, prec, ddd)
    s = neg if deg < 0 else pos
    return s + t

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
