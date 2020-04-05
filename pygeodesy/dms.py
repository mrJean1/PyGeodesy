
# -*- coding: utf-8 -*-

u'''Functions to parse and format bearing, compass, lat- and longitudes
in various forms of degrees, minutes and seconds.

After I{(C) Chris Veness 2011-2015} published under the same MIT Licence**, see
U{Latitude/Longitude<https://www.Movable-Type.co.UK/scripts/latlong.html>} and
U{Vector-based geodesy<https://www.Movable-Type.co.UK/scripts/latlong-vectors.html>}.

@newfield example: Example, Examples
'''

from pygeodesy.basics import isint
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.named import LatLon2Tuple, LatLon3Tuple
from pygeodesy.streprs import fstr, fstrzs

from math import copysign, radians
try:
    from string import letters as _LETTERS
except ImportError:  # Python 3+
    from string import ascii_letters as _LETTERS

# all public contants, classes and functions
__all__ = _ALL_LAZY.dms
__version__ = '20.04.04'

F_D   = 'd'    #: Format degrees as unsigned "deg°" plus suffix (C{str}).
F_DM  = 'dm'   #: Format degrees as unsigned "deg°min′" plus suffix (C{str}).
F_DMS = 'dms'  #: Format degrees as unsigned "deg°min′sec″" plus suffix (C{str}).
F_DEG = 'deg'  #: Format degrees as unsigned "[D]DD" plus suffix without symbol (C{str}).
F_MIN = 'min'  #: Format degrees as unsigned "[D]DDMM" plus suffix without symbols (C{str}).
F_SEC = 'sec'  #: Format degrees as unsigned "[D]DDMMSS" plus suffix without symbols (C{str}).
F_RAD = 'rad'  #: Convert degrees to radians and format as unsigned "RR" plus suffix (C{str}).

F_D_   = '-d'    #: Format degrees as signed "-/deg°" without suffix (C{str}).
F_DM_  = '-dm'   #: Format degrees as signed "-/deg°min′" without suffix (C{str}).
F_DMS_ = '-dms'  #: Format degrees as signed "-/deg°min′sec″" without suffix (C{str}).
F_DEG_ = '-deg'  #: Format degrees as signed "-/[D]DD" without suffix and symbol (C{str}).
F_MIN_ = '-min'  #: Format degrees as signed "-/[D]DDMM" without suffix and symbols (C{str}).
F_SEC_ = '-sec'  #: Format degrees as signed "-/[D]DDMMSS" without suffix and symbols (C{str}).
F_RAD_ = '-rad'  #: Convert degrees to radians and format as signed "-/RR" without suffix (C{str}).

F_D__   = '+d'    #: Format degrees as signed "-/+deg°" without suffix (C{str}).
F_DM__  = '+dm'   #: Format degrees as signed "-/+deg°min′" without suffix (C{str}).
F_DMS__ = '+dms'  #: Format degrees as signed "-/+deg°min′sec″" without suffix (C{str}).
F_DEG__ = '+deg'  #: Format degrees as signed "-/+[D]DD" without suffix and symbol (C{str}).
F_MIN__ = '+min'  #: Format degrees as signed "-/+[D]DDMM" without suffix and symbols (C{str}).
F_SEC__ = '+sec'  #: Format degrees as signed "-/+[D]DDMMSS" without suffix and symbols (C{str}).
F_RAD__ = '+rad'  #: Convert degrees to radians and format as signed "-/+RR" without suffix (C{str}).

S_DEG = '°'  #: Degrees "°" symbol (C{str}).
S_MIN = '′'  #: Minutes "′" symbol (C{str}).
S_SEC = '″'  #: Seconds "″" symbol (C{str}).
S_RAD = ''   #: Radians symbol "" (C{str}).
S_SEP = ''   #: Separator between deg, min and sec "" (C{str}).

_F_prec = {F_D:   6, F_DM:  4, F_DMS: 2,  #: (INTERNAL) default precs.
           F_DEG: 6, F_MIN: 4, F_SEC: 2, F_RAD: 5}

_S_norm = {'^': S_DEG, '˚': S_DEG,  #: (INTERNAL) normalized DMS.
           "'": S_MIN, '’': S_MIN, '′': S_MIN,
           '"': S_SEC, '″': S_SEC, '”': S_SEC}
_S_ALL  = (S_DEG, S_MIN, S_SEC) + tuple(_S_norm.keys())  #: (INTERNAL) alternates.

_rangerrors = True


class RangeError(ValueError):
    '''Error raised for lat- or longitude values outside the B{C{clip}},
       B{C{clipLat}}, B{C{clipLon}} or B{C{limit}} range in function L{clipDMS},
       L{parse3llh}, L{parseDMS} or L{parseDMS2}.

       @see: Function L{rangerrors}.
    '''
    pass


def _toDMS(deg, form, prec, sep, ddd, suff):  # MCCABE 14
    '''(INTERNAL) Convert degrees to C{str}, with/-out sign and/or suffix.
    '''
    try:
        d = abs(float(deg))
    except ValueError:
        raise ValueError('%s invalid: %r' % ('deg', deg))

    form = form.lower()
    sign = form[:1]
    if sign in '-+':
        form = form[1:]
    else:
        sign = ''

    if prec is None:
        z = p = _F_prec.get(form, 6)
    else:
        z = int(prec)
        p = abs(z)
    w = p + (1 if p else 0)

    if form in (F_DEG, F_MIN, F_SEC):
        s_deg = s_min = s_sec = ''  # no symbols
    else:
        s_deg, s_min, s_sec = S_DEG, S_MIN, S_SEC

    if form in (F_D, F_DEG, 'degrees'):  # deg°, degrees
        t = '%0*.*f' % (ddd+w,p,d)
        s = s_deg

    elif form in (F_RAD, 'radians'):
        t = '%.*F' % (p,radians(d))
        s = S_RAD

    elif form in (F_DM, F_MIN, 'deg+min'):
        d, m = divmod(d * 60, 60)
        d, m = int(d), round(m, p)
        if m >= 60:  # corner case
            m -= 60
            d += 1
        t = "%0*d%s%s%0*.*f" % (ddd,d,s_deg, sep, w+2,p,m)
        s = s_min

    else:  # F_DMS, F_SEC, 'deg+min+sec'
        d, s = divmod(d * 3600, 3600)
        m, s = divmod(s, 60)
        d, m, s = int(d), int(m), round(s, p)
        if s >= 60:  # corner case
            s -= 60
            m += 1
            if m >= 60:
                m -= 60
                d += 1
        t = "%0*d%s%s%02d%s%s%0*.*f" % (ddd,d,s_deg, sep,
                                            m,s_min, sep,
                                      w+2,p,s)
        s = s_sec

    if z > 1:
        t = fstrzs(t)

    if sign:
        if deg < 0:
            t = '-' + t
        elif deg > 0 and sign == '+':
            t = '+' + t
    elif suff:
        s += sep + suff
    return t + s


def bearingDMS(bearing, form=F_D, prec=None, sep=S_SEP):
    '''Convert bearing to a string.

       @arg bearing: Bearing from North (compass C{degrees360}).
       @kwarg form: Optional B{C{bearing}} format (C{str} or L{F_D},
                    L{F_DM}, L{F_DMS}, L{F_DEG}, L{F_MIN}, L{F_SEC},
                    L{F_RAD}, L{F_D_}, L{F_DM_}, L{F_DMS_}, L{F_DEG_},
                    L{F_MIN_}, L{F_SEC_}, L{F_RAD_}, L{F_D__},
                    L{F_DM__}, L{F_DMS__}, L{F_DEG__}, L{F_MIN__},
                    L{F_SEC__} or L{F_RAD__}).
       @kwarg prec: Optional number of decimal digits (0..9 or
                    C{None} for default).  Trailing zero decimals
                    are stripped for B{C{prec}} values of 1 and
                    above, but kept for negative B{C{prec}}.
       @kwarg sep: Optional separator (C{str}).

       @return: Compass degrees per the specified B{C{form}} (C{str}).

       @JSname: I{toBrng}.
    '''
    return _toDMS(bearing % 360, form, prec, sep, 1, '')


def clipDMS(deg, limit):
    '''Clip a lat- or longitude to the given range.

       @arg deg: Unclipped lat- or longitude (C{degrees}).
       @arg limit: Valid B{C{-limit..+limit}} range (C{degrees}).

       @return: Clipped value (C{degrees}).

       @raise RangeError: If B{C{abs(deg)}} beyond B{C{limit}} and
                          L{rangerrors} set to C{True}.
    '''
    if limit > 0:
        c = min(limit, max(-limit, deg))
        if _rangerrors and deg != c:
            raise RangeError('%s beyond %s degrees' % (fstr(deg, prec=6),
                             fstr(copysign(limit, deg), prec=3, ints=True)))
        deg = c
    return deg


def compassDMS(bearing, form=F_D, prec=None, sep=S_SEP):
    '''Convert bearing to a string suffixed with compass point.

       @arg bearing: Bearing from North (compass C{degrees360}).
       @kwarg form: Optional B{C{bearing}} format (C{str} or L{F_D},
                    L{F_DM}, L{F_DMS}, L{F_DEG}, L{F_MIN}, L{F_SEC},
                    L{F_RAD}, L{F_D_}, L{F_DM_}, L{F_DMS_}, L{F_DEG_},
                    L{F_MIN_}, L{F_SEC_}, L{F_RAD_}, L{F_D__},
                    L{F_DM__}, L{F_DMS__}, L{F_DEG__}, L{F_MIN__},
                    L{F_SEC__} or L{F_RAD__}).
       @kwarg prec: Optional number of decimal digits (0..9 or
                    C{None} for default).  Trailing zero decimals
                    are stripped for B{C{prec}} values of 1 and
                    above, but kept for negative B{C{prec}}.
       @kwarg sep: Optional separator (C{str}).

       @return: Compass degrees and point in the specified form (C{str}).
    '''
    return _toDMS(bearing % 360, form, prec, sep, 1, compassPoint(bearing))


def compassPoint(bearing, prec=3):
    '''Convert bearing to a compass point.

       @arg bearing: Bearing from North (compass C{degrees360}).
       @kwarg prec: Optional precision (1 for cardinal or basic winds,
                    2 for intercardinal or ordinal or principal winds,
                    3 for secondary-intercardinal or half-winds or
                    4 for quarter-winds).

       @return: Compass point (1-, 2-, 3- or 4-letter C{str}).

       @raise ValueError: Invalid B{C{prec}}.

       @see: U{Dms.compassPoint
             <https://GitHub.com/chrisveness/geodesy/blob/master/dms.js>}
             and U{Compass rose<https://WikiPedia.org/wiki/Compass_rose>}.

       @example:

       >>> p = compassPoint(24, 1)  # 'N'
       >>> p = compassPoint(24, 2)  # 'NE'
       >>> p = compassPoint(24, 3)  # 'NNE'
       >>> p = compassPoint(24)     # 'NNE'
       >>> p = compassPoint(11, 4)  # 'NbE'
       >>> p = compassPoint(30, 4)  # 'NEbN'

       >>> p = compassPoint(11.249)  # 'N'
       >>> p = compassPoint(11.25)   # 'NNE'
       >>> p = compassPoint(-11.25)  # 'N'
       >>> p = compassPoint(348.749) # 'NNW'
    '''
    try:  # m = 2 << prec; x = 32 // m
        m, x = _MOD_X[prec]
    except KeyError:
        raise ValueError('%s invalid: %r' % ('prec', prec))
    # not round(), i.e. half-even rounding in Python 3,
    # but round-away-from-zero as int(b + 0.5) iff b is
    # non-negative, otherwise int(b + copysign(0.5, b))
    q = int((bearing % 360) * m / 360.0 + 0.5) % m
    return _WINDS[q * x]


_MOD_X = {1: (4, 8), 2: (8, 4), 3: (16, 2), 4: (32, 1)}  #: (INTERNAL) [prec]
_WINDS = ('N', 'NbE', 'NNE', 'NEbN', 'NE', 'NEbE', 'ENE', 'EbN',
          'E', 'EbS', 'ESE', 'SEbE', 'SE', 'SEbS', 'SSE', 'SbE',
          'S', 'SbW', 'SSW', 'SWbS', 'SW', 'SWbW', 'WSW', 'WbS',
          'W', 'WbN', 'WNW', 'NWbW', 'NW', 'NWbN', 'NNW', 'NbW')  #: (INTERNAL) cardinals


def degDMS(deg, prec=6, s_D=S_DEG, s_M=S_MIN, s_S=S_SEC, neg='-', pos=''):
    '''Convert degrees to a string in degrees, minutes B{I{or}} seconds.

       @arg deg: Value in degrees (C{scalar}).
       @kwarg prec: Optional number of decimal digits (0..9 or
                    C{None} for default).  Trailing zero decimals
                    are stripped for B{C{prec}} values of 1 and
                    above, but kept for negative B{C{prec}}.
       @kwarg s_D: Symbol for degrees (C{str}).
       @kwarg s_M: Symbol for minutes (C{str}) or C{""}.
       @kwarg s_S: Symbol for seconds (C{str}) or C{""}.
       @kwarg neg: Optional sign for negative (C{'-'}).
       @kwarg pos: Optional sign for positive (C{''}).

       @return: I{Either} degrees, minutes B{I{or}} seconds (C{str}).
    '''
    d, s = abs(deg), s_D
    if d < 1:
        if s_M:
            d *= 60
            if d < 1 and s_S:
                d *= 60
                s = s_S
            else:
                s = s_M
        elif s_S:
            d *= 3600
            s = s_S

    n = neg if deg < 0 else pos
    z = int(prec)
    t = '%s%.*F' % (n, abs(z),d)
    if z > 1:
        t = fstrzs(t)
    return t + s


def latDMS(deg, form=F_DMS, prec=2, sep=S_SEP):
    '''Convert latitude to a string, optionally suffixed with N or S.

       @arg deg: Latitude to be formatted (C{degrees}).
       @kwarg form: Optional B{C{deg}} format (C{str} or L{F_D},
                    L{F_DM}, L{F_DMS}, L{F_DEG}, L{F_MIN}, L{F_SEC},
                    L{F_RAD}, L{F_D_}, L{F_DM_}, L{F_DMS_}, L{F_DEG_},
                    L{F_MIN_}, L{F_SEC_}, L{F_RAD_}, L{F_D__},
                    L{F_DM__}, L{F_DMS__}, L{F_DEG__}, L{F_MIN__},
                    L{F_SEC__} or L{F_RAD__}).
       @kwarg prec: Optional number of decimal digits (0..9 or
                    C{None} for default).  Trailing zero decimals
                    are stripped for B{C{prec}} values of 1 and
                    above, but kept for negative B{C{prec}}.
       @kwarg sep: Optional separator (C{str}).

       @return: Degrees in the specified form (C{str}).

       @JSname: I{toLat}.
    '''
    return _toDMS(deg, form, prec, sep, 2, 'S' if deg < 0 else 'N')


def lonDMS(deg, form=F_DMS, prec=2, sep=S_SEP):
    '''Convert longitude to a string, optionally suffixed with E or W.

       @arg deg: Longitude to be formatted (C{degrees}).
       @kwarg form: Optional B{C{deg}} format (C{str} or L{F_D},
                    L{F_DM}, L{F_DMS}, L{F_DEG}, L{F_MIN}, L{F_SEC},
                    L{F_RAD}, L{F_D_}, L{F_DM_}, L{F_DMS_}, L{F_DEG_},
                    L{F_MIN_}, L{F_SEC_}, L{F_RAD_}, L{F_D__},
                    L{F_DM__}, L{F_DMS__}, L{F_DEG__}, L{F_MIN__},
                    L{F_SEC__} or L{F_RAD__}).
       @kwarg prec: Optional number of decimal digits (0..9 or
                    C{None} for default).  Trailing zero decimals
                    are stripped for B{C{prec}} values of 1 and
                    above, but kept for negative B{C{prec}}.
       @kwarg sep: Optional separator (C{str}).

       @return: Degrees in the specified form (C{str}).

       @JSname: I{toLon}.
    '''
    return _toDMS(deg, form, prec, sep, 3, 'W' if deg < 0 else 'E')


def normDMS(strDMS, norm=''):
    '''Normalize all degree ˚, minute ' and second " symbols in a
       string to the default symbols %s, %s and %s.

       @arg strDMS: DMS (C{str}).
       @kwarg norm: Optional replacement symbol, default symbol
                    otherwise (C{str}).

       @return: Normalized DMS (C{str}).
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

       @arg strll: Latitude, longitude[, height] (C{str}, ...).
       @kwarg height: Optional, default height (C{meter}).
       @kwarg sep: Optional separator (C{str}).
       @kwarg clipLat: Keep latitude in B{C{-clipLat..+clipLat}} (C{degrees}).
       @kwarg clipLon: Keep longitude in B{C{-clipLon..+clipLon}} range (C{degrees}).

       @return: A L{LatLon3Tuple}C{(lat, lon, height)} in
                C{degrees}, C{degrees} and C{float}.

       @raise RangeError: Lat- or longitude value of B{C{strll}} outside
                          valid range and L{rangerrors} set to C{True}.

       @raise ValueError: Invalid B{C{strll}}.

       @see: Functions L{parseDMS} and L{parseDMS2} for more details
             on the forms and symbols accepted.

       @example:

       >>> parse3llh('000°00′05.31″W, 51° 28′ 40.12″ N')
       (51.4778°N, 000.0015°W, 0)
    '''
    try:
        ll = strll.strip().split(sep)
        if len(ll) > 2:  # XXX interpret height unit
            h = float(ll.pop(2).strip().rstrip(_LETTERS).rstrip())
        else:
            h = height
        if len(ll) != 2:
            raise ValueError
    except (AttributeError, TypeError, ValueError):
        return ValueError('parsing %r failed' % (strll,))

    a, b = [_.strip() for _ in ll]
    if a[-1:] in 'EW' or b[-1:] in 'NS':
        a, b = b, a
    a, b = parseDMS2(a, b, clipLat=clipLat, clipLon=clipLon)  # PYCHOK LatLon2Tuple
    return LatLon3Tuple(a, b, h)


def parseDMS(strDMS, suffix='NSEW', sep=S_SEP, clip=0):
    '''Parse a string representing deg°min′sec″ to degrees.

       This is very flexible on formats, allowing signed decimal
       degrees, degrees and minutes or degrees minutes and seconds
       optionally suffixed by compass direction NSEW.

       A variety of symbols, separators and suffixes are accepted,
       for example 3° 37′ 09″W.  Minutes and seconds may be omitted.

       @arg strDMS: Degrees in any of several forms (C{str} or C{degrees}).
       @kwarg suffix: Optional, valid compass directions (NSEW).
       @kwarg sep: Optional separator between deg°, min′ and sec″ ('').
       @kwarg clip: Optionally, limit value to -clip..+clip (C{degrees}).

       @return: Degrees (C{float}).

       @raise RangeError: Value of B{C{strDMS}} outside the valid range
                          and L{rangerrors} set to C{True}.

       @raise ValueError: Invalid B{C{strDMS}}.

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

        except (AttributeError, IndexError, TypeError, ValueError):
            raise ValueError('parsing %r failed' % (strDMS,))

    return clipDMS(d, float(clip))


def parseDMS2(strLat, strLon, sep=S_SEP, clipLat=90, clipLon=180):
    '''Parse lat- and longitude representions.

       @arg strLat: Latitude in any of several forms (C{str} or C{degrees}).
       @arg strLon: Longitude in any of several forms (C{str} or C{degrees}).
       @kwarg sep: Optional separator between deg°, min′ and sec″ ('').
       @kwarg clipLat: Keep latitude in B{C{-clipLat..+clipLat}} range (C{degrees}).
       @kwarg clipLon: Keep longitude in B{C{-clipLon..+clipLon}} range (C{degrees}).

       @return: A L{LatLon2Tuple}C{(lat, lon)} in C{degrees}.

       @raise RangeError: Value of B{C{strLat}} or B{C{strLon}} outside the
                          valid range and L{rangerrors} set to C{True}.

       @raise ValueError: Invalid B{C{strLat}} or B{C{strLon}}.

       @see: Function L{parse3llh} to parse a string with lat-,
             longitude and height values and function L{parseDMS}
             to parse individual lat- or longitudes.
    '''
    return LatLon2Tuple(parseDMS(strLat, suffix='NS', sep=sep, clip=clipLat),
                        parseDMS(strLon, suffix='EW', sep=sep, clip=clipLon))


def _parseUTMUPS(strUTMUPS, band=''):  # see .utm.py
    '''(INTERNAL) Parse a string representing a UTM or UPS coordinate
       consisting of C{"zone[band] hemisphere/pole easting northing"}.

       @arg strUTMUPS: A UTM or UPS coordinate (C{str}).
       @kwarg band: Optional, default Band letter (C{str}).

       @return: 5-Tuple (C{zone, hemisphere/pole, easting, northing,
                band}).

       @raise Value: Invalid B{C{strUTMUPS}}.
    '''
    try:
        u = strUTMUPS.strip().replace(',', ' ').split()
        if len(u) < 4:
            raise ValueError

        z, h = u[:2]
        if h[:1] not in 'NnSs':
            raise ValueError

        if z.isdigit():
            z, B = int(z), band
        else:
            for i in range(len(z)):
                if not z[i].isdigit():
                    # int('') raises ValueError
                    z, B = int(z[:i]), z[i:]
                    break
            else:
                raise ValueError

        e, n = map(float, u[2:4])

    except (AttributeError, TypeError, ValueError):
        raise ValueError('%s invalid: %r' % ('strUTMUPS', strUTMUPS))

    return z, h.upper(), e, n, B.upper()


def precision(form, prec=None):
    '''Set the default precison for a given F_ form.

       @arg form: L{F_D}, L{F_DM}, L{F_DMS}, L{F_DEG}, L{F_MIN},
                  L{F_SEC} or L{F_RAD} (C{str}).
       @kwarg prec: Optional number of decimal digits (0..9 or
                    C{None} for default).  Trailing zero decimals
                    are stripped for B{C{prec}} values of 1 and
                    above, but kept for negative B{C{prec}}.

       @return: Previous precision (C{int}).

       @raise ValueError: Invalid B{C{form}} or B{C{prec}} or beyond valid range.
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
    '''Get/set raising of range errors.

       @kwarg raiser: Choose C{True} to raise or C{False} to ignore
                      L{RangeError} exceptions.  Use C{None} to leave
                      the setting unchanged.

       @return: Previous setting (C{bool}).

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

       @arg deg: Degrees to be formatted (C{degrees}).
       @kwarg form: Optional B{C{deg}} format (C{str} or L{F_D},
                    L{F_DM}, L{F_DMS}, L{F_DEG}, L{F_MIN}, L{F_SEC},
                    L{F_RAD} without suffix, L{F_D_}, L{F_DM_},
                    L{F_DMS_}, L{F_DEG_}, L{F_MIN_}, L{F_SEC_},
                    L{F_RAD_}, L{F_D__}, L{F_DM__}, L{F_DMS__},
                    L{F_DEG__}, L{F_MIN__}, L{F_SEC__} or L{F_RAD__}).
       @kwarg prec: Optional number of decimal digits (0..9 or
                    C{None} for default).  Trailing zero decimals
                    are stripped for B{C{prec}} values of 1 and
                    above, but kept for negative B{C{prec}}.
       @kwarg sep: Optional separator (C{str}).
       @kwarg ddd: Optional number of digits for deg° (2 or 3).
       @kwarg neg: Optional sign for negative degrees ('-').
       @kwarg pos: Optional sign for positive degrees ('').

       @return: Degrees in the specified form (C{str}).
    '''
    t = _toDMS(deg, form, prec, sep, ddd, '')
    if form[:1] not in '-+':
        t = (neg if deg < 0 else pos) + t
    return t

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
