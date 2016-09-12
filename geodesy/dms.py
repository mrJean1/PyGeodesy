
# -*- coding: utf-8 -*-

# Common base class and functions for modules geodesySphere and -Vector.

# After (C) Chris Veness 2011-2015 published under the same MIT Licence,
# see <http://www.movable-type.co.uk/scripts/latlong.html>
# and <http://www.movable-type.co.uk/scripts/latlong-vectors.html>

from math import radians
try:
    from string import letters as LETTERS
except ImportError:  # Python 3+
    from string import ascii_letters as LETTERS

# all public contants, classes and functions
__all__ = ('F_D', 'F_DM', 'F_DMS', 'F_RAD',  # format contants
           'S_DEG', 'S_MIN', 'S_SEC', 'S_SEP',  # symbols
           'bearingDMS', 'compassDMS', 'compassPoint',  # functions
           'latDMS', 'lonDMS', 'normDMS',
           'parseDMS', 'parse3llh', 'precision', 'toDMS')
__version__ = '16.09.05'

F_D   = 'd'    # format degrees as deg°
F_DM  = 'dm'   # format degrees as deg°min'
F_DMS = 'dms'  # format degrees as deg°min'sec"
F_RAD = 'rad'  # convert degrees to radians and format

S_DEG = '°'  # degree symbol
S_MIN = '′'  # minutes symbol
S_SEC = '″'  # seconds symbol
S_SEP = ''   # separator between deg°, min′ and sec″

_F_prec = {F_D: 6, F_DM: 4, F_DMS: 2}  # default precs

_S_norm = {'^': S_DEG, '˚': S_DEG,  # normalized DMS
           "'": S_MIN, '’': S_MIN, '′': S_MIN,
           '"': S_SEC, '″': S_SEC, '”': S_SEC}
_S_ALL  = (S_DEG, S_MIN, S_SEC) + tuple(_S_norm.keys())


def _toDMS(deg, form=F_DMS, prec=None, ddd=3):
    '''Convert numeric degrees, without sign or suffix.
    '''
    try:
        d = abs(float(deg))
    except ValueError:
        raise ValueError('%s not numeric: %r' % ('deg', deg))

    if prec is None:
        z = _F_prec.get(form, 6)
    else:
        z = int(prec) or 0
    p = abs(z)
    w = p + (1 if p else 0)

    f = form.lower()
    if f in (F_D, 'deg'):
        t = '%0*.*f' % (ddd+w,p,d)
        s = S_DEG

    elif f in (F_RAD, 'radians'):
        t = '%.*f' % (p, radians(d))
        s = ''

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

    if z > 1 and t.endswith('0'):
        # strip trailing decimal zeros, except one
        z = len(t) - z + 1
        t = t[:z] + t[z:].rstrip('0')

    return t + s


def bearingDMS(deg, form=F_D, prec=6):
    '''Convert numeric degrees to a bearing 0°..360°.

       @param {number} deg - Degrees to be formatted as specified.
       @param {string} [form=F_D] - Use F_D, F_DM, F_DMS for deg°, deg°min', deg°min'sec".
       @param {number} [prec=0..8] - Number of decimal digits.

       @returns {string} Degrees formatted per the specified form.
    '''
    return _toDMS(deg % 360, form=form, prec=prec, ddd=3)

toBrng = bearingDMS  # XXX original name

_COMPASS = ('N', 'NNE', 'NE', 'ENE',
            'E', 'ESE', 'SE', 'SSE',
            'S', 'SSW', 'SW', 'WSW',
            'W', 'WNW', 'NW', 'NNW')

_M_X = {1: (4, 4), 2: (8, 2), 3: (16, 1)}


def compassDMS(deg, form=F_D, prec=None):
    '''Convert numeric degrees to a bearing 0°..360° and compass point.

       @param {number} deg - Degrees to be formatted as specified.
       @param {string} [form=F_D] - Use F_D, F_DM, F_DMS for deg°, deg°min', deg°min'sec".
       @param {number} [prec=0..8] - Number of decimal digits.

       @returns {string} Formatted bearing and compass point.
    '''
    return bearingDMS(deg , form=form, prec=prec) + S_SEP + compassPoint(deg)


def compassPoint(bearing, prec=3):
    '''Returns compass point (to given precision) for supplied bearing.

       @param {number} bearing - Bearing in degrees from North.
       @param {number} [prec=3] - Precision (1: cardinal,
                                             2: intercardinal or
                                             3: secondary-intercardinal).

       @returns {string} Compass point for supplied bearing.

       @example
       p = compass(24)     # 'NNE'
       p = compass(24, 1)  # 'N'
    '''
    try:
        m, x = _M_X[prec]
    except KeyError:
        raise ValueError('invalid prec: %r' % (prec,))

    q = int(round((bearing % 360) * m / 360.0)) % m
    return _COMPASS[q * x]


def latDMS(deg, form=F_DMS, prec=2):
    '''Convert latitude degrees, suffixed with N or S.

       @param {number} deg - Degrees to be formatted as specified.
       @param {string} [form=F_DMS] - Use F_D, F_DM, F_DMS for deg°, deg°min', deg°min'sec".
       @param {number} [prec=6] - Number of decimal digits (0..8).

       @returns {string} Degrees formatted per the specified form.
    '''
    return _toDMS(deg, form=form, prec=prec, ddd=2) + ('S' if deg < 0 else 'N')

toLat = latDMS  # XXX original name


def lonDMS(deg, form=F_DMS, prec=2):
    '''Convert longitude degrees, suffixed with E or W.

       @param {number} deg - Degrees to be formatted as specified.
       @param {string} [form=F_DMS] - Use F_D, F_DM, F_DMS for deg°, deg°min', deg°min'sec".
       @param {number} [prec=6] - Number of decimal digits (0..8).

       @returns {string} Degrees formatted per the specified form.
    '''
    return _toDMS(deg, form=form, prec=prec, ddd=3) + ('W' if deg < 0 else 'E')

toLon = lonDMS  # XXX original name


def normDMS(strDMS):
    '''Normalize the ˚, ' and " symbols in a string
       to the defaults %s, %s and %s.

       @param {string} strDMS - String.

       @returns {string} Normalized DMS string.
    '''
    for s, S in _S_norm.items():
        strDMS = strDMS.replace(s, S)
    return strDMS

normDMS.__doc__ %= S_DEG, S_MIN, S_SEC


def parse3llh(strll, height=0, sep=','):
    '''Parse a string representing lat-/longitude point and return a
       3-tuple (lat, lon, height).

       The lat- and longitude must be separated by a comma.  If height
       is present it must follow the longitude, separated by another
       comma.

        The lat- and longitude must be separated by a sep[arator]
        character.  If height is present it must follow and be
        separated by another sep[arator].  Lat- and longitude may
        be swapped, provided at least one ends with the proper
        compass direction.

       See parseDMS for more details on the forms accepted.

       @param {string} strll - Latitude, longitude [, height].
       @param {meter} [height=0] - Default height in meter.
       @param {string} [sep=','] - Separator.

       @returns {(degrees, degrees, float)} 3-Tuple (lat, lon, height).

       @throws {ValueError} Invalid strll.

       @example
       t = parse3llh('000°00′05.31″W, 51° 28′ 40.12″ N')  # (51.4778°N, 000.0015°W, 0)
    '''
    ll = strll.strip().split(sep)
    if len(ll) > 2:  # XXX interpret height unit
        h = float(ll.pop(2).strip().rstrip(LETTERS).rstrip())
    else:
        h = height
    if len(ll) != 2:
        return ValueError('parsing %r' % (strll,))

    a, b = [_.strip() for _ in ll]
    if a[-1:] in 'EW' or b[-1:] in 'NS':
        a, b = b, a
    return parseDMS(a, suffix='NS'), parseDMS(b, suffix='EW'), h


def parseDMS(strDMS, suffix='NSEW'):
    '''Parse a string representing deg°min'sec" into degrees.

       This is very flexible on formats, allowing signed decimal
       degrees, degrees and minutes or degrees minutes and seconds
       optionally suffixed by compass direction NSEW.  A variety of
       separators are accepted, for example 3° 37' 09"W.  Minutes
       and seconds may be omitted.

       @param {string|number} strDMS - Degrees in one of several forms.
       @param {string} [suffix='NSEW'] - Valid compass directions.

       @returns {degrees} Value as a float.

       @throws {ValueError} Invalid strDMS.
    '''
    try:  # signed decimal degrees without NSEW
        return float(strDMS)
    except ValueError:
        pass

    try:
        strDMS = strDMS.strip()

        t = strDMS.lstrip('-+').rstrip(suffix.upper())
        if S_SEP:
            t = t.replace(S_SEP, ' ')
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
        raise ValueError('parsing %r' % (strDMS,))

    return d


def precision(form, prec=None):
    '''Set the default precison for a given format.

       @param {string} form - F_D, F_DM, F_DMS.
       @param {number} [prec=None] - Number of decimal digits
                                  or None to remain unchanged.

       @returns {number} Previous precision.
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
    '''Convert signed degrees, without suffix.

       @param {number} deg - Degrees to be formatted as specified.
       @param {string} [form=F_DMS] - F_D, F_DM, F_DMS, F_RAD for deg°, deg°min', deg°min'sec".
       @param {number} [prec=6] - Number of decimal digits (0..8).
       @param {number} [ddd=3] - Number of digits for deg° (2|3).
       @param {string} [neg='-'] - Sign for negative degrees.
       @param {string} [pos=''] - Sign for positive degrees.

       @returns {string} Degrees formatted per the specified form.
    '''
    t = _toDMS(deg, form=form, prec=prec, ddd=ddd)
    s = neg if deg < 0 else pos
    return s + t


if __name__ == '__main__':

    from tests import Tests as _Tests

    class Tests(_Tests):

        def testDMS(self):
            self.test('parseDMS', parseDMS(  '0.0°'), '0.0')
            self.test('parseDMS', parseDMS(    '0°'), '0.0')
            self.test('parseDMS', parseDMS('''000°00'00"'''),   '0.0')
            self.test('parseDMS', parseDMS('''000°00'00.0"'''), '0.0')
            self.test('parseDMS', parseDMS('''000° 00'00"'''),    '0.0')
            self.test('parseDMS', parseDMS('''000°00 ' 00.0"'''), '0.0')

            x = parse3llh('000° 00′ 05.31″W, 51° 28′ 40.12″ N')
            x = ', '.join('%.6f' % a for a in x)  # XXX fStr
            self.test('parse3llh', x, '51.477811, -0.001475, 0.000000')

            for a, x in (((),            '''45°45'45.36"'''),
                         ((F_D, None),     '45.7626°'),
                         ((F_DM, None),    "45°45.756'"),
                         ((F_DMS, None), '''45°45'45.36"'''),
                         ((F_D, 6),     '45.7626°'),
                         ((F_DM, -4),   "45°45.7560'"),
                         ((F_DMS, 2), '''45°45'45.36"''')):
                self.test('toDMS', toDMS(45.76260, *a), x)

            for a, x in (((1,),   'N'),
                         ((0,),   'N'),
                         ((-1,),  'N'),
                         ((359,), 'N'),
                         ((24,),   'NNE'),
                         ((24, 1), 'N'),
                         ((24, 2), 'NE'),
                         ((24, 3), 'NNE'),
                         ((226,),   'SW'),
                         ((226, 1), 'W'),
                         ((226, 2), 'SW'),
                         ((226, 3), 'SW'),
                         ((237,),   'WSW'),
                         ((237, 1), 'W'),
                         ((237, 2), 'SW'),
                         ((237, 3), 'WSW')):
                self.test('compassPoint', compassPoint(*a), x)

    t = Tests(__file__, __version__)
    t.testDMS()
    t.results()

    # Typical test results (on MacOS X):

    # testing dms.py version 16.09.05
    # test 1 parseDMS: 0.0
    # test 2 parseDMS: 0.0
    # test 3 parseDMS: 0.0
    # test 4 parseDMS: 0.0
    # test 5 parseDMS: 0.0
    # test 6 parseDMS: 0.0
    # test 7 parse3llh: 51.477811, -0.001475, 0.000000
    # test 8 toDMS: 45°45′45.36″
    # test 9 toDMS: 45.7626°
    # test 10 toDMS: 45°45.756′
    # test 11 toDMS: 45°45′45.36″
    # test 12 toDMS: 45.7626°
    # test 13 toDMS: 45°45.7560′
    # test 14 toDMS: 45°45′45.36″
    # test 15 compassPoint: N
    # test 16 compassPoint: N
    # test 17 compassPoint: N
    # test 18 compassPoint: N
    # test 19 compassPoint: NNE
    # test 20 compassPoint: N
    # test 21 compassPoint: NE
    # test 22 compassPoint: NNE
    # test 23 compassPoint: SW
    # test 24 compassPoint: W
    # test 25 compassPoint: SW
    # test 26 compassPoint: SW
    # test 27 compassPoint: WSW
    # test 28 compassPoint: W
    # test 29 compassPoint: SW
    # test 30 compassPoint: WSW
    # all dms.py tests passed (Python 2.7.10)

    # testing dms.py version 16.09.05
    # test 1 parseDMS: 0.0
    # test 2 parseDMS: 0.0
    # test 3 parseDMS: 0.0
    # test 4 parseDMS: 0.0
    # test 5 parseDMS: 0.0
    # test 6 parseDMS: 0.0
    # test 7 parse3llh: 51.477811, -0.001475, 0.000000
    # test 8 toDMS: 45°45′45.36″
    # test 9 toDMS: 45.7626°
    # test 10 toDMS: 45°45.756′
    # test 11 toDMS: 45°45′45.36″
    # test 12 toDMS: 45.7626°
    # test 13 toDMS: 45°45.7560′
    # test 14 toDMS: 45°45′45.36″
    # test 15 compassPoint: N
    # test 16 compassPoint: N
    # test 17 compassPoint: N
    # test 18 compassPoint: N
    # test 19 compassPoint: NNE
    # test 20 compassPoint: N
    # test 21 compassPoint: NE
    # test 22 compassPoint: NNE
    # test 23 compassPoint: SW
    # test 24 compassPoint: W
    # test 25 compassPoint: SW
    # test 26 compassPoint: SW
    # test 27 compassPoint: WSW
    # test 28 compassPoint: W
    # test 29 compassPoint: SW
    # test 30 compassPoint: WSW
    # all dms.py tests passed (Python 3.5.1)
